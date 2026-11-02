#!/usr/bin/env python3
"""Aggregate the C2 hydration analysis to replica-level statistics.

Reads campaign_100x40ns/analysis/c2_hydration/<cond>/<sys>/replica<nnn>/*.dat and reports,
for each observable: the per-condition mean with a REPLICA-level standard error, the
OPC - TIP3P difference with a Welch standard error, and -- from the observed between-replica
SD -- the number of replicas per arm that an 80% powered test of that difference would need.

Replica-level errors are used deliberately: pooling frames would understate the uncertainty
in exactly the way Reviewer 1 Comment 5 criticised. 400 frames at 100 ps are independent
solvent configurations, so there is no within-replica correlation to exploit anyway.
"""
import os, glob, math
import numpy as np

A = '/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/c2_hydration'
CONDS = ['opc-salt', 'tip3p']
SYSTEMS = ['lats1_kinase', 'ndr2_kinase']

def cols(path):
    """Read a cpptraj data file into (header_names, array).

    Only the leading run of numeric fields is kept. hbond.dat ends with a PW[ID] column
    holding strings like "1073(37+38+),1276(15+18+)", which must not reach float().
    """
    if not os.path.exists(path):
        return None, None
    names, rows, width = None, [], None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#'):
                if names is None:
                    names = ln.lstrip('#').split()
                continue
            p = ln.split()
            if not p:
                continue
            vals = []
            for x in p:
                try:
                    vals.append(float(x))
                except ValueError:
                    break
            if not vals:
                continue
            if width is None:
                width = len(vals)
            rows.append(vals[:width])
    if not rows:
        return names, None
    return names, np.asarray(rows)


def col_by_name(names, arr, suffix):
    """Column whose header name ends with `suffix` (e.g. '[UV]'), or None."""
    if names is None or arr is None:
        return None
    for i, nm in enumerate(names):
        if nm.endswith(suffix):
            return arr[:, i] if i < arr.shape[1] else None
    return None

def peak(path, rmin=1.5, rmax=5.0):
    """First-peak height and position of an RDF within [rmin, rmax]."""
    _, a = cols(path)
    if a is None or a.shape[1] < 2:
        return np.nan, np.nan
    r, g = a[:, 0], a[:, 1]
    m = (r >= rmin) & (r <= rmax)
    if not m.any() or not np.isfinite(g[m]).any():
        return np.nan, np.nan
    i = np.nanargmax(g[m])
    return float(g[m][i]), float(r[m][i])

def replica(cond, sysn, n):
    """Per-replica means of every observable. Returns {} if the replica is incomplete."""
    d = f'{A}/{cond}/{sysn}/replica{n}'
    if not os.path.exists(f'{d}/sasa.dat'):
        return {}
    o = {}
    nm, hb = cols(f'{d}/hbond.dat')
    if hb is not None:
        # header is: #Frame PW[UU] PW[UV] PW[Bridge] PW[ID]
        # UU = solute-solute, UV = solute-solvent, Bridge = waters bridging two solute sites.
        # Select by name: PW[ID] is non-numeric and the column count is not fixed.
        for tag, suf in [('hbond_solute_solute', '[UU]'),
                         ('hbond_solute_solvent', '[UV]'),
                         ('water_bridges', '[Bridge]')]:
            c = col_by_name(nm, hb, suf)
            if c is not None:
                o[tag] = c.mean()
    for tag, f in [('protein', 'ws_protein'), ('cationic', 'ws_cationic'), ('anionic', 'ws_anionic')]:
        nmw, w = cols(f'{d}/{f}.dat')
        if w is not None and w.shape[1] >= 3:
            # header: #Frame WS_000nn[lower] WS_000nn[upper]
            lo = col_by_name(nmw, w, '[lower]')
            hi = col_by_name(nmw, w, '[upper]')
            o[f'shell1_{tag}'] = (lo if lo is not None else w[:, 1]).mean()
            o[f'shell2_{tag}'] = (hi if hi is not None else w[:, 2]).mean()
    _, s = cols(f'{d}/sasa.dat')
    if s is not None:
        o['sasa'] = s[:, 1].mean()
    for tag, f, lo, hi in [('rdf_cat_watO', 'rdf_cat_watO', 1.5, 5.0),
                           ('rdf_ani_watO', 'rdf_ani_watO', 1.5, 5.0),
                           ('rdf_cat_Cl',  'rdf_cat_Cl',  2.0, 6.0),
                           ('rdf_ani_Na',  'rdf_ani_Na',  1.5, 5.0)]:
        h, p = peak(f'{d}/{f}.dat', lo, hi)
        o[f'{tag}_peak'] = h
        o[f'{tag}_rpeak'] = p
    return o

def main():
    data = {}
    for cond in CONDS:
        for sysn in SYSTEMS:
            for d in sorted(glob.glob(f'{A}/{cond}/{sysn}/replica*')):
                n = d.rsplit('replica', 1)[1]
                r = replica(cond, sysn, n)
                if r:
                    data.setdefault((cond, sysn), {})[n] = r
    if not data:
        print('no completed replicas yet'); return
    print('C2 hydration: replica-level statistics\n')
    for (cond, sysn), reps in sorted(data.items()):
        print(f'  {cond:9s} {sysn:13s} replicas: {len(reps)}  ({",".join(sorted(reps))})')
    keys = sorted({k for reps in data.values() for r in reps.values() for k in r})

    for sysn in SYSTEMS:
        have = [(c, sysn) in data for c in CONDS]
        if not all(have):
            continue
        print(f'\n{"="*104}\n{sysn}  --  OPC vs TIP3P, replica-level\n{"="*104}')
        print(f'{"observable":26s} {"OPC mean+/-SE":>20s} {"TIP3P mean+/-SE":>20s} '
              f'{"diff+/-SE":>17s} {"t":>6s} {"n/arm 80%":>10s}')
        for k in keys:
            vals = {}
            for c in CONDS:
                v = np.array([r[k] for r in data[(c, sysn)].values() if k in r and np.isfinite(r[k])])
                vals[c] = v
            o, t = vals['opc-salt'], vals['tip3p']
            if len(o) < 2 or len(t) < 2:
                continue
            so, st = o.std(ddof=1), t.std(ddof=1)
            seo, set_ = so / math.sqrt(len(o)), st / math.sqrt(len(t))
            diff = o.mean() - t.mean()
            se = math.hypot(seo, set_)
            tstat = diff / se if se > 0 else float('nan')
            sp = math.sqrt((so**2 + st**2) / 2)
            if sp > 0 and diff != 0:
                # floor at 2: a two-sample test needs at least two replicas per arm
                n80 = max(2, math.ceil(2 * (1.96 + 0.8416)**2 / (diff / sp)**2))
                n80s = str(n80)
            elif sp == 0 and diff != 0:
                # both arms landed in a single RDF bin; the 0.1 A binning, not sampling,
                # limits this comparison -- report it as bin-resolution-limited
                n80s = 'bin-ltd'
            else:
                n80s = '-'
            tstr = f'{tstat:+6.2f}' if math.isfinite(tstat) else '   n/a'
            print(f'{k:26s} {o.mean():12.3f}+/-{seo:6.3f} {t.mean():12.3f}+/-{set_:6.3f} '
                  f'{diff:+10.3f}+/-{se:5.3f} {tstr} {n80s:>10}')
        print('\n  n/arm 80% is computed from the OBSERVED between-replica SD, so it is the empirical')
        print('  answer to the sample-size question -- not the ENPOLAR proxy used before C2 ran.')

if __name__ == '__main__':
    main()
