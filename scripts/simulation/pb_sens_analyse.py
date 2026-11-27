#!/usr/bin/env python3
"""Analyse the PB istrng sensitivity runs.

Parses MMPBSA.py's raw -eo per_frame.csv, which contains four blocks in order:
Complex, Receptor, Ligand, DELTA. (The campaign's deposited mmpbsa_perframe_*.csv.gz files are
a TIDIED form of this produced by mmpbsa_tidy.py, with cpx_/rec_/lig_ columns; the raw files
here are not in that layout.)

The two protocols are handled differently, because the receptor term differs between them.

  ST  (the number quoted in the paper's Methods): MOB1A is taken from each complex's own
      snapshots, so it is system-specific and does NOT cancel between LATS1 and NDR2. The
      binding energy of a replica is MMPBSA.py's own DELTA block (complex - receptor - ligand).
  MT  the unbound MOB1A reference is POOLED over both systems, so the shared term cancels in
      the LATS1-minus-NDR2 difference. The MT arm deposited here was scored on the unbound
      runs that preceded the matched-salt correction; it is historical, is not used in the
      paper, and is kept only for completeness.

(Release v2.0.0 of this script applied the pooled-reference cancellation to the ST arm as well,
which gave -0.015 +- 0.011; the correct ST value, printed below, is -0.021 +- 0.010 kcal/mol.)

The headline number is the PAIRED per-replica difference between istrng arms -- the same frames
are scored in both, so frame-to-frame noise cancels. Absolute ddH at 10 snapshots/replica is
noisy and shown only for orientation.
"""
import csv, glob, math, os, sys, collections

def blocks(path):
    """return {'complex':[totals], 'receptor':[...], 'ligand':[...], 'delta':[...]}"""
    names = ['complex', 'receptor', 'ligand', 'delta']
    out, cur, idx = {}, None, -1
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row: continue
            if row[0].startswith('Frame #'):
                idx += 1
                if idx < len(names): cur = names[idx]; out[cur] = []
                ti = len(row) - 1          # TOTAL / DELTA TOTAL is the last column
                continue
            if cur is None: continue
            try: out[cur].append(float(row[-1]))
            except (ValueError, IndexError): pass
    return out

def mse(v):
    n = len(v); m = sum(v) / n
    if n < 2: return m, float('nan')
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) / math.sqrt(n)

def load(arm):
    """per-replica means, keyed (system, replica)"""
    D = {}
    for d in sorted(glob.glob(f'{arm}/*_[0-9][0-9][0-9]')):
        sysn, rep = os.path.basename(d).rsplit('_', 1)
        f = f'{d}/per_frame.csv'
        if not os.path.exists(f): continue
        b = blocks(f)
        if not all(k in b and b[k] for k in ('complex', 'receptor', 'ligand')): continue
        D[(sysn, int(rep))] = {k: sum(b[k]) / len(b[k]) for k in ('complex', 'receptor', 'ligand', 'delta') if b.get(k)}
    return D

def ddh_own(D):
    """ST: each complex's own bound MOB1A (the DELTA block)."""
    per = {s: mse([v['delta'] for (ss, r), v in D.items() if ss == s]) for s in ('lats1', 'ndr2')}
    return per['lats1'][0] - per['ndr2'][0], math.hypot(per['lats1'][1], per['ndr2'][1]), len(D)

def ddh_pooled(D):
    ref = [v['receptor'] for v in D.values()]
    rm, _ = mse(ref)
    per = {}
    for s in ('lats1', 'ndr2'):
        vals = [v['complex'] - v['ligand'] - rm for (ss, r), v in D.items() if ss == s]
        if not vals: return None
        per[s] = mse(vals)
    return per['lats1'][0] - per['ndr2'][0], math.hypot(per['lats1'][1], per['ndr2'][1]), len(ref)

root = sys.argv[1] if len(sys.argv) > 1 else 'pb_sensitivity/opc-salt'
for mode in ('st', 'mt'):
    arms = {}
    for d in sorted(glob.glob(f'{root}/{mode}/istrng_*')):
        D = load(d)
        if D: arms[os.path.basename(d).replace('istrng_', '')] = D
    if len(arms) < 2:
        print(f'{mode.upper()}: {len(arms)} arm(s) with data -- not ready\n'); continue
    lo, hi = sorted(arms)[0], sorted(arms)[-1]
    print('=' * 80)
    ddh = ddh_own if mode == 'st' else ddh_pooled
    label = ("receptor+ligand from the COMPLEX frames; each complex's own MOB1A (DELTA block)"
             if mode == 'st' else
             'pooled MOB1A reference; HISTORICAL: unbound runs from before the matched-salt correction, not used in the paper')
    print(f'{mode.upper()}\n  {label}')
    print('=' * 80)
    for k in sorted(arms):
        r = ddh(arms[k])
        print(f'  istrng={k}   ddH = {r[0]:+7.2f} +- {r[1]:5.2f} kcal/mol   ({r[2]} replicas)')
    common = sorted(set(arms[lo]) & set(arms[hi]))
    pair = {}
    for s in ('lats1', 'ndr2'):
        if mode == 'st':
            pair[s] = mse([arms[hi][k]['delta'] - arms[lo][k]['delta'] for k in common if k[0] == s])
        else:
            pair[s] = mse([(arms[hi][k]['complex'] - arms[hi][k]['ligand'])
                           - (arms[lo][k]['complex'] - arms[lo][k]['ligand'])
                           for k in common if k[0] == s])
    dref = mse([arms[hi][k]['receptor'] - arms[lo][k]['receptor'] for k in common])
    d = pair['lats1'][0] - pair['ndr2'][0]; se = math.hypot(pair['lats1'][1], pair['ndr2'][1])
    print(f'\n  PAIRED {lo} -> {hi}  (same frames, so frame noise cancels)')
    print(f'    d(ddH)      = {d:+.3f} +- {se:.3f} kcal/mol'
          f'   [{abs(d)/se:.1f} sigma]' if se == se and se > 0 else '')
    print(f'    LATS1 shift = {pair["lats1"][0]:+.2f} +- {pair["lats1"][1]:.2f}')
    print(f'    NDR2  shift = {pair["ndr2"][0]:+.2f} +- {pair["ndr2"][1]:.2f}')
    if mode == 'mt':
        print(f'    shared MOB1A reference shift = {dref[0]:+.2f} +- {dref[1]:.2f}  (cancels in ddH)')
    for k in sorted(arms):
        r = ddh(arms[k])
        print(f'    sign at istrng={k}: {"POSITIVE (experimental order)" if r[0] > 0 else "NEGATIVE"}')
    print()
