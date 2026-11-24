#!/usr/bin/env python3
"""Energetic origin of the water-model dependence of MT ddH_bind (LATS1 minus NDR2).

Reads ONLY the corrected matched-salt MM-PBSA results:
    campaign_saltmatched/mmpbsa/{opc-salt-matched,tip3p-matched}/{st,mt}/{lats1,ndr2}_NNN/per_frame.csv

Conventions (match manuscript_rev/template.tex Tables 1-3):
  groups   Covalent = BOND+ANGLE+DIHED+UB+IMP+CMAP+1-4 VDW+1-4 EEL
           Electrostatic = EEL+EPB ; Dispersion = VDWAALS+EDISPER ; Hydrophobic = ENPOLAR
  MT       dH_n = C_n - K_n - Rbar      (Rbar = pooled unbound MOB1A, 200 replicas)
  ST       dH_n = C_n - R_n - K_n       (all from the complex trajectory)
  SE       replica-level: s_replica/sqrt(N), N=100
  ddH      (C_L-K_L) - (C_N-K_N): Rbar cancels EXACTLY and is excluded from its SE

Bootstrap: resamples REPLICA INDICES, never frames. Each (condition, system) is an
independent simulation set and is resampled independently. For reference-free quantities
(any LATS1-minus-NDR2 contrast, and the OPC-minus-TIP3P contrast of those) Rbar is not
resampled, because it cancels algebraically. For single-system absolute quantities, which
do contain Rbar, the 200 pooled reference replicas ARE resampled.
"""
import csv, glob, json, math, os, random, sys
from collections import defaultdict

ROOT = '/home3/jocelyn/projects/hippo/campaign_saltmatched'
OUT  = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
CONDS = {'OPC': 'opc-salt-matched', 'TIP3P': 'tip3p-matched'}
SYSTEMS = ['lats1', 'ndr2']
GROUPS = {
    'Covalent':      ['BOND','ANGLE','DIHED','UB','IMP','CMAP','1-4 VDW','1-4 EEL'],
    'Electrostatic': ['EEL','EPB'],
    'Dispersion':    ['VDWAALS','EDISPER'],
    'Hydrophobic':   ['ENPOLAR'],
}
GORDER = ['Covalent','Electrostatic','Dispersion','Hydrophobic']
BLOCKS = ['complex','receptor','ligand','delta']
NBOOT = 10000
SEED = 20261122

# ---------------------------------------------------------------- parsing
def read_per_frame(path):
    """-> {block: {col: mean_over_frames}}, plus frame count. Raw MMPBSA.py -eo layout."""
    out, cur, idx, hdr = {}, None, -1, None
    nfr = {}
    acc = defaultdict(lambda: defaultdict(float))
    cnt = defaultdict(int)
    for row in csv.reader(open(path)):
        if not row:
            continue
        if row[0].startswith('Frame #'):
            hdr = row[1:]; idx += 1
            cur = BLOCKS[idx] if idx < len(BLOCKS) else None
            continue
        if cur is None or len(row) < 10:
            continue
        try:
            vals = [float(x) for x in row[1:]]
        except (ValueError, IndexError):
            continue
        for h, v in zip(hdr, vals):
            acc[cur][h] += v
        cnt[cur] += 1
    for b in cnt:
        nfr[b] = cnt[b]
        out[b] = {h: acc[b][h] / cnt[b] for h in acc[b]}
    return out, nfr

def grouped(block):
    """Collapse raw columns into the four manuscript groups."""
    g = {}
    for name, cols in GROUPS.items():
        g[name] = sum(block.get(c, 0.0) for c in cols)
    g['TOTAL'] = block.get('TOTAL', block.get('DELTA TOTAL', float('nan')))
    for c in ('BOND','ANGLE','DIHED','1-4 VDW','1-4 EEL'):
        g[c] = block.get(c, 0.0)
    return g

def load(cond_dir, mode):
    """-> per[(system, replica)] = {block: grouped dict}"""
    per = {}
    pat = f'{ROOT}/mmpbsa/{cond_dir}/{mode}/*_[0-9][0-9][0-9]'
    for d in sorted(glob.glob(pat)):
        base = os.path.basename(d)
        sysn, rep = base.rsplit('_', 1)
        f = os.path.join(d, 'per_frame.csv')
        if not os.path.exists(f):
            continue
        blocks, nfr = read_per_frame(f)
        if not all(b in blocks for b in ('complex','receptor','ligand')):
            continue
        rec = {b: grouped(blocks[b]) for b in blocks}
        rec['_nframes'] = nfr
        per[(sysn, int(rep))] = rec
    return per

# ---------------------------------------------------------------- stats
def mean(v):
    return sum(v) / len(v)

def se(v):
    n = len(v)
    if n < 2:
        return float('nan')
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) / math.sqrt(n)

def ci95(samples):
    s = sorted(samples)
    n = len(s)
    return s[int(0.025 * n)], s[int(0.975 * n)]

# ---------------------------------------------------------------- Step 1: validation
def validate(data):
    """Returns (list_of_report_lines, ok_flag)."""
    L, ok = [], True
    L.append('## Step 1 validation\n')

    # 1.2 replica counts and frame counts
    L.append('### 1.2 Replica and frame counts')
    L.append('| condition | mode | system | replicas | frames/replica |')
    L.append('|---|---|---|---|---|')
    for wm, cd in CONDS.items():
        for mode in ('st','mt'):
            for s in SYSTEMS:
                reps = [r for (ss, r) in data[(wm, mode)] if ss == s]
                nf = {tuple(sorted(v['_nframes'].items()))
                      for (ss, r), v in data[(wm, mode)].items() if ss == s}
                fr = sorted({n for t in nf for _, n in t})
                L.append(f'| {wm} | {mode.upper()} | {s} | {len(reps)} | {fr} |')
                if len(reps) != 100:
                    L.append(f'> **DISCREPANCY**: expected 100 replicas, found {len(reps)}')
                    ok = False
    L.append('')

    # 1.3 reconstruct DELTA TOTAL from raw terms
    L.append('### 1.3 Reconstruction of DELTA TOTAL from the four groups')
    L.append('ST files carry an explicit DELTA block, so the four groups can be checked')
    L.append('against a value MMPBSA.py wrote itself. MT files carry only three blocks,')
    L.append('so there the check is that the groups sum to the TOTAL-column arithmetic.\n')
    L.append('| condition | mode | max abs residual (kcal/mol) | basis |')
    L.append('|---|---|---|---|')
    for wm, cd in CONDS.items():
        for mode in ('st','mt'):
            per = data[(wm, mode)]
            worst = 0.0
            if mode == 'st':
                for (s, r), v in per.items():
                    if 'delta' not in v:
                        continue
                    lhs = sum(v['delta'][g] for g in GORDER)
                    rhs = v['delta']['TOTAL']
                    worst = max(worst, abs(lhs - rhs))
                basis = 'written DELTA block'
            else:
                pool = mean([v['receptor']['TOTAL'] for v in per.values()])
                poolg = {g: mean([v['receptor'][g] for v in per.values()]) for g in GORDER}
                for (s, r), v in per.items():
                    lhs = sum(v['complex'][g] - v['ligand'][g] - poolg[g] for g in GORDER)
                    rhs = v['complex']['TOTAL'] - v['ligand']['TOTAL'] - pool
                    worst = max(worst, abs(lhs - rhs))
                basis = 'TOTAL-column arithmetic'
            flag = '' if worst < 5e-3 else '  <-- **DISCREPANCY**'
            if worst >= 5e-3:
                ok = False
            L.append(f'| {wm} | {mode.upper()} | {worst:.2e} | {basis}{flag} |')
    L.append('')

    # 1.4 algebraic cancellation of the pooled reference
    L.append('### 1.4 Pooled MOB1A reference cancels in LATS1-minus-NDR2')
    L.append('ddH computed two ways: (a) each system against the pooled reference, then')
    L.append('differenced; (b) reference omitted entirely. These must be bit-identical.\n')
    L.append('| condition | with Rbar | without Rbar | difference |')
    L.append('|---|---|---|---|')
    for wm, cd in CONDS.items():
        per = data[(wm, 'mt')]
        poolg = {g: mean([v['receptor'][g] for v in per.values()]) for g in GORDER}
        a, b = {}, {}
        for s in SYSTEMS:
            vs = [v for (ss, r), v in per.items() if ss == s]
            a[s] = mean([sum(v['complex'][g] - v['ligand'][g] - poolg[g] for g in GORDER) for v in vs])
            b[s] = mean([sum(v['complex'][g] - v['ligand'][g] for g in GORDER) for v in vs])
        da, db = a['lats1'] - a['ndr2'], b['lats1'] - b['ndr2']
        resid = abs(da - db)
        flag = '' if resid < 1e-9 else '  <-- **DISCREPANCY**'
        if resid >= 1e-9:
            ok = False
        L.append(f'| {wm} | {da:+.10f} | {db:+.10f} | {resid:.2e}{flag} |')
    L.append('')
    return L, ok

# ---------------------------------------------------------------- Step 2: quantities
def series(per, mode, group):
    """Per-replica series for each system, plus the pooled reference series.
    Returns (dict system -> {'C':[],'K':[],'dH_nr':[]}, reference_list)
    dH_nr = C - K  (reference-free part of MT dH); for ST, C - R - K."""
    ref = [v['receptor'][group] for v in per.values()]
    outd = {}
    for s in SYSTEMS:
        reps = sorted(r for (ss, r) in per if ss == s)
        C = [per[(s, r)]['complex'][group] for r in reps]
        K = [per[(s, r)]['ligand'][group] for r in reps]
        if mode == 'st':
            R = [per[(s, r)]['receptor'][group] for r in reps]
            dh = [c - rr - k for c, rr, k in zip(C, R, K)]
        else:
            dh = [c - k for c, k in zip(C, K)]
        outd[s] = {'C': C, 'K': K, 'dH_nr': dh}
    return outd, ref

def main():
    random.seed(SEED)
    data = {}
    for wm, cd in CONDS.items():
        for mode in ('st','mt'):
            data[(wm, mode)] = load(cd, mode)

    vlines, vok = validate(data)

    rows = []           # CSV rows
    boot_store = {}     # for reuse in the report

    # ---- build per-replica series per (condition, mode, group)
    S = {}
    for wm in CONDS:
        for mode in ('st','mt'):
            for g in GORDER:
                S[(wm, mode, g)] = series(data[(wm, mode)], mode, g)

    def dH_total(wm, mode, s):
        """per-replica reference-free dH summed over groups"""
        n = len(S[(wm, mode, GORDER[0])][0][s]['dH_nr'])
        return [sum(S[(wm, mode, g)][0][s]['dH_nr'][i] for g in GORDER) for i in range(n)]

    def ref_total(wm, mode):
        n = len(S[(wm, mode, GORDER[0])][1])
        return [sum(S[(wm, mode, g)][1][i] for g in GORDER) for i in range(n)]

    # ================= Panel A / Step 2: ddH per condition, per group and total
    def ddh_stats(wm, mode, group=None):
        if group is None:
            L, N = dH_total(wm, mode, 'lats1'), dH_total(wm, mode, 'ndr2')
        else:
            L = S[(wm, mode, group)][0]['lats1']['dH_nr']
            N = S[(wm, mode, group)][0]['ndr2']['dH_nr']
        d = mean(L) - mean(N)
        s = math.hypot(se(L), se(N))
        bs = []
        for _ in range(NBOOT):
            a = mean([random.choice(L) for _ in L])
            b = mean([random.choice(N) for _ in N])
            bs.append(a - b)
        return d, s, ci95(bs), bs

    for mode in ('mt','st'):
        for wm in CONDS:
            d, s, ci, bs = ddh_stats(wm, mode, None)
            boot_store[('ddh', wm, mode, 'TOTAL')] = bs
            rows.append(dict(panel='A', quantity='ddH_bind', mode=mode.upper(), condition=wm,
                             group='TOTAL', system='LATS1-NDR2', mean=d, se=s,
                             ci_lo=ci[0], ci_hi=ci[1], sigma=abs(d)/s if s else float('nan')))
            for g in GORDER:
                d, s, ci, bs = ddh_stats(wm, mode, g)
                boot_store[('ddh', wm, mode, g)] = bs
                rows.append(dict(panel='A', quantity='ddH_bind', mode=mode.upper(), condition=wm,
                                 group=g, system='LATS1-NDR2', mean=d, se=s,
                                 ci_lo=ci[0], ci_hi=ci[1], sigma=abs(d)/s if s else float('nan')))

    # ================= Panel B / Step 2: OPC-minus-TIP3P of the ddH contrast
    for mode in ('mt','st'):
        for g in ['TOTAL'] + GORDER:
            bo = boot_store[('ddh','OPC',mode,g)]
            bt = boot_store[('ddh','TIP3P',mode,g)]
            do, so, _, _ = ddh_stats('OPC', mode, None if g=='TOTAL' else g)
            dt, st_, _, _ = ddh_stats('TIP3P', mode, None if g=='TOTAL' else g)
            d = do - dt
            s = math.hypot(so, st_)
            bs = [x - y for x, y in zip(bo, bt)]     # independent sim sets
            ci = ci95(bs)
            rows.append(dict(panel='B', quantity='watermodel_effect_on_contrast', mode=mode.upper(),
                             condition='OPC-TIP3P', group=g, system='LATS1-NDR2', mean=d, se=s,
                             ci_lo=ci[0], ci_hi=ci[1], sigma=abs(d)/s if s else float('nan')))

    # ================= Panel C: per-system OPC-minus-TIP3P, per group
    # single-system dH DOES contain Rbar -> resample the 200 reference replicas too
    def sys_dH_stats(wm, mode, s, group=None):
        if group is None:
            X, R = dH_total(wm, mode, s), ref_total(wm, mode)
        else:
            X = S[(wm, mode, group)][0][s]['dH_nr']
            R = S[(wm, mode, group)][1]
        if mode == 'st':
            m, sd = mean(X), se(X)
            bs = [mean([random.choice(X) for _ in X]) for _ in range(NBOOT)]
            return m, sd, bs
        m = mean(X) - mean(R)
        sd = math.hypot(se(X), se(R))
        bs = [mean([random.choice(X) for _ in X]) - mean([random.choice(R) for _ in R])
              for _ in range(NBOOT)]
        return m, sd, bs

    for mode in ('mt',):
        for s in SYSTEMS:
            for g in ['TOTAL'] + GORDER:
                gg = None if g == 'TOTAL' else g
                mo, so, bo = sys_dH_stats('OPC', mode, s, gg)
                mt_, st_, bt = sys_dH_stats('TIP3P', mode, s, gg)
                for wm, m, sd in (('OPC', mo, so), ('TIP3P', mt_, st_)):
                    rows.append(dict(panel='C', quantity='dH_bind', mode=mode.upper(),
                                     condition=wm, group=g, system=s.upper(), mean=m, se=sd,
                                     ci_lo=float('nan'), ci_hi=float('nan'),
                                     sigma=float('nan')))
                d = mo - mt_
                sd = math.hypot(so, st_)
                bs = [x - y for x, y in zip(bo, bt)]
                ci = ci95(bs)
                rows.append(dict(panel='C', quantity='dH_bind_OPC_minus_TIP3P', mode=mode.upper(),
                                 condition='OPC-TIP3P', group=g, system=s.upper(), mean=d, se=sd,
                                 ci_lo=ci[0], ci_hi=ci[1], sigma=abs(d)/sd if sd else float('nan')))

    # ================= Panel D: bound vs unbound decomposition
    # ddH = (C_L - C_N) - (K_L - K_N).  Rbar cancels, so it does not appear at all.
    # water-model effect = [dC(OPC)-dC(TIP3P)] - [dK(OPC)-dK(TIP3P)]
    def species_contrast(wm, group, which):
        """which in {'C','K'}: LATS1-minus-NDR2 of that species block."""
        if group is None:
            Lv = [sum(S[(wm,'mt',g)][0]['lats1'][which][i] for g in GORDER)
                  for i in range(len(S[(wm,'mt',GORDER[0])][0]['lats1'][which]))]
            Nv = [sum(S[(wm,'mt',g)][0]['ndr2'][which][i] for g in GORDER)
                  for i in range(len(S[(wm,'mt',GORDER[0])][0]['ndr2'][which]))]
        else:
            Lv = S[(wm,'mt',group)][0]['lats1'][which]
            Nv = S[(wm,'mt',group)][0]['ndr2'][which]
        d = mean(Lv) - mean(Nv)
        s = math.hypot(se(Lv), se(Nv))
        bs = [mean([random.choice(Lv) for _ in Lv]) - mean([random.choice(Nv) for _ in Nv])
              for _ in range(NBOOT)]
        return d, s, bs

    for g in ['TOTAL'] + GORDER:
        gg = None if g == 'TOTAL' else g
        store = {}
        for which, sign, label in (('C', +1, 'bound complex'), ('K', -1, 'unbound kinase')):
            dco, sco, bco = species_contrast('OPC', gg, which)
            dct, sct, bct = species_contrast('TIP3P', gg, which)
            # contribution to the OPC-minus-TIP3P water-model effect, with algebraic sign
            contrib = sign * (dco - dct)
            scon = math.hypot(sco, sct)
            bs = [sign * (x - y) for x, y in zip(bco, bct)]
            ci = ci95(bs)
            store[which] = (contrib, scon, bs)
            rows.append(dict(panel='D', quantity=f'watermodel_contrib_{label.replace(" ","_")}',
                             mode='MT', condition='OPC-TIP3P', group=g,
                             system=f'sign={sign:+d}', mean=contrib, se=scon,
                             ci_lo=ci[0], ci_hi=ci[1],
                             sigma=abs(contrib)/scon if scon else float('nan')))
        # sum check: bound + unbound must equal the Panel B water-model effect
        tot = store['C'][0] + store['K'][0]
        bs = [a + b for a, b in zip(store['C'][2], store['K'][2])]
        ci = ci95(bs)
        rows.append(dict(panel='D', quantity='watermodel_contrib_SUM_CHECK', mode='MT',
                         condition='OPC-TIP3P', group=g, system='bound+unbound',
                         mean=tot, se=math.hypot(store['C'][1], store['K'][1]),
                         ci_lo=ci[0], ci_hi=ci[1], sigma=float('nan')))

    # ---- pooled reference absolute values, per condition (context for Step 2)
    for wm in CONDS:
        for g in ['TOTAL'] + GORDER:
            R = ref_total(wm,'mt') if g=='TOTAL' else S[(wm,'mt',g)][1]
            rows.append(dict(panel='ref', quantity='pooled_unbound_MOB1A', mode='MT',
                             condition=wm, group=g, system='pooled_200', mean=mean(R), se=se(R),
                             ci_lo=float('nan'), ci_hi=float('nan'), sigma=float('nan')))

    # ---- write CSV
    os.makedirs(OUT, exist_ok=True)
    cols = ['panel','quantity','mode','condition','group','system','mean','se','ci_lo','ci_hi','sigma']
    with open(f'{OUT}/watermodel_origin.csv','w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{r[k]:.6f}' if isinstance(r[k], float) else r[k]) for k in cols})

    with open(f'{OUT}/validation.md','w') as fh:
        fh.write('\n'.join(vlines))

    print('\n'.join(vlines))
    print(f'\nwrote {OUT}/watermodel_origin.csv  ({len(rows)} rows)')
    print(f'validation: {"ALL CHECKS PASS" if vok else "DISCREPANCIES FOUND -- see validation.md"}')
    return 0 if vok else 1

if __name__ == '__main__':
    sys.exit(main())
