#!/usr/bin/env python3
"""Per-species absolute block energies, OPC vs TIP3P, to separate
'unbound kinases did not change' from 'both changed equally and cancelled'."""
import math, random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_watermodel_origin import load, CONDS, GORDER, SYSTEMS, mean, se, ci95, NBOOT, SEED

random.seed(SEED)
D = {wm: load(cd, 'mt') for wm, cd in CONDS.items()}

def blockseries(wm, sysn, blk):
    per = D[wm]
    if blk == 'ref':
        return [sum(v['receptor'][g] for g in GORDER) for v in per.values()]
    reps = sorted(r for (ss, r) in per if ss == sysn)
    key = {'complex': 'complex', 'kinase': 'ligand'}[blk]
    return [sum(per[(sysn, r)][key][g] for g in GORDER) for r in reps]

print("="*100)
print("ABSOLUTE block energies by species (sum of the four groups), MT")
print("="*100)
print(f"{'species/block':22s} {'OPC':>18s} {'TIP3P':>18s} {'OPC-TIP3P':>18s} {'95% CI':>20s}")
res = {}
items = [('lats1','complex'),('ndr2','complex'),('lats1','kinase'),('ndr2','kinase'),(None,'ref')]
for sysn, blk in items:
    o = blockseries('OPC', sysn, blk); t = blockseries('TIP3P', sysn, blk)
    mo, so = mean(o), se(o); mt_, st_ = mean(t), se(t)
    d = mo - mt_; sd = math.hypot(so, st_)
    bs = [mean([random.choice(o) for _ in o]) - mean([random.choice(t) for _ in t])
          for _ in range(NBOOT)]
    lo, hi = ci95(bs)
    lbl = f"{sysn or 'pooled'}_{blk}"
    res[(sysn, blk)] = (d, sd)
    print(f"{lbl:22s} {mo:11.2f}+-{so:4.2f} {mt_:11.2f}+-{st_:4.2f} "
          f"{d:+11.2f}+-{sd:4.2f} [{lo:+7.2f},{hi:+7.2f}]")

print()
print("="*100)
print("KEY QUESTION: do the two unbound kinases shift by the SAME amount between water models?")
print("="*100)
dl, sl = res[('lats1','kinase')]
dn, sn = res[('ndr2','kinase')]
diff = dl - dn; sdif = math.hypot(sl, sn)
print(f"  LATS1 kinase shift  {dl:+8.2f} +- {sl:.2f}")
print(f"  NDR2  kinase shift  {dn:+8.2f} +- {sn:.2f}")
print(f"  difference          {diff:+8.2f} +- {sdif:.2f}   ({abs(diff)/sdif:.1f} sigma)")
print(f"  -> {'SAME within error: common shift, cancels in the contrast' if abs(diff)/sdif < 2 else 'DIFFERENT: contributes to the contrast'}")
print()
cl, scl = res[('lats1','complex')]
cn, scn = res[('ndr2','complex')]
cdiff = cl - cn; scdif = math.hypot(scl, scn)
print(f"  LATS1 complex shift {cl:+8.2f} +- {scl:.2f}")
print(f"  NDR2  complex shift {cn:+8.2f} +- {scn:.2f}")
print(f"  difference          {cdiff:+8.2f} +- {scdif:.2f}   ({abs(cdiff)/scdif:.1f} sigma)")
print(f"  -> {'SAME within error' if abs(cdiff)/scdif < 2 else 'DIFFERENT: this is what drives the contrast'}")
