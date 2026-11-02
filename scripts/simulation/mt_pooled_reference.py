#!/usr/bin/env python3
"""mt_pooled_reference.py [window]   window = "" (0-40 ns) or "_late" (40-100 ns)

Writes a POOLED-MOB1A-reference MT table beside the own-reference one.

The MT protocol as executed references each complex to its OWN unbound MOB1A ensemble
(-yr <sys>_mob1a_<n>.nc). Because MOB1A is the same molecule in both systems (181 residues,
2573 atoms, identical residue labels and atom names -- verified), that choice injects a pure
preparation artefact: the MOB1A internal control measures it at +11.19 +/- 2.56 (0-40 ns) and
+7.71 +/- 2.55 (40-100 ns), which is 4-5x the binding contrast and large enough to reverse its
sign. Pooling the reference over both systems and all replicas removes it, and in the
LATS1-NDR2 contrast the reference term cancels exactly.

Every energy term is recomputed as  cpx_<term> - <pooled rec_<term>> - lig_<term>.
The own-reference table is left untouched; it is the diagnostic, this is the result.
"""
import sys
import pandas as pd

C = "/home3/jocelyn/projects/hippo/campaign_100x40ns"
win = sys.argv[1] if len(sys.argv) > 1 else ""
src = f"{C}/analysis/mmpbsa_perframe_opc-salt_mt{win}.csv"
dst = f"{C}/analysis/mmpbsa_perframe_opc-salt_mt{win}_pooledref.csv"

d = pd.read_csv(src)
terms = [c[4:] for c in d.columns if c.startswith("cpx_")]

# pooled reference: mean of each rec_ term over every system and replica in this window
pooled = {t: d[f"rec_{t}"].astype(float).mean() for t in terms if f"rec_{t}" in d.columns}

for t, ref in pooled.items():
    d[t] = d[f"cpx_{t}"].astype(float) - ref - d[f"lig_{t}"].astype(float)
    d[f"rec_{t}"] = ref                      # record the reference actually used
d["delta_method"] = "pooled_mob1a_reference"

d.to_csv(dst, index=False)
print(f"wrote {dst}")
print(f"  pooled reference, {len(pooled)} terms; rec_TOTAL = {pooled.get('TOTAL', float('nan')):.2f}")
g = d[d.replica <= 20].groupby(["system", "replica"])["TOTAL"].mean()
for s in sorted({i[0] for i in g.index}):
    v = g.loc[s]
    print(f"  {s}: {v.mean():+.2f} +/- {v.std(ddof=1)/len(v)**0.5:.2f}  (n={len(v)})")
