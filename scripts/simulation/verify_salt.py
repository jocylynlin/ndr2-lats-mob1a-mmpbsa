#!/usr/bin/env python3
"""Check every salt-matched build against the complex box it is meant to match.

Criterion: NaCl pairs per water molecule.  That is what fixes the salt
concentration in the aqueous phase, and unlike pairs-per-box-volume it is
unaffected by the box shrinking during NPT equilibration.
"""
import collections, glob, os, re, sys

TO = 4 / (3 * 3 ** 0.5)          # truncated-octahedron volume factor
# the complex boxes we are KEEPING (published campaign, 80 pairs each)
TARGET = {('lats1', 'opc-salt-matched'): (80, 13930), ('ndr2', 'opc-salt-matched'): (80, 13574),
          ('lats1', 'tip3p-matched'):    (80, 13913), ('ndr2', 'tip3p-matched'):    (80, 13593)}

def parse(prm):
    t = open(prm).read()
    rl = re.search(r'%FLAG RESIDUE_LABEL.*?%FORMAT\([^)]*\)\s*\n(.*?)%FLAG', t, re.S).group(1).split()
    h = collections.Counter(rl)
    a = float(re.search(r'%FLAG BOX_DIMENSIONS.*?%FORMAT\([^)]*\)\s*\n(.*?)(?:%FLAG|$)',
                        t, re.S).group(1).split()[1])
    return h['Na+'], h['Cl-'], h['WAT'], a

root = os.path.dirname(os.path.abspath(__file__))
print(f"{'condition':18s} {'species':13s} {'Na':>4s} {'Cl':>4s} {'pairs':>6s} {'waters':>7s} "
      f"{'pairs/water':>12s} {'target':>11s} {'dev':>7s}")
bad = 0
for prm in sorted(glob.glob(f'{root}/inputs/*/*/*.prmtop')):
    cond, sp = prm.split('/')[-3], prm.split('/')[-2]
    na, cl, wat, a = parse(prm)
    pairs = min(na, cl)
    sysn = 'lats1' if sp.startswith('lats1') else 'ndr2'
    tp, tw = TARGET[(sysn, cond)]
    r, tr = pairs / wat, tp / tw
    dev = (r / tr - 1) * 100
    flag = '' if abs(dev) <= 2.0 else '   <-- OFF'
    if abs(dev) > 2.0: bad += 1
    print(f"{cond:18s} {sp:13s} {na:4d} {cl:4d} {pairs:6d} {wat:7d} "
          f"{r*1e3:10.4f}e-3 {tr*1e3:9.4f}e-3 {dev:+6.1f}%{flag}")
print()
# MOB1A is POOLED across both systems, so its two boxes must also match EACH OTHER,
# not just their respective complexes. That is why the NDR2 side carries 40 pairs
# rather than the 41 that would track its own complex most closely.
mob = {}
for prm in sorted(glob.glob(f'{root}/inputs/*/*/*.prmtop')):
    cond, sp = prm.split('/')[-3], prm.split('/')[-2]
    if not sp.endswith('mob1a'): continue
    na, cl, wat, a = parse(prm)
    mob[(cond, sp)] = min(na, cl) / wat
for cond in sorted({c for c, _ in mob}):
    rl = mob[(cond, 'lats1_mob1a')]; rn = mob[(cond, 'ndr2_mob1a')]
    spread = abs(rl / rn - 1) * 100
    ok = spread <= 0.5
    if not ok: bad += 1
    print(f"  {cond:18s} MOB1A-to-MOB1A spread {spread:5.2f}%  "
          f"{'ok (poolable)' if ok else '<-- TOO WIDE TO POOL'}")
print()
print(f"  {'all checks pass' if not bad else str(bad)+' check(s) failed'}")
sys.exit(1 if bad else 0)
