#!/usr/bin/env python3
"""check_topologies.py complex.prmtop receptor.prmtop ligand.prmtop

Verifies the MM-PBSA topology split, rather than trusting the residue numbers:
  1. the complex splits into exactly two solute molecules; molecule 1 is the receptor (MOB1A + Zn), molecule 2
     the kinase. The receptor/ligand residue ranges are DERIVED from that boundary.
  2. the zinc is a single residue and sits at the end of the receptor, bonded into molecule 1
  3. receptor and ligand atom counts sum to the complex solute atom count
  4. every receptor and ligand atom matches the corresponding complex atom in name, type, mass, charge and
     Lennard-Jones parameters. parmstrip cuts from the complex, so these must be identical by construction;
     a difference means the wrong topology was used.
Exit status 0 only if all pass.
"""
import re, sys
import numpy as np
from pathlib import Path

def P(p):
    t = Path(p).read_text(errors="replace")
    F = {m.group(1): m.group(3) for m in re.finditer(r"%FLAG (\S+)\s*\n%FORMAT\(([^)]+)\)\s*\n(.*?)(?=%FLAG|\Z)", t, re.S)}
    def s4(k, w=4):
        x = F[k].replace("\n", ""); return [x[i:i+w].strip() for i in range(0, len(x), w)]
    g = lambda k, ty=float: np.array(F[k].split(), ty)
    ptr = g("POINTERS", int); nt = ptr[1]; ti = g("ATOM_TYPE_INDEX", int)
    ico = g("NONBONDED_PARM_INDEX", int)[nt*(ti-1)+(ti-1)] - 1
    A, B = g("LENNARD_JONES_ACOEF")[ico], g("LENNARD_JONES_BCOEF")[ico]
    return dict(nat=ptr[0], name=s4("ATOM_NAME"), type=s4("AMBER_ATOM_TYPE"), res=s4("RESIDUE_LABEL"),
                rptr=g("RESIDUE_POINTER", int),
                bonds=np.concatenate([g("BONDS_INC_HYDROGEN", int), g("BONDS_WITHOUT_HYDROGEN", int)]).reshape(-1, 3),
                mass=g("MASS"), chg=g("CHARGE")/18.2223, A=A, B=B)

argv = [a for a in sys.argv[1:] if a != "--ranges"]
ranges_only = "--ranges" in sys.argv
cpx = P(argv[0])
fail = []
def molecules(d):
    """connected components of the bond graph, as sorted atom-index arrays (works on stripped topologies,
    which carry no ATOMS_PER_MOLECULE flag)"""
    n = d["nat"]; parent = list(range(n))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i, j, _ in d["bonds"]:
        a, b = find(abs(i)//3), find(abs(j)//3)
        if a != b: parent[a] = b
    groups = {}
    for a in range(n): groups.setdefault(find(a), []).append(a)
    return sorted((np.array(v) for v in groups.values()), key=lambda v: -len(v))
comps = molecules(cpx)
mols = [len(c) for c in comps if len(c) > 3]                 # solute molecules (ions/waters are 1-3 atoms)
if len(mols) >= 2 and not (comps[0].min() == 0 and comps[1].min() == comps[0].max() + 1):
    fail.append("the two largest molecules are not contiguous blocks starting at atom 1")
if len(mols) != 2: fail.append(f"complex has {len(mols)} solute molecules, expected 2 (sizes {mols[:4]})")
n_rec = int(mols[0]) if mols else 0
n_lig = int(mols[1]) if len(mols) > 1 else 0
if ranges_only:
    n_rec_ = int(mols[0]); last_ = int(np.searchsorted(cpx["rptr"] - 1, n_rec_))
    print(last_, n_rec_, int(mols[1])); sys.exit(0 if not fail else 1)
rec, lig = P(argv[1]), P(argv[2])
zn = [i+1 for i, r in enumerate(cpx["res"]) if r.startswith("ZN")]
if len(zn) != 1: fail.append(f"expected exactly one ZN residue, found {zn}")
elif cpx["rptr"][zn[0]-1] - 1 >= n_rec: fail.append(f"ZN residue {zn[0]} lies outside molecule 1 (the receptor)")
last_rec_res = int(np.searchsorted(cpx["rptr"] - 1, n_rec))   # residues 1..last_rec_res are the receptor
if zn and zn[0] != last_rec_res: fail.append(f"ZN is residue {zn[0]} but the receptor ends at residue {last_rec_res}")
if rec["nat"] != n_rec: fail.append(f"receptor has {rec['nat']} atoms, molecule 1 has {n_rec}")
if lig["nat"] != n_lig: fail.append(f"ligand has {lig['nat']} atoms, molecule 2 has {n_lig}")
if rec["nat"] + lig["nat"] != n_rec + n_lig: fail.append("receptor + ligand != complex solute atoms")

# parameters must be identical to the corresponding complex atoms
for lab, sub, off in (("receptor", rec, 0), ("ligand", lig, n_rec)):
    n = sub["nat"]
    if off + n > cpx["nat"]: fail.append(f"{lab}: range exceeds complex"); continue
    sl = slice(off, off+n)
    if cpx["name"][sl] != sub["name"][:n]: fail.append(f"{lab}: atom names differ from complex")
    if cpx["type"][sl] != sub["type"][:n]: fail.append(f"{lab}: atom types differ from complex")
    for key, tol in (("mass", 1e-9), ("chg", 1e-9), ("A", 1e-6), ("B", 1e-6)):
        d = float(np.abs(cpx[key][sl] - sub[key][:n]).max())
        if d > tol: fail.append(f"{lab}: {key} differs from complex by up to {d:.3e}")

print(f"complex {cpx['nat']} atoms -> receptor {rec['nat']} (residues 1-{last_rec_res}, ends at {cpx['res'][last_rec_res-1]}) "
      f"+ ligand {lig['nat']} (residues {last_rec_res+1}-{len(cpx['res'])}); "
      + ("all parameters identical to the complex -> OK" if not fail else "FAIL"))
for f in fail: print("   - " + f)
sys.exit(0 if not fail else 1)
