#!/usr/bin/env python3
"""Stereochemistry check for MOB1A (residues 1-181) in one or more PDB files.

For each file prints: residues with D-configured C-alpha (improper N-CA-C-CB > 0; L ~ -120 deg, D ~ +120 deg),
non-proline cis peptide bonds (|omega| < 30 deg), heavy-atom contacts < 2.5 A involving loop residues 60-78,
and the improper angle of loop residues 64-74 (simulation numbering; MOB1A = +32).
Usage: stereo_check.py label1=file1.pdb label2=file2.pdb ...
"""
import sys
import numpy as np

def load(path):
    R = {}
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            r = int(l[22:26])
            if r > 181:
                continue
            R.setdefault(r, {"name": l[17:20]})[l[12:16].strip()] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return R

def dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))

def check(path):
    R = load(path)
    rs = sorted(R)
    imp = {r: dih(R[r]["N"], R[r]["CA"], R[r]["C"], R[r]["CB"]) for r in rs
           if all(k in R[r] for k in ("N", "CA", "C", "CB"))}
    D = [r for r, v in imp.items() if v > 0]
    cis = [f"{r}-{r+1}" for r in rs if r + 1 in R and "C" in R[r] and "CA" in R[r] and "N" in R[r + 1]
           and "CA" in R[r + 1] and R[r + 1]["name"] != "PRO"
           and abs(dih(R[r]["CA"], R[r]["C"], R[r + 1]["N"], R[r + 1]["CA"])) < 30]
    H = [(r, x) for r in rs for k, x in R[r].items() if k != "name" and not k.startswith("H") and k != "ZN"
         and not (k[0].isdigit() and k[1] == "H")]
    X = np.array([h[1] for h in H]); rr = np.array([h[0] for h in H])
    loop = (rr >= 60) & (rr <= 78)
    Dm = np.sqrt(((X[loop][:, None] - X[None]) ** 2).sum(-1))
    far = np.abs(rr[loop][:, None] - rr[None]) >= 2
    clashes = int(((Dm < 2.5) & far).sum() // 1)
    return D, cis, clashes, (Dm[far].min() if far.any() else np.nan), imp

if __name__ == "__main__":
    rows = []
    for arg in sys.argv[1:]:
        label, path = arg.split("=", 1)
        D, cis, ncl, dmin, imp = check(path)
        loop = " ".join(f"{imp[r]:+.0f}" if r in imp else "  G " for r in range(64, 75))
        print(f"{label:34s} D: {str(D or '-'):<18s} cis: {str(cis or '-'):<34s} clashes<2.5: {ncl:3d} (min {dmin:.2f})  | {loop}")
