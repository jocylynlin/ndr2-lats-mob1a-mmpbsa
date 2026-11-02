#!/usr/bin/env python3
"""verify_strip.py raw.nc stripped.nc raw.prmtop expected_frames [solute_range]

Gate for deleting a raw trajectory. Every check must pass:
  1. the stripped file has exactly expected_frames frames
  2. its atom count equals the solute atom count implied by solute_range (when given)
  3. every frame is readable (full traversal) and contains no NaN or infinite coordinate
  4. coordinates agree with the raw file at the first, last and several deterministic intermediate frames,
     compared per molecule after removing each molecule's centroid, because autoimage shifts whole molecules
  5. the stripped file has no box information (it was written with `nobox`) and the raw file does
Exit status 0 only if all pass.
"""
import re, sys
import numpy as np
from pathlib import Path
from scipy.io import netcdf_file

raw, strip, prm, want = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4])
rng = sys.argv[5] if len(sys.argv) > 5 else None
t = prm.read_text(errors="replace")
F = {m.group(1): m.group(3) for m in re.finditer(r"%FLAG (\S+)\s*\n%FORMAT\(([^)]+)\)\s*\n(.*?)(?=%FLAG|\Z)", t, re.S)}
apm = np.array(F["ATOMS_PER_MOLECULE"].split(), dtype=int)
rptr = np.array(F["RESIDUE_POINTER"].split(), dtype=int)
natom_raw = int(F["POINTERS"].split()[0])
fail = []

def expected_solute_atoms(rng):
    if not rng: return None
    last = int(rng.split("-")[1])
    return (rptr[last] - 1) if last < len(rptr) else natom_raw

fs = netcdf_file(str(strip), "r", mmap=True)
cs = fs.variables["coordinates"]
ns, na_s = cs.shape[0], cs.shape[1]
has_box_s = "cell_lengths" in fs.variables
if ns != want: fail.append(f"frame count {ns} != {want}")
exp = expected_solute_atoms(rng)
if exp is not None and na_s != exp: fail.append(f"atom count {na_s} != expected solute {exp}")
if has_box_s: fail.append("stripped file still carries box information")

# full traversal, in blocks, checking every frame is readable and finite
bad = 0
for i0 in range(0, ns, 250):
    blk = np.array(cs[i0:i0+250], dtype=np.float32)
    if not np.isfinite(blk).all(): bad += int((~np.isfinite(blk)).any(axis=(1, 2)).sum())
if bad: fail.append(f"{bad} frames contain NaN or infinite coordinates")
fs.close()

fr = netcdf_file(str(raw), "r", mmap=True)
cr = fr.variables["coordinates"]
nr = cr.shape[0]
if nr < want: fail.append(f"raw file has {nr} frames, fewer than {want}")
if "cell_lengths" not in fr.variables: fail.append("raw file has no box information (unexpected)")
usable = min(ns, nr, want)
checks = sorted({0, want // 4, want // 2, (3 * want) // 4, want - 1} & set(range(usable)))
edges = np.cumsum(apm); bounds = np.concatenate([[0], edges[edges <= na_s]])
dev = 0.0
fs = netcdf_file(str(strip), "r", mmap=True); cs = fs.variables["coordinates"]
for f in checks:
    rr = np.array(cr[f], float)[:na_s]; ss = np.array(cs[f], float)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        dev = max(dev, float(np.abs((rr[lo:hi] - rr[lo:hi].mean(0)) - (ss[lo:hi] - ss[lo:hi].mean(0))).max()))
fs.close(); fr.close()
if dev >= 0.01: fail.append(f"coordinate deviation {dev:.5f} A at frames {checks}")

status = "OK" if not fail else "FAIL"
print(f"{strip.name}: frames {ns}/{want}, atoms {na_s}"
      f"{'' if exp is None else '/'+str(exp)}, checked frames {checks}, max deviation {dev:.5f} A -> {status}"
      + ("" if not fail else "  [" + "; ".join(fail) + "]"))
sys.exit(0 if not fail else 1)
