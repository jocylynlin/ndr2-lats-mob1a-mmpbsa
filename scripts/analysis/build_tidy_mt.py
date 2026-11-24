#!/usr/bin/env python3
"""build_tidy_mt.py <mmpbsa_root> <condition_dir> <condition_label> <water_model> <out.csv>

Builds the tidy multi-trajectory table energies/mmpbsa_perframe/mmpbsa_perframe_<condition>_mt.csv
from the per-replica MMPBSA.py outputs (<mmpbsa_root>/<condition_dir>/mt/<system>_<NNN>/per_frame.csv).

One row per frame per replica. cpx_*, rec_* and lig_* are the raw energies of the complex,
the unbound MOB1A (receptor) and the unbound kinase segment (ligand). The unprefixed terms are
the index-matched difference cpx - rec - lig, which pairs each complex replica with the
unbound replicas of the same index. The published multi-trajectory values do NOT use that
difference: they use the unbound MOB1A reference pooled over both systems, computed from the
rec_* columns by scripts/analysis/reproduce_tables.py.

Frames are 100 ps apart (every 10th saved frame); production starts 51.01 ns after the start
of heating, as in the single-trajectory tables.
"""
import csv, glob, os, sys

root, cdir, label, water, out = sys.argv[1:6]
SEC = {"Complex Energy Terms": "cpx", "Receptor Energy Terms": "rec", "Ligand Energy Terms": "lig"}

def parse(p):
    blocks, cur, hdr = {}, None, None
    for row in csv.reader(open(p)):
        if not row or not row[0].strip():
            continue
        key = row[0].strip()
        if key in SEC:
            cur = SEC[key]; blocks[cur] = {}; hdr = None; continue
        if key.startswith("DELTA"):
            cur = None; continue
        if key.startswith("Frame"):
            hdr = [c.strip() for c in row[1:]]; terms = hdr; continue
        if cur and hdr and key.isdigit():
            blocks[cur][int(key)] = [float(v) for v in row[1:]]
    return blocks, terms

rows, terms0 = [], None
for d in sorted(glob.glob(f"{root}/{cdir}/mt/*_[0-9][0-9][0-9]")):
    sysn, rep = os.path.basename(d).rsplit("_", 1)
    b, terms = parse(os.path.join(d, "per_frame.csv"))
    terms0 = terms0 or terms
    assert terms == terms0 and set(b) == {"cpx", "rec", "lig"}, d
    n = min(len(b[k]) for k in b)
    for i in range(n):
        c, r, l = b["cpx"][i], b["rec"][i], b["lig"][i]
        delta = [x - y - z for x, y, z in zip(c, r, l)]
        rows.append([sysn, water, label, "MT", "index_matched_cpx_minus_rec_minus_lig", int(rep), i,
                     1 + 10 * i, round(0.1 * i, 1), round(51.01 + 0.1 * i, 2)] + delta + c + r + l)

hdr = ["system", "water_model", "salt_condition", "protocol", "delta_method", "replica", "frame",
       "traj_frame", "time_ns", "time_ns_abs"] + terms0 + [f"{p}_{t}" for p in ("cpx", "rec", "lig") for t in terms0]
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
print(out, len(rows), "rows")
