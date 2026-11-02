#!/usr/bin/env python3
"""mmpbsa_tidy.py [condition] [mode]

Builds ONE tidy table with a row per frame per replica, from the preserved per_frame.csv files.
Nothing is aggregated - every MM-PBSA snapshot is kept, so any averaging or convergence analysis can be done
downstream. The original MM-PBSA outputs are left untouched.

Columns: system, water_model, salt_condition, protocol, replica, frame, traj_frame, time_ns, time_ns_abs,
then the DELTA (binding) terms unprefixed, and the same terms per species as cpx_*, rec_*, lig_*.
  frame       0-based index within the MM-PBSA calculation (0..399)
  traj_frame  1-based frame in the stripped trajectory (1, 11, 21, ...), i.e. startframe + frame*interval
  time_ns     production time of that frame, relative to the first analysed frame
  time_ns_abs absolute time on the trajectory clock (production starts at ~51 ns after heating+equilibration)
Writes analysis/mmpbsa_perframe_<condition>_<mode>.csv
"""
import csv, sys
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

C = Path("/home3/jocelyn/projects/hippo/campaign_100x40ns")
cond = sys.argv[1] if len(sys.argv) > 1 else "opc-salt"
mode = sys.argv[2] if len(sys.argv) > 2 else "st"
# modes ending in _late are the 40-100 ns window: MM-PBSA startframe=4001, not 1.
STARTFRAME = 4001 if mode.endswith("_late") else 1
META = {"opc-salt": ("OPC", "150 mM NaCl"), "tip3p": ("TIP3P", "150 mM NaCl"),
        "counterion": ("OPC", "counterion only")}
water, salt = META[cond]
SEC = {"Complex Energy Terms": "cpx", "Receptor Energy Terms": "rec",
       "Ligand Energy Terms": "lig", "DELTA Energy Terms": ""}

def parse(p):
    """-> {prefix: {frame: {term: value}}}, and the term order"""
    out, terms, cur, hdr = {}, None, None, None
    for row in csv.reader(open(p)):
        if not row or not row[0].strip():
            continue
        key = row[0].strip()
        if key in SEC:
            cur = SEC[key]; out[cur] = {}; hdr = None; continue
        if key.startswith("Frame"):
            hdr = [c.strip() for c in row[1:]]
            if cur == "": terms = hdr          # take the term list from the DELTA section
            elif terms is None: terms = hdr
            continue
        if cur is not None and hdr and key.lstrip("-").isdigit():
            # the DELTA section names its last columns "DELTA TOTAL" / "DELTA G gas" / "DELTA G solv";
            # normalise so every section uses the same term names
            out[cur][int(key)] = {t[6:] if t.startswith("DELTA ") else t: float(v)
                                  for t, v in zip(hdr, row[1:])}
    return out, [t[6:] if t.startswith("DELTA ") else t for t in (terms or [])]

def frame_times(traj, n_sel, interval=10, startframe=1):
    """startframe is the 1-based MM-PBSA startframe. Default 1 reproduces the original
    behaviour exactly; the 40-100 ns late window uses startframe=4001."""
    f = netcdf_file(str(traj), "r", mmap=True)
    if "time" not in f.variables:
        f.close(); return None, None
    t = np.array(f.variables["time"][:], float); f.close()
    idx = (startframe - 1) + np.arange(0, n_sel * interval, interval)
    idx = idx[idx < len(t)]
    return t[idx] / 1000.0, idx + 1          # ns (absolute), 1-based trajectory frame

rows, terms_ref, missing_time = [], None, 0
base = C/"mmpbsa"/cond/mode
for d in sorted(base.glob("*_[0-9][0-9][0-9]")):
    pf = d/"per_frame.csv"
    if not pf.exists() or pf.stat().st_size == 0:
        continue
    sysn, slot = d.name.rsplit("_", 1)
    data, terms = parse(pf)
    # Single-trajectory files carry a DELTA section. Multi-trajectory files do NOT: the three species come
    # from independent simulations, so MM-PBSA writes no frame-wise difference. For MT we form the per-frame
    # delta ourselves as complex - receptor - ligand at the same frame INDEX. That pairing is bookkeeping,
    # not physical pairing, and is recorded in the delta_method column.
    if "" in data and data[""]:
        method = "mmpbsa_delta_section"
    elif all(k in data and data[k] for k in ("cpx", "rec", "lig")):
        method = "index_matched_cpx_minus_rec_minus_lig"
        terms = terms or [t for t in data["cpx"][sorted(data["cpx"])[0]]]
        data[""] = {k: {t: data["cpx"][k][t] - data["rec"][k].get(t, 0.0) - data["lig"][k].get(t, 0.0)
                        for t in data["cpx"][k]}
                    for k in sorted(set(data["cpx"]) & set(data["rec"]) & set(data["lig"]))}
    else:
        continue
    terms_ref = terms_ref or terms
    n = len(data[""])
    traj = C/"trajectories"/cond/f"{sysn}_complex"/f"{sysn}_complex_{slot}.nc"
    tabs, tfr = frame_times(traj, n, startframe=STARTFRAME) if traj.exists() else (None, None)
    if tabs is None: missing_time += 1
    for k in sorted(data[""]):
        r = dict(system=sysn, water_model=water, salt_condition=salt, protocol=mode.upper(),
                 delta_method=method, replica=int(slot), frame=k,
                 traj_frame=int(tfr[k]) if tfr is not None and k < len(tfr) else "",
                 time_ns=round(float(tabs[k] - tabs[0]), 4) if tabs is not None and k < len(tabs) else "",
                 time_ns_abs=round(float(tabs[k]), 4) if tabs is not None and k < len(tabs) else "")
        for t in terms_ref:
            r[t] = data[""][k].get(t, "")
        for pre in ("cpx", "rec", "lig"):
            for t in terms_ref:
                r[f"{pre}_{t}"] = data.get(pre, {}).get(k, {}).get(t, "")
        rows.append(r)

if not rows:
    sys.exit(f"no completed per_frame.csv found under {base}")
out = C/"analysis"/f"mmpbsa_perframe_{cond}_{mode}.csv"
out.parent.mkdir(exist_ok=True)
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
reps = sorted({(r["system"], r["replica"]) for r in rows})
per = {}
for r in rows: per.setdefault((r["system"], r["replica"]), 0); per[(r["system"], r["replica"])] += 1
bad = [k for k, v in per.items() if v != 400]
print(f"wrote {out}")
print(f"  {len(rows)} rows, {len(reps)} replicas, {len(terms_ref)} energy terms x 4 species")
print(f"  frames per replica: {'all 400' if not bad else f'IRREGULAR: {bad[:5]}'}")
print(f"  replicas missing trajectory time values: {missing_time}")
for s in sorted({r['system'] for r in rows}):
    n = len({r['replica'] for r in rows if r['system'] == s})
    print(f"  {s}: {n} replicas")
