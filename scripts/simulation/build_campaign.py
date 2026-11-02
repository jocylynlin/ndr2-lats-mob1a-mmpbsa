#!/usr/bin/env python3
"""Build the organised directory for the 100 x 40 ns campaign."""
import csv, hashlib, os, re, shutil, sys
from pathlib import Path
from scipy.io import netcdf_file
import numpy as np

H = Path("/home3/jocelyn/projects/hippo/hippo2")
C = Path("/home3/jocelyn/projects/hippo/campaign_100x40ns")
N_TARGET, PROD_NS, DT_PS, FRAME_PS = 100, 40, 2.0, 10.0

# condition -> system -> species -> (setup dir, stripped dir, file stem, reusable indices, campaign prmtop name)
SRC = {
 "opc-salt": {
  "lats1": {"complex": (H/"water_test/setup/moblats-opc-2", H/"water_test/mmpbsa/opc/moblats/trajectory_moblats", "moblats", range(21,41)),
            "mob1a":   (H/"water_test/setup/mobl-opc-2",    H/"water_test/mmpbsa/opc/moblats/trajectory_mob",      "mob",     range(21,41)),
            "kinase":  (H/"water_test/setup/lats-opc-2",    H/"water_test/mmpbsa/opc/moblats/trajectory_lats",     "lats",    range(21,41))},
  "ndr2":  {"complex": (H/"water_test/setup/mobstk-opc-chainB", H/"water_test/mmpbsa/opc/mobndr-chainB/traj_mobstk-chainB", "mobstk", range(1,41)),
            "mob1a":   (H/"water_test/setup/mobs-opc-chainB",   H/"water_test/mmpbsa/opc/mobndr-chainB/traj_mobs-chainB",   "mobS",   range(1,41)),
            "kinase":  (H/"water_test/setup/stk-opc-chainD",    H/"water_test/mmpbsa/opc/mobndr-chainB/traj_ndr-chainD",    "ndrD",   range(1,41))}},
 "tip3p": {
  "lats1": {"complex": (H/"water_test/setup/moblats-tip3p-2", H/"water_test/mmpbsa/tip3p/moblats/trajectory_moblats", "moblats", range(21,41)),
            "mob1a":   (H/"water_test/setup/mobL-tip3p-2",    H/"water_test/mmpbsa/tip3p/moblats/trajectory_mob",     "mob",     range(21,41)),
            "kinase":  (H/"water_test/setup/lats-tip3p-2",    H/"water_test/mmpbsa/tip3p/moblats/trajectory_lats",    "lats",    range(21,41))},
  "ndr2":  {"complex": (H/"water_test/setup/mobstk-tip3p-chainB", H/"water_test/mmpbsa/tip3p/mobndr-chainB/traj_mobstk-chainB", "mobstk", range(1,21)),
            "mob1a":   (H/"water_test/setup/mobs-tip3p-chainB",   H/"water_test/mmpbsa/tip3p/mobndr-chainB/traj_mobs-chainB",   "mobS",   range(1,21)),
            "kinase":  (H/"water_test/setup/stk-tip3p-chainD",    H/"water_test/mmpbsa/tip3p/mobndr-chainB/traj_ndr-chainD",    "ndrD",   range(1,21))}},
 "counterion": {
  "lats1": {"complex": (H/"ion_test/setup/moblats-no-ions-2", H/"ion_test/mmpbsa/moblats-no-ions/traj_moblats", "moblats", range(21,41)),
            "mob1a":   (H/"ion_test/setup/mobL-no-ions-2",    H/"ion_test/mmpbsa/moblats-no-ions/traj_mob",     "mob",     range(21,41)),
            "kinase":  (H/"ion_test/setup/lats-no-ions-2",    H/"ion_test/mmpbsa/moblats-no-ions/traj_lats",    "lats",    range(21,41))},
  "ndr2":  {"complex": (H/"water_test/setup/mobstk-noions-chainB", H/"ion_test/mmpbsa/mobndr-chainB/traj_mobstk-chainB", "mobstk", range(1,21)),
            "mob1a":   (H/"water_test/setup/mobs-noions-chainB",   H/"ion_test/mmpbsa/mobndr-chainB/traj_mobs-chainB",   "mobS",   range(1,21)),
            "kinase":  (H/"water_test/setup/stk-noions-chainD",    None, None, [])}}}

md5 = lambda p: hashlib.md5(Path(p).read_bytes()).hexdigest()
def prm_info(p):
    t = Path(p).read_text(errors="replace")
    F = {m.group(1): m.group(3) for m in re.finditer(r"%FLAG (\S+)\s*\n%FORMAT\(([^)]+)\)\s*\n(.*?)(?=%FLAG|\Z)", t, re.S)}
    ptr = F["POINTERS"].split()
    chg = np.array(F["CHARGE"].split(), dtype=float).sum() / 18.2223
    s = F["RESIDUE_LABEL"].replace("\n", "")
    lab = [s[i:i+4].strip() for i in range(0, len(s), 4)]
    return int(ptr[0]), len(lab), round(float(chg), 4), sum(1 for x in lab if x == "WAT"), \
           sum(1 for x in lab if x == "Na+"), sum(1 for x in lab if x == "Cl-")

rows, problems = [], []
for cond, systems in SRC.items():
    for sysn, species in systems.items():
        for sp, (setup, strip_dir, stem, idxs) in species.items():
            tag = f"{sysn}_{sp}"
            inp = C/"inputs"/cond/tag; inp.mkdir(parents=True, exist_ok=True)
            (C/"runs"/cond/tag).mkdir(parents=True, exist_ok=True)
            trj = C/"trajectories"/cond/tag; trj.mkdir(parents=True, exist_ok=True)
            if not setup.exists(): problems.append(f"missing setup {setup}"); continue
            prm = sorted(setup.glob("*.prmtop"))
            if not prm or not (setup/"min2.rst").exists():
                problems.append(f"{cond}/{tag}: no prmtop or min2.rst in {setup}"); continue
            prm = prm[0]
            # copy the run inputs
            copied = []
            for f in [prm.name, "min2.rst", "min2.out", "min2.in", "min1.in", "tleap.in", "ZAFF.frcmod", "ZAFF.prep",
                      "heating_and_eq.sh", "run_prod.sh", "leap.log"]:
                if (setup/f).exists(): shutil.copy2(setup/f, inp/f); copied.append(f)
            na, nres, chg, wat, na_, cl = prm_info(inp/prm.name)
            # reusable stripped trajectories -> campaign numbering 001..n
            reuse, missing = [], 0
            if strip_dir is not None:
                for slot, i in enumerate(idxs, start=1):
                    src = strip_dir/f"{stem}{i}.nc"
                    if not src.exists(): missing += 1; continue
                    link = trj/f"{tag}_{slot:03d}.nc"
                    if link.is_symlink() or link.exists(): link.unlink()
                    link.symlink_to(src)
                    reuse.append((slot, i, str(src)))
            if missing: problems.append(f"{cond}/{tag}: {missing} expected stripped files missing in {strip_dir}")
            # verify one reused file has enough frames for 40 ns
            frames_ok = ""
            if reuse:
                f = netcdf_file(str(trj/f"{tag}_001.nc"), "r", mmap=True)
                nf = f.variables["coordinates"].shape[0]; f.close()
                need = int(PROD_NS*1000/FRAME_PS)
                frames_ok = f"{nf} frames (need {need})" + ("" if nf >= need else "  TOO SHORT")
                if nf < need: problems.append(f"{cond}/{tag}: reused trajectories have only {nf} frames")
            with (inp/"PROVENANCE.txt").open("w") as fh:
                fh.write(f"campaign slot: {cond} / {tag}\nsource setup: {setup}\n")
                fh.write(f"copied: {', '.join(copied)}\n")
                fh.write(f"{prm.name} md5 {md5(inp/prm.name)}\nmin2.rst md5 {md5(inp/'min2.rst')}\n")
                fh.write(f"topology: {na} atoms, {nres} residues, charge {chg:+.4f} e, {wat} WAT, {na_} Na+, {cl} Cl-\n")
                fh.write(f"reused stripped trajectories: {len(reuse)} (source indices {min(i for _,i,_ in reuse) if reuse else '-'}"
                         f"-{max(i for _,i,_ in reuse) if reuse else '-'})\n")
                for slot, i, src in reuse: fh.write(f"  slot {slot:03d} <- {src}\n")
            if reuse:
                with (trj/"REUSE.csv").open("w", newline="") as fh:
                    w = csv.writer(fh); w.writerow(["campaign_slot", "source_index", "source_path"])
                    w.writerows(reuse)
            rows.append(dict(condition=cond, system=sysn, species=sp, tag=tag, setup=str(setup),
                             prmtop=prm.name, atoms=na, residues=nres, charge=chg, waters=wat, na=na_, cl=cl,
                             reused=len(reuse), new_needed=N_TARGET-len(reuse), reuse_frames=frames_ok,
                             prmtop_md5=md5(inp/prm.name), min2_md5=md5(inp/"min2.rst")))
C.mkdir(exist_ok=True)
with (C/"MANIFEST.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"{'condition':12s}{'system':7s}{'species':9s}{'atoms':>8s}{'charge':>9s}{'reused':>8s}{'new':>6s}  {'reuse frames'}")
for r in rows:
    print(f"{r['condition']:12s}{r['system']:7s}{r['species']:9s}{r['atoms']:8d}{r['charge']:+9.4f}"
          f"{r['reused']:8d}{r['new_needed']:6d}  {r['reuse_frames']}")
print(f"\ntotal new replicas to run: {sum(r['new_needed'] for r in rows)}")
print("problems:" if problems else "no problems found")
for p in problems: print("  " + p)
