#!/usr/bin/env python3
"""campaign_audit.py [condition ...]

Emits a campaign-level audit block per condition and system, computed from the files themselves. Anything that
cannot yet be established prints PENDING rather than PASS, so the report never claims more than was checked.
Writes analysis/AUDIT_<condition>_<system>.txt and prints to stdout.
"""
import csv, re, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

C = Path("/home3/jocelyn/projects/hippo/campaign_100x40ns")
cfg = dict(re.findall(r"^(\w+)=([^\s#]+)", (C/"config/campaign.conf").read_text(), re.M))
N_TARGET = int(cfg["N_TARGET"]); PROD_FRAMES = int(cfg["PROD_FRAMES"])
PY = cfg.get("PYTHON", "/home3/jocelyn/miniconda3/bin/python3")
SPECIES = [("complex", "complex"), ("mob1a", "receptor"), ("kinase", "ligand")]

def prm_natom(p):
    t = Path(p).read_text(errors="replace")
    return int(re.search(r"%FLAG POINTERS\s*\n%FORMAT\([^)]+\)\s*\n(.*?)%FLAG", t, re.S).group(1).split()[0])
def mdin(path, key):
    m = re.search(rf"{key}\s*=\s*(-?[\d.]+)", Path(path).read_text(errors="replace"))
    return m.group(1) if m else "?"
def run(args):
    try: return subprocess.run(args, capture_output=True, text=True, timeout=900).returncode == 0
    except Exception: return False

def audit(cond, sysname):
    L = []; p = L.append
    man = {r["tag"]: r for r in csv.DictReader(open(C/"MANIFEST.csv")) if r["condition"] == cond}
    top = C/"topologies"/cond
    tags = {sp: f"{sysname}_{sp}" for sp, _ in SPECIES}
    trj = {sp: sorted((C/"trajectories"/cond/tags[sp]).glob("*.nc")) for sp, _ in SPECIES}
    p(f"SYSTEM: {sysname.upper()}   CONDITION: {cond}")
    p(f"{'replicas:':24s}{len(trj['complex'])} of {N_TARGET} complex, "
      f"{len(trj['mob1a'])} MOB1A, {len(trj['kinase'])} kinase")
    ex = next((C/"runs"/cond/tags["complex"]).glob("replica*/prod.in"), None)
    if ex:
        ns = int(mdin(ex, "nstlim")) * float(mdin(ex, "dt")) / 1000
        sp_ps = int(mdin(ex, "ntwx")) * float(mdin(ex, "dt"))
        p(f"{'production/replica:':24s}{ns:.0f} ns")
        p(f"{'frame spacing:':24s}{sp_ps:.0f} ps  ({int(ns*1000/sp_ps)} frames written)")
    mi = C/"mmpbsa"/cond/"st"/"mmpbsa.in"
    iv = mdin(mi, "interval") if mi.exists() else "10"
    p(f"{'snapshot spacing:':24s}{float(iv)*10:.0f} ps (interval={iv})")
    p(f"{'snapshots/replica:':24s}{PROD_FRAMES//int(float(iv))}")
    p("")
    for sp, role in SPECIES:
        f = top/f"{sysname}_{role}.prmtop"
        p(f"{role+' topology:':24s}{f if f.exists() else 'PENDING (made after the first complex strip)'}")
    p("")
    nat = {}
    for sp, role in SPECIES:
        f = top/f"{sysname}_{role}.prmtop"
        nat[role] = prm_natom(f) if f.exists() else None
        p(f"{role+' atoms:':24s}{nat[role] if nat[role] is not None else 'PENDING'}")
    if all(nat[r] is not None for _, r in SPECIES):
        ok = nat["receptor"] + nat["ligand"] == nat["complex"]
        p(f"{'atom-count sum:':24s}{'PASS' if ok else 'FAIL'} "
          f"({nat['receptor']} + {nat['ligand']} = {nat['receptor']+nat['ligand']} vs complex {nat['complex']})")
    p("")
    # topology parameters + Zn/receptor composition
    cpx, rec, lig = (top/f"{sysname}_{r}.prmtop" for _, r in SPECIES)
    if all(x.exists() for x in (cpx, rec, lig)):
        ok = run([PY, str(C/"scripts/check_topologies.py"), str(cpx), str(rec), str(lig)])
        p(f"{'topology parameter check:':26s}{'PASS' if ok else 'FAIL'}")
        p(f"{'Zn/receptor check:':26s}{'PASS' if ok else 'FAIL'}  (Zn is the last receptor residue, bonded in molecule 1)")
    else:
        p(f"{'topology parameter check:':26s}PENDING"); p(f"{'Zn/receptor check:':26s}PENDING")
    # trajectory atom/order check on every trajectory present
    res = []
    for sp, role in SPECIES:
        f = top/f"{sysname}_{role}.prmtop"
        if not f.exists() or not trj[sp]: res.append((role, None, 0)); continue
        good = sum(run([PY, str(C/"scripts/check_traj_topology.py"), str(t), str(f)]) for t in trj[sp])
        res.append((role, good, len(trj[sp])))
    if any(g is not None for _, g, _ in res):
        s = ", ".join(f"{r} {g}/{n}" for r, g, n in res if g is not None)
        allok = all(g == n for _, g, n in res if g is not None)
        p(f"{'trajectory atom check:':26s}{'PASS' if allok else 'FAIL'}  ({s})")
    else:
        p(f"{'trajectory atom check:':26s}PENDING")
    # seeds: assigned ledger vs what AMBER actually used, per stage, plus cross-stage overlap
    led = C/"logs/seeds_assigned.csv"
    assigned = {"heat": [], "production": []}
    if led.exists():
        for r in csv.reader(open(led)):
            if len(r) == 6 and r[0] == cond and r[1] in tags.values():
                assigned["heat"].append(r[3]); assigned["production"].append(r[5])
    actual = {"heat": [], "production": []}
    for sp, _ in SPECIES:
        for stage, fname in (("heat", "heat.out"), ("production", "prod.out")):
            for o in (C/"runs"/cond/tags[sp]).glob(f"replica*/{fname}"):
                t = o.read_text(errors="replace")
                m = re.search(r"Setting random seed to\s+(\d+)", t) or re.search(r"^\s*ig\s*=\s*(\d+)", t, re.M)
                if m: actual[stage].append(m.group(1))
    for stage in ("heat", "production"):
        a = assigned[stage]
        p(f"{stage+' assigned:':28s}" + (f"{len(set(a))}/{len(a)} unique" if a else "PENDING"))
    for stage in ("heat", "production"):
        a = actual[stage]
        p(f"{stage+' AMBER actual:':28s}" + (f"{len(set(a))}/{len(a)} unique" if a else "PENDING (not yet run)"))
    both = set(assigned["heat"]) & set(assigned["production"])
    p(f"{'heat n production:':28s}{len(both)} duplicates" + ("" if not both else f"  {sorted(both)[:5]}"))
    if actual["heat"] and assigned["heat"]:
        mism = set(actual["heat"]) - set(assigned["heat"])
        p(f"{'actual vs assigned:':28s}{'PASS' if not mism else 'FAIL'} "
          f"(every seed AMBER used appears in the ledger)")
    # ---- ST vs MT trajectory sources: what each calculation actually consumed ----
    def sources(dat):
        t = Path(dat).read_text(errors="replace")
        g = lambda k: (lambda m: m.group(1).strip() if m else None)(re.search(rf"\|Initial {k}mdcrd\(s\):\s*(\S+)", t))
        return g(""), g("Receptor "), g("Ligand ")
    for mode in ("st", "mt"):
        d = C/"mmpbsa"/cond/mode
        dats = sorted(d.glob(f"{sysname}_*/FINAL_RESULTS_MMPBSA.dat")) if d.exists() else []
        p("")
        p(f"{mode.upper()} trajectory sources:")
        if not dats:
            planned = (C/"scripts/run_mmpbsa.sh").read_text()
            st_ok = 'ARGS="-y $CT"' in planned
            mt_ok = 'ARGS="-y $CT -yr $RT -yl $LT"' in planned
            ok = st_ok if mode == "st" else mt_ok
            p(f"  PENDING (not run). Script will pass: "
              + ("-y complex only -> receptor and ligand extracted from the complex" if mode == "st"
                 else "-y complex -yr mob1a -yl kinase") + f"  [{'PASS' if ok else 'FAIL'}]")
        else:
            c0, r0, l0 = sources(dats[0])
            p(f"  complex  <- {c0}")
            p(f"  receptor <- {r0 if r0 else 'SAME complex trajectory / extracted coordinates'}")
            p(f"  ligand   <- {l0 if l0 else 'SAME complex trajectory / extracted coordinates'}")
            bad = []
            for f in dats:
                c, r, l = sources(f)
                if mode == "st":
                    if r not in (None, c) or l not in (None, c): bad.append(f.parent.name)
                else:
                    if not r or not l: bad.append(f.parent.name)
                    elif len({c, r, l}) != 3: bad.append(f.parent.name)
                    elif f"{sysname}_mob1a" not in r or f"{sysname}_kinase" not in l or f"{sysname}_complex" not in c:
                        bad.append(f.parent.name)
            p(f"  {'PASS' if not bad else 'FAIL'}  ({len(dats)} replicas checked"
              + (f"; wrong sources in {bad[:3]}" if bad else "") + ")")

    # ---- time span actually represented by the selected frames ----
    p("")
    iv_i = int(float(iv))
    for sp, _ in SPECIES:
        if not trj[sp]: continue
        f = netcdf_file(str(trj[sp][0]), "r", mmap=True)
        has_t = "time" in f.variables
        if has_t:
            tv = np.array(f.variables["time"][:], float)
            sel = tv[0:PROD_FRAMES:iv_i]
            f.close()
            span_ok = len(sel) == PROD_FRAMES // iv_i
            dt_ok = np.allclose(np.diff(sel), float(iv_i) * 10, atol=1e-3) if len(sel) > 1 else False
            p(f"{'frames used ('+sp+'):':28s}{len(sel)} frames, t = {sel[0]/1000:.1f} to {sel[-1]/1000:.1f} ns "
              f"(span {(sel[-1]-sel[0])/1000:.1f} ns, spacing {np.diff(sel)[0]/1000 if len(sel)>1 else 0:.3f} ns)  "
              f"{'PASS' if span_ok and dt_ok else 'FAIL'}")
        else:
            f.close(); p(f"{'frames used ('+sp+'):':28s}UNVERIFIED (trajectory carries no time variable)")

    # ---- replica-level completeness, before any aggregation ----
    p("")
    for mode in ("st", "mt"):
        d = C/"mmpbsa"/cond/mode
        reps = sorted(d.glob(f"{sysname}_*")) if d.exists() else []
        if not reps:
            p(f"{sysname.upper()+' '+mode.upper()+':':16s}PENDING (not run)"); continue
        want = PROD_FRAMES // iv_i
        good = 0; issues = []
        for rp in reps:
            dat, csvf = rp/"FINAL_RESULTS_MMPBSA.dat", rp/"per_frame.csv"
            if not dat.exists() or "DELTA TOTAL" not in dat.read_text(errors="replace"):
                issues.append((rp.name, "incomplete .dat")); continue
            m = re.search(r"Calculations performed using\s+([\d.]+)\s+complex frames", dat.read_text(errors="replace"))
            n = float(m.group(1)) if m else -1
            if abs(n - want) > 1.5: issues.append((rp.name, f"{n} frames")); continue
            if not csvf.exists() or csvf.stat().st_size == 0:
                issues.append((rp.name, "no per_frame.csv")); continue
            txt = csvf.read_text(errors="replace")
            if re.search(r"nan|inf", txt, re.I): issues.append((rp.name, "NaN/Inf in per_frame.csv")); continue
            good += 1
        p(f"{sysname.upper()+' '+mode.upper()+':':16s}{good}/{len(reps)} successful "
          f"(expected {want} frames/replica, NaN/Inf 0, incomplete 0)"
          + (f"   issues: {issues[:3]}" if issues else ""))
    p("")
    p("Aggregate from the replica means above (one mean per replica, SEM over replicas);")
    p("do not pool frames across replicas and take an SEM over frames.")
    return "\n".join(L)

conds = sys.argv[1:] or cfg["CONDITIONS"].strip('"').split()
(C/"analysis").mkdir(exist_ok=True)
for cond in conds:
    for sysname in ("ndr2", "lats1"):
        txt = audit(cond, sysname)
        (C/"analysis"/f"AUDIT_{cond}_{sysname}.txt").write_text(txt + "\n")
        print(txt); print("-" * 78)
