#!/usr/bin/env python3
"""reproduce_tables.py  -- run from the deposit root:  python3 scripts/analysis/reproduce_tables.py

Recomputes the headline numbers of the paper from the files in energies/ alone and compares
each with the value printed in the paper. No absolute paths, no trajectories.

Conventions (Methods, MM-PBSA Calculations):
  * the replica is the statistical unit; SE = s(replica means) / sqrt(N);
  * multi-trajectory (MT) values use the unbound MOB1A reference pooled over the two MOB1A
    preparations within a condition (200 replicas). Per-system SEs include the pooled
    reference; LATS1-minus-NDR2 differences exclude it, because it cancels exactly;
  * the complex, kinase-segment and MOB1A simulations are independent, so their
    replica-level SEs are combined in quadrature; replicas are never paired by index.
"""
import glob, gzip, io, math, os, re, sys, tarfile
import numpy as np, pandas as pd

E = "energies"
se = lambda v: float(np.std(np.asarray(v, float), ddof=1) / math.sqrt(len(v)))
GROUPS = {"Covalent": ["BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"], "Electrostatic": ["EEL", "EPB"],
          "Dispersion": ["VDWAALS", "EDISPER"], "Hydrophobic": ["ENPOLAR"], "dH": ["TOTAL"]}
SYS = ("lats1", "ndr2")
fails = []

def check(label, got, want, tol=0.006):
    ok = all(abs(g - w) <= tol for g, w in zip(got, want))
    if not ok: fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:46s} " + "  ".join(f"{g:+.2f}" for g in got)
          + ("" if ok else "   paper: " + "  ".join(f"{w:+.2f}" for w in want)))

def table(cond, proto):
    return pd.read_csv(f"{E}/mmpbsa_perframe/mmpbsa_perframe_{cond}_{proto}.csv.gz")

def rep_means(d, cols):
    """-> {system: array of replica means of sum(cols)}"""
    d = d.assign(_v=d[cols].sum(axis=1))
    return {s: d[d.system == s].groupby("replica")._v.mean().values for s in SYS}

def st(cond, terms):
    v = rep_means(table(cond, "st"), terms)
    return {s: (v[s].mean(), se(v[s])) for s in SYS} | {
        "diff": (v["lats1"].mean() - v["ndr2"].mean(), math.hypot(se(v["lats1"]), se(v["ndr2"])))}

_mt = {}
def mt_parts(cond, terms):
    d = _mt.setdefault(cond, table(cond, "mt"))
    c = rep_means(d, ["cpx_" + t for t in terms]); r = rep_means(d, ["rec_" + t for t in terms])
    l = rep_means(d, ["lig_" + t for t in terms])
    return c, r, l

def mt(cond, terms):
    c, r, l = mt_parts(cond, terms)
    pool = np.concatenate([r[s] for s in SYS])
    q = lambda *v: math.sqrt(sum(se(x) ** 2 for x in v))
    out = {s: (c[s].mean() - l[s].mean() - pool.mean(), q(c[s], l[s], pool)) for s in SYS}
    out["diff"] = ((c["lats1"].mean() - l["lats1"].mean()) - (c["ndr2"].mean() - l["ndr2"].mean()),
                   q(c["lats1"], l["lats1"], c["ndr2"], l["ndr2"]))
    return out

flat = lambda r: [r["lats1"][0], r["lats1"][1], r["ndr2"][0], r["ndr2"][1], r["diff"][0], r["diff"][1]]

print("Table 1 (reference condition: OPC, added NaCl)")
check("ST dH   L, se, N, se, diff, se", flat(st("opc-salt", ["TOTAL"])), [-54.74, .51, -56.65, .47, 1.91, .69])
ref = {g: mt("opc-salt", t) for g, t in GROUPS.items()}
want = {"Covalent": [-3.14, 0.70, -5.63, 0.90, 2.49, 1.06], "Electrostatic": [-15.62, .21, -14.32, .20, -1.30, .27],
        "Dispersion": [54.50, 0.65, 52.48, 0.67, 2.02, 0.85], "Hydrophobic": [-88.67, 1.06, -87.38, 1.08, -1.29, 1.36],
        "dH": [-52.93, 1.19, -54.85, 1.33, 1.92, 1.61]}
for g in GROUPS: check(f"MT {g}", flat(ref[g]), want[g])

print("Table 2 (TIP3P)")
tip = {g: mt("tip3p", t) for g, t in GROUPS.items()}
want = {"Covalent": [-8.49, .71, -6.13, .82, -2.36, .99], "Electrostatic": [-15.71, 0.18, -15.05, 0.19, -0.67, 0.24],
        "Dispersion": [56.17, 0.62, 51.12, 0.79, 5.06, 0.93], "Hydrophobic": [-89.35, 0.87, -82.73, 0.96, -6.62, 1.18],
        "dH": [-57.38, 1.13, -52.79, 1.17, -4.59, 1.45]}
for g in GROUPS: check(f"MT {g}", flat(tip[g]), want[g])
check("ST dH diff, se (TIP3P)", list(st("tip3p", ["TOTAL"])["diff"]), [0.32, .68])
d1, d2 = ref["dH"]["diff"], tip["dH"]["diff"]
check("MT water-model contrast, se", [d1[0] - d2[0], math.hypot(d1[1], d2[1])], [6.51, 2.17])

print("Table 3 (counterion only)")
cnt = {g: mt("counterion", t) for g, t in GROUPS.items()}
want = {"Covalent": [-5.99, 0.79, -5.47, 0.89, -0.52, 1.10], "Electrostatic": [-14.94, 0.18, -14.60, 0.18, -0.34, 0.23],
        "Dispersion": [57.49, 0.60, 51.38, 0.60, 6.12, 0.77], "Hydrophobic": [-94.36, 0.99, -90.39, 1.08, -3.98, 1.32],
        "dH": [-57.80, 1.29, -59.08, 1.28, 1.28, 1.63]}
for g in GROUPS: check(f"MT {g}", flat(cnt[g]), want[g])
d3 = cnt["dH"]["diff"]
check("ionic-condition contrast, se", [d1[0] - d3[0], math.hypot(d1[1], d3[1])], [0.64, 2.29])

print("Table S1 (covalent components, LATS1 - NDR2, reference / TIP3P)")
for cond, w in (("opc-salt", {"BOND": (-.17, .21), "ANGLE": (-0.27, 0.48), "DIHED": (2.88, 0.84), "1-4 VDW": (-0.96, 0.21), "1-4 EEL": (1.01, 0.14)}),
                ("tip3p", {"BOND": (-0.08, 0.20), "ANGLE": (-0.99, 0.45), "DIHED": (-0.43, 0.76), "1-4 VDW": (-1.59, .20), "1-4 EEL": (.73, .13)})):
    for t, v in w.items(): check(f"{cond} {t}", list(mt(cond, [t])["diff"]), list(v))

print("Unbound MOB1A preparations and own-reference estimate (Section 2.1)")
c, r, l = mt_parts("opc-salt", ["TOTAL"])
ctrl = r["lats1"].mean() - r["ndr2"].mean()
check("preparation difference, se", [ctrl, math.hypot(se(r["lats1"]), se(r["ndr2"]))], [5.74, 1.01])
own = d1[0] - ctrl
q = math.sqrt(sum(se(x[s]) ** 2 for x in (c, r, l) for s in SYS))
check("own-reference ddH, se (quadrature)", [own, q], [-3.82, 1.90])

print("Table S5 (first vs second half of the 40 ns window, reference condition)")
def halves(d, cols):
    d = d.assign(_v=d[cols].sum(axis=1), _h=(d.frame >= 200))
    g = d.groupby(["system", "replica", "_h"])._v.mean().unstack()
    return {s: (g.loc[s][True] - g.loc[s][False]).values for s in SYS}
w = halves(table("opc-salt", "st"), ["TOTAL"])
check("ST late-early L, se, N, se", [w["lats1"].mean(), se(w["lats1"]), w["ndr2"].mean(), se(w["ndr2"])], [-.19, .38, -.20, .33])
check("ST contrast, se", [w["lats1"].mean() - w["ndr2"].mean(), math.hypot(se(w["lats1"]), se(w["ndr2"]))], [.01, .50])
dm = _mt["opc-salt"]
wc, wr, wl = halves(dm, ["cpx_TOTAL"]), halves(dm, ["rec_TOTAL"]), halves(dm, ["lig_TOTAL"])
pool = np.concatenate([wr[s] for s in SYS])
q = lambda *v: math.sqrt(sum(se(x_) ** 2 for x_ in v))
check("pooled MOB1A drift, se", [pool.mean(), se(pool)], [.96, .59])
check("MT late-early L, se, N, se", [wc["lats1"].mean() - wl["lats1"].mean() - pool.mean(), q(wc["lats1"], wl["lats1"], pool),
                                     wc["ndr2"].mean() - wl["ndr2"].mean() - pool.mean(), q(wc["ndr2"], wl["ndr2"], pool)], [-3.84, 1.30, -2.01, 1.24])
check("MT contrast, se", [(wc["lats1"].mean() - wl["lats1"].mean()) - (wc["ndr2"].mean() - wl["ndr2"].mean()),
                          q(wc["lats1"], wl["lats1"], wc["ndr2"], wl["ndr2"])], [-1.83, 1.60])

print("Table S3 (normal-mode entropy, -T dS at 310.15 K)")
def nm(fn):
    txt = gzip.open(fn, "rt").read().split("POISSON BOLTZMANN:")[0]
    out = {}
    for key, lab in (("Complex", "cpx"), ("Receptor", "rec"), ("Ligand", "lig")):
        m = re.search(key + r" Energy Terms\nFrame #,Translational,Rotational,Vibrational,Total\n((?:[-\d].*\n)+)", txt)
        if m is None: return None
        out[lab] = np.array([[float(v) for v in ln.split(",")] for ln in m.group(1).strip().split("\n")])
    n = min(len(v) for v in out.values())
    return {k: v[:n, 1:5].mean(axis=0) for k, v in out.items()}
N = {s: [o for o in (nm(f) for f in sorted(glob.glob(f"{E}/nmode/opc-salt/{s}/traj_*/per_frame.csv.gz"))) if o] for s in SYS}
print("  replicas:", {s: len(N[s]) for s in SYS})
k = 310.15 / 298.15
cx = {s: np.array([o["cpx"] for o in N[s]]) for s in SYS}; lg = {s: np.array([o["lig"] for o in N[s]]) for s in SYS}
pool = np.array([o["rec"] for s in SYS for o in N[s]])
ent = {}
for j, name in enumerate(("trans", "rot", "vib", "total")):
    row = []
    for s in SYS:
        row += [k * (-(cx[s][:, j].mean() - lg[s][:, j].mean()) + pool[:, j].mean()), k * q(cx[s][:, j], lg[s][:, j], pool[:, j])]
    row += [k * (-(cx["lats1"][:, j].mean() - lg["lats1"][:, j].mean()) + (cx["ndr2"][:, j].mean() - lg["ndr2"][:, j].mean())),
            k * q(cx["lats1"][:, j], lg["lats1"][:, j], cx["ndr2"][:, j], lg["ndr2"][:, j])]
    ent[name] = row
want = {"trans": [16.00, 0, 15.99, 0, 0, 0], "rot": [15.66, .01, 15.64, 0, .02, .01], "vib": [33.28, 0.82, 30.72, 0.76, 2.56, 1.02], "total": [64.94, 0.82, 62.35, 0.76, 2.59, 1.02]}
for name in ent: check(f"-TdS {name}", ent[name], want[name])
ds = ent["total"][4:]
check("dG L, se, N, se", [ref["dH"]["lats1"][0] + ent["total"][0], math.hypot(ref["dH"]["lats1"][1], ent["total"][1]),
                          ref["dH"]["ndr2"][0] + ent["total"][2], math.hypot(ref["dH"]["ndr2"][1], ent["total"][3])], [12.01, 1.44, 7.50, 1.53])
check("ddG, se", [d1[0] + ds[0], math.hypot(d1[1], ds[1])], [4.51, 1.91])

print("SI Section S5 (terminal / core partition of the kinase bonded energy)")
COV = ["bond", "angle", "dih", "vdw14", "elec14"]
tc = {}
with tarfile.open(f"{E}/terminal_core/terminal_core_opc-salt.tar.gz") as tf:
    for m in tf.getmembers():
        g = re.match(r"(?:\./)?(lats1|ndr2)_(bound|unbound)_(\d{3})_(full|core)\.dat$", m.name)
        if not g: continue
        txt = tf.extractfile(m).read().decode()
        hdr = txt.split("\n", 1)[0].lstrip("#").split()
        a = np.loadtxt(io.StringIO(txt), comments="#", ndmin=2)
        idx = [hdr.index(f"{g.group(4).upper()}[{cc}]") for cc in COV]
        tc[g.groups()] = a[:, idx].sum(axis=1).mean()
def part(s, state, p):
    reps = sorted({kk[2] for kk in tc if kk[0] == s and kk[1] == state})
    f = np.array([tc[(s, state, r_, "full")] for r_ in reps]); co = np.array([tc[(s, state, r_, "core")] for r_ in reps])
    return {"full": f, "core": co, "termini": f - co}[p]
res = {}
for s in SYS:
    for p in ("full", "core", "termini"):
        b, u = part(s, "bound", p), part(s, "unbound", p)
        res[(s, p)] = (b.mean() - u.mean(), math.hypot(se(b), se(u)))
check("LATS1 total, core, termini", [res[("lats1", p)][0] for p in ("full", "core", "termini")], [-6.65, -3.37, -3.27])
check("NDR2 total, core, termini", [res[("ndr2", p)][0] for p in ("full", "core", "termini")], [-14.67, -7.04, -7.64])
for p, wv in (("termini", (4.36, .73)), ("core", (3.66, 1.73))):
    check(f"LATS1 - NDR2 {p}, se", [res[("lats1", p)][0] - res[("ndr2", p)][0], math.hypot(res[("lats1", p)][1], res[("ndr2", p)][1])], list(wv))

print("Interface (Sections 2.2 and 2.3): buried area and salt bridges, per-replica means")
def tar_means(path, col):
    out = {s: [] for s in SYS}
    with tarfile.open(path) as tf:
        for m in sorted(tf.getmembers(), key=lambda m_: m_.name):
            g = re.match(r"(?:\./)?(lats1|ndr2)_\d{3}\.csv$", m.name)
            if g: out[g.group(1)].append(pd.read_csv(tf.extractfile(m))[col].mean())
    return {s: np.array(v) for s, v in out.items()}
B = {c_: tar_means(f"{E}/interface/buried_area_{c_}.tar.gz", "buried") for c_ in ("opc-salt", "tip3p", "counterion")}
SB = {c_: tar_means(f"{E}/interface/saltbridge_{c_}.tar.gz", "total") for c_ in ("opc-salt", "tip3p", "counterion")}
dd = lambda a, b: [a.mean() - b.mean(), math.hypot(se(a), se(b))]
check("buried area TIP3P - OPC, L and N", dd(B["tip3p"]["lats1"], B["opc-salt"]["lats1"]) + dd(B["tip3p"]["ndr2"], B["opc-salt"]["ndr2"]), [99.4, 19.9, 48.7, 20.4], tol=0.06)
check("salt bridges TIP3P - OPC, L and N", dd(SB["tip3p"]["lats1"], SB["opc-salt"]["lats1"]) + dd(SB["tip3p"]["ndr2"], SB["opc-salt"]["ndr2"]), [.73, .09, .50, .09])
check("salt bridges counterion - added, L and N", dd(SB["counterion"]["lats1"], SB["opc-salt"]["lats1"]) + dd(SB["counterion"]["ndr2"], SB["opc-salt"]["ndr2"]), [.40, .10, .24, .10])
check("buried area counterion - added, L and N", dd(B["counterion"]["lats1"], B["opc-salt"]["lats1"]) + dd(B["counterion"]["ndr2"], B["opc-salt"]["ndr2"]), [-5.8, 17.3, 32.7, 18.3], tol=0.06)

print(f"\n{'ALL CHECKS PASSED' if not fails else str(len(fails)) + ' CHECK(S) FAILED: ' + '; '.join(fails)}")
sys.exit(1 if fails else 0)
