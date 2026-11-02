#!/usr/bin/env python3
"""Reviewer 2, Comment 2: are 10 normal-mode snapshots per replica enough?
Uses the existing nmode runs (hippo-claude/nmode): separate runs with 2, 9 and 10 snapshots per replica, and the
per-frame values of the 10-snapshot run. -TdS reported at 310.15 K (S from Amber at its default 298.15 K, rescaled by
T, as in Table S3). Replica-level SE = SD of replica means / sqrt(n). NOTE: NDR2 = original (chain A) model.
Output: results/nmode_snapshot_convergence.csv and a text summary."""
from pathlib import Path
import numpy as np, pandas as pd
N = Path("/home3/jocelyn/projects/hippo-claude/nmode"); RES = Path(__file__).resolve().parents[1] / "results"
F = 310.15 / 298.15
se = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
rows = []
print("-TdS at 310.15 K (kcal/mol), mean ± replica-level SE [frame-level SE for comparison]")
for k in (2, 9, 10):
    res = {}
    for sysn, d in (("LATS1", "moblats"), ("NDR2", "mobndr")):
        x = pd.read_csv(N / f"{d}_opc_salt_{k}frames/nmode_entropy_{k}frames_parsed.csv").neg_TdS.values * F
        res[sysn] = x; rows.append(dict(source=f"{k}-snapshot run", system=sysn, mean=x.mean(), se=se(x), n=len(x)))
    diff = res["LATS1"].mean() - res["NDR2"].mean(); dse = np.hypot(se(res["LATS1"]), se(res["NDR2"]))
    rows.append(dict(source=f"{k}-snapshot run", system="LATS1-NDR2", mean=diff, se=dse, n=None))
    print(f"  {k:2d} snapshots/replica: LATS1 {res['LATS1'].mean():6.2f} ± {se(res['LATS1']):.2f}   NDR2 {res['NDR2'].mean():6.2f} ± {se(res['NDR2']):.2f}   difference {diff:+.2f} ± {dse:.2f}")
print("\nFrom the per-frame values of the 10-snapshot run:")
pf = {s: pd.read_csv(N / f"{d}_opc_salt_10frames/entropy_perframe_all.csv") for s, d in (("LATS1", "moblats"), ("NDR2", "mobndr"))}
for s, d in pf.items():
    v = -d.TdS_total
    rep = v.groupby(d.replica)
    within = rep.std(ddof=1).mean(); between = rep.mean().std(ddof=1)
    print(f"  {s}: frames {len(v)}, replicas {rep.ngroups}; mean {v.mean():.2f}; frame-level SE {v.std(ddof=1)/np.sqrt(len(v)):.2f}; "
          f"replica-level SE {rep.mean().std(ddof=1)/np.sqrt(rep.ngroups):.2f}; SD within replica {within:.2f}, SD of replica means {between:.2f}")
print("\nCumulative estimate using the first k snapshots of each replica (10-snapshot run):")
for k in (2, 4, 6, 8, 10):
    m = {}
    for s, d in pf.items():
        v = d[d.frame <= k].groupby("replica").apply(lambda g: -g.TdS_total.mean()); m[s] = v
    diff = m["LATS1"].mean() - m["NDR2"].mean(); dse = np.hypot(se(m["LATS1"]), se(m["NDR2"]))
    rows.append(dict(source=f"first {k} of 10", system="LATS1-NDR2", mean=diff, se=dse, n=None))
    print(f"  k={k:2d}: LATS1 {m['LATS1'].mean():6.2f} ± {se(m['LATS1']):.2f}   NDR2 {m['NDR2'].mean():6.2f} ± {se(m['NDR2']):.2f}   difference {diff:+.2f} ± {dse:.2f}")
pd.DataFrame(rows).to_csv(RES / "nmode_snapshot_convergence.csv", index=False)
