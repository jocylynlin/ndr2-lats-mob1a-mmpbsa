#!/usr/bin/env python3
r"""Main-text covalent decomposition figure, regenerated at n = 100.

Same layout as the submitted figure: grouped bars, one pair per internal-coordinate component,
replica-level error bars. Values are recomputed from the per-frame table rather than transcribed,
using the POOLED unbound MOB1A reference so the bars sum to the covalent totals in Table 1
(-2.33 for LATS1, -5.66 for NDR2) and their difference to +3.34.

Palette: categorical slots 1 and 2, validated (CVD dE 24.7, normal-vision dE 33.6, all checks
pass) -- the same pair used by the hydration figures.
"""
import math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis"
OUT = "/home3/jocelyn/projects/hippo/reviewer_comments/manuscript_rev"
LATS1, NDR2 = "#2a78d6", "#eb6834"
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#8a8a85", "#fcfcfb"
COV = ["BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"]

mt = pd.read_csv(f"{A}/mmpbsa_perframe_opc-salt_mt.csv")
se = lambda v: np.std(v, ddof=1) / math.sqrt(len(v))

def per_replica(system, prefix, term):
    d = mt[mt.system == system]
    return d.groupby("replica")[prefix + term].mean().values

vals, errs = {}, {}
for t in COV:
    # pooled MOB1A reference: both systems' unbound ensembles combined (200 replicas)
    pool = np.concatenate([per_replica(s, "rec_", t) for s in ("lats1", "ndr2")])
    for s in ("lats1", "ndr2"):
        cpx, lig = per_replica(s, "cpx_", t), per_replica(s, "lig_", t)
        vals[(s, t)] = (cpx - pool.mean() - lig).mean()
        errs[(s, t)] = math.hypot(se(cpx), se(pool), se(lig))

print("per-system delta components (pooled reference), kcal/mol")
for s in ("lats1", "ndr2"):
    tot = sum(vals[(s, t)] for t in COV)
    print(f"  {s}: " + "  ".join(f"{t} {vals[(s,t)]:+.2f}" for t in COV) + f"   sum {tot:+.2f}")
print(f"  difference of sums: {sum(vals[('lats1',t)] for t in COV) - sum(vals[('ndr2',t)] for t in COV):+.2f}"
      "   (Table 1 covalent difference: +3.34)")

x = np.arange(len(COV)); w = 0.36
fig, ax = plt.subplots(figsize=(7.4, 5.0))
fig.patch.set_alpha(0.0); ax.patch.set_alpha(0.0)
for off, s, col, lab in [(-w/2, "lats1", LATS1, "MOB1A--LATS1"), (w/2, "ndr2", NDR2, "MOB1A--NDR2")]:
    y = [vals[(s, t)] for t in COV]; e = [errs[(s, t)] for t in COV]
    ax.bar(x + off, y, w, color=col, label=lab.replace("--", "–"), zorder=2,
           edgecolor=SURFACE, linewidth=1.2)
    ax.errorbar(x + off, y, yerr=e, fmt="none", ecolor=INK, elinewidth=1.3,
                capsize=3.5, capthick=1.3, zorder=3)
ax.axhline(0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
ax.set_xticks(x); ax.set_xticklabels(COV, fontsize=10.5, color=INK)
ax.set_ylabel(r"Grand mean $\Delta$Component (kcal/mol)", fontsize=11, color=INK2)
ax.grid(True, axis="y", color="#e8e8e4", lw=0.6); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.8)
ax.tick_params(colors=INK2, labelsize=10, length=3, width=0.8)
leg = ax.legend(loc="lower right", frameon=False, fontsize=10.5)
for t in leg.get_texts(): t.set_color(INK2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/covalent_decomposition_n100.{ext}", dpi=300, transparent=True,
                bbox_inches="tight")
print(f"wrote {OUT}/covalent_decomposition_n100.png and .pdf")
