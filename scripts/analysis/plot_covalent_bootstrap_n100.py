#!/usr/bin/env python3
r"""SI bootstrap forest for the covalent difference, rebuilt at n = 100.

Same layout as the submitted figure: one row block per condition, TOTAL plus the five
components, observed difference with a 95 per cent percentile interval, and the fraction of
resamples preserving the observed sign printed beside each row. The replica is the unit of
resampling, 10,000 resamples.

Two changes from the submitted version. Values use the POOLED unbound MOB1A reference within
each condition, so the TOTAL rows reproduce Tables 1-3. And a row is marked as sign-stable only
when its interval excludes zero, which is the same threshold the percentage reports, rather than
by eye.
"""
import math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis"
OUT = "/home3/jocelyn/projects/hippo/reviewer_comments/manuscript_rev/Supplemental Information"
TERMS = ["BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"]
ROWS = ["TOTAL"] + TERMS
COND = [("opc-salt", "Reference"), ("tip3p", "TIP3P"), ("counterion", "Counterion only")]
STABLE, UNSTABLE = "#1c6b2f", "#8a8a85"      # status: interval excludes zero / does not
INK, INK2, MUTED, SURFACE, ZERO = "#0b0b0b", "#52514e", "#8a8a85", "#fcfcfb", "#b4443f"
NB = 10000

def per_replica(d, system, prefix, terms, pool=None):
    g = d[d.system == system]
    return sum(g.groupby("replica")[f"{prefix}{t}"].mean().values for t in terms)

fig, axes = plt.subplots(len(COND), 1, figsize=(9.2, 10.4), sharex=True)
fig.patch.set_alpha(0.0)
rng = np.random.default_rng(0)
print("Bootstrap sign stability of the LATS1-minus-NDR2 covalent difference, n=100\n")

for ax, (cond, clab) in zip(axes, COND):
    d = pd.read_csv(f"{A}/mmpbsa_perframe_{cond}_mt.csv")
    ax.patch.set_alpha(0.0)
    print(f"  --- {clab} ---")
    for r, term in enumerate(ROWS):
        tl = TERMS if term == "TOTAL" else [term]
        pool = np.concatenate([per_replica(d, s, "rec_", tl) for s in ("lats1", "ndr2")]).mean()
        v = {s: per_replica(d, s, "cpx_", tl) - pool - per_replica(d, s, "lig_", tl)
             for s in ("lats1", "ndr2")}
        obs = v["lats1"].mean() - v["ndr2"].mean()
        boot = np.array([rng.choice(v["lats1"], 100, replace=True).mean()
                         - rng.choice(v["ndr2"], 100, replace=True).mean() for _ in range(NB)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        frac = (boot > 0).mean() if obs > 0 else (boot < 0).mean()
        col = STABLE if lo * hi > 0 else UNSTABLE          # interval excludes zero
        y = len(ROWS) - 1 - r
        ax.plot([lo, hi], [y, y], color=col, lw=2.0, solid_capstyle="round", zorder=2)
        ax.plot([lo, lo], [y - .12, y + .12], color=col, lw=1.6, zorder=2)
        ax.plot([hi, hi], [y - .12, y + .12], color=col, lw=1.6, zorder=2)
        ax.plot([obs], [y], "o", ms=9 if term == "TOTAL" else 6.5, color=col, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=1.1)
        ax.text(hi + 0.45, y, f"{100*frac:.1f}%", va="center", ha="left", fontsize=9,
                color=col, fontweight="semibold" if col == STABLE else "normal")
        print(f"    {term:9s} {obs:+6.2f}  [{lo:+6.2f}, {hi:+6.2f}]  sign preserved {100*frac:5.1f}%")
    ax.axvline(0, color=ZERO, lw=1.1, ls=(0, (2, 2)), zorder=1)
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS[::-1], fontsize=10, color=INK)
    ax.set_ylim(-0.6, len(ROWS) - 0.4)
    ax.set_ylabel(clab, fontsize=11, color=INK, fontweight="semibold", labelpad=10)
    ax.grid(True, axis="x", color="#e8e8e4", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=9.5, length=3, width=0.8)

axes[-1].set_xlabel(r"LATS1 $-$ NDR2 covalent difference (kcal/mol)", fontsize=11, color=INK2)
# one shared x-range across the blocks, with just enough headroom for the percentage labels
lo = min(ax.get_xlim()[0] for ax in axes)
hi = max(max(l.get_xdata()) for ax in axes for l in ax.lines if len(l.get_xdata()))
for ax in axes: ax.set_xlim(lo, hi + 0.18 * (hi - lo))
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/covalent_bootstrap_forest_n100.{ext}", dpi=300, transparent=True,
                bbox_inches="tight")
print(f"\nwrote {OUT}/covalent_bootstrap_forest_n100.png and .pdf")
