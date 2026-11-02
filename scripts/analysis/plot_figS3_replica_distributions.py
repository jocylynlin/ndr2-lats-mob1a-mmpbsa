#!/usr/bin/env python3
r"""Supplementary Figure S3: replica-level distributions, as promised in Reviewer 1 Response 4.

Panels A and B: single- and multi-trajectory dHbind per replica. Panel C: the multi-trajectory
covalent term. Panel D: bootstrap distributions of the LATS1-minus-NDR2 ddHbind from 4,000
replica-level resamples. Violins carry every replica as a point, with the mean and its
replica-level standard error marked.

Multi-trajectory values use the POOLED unbound MOB1A reference within each condition, matching
Tables 1-3, so the panel means reproduce the table entries.

Palettes validated separately per panel, since panels are never compared directly:
  A-C systems  blue/orange  CVD dE 24.7, normal 33.6, contrast 3.12
  D  protocols green/violet CVD dE 27.9, normal 34.1, contrast 4.82
"""
import math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis"
OUT = "/home3/jocelyn/projects/hippo/reviewer_comments/manuscript_rev/Supplemental Information"
COND = [("opc-salt", "OPC"), ("tip3p", "TIP3P"), ("counterion", "Counterion\nonly")]
SYS = [("lats1", "LATS1"), ("ndr2", "NDR2")]
C_SYS = {"lats1": "#2a78d6", "ndr2": "#eb6834"}
C_PROT = {"ST": "#008300", "MT": "#4a3aa7"}
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#8a8a85", "#fcfcfb"
se = lambda v: np.std(v, ddof=1) / math.sqrt(len(v))

def per_replica(cond, proto, term):
    """Per-replica values. MT uses the pooled MOB1A reference, as the tables do."""
    d = pd.read_csv(f"{A}/mmpbsa_perframe_{cond}_{proto}.csv")
    out = {}
    if proto == "st":
        for s, _ in SYS:
            out[s] = d[d.system == s].groupby("replica")[term].mean().values
    else:
        pool = np.concatenate([d[d.system == s].groupby("replica")[f"rec_{term}"].mean().values
                               for s, _ in SYS]).mean()
        for s, _ in SYS:
            g = d[d.system == s]
            cpx = g.groupby("replica")[f"cpx_{term}"].mean().values
            lig = g.groupby("replica")[f"lig_{term}"].mean().values
            out[s] = cpx - pool - lig
    return out

def violin_panel(ax, data, title, ylab):
    for k, (cond, clab) in enumerate(COND):
        for j, (s, slab) in enumerate(SYS):
            v = data[cond][s]
            x = k + (j - 0.5) * 0.34
            parts = ax.violinplot([v], positions=[x], widths=0.28, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(C_SYS[s]); b.set_alpha(0.22); b.set_edgecolor("none")
            rng = np.random.default_rng(j * 7 + k)
            ax.scatter(x + rng.uniform(-0.055, 0.055, len(v)), v, s=5, color=C_SYS[s],
                       alpha=0.30, linewidths=0, zorder=2)
            ax.errorbar(x, v.mean(), yerr=se(v), fmt="o", ms=6, color=C_SYS[s], ecolor=INK,
                        elinewidth=1.4, capsize=3.5, capthick=1.4, zorder=3,
                        markeredgecolor=SURFACE, markeredgewidth=1.1,
                        label=slab if k == 0 else None)
    ax.set_xticks(range(len(COND))); ax.set_xticklabels([c[1] for c in COND], fontsize=9, color=INK)
    ax.set_xlim(-0.55, len(COND) - 0.45)
    ax.set_ylabel(ylab, fontsize=9.5, color=INK2)
    ax.set_title(title, fontsize=10, color=INK, loc="left", pad=6)

fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
fig.patch.set_alpha(0.0)
for ax in axes.ravel():
    ax.patch.set_alpha(0.0)
    ax.grid(True, axis="y", color="#e8e8e4", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=9, length=3, width=0.8)

st = {c: per_replica(c, "st", "TOTAL") for c, _ in COND}
mt = {c: per_replica(c, "mt", "TOTAL") for c, _ in COND}
cov = {}
for c, _ in COND:
    d = pd.read_csv(f"{A}/mmpbsa_perframe_{c}_mt.csv")
    terms = ["BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"]
    pool = np.concatenate([sum(d[d.system == s].groupby("replica")[f"rec_{t}"].mean().values
                               for t in terms) for s, _ in SYS]).mean()
    cov[c] = {}
    for s, _ in SYS:
        g = d[d.system == s]
        cpx = sum(g.groupby("replica")[f"cpx_{t}"].mean().values for t in terms)
        lig = sum(g.groupby("replica")[f"lig_{t}"].mean().values for t in terms)
        cov[c][s] = cpx - pool - lig

violin_panel(axes[0, 0], st,  r"(A)  Single-trajectory $\Delta H_{\mathrm{bind}}$", "kcal/mol")
violin_panel(axes[0, 1], mt,  r"(B)  Multi-trajectory $\Delta H_{\mathrm{bind}}$",  "kcal/mol")
violin_panel(axes[1, 0], cov, r"(C)  Multi-trajectory covalent term",               "kcal/mol")

# ---- Panel D: bootstrap of the LATS1 - NDR2 difference -----------------------------------
ax = axes[1, 1]
rng = np.random.default_rng(0)
NB = 4000
print("bootstrap ddH (LATS1 - NDR2), 4000 replica-level resamples")
for k, (cond, clab) in enumerate(COND):
    for j, (proto, src) in enumerate((("ST", st), ("MT", mt))):
        a, b = src[cond]["lats1"], src[cond]["ndr2"]
        boot = np.array([rng.choice(a, len(a), replace=True).mean()
                         - rng.choice(b, len(b), replace=True).mean() for _ in range(NB)])
        x = k + (j - 0.5) * 0.34
        parts = ax.violinplot([boot], positions=[x], widths=0.28, showextrema=False)
        for bb in parts["bodies"]:
            bb.set_facecolor(C_PROT[proto]); bb.set_alpha(0.30); bb.set_edgecolor("none")
        lo, hi = np.percentile(boot, [2.5, 97.5])
        ax.plot([x, x], [lo, hi], color=C_PROT[proto], lw=1.6, solid_capstyle="round", zorder=3)
        ax.plot([x], [boot.mean()], "o", ms=6, color=C_PROT[proto], zorder=4,
                markeredgecolor=SURFACE, markeredgewidth=1.1,
                label=proto if k == 0 else None)
        print(f"  {clab.replace(chr(10),' '):16s} {proto}: {boot.mean():+6.2f}  95% CI [{lo:+.2f}, {hi:+.2f}]")
ax.axhline(0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
ax.set_xticks(range(len(COND))); ax.set_xticklabels([c[1] for c in COND], fontsize=9, color=INK)
ax.set_xlim(-0.55, len(COND) - 0.45)
ax.set_ylabel("kcal/mol", fontsize=9.5, color=INK2)
ax.set_title(r"(D)  Bootstrap $\Delta\Delta H_{\mathrm{bind}}$ (LATS1 $-$ NDR2)",
             fontsize=10, color=INK, loc="left", pad=6)

for a_, loc in ((axes[0, 0], "lower left"), (axes[1, 1], "upper left")):
    leg = a_.legend(loc=loc, frameon=False, fontsize=9)
    for t in leg.get_texts(): t.set_color(INK2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/replica_distributions_n100.{ext}", dpi=300, transparent=True,
                bbox_inches="tight")
print(f"\nwrote {OUT}/replica_distributions_n100.png and .pdf")
