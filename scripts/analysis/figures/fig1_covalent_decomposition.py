#!/usr/bin/env python3
"""Figure 1 -- covalent decomposition, regenerated from the final matched-salt dataset.

Data   campaign_saltmatched/mmpbsa/opc-salt-matched/mt/{lats1,ndr2}_NNN/per_frame.csv
       100 replicas per system; unbound MOB1A pooled over both systems (200 replicas).
Checks the five plotted components sum to the covalent total and reproduce Table S1
(SI tab:covalent_decomp, Reference block) before writing anything.

Palette: categorical slots 1 and 2 of the validated default set (CVD dE 24.7,
normal-vision dE 33.6 -- all six checks pass). validate_palette.js cannot be executed
on this host (node needs GLIBC >= 2.28; login node is CentOS 7), so the pair is reused
unmodified from the pre-validated set rather than re-derived.
Layout follows ray_comment/scripts/plot_covalent_decomp_n100.py, the submitted figure's
script, so the revised figure is visually continuous with the rest of the manuscript.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from covalent_components import components, COV

OUT  = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
assert os.path.exists(os.path.join(OUT,'template.tex')), f'wrong output dir: {OUT}'
LATS1, NDR2 = "#2a78d6", "#eb6834"
INK, INK2, MUTED, SURFACE, GRID = "#0b0b0b", "#52514e", "#8a8a85", "#fcfcfb", "#e8e8e4"

out, diff, nl, nn = components("opc-salt-matched", "mt")
assert nl == nn == 100, f"expected 100 replicas per system, got {nl}/{nn}"
for s in ("lats1", "ndr2"):
    parts, tot = sum(out[(s, t)][0] for t in COV), out[(s, "Sum")][0]
    assert abs(parts - tot) < 1e-9, f"{s}: components {parts} do not close on total {tot}"
print("closure verified; components reconstruct the covalent total to < 1e-9 kcal/mol")
for s, lab in (("lats1", "LATS1"), ("ndr2", "NDR2")):
    print(f"  {lab:5s} " + "  ".join(f"{t} {out[(s,t)][0]:+.2f}" for t in COV)
          + f"   total {out[(s,'Sum')][0]:+.2f}")
print("  LATS1-NDR2 " + "  ".join(f"{t} {diff[t][0]:+.2f}+-{diff[t][1]:.2f}" for t in COV)
      + f"   total {diff['Sum'][0]:+.2f}+-{diff['Sum'][1]:.2f}")

x, w = np.arange(len(COV)), 0.36
fig, ax = plt.subplots(figsize=(7.4, 5.0))
fig.patch.set_alpha(0.0); ax.patch.set_alpha(0.0)
for off, s, col, lab in [(-w/2, "lats1", LATS1, "LATS1–MOB1A"),
                         ( w/2, "ndr2",  NDR2,  "NDR2–MOB1A")]:
    y = [out[(s, t)][0] for t in COV]; e = [out[(s, t)][1] for t in COV]
    ax.bar(x + off, y, w, color=col, label=lab, zorder=2,
           edgecolor=SURFACE, linewidth=1.2)
    ax.errorbar(x + off, y, yerr=e, fmt="none", ecolor=INK, elinewidth=1.3,
                capsize=3.5, capthick=1.3, zorder=3)
# direct label on the one component that carries the between-system difference
ax.text(x[COV.index("DIHED")], 0.45,
        f"$\\Delta$ = {diff['DIHED'][0]:+.2f} $\\pm$ {diff['DIHED'][1]:.2f}",
        ha="center", va="bottom", fontsize=9.5, color=INK2, zorder=4)
ax.axhline(0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
ax.set_xticks(x); ax.set_xticklabels(COV, fontsize=10.5, color=INK)
ax.set_ylabel(r"Grand mean $\Delta$Component (kcal/mol)", fontsize=11, color=INK2)
ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.8)
ax.tick_params(colors=INK2, labelsize=10, length=3, width=0.8)
leg = ax.legend(loc="lower left", frameon=False, fontsize=10.5)
for t in leg.get_texts(): t.set_color(INK2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/covalent_decomposition_n100.{ext}", dpi=300,
                transparent=True, bbox_inches="tight")
print(f"wrote {OUT}/covalent_decomposition_n100.png and .pdf")
