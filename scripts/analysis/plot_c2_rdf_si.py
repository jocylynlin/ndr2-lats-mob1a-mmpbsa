#!/usr/bin/env python3
"""plot_c2_rdf.py -- RDF figure for the C2 solvent-retaining set (Reviewer 1 Comment 3).

Form: radial distribution functions are a continuous function of distance, so lines.
Two series only (OPC, TIP3P) -> categorical slots 1 and 2, validated with
scratchpad/validate_palette.py (CVD dE 24.7, normal-vision dE 33.6, all six checks PASS).
Bands are +/- 1 SE across the 5 replicas -- replica-level, never frame-pooled.
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/c2_hydration"
OUT = "/home3/jocelyn/projects/hippo/reviewer_comments/manuscript_rev/Supplemental Information"

OPC, TIP3P = "#2a78d6", "#eb6834"          # categorical slots 1 and 2
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
SURFACE = "#fcfcfb"

# (tag, panel title, x range, draw the g=1 bulk reference line)
# The ion panels peak near 0.16, so forcing g=1 onto their y-axis would flatten the
# very difference the panel exists to show. They autoscale and carry a note instead.
GROUPS = {
    "rdf_water": [
        ("rdf_ani_watO", "Water O around\ncarboxylate O", (1.5, 6.0), True),
        ("rdf_cat_watO", "Water O around\nArg/Lys N",     (1.5, 6.0), True),
    ],
    "rdf_ions": [
        ("rdf_ani_Na",   "Na$^+$ around\ncarboxylate O",  (1.5, 8.0), True),
        ("rdf_cat_Cl",   "Cl$^-$ around\nArg/Lys N",      (2.0, 8.0), True),
    ],
}
SYSTEMS = [("lats1_kinase", "LATS1 kinase"), ("ndr2_kinase", "NDR2 kinase")]



# cpptraj's `radial` was run without `density`, so it normalised every RDF by its default
# 0.033456 molecules/Ang^3 (bulk water). That is right for the water panels and wrong for the
# ion panels, where the correct reference is the ion number density N_ion/V. Since
# rho_water ~ N_wat/V, the correction is the pure constant N_wat/N_ion. It is identical for OPC
# and TIP3P to within 0.1%, so it does not touch the TIP3P/OPC ratio - only the absolute scale.
NWAT = {("opc-salt", "lats1_kinase"): 8970, ("tip3p", "lats1_kinase"): 8980,
        ("opc-salt", "ndr2_kinase"):  8653, ("tip3p", "ndr2_kinase"):  8664}
NION = {"rdf_ani_Na": 80, "rdf_cat_Cl": 86}

def ion_scale(cond, sysn, tag):
    return NWAT[(cond, sysn)] / NION[tag] if tag in NION else 1.0

def curves(cond, sysn, tag):
    """-> r, mean g(r), SE g(r) over replicas."""
    g = []
    r = None
    for f in sorted(glob.glob(f"{A}/{cond}/{sysn}/replica*/{tag}.dat")):
        d = np.loadtxt(f, comments="#")
        if d.ndim != 2 or d.shape[1] < 2:
            continue
        r = d[:, 0]
        g.append(d[:, 1])
    if not g:
        return None, None, None
    g = np.vstack(g) * ion_scale(cond, sysn, tag)
    return r, g.mean(0), g.std(0, ddof=1) / np.sqrt(len(g))



for gname, PANELS in GROUPS.items():
    fig, axes = plt.subplots(len(PANELS), 2, figsize=(6.5, 5.4), sharex=False)
    fig.patch.set_alpha(0.0)

    fig.patch.set_alpha(0.0)

    for i, (sysn, syslab) in enumerate(SYSTEMS):
        for j, (tag, panlab, xlim, bulkline) in enumerate(PANELS):
            ax = axes[j, i]        # rows = RDF pair, cols = system
            ax.patch.set_alpha(0.0)
            for cond, col, name in [("opc-salt", OPC, "OPC"), ("tip3p", TIP3P, "TIP3P")]:
                r, m, s = curves(cond, sysn, tag)
                if r is None:
                    continue
                ax.fill_between(r, m - s, m + s, color=col, alpha=0.20, linewidth=0)
                ax.plot(r, m, color=col, lw=1.8, solid_capstyle="round", label=name)
            ax.set_xlim(*xlim)
            if bulkline:
                ax.axhline(1.0, color=MUTED, lw=0.7, ls=(0, (4, 3)), zorder=0)
                ax.text(xlim[1] - 0.08, 1.0, "bulk", fontsize=7.5, color=MUTED,
                        va="bottom", ha="right")
            else:
                ax.set_ylim(bottom=0)
                ax.margins(y=0.14)
            ax.grid(True, color="#e8e8e4", lw=0.6, zorder=0)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(MUTED)
                ax.spines[sp].set_linewidth(0.8)
            ax.tick_params(colors=INK2, labelsize=9, length=3, width=0.8)
            if j == 0:
                ax.set_title(syslab, fontsize=10.5, color=INK, pad=8)
            if not bulkline and i == 1:
                ax.text(0.97, 0.94, "bulk $g=1$ off scale", transform=ax.transAxes,
                        fontsize=7.5, color=MUTED, ha="right", va="top")
            if j == len(PANELS) - 1:
                ax.set_xlabel(r"distance ($\mathrm{\AA}$)", fontsize=9.5, color=INK2)
            if i == 0:
                ax.set_ylabel(f"{panlab}\n\n$g(r)$", fontsize=9, color=INK)
            else:
                ax.set_ylabel("")

    # direct-label the two series on the first panel of this group
    ax0 = axes[0, 0]
    r0, m0, _ = curves("opc-salt", "lats1_kinase", PANELS[0][0])
    k = int(np.nanargmax(m0))
    for cond, col, name, dy in [("opc-salt", OPC, "OPC", 0.14), ("tip3p", TIP3P, "TIP3P", -0.04)]:
        r, m, _ = curves(cond, "lats1_kinase", PANELS[0][0])
        ax0.annotate(name, xy=(r[k], m[k]), xytext=(r[k] + 0.9, m[k] * (0.80 + dy)),
                     fontsize=9, color=col, fontweight="semibold", va="center", ha="left")

    h, l = axes[0, 0].get_legend_handles_labels()
    leg = fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=10,
                     bbox_to_anchor=(0.5, -0.005), handlelength=1.8)
    for t in leg.get_texts():
        t.set_color(INK2)


    fig.tight_layout(rect=(0, 0.055, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{gname}_si.{ext}", dpi=300, transparent=True, bbox_inches="tight")
    print(f"wrote {OUT}/{gname}_si.png")
    plt.close(fig)
