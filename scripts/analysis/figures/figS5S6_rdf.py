#!/usr/bin/env python3
"""SI Figures S5 and S6 -- hydration and ionic atmosphere, final dataset.

Layout, palette, panel order and normalization logic are taken unchanged from
ray_comment/scripts/plot_c2_rdf_si.py, the script behind the submitted figures.
Two things change.

  Data   analysis_hydration_matched/data/{cond}/{kinase}/replicaNNN -- the matched-salt
         solvent-retaining set, 20 replicas per condition and kinase instead of 5.
  Ions   the actual ion content of the matched boxes, read from verify_salt.py rather
         than the fixed 80 Na+ / 86 Cl- of the superseded boxes.

Normalization: cpptraj `radial` was run without `density`, so it divides by its default
0.033456 molecules/Ang^3 (bulk water). That is right for the water panels and wrong for the
ion panels, where the reference is the ion number density N_ion/V. Since rho_water ~ N_wat/V,
the correction is the constant N_wat/N_ion. It is the same for OPC and TIP3P within each
kinase to 0.2%, so it does not affect TIP3P/OPC peak-height ratios.
"""
import glob, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home3/jocelyn/projects/hippo/reviewer_comments/analysis_hydration_matched/data"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                   "Supplemental Information")
assert os.path.exists(os.path.join(OUT, "manuscript_SI_polished.tex")), f"wrong output dir: {OUT}"
OPC, TIP3P = "#2a78d6", "#eb6834"
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#8a8a85", "#fcfcfb"

# Matched-salt box composition (campaign_saltmatched/verify_salt.py).
NWAT = {("opc-salt-matched","lats1_kinase"): 9002, ("tip3p-matched","lats1_kinase"): 9016,
        ("opc-salt-matched","ndr2_kinase"):  8693, ("tip3p-matched","ndr2_kinase"):  8706}
NION = {("lats1_kinase","rdf_ani_Na"): 52, ("lats1_kinase","rdf_cat_Cl"): 58,
        ("ndr2_kinase","rdf_ani_Na"):  51, ("ndr2_kinase","rdf_cat_Cl"):  57}

GROUPS = {
    "rdf_water": [("rdf_ani_watO", "Water O around\ncarboxylate O", (1.5, 6.0), True),
                  ("rdf_cat_watO", "Water O around\nArg/Lys N",     (1.5, 6.0), True)],
    "rdf_ions":  [("rdf_ani_Na",   "Na$^+$ around\ncarboxylate O",  (1.5, 8.0), True),
                  ("rdf_cat_Cl",   "Cl$^-$ around\nArg/Lys N",      (2.0, 8.0), True)],
}
SYSTEMS = [("lats1_kinase", "LATS1 kinase segment"), ("ndr2_kinase", "NDR2 kinase segment")]
CONDS = [("opc-salt-matched", OPC, "OPC"), ("tip3p-matched", TIP3P, "TIP3P")]

def ion_scale(cond, sysn, tag):
    return NWAT[(cond,sysn)]/NION[(sysn,tag)] if (sysn,tag) in NION else 1.0

def curves(cond, sysn, tag):
    g=[]; r=None
    for f in sorted(glob.glob(f"{A}/{cond}/{sysn}/replica*/{tag}.dat")):
        d=np.loadtxt(f, comments="#")
        if d.ndim!=2 or d.shape[1]<2: continue
        r=d[:,0]; g.append(d[:,1])
    if not g: return None,None,None,0
    g=np.vstack(g)*ion_scale(cond,sysn,tag)
    return r, g.mean(0), g.std(0,ddof=1)/np.sqrt(len(g)), len(g)

def peak(cond,sysn,tag,lo,hi):
    """first-peak height and position, with replica-level SE on the height."""
    gs=[]; r=None
    for f in sorted(glob.glob(f"{A}/{cond}/{sysn}/replica*/{tag}.dat")):
        d=np.loadtxt(f,comments="#")
        if d.ndim!=2 or d.shape[1]<2: continue
        r=d[:,0]; gs.append(d[:,1])
    if not gs: return None
    G=np.vstack(gs)*ion_scale(cond,sysn,tag)
    m=G.mean(0); w=(r>=lo)&(r<=hi); k=np.argmax(np.where(w,m,-np.inf))
    h=G[:,k]
    return float(m[k]), float(np.std(h,ddof=1)/math.sqrt(len(h))), float(r[k]), len(gs)

if __name__=="__main__":
    print("replica counts and first peaks (final dataset)")
    for sysn,_ in SYSTEMS:
        for cond,_,name in CONDS:
            for tag,lo,hi in (("rdf_ani_watO",2.4,3.2),("rdf_ani_Na",2.0,3.0),("rdf_cat_Cl",3.0,4.2)):
                p=peak(cond,sysn,tag,lo,hi)
                if p: print(f"  {sysn:13s} {name:6s} {tag:14s} peak {p[0]:7.3f}+-{p[1]:.3f} at {p[2]:.2f} A  (n={p[3]})")
    print("\nTIP3P/OPC first-peak ratios")
    for sysn,lab in SYSTEMS:
        for tag,lo,hi in (("rdf_ani_Na",2.0,3.0),):
            o=peak("opc-salt-matched",sysn,tag,lo,hi); t=peak("tip3p-matched",sysn,tag,lo,hi)
            if o and t:
                ratio=t[0]/o[0]; se=ratio*math.hypot(t[1]/t[0],o[1]/o[0])
                print(f"  {lab:13s} {tag}: TIP3P/OPC {ratio:.2f}+-{se:.2f}  "
                      f"(OPC peak at {o[2]:.2f} A, TIP3P at {t[2]:.2f} A)")
    for sysn,lab in SYSTEMS:
        o=peak("opc-salt-matched",sysn,"rdf_cat_Cl",3.0,4.2); t=peak("tip3p-matched",sysn,"rdf_cat_Cl",3.0,4.2)
        if o and t: print(f"  {lab:13s} Cl- first-peak position: OPC {o[2]:.2f} A, TIP3P {t[2]:.2f} A, shift {t[2]-o[2]:+.2f} A")

    for gname,PANELS in GROUPS.items():
        fig,axes=plt.subplots(len(PANELS),2,figsize=(6.5,5.4),sharex=False)
        fig.patch.set_alpha(0.0)
        for i,(sysn,syslab) in enumerate(SYSTEMS):
            for j,(tag,panlab,xlim,bulkline) in enumerate(PANELS):
                ax=axes[j,i]; ax.patch.set_alpha(0.0)
                for cond,col,name in CONDS:
                    r,m,s,n=curves(cond,sysn,tag)
                    if r is None: continue
                    ax.fill_between(r,m-s,m+s,color=col,alpha=0.20,linewidth=0)
                    ax.plot(r,m,color=col,lw=1.8,solid_capstyle="round",label=name)
                ax.set_xlim(*xlim)
                if bulkline:
                    ax.axhline(1.0,color=MUTED,lw=0.7,ls=(0,(4,3)),zorder=0)
                    ax.text(xlim[1]-0.08,1.0,"bulk",fontsize=7.5,color=MUTED,va="bottom",ha="right")
                ax.grid(True,color="#e8e8e4",lw=0.6,zorder=0); ax.set_axisbelow(True)
                for sp in ("top","right"): ax.spines[sp].set_visible(False)
                for sp in ("left","bottom"):
                    ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.8)
                ax.tick_params(colors=INK2,labelsize=9,length=3,width=0.8)
                if j==0: ax.set_title(syslab,fontsize=10.5,color=INK,pad=8)
                if j==len(PANELS)-1: ax.set_xlabel(r"distance ($\mathrm{\AA}$)",fontsize=9.5,color=INK2)
                ax.set_ylabel(f"{panlab}\n\n$g(r)$",fontsize=9,color=INK) if i==0 else ax.set_ylabel("")
        ax0=axes[0,0]
        r0,m0,_,_=curves("opc-salt-matched","lats1_kinase",PANELS[0][0])
        k=int(np.nanargmax(m0))
        for cond,col,name,dy in [("opc-salt-matched",OPC,"OPC",0.14),("tip3p-matched",TIP3P,"TIP3P",-0.04)]:
            r,m,_,_=curves(cond,"lats1_kinase",PANELS[0][0])
            ax0.annotate(name,xy=(r[k],m[k]),xytext=(r[k]+0.9,m[k]*(0.80+dy)),
                         fontsize=9,color=col,fontweight="semibold",va="center",ha="left")
        h,l=axes[0,0].get_legend_handles_labels()
        leg=fig.legend(h,l,loc="lower center",ncol=2,frameon=False,fontsize=10,
                       bbox_to_anchor=(0.5,-0.005),handlelength=1.8)
        for t in leg.get_texts(): t.set_color(INK2)
        fig.tight_layout(rect=(0,0.055,1,1))
        for ext in ("png","pdf"):
            fig.savefig(f"{OUT}/{gname}_si.{ext}",dpi=300,transparent=True,bbox_inches="tight")
        print(f"wrote {OUT}/{gname}_si.png")
        plt.close(fig)
