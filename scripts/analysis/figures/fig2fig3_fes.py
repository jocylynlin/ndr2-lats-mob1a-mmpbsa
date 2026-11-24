#!/usr/bin/env python3
"""Main-text Figures 2 and 3 -- conformational free-energy surfaces, final dataset.

Grid, contour levels, colours and layout are taken unchanged from
revision_analysis/.../scripts/fes_physcv.py, the script behind the submitted figures:
one GLOBAL (angle, Rg) grid built from every condition and state, so all eight panels of
both figures are directly comparable, 70x70 bins, Gaussian smoothing sigma = 1 bin,
F = -ln P in kT, contours every 0.5 kT to 4 kT.
Unbound CVs come from the matched-salt trajectories; bound CVs are unchanged.
"""
import os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from physcv_data import all_data, CONDS, SYS

OUT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
assert os.path.exists(os.path.join(OUT,'template.tex')), f'wrong output dir: {OUT}'
BOUND='#0072B2'; UNB='#D55E00'
D={k:(v[0].ravel(), v[1].ravel()) for k,v in all_data().items()}

AA=np.concatenate([v[0] for v in D.values()]); GG=np.concatenate([v[1] for v in D.values()])
xlo,xhi=np.percentile(AA,0.05),np.percentile(AA,99.95)
ylo,yhi=np.percentile(GG,0.05),np.percentile(GG,99.95)
xe=np.linspace(xlo,xhi,70); ye=np.linspace(ylo,yhi,70)
xc,yc=(xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2

def fes(a,g):
    H,_,_=np.histogram2d(a,g,bins=[xe,ye]); H=gaussian_filter(H,1.0)
    Pr=H/H.sum(); F=np.full(Pr.shape,np.nan)
    m=Pr>0; F[m]=-np.log(Pr[m]); F-=np.nanmin(F)
    return F.T
LV=np.arange(0.5,4.01,0.5)

def build(cols,fname,tag):
    fig,axes=plt.subplots(2,2,figsize=(11,9.2),dpi=300,sharex=True,sharey=True)
    for i,s in enumerate(SYS):
        for j,(cond,title) in enumerate(cols):
            ax=axes[i,j]
            for st,col in (('unbound',UNB),('bound',BOUND)):
                a,g=D[(cond,s,st)]
                ax.contour(xc,yc,fes(a,g),levels=LV,colors=col,linewidths=0.9,alpha=0.95)
            ax.grid(True,color='#ededed',lw=0.6); ax.set_axisbelow(True)
            for sp in ('top','right'): ax.spines[sp].set_visible(False)
            for sp in ('left','bottom'): ax.spines[sp].set_color('#444'); ax.spines[sp].set_linewidth(0.9)
            ax.tick_params(labelsize=10,color='#444')
            if i==0: ax.set_title(title,fontsize=13,fontweight='bold',pad=10)
            if i==1: ax.set_xlabel('Interhelical angle (deg)',fontsize=11.5)
            if j==0:
                ax.set_ylabel('$R_g$ ($\\AA$)',fontsize=11.5)
                ax.annotate({'lats1':'LATS1 segment','ndr2':'NDR2 segment'}[s],xy=(-0.19,0.5),
                            xycoords='axes fraction',rotation=90,ha='center',va='center',
                            fontsize=13,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=BOUND,lw=2,label='Bound'),
                        Line2D([],[],color=UNB,lw=2,label='Unbound')],
               loc='lower center',ncol=2,frameon=False,fontsize=12,bbox_to_anchor=(0.5,0.008))
    fig.tight_layout(rect=[0.025,0.045,1,1])
    for ext in ('png','pdf'): fig.savefig(f'{OUT}/{fname}.{ext}',bbox_inches='tight')
    plt.close(fig); print(f'  wrote {fname}  ({tag})')

build([('opc-salt','OPC'),('tip3p','TIP3P')],'fes_physcv_water','water model')
build([('opc-salt','Added NaCl'),('counterion','Counterion only')],'fes_physcv_salt','ionic condition')
print(f'  shared axes: angle {xlo:.1f}-{xhi:.1f} deg, Rg {ylo:.2f}-{yhi:.2f} A  (all 8 panels)')
