import pandas as pd, numpy as np, glob, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
P='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/physcv/'
OUT='/home3/jocelyn/projects/hippo/reviewer_comments/revision_analysis/results/campaign_n100/'
BOUND='#0072B2'; UNB='#D55E00'
import os
_here = os.path.dirname(os.path.abspath(__file__))
_cand = [os.path.join(_here, 'cv_metrics.json'),
         os.path.join(os.path.dirname(_here), 'cv_metrics.json'),
         'cv_metrics.json']
_mf = next(f for f in _cand if os.path.exists(f))
M = json.load(open(_mf))
def load(c,s,st):
    fs=sorted(glob.glob(f'{P}{c}/{st}/{s}_*.csv')); d=[pd.read_csv(f) for f in fs]
    return (np.concatenate([x['angle'].values for x in d]),
            np.concatenate([x['rg'].values for x in d]))
CONDS=['opc-salt','tip3p','counterion']; SYS=['lats1','ndr2']
D={(c,s,st):load(c,s,st) for c in CONDS for s in SYS for st in ['bound','unbound']}
# one global grid so every panel of both figures is directly comparable
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
def build(cols,fname,caption_tag):
    fig,axes=plt.subplots(2,2,figsize=(11,9.2),dpi=300,sharex=True,sharey=True)
    for i,s in enumerate(SYS):
        for j,(cond,title) in enumerate(cols):
            ax=axes[i,j]
            for st,col in [('unbound',UNB),('bound',BOUND)]:
                a,g=D[(cond,s,st)]
                ax.contour(xc,yc,fes(a,g),levels=LV,colors=col,linewidths=0.9,alpha=0.95)
            # quantitative metrics are tabulated in the Supporting Information, not on the panels
            ax.grid(True,color='#ededed',lw=0.6); ax.set_axisbelow(True)
            for sp in ('top','right'): ax.spines[sp].set_visible(False)
            for sp in ('left','bottom'): ax.spines[sp].set_color('#444'); ax.spines[sp].set_linewidth(0.9)
            ax.tick_params(labelsize=10,color='#444')
            if i==0: ax.set_title(title,fontsize=13,fontweight='bold',pad=10)
            if i==1: ax.set_xlabel('Interhelical angle (deg)',fontsize=11.5)
            if j==0:
                ax.set_ylabel('$R_g$ ($\\AA$)',fontsize=11.5)
                ax.annotate({'lats1':'LATS1-MBD','ndr2':'NDR2-MBD'}[s],xy=(-0.19,0.5),
                            xycoords='axes fraction',rotation=90,ha='center',va='center',
                            fontsize=13,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=BOUND,lw=2,label='Bound'),
                        Line2D([],[],color=UNB,lw=2,label='Unbound')],
               loc='lower center',ncol=2,frameon=False,fontsize=12,bbox_to_anchor=(0.5,0.008))
    fig.tight_layout(rect=[0.025,0.045,1,1])
    for ext in ('png','pdf'): fig.savefig(f'{OUT}{fname}.{ext}',bbox_inches='tight')
    plt.close(fig); print(f'  wrote {fname}  ({caption_tag})')
build([('opc-salt','OPC'),('tip3p','TIP3P')],'fes_physcv_water','water model')
build([('opc-salt','Added NaCl (80 ion pairs)'),('counterion','Counterion only')],'fes_physcv_salt','ionic condition')
print(f'  shared axes: angle {xlo:.1f}-{xhi:.1f} deg, Rg {ylo:.2f}-{yhi:.2f} A  (all 8 panels)')
