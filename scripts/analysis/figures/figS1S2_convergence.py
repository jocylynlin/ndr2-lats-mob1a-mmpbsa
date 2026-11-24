#!/usr/bin/env python3
"""SI Figures S1 and S2 -- time-series stability of dH_bind and of the covalent term.

Layout, smoothing, CI definition, colours and row order are taken unchanged from
revision_analysis/.../scripts/convplot.py, the script behind the submitted figures.
Only the data source changes:

  Reference / Single-traj / TIP3P   ->  matched-salt per-frame data via frames_cache.py
  Counterion only                   ->  campaign_100x40ns tidied CSV, unchanged, because
                                        that condition contains no added NaCl and was
                                        not resimulated.

Multi-trajectory rows use the unbound MOB1A reference pooled over both systems
(200 replicas); its SE enters the per-system band in quadrature, as in Table 1.
"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE=os.path.dirname(os.path.abspath(__file__))
OUT=os.path.join(os.path.dirname(os.path.dirname(HERE)),'manuscript_working_20261122','Supplemental Information')
OUT=os.path.join(os.path.dirname(os.path.dirname(HERE)),'Supplemental Information')
assert os.path.exists(os.path.join(OUT,'manuscript_SI_polished.tex')), f'wrong output dir: {OUT}'
OLD='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
CLR='#7b2d5e'; FILL='#e3cfdc'
COV=['BOND','ANGLE','DIHED','UB','IMP','CMAP','1-4 VDW','1-4 EEL']
ROWS=[('Reference','opc-salt-matched','mt'),('Single-traj','opc-salt-matched','st'),
      ('TIP3P','tip3p-matched','mt'),('Counterion only','counterion','mt')]
SYS=[('lats1','LATS1–MOB1A'),('ndr2','NDR2–MOB1A')]
NF=400; T=np.arange(1,NF+1)*0.1

def series_matched(cond,prot,which):
    D=np.load(os.path.join(HERE,f'frames_{cond}_{prot}.npz'))
    q={'dH':'TOTAL','cov':'COV'}[which]
    out={}
    if prot=='st':
        for s,_ in SYS:
            X=D[f'cpx_{s}_{q}']-D[f'rec_{s}_{q}']-D[f'lig_{s}_{q}']
            out[s]=(T[:X.shape[1]], X.mean(0), X.std(0,ddof=1)/np.sqrt(X.shape[0]))
    else:
        R=np.vstack([D[f'rec_{s}_{q}'] for s,_ in SYS])      # pooled, 200 replicas
        rm=R.mean(0); rse=R.std(0,ddof=1)/np.sqrt(R.shape[0])
        for s,_ in SYS:
            Cx=D[f'cpx_{s}_{q}']; Lx=D[f'lig_{s}_{q}']; X=Cx-Lx
            sem=lambda A: A.std(0,ddof=1)/np.sqrt(A.shape[0])     # independent simulations: quadrature
            out[s]=(T[:X.shape[1]], X.mean(0)-rm, np.sqrt(sem(Cx)**2+sem(Lx)**2+rse**2))
    return out

def series_old(cond,prot,which):
    df=pd.read_csv(f'{OLD}mmpbsa_perframe_{cond}_{prot}.csv')
    for pre in ['','cpx_','rec_','lig_']:
        df[pre+'COV']=sum(df[pre+c] for c in COV)
    col={'dH':'TOTAL','cov':'COV'}[which]
    out={}
    rec=[df[df.system==s].pivot_table(index='replica',columns='time_ns',values='rec_'+col).values
         for s,_ in SYS]
    R=np.vstack(rec); rm=R.mean(0); rse=R.std(0,ddof=1)/np.sqrt(R.shape[0])
    for s,_ in SYS:
        d=df[df.system==s]
        Pc=d.pivot_table(index='replica',columns='time_ns',values='cpx_'+col)
        Pl=d.pivot_table(index='replica',columns='time_ns',values='lig_'+col)
        X=Pc.values-Pl.values
        sem=lambda A: A.std(0,ddof=1)/np.sqrt(A.shape[0])         # independent simulations: quadrature
        out[s]=(Pc.columns.values, X.mean(0)-rm, np.sqrt(sem(Pc.values)**2+sem(Pl.values)**2+rse**2))
    return out

series=lambda c,p,w: series_old(c,p,w) if c=='counterion' else series_matched(c,p,w)
roll=lambda y,w=20: pd.Series(y).rolling(w,center=True,min_periods=1).mean().values

def build(which,ylab,fname,tag):
    DATA={lab:series(c,p,which) for lab,c,p in ROWS}
    lo=min(np.nanmin(m-1.96*se) for d in DATA.values() for _,m,se in d.values())
    hi=max(np.nanmax(m+1.96*se) for d in DATA.values() for _,m,se in d.values())
    pad=0.06*(hi-lo); ylim=(lo-pad,hi+pad)
    fig,axes=plt.subplots(4,2,figsize=(9.9,12.6),dpi=300,sharex=True,sharey=True)
    for i,(lab,_,_) in enumerate(ROWS):
        for j,(s,title) in enumerate(SYS):
            ax=axes[i,j]; t,m,se=DATA[lab][s]
            ax.plot(t,m,color=CLR,lw=0.5,alpha=0.30,zorder=1)
            ax.fill_between(t,roll(m-1.96*se),roll(m+1.96*se),color=FILL,alpha=0.85,lw=0,zorder=2)
            ax.plot(t,roll(m),color=CLR,lw=2.0,zorder=4)
            ax.axhline(m.mean(),color=CLR,ls='--',lw=1.2,alpha=0.85,zorder=3)
            ax.text(0.975,0.94,f'{tag}: {m.mean():.2f}',transform=ax.transAxes,
                    ha='right',va='top',color=CLR,fontsize=11)
            ax.grid(True,color='#e5e5e5',lw=0.7,zorder=0); ax.set_axisbelow(True)
            for sp in ('top','right'): ax.spines[sp].set_visible(False)
            for sp in ('left','bottom'): ax.spines[sp].set_color('#444444'); ax.spines[sp].set_linewidth(0.9)
            ax.tick_params(labelsize=9,color='#444444')
            ax.set_ylim(*ylim); ax.set_xlim(t.min(),t.max())
            if i==0: ax.set_title(title,fontsize=13,fontweight='bold',pad=12)
            if j==0:
                ax.set_ylabel(ylab,fontsize=10)
                ax.annotate(lab,xy=(-0.30,0.5),xycoords='axes fraction',rotation=90,
                            ha='center',va='center',fontsize=13,fontweight='bold')
            if i==3: ax.set_xlabel('Simulation Time (ns)',fontsize=11)
    fig.legend(handles=[Line2D([],[],color=CLR,lw=2.0,label='Rolling mean (2 ns)'),
                        Patch(facecolor=FILL,label='95% CI (smoothed)'),
                        Line2D([],[],color=CLR,lw=0.8,alpha=0.45,label='Raw per-frame mean'),
                        Line2D([],[],color=CLR,lw=1.2,ls='--',label='Overall mean')],
               loc='lower center',ncol=4,frameon=False,fontsize=10.5,bbox_to_anchor=(0.5,0.005))
    fig.tight_layout(rect=[0.035,0.035,1,1])
    for ext in ('png','pdf'): fig.savefig(f'{OUT}/{fname}.{ext}',bbox_inches='tight')
    plt.close(fig)
    print(f'  {fname}: ylim ({ylim[0]:.1f}, {ylim[1]:.1f})')
    for lab,_,_ in ROWS:
        print(f'    {lab:16s} '+'  '.join(f'{t2} {DATA[lab][s][1].mean():8.2f}' for s,t2 in SYS))

build('dH', r'$\Delta H_{\mathrm{bind}}$  (kcal/mol)','convergence_dHbind_n100','Mean $\\Delta$H')
build('cov',r'$\Delta E_{\mathrm{cov}}$  (kcal/mol)','covalent_convergence_n100','Mean $\\Delta E$')
