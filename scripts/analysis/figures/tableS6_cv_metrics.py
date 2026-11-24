#!/usr/bin/env python3
"""SI Table S6 -- conformational-distribution metrics, final dataset.

Definitions are byte-identical to revision_analysis/.../scripts/cv_metrics.py:
  shared 61x61 (angle, Rg) grid PER SYSTEM, spanning every condition and state, so values
  are comparable across conditions;
  Area    = fraction of occupied grid cells within 3 kT of the free-energy minimum
            (a dimensionless grid fraction, not a physical area);
  Overlap = sum_i min(p_i, q_i) between the bound and unbound distributions;
  JS      = Jensen-Shannon divergence in bits;
  uncertainties = replica-level bootstrap SE over 200 resamples of the 100 replicas,
            bound and unbound resampled independently.
"""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from physcv_data import all_data, CONDS, SYS

rng=np.random.default_rng(0)
D=all_data()

def grid(s):
    A=np.concatenate([D[(c,s,st)][0].ravel() for c in CONDS for st in ('bound','unbound')])
    G=np.concatenate([D[(c,s,st)][1].ravel() for c in CONDS for st in ('bound','unbound')])
    return np.linspace(A.min(),A.max(),61), np.linspace(G.min(),G.max(),61)
GR={s:grid(s) for s in SYS}

def hist(a,g,s):
    xe,ye=GR[s]
    H,_,_=np.histogram2d(a.ravel(),g.ravel(),bins=[xe,ye])
    return H/H.sum()

def metrics(Hb,Hu):
    p,q=Hb.ravel(),Hu.ravel()
    ov=np.minimum(p,q).sum()
    m=0.5*(p+q)
    nz=lambda x,y: np.sum(x[(x>0)&(y>0)]*np.log2(x[(x>0)&(y>0)]/y[(x>0)&(y>0)]))
    js=0.5*nz(p,m)+0.5*nz(q,m)
    def area(H):
        h=H[H>0]; F=-np.log(h); return (F<F.min()+3).sum()/H.size
    return dict(ab=area(Hb),au=area(Hu),ov=ov,js=js)

LBL={'lats1':'LATS1-MBD','ndr2':'NDR2-MBD'}
TITLE={'opc-salt':'Reference (OPC, added NaCl)','tip3p':'TIP3P','counterion':'Counterion only'}
for c in CONDS:
    print(f"\n% --- {TITLE[c]} ---")
    for s in SYS:
        ab,gb=D[(c,s,'bound')]; au,gu=D[(c,s,'unbound')]
        pt=metrics(hist(ab,gb,s),hist(au,gu,s))
        bs=[]
        for _ in range(200):
            i=rng.integers(0,len(ab),len(ab)); j=rng.integers(0,len(au),len(au))
            bs.append(metrics(hist(ab[i],gb[i],s),hist(au[j],gu[j],s)))
        se={k:float(np.std([b[k] for b in bs],ddof=1)) for k in pt}
        cells=[f"\\saltvalue{{${pt[k]:.3f} \\pm {se[k]:.3f}$}}" for k in ('ab','au','ov','js')]
        print(f"{LBL[s]} & " + " & ".join(cells) + r" \\")
