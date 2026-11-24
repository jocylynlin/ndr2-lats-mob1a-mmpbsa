#!/usr/bin/env python3
"""Split-half null for the conformational free-energy surfaces.

The overlap and Jensen-Shannon values of Table S6 compare the BOUND and UNBOUND ensembles.
They do not say whether two CONDITIONS are distinguishable from one another, which is the
claim Figures 2 and 3 are cited for. This script answers that directly:

  real  = JS between the unbound ensembles of two conditions
  null  = JS between two random halves of the SAME condition's 100 replicas,
          i.e. the floor set by sampling noise alone (200 random splits)

A real value inside the null band means the two conditions are not distinguishable at this
sample size, however different their contours look. Same per-system 61x61 grid as Table S6.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from physcv_data import all_data, CONDS, SYS

NSPLIT = 200
rng = np.random.default_rng(0)
D = all_data()

def grid(s):
    A=np.concatenate([D[(c,s,st)][0].ravel() for c in CONDS for st in ('bound','unbound')])
    G=np.concatenate([D[(c,s,st)][1].ravel() for c in CONDS for st in ('bound','unbound')])
    return np.linspace(A.min(),A.max(),61), np.linspace(G.min(),G.max(),61)
GR={s:grid(s) for s in SYS}

def hist(a,g,s):
    xe,ye=GR[s]; H,_,_=np.histogram2d(a.ravel(),g.ravel(),bins=[xe,ye]); return H/H.sum()
def js(p,q):
    p,q=p.ravel(),q.ravel(); m=0.5*(p+q)
    nz=lambda x,y: np.sum(x[(x>0)&(y>0)]*np.log2(x[(x>0)&(y>0)]/y[(x>0)&(y>0)]))
    return 0.5*nz(p,m)+0.5*nz(q,m)
def ov(p,q): return np.minimum(p.ravel(),q.ravel()).sum()

def null(c,s,st):
    a,g=D[(c,s,st)]; N=len(a); o=[]; j=[]
    for _ in range(NSPLIT):
        i=rng.permutation(N); A,B=i[:N//2],i[N//2:]
        HA,HB=hist(a[A],g[A],s),hist(a[B],g[B],s)
        o.append(ov(HA,HB)); j.append(js(HA,HB))
    return np.array(o),np.array(j)

PAIRS=[('opc-salt','tip3p','OPC vs TIP3P'),('opc-salt','counterion','added NaCl vs counterion')]
if __name__=='__main__':
    print(f"Split-half null, {NSPLIT} random splits of the 100 replicas, unbound ensembles\n")
    for s in SYS:
        print(f"--- {s} ---")
        nj={}
        for c in CONDS:
            o,j=null(c,s,'unbound'); nj[c]=(j.mean(),j.std(ddof=1))
            print(f"   null {c:11s}  overlap {o.mean():.3f}+-{o.std(ddof=1):.3f}   JS {j.mean():.3f}+-{j.std(ddof=1):.3f}")
        for c1,c2,lab in PAIRS:
            H1=hist(*D[(c1,s,'unbound')],s); H2=hist(*D[(c2,s,'unbound')],s)
            r=js(H1,H2)
            z=min(abs(r-nj[c][0])/nj[c][1] for c in (c1,c2))
            print(f"   REAL {lab:26s} overlap {ov(H1,H2):.3f}  JS {r:.3f}   "
                  f"-> {z:.1f} SD above the nearer null")
        print()
    print("Bound-vs-unbound, for scale (the contrast the surfaces DO resolve):")
    for s in SYS:
        for c in CONDS:
            print(f"   {s:6s} {c:11s} JS {js(hist(*D[(c,s,'bound')],s),hist(*D[(c,s,'unbound')],s)):.3f}")
