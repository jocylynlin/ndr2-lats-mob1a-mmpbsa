#!/usr/bin/env python3
"""SI Table S4 -- sensitivity of the normal-mode entropy to snapshots per replica.

Same source and conventions as tableS3_entropy.py. Subsets are EVENLY SPACED across each
replica's 10 snapshots, as the caption states ("every subset spans the full sampled time
range"), via linspace over the available indices -- not the first k. The 10-snapshot row
is by construction the Table S3 total row.
Pooled unbound MOB1A reference; per-system SE includes it, the difference excludes it.
"""
import glob, re, math
import numpy as np
ROOT='/home3/jocelyn/projects/hippo/campaign_saltmatched/nmode/opc-salt-matched'
SCALE=310.15/298.15
sem=lambda x: float(np.std(np.asarray(x,float),ddof=1)/math.sqrt(len(x)))

def parse(fn):
    ent=open(fn).read().split('POISSON BOLTZMANN:')[0]; out={}
    for key,lab in [('Complex','cpx'),('Receptor','rec'),('Ligand','lig')]:
        m=re.search(key+r' Energy Terms\nFrame #,Translational,Rotational,Vibrational,Total\n((?:[-\d].*\n)+)',ent)
        if m is None: return None
        out[lab]=np.array([[float(x) for x in l.split(',')] for l in m.group(1).strip().split('\n')])
    return out

S={s:[] for s in ('lats1','ndr2')}
for s in S:
    for d in sorted(glob.glob(f'{ROOT}/{s}/traj_[0-9][0-9][0-9]')):
        o=parse(f'{d}/per_frame.csv')
        if o is None: continue
        n=min(len(o['cpx']),len(o['rec']),len(o['lig']))
        S[s].append({k:o[k][:n,4] for k in ('cpx','rec','lig')})     # col 4 = Total
print(f"replicas: lats1 {len(S['lats1'])}, ndr2 {len(S['ndr2'])}")

# single-snapshot and replica-mean spread, for the caption's SD sentence
for s in S:
    allsnap=np.concatenate([-SCALE*(r['cpx']-r['rec']-r['lig']) for r in S[s]])
    repmean=np.array([(-SCALE*(r['cpx']-r['rec']-r['lig'])).mean() for r in S[s]])
    print(f"  {s}: SD of single-snapshot values {allsnap.std(ddof=1):.1f}, "
          f"SD of replica means {repmean.std(ddof=1):.1f} kcal/mol")

print("\nSnapshots & LATS1--MOB1A & NDR2--MOB1A & LATS1 - NDR2")
for k in (2,4,6,8,10):
    x={}; rec={}; C={}; L={}
    for s in S:
        v=[]; rv=[]; cv=[]; lv=[]
        for r in S[s]:
            idx=np.unique(np.linspace(0,len(r['cpx'])-1,k).round().astype(int))
            v.append(-(r['cpx'][idx]-r['lig'][idx]).mean())
            cv.append(r['cpx'][idx].mean()); lv.append(r['lig'][idx].mean())
            rv.append(r['rec'][idx].mean())
        x[s]=np.array(v); rec[s]=np.array(rv); C[s]=np.array(cv); L[s]=np.array(lv)
    pool=np.concatenate([rec['lats1'],rec['ndr2']])
    cells=[]
    for s in ('lats1','ndr2'):
        m=SCALE*(x[s].mean()+pool.mean()); e=SCALE*math.hypot(sem(C[s]),sem(L[s]),sem(pool))
        cells.append(f"\\saltvalue{{${m:.2f} \\pm {e:.2f}$}}")
    dm=SCALE*(x['lats1'].mean()-x['ndr2'].mean()); de=SCALE*math.hypot(sem(C['lats1']),sem(L['lats1']),sem(C['ndr2']),sem(L['ndr2']))
    cells.append(f"\\saltvalue{{${dm:+.2f} \\pm {de:.2f}$}}")
    print(f"{k:<2d} & " + " & ".join(cells) + r" \\")
print("\ncurrently printed: 2 65.51/64.95/+0.57 | 4 66.67/65.42/+1.25 | 6 67.20/64.01/+3.19 |"
      " 8 67.31/64.94/+2.37 | 10 67.27/64.63/+2.64")
print("Table S3 total row (must equal the k=10 row): 64.94+-0.78 / 62.35+-0.79 / +2.59+-1.01")
