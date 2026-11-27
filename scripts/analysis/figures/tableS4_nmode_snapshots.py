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
        S[s].append({k:o[k][:,4] for k in ('cpx','rec','lig')})      # col 4 = Total; all converged snapshots
print(f"replicas: lats1 {len(S['lats1'])}, ndr2 {len(S['ndr2'])}")

# single-snapshot and replica-mean spread, for the caption's SD sentence
for s in S:
    # descriptive only: snapshots combined by index, over the indices all three species have
    def comb(r):
        n=min(len(r['cpx']),len(r['rec']),len(r['lig'])); return -SCALE*(r['cpx'][:n]-r['rec'][:n]-r['lig'][:n])
    allsnap=np.concatenate([comb(r) for r in S[s]])
    repmean=np.array([-SCALE*(r['cpx'].mean()-r['rec'].mean()-r['lig'].mean()) for r in S[s]])
    print(f"  {s}: SD of single-snapshot values {allsnap.std(ddof=1):.1f}, "
          f"SD of replica means {repmean.std(ddof=1):.1f} kcal/mol")

print("\nSnapshots & LATS1--MOB1A & NDR2--MOB1A & LATS1 - NDR2")
for k in (2,4,6,8,10):
    x={}; rec={}; C={}; L={}
    for s in S:
        v=[]; rv=[]; cv=[]; lv=[]
        for r in S[s]:
            pick=lambda a: a[np.unique(np.linspace(0,len(a)-1,k).round().astype(int))]   # per species
            cv.append(pick(r['cpx']).mean()); lv.append(pick(r['lig']).mean())
            rv.append(pick(r['rec']).mean()); v.append(-(cv[-1]-lv[-1]))
        x[s]=np.array(v); rec[s]=np.array(rv); C[s]=np.array(cv); L[s]=np.array(lv)
    pool=np.concatenate([rec['lats1'],rec['ndr2']])
    cells=[]
    for s in ('lats1','ndr2'):
        m=SCALE*(x[s].mean()+pool.mean()); e=SCALE*math.hypot(sem(C[s]),sem(L[s]),sem(pool))
        cells.append(f"\\saltvalue{{${m:.2f} \\pm {e:.2f}$}}")
    dm=SCALE*(x['lats1'].mean()-x['ndr2'].mean()); de=SCALE*math.hypot(sem(C['lats1']),sem(L['lats1']),sem(C['ndr2']),sem(L['ndr2']))
    cells.append(f"\\saltvalue{{${dm:+.2f} \\pm {de:.2f}$}}")
    print(f"{k:<2d} & " + " & ".join(cells) + r" \\")
