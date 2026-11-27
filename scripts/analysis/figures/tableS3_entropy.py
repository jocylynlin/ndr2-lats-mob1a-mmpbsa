#!/usr/bin/env python3
"""SI Table S3 -- normal-mode entropy decomposition, final dataset.

Source  campaign_saltmatched/nmode/opc-salt-matched/{lats1,ndr2}/traj_NNN/per_frame.csv
        100 replicas per system, 10 snapshots each, harmonic approximation.

Follows revision_analysis/.../scripts/ent_pool.py, the script behind the published
entropy numbers, exactly:
  * snapshot-level per_frame.csv, NOT the replica-averaged FINAL_RESULTS file;
  * all converged snapshots of each species are used; the single non-converged ligand
    snapshot (ndr2 traj_067) is the only one missing;
  * replica means are the statistical unit; SE = s_replica / sqrt(100).
Pooled unbound MOB1A reference (200 replicas): per-system SEs INCLUDE it, the
LATS1-minus-NDR2 difference EXCLUDES it because it cancels exactly.
MMPBSA.py reports T*S at 298.15 K and the manuscript reports -T*dS at 310.15 K.
"""
import glob, os, re, math, sys
import numpy as np
ROOT='/home3/jocelyn/projects/hippo/campaign_saltmatched/nmode/opc-salt-matched'
T_RUN, T_REP = 298.15, 310.15
COLS={'Translational':1,'Rotational':2,'Vibrational':3,'Total':4}
TERMS=list(COLS)
sem=lambda x: float(np.std(np.asarray(x,float),ddof=1)/math.sqrt(len(x)))

def parse(fn):
    ent=open(fn).read().split('POISSON BOLTZMANN:')[0]
    out={}
    for key,lab in [('Complex','cpx'),('Receptor','rec'),('Ligand','lig')]:
        m=re.search(key+r' Energy Terms\nFrame #,Translational,Rotational,Vibrational,Total\n((?:[-\d].*\n)+)',ent)
        if m is None: return None
        out[lab]=np.array([[float(x) for x in l.split(',')] for l in m.group(1).strip().split('\n')])
    return out

S={s:{k:[] for k in ('cpx','rec','lig')} for s in ('lats1','ndr2')}
ndrop=0
for s in ('lats1','ndr2'):
    for d in sorted(glob.glob(f'{ROOT}/{s}/traj_[0-9][0-9][0-9]')):
        o=parse(f'{d}/per_frame.csv')
        if o is None: continue
        # every converged snapshot of every species is kept; the species are independent
        # simulations, so a failed ligand snapshot does not remove complex or receptor rows
        ndrop+=sum(10-len(o[k]) for k in ('cpx','rec','lig'))
        for k in ('cpx','rec','lig'): S[s][k].append(o[k])
for s in S:
    print(f"{s}: {len(S[s]['cpx'])} replicas")
print(f"snapshots missing (non-converged): {ndrop}\n")

SCALE=T_REP/T_RUN
rows={}
for t,ci in COLS.items():
    rec={s:np.array([r[:,ci].mean() for r in S[s]['rec']]) for s in S}
    cpx={s:np.array([r[:,ci].mean() for r in S[s]['cpx']]) for s in S}
    lig={s:np.array([r[:,ci].mean() for r in S[s]['lig']]) for s in S}
    pool=np.concatenate([rec['lats1'],rec['ndr2']])
    x={s:-(cpx[s]-lig[s]) for s in S}                       # -(C - L); reference added below
    for s in S:
        rows[(s,t)]=(SCALE*(x[s].mean()+pool.mean()), SCALE*math.hypot(sem(cpx[s]),sem(lig[s]),sem(pool)))
    rows[('diff',t)]=(SCALE*(x['lats1'].mean()-x['ndr2'].mean()),
                      SCALE*math.hypot(sem(cpx['lats1']),sem(lig['lats1']),sem(cpx['ndr2']),sem(lig['ndr2'])))

LBL={'Translational':r'$-T\Delta S_\text{trans}$','Rotational':r'$-T\Delta S_\text{rot}$',
     'Vibrational':r'$-T\Delta S_\text{vib}$','Total':r'$-T\Delta S_\text{total}$'}
print(f"-T*dS rescaled from {T_RUN} K to {T_REP} K (factor {SCALE:.6f}):\n")
for t in TERMS:
    c=[]
    for key in (('lats1',t),('ndr2',t),('diff',t)):
        m,e=rows[key]
        c.append(f"$+{m:.2f} \\pm {e:.2f}$" if key[0]!='diff' else f"${m:+.2f} \\pm {e:.2f}$")
    print(f"{LBL[t]:34s} & " + " & ".join(f"\\saltvalue{{{x}}}" for x in c) + r" \\")
chk=sum(rows[('lats1',t)][0] for t in ('Translational','Rotational','Vibrational'))
print(f"closure: trans+rot+vib for LATS1 = {chk:.4f}  vs total {rows[('lats1','Total')][0]:.4f}")
