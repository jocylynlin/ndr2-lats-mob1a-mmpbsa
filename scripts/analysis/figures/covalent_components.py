#!/usr/bin/env python3
"""Per-replica covalent components from the final matched-salt dataset.

Shared by the main-text covalent decomposition figure and its verification.
Provenance: campaign_saltmatched/mmpbsa/{cond}/{mode}/{sys}_NNN/per_frame.csv
Conventions (manuscript Tables 1-3):
  MT  dE_n = C_n - K_n - Rbar,  Rbar = unbound MOB1A pooled over BOTH systems (200 replicas)
  SE  replica-level, s_replica/sqrt(N); per-system SEs INCLUDE the pooled reference,
      LATS1-minus-NDR2 differences EXCLUDE it because it cancels exactly.
"""
import math, sys
import numpy as np
sys.path.insert(0,'/home3/jocelyn/projects/hippo/reviewer_comments/analysis_watermodel_origin/scripts')
from analyze_watermodel_origin import load

COV   = ['BOND','ANGLE','DIHED','1-4 VDW','1-4 EEL']
CONDS = {'Reference':'opc-salt-matched', 'TIP3P':'tip3p-matched'}
SYS   = ['lats1','ndr2']
se    = lambda a: float(np.std(a, ddof=1)/math.sqrt(len(a)))

def components(cond_dir, mode='mt'):
    """-> {(system, term): (mean, se)}, {(term,): (diff, se)}, and the per-replica arrays."""
    per  = load(cond_dir, mode)
    reps = {s: sorted(r for (ss,r) in per if ss==s) for s in SYS}
    arr  = {}
    for s in SYS:
        for b in ('complex','ligand','receptor'):
            for t in COV:
                arr[(s,b,t)] = np.array([per[(s,r)][b][t] for r in reps[s]])
    out, diff = {}, {}
    for t in COV:
        pool = np.concatenate([arr[(s,'receptor',t)] for s in SYS])   # 200 replicas
        rbar, rse = pool.mean(), se(pool)
        for s in SYS:
            d = arr[(s,'complex',t)] - arr[(s,'ligand',t)]
            out[(s,t)] = (float(d.mean() - rbar),
                          math.hypot(se(arr[(s,'complex',t)]), rse, se(arr[(s,'ligand',t)])))
        dl = arr[('lats1','complex',t)] - arr[('lats1','ligand',t)]
        dn = arr[('ndr2','complex',t)]  - arr[('ndr2','ligand',t)]
        diff[t] = (float(dl.mean() - dn.mean()),                                 # Rbar cancels
                   math.hypot(*(se(arr[(s_,b_,t)]) for s_ in SYS for b_ in ('complex','ligand'))))
    # totals over the five components, from the per-replica sums (not from rounded parts)
    poolT = np.concatenate([sum(arr[(s,'receptor',t)] for t in COV) for s in SYS])
    for s in SYS:
        dT = sum(arr[(s,'complex',t)] - arr[(s,'ligand',t)] for t in COV)
        out[(s,'Sum')] = (float(dT.mean() - poolT.mean()),
                          math.hypot(se(sum(arr[(s,'complex',t)] for t in COV)), se(poolT),
                                     se(sum(arr[(s,'ligand',t)] for t in COV))))
    dlT = sum(arr[('lats1','complex',t)] - arr[('lats1','ligand',t)] for t in COV)
    dnT = sum(arr[('ndr2','complex',t)]  - arr[('ndr2','ligand',t)]  for t in COV)
    diff['Sum'] = (float(dlT.mean() - dnT.mean()),
                   math.hypot(*(se(sum(arr[(s_,b_,t)] for t in COV)) for s_ in SYS for b_ in ('complex','ligand'))))
    return out, diff, len(reps['lats1']), len(reps['ndr2'])

if __name__ == '__main__':
    # Table S1 as currently printed, for verification
    S1 = {
     'Reference': {('lats1','BOND'):(-0.39,0.16),('ndr2','BOND'):(-0.22,0.17),('BOND',):(-0.17,0.21),
                   ('lats1','ANGLE'):(0.04,0.36),('ndr2','ANGLE'):(0.30,0.38),('ANGLE',):(-0.27,0.49),
                   ('lats1','DIHED'):(-2.18,0.51),('ndr2','DIHED'):(-5.06,0.75),('DIHED',):(2.88,0.85),
                   ('lats1','1-4 VDW'):(-0.13,0.15),('ndr2','1-4 VDW'):(0.83,0.17),('1-4 VDW',):(-0.96,0.22),
                   ('lats1','1-4 EEL'):(-0.48,0.10),('ndr2','1-4 EEL'):(-1.48,0.11),('1-4 EEL',):(1.01,0.13),
                   ('lats1','Sum'):(-3.14,0.70),('ndr2','Sum'):(-5.63,0.90),('Sum',):(2.49,1.07)},
     'TIP3P':     {('lats1','BOND'):(-0.54,0.16),('ndr2','BOND'):(-0.45,0.17),('BOND',):(-0.08,0.21),
                   ('lats1','ANGLE'):(-1.57,0.35),('ndr2','ANGLE'):(-0.59,0.36),('ANGLE',):(-0.99,0.46),
                   ('lats1','DIHED'):(-5.20,0.46),('ndr2','DIHED'):(-4.77,0.67),('DIHED',):(-0.43,0.75),
                   ('lats1','1-4 VDW'):(-0.50,0.15),('ndr2','1-4 VDW'):(1.09,0.16),('1-4 VDW',):(-1.59,0.20),
                   ('lats1','1-4 EEL'):(-0.69,0.09),('ndr2','1-4 EEL'):(-1.41,0.11),('1-4 EEL',):(0.73,0.13),
                   ('lats1','Sum'):(-8.49,0.71),('ndr2','Sum'):(-6.13,0.82),('Sum',):(-2.36,0.99)},
    }
    bad=0
    for cond, cdir in CONDS.items():
        out, diff, nl, nn = components(cdir)
        print(f"\n=== {cond}  ({cdir}, MT, n={nl}/{nn} replicas) ===")
        print(f"{'component':10s} {'LATS1':>17s} {'NDR2':>17s} {'LATS1-NDR2':>17s}")
        for t in COV+['Sum']:
            cells=[]
            for key in [('lats1',t),('ndr2',t),(t,)]:
                m,e = out[key] if len(key)==2 else diff[t]
                ref = S1[cond].get(key)
                flag=''
                if ref:
                    if abs(round(m,2)-ref[0])>5e-3 or abs(round(e,2)-ref[1])>5e-3:
                        flag='  <-- MISMATCH vs Table S1 %+.2f+-%.2f'%ref; bad+=1
                cells.append(f"{m:+8.2f}+-{e:.2f}{flag}")
            print(f"{t:10s} " + " ".join(f"{c:>17s}" for c in cells))
        # closure: do the five components reconstruct the total?
        for s in SYS:
            parts=sum(out[(s,t)][0] for t in COV); tot=out[(s,'Sum')][0]
            print(f"  closure {s:5s}: sum of components {parts:+.6f}  vs total {tot:+.6f}  |diff| {abs(parts-tot):.2e}")
        parts=sum(diff[t][0] for t in COV); tot=diff['Sum'][0]
        print(f"  closure diff : sum of components {parts:+.6f}  vs total {tot:+.6f}  |diff| {abs(parts-tot):.2e}")
    print(f"\nmismatches vs Table S1: {bad}")
