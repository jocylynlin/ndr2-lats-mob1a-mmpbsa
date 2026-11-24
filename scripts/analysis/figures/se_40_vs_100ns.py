#!/usr/bin/env python3
"""Section 2.1: does lengthening a trajectory from 40 to 100 ns reduce the replica-level SE?

Single-trajectory dH_bind for replicas 1-20 of each complex under the reference condition, the
replicas whose production was continued to 100 ns. Complex trajectories are unaffected by the
matched-salt preparation, so the campaign_100x40ns per-frame tables are current.
"""
import numpy as np, pandas as pd
A='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
a=pd.read_csv(A+'mmpbsa_perframe_opc-salt_st.csv',usecols=['system','replica','TOTAL'])
b=pd.read_csv(A+'mmpbsa_perframe_opc-salt_st_late.csv',usecols=['system','replica','TOTAL'])
se=lambda v: v.std(ddof=1)/np.sqrt(len(v))
for s in ('lats1','ndr2'):
    reps=sorted(b[b.system==s].replica.unique())
    e=a[(a.system==s)&a.replica.isin(reps)].groupby('replica').TOTAL; l=b[b.system==s].groupby('replica').TOTAL
    m40=e.mean(); m100=(e.sum()+l.sum())/(e.count()+l.count())
    print(f"{s}: replicas {reps[0]}-{reps[-1]} (n={len(reps)})  SE 40 ns {se(m40):.3f}  100 ns {se(m100):.3f}  reduction {100*(1-se(m100)/se(m40)):.0f}%")
