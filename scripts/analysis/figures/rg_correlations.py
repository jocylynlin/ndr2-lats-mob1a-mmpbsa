#!/usr/bin/env python3
"""Per-replica Rg of the unbound kinase against its surface area and covalent energy (Section 2.2).

Rg       analysis_physcv_matched (heavy atoms, replica mean over 400 frames)
surface  unbound-kinase ENPOLAR / gamma, gamma = 0.0378 kcal/mol/A^2 (inp=2 cavity term), i.e. the
         solvent-accessible area the MM-PBSA step itself uses
covalent unbound-kinase covalent group from the same MM-PBSA evaluation
Replicas pooled across the two water models, 200 per kinase, as in the text.
"""
import sys, os, glob
import numpy as np, pandas as pd
sys.path.insert(0,'/home3/jocelyn/projects/hippo/reviewer_comments/analysis_watermodel_origin/scripts')
from analyze_watermodel_origin import load
P='/home3/jocelyn/projects/hippo/reviewer_comments/analysis_physcv_matched/data'
GAMMA=0.0378
for s in ('lats1','ndr2'):
    rg=[];sa=[];cv=[]
    for cond in ('opc-salt-matched','tip3p-matched'):
        per=load(cond,'mt')
        for r in sorted(rr for (ss,rr) in per if ss==s):
            f=f'{P}/{cond}/unbound/{s}_{r:03d}.csv'
            if not os.path.exists(f): continue
            rg.append(pd.read_csv(f)['rg'].mean())
            sa.append(per[(s,r)]['ligand']['Hydrophobic']/GAMMA)
            cv.append(per[(s,r)]['ligand']['Covalent'])
    rg,sa,cv=map(np.array,(rg,sa,cv))
    r1=np.corrcoef(rg,sa)[0,1]; r2=np.corrcoef(rg,cv)[0,1]
    print(f"{s}: n={len(rg)}  r(Rg, surface area) = {r1:+.2f}   r(Rg, covalent) = {r2:+.2f}  ({100*r2*r2:.0f}% of variance)")
print("text currently: +0.83 / +0.82 and -0.14 / -0.15 (about 2%)")
