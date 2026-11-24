#!/usr/bin/env python3
"""SI Table S2 -- absolute bound and unbound covalent energies, final dataset.

Bound   = complex-block covalent components (complex trajectory; unchanged by the
          matched-salt preparation, since the complex boxes were kept as built).
Unbound = pooled unbound MOB1A reference + unbound kinase, both resimulated.
Delta   = bound - unbound, which is identically the per-system column of Table S1,
          and its SE includes the pooled reference, as Table S1's does.
Counterion-only is NOT recomputed and is left exactly as published: that condition
contains no added NaCl and was not resimulated.
"""
import sys, math
import numpy as np
sys.path.insert(0,'/home3/jocelyn/projects/hippo/reviewer_comments/analysis_watermodel_origin/scripts')
from analyze_watermodel_origin import load
COV=['BOND','ANGLE','DIHED','1-4 VDW','1-4 EEL']
LBL={'BOND':'BOND','ANGLE':'ANGLE','DIHED':'DIHED','1-4 VDW':'1--4 VDW','1-4 EEL':'1--4 EEL'}
se=lambda a: float(np.std(a,ddof=1)/math.sqrt(len(a)))

def block(cond_dir):
    per=load(cond_dir,'mt')
    reps={s:sorted(r for (ss,r) in per if ss==s) for s in ('lats1','ndr2')}
    A=lambda s,b,t: np.array([per[(s,r)][b][t] for r in reps[s]])
    rows={}
    for t in COV+['TOTAL']:
        cols=COV if t=='TOTAL' else [t]
        pool=np.concatenate([sum(A(s,'receptor',c) for c in cols) for s in ('lats1','ndr2')])
        for s in ('lats1','ndr2'):
            cpx=sum(A(s,'complex',c) for c in cols); lig=sum(A(s,'ligand',c) for c in cols)
            rows[(s,t)]=(cpx.mean(), pool.mean()+lig.mean(),
                         cpx.mean()-pool.mean()-lig.mean(),
                         math.hypot(se(cpx),se(pool),se(lig)))
    return rows

for cond,cdir in (('Reference','opc-salt-matched'),('TIP3P','tip3p-matched')):
    r=block(cdir)
    print(f"\n% ---- {cond} ----")
    for t in COV:
        c=[]
        for s in ('lats1','ndr2'):
            b,u,d,e=r[(s,t)]
            c.append(f"\\saltvalue{{{b:.1f}}} & \\saltvalue{{{u:.1f}}} & \\saltvalue{{${d:+.2f} \\pm {e:.2f}$}}")
        print(f"{LBL[t]:9s} & " + " & ".join(c) + r" \\")
    c=[]
    for s in ('lats1','ndr2'):
        b,u,d,e=r[(s,'TOTAL')]
        c.append(f"\\saltvalue{{{b:.1f}}} & \\saltvalue{{{u:.1f}}} & \\saltvalue{{${d:+.2f} \\pm {e:.2f}$}}")
    print("Total     & " + " & ".join(c) + r" \\")
