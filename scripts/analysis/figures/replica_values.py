#!/usr/bin/env python3
"""Per-replica MM-PBSA values for the SI distribution figures (S3, S4).

One accessor over two data sources:
  opc-salt-matched / tip3p-matched -> campaign_saltmatched per_frame.csv (final dataset)
  counterion                       -> campaign_100x40ns tidied CSV, unchanged, because
                                      that condition has no added NaCl and was not resimulated.

MT uses the unbound MOB1A reference POOLED over both systems within each condition, as
Tables 1-3 do, so panel means reproduce the table entries exactly.
"""
import sys, math
import numpy as np, pandas as pd
sys.path.insert(0,'/home3/jocelyn/projects/hippo/reviewer_comments/analysis_watermodel_origin/scripts')
from analyze_watermodel_origin import load

OLD='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis'
COV5=['BOND','ANGLE','DIHED','1-4 VDW','1-4 EEL']
SYS=['lats1','ndr2']
MATCHED={'opc-salt':'opc-salt-matched','tip3p':'tip3p-matched'}
se=lambda v: float(np.std(v,ddof=1)/math.sqrt(len(v)))
_cache={}

def per_replica(cond, proto, term='TOTAL'):
    """term: 'TOTAL', 'COV5', or a single raw column such as 'DIHED'.
    -> {system: np.array(n_replicas)}"""
    key=(cond,proto,term)
    if key in _cache: return _cache[key]
    if cond in MATCHED:
        per=_cache.setdefault(('_load',cond,proto), load(MATCHED[cond],proto))
        reps={s:sorted(r for (ss,r) in per if ss==s) for s in SYS}
        def grab(s,blk):
            cols = COV5 if term=='COV5' else [term]
            return np.array([sum(per[(s,r)][blk][t] for t in cols) for r in reps[s]])
        out={}
        if proto=='st':
            for s in SYS: out[s]=grab(s,'complex')-grab(s,'receptor')-grab(s,'ligand')
        else:
            pool=np.concatenate([grab(s,'receptor') for s in SYS]).mean()
            for s in SYS: out[s]=grab(s,'complex')-pool-grab(s,'ligand')
    else:
        d=_cache.setdefault(('_csv',cond,proto), pd.read_csv(f'{OLD}/mmpbsa_perframe_{cond}_{proto}.csv'))
        def grab(s,pre):
            g=d[d.system==s]
            cols = COV5 if term=='COV5' else [term]
            return sum(g.groupby('replica')[pre+t].mean().values for t in cols)
        out={}
        if proto=='st':
            for s in SYS: out[s]=grab(s,'cpx_')-grab(s,'rec_')-grab(s,'lig_')
        else:
            pool=np.concatenate([grab(s,'rec_') for s in SYS]).mean()
            for s in SYS: out[s]=grab(s,'cpx_')-pool-grab(s,'lig_')
    _cache[key]=out
    return out

def mt_parts(cond, term='TOTAL'):
    """Independent replica means for the multi-trajectory protocol.
    -> ({system: complex}, {system: kinase segment}, pooled MOB1A array of 200).
    Complex, kinase-segment and MOB1A simulations are independent, so uncertainties are
    combined in quadrature from these arrays; replicas are never paired by index."""
    cols = COV5 if term=='COV5' else [term]
    if cond in MATCHED:
        per=_cache.setdefault(('_load',cond,'mt'), load(MATCHED[cond],'mt'))
        reps={s:sorted(r for (ss,r) in per if ss==s) for s in SYS}
        g=lambda s,blk: np.array([sum(per[(s,r)][blk][t] for t in cols) for r in reps[s]])
        c={s:g(s,'complex') for s in SYS}; l={s:g(s,'ligand') for s in SYS}; r={s:g(s,'receptor') for s in SYS}
    else:
        d=_cache.setdefault(('_csv',cond,'mt'), pd.read_csv(f'{OLD}/mmpbsa_perframe_{cond}_mt.csv'))
        g=lambda s,pre: sum(d[d.system==s].groupby('replica')[pre+t].mean().values for t in cols)
        c={s:g(s,'cpx_') for s in SYS}; l={s:g(s,'lig_') for s in SYS}; r={s:g(s,'rec_') for s in SYS}
    return c, l, np.concatenate([r[s] for s in SYS])

def quad(*arrs):
    return math.sqrt(sum(se(a)**2 for a in arrs))

def mt_summary(cond, term='TOTAL'):
    """-> {'lats1':(mean,se), 'ndr2':(mean,se), 'diff':(mean,se)}; per-system SEs include the
    pooled reference, the difference excludes it (it cancels exactly)."""
    c,l,pool=mt_parts(cond,term)
    out={s:(c[s].mean()-l[s].mean()-pool.mean(), quad(c[s],l[s],pool)) for s in SYS}
    out['diff']=((c['lats1'].mean()-l['lats1'].mean())-(c['ndr2'].mean()-l['ndr2'].mean()),
                 quad(c['lats1'],l['lats1'],c['ndr2'],l['ndr2']))
    return out

def mt_bootstrap(cond, term, nb, rng):
    """Bootstrap of the LATS1-minus-NDR2 difference with complex and kinase-segment replicas
    resampled independently (the pooled reference cancels)."""
    c,l,_=mt_parts(cond,term)
    m=lambda a: a[rng.integers(0,len(a),(nb,len(a)))].mean(1)
    return (m(c['lats1'])-m(l['lats1']))-(m(c['ndr2'])-m(l['ndr2']))

if __name__=='__main__':
    print("per-replica means (should reproduce Tables 1-3 and S1)")
    for cond in ('opc-salt','tip3p','counterion'):
        for proto in ('st','mt'):
            for term in ('TOTAL','COV5'):
                try: v=per_replica(cond,proto,term)
                except Exception as e: print(f"  {cond}/{proto}/{term}: {e}"); continue
                if proto=='mt':
                    q=mt_summary(cond,term)
                    print(f"  {cond:11s} mt {term:6s}  LATS1 {q['lats1'][0]:+8.2f}+-{q['lats1'][1]:.2f}  "
                          f"NDR2 {q['ndr2'][0]:+8.2f}+-{q['ndr2'][1]:.2f}  diff {q['diff'][0]:+6.2f}+-{q['diff'][1]:.2f}")
                    continue
                print(f"  {cond:11s} {proto} {term:6s}  "
                      f"LATS1 {v['lats1'].mean():+8.2f}+-{se(v['lats1']):.2f}  "
                      f"NDR2 {v['ndr2'].mean():+8.2f}+-{se(v['ndr2']):.2f}  "
                      f"diff {v['lats1'].mean()-v['ndr2'].mean():+6.2f}+-"
                      f"{math.hypot(se(v['lats1']),se(v['ndr2'])):.2f}")
