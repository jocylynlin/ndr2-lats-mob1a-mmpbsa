#!/usr/bin/env python3
"""Physical collective variables for the FES figures and SI Table S6.

Bound    campaign_100x40ns/analysis/physcv/{cond}/bound  -- unchanged. The complex
         topologies and trajectories are byte-identical to the originals (md5 verified),
         because the complex boxes were kept as built.
Unbound  reviewer_comments/analysis_physcv_matched/data/{matched cond}/unbound for the two
         added-salt conditions (recomputed on the matched-salt kinase trajectories), and the
         original directory for counterion-only, which has no added NaCl and was not resimulated.

CV definitions are byte-identical to physcv_one.sh: interhelical angle between the
18->31 and 38->57 CA axes, and radius of gyration of residues 1-59 over HEAVY ATOMS.
"""
import glob, os
import numpy as np, pandas as pd

OLD='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/physcv'
NEW='/home3/jocelyn/projects/hippo/reviewer_comments/analysis_physcv_matched/data'
MATCHED={'opc-salt':'opc-salt-matched','tip3p':'tip3p-matched'}
CONDS=['opc-salt','tip3p','counterion']; SYS=['lats1','ndr2']

def _dir(cond,state):
    if state=='unbound' and cond in MATCHED:
        return f'{NEW}/{MATCHED[cond]}/unbound'
    return f'{OLD}/{cond}/{state}'

def load(cond,sysn,state):
    fs=sorted(glob.glob(f'{_dir(cond,state)}/{sysn}_[0-9][0-9][0-9].csv'))
    if not fs: raise FileNotFoundError(_dir(cond,state))
    d=[pd.read_csv(f) for f in fs]
    n=min(len(x) for x in d)
    return (np.array([x['angle'].values[:n] for x in d]),
            np.array([x['rg'].values[:n] for x in d]))

def all_data():
    return {(c,s,st):load(c,s,st) for c in CONDS for s in SYS for st in ('bound','unbound')}

if __name__=='__main__':
    D=all_data()
    print(f"{'cond':11s} {'sys':6s} {'state':8s} {'replicas':>9s} {'frames':>7s} "
          f"{'angle mean':>11s} {'Rg mean':>8s}   source")
    for c in CONDS:
        for s in SYS:
            for st in ('bound','unbound'):
                a,g=D[(c,s,st)]
                src='matched' if (st=='unbound' and c in MATCHED) else 'original'
                print(f"{c:11s} {s:6s} {st:8s} {a.shape[0]:9d} {a.shape[1]:7d} "
                      f"{a.mean():11.2f} {g.mean():8.3f}   {src}")
