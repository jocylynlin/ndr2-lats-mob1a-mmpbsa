#!/usr/bin/env python3
"""Frame-resolved MM-PBSA cache for the convergence figures (SI S1, S2).

The replica-mean loader in analyze_watermodel_origin.py collapses frames; the
time-series figures need the 400 individual snapshots, so this reads the same
per_frame.csv files and keeps the frame axis.

Output: tools/figures/frames_<condition>_<mode>.npz with arrays (n_replicas, 400)
  {blk}_{sys}_{quantity}   blk in cpx/rec/lig, sys in lats1/ndr2,
                           quantity in TOTAL / COV
COV = BOND+ANGLE+DIHED+UB+IMP+CMAP+1-4 VDW+1-4 EEL, the manuscript's covalent group.
"""
import csv, glob, os, sys
from collections import defaultdict
import numpy as np

ROOT='/home3/jocelyn/projects/hippo/campaign_saltmatched/mmpbsa'
OUTD=os.path.dirname(os.path.abspath(__file__))
COV=['BOND','ANGLE','DIHED','UB','IMP','CMAP','1-4 VDW','1-4 EEL']
BLOCKS=['complex','receptor','ligand','delta']
SHORT={'complex':'cpx','receptor':'rec','ligand':'lig'}

def frames(path):
    """-> {block: {col: np.array(nframes)}}"""
    out=defaultdict(lambda: defaultdict(list)); cur=None; idx=-1; hdr=None
    for row in csv.reader(open(path)):
        if not row: continue
        if row[0].startswith('Frame #'):
            hdr=row[1:]; idx+=1; cur=BLOCKS[idx] if idx<len(BLOCKS) else None; continue
        if cur is None or len(row)<10: continue
        try: vals=[float(x) for x in row[1:]]
        except (ValueError,IndexError): continue
        for h,v in zip(hdr,vals): out[cur][h].append(v)
    return {b:{k:np.array(v) for k,v in d.items()} for b,d in out.items()}

def build(cond, mode):
    store={}
    for sysn in ('lats1','ndr2'):
        dirs=sorted(glob.glob(f'{ROOT}/{cond}/{mode}/{sysn}_[0-9][0-9][0-9]'))
        acc=defaultdict(list)
        for d in dirs:
            f=os.path.join(d,'per_frame.csv')
            if not os.path.exists(f): continue
            B=frames(f)
            if not all(b in B for b in ('complex','receptor','ligand')): continue
            for b in ('complex','receptor','ligand'):
                acc[(b,'TOTAL')].append(B[b].get('TOTAL', B[b].get('DELTA TOTAL')))
                acc[(b,'COV')].append(sum(B[b][c] for c in COV if c in B[b]))
        for (b,q),v in acc.items():
            n=min(len(x) for x in v)
            store[f'{SHORT[b]}_{sysn}_{q}']=np.array([x[:n] for x in v])
        print(f"  {cond}/{mode}/{sysn}: {len(dirs)} replicas, "
              f"{store[f'cpx_{sysn}_TOTAL'].shape[1]} frames")
    p=os.path.join(OUTD,f'frames_{cond}_{mode}.npz')
    np.savez_compressed(p, **store); print(f"  wrote {p}")

if __name__=='__main__':
    for cond in ('opc-salt-matched','tip3p-matched'):
        for mode in ('mt','st'):
            build(cond,mode)
