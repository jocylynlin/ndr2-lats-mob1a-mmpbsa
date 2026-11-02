import numpy as np, glob, os, re
from scipy import stats
base='/home3/jocelyn/projects/hippo/campaign_100x40ns/nmode/opc-salt'
def parse(fn):
    ent=open(fn).read().split('POISSON BOLTZMANN:')[0]
    out={}
    for key,lab in [('Complex','cpx'),('Receptor','rec'),('Ligand','lig')]:
        m=re.search(key+r' Energy Terms\nFrame #,Translational,Rotational,Vibrational,Total\n((?:[-\d].*\n)+)',ent)
        out[lab]=np.array([[float(x) for x in l.split(',')] for l in m.group(1).strip().split('\n')])
    return out
def sem(x): x=np.asarray(x,float); return x.std(ddof=1)/np.sqrt(len(x))
S={}
for s in ['lats1','ndr2']:
    c,r,l=[],[],[]
    for d in sorted(glob.glob(f'{base}/{s}/traj_*')):
        o=parse(f'{d}/per_frame.csv'); n=min(len(o['cpx']),len(o['rec']),len(o['lig']))
        c.append(o['cpx'][:n,4].mean()); r.append(o['rec'][:n,4].mean()); l.append(o['lig'][:n,4].mean())
    S[s]=dict(cpx=np.array(c),rec=np.array(r),lig=np.array(l))
rL,rN=S['lats1']['rec'],S['ndr2']['rec']
gap=rL.mean()-rN.mean()
print('TS of unbound MOB1A (nmode, OPC, kcal/mol at 298.15 K):')
print(f'  LATS1-side {rL.mean():9.2f} ± {sem(rL):.2f} (SD {rL.std(ddof=1):.2f})   NDR2-side {rN.mean():9.2f} ± {sem(rN):.2f} (SD {rN.std(ddof=1):.2f})')
print(f'  gap {gap:+.2f} ± {np.hypot(sem(rL),sem(rN)):.2f}  (Welch p = {stats.ttest_ind(rL,rN,equal_var=False)[1]:.2g})')
pool=np.concatenate([rL,rN])
print(f'  pooled {pool.mean():9.2f} ± {sem(pool):.2f} (n=200)\n')
print('-T*dS (kcal/mol, 298.15 K):')
rows={}
for s in ['lats1','ndr2']:
    own=-(S[s]['cpx']-S[s]['rec']-S[s]['lig'])
    x=-(S[s]['cpx']-S[s]['lig']); pooledv=x+pool.mean()
    rows[s]=(own,x)
    print(f'  {s:6s} own reference {own.mean():7.2f} ± {sem(own):4.2f}    pooled reference {pooledv.mean():7.2f} ± {np.hypot(sem(x),sem(pool)):4.2f}')
dO=rows['lats1'][0].mean()-rows['ndr2'][0].mean(); sO=np.hypot(sem(rows['lats1'][0]),sem(rows['ndr2'][0]))
dP=rows['lats1'][1].mean()-rows['ndr2'][1].mean(); sP=np.hypot(sem(rows['lats1'][1]),sem(rows['ndr2'][1]))
print(f'\n  d(-TdS) LATS1-NDR2: own reference {dO:+.2f} ± {sO:.2f}   pooled reference {dP:+.2f} ± {sP:.2f}')
print('\n--- dG = dH + (-TdS), OPC MT, n=100, entropy at 298.15 K ---')
import pandas as pd
A='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
df=pd.read_csv(f'{A}mmpbsa_perframe_opc-salt_mt.csv',usecols=['system','replica','TOTAL','cpx_TOTAL','rec_TOTAL','lig_TOTAL'])
H={}
for s in ['lats1','ndr2']:
    g=df[df.system==s].groupby('replica')[['TOTAL','cpx_TOTAL','rec_TOTAL','lig_TOTAL']].mean()
    H[s]=g
mob_pool=np.concatenate([H['lats1']['rec_TOTAL'].values,H['ndr2']['rec_TOTAL'].values])
for ref in ['own','pooled']:
    vals={}
    for s in ['lats1','ndr2']:
        if ref=='own':
            dh=H[s]['TOTAL'].values; sdh=sem(dh)
            tds=rows[s][0]; stds=sem(tds)
        else:
            dh=(H[s]['cpx_TOTAL']-H[s]['lig_TOTAL']).values-mob_pool.mean(); sdh=np.hypot(sem(H[s]['cpx_TOTAL']-H[s]['lig_TOTAL']),sem(mob_pool))
            tds=rows[s][1]+pool.mean(); stds=np.hypot(sem(rows[s][1]),sem(pool))
        g=dh.mean()+tds.mean(); sg=np.hypot(sdh,stds)
        vals[s]=(g,sg)
        print(f'  {ref:6s} {s:6s} dG {g:+7.2f} ± {sg:4.2f}')
    d=vals['lats1'][0]-vals['ndr2'][0]; sd=np.hypot(vals['lats1'][1],vals['ndr2'][1])
    print(f'  {ref:6s} ddG (LATS1-NDR2) {d:+7.2f} ± {sd:4.2f}\n')
