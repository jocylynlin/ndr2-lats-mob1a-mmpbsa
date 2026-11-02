import pandas as pd, numpy as np, glob, os, re
A='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
COV=['BOND','ANGLE','DIHED','UB','IMP','CMAP','1-4 VDW','1-4 EEL']
def sem(x): x=np.asarray(x,float); return x.std(ddof=1)/np.sqrt(len(x))
def grp(df):
    return pd.DataFrame({'Covalent':df[COV].sum(axis=1),'Electrostatic':df['EEL']+df['EPB'],
        'Dispersion':df['VDWAALS']+df['EDISPER'],'Hydrophobic':df['ENPOLAR'],'dH':df['TOTAL'],
        'system':df['system'],'replica':df['replica']})
cols=['Covalent','Electrostatic','Dispersion','Hydrophobic','dH']
V={}
for prot in ['st','mt']:
    g=grp(pd.read_csv(f'{A}mmpbsa_perframe_opc-salt_{prot}.csv')).groupby(['system','replica']).mean()
    for s in ['lats1','ndr2']:
        for c in cols: V[(prot,s,c)]=g.loc[s,c].values

# --- entropy: nmode at 298.15 K, evaluated at the simulation temperature 310.15 K (manuscript convention)
base='/home3/jocelyn/projects/hippo/campaign_100x40ns/nmode/opc-salt'
SC=310.15/298.15
def parse(fn):
    ent=open(fn).read().split('POISSON BOLTZMANN:')[0]; o={}
    for k,l in [('Complex','cpx'),('Receptor','rec'),('Ligand','lig')]:
        m=re.search(k+r' Energy Terms\nFrame #,Translational,Rotational,Vibrational,Total\n((?:[-\d].*\n)+)',ent)
        o[l]=np.array([[float(x) for x in r.split(',')] for r in m.group(1).strip().split('\n')])
    return o
TDS={}
for s in ['lats1','ndr2']:
    a=[]
    for d in sorted(glob.glob(f'{base}/{s}/traj_*')):
        o=parse(f'{d}/per_frame.csv'); n=min(len(o['cpx']),len(o['rec']),len(o['lig']))
        a.append(-(o['cpx'][:n,4]-o['rec'][:n,4]-o['lig'][:n,4]).mean())
    TDS[s]=np.array(a)*SC

print('=== Table 1 replacement values, OPC, n=100 (replica-level SE) ===\n')
for prot,lab in [('st','Single-trajectory'),('mt','Multi-trajectory')]:
    print(f'-- {lab} --')
    for s in ['lats1','ndr2']:
        r=[]
        for c in cols:
            v=V[(prot,s,c)]
            r.append(f'{c}: {v.mean():8.2f} ± {sem(v):4.2f}')
        line='   '.join(r)
        if prot=='mt':
            dG=V[(prot,s,'dH')].mean()+TDS[s].mean()
            sdG=np.hypot(sem(V[(prot,s,'dH')]),sem(TDS[s]))
            line+=f'   dG: {dG:+7.2f} ± {sdG:4.2f}'
        print(f'  {s.upper():6s} {line}')
    r=[]
    for c in cols:
        L,N=V[(prot,'lats1',c)],V[(prot,'ndr2',c)]
        d=L.mean()-N.mean(); s_=np.hypot(sem(L),sem(N))
        r.append(f'{c}: {d:+8.2f} ± {s_:4.2f}')
    line='   '.join(r)
    if prot=='mt':
        gL=V[(prot,'lats1','dH')].mean()+TDS['lats1'].mean(); gN=V[(prot,'ndr2','dH')].mean()+TDS['ndr2'].mean()
        sL=np.hypot(sem(V[(prot,'lats1','dH')]),sem(TDS['lats1'])); sN=np.hypot(sem(V[(prot,'ndr2','dH')]),sem(TDS['ndr2']))
        line+=f'   dG: {gL-gN:+7.2f} ± {np.hypot(sL,sN):4.2f}'
    print(f'  DIFF   {line}')
    print()
print(f'-TdS at 310.15 K:  LATS1 {TDS["lats1"].mean():.2f} ± {sem(TDS["lats1"]):.2f}   '
      f'NDR2 {TDS["ndr2"].mean():.2f} ± {sem(TDS["ndr2"]):.2f}   '
      f'diff {TDS["lats1"].mean()-TDS["ndr2"].mean():+.2f} ± {np.hypot(sem(TDS["lats1"]),sem(TDS["ndr2"])):.2f}')
print(f'  (at 298.15 K unscaled: {TDS["lats1"].mean()/SC:.2f} / {TDS["ndr2"].mean()/SC:.2f}, diff {(TDS["lats1"].mean()-TDS["ndr2"].mean())/SC:+.2f})')
# shared-reference MT block
print('\n=== optional shared-MOB1A-reference MT block ===')
df=pd.read_csv(f'{A}mmpbsa_perframe_opc-salt_mt.csv')
S={}
for s in ['lats1','ndr2']:
    d=df[df.system==s]
    for p in ['cpx_','rec_','lig_']:
        gg=pd.DataFrame({'Covalent':d[[p+c for c in COV]].sum(axis=1),'Electrostatic':d[p+'EEL']+d[p+'EPB'],
            'Dispersion':d[p+'VDWAALS']+d[p+'EDISPER'],'Hydrophobic':d[p+'ENPOLAR'],'dH':d[p+'TOTAL'],'replica':d['replica']})
        S[(s,p)]=gg.groupby('replica').mean()
for c in cols:
    pool=np.concatenate([S[('lats1','rec_')][c].values,S[('ndr2','rec_')][c].values]); pm=pool.mean()
    out=[]
    for s in ['lats1','ndr2']:
        x=(S[(s,'cpx_')][c]-S[(s,'lig_')][c]).values
        out.append((x.mean()-pm, np.hypot(sem(x),sem(pool)), x))
    d=out[0][2].mean()-out[1][2].mean(); sd=np.hypot(sem(out[0][2]),sem(out[1][2]))
    print(f'  {c:14s} LATS1 {out[0][0]:8.2f} ± {out[0][1]:4.2f}   NDR2 {out[1][0]:8.2f} ± {out[1][1]:4.2f}   diff {d:+7.2f} ± {sd:4.2f}')

print('\n=== shared-reference dG (310.15 K) ===')
ENT={}
for s in ['lats1','ndr2']:
    c,r,l=[],[],[]
    for d in sorted(glob.glob(f'{base}/{s}/traj_*')):
        o=parse(f'{d}/per_frame.csv'); n=min(len(o['cpx']),len(o['rec']),len(o['lig']))
        c.append(o['cpx'][:n,4].mean()); r.append(o['rec'][:n,4].mean()); l.append(o['lig'][:n,4].mean())
    ENT[s]=(np.array(c),np.array(r),np.array(l))
rpool=np.concatenate([ENT['lats1'][1],ENT['ndr2'][1]])
G={}
for s in ['lats1','ndr2']:
    tds=(-(ENT[s][0]-ENT[s][2])+rpool.mean())*SC
    stds=np.hypot(sem(-(ENT[s][0]-ENT[s][2])*SC),sem(rpool*SC))
    hpool=np.concatenate([S[('lats1','rec_')]['dH'].values,S[('ndr2','rec_')]['dH'].values])
    x=(S[(s,'cpx_')]['dH']-S[(s,'lig_')]['dH']).values
    dh=x.mean()-hpool.mean(); sdh=np.hypot(sem(x),sem(hpool))
    G[s]=(dh+tds.mean(), np.hypot(sdh,stds), x, tds)
    print(f'  {s:6s} -TdS {tds.mean():6.2f} ± {stds:4.2f}   dH {dh:+7.2f} ± {sdh:4.2f}   dG {dh+tds.mean():+7.2f} ± {np.hypot(sdh,stds):4.2f}')
dg=G['lats1'][0]-G['ndr2'][0]
sdg=np.hypot(np.hypot(sem(G['lats1'][2]),sem(G['lats1'][3])),np.hypot(sem(G['ndr2'][2]),sem(G['ndr2'][3])))
tl=(-(ENT['lats1'][0]-ENT['lats1'][2]))*SC; tn=(-(ENT['ndr2'][0]-ENT['ndr2'][2]))*SC
print(f'  diff -TdS {tl.mean()-tn.mean():+.2f} ± {np.hypot(sem(tl),sem(tn)):.2f}   dG diff {dg:+.2f} ± {sdg:.2f}')
