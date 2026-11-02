import pandas as pd, numpy as np
A='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
OUT='/home3/jocelyn/projects/hippo/reviewer_comments/revision_analysis/results/campaign_n100/'
COMP=[('BOND','BOND'),('ANGLE','ANGLE'),('DIHED','DIHED'),('1-4 VDW','1-4 VDW'),('1-4 EEL','1-4 EEL')]
BLOCKS=[('opc-salt','Reference'),('tip3p','TIP3P'),('counterion','Counterion-only')]
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))
L=[]
print('Covalent decomposition, n=100, pooled MOB1A reference, replica-level SE')
for cond,lab in BLOCKS:
    mt=pd.read_csv(f'{A}mmpbsa_perframe_{cond}_mt.csv')
    M={}
    for s in ['lats1','ndr2']:
        d=mt[mt.system==s]
        M[s]={c:(d.groupby('replica')['cpx_'+c].mean().values,
                 d.groupby('replica')['rec_'+c].mean().values,
                 d.groupby('replica')['lig_'+c].mean().values) for c,_ in COMP}
    print(f'\n--- {lab} ---')
    L.append(f'\\multicolumn{{4}}{{l}}{{\\textit{{{lab}}}}} \\\\')
    sums={'lats1':[],'ndr2':[]}
    for c,disp in COMP:
        pool=np.concatenate([M['lats1'][c][1],M['ndr2'][c][1]])
        vals={}
        for s in ['lats1','ndr2']:
            x=M[s][c][0]-M[s][c][2]
            vals[s]=(x.mean()-pool.mean(), np.hypot(sem(x),sem(pool)), x)
            sums[s].append(x)
        d=vals['lats1'][2].mean()-vals['ndr2'][2].mean()
        e=np.hypot(sem(vals['lats1'][2]),sem(vals['ndr2'][2]))
        print(f'  {disp:9s} LATS1 {vals["lats1"][0]:+6.2f}+/-{vals["lats1"][1]:4.2f}   '
              f'NDR2 {vals["ndr2"][0]:+6.2f}+/-{vals["ndr2"][1]:4.2f}   '
              f'L-N {d:+6.2f}+/-{e:4.2f} ({abs(d)/e:4.1f}s)')
        f=lambda v,s_: f'${v:.2f} \\pm {s_:.2f}$'
        L.append(f'{disp.replace("1-4","1-4"):9s} & {f(*vals["lats1"][:2])} & {f(*vals["ndr2"][:2])} & {f(d,e)} \\\\')
    tot={}
    for s in ['lats1','ndr2']:
        tot[s]=np.sum(sums[s],axis=0)
    poolT=np.concatenate([np.sum([M['lats1'][c][1] for c,_ in COMP],axis=0),
                          np.sum([M['ndr2'][c][1] for c,_ in COMP],axis=0)])
    dt=tot['lats1'].mean()-tot['ndr2'].mean(); et=np.hypot(sem(tot['lats1']),sem(tot['ndr2']))
    print(f'  {"Sum":9s} LATS1 {tot["lats1"].mean()-poolT.mean():+6.2f}+/-'
          f'{np.hypot(sem(tot["lats1"]),sem(poolT)):4.2f}   '
          f'NDR2 {tot["ndr2"].mean()-poolT.mean():+6.2f}+/-{np.hypot(sem(tot["ndr2"]),sem(poolT)):4.2f}   '
          f'L-N {dt:+6.2f}+/-{et:4.2f} ({abs(dt)/et:4.1f}s)')
    L.append('\\midrule')
    L.append(f'Sum       & ${tot["lats1"].mean()-poolT.mean():.2f} \\pm '
             f'{np.hypot(sem(tot["lats1"]),sem(poolT)):.2f}$ & '
             f'${tot["ndr2"].mean()-poolT.mean():.2f} \\pm '
             f'{np.hypot(sem(tot["ndr2"]),sem(poolT)):.2f}$ & ${dt:.2f} \\pm {et:.2f}$ \\\\')
    if cond!='counterion': L.append('\\midrule')
open(OUT+'tableS_covalent_decomp_n100.tex','w').write('\n'.join(L)+'\n')
print('\nwrote tableS_covalent_decomp_n100.tex (body rows only, drop into the existing table)')
