import pandas as pd, numpy as np
PATH='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
ROWS=[('Reference (OPC MT)','opc-salt','mt'),('Single-traj (OPC ST)','opc-salt','st'),
      ('TIP3P (MT)','tip3p','mt'),('Counterion only (MT)','counterion','mt')]
WIN=[('0-5',0,5),('5-10',5,10),('10-20',10,20),('20-40',20,40),('0-40 (reported)',0,40)]
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))
print('ddH = LATS1 minus NDR2, by analysis window. Pooled MOB1A reference cancels in the')
print('difference, so SE = hypot of the two per-replica SEMs of (cpx - lig).\n')
print(f'{"condition":22s} '+'  '.join(f'{w[0]:>16s}' for w in WIN))
for lab,c,p in ROWS:
    df=pd.read_csv(f'{PATH}mmpbsa_perframe_{c}_{p}.csv')
    col='TOTAL'
    cells=[]
    for wn,a,b in WIN:
        d=df[(df.time_ns>=a)&(df.time_ns<b)]
        vals={}
        for s in ['lats1','ndr2']:
            e=d[d.system==s]
            if p=='st': vals[s]=e.groupby('replica')[col].mean().values
            else: vals[s]=(e.groupby('replica')['cpx_'+col].mean()
                           -e.groupby('replica')['lig_'+col].mean()).values
        v=vals['lats1'].mean()-vals['ndr2'].mean()
        e_=np.hypot(sem(vals['lats1']),sem(vals['ndr2']))
        cells.append(f'{v:+6.2f}+/-{e_:4.2f}')
    print(f'{lab:22s} '+'  '.join(f'{x:>16s}' for x in cells))
print()
print('Significance of the reported (0-40) value vs the late-window (20-40) value:')
