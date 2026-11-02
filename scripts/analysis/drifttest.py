import pandas as pd, numpy as np
PATH='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/'
COV=['BOND','ANGLE','DIHED','UB','IMP','CMAP','1-4 VDW','1-4 EEL']
ROWS=[('Reference','opc-salt','mt'),('Single-traj','opc-salt','st'),
      ('TIP3P','tip3p','mt'),('Counterion only','counterion','mt')]
SYS=['lats1','ndr2']
def series(cond,prot,which):
    df=pd.read_csv(f'{PATH}mmpbsa_perframe_{cond}_{prot}.csv')
    for pre in ['','cpx_','rec_','lig_']:
        df[pre+'COV']=sum(df[pre+c] for c in COV)
    col={'dH':'TOTAL','cov':'COV'}[which]
    out={}
    if prot=='st':
        for s in SYS:
            P=df[df.system==s].pivot_table(index='replica',columns='time_ns',values=col)
            out[s]=(P.columns.values,P.values.mean(0))
    else:
        R=np.vstack([df[df.system==s].pivot_table(index='replica',columns='time_ns',
                     values='rec_'+col).values for s in SYS])
        rm=R.mean(0)
        for s in SYS:
            d=df[df.system==s]
            Pc=d.pivot_table(index='replica',columns='time_ns',values='cpx_'+col)
            Pl=d.pivot_table(index='replica',columns='time_ns',values='lig_'+col)
            out[s]=(Pc.columns.values,(Pc.values-Pl.values).mean(0)-rm)
    return out
def trend(t,m,t0):
    k=t>=t0; X=np.vstack([np.ones(k.sum()),t[k]]).T; y=m[k]
    b=np.linalg.lstsq(X,y,rcond=None)[0]; r=y-X@b
    se=np.sqrt((r@r/(k.sum()-2))/(((t[k]-t[k].mean())**2).sum()))
    return b[1],se
for which,lab2 in [('dH','dH_bind'),('cov','dE_cov')]:
    print(f'\n=== {lab2}: early-vs-late on the per-time-point mean (kcal/mol) ===')
    print(f'{"condition":17s} {"sys":6s} {"0-5ns":>8s} {"5-10":>8s} {"10-20":>8s} {"20-40":>8s} '
          f'{"late-early":>11s} {"slope 10-40 ns":>22s}')
    for lab,c,p in ROWS:
        d=series(c,p,which)
        for s in SYS:
            t,m=d[s]; w=lambda a,b: m[(t>=a)&(t<b)].mean()
            sl,sse=trend(t,m,10.0)
            flag=' <== DRIFT' if abs(sl)>2*sse and abs(w(20,40)-w(0,5))>2.0 else ''
            print(f'{lab:17s} {s:6s} {w(0,5):8.2f} {w(5,10):8.2f} {w(10,20):8.2f} {w(20,40):8.2f} '
                  f'{w(20,40)-w(0,5):11.2f} {sl:+10.4f} +/- {sse:.4f}{flag}')
