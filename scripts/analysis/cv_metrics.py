import pandas as pd, numpy as np, glob, json
P='/home3/jocelyn/projects/hippo/campaign_100x40ns/analysis/physcv/'
CONDS=['opc-salt','tip3p','counterion']; SYS=['lats1','ndr2']
rng=np.random.default_rng(0)
def load(c,s,st):
    fs=sorted(glob.glob(f'{P}{c}/{st}/{s}_*.csv')); d=[pd.read_csv(f) for f in fs]
    return np.array([x['angle'].values for x in d]), np.array([x['rg'].values for x in d])
D={(c,s,st):load(c,s,st) for c in CONDS for s in SYS for st in ['bound','unbound']}
def sem(x): return np.std(x,ddof=1)/np.sqrt(len(x))

print('='*98)
print('COUNTERION vs OPC (both OPC water): is there ANY sampling difference? (unbound)')
print('='*98)
for s in SYS:
    print(f'\n--- {s} ---')
    for lab,i in [('interhelical angle',0),('Rg',1)]:
        o=D[('opc-salt',s,'unbound')][i]; c=D[('counterion',s,'unbound')][i]
        so,sc=o.std(1),c.std(1); d=so.mean()-sc.mean(); e=np.hypot(sem(so),sem(sc))
        print(f'  {lab:19s} within-replica SD: OPC {so.mean():6.2f}  counterion {sc.mean():6.2f}  '
              f'diff {d:+6.2f}+/-{e:4.2f}  -> {"DIFFERENT" if abs(d)>2*e else "NO DIFFERENCE"}')
        print(f'  {"":19s} pooled SD: OPC {o.std():6.2f}  counterion {c.std():6.2f}  ratio {o.std()/c.std():.2f}')

# ---- R1.6 metrics in a SHARED (angle, Rg) grid per system ----
def grid(s):
    A=np.concatenate([D[(c,s,st)][0].ravel() for c in CONDS for st in ['bound','unbound']])
    G=np.concatenate([D[(c,s,st)][1].ravel() for c in CONDS for st in ['bound','unbound']])
    return (np.linspace(A.min(),A.max(),61), np.linspace(G.min(),G.max(),61))
GR={s:grid(s) for s in SYS}
def hist(a,g,s):
    xe,ye=GR[s]
    H,_,_=np.histogram2d(a.ravel(),g.ravel(),bins=[xe,ye])
    return H/H.sum()
def metrics(Hb,Hu):
    p,q=Hb.ravel(),Hu.ravel()
    ov=np.minimum(p,q).sum()
    m=0.5*(p+q); nz=lambda x,y: np.sum(x[(x>0)&(y>0)]*np.log2(x[(x>0)&(y>0)]/y[(x>0)&(y>0)]))
    js=0.5*nz(p,m)+0.5*nz(q,m)
    kl_ub=nz(q,p)
    # 3 kT occupied area: cells with F < 3kT above min, F = -ln P
    def area(H):
        h=H[H>0]; F=-np.log(h); return (F<F.min()+3).sum()/H.size
    return dict(overlap=ov,JS=js,KL_unbound_bound=kl_ub,
                area3kT_bound=area(Hb),area3kT_unbound=area(Hu))
print()
print('='*98)
print('REVIEWER 1 COMMENT 6 METRICS -- bound vs unbound, shared (angle, Rg) grid per system')
print('61x61 bins spanning all conditions and states, so values are comparable ACROSS conditions.')
print('Errors: bootstrap over 100 replicas, 200 resamples.')
print('='*98)
out={}
print(f'\n{"cond":11s} {"sys":6s} {"overlap":>14s} {"JS (bits)":>14s} {"KL(u||b)":>14s} '
      f'{"area3kT bnd":>13s} {"area3kT unb":>13s}')
for c in CONDS:
    for s in SYS:
        ab,gb=D[(c,s,'bound')]; au,gu=D[(c,s,'unbound')]
        pt=metrics(hist(ab,gb,s),hist(au,gu,s))
        bs=[]
        for _ in range(200):
            i=rng.integers(0,len(ab),len(ab)); j=rng.integers(0,len(au),len(au))
            bs.append(metrics(hist(ab[i],gb[i],s),hist(au[j],gu[j],s)))
        se={k:np.std([b[k] for b in bs],ddof=1) for k in pt}
        out[f'{c}|{s}']={'point':pt,'se':se}
        print(f'{c:11s} {s:6s} '+'  '.join(
            f'{pt[k]:7.3f}+/-{se[k]:5.3f}' for k in
            ['overlap','JS','KL_unbound_bound','area3kT_bound','area3kT_unbound']))
json.dump(out,open(f'{S}/cv_metrics.json' if (S:=__import__("os").path.dirname(__file__)) else 'cv_metrics.json','w'),indent=1)
print()
print('Water-model contrast (OPC minus TIP3P), with bootstrap SEs combined in quadrature:')
for s in SYS:
    a=out[f'opc-salt|{s}']; b=out[f'tip3p|{s}']
    for k in ['overlap','JS','area3kT_unbound']:
        d=a['point'][k]-b['point'][k]; e=np.hypot(a['se'][k],b['se'][k])
        print(f'  {s:6s} {k:17s} {d:+7.3f} +/- {e:5.3f}  ({abs(d)/e:4.1f}s) '
              f'{"RESOLVED" if abs(d)>2*e else "not resolved"}')
