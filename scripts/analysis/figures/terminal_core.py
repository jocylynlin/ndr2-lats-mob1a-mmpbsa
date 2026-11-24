#!/usr/bin/env python3
"""Terminal-vs-core covalent decomposition of the kinase, final dataset.

Data   reviewer_comments/analysis_terminal_core/data/{sys}_{state}_{NNN}_{full,core}.dat
       produced by analysis_terminal_core/scripts/run_terminal.sh (see that file for the
       trajectory and topology provenance: bound from the unchanged complex boxes, unbound
       from the matched-salt kinase boxes).

Partition   FULL = residues 1-59, CORE = residues 6-54, TERMINI = FULL - CORE
            (residues 1-5 and 55-59; ten residues in all).
            Bonds crossing the core boundary are missing from CORE, but that offset is the
            same in bound and unbound, so it cancels in the bound-minus-unbound difference.
Covalent    bond + angle + dih + vdw14 + elec14, the manuscript's covalent group.
Statistics  replica is the unit; 100 bound and 100 unbound replicas are INDEPENDENT
            simulation sets, so each quantity's SE is the two sets' replica-level SEs in
            quadrature, never a paired difference. Each of FULL, CORE and TERMINI is
            computed separately rather than by subtracting rounded parts, so SE(FULL) is
            not SE(CORE) + SE(TERMINI): the two shares are correlated within a frame.
"""
import glob, math, os, re, sys
import numpy as np

DATA='/home3/jocelyn/projects/hippo/reviewer_comments/analysis_terminal_core/data'
COVCOLS=['bond','angle','dih','vdw14','elec14']
SYS=['lats1','ndr2']
se=lambda a: float(np.std(a,ddof=1)/math.sqrt(len(a)))

def cov_of(path, tag):
    """Mean covalent energy over the frames of one replica."""
    hdr=None; tot=None; n=0
    for line in open(path):
        if line.startswith('#'):
            hdr=line.lstrip('#').split(); continue
        v=line.split()
        if len(v)<2: continue
        if tot is None: tot=0.0
        idx=[hdr.index(f'{tag}[{c}]') for c in COVCOLS]
        tot+=sum(float(v[i]) for i in idx); n+=1
    return tot/n if n else float('nan')

def collect():
    """-> {(sys,state,part): array over replicas}"""
    out={}
    for s in SYS:
        for st in ('bound','unbound'):
            full=[]; core=[]; reps=[]
            for f in sorted(glob.glob(f'{DATA}/{s}_{st}_[0-9][0-9][0-9]_full.dat')):
                n=re.search(r'_(\d{3})_full',f).group(1)
                c=f.replace('_full.dat','_core.dat')
                if not os.path.exists(c): continue
                full.append(cov_of(f,'FULL')); core.append(cov_of(c,'CORE')); reps.append(n)
            full=np.array(full); core=np.array(core)
            out[(s,st,'full')]=full; out[(s,st,'core')]=core
            out[(s,st,'termini')]=full-core
            out[(s,st,'_n')]=len(reps)
    return out

def report(D):
    res={}
    print(f"{'system':6s} {'part':8s} {'bound':>9s} {'unbound':>9s} {'bound-unbound':>22s}")
    for s in SYS:
        for part in ('full','core','termini'):
            b,u=D[(s,'bound',part)], D[(s,'unbound',part)]
            m=b.mean()-u.mean(); e=math.hypot(se(b),se(u))
            res[(s,part)]=(m,e)
            print(f"{s:6s} {part:8s} {b.mean():9.2f} {u.mean():9.2f}   {m:+8.2f} +- {e:.2f}")
    print()
    for part in ('full','core','termini'):
        bl,ul=D[('lats1','bound',part)],D[('lats1','unbound',part)]
        bn,un=D[('ndr2','bound',part)], D[('ndr2','unbound',part)]
        m=(bl.mean()-ul.mean())-(bn.mean()-un.mean())
        e=math.hypot(se(bl),se(ul),se(bn),se(un))
        res[('diff',part)]=(m,e)
        print(f"LATS1 minus NDR2, {part:8s}: {m:+.2f} +- {e:.2f}")
    print()
    for s in SYS:
        tot=res[(s,'full')][0]
        if tot: print(f"{s}: termini are {100*res[(s,'termini')][0]/tot:.0f}% of the kinase covalent change")
    return res

if __name__=='__main__':
    D=collect()
    for s in SYS:
        for st in ('bound','unbound'):
            print(f"  {s}_{st}: {D[(s,st,'_n')]} replicas")
    print()
    report(D)
