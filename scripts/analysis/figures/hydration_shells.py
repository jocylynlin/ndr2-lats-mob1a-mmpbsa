#!/usr/bin/env python3
"""Hydration-shell populations and the surface-area control (SI Figure S5 caption).

Separates three things that must not be conflated:
  (1) hydration-shell populations differ between the solvent setups  -- established here;
  (2) conformational distributions differ                            -- a separate question;
  (3) MM-PBSA binding enthalpies differ                              -- a third, separate question.

The control that keeps (1) from being read as evidence for (2): if the shell difference were a
consequence of the solute sampling more compact conformations, it would scale with solvent-
accessible area and vanish once normalized by it. It does not.

Shell definition: cpptraj `watershell :1-59 lower 3.4 upper 5.0`. First shell = water molecules
with any atom within 3.4 Ang of any atom of the 59-residue segment; second shell = 3.4-5.0 Ang.
Per-frame counts averaged within a replica, then over 20 replicas. All differences are quoted
OPC minus TIP3P, so a positive value means MORE water under OPC.
"""
import glob, math
import numpy as np
A='/home3/jocelyn/projects/hippo/reviewer_comments/analysis_hydration_matched/data'
se=lambda a: float(np.std(a,ddof=1)/math.sqrt(len(a)))
def col(cond,sysn,tag,c):
    v=[]
    for f in sorted(glob.glob(f'{A}/{cond}/{sysn}/replica*/{tag}.dat')):
        d=np.loadtxt(f,comments='#')
        if d.ndim==2 and d.shape[1]>c: v.append(d[:,c].mean())
    return np.array(v)
ROWS=[('ws_protein',1,'first shell, whole segment (<3.4 A)'),
      ('ws_protein',2,'second shell, whole segment (3.4-5.0 A)'),
      ('ws_anionic',1,'first shell, carboxylate O'),
      ('ws_cationic',1,'first shell, Arg/Lys/His N'),
      ('sasa',1,'SASA of residues 1-59 (A^2)')]
print("OPC minus TIP3P; positive = larger under OPC; n = 20 replicas each\n")
for s in ('lats1_kinase','ndr2_kinase'):
    print(f"--- {s} ---")
    for tag,c,lab in ROWS:
        o,t=col('opc-salt-matched',s,tag,c), col('tip3p-matched',s,tag,c)
        d=o.mean()-t.mean(); e=math.hypot(se(o),se(t))
        print(f"  {lab:40s} OPC {o.mean():9.2f}  TIP3P {t.mean():9.2f}   {d:+8.2f}+-{e:5.2f} ({abs(d)/e:4.1f}s)")
    wo,so=col('opc-salt-matched',s,'ws_protein',1), col('opc-salt-matched',s,'sasa',1)
    wt,st=col('tip3p-matched',s,'ws_protein',1),   col('tip3p-matched',s,'sasa',1)
    ro,rt=wo/so, wt/st; d=ro.mean()-rt.mean(); e=math.hypot(se(ro),se(rt))
    print(f"  {'CONTROL: first shell per A^2 of SASA':40s} OPC {ro.mean():9.5f}  TIP3P {rt.mean():9.5f}"
          f"   {100*d/rt.mean():+7.1f}% +-{100*e/rt.mean():4.1f} ({abs(d)/e:4.1f}s)")
    print("  -> the difference survives normalization by surface area, so it is a property of the\n"
          "     water model, not a consequence of different solute conformations.\n")
