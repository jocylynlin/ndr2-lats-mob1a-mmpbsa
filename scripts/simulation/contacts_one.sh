#!/bin/bash
# SUPERSEDED 2026-10-20 — resoffset is silently ignored without `byresidue`; this metric has
# NO sequence-separation filter and is 94-96% local contacts. Use contacts_tertiary.sh instead.
# Retained for provenance only. See analysis/contacts/INVALID_DO_NOT_USE.txt
# contacts_one.sh <condition> <system> <slot>
# Intramolecular heavy-atom contacts in the UNBOUND kinase, two shells, |i-j|>=3:
#   N(<5.0 A) and N(<3.4 A).  packing shell = N(5.0)-N(3.4) (attractive region)
#                             close contacts = N(3.4)       (potentially repulsive)
. "$(dirname "$0")/lib.sh"
COND=$1; SYS=$2; S=$(printf %03d $3)
OUT=$CAMPAIGN/analysis/contacts/$COND; mkdir -p "$OUT"
CSV=$OUT/${SYS}_${S}.csv
[ -s "$CSV" ] && exit 0
TRJ=$CAMPAIGN/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${S}.nc
PRM=$CAMPAIGN/topologies/$COND/${SYS}_kinase.prmtop
[ -f "$TRJ" ] || { echo "missing $TRJ"; exit 1; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$CPPTRAJ -p "$PRM" >/dev/null 2>&1 <<CP
trajin $TRJ 1 $PROD_FRAMES 10
nativecontacts !@H=,H out $T/far.dat  distance 5.0 resoffset 3
nativecontacts !@H=,H out $T/near.dat distance 3.4 resoffset 3
go
CP
$PYTHON - "$T" "$CSV" <<'PY'
import sys,numpy as np
T,out=sys.argv[1:3]
f=np.loadtxt(f"{T}/far.dat",skiprows=1); n=np.loadtxt(f"{T}/near.dat",skiprows=1)
tf=f[:,1]+f[:,2]; tn=n[:,1]+n[:,2]
with open(out,"w") as fh:
    fh.write("frame,close_lt3.4,shell_3.4to5.0,total_lt5.0\n")
    for i,(a,b) in enumerate(zip(tn,tf),1): fh.write(f"{i},{a:.0f},{b-a:.0f},{b:.0f}\n")
PY
