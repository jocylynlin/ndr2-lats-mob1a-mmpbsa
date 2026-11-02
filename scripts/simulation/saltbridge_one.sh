#!/bin/bash
# saltbridge_one.sh <condition> <system> <slot>
# Interface salt bridges from the complex trajectory: MOB1A (:1-181) acidic O vs kinase (:182-240)
# basic N, and the reverse. Total contacts per frame = native + nonnative at a 4.0 A cutoff.
# GLH is protonated and therefore neutral, so it is excluded. Frames match MM-PBSA (1..4000 stride 10).
. "$(dirname "$0")/lib.sh"
COND=$1; SYS=$2; S=$(printf %03d $3)
TRAJ=$CAMPAIGN/trajectories/$COND/${SYS}_complex/${SYS}_complex_${S}.nc
PRM=$CAMPAIGN/topologies/$COND/${SYS}_complex.prmtop
OUT=$CAMPAIGN/analysis/saltbridge/$COND; mkdir -p "$OUT"
CSV=$OUT/${SYS}_${S}.csv
[ -s "$CSV" ] && exit 0
[ -f "$TRAJ" ] && [ -f "$PRM" ] || { echo "missing $TRAJ"; exit 1; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$CPPTRAJ -p "$PRM" > "$T/log" 2>&1 <<CP
trajin $TRAJ 1 $PROD_FRAMES 10
nativecontacts :1-181&:ASP,GLU@OD1,OD2,OE1,OE2 :182-240&:ARG,LYS@NE,NH1,NH2,NZ distance 4.0 out $T/mob_neg.dat resout $T/mob_neg_res.dat
nativecontacts :1-181&:ARG,LYS@NE,NH1,NH2,NZ :182-240&:ASP,GLU@OD1,OD2,OE1,OE2 distance 4.0 out $T/mob_pos.dat resout $T/mob_pos_res.dat
go
CP
$PYTHON - "$T" "$CSV" <<'PY'
import sys, numpy as np
T, out = sys.argv[1:3]
def tot(f):
    v = np.loadtxt(f, skiprows=1)
    return v[:,1] + v[:,2]          # native + nonnative = all contacts at the cutoff
a, b = tot(f"{T}/mob_neg.dat"), tot(f"{T}/mob_pos.dat")
with open(out,"w") as fh:
    fh.write("frame,mob1a_acidic_kinase_basic,mob1a_basic_kinase_acidic,total\n")
    for i,(x,y) in enumerate(zip(a,b),1):
        fh.write(f"{i},{x:.0f},{y:.0f},{x+y:.0f}\n")
print(f"{out}: {len(a)} frames, mean total {np.mean(a+b):.2f}")
PY
