#!/bin/bash
# dihedral_one.sh <condition> <system> <state> <slot>
#   state = bound   -> kinase phi/psi taken from the complex trajectory (strip :1-181)
#         = unbound -> kinase phi/psi from its own trajectory
# Frames 1..4000 stride 10 (400 per replica), matching the MM-PBSA sampling.
. "$(dirname "$0")/lib.sh"
COND=$1; SYS=$2; STATE=$3; S=$(printf %03d $4)
OUT=$CAMPAIGN/analysis/dihedral/$COND/$STATE; mkdir -p "$OUT"
DAT=$OUT/${SYS}_${S}.dat
[ -s "$DAT" ] && exit 0
if [ "$STATE" = bound ]; then
  PRM=$CAMPAIGN/topologies/$COND/${SYS}_complex.prmtop
  TRJ=$CAMPAIGN/trajectories/$COND/${SYS}_complex/${SYS}_complex_${S}.nc
  STRIP="strip :1-181"
else
  PRM=$CAMPAIGN/topologies/$COND/${SYS}_kinase.prmtop
  TRJ=$CAMPAIGN/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${S}.nc
  STRIP=""
fi
[ -f "$TRJ" ] || { echo "missing $TRJ"; exit 1; }
$CPPTRAJ -p "$PRM" >/dev/null 2>&1 <<CP
trajin $TRJ 1 $PROD_FRAMES 10
$STRIP
multidihedral phi psi out $DAT
go
CP
[ -s "$DAT" ] || echo "FAILED $COND/$SYS/$STATE/$S"
