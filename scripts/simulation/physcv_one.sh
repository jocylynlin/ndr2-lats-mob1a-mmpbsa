#!/bin/bash
# physcv_one.sh <condition> <system> <state: bound|unbound> <slot>
# Physical collective variables for the kinase fragment, no fitting involved:
#   angle : interhelical angle between helix A (res 18-31) and helix B (res 38-57) CA axes
#   rg    : radius of gyration of the fragment, heavy atoms
#   rmsd  : CA RMSD to a single common bound-state reference per system
. "$(dirname "$0")/lib.sh"
COND=$1; SYS=$2; STATE=$3; S=$(printf %03d $4)
OUT=$CAMPAIGN/analysis/physcv/$COND/$STATE; mkdir -p "$OUT"
CSV=$OUT/${SYS}_${S}.csv
[ -s "$CSV" ] && exit 0
REF=$CAMPAIGN/refs/ref_${SYS}.pdb
if [ "$STATE" = unbound ]; then
  TRJ=$CAMPAIGN/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${S}.nc
  PRM=$CAMPAIGN/topologies/$COND/${SYS}_kinase.prmtop; STRIP=""
else
  TRJ=$CAMPAIGN/trajectories/$COND/${SYS}_complex/${SYS}_complex_${S}.nc
  PRM=$CAMPAIGN/topologies/$COND/${SYS}_complex.prmtop; STRIP="strip :1-181"
fi
[ -f "$TRJ" ] || exit 1
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$CPPTRAJ -p "$PRM" >/dev/null 2>&1 <<CP
parm $REF name REF
reference $REF parm REF name BOUND
trajin $TRJ 1 $PROD_FRAMES 10
$STRIP
vector v1 :18@CA :31@CA
vector v2 :38@CA :57@CA
vectormath vec1 v1 vec2 v2 dotangle out $T/ang.dat
radgyr :1-59&!@H= out $T/rg.dat
rms toBound ref BOUND :1-59@CA out $T/rms.dat
go
CP
paste <(awk 'NR>1{print $2}' $T/ang.dat) <(awk 'NR>1{print $2}' $T/rg.dat) \
      <(awk 'NR>1{print $2}' $T/rms.dat) | awk 'BEGIN{print "frame,angle,rg,rmsd"}
      {printf "%d,%.3f,%.4f,%.4f\n",NR,$1,$2,$3}' > "$CSV"
