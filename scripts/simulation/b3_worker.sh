#!/bin/bash
# b3_worker.sh SY ST N -- bonded energy of the full fragment and the core, one replica.
set -u
SY=$1; ST=$2; n=$3
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:$C/lib
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
O=$C/analysis/b3_terminal_all; mkdir -p "$O"
[ -s "$O/${SY}_${ST}_${n}.dat" ] && [ -s "$O/${SY}_${ST}_core_${n}.dat" ] && exit 0
if [ "$ST" = unbound ]; then
  PRM=$C/topologies/opc-salt/${SY}_kinase.prmtop; STRIP=""
  TRJD=$C/trajectories/opc-salt/${SY}_kinase; PFX=${SY}_kinase
else
  PRM=$C/topologies/opc-salt/${SY}_complex.prmtop; STRIP="strip :1-181"
  TRJD=$C/trajectories/opc-salt/${SY}_complex; PFX=${SY}_complex
fi
[ -s "$TRJD/${PFX}_${n}.nc" ] || exit 0
$CP -p "$PRM" >/dev/null 2>&1 <<CPT
trajin $TRJD/${PFX}_${n}.nc 1 4000 10
$STRIP
energy FULL :1-59 out $O/${SY}_${ST}_$n.dat
energy CORE :6-54 out $O/${SY}_${ST}_core_$n.dat
go
CPT
