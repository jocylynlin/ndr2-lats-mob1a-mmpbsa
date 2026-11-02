#!/bin/bash
# late_mmpbsa_worker.sh MODE SYS N -- one MM-PBSA calculation of the 40-100 ns window.
# Idempotent: returns immediately if per_frame.csv already exists.
set -u
MODE=$1; SYS=$2; n=$3
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
# MMPBSA.py shells out to amber25 cpptraj, which links libreadline.so.6 / libncurses.so.5 /
# libtinfo.so.5 -- present only in $C/lib on this image. Omitting this makes every call fail
# with "TrajError: cpptraj failed when querying <traj>".
export LD_LIBRARY_PATH=$C/lib:${LD_LIBRARY_PATH:-}
TOP=$C/topologies/opc-salt; T100=$C/trajectories100/opc-salt
OUT=$C/mmpbsa/opc-salt/${MODE}_late; d=$OUT/${SYS}_$n
mkdir -p "$d"
[ -s "$d/per_frame.csv" ] && { echo "  SKIP $MODE $SYS $n (already done)"; exit 0; }
CT=$C/trajectories/opc-salt/${SYS}_complex/${SYS}_complex_${n}.nc
RT=$C/trajectories/opc-salt/${SYS}_mob1a/${SYS}_mob1a_${n}.nc
LT=$T100/${SYS}_kinase/${SYS}_kinase_${n}.nc
[ -s "$CT" ] || { echo "  NO complex $SYS $n"; exit 1; }
if [ "$MODE" = mt ]; then
  [ -s "$RT" ] && [ -s "$LT" ] || { echo "  NO unbound $SYS $n"; exit 1; }
  ARGS="-y $CT -yr $RT -yl $LT"
else ARGS="-y $CT"; fi
W=/tmp/$USER/latepbsa_${MODE}_${SYS}_${n}_$$; rm -rf $W; mkdir -p $W
( cd $W && MMPBSA.py -O -i $OUT/mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -eo per_frame.csv \
    -cp $TOP/${SYS}_complex.prmtop -rp $TOP/${SYS}_receptor.prmtop -lp $TOP/${SYS}_ligand.prmtop \
    $ARGS > mmpbsa.stdout 2>&1 )
for f in FINAL_RESULTS_MMPBSA.dat per_frame.csv mmpbsa.stdout; do [ -f "$W/$f" ] && cp -p "$W/$f" "$d/"; done
rm -rf $W
[ -s "$d/per_frame.csv" ] && echo "  $MODE $SYS $n done" || echo "  $MODE $SYS $n FAILED"
