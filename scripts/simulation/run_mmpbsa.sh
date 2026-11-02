#!/bin/bash
# run_mmpbsa.sh <condition> <mode: st|mt> [first] [last] [systems]
#   PARALLEL=n  number of replicas computed at once (default 4)
#
# ST: complex, receptor and ligand energies all evaluated from the complex trajectory.
# MT: complex from the complex trajectory, receptor from the unbound MOB1A run, ligand from the unbound kinase.
# Always the first PROD_FRAMES frames with interval=10 -> 400 snapshots at 100 ps, so reused 100 ns
# trajectories contribute exactly their first 40 ns. Per-frame energies are written for every replica
# (-eo per_frame.csv, verbose=2) and every original MM-PBSA output file is kept under
# mmpbsa/<condition>/<mode>/<system>_<replica>/.
. "$(dirname "$0")/lib.sh"
COND=$1; MODE=$2; FIRST=${3:-1}; LAST=${4:-$N_TARGET}; SYSLIST=${5:-"lats1 ndr2"}
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
PARALLEL=${PARALLEL:-4}
TOP=$CAMPAIGN/topologies/$COND; OUT=$CAMPAIGN/mmpbsa/$COND/$MODE; mkdir -p "$OUT"
LOG=$CAMPAIGN/logs/mmpbsa_${COND}_${MODE}_$(date +%Y%m%d_%H%M).log
cat > "$OUT/mmpbsa.in" <<IN
Input file for PB per-frame calculations ($MODE), $COND, first $PROD_FRAMES frames = $PROD_NS ns
&general
startframe=1, endframe=$PROD_FRAMES, interval=10,
verbose=2, keep_files=0,strip_mask=":WAT,CL,CIO,CS,IB,K,LI,MG,RB"
/
&pb
istrng=0.150, indi=20.0,exdi=80.0,
/
IN

run_one() {                     # $1 system  $2 slot
  local SYS=$1 s=$2 n d CT RT LT ARGS
  n=$(printf %03d "$s"); d=$OUT/${SYS}_$n; mkdir -p "$d"
  if [ -f "$d/FINAL_RESULTS_MMPBSA.dat" ] && grep -aq "DELTA TOTAL" "$d/FINAL_RESULTS_MMPBSA.dat" \
     && [ -s "$d/per_frame.csv" ]; then return 0; fi
  rm -f "$d/FINAL_RESULTS_MMPBSA.dat"
  CT=$CAMPAIGN/trajectories/$COND/${SYS}_complex/${SYS}_complex_${n}.nc
  RT=$CAMPAIGN/trajectories/$COND/${SYS}_mob1a/${SYS}_mob1a_${n}.nc
  LT=$CAMPAIGN/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${n}.nc
  [ -f "$CT" ] || { echo "$(date +%T) $SYS $n: no complex trajectory" >> "$LOG"; return 1; }
  if [ "$MODE" = mt ]; then
    [ -f "$RT" ] && [ -f "$LT" ] || { echo "$(date +%T) $SYS $n: unbound trajectories missing" >> "$LOG"; return 1; }
    ARGS="-y $CT -yr $RT -yl $LT"
  else
    ARGS="-y $CT"
  fi
  # pre-flight: atom count AND ordering must match the topology each trajectory will be scored with
  $PYTHON "$CAMPAIGN/scripts/check_traj_topology.py" "$CT" "$TOP/${SYS}_complex.prmtop" >> "$LOG" 2>&1 || {
    echo "$(date +%T) $SYS $n: PRE-FLIGHT FAILED (complex)" >> "$LOG"; return 1; }
  if [ "$MODE" = mt ]; then
    $PYTHON "$CAMPAIGN/scripts/check_traj_topology.py" "$RT" "$TOP/${SYS}_receptor.prmtop" >> "$LOG" 2>&1 || {
      echo "$(date +%T) $SYS $n: PRE-FLIGHT FAILED (receptor)" >> "$LOG"; return 1; }
    $PYTHON "$CAMPAIGN/scripts/check_traj_topology.py" "$LT" "$TOP/${SYS}_ligand.prmtop" >> "$LOG" 2>&1 || {
      echo "$(date +%T) $SYS $n: PRE-FLIGHT FAILED (ligand)" >> "$LOG"; return 1; }
  fi
  # run in node-local scratch: MMPBSA writes large temp files (_MMPBSA_*.mdcrd) and a ~350 MB reference.frc
  # force dump that nothing reads. Keeping them off NFS removes the main I/O bottleneck and the wasted space.
  local T=/tmp/$USER/mmpbsa_${COND}_${MODE}_${SYS}_${n}_$$
  mkdir -p "$T" || { echo "$(date +%T) $SYS $n: cannot create scratch $T" >> "$LOG"; return 1; }
  ( cd "$T" && MMPBSA.py -O -i "$OUT/mmpbsa.in" -o FINAL_RESULTS_MMPBSA.dat -eo per_frame.csv \
      -cp "$TOP/${SYS}_complex.prmtop" -rp "$TOP/${SYS}_receptor.prmtop" -lp "$TOP/${SYS}_ligand.prmtop" \
      $ARGS >> "$LOG" 2>&1 )
  # copy back only what is needed for analysis; the force dump and mdcrd scratch stay behind and are removed
  for f in FINAL_RESULTS_MMPBSA.dat per_frame.csv; do [ -f "$T/$f" ] && cp -p "$T/$f" "$d/"; done
  rm -rf "$T"
  if grep -aq "DELTA TOTAL" "$d/FINAL_RESULTS_MMPBSA.dat" 2>/dev/null && [ -s "$d/per_frame.csv" ]; then
    echo "$(date +%T) $SYS $n done ($(wc -l < "$d/per_frame.csv") lines in per_frame.csv)" >> "$LOG"
  else
    echo "$(date +%T) $SYS $n FAILED or incomplete" >> "$LOG"
    mv -f "$d/FINAL_RESULTS_MMPBSA.dat" "$d/FINAL_RESULTS_MMPBSA.dat.incomplete" 2>/dev/null
  fi
}
export -f run_one
export OUT TOP LOG CAMPAIGN COND MODE PYTHON PROD_FRAMES AMBERHOME PATH

echo "$(date '+%F %T')  ===== $MODE MM-PBSA, $COND, systems: $SYSLIST, replicas $FIRST-$LAST, $PARALLEL at a time" | tee -a "$LOG"
for SYS in $SYSLIST; do
  seq "$FIRST" "$LAST" | xargs -P "$PARALLEL" -I{} bash -c "run_one $SYS {}"
  done_n=$(ls -d $OUT/${SYS}_*/ 2>/dev/null | wc -l)
  ok_n=$(grep -alc "DELTA TOTAL" $OUT/${SYS}_*/FINAL_RESULTS_MMPBSA.dat 2>/dev/null | wc -l)
  echo "$(date '+%F %T')  $SYS: $ok_n of $done_n replica directories have complete output" | tee -a "$LOG"
done
