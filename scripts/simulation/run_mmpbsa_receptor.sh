#!/bin/bash
# run_mmpbsa_receptor.sh <condition> [first] [last] [systems]
#   PARALLEL=n   replicas computed at once (default 4)
#
# Absolute (stability) PB energy of the UNBOUND MOB1A ensemble, scored with the complex-derived
# receptor topology - exactly the receptor leg of the MT calculation, run on its own so the two
# systems' unbound MOB1A ensembles can be compared without waiting for the kinase trajectories.
# Settings are byte-identical to scripts/run_mmpbsa.sh so the numbers are directly comparable to
# the rec_* columns already produced by ST and MT.
. "$(dirname "$0")/lib.sh"
COND=$1; FIRST=${2:-1}; LAST=${3:-$N_TARGET}; SYSLIST=${4:-"ndr2 lats1"}
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
PARALLEL=${PARALLEL:-4}
TOP=$CAMPAIGN/topologies/$COND; OUT=$CAMPAIGN/mmpbsa/$COND/rec; mkdir -p "$OUT"
LOG=$CAMPAIGN/logs/mmpbsa_${COND}_rec_$(date +%Y%m%d_%H%M).log
cat > "$OUT/mmpbsa.in" <<IN
Stability (unbound MOB1A) PB per-frame calculation, $COND, first $PROD_FRAMES frames = $PROD_NS ns
&general
startframe=1, endframe=$PROD_FRAMES, interval=10,
verbose=2, keep_files=0,strip_mask=":WAT,CL,CIO,CS,IB,K,LI,MG,RB"
/
&pb
istrng=0.150, indi=20.0,exdi=80.0,
/
IN

run_one() {
  local SYS=$1 s=$2 n d RT
  n=$(printf %03d "$s"); d=$OUT/${SYS}_$n; mkdir -p "$d"
  # already complete? (stability output has no DELTA section, so key on TOTAL + per-frame file)
  if grep -aq "^TOTAL" "$d/FINAL_RESULTS_MMPBSA.dat" 2>/dev/null && [ -s "$d/per_frame.csv" ]; then return 0; fi
  rm -f "$d/FINAL_RESULTS_MMPBSA.dat"
  RT=$CAMPAIGN/trajectories/$COND/${SYS}_mob1a/${SYS}_mob1a_${n}.nc
  [ -f "$RT" ] || { echo "$(date +%T) $SYS $n: no mob1a trajectory" >> "$LOG"; return 1; }
  # pre-flight: the unbound trajectory must match the complex-derived receptor topology exactly
  $PYTHON "$CAMPAIGN/scripts/check_traj_topology.py" "$RT" "$TOP/${SYS}_receptor.prmtop" >> "$LOG" 2>&1 || {
    echo "$(date +%T) $SYS $n: PRE-FLIGHT FAILED (receptor)" >> "$LOG"; return 1; }
  local T=/tmp/$USER/mmpbsa_${COND}_rec_${SYS}_${n}_$$
  mkdir -p "$T" || { echo "$(date +%T) $SYS $n: cannot create scratch $T" >> "$LOG"; return 1; }
  ( cd "$T" && MMPBSA.py -O -i "$OUT/mmpbsa.in" -o FINAL_RESULTS_MMPBSA.dat -eo per_frame.csv \
      -cp "$TOP/${SYS}_receptor.prmtop" -y "$RT" >> "$LOG" 2>&1 )
  for f in FINAL_RESULTS_MMPBSA.dat per_frame.csv; do [ -f "$T/$f" ] && cp -p "$T/$f" "$d/"; done
  rm -rf "$T"
  if grep -aq "^TOTAL" "$d/FINAL_RESULTS_MMPBSA.dat" 2>/dev/null && [ -s "$d/per_frame.csv" ]; then
    echo "$(date +%T) $SYS $n done ($(wc -l < "$d/per_frame.csv") lines)" >> "$LOG"
  else
    echo "$(date +%T) $SYS $n FAILED or incomplete" >> "$LOG"
    mv -f "$d/FINAL_RESULTS_MMPBSA.dat" "$d/FINAL_RESULTS_MMPBSA.dat.incomplete" 2>/dev/null
  fi
}
export -f run_one
export OUT TOP LOG CAMPAIGN COND PYTHON PROD_FRAMES AMBERHOME PATH

echo "$(date '+%F %T')  ===== unbound-MOB1A stability MM-PBSA, $COND, systems: $SYSLIST, replicas $FIRST-$LAST, $PARALLEL at a time" | tee -a "$LOG"
for SYS in $SYSLIST; do
  seq "$FIRST" "$LAST" | xargs -P "$PARALLEL" -I{} bash -c "run_one $SYS {}"
  ok=$(grep -alc "^TOTAL" $OUT/${SYS}_*/FINAL_RESULTS_MMPBSA.dat 2>/dev/null | wc -l)
  echo "$(date '+%F %T')  $SYS: $ok replica directories have complete output" | tee -a "$LOG"
done
