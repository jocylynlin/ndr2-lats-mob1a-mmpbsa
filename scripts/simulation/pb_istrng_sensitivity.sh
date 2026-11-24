#!/bin/bash
# pb_istrng_sensitivity.sh <condition> <istrng> [nrep] [interval]
#
# PB rescoring only -- no new MD. Re-evaluates existing stripped trajectories at a different PB
# ionic strength so the istrng=0.150 choice can be shown not to drive the result. Paired design:
# the SAME frames are scored at every istrng, so frame-to-frame noise cancels in the comparison.
#
# Mirrors campaign_100x40ns/scripts/run_mmpbsa.sh exactly except for (a) istrng, (b) a coarser
# frame interval, (c) a separate output tree so nothing published is touched.
#
# REFERENCE CONVENTION: like run_mmpbsa.sh, the MT arm passes the OWN-system MOB1A as -yr, so the
# raw output is index-matched, NOT pooled. The pooled reference the paper uses is applied in
# post-processing from the rec_* columns of per_frame.csv (grand mean over both systems), which is
# how Tables 1-3 were reproduced exactly from the index-matched per-frame files. Report ddH under
# the pooled convention; mixing conventions can flip the sign.
set -u
COND=${1:?condition}; ISTR=${2:?istrng}; NREP=${3:-10}; INTERVAL=${4:-40}; MODE=${5:-mt}
SRC=/home3/jocelyn/projects/hippo/campaign_100x40ns
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
PARALLEL=${PARALLEL:-8}
TOP=$SRC/topologies/$COND
OUT=/home3/jocelyn/projects/hippo/campaign_saltmatched/pb_sensitivity/$COND/$MODE/istrng_$ISTR
mkdir -p "$OUT"
LOG=$OUT/run.log
cat > "$OUT/mmpbsa.in" <<IN
PB rescoring sensitivity, $COND, istrng=$ISTR, interval=$INTERVAL
&general
startframe=1, endframe=400, interval=$INTERVAL,
verbose=2, keep_files=0,strip_mask=":WAT,CL,CIO,CS,IB,K,LI,MG,RB"
/
&pb
istrng=$ISTR, indi=20.0,exdi=80.0,
/
IN

run_one() {
  local SYS=$1 s=$2 n d CT RT LT T
  n=$(printf %03d "$s"); d=$OUT/${SYS}_$n; mkdir -p "$d"
  [ -s "$d/per_frame.csv" ] && return 0
  CT=$SRC/trajectories/$COND/${SYS}_complex/${SYS}_complex_${n}.nc
  RT=$SRC/trajectories/$COND/${SYS}_mob1a/${SYS}_mob1a_${n}.nc
  LT=$SRC/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${n}.nc
  local NEED="$CT"; [ "$MODE" = mt ] && NEED="$CT $RT $LT"
  for f in $NEED; do [ -f "$f" ] || { echo "$SYS $n missing $(basename $f)" >>"$LOG"; return 1; }; done
  # ST: receptor and ligand are taken from the complex trajectory, so the unbound runs -- and the
  # salt bug that affects only them -- play no part. ST sensitivity is therefore final-quality now.
  # MT: receptor from the unbound MOB1A run, ligand from the unbound kinase. Those unbound
  # trajectories are the over-salted ones being replaced, so MT here is a first look only.
  local ARGS="-y $CT"; [ "$MODE" = mt ] && ARGS="-y $CT -yr $RT -yl $LT"
  T=/tmp/$USER/pbsens_${COND}_${MODE}_${ISTR}_${SYS}_${n}_$$; mkdir -p "$T"
  ( cd "$T" && MMPBSA.py -O -i "$OUT/mmpbsa.in" -o FINAL_RESULTS_MMPBSA.dat -eo per_frame.csv \
      -cp "$TOP/${SYS}_complex.prmtop" -rp "$TOP/${SYS}_receptor.prmtop" -lp "$TOP/${SYS}_ligand.prmtop" \
      $ARGS >>"$LOG" 2>&1 )
  cp -f "$T/per_frame.csv" "$T/FINAL_RESULTS_MMPBSA.dat" "$d/" 2>/dev/null
  rm -rf "$T"
  [ -s "$d/per_frame.csv" ] && echo "$(date +%T) $SYS $n done (istrng=$ISTR)" >>"$LOG" \
    || echo "$(date +%T) $SYS $n FAILED" >>"$LOG"
}

j=0
for SYS in lats1 ndr2; do
  for ((s=1; s<=NREP; s++)); do
    run_one "$SYS" "$s" &
    j=$((j+1)); [ $((j % PARALLEL)) -eq 0 ] && wait
  done
done
wait
n=$(find "$OUT" -name per_frame.csv -size +0 2>/dev/null | wc -l)
want=$((NREP*2))
if [ "$n" -eq "$want" ]; then
  echo "PB_SENS_DONE $COND $MODE istrng=$ISTR ($n/$want per_frame.csv)"
else
  echo "PB_SENS_INCOMPLETE $COND $MODE istrng=$ISTR ($n/$want per_frame.csv) -- see $OUT/run.log"
  exit 1
fi
