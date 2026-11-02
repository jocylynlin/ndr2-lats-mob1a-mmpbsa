#!/bin/bash
#SBATCH --job-name=latepar
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/late_par_%j.out
#SBATCH --error=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/late_par_%j.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=96 --time=1-00:00:00 --partition=cpu-long
# Whole node (96 cores, 128 GB): there are exactly 80 tasks, so all run concurrently and the
# wall clock is one PB calculation (~2.5 h for 600 frames on one core) rather than 80 of them.
# Taking the full node also means no other job competes for memory with 80 PB grids.
# Parallel replacement for STEP 2 of late_mmpbsa_pipeline.sh. STEP 1 (trajectory concatenation)
# is already complete in trajectories100/opc-salt -- 40 files, 10000 frames each, verified.
# The original ran 80 PB calculations serially on a 24-core allocation; this runs 24 at a time.
# Fully idempotent: any per_frame.csv already present is skipped, so this is safe to re-run.
set -u
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
S=$C/scripts
NPROC=${NPROC_OVERRIDE:-80}   # 80 tasks, run them all at once
echo "$(date '+%F %T')  STEP 2 parallel: MM-PBSA frames 4001-10000, $NPROC concurrent"
for MODE in st mt; do
  OUT=$C/mmpbsa/opc-salt/${MODE}_late; mkdir -p $OUT
  cat > $OUT/mmpbsa.in <<IN
Input file for PB per-frame calculations ($MODE), opc-salt, LATE window frames 4001-10000 = 40-100 ns
&general
startframe=4001, endframe=10000, interval=10,
verbose=2, keep_files=0,strip_mask=":WAT,CL,CIO,CS,IB,K,LI,MG,NA,RB"
/
&pb
istrng=0.150, indi=20.0,exdi=80.0,
/
IN
done
# build the full task list, both protocols together so all 24 cores stay busy to the end
for MODE in st mt; do for SYS in lats1 ndr2; do for i in $(seq 1 20); do
  printf '%s %s %03d\n' "$MODE" "$SYS" "$i"
done; done; done > $C/logs/late_tasks.txt
echo "tasks: $(wc -l < $C/logs/late_tasks.txt)"
xargs -a $C/logs/late_tasks.txt -P "$NPROC" -n 3 $S/late_mmpbsa_worker.sh
for MODE in st mt; do
  echo "  $MODE complete: $(find $C/mmpbsa/opc-salt/${MODE}_late -name per_frame.csv -size +1k | wc -l)/40"
done
echo "$(date '+%F %T')  LATEPBSADONE"
