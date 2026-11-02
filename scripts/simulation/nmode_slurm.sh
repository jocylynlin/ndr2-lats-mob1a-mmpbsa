#!/bin/bash
#SBATCH --job-name=nmode
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/nmode_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40   # override with sbatch -c N. MEMORY: the complex Hessian alone is
# (3N)^2 doubles = 1.12 GB for a 3937-atom complex, plus diagonalisation workspace; 40-way parallelism
# exhausted a 128 GB node and lost 4 of 40 replicas on 2026-09-28. Use NMODE_PARALLEL<=20.
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-long
# nmode_slurm.sh <condition> <system> <first> <last>
# Normal-mode entropy, multi-trajectory, 10 evenly spaced frames per replica ending at the last frame -
# the same protocol as nmode/moblats_opc_salt_10frames, applied to the campaign's 4000-frame trajectories
# (interval 400, startframe 400) and the campaign's complex-derived topologies.
cd /home3/jocelyn/projects/hippo/campaign_100x40ns || exit 1
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_100x40ns/lib:$AMBERHOME/lib:$LD_LIBRARY_PATH
COND=$1; SYS=$2; FIRST=${3:-1}; LAST=${4:-20}
CAMP=/home3/jocelyn/projects/hippo/campaign_100x40ns
TOP=$CAMP/topologies/$COND; OUT=$CAMP/nmode/$COND/$SYS; mkdir -p "$OUT"
LOG=$CAMP/logs/nmode_${COND}_${SYS}_$(date +%Y%m%d_%H%M).log
echo "node $(hostname), $(nproc) cores | nmode $COND/$SYS replicas $FIRST-$LAST" | tee -a "$LOG"

nmode_one() {
  local s=$1 n d ct rt lt
  n=$(printf %03d $s); d=$OUT/traj_$n; mkdir -p "$d"
  [ -f "$d/FINAL_RESULTS_MMPBSA.dat" ] && grep -aq "Entropy" "$d/FINAL_RESULTS_MMPBSA.dat" && return 0
  ct=$CAMP/trajectories/$COND/${SYS}_complex/${SYS}_complex_${n}.nc
  rt=$CAMP/trajectories/$COND/${SYS}_mob1a/${SYS}_mob1a_${n}.nc
  lt=$CAMP/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${n}.nc
  for f in "$ct" "$rt" "$lt"; do [ -f "$f" ] || { echo "$(date +%T) $n missing $f" >> "$LOG"; return 1; }; done
  cat > "$d/nmode.in" <<IN
Input file for nmode entropy calculation (Multi-Trajectory, $COND, 10 frames)
&general
startframe=400, endframe=4000, interval=400,
verbose=1, keep_files=1,
/
&pb
istrng=0.150, indi=20.0, exdi=80.0,
/
&nmode
nmstartframe=1, nmendframe=10, nminterval=1,
maxcyc=10000, drms=0.001,
/
IN
  local T=/tmp/$USER/nmode_${COND}_${SYS}_${n}_$$; mkdir -p "$T"
  ( cd "$T" && MMPBSA.py -O -i "$d/nmode.in" -o FINAL_RESULTS_MMPBSA.dat -eo per_frame.csv \
      -cp "$TOP/${SYS}_complex.prmtop" -rp "$TOP/${SYS}_receptor.prmtop" -lp "$TOP/${SYS}_ligand.prmtop" \
      -y "$ct" -yr "$rt" -yl "$lt" >> "$LOG" 2>&1 )
  # keep everything the entropy analysis needs: the per-replica summary, the per-frame energies, and the
  # per-species normal-mode output, which is where the PER-FRAME entropies live (entropy_perframe_SE.py
  # parses _MMPBSA_{complex,receptor,ligand}_nm.out.0). Without these only replica averages survive and
  # snapshot-convergence analysis becomes impossible.
  for f in FINAL_RESULTS_MMPBSA.dat per_frame.csv \
           _MMPBSA_complex_nm.out.0 _MMPBSA_receptor_nm.out.0 _MMPBSA_ligand_nm.out.0; do
    [ -f "$T/$f" ] && cp -p "$T/$f" "$d/"
  done
  rm -rf "$T"
  if grep -aq "Entropy" "$d/FINAL_RESULTS_MMPBSA.dat" 2>/dev/null; then
    echo "$(date +%T) $SYS traj_$n done" >> "$LOG"
  else
    echo "$(date +%T) $SYS traj_$n FAILED or incomplete" >> "$LOG"
  fi
}
export -f nmode_one; export OUT TOP CAMP COND SYS LOG AMBERHOME PATH LD_LIBRARY_PATH
P=${NMODE_PARALLEL:-${SLURM_CPUS_PER_TASK:-20}}
echo "running $P replicas at a time" | tee -a "$LOG"
seq "$FIRST" "$LAST" | xargs -P "$P" -I{} bash -c "nmode_one {}"
echo "finished: $(ls -d $OUT/traj_*/ 2>/dev/null | wc -l) replica directories" | tee -a "$LOG"
