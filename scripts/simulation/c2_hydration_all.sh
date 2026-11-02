#!/bin/bash
#SBATCH --job-name=c2hydall
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/c2_hyd_%j.out
#SBATCH --error=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/c2_hyd_%j.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=20 --time=4:00:00 --partition=cpu-short
# Hydration analysis of every completed C2 solvent-retaining run, 20 cpptraj jobs concurrent.
# Idempotent: a replica whose sasa.dat already exists is skipped, so re-running after the last
# run lands only processes the new one.
set -u
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
S=$C/scripts; NPROC=${SLURM_CPUS_PER_TASK:-20}
TASKS=$C/logs/c2_hyd_tasks.txt; : > $TASKS
for COND in opc-salt tip3p; do
  case $COND in opc-salt) SL="021 022 023 024 025";; tip3p) SL="001 002 003 004 005";; esac
  for SYS in lats1_kinase ndr2_kinase; do for n in $SL; do
    [ -s "$C/runs/$COND/$SYS/replica$n/wat.nc" ] || { echo "PENDING $COND $SYS $n (no wat.nc)"; continue; }
    [ -s "$C/analysis/c2_hydration/$COND/$SYS/replica$n/sasa.dat" ] && { echo "SKIP $COND $SYS $n"; continue; }
    echo "$COND $SYS $n" >> $TASKS
  done; done
done
echo "$(date '+%F %T')  tasks to run: $(wc -l < $TASKS)"
[ -s "$TASKS" ] && xargs -a $TASKS -P "$NPROC" -n 3 $S/c2_hydration_one.sh
echo "$(date '+%F %T')  C2HYDDONE  complete replicas: $(find $C/analysis/c2_hydration -name sasa.dat -size +1k | wc -l)/20"
