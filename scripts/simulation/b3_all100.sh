#!/bin/bash
#SBATCH --job-name=b3_100
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/b3_100_%j.out
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=04:00:00 --partition=cpu-short
# Extend the terminal-vs-core bonded decomposition from 20 to all 100 replicas.
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
for SY in lats1 ndr2; do for ST in unbound bound; do for i in $(seq 1 100); do
  printf '%s %s %03d\n' "$SY" "$ST" "$i"
done; done; done > $C/logs/b3_tasks.txt
echo "tasks: $(wc -l < $C/logs/b3_tasks.txt)"
xargs -a $C/logs/b3_tasks.txt -P 32 -n 3 $C/scripts/b3_worker.sh
echo "B3_100_DONE  dat files: $(ls $C/analysis/b3_terminal_all/*.dat 2>/dev/null | wc -l)"
