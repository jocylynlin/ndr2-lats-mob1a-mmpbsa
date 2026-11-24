#!/bin/bash
#SBATCH --job-name=pbsens
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_saltmatched/logs/pbsens_%j.out
#SBATCH --error=/home3/jocelyn/projects/hippo/campaign_saltmatched/logs/pbsens_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --partition=cpu-long
#SBATCH --exclude=cpu-2-0,cpu-2-1
# pb_sens_slurm.sh <condition> <istrng> [nrep] [interval]
cd /home3/jocelyn/projects/hippo/campaign_saltmatched
# The compute nodes lack libreadline.so.6 / libncurses.so.5 / libtinfo.so.5, which cpptraj
# needs; only the login node has them. The campaign stages copies in <campaign>/lib and puts
# that on LD_LIBRARY_PATH (see campaign_100x40ns/scripts/mmpbsa_slurm.sh). Same here.
export AMBERHOME=/home3/jocelyn/apps/amber25
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_saltmatched/lib:$AMBERHOME/lib:${LD_LIBRARY_PATH:-}
echo "host $(hostname)  cond=$1 istrng=$2 nrep=${3:-10} interval=${4:-40} mode=${5:-mt}"
PARALLEL=12 bash scripts/pb_istrng_sensitivity.sh "$1" "$2" "${3:-10}" "${4:-40}" "${5:-mt}"
