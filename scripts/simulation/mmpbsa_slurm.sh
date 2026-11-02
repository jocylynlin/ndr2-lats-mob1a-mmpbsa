#!/bin/bash
#SBATCH --job-name=mmpbsa
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/mmpbsa_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-long
# mmpbsa_slurm.sh <condition> <st|mt> <first> <last> <systems...>
# Verifies the AMBER environment on this node first (compute nodes lack libreadline.so.6; the campaign ships
# a copy), then runs MM-PBSA with 24 replicas in parallel, each in node-local scratch.
cd /home3/jocelyn/projects/hippo/campaign_100x40ns || exit 1
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_100x40ns/lib:$AMBERHOME/lib:$LD_LIBRARY_PATH
COND=$1; MODE=$2; FIRST=$3; LAST=$4; shift 4; SYSTEMS="$*"
echo "node $(hostname), $(nproc) cores | $MODE $COND replicas $FIRST-$LAST | systems: $SYSTEMS"
T=$(mktemp -d); ( cd $T && cpptraj -i /dev/stdin > cp.log 2>&1 <<CP
parm /home3/jocelyn/projects/hippo/campaign_100x40ns/topologies/$COND/ndr2_complex.prmtop
trajin /home3/jocelyn/projects/hippo/campaign_100x40ns/trajectories/$COND/ndr2_complex/ndr2_complex_001.nc 1 2
surf t out t.dat
go
CP
)
[ -s $T/t.dat ] || { echo "ENVIRONMENT CHECK FAILED on $(hostname)"; tail -3 $T/cp.log; rm -rf $T; exit 1; }
rm -rf $T; echo "environment check passed"
for SYS in $SYSTEMS; do PARALLEL=24 bash scripts/run_mmpbsa.sh $COND $MODE $FIRST $LAST $SYS; done
echo "finished: $MODE $COND $SYSTEMS"
