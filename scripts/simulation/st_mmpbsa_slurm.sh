#!/bin/bash
#SBATCH --job-name=st_mmpbsa
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/st_mmpbsa_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --partition=cpu-long
# st_mmpbsa_slurm.sh <condition> <systems...>
# Runs ST MM-PBSA on a compute node, 24 replicas at a time. Verifies the AMBER environment works on this
# node before starting: the earlier SASA attempt failed on cpu-long because cpptraj produced no output there.
cd /home3/jocelyn/projects/hippo/campaign_100x40ns || exit 1
export AMBERHOME=/home3/jocelyn/apps/amber25; export PATH=$AMBERHOME/bin:$PATH
# compute nodes lack libreadline.so.6 (the login node has it); ship a copy with the campaign so the SAME
# AmberTools build runs in both places and results stay comparable
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_100x40ns/lib:$AMBERHOME/lib:$LD_LIBRARY_PATH
COND=${1:-opc-salt}; shift; SYSTEMS=${*:-ndr2}
echo "node $(hostname), $(nproc) cores, systems: $SYSTEMS"
T=$(mktemp -d); cd $T
cpptraj -i /dev/stdin > cp.log 2>&1 <<CP
parm /home3/jocelyn/projects/hippo/campaign_100x40ns/topologies/$COND/ndr2_complex.prmtop
trajin /home3/jocelyn/projects/hippo/campaign_100x40ns/trajectories/$COND/ndr2_complex/ndr2_complex_001.nc 1 2
surf test out t.dat
go
CP
if [ ! -s t.dat ]; then echo "ENVIRONMENT CHECK FAILED on $(hostname): cpptraj produced no output"; cat cp.log | tail -5; exit 1; fi
echo "environment check passed (cpptraj works on $(hostname))"
cd /home3/jocelyn/projects/hippo/campaign_100x40ns; rm -rf $T
for SYS in $SYSTEMS; do
  PARALLEL=24 bash scripts/run_mmpbsa.sh $COND st 1 100 $SYS
done
echo "ST MM-PBSA finished for: $SYSTEMS"
