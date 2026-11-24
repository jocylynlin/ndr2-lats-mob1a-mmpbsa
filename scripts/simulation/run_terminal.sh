#!/bin/bash
#SBATCH --job-name=termcore
#SBATCH --output=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_terminal_core/logs/t_%A_%a.out
#SBATCH --error=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_terminal_core/logs/t_%A_%a.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --time=02:00:00 --partition=cpu-short
#SBATCH --array=0-399%60
#
# Terminal-vs-core covalent decomposition on the FINAL dataset.
#
#   bound   : complex trajectory, MOB1A (:1-181) stripped -> kinase renumbered 1-59.
#             Complex topologies and trajectories are byte-identical to the originals
#             (md5 verified), because the complex boxes were kept as built.
#   unbound : matched-salt kinase trajectory and its own topology (NOT the superseded
#             fixed-80-pair kinase boxes).
#
#   FULL = :1-59 (whole fragment)   CORE = :6-54 (49 residues)
#   termini = FULL - CORE = residues 1-5 and 55-59, ten residues in all.
#   Bonds crossing the core boundary are absent from CORE, but that offset is identical
#   in bound and unbound, so it cancels in the bound-minus-unbound reorganization.
#
# cpptraj `energy` with no keyword list emits bond, angle, dih, vdw14, elec14, vdw, elec.
# The manuscript's covalent term is bond+angle+dih+vdw14+elec14; vdw/elec are written
# but not used, so the raw files also support the nonbonded check if ever wanted.
set -u
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_100x40ns/lib:$LD_LIBRARY_PATH
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
NEW=/home3/jocelyn/projects/hippo/campaign_saltmatched
O=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_terminal_core/data

i=$SLURM_ARRAY_TASK_ID
rep=$(( i % 100 + 1 )); k=$(( i / 100 ))
case $k in
  0) SY=lats1; ST=bound   ;; 1) SY=lats1; ST=unbound ;;
  2) SY=ndr2;  ST=bound   ;; 3) SY=ndr2;  ST=unbound ;;
esac
n=$(printf %03d $rep)

if [ "$ST" = unbound ]; then
  PRM=$NEW/topologies/opc-salt-matched/${SY}_kinase.prmtop
  TRJ=$NEW/trajectories/opc-salt-matched/${SY}_kinase/${SY}_kinase_${n}.nc
  STRIP=""
else
  PRM=$NEW/topologies/opc-salt-matched/${SY}_complex.prmtop
  TRJ=$NEW/trajectories/opc-salt-matched/${SY}_complex/${SY}_complex_${n}.nc
  STRIP="strip :1-181"
fi
[ -s "$TRJ" ] || { echo "MISSING TRAJ $TRJ"; exit 1; }

$CP -p "$PRM" <<CPT
trajin $TRJ 1 4000 10
$STRIP
energy FULL :1-59 out $O/${SY}_${ST}_${n}_full.dat
energy CORE :6-54 out $O/${SY}_${ST}_${n}_core.dat
go
CPT
rc=$?
[ -s "$O/${SY}_${ST}_${n}_full.dat" ] && [ -s "$O/${SY}_${ST}_${n}_core.dat" ] || { echo "EMPTY OUTPUT $SY $ST $n"; exit 1; }
echo "OK $SY $ST $n rc=$rc"
