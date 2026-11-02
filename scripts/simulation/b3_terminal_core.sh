#!/bin/bash
#SBATCH --job-name=b3term
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/b3_%j.out
#SBATCH --error=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/b3_%j.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=04:00:00 --partition=cpu-short
# B3: is the negative covalent term carried by the frayed termini or by the 49-residue core?
# Bonded energy of the FULL fragment (1-59) and of the CORE only (6-54), bound and unbound,
# reference condition. Bonds spanning the core boundary are absent from the core calculation,
# but that offset is identical in bound and unbound, so it cancels in the reorganization.
export LD_LIBRARY_PATH=/home3/jocelyn/projects/hippo/campaign_100x40ns/lib:$LD_LIBRARY_PATH
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
O=$C/analysis/b3_terminal; mkdir -p "$O"
for SY in lats1 ndr2; do
  for ST in unbound bound; do
    if [ "$ST" = unbound ]; then
      PRM=$C/topologies/opc-salt/${SY}_kinase.prmtop; STRIP=""
      TRJD=$C/trajectories/opc-salt/${SY}_kinase; PFX=${SY}_kinase
    else
      PRM=$C/topologies/opc-salt/${SY}_complex.prmtop; STRIP="strip :1-181"
      TRJD=$C/trajectories/opc-salt/${SY}_complex; PFX=${SY}_complex
    fi
    for i in $(seq 1 20); do
      n=$(printf %03d $i)
      [ -s $O/${SY}_${ST}_$n.dat ] && continue
      $CP -p "$PRM" <<CPT
trajin $TRJD/${PFX}_${n}.nc 1 4000 10
$STRIP
energy FULL :1-59 out $O/${SY}_${ST}_$n.dat bond angle dihedral
energy CORE :6-54 out $O/${SY}_${ST}_core_$n.dat bond angle dihedral
go
CPT
    done
  done
done
echo B3DONE
