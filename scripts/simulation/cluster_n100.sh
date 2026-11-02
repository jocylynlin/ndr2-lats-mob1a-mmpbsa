#!/bin/bash
#SBATCH --job-name=clust100
#SBATCH --output=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/cluster_%j.out
#SBATCH --error=/home3/jocelyn/projects/hippo/campaign_100x40ns/logs/cluster_%j.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=6:00:00 --partition=cpu-short
# K-means clustering of the kinase CA coordinates, bound and unbound, reference condition.
# Same settings as the submitted analysis: 5 clusters, kmeans, maxit 500.
# All 100 replicas per ensemble, strided to every 10th frame (400 -> 40 per replica,
# 4,000 frames per ensemble) to keep the pairwise work tractable.
set -u
C=/home3/jocelyn/projects/hippo/campaign_100x40ns
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:$C/lib
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
O=$C/analysis/cluster_n100; mkdir -p $O

for SY in lats1 ndr2; do
  # ---- unbound kinase: its own trajectories ----
  {
    echo "parm $C/topologies/opc-salt/${SY}_ligand.prmtop"
    for i in $(seq 1 100); do
      n=$(printf %03d $i)
      f=$C/trajectories/opc-salt/${SY}_kinase/${SY}_kinase_${n}.nc
      [ -s "$f" ] && echo "trajin $f 1 4000 10"
    done
    echo "rms first :1-59@CA"
    echo "cluster C0 kmeans clusters 5 maxit 500 rms :1-59@CA \\"
    echo "  out $O/${SY}_unbound_cnum.dat summary $O/${SY}_unbound_summary.dat \\"
    echo "  repout $O/${SY}_unbound_rep repfmt pdb"
    echo "run"
  } > $O/${SY}_unbound.in
  env -u PERL5LIB $CP -i $O/${SY}_unbound.in > $O/${SY}_unbound.log 2>&1
  echo "  $SY unbound done"

  # ---- bound kinase: taken from the complex, MOB1A stripped ----
  {
    echo "parm $C/topologies/opc-salt/${SY}_complex.prmtop"
    for i in $(seq 1 100); do
      n=$(printf %03d $i)
      f=$C/trajectories/opc-salt/${SY}_complex/${SY}_complex_${n}.nc
      [ -s "$f" ] && echo "trajin $f 1 4000 10"
    done
    echo "strip :1-181"
    echo "rms first :1-59@CA"
    echo "cluster C0 kmeans clusters 5 maxit 500 rms :1-59@CA \\"
    echo "  out $O/${SY}_bound_cnum.dat summary $O/${SY}_bound_summary.dat \\"
    echo "  repout $O/${SY}_bound_rep repfmt pdb"
    echo "run"
  } > $O/${SY}_bound.in
  env -u PERL5LIB $CP -i $O/${SY}_bound.in > $O/${SY}_bound.log 2>&1
  echo "  $SY bound done"
done
echo CLUSTERDONE
