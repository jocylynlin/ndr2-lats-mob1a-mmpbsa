#!/bin/bash
#SBATCH --job-name=physcv
#SBATCH --output=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_physcv_matched/logs/p_%A_%a.out
#SBATCH --error=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_physcv_matched/logs/p_%A_%a.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=4G --time=02:00:00 --partition=cpu-short
#SBATCH --array=0-399%60
#
# Physical collective variables for the UNBOUND kinase, recomputed on the matched-salt
# trajectories. Definitions are byte-identical to physcv_one.sh of the submitted analysis:
#   angle : dot-product angle between the helix A (18->31 CA) and helix B (38->57 CA) axes
#   rg    : radius of gyration of residues 1-59, HEAVY ATOMS (:1-59&!@H=)
#   rmsd  : CA RMSD to the same per-system bound reference, campaign_100x40ns/refs/ref_SYS.pdb
# Same frame selection (1-4000 stride 10 -> 400 frames) and the same reference files, so the
# recomputed unbound CVs sit on exactly the grid the retained bound CVs were built on.
#
# BOUND panels are NOT recomputed: the complex topologies and trajectories are byte-identical
# to the originals (md5 verified), so campaign_100x40ns/analysis/physcv/*/bound is current.
set -eu
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/jocelyn/projects/hippo/campaign_100x40ns/lib
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
NEW=/home3/jocelyn/projects/hippo/campaign_saltmatched
REFD=/home3/jocelyn/projects/hippo/campaign_100x40ns/refs
OUT=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_physcv_matched/data

i=$SLURM_ARRAY_TASK_ID
rep=$(( i % 100 + 1 )); k=$(( i / 100 ))
case $k in
  0) COND=opc-salt-matched; SYS=lats1 ;;  1) COND=opc-salt-matched; SYS=ndr2 ;;
  2) COND=tip3p-matched;    SYS=lats1 ;;  3) COND=tip3p-matched;    SYS=ndr2 ;;
esac
S=$(printf %03d $rep)
O=$OUT/$COND/unbound; mkdir -p "$O"
CSV=$O/${SYS}_${S}.csv
[ -s "$CSV" ] && { echo "exists $CSV"; exit 0; }

TRJ=$NEW/trajectories/$COND/${SYS}_kinase/${SYS}_kinase_${S}.nc
PRM=$NEW/topologies/$COND/${SYS}_kinase.prmtop
REF=$REFD/ref_${SYS}.pdb
[ -f "$TRJ" ] || { echo "NO TRAJ $TRJ"; exit 1; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$CP -p "$PRM" >/dev/null 2>&1 <<CPT
parm $REF name REF
reference $REF parm REF name BOUND
trajin $TRJ 1 4000 10
vector v1 :18@CA :31@CA
vector v2 :38@CA :57@CA
vectormath vec1 v1 vec2 v2 dotangle out $T/ang.dat
radgyr :1-59&!@H= out $T/rg.dat
rms toBound ref BOUND :1-59@CA out $T/rms.dat
go
CPT
paste <(awk 'NR>1{print $2}' $T/ang.dat) <(awk 'NR>1{print $2}' $T/rg.dat) \
      <(awk 'NR>1{print $2}' $T/rms.dat) | awk 'BEGIN{print "frame,angle,rg,rmsd"}
      {printf "%d,%.3f,%.4f,%.4f\n",NR,$1,$2,$3}' > "$CSV"
n=$(awk 'END{print NR}' "$CSV")
[ "$n" -eq 401 ] || { echo "BAD ROWS $n in $CSV"; exit 1; }
echo "OK $COND $SYS $S rows=$n"
