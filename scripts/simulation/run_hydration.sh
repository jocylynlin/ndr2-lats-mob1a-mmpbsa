#!/bin/bash
#SBATCH --job-name=hydmatch
#SBATCH --output=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_hydration_matched/logs/h_%A_%a.out
#SBATCH --error=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_hydration_matched/logs/h_%A_%a.err
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=04:00:00 --partition=cpu-short
#SBATCH --array=0-79%40
#
# Hydration and ion-distribution analysis on the FINAL matched-salt solvent-retaining set.
# 20 replicas per (condition, kinase) -- four times the 5 used for the submitted figures.
# prod.nc keeps solvent at 4000 frames / 10 ps; strided by 10 to 400 frames, matching the
# 100 ps sampling of the original analysis so the two are directly comparable.
#
# Charged-sidechain masks are identical to the submitted analysis; they were re-derived
# from the matched topologies' RESIDUE_LABEL and agree residue for residue.
set -eu
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/jocelyn/projects/hippo/campaign_100x40ns/lib
CP=/home3/jocelyn/apps/amber25/bin/cpptraj
NEW=/home3/jocelyn/projects/hippo/campaign_saltmatched
OUT=/home3/jocelyn/projects/hippo/reviewer_comments/analysis_hydration_matched/data

i=$SLURM_ARRAY_TASK_ID
rep=$(( i % 20 + 1 )); k=$(( i / 20 ))
case $k in
  0) COND=opc-salt-matched; SYS=lats1_kinase ;;
  1) COND=opc-salt-matched; SYS=ndr2_kinase  ;;
  2) COND=tip3p-matched;    SYS=lats1_kinase ;;
  3) COND=tip3p-matched;    SYS=ndr2_kinase  ;;
esac
N=$(printf %03d $rep)
case $SYS in
  lats1_kinase) PRM=lats.prmtop;       CAT="1,13,18,21,22,23,31,43,44,49,55,57,58"; ANI="5,9,26,28,37,40,50"   ;;
  ndr2_kinase)  PRM=stk_chainD.prmtop; CAT="1,18,21,23,24,39,40,42,43,48,49,55,57,58"; ANI="5,17,26,30,32,36,50,52" ;;
esac
D=$NEW/runs/$COND/$SYS/replica$N
O=$OUT/$COND/$SYS/replica$N
[ -f "$D/prod.nc" ] || { echo "NO prod.nc: $COND $SYS $N"; exit 1; }
mkdir -p "$O"
CATM=":$CAT@NZ,NE,NH1,NH2,ND1,NE2"
ANIM=":$ANI@OD1,OD2,OE1,OE2"
$CP -p "$D/$PRM" <<CPP > "$O/cpptraj.log" 2>&1
trajin $D/prod.nc 1 4000 10
autoimage
hbond PW :1-59 solventdonor :WAT solventacceptor :WAT@O out $O/hbond.dat avgout $O/hbond_avg.dat
watershell :1-59 out $O/ws_protein.dat lower 3.4 upper 5.0
watershell $CATM out $O/ws_cationic.dat lower 3.4 upper 5.0
watershell $ANIM out $O/ws_anionic.dat lower 3.4 upper 5.0
radial $O/rdf_cat_watO.dat 0.1 10.0 :WAT@O $CATM
radial $O/rdf_ani_watO.dat 0.1 10.0 :WAT@O $ANIM
radial $O/rdf_cat_Cl.dat 0.1 12.0 :Cl- $CATM
radial $O/rdf_ani_Na.dat 0.1 12.0 :Na+ $ANIM
surf :1-59 out $O/sasa.dat
run
CPP
for f in rdf_cat_watO rdf_ani_watO rdf_cat_Cl rdf_ani_Na ws_protein; do
  [ -s "$O/$f.dat" ] || { echo "EMPTY $f : $COND $SYS $N"; exit 1; }
done
echo "OK $COND $SYS $N"
