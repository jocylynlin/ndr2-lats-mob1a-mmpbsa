#!/bin/bash
# c2_hydration_one.sh COND SYS REPLICA
# Hydration analysis of one C2 solvent-retaining run (wat.nc, 400 frames at 100 ps).
# Measurable at this output interval: time-averaged protein-water H-bond counts, first/second
# shell water occupancies (whole protein and charged sidechains), RDFs of water O and of ions
# around charged sidechains. NOT measurable: residence times / H-bond lifetimes (100 ps spacing).
set -eu
COND=$1; SYS=$2; N=$3
CAMPAIGN=/home3/jocelyn/projects/hippo/campaign_100x40ns
# gcc 13.2 MUST precede the node default gcc 11.3: amber24 libemil.so needs GLIBCXX_3.4.32.
# $CAMPAIGN/lib supplies libreadline.so.6 / libncurses.so.5 / libtinfo.so.5 that cpptraj links.
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:$CAMPAIGN/lib
CPPTRAJ=/home3/apps/amber24/bin/cpptraj
case $SYS in lats1_kinase) PRM=lats.prmtop;  CAT="1,13,18,21,22,23,31,43,44,49,55,57,58"; ANI="5,9,26,28,37,40,50";;
             ndr2_kinase)  PRM=stk_chainD.prmtop; CAT="1,18,21,23,24,39,40,42,43,48,49,55,57,58"; ANI="5,17,26,30,32,36,50,52";; esac
D=$CAMPAIGN/runs/$COND/$SYS/replica$N
O=$CAMPAIGN/analysis/c2_hydration/$COND/$SYS/replica$N
[ -f "$D/wat.nc" ] || { echo "NO wat.nc: $COND $SYS $N"; exit 1; }
mkdir -p "$O"
CATM=":$CAT@NZ,NE,NH1,NH2,ND1,NE2"
ANIM=":$ANI@OD1,OD2,OE1,OE2"
$CPPTRAJ -p "$D/$PRM" <<CPP > "$O/cpptraj.log" 2>&1
trajin $D/wat.nc
autoimage
# --- protein-water hydrogen bonds (solute-solvent) ---
hbond PW :1-59 solventdonor :WAT solventacceptor :WAT@O out $O/hbond.dat avgout $O/hbond_avg.dat
# --- hydration shell occupancies ---
watershell :1-59 out $O/ws_protein.dat lower 3.4 upper 5.0
watershell $CATM out $O/ws_cationic.dat lower 3.4 upper 5.0
watershell $ANIM out $O/ws_anionic.dat lower 3.4 upper 5.0
# --- radial distribution functions, water oxygen around charged sidechains ---
radial $O/rdf_cat_watO.dat 0.1 10.0 :WAT@O $CATM
radial $O/rdf_ani_watO.dat 0.1 10.0 :WAT@O $ANIM
# --- ion distributions around charged sidechains ---
radial $O/rdf_cat_Cl.dat 0.1 12.0 :Cl- $CATM
radial $O/rdf_ani_Na.dat 0.1 12.0 :Na+ $ANIM
# --- control: solvent-accessible surface, the ENPOLAR driver ---
surf :1-59 out $O/sasa.dat
run
CPP
nf=$(awk 'END{print NR}' "$O/ws_protein.dat" 2>/dev/null || echo 0)
echo "done $COND $SYS $N  (frames in ws_protein.dat incl header: $nf)"
