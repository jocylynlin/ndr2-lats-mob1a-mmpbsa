#!/bin/bash
# setup_c2.sh  -- C2: solvent-retaining runs for Reviewer 1 Comment 3.
# Unbound LATS1 and NDR2 kinase, OPC and TIP3P, 5 replicas each, 40 ns, WATER RETAINED.
# Continues from prod.rst so the conformational ensemble is the production one.
# Output every 100 ps (ntwx=50000) -> 400 frames/replica, ~0.17 GB each, ~3.3 GB total.
# OPC uses slots 21-25 to avoid any clash with the C1 extensions (slots 1-20).
CAMPAIGN=/home3/jocelyn/projects/hippo/campaign_100x40ns
STEPS=20000000          # 40 ns at dt=0.002
NTWX=50000              # 100 ps
declare -A PRM=( [lats1_kinase]=lats.prmtop [ndr2_kinase]=stk_chainD.prmtop )
declare -A SLOTS=( [opc-salt]="21 22 23 24 25" [tip3p]="1 2 3 4 5" )
for COND in opc-salt tip3p; do
  for SYS in lats1_kinase ndr2_kinase; do
    for i in ${SLOTS[$COND]}; do
      n=$(printf %03d $i); D=$CAMPAIGN/runs/$COND/$SYS/replica$n
      [ -f "$D/prod.rst" ] || { echo "SKIP $COND $SYS $n: no prod.rst"; continue; }
      [ -f "$D/wat.nc" ]   && { echo "SKIP $COND $SYS $n: wat.nc exists"; continue; }
      sed -e "s/nstlim=[0-9]*/nstlim=$STEPS/" \
          -e "s/ntwx=[0-9]*/ntwx=$NTWX/" \
          -e "s/! 40 ns for production/! 40 ns, solvent retained (C2)/" "$D/prod.in" > "$D/wat.in"
      cat > "$D/submit_wat.sh" <<EOF
#!/bin/bash
#SBATCH --job-name=W${COND:0:1}${SYS:0:1}$n
#SBATCH --output=wat_%j.out
#SBATCH --error=wat_%j.err
#SBATCH --nodes=1 --ntasks=1 --time=24:00:00 --partition=gpu
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/apps/cuda-12.1/lib64:\$LD_LIBRARY_PATH
WORK=\$(pwd); T=/tmp/\$USER/c2_${COND}_${SYS}_${n}_\$SLURM_JOB_ID
mkdir -p \$T && cp \$WORK/${PRM[$SYS]} \$WORK/prod.rst \$WORK/wat.in \$T/ && cd \$T
echo "C2 solvent-retaining run $COND $SYS replica $n on \$(hostname)"
/home3/apps/amber24/bin/pmemd.cuda -O -i wat.in -o wat.out -p ${PRM[$SYS]} \\
    -c prod.rst -r wat.rst -x wat.nc || { cp wat.out \$WORK/ 2>/dev/null; exit 1; }
cp wat.out wat.rst wat.nc \$WORK/
cd \$WORK && rm -rf \$T
echo "$COND $SYS replica $n solvent run done"
EOF
      echo "prepared $COND $SYS $n"
    done
  done
done
