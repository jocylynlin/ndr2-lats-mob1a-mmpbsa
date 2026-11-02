#!/bin/bash
# extend_c1.sh [NREP]  -- C1: extend unbound kinase replicas from 40 ns to 100 ns (OPC reference).
# Restarts from prod.rst; writes prod2.{in,out,rst,nc}. Originals are never touched.
# Complexes and unbound MOB1A already have 100 ns for these slots, so only the kinases are extended.
CAMPAIGN=/home3/jocelyn/projects/hippo/campaign_100x40ns
NREP=${1:-20}
STEPS=30000000          # 60 ns at dt=0.002
declare -A PRM=( [lats1_kinase]=lats.prmtop [ndr2_kinase]=stk_chainD.prmtop )
for SYS in lats1_kinase ndr2_kinase; do
  for i in $(seq 1 $NREP); do
    n=$(printf %03d $i); D=$CAMPAIGN/runs/opc-salt/$SYS/replica$n
    [ -f "$D/prod.rst" ] || { echo "SKIP $SYS $n: no prod.rst"; continue; }
    [ -f "$D/prod2.nc" ] && { echo "SKIP $SYS $n: prod2.nc exists"; continue; }
    sed -e "s/nstlim=20000000/nstlim=$STEPS/" \
        -e "s/! 40 ns for production/! 60 ns extension (40 -> 100 ns)/" \
        -e "s/- 40ns/- extension 40-100ns/" "$D/prod.in" > "$D/prod2.in"
    cat > "$D/submit_prod2.sh" <<EOF
#!/bin/bash
#SBATCH --job-name=X${SYS:0:1}$n
#SBATCH --output=prod2_%j.out
#SBATCH --error=prod2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/apps/cuda-12.1/lib64:\$LD_LIBRARY_PATH
WORK=\$(pwd); T=/tmp/\$USER/c1_${SYS}_${n}_\$SLURM_JOB_ID
mkdir -p \$T && cp \$WORK/${PRM[$SYS]} \$WORK/prod.rst \$WORK/prod2.in \$T/ && cd \$T
echo "C1 extension $SYS replica $n on \$(hostname)"
/home3/apps/amber24/bin/pmemd.cuda -O -i prod2.in -o prod2.out -p ${PRM[$SYS]} \\
    -c prod.rst -r prod2.rst -x prod2.nc || { cp prod2.out \$WORK/ 2>/dev/null; exit 1; }
cp prod2.out prod2.rst prod2.nc \$WORK/
cd \$WORK && rm -rf \$T
echo "$SYS replica $n extension done"
EOF
    ( cd "$D" && sbatch submit_prod2.sh >/dev/null 2>&1 && echo "submitted $SYS $n" )
  done
done
