#!/bin/bash
# generate_replicas.sh <condition> <tag> [first_slot] [last_slot]
# Creates runs/<condition>/<tag>/replicaNNN with inputs and submit scripts. Submits nothing.
. "$(dirname "$0")/lib.sh"
COND=$1; TAG=$2
REUSED=$(reused_count "$COND" "$TAG"); FIRST=${3:-$((REUSED+1))}; LAST=${4:-$N_TARGET}
IN=$(inputs_dir "$COND" "$TAG"); RUN=$(runs_dir "$COND" "$TAG"); PRM=$(basename "$(prmtop_of "$COND" "$TAG")")
[ -f "$IN/$PRM" ] && [ -f "$IN/min2.rst" ] || { echo "missing inputs in $IN"; exit 1; }
mkdir -p "$RUN"
for ((s=FIRST; s<=LAST; s++)); do
  d=$(printf "%s/replica%03d" "$RUN" "$s"); mkdir -p "$d"
  cp -n "$IN/$PRM" "$IN/min2.rst" "$d/"
  # explicit, unique, reproducible seeds: index = condition, tag and slot -> one seed per stage
  ci=0; for c in $CONDITIONS; do [ "$c" = "$COND" ] && break; ci=$((ci+1)); done
  ti=0; for t in $TAGS;       do [ "$t" = "$TAG"  ] && break; ti=$((ti+1)); done
  idx=$(( ci*600 + ti*100 + s ))
  SEED_H=$(( SEED_BASE + idx*3 )); SEED_E=$(( SEED_H + 1 )); SEED_P=$(( SEED_H + 2 ))
  sed -e "s/REPLICA_ID/$s/" -e "s/ig=-1/ig=$SEED_H/" "$IN/heat.in" > "$d/heat.in"
  sed -e "s/REPLICA_ID/$s/" -e "s/ig=-1/ig=$SEED_E/" "$IN/eq.in"   > "$d/eq.in"
  echo "$COND,$TAG,$s,$SEED_H,$SEED_E,$SEED_P" >> "$CAMPAIGN/logs/seeds_assigned.csv"
  sed -e "s/REPLICA_ID/$s/" -e "s/ig=-1/ig=$SEED_P/" "$IN/prod.in" > "$d/prod.in"
  for ph in eq prod; do
    if [ "$ph" = eq ]; then jn="E"; steps="heat+eq"; else jn="P"; steps="prod"; fi
    cat > "$d/submit_${ph}.sh" <<SB
#!/bin/bash
#SBATCH --job-name=${jn}$(printf %03d $s)
#SBATCH --output=${ph}_%j.out
#SBATCH --error=${ph}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=$PARTITION
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/apps/cuda-12.1/lib64:\$LD_LIBRARY_PATH
WORK=\$(pwd); T=/tmp/\$USER/${COND}_${TAG}_$(printf %03d $s)_\$SLURM_JOB_ID
mkdir -p \$T && cp \$WORK/* \$T/ 2>/dev/null; cd \$T
if [ "$ph" = eq ]; then
  echo "heating replica $s on \$(hostname)"
  $AMBER -O -i heat.in -o heat.out -p $PRM -c min2.rst -r heat.rst -x heat.nc || { cp heat.out \$WORK/; exit 1; }
  echo "equilibrating replica $s"
  $AMBER -O -i eq.in -o eq.out -p $PRM -c heat.rst -r eq.rst -x eq.nc || { cp eq.out \$WORK/; exit 1; }
  cp heat.out heat.rst eq.out eq.rst \$WORK/
else
  cp \$WORK/eq.rst \$T/ 2>/dev/null
  echo "production replica $s on \$(hostname)"
  $AMBER -O -i prod.in -o prod.out -p $PRM -c eq.rst -r prod.rst -x prod.nc || { cp prod.out \$WORK/; exit 1; }
  cp prod.out prod.rst prod.nc \$WORK/
fi
cd \$WORK && rm -rf \$T
echo "replica $s $steps done"
SB
    chmod +x "$d/submit_${ph}.sh"
  done
done
echo "generated $COND/$TAG slots $FIRST-$LAST in $RUN"
