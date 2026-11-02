#!/bin/bash
# minimize.sh <condition> <tag> [more cond/tag pairs via repeated calls]
# Two-stage minimisation (steepest descent then conjugate gradient, 30,000 cycles each) of a rebuilt system,
# producing the min2.rst that heating starts from. Inputs are taken from the slot's pre_zaff_archive so the
# protocol matches the original build exactly; the engine is pmemd.cuda rather than the original CPU
# pmemd.MPI (recorded in PROVENANCE.txt).
. "$(dirname "$0")/lib.sh"
COND=$1; TAG=$2
IN=$(inputs_dir "$COND" "$TAG"); PRM=$(basename "$(prmtop_of "$COND" "$TAG")"); STEM=${PRM%.prmtop}
[ -f "$IN/$STEM.inpcrd" ] || { echo "no $STEM.inpcrd in $IN"; exit 1; }
cp -n "$IN/pre_zaff_archive/min1.in" "$IN/pre_zaff_archive/min2.in" "$IN/" 2>/dev/null
cat > "$IN/submit_min.sh" <<SB
#!/bin/bash
#SBATCH --job-name=min_${TAG}
#SBATCH --output=min_%j.out
#SBATCH --error=min_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --partition=$PARTITION
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/apps/cuda-12.1/lib64:\$LD_LIBRARY_PATH
cd \$SLURM_SUBMIT_DIR
echo "minimising $COND/$TAG on \$(hostname)"
$AMBER -O -i min1.in -o min1.out -p $PRM -c $STEM.inpcrd -r min1.rst || exit 1
$AMBER -O -i min2.in -o min2.out -p $PRM -c min1.rst      -r min2.rst || exit 1
echo "done"
SB
chmod +x "$IN/submit_min.sh"
id=$(cd "$IN" && sbatch --parsable submit_min.sh) && echo "$COND/$TAG submitted as job $id" || echo "$COND/$TAG SUBMIT FAILED"
