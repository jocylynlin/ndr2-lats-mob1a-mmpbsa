#!/bin/bash
#SBATCH --job-name=min_lats1_kinase
#SBATCH --output=min_%j.out
#SBATCH --error=min_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --partition=gpu
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:/home8/apps/gcc/11.3.0/lib64:/home3/apps/cuda-12.1/lib64:$LD_LIBRARY_PATH
cd $SLURM_SUBMIT_DIR
echo "minimising tip3p/lats1_kinase on $(hostname)"
/home3/apps/amber24/bin/pmemd.cuda -O -i min1.in -o min1.out -p lats.prmtop -c lats.inpcrd -r min1.rst || exit 1
/home3/apps/amber24/bin/pmemd.cuda -O -i min2.in -o min2.out -p lats.prmtop -c min1.rst      -r min2.rst || exit 1
echo "done"
