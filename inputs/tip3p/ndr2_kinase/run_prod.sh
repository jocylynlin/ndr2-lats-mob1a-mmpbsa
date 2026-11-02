#!/bin/bash

# Number of replicas
NUM_REPLICAS=20

# MD production length (ns)
PROD_LENGTH=100
# Steps for production (assuming 2 fs timestep)
PROD_STEPS=$((PROD_LENGTH * 500000))

# Function to create production input file
create_prod_input() {
   local replica=$1
   local random_seed=$2
   cat << EOF > replica${replica}/prod.in
Production MD at 310.15K - Replica $replica - ${PROD_LENGTH}ns
&cntrl
  imin=0,        ! No minimization
  ntx=5,         ! Read coords and velocities
  irest=1,       ! Restart
  nstlim=${PROD_STEPS}, ! ${PROD_LENGTH} ns for production
  dt=0.002,      ! Time step (ps)
  ntf=2,         ! SHAKE on H bonds
  ntc=2,         ! SHAKE on H
  temp0=310.15,  ! Target temp
  ntt=3,         ! Langevin thermostat
  gamma_ln=1.0,  ! Collision frequency
  ntpr=5000,     ! Print every 5000 steps
  ntwx=5000,     ! Write coords every 5000 steps
  ntwr=50000,    ! Write restarts every 50000 steps
  cut=10.0,      ! Nonbonded cutoff
  ntb=2,         ! Constant pressure
  ntp=1,         ! Isotropic pressure scaling
  pres0=1.0,     ! Target pressure (bar)
  taup=2.0,      ! Pressure relaxation time
  iwrap=1,       ! Wrap coords in output
  ig=-1, ! Same random seed pattern as before
  ntr=0,         ! No restraints
/
EOF
}

# Function to create submission script for production
create_prod_submission_script() {
   local replica=$1
   
   cat << EOF > replica${replica}/submit_prod_replica${replica}.sh
#!/bin/bash
#SBATCH --job-name=stkt1_R${replica}
#SBATCH --output=PROD_R${replica}_%j.out
#SBATCH --error=PROD_R${replica}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
# No node exclusions (2026-09-20)
# No nodelist - let scheduler decide

export PATH=/home3/apps/gcc/13.2.0/bin:$PATH
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:$LD_LIBRARY_PATH
export CUDA_HOME=/home8/apps/cuda-12.8
export CUDA_TOOLKIT_ROOT_DIR=/home8/apps/cuda-12.8
amber="/home3/apps/amber24/bin/pmemd.cuda"

# Store original working directory
WORKDIR=\$(pwd)

# Create unique temporary directory
TMPDIR="/tmp/\$USER/prod_replica${replica}_\${SLURM_JOB_ID}"
mkdir -p \$TMPDIR

# Copy all necessary files to temporary directory
cp * \$TMPDIR/

# Change to temporary directory
cd \$TMPDIR

echo "Starting production for replica${replica} on node \$(hostname)"
\$amber -O \\
  -i prod.in \\
  -o prod.out \\
  -p stk_chainD.prmtop \\
  -c eq.rst \\
  -r prod.rst \\
  -x prod.nc

# Check if production completed successfully
if [ \$? -ne 0 ]; then
  echo "ERROR: Production failed for replica${replica}"
  # Copy logs back even on failure
  cp prod.out \$WORKDIR/
  exit 1
fi

# Copy results back to working directory
echo "Copying results back to \$WORKDIR..."
cp prod.* \$WORKDIR/

# Clean up temporary directory
cd \$WORKDIR
rm -rf \$TMPDIR

echo "Production for replica${replica} completed successfully"
EOF
   chmod +x replica${replica}/submit_prod_replica${replica}.sh
}

# Create check dependency script
create_dependency_check() {
   cat << EOF > check_before_prod.sh
#!/bin/bash
echo "Checking if replicas are ready for production..."
echo "-----------------------------------------------"

READY=0
NOT_READY=0
READY_REPLICAS=()
NOT_READY_REPLICAS=()

for replica in \$(seq 1 $NUM_REPLICAS); do
  if [ -f "replica\${replica}/eq.rst" ]; then
    echo "Replica \$replica: READY for production"
    READY=\$((READY+1))
    READY_REPLICAS+=(\$replica)
  else
    echo "Replica \$replica: NOT READY (eq.rst not found)"
    NOT_READY=\$((NOT_READY+1))
    NOT_READY_REPLICAS+=(\$replica)
  fi
done

echo ""
echo "Summary: \$READY replicas ready, \$NOT_READY replicas not ready"
echo ""
echo "Ready replicas: \${READY_REPLICAS[@]}"
echo "Not ready replicas: \${NOT_READY_REPLICAS[@]}"
EOF
   chmod +x check_before_prod.sh
}

# Create production files and scripts for all replicas
echo "Creating production files for $NUM_REPLICAS replicas..."
for ((replica=1; replica<=NUM_REPLICAS; replica++)); do
   if [ ! -d "replica${replica}" ]; then
      echo "ERROR: replica${replica} directory not found!"
      continue
   fi
   
   # Create random seed for production (different than previous phases)
   random_seed=$((replica * 100 + 789))
   
   # Create production input files
   create_prod_input $replica $random_seed
   create_prod_submission_script $replica
   
   echo "Configured production for replica${replica}"
done

# Create dependency check script
create_dependency_check

# Submit all jobs
echo "Checking replica status and submitting production jobs..."
submitted_jobs=()
skipped_replicas=()

for ((replica=1; replica<=NUM_REPLICAS; replica++)); do
   if [ -f "replica${replica}/eq.rst" ]; then
      cd replica${replica}
      job_id=$(sbatch --parsable submit_prod_replica${replica}.sh)
      submitted_jobs+=("$job_id")
      cd ..
      echo "Submitted production for replica${replica} (Job ID: $job_id)"
   else
      echo "Skipping replica${replica} - not ready for production (eq.rst missing)"
      skipped_replicas+=("$replica")
   fi
done

# Save job IDs to a file for reference
echo "${submitted_jobs[@]}" > production_job_ids.txt
echo "All production job IDs saved to production_job_ids.txt"

if [ ${#skipped_replicas[@]} -gt 0 ]; then
   echo "Skipped replicas: ${skipped_replicas[@]}"
   echo "Run './check_before_prod.sh' for detailed status"
fi

# Create a summary script to check production status
cat << EOF > check_production_status.sh
#!/bin/bash
echo "Production Status Summary:"
echo "-------------------------"
for replica in \$(seq 1 $NUM_REPLICAS); do
  # Check job status using sacct
  job_id=\$(grep -o "\b[0-9]\+" replica\${replica}/PROD_R\${replica}_*.out 2>/dev/null | head -1)
  if [ -z "\$job_id" ]; then
    job_status="UNKNOWN"
  else
    job_status=\$(sacct -j \$job_id --format=State --noheader | head -1 | tr -d ' ')
  fi
  
  # Check file status
  if [ -f "replica\${replica}/prod.rst" ]; then
    file_status="COMPLETE"
  elif [ -f "replica\${replica}/prod.out" ]; then
    file_status="IN PROGRESS"
  elif [ -f "replica\${replica}/eq.rst" ]; then
    file_status="READY BUT NOT STARTED"
  else
    file_status="NOT READY"
  fi
  
  echo "Replica \$replica: Job \$job_status, Files \$file_status"
done
EOF
chmod +x check_production_status.sh

echo "Created check_production_status.sh to monitor production progress"
echo "To monitor all jobs, run: squeue -u \$USER"
