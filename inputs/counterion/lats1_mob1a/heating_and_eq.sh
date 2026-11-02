#!/bin/bash

# Number of replicas
NUM_REPLICAS=20

# Function to create heating input file
create_heating_input() {
   local replica=$1
   local random_seed=$2
   cat << EOF > replica${replica}/heat.in
Heating from 0K to 310.15K - Replica $replica
&cntrl
  imin=0,        ! No minimization
  ntx=1,         ! Read coordinates only
  irest=0,       ! No restart
  nstlim=500000, ! 1ns for heating
  dt=0.002,      ! Time step (ps)
  ntf=2,         ! SHAKE on bonds with H
  ntc=2,         ! SHAKE on H
  tempi=0.0,     ! Initial temp
  temp0=310.15,  ! Final temp
  ntt=3,         ! Langevin thermostat
  gamma_ln=2.0,  ! Collision frequency
  ntpr=1000,     ! Print every 1000 steps
  ntwx=1000,     ! Write coords every 1000 steps
  cut=10.0,      ! Nonbonded cutoff
  ntb=1,         ! Constant volume
  ntp=0,         ! No pressure control
  iwrap=1,       ! Wrap coords in output
  ig=-1! Different random seed for each replica
/
EOF
}

# Function to create equilibration input file
create_eq_input() {
   local replica=$1
   local random_seed=$2
   cat << EOF > replica${replica}/eq.in
NPT equilibration at 310.15K - Replica $replica
&cntrl
  imin=0,        ! No minimization
  ntx=5,         ! Read coords and velocities
  irest=1,       ! Restart
  nstlim=25000000, ! 50 ns for eq
  dt=0.002,      ! Time step (ps)
  ntf=2,         ! SHAKE on H bonds
  ntc=2,         ! SHAKE on H
  temp0=310.15,  ! Target temp
  ntt=3,         ! Langevin thermostat
  gamma_ln=1.0,  ! Collision frequency
  ntpr=1000,     ! Print every 1000 steps
  ntwx=1000,     ! Write coords every 1000 steps
  cut=10.0,      ! Nonbonded cutoff
  ntb=2,         ! Constant pressure
  ntp=1,         ! Isotropic pressure scaling
  pres0=1.0,     ! Target pressure (bar)
  taup=2.0,      ! Pressure relaxation time
  iwrap=1,       ! Wrap coords in output
  ig=-1, ! Same random seed as heating
  ntr=0,         ! No restraints
/
EOF
}

# Function to create submission script
create_submission_script() {
   local replica=$1
   
   cat << EOF > replica${replica}/submit_replica${replica}.sh
#!/bin/bash
#SBATCH --job-name=MD_R${replica}
#SBATCH --output=MD_R${replica}_%j.out
#SBATCH --error=MD_R${replica}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-0-1,gpu-0-2,gpu-1-5,gpu-1-6

export PATH=/home3/apps/gcc/13.2.0/bin:$PATH
export LD_LIBRARY_PATH=/home3/apps/gcc/13.2.0/lib64:$LD_LIBRARY_PATH
export CUDA_HOME=/home8/apps/cuda-12.8
export CUDA_TOOLKIT_ROOT_DIR=/home8/apps/cuda-12.8
amber="/home3/apps/amber24/bin/pmemd.cuda"

# Store original working directory
WORKDIR=\$(pwd)

# Create unique temporary directory
TMPDIR="/tmp/\$USER/replica${replica}_\${SLURM_JOB_ID}"
mkdir -p \$TMPDIR

# Copy all necessary files to temporary directory
cp * \$TMPDIR/

# Change to temporary directory
cd \$TMPDIR

echo "Starting heating for replica${replica} on node \$(hostname)"
\$amber -O \\
  -i heat.in \\
  -o heat.out \\
  -p mobL.prmtop \\
  -c min2.rst \\
  -r heat.rst \\
  -x heat.nc

# Check if heating completed successfully
if [ \$? -ne 0 ]; then
  echo "ERROR: Heating failed for replica${replica}"
  # Copy logs back even on failure
  cp heat.out \$WORKDIR/
  exit 1
fi

echo "Starting equilibration for replica${replica}"
\$amber -O \\
  -i eq.in \\
  -o eq.out \\
  -p mobL.prmtop \\
  -c heat.rst \\
  -r eq.rst \\
  -x eq.nc

# Copy results back to working directory
echo "Copying results back to \$WORKDIR..."
cp heat.* eq.* \$WORKDIR/

# Clean up temporary directory
cd \$WORKDIR
rm -rf \$TMPDIR

echo "Replica${replica} completed successfully"
EOF
   chmod +x replica${replica}/submit_replica${replica}.sh
}

# Create replicas
for ((replica=1; replica<=NUM_REPLICAS; replica++)); do
   # Create random seed for this replica
   random_seed=$((replica * 100 + 123))
   
   # Create directory
   mkdir -p replica${replica}
   
   # Copy topology and coordinate files to replica directory
   cp mobL.prmtop min2.rst replica${replica}/
   
   # Create input files
   create_heating_input $replica $random_seed
   create_eq_input $replica $random_seed
   create_submission_script $replica
   
   echo "Configured replica${replica}"
done

# Submit all jobs
echo "Submitting jobs..."
submitted_jobs=()
for ((replica=1; replica<=NUM_REPLICAS; replica++)); do
   cd replica${replica}
   job_id=$(sbatch --parsable submit_replica${replica}.sh)
   submitted_jobs+=("$job_id")
   cd ..
   echo "Submitted replica${replica} (Job ID: $job_id)"
done

# Save job IDs to a file for reference
echo "${submitted_jobs[@]}" > replica_job_ids.txt
echo "All job IDs saved to replica_job_ids.txt"

# Create a summary script to check status later
cat << EOF > check_replica_status.sh
#!/bin/bash
echo "Replica Status Summary:"
echo "----------------------"
for replica in \$(seq 1 $NUM_REPLICAS); do
  # Check job status using sacct
  job_id=\$(grep -o "\b[0-9]\+" replica\${replica}/MD_R\${replica}_*.out 2>/dev/null | head -1)
  if [ -z "\$job_id" ]; then
    job_status="UNKNOWN"
  else
    job_status=\$(sacct -j \$job_id --format=State --noheader | head -1 | tr -d ' ')
  fi
  
  # Check file status
  if [ -f "replica\${replica}/eq.rst" ]; then
    file_status="COMPLETE"
  elif [ -f "replica\${replica}/heat.rst" ]; then
    file_status="HEATING DONE"
  else
    file_status="NOT STARTED/FAILED"
  fi
  
  echo "Replica \$replica: Job \$job_status, Files \$file_status"
done
EOF
chmod +x check_replica_status.sh

echo "Created check_replica_status.sh to monitor progress"
echo "To monitor all jobs, run: squeue -u \$USER"
