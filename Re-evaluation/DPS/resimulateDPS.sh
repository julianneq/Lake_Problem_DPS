#!/bin/bash
#SBATCH -N 5            																# Total number of nodes to request (up to 120)
#SBATCH --ntasks-per-node 20           													# Number of processors per node (up to 20)
#SBATCH -p parallel           															# Queue name "parallel"
#SBATCH -A quinnlab       																# allocation name
#SBATCH -t 0:30:00       											 					# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=js4yd@virginia.edu             										# address for email notification
#SBATCH --mail-type=ALL                  												# email at Begin and End of job

cd $SLURM_SUBMIT_DIR

module load gcc/system openmpi/2.1.5 python/2.7.14 mpi4py/3.0.0-py2.7
mpirun python resimulateDPS.py