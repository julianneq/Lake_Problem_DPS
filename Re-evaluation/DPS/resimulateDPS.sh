#!/bin/bash
#SBATCH -D /scratch/jdq6nn/Lake_Problem_DPS/Re-evaluation/DPS/		# working directory
#SBATCH -o /scratch/jdq6nn/Lake_Problem_DPS/Re-evaluation/DPS/output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 5            					# Total number of nodes to request (up to 120)
#SBATCH --ntasks-per-node 20           		# Number of processors per node (up to 20)
#SBATCH -p parallel           				# Queue name "parallel"
#SBATCH -A quinnlab       					# allocation name
#SBATCH -t 0:30:00       					# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=jdq6nn@virginia.edu     # address for email notification
#SBATCH --mail-type=ALL 					# email at Begin and End of job

module load gcc openmpi python mpi4py
mpirun python resimulateDPS.py
