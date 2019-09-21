#!/bin/bash
#SBATCH -D /scratch/jdq6nn/Lake_Problem_DPS/Optimization/DPS/		# working directory
#SBATCH -o /scratch/jdq6nn/Lake_Problem_DPS/Optimization/DPS/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 3            					# Total number of nodes to request (up to 120)
#SBATCH --ntasks-per-node 20           		# Number of processors per node (up to 20)
#SBATCH -p parallel           				# Queue name "parallel"
#SBATCH -A quinnlab       					# allocation name
#SBATCH -t 1:00:00       					# Run time (hh:mm:ss) - up to 36 hours
#SBATCH --mail-user=jdq6nn@virginia.edu     # address for email notification
#SBATCH --mail-type=ALL                  	# email at Begin and End of job

module load gcc openmpi boost/1.54.0

# Your commands go here
# arguments are <seed> <NFE>
for i in {1..50}
do
  mpirun ./LakeDPSparallel $i 200000
done
