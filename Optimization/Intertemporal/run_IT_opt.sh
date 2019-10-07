#!/bin/bash
#SBATCH --nodes=4             # specify number of nodes
#SBATCH --ntasks-per-node=16  # specify number of core per node
#SBATCH --export=ALL
#SBATCH -t 1:00:00            # set max wallclock time
#SBATCH --job-name="IT_opt" # name your job 
#SBATCH --output="IT_out.out"

# Your commands go here
# arguments are <seed> <NFE>
for i in {1..50}
do
  mpirun ./LakeITparallel $i 200000
done

