#!/bin/bash
#SBATCH --nodes=17             # specify number of nodes
#SBATCH --ntasks-per-node=16  # specify number of core per node
#SBATCH --export=ALL
#SBATCH -t 1:00:00            # set max wallclock time
#SBATCH --job-name="IT_reevaluation" # name your job 
#SBATCH --output="IT_reevaluation.out"

module load python
mpirun python3 resimulateIT.py
