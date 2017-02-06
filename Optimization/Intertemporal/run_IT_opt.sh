#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=4:ppn=16
#PBS -j oe

cd $PBS_O_WORKDIR

# Your commands go here
# arguments are <seed> <NFE>
for i in {1..50}
do
  mpirun ./LakeITparallel $i 200000
done
