#!/bin/bash
#SBATCH -D /scratch/jdq6nn/Lake_Problem_DPS/Optimization	 # working directory
#SBATCH -o /scratch/jdq6nn/Lake_Problem_DPS/Optimization/output/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -A quinnlab				# allocation name
#SBATCH -t 1:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem-per-cpu=12288		# Memory per cpu (bytes)
#SBATCH --array=1-50			# Array of jobs to loop through

java -cp MOEAFramework-2.12-Demo.jar org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
	-d 4 -i ./DPS/objs/LakeDPS_S${SLURM_ARRAY_TASK_ID}.obj -r Overall.reference \
	-o ./DPS/metrics/LakeDPS_S${SLURM_ARRAY_TASK_ID}.metrics

java -cp MOEAFramework-2.12-Demo.jar org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
	-d 4 -i ./Intertemporal/objs/LakeIT_S${SLURM_ARRAY_TASK_ID}.obj -r Overall.reference \
	-o ./Intertemporal/metrics/LakeIT_S${SLURM_ARRAY_TASK_ID}.metrics
