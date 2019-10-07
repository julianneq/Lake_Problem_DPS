NSEEDS=50
SEEDS=$(seq 1 ${NSEEDS})
JAVA_ARGS="-cp MOEAFramework-2.12-Demo.jar"

for SEED in ${SEEDS}
do
	NAME=Runtime_Metrics_S${SEED}
	SLURM="
  #!/bin/bash\n\
  #SBATCH --job-name=${NAME}\n\
	#SBATCH --nodes=1\n\
  #SBATCH --ntasks-per-node=1\n\
	#SBATCH -t 1:00:00\n\
	#SBATCH --output=/output/${NAME}\n\
	#SBATCH --error=/error/${NAME}\n\
	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
		-d 4 -i ./DPS/objs/LakeDPS_S${SEED}.obj -r Overall.reference \
		-o ./DPS/metrics/LakeDPS_S${SEED}.metrics\n\
	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
		-d 4 -i ./Intertemporal/objs/LakeIT_S${SEED}.obj -r Overall.reference \
		-o ./Intertemporal/metrics/LakeIT_S${SEED}.metrics"
	echo -e $SLURM | sbatch
done