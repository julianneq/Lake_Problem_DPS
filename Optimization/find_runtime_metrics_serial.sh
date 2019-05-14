NSEEDS=50
SEEDS=$(seq 1 ${NSEEDS})
JAVA_ARGS="-cp MOEAFramework-2.12-Demo.jar"

for SEED in ${SEEDS}
do
	NAME=Runtime_Metrics_S${SEED}
	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
		-d 4 -i ./DPS/objs/LakeDPS_S${SEED}.obj -r Overall.reference \
		-o ./DPS/metrics/LakeDPS_S${SEED}.metrics
	java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator \
		-d 4 -i ./Intertemporal/objs/LakeIT_S${SEED}.obj -r Overall.reference \
		-o ./Intertemporal/metrics/LakeIT_S${SEED}.metrics
done
