#!/bin/bash
NSEEDS=50
SEEDS=$(seq 1 ${NSEEDS})

for SEED in ${SEEDS}
do
	awk 'BEGIN {FS=" "}; /^#/ {print $0}; /^[^#/]/ {printf("%s %s %s %s\n",$7,$8,$9,$10)}' ./DPS/runtime/LakeDPS_S${SEED}.runtime \
	 	>./DPS/objs/LakeDPS_S${SEED}.obj

	awk 'BEGIN {FS=" "}; /^#/ {print $0}; /^[^#/]/ {printf("%s %s %s %s\n",$101,$102,$103,$104)}' ./Intertemporal/runtime/LakeIT_S${SEED}.runtime \
		>./Intertemporal/objs/LakeIT_S${SEED}.obj
done