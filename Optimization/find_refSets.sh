#/bin/bash

source /etc/profile.d/modules.sh
module load python-2.7.5

python pareto.py ./DPS/sets/*.set -o 6-9 -e 0.01 0.01 0.001 0.001 --output DPS.resultfile --delimiter=" " --comment="#"
cut -d ' ' -f 7-10 DPS.resultfile >DPS.reference

python pareto.py ./Intertemporal/sets/*.set -o 100-103 -e 0.01 0.01 0.001 0.001 --output Intertemporal.resultfile --delimiter=" " --comment="#"
cut -d ' ' -f 101-104 Intertemporal.resultfile >Intertemporal.reference

python pareto.py ./*.reference -o 0-3 -e 0.01 0.01 0.001 0.001 --output Overall.reference --delimiter=" " --comment="#"