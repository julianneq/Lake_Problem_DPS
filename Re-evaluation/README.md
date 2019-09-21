### Lake Problem DPS vs. Intertemporal Robustness Analysis

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

Contents:
`DPS`: directory containing `resimulateDPS.py` and `LakeModel_DPS.py` (Python codes used to re-evalute the DPS-optimized policies for the four objective lake management problem), as well as a job submission script, `resimulateDPS.sh`

`Intertemporal`: directory containing `resimulateIT.py` and `LakeModel_IT.py` (Python codes used to re-evalute the intertemporal-optimized policies for the four objective lake management problem), as well as a job submission script, `resimulateIT.sh`

`calcRobustness.py`: Python code used to calculate the domain criterion satisficing metrics for each policy after re-evaluating them on alternative SOWs

`Lake_params.txt`: file listing the uncertain parameters and ranges they are sampled over to generate alternative SOWs

`sample_parameters.sh`: job script to generate 1000 alternative SOWs for the re-evaluation. This should create a file called `LHsamples.txt`.

To run the re-evaluation:
* Generate alternative SOWs (`LHsamples.txt`) in which to evaluate the different policies by running the following command (you may need to change the version of the MOEAFramework given on line 3):   
`sh sample_parameters.sh`

* Re-evaluate the DPS and intertemporal policies. You will need to make a directory for the output. This is set up to run in parallel using mpi4py. From this directory, run the following commands:   
`mkdir DPS/output`,   
`mkdir Intertemporal/output`,   
`cd DPS && sbatch resimulateDPS.sh`,   
`cd ./../Intertemporal && sbatch resimulateIT.sh`.   
Within `resimulateDPS.sh` and `resimulateIT.sh` you should change the directory names on lines 2 and 3, and your email address on line 9. You can also change the number of nodes and processors on lines 4 and 5, but make sure to then scale the walltime on line 8 up or down, accordingly. If necessary, change line 12 for your machine.

* Next, calculate the domain satisficing criterion for the policies found by each solution strategy. From this directory run the following commands:   
`module load gcc/5.4.0 openmpi/3.1.4 python/2.7.16`   
`python calcRobustness.py`.   
This should write `DPSrobustness.txt` and `ITrobustness.txt` to this directory.

Next, go to the FigureGeneration directory (`cd ./../FigureGeneration`) to generate the plots found in [Quinn et al. (2017)](https://doi.org/10.1016/j.envsoft.2017.02.017)
