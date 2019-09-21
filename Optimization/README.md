### Lake Problem DPS vs. Intertemporal Optimization

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

Contents:
`DPS`: directory containing main-lake-DPS.cpp (the C++ code used to optimize the four objective lake management problem using direct policy search), as well as its makefile and a job submission script

`Intertemporal`: contains main-lake-IT.cpp (the C++ code used to optimize the four objective lake management problem using an intertemporal, open loop solution strategy), as well as its makefile and a job submission script

`SOWs_Type6.txt`: file containing random natural P inflows generated from a log-normal distribution with mean 0.03 and variance 10^-5

`boostutil.h`: file containing utility functions for boost matrices and vectors

`moeaframework.c`,`moeaframework.h`: provides methods in the lake problem C++ code for use with the MOEAFramework, which are necessary to compile the code as written

`find_refSets.sh`, `find_runtime_metrics.sh`, `get_objs.sh`: scripts for processing the MOEA output

To compile and run:
* Make sure you have the [MOEAFramework](http://www.moeaframework.org). Download the MOEAFramework-\*-Demo.jar file and copy it to this directory (`Lake_Problem_DPS/Optimization`).

* Download the [Borg](http://borgmoea.org/) source code. Make a new directory in `Lake_Problem_DPS/Optimization/` called `borg` and copy the Borg source code to it. Also add `moeaframework.c` and `moeaframework.h` to the `/borg` directory. To compile Borg and the Lake Problem (two steps down for DPS and three steps down for Intertemporal), you may need to edit line 260 of borgms.c from `BORG_Master_initialize_preseed() {` to `void BORG_Master_initialize_preseed() {`

* Download [Pareto.py](https://github.com/matthewjwoodruff/pareto.py) and put it in this directory (`Lake_Problem_DPS/Optimization/`).

* From this directory (`Lake_Problem_DPS/Optimization/`), compile and run the DPS optimization code. This is set up to run in parallel using MPI. You will need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`cd DPS`,   
`module load gcc openmpi boost/1.54.0`,   
`make`,   
`mkdir runtime`,   
`mkdir sets`,   
`sbatch run_DPS_opt.sh`.   
You may need to change the include statement on line 4 of the makefile to match the location of the Boost C++ library on your machine. Within `run_DPS_opt.sh` you should change the directory names on lines 2 and 3, and your email address on line 9. You can also change the number of nodes and processors on lines 4 and 5, but make sure to then scale the walltime on line 8 up or down, accordingly. Note it may take awhile to run.

* From this directory (`Lake_Problem_DPS/Optimization/`), compile and run the intertemporal optimization code. This is set up to run in parallel using MPI. You will need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`cd Intertemporal && make`,   
`mkdir runtime`,   
`mkdir sets`,   
`sbatch run_IT_opt.sh`.   
You may need to change the include statement on line 4 of the makefile to match the location of the Boost C++ library on your machine. Within `run_IT_opt.sh` you should change the directory names on lines 2 and 3, and your email address on line 9. You can also change the number of nodes and processors on lines 4 and 5, but make sure to then scale the walltime on line 8 up or down, accordingly. Note it may take awhile to run.

* Once the DPS and intertemporal optimization jobs have finished running, return to this directory (`Lake_Problem_DPS/Optimization/`). Next, find the reference sets across the 50 seeds of each solution strategy, and calculate runtime metrics. You will again need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`mkdir DPS/metrics`,   
`mkdir DPS/objs`,   
`mkdir Intertemporal/metrics`,   
`mkdir Intertemporal/objs`,    
`sh get_objs.sh`,   
`sh find_refSets.sh`,  
`mkdir output`, 
`sbatch find_runtime_metrics.sh`.   
`find_runtime_metrics.sh` runs multiple jobs in parallel. You should change the directory names on lines 2 and 3. Lines 11 and 15 may also need to be changed to correspond to the version of the MOEAFramework you are using. Finally, line 3 of `find_refSets.sh` may be unnecessary, or need to be changed, depending on the machine you are using.

Next, go to the Re-evaluation directory (`cd ./../Re-evaluation`) to re-evaluate the optimized policies on alternative SOWs.
