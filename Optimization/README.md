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

* Download the [Borg](http://borgmoea.org/) source code. Make a new directory in `Lake_Problem_DPS/Optimization/` called `borg` and copy the Borg source code to it. Also add `moeaframework.c` and `moeaframework.h` to the `/borg` directory.

* Download [Pareto.py](https://github.com/matthewjwoodruff/pareto.py) and put it in this directory (`Lake_Problem_DPS/Optimization/`).

* From this directory (`Lake_Problem_DPS/Optimization/`), compile and run the DPS optimization code. This is set up to run in parallel using MPI. You will need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`cd DPS && make`,   
`mkdir runtime`,   
`mkdir sets`,   
`qsub run_DPS_opt.sh`.   
You may need to change the include statement on line 4 of the makefile to match the location of the Boost C++ library on your machine. You can also change the number of nodes and processors on line 3 of `run_DPS_opt.sh`. Make sure to also scale the walltime on line 2 up or down, accordingly. Note it may take awhile to run.

* From this directory (`Lake_Problem_DPS/Optimization/`), compile and run the intertemporal optimization code. This is set up to run in parallel using MPI. You will need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`cd Intertemporal && make`,   
`mkdir runtime`,   
`mkdir sets`,   
`qsub run_IT_opt.sh`.   
You may need to change the include statement on line 4 of the makefile to match the location of the Boost C++ library on your machine. You can also change the number of nodes and processors on line 3 of `run_IT_opt.sh`. Make sure to also scale the walltime on line 2 up or down, accordingly. Note it may take awhile to run. If there are enough processors available on your machine, you can run the DPS and intertemporal optimization at the same time.

* Once the DPS and intertemporal optimization jobs have finished running, return to this directory (`Lake_Problem_DPS/Optimization/`). Next, find the reference sets across the 50 seeds of each solution strategy, and calculate runtime metrics. You will again need to make directories for the output by running the following commands from `Lake_Problem_DPS/Optimization/`:   
`mkdir DPS/metrics`,   
`mkdir DPS/objs`,   
`mkdir Intertemporal/metrics`,   
`mkdir Intertemporal/objs`,   
`mkdir output`,   
`mkdir error`,   
`sh get_objs.sh`,   
`sh find_refSets.sh`,   
`sh find_runtime_metrics.sh`.   
`find_runtime_metrics.sh` runs multiple jobs in parallel. Line 3 may need to be changed to correspond to the version of the MOEAFramework you are using. Lines 3 and 4 of `find_refSets.sh` may be unnecessary, or need to be changed, depending on the machine you are using.

Next, go to the Re-evaluation directory (`cd ./../Re-evaluation`) to re-evaluate the optimized policies on alternative SOWs.
