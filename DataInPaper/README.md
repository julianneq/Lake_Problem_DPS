### Lake Problem DPS vs. Intertemporal Study (C++, Python, Bash)
This repository contains all of the code used for the study described in the following submitted paper:  

Quinn, J. D., P. M. Reed, and K. Keller, "Direct Policy Search for Robust Multi-Objective Management of Deeply Uncertain Socio-Ecological Tipping Points", Environmental Modelling and Software, In Review.

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

Contents:
The actual data used in Quinn et al. (In Review)

* `DPS`: Directory containing 
1) Reference set (`DPS.reference`) and result file (`DPS.resultfile`) found for the 4-objective version of the lake problem using DPS,
2) Metrics of the 50 DPS seeds in the subdirectory `metrics`,
3) The re-evaluated objective function values in the 1000 alternative SOWs in the subdirectory `re-eval`, and
4) Values of the domain satisficing criterion for each policy on the 3 criteria considered in the paper (`DPSrobustness.txt`).

* `Intertemporal`: Directory containing 
1) Reference set (`Intertemporal.reference`) and result file (`Intertemporal.resultfile`) found for the 4-objective version of the lake problem using DPS,
2) Metrics of the 50 intertemporal seeds in the subdirectory `metrics`,
3) The re-evaluated objective function values in the 1000 alternative SOWs in the subdirectory `re-eval`, and
4) Values of the domain satisficing criterion for each policy on the 3 criteria considered in the paper (`ITrobustness.txt`).

* `LHsamples.txt`: file containing the alternative SOWs used for the re-evaluation. Columns are b, q, mu, sigma, delta and X0, in that order.

* `Overall.reference`: reference set across all DPS and Intertemporal Seeds
