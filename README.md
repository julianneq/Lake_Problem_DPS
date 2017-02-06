### Lake Problem DPS vs. Intertemporal Study (C++, Python, Bash)
This repository contains all of the code used for the study described in the following submitted paper: Quinn, J. D., P. M. Reed, and K. Keller, "Direct Policy Search for Robust Multi-Objective Management of Deeply Uncertain Socio-Ecological Tipping Points", Environmental Modelling and Software, In Review.

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

Contents:

* `DataInPaper`: Directory containing all of the data generated in this paper, including reference sets from the optimizaiton, SOWs used for re-evaluation, and the domain satisficing metric for each policy

* `Optimization`: Directory containing the code for optimizing the lake problem using DPS and the intertemporal open loop control strategy

* `Re-evaluation`: Directory containing the code for re-evaluating the policies found from the above optimization on a set of alternative SOWs and calculating the satisficing metric for each policy

* `FigureGeneration`: Code for generating the figures found in Quinn et al. (In Review)

To reproduce the results from this study, first follow the steps given in the README of the `Optimization` directory, then those in `Re-evaluation`, and finally those in `FigureGeneration.`
