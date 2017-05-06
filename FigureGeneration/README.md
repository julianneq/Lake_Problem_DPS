### Lake Problem DPS vs. Intertemporal Figure Generation

Intended for use with the MOEAFramework, Borg and pareto.py. Licensed under the GNU Lesser General Public License.

Contents:
Separate Python scripts to generate each of the 11 figures in Quinn et al. (2017).

To make the figures run `python makeAllFigures.py`. This should generate 11 pdfs of the figures in the paper. 
There will be slight differences if you run it on your own optimization set. The seed used to sample synthetic natural P inflows in the C++ optimization code is based on the time of day, and therefore different reference sets will be found depending on when the optimizaiton is run.

If you want to generate the same figures as in the paper, move the data from `./../DataInPaper` to the following directories:
`mv DPS/DPS.reference ./../Optimization`,   
`mv DPS/DPS.resultfile ./../Optimization`,   
`mv DPS/DPSrobustness.txt ./../Re-evaluation`,   
`mv DPS/metrics/* ./../Optimizaiton/DPS/metrics`,   
`mv DPS/re-eval/* ./../Re-evaluation/DPS/output`, and   
   
`mv Intertemporal/Intertemporal.reference ./../Optimization`,   
`mv Intertemporal/Intertemporal.resultfile ./../Optimization`,   
`mv Intertemporal/ITrobustness.txt ./../Re-evaluation`,   
`mv Intertemporal/metrics/* ./../Optimizaiton/Intertemporal/metrics`,   
`mv Intertemporal/re-eval/* ./../Re-evaluation/Intertemporal/output`.   

Additionally, both methods find several solutions with perfect reliability, but only one was plotted in the paper. To plot the same policies used in the figure, make the following changes to the codes in this directory:   
To `makeFigure2.py`, change line 15 to `DPSmostRel = 26`,   
to `makeFigure4.py`, change line 90 to `IT_best_rel = 8` and line 98 to `DPS_best_rel = 26`, and   
to `makeFigure6.py` and `makeFiure7.py`, change line 33 to `ITmostRel = 8` and line 34 to `DPSmostRel = 26`.   
