# Python implementation of the 4 objective lake problem formulation

# Julianne Quinn
# Cornell University
# jdq2101@gmail.com
# July 9, 2015

# Lake parameters
# b    : proportion of phosphorous retained in the lake each year;
# q    : steepness of the sigmoid curve, large values give a steeper slope;

# Parameters related to utiltiy
# delta : discount rate, set at 0.98
# alpha : utility from pollution, set at 0.4

# Lake state variabiles
# lake_state : phosphorous concentration in the lake at a given time step;
# initially set to 0

# Stochasticity is introduced by natural variability around anthropogenic pollution flow
# which is generated form a log-normal distribution

# nYears : the time horizon for planning (100 years in this example)

# Decision variable
# vars : vector of length nYears giving the amount of pollution emitted each year

# Objectives
# 1: average discounted net benefit
# 2: maximum average P concentration
# 3: average inertia
# 4: average reliability

# Constraints
# Reliability must be > 85%

import numpy as np
from scipy.optimize import root
import scipy.stats as ss

# Set the number of decision variables and objectives
nvars = 100
nobjs = 4

# Set lake model parameters, number of years simulated, and number of samples generated
alpha = 0.4
nYears = 100
nSamples = 100

# Establish reliability and inertia thresholds
reliability_threshold = 0.85
inertia_threshold = -0.02

def LakeModel_IT(seed, vars, b=0.42, q=2.0, mu=0.03, sigma=np.sqrt(10**(-5)), \
    delta=0.98, lake0=0):
    # Compute Pcrit given b and q
    def fun(x):
        return [(x[0]**q)/(1+x[0]**q) - b*x[0]]
    soln = root(fun, 0.75)
    Pcrit = soln.x
    
    # Initialize arrays to store objective function values and constraints
    objs = [0.0]*nobjs
    
    # Set inflow distribution parameters
    log_sigma = np.sqrt(np.log(1+sigma**2/mu**2))
    log_mu = np.log(mu) - 0.5*(log_sigma**2)
    
    # Initialize arrays to store average daily P and performance of solution in each generated sample
    average_annual_P = np.zeros([nYears])
    discounted_benefit = np.zeros([nSamples])
    yrs_inertia_met = np.zeros([nSamples])
    yrs_Pcrit_met = np.zeros([nSamples])
    lake_state = np.zeros([nYears+1])
    
    # Randomly generate nSamples of nYears of natural P inflows
    natFlow = np.zeros([nSamples,nYears])
    for i in range(nSamples):
        np.random.seed(seed+i)
        natFlow[i,:] = np.exp(ss.norm.rvs(log_mu, log_sigma, nYears))
    
    # Run lake model simulation
    for s in range(nSamples):
        lake_state[0] = lake0
        
        for i in range(nYears):
            lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]**q)/(1+(lake_state[i]**q)) + vars[i] + natFlow[s,i]
            average_annual_P[i] = average_annual_P[i] + lake_state[i+1]/nSamples
            discounted_benefit[s] = discounted_benefit[s] + alpha*vars[i]*delta**i
            
            if i>=1 and ((vars[i] - vars[i-1]) > inertia_threshold):
                yrs_inertia_met[s] = yrs_inertia_met[s] + 1
                
            if lake_state[i+1] < Pcrit:
                yrs_Pcrit_met[s] = yrs_Pcrit_met[s] + 1
                
    # Calculate minimization objectives (defined in comments at beginning of file)
    objs[0] = np.mean(discounted_benefit) #average economic benefit
    objs[2] = np.max(average_annual_P) #maximum average annual P concentration
    objs[3] = np.sum(yrs_inertia_met)/((nYears-1)*nSamples) #average pct of transitions meeting inertia thershold
    objs[3] = np.sum(yrs_Pcrit_met)/(nYears*nSamples) #average reliability
    
    return objs
