# Python implementation of the 4 objective lake problem formulation
# using Direct Policy Search, optimized by Borg

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

# Stochasticity is introduced by natural variability around anthropogenic pollution flow, 
# which is generated form a log-normal distribution 

# nYears : the time horizon for planning (100 years in this example)

# Decision variable
# vars : vector of length 3n describing the policy parameters
# n is the number of RBFs, each with a weight, center and radius
# format of vars is [c1, b1, w1, c2, b2, w2, ..., cn, bn, wn]
# all variables must be between 0 and 1, and all weights must sum to 1
# The weights determined by Borg aren't constrained to sum to 1,
# but are normalized to do so in the model simulation
# The actual decision, the amount of pollution emitted by the town, is Y,
# a function of vars and the current state of the lake (lake_state)

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

# Set the number of RBFs (n), decision variables, objectives and constraints
n = 2
nvars = 3*n
nobjs = 4

# Set lake model parameters, number of years simulated, and number of samples generated
alpha = 0.4
nYears = 100
nSamples = 100

# Establish reliability and inertia thresholds
reliability_threshold = 0.85
inertia_threshold = -0.02

def LakeModel_DPS(seed, vars, b=0.42, q=2.0, mu=0.03, sigma=np.sqrt(10**(-5)), \
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
        
    # Determine centers, radii and weights of RBFs
    C = vars[0::3]
    R = vars[1::3]
    W = vars[2::3]
    newW = np.zeros(len(W))
    
    # Normalize weights to sum to 1
    total = sum(W)
    if total != 0.0:
        for i in range(len(W)):
            newW[i] = W[i]/total
    else:
        for i in range(len(W)):
            newW[i] = 1/n
    
    # Run lake model simulation
    for s in range(nSamples):
        lake_state[0] = lake0
        Y = np.zeros([nYears])
        #find policy-derived emission
        Y[0] = RBFpolicy(lake_state[0], C, R, newW)

        for i in range(nYears):
            lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]**q)/(1+(lake_state[i]**q)) + Y[i] + natFlow[s,i]
            average_annual_P[i] = average_annual_P[i] + lake_state[i+1]/nSamples
            discounted_benefit[s] = discounted_benefit[s] + alpha*Y[i]*delta**i
                
            if i>=1 and ((Y[i] - Y[i-1]) > inertia_threshold):
                yrs_inertia_met[s] = yrs_inertia_met[s] + 1
                
            if lake_state[i+1] < Pcrit:
                yrs_Pcrit_met[s] = yrs_Pcrit_met[s] + 1
                
            if i<(nYears-1):
                #find policy-derived emission
                Y[i+1] = RBFpolicy(lake_state[i+1], C, R, newW)
                
    # Calculate minimization objectives (defined in comments at beginning of file)
    objs[0] = np.mean(discounted_benefit) #average economic benefit
    objs[1] = np.max(average_annual_P) #maximum average annual P concentration
    objs[2] = np.sum(yrs_inertia_met)/((nYears-1)*nSamples) #average pct of transitions meeting inertia thershold
    objs[3] = np.sum(yrs_Pcrit_met)/(nYears*nSamples) #average reliability

    return objs

def RBFpolicy(lake_state, C, R, W):    
    # Determine pollution emission decision, Y
    Y = 0
    for i in range(len(C)):
        if R[i] != 0:
            Y = Y + W[i]*((np.absolute(lake_state-C[i])/R[i])**3)

    Y = min(0.1,max(Y,0.01))
    
    return Y
