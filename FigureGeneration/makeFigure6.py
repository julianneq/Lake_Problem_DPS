import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
import scipy.stats as ss

# Set the number of RBFs (n), decision variables, objectives and constraints, and random seeds
n = 2
nSeeds = 1

# Set lake model parameters, number of years simulated, and number of samples generated
q = 2
b = 0.42
alpha = 0.4
delta = 0.98
mu = 0.03
sigma = np.sqrt(10**(-5.0))
lake0 = 0
nYears = 100
nSamples = 100

# Establish reliability and inertia thresholds
reliability_threshold = 0.85

def fun(x):
    return [(x[0]**q)/(1+x[0]**q) - b*x[0]]
soln = root(fun, 0.75)
pCrit = soln.x

def makeFigure6():
    IT = np.loadtxt('./../Optimization/Intertemporal.resultfile')
    DPS = np.loadtxt('./../Optimization/DPS.resultfile')
    
    ITmostRel = np.argmin(IT[:,103]) # one of several solutions with perfect reliability; in paper, row 8 of IT.resultfile was plotted
    DPSmostRel = np.argmin(DPS[:,9]) # one of several solutions with perfect reliability; in paper, row 26 of DPS.resultfile was plotted
    
    # plot lake P concentration vs. time
    seed=1
    lake_state3 = LakeModel_IT(seed,IT[ITmostRel,0:100]) # intertemporal best reliability
    lake_state4 = LakeModel_IT(seed,IT[np.argmin(IT[:,100]),0:100]) # intertemporal best benefits
       
    # plot P discharge vs. P concentration in the lake
    lake_state = np.arange(0,2.5,0.01)
    Y1 = np.zeros(len(lake_state))
    Y2 = np.zeros(len(lake_state))
    for i in range(len(lake_state)):
        Y1[i] = DPSpolicy(lake_state[i], DPS[DPSmostRel,0:6]) # DPS best reliability
        Y2[i] = DPSpolicy(lake_state[i], DPS[np.argmin(DPS[:,6]),0:6]) # DPS best benefits
        
    fig = plt.figure()
    ax = fig.add_subplot(1,3,1)
    line1, = ax.plot(lake_state, Y1, c="#08519c",linewidth=2) # DPS best reliability
    line2, = ax.plot(lake_state, Y2, c="#006d2c",linewidth=2) # DPS best benefits
    line5, = ax.plot([pCrit,pCrit],[0,0.1],c='#a50f15',linewidth=2) # critical P threshold
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.15)
    ax.set_ylabel('Anthropogenic P Release, $a_t$',fontsize=16)
    ax.set_yticks(np.arange(0,0.16,0.04))
    ax.tick_params(axis='both',labelsize=14)
    ax.text(0.03,0.125,'a) Reliability & Benefits-Maximizing\nDPS Policies', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height*0.8])
    
    ax = fig.add_subplot(1,3,2)
    line1, = ax.plot(lake_state, Y1, c="#08519c",linewidth=2) # DPS best reliability
    line3 = ax.scatter(lake_state3[0:-1],IT[ITmostRel,0:100],c="#08519c",s=20,linewidth=0) # intertemporal best reliability
    line5, = ax.plot([pCrit,pCrit],[0,0.1],c='#a50f15',linewidth=2) # critical P threshold
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.15)
    ax.set_xlabel('Lake P Concentration, $X_t$',fontsize=16)
    ax.tick_params(axis='both',labelsize=14)
    ax.text(0.03,0.13,'b) Reliability-Maximizing Policies', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height*0.8])
    
    ax = fig.add_subplot(1,3,3)
    line2, = ax.plot(lake_state, Y2, c="#006d2c",linewidth=2) # DPS best benefits
    line4 = ax.scatter(lake_state4[0:-1],IT[np.argmin(IT[:,100]),0:100],c="#006d2c",s=20,linewidth=0) # intertemporal best benefits
    line5, = ax.plot([pCrit,pCrit],[0,0.1],c='#a50f15',linewidth=2) # critical P threshold
    ax.set_xlim(0,np.max(lake_state4))
    ax.set_ylim(0,0.15)
    ax.tick_params(axis='both',labelsize=14)
    ax.text(0.05,0.13,'c) Benefits-Maximizing Policies', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height*0.8])
    
    plt.suptitle('Optimized DPS Policies vs.Simulated Intertemporal Policies',fontsize=16)
    plt.figlegend([line5,line1,line2,line3,line4],['Critical P Threshold','DPS Highest Reliability',\
        'DPS Highest Benefits','Intertemporal Highest Reliability','Intertemporal Highest Benefits'], \
        loc = 'lower center', ncol=2)
    fig.set_size_inches([18.25, 6.2875])
    fig.savefig('Figure6.pdf')
    fig.clf()
    
    return None
    
def LakeModel_IT(seed, vars):
    # Set inflow distribution parameters
    log_std = np.sqrt(np.log(1+sigma**2/mu**2))
    log_mu = np.log(mu) - 0.5*(log_std**2)
    
    # Initialize arrays to store P level in the lake at each time step
    lake_state = np.zeros([nYears+1])

    # Randomly generate nSamples of nYears of natural P inflows
    natFlow = np.zeros([nYears])
    np.random.seed(seed)
    natFlow= np.exp(ss.norm.rvs(log_mu, log_std, nYears))
    
    # Run lake model simulation
    lake_state[0] = lake0
    for i in range(nYears):
        lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]**q)/(1+(lake_state[i]**q)) + vars[i] + natFlow[i]
                
    return lake_state
    
def DPSpolicy(lake_state, vars):
    # Determine centers, radii and weights of RBFs
    C = vars[0::3]
    B = vars[1::3]
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
    
    # Determine pollution emission decision, Y
    Y = 0
    for i in range(len(C)):
        if B[i] != 0:
            Y = Y + W[i]*((np.absolute(lake_state-C[i])/B[i])**3)
            
    Y = min(0.1,max(Y,0.01))
    
    return Y
    
makeFigure6()