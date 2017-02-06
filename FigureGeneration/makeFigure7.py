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

def makeFigure7():
    IT = np.loadtxt('./../Optimization/Intertemporal.resultfile')
    DPS = np.loadtxt('./../Optimization/DPS.resultfile')
    
    ITmostRel = np.argmin(IT[:,103]) # one of several solutions with perfect reliability; in paper, row 8 of IT.resultfile was plotted
    DPSmostRel = np.argmin(DPS[:,9]) # one of several solutions with perfect reliability; in paper, row 26 of DPS.resultfile was plotted
    
    # plot lake P concentration vs. time
    seed=1
    lake_state1, Y1 = LakeModel_DPS(seed,DPS[DPSmostRel,0:6]) # DPS best reliability
    lake_state2, Y2 = LakeModel_DPS(seed,DPS[np.argmin(DPS[:,6]),0:6]) # DPS best benefits
    lake_state3 = LakeModel_IT(seed,IT[ITmostRel,0:100]) # intertemporal best reliability
    lake_state4 = LakeModel_IT(seed,IT[np.argmin(IT[:,100]),0:100]) # intertemporal best benefits
    
    fig = plt.figure()
    ax = fig.add_subplot(2,3,1)
    time = np.arange(0,len(lake_state1),1)
    line1, = ax.plot(time, lake_state1, c='#08519c',linewidth=2) # DPS best reliability
    line2, = ax.plot(time, lake_state2, c='#006d2c', linewidth=2) # DPS best benefits
    line5, = ax.plot(time, np.ones(len(time))*pCrit, c='#a50f15', linewidth=2) # critical threshold
    ax.set_ylabel('Lake P Concentration, $X_t$',fontsize=16)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.8)
    ax.set_yticks(np.arange(0,1.0,0.2))
    ax.text(5,0.65,'a) DPS Reliability & Benefits-\nMaximizing Lake P Trajectories', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
    
    ax = fig.add_subplot(2,3,2)
    line1, = ax.plot(time, lake_state1, c='#08519c',linewidth=2) # DPS best reliability
    line3, = ax.plot(time, lake_state3, c='#08519c', linewidth=2, linestyle='--') # intertemporal best reliability
    line5, = ax.plot(time, np.ones(len(time))*pCrit, c='#a50f15', linewidth=2) #critical threshold
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.8)
    ax.set_yticks(np.arange(0,1.0,0.2))
    ax.text(5,0.65,'b) Reliability-Maximizing\nLake P Trajectories', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
    
    ax = fig.add_subplot(2,3,3)
    line2, = ax.plot(time, lake_state2, c='#006d2c', linewidth=2) # DPS best benefits
    line4, = ax.plot(time, lake_state4, c='#006d2c', linewidth=2, linestyle='--') # intertemporal best benefits
    line5, = ax.plot(time, np.ones(len(time))*pCrit, c='#a50f15', linewidth=2) # critical threshold
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.8)
    ax.set_yticks(np.arange(0,1.0,0.2))
    ax.text(5,0.65,'c) Benefits-Maximizing\nLake P Trajectories', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
    
    ax = fig.add_subplot(2,3,4)
    line1, = ax.plot(time[1::], Y1, c='#08519c', linewidth=2) # DPS best reliability
    line2, = ax.plot(time[1::], Y2, c='#006d2c', linewidth=2) # DPS best benefits
    ax.set_ylabel('Anthropogenic P Release, $a_t$',fontsize=16)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.11)
    ax.text(5,0.09,'d) DPS Reliability & Benefits-\nMaximizing P Release Decisions',fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.2, box.width, box.height*0.9])
    
    ax = fig.add_subplot(2,3,5)
    line1, = ax.plot(time[1::], Y1, c='#08519c', linewidth=2) # DPS best reliability
    line3, = ax.plot(time[1::], IT[ITmostRel,0:100], c='#08519c', linewidth=2, linestyle='--') # intertemporal best reliability
    ax.set_xlabel('Year, $t$',fontsize=16)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.11)
    ax.text(5,0.09,'e) Reliability-Maximizing\nP Release Decisions', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.2, box.width, box.height*0.9])
    
    ax = fig.add_subplot(2,3,6)
    line2, = ax.plot(time[1::], Y2, c='#006d2c', linewidth=2) # DPS best benefits
    line4, = ax.plot(time[1::], IT[np.argmin(IT[:,101]),0:100], c='#006d2c', linewidth=2, linestyle='--') # intertemporal best benefits
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,100)
    ax.set_ylim(0,0.11)
    ax.text(5,0.09,'f) Benefits-Maximizing\nP Release Decisions', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.2, box.width, box.height*0.9])
    
    plt.figlegend([line5,line1,line2,line3,line4],['Critical P Threshold','DPS Highest Reliability',\
        'DPS Highest Benefits','Intertemporal Highest Reliability','Intertemporal Highest Benefits'], \
        loc = 'lower center', ncol=2)
    plt.suptitle('Lake P Dynamics',fontsize=18)
    fig.set_size_inches([19.825, 10.225])
    fig.savefig('Figure7.pdf')
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

def LakeModel_DPS(seed, vars):
    # Set inflow distribution parameters
    log_std = np.sqrt(np.log(1+sigma**2/mu**2))
    log_mu = np.log(mu) - 0.5*(log_std**2)
    
    # Initialize arrays to store P level in the lake at each time step
    lake_state = np.zeros([nYears+1])

    # Randomly generate nSamples of nYears of natural P inflows
    natFlow = np.zeros([nYears])
    np.random.seed(seed)
    natFlow= np.exp(ss.norm.rvs(log_mu, log_std, nYears))
        
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
    lake_state[0] = lake0
    Y = np.zeros([nYears])
    #find policy-derived emission
    Y[0] = RBFpolicy(lake_state[0], C, R, newW)

    for i in range(nYears):
        lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]**q)/(1+(lake_state[i]**q)) + Y[i] + natFlow[i]
        if i<(nYears-1):
            #find policy-derived emission
            Y[i+1] = RBFpolicy(lake_state[i+1], C, R, newW)
                
    return lake_state, Y

def RBFpolicy(lake_state, C, R, W):    
    # Determine pollution emission decision, Y
    Y = 0
    for i in range(len(C)):
        if R[i] != 0:
            Y = Y + W[i]*((np.absolute(lake_state-C[i])/R[i])**3)
            
    Y = min(0.1,max(Y,0.01))
    
    return Y
    
makeFigure7()