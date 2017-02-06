import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

def makeFigure11():
    # read in Latin hypercube samples of uncertain inputs
    LHsamples = np.loadtxt('./../Re-evaluation/LHsamples.txt')

    # read in satisficing measure of DPS solutions
    DPS = np.loadtxt('./../Re-evaluation/DPSrobustness.txt',delimiter=' ')

    # find most robust DPS solution
    DPSbest = np.argmax(DPS[:,2])

    # read in objective values of most robust DPS solution in alternative SOWs
    mostRobustDPS = np.loadtxt('./../Re-evaluation/DPS/output/DPSobjs_' + str(int(DPSbest)) + '.txt', delimiter=' ')

    # determine in which SOWs the most robust DPS solution fails
    successes = [k for k in range(np.shape(LHsamples)[0]) if mostRobustDPS[k,0] > 0.2 and mostRobustDPS[k,3] > 0.95]
    failures = [k for k in range(np.shape(LHsamples)[0]) if mostRobustDPS[k,0] <= 0.2 or mostRobustDPS[k,3] <=0.95]
    DPS_success = LHsamples[successes,:]
    DPS_fail = LHsamples[failures,:]
    
    def fun(x):
        return [(x[0]**q)/(1+x[0]**q) - b*x[0]]
    
    B = np.arange(0.1,0.46,(0.46-0.1)/50)
    Q = np.arange(2.0,4.6,(4.6-2.0)/50)
    X, Y = np.meshgrid(B,Q)
    Z = np.zeros(np.shape(X))
    for i in range(len(B)):
        for j in range(len(Q)):
            b = X[i,j]
            q = Y[i,j]
            soln = root(fun, 0.75)
            Z[i,j] = soln.x
            
    levels = [0,0.25,0.5,0.75,1.0]
    plt.contourf(X, Y, Z, levels, cmap='RdYlGn')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)
    plt.scatter(DPS_success[:,0],DPS_success[:,1],facecolor='k',edgecolor='none')
    plt.scatter(DPS_fail[:,0],DPS_fail[:,1],facecolor='0.25',edgecolor='none',alpha=0.5)
    plt.xlim(0.1,0.45)
    plt.ylim(2.0,4.5)
    plt.xticks(np.arange(0.1,0.5,0.1))
    plt.yticks(np.arange(2.0,5.0,1.0))
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('b',fontsize=16)
    plt.ylabel('q',fontsize=16,rotation='horizontal')
    plt.title('Critical P Threshold',fontsize=16)
    plt.savefig('Figure11.pdf')
    plt.clf()
    
    return None
    
makeFigure11()