import numpy as np
from scipy.optimize import root
from matplotlib import pyplot as plt

n = 2
q = 2
b = 0.42
def fun(x):
    return [(x[0]**q)/(1+x[0]**q) - b*x[0]]
soln = root(fun, 0.75)
pCrit = soln.x

def makeFigure2():
    DPS = np.loadtxt('./../Optimization/DPS.resultfile')
    DPSmostRel = np.argmin([DPS[:,9]]) # one of several solutions with perfect reliability; in paper, row 26 of DPS.resultfile was plotted
    lake_state = np.arange(0.0,2.01,0.01)
    Y = np.zeros([len(lake_state),2])
    
    for j in range(len(lake_state)):
        Y[j,0] = RBFpolicy([lake_state[j]], DPS[DPSmostRel,0:6]) # best reliability
        
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(lake_state,Y[:,0],c='#08519c', linewidth=2, label='Sample P Release Policy')
    ax.plot([pCrit, pCrit],[0.0,0.1],c='#a50f15', linewidth=2, label='Critical P Threshold')
    ax.set_xlabel('Lake P Concentration, $X_t$',fontsize=16)
    ax.set_ylabel('Anthropogenic P Release, $a_t$',fontsize=16)
    ax.set_xlim(0,1)
    ax.set_ylim(0,0.15)
    ax.tick_params(axis='both',labelsize=14)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
    ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.2), fancybox=True)
    fig.set_size_inches([6.2, 5.6])
    fig.savefig('Figure2.pdf')
    fig.clf()

    return None

def RBFpolicy(inputs, vars):
    # Determine centers, radii and weights of RBFs
    m = len(inputs)
    C = np.zeros([n,m])
    B = np.zeros([n,m])
    W = np.zeros(n)
    newW = np.zeros(len(W))

    for i in range(n):
        W[i] = vars[(i+1)*(2*m+1)-1]
        for j in range(m):
            C[i,j] = vars[(2*m+1)*i + 2*j]
            B[i,j] = vars[(2*m+1)*i + 2*j + 1]
    
    # Normalize weights to sum to 1
    total = sum(W)
    if total != 0.0:
        for i in range(len(W)):
            newW[i] = W[i]/total
    else:
        for i in range(len(W)):
            newW[i] = 1/n
            
    Y = 0
    for i in range(n):
        for j in range(m):
            if B[i,j] != 0:
                Y = Y + W[i]*((np.absolute(inputs[j]-C[i,j])/B[i,j])**3)

    Y = min(0.1,max(Y,0.01))
    
    return Y
    
makeFigure2()