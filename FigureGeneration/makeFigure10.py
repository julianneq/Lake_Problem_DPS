import numpy as np
import matplotlib.pyplot as plt

def makeFigure10():
    # read in Latin hypercube samples of uncertain inputs
    LHsamples = np.loadtxt('./../Re-evaluation/LHsamples.txt')

    # read in satisficing measure of IT and DPS solutions
    IT = np.loadtxt('./../Re-evaluation/ITrobustness.txt',delimiter=' ')
    DPS = np.loadtxt('./../Re-evaluation/DPSrobustness.txt',delimiter=' ')

    # find most robust IT and DPS solutions
    ITbest = np.argmax(IT[:,2])
    DPSbest = np.argmax(DPS[:,2])

    # read in objective values of most robust IT and DPS solutions in alternative SOWs
    mostRobustIT = np.loadtxt('./../Re-evaluation/Intertemporal/output/ITobjs_' + str(int(ITbest)) + '.txt', delimiter=' ')
    mostRobustDPS = np.loadtxt('./../Re-evaluation/DPS/output/DPSobjs_' + str(int(DPSbest)) + '.txt', delimiter=' ')

    # determine in which SOWs the most robust IT and DPS solutions fails
    successes = [k for k in range(np.shape(LHsamples)[0]) if mostRobustIT[k,0] > 0.2 and mostRobustIT[k,3] > 0.95]
    failures = [k for k in range(np.shape(LHsamples)[0]) if mostRobustIT[k,0] <= 0.2 or mostRobustIT[k,3] <= 0.95]
    IT_success = LHsamples[successes,:]
    IT_fail = LHsamples[failures,:]
    
    successes = [k for k in range(np.shape(LHsamples)[0]) if mostRobustDPS[k,0] > 0.2 and mostRobustDPS[k,3] > 0.95]
    failures = [k for k in range(np.shape(LHsamples)[0]) if mostRobustDPS[k,0] <= 0.2 or mostRobustDPS[k,3] <= 0.95]
    DPS_success = LHsamples[successes,:]
    DPS_fail = LHsamples[failures,:]
    
    fig = plt.figure()
    ax = fig.add_subplot(2,2,1)
    successPts = ax.scatter(DPS_success[:,0],DPS_success[:,1],facecolor='#006d2c',edgecolor='none')
    failPts = ax.scatter(DPS_fail[:,0],DPS_fail[:,1],facecolor='0.25',edgecolor='none',alpha=0.5)
    ax.set_xlim(0.1,0.45)
    ax.set_xticks(np.arange(0.1,0.5,0.1))
    ax.set_yticks(np.arange(2.0,5.0,1.0))
    ax.set_ylim(2.0,4.5)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_ylabel('q',fontsize=16,rotation='horizontal')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.05, box.width, box.height*0.95])
    ax.set_title('a) DPS: Effect of b and q',loc='left',fontsize=12)
    
    ax = fig.add_subplot(2,2,2)
    ax.scatter(DPS_success[:,4],DPS_success[:,1],facecolor='#006d2c',edgecolor='none')
    ax.scatter(DPS_fail[:,4],DPS_fail[:,1],facecolor='0.25',edgecolor='none',alpha=0.5)
    ax.set_xlim(0.93,0.99)
    ax.set_xticks(np.arange(0.93,0.99,0.02))
    ax.set_yticks(np.arange(2.0,5.0,1.0))
    ax.set_ylim(2.0,4.5)
    ax.tick_params(axis='both',labelsize=14)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.05, box.width, box.height*0.95])
    ax.set_title('b) DPS: Effect of $\delta$',loc='left',fontsize=12)
    
    ax = fig.add_subplot(2,2,3)
    ax.scatter(IT_success[:,0],IT_success[:,1],facecolor='#006d2c',edgecolor='none')
    ax.scatter(IT_fail[:,0],IT_fail[:,1],facecolor='0.25',edgecolor='none',alpha=0.5)
    ax.set_xlim(0.1,0.45)
    ax.set_xticks(np.arange(0.1,0.5,0.1))
    ax.set_yticks(np.arange(2.0,5.0,1.0))
    ax.set_ylim(2.0,4.5)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel('b',fontsize=16)
    ax.set_ylabel('q',fontsize=16,rotation='horizontal')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.95])
    ax.set_title('a) Intertemporal: Effect of b and q',loc='left',fontsize=12)
    
    ax = fig.add_subplot(2,2,4)
    ax.scatter(IT_success[:,4],IT_success[:,1],facecolor='#006d2c',edgecolor='none')
    ax.scatter(IT_fail[:,4],IT_fail[:,1],facecolor='0.25',edgecolor='none',alpha=0.5)
    ax.set_xlim(0.93,0.99)
    ax.set_xticks(np.arange(0.93,0.99,0.02))
    ax.set_yticks(np.arange(2.0,5.0,1.0))
    ax.set_ylim(2.0,4.5)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel(r'$\delta$',fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.95])
    ax.set_title('a) Intertemporal: Effect of $\delta$',loc='left',fontsize=12)
    
    fig.suptitle('Parameter Combinations Leading to Failure',fontsize=16)
    plt.figlegend([successPts, failPts],['Meets Criteria','Fails to Meet Criteria'],loc='lower center',ncol=2)
    fig.set_size_inches([8.7375, 7.2])
    fig.savefig('Figure10.pdf')
    fig.clf()
    
    return None
    
makeFigure10()