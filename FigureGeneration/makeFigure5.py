import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os

def makeFigure5():

    subfolders = ['./../Optimization/DPS','./../Optimization/Intertemporal']
    colors = ['#08519c','#a50f15']
    steps = 200
    
    #normalized between 0 and 1
    
    All_DPS_HV = np.zeros([50,steps+1])
    All_Trad_HV = np.zeros([50,steps+1])
    NFE = np.arange(0,steps+1,1)*1000
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i in range(len(subfolders)):
        if i == 1:
            metricsList = [f for f in os.listdir( subfolders[i] + "/metrics") if f[0:7] == 'LakeIT_']
        else:
            metricsList = [f for f in os.listdir(subfolders[i] + "/metrics") if f[0:8] == 'LakeDPS_']
        
        for j in range(1,len(metricsList)+1):
            if i == 1:
                HV = np.loadtxt(subfolders[i] + '/metrics/LakeIT_S%d.metrics' % j, skiprows=1, usecols=[0])
            else:
                HV = np.loadtxt(subfolders[i] + '/metrics/LakeDPS_S%d.metrics' % j, skiprows=1, usecols=[0])
            
            if len(HV) < steps+1:
                HV = np.concatenate((np.zeros(steps+1-len(HV)),HV),0)
                
            if i == 1:
                All_Trad_HV[j-1,:] = HV
            else:
                All_DPS_HV[j-1,:] = HV
                
    maxHV = np.max([np.max(All_Trad_HV[:,-1]), np.max(All_DPS_HV[:,-1])])
    All_Trad_HV = (1/maxHV)*All_Trad_HV
    All_DPS_HV = (1/maxHV)*All_DPS_HV
    
    for j in range(np.shape(All_DPS_HV)[0]):
        line1, = ax.plot(NFE, All_DPS_HV[j,:], color=colors[0], linewidth=2)
        line2, = ax.plot(NFE, All_Trad_HV[j,:], color=colors[1], linewidth=2)
            
    ax.set_ylabel('Relative Hypervolume Performance', fontsize=16)
    ax.set_xlabel('NFE', fontsize=16)
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.set_xlim(0,200000)
    ax.set_ylim(0,1.05)
    ax.set_yticks(np.arange(0,1.2,0.2))
    ax.tick_params(axis='both',labelsize=14)
    legend = plt.legend([line1, line2],['DPS','Intertemporal'],loc='lower right',title='Solution Strategy',fontsize=16)#bbox_to_anchor=(0.72,1.04,1.0,0.1),loc=2,borderaxespad=0.)
    plt.setp(legend.get_title(),fontsize=16)
    fig.set_size_inches([8, 5.975])
    plt.savefig('Figure5.pdf')
    plt.clf()
    
    return None
    
makeFigure5()