import numpy as np
import matplotlib.pyplot as plt

def makeFigure8():
    IT = 100*np.loadtxt('./../Re-evaluation/ITrobustness.txt',delimiter=' ')
    DPS = 100*np.loadtxt('./../Re-evaluation/DPSrobustness.txt',delimiter=' ')
    
    titles = ['a) Economic Benefits > 0.2','b) Reliability > 95%','c) Economic Benefits > 0.2 & Reliability > 95%']
    p1 = plt.Rectangle((0, 0), 1, 1, fc='#08519c', edgecolor='none') # DPS color
    p2 = plt.Rectangle((0, 0), 1, 1, fc='#a50f15', edgecolor='none') # intertemporal color
    multiplier = [0.06, 0.03, 0.0]
    
    fig = plt.figure()
    for i in range(len(titles)):
        ax = fig.add_subplot(3,1,i+1)
        ax.plot(range(np.shape(DPS)[0]+1),np.append(np.sort(DPS[:,i])[::-1],0),color='#08519c', linewidth=2)
        ax.plot(range(np.shape(IT)[0]+1),np.append(np.sort(IT[:,i])[::-1],0),color='#a50f15', linewidth=2)
        ax.fill_between(range(np.shape(DPS)[0]+1),np.append(np.sort(DPS[:,i])[::-1],0),color='#08519c')
        ax.fill_between(range(np.shape(IT)[0]+1),np.append(np.sort(IT[:,i])[::-1],0),color='#a50f15')
        ax.tick_params(axis='both',labelsize=14)
        ax.set_xlim([0,np.shape(DPS)[0]+1])
        ax.set_ylim([0,100])
        ax.set_title(titles[i],fontsize=16,loc='left')
        box = ax.get_position()
        ax.set_position([box.x0,box.y0+box.height*multiplier[i], box.width, box.height*0.97])
        if i == 2:
            ax.set_xlabel('Solution # (sorted by rank)',fontsize=16)
            
    fig.text(0.02, 0.5, 'Percent of Sampled SOWs in which Criteria are Met', va='center', rotation='vertical',fontsize=14)
    plt.figlegend([p1,p2],['DPS','Intertemporal'], loc='upper center', ncol=2)
    fig.set_size_inches([6.1625, 12.35])
    fig.savefig('Figure8.pdf')
    fig.clf()
    
    return None
    
makeFigure8()
