import matplotlib.pyplot as plt
from scipy.optimize import root
import matplotlib
import numpy as np

def makeFigure1():
    def fun(x):
        return [(x[0]**qq)/(1+x[0]**qq) - bb*x[0]]
        
    b = [0.4,0.3,0.2,0.1]
    q = [2.5,3,3.5,4]
    
    x = np.arange(0,2.6,0.1)
    y = np.zeros(len(x))
    
    cmap = matplotlib.cm.get_cmap('RdYlGn')
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    # for b = 0.4, plot different values of q
    colors = []
    for i in range(len(b)):
        colors.append(cmap(0.25*(i+1)))
        
    lines=[]
    line1, = ax.plot(x,b[0]*x,c='k', linewidth=2)
    lines.append(line1)
    for j in range(len(q)):
        for i in range(len(x)):
            y[i] = x[i]**q[j]/(1+x[i]**q[j])
        
        line1, = ax.plot(x, y, c = colors[j], linewidth=2)
        lines.append(line1)
        bb = b[0]
        qq = q[j]
        soln = root(fun,1.0)
        lines.append(ax.scatter(soln.x,b[0]*soln.x,facecolor='none',edgecolor='k',s=30))
        soln = root(fun,2.5)
        lines.append(ax.scatter(soln.x,b[0]*soln.x,facecolor='k',edgecolor='k',s=30))
        
    lines.append(ax.scatter(0,0,facecolor='k',edgecolor='k',s=30))
    legend1 = plt.legend([lines[0], lines[1], lines[4], lines[7], lines[10]],\
        ['b = 0.4', 'q = 2.5', 'q = 3.0', 'q = 3.5', 'q = 4.0'], loc='lower right')
    plt.setp(legend1.get_title(),fontsize=14)
    plt.gca().add_artist(legend1)
    plt.legend([lines[3], lines[2]],['Stable Equilibria','Unstable Equilibria'],loc='upper left')
    ax.set_ylabel('Fluxes of P',fontsize=16)
    ax.tick_params(axis='y',labelsize=14)
    ax.set_xlim(0,2.5)
    ax.set_ylim(0,1)
    ax.set_title('a) Effect of q on Lake Dynamics',loc='left')
    
    ax = fig.add_subplot(2,1,2)
    colors = []
    for i in range(len(b)):
        colors.append(cmap(1-(0.25*i)))
        
    #for q = 2.5, plot different values of b
    for i in range(len(x)):
        y[i] = x[i]**q[0]/(1+x[i]**q[0])
            
    lines = []
    line1, = ax.plot(x,y,c='k',label='q = ' + str(q[0]),linewidth=2)
    lines.append(line1)
    for i in range(len(b)):        
        line1, = ax.plot(x,b[i]*x,c=colors[i],label='b = ' + str(b[i]),linewidth=2)
        lines.append(line1)
        bb = b[i]
        qq = q[0]
        soln = root(fun,1.0)
        lines.append(ax.scatter(soln.x,b[i]*soln.x,facecolor='none',edgecolor='k',s=30))
        soln = root(fun,2.5)
        lines.append(ax.scatter(soln.x,b[i]*soln.x,facecolor='k',edgecolor='k',s=30))
        
    lines.append(ax.scatter(0,0,facecolor='k',edgecolor='k',s=30))
    ax.legend([lines[0], lines[1], lines[4], lines[7], lines[10]],\
        ['q = 2.5', 'b = 0.4', 'b = 0.3', 'b = 0.2', 'b = 0.1'],\
        scatterpoints = 1, loc='upper left')
    ax.set_xlabel('Lake P Concentration,$X_t$',fontsize=16)
    ax.set_ylabel('Fluxes of P',fontsize=16)
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlim(0,2.5)
    ax.set_ylim(0,1)
    ax.set_title('b) Effect of b on Lake Dynamics',loc='left')
    fig.set_size_inches([8,11.85])
    fig.savefig('Figure1.pdf')
    fig.clf()
    
    return None
    
makeFigure1()