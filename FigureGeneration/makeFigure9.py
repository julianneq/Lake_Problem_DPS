'''
Created on Jul 11, 2014
@author: Victoria Lynn Ward
vlw27@cornell.edu
Cornell University
python script for a 4 objective 3d scatter plot

Modified by Julianne Quinn
June 16, 2015
'''
from matplotlib import pyplot as plt
#from matplotlib.backends import backend_agg as agg #raster backend
from mpl_toolkits.mplot3d import Axes3D
#import pandas # data analysis library
import numpy as np

def makeFigure9():
    # read in reference sets and negate all objectives
    IT = -1*np.loadtxt('./../Optimization/Intertemporal.reference')
    DPS = -1*np.loadtxt('./../Optimization/DPS.reference')
    
    # find positive values of minimization objectives
    IT[:,1] = -IT[:,1]
    DPS[:,1] = -DPS[:,1]
    
    #read in robustness
    ITrobustness = np.loadtxt('./../Re-evaluation/ITrobustness.txt',delimiter=' ')
    ITrobustness = ITrobustness*100
    DPSrobustness = np.loadtxt('./../Re-evaluation/DPSrobustness.txt',delimiter=' ')
    DPSrobustness = DPSrobustness*100
    
    objs_min = np.zeros(4)
    objs_max = np.zeros(4)
    for i in range(4):
        objs_min[i] = np.min([np.min(IT[:,i]), np.min(DPS[:,i])])
        objs_max[i] = np.max([np.max(IT[:,i]), np.max(DPS[:,i])])
        
    fig = plt.figure(figsize=plt.figaspect(0.5))
    
    ax = fig.add_subplot(1,1,1,projection='3d')
    IT_pts, DPS_pts, ideal_pt = makeSubPlot(IT, ITrobustness[:,2], DPS, \
            DPSrobustness[:,2], objs_min, objs_max, ax, fig)
    ax.view_init(elev=25.3, azim=42.7)
    ax.dist = 10
    fig.colorbar(IT_pts, ax=ax, shrink=0.75)
    fig.colorbar(DPS_pts, ax=ax, shrink=0.75)
    fig.axes[-1].set_ylabel("Percent of SOWs in which\nReliability > 95% and\nEconomic Benefits > 0.1")
        
    l1 = ax.scatter([],[],[], s=20, color='k', linewidth=0, depthshade=False)
    l2 = ax.scatter([],[],[], s=80, color='k', linewidth=0, depthshade=False)
        
    legend1 = fig.legend([IT_pts, DPS_pts, ideal_pt],["Intertemporal", "DPS", "Ideal Point"], \
        loc='lower left',scatterpoints=1)
    plt.setp(legend1.get_title(), fontsize=14)
    ax = plt.gca().add_artist(legend1)
    fig.legend([l1, l2],[str(round(objs_min[3],2)),str(round(objs_max[3],2))], scatterpoints=1, \
        title="Reliability", fontsize=14, loc='lower right')
    fig.set_size_inches([10.2125, 5.975])
    fig.show()
    fig.savefig('Figure9.pdf')
    fig.clf()
    
    return None
    
def makeSubPlot(IT, ITrobustness, DPS, DPSrobustness, objs_min, objs_max, ax, fig):
    pts_IT = ax.scatter(IT[:,2], IT[:,1], IT[:,0], \
        s=20+60*(IT[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), c=ITrobustness, \
        cmap=plt.cm.get_cmap("Reds"), facecolors='none', linewidth=1, depthshade=False)
    IT_index = np.argmax(ITrobustness)
    cmap = plt.get_cmap("Reds")
    ax.scatter(IT[IT_index,2], IT[IT_index,1], IT[IT_index,0],\
        s=800, c=cmap(1.0), facecolors='none', linewidth=1, depthshade=False)
    pts_DPS = ax.scatter(DPS[:,2], DPS[:,1], DPS[:,0], \
        s=20+60*(DPS[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), c=DPSrobustness, \
        cmap=plt.cm.get_cmap('Blues'), linewidth=0, depthshade=False)
    DPS_index = np.argmax(DPSrobustness)
    cmap = plt.get_cmap("Blues")
    ax.scatter(DPS[DPS_index,2], DPS[DPS_index,1], DPS[DPS_index,0],\
        s=800, c=cmap(1.0), linewidth=0, depthshade=False)
    

    pt_ideal = ax.scatter(1.0, 0.0, objs_max[0], c='black', s=500, linewidth=0, marker='*',)
    ax.set_xlabel("Inertia")
    ax.set_ylabel("P Concentration")
    ax.set_zlabel("Economic Benefit")
    ax.set_xlim(objs_min[2],1.0)
    ax.set_ylim(0.0,objs_max[1])
    ax.set_zlim(objs_min[0],objs_max[0])
    
    return pts_IT, pts_DPS, pt_ideal
    
makeFigure9()