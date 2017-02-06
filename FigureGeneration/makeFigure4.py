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

def makeFigure4():
    # read in reference sets and negate all objectives
    IT = -1*np.loadtxt('./../Optimization/Intertemporal.reference')
    DPS = -1*np.loadtxt('./../Optimization/DPS.reference')
    
    # find positive values of minimization objectives
    IT[:,1] = -IT[:,1]
    DPS[:,1] = -DPS[:,1]
    
    # find minimum and maximum values of each objective across the two sets
    objs_min = np.zeros(4)
    objs_max = np.zeros(4)
    for i in range(4):
        objs_min[i] = np.min([np.min(IT[:,i]), np.min(DPS[:,i])])
        objs_max[i] = np.max([np.max(IT[:,i]), np.max(DPS[:,i])])
        
    # figure with all points in color
    fig = plt.figure(figsize=plt.figaspect(0.5))
    
    ax = fig.add_subplot(1,2,1,projection='3d')
    IT_pts, DPS_pts, ideal_pt = makeSubPlot(IT, DPS, objs_min, objs_max, ax, fig)
    ax.view_init(elev=25.3, azim=42.7)
    ax.dist = 10

    # gray points except for best economic benefits and reliability
    ax = fig.add_subplot(1,2,2,projection='3d')
    IT_pts, DPS_pts, ideal_pt = BrushedPlot(IT, DPS, objs_min, objs_max, ax, fig)
    ax.view_init(elev=25.3, azim=42.7)
    ax.dist = 10
        
    l1 = ax.scatter([],[],[], s=20, color='k', linewidth=0, depthshade=False)
    l2 = ax.scatter([],[],[], s=80, color='k', linewidth=0, depthshade=False)
        
    legend1 = fig.legend([IT_pts, DPS_pts, ideal_pt],["Intertemporal", "DPS", "Ideal Point"], \
        loc='lower left',scatterpoints=1, title='Solution Strategy')
    plt.setp(legend1.get_title(), fontsize=14)
    ax = plt.gca().add_artist(legend1)
    fig.legend([l1, l2],[str(int(100*objs_min[3])) + '%',str(int(100*objs_max[3])) + '%'], \
        scatterpoints=1, title="Reliability", fontsize=14, loc='lower right')
    fig.set_size_inches([12,5.975])
    fig.savefig('Figure4.pdf')
    fig.clf()
    
    return None
    
def makeSubPlot(IT, DPS, objs_min, objs_max, ax, fig):
    pts_IT = ax.scatter(IT[:,2], IT[:,1], IT[:,0], \
        s=20+60*(IT[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), edgecolor='#a50f15', \
        facecolor='none', linewidth=1, depthshade=False)
    pts_DPS = ax.scatter(DPS[:,2], DPS[:,1], DPS[:,0], \
        s=20+60*(DPS[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), color='#08519c', \
        linewidth=0, depthshade=False)

    pt_ideal = ax.scatter(1.0, 0.0, objs_max[0], c='black', s=500, linewidth=0, marker='*',)
    ax.set_xticks(np.arange(0.92, 1.02, 0.02))
    ax.set_yticks(np.arange(0.0, 2.5, 0.5))
    ax.set_zticks(np.arange(0.25, 0.6, 0.05))
    ax.set_xlabel("Inertia")
    ax.set_ylabel("P Concentration")
    ax.set_zlabel("Economic Benefit")
    ax.set_xlim(objs_min[2],1.0)
    ax.set_ylim(0.0,objs_max[1])
    ax.set_zlim(objs_min[0],objs_max[0])
    
    return pts_IT, pts_DPS, pt_ideal
    
def BrushedPlot(IT, DPS, objs_min, objs_max, ax, fig):
    ax.scatter(IT[:,2], IT[:,1], IT[:,0], \
        s=20+60*(IT[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), edgecolor='#969696', \
        facecolor='none', linewidth=1, depthshade=False, alpha=0.3)
    ax.scatter(DPS[:,2], DPS[:,1], DPS[:,0], \
        s=20+60*(DPS[:,3]-objs_min[3])/(objs_max[3]-objs_min[3]), color='#969696', \
        linewidth=0, depthshade=False, alpha=0.3)
    IT_best_rel = np.argmax(IT[:,3]) # one of several solutions with perfect reliability; in paper, row 8 of IT.resultfile was plotted
    IT_best_ben = np.argmax(IT[:,0])
    pt_IT_rel = ax.scatter(IT[IT_best_rel,2], IT[IT_best_rel,1], IT[IT_best_rel,0], \
        s=20+60*(IT[IT_best_rel,3]-objs_min[3])/(objs_max[3]-objs_min[3]), edgecolor="#a50f15", \
        facecolor='none', linewidth=2, depthshade=False)
    ax.scatter(IT[IT_best_ben,2], IT[IT_best_ben,1], IT[IT_best_ben,0], \
        s=20+60*(IT[IT_best_ben,3]-objs_min[3])/(objs_max[3]-objs_min[3]), edgecolor="#a50f15", \
        facecolor='none', linewidth=2, depthshade=False)
    DPS_best_rel = np.argmax(DPS[:,3]) # one of several solutions with perfect reliability; in paper, row 26 of DPS.resultfile was plotted
    DPS_best_ben = np.argmax(DPS[:,0])
    pt_DPS_rel = ax.scatter(DPS[DPS_best_rel,2], DPS[DPS_best_rel,1], DPS[DPS_best_rel,0], \
        s=20+60*(DPS[DPS_best_rel,3]-objs_min[3])/(objs_max[3]-objs_min[3]), color="#08519c", \
        linewidth=0, depthshade=False)
    ax.scatter(DPS[DPS_best_ben,2], DPS[DPS_best_ben,1], DPS[DPS_best_ben,0], \
        s=20+60*(DPS[DPS_best_ben,3]-objs_min[3])/(objs_max[3]-objs_min[3]), color="#08519c", \
        linewidth=0, depthshade=False)

    pt_ideal = ax.scatter(1.0, 0.0, objs_max[0], c='black', s=500, linewidth=0, marker='*',)
    ax.set_xticks(np.arange(0.92, 1.02, 0.02))
    ax.set_yticks(np.arange(0.0, 2.5, 0.5))
    ax.set_zticks(np.arange(0.25, 0.6, 0.05))
    ax.set_xlabel("Inertia")
    ax.set_ylabel("P Concentration")
    ax.set_zlabel("Economic Benefit")
    ax.set_xlim(objs_min[2],1.0)
    ax.set_ylim(0.0,objs_max[1])
    ax.set_zlim(objs_min[0],objs_max[0])
    
    return pt_IT_rel, pt_DPS_rel, pt_ideal
    
makeFigure4()