import math
import os
import numpy as np

def reformatData(method, abbrev, nSamples, nObjs):
    numPts = (len(os.listdir('./' + method + '/output')))
    objs = np.zeros([numPts, nSamples, nObjs])
    
    for i in range(numPts):
        objs[i,:,:] = np.loadtxt('./' + method + '/output/' + abbrev + 'objs_' + str(int(i)) + '.txt')

    return objs

def calcSatisfaction(Objs):
    '''Calculates the percent of SOWs in which a specified Satisficing Criteria\n
    is met.'''
    SatisfyMatrix = np.zeros([np.shape(Objs)[0],np.shape(Objs)[1],3])
    for i in range(np.shape(SatisfyMatrix)[0]):
        for j in range(np.shape(SatisfyMatrix)[1]):
            if Objs[i,j,0] > 0.2:
                SatisfyMatrix[i,j,0] = 1
            if Objs[i,j,3] > 0.95:
                SatisfyMatrix[i,j,1] = 1
            if Objs[i,j,0] > 0.2 and Objs[i,j,3] > 0.95:
                SatisfyMatrix[i,j,2] = 1
                
    satisfaction = np.zeros([np.shape(SatisfyMatrix)[0],np.shape(SatisfyMatrix)[2]])
    for i in range(np.shape(SatisfyMatrix)[0]):
        for j in range(np.shape(SatisfyMatrix)[2]):
            satisfaction[i,j] = np.mean(SatisfyMatrix[i,:,j])
    
    return satisfaction

nSamples=1000
nObjs=4
directory = './output'

methods = ['DPS','Intertemporal']
abbrevs=['DPS','IT']
for i in range(len(methods)):
	objs = reformatData(methods[i], abbrevs[i], nSamples, nObjs)
	satisfaction = calcSatisfaction(objs)
	np.savetxt(abbrevs[i] + 'robustness.txt', satisfaction, delimiter=' ')

