# -*- coding: utf-8 -*-
#pylint: disable=pointless-string-statement

from __future__ import division
"""
    Launcher, entry point for the pleiotropy experiments
"""
#from ea.genome import Genome
#TODO GITify everything. Create repository on github

from pybrain.structure import FeedForwardNetwork, LinearLayer
from pybrain.structure import FullConnection, SigmoidLayer, BiasUnit
#from pybrain.structure import RecurrentNetwork
import numpy as np
#from pybrain.structure.evolvables.evolvable import Evolvable
#from numpy.random import multivariate_normal
#TODO lots of sanity checks, such as tournSize < popSize
#TODO Final integration test. Migrate 

###############################################################################
####################################Plotting###################################
#import os, re, sys
#from os import listdir
#from os.path import isfile, join
#from scipy.stats import *
#from pylab import *
import brewer2mpl


bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors   

###############################################################################
###############################################################################
def perc(data_l):    
    data = np.asarray(data_l)
    median = np.zeros(data.shape[1])
    perc_25 = np.zeros(data.shape[1])
    perc_75 = np.zeros(data.shape[1])
    for i in range(0, len(median)):
        median[i] = np.median(data[:, i])
        perc_25[i] = np.percentile(data[:, i], 25)
        perc_75[i] = np.percentile(data[:, i], 75)
    return median, perc_25, perc_75


def plot_one_curve(data, color, axis, label, quartiles=False):
    med, perc_25, perc_75 =  perc(data)
    if quartiles :
        axis.fill_between(np.arange(0, len(med)), perc_25, perc_75,
                          alpha=0.25, linewidth=0, color=color)
    axis.plot(med, lw=2, label=label, color=color)
    
      
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    for spine in axis.spines.values():
        spine.set_position(('outward', 5))
    axis.set_axisbelow(True)

def fullConnectionNet(nbIn, nbOut, nbHid, neurPerHid, bias):
    
    net = FeedForwardNetwork()
    
    for i in xrange(nbIn):
        net.addInputModule(LinearLayer(1,name="i" + str(i)))        
    
    for i in xrange(nbOut):
        net.addOutputModule(SigmoidLayer(1,name="o" + str(i)))
    
    for i in xrange(nbHid):
        for j in xrange(neurPerHid):
            net.addModule(
            SigmoidLayer(1,
                         name="h(" + str(i) + "-" + str(j) + ")"))

    if nbHid == 0 :        
        for i in xrange(nbIn):
            for j in xrange(nbOut):
                net.addConnection(
                FullConnection(
                net["i" + str(i)],
                net["o" + str(j)],
                name = "i" + str(i) + "-o" + str(j)))
    else:
        for i in xrange(nbIn):
            for j in xrange(neurPerHid):                
                net.addConnection(
                FullConnection(
                net["i" + str(i)],
                net["h(0-" + str(j) + ")"] ,
                name = "i" + str(i) + "-h(0-" + str(j) + ")"))
        
        for k in xrange(nbHid - 1):
            for i in xrange(neurPerHid):
                for j in xrange(neurPerHid):    
                    net.addConnection(
                    FullConnection(
                    net["h(" + str(k) + "-" + str(i) + ")"], 
                    net["h(" + str(k + 1) + "-" + str(j) + ")"],
                    name = "h(" + str(k) + "-" + str(i) + ")-h(" 
                        + str(k + 1) + "-" + str(j) + ")"))
            
        for i in xrange(neurPerHid):
            for j in xrange(nbOut):
                net.addConnection(
                FullConnection(
                    net["h(" + str(nbHid - 1) + "-" + str(i) + ")"],
                    net["o" + str(j)],
                    name = "h(" + str(nbHid - 1) + "-" + str(i) + ")-o" + str(j)))
    if bias:
        b = BiasUnit(name = "b")
        net.addModule(b)
        for i in xrange(nbHid):
            for j in xrange(neurPerHid):
                net.addConnection(
                FullConnection(b, 
                               net["h(" + str(i) + "-" + str(j) + ")"],
                               name = "b-h(" + str(i) + "-" + str(j) + ")"))
        for i in xrange(nbOut):
            net.addConnection(
            FullConnection(b,
                           net["o" + str(i)],
                           name = "b-o" + str(i)))
    
    net.sortModules()    
    
    return net

###############################################################################
    
def evaluate(indiv, pbDb, labelDb):
    fit = 0
    for i in xrange(len(pbDb)):
        answer = round(indiv.activate(pbDb[i]),0)
        
        fit += (answer == labelDb[i])    
    #Balance fitness w.r.t. proportion of same-valued answers:
    #e.g. 99 True, 1 False
    return fit/float(len(pbDb)) # accuracy

#TODO Check all selection methods
def select(p, f, selM, nbIndiv, tSize, elit):
    #Pass into global.py or config.py
    selMethods = {"RW": roulette,
              "RB": rankBased,
              "Tourn": tournament,
              "Best": best,
              }
    #TODO try factorised version like: (sortedPop, sortedF) = 
              #[(x, y) for (y,x) in sorted(zip(f,p), key=lambda pair: pair[0])] 
    sortedPop= [x 
        for (y,x) in sorted(zip(f,p), key=lambda pair: pair[0], reverse=True)]
    sortedF = [y 
        for (y,x) in sorted(zip(f,p), key=lambda pair: pair[0], reverse=True)]  
    
    parents = []

    if elit:
        parents.append(sortedPop[0])
        sortedPop.remove(sortedPop[0])
    
    parents.extend(selMethods[selM](sortedPop, sortedF, nbIndiv, tSize)) 
           
    return parents

def roulette(sPop, sF, nb, tSz):
    result = []
    total =  sum(sF)
    prob =  sF / total
    print prob
    
    n=0
    #TODO debug random.choice for roulette
    while n < nb:
        result.append(np.random.choice(sPop, replace=True, p = prob))
    return result

def rankBased(sPop, sF, nb, tSz):
    result = []
    total = len(sPop)*(len(sPop) + 1) // 2 #n(n + 1)/2 
    
    rankWeights = range(1, len(sPop) + 1)
    rankWeights = list(reversed(rankWeights))
    
    #Another version to use with random.choice()
    #prob =  [x / total for x in rankWeights]
    
    n = 0
    while n < nb:
        dice = np.random.random()* total
        index = 1
        while dice > 0 and index <= len(sPop):
            dice -= index 
            index += 1
        result.append(sPop[len(sPop) - (index - 1)])
        n += 1
    
    return result

def tournament(sPop, sF, nb, tSz):
    result = []
    #TODO Tournament
    for i in xrange(nb):
        #No replacement when sampling for tournament        
        print sorted(np.random.sample(xrange(len(sPop)), tSz, replace = False))
        #[sPop[i] 
        #for i in 
        #sorted(np.random.sample(xrange(len(sPop)), tSz, replace = False))][0]
        #result.append()
        
    return result

def best(sPop, sF, nb, tSz):
    #TODO test
    result = []
    result.extend(sPop[0 : nb - 1]) 
    return result


def gnrtOffspring(par, pSz, pSurv, sig):
    children = []    

    if pSurv:
        #TODO test
        children.extend(par)

    i = 0        
    while len(children) < pSz:
        children.append(mutate(par[i], sig))
        #Loop through the parents in order. 
        #This could be randomized or done otherwise
        i = (i + 1) % len(par) 
        
    return children
    
def mutate(indiv, sigma):
    result = indiv.copy()

    for i in xrange(result.paramdim):
        result.params[i] = np.random.normal(result.params[i], sigma)
    
    return result

def task(vector,taskID):
    #TODO Pass into global.py or config.py
    tasks= {"AND": all,
              "OR": any,
              "MAJ": majority,
              "MIN": minority,
              "MUX1": mux1,
              }
    return tasks[taskID](vector)
    
def majority(x):
    total = sum(x)
    if total >= len(x):
        return True
    else:
        return False

def minority(x):
    return not majority(x)
    
def mux1(x): 
    '''Multiplexer type 1: 
        if x[0] is False, then return x[1:ceil((len(x)-1)/2)]
        else if x[0] is True, then returnx[ceil((len(x)-1)/2) + 1: len(x) - 1]
    '''
    if x[0] == False:
        return x[1:np.ceil((len(x) - 1)/2)]
    else:
        return x[np.ceil((len(x)-1)/2) + 1: len(x) - 1]            

###############################################################################

percepts = 6
bias = False
out = 1
nbHLayers = 2
neurPerLayer = 8

###############################################################################

pbDb = []
labelDb = []
pbSz = percepts
dbSz = 200
problemID = "MAJ" # "MIN", "AND", "OR", "MAJ". 
#TODO TOFIX: "MUX1" (distance definitions in N-dimensional spaces needed)

for i in xrange(dbSz):
    v = np.random.rand(pbSz)
    for i in xrange(pbSz):
        # Or random.choice( [0,1] )
        v[i] = round(v[i], 0) 
    pbDb.append(v)
    #Different tasks. 
    #TODO Refactor into independent "task" module    
    #labelDb.append(all(v)) #and function
    #Learn random associations (impossible task)
    #labelDb.append(round(np.random.rand(),0)) 
    
    labelDb.append(task(v, problemID))


###############################################################################

#Evo algo
sigma = 0.5
popSz = 4

maxGen = 20
#TODO Implemented: RB, Best; TOTEST/FIX: RW, Tourn (N)
selection = "RW" 
N = 2 # Tourn. size if Tourn(N) selected

#Choose the best as parent #Maybe just make it survive as it is
elitist = False 
# Keep the parents as they are for the next population, 
#then generated the remaining popSz - |parents| individuals to fill the pop.
parSurv = False 

#Percentage of individuals in population that are selected as parents
#Other selection strategies can be considered
selFactor = 0.4 
nbPar = int(np.ceil(selFactor * popSz))

###############################################################################
#TODO Other methods for pop init. 
#Implemented: fully connected MLP, bias allowed
pop = []
for i in xrange(popSz):
    pop.append(
        fullConnectionNet( percepts, out, nbHLayers, neurPerLayer, bias))

    
fitness = [[0]*popSz for i in range(maxGen + 1)]

#TODO TOCHECK termination conditions and data structures boundaries
###############################################################################

#Start Evo loop
for g in xrange(maxGen):   
    for ind in xrange(popSz):
        fitness[g][ind] = evaluate(pop[ind], pbDb, labelDb)
    print "Generation: " + str(g)
    #ToFinish: parent selection, offspring generation
    parents = select(pop, fitness[g], selection, nbPar, N, elitist)
    
    offspring = gnrtOffspring(parents, popSz, parSurv, sigma)    
    pop = offspring
            

#Last generation eval
for ind in xrange(popSz):    
    fitness[maxGen][ind] = evaluate(pop[ind], pbDb, labelDb)

print fitness

axis = subplot2grid((1,1), (0,0))

plot_one_curve(zip(*fitness), colors[0], axis, "Fitness", True)

print "Max final fitness: " + str(np.max(fitness[maxGen])) + ", "
print "Mean final fitness: " + str(np.mean(fitness[maxGen])) + ", "
#print "Reachable fitness with systematic \"True\" answer: " 
        #+ str((2**pbSz-1)/2**pbSz)

###############################################################################

"""
g1 = Genome(9,2,0,0,100)

g2 = Genome(9,2,1,5,100)

g3 = Genome(9,2,2,5,100)

print g1.getMaxConnections()

print g2.getMaxConnections()

print g3.getMaxConnections()

n. """



































'''
means = [(-1,0),(2,4),(3,1)]
cov = [np.diag([1,1]), np.diag([0.5,1.2]), np.diag([1.5,0.7])]
alldata = ClassificationDataSet(2, 1, nb_classes=3)

for n in xrange(400):
    for klass in range(3):
        input = multivariate_normal(means[klass],cov[klass])
        alldata.addSample(input, [klass])

tstdata, trndata = alldata.splitWithProportion( 0.25 )

#trndata._convertToOneOfMany( )
#tstdata._convertToOneOfMany( )

print "Number of training patterns: ", len(trndata)
print "Input and output dimensions: ", trndata.indim, trndata.outdim
print "First sample (input, target, class):"
print trndata['input'][0], trndata['target'][0]#, trndata['class'][0]

fnn = buildNetwork( trndata.indim, 5, trndata.outdim, outclass=SoftmaxLayer )
'''

#testNN = buildNetwork( 16, 5, 4, 2, outclass=SoftmaxLayer, bias=True )
#print testNN



