# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
#pylint: disable=pointless-string-statement
from __future__ import division
"""
    Launcher, entry point for the pleiotropy experiments
"""
#TODO (lambda + mu) - ES
#TODO check for memory leaks (due to NN/genomes copy)
#TODO check imports are useful. Static code analysis.
#TODO To uncomment following lines
#from ea.genome import Genome
#TODO check GITification. 
#TODO Ssh Keys not working correctly
from pybrain.structure import FeedForwardNetwork, LinearLayer
from pybrain.structure import FullConnection, SigmoidLayer, BiasUnit
#from pybrain.structure import RecurrentNetwork
import numpy as np
import matplotlib.pyplot as plt
#from pybrain.structure.evolvables.evolvable import Evolvable
#TODO lots of sanity checks, such as tournSize < popSize
#TODO Final integration test. Refactor code into different classes

###############################################################################
####################################Plotting###################################
import sys # os, re
#from os import listdir
#from os.path import isfile, join
#from scipy.stats import *
#from pylab import *

########################Plotting params and functions##########################
import brewer2mpl
bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
colors = bmap.mpl_colors   
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
#########################NN structure definitions##############################
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
    
def decode(g, i, o, nbH, nPerH, b):
    n = fullConnectionNet(i, o, nbH, nPerH, b)    
    for i in xrange(n.paramdim):
        n.params[i] = g[i]
    n.reset()
    return n
    
###############################################################################
    
def evaluate(indiv, pbDb, labelDb):
    fit = 0
    for i in xrange(len(pbDb)):
        answer = indiv.activate(pbDb[i])
        #Quadratic error
        fit += (answer - labelDb[i])**2    
    #?Balance fitness w.r.t. proportion of same-valued answers:
    #e.g. 99 True, 1 False
    return -fit

def select(p, f, nbSel):
    sortedPop= [x 
        for (y,x) in sorted(zip(f,p), key=lambda pair: pair[0], reverse=True)]
    #Truncation selection ","
    parents = sortedPop[:mu]
           
    return np.mean(parents, axis=0)

def sampleOffspring(centr, l, sig, n):
    return np.random.multivariate_normal(centr, np.identity(n) * sig, l)
'''
Deprecated    
def mutate(indiv, sigma):
    result = indiv.copy()

    for i in xrange(result.paramdim):
        result.params[i] = np.random.normal(result.params[i], sigma)
    
    return result
'''
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

##############NN fixed structure definition - fully connected MLP##############
percepts = 5; out = 1; nbHLayers = 1; neurPerLayer = 3; bias = True
#Unnecessary if initialized centroid:
#Boundaries of initial weights. minWVal = -10; maxWVal = +10
#########################Training DB construction##############################
pbDb = []; labelDb = []; pbSz = percepts; dbSz = 100
problemID = "MIN"#"MIN", "AND", "OR", "MAJ"
#TODO TOFIX: "MUX1" (distance definitions in N-dimensional spaces needed)
for i in xrange(dbSz):
    v = np.random.rand(pbSz)
    for i in xrange(pbSz):
        # Or random.choice( [0,1] )
        v[i] = round(v[i], 0) 
    pbDb.append(v)    
    #TODO Refactor into independent "task" module    
    #Learn random associations (impossible task)
    #labelDb.append(round(np.random.rand(),0)) 
    labelDb.append(task(v, problemID))
    
##################Evolutionary algorithm (mu/mu, lambda)-ES####################
sigma0 = 0.5; sigma = sigma0; lmbda = 20; mu = int(np.floor(lmbda / 2  ));
maxGen = 10; bestEver = None #Best individual
sigmas = [0]*(maxGen + 1); sigmas[0] = sigma
########################Population initialization##############################
#TODO Other connectivity patterns
#Implemented: fully connected MLP with 0-* hidden layers with the same number
#of neurons per layer optional bias 
auxNN = fullConnectionNet(percepts, out, nbHLayers, neurPerLayer, bias)
spaceDim = auxNN.paramdim
#Initial centroid of the distribution
m0 = [0] * spaceDim; m = m0
centroids = [[0]*spaceDim for i in range(maxGen + 1)]; centroids[0] = m
#parents = (maxWVal - minWVal) * np.random_sample((mu, spaceDim)) - minWVal
fitness = [[0]*lmbda for i in range(maxGen + 1)]
print "Starting evolution of neural nets for " + problemID + " problem with "
print str(pbSz) + " boolean inputs"
print "Search space : R^" + str(spaceDim)

#TODO TOCHECK termination conditions and data structures boundaries
#TODO other stopping criteria    
#############################Start Evo loop####################################
for g in xrange(maxGen):       
    print "Generation: " + str(g)
    offspring = sampleOffspring(m, lmbda, sigma, spaceDim)
    for ind in xrange(lmbda):
        fitness[g][ind] = evaluate(
                            decode(offspring[ind],
                                   percepts, out, nbHLayers,
                                       neurPerLayer, bias),
                            pbDb, labelDb)
    #TODO sort to keep bestEver        
    #TODO check selection and mutation/sampleOffspring methods
    m = select(offspring, fitness[g], mu)
    centroids[g + 1] = m
    #TODO Update sigma and store it    
    
print "Generation: " + str(maxGen)         
#Last generation eval
offspring = sampleOffspring(m, lmbda, sigma, spaceDim)
for ind in xrange(lmbda):    
    fitness[maxGen][ind] = evaluate(decode(offspring[ind], percepts, out, 
                                nbHLayers, neurPerLayer, bias), pbDb, labelDb)

axis = plt.subplot2grid((1,1), (0,0))
plot_one_curve(zip(*fitness), colors[0], axis, "Fitness", True)