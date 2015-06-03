import matplotlib.pyplot as plt
import random
import numpy as np
import copy
import sys

from scipy import misc

###Aux object-oriented dict class###
class objdict(dict):
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
	    self[name] = value
	    
    def __delattr__(self, name):
	    if name in self:
		    del self[name]
	    else:
		    raise AttributeError("No such attribute: " + name)

########GENOME#########

class Genome:
	def __init__(self,nbIn,nbOut, nbLayers,neuronsPerLayer,genomeSize):
		#initialization
		self.nbI = nbIn
		self.nbO = nbOut

		d =  dict()
		d["nblay"] = nbLayers
		d["neurPerLay"] = neuronsPerLayer
		self.nnSpec = []
		od = objdict(d)
		self.nnSpec = od
		self.nnSpec["nblay"] = nbLayers
		self.nnSpec["neurPerLay"] = neuronsPerLayer
		
		#Counting maximum number of connections 
		if nbLayers == 0:
			self.maxConnections = nbIn * nbOut
		else:
			self.maxConnections = nbIn * neuronsPerLayer
			for i in range(nbLayers-1):
				self.maxConnections += neuronsPerLayer * neuronsPerLayer
			self.maxConnections += neuronsPerLayer * nbOut
			
		self.genome = [True] * genomeSize
		self.nn = None
		#Fully random initialization
		for i in range(genomeSize):
			if random.random() < 0.5:
				self.genome[i] = False
				
		#Possibility of explicitly introducing valid genes, 
		#codons or any other initial blueprint
		#HERE

	#####getters
	
	def getNbLay(self):
		return self.nnSpec["nblay"]
	
	def getNeurPerLay(self):
		return self.nnSpec["neurPerLay"]

	def getGenome(self):
		return self.genome

	def getNN(self):
		return self.nn
	
	def getMaxConnections(self):
		return self.maxConnections



#Test miniscripts
'''
import genome
from genome import Genome


g1 = Genome(9,2,0,0,100)

g2 = Genome(9,2,1,5,100)

g3 = Genome(9,2,2,5,100)


print g1.getMaxConnections()

print g2.getMaxConnections()

print g3.getMaxConnections()

'''
