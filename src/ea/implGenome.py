# -*- coding: utf-8 -*-
#TODO General list: config global params(probas, start and end codons, other)
#test mutations and boundary conditions on fragments etc.
import time
import random, sys, math
import numpy as np
from pybrain.structure import FeedForwardNetwork, LinearLayer
from pybrain.structure import FullConnection, SigmoidLayer
import ConfigParser, itertools
sys.path.insert(0, '../plots')
import plot
import colorama as clr
import brewer2mpl
from pylab import subplot2grid

#Data structures and functions associated with the creation
#mapping to neural net (pybrain), mutation, recombination
#TODO Learn random associations (~impossible (overfitted) task)
#labelDb.append(round(np.random.rand(),0))         

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def mapConfig(section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def doNothing(g):
    pass
    return g

def randomGenome(size, symbols = ['0','1']):
    g = [];    
    for i in xrange(size):
        g.append(random.choice(symbols));
    return g;    

def task(vector,taskID):
    tasks= {"AND": all,
              "OR": any,
              "MAJ": majority,
              "MIN": minority,
             # "MUX1": mux1,
              "RETAND": retinaAnd,
              "RETOR": retinaOr
              }    
    intVector= [int(round(x)) for x in vector]              
    return tasks[taskID](intVector)

#retina "and" problem
def retinaAnd(x):
    #input vector is arranged in [left1,l2,r1,r2, (2nd row) l3,l4,r3,r4]
    lInput = [v for idx,v in enumerate(x) if idx in [0,1,4,5]]
    rInput = [v for idx,v in enumerate(x) if idx in [2,3,6,7]]
    return isLObject(lInput) and isRObject(rInput)

#retina "or" problem
def retinaOr(x):
    #input vector is arranged in [left1,l2,r1,r2, (2nd row) l3,l4,r3,r4]
    lInput = [v for idx,v in enumerate(x) if idx in [0,1,4,5]]
    rInput = [v for idx,v in enumerate(x) if idx in [2,3,6,7]]
    return isLObject(lInput) or isRObject(rInput)

def isLObject(lX):   
    if lX.count(1) >= 3:
        return True
    else:
        #Only one or two pixeles in left column, none in right column
        if (lX[1] == 0) and (lX[3] == 0):
            if (lX[0] == 1) or (lX[2] == 1):
                return True;
    return False
    
def isRObject(rX):
    if rX.count(1) >= 3:
        return True
    else:
        #Only one or two pixeles in right column, none in leftcolumn
        if (rX[0] == 0) and (rX[2] == 0):
            if (rX[1] == 1) or (rX[3] == 1):
                return True;
    return False  
#Binary majority
def majority(x):
    total = sum(x)
    if total >= len(x)/2.0:
        return True
    else:
        return False

#Binary minority
def minority(x):
    return not majority(x)
            
###############################################################################
'''
def mux1(x): 
    #Multiplexer type 1: 
     #   if x[0] is False, then return x[1:ceil((len(x)-1)/2)]
      #  else if x[0] is True, then returnx[ceil((len(x)-1)/2) + 1: len(x) - 1]
    
    if x[0] == False:
        return x[1:np.ceil((len(x) - 1)/2)]
    else:
        return x[np.ceil((len(x)-1)/2) + 1: len(x) - 1]            

def mate(g1,g2):
    #one or two point crossover? something else?
    offspring = None    
    
    return offspring;

def extractExpressedCodons(genome, totalNbLinks, sC, eC, symbols = ['0','1']):

    allCodons, outStr = extractCodons(genome, totalNbLinks, sC, eC, symbols = ['0','1']);
    expressedLinks = [];
    eCodons = []    
    
    if config.DominanceScheme == "first":
       for c in allCodons:
           if not(c["target"] in expressedLinks):
               eCodons.append(c);
               expressedLinks.append(c["target"]);
    else:
        print "Other dominance schemes to be implemented."
    
    return eCodons;

'''

def fragCopy(g, symbols = ['0','1']):
    offspring = [];
    size = len(g);    
    #Draw a first position
    first = random.choice(xrange(size + 1))    
    #Draw a second position (after "first" position)
    secondRange = range(first, size + 1)#[first:size + 1]    
    second = random.choice(secondRange)    
    at = random.choice(xrange(size + 1))       
    
    #TODO test boundary conditions
    for i in xrange(size):
        if i == at:
            offspring.extend(g[first:second])
            
        offspring.append(g[i])
        
    if at == size + 1:        
        offspring.extend(g[first:second])
        
    return offspring;    
    
def fragMove(g, symbols = ['0','1']):
    offspring = [];
    size = len(g);
    #TODO IMPORTANT TOTEST    
    #Draw a first position
    first = random.choice(xrange(size + 1))    
    #Draw a second position (after "first" position)
    secondRange = range(first, size + 1)    
    second = random.choice(secondRange)
    frag = g[first:second]
    #Insert drawn fragment at random position in the pruned offspring    
    at = random.choice(xrange(len(offspring) + 1))              

    offspring[:] = [item for i,item in enumerate(g) if not(i in xrange(first, second))]
    
    beginning = offspring[0:at]
    end = offspring[at:]
    offspring = beginning + frag + end
            
    return offspring;    
        
def fragDel(g, symbols = ['0','1']):
    offspring = [];

    #TODO IMPORTANT TOTEST
    size = len(g);
    
    #Draw interval to delete
    first = random.choice(xrange(size + 1))    
    secondRange = range(size + 1)[first:size + 1]    
    second = random.choice(secondRange)

    offspring[:] = [item for i,item in enumerate(g) if not(item in xrange(first, second))]

    return offspring;    

#For binary genomes only. Other mutation operators to define in non-binary case
def bitFlip(g, p):
    offspring = g
    
    for i,c in enumerate(offspring):
        if np.random.uniform() < p:
            if c == '0':
                offspring[i] = '1'
            else:
                if c == '1':
                    offspring[i] = '0'
                else:
                    sys.exit("Wrong char in genome")
    
    return offspring

def randomGenomeWithGenes(junk, nbLinks, sC, eC,\
                             nuclPerGene, symbols = ['0','1']):    
    g = []
    
    interstices = nbLinks + 1;
    junkPerInterstice = int(math.floor(junk /interstices))
    tgtSize = max(1, math.ceil(math.log(nbLinks,len(symbols))));
    
    g = g + list(np.random.choice(symbols, size = junkPerInterstice))
    #loop for each link gene
    for i in xrange(nbLinks):
        #Start codon
        g = g + list(sC);
        #Targeted link ID ' + str(int(tgtSize)) + "
        target = list(("{0:0%sb}" % str(int(tgtSize))).format(i)) 

        g = g + target;
        #Some random nucleotids for weight field        
        for j in xrange(nuclPerGene):
            g = g + list(random.choice(symbols))
        #End codon
        g = g + list(eC);    
        g = g + list(np.random.choice(symbols, size = junkPerInterstice))

    return g;    
    
def createGenome(nLinks, fractionJunkGenes, nucleotPerGene, \
             startCodon =  None, endCodon = None,  symbols = ['0','1']):
                 
    targetSize = int(math.ceil(max(1, math.log(nLinks,len(symbols)))))
    sze =  nLinks * (targetSize + len(startCodon) + len(endCodon)  + nucleotPerGene)
    effSize = int(round((1.0/(1.0 - fractionJunkGenes)) * sze))

    g = randomGenomeWithGenes(int(math.floor(effSize * fractionJunkGenes)), 
                                  nLinks, startCodon, endCodon, nucleotPerGene)
    #TODO TOCHECK Robustness through redundancy?
    return g;


def mapToStdFFNN(codons, targetSize, nIn, nOut, nHid = 0, neurPerHid = 0,
                      bias = False, polygeneStrategy = "avg"):
    net = None    
    
    if nHid <= 0:
        nLinks = nIn * nOut;
    else:
        nLinks = nIn * neurPerHid + neurPerHid * nOut + \
                    (nHid - 1) * (neurPerHid * neurPerHid)              
    #Neurons' IDs starts counting from 0 in the sequence:
    #Input - Bias - Output - (Hidden)
    #Neuron IDs are fixed and topology does not evolve for now
    #Link order : cf. below

    #TODO evolution of topology?
    
    net = FeedForwardNetwork()
    
    for i in xrange(nIn):
        net.addInputModule(LinearLayer(1,name="i" + str(i)))        
    
    for i in xrange(nOut):
        net.addOutputModule(SigmoidLayer(1,name="o" + str(i)))

    for i in xrange(nHid):
        for j in xrange(neurPerHid):
            net.addModule(SigmoidLayer(1, name="h" + str(i) + "-" + str(j)))
        
    #Loop link codons
    #Link ordering: 
    #    No hidden layer: 
    #              0        1     ....
    #           [i0,o0], [i1,o0], .... , [i0,o1], [i1,o1], ...., [iN,oN'] 
    #   W. hidden layer:
    #           [i0,h00],[i1,h00],...,[i0,h01],...
    #           [h00,h10],[h01,h10],..., [hK0,o0], [hK1,o0],..., [hKN'',oN']
    i = 0; o = 0; h = 0; 
    
    #TODO attention on multiple codons 
    targetSet = set(xrange(nLinks))
    groupedCodons = [[x, [y[1] for y in codons if y[0]==x]]  for x in targetSet]
    for l in xrange(nLinks):
        geneWeights = [y[1] for y in groupedCodons if y[0] == l]
        geneWeights = list(itertools.chain.from_iterable(geneWeights))
        
        if len(geneWeights) == 0:
            w = 0.0
        else:
            #Different polygenic schemes, dominance etc 
            #Average
            w = sum(geneWeights)/float(len(geneWeights))
            #TODO other schemes
        if nHid == 0:
            inNeur = net["i" + str(i)]
            outNeur =  net["o" + str(o)]            
            if i == nIn - 1:
                #To first in neuron for next out neuron                
                i = 0;
                if o == nOut - 1:
                    #Finished setting up all links
                    o = 0
                    conn = FullConnection(inNeur, outNeur, name = str(l))
                    net.addConnection(conn)
                    conn.params[0] = w
                    break
                else:
                    #Next out neuron
                    o = o + 1;
            else:
                #Next in neuron
                i = i + 1;
        else:
            if h == 0:
                inNeur = net["i" + str(i)]
                outNeur =  net["h" + str(h) + "-" + str(o)]
                if i == nIn - 1:
                    i = 0
                    if o == neurPerHid - 1:
                        o = 0
                        h = h + 1
                    else:
                        o = o + 1
                else:
                    i = i + 1
            else:
                if h < nHid:
                    inNeur = net["h" + str(h - 1) + "-" + str(i)]
                    outNeur = net["h" + str(h) + "-" + str(o)] 
                    if i == neurPerHid - 1:
                        i = 0
                        if o == neurPerHid - 1:
                            o = 0
                            h = h + 1
                        else:
                            o = o + 1
                    else:
                        i = i + 1
                else:
                    inNeur = net["h" + str(h - 1) + "-" + str(i)]
                    outNeur = net["o" + str(o)]                     
                    if i == neurPerHid - 1:
                        i = 0
                        if o == nOut - 1:
                            o = 0
                            conn = FullConnection(inNeur, outNeur, name = str(l))
                            net.addConnection(conn)
                            conn.params[0] = w
                            break
                        else:
                            o = o + 1
                    else:
                        i = i + 1        
        conn = FullConnection(inNeur, outNeur, name = str(l))
        net.addConnection(conn)
        conn.params[0] = w
        
    '''
        if bias:
            b = BiasUnit(name = "b")
            net.addModule(b)
                Bias??
    '''
    net.sortModules()

    return net;

def extractCodons(genome, tgtSize, sC, eC, symbols = ['0','1'], maxW = 5.0):

    Codons = [];  i = 0; outStr = ""; spacing = " \n---\n "
    
    while i in xrange(len(genome)):
        if "".join(genome[i:]).startswith(sC):
            i = i + len(sC);
            outStr = outStr + spacing + clr.Back.MAGENTA + sC + clr.Back.RESET;
            
            tgt = genome[i:i + tgtSize];
            outStr = outStr + clr.Back.CYAN + "".join(tgt) + clr.Back.RESET;
            i = i + tgtSize;

            j = i;
            while not("".join(genome[j:]).startswith(eC)) and j < len(genome):
                j = j + 1;
                #TODO condition if end genome and not finished gene
            if j <= len(genome):
                weightStr = genome[i:j];
                outStr = outStr + clr.Back.GREEN + "".join(weightStr) + clr.Back.RESET;
                if len(weightStr) == 0:
                    Codons.append([int("".join(tgt), 2),0.0])            
                else:
                    Codons.append([int("".join(tgt), 2),
                               2 * maxW * weightStr.count("1")/len(weightStr) - maxW])            
            outStr = outStr + clr.Back.MAGENTA + eC + clr.Back.RESET + spacing;
            i = j + len(eC);
        else:
            outStr = outStr + str(genome[i])
            i = i + 1;
    return Codons, outStr;
###############################################################################       
def mutate(g):
    #pFrCpy = 0.12;  pFrMove = 0.12; pFrDel = 0.12; pBitFlip = 0.01 ; #pGeneIns= 0.001;
    pFrCpy = 0.05;  pFrMove = 0.10; pFrDel = 0.05; pBitFlip = 0.03 ; #pGeneIns= 0.001;
        
    muts = {0 : fragCopy,
            1 : fragMove, 
            2 : fragDel,
            3 : doNothing} #4 : [geneIns, pGeneIns]};

    pNothing = 1 - pFrCpy - pFrMove - pFrDel
    
    mut = np.random.choice(xrange(0,4),p = [pFrCpy, pFrMove, pFrDel, pNothing]);
    #First, copy, move or delete a random fragment of the genome at random pos.
    #Second, bit-flip for each bit of the genome with small probability
    offspring = bitFlip(muts[mut](g), pBitFlip)           
       
    return offspring; 

def evaluate(indiv, pbDb, labelDb, nbInstances = -1):
    fit = 0
    if nbInstances == -1:        
        instances = xrange(len(pbDb))
    else:
        #Randomly sample "nbInstances" instances without replacement
        instances = np.random.choice(xrange(len(pbDb)), nbInstances, replace=False)
    for i in instances:
        answer = indiv.activate(pbDb[i])
        #Error on each instance (euclidian distance to right answer) 
        #is normalized between 0 and 1
        error = np.linalg.norm(np.subtract(answer, labelDb[i])) / math.sqrt(len(pbDb[i]))
        fit += -error
    #Normalized for the size of the database: average euclidian error over instances
    return fit/float(len(pbDb)) 

def select(fitness, nbParents):        
    parents = [] 
    sortedIdx = [i[0] for i in sorted(enumerate(fitness), key=lambda x:x[1])]
    rankWeight = [i[0] for i in enumerate(sortedIdx)]    
    triangNumber = len(fitness) * (len(fitness) + 1)/2
    rankProba = [float(i + 1)/float(triangNumber) for i in rankWeight]

    for i in xrange(nbParents):
        c = np.random.choice(sortedIdx, p = rankProba)                
        parents.append(c)

    condition = (nbParents == len(parents))
    assert condition, "Wrong number of selected parents" 
    return parents
    
def survive(f, mu, fullPop):
    sortedIdx = [i[0] for i in sorted(enumerate(f), key=lambda x:x[1])]
    survivorsIdx = sortedIdx[-mu:]
    
    genomes = [fullPop[i] for i in survivorsIdx]

    fitP = [f[i] for i in survivorsIdx]
    
    return survivorsIdx, genomes, fitP

if __name__ == "__main__":
    clr.init()
    '''    
    g  = list(np.random.choice(['0','1'], size = 100 ))  
    print g, "\n", bitFlip(g,1.0)        
    '''
    '''
    n = 10000
    k = 1000000
    indexes = xrange(1, n + 1)
    res = select(indexes,k)
    test = [float(i) / k for i in np.histogram(res,n)[0]]
    teor =  [float(i) / (n * (n + 1) / 2) for i in xrange(1,n + 1)]
    print np.linalg.norm(np.subtract(test, teor))
    '''    
    '''
    inp = 8; o = 1; h = 0; neurH = 0;
    st = "1111111"; end = "0000000"; nuclW = 5; #st = "START"; end = "END";
    symbols = ['0','1']
    if h == 0:
        nLinks = inp * o
    else:
        nLinks = inp * neurH + (h - 1) * neurH ** 2 + neurH * o

    g = createGenome(nLinks, 0.0, nuclW, st, end)
    g  = list(np.random.choice(['0','1'], size = 100 ))  
    print "\n\n\n\n"
    print "".join(g)
    print "\n\n\n\n"
    tgtSize  = int(math.ceil(max(1, math.log(nLinks,len(symbols)))));
    cod, prettyG = extractCodons(g,tgtSize, st, end)
    print prettyG
    
    pbDb = []; labelDb = []; pbSz = inp; dbSz = -1 # dbSz = -1 for whole problem
    problemID = "RETAND"#"RETAND" # "MIN", "AND", "OR", "MAJ", "RETOR"        

    #Generate all instances 
    for i in xrange(len(symbols)**pbSz):            
        #binary in this case
        v = [float(k) for k in list(("{0:0%sb}" % str(pbSz)).format(i))] 
        pbDb.append(v)            
        labelDb.append(task(v, problemID))

    net = mapToStdFFNN(cod, tgtSize, inp, o, nHid = h, neurPerHid = neurH)
    
    print evaluate(net, pbDb, labelDb, nbInstances=dbSz)    
    for c in [connection for connections in net.connections.values() for connection in connections]:
        print("{} -> {} => {}".format(c.inmod.name, c.outmod.name, c.params))    
    print net.connections
    sys.exit("End test")
    '''
    c = "config.ini"    
    config = ConfigParser.ConfigParser()
    config.read(c)
    a = AutoVivification()
    C = []
    for s in config.sections():
         a[s] = mapConfig(s)

###############################################################################
    #TODO different symbol alphabet (different base for weight, different 
    #numbering of links)
    #One output for logical binary output
    nO = 1; nI = 8 #8 for retina pb
    nHid = 2; neurPerHiddenLayer = 5    
    #Proportion of junk genes in-between genes on initialization of the genome
    fracJunkGenes = 0.2;
    symbols = ['0','1']
    #Number of nucleotids for weight field on initialization
    #(the more nucleotids, the finer the resolution of initial weights)
    nbWeightValues = 15;
    nuclPerW = nbWeightValues - 1;    
    #How to select start and end codons? Shorter (i.e. easier to randomly draw)
    #codons are more prone to disruptive variation? (are offspring viable?)
    #Grey-coding for links? minimizing "distances" between links?
    #Look for a symmetry in codons?
    #sC = "111010010100"; eC = "010100011011" ;
    #sC = "11111"; eC = "00000"
    #sC = "11010"; eC = "00110"
    sC = "111100111000110010010010"; eC = "101011101101101100001111";    
    if nHid == 0:
        nLinks = nI * nO
    else:
        nLinks = nI * neurPerHiddenLayer + (nHid - 1)*neurPerHiddenLayer**2 + \
                            neurPerHiddenLayer * nO
    tgtSize  = int(math.ceil(max(1, math.log(nLinks,len(symbols)))));
###############################################################################
    pbDb = []; labelDb = []; pbSz = nI; dbSz = -1 # dbSz = -1 for whole problem
    problemID = "RETAND"#"RETAND" # "MIN", "AND", "OR", "MAJ", "RETOR"        
    #Generate all instances
    for i in xrange(len(symbols)**pbSz):            
        #v = list(("{0:0%sb}" % str(pbSz)).format(i))
        #binary in this case, different base for other symbol sets
        v = [float(k) for k in list(("{0:0%sb}" % str(pbSz)).format(i))] 
        pbDb.append(v)            
        labelDb.append(task(v, problemID))
###############################################################################
###############################################################################
##########################Evolutionary algorithm###############################
    mu = 20; nbGenerations = 80; 
    parentChildrenRatio = 1.0 # number of children per parent
    # selectionRatio * mu == lambda nb of children
    lbda = int(round(mu * parentChildrenRatio)); 
    #mu parents in the beginning and lbda children every generation 
    maxEval = mu + lbda * nbGenerations; 
    genomes = []; nets = []; nbEval = 0; outStr = ""
    print "Init population"
    #Initial population: mu individuals 
    #Valid ones (all random links correctly encoded) OR random binary string
    for i in xrange(mu):
        #codons = [] #while len(codons) == 0:
        #g  = list(np.random.choice(symbols, size = 150 ))        
        g = createGenome(nLinks, fracJunkGenes, nuclPerW, sC, eC) 
        codons, outStr = extractCodons(g, tgtSize, sC, eC)
        net = mapToStdFFNN(codons, tgtSize, nI, nO, nHid=nHid, \
            neurPerHid=neurPerHiddenLayer)
        genomes.append(g); nets.append(net)
    
    fitP = [0.0] * mu ; it = 0        
    #Evaluate initial population
    for i, net in enumerate(nets):
        fitP[i] = evaluate(net,pbDb,labelDb, nbInstances=dbSz)
        
    sys.stdout.write(str(it) + " - ")    
    log = []; log.append(fitP); timeLog = []; codonLog = []
    
    it = it + 1
    #Total Time
    timeSt = time.time()

    #Evolutionary loop    
    for it in xrange(1, nbGenerations + 1):
        #Time per generation
        timeGSt = time.time()
        #Select lambda (possibly repeated) parents  with rank-based selection
        parents = select(fitP, lbda)           
        children = []; nets = []        
        #Generate lambda children        
        #TODO  look for disruptive mutations (e.g. by measuring valid codons)
        #Note: it is easier to break a link gene than creating a new one
        for i in parents:        
            child = mutate(genomes[i])
            codons, outStr = extractCodons(child, tgtSize, sC, eC)
            net = mapToStdFFNN(codons, tgtSize, nI, nO, nHid=nHid, \
                        neurPerHid=neurPerHiddenLayer)
            children.append(child); nets.append(net)
            
        fitChild = [0.0] * lbda
        #Evaluate children
        for i, net in enumerate(nets):
            fitChild[i] = evaluate(net,pbDb,labelDb, nbInstances=dbSz)
            
        #Truncation "plus" survivor sel./replacement [parents + children]
        survIdx, genomes, fitP = survive(fitP + fitChild, mu, genomes + children)
        #TODO TO TEST Re-evaluate genomes     
        fitParents = []
        for g in genomes:
            codons, outStr = extractCodons(g, tgtSize, sC, eC)
            net = mapToStdFFNN(codons, tgtSize, nI, nO, \
                           nHid=nHid, neurPerHid=neurPerHiddenLayer)
            fitParents.append(evaluate(net,pbDb,labelDb, nbInstances=dbSz))
            
        #print fitParents
        #print fitP
        #print fitP == fitParents
        #For codon counting only        
        nbCodons = []          
        for i in genomes:        
            codons, outStr = extractCodons(i, tgtSize, sC, eC)
            nbCodons.append(len(codons))
        codonLog.append(nbCodons)
        #Logging operations at the end of generation                
        sys.stdout.write(str(it) + " - ")        
        it = it + 1        
        log.append(fitP)                    
        timeGEnd = time.time()
        timeLog.append([timeGEnd - timeGSt])
        
###############################################################################        
    timeEnd = time.time()
  
    print "\nIt took: ", str(timeEnd - timeSt), " seconds"
    print "Time per generation (seconds)"    
    plot.draw_data([["Time",[list(x) for x in zip(*timeLog)]]])   
    print "For mu = ", str(mu), ", lambda = ", str(lbda), \
        ", number of generations = ", nbGenerations, ", problem = ", problemID
    print "\nFitness through generations ( -error )"    
    plot.draw_data([["Test",[list(x) for x in zip(*log)]]])
    print "Number of codons (of the mu survirving individuals) per generation"
    plot.draw_data([["Codons",[list(x) for x in zip(*codonLog)]]])
###############################################################################    
    #'''
    #Test (& generalization if tested on more instances) 
    popNets = []
    for g in genomes:                    
        codons, outStr = extractCodons(g, tgtSize, sC, eC)
        net = mapToStdFFNN(codons, tgtSize, nI, nO, \
                           nHid=nHid, neurPerHid=neurPerHiddenLayer)
        popNets.append(net)         
    testResults = []    
    for net in popNets:
        accuracy = 0.0
        for i in xrange(2**pbSz):
            answer = net.activate(pbDb[i])        
            if int(round(answer)) == labelDb[i]:
                accuracy = accuracy + 1.0
        accuracy = accuracy / (2 ** pbSz)
        testResults.append(accuracy)
    print "Test results (accuracy of last generation)"
    axis = subplot2grid((1,1), (0,0))
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    plot.plot_boxplot([testResults], bmap.mpl_colors, axis, ["Test"])        
    #'''
###############################################################################    
''' 
#for c in [connection for connections in net.connections.values() for connection in connections]:
#    print("{} -> {} => {}".format(c.inmod.name, c.outmod.name, c.params))    
#TODO bases other than binary randomGenome(10, list(string.ascii_lowercase))
'''
###############################################################################