# -*- coding: utf-8 -*-

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------

import re
from os import listdir
from os.path import isfile, join
from scipy.stats import *
from pylab import *
import brewer2mpl

# find datalog in a directory coresponfing to one experiment
# in  : a directory path
# out : a tuple of name and a list of file name (name, [datalog, ...])
def list_datalogs(path):
    logs = [ f for f in listdir(path)
             if isfile(join(path,f)) and f.startswith('datalog') ]

    datalogs = map(lambda f: join(path,f), logs)
    return datalogs

def list_logfiles(path):
    files = [ f for f in listdir(path)
              if isfile(join(path,f)) and (f.endswith('Nav') or f.endswith('For') or f.endswith('ForWithOdN') or f.endswith('NavWithOdN'))]
    logs = map(lambda f: join(path,f), files)
    return logs

def list_species(path):
    files = [ f for f in listdir(path)
              if isfile(join(path,f)) and (f.endswith('.species'))]
    logs = map(lambda f: join(path,f), files)
    return logs

def list_sizeNN(path):
    files = [ f for f in listdir(path)
              if isfile(join(path,f)) and (f.endswith('.sizeNN'))]
    logs = map(lambda f: join(path,f), files)
    return logs
def list_srfiles(path):
    files = [ f for f in listdir(path)
              if isfile(join(path,f)) and f.endswith('.log.sr2') ]

    logs = map(lambda f: join(path,f), files)
    return logs

# read a datalog and extract the [fit:##.###] component
# in  : a file name
# out : returns a list of fitness values.
def process_datalog(fname):
    d=[]
    fh = open(fname, 'r')
    for line in fh :
        data = line.split()
        for o in data :
            if o.startswith('fit', 1, 4):
                f = re.sub('[:\[\]]', ' ', o)
                f = float(f.split()[1])
                d.append(f)
    fh.close()
    return d
# read a logfile and extract the [fitness:##.### component
# in  : a file name
# out : returns a list of fitness values.
def process_logfile(fname):
    d=[]
    fh = open(fname, 'r')
    for line in fh :
        data = line.split()
        for o in data :
            if o.startswith('[fitness', 0, 8):
                f = re.sub('[:]', ' ', o)
                f = float(f.split()[1])
                d.append(f)
    fh.close()
    return d

# read a logfile and extract the [species:##.### component
# in  : a file name
# out : returns a list of number of species.
def process_species(fname):
    d=[]
    fh = open(fname, 'r')
    for line in fh :
        data = line.split()
        for o in data :
            if o.startswith('[species', 0, 8):
                f = re.sub('[:]', ' ', o)
                f = float(f.split()[1])
                d.append(f)
    fh.close()
    return d

# read a logfile and extract the [sizeNN:##.### component
# in  : a file name
# out : returns a list of sizes of Neural Network.
def process_sizeNN(fname):
    d=[]
    fh = open(fname, 'r')
    for line in fh :
        data = line.split()
        for o in data :
            if o.startswith('[sizeNN', 0, 8):
                f = re.sub('[:]', ' ', o)
                f = float(f.split()[1])
                d.append(f)
    fh.close()
    return d

# read the survival rate file
def process_srfile(fname):
    d=[]
    fh = open(fname, 'r')
    for line in fh :
        d.append(float(line))
    fh.close()
    return d

# statistics measures
def ave_accu_sf(data, cut=0.1):
    gen = int(cut * len(data[0]))
    result=[]
    for d in data:
        result.append(average(d[-gen:]))
    return result

def fix_budg_sf(data, cut=0.9):
    gen = int(cut * len(data[0]))
    result=[]
    for d in data:
        result.append(d[gen])
    return result

def time_reach_target(data, maxValue,pers=0.9):
    target = pers*maxValue
    result=[]
    for d in data:
        for g in xrange(len(d)):
            if d[g] >= target :
                break
        result.append(g)
    return result

def acc_above_target(data, maxValue):
    result=[]

    for d in data :
        aux = 0.0
        for fit in d :
            if fit > maxValue :
                aux = aux + (fit - maxValue)
        result.append(aux)
    return result

def stat_test(x,y, ks=False):
    if ks:
        return ks_2samp(x,y)
    else :
        return mannwhitneyu(x,y)

def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "-"

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

#returns the maximum swarm fitness value
#in the experiment
def extract_max(path):
    logfiles = list_logfiles(path)

    L = []
    for l in logfiles:
        L.append(process_logfile(l))
    return max(max(L))

# process multiple datalogs
def process_experiment(path, maxFit):
    sizeLogs = list_sizeNN(path)
    speciesLogs = list_species(path)

	##Inaki
    logfiles = list_logfiles(path)
    L = []
    NN = []
    Sp = []

    for l in logfiles:
        L.append(process_logfile(l))

    for nn in sizeLogs:
        NN.append(process_sizeNN(nn))

    for sp in speciesLogs:
        Sp.append(process_species(sp))

    S = {}
    S['aasf'] = ave_accu_sf(L)
    S['fbsf'] = fix_budg_sf(L)
    S['trt']  = time_reach_target(L,maxFit)
    S['aat']  = acc_above_target(L, maxFit)

    statsNN = {}
    statsNN['avgAcc'] = ave_accu_sf(NN)

    return L,S,Sp,NN,statsNN #,R

# draw a figure

def plot_one_curve(data, color, axis, label, quartiles=False):
    print data
    med, perc_25, perc_75 = perc(data)
    if quartiles :
        axis.fill_between(np.arange(0, len(med)), perc_25, perc_75,
                          alpha=0.25, linewidth=0, color=color)
    axis.plot(med, lw=5, label=label, color=color)


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

def plot_boxplot(stats, colors, axis, labels, sig=False):

    bp = axis.boxplot(stats)

    for i in range(0, len(bp['boxes'])):
        bp['boxes'][i].set_color(colors[i])
        # we have two whiskers!
        bp['whiskers'][i*2].set_color(colors[i])
        bp['whiskers'][i*2 + 1].set_color(colors[i])
        bp['whiskers'][i*2].set_linewidth(2)
        bp['whiskers'][i*2 + 1].set_linewidth(2)
        # top and bottom fliers
        # (set allows us to set many parameters at once)

        #bp['fliers'][i * 2].set(markerfacecolor=colors[i],
        #                        marker='o', alpha=0.75, markersize=6,
        #                        markeredgecolor='none')

        #bp['fliers'][i * 2 + 1].set(markerfacecolor=colors[i],
        #                            marker='o', alpha=0.75, markersize=6,
        #                            markeredgecolor='none')

        bp['medians'][i].set_color('black')
        bp['medians'][i].set_linewidth(3)
        # and 4 caps to remove
        for c in bp['caps']:
            c.set_linewidth(0)

    for i in range(0, len(bp['fliers'])):
        bp['fliers'][i].set(markerfacecolor=colors[i], marker='o', alpha=0.75,
                                markersize=6, markeredgecolor='none')

    for i in range(len(bp['boxes'])):
        box = bp['boxes'][i]
        box.set_linewidth(0)
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX,boxY)
            boxPolygon = Polygon(boxCoords, facecolor = colors[i], linewidth=0)
            axis.add_patch(boxPolygon)



    axis.set_xticklabels(labels)
    for tick in axis.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
                # specify integer or one of preset strings, e.g.
                tick.label.set_fontsize(30)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.set_axisbelow(True)


    # stat test
    if sig :
        for i in xrange(len(stats)) :
            for j in xrange(i+1,len(stats)) :
                y_max = max(concatenate((stats[i], stats[j])))
                y_min = min(concatenate((stats[i], stats[j])))
                z,p = stat_test(stats[i], stats[j])

                axis.annotate("", xy=(i+1, y_max), xycoords='data',
                             xytext=(j+1, y_max), textcoords='data',
                             arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                             connectionstyle="bar,fraction=0.1"))
                axis.text((j-i)/2.0 + j, y_max + abs(y_max - y_min)*0.1, stars(p*2.0),
                         horizontalalignment='center',
                         verticalalignment='center')






def draw_data(exp, text = "", runs=False, tex=False):
    font = {'family' : 'sans', 'size'   : 20}
    if tex :
        matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', **font)
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors

    figure(num=None, figsize=(5, 2.5), dpi=100)
    clf()

    ax1 = subplot2grid((1,1), (0,0))
    c=0

    for e in exp:
        (n, data) = e
        plot_one_curve(data, colors[(c-5) % len(colors)], ax1,  re.sub('[_/]', '', n), True)
        c=c+1
    #ax1.set_title('Avg. Swarm Fitness over time (%d runs)'%(len(data)))
    #ax1.legend(loc='lower right',prop={'size':25})
    #ax1.set_xlabel('Time (x100)')
    #ax1.set_ylabel('Fitness')


    #figtext(.08, .02, text, fontsize=10)

    #plt.tight_layout()

    draw()
    show()
