# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:47:03 2015

@author: fernandi
"""
import re
import sys
sys.path.insert(0, '../plots')
import plot

filename = sys.argv[1]

f = open(filename, "r")

list_lines = f.readlines()

experiments = []
for exp in list_lines:
    list_str = []
    for elt in exp[0:-1].split("], ["):
        list_str.append(elt.replace('[', "").replace(']', ""))

    floats = []

    for gen in list_str:
        gener = []
        for ind in re.findall(r"'([^']*)'", gen):
            gener.append(float(ind))
        floats.append(gener)
    experiments.append(floats)

#for x in experiments:
#    plot.draw_data([("test", x)])

result = []

for ind_exp in xrange(len(experiments)):
    for x in xrange(len(experiments[ind_exp])):
        result.append(experiments[ind_exp][x])

plot.draw_data([("test", zip(*result))])
