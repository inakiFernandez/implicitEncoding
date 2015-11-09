# -*- coding: utf-8 -*-
"""
Read a set of independent runs of the same experiment in the same folder
and plots the corresponding fitness values
@author: fernandi
"""
import re
import sys
sys.path.insert(0, '../plots')
import plot
import os
import numpy as np

folder = sys.argv[1]
result_avg = []
for subdir, dirs, files in os.walk(folder):
    for fname in files:
        if fname.endswith(".log") and fname.startswith("fitness"):

            f = open(os.path.join(subdir, fname), "r")

            list_lines = f.readlines()

            experiments = []
            average_list = []
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
                    #print np.average(gener)
                    average_list.append(np.average(gener))
                experiments.append(floats)
            #for x in experiments:
            #    plot.draw_data([("test", x)])

            result = []

            for ind_exp in xrange(len(experiments)):
                for x in xrange(len(experiments[ind_exp])):
                    result.append(experiments[ind_exp][x])

            #plot.draw_data([("test", zip(*result))])
            #print result_avg
            result_avg.append(average_list)
#result_avg = zip(*result_avg)
#print result_avg
print len(result_avg)
print len(result_avg[0])
plot.draw_data([("test", result_avg)])
