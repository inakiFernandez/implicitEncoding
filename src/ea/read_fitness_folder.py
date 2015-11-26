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
type_file = sys.argv[2]
out_folder_name = folder.split("/")[1] + "-" + type_file
if not os.path.exists(out_folder_name):
    os.mkdir(out_folder_name)
for subdir, dirs, files in os.walk(folder):
    for fname in files:
        if fname.endswith(".log") and fname.startswith(type_file):

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
                    average_list.append(np.median(gener))
                experiments.append(floats)
            for index, x in enumerate(experiments):
                plot.draw_data([(subdir.split("/")[2] + "-" +
                                 str(index) + "-" + type_file, zip(*x))],
                               out_folder=out_folder_name)

            result = []

            for ind_exp in xrange(len(experiments)):
                for x in xrange(len(experiments[ind_exp])):
                    result.append(experiments[ind_exp][x])

            plot.draw_data([(subdir.split("/")[2] + "-allIndiv-" + type_file,
                             zip(*result))],
                           out_folder=out_folder_name)
            result_avg.append(average_list)

plot.draw_data([("median-" + type_file, result_avg)],
               out_folder=out_folder_name)
