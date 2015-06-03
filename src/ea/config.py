#!/usr/bin/python
import matplotlib.pyplot as plt
import random
import numpy as np
import copy
import sys
from scipy import misc
import ConfigParser, os
import getopt

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

'''

'''
def mapConfig(section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"l:",["config="])
   except getopt.GetoptError:
      print 'test.py -l <configfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
          print 'test.py -l <configfile>'
          sys.exit()
      elif opt in ("-l", "--config"):
          c = arg
   #print 'Config file is "', c
   return [c]

if __name__ == "__main__":

    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)
    [c] = main(sys.argv[1:])
    config = ConfigParser.ConfigParser()
    config.read(c)
    
    #print c
    #print config.sections()
    #print config
    
    a = AutoVivification()
    C = []
    for s in config.sections():
         a[s] = mapConfig(s)

    print a


