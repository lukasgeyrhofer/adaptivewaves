#!/usr/bin/env python

import numpy as np
import argparse
import math,sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-v","--speed",type=float,default=1.)
parser.add_argument("-s","--sigma",type=float,default=1e-3)
parser.add_argument("-d","--dx",type=float,default=1.)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open files"
  exit(1)

re = data[:,1]
im = data[:,2]

for i in range(len(re)):
  print re[i]/args.sigma,im[i]*args.dx/args.speed