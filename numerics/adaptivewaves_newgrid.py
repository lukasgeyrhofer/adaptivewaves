#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",required=True)
parser.add_argument("-s","--space",type=int,default=2000)
parser.add_argument("-z","--space0",type=int,default=1000)
parser.add_argument("-d","--dx",type=float,default=5e-2)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open files"
  exit(1)

x = data[:,0]
c = data[:,1]

nx = (np.arange(args.space)-args.space0)*args.dx
nc = np.interp(nx,x,c)

for i in range(args.space):
  print nx[i],nc[i]