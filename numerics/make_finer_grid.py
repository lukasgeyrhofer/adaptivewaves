#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",required=True)
parser.add_argument("-r","--ratio",type=int,default=2)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)


x = data[:,0]
u = data[:,1]

au = np.append(u,2*u[-1]-u[-2])

dx = x[1]-x[0]
space=len(x)
space0=(x*x).argmin()

newdx = dx/args.ratio
newspace = space*args.ratio
newspace0 = space0*args.ratio

newx = (np.arange(newspace)-newspace0)*newdx
newu = np.zeros(newspace)

for i in range(space):
  for j in range(args.ratio):
    newu[args.ratio*i+j] = au[i] + float(j)/float(args.ratio)*(au[i+1]-au[i])

for i in range(newspace):
  print newx[i],newu[i]