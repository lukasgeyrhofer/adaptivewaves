#!/usr/bin/env python


import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile")
parser.add_argument("-V","--variance",default=-1.,type=float)
args = parser.parse_args()


try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)


x = udata[:,0]
u = udata[:,1]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

if args.variance < 0:
  v = np.sum(x[-20:-10]**2 - 2*x[-20:-10]*u[-20:-10])/10.
else:
  v = args.variance
  
  
c = np.exp(-0.5*x*x/v)
norm = np.sum(u*c)*dx
c/=norm

for i in range(space):
  print >> sys.stderr,x[i],c[i],u[i]
  
  
popsize = np.sum(c)*dx
print popsize,v
  



