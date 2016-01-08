#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-n","--extension",type=int,default=2)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

x = data[:,0]
c = data[:,1]

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1]-x[0]

ndx = dx/args.extension
nspace = space*args.extension
nspace0 = space0*args.extension

nx = (np.arange(nspace)-nspace0)*ndx

nc = np.interp(nx,x,c)

for i in range(nspace):
  print nx[i],nc[i]

