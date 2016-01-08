#!/usr/bin/env python

import numpy as np
import sys,math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-C","--c2file")
parser.add_argument("-c","--cfile")
parser.add_argument("-u","--ufile")
args = parser.parse_args()


try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

x = udata[:,0]
u = udata[:,1]
space = len(x)
space0=(x*x).argmin()
dx = x[1] - x[0]

try:
  c2data = np.genfromtxt(args.c2file)
  c2 = c2data[:,2].reshape((space,space))
  c = np.dot(c2,u)
except:
  try:
    cdata = np.genfromtxt(args.cfile)
    c = cdata[:,1]
  except:
    print >> sys.stderr,"could not open file"
    exit(1)

print np.dot(u*u,c)/np.dot(u,c)