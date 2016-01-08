#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-T","--sizethreshold",type=float,default=-1)
parser.add_argument("-p","--printthreshold",type=float,default=1e-20)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

t = data[:,0]
x = data[:,1]
c = data[:,2]

dx = x[1]-x[0]

xindex = np.array(x/dx,dtype=int)

x0 = np.unique(xindex)
i = 1

if args.sizethreshold > 0:
  alltrajectories = False
else:
  alltrajectories = True

for xx in x0:
  curt = t[xindex==xx]
  curx = x[xindex==xx]
  curc = c[xindex==xx]
  if not alltrajectories:
    if len(curc[curc >= args.sizethreshold]) == 0:
      continue
  for j in range(len(curt)):
    if curc[j] > args.printthreshold:
      print i,curt[j],curx[j],curc[j]
  print
  i += 1