#!/usr/bin/env python

import numpy as np
import scipy.optimize as opt
import sys,math
import argparse


#mutationkernel = np.zeros(1)
#xm4um = np.zeros(1)

def f(cc2):
  tmp = np.zeros((space,space))
  c1_2dx = np.dot(cc2,u)*2.
  rescale = np.dot(c1_2dx,u)*.5*dx*dx
  cc2 /= rescale
  for i in range(1,space-1):
    for j in range(1,space-1):
      tmp[i][j]  = speedprefactor*(cc2[i+1,j]-cc2[i-1,j]+cc2[i,j+1]-cc2[i,j-1])
      tmp[i][j] += xm4um[i] + xm4um[j]
      for k in range(1,min(i,j)):
	tmp[i][j] += mutationkernel[k]*(cc2[i-k,j]+cc2[i,j-k])
      if i<j:
	for k in range(i,j):
	  tmp[i,j] += mutationkernel[k]*cc2[i,j-k]
      elif j<i:
	for k in range(j,i):
	  tmp[i,j] += mutationkernel[k]*cc2[i-k,j]
      else:
	tmp[i,j] += c1_2dx[i]
  return tmp

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-c","--cfile")
parser.add_argument("-C","--c2file")
parser.add_argument("-s","--mutationsigma",type=float,default=1e-3)
parser.add_argument("-m","--mutationrate",type=float,default=1e-5)
parser.add_argument("-v","--wavespeed",type=float,default=1e-5)
parser.add_argument("-S","--maxSteps",type=int,default=None)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)

x = udata[:,0]
u = udata[:,1]/1.5

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

mutationkernel = np.exp(-np.arange(0,space)*dx/args.mutationsigma)/args.mutationsigma*args.mutationrate*dx
xm4um = x-4.*u - args.mutationrate
speedprefactor = .5*args.wavespeed/dx

if args.c2file != None:
  try:
    c2data = np.genfromtxt(args.c2file)
  except:
    print >> sys.stderr,"could not open c2file"
    exit(1)
  c2 = c2data[:,2]
  c2.reshape((space,space))
elif args.cfile != None:
  try:
    cdata = np.genfromtxt(args.cfile)
  except:
    print >> sys.stderr,"could not open cfile"
    exit(1)
  c1 = cdata[:,1]
  c2 = np.outer(c1,c1)
else:
  c2 = np.outer(np.exp(-x*x/(2*args.wavespeed)),np.exp(-x*x/(2*args.wavespeed)))

final = opt.newton_krylov(f,c2,iter=args.maxSteps)

rescale = np.dot(np.dot(final,u),u)
final /= rescale

for i in range(space):
  for j in range(space):
    print x[i],x[j],final[i][j]
  print
