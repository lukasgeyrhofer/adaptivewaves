#!/usr/bin/env python

import numpy as np
import math,sys
import argparse

def derivative(ff,order):
  # boundary terms:
  # population density falls off exponentially, otherwise set to zero
  if ff[-2] > 0:
    n = ff[-1]*ff[-1]/ff[-2]
  else:
    n = 0.
  if ff[1] > 0:
    l = ff[0]*ff[0]/ff[1]
  else:
    l = 0.
  if order == 1:
    return 0.5/args.dx*(np.diff(np.append(ff,n))+np.diff(np.insert(ff,0,l)))
  if order == 2:
    return 1./(dx*dx)*np.diff(np.append(np.insert(ff,0,l),n),n=2)

def ders(ff):
  n = (2.*ff[-1]**-1-ff[-2]**-1)**-1  # effective selection (X-2u) = v/X (for large X)
  l = 2*ff[1]-ff[2]                   # effective selection (X-2u) = X   (for negative X)
  return 0.5/args.dx*(np.diff(np.append(ff,n))+np.diff(np.insert(ff,0,l)))

parser = argparse.ArgumentParser()

parser.add_argument("-s","--space",type=int,default=4000)
parser.add_argument("-z","--space0",type=int,default=1600)
parser.add_argument("-d","--dx",type=float,default=1/20.)
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-u","--ufile",required=True)

parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-O","--outputstep",type=int,default=1000)

parser.add_argument("-o","--outfile")
parser.add_argument("-i","--infile")

parser.add_argument("-q","--quiet",default=False,action="store_true")

args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
  x = udata[:,0]
  u = udata[:,1]
  space = len(x)
  space0 = (x*x).argmin()
  dx = x[1] - x[0]
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)

try:
  cdata = np.genfromtxt(args.infile)
  c = cdata[:,1]
  speed = args.speed
  mutationrate = args.mutationrate
except:
  print >> sys.stderr,"starting from scratch"
  speed = args.speed
  mutationrate= args.mutationrate
  c = np.exp(-x*x/(2*speed))
  xcrossover = (derivative(u,2)).argmin()
  c[xcrossover:] = np.exp(-(speed-mutationrate)*(x[xcrossover:]-x[xcrossover])/speed)*c[xcrossover]
  norm = np.dot(u,c)*dx
  c/=norm

assert speed > mutationrate,"note that (speed = variance + mutationrate) ! check your parameters!"

s = x-2*u
coeff2 = speed
coeff1 = speed - mutationrate + s
coeff0 = s + ders(s)

for i in range(args.maxsteps):
  c2  = derivative(c,2)
  c1  = derivative(c,1)
  
  f  = coeff2 * c2 + coeff1 * c1 + coeff0 * c
  fc = -2.*speed/(dx*dx) + coeff0
  
  c -= args.alpha*f/fc
  
  c /= np.dot(c,u)*dx
  
  if (not args.quiet) and (i%args.outputstep == 0):
    print >> sys.stderr,i,np.dot(u*u,c)

for i in range(space):
  print "%lf %.14e"%(x[i],c[i])

