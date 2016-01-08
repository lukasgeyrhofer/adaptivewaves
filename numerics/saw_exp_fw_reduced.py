#!/usr/bin/env python

import numpy as np
import math,sys
import argparse

def derivative(f,order):
  n = 2*f[-1] - f[-2] # linear increase
  l = f[0]*f[0]/f[1]  # exponential decrease
  if order == 1:
    return 0.5/dx*(np.diff(np.append(f,n))+np.diff(np.insert(f,0,l)))
  if order == 2:
    return 1./(dx*dx)*np.diff(np.append(np.insert(f,0,l),n),n=2)

parser = argparse.ArgumentParser()
parser.add_argument("-s","--space",type=int,default=4000)
parser.add_argument("-z","--space0",type=int,default=1600)
parser.add_argument("-d","--dx",type=float,default=5e-2)
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-S","--maxsteps",type=int,default=1000)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-i","--infile")
parser.add_argument("-R","--redblack",action="store_true",default=False)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.infile)
  x = udata[:,0]
  u = udata[:,1]
  space = len(x)
  space0 = (x*x).argmin()
  dx = x[1]-x[0]
  speed = args.speed
  mutationrate = args.mutationrate
except:
  print >> sys.stderr,"starting from scratch"
  space=args.space
  space0=args.space0
  dx = args.dx
  speed = args.speed
  mutationrate= args.mutationrate
  x = (np.arange(space)-space0)*dx
  u = np.zeros(space)
  u[x< 0] = np.exp((speed-mutationrate)/speed*x[x<0])
  u[x>=0] = 1

#assert speed > mutationrate,"note that condition (speed > mutationrate) has to be fulfilled! check your parameters!"


growth            = np.zeros(space)
growth[x>=0]      = 1.
deltagrowth       = np.zeros(space)
deltagrowth[x==0] = 1./dx

for i in range(args.maxsteps):
  u2  = derivative(u,2)
  u1  = derivative(u,1)
  uu1 = derivative(u*u,1)
    
  f  = speed*u2 - (speed - mutationrate +  growth )*u1 + (growth - deltagrowth)*u - u*u + uu1
  fu = -2.*speed/(dx*dx) + (growth-deltagrowth) - 2*u + 2*u1
  
  if args.redblack:
    # update only every second element, start with 1 if step i is odd, start with 0 if step i is even
    f[i%2::2] = 0
  
  u -= args.alpha*f/fu

for i in range(space):
  print "%lf %.14e"%(x[i],u[i])
