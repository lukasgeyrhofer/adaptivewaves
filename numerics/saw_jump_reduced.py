#!/usr/bin/env python

import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-v","--speed",type=float,default=1.)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-s","--space",type=int,default=1000)
parser.add_argument("-z","--space0",type=int,default=500)
parser.add_argument("-d","--dx",type=float,default=5e-2)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-S","--maxsteps",type=int,default=20)
args = parser.parse_args()


try:
  data = np.genfromtxt(args.infile)
  x = data[:,0]
  u = data[:,1]
  space = len(x)
  space0 = (x*x).argmin()
  dx = x[1] - x[0]
except:
  x = (np.arange(args.space)-args.space0)*args.dx
  space = args.space
  space0 = args.space0
  dx = args.dx
  u = np.exp(x)*0.5*dx
  u[x>0] = 0.5*x[x>0]

fu = np.zeros((space,space))

fu += np.diagflat(x-args.mutationrate)
fu += np.diagflat(np.ones(space-int(1./dx)),k=int(1./dx))*args.mutationrate
jacbase = fu[:,:]

#fu -= args.speed*(np.diagflat(np.ones(space-2),k=-2) - 8.*np.diagflat(np.ones(space-1),k=-1) + 8.*np.diagflat(np.ones(space-1),k=1) - np.diagflat(np.ones(space-2),k=2))/(12.*dx)

fu -= (-np.diagflat(np.ones(space-1),k=-1) + np.diagflat(np.ones(space-1),k=1))/(2.*dx)

# speed correction, assume linear extrapolation of u
fu[space-1,space-1] += 2/(2.*dx)
fu[space-1,space-2] -= 1/(2.*dx)

for o in range(args.maxsteps):
  print >> sys.stderr,o
  unew = np.dot(fu,u)-2*u*u
  unew[space-1] -= args.speed/(2*dx)*(2*u[space-1]-u[space-2])
  if u[1] > 1e-300:
    unew[0] -= args.speed/(2*dx)*(u[0]**2/u[1])
  jac = jacbase - 4*np.diagflat(u)
  ijac = np.linalg.inv(jac)
  u -= args.alpha*np.dot(ijac,unew)
  u[u<1e-300] = 1e-300


for i in range(space):
  print x[i],u[i]