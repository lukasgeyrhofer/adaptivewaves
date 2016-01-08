#!/usr/bin/env python

import numpy as np
import math,sys
import argparse


def update_constraint(cc,uu,ddxx):
  return cc/(ddxx*np.dot(cc,uu))

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-c","--cfile",default=None)
parser.add_argument("-v","--wavespeed",type=float,default=1e-6)
parser.add_argument("-s","--sigma",type=float,default=1e-3)
parser.add_argument("-m","--mutationrate",type=float,default=1e-6)
parser.add_argument("-S","--maxsteps",type=int,default=1000)
parser.add_argument("-e","--epsilon",type=float,default=1e-2)
parser.add_argument("-C","--applyconstraint",default=False,action='store_true')
parser.add_argument("-I","--showintegrated",default=False,action='store_true')
parser.add_argument("-O","--outputstep",type=int,default=100)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

x = udata[:,0]
u = udata[:,1]

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]

if args.cfile != None:
  try:
    cdata = np.genfromtxt(args.cfile)
  except:
    print >> sys.stderr,"could not open file"
    exit(1)
  c = cdata[:,1]
  assert space == len(c),"lattice does not match"
else:
  c = np.exp(-x*x/(2*(args.wavespeed + args.mutationrate*args.sigma)))

c = update_constraint(c,u,dx)

lmatrix = np.diagflat(x-2*u-args.mutationrate*np.exp(-dx/args.sigma))
for i in range(space):
  if i>0:
    lmatrix[i,i-1] -= 0.5*args.wavespeed/dx
  if i<space-1:
    lmatrix[i,i+1] += 0.5*args.wavespeed/dx
  for j in range(i):
    lmatrix[i,i-j] += args.mutationrate * (1. - np.exp(-dx/args.sigma)) * np.exp(-j*dx/args.sigma)

for i in range(space):
  print "%10.6lf %10.6lf %18.10e"%(0,x[i],c[i])

for o in range(1,args.maxsteps):
  c+= args.epsilon*np.dot(lmatrix,c)
  if args.applyconstraint:
    c = update_constraint(c,u,dx)
  if o % args.outputstep == 0:
    for i in range(space):
      print "%10.6lf %10.6lf %18.10e"%(o*args.epsilon,x[i],c[i])
    print
    if args.showintegrated:
      ps = np.sum(c)
      xc = np.dot(x,c)/ps
      xxc = np.dot(x*x,c)/ps
      print >> sys.stderr,o*args.epsilon,ps*dx,xc,xxc-xc*xc


