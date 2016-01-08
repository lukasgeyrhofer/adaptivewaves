#!/usr/bin/env python

import numpy as np
from scipy.special import airy
import sys,math
import argparse


def derivative(f,order):
  n = 2*f[-1] - f[-2]   # next point after lattice: linear increase
  if f[1]>0:
    l = f[0]*f[0]/f[1]  # last point before lattice: exponential decrease
  else:
    l = 0
  if order == 1:
    return 0.5/dx*(np.diff(np.append(f,n))+np.diff(np.insert(f,0,l)))
  if order == 2:
    return 1./(dx*dx)*np.diff(np.append(np.insert(f,0,l),n),n=2)


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-v","--speed",type=float,default=1.)

parser.add_argument("-s","--space",type=int,default=2000)
parser.add_argument("-z","--zero",type=int,default=1000)
parser.add_argument("-d","--dx",type=float,default=5e-2)

parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-R","--redblack",action="store_true",default=False)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-B","--enforceboundaries",action="store_true",default=False)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.infile)
  x = udata[:,0]
  u = udata[:,1]
  space = len(x)
  space0 = (x*x).argmin()
  dx = x[1]-x[0]
  speed=args.speed
except:
  print >> sys.stderr,"starting from scratch"
  dx = args.dx
  space= args.space
  space0 = args.zero
  speed = args.speed
  x = (np.arange(space)-space0)*dx
  uairy =  np.exp(speed*x/2.)*airy(speed**2/4.-x)[0]
  u = uairy
  firstAiZero = min(int((2.3381-speed**2/4)/dx)+space0,space-1)
  idxmax = u[:firstAiZero].argmax()
  u = u/u[idxmax]*x[idxmax]/2
  u[idxmax:] = x[idxmax:]/2

if args.enforceboundaries:
  maxvalue = (space-space0)*dx

cm = np.exp(speed*dx/2.)/dx**2
cp = np.exp(-speed*dx/2.)/dx**2
c0 = - 2.*np.cosh(-speed*dx/2.)/dx**2+x

for i in range(args.maxsteps):

    um = np.concatenate((np.array([u[0]*u[0]/u[1]]) if u[1]>0 else np.zeros(1),u[:-1]))
    up = np.concatenate((u[1:],np.array([2*u[-1]-u[-2]])))
    f  = cm*um + c0*u + cp*up - 2.*u*u
    fu = c0 - 4*u
    
    if args.redblack:
        # update only every second element, start with 1 if step i is odd, start with 0 if step i is even
        f[i%2::2] = 0
    
    u -= args.alpha*f/fu
    
    if args.enforceboundaries:
        u[u<0] = 0
        u[u>maxvalue] = maxvalue

for i in range(space):
  print "%lf %.14e"%(x[i],u[i])
