#!/usr/bin/env python

import numpy as np
import sys,math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-i","--infile")
parser.add_argument("-v","--speed",type=float,default=1.)

parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-O","--outputstep",type=int,default=100)
parser.add_argument("-R","--redblack",action="store_true",default=False)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-Q","--quiet",action="store_true",default=False)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)

x = udata[:,0]
u = udata[:,1]
space = len(x)
space0 = (x*x).argmin()
dx = x[1]-x[0]
speed=args.speed


try:
  cdata = np.genfromtxt(args.infile)
  cx = cdata[:,0]
  cspace = len(cx)
  cspace0 = (cx*cx).argmin()
  cdx = cx[1]-cx[0]
  assert (space == cspace) and (space0 == cspace0) and (dx - cdx < 1e-6), "lattice parameters do not match!"
  c = cdata[:,1]
except:
  print >> sys.stderr,"# starting from scratch"
  c = np.exp(-x*x/(2*speed))

c /= np.dot(u,c)*dx

coeffm = 1./dx**2 - 0.5*speed/dx
coeff0 = -2./dx**2 + x - 2.*u
coeffp = 1./dx**2 + 0.5*speed/dx

for i in range(args.maxsteps):
  cm  = np.concatenate((np.array([c[0]**2/c[1] if c[1]>0 else 0]),c[:-1]))
  cp  = np.concatenate((c[1:],np.array([c[-1]**2/c[-2] if c[-2]>0 else 0])))
    
  f  = coeffm * cm + coeff0 * c + coeffp * cp
  fu = coeff0
  
  if args.redblack:
    # update only every second element, start with 1 if step i is odd, start with 0 if step i is even
    f[i%2::2] = 0
  
  c -= args.alpha*f/fu
  c[c<0] = 0
  
  c /= np.dot(u,c)*dx
  
  if not args.quiet:
    if args.outputstep>0 and (i%args.outputstep == 0):
      print >> sys.stderr,i,np.dot(u*u,c)*dx

for i in range(space):
  print "%lf %.14e %.14e"%(x[i],c[i],u[i])
