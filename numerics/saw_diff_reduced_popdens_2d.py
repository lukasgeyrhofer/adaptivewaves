#!/usr/bin/env python

import numpy as np
import sys,math
import argparse


def derivative(f,order,direction):
  if direction == 0:
    a1 = f[:,1]
    a2 = f[:,2]
    b1 = f[:,-1]
    b2 = f[:,-2]
  else:
    a1 = f[1,:]
    a2 = f[2,:]
    b1 = f[-1,:]
    b2 = f[-2,:]
  n = np.zeros(space)
  l = np.zeros(space)
  n[b2>0] = b1[b2>0]*b1[b2>0]/b2[b2>0]
  l[a2>0] = a1[a2>0]*a1[a2>0]/a2[b2>0]
  if order == 1:
    return 0.5/dx*(np.diff(np.append(f,n,axis=direction),axis=direction)+np.diff(np.insert(f,0,l,axis=direction),axis=direction))
  if order == 2:
    return 1./(dx*dx)*np.diff(np.append(np.insert(f,0,l,axis=direction),n,axis=direction),n=2,axis=direction)


parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-C","--c2file")
parser.add_argument("-c","--cfile")
parser.add_argument("-v","--speed",type=float,default=1.)

parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-O","--outputstep",type=int,default=100)
parser.add_argument("-a","--alpha",type=float,default=1.)
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
  cdata = np.genfromtxt(args.c2file)
  cx = cdata[:,0]
  cspace = len(cx)
  cspace0 = (cx*cx).argmin()
  cdx = cx[1]-cx[0]
  assert (space == cspace) and (space0 == cspace0) and (dx - cdx < 1e-6), "lattice parameters do not match!"
  c2 = cdata[:,2].reshape((space,space))
except:
  try:
    cdata = np.genfromtxt(args.cfile)
    cx = cdata[:,0]
    c1 = cdata[:,1]
    cspace = len(cx)
    cspace0 = (cx*cx).argmin()
    cdx = cx[1]-cx[0]
    c2 = np.outer(c1,c1)
  except:
    print >> sys.stderr,"# starting from scratch"
    c1 = np.exp(-x*x/(2*speed))
    c2 = np.outer(c1,c1)

normalize = np.dot(u,np.dot(u,c2))*dx*dx
c2 /= normalize

xu2 = np.zeros((space,space))
for i in range(space):
  for j in range(space):
    xu2[i,j] = x[i]-4.*u[i] +x[j]- 4.*u[j]

for i in range(args.maxsteps):
  cdx2  = derivative(c2,2,0)
  cdx1  = derivative(c2,1,0)
  cdy2  = derivative(c2,2,1)
  cdy1  = derivative(c2,1,1)
  
  c1 = np.dot(u,c2)         # *dx
  diagterm = 2.*np.diag(c1) # /dx
                            # a "discrete Dirac-delta" has a finite height at 0!

  f  = cdx2 + speed*cdx1 + cdy2 + speed*cdy1 + xu2*c2 + diagterm
  fu = -4./(dx*dx) + xu2
  
  c2 -= args.alpha*f/fu
  
  normalize = np.dot(u,np.dot(u,c2))*dx*dx
  c /= normalize
  
  if args.outputstep>0 and (i%args.outputstep == 0):
    print >> sys.stderr,i,np.dot(u,np.dot(u*u,c2)*dx)*dx

for i in range(space):
  for j in range(space):
    print "%lf %lf %.14e %.14e %.14e"%(x[i],x[j],c2[i,j],u[i],u[j])
  print
