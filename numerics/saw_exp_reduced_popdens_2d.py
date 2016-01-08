#!/usr/bin/env python

import numpy as np
import argparse
import math,sys


def derivativeu(ff):
  n = 2*ff[-1]-ff[-2]
  if ff[1] > 0:
    l = ff[0]*ff[0]/ff[1]
  else:
    l = 0.
  return 0.5/args.dx*(np.diff(np.append(ff,n))+np.diff(np.insert(ff,0,l)))


def derivative2(ff,direction):
  if direction == 0:
    ff[ff[:,-1]<1e-200] = 1e-200
    ff[ff[:,-2]<1e-200] = 1e-200
    ff[ff[:,0]<1e-200] = 1e-200
    ff[ff[:,1]<1e-200] = 1e-200

    n = [ff[:,-1]*ff[:,-1]/ff[:,-2]]
    l = [ff[:,0]*ff[:,0]/ff[:,1]]
    
  else:
    ff[ff[-1,:]<1e-200] = 1e-200
    ff[ff[-2,:]<1e-200] = 1e-200
    ff[ff[0,:]<1e-200] = 1e-200
    ff[ff[1,:]<1e-200] = 1e-200

    n = np.transpose([ff[-1,:]*ff[-1,:]/ff[-2,:]])
    l = np.transpose([ff[0,:]*ff[0,:]/ff[1,:]])
    
  d = 0.5/dx * ( np.diff(np.append(l,ff,axis=direction),axis=direction) + np.diff(np.append(ff,n,axis=direction),axis=direction) )
  
  return d
    


parser = argparse.ArgumentParser()

parser.add_argument("-v","--speed",type=float,default=2)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-u","--ufile",required=True)

parser.add_argument("-S","--maxsteps",type=int,default=1000)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-O","--outputstep",type=int,default=1000)

parser.add_argument("-i","--infile")
parser.add_argument("-c","--c1file")


args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
  x = udata[:,0]
  u = 2.*udata[:,1]/3.
  space = len(x)
  space0 = (x*x).argmin()
  dx = x[1] - x[0]
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)


try:
  c2data = np.genfromtxt(args.infile)
  c2 = (c2data[:,2]).reshape((space,space))
  speed = args.speed
  mutationrate = args.mutationrate
except:
  try:
    c1data = np.genfromtxt(args.c1file)
    c1 = c1data[:,1]
    c2 = np.outer(c1,c1)
    normalize = np.dot(np.dot(c2,u),u)*dx*dx
    c2 /= normalize
    speed = args.speed
    mutationrate = args.mutationrate
  except:
    print >> sys.stderr,"starting from scratch"
    speed = args.speed
    mutationrate = args.mutationrate
    c1 = np.exp(-x*x/(2*speed))
    c2 = np.outer(c1,c1)
    normalize = np.dot(np.dot(c2,u),u)*dx*dx
    c2 /= normalize


s  = x - 4*u
sx0 = np.repeat([s],space,axis=0)
sy0 = np.repeat(np.transpose([s]),space,axis=1)
sx1 = derivative2(sx0,0)
sy1 = derivative2(sy0,1)

ux0 = np.repeat([u],space,axis=0)
uy0 = np.repeat(np.transpose([u]),space,axis=1)

fc = sx0 + sx1 + sy0 + sy1 
fc += - 4*speed - (1+2/dx**2)*(ux0+uy0)

for o in range(args.maxsteps):
  cd00 = c2
  cd10 = derivative2(c2,0)
  cd01 = derivative2(c2,1)
  cd20 = derivative2(cd10,0)
  cd02 = derivative2(cd01,1)
  cd11 = derivative2(cd10,1)
  cd21 = derivative2(cd11,0)
  cd12 = derivative2(cd11,1)
  
  c1  = 2*np.dot(c2,u)
  c1m = 0.5*(np.delete(c1,0)+np.delete(c1,-1))
  diagterm = np.diag(c1)*(1.+2/dx**2) - np.diag(c1m,k=1)/dx**2 - np.diag(c1m,k=-1)/dx**2
  
  f  = cd00*(sx0+sx1+sy0+sy1) 
  f += cd10*(speed - mutationrate + sx0 + sy0 + sy1)
  f += cd01*(speed - mutationrate + sy0 + sx0 + sx1)
  f += cd11*(2*speed - 2*mutationrate + sx0 + sy0)
  f += speed*(cd20+cd21+cd02+cd12)
  f += diagterm
  
  c2 -= args.alpha* f / fc
  #c2[c2<0] = 0

  normalize = np.dot(np.dot(c2,u),u)*dx**2
  
  c2 /= normalize

for i in range(space):
  for j in range(space):
    print x[i],x[j],c2[i,j]
  print

