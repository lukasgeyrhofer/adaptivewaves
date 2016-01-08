#!/usr/bin/env python


import numpy as np
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-m","--mutationrate",type=float,default=1e-5)
parser.add_argument("-s","--sigma",type=float,default=1e-3)
parser.add_argument("-v","--wavespeed",type=float,default=2e-6)
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)
  
u = udata[:,1]
x = udata[:,0]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]


newmat = np.zeros((space*space,space*space+1))
b = np.zeros(space*space+1)
b[-1] = 1/(dx*dx)

for i in range(space):
  for j in range(space):
    newmat[i*space+j,i*space+j] += x[i]+x[j]-2*args.mutationrate - 4*u[i] - 4*u[j]
    if i < space-1:
      newmat[i*space+j,(i+1)*space+j] +=  args.wavespeed/(2.*dx)
    if i > 0:
      newmat[i*space+j,(i-1)*space+j] += -args.wavespeed/(2.*dx)
    if j < space-1:
      newmat[i*space+j,i*space+j+1]   +=  args.wavespeed/(2.*dx)
    if j > 0:
      newmat[i*space+j,i*space+j-1]   += -args.wavespeed/(2.*dx)
    
    for k in range(min(i,j)):
      newmat[i*space+j,k*space+j] += dx/args.sigma * args.mutationrate * np.exp(-(i-k)*dx/args.sigma)
      newmat[i*space+j,i*space+k] += dx/args.sigma * args.mutationrate * np.exp(-(j-k)*dx/args.sigma)
    if i<j:
      for k in range(i,j):
	newmat[i*space+j,i*space+k] += dx/args.sigma * args.mutationrate * np.exp(-(j-k)*dx/args.sigma)
    elif j<i:
      for k in range(j,i):
	newmat[i*space+j,k*space+j] += dx/args.sigma * args.mutationrate * np.exp(-(i-k)*dx/args.sigma)
    else:
      for k in range(space):
	newmat[i*space+j,k*space+j] += u[k]/dx
	newmat[i*space+j,i*space+k] += u[k]/dx
    newmat[i*space+j,-1] = u[i]*u[j]

c = np.dot(np.linalg.pinv(newmat),b)
c[c<0]=0


for i in range(space):
  for j in range(space):
    print x[i],x[j],c[i*space+j]
  print
   