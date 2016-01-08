#!/usr/bin/env python

import numpy as np
import sys,math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-v","--wavespeed",type=float,default=1e-5)
parser.add_argument("-D","--mutationrate",type=float,default=1e-5)
parser.add_argument("-m","--mutationmodel",choices=["diffusion","jumps","expdecay","dexp"],default="expdecay")
parser.add_argument("-s","--mutationsigma",type=float,default=1e-3)
parser.add_argument("-K","--mutationsigman",type=float,default=1e-3)
parser.add_argument("-k","--ratiodelmut",type=float,default=0.5)
parser.add_argument("-j","--jumpwidth",type=int,default=100)
parser.add_argument("-u","--ufile")
args = parser.parse_args()

try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)
  
x = udata[:,0]
u = udata[:,1]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

fc = np.zeros(shape=(space+1,space))

b = np.zeros(space+1)
b[space] = 1



for i in range(space):
  
  # selection and constraint
  fc[i,i] = x[i]-2*u[i]
  # constraint again
  fc[space,i] = u[i]*dx

  # comoving frame
  if i>0:
    fc[i,i-1] -= 0.5*args.wavespeed/dx
  if i<space-1:
    fc[i,i+1] += 0.5*args.wavespeed/dx

  # mutations 
  if args.mutationmodel == "expdecay":
    fc[i,i] -= args.mutationrate
    for j in range(i):
      fc[i,j] += np.exp((j-i)*dx/args.mutationsigma)*args.mutationrate/args.mutationsigma*dx
  elif args.mutationmodel == "dexp":
    for j in range(i):
      fc[i,j] += np.exp((j-i)*dx/args.mutationsigma)*args.mutationrate*(1-np.exp(-dx/args.mutationsigma))
    fc[i,i] -= args.mutationrate*(np.exp(-dx/args.mutationsigma) + args.ratiodelmut*np.exp(-dx/args.mutationsigman))
    for j in range(i+1,space):
      fc[i,j] += np.exp((i-j)*dx/args.mutationsigman)*args.mutationrate*args.ratiodelmut*(1-np.exp(-dx/args.mutationsigman))
  elif args.mutationmodel == "diffusion":
    fc[i,i] -= 2*args.mutationrate/(dx*dx)
    if i>0:
      fc[i,i-1] += args.mutationrate/(dx*dx)
    if i<space-1:
      fc[i,i+1] += args.mutationrate/(dx*dx)
  elif args.mutationmodel == "jumps":
    fc[i,i] -= args.mutationrate
    if i >= args.jumpwidth:
      fc[i,i-args.jumpwidth] += args.mutationrate


c = np.dot(np.linalg.pinv(fc),b)
c[c<0]=0

rescale = np.dot(c,u)*dx
c /= rescale


for i in range(space):
  print x[i],c[i]