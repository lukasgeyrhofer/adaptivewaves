#!/usr/bin/env python

import argparse
import numpy as np
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-c","--cfile")
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=0.)
parser.add_argument("-T","--maxtime",type=int,default=2000)
parser.add_argument("-t","--dt",type=float,default=1e-2)
args = parser.parse_args()


try:
  udata = np.genfromtxt(args.ufile)
  cdata = np.genfromtxt(args.cfile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

x     = udata[:,0]
u     = udata[:,1]
c2raw = cdata[:,2]

space = len(x)
dx    = x[1] - x[0]

c2 = c2raw.reshape((space,space))

c1 = np.dot(c2,u)*dx

print np.dot(c1,u)*dx

#for i in range(space):
  #print >> sys.stderr,x[i],u[i],c1[i]

#exit(1)

c2y = np.zeros((space,space))

for i in range(space):
  cd = np.diag(c2,k=i)
  #print i,len(cd)
  if i==0:
    y = cd[:]
  elif (i%2 == 0) and (i>0):
    z = np.zeros(i/2)
    y = np.concatenate((z,cd,z))
  else:
    zm = np.zeros((i-1)/2)
    zp = np.zeros((i+1)/2)
    y = 0.5*(np.concatenate((zm,cd,zp)) + np.concatenate((zp,cd,zm)))
  c2y[i,:] = y[:]


ps   = np.sum(c1)
xxc  = np.sum(x*x*c1)/ps
xc   = np.sum(x*c1)/ps

fv   = xxc-xc*xc

print fv,xxc,xc

ps2  = np.sum(np.sum(c2))

vfv = np.dot((x-xc)*(x-xc),np.dot((x-xc)*(x-xc),c2))/ps2

print vfv
