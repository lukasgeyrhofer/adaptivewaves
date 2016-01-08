#!/usr/bin/env python


import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
#parser.add_argument("-s","--space",type=int,default=200)
#parser.add_argument("-z","--zero",type=int,default=50)
#parser.add_argument("-d","--dx",type=float,default=1e-3)
parser.add_argument("-D","--mutationrate",type=float,default=1e-5)
parser.add_argument("-M","--mutationsigma",type=float,default=1e-3)
parser.add_argument("-v","--wavespeed",type=float,default=1e-5)
parser.add_argument("-S","--maxiterations",type=int,default=10)

args = parser.parse_args()

mutationrate = args.mutationrate
mutationsigma = args.mutationsigma
wavespeed = args.wavespeed
maxiterations = args.maxiterations

try:
  initialdata = np.genfromtxt(args.inputfile)
except:
  exit(1)

x = initialdata[:,0]
dx = x[1]-x[0]
u = initialdata[:,1]
w = np.sqrt(abs(u))
space = len(x)
space0 = np.argmin(x*x)
mutationmatrix = np.zeros([space,space])
mutationcorrection = np.zeros(space)


for i in range(space):
  for j in range(i,space):
    mutationmatrix[i,j] = np.exp((i-j)*dx/mutationsigma)/mutationsigma*mutationrate*dx
  mutationmatrix[i,i] -= mutationrate
  mutationcorrection[i] = 0.5*mutationrate*(mutationsigma+(space-space0)*dx)*np.exp(-(space-i)*dx/mutationsigma)


for m in range(maxiterations):
  print >> sys.stderr,m
  jacobian = np.zeros([space,space])
  fw = np.zeros(space)

  for i in range(space):
    fw[i] = x[i]*w[i]**2. - 2.*w[i]**4.
    jacobian[i,i] = 2.*x[i]*w[i]-8.*w[i]**3.
    if i>0:
      fw[i] += wavespeed/dx*w[i-1]*w[i]
      jacobian[i,i-1] += wavespeed/dx*w[i]
      jacobian[i,i]   += wavespeed/dx*w[i-1]
    if i<space-1:
      fw[i] -= wavespeed/dx*w[i+1]*w[i]
      jacobian[i,i+1] -= wavespeed/dx*w[i]
      jacobian[i,i]   -= wavespeed/dx*w[i+1]
    else:
      fw[i] -= wavespeed/dx*np.sqrt((space-space0)*dx*0.5)*w[i]
      jacobian[i,i]   -= wavespeed/dx*np.sqrt((space-space0)*dx*0.5)
    for j in range(space):
      jacobian[i,j] += 2.*mutationmatrix[i,j]*w[j]
      fw[i] += mutationmatrix[i,j]*w[j]**2.
      
  fw += mutationcorrection
  
  
  ijac=np.linalg.inv(jacobian)

  for i in range(space):
    for j in range(space):
      w[i] -= ijac[i,j]*fw[j]
    #if w[i] < 0: w[i]=0.

for i in range(space):
  print "%10lf\t%20.14e\t%20.14e"%(x[i],w[i]**2,w[i])
  
