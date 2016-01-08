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
  fu = np.zeros(space)

  for i in range(space):
    fu[i] = x[i]*u[i] - 2.*u[i]**2.
    jacobian[i,i] = x[i]-4.*u[i]
    if i>0:
      fu[i] += 0.5*wavespeed/dx*u[i-1]
      jacobian[i,i-1] += 0.5*wavespeed/dx
    if i<space-1:
      fu[i] -= 0.5*wavespeed/dx*u[i+1]
      jacobian[i,i+1] -= 0.5*wavespeed/dx
    else:
      #print >> sys.stderr,i
      fu[i] -= wavespeed/dx*(space-space0)*dx*0.25
    for j in range(space):
      jacobian[i,j] += mutationmatrix[i,j]
      fu[i] += mutationmatrix[i,j]*u[j]
      
  fu += mutationcorrection
  
  
  ijac=np.linalg.inv(jacobian)

  for i in range(space):
    for j in range(space):
      u[i] -= ijac[i,j]*fu[j]
    #if u[i]<0:u[i]=0.

for i in range(space):
  print "%10lf\t%20.14e"%(x[i],u[i])
  
  
  
  