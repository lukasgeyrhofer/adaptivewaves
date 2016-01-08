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
m = np.log(abs(u))
space = len(x)
space0 = np.argmin(x*x)
mutationmatrix = np.zeros([space,space])
mutationcorrection = np.zeros(space)

for i in range(space):
  for j in range(i,space):
    mutationmatrix[i,j] = np.exp((i-j)*dx/mutationsigma)/mutationsigma*mutationrate*dx
  if x[i] > 0:
    mutationcorrection[i] = mutationrate*(mutationsigma+(space-space0)*dx)*np.exp(-(space-i)*dx/mutationsigma)/x[i]


for o in range(maxiterations):
  print >> sys.stderr,o
  jacobian = np.zeros([space,space])
  fm = np.zeros(space)

  for i in range(space):
    fm[i] = x[i] - 2.*np.exp(m[i])-mutationrate
    jacobian[i,i] = -2.*m[i]*np.exp(m[i])
    if i>0:
      fm[i] += 0.5*wavespeed/dx*m[i-1]
      jacobian[i,i-1] += 0.5*wavespeed/dx
    if i<space-1:
      fm[i] -= 0.5*wavespeed/dx*m[i+1]
      jacobian[i,i+1] -= 0.5*wavespeed/dx
    else:
      #print >> sys.stderr,i
      fm[i] -= 0.5*wavespeed/dx*np.exp((space-space0)*dx*0.5)
    for j in range(space):
      jacobian[i,j] += mutationmatrix[i,j]*np.exp(m[j]+m[i])*m[j]
      fm[i] += mutationmatrix[i,j]*np.exp(m[j]+m[i])
      
  fm += mutationcorrection
  
  
  ijac=np.linalg.inv(jacobian)

  for i in range(space):
    for j in range(space):
      m[i] -= ijac[i,j]*fm[j]
    #if u[i]<0:u[i]=0.

for i in range(space):
  print "%10lf\t%20.14e"%(x[i],np.exp(m[i]))
  
  
  
  