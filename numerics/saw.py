#!/usr/bin/env python


import numpy as np
import sys
import math


mutationrate = 1e-5
mutationsigma = 1e-3
wavespeed = 1e-5

maxiterations = 10

try:
  initialdata = np.genfromtxt(sys.argv[1])
except:
  exit(1)


x = initialdata[:,0]
dx = x[1]-x[0]
w = np.sqrt(abs(initialdata[:,1]))
space = len(x)

mutationmatrix = np.zeros([space,space])

for i in range(space):
  for j in range(i,space):
    mutationmatrix[i,j] = np.exp(-(i-j)*dx/mutationsigma)/mutationsigma*mutationrate*dx
  mutationmatrix[i,i] -= mutationrate

for m in range(maxiterations):
  jacobian = np.zeros([space,space])
  fw = np.zeros(space)

  for i in range(space):
    fw[i] = x[i]*w[i]**2. - 2.*w[i]**4.
    jacobian[i,i] = 2.*x[i]*w[i]-8.*w[i]**3.
    if i>0:
      fw[i] += wavespeed/dx*w[i-1]*w[i]
      jacobian[i,i]   += wavespeed/dx*w[i-1]
      jacobian[i,i-1] += wavespeed/dx*w[i]
    if i<space-1:
      fw[i] -= wavespeed/dx*w[i+1]*w[i]
      jacobian[i,i]   -= wavespeed/dx*w[i+1]
      jacobian[i,i+1] -= wavespeed/dx*w[i]
    for j in range(space):
      jacobian[i,j] += 2.*mutationmatrix[i,j]*w[j]
      fw[i] += mutationmatrix[i,j]*w[j]**2
      


  ijac=np.linalg.inv(jacobian)

  w -= ijac.dot(fw)

for i in range(space):
  print "%10lf\t%20.14e"%(x[i],w[i]**2)
  
  
  
  