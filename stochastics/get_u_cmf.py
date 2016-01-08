#!/usr/bin/env python

import numpy as np
import sys,math
import argparse


argparser = argparse.ArgumentParser()
argparser.add_argument("-u","--ufile")
argparser.add_argument("-U","--densustarlatticeratio",type=int,default=1)
argparser.add_argument("-t","--time",type=float,default=0.)
argparser.add_argument("-v","--wavespeed",type=float,default=0.)
argparser.add_argument("-H","--nonhistomode",action='store_false',default=True)
args = argparser.parse_args()

try:
  data = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)
  
x = data[:,0]
u_read = data[:,1]

space_u = len(x)
dx_u = x[1]-x[0]
space0_u = (x*x).argmin()

assert space_u%args.densustarlatticeratio == 0,"lattice ratio mismatch"

space_n = int(space_u/args.densustarlatticeratio)
space0_n = int(space0_u/args.densustarlatticeratio)
dx_n = dx_u * args.densustarlatticeratio


fitnessoffset = args.time * args.wavespeed

baseoffset = int(np.floor(fitnessoffset/dx_n))
uoffset = int(np.floor(fitnessoffset/dx_u)) - baseoffset*args.densustarlatticeratio

fracoffset = (fitnessoffset - uoffset*dx_u - baseoffset*dx_n)/dx_u

xn = (np.arange(space_n)-space0_n + baseoffset)*dx_n

u = np.zeros(space_n)

if args.densustarlatticeratio > 1:
  i = 0
  j = 0
  while i<space_n and j<space_u :
    if((j+uoffset)%args.densustarlatticeratio == 0):
      u[i] = (fracoffset) * u_read[j]
    elif((j+uoffset)%args.densustarlatticeratio == args.densustarlatticeratio - 1):
      u[i] += (1.-fracoffset) * u_read[j]
      u[i] /= args.densustarlatticeratio
      i+=1
    else:
      u[i] += u_read[j]
    j+=1
else:
  for i in range(1,space_n):
    u[i] = fracoffset * u_read[i-1] + (1.-fracoffset)*u_read[i]


for i in range(space_n):
  print xn[i],u[i]
  if args.nonhistomode:
    print xn[i]+dx_n,u[i]
