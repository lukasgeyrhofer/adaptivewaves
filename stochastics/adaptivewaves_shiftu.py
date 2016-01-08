#!/usr/bin/env python


import numpy as np
import math,sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-v","--speed",type=float,default=0)
parser.add_argument("-t","--time",type=float,default=0)
parser.add_argument("-U","--latticeratio",type=int,default=1)
parser.add_argument("-H","--printhistomode",action="store_true",default=False)
args = parser.parse_args()


try:
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)


xu = udata[:,0]
uu = udata[:,1]

spaceu  = len(xu)
space0u = (xu*xu).argmin()
dxu     = xu[1] - xu[0]

space = int(spaceu/args.latticeratio)
space0 = int(spaceu/args.latticeratio)
dx = dxu*args.latticeratio

u = np.zeros(space+1)

fracshift = args.time*args.speed/dxu - int(args.time*args.speed/dxu)
#print fracshift
if args.latticeratio > 1:
  i = 0
  j = 0
  while i<space and j<spaceu:
    if j%args.latticeratio == 0:
      u[i] += (1.-fracshift)*uu[j]
      u[i] /= args.latticeratio
      i += 1
      u[i] = fracshift * uu[j]
    else:
      u[i] += uu[j]
    j += 1
else:
  for i in range(1,space):
    u[i] = fracshift*uu[i-1] + (1.-fracshift)*uu[i]

for i in range(space):
  print u[i]
  if args.printhistomode:
    print u[i]