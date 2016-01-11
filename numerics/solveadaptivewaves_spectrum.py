#!/usr/bin/env python

import numpy as np
from scipy import linalg
import argparse
import sys,math

parser = argparse.ArgumentParser()
#parser.add_argument("-s","--space",type=int,default=2000)
#parser.add_argument("-z","--space0",type=int,default=1000)
#parser.add_argument("-d","--dx",type=float,default=5e-2)

parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-M","--mutationmodel",choices=("diff","exp"),default="diff")
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-G","--growthterm",choices=("selection","step"),default="selection")
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
speed = args.speed
mutationrate = args.mutationrate

if args.growthterm == "selection":
    growth = x
elif args.growthterm == "step":
    growth = np.zeros(space)
    growth[x>0] = 1.

o0 = np.ones(space)
o1 = np.ones(space-1)
f  = np.diag(growth-2*u)
f += 0.5*speed/dx*(np.diag(o1,k=1) - np.diag(o1,k=-1))

if args.mutationmodel == "diff":
    f += (np.diag(o1,k=-1) -2*np.diag(o0,k=0) + np.diag(o1,k=1))/(dx*dx)
if args.mutationmodel == "exp":
    mut_outflow = 0
    for i in range(space):
	oi = np.ones(space-i)
	f += mutationrate*np.exp(-i)*np.diag(oi,k=i)
	mut_outflow += np.exp(-i)
    f -= mutationrate*mut_outflow*np.diag(o0)

ev = linalg.eigvals(f)
for i in range(len(ev)):
    print i,np.real(ev[i]),np.imag(ev[i])
