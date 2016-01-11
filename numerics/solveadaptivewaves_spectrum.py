#!/usr/bin/env python

import numpy as np
from scipy import linalg
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-M","--mutationmodel",choices=("diff","exp"),default="diff")
parser.add_argument("-G","--growthterm",choices=("selection","step"),default="selection")
args = parser.parse_args()

try:
    udata = np.genfromtxt(args.ufile)
except:
    # an error occurred during loading. aborting
    print >> sys.stderr,"could not open file"
    exit(1)

# generate necessary variables and profiles
x = udata[:,0]
u = udata[:,1]
space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]
speed  = args.speed
mutationrate = args.mutationrate



if args.growthterm == "selection":
    growth = x
elif args.growthterm == "step":
    growth = np.zeros(space)
    growth[x>0] = 1.

# construct matrix f
# growth in tuned models
f  = np.diag(growth-2*u)
# wavespeed
f += 0.5*speed/dx*(np.diag(np.ones(space-1),k=1) - np.diag(np.ones(space-1),k=-1))

# mutation model
if args.mutationmodel == "diff":
    f += (np.diag(np.ones(space-1),k=-1) -2*np.diag(np.ones(space),k=0) + np.diag(np.ones(space-1),k=1))/(dx*dx)
if args.mutationmodel == "exp":
    mut_outflow = 0
    for i in range(space):
	f += mutationrate*np.exp(-i)*(1-np.exp(-dx))*np.diag(np.ones(space-i),k=i) 
	mut_outflow += np.exp(-i)
    f -= mutationrate*mut_outflow*(1-np.exp(-dx))*np.diag(np.ones(space))

# compute eigenvalues
ev = linalg.eigvals(f)

# output
for i in range(len(ev)):
    print i,np.real(ev[i]),np.imag(ev[i])


