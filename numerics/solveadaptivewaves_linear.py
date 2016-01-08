#!/usr/bin/env python

import numpy as np
from scipy.special import airy
import sys,math
import argparse


parser = argparse.ArgumentParser(description="Numerical solution for fixation probabilities (n=1) in traveling wave models of adaptation with diffusion and exponential mutation kernels")
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-v","--speed",type=float,default=1.)
parser.add_argument("-M","--mutationmodel",choices=("diff","exp"),default="diff")
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)

parser.add_argument("-s","--space",type=int,default=2000)
parser.add_argument("-z","--zero",type=int,default=1000)
parser.add_argument("-d","--dx",type=float,default=5e-2)

parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-R","--redblack",action="store_true",default=False)
parser.add_argument("-a","--alpha",type=float,default=1.)
parser.add_argument("-B","--enforceboundaries",action="store_true",default=False)
args = parser.parse_args()

try:
    udata = np.genfromtxt(args.infile)
    x = udata[:,0]
    u = udata[:,1]
    space        = len(x)
    space0       = (x*x).argmin()
    dx           = x[1]-x[0]
    speed        = args.speed
    mutationrate = args.mutationrate
except:
    print >> sys.stderr,"# starting from scratch"
    dx           = args.dx
    space        = args.space
    space0       = args.zero
    speed        = args.speed
    mutationrate = args.mutationrate
    x = (np.arange(space)-space0)*dx
    
    # use semianalytical approximations as starting conditions
    if args.mutationmodel == "diff":
        uairy =  np.exp(speed*x/2.)*airy(speed**2/4.-x)[0]
        u = uairy
        firstAiZero = min(int((2.3381-speed**2/4)/dx)+space0,space-1)
        idxmax = u[:firstAiZero].argmax()
        u = u/u[idxmax]*x[idxmax]/2
        u[idxmax:] = x[idxmax:]/2
    elif args.mutationmodel == "exp":
        popsize = np.exp(np.sqrt(2.*np.log(1./mutationrate)*speed))/mutationrate # inverting relation in Good et al. (2012)
        idxcrossover1 = ((x-speed)**2).argmin()
        idxcrossover2 = ((x-np.sqrt(2*speed*np.log(popsize*np.sqrt(speed))))**2).argmin()
        u = np.zeros(space)
        if idxcrossover1 < idxcrossover2:
            u[idxcrossover2:] = x[idxcrossover2:]/2 - speed/(x[idxcrossover2:]*2)
            u[idxcrossover1:idxcrossover2] = u[idxcrossover2]*np.exp(x[idxcrossover1:idxcrossover2]**2/(2*speed))
            u[:idxcrossover1] = np.exp(x[:idxcrossover1])
        else:
            u[idxcrossover2:] = x[idxcrossover2:]/2 - speed/(x[idxcrossover2:]*2)
            u[:idxcrossover2] = u[idxcrossover2]*np.exp(x[:idxcrossover2])

if args.enforceboundaries:
  maxvalue = (space-space0)*dx

if args.mutationmodel == "diff":
    coeff_prev =  1./(dx*dx) + 0.5*speed/dx
    coeff_next =  1./(dx*dx) - 0.5*speed/dx
    coeff_0    = -2./(dx*dx) + x
elif args.mutationmodel == "exp":
    coeff_prev = speed/(dx*dx) + 0.5*(speed - mutationrate + x)/dx
    coeff_next = speed/(dx*dx) - 0.5*(speed - mutationrate + x)/dx
    coeff_0    = -2*speed/(dx*dx) + x - 1

# coupled Newton-Raphson iterations for each lattice point
for i in range(args.maxsteps):
    # shift profile for terms with derivatives
    u_prev = np.concatenate((np.array([u[0]*u[0]/u[1]]) if u[1]>0 else np.zeros(1),u[:-1])) # exponential decay before lattice
    u_next = np.concatenate((u[1:],np.array([2*u[-1]-u[-2]])))                              # linear increase after lattice
    
    f  = coeff_prev * u_prev + coeff_0 * u + coeff_next * u_next - 2.*u*u
    fu = coeff_0 - 4*u
    
    # more nonlinear terms with exponential kernel
    if args.mutationmodel == "exp":
        f  += 2*(u_next - u_prev)*u/dx
        fu += 2*(u_next - u_prev)/dx
    
    if args.redblack:
        # update only every second element, start with 1 if step i is odd, start with 0 if step i is even
        f[i%2::2] = 0
    
    # NR step
    u -= args.alpha*f/fu
    
    if args.enforceboundaries:
        u[u<0] = 0
        u[u>maxvalue] = maxvalue


# output
if args.outfile == None:    fp = sys.stdout
else:                       fp = open(args.outfile,"w")
for i in range(space):
  print >> fp,"%lf %.14e"%(x[i],u[i])

fp.close()

