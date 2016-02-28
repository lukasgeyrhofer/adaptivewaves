#!/usr/bin/env python

import numpy as np
from scipy.special import airy
import sys,math
import argparse

parser = argparse.ArgumentParser(description="Numerical solution for fixation probabilities (n=1) in traveling wave models of asexual adaptation with diffusion and exponential mutation kernels (using reduced units)")

parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser_alg.add_argument("-i","--infile",help="Start iterations from profile specified in INFILE. If parameter is not used or errors occur during loading, start from scratch.")
parser_alg.add_argument("-o","--outfile",default=None,help="Write output to OUTFILE instead of stdout")
parser_alg.add_argument("-S","--maxsteps",type=int,default=10000,help="Number of iteration steps [default: 10000]")
parser_alg.add_argument("-R","--redblack",action="store_true",default=False,help="Iterate even or odd lattice points alternatively, helps with stability. Should be used for profiles of fixation probabilities [default: OFF]")
parser_alg.add_argument("-a","--alpha",type=float,default=1.,help="'Speed' of Newton-Raphson iteration: u/c -= alpha f/f'. Classical NR: alpha=1. Slower convergence but more stability for alpha<1 [default: 1]")
parser_alg.add_argument("-B","--enforceboundaries",action="store_true",default=False,help="Force solution to positive (and bounded) values [default: OFF]")

parser_lattice = parser.add_argument_group(description="####   Lattice parameters   ####")
parser_lattice.add_argument("-s","--space",type=int,default=2000,help="Number of lattice points [default: 2000]")
parser_lattice.add_argument("-z","--zero",type=int,default=1000,help="Position of Zero on lattice [default: 1000]")
parser_lattice.add_argument("-d","--dx",type=float,default=5e-2,help="Lattice spacing [default: 0.05]")

parser_params = parser.add_argument_group(description="####   Profile parameters (in reduced units)  ####")
parser_params.add_argument("-v","--speed",type=float,default=1.,help="Adaptation speed [default: 1]")
parser_params.add_argument("-M","--mutationmodel",choices=("diff","exp"),default=None,help="Mutation kernel")
parser_params.add_argument("-m","--mutationrate",type=float,default=None,help="Mutation rate used in exponential mutation kernel. Value is ignored for diffusion mutation kernel. [default: None]")
parser_params.add_argument("-G","--growthterm",choices=("selection","step"),default="selection",help="Growth is either given by linear gradient for adaptation (\"selection\") or as step function for Fisher waves (\"step\") [default: selection]") 

args = parser.parse_args()

if args.mutationmodel == None:
    if args.mutationrate == None:
        mutationmodel = "diff"
    else:
        mutationmodel = "exp"
else:
    mutationmodel = "diff"

try:
    udata = np.genfromtxt(args.infile)
    x = udata[:,0]
    u = udata[:,1]
    space        = len(x)
    space0       = (x*x).argmin()
    dx           = x[1]-x[0]
    speed        = args.speed
    if mutationmodel == "exp":
        mutationrate = args.mutationrate
except:
    # could not load file, either because input file was not provided (no -i option), did not exist or some other error while loading
    # thus, start from scratch
    print >> sys.stderr,"# starting from scratch"
    dx           = args.dx
    space        = args.space
    space0       = args.zero
    speed        = args.speed
    if mutationmodel == "exp":
        mutationrate = args.mutationrate
    x = (np.arange(space)-space0)*dx
    
    # use semianalytical approximations as starting conditions
    if args.growthterm == "selection":
        try:
            if mutationmodel == "diff":
                # solution is proportial to exp(vx/2)Airy(v^2/4-x) for small fitness (by neglecting nonlinear term)
                # Airy function is oscillating, however, thus find first zero of it, and approximate u by its
                # linear branch starting at the maximum of the asymptotic solution above (before its first zero)
                uairy =  np.exp(speed*x/2.)*airy(speed**2/4.-x)[0]
                u = uairy
                firstAiZero = min(int((2.3381-speed**2/4)/dx)+space0,space-1)
                idxmax = u[:firstAiZero].argmax()
                u = u/u[idxmax]*x[idxmax]/2
                u[idxmax:] = x[idxmax:]/2
            elif mutationmodel == "exp":
                # from semianalytic considerations we know u has three regimes (see [Geyrhofer, 2014], also defined in [Good et al., 2012])
                popsize = np.exp(np.sqrt(2.*np.log(1./mutationrate)*speed))/mutationrate          # inverting relation in [Good et al., 2012]
                idxcrossover1 = ((x-speed)**2).argmin()                                           # crossover at v
                idxcrossover2 = ((x-np.sqrt(2*speed*np.log(popsize*np.sqrt(speed))))**2).argmin() # crossover at x_c
                u = np.zeros(space)
                if idxcrossover1 < idxcrossover2:
                    u[idxcrossover2:]              = x[idxcrossover2:]/2
                    u[idxcrossover1:idxcrossover2] = u[idxcrossover2]*np.exp(x[idxcrossover1:idxcrossover2]**2/(2*speed))
                    u[:idxcrossover1]              = u[idxcrossover1]*np.exp(x[:idxcrossover1])
                else:
                    u[idxcrossover2:] = x[idxcrossover2:]/2
                    u[:idxcrossover2] = u[idxcrossover2]*np.exp(x[:idxcrossover2])
        except:
            # if anything in the more elaborate approximations fails, fall back to most basic approximation
            u       = x/2
            u[x<=0] = np.exp(x[x<=0])*u[space0+1]
    elif args.growthterm == "step":
	u      = np.exp(x)
	u[x>0] = 1

if args.enforceboundaries:
  maxvalue = (space-space0)*dx


if args.growthterm == "selection":
    growth = x
    dgrowth = np.ones(space)
elif args.growthterm == "step":
    growth = np.zeros(space)
    growth[x>=0] = 1
    dgrowth = np.zeros(space)
    dgrowth[space0] = 1/dx

if mutationmodel == "diff":
    coeff_prev =  1./(dx*dx) + 0.5*speed/dx
    coeff_next =  1./(dx*dx) - 0.5*speed/dx
    coeff_0    = -2./(dx*dx) + growth
elif mutationmodel == "exp":
    coeff_prev = speed/(dx*dx) + 0.5*(speed - mutationrate + x)/dx
    coeff_next = speed/(dx*dx) - 0.5*(speed - mutationrate + x)/dx
    coeff_0    = -2*speed/(dx*dx) + growth - dgrowth

# coupled Newton-Raphson iterations for each lattice point
for i in range(args.maxsteps):
    # shift profile for terms with derivatives
    u_prev = np.concatenate((np.array([u[0]*u[0]/u[1]]) if u[1]>0 else np.zeros(1),u[:-1])) # exponential decay before lattice
    u_next = np.concatenate((u[1:],np.array([2*u[-1]-u[-2]])))                              # linear increase after lattice
    
    f  = coeff_prev * u_prev + coeff_0 * u + coeff_next * u_next - 2.*u*u
    fu = coeff_0 - 4*u
    
    # more nonlinear terms with exponential kernel
    if mutationmodel == "exp":
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

if args.outfile != None:    fp.close()

