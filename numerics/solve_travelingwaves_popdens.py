#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()

parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser_alg.add_argument("-i","--infile",help="start from profile given by INFILE")
parser_alg.add_argument("-o","--outfile",default=None,help="write output to OUTFILE instead of stdout")
parser_alg.add_argument("-u","--ufile",help="File with profile for fixation probability (n=1)")
parser_alg.add_argument("-S","--maxsteps",type=int,default=10000)
parser_alg.add_argument("-a","--alpha",type=float,default=1.,help="Slower convergence but more stability for alpha<1")
parser_alg.add_argument("-O","--outputstep",type=int,default=0,help="print expression <uu|c> at each OUTPUTSTEP steps to show convergence (default: 0 [=OFF])")

parser_params = parser.add_argument_group(description = "####   Profile parameters    ####")
parser_params.add_argument("-v","--speed",type=float,default=1.,help="Adaptation speed (default: 1)")
parser_params.add_argument("-M","--mutationmodel",choices=("diff","exp"),default="diff",help="Mutation kernel (default: \"diff\")")
parser_params.add_argument("-m","--mutationrate",type=float,default=1e-2,help="Mutation rate (only used for exp kernel, default: 1e-2)")
parser_params.add_argument("-G","--growthterm",choices=("selection","step"),default="selection",help="Growth is either given by linear gradient for adaptation (\"selection\") or as step function for Fisher waves (\"step\")") 

args = parser.parse_args()


try:
    udata = np.genfromtxt(args.ufile)
    x = udata[:,0]
    u = udata[:,1]
    space  = len(x)
    space0 = (x*x).argmin()
    dx     = x[1] - x[0]
except:
    print >> sys.stderr,"could not open ufile"
    exit(1)

try:
    cdata = np.genfromtxt(args.infile)
    c = cdata[:,1]
    speed = args.speed
    mutationrate = args.mutationrate
except:
    print >> sys.stderr,"# starting from scratch"
    speed = args.speed
    mutationrate= args.mutationrate
    c = np.exp(-x*x/(2*speed))
    norm = np.dot(u,c)*dx
    c/=norm

assert speed > mutationrate,"note that (speed ~ variance + mutationrate) ! check your parameters!"

if args.growthterm == "selection":
    growth = x
    dgrowth = np.ones(space)
elif args.growthterm == "step":
    growth = np.zeros(space)
    growth[x>=0] = 1
    dgrowth = np.zeros(space)
    dgrowth[space0] = 1/dx

if args.mutationmodel == "diff":
    coeff_prev = 1./(dx*dx) - 0.5*speed/dx
    coeff_next = 1./(dx*dx) + 0.5*speed/dx
    coeff_0    = -2./(dx*dx) + growth - 2*u
if args.mutationmodel == "exp":
    u_prev = np.concatenate((np.array([u[0]*u[0]/u[1]]) if u[1]>0 else np.zeros(1),u[:-1]))
    u_next = np.concatenate((u[1:],np.array([2*u[-1]-u[-2]])))
    du = 0.5*(u_next - u_prev)/dx
    coeff_prev = speed/(dx*dx) - 0.5*(speed-mutationrate+growth-2*u)/dx
    coeff_next = speed/(dx*dx) + 0.5*(speed-mutationrate+growth-2*u)/dx
    coeff_0    = -2*speed/(dx*dx) + growth - 2*u + dgrowth - 2*du
    
for i in range(args.maxsteps):
    c_prev = np.concatenate((np.array([c[0]*c[0]/c[1]]) if c[1]>0 else np.zeros(1),c[:-1]))
    c_next = np.concatenate((c[1:],np.array([c[-1]*c[-1]/c[-2]]) if c[-2]>0 else np.zeros(1)))
  
    f  = coeff_prev * c_prev + coeff_0 * c + coeff_next * c_next
    fc = coeff_0
    
    c -= args.alpha*f/fc
    
    c /= np.dot(c,u)*dx
    
    if args.outputstep > 0:
	if i%args.outputstep == 0:
	    print >> sys.stderr,i,np.dot(u*u,c)*dx

if args.outfile == None:    fp = sys.stdout
else:                       fp = open(args.outfile,"w")
for i in range(space):
    print >>fp, "%lf %.14e"%(x[i],c[i])
fp.close()
