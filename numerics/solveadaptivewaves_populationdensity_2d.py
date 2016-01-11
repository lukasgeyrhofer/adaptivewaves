#!/usr/bin/env python

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()

parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser_alg_c = parser_alg.add_mutually_exclusive_group()
parser_alg_c.add_argument("-i","--c2file",default=None,description="Start with previous 2d profile")
parser_alg_c.add_argument("-c","--c1file",default=None,description="Start by expanding 1d profile")
parser_alg.add_argument("-o","--outfile",default=None)
parser_alg.add_argument("-u","--ufile",default=None,description="File with profile for fixation probability (n=1)")
parser_alg.add_argument("-a","--alpha",type=float,default=1)
parser_alg.add_argument("-S","--maxsteps",type=int,default=1000)
parser_alg.add_argument("-O","--outputstep",type=int,default=0,help="print expression <uu|c> at each OUTPUTSTEP steps to show convergence (default: 0 [=OFF])")

parser_params = parser.add_argument_group(description="####   Profile parameters   ####")
parser_params.add_argument("-v","--speed",default=1,type=float)

args = parser.parse_args()


try:
    udata = np.genfromtxt(args.ufile)
    x = udata[:,0]
    u = udata[:,1] * 2/3. # numerical profiles are usually for n=1, have to rescale to comply with closure at n=2
except:
    print >> sys.stderr,"could not open ufile"
    exit(1)

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]
speed = args.speed

try:
    c2data = np.genfromtxt(args.c2file)
    c2 = np.reshape(c2data[:,2],(space,space))
    print >> sys.stderr,"# starting from c2 profile (file: '%s')"%args.c2file
except:
    try:
        c1data = np.genfromtxt(args.c1file)
        c1 = c1data[:,1]
        c2 = np.outer(c1,c1)
        print >> sys.stderr,"# starting from c1 profile (file: '%s')"%args.c1file
    except:
        c1 = np.exp(-0.5*x*x/speed)
        c2 = np.outer(c1,c1)
        print >> sys.stderr,"# starting from scratch"
        
# rescale profile to comply with constraint initially
c2 /= np.dot(np.dot(c2,u),u)*dx*dx


# generate coefficient matrices
coeff_xprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_xnext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_yprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_ynext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_00    = -4*np.ones((space,space))/(dx*dx) + np.outer(x-3*u,np.ones(space)) + np.outer(np.ones(space),x-3*u)

# Newton-Raphson iterations
for i in range(args.maxsteps):
    # generate c2 profiles shifted by 1 lattice point in each direction
    c2_xprev = np.concatenate(( np.reshape(2*c2[:,0]-c2[:,1],(space,1)),  c2[:,:-1] ),axis=1)   # assume profile extends linearly,
    c2_xnext = np.concatenate(( c2[:,1:], np.reshape(2*c2[:,-1]-c2[:,-2],(space,1)) ),axis=1)   # which is much faster than an
    c2_yprev = np.concatenate(( np.reshape(2*c2[0,:]-c2[1,:],(1,space)),  c2[:-1,:] ),axis=0)   # exponential decay, where one has
    c2_ynext = np.concatenate(( c2[1:,:], np.reshape(2*c2[-1,:]-c2[-2,:],(1,space)) ),axis=0)   # to check for each element == 0
                                                                                                # to be allowed to divide elements
    
    f  = coeff_xnext * c2_xnext + coeff_xprev * c2_xprev
    f += coeff_ynext * c2_ynext + coeff_yprev * c2_yprev
    f += coeff_00 * c2 + 2*np.diag(np.dot(c2,u))
    
    fc = coeff_00 + 2*np.diag(u)
    
    # NR step
    c2 -= args.alpha * f/fc
    
    # rescaling
    c2 /= np.dot(np.dot(c2,u),u)*dx*dx
    
    # output?
    if args.outputstep > 0:
        if i % args.outputstep == 0:
            uuc = np.dot(np.dot(c2,u),u*u)*dx*dx
            print >> sys.stderr,i,uuc

# ... and final output
if args.outfile == None:    fp = sys.stdout
else:                       fp = open(args.outfile,"w")
for i in range(space):
    for j in range(space):
        print >>fp,"%lf %lf %.14e",x[i],x[j],c2[i,j]
    print >>fp
fp.close()
