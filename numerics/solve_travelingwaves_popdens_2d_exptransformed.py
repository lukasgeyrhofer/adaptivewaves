#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()

parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")

parser_alg_c = parser_alg.add_mutually_exclusive_group()
parser_alg_c.add_argument("-i","--c2file",default=None,help="Start with previous 2d profile")
parser_alg_c.add_argument("-c","--c1file",default=None,help="Start by expanding 1d profile")

parser_alg.add_argument("-o","--outfile",default=None,help="Write output to OUTFILE instead of stdout")
parser_alg.add_argument("-u","--ufile",default=None,help="File with profile for fixation probability (n=1)")
parser_alg.add_argument("-a","--alpha",type=float,default=1,help="'Speed' of Newton-Raphson iteration: c2 -= alpha f/f'. Classical NR: alpha=1. Slower convergence but more stability for alpha<1 [default: 1]")
parser_alg.add_argument("-S","--maxsteps",type=int,default=1000,help="Number of iteration steps [default: 1000]")
parser_alg.add_argument("-O","--outputstep",type=int,default=0,help="print expression <uu|c> at each OUTPUTSTEP steps to show convergence (default: 0 [=OFF])")

parser_params = parser.add_argument_group(description="####   Profile parameters   ####")
parser_params.add_argument("-v","--speed",default=1,type=float,help="Adaptation speed [default: 1]")
parser_params.add_argument("-m","--mutationrate",default=1e-2,type=float,help="Mutation rate used in exponential mutation kernel [default: 0.01]")

args = parser.parse_args()


try:
    udata = np.genfromtxt(args.ufile)
    x = udata[:,0]
    w = udata[:,1] * 2. # numerical profiles are usually for n=1, have to rescale to comply with closure at n=2
    u = w/3.            
except:
    print >> sys.stderr,"could not open ufile"
    exit(1)

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]
speed  = args.speed
mutationrate = args.mutationrate

# make array larger such that boundary conditions can be directly implemented in it
c2 = np.zeros((space+2,space+2))

# starting conditions
try:
    c2data = np.genfromtxt(args.c2file)
    c2[1:space+1,1:space+1] = np.reshape(c2data[:,2],(space,space))
    print >> sys.stderr,"# starting from c2 profile (file: '%s')"%args.c2file
except:
    try:
        c1data = np.genfromtxt(args.c1file)
        c1 = c1data[:,1]
        c2[1:space+1,1:space+1] = np.outer(c1,c1)
        print >> sys.stderr,"# starting from c1 profile (file: '%s')"%args.c1file
    except:
        c1 = np.exp(-0.5*x*x/speed)
        c2[1:space+1,1:space+1] = np.outer(c1,c1)
        print >> sys.stderr,"# starting from scratch"
        
# rescale profile to comply with constraint initially
c2 /= np.dot(np.dot(c2[1:space+1,1:space+1],u),u)*dx*dx


# generate coefficient matrices
s  = x - 4*u
ds = 0.5*np.diff( np.concatenate([np.array([2*s[0] - s[1]]),s]) + np.concatenate([s,np.array([2*s[-1] - s[-2]])]) )/dx

mat_ones = np.ones((space,space))
lin_ones = np.ones(space)

idx1 = 1./(dx)
idx2 = 1./(dx*dx)
idx3 = 1./(dx*dx*dx)

const_sdx2 = idx2/6.
const_tdx  = idx1/3.
const_tt   = 2./3.

coeff_xp_yp = ( speed*idx3 + 0.5*(speed-mutationrate)*idx2             ) * mat_ones + 0.25*(np.outer(lin_ones,s)    + np.outer(s,lin_ones))*idx2
coeff_xp_y0 = (-speed*idx3 + speed*idx2 + 0.5*(speed-mutationrate)*idx1) * mat_ones + 0.50*(np.outer(lin_ones,s+ds) + np.outer(s,lin_ones))*idx1
coeff_xp_ym = (            - 0.5*(speed-mutationrate)*idx2             ) * mat_ones - 0.25*(np.outer(lin_ones,s)    + np.outer(s,lin_ones))*idx2

coeff_x0_yp = (-speed*idx3 + speed*idx2 + 0.5*(speed-mutationrate)*idx1) * mat_ones + 0.50*(np.outer(lin_ones,s)    + np.outer(s+ds,lin_ones))*idx1
coeff_x0_y0 = (-4.*speed*idx2                                          ) * mat_ones +       np.outer(lin_ones,s+ds) + np.outer(s+ds,lin_ones)
coeff_x0_ym = ( speed*idx3 + speed*idx2 - 0.5*(speed-mutationrate)*idx1) * mat_ones - 0.50*(np.outer(lin_ones,s)    + np.outer(s+ds,lin_ones))*idx1

coeff_xm_yp = (            - 0.5*(speed-mutationrate)*idx2             ) * mat_ones - 0.25*(np.outer(lin_ones,s)    + np.outer(s,lin_ones))*idx2
coeff_xm_y0 = ( speed*idx3 + speed*idx2 - 0.5*(speed-mutationrate)*idx1) * mat_ones - 0.50*(np.outer(lin_ones,s+ds) + np.outer(s,lin_ones))*idx1
coeff_xm_ym = (-speed*idx3 + 0.5*(speed-mutationrate)*idx2             ) * mat_ones + 0.25*(np.outer(lin_ones,s)    + np.outer(s,lin_ones))*idx2

# completely symmetrized version
#fc = coeff_x0_y0 + np.diag(2.*w/3.) + np.diag( (w[:space-1] - w[1:])/(6.*dx) ,k=1) + np.diag( (w[:space-1] - w[1:])/(6.*dx) ,k=-1)
# ignore symmetry (??)
fc = coeff_x0_y0 + np.diag(2.*w/3.) + np.diag( const_tdx * w[:space-1],k=1) - np.diag( const_tdx * w[:space-1], k = -1)

# Newton-Raphson iterations
for i in range(args.maxsteps):
    
    f  = coeff_xp_yp * c2[2:space+2,2:space+2] + coeff_xp_y0 * c2[2:space+2,1:space+1] + coeff_xp_ym * c2[2:space+2,0:space]
    f += coeff_x0_yp * c2[1:space+1,2:space+2] + coeff_x0_y0 * c2[1:space+1,1:space+1] + coeff_x0_ym * c2[1:space+1,0:space]
    f += coeff_xm_yp * c2[0:space  ,2:space+2] + coeff_xm_y0 * c2[0:space  ,1:space+1] + coeff_xm_ym * c2[0:space  ,0:space]
    
    # compute c1 for terms on diagonal, assume linear extrapolation
    b  = np.dot(c2[1:space+1,1:space+1],w)
    bp = np.concatenate([b[:space-1],np.array([2*b[space-1] - b[space-2]])])
    bm = np.concatenate([np.array([2*b[0] - b[1]]),b[1:]])
    
    # apply diagonal terms arising from delta
    f += np.diag(b*const_tt + (bm + bp)*const_sdx2 ) 
    f += np.diag(const_tdx * (b[1:]-b[:space-1]),k= 1) + np.diag(-const_sdx2*b[1:space-1],k= 2)
    f += np.diag(const_tdx * (b[1:]-b[:space-1]),k=-1) + np.diag(-const_sdx2*b[1:space-1],k=-2)
    
    # NR step
    c2[1:space+1,1:space+1] -= args.alpha * f/fc
    
    # profile should never be smaller than zero
    c2[c2<0] = 0
    
    # rescaling to enforce constraint
    c2 /= np.dot(np.dot(c2[1:space+1,1:space+1],u),u)*dx*dx
    
    # output of <uu|c>?
    if args.outputstep > 0:
        if i % args.outputstep == 0:
            uuc = np.dot(np.dot(c2[1:space+1,1:space+1],u),u*u)*dx*dx
            print >> sys.stderr,i,uuc,1.-constraint

# ... and final output
if args.outfile == None:    fp = sys.stdout
else:                       fp = open(args.outfile,"w")
for i in range(space):
    for j in range(space):
        print >>fp,"{:.2f} {:.2f} {:.10e}".format(x[i],x[j],c2[i+1,j+1])
    print >>fp
fp.close()










