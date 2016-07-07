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

parser_alg.add_argument("-o","--outfile",default=None)
parser_alg.add_argument("-u","--ufile",default=None,help="File with profile for fixation probability (n=1)")
parser_alg.add_argument("-a","--alpha",type=float,default=1)
parser_alg.add_argument("-S","--maxsteps",type=int,default=1000)
parser_alg.add_argument("-O","--outputstep",type=int,default=0,help="print expression <uu|c> at each OUTPUTSTEP steps to show convergence (default: 0 [=OFF])")

parser_params = parser.add_argument_group(description="####   Profile parameters   ####")
parser_params.add_argument("-v","--speed",default=1,type=float)
parser_params.add_argument("-m","--mutationrate",default=1e-2,type=float)

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
speed = args.speed
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
s = x-4*u
ds = 0.5*np.diff( np.concatenate([np.array([2*s[0] - s[1]]),s]) + np.concatenate([s,np.array([2*s[-1] - s[-2]])]) )/dx
ones = np.ones((space,space))
lin = np.ones(space)

coeff_xp_yp = (speed/dx**3 + 0.5*(speed-mutationrate)/dx**2)             * ones + 0.25*(np.outer(lin,s)    + np.outer(s,lin))/dx**2
coeff_xp_y0 = (-speed/dx**3 + speed/dx**2 + 0.5*(speed-mutationrate)/dx) * ones + 0.50*(np.outer(lin,s+ds) + np.outer(s,lin))/dx
coeff_xp_ym = (-0.50*(speed-mutationrate)/dx**2)                         * ones - 0.25*(np.outer(lin,s)    + np.outer(s,lin))/dx**2

coeff_x0_yp = (-speed/dx**3 + speed/dx**2 + 0.5*(speed-mutationrate)/dx) * ones + 0.50*(np.outer(lin,s)    + np.outer(s+ds,lin))/dx
coeff_x0_y0 = (-4.*speed/dx**2)                                          * ones +       np.outer(lin,s+ds) + np.outer(s+ds,lin)
coeff_x0_ym = (speed/dx**3 + speed/dx**2 - 0.5*(speed-mutationrate)/dx)  * ones - 0.50*(np.outer(lin,s)    + np.outer(s+ds,lin))/dx

coeff_xm_yp = (-0.50*(speed-mutationrate)/dx**2)                         * ones - 0.25*(np.outer(lin,s)    + np.outer(s,lin))/dx**2
coeff_xm_y0 = (speed/dx**3 + speed/dx**2 - 0.5*(speed-mutationrate)/dx)  * ones - 0.50*(np.outer(lin,s+ds) + np.outer(s,lin))/dx
coeff_xm_ym = (-speed/dx**3 + 0.5*(speed-mutationrate)/dx**2)            * ones + 0.25*(np.outer(lin,s)    + np.outer(s,lin))/dx**2

fc = coeff_x0_y0 + 2*np.diag(w)/3. + np.diag(w[:space-1] - w[1:],k=1)/(6.*dx) + np.diag(w[1:] - w[:space-1],k=-1)/(6.*dx)

const_sdx2 = 1./(6.*dx**2)
const_tdx  = 1./(3.*dx)
const_ft   = 4./3.

# Newton-Raphson iterations
for i in range(args.maxsteps):
    
    f  = coeff_xp_yp * c2[2:space+2,2:space+2] + coeff_xp_y0 * c2[2:space+2,1:space+1] + coeff_xp_ym * c2[2:space+2,0:space]
    f += coeff_x0_yp * c2[1:space+1,2:space+2] + coeff_x0_y0 * c2[1:space+1,1:space+1] + coeff_x0_ym * c2[1:space+1,0:space]
    f += coeff_xm_yp * c2[0:space  ,2:space+2] + coeff_xm_y0 * c2[0:space  ,1:space+1] + coeff_xm_ym * c2[0:space  ,0:space]
    
    b = np.dot(c2[1:space+1,1:space+1],w)
    bp = np.concatenate([b[:space-1],np.zeros(1)])
    bm = np.concatenate([np.zeros(1),b[1:]])
    
    f += np.diag(bp*const_sdx2 + b*const_ft + bm*const_sdx2) 
    f += np.diag(const_tdx * (b[1:]-b[:space-1]),k= 1) + np.diag(-const_sdx2*b[1:space-1],k= 2)
    f += np.diag(const_tdx * (b[1:]-b[:space-1]),k=-1) + np.diag(-const_sdx2*b[1:space-1],k=-2)
    
    # NR step
    c2[1:space+1,1:space+1] -= args.alpha * f/fc
    
    c2[c2<0] = 0
    
    # rescaling
    c2 /= np.dot(np.dot(c2[1:space+1,1:space+1],u),u)*dx*dx
    
    # output?
    if args.outputstep > 0:
        if i % args.outputstep == 0:
            uuc = np.dot(np.dot(c2[1:space+1,1:space+1],u),u*u)*dx*dx
            print >> sys.stderr,i,uuc

# ... and final output
if args.outfile == None:    fp = sys.stdout
else:                       fp = open(args.outfile,"w")
for i in range(space):
    for j in range(space):
        print >>fp,"{:.2f} {:.2f} {:.10e}".format(x[i],x[j],c2[i+1,j+1])
    print >>fp
fp.close()










