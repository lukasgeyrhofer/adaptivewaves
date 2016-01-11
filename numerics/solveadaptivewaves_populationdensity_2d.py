#!/usr/bin/env python

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()

parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser_alg.add_argument("-u","--ufile",default=None,description="File with profile for fixation probability (n=1)")
parser_alg_c = parser_alg.add_mutually_exclusive_group()
parser_alg_c.add_argument("-i","--c2file",default=None,description="Start with previous 2d profile")
parser_alg_c.add_argument("-c","--c1file",default=None,description="Start by expanding 1d profile")
parser_alg.add_argument("-a","--alpha",type=float,default=1)
parser_alg.add_argument("-S","--maxsteps",type=int,default=10)

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
except:
    try:
        c1data = np.genfromtxt(args.c1file)
        c1 = c1data[:,1]
        c2 = np.outer(c1,c1)
    except:
        c1 = np.exp(-0.5*x*x/speed)
        c2 = np.outer(c1,c1)
        

c2 /= np.dot(np.dot(c2,u),u)*dx*dx


coeff_xprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_xnext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_yprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_ynext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_00    = -4*np.ones((space,space))/(dx*dx) + np.outer(x-3*u,np.ones(space)) + np.outer(np.ones(space),x-3*u)


for i in range(args.maxsteps):
    print >> sys.stderr,i
    c2_xprev = np.concatenate(( np.reshape(2*c2[:,0]-c2[:,1],(space,1)),  c2[:,:-1] ),axis=1)
    c2_xnext = np.concatenate(( c2[:,1:], np.reshape(2*c2[:,-1]-c2[:,-2],(space,1)) ),axis=1)
    c2_yprev = np.concatenate(( np.reshape(2*c2[0,:]-c2[1,:],(1,space)),  c2[:-1,:] ),axis=0)
    c2_ynext = np.concatenate(( c2[1:,:], np.reshape(2*c2[-1,:]-c2[-2,:],(1,space)) ),axis=0)
    
    f  = coeff_xnext * c2_xnext + coeff_xprev * c2_xprev
    f += coeff_ynext * c2_ynext + coeff_yprev * c2_yprev
    f += coeff_00 * c2 + 2*np.diag(np.dot(c2,u))
    
    fc = coeff_00 + 2*np.diag(u)
    
    c2 -= args.alpha * f/fc
    
    c2 /= np.dot(np.dot(c2,u),u)*dx*dx
    
for i in range(space):
    for j in range(space):
        print x[i],x[j],c2[i,j]
    print
