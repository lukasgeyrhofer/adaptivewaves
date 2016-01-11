#!/usr/bin/env python

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",default=1,type=float)
parser.add_argument("-S","--maxsteps",type=int,default=10)
parser.add_argument("-a","--alpha",type=float,default=1)
args = parser.parse_args()


try:
    udata = np.genfromtxt(args.ufile)
    x = udata[:,0]
    u = udata[:,1]
except:
    print >> sys.stderr,"could not open ufile"
    exit(1)

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]

speed = args.speed

c1 = np.exp(-0.5*x*x/speed)

c1 /= np.dot(c1,u)*dx

c2 = np.outer(c1,c1)

o = np.ones(space)

coeff_xprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_xnext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_yprev = (1/(dx*dx) - 0.5*speed/dx)*np.ones((space,space))
coeff_ynext = (1/(dx*dx) + 0.5*speed/dx)*np.ones((space,space))
coeff_00    = -4*np.ones((space,space))/(dx*dx) + np.outer(x-2*u,o) + np.outer(o,x-2*u)

uu = np.outer(u,o) * np.outer(o,u) * dx * dx


for i in range(args.maxsteps):
    print i
    c2_xprev = np.concatenate(( np.reshape(2*c2[:,0]-c2[:,1],(space,1)),  c2[:,:-1] ),axis=1)
    c2_xnext = np.concatenate(( c2[:,1:], np.reshape(2*c2[:,-1]-c2[:,-2],(space,1)) ),axis=1)
    c2_yprev = np.concatenate(( np.reshape(2*c2[0,:]-c2[1,:],(1,space)),  c2[:-1,:] ),axis=0)
    c2_ynext = np.concatenate(( c2[1:,:], np.reshape(2*c2[-1,:]-c2[-2,:],(1,space)) ),axis=0)
    
    f  = coeff_xnext * c2_xnext + coeff_xprev * c2_xprev + coeff_ynext * c2_ynext + coeff_yprev * c2_yprev + coeff_00 * c2 + 2*np.diag(np.dot(c2,u))
    fc = coeff_00 + 2*np.diag(u)
    
    c2 -= args.alpha * f/fc
    
    c2 /= np.dot(np.dot(c2,u),u)*dx*dx
    
for i in range(space):
    for j in range(space):
        print x[i],x[j],c2[i,j]
    print