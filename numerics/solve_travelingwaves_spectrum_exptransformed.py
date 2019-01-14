#!/usr/bin/env python
# -*- coding: utf-8 -*-


# solve spectrum of tuned adaptation operator with an exponential mutation kernel

import numpy as np
from scipy import linalg
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=0)
parser.add_argument("-m","--mutationrate",type=float,default=0)
parser.add_argument("-n","--closurelevel",type=int,default=1)
parser.add_argument("-R","--reduction",type=int,default=1)
parser.add_argument("-V","--compute_eigenvectors",default=False,action="store_true")
args = parser.parse_args()


try:
    udata = np.genfromtxt(args.ufile)
except:
    raise IOError("could not open file '{}'".format(args.ufile))
    exit(1)

x = udata[::args.reduction,0]
un = 4. * args.closurelevel/(args.closurelevel + 1. ) * udata[::args.reduction,1]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

# state generalized eigenvalue equation as 
#   A c_l = l B c_l
# with eigenvalue l to eigenvector c_l
# overall, we want to solve the equation
#   v c'' + (S-m-v)c' + (S+S')c = l c + l c'

s = x - un
ds = 0.5/dx * (np.diff(np.concatenate((s,np.array([2*s[-1] - s[-2]])))) + np.diff(np.concatenate((np.array([2*s[0]-s[1]]),s))))

dxm = 0.5/dx* ( np.diag(np.ones(space-1),k=1) - np.diag(np.ones(space-1),k=-1))

A  = args.speed/(dx*dx)*(np.diag(np.ones(space-1),k=1) - 2.*np.diag(np.ones(space)) + np.diag(np.ones(space-1),k=-1)) # speed times second derivative
A += np.dot(np.diag(s - args.mutationrate + args.speed),dxm)                                                            # 
A += np.diag(s+ds)

B = np.eye(space) + dxm


if args.compute_eigenvectors:
    w,v = linalg.eig(a = A,b = B,overwrite_a = True,overwrite_b = True,left = False,right = True )
else:
    w   = linalg.eig(a = A,b = B,overwrite_a = True,overwrite_b = True,left = False,right = False)

for i in range(len(w)):
    print i,np.real(w[i]), np.imag(w[i])

if args.compute_eigenvectors:
    for i in range(len(v)):
        print >> sys.stderr,x[i],
        for j in range(len(v)):
            print >> sys.stderr,np.real(v[i,j]),np.imag(v[i,j]),
        print >> sys.stderr
