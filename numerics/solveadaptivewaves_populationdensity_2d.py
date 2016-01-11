#!/usr/bin/env pyton

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",default=1,type=float)
parser.add_argument("-S","--maxsteps",type=int,default=100)
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

coeff_0 = np.outer(x-2*u,o) + np.outer(o,x-2*u)

print coeff_0




