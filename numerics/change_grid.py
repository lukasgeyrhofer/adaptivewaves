#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-s","--space",type=int,default=None)
parser.add_argument("-z","--space0",type=int,default=None)
parser.add_argument("-d","--dx",type=float,default=None)
parser.add_argument("-e","--extractionlength",type=int,default=100)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

x = data[:,0]
u = data[:,1]

if args.space == None:  space = len(x)
else:                   space = args.space

if args.space0 == None: space0 = (x*x).argmin()
else:                   space0 = args.space0

if args.dx == None:     dx = x[1]-x[0]
else:                   dx = args.dx

nx = (np.arange(space)-space0)*dx
nu = np.interp(nx,x,u)

if nx[0] < x[0]:
    idx = ((nx-x[0])**2).argmin()
    e   = args.extractionlength
    sx  = np.sum(x[:e])
    sxx = np.dot(x[:e],x[:e])
    try:
        sy    = np.sum(np.log(u[:e]))
        sxy   = np.dot(np.log(u[:e]),x[:e])
        decay = (e*sxy-sy*sx)/(e*sxx-sx*sx)
        nu[:idx] = np.exp(decay*(nx[:idx]-nx[idx]))*nu[idx]
    except:
        nu[:idx] = 0

if nx[-1] > x[-1]:
    idx   = ((nx-x[-1])**2).argmin()
    e     = args.extractionlength
    sx    = np.sum(x[-e:])
    sxx   = np.dot(x[-e:],x[-e:])
    sy    = np.sum(u[-e:])
    sxy   = np.dot(u[-e:],x[-e:])
    increase = (e*sxy-sy*sx)/(e*sxx-sx*sx)
    nu[idx:] = nu[idx] + (nx[idx:]-nx[idx])*increase


if args.outfile == None:    fp = sys.stderr
else:                       fp = open(args.outfile,"w")
for i in range(space):
    print >> fp, "%.6lf %.14e"%(nx[i],nu[i])
if args.outfile != None:    fp.close()

