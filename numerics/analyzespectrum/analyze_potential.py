#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

def d(y,dx = 1e-2):
    return 0.5 * np.diff(np.concatinate([np.array([2*y[0] - y[1]]),y]) + np.concatinate([y,np.array([2*y[-1] - y[-2]])]))/dx


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)


x = data[:,0]
space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

pot0r = data[:,6]
pot0i = data[:,7]
pot1r = data[:,12]
pot1i = data[:,13]

pot0 = pot0r + 1j * pot0i
pot1 = pot1r + 1j * pot1i


a = (pot0r[2:] - pot0r[1:space-1]) * (pot0r[1:space-1] - pot0r[:space-2])
b = (pot1r[2:] - pot1r[1:space-1]) * (pot1r[1:space-1] - pot1r[:space-2])

for i in np.where(a<0)[0]:
    print "{} {}".format(x[i+1],pot0r[i+1]),
print

for i in np.where(b<0)[0]:
    print "{} {}".format(x[i+1],pot1r[i+1]),
print
