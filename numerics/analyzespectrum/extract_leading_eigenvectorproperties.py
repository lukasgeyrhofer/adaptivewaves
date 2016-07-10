#!/usr/bin/env python

#   adaptive waves with exponential mutation kernel, 
#   get the first two eigenvectors from already sorted file
#   and compute their integrating factor, the "schroedinger" potential and the remaining eigenfunction


import numpy as np
import argparse
import sys,math

def d(y,dx):
    return 0.5/dx*np.diff(np.concatenate([y,np.array([2*y[-1] - y[-2]])]) + np.concatenate([np.array([2*y[0]-y[1]]),y]))


parser = argparse.ArgumentParser()
parser.add_argument("-V","--eigenvectorfile")
parser.add_argument("-E","--eigenvaluefile")
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=1)
parser.add_argument("-m","--mutationrate",type=float,default=1e-2)
parser.add_argument("-T","--numericnoisethreshold",type=float,default=None)
parser.add_argument("-P","--reducePhase",action="store_true",default=False)
args = parser.parse_args()


try:
    eigvecdata = np.genfromtxt(args.eigenvectorfile)
    eigvaldata = np.genfromtxt(args.eigenvaluefile)
    udata = np.genfromtxt(args.ufile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)
    
    
x = udata[:,0]
u = udata[:,1]
s = x-2*u

v = args.speed
m = args.mutationrate

dx     = x[1] - x[0]
space  = len(x)
space0 = (x*x).argmin()

ds = d(s,dx)

eigvec0 = eigvecdata[:,1] + 1j * eigvecdata[:,2]
eigvec1 = eigvecdata[:,3] + 1j * eigvecdata[:,4]
eigval  = np.array([eigvaldata[0,0] + 1j * eigvaldata[0,1],eigvaldata[1,0] + 1j * eigvaldata[1,1]])

if not (args.numericnoisethreshold is None):
    eigvec0[np.absolute(eigvec0) < args.numericnoisethreshold] = 0. + 1j * 0.
    eigvec1[np.absolute(eigvec1) < args.numericnoisethreshold] = 0. + 1j * 0.

if0expint = (s + v - m - eigval[0])/(2*v)
if0       = np.exp(-np.cumsum(if0expint*dx))
if0      /= if0[space0]

if1expint = (s + v - m - eigval[1])/(2*v)
if1       = np.exp(-np.cumsum(if1expint*dx))
if1      /= if1[space0]

phi0 = eigvec0/if0
phi1 = eigvec1/if1

potential0 = (s - v - m - eigval[0])**2/(4*v**2) - 0.5*ds/v - m/v
potential1 = (s - v - m - eigval[1])**2/(4*v**2) - 0.5*ds/v - m/v

if args.reducePhase:
    eigvec0 *= np.exp(-1j*np.imag(eigval[0])/v*x)
    eigvec1 *= np.exp(-1j*np.imag(eigval[1])/v*x)
    phi0    *= np.exp(-0.5*1j*np.imag(eigval[0])/v*x)
    phi0    *= np.exp(-0.5*1j*np.imag(eigval[1])/v*x)
    
    eigvec0 /= eigvec0[np.absolute(eigvec0).argmax()]
    eigvec1 /= eigvec1[np.absolute(eigvec1).argmax()]
    phi0    /= phi0[np.absolute(phi0).argmax()]
    phi1    /= phi1[np.absolute(phi1).argmax()]


for i in range(space):
    print "{:6.2f} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e} {:13.6e}"\
        .format(x[i],u[i],\
            np.real(eigvec0)[i],    np.imag(eigvec0)[i],\
            np.real(phi0)[i],       np.imag(phi0)[i],\
            np.real(potential0)[i], np.imag(potential0)[i],\
            np.real(eigvec1)[i],    np.imag(eigvec1)[i],\
            np.real(phi1)[i],       np.imag(phi1)[i],\
            np.real(potential1)[i], np.imag(potential1)[i])
