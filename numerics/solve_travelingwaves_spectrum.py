#!/usr/bin/env python
# -*- coding: utf-8 -*-


# =========================================================================== #
#   compute eigenvalues and eigenvectors in tuned traveling wave models       #
#   using reduced units.                                                      #
# =========================================================================== #
#                                                                             #
# simplest usage:                                                             #
#   ./solve_travelingwaves_spectrum.py -u fixationprofile.txt -v 2            #
# to compute eigenvalues for a fixation probability profile given in the      #
# file 'fixationprofile.txt' and for a speed of 2                             #
#                                                                             #
# see option --help for available parameters                                  #
#                                                                             #
# Code can be used as is. If reference is needed then cite                    #
#   Hallatschek, Geyrhofer, Genetics 202 (3), pp1201-1227, 2016               #
#   DOI: 10.1534/genetics.115.181271                                          #
#                                                                             #
#                                                                             #
# Lukas Geyrhofer, 2015 - 2016                                                #
#                                                                             #
# *************************************************************************** #


import numpy as np
from scipy import linalg
import argparse
import sys,math

parser = argparse.ArgumentParser(description = "Compute eigenvalues (and eigenvectors) in tuned traveling wave models")
parser_alg = parser.add_argument_group(description="####   Algorithm and IO parameters   ####")
parser_alg.add_argument("-u","--ufile")
parser_alg.add_argument("-n","--closurelevel",type=int,default=1)
parser_alg.add_argument("-R","--reduction",type=int,default=1)
parser_alg.add_argument("-V","--compute_eigenvectors",default=False,action="store_true")

parser_params = parser.add_argument_group(description="####   Profile parameters   ####")
parser_params.add_argument("-v","--speed",type=float,default=1)
parser_params.add_argument("-m","--mutationrate",type=float,default=None)
parser_params.add_argument("-M","--mutationmodel",choices=("diff","exp"),default=None)
parser_params.add_argument("-G","--growthterm",choices=("selection","step"),default="selection")
args = parser.parse_args()

try:
    udata = np.genfromtxt(args.ufile)
except:
    # an error occurred during loading. aborting
    print >> sys.stderr,"could not open file for fixation probability profile '%s'"%args.ufile
    parser.print_help()
    exit(1)


# generate necessary variables and profiles
x = udata[::args.reduction,0]
w = 2*2*args.closurelevel/(args.closurelevel + 1 ) * udata[::args.reduction,1] # assume profile was computed for n=1, thus already rescale with additional factor 2
space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]

# set speed
speed  = args.speed
# set mutationrate and mutationmodel
if args.mutationmodel == None:
    if args.mutationrate == None:
        mutationmodel = "diff"
    else:
        mutationmodel = "exp"
        mutationrate = args.mutationrate
else:
    if args.mutationrate != None and args.mutationmodel == "diff":
        print >> sys.stderr,"Cannot set mutation rate (option -m) in diffusion mutation kernel"
        parser.print_help()
        exit(1)
    mutationmodel = args.mutationmodel
    if mutationmodel == "exp":
        if args.mutationrate == None:
            mutationrate = 0.1 # set default value
        else:
            mutationrate = args.mutationrate

# set growth term
if args.growthterm == "selection": # adaptation waves
    growth = x
elif args.growthterm == "step":    # invasion waves
    growth = np.zeros(space)
    growth[x>0] = 1.

# effective selection
s = growth - w

# generate matrices for first and second derivative
dx1 = 0.5/dx*(np.diag(np.ones(space-1),k=1) - np.diag(np.ones(space-1),k=-1))
dx2 = (np.diag(np.ones(space-1),k=-1) -2*np.diag(np.ones(space),k=0) + np.diag(np.ones(space-1),k=1))/(dx*dx)

# construct linear operator for both, diffusion and exponential, mutation models
if mutationmodel == "diff":
    A  = np.diag(growth-2*u) # linear term
    A += speed*dx1           # first derivative
    A += dx2                 # second derivative
    
    # actual computation of eigenvalues
    if args.compute_eigenvectors:
        w,v = linalg.eig(a = A,overwrite_a = True, left = False, right = True )
    else:
        w   = linalg.eig(a = A,overwrite_a = True, left = False, right = False)
    
    
elif mutationmodel == "exp":
    # need derivative of effective selection term
    ds = 0.5/dx*np.diff(np.concatenate([np.array([2*s[0]-s[1]]),s[:-1]]) + np.concatenate([s[1:],np.array([2*s[-1]-s[-2]])]))
    
    A  = np.diag(s+ds)       # linear term
    A += np.dot(np.diag(s+speed-mutationrate),dx1)
                             # first derivative
    A += speed*dx2           # second derative
    
    B = np.eye(space) + dx1  # exponential mutation kernel also need second matrix for generalized eigenvalue problem

    # actual computation of eigenvalues
    if args.compute_eigenvectors:
        w,v = linalg.eig(a = A,b = B,overwrite_a = True,overwrite_b = True,left = False,right = True )
    else:
        w   = linalg.eig(a = A,b = B,overwrite_a = True,overwrite_b = True,left = False,right = False)


# output of eigenvalues
for i in range(len(w)):
    print i,np.real(ev[i]),np.imag(ev[i])

# if eigenvectors are computed, also print them
if args.compute_eigenvectors:
    for i in range(len(w)):
        print >> sys.stderr,x[i],
        for j in range(len(v)):
            print >> sys.stderr,np.real(v[i,j]),np.imag(v[i,j]),
        print >> sys.stderr
    
