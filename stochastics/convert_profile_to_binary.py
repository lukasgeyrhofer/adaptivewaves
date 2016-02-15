#!/usr/bin/python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser(description="Converting numerics profiles for use with stochastic simulations (which uses units and not reduced observables).")
parser_io = parser.add_argument_group(description="Input and Output")
parser_io.add_argument("-i","--infile",required=True,help="Input file obtained from numerics, solving reduced equations in TXT format")
parser_io.add_argument("-o","--outfile",required=True,help="Binary output file, using non-reduced units")
parser_model = parser.add_argument_group(description="General profile options")
mode = parser_model.add_mutually_exclusive_group()
mode.add_argument("-C","--densityprofile",action="store_true",default=False,help="Profiles are rescaled with SCALE^-2 (default)")
mode.add_argument("-U","--fixationprofile",action="store_true",default=False,help="Profiles are rescaled with SCALE")
parser_model.add_argument("-v","--speed",type=float,default=1,help="Reduced adaptation speed (default v=1)")
parser_model.add_argument("-M","--mutationmodel",choices=("diff","exp"),required=True)
parser_exp = parser.add_argument_group(description="Parameters for exponential mutation kernel")
parser_exp.add_argument("-m","--mutationrate",type=float,default=None,help="Reduced mutation rate")
parser_exp.add_argument("-s","--mutationeffect",type=float,default=None,help="Scale in mutation kernel")
parser_diff = parser.add_argument_group(description="Parameters for diffusion mutation kernel")
parser_diff.add_argument("-D","--diffusionconstant",type=float,default=None,help="Diffusion constant sets scale via D^1/3")
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"Could not open input-file '%s'"%args.infile
    exit(1)
    
x = data[:,0]
f = data[:,1]

space  = len(x)
space0 = (x*x).argmin()
dx     = x[1] - x[0]

icount = 0
if args.mutationmodel == "exp":
    scale = args.mutationeffect
    if not args.fixationprofile:
        h04 = np.array([args.mutationrate*scale,args.speed*scale**2,scale,0.],dtype = np.float64)
        dcount = 4
        rescaleprofile = 1./scale**2
    else:
        h04 = np.array([args.mutationrate*scale,args.speed*scale**2,scale],dtype = np.float64)
        dcount = 3 
        rescaleprofile = scale
    if args.diffusionconstant != None:
        print >> sys.stderr,"cannot use diffusionconstant in exponential kernel"
        exit(1)
        
elif args.mutationmodel == "diff":
    scale = np.power(args.diffusionconstant,1/3.)
    if not args.fixationprofile:
        h04 = np.array([args.diffusionconstant,args.speed*scale**2,0.],dtype=np.float64)
        dcount = 3
        rescaleprofile = 1./scale**2
    else:
        h04 = np.array([args.diffusionconstant,args.speed*scale**2],dtype=np.float64)
        dcount = 2
        rescaleprofile = scale
    if args.mutationrate != None or args.mutationeffect != None:
        print >> sys.stderr,"cannot use mutationeffect and/or mutationrate in diffusion kernel"
        exit(1)
else:
    print >> sys.stderr,"mutation kernel has to be specified"
    exit(1)


h01 = np.array([icount,dcount],dtype=np.int32)
h02 = np.array([dx*scale],dtype=np.float64)
h03 = np.array([space,space0],dtype=np.int32)

prof = np.array(f*rescaleprofile,dtype=np.float64)


fp = open(args.outfile,"wb")

h01.tofile(fp)
h02.tofile(fp)
h03.tofile(fp)
h04.tofile(fp)
prof.tofile(fp)

fp.close()





    