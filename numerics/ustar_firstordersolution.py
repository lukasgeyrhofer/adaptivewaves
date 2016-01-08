#!/usr/bin/env python


from scipy.special import erfi
import numpy as np
import math
import argparse


def ustar(x,ps,fv):
  return np.exp(x**2/(2.*fv))/(ps+np.sqrt(2*math.pi/fv)*erfi(x/np.sqrt(2.*fv)))



parser = argparse.ArgumentParser()
parser.add_argument("-N","--populationsize",type=float,default=1e7)
parser.add_argument("-V","--fitnessvariance",type=float,default=1e-5)
parser.add_argument("-s","--space",type=int,default=100)
parser.add_argument("-z","--zero",type=int,default=50)
parser.add_argument("-d","--latticespacing",type=float,default=1e-3)

args = parser.parse_args()



x = np.array([(i-args.zero) * args.latticespacing for i in range(args.space)])
u = ustar(x,args.populationsize,args.fitnessvariance)

for i in range(len(x)):
  print "%9.6lf %21.14e"%(x[i],u[i])




