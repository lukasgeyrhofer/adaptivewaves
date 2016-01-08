#!/usr/bin/env python


import numpy as np
import sys,math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)


x    = data[:,0]
rho  = data[:,5]
rhof = data[:,6]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

nr = np.sum(rho)*dx
rho /= nr
nrf = np.sum(rhof)*dx
rhof /= nrf



Dkl_rho_rhof = 0. #np.sum(rho*np.log(rho/rhof))*dx
Dkl_rhof_rho = 0. #np.sum(rhof*np.log(rhof/rho))*dx
for i in range(space):
  if (rhof[i] > 1e-300) and (rho[i] > 1e-300):
    Dkl_rho_rhof += rhof[i] * np.log(rhof[i]/rho[i]) * dx
    Dkl_rhof_rho += rho[i]  * np.log(rho[i]/rhof[i]) * dx
    

print Dkl_rho_rhof,Dkl_rhof_rho