#!/usr/bin/env python

deviationthreshold = 1e-2

from scipy.special import erfi
import scipy.optimize as optimization
import numpy as np
import math
import argparse
import sys

def ustar(x,ps,fv):
  return np.exp(x**2/(2.*fv))/(ps+np.sqrt(2*math.pi/fv)*erfi(x/np.sqrt(2.*fv)))

parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
parser.add_argument("-N","--populationsize",type=float,default=1e7)
parser.add_argument("-V","--fitnessvariance",type=float,default=1e-5)
parser.add_argument("-f","--fractionoflinearregime",type=float,default=-1.)
parser.add_argument("-s","--outputstring",type=str,default="")

args = parser.parse_args()

try:
  data = np.genfromtxt(args.inputfile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

xdata = data[:,0]
ydata = data[:,1]

try:
  startindex = np.where(ydata < 0)[-1][-1]+1
except:
  startindex = (xdata*xdata).argmin()
stopindex = len(xdata)-1

if 1 > args.fractionoflinearregime > 0:
  xcindex = startindex
  for i in range(len(xdata)-1,startindex,-1):
    if ydata[i]<0.5*xdata[i]-deviationthreshold:
      xcindex = i
      break
  if xcindex-startindex > 0:
    stopindex = int( (-args.fractionoflinearregime*startindex + xcindex)/(1-args.fractionoflinearregime))

xused = xdata[startindex:stopindex]
yused = ydata[startindex:stopindex]

sigma = np.ones(len(xused))
startvalues = (args.populationsize,args.fitnessvariance)

parameters,covariance = optimization.curve_fit(ustar,xused,yused,startvalues,sigma)

yapprox = ustar(xdata,parameters[0],parameters[1])

for i in range(len(xdata)):
  print >> sys.stderr,"%9.6lf\t%21.14e\t%21.14e"%(xdata[i],ydata[i],yapprox[i])

print >> sys.stdout,"%s%21.14e\t%21.14e"%(args.outputstring,parameters[0],parameters[1])
