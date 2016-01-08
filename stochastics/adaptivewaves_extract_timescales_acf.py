#!/usr/bin/env python

import numpy as np
import sys,math
import argparse
from scipy.optimize import curve_fit

def acf(time,taudec,tauosc):
  return np.exp(-time/taudec)*np.cos(2.*math.pi*time/tauosc)

def acf_phase(time,taudec,tauosc,phi):
  return np.exp(-time/taudec)*np.cos(2.*math.pi*(time+phi)/tauosc)/np.cos(2.*math.pi*phi/tauosc)


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",required=True)
parser.add_argument("-P","--phase",default=False,action="store_true")
parser.add_argument("-M","--maxtime",default=-1,type=float)
args = parser.parse_args()

try:
  data = np.genfromtxt(args.infile)
except:
  print >> sys.stderr,"could not open file"
  exit(1)


t0 = data[:,0]
a0 = data[:,1]

if args.maxtime > 0:
  t = t0[t0<args.maxtime]
  a = a0[t0<args.maxtime]
else:
  t = t0
  a = a0

if args.phase:
  startvalues = np.array((200,200,1))
  popt,pcov = curve_fit(acf_phase,t,a,p0=startvalues)
else:
  startvalues = np.array((200,200))
  popt,pcov = curve_fit(acf,t,a,p0=startvalues)


print popt[0],popt[1],
if len(popt) > 2:
  print popt[2],
else:
  print 0,
print pcov[0,1]/np.sqrt(pcov[0,0]*pcov[1,1])
