#!/usr/bin/env python

# calculates X_c in u* profiles
# various methods are implemented, chosen via option "-m METHOD"
# METHOD == 90		x-value where u* profile is less than 90% of the asymptotic solution X/2
# METHOD == 95		x-value where u* profile is less than 95% of the asymptotic solution X/2
# METHOD == diff	x-value where first derivative of u* is maximal


import numpy as np
import sys,math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",default="-")
parser.add_argument("-m","--method",choices=["90","95","diff","2nd"],default="diff")
args = parser.parse_args()

if (args.infile == "-"):
  f = sys.stdin
else:
  f = args.infile
  
try:
  data = np.genfromtxt(f)
except:
  print >> sys.stderr,"could not open file"
  exit(1)

try:
  x = data[:,0]
  u = data[:,1]
except:
  print >> sys.stderr,"not enough columns in data-file"
  exit(1)

if args.method == "diff":
  xc = x[(np.diff(u)).argmax()]
elif args.method == "90":
  xc = max(x[(u-0.90*.5*x)<0])
elif args.method == "95":
  xc = max(x[(u-0.95*.5*x)<0])
elif args.method == "2nd":
  xc = x[np.diff(u,n=2).argmin()+1]
  #for i in range(1,len(u)-1):
    #print >> sys.stderr,x[i],np.diff(u,n=2)[i-1]
  
else:
  print >> sys.stderr,"could not find method to obtain xc"
  exit(1)

print xc