#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import os.path as op


def d(y,dx = 1e-2):
    return 0.5 * np.diff(np.concatenate([np.array([2*y[0] - y[1]]),y]) + np.concatenate([y,np.array([2*y[-1] - y[-2]])]))/dx


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-v","--speedfile")
parser.add_argument("-m","--mutationrate",type=float,default=None)
parser.add_argument("-s","--prependstring",default=None)
parser.add_argument("-o","--outfile",default="dpot0")
args = parser.parse_args()


try:
    speedlist = {}
    fpspeed = open(args.speedfile,"r")
    for line in fpspeed.readlines():
        values = line.split()
        if len(values) > 1:
            if not speedlist.has_key(values[0]):
                speedlist[values[0]] = {}
            speedlist[values[0]][values[1]] = float(values[2])
    fpspeed.close()
except:
    print >> sys.stderr,"could not open speed file"
    exit(1)

for fn in args.infiles:
    #print fn
    try:
        data = np.genfromtxt(fn)
        speed = speedlist[fn[1:3]][fn[5:7]]
    except:
        continue



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

    
    dpot0r = d(pot0r)
    outfilename = op.split(fn)[0] + "/" + args.outfile
    outfile = open(outfilename,"w")
    #print op.split(fn),outfilename
    #exit(1)
    #for i in range(space):
        #print >>outfile,x[i],pot0r[i],dpot0r[i]
    #outfile.close()

    positionsextrema = np.zeros(10)
    countextrema = 0
    for i in range(space-1):
        if dpot0r[i] * dpot0r[i+1] < 0:
            positionsextrema[countextrema] = (abs(dpot0r[i]) * x[i] + abs(dpot0r[i+1]) * x[i+1]) / (abs(dpot0r[i]) + abs(dpot0r[i+1]))
            countextrema+=1
    extrema = np.interp(positionsextrema[:countextrema],x,pot0r)
    print "{} {} {:10.6f}".format(fn[1:3],fn[5:7],speed),
    for i in range(countextrema):
        print "{:10.6f} {:15.6e} ".format(positionsextrema[i],extrema[i]),
    print

    #exit(1)

    #a = (pot0r[2:] - pot0r[1:space-1]) * (pot0r[1:space-1] - pot0r[:space-2])
    #b = (pot1r[2:] - pot1r[1:space-1]) * (pot1r[1:space-1] - pot1r[:space-2])


    #print args.prependstring,args.speed,args.mutationrate,
    #for i in np.where(a<0)[0]:
        #print "{} {}".format(x[i+1],pot0r[i+1]),
    #print

    #for i in np.where(b<0)[0]:
        #print "{} {}".format(x[i+1],pot1r[i+1]),
    #print
