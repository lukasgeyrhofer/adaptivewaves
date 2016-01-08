#!/usr/bin/env python


import numpy as np
import sys,math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u","--ufile")
parser.add_argument("-v","--speed",type=float,default=1e-5)
parser.add_argument("-E","--notextend",action='store_true',default=False)
#parser.add_argument("-V","--fitnessvariance",type=float,default=1e-5)
#parser.add_argument("-M","--mutationsigma",type=float,default=1e-3)
#parser.add_argument("-D","--mutationrate",type=float,default=1e-5)
args = parser.parse_args()


try:
  data = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open files"
  exit(1)


x = data[:,0]
u = data[:,1]

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]


if args.notextend:
  xall = x[:]
  uall = u[:]
else:  
  grid = np.arange(10*space)


  xpos = grid*dx+x[-1]+dx
  xneg = grid*dx+x[0]-len(grid)*dx


  xfitneg = x[:100]
  lufitneg = np.log(u[:100])

  expneg = (100*np.dot(xfitneg,lufitneg) - np.sum(xfitneg)*np.sum(lufitneg))/(100*np.dot(xfitneg,xfitneg)-np.sum(xfitneg)*np.sum(xfitneg))

  uneg = np.exp(expneg*(xneg-x[0]))*u[0]
  upos = 0.5*(xpos - args.speed/xpos)

  xall = np.concatenate((xneg,x,xpos))
  uall = np.concatenate((uneg,u,upos))


gaussall = np.exp(-xall*xall/(2*args.speed))/np.sqrt(2*math.pi*args.speed)

popsize = 1/(np.dot(uall,gaussall)*dx)

call = gaussall[:]*popsize

gall = uall*call


g = np.sum(gall)
xg = np.dot(xall,gall)/g
xxg = np.dot(xall*xall,gall)/g


for i in range(len(xall)):
  print >> sys.stderr,xall[i],call[i],uall[i],gall[i]

print popsize,xg,xxg-xg*xg,x[-1]



