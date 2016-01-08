#!/usr/bin/env python


import numpy as np
import argparse
import sys,math

currentshift = 0


def reprod(nn,xx,time,speed):
  mutkernel = args.epsilon*args.mutationrate/args.mutationsigma*np.exp(-(space-np.arange(space))*dx/args.mutationsigma)
  speedprefactor = args.epsilon*args.wavespeed/dx
  twoepssqrt = np.sqrt(args.epsilon*2.)

  nn[0] = nn[-1] = 0

  for i in range(1,space-1):
    tmpn = nn[i]+args.epsilon*(xx[i]-time*speed-args.mutationrate)*nn[i] + np.dot(mutkernel[i-1::-1],nn[:i])
    if tmpn < 0:
      tmpn = 0
    else:
      tmpn += twoepssqrt*(np.random.poisson(tmpn)-tmpn)
    nn[i] = tmpn

  return nn

def populationcontrol(nn,uu):
  nn /= np.dot(nn,uu)
  return nn


def averagestuff(nn,xx):
  popsize = np.dot(np.ones(space),nn)
  meanfit = np.dot(xx,nn)/popsize
  varfit  = np.dot(xx*xx,nn)/popsize - meanfit**2
  return (meanfit,varfit,popsize)


def print_dens(nn,xx,ddx,time):
  for i in range(len(xx)):
    print >> sys.stderr,time,xx[i],nn[i]/ddx
  print >> sys.stderr


def update_u(uu,time,speed,ddx,sspace):
  baseshift = int(time*speed/ddx)
  fracshift = time*speed/ddx - 1.*baseshift
  tmpu = np.zeros(len(uu))
  for i in range(baseshift,sspace-1):
    tmpu[i] = fracshift * uu[i-baseshift+1] + (1-fracshift)*uu[i-baseshift]
  tmpu[sspace-1] = tmpu[sspace-2]+ddx/2.
  return tmpu


parser = argparse.ArgumentParser()
parser.add_argument("-c","--cfile",type=str)
parser.add_argument("-u","--ufile",type=str,required=True)
parser.add_argument("-v","--wavespeed",type=float,default=1e-5)
parser.add_argument("-D","--mutationrate",type=float,default=1e-5)
parser.add_argument("-M","--mutationsigma",type=float,default=1e-2)
parser.add_argument("-e","--epsilon",type=float,default=1e-2)
parser.add_argument("-S","--maxsteps",type=int,default=10000)
parser.add_argument("-O","--outputsteps",type=int,default=100)
parser.add_argument("-o","--densoutputsteps",type=int,default=1000)
args = parser.parse_args()

try:
  datau = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open ufile"
  exit(1)


x = datau[:,0]
u = datau[:,1]


try:
  datac = np.genfromtxt(args.cfile)
  c=datac[:,1]
except:
  c = np.exp(-x*x/(2*args.wavespeed))
  


dx = x[1] - x[0]
space = len(x)
space0 = (x*x).argmin()

n = c*dx

n = populationcontrol(n,u)

for i in range(1,args.maxsteps+1):
  n = reprod(n,x,i*args.epsilon,args.wavespeed)
  curu = update_u(u,i*args.epsilon,args.wavespeed,dx,space)
  n = populationcontrol(n,curu)
  if i%args.outputsteps == 0:
    (mf,vf,ps) = averagestuff(n,x-i*args.epsilon*args.wavespeed)
    print i*args.epsilon,mf,vf,ps
  #if i%args.densoutputsteps == 0:
  #  print_dens(n,x-i*args.epsilon*args.wavespeed,dx,i*args.epsilon)


for i in range(space):
  print >> sys.stderr,x[i]-args.maxsteps*args.epsilon*args.wavespeed,n[i]/dx






