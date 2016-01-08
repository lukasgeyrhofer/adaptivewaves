#!/usr/bin/env python


import numpy as np
import math
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("-c","--cfile")
parser.add_argument("-u","--ufile")
parser.add_argument("-M","--mutationexpdecay_p",type=float,default=1e-3)
parser.add_argument("-K","--mutationexpdecay_n",type=float,default=1e-3)
parser.add_argument("-k","--ratiodelmut",type=float,default=0.5)
parser.add_argument("-N","--normalize",action="store_true")
args = parser.parse_args()

try:
  urawdata = np.genfromtxt(args.ufile)
  crawdata = np.genfromtxt(args.cfile)
except:
  print >> sys.stderr,"could not open files"
  exit(1)
  
x = urawdata[:,0]
u = urawdata[:,1]
c = crawdata[:,1]

dx = x[1]-x[0]
space=len(x)
space0=(x*x).argmin()

xneg = x[x<0]
xpos = x[x>=0]
mkn = np.exp(xneg/args.mutationexpdecay_n)*args.ratiodelmut/args.mutationexpdecay_n
mkp = np.exp(-xpos/args.mutationexpdecay_p)/args.mutationexpdecay_p

mknorm = np.sum(mkn)+np.sum(mkp)

mk = np.concatenate((mkn,mkp))/mknorm

g = c*u

#for convolution we need the u array (3 times) as long
#extend it with u*(x) ~ exp(x \sigma) for negative x, and with u*(x) ~ x/2 for positive x
dropborders = 50
ucrop = u[dropborders:-dropborders]
ucropp = ucrop[-dropborders:]+ucrop[-1]-ucrop[-dropborders-1]
ucropm = np.exp((-x[dropborders]+x[:dropborders])/args.mutationexpdecay_p)*u[dropborders]
uextn = np.exp((x-x[-1]-dx)/args.mutationexpdecay_p)*ucropm[0]
xextp = x+(x[-1]-x[0])+dx
uextp = xextp*0.5
uext = np.concatenate((uextn,ucropm,ucrop,ucropp,uextp))
xext = np.concatenate((x-(x[-1]-x[0])-dx,x,xextp))


pi = np.zeros(space)
for i in range(space):
  pi[i] = np.dot(c,uext[space-space0+i:2*space-space0+i])


rhof = pi*mk

normalize = np.ones(6)
if args.normalize:
  normalize[0] = max(c)
  normalize[1] = max(u)
  normalize[2] = max(g)
  normalize[3] = max(pi)
  normalize[4] = max(mk)
  normalize[5] = max(rhof)
  

# get critical parameters

for i in range(space):
  if u[space-1-i] < 0.5*0.75*x[space-1-i]:
    xc = x[space-1-i]
    break

# mean fitness of g and \rho_f, respectively (fixating clones, fixating mutations)
xg = np.dot(x,g)/np.dot(np.ones(space),g)
xf = np.dot(x,rhof)/np.dot(np.ones(space),rhof)

vg = np.dot(x*x,g)/np.dot(np.ones(space),g)-xg**2

# maximum of L, maybe needed for laplace approximation (??)
xu = x[(x-2*u).argmax()]

# inverse coalescence time (??)
xn2 = np.dot(c,u*u)/np.dot(c,u)
xn3 = np.dot(c,u*u*u)/np.dot(c,u)
xn4 = np.dot(c,u*u*u*u)/np.dot(c,u)
xn5 = np.dot(c,u*u*u*u*u)/np.dot(c,u)


#population size
popsize = np.dot(c,np.ones(space))*dx

# 1 x
# 2 c
# 3 u*
# 4 g
# 5 pi ( = c conv u*)
# 6 mutkernel
# 7 rhof
for i in range(space):
  print >> sys.stderr,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e"%(x[i],c[i]/normalize[0],u[i]/normalize[1],g[i]/normalize[2],pi[i]/normalize[3],mk[i]/normalize[4],rhof[i]/normalize[5])


print "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e"%(xc,xg,xf,xu,xn2,xn3,xn4,xn5,vg,popsize)