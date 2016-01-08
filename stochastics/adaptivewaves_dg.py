#!/usr/bin/env python


import numpy as np
import sys,math
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-c","--cfile",required=True)
parser.add_argument("-u","--ufile",required=True)
parser.add_argument("-m","--mutationmodel",choices=["diffusion","expdecay","jump","dexp"],default="diffusion")
parser.add_argument("-D","--mutationrate",type=float,default=1.)
parser.add_argument("-M","--mutationoption",type=float,default=1.)
args = parser.parse_args()


try:
  cdata = np.genfromtxt(args.cfile)
  udata = np.genfromtxt(args.ufile)
except:
  print >> sys.stderr,"could not open files"
  exit(1)
  
  
x = udata[:,0]
u = udata[:,1]
c = cdata[:,1]

assert len(u) == len(c)

space = len(x)
space0 = (x*x).argmin()
dx = x[1] - x[0]

g = u*c

if args.mutationmodel == "diffusion":
  uext = np.zeros(space+2)
  uext[1:space+1] = u
  uext[0] = u[0]*u[0]/u[1]
  uext[space+1] = 2*u[space-1]-u[space-2]

  cext = np.zeros(space+2)
  cext[1:space+1] = c
  cext[0] = c[0]*c[0]/c[1]
  cext[space+1] = c[space-1]*c[space-1]/c[space-2]

  dg = args.mutationrate*(u*np.diff(cext,n=2) - c*np.diff(uext,n=2))/(dx*dx)
  
elif args.mutationmodel == "expdecay":
  dropborders = 50
  ucrop = u[dropborders:-dropborders]
  ucropp = ucrop[-dropborders:]+ucrop[-1]-ucrop[-dropborders-1]
  ucropm = np.exp((-x[dropborders]+x[:dropborders])/args.mutationoption)*u[dropborders]
  uextn = np.exp((x-x[-1]-dx)/args.mutationoption)*ucropm[0]
  xextp = x+(x[-1]-x[0])+dx
  uextp = xextp*0.5
  uext = np.concatenate((uextn,ucropm,ucrop,ucropp,uextp))
  xext = np.concatenate((x-(x[-1]-x[0])-dx,x,xextp))
  
  xneg = x[x<0]
  xpos = x[x>=0]
  rhoneg = np.zeros(len(xneg))
  rhopos = np.exp(-xpos/args.mutationoption)/args.mutationoption
  rho = np.concatenate((rhoneg,rhopos))
  
  rhoext = np.concatenate((np.zeros(space),rho,np.exp(-(np.arange(space)*dx+x[-1]+dx)/args.mutationoption)/args.mutationoption))
  
  conv_rho_u = np.zeros(space)
  conv_rho_c = np.zeros(space)
  conv_rho_uc = np.zeros(space)
  pi = np.zeros(space)

  for i in range(space):
    pi[i] = np.dot(c,uext[space-space0+i:2*space-space0+i])*dx
    conv_rho_u[i] = np.dot(rho,uext[space-space0+i:2*space-space0+i])*dx
    conv_rho_c[i] = np.dot(c,rhoext[i+space+space0:i+space0:-1])*dx
    conv_rho_uc[i] = np.dot(u*c,rhoext[i+space+space0:i+space0:-1])*dx

  rhof   = pi*rho
  rhobg  = conv_rho_u*c
  rhowin = conv_rho_c*u
  
  dg = args.mutationrate*(rhowin - rhobg)

elif args.mutationmodel == "jump":
  ushift = np.zeros(space)
  ushift[0:space-int(args.mutationoption)] = u[int(args.mutationoption):space]
  ushift[space-int(args.mutationoption)] = np.arange(int(args.mutationoption))*0.5*dx+u[space-1]
  
  cshift = np.zeros(space)
  cshift[int(args.mutationoption):space] = c[0:space-int(args.mutationoption)]
  # c below the jump-width is anyway almost zero...
  
  dg = args.mutationrate * (u*cshift - ushift*c)

else:
  print >> sys.stderr,"mutationmodel not yet implemented"
  exit(1)


xmi = x[dg.argmin()]
xma = x[dg.argmax()]

mi = dg[dg.argmin()]
ma = dg[dg.argmax()]

xg = np.dot(x,g)/np.sum(g)

for i in range(space):
  print >>sys.stderr, x[i],g[i],dg[i]
  
print "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e"%(xmi,xma,xma-xmi,mi,ma,ma-mi,xg)

