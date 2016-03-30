#!/usr/bin/env python


import numpy as np
import math
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("-c","--cfile")
parser.add_argument("-u","--ufile")
parser.add_argument("-s","--mutationexpdecay",type=float,default=1)
parser.add_argument("-U","--mutationrate",type=float,default=0.1)
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
rhoneg = np.zeros(len(xneg))
rhopos = np.exp(-xpos/args.mutationexpdecay)/args.mutationexpdecay
rho = np.concatenate((rhoneg,rhopos))
mu = args.mutationrate/args.mutationexpdecay
mrhopos = np.exp(-(1-mu)*xpos/args.mutationexpdecay)/args.mutationexpdecay
mrho = np.concatenate((rhoneg,mrhopos))


g = c*u

#for convolution we need the u array (3 times) as long
#extend it with u*(x) ~ exp(x \sigma) for negative x, and with u*(x) ~ x/2 for positive x
dropborders = 50
ucrop = u[dropborders:-dropborders]
ucropp = ucrop[-dropborders:]+ucrop[-1]-ucrop[-dropborders-1]
ucropm = np.exp((-x[dropborders]+x[:dropborders])/args.mutationexpdecay)*u[dropborders]
uextn = np.exp((x-x[-1]-dx)/args.mutationexpdecay)*ucropm[0]
xextp = x+(x[-1]-x[0])+dx
uextp = xextp*0.5
uext = np.concatenate((uextn,ucropm,ucrop,ucropp,uextp))
xext = np.concatenate((x-(x[-1]-x[0])-dx,x,xextp))

rhoext = np.concatenate((np.zeros(space),rho,np.exp(-(np.arange(space)*dx+x[-1]+dx)/args.mutationexpdecay)/args.mutationexpdecay))

conv_rho_u = np.zeros(space)
conv_rho_c = np.zeros(space)
conv_rho_uc = np.zeros(space)
pi = np.zeros(space)

for i in range(space):
  pi[i] = np.dot(c,uext[space-space0+i:2*space-space0+i])
  conv_rho_u[i] = np.dot(rho,uext[space-space0+i:2*space-space0+i])
  conv_rho_c[i] = np.dot(c,rhoext[i+space+space0:i+space0:-1])
  conv_rho_uc[i] = np.dot(u*c,rhoext[i+space+space0:i+space0:-1])

rhof   = pi*rho
rhomf  = pi*mrho
rhobg  = conv_rho_u*c
rhowin = conv_rho_c*u


dg = rhowin - rhobg

normalize = np.ones(10)
if args.normalize:
  normalize[0] = max(c)
  normalize[1] = max(u)
  normalize[2] = max(g)
  normalize[3] = max(pi)
  normalize[4] = max(rho)
  normalize[5] = max(rhof)
  normalize[6] = max(rhobg)
  normalize[7] = max(rhowin)
  normalize[8] = max(conv_rho_uc)
  normalize[9] = max(rhomf)
else:
  normalize[4] = np.sum(rho)*dx
  normalize[5] = np.sum(rhof)*dx
  normalize[6] = np.sum(rhobg)*dx
  normalize[7] = np.sum(rhowin)*dx
  normalize[8] = np.sum(conv_rho_uc)*dx
  normalize[9] = np.sum(rhomf)*dx
maxg = max(g)  

# get crossover
xcd   = x[(np.diff(u)).argmax()]
xc90  = max(x[(u-0.90*.5*x)<0])
xc95  = max(x[(u-0.95*.5*x)<0])
xc2nd = x[np.diff(u,n=2).argmin()+1]



xc=0
for i in range(space):
  if u[space-1-i] < 0.5*0.9*x[space-1-i]:
    xc = x[space-1-i]
    break

# mean fitness of g and \rho_f, respectively (fixating clones, fixating mutations)
xg = np.dot(x,g)/np.dot(np.ones(space),g)
xf = np.dot(x,rhof)/np.dot(np.ones(space),rhof)

# width of g
vg = np.dot(x*x,g)/np.dot(np.ones(space),g)-xg**2
startindexg = ((x-xg)*(x-xg)).argmin()
foundlower = 0
foundupper = 0
i = 0
while foundlower + foundupper < 2:
  if foundlower == 0 and g[startindexg - i] < 0.5*maxg:
    foundlower = 1
    glower = startindexg - i
  if foundupper == 0 and g[startindexg + i] < 0.5*maxg:
    foundupper = 1
    gupper = startindexg + i
  i += 1
fwhmg = x[gupper] - x[glower]
    
  


# maximum of L, maybe needed for laplace approximation (??)
xu = x[(x-2*u).argmax()]

# inverse coalescence time (??)
idcu = 1./np.dot(c,u)
xn2 = np.dot(c,u*u)*idcu
xn3 = np.dot(c,u*u*u)*idcu
xn4 = np.dot(c,u*u*u*u)*idcu
xn5 = np.dot(c,u*u*u*u*u)*idcu


#population size
popsize = np.dot(c,np.ones(space))*dx

#properties of density
csum = np.dot(c,np.ones(space))
xdens = np.dot(x,c)/csum
x2dens = np.dot(x*x,c)/csum
vardens = x2dens - xdens**2


#properties of change in g
xdgma = x[dg.argmax()]
xdgmi = x[dg.argmin()]
xwin = np.dot(x,rhowin)/np.dot(rhowin,np.ones(space))
xbg  = np.dot(x,rhobg)/np.dot(rhobg,np.ones(space))
xmf  = np.dot(x,rhomf)/np.dot(rhomf,np.ones(space))

for i in range(space):
  print >> sys.stderr,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e"%(x[i],c[i]/normalize[0],u[i]/normalize[1],g[i]/normalize[2],pi[i]/normalize[3],rho[i]/normalize[4],rhof[i]/normalize[5],rhobg[i]/normalize[6],rhowin[i]/normalize[7],conv_rho_uc[i]/normalize[8],rhomf[i]/normalize[9])
#                                                                                          1    2                 3                 4                 5                  6                   7                    8                     9                      10                          11

print "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e"%(xc90,xc95,xcd,xc2nd,xg,xf,xu,xn2,xn3,xn4,xn5,vg,fwhmg,popsize,xdens,vardens,xdgma,xdgmi,xwin,xbg,xmf)
#                                                                                                                                      1    2    3   4     5  6  7  8   9   10  11  12 13    14      15    16      17    18    19   20  21

