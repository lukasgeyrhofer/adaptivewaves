#!/usr/bin/env python

import sys


popsize = {}




def switch_level(v1,r1,v2,r2,levelc):
  p1 = popsize[v1][r1]
  p2 = popsize[v2][r2]
  if ((p1 < levelc) and (p2 > levelc)) or ((p1 > levelc) and (p2 < levelc)):
    retval = True
    #print "check",p1,p2
  else:
    retval = False
  return retval

def get_level(v1,r1,v2,r2,levelc):
  l1 = popsize[v1][r1]
  l2 = popsize[v2][r2]
  
  diff_l1_l2 = l2-l1
  diff_l1_level = levelc-l1
  try:
    ratio = diff_l1_level/diff_l1_l2
  except:
    ratio = 0
  if not switch_level(v1,r1,v2,r2,levelc):
    print >> sys.stderr,"level not switched"
    exit(1)
  return (v1+(float(v2)-float(v1))*ratio,r1+(float(r2)-float(r1))*ratio,level)

  
  
def get_next_step(v1,r1,v2,r2,levelc):
  (dv,dr)=(v2-v1,r2-r1)
  
  if (dv,dr) == (1,0):
    print >> sys.stderr,"(1,0)",v2,r2
    if switch_level(v2,r2+1,v2+1,r2+1,levelc):
      return (v2,r2,v2,r2+1,v2,r2+1,v2+1,r2+1)
    elif switch_level(v2+1,r2+1,v2+1,r2,levelc):
      return (v2,r2,v2+1,r2,v2+1,r2+1,v2+1,r2)
    elif switch_level(v2,r2,v2+1,r2,levelc):
      return (v2,r2,v2,r2-1,v2,r2,v2+1,r2)
    else:
      return (0,0,0,0,0,0,0,0)
  elif (dv,dr) == (0,1):
    print >> sys.stderr,"(0,1)",v2,r2
    if switch_level(v2,r2,v2,r2+1,levelc):
      return (v2,r2,v2-1,r2,v2,r2,v2,r2+1)
    elif switch_level(v2,r2+1,v2+1,r2+1,levelc):
      return (v2,r2,v2,r2+1,v2,r2+1,v2+1,r2+1)
    elif switch_level(v2+1,r2+1,v2+1,r2,levelc):
      return (v2,r2,v2+1,r2,v2+1,r2+1,v2+1,r2)
    else:
      return (0,0,0,0,0,0,0,0)
  elif (dv,dr) == (-1,0):
    print >> sys.stderr,"(-1,0)",v2,r2
    if switch_level(v2,r2,v2+1,r2,levelc):
      return (v2,r2,v2,r2-1,v2,r2,v2+1,r2)
    elif switch_level(v2,r2,v2,r2+1,levelc):
      return (v2,r2,v2-1,r2,v2,r2,v2,r2+1)
    elif switch_level(v2,r2+1,v2-1,r2+1,levelc):
      return (v2,r2,v2,r2+1,v2,r2+1,v2-1,r2+1)
    else:
      return (0,0,0,0,0,0,0,0)
  elif (dv,dr) == (0,-1):
    print >> sys.stderr,"(0,-1)",v2,r2
    #print v2+1,r2,v2+1,r2-1
    if switch_level(v2+1,r2,v2+1,r2+1,levelc):
      return (v2,r2,v2+1,r2,v2+1,r2+1,v2+1,r2)
    elif switch_level(v2,r2,v2+1,r2,levelc):
      return (v2,r2,v2,r2-1,v2,r2,v2+1,r2)
    elif switch_level(v2,r2,v2,r2+1,levelc):
      return (v2,r2,v2-1,r2,v2,r2,v2,r2+1)
    else:
      return (0,0,0,0,0,0,0,0)
  else:
    return (0,0,0,0,0,0,0,0)
       
  
  
fp = open(sys.argv[1])

try:
  level = float(sys.argv[2])
except:
  print >> sys.stderr,"need 2 parameters: FILENAME LEVEL"
  exit(1)



for line in fp.readlines():
  values = line.split()
  if len(values) == 3:
    vv = float(values[0])
    rr = float(values[1])
    v = int(round(100.*vv/0.1))
    r = int(round(100.*rr/0.0025))
    #print v,r,vv,rr
    try:
      n = float(values[2])
      if not popsize.has_key(v):
	popsize[v] = {}
      popsize[v][r] = n
    except:
      continue

#print popsize      



# get level curve

r = 100
while not switch_level(0,r,0,r-1,level):
  #print popsize[v][r], popsize[v][r-1]
  r = r-1

(pv,pr,pn) = get_level(0,r,0,r-1,level)
print "%.10lf %.10lf %lf"%(float(pv)*0.001,float(pr)*0.000025,pn)

offlattice_lastv = -1
offlattice_lastr = r-1

offlattice_newv = 0
offlattice_newr = r-1

while (offlattice_lastv,offlattice_lastr,offlattice_newv,offlattice_newr) != (0,0,0,0):
  (offlattice_lastv,offlattice_lastr,offlattice_newv,offlattice_newr,v1,r1,v2,r2) = get_next_step(offlattice_lastv,offlattice_lastr,offlattice_newv,offlattice_newr,level)
  (v,r,n) = get_level(v1,r1,v2,r2,level)
  print "%.10lf %.10lf %lf"%(float(v)*0.001,float(r)*0.000025,level)





# ================================
# get minimum curve
# ================================

#for vv in range(101):
  #v = float(vv)/1000.
  #min_for_v_N = 1.e20
  #min_for_v_r = 0;
  #for rr in range(101):
    #r = float(rr)/100.*0.0025
    #try:
      #if popsize[vv][rr] < min_for_v_N:
	#min_for_v_N = popsize[vv][rr]
	#min_for_v_r = rr
	##print "found min", r
    #except:
      #continue
	
      
  #print vv,min_for_v_r,min_for_v_N
    
  
  
      
#print popsize	