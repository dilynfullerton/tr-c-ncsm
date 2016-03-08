#!/usr/bin/python
"""
$ python GetSmallerInteraction.py <fname> <n1> <n2>

Reads the TBME file and removes elements beyond the
cutoff given by n1. 
NOTE: These are printed out, so the
result should be piped to the desired file. Then,
the final line needs to be moved to the top.
"""

from sys import argv

fname = argv[1]
n1,n2 = [int(x) for x in argv[2:4]]
f = open(fname)

line=f.readline()
ldat = line.split()
ntbme_old,n1_old,n2_old = [int(x) for x in ldat[:3]]
hw,Atbme = [float(x) for x in ldat[3:]]

index_max = (n1+1)*(n1+2)/2
def GetN(indx):
  n = 0
  while (n+1)*(n+2)/2 < indx: 
	n+=1
  return n

print '%d  %d  %d  %.4f  %.4f' % (ntbme_old,n1,n2,hw,Atbme)

ntbme = 0

for line in f:
  ldat = line.split()
  a,b,c,d,J,T = [ int(x) for x in ldat[:6] ]
  if a>index_max or b>index_max or c>index_max or d>index_max: 
	continue
  if GetN(a)+GetN(b) > n2 or GetN(c)+GetN(d) > n2: 
	continue
  print line.rstrip('\n')
  ntbme += 1

print ntbme



