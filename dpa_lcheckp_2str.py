#!/usr/bin/python

# Copyright (C) David Powell <david@drp.id.au>
#  This program comes with ABSOLUTELY NO WARRANTY; and is provided
#  under the GNU Public License v2
#
#
# Implementation of the DPA algorithm for linear (or affine)
# gap costs using check-pointing(1) to recover the
# alignment.  It has time complexity of O(n*n), and space
# complexity of O(n)
#
#1:  D. R. Powell, L. Allison and T. I. Dix,
#    "A Versatile Divide and Conquer Technique for Optimal
#    String Alignment", Information Processing Letters,
#    1999, 70:3, pp 127-139


import sys

BIG_NUM = 1<<30

matchCost    = 0
mismatchCost = 1
startGap     = 3
contGap      = 1

diag,vert,horz = (0,1,2)

def minIndex(a,b,c):
    if a<=b and a<=c: return 0
    elif        b<=c: return 1
    else:             return 2



def createRow(len):
    r = []
    for i in range(len):
      r.append( [State(), State(), State()] )
    return r
                

def cost(fromState, toState, c1, c2):
    if toState == diag:
      return (mismatchCost, matchCost)[c1 == c2]

    if toState == fromState:
      return contGap           # Same state, so just continue a gap

    return startGap + contGap  # Starting a new gap


class State:
    def __init__(self):
        self.reset()
        
    def reset(self):
        self.cost      = BIG_NUM
        self.fromCol   = BIG_NUM
        self.fromState = BIG_NUM

    def __str__(self):
        return "c=%d fromCol=%d fromState=%d\n"%(self.cost,self.fromCol,self.fromState)

    def minWithFrom(self, cell, toState, c1 = None, c2 = None):
        costs = []
        for fromState in range(3):
          costs.append( cell[fromState].cost + cost(fromState, toState, c1, c2) )

        minI = minIndex(*costs)
        self.cost      = costs[minI]
        self.fromCol   = cell[minI].fromCol
        self.fromState = cell[minI].fromState


class DPA_CP:
    def __init__(self, str1, str2):
        self.str1 = str1
        self.str2 = str2
        self.len1 = len(str1)
        self.len2 = len(str2)
        self.r1 = createRow(self.len2+1)
        self.r2 = createRow(self.len2+1)

    def getAlignment(self):
        (cost, a1, a2) = self.align(0, self.len1, 0, self.len2)

        #print cost
        #print a1[::-1]
        #print a2[::-1]
        a1 = list(a1)
        a1.reverse()
        a2 = list(a2)
        a2.reverse()
        s = map( lambda x,y : "<"+x+","+y+">", a1, a2)
        print ' ', ' '.join(s), "Edit cost = %d"%cost

    def align(self, i1, i2, j1, j2, sState = 0, eState = -1):
        #print "Doing (%d,%d) -> (%d,%d)   (s=%d -> e=%d)"%(i1,j1,i2,j2,sState,eState)
        r1 = self.r1
        r2 = self.r2

        CProw = (i1+i2)/2

        for i in range(i1, i2+1):
            (r1, r2) = (r2, r1)

            for j in range(j1, j2+1):
                c1 = (i==0 and '-' or self.str1[i-1])
                c2 = (j==0 and '-' or self.str2[j-1])
                
                if i==i1 and j==j1:
                    for s in range(3):
                      r2[j][s].cost = (BIG_NUM, 0)[ s == sState ]
                    continue

                if i==i1 or j==j1:
                    r2[j][diag].cost = BIG_NUM
                else:
                    r2[j][diag].minWithFrom(r1[j-1], diag, c1, c2)
                
                if i==i1:
                    r2[j][vert].cost = BIG_NUM
                else:
                    r2[j][vert].minWithFrom(r1[j], vert)
                
                if j==j1:
                    r2[j][horz].cost = BIG_NUM
                else:
                    r2[j][horz].minWithFrom(r2[j-1], horz)

                if i==CProw:
                    for s in range(3):
                        setattr( r2[j][s], 'fromCol',   j)
                        setattr( r2[j][s], 'fromState', s)

        if eState<0:
          eState = minIndex( r2[j][diag].cost, r2[j][vert].cost, r2[j][horz].cost)
        minCost = r2[j][eState].cost

        if (i2 - i1 > 1):
            jSplit = r2[j][eState].fromCol
            sSplit = r2[j][eState].fromState

            (_, strA, strB) = self.align(CProw, i2, jSplit, j2, sSplit, eState)
            (_, sA, sB)     = self.align(i1, CProw, j1, jSplit, sState, sSplit)
            strA += sA
            strB += sB
        else:
          (i,j,s) = (i2, j2, eState)
          strA = ''
          strB = ''
          while i!=i1 or j!=j1:
              if s == diag:
                  (c1, c2) = (self.str1[i-1], self.str2[j-1])
                  (iN,jN) = (i-1, j-1)
                  v2 = r2
                  v1 = r1

              if s == vert:
                  (c1, c2) = (self.str1[i-1], '-')
                  (iN,jN) = (i-1, j)
                  v2 = r2
                  v1 = r1

              if s == horz:
                  (c1, c2) = ('-', self.str2[j-1])
                  (iN,jN) = (i, j-1)
                  v1 = v2 = (r1,r2)[ i == i2]

              v = v2[j][s].cost

              sN = None
              for fromState in range(3):
                  if v == v1[jN][fromState].cost + cost( fromState, s, c1, c2):
                      sN = fromState
                      break

              (i,j,s) = (iN, jN, sN)

              strA += c1
              strB += c2

          if s != sState:
              print "BAD sState != s : %d %d"%(sState,s)
          
        return (minCost, strA, strB)


if len(sys.argv)!=1 and len(sys.argv)!=5:
    print "Usage: %s [matchCost mismatchCost startGap contGap]"%sys.argv[0]
    sys.exit(-1)

if  len(sys.argv)==5:
    matchCost    = int(sys.argv[1])
    mismatchCost = int(sys.argv[2])
    startGap     = int(sys.argv[3])
    contGap      = int(sys.argv[4])

print "Enter string A : ",
strA = sys.stdin.readline().rstrip("\n")
print
print "Enter string B : ",
strB = sys.stdin.readline().rstrip("\n")
print

d = DPA_CP(strA, strB)

d.getAlignment()
