#!/usr/bin/env python

import math
from scipy import *
from scipy.integrate import quad, dblquad

func = lambda x: x**4

class OnlineMeanVarianceCalc(object):
  def __init__(self):
    self.N = 0.
    self.mean = 0.
    self.M2 = 0.

  def addPoint(self,x):
    self.N += 1
    delta = x - self.mean
    self.mean += delta/self.N
    self.M2 += delta*(x-self.mean)

  def getMean(self):
    if self.N <1:
      return float('nan')
    else:
      return self.mean
    
  def getVariance(self):
    if self.N < 2:
      return float('nan')
    else:
      return self.M2 / (self.N - 1)
  def getStdDev(self):
    return math.sqrt(self.getVariance())

def mcInt(N,xmin,xmax):
  totalvolume = xmax-xmin

  nAccept = 0
  fvalSum = 0.
  meanVarCalc = OnlineMeanVarianceCalc()
  energies = []
  thetas = []
  while nAccept < N:
    x = rand()*(xmax-xmin)+xmin
    fval = func(x)
    fvalSum += fval
    meanVarCalc.addPoint(fval)
    nAccept += 1

  integralEstimate = totalvolume*fvalSum/N
  integralVariance = totalvolume**2*meanVarCalc.getVariance()/N
  integralUnc = math.sqrt(integralVariance)
  return integralEstimate, integralUnc

def directInt(xmin,xmax):
  result, uncertainty = quad(func,xmin,xmax) 
  return result, uncertainty

if __name__ == "__main__":
  xlims = [
    [0,1],
    [1,2],
    [2,3],
    [3,4],
  ]
  for xlim in xlims:
    mc, mcUnc = mcInt(40000,xlim[0],xlim[1])
    direct, directUnc = directInt(xlim[0],xlim[1])
    z = (mc - direct)/math.sqrt(directUnc**2 + mcUnc**2)
    print "{0:3} {1:3} {2:9.3g} {3:9.3g} {4:9.3g} {5:9.3g} {6:9.3g}".format(xlim[0],xlim[1],direct,mc,mc/direct,mcUnc,z)


