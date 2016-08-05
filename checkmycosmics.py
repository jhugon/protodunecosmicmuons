#!/usr/bin/env python

import math
from scipy import *
from mycosmics import mcInt
from mycosmics_differentialplots import getTotalFlux

####def mcInt(N,emin,emax,thetamin,thetamax):
####  return energies, thetas, integralEstimate, integralUnc
####def getTotalFlux(eMin,eMax,costhetamin,costhetamax):
####  """
####  In units of Hz cm^-2
####  """
####  return result

if __name__ == "__main__":

  eRanges = [
    [0.0001,1.],
    [1.,2.],
    [2.,3.],
    [2.5,3.5],
    [3,4],
    [4,5],
    [5,6],
    [6,7],
    [7,8],
    #[5,10],
    #[10,50],
    #[50,150],
    #[150,300],
    #[300,500],
    #[500,1000],
    #[1000,5000],
  ]

  thetaRanges = [
    [0,90.],
    #[0,70.],
    [0,45.],
    [45.,70],
  ]
    
  for thetaRange in thetaRanges:
    thetaRangeRad = [theta*pi/180 for theta in thetaRange]
    costhetaRange = [cos(theta) for theta in thetaRangeRad]
    for eRange in eRanges:
      outStr = "theta [deg]: {0:4.1f}, {1:4.1f}, energy [GeV]: {2:6.2f}, {3:6.2f} ".format(thetaRange[0],thetaRange[1],eRange[0],eRange[1])
      energies, thetas, fluxMC, fluxMCError = mcInt(1000,eRange[0],eRange[1],thetaRangeRad[0],thetaRangeRad[1])
      flux = getTotalFlux(eRange[0],eRange[1],costhetaRange[1],costhetaRange[0])
      outStr += "Direct Int Flux: {0:<8.3e}, MC Int Flux: {1:<8.3e} +/-  {2:<8.3g}, MC/Direct Ratio: {3:<5.3f}".format(flux,fluxMC,fluxMCError,fluxMC/flux)
      print outStr
