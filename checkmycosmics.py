#!/usr/bin/env python

import math
from scipy import *
from mycosmics import mcInt, MUONMASS
from mycosmics_differentialplots import getTotalFlux, dIdtheta_only

####def mcInt(N,emin,emax,thetamin,thetamax):
####  """
####  args:
####    N: n events to generate
####    emin, emax: muon energy bounds in GeV
####    thetamin,thetamax: azimuthal angle bounds in rad
####  returns:
####    muon energies in GeV (ndarray)
####    muon thetas in rad (ndarray)
####    integral flux estimate in Hz cm^-2
####  """
####  return energies, thetas, integralEstimate
####def dIdtheta_only(costheta,eMin,eMax):
####  """
####  In units of Hz cm^-2 sr^-1
####  """
####  return result
####def getTotalFlux(eMin,eMax,costhetamin,costhetamax):
####  """
####  In units of Hz cm^-2
####  """
####  return result

if __name__ == "__main__":

  IforCostheta0 = dIdtheta_only(1,0,10000)
  IforCostheta0 *= 100**2
  print("Flux in Hz m^-2 sr^-1 for vertical muons: {0:.3f}, should be ~70".format(IforCostheta0))
  
  print "##############################################################"

  eRanges = [
    #[0.0001,1.],
    #[1.,2.],
    #[2.,3.],
    #[2.5,3.5],
    #[3,4],
    #[4,5],
    #[5,6],
    #[6,7],
    #[7,8],
    #[5,10],
    #[10,50],
    #[50,150],
    #[150,300],
    #[300,500],
    #[500,1000],
    #[1000,5000],

    [MUONMASS,100],
    [1,100],
  ]

  thetaRanges = [
    [0,90.],
    [0,70.],
    [85,90],
  ]
    
  print("Flux in cm^-2 min^-1, should be ~1 for all angles and E>1GeV")
  for thetaRange in thetaRanges:
    thetaRangeRad = [theta*pi/180 for theta in thetaRange]
    costhetaRange = [cos(theta) for theta in thetaRangeRad]
    for eRange in eRanges:
      outStr = "theta [deg]: {0:4.1f}, {1:4.1f}, energy [GeV]: {2:6.2f}, {3:6.2f} ".format(thetaRange[0],thetaRange[1],eRange[0],eRange[1])
      energies, thetas, fluxMC = mcInt(10000,eRange[0],eRange[1],thetaRangeRad[0],thetaRangeRad[1])
      flux = getTotalFlux(eRange[0],eRange[1],costhetaRange[1],costhetaRange[0])
      outStr += "Direct Int Flux: {0:<8.3f}, MC Int Flux: {1:<8.3f}, MC/Direct Ratio: {2:<5.3f}".format(flux*60,fluxMC*60,fluxMC/flux)
      print outStr
