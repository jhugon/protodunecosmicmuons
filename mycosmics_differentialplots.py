#!/usr/bin/env python

from math import acos,asin,atan2
from scipy import *
from scipy.integrate import quad, dblquad
from matplotlib import pyplot as mpl
from mycosmics import differentialFlux, MUONMASS

def dIdE_only(energy,costhetaMin,costhetaMax):
  """
  In units of Hz cm^-2 GeV^-1
  """
  assert(costhetaMin<costhetaMax)
  result, uncertainty = quad(lambda x: differentialFlux(energy,x),costhetaMin,costhetaMax)
  result *= 2*pi
  uncertainty *= 2*pi
  return result

def dIdtheta_only(costheta,eMin,eMax):
  """
  In units of Hz cm^-2 sr^-1
  """
  assert(eMin<eMax)
  result, uncertainty = quad(lambda x: differentialFlux(x,costheta),eMin,eMax)
  return result

def getTotalFlux(eMin,eMax,costhetamin,costhetamax):
  """
  In units of Hz cm^-2
  """
  assert(costhetamin<costhetamax)

  # y,x -> energy, costheta
  result, uncertainty = dblquad(differentialFlux,costhetamin,costhetamax,lambda x: eMin, lambda x: eMax)
  result *= 2*pi
  uncertainty *= 2*pi
  return result

if __name__ == "__main__":

  fig, ax = mpl.subplots()
  
  thetas = linspace(0,pi/2.)
  costhetas = cos(thetas)
  
  for energy in arange(1,11):
    f = [differentialFlux(energy,costheta) for costheta in costhetas]
    ax.plot(thetas*180/pi,f/f[0],label="%s GeV" % energy) 
  ax.plot(thetas*180/pi,cos(thetas)**2,ls="--",label=r"cos($\theta$)^2") 
  #ax.plot(thetas[:len(thetas)/2]*180/pi,1./cos(thetas[:len(thetas)/2]),label=r"sec($\theta$)") 
  ax.legend()
  ax.set_xlabel(r"$\theta$ [deg]")
  ax.set_ylabel(r"dI/(d$\theta$ dE) [Hz cm$^{-2}$ sr$^{-1}$ GeV$^{-1}$]")
  fig.savefig("dIdtheta_multiE.png")

  ###########################################################
  
  energies = logspace(0,4)
  fig, ax = mpl.subplots()
  
  for thetaDeg in [0,56,75]:
    costheta = math.cos(thetaDeg*pi/180)
    f = [differentialFlux(energy,costheta) for energy in energies]
    f = array(f)
    f *= energies**2.7
    ax.loglog(energies,f,label=r"$\theta$ = %s deg" % thetaDeg)
  
  ax.legend()
  ax.set_xlim(7,7000)
  ax.set_ylim(0.003,0.2)
  ax.set_xlabel(r"E [GeV]")
  ax.set_ylabel(r"dI/(d$\theta$ dE) [Hz cm$^{-2}$ sr$^{-1}$ GeV$^{-1}$]")
  fig.savefig("dIdE_multiTheta.png")
  
  ###########################################################
  
  costhetamin = 0.03
  costhetamax = 1.
  dIdE = [dIdE_only(energy,costhetamin,costhetamax)*2*pi for energy in energies]
  
  fig, ax = mpl.subplots()
  ax.loglog(energies,dIdE)
  ax.set_xlabel("E [GeV]")
  ax.set_ylabel("dI/dE [Hz cm$^{-2}$ GeV$^{-1}$]")
  ax.text(0.95,0.95,r"{0:.2f} < cos($\theta$) < {1:.2f}".format(costhetamin,costhetamax),transform=ax.transAxes,verticalalignment="top",horizontalalignment="right")
  fig.savefig("dIdE_thetaInt.png")
  
  ###########################################################
  
  emin = MUONMASS
  emax = 1000.
  
  dIdtheta = [dIdtheta_only(costheta,emin,emax) for costheta in costhetas]
  
  fig, ax = mpl.subplots()
  ax.plot(thetas*180/pi,dIdtheta,label="Formula")
  ax.plot(thetas*180/pi,cos(thetas)**2*dIdtheta[0],"r--",label= r"cos($\theta$)$^2$")
  ax.set_xlabel(r"$\theta$ [deg]")
  ax.set_ylabel(r"dI/d$\theta$ [Hz cm$^{-2}$ sr$^{-1}$]")
  ax.legend()
  ax.text(0.95,0.80,r"$m_\mu c^{{-2}}$ < E < {1:.2f} GeV".format(emin,emax),transform=ax.transAxes,verticalalignment="top",horizontalalignment="right")
  fig.savefig("dIdtheta_EInt.png")
  
  ###########################################################
  
  fig, ax = mpl.subplots()
  energies = logspace(0,3)
  for thetadeg in [0.,70.]:
    f = [differentialFlux(energy,cos(thetadeg*pi/180)) * energies**2.7 for energy in energies]
    ax.loglog(energies,f,label=r"$\theta$ = %s deg" % thetadeg)
  ax.set_xlabel("E [GeV]")
  ax.set_ylabel("E$^{2.7}$ dI/dE [Hz cm$^{-2}$ sr$^{-1}$ GeV$^{1.7}$]")
  ax.set_xlim(1,2000)
  ax.set_ylim(0.002,0.2)
  ax.legend()
  fig.savefig("dIdE_logx.png")
  
  ##############################################################
  
  thetamin=0.0
  thetamax=1.
  emin=0.
  emax=1000.
  I = getTotalFlux(emin,emax,thetamin,thetamax)
  print("I = {0:10.4f} Hz cm^-2 = {5} cm^-2/min for theta in {1:.2f}, {2:.2f} and E in {3:.2f}, {4:.2f} GeV".format(I,acos(thetamax)*180/pi,acos(thetamin)*180/pi,emin,emax,I*60))

  thetamin=0.
  thetamax=1.
  emin=1.
  emax=1000.
  I = getTotalFlux(emin,emax,thetamin,thetamax)
  print("I = {0:10.4f} Hz cm^-2 = {5} cm^-2/min for theta in {1:.2f}, {2:.2f} and E in {3:.2f}, {4:.2f} GeV".format(I,acos(thetamax)*180/pi,acos(thetamin)*180/pi,emin,emax,I*60.))

  ##############################################################

  thetamin=cos(5.*pi/180.)
  thetamax=1.
  energyBins = linspace(0.,1000.,21)
  binWidth = energyBins[1]-energyBins[0]
  Is = []
  for iBin in range(len(energyBins)-1):
    Is.append(10**4*getTotalFlux(energyBins[iBin],energyBins[iBin+1],thetamin,thetamax))
  fig, ax = mpl.subplots()
  ax.bar(energyBins[:-1],Is,binWidth,label=r"$\theta$ = 0 to 5 deg".format(0,5))
  ax.plot(energyBins[1:],1e3*energyBins[1:]**(-1.7),"g-")
  ax.set_yscale('log')
  ax.set_title(r"For $\theta$ = 0 to 5 deg".format(0,5))
  ax.set_xlabel("E [GeV]")
  ax.set_ylabel("I [Hz m$^{{-2}}$ per {0:.1f} GeV]".format(binWidth))
  #ax.set_xlim(1,2000)
  #ax.set_ylim(0.002,0.2)
  fig.savefig("I_1000.png")
  
  energyBins = linspace(0.,10.,101)
  binWidth = energyBins[1]-energyBins[0]
  Is = []
  for iBin in range(len(energyBins)-1):
    Is.append(10**4*getTotalFlux(energyBins[iBin],energyBins[iBin+1],thetamin,thetamax))
  fig, ax = mpl.subplots()
  ax.bar(energyBins[:-1],Is,binWidth,label=r"$\theta$ = 0 to 5 deg".format(0,5))
  ax.plot(energyBins[1:],2e-3*energyBins[1:]**(-2.7),"g-")
  ax.set_yscale('linear')
  ax.set_title(r"For $\theta$ = 0 to 5 deg".format(0,5))
  ax.set_xlabel("E [GeV]")
  ax.set_ylabel("I [Hz m$^{{-2}}$ per {0:.1f} GeV]".format(binWidth))
  #ax.set_xlim(1,2000)
  #ax.set_ylim(0.002,0.2)
  fig.savefig("I_10.png")
  
  energyBins = linspace(0.,2.,21)
  binWidth = energyBins[1]-energyBins[0]
  Is = []
  for iBin in range(len(energyBins)-1):
    Is.append(10**4*getTotalFlux(energyBins[iBin],energyBins[iBin+1],thetamin,thetamax))
  fig, ax = mpl.subplots()
  ax.bar(energyBins[:-1],Is,binWidth,label=r"$\theta$ = 0 to 5 deg".format(0,5))
  ax.plot(energyBins[1:],2e-3*energyBins[1:]**(-1.7),"g-")
  ax.set_yscale('linear')
  ax.set_title(r"For $\theta$ = 0 to 5 deg".format(0,5))
  ax.set_xlabel("E [GeV]")
  ax.set_ylabel("I [Hz m$^{{-2}}$ per {0:.1f} GeV]".format(binWidth))
  #ax.set_xlim(1,2000)
  #ax.set_ylim(0.002,0.2)
  fig.savefig("I_2.png")
