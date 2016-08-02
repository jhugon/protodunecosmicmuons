#!/usr/bin/env python

from scipy import *
from scipy.integrate import dblquad

#Double_t fluxformula(Double_t *x, Double_t *par)
#{
#  //x = energy of the muons, par = zenith angle
#  Double_t xx = x[0];
#  Double_t p1 = 0.102573;
#  Double_t p2 = -0.068287;
#  Double_t p3 = 0.958633;
#  Double_t p4 = 0.0407253;
#  Double_t p5 = 0.817285;
#  Double_t costheta = cos((par[0]*TMath::Pi)/180);
#  Double_t workaround = TMath::Power(costheta,p5);
#  //the term on the above line was acting up when put into the parstar formula
#
#  Double_t parstar = TMath::Sqrt((costheta*costheta+p1*p1+p2*TMath::Power(costheta,p3)+p4*workaround)/(1+p1*p1+p2+p4));
#
#  Double_t twopointseven = xx*(1+3.64/(xx*TMath::Power(parstar,1.29)));
#  //please excuse my naming conventions, but they help me to remember!
#
#  Double_t f = .140*TMath::Power(twopointseven,-2.7)*((1/(1+1.1*xx*parstar/115))+0.054/(1+1.1*xx*parstar/850));
#  return f;
#}

def differentialFlux(energy,theta):
  """
    dI
  -----
  dE d(costheta)

  in units of Hz m^-2 sr^-1 GeV^-1
  """
  p1 = 0.102573
  p2 = -0.068287
  p3 = 0.958633
  p4 = 0.0407253
  p5 = 0.817285
  costheta = cos(theta)
  
  costhetastar = sqrt((costheta**2+p1**2+p2*(costheta**p3)+p4*(costheta**p5))/(1.+p1**2+p2+p4))
  twopointseven = energy*(1+3.64/(energy*(costhetastar**1.29)))
  f = .140*(twopointseven**(-2.7))*((1./(1.+1.1*energy*costhetastar/115.))+0.054/(1.+1.1*energy*costhetastar/850.))
  return f

#######################################3
#
# dI/dx = f(x)
# I = Int f(x) dx
#
# Find f max
#
# draw x & fsample. See which fsample values are less than f(x)
# integral is then frac less than f(x) * fmax*(xdomain)
#
#######################################3

def sample(N,emin,emax,thetamin,thetamax):
  fmax = max(differentialFlux(emin,thetamin),differentialFlux(emin,thetamax))
  totalvolume = fmax*(thetamax-thetamin)*(emax-emin)
  thetas = rand(N)*(thetamax-thetamin)+thetamin
  energies = rand(N)*(emax-emin)+emin
  ftrials = rand(N)*fmax

  fvals = differentialFlux(energies,thetas)

  passes = ftrials < fvals

  integralEstimate = totalvolume*len(passes)/N

  return energies[passes], thetas[passes], integralEstimate

energies, thetas, integralEst = sample(10000,0.1,10,0.01,0.99)
print integralEst
print len(energies)
print energies, thetas

#from matplotlib import pyplot as mpl
#
#fig, ax = mpl.subplots()
#
#thetas = linspace(0,pi)
#costhetas = cos(thetas)
#
#for energy in arange(1,11):
#  f = differentialFlux(energy,costhetas)
#  ax.plot(thetas*180/pi,f,label="%s GeV" % energy) 
#
#ax.plot(thetas*180/pi,cos(thetas)**2,'k--')
#
#ax.legend()
#fig.savefig("test.png")
#
#fig, ax = mpl.subplots()
#
#energies = logspace(-2,1)
#f = differentialFlux(energies,1.)
#ax.plot(energies,f) 
#fig.savefig("test2.png")
#
#import ROOT as root
#from helpers import *
#root.gROOT.SetBatch(True)
#
#thetaHist = root.TH1F("theta","",100,0,180)
#
#theta
