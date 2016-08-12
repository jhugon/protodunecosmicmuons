#!/usr/bin/env python

import math
from scipy import *
from scipy.integrate import dblquad

MUONMASS = 0.1134289267

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

def differentialFlux(energy,costheta):
  """
    dI
  -----
    dE

  in units of Hz m^-2 sr^-1 GeV^-1
  """
  p1 = 0.102573
  p2 = -0.068287
  p3 = 0.958633
  p4 = 0.0407253
  p5 = 0.817285
  
  costhetastar = sqrt((costheta**2+p1**2+p2*(costheta**p3)+p4*(costheta**p5))/(1.+p1**2+p2+p4))
  twopointseven = energy*(1+3.64/(energy*(costhetastar**1.29)))
  f = .140*(twopointseven**(-2.7))*((1./(1.+1.1*energy*costhetastar/115.))+0.054/(1.+1.1*energy*costhetastar/850.))
  return f

#######################################3
#
# dI/dx = f(x)
# I = Int f(x) dx
#
# draw x samples. find f(x) for the samples
# integral is then  V * sum(f(x)) / Nsamples,
# where V is the volume you are integrating 
# over in x.
#
#######################################3

def mcInt(N,emin,emax,thetamin,thetamax):
  #print N,emin,emax,thetamin,thetamax
  costhetamin = cos(thetamin)
  costhetamax = cos(thetamax)
  sinthetamin = sin(thetamin)
  sinthetamax = sin(thetamax)
  if costhetamax < costhetamin:
    tmp = costhetamax
    costhetamax = costhetamin
    costhetamin = tmp
  totalvolume = abs(costhetamax-costhetamin)*(emax-emin)
  fmax = max(differentialFlux(emin,costhetamin),differentialFlux(emin,costhetamax))

  nAccept = 0
  fvalSum = 0.
  energies = []
  thetas = []
  while nAccept < N:
    #theta = rand()*(thetamax-thetamin)+thetamin
    #costheta = cos(theta)
    costheta = rand()*(costhetamax-costhetamin)+costhetamin
    energy = rand()*(emax-emin)+emin
    fval = differentialFlux(energy,costheta)
    fvalSum += fval

    randfval = rand()*fmax
    if randfval <= fval:
      theta = math.acos(costheta)
      thetas.append(theta)
      energies.append(energy)
      nAccept += 1

  integralEstimate = totalvolume*fvalSum/N
  energies = array(energies)
  thetas = array(thetas)

  return energies, thetas, integralEstimate

class Muon(object):
  def __init__(self,x,y,z,t,px,py,pz,e):
    self.px = px
    self.py = py
    self.pz = pz
    self.p = (self.px**2+self.py**2+self.pz**2)**(0.5)
    self.e =  e
    self.m =  e**2-self.p**2
    self.x =  x
    self.y =  y
    self.z =  z
    self.t =  t
    self.thetaz = math.atan2((self.px**2+self.pz**2)**0.5,-self.py)
    self.phiz = math.atan2(self.pz,self.px)
    self.theta = math.atan2((self.px**2+self.py**2)**0.5,self.pz)
    self.phi = math.atan2(self.py,self.px)

  def __str__(self):
      return "Particle: px: {} py: {} pz: {} e: {} x: {} y: {} z: {} t: {}".format(self.px,self.py,self.pz,self.e,self.x,self.y,self.z,self.t)
  def hepevt(self):
    """
    <Status> <PDG ID> <1st Mother> <2nd Mother> <1st Daughter> <2nd Daughter> <Px> <Py> <Px> <E> <Mass> <x> <y> <z> <t>
    """
    result = "1 -13 0 0 0 0"
    result += " {p.px} {p.py} {p.pz} {p.e} {mass} {p.x} {p.y} {p.z} {p.t}".format(p=self,mass=MUONMASS)
    return result

def sample(N,emin,emax,thetamin,thetamax,xmin,xmax,ymin,ymax,zmin,zmax):
  m = MUONMASS
  assert(emin>=m)
  energies, thetas, integralEstimate = mcInt(N,emin,emax,thetamin,thetamax)
  xs = rand(N)*(xmax-xmin) + xmin
  ys = rand(N)*(ymax-ymin) + ymin
  zs = rand(N)*(zmax-zmin) + zmin
  ts = zeros(N)

  p2s = energies**2 - m**2
  ps = sqrt(p2s)

  phis = rand(N)*2*math.pi - math.pi

  pys = -ps*cos(thetas)
  pxzs = ps*sin(thetas)
  pxs = pxzs*cos(phis)
  pzs = pxzs*sin(phis)

  particles = []
  for i in range(N):
    tmp = Muon(xs[i],ys[i],zs[i],ts[i],pxs[i],pys[i],pzs[i],energies[i])
    particles.append(tmp)
  return particles, integralEstimate

if __name__ == "__main__":

  import argparse

  minenergy_default = MUONMASS
  maxenergy_default = 1000.
  mintheta_default = 0.
  maxtheta_default = 90.

  parser = argparse.ArgumentParser(description='Generate cosmic ray muon events in hepevt format.')
  parser.add_argument('outfilename',
                      help='Output file name')
  parser.add_argument('--nEvents','-n', type=int,
                      required=True,
                      help='Number of events to generate')
  parser.add_argument('--minenergy', type=float,
                      default=minenergy_default,
                      help='Minimum energy [GeV], default={0}'.format(minenergy_default))
  parser.add_argument('--maxenergy', type=float,
                      default=maxenergy_default,
                      help='Maximum energy [GeV], default={0}'.format(maxenergy_default))
  parser.add_argument('--mintheta', type=float,
                      default=mintheta_default,
                      help='Minimum zenith angle [degrees], default={0}'.format(mintheta_default))
  parser.add_argument('--maxtheta', type=float,
                      default=maxtheta_default,
                      help='Maximum zenith angle [degrees], default={0}'.format(maxtheta_default))
  parser.add_argument('--diagnostics','-d', action="store_true",
                      help='Make diagnostic plots')
  parser.add_argument('--debug','-D', action="store_true",
                      help='Print debug messages')
  parser.add_argument('--roottree','-R', action="store_true",
                      help='Make ROOT Tree')
  
  
  args = parser.parse_args()

  print "N events: {0}".format(args.nEvents)
  print "Energy range: {0:.3f}-{1:.3f} GeV".format(args.minenergy,args.maxenergy)
  print "Theta range: {0:.3f}-{1:.3f} degrees".format(args.mintheta,args.maxtheta)
  muons, integralEst = sample(args.nEvents,args.minenergy,args.maxenergy,args.mintheta*math.pi/180,args.maxtheta*math.pi/180,-1,1,-1,1,-1,1)
  print "Flux: {0:.3g} Hz".format(integralEst)
  print "Outputing HEPEVT file '{0}'".format(args.outfilename)
  with open(args.outfilename,'w') as outfile:
    for i, muon in enumerate(muons):
      outfile.write("{0} {1}".format(i,1))
      outfile.write('\n')
      outfile.write(muon.hepevt())
      outfile.write('\n')

  if args.debug:
    print len(muons)
    for muon in muons:
      print muon
      print muon.e, muon.thetaz

  if args.diagnostics:
    
    import ROOT as root
    from helpers import *
    root.gROOT.SetBatch(True)
    c = root.TCanvas()

    thetaHist = root.TH1F("theta","",10,0,90)
    phiHist = root.TH1F("phi","",30,-180,180)
    energyHist = root.TH1F("energy","",100,0,10)
    
    setHistTitles(thetaHist,"#theta_{zenith} [degrees]","Events/bin")
    setHistTitles(phiHist,"#phi_{azimuth} [degrees]","Events/bin")
    setHistTitles(energyHist,"E_{#mu} [GeV]","Events/bin")
    
    for muon in muons:
      thetaHist.Fill(muon.thetaz*180/math.pi)
      phiHist.Fill(muon.phiz*180/math.pi)
      energyHist.Fill(muon.e)
    
    phiHist.Draw()
    c.SaveAs("phiHist.png")
    c.SaveAs("phiHist.pdf")
    energyHist.Draw()
    c.SaveAs("energyHist.png")
    c.SaveAs("energyHist.pdf")
  
    thetaHistIntegral = thetaHist.Integral()
  
    thetaHist.Draw()
    c.SaveAs("thetaHist.png")
    c.SaveAs("thetaHist.pdf")

  if args.roottree:
    from makeRootTree import makeRootTree
    import os.path
    outfileNameRoot = os.path.splitext(args.outfilename)[0] + ".root"
    print "Outputing ROOT file '{0}' with tree name 'tree'".format(outfileNameRoot)
    makeRootTree(args.outfilename,outfileNameRoot,-1,False)
