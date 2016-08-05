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

  assert(costheta <= 1.)
  assert(costheta >= 0.)
  assert(energy >= 0.)
  
  p1 = 0.102573
  p2 = -0.068287
  p3 = 0.958633
  p4 = 0.0407253
  p5 = 0.817285
  
  costhetastar = sqrt((costheta**2+p1**2+p2*(costheta**p3)+p4*(costheta**p5))/(1.+p1**2+p2+p4))
  twopointseven = energy*(1+3.64/(energy*(costhetastar**1.29)))
  f = .140*(twopointseven**(-2.7))*((1./(1.+1.1*energy*costhetastar/115.))+0.054/(1.+1.1*energy*costhetastar/850.))
  return f

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
  assert(emax>emin)
  assert(thetamax>thetamin)
  costhetamin = cos(thetamin)
  costhetamax = cos(thetamax)
  assert(costhetamax<costhetamin)
  sinthetamin = sin(thetamin)
  sinthetamax = sin(thetamax)
  totalvolume = abs(costhetamax-costhetamin)*(emax-emin)
  fmax = differentialFlux(emin,costhetamin)

  nAccept = 0
  meanVarCalc = OnlineMeanVarianceCalc()
  energies = []
  thetas = []
  while nAccept < N:
    #costheta = rand()*(costhetamax-costhetamin)+costhetamin
    sintheta = rand()*(sinthetamax-sinthetamin)+sinthetamin
    theta = math.asin(sintheta)
    costheta = cos(theta)
    energy = rand()*(emax-emin)+emin
    fval = differentialFlux(energy,costheta)
    meanVarCalc.addPoint(fval)

    if fval >fmax:
      raise Error("Got a larger fval than fmax for energy: %s costheta %s",energy,costheta)
    randfval = rand()*fmax
    if randfval <= fval:
      thetas.append(theta)
      energies.append(energy)
      nAccept += 1

  integralEstimate = totalvolume*meanVarCalc.getMean()
  integralVariance = totalvolume**2*meanVarCalc.getVariance()/N
  integralUnc = math.sqrt(integralVariance)

  energies = array(energies)
  thetas = array(thetas)

  return energies, thetas, integralEstimate, integralUnc

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
  assert(emax>emin)
  assert(thetamax>thetamin)
  energies, thetas, integralEstimate, integralUnc = mcInt(N,emin,emax,thetamin,thetamax)
  xs = rand(N)*(xmax-xmin) + xmin
  ys = rand(N)*(ymax-ymin) + ymin
  zs = rand(N)*(zmax-zmin) + zmin
  ts = zeros(N)

  p2s = energies**2 - m**2
  ps = sqrt(p2s)

  phis = rand(N)*2*math.pi - math.pi

  pys = ps*cos(thetas)
  pxzs = ps*sin(thetas)
  pxs = pxzs*cos(phis)
  pzs = pxzs*sin(phis)

  particles = []
  for i in range(N):
    tmp = Muon(xs[i],ys[i],zs[i],ts[i],pxs[i],pys[i],pzs[i],energies[i])
    particles.append(tmp)
  return particles, integralEstimate

if __name__ == "__main__":

  muons, integralEst = sample(10000,MUONMASS,1000,0.000,70.*math.pi/180,-1,1,-1,1,-1,1)
  print integralEst
  print len(muons)
  #for muon in muons:
  #  print muon.e, muon.thetaz
  
  import ROOT as root
  from helpers import *
  root.gROOT.SetBatch(True)
  c = root.TCanvas()
  
  thetaHist = root.TH1F("theta","",30,0,90)
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
  energyHist.Draw()
  c.SaveAs("energyHist.png")

  thetaHistIntegral = thetaHist.Integral()

  cos2ThetaGraph = root.TGraph()
  cos2ThetaGraph.SetLineColor(root.kBlue)
  #funcNormalization = (thetaHistIntegral/math.pi/0.25)
  #funcNormalization = (thetaHistIntegral/(math.pi/4.))/(thetaHist.GetNbinsX())
  funcNormalization = thetaHist.GetBinContent(1)/cos(thetaHist.GetXaxis().GetBinCenter(1)*math.pi/180.)**2
  print funcNormalization
  print thetaHist.GetBinContent(1)
  print thetaHist.GetBinContent(1)/funcNormalization
  for i in range(1,thetaHist.GetNbinsX()):
    tmpTheta = thetaHist.GetXaxis().GetBinCenter(i)
    tmpCos2 = funcNormalization*cos(tmpTheta*math.pi/180.)**2
    cos2ThetaGraph.SetPoint(i-1,tmpTheta,tmpCos2) 

  thetaHist.Draw()
  cos2ThetaGraph.Draw("l")
  c.SaveAs("thetaHist.png")

  #############################################
  
  thetaIntHist = getIntegralHist(thetaHist,False)
  energyIntHist = getIntegralHist(energyHist,False)
  setHistTitles(thetaIntHist,thetaIntHist.GetXaxis().GetTitle(),"Events #geq X")
  setHistTitles(energyIntHist,energyIntHist.GetXaxis().GetTitle(),"Events #geq X")
  
  thetaIntHist.Draw()
  c.SaveAs("thetaIntHist.png")
  energyIntHist.Draw()
  c.SaveAs("energyIntHist.png")

  with open("cosmics.txt",'w') as outfile:
    for i, muon in enumerate(muons):
      outfile.write("{0} {1}".format(i,1))
      outfile.write('\n')
      outfile.write(muon.hepevt())
      outfile.write('\n')
