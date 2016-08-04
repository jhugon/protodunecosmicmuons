#!/usr/bin/env python

import re
import math

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)

class Particle(object):
  def __init__(self,line,iEvent=None):
    self.bad = False
    reIntParts = r"([0-9-]+)\s+"*6
    reFloatParts = r"([0-9.-]+)\s+"*9
    reParticleLine = r"^"+reIntParts+reFloatParts[:-3]+r"$"
    particleLineMatch = re.match(reParticleLine,line)
    if not particleLineMatch:
      self.bad = True
      return
    self.status = int(particleLineMatch.group(1))
    self.pdgid = int(particleLineMatch.group(2))
    self.px = float(particleLineMatch.group(7))
    self.py = -abs(float(particleLineMatch.group(8)))
    self.pz = float(particleLineMatch.group(9))
    self.p = math.sqrt(self.px**2+self.py**2+self.pz**2)
    self.e =  float(particleLineMatch.group(10))
    self.m =  float(particleLineMatch.group(11))
    self.x =  float(particleLineMatch.group(12))
    self.y =  float(particleLineMatch.group(13))
    self.z =  float(particleLineMatch.group(14))
    self.t =  float(particleLineMatch.group(15))
    self.thetaz = math.atan2((self.px**2+self.pz**2)**0.5,-self.py)
    self.phiz = math.atan2(self.pz,self.px)
    self.theta = math.atan2((self.px**2+self.py**2)**0.5,self.pz)
    self.phi = math.atan2(self.py,self.px)
    print self.thetaz, self.py, (self.px**2+self.pz**2)**0.5

  def __nonzero__(self):
    return not self.bad

  def __str__(self):
    if self:
      return "Particle: status: {} pdgid: {} px: {} py: {} pz: {} e: {} x: {} y: {} z: {} t: {}".format(self.status,self.pdgid,self.px,self.py,self.pz,self.e,self.x,self.y,self.z,self.t)
    else:
      return "Particle: Invalid"
  
def makeParticles(fname,limit=1000000000):
  print fname
  result = []
  f = open(fname)
  
  iEvent = None
  for iLine, line in enumerate(f):
    if iLine > limit*2:
      break
    eventLineMatch = re.match(r"^([0-9]+)\s+[0-9]+$",line)
    if eventLineMatch:
      iEvent = int(eventLineMatch.group(1))
    p = Particle(line,iEvent)
    if p:
      result.append(p)
  return result

if __name__ == "__main__":

  c = root.TCanvas()
  
  xb = root.TH1F("xb","",100,-1000,1000)
  yb = root.TH1F("yb","",100,-1000,1000)
  zb = root.TH1F("zb","",100,-1000,1000)
  thetaz = root.TH1F("thetaz","",30,0,90)
  phiz = root.TH1F("phiz","",30,-math.pi,math.pi)
  theta = root.TH1F("theta","",30,0,90)
  phi = root.TH1F("phi","",30,-math.pi,math.pi)
  
  xb2 = root.TH1F("xb2","",100,-1000,1000)
  yb2 = root.TH1F("yb2","",100,-1000,1000)
  zb2 = root.TH1F("zb2","",100,-1000,1000)
  thetaz2 = root.TH1F("thetaz2","",30,0,90)
  phiz2 = root.TH1F("phiz2","",30,-math.pi,math.pi)
  xb2.SetLineColor(root.kRed-1)
  yb2.SetLineColor(root.kRed-1)
  zb2.SetLineColor(root.kRed-1)
  thetaz2.SetLineColor(root.kRed-1)
  phiz2.SetLineColor(root.kRed-1)
  
  xb3 = root.TH1F("xb3","",100,-1000,1000)
  yb3 = root.TH1F("yb3","",100,-1000,1000)
  zb3 = root.TH1F("zb3","",100,-1000,1000)
  thetaz3 = root.TH1F("thetaz3","",30,0,90)
  phiz3 = root.TH1F("phiz3","",30,-math.pi,math.pi)
  xb3.SetLineColor(root.kBlue)
  yb3.SetLineColor(root.kBlue)
  zb3.SetLineColor(root.kBlue)
  thetaz3.SetLineColor(root.kBlue)
  phiz3.SetLineColor(root.kBlue)
  
  xb4 = root.TH1F("xb4","",100,-1000,1000)
  yb4 = root.TH1F("yb4","",100,-1000,1000)
  zb4 = root.TH1F("zb4","",100,-1000,1000)
  thetaz4 = root.TH1F("thetaz4","",30,0,90)
  phiz4 = root.TH1F("phiz4","",30,-math.pi,math.pi)
  xb4.SetLineColor(root.kCyan)
  yb4.SetLineColor(root.kCyan)
  zb4.SetLineColor(root.kCyan)
  thetaz4.SetLineColor(root.kCyan)
  phiz4.SetLineColor(root.kCyan)
  
  setHistTitles(xb,"Muon Starting x position [cm]","Events/bin")
  setHistTitles(yb,"Muon Starting y position [cm]","Events/bin")
  setHistTitles(zb,"Muon Starting z position [cm]","Events/bin")
  setHistTitles(thetaz,"Muon Starting zenith angle [degrees]","Events/bin")
  setHistTitles(phiz,"Muon Starting azimuthal angle with zenith [radians]","Events/bin")
  
  particles = makeParticles("jti3/AntiMuonCutEvents_1000000.txt",100)
  particles2 = makeParticles("jti3/DATAte",100)
  particles3 = makeParticles("MillionEvents_RawData_3/AntiMuonCutEvents_1000000.txt",100)
  particles4 = makeParticles("hamlet.txt",100)
  print 1
  for particle in particles:
    xb.Fill(particle.x)
    yb.Fill(particle.y)
    zb.Fill(particle.z)
    thetaz.Fill(particle.thetaz*180/math.pi)
    phiz.Fill(particle.phiz)
    theta.Fill(particle.theta*180/math.pi)
    phi.Fill(particle.phi)
    #print particle.thetaz, particle.px, particle.py, particle.pz
  print 2
  for particle in particles2:
    xb2.Fill(particle.x)
    yb2.Fill(particle.y)
    zb2.Fill(particle.z)
    thetaz2.Fill(particle.thetaz*180/math.pi)
    phiz2.Fill(particle.phiz)
  print 3
  for particle in particles3:
    xb3.Fill(particle.x)
    yb3.Fill(particle.y)
    zb3.Fill(particle.z)
    thetaz3.Fill(particle.thetaz*180/math.pi)
    phiz3.Fill(particle.phiz)
  print 4
  for particle in particles4:
    xb4.Fill(particle.x)
    yb4.Fill(particle.y)
    zb4.Fill(particle.z)
    thetaz4.Fill(particle.thetaz*180/math.pi)
    phiz4.Fill(particle.phiz)
  
  xb.Draw()
  xb2.Draw("same")
  xb3.Draw("same")
  xb4.Draw("same")
  c.SaveAs("xb.png")
  yb.Draw()
  yb2.Draw("same")
  yb3.Draw("same")
  yb4.Draw("same")
  c.SaveAs("yb.png")
  zb.Draw()
  zb2.Draw("same")
  zb3.Draw("same")
  zb4.Draw("same")
  c.SaveAs("zb.png")
  thetaz.Draw()
  thetaz2.Draw("same")
  thetaz3.Draw("same")
  thetaz4.Draw("same")
  c.SaveAs("thetaz.png")
  phiz.Draw()
  phiz2.Draw("same")
  phiz3.Draw("same")
  phiz4.Draw("same")
  c.SaveAs("phiz.png")
  theta.Draw()
  c.SaveAs("theta.png")
  phi.Draw()
  c.SaveAs("phi.png")
