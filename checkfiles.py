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
    #print self.thetaz, self.py, (self.px**2+self.pz**2)**0.5

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

  fileConfigs = [
    {
        "fn": "hamlet_1.txt",
        "title": "E < 1 GeV",
        "color": root.kCyan,
    },
    {
        "fn": "hamlet_2.txt",
        "title": "E < 2 GeV",
        "color": root.kMagenta,
    },
    {
        "fn": "hamlet_3.txt",
        "title": "E < 3 GeV",
        "color": root.kRed,
    },
    {
        "fn": "hamlet_4.txt",
        "title": "E < 4 GeV",
        "color": root.kGreen+1,
    },
    {
        "fn": "hamlet_5.txt",
        "title": "E < 5 GeV",
        "color": root.kBlue,
    },
  ]
  histConfigs = [
    {
        "name": "xb",
        "attr": "x",
        "title": "Muon starting x position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "yb",
        "attr": "y",
        "title": "Muon starting y position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "zb",
        "attr": "z",
        "title": "Muon starting z position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "phiz",
        "attr": "phiz",
        "scalexby": 180/math.pi,
        "title": "Muon #phi with respect to zenith [degrees]",
        "binning": [30,-90,90],
    },
    {
        "name": "thetaz",
        "attr": "thetaz",
        "scalexby": 180/math.pi,
        "title": "Muon #theta with respect to zenith [degrees]",
        "binning": [45,0.,90.],
        "plottf1": "450*cos(x/180*3.14159)*cos(x/180*3.14159)",
        "titletf1": "cos^{2}(#theta)",
        "colortf1": root.kGreen+1,
    },
    {
        "name": "phi",
        "attr": "phi",
        "scalexby": 180/math.pi,
        "title": "Muon #phi with respect to beam direction [degrees]",
        "binning": [30,-90,90],
    },
    {
        "name": "theta",
        "attr": "theta",
        "scalexby": 180/math.pi,
        "title": "Muon #theta with respect to beam direction [degrees]",
        "binning": [30,0.,180],
    },
    {
        "name": "energy",
        "attr": "e",
        "title": "Muon starting energy [GeV]",
        "binning": [100,0,10],
        #"binning": [100,array.array('f',getLogBins(100,0.1,1000))],
        "logy": True,
        #"logx": True,
    },
  ]

  hists = {}
  for histConfig in histConfigs:
    hists[histConfig['name']] = {}
  for fileConfig in fileConfigs:
    for histConfig in histConfigs:
      histArgs = [histConfig['name']+"_"+fileConfig['fn'],""] + histConfig['binning']
      hists[histConfig['name']][fileConfig['fn']] = root.TH1F(*histArgs)
      hists[histConfig['name']][fileConfig['fn']].SetLineColor(fileConfig['color'])
    particles = makeParticles(fileConfig['fn'],100000)
    for particle in particles:
      for histConfig in histConfigs:
        x = getattr(particle,histConfig['attr'])
        if "scalexby" in histConfig:
          x *= histConfig["scalexby"]
        hists[histConfig['name']][fileConfig['fn']].Fill(x)


  for histConfig in histConfigs:
    theseHists = [ hists[histConfig['name']][fileConfig['fn']] for fileConfig in fileConfigs]
    axisHist = makeStdAxisHist(theseHists,freeTopSpace=0.45)
    if 'logx' in histConfig and histConfig['logx']:
      c.SetLogx(True)
    else:
      c.SetLogx(False)
    if 'logy' in histConfig and histConfig['logy']:
      c.SetLogy(True)
      axisHist = makeStdAxisHist(theseHists,True,freeTopSpace=0.45)
    else:
      c.SetLogy(False)
    setHistTitles(axisHist,histConfig['title'],"Events/bin")
    axisHist.Draw()
    for hist in theseHists:
      hist.Draw("same")
    ### plot tf1 time ###
    plots = []
    if 'plottf1' in histConfig:
      tmpPlot = root.TF1(
                            "plot"+str(random.randint(1000,1000000)),histConfig['plottf1'],
                            axisHist.GetXaxis().GetBinLowEdge(1),
                            axisHist.GetXaxis().GetBinUpEdge(axisHist.GetNbinsX())
                        )
      tmpPlot.SetLineStyle(2)
      if 'colortf1' in histConfig:
        tmpPlot.SetLineColor(histConfig["colortf1"])
      tmpPlot.Draw('LSAME')
      plots.append(tmpPlot)
    ### Legend time
    theseLabels = [ fileConfig['title'] for fileConfig in fileConfigs]
    theseHistsForLeg = [x for x in theseHists]
    if 'titletf1' in histConfig and len(plots)> 0:
      theseLabels.append(histConfig['titletf1'])
      theseHistsForLeg.append(plots[0])
    leg = drawNormalLegend(theseHistsForLeg,theseLabels)
    leg.Draw()
    ###
    c.RedrawAxis()
    c.SaveAs(histConfig['name']+".png")
    c.SaveAs(histConfig['name']+".pdf")
  
#
#  cos2Graph = root.TGraph()
#  cos2Graph.SetLineStyle(2)
#  Npoints = 1000
#  for i in range(Npoints):
#    cos2Graph.SetPoint(i,i*90./Npoints,math.cos(i*(math.pi/2.)/Npoints)**2*7000)
#  
#  xb.GetYaxis().SetRangeUser(0,14000)
#  xb.Draw()
#  #xb2.Draw("same")
#  #xb3.Draw("same")
#  xb4.Draw("same")
#  xb5.Draw("same")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("xb.png")
#  c.SaveAs("xb.pdf")
#  yb.GetYaxis().SetRangeUser(0,2e5)
#  yb.Draw()
#  #yb2.Draw("same")
#  #yb3.Draw("same")
#  yb4.Draw("same")
#  yb5.Draw("same")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("yb.png")
#  c.SaveAs("yb.pdf")
#  zb.GetYaxis().SetRangeUser(0,14000)
#  zb.Draw()
#  #zb2.Draw("same")
#  #zb3.Draw("same")
#  zb4.Draw("same")
#  zb5.Draw("same")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("zb.png")
#  c.SaveAs("zb.pdf")
#  thetaz.GetYaxis().SetRangeUser(0,8000)
#  thetaz.Draw()
#  #thetaz2.Draw("same")
#  #thetaz3.Draw("same")
#  thetaz4.Draw("same")
#  thetaz5.Draw("same")
#  cos2Graph.Draw("l")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("thetaz.png")
#  c.SaveAs("thetaz.pdf")
#  phiz.GetYaxis().SetRangeUser(0,7000)
#  phiz.Draw()
#  #phiz2.Draw("same")
#  #phiz3.Draw("same")
#  phiz4.Draw("same")
#  phiz5.Draw("same")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("phiz.png")
#  c.SaveAs("phiz.pdf")
#  theta.Draw()
#  c.RedrawAxis()
#  c.SaveAs("theta.png")
#  c.SaveAs("theta.pdf")
#  phi.Draw()
#  c.RedrawAxis()
#  c.SaveAs("phi.png")
#  c.SaveAs("phi.pdf")
#  c.SetLogy()
#  c.SetLogx()
#  energy.Draw()
#  energy4.Draw("same")
#  energy5.Draw("same")
#  leg.Draw()
#  c.RedrawAxis()
#  c.SaveAs("energy.png")
#  c.SaveAs("energy.pdf")
#
#  energy27 = energy4.Clone("energy24")
#  for iBin in range(1,energy27.GetNbinsX()+1):
#    x = energy27.GetBinCenter(iBin)
#    y = energy27.GetBinContent(iBin)
#    energy27.SetBinContent(iBin,y*x**2.7)
#  energy27.Draw()
#  c.SaveAs("energy27.png")
#  c.SaveAs("energy27.pdf")
