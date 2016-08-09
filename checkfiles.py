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

def normToBinWidth(hist):
  xaxis = hist.GetXaxis()
  nBins = xaxis.GetNbins()
  for i in range(1,nBins+1):
    binContent = hist.GetBinContent(i)
    binWidth = hist.GetBinWidth(i)
    hist.SetBinContent(i,binContent/binWidth)
  return hist

if __name__ == "__main__":

  c = root.TCanvas()
  fileConfigs = [
    {
        "fn": "jti3/AntiMuonCutEvents_1000000.txt",
        "title": "Original File",
        "color": root.kBlack,
    },
    {
        "fn": "hamlet.txt",
        "title": "SamplingProgram.C",
        "color": root.kBlue,
    },
    {
        "fn": "cosmics.txt",
        "title": "mycosmics.py",
        "color": root.kRed,
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
        "plottf1": "4700*cos(x/180*3.14159)*cos(x/180*3.14159)",
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
        #"binning": [1000,0,100],
        "binning": [100,array.array('f',getLogBins(100,0.1,1000))],
        "logy": True,
        "logx": True,
        "normToBinWidth":True
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
    if "normToBinWidth" in histConfig and histConfig["normToBinWidth"]:
      for hist in theseHists:
        normToBinWidth(hist)
    axisHist = makeStdAxisHist(theseHists,freeTopSpace=0.35)
    if 'logx' in histConfig and histConfig['logx']:
      c.SetLogx(True)
    else:
      c.SetLogx(False)
    if 'logy' in histConfig and histConfig['logy']:
      c.SetLogy(True)
      axisHist = makeStdAxisHist(theseHists,True,freeTopSpace=0.35)
    else:
      c.SetLogy(False)
    ylabel = "Events/bin"
    xlabel = histConfig['title']
    if "normToBinWidth" in histConfig and histConfig["normToBinWidth"]:
      ylabel = "Events/GeV"
    setHistTitles(axisHist,xlabel,ylabel)
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
  
