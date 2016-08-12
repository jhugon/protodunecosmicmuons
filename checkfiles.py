#!/usr/bin/env python

import re
import math

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)

from rootpy.plotting import Hist, Hist2D, Canvas, Graph
from rootpy.tree import Tree
from rootpy.io import root_open

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
        "fn": "testcosmics.root",
        "title": "testcosmics.root",
        "color": root.kRed,
    },
    #{
    #    "fn": "jti3/AntiMuonCutEvents_1000000.root",
    #    "title": "Original File",
    #    "color": root.kBlack,
    #},
    #{
    #    "fn": "jti3/AntiMuonCutEvents_1000000.root",
    #    "title": "Original File, E>4 GeV",
    #    "color": root.kRed,
    #    "cuts": "E > 5",
    #},
    #{
    #    "fn": "jti3/AntiMuonCutEvents_1000000.root",
    #    "title": "Original File, 2 GeV < E < 4 GeV",
    #    "color": root.kGreen+1,
    #    "cuts": "E > 2 && e < 4", 
    #},
    #{
    #    "fn": "jti3/AntiMuonCutEvents_1000000.root",
    #    "title": "Original File, E<2 GeV",
    #    "color": root.kBlue,
    #    "cuts": "E < 2",
    #},
    #{
    #    "fn": "hamlet.root",
    #    "title": "SamplingProgram.C",
    #    "color": root.kBlue,
    #},
    #{
    #    "fn": "cosmics.root",
    #    "title": "mycosmics.py",
    #    "color": root.kRed,
    #},
    #{
    #    "fn": "cosmics.root",
    #    "title": "mycosmics.py, E>4 GeV",
    #    "color": root.kRed,
    #    "cuts": "E > 5",
    #},
    #{
    #    "fn": "cosmics.root",
    #    "title": "mycosmics.py, 2 GeV < E < 4 GeV",
    #    "color": root.kGreen+1,
    #    "cuts": "E > 2 && e < 4", 
    #},
    #{
    #    "fn": "cosmics.root",
    #    "title": "mycosmics.py, E<2 GeV",
    #    "color": root.kBlue,
    #    "cuts": "E < 2",
    #},
  ]
  histConfigs = [
    {
        "name": "xb",
        "var": "x",
        "title": "Muon starting x position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "yb",
        "var": "y",
        "title": "Muon starting y position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "zb",
        "var": "z",
        "title": "Muon starting z position [cm]",
        "binning": [100,-1000,1000],
    },
    {
        "name": "phiz",
        "var": "phiz*180/pi",
        "title": "Muon #phi with respect to zenith [degrees]",
        "binning": [30,-90,90],
    },
    {
        "name": "thetaz",
        "var": "thetaz*180/pi",
        "title": "Muon #theta with respect to zenith [degrees]",
        "binning": [45,0.,90.],
        "plottf1": "4700*cos(x/180*3.14159)*cos(x/180*3.14159)",
        "titletf1": "cos^{2}(#theta)",
        "colortf1": root.kGreen+1,
    },
    {
        "name": "thetazNorm",
        "var": "thetaz*180/pi",
        "title": "Muon #theta with respect to zenith [degrees]",
        "binning": [45,0.,90.],
        "plottf1": "0.05*cos(x/180*3.14159)*cos(x/180*3.14159)",
        "titletf1": "cos^{2}(#theta)",
        "colortf1": root.kGreen+1,
        "normalize": True
    },
    {
        "name": "phi",
        "var": "phi*180/pi",
        "title": "Muon #phi with respect to beam direction [degrees]",
        "binning": [30,-90,90],
    },
    {
        "name": "theta",
        "var": "theta*180/pi",
        "title": "Muon #theta with respect to beam direction [degrees]",
        "binning": [30,0.,180],
    },
    {
        "name": "energy",
        "var": "E",
        "title": "Muon starting energy [GeV]",
        #"binning": [1000,0,100],
        "binning": getLogBins(100,0.1,1000),
        "logy": True,
        "logx": True,
        "normToBinWidth":True
    },
  ]

  hists = []
  for histConfig in histConfigs:
    hists.append([])
  for iFile,fileConfig in enumerate(fileConfigs):
    for iHist,histConfig in enumerate(histConfigs):
      histArgs = histConfig['binning']
      hist = None
      if len(histArgs) == 3:
        hist = Hist(*histArgs)
      else:
        hist = Hist(histArgs)
      hist.UseCurrentStyle()
      hist.SetLineColor(fileConfig['color'])
      hist.SetMarkerStyle(0)
      hists[iHist].append(hist)
    with root_open(fileConfig['fn']) as infile:
      tree = infile.Get("tree")
      for iHist,histConfig in enumerate(histConfigs):
        var = histConfig['var']
        cuts = ""
        if "cuts" in histConfig:
          cuts = histConfig['cuts']
        tree.Draw(var,selection=cuts,hist=hists[iHist][iFile])

  for iHist,histConfig in enumerate(histConfigs):
    theseHists = [ hists[iHist][iFile] for iFile in range(len(fileConfigs))]
    if "normToBinWidth" in histConfig and histConfig["normToBinWidth"]:
      for hist in theseHists:
        normToBinWidth(hist)
    elif "normalize" in histConfig and histConfig["normalize"]:
      for hist in theseHists:
        normalizeHist(hist)
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
    elif "normalize" in histConfig and histConfig["normalize"]:
      ylabel = "Normalized Events/bin"
    setHistTitles(axisHist,xlabel,ylabel)
    axisHist.Draw("")
    for hist in theseHists:
      hist.Draw("histsame")
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
  
