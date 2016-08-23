#!/usr/bin/env python

import re
import math

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)

from scipy import *
import matplotlib.pyplot as mpl

from rootpy.plotting import Hist, Hist2D, Canvas, Graph
from rootpy.tree import Tree
from rootpy.io import root_open

from throughgoingstudy import Projector

if __name__ == "__main__":
  c = root.TCanvas()

  minx = -5
  maxx = 5
  miny = -5
  maxy = 5
  minz = -5
  maxz = 5

  projector = Projector(minx,maxx,miny,maxy,minz,maxz)

  infilename = "cosmicstheta70xyz5.root"
  sampleRate = 169.89438884
  sampleArea = 1. # meters^2
  detectorArea = (maxx-minx)*(maxz-minz)/100.**2 # meters^2
  detectorRate = sampleRate*detectorArea/sampleArea
  print "Sample: "+infilename
  print "Sample Rate: {} Hz".format(sampleRate)
  print "Sample Area: {} m^2".format(sampleArea)
  print "Detector Rate: {} Hz".format(detectorRate)
  print "Detector Area: {} m^2".format(detectorArea)
  
  infile = root.TFile(infilename)
  tree = infile.Get("tree")

  nEventsTotal = 0
  nEventsEntered = 0

  thetaVthetax = Hist2D(18,0,180,18,0,180)
  
  for iEvent in range(tree.GetEntries()):
    if iEvent > 10000:
      break
      pass
    tree.GetEntry(iEvent)
    nEventsTotal += 1
    ## Getting entry/exit points ordered by z
    minpoint, maxpoint = projector.getExtremaZPoints(tree)
    if minpoint:
      nEventsEntered += 1
      theta = tree.theta*180/pi
      p = sqrt(tree.px**2+tree.py**2+tree.pz**2)
      thetax = math.acos(tree.px/p)*180/pi
      thetaVthetax.Fill(thetax,theta)
    
    #print("Event: {:4} x y z: {:6.1f} {:6.1f} {:6.1f} px py pz E: {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(iEvent,tree.x,tree.y,tree.z,tree.px,tree.py,tree.pz,tree.E))
    #print(" {}    {}".format(minpointZ,maxpointZ))
    #print(" {}    {}".format(hitMin,hitMax))

  eventWeight = sampleRate/float(nEventsTotal)
  print("eventWeight: {} Hz".format(eventWeight))
  print("Rate going through detector: {} Hz".format(nEventsEntered*eventWeight))


  setupCOLZFrame(c)
  thetaVthetax.Scale(eventWeight)
  setHistTitles(thetaVthetax,"#theta w.r.t. x-axis [deg]","#theta w.r.t. z-axis [deg]")
  thetaVthetax.Draw("colz")
  drawStandardCaptions(c,"Going through 10 #times 10 #times 10 cm^{3} voxel")
  c.SaveAs("thetaVthetax.png")
  
  
