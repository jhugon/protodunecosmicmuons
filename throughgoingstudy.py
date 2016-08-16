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

class Projector(object):

  def __init__(self,minx,maxx,miny,maxy,minz,maxz):
    """
    creates Projector object
    arguments are boundaries of box to project onto
    """
    self.minx = minx
    self.miny = miny
    self.minz = minz
    self.maxx = maxx
    self.maxy = maxy
    self.maxz = maxz

  def projectX(self,tree,xplane):
    # first y coord
    my = tree.py/tree.px
    by = tree.y - my*tree.x
    y = my*xplane + by

    # then z coord
    mz = tree.pz/tree.px
    bz = tree.z - mz*tree.x
    z = mz*xplane + bz
    return [xplane,y,z]
    
  def projectY(self,tree,yplane):
    # first x coord
    mx = tree.px/tree.py
    bx = tree.x - mx*tree.y
    x = mx*yplane + bx

    # then z coord
    mz = tree.pz/tree.py
    bz = tree.z - mz*tree.y
    z = mz*yplane + bz
    return [x,yplane,z]
    
  def projectZ(self,tree,zplane):
    # first x coord
    mx = tree.px/tree.pz
    bx = tree.x - mx*tree.z
    x = mx*zplane + bx

    # then y coord
    my = tree.py/tree.pz
    by = tree.y - my*tree.z
    y = my*zplane + by
    return [x,y,zplane]

  def projectXBounds(self,tree):
    return self.projectX(tree,self.minx), self.projectX(tree,self.maxx)

  def projectYBounds(self,tree):
    return self.projectY(tree,self.miny), self.projectY(tree,self.maxy)

  def projectZBounds(self,tree):
    return self.projectZ(tree,self.minz), self.projectZ(tree,self.maxz)

  def hitsXFaces(self,tree):
    minXPoint, maxXPoint = self.projectXBounds(tree)
    minXHits = False
    if minXPoint[1] < self.maxy and minXPoint[1] > self.miny:
      if minXPoint[2] < self.maxz and minXPoint[2] > self.minz:
        minXHits = True
    maxXHits = False
    if maxXPoint[1] < self.maxy and maxXPoint[1] > self.miny:
      if maxXPoint[2] < self.maxz and maxXPoint[2] > self.minz:
        maxXHits = True
    return minXHits, maxXHits

  def hitsYFaces(self,tree):
    minYPoint, maxYPoint = self.projectYBounds(tree)
    minYHits = False
    if minYPoint[0] < self.maxx and minYPoint[0] > self.minx:
      if minYPoint[2] < self.maxz and minYPoint[2] > self.minz:
        minYHits = True
    maxYHits = False
    if maxYPoint[0] < self.maxx and maxYPoint[0] > self.minx:
      if maxYPoint[2] < self.maxz and maxYPoint[2] > self.minz:
        maxYHits = True
    return minYHits, maxYHits

  def hitsZFaces(self,tree):
    minZPoint, maxZPoint = self.projectZBounds(tree)
    minZHits = False
    if minZPoint[1] < self.maxy and minZPoint[1] > self.miny:
      if minZPoint[0] < self.maxx and minZPoint[0] > self.minx:
        minZHits = True
    maxZHits = False
    if maxZPoint[1] < self.maxy and maxZPoint[1] > self.miny:
      if maxZPoint[0] < self.maxx and maxZPoint[0] > self.minx:
        maxZHits = True
    return minZHits, maxZHits

def plotFaceHits(tree,projector):
  hitPointsMinX = []
  notHitPointsMinX = []
  hitPointsMaxX = []
  notHitPointsMaxX = []

  hitPointsMinY = []
  notHitPointsMinY = []
  hitPointsMaxY = []
  notHitPointsMaxY = []

  hitPointsMinZ = []
  notHitPointsMinZ = []
  hitPointsMaxZ = []
  notHitPointsMaxZ = []
  for iEvent in range(tree.GetEntries()):
    if iEvent > 10000:
      break
    tree.GetEntry(iEvent)
    minXpoint, maxXpoint = projector.projectXBounds(tree)
    hitMinX,hitMaxX = projector.hitsXFaces(tree)
    minYpoint, maxYpoint = projector.projectYBounds(tree)
    hitMinY,hitMaxY = projector.hitsYFaces(tree)
    minZpoint, maxZpoint = projector.projectZBounds(tree)
    hitMinZ,hitMaxZ = projector.hitsZFaces(tree)
    if hitMinX:
      hitPointsMinX.append(minXpoint)
    else:
      notHitPointsMinX.append(minXpoint)
    if hitMaxX:
      hitPointsMaxX.append(maxXpoint)
    else:
      notHitPointsMaxX.append(maxXpoint)
    if hitMinY:
      hitPointsMinY.append(minYpoint)
    else:
      notHitPointsMinY.append(minYpoint)
    if hitMaxY:
      hitPointsMaxY.append(maxYpoint)
    else:
      notHitPointsMaxY.append(maxYpoint)
    if hitMinZ:
      hitPointsMinZ.append(minZpoint)
    else:
      notHitPointsMinZ.append(minZpoint)
    if hitMaxZ:
      hitPointsMaxZ.append(maxZpoint)
    else:
      notHitPointsMaxZ.append(maxZpoint)
  hitPointsMinX = array(hitPointsMinX)
  notHitPointsMinX = array(notHitPointsMinX)
  hitPointsMaxX = array(hitPointsMaxX)
  notHitPointsMaxX = array(notHitPointsMaxX)
  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.axvspan(minz,maxz,fc='0.8',ec='0.8')
  ax.plot(hitPointsMinX[:,2],hitPointsMinX[:,1],".b",label="Hit")
  ax.plot(notHitPointsMinX[:,2],notHitPointsMinX[:,1],'.r',label="Not Hit")
  ax.set_xlim(-200,900)
  ax.set_ylim(-200,800)
  ax.set_xlabel("z [cm]")
  ax.set_ylabel("y [cm]")
  fig.savefig("hitsMinX.png")
  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.axvspan(minz,maxz,fc='0.8',ec='0.8')
  ax.plot(hitPointsMaxX[:,2],hitPointsMaxX[:,1],".b",label="Hit")
  ax.plot(notHitPointsMaxX[:,2],notHitPointsMaxX[:,1],'.r',label="Not Hit")
  ax.set_xlim(-200,900)
  ax.set_ylim(-200,800)
  ax.set_xlabel("z [cm]")
  ax.set_ylabel("y [cm]")
  fig.savefig("hitsMaxX.png")

  hitPointsMinY = array(hitPointsMinY)
  notHitPointsMinY = array(notHitPointsMinY)
  hitPointsMaxY = array(hitPointsMaxY)
  notHitPointsMaxY = array(notHitPointsMaxY)
  fig, ax = mpl.subplots()
  ax.axhspan(minz,maxz,fc='0.8',ec='0.8')
  ax.axvspan(minx,maxx,fc='0.8',ec='0.8')
  ax.plot(hitPointsMinY[:,0],hitPointsMinY[:,2],".b",label="Hit")
  ax.plot(notHitPointsMinY[:,0],notHitPointsMinY[:,2],'.r',label="Not Hit")
  ax.set_xlim(-1000,1000)
  ax.set_ylim(-200,800)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("z [cm]")
  fig.savefig("hitsMinY.png")
  fig, ax = mpl.subplots()
  ax.axhspan(minz,maxz,fc='0.8',ec='0.8')
  ax.axvspan(minx,maxx,fc='0.8',ec='0.8')
  ax.plot(hitPointsMaxY[:,0],hitPointsMaxY[:,2],".b",label="Hit")
  ax.plot(notHitPointsMaxY[:,0],notHitPointsMaxY[:,2],'.r',label="Not Hit")
  ax.set_xlim(-1000,1000)
  ax.set_ylim(-200,800)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("z [cm]")
  fig.savefig("hitsMaxY.png")

  hitPointsMinZ = array(hitPointsMinZ)
  notHitPointsMinZ = array(notHitPointsMinZ)
  hitPointsMaxZ = array(hitPointsMaxZ)
  notHitPointsMaxZ = array(notHitPointsMaxZ)
  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.axvspan(minx,maxx,fc='0.8',ec='0.8')
  ax.plot(hitPointsMinZ[:,0],hitPointsMinZ[:,1],".b",label="Hit")
  ax.plot(notHitPointsMinZ[:,0],notHitPointsMinZ[:,1],'.r',label="Not Hit")
  ax.set_xlim(-1000,1000)
  ax.set_ylim(-200,800)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("y [cm]")
  fig.savefig("hitsMinZ.png")
  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.axvspan(minx,maxx,fc='0.8',ec='0.8')
  ax.plot(hitPointsMaxZ[:,0],hitPointsMaxZ[:,1],".b",label="Hit")
  ax.plot(notHitPointsMaxZ[:,0],notHitPointsMaxZ[:,1],'.r',label="Not Hit")
  ax.set_xlim(-1000,1000)
  ax.set_ylim(-200,800)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("y [cm]")
  fig.savefig("hitsMaxZ.png")

def plotCoordVTheta(tree,projector):
  hitPointsMinZ = []
  hitPointsMaxZ = []
  thetazMinZ = []
  thetazMaxZ = []
  for iEvent in range(tree.GetEntries()):
    if iEvent > 10000:
      break
    tree.GetEntry(iEvent)
    minZpoint, maxZpoint = projector.projectZBounds(tree)
    hitMinZ,hitMaxZ = projector.hitsZFaces(tree)
    if hitMinZ:
      hitPointsMaxZ.append(maxZpoint)
      thetazMaxZ.append(tree.thetaz)
    if hitMaxZ:
      hitPointsMinZ.append(minZpoint)
      thetazMinZ.append(tree.thetaz)
  hitPointsMinZ = array(hitPointsMinZ)
  hitPointsMaxZ = array(hitPointsMaxZ)
  thetazMinZ = array(thetazMinZ)*180/pi
  thetazMaxZ = array(thetazMaxZ)*180/pi

  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.plot(thetazMinZ,hitPointsMinZ[:,1],".b")
  #ax.set_xlim(-200,900)
  ax.set_ylim(-200,800)
  ax.set_xlabel(r"$\theta$ [deg]")
  ax.set_ylabel("y [cm]")
  fig.savefig("minZyVthetaz.png")

  fig, ax = mpl.subplots()
  ax.axhspan(miny,maxy,fc='0.8',ec='0.8')
  ax.plot(thetazMaxZ,hitPointsMaxZ[:,1],".b")
  #ax.set_xlim(-200,900)
  ax.set_ylim(-200,800)
  ax.set_xlabel(r"$\theta$ [deg]")
  ax.set_ylabel("y [cm]")
  fig.savefig("maxZyVthetaz.png")

if __name__ == "__main__":
  c = root.TCanvas()

  minx = -360
  maxx = 360
  miny = 0.
  maxy = 607.5
  minz = -0.5
  maxz = 695

  projector = Projector(minx,maxx,miny,maxy,minz,maxz)
  projectorLeft = Projector(minx,0.,miny,maxy,minz,maxz)
  projectorRight = Projector(0.,maxx,miny,maxy,minz,maxz)

  infilename = "cosmics.root" 
  sampleRate = 8964.12
  sampleArea = sampleRate/0.014/100**2
  detectorArea = (maxx-minx)*(maxz-minz)/100**2
  detectorRate = sampleRate*detectorArea/sampleArea
  print "Sample: "+infilename
  print "Sample Rate: {} Hz".format(sampleRate)
  print "Sample Area: {} m^2".format(sampleArea)
  print "Detector Rate: {} Hz".format(detectorRate)
  print "Detector Area: {} m^2".format(detectorArea)
  

  infile = root.TFile(infilename)
  tree = infile.Get("tree")

  #plotFaceHits(tree,projector)
  #plotCoordVTheta(tree,projector)
  
  nEventsTotal = 0
  nEventsHitFrontBackZ = 0

  
  dxHitFrontBackHist = Hist(20,0,720)
  dxHitFrontBackLeftHist = Hist(10,0,360)
  dxHitFrontBackRightHist = Hist(10,0,360)

  for iEvent in range(tree.GetEntries()):
    if iEvent > 10000:
      #break
      pass
    tree.GetEntry(iEvent)
    nEventsTotal += 1
    minpoint, maxpoint = projector.projectZBounds(tree)
    hitMin,hitMax = projector.hitsZFaces(tree)
    if hitMin and hitMax:
      nEventsHitFrontBackZ += 1
      dxHitFrontBackHist.Fill(abs(maxpoint[0]-minpoint[0]))
    hitMinLeft, hitMaxLeft = projectorLeft.hitsZFaces(tree)
    if hitMinLeft and hitMaxLeft:
      dxHitFrontBackLeftHist.Fill(abs(maxpoint[0]-minpoint[0]))
    hitMinRight, hitMaxRight = projectorRight.hitsZFaces(tree)
    if hitMinRight and hitMaxRight:
      dxHitFrontBackRightHist.Fill(abs(maxpoint[0]-minpoint[0]))
    
    #print("Event: {:4} x y z: {:6.1f} {:6.1f} {:6.1f} px py pz E: {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(iEvent,tree.x,tree.y,tree.z,tree.px,tree.py,tree.pz,tree.E))
    #print(" {}    {}".format(minpoint,maxpoint))
    #print(" {}    {}".format(hitMin,hitMax))

  eventWeight = sampleRate/float(nEventsTotal)
  print("eventWeight: {} Hz".format(eventWeight))
  print("Rate going through front and back: {} Hz".format(nEventsHitFrontBackZ*eventWeight))

  dxHitFrontBackHist.Scale(eventWeight)
  setHistTitles(dxHitFrontBackHist,"|x_{front}-x_{rear}| [cm]","Rate/bin [Hz]")
  #normToBinWidth(dxHitFrontBackHist)
  #setHistTitles(dxHitFrontBackHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dxHitFrontBackHist.UseCurrentStyle()
  dxHitFrontBackHist.Draw("hist")
  drawStandardCaptions(c,"",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBack.pdf")

  dxHitFrontBackIntegralHist = getIntegralHist(dxHitFrontBackHist)
  dxHitFrontBackIntegralHist.GetYaxis().SetTitle("Rate for X-axis value #geq |x_{front}-x_{rear}| [Hz]")
  #setHistTitles(dxHitFrontBackHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dxHitFrontBackIntegralHist.Draw("hist")
  drawStandardCaptions(c,"",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackIntegral.pdf")

  #######

  dxHitFrontBackLeftHist.Scale(eventWeight)
  setHistTitles(dxHitFrontBackLeftHist,"|x_{front}-x_{rear}| [cm]","Rate/bin [Hz]")
  dxHitFrontBackLeftHist.UseCurrentStyle()
  dxHitFrontBackLeftHist.Draw("hist")
  drawStandardCaptions(c,"Left TPC",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackLeft.pdf")

  dxHitFrontBackRightHist.Scale(eventWeight)
  setHistTitles(dxHitFrontBackRightHist,"|x_{front}-x_{rear}| [cm]","Rate/bin [Hz]")
  dxHitFrontBackRightHist.UseCurrentStyle()
  dxHitFrontBackRightHist.Draw("hist")
  drawStandardCaptions(c,"Right TPC",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackRight.pdf")
