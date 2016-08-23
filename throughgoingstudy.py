#!/usr/bin/env python

import re
import math
import sys

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

  def getExtremaZPoints(self,tree):
    minXHits, maxXHits = self.hitsXFaces(tree)
    minYHits, maxYHits = self.hitsYFaces(tree)
    minZHits, maxZHits = self.hitsZFaces(tree)
    minXPoint, maxXPoint = self.projectXBounds(tree)
    minYPoint, maxYPoint = self.projectYBounds(tree)
    minZPoint, maxZPoint = self.projectZBounds(tree)
    points = [minZPoint,maxZPoint,minYPoint,maxYPoint,minXPoint,maxXPoint]
    hits = [minZHits,maxZHits,minYHits,maxYHits,minXHits,maxXHits]
    minPoint = [0.,0.,1.e15]
    maxPoint = [0.,0.,-1.e15]
    enteredDetector = False
    for point,hit in zip(points,hits):
      if hit:
        enteredDetector = True
        if point[2] < minPoint[2]:
          minPoint = point
        if point[2] > maxPoint[2]:
          maxPoint = point
    if not enteredDetector:
      return None, None
    assert(minPoint[2]< 1.e15)
    assert(maxPoint[2]> -1.e15)
    return minPoint, maxPoint

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

def estimatePerSrForVerticalAndEgt1GeV(tree,nmax=10000000000):
  canvas = root.TCanvas("c2")

  func = root.TF1("fitfunc","[0]*cos(x)*cos(x)",0,pi)

  cuts = "E>1."
  theta = Hist(100,0,0.5*pi)
  theta.Sumw2()
  theta.UseCurrentStyle()
  tree.Draw("thetaz >> {}".format(theta.GetName()),cuts,"hist",nmax)
  setHistTitles(theta,"#theta","Events/Sr")
  drawStandardCaptions(canvas,"E_{#mu} > 1 GeV")
  xAxis = theta.GetXaxis()
  for iBinX in range(1,xAxis.GetNbins()+1):
    xLow = xAxis.GetBinLowEdge(iBinX)
    xHigh = xAxis.GetBinUpEdge(iBinX)
    solidAngle = 2*pi*(-cos(xHigh)+cos(xLow))
    binContent = theta.GetBinContent(iBinX)
    theta.SetBinContent(iBinX,binContent/solidAngle)
  theta.Draw("E")
  fitResult = theta.Fit(func,"WLMSQ",'',0.05,0.3)
  canvas.SaveAs("normalizationFit.png")

  normalization = fitResult.Parameter(0)
  normalizationUnc = fitResult.ParError(0)
  return normalization, normalizationUnc

if __name__ == "__main__":
  c = root.TCanvas()

  NMAX = 10000000
  #NMAX = 2000
  minx = -360
  maxx = 360
  miny = 0.
  maxy = 607.5
  minz = -0.5
  maxz = 695

  distanceToZFaceCut = 35.
  distanceToXFaceCut = 18.

  projector = Projector(minx,maxx,miny,maxy,minz,maxz)
  projectorLeft = Projector(minx,0.,miny,maxy,minz,maxz)
  projectorRight = Projector(0.,maxx,miny,maxy,minz,maxz)

  infilename = "cosmicsJustInProtoDUNE.root" 
  sampleRate = 8887.949
  #sampleArea = sampleRate/0.014/100**2
  sampleArea = (maxx-minx)*(maxz-minz)/100**2
  detectorArea = (maxx-minx)*(maxz-minz)/100**2
  detectorRate = sampleRate*detectorArea/sampleArea
  print "Sample: "+infilename
  print "Sample Rate: {} Hz".format(sampleRate)
  print "Sample Area: {} m^2".format(sampleArea)
  print "Detector Rate: {} Hz".format(detectorRate)
  print "Detector Area: {} m^2".format(detectorArea)

  infile = root.TFile(infilename)
  tree = infile.Get("tree")
  nEntries = min(NMAX,tree.GetEntries())
  print "N events in sample: {}".format(nEntries)

  nEvtsVertPerSr, nEvtsVertPerSrUnc = estimatePerSrForVerticalAndEgt1GeV(tree,NMAX)
  print "N events per Sr vertical and E>1GeV: {} +/- {}".format(nEvtsVertPerSr,nEvtsVertPerSrUnc)
  print "N events per Sr per m^2 vertical and E>1GeV: {}".format(nEvtsVertPerSr/sampleArea)
  print "Vertical rate per Sr per m^2 and E>1GeV should be: {} Hz".format(70.)
  print "Vertical rate per Sr and E>1GeV should be: {} Hz".format(70.*sampleArea)
  print "Sample corresponds to {} s".format(nEvtsVertPerSr/sampleArea/70.)
  eventWeight = 70.*sampleArea/float(nEvtsVertPerSr)
  print "eventWeight: {} Hz".format(eventWeight)
  print "Real Sample Rate: {} Hz".format(eventWeight*nEntries)

  #plotFaceHits(tree,projector)
  #plotCoordVTheta(tree,projector)
  
  nEventsTotal = 0
  nEventsHitFrontBackZ = 0
  nEventsEntered = 0
  nEventsEnteredLeft = 0
  nEventsEnteredRight = 0
  nEventsCuts = 0
  nEventsCutsLeft = 0
  nEventsCutsRight = 0
  
  dxHitFrontBackHist = Hist(20,0,720)
  dxHitFrontBackLeftHist = Hist(10,0,360)
  dxHitFrontBackRightHist = Hist(10,0,360)
  dxHitCutsHist = Hist(20,0,720)
  dxHitCutsLeftHist = Hist(10,0,360)
  dxHitCutsRightHist = Hist(10,0,360)
  dzHitCutsHist = Hist(20,maxz-minz-distanceToZFaceCut,maxz-minz)
  dzHitCutsLeftHist = Hist(20,maxz-minz-distanceToZFaceCut,maxz-minz)
  dzHitCutsRightHist = Hist(20,maxz-minz-distanceToZFaceCut,maxz-minz)
  thetaFrontBackHist = Hist(30,0,90)
  thetaCutsHist = Hist(30,0,90)

  #dxdzHitLeftHist = Hist2D(90,0.,maxx,50,0.,maxz-minz)
  #dxdzHitRightHist = Hist2D(90,0.,maxx,50,0.,maxz-minz)
  dxdzHitLeftHist = Hist2D(90,0.,maxx,[0,200,400,600,650,maxz-minz-25,maxz-minz-10])
  dxdzHitRightHist = Hist2D(90,0.,maxx,[0,200,400,600,650,maxz-minz-25,maxz-minz-10])

  for iEvent in range(nEntries):
    tree.GetEntry(iEvent)
    nEventsTotal += 1
    #minpointX, maxpointX = projector.projectXBounds(tree)
    #minpointY, maxpointY = projector.projectYBounds(tree)
    minpointZ, maxpointZ = projector.projectZBounds(tree)
    #minpointLeftX, maxpointLeftX = projectorLeft.projectXBounds(tree)
    #minpointLeftY, maxpointLeftY = projectorLeft.projectYBounds(tree)
    minpointLeftZ, maxpointLeftZ = projectorLeft.projectZBounds(tree)
    #minpointRightX, maxpointRightX = projectorRight.projectXBounds(tree)
    #minpointRightY, maxpointRightY = projectorRight.projectYBounds(tree)
    minpointRightZ, maxpointRightZ = projectorRight.projectZBounds(tree)
    hitMinX, hitMaxX = projector.hitsXFaces(tree)
    #hitMinY, hitMaxY = projector.hitsYFaces(tree)
    hitMinZ, hitMaxZ = projector.hitsZFaces(tree)
    hitMinXLeft, hitMaxXLeft = projectorLeft.hitsXFaces(tree)
    #hitMinYLeft, hitMaxYLeft = projectorLeft.hitsYFaces(tree)
    hitMinZLeft, hitMaxZLeft = projectorLeft.hitsZFaces(tree)
    hitMinXRight, hitMaxXRight = projectorRight.hitsXFaces(tree)
    #hitMinYRight, hitMaxYRight = projectorRight.hitsYFaces(tree)
    hitMinZRight, hitMaxZRight = projectorRight.hitsZFaces(tree)
    ## Just hitting front and back of TPC(s)
    if hitMinZ and hitMaxZ:
      nEventsHitFrontBackZ += 1
      dxHitFrontBackHist.Fill(abs(maxpointZ[0]-minpointZ[0]),eventWeight)
      thetaFrontBackHist.Fill(tree.thetaz*180./pi,eventWeight)
    if hitMinZLeft and hitMaxZLeft:
      dxHitFrontBackLeftHist.Fill(abs(maxpointZ[0]-minpointZ[0]),eventWeight)
    if hitMinZRight and hitMaxZRight:
      dxHitFrontBackRightHist.Fill(abs(maxpointZ[0]-minpointZ[0]),eventWeight)
    ## Getting entry/exit points ordered by z
    minpoint, maxpoint = projector.getExtremaZPoints(tree)
    minpointLeft, maxpointLeft = projectorLeft.getExtremaZPoints(tree)
    minpointRight, maxpointRight = projectorRight.getExtremaZPoints(tree)
    if minpoint:
      nEventsEntered += 1
      if minpoint[2] < minz + distanceToZFaceCut and maxpoint[2] > maxz - distanceToZFaceCut:
        #if minpoint[0] > maxpoint[0]  and minpoint[0] > maxx - distanceToXFaceCut and maxpoint[0] < minx + distanceToXFaceCut:
        #    dxHitCutsHist.Fill(abs(maxpoint[0]-minpoint[0]),eventWeight)
        #    dzHitCutsHist.Fill(abs(maxpoint[2]-minpoint[2]),eventWeight)
        #    nEventsCuts += 1
        #elif maxpoint[0] > minpoint[0] and maxpoint[0] > maxx - distanceToXFaceCut and minpoint[0] < minx + distanceToXFaceCut:
            dxHitCutsHist.Fill(abs(maxpoint[0]-minpoint[0]),eventWeight)
            dzHitCutsHist.Fill(abs(maxpoint[2]-minpoint[2]),eventWeight)
            thetaCutsHist.Fill(tree.thetaz*180/pi,eventWeight)
            nEventsCuts += 1
    if minpointLeft:
      nEventsEnteredLeft += 1
      dxdzHitLeftHist.Fill(abs(maxpointLeft[0]-minpointLeft[0]),abs(maxpointLeft[2]-minpointLeft[2]),eventWeight)
      if minpointLeft[2] < minz + distanceToZFaceCut and maxpointLeft[2] > maxz - distanceToZFaceCut:
        #if minpointLeft[0] > maxpointLeft[0] and minpointLeft[0] > 0. - distanceToXFaceCut and maxpointLeft[0] < minx + distanceToXFaceCut:
        #  nEventsCutsLeft += 1
        #  dxHitCutsLeftHist.Fill(abs(maxpointLeft[0]-minpointLeft[0]),eventWeight)
        #  dzHitCutsLeftHist.Fill(abs(maxpointLeft[2]-minpointLeft[2]),eventWeight)
        #elif maxpointLeft[0] > minpointLeft[0] and maxpointLeft[0] > 0. - distanceToXFaceCut and minpointLeft[0] < minx + distanceToXFaceCut:
          dxHitCutsLeftHist.Fill(abs(maxpointLeft[0]-minpointLeft[0]),eventWeight)
          dzHitCutsLeftHist.Fill(abs(maxpointLeft[2]-minpointLeft[2]),eventWeight)
          nEventsCutsLeft += 1
    if minpointRight:
      nEventsEnteredRight += 1
      dxdzHitRightHist.Fill(abs(maxpointRight[0]-minpointRight[0]),abs(maxpointRight[2]-minpointRight[2]),eventWeight)
      if minpointRight[2] < minz + distanceToZFaceCut and maxpointRight[2] > maxz - distanceToZFaceCut:
        #if minpointRight[0] > maxpointRight[0] and minpointRight[0] > maxx - distanceToXFaceCut and maxpointRight[0] < 0. + distanceToXFaceCut:
        #  nEventsCutsRight += 1
        #  dxHitCutsRightHist.Fill(abs(maxpointRight[0]-minpointRight[0]),eventWeight)
        #  dzHitCutsRightHist.Fill(abs(maxpointRight[2]-minpointRight[2]),eventWeight)
        #elif maxpointRight[0] > minpointRight[0] and maxpointRight[0] > maxx - distanceToXFaceCut and minpointRight[0] < 0. + distanceToXFaceCut:
          nEventsCutsRight += 1
          dxHitCutsRightHist.Fill(abs(maxpointRight[0]-minpointRight[0]),eventWeight)
          dzHitCutsRightHist.Fill(abs(maxpointRight[2]-minpointRight[2]),eventWeight)
    
    #print("Event: {:4} x y z: {:6.1f} {:6.1f} {:6.1f} px py pz E: {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(iEvent,tree.x,tree.y,tree.z,tree.px,tree.py,tree.pz,tree.E))
    #print(" {}    {}".format(minpointZ,maxpointZ))
    #print(" {}    {}".format(hitMin,hitMax))

  print("Rate going through front and back: {} Hz".format(nEventsHitFrontBackZ*eventWeight))
  print("Rate going through detector: {} Hz".format(nEventsEntered*eventWeight))
  print("Rate going through left detector: {} Hz".format(nEventsEnteredLeft*eventWeight))
  print("Rate going through right detector: {} Hz".format(nEventsEnteredRight*eventWeight))
  print("Rate passing cuts, whole detector: {} Hz".format(nEventsCuts*eventWeight))
  print("Rate passing cuts, left detector: {} Hz".format(nEventsCutsLeft*eventWeight))
  print("Rate passing cuts, right detector: {} Hz".format(nEventsCutsRight*eventWeight))

  setHistTitles(dxHitFrontBackHist,"|#Delta x| [cm]","Rate/bin [Hz]")
  #normToBinWidth(dxHitFrontBackHist)
  #setHistTitles(dxHitFrontBackHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dxHitFrontBackHist.UseCurrentStyle()
  dxHitFrontBackHist.Draw("hist")
  drawStandardCaptions(c,"",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBack.png")
  c.SaveAs("dxHitFrontBack.pdf")

  dxHitFrontBackIntegralHist = getIntegralHist(dxHitFrontBackHist)
  dxHitFrontBackIntegralHist.GetYaxis().SetTitle("Rate for X-axis value #geq |#Delta x| [Hz]")
  #setHistTitles(dxHitFrontBackHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dxHitFrontBackIntegralHist.Draw("hist")
  drawStandardCaptions(c,"Includes both left and right TPCs",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackIntegral.png")
  c.SaveAs("dxHitFrontBackIntegral.pdf")

  #######

  setHistTitles(dxHitFrontBackLeftHist,"|#Delta x| [cm]","Rate/bin [Hz]")
  dxHitFrontBackLeftHist.UseCurrentStyle()
  dxHitFrontBackLeftHist.Draw("hist")
  drawStandardCaptions(c,"Left TPC",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackLeft.png")
  c.SaveAs("dxHitFrontBackLeft.pdf")

  setHistTitles(dxHitFrontBackRightHist,"|#Delta x| [cm]","Rate/bin [Hz]")
  dxHitFrontBackRightHist.UseCurrentStyle()
  dxHitFrontBackRightHist.Draw("hist")
  drawStandardCaptions(c,"Right TPC",captionright1="For events passing through front",captionright2="& rear z-faces of detector")
  c.SaveAs("dxHitFrontBackRight.png")
  c.SaveAs("dxHitFrontBackRight.pdf")

  #######
  ## Cuts

  setHistTitles(dxHitCutsHist,"|#Delta x| [cm]","Rate/bin [Hz]")
  #normToBinWidth(dxHitCutsHist)
  #setHistTitles(dxHitCutsHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dxHitCutsHist.UseCurrentStyle()
  dxHitCutsHist.Draw("hist")
  drawStandardCaptions(c,"Across left and right TPCs",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dxHitCuts.png")
  c.SaveAs("dxHitCuts.pdf")

  setHistTitles(dxHitCutsLeftHist,"|#Delta x| [cm]","Rate/bin [Hz]")
  dxHitCutsLeftHist.UseCurrentStyle()
  dxHitCutsLeftHist.Draw("hist")
  drawStandardCaptions(c,"Left TPC",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dxHitCutsLeft.png")
  c.SaveAs("dxHitCutsLeft.pdf")

  setHistTitles(dxHitCutsRightHist,"|#Delta z| [cm]","Rate/bin [Hz]")
  dxHitCutsRightHist.UseCurrentStyle()
  dxHitCutsRightHist.Draw("hist")
  drawStandardCaptions(c,"Right TPC",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dxHitCutsRight.png")
  c.SaveAs("dxHitCutsRight.pdf")

  setHistTitles(dzHitCutsHist,"|#Delta z| [cm]","Rate/bin [Hz]")
  #normToBinWidth(dzHitCutsHist)
  #setHistTitles(dzHitCutsHist,"|x_{front}-x_{rear}| [cm]","Rate/|x_{front}-x_{rear}| [Hz/cm]")
  dzHitCutsHist.UseCurrentStyle()
  dzHitCutsHist.Draw("hist")
  drawStandardCaptions(c,"Across left and right TPCs",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dzHitCuts.png")
  c.SaveAs("dzHitCuts.pdf")

  setHistTitles(dzHitCutsLeftHist,"|#Delta z| [cm]","Rate/bin [Hz]")
  dzHitCutsLeftHist.UseCurrentStyle()
  dzHitCutsLeftHist.Draw("hist")
  drawStandardCaptions(c,"Left TPC",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dzHitCutsLeft.png")
  c.SaveAs("dzHitCutsLeft.pdf")

  setHistTitles(dzHitCutsRightHist,"|#Delta z| [cm]","Rate/bin [Hz]")
  dzHitCutsRightHist.UseCurrentStyle()
  dzHitCutsRightHist.Draw("hist")
  drawStandardCaptions(c,"Right TPC",captionright1="Track within {} cm of z-faces".format(distanceToXFaceCut))
  c.SaveAs("dzHitCutsRight.png")
  c.SaveAs("dzHitCutsRight.pdf")

  thetaFrontBackHist.UseCurrentStyle()
  setHistTitles(thetaFrontBackHist,"#theta with respect to zenith [deg]","Events/bin")
  thetaFrontBackHist.Draw("hist")
  c.SaveAs("thetaHists.png")
  c.SaveAs("thetaHists.pdf")
  print("thetaFrontBackHist integral: %s" % thetaFrontBackHist.Integral())

  ############
  ### 2D plots of Delta z v Delta x
  setupCOLZFrame(c)
  #c.SetLogz()
  setHistTitles(dxdzHitLeftHist,"|#Delta x| [cm]", "|#Delta z| [cm]","Rate/bin [Hz]")
  #setHistRange(dxdzHitLeftHist,300,maxx,400,maxz-minz)
  dxdzHitLeftHist.Draw("colz")
  drawStandardCaptions(c,"Left TPC")
  c.SaveAs("dxdzHitLeft.png")
  c.SaveAs("dxdzHitLeft.pdf")

  dxdzHitLeftIntegralHist = getIntegralHist(dxdzHitLeftHist)
  dxdzHitLeftIntegralHist.GetZaxis().SetTitle("Rate for X #geq |#Delta x| and Y #geq |#Delta z| [Hz]")
  #setHistRange(dxdzHitLeftIntegralHist,300,maxx,400,maxz-minz)
  dxdzHitLeftIntegralHist.Draw("colz")
  drawStandardCaptions(c,"Left TPC")
  c.SaveAs("dxdzHitLeftIntegral.png")
  c.SaveAs("dxdzHitLeftIntegral.pdf")

  setupCOLZFrame(c,True)
  c.SetLogy(True)
  hists = []
  labels = []
  for i, xBin in enumerate([4,5,6,7]):
    deltaZMin = dxdzHitLeftIntegralHist.GetYaxis().GetBinLowEdge(xBin)
    hist = getYBinHist(dxdzHitLeftIntegralHist,xBin)
    hist.UseCurrentStyle()
    hist.SetLineColor(i+1)
    hists.append(hist)
    labels.append(r"|#Delta z| #geq {:.0f} cm".format(deltaZMin))
  axisHist = makeStdAxisHist(hists,logy=True,freeTopSpace=0.35,xlim=[200,maxx])
  axisHist.Draw()
  setHistTitles(axisHist,"|#Delta x| [cm]","Rate for X #geq |#Delta x| [Hz]")
  for hist in hists:
    hist.Draw("histsame")
  leg = drawNormalLegend(hists,labels)
  c.SaveAs("dxdzHitLeftIntegralScan.png")
  c.SaveAs("dxdzHitLeftIntegralScan.pdf")
  c.SetLogy(False)
  setupCOLZFrame(c,True)

  setHistTitles(dxdzHitRightHist,"|#Delta x| [cm]", "|#Delta z| [cm]","Rate/bin [Hz]")
  #setHistRange(dxdzHitRightHist,300,maxx,400,maxz-minz)
  drawStandardCaptions(c,"Right TPC")
  dxdzHitRightHist.Draw("colz")
  c.SaveAs("dxdzHitRight.png")
  c.SaveAs("dxdzHitRight.pdf")

  dxdzHitRightIntegralHist = getIntegralHist(dxdzHitRightHist)
  dxdzHitRightIntegralHist.GetZaxis().SetTitle("Rate for X #geq |#Delta x| and Y #geq |#Delta z| [Hz]")
  #setHistRange(dxdzHitRightIntegralHist,300,maxx,400,maxz-minz)
  dxdzHitRightIntegralHist.Draw("colz")
  drawStandardCaptions(c,"Right TPC")
  c.SaveAs("dxdzHitRightIntegral.png")
  c.SaveAs("dxdzHitRightIntegral.pdf")



