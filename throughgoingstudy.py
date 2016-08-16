#!/usr/bin/env python

import re
import math

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)

from scipy import *

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
    return array([xplane,y,z])
    
  def projectY(self,tree,yplane):
    # first x coord
    mx = tree.px/tree.py
    bx = tree.x - mx*tree.y
    x = mx*yplane + bx

    # then z coord
    mz = tree.pz/tree.py
    bz = tree.z - mz*tree.y
    z = mz*yplane + bz
    return array([x,yplane,z])
    
  def projectZ(self,tree,zplane):
    # first x coord
    mx = tree.px/tree.pz
    bx = tree.x - mx*tree.z
    x = mx*zplane + bx

    # then y coord
    my = tree.py/tree.pz
    by = tree.y - my*tree.z
    y = my*zplane + by
    return array([x,y,zplane])

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

if __name__ == "__main__":

  minx = -360
  maxx = 360
  miny = 0.
  maxy = 607.5
  minz = -0.5
  maxz = 695

  projector = Projector(minx,maxx,miny,maxy,minz,maxz)

  infilename = "cosmics.root" 
  infile = root.TFile(infilename)
  tree = infile.Get("tree")
  for iEvent in range(tree.GetEntries()):
    if iEvent > 10:
      break
    tree.GetEntry(iEvent)
    minpoint, maxpoint = projector.projectZBounds(tree)
    hitMin,hitMax = projector.hitsZFaces(tree)
    
    print("Event: {:4} x y z: {:6.1f} {:6.1f} {:6.1f} px py pz E: {:7.2f} {:7.2f} {:7.2f} {:7.2f}".format(iEvent,tree.x,tree.y,tree.z,tree.px,tree.py,tree.pz,tree.E))
    print(" {}    {}".format(minpoint,maxpoint))
    print(" {}    {}".format(hitMin,hitMax))
