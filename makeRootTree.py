#!/usr/bin/env python

import math
from rootpy.plotting import Hist, Hist2D, Canvas, Graph
from rootpy.tree import Tree
from rootpy.io import root_open
from rootpy import ROOT as root
import helpers
root.gROOT.SetBatch(True)

canvas = Canvas()

def makeRootTree(infilename,outfilename,maxEvents,diagnostics):

  with root_open(outfilename, "recreate") as outfile:
    tree = Tree("tree")
    tree.create_branches(
        {
         # HEPEVT info
         'status': 'I',
         'pdgid': 'I',
         'x': 'F',
         'y': 'F',
         'z': 'F',
         't': 'F',
         'px': 'F',
         'py': 'F',
         'pz': 'F',
         'E': 'F',
         'm': 'F',
         # Additional info
         'theta': 'F',
         'phi': 'F',
         'thetaz': 'F',
         'phiz': 'F',
        })
    
    with open(infilename) as infile:
      for iline, line in enumerate(infile):
        if maxEvents >0:
          if iline*2 >= maxEvents:
            break
        """
        <Event number> <Number of entries in this event> # lines after this line are particle entries
        <Status> <PDG ID> <1st Mother> <2nd Mother> <1st Daughter> <2nd Daughter> <Px> <Py> <Px> <E> <Mass> <x> <y> <z> <t>
        """
        line = line.split(" ")
        if len(line) < 15:
          continue
        tree.status = line[0]
        tree.pdgid = line[1]
        tree.px = line[6]
        tree.py = line[7]
        tree.pz = line[8]
        tree.E = line[9]
        tree.m = line[10]
        tree.x = line[11]
        tree.y = line[12]
        tree.z = line[13]
        tree.t = line[14]
        tree.thetaz = math.atan2((tree.px**2+tree.pz**2)**0.5,-tree.py)
        tree.phiz = math.atan2(tree.pz,tree.px)
        tree.theta = math.atan2((tree.px**2+tree.py**2)**0.5,tree.pz)
        tree.phi = math.atan2(tree.py,tree.px)
        tree.Fill()
    tree.write()
  
    if diagnostics:
      h = tree.Draw("px")
      canvas.SaveAs("px.png")
      h = tree.Draw("py")
      canvas.SaveAs("py.png")
      h = tree.Draw("pz")
      canvas.SaveAs("pz.png")
      h = tree.Draw("E")
      canvas.SaveAs("E.png")
      h = tree.Draw("m")
      canvas.SaveAs("m.png")
      h = tree.Draw("x")
      canvas.SaveAs("x.png")
      h = tree.Draw("y")
      canvas.SaveAs("y.png")
      h = tree.Draw("z")
      canvas.SaveAs("z.png")
      h = tree.Draw("t")
      canvas.SaveAs("t.png")
      h = tree.Draw("theta")
      canvas.SaveAs("theta.png")
      h = tree.Draw("phi")
      canvas.SaveAs("phi.png")
      h = tree.Draw("thetaz")
      canvas.SaveAs("thetaz.png")
      h = tree.Draw("phiz")
      canvas.SaveAs("phiz.png")

if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser(description='Makes a root tree from a HEPEVT file. Assumes one particle per event.')
  parser.add_argument('infilename',
                      help='Input file name')
  parser.add_argument('outfilename',
                      help='Output file name')
  parser.add_argument('--nEvents','-n', type=int, default=-1,
                      help='Number of events to process')
  parser.add_argument('--diagnostics','-d', action="store_true",
                      help='Make diagnostic plots')
  
  args = parser.parse_args()
  makeRootTree(args.infilename,args.outfilename,args.nEvents,args.diagnostics)
