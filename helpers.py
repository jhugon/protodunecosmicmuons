
import ROOT as root
from ROOT import gStyle as gStyle
#from ROOT import RooRealVar, RooGaussian, RooArgList, RooDataHist
import re
import csv
import glob
from math import exp
from math import sqrt
from math import log10
import math
import array
import os
import sys
import time
import datetime
import random
#import matplotlib.pyplot as mpl

def getOrdinalStr(inInt):
  result = str(inInt)
  if result[-1] == "1":
    result += "st"
  elif result[-1] == "2":
    result += "nd"
  elif result[-1] == "3":
    result += "rd"
  else:
    result += "th"
  return result

# calculate FWHM
def calcFWHM(pdf,obs,min,max,step):

  var = pdf.getObservables(root.RooArgSet(obs)).first();

  ymaxVal = float(0)
  xmaxVal = float(0)

  # find the maximum value
  for x in drange(min,max,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    if (pdfVal > ymaxVal):
       xmaxVal = x
       ymaxVal = pdfVal
       
    #print "x=%s, pdfVal=%s" % (x,pdfVal)

  #print "xMax=%s, ymaxVal=%s\n\n\n\n" % (xmaxVal,ymaxVal)


  # find lower boundary with y=max/2
  xLow = float(0)
  for x in drange(min,max,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    #print "x=%s, pdfVal=%s, ymaxVal/2.=%s" % (x,pdfVal, ymaxVal/2.)
    if (pdfVal > ymaxVal/2. and xLow==0):
       xLow = x

  #print "xLow=%s" % xLow
  

  # find higher boundary with y=max/2
  xHigh = float(0)
  for x in revdrange(max,min,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    if (pdfVal > ymaxVal/2. and xHigh==0):
       xHigh = x

  #print "xHigh=%s" % xHigh
  
  return (xHigh-xLow)

def doubleGauss(x,par):
  meanG1  = par[0]
  widthG1 = par[1]
  meanG2  = par[2]
  widthG2 = par[3]
  mixGG   = par[4]
  scale   = par[5]
  
  #if (par[1] != 0.0):
  
  arg1 = (x[0]-meanG1)/widthG1
  arg2 = (x[0]-meanG2)/widthG2
  
  gauss1 = exp(-0.5*arg1*arg1)
  gauss2 = exp(-0.5*arg2*arg2)
  dgauss = (1-mixGG)*gauss1 + mixGG*gauss2 
  
  return scale*dgauss
  #return meanG1 + widthG1*x[0]
  
def getXBinHist(inHist, xBin):
  """
  Makes a TH1 hisogram from a TH2
  A vertical slice of a 2D histo
  """
  outHist = inHist.ProjectionY()
  outHist.Reset()
  outHist.SetName(inHist.GetName()+"XSliceBin"+str(xBin))
  outHist.Sumw2()
  nBins = outHist.GetXaxis().GetNbins()
  for i in range(0,nBins+2):
    outHist.SetBinContent(i,inHist.GetBinContent(xBin,i))
    outHist.SetBinError(i,inHist.GetBinError(xBin,i))
  return outHist

def getYBinHist(inHist, yBin):
  """
  Makes a TH1 hisogram from a TH2
  A horizontal slice of a 2D histo
  """
  outHist = inHist.ProjectionX()
  outHist.Reset()
  outHist.SetName(inHist.GetName()+"YSliceBin"+str(yBin))
  outHist.Sumw2()
  nBins = outHist.GetXaxis().GetNbins()
  for i in range(0,nBins+2):
    outHist.SetBinContent(i,inHist.GetBinContent(i,yBin))
    outHist.SetBinError(i,inHist.GetBinError(i,yBin))
  return outHist

def divideYValByXVal(hist):
    nBinsX = hist.GetXaxis().GetNbins()
    for iBinX in range(1,nBinsX+1):
	binVal = hist.GetBinContent(iBinX)
	binErrVal = hist.GetBinError(iBinX)
	xVal = hist.GetXaxis().GetBinCenter(iBinX)
	hist.SetBinContent(iBinX,binVal/xVal)
	hist.SetBinError(iBinX,binErrVal/xVal)

def setNormalColorTable():
  rArray = array.array('d',[0.0,1.0,1.0])
  gArray = array.array('d',[1.0,1.0,0.0])
  bArray = array.array('d',[0.0,0.0,0.0])
  stopArray = array.array('d',[0.,0.5,1.])
  nTabColors = 500
  root.TColor.CreateGradientColorTable(len(stopArray),
            stopArray,rArray,gArray,bArray,nTabColors
         )
def setInvertColorTable():
  rArray = array.array('d',[1.0,1.0,0.0])
  gArray = array.array('d',[0.0,1.0,1.0])
  bArray = array.array('d',[0.0,0.0,0.0])
  stopArray = array.array('d',[0.,0.5,1.])
  nTabColors = 500
  root.TColor.CreateGradientColorTable(len(stopArray),
            stopArray,rArray,gArray,bArray,nTabColors
         )

def setStyle():
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderSize(10)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasDefH(700)
  gStyle.SetCanvasDefW(700)

  gStyle.SetPadColor       (0)
  gStyle.SetPadBorderSize  (10)
  gStyle.SetPadBorderMode  (0)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadTopMargin   (0.08)
  gStyle.SetPadLeftMargin  (0.15)
  gStyle.SetPadRightMargin (0.05)
  gStyle.SetPadGridX       (0)
  gStyle.SetPadGridY       (0)
  gStyle.SetPadTickX       (1)
  gStyle.SetPadTickY       (1)

  gStyle.SetFrameFillStyle ( 0)
  gStyle.SetFrameFillColor ( 0)
  gStyle.SetFrameLineColor ( 1)
  gStyle.SetFrameLineStyle ( 0)
  gStyle.SetFrameLineWidth ( 1)
  gStyle.SetFrameBorderSize(10)
  gStyle.SetFrameBorderMode( 0)

  gStyle.SetNdivisions(505)

  gStyle.SetLineWidth(2)
  gStyle.SetHistLineWidth(2)
  gStyle.SetFrameLineWidth(2)
  gStyle.SetLegendFillColor(root.kWhite)
  gStyle.SetLegendFont(42)
  gStyle.SetMarkerSize(1.2)
  gStyle.SetMarkerStyle(20)
  gStyle.SetHistLineColor(1)
 
  gStyle.SetLabelSize(0.040,"X")
  gStyle.SetLabelSize(0.040,"Y")

  gStyle.SetLabelOffset(0.010,"X")
  gStyle.SetLabelOffset(0.010,"Y")
 
  gStyle.SetLabelFont(42,"X")
  gStyle.SetLabelFont(42,"Y")
 
  gStyle.SetTitleBorderSize(0)
  gStyle.SetTitleFont(42)
  gStyle.SetTitleFont(42,"X")
  gStyle.SetTitleFont(42,"Y")

  gStyle.SetTitleSize(0.045,"X")
  gStyle.SetTitleSize(0.045,"Y")
 
  gStyle.SetTitleOffset(1.4,"X")
  gStyle.SetTitleOffset(1.6,"Y")
 
  gStyle.SetTextSize(0.055)
  gStyle.SetTextFont(42)
 
  gStyle.SetOptStat(0)
  setNormalColorTable()
  #gStyle.SetPalette(53)
  
setStyle()

def setHistTitles(hist,xlabel,ylabel,zlabel=None):
    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle(ylabel)
    if zlabel:
      hist.GetZaxis().SetTitle(zlabel)

def setHistRange(hist,xMin,xMax,yMin,yMax):
    hist.GetXaxis().SetRangeUser(xMin,xMax)
    hist.GetYaxis().SetRangeUser(yMin,yMax)

def makeWeightHist(f1,canvas,leg):
  firstHist = True
  canvas.cd()
  canvas.SetLogy()
  colorsList = [1,2,3,4,5,6,7,8]
  nColors = len(colorsList)
  iDir = 0
  leg.Clear()
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  tmpList = []
  for dirName in f1.GetListOfKeys():
    tmpList.append(dirName)
  tmpList.reverse()
  for dirName in tmpList:
    print(dirName.GetName())
    if(re.search(r"data",dirName.GetName())):
	continue
    directory = dirName.ReadObj()
    for histKey in directory.GetListOfKeys():
      if(histKey.GetName()=="hWeight"):
        hist = histKey.ReadObj()
	hist.UseCurrentStyle()
	hist.SetLineColor(colorsList[iDir % nColors])
	hist.SetMarkerColor(colorsList[iDir % nColors])
	allIntegral = hist.Integral(0,hist.GetNbinsX()+1)
	integral = hist.Integral()
	if integral > 0.0:
	  print("Fraction Outside of bounds: %f" % (allIntegral/integral-1.0))
	  #hist.Scale(1.0/allIntegral)
	  hist.Scale(1.0/integral)
	else:
	  leg.AddEntry(hist,dirName.GetName(),"lep")
	if(firstHist):
	  firstHist=False
	  hist.GetYaxis().SetTitle("Fraction of Events")
	  hist.GetXaxis().SetTitle("Event Weight")
	  #hist.GetXaxis().SetRangeUser(0.0,1.0)
	  hist.Draw()
	else:
	  hist.Draw("same")
    iDir += 1
  leg.Draw("same")

class DataMCStack:
  def __init__(self, mcHistList, dataHist, canvas, xtitle, ytitle="", drawStack=True,nDivX=7,xlimits=[],showOverflow=False,lumi=5.0,logy=False,signalsNoStack=[],showCompatabilityTests=True,integralPlot=False,energyStr="8TeV",ylimits=[],ylimitsRatio=[],pullType="",doMCErrors=False,showPullStats=False,yMaxVals=[],yMaxXRanges=[],mcVariations=None,scaleMC2Data=False):
    nBinsX = dataHist.GetNbinsX()
    self.xlimits = xlimits
    self.ylimits = ylimits
    self.logy = logy
    self.nBinsX = nBinsX
    self.dataHist = dataHist
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)
    self.mcVarHist = None
    setYLimitsAuto = getattr(self,"setYLimitsAuto")
    if ytitle=="":
      ytitle="Events/%s" % (getBinWidthStr(dataHist))
    for mcHist in mcHistList:
      #print("nBinsX data: %i, mc: %i" % (nBinsX,mcHist.GetNbinsX()))
      assert(nBinsX == mcHist.GetNbinsX())
    for sigHist in signalsNoStack:
      assert(nBinsX == sigHist.GetNbinsX())

    if integralPlot:
      dataHist = getIntegralHist(dataHist,True)
      self.dataHist = dataHist
      newMcHistList = []
      for i in mcHistList:
        newMcHistList.append(getIntegralHist(i))
      mcHistList = newMcHistList
      newSigHistList = []
      for i in signalsNoStack:
        newSigHistList.append(getIntegralHist(i))
      signalsNoStack = newSigHistList
      ytitle = "Integral of "+ytitle+" #geq X"
    self.signalsNoStack = signalsNoStack
    self.mcHistList = mcHistList
    self.dataHist = dataHist

    self.nDataEvents = dataHist.Integral(0,dataHist.GetNbinsX()+1)
    self.mc2DataSF = 1.
    if scaleMC2Data:
      tmpMCSum = 0.
      for mcHist in mcHistList:
        tmpMCSum += mcHist.Integral(0,mcHist.GetNbinsX()+1)
      self.mc2DataSF = float(self.nDataEvents)/tmpMCSum
      print("DataMC SF: %.2f" % self.mc2DataSF)

    # Make MC Stack/sumHist
    self.stack = root.THStack()
    self.mcSumHist = dataHist.Clone("mcSumHist"+dataHist.GetName())
    self.mcSumHist.Reset()
    for mcHist in mcHistList:
      mcHist.SetMaximum(1e12)
      mcHist.SetMinimum(1e-12)
      mcHist.SetLineColor(mcHist.GetFillColor())
      if showOverflow:
        showHistOverflow(mcHist)
      mcHist.Scale(self.mc2DataSF)
      self.mcSumHist.Add(mcHist)
      self.stack.Add(mcHist)

    if showOverflow:
        showHistOverflow(dataHist)

    self.doMCVariations(mcVariations)

    self.mcSumHist.SetFillColor(root.kGray+3)
    self.mcSumHist.SetFillStyle(3254)
    self.mcSumHist.SetMarkerSize(0)
    if doMCErrors and drawStack:
        self.mcSumHist.SetLineStyle(0)

    self.nMCEvents = self.mcSumHist.Integral(0,self.mcSumHist.GetNbinsX()+1)

    # Get chi^2 Prob Data/MC
    self.normchi2 = dataHist.Chi2Test(self.mcSumHist,"UW CHI2/NDF")
    self.chi2Prob = dataHist.Chi2Test(self.mcSumHist,"UW")
    self.KSProb = dataHist.KolmogorovTest(self.mcSumHist)
    if self.mcVarHist != None:
      self.normchi2 = dataHist.Chi2Test(self.mcVarHist,"UW CHI2/NDF")
      self.chi2Prob = dataHist.Chi2Test(self.mcVarHist,"UW")
      self.KSProb = dataHist.KolmogorovTest(self.mcVarHist)
    if self.chi2Prob < 1e-20:
        self.chi2Prob = 0.0
    if self.KSProb < 1e-20:
        self.KSProb = 0.0

    # Make Pull Hist
    self.pullList = []
    self.pullHist = dataHist.Clone("pullHist"+dataHist.GetName())
    self.pullHist.Reset()
    self.oneGraph = root.TGraph()
    self.oneGraph.SetLineWidth(2)
    self.oneGraph.SetLineStyle(2)
    iGraph = 0
    for i in range(0,self.pullHist.GetNbinsX()+2):
      nData = dataHist.GetBinContent(i)
      nMC = self.mcSumHist.GetBinContent(i)
      error = dataHist.GetBinError(i)
      errorMC = self.mcSumHist.GetBinError(i)
      if self.mcVarHist != None:
        errorMC = self.mcVarHist.GetBinError(i)
      pull = 0.0
      ratio = 0.0
      ratioErr = 0.0
      self.oneGraph.SetPoint(iGraph,dataHist.GetXaxis().GetBinCenter(i),1.0)
      iGraph += 1
      if error != 0.0:
        if pullType=="adrian1":
          pull = (nData -nMC)/nData
        else:
          pull = (nData -nMC)/error
      if pullType=="pullMC":
        if errorMC != 0.0:
          pull = (nData -nMC)/errorMC
        else:
          pull = 0.0
      if nMC != 0.0:
        ratio = nData/nMC
        ratioErr = error/nMC
      if pullType=="ratio":
        self.pullHist.SetBinContent(i,ratio)
        self.pullHist.SetBinError(i,ratioErr)
        #print("nData: {0:.2f} +/- {1:.2f}, nMC: {2:.2f}, ratio: {3:.2f} +/- {4:.2f}".format(nData,error,nMC,ratio,ratioErr))
      else:
        self.pullHist.SetBinContent(i,pull)
        #print("nData: %f, nMC: %f, error: %f, pull: %f" % (nData,nMC,error,pull))
      #pullDistribution
      if pullType == "pullMC":
        if errorMC != 0.0:
          self.pullList.append((nData -nMC)/errorMC)
      else:
        if error != 0.0:
          self.pullList.append((nData -nMC)/error)
    #print getattr(self,"getPullDistributionParams")(self.pullList)

    #Find Maximum y-value
    if xlimits != []:
      self.mcSumHist.GetXaxis().SetRangeUser(*xlimits)
      self.dataHist.GetXaxis().SetRangeUser(*xlimits)
    mcMax = self.mcSumHist.GetMaximum()
    if self.mcVarHist != None:
      mcMax = self.mcSumHist.GetMaximum()
    dataMaxBin = self.dataHist.GetMaximumBin()
    dataMax = dataHist.GetBinContent(dataMaxBin)+dataHist.GetBinError(dataMaxBin)
    ymax = 0.0
    if mcMax > dataMax:
       ymax = mcMax
    else:
       ymax = dataMax
    self.ymax = ymax
  
    #Setup Canvas
    canvas.cd()
    self.pad1Top = 0.98
    self.pad1Bot = 0.30
    self.pad1Right = 0.98
    self.pad1Left = 0.02
    pad1 = root.TPad("pad1"+dataHist.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+dataHist.GetName(),"",0.02,0.01,0.98,0.29,0)
    self.pad1 = pad1
    self.pad2 = pad2
  
    pad1.SetBottomMargin(0.005);
    pad2.SetTopMargin   (0.005);
    pad2.SetBottomMargin(0.33);
    """
    pad1.SetBottomMargin(0.01);
    pad2.SetTopMargin   (0.3);
    pad2.SetBottomMargin(0.33);
    """
    canvas.SetLogy(0)
    pad2.SetLogy(0)
    if logy:
        pad1.SetLogy(1)
    else:
        pad1.SetLogy(0)
  
    pad1.Draw() # Projections pad
    pad2.Draw() # Residuals   pad

    pad1Width = pad1.XtoPixel(pad1.GetX2())
    pad1Height = pad1.YtoPixel(pad1.GetY1())
    pad2Height = pad2.YtoPixel(pad2.GetY1())
    #pad1ToPad2FontScalingFactor = float(pad1Width)/pad2Height
    pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
  
    # Main Pad
    pad1.cd();
    xAxis = None
    yAxis = None
    histForAxis = None
    if len(self.ylimits)==2:
      ylimits[0] += 1e-3
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,self.ylimits[0],self.ylimits[1])
    elif self.logy:
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,0.1,ymax*2.0)
    else:
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,1e-3,ymax*1.05)
    self.histForAxis = histForAxis
    self.histForAxis.Draw()
    self.mcSumHist.Draw("e1same")
    #self.canvas.SaveAs("debug.png")
    if len(self.ylimits)!=2:
      setYLimitsAuto(yMaxXRanges,yMaxVals,self.ymax)
    self.histForAxis.Draw()
    self.histForAxis.GetXaxis().SetTitle("")
    self.histForAxis.GetXaxis().SetLabelSize(0)
    self.histForAxis.GetYaxis().SetTitle(ytitle)
    self.histForAxis.GetYaxis().SetLabelSize(0.050)
    self.histForAxis.GetYaxis().SetTitleSize(0.055)
    self.histForAxis.GetXaxis().SetNdivisions(nDivX)
    self.histForAxis.GetXaxis().SetTitleColor(0)
    self.histForAxis.GetXaxis().SetLabelColor(0)
    if drawStack:
      self.stack.Draw("hist same")
      if doMCErrors:
        if self.mcVarHist != None:
          self.mcVarHist.Draw("e2same")
        self.mcSumHist.Draw("e2same")
      pad1.Update()
    else:
      self.mcSumHist.SetFillColor(856)
      self.mcSumHist.SetLineColor(856)
      self.mcSumHist.SetMarkerColor(856)
      self.mcSumHist.SetFillStyle(1001)
      self.mcSumHist.Draw("histo b")
    for sigHist in signalsNoStack:
      sigHist.Draw("histo same")
    dataHist.Draw("pe same")

    pad1.RedrawAxis() # Updates Axis Lines
  
    # Pulls Pad
    pad2.cd()
    self.pullHist.SetTitle("")
    if xlimits != []:
      self.pullHist.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist.GetXaxis().SetTitle(xtitle)
    self.pullHist.GetXaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetNdivisions(nDivX)
    self.pullHist.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.SetLineColor(root.kBlue)
    self.pullHist.SetLineStyle(1)
    self.pullHist.SetLineWidth(2)
    if pullType=="adrian1":
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{Data}")
    elif pullType=="pullMC":
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{\sigma_{MC}}")
    else:
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{\sigma_{Data}}")
    self.pullHist.GetYaxis().SetTitleSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetTitleOffset(0.75*self.pullHist.GetXaxis().GetTitleOffset())
    self.pullHist.GetYaxis().SetTitleOffset(0.70)
    self.pullHist.SetFillColor(856)
    self.pullHist.SetFillStyle(1001)
    if len(ylimitsRatio) == 2:
      ylimitsRatio[0] += 1e-3
      ylimitsRatio[1] -= 1e-3
      self.pullHist.GetYaxis().SetRangeUser(*ylimitsRatio)

    if pullType=="ratio":
      #pad2.SetGridy(1)
      self.pullHist.GetYaxis().SetTitle("#frac{Data}{MC}")
      self.pullHist.Draw("")
      self.oneGraph.Draw()
      self.pullHist.Draw("same")
    else:
      self.pullHist.Draw("histo")

    if showCompatabilityTests:
      self.problatex = root.TLatex()
      self.problatex.SetNDC()
      self.problatex.SetTextFont(root.gStyle.GetLabelFont())
      self.problatex.SetTextSize(self.pullHist.GetYaxis().GetLabelSize())
      self.problatex.SetTextAlign(12)
      yToDraw = 0.41 #bottom
      yToDraw = 0.92 #top
      #self.problatex.DrawLatex(0.18,yToDraw,"KS Prob: {0:.3g}".format(self.KSProb))
      self.problatex.DrawLatex(0.18,yToDraw,"#chi^{2}/NDF: %.3g" % (self.normchi2))
      self.problatex.DrawLatex(0.18,yToDraw-0.08,"#chi^{2}  Prob: %.3g" % (self.chi2Prob))

    pad2.Update()
    pad2.GetFrame().DrawClone()
    pad2.RedrawAxis() # Updates Axis Lines
  
    canvas.cd()
    #self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

  def getPullDistributionParams(self,pullList):
    pull = root.RooRealVar("pull","pull",-20,20)
    mean = root.RooRealVar("mean","pull Mean",0.0,-20,20)
    sigma = root.RooRealVar("sigma","pull sigma",1.0,0.01,20)
    self.pullGaus = root.RooGaussian("pullGaus","pullGaus",pull,mean,sigma)
    self.pullDS = root.RooDataSet("pullDS","pullDS",root.RooArgSet(pull))
    for i in pullList:
      pull.setVal(i)
      self.pullDS.add(root.RooArgSet(pull))
    self.pullFR = self.pullGaus.fitTo(self.pullDS,PRINTLEVEL)
    self.pullMean = mean
    self.pullSigma = sigma
    meanStr = "<Pull> = %.2f #pm %.2f" % (mean.getVal(), mean.getError())
    sigmaStr = "#sigma(Pull) = %.2f #pm %.2f" % (sigma.getVal(), sigma.getError())

    frame = pull.frame(root.RooFit.Bins(20))
    self.pullDS.plotOn(frame)
    self.pullGaus.plotOn(frame)
    frame.Draw()
    self.canvas.SaveAs("pullDist"+self.dataHist.GetName()+".png")
    return meanStr, sigmaStr

  def getXNDC(self,x):
    minX = self.pad1.GetX1()
    maxX = self.pad1.GetX2()
    result=(x-minX)/(maxX-minX)
    return result
  def getYNDC(self,y):
    minY = self.pad1.GetY1()
    maxY = self.pad1.GetY2()
    result=(y-minY)/(maxY-minY)
    return result
  def getXUser(self,x):
    minX = self.pad1.GetX1()
    maxX = self.pad1.GetX2()
    result=x*(maxX-minX)+minX
    return result
  def getYUser(self,y):
    minY = self.pad1.GetY1()
    maxY = self.pad1.GetY2()
    result=y*(maxY-minY)+minY
    #print "running getYUser with: %.2f" % y
    #print "  minY: %.2f" % minY
    #print "  maxY: %.2f" % maxY
    #print "  result: %.2f" % result
    return result
  def setYLimitsAuto(self,rangesNDC,yNDCLimits,yMaxCurrent):
    #self.canvas.SaveAs("before_"+str(int(time.time()*100))+".png")
    #print("Running setYLimitsAuto...")
    self.pad1.Update()
    self.canvas.Update()
    getXUser = getattr(self,"getXUser")
    getYUser = getattr(self,"getYUser")
    setYLimitsAuto = getattr(self,"setYLimitsAuto")
    self.pad1.cd()
    ranges = [[getXUser(i[0]),getXUser(i[1])] for i in rangesNDC]
    yLimitsScaleFactor = 1.0
    if self.logy:
      yLimitsScaleFactor = 0.75
    yLimits = [getYUser(i)*yLimitsScaleFactor for i in yNDCLimits]
    maxPoints = []
    xAxis = self.mcSumHist.GetXaxis()
    #print("yMaxCurrent: %.2f " % (yMaxCurrent))
    for r,yLim in zip(ranges,yLimits):
      maxY = 0.0
      for i in range(1,xAxis.GetNbins()+1):
        if xAxis.GetBinUpEdge(i) >= r[0] and xAxis.GetBinLowEdge(i) <= r[1]:
          y = self.mcSumHist.GetBinContent(i)
          yErrTmp = self.mcSumHist.GetBinError(i)
          yErr2Tmp = 0.
          if self.mcVarHist != None:
            yErr2Tmp = self.mcVarHist.GetBinError(i)
          y += max(yErrTmp,yErr2Tmp)
          maxY = max(y,maxY)
      maxPoints += [maxY]
    rescale = 0.0
    if self.logy:
      newMaxPoints = []
      for x in maxPoints:
        if x>0.:
          newMaxPoints += [log10(x)]
        else:
          newMaxPoints += [0.]
      maxPoints = newMaxPoints
    for yLim,maxY in zip(yLimits,maxPoints):
      #print("yLim: %.2f maxY: %.2f" % (yLim, maxY))
      if maxY > yLim:
        rescaleTmp = (maxY/yLim)
        if rescaleTmp > rescale:
          rescale = rescaleTmp
    if rescale == 0.0:
        self.ymax = yMaxCurrent*1.1
        return
    if self.logy:
      rescale = 10**rescale*5.
    #print(rescale)
    newYMax = yMaxCurrent*rescale*1.5
    newYMin = 1e-3
    if self.logy:
      newYMin = 0.1
    self.histForAxis = root.TH2F(self.histForAxis.GetName()+"ForAxis","",1,self.xlimits[0],self.xlimits[1],1,newYMin,newYMax)
    self.histForAxis.Draw("")
    self.mcSumHist.Draw("e1 same")
    #self.canvas.SaveAs("after_"+str(int(time.time()*100))+".png")
    setYLimitsAuto(rangesNDC,yNDCLimits,newYMax)

  def doMCVariations(self,mcVariations):
    self.mcVarHist = None
    if mcVariations==None:
      return
    for key in mcVariations:
      for hist in mcVariations[key]:
        hist.Scale(self.mc2DataSF)
    errorTypes = set()
    for key in mcVariations:
      key = re.sub("Up$","",key)
      key = re.sub("Down$","",key)
      if not key in errorTypes:
        errorTypes.add(key)
    mcSumVariations = {}
    for key in mcVariations:
      if len(mcVariations[key])==0:
        continue
      sumHist = mcVariations[key][0].Clone()
      sumHist.Reset()
      for h in mcVariations[key]:
        sumHist.Add(h)
      mcSumVariations[key] = sumHist
    self.mcVarHist = self.mcSumHist.Clone(self.mcSumHist.GetName()+"_mcVariations")
    for iBin in range(1,self.mcVarHist.GetNbinsX()+1):
      nom = self.mcVarHist.GetBinContent(iBin)
      err2 = self.mcVarHist.GetBinError(iBin)**2
      for eBase in errorTypes:
        errUp = mcSumVariations[eBase+"Up"].GetBinContent(iBin)
        errDown = mcSumVariations[eBase+"Down"].GetBinContent(iBin)
        errUp = abs(nom-errUp)
        errDown = abs(nom-errDown)
        if errUp > errDown:
            err2 += errUp**2
        else:
            err2 += errDown**2
      err = sqrt(err2)
      self.mcVarHist.SetBinError(iBin,err)
    self.mcVarHist.SetFillColor(root.kRed)
    self.mcVarHist.SetFillStyle(3245)
    self.mcVarHist.SetMarkerSize(0)
    self.mcVarHist.SetLineStyle(0)


class CompareTwoHists:
  def __init__(self, hist1,hist2, canvas, xtitle, ytitle="Events",nDivX=7,nDivPullY=5,xlimits=[],ylimits=[],pullHistRangeY=[0.0,2.0],energyStr="8TeV",lumi=19.4):
    nBinsX = hist1.GetNbinsX()
    assert(nBinsX == hist2.GetNbinsX())
    self.nBinsX = nBinsX
    self.hist1 = hist1
    self.hist2 = hist2
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)

    if xlimits != []:
      self.hist1.GetXaxis().SetRangeUser(*xlimits)
      self.hist2.GetXaxis().SetRangeUser(*xlimits)
  
    # Make Pull Hist
    self.pullHist = hist1.Clone("pullHist"+hist1.GetName())
    self.pullHist.Reset()
    self.pullHist.SetLineColor(hist2.GetLineColor())
    self.pullHist.SetMarkerColor(hist2.GetMarkerColor())
    self.pullErrorBand = root.TGraphAsymmErrors()
    self.pullErrorBand.SetFillColor(856)
    self.pullErrorBand.SetFillStyle(1001)
    self.pullErrorBand.SetLineStyle(2)
    self.pullErrorBand.SetLineColor(root.kBlack)
    for i in range(0,nBinsX+2):
      nhist1 = hist1.GetBinContent(i)
      nhist2 = hist2.GetBinContent(i)
      nhist1Err = hist1.GetBinError(i)
      nhist2Err = hist2.GetBinError(i)
      ratio = 0.0
      ratioErr = 0.0
      if nhist1 != 0.0 and nhist2 != 0.0:
        ratio = nhist2/nhist1
        ratioErr = nhist2/nhist1 * sqrt((nhist1Err/nhist1)**2+(nhist2Err/nhist2)**2)
      self.pullHist.SetBinContent(i,ratio)
      self.pullHist.SetBinError(i,ratioErr)
      tmpAxis = hist1.GetXaxis()
      self.pullErrorBand.SetPoint(i,tmpAxis.GetBinCenter(i),1.0)
      if nhist1 != 0.0:
        self.pullErrorBand.SetPointError(i,tmpAxis.GetBinLowEdge(i),tmpAxis.GetBinUpEdge(i),
                                                            nhist1Err/nhist1,nhist1Err/nhist1)
      else:
        self.pullErrorBand.SetPointError(i,tmpAxis.GetBinLowEdge(i),tmpAxis.GetBinUpEdge(i),0.,0.)
      #print("nData: %f, nMC: %f, error: %f, pull: %f" % (nData,nMC,error,pull))

    firstVizBin = self.pullHist.GetXaxis().GetFirst()
    lastVizBin = self.pullHist.GetXaxis().GetLast()
    for i in range(0,firstVizBin):
      self.pullErrorBand.SetPointEYhigh(i,self.pullErrorBand.GetErrorYhigh(firstVizBin))
      self.pullErrorBand.SetPointEYlow(i,self.pullErrorBand.GetErrorYlow(firstVizBin))
    for i in range(lastVizBin+1,nBinsX+2):
      self.pullErrorBand.SetPointEYhigh(i,self.pullErrorBand.GetErrorYhigh(lastVizBin))
      self.pullErrorBand.SetPointEYlow(i,self.pullErrorBand.GetErrorYlow(lastVizBin))

    #Find Maximum y-value
    max1 = hist1.GetMaximum()
    max2 = hist2.GetMaximum()
    max1Bin = hist1.GetMaximumBin()
    max2Bin = hist2.GetMaximumBin()
    max1 = hist1.GetBinContent(max1Bin)+hist1.GetBinError(max1Bin)
    max2 = hist2.GetBinContent(max2Bin)+hist2.GetBinError(max2Bin)
    ymax = 0.0
    if max1 > max2:
       ymax = max1
    else:
       ymax = max2
    self.ymax = ymax

    #Setup Canvas
    canvas.cd()
    canvas.Clear()
    pad1 = root.TPad("pad1"+hist1.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+hist1.GetName(),"",0.02,0.01,0.98,0.29,0)
    self.pad1 = pad1
    self.pad2 = pad2
  
    pad1.SetBottomMargin(0.005);
    pad2.SetTopMargin   (0.005);
    pad2.SetBottomMargin(0.33);
  
    pad1.Draw() # Projections pad
    pad2.Draw() # Residuals   pad

    pad1Width = pad1.XtoPixel(pad1.GetX2())
    pad1Height = pad1.YtoPixel(pad1.GetY1())
    pad2Height = pad2.YtoPixel(pad2.GetY1())
    #pad1ToPad2FontScalingFactor = float(pad1Width)/pad2Height
    pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
  
    # Main Pad
    pad1.cd();
    self.hist2.GetXaxis().SetTitle("")
    self.hist2.GetXaxis().SetLabelSize(0)
    self.hist2.GetYaxis().SetTitle(ytitle)
    self.hist2.GetYaxis().SetLabelSize(0.050)
    self.hist2.GetYaxis().SetTitleSize(0.055)
    self.hist2.GetXaxis().SetNdivisions(nDivX)
    if ylimits==[]:
      self.hist2.GetYaxis().SetRangeUser(0.0,ymax*1.04)
    else:
      self.hist2.GetYaxis().SetRangeUser(*ylimits)
    self.hist2.Draw("")
    self.hist1.Draw("same")
  
    # Pulls Pad
    pad2.cd()
    self.pullHist.SetTitle("")
    if xlimits != []:
      self.pullHist.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist.GetXaxis().SetTitle(xtitle)
    self.pullHist.GetXaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetNdivisions(nDivX)
    self.pullHist.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetTitle("Ratio")
    self.pullHist.GetYaxis().SetTitleSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().CenterTitle(1)
    self.pullHist.GetYaxis().SetTitleOffset(0.5)
    self.pullHist.GetYaxis().SetNdivisions(nDivPullY)
    self.pullHist.GetYaxis().SetRangeUser(*pullHistRangeY)
    self.pullHist.GetXaxis().SetTitleOffset(0.75*self.pullHist.GetXaxis().GetTitleOffset())
    self.pullHist.GetYaxis().SetTitleOffset(0.70)
    self.pullHist.SetFillColor(856)
    self.pullHist.SetFillStyle(1001)
    self.pullHist.Draw("")
    self.pullErrorBand.Draw("3")
    self.pullErrorBand.Draw("LX")
    #self.pullErrorBandLine = self.pullErrorBand.Clone(self.pullErrorBand.GetName()+"Line")
    #self.pullErrorBandLine.SetFillStyle(0)
    #self.pullErrorBandLine.Draw("same HIST L")
    self.pullHist.Draw("same")

    pad1.RedrawAxis() # Updates Axis Lines
    pad2.RedrawAxis() # Updates Axis Lines
  
    #canvas.cd()
    pad1.cd()
    #self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

class CompareTwoHistsAndData:
  def __init__(self, hist1,hist2, data, canvas, xtitle, ytitle="Events",nDivX=7,nDivPullY=5,xlimits=[],ylimits=[],pullHistRangeY=[0.0,2.0],isPreliminary=True,is7TeV=False,lumi=5.0,logy=False,integralPlot=False,energyStr="8TeV"):
    nBinsX = hist1.GetNbinsX()
    assert(nBinsX == hist2.GetNbinsX())
    assert(nBinsX == data.GetNbinsX())
    self.nBinsX = nBinsX
    self.hist1 = hist1
    self.hist2 = hist2
    self.data = data
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)

    if xlimits != []:
      self.hist1.GetXaxis().SetRangeUser(*xlimits)
      self.hist2.GetXaxis().SetRangeUser(*xlimits)
  
    # Make Pull Hist
    self.pullHist1 = hist1.Clone("pullHist1"+data.GetName())
    self.pullHist1.Reset()
    self.pullHist2 = hist2.Clone("pullHist2"+data.GetName())
    self.pullHist2.Reset()
    for i in range(1,self.pullHist1.GetNbinsX()):
      nData = data.GetBinContent(i)
      nMC1 = hist1.GetBinContent(i)
      nMC2 = hist2.GetBinContent(i)
      error = data.GetBinError(i)
      pull1 = 0.0
      pull2 = 0.0
      if error != 0.0:
        pull1 = (nMC1 - nData)/error
        pull2 = (nMC2 - nData)/error
      self.pullHist1.SetBinContent(i,pull1)
      self.pullHist2.SetBinContent(i,pull2)

    #Find Maximum y-value
    max1 = hist1.GetMaximum()
    max2 = hist2.GetMaximum()
    max1Bin = hist1.GetMaximumBin()
    max2Bin = hist2.GetMaximumBin()
    max1 = hist1.GetBinContent(max1Bin)+hist1.GetBinError(max1Bin)
    max2 = hist2.GetBinContent(max2Bin)+hist2.GetBinError(max2Bin)
    ymax = 0.0
    if max1 > max2:
       ymax = max1
    else:
       ymax = max2
    self.ymax = ymax

    #Find min/max pulls
    max1 = self.pullHist1.GetMaximum()
    max2 = self.pullHist2.GetMaximum()
    pullmax = 0.0
    if max1 > max2:
       pullmax = max1
    else:
       pullmax = max2
    self.pullmax = pullmax
    min1 = self.pullHist1.GetMinimum()
    min2 = self.pullHist2.GetMinimum()
    pullmin = 0.0
    if min1 < min2:
       pullmin = min1
    else:
       pullmin = min2
    self.pullmin = pullmin

    #Setup Canvas
    canvas.cd()
    canvas.Clear()
    pad1 = root.TPad("pad1"+hist1.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+hist1.GetName(),"",0.02,0.01,0.98,0.29,0)
  
    pad1.SetBottomMargin(0.005);
    pad2.SetTopMargin   (0.005);
    pad2.SetBottomMargin(0.33);
  
    pad1.Draw() # Projections pad
    pad2.Draw() # Residuals   pad

    pad1Width = pad1.XtoPixel(pad1.GetX2())
    pad1Height = pad1.YtoPixel(pad1.GetY1())
    pad2Height = pad2.YtoPixel(pad2.GetY1())
    #pad1ToPad2FontScalingFactor = float(pad1Width)/pad2Height
    pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
  
    # Main Pad
    pad1.cd();
    self.hist2.SetTitle("")
    self.hist2.GetXaxis().SetTitle("")
    self.hist2.GetXaxis().SetLabelSize(0)
    self.hist2.GetYaxis().SetTitle(ytitle)
    self.hist2.GetYaxis().SetLabelSize(0.050)
    self.hist2.GetYaxis().SetTitleSize(0.055)
    self.hist2.GetXaxis().SetNdivisions(nDivX)
    if ylimits==[]:
      self.hist2.GetYaxis().SetRangeUser(0.0,ymax*1.04)
    else:
      self.hist2.GetYaxis().SetRangeUser(*ylimits)
    self.hist2.SetFillStyle(0)
    self.hist1.SetFillStyle(0)
    self.hist2.Draw("hist")
    self.hist1.Draw("hist same")
    self.data.Draw("pe same")
  
    # Pulls Pad
    pad2.cd()
    self.pullHist1.SetTitle("")
    if xlimits != []:
      self.pullHist1.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist0 = self.pullHist1.Clone("pullHist0")
    self.pullHist0.Reset()
    self.pullHist0.SetLineColor(1)
    self.pullHist0.SetLineStyle(2)
    self.pullHist0.SetFillStyle(0)
    self.pullHist1.GetXaxis().SetTitle(xtitle)
    self.pullHist1.GetXaxis().CenterTitle(1)
    self.pullHist1.GetXaxis().SetNdivisions(nDivX)
    self.pullHist1.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().SetTitle("#frac{MC-Data}{Error}")
    self.pullHist1.GetYaxis().SetTitleSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().CenterTitle(1)
    self.pullHist1.GetYaxis().SetTitleOffset(0.5)
    self.pullHist1.GetYaxis().SetNdivisions(nDivPullY)
    self.pullHist1.GetYaxis().SetRangeUser(pullmin*0.90,pullmax*1.1)
    self.pullHist1.GetXaxis().SetTitleOffset(0.75*self.pullHist1.GetXaxis().GetTitleOffset())
    #self.pullHist1.GetYaxis().SetTitleOffset(0.70)
    self.pullHist1.SetFillStyle(0)
    self.pullHist2.SetFillStyle(0)
    self.pullHist1.Draw("hist")
    self.pullHist0.Draw("hist same")
    self.pullHist1.Draw("hist same")
    self.pullHist2.Draw("hist same")
    pad2.Update()
    pad2.GetFrame().DrawClone()
  
    canvas.cd()
    #self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

def makeBootstrapHist(hist,outHist,entries=None):
 outHist.Reset()
 samples = entries
 if samples == None:
   integral = hist.Integral()
 for i in range(samples):
   outHist.Fill(hist.GetRandom())

def sqrtThisHistogram(hist):
    """
        Sqrt's bin contents 
        of the input hist bin contents, properly treating the errors.
    """
    nBins = hist.GetNbinsX()

    for i in range(nBins+2):
      y = hist.GetBinContent(i)
      yErr = hist.GetBinError(i)
      if y < 0.0:
        print("Warning sqrtThisHIstogram: hist named %s bin %i has negative y value %f" % (hist.GetName(),i,y))
        hist.SetBinContent(i,0.0)
        hist.SetBinError(i,0.0)
        continue
      if y == 0.0:
        print("Warning sqrtThisHIstogram: hist named %s bin %i has zero y value" % (hist.GetName(),i))
        hist.SetBinContent(i,0.0)
        hist.SetBinError(i,0.0)
        continue
      hist.SetBinContent(i,sqrt(y))
      hist.SetBinError(i,yErr/(sqrt(2*y)))

def getSqrtCopyOfHistogram(hist):
    """
        Reterns a histogram of where the bin contents are the sqrt
        of the input hist bin contents, properly treating the errors.
    """
    outHist = hist.Clone(hist.GetName()+"SqrtHist")
    sqrtThisHistogram(outHist)
    return outHist

def drawSilly(isPreliminary=True,is7TeV=False):
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(22)
    if isPreliminary:
      tlatex.DrawLatex(0.33,0.96,"CMS Preliminary")
    if is7TeV:
      tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=8 TeV, L=4.7 fb^{-1}")

def normalizeHist(hist):
  integral = hist.Integral(0,hist.GetNbinsX()+1)
  if integral != 0.0:
    hist.Scale(1.0/integral)

def showHistOverflow(hist):
  nBins = hist.GetNbinsX()

  overflow = hist.GetBinContent(nBins+1)
  overflowErr = hist.GetBinError(nBins+1)
  lastBin = hist.GetBinContent(nBins)
  lastBinErr = hist.GetBinError(nBins)

  hist.SetBinContent(nBins,lastBin+overflow)
  hist.SetBinError(nBins,sqrt(lastBinErr**2+overflowErr**2))

  underflow = hist.GetBinContent(0)
  underflowErr = hist.GetBinError(0)
  firstBin = hist.GetBinContent(1)
  firstBinErr = hist.GetBinError(1)

  hist.SetBinContent(1,firstBin+underflow)
  hist.SetBinError(1,sqrt(firstBinErr**2+underflowErr**2))

class PlotOfSlices:
  def __init__(self, hist2D, xtitle, ytitle, canvas, xlimits=[], ylimits=[],sliceLabelPrefix="",isPreliminary=True,is7TeV=False):
    canvas.cd(0)
    canvas.Clear()
    nBinsX = hist2D.GetNbinsX()
    nBinsY = hist2D.GetNbinsY()
    self.nBinsX = nBinsX
    self.nBinsY = nBinsY
    self.hist2D = hist2D
    self.canvas = canvas
    self.sliceLabelPrefix = sliceLabelPrefix
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)
    self.histList = []

    colorsListTmp = [root.kRed+1,root.kBlue+1,root.kGreen+1,root.kCyan,root.kMagenta+1]
    self.colorsList=[]
    for i in [0,-11,+2]:
        for j in range(len(colorsListTmp)):
            self.colorsList.append(colorsListTmp[j]+i)

    if xlimits != []:
      self.hist2D.GetXaxis().SetRangeUser(*xlimits)
    if ylimits != []:
      self.hist2D.GetYaxis().SetRangeUser(*ylimits)

    ymax = 0.0
    for i in range(nBinsX+2):
        tmpHist = root.TH1F(hist2D.GetName()+"_slice"+str(i),"",
                            nBinsY,hist2D.GetYaxis().GetXbins().GetArray())
        for j in range(nBinsY+2):
            tmpHist.SetBinContent(j,hist2D.GetBinContent(i,j))
        tmpMax = tmpHist.GetMaximum()
        if tmpMax > ymax:
            ymax = tmpMax
        tmpHist.SetLineColor(self.colorsList[i])
        self.histList.append(tmpHist)
    
    firstHist = self.histList[0]
    firstHist.SetTitle("")
    firstHist.GetXaxis().SetTitle(xtitle)
    firstHist.GetYaxis().SetTitle(ytitle)
    if xlimits != []:
        firstHist.GetXaxis().SetRangeUser(*xlimits)
    if ylimits==[]:
        firstHist.GetYaxis().SetRangeUser(0.0,ymax*1.05)
    else:
        firstHist.GetYaxis().SetRangeUser(*ylimits)

    firstHist.Draw("")
    for hist in self.histList[1:]:
        hist.Draw("same")

    if isPreliminary:
      self.tlatex.DrawLatex(0.33,0.96,"CMS Preliminary")
    if is7TeV:
      self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=8 TeV, L=4.7 fb^{-1}")

    ## Lgend Part

    leg = root.TLegend(0.6,0.3,0.9,0.9)
    leg.SetLineColor(root.kWhite)
    leg.SetFillColor(root.kWhite)
    self.leg = leg
    xAxis = self.hist2D.GetXaxis()
    xBin = 0
    for hist in self.histList:
      tmpLabel = ""
      if xBin == 0:
        tmpLabel = "%s [0.0,%.1f]" % (sliceLabelPrefix,xAxis.GetBinUpEdge(xBin))
      elif xBin == nBinsX+1:
        tmpLabel = "%s [%.1f,#infty]" % (sliceLabelPrefix,xAxis.GetBinLowEdge(xBin))
      else:
        tmpLabel = "%s [%.1f,%.1f]" % (sliceLabelPrefix,xAxis.GetBinLowEdge(xBin),xAxis.GetBinUpEdge(xBin))
      leg.AddEntry(hist,tmpLabel,"l")
      xBin += 1
    leg.Draw("same")

def getIntegralHist(hist,setErrors=True):
  result = hist.Clone(hist.GetName()+"_Integral")
  if hist.InheritsFrom("TH2"):
    nBinsX = result.GetNbinsX()
    nBinsY = result.GetNbinsY()
    for iX in range(nBinsX+1):
      for iY in range(nBinsY+1):
        sumw = 0.0
        sumw2 = 0.0
        for jX in range(iX,nBinsX+2):
          for jY in range(iY,nBinsY+2):
            sumw += result.GetBinContent(jX,jY)
            sumw2 += (result.GetBinError(jX,jY))**2
        result.SetBinContent(iX,iY,sumw)
        if setErrors:
            result.SetBinError(iX,iY,sumw2**0.5)
  else:
    nBins = result.GetNbinsX()
    for i in range(nBins+1):
      sumw = 0.0
      sumw2 = 0.0
      for j in range(i,nBins+2):
        sumw += result.GetBinContent(j)
        sumw2 += (result.GetBinError(j))**2
      result.SetBinContent(i,sumw)
      if setErrors:
          result.SetBinError(i,sumw2**0.5)
  return result

def hist2to1(hist):
  assert(hist.InheritsFrom("TH1"))
  result = None
  nBinsX = hist.GetNbinsX()
  nBinsY = hist.GetNbinsY()
  totalBins = (nBinsX+2)*(nBinsY+2) - 2 #include underflow/overflow
  if hist.InheritsFrom("TH2F"):
    result = root.TH1F(hist.GetName()+"_1d","",totalBins,0,totalBins)
  elif hist.InheritsFrom("TH2D"):
    result = root.TH1D(hist.GetName()+"_1d","",totalBins,0,totalBins)
  else:
    print("Error: hist2to1: Input hist must be TH2F or TH2D, exiting.")
    sys.exit(1)
  k = 0
  for i in range(nBinsX+2):
    for j in range(nBinsY+2):
      tmp = hist.GetBinContent(i,j)
      tmpErr = hist.GetBinError(i,j)
      result.SetBinContent(k,tmp)
      result.SetBinError(k,tmpErr)
      k += 1
  return result

def hist2to1CollapseY(hist,xcuts=[]):
  assert(hist.InheritsFrom("TH1"))
  result = None
  nBinsX = hist.GetNbinsX()
  nBinsY = hist.GetNbinsY()
  ymin = hist.GetYaxis().GetXmin()
  ymax = hist.GetYaxis().GetXmax()
  totalBins = (nBinsX+2)*(nBinsY+2) - 2 #include underflow/overflow
  if hist.InheritsFrom("TH2F"):
    result = root.TH1F(hist.GetName()+"_1d","",nBinsY,ymin,ymax)
  elif hist.InheritsFrom("TH2D"):
    result = root.TH1D(hist.GetName()+"_1d","",nBinsY,ymin,ymax)
  else:
    print("Error: hist2to1CollapseY: Input hist must be TH2F or TH2D, exiting.")
    sys.exit(1)
  minBinX = 0
  maxBinX = nBinsX+2
  if len(xcuts)==2:
    minBinX = hist.GetXaxis().FindBin(xcuts[0])
    maxBinX = hist.GetXaxis().FindBin(xcuts[1])
    if hist.GetXaxis().GetBinCenter(maxBinX)> xcuts[1]:
        maxBinX -= 1
  for j in range(nBinsY+2):
    tmpSum = 0.0
    tmpSumErr2 = 0.0
    for i in range(minBinX,maxBinX):
      tmp = hist.GetBinContent(i,j)
      tmpErr = hist.GetBinError(i,j)
      tmpSum += tmp
      tmpSumErr2 += tmpErr*tmpErr
    result.SetBinContent(j,tmpSum)
    result.SetBinError(j,sqrt(tmpSumErr2))
  return result

def shrinkTH1(hist,xlow,xhigh,deleteOld=False):
  assert(hist.InheritsFrom("TH1"))
  taxis=hist.GetXaxis()
  oldXlow=taxis.GetXmin()
  oldXhigh=taxis.GetXmax()
  assert(xlow >= oldXlow)
  assert(xhigh <= oldXhigh)
  lowBin = taxis.FindBin(xlow)
  highBin = taxis.FindBin(xhigh)
  if taxis.GetBinLowEdge(highBin)==float(xhigh):
    highBin -= 1
  xlow = taxis.GetBinLowEdge(lowBin)
  xhigh = taxis.GetBinUpEdge(highBin)
  oldN = hist.GetNbinsX()
  newN = int((xhigh-xlow)/(oldXhigh-oldXlow)*oldN)
  name = hist.GetName()
  title = hist.GetTitle()
  hist.SetName(name+"_Old")
  newHist = root.TH1F(name,title,newN,xlow,xhigh)
  newHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  newHist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
  for iOld,iNew in zip(range(lowBin,highBin+1),range(1,newN+1)):
    newHist.SetBinContent(iNew,hist.GetBinContent(iOld))
    newHist.SetBinError(iNew,hist.GetBinError(iOld))
  if deleteOld:
    hist.Delete()
  return newHist

def toyHistogram(hist):
  nBins = hist.GetNbinsX()
  random = root.TRandom3()
  for i in range(nBins+2):
    mean = hist.GetBinContent(i)
    n = random.Poisson(mean)
    err = sqrt(n)
    hist.SetBinContent(i,n)
    hist.SetBinError(i,err)

def getXbinsHighLow(hist,low,high):
  axis = hist.GetXaxis()
  xbinLow = axis.FindBin(low)
  xbinHigh = axis.FindBin(high)
  #print("xbinhigh: {0}, {1}, {2}".format(xbinHigh,axis.GetBinLowEdge(xbinHigh),float(high)))
  if axis.GetBinLowEdge(xbinHigh)==float(high):
    xbinHigh -= 1
  return xbinLow, xbinHigh

def getIntegralAll(hist,boundaries=[]):
  xbinLow = None
  xbinHigh = None
  if len(boundaries)==0:
    xbinLow = 0
    xbinHigh = hist.GetXaxis().GetNbins()
  elif len(boundaries)==2:
    xbinLow, xbinHigh = getXbinsHighLow(hist,boundaries[0],boundaries[1])
  else:
    return -1
  if hist.InheritsFrom("TH2"):
    nBinsY = hist.GetYaxis().GetNbins()
    return hist.Integral(xbinLow,xbinHigh,0,nBinsY+1)
  elif hist.InheritsFrom("TH1"):
    return hist.Integral(xbinLow,xbinHigh)
  else:
    return -1

def getIntegralLowHigh(hist,lowBoundaries,highBoundaries):
  lowInt = getIntegralAll(hist,lowBoundaries)
  highInt = getIntegralAll(hist,highBoundaries)
  return lowInt+highInt

def sqrtTH1(hist):
  nBins = hist.GetNbinsX()
  for i in range(nBins+2):
    n = hist.GetBinContent(i)
    nErr = hist.GetBinError(i)
    if n < 0.0:
      n = 0.0
    hist.SetBinContent(i,sqrt(n))
    hist.SetBinError(i,sqrt(nErr))

def saveAs(canvas,name):
  canvas.SaveAs(name+".png")
  canvas.SaveAs(name+".pdf")
  canvas.SaveAs(name+".eps")
  canvas.SaveAs(name+".root")
  canvas.SaveAs(name+".C")

def setLegPos(leg,legPos):
  leg.SetX1NDC(legPos[0])
  leg.SetX2NDC(legPos[2])
  leg.SetY1NDC(legPos[1])
  leg.SetY2NDC(legPos[3])

def getBinWidthStr(hist):
    binWidth = (hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin())/hist.GetXaxis().GetNbins()
    binWidthPrec = "0"
    if binWidth % 1 > 0.0:
      binWidthPrec = "1"
      if binWidth*10 % 1 > 0.0:
        binWidthPrec = "2"
    return ("%."+binWidthPrec+"f") % (binWidth)

def getEfficiencyInterval(passed,total):
  eff = root.TEfficiency()
  nom = 0.
  quant=0.682689492137
  if total>0:
    nom = float(passed)/total
  low = eff.ClopperPearson(int(total),int(passed),quant,False)
  high = eff.ClopperPearson(int(total),int(passed),quant,True)
  return [low,nom,high]

def drawStandardCaptions(canvas,caption,captionleft1="",captionleft2="",captionleft3="",captionright1="",captionright2="",captionright3="",preliminaryString=""):
  tlatex = root.TLatex()
  tlatex.SetNDC()

  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,preliminaryString)

  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(1.0-canvas.GetRightMargin(),0.96,caption)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+canvas.GetLeftMargin(),0.88,captionleft1)
  tlatex.DrawLatex(0.02+canvas.GetLeftMargin(),0.82,captionleft2)
  tlatex.DrawLatex(0.02+canvas.GetLeftMargin(),0.76,captionleft3)
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-canvas.GetRightMargin(),0.88,captionright1)
  tlatex.DrawLatex(0.97-canvas.GetRightMargin(),0.82,captionright2)
  tlatex.DrawLatex(0.97-canvas.GetRightMargin(),0.76,captionright3)
  return tlatex

def copyTreeBranchToNewNameTree(tree,oldBranchName,newBranchName):
  """
  Returns a new tree with the contents of oldBranchName in the old tree, but with
  the branch name newBranchName

  Assumes both branches just contain floats!!
  """
  result = root.TTree(tree.GetName()+"_newNames","")

  newVal = array.array( 'f', [ 0. ] ) # one element array so we get a pointer to the value
  newBranch = result.Branch( newBranchName, newVal, newBranchName+'/F' )

  nEntries = tree.GetEntries()

  for i in range(nEntries):
    tree.GetEntry(i)
    oldVal = getattr(tree,oldBranchName)
    newVal[0] = oldVal
    #print i,newVal[0]
    result.Fill()
  newBranch.SetAddress(0)
  return result

def getHistMax(hist):
  iBin = hist.GetMaximumBin()
  result = hist.GetBinContent(iBin)
  return result

def makeStdAxisHist(histList,logy=False,freeTopSpace=0.5,xlim=[],ylim=[]):
  assert(len(histList)>0)
  assert(len(xlim)==0 or len(xlim)==2)
  assert(len(ylim)==0 or len(ylim)==2)
  multiplier = 1./(1.-freeTopSpace)
  yMin = 0.
  yMax = 0.
  xMin = 1e15
  xMax = -1e15
  for hist in histList:
    histMax = getHistMax(hist)
    yMax = max(yMax,histMax)
    nBins = hist.GetNbinsX()
    xMax = max(xMax,hist.GetXaxis().GetBinUpEdge(nBins))
    xMin = min(xMin,hist.GetBinLowEdge(1))
  if yMax == 0.:
    yMax = 1.
  if logy:
    yMin = 10**(-1)
    yMax = (math.log10(yMax) + 1.)*multiplier - 1.
    yMax = 10**yMax
  else:
    yMax = yMax*multiplier
    if yMax == 0.:
      yMax = 1.
  if len(xlim)==2:
    xMin = xlim[0]
    xMax = xlim[1]
  if len(ylim)==2:
    yMin = ylim[0]
    yMax = ylim[1]
  axisHist = root.TH2F("axisHist"+str(random.randint(1000,1000000)),"",1,xMin,xMax,1,yMin,yMax)
  return axisHist

def getLogBins(nBins,xMin,xMax):
  xMinLog = math.log10(xMin)
  delta = (math.log10(xMax)-xMinLog)/nBins
  return [10**(xMinLog + x*delta) for x in range(nBins+1)]

def drawNormalLegend(hists,labels,option="l"):
  assert(len(hists)==len(labels))
  #leg = root.TLegend(0.55,0.6,0.91,0.89)
  #leg = root.TLegend(0.35,0.6,0.91,0.89)
  leg = root.TLegend(0.40,0.7,0.91,0.89)
  leg.SetLineColor(root.kWhite)
  for hist,label in zip(hists,labels):
    leg.AddEntry(hist,label,option)
  leg.Draw()
  return leg

def setupCOLZFrame(pad,reset=False):
   if reset:
     pad.SetRightMargin(gStyle.GetPadRightMargin())
   else:
     pad.SetRightMargin(0.15)

def normToBinWidth(hist):
  xaxis = hist.GetXaxis()
  nBins = xaxis.GetNbins()
  for i in range(1,nBins+1):
    binContent = hist.GetBinContent(i)
    binWidth = hist.GetBinWidth(i)
    hist.SetBinContent(i,binContent/binWidth)
  return hist

if __name__ == "__main__":

  root.gROOT.SetBatch(True)
  print("Running helpers.py")
  
