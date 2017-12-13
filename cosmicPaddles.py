#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as mpl
import numpy
import random
import scipy.spatial

class CosmicPaddle(object):
  def __init__(self,bounds):
    self.boundList = bounds
    self.xminTop = bounds[0]
    self.xmaxTop = bounds[1]
    self.xminBot = bounds[2]
    self.xmaxBot = bounds[3]
    self.ymin = bounds[4]
    self.ymax = bounds[5]
    self.zmin = bounds[6]
    self.zmax = bounds[7]

    self.zmid = 0.5*(self.zmin+self.zmax)
    self.ymid = 0.5*(self.ymin+self.ymax)
    self.xmid = 0.5*(self.xminTop+self.xmaxTop)

    self.xmax = max(self.xmaxTop,self.xmaxBot)
    self.xmin = min(self.xminTop,self.xminBot)

    self.xmaxSlope = float(self.ymax-self.ymin)/(self.xmaxTop-self.xmaxBot)
    self.xminSlope = float(self.ymax-self.ymin)/(self.xminTop-self.xminBot)
    assert(self.xmaxSlope > 0)
    assert(self.xminSlope < 0)

    self.xmaxYIntercept = self.ymin - self.xmaxSlope * self.xmaxBot
    self.xminYIntercept = self.ymin - self.xminSlope * self.xminBot

  def getCorners(self):
    result = []
    y = self.ymax
    for x in [self.xminTop, self.xmaxTop]:
        for z in [self.zmin,self.zmax]:
            result.append([x,y,z])
    y = self.ymin
    for x in [self.xminBot, self.xmaxBot]:
        for z in [self.zmin,self.zmax]:
            result.append([x,y,z])
    return result

  def getMid(self):
    return [self.xmid,self.ymid,self.zmid]

  def checkWentThrough(self,position,direction,fast=False):
    for zDet in [self.zmin,self.zmax]:
      dirScaleFactor = (zDet - position[2])/float(direction[2])
      x = position[0] + dirScaleFactor * direction[0]
      y = position[1] + dirScaleFactor * direction[1]
      if y > self.ymin and y < self.ymax:
        if x > self.xmin and x < self.xmax:
          if fast:
            return True
          else:
            xmaxBoundYmin = x * self.xmaxSlope + self.xmaxYIntercept
            xminBoundYmin = x * self.xminSlope + self.xminYIntercept
            if y > xmaxBoundYmin and y > xminBoundYmin:
                return True
    return False

  def checkWentThroughArray(self,x,y,z,px,py,pz,fast=False):
    dirScaleFactor = (self.zmid - z) / pz
    thisx = x + dirScaleFactor * px
    thisy = y + dirScaleFactor * py

    inX = numpy.logical_and(thisx > self.xmin, thisx < self.xmax)
    inY = numpy.logical_and(thisy > self.ymin, thisy < self.ymax)
    inXY = numpy.logical_and(inX,inY)
    if fast:
        return inXY
    else:
        xmaxBoundYmin = thisx * self.xmaxSlope + self.xmaxYIntercept
        xminBoundYmin = thisx * self.xminSlope + self.xminYIntercept
        inBounds = numpy.logical_and(thisy > xmaxBoundYmin, thisy > xminBoundYmin)
        result = numpy.logical_and(inXY,inBounds)
        return result

  def plot(self,ax,lc='-b',swapYZ=False):

    if swapYZ:
        ax.plot3D([self.xminTop,self.xmaxTop],[self.zmin,self.zmin],[self.ymax,self.ymax],lc)
        ax.plot3D([self.xminTop,self.xmaxTop],[self.zmax,self.zmax],[self.ymax,self.ymax],lc)
                                                                                         
        ax.plot3D([self.xminBot,self.xmaxBot],[self.zmin,self.zmin],[self.ymin,self.ymin],lc)
        ax.plot3D([self.xminBot,self.xmaxBot],[self.zmax,self.zmax],[self.ymin,self.ymin],lc)
                                                                                         
        ax.plot3D([self.xminBot,self.xminTop],[self.zmin,self.zmin],[self.ymin,self.ymax],lc)
        ax.plot3D([self.xminBot,self.xminTop],[self.zmax,self.zmax],[self.ymin,self.ymax],lc)
                                                                                         
        ax.plot3D([self.xmaxBot,self.xmaxTop],[self.zmin,self.zmin],[self.ymin,self.ymax],lc)
        ax.plot3D([self.xmaxBot,self.xmaxTop],[self.zmax,self.zmax],[self.ymin,self.ymax],lc)
                                                                                         
        ax.plot3D([self.xminTop,self.xminTop],[self.zmin,self.zmax],[self.ymax,self.ymax],lc)
        ax.plot3D([self.xminBot,self.xminBot],[self.zmin,self.zmax],[self.ymin,self.ymin],lc)
                                                                                         
        ax.plot3D([self.xmaxTop,self.xmaxTop],[self.zmin,self.zmax],[self.ymax,self.ymax],lc)
        ax.plot3D([self.xmaxBot,self.xmaxBot],[self.zmin,self.zmax],[self.ymin,self.ymin],lc)
    else:
        ax.plot3D([self.xminTop,self.xmaxTop],[self.ymax,self.ymax],[self.zmin,self.zmin],lc)
        ax.plot3D([self.xminTop,self.xmaxTop],[self.ymax,self.ymax],[self.zmax,self.zmax],lc)

        ax.plot3D([self.xminBot,self.xmaxBot],[self.ymin,self.ymin],[self.zmin,self.zmin],lc)
        ax.plot3D([self.xminBot,self.xmaxBot],[self.ymin,self.ymin],[self.zmax,self.zmax],lc)

        ax.plot3D([self.xminBot,self.xminTop],[self.ymin,self.ymax],[self.zmin,self.zmin],lc)
        ax.plot3D([self.xminBot,self.xminTop],[self.ymin,self.ymax],[self.zmax,self.zmax],lc)

        ax.plot3D([self.xmaxBot,self.xmaxTop],[self.ymin,self.ymax],[self.zmin,self.zmin],lc)
        ax.plot3D([self.xmaxBot,self.xmaxTop],[self.ymin,self.ymax],[self.zmax,self.zmax],lc)

        ax.plot3D([self.xminTop,self.xminTop],[self.ymax,self.ymax],[self.zmin,self.zmax],lc)
        ax.plot3D([self.xminBot,self.xminBot],[self.ymin,self.ymin],[self.zmin,self.zmax],lc)

        ax.plot3D([self.xmaxTop,self.xmaxTop],[self.ymax,self.ymax],[self.zmin,self.zmax],lc)
        ax.plot3D([self.xmaxBot,self.xmaxBot],[self.ymin,self.ymin],[self.zmin,self.zmax],lc)

def plotRectangle(ax,xmin,xmax,ymin,ymax,zmin,zmax,lc='-g',swapYZ=False):
    if swapYZ:
        ax.plot3D([xmin,xmin],[zmin,zmax],[ymin,ymin],lc)
        ax.plot3D([xmin,xmin],[zmin,zmax],[ymax,ymax],lc)
        ax.plot3D([xmax,xmax],[zmin,zmax],[ymin,ymin],lc)
        ax.plot3D([xmax,xmax],[zmin,zmax],[ymax,ymax],lc)
                                                     
        ax.plot3D([xmin,xmin],[zmin,zmin],[ymin,ymax],lc)
        ax.plot3D([xmin,xmin],[zmax,zmax],[ymin,ymax],lc)
        ax.plot3D([xmax,xmax],[zmin,zmin],[ymin,ymax],lc)
        ax.plot3D([xmax,xmax],[zmax,zmax],[ymin,ymax],lc)
                                                     
        ax.plot3D([xmin,xmax],[zmin,zmin],[ymin,ymin],lc)
        ax.plot3D([xmin,xmax],[zmax,zmax],[ymin,ymin],lc)
        ax.plot3D([xmin,xmax],[zmin,zmin],[ymax,ymax],lc)
        ax.plot3D([xmin,xmax],[zmax,zmax],[ymax,ymax],lc)
    else:
        ax.plot3D([xmin,xmin],[ymin,ymin],[zmin,zmax],lc)
        ax.plot3D([xmin,xmin],[ymax,ymax],[zmin,zmax],lc)
        ax.plot3D([xmax,xmax],[ymin,ymin],[zmin,zmax],lc)
        ax.plot3D([xmax,xmax],[ymax,ymax],[zmin,zmax],lc)

        ax.plot3D([xmin,xmin],[ymin,ymax],[zmin,zmin],lc)
        ax.plot3D([xmin,xmin],[ymin,ymax],[zmax,zmax],lc)
        ax.plot3D([xmax,xmax],[ymin,ymax],[zmin,zmin],lc)
        ax.plot3D([xmax,xmax],[ymin,ymax],[zmax,zmax],lc)

        ax.plot3D([xmin,xmax],[ymin,ymin],[zmin,zmin],lc)
        ax.plot3D([xmin,xmax],[ymin,ymin],[zmax,zmax],lc)
        ax.plot3D([xmin,xmax],[ymax,ymax],[zmin,zmin],lc)
        ax.plot3D([xmin,xmax],[ymax,ymax],[zmax,zmax],lc)

def plotTrack(ax,position,direction,lc='-r',swapYZ=False):
    def inBoundaries(x,y,z):
        if x > 150. or x < -150.:
          return False
        if y > 150. or y < -150.:
          return False
        if x > 200. or z < -200.:
          return False
        return True

    assert(len(position)==len(direction))
    assert(len(position)==3)
    endPos = [0.,0.,0.]
    for val,i in [(-150,1),(150,1),(200,2),(-200,2),(150,0),(-150,0)]:
        dirScaleFactor = (val - position[i])/float(direction[i])
        endPos = [position[i] + dirScaleFactor*direction[i] for i in range(3)]
        if dirScaleFactor > 0 and inBoundaries(*endPos):
            break
    if swapYZ:
        ax.plot3D([position[0],endPos[0]],[position[2],endPos[2]],[position[1],endPos[1]],lc)
    else:
        ax.plot3D([position[0],endPos[0]],[position[1],endPos[1]],[position[2],endPos[2]],lc)

def findLinePoints(paddle1,paddle2,y=100.):
    corners1 = paddle1.getCorners()
    corners2 = paddle2.getCorners()
    result = []
    angles = []
    phis = []
    xs = []
    zs = []
    for c1 in corners1:
        for c2 in corners2:
            d = [c2[i] - c1[i] for i in range(3)]
            p = c1
            dirScaleFactor = (y - p[1])/float(d[1])
            point = [c1[i]+dirScaleFactor*d[i] for i in range(3)]
            result.append(point)
            theta = numpy.arctan2((d[0]**2+d[2]**2)**0.5,d[1])*180/numpy.pi
            phi = numpy.arctan2(d[2],d[0])*180/numpy.pi
            if theta < 90.:
                theta = 180 - theta
                phi += 180
            angle = 180 - theta # want from angle going down
            
            angles.append(angle)
            phis.append(phi)
            xs.append(point[0])
            zs.append(point[2])
            #print "  position: ({:5.1f},{:5.1f},{:5.1f}) zenith angle: {angle:4.1f} deg".format(*point,angle=zenithAngleDeg)

    #print("  xmin/max: {:5.1f},{:5.1f} zmin/max: {:5.1f},{:5.1f} angle min/max: {:5.1f},{:5.1f}".format(
    #                                    min(xs),max(xs),min(zs),max(zs),min(angles),max(angles)))
    return numpy.array(result), numpy.array(angles), numpy.array(phis)


#,xmin Top [cm],xmax Top [cm],xmin Bottom [cm],xmax Bottom [cm],ymin [cm],ymax [cm],zmin [cm],zmax [cm]
Bounds_Cosmic10x1 = [-80.4911,-48.2911,-77.7411,-51.0411,-108.9835,-47.9835,-131.502,-128.482]
Bounds_Cosmic20x1 = [70.9371,103.1371,73.6871,100.3871,17.2977,78.2977,152.904,155.924]
Bounds_Cosmic30x1 = [-80.5593,-47.3593,-77.4593,-50.4593,24.7269,95.5269,-131.995,-128.975]
Bounds_Cosmic40x1 = [67.7564,100.9564,70.8564,97.8564,-82.0192,-11.2192,153.008,156.028]

cosmic1 = CosmicPaddle(Bounds_Cosmic10x1)
cosmic2 = CosmicPaddle(Bounds_Cosmic20x1)
cosmic3 = CosmicPaddle(Bounds_Cosmic30x1)
cosmic4 = CosmicPaddle(Bounds_Cosmic40x1)

paddles = [cosmic1,cosmic2,cosmic3,cosmic4]

tpcBoundaries = [-0.8,49.17,-25,25,-5,95]
tpcActiveBoundaries = [0.4,47.9,-20,20,0,90]

def eventViewer(tracks,listsOfPoints=[]):
    fig = mpl.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    
    for iPaddle, paddle in enumerate(paddles):
      paddle.plot(ax,['c','m','y','k'][iPaddle],swapYZ=True)
    
    plotRectangle(ax,*tpcBoundaries,swapYZ=True)
    plotRectangle(ax,*tpcActiveBoundaries,lc='b',swapYZ=True)

    for itrk, track in enumerate(tracks):
        print("Track: {:2d}".format(itrk))
        plotTrack(ax,track[0],track[1], swapYZ=True)
        for ipad, paddle in enumerate(paddles):
            print("  {:2d} {}".format(ipad,paddle.checkWentThrough(*track,fast=True)))

    for iList, listOfPoints in enumerate(listsOfPoints):
        ax.scatter(listOfPoints[:,0],listOfPoints[:,2],listOfPoints[:,1], c=['r','g','b','k','o'][iList], marker='x')
    
    ax.invert_xaxis()
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('z [cm]')
    ax.set_zlabel('y [cm]')
    
    mpl.show()

def genPositionsAngles(paddleSets,debugPlots=False,nScaleFactor=0.01):
    def genUniformPoints(N,xmin,xmax,ymin,ymax):
        randomPoints = numpy.random.rand(N,2)
        randomPoints[:,0] *= abs(xmax-xmin)
        randomPoints[:,1] *= abs(ymax-ymin)
        randomPoints[:,0] += xmin
        randomPoints[:,1] += ymin
        return randomPoints

    pointSets = []
    angleSets = []
    phiSets = []
    convexHulls = []
    delaunays = []
    xExtremas = []
    zExtremas = []
    areas = []
    thetaExtremas = []
    thetaWidths = []
    thetaMids = []
    phiExtremas = []
    phiWidths = []
    phiMids = []
    for paddleSet in paddleSets:
        pointSet, angleSet, phiSet = findLinePoints(*paddleSet)
        convexHull = scipy.spatial.ConvexHull(pointSet[:,[0,2]])
        delaunay = scipy.spatial.Delaunay(pointSet[convexHull.vertices][:,[0,2]])
        pointSets.append(pointSet)
        angleSets.append(angleSet)
        phiSets.append(phiSet)
        convexHulls.append(convexHull)
        delaunays.append(delaunay)
        xExtrema = [min(pointSet[:,0]),max(pointSet[:,0])]
        zExtrema = [min(pointSet[:,2]),max(pointSet[:,2])]
        area = abs(xExtrema[0]-xExtrema[1])*abs(zExtrema[0]-zExtrema[1])
        xExtremas.append(xExtrema)
        zExtremas.append(zExtrema)
        areas.append(area)
        thetaExtrema = [min(angleSet),max(angleSet)]
        thetaWidth = thetaExtrema[1]-thetaExtrema[0]
        thetaMid = 0.5*(thetaExtrema[1]+thetaExtrema[0])
        phiExtrema = [min(phiSet),max(phiSet)]
        phiWidth = phiExtrema[1]-phiExtrema[0]
        phiMid = 0.5*(phiExtrema[1]+phiExtrema[0])
        thetaExtremas.append(thetaExtrema)
        thetaWidths.append(thetaWidth)
        thetaMids.append(thetaMid)
        phiExtremas.append(phiExtrema)
        phiWidths.append(phiWidth)
        phiMids.append(phiMid)
        if debugPlots:
            print "thetaExtrema: {:3.0f}, {:3.0f} deg".format(*thetaExtrema)
            print "phiExtrema: {:3.0f}, {:3.0f} deg, mid: {phiMid:5.1f} deg, width: {phiWidth:5.1f} deg".format(*phiExtrema,phiWidth=phiWidth,phiMid=phiMid)

    maxPhiWidth = max(phiWidths)*1.02
    halfMaxPhiWidth = 0.5*maxPhiWidth
    maxThetaWidth = max(thetaWidths)*1.02
    halfMaxThetaWidth = 0.5*maxThetaWidth

    result = None
    resultPhis = None
    resultThetas = None
    nTries = 2
    for iTry in range(nTries):
        for delaunay, area, xExtrema, zExtrema, thetaMid, phiMid in zip(delaunays,areas,xExtremas,zExtremas,thetaMids,phiMids):
            Ntry = int(area*nScaleFactor)
            randomPoints = genUniformPoints(Ntry,xExtrema[0],xExtrema[1],zExtrema[0],zExtrema[1])
            randomPointsInDelaunay = delaunay.find_simplex(randomPoints) >= 0
            randomPoints = randomPoints[randomPointsInDelaunay]
            randomThetas = numpy.random.rand(len(randomPoints))*maxThetaWidth+thetaMid-halfMaxThetaWidth
            randomPhis = numpy.random.rand(len(randomPoints))*maxPhiWidth+phiMid-halfMaxPhiWidth
            if result is None:
                result = randomPoints
                resultThetas = randomThetas
                resultPhis = randomPhis
            else:
                result = numpy.concatenate((result,randomPoints),0)
                resultThetas = numpy.concatenate((resultThetas,randomThetas),0)
                resultPhis = numpy.concatenate((resultPhis,randomPhis),0)

    if debugPlots:
        fig, (ax1,ax2) = mpl.subplots(2)
        ax1.set_aspect("equal")
#        ax2.set_aspect("equal")
        for iColor, points, convexHull in zip(range(len(pointSets)),pointSets,convexHulls):
            for simplex in convexHull.simplices:
                ax1.plot(points[:,[0,2]][simplex,0],points[:,[0,2]][simplex,1],'-'+['r','g','b','k'][iColor])
            #ax1.scatter(points.T[0],points.T[2],c=['r','g','b','k'][iColor])

        ax1.scatter(result[:,0],result[:,1],s=10,c='c',edgecolors='none',cmap="brg")
        ax1.set_xlabel('x [cm]')
        ax1.set_ylabel('z [cm]')

        #for phi in resultPhis:
        #    ax2.axvline(phi,ymax=0.1,c='c')
        ax2.scatter(resultPhis,resultThetas,s=10,c='c',edgecolors='none')
        for iColor, angleSet, phiSet in zip(range(len(pointSets)),angleSets,phiSets):
            ax2.scatter(phiSet,angleSet,s=10,c=['r','g','b','k'][iColor],edgecolors='none')
        ax2.set_xlabel(r'$\phi_{zx}$ [deg]')
        ax2.set_ylabel(r'$\theta_z$ [deg]')


        #fig.savefig("projPoints.png")
        mpl.show()

    return result, resultThetas, resultPhis

if __name__ == "__main__":

    import sys

    randomPoints, randomThetas, randomPhis = genPositionsAngles([[cosmic1,cosmic2],[cosmic3,cosmic4]])

    #################################

    tracks = [
                [[-120.,80.,-190.],[1.5,-1.0,2.5]],
                [[-80.,80.,-190.],[1.,-0.8,2.]],
            ]
    #eventViewer(tracks)

    
    
