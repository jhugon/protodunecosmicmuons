#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as mpl
import numpy
import random

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
    for x in [self.xmin, self.xmax]:
        for y in [self.ymin, self.ymax]:
            result.append([x,y,self.zmid])
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

  def plot(self,ax,lc='-b'):

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

def plotRectangle(ax,xmin,xmax,ymin,ymax,zmin,zmax,lc='-g'):
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

def plotTrack(ax,position,direction,lc='-r'):
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
    ax.plot3D([position[0],endPos[0]],[position[1],endPos[1]],[position[2],endPos[2]],lc)

def findLinePoints(paddle1,paddle2,y=100.):
    corners1 = paddle1.getCorners()
    corners2 = paddle2.getCorners()
    result = []
    angles = []
    xs = []
    zs = []
    for c1 in corners1:
        for c2 in corners2:
            d = [c2[i] - c1[i] for i in range(3)]
            p = c1
            dirScaleFactor = (y - p[1])/float(d[1])
            point = [c1[i]+dirScaleFactor*d[i] for i in range(3)]
            result.append(point)
            angle = abs(numpy.arctan(d[1]/(d[0]**2+d[2]**2)**0.5))
            angleDeg = angle*180/numpy.pi
            zenithAngleDeg = 90-angleDeg
            angles.append(zenithAngleDeg)
            xs.append(point[0])
            zs.append(point[2])
            print "  position: ({:5.1f},{:5.1f},{:5.1f}) zenith angle: {angle:4.1f} deg".format(*point,angle=zenithAngleDeg)

    print("  xmin/max: {:5.1f},{:5.1f} zmin/max: {:5.1f},{:5.1f} angle min/max: {:5.1f},{:5.1f}".format(
                                        min(xs),max(xs),min(zs),max(zs),min(angles),max(angles)))
    return result
    return angles

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


if __name__ == "__main__":

    print "Paddles 1 and 2"
    findLinePoints(cosmic1,cosmic2)
    print "Paddles 3 and 4"
    findLinePoints(cosmic3,cosmic4)

    ################################

    fig = mpl.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    
    #xs = [paddle.getMid()[0] for paddle in paddles]
    #ys = [paddle.getMid()[1] for paddle in paddles]
    #zs = [paddle.getMid()[2] for paddle in paddles]
    #ax.scatter(xs, ys, zs, c='r', marker='o')

    for iPaddle, paddle in enumerate(paddles):
      paddle.plot(ax,['c','m','y','k'][iPaddle])
    
    #plotRectangle(ax,*tpcBoundaries)
    plotRectangle(ax,*tpcActiveBoundaries,lc='b')

    tracks = [
                [[-120.,80.,-190.],[1.5,-1.0,2.5]],
                [[-80.,80.,-190.],[1.,-0.8,2.]],
            ]
    for i in range(len(tracks)):
        s = (sum([j**2 for j in tracks[i][1] ]))**0.5
        print(s)
        tracks[i][1] = [j/s for j in tracks[i][1]]
    for itrk, track in enumerate(tracks):
        print("Track: {:2d}".format(itrk))
        plotTrack(ax,track[0],track[1])
        for ipad, paddle in enumerate(paddles):
            print("  {:2d} {}".format(ipad,paddle.checkWentThrough(*track,fast=True)))
    
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_zlabel('z [cm]')
    
    mpl.show()

