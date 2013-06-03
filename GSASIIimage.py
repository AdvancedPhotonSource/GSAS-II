# -*- coding: utf-8 -*-
#GSASII image calculations: ellipse fitting & image integration        
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIimage: Image calc module*
================================

Ellipse fitting & image integration

needs some minor refactoring to remove wx code
'''

import math
import wx
import time
import numpy as np
import numpy.linalg as nl
from scipy.optimize import leastsq
import copy
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIplot as G2plt
import GSASIIlattice as G2lat
import GSASIIpwd as G2pwd
import fellipse as fel

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
tand = lambda x: math.tan(x*math.pi/180.)
atand = lambda x: 180.*math.atan(x)/math.pi
atan2d = lambda y,x: 180.*math.atan2(y,x)/math.pi
cosd = lambda x: math.cos(x*math.pi/180.)
acosd = lambda x: 180.*math.acos(x)/math.pi
rdsq2d = lambda x,p: round(1.0/math.sqrt(x),p)
#numpy versions
npsind = lambda x: np.sin(x*np.pi/180.)
npasind = lambda x: 180.*np.arcsin(x)/np.pi
npcosd = lambda x: np.cos(x*np.pi/180.)
nptand = lambda x: np.tan(x*np.pi/180.)
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
def pointInPolygon(pXY,xy):
    'Needs a doc string'
    #pXY - assumed closed 1st & last points are duplicates
    Inside = False
    N = len(pXY)
    p1x,p1y = pXY[0]
    for i in range(N+1):
        p2x,p2y = pXY[i%N]
        if (max(p1y,p2y) >= xy[1] > min(p1y,p2y)) and (xy[0] <= max(p1x,p2x)):
            if p1y != p2y:
                xinters = (xy[1]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x or xy[0] <= xinters:
                Inside = not Inside
        p1x,p1y = p2x,p2y
    return Inside
    
def peneCorr(tth,dep):
    'Needs a doc string'
    return dep*(1.-npcosd(tth))         #best one
#    return dep*npsind(tth)             #not as good as 1-cos2Q
        
def makeMat(Angle,Axis):
    '''Make rotation matrix from Angle and Axis

    :param float Angle: in degrees
    :param int Axis: 0 for rotation about x, 1 for about y, etc.
    '''
    cs = cosd(Angle)
    ss = sind(Angle)
    M = np.array(([1.,0.,0.],[0.,cs,-ss],[0.,ss,cs]),dtype=np.float32)
    return np.roll(np.roll(M,Axis,axis=0),Axis,axis=1)
                    
def FitRing(ring,delta):
    'Needs a doc string'
    parms = []
    if delta:
        err,parms = FitEllipse(ring)
        errc,parmsc = FitCircle(ring)
        errc = errc[0]/(len(ring)*parmsc[2][0]**2)
        if not parms or errc < .1:
            parms = parmsc
    else:
        err,parms = FitCircle(ring)
    return parms
        
def FitCircle(ring):
    'Needs a doc string'
    
    def makeParmsCircle(B):
        cent = [-B[0]/2,-B[1]/2]
        phi = 0.
        sr1 = sr2 = math.sqrt(cent[0]**2+cent[1]**2-B[2])
        return cent,phi,[sr1,sr2]
        
    ring = np.array(ring)
    x = np.asarray(ring.T[0])
    y = np.asarray(ring.T[1])
    
    M = np.array((x,y,np.ones_like(x)))
    B = np.array(-(x**2+y**2))
    result = nl.lstsq(M.T,B)
    return result[1],makeParmsCircle(result[0])
        
def FitEllipse(ring):
    'Needs a doc string'
            
    def makeParmsEllipse(B):
        det = 4.*(1.-B[0]**2)-B[1]**2
        if det < 0.:
            print 'hyperbola!'
            return 0
        elif det == 0.:
            print 'parabola!'
            return 0
        cent = [(B[1]*B[3]-2.*(1.-B[0])*B[2])/det, \
            (B[1]*B[2]-2.*(1.+B[0])*B[3])/det]
        phi = 0.5*atand(0.5*B[1]/B[0])
        
        a = (1.+B[0])/cosd(2*phi)
        b = 2.-a
        f = (1.+B[0])*cent[0]**2+(1.-B[0])*cent[1]**2+B[1]*cent[0]*cent[1]-B[4]
        if f/a < 0 or f/b < 0:
            return 0
        sr1 = math.sqrt(f/a)
        sr2 = math.sqrt(f/b)
        if sr1 > sr2:
            sr1,sr2 = sr2,sr1
            phi -= 90.
            if phi < -90.:
                phi += 180.
        return cent,phi,[sr1,sr2]
                
    ring = np.array(ring)
    x = np.asarray(ring.T[0])
    y = np.asarray(ring.T[1])
    M = np.array((x**2-y**2,x*y,x,y,np.ones_like(x)))
    B = np.array(-(x**2+y**2))
    bb,err = fel.fellipse(len(x),x,y,1.E-7)
#    print nl.lstsq(M.T,B)[0]
#    print bb
    return err,makeParmsEllipse(bb)
    
def FitDetector(rings,varyList,parmDict):
    'Needs a doc string'
        
    def CalibPrint(ValSig):
        print 'Image Parameters:'
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for name,fmt,value,sig in ValSig:
            ptlbls += "%s" % (name.rjust(12))
            if name == 'rotate':
                ptstr += fmt % (value-90.)      #kluge to get rotation from vertical - see GSASIIimgGUI
            else:
                ptstr += fmt % (value)
            if sig:
                sigstr += fmt % (sig)
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr        
        
    def ellipseCalcD(B,xyd,varyList,parmDict):
        x = xyd[0]
        y = xyd[1]
        dsp = xyd[2]
        wave = parmDict['wave']
        if 'dep' in varyList:
            dist,x0,y0,tilt,phi,dep = B[:6]
        else:
            dist,x0,y0,tilt,phi = B[:5]
            dep = parmDict['dep']
        if 'wave' in varyList:
            wave = B[-1]
        tth = 2.0*npasind(wave/(2.*dsp))
        dxy = peneCorr(tth,dep)
        ttth = nptand(tth)
        radius = (dist+dxy)*ttth
        stth = npsind(tth)
        cosb = npcosd(tilt)
        R1 = (dist+dxy)*stth*npcosd(tth)*cosb/(cosb**2-stth**2)
        R0 = np.sqrt(R1*radius*cosb)
        zdis = R1*ttth*nptand(tilt)
        X = x-x0+zdis*npsind(phi)
        Y = y-y0-zdis*npcosd(phi)
        XR = X*npcosd(phi)-Y*npsind(phi)
        YR = X*npsind(phi)+Y*npcosd(phi)
        return (XR/R0)**2+(YR/R1)**2-1
        
    names = ['dist','det-X','det-Y','tilt','phi','dep','wave']
    fmt = ['%12.2f','%12.2f','%12.2f','%12.2f','%12.2f','%12.3f','%12.5f']
    p0 = [parmDict[key] for key in varyList]
    result = leastsq(ellipseCalcD,p0,args=(rings.T,varyList,parmDict),full_output=True)
    parmDict.update(zip(varyList,result[0]))
    vals = list(result[0])
    chi = np.sqrt(np.sum(ellipseCalcD(result[0],rings.T,varyList,parmDict)**2))
    sig = list(chi*np.sqrt(np.diag(result[1])))
    sigList = np.zeros(7)
    for i,name in enumerate(names):
        if name in varyList:
            sigList[i] = sig[varyList.index(name)]
    ValSig = zip(names,fmt,vals,sig)
    CalibPrint(ValSig)
#    try:
#        print 'Trial refinement of wavelength - not used for calibration'
#        p0 = result[0]
#        print p0
#        print parms
#        p0 = np.append(p0,parms[0])
#        resultW = leastsq(ellipseCalcW,p0,args=(rings.T,parms[1:]),full_output=True)
#        resultW[0][3] = np.mod(result[0][3],360.0)          #remove random full circles
#        sig = np.sqrt(np.sum(ellipseCalcW(resultW[0],rings.T)**2))
#        ValSig = zip(names,fmt,resultW[0],sig*np.sqrt(np.diag(resultW[1])))
#        CalibPrint(ValSig)
#        return result[0],resultW[0][-1]        
#    except ValueError:
#        print 'Bad refinement - no result'
#        return result[0],wave
                    
def ImageLocalMax(image,w,Xpix,Ypix):
    'Needs a doc string'
    w2 = w*2
    sizey,sizex = image.shape
    xpix = int(Xpix)            #get reference corner of pixel chosen
    ypix = int(Ypix)
    if (w < xpix < sizex-w) and (w < ypix < sizey-w) and image[ypix,xpix]:
        Z = image[ypix-w:ypix+w,xpix-w:xpix+w]
        Zmax = np.argmax(Z)
        Zmin = np.argmin(Z)
        xpix += Zmax%w2-w
        ypix += Zmax/w2-w
        return xpix,ypix,np.ravel(Z)[Zmax],max(0.0001,np.ravel(Z)[Zmin])   #avoid neg/zero minimum
    else:
        return 0,0,0,0      
    
def makeRing(dsp,ellipse,pix,reject,scalex,scaley,image):
    'Needs a doc string'
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    amin = 0
    amax = 360
    for a in range(0,360,1):
        x = radii[0]*cosd(a)
        y = radii[1]*sind(a)
        X = (cphi*x-sphi*y+cent[0])*scalex      #convert mm to pixels
        Y = (sphi*x+cphi*y+cent[1])*scaley
        X,Y,I,J = ImageLocalMax(image,pix,X,Y)
        if I and J and I/J > reject:
            X += .5                             #set to center of pixel
            Y += .5
            X /= scalex                         #convert to mm
            Y /= scaley
            amin = min(amin,a)
            amax = max(amax,a)
            if [X,Y,dsp] not in ring:
                ring.append([X,Y,dsp])
    delt = amax-amin
    if len(ring) < 10:             #want more than 20 deg
        return [],delt > 90
    return ring,delt > 90
    
def makeIdealRing(ellipse,azm=None):
    'Needs a doc string'
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    if azm:
        aR = azm[0]-90,azm[1]-90,azm[1]-azm[0]
        if azm[1]-azm[0] > 180:
            aR[2] /= 2
    else:
        aR = 0,362,181
        
    a = np.linspace(aR[0],aR[1],aR[2])
    x = radii[0]*npcosd(a-phi)
    y = radii[1]*npsind(a-phi)
    X = (cphi*x-sphi*y+cent[0])
    Y = (sphi*x+cphi*y+cent[1])
    return zip(X,Y)
                
def calcDist(radii,tth):
    'Needs a doc string'
    stth = sind(tth)
    ctth = cosd(tth)
    ttth = tand(tth)
    return math.sqrt(radii[0]**4/(ttth**2*((radii[0]*ctth)**2+(radii[1]*stth)**2)))
    
def calcZdisCosB(radius,tth,radii):
    'Needs a doc string'
    cosB = sinb = radii[0]**2/(radius*radii[1])
    if cosB > 1.:
        return 0.,1.
    else:
        cosb = math.sqrt(1.-sinb**2)
        ttth = tand(tth)
        zdis = radii[1]*ttth*cosb/sinb
        return zdis,cosB
    
def GetEllipse(dsp,data):
    'Needs a doc string'
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    dep = data['DetDepth']
    radii = [0,0]
    tth = 2.0*asind(data['wavelength']/(2.*dsp))
    ttth = tand(tth)
    stth = sind(tth)
    ctth = cosd(tth)
    cosb = cosd(tilt)
    dxy = peneCorr(tth,dep)
    radius = ttth*(dist+dxy)
    radii[1] = (dist+dxy)*stth*ctth*cosb/(cosb**2-stth**2)
    if radii[1] > 0:
        radii[0] = math.sqrt(radii[1]*radius*cosb)
        zdis = radii[1]*ttth*tand(tilt)
        elcent = [cent[0]-zdis*sind(phi),cent[1]+zdis*cosd(phi)]
        return elcent,phi,radii
    else:
        print 'bad ellipse - radii:',radii
        return False
        
def GetDetectorXY(dsp,azm,data):
    'Needs a doc string'
    from scipy.optimize import fsolve
    def func(xy,*args):
       azm,phi,R0,R1,A,B = args
       cp = cosd(phi)
       sp = sind(phi)
       x,y = xy
       out = []
       out.append(y-x*tand(azm))
       out.append(R0**2*((x+A)*sp-(y+B)*cp)**2+R1**2*((x+A)*cp+(y+B)*sp)**2-(R0*R1)**2)
       return out
    elcent,phi,radii = GetEllipse(dsp,data)
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    wave = data['wavelength']
    dist = data['distance']
    dep = data['DetDepth']
    tth = 2.0*asind(wave/(2.*dsp))
    dxy = peneCorr(tth,dep)
    ttth = tand(tth)
    radius = (dist+dxy)*ttth
    stth = sind(tth)
    cosb = cosd(tilt)
    R1 = (dist+dxy)*stth*cosd(tth)*cosb/(cosb**2-stth**2)
    R0 = math.sqrt(R1*radius*cosb)
    zdis = R1*ttth*tand(tilt)
    A = zdis*sind(phi)
    B = -zdis*cosd(phi)
    xy0 = [radius*cosd(azm),radius*sind(azm)]
    xy = fsolve(func,xy0,args=(azm,phi,R0,R1,A,B))+cent
    return xy
    
def GetDetXYfromThAzm(Th,Azm,data):
    'Needs a doc string'
    dsp = data['wavelength']/(2.0*npsind(Th))
    
    return GetDetectorXY(dsp,azm,data)
                    
def GetTthAzmDsp(x,y,data):
    'Needs a doc string'
    wave = data['wavelength']
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    dep = data['DetDepth']
    LRazim = data['LRazimuth']
    azmthoff = data['azmthOff']
    dx = np.array(x-cent[0],dtype=np.float32)
    dy = np.array(y-cent[1],dtype=np.float32)
    X = np.array(([dx,dy,np.zeros_like(dx)]),dtype=np.float32).T
    X = np.dot(X,makeMat(phi,2))
    Z = np.dot(X,makeMat(tilt,0)).T[2]
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z))
    dxy = peneCorr(tth,dep)
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z+dxy))
    dsp = wave/(2.*npsind(tth/2.))
    azm = (npatan2d(dx,-dy)+azmthoff+720.)%360.
    return tth,azm,dsp
    
def GetTth(x,y,data):
    'Needs a doc string'
    return GetTthAzmDsp(x,y,data)[0]
    
def GetTthAzm(x,y,data):
    'Needs a doc string'
    return GetTthAzmDsp(x,y,data)[0:2]
    
def GetDsp(x,y,data):
    'Needs a doc string'
    return GetTthAzmDsp(x,y,data)[2]
       
def GetAzm(x,y,data):
    'Needs a doc string'
    return GetTthAzmDsp(x,y,data)[1]
       
def ImageCompress(image,scale):
    'Needs a doc string'
    if scale == 1:
        return image
    else:
        return image[::scale,::scale]
        
def checkEllipse(Zsum,distSum,xSum,ySum,dist,x,y):
    'Needs a doc string'
    avg = np.array([distSum/Zsum,xSum/Zsum,ySum/Zsum])
    curr = np.array([dist,x,y])
    return abs(avg-curr)/avg < .02

def EdgeFinder(image,data):
    '''this makes list of all x,y where I>edgeMin suitable for an ellipse search?
    '''
    import numpy.ma as ma
    Nx,Ny = data['size']
    pixelSize = data['pixelSize']
    edgemin = data['edgemin']
    scalex = pixelSize[0]/1000.
    scaley = pixelSize[1]/1000.    
    tay,tax = np.mgrid[0:Nx,0:Ny]
    tax = np.asfarray(tax*scalex,dtype=np.float32)
    tay = np.asfarray(tay*scaley,dtype=np.float32)
    tam = ma.getmask(ma.masked_less(image.flatten(),edgemin))
    tax = ma.compressed(ma.array(tax.flatten(),mask=tam))
    tay = ma.compressed(ma.array(tay.flatten(),mask=tam))
    return zip(tax,tay)
    
def ImageRecalibrate(self,data):
    'Needs a doc string'
    import ImageCalibrants as calFile
    print 'Image recalibration:'
    time0 = time.time()
    pixelSize = data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    pixLimit = data['pixLimit']
    cutoff = data['cutoff']
    data['rings'] = []
    data['ellipses'] = []
    if not data['calibrant']:
        print 'no calibration material selected'
        return True
    
    skip = data['calibskip']
    dmin = data['calibdmin']
    Bravais,Cells = calFile.Calibrants[data['calibrant']][:2]
    HKL = []
    for bravais,cell in zip(Bravais,Cells):
        A = G2lat.cell2A(cell)
        hkl = G2lat.GenHBravais(dmin,bravais,A)[skip:]
        HKL += hkl
    HKL = G2lat.sortHKLd(HKL,True,False)
    varyList = ['dist','det-X','det-Y','tilt','phi']
    parmDict = {'dist':data['distance'],'det-X':data['center'][0],'det-Y':data['center'][1],
        'tilt':data['tilt'],'phi':data['rotation'],'wave':data['wavelength'],'dep':data['DetDepth']}
    for H in HKL: 
        dsp = H[3]
        ellipse = GetEllipse(dsp,data)
        Ring,delt = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
#            numZ = len(Ring)
            data['rings'].append(np.array(Ring))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
        else:
            continue
    rings = np.concatenate((data['rings']),axis=0)
    if data['DetDepthRef']:
        varyList.append('dep')
    FitDetector(rings,varyList,parmDict)
    data['distance'] = parmDict['dist']
    data['center'] = [parmDict['det-X'],parmDict['det-Y']]
    data['rotation'] = np.mod(parmDict['phi'],360.0)
    data['tilt'] = parmDict['tilt']
    data['DetDepth'] = parmDict['dep']
    N = len(data['ellipses'])
    data['ellipses'] = []           #clear away individual ellipse fits
    for H in HKL[:N]:
        ellipse = GetEllipse(H[3],data)
        data['ellipses'].append(copy.deepcopy(ellipse+('b',)))
    
    print 'calibration time = ',time.time()-time0
    G2plt.PlotImage(self,newImage=True)        
    return True
            
def ImageCalibrate(self,data):
    'Needs a doc string'
    import copy
    import ImageCalibrants as calFile
    print 'Image calibration:'
    time0 = time.time()
    ring = data['ring']
    pixelSize = data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    pixLimit = data['pixLimit']
    cutoff = data['cutoff']
    if len(ring) < 5:
        print 'not enough inner ring points for ellipse'
        return False
        
    #fit start points on inner ring
    data['ellipses'] = []
    outE = FitRing(ring,True)
    if outE:
        print 'start ellipse:',outE
        ellipse = outE
    else:
        return False
        
    #setup 360 points on that ring for "good" fit
    Ring,delt = makeRing(1.0,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
    if Ring:
        ellipse = FitRing(Ring,delt)
        Ring,delt = makeRing(1.0,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)    #do again
        ellipse = FitRing(Ring,delt)
    else:
        print '1st ring not sufficiently complete to proceed'
        return False
    print 'inner ring:',ellipse     #cent,phi,radii
    data['center'] = copy.copy(ellipse[0])           #not right!! (but useful for now)
    data['ellipses'].append(ellipse[:]+('r',))
    G2plt.PlotImage(self,newImage=True)
    
    #setup for calibration
    data['rings'] = []
    data['ellipses'] = []
    if not data['calibrant']:
        print 'no calibration material selected'
        return True
    
    skip = data['calibskip']
    dmin = data['calibdmin']
    Bravais,Cells = calFile.Calibrants[data['calibrant']][:2]
    HKL = []
    for bravais,cell in zip(Bravais,Cells):
        A = G2lat.cell2A(cell)
        hkl = G2lat.GenHBravais(dmin,bravais,A)[skip:]
        HKL += hkl
    HKL = G2lat.sortHKLd(HKL,True,False)
    wave = data['wavelength']
    cent = data['center']
    elcent,phi,radii = ellipse
    dsp = HKL[0][3]
    tth = 2.0*asind(wave/(2.*dsp))
    ttth = tand(tth)
    data['distance'] = dist = calcDist(radii,tth)
    radius = dist*tand(tth)
    zdis,cosB = calcZdisCosB(radius,tth,radii)
    cent1 = []
    cent2 = []
    xSum = 0
    ySum = 0
    zxSum = 0
    zySum = 0
    phiSum = 0
    tiltSum = 0
    distSum = 0
    Zsum = 0
    for i,H in enumerate(HKL):
        dsp = H[3]
        tth = 2.0*asind(0.5*wave/dsp)
        stth = sind(tth)
        ctth = cosd(tth)
        ttth = tand(tth)
        radius = dist*ttth
        elcent,phi,radii = ellipse
        if i:
            radii[1] = dist*stth*ctth*cosB/(cosB**2-stth**2)
            radii[0] = math.sqrt(radii[1]*radius*cosB)
            zdis,cosB = calcZdisCosB(radius,tth,radii)
            zsinp = zdis*sind(phi)
            zcosp = zdis*cosd(phi)
            cent = data['center']
            elcent = [cent[0]-zsinp,cent[1]+zcosp]
            ellipse = (elcent,phi,radii)
        ratio = radii[1]/radii[0]
        Ring,delt = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
            numZ = len(Ring)
            data['rings'].append(np.array(Ring))
            newellipse = FitRing(Ring,delt)
            elcent,phi,radii = newellipse                
            if abs(phi) > 45. and phi < 0.:
                phi += 180.
            dist = calcDist(radii,tth)
            distR = 1.-dist/data['distance']
            if abs(distR) > 0.1:
#                print dsp,dist,data['distance'],distR,len(Ring),delt
                break
            if distR > 0.001:
                print 'Wavelength too large?'
            elif distR < -0.001:
                print 'Wavelength too small?'
            else:
                ellipse = newellipse
            zdis,cosB = calcZdisCosB(radius,tth,radii)
            Tilt = acosd(cosB)          # 0 <= tilt <= 90
            zsinp = zdis*sind(ellipse[1])
            zcosp = zdis*cosd(ellipse[1])
            cent1.append(np.array([elcent[0]+zsinp,elcent[1]-zcosp]))
            cent2.append(np.array([elcent[0]-zsinp,elcent[1]+zcosp]))
            if i:
                d1 = cent1[-1]-cent1[-2]        #get shift of 2 possible center solutions
                d2 = cent2[-1]-cent2[-2]
                if np.dot(d2,d2) > np.dot(d1,d1):  #right solution is the larger shift
                    data['center'] = cent1[-1]
                else:
                    data['center'] = cent2[-1]
                Zsum += numZ
                phiSum += numZ*phi
                distSum += numZ*dist
                xSum += numZ*data['center'][0]
                ySum += numZ*data['center'][1]
                tiltSum += numZ*abs(Tilt)
                if not np.all(checkEllipse(Zsum,distSum,xSum,ySum,dist,data['center'][0],data['center'][1])):
                    print 'Bad fit for ring # %i. Try reducing Pixel search range'%(i) 
            cent = data['center']
#            print ('for ring # %2i @ d-space %.4f: dist %.3f rotate %6.2f tilt %6.2f Xcent %.3f Ycent %.3f Npts %d' 
#                %(i,dsp,dist,phi-90.,Tilt,cent[0],cent[1],numZ))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
        else:
            break
    G2plt.PlotImage(self,newImage=True)
    fullSize = len(self.ImageZ)/scalex
    if 2*radii[1] < .9*fullSize:
        print 'Are all usable rings (>25% visible) used? Try reducing Min ring I/Ib'
    if not Zsum:
        print 'Only one ring fitted. Check your wavelength.'
        return False
    data['center'] = [xSum/Zsum,ySum/Zsum]
    data['distance'] = distSum/Zsum
    
    #possible error if no. of rings < 3! Might need trap here
    d1 = cent1[-1]-cent1[1]             #compare last ring to 2nd ring
    d2 = cent2[-1]-cent2[1]
    Zsign = 1
    len1 = math.sqrt(np.dot(d1,d1))
    len2 = math.sqrt(np.dot(d2,d2))
    t1 = d1/len1
    t2 = d2/len2
    if len2 > len1:
        if -135. < atan2d(t2[1],t2[0]) < 45.:
            Zsign = -1
    else:
        if -135. < atan2d(t1[1],t1[0]) < 45.:
            Zsign = -1
    
    data['tilt'] = Zsign*tiltSum/Zsum
    data['rotation'] = phiSum/Zsum
    varyList = ['dist','det-X','det-Y','tilt','phi']
    parmDict = {'dist':data['distance'],'det-X':data['center'][0],'det-Y':data['center'][1],
        'tilt':data['tilt'],'phi':data['rotation'],'wave':data['wavelength'],'dep':data['DetDepth']}
    rings = np.concatenate((data['rings']),axis=0)
    if data['DetDepthRef']:
        varyList.append('dep')
        FitDetector(rings,varyList,parmDict)
    data['distance'] = parmDict['dist']
    data['center'] = [parmDict['det-X'],parmDict['det-Y']]
    data['rotation'] = np.mod(parmDict['phi'],360.0)
    data['tilt'] = parmDict['tilt']
    data['DetDepth'] = parmDict['dep']
    N = len(data['ellipses'])
    data['ellipses'] = []           #clear away individual ellipse fits
    for H in HKL[:N]:
        ellipse = GetEllipse(H[3],data)
        data['ellipses'].append(copy.deepcopy(ellipse+('b',)))
    print 'calibration time = ',time.time()-time0
    G2plt.PlotImage(self,newImage=True)        
    return True
    
def Make2ThetaAzimuthMap(data,masks,iLim,jLim):
    'Needs a doc string'
    import numpy.ma as ma
    import polymask as pm
    #transforms 2D image from x,y space to 2-theta,azimuth space based on detector orientation
    pixelSize = data['pixelSize']
    scalex = pixelSize[0]/1000.
    scaley = pixelSize[1]/1000.
    
    tay,tax = np.mgrid[iLim[0]+0.5:iLim[1]+.5,jLim[0]+.5:jLim[1]+.5]         #bin centers not corners
    tax = np.asfarray(tax*scalex,dtype=np.float32)
    tay = np.asfarray(tay*scaley,dtype=np.float32)
    nI = iLim[1]-iLim[0]
    nJ = jLim[1]-jLim[0]
    #make position masks here
    spots = masks['Points']
    polygons = masks['Polygons']
    tam = ma.make_mask_none((nI,nJ))
    for polygon in polygons:
        if polygon:
            tamp = ma.make_mask_none((nI*nJ))
            tamp = ma.make_mask(pm.polymask(nI*nJ,tax.flatten(),
                tay.flatten(),len(polygon),polygon,tamp))
            tam = ma.mask_or(tam.flatten(),tamp)
    if tam.shape: tam = np.reshape(tam,(nI,nJ))
    for X,Y,diam in spots:
        tamp = ma.getmask(ma.masked_less((tax-X)**2+(tay-Y)**2,(diam/2.)**2))
        tam = ma.mask_or(tam,tamp)
    TA = np.array(GetTthAzm(tax,tay,data))
    TA[1] = np.where(TA[1]<0,TA[1]+360,TA[1])
    return np.array(TA),tam           #2-theta & azimuth arrays & position mask

def Fill2ThetaAzimuthMap(masks,TA,tam,image):
    'Needs a doc string'
    import numpy.ma as ma
    Zlim = masks['Thresholds'][1]
    rings = masks['Rings']
    arcs = masks['Arcs']
    TA = np.dstack((ma.getdata(TA[1]),ma.getdata(TA[0])))    #azimuth, 2-theta
    tax,tay = np.dsplit(TA,2)    #azimuth, 2-theta
    for tth,thick in rings:
        tam = ma.mask_or(tam.flatten(),ma.getmask(ma.masked_inside(tay.flatten(),max(0.01,tth-thick/2.),tth+thick/2.)))
    for tth,azm,thick in arcs:
        tamt = ma.getmask(ma.masked_inside(tay.flatten(),max(0.01,tth-thick/2.),tth+thick/2.))
        tama = ma.getmask(ma.masked_inside(tax.flatten(),azm[0],azm[1]))
        tam = ma.mask_or(tam.flatten(),tamt*tama)
    taz = ma.masked_outside(image.flatten(),int(Zlim[0]),Zlim[1])
    tam = ma.mask_or(tam.flatten(),ma.getmask(taz))
    tax = ma.compressed(ma.array(tax.flatten(),mask=tam))
    tay = ma.compressed(ma.array(tay.flatten(),mask=tam))
    taz = ma.compressed(ma.array(taz.flatten(),mask=tam))
    del(tam)
    return tax,tay,taz
    
def ImageIntegrate(image,data,masks):
    'Needs a doc string'
    import histogram2d as h2d
    print 'Begin image integration'
    LUtth = data['IOtth']
    LRazm = np.array(data['LRazimuth'],dtype=np.float64)
    numAzms = data['outAzimuths']
    numChans = data['outChannels']
    Dtth = (LUtth[1]-LUtth[0])/numChans
    Dazm = (LRazm[1]-LRazm[0])/numAzms
    if data['centerAzm']:
        LRazm += Dazm/2.
    NST = np.zeros(shape=(numAzms,numChans),order='F',dtype=np.float32)
    H0 = np.zeros(shape=(numAzms,numChans),order='F',dtype=np.float32)
    imageN = len(image)
    Nx,Ny = data['size']
    nXBlks = (Nx-1)/1024+1
    nYBlks = (Ny-1)/1024+1
    Nup = nXBlks*nYBlks*3+3
    dlg = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
    try:
        t0 = time.time()
        Nup = 0
        dlg.Update(Nup)
        for iBlk in range(nYBlks):
            iBeg = iBlk*1024
            iFin = min(iBeg+1024,Ny)
            for jBlk in range(nXBlks):
                jBeg = jBlk*1024
                jFin = min(jBeg+1024,Nx)                
                print 'Process map block:',iBlk,jBlk,' limits:',iBeg,iFin,jBeg,jFin
                TA,tam = Make2ThetaAzimuthMap(data,masks,(iBeg,iFin),(jBeg,jFin))           #2-theta & azimuth arrays & create position mask
                
                Nup += 1
                dlg.Update(Nup)
                Block = image[iBeg:iFin,jBeg:jFin]
                tax,tay,taz = Fill2ThetaAzimuthMap(masks,TA,tam,Block)    #and apply masks
                del TA,tam
                Nup += 1
                dlg.Update(Nup)
                tax = np.where(tax > LRazm[1],tax-360.,tax)                 #put azm inside limits if possible
                tax = np.where(tax < LRazm[0],tax+360.,tax)
                NST,H0 = h2d.histogram2d(len(tax),tax,tay,taz,numAzms,numChans,LRazm,LUtth,Dazm,Dtth,NST,H0)
                print 'block done'
                del tax,tay,taz
                Nup += 1
                dlg.Update(Nup)
        NST = np.array(NST)
        H0 = np.divide(H0,NST)
        H0 = np.nan_to_num(H0)
        del NST
        if Dtth:
            H2 = np.array([tth for tth in np.linspace(LUtth[0],LUtth[1],numChans+1)])
        else:
            H2 = np.array(LUtth)
        if Dazm:        
            H1 = np.array([azm for azm in np.linspace(LRazm[0],LRazm[1],numAzms+1)])
        else:
            H1 = LRazm
        H0[0] /= npcosd(H2[:-1])**2
        if data['Oblique'][1]:
            H0[0] /= G2pwd.Oblique(data['Oblique'][0],H2[:-1])
        Nup += 1
        dlg.Update(Nup)
        t1 = time.time()
    finally:
        dlg.Destroy()
    print 'Integration complete'
    print "Elapsed time:","%8.3f"%(t1-t0), "s"
    return H0,H1,H2
    
def FitStrSta(Image,StrSta,Controls,Masks):
    'Needs a doc string'
    
#    print 'Masks:',Masks
    StaControls = copy.deepcopy(Controls)
    phi = StrSta['Sample phi']
    wave = Controls['wavelength']
    StaControls['distance'] += StrSta['Sample z']*cosd(phi)
    pixSize = StaControls['pixelSize']
    scalex = 1000./pixSize[0]
    scaley = 1000./pixSize[1]
    rings = []

    for ring in StrSta['d-zero']:       #get observed x,y,d points for the d-zeros
        ellipse = GetEllipse(ring['Dset'],StaControls)
        Ring,delt = makeRing(ring['Dset'],ellipse,ring['pixLimit'],ring['cutoff'],scalex,scaley,Image)
        Ring = np.array(Ring).T
        ring['ImxyObs'] = np.array(Ring[:2])      #need to apply masks to this to eliminate bad points
        Ring[:2] = GetTthAzm(Ring[0],Ring[1],StaControls)       #convert x,y to tth,azm
        Ring[0] /= 2.                                           #convert to theta
        if len(rings):
            rings = np.concatenate((rings,Ring),axis=1)
        else:
            rings = np.array(Ring)
    E = StrSta['strain']
    p0 = [E[0][0],E[1][1],E[2][2],E[0][1],E[0][2],E[1][2]]
    E = FitStrain(rings,p0,wave,phi)
    StrSta['strain'] = E

def calcFij(omg,phi,azm,th):
    '''Does something...

    Uses parameters as defined by Bob He & Kingsley Smith, Adv. in X-Ray Anal. 41, 501 (1997)

    :param omg: his omega = sample omega rotation; 0 when incident beam || sample surface, 90 when perp. to sample surface
    :param phi: his phi = sample phi rotation; usually = 0, axis rotates with omg.
    :param azm: his chi = azimuth around incident beam
    :param th:  his theta = theta
    '''
    a = npsind(th)*npcosd(omg)+npsind(azm)*npcosd(th)*npsind(omg)
    b = -npcosd(azm)*npcosd(th)
    c = npsind(th)*npsind(omg)-npsind(azm)*npcosd(th)*npcosd(omg)
    d = a*npcosd(phi)-b*npsind(phi)
    e = a*npsind(phi)+b*npcosd(phi)
    Fij = np.array([
        [d**2,d*e,c*d],
        [d*e,e**2,c*e],
        [c*d,c*e,c**2]])
    return -Fij*nptand(th)

def FitStrain(rings,p0,wave,phi):
    'Needs a doc string'
    def StrainPrint(ValSig):
        print 'Strain tensor:'
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for name,fmt,value,sig in ValSig:
            ptlbls += "%s" % (name.rjust(12))
            ptstr += fmt % (value)
            if sig:
                sigstr += fmt % (sig)
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr
        
    def strainCalc(p,xyd,wave,phi):
#        E = np.array([[p[0],p[3],p[4]],[p[3],p[1],p[5]],[p[4],p[5],p[2]]])
        E = np.array([[p[0],0,0],[0,p[1],0],[0,0,0]])
        th,azm,dsp = xyd
        th0 = npasind(wave/(2.*dsp))
        dth = 180.*np.sum(E*calcFij(phi,0.,azm,th).T)/(np.pi*1.e6) #in degrees & microstrain units
        th0 += dth
        return (th-th0)**2
        
    names = ['e11','e22','e33','e12','e13','e23']
    fmt = ['%12.2f','%12.2f','%12.2f','%12.2f','%12.2f','%12.5f']
    p1 = [p0[0],p0[1]]   
    result = leastsq(strainCalc,p1,args=(rings,wave,phi),full_output=True)
    vals = list(result[0])
    chi = np.sqrt(np.sum(strainCalc(result[0],rings,wave,phi)**2))
    sig = list(chi*np.sqrt(np.diag(result[1])))
    ValSig = zip(names,fmt,vals,sig)
    StrainPrint(ValSig)
#    return np.array([[vals[0],vals[3],vals[4]],[vals[3],vals[1],vals[5]],[vals[4],vals[5],vals[2]]])
    return np.array([[vals[0],0,0],[0,vals[1],0],[0,0,0]])
    
        
