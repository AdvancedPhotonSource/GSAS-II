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
actually not easy in this case wx.ProgressDialog needs #of blocks to process when started
not known in G2imageGUI.
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
npacosd = lambda x: 180.*np.arccos(x)/np.pi
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
    return parms,err
        
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
            print B
            return 0
        elif det == 0.:
            print 'parabola!'
            print B
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
    
def FitDetector(rings,varyList,parmDict,Print=True):
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
        x,y,dsp = xyd
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
        stth = npsind(tth)
        cosb = npcosd(tilt)
        tanb = nptand(tilt)        
        tbm = nptand((tth-tilt)/2.)
        tbp = nptand((tth+tilt)/2.)
        sinb = npsind(tilt)
        d = dist+dxy
        fplus = d*tanb*stth/(cosb+stth)
        fminus = d*tanb*stth/(cosb-stth)
        vplus = d*(tanb+(1+tbm)/(1-tbm))*stth/(cosb+stth)
        vminus = d*(tanb+(1-tbp)/(1+tbp))*stth/(cosb-stth)
        R0 = np.sqrt((vplus+vminus)**2-(fplus+fminus)**2)/2.      #+minor axis
        R1 = (vplus+vminus)/2.                                    #major axis
        zdis = (fplus-fminus)/2.
        X = x-x0-zdis*npsind(phi)
        Y = y-y0-zdis*npcosd(phi)
        XR = X*npcosd(phi)+Y*npsind(phi)
        YR = X*npsind(phi)+Y*npcosd(phi)
        return (XR/R0)**2+(YR/R1)**2-1
        
    names = ['dist','det-X','det-Y','tilt','phi','dep','wave']
    fmt = ['%12.2f','%12.2f','%12.2f','%12.2f','%12.2f','%12.3f','%12.5f']
    p0 = [parmDict[key] for key in varyList]
    result = leastsq(ellipseCalcD,p0,args=(rings.T,varyList,parmDict),full_output=True,ftol=1.e-8)
    chisq = np.sum(result[2]['fvec']**2)
    parmDict.update(zip(varyList,result[0]))
    vals = list(result[0])
    chi = np.sqrt(np.sum(ellipseCalcD(result[0],rings.T,varyList,parmDict)**2))
    sig = list(chi*np.sqrt(np.diag(result[1])))
    sigList = np.zeros(7)
    for i,name in enumerate(names):
        if name in varyList:
            sigList[i] = sig[varyList.index(name)]
    ValSig = zip(names,fmt,vals,sig)
    if Print:
        CalibPrint(ValSig)
    return chi,chisq
                    
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
    def ellipseC():
        'compute estimate of ellipse circumference'
        if radii[0] < 0:        #hyperbola
            theta = npacosd(1./np.sqrt(1.+(radii[0]/radii[1])**2))
            print theta
            return 0
        apb = radii[1]+radii[0]
        amb = radii[1]-radii[0]
        return np.pi*apb*(1+3*(amb/apb)**2/(10+np.sqrt(4-3*(amb/apb)**2)))
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    sumI = 0
    amin = 0
    amax = 360
    for a in range(0,int(ellipseC()),1):
        x = radii[0]*cosd(a)
        y = radii[1]*sind(a)
        X = (cphi*x-sphi*y+cent[0])*scalex      #convert mm to pixels
        Y = (sphi*x+cphi*y+cent[1])*scaley
        X,Y,I,J = ImageLocalMax(image,pix,X,Y)
        if I and J and I/J > reject:
            sumI += I/J
            X += .5                             #set to center of pixel
            Y += .5
            X /= scalex                         #convert to mm
            Y /= scaley
            amin = min(amin,a)
            amax = max(amax,a)
            if [X,Y,dsp] not in ring:
                ring.append([X,Y,dsp])
    delt = amax-amin
    if len(ring) < 10:             #want more than 10 deg
        return [],(delt > 90),0
    return ring,(delt > 90),sumI/len(ring)
    
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
    slr = radii[0]**2/radii[1]
    rat2 = (radii[0]/radii[1])**2
    if radii[0] > 0.:       #ellipse
        ecc = np.sqrt(max(1.-rat2,0.))
    else:                   #hyperbola
        ecc = np.sqrt(1.+rat2)
    rad = slr/(1.+ecc*npcosd(a))
    xy = np.array([rad*npsind(a)+cent[0],rad*npcosd(a)+cent[1]])
    return xy
                    
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
        
def GetEllipse2(tth,dxy,dist,cent,tilt,phi):
    '''uses Dandelin spheres to find ellipse or hyperbola parameters from detector geometry
    on output
    radii[0] (b-minor axis) set < 0. for hyperbola
    
    '''
    radii = [0,0]
    ttth = tand(tth)
    stth = sind(tth)
    ctth = cosd(tth)
    cosb = cosd(tilt)
    tanb = tand(tilt)
    tbm = tand((tth-tilt)/2.)
    tbp = tand((tth+tilt)/2.)
    sinb = sind(tilt)
    d = dist+dxy
    if tth+tilt < 90.:      #ellipse
        fplus = d*tanb*stth/(cosb+stth)
        fminus = d*tanb*stth/(cosb-stth)
        vplus = d*(tanb+(1+tbm)/(1-tbm))*stth/(cosb+stth)
        vminus = d*(tanb+(1-tbp)/(1+tbp))*stth/(cosb-stth)
        radii[0] = np.sqrt((vplus+vminus)**2-(fplus+fminus)**2)/2.      #+minor axis
        radii[1] = (vplus+vminus)/2.                                    #major axis
        zdis = (fplus-fminus)/2.
    else:   #hyperbola!
        f = d*tanb*stth/(cosb+stth)
        v = d*(tanb+tand(tth-tilt))
        delt = d*stth*(1+stth*cosb)/(sinb*cosb*(stth+cosb))
        eps = (v-f)/(delt-v)
        radii[0] = -eps*(delt-f)/np.sqrt(eps**2-1.)                     #-minor axis
        radii[1] = eps*(delt-f)/(eps**2-1.)                             #major axis
        zdis = f+radii[1]*eps
    elcent = [cent[0]+zdis*sind(phi),cent[1]+zdis*cosd(phi)]
    return elcent,phi,radii
    
def GetEllipse(dsp,data):
    '''uses Dandelin spheres to find ellipse or hyperbola parameters from detector geometry
    as given in image controls dictionary (data) and a d-spacing (dsp)
    '''
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    dep = data['DetDepth']
    tth = 2.0*asind(data['wavelength']/(2.*dsp))
    dxy = peneCorr(tth,dep)
    dist = data['distance']
    return GetEllipse2(tth,dxy,dist,cent,tilt,phi)
        
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
    phi = data['rotation']
    wave = data['wavelength']
    tilt = data['tilt']
    dist = data['distance']/cosd(tilt)
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
    'Needs a doc string - checked OK for ellipses dont know about hyperbola'
    wave = data['wavelength']
    cent = data['center']
    tilt = data['tilt']
    dist = data['distance']/cosd(tilt)
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
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z+dxy)) #depth corr not correct for tilted detector
    dsp = wave/(2.*npsind(tth/2.))
    azm = (npatan2d(dy,dx)+azmthoff+720.)%360.
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
    Not currently used but might be useful in future?
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
        hkl = G2lat.GenHBravais(dmin,bravais,A)
        HKL += hkl
    HKL = G2lat.sortHKLd(HKL,True,False)
    varyList = ['dist','det-X','det-Y','tilt','phi']
    parmDict = {'dist':data['distance'],'det-X':data['center'][0],'det-Y':data['center'][1],
        'tilt':data['tilt'],'phi':data['rotation'],'wave':data['wavelength'],'dep':data['DetDepth']}
    Found = False
    for H in HKL: 
        dsp = H[3]
        ellipse = GetEllipse(dsp,data)
        Ring,delt,sumI = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
            data['rings'].append(np.array(Ring))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
            Found = True
        elif not Found:         #skipping inner rings, keep looking until ring found 
            continue
        else:                   #no more rings beyond edge of detector
            break
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
    outE,err = FitRing(ring,True)
    fmt  = '%s X: %.3f, Y: %.3f, phi: %.3f, R1: %.3f, R2: %.3f'
    fmt2 = '%s X: %.3f, Y: %.3f, phi: %.3f, R1: %.3f, R2: %.3f, chi^2: %.3f, N: %d'
    if outE:
        print fmt%('start ellipse: ',outE[0][0],outE[0][1],outE[1],outE[2][0],outE[2][1])
        ellipse = outE
    else:
        return False
        
    #setup 360 points on that ring for "good" fit
    Ring,delt,sumI = makeRing(1.0,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
    if Ring:
        ellipse,err = FitRing(Ring,delt)
        Ring,delt,sumI = makeRing(1.0,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)    #do again
        ellipse,err = FitRing(Ring,delt)
    else:
        print '1st ring not sufficiently complete to proceed'
        return False
    print fmt2%('inner ring:    ',ellipse[0][0],ellipse[0][1],ellipse[1],ellipse[2][0],ellipse[2][1],err,len(Ring))     #cent,phi,radii
    data['ellipses'].append(ellipse[:]+('r',))
    data['rings'].append(np.array(Ring))
    G2plt.PlotImage(self,newImage=True)
    
    #setup for calibration
    data['rings'] = []
    data['ellipses'] = []
    if not data['calibrant']  or 'None' in data['calibrant']:
        print 'no calibration material selected'
        return True
    
    skip = data['calibskip']
    dmin = data['calibdmin']
#generate reflection set
    Bravais,Cells = calFile.Calibrants[data['calibrant']][:2]
    HKL = []
    for bravais,cell in zip(Bravais,Cells):
        A = G2lat.cell2A(cell)
        hkl = G2lat.GenHBravais(dmin,bravais,A)[skip:]
        HKL += hkl
    HKL = G2lat.sortHKLd(HKL,True,False)
    wave = data['wavelength']
#set up 1st ring
    elcent,phi,radii = ellipse              #from fit of 1st ring
    dsp = HKL[0][3]
    tth = 2.0*asind(wave/(2.*dsp))
    ttth = nptand(tth)
    stth = npsind(tth)
    ctth = npcosd(tth)
#1st estimate of tilt; assume ellipse - don't know sign though
    tilt = npasind(np.sqrt(max(0.,1.-(radii[0]/radii[1])**2))*ctth)
#1st estimate of dist: sample to detector normal to plane
    data['distance'] = dist = radii[0]**2/(ttth*radii[1])
#ellipse to cone axis (x-ray beam); 2 choices depending on sign of tilt
    zdisp = radii[1]*ttth*tand(tilt)
    zdism = radii[1]*ttth*tand(-tilt)
#cone axis position; 2 choices. Which is right?
    centp = [elcent[0]+zdisp*sind(phi),elcent[1]+zdisp*cosd(phi)]
    centm = [elcent[0]+zdism*sind(phi),elcent[1]+zdism*cosd(phi)]
#check get same ellipse parms either way
#now do next ring; estimate either way & check sum of Imax/Imin in 3x3 block around each point
    dsp = HKL[1][3]
    tth = 2.0*asind(wave/(2.*dsp))
    ellipsep = GetEllipse2(tth,0.,dist,centp,tilt,phi)
    Ringp,delt,sumIp = makeRing(dsp,ellipsep,2,cutoff,scalex,scaley,self.ImageZ)
    outEp,errp = FitRing(Ringp,True)
    print '+',ellipsep,errp
    ellipsem = GetEllipse2(tth,0.,dist,centm,-tilt,phi)
    Ringm,delt,sumIm = makeRing(dsp,ellipsem,2,cutoff,scalex,scaley,self.ImageZ)
    outEm,errm = FitRing(Ringm,True)
    print '-',ellipsem,errm
    if errp < errm:
        data['tilt'] = tilt
        data['center'] = centp
        negTilt = 1
    else:
        data['tilt'] = -tilt
        data['center'] = centm
        negTilt = -1
    data['rotation'] = phi
    parmDict = {'dist':data['distance'],'det-X':data['center'][0],'det-Y':data['center'][1],
        'tilt':data['tilt'],'phi':data['rotation'],'wave':data['wavelength'],'dep':data['DetDepth']}
    varyList = ['dist','det-X','det-Y','tilt','phi']
    if data['DetDepthRef']:
        varyList.append('dep')
    for i,H in enumerate(HKL):
        dsp = H[3]
        tth = 2.0*asind(wave/(2.*dsp))
        print 'HKLD:',H[:4],'2-theta: %.4f'%(tth)
        elcent,phi,radii = ellipse = GetEllipse(dsp,data)
        data['ellipses'].append(copy.deepcopy(ellipse+('g',)))
        print fmt%('predicted ellipse:',elcent[0],elcent[1],phi,radii[0],radii[1])
        Ring,delt,sumI = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
            data['rings'].append(np.array(Ring))
            rings = np.concatenate((data['rings']),axis=0)
            if i:
                chi,chisq = FitDetector(rings,varyList,parmDict,False)
                data['distance'] = parmDict['dist']
                data['center'] = [parmDict['det-X'],parmDict['det-Y']]
                data['rotation'] = np.mod(parmDict['phi'],360.0)
                data['tilt'] = parmDict['tilt']
                data['DetDepth'] = parmDict['dep']
                elcent,phi,radii = ellipse = GetEllipse(dsp,data)
                print fmt2%('fitted ellipse:   ',elcent[0],elcent[1],phi,radii[0],radii[1],chisq,len(rings))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
        else:
            print 'insufficient number of points in this ellipse to fit'
            break
    G2plt.PlotImage(self,newImage=True)
    fullSize = len(self.ImageZ)/scalex
    if 2*radii[1] < .9*fullSize:
        print 'Are all usable rings (>25% visible) used? Try reducing Min ring I/Ib'
    N = len(data['ellipses'])
    if N > 2:
        FitDetector(rings,varyList,parmDict)
        data['distance'] = parmDict['dist']
        data['center'] = [parmDict['det-X'],parmDict['det-Y']]
        data['rotation'] = np.mod(parmDict['phi'],360.0)
        data['tilt'] = parmDict['tilt']
        data['DetDepth'] = parmDict['dep']
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
    frame = masks['Frames']
    tam = ma.make_mask_none((nI,nJ))
    if frame:
        tamp = ma.make_mask_none((1024*1024))
        tamp = ma.make_mask(pm.polymask(nI*nJ,tax.flatten(),
            tay.flatten(),len(frame),frame,tamp)[:nI*nJ])-True  #switch to exclude around frame
        tam = ma.mask_or(tam.flatten(),tamp)
    polygons = masks['Polygons']
    for polygon in polygons:
        if polygon:
            tamp = ma.make_mask_none((1024*1024))
            tamp = ma.make_mask(pm.polymask(nI*nJ,tax.flatten(),
                tay.flatten(),len(polygon),polygon,tamp)[:nI*nJ])
            tam = ma.mask_or(tam.flatten(),tamp)
    if tam.shape: tam = np.reshape(tam,(nI,nJ))
    spots = masks['Points']
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
    blkSize = 128   #this seems to be optimal; will break in polymask if >1024
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
    nXBlks = (Nx-1)/blkSize+1
    nYBlks = (Ny-1)/blkSize+1
    Nup = nXBlks*nYBlks*3+3
    dlg = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
    try:
        t0 = time.time()
        Nup = 0
        dlg.Update(Nup)
        for iBlk in range(nYBlks):
            iBeg = iBlk*blkSize
            iFin = min(iBeg+blkSize,Ny)
            for jBlk in range(nXBlks):
                jBeg = jBlk*blkSize
                jFin = min(jBeg+blkSize,Nx)                
#                print 'Process map block:',iBlk,jBlk,' limits:',iBeg,iFin,jBeg,jFin
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
                if any([tax.shape[0],tay.shape[0],taz.shape[0]]):
                    NST,H0 = h2d.histogram2d(len(tax),tax,tay,taz,numAzms,numChans,LRazm,LUtth,Dazm,Dtth,NST,H0)
#                print 'block done'
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
        H0 /= npcosd(H2[:-1])**2
        if data['Oblique'][1]:
            H0 /= G2pwd.Oblique(data['Oblique'][0],H2[:-1])
        if 'SASD' in data['type'] and data['PolaVal'][1]:
            H0 /= np.array([G2pwd.Polarization(data['PolaVal'][0],H2[:-1],Azm=azm)[0] for azm in H1[:-1]+np.diff(H1)/2.])
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
        Ring,delt,sumI = makeRing(ring['Dset'],ellipse,ring['pixLimit'],ring['cutoff'],scalex,scaley,Image)
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
    
        
