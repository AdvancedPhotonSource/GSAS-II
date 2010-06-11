#GSASII image calculations: ellipse fitting & image integration        
import math
import wx
import time
import numpy as np
import numpy.linalg as nl
import GSASIIplot as G2plt
import GSASIIlattice as G2lat

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
        
def makeMat(Angle,Axis):
    #Make rotation matrix from Angle in degrees,Axis =0 for rotation about x, =1 for about y, etc.
    cs = cosd(Angle)
    ss = sind(Angle)
    M = np.array(([1.,0.,0.],[0.,cs,-ss],[0.,ss,cs]),dtype=np.float32)
    return np.roll(np.roll(M,Axis,axis=0),Axis,axis=1)
                    
def FitRing(ring):
    Err,parms = FitCircle(ring)
    Err /= len(ring)
#    print 'circle error:','%8f'%(Err)
    if Err > 20000.:
        eparms = FitEllipse(ring)
        if eparms:
            parms = eparms
    return parms
        
def FitCircle(ring):
    
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
    result = nl.lstsq(M.T,B)
    return makeParmsEllipse(result[0])
    
def FitDetector(rings,p0,wave):
    from scipy.optimize import leastsq
    def ellipseCalc(B,xyd,wave):
        x = xyd[0]
        y = xyd[1]
        dsp = xyd[2]
        dist,x0,y0,phi,tilt = B
        tth = 2.0*npasind(wave/(2.*dsp))
        ttth = nptand(tth)
        radius = dist*ttth
        stth = npsind(tth)
        cosb = npcosd(tilt)
        R1 = dist*stth*npcosd(tth)*cosb/(cosb**2-stth**2)
        R0 = np.sqrt(R1*radius*cosb)
        zdis = R1*ttth*nptand(tilt)
        X = x-x0+zdis*npsind(phi)
        Y = y-y0-zdis*npcosd(phi)
        XR = X*npcosd(phi)-Y*npsind(phi)
        YR = X*npsind(phi)+Y*npcosd(phi)
        return (XR/R0)**2+(YR/R1)**2-1
    result = leastsq(ellipseCalc,p0,args=(rings.T,wave))
    return result[0]
            
def ImageLocalMax(image,w,Xpix,Ypix):
    w2 = w*2
    size = len(image)
    xpix = int(Xpix)            #get reference corner of pixel chosen
    ypix = int(Ypix)
    if (w < xpix < size-w) and (w < ypix < size-w) and image[ypix,xpix]:
        Z = image[ypix-w:ypix+w,xpix-w:xpix+w]
        Zmax = np.argmax(Z)
        Zmin = np.argmin(Z)
        xpix += Zmax%w2-w
        ypix += Zmax/w2-w
        return xpix,ypix,np.ravel(Z)[Zmax],np.ravel(Z)[Zmin]
    else:
        return 0,0,0,0
    
def makeRing(dsp,ellipse,pix,reject,scalex,scaley,image):
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    for a in range(-180,180,2):
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
            ring.append([X,Y,dsp])
    if len(ring) < 45:             #want more than 1/4 of a circle
        return []
    return ring
    
def makeIdealRing(ellipse):
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    for a in range(0,360,2):
        x = radii[0]*cosd(a)
        y = radii[1]*sind(a)
        X = (cphi*x-sphi*y+cent[0])
        Y = (sphi*x+cphi*y+cent[1])
        ring.append([X,Y])
    return ring
    
def calcDist(radii,tth):
    stth = sind(tth)
    ctth = cosd(tth)
    ttth = tand(tth)
    return math.sqrt(radii[0]**4/(ttth**2*((radii[0]*ctth)**2+(radii[1]*stth)**2)))
    
def calcZdisCosB(radius,tth,radii):
    cosB = sinb = radii[0]**2/(radius*radii[1])
    if cosB > 1.:
        return 0.,1.
    else:
        cosb = math.sqrt(1.-sinb**2)
        ttth = tand(tth)
        zdis = radii[1]*ttth*cosb/sinb
        return zdis,cosB
    
def GetEllipse(dsp,data):
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    radii = [0,0]
    tth = 2.0*asind(data['wavelength']/(2.*dsp))
    ttth = tand(tth)
    stth = sind(tth)
    ctth = cosd(tth)
    cosb = cosd(tilt)
    radius = dist*ttth
    radii[1] = dist*stth*ctth*cosb/(cosb**2-stth**2)
    if radii[1] > 0:
        radii[0] = math.sqrt(radii[1]*radius*cosb)
        zdis = radii[1]*ttth*tand(tilt)
        elcent = [cent[0]-zdis*sind(phi),cent[1]+zdis*cosd(phi)]
        return elcent,phi,radii
    else:
        return False
        
def GetDetectorXY(dsp,azm,data):
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
    tth = 2.0*asind(wave/(2.*dsp))
    ttth = tand(tth)
    radius = dist*ttth
    stth = sind(tth)
    cosb = cosd(tilt)
    R1 = dist*stth*cosd(tth)*cosb/(cosb**2-stth**2)
    R0 = math.sqrt(R1*radius*cosb)
    zdis = R1*ttth*tand(tilt)
    A = zdis*sind(phi)
    B = -zdis*cosd(phi)
    xy0 = [radius*cosd(azm),radius*sind(azm)]
    xy = fsolve(func,xy0,args=(azm,phi,R0,R1,A,B))+cent
    return xy
                    
def GetTthAzmDsp(x,y,data):
    wave = data['wavelength']
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    dx = np.array(x-cent[0],dtype=np.float32)
    dy = np.array(y-cent[1],dtype=np.float32)
    X = np.array(([dx,dy,np.zeros_like(dx)]),dtype=np.float32).T
    X = np.dot(X,makeMat(phi,2))
    Z = np.dot(X,makeMat(tilt,0)).T[2]
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z))
    dsp = wave/(2.*npsind(tth/2.))
    azm = npatan2d(dy,dx)
    return tth,azm,dsp
    
def GetTth(x,y,data):
    return GetTthAzmDsp(x,y,data)[0]
    
def GetTthAzm(x,y,data):
    return GetTthAzmDsp(x,y,data)[0:2]
    
def GetDsp(x,y,data):
    return GetTthAzmDsp(x,y,data)[2]
       
def ImageCompress(image,scale):
    if scale == 1:
        return image
    else:
        return image[::scale,::scale]
        
def ImageCalibrate(self,data):
    import copy
    import ImageCalibrants as calFile
    print 'image calibrate'
    ring = data['ring']
    pixelSize = data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    cutoff = data['cutoff']
    if len(ring) < 5:
        print 'not enough inner ring points for ellipse'
        return False
        
    #fit start points on inner ring
    data['ellipses'] = []
    outE = FitRing(ring)
    if outE:
        print 'start ellipse:',outE
        ellipse = outE
    else:
        return False
        
    #setup 180 points on that ring for "good" fit
    Ring = makeRing(1.0,ellipse,20,cutoff,scalex,scaley,self.ImageZ)
    if Ring:
        ellipse = FitRing(Ring)
        Ring = makeRing(1.0,ellipse,20,cutoff,scalex,scaley,self.ImageZ)    #do again
        ellipse = FitRing(Ring)
    else:
        print '1st ring not sufficiently complete to proceed'
        return False
    print 'inner ring:',ellipse
    data['center'] = copy.copy(ellipse[0])           #not right!! (but useful for now)
    data['ellipses'].append(ellipse[:]+('r',))
    G2plt.PlotImage(self)
    
    #setup for calibration
    data['rings'] = []
    data['ellipses'] = []
    if not data['calibrant']:
        print 'no calibration material selected'
        return True
        
    Bravais,cell,skip = calFile.Calibrants[data['calibrant']]
    A = G2lat.cell2A(cell)
    wave = data['wavelength']
    cent = data['center']
    pixLimit = data['pixLimit']
    elcent,phi,radii = ellipse
    HKL = G2lat.GenHBravais(0.5,Bravais,A)[skip:]
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
        radii[1] = dist*stth*ctth*cosB/(cosB**2-stth**2)
        radii[0] = math.sqrt(radii[1]*radius*cosB)
        zdis,cosB = calcZdisCosB(radius,tth,radii)
        zsinp = zdis*sind(phi)
        zcosp = zdis*cosd(phi)
        cent = data['center']
        elcent = [cent[0]+zsinp,cent[1]-zcosp]
        ratio = radii[1]/radii[0]
        Ring = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
            numZ = len(Ring)
            data['rings'].append(np.array(Ring))
            ellipse = FitRing(Ring)
            elcent,phi,radii = ellipse                
            if abs(phi) > 45. and phi < 0.:
                phi += 180.
            dist = calcDist(radii,tth)
            distR = 1.-dist/data['distance']
            if distR > 0.001:
                print 'Wavelength too large?'
            elif distR < -0.001:
                print 'Wavelength too small?'
            else:
                if abs((radii[1]/radii[0]-ratio)/ratio) > 0.01:
                    print 'Bad fit for ring # %i. Try reducing Pixel search range'%(i)
                    return False
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
            cent = data['center']
            print ('for ring # %2i dist %.3f rotate %6.2f tilt %6.2f Xcent %.3f Ycent %.3f Npts %d' 
                %(i,dist,phi,Tilt,cent[0],cent[1],numZ))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
            G2plt.PlotImage(self)
        else:
            break
    fullSize = len(self.ImageZ)/scalex
    if 2*radii[1] < .9*fullSize:
        print 'Are all usable rings (>25% visible) used? Try reducing Min ring I/Ib'
    if not Zsum:
        print 'Only one ring fitted. Check your wavelength.'
        return False
    cent = data['center'] = [xSum/Zsum,ySum/Zsum]
    wave = data['wavelength']
    dist = data['distance'] = distSum/Zsum
    
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
    
    tilt = data['tilt'] = Zsign*tiltSum/Zsum
    phi = data['rotation'] = phiSum/Zsum
    rings = np.concatenate((data['rings']),axis=0)
    p0 = [dist,cent[0],cent[1],phi,tilt]
    result = FitDetector(rings,p0,wave)
    data['distance'] = result[0]
    data['center'] = result[1:3]
    data['rotation'] = np.mod(result[3],360.0)
    data['tilt'] = result[4]
    N = len(data['ellipses'])
    data['ellipses'] = []           #clear away individual ellipse fits
    for H in HKL[:N]:
        ellipse = GetEllipse(H[3],data)
        data['ellipses'].append(copy.deepcopy(ellipse+('b',)))
    G2plt.PlotImage(self)        
    return True
    
def Make2ThetaAzimuthMap(data,imageN):
    #transforms 2D image from x,y space to 2-theta,azimuth space based on detector orientation
    pixelSize = data['pixelSize']
    scalex = pixelSize[0]/1000.
    scaley = pixelSize[1]/1000.
    tax,tay = np.mgrid[0.5:imageN+.5,0.5:imageN+.5]         #bin centers not corners
    tax = np.asfarray(tax*scalex,dtype=np.float32)
    tay = np.asfarray(tay*scaley,dtype=np.float32)
    return GetTthAzm(tay,tax,data)           #2-theta & azimuth arrays

def Fill2ThetaAzimuthMap(masks,TA,image):
    import numpy.ma as ma
    Zlim = masks['Thresholds'][1]
    imageN = len(image)
    TA = np.reshape(TA,(2,imageN,imageN))
    TA = np.dstack((ma.getdata(TA[1]),ma.getdata(TA[0])))    #azimuth, 2-theta
    tax,tay = np.dsplit(TA,2)    #azimuth, 2-theta
    taz = ma.masked_greater(ma.masked_less(image,Zlim[0]),Zlim[1]).flatten()
    tam = ma.getmask(taz)
    tax = ma.compressed(ma.array(tax.flatten(),mask=tam))
    tay = ma.compressed(ma.array(tay.flatten(),mask=tam))
    taz = ma.compressed(taz)
    del(tam)
    return tax,tay,taz
    
def Bin2ThetaAzimuthMap(data,tax,tay,taz):
    import numpy.ma as ma
    LUtth = data['IOtth']
    if data['fullIntegrate']:
        LRazm = [-180,180]
    else:
        LRazm = data['LRazimuth']
    numAzms = data['outAzimuths']
    numChans = data['outChannels']
    NST = np.histogram2d(tax,tay,normed=False,bins=(numAzms,numChans),range=[LRazm,LUtth])
    HST = np.histogram2d(tax,tay,normed=False,bins=(numAzms,numChans),range=[LRazm,LUtth],weights=taz)
    return NST,HST

def ImageIntegrate(self,data,masks):
    dlg = wx.ProgressDialog("Elapsed time","2D image integration",5,
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
    try:
        print 'Begin image integration'
        print 'Create 2-theta,azimuth map'
        t0 = time.time()
        dlg.Update(0)
        imageN = len(self.ImageZ)
        TA = Make2ThetaAzimuthMap(data,imageN)           #2-theta & azimuth arrays
        dlg.Update(1)
        print 'Fill map with 2-theta/azimuth values'
        tax,tay,taz = Fill2ThetaAzimuthMap(masks,TA,self.ImageZ)
        del TA
        dlg.Update(2)
        print 'Bin image by 2-theta/azimuth intervals'
        NST,HST = Bin2ThetaAzimuthMap(data,tax,tay,taz)
        del tax,tay,taz
        dlg.Update(3)
        print 'Form normalized 1D pattern(s)'
        self.Integrate = [HST[0]/NST[0],HST[1],HST[2]]
        del NST,HST
        dlg.Update(4)
        t1 = time.time()
        print 'Integration complete'
        print "Elapsed time:","%8.3f"%(t1-t0), "s"
    finally:
        dlg.Destroy()
