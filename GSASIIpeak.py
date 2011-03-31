#GSASII powder calculation module
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import sys
import math
import wx
import time
import numpy as np
import numpy.linalg as nl
import GSASIIpath
import pypowder as pyp              #assumes path has been amended to include correctr bin directory
import GSASIIplot as G2plt

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
npasind = lambda x: 180.*np.arcsin(x)/math.pi
npcosd = lambda x: np.cos(x*math.pi/180.)
nptand = lambda x: np.tan(x*math.pi/180.)
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

def factorize(num):
    ''' Provide prime number factors for integer num
    Returns dictionary of prime factors (keys) & power for each (data)
    '''
    factors = {}
    orig = num

    # we take advantage of the fact that (i +1)**2 = i**2 + 2*i +1
    i, sqi = 2, 4
    while sqi <= num:
        while not num%i:
            num /= i
            factors[i] = factors.get(i, 0) + 1

        sqi += 2*i + 1
        i += 1

    if num != 1 and num != orig:
        factors[num] = factors.get(num, 0) + 1

    if factors:
        return factors
    else:
        return {num:1}          #a prime number!
            
def makeFFTsizeList(nmin=1,nmax=1023,thresh=15):
    ''' Provide list of optimal data sizes for FFT calculations
    Input:
        nmin: minimum data size >= 1
        nmax: maximum data size > nmin
        thresh: maximum prime factor allowed
    Returns:
        list of data sizes where the maximum prime factor is < thresh
    ''' 
    plist = []
    nmin = max(1,nmin)
    nmax = max(nmin+1,nmax)
    for p in range(nmin,nmax):
        if max(factorize(p).keys()) < thresh:
            plist.append(p)
    return plist

def Transmission(Geometry,Abs,Diam):
#Calculate sample transmission
#   Geometry: one of 'Cylinder','Bragg-Brentano','Tilting Flat Plate in transmission','Fixed flat plate'
#   Abs: absorption coeff in cm-1
#   Diam: sample thickness/diameter in mm
    if 'Cylinder' in Geometry:      #Lobanov & Alte da Veiga for 2-theta = 0; beam fully illuminates sample
        MuR = Abs*Diam/5.0
        if MuR <= 3.0:
            T0 = 16/(3.*math.pi)
            T1 = -0.045780
            T2 = -0.02489
            T3 = 0.003045
            T = -T0*MuR-T1*MuR**2-T2*MuR**3-T3*MuR**4
            if T < -20.:
                return 2.06e-9
            else:
                return math.exp(T)
        else:
            T1 = 1.433902
            T2 = 0.013869+0.337894
            T3 = 1.933433+1.163198
            T4 = 0.044365-0.04259
            T = (T1-T4)/(1.0+T2*(MuR-3.0))**T3+T4
            return T/100.

def Absorb(Geometry,Abs,Diam,Tth,Phi=0,Psi=0):
#Calculate sample absorption
#   Geometry: one of 'Cylinder','Bragg-Brentano','Tilting Flat Plate in transmission','Fixed flat plate'
#   Abs: absorption coeff in cm-1
#   Diam: sample thickness/diameter in mm
#   Tth: 2-theta scattering angle - can be numpy array
#   Phi: flat plate tilt angle - future
#   Psi: flat plate tilt axis - future
    MuR = Abs*Diam/5.0
    Sth2 = npsind(Tth/2.0)**2
    Cth2 = 1.-Sth2
    if 'Cylinder' in Geometry:      #Lobanov & Alte da Veiga for 2-theta = 0; beam fully illuminates sample
        if MuR < 3.0:
            T0 = 16.0/(3*np.pi)
            T1 = (25.99978-0.01911*Sth2**0.25)*np.exp(-0.024551*Sth2)+ \
                0.109561*np.sqrt(Sth2)-26.04556
            T2 = -0.02489-0.39499*Sth2+1.219077*Sth2**1.5- \
                1.31268*Sth2**2+0.871081*Sth2**2.5-0.2327*Sth2**3
            T3 = 0.003045+0.018167*Sth2-0.03305*Sth2**2
            Trns = -T0*MuR-T1*MuR**2-T2*MuR**3-T3*MuR**4
            return np.exp(Trns)
        else:
            T1 = 1.433902+11.07504*Sth2-8.77629*Sth2*Sth2+ \
                10.02088*Sth2**3-3.36778*Sth2**4
            T2 = (0.013869-0.01249*Sth2)*np.exp(3.27094*Sth2)+ \
                (0.337894+13.77317*Sth2)/(1.0+11.53544*Sth2)**1.555039
            T3 = 1.933433/(1.0+23.12967*Sth2)**1.686715- \
                0.13576*np.sqrt(Sth2)+1.163198
            T4 = 0.044365-0.04259/(1.0+0.41051*Sth2)**148.4202
            Trns = (T1-T4)/(1.0+T2*(MuR-3.0))**T3+T4
            return Trns/100.
    elif 'Bragg' in Geometry:
        return 1.0
    elif 'Fixed' in Geometry: #assumes sample plane is perpendicular to incident beam
        # and only defined for 2theta < 90
        T1 = np.exp(-MuR)
        T2 = np.exp(-MuR/(1.-2.*Sth2))
        Tb = -2.*Abs*Sth2
        return (T1-T2)/Tb
    elif 'Tilting' in Geometry: #assumes symmetric tilt so sample plane is parallel to diffraction vector
        cth = npcosd(Tth/2.0)
        return (Diam/cth)*np.exp(-MuR/cth)
        
def Polarization(Pola,Azm,Tth):
#   Calculate x-ray polarization correction
#   Pola: polarization coefficient e.g 1.0 fully polarized, 0.5 unpolarized
#   Azm: azimuthal angle e.g. 0.0 in plane of polarization(?)
#   Tth: 2-theta scattering angle - can be numpy array
    pass
    
        

def ValEsd(value,esd=0,nTZ=False):                  #NOT complete - don't use
    # returns value(esd) string; nTZ=True for no trailing zeros
    # use esd < 0 for level of precision shown e.g. esd=-0.01 gives 2 places beyond decimal
    #get the 2 significant digits in the esd 
    edig = lambda esd: int(round(10**(math.log10(esd) % 1+1)))
    #get the number of digits to represent them 
    epl = lambda esd: 2+int(1.545-math.log10(10*edig(esd)))
    
    mdec = lambda esd: -int(math.log10(abs(esd)))
    ndec = lambda esd: int(1.545-math.log10(abs(esd)))
    if esd > 0:
        fmt = '"%.'+str(ndec(esd))+'f(%d)"'
        print fmt,ndec(esd),esd*10**(mdec(esd)+1)
        return fmt%(value,int(esd*10**(mdec(esd)+1)))
    elif esd < 0:
         return str(round(value,mdec(esd)))
    else:
        text = "%F"%(value)
        if nTZ:
            return text.rstrip('0')
        else:
            return text

        
#GSASII peak fitting routine: Thompson, Cox & Hastings; Finger, Cox & Jephcoat model        

def DoPeakFit(peaks,background,limits,inst,data):
    
    def backgroundPrint(background,sigback):
        if background[1]:
            print 'Background coefficients for',background[0],'function'
            ptfmt = "%12.5f"
            ptstr =  'values:'
            sigstr = 'esds  :'
            for i,back in enumerate(background[3:]):
                ptstr += ptfmt % (back)
                sigstr += ptfmt % (sigback[i+3])
            print ptstr
            print sigstr
        else:
            print 'Background not refined'
    
    def instPrint(instVal,siginst,insLabels):
        print 'Instrument Parameters:'
        ptfmt = "%12.6f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for i,value in enumerate(instVal):
            ptlbls += "%s" % (insLabels[i].center(12))
            ptstr += ptfmt % (value)
            if siginst[i]:
                sigstr += ptfmt % (siginst[i])
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr
    
    def peaksPrint(peaks,sigpeaks):
        print 'Peak list='

    begin = time.time()
    Np = len(peaks[0])
    DataType = inst[1][0]
    instVal = inst[1][1:]
    Insref = inst[2][1:]
    insLabels = inst[3][1:]
    Ka2 = False
    Ioff = 3
    if len(instVal) == 12:
        lamratio = instVal[1]/instVal[0]
        Ka2 = True
        Ioff = 5
    insref = Insref[len(Insref)-7:-1]               #just U,V,W,X,Y,SH/L
    for peak in peaks:
        dip = []
        dip.append(tand(peak[0]/2.0)**2)
        dip.append(tand(peak[0]/2.0))
        dip.append(1.0/cosd(peak[0]/2.0))
        dip.append(tand(peak[0]/2.0))
        peak.append(dip)
    B = background[2]
    bcof = background[3:3+B]
    Bv = 0
    if background[1]:
        Bv = B
    x,y,w,yc,yb,yd = data               #these are numpy arrays!
    V = []
    A = []
    swobs = 0.0
    smin = 0.0
    first = True
    GoOn = True
    Go = True
    dlg = wx.ProgressDialog("Elapsed time","Fitting peaks to pattern",len(x), \
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
    screenSize = wx.DisplaySize()
    Size = dlg.GetSize()
    dlg.SetPosition(wx.Point(screenSize[0]-Size[0]-300,0))
    try:
        i = 0
        for xi in x :
            Go = dlg.Update(i)[0]
            if GoOn:
                GoOn = Go
            if limits[0] <= xi <= limits[1]:
                yb[i] = 0.0
                dp = []
                for j in range(B):
                    t = (xi-limits[0])**j
                    yb[i] += t*bcof[j]
                    if background[1]:
                        dp.append(t)
                yc[i] = yb[i]
                Iv = 0
                for j in range(6):
                    if insref[j]:
                        dp.append(0.0)
                        Iv += 1
                for peak in peaks:
                    dip = peak[-1]
                    f = pyp.pypsvfcj(peak[2],xi-peak[0],peak[0],peak[4],peak[6],instVal[-2],0.0)
                    yc[i] += f[0]*peak[2]
                    if f[0] > 0.0:
                        j = 0
                        if insref[0]:              #U
                            dp[Bv+j] += f[3]*dip[0]
                            j += 1
                        if insref[1]:              #V
                            dp[Bv+j] += f[3]*dip[1]
                            j += 1
                        if insref[2]:              #W
                            dp[Bv+j] += f[3]
                            j += 1
                        if insref[3]:              #X
                            dp[Bv+j] += f[4]*dip[2]
                            j += 1
                        if insref[4]:              #Y
                            dp[Bv+j] += f[4]*dip[3]
                            j += 1
                        if insref[5]:              #SH/L
                            dp[Bv+j] += f[5]
                    if Ka2:
                       pos2 = 2.0*asind(lamratio*sind(peak[0]/2.0))
                       f2 = pyp.pypsvfcj(peak[2],xi-pos2,peak[0],peak[4],peak[6],instVal[-2],0.0)
                       yc[i] += f2[0]*peak[2]*instVal[3]
                       if f[0] > 0.0:
                           j = 0
                           if insref[0]:              #U
                               dp[Bv+j] += f2[3]*dip[0]*instVal[3]
                               j += 1
                           if insref[1]:              #V
                               dp[Bv+j] += f2[3]*dip[1]*instVal[3]
                               j += 1
                           if insref[2]:              #W
                               dp[Bv+j] += f2[3]*instVal[3]
                               j += 1
                           if insref[3]:              #X
                               dp[Bv+j] += f2[4]*dip[2]*instVal[3]
                               j += 1
                           if insref[4]:              #Y
                               dp[Bv+j] += f2[4]*dip[3]*instVal[3]
                               j += 1
                           if insref[5]:              #SH/L
                               dp[Bv+j] += f2[5]*instVal[3]                       
                    for j in range(0,Np,2):
                        if peak[j+1]: dp.append(f[j/2+1])
                yd[i] = y[i]-yc[i]
                swobs += w[i]*y[i]**2
                t2 = w[i]*yd[i]
                smin += t2*yd[i]
                if first:
                    first = False
                    M = len(dp)
                    A = np.zeros(shape=(M,M))
                    V = np.zeros(shape=(M))
                A,V = pyp.buildmv(t2,w[i],M,dp,A,V)
            i += 1
    finally:
        dlg.Destroy()
    Rwp = smin/swobs
    Rwp = math.sqrt(Rwp)*100.0
    norm = np.diag(A)
    for i,elm in enumerate(norm):
        if elm <= 0.0:
            print norm
            return False,0,0,0,False
    for i in xrange(len(V)):
        norm[i] = 1.0/math.sqrt(norm[i])
        V[i] *= norm[i]
        a = A[i]
        for j in xrange(len(V)):
            a[j] *= norm[i]
    A = np.transpose(A)
    for i in xrange(len(V)):
        a = A[i]
        for j in xrange(len(V)):
            a[j] *= norm[i]
    b = nl.solve(A,V)
    A = nl.inv(A)
    sig = np.diag(A)
    for i in xrange(len(V)):
        b[i] *= norm[i]
        sig[i] *= norm[i]*norm[i]
        sig[i] = math.sqrt(abs(sig[i]))
    sigback = [0,0,0]
    for j in range(Bv):
        background[j+3] += b[j]
        sigback.append(sig[j])
    backgroundPrint(background,sigback)
    k = 0
    delt = []
    if Ka2:
        siginst = [0,0,0,0,0]
    else:
        siginst = [0,0,0]
    for j in range(6):
        if insref[j]:
            instVal[j+Ioff] += b[Bv+k]
            siginst.append(sig[Bv+k])
            delt.append(b[Bv+k])
            k += 1
        else:
            delt.append(0.0)
            siginst.append(0.0)
    delt.append(0.0)                    #dummies for azm
    siginst.append(0.0)
    instPrint(instVal,siginst,insLabels)
    inst[1] = [DataType,]
    for val in instVal:
        inst[1].append(val)
    B = Bv+Iv
    for peak in peaks:
        del peak[-1]                        # remove dip from end
        delsig = delt[0]*tand(peak[0]/2.0)**2+delt[1]*tand(peak[0]/2.0)+delt[2]
        delgam = delt[3]/cosd(peak[0]/2.0)+delt[4]*tand(peak[0]/2.0)
        for j in range(0,len(peak[:-1]),2):
            if peak[j+1]: 
                peak[j] += b[B]
                B += 1
        peak[4] += delsig
        if peak[4] < 0.0:
            print 'ERROR - negative sigma'
            return False,0,0,0,False            
        peak[6] += delgam
        if peak[6] < 0.0:
            print 'ERROR - negative gamma'
            return False,0,0,0,False
    runtime = time.time()-begin    
    data = [x,y,w,yc,yb,yd]
    return True,smin,Rwp,runtime,GoOn

def ComputePDF(data,xydata):
    for key in data:
        print key,data[key]
    #subtract backgrounds - if any
    xydata['Sample corrected'] = xydata['Sample']
    if 'Sample Bkg.' in xydata:
        xydata['Sample corrected'][1][1] -= (xydata['Sample Bkg.'][1][1]+
            data['Sample Bkg.']['Add'])*data['Sample Bkg.']['Mult']
    if 'Container' in xydata:    
        xydata['Sample corrected'][1][1] -= (xydata['Container'][1][1]+
            data['Container']['Add'])*data['Container']['Mult']
    if 'Container Bkg.' in xydata:
        xydata['Sample corrected'][1][1] += (xydata['Container Bkg.'][1][1]+
            data['Container Bkg.']['Add'])*data['Container Bkg.']['Mult']
    
           
        
    return xydata['Sample corrected'],[]
        
