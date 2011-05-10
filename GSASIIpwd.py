#GSASII powder calculation module
########### SVN repository information ###################
# $Date: 2011-04-20 13:09:53 -0500 (Wed, 20 Apr 2011) $
# $Author: vondreele $
# $Revision: 267 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIIpwd.py $
# $Id: GSASIIpwd.py 267 2011-04-20 18:09:53Z vondreele $
########### SVN repository information ###################
import sys
import math
import wx
import time
import numpy as np
import scipy as sp
import numpy.linalg as nl
import scipy.interpolate as si
import scipy.stats as st
import GSASIIpath
import pypowder as pyp              #assumes path has been amended to include correctr bin directory
import GSASIIplot as G2plt
import GSASIIlattice as G2lat
import GSASIIElem as G2elem
import GSASIIgrid as G2gd

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
npacosd = lambda x: 180.*np.arccos(x)/math.pi
nptand = lambda x: np.tan(x*math.pi/180.)
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
npT2stl = lambda tth, wave: 2.0*npsind(tth/2.0)/wave
npT2q = lambda tth,wave: 2.0*np.pi*npT2stl(tth,wave)
    
np.seterr(divide='ignore')      #this is for the FCJ functions

#Peak shape definitions
# Finger-Cox_Jephcoat D(2phi,2th) function for S/L = H/L

class fcjde_gen(st.rv_continuous):
    
    """
    Finger-Cox-Jephcoat D(2phi,2th) function for S/L = H/L
    Ref: J. Appl. Cryst. (1994) 27, 892-900.
    Parameters
    -----------------------------------------
    x: array like 2-theta positions
    t: 2-theta position of peak
    s: sum(S/L,H/L); S: sample height, H: detector opening, 
        L: sample to detector opening distance
    Result for fcj.pdf
    -----------------------------------------
    if x < t & s = S/L+H/L: 
        fcj.pdf = [1/sqrt({cos(x)**2/cos(t)**2}-1) - 1/s]/cos(x)    
    if x >= t:
        fcj.pdf = 0    
    """
    def _pdf(self,x,t,s):
        ax = npcosd(x)**2
        bx = npcosd(t)**2
        bx = np.where(ax>bx,bx,ax)
        fx = np.where(ax>bx,(1./np.sqrt((ax/bx)-1.)-1./s)/npcosd(x),0.0)
        return np.where(((x < t) & (fx > 0)),fx,0.0)
#    def _cdf(self, x):
#    def _ppf(self, q):
#    def _sf(self, x):
#    def _isf(self, q):
#    def _stats(self):
#    def _entropy(self):
        
fcjde = fcjde_gen(name='fcjde')
                
        
# Finger-Cox_Jephcoat D(2phi,2th) function for S/L != H/L

class fcjd_gen(st.rv_continuous):
    """
    Finger-Cox-Jephcoat D(2phi,2th) function for S/L != H/L
    Ref: J. Appl. Cryst. (1994) 27, 892-900.
    Parameters
    -----------------------------------------
    x: array like 2-theta positions
    t: 2-theta position of peak
    h: detector opening height/sample to detector opening distance
    s: sample height/sample to detector opening distance
    Result for fcj2.pdf
    -----------------------------------------
    infl = acos(cos(t)*sqrt((h-s)**2+1))
    if x < infl:    
        fcj.pdf = [1/sqrt({cos(x)**2/cos(t)**2}-1) - 1/shl]/cos(2phi)
    
    for 2phi < 2tth & shl = S/L+H/L
    
    fcj.pdf(x,tth,shl) = 0
    
    for 2phi >= 2th
    """
    def _pdf(self,x,t,h,s):
        a = npcosd(t)*(np.sqrt((h-s)**2+1.))
        infl = np.where((a <= 1.),npacosd(a),t)
        ax = npcosd(x)**2
        bx = npcosd(t)**2
        bx = np.where(ax>bx,bx,ax)
        H = np.where(ax>bx,np.sqrt((ax/bx)-1.),0.0)
        W1 = h+s-H
        W2 = np.where ((h > s),2.*s,2.*h)
        fx = 2.*h*np.sqrt((ax/bx)-1.)*npcosd(x)
        fx = np.where(fx>0.0,1./fx,0.0)
        fx = np.where((x < infl),fx*W1,fx*W2)
        return np.where((fx > 0.),fx,0.0)
#    def _cdf(self, x):
#    def _ppf(self, q):
#    def _sf(self, x):
#    def _isf(self, q):
#    def _stats(self):
#    def _entropy(self):
        
fcjd = fcjd_gen(name='fcjd')
                
# Finger-Cox_Jephcoat D(2phi,2th) function for S/L != H/L using sum & difference

class fcjdsd_gen(st.rv_continuous):
    """
    Finger-Cox-Jephcoat D(2phi,2th) function for S/L != H/L using sum & difference
    
    fcj.pdf(x,tth,shl) = [1/sqrt({cos(2phi)**2/cos(2th)**2}-1) - 1/shl]/cos(2phi)
    
    for 2phi < 2tth & shl = S/L+H/L
    
    fcj.pdf(x,tth,shl) = 0
    
    for 2phi >= 2th
    """
    def _argcheck(self,t,s,d):
        return (t > 0)&(s > 0)&(abs(d) < s)
    def _pdf(self,x,t,s,d):
        a = npcosd(t)*np.sqrt(d**2+1.)
        infl = np.where((a < 1.),npacosd(a),t)
        ax = npcosd(x)**2
        bx = npcosd(t)**2
        bx = np.where(ax>bx,bx,ax)
        H = np.where(ax>bx,np.sqrt((ax/bx)-1.),0.0)
        W1 = s-H
        W2 = np.where ((d > 0),s-d,s+d)
        fx = np.where(ax>bx,1./((s+d)*np.sqrt((ax/bx)-1.)*npcosd(x)),0.0)
        fx = np.where((x < infl),fx*W1,fx*W2)
        return np.where((fx > 0.),fx,0.0)
#    def _cdf(self, x):
#    def _ppf(self, q):
#    def _sf(self, x):
#    def _isf(self, q):
#    def _stats(self):
#    def _entropy(self):
        
fcjdsd = fcjdsd_gen(name='fcjdsd')
                
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
#   Geometry: one of 'Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate'
#   Abs: absorption coeff in cm-1
#   Diam: sample thickness/diameter in mm
    if 'Cylinder' in Geometry:      #Lobanov & Alte da Veiga for 2-theta = 0; beam fully illuminates sample
        MuR = Abs*Diam/20.0
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
    elif 'plate' in Geometry:
        MuR = Abs*Diam/10.
        return math.exp(-MuR)
    elif 'Bragg' in Geometry:
        return 0.0

def Absorb(Geometry,Abs,Diam,Tth,Phi=0,Psi=0):
#Calculate sample absorption
#   Geometry: one of 'Cylinder','Bragg-Brentano','Tilting Flat Plate in transmission','Fixed flat plate'
#   Abs: absorption coeff in cm-1
#   Diam: sample thickness/diameter in mm
#   Tth: 2-theta scattering angle - can be numpy array
#   Phi: flat plate tilt angle - future
#   Psi: flat plate tilt axis - future
    Sth2 = npsind(Tth/2.0)**2
    Cth2 = 1.-Sth2
    if 'Cylinder' in Geometry:      #Lobanov & Alte da Veiga for 2-theta = 0; beam fully illuminates sample
        MuR = Abs*Diam/20.0
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
        MuR = Abs*Diam/10.0
        T1 = np.exp(-MuR)
        T2 = np.exp(-MuR/npcosd(Tth))
        Tb = MuR-MuR/npcosd(Tth)
        return (T2-T1)/Tb
    elif 'Tilting' in Geometry: #assumes symmetric tilt so sample plane is parallel to diffraction vector
        MuR = Abs*Diam/10.0
        cth = npcosd(Tth/2.0)
        return np.exp(-MuR/cth)/cth
        
def Polarization(Pola,Tth,Azm=0.0):
#   Calculate angle dependent x-ray polarization correction (not scaled correctly!)
#   Pola: polarization coefficient e.g 1.0 fully polarized, 0.5 unpolarized
#   Azm: azimuthal angle e.g. 0.0 in plane of polarization
#   Tth: 2-theta scattering angle - can be numpy array
#       which (if either) of these is "right"?
#    return (Pola*npcosd(Azm)**2+(1.-Pola)*npsind(Azm)**2)*npcosd(Tth)**2+ \
#        Pola*npsind(Azm)**2+(1.-Pola)*npcosd(Azm)**2
    return ((1.0-Pola)*npcosd(Azm)**2+Pola*npsind(Azm)**2)*npcosd(Tth)**2+   \
        (1.0-Pola)*npsind(Azm)**2+Pola*npcosd(Azm)**2
    
def Oblique(ObCoeff,Tth):
    if ObCoeff:
        return (1.-ObCoeff)/(1.0-np.exp(np.log(ObCoeff)/npcosd(Tth)))
    else:
        return 1.0
                
def Ruland(RulCoff,wave,Q,Compton):
    C = 2.9978e8
    D = 1.5e-3
    hmc = 0.024262734687
    sinth2 = (Q*wave/(4.0*np.pi))**2
    dlam = (wave**2)*Compton*Q/C
    dlam_c = 2.0*hmc*sinth2-D*wave**2
    return 1.0/((1.0+dlam/RulCoff)*(1.0+(np.pi*dlam_c/(dlam+RulCoff))**2))
    
def LorchWeight(Q):
    return np.sin(np.pi*(Q[-1]-Q)/(2.0*Q[-1]))
            
def GetAsfMean(ElList,Sthl2):
#   Calculate various scattering factor terms for PDF calcs
#   ElList: element dictionary contains scattering factor coefficients, etc.
#   Sthl2: numpy array of sin theta/lambda squared values
#   returns: mean(f^2), mean(f)^2, mean(compton)
    sumNoAtoms = 0.0
    FF = np.zeros_like(Sthl2)
    FF2 = np.zeros_like(Sthl2)
    CF = np.zeros_like(Sthl2)
    for El in ElList:
        sumNoAtoms += ElList[El]['FormulaNo']
    for El in ElList:
        el = ElList[El]
        ff2 = (G2elem.ScatFac(el,Sthl2)+el['fp'])**2+el['fpp']**2
        cf = G2elem.ComptonFac(el,Sthl2)
        FF += np.sqrt(ff2)*el['FormulaNo']/sumNoAtoms
        FF2 += ff2*el['FormulaNo']/sumNoAtoms
        CF += cf*el['FormulaNo']/sumNoAtoms
    return FF2,FF**2,CF
    
def GetNumDensity(ElList,Vol):
    sumNoAtoms = 0.0
    for El in ElList:
        sumNoAtoms += ElList[El]['FormulaNo']
    return sumNoAtoms/Vol
           
def MultGetQ(Tth,MuT,Geometry,b=88.0,a=0.01):
    NS = 500
    Gama = np.linspace(0.,np.pi/2.,NS,False)[1:]
    Cgama = np.cos(Gama)[:,np.newaxis]
    Sgama = np.sin(Gama)[:,np.newaxis]
    CSgama = 1.0/Sgama
    Delt = Gama[1]-Gama[0]
    emc = 7.94e-26
    Navo = 6.023e23
    Cth = npcosd(Tth/2.0)
    CTth = npcosd(Tth)
    Sth = npcosd(Tth/2.0)
    STth = npsind(Tth)
    CSth = 1.0/Sth
    CSTth = 1.0/STth
    SCth = 1.0/Cth
    SCTth = 1.0/CTth
    if 'Bragg' in Geometry:
        G = 8.0*Delt*Navo*emc*Sth/((1.0-CTth**2)*(1.0-np.exp(-2.0*MuT*CSth)))
        Y1 = np.pi
        Y2 = np.pi/2.0
        Y3 = 3.*np.pi/8. #3pi/4?
        W = 1.0/(Sth+np.fabs(Sgama))
        W += np.exp(-MuT*CSth)*(2.0*np.fabs(Sgama)*np.exp(-MuT*np.fabs(CSgama))-
            (Sth+np.fabs(Sgama))*np.exp(-MuT*CSth))/(Sth**2-Sgama**2)
        Fac0 = Sth**2*Sgama**2
        X = Fac0+(Fac0+CTth)**2/2
        Y = Cgama**2*Cth**2*(1.0-Fac0-CTth)
        Z = Cgama**4*Cth**4/2.0
        E = 2.0*(1.0-a)/(b*Cgama/Cth)
        F1 = (2.0+b*(1.0+Sth*Sgama))/(b*Cth*Cgama) #trouble if < 1
        F2 = (2.0+b*(1.0-Sth*Sgama))/(b*Cth*Cgama)
        T1 = np.pi/np.sqrt(F1**2-1.0)
        T2 = np.pi/np.sqrt(F2**2-1.0)
        Y4 = T1+T2
        Y5 = F1**2*T1+F2**2*T2-np.pi*(F1+F2)
        Y6 = F1**4*T1+F2**4*T2-np.pi*(F1+F2)/2.0-np.pi*(F1**3+F2**3)
        Y7 = (T2-T1)/(F1-F2)
        YT = F2**2*T2-F1**2*T1
        Y8 = Y1+YT/(F1-F2)
        Y9 = Y2+(F2**4*T2-F1**4*T1)/(F1-F2)+Y1*((F1+F2)**2-F1*F2)
        M = (a**2*(X*Y1+Y*Y2+Z*Y3)+a*E*(X*Y4+Y*Y5+Z*Y6)+E**2*(X*Y7+Y*Y8+Z*Y9))*Cgama
        
        Q = np.sum(W*M,axis=0)
        return Q*G*NS/(NS-1.)
#
#      cos2th=2.0d*costh^2 - 1.0d
#      G= delta * 8.0d * Na * emc * sinth/(1.0d + cos2th^2)/(1.0d - exp(-2.0d*mut*cscth))
#      Y1=3.1415926d
#      Y2=Y1*0.5d
#      Y3=Y2*0.75d
#      for i=1,num_steps-1 do begin
#         cosgama=double(cos(gama[i]))
#         singama=double(sin(gama[i]))
#         cscgama=1.0d / singama
#
#         W=1.0d/(sinth+abs(singama))
#         W=W+exp(-1.0*mut*cscth)*(2.0d*abs(singama)*exp(-1.0d*mut*abs(cscgama))- $
#                 (sinth+abs(singama))*exp(-1.0d*mut*cscth))/(sinth^2-singama^2)
#
#         factor0=sinth^2*singama^2
#         X=factor0+(factor0+cos2th)^2/2.0d
#         Y=cosgama^2*(1.0d - factor0-cos2th)*costh^2
#         Z=cosgama^4/2.0d*costh^4
#         E=2.0d*(1.0-a)/b/cosgama/costh
#
#         F1=1.0d/b/cosgama*(2.0d + b*(1.0+sinth*singama))/costh
#         F2=1.0d/b/cosgama*(2.0d + b*(1.0-sinth*singama))/costh
#
#         T1=3.14159/sqrt(F1^2-1.0d)
#         T2=3.14159/sqrt(F2^2-1.0d)
#         Y4=T1+T2
#         Y5=F1^2*T1+F2^2*T2-3.14159*(F1+F2)
#         Y6=F1^4*T1+F2^4*T2-3.14159*(F1+F2)/2.0-3.14159*(F1^3+F2^3)
#         Y7=(T2-T1)/(F1-F2)
#         Y8=Y1+(F2^2*T2-F1^2*T1)/(F1-F2)
#         Y9=Y2+(F2^4*T2-F1^4*T1)/(F1-F2)+Y1*((F1+F2)^2-F1*F2)
#         M=(a^2*(X*Y1+Y*Y2+Z*Y3)+a*E*(X*Y4+Y*Y5+Z*Y6)+E^2* $
#                      (X*Y7+Y*Y8+Z*Y9))*cosgama
#
#         Q=Q+W*M
#
#      endfor
#      Q=double(num_steps)/Double(num_steps-1)*Q*G
#      end
    elif 'Cylinder' in Geometry:
        Q = 0.
    elif 'Fixed' in Geometry:   #Dwiggens & Park, Acta Cryst. A27, 264 (1971) with corrections
        EMA = np.exp(-MuT*(1.0-SCTth))
        Fac1 = (1.-EMA)/(1.0-SCTth)
        G = 2.0*Delt*Navo*emc/((1.0+CTth**2)*Fac1)
        Fac0 = Cgama/(1-Sgama*SCTth)
        Wp = Fac0*(Fac1-(EMA-np.exp(-MuT*(CSgama-SCTth)))/(CSgama-1.0))
        Fac0 = Cgama/(1.0+Sgama*SCTth)
        Wm = Fac0*(Fac1+(np.exp(-MuT*(1.0+CSgama))-1.0)/(CSgama+1.0))
        X = (Sgama**2+CTth**2*(1.0-Sgama**2+Sgama**4))/2.0
        Y = Sgama**3*Cgama*STth*CTth
        Z = Cgama**2*(1.0+Sgama**2)*STth**2/2.0
        Fac2 = 1.0/(b*Cgama*STth)
        U = 2.0*(1.0-a)*Fac2
        V = (2.0+b*(1.0-CTth*Sgama))*Fac2
        Mp = 2.0*np.pi*(a+2.0*(1.0-a)/(2.0+b*(1.0-Sgama)))*(a*X+a*Z/2.0-U*Y+U*(X+Y*V+Z*V**2)/np.sqrt(V**2-1.0)-U*Z*V)
        V = (2.0+b*(1.0+CTth*Sgama))*Fac2
        Y = -Y
        Mm = 2.0*np.pi*(a+2.0*(1.0-a)/(2.0+b*(1.0+Sgama)))*(a*X+a*Z/2.0-U*Y+U*(X+Y*V+Z*V**2)/np.sqrt(V**2-1.0)-U*Z*V)
        Q = np.sum(Wp*Mp+Wm*Mm,axis=0)
        return Q*G*NS/(NS-1.)
    elif 'Tilting' in Geometry:
        EMA = np.exp(-MuT*SCth)
        G = 2.0*Delt*Navo*emc/((1.0+CTth**2)*EMA)
#          Wplus = expmutsecth/(1.0d - singama*secth) + singama/mut/(1.0 -singama*secth)/(1.0-singama*secth)* $
#                                                       (Exp(-1.0*mut*cscgama) - expmutsecth)
#          Wminus = expmutsecth/(1.0d + singama*secth) - singama/mut/(1.0 +singama*secth)/(1.0+singama*secth)* $
#                                                        expmutsecth*(1.0d - expmutsecth*Exp(-1.0*mut*cscgama))
        Wp = EMA/(1.0-Sgama*SCth)+Sgama/MuT/(1.0-Sgama*SCth)/(1.0-Sgama*SCth)*(np.exp(-MuT*CSgama)-EMA)
#        Wp = EMA/(1.0-Sgama*SCth)+Sgama/MuT/(1.0-Sgama*SCth)**2*(np.exp(-MuT*CSgama)-EMA)
        Wm = EMA/(1.0+Sgama*SCth)-Sgama/MuT/(1.0+Sgama*SCth)/(1.0+Sgama*SCth)*EMA*(1.0-EMA*np.exp(-MuT*CSgama))
#        Wm = EMA/(1.0+Sgama*SCth)-Sgama/MuT/(1.0+Sgama*SCth)**2*EMA*(1.0-EMA*np.exp(-MuT*CSgama))
        X = 0.5*(Cth**2*(Cth**2*Sgama**4-4.0*Sth**2*Cgama**2)+1.0)
        Y = Cgama**2*(1.0+Cgama**2)*Cth**2*Sth**2
        Z = 0.5*Cgama**4*Sth**4
#          X = 0.5*(costh*costh*(costh*costh*singama*singama*singama*singama - $
#                           4.0*sinth*sinth*cosgama*cosgama) +1.0d)
#
#          Y = cosgama*cosgama*(1.0 + cosgama*cosgama)*costh*costh*sinth*sinth
#
#          Z= 0.5 *cosgama*cosgama*cosgama*cosgama* (sinth^4)
#
        AlP = 2.0+b*(1.0-Cth*Sgama)
        AlM = 2.0+b*(1.0+Cth*Sgama)
#          alphaplus = 2.0 + b*(1.0 - costh*singama)
#          alphaminus = 2.0 + b*(1.0 + costh*singama)
        BeP = np.sqrt(np.fabs(AlP**2-(b*Cgama*Sth)**2))
        BeM = np.sqrt(np.fabs(AlM**2-(b*Cgama*Sth)**2))
#          betaplus = Sqrt(Abs(alphaplus*alphaplus - b*b*cosgama*cosgama*sinth*sinth))
#          betaminus = Sqrt(Abs(alphaminus*alphaminus - b*b*cosgama*cosgama*sinth*sinth))
        Mp = Cgama*(np.pi*a**2*(2.0*X+Y+0.75*Z)+(2.0*np.pi*(1.0-a))*(1.0-a+a*AlP)* \
            (4.0*X/AlP/BeP+(4.0*(1.0+Cgama**2)/b/b*Cth**2)*(AlP/BeP-1.0)+
            2.0/b**4*AlP/BeP*AlP**2-2.0/b**4*AlP**2-Cgama**2/b/b*Sth*2))
#          Mplus = cosgama*(!DPI * a * a * (2.0*x + y + 0.75*z) + $
#                   (2.0*!DPI*(1.0 - a)) *(1.0 - a + a*alphaplus)* $
#                   (4.0*x/alphaplus/betaplus + (4.0*(1.0+cosgama*cosgama)/b/b*costh*costh)*(alphaplus/betaplus -1.0) + $
#                   2.0/(b^4)*alphaplus/betaplus*alphaplus*alphaplus - 2.0/(b^4)*alphaplus*alphaplus - $
#                   cosgama*cosgama/b/b*sinth*sinth))
        Mm =Cgama*(np.pi*a**2*(2.0*X+Y+0.75*Z)+(2.0*np.pi*(1.0-a))*(1.0-a+a*AlM)* \
            (4.0*X/AlM/BeM+(4.0*(1.0+Cgama**2)/b/b*Cth**2)*(AlM/BeM-1.0)+
            2.0/b**4*AlM/BeM*AlM**2-2.0/b**4*AlM**2-Cgama**2/b/b*Sth*2))
#          Mminus = cosgama*(!DPI * a * a * (2.0*x + y + 0.75*z) + $
#                   (2.0*!DPI*(1.0 - a)) *(1.0 - a + a*alphaminus)* $
#                   (4.0*x/alphaminus/betaminus + (4.0*(1.0+cosgama*cosgama)/b/b*costh*costh)*(alphaminus/betaminus -1.0) + $
#                   2.0/(b^4)*alphaminus/betaminus*alphaminus*alphaminus - 2.0/(b^4)*alphaminus*alphaminus - $
#                   cosgama*cosgama/b/b*sinth*sinth))
        Q = np.sum(Wp*Mp+Wm*Mm,axis=0)
        return Q*G*NS/(NS-1.)
#       expmutsecth = Exp(-1.0*mut*secth)
#       G= delta * 2.0 * Na * emc /(1.0+costth^2)/expmutsecth
#       for i=1, num_steps-1 do begin
#          cosgama=double(cos(gama[i]))
#          singama=double(sin(gama[i]))
#          cscgama=1.0d/singama
#
#
#     ; print, "W", min(wplus), max(wplus), min(wminus), max(wminus)
#
#
#
#
#    ;               print, a,b
#  ; print, "M", min(mplus), max(mplus), min(mminus), max(mminus)
#          Q=Q+ Wplus*Mplus + Wminus*Mminus
#      endfor
#      Q=double(num_steps)/double(num_steps-1)*Q*G
#   ;   print, min(q), max(q)
#      end

def MultiScattering(Geometry,ElList,Tth):
    BN = BD = 0.0
    Amu = 0.0
    for El in ElList:
        el = ElList[El]
        BN += el['Z']*el['FormulaNo']
        BD += el['FormulaNo']
        Amu += el['FormulaNo']*el['mu']
        

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

def CalcPDF(data,inst,xydata):
    auxPlot = []
    import copy
    import scipy.fftpack as ft
    #subtract backgrounds - if any
    xydata['IofQ'] = copy.deepcopy(xydata['Sample'])
    if data['Sample Bkg.']['Name']:
        xydata['IofQ'][1][1] += (xydata['Sample Bkg.'][1][1]+
            data['Sample Bkg.']['Add'])*data['Sample Bkg.']['Mult']
    if data['Container']['Name']:
        xycontainer = (xydata['Container'][1][1]+data['Container']['Add'])*data['Container']['Mult']
        if data['Container Bkg.']['Name']:
            xycontainer += (xydata['Container Bkg.'][1][1]+
                data['Container Bkg.']['Add'])*data['Container Bkg.']['Mult']
        xydata['IofQ'][1][1] += xycontainer
    #get element data & absorption coeff.
    ElList = data['ElList']
    Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
    #Apply angle dependent corrections
    Tth = xydata['Sample'][1][0]
    dt = (Tth[1]-Tth[0])
    xydata['IofQ'][1][1] /= Absorb(data['Geometry'],Abs,data['Diam'],Tth)
    xydata['IofQ'][1][1] /= Polarization(inst['Polariz.'],Tth,Azm=inst['Azimuth'])
    if data['DetType'] == 'Image plate':
        xydata['IofQ'][1][1] *= Oblique(data['ObliqCoeff'],Tth)
    XY = xydata['IofQ'][1]    
    #convert to Q
    hc = 12.397639
    if 'Lam' in inst:
        wave = inst['Lam']
    else:
        wave = inst['Lam1']
    keV = hc/wave
    minQ = npT2q(Tth[0],wave)
    maxQ = npT2q(Tth[-1],wave)    
    Qpoints = np.linspace(0.,maxQ,len(XY[0]),endpoint=True)
    dq = Qpoints[1]-Qpoints[0]
    XY[0] = npT2q(XY[0],wave)    
#    Qdata = np.nan_to_num(si.griddata(XY[0],XY[1],Qpoints,method='linear')) #only OK for scipy 0.9!
    T = si.interp1d(XY[0],XY[1],bounds_error=False,fill_value=0.0)      #OK for scipy 0.8
    Qdata = T(Qpoints)
    
    qLimits = data['QScaleLim']
    minQ = np.searchsorted(Qpoints,qLimits[0])
    maxQ = np.searchsorted(Qpoints,qLimits[1])
    newdata = []
    xydata['IofQ'][1][0] = Qpoints
    xydata['IofQ'][1][1] = Qdata
    for item in xydata['IofQ'][1]:
        newdata.append(item[:maxQ])
    xydata['IofQ'][1] = newdata
    

    xydata['SofQ'] = copy.deepcopy(xydata['IofQ'])
    FFSq,SqFF,CF = GetAsfMean(ElList,(xydata['SofQ'][1][0]/(4.0*np.pi))**2)  #these are <f^2>,<f>^2,Cf
    Q = xydata['SofQ'][1][0]
    ruland = Ruland(data['Ruland'],wave,Q,CF)
#    auxPlot.append([Q,ruland,'Ruland'])      
    CF *= ruland
#    auxPlot.append([Q,CF,'CF-Corr'])
    scale = np.sum((FFSq+CF)[minQ:maxQ])/np.sum(xydata['SofQ'][1][1][minQ:maxQ])
    xydata['SofQ'][1][1] *= scale
    xydata['SofQ'][1][1] -= CF
    xydata['SofQ'][1][1] = xydata['SofQ'][1][1]/SqFF
    scale = len(xydata['SofQ'][1][1][minQ:maxQ])/np.sum(xydata['SofQ'][1][1][minQ:maxQ])
    xydata['SofQ'][1][1] *= scale
    
    xydata['FofQ'] = copy.deepcopy(xydata['SofQ'])
    xydata['FofQ'][1][1] = xydata['FofQ'][1][0]*(xydata['SofQ'][1][1]-1.0)
    if data['Lorch']:
        xydata['FofQ'][1][1] *= LorchWeight(Q)
    
    xydata['GofR'] = copy.deepcopy(xydata['FofQ'])
    nR = len(xydata['GofR'][1][1])
    xydata['GofR'][1][1] = -dq*np.imag(ft.fft(xydata['FofQ'][1][1],4*nR)[:nR])
    xydata['GofR'][1][0] = 0.5*np.pi*np.linspace(0,nR,nR)/qLimits[1]
    
        
    return auxPlot
        
