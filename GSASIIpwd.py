#/usr/bin/env python
# -*- coding: utf-8 -*-
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
from numpy.fft import ifft, fft, fftshift
import scipy.interpolate as si
import scipy.stats as st
import scipy.optimize as so

import GSASIIpath
#import pypowder as pyp
import GSASIIplot as G2plt
import GSASIIlattice as G2lat
import GSASIIElem as G2elem
import GSASIIgrid as G2gd
import GSASIIIO as G2IO

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

def ValuesOut(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def ValuesIn(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
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
        
#GSASII peak fitting routine: Finger, Cox & Jephcoat model        

np.seterr(divide='ignore')

# Voigt function = convolution (norm,cauchy)

class voigt_gen(st.rv_continuous):
    '''
    Voigt function
    Parameters
    -----------------------------------------
    '''
    def _pdf(self,x,s,g):
        z = np.empty([2,len(x)])
        z[0] = st.norm.pdf(x,scale=s)
        z[1] = st.cauchy.pdf(x,scale=g)
        Z = fft(z)
        return fftshift(ifft(Z.prod(axis=0))).real
#    def _cdf(self, x):
#    def _ppf(self, q):
#    def _sf(self, x):
#    def _isf(self, q):
#    def _stats(self):
#    def _entropy(self):

voigt = voigt_gen(name='voigt',shapes='ts,g')
       
# Finger-Cox_Jephcoat D(2phi,2th) function for S/L = H/L

class fcjde_gen(st.rv_continuous):
    """
    Finger-Cox-Jephcoat D(2phi,2th) function for S/L = H/L
    Ref: J. Appl. Cryst. (1994) 27, 892-900.
    Parameters
    -----------------------------------------
    x: array -1 to 1
    t: 2-theta position of peak
    s: sum(S/L,H/L); S: sample height, H: detector opening, 
        L: sample to detector opening distance
    dx: 2-theta step size in deg
    Result for fcj.pdf
    -----------------------------------------
    T = x*dx+t
    s = S/L+H/L
    if x < 0: 
        fcj.pdf = [1/sqrt({cos(T)**2/cos(t)**2}-1) - 1/s]/|cos(T)|    
    if x >= 0:
        fcj.pdf = 0    
    """
    def _pdf(self,x,t,s,dx):
        T = dx*x+t
        ax = npcosd(T)**2
        bx = npcosd(t)**2
        bx = np.where(ax>bx,bx,ax)
        fx = np.where(ax>bx,(np.sqrt(bx/(ax-bx))-1./s)/np.sqrt(ax),0.0)
        fx = np.where(fx > 0.,fx,0.0)
        return fx
#    def _cdf(self, x):
#    def _ppf(self, q):
#    def _sf(self, x):
#    def _isf(self, q):
#    def _stats(self):
#    def _entropy(self):
        
fcjde = fcjde_gen(name='fcjde',shapes='t,s,dx')
                
def getBackground(parmDict,bakType,xdata):
    if bakType == 'chebyschev':
        yb = np.zeros_like(xdata)
        iBak = 0
        while True:
            key = 'Back'+str(iBak)
            try:
                yb += parmDict[key]*(xdata-xdata[0])**iBak
                iBak += 1
            except KeyError:
                return yb
                
def getPeakProfile(parmDict,xdata,varyList,bakType):
    
    def calcPeakFFT(x,fxns,widths,pos,args):
        dx = x[1]-x[0]
        z = np.empty([len(fxns),len(x)])
        for i,(fxn,width,arg) in enumerate(zip(fxns,widths,args)):
            z[i] = fxn.pdf(x,loc=pos,scale=width,*arg)*dx
        Z = fft(z)
        if len(fxns)%2:
            return ifft(Z.prod(axis=0)).real
        else:
            return fftshift(ifft(Z.prod(axis=0))).real
                        
    def getFCJVoigt(pos,intens,sig,gam,shl,xdata):
        
        DX = xdata[1]-xdata[0]
        widths,fmin,fmax = getWidths(pos,sig,gam,shl)
        x = np.linspace(pos-fmin,pos+fmin,1024)
        if pos > 90:
            fmin,fmax = [fmax,fmin]          
        dx = x[1]-x[0]
        arg = [pos,shl/10.,dx,]
        Df = fcjde.pdf(x,*arg,loc=pos,scale=1.0)
        if len(np.nonzero(Df)[0])>5:
            fxns = [st.norm,st.cauchy,fcjde]
            args = [[],[],arg]
        else:
            fxns = [st.norm,st.cauchy]
            args = [[],[]]
        Df = calcPeakFFT(x,fxns,widths,pos,args)
        Df /= np.sum(Df)
        Df = si.interp1d(x,Df,bounds_error=False,fill_value=0.0)
        return intens*Df(xdata)*DX/dx       #*10 to get close to old fxn???
                    
    yb = getBackground(parmDict,bakType,xdata)
    yc = np.zeros_like(yb)
    U = parmDict['U']
    V = parmDict['V']
    W = parmDict['W']
    X = parmDict['X']
    Y = parmDict['Y']
    shl = max(parmDict['SH/L'],0.001)
    Ka2 = False
    if 'Lam1' in parmDict.keys():
        Ka2 = True
        lamRatio = parmDict['Lam2']/parmDict['Lam1']
        kRatio = parmDict['I(L2)/I(L1)']
    iPeak = 0
    while True:
        try:
            pos = parmDict['pos'+str(iPeak)]
            intens = parmDict['int'+str(iPeak)]
            sigName = 'sig'+str(iPeak)
            if sigName in varyList:
                sig = parmDict[sigName]
            else:
                sig = U*tand(pos/2.0)**2+V*tand(pos/2.0)+W
            sig = max(sig,0.001)          #avoid neg sigma
            gamName = 'gam'+str(iPeak)
            if gamName in varyList:
                gam = parmDict[gamName]
            else:
                gam = X/cosd(pos/2.0)+Y*tand(pos/2.0)
            gam = max(gam,0.01)             #avoid neg gamma
            Wd,fmin,fmax = getWidths(pos,sig,gam,shl)
            if pos > 90:
                fmin,fmax = [fmax,fmin]          
            iBeg = np.searchsorted(xdata,pos-fmin)
            iFin = np.searchsorted(xdata,pos+fmax)
            if not iBeg+iFin:       #peak below low limit
                iPeak += 1
                continue
            elif not iBeg-iFin:     #peak above high limit
                return yb+yc
            yc[iBeg:iFin] += getFCJVoigt(pos,intens,sig,gam,shl,xdata[iBeg:iFin])
            if Ka2:
                pos2 = 2.0*asind(lamRatio*sind(pos/2.0))
                Wd,fmin,fmax = getWidths(pos2,sig,gam,shl)
                if pos > 90:
                    fmin,fmax = [fmax,fmin]          
                iBeg = np.searchsorted(xdata,pos2-fmin)
                iFin = np.searchsorted(xdata,pos2+fmax)
                yc[iBeg:iFin] += getFCJVoigt(pos2,intens*kRatio,sig,gam,shl,xdata[iBeg:iFin])
            iPeak += 1
        except KeyError:        #no more peaks to process
            return yb+yc
        
def getWidths(pos,sig,gam,shl):
    if shl:
        widths = [np.sqrt(sig)/100.,gam/200.,.1]
    else:
        widths = [np.sqrt(sig)/100.,gam/200.]
    fwhm = 2.355*widths[0]+2.*widths[1]
    fmin = 10.*(fwhm+shl*abs(npcosd(pos)))
    fmax = 15.0*fwhm
    return widths,fmin,fmax
                
def DoPeakFit(FitPgm,Peaks,Background,Limits,Inst,data):
    
    def SetBackgroundParms(Background):
        bakType,bakFlag = Background[:2]
        backVals = Background[3:]
        backNames = ['Back'+str(i) for i in range(len(backVals))]
        if bakFlag: #returns backNames as varyList = backNames
            return bakType,dict(zip(backNames,backVals)),backNames
        else:       #no background varied; varyList = []
            return bakType,dict(zip(backNames,backVals)),[]
        
    def GetBackgroundParms(parmList,Background):
        iBak = 0
        while True:
            try:
                bakName = 'Back'+str(iBak)
                Background[iBak+3] = parmList[bakName]
                iBak += 1
            except KeyError:
                break
                
    def BackgroundPrint(Background,sigDict):
        if Background[1]:
            print 'Background coefficients for',Background[0],'function'
            ptfmt = "%12.5f"
            ptstr =  'values:'
            sigstr = 'esds  :'
            for i,back in enumerate(Background[3:]):
                ptstr += ptfmt % (back)
                sigstr += ptfmt % (sigDict['Back'+str(i)])
            print ptstr
            print sigstr
        else:
            print 'Background not refined'
            
    def SetInstParms(Inst):
        insVals,insFlags,insNames = Inst[1:4]
        dataType = insVals[0]
        insVary = []
        for i,flag in enumerate(insFlags):
            if flag:
                insVary.append(insNames[i])
        instDict = dict(zip(insNames,insVals))
        instDict['X'] = max(instDict['X'],0.1)
        instDict['Y'] = max(instDict['Y'],0.1)
        instDict['SH/L'] = max(instDict['SH/L'],0.001)
        return dataType,instDict,insVary
        
    def GetInstParms(parmDict,Inst,varyList,Peaks):
        instNames = Inst[3]
        for i,name in enumerate(instNames):
            Inst[1][i] = parmDict[name]
        iPeak = 0
        while True:
            try:
                sigName = 'sig'+str(iPeak)
                pos = parmDict['pos'+str(iPeak)]
                if sigName not in varyList:
                    parmDict[sigName] = parmDict['U']*tand(pos/2.0)**2+parmDict['V']*tand(pos/2.0)+parmDict['W']
                gamName = 'gam'+str(iPeak)
                if gamName not in varyList:
                    parmDict[gamName] = parmDict['X']/cosd(pos/2.0)+parmDict['Y']*tand(pos/2.0)
                iPeak += 1
            except KeyError:
                break
        
    def InstPrint(Inst,sigDict):
        print 'Instrument Parameters:'
        ptfmt = "%12.6f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        instNames = Inst[3][1:]
        for i,name in enumerate(instNames):
            ptlbls += "%s" % (name.center(12))
            ptstr += ptfmt % (Inst[1][i+1])
            if name in sigDict:
                sigstr += ptfmt % (sigDict[name])
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr

    def setPeaksParms(Peaks):
        peakNames = []
        peakVary = []
        peakVals = []
        names = ['pos','int','sig','gam']
        for i,peak in enumerate(Peaks):
            for j in range(4):
                peakVals.append(peak[2*j])
                parName = names[j]+str(i)
                peakNames.append(parName)
                if peak[2*j+1]:
                    peakVary.append(parName)
        return dict(zip(peakNames,peakVals)),peakVary
                
    def GetPeaksParms(parmDict,Peaks,varyList):
        names = ['pos','int','sig','gam']
        for i,peak in enumerate(Peaks):
            for j in range(4):
                pos = parmDict['pos'+str(i)]
                parName = names[j]+str(i)
                if parName in varyList:
                    peak[2*j] = parmDict[parName]
                elif 'sig' in parName:
                    peak[2*j] = parmDict['U']*tand(pos/2.0)**2+parmDict['V']*tand(pos/2.0)+parmDict['W']
                elif 'gam' in parName:
                    peak[2*j] = parmDict['X']/cosd(pos/2.0)+parmDict['Y']*tand(pos/2.0)
                        
    def PeaksPrint(parmDict,sigDict,varyList):
        print 'Peak coefficients:'
        names = ['pos','int','sig','gam']
        head = 15*' '
        for name in names:
            head += name.center(12)+'esd'.center(12)
        print head
        ptfmt = {'pos':"%12.5f",'int':"%12.1f",'sig':"%12.3f",'gam':"%12.3f"}
        for i,peak in enumerate(Peaks):
            ptstr =  ':'
            for j in range(4):
                name = names[j]
                parName = name+str(i)
                ptstr += ptfmt[name] % (parmDict[parName])
                if parName in varyList:
#                    ptstr += G2IO.ValEsd(parmDict[parName],sigDict[parName])
                    ptstr += ptfmt[name] % (sigDict[parName])
                else:
#                    ptstr += G2IO.ValEsd(parmDict[parName],0.0)
                    ptstr += 12*' '
            print '%s'%(('Peak'+str(i+1)).center(8)),ptstr
            
    def errPeakProfile(values, xdata, ydata, weights, parmdict, varylist,bakType,dlg):        
        parmdict.update(zip(varylist,values))
        M = np.sqrt(weights)*(ydata-getPeakProfile(parmdict,xdata,varylist,bakType))
        Rwp = min(100.,np.sqrt(np.sum(M**2)/np.sum(weights*ydata**2))*100.)
        if dlg:
            GoOn = dlg.Update(Rwp,newmsg='%s%8.3f%s'%('Peak fit Rwp =',Rwp,'%'))[0]
            if not GoOn:
                return -M           #abort!!
        return M
        
    x,y,w,yc,yb,yd = data               #these are numpy arrays!
    xBeg = np.searchsorted(x,Limits[0])
    xFin = np.searchsorted(x,Limits[1])
    bakType,bakDict,bakVary = SetBackgroundParms(Background)
    dataType,insDict,insVary = SetInstParms(Inst)
    peakDict,peakVary = setPeaksParms(Peaks)
    parmDict = {}
    parmDict.update(bakDict)
    parmDict.update(insDict)
    parmDict.update(peakDict)
    varyList = bakVary+insVary+peakVary
    while True:
        begin = time.time()
        values =  np.array(ValuesOut(parmDict, varyList))
        if FitPgm == 'LSQ':
            dlg = wx.ProgressDialog('Residual','Peak fit Rwp = ',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
            screenSize = wx.ClientDisplayRect()
            Size = dlg.GetSize()
            dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
            try:
                result = so.leastsq(errPeakProfile,values,
                    args=(x[xBeg:xFin],y[xBeg:xFin],w[xBeg:xFin],parmDict,varyList,bakType,dlg),full_output=True)
            finally:
                dlg.Destroy()
            runtime = time.time()-begin    
            print 'Number of function calls:',result[2]['nfev'],' Number of observations: ',xFin-xBeg,' Number of parameters: ',len(varyList)
            print "%s%8.3f%s " % ('fitpeak time =',runtime,'s')
            ValuesIn(parmDict, varyList, result[0])
            chisq = np.sum(errPeakProfile(result[0],x[xBeg:xFin],y[xBeg:xFin],w[xBeg:xFin],parmDict,varyList,bakType,0)**2)
            Rwp = np.sqrt(chisq/np.sum(w[xBeg:xFin]*y[xBeg:xFin]**2))*100.      #to %
            GOF = chisq/(xFin-xBeg-len(varyList))
            print "%s%7.2f%s%12.6g%s%6.2f" % ('Rwp = ',Rwp,'%, chi**2 = ',chisq,' reduced chi**2 = ',GOF)
            try:
                sig = np.sqrt(np.diag(result[1])*chisq)
                if np.any(np.isnan(sig)):
                    print '*** Least squares aborted - some invalid esds possible ***'
                break                   #refinement succeeded - finish up!
            except ValueError:          #result[1] is None on singular matrix
                print 'Refinement failed - singular matrix'
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        print 'Removing parameter: ',varyList[ipvt-1]
                        del(varyList[ipvt-1])
                        break
        elif FitPgm == 'BFGS':
            print 'Other program here'
            return
        
    sigDict = dict(zip(varyList,sig))
    yb[xBeg:xFin] = getBackground(parmDict,bakType,x[xBeg:xFin])
    yc[xBeg:xFin] = getPeakProfile(parmDict,x[xBeg:xFin],varyList,bakType)
    yd[xBeg:xFin] = y[xBeg:xFin]-yc[xBeg:xFin]
    GetBackgroundParms(parmDict,Background)
    BackgroundPrint(Background,sigDict)
    GetInstParms(parmDict,Inst,varyList,Peaks)
    InstPrint(Inst,sigDict)
    GetPeaksParms(parmDict,Peaks,varyList)    
    PeaksPrint(parmDict,sigDict,varyList)
    
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
        
#testing data
import plot
NeedTestData = True
def TestData():
#    global NeedTestData
    NeedTestData = False
    global bakType
    bakType = 'chebyschev'
    global xdata
    xdata = np.linspace(4.0,40.0,36000)
    global parmDict0
    parmDict0 = {
        'pos0':5.6964,'int0':8835.8,'sig0':1.0,'gam0':1.0,
        'pos1':11.4074,'int1':3922.3,'sig1':1.0,'gam1':1.0,
        'pos2':20.6426,'int2':1573.7,'sig2':1.0,'gam2':1.0,
        'pos3':26.9568,'int3':925.1,'sig3':1.0,'gam3':1.0,
        'U':1.163,'V':-0.605,'W':0.093,'X':0.0,'Y':2.183,'SH/L':0.002,
        'Back0':5.384,'Back1':-0.015,'Back2':.004,
        }
    global parmDict1
    parmDict1 = {
        'pos0':13.4924,'int0':48697.6,'sig0':1.0,'gam0':1.0,
        'pos1':23.4360,'int1':43685.5,'sig1':1.0,'gam1':1.0,
        'pos2':27.1152,'int2':123712.6,'sig2':1.0,'gam2':1.0,
        'pos3':33.7196,'int3':65349.4,'sig3':1.0,'gam3':1.0,
        'pos4':36.1119,'int4':115829.8,'sig4':1.0,'gam4':1.0,
        'pos5':39.0122,'int5':6916.9,'sig5':1.0,'gam5':1.0,
        'U':22.75,'V':-17.596,'W':10.594,'X':1.577,'Y':5.778,'SH/L':0.002,
        'Back0':36.897,'Back1':-0.508,'Back2':.006,
        'Lam1':1.540500,'Lam2':1.544300,'I(L2)/I(L1)':0.5,
        }
    global varyList
    varyList = []

def test0():
    if NeedTestData: TestData()
    msg = 'test '
    gplot = plotter.add('FCJ-Voigt, 11BM').gca()
    gplot.plot(xdata,getBackground(parmDict0,bakType,xdata))   
    gplot.plot(xdata,getPeakProfile(parmDict0,xdata,varyList,bakType))
    fplot = plotter.add('FCJ-Voigt, Ka1+2').gca()
    fplot.plot(xdata,getBackground(parmDict1,bakType,xdata))   
    fplot.plot(xdata,getPeakProfile(parmDict1,xdata,varyList,bakType))
    time0 = time.time()
    for i in range(100):
        y = getPeakProfile(parmDict1,xdata,varyList,bakType)
    print '100+6*Ka1-2 peaks=3600 peaks',time.time()-time0
    

if __name__ == '__main__':
    global plotter
    plotter = plot.PlotNotebook()
    test0()
    print "OK"
    plotter.StartEventLoop()
