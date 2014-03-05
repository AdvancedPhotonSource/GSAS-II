#/usr/bin/env python
# -*- coding: utf-8 -*-
'''
*GSASII small angle calculation module*
==================================

'''
########### SVN repository information ###################
# $Date: 2014-01-09 11:09:53 -0600 (Thu, 09 Jan 2014) $
# $Author: vondreele $
# $Revision: 1186 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIsasd.py $
# $Id: GSASIIsasd.py 1186 2014-01-09 17:09:53Z vondreele $
########### SVN repository information ###################
import sys
import math
import time

import numpy as np
import scipy as sp
import numpy.linalg as nl
from numpy.fft import ifft, fft, fftshift
import scipy.special as scsp
import scipy.interpolate as si
import scipy.stats as st
import scipy.optimize as so

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1186 $")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIgrid as G2gd
import GSASIIIO as G2IO
import GSASIImath as G2mth
import pypowder as pyd

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
    
###############################################################################
#### Particle form factors & volumes as class definitions
###############################################################################

class SASDParticles(object):
    def __init__(self,name=None,parNames=None):
        self.name = name
        self.ParmNames = parmNames
        
    def FormFactor():
        FF = 0.
        return FF
        
    def Volume():
        Vol = 0.
        return Vol
        
class Sphere(SASDParticles):
    def __init__(self,name=None,parmNames=None):
        self.Name = name
        if self.Name == None:
            self.Name = 'Sphere'
        self.ParmNames = parmNames
        if self.ParmNames == None:
            self.ParmNames = ['Radius',]
        
    def FormFactor(self,Q,Parms):
        ''' Compute hard sphere form factor - can use numpy arrays
        param float:Q Q value array (usually in A-1)
        param float:R sphere radius (Usually in A - must match Q-1 units)
        returns float: form factors as array as needed
        '''
        QR = Q*Parms[0]
        return (3./(QR**3))*(np.sin(QR)-(QR*np.cos(QR)))
        
    def Volume(self,Parms):
        ''' Compute volume of sphere
        - numpy array friendly
        param float:R sphere radius
        returns float: volume
        '''
        return (4./3.)*np.pi*Parms[0]**3
        
class Spheroid(SASDParticles):
    def __init__(self,name=None,parmNames=None):
        self.Name = name
        if self.Name == None:
            self.Name = 'Spheroid'
        self.ParmNames = parmNames
        if self.ParmNames == None:
            self.ParmNames = ['Radius','Aspect ratio']
    
    def FormFactor(self,Q,Parms):
        ''' Compute form factor of cylindrically symmetric ellipsoid (spheroid) 
        - can use numpy arrays for R & AR; will return corresponding numpy array
        param float:Q Q value array (usually in A-1)
        param float R: radius along 2 axes of spheroid
        param float AR: aspect ratio so 3rd axis = R*AR
        returns float: form factors as array as needed
        '''
        R,AR = Parms
        NP = 50
        if 0.99 < AR < 1.01:
            return SphereFF(Q,R,0)
        else:
            cth = np.linspace(0,1.,NP)
            Rct = R*np.sqrt(1.+(AR**2-1.)*cth**2)
            return np.sqrt(np.sum(SphereFF(Q[:,np.newaxis],Rct,0)**2,axis=1)/NP)
        
    def Volume(self,Parms):
        ''' Compute volume of cylindrically symmetric ellipsoid (spheroid) 
        - numpy array friendly
        param float R: radius along 2 axes of spheroid
        param float AR: aspect ratio so radius of 3rd axis = R*AR
        returns float: volume
        '''
        R,AR = Parms
        return AR*(4./3.)*np.pi*R**3
        
###############################################################################
#### Particle form factors
###############################################################################

def SphereFF(Q,R,dummy=0):
    ''' Compute hard sphere form factor - can use numpy arrays
    param float:Q Q value array (usually in A-1)
    param float:R sphere radius (Usually in A - must match Q-1 units)
    returns float: form factors as array as needed
    '''
    QR = Q*R
    return (3./(QR**3))*(np.sin(QR)-(QR*np.cos(QR)))
    
def SpheroidFF(Q,R,AR):
    ''' Compute form factor of cylindrically symmetric ellipsoid (spheroid) 
    - can use numpy arrays for R & AR; will return corresponding numpy array
    param float:Q Q value array (usually in A-1)
    param float R: radius along 2 axes of spheroid
    param float AR: aspect ratio so 3rd axis = R*AR
    returns float: form factors as array as needed
    '''
    NP = 50
    if 0.99 < AR < 1.01:
        return SphereFF(Q,R,0)
    else:
        cth = np.linspace(0,1.,NP)
        Rct = R*np.sqrt(1.+(AR**2-1.)*cth**2)
        return np.sqrt(np.sum(SphereFF(Q[:,np.newaxis],Rct,0)**2,axis=1)/NP)
            
def CylinderFF(Q,R,L):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: R cylinder radius (A)
    param float: L cylinder length (A)
    returns float: form factor
    '''
    NP = 200
    alp = np.linspace(0,np.pi/2.,NP)
    LBessArg = 0.5*Q[:,np.newaxis]*L*np.cos(alp)
    LBess = np.where(LBessArg<1.e-6,1.,np.sin(LBessArg)/LBessArg)
    SBessArg = Q[:,np.newaxis]*R*np.sin(alp)
    SBess = np.where(SBessArg<1.e-6,0.5,scsp.jv(1,SBessArg)/SBessArg)
    return np.sqrt(2.*np.pi*np.sum(np.sin(alp)*(LBess*SBess)**2,axis=1)/NP)
    
def CylinderDFF(Q,L,D):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: L cylinder length (A)
    param float: D cylinder diameter (A)
    returns float: form factor
    '''
    return CylinderFF(Q,D/2.,L)    
    
def CylinderARFF(Q,R,AR): 
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: R cylinder radius (A)
    param float: AR cylinder aspect ratio = L/D = L/2R
    returns float: form factor
    '''
    return CylinderFF(Q,R,2.*R*AR)    
    
def UniSphereFF(Q,R,dummy=0):
    Rg = np.sqrt(3./5.)*R
    B = 1.62/(Rg**4)    #are we missing *np.pi? 1.62 = 6*(3/5)**2/(4/3) sense?
    QstV = Q/(scsp.erf(Q*Rg/np.sqrt(6)))**3
    return np.sqrt(np.exp((-Q**2*Rg**2)/3.)+(B/QstV**4))
    
def UniRodFF(Q,R,L):
    Rg2 = np.sqrt(R**2/2+L**2/12)
    B2 = np.pi/L
    Rg1 = np.sqrt(3.)*R/2.
    G1 = (2./3.)*R/L
    B1 = 4.*(L+R)/(R**3*L**2)
    QstV = Q/(scsp.erf(Q*Rg2/np.sqrt(6)))**3
    FF = np.exp(-Q**2*Rg2**2/3.)+(B2/QstV)*np.exp(-Rg1**2*Q**2/3.)
    QstV = Q/(scsp.erf(Q*Rg1/np.sqrt(6)))**3
    FF += G1*np.exp(-Q**2*Rg1**2/3.)+(B1/QstV**4)
    return np.sqrt(FF)
    
def UniRodARFF(Q,R,AR):
    return UniRodFF(Q,R,AR*R)
    
def UniDiskFF(Q,R,T):
    Rg2 = np.sqrt(R**2/2.+T**2/12.)
    B2 = 2./R**2
    Rg1 = np.sqrt(3.)*T/2.
    RgC2 = 1.1*Rg1
    G1 = (2./3.)*(T/R)**2
    B1 = 4.*(T+R)/(R**3*T**2)
    QstV = Q/(scsp.erf(Q*Rg2/np.sqrt(6)))**3
    FF = np.exp(-Q**2*Rg2**2/3.)+(B2/QstV**2)*np.exp(-RgC2**2*Q**2/3.)
    QstV = Q/(scsp.erf(Q*Rg1/np.sqrt(6)))**3
    FF += G1*np.exp(-Q**2*Rg1**2/3.)+(B1/QstV**4)
    return np.sqrt(FF)
    
def UniTubeFF(Q,R,L,T):
    Ri = R-T
    DR2 = R**2-Ri**2
    Vt = np.pi*DR2*L
    Rg3 = np.sqrt(DR2/2.+L**2/12.)
    B1 = 4.*np.pi**2*(DR2+L*(R+Ri))/Vt**2
    B2 = np.pi**2*T/Vt
    B3 = np.pi/L
    QstV = Q/(scsp.erf(Q*Rg3/np.sqrt(6)))**3
    FF = np.exp(-Q**2*Rg3**2/3.)+(B3/QstV)*np.exp(-Q**2*R**2/3.)
    QstV = Q/(scsp.erf(Q*R/np.sqrt(6)))**3
    FF += (B2/QstV**2)*np.exp(-Q**2*T**2/3.)
    QstV = Q/(scsp.erf(Q*T/np.sqrt(6)))**3
    FF += B1/QstV**4
    return np.sqrt(FF)
    
    
    
    

###############################################################################
#### Particle volumes
###############################################################################

def SphereVol(R):
    ''' Compute volume of sphere
    - numpy array friendly
    param float:R sphere radius
    returns float: volume
    '''
    return (4./3.)*np.pi*R**3

def SpheroidVol(R,AR):
    ''' Compute volume of cylindrically symmetric ellipsoid (spheroid) 
    - numpy array friendly
    param float R: radius along 2 axes of spheroid
    param float AR: aspect ratio so radius of 3rd axis = R*AR
    returns float: volume
    '''
    return AR*SphereVol(R)
    
def CylinderVol(R,L):
    ''' Compute cylinder volume for radius & length
    - numpy array friendly
    param float: R diameter (A)
    param float: L length (A)
    returns float:volume (A^3)
    '''
    return np.pi*L*R**2
    
def CylinderDVol(L,D):
    ''' Compute cylinder volume for length & diameter
    - numpy array friendly
    param float: L length (A)
    param float: D diameter (A)
    returns float:volume (A^3)
    '''
    return CylinderVol(D/2.,L)
    
def CylinderARVol(R,AR):
    ''' Compute cylinder volume for radius & aspect ratio = L/D
    - numpy array friendly
    param float: R radius (A)
    param float: AR=L/D=L/2R aspect ratio
    returns float:volume
    '''
    return CylinderVol(R,2.*R*AR)
        
    
    
###############################################################################
#### Size distribution
###############################################################################

def SizeDistribution(Profile,ProfDict,Limits,Substances,Sample,data):
    if data['Size']['logBins']:
        Bins = np.logspace(np.log10(data['Size']['MinDiam']),np.log10(data['Size']['MaxDiam']),
            data['Size']['Nbins']+1,True)
    else:
        Bins = np.linspace(data['Size']['MinDiam'],data['Size']['MaxDiam'],
            data['Size']['Nbins']+1,True)
    Dbins = np.diff(Bins)
    BinMag = Dbins*np.ones_like(Dbins)
    print np.sum(BinMag)
        
    print data['Size']
#    print Limits
#    print Substances
#    print Sample
#    print Profile
    
