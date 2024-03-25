#/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
Classes and routines defined in :mod:`GSASIIsasd` follow. 
'''
from __future__ import division, print_function
import os
import math

import numpy as np
import scipy.special as scsp
import scipy.optimize as so
#import pdb

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIpwd as G2pwd

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
# npT2stl = lambda tth, wave: 2.0*npsind(tth/2.0)/wave
# npT2q = lambda tth,wave: 2.0*np.pi*npT2stl(tth,wave)
    
###############################################################################
#### Particle form factors
###############################################################################

def SphereFF(Q,R,args=()):
    ''' Compute hard sphere form factor - can use numpy arrays
    :param float Q: Q value array (usually in A-1)
    :param float R: sphere radius (Usually in A - must match Q-1 units)
    :param array args: ignored
    :returns: form factors as array as needed (float)
    '''
    QR = Q[:,np.newaxis]*R
    return (3./(QR**3))*(np.sin(QR)-(QR*np.cos(QR)))
    
def SphericalShellFF(Q,R,args=()):
    ''' Compute spherical shell form factor - can use numpy arrays
    :param float Q: Q value array (usually in A-1)
    :param float R: sphere radius (Usually in A - must match Q-1 units)
    :param array args: [float r]: controls the shell thickness: R_inner = min(r*R,R), R_outer = max(r*R,R)
    :returns float: form factors as array as needed
    
	Contributed by: L.A. Avakyan, Southern Federal University, Russia
    '''
    r = args[0]
    if r < 0: # truncate to positive value
        r = 0
    if r < 1:  # r controls inner sphere radius
        return SphereFF(Q,R) - SphereFF(Q,R*r)
    else:      # r controls outer sphere radius
        return SphereFF(Q,R*r) - SphereFF(Q,R)

def SpheroidFF(Q,R,args):
    ''' Compute form factor of cylindrically symmetric ellipsoid (spheroid) 
    - can use numpy arrays for R & AR; will return corresponding numpy array
    param float Q : Q value array (usually in A-1)
    param float R: radius along 2 axes of spheroid
    param array args: [float AR]: aspect ratio so 3rd axis = R*AR
    returns float: form factors as array as needed
    '''
    NP = 50
    AR = args[0]
    if 0.99 < AR < 1.01:
        return SphereFF(Q,R,0)
    else:
        cth = np.linspace(0,1.,NP)
        try:
            Rct = R[:,np.newaxis]*np.sqrt(1.+(AR**2-1.)*cth**2)
        except:
            Rct = R*np.sqrt(1.+(AR**2-1.)*cth**2)
        return np.sqrt(np.sum(SphereFF(Q[:,np.newaxis],Rct,0)**2,axis=2)/NP)
            
def CylinderFF(Q,R,args):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float L]: cylinder length (A)
    returns float: form factor
    '''
    L = args[0]
    NP = 200
    alp = np.linspace(0,np.pi/2.,NP)
    QL = Q[:,np.newaxis]*L
    LBessArg = 0.5*QL[:,:,np.newaxis]**np.cos(alp)
    LBess = np.where(LBessArg<1.e-6,1.,np.sin(LBessArg)/LBessArg)
    QR = Q[:,np.newaxis]*R
    SBessArg = QR[:,:,np.newaxis]*np.sin(alp)
    SBess = np.where(SBessArg<1.e-6,0.5,scsp.jv(1,SBessArg)/SBessArg)
    LSBess = LBess*SBess
    return np.sqrt(2.*np.pi*np.sum(np.sin(alp)*LSBess**2,axis=2)/NP)
    
def CylinderDFF(Q,L,args):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float Q: Q value array (A-1)
    param float L: cylinder half length (A)
    param array args: [float R]: cylinder radius (A)
    returns float: form factor
    '''
    R = args[0]
    return CylinderFF(Q,R,args=[2.*L])    
    
def CylinderARFF(Q,R,args): 
    ''' Compute form factor for cylinders - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float AR]: cylinder aspect ratio = L/D = L/2R
    returns float: form factor
    '''
    AR = args[0]
    return CylinderFF(Q,R,args=[2.*R*AR])    
    
def UniSphereFF(Q,R,args=0):
    ''' Compute form factor for unified sphere - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: ignored
    returns float: form factor
    '''
    Rg = np.sqrt(3./5.)*R
    B = np.pi*1.62/(Rg**4)    #are we missing *np.pi? 1.62 = 6*(3/5)**2/(4/3) sense?
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg/np.sqrt(6)))**3
    return np.sqrt(np.exp((-Q[:,np.newaxis]**2*Rg**2)/3.)+(B/QstV**4))
    
def UniRodFF(Q,R,args):
    ''' Compute form factor for unified rod - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float R]: cylinder radius (A)
    returns float: form factor
    '''
    L = args[0]
    Rg2 = np.sqrt(R**2/2+L**2/12)
    B2 = np.pi/L
    Rg1 = np.sqrt(3.)*R/2.
    G1 = (2./3.)*R/L
    B1 = 4.*(L+R)/(R**3*L**2)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg2/np.sqrt(6)))**3
    FF = np.exp(-Q[:,np.newaxis]**2*Rg2**2/3.)+(B2/QstV)*np.exp(-Rg1**2*Q[:,np.newaxis]**2/3.)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg1/np.sqrt(6)))**3
    FF += G1*np.exp(-Q[:,np.newaxis]**2*Rg1**2/3.)+(B1/QstV**4)
    return np.sqrt(FF)
    
def UniRodARFF(Q,R,args):
    ''' Compute form factor for unified rod of fixed aspect ratio - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float AR]: cylinder aspect ratio = L/D = L/2R
    returns float: form factor
    '''
    AR = args[0]
    return UniRodFF(Q,R,args=[2.*AR*R,])
    
def UniDiskFF(Q,R,args):
    ''' Compute form factor for unified disk - can use numpy arrays
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float T]: disk thickness (A)
    returns float: form factor
    '''
    T = args[0]
    Rg2 = np.sqrt(R**2/2.+T**2/12.)
    B2 = 2./R**2
    Rg1 = np.sqrt(3.)*T/2.
    RgC2 = 1.1*Rg1
    G1 = (2./3.)*(T/R)**2
    B1 = 4.*(T+R)/(R**3*T**2)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg2/np.sqrt(6)))**3
    FF = np.exp(-Q[:,np.newaxis]**2*Rg2**2/3.)+(B2/QstV**2)*np.exp(-RgC2**2*Q[:,np.newaxis]**2/3.)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg1/np.sqrt(6)))**3
    FF += G1*np.exp(-Q[:,np.newaxis]**2*Rg1**2/3.)+(B1/QstV**4)
    return np.sqrt(FF)
    
def UniTubeFF(Q,R,args):
    ''' Compute form factor for unified tube - can use numpy arrays
    assumes that core of tube is same as the matrix/solvent so contrast
    is from tube wall vs matrix
    param float Q: Q value array (A-1)
    param float R: cylinder radius (A)
    param array args: [float L,T]: tube length & wall thickness(A)
    returns float: form factor
    '''
    L,T = args[:2]
    Ri = R-T
    DR2 = R**2-Ri**2
    Vt = np.pi*DR2*L
    Rg3 = np.sqrt(DR2/2.+L**2/12.)
    B1 = 4.*np.pi**2*(DR2+L*(R+Ri))/Vt**2
    B2 = np.pi**2*T/Vt
    B3 = np.pi/L
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*Rg3/np.sqrt(6)))**3
    FF = np.exp(-Q[:,np.newaxis]**2*Rg3**2/3.)+(B3/QstV)*np.exp(-Q[:,np.newaxis]**2*R**2/3.)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*R/np.sqrt(6)))**3
    FF += (B2/QstV**2)*np.exp(-Q[:,np.newaxis]**2*T**2/3.)
    QstV = Q[:,np.newaxis]/(scsp.erf(Q[:,np.newaxis]*T/np.sqrt(6)))**3
    FF += B1/QstV**4
    return np.sqrt(FF)

###############################################################################
#### Particle volumes
###############################################################################

def SphereVol(R,args=()):
    ''' Compute volume of sphere
    - numpy array friendly
    param float R: sphere radius
    param array args: ignored
    returns float: volume
    '''
    return (4./3.)*np.pi*R**3

def SphericalShellVol(R,args):
    ''' Compute volume of spherical shell
    - numpy array friendly
    param float R: sphere radius
    param array args: [float r]: controls shell thickness, see SphericalShellFF description
    returns float: volume
    '''
    r = args[0]
    if r < 0:
        r = 0
    if r < 1:
        return SphereVol(R) - SphereVol(R*r)
    else:
        return SphereVol(R*r) - SphereVol(R)

def SpheroidVol(R,args):
    ''' Compute volume of cylindrically symmetric ellipsoid (spheroid) 
    - numpy array friendly
    param float R: radius along 2 axes of spheroid
    param array args: [float AR]: aspect ratio so radius of 3rd axis = R*AR
    returns float: volume
    '''
    AR = args[0]
    return AR*SphereVol(R)
    
def CylinderVol(R,args):
    ''' Compute cylinder volume for radius & length
    - numpy array friendly
    param float R: diameter (A)
    param array args: [float L]: length (A)
    returns float:volume (A^3)
    '''
    L = args[0]
    return np.pi*L*R**2
    
def CylinderDVol(L,args):
    ''' Compute cylinder volume for length & diameter
    - numpy array friendly
    param float: L half length (A)
    param array args: [float D]: diameter (A)
    returns float:volume (A^3)
    '''
    D = args[0]
    return CylinderVol(D/2.,[2.*L,])
    
def CylinderARVol(R,args):
    ''' Compute cylinder volume for radius & aspect ratio = L/D
    - numpy array friendly
    param float: R radius (A)
    param array args: [float AR]: =L/D=L/2R aspect ratio
    returns float:volume
    '''
    AR = args[0]
    return CylinderVol(R,[2.*R*AR,])
    
def UniSphereVol(R,args=()):
    ''' Compute volume of sphere
    - numpy array friendly
    param float R: sphere radius
    param array args: ignored
    returns float: volume
    '''
    return SphereVol(R)
    
def UniRodVol(R,args):
    ''' Compute cylinder volume for radius & length
    - numpy array friendly
    param float R: diameter (A)
    param array args: [float L]: length (A)
    returns float:volume (A^3)
    '''
    L = args[0]
    return CylinderVol(R,[L,])
    
def UniRodARVol(R,args):
    ''' Compute rod volume for radius & aspect ratio
    - numpy array friendly
    param float R: diameter (A)
    param array args: [float AR]: =L/D=L/2R aspect ratio
    returns float:volume (A^3)
    '''
    AR = args[0]
    return CylinderARVol(R,[AR,])
    
def UniDiskVol(R,args):
    ''' Compute disk volume for radius & thickness
    - numpy array friendly
    param float R: diameter (A)
    param array args: [float T]: thickness
    returns float:volume (A^3)
    '''
    T = args[0]
    return CylinderVol(R,[T,])
    
def UniTubeVol(R,args):
    ''' Compute tube volume for radius, length & wall thickness
    - numpy array friendly
    param float R: diameter (A)
    param array args: [float L,T]: tube length & wall thickness(A)
    returns float: volume (A^3) of tube wall
    '''
    L,T = args[:2]
    return CylinderVol(R,[L,])-CylinderVol(R-T,[L,])
    
################################################################################
#### Distribution functions & their cumulative fxns
################################################################################

def LogNormalDist(x,pos,args):
    ''' Standard LogNormal distribution - numpy friendly on x axis
    ref: http://www.itl.nist.gov/div898/handbook/index.htm 1.3.6.6.9
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (m)
    param float shape: shape - (sigma of log(LogNormal))
    returns float: LogNormal distribution
    '''
    scale,shape = args
    return np.exp(-np.log((x-pos)/scale)**2/(2.*shape**2))/(np.sqrt(2.*np.pi)*(x-pos)*shape)
    
def GaussDist(x,pos,args):
    ''' Standard Normal distribution - numpy friendly on x axis
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (sigma)
    param float shape: not used
    returns float: Normal distribution
    '''
    scale = args[0]
    return (1./(scale*np.sqrt(2.*np.pi)))*np.exp(-(x-pos)**2/(2.*scale**2))
    
def LSWDist(x,pos,args=[]):
    ''' Lifshitz-Slyozov-Wagner Ostwald ripening distribution - numpy friendly on x axis
    ref: 
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: not used
    param float shape: not used
    returns float: LSW distribution    
    '''
    redX = x/pos
    result = (81.*2**(-5/3.))*(redX**2*np.exp(-redX/(1.5-redX)))/((1.5-redX)**(11/3.)*(3.-redX)**(7/3.))
    
    return np.nan_to_num(result/pos)
    
def SchulzZimmDist(x,pos,args):
    ''' Schulz-Zimm macromolecule distribution - numpy friendly on x axis
    ref: http://goldbook.iupac.org/S05502.html
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (sigma)
    param float shape: not used
    returns float: Schulz-Zimm distribution
    '''
    scale = args[0]
    b = (2.*pos/scale)**2
    a = b/pos
    if b<70:    #why bother?
        return (a**(b+1.))*x**b*np.exp(-a*x)/scsp.gamma(b+1.)
    else:
        return np.exp((b+1.)*np.log(a)-scsp.gammaln(b+1.)+b*np.log(x)-(a*x))
           
def LogNormalCume(x,pos,args):
    ''' Standard LogNormal cumulative distribution - numpy friendly on x axis
    ref: http://www.itl.nist.gov/div898/handbook/index.htm 1.3.6.6.9
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (sigma)
    param float shape: shape parameter
    returns float: LogNormal cumulative distribution
    '''
    scale,shape = args
    lredX = np.log((x-pos)/scale)
    return (scsp.erf((lredX/shape)/np.sqrt(2.))+1.)/2.
    
def GaussCume(x,pos,args):
    ''' Standard Normal cumulative distribution - numpy friendly on x axis
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (sigma)
    param float shape: not used
    returns float: Normal cumulative distribution
    '''
    scale = args[0]
    redX = (x-pos)/scale
    return (scsp.erf(redX/np.sqrt(2.))+1.)/2.
    
def LSWCume(x,pos,args=[]):
    ''' Lifshitz-Slyozov-Wagner Ostwald ripening cumulative distribution - numpy friendly on x axis
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: not used
    param float shape: not used
    returns float: LSW cumulative distribution
    '''
    nP = 500
    xMin,xMax = [0.,2.*pos]
    X = np.linspace(xMin,xMax,nP,True)
    fxn = LSWDist(X,pos)
    mat = np.outer(np.ones(nP),fxn)
    cume = np.sum(np.tril(mat),axis=1)/np.sum(fxn)
    return np.interp(x,X,cume,0,1)
    
def SchulzZimmCume(x,pos,args):
    ''' Schulz-Zimm cumulative distribution - numpy friendly on x axis
    param float x: independent axis (can be numpy array)
    param float pos: location of distribution
    param float scale: width of distribution (sigma)
    param float shape: not used
    returns float: Normal distribution
    '''
    scale = args[0]
    nP = 500
    xMin = np.fmax(0.,pos-4.*scale)
    xMax = np.fmin(pos+4.*scale,1.e5)
    X = np.linspace(xMin,xMax,nP,True)
    fxn = LSWDist(X,pos)
    mat = np.outer(np.ones(nP),fxn)
    cume = np.sum(np.tril(mat),axis=1)/np.sum(fxn)
    return np.interp(x,X,cume,0,1)
    
    return []
    
################################################################################
#### Structure factors for condensed systems
################################################################################

def DiluteSF(Q,args=[]):
    ''' Default: no structure factor correction for dilute system
    '''
    return np.ones_like(Q)  #or return 1.0

def HardSpheresSF(Q,args):
    ''' Computes structure factor for not dilute monodisperse hard spheres
    Refs.: PERCUS,YEVICK PHYS. REV. 110 1 (1958),THIELE J. CHEM PHYS. 39 474 (1968),
    WERTHEIM  PHYS. REV. LETT. 47 1462 (1981)
    
    param float Q: Q value array (A-1)
    param array args: [float R, float VolFrac]: interparticle distance & volume fraction
    returns numpy array S(Q)
    '''
    
    R,VolFr = args
    denom = (1.-VolFr)**4
    num = (1.+2.*VolFr)**2
    alp = num/denom
    bet = -6.*VolFr*(1.+VolFr/2.)**2/denom
    gamm = 0.5*VolFr*num/denom        
    A = 2.*Q*R
    A2 = A**2
    A3 = A**3
    A4 = A**4
    Rca = np.cos(A)
    Rsa = np.sin(A)
    calp = alp*(Rsa/A2-Rca/A)
    cbet = bet*(2.*Rsa/A2-(A2-2.)*Rca/A3-2./A3)
    cgam = gamm*(-Rca/A+(4./A)*((3.*A2-6.)*Rca/A4+(A2-6.)*Rsa/A3+6./A4))
    pfac = -24.*VolFr/A
    C = pfac*(calp+cbet+cgam)
    return 1./(1.-C)
        
def SquareWellSF(Q,args):
    '''Computes structure factor for not dilute monodisperse hard sphere with a
    square well potential interaction. 
    Refs.: SHARMA,SHARMA, PHYSICA 89A,(1977),213-
    
    :param float Q: Q value array (A-1)
    :param array args: [float R, float VolFrac, float depth, float width]: 
        interparticle distance, volume fraction (<0.08), well depth (e/kT<1.5kT),
        well width
    :returns: numpy array S(Q)
      well depth > 0 attractive & values outside above limits nonphysical cf. 
      Monte Carlo simulations 
    '''
    R,VolFr,Depth,Width = args 
    eta = VolFr
    eta2 = eta*eta
    eta3 = eta*eta2
    eta4 = eta*eta3       
    etam1 = 1. - eta 
    etam14 = etam1**4
    alp = (  ( (1. + 2.*eta)**2 ) + eta3*( eta-4.0 )  )/etam14
    bet = -(eta/3.0) * ( 18. + 20.*eta - 12.*eta2 + eta4 )/etam14
    gam = 0.5*eta*( (1. + 2.*eta)**2 + eta3*(eta-4.) )/etam14

    SK = 2.*Q*R
    SK2 = SK*SK
    SK3 = SK*SK2
    SK4 = SK3*SK
    T1 = alp * SK3 * ( np.sin(SK) - SK * np.cos(SK) )
    T2 = bet * SK2 * ( 2.*SK*np.sin(SK) - (SK2-2.)*np.cos(SK) - 2.0 )
    T3 =   ( 4.0*SK3 - 24.*SK ) * np.sin(SK)  
    T3 = T3 - ( SK4 - 12.0*SK2 + 24.0 )*np.cos(SK) + 24.0    
    T3 = gam*T3
    T4 = -Depth*SK3*(np.sin(Width*SK) - Width*SK*np.cos(Width*SK)+ SK*np.cos(SK) - np.sin(SK) )
    CK =  -24.0*eta*( T1 + T2 + T3 + T4 )/SK3/SK3
    return 1./(1.-CK)
    
def StickyHardSpheresSF(Q,args):
    ''' Computes structure factor for not dilute monodisperse hard spheres
    Refs.: PERCUS,YEVICK PHYS. REV. 110 1 (1958),THIELE J. CHEM PHYS. 39 474 (1968),
    WERTHEIM  PHYS. REV. LETT. 47 1462 (1981)
    
    param float Q: Q value array (A-1)
    param array args: [float R, float VolFrac]: sphere radius & volume fraction
    returns numpy array S(Q)
    '''
    R,VolFr,epis,sticky = args
    eta = VolFr/(1.0-epis)/(1.0-epis)/(1.0-epis)	
    sig = 2.0 * R
    aa = sig/(1.0 - epis)
    etam1 = 1.0 - eta
#  SOLVE QUADRATIC FOR LAMBDA
    qa = eta/12.0
    qb = -1.0*(sticky + eta/(etam1))
    qc = (1.0 + eta/2.0)/(etam1*etam1)
    radic = qb*qb - 4.0*qa*qc
#   KEEP THE SMALLER ROOT, THE LARGER ONE IS UNPHYSICAL
    lam1 = (-1.0*qb-np.sqrt(radic))/(2.0*qa)
    lam2 = (-1.0*qb+np.sqrt(radic))/(2.0*qa)
    lam = min(lam1,lam2)
    mu = lam*eta*etam1
    alp = (1.0 + 2.0*eta - mu)/(etam1*etam1)
    bet = (mu - 3.0*eta)/(2.0*etam1*etam1)
#   CALCULATE THE STRUCTURE FACTOR<P></P>
    kk = Q*aa
    k2 = kk*kk
    k3 = kk*k2
    ds = np.sin(kk)
    dc = np.cos(kk)
    aq1 = ((ds - kk*dc)*alp)/k3
    aq2 = (bet*(1.0-dc))/k2
    aq3 = (lam*ds)/(12.0*kk)
    aq = 1.0 + 12.0*eta*(aq1+aq2-aq3)

    bq1 = alp*(0.5/kk - ds/k2 + (1.0 - dc)/k3)
    bq2 = bet*(1.0/kk - ds/k2)
    bq3 = (lam/12.0)*((1.0 - dc)/kk)
    bq = 12.0*eta*(bq1+bq2-bq3)
    sq = 1.0/(aq*aq +bq*bq)

    return sq

#def HayterPenfoldSF(Q,args): #big & ugly function - do later (if ever)
#    pass
    
def InterPrecipitateSF(Q,args):
    ''' Computes structure factor for precipitates in a matrix
    Refs.: E-Wen Huang, Peter K. Liaw, Lionel Porcar, Yun Liu, Yee-Lang Liu, 
    Ji-Jung Kai, and Wei-Ren Chen,APPLIED PHYSICS LETTERS 93, 161904 (2008)
    R. Giordano, A. Grasso, and J. Teixeira, Phys. Rev. A 43, 6894    1991    
    param float Q: Q value array (A-1)
    param array args: [float R, float VolFr]: "radius" & volume fraction
    returns numpy array S(Q)
    '''
    R,VolFr = args
    QV2 = Q**2*VolFr**2
    top = 1 - np.exp(-QV2/4)*np.cos(2.*Q*R)
    bot = 1-2*np.exp(-QV2/4)*np.cos(2.*Q*R) + np.exp(-QV2/2)
    return 2*(top/bot) - 1
          
################################################################################
##### SB-MaxEnt
################################################################################

def G_matrix(q,r,contrast,FFfxn,Volfxn,args=()):
    '''Calculates the response matrix :math:`G(Q,r)` 
    
    :param float q: :math:`Q`
    :param float r: :math:`r`
    :param float contrast: :math:`|\\Delta\\rho|^2`, the scattering contrast
    :param function FFfxn: form factor function FF(q,r,args)
    :param function Volfxn: volume function Vol(r,args)
    :returns float: G(Q,r)
    '''
    FF = FFfxn(q,r,args)
    Vol = Volfxn(r,args)
    return 1.e-4*(contrast*Vol*FF**2).T     #10^-20 vs 10^-24
    
'''
sbmaxent

Entropy maximization routine as described in the article
J Skilling and RK Bryan; MNRAS 211 (1984) 111 - 124.
("MNRAS": "Monthly Notices of the Royal Astronomical Society")

:license: Copyright (c) 2013, UChicago Argonne, LLC
:license: This file is distributed subject to a Software License Agreement found
     in the file LICENSE that is included with this distribution. 

References:

1. J Skilling and RK Bryan; MON NOT R ASTR SOC 211 (1984) 111 - 124.
2. JA Potton, GJ Daniell, and BD Rainford; Proc. Workshop
   Neutron Scattering Data Analysis, Rutherford
   Appleton Laboratory, UK, 1986; ed. MW Johnson,
   IOP Conference Series 81 (1986) 81 - 86, Institute
   of Physics, Bristol, UK.
3. ID Culverwell and GP Clarke; Ibid. 87 - 96.
4. JA Potton, GK Daniell, & BD Rainford,
   J APPL CRYST 21 (1988) 663 - 668.
5. JA Potton, GJ Daniell, & BD Rainford,
   J APPL CRYST 21 (1988) 891 - 897.

'''

class MaxEntException(Exception): 
    '''Any exception from this module'''
    pass

def MaxEnt_SB(datum, sigma, G, base, IterMax, image_to_data=None, data_to_image=None, report=False):
    '''
    do the complete Maximum Entropy algorithm of Skilling and Bryan
    
    :param float datum[]:
    :param float sigma[]:
    :param float[][] G: transformation matrix
    :param float base[]:
    :param int IterMax:
    :param obj image_to_data: opus function (defaults to opus)
    :param obj data_to_image: tropus function (defaults to tropus)
    
    :returns float[]: :math:`f(r) dr`
    '''
        
    TEST_LIMIT        = 0.05                    # for convergence
    CHI_SQR_LIMIT     = 0.01                    # maximum difference in ChiSqr for a solution
    SEARCH_DIRECTIONS = 3                       # <10.  This code requires value = 3
    RESET_STRAYS      = 1                       # was 0.001, correction of stray negative values
    DISTANCE_LIMIT_FACTOR = 0.1                 # limitation on df to constrain runaways
    
    MAX_MOVE_LOOPS    = 5000                     # for no solution in routine: move, 
    MOVE_PASSES       = 0.001                   # convergence test in routine: move

    def tropus (data, G):
        '''
        tropus: transform data-space -> solution-space:  [G] * data
        
        default definition, caller can use this definition or provide an alternative
        
        :param float[M] data: observations, ndarray of shape (M)
        :param float[M][N] G: transformation matrix, ndarray of shape (M,N)
        :returns float[N]: calculated image, ndarray of shape (N)
        '''
        return G.dot(data)

    def opus (image, G):
        '''
        opus: transform solution-space -> data-space:  [G]^tr * image
        
        default definition, caller can use this definition or provide an alternative
        
        :param float[N] image: solution, ndarray of shape (N)
        :param float[M][N] G: transformation matrix, ndarray of shape (M,N)
        :returns float[M]: calculated data, ndarray of shape (M)
        '''
        return np.dot(G.T,image)    #G.transpose().dot(image)

    def Dist(s2, beta):
        '''measure the distance of this possible solution'''
        w = 0
        n = beta.shape[0]
        for k in range(n):
            z = -sum(s2[k] * beta)
            w += beta[k] * z
        return w
    
    def ChiNow(ax, c1, c2, s1, s2):
        '''
        ChiNow
        
        :returns tuple: (ChiNow computation of ``w``, beta)
        '''
        
        bx = 1 - ax
        a =   bx * c2  -  ax * s2
        b = -(bx * c1  -  ax * s1)
    
        beta = ChoSol(a, b)
        w = 1.0
        for k in range(SEARCH_DIRECTIONS):
            w += beta[k] * (c1[k] + 0.5*sum(c2[k] * beta))
        return w, beta
    
    def ChoSol(a, b):
        '''
        ChoSol: ? chop the solution vectors ?
        
        :returns: new vector beta
        '''
        n = b.shape[0]
        fl = np.zeros((n,n))
        bl = np.zeros_like(b)
        
        #print_arr("ChoSol: a", a)
        #print_vec("ChoSol: b", b)
    
        if (a[0][0] <= 0):
            msg = "ChoSol: a[0][0] = " 
            msg += str(a[0][0])
            msg += '  Value must be positive'
            raise MaxEntException(msg)
    
        # first, compute fl from a
        # note fl is a lower triangular matrix
        fl[0][0] = math.sqrt (a[0][0])
        for i in (1, 2):
            fl[i][0] = a[i][0] / fl[0][0]
            for j in range(1, i+1):
                z = 0.0
                for k in range(j):
                    z += fl[i][k] * fl[j][k]
                    #print "ChoSol: %d %d %d  z = %lg" % ( i, j, k, z)
                z = a[i][j] - z
                if j == i:
                    y = math.sqrt(max(0.,z))
                else:
                    y = z / fl[j][j]
                fl[i][j] = y
        #print_arr("ChoSol: fl", fl)
    
        # next, compute bl from fl and b
        bl[0] = b[0] / fl[0][0]
        for i in (1, 2):
            z = 0.0
            for k in range(i):
                z += fl[i][k] * bl[k]
                #print "\t", i, k, z
            bl[i] = (b[i] - z) / fl[i][i]
        #print_vec("ChoSol: bl", bl)
    
        # last, compute beta from bl and fl
        beta = np.empty((n))
        beta[-1] = bl[-1] / fl[-1][-1]
        for i in (1, 0):
            z = 0.0
            for k in range(i+1, n):
                z += fl[k][i] * beta[k]
                #print "\t\t", i, k, 'z=', z
            beta[i] = (bl[i] - z) / fl[i][i]
        #print_vec("ChoSol: beta", beta)
    
        return beta

    def MaxEntMove(fSum, blank, chisq, chizer, c1, c2, s1, s2):
        '''
        move beta one step closer towards the solution
        '''
        a_lower, a_upper = 0., 1.          # bracket  "a"
        cmin, beta = ChiNow (a_lower, c1, c2, s1, s2)
        #print "MaxEntMove: cmin = %g" % cmin
        if cmin*chisq > chizer:
            ctarg = (1.0 + cmin)/2
        else:
            ctarg = chizer/chisq
        f_lower = cmin - ctarg
        c_upper, beta = ChiNow (a_upper, c1, c2, s1, s2)
        f_upper = c_upper - ctarg
    
        fx = 2*MOVE_PASSES      # just to start off
        loop = 1
        while abs(fx) >= MOVE_PASSES and loop <= MAX_MOVE_LOOPS:
            a_new = (a_lower + a_upper) * 0.5           # search by bisection
            c_new, beta = ChiNow (a_new, c1, c2, s1, s2)
            fx = c_new - ctarg
            # tighten the search range for the next pass
            if f_lower*fx > 0:
                a_lower, f_lower = a_new, fx
            if f_upper*fx > 0:
                a_upper, f_upper = a_new, fx
            loop += 1
    
        if abs(fx) >= MOVE_PASSES or loop > MAX_MOVE_LOOPS:
            msg = "MaxEntMove: Loop counter = " 
            msg += str(MAX_MOVE_LOOPS)
            msg += '  No convergence in alpha chop'
            raise MaxEntException(msg)
    
        w = Dist (s2, beta);
        m = SEARCH_DIRECTIONS
        if (w > DISTANCE_LIMIT_FACTOR*fSum/blank):        # invoke the distance penalty, SB eq. 17
            for k in range(m):
                beta[k] *= math.sqrt (fSum/(blank*w))
        chtarg = ctarg * chisq
        return w, chtarg, loop, a_new, fx, beta
        
#MaxEnt_SB starts here    

    if image_to_data == None:
        image_to_data = opus
    if data_to_image == None:
        data_to_image = tropus
    n   = len(base)
    npt = len(datum)

    # Note that the order of subscripts for
    # "xi" and "eta" has been reversed from
    # the convention used in the FORTRAN version
    # to enable parts of them to be passed as
    # as vectors to "image_to_data" and "data_to_image".
    xi      = np.zeros((SEARCH_DIRECTIONS, n))
    eta     = np.zeros((SEARCH_DIRECTIONS, npt))
    beta    = np.zeros((SEARCH_DIRECTIONS))
    # s1      = np.zeros((SEARCH_DIRECTIONS))
    # c1      = np.zeros((SEARCH_DIRECTIONS))
    s2      = np.zeros((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))
    c2      = np.zeros((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))

    # TODO: replace blank (scalar) with base (vector)
    blank = sum(base) / len(base)   # use the average value of base

    chizer, chtarg = npt*1.0, npt*1.0
    f = base * 1.0                              # starting distribution is base

    fSum  = sum(f)                              # find the sum of the f-vector
    z = (datum - image_to_data (f, G)) / sigma  # standardized residuals, SB eq. 3
    chisq = sum(z*z)                            # Chi^2, SB eq. 4

    for iter in range(IterMax):
        ox = -2 * z / sigma                        # gradient of Chi^2

        cgrad = data_to_image (ox, G)              # cgrad[i] = del(C)/del(f[i]), SB eq. 8
        sgrad = -np.log(f/base) / (blank*math.exp (1.0))  # sgrad[i] = del(S)/del(f[i])
        snorm = math.sqrt(sum(f * sgrad*sgrad))    # entropy term, SB eq. 22
        cnorm = math.sqrt(sum(f * cgrad*cgrad))    # ChiSqr term, SB eq. 22
        tnorm = sum(f * sgrad * cgrad)             # norm for gradient term TEST 

        a = 1.0
        b = 1.0 / cnorm
        if iter == 0:
            test = 0.0     # mismatch between entropy and ChiSquared gradients
        else:
            test = math.sqrt( ( 1.0 - tnorm/(snorm*cnorm) )/2 ) # SB eq. 37?
            a = 0.5 / (snorm * test)
            b *= 0.5 / test
        xi[0] = f * cgrad / cnorm
        xi[1] = f * (a * sgrad - b * cgrad)

        eta[0] = image_to_data (xi[0], G);          # image --> data
        eta[1] = image_to_data (xi[1], G);          # image --> data
        ox = eta[1] / (sigma * sigma)
        xi[2] = data_to_image (ox, G);              # data --> image
        a = 1.0 / math.sqrt(sum(f * xi[2]*xi[2]))
        xi[2] = f * xi[2] * a
        eta[2] = image_to_data (xi[2], G)           # image --> data
        
#         print_arr("MaxEnt: eta.transpose()", eta.transpose())
#         print_arr("MaxEnt: xi.transpose()", xi.transpose())

        # prepare the search directions for the conjugate gradient technique
        c1 = xi.dot(cgrad) / chisq		            # C_mu, SB eq. 24
        s1 = xi.dot(sgrad)                          # S_mu, SB eq. 24
#         print_vec("MaxEnt: c1", c1)
#         print_vec("MaxEnt: s1", s1)

        for k in range(SEARCH_DIRECTIONS):
            for l in range(k+1):
                c2[k][l] = 2 * sum(eta[k] * eta[l] / sigma/sigma) / chisq
                s2[k][l] = -sum(xi[k] * xi[l] / f) / blank
#         print_arr("MaxEnt: c2", c2)
#         print_arr("MaxEnt: s2", s2)

        # reflect across the body diagonal
        for k, l in ((0,1), (0,2), (1,2)):
            c2[k][l] = c2[l][k] 		    #  M_(mu,nu)
            s2[k][l] = s2[l][k] 		    #  g_(mu,nu)
 
        beta[0] = -0.5 * c1[0] / c2[0][0]
        beta[1] = 0.0
        beta[2] = 0.0
        if (iter > 0):
            w, chtarg, loop, a_new, fx, beta = MaxEntMove(fSum, blank, chisq, chizer, c1, c2, s1, s2)

        f_old = f.copy()    # preserve the last image
        f += xi.transpose().dot(beta)   # move the image towards the solution, SB eq. 25
        
        # As mentioned at the top of p.119,
        # need to protect against stray negative values.
        # In this case, set them to RESET_STRAYS * base[i]
        #f = f.clip(RESET_STRAYS * blank, f.max())
        for i in range(n):
            if f[i] <= 0.0:
                f[i] = RESET_STRAYS * base[i]
        df = f - f_old
        fSum = sum(f)
        fChange = sum(df)

        # calculate the normalized entropy
        S = sum((f/fSum) * np.log(f/fSum))      # normalized entropy, S&B eq. 1
        z = (datum - image_to_data (f, G)) / sigma  # standardized residuals
        chisq = sum(z*z)                            # report this ChiSq

        if report:
            print (" MaxEnt trial/max: %3d/%3d" % ((iter+1), IterMax))
            print (" Residual: %5.2lf%% Entropy: %8lg" % (100*test, S))
            print (" Function sum: %.6lg Change from last: %.2lf%%\n" % (fSum,100*fChange/fSum))

        # See if we have finished our task.
        # do the hardest test first
        if (abs(chisq/chizer-1.0) < CHI_SQR_LIMIT) and  (test < TEST_LIMIT):
            print (' Convergence achieved.')
            return chisq,f,image_to_data(f, G)     # solution FOUND returns here
    print (' No convergence! Try increasing Error multiplier.')
    return chisq,f,image_to_data(f, G)       # no solution after IterMax iterations

    
###############################################################################
#### IPG/TNNLS Routines
###############################################################################

def IPG(datum,sigma,G,Bins,Dbins,IterMax,Qvec=[],approach=0.8,Power=-1,report=False):
    ''' An implementation of the Interior-Point Gradient method of 
    Michael Merritt & Yin Zhang, Technical Report TR04-08, Dept. of Comp. and 
    Appl. Math., Rice Univ., Houston, Texas 77005, U.S.A. found on the web at
    http://www.caam.rice.edu/caam/trs/2004/TR04-08.pdf
    Problem addressed: Total Non-Negative Least Squares (TNNLS)
    :param float datum[]:
    :param float sigma[]:
    :param float[][] G: transformation matrix
    :param int IterMax:
    :param float Qvec: data positions for Power = 0-4
    :param float approach: 0.8 default fitting parameter
    :param int Power: 0-4 for Q^Power weighting, -1 to use input sigma
    
    '''
    if Power < 0: 
        GmatE = G/sigma[:np.newaxis]
        IntE = datum/sigma
        pwr = 0
        QvecP = np.ones_like(datum)
    else:
        GmatE = G[:]
        IntE = datum[:]
        pwr = Power
        QvecP = Qvec**pwr
    Amat = GmatE*QvecP[:np.newaxis]
    AAmat = np.inner(Amat,Amat)
    Bvec = IntE*QvecP
    Xw = np.ones_like(Bins)*1.e-32
    calc = np.dot(G.T,Xw)
    nIter = 0
    err = 10.
    while (nIter<IterMax) and (err > 1.):
        #Step 1 in M&Z paper:
        Qk = np.inner(AAmat,Xw)-np.inner(Amat,Bvec)
        Dk = Xw/np.inner(AAmat,Xw)
        Pk = -Dk*Qk
        #Step 2 in M&Z paper:
        alpSt = -np.inner(Pk,Qk)/np.inner(Pk,np.inner(AAmat,Pk))
        alpWv = -Xw/Pk
        AkSt = [np.where(Pk[i]<0,np.min((approach*alpWv[i],alpSt)),alpSt) for i in range(Pk.shape[0])]
        #Step 3 in M&Z paper:
        shift = AkSt*Pk
        Xw += shift
        #done IPG; now results
        nIter += 1
        calc = np.dot(G.T,Xw)
        chisq = np.sum(((datum-calc)/sigma)**2)
        err = chisq/len(datum)
        if report:
            print (' Iteration: %d, chisq: %.3g, sum(shift^2): %.3g'%(nIter,chisq,np.sum(shift**2)))
    return chisq,Xw,calc

###############################################################################
#### SASD Utilities
###############################################################################

def SetScale(Data,refData):
    Profile,Limits,Sample = Data
    refProfile,refLimits,refSample = refData
    x,y = Profile[:2]
    rx,ry = refProfile[:2]
    Beg = np.fmax([rx[0],x[0],Limits[1][0],refLimits[1][0]])
    Fin = np.fmin([rx[-1],x[-1],Limits[1][1],refLimits[1][1]])
    iBeg = np.searchsorted(x,Beg)
    iFin = np.searchsorted(x,Fin)+1        #include last point
    sum = np.sum(y[iBeg:iFin])
    refsum = np.sum(np.interp(x[iBeg:iFin],rx,ry,0,0))
    Sample['Scale'][0] = refSample['Scale'][0]*refsum/sum
    
def Bestimate(G,Rg,P):
    return (G*P/Rg**P)*np.exp(scsp.gammaln(P/2))
    
def SmearData(Ic,Q,slitLen,Back):
    Np = Q.shape[0]
    Qtemp = np.concatenate([Q,Q[-1]+20*Q])
    Ictemp = np.concatenate([Ic,Ic[-1]*(1-(Qtemp[Np:]-Qtemp[Np])/(20*Qtemp[Np-1]))])
    Icsm = np.zeros_like(Q)
    Qsm = 2*slitLen*(np.interp(np.arange(2*Np)/2.,np.arange(Np),Q)-Q[0])/(Q[-1]-Q[0])
    Sp = np.searchsorted(Qsm,slitLen)
    DQsm = np.diff(Qsm)[:Sp]
    Ism = np.interp(np.sqrt(Q[:,np.newaxis]**2+Qsm**2),Qtemp,Ictemp)
    Icsm = np.sum((Ism[:,:Sp]*DQsm),axis=1)
    Icsm /= slitLen
    return Icsm
    
###############################################################################
#### Size distribution
###############################################################################

def SizeDistribution(Profile,ProfDict,Limits,Sample,data):
    shapes = {'Spheroid':[SpheroidFF,SpheroidVol],'Cylinder':[CylinderDFF,CylinderDVol],
        'Cylinder AR':[CylinderARFF,CylinderARVol],'Unified sphere':[UniSphereFF,UniSphereVol],
        'Unified rod':[UniRodFF,UniRodVol],'Unified rod AR':[UniRodARFF,UniRodARVol],
        'Unified disk':[UniDiskFF,UniDiskVol],'Sphere':[SphereFF,SphereVol],
        'Cylinder diam':[CylinderDFF,CylinderDVol],
        'Spherical shell': [SphericalShellFF, SphericalShellVol]}
    Shape = data['Size']['Shape'][0]
    Parms = data['Size']['Shape'][1:]
    if data['Size']['logBins']:
        Bins = np.logspace(np.log10(data['Size']['MinDiam']),np.log10(data['Size']['MaxDiam']),
            data['Size']['Nbins']+1,True)/2.        #make radii
    else:
        Bins = np.linspace(data['Size']['MinDiam'],data['Size']['MaxDiam'],
            data['Size']['Nbins']+1,True)/2.        #make radii
    Dbins = np.diff(Bins)
    Bins = Bins[:-1]+Dbins/2.
    Contrast = Sample['Contrast'][1]
    Scale = Sample['Scale'][0]
    Sky = 10**data['Size']['MaxEnt']['Sky']
    BinsBack = np.ones_like(Bins)*Sky*Scale/Contrast
    Back = data['Back']
    Q,Io,wt,Ic,Ib = Profile[:5]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    wtFactor = ProfDict['wtFactor']
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)+1        #include last point
    BinMag = np.zeros_like(Bins)
    Ic[:] = 0.
    Gmat = G_matrix(Q[Ibeg:Ifin],Bins,Contrast,shapes[Shape][0],shapes[Shape][1],args=Parms)
    if 'MaxEnt' == data['Size']['Method']:
        chisq,BinMag,Ic[Ibeg:Ifin] = MaxEnt_SB(Scale*Io[Ibeg:Ifin]-Back[0],
            Scale/np.sqrt(wtFactor*wt[Ibeg:Ifin]),Gmat,BinsBack,
            data['Size']['MaxEnt']['Niter'],report=True)
    elif 'IPG' == data['Size']['Method']:
        chisq,BinMag,Ic[Ibeg:Ifin] = IPG(Scale*Io[Ibeg:Ifin]-Back[0],Scale/np.sqrt(wtFactor*wt[Ibeg:Ifin]),
            Gmat,Bins,Dbins,data['Size']['IPG']['Niter'],Q[Ibeg:Ifin],approach=0.8,
            Power=data['Size']['IPG']['Power'],report=True)
    Ib[:] = Back[0]
    Ic[Ibeg:Ifin] += Back[0]
    print (' Final chi^2: %.3f'%(chisq))
    data['Size']['Distribution'] = [Bins,Dbins,BinMag/(2.*Dbins)]
    
################################################################################
#### Modelling
################################################################################

def PairDistFxn(Profile,ProfDict,Limits,Sample,data):
    
    def CalcMoore():
        
        def MoorePOR(cw,r,dmax):    #lines 1417-1428
            n = 0
            nmax = len(cw)
            POR = np.zeros(len(r))
            while n < nmax:
                POR += 4.*r*cw[n]*np.sin((n+1.)*np.pi*r/dmax)
                n += 1
            return POR
        
        def MooreIOREFF(cw,q,dmax):      #lines 1437-1448
            n = 0
            nmax = len(cw)
            POR = np.zeros(len(q))
            dq = dmax*q
            fpd = 8.*(np.pi**2)*dmax*np.sin(dq)/q
            while n < nmax:
                POR += cw[n]*(n+1.)*((-1)**n)*fpd/(((n+1.)*np.pi)**2-dq**2)
                n += 1
            return POR
        
        def calcSASD(values,Q,Io,wt,Ifb,dmax,ifBack):
            if ifBack:
                M = np.sqrt(wt)*(MooreIOREFF(values[:-1],Q,dmax)+Ifb+values[-1]-Io)
            else:
                M = np.sqrt(wt)*(MooreIOREFF(values,Q,dmax)+Ifb-Io)
            return M

        Q,Io,wt,Ic,Ib,Ifb = Profile[:6]
        Qmin = Limits[1][0]
        Qmax = Limits[1][1]
        wtFactor = ProfDict['wtFactor']
        Back,ifBack = data['Back']
        Ibeg = np.searchsorted(Q,Qmin)
        Ifin = np.searchsorted(Q,Qmax)+1    #include last point
        Ic[Ibeg:Ifin] = 0
        Bins = np.linspace(0.,pairData['MaxRadius'],pairData['NBins']+1,True)
        Dbins = np.diff(Bins)
        Bins = Bins[:-1]+Dbins/2.
        N = pairData['Moore']
        if ifBack:
            N += 1
        MPV = np.ones(N)*1.e-5
        dmax = pairData['MaxRadius']
        if 'User' in pairData['Errors']:
            Wt = wt[Ibeg:Ifin]
        elif 'Sqrt' in pairData['Errors']:
            Wt = 1./Io[Ibeg:Ifin]
        elif 'Percent' in pairData['Errors']:
            Wt = 1./(pairData['Percent error']*Io[Ibeg:Ifin])
        result = so.leastsq(calcSASD,MPV,full_output=True,epsfcn=1.e-8,   #ftol=Ftol,
            args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*Wt,Ifb[Ibeg:Ifin],dmax,ifBack))
        if ifBack:
            MPVR = result[0][:-1]
            data['Back'][0] = result[0][-1]
            Back = data['Back'][0]
        else:        
            MPVR = result[0]
            Back = 0.
        chisq = np.sum(result[2]['fvec']**2)
        covM = result[1]
        Ic[Ibeg:Ifin] = MooreIOREFF(MPVR,Q[Ibeg:Ifin],dmax)+Ifb[Ibeg:Ifin]+Back
        ncalc = result[2]['nfev']
        GOF = chisq/(Ifin-Ibeg-N)
        Rwp = np.sqrt(chisq/np.sum(wt[Ibeg:Ifin]*Io[Ibeg:Ifin]**2))*100.      #to %
        print (' Results of small angle data modelling fit of P(R):')
        print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(ncalc,Ifin-Ibeg,N))
        print ('Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rwp,chisq,GOF))
        if len(covM):
            sig = np.sqrt(np.diag(covM)*GOF)
            for val,esd in zip(result[0],sig):
                print(' parameter: %.4g esd: %.4g'%(val,esd))            
        BinMag = MoorePOR(MPVR,Bins,dmax)/2.
        return Bins,Dbins,BinMag
    
    pairData = data['Pair']
    
    if pairData['Method'] == 'Regularization':      #not used
        print('Regularization calc; dummy Gaussian - TBD')
        pairData['Method'] = 'Moore'
        
        
    elif pairData['Method'] == 'Moore':
        Bins,Dbins,BinMag = CalcMoore()
        BinSum = np.sum(BinMag*Dbins)
        BinMag /= BinSum        
    
    data['Pair']['Distribution'] = [Bins,Dbins,BinMag]      #/(2.*Dbins)]
    if 'Pair Calc' in data['Pair']:
        del data['Pair']['Pair Calc']
    
    
        
################################################################################
#### Modelling
################################################################################

def ModelFit(Profile,ProfDict,Limits,Sample,Model):
    shapes = {'Spheroid':[SpheroidFF,SpheroidVol],'Cylinder':[CylinderDFF,CylinderDVol],
        'Cylinder AR':[CylinderARFF,CylinderARVol],'Unified sphere':[UniSphereFF,UniSphereVol],
        'Unified rod':[UniRodFF,UniRodVol],'Unified rod AR':[UniRodARFF,UniRodARVol],
        'Unified disk':[UniDiskFF,UniDiskVol],'Sphere':[SphereFF,SphereVol],
        'Unified tube':[UniTubeFF,UniTubeVol],'Cylinder diam':[CylinderDFF,CylinderDVol],
        'Spherical shell':[SphericalShellFF,SphericalShellVol]}
            
    sfxns = {'Dilute':DiluteSF,'Hard sphere':HardSpheresSF,'Square well':SquareWellSF,
            'Sticky hard sphere':StickyHardSpheresSF,'InterPrecipitate':InterPrecipitateSF,}
                
    parmOrder = ['Volume','Radius','Mean','StdDev','MinSize','G','Rg','B','P','Cutoff',
        'PkInt','PkPos','PkSig','PkGam',]
        
    FFparmOrder = ['Aspect ratio','Length','Diameter','Thickness','Shell thickness']
    
    SFparmOrder = ['Dist','VolFr','epis','Sticky','Depth','Width']

    def GetModelParms():
        parmDict = {'Scale':Sample['Scale'][0],'SlitLen':Sample.get('SlitLen',0.0),}
        varyList = []
        values = []
        levelTypes = []
        Back = Model['Back']
        if Back[1]:
            varyList += ['Back',]
            values.append(Back[0])
        parmDict['Back'] = Back[0]
        partData = Model['Particle']
        for i,level in enumerate(partData['Levels']):
            cid = str(i)+';'
            controls = level['Controls']
            Type = controls['DistType']
            levelTypes.append(Type)
            if Type in ['LogNormal','Gaussian','LSW','Schulz-Zimm','Monodisperse']:
                if Type not in ['Monodosperse',]:
                    parmDict[cid+'NumPoints'] = controls['NumPoints']
                    parmDict[cid+'Cutoff'] = controls['Cutoff']
                parmDict[cid+'FormFact'] = shapes[controls['FormFact']][0]
                parmDict[cid+'FFVolume'] = shapes[controls['FormFact']][1]
                parmDict[cid+'StrFact'] = sfxns[controls['StrFact']]
                parmDict[cid+'Contrast'] = controls['Contrast']
                for item in FFparmOrder:
                    if item in controls['FFargs']:
                        parmDict[cid+item] = controls['FFargs'][item][0]
                        if controls['FFargs'][item][1]:
                            varyList.append(cid+item)
                            values.append(controls['FFargs'][item][0])
                for item in SFparmOrder:
                    if item in controls.get('SFargs',{}):
                        parmDict[cid+item] = controls['SFargs'][item][0]
                        if controls['SFargs'][item][1]:
                            varyList.append(cid+item)
                            values.append(controls['SFargs'][item][0])
            distDict = controls['DistType']
            for item in parmOrder:
                if item in level[distDict]:
                    parmDict[cid+item] = level[distDict][item][0]
                    if level[distDict][item][1]:
                        values.append(level[distDict][item][0])
                        varyList.append(cid+item)
        return levelTypes,parmDict,varyList,values
        
    def SetModelParms():
        print (' Refined parameters: Histogram scale: %.4g'%(parmDict['Scale']))
        if 'Back' in varyList:
            Model['Back'][0] = parmDict['Back']
            print ('  %15s %15.4f esd: %15.4g'%('Background:',parmDict['Back'],sigDict['Back']))
        partData = Model['Particle']
        for i,level in enumerate(partData['Levels']):
            controls = level['Controls']
            Type = controls['DistType']
            if Type in ['LogNormal','Gaussian','LSW','Schulz-Zimm','Monodisperse']:
                print (' Component %d: Type: %s: Structure Factor: %s Contrast: %12.3f'  \
                    %(i,Type,controls['StrFact'],controls['Contrast']))              
            else:
                print (' Component %d: Type: %s: '%(i,Type,))
            cid = str(i)+';'
            if Type in ['LogNormal','Gaussian','LSW','Schulz-Zimm','Monodisperse']:
                for item in FFparmOrder:
                    if cid+item in varyList:
                        controls['FFargs'][item][0] = parmDict[cid+item]
                        print (' %15s: %15.4g esd: %15.4g'%(cid+item,parmDict[cid+item],sigDict[cid+item]))
                for item in SFparmOrder:
                    if cid+item in varyList:
                        controls['SFargs'][item][0] = parmDict[cid+item]
                        print (' %15s: %15.4g esd: %15.4g'%(cid+item,parmDict[cid+item],sigDict[cid+item]))
            distDict = controls['DistType']
            for item in level[distDict]:
                if cid+item in varyList:
                    level[distDict][item][0] = parmDict[cid+item]
                    print (' %15s: %15.4g esd: %15.4g'%(cid+item,parmDict[cid+item],sigDict[cid+item]))
                    
    def calcSASD(values,Q,Io,wt,Ifb,levelTypes,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        M = np.sqrt(wt)*(getSASD(Q,levelTypes,parmDict)+Ifb-parmDict['Scale']*Io)
        return M
        
    def getSASD(Q,levelTypes,parmDict):
        Ic = np.zeros_like(Q)
        for i,Type in enumerate(levelTypes):
            cid = str(i)+';'
            if Type in ['LogNormal','Gaussian','LSW','Schulz-Zimm']:
                FFfxn = parmDict[cid+'FormFact']
                Volfxn = parmDict[cid+'FFVolume']
                SFfxn = parmDict[cid+'StrFact']
                FFargs = []
                SFargs = []
                for item in [cid+'Aspect ratio',cid+'Length',cid+'Thickness',cid+'Diameter',cid+'Shell thickness']:
                    if item in parmDict: 
                        FFargs.append(parmDict[item])
                for item in [cid+'Dist',cid+'VolFr',cid+'epis',cid+'Sticky',cid+'Depth',cid+'Width']:
                    if item in parmDict: 
                        SFargs.append(parmDict[item])
                distDict = {}
                for item in [cid+'Volume',cid+'Mean',cid+'StdDev',cid+'MinSize',]:
                    if item in parmDict:
                        distDict[item.split(';')[1]] = parmDict[item]
                contrast = parmDict[cid+'Contrast']
                rBins,dBins,dist = MakeDiamDist(Type,parmDict[cid+'NumPoints'],parmDict[cid+'Cutoff'],distDict)
                Gmat = G_matrix(Q,rBins,contrast,FFfxn,Volfxn,FFargs).T
                dist *= parmDict[cid+'Volume']
                Ic += np.dot(Gmat,dist)*SFfxn(Q,args=SFargs)
            elif 'Unified' in Type:
                Rg,G,B,P,Rgco = parmDict[cid+'Rg'],parmDict[cid+'G'],parmDict[cid+'B'], \
                    parmDict[cid+'P'],parmDict[cid+'Cutoff']
                Qstar = Q/(scsp.erf(Q*Rg/np.sqrt(6)))**3
                Guin = G*np.exp(-(Q*Rg)**2/3)
                Porod = (B/Qstar**P)*np.exp(-(Q*Rgco)**2/3)
                Ic += Guin+Porod
            elif 'Porod' in Type:
                B,P,Rgco = parmDict[cid+'B'],parmDict[cid+'P'],parmDict[cid+'Cutoff']
                Porod = (B/Q**P)*np.exp(-(Q*Rgco)**2/3)
                Ic += Porod
            elif 'Mono' in Type:
                FFfxn = parmDict[cid+'FormFact']
                Volfxn = parmDict[cid+'FFVolume']
                SFfxn = parmDict[cid+'StrFact']
                FFargs = []
                SFargs = []
                for item in [cid+'Aspect ratio',cid+'Length',cid+'Thickness',cid+'Diameter',cid+'Shell thickness']:
                    if item in parmDict: 
                        FFargs.append(parmDict[item])
                for item in [cid+'Dist',cid+'VolFr',cid+'epis',cid+'Sticky',cid+'Depth',cid+'Width']:
                    if item in parmDict: 
                        SFargs.append(parmDict[item])
                contrast = parmDict[cid+'Contrast']
                R = parmDict[cid+'Radius']
                Gmat = G_matrix(Q,R,contrast,FFfxn,Volfxn,FFargs)             
                Ic += Gmat[0]*parmDict[cid+'Volume']*SFfxn(Q,args=SFargs)
            elif 'Bragg' in Type:
                Ic += parmDict[cid+'PkInt']*G2pwd.getPsVoigt(parmDict[cid+'PkPos'],
                    parmDict[cid+'PkSig'],parmDict[cid+'PkGam'],Q)[0]
        Ic += parmDict['Back']  #/parmDict['Scale']
        slitLen = Sample['SlitLen']
        if slitLen:
            Ic = SmearData(Ic,Q,slitLen,parmDict['Back'])
        return Ic
        
    Q,Io,wt,Ic,Ib,Ifb = Profile[:6]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    wtFactor = ProfDict['wtFactor']
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)+1    #include last point
    Ic[:] = 0
    levelTypes,parmDict,varyList,values = GetModelParms()
    if varyList:
        result = so.leastsq(calcSASD,values,full_output=True,epsfcn=1.e-8,   #ftol=Ftol,
            args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Ifb[Ibeg:Ifin],levelTypes,parmDict,varyList))
        parmDict.update(zip(varyList,result[0]))
        chisq = np.sum(result[2]['fvec']**2)
        ncalc = result[2]['nfev']
        covM = result[1]
    else:   #nothing varied
        M = calcSASD(values,Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Ifb[Ibeg:Ifin],levelTypes,parmDict,varyList)
        chisq = np.sum(M**2)
        ncalc = 0
        covM = []
        sig = []
        sigDict = {}
        result = []
    Rvals = {}
    Rvals['Rwp'] = np.sqrt(chisq/np.sum(wt[Ibeg:Ifin]*Io[Ibeg:Ifin]**2))*100.      #to %
    Rvals['GOF'] = chisq/(Ifin-Ibeg-len(varyList))       #reduced chi^2
    Ic[Ibeg:Ifin] = getSASD(Q[Ibeg:Ifin],levelTypes,parmDict)
    Msg = 'Failed to converge'
    try:
        Nans = np.isnan(result[0])
        if np.any(Nans):
            name = varyList[Nans.nonzero(True)[0]]
            Msg = 'Nan result for '+name+'!'
            raise ValueError
        Negs = np.less_equal(result[0],0.)
        if np.any(Negs):
            name = varyList[Negs.nonzero(True)[0]]
            Msg = 'negative coefficient for '+name+'!'
            raise ValueError
        if len(covM):
            sig = np.sqrt(np.diag(covM)*Rvals['GOF'])
            sigDict = dict(zip(varyList,sig))
        print (' Results of small angle data modelling fit:')
        print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(ncalc,Ifin-Ibeg,len(varyList)))
        print ('Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],chisq,Rvals['GOF']))
        SetModelParms()
        covMatrix = covM*Rvals['GOF']
        return True,result,varyList,sig,Rvals,covMatrix,parmDict,''
    except (ValueError,TypeError):      #when bad LS refinement; covM missing or with nans
        return False,0,0,0,0,0,0,Msg

def getSASDRg(Q,parmDict):
    Ic = np.zeros_like(Q)
    Rg,G,B = parmDict['Rg'],parmDict['G'],parmDict['B']
    Qstar = Q/(scsp.erf(Q*Rg/np.sqrt(6)))**3
    Guin = G*np.exp(-(Q*Rg)**2/3)
    Porod = (B/Qstar**4)        #*np.exp(-(Q*B)**2/3)
    Ic += Guin+Porod+parmDict['Back']
    return Ic
        
def RgFit(Profile,ProfDict,Limits,Sample,Model):
    print('unified fit single Rg to data to estimate a size')
    
    def GetModelParms():
        parmDict = {}
        varyList = []
        values = []
        Back = Model['Back']
        if Back[1]:
            varyList += ['Back',]
            values.append(Back[0])
        parmDict['Back'] = Back[0]
        pairData = Model['Pair']
        parmDict['G'] = pairData.get('Dist G',Io[Ibeg])
        parmDict['Rg'] = pairData['MaxRadius']/2.5
        parmDict['B'] = pairData.get('Dist B',Io[Ifin-1]*Q[Ifin-1]**4)
        varyList  += ['G','Rg','B']
        values += [parmDict['G'],parmDict['Rg'],parmDict['B']]
        return parmDict,varyList,values

    def calcSASD(values,Q,Io,wt,Ifb,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        M = np.sqrt(wt)*(getSASDRg(Q,parmDict)-Io)
        return M
    
    def SetModelParms():
        print (' Refined parameters: ')
        if 'Back' in varyList:
            Model['Back'][0] = parmDict['Back']
            print ('  %15s %15.4f esd: %15.4g'%('Background:',parmDict['Back'],sigDict['Back']))
        pairData = Model['Pair']
        pairData['Dist G'] = parmDict['G']
        pairData['MaxRadius'] = parmDict['Rg']*2.5
        pairData['Dist B'] = parmDict['B']
        for item in varyList:
            print (' %15s: %15.4g esd: %15.4g'%(item,parmDict[item],sigDict[item]))

    Q,Io,wt,Ic,Ib,Ifb = Profile[:6]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    wtFactor = ProfDict['wtFactor']
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)+1    #include last point
    Ic[:] = 0
    pairData = Model['Pair']
    if 'User' in pairData['Errors']:
        Wt = wt[Ibeg:Ifin]
    elif 'Sqrt' in pairData['Errors']:
        Wt = 1./Io[Ibeg:Ifin]
    elif 'Percent' in pairData['Errors']:
        Wt = 1./(pairData['Percent error']*Io[Ibeg:Ifin])
    parmDict,varyList,values = GetModelParms()
    result = so.leastsq(calcSASD,values,full_output=True,epsfcn=1.e-12,factor=0.1,  #ftol=Ftol,
        args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*Wt,Ifb[Ibeg:Ifin],parmDict,varyList))
    parmDict.update(dict(zip(varyList,result[0])))
    chisq = np.sum(result[2]['fvec']**2)
    ncalc = result[2]['nfev']
    covM = result[1]
    Rvals = {}
    Rvals['Rwp'] = np.sqrt(chisq/np.sum(wt[Ibeg:Ifin]*Io[Ibeg:Ifin]**2))*100.      #to %
    Rvals['GOF'] = chisq/(Ifin-Ibeg-len(varyList))       #reduced chi^2
    Ic[Ibeg:Ifin] = getSASDRg(Q[Ibeg:Ifin],parmDict)
    Msg = 'Failed to converge'
    try:
        Nans = np.isnan(result[0])
        if np.any(Nans):
            name = varyList[Nans.nonzero(True)[0]]
            Msg = 'Nan result for '+name+'!'
            raise ValueError
        Negs = np.less_equal(result[0],0.)
        if np.any(Negs):
            name = varyList[Negs.nonzero(True)[0]]
            Msg = 'negative coefficient for '+name+'!'
            raise ValueError
        if len(covM):
            sig = np.sqrt(np.diag(covM)*Rvals['GOF'])
            sigDict = dict(zip(varyList,sig))
        print (' Results of Rg fit:')
        print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(ncalc,Ifin-Ibeg,len(varyList)))
        print ('Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],chisq,Rvals['GOF']))
        SetModelParms()
        covMatrix = covM*Rvals['GOF']
        return True,result,varyList,sig,Rvals,covMatrix,parmDict,''
    except (ValueError,TypeError):      #when bad LS refinement; covM missing or with nans
        return False,0,0,0,0,0,0,Msg

    return [None,]

    
def ModelFxn(Profile,ProfDict,Limits,Sample,sasdData):
    
    shapes = {'Spheroid':[SpheroidFF,SpheroidVol],'Cylinder':[CylinderDFF,CylinderDVol],
        'Cylinder AR':[CylinderARFF,CylinderARVol],'Unified sphere':[UniSphereFF,UniSphereVol],
        'Unified rod':[UniRodFF,UniRodVol],'Unified rod AR':[UniRodARFF,UniRodARVol],
        'Unified disk':[UniDiskFF,UniDiskVol],'Sphere':[SphereFF,SphereVol],
        'Unified tube':[UniTubeFF,UniTubeVol],'Cylinder diam':[CylinderDFF,CylinderDVol],
        'Spherical shell':[SphericalShellFF,SphericalShellVol]}
    sfxns = {'Dilute':DiluteSF,'Hard sphere':HardSpheresSF,'Square well':SquareWellSF,
            'Sticky hard sphere':StickyHardSpheresSF,'InterPrecipitate':InterPrecipitateSF,}

#    pdb.set_trace()
    partData = sasdData['Particle']
    Back = sasdData['Back']
    Q,Io,wt,Ic,Ib,Ifb = Profile[:6]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)+1    #include last point
    Ib[:] = Back[0]
    Ic[:] = 0
    Rbins = []
    Dist = []
    for level in partData['Levels']:
        controls = level['Controls']
        distFxn = controls['DistType']
        if distFxn in ['LogNormal','Gaussian','LSW','Schulz-Zimm']:
            parmDict = level[controls['DistType']]
            FFfxn = shapes[controls['FormFact']][0]
            Volfxn = shapes[controls['FormFact']][1]
            SFfxn = sfxns[controls['StrFact']]
            FFargs = []
            SFargs = []
            for item in ['Dist','VolFr','epis','Sticky','Depth','Width',]:
                if item in controls.get('SFargs',{}):
                    SFargs.append(controls['SFargs'][item][0])
            for item in ['Aspect ratio','Length','Thickness','Diameter','Shell thickness']:
                if item in controls['FFargs']: 
                    FFargs.append(controls['FFargs'][item][0])
            contrast = controls['Contrast']
            distDict = {}
            for item in parmDict:
                distDict[item] = parmDict[item][0]
            rBins,dBins,dist = MakeDiamDist(controls['DistType'],controls['NumPoints'],controls['Cutoff'],distDict)
            Gmat = G_matrix(Q[Ibeg:Ifin],rBins,contrast,FFfxn,Volfxn,FFargs).T
            dist *= level[distFxn]['Volume'][0]
            Ic[Ibeg:Ifin] += np.dot(Gmat,dist)*SFfxn(Q[Ibeg:Ifin],args=SFargs)
            Rbins.append(rBins)
            Dist.append(dist/(4.*dBins))
        elif 'Unified' in distFxn:
            parmDict = level[controls['DistType']]
            Rg,G,B,P,Rgco = parmDict['Rg'][0],parmDict['G'][0],parmDict['B'][0],    \
                parmDict['P'][0],parmDict['Cutoff'][0]
            Qstar = Q[Ibeg:Ifin]/(scsp.erf(Q[Ibeg:Ifin]*Rg/np.sqrt(6)))**3
            Guin = G*np.exp(-(Q[Ibeg:Ifin]*Rg)**2/3)
            Porod = (B/Qstar**P)*np.exp(-(Q[Ibeg:Ifin]*Rgco)**2/3)
            Ic[Ibeg:Ifin] += Guin+Porod
            Rbins.append([])
            Dist.append([])
        elif 'Porod' in distFxn:
            parmDict = level[controls['DistType']]
            B,P,Rgco = parmDict['B'][0],parmDict['P'][0],parmDict['Cutoff'][0]
            Porod = (B/Q[Ibeg:Ifin]**P)*np.exp(-(Q[Ibeg:Ifin]*Rgco)**2/3)
            Ic[Ibeg:Ifin] += Porod
            Rbins.append([])
            Dist.append([])
        elif 'Mono' in distFxn:
            parmDict = level[controls['DistType']]
            R = level[controls['DistType']]['Radius'][0]
            FFfxn = shapes[controls['FormFact']][0]
            Volfxn = shapes[controls['FormFact']][1]
            SFfxn = sfxns[controls['StrFact']]
            FFargs = []
            SFargs = []
            for item in ['Dist','VolFr','epis','Sticky','Depth','Width',]:
                if item in controls.get('SFargs',{}):
                    SFargs.append(controls['SFargs'][item][0])
            for item in ['Aspect ratio','Length','Thickness','Diameter','Shell thickness']:
                if item in controls['FFargs']: 
                    FFargs.append(controls['FFargs'][item][0])
            contrast = controls['Contrast']
            Gmat = G_matrix(Q[Ibeg:Ifin],R,contrast,FFfxn,Volfxn,FFargs)             
            Ic[Ibeg:Ifin] += Gmat[0]*level[distFxn]['Volume'][0]*SFfxn(Q[Ibeg:Ifin],args=SFargs)
            Rbins.append([])
            Dist.append([])
        elif 'Bragg' in distFxn:
            parmDict = level[controls['DistType']]
            Ic[Ibeg:Ifin] += parmDict['PkInt'][0]*G2pwd.getPsVoigt(parmDict['PkPos'][0],
                parmDict['PkSig'][0],parmDict['PkGam'][0],Q[Ibeg:Ifin])[0]
            Rbins.append([])
            Dist.append([])
    Ic[Ibeg:Ifin] += Back[0]
    slitLen = Sample.get('SlitLen',0.)
    if slitLen:
        Ic[Ibeg:Ifin] = SmearData(Ic,Q,slitLen,Back[0])[Ibeg:Ifin]
    sasdData['Size Calc'] = [Rbins,Dist]
    
def MakeDiamDist(DistName,nPoints,cutoff,distDict):
    
    if 'LogNormal' in DistName:
        distFxn = 'LogNormalDist'
        cumeFxn = 'LogNormalCume'
        pos = distDict['MinSize']
        args = [distDict['Mean'],distDict['StdDev']]
        step = 4.*np.sqrt(np.exp(distDict['StdDev']**2)*(np.exp(distDict['StdDev']**2)-1.))
        mode = distDict['MinSize']+distDict['Mean']/np.exp(distDict['StdDev']**2)
        minX = 1. #pos
        maxX = 1.e4 #np.min([mode+nPoints*step*40,1.e5])
    elif 'Gauss' in DistName:
        distFxn = 'GaussDist'
        cumeFxn = 'GaussCume'
        pos = distDict['Mean']
        args = [distDict['StdDev']]
        mode = distDict['Mean']
        minX = np.max([mode-4.*distDict['StdDev'],1.])
        maxX = np.min([mode+4.*distDict['StdDev'],1.e5])
    elif 'LSW' in DistName:
        distFxn = 'LSWDist'
        cumeFxn = 'LSWCume'
        pos = distDict['Mean']
        args = []
        minX,maxX = [0.,2.*pos]
    elif 'Schulz' in DistName:
        distFxn = 'SchulzZimmDist'
        cumeFxn = 'SchulzZimmCume'
        pos = distDict['Mean']
        args = [distDict['StdDev']]
        minX = np.max([1.,pos-4.*distDict['StdDev']])
        maxX = np.min([pos+4.*distDict['StdDev'],1.e5])
    nP = 500
    Diam = np.logspace(0.,5.,nP,True)
    TCW = eval(cumeFxn+'(Diam,pos,args)')
    CumeTgts = np.linspace(cutoff,(1.-cutoff),nPoints+1,True)
    rBins = np.interp(CumeTgts,TCW,Diam,0,0)
    dBins = np.diff(rBins)
    rBins = rBins[:-1]+dBins/2.
    return rBins,dBins,eval(distFxn+'(rBins,pos,args)')

################################################################################
#### MaxEnt testing stuff
################################################################################

def print_vec(text, a):
    '''print the contents of a vector to the console'''
    n = a.shape[0]
    print ("%s[ = (" % text,end='')
    for i in range(n):
        s = " %g, " % a[i]
        print (s,end='')
    print (")")

def print_arr(text, a):
    '''print the contents of an array to the console'''
    n, m = a.shape
    print ("%s[][] = (" % text)
    for i in range(n):
        print (" (",end='')
        for j in range(m):
            print (" %g, " % a[i][j],end='')
        print ("),")
    print (")")

def test_MaxEnt_SB(report=True):
    def readTextData(filename):
        '''return q, I, dI from a 3-column text file'''
        if not os.path.exists(filename):
            raise Exception("file not found: " + filename)
        buf = [line.split() for line in open(filename, 'r').readlines()]
        buf = zip(*buf)         # transpose rows and columns
        q  = np.array(buf[0], dtype=np.float64)
        I  = np.array(buf[1], dtype=np.float64)
        dI = np.array(buf[2], dtype=np.float64)
        return q, I, dI
    print ("MaxEnt_SB: ")
    test_data_file = os.path.join( 'testinp', 'test.sas')
    rhosq = 100     # scattering contrast, 10^20 1/cm^-4
    bkg   = 0.1     #   I = I - bkg
    dMin, dMax, nRadii = 25, 9000, 40
    defaultDistLevel = 1.0e-6
    IterMax = 40
    errFac = 1.05
    
    r    = np.logspace(math.log10(dMin), math.log10(dMax), nRadii)/2
    dr   = r * (r[1]/r[0] - 1)          # step size
    f_dr = np.ndarray((nRadii)) * 0  # volume fraction histogram
    b    = np.ndarray((nRadii)) * 0 + defaultDistLevel  # MaxEnt "sky background"
    
    qVec, I, dI = readTextData(test_data_file)
    G = G_matrix(qVec,r,rhosq,SphereFF,SphereVol,args=())
    
    chisq,f_dr,Ic = MaxEnt_SB(I - bkg, dI*errFac, b, IterMax, G, report=report)
    if f_dr is None:
        print ("no solution")
        return
    
    print ("solution reached")
    for a,b,c in zip(r.tolist(), dr.tolist(), f_dr.tolist()):
        print ('%10.4f %10.4f %12.4g'%(a,b,c))

def tests():
    test_MaxEnt_SB(report=True)

if __name__ == '__main__':
    tests()

