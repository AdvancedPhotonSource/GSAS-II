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
#### Particle form factors & volumes
###############################################################################

def SpheroidFF(Q,R,AR):
    '''
    '''
    QR = Q*R
    if 0.98 < AR < 1.02:
        return (3./(QR**3))*(np.sin(QR)-(QR*np.cos(QR)))
        
    
def SpheroidVol(R,AR):
    ''' Compute volume of cylindrically symmetric ellipsoid (spheroid) 
    - can use numpy arrays for R & AR; will return corresponding numpy array
    param float R: radius along 2 axes of spheroid
    param float AR: aspect ratio so 3rd axis = R*AR
    returns float: volume
    '''
    return AR*(4./3.)*np.pi*R**3
    
def SizeDistribution(Profile,Limits,Substances,Sample,data):
    print data['Size']
#    print Limits
#    print Substances
#    print Sample
#    print Profile
    
