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
#### Particle form factors
###############################################################################

def SphereFF(Q,R,args=()):
    ''' Compute hard sphere form factor - can use numpy arrays
    param float:Q Q value array (usually in A-1)
    param float:R sphere radius (Usually in A - must match Q-1 units)
    returns float: form factors as array as needed
    '''
    QR = Q[:,np.newaxis]*R
    return (3./(QR**3))*(np.sin(QR)-(QR*np.cos(QR)))
    
def SpheroidFF(Q,R,args):
    ''' Compute form factor of cylindrically symmetric ellipsoid (spheroid) 
    - can use numpy arrays for R & AR; will return corresponding numpy array
    param float:Q Q value array (usually in A-1)
    param float R: radius along 2 axes of spheroid
    param float AR: aspect ratio so 3rd axis = R*AR
    returns float: form factors as array as needed
    '''
    NP = 50
    AR = args[0]
    if 0.99 < AR < 1.01:
        return SphereFF(Q,R,0)
    else:
        cth = np.linspace(0,1.,NP)
        Rct = R*np.sqrt(1.+(AR**2-1.)*cth**2)
        return np.sqrt(np.sum(SphereFF(Q[:,np.newaxis],Rct,0)**2,axis=1)/NP)
            
def CylinderFF(Q,R,args):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: R cylinder radius (A)
    param float: L cylinder length (A)
    returns float: form factor
    '''
    L = args[0]
    NP = 200
    alp = np.linspace(0,np.pi/2.,NP)
    LBessArg = 0.5*Q[:,np.newaxis]*L*np.cos(alp)
    LBess = np.where(LBessArg<1.e-6,1.,np.sin(LBessArg)/LBessArg)
    SBessArg = Q[:,np.newaxis]*R*np.sin(alp)
    SBess = np.where(SBessArg<1.e-6,0.5,scsp.jv(1,SBessArg)/SBessArg)
    return np.sqrt(2.*np.pi*np.sum(np.sin(alp)*(LBess*SBess)**2,axis=1)/NP)
    
def CylinderDFF(Q,L,args):
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: L cylinder length (A)
    param float: D cylinder diameter (A)
    returns float: form factor
    '''
    D = args[0]
    return CylinderFF(Q,D/2.,L)    
    
def CylinderARFF(Q,R,args): 
    ''' Compute form factor for cylinders - can use numpy arrays
    param float: Q Q value array (A-1)
    param float: R cylinder radius (A)
    param float: AR cylinder aspect ratio = L/D = L/2R
    returns float: form factor
    '''
    AR = args[0]
    return CylinderFF(Q,R,2.*R*AR)    
    
def UniSphereFF(Q,R,args=0):
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

def SphereVol(R,arg=()):
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
    
def UniSphereVol(R,arg=()):
    ''' Compute volume of sphere
    - numpy array friendly
    param float:R sphere radius
    returns float: volume
    '''
    return SphereVol(R)
    
def UniRodVol(R,L):
    ''' Compute cylinder volume for radius & length
    - numpy array friendly
    param float: R diameter (A)
    param float: L length (A)
    returns float:volume (A^3)
    '''
    return CylinderVol(R,L)
    
def UniRodARVol(R,AR):
    return CylinderARVol(R,AR)
    
def UniDiskVol(R,T):
    return CylinderVol(R,T)
    
def UniTubeVol(R,L,T):
    ''' Compute tube volume for radius, length & wall thickness
    - numpy array friendly
    param float: R diameter (A)
    param float: L length (A)
    param float: T tube wall thickness (A)
    returns float: volume (A^3) of tube wall
    '''
    return CylinderVol(R,L)-CylinderVol(R-T,L)
          
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
    return 1.e-4*(contrast*Vol*FF**2)     #10^-20 vs 10^-24
    
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

import os
import sys
import math
import numpy

class MaxEntException(Exception): 
    '''Any exception from this module'''
    pass

def MaxEnt_SB(datum, sigma, base, IterMax, G, image_to_data=None, data_to_image=None, report=False):
    '''
    do the complete Maximum Entropy algorithm of Skilling and Bryan
    
    :param float datum[]:
    :param float sigma[]:
    :param float base[]:
    :param int IterMax:
    :param float[][] G: transformation matrix
    :param obj image_to_data: opus function (defaults to opus)
    :param obj data_to_image: tropus function (defaults to tropus)
    
    :returns float[]: :math:`f(r) dr`
    '''
        
    TEST_LIMIT        = 0.10                    # for convergence
    CHI_SQR_LIMIT     = 0.01                    # maximum difference in ChiSqr for a solution
    SEARCH_DIRECTIONS = 3                       # <10.  This code requires value = 3
    RESET_STRAYS      = 1                       # was 0.001, correction of stray negative values
    DISTANCE_LIMIT_FACTOR = 0.1                 # limitation on df to constrain runaways
    
    MAX_MOVE_LOOPS    = 500                     # for no solution in routine: move, 
    MOVE_PASSES       = 0.001                   # convergence test in routine: move

    def opus (data, G):
        '''
        opus: transform data-space -> solution-space:  [G] * data
        
        default definition, caller can use this definition or provide an alternative
        
        :param float[M] data: observations, ndarray of shape (M)
        :param float[M][N] G: transformation matrix, ndarray of shape (M,N)
        :returns float[N]: calculated image, ndarray of shape (N)
        '''
        return G.dot(data)

    def tropus (image, G):
        '''
        tropus: transform solution-space -> data-space:  [G]^tr * image
        
        default definition, caller can use this definition or provide an alternative
        
        :param float[N] image: solution, ndarray of shape (N)
        :param float[M][N] G: transformation matrix, ndarray of shape (M,N)
        :returns float[M]: calculated data, ndarray of shape (M)
        '''
        return G.transpose().dot(image)

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
        fl = numpy.ndarray((n, n))*0
        bl = numpy.ndarray((n))*0
        
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
                    y = math.sqrt(z)
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
        beta = numpy.ndarray((n))
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
    xi      = 0*numpy.ndarray((SEARCH_DIRECTIONS, n))
    eta     = 0*numpy.ndarray((SEARCH_DIRECTIONS, npt))
    beta    = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    # s1      = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    # c1      = 0*numpy.ndarray((SEARCH_DIRECTIONS))
    s2      = 0*numpy.ndarray((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))
    c2      = 0*numpy.ndarray((SEARCH_DIRECTIONS, SEARCH_DIRECTIONS))

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
        sgrad = -numpy.log(f/base) / (blank*math.exp (1.0))  # sgrad[i] = del(S)/del(f[i])
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
        S = -sum((f/fSum) * numpy.log(f/fSum))      # normalized entropy, S&B eq. 1
        z = (datum - image_to_data (f, G)) / sigma  # standardized residuals
        chisq = sum(z*z)                            # report this ChiSq

        if report:
            print "%3d/%3d" % ((iter+1), IterMax)
            print " %5.2lf%% %8lg" % (100*test, S)
            if iter > 0:
                value = 100*( math.sqrt(chisq/chtarg)-1)
            else:
                value = 0
            print " %12.5lg %10.4lf" % ( math.sqrt (chtarg/npt), value )
            print "%12.6lg %8.2lf\n" % (fSum, 100*fChange/fSum)

        # See if we have finished our task.
        # do the hardest test first
        if (abs(chisq/chizer-1.0) < CHI_SQR_LIMIT) and  (test < TEST_LIMIT):
            return f,image_to_data (f, G)     # solution FOUND returns here
    
    return f,image_to_data (f, G)       # no solution after IterMax iterations

    
################################################################################
#### MaxEnt testing stuff
################################################################################

def print_vec(text, a):
    '''print the contents of a vector to the console'''
    n = a.shape[0]
    print "%s[ = (" % text,
    for i in range(n):
        s = " %g, " % a[i]
        print s,
    print ")"

def print_arr(text, a):
    '''print the contents of an array to the console'''
    n, m = a.shape
    print "%s[][] = (" % text
    for i in range(n):
        print " (",
        for j in range(m):
            print " %g, " % a[i][j],
        print "),"
    print ")"

def test_MaxEnt_SB(report=True):
    def readTextData(filename):
        '''return q, I, dI from a 3-column text file'''
        if not os.path.exists(filename):
            raise Exception("file not found: " + filename)
        buf = [line.split() for line in open(filename, 'r').readlines()]
        M = len(buf)
        buf = zip(*buf)         # transpose rows and columns
        q  = numpy.array(buf[0], dtype=numpy.float64)
        I  = numpy.array(buf[1], dtype=numpy.float64)
        dI = numpy.array(buf[2], dtype=numpy.float64)
        return q, I, dI
    print "MaxEnt_SB: "
    test_data_file = os.path.join( 'testinp', 'test.sas')
    rhosq = 100     # scattering contrast, 10^20 1/cm^-4
    bkg   = 0.1     #   I = I - bkg
    dMin, dMax, nRadii = 25, 9000, 40
    defaultDistLevel = 1.0e-6
    IterMax = 40
    errFac = 1.05
    
    r    = numpy.logspace(math.log10(dMin), math.log10(dMax), nRadii)/2
    dr   = r * (r[1]/r[0] - 1)          # step size
    f_dr = numpy.ndarray((nRadii)) * 0  # volume fraction histogram
    b    = numpy.ndarray((nRadii)) * 0 + defaultDistLevel  # MaxEnt "sky background"
    
    qVec, I, dI = readTextData(test_data_file)
    G = G_matrix(qVec,r,rhosq,SphereFF,SphereVol,args=())
    
    f_dr,Ic = MaxEnt_SB(I - bkg, dI*errFac, b, IterMax, G, report=report)
    if f_dr is None:
        print "no solution"
        return
    
    print "solution reached"
    for a,b,c in zip(r.tolist(), dr.tolist(), f_dr.tolist()):
        print '%10.4f %10.4f %12.4g'%(a,b,c)

def tests():
    test_MaxEnt_SB(report=True)

if __name__ == '__main__':
    tests()
    
###############################################################################
#### Size distribution
###############################################################################

def SizeDistribution(Profile,ProfDict,Limits,Substances,Sample,data):
    shapes = {'Spheroid':[SpheroidFF,SpheroidVol],'Cylinder':[CylinderDFF,CylinderDVol],
        'Cylinder AR':[CylinderARFF,CylinderARVol],'Unified sphere':[UniSphereFF,UniSphereVol],
        'Unified rod':[UniRodFF,UniRodVol],'Unified rod AR':[UniRodARFF,UniRodARVol],
        'Unified disk':[UniDiskFF,UniDiskVol]}
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
    BinsBack = np.ones_like(Bins)*1.e-6
    Back = data['Back']
    Q,Io,wt,Ic,Ib = Profile[:5]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    Contrast = Sample['Contrast'][1]
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)
    if Back[1]:
        Ib = Back[0]
        Ic[Ibeg:Ifin] = Back[0]
    Gmat = G_matrix(Q[Ibeg:Ifin],Bins,Contrast,shapes[Shape][0],shapes[Shape][1],args=Parms)
    BinMag,Ic[Ibeg:Ifin] = MaxEnt_SB(Io[Ibeg:Ifin]-Back[0],1./np.sqrt(wt[Ibeg:Ifin]),BinsBack,
        data['Size']['MaxEnt']['Niter'],Gmat)
    print BinMag.shape
    data['Size']['Distribution'] = [Bins,Dbins,BinMag]
    print np.sum(BinMag)
        
    
