# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2024-06-13 07:33:46 -0500 (Thu, 13 Jun 2024) $
# $Author: toby $
# $Revision: 5790 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIImpsubs.py $
# $Id: GSASIImpsubs.py 5790 2024-06-13 12:33:46Z toby $
########### SVN repository information ###################
'''
The routines here are called either directly when GSAS-II is used without multiprocessing
or in separate cores when multiprocessing is used.

These routines are designed to be used in one of two ways:

 * when multiprocessing is
   enabled (see global variable useMP) the computational routines are called in
   separate Python interpreter that is created and then deleted after use.

 * when useMP is False, these routines are called directly from the main "thread".

Note that :func:`GSASIImpsubs.InitMP` should be called before any of the other routines
in this module are used. 
'''
from __future__ import division, print_function
import multiprocessing as mp
import numpy as np
import numpy.ma as ma
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5790 $")
import GSASIIpwd as G2pwd
import GSASIIfiles as G2fil

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
#asind = lambda x: 180.*np.arcsin(x)/np.pi
#acosd = lambda x: 180.*np.arccos(x)/np.pi
#atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ncores = None

def ResetMP():
    '''Call after changing Config var 'Multiprocessing_cores' to force a resetting
    of the useMP from the parameter.
    '''
    global ncores
    ncores = None
    
def InitMP(allowMP=True):
    '''Called to initialize use of Multiprocessing
    '''
    global useMP,ncores
    if ncores is not None: return useMP,ncores
    useMP = False
    if not allowMP:
        G2fil.G2Print('Multiprocessing disabled')
        ncores = 0
        return useMP,ncores
    ncores = GSASIIpath.GetConfigValue('Multiprocessing_cores',0)
    if ncores < 0: ncores = mp.cpu_count()//2
    if ncores > 1:
        useMP = True
    if useMP:
        G2fil.G2Print('Multiprocessing with {} cores enabled'.format(ncores))
    return useMP,ncores

################################################################################
# Fobs Squared computation
################################################################################        
def InitFobsSqGlobals(x1,ratio1,shl1,xB1,xF1,im1,lamRatio1,kRatio1,xMask1,Ka21):
    '''Initialize for the computation of Fobs Squared for powder histograms.
    Puts lots of junk into the global namespace in this module.
    '''
    global x,ratio,shl,xB,xF,im,lamRatio,kRatio,xMask,Ka2,cw
    x = ma.getdata(x1)
    cw = np.diff(x)
    cw = np.append(cw,cw[-1])
    ratio = ratio1
    shl = shl1
    xB = xB1
    xF = xF1
    im = im1
    lamRatio = lamRatio1
    kRatio = kRatio1
    xMask = xMask1
    Ka2 = Ka21

def ComputeFobsSqCWbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqCW(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            resList.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList

def ComputeFobsSqTOFbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqTOF(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            resList.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList
        
def ComputeFobsSqCWBbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqCWB(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            resList.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList

def ComputeFobsSqCWAbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqCWA(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            resList.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList

def ComputeFobsSqEDbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqED(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            resList.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList

def ComputeFobsSqCW(refl,iref):
    yp = np.zeros(len(x)) # not masked
    sInt = 0
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5+im],refl[6+im],refl[7+im],shl)
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    iFin2 = iFin
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    elif iBeg < iFin:
        fp,sumfp = G2pwd.getFCJVoigt3(refl[5+im],refl[6+im],refl[7+im],shl,x[iBeg:iFin])
        yp[iBeg:iFin] = 100.*refl[11+im]*refl[9+im]*fp*cw[iBeg:iFin]/sumfp
        sInt = refl[11+im]*refl[9+im]
        if Ka2:
            pos2 = refl[5+im]+lamRatio*tand(refl[5+im]/2.0)       # + 360/pi * Dlam/lam * tan(th)
            Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6+im],refl[7+im],shl)
            iBeg2 = max(xB,np.searchsorted(x,pos2-fmin))
            iFin2 = min(np.searchsorted(x,pos2+fmax),xF)
            if iFin2 > iBeg2:
                fp2,sumfp2 = G2pwd.getFCJVoigt3(pos2,refl[6+im],refl[7+im],shl,x[iBeg2:iFin2])
                yp[iBeg2:iFin2] += 100.*refl[11+im]*refl[9+im]*kRatio*fp2*cw[iBeg2:iFin2]/sumfp2
                sInt *= 1.+kRatio
    refl8im = np.sum(np.where(ratio[iBeg:iFin2]>0.,yp[iBeg:iFin2]*ratio[iBeg:iFin2]/(refl[11+im]*(1.+kRatio)),0.0))
    return refl8im,sInt

def ComputeFobsSqCWB(refl,iref):
    yp = np.zeros(len(x)) # not masked
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.)
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    if iBeg < iFin:
        fp,sumfp = G2pwd.getEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.,x[iBeg:iFin])
        yp[iBeg:iFin] = refl[11+im]*refl[9+im]*fp*cw[iBeg:iFin]/sumfp
    refl8im = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin]/refl[11+im],0.0))
    return refl8im,refl[11+im]*refl[9+im]

def ComputeFobsSqCWA(refl,iref):
    yp = np.zeros(len(x)) # not masked
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im])
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    if iBeg < iFin:
        fp,sumfp = G2pwd.getExpFCJVoigt3(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],shl,x[iBeg:iFin])
        yp[iBeg:iFin] = refl[11+im]*refl[9+im]*fp*cw[iBeg:iFin]/sumfp
    refl8im = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin]/refl[11+im],0.0))
    return refl8im,refl[11+im]*refl[9+im]

def ComputeFobsSqTOF(refl,iref):
    yp = np.zeros(len(x)) # not masked
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im])
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    if iBeg < iFin:
        fp,sumfp = G2pwd.getEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],x[iBeg:iFin])
        yp[iBeg:iFin] = refl[11+im]*refl[9+im]*fp*cw[iBeg:iFin]/sumfp
    refl8im = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin]/refl[11+im],0.0))
    return refl8im,refl[11+im]*refl[9+im]

def ComputeFobsSqED(refl,iref):
    yp = np.zeros(len(x)) # not masked
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsED(refl[5+im],refl[6+im],refl[7+im])
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    if iBeg < iFin:
        fp,sumfp = G2pwd.getPsVoigt(refl[5+im],refl[6+im]*1.e4,refl[6+im]*100.,x[iBeg:iFin])
        yp[iBeg:iFin] = 100.*refl[9+im]*fp*cw[iBeg:iFin]/sumfp
    refl8im = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin],0.0))
    return refl8im,refl[9+im]
################################################################################
# Powder Profile computation
################################################################################        
def InitPwdrProfGlobals(im1,shl1,x1):
    '''Initialize for the computation of Fobs Squared for powder histograms.
    Puts lots of junk into the global namespace in this module.
    '''
    global im,shl,x,cw,yc
    im = im1
    shl = shl1
    x = ma.getdata(x1)
    cw = np.diff(x)
    cw = np.append(cw,cw[-1])
    # create local copies of ycalc array
    yc = np.zeros_like(x)

def ComputePwdrProfCW(profList):
    'Compute the peaks profile for a set of CW peaks and add into the yc array'
    for pos,refl,iBeg,iFin,kRatio in profList:
        fp = G2pwd.getFCJVoigt3(pos,refl[6+im],refl[7+im],shl,x[iBeg:iFin])[0]
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*kRatio*fp/cw[iBeg:iFin]
    return yc

def ComputePwdrProfTOF(profList):
    'Compute the peaks profile for a set of TOF peaks and add into the yc array'
    for pos,refl,iBeg,iFin in profList:
        fp = G2pwd.getEpsVoigt(pos,refl[12+im],refl[13+im],refl[6+im],refl[7+im],x[iBeg:iFin])[0]
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*fp*cw[iBeg:iFin]
    
def ComputePwdrProfCWB(profList):
    'Compute the peaks profile for a set of TOF peaks and add into the yc array'
    for pos,refl,iBeg,iFin in profList:
        fp = G2pwd.getEpsVoigt(pos,refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.,x[iBeg:iFin])[0]/100.
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*fp
    return yc

def ComputePwdrProfCWA(profList):
    'Compute the peaks profile for a set of TOF peaks and add into the yc array'
    for pos,refl,iBeg,iFin in profList:
        fp = G2pwd.getEpsVoigt(pos,refl[12+im],refl[13+im],refl[6+im],refl[7+im],shl,x[iBeg:iFin])[0]
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*fp
    return yc

def ComputePwdrProfED(profList):
    'Compute the peaks profile for a set of TOF peaks and add into the yc array'
    for pos,refl,iBeg,iFin in profList:
        fp = G2pwd.getPsVoigt(pos,refl[6+im]*1.e4,0.001,x[iBeg:iFin])[0]
        yc[iBeg:iFin] += refl[9+im]*fp
    return yc
