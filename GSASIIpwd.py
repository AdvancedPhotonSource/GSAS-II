#/usr/bin/env python
# -*- coding: utf-8 -*-
'''
*GSASII powder calculation module*
==================================

'''
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
from __future__ import division, print_function
import sys
import math
import time
import os
import os.path
import subprocess as subp
import datetime as dt
import copy

import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import random as rand
import numpy.fft as fft
import scipy.interpolate as si
import scipy.stats as st
import scipy.optimize as so
import scipy.special as sp

import GSASIIpath
filversion = "$Revision$"
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIImath as G2mth
try:
    import pypowder as pyd
except ImportError:
    print ('pypowder is not available - profile calcs. not allowed')
try:
    import pydiffax as pyx
except ImportError:
    print ('pydiffax is not available for this platform')
import GSASIIfiles as G2fil

    
# trig functions in degrees
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
npT2stl = lambda tth, wave: 2.0*npsind(tth/2.0)/wave    #=d*
npT2q = lambda tth,wave: 2.0*np.pi*npT2stl(tth,wave)    #=2pi*d*
ateln2 = 8.0*math.log(2.0)
sateln2 = np.sqrt(ateln2)
nxs = np.newaxis

################################################################################
#### Powder utilities
################################################################################

def PhaseWtSum(G2frame,histo):
    '''
    Calculate sum of phase mass*phase fraction for PWDR data (exclude magnetic phases)
    
    :param G2frame: GSASII main frame structure
    :param str histo: histogram name
    :returns: sum(scale*mass) for phases in histo
    '''
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    wtSum = 0.0
    for phase in Phases:
        if Phases[phase]['General']['Type'] != 'magnetic':
            if histo in Phases[phase]['Histograms']:
                if not Phases[phase]['Histograms'][histo]['Use']: continue
                mass = Phases[phase]['General']['Mass']
                phFr = Phases[phase]['Histograms'][histo]['Scale'][0]
                wtSum += mass*phFr
    return wtSum
    
################################################################################
#### GSASII pwdr & pdf calculation routines
################################################################################
        
def Transmission(Geometry,Abs,Diam):
    '''
    Calculate sample transmission

    :param str Geometry: one of 'Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate'
    :param float Abs: absorption coeff in cm-1
    :param float Diam: sample thickness/diameter in mm
    '''
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
        
def SurfaceRough(SRA,SRB,Tth):
    ''' Suortti (J. Appl. Cryst, 5,325-331, 1972) surface roughness correction
    :param float SRA: Suortti surface roughness parameter
    :param float SRB: Suortti surface roughness parameter
    :param float Tth: 2-theta(deg) - can be numpy array
    
    '''
    sth = npsind(Tth/2.)
    T1 = np.exp(-SRB/sth)
    T2 = SRA+(1.-SRA)*np.exp(-SRB)
    return (SRA+(1.-SRA)*T1)/T2
    
def SurfaceRoughDerv(SRA,SRB,Tth):
    ''' Suortti surface roughness correction derivatives
    :param float SRA: Suortti surface roughness parameter (dimensionless)
    :param float SRB: Suortti surface roughness parameter (dimensionless)
    :param float Tth: 2-theta(deg) - can be numpy array
    :return list: [dydSRA,dydSRB] derivatives to be used for intensity derivative
    '''
    sth = npsind(Tth/2.)
    T1 = np.exp(-SRB/sth)
    T2 = SRA+(1.-SRA)*np.exp(-SRB)
    Trans = (SRA+(1.-SRA)*T1)/T2
    dydSRA = ((1.-T1)*T2-(1.-np.exp(-SRB))*Trans)/T2**2
    dydSRB = ((SRA-1.)*T1*T2/sth-Trans*(SRA-T2))/T2**2
    return [dydSRA,dydSRB]

def Absorb(Geometry,MuR,Tth,Phi=0,Psi=0):
    '''Calculate sample absorption
    :param str Geometry: one of 'Cylinder','Bragg-Brentano','Tilting Flat Plate in transmission','Fixed flat plate'
    :param float MuR: absorption coeff * sample thickness/2 or radius
    :param Tth: 2-theta scattering angle - can be numpy array
    :param float Phi: flat plate tilt angle - future
    :param float Psi: flat plate tilt axis - future
    '''
    
    def muRunder3(MuR,Sth2):
        T0 = 16.0/(3.*np.pi)
        T1 = (25.99978-0.01911*Sth2**0.25)*np.exp(-0.024551*Sth2)+ \
            0.109561*np.sqrt(Sth2)-26.04556
        T2 = -0.02489-0.39499*Sth2+1.219077*Sth2**1.5- \
            1.31268*Sth2**2+0.871081*Sth2**2.5-0.2327*Sth2**3
        T3 = 0.003045+0.018167*Sth2-0.03305*Sth2**2
        Trns = -T0*MuR-T1*MuR**2-T2*MuR**3-T3*MuR**4
        return np.exp(Trns)
    
    def muRover3(MuR,Sth2):
        T1 = 1.433902+11.07504*Sth2-8.77629*Sth2*Sth2+ \
            10.02088*Sth2**3-3.36778*Sth2**4
        T2 = (0.013869-0.01249*Sth2)*np.exp(3.27094*Sth2)+ \
            (0.337894+13.77317*Sth2)/(1.0+11.53544*Sth2)**1.555039
        T3 = 1.933433/(1.0+23.12967*Sth2)**1.686715- \
            0.13576*np.sqrt(Sth2)+1.163198
        T4 = 0.044365-0.04259/(1.0+0.41051*Sth2)**148.4202
        Trns = (T1-T4)/(1.0+T2*(MuR-3.0))**T3+T4
        return Trns/100.
        
    Sth2 = npsind(Tth/2.0)**2
    if 'Cylinder' in Geometry:      #Lobanov & Alte da Veiga for 2-theta = 0; beam fully illuminates sample
        if 'array' in str(type(MuR)):
            MuRSTh2 = np.vstack((MuR,Sth2))
            AbsCr = np.where(MuRSTh2[0]<=3.0,muRunder3(MuRSTh2[0],MuRSTh2[1]),muRover3(MuRSTh2[0],MuRSTh2[1]))
            return AbsCr
        else:
            if MuR <= 3.0:
                return muRunder3(MuR,Sth2)
            else:
                return muRover3(MuR,Sth2)
    elif 'Bragg' in Geometry:
        return 1.0
    elif 'Fixed' in Geometry: #assumes sample plane is perpendicular to incident beam
        # and only defined for 2theta < 90
        MuT = 2.*MuR
        T1 = np.exp(-MuT)
        T2 = np.exp(-MuT/npcosd(Tth))
        Tb = MuT-MuT/npcosd(Tth)
        return (T2-T1)/Tb
    elif 'Tilting' in Geometry: #assumes symmetric tilt so sample plane is parallel to diffraction vector
        MuT = 2.*MuR
        cth = npcosd(Tth/2.0)
        return np.exp(-MuT/cth)/cth
        
def AbsorbDerv(Geometry,MuR,Tth,Phi=0,Psi=0):
    'needs a doc string'
    dA = 0.001
    AbsP = Absorb(Geometry,MuR+dA,Tth,Phi,Psi)
    if MuR:
        AbsM = Absorb(Geometry,MuR-dA,Tth,Phi,Psi)
        return (AbsP-AbsM)/(2.0*dA)
    else:
        return (AbsP-1.)/dA
        
def Polarization(Pola,Tth,Azm=0.0):
    """   Calculate angle dependent x-ray polarization correction (not scaled correctly!)

    :param Pola: polarization coefficient e.g 1.0 fully polarized, 0.5 unpolarized
    :param Azm: azimuthal angle e.g. 0.0 in plane of polarization
    :param Tth: 2-theta scattering angle - can be numpy array
      which (if either) of these is "right"?
    :return: (pola, dpdPola)
      * pola = ((1-Pola)*npcosd(Azm)**2+Pola*npsind(Azm)**2)*npcosd(Tth)**2+ \
        (1-Pola)*npsind(Azm)**2+Pola*npcosd(Azm)**2
      * dpdPola: derivative needed for least squares

    """
    cazm = npcosd(Azm)**2
    sazm = npsind(Azm)**2
    pola = ((1.0-Pola)*cazm+Pola*sazm)*npcosd(Tth)**2+(1.0-Pola)*sazm+Pola*cazm
    dpdPola = -npsind(Tth)**2*(sazm-cazm)
    return pola,dpdPola
    
def Oblique(ObCoeff,Tth):
    'currently assumes detector is normal to beam'
    if ObCoeff:
        return (1.-ObCoeff)/(1.0-np.exp(np.log(ObCoeff)/npcosd(Tth)))
    else:
        return 1.0
                
def Ruland(RulCoff,wave,Q,Compton):
    'needs a doc string'
    C = 2.9978e8
    D = 1.5e-3
    hmc = 0.024262734687
    sinth2 = (Q*wave/(4.0*np.pi))**2
    dlam = (wave**2)*Compton*Q/C
    dlam_c = 2.0*hmc*sinth2-D*wave**2
    return 1.0/((1.0+dlam/RulCoff)*(1.0+(np.pi*dlam_c/(dlam+RulCoff))**2))
    
def LorchWeight(Q):
    'needs a doc string'
    return np.sin(np.pi*(Q[-1]-Q)/(2.0*Q[-1]))
            
def GetAsfMean(ElList,Sthl2):
    '''Calculate various scattering factor terms for PDF calcs

    :param dict ElList: element dictionary contains scattering factor coefficients, etc.
    :param np.array Sthl2: numpy array of sin theta/lambda squared values
    :returns: mean(f^2), mean(f)^2, mean(compton)
    '''
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
    'needs a doc string'
    sumNoAtoms = 0.0
    for El in ElList:
        sumNoAtoms += ElList[El]['FormulaNo']
    return sumNoAtoms/Vol
           
def CalcPDF(data,inst,limits,xydata):
    '''Computes I(Q), S(Q) & G(r) from Sample, Bkg, etc. diffraction patterns loaded into
    dict xydata; results are placed in xydata.
    Calculation parameters are found in dicts data and inst and list limits.
    The return value is at present an empty list.
    '''
    auxPlot = []
    if 'T' in inst['Type'][0]:
        Ibeg = 0
        Ifin = len(xydata['Sample'][1][0])
    else:
        Ibeg = np.searchsorted(xydata['Sample'][1][0],limits[0])
        Ifin = np.searchsorted(xydata['Sample'][1][0],limits[1])+1
    #subtract backgrounds - if any & use PWDR limits
    IofQ = copy.deepcopy(xydata['Sample'])
    IofQ[1] = np.array([I[Ibeg:Ifin] for I in IofQ[1]])
    if data['Sample Bkg.']['Name']:
        IofQ[1][1] += xydata['Sample Bkg.'][1][1][Ibeg:Ifin]*data['Sample Bkg.']['Mult']
    if data['Container']['Name']:
        xycontainer = xydata['Container'][1][1]*data['Container']['Mult']
        if data['Container Bkg.']['Name']:
            xycontainer += xydata['Container Bkg.'][1][1][Ibeg:Ifin]*data['Container Bkg.']['Mult']
        IofQ[1][1] += xycontainer[Ibeg:Ifin]
    data['IofQmin'] = IofQ[1][1][-1]
    IofQ[1][1] -= data.get('Flat Bkg',0.)
    #get element data & absorption coeff.
    ElList = data['ElList']
    Tth = IofQ[1][0]    #2-theta or TOF!
    if 'X' in inst['Type'][0]:
        Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
        #Apply angle dependent corrections
        MuR = Abs*data['Diam']/20.0
        IofQ[1][1] /= Absorb(data['Geometry'],MuR,Tth)
        IofQ[1][1] /= Polarization(inst['Polariz.'][1],Tth,Azm=inst['Azimuth'][1])[0]
        if data['DetType'] == 'Image plate':
            IofQ[1][1] *= Oblique(data['ObliqCoeff'],Tth)
    elif 'T' in inst['Type'][0]:    #neutron TOF normalized data - needs wavelength dependent absorption
        wave = 2.*G2lat.TOF2dsp(inst,IofQ[1][0])*npsind(inst['2-theta'][1]/2.)
        Els = ElList.keys()
        Isotope = {El:'Nat. abund.' for El in Els}
        GD = {'AtomTypes':ElList,'Isotope':Isotope}
        BLtables = G2elem.GetBLtable(GD)
        FP,FPP = G2elem.BlenResTOF(Els,BLtables,wave)
        Abs = np.zeros(len(wave))
        for iel,El in enumerate(Els):
            BL = BLtables[El][1]
            SA = BL['SA']*wave/1.798197+4.0*np.pi*FPP[iel]**2 #+BL['SL'][1]?
            SA *= ElList[El]['FormulaNo']/data['Form Vol']
            Abs += SA
        MuR = Abs*data['Diam']/2.
        IofQ[1][1] /= Absorb(data['Geometry'],MuR,inst['2-theta'][1]*np.ones(len(wave)))        
    XY = IofQ[1]    
    #convert to Q
#    nQpoints = len(XY[0])     #points for Q interpolation
    nQpoints = 5000
    if 'C' in inst['Type'][0]:
        wave = G2mth.getWave(inst)
        minQ = npT2q(Tth[0],wave)
        maxQ = npT2q(Tth[-1],wave)    
        Qpoints = np.linspace(0.,maxQ,nQpoints,endpoint=True)
        dq = Qpoints[1]-Qpoints[0]
        XY[0] = npT2q(XY[0],wave)
        Qdata = si.griddata(XY[0],XY[1],Qpoints,method='linear',fill_value=XY[1][0])    #interpolate I(Q)
    elif 'T' in inst['Type'][0]:
        difC = inst['difC'][1]
        minQ = 2.*np.pi*difC/Tth[-1]
        maxQ = 2.*np.pi*difC/Tth[0]
        Qpoints = np.linspace(0.,maxQ,nQpoints,endpoint=True)
        dq = Qpoints[1]-Qpoints[0]
        XY[0] = 2.*np.pi*difC/XY[0]
        Qdata = si.griddata(XY[0],XY[1],Qpoints,method='linear',fill_value=XY[1][-1])    #interpolate I(Q)
    Qdata -= np.min(Qdata)*data['BackRatio']
    
    qLimits = data['QScaleLim']
    maxQ = np.searchsorted(Qpoints,min(Qpoints[-1],qLimits[1]))+1
    minQ = np.searchsorted(Qpoints,min(qLimits[0],0.90*Qpoints[-1]))
    qLimits = [Qpoints[minQ],Qpoints[maxQ-1]]
    newdata = []
    if len(IofQ) < 3:
        xydata['IofQ'] = [IofQ[0],[Qpoints,Qdata],'']
    else:
        xydata['IofQ'] = [IofQ[0],[Qpoints,Qdata],IofQ[2]]
    for item in xydata['IofQ'][1]:
        newdata.append(item[:maxQ])
    xydata['IofQ'][1] = newdata
    
    xydata['SofQ'] = copy.deepcopy(xydata['IofQ'])
    if 'XC' in inst['Type'][0]:
        FFSq,SqFF,CF = GetAsfMean(ElList,(xydata['SofQ'][1][0]/(4.0*np.pi))**2)  #these are <f^2>,<f>^2,Cf
    else: #TOF
        CF = np.zeros(len(xydata['SofQ'][1][0]))
        FFSq = np.ones(len(xydata['SofQ'][1][0]))
        SqFF = np.ones(len(xydata['SofQ'][1][0]))
    Q = xydata['SofQ'][1][0]
#    auxPlot.append([Q,np.copy(CF),'CF-unCorr'])
    if 'XC' in inst['Type'][0]:
        ruland = Ruland(data['Ruland'],wave,Q,CF)
#        auxPlot.append([Q,ruland,'Ruland'])      
        CF *= ruland
#    auxPlot.append([Q,CF,'CF-Corr'])
    scale = np.sum((FFSq+CF)[minQ:maxQ])/np.sum(xydata['SofQ'][1][1][minQ:maxQ])
    xydata['SofQ'][1][1] *= scale
    if 'XC' in inst['Type'][0]:
        xydata['SofQ'][1][1] -= CF
    xydata['SofQ'][1][1] = xydata['SofQ'][1][1]/SqFF
    scale = len(xydata['SofQ'][1][1][minQ:maxQ])/np.sum(xydata['SofQ'][1][1][minQ:maxQ])
    xydata['SofQ'][1][1] *= scale
    xydata['FofQ'] = copy.deepcopy(xydata['SofQ'])
    xydata['FofQ'][1][1] = xydata['FofQ'][1][0]*(xydata['SofQ'][1][1]-1.0)
    if data['Lorch']:
        xydata['FofQ'][1][1] *= LorchWeight(Q)    
    xydata['GofR'] = copy.deepcopy(xydata['FofQ'])
    xydata['gofr'] = copy.deepcopy(xydata['FofQ'])
    nR = len(xydata['GofR'][1][1])
    Rmax = GSASIIpath.GetConfigValue('PDF_Rmax',100.)
    mul = int(round(2.*np.pi*nR/(Rmax*qLimits[1])))
#    mul = int(round(2.*np.pi*nR/(data.get('Rmax',100.)*qLimits[1])))
    R = 2.*np.pi*np.linspace(0,nR,nR,endpoint=True)/(mul*qLimits[1])
    xydata['GofR'][1][0] = R
    xydata['gofr'][1][0] = R
    GR = -dq*np.imag(fft.fft(xydata['FofQ'][1][1],mul*nR)[:nR])
    xydata['GofR'][1][1] = GR
    gr = GR/(np.pi*R)
    xydata['gofr'][1][1] = gr
    numbDen = 0.
    if 'ElList' in data:
        numbDen = GetNumDensity(data['ElList'],data['Form Vol'])
    if data.get('noRing',True):
        Rmin = data['Rmin']
        xydata['gofr'][1][1] = np.where(R<Rmin,-4.*numbDen,xydata['gofr'][1][1])
        xydata['GofR'][1][1] = np.where(R<Rmin,-4.*R*np.pi*numbDen,xydata['GofR'][1][1])
    return auxPlot
    
def PDFPeakFit(peaks,data):
    rs2pi = 1./np.sqrt(2*np.pi)
    
    def MakeParms(peaks):
        varyList = []
        parmDict = {'slope':peaks['Background'][1][1]}
        if peaks['Background'][2]:
            varyList.append('slope')
        for i,peak in enumerate(peaks['Peaks']):
            parmDict['PDFpos;'+str(i)] = peak[0]
            parmDict['PDFmag;'+str(i)] = peak[1]
            parmDict['PDFsig;'+str(i)] = peak[2]
            if 'P' in peak[3]:
                varyList.append('PDFpos;'+str(i))
            if 'M' in peak[3]:
                varyList.append('PDFmag;'+str(i))
            if 'S' in peak[3]:
                varyList.append('PDFsig;'+str(i))
        return parmDict,varyList
        
    def SetParms(peaks,parmDict,varyList):
        if 'slope' in varyList:
            peaks['Background'][1][1] = parmDict['slope']
        for i,peak in enumerate(peaks['Peaks']):
            if 'PDFpos;'+str(i) in varyList:
                peak[0] = parmDict['PDFpos;'+str(i)]
            if 'PDFmag;'+str(i) in varyList:
                peak[1] = parmDict['PDFmag;'+str(i)]
            if 'PDFsig;'+str(i) in varyList:
                peak[2] = parmDict['PDFsig;'+str(i)]
        
    
    def CalcPDFpeaks(parmdict,Xdata):
        Z = parmDict['slope']*Xdata
        ipeak = 0
        while True:
            try:
                pos = parmdict['PDFpos;'+str(ipeak)]
                mag = parmdict['PDFmag;'+str(ipeak)]
                wid = parmdict['PDFsig;'+str(ipeak)]
                wid2 = 2.*wid**2
                Z += mag*rs2pi*np.exp(-(Xdata-pos)**2/wid2)/wid
                ipeak += 1
            except KeyError:        #no more peaks to process
                return Z
                
    def errPDFProfile(values,xdata,ydata,parmdict,varylist):        
        parmdict.update(zip(varylist,values))
        M = CalcPDFpeaks(parmdict,xdata)-ydata
        return M
            
    newpeaks = copy.copy(peaks)
    iBeg = np.searchsorted(data[1][0],newpeaks['Limits'][0])
    iFin = np.searchsorted(data[1][0],newpeaks['Limits'][1])+1
    X = data[1][0][iBeg:iFin]
    Y = data[1][1][iBeg:iFin]
    parmDict,varyList = MakeParms(peaks)
    if not len(varyList):
        G2fil.G2Print (' Nothing varied')
        return newpeaks,None,None,None,None,None
    
    Rvals = {}
    values =  np.array(Dict2Values(parmDict, varyList))
    result = so.leastsq(errPDFProfile,values,full_output=True,ftol=0.0001,
           args=(X,Y,parmDict,varyList))
    chisq = np.sum(result[2]['fvec']**2)
    Values2Dict(parmDict, varyList, result[0])
    SetParms(peaks,parmDict,varyList)
    Rvals['Rwp'] = np.sqrt(chisq/np.sum(Y**2))*100.      #to %
    chisq = np.sum(result[2]['fvec']**2)/(len(X)-len(values))   #reduced chi^2 = M/(Nobs-Nvar)
    sigList = list(np.sqrt(chisq*np.diag(result[1])))    
    Z = CalcPDFpeaks(parmDict,X)
    newpeaks['calc'] = [X,Z]
    return newpeaks,result[0],varyList,sigList,parmDict,Rvals    
    
def MakeRDF(RDFcontrols,background,inst,pwddata):
    import scipy.signal as signal
    auxPlot = []
    if 'C' in inst['Type'][0]:
        Tth = pwddata[0]
        wave = G2mth.getWave(inst)
        minQ = npT2q(Tth[0],wave)
        maxQ = npT2q(Tth[-1],wave)
        powQ = npT2q(Tth,wave)  
    elif 'T' in inst['Type'][0]:
        TOF = pwddata[0]
        difC = inst['difC'][1]
        minQ = 2.*np.pi*difC/TOF[-1]
        maxQ = 2.*np.pi*difC/TOF[0]
        powQ = 2.*np.pi*difC/TOF
    piDQ = np.pi/(maxQ-minQ)
    Qpoints = np.linspace(minQ,maxQ,len(pwddata[0]),endpoint=True)
    if RDFcontrols['UseObsCalc'] == 'obs-calc':
        Qdata = si.griddata(powQ,pwddata[1]-pwddata[3],Qpoints,method=RDFcontrols['Smooth'],fill_value=0.)
    elif RDFcontrols['UseObsCalc'] == 'obs-back':
        Qdata = si.griddata(powQ,pwddata[1]-pwddata[4],Qpoints,method=RDFcontrols['Smooth'],fill_value=pwddata[1][0])
    elif RDFcontrols['UseObsCalc'] == 'calc-back':
        Qdata = si.griddata(powQ,pwddata[3]-pwddata[4],Qpoints,method=RDFcontrols['Smooth'],fill_value=pwddata[1][0])
    Qdata *= np.sin((Qpoints-minQ)*piDQ)/piDQ
    Qdata *= 0.5*np.sqrt(Qpoints)       #Qbin normalization
#    GSASIIpath.IPyBreak()
    dq = Qpoints[1]-Qpoints[0]
    nR = len(Qdata)
    R = 0.5*np.pi*np.linspace(0,nR,nR)/(4.*maxQ)
    iFin = np.searchsorted(R,RDFcontrols['maxR'])+1
    bBut,aBut = signal.butter(4,0.01)
    Qsmooth = signal.filtfilt(bBut,aBut,Qdata)
#    auxPlot.append([Qpoints,Qdata,'interpolate:'+RDFcontrols['Smooth']])
#    auxPlot.append([Qpoints,Qsmooth,'interpolate:'+RDFcontrols['Smooth']])
    DofR = dq*np.imag(fft.fft(Qsmooth,16*nR)[:nR])
#    DofR = dq*np.imag(ft.fft(Qsmooth,16*nR)[:nR])
    auxPlot.append([R[:iFin],DofR[:iFin],'D(R) for '+RDFcontrols['UseObsCalc']])    
    return auxPlot

# PDF optimization =============================================================
def OptimizePDF(data,xydata,limits,inst,showFit=True,maxCycles=5):
    import scipy.optimize as opt
    numbDen = GetNumDensity(data['ElList'],data['Form Vol'])
    Min,Init,Done = SetupPDFEval(data,xydata,limits,inst,numbDen)
    xstart = Init()
    bakMul = data['Sample Bkg.']['Mult']
    if showFit:
        rms = Min(xstart)
        G2fil.G2Print('  Optimizing corrections to improve G(r) at low r')
        if data['Sample Bkg.'].get('Refine',False):
#            data['Flat Bkg'] = 0.
            G2fil.G2Print('  start: Ruland={:.3f}, Sample Bkg mult={:.3f} (RMS:{:.4f})'.format(
                data['Ruland'],data['Sample Bkg.']['Mult'],rms))
        else:
            G2fil.G2Print('  start: Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f} (RMS:{:.4f})'.format(
                data['Flat Bkg'],data['BackRatio'],data['Ruland'],rms))
    if data['Sample Bkg.'].get('Refine',False):
        res = opt.minimize(Min,xstart,bounds=([0.01,1],[1.2*bakMul,0.8*bakMul]),
                    method='L-BFGS-B',options={'maxiter':maxCycles},tol=0.001)
    else:
        res = opt.minimize(Min,xstart,bounds=([0,None],[0,1],[0.01,1]),
                    method='L-BFGS-B',options={'maxiter':maxCycles},tol=0.001)
    Done(res['x'])
    if showFit:
        if res['success']:
            msg = 'Converged'
        else:
            msg = 'Not Converged'
        if data['Sample Bkg.'].get('Refine',False):
            G2fil.G2Print('  end:   Ruland={:.3f}, Sample Bkg mult={:.3f} (RMS:{:.4f}) *** {} ***\n'.format(
                data['Ruland'],data['Sample Bkg.']['Mult'],res['fun'],msg))
        else:
            G2fil.G2Print('  end:   Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f}) *** {} ***\n'.format(
                data['Flat Bkg'],data['BackRatio'],data['Ruland'],res['fun'],msg))
    return res

def SetupPDFEval(data,xydata,limits,inst,numbDen):
    Data = copy.deepcopy(data)
    BkgMax = 1.
    def EvalLowPDF(arg):
        '''Objective routine -- evaluates the RMS deviations in G(r)
        from -4(pi)*#density*r for for r<Rmin
        arguments are ['Flat Bkg','BackRatio','Ruland'] scaled so that
        the min & max values are between 0 and 1. 
        '''
        if Data['Sample Bkg.'].get('Refine',False):
            R,S = arg
            Data['Sample Bkg.']['Mult'] = S
        else:
            F,B,R = arg
            Data['Flat Bkg'] = F*BkgMax
            Data['BackRatio'] = B
        Data['Ruland'] = R/10.
        CalcPDF(Data,inst,limits,xydata)
        # test low r computation
        g = xydata['GofR'][1][1]
        r = xydata['GofR'][1][0]
        g0 = g[r < Data['Rmin']] + 4*np.pi*r[r < Data['Rmin']]*numbDen
        M = sum(g0**2)/len(g0)
        return M
    def GetCurrentVals():
        '''Get the current ['Flat Bkg','BackRatio','Ruland'] with scaling
        '''
        if data['Sample Bkg.'].get('Refine',False):
                return [max(10*data['Ruland'],.05),data['Sample']['Mult']]
        try:
            F = data['Flat Bkg']/BkgMax
        except:
            F = 0
        return [F,data['BackRatio'],max(10*data['Ruland'],.05)]
    def SetFinalVals(arg):
        '''Set the 'Flat Bkg', 'BackRatio' & 'Ruland' values from the
        scaled, refined values and plot corrected region of G(r) 
        '''
        if data['Sample Bkg.'].get('Refine',False):
            R,S = arg
            data['Sample Bkg.']['Mult'] = S
        else:
            F,B,R = arg
            data['Flat Bkg'] = F*BkgMax
            data['BackRatio'] = B
        data['Ruland'] = R/10.
        CalcPDF(data,inst,limits,xydata)
    EvalLowPDF(GetCurrentVals())
    BkgMax = max(xydata['IofQ'][1][1])/50.
    return EvalLowPDF,GetCurrentVals,SetFinalVals

################################################################################        
#### GSASII peak fitting routines: Finger, Cox & Jephcoat model        
################################################################################

def factorize(num):
    ''' Provide prime number factors for integer num
    :returns: dictionary of prime factors (keys) & power for each (data)
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

    :param int nmin: minimum data size >= 1
    :param int nmax: maximum data size > nmin
    :param int thresh: maximum prime factor allowed
    :Returns: list of data sizes where the maximum prime factor is < thresh
    ''' 
    plist = []
    nmin = max(1,nmin)
    nmax = max(nmin+1,nmax)
    for p in range(nmin,nmax):
        if max(list(factorize(p).keys())) < thresh:
            plist.append(p)
    return plist

np.seterr(divide='ignore')

# Normal distribution

# loc = mu, scale = std
_norm_pdf_C = 1./math.sqrt(2*math.pi)
class norm_gen(st.rv_continuous):
    'needs a doc string'
      
    def pdf(self,x,*args,**kwds):
        loc,scale=kwds['loc'],kwds['scale']
        x = (x-loc)/scale
        return np.exp(-x**2/2.0) * _norm_pdf_C / scale
        
norm = norm_gen(name='norm',longname='A normal',extradoc="""

Normal distribution

The location (loc) keyword specifies the mean.
The scale (scale) keyword specifies the standard deviation.

normal.pdf(x) = exp(-x**2/2)/sqrt(2*pi)
""")

## Cauchy

# median = loc

class cauchy_gen(st.rv_continuous):
    'needs a doc string'

    def pdf(self,x,*args,**kwds):
        loc,scale=kwds['loc'],kwds['scale']
        x = (x-loc)/scale
        return 1.0/np.pi/(1.0+x*x) / scale
        
cauchy = cauchy_gen(name='cauchy',longname='Cauchy',extradoc="""

Cauchy distribution

cauchy.pdf(x) = 1/(pi*(1+x**2))

This is the t distribution with one degree of freedom.
""")
    
    
#GSASII peak fitting routine: Finger, Cox & Jephcoat model        


class fcjde_gen(st.rv_continuous):
    """
    Finger-Cox-Jephcoat D(2phi,2th) function for S/L = H/L
    Ref: J. Appl. Cryst. (1994) 27, 892-900.

    :param x: array -1 to 1
    :param t: 2-theta position of peak
    :param s: sum(S/L,H/L); S: sample height, H: detector opening, 
      L: sample to detector opening distance
    :param dx: 2-theta step size in deg

    :returns: for fcj.pdf

     * T = x*dx+t
     * s = S/L+H/L
     * if x < 0::

        fcj.pdf = [1/sqrt({cos(T)**2/cos(t)**2}-1) - 1/s]/|cos(T)| 

     * if x >= 0: fcj.pdf = 0    
    """
    def _pdf(self,x,t,s,dx):
        T = dx*x+t
        ax2 = abs(npcosd(T))
        ax = ax2**2
        bx = npcosd(t)**2
        bx = np.where(ax>bx,bx,ax)
        fx = np.where(ax>bx,(np.sqrt(bx/(ax-bx))-1./s)/ax2,0.0)
        fx = np.where(fx > 0.,fx,0.0)
        return fx
             
    def pdf(self,x,*args,**kwds):
        loc=kwds['loc']
        return self._pdf(x-loc,*args)
        
fcjde = fcjde_gen(name='fcjde',shapes='t,s,dx')
                
def getWidthsCW(pos,sig,gam,shl):
    '''Compute the peak widths used for computing the range of a peak
    for constant wavelength data. On low-angle side, 50 FWHM are used, 
    on high-angle side 75 are used, low angle side extended for axial divergence
    (for peaks above 90 deg, these are reversed.)
    '''
    widths = [np.sqrt(sig)/100.,gam/100.]
    fwhm = 2.355*widths[0]+widths[1]
    fmin = 50.*(fwhm+shl*abs(npcosd(pos)))
    fmax = 75.0*fwhm
    if pos > 90:
        fmin,fmax = [fmax,fmin]          
    return widths,fmin,fmax
    
def getWidthsTOF(pos,alp,bet,sig,gam):
    '''Compute the peak widths used for computing the range of a peak
    for constant wavelength data. 50 FWHM are used on both sides each 
    extended by exponential coeff.
    '''
    widths = [np.sqrt(sig),gam]
    fwhm = 2.355*widths[0]+2.*widths[1]
    fmin = 50.*fwhm*(1.+1./alp)    
    fmax = 50.*fwhm*(1.+1./bet)
    return widths,fmin,fmax
    
def getFWHM(pos,Inst):
    '''Compute total FWHM from Thompson, Cox & Hastings (1987) , J. Appl. Cryst. 20, 79-83
    via getgamFW(g,s).
    
    :param pos: float peak position in deg 2-theta or tof in musec
    :param Inst: dict instrument parameters
    
    :returns float: total FWHM of pseudoVoigt in deg or musec
    ''' 
    
    sig = lambda Th,U,V,W: np.sqrt(max(0.001,U*tand(Th)**2+V*tand(Th)+W))
    sigTOF = lambda dsp,S0,S1,S2,Sq: np.sqrt(S0+S1*dsp**2+S2*dsp**4+Sq*dsp)
    gam = lambda Th,X,Y,Z: Z+X/cosd(Th)+Y*tand(Th)
    gamTOF = lambda dsp,X,Y,Z: Z+X*dsp+Y*dsp**2
    alpTOF = lambda dsp,alp: alp/dsp
    betTOF = lambda dsp,bet0,bet1,betq: bet0+bet1/dsp**4+betq/dsp**2
    alpPink = lambda pos,alp0,alp1: alp0+alp1*tand(pos/2.)
    betPink = lambda pos,bet0,bet1: bet0+bet1*tand(pos/2.)
    if 'T' in Inst['Type'][0]:
        dsp = pos/Inst['difC'][1]
        alp = alpTOF(dsp,Inst['alpha'][1])
        bet = betTOF(dsp,Inst['beta-0'][1],Inst['beta-1'][1],Inst['beta-q'][1])
        s = sigTOF(dsp,Inst['sig-0'][1],Inst['sig-1'][1],Inst['sig-2'][1],Inst['sig-q'][1])
        g = gamTOF(dsp,Inst['X'][1],Inst['Y'][1],Inst['Z'][1])
        return getgamFW(g,s)+np.log(2.0)*(alp+bet)/(alp*bet)
    elif 'C' in Inst['Type'][0]:
        s = sig(pos/2.,Inst['U'][1],Inst['V'][1],Inst['W'][1])
        g = gam(pos/2.,Inst['X'][1],Inst['Y'][1],Inst['Z'][1])
        return getgamFW(g,s)/100.  #returns FWHM in deg
    else:   #'B'
        alp = alpPink(pos,Inst['alpha-0'][1],Inst['alpha-1'][1])
        bet = betPink(pos,Inst['beta-0'][1],Inst['beta-1'][1])
        s = sig(pos/2.,Inst['U'][1],Inst['V'][1],Inst['W'][1])
        g = gam(pos/2.,Inst['X'][1],Inst['Y'][1],Inst['Z'][1])
        return getgamFW(g,s)/100.+np.log(2.0)*(alp+bet)/(alp*bet)  #returns FWHM in deg
    
def getgamFW(g,s):
    '''Compute total FWHM from Thompson, Cox & Hastings (1987), J. Appl. Cryst. 20, 79-83
    lambda fxn needs FWHM for both Gaussian & Lorentzian components
    
    :param g: float Lorentzian gamma = FWHM(L)
    :param s: float Gaussian sig
    
    :returns float: total FWHM of pseudoVoigt
    ''' 
    gamFW = lambda s,g: np.exp(np.log(s**5+2.69269*s**4*g+2.42843*s**3*g**2+4.47163*s**2*g**3+0.07842*s*g**4+g**5)/5.)
    return gamFW(2.35482*s,g)   #sqrt(8ln2)*sig = FWHM(G)
                
def getFCJVoigt(pos,intens,sig,gam,shl,xdata):    
    '''Compute the Finger-Cox-Jepcoat modified Voigt function for a
    CW powder peak by direct convolution. This version is not used.
    '''
    DX = xdata[1]-xdata[0]
    widths,fmin,fmax = getWidthsCW(pos,sig,gam,shl)
    x = np.linspace(pos-fmin,pos+fmin,256)
    dx = x[1]-x[0]
    Norm = norm.pdf(x,loc=pos,scale=widths[0])
    Cauchy = cauchy.pdf(x,loc=pos,scale=widths[1])
    arg = [pos,shl/57.2958,dx,]
    FCJ = fcjde.pdf(x,*arg,loc=pos)
    if len(np.nonzero(FCJ)[0])>5:
        z = np.column_stack([Norm,Cauchy,FCJ]).T
        Z = fft.fft(z)
        Df = fft.ifft(Z.prod(axis=0)).real
    else:
        z = np.column_stack([Norm,Cauchy]).T
        Z = fft.fft(z)
        Df = fft.fftshift(fft.ifft(Z.prod(axis=0))).real
    Df /= np.sum(Df)
    Df = si.interp1d(x,Df,bounds_error=False,fill_value=0.0)
    return intens*Df(xdata)*DX/dx

def getBackground(pfx,parmDict,bakType,dataType,xdata,fixedBkg={}):
    '''Computes the background from vars pulled from gpx file or tree.
    '''
    if 'T' in dataType:
        q = 2.*np.pi*parmDict[pfx+'difC']/xdata
    else:
        wave = parmDict.get(pfx+'Lam',parmDict.get(pfx+'Lam1',1.0))
        q = npT2q(xdata,wave)
    yb = np.zeros_like(xdata)
    nBak = 0
    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])
    sumBk = [0.,0.,0]
    while True:
        key = pfx+'Back;'+str(nBak)
        if key in parmDict:
            nBak += 1
        else:
            break
#empirical functions
    if bakType in ['chebyschev','cosine','chebyschev-1']:
        dt = xdata[-1]-xdata[0]    
        for iBak in range(nBak):
            key = pfx+'Back;'+str(iBak)
            if bakType == 'chebyschev':
                ybi = parmDict[key]*(-1.+2.*(xdata-xdata[0])/dt)**iBak
            elif bakType == 'chebyschev-1':
                xpos = -1.+2.*(xdata-xdata[0])/dt
                ybi = parmDict[key]*np.cos(iBak*np.arccos(xpos))
            elif bakType == 'cosine':
                ybi = parmDict[key]*npcosd(180.*xdata*iBak/xdata[-1])
            yb += ybi
        sumBk[0] = np.sum(yb)
    elif bakType in ['Q^2 power series','Q^-2 power series']:
        QT = 1.
        yb += np.ones_like(yb)*parmDict[pfx+'Back;0']
        for iBak in range(nBak-1):
            key = pfx+'Back;'+str(iBak+1)
            if '-2' in bakType:
                QT *= (iBak+1)*q**-2
            else:
                QT *= q**2/(iBak+1)
            yb += QT*parmDict[key]
        sumBk[0] = np.sum(yb)
    elif bakType in ['lin interpolate','inv interpolate','log interpolate',]:
        if nBak == 1:
            yb = np.ones_like(xdata)*parmDict[pfx+'Back;0']
        elif nBak == 2:
            dX = xdata[-1]-xdata[0]
            T2 = (xdata-xdata[0])/dX
            T1 = 1.0-T2
            yb = parmDict[pfx+'Back;0']*T1+parmDict[pfx+'Back;1']*T2
        else:
            xnomask = ma.getdata(xdata)
            xmin,xmax = xnomask[0],xnomask[-1]
            if bakType == 'lin interpolate':
                bakPos = np.linspace(xmin,xmax,nBak,True)
            elif bakType == 'inv interpolate':
                bakPos = 1./np.linspace(1./xmax,1./xmin,nBak,True)
            elif bakType == 'log interpolate':
                bakPos = np.exp(np.linspace(np.log(xmin),np.log(xmax),nBak,True))
            bakPos[0] = xmin
            bakPos[-1] = xmax
            bakVals = np.zeros(nBak)
            for i in range(nBak):
                bakVals[i] = parmDict[pfx+'Back;'+str(i)]
            bakInt = si.interp1d(bakPos,bakVals,'linear')
            yb = bakInt(ma.getdata(xdata))
        sumBk[0] = np.sum(yb)
#Debye function        
    if pfx+'difC' in parmDict:
        ff = 1.
    else:        
        try:
            wave = parmDict[pfx+'Lam']
        except KeyError:
            wave = parmDict[pfx+'Lam1']
        SQ = (q/(4.*np.pi))**2
        FF = G2elem.GetFormFactorCoeff('Si')[0]
        ff = np.array(G2elem.ScatFac(FF,SQ)[0])**2
    iD = 0        
    while True:
        try:
            dbA = parmDict[pfx+'DebyeA;'+str(iD)]
            dbR = parmDict[pfx+'DebyeR;'+str(iD)]
            dbU = parmDict[pfx+'DebyeU;'+str(iD)]
            ybi = ff*dbA*np.sin(q*dbR)*np.exp(-dbU*q**2)/(q*dbR)
            yb += ybi
            sumBk[1] += np.sum(ybi)
            iD += 1       
        except KeyError:
            break
#peaks
    iD = 0
    while True:
        try:
            pkP = parmDict[pfx+'BkPkpos;'+str(iD)]
            pkI = max(parmDict[pfx+'BkPkint;'+str(iD)],0.1)
            pkS = max(parmDict[pfx+'BkPksig;'+str(iD)],1.)
            pkG = max(parmDict[pfx+'BkPkgam;'+str(iD)],0.1)
            if 'C' in dataType:
                Wd,fmin,fmax = getWidthsCW(pkP,pkS,pkG,.002)
            else: #'T'OF
                Wd,fmin,fmax = getWidthsTOF(pkP,1.,1.,pkS,pkG)
            iBeg = np.searchsorted(xdata,pkP-fmin)
            iFin = np.searchsorted(xdata,pkP+fmax)
            lenX = len(xdata)
            if not iBeg:
                iFin = np.searchsorted(xdata,pkP+fmax)
            elif iBeg == lenX:
                iFin = iBeg
            else:
                iFin = np.searchsorted(xdata,pkP+fmax)
            if 'C' in dataType:
                ybi = pkI*getFCJVoigt3(pkP,pkS,pkG,0.002,xdata[iBeg:iFin])
                yb[iBeg:iFin] += ybi
            elif 'T' in dataType:
                ybi = pkI*getEpsVoigt(pkP,1.,1.,pkS,pkG,xdata[iBeg:iFin])
                yb[iBeg:iFin] += ybi
            elif 'B' in dataType:
                ybi = pkI*getEpsVoigt(pkP,1.,1.,pkS/100.,pkG/1.e4,xdata[iBeg:iFin])
                yb[iBeg:iFin] += ybi
            sumBk[2] += np.sum(ybi)
            iD += 1       
        except KeyError:
            break
        except ValueError:
            G2fil.G2Print ('**** WARNING - backround peak '+str(iD)+' sigma is negative; fix & try again ****')
            break
    # fixed background from file
    if len(fixedBkg) >= 3:
        mult = fixedBkg.get('_fixedMult',0.0)
        if len(fixedBkg.get('_fixedValues',[])) != len(yb):
            G2fil.G2Print('Lengths of backgrounds do not agree: yb={}, fixed={}'.format(
                len(yb),len(fixedBkg.get('_fixedValues',[]))))
        elif mult: 
            yb -= mult*fixedBkg.get('_fixedValues',[]) # N.B. mult is negative
            sumBk[0] = sum(yb)
    return yb,sumBk
    
def getBackgroundDerv(hfx,parmDict,bakType,dataType,xdata):
    'needs a doc string'
    if 'T' in dataType:
        q = 2.*np.pi*parmDict[hfx+'difC']/xdata
    else:
        wave = parmDict.get(hfx+'Lam',parmDict.get(hfx+'Lam1',1.0))
        q = 2.*np.pi*npsind(xdata/2.)/wave
    nBak = 0
    while True:
        key = hfx+'Back;'+str(nBak)
        if key in parmDict:
            nBak += 1
        else:
            break
    dydb = np.zeros(shape=(nBak,len(xdata)))
    dyddb = np.zeros(shape=(3*parmDict[hfx+'nDebye'],len(xdata)))
    dydpk = np.zeros(shape=(4*parmDict[hfx+'nPeaks'],len(xdata)))
    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])

    if bakType in ['chebyschev','cosine','chebyschev-1']:
        dt = xdata[-1]-xdata[0]    
        for iBak in range(nBak):    
            if bakType == 'chebyschev':
                dydb[iBak] = (-1.+2.*(xdata-xdata[0])/dt)**iBak
            elif bakType == 'chebyschev-1':
                xpos = -1.+2.*(xdata-xdata[0])/dt
                dydb[iBak] = np.cos(iBak*np.arccos(xpos))
            elif bakType == 'cosine':
                dydb[iBak] = npcosd(180.*xdata*iBak/xdata[-1])
    elif bakType in ['Q^2 power series','Q^-2 power series']:
        QT = 1.
        dydb[0] = np.ones_like(xdata)
        for iBak in range(nBak-1):
            if '-2' in bakType:
                QT *= (iBak+1)*q**-2
            else:
                QT *= q**2/(iBak+1)
            dydb[iBak+1] = QT
    elif bakType in ['lin interpolate','inv interpolate','log interpolate',]:
        if nBak == 1:
            dydb[0] = np.ones_like(xdata)
        elif nBak == 2:
            dX = xdata[-1]-xdata[0]
            T2 = (xdata-xdata[0])/dX
            T1 = 1.0-T2
            dydb = [T1,T2]
        else:
            xnomask = ma.getdata(xdata)
            xmin,xmax = xnomask[0],xnomask[-1]
            if bakType == 'lin interpolate':
                bakPos = np.linspace(xmin,xmax,nBak,True)
            elif bakType == 'inv interpolate':
                bakPos = 1./np.linspace(1./xmax,1./xmin,nBak,True)
            elif bakType == 'log interpolate':
                bakPos = np.exp(np.linspace(np.log(xmin),np.log(xmax),nBak,True))
            bakPos[0] = xmin
            bakPos[-1] = xmax
            for i,pos in enumerate(bakPos):
                if i == 0:
                    dydb[0] = np.where(xdata<bakPos[1],(bakPos[1]-xdata)/(bakPos[1]-bakPos[0]),0.)
                elif i == len(bakPos)-1:
                    dydb[i] = np.where(xdata>bakPos[-2],(bakPos[-1]-xdata)/(bakPos[-1]-bakPos[-2]),0.)
                else:
                    dydb[i] = np.where(xdata>bakPos[i],
                        np.where(xdata<bakPos[i+1],(bakPos[i+1]-xdata)/(bakPos[i+1]-bakPos[i]),0.),
                        np.where(xdata>bakPos[i-1],(xdata-bakPos[i-1])/(bakPos[i]-bakPos[i-1]),0.))
    if hfx+'difC' in parmDict:
        ff = 1.
    else:
        wave = parmDict.get(hfx+'Lam',parmDict.get(hfx+'Lam1',1.0))
        q = npT2q(xdata,wave)
        SQ = (q/(4*np.pi))**2
        FF = G2elem.GetFormFactorCoeff('Si')[0]
        ff = np.array(G2elem.ScatFac(FF,SQ)[0])*np.pi**2    #needs pi^2~10. for cw data (why?)
    iD = 0        
    while True:
        try:
            if hfx+'difC' in parmDict:
                q = 2*np.pi*parmDict[hfx+'difC']/xdata
            dbA = parmDict[hfx+'DebyeA;'+str(iD)]
            dbR = parmDict[hfx+'DebyeR;'+str(iD)]
            dbU = parmDict[hfx+'DebyeU;'+str(iD)]
            sqr = np.sin(q*dbR)/(q*dbR)
            cqr = np.cos(q*dbR)
            temp = np.exp(-dbU*q**2)
            dyddb[3*iD] = ff*sqr*temp
            dyddb[3*iD+1] = ff*dbA*temp*(cqr-sqr)/(dbR)
            dyddb[3*iD+2] = -ff*dbA*sqr*temp*q**2
            iD += 1
        except KeyError:
            break
    iD = 0
    while True:
        try:
            pkP = parmDict[hfx+'BkPkpos;'+str(iD)]
            pkI = max(parmDict[hfx+'BkPkint;'+str(iD)],0.1)
            pkS = max(parmDict[hfx+'BkPksig;'+str(iD)],1.0)
            pkG = max(parmDict[hfx+'BkPkgam;'+str(iD)],0.1)
            if 'C' in dataType:
                Wd,fmin,fmax = getWidthsCW(pkP,pkS,pkG,.002)
            else: #'T' or 'B'
                Wd,fmin,fmax = getWidthsTOF(pkP,1.,1.,pkS,pkG)
            iBeg = np.searchsorted(xdata,pkP-fmin)
            iFin = np.searchsorted(xdata,pkP+fmax)
            lenX = len(xdata)
            if not iBeg:
                iFin = np.searchsorted(xdata,pkP+fmax)
            elif iBeg == lenX:
                iFin = iBeg
            else:
                iFin = np.searchsorted(xdata,pkP+fmax)
            if 'C' in dataType:
                Df,dFdp,dFds,dFdg,x = getdFCJVoigt3(pkP,pkS,pkG,.002,xdata[iBeg:iFin])
                dydpk[4*iD][iBeg:iFin] += 100.*cw[iBeg:iFin]*pkI*dFdp
                dydpk[4*iD+1][iBeg:iFin] += 100.*cw[iBeg:iFin]*Df
                dydpk[4*iD+2][iBeg:iFin] += 100.*cw[iBeg:iFin]*pkI*dFds
                dydpk[4*iD+3][iBeg:iFin] += 100.*cw[iBeg:iFin]*pkI*dFdg
            else:   #'T'OF
                Df,dFdp,x,x,dFds,dFdg = getdEpsVoigt(pkP,1.,1.,pkS,pkG,xdata[iBeg:iFin])
                dydpk[4*iD][iBeg:iFin] += pkI*dFdp
                dydpk[4*iD+1][iBeg:iFin] += Df
                dydpk[4*iD+2][iBeg:iFin] += pkI*dFds
                dydpk[4*iD+3][iBeg:iFin] += pkI*dFdg
            iD += 1       
        except KeyError:
            break
        except ValueError:
            G2fil.G2Print ('**** WARNING - backround peak '+str(iD)+' sigma is negative; fix & try again ****')
            break        
    return dydb,dyddb,dydpk

#use old fortran routine
def getFCJVoigt3(pos,sig,gam,shl,xdata):
    '''Compute the Finger-Cox-Jepcoat modified Pseudo-Voigt function for a
    CW powder peak in external Fortran routine
    '''
    Df = pyd.pypsvfcj(len(xdata),xdata-pos,pos,sig,gam,shl)
#    Df = pyd.pypsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)
    Df /= np.sum(Df)
    return Df

def getdFCJVoigt3(pos,sig,gam,shl,xdata):
    '''Compute analytic derivatives the Finger-Cox-Jepcoat modified Pseudo-Voigt
    function for a CW powder peak
    '''
    Df,dFdp,dFds,dFdg,dFdsh = pyd.pydpsvfcj(len(xdata),xdata-pos,pos,sig,gam,shl)
#    Df,dFdp,dFds,dFdg,dFdsh = pyd.pydpsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)
    return Df,dFdp,dFds,dFdg,dFdsh

def getPsVoigt(pos,sig,gam,xdata):
    'needs a doc string'
    
    Df = pyd.pypsvoigt(len(xdata),xdata-pos,sig,gam)
    Df /= np.sum(Df)
    return Df

def getdPsVoigt(pos,sig,gam,xdata):
    'needs a doc string'
    
    Df,dFdp,dFds,dFdg = pyd.pydpsvoigt(len(xdata),xdata-pos,sig,gam)
    return Df,dFdp,dFds,dFdg

def getEpsVoigt(pos,alp,bet,sig,gam,xdata):
    'needs a doc string'
    Df = pyd.pyepsvoigt(len(xdata),xdata-pos,alp,bet,sig,gam)
    Df /= np.sum(Df)
    return Df  
    
def getdEpsVoigt(pos,alp,bet,sig,gam,xdata):
    'needs a doc string'
    Df,dFdp,dFda,dFdb,dFds,dFdg = pyd.pydepsvoigt(len(xdata),xdata-pos,alp,bet,sig,gam)
    return Df,dFdp,dFda,dFdb,dFds,dFdg   

def ellipseSize(H,Sij,GB):
    'Implements r=1/sqrt(sum((1/S)*(q.v)^2) per note from Alexander Brady'
    HX = np.inner(H.T,GB)
    lenHX = np.sqrt(np.sum(HX**2))
    Esize,Rsize = nl.eigh(G2lat.U6toUij(Sij))            
    R = np.inner(HX/lenHX,Rsize)**2*Esize         #want column length for hkl in crystal
    lenR = 1./np.sqrt(np.sum(R))
    return lenR

def ellipseSizeDerv(H,Sij,GB):
    'needs a doc string'
    lenR = ellipseSize(H,Sij,GB)
    delt = 0.001
    dRdS = np.zeros(6)
    for i in range(6):
        Sij[i] -= delt
        lenM = ellipseSize(H,Sij,GB)
        Sij[i] += 2.*delt
        lenP = ellipseSize(H,Sij,GB)
        Sij[i] -= delt
        dRdS[i] = (lenP-lenM)/(2.*delt)
    return lenR,dRdS

def getMustrain(HKL,G,SGData,muStrData):
    if muStrData[0] == 'isotropic':
        return np.ones(HKL.shape[1])*muStrData[1][0]
    elif muStrData[0] == 'uniaxial':
        H = np.array(HKL)
        P = np.array(muStrData[3])
        cosP,sinP = np.array([G2lat.CosSinAngle(h,P,G) for h in H.T]).T
        Si = muStrData[1][0]
        Sa = muStrData[1][1]
        return Si*Sa/(np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
    else:       #generalized - P.W. Stephens model
        H = np.array(HKL)
        rdsq = np.array([G2lat.calc_rDsq2(h,G) for h in H.T])
        Strms = np.array(G2spc.MustrainCoeff(H,SGData))
        Sum = np.sum(np.array(muStrData[4])[:,nxs]*Strms,axis=0)
        return np.sqrt(Sum)/rdsq
    
def getCrSize(HKL,G,GB,sizeData):
    if sizeData[0] == 'isotropic':
        return np.ones(HKL.shape[1])*sizeData[1][0]
    elif sizeData[0] == 'uniaxial':
        H = np.array(HKL)
        P = np.array(sizeData[3])
        cosP,sinP = np.array([G2lat.CosSinAngle(h,P,G) for h in H.T]).T
        Si = sizeData[1][0]
        Sa = sizeData[1][1]
        return Si*Sa/(np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
    else:
        Sij =[sizeData[4][i] for i in range(6)]
        H = np.array(HKL)
        return 1./np.array([ellipseSize(h,Sij,GB) for h in H.T])**2

def getHKLpeak(dmin,SGData,A,Inst=None,nodup=False):
    '''
    Generates allowed by symmetry reflections with d >= dmin
    NB: GenHKLf & checkMagextc return True for extinct reflections

    :param dmin:  minimum d-spacing 
    :param SGData: space group data obtained from SpcGroup
    :param A: lattice parameter terms A1-A6
    :param Inst: instrument parameter info
    :returns: HKLs: np.array hkl, etc for allowed reflections

    '''
    HKL = G2lat.GenHLaue(dmin,SGData,A)        
    HKLs = []
    ds = []
    for h,k,l,d in HKL:
        ext = G2spc.GenHKLf([h,k,l],SGData)[0]
        if ext and 'MagSpGrp' in SGData:
            ext = G2spc.checkMagextc([h,k,l],SGData)
        if not ext:
            if nodup and int(10000*d) in ds:
                continue
            ds.append(int(10000*d))
            if Inst == None:
                HKLs.append([h,k,l,d,0,-1])
            else:
                HKLs.append([h,k,l,d,G2lat.Dsp2pos(Inst,d),-1])
    return np.array(HKLs)

def getHKLMpeak(dmin,Inst,SGData,SSGData,Vec,maxH,A):
    'needs a doc string'
    HKLs = []
    vec = np.array(Vec)
    vstar = np.sqrt(G2lat.calc_rDsq(vec,A))     #find extra needed for -n SS reflections
    dvec = 1./(maxH*vstar+1./dmin)
    HKL = G2lat.GenHLaue(dvec,SGData,A)        
    SSdH = [vec*h for h in range(-maxH,maxH+1)]
    SSdH = dict(zip(range(-maxH,maxH+1),SSdH))
    ifMag = False
    if 'MagSpGrp' in SGData:
        ifMag = True
    for h,k,l,d in HKL:
        ext = G2spc.GenHKLf([h,k,l],SGData)[0]
        if not ext and d >= dmin:
            HKLs.append([h,k,l,0,d,G2lat.Dsp2pos(Inst,d),-1])
        for dH in SSdH:
            if dH:
                DH = SSdH[dH]
                H = [h+DH[0],k+DH[1],l+DH[2]]
                d = float(1/np.sqrt(G2lat.calc_rDsq(H,A)))
                if d >= dmin:
                    HKLM = np.array([h,k,l,dH])
                    if G2spc.checkSSextc(HKLM,SSGData) or ifMag:
                        HKLs.append([h,k,l,dH,d,G2lat.Dsp2pos(Inst,d),-1])    
    return G2lat.sortHKLd(HKLs,True,True,True)

def getPeakProfile(dataType,parmDict,xdata,varyList,bakType):
    'Computes the profile for a powder pattern'
    
    yb = getBackground('',parmDict,bakType,dataType,xdata)[0]
    yc = np.zeros_like(yb)
    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])
    if 'C' in dataType:
        shl = max(parmDict['SH/L'],0.002)
        Ka2 = False
        if 'Lam1' in parmDict.keys():
            Ka2 = True
            lamRatio = 360*(parmDict['Lam2']-parmDict['Lam1'])/(np.pi*parmDict['Lam1'])
            kRatio = parmDict['I(L2)/I(L1)']
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                tth = (pos-parmDict['Zero'])
                intens = parmDict['int'+str(iPeak)]
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                sig = max(sig,0.001)          #avoid neg sigma^2
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsCW(pos,sig,gam,shl)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    return yb+yc
                yc[iBeg:iFin] += intens*getFCJVoigt3(pos,sig,gam,shl,xdata[iBeg:iFin])
                if Ka2:
                    pos2 = pos+lamRatio*tand(pos/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    iBeg = np.searchsorted(xdata,pos2-fmin)
                    iFin = np.searchsorted(xdata,pos2+fmin)
                    if iBeg-iFin:
                        yc[iBeg:iFin] += intens*kRatio*getFCJVoigt3(pos2,sig,gam,shl,xdata[iBeg:iFin])
                iPeak += 1
            except KeyError:        #no more peaks to process
                return yb+yc
    elif 'B' in dataType:
        iPeak = 0
        dsp = 1.0 #for now - fix later
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                tth = (pos-parmDict['Zero'])
                intens = parmDict['int'+str(iPeak)]
                alpName = 'alp'+str(iPeak)
                if alpName in varyList:
                    alp = parmDict[alpName]
                else:
                    alp = G2mth.getPinkalpha(parmDict,tth)
                alp = max(0.1,alp)
                betName = 'bet'+str(iPeak)
                if betName in varyList:
                    bet = parmDict[betName]
                else:
                    bet = G2mth.getPinkbeta(parmDict,tth)
                bet = max(0.0001,bet)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                sig = max(sig,0.001)          #avoid neg sigma^2
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsTOF(pos,alp,bet,sig,gam)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    return yb+yc
                yc[iBeg:iFin] += intens*getEpsVoigt(pos,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])
                iPeak += 1
            except KeyError:        #no more peaks to process
                return yb+yc        
    else:
        Pdabc = parmDict['Pdabc']
        difC = parmDict['difC']
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]                
                tof = pos-parmDict['Zero']
                dsp = tof/difC
                intens = parmDict['int'+str(iPeak)]
                alpName = 'alp'+str(iPeak)
                if alpName in varyList:
                    alp = parmDict[alpName]
                else:
                    if len(Pdabc):
                        alp = np.interp(dsp,Pdabc[0],Pdabc[1])
                    else:
                        alp = G2mth.getTOFalpha(parmDict,dsp)
                alp = max(0.1,alp)
                betName = 'bet'+str(iPeak)
                if betName in varyList:
                    bet = parmDict[betName]
                else:
                    if len(Pdabc):
                        bet = np.interp(dsp,Pdabc[0],Pdabc[2])
                    else:
                        bet = G2mth.getTOFbeta(parmDict,dsp)
                bet = max(0.0001,bet)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getTOFsig(parmDict,dsp)
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                else:
                    gam = G2mth.getTOFgamma(parmDict,dsp)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsTOF(pos,alp,bet,sig,gam)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmax)
                lenX = len(xdata)
                if not iBeg:
                    iFin = np.searchsorted(xdata,pos+fmax)
                elif iBeg == lenX:
                    iFin = iBeg
                else:
                    iFin = np.searchsorted(xdata,pos+fmax)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    return yb+yc
                yc[iBeg:iFin] += intens*getEpsVoigt(pos,alp,bet,sig,gam,xdata[iBeg:iFin])
                iPeak += 1
            except KeyError:        #no more peaks to process
                return yb+yc
            
def getPeakProfileDerv(dataType,parmDict,xdata,varyList,bakType):
    'needs a doc string'
# needs to return np.array([dMdx1,dMdx2,...]) in same order as varylist = backVary,insVary,peakVary order
    dMdv = np.zeros(shape=(len(varyList),len(xdata)))
    dMdb,dMddb,dMdpk = getBackgroundDerv('',parmDict,bakType,dataType,xdata)
    if 'Back;0' in varyList:            #background derivs are in front if present
        dMdv[0:len(dMdb)] = dMdb
    names = ['DebyeA','DebyeR','DebyeU']
    for name in varyList:
        if 'Debye' in name:
            parm,Id = name.split(';')
            ip = names.index(parm)
            dMdv[varyList.index(name)] = dMddb[3*int(Id)+ip]
    names = ['BkPkpos','BkPkint','BkPksig','BkPkgam']
    for name in varyList:
        if 'BkPk' in name:
            parm,Id = name.split(';')
            ip = names.index(parm)
            dMdv[varyList.index(name)] = dMdpk[4*int(Id)+ip]
    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])
    if 'C' in dataType:
        shl = max(parmDict['SH/L'],0.002)
        Ka2 = False
        if 'Lam1' in parmDict.keys():
            Ka2 = True
            lamRatio = 360*(parmDict['Lam2']-parmDict['Lam1'])/(np.pi*parmDict['Lam1'])
            kRatio = parmDict['I(L2)/I(L1)']
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                tth = (pos-parmDict['Zero'])
                intens = parmDict['int'+str(iPeak)]
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                    dsdU = dsdV = dsdW = 0
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                    dsdU,dsdV,dsdW = G2mth.getCWsigDeriv(tth)
                sig = max(sig,0.001)          #avoid neg sigma
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                    dgdX = dgdY = dgdZ = 0
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                    dgdX,dgdY,dgdZ = G2mth.getCWgamDeriv(tth)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsCW(pos,sig,gam,shl)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    break
                dMdpk = np.zeros(shape=(6,len(xdata)))
                dMdipk = getdFCJVoigt3(pos,sig,gam,shl,xdata[iBeg:iFin])
                for i in range(1,5):
                    dMdpk[i][iBeg:iFin] += 100.*cw[iBeg:iFin]*intens*dMdipk[i]
                dMdpk[0][iBeg:iFin] += 100.*cw[iBeg:iFin]*dMdipk[0]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4]}
                if Ka2:
                    pos2 = pos+lamRatio*tand(pos/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    iBeg = np.searchsorted(xdata,pos2-fmin)
                    iFin = np.searchsorted(xdata,pos2+fmin)
                    if iBeg-iFin:
                        dMdipk2 = getdFCJVoigt3(pos2,sig,gam,shl,xdata[iBeg:iFin])
                        for i in range(1,5):
                            dMdpk[i][iBeg:iFin] += 100.*cw[iBeg:iFin]*intens*kRatio*dMdipk2[i]
                        dMdpk[0][iBeg:iFin] += 100.*cw[iBeg:iFin]*kRatio*dMdipk2[0]
                        dMdpk[5][iBeg:iFin] += 100.*cw[iBeg:iFin]*dMdipk2[0]
                        dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':dMdpk[5]*intens}
                for parmName in ['pos','int','sig','gam']:
                    try:
                        idx = varyList.index(parmName+str(iPeak))
                        dMdv[idx] = dervDict[parmName]
                    except ValueError:
                        pass
                if 'U' in varyList:
                    dMdv[varyList.index('U')] += dsdU*dervDict['sig']
                if 'V' in varyList:
                    dMdv[varyList.index('V')] += dsdV*dervDict['sig']
                if 'W' in varyList:
                    dMdv[varyList.index('W')] += dsdW*dervDict['sig']
                if 'X' in varyList:
                    dMdv[varyList.index('X')] += dgdX*dervDict['gam']
                if 'Y' in varyList:
                    dMdv[varyList.index('Y')] += dgdY*dervDict['gam']
                if 'Z' in varyList:
                    dMdv[varyList.index('Z')] += dgdZ*dervDict['gam']
                if 'SH/L' in varyList:
                    dMdv[varyList.index('SH/L')] += dervDict['shl']         #problem here
                if 'I(L2)/I(L1)' in varyList:
                    dMdv[varyList.index('I(L2)/I(L1)')] += dervDict['L1/L2']
                iPeak += 1
            except KeyError:        #no more peaks to process
                break
    elif 'B' in dataType:
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                tth = (pos-parmDict['Zero'])
                intens = parmDict['int'+str(iPeak)]
                alpName = 'alp'+str(iPeak)
                if alpName in varyList:
                    alp = parmDict[alpName]
                    dada0 = dada1 = 0.0
                else:
                    alp = G2mth.getPinkalpha(parmDict,tth)
                    dada0,dada1 = G2mth.getPinkalphaDeriv(tth)
                alp = max(0.0001,alp)
                betName = 'bet'+str(iPeak)
                if betName in varyList:
                    bet = parmDict[betName]
                    dbdb0 = dbdb1 = 0.0
                else:
                    bet = G2mth.getPinkbeta(parmDict,tth)
                    dbdb0,dbdb1 = G2mth.getPinkbetaDeriv(tth)
                bet = max(0.0001,bet)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                    dsdU = dsdV = dsdW = 0
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                    dsdU,dsdV,dsdW = G2mth.getCWsigDeriv(tth)
                sig = max(sig,0.001)          #avoid neg sigma
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                    dgdX = dgdY = dgdZ = 0
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                    dgdX,dgdY,dgdZ = G2mth.getCWgamDeriv(tth)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsTOF(pos,alp,bet,sig/1.e4,gam/100.)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    break
                dMdpk = np.zeros(shape=(7,len(xdata)))
                dMdipk = getdEpsVoigt(pos,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])
                for i in range(1,6):
                    dMdpk[i][iBeg:iFin] += cw[iBeg:iFin]*intens*dMdipk[i]
                dMdpk[0][iBeg:iFin] += cw[iBeg:iFin]*dMdipk[0]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4]/1.e4,'gam':dMdpk[5]/100.}
                for parmName in ['pos','int','alp','bet','sig','gam']:
                    try:
                        idx = varyList.index(parmName+str(iPeak))
                        dMdv[idx] = dervDict[parmName]
                    except ValueError:
                        pass
                if 'U' in varyList:
                    dMdv[varyList.index('U')] += dsdU*dervDict['sig']
                if 'V' in varyList:
                    dMdv[varyList.index('V')] += dsdV*dervDict['sig']
                if 'W' in varyList:
                    dMdv[varyList.index('W')] += dsdW*dervDict['sig']
                if 'X' in varyList:
                    dMdv[varyList.index('X')] += dgdX*dervDict['gam']
                if 'Y' in varyList:
                    dMdv[varyList.index('Y')] += dgdY*dervDict['gam']
                if 'Z' in varyList:
                    dMdv[varyList.index('Z')] += dgdZ*dervDict['gam']
                if 'alpha-0' in varyList:
                    dMdv[varyList.index('alpha-0')] += dada0*dervDict['alp']
                if 'alpha-1' in varyList:
                    dMdv[varyList.index('alpha-1')] += dada1*dervDict['alp']
                if 'beta-0' in varyList:
                    dMdv[varyList.index('beta-0')] += dbdb0*dervDict['bet']
                if 'beta-1' in varyList:
                    dMdv[varyList.index('beta-1')] += dbdb1*dervDict['bet']
                iPeak += 1
            except KeyError:        #no more peaks to process
                break        
    else:
        Pdabc = parmDict['Pdabc']
        difC = parmDict['difC']
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]                
                tof = pos-parmDict['Zero']
                dsp = tof/difC
                intens = parmDict['int'+str(iPeak)]
                alpName = 'alp'+str(iPeak)
                if alpName in varyList:
                    alp = parmDict[alpName]
                else:
                    if len(Pdabc):
                        alp = np.interp(dsp,Pdabc[0],Pdabc[1])
                        dada0 = 0
                    else:
                        alp = G2mth.getTOFalpha(parmDict,dsp)
                        dada0 = G2mth.getTOFalphaDeriv(dsp)
                betName = 'bet'+str(iPeak)
                if betName in varyList:
                    bet = parmDict[betName]
                else:
                    if len(Pdabc):
                        bet = np.interp(dsp,Pdabc[0],Pdabc[2])
                        dbdb0 = dbdb1 = dbdb2 = 0
                    else:
                        bet = G2mth.getTOFbeta(parmDict,dsp)
                        dbdb0,dbdb1,dbdb2 = G2mth.getTOFbetaDeriv(dsp)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList:
                    sig = parmDict[sigName]
                    dsds0 = dsds1 = dsds2 = dsds3 = 0
                else:
                    sig = G2mth.getTOFsig(parmDict,dsp)
                    dsds0,dsds1,dsds2,dsds3 = G2mth.getTOFsigDeriv(dsp)
                gamName = 'gam'+str(iPeak)
                if gamName in varyList:
                    gam = parmDict[gamName]
                    dsdX = dsdY = dsdZ = 0
                else:
                    gam = G2mth.getTOFgamma(parmDict,dsp)
                    dsdX,dsdY,dsdZ = G2mth.getTOFgammaDeriv(dsp)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsTOF(pos,alp,bet,sig,gam)
                iBeg = np.searchsorted(xdata,pos-fmin)
                lenX = len(xdata)
                if not iBeg:
                    iFin = np.searchsorted(xdata,pos+fmax)
                elif iBeg == lenX:
                    iFin = iBeg
                else:
                    iFin = np.searchsorted(xdata,pos+fmax)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    break
                dMdpk = np.zeros(shape=(7,len(xdata)))
                dMdipk = getdEpsVoigt(pos,alp,bet,sig,gam,xdata[iBeg:iFin])
                for i in range(1,6):
                    dMdpk[i][iBeg:iFin] += intens*cw[iBeg:iFin]*dMdipk[i]
                dMdpk[0][iBeg:iFin] += cw[iBeg:iFin]*dMdipk[0]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5]}
                for parmName in ['pos','int','alp','bet','sig','gam']:
                    try:
                        idx = varyList.index(parmName+str(iPeak))
                        dMdv[idx] = dervDict[parmName]
                    except ValueError:
                        pass
                if 'alpha' in varyList:
                    dMdv[varyList.index('alpha')] += dada0*dervDict['alp']
                if 'beta-0' in varyList:
                    dMdv[varyList.index('beta-0')] += dbdb0*dervDict['bet']
                if 'beta-1' in varyList:
                    dMdv[varyList.index('beta-1')] += dbdb1*dervDict['bet']
                if 'beta-q' in varyList:
                    dMdv[varyList.index('beta-q')] += dbdb2*dervDict['bet']
                if 'sig-0' in varyList:
                    dMdv[varyList.index('sig-0')] += dsds0*dervDict['sig']
                if 'sig-1' in varyList:
                    dMdv[varyList.index('sig-1')] += dsds1*dervDict['sig']
                if 'sig-2' in varyList:
                    dMdv[varyList.index('sig-2')] += dsds2*dervDict['sig']
                if 'sig-q' in varyList:
                    dMdv[varyList.index('sig-q')] += dsds3*dervDict['sig']
                if 'X' in varyList:
                    dMdv[varyList.index('X')] += dsdX*dervDict['gam']
                if 'Y' in varyList:
                    dMdv[varyList.index('Y')] += dsdY*dervDict['gam']
                if 'Z' in varyList:
                    dMdv[varyList.index('Z')] += dsdZ*dervDict['gam']
                iPeak += 1
            except KeyError:        #no more peaks to process
                break
    return dMdv
        
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def SetBackgroundParms(Background):
    'Loads background parameters into dicts/lists to create varylist & parmdict'
    if len(Background) == 1:            # fix up old backgrounds
        Background.append({'nDebye':0,'debyeTerms':[]})
    bakType,bakFlag = Background[0][:2]
    backVals = Background[0][3:]
    backNames = ['Back;'+str(i) for i in range(len(backVals))]
    Debye = Background[1]           #also has background peaks stuff
    backDict = dict(zip(backNames,backVals))
    backVary = []
    if bakFlag:
        backVary = backNames

    backDict['nDebye'] = Debye['nDebye']
    debyeDict = {}
    debyeList = []
    for i in range(Debye['nDebye']):
        debyeNames = ['DebyeA;'+str(i),'DebyeR;'+str(i),'DebyeU;'+str(i)]
        debyeDict.update(dict(zip(debyeNames,Debye['debyeTerms'][i][::2])))
        debyeList += zip(debyeNames,Debye['debyeTerms'][i][1::2])
    debyeVary = []
    for item in debyeList:
        if item[1]:
            debyeVary.append(item[0])
    backDict.update(debyeDict)
    backVary += debyeVary 

    backDict['nPeaks'] = Debye['nPeaks']
    peaksDict = {}
    peaksList = []
    for i in range(Debye['nPeaks']):
        peaksNames = ['BkPkpos;'+str(i),'BkPkint;'+str(i),'BkPksig;'+str(i),'BkPkgam;'+str(i)]
        peaksDict.update(dict(zip(peaksNames,Debye['peaksList'][i][::2])))
        peaksList += zip(peaksNames,Debye['peaksList'][i][1::2])
    peaksVary = []
    for item in peaksList:
        if item[1]:
            peaksVary.append(item[0])
    backDict.update(peaksDict)
    backVary += peaksVary
    return bakType,backDict,backVary
    
def DoCalibInst(IndexPeaks,Inst):
    
    def SetInstParms():
        dataType = Inst['Type'][0]
        insVary = []
        insNames = []
        insVals = []
        for parm in Inst:
            insNames.append(parm)
            insVals.append(Inst[parm][1])
            if parm in ['Lam','difC','difA','difB','Zero',]:
                if Inst[parm][2]:
                    insVary.append(parm)
        instDict = dict(zip(insNames,insVals))
        return dataType,instDict,insVary
        
    def GetInstParms(parmDict,Inst,varyList):
        for name in Inst:
            Inst[name][1] = parmDict[name]
        
    def InstPrint(Inst,sigDict):
        print ('Instrument Parameters:')
        if 'C' in Inst['Type'][0] or 'B' in Inst['Type'][0]:
            ptfmt = "%12.6f"
        else:
            ptfmt = "%12.3f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for parm in Inst:
            if parm in  ['Lam','difC','difA','difB','Zero',]:
                ptlbls += "%s" % (parm.center(12))
                ptstr += ptfmt % (Inst[parm][1])
                if parm in sigDict:
                    sigstr += ptfmt % (sigDict[parm])
                else:
                    sigstr += 12*' '
        print (ptlbls)
        print (ptstr)
        print (sigstr)
        
    def errPeakPos(values,peakDsp,peakPos,peakWt,dataType,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        return np.sqrt(peakWt)*(G2lat.getPeakPos(dataType,parmDict,peakDsp)-peakPos)

    peakPos = []
    peakDsp = []
    peakWt = []
    for peak,sig in zip(IndexPeaks[0],IndexPeaks[1]):
        if peak[2] and peak[3] and sig > 0.:
            peakPos.append(peak[0])
            peakDsp.append(peak[-1])    #d-calc
#            peakWt.append(peak[-1]**2/sig**2)   #weight by d**2
            peakWt.append(1./(sig*peak[-1]))   #
    peakPos = np.array(peakPos)
    peakDsp = np.array(peakDsp)
    peakWt = np.array(peakWt)
    dataType,insDict,insVary = SetInstParms()
    parmDict = {}
    parmDict.update(insDict)
    varyList = insVary
    if not len(varyList):
        G2fil.G2Print ('**** ERROR - nothing to refine! ****')
        return False
    while True:
        begin = time.time()
        values =  np.array(Dict2Values(parmDict, varyList))
        result = so.leastsq(errPeakPos,values,full_output=True,ftol=0.000001,
            args=(peakDsp,peakPos,peakWt,dataType,parmDict,varyList))
        ncyc = int(result[2]['nfev']/2)
        runtime = time.time()-begin    
        chisq = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
        GOF = chisq/(len(peakPos)-len(varyList))       #reduced chi^2
        G2fil.G2Print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(result[2]['nfev'],len(peakPos),len(varyList)))
        G2fil.G2Print ('calib time = %8.3fs, %8.3fs/cycle'%(runtime,runtime/ncyc))
        G2fil.G2Print ('chi**2 = %12.6g, reduced chi**2 = %6.2f'%(chisq,GOF))
        try:
            sig = np.sqrt(np.diag(result[1])*GOF)
            if np.any(np.isnan(sig)):
                G2fil.G2Print ('*** Least squares aborted - some invalid esds possible ***')
            break                   #refinement succeeded - finish up!
        except ValueError:          #result[1] is None on singular matrix
            G2fil.G2Print ('**** Refinement failed - singular matrix ****')
        
    sigDict = dict(zip(varyList,sig))
    GetInstParms(parmDict,Inst,varyList)
    InstPrint(Inst,sigDict)
    return True
            
def DoPeakFit(FitPgm,Peaks,Background,Limits,Inst,Inst2,data,fixback=None,prevVaryList=[],oneCycle=False,controls=None,wtFactor=1.0,dlg=None):
    '''Called to perform a peak fit, refining the selected items in the peak
    table as well as selected items in the background.

    :param str FitPgm: type of fit to perform. At present this is ignored.
    :param list Peaks: a list of peaks. Each peak entry is a list with 8 values:
      four values followed by a refine flag where the values are: position, intensity,
      sigma (Gaussian width) and gamma (Lorentzian width). From the Histogram/"Peak List"
      tree entry, dict item "peaks"
    :param list Background: describes the background. List with two items.
      Item 0 specifies a background model and coefficients. Item 1 is a dict.
      From the Histogram/Background tree entry.
    :param list Limits: min and max x-value to use
    :param dict Inst: Instrument parameters
    :param dict Inst2: more Instrument parameters
    :param numpy.array data: a 5xn array. data[0] is the x-values,
      data[1] is the y-values, data[2] are weight values, data[3], [4] and [5] are
      calc, background and difference intensities, respectively. 
    :param array fixback: fixed background values
    :param list prevVaryList: Used in sequential refinements to override the
      variable list. Defaults as an empty list.
    :param bool oneCycle: True if only one cycle of fitting should be performed
    :param dict controls: a dict specifying two values, Ftol = controls['min dM/M']
      and derivType = controls['deriv type']. If None default values are used. 
    :param float wtFactor: weight multiplier; = 1.0 by default
    :param wx.Dialog dlg: A dialog box that is updated with progress from the fit.
      Defaults to None, which means no updates are done.
    '''
    def GetBackgroundParms(parmList,Background):
        iBak = 0
        while True:
            try:
                bakName = 'Back;'+str(iBak)
                Background[0][iBak+3] = parmList[bakName]
                iBak += 1
            except KeyError:
                break
        iDb = 0
        while True:
            names = ['DebyeA;','DebyeR;','DebyeU;']
            try:
                for i,name in enumerate(names):
                    val = parmList[name+str(iDb)]
                    Background[1]['debyeTerms'][iDb][2*i] = val
                iDb += 1
            except KeyError:
                break
        iDb = 0
        while True:
            names = ['BkPkpos;','BkPkint;','BkPksig;','BkPkgam;']
            try:
                for i,name in enumerate(names):
                    val = parmList[name+str(iDb)]
                    Background[1]['peaksList'][iDb][2*i] = val
                iDb += 1
            except KeyError:
                break
                
    def BackgroundPrint(Background,sigDict):
        print ('Background coefficients for '+Background[0][0]+' function')
        ptfmt = "%12.5f"
        ptstr =  'value: '
        sigstr = 'esd  : '
        for i,back in enumerate(Background[0][3:]):
            ptstr += ptfmt % (back)
            if Background[0][1]:
                prm = 'Back;'+str(i)
                if prm in sigDict:
                    sigstr += ptfmt % (sigDict[prm])
                else:
                    sigstr += " "*12
            if len(ptstr) > 75:
                print (ptstr)
                if Background[0][1]: print (sigstr)
                ptstr =  'value: '
                sigstr = 'esd  : '
        if len(ptstr) > 8:
            print (ptstr)
            if Background[0][1]: print (sigstr)

        if Background[1]['nDebye']:
            parms = ['DebyeA;','DebyeR;','DebyeU;']
            print ('Debye diffuse scattering coefficients')
            ptfmt = "%12.5f"
            print (' term       DebyeA       esd        DebyeR       esd        DebyeU        esd')
            for term in range(Background[1]['nDebye']):
                line = ' term %d'%(term)
                for ip,name in enumerate(parms):
                    line += ptfmt%(Background[1]['debyeTerms'][term][2*ip])
                    if name+str(term) in sigDict:
                        line += ptfmt%(sigDict[name+str(term)])
                    else:
                        line += " "*12
                print (line)
        if Background[1]['nPeaks']:
            print ('Coefficients for Background Peaks')
            ptfmt = "%15.3f"
            for j,pl in enumerate(Background[1]['peaksList']):
                names =  'peak %3d:'%(j+1)
                ptstr =  'values  :'
                sigstr = 'esds    :'
                for i,lbl in enumerate(['BkPkpos','BkPkint','BkPksig','BkPkgam']):
                    val = pl[2*i]
                    prm = lbl+";"+str(j)
                    names += '%15s'%(prm)
                    ptstr += ptfmt%(val)
                    if prm in sigDict:
                        sigstr += ptfmt%(sigDict[prm])
                    else:
                        sigstr += " "*15
                print (names)
                print (ptstr)
                print (sigstr)
                            
    def SetInstParms(Inst):
        dataType = Inst['Type'][0]
        insVary = []
        insNames = []
        insVals = []
        for parm in Inst:
            insNames.append(parm)
            insVals.append(Inst[parm][1])
            if parm in ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','alpha-0','alpha-1'] and Inst[parm][2]:
                    insVary.append(parm)
        instDict = dict(zip(insNames,insVals))
#        instDict['X'] = max(instDict['X'],0.01)
#        instDict['Y'] = max(instDict['Y'],0.01)
        if 'SH/L' in instDict:
            instDict['SH/L'] = max(instDict['SH/L'],0.002)
        return dataType,instDict,insVary
        
    def GetInstParms(parmDict,Inst,varyList,Peaks):
        for name in Inst:
            Inst[name][1] = parmDict[name]
        iPeak = 0
        while True:
            try:
                sigName = 'sig'+str(iPeak)
                pos = parmDict['pos'+str(iPeak)]
                if sigName not in varyList:
                    if 'T' in Inst['Type'][0]:
                        dsp = G2lat.Pos2dsp(Inst,pos)
                        parmDict[sigName] = G2mth.getTOFsig(parmDict,dsp)
                    else:
                        parmDict[sigName] = G2mth.getCWsig(parmDict,pos)
                gamName = 'gam'+str(iPeak)
                if gamName not in varyList:
                    if 'T' in Inst['Type'][0]:
                        dsp = G2lat.Pos2dsp(Inst,pos)
                        parmDict[gamName] = G2mth.getTOFgamma(parmDict,dsp)
                    else:
                        parmDict[gamName] = G2mth.getCWgam(parmDict,pos)
                iPeak += 1
            except KeyError:
                break
        
    def InstPrint(Inst,sigDict):
        print ('Instrument Parameters:')
        ptfmt = "%12.6f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for parm in Inst:
            if parm in  ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','alpha-0','alpha-1']:
                ptlbls += "%s" % (parm.center(12))
                ptstr += ptfmt % (Inst[parm][1])
                if parm in sigDict:
                    sigstr += ptfmt % (sigDict[parm])
                else:
                    sigstr += 12*' '
        print (ptlbls)
        print (ptstr)
        print (sigstr)

    def SetPeaksParms(dataType,Peaks):
        peakNames = []
        peakVary = []
        peakVals = []
        if 'C' in dataType:
            names = ['pos','int','sig','gam']
        else:   #'T' and 'B'
            names = ['pos','int','alp','bet','sig','gam']
        for i,peak in enumerate(Peaks):
            for j,name in enumerate(names):
                peakVals.append(peak[2*j])
                parName = name+str(i)
                peakNames.append(parName)
                if peak[2*j+1]:
                    peakVary.append(parName)
        return dict(zip(peakNames,peakVals)),peakVary
                
    def GetPeaksParms(Inst,parmDict,Peaks,varyList):
        if 'C' in Inst['Type'][0]:
            names = ['pos','int','sig','gam']
        else:   #'T' & 'B'
            names = ['pos','int','alp','bet','sig','gam']
        for i,peak in enumerate(Peaks):
            pos = parmDict['pos'+str(i)]
            if 'difC' in Inst:
                dsp = pos/Inst['difC'][1]
            for j in range(len(names)):
                parName = names[j]+str(i)
                if parName in varyList:
                    peak[2*j] = parmDict[parName]
                elif 'alp' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j] = G2mth.getTOFalpha(parmDict,dsp)
                    else: #'B'
                        peak[2*j] = G2mth.getPinkalpha(parmDict,pos)
                elif 'bet' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j] = G2mth.getTOFbeta(parmDict,dsp)
                    else:   #'B'
                        peak[2*j] = G2mth.getPinkbeta(parmDict,pos)
                elif 'sig' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j] = G2mth.getTOFsig(parmDict,dsp)
                    else:   #'C' & 'B'
                        peak[2*j] = G2mth.getCWsig(parmDict,pos)
                elif 'gam' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j] = G2mth.getTOFgamma(parmDict,dsp)
                    else:   #'C' & 'B'
                        peak[2*j] = G2mth.getCWgam(parmDict,pos)
                        
    def PeaksPrint(dataType,parmDict,sigDict,varyList,ptsperFW):
        print ('Peak coefficients:')
        if 'C' in dataType:
            names = ['pos','int','sig','gam']
        else:   #'T' & 'B'
            names = ['pos','int','alp','bet','sig','gam']            
        head = 13*' '
        for name in names:
            if name in ['alp','bet']:
                head += name.center(8)+'esd'.center(8)
            else:
                head += name.center(10)+'esd'.center(10)
        head += 'bins'.center(8)
        print (head)
        if 'C' in dataType:
            ptfmt = {'pos':"%10.5f",'int':"%10.1f",'sig':"%10.3f",'gam':"%10.3f"}
        elif 'T' in dataType:
            ptfmt = {'pos':"%10.2f",'int':"%10.4f",'alp':"%8.3f",'bet':"%8.5f",'sig':"%10.3f",'gam':"%10.3f"}
        else: #'B'
            ptfmt = {'pos':"%10.5f",'int':"%10.1f",'alp':"%8.2f",'bet':"%8.4f",'sig':"%10.3f",'gam':"%10.3f"}
        for i,peak in enumerate(Peaks):
            ptstr =  ':'
            for j in range(len(names)):
                name = names[j]
                parName = name+str(i)
                ptstr += ptfmt[name] % (parmDict[parName])
                if parName in varyList:
                    ptstr += ptfmt[name] % (sigDict[parName])
                else:
                    if name in ['alp','bet']:
                        ptstr += 8*' '
                    else:
                        ptstr += 10*' '
            ptstr += '%9.2f'%(ptsperFW[i])
            print ('%s'%(('Peak'+str(i+1)).center(8)),ptstr)
                
    def devPeakProfile(values,xdata,ydata, weights,dataType,parmdict,varylist,bakType,dlg):
        parmdict.update(zip(varylist,values))
        return np.sqrt(weights)*getPeakProfileDerv(dataType,parmdict,xdata,varylist,bakType)
            
    def errPeakProfile(values,xdata,ydata,weights,dataType,parmdict,varylist,bakType,dlg):        
        parmdict.update(zip(varylist,values))
        M = np.sqrt(weights)*(getPeakProfile(dataType,parmdict,xdata,varylist,bakType)-ydata)
        Rwp = min(100.,np.sqrt(np.sum(M**2)/np.sum(weights*ydata**2))*100.)
        if dlg:
            dlg.Raise()
            GoOn = dlg.Update(Rwp,newmsg='%s%8.3f%s'%('Peak fit Rwp =',Rwp,'%'))[0]
            if not GoOn:
                return -M           #abort!!
        return M
        
    if controls:
        Ftol = controls['min dM/M']
    else:
        Ftol = 0.0001
    if oneCycle:
        Ftol = 1.0
    x,y,w,yc,yb,yd = data   #these are numpy arrays - remove masks!
    if fixback is None:
        fixback = np.zeros_like(y)
    yc *= 0.                            #set calcd ones to zero
    yb *= 0.
    yd *= 0.
    cw = x[1:]-x[:-1]
    xBeg = np.searchsorted(x,Limits[0])
    xFin = np.searchsorted(x,Limits[1])+1
    bakType,bakDict,bakVary = SetBackgroundParms(Background)
    dataType,insDict,insVary = SetInstParms(Inst)
    peakDict,peakVary = SetPeaksParms(Inst['Type'][0],Peaks)
    parmDict = {}
    parmDict.update(bakDict)
    parmDict.update(insDict)
    parmDict.update(peakDict)
    parmDict['Pdabc'] = []      #dummy Pdabc
    parmDict.update(Inst2)      #put in real one if there
    if prevVaryList:
        varyList = prevVaryList[:]
    else:
        varyList = bakVary+insVary+peakVary
    fullvaryList = varyList[:]
    while True:
        begin = time.time()
        values =  np.array(Dict2Values(parmDict, varyList))
        Rvals = {}
        badVary = []
        result = so.leastsq(errPeakProfile,values,Dfun=devPeakProfile,full_output=True,ftol=Ftol,col_deriv=True,
               args=(x[xBeg:xFin],(y+fixback)[xBeg:xFin],wtFactor*w[xBeg:xFin],dataType,parmDict,varyList,bakType,dlg))
        ncyc = int(result[2]['nfev']/2)
        runtime = time.time()-begin    
        chisq = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
        Rvals['Rwp'] = np.sqrt(chisq/np.sum(wtFactor*w[xBeg:xFin]*(y+fixback)[xBeg:xFin]**2))*100.      #to %
        Rvals['GOF'] = chisq/(xFin-xBeg-len(varyList))       #reduced chi^2
        G2fil.G2Print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(result[2]['nfev'],xFin-xBeg,len(varyList)))
        if ncyc:
            G2fil.G2Print ('fitpeak time = %8.3fs, %8.3fs/cycle'%(runtime,runtime/ncyc))
        G2fil.G2Print ('Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],chisq,Rvals['GOF']))
        sig = [0]*len(varyList)
        if len(varyList) == 0: break  # if nothing was refined
        try:
            sig = np.sqrt(np.diag(result[1])*Rvals['GOF'])
            if np.any(np.isnan(sig)):
                G2fil.G2Print ('*** Least squares aborted - some invalid esds possible ***')
            break                   #refinement succeeded - finish up!
        except ValueError:          #result[1] is None on singular matrix
            G2fil.G2Print ('**** Refinement failed - singular matrix ****')
            Ipvt = result[2]['ipvt']
            for i,ipvt in enumerate(Ipvt):
                if not np.sum(result[2]['fjac'],axis=1)[i]:
                    G2fil.G2Print ('Removing parameter: '+varyList[ipvt-1])
                    badVary.append(varyList[ipvt-1])
                    del(varyList[ipvt-1])
                    break
            else: # nothing removed
                break
    if dlg: dlg.Destroy()
    sigDict = dict(zip(varyList,sig))
    yb[xBeg:xFin] = getBackground('',parmDict,bakType,dataType,x[xBeg:xFin])[0]-fixback[xBeg:xFin]
    yc[xBeg:xFin] = getPeakProfile(dataType,parmDict,x[xBeg:xFin],varyList,bakType)-fixback[xBeg:xFin]
    yd[xBeg:xFin] = y[xBeg:xFin]-yc[xBeg:xFin]
    GetBackgroundParms(parmDict,Background)
    if bakVary: BackgroundPrint(Background,sigDict)
    GetInstParms(parmDict,Inst,varyList,Peaks)
    if insVary: InstPrint(Inst,sigDict)
    GetPeaksParms(Inst,parmDict,Peaks,varyList)
    binsperFWHM = []
    for peak in Peaks:
        FWHM = getFWHM(peak[0],Inst)
        try:
            binsperFWHM.append(FWHM/cw[x.searchsorted(peak[0])])
        except IndexError:
            binsperFWHM.append(0.)
    if peakVary: PeaksPrint(dataType,parmDict,sigDict,varyList,binsperFWHM)
    if len(binsperFWHM):
        if min(binsperFWHM) < 1.:
            G2fil.G2Print ('*** Warning: calculated peak widths are too narrow to refine profile coefficients ***')
            if 'T' in Inst['Type'][0]:
                G2fil.G2Print (' Manually increase sig-0, 1, or 2 in Instrument Parameters')
            else:
                G2fil.G2Print (' Manually increase W in Instrument Parameters')
        elif min(binsperFWHM) < 4.:
            G2fil.G2Print ('*** Warning: data binning yields too few data points across peak FWHM for reliable Rietveld refinement ***')
            G2fil.G2Print ('*** recommended is 6-10; you have %.2f ***'%(min(binsperFWHM)))
    return sigDict,result,sig,Rvals,varyList,parmDict,fullvaryList,badVary
    
def calcIncident(Iparm,xdata):
    'needs a doc string'

    def IfunAdv(Iparm,xdata):
        Itype = Iparm['Itype']
        Icoef = Iparm['Icoeff']
        DYI = np.ones((12,xdata.shape[0]))
        YI = np.ones_like(xdata)*Icoef[0]
        
        x = xdata/1000.                 #expressions are in ms
        if Itype == 'Exponential':
            for i in [1,3,5,7,9]:
                Eterm = np.exp(-Icoef[i+1]*x**((i+1)/2))
                YI += Icoef[i]*Eterm
                DYI[i] *= Eterm
                DYI[i+1] *= -Icoef[i]*Eterm*x**((i+1)/2)            
        elif 'Maxwell'in Itype:
            Eterm = np.exp(-Icoef[2]/x**2)
            DYI[1] = Eterm/x**5
            DYI[2] = -Icoef[1]*DYI[1]/x**2
            YI += (Icoef[1]*Eterm/x**5)
            if 'Exponential' in Itype:
                for i in range(3,11,2):
                    Eterm = np.exp(-Icoef[i+1]*x**((i+1)/2))
                    YI += Icoef[i]*Eterm
                    DYI[i] *= Eterm
                    DYI[i+1] *= -Icoef[i]*Eterm*x**((i+1)/2)
            else:   #Chebyschev
                T = (2./x)-1.
                Ccof = np.ones((12,xdata.shape[0]))
                Ccof[1] = T
                for i in range(2,12):
                    Ccof[i] = 2*T*Ccof[i-1]-Ccof[i-2]
                for i in range(1,10):
                    YI += Ccof[i]*Icoef[i+2]
                    DYI[i+2] =Ccof[i]
        return YI,DYI
        
    Iesd = np.array(Iparm['Iesd'])
    Icovar = Iparm['Icovar']
    YI,DYI = IfunAdv(Iparm,xdata)
    YI = np.where(YI>0,YI,1.)
    WYI = np.zeros_like(xdata)
    vcov = np.zeros((12,12))
    k = 0
    for i in range(12):
        for j in range(i,12):
            vcov[i][j] = Icovar[k]*Iesd[i]*Iesd[j]
            vcov[j][i] = Icovar[k]*Iesd[i]*Iesd[j]
            k += 1
    M = np.inner(vcov,DYI.T)
    WYI = np.sum(M*DYI,axis=0)
    WYI = np.where(WYI>0.,WYI,0.)
    return YI,WYI

################################################################################
#### RMCutilities
################################################################################
   
def MakeInst(PWDdata,Name,Size,Mustrain,useSamBrd):
    inst = PWDdata['Instrument Parameters'][0]
    Xsb = 0.
    Ysb = 0.
    if 'T' in inst['Type'][1]:
        difC = inst['difC'][1]
        if useSamBrd[0]:
            if 'ellipsoidal' not in Size[0]:    #take the isotropic term only
                Xsb = 1.e-4*difC/Size[1][0]
        if useSamBrd[1]:
            if 'generalized' not in Mustrain[0]:    #take the isotropic term only
                Ysb = 1.e-6*difC*Mustrain[1][0]
        prms = ['Bank',
                'difC','difA','Zero','2-theta',
                'alpha','beta-0','beta-1',
                'sig-0','sig-1','sig-2',
                'Z','X','Y']
        fname = Name+'.inst'
        fl = open(fname,'w')
        fl.write('1\n')
        fl.write('%d\n'%int(inst[prms[0]][1]))
        fl.write('%19.11f%19.11f%19.11f%19.11f\n'%(inst[prms[1]][1],inst[prms[2]][1],inst[prms[3]][1],inst[prms[4]][1]))
        fl.write('%12.6e%14.6e%14.6e\n'%(inst[prms[5]][1],inst[prms[6]][1],inst[prms[7]][1]))
        fl.write('%12.6e%14.6e%14.6e\n'%(inst[prms[8]][1],inst[prms[9]][1],inst[prms[10]][1]))    
        fl.write('%12.6e%14.6e%14.6e%14.6e%14.6e\n'%(inst[prms[11]][1],inst[prms[12]][1]+Ysb,inst[prms[13]][1]+Xsb,0.0,0.0))
        fl.close()
    else:
        if useSamBrd[0]:
            wave = G2mth.getWave(inst)
            if 'ellipsoidal' not in Size[0]:    #take the isotropic term only
                Xsb = 1.8*wave/(np.pi*Size[1][0])
        if useSamBrd[1]:
            if 'generalized' not in Mustrain[0]:    #take the isotropic term only
                Ysb = 0.0180*Mustrain[1][0]/np.pi
        prms = ['Bank',
                'Lam','Zero','Polariz.',
                'U','V','W',
                'X','Y']
        fname = Name+'.inst'
        fl = open(fname,'w')
        fl.write('1\n')
        fl.write('%d\n'%int(inst[prms[0]][1]))
        fl.write('%10.5f%10.5f%10.4f%10d\n'%(inst[prms[1]][1],-100.*inst[prms[2]][1],inst[prms[3]][1],0))
        fl.write('%10.3f%10.3f%10.3f\n'%(inst[prms[4]][1],inst[prms[5]][1],inst[prms[6]][1]))
        fl.write('%10.3f%10.3f%10.3f\n'%(inst[prms[7]][1]+Xsb,inst[prms[8]][1]+Ysb,0.0))    
        fl.write('%10.3f%10.3f%10.3f\n'%(0.0,0.0,0.0))
        fl.close()
    return fname
    
def MakeBack(PWDdata,Name):
    Back = PWDdata['Background'][0]
    inst = PWDdata['Instrument Parameters'][0]
    if 'chebyschev-1' != Back[0]:
        return None
    Nback = Back[2]
    BackVals = Back[3:]
    fname = Name+'.back'
    fl = open(fname,'w')
    fl.write('%10d\n'%Nback)
    for val in BackVals:
        if 'T' in inst['Type'][1]:
            fl.write('%12.6g\n'%(float(val)))
        else:
            fl.write('%12.6g\n'%val)
    fl.close()
    return fname

def findDup(Atoms):
    Dup = []
    Fracs = []
    for iat1,at1 in enumerate(Atoms):
        if any([at1[0] in dup for dup in Dup]):
            continue
        else:
            Dup.append([at1[0],])
            Fracs.append([at1[6],])
        for iat2,at2 in enumerate(Atoms[(iat1+1):]):
            if np.sum((np.array(at1[3:6])-np.array(at2[3:6]))**2) < 0.00001:
                Dup[-1] += [at2[0],]
                Fracs[-1]+= [at2[6],]
    return Dup,Fracs

def MakeRMC6f(PWDdata,Name,Phase,RMCPdict):    
    
    Meta = RMCPdict['metadata']
    Atseq = RMCPdict['atSeq']
    Supercell =  RMCPdict['SuperCell']
    generalData = Phase['General']
    Dups,Fracs = findDup(Phase['Atoms'])
    Sfracs = [np.cumsum(fracs) for fracs in Fracs]
    ifSfracs = any([np.any(sfracs-1.) for sfracs in Sfracs])
    Sample = PWDdata['Sample Parameters']
    Meta['temperature'] = Sample['Temperature']
    Meta['pressure'] = Sample['Pressure']
    Cell = generalData['Cell'][1:7]
    Trans = np.eye(3)*np.array(Supercell)
    newPhase = copy.deepcopy(Phase)
    newPhase['General']['SGData'] = G2spc.SpcGroup('P 1')[1]
    newPhase['General']['Cell'][1:] = G2lat.TransformCell(Cell,Trans)
    GB = G2lat.cell2Gmat( newPhase['General']['Cell'][1:7])[0]
    RMCPdict['Rmax'] = np.min(np.sqrt(np.array([1./G2lat.calc_rDsq2(H,GB) for H in [[1,0,0],[0,1,0],[0,0,1]]])))/2.
    newPhase,Atcodes = G2lat.TransformPhase(Phase,newPhase,Trans,np.zeros(3),np.zeros(3),ifMag=False,Force=True)
    Natm = np.core.defchararray.count(np.array(Atcodes),'+')    #no. atoms in original unit cell
    Natm = np.count_nonzero(Natm-1)
    Atoms = newPhase['Atoms']
    reset = False
    
    if ifSfracs:
        Natm = np.core.defchararray.count(np.array(Atcodes),'+')    #no. atoms in original unit cell
        Natm = np.count_nonzero(Natm-1)
        Satoms = []
        for i in range(len(Atoms)//Natm):
            ind = i*Natm
            Satoms.append(G2mth.sortArray(G2mth.sortArray(G2mth.sortArray(Atoms[ind:ind+Natm],5),4),3))
        Natoms = []
        for satoms in Satoms:
            for idup,dup in enumerate(Dups):
                ldup = len(dup)
                natm = len(satoms)
                i = 0
                while i < natm:
                    if satoms[i][0] in dup:
                        atoms = satoms[i:i+ldup]
                        try:
                            atom = atoms[np.searchsorted(Sfracs[idup],rand.random())]
                            Natoms.append(atom)
                        except IndexError:      #what about vacancies?
                            if 'Va' not in Atseq:
                                reset = True
                                Atseq.append('Va')
                                RMCPdict['aTypes']['Va'] = 0.0
                            atom = atoms[0]
                            atom[1] = 'Va'
                            Natoms.append(atom)
                        i += ldup
                    else:
                       i += 1
    else:
        Natoms = Atoms
    
    NAtype = np.zeros(len(Atseq))
    for atom in Natoms:
        NAtype[Atseq.index(atom[1])] += 1
    NAstr = ['%6d'%i for i in NAtype]
    Cell = newPhase['General']['Cell'][1:7]
    if os.path.exists(Name+'.his6f'):
        os.remove(Name+'.his6f')
    if os.path.exists(Name+'.neigh'):
        os.remove(Name+'.neigh')
    fname = Name+'.rmc6f'
    fl = open(fname,'w')
    fl.write('(Version 6f format configuration file)\n')
    for item in Meta:
        fl.write('%-20s%s\n'%('Metadata '+item+':',Meta[item]))
    fl.write('Atom types present:                 %s\n'%'    '.join(Atseq))
    fl.write('Number of each atom type:       %s\n'%''.join(NAstr))
    fl.write('Number of atoms:                %d\n'%len(Natoms))
    fl.write('%-35s%4d%4d%4d\n'%('Supercell dimensions:',Supercell[0],Supercell[1],Supercell[2]))
    fl.write('Cell (Ang/deg): %12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n'%(
            Cell[0],Cell[1],Cell[2],Cell[3],Cell[4],Cell[5]))
    A,B = G2lat.cell2AB(Cell,True)
    fl.write('Lattice vectors (Ang):\n')   
    for i in [0,1,2]:
        fl.write('%12.6f%12.6f%12.6f\n'%(A[i,0],A[i,1],A[i,2]))
    fl.write('Atoms (fractional coordinates):\n')
    nat = 0
    for atm in Atseq:
        for iat,atom in enumerate(Natoms):
            if atom[1] == atm:
                nat += 1
                atcode = Atcodes[iat].split(':')
                cell = [0,0,0]
                if '+' in atcode[1]:
                    cell = eval(atcode[1].split('+')[1])
                fl.write('%6d%4s  [%s]%19.15f%19.15f%19.15f%6d%4d%4d%4d\n'%(       
                        nat,atom[1].strip(),atcode[0],atom[3],atom[4],atom[5],(iat)%Natm+1,cell[0],cell[1],cell[2]))
    fl.close()
    return fname,reset

def MakeBragg(PWDdata,Name,Phase):
    generalData = Phase['General']
    Vol = generalData['Cell'][7]
    Data = PWDdata['Data']
    Inst = PWDdata['Instrument Parameters'][0]
    Bank = int(Inst['Bank'][1])
    Sample = PWDdata['Sample Parameters']
    Scale = Sample['Scale'][0]
    if 'X' in Inst['Type'][1]:
        Scale *= 2.
    Limits = PWDdata['Limits'][1]
    Ibeg = np.searchsorted(Data[0],Limits[0])
    Ifin = np.searchsorted(Data[0],Limits[1])+1
    fname = Name+'.bragg'
    fl = open(fname,'w')
    fl.write('%12d%6d%15.7f%15.4f\n'%(Ifin-Ibeg-2,Bank,Scale,Vol))
    if 'T' in Inst['Type'][0]:
        fl.write('%12s%12s\n'%('   TOF,ms','  I(obs)'))
        for i in range(Ibeg,Ifin-1):
            fl.write('%12.8f%12.6f\n'%(Data[0][i]/1000.,Data[1][i]))
    else:
        fl.write('%12s%12s\n'%('   2-theta, deg','  I(obs)'))
        for i in range(Ibeg,Ifin-1):
            fl.write('%11.6f%15.2f\n'%(Data[0][i],Data[1][i]))        
    fl.close()
    return fname

def MakeRMCPdat(PWDdata,Name,Phase,RMCPdict):
    Meta = RMCPdict['metadata']
    Times = RMCPdict['runTimes']
    Atseq = RMCPdict['atSeq']
    Atypes = RMCPdict['aTypes']
    atPairs = RMCPdict['Pairs']
    Files = RMCPdict['files']
    BraggWt = RMCPdict['histogram'][1]
    inst = PWDdata['Instrument Parameters'][0]
    try:
        refList = PWDdata['Reflection Lists'][Name]['RefList']
    except KeyError:
        return 'Error - missing reflection list; you must do Refine first'
    dMin = refList[-1][4]
    gsasType = 'xray2'
    if 'T' in inst['Type'][1]:
        gsasType = 'gsas3'
    elif 'X' in inst['Type'][1]:
        XFF = G2elem.GetFFtable(Atseq)
        Xfl = open(Name+'.xray','w')
        for atm in Atseq:
            fa = XFF[atm]['fa']
            fb = XFF[atm]['fb']
            fc = XFF[atm]['fc']
            Xfl.write('%2s  %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n'%(
                    atm.upper(),fa[0],fb[0],fa[1],fb[1],fa[2],fb[2],fa[3],fb[3],fc))
        Xfl.close()
    lenA = len(Atseq)
    Pairs = []
    for pair in [[' %s-%s'%(Atseq[i],Atseq[j]) for j in range(i,lenA)] for i in range(lenA)]:
        Pairs += pair
    pairMin = [atPairs[pair]for pair in Pairs if pair in atPairs]
    maxMoves = [Atypes[atm] for atm in Atseq if atm in Atypes]
    fname = Name+'.dat'
    fl = open(fname,'w')
    fl.write(' %% Hand edit the following as needed\n')
    fl.write('TITLE :: '+Name+'\n')
    fl.write('MATERIAL :: '+Meta['material']+'\n')
    fl.write('PHASE :: '+Meta['phase']+'\n')
    fl.write('TEMPERATURE :: '+str(Meta['temperature'])+'\n')
    fl.write('INVESTIGATOR :: '+Meta['owner']+'\n')
    minHD = ' '.join(['%6.3f'%dist[0] for dist in pairMin])
    minD = ' '.join(['%6.3f'%dist[1] for dist in pairMin])
    maxD = ' '.join(['%6.3f'%dist[2] for dist in pairMin])
    fl.write('MINIMUM_DISTANCES ::   %s  Angstrom\n'%minHD)
    maxMv = ' '.join(['%6.3f'%mov for mov in maxMoves])
    fl.write('MAXIMUM_MOVES ::   %s Angstrom\n'%maxMv)
    fl.write('R_SPACING ::  0.0200 Angstrom\n')
    fl.write('PRINT_PERIOD :: 100\n')
    fl.write('TIME_LIMIT ::     %.2f MINUTES\n'%Times[0])
    fl.write('SAVE_PERIOD ::    %.2f MINUTES\n'%Times[1])
    fl.write('\n')
    fl.write('ATOMS :: '+' '.join(Atseq)+'\n')
    fl.write('\n')
    fl.write('FLAGS ::\n')
    fl.write('  > NO_MOVEOUT\n')
    fl.write('  > NO_SAVE_CONFIGURATIONS\n')
    fl.write('  > NO_RESOLUTION_CONVOLUTION\n')
    fl.write('\n')
    fl.write('INPUT_CONFIGURATION_FORMAT ::  rmc6f\n')
    fl.write('SAVE_CONFIGURATION_FORMAT  ::  rmc6f\n')
    fl.write('IGNORE_HISTORY_FILE ::\n')
    fl.write('\n')
    fl.write('DISTANCE_WINDOW ::\n')
    fl.write('  > MNDIST :: %s\n'%minD)
    fl.write('  > MXDIST :: %s\n'%maxD)
    if len(RMCPdict['Potentials']['Stretch']) or len(RMCPdict['Potentials']['Stretch']):
        fl.write('\n')
        fl.write('POTENTIALS ::\n')
        fl.write('  > TEMPERATURE :: %.1f K\n'%RMCPdict['Potentials']['Pot. Temp.'])
        fl.write('  > PLOT :: pixels=400, colour=red, zangle=90, zrotation=45 deg\n')
        if len(RMCPdict['Potentials']['Stretch']):
            fl.write('  > STRETCH_SEARCH :: %.1f%%\n'%RMCPdict['Potentials']['Stretch search'])
            for bond in RMCPdict['Potentials']['Stretch']:
                fl.write('  > STRETCH :: %s %s %.2f eV %.2f Ang\n'%(bond[0],bond[1],bond[3],bond[2]))        
        if len(RMCPdict['Potentials']['Angles']):
            fl.write('  > ANGLE_SEARCH :: %.1f%%\n'%RMCPdict['Potentials']['Angle search'])
            for angle in RMCPdict['Potentials']['Angles']:
                fl.write('  > ANGLE :: %s %s %s %.2f eV %.2f deg %.2f %.2f Ang\n'%
                    (angle[1],angle[0],angle[2],angle[6],angle[3],angle[4],angle[5]))
    if RMCPdict['useBVS']:
        fl.write('BVS ::\n')
        fl.write('  > ATOM :: '+' '.join(Atseq)+'\n')
        fl.write('  > WEIGHTS :: %s\n'%' '.join(['%6.3f'%RMCPdict['BVS'][bvs][2] for bvs in RMCPdict['BVS']]))
        oxid = []
        for val in RMCPdict['Oxid']:
            if len(val) == 3:
                oxid.append(val[0][1:])
            else:
                oxid.append(val[0][2:])
        fl.write('  > OXID :: %s\n'%' '.join(oxid))
        fl.write('  > RIJ :: %s\n'%' '.join(['%6.3f'%RMCPdict['BVS'][bvs][0] for bvs in RMCPdict['BVS']]))
        fl.write('  > BVAL :: %s\n'%' '.join(['%6.3f'%RMCPdict['BVS'][bvs][1] for bvs in RMCPdict['BVS']]))
        fl.write('  > CUTOFF :: %s\n'%' '.join(['%6.3f'%RMCPdict['BVS'][bvs][2] for bvs in RMCPdict['BVS']]))        
        fl.write('  > SAVE :: 100000\n')
        fl.write('  > UPDATE :: 100000\n')
        if len(RMCPdict['Swap']):
            fl.write('\n')
            fl.write('SWAP_MULTI ::\n')
            for swap in RMCPdict['Swap']:
                try:
                    at1 = Atseq.index(swap[0])
                    at2 = Atseq.index(swap[1])
                except ValueError:
                    break
                fl.write('  > SWAP_ATOMS :: %d %d %.2f\n'%(at1,at2,swap[2]))
        
    if len(RMCPdict['FxCN']):
        fl.write('FIXED_COORDINATION_CONSTRAINTS ::  %d\n'%len(RMCPdict['FxCN']))        
        for ifx,fxcn in enumerate(RMCPdict['FxCN']):
            try:
                at1 = Atseq.index(fxcn[0])
                at2 = Atseq.index(fxcn[1])
            except ValueError:
                break
            fl.write('  > CSTR%d ::   %d %d %.2f %.2f %.2f %.2f %.6f\n'%(ifx+1,at1+1,at2+1,fxcn[2],fxcn[3],fxcn[4],fxcn[5],fxcn[6]))
    if len(RMCPdict['AveCN']):
        fl.write('AVERAGE_COORDINATION_CONSTRAINTS ::  %d\n'%len(RMCPdict['AveCN']))
        for iav,avcn in enumerate(RMCPdict['AveCN']):
            try:
                at1 = Atseq.index(avcn[0])
                at2 = Atseq.index(avcn[1])
            except ValueError:
                break
            fl.write('  > CAVSTR%d ::   %d %d %.2f %.2f %.2f %.6f\n'%(iav+1,at1+1,at2+1,avcn[2],avcn[3],avcn[4],avcn[5]))
    for File in Files:
        if Files[File][0] and Files[File][0] != 'Select':
            if 'Xray' in File and 'F(Q)' in File:
                fqdata = open(Files[File][0],'r')
                lines = int(fqdata.readline()[:-1])
            fl.write('\n')
            fl.write('%s ::\n'%File.split(';')[0].upper().replace(' ','_'))
            fl.write('  > FILENAME :: %s\n'%Files[File][0])
            fl.write('  > DATA_TYPE :: %s\n'%Files[File][2])
            fl.write('  > FIT_TYPE :: %s\n'%Files[File][2])
            if 'Xray' not in File:
                fl.write('  > START_POINT :: 1\n')
                fl.write('  > END_POINT :: 3000\n')
                fl.write('  > WEIGHT :: %.4f\n'%Files[File][1])
            fl.write('  > CONSTANT_OFFSET 0.000\n')
            fl.write('  > NO_FITTED_OFFSET\n')
            if RMCPdict['FitScale']:
                fl.write('  > FITTED_SCALE\n')
            else:
                fl.write('  > NO_FITTED_SCALE\n')
            if Files[File][3] !='RMC':
                fl.write('  > %s\n'%Files[File][3])
            if 'reciprocal' in File:
                fl.write('  > CONVOLVE ::\n')
                if 'Xray' in File:
                    fl.write('  > RECIPROCAL_SPACE_FIT :: 1 %d 1\n'%lines)
                    fl.write('  > RECIPROCAL_SPACE_PARAMETERS :: 1 %d %.4f\n'%(lines,Files[File][1]))
                    fl.write('  > REAL_SPACE_FIT :: 1 %d 1\n'%(3*lines//2))
                    fl.write('  > REAL_SPACE_PARAMETERS :: 1 %d %.4f\n'%(3*lines//2,1./Files[File][1]))
    fl.write('\n')
    fl.write('BRAGG ::\n')
    fl.write('  > BRAGG_SHAPE :: %s\n'%gsasType)
    fl.write('  > RECALCUATE\n')
    fl.write('  > DMIN :: %.2f\n'%(dMin-0.02))
    fl.write('  > WEIGHT :: %10.3f\n'%BraggWt)
    fl.write('\n')
    fl.write('END  ::\n')
    fl.close()
    return fname

# def FindBonds(Phase,RMCPdict):
#     generalData = Phase['General']
#     cx,ct,cs,cia = generalData['AtomPtrs']
#     atomData = Phase['Atoms']
#     Res = 'RMC'
#     if 'macro' in generalData['Type']:
#         Res = atomData[0][ct-3]
#     AtDict = {atom[ct-1]:atom[ct] for atom in atomData}
#     Pairs = RMCPdict['Pairs']   #dict!
#     BondList = []
#     notNames = []
#     for FrstName in AtDict:
#         nbrs = G2mth.FindAllNeighbors(Phase,FrstName,list(AtDict.keys()),notName=notNames,Short=True)[0]
#         Atyp1 = AtDict[FrstName]
#         if 'Va' in Atyp1:
#             continue
#         for nbr in nbrs:
#             Atyp2 = AtDict[nbr[0]]
#             if 'Va' in Atyp2:
#                 continue
#             try:
#                 bndData = Pairs[' %s-%s'%(Atyp1,Atyp2)][1:]
#             except KeyError:
#                 bndData = Pairs[' %s-%s'%(Atyp2,Atyp1)][1:]
#             if any(bndData):
#                 if bndData[0] <= nbr[1] <= bndData[1]:
#                     bondStr = str((FrstName,nbr[0])+tuple(bndData))+',\n'
#                     revbondStr = str((nbr[0],FrstName)+tuple(bndData))+',\n'
#                     if bondStr not in BondList and revbondStr not in BondList:
#                         BondList.append(bondStr)
#         notNames.append(FrstName)
#     return Res,BondList

# def FindAngles(Phase,RMCPdict):
#     generalData = Phase['General']
#     Cell = generalData['Cell'][1:7]
#     Amat = G2lat.cell2AB(Cell)[0]
#     cx,ct,cs,cia = generalData['AtomPtrs']
#     atomData = Phase['Atoms']
#     AtLookup = G2mth.FillAtomLookUp(atomData,cia+8)
#     AtDict = {atom[ct-1]:atom[ct] for atom in atomData}
#     Angles = RMCPdict['Angles']
#     AngDict = {'%s-%s-%s'%(angle[0],angle[1],angle[2]):angle[3:] for angle in Angles}
#     AngleList = []
#     for MidName in AtDict:
#         nbrs,nbrIds = G2mth.FindAllNeighbors(Phase,MidName,list(AtDict.keys()),Short=True)
#         if len(nbrs) < 2: #need 2 neighbors to make an angle
#             continue
#         Atyp2 = AtDict[MidName]
#         for i,nbr1 in enumerate(nbrs):
#             Atyp1 = AtDict[nbr1[0]]
#             for j,nbr3 in enumerate(nbrs[i+1:]):
#                 Atyp3 = AtDict[nbr3[0]]
#                 IdList = [nbrIds[1][i],nbrIds[0],nbrIds[1][i+j+1]]
#                 try:
#                     angData = AngDict['%s-%s-%s'%(Atyp1,Atyp2,Atyp3)]
#                 except KeyError:
#                     try:
#                         angData = AngDict['%s-%s-%s'%(Atyp3,Atyp2,Atyp1)]
#                     except KeyError:
#                         continue
#                 XYZ = np.array(G2mth.GetAtomItemsById(atomData,AtLookup,IdList,cx,numItems=3))
#                 calAngle = G2mth.getRestAngle(XYZ,Amat)
#                 if angData[0] <= calAngle <= angData[1]:
#                     angStr = str((MidName,nbr1[0],nbr3[0])+tuple(angData))+',\n'
#                     revangStr = str((MidName,nbr3[0],nbr1[0])+tuple(angData))+',\n'
#                     if angStr not in AngleList and revangStr not in AngleList:
#                         AngleList.append(angStr)
#     return AngleList

# def GetSqConvolution(XY,d):

#     n = XY.shape[1]
#     snew = np.zeros(n)
#     dq = np.zeros(n)
#     sold = XY[1]
#     q = XY[0]
#     dq[1:] = np.diff(q)
#     dq[0] = dq[1]
    
#     for j in range(n):
#         for i in range(n):
#             b = abs(q[i]-q[j])
#             t = q[i]+q[j]
#             if j == i:
#                 snew[j] += q[i]*sold[i]*(d-np.sin(t*d)/t)*dq[i]
#             else:
#                 snew[j] += q[i]*sold[i]*(np.sin(b*d)/b-np.sin(t*d)/t)*dq[i]
#         snew[j] /= np.pi*q[j]
    
#     snew[0] = snew[1]
#     return snew

# def GetMaxSphere(pdbName):
#     try:
#         pFil = open(pdbName,'r')
#     except FileNotFoundError:
#         return None
#     while True:
#         line = pFil.readline()
#         if 'Boundary' in line:
#             line = line.split()[3:]
#             G = np.array([float(item) for item in line])
#             G = np.reshape(G,(3,3))**2
#             G = nl.inv(G)
#             pFil.close()
#             break
#     dspaces = [0.5/np.sqrt(G2lat.calc_rDsq2(H,G)) for H in np.eye(3)]
#     return min(dspaces)
    
def MakefullrmcRun(pName,Phase,RMCPdict):
    BondList = {}
    for k in RMCPdict['Pairs']:
        if RMCPdict['Pairs'][k][1]+RMCPdict['Pairs'][k][2]>0:
            BondList[k] = (RMCPdict['Pairs'][k][1],RMCPdict['Pairs'][k][2])
    AngleList = []
    for angle in RMCPdict['Angles']:
        if angle[3] == angle[4] or angle[5] >= angle[6] or angle[6] <= 0:
            continue
        for i in (0,1,2):
            angle[i] = angle[i].strip()
        AngleList.append(angle)
    rmin = RMCPdict['min Contact']
    cell = Phase['General']['Cell'][1:7]
    SymOpList = G2spc.AllOps(Phase['General']['SGData'])[0]
    cx,ct,cs,cia = Phase['General']['AtomPtrs']
    atomsList = []
    for atom in Phase['Atoms']:
        el = ''.join([i for i in atom[ct] if i.isalpha()])
        atomsList.append([el] + atom[cx:cx+4])
    rname = pName+'-fullrmc.py'
    restart = '%s_restart.pdb'%pName
    Files = RMCPdict['files']
    rundata = ''
    rundata += '#### fullrmc %s file; edit by hand if you so choose #####\n'%rname
    rundata += '# created in '+__file__+" v"+filversion.split()[1]
    rundata += dt.datetime.strftime(dt.datetime.now()," at %Y-%m-%dT%H:%M\n")
    rundata += '''
# fullrmc imports (all that are potentially useful)
import os,glob
import time
#import matplotlib.pyplot as plt
import numpy as np
from fullrmc.Core import Collection
from fullrmc.Engine import Engine
import fullrmc.Constraints.PairDistributionConstraints as fPDF
from fullrmc.Constraints.StructureFactorConstraints import ReducedStructureFactorConstraint, StructureFactorConstraint
from fullrmc.Constraints.DistanceConstraints import DistanceConstraint
from fullrmc.Constraints.BondConstraints import BondConstraint
from fullrmc.Constraints.AngleConstraints import BondsAngleConstraint
from fullrmc.Constraints.DihedralAngleConstraints import DihedralAngleConstraint
from fullrmc.Generators.Swaps import SwapPositionsGenerator

### When True, erases an existing enging to provide a fresh start
FRESH_START = {}
time0 = time.time()
'''.format(RMCPdict['ReStart'][0])
    
    rundata += '# setup structure\n'
    rundata += 'cell = ' + str(cell) + '\n'
    rundata += "SymOpList = "+str([i.lower() for i in SymOpList]) + '\n'
    rundata += 'atomList = ' + str(atomsList).replace('],','],\n  ') + '\n'
    rundata += 'supercell = ' + str(RMCPdict['SuperCell']) + '\n'

    rundata += '\n# initialize engine\n'
    rundata += 'engineFileName = "%s.rmc"\n'%pName

    rundata += '''\n# check Engine exists if so (and not FRESH_START) load it
# otherwise build it 
ENGINE = Engine(path=None)
if not ENGINE.is_engine(engineFileName) or FRESH_START:
    ## create structure
    ENGINE = Engine(path=engineFileName, freshStart=True)
    ENGINE.build_crystal_set_pdb(symOps     = SymOpList,
                                 atoms      = atomList,
                                 unitcellBC = cell,
                                 supercell  = supercell)
    rmax = min( [ENGINE.boundaryConditions.get_a(), ENGINE.boundaryConditions.get_b(), ENGINE.boundaryConditions.get_c()] ) /2.
'''    
    import atmdata
    rundata += '# conversion factors (may be needed)\n'
    rundata += '    sumCiBi2 = 0.\n'
    for elem in Phase['General']['AtomTypes']:
        rundata += '    Ci = ENGINE.numberOfAtomsPerElement["{}"]/len(ENGINE.allElements)\n'.format(elem)
        rundata += '    sumCiBi2 += (Ci*{})**2\n'.format(atmdata.AtmBlens[elem+'_']['SL'][0])
    rundata += '    rho0 = len(ENGINE.allNames)/ENGINE.volume\n'
    # settings that require a new Engine
    for File in Files:
        filDat = RMCPdict['files'][File]
        if not os.path.exists(filDat[0]): continue
        sfwt = 'neutronCohb'
        if 'Xray' in File:
            sfwt = 'atomicNumber'
        if 'G(r)' in File:
            rundata += '    GR = np.loadtxt("%s").T\n'%filDat[0]
            if filDat[3] == 0:
                rundata += '''    # read and xform G(r) as defined in RMCProfile
    # see eq. 44 in Keen, J. Appl. Cryst. (2001) 34 172-177\n'''
                rundata += '    GR[1] *= 4 * np.pi * GR[0] * rho0 / sumCiBi2\n'
                rundata += '    GofR = fPDF.PairDistributionConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            elif filDat[3] == 1:
                rundata += '    # This is G(r) as defined in PDFFIT\n'
                rundata += '    GofR = fPDF.PairDistributionConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            elif filDat[3] == 2:
                rundata += '    # This is g(r)\n'
                rundata += '    GofR = fPDF.PairCorrelationConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            else:
                raise ValueError('Invalid G(r) type: '+str(filDat[3]))
            rundata += '    ENGINE.add_constraints([GofR])\n'
            rundata += '    GofR.set_limits((None, rmax))\n'
        elif '(Q)' in File:
            rundata += '    SOQ = np.loadtxt("%s").T\n'%filDat[0]
            if filDat[3] == 0:
                rundata += '    # Read & xform F(Q) as defined in RMCProfile to S(Q)-1\n'
                rundata += '    SOQ[1] *= 1 / sumCiBi2\n'
            elif filDat[3] == 1:
                rundata += '    # This is S(Q) as defined in PDFFIT\n'
                rundata += '    SOQ[1] -= 1\n'
            if filDat[4]:
                rundata += '    SOQ[1] = Collection.sinc_convolution(q=SOQ[0],sq=SOQ[1],rmax=rmax)\n'
            rundata += '    SofQ = ReducedStructureFactorConstraint(experimentalData=SOQ.T, weighting="%s")\n'%sfwt
            rundata += '    ENGINE.add_constraints([SofQ])\n'
        else:
            print('What is this?')
    rundata += '    ENGINE.add_constraints(DistanceConstraint(defaultLowerDistance={}))\n'.format(RMCPdict['min Contact'])
    if BondList:
        rundata += '''    B_CONSTRAINT   = BondConstraint()
    ENGINE.add_constraints(B_CONSTRAINT)
    B_CONSTRAINT.create_supercell_bonds(bondsDefinition=[
'''
        for pair in BondList:
            e1,e2 = pair.split('-')
            rundata += '            ("element","{}","{}",{},{}),\n'.format(
                                        e1.strip(),e2.strip(),*BondList[pair])
        rundata += '             ])\n'
    if AngleList:
        rundata += '''    A_CONSTRAINT   = BondsAngleConstraint()
    ENGINE.add_constraints(A_CONSTRAINT)
    A_CONSTRAINT.create_supercell_angles(anglesDefinition=[
'''
        for item in AngleList:
            rundata += ('            '+
               '("element","{1}","{0}","{2}",{5},{6},{5},{6},{3},{4}),\n'.format(*item))
        rundata += '             ])\n'
    rundata += '    for f in glob.glob("{}_*.log"): os.remove(f)\n'.format(pName)
    rundata += '''
    ENGINE.save()
else:
    ENGINE = ENGINE.load(path=engineFileName)
'''
    rundata += 'ENGINE.set_log_file("{}")\n'.format(pName)
    if RMCPdict['Swaps']:
        rundata += '\n#set up for site swaps\n'
        rundata += 'aN = ENGINE.allNames\n'
        rundata += 'SwapGen = {}\n'
        for swap in RMCPdict['Swaps']:
            rundata += 'SwapA = [[idx] for idx in range(len(aN)) if aN[idx]=="%s"]\n'%swap[0]
            rundata += 'SwapB = [[idx] for idx in range(len(aN)) if aN[idx]=="%s"]\n'%swap[1]
            rundata += 'SwapGen["%s-%s"] = [SwapPositionsGenerator(swapList=SwapA),SwapPositionsGenerator(swapList=SwapB),%.2f]\n'%(swap[0],swap[1],swap[2])
        rundata += '    for swaps in SwapGen:\n'
        rundata += '        AB = swaps.split("-")\n'
        rundata += '        ENGINE.set_groups_as_atoms()\n'
        rundata += '        for g in ENGINE.groups:\n'
        rundata += '            if aN[g.indexes[0]]==AB[0]:\n'
        rundata += '                g.set_move_generator(SwapGen[swaps][0])\n'
        rundata += '            elif aN[g.indexes[0]]==AB[1]:\n'
        rundata += '                g.set_move_generator(SwapGen[swaps][1])\n'
        rundata += '            sProb = SwapGen[swaps][2]\n'
    rundata += '\n# set weights -- do this now so values can be changed without a restart\n'
    rundata += 'wtDict = {}\n'
    for File in Files:
        filDat = RMCPdict['files'][File]
        if not os.path.exists(filDat[0]): continue
        if 'Xray' in File:
            sfwt = 'atomicNumber'
        else:
            sfwt = 'neutronCohb'
        if 'G(r)' in File:
            typ = 'Pair'
        elif '(Q)' in File:
            typ = 'Struct'
        rundata += 'wtDict["{}-{}"] = {}\n'.format(typ,sfwt,filDat[1])
    rundata += 'for c in ENGINE.constraints:  # loop over predefined constraints\n'
    rundata += '    if type(c) is fPDF.PairDistributionConstraint:\n'
    rundata += '        c.set_variance_squared(1./wtDict["Pair-"+c.weighting])\n'
    rundata += '        c.set_limits((None,rmax))\n'
    if RMCPdict['FitScale']:
        rundata += '        c.set_adjust_scale_factor((10, 0.01, 100.))\n'
    rundata += '    elif type(c) is ReducedStructureFactorConstraint:\n'
    rundata += '        c.set_variance_squared(1./wtDict["Struct-"+c.weighting])\n'
    if RMCPdict['FitScale']:
        rundata += '        c.set_adjust_scale_factor((10, 0.01, 100.))\n'
    # torsions difficult to implement, must be internal to cell & named with
    # fullrmc atom names
    # if len(RMCPdict['Torsions']):         # Torsions currently commented out in GUI
    #     rundata += 'for c in ENGINE.constraints:  # look for Dihedral Angle Constraints\n'
    #     rundata += '    if type(c) is DihedralAngleConstraint:\n'
    #     rundata += '        c.set_variance_squared(%f)\n'%RMCPdict['Torsion Weight']
    #     rundata += '        c.create_angles_by_definition(anglesDefinition={"%s":[\n'%Res
    #     for torsion in RMCPdict['Torsions']:
    #         rundata += '    %s\n'%str(tuple(torsion))
    #     rundata += '        ]})\n'            
    rundata += '\n# setup runs for fullrmc\n'

    rundata += 'steps = 10000\n'
    rundata += 'for _ in range({}):\n'.format(RMCPdict['Cycles'])
    rundata += '    ENGINE.set_groups_as_atoms()\n'
    rundata += '    expected = ENGINE.generated+steps\n'
    
    rundata += '    ENGINE.run(restartPdb="%s",numberOfSteps=steps, saveFrequency=1000)\n'%restart
    rundata += '    if ENGINE.generated != expected: break # run was stopped\n'
    rundata += 'print("ENGINE run time %.2f s"%(time.time()-time0))\n'
    rfile = open(rname,'w')
    rfile.writelines(rundata)
    rfile.close()
    return rname
    
def GetRMCBonds(general,RMCPdict,Atoms,bondList):
    bondDist = []
    Cell = general['Cell'][1:7]
    Supercell =  RMCPdict['SuperCell']
    Trans = np.eye(3)*np.array(Supercell)
    Cell = G2lat.TransformCell(Cell,Trans)[:6]
    Amat,Bmat = G2lat.cell2AB(Cell)
    indices = (-1,0,1)
    Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    for bonds in bondList:
        Oxyz = np.array(Atoms[bonds[0]][1:])
        Txyz = np.array([Atoms[tgt-1][1:] for tgt in bonds[1]])        
        Dx = np.array([Txyz-Oxyz+unit for unit in Units])
        Dx = np.sqrt(np.sum(np.inner(Dx,Amat)**2,axis=2))
        for dx in Dx.T:
            bondDist.append(np.min(dx))
    return np.array(bondDist)
    
def GetRMCAngles(general,RMCPdict,Atoms,angleList):
    bondAngles = []
    Cell = general['Cell'][1:7]
    Supercell =  RMCPdict['SuperCell']
    Trans = np.eye(3)*np.array(Supercell)
    Cell = G2lat.TransformCell(Cell,Trans)[:6]
    Amat,Bmat = G2lat.cell2AB(Cell)
    indices = (-1,0,1)
    Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    for angle in angleList:
        Oxyz = np.array(Atoms[angle[0]][1:])
        TAxyz = np.array([Atoms[tgt-1][1:] for tgt in angle[1].T[0]])        
        TBxyz = np.array([Atoms[tgt-1][1:] for tgt in angle[1].T[1]])        
        DAxV = np.inner(np.array([TAxyz-Oxyz+unit for unit in Units]),Amat)
        DAx = np.sqrt(np.sum(DAxV**2,axis=2))
        DBxV = np.inner(np.array([TBxyz-Oxyz+unit for unit in Units]),Amat)
        DBx = np.sqrt(np.sum(DBxV**2,axis=2))
        iDAx = np.argmin(DAx,axis=0)
        iDBx = np.argmin(DBx,axis=0)
        for i,[iA,iB] in enumerate(zip(iDAx,iDBx)):
            DAv = DAxV[iA,i]/DAx[iA,i]
            DBv = DBxV[iB,i]/DBx[iB,i]
            bondAngles.append(npacosd(np.sum(DAv*DBv)))
    return np.array(bondAngles)
    
################################################################################
#### Reflectometry calculations
################################################################################

def REFDRefine(Profile,ProfDict,Inst,Limits,Substances,data):
    G2fil.G2Print ('fit REFD data by '+data['Minimizer']+' using %.2f%% data resolution'%(data['Resolution'][0]))
    
    class RandomDisplacementBounds(object):
        """random displacement with bounds"""
        def __init__(self, xmin, xmax, stepsize=0.5):
            self.xmin = xmin
            self.xmax = xmax
            self.stepsize = stepsize
    
        def __call__(self, x):
            """take a random step but ensure the new position is within the bounds"""
            while True:
                # this could be done in a much more clever way, but it will work for example purposes
                steps = self.xmax-self.xmin
                xnew = x + np.random.uniform(-self.stepsize*steps, self.stepsize*steps, np.shape(x))
                if np.all(xnew < self.xmax) and np.all(xnew > self.xmin):
                    break
            return xnew
    
    def GetModelParms():
        parmDict = {}
        varyList = []
        values = []
        bounds = []
        parmDict['dQ type'] = data['dQ type']
        parmDict['Res'] = data['Resolution'][0]/(100.*sateln2)     #% FWHM-->decimal sig
        for parm in ['Scale','FltBack']:
            parmDict[parm] = data[parm][0]
            if data[parm][1]:
                varyList.append(parm)
                values.append(data[parm][0])
                bounds.append(Bounds[parm])
        parmDict['Layer Seq'] = np.array(['0',]+data['Layer Seq'].split()+[str(len(data['Layers'])-1),],dtype=int)
        parmDict['nLayers'] = len(parmDict['Layer Seq'])
        for ilay,layer in enumerate(data['Layers']):
            name = layer['Name']
            cid = str(ilay)+';'
            parmDict[cid+'Name'] = name
            for parm in ['Thick','Rough','DenMul','Mag SLD','iDenMul']:
                parmDict[cid+parm] = layer.get(parm,[0.,False])[0]
                if layer.get(parm,[0.,False])[1]:
                    varyList.append(cid+parm)
                    value = layer[parm][0]
                    values.append(value)
                    if value:
                        bound = [value*Bfac,value/Bfac]
                    else:
                        bound = [0.,10.]
                    bounds.append(bound)
            if name not in ['vacuum','unit scatter']:
                parmDict[cid+'rho'] = Substances[name]['Scatt density']
                parmDict[cid+'irho'] = Substances[name].get('XImag density',0.)
        return parmDict,varyList,values,bounds
    
    def SetModelParms():
        line = ' Refined parameters: Histogram scale: %.4g'%(parmDict['Scale'])
        if 'Scale' in varyList:
            data['Scale'][0] = parmDict['Scale']
            line += ' esd: %.4g'%(sigDict['Scale'])                                                             
        G2fil.G2Print (line)
        line = ' Flat background: %15.4g'%(parmDict['FltBack'])
        if 'FltBack' in varyList:
            data['FltBack'][0] = parmDict['FltBack']
            line += ' esd: %15.3g'%(sigDict['FltBack'])
        G2fil.G2Print (line)
        for ilay,layer in enumerate(data['Layers']):
            name = layer['Name']
            G2fil.G2Print (' Parameters for layer: %d %s'%(ilay,name))
            cid = str(ilay)+';'
            line = ' '
            line2 = ' Scattering density: Real %.5g'%(Substances[name]['Scatt density']*parmDict[cid+'DenMul'])
            line2 += ' Imag %.5g'%(Substances[name].get('XImag density',0.)*parmDict[cid+'DenMul'])
            for parm in ['Thick','Rough','DenMul','Mag SLD','iDenMul']:
                if parm in layer:
                    if parm == 'Rough':
                        layer[parm][0] = abs(parmDict[cid+parm])    #make positive
                    else:
                        layer[parm][0] = parmDict[cid+parm]
                    line += ' %s: %.3f'%(parm,layer[parm][0])
                    if cid+parm in varyList:
                        line += ' esd: %.3g'%(sigDict[cid+parm])
            G2fil.G2Print (line)
            G2fil.G2Print (line2)
    
    def calcREFD(values,Q,Io,wt,Qsig,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        M = np.sqrt(wt)*(getREFD(Q,Qsig,parmDict)-Io)
        return M
    
    def sumREFD(values,Q,Io,wt,Qsig,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        M = np.sqrt(wt)*(getREFD(Q,Qsig,parmDict)-Io)
        return np.sum(M**2)
    
    def getREFD(Q,Qsig,parmDict):
        Ic = np.ones_like(Q)*parmDict['FltBack']
        Scale = parmDict['Scale']
        Nlayers = parmDict['nLayers']
        Res = parmDict['Res']
        depth = np.zeros(Nlayers)
        rho = np.zeros(Nlayers)
        irho = np.zeros(Nlayers)
        sigma = np.zeros(Nlayers)
        for ilay,lay in enumerate(parmDict['Layer Seq']):
            cid = str(lay)+';'
            depth[ilay] = parmDict[cid+'Thick']
            sigma[ilay] = parmDict[cid+'Rough']
            if parmDict[cid+'Name'] == u'unit scatter':
                rho[ilay] = parmDict[cid+'DenMul']
                irho[ilay] = parmDict[cid+'iDenMul']
            elif 'vacuum' != parmDict[cid+'Name']:
                rho[ilay] = parmDict[cid+'rho']*parmDict[cid+'DenMul']
                irho[ilay] = parmDict[cid+'irho']*parmDict[cid+'DenMul']
            if cid+'Mag SLD' in parmDict:
                rho[ilay] += parmDict[cid+'Mag SLD']
        if parmDict['dQ type'] == 'None':
            AB = abeles(0.5*Q,depth,rho,irho,sigma[1:])     #Q --> k, offset roughness for abeles
        elif 'const' in parmDict['dQ type']:
            AB = SmearAbeles(0.5*Q,Q*Res,depth,rho,irho,sigma[1:])
        else:       #dQ/Q in data
            AB = SmearAbeles(0.5*Q,Qsig,depth,rho,irho,sigma[1:])
        Ic += AB*Scale
        return Ic
        
    def estimateT0(takestep):
        Mmax = -1.e-10
        Mmin = 1.e10
        for i in range(100):
            x0 = takestep(values)
            M = sumREFD(x0,Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList)
            Mmin = min(M,Mmin)
            MMax = max(M,Mmax)
        return 1.5*(MMax-Mmin)

    Q,Io,wt,Ic,Ib,Qsig = Profile[:6]
    if data.get('2% weight'):
        wt = 1./(0.02*Io)**2
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    wtFactor = ProfDict['wtFactor']
    Bfac = data['Toler']
    Ibeg = np.searchsorted(Q,Qmin)
    Ifin = np.searchsorted(Q,Qmax)+1    #include last point
    Ic[:] = 0
    Bounds = {'Scale':[data['Scale'][0]*Bfac,data['Scale'][0]/Bfac],'FltBack':[0.,1.e-6],
              'DenMul':[0.,1.],'Thick':[1.,500.],'Rough':[0.,10.],'Mag SLD':[-10.,10.],'iDenMul':[-1.,1.]}
    parmDict,varyList,values,bounds = GetModelParms()
    Msg = 'Failed to converge'
    if varyList:
        if data['Minimizer'] == 'LMLS': 
            result = so.leastsq(calcREFD,values,full_output=True,epsfcn=1.e-8,ftol=1.e-6,
                args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList))
            parmDict.update(zip(varyList,result[0]))
            chisq = np.sum(result[2]['fvec']**2)
            ncalc = result[2]['nfev']
            covM = result[1]
            newVals = result[0]
        elif data['Minimizer'] == 'Basin Hopping':
            xyrng = np.array(bounds).T
            take_step = RandomDisplacementBounds(xyrng[0], xyrng[1])
            T0 = estimateT0(take_step)
            G2fil.G2Print (' Estimated temperature: %.3g'%(T0))
            result = so.basinhopping(sumREFD,values,take_step=take_step,disp=True,T=T0,stepsize=Bfac,
                interval=20,niter=200,minimizer_kwargs={'method':'L-BFGS-B','bounds':bounds,
                'args':(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList)})
            chisq = result.fun
            ncalc = result.nfev
            newVals = result.x
            covM = []
        elif data['Minimizer'] == 'MC/SA Anneal':
            xyrng = np.array(bounds).T
            result = G2mth.anneal(sumREFD, values, 
                args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList),
                schedule='log', full_output=True,maxeval=None, maxaccept=None, maxiter=10,dwell=1000,
                boltzmann=10.0, feps=1e-6,lower=xyrng[0], upper=xyrng[1], slope=0.9,ranStart=True,
                ranRange=0.20,autoRan=False,dlg=None)
            newVals = result[0]
            parmDict.update(zip(varyList,newVals))
            chisq = result[1]
            ncalc = result[3]
            covM = []
            G2fil.G2Print (' MC/SA final temperature: %.4g'%(result[2]))
        elif data['Minimizer'] == 'L-BFGS-B':
            result = so.minimize(sumREFD,values,method='L-BFGS-B',bounds=bounds,   #ftol=Ftol,
                args=(Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList))
            parmDict.update(zip(varyList,result['x']))
            chisq = result.fun
            ncalc = result.nfev
            newVals = result.x
            covM = []
    else:   #nothing varied
        M = calcREFD(values,Q[Ibeg:Ifin],Io[Ibeg:Ifin],wtFactor*wt[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict,varyList)
        chisq = np.sum(M**2)
        ncalc = 0
        covM = []
        sig = []
        sigDict = {}
        result = []
    Rvals = {}
    Rvals['Rwp'] = np.sqrt(chisq/np.sum(wt[Ibeg:Ifin]*Io[Ibeg:Ifin]**2))*100.      #to %
    Rvals['GOF'] = chisq/(Ifin-Ibeg-len(varyList))       #reduced chi^2
    Ic[Ibeg:Ifin] = getREFD(Q[Ibeg:Ifin],Qsig[Ibeg:Ifin],parmDict)
    Ib[Ibeg:Ifin] = parmDict['FltBack']
    try:
        if not len(varyList):
            Msg += ' - nothing refined'
            raise ValueError
        Nans = np.isnan(newVals)
        if np.any(Nans):
            name = varyList[Nans.nonzero(True)[0]]
            Msg += ' Nan result for '+name+'!'
            raise ValueError
        Negs = np.less_equal(newVals,0.)
        if np.any(Negs):
            indx = Negs.nonzero()
            name = varyList[indx[0][0]]
            if name != 'FltBack' and name.split(';')[1] in ['Thick',]:
                Msg += ' negative coefficient for '+name+'!'
                raise ValueError
        if len(covM):
            sig = np.sqrt(np.diag(covM)*Rvals['GOF'])
            covMatrix = covM*Rvals['GOF']
        else:
            sig = np.zeros(len(varyList))
            covMatrix = []
        sigDict = dict(zip(varyList,sig))
        G2fil.G2Print (' Results of reflectometry data modelling fit:')
        G2fil.G2Print ('Number of function calls: %d Number of observations: %d Number of parameters: %d'%(ncalc,Ifin-Ibeg,len(varyList)))
        G2fil.G2Print ('Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],chisq,Rvals['GOF']))
        SetModelParms()
        return True,result,varyList,sig,Rvals,covMatrix,parmDict,''
    except (ValueError,TypeError):      #when bad LS refinement; covM missing or with nans
        G2fil.G2Print (Msg)
        return False,0,0,0,0,0,0,Msg
        
def makeSLDprofile(data,Substances):
    
    sq2 = np.sqrt(2.)
    laySeq = ['0',]+data['Layer Seq'].split()+[str(len(data['Layers'])-1),]
    Nlayers = len(laySeq)
    laySeq = np.array(laySeq,dtype=int)
    interfaces = np.zeros(Nlayers)
    rho = np.zeros(Nlayers)
    sigma = np.zeros(Nlayers)
    name = data['Layers'][0]['Name']
    thick = 0.
    for ilay,lay in enumerate(laySeq):
        layer = data['Layers'][lay]
        name = layer['Name']
        if 'Thick' in layer:
            thick += layer['Thick'][0]
            interfaces[ilay] = layer['Thick'][0]+interfaces[ilay-1]
        if 'Rough' in layer:
            sigma[ilay] = max(0.001,layer['Rough'][0])
        if name != 'vacuum':
            if name == 'unit scatter':
                rho[ilay] = np.sqrt(layer['DenMul'][0]**2+layer['iDenMul'][0]**2)
            else:
                rrho = Substances[name]['Scatt density']
                irho = Substances[name]['XImag density']
                rho[ilay] = np.sqrt(rrho**2+irho**2)*layer['DenMul'][0]
        if 'Mag SLD' in layer:
            rho[ilay] += layer['Mag SLD'][0]
    name = data['Layers'][-1]['Name']
    x = np.linspace(-0.15*thick,1.15*thick,1000,endpoint=True)
    xr = np.flipud(x)
    interfaces[-1] = x[-1]
    y = np.ones_like(x)*rho[0]
    iBeg = 0
    for ilayer in range(Nlayers-1):
        delt = rho[ilayer+1]-rho[ilayer]
        iPos = np.searchsorted(x,interfaces[ilayer])
        y[iBeg:] += (delt/2.)*sp.erfc((interfaces[ilayer]-x[iBeg:])/(sq2*sigma[ilayer+1]))
        iBeg = iPos
    return x,xr,y    

def REFDModelFxn(Profile,Inst,Limits,Substances,data):
    
    Q,Io,wt,Ic,Ib,Qsig = Profile[:6]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    iBeg = np.searchsorted(Q,Qmin)
    iFin = np.searchsorted(Q,Qmax)+1    #include last point
    Ib[:] = data['FltBack'][0]
    Ic[:] = 0
    Scale = data['Scale'][0]
    if data['Layer Seq'] == []:
        return
    laySeq = ['0',]+data['Layer Seq'].split()+[str(len(data['Layers'])-1),]
    Nlayers = len(laySeq)
    depth = np.zeros(Nlayers)
    rho = np.zeros(Nlayers)
    irho = np.zeros(Nlayers)
    sigma = np.zeros(Nlayers)
    for ilay,lay in enumerate(np.array(laySeq,dtype=int)):
        layer = data['Layers'][lay]
        name = layer['Name']
        if 'Thick' in layer:    #skips first & last layers
            depth[ilay] = layer['Thick'][0]
        if 'Rough' in layer:    #skips first layer
            sigma[ilay] = layer['Rough'][0]
        if 'unit scatter' == name:
            rho[ilay] = layer['DenMul'][0]
            irho[ilay] = layer['iDenMul'][0]
        else:
            rho[ilay] = Substances[name]['Scatt density']*layer['DenMul'][0]
            irho[ilay] = Substances[name].get('XImag density',0.)*layer['DenMul'][0]
        if 'Mag SLD' in layer:
            rho[ilay] += layer['Mag SLD'][0]
    if data['dQ type'] == 'None':
        AB = abeles(0.5*Q[iBeg:iFin],depth,rho,irho,sigma[1:])     #Q --> k, offset roughness for abeles
    elif 'const' in data['dQ type']:
        res = data['Resolution'][0]/(100.*sateln2)
        AB = SmearAbeles(0.5*Q[iBeg:iFin],res*Q[iBeg:iFin],depth,rho,irho,sigma[1:])
    else:       #dQ/Q in data
        AB = SmearAbeles(0.5*Q[iBeg:iFin],Qsig[iBeg:iFin],depth,rho,irho,sigma[1:])
    Ic[iBeg:iFin] = AB*Scale+Ib[iBeg:iFin]

def abeles(kz, depth, rho, irho=0, sigma=0):
    """
    Optical matrix form of the reflectivity calculation.
    O.S. Heavens, Optical Properties of Thin Solid Films
    
    Reflectometry as a function of kz for a set of slabs.

    :param kz: float[n] (1/Ang). Scattering vector, :math:`2\\pi\\sin(\\theta)/\\lambda`.
        This is :math:`\\tfrac12 Q_z`.        
    :param depth:  float[m] (Ang).
        thickness of each layer.  The thickness of the incident medium
        and substrate are ignored.
    :param rho:  float[n,k] (1e-6/Ang^2)
        Real scattering length density for each layer for each kz
    :param irho:  float[n,k] (1e-6/Ang^2)
        Imaginary scattering length density for each layer for each kz
        Note: absorption cross section mu = 2 irho/lambda for neutrons
    :param sigma: float[m-1] (Ang)
        interfacial roughness.  This is the roughness between a layer
        and the previous layer. The sigma array should have m-1 entries.

    Slabs are ordered with the surface SLD at index 0 and substrate at
    index -1, or reversed if kz < 0.
    """
    def calc(kz, depth, rho, irho, sigma):
        if len(kz) == 0: return kz
    
        # Complex index of refraction is relative to the incident medium.
        # We can get the same effect using kz_rel^2 = kz^2 + 4*pi*rho_o
        # in place of kz^2, and ignoring rho_o
        kz_sq = kz**2 + 4e-6*np.pi*rho[:,0]
        k = kz
    
        # According to Heavens, the initial matrix should be [ 1 F; F 1],
        # which we do by setting B=I and M0 to [1 F; F 1].  An extra matrix
        # multiply versus some coding convenience.
        B11 = 1
        B22 = 1
        B21 = 0
        B12 = 0
        for i in range(0, len(depth)-1):
            k_next = np.sqrt(kz_sq - 4e-6*np.pi*(rho[:,i+1] + 1j*irho[:,i+1]))
            F = (k - k_next) / (k + k_next)
            F *= np.exp(-2*k*k_next*sigma[i]**2)
            #print "==== layer",i
            #print "kz:", kz
            #print "k:", k
            #print "k_next:",k_next
            #print "F:",F
            #print "rho:",rho[:,i+1]
            #print "irho:",irho[:,i+1]
            #print "d:",depth[i],"sigma:",sigma[i]
            M11 = np.exp(1j*k*depth[i]) if i>0 else 1
            M22 = np.exp(-1j*k*depth[i]) if i>0 else 1
            M21 = F*M11
            M12 = F*M22
            C1 = B11*M11 + B21*M12
            C2 = B11*M21 + B21*M22
            B11 = C1
            B21 = C2
            C1 = B12*M11 + B22*M12
            C2 = B12*M21 + B22*M22
            B12 = C1
            B22 = C2
            k = k_next
    
        r = B12/B11
        return np.absolute(r)**2

    if np.isscalar(kz): kz = np.asarray([kz], 'd')

    m = len(depth)

    # Make everything into arrays
    depth = np.asarray(depth,'d')
    rho = np.asarray(rho,'d')
    irho = irho*np.ones_like(rho) if np.isscalar(irho) else np.asarray(irho,'d')
    sigma = sigma*np.ones(m-1,'d') if np.isscalar(sigma) else np.asarray(sigma,'d')

    # Repeat rho,irho columns as needed
    if len(rho.shape) == 1:
        rho = rho[None,:]
        irho = irho[None,:]

    return calc(kz, depth, rho, irho, sigma)
    
def SmearAbeles(kz,dq, depth, rho, irho=0, sigma=0):
    y = abeles(kz, depth, rho, irho, sigma)
    s = dq/2.
    y += 0.1354*(abeles(kz+2*s, depth, rho, irho, sigma)+abeles(kz-2*s, depth, rho, irho, sigma))
    y += 0.24935*(abeles(kz-5*s/3., depth, rho, irho, sigma)+abeles(kz+5*s/3., depth, rho, irho, sigma)) 
    y += 0.4111*(abeles(kz-4*s/3., depth, rho, irho, sigma)+abeles(kz+4*s/3., depth, rho, irho, sigma)) 
    y += 0.60653*(abeles(kz-s, depth, rho, irho, sigma) +abeles(kz+s, depth, rho, irho, sigma))
    y += 0.80074*(abeles(kz-2*s/3., depth, rho, irho, sigma)+abeles(kz-2*s/3., depth, rho, irho, sigma))
    y += 0.94596*(abeles(kz-s/3., depth, rho, irho, sigma)+abeles(kz-s/3., depth, rho, irho, sigma))
    y *= 0.137023
    return y
        
def makeRefdFFT(Limits,Profile):
    G2fil.G2Print ('make fft')
    Q,Io = Profile[:2]
    Qmin = Limits[1][0]
    Qmax = Limits[1][1]
    iBeg = np.searchsorted(Q,Qmin)
    iFin = np.searchsorted(Q,Qmax)+1    #include last point
    Qf = np.linspace(0.,Q[iFin-1],5000)
    QI = si.interp1d(Q[iBeg:iFin],Io[iBeg:iFin],bounds_error=False,fill_value=0.0)
    If = QI(Qf)*Qf**4
    R = np.linspace(0.,5000.,5000)
    F = fft.rfft(If)
    return R,F

    
################################################################################
#### Stacking fault simulation codes
################################################################################

def GetStackParms(Layers):
    
    Parms = []
#cell parms
    if Layers['Laue'] in ['-3','-3m','4/m','4/mmm','6/m','6/mmm']:
        Parms.append('cellA')
        Parms.append('cellC')
    else:
        Parms.append('cellA')
        Parms.append('cellB')
        Parms.append('cellC')
        if Layers['Laue'] != 'mmm':
            Parms.append('cellG')
#Transition parms
    for iY in range(len(Layers['Layers'])):
        for iX in range(len(Layers['Layers'])):
            Parms.append('TransP;%d;%d'%(iY,iX))
            Parms.append('TransX;%d;%d'%(iY,iX))
            Parms.append('TransY;%d;%d'%(iY,iX))
            Parms.append('TransZ;%d;%d'%(iY,iX))
    return Parms

def StackSim(Layers,ctrls,scale=0.,background={},limits=[],inst={},profile=[]):
    '''Simulate powder or selected area diffraction pattern from stacking faults using DIFFaX
    
    :param dict Layers: dict with following items

      ::

       {'Laue':'-1','Cell':[False,1.,1.,1.,90.,90.,90,1.],
       'Width':[[10.,10.],[False,False]],'Toler':0.01,'AtInfo':{},
       'Layers':[],'Stacking':[],'Transitions':[]}
        
    :param str ctrls: controls string to be written on DIFFaX controls.dif file
    :param float scale: scale factor
    :param dict background: background parameters
    :param list limits: min/max 2-theta to be calculated
    :param dict inst: instrument parameters dictionary
    :param list profile: powder pattern data
    
    Note that parameters all updated in place    
    '''
    import atmdata
    path = sys.path
    for name in path:
        if 'bin' in name:
            DIFFaX = name+'/DIFFaX.exe'
            G2fil.G2Print (' Execute '+DIFFaX)
            break
    # make form factor file that DIFFaX wants - atom types are GSASII style
    sf = open('data.sfc','w')
    sf.write('GSASII special form factor file for DIFFaX\n\n')
    atTypes = list(Layers['AtInfo'].keys())
    if 'H' not in atTypes:
        atTypes.insert(0,'H')
    for atType in atTypes:
        if atType == 'H': 
            blen = -.3741
        else:
            blen = Layers['AtInfo'][atType]['Isotopes']['Nat. Abund.']['SL'][0]
        Adat = atmdata.XrayFF[atType]
        text = '%4s'%(atType.ljust(4))
        for i in range(4):
            text += '%11.6f%11.6f'%(Adat['fa'][i],Adat['fb'][i])
        text += '%11.6f%11.6f'%(Adat['fc'],blen)
        text += '%3d\n'%(Adat['Z'])
        sf.write(text)
    sf.close()
    #make DIFFaX control.dif file - future use GUI to set some of these flags
    cf = open('control.dif','w')
    if ctrls == '0\n0\n3\n' or ctrls == '0\n1\n3\n': 
        x0 = profile[0]
        iBeg = np.searchsorted(x0,limits[0])
        iFin = np.searchsorted(x0,limits[1])+1
        if iFin-iBeg > 20000:
            iFin = iBeg+20000
        Dx = (x0[iFin]-x0[iBeg])/(iFin-iBeg)
        cf.write('GSASII-DIFFaX.dat\n'+ctrls)
        cf.write('%.6f %.6f %.6f\n1\n1\nend\n'%(x0[iBeg],x0[iFin],Dx))
    else:
        cf.write('GSASII-DIFFaX.dat\n'+ctrls)
        inst = {'Type':['XSC','XSC',]}
    cf.close()
    #make DIFFaX data file
    df = open('GSASII-DIFFaX.dat','w')
    df.write('INSTRUMENTAL\n')
    if 'X' in inst['Type'][0]:
        df.write('X-RAY\n')
    elif 'N' in inst['Type'][0]:
        df.write('NEUTRON\n')
    if ctrls == '0\n0\n3\n' or ctrls == '0\n1\n3\n': 
        df.write('%.4f\n'%(G2mth.getMeanWave(inst)))
        U = ateln2*inst['U'][1]/10000.
        V = ateln2*inst['V'][1]/10000.
        W = ateln2*inst['W'][1]/10000.
        HWHM = U*nptand(x0[iBeg:iFin]/2.)**2+V*nptand(x0[iBeg:iFin]/2.)+W
        HW = np.sqrt(np.mean(HWHM))
    #    df.write('PSEUDO-VOIGT 0.015 -0.0036 0.009 0.605 TRIM\n')
        if 'Mean' in Layers['selInst']:
            df.write('GAUSSIAN %.6f TRIM\n'%(HW))     #fast option - might not really matter
        elif 'Gaussian' in Layers['selInst']:
            df.write('GAUSSIAN %.6f %.6f %.6f TRIM\n'%(U,V,W))    #slow - make a GUI option?
        else:
            df.write('None\n')
    else:
        df.write('0.10\nNone\n')
    df.write('STRUCTURAL\n')
    a,b,c = Layers['Cell'][1:4]
    gam = Layers['Cell'][6]
    df.write('%.4f %.4f %.4f %.3f\n'%(a,b,c,gam))
    laue = Layers['Laue']
    if laue == '2/m(ab)':
        laue = '2/m(1)'
    elif laue == '2/m(c)':
        laue = '2/m(2)'
    if 'unknown' in Layers['Laue']:
        df.write('%s %.3f\n'%(laue,Layers['Toler']))
    else:    
        df.write('%s\n'%(laue))
    df.write('%d\n'%(len(Layers['Layers'])))
    if Layers['Width'][0][0] < 1. or Layers['Width'][0][1] < 1.:
        df.write('%.1f %.1f\n'%(Layers['Width'][0][0]*10000.,Layers['Width'][0][0]*10000.))    #mum to A
    layerNames = []
    for layer in Layers['Layers']:
        layerNames.append(layer['Name'])
    for il,layer in enumerate(Layers['Layers']):
        if layer['SameAs']:
            df.write('LAYER %d = %d\n'%(il+1,layerNames.index(layer['SameAs'])+1))
            continue
        df.write('LAYER %d\n'%(il+1))
        if '-1' in layer['Symm']:
            df.write('CENTROSYMMETRIC\n')
        else:
            df.write('NONE\n')
        for ia,atom in enumerate(layer['Atoms']):
            [name,atype,x,y,z,frac,Uiso] = atom
            if '-1' in layer['Symm'] and [x,y,z] == [0.,0.,0.]:
                frac /= 2.
            df.write('%4s %3d %.5f %.5f %.5f %.4f %.2f\n'%(atype.ljust(6),ia,x,y,z,78.9568*Uiso,frac))
    df.write('STACKING\n')
    df.write('%s\n'%(Layers['Stacking'][0]))
    if 'recursive' in Layers['Stacking'][0]:
        df.write('%s\n'%Layers['Stacking'][1])
    else:
        if 'list' in Layers['Stacking'][1]:
            Slen = len(Layers['Stacking'][2])
            iB = 0
            iF = 0
            while True:
                iF += 68
                if iF >= Slen:
                    break
                iF = min(iF,Slen)
                df.write('%s\n'%(Layers['Stacking'][2][iB:iF]))
                iB = iF
        else:
            df.write('%s\n'%Layers['Stacking'][1])    
    df.write('TRANSITIONS\n')
    for iY in range(len(Layers['Layers'])):
        sumPx = 0.
        for iX in range(len(Layers['Layers'])):
            p,dx,dy,dz = Layers['Transitions'][iY][iX][:4]
            p = round(p,3)
            df.write('%.3f %.5f %.5f %.5f\n'%(p,dx,dy,dz))
            sumPx += p
        if sumPx != 1.0:    #this has to be picky since DIFFaX is.
            G2fil.G2Print ('ERROR - Layer probabilities sum to %.3f DIFFaX will insist it = 1.0'%sumPx)
            df.close()
            os.remove('data.sfc')
            os.remove('control.dif')
            os.remove('GSASII-DIFFaX.dat')
            return       
    df.close()    
    time0 = time.time()
    try:
        subp.call(DIFFaX)
    except OSError:
        G2fil.G2Print('DIFFax.exe is not available for this platform',mode='warn')
    G2fil.G2Print (' DIFFaX time = %.2fs'%(time.time()-time0))
    if os.path.exists('GSASII-DIFFaX.spc'):
        Xpat = np.loadtxt('GSASII-DIFFaX.spc').T
        iFin = iBeg+Xpat.shape[1]
        bakType,backDict,backVary = SetBackgroundParms(background)
        backDict['Lam1'] = G2mth.getWave(inst)
        profile[4][iBeg:iFin] = getBackground('',backDict,bakType,inst['Type'][0],profile[0][iBeg:iFin])[0]    
        profile[3][iBeg:iFin] = Xpat[-1]*scale+profile[4][iBeg:iFin]
        if not np.any(profile[1]):                   #fill dummy data x,y,w,yc,yb,yd
            rv = st.poisson(profile[3][iBeg:iFin])
            profile[1][iBeg:iFin] = rv.rvs()
            Z = np.ones_like(profile[3][iBeg:iFin])
            Z[1::2] *= -1
            profile[1][iBeg:iFin] = profile[3][iBeg:iFin]+np.abs(profile[1][iBeg:iFin]-profile[3][iBeg:iFin])*Z
            profile[2][iBeg:iFin] = np.where(profile[1][iBeg:iFin]>0.,1./profile[1][iBeg:iFin],1.0)
        profile[5][iBeg:iFin] = profile[1][iBeg:iFin]-profile[3][iBeg:iFin]
    #cleanup files..
        os.remove('GSASII-DIFFaX.spc')
    elif os.path.exists('GSASII-DIFFaX.sadp'):
        Sadp = np.fromfile('GSASII-DIFFaX.sadp','>u2')
        Sadp = np.reshape(Sadp,(256,-1))
        Layers['Sadp']['Img'] = Sadp
        os.remove('GSASII-DIFFaX.sadp')
    os.remove('data.sfc')
    os.remove('control.dif')
    os.remove('GSASII-DIFFaX.dat')
    
def SetPWDRscan(inst,limits,profile):
    
    wave = G2mth.getMeanWave(inst)
    x0 = profile[0]
    iBeg = np.searchsorted(x0,limits[0])
    iFin = np.searchsorted(x0,limits[1])
    if iFin-iBeg > 20000:
        iFin = iBeg+20000
    Dx = (x0[iFin]-x0[iBeg])/(iFin-iBeg)
    pyx.pygetinst(wave,x0[iBeg],x0[iFin],Dx)
    return iFin-iBeg
       
def SetStackingSF(Layers,debug):
# Load scattering factors into DIFFaX arrays
    import atmdata
    atTypes = Layers['AtInfo'].keys()
    aTypes = []
    for atype in atTypes:
        aTypes.append('%4s'%(atype.ljust(4)))
    SFdat = []
    for atType in atTypes:
        Adat = atmdata.XrayFF[atType]
        SF = np.zeros(9)
        SF[:8:2] = Adat['fa']
        SF[1:8:2] = Adat['fb']
        SF[8] = Adat['fc']
        SFdat.append(SF)
    SFdat = np.array(SFdat)
    pyx.pyloadscf(len(atTypes),aTypes,SFdat.T,debug)
    
def SetStackingClay(Layers,Type):
# Controls
    rand.seed()
    ranSeed = rand.randint(1,2**16-1)
    try:    
        laueId = ['-1','2/m(ab)','2/m(c)','mmm','-3','-3m','4/m','4/mmm',
            '6/m','6/mmm'].index(Layers['Laue'])+1
    except ValueError:  #for 'unknown'
        laueId = -1
    if 'SADP' in Type:
        planeId = ['h0l','0kl','hhl','h-hl'].index(Layers['Sadp']['Plane'])+1
        lmax = int(Layers['Sadp']['Lmax'])
    else:
        planeId = 0
        lmax = 0
# Sequences
    StkType = ['recursive','explicit'].index(Layers['Stacking'][0])
    try:
        StkParm = ['infinite','random','list'].index(Layers['Stacking'][1])
    except ValueError:
        StkParm = -1
    if StkParm == 2:    #list
        StkSeq = [int(val) for val in Layers['Stacking'][2].split()]
        Nstk = len(StkSeq)
    else:
        Nstk = 1
        StkSeq = [0,]
    if StkParm == -1:
        StkParm = int(Layers['Stacking'][1])
    Wdth = Layers['Width'][0]
    mult = 1
    controls = [laueId,planeId,lmax,mult,StkType,StkParm,ranSeed]
    LaueSym = Layers['Laue'].ljust(12)
    pyx.pygetclay(controls,LaueSym,Wdth,Nstk,StkSeq)
    return laueId,controls
    
def SetCellAtoms(Layers):
    Cell = Layers['Cell'][1:4]+Layers['Cell'][6:7]
# atoms in layers
    atTypes = list(Layers['AtInfo'].keys())
    AtomXOU = []
    AtomTp = []
    LayerSymm = []
    LayerNum = []
    layerNames = []
    Natm = 0
    Nuniq = 0
    for layer in Layers['Layers']:
        layerNames.append(layer['Name'])
    for il,layer in enumerate(Layers['Layers']):
        if layer['SameAs']:
            LayerNum.append(layerNames.index(layer['SameAs'])+1)
            continue
        else:
            LayerNum.append(il+1)
            Nuniq += 1
        if '-1' in layer['Symm']:
            LayerSymm.append(1)
        else:
            LayerSymm.append(0)
        for ia,atom in enumerate(layer['Atoms']):
            [name,atype,x,y,z,frac,Uiso] = atom
            Natm += 1
            AtomTp.append('%4s'%(atype.ljust(4)))
            Ta = atTypes.index(atype)+1
            AtomXOU.append([float(Nuniq),float(ia+1),float(Ta),x,y,z,frac,Uiso*78.9568])
    AtomXOU = np.array(AtomXOU)
    Nlayers = len(layerNames)
    pyx.pycellayer(Cell,Natm,AtomTp,AtomXOU.T,Nuniq,LayerSymm,Nlayers,LayerNum)
    return Nlayers
    
def SetStackingTrans(Layers,Nlayers):
# Transitions
    TransX = []
    TransP = []
    for Ytrans in Layers['Transitions']:
        TransP.append([trans[0] for trans in Ytrans])   #get just the numbers
        TransX.append([trans[1:4] for trans in Ytrans])   #get just the numbers
    TransP = np.array(TransP,dtype='float').T
    TransX = np.array(TransX,dtype='float')
#    GSASIIpath.IPyBreak()
    pyx.pygettrans(Nlayers,TransP,TransX)
    
def CalcStackingPWDR(Layers,scale,background,limits,inst,profile,debug):
# Scattering factors
    SetStackingSF(Layers,debug)
# Controls & sequences
    laueId,controls = SetStackingClay(Layers,'PWDR')
# cell & atoms
    Nlayers = SetCellAtoms(Layers)
    Volume = Layers['Cell'][7]    
# Transitions
    SetStackingTrans(Layers,Nlayers)
# PWDR scan
    Nsteps = SetPWDRscan(inst,limits,profile)
# result as Spec
    x0 = profile[0]
    profile[3] = np.zeros(len(profile[0]))
    profile[4] = np.zeros(len(profile[0]))
    profile[5] = np.zeros(len(profile[0]))
    iBeg = np.searchsorted(x0,limits[0])
    iFin = np.searchsorted(x0,limits[1])+1
    if iFin-iBeg > 20000:
        iFin = iBeg+20000
    Nspec = 20001       
    spec = np.zeros(Nspec,dtype='double')    
    time0 = time.time()
    pyx.pygetspc(controls,Nspec,spec)
    G2fil.G2Print (' GETSPC time = %.2fs'%(time.time()-time0))
    time0 = time.time()
    U = ateln2*inst['U'][1]/10000.
    V = ateln2*inst['V'][1]/10000.
    W = ateln2*inst['W'][1]/10000.
    HWHM = U*nptand(x0[iBeg:iFin]/2.)**2+V*nptand(x0[iBeg:iFin]/2.)+W
    HW = np.sqrt(np.mean(HWHM))
    BrdSpec = np.zeros(Nsteps)
    if 'Mean' in Layers['selInst']:
        pyx.pyprofile(U,V,W,HW,1,Nsteps,BrdSpec)
    elif 'Gaussian' in Layers['selInst']:
        pyx.pyprofile(U,V,W,HW,4,Nsteps,BrdSpec)
    else:
        BrdSpec = spec[:Nsteps]
    BrdSpec /= Volume
    iFin = iBeg+Nsteps
    bakType,backDict,backVary = SetBackgroundParms(background)
    backDict['Lam1'] = G2mth.getWave(inst)
    profile[4][iBeg:iFin] = getBackground('',backDict,bakType,inst['Type'][0],profile[0][iBeg:iFin])[0]    
    profile[3][iBeg:iFin] = BrdSpec*scale+profile[4][iBeg:iFin]
    if not np.any(profile[1]):                   #fill dummy data x,y,w,yc,yb,yd
        try:
            rv = st.poisson(profile[3][iBeg:iFin])
            profile[1][iBeg:iFin] = rv.rvs()
        except ValueError:
            profile[1][iBeg:iFin] = profile[3][iBeg:iFin]
        Z = np.ones_like(profile[3][iBeg:iFin])
        Z[1::2] *= -1
        profile[1][iBeg:iFin] = profile[3][iBeg:iFin]+np.abs(profile[1][iBeg:iFin]-profile[3][iBeg:iFin])*Z
        profile[2][iBeg:iFin] = np.where(profile[1][iBeg:iFin]>0.,1./profile[1][iBeg:iFin],1.0)
    profile[5][iBeg:iFin] = profile[1][iBeg:iFin]-profile[3][iBeg:iFin]
    G2fil.G2Print (' Broadening time = %.2fs'%(time.time()-time0))
    
def CalcStackingSADP(Layers,debug):
    
# Scattering factors
    SetStackingSF(Layers,debug)
# Controls & sequences
    laueId,controls = SetStackingClay(Layers,'SADP')
# cell & atoms
    Nlayers = SetCellAtoms(Layers)    
# Transitions
    SetStackingTrans(Layers,Nlayers)
# result as Sadp
    Nspec = 20001       
    spec = np.zeros(Nspec,dtype='double')    
    time0 = time.time()
    hkLim,Incr,Nblk = pyx.pygetsadp(controls,Nspec,spec)
    Sapd = np.zeros((256,256))
    iB = 0
    for i in range(hkLim):
        iF = iB+Nblk
        p1 = 127+int(i*Incr)
        p2 = 128-int(i*Incr)
        if Nblk == 128:
            if i:
                Sapd[128:,p1] = spec[iB:iF]
                Sapd[:128,p1] = spec[iF:iB:-1]
            Sapd[128:,p2] = spec[iB:iF]
            Sapd[:128,p2] = spec[iF:iB:-1]
        else:
            if i:
                Sapd[:,p1] = spec[iB:iF]
            Sapd[:,p2] = spec[iB:iF]
        iB += Nblk
    Layers['Sadp']['Img'] = Sapd
    G2fil.G2Print (' GETSAD time = %.2fs'%(time.time()-time0))
    
###############################################################################
#### Maximum Entropy Method - Dysnomia
###############################################################################
    
def makePRFfile(data,MEMtype):
    ''' makes Dysnomia .prf control file from Dysnomia GUI controls
    
    :param dict data: GSAS-II phase data
    :param int MEMtype: 1 for neutron data with negative scattering lengths
                        0 otherwise
    :returns str: name of Dysnomia control file
    '''

    generalData = data['General']
    pName = generalData['Name'].replace(' ','_')
    DysData = data['Dysnomia']
    prfName = pName+'.prf'
    prf = open(prfName,'w')
    prf.write('$PREFERENCES\n')
    prf.write(pName+'.mem\n') #or .fos?
    prf.write(pName+'.out\n')
    prf.write(pName+'.pgrid\n')
    prf.write(pName+'.fba\n')
    prf.write(pName+'_eps.raw\n')
    prf.write('%d\n'%MEMtype)
    if DysData['DenStart'] == 'uniform':
        prf.write('0\n')
    else:
        prf.write('1\n')
    if DysData['Optimize'] == 'ZSPA':
        prf.write('0\n')
    else:
        prf.write('1\n')
    prf.write('1\n')
    if DysData['Lagrange'][0] == 'user':
        prf.write('0\n')
    else:
        prf.write('1\n')
    prf.write('%.4f %d\n'%(DysData['Lagrange'][1],DysData['wt pwr']))
    prf.write('%.3f\n'%DysData['Lagrange'][2])
    prf.write('%.2f\n'%DysData['E_factor'])
    prf.write('1\n')
    prf.write('0\n')
    prf.write('%d\n'%DysData['Ncyc'])
    prf.write('1\n')
    prf.write('1 0 0 0 0 0 0 0\n')
    if DysData['prior'] == 'uniform':
        prf.write('0\n')
    else:
        prf.write('1\n')
    prf.close()
    return prfName

def makeMEMfile(data,reflData,MEMtype,DYSNOMIA):
    ''' make Dysnomia .mem file of reflection data, etc.

    :param dict data: GSAS-II phase data
    :param list reflData: GSAS-II reflection data
    :param int MEMtype: 1 for neutron data with negative scattering lengths
                        0 otherwise
    :param str DYSNOMIA: path to dysnomia.exe
    '''
    
    DysData = data['Dysnomia']
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    A = G2lat.cell2A(cell)
    SGData = generalData['SGData']
    pName = generalData['Name'].replace(' ','_')
    memName = pName+'.mem'
    Map = generalData['Map']
    Type = Map['Type']
    UseList = Map['RefList']
    mem = open(memName,'w')
    mem.write('%s\n'%(generalData['Name']+' from '+UseList[0]))
    a,b,c,alp,bet,gam = cell
    mem.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'%(a,b,c,alp,bet,gam))
    mem.write('      0.0000000      0.0000000     -1    0    0    0     P\n')   #dummy PO stuff
    SGSym = generalData['SGData']['SpGrp']
    try:
        SGId = G2spc.spgbyNum.index(SGSym)
    except ValueError:
        return False
    org = 1
    if SGSym in G2spc.spg2origins:
        org = 2
    mapsize = Map['rho'].shape
    sumZ = 0.
    sumpos = 0.
    sumneg = 0.
    mem.write('%5d%5d%5d%5d%5d\n'%(SGId,org,mapsize[0],mapsize[1],mapsize[2]))
    for atm in generalData['NoAtoms']:
        Nat = generalData['NoAtoms'][atm]
        AtInfo = G2elem.GetAtomInfo(atm)
        sumZ += Nat*AtInfo['Z']
        isotope = generalData['Isotope'][atm]
        blen = generalData['Isotopes'][atm][isotope]['SL'][0]
        if blen < 0.:
            sumneg += blen*Nat
        else:
            sumpos += blen*Nat
    if 'X' in Type:
        mem.write('%10.2f  0.001\n'%sumZ)
    elif 'N' in Type and MEMtype:
        mem.write('%10.3f%10.3f 0.001\n'%(sumpos,sumneg))
    else:
        mem.write('%10.3f 0.001\n'%sumpos)
        
    dmin = DysData['MEMdmin']
    TOFlam = 2.0*dmin*npsind(80.0)
    refSet = G2lat.GenHLaue(dmin,SGData,A)      #list of h,k,l,d
    refDict = {'%d %d %d'%(ref[0],ref[1],ref[2]):ref for ref in refSet}
        
    refs = []
    prevpos = 0.
    for ref in reflData:
        if ref[3] < 0:
            continue
        if 'T' in Type:
            h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,x,x,x,prfo = ref[:16]
            s = np.sqrt(max(sig,0.0001))   #var -> sig in deg
            FWHM = getgamFW(gam,s)
            if dsp < dmin:
                continue
            theta = npasind(TOFlam/(2.*dsp))
            FWHM *= nptand(theta)/pos
            pos = 2.*theta
        else:
            h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,prfo = ref[:13]
            g = gam/100.    #centideg -> deg
            s = np.sqrt(max(sig,0.0001))/100.   #var -> sig in deg
            FWHM = getgamFW(g,s)
        delt = pos-prevpos
        refs.append([h,k,l,mult,pos,FWHM,Fobs,phase,delt])
        prevpos = pos
            
    ovlp = DysData['overlap']
    refs1 = []
    refs2 = []
    nref2 = 0
    iref = 0
    Nref = len(refs)
    start = False
    while iref < Nref-1:
        if refs[iref+1][-1] < ovlp*refs[iref][5]:
            if refs[iref][-1] > ovlp*refs[iref][5]:
                refs2.append([])
                start = True
            if nref2 == len(refs2):
                refs2.append([])
            refs2[nref2].append(refs[iref])
        else:
            if start:
                refs2[nref2].append(refs[iref])
                start = False
                nref2 += 1
            else:
                refs1.append(refs[iref])
        iref += 1
    if start:
        refs2[nref2].append(refs[iref])
    else:
        refs1.append(refs[iref])
    
    mem.write('%5d\n'%len(refs1))
    for ref in refs1:
        h,k,l = ref[:3]
        hkl = '%d %d %d'%(h,k,l)
        if hkl in refDict:
            del refDict[hkl]
        Fobs = np.sqrt(ref[6])
        mem.write('%5d%5d%5d%10.3f%10.3f%10.3f\n'%(h,k,l,Fobs*npcosd(ref[7]),Fobs*npsind(ref[7]),max(0.01*Fobs,0.1)))
    while True and nref2:
        if not len(refs2[-1]):
            del refs2[-1]
        else:
            break
    mem.write('%5d\n'%len(refs2))
    for iref2,ref2 in enumerate(refs2):
        mem.write('#%5d\n'%iref2)
        mem.write('%5d\n'%len(ref2))
        Gsum = 0.
        Msum = 0
        for ref in ref2:
            Gsum += ref[6]*ref[3]
            Msum += ref[3]
        G = np.sqrt(Gsum/Msum)
        h,k,l = ref2[0][:3]
        hkl = '%d %d %d'%(h,k,l)
        if hkl in refDict:
            del refDict[hkl]
        mem.write('%5d%5d%5d%10.3f%10.3f%5d\n'%(h,k,l,G,max(0.01*G,0.1),ref2[0][3]))
        for ref in ref2[1:]:
            h,k,l,m = ref[:4]
            mem.write('%5d%5d%5d%5d\n'%(h,k,l,m))
            hkl = '%d %d %d'%(h,k,l)
            if hkl in refDict:
                del refDict[hkl]
    if len(refDict):
        mem.write('%d\n'%len(refDict))
        for hkl in list(refDict.keys()):
            h,k,l = refDict[hkl][:3]
            mem.write('%5d%5d%5d\n'%(h,k,l))
    else:
        mem.write('0\n')
    mem.close()
    return True

def MEMupdateReflData(prfName,data,reflData):
    ''' Update reflection data with new Fosq, phase result from Dysnomia

    :param str prfName: phase.mem file name
    :param list reflData: GSAS-II reflection data
    '''
    
    generalData = data['General']
    Map = generalData['Map']
    Type = Map['Type']
    cell = generalData['Cell'][1:7]
    A = G2lat.cell2A(cell)
    reflDict = {}
    newRefs = []
    for iref,ref in enumerate(reflData):
        if ref[3] > 0:
            newRefs.append(ref)
            reflDict[hash('%5d%5d%5d'%(ref[0],ref[1],ref[2]))] = iref
    fbaName = os.path.splitext(prfName)[0]+'.fba'
    if os.path.isfile(fbaName):
        fba = open(fbaName,'r')
    else:
        return False
    fba.readline()
    Nref = int(fba.readline()[:-1])
    fbalines = fba.readlines()
    for line in fbalines[:Nref]:
        info = line.split()
        h = int(info[0])
        k = int(info[1])
        l = int(info[2])
        FoR = float(info[3])
        FoI = float(info[4])
        Fosq = FoR**2+FoI**2
        phase = npatan2d(FoI,FoR)
        try:
            refId = reflDict[hash('%5d%5d%5d'%(h,k,l))]
        except KeyError:    #added reflections at end skipped
            d = float(1/np.sqrt(G2lat.calc_rDsq([h,k,l],A)))
            if 'T' in Type:
                newRefs.append([h,k,l,-1,d,0.,0.01,1.0,Fosq,Fosq,phase,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
            else:
                newRefs.append([h,k,l,-1,d,0.,0.01,1.0,Fosq,Fosq,phase,1.0,1.0,1.0,1.0])
            continue
        newRefs[refId][8] = Fosq
        newRefs[refId][10] = phase
    newRefs = np.array(newRefs)
    return True,newRefs
    
#### testing data
NeedTestData = True
def TestData():
    'needs a doc string'
#    global NeedTestData
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
        'U':1.163,'V':-0.605,'W':0.093,'X':0.0,'Y':2.183,'Z':0.0,'SH/L':0.002,
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
        'U':22.75,'V':-17.596,'W':10.594,'X':1.577,'Y':5.778,'Z':0.0,'SH/L':0.002,
        'Back0':36.897,'Back1':-0.508,'Back2':.006,
        'Lam1':1.540500,'Lam2':1.544300,'I(L2)/I(L1)':0.5,
        }
    global parmDict2
    parmDict2 = {
        'pos0':5.7,'int0':1000.0,'sig0':0.5,'gam0':0.5,
        'U':2.,'V':-2.,'W':5.,'X':0.5,'Y':0.5,'Z':0.0,'SH/L':0.02,
        'Back0':5.,'Back1':-0.02,'Back2':.004,
#        'Lam1':1.540500,'Lam2':1.544300,'I(L2)/I(L1)':0.5,
        }
    global varyList
    varyList = []

def test0():
    if NeedTestData: TestData()
    gplot = plotter.add('FCJ-Voigt, 11BM').gca()
    gplot.plot(xdata,getBackground('',parmDict0,bakType,'PXC',xdata)[0])   
    gplot.plot(xdata,getPeakProfile(parmDict0,xdata,varyList,bakType))
    fplot = plotter.add('FCJ-Voigt, Ka1+2').gca()
    fplot.plot(xdata,getBackground('',parmDict1,bakType,'PXC',xdata)[0])   
    fplot.plot(xdata,getPeakProfile(parmDict1,xdata,varyList,bakType))
   
def test1():
    if NeedTestData: TestData()
    time0 = time.time()
    for i in range(100):
        getPeakProfile(parmDict1,xdata,varyList,bakType)
    G2fil.G2Print ('100+6*Ka1-2 peaks=1200 peaks %.2f'%time.time()-time0)
    
def test2(name,delt):
    if NeedTestData: TestData()
    varyList = [name,]
    xdata = np.linspace(5.6,5.8,400)
    hplot = plotter.add('derivatives test for '+name).gca()
    hplot.plot(xdata,getPeakProfileDerv(parmDict2,xdata,varyList,bakType)[0])
    y0 = getPeakProfile(parmDict2,xdata,varyList,bakType)
    parmDict2[name] += delt
    y1 = getPeakProfile(parmDict2,xdata,varyList,bakType)
    hplot.plot(xdata,(y1-y0)/delt,'r+')
    
def test3(name,delt):
    if NeedTestData: TestData()
    names = ['pos','sig','gam','shl']
    idx = names.index(name)
    myDict = {'pos':parmDict2['pos0'],'sig':parmDict2['sig0'],'gam':parmDict2['gam0'],'shl':parmDict2['SH/L']}
    xdata = np.linspace(5.6,5.8,800)
    dx = xdata[1]-xdata[0]
    hplot = plotter.add('derivatives test for '+name).gca()
    hplot.plot(xdata,100.*dx*getdFCJVoigt3(myDict['pos'],myDict['sig'],myDict['gam'],myDict['shl'],xdata)[idx+1])
    y0 = getFCJVoigt3(myDict['pos'],myDict['sig'],myDict['gam'],myDict['shl'],xdata)
    myDict[name] += delt
    y1 = getFCJVoigt3(myDict['pos'],myDict['sig'],myDict['gam'],myDict['shl'],xdata)
    hplot.plot(xdata,(y1-y0)/delt,'r+')

if __name__ == '__main__':
    import GSASIItestplot as plot
    global plotter
    plotter = plot.PlotNotebook()
#    test0()
#    for name in ['int0','pos0','sig0','gam0','U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)']:
    for name,shft in [['int0',0.1],['pos0',0.0001],['sig0',0.01],['gam0',0.00001],
        ['U',0.1],['V',0.01],['W',0.01],['X',0.0001],['Y',0.0001],['Z',0.0001],['SH/L',0.00005]]:
        test2(name,shft)
    for name,shft in [['pos',0.0001],['sig',0.01],['gam',0.0001],['shl',0.00005]]:
        test3(name,shft)
    G2fil.G2Print ("OK")
    plotter.StartEventLoop()
