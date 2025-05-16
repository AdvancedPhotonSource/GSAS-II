# -*- coding: utf-8 -*-
'''
Classes and routines defined in :mod:`GSASIIpwd` follow.
'''

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
import scipy.signal as signal

from . import GSASIIpath
GSASIIpath.SetBinaryPath()

filversion = str(GSASIIpath.GetVersionNumber())
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIElem as G2elem
from . import GSASIImath as G2mth
from . import GSASIIfiles as G2fil
try:
    if GSASIIpath.binaryPath:
        import  pypowder as pyd
    else:
        from . import pypowder as pyd
except ImportError:
    print ('pypowder is not available - profile calcs. not allowed')
try:
    if GSASIIpath.binaryPath:
        import pydiffax as pyx
    else:
        from . import pydiffax as pyx
except ImportError:
    print ('pydiffax is not available for this platform')

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
npq2T = lambda Q,wave: 2.0*npasind(0.25*Q*wave/np.pi)
ateln2 = 8.0*np.log(2.0)
sateln2 = np.sqrt(ateln2)
nxs = np.newaxis
is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)

#### Powder utilities ################################################################################
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

#### GSASII pwdr & pdf calculation routines ################################################################################
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
    :param Azm: azimuthal angle e.g. 0.0 in plane of polarization - can be numpy array
    :param Tth: 2-theta scattering angle - can be numpy array
      which (if either) of these is "right"?
    :return: (pola, dpdPola) - both 2-d arrays
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
        K = (1.-ObCoeff)/(1.0-np.exp(np.log(ObCoeff)/npcosd(Tth)))
        return K
    else:
        return 1.0

def Ruland(RulCoff,wave,Q,Compton):
    'needs a doc string'
    C = 2.9978e8
    D = 1.5e-3
    hmc = 0.024262734687    #Compton wavelength in A
    sinth2 = (Q*wave/(4.0*np.pi))**2
    dlam = (wave**2)*Compton*Q/C
    dlam_c = 2.0*hmc*sinth2-D*wave**2
    return 1.0/((1.0+dlam/RulCoff)*(1.0+(np.pi*dlam_c/(dlam+RulCoff))**2))

def KleinNishina(wave,Q):
    hmc = 0.024262734687    #Compton wavelength in A
    TTh = npq2T(Q,wave)
    P = 1./(1.+(1.-npcosd(TTh)*(hmc/wave)))
    KN = (P**3-(P*npsind(TTh))**2+P)/(1.+npcosd(TTh)**2)
    return KN

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
        try:   # fails if background differs in number of points
            IofQ[1][1] += xydata['Sample Bkg.'][1][1][Ibeg:Ifin]*data['Sample Bkg.']['Mult']
        except ValueError:
            print("Interpolating Sample background since points don't match")
            interpF = si.interp1d(xydata['Sample Bkg.'][1][0],xydata['Sample Bkg.'][1][1],
                                  fill_value='extrapolate')
            IofQ[1][1] += interpF(IofQ[1][0]) * data['Sample Bkg.']['Mult']
    if data['Container']['Name']:
        xycontainer = xydata['Container'][1][1]*data['Container']['Mult']
        if data['Container Bkg.']['Name']:
            try:
                xycontainer += xydata['Container Bkg.'][1][1][Ibeg:Ifin]*data['Container Bkg.']['Mult']
            except ValueError:
                print('Number of points do not agree between Container and Container Bkg.')
                return
        try:   # fails if background differs in number of points
            IofQ[1][1] += xycontainer[Ibeg:Ifin]
        except ValueError:
            print("Interpolating Container background since points don't match")
            interpF = si.interp1d(xydata['Container'][1][0],xycontainer,fill_value='extrapolate')
            IofQ[1][1] += interpF(IofQ[1][0])

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
        if data['DetType'] == 'Area detector':
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
    # improves look of F(Q) but no impact on G(R)
    # bBut,aBut = signal.butter(8,.5,"lowpass")
    # IofQ[1][1] = signal.filtfilt(bBut,aBut,IofQ[1][1])
    XY = IofQ[1]
    #convert to Q
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
#        CF *= KleinNishina(wave,Q)
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
    GR = -(2./np.pi)*dq*np.imag(fft.fft(xydata['FofQ'][1][1],mul*nR)[:nR])*data.get('GR Scale',1.0)
#    GR = -dq*np.imag(fft.fft(xydata['FofQ'][1][1],mul*nR)[:nR])*data.get('GR Scale',1.0)
    xydata['GofR'][1][1] = GR
    numbDen = 0.
    if 'ElList' in data:
        numbDen = GetNumDensity(data['ElList'],data['Form Vol'])
        gr = GR/(4.*np.pi*numbDen*R)+1.
#        gr = GR/(np.pi*R) ##mising numberdensity
        xydata['gofr'][1][1] = gr
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
    auxPlot = []
    if 'C' in inst['Type'][0] or 'B' in inst['Type'][0]:
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
    elif RDFcontrols['UseObsCalc'] == 'auto-back':
        auto = autoBkgCalc(background[1],pwddata[1])
        Qdata = si.griddata(powQ,auto-pwddata[4],Qpoints,method=RDFcontrols['Smooth'],fill_value=0.)
    Qdata *= np.sin((Qpoints-minQ)*piDQ)/piDQ
    Qdata *= 0.5*np.sqrt(Qpoints)       #Qbin normalization
    dq = Qpoints[1]-Qpoints[0]
    nR = len(Qdata)
    R = 0.5*np.pi*np.linspace(0,nR,nR)/(4.*maxQ)
    iFin = np.searchsorted(R,RDFcontrols['maxR'])+1
    bBut,aBut = signal.butter(4,0.01)
    Qsmooth = signal.filtfilt(bBut,aBut,Qdata)
#    auxPlot.append([Qpoints,Qdata,'interpolate:'+RDFcontrols['Smooth']])
#    auxPlot.append([Qpoints,Qsmooth,'interpolate:'+RDFcontrols['Smooth']])
    DofR = dq*np.imag(fft.fft(Qsmooth,16*nR)[:nR])
    auxPlot.append([R[:iFin],DofR[:iFin],'D(R) for '+RDFcontrols['UseObsCalc']])
    return auxPlot

# PDF optimization =============================================================
def OptimizePDF(data,xydata,limits,inst,showFit=True,maxCycles=25):
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
        res = opt.minimize(Min,xstart,bounds=([0.01,1.],[1.2*bakMul,0.8*bakMul]),
                    method='L-BFGS-B',options={'maxiter':maxCycles},tol=0.001)
    else:
        res = opt.minimize(Min,xstart,bounds=([0.,None],[0,1],[0.01,1.]),
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
            G2fil.G2Print('  end:   Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f} RMS:{:.4f}) *** {} ***\n'.format(
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
            Data['Flat Bkg'] = BkgMax*(2.*F-1.)
            Data['BackRatio'] = B
        Data['Ruland'] = R
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
                return [max(data['Ruland'],.05),data['Sample']['Mult']]
        try:
            F = 0.5+0.5*data['Flat Bkg']/BkgMax
        except:
            F = 0
        return [F,data['BackRatio'],max(data['Ruland'],.05)]
    def SetFinalVals(arg):
        '''Set the 'Flat Bkg', 'BackRatio' & 'Ruland' values from the
        scaled, refined values and plot corrected region of G(r)
        '''
        if data['Sample Bkg.'].get('Refine',False):
            R,S = arg
            data['Sample Bkg.']['Mult'] = S
        else:
            F,B,R = arg
            data['Flat Bkg'] = BkgMax*(2.*F-1.)
            data['BackRatio'] = B
        data['Ruland'] = R
        CalcPDF(data,inst,limits,xydata)
    EvalLowPDF(GetCurrentVals())
    BkgMax = max(xydata['IofQ'][1][1])/50.
    return EvalLowPDF,GetCurrentVals,SetFinalVals

#### GSASII convolution peak fitting routines: Finger, Cox & Jephcoat model
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
    '''
    Normal distribution

The location (loc) keyword specifies the mean.
The scale (scale) keyword specifies the standard deviation.

normal.pdf(x) = exp(-x**2/2)/sqrt(2*pi)

    '''

    def pdf(self,x,*args,**kwds):
        loc,scale=kwds['loc'],kwds['scale']
        x = (x-loc)/scale
        return np.exp(-x**2/2.0) * _norm_pdf_C / scale

norm = norm_gen(name='norm')

## Cauchy

# median = loc

class cauchy_gen(st.rv_continuous):
    '''
Cauchy distribution

cauchy.pdf(x) = 1/(pi*(1+x**2))

This is the t distribution with one degree of freedom.
    '''

    def pdf(self,x,*args,**kwds):
        loc,scale=kwds['loc'],kwds['scale']
        x = (x-loc)/scale
        return 1.0/np.pi/(1.0+x*x) / scale

cauchy = cauchy_gen(name='cauchy')


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

#### GSASII peak fitting routine: Finger, Cox & Jephcoat model

def getWidthsCW(pos,sig,gam,shl):
    '''Compute the peak widths used for computing the range of a peak
    for constant wavelength data. On low-angle side, 50 FWHM are used,
    on high-angle side 75 are used, low angle side extended for axial divergence
    (for peaks above 90 deg, these are reversed.)

    :param pos: peak position; 2-theta in degrees
    :param sig: Gaussian peak variance in centideg^2
    :param gam: Lorentzian peak width in centidegrees
    :param shl: axial divergence parameter (S+H)/L

    :returns: widths; [Gaussian sigma, Lorentzian gamma] in degrees, and
        low angle, high angle ends of peak; 50 FWHM & 75 FWHM from position
        reversed for 2-theta > 90 deg.
    '''
    widths = [np.sqrt(sig)/100.,gam/100.]
    fwhm = max(2.355*widths[0]+widths[1],0.0001)
    fmin = 50.*(fwhm+shl*abs(npcosd(pos)))
    fmax = 75.0*fwhm
    if pos > 90.:
        fmin,fmax = [fmax,fmin]
    return widths,fmin,fmax

def getWidthsED(pos,sig,gam):
    '''Compute the peak widths used for computing the range of a peak
    for energy dispersive data. On low-energy side, 20 FWHM are used,
    on high-energy side 20 are used

    :param pos: peak position; energy in keV (not used)
    :param sig: Gaussian peak variance in keV^2
    :param gam: Lorentzian peak width in keV

    :returns: widths; [Gaussian sigma, Lorentzian gamma] in keV, and
        low angle, high angle ends of peak; 5 FWHM & 5 FWHM from position
    '''
    widths = [np.sqrt(sig),gam]
    fwhm = 2.355*widths[0]+widths[1]
    fmin = 5.*fwhm
    fmax = 5.*fwhm
    return widths,fmin,fmax

def getWidthsTOF(pos,alp,bet,sig,gam):
    '''Compute the peak widths used for computing the range of a peak
    for TOF data. 50 FWHM are used on both sides each
    extended by exponential coeff.

    param pos: peak position; TOF in musec (not used)
    param alp,bet: TOF peak exponential rise & decay parameters
    param sig: Gaussian peak variance in musec^2
    param gam: Lorentzian peak width in musec

    returns: widths; [Gaussian sigma, Lornetzian gamma] in musec
    returns: low TOF, high TOF ends of peak; 50FWHM from position
    '''
    widths = [np.sqrt(sig),gam]
    fwhm = 2.355*widths[0]+2.*widths[1]
    fmin = 50.*fwhm*(1.+1./alp)
    fmax = 50.*fwhm*(1.+1./bet)
    return widths,fmin,fmax

def getWidthsCWA(pos,alp,bet,sig,gam,shl):
    '''Compute the peak widths used for computing the range of a peak
    for constant wavelength data with axial divergence. 50 & 75 FWHM are used on
    each side each extended by exponential coeff.

    :param pos: peak position; 2-theta in degrees
    :param alp,bet: TOF peak exponential rise & decay parameters
    :param sig: Gaussian peak variance in centideg^2
    :param gam: Lorentzian peak width in centidegrees
    :param shl: axial divergence parameter (S+H)/L

    :returns: widths; [Gaussian sigma, Lorentzian gamma] in degrees, and
        low angle, high angle ends of peak; 50 FWHM & 75 FWHM from position
        reversed for 2-theta > 90 deg.
    '''
    widths = [np.sqrt(sig)/100.,gam/100.]
    fwhm = 2.355*widths[0]+2.*widths[1]
    if pos < 90.:
        fmin = 50.0*(fwhm*(1.+1./alp)+shl*abs(npcosd(pos)))
        fmax = 75.0*fwhm*(1.+1./bet)
    else:
        fmin = 75.0*fwhm*(1.+1./alp)
        fmax = 50.0*(fwhm*(1.+1./bet)+shl*abs(npcosd(pos)))
    return widths,fmin,fmax

def getWidthsCWB(pos,alp,bet,sig,gam):
    '''Compute the peak widths used for computing the range of a peak
    for constant wavelength data without axial divergence. 50 FWHM are used on
    both sides each extended by exponential coeff.

    :param pos: peak position; 2-theta in degrees
    :param alp,bet: TOF peak exponential rise & decay parameters
    :param sig: Gaussian peak variance in centideg^2
    :param gam: Lorentzian peak width in centidegrees

    returns: widths; [Gaussian sigma, Lornetzian gamma] in degrees
    returns: low angle, high angle ends of peak; 50FWHM from position
    '''
    widths = [np.sqrt(sig)/100.,gam/100.]
    fwhm = 2.355*widths[0]+2.*widths[1]
    fmin = 50.*fwhm*(1.+1./alp)
    fmax = 50.*fwhm*(1.+1./bet)
    return widths,fmin,fmax

def getFWHM(pos,Inst,N=1):
    '''Compute total FWHM from Thompson, Cox & Hastings (1987) , J. Appl. Cryst. 20, 79-83
    via getgamFW(g,s).

    :param pos: float peak position in deg 2-theta or tof in musec
    :param Inst: dict instrument parameters
    :param N: int Inst index (0 for input, 1 for fitted)

    :returns float: total FWHM of pseudoVoigt in deg or musec
    '''

    sig = lambda Th,U,V,W: np.sqrt(max(0.001,U*tand(Th)**2+V*tand(Th)+W))
    sigED = lambda E,A,B,C: np.sqrt(max(0.001,A*E**2+B*E+C))
    sigTOF = lambda dsp,S0,S1,S2,Sq: np.sqrt(S0+S1*dsp**2+S2*dsp**4+Sq*dsp)
    gam = lambda Th,X,Y,Z: Z+X/cosd(Th)+Y*tand(Th)
    gamED = lambda E,X,Y,Z: max(0.001,X*E**2+Y*E+Z)
    gamTOF = lambda dsp,X,Y,Z: Z+X*dsp+Y*dsp**2
    alpTOF = lambda dsp,alp: alp/dsp
    betTOF = lambda dsp,bet0,bet1,betq: bet0+bet1/dsp**4+betq/dsp**2
    alpPinkX = lambda pos,alp0,alp1: alp0+alp1*nptand(pos/2.)
    betPinkX = lambda pos,bet0,bet1: bet0+bet1*nptand(pos/2.)
    alpPinkN = lambda pos,alp0,alp1: alp0+alp1*npsind(pos/2.)
    betPinkN = lambda pos,bet0,bet1: bet0+bet1*npsind(pos/2.)
    if 'LF' in Inst['Type'][0]:
        return 3
    elif 'T' in Inst['Type'][0]:
        dsp = pos/Inst['difC'][N]
        alp = alpTOF(dsp,Inst['alpha'][N])
        bet = betTOF(dsp,Inst['beta-0'][1],Inst['beta-1'][N],Inst['beta-q'][N])
        s = sigTOF(dsp,Inst['sig-0'][N],Inst['sig-1'][N],Inst['sig-2'][N],Inst['sig-q'][N])
        g = gamTOF(dsp,Inst['X'][N],Inst['Y'][N],Inst['Z'][N])
        return getgamFW(g,s)+np.log(2.0)*(alp+bet)/(alp*bet)
    elif 'C' in Inst['Type'][0]:
        s = sig(pos/2.,Inst['U'][N],Inst['V'][N],Inst['W'][N])
        g = gam(pos/2.,Inst['X'][N],Inst['Y'][N],Inst['Z'][N])
        return getgamFW(g,s)/100.  #returns FWHM in deg
    elif 'E' in Inst['Type'][0]:
        s = sigED(pos,Inst['A'][N],Inst['B'][N],Inst['C'][N])
        g = gamED(pos,Inst['X'][N],Inst['Y'][N],Inst['Z'][N])
        return getgamFW(g,s)
    else:   #'B'
        if 'X' in Inst['Type'][0]:
            alp = alpPinkX(pos,Inst['alpha-0'][N],Inst['alpha-1'][N])
            bet = betPinkX(pos,Inst['beta-0'][N],Inst['beta-1'][N])
        else:
            alp = alpPinkN(pos,Inst['alpha-0'][N],Inst['alpha-1'][N])
            bet = betPinkN(pos,Inst['beta-0'][N],Inst['beta-1'][N])
        s = sig(pos/2.,Inst['U'][N],Inst['V'][N],Inst['W'][N])
        g = gam(pos/2.,Inst['X'][N],Inst['Y'][N],Inst['Z'][N])
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

def getBackground(pfx,parmDict,bakType,dataType,xdata,fixback=None):
    '''Computes the background based on parameters that may be taken from
    a gpx file or the data tree.

    :param str pfx: histogram prefix (:h:)
    :param dict parmDict: Refinement parameter values
    :param str bakType: defines background function to be used. Should be
      one of these: 'chebyschev', 'cosine', 'chebyschev-1',
      'Q^2 power series', 'Q^-2 power series', 'lin interpolate',
      'inv interpolate', 'log interpolate'
    :param str dataType: Code to indicate histogram type (PXC, PNC, PNT,...)
    :param MaskedArray xdata: independent variable, 2theta (deg*100) or
      TOF (microsec?)
    :param numpy.array fixback: Array of fixed background points (length
      matching xdata) or None

    :returns: yb,sumBK where yp is an array of background values (length
      matching xdata) and sumBK is a list with three values. The sumBK[0] is
      the sum of all yb values, sumBK[1] is the sum of Debye background terms
      and sumBK[2] is the sum of background peaks.
    '''
    if 'T' in dataType:
        q = 2.*np.pi*parmDict[pfx+'difC']/xdata
    elif 'E' in dataType:
        const = 4.*np.pi*npsind(parmDict[pfx+'2-theta']/2.0)
        q = const*xdata
    else:
        wave = parmDict.get(pfx+'Lam',parmDict.get(pfx+'Lam1',1.0))
        q = npT2q(xdata,wave)
    yb = np.zeros_like(xdata)
    nBak = 0
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
    if pfx+'difC' in parmDict or 'E' in dataType:
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
            pkS = max(parmDict[pfx+'BkPksig;'+str(iD)],0.01)
            pkG = max(parmDict[pfx+'BkPkgam;'+str(iD)],0.1)
            if 'C' in dataType:
                shl = parmDict[pfx+'SH/L']
                Wd,fmin,fmax = getWidthsCW(pkP,pkS,pkG,shl)
            elif 'B' in dataType:
                sinPos = npsind(pkP/2.0)
                alp = max(0.1,parmDict[pfx+'alpha-0']+parmDict[pfx+'alpha-1']*sinPos)
                bet = max(0.001,parmDict[pfx+'beta-0']+parmDict[pfx+'beta-1']*sinPos)
                Wd,fmin,fmax = getWidthsCWB(pkP,alp,bet,pkS,pkG)
            elif 'A'' in dataType':
                shl = parmDict[pfx+'SH/L']
                sinPos = npsind(pkP/2.0)
                alp = max(0.1,parmDict[pfx+'alpha-0']+parmDict[pfx+'alpha-1']*sinPos)
                bet = max(0.001,parmDict[pfx+'beta-0']+parmDict[pfx+'beta-1']*sinPos)
                Wd,fmin,fmax = getWidthsCWA(pkP,alp,bet,pkS,pkG,shl)
            elif 'E' in dataType:
                Wd,fmin,fmax = getWidthsED(pkP,pkS)
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
                ybi = pkI*getFCJVoigt3(pkP,pkS,pkG,shl,xdata[iBeg:iFin])[0]
            elif 'B' in dataType:
                ybi = pkI*getEpsVoigt(pkP,alp,bet,pkS/1.e4,pkG/100.,xdata[iBeg:iFin])[0]/100.
            elif 'A' in dataType:
                ybi = pkI*getExpFCJVoigt3(pkP,alp,bet,pkS,pkG,shl,xdata[iBeg:iFin])[0]
            elif 'T' in dataType:
                ybi = pkI*getEpsVoigt(pkP,1.,1.,pkS,pkG,xdata[iBeg:iFin])[0]
            elif 'E' in dataType:
                ybi = pkI*getPsVoigt(pkP,pkS*10.**4,pkG*100.,xdata[iBeg:iFin])[0]
            else:
                raise Exception('dataType of {:} should not happen!'.format(dataType))
            yb[iBeg:iFin] += ybi
            sumBk[2] += np.sum(ybi)
            iD += 1
        except KeyError:
            break
        except ValueError:
            G2fil.G2Print ('**** WARNING - backround peak '+str(iD)+' sigma is negative; fix & try again ****')
            break
    if fixback is not None:
        yb += parmDict.get(pfx+'BF mult',1.0)*fixback
        sumBk[0] = sum(yb)
    return yb,sumBk

def getBackgroundDerv(hfx,parmDict,bakType,dataType,xdata,fixback=None):
    '''Computes the derivatives of the background
    Parameters passed to this may be pulled from gpx file or data tree.
    See :func:`getBackground` for parameter definitions.

    :returns: dydb,dyddb,dydpk,dydfb where the first three are 2-D arrays
      of derivatives with respect to the background terms, the Debye terms and
      the background peak terms vs. the points in the diffracton pattern. The
      final 1D array is the derivative with respect to the fixed-background
      multiplier (= the fixed background values).
    '''
    if 'T' in dataType:
        q = 2.*np.pi*parmDict[hfx+'difC']/xdata
    elif 'E' in dataType:
        const = 4.*np.pi*npsind(parmDict[hfx+'2-theta']/2.0)
        q = const*xdata
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
    dydfb = []

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
            pkS = max(parmDict[hfx+'BkPksig;'+str(iD)],0.01)
            pkG = max(parmDict[hfx+'BkPkgam;'+str(iD)],0.1)
            if 'C' in dataType:
                shl = parmDict[hfx+'SH/L']
                Wd,fmin,fmax = getWidthsCW(pkP,pkS,pkG,shl)
            elif 'B' in dataType:
                sinPos = npsind(pkP/2.0)
                alp = max(0.1,parmDict[hfx+'alpha-0']+parmDict[hfx+'alpha-1']*sinPos)
                bet = max(0.001,parmDict[hfx+'beta-0']+parmDict[hfx+'beta-1']*sinPos)
                Wd,fmin,fmax = getWidthsCWB(pkP,alp,bet,pkS,pkG)
            elif 'A'' in dataType':
                shl = parmDict[hfx+'SH/L']
                sinPos = npsind(pkP/2.0)
                alp = max(0.1,parmDict[hfx+'alpha-0']+parmDict[hfx+'alpha-1']*sinPos)
                bet = max(0.001,parmDict[hfx+'beta-0']+parmDict[hfx+'beta-1']*sinPos)
                Wd,fmin,fmax = getWidthsCWA(pkP,alp,bet,pkS,pkG,shl)
            elif 'E' in dataType:
                Wd,fmin,fmax = getWidthsED(pkP,pkS)
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
                Df,dFdp,dFds,dFdg,x = getdFCJVoigt3(pkP,pkS,pkG,shl,xdata[iBeg:iFin])
            elif 'B' in dataType:
                Df,dFdp,x,x,dFds,dFdg = getdEpsVoigt(pkP,alp,bet,pkS/1.e4,pkG/100.,xdata[iBeg:iFin])
                dFdp /= 100.
                dFds /= 1.e6
                dFdg /= 1.e4
                Df /= 100.
            elif 'A' in dataType:
                Df,dFdp,x,x,dFds,dFdg,x = getdExpFCJVoigt3(pkP,alp,bet,pkS,pkG,shl,xdata[iBeg:iFin])
            elif 'E' in dataType:
                Df,dFdp,dFds,dFdg = getdPsVoigt(pkP,pkS*10.**4,pkG*100.,xdata[iBeg:iFin])
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
    # fixed background from file
    if  fixback is not None:
        dydfb = fixback
    return dydb,dyddb,dydpk,dydfb

#### Using old gsas fortran routines for powder peak shapes & derivatives
def getFCJVoigt3(pos,sig,gam,shl,xdata):
    '''Compute the Finger-Cox-Jepcoat modified Pseudo-Voigt function for a
    CW powder peak in external Fortran routine

    param pos: peak position in degrees
    param sig: Gaussian variance in centideg^2
    param gam: Lorentzian width in centideg
    param shl: axial divergence parameter (S+H)/L
    param xdata: array; profile points for peak to be calculated; bounded by 20FWHM to 50FWHM (or vv)

    returns: array: calculated peak function at each xdata
    returns: integral of peak; nominally = 1.0
    '''
    if len(xdata):
        cw = np.diff(xdata)
        cw = np.append(cw,cw[-1])
        Df = pyd.pypsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)
        return Df,np.sum(100.*Df*cw)
    else:
        return 0.,1.

def getdFCJVoigt3(pos,sig,gam,shl,xdata):
    '''Compute analytic derivatives the Finger-Cox-Jepcoat modified Pseudo-Voigt
    function for a CW powder peak

    param pos: peak position in degrees
    param sig: Gaussian variance in centideg^2
    param gam: Lorentzian width in centideg
    param shl: axial divergence parameter (S+H)/L
    param xdata: array; profile points for peak to be calculated; bounded by 20FWHM to 50FWHM (or vv)

    returns: arrays: function and derivatives of pos, sig, gam, & shl
    '''
    Df,dFdp,dFds,dFdg,dFdsh = pyd.pydpsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)
    return Df,dFdp,dFds,dFdg,dFdsh

def getExpFCJVoigt3(pos,alp,bet,sig,gam,shl,xdata):
    '''Compute the Finger-Cox-Jepcoat modified double exponential Pseudo-Voigt
    convolution function for a CW powder peak in external Fortran routine

    param pos: peak position in degrees
    param sig: Gaussian variance in centideg^2
    param alp: Rise exponential coefficient
    param bet: Decay exponential coefficient
    param gam: Lorentzian width in centideg
    param shl: axial divergence parameter (S+H)/L
    param xdata: array; profile points for peak to be calculated; bounded by 20FWHM to 50FWHM (or vv)

    returns: array: calculated peak function at each xdata
    returns: integral of peak; nominally = 1.0
    '''
    if len(xdata):
        cw = np.diff(xdata)
        cw = np.append(cw,cw[-1])
        Df = pyd.pypsvfcjexpo(len(xdata),xdata-pos,pos,alp,bet,sig,gam,shl)
        return Df,np.sum(100.*Df*cw)
    else:
        return 0.,1.

def getdExpFCJVoigt3(pos,alp,bet,sig,gam,shl,xdata):
    '''Compute analytic derivatives the Finger-Cox-Jepcoat modified double
    exponential Pseudo-Voigt convolution function for a CW powder peak

    param pos: peak position in degrees
    param alp: Rise exponential coefficient
    param bet: Decay exponential coefficient
    param sig: Gaussian variance in centideg^2
    param gam: Lorentzian width in centideg
    param shl: axial divergence parameter (S+H)/L
    param xdata: array; profile points for peak to be calculated; bounded by 20FWHM to 50FWHM (or vv)

    returns: arrays: function and derivatives of pos, alp, bet, sig, gam, & shl
    '''
    Df,dFdp,dFda,dFdb,dFds,dFdg,dFdsh = pyd.pydpsvfcjexpo(len(xdata),xdata-pos,pos,alp,bet,sig,gam,shl)
    return Df,dFdp,dFda,dFdb,dFds,dFdg,dFdsh

def getPsVoigt(pos,sig,gam,xdata):
    '''Compute the simple Pseudo-Voigt function for a
    small angle Bragg peak in external Fortran routine

    param pos: peak position in degrees
    param sig: Gaussian variance in centideg^2
    param gam: Lorentzian width in centideg
    param xdata: array; profile points for peak to be calculated

    returns: array: calculated peak function at each xdata
    returns: integral of peak; nominally = 1.0
    '''

    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])
    Df = pyd.pypsvoigt(len(xdata),xdata-pos,sig,gam)
    return Df,np.sum(100.*Df*cw)

def getdPsVoigt(pos,sig,gam,xdata):
    '''Compute the simple Pseudo-Voigt function derivatives for a
    small angle Bragg peak peak in external Fortran routine

    param pos: peak position in degrees
    param sig: Gaussian variance in centideg^2
    param gam: Lorentzian width in centideg
    param xdata: array; profile points for peak to be calculated

    returns: arrays: function and derivatives of pos, sig & gam
    NB: the pos derivative has the opposite sign compared to that in other profile functions
    '''

    Df,dFdp,dFds,dFdg = pyd.pydpsvoigt(len(xdata),xdata-pos,sig,gam)
    return Df,dFdp,dFds,dFdg

def getEpsVoigt(pos,alp,bet,sig,gam,xdata):
    '''Compute the double exponential Pseudo-Voigt convolution function for a
    neutron TOF & CW pink peak in external Fortran routine
    '''

    cw = np.diff(xdata)
    cw = np.append(cw,cw[-1])
    Df = pyd.pyepsvoigt(len(xdata),xdata-pos,alp,bet,sig,gam)
    return Df,np.sum(Df*cw)

def getdEpsVoigt(pos,alp,bet,sig,gam,xdata):
    '''Compute the double exponential Pseudo-Voigt convolution function derivatives for a
    neutron TOF & CW pink peak in external Fortran routine
    '''

    Df,dFdp,dFda,dFdb,dFds,dFdg = pyd.pydepsvoigt(len(xdata),xdata-pos,alp,bet,sig,gam)
    return Df,dFdp,dFda,dFdb,dFds,dFdg

def ellipseSize(H,Sij,GB):
    '''Implements r=1/sqrt(sum((1/S)*(q.v)^2) per note from Alexander Brady
    '''

    HX = np.inner(H.T,GB)
    lenHX = np.sqrt(np.sum(HX**2))
    Esize,Rsize = nl.eigh(G2lat.U6toUij(Sij))
    R = np.inner(HX/lenHX,Rsize)**2*Esize         #want column length for hkl in crystal
    lenR = 1./np.sqrt(np.sum(R))
    return lenR

def ellipseSizeDerv(H,Sij,GB):
    '''Implements r=1/sqrt(sum((1/S)*(q.v)^2) derivative per note from Alexander Brady
    '''

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

peakInstPrmMode = True
'''Determines the mode used for peak fitting. When peakInstPrmMode=True peak
width parameters are computed from the instrument parameters (UVW,... or
alpha,... etc) unless the individual parameter is refined. This allows the
instrument parameters to be refined. When peakInstPrmMode=False, the instrument
parameters are not used and cannot be refined.
The default is peakFitMode=True. This is changed only in
:func:`setPeakInstPrmMode`, which is called from :mod:`GSASIIscriptable`
or GSASIIphsGUI.OnSetPeakWidMode ('Gen unvaried widths' menu item).
'''

def setPeakInstPrmMode(normal=True):
    '''Determines the mode used for peak fitting. If normal=True (default)
    peak width parameters are computed from the instrument parameters
    unless the individual parameter is refined. If normal=False,
    peak widths are used as supplied for each peak.

    Note that normal=True unless this routine is called. Also,
    instrument parameters can only be refined with normal=True.

    :param bool normal: setting to apply to global variable
      :data:`peakInstPrmMode`
    '''
    global peakInstPrmMode
    peakInstPrmMode = normal

def getPeakProfile(dataType,parmDict,xdata,fixback,varyList,bakType):
    '''Computes the profiles from multiple peaks for individual peak fitting
    for powder patterns.
    NB: not used for Rietveld refinement
    '''
    yb = getBackground('',parmDict,bakType,dataType,xdata,fixback)[0]
    yc = np.zeros_like(yb)
    if 'LF' in dataType:
        if 'Lam1' in parmDict.keys():
            lam = parmDict['Lam1']
            lam2 = parmDict['Lam2']
            Ka2 = True
            lamRatio = 360*(lam2-lam)/(np.pi*lam)
            kRatio = parmDict['I(L2)/I(L1)']
        else:
            lam = parmDict['Lam']
            Ka2 = False
        shol = 0
        # loop over peaks
        iPeak = -1
        try:
            ncells =  parmDict['ncell']
            clat =  parmDict['clat']
        except KeyError: # no Laue info must be bkg fit
            print('Laue Fit: no params, assuming bkg fit')
            return yb
        while True:
            iPeak += 1
            try:
                #Qcen = 2 * np.pi * lam * (iPeak+1) / parmDict['clat']
                l = parmDict['l'+str(iPeak)]
                pos = 360 * np.arcsin(0.5 * lam * l / parmDict['clat']) / np.pi
                parmDict['pos'+str(iPeak)] = pos
                #tth = (pos-parmDict['Zero'])
                intens = parmDict['int'+str(iPeak)]
                dampM =  parmDict['dampM'+str(iPeak)]
                dampP =  parmDict['dampP'+str(iPeak)]
                sig =  parmDict['sig'+str(iPeak)]
                gam =  parmDict['gam'+str(iPeak)]
                fmin = parmDict.get('fitRange',8.0) # width for peak computation: defaults to 8 deg.
                fmin = min(0.9*abs(xdata[-1] - xdata[0]),fmin) # unless the data range is smaller
                fitPowerM = parmDict.get('fitPowerM',2.0)
                fitPowerP = parmDict.get('fitPowerP',2.0)
                iBeg = np.searchsorted(xdata,pos-fmin/2)
                iFin = np.searchsorted(xdata,pos+fmin/2)
                if not iBeg+iFin:       # skip peak below low limit
                    continue
                elif not iBeg-iFin:     # got peak above high limit (peaks sorted, so we can stop)
                    break
                LaueFringePeakCalc(xdata,yc,lam,pos,intens,sig,gam,shol,ncells,clat,dampM,dampP,fmin,fitPowerM,fitPowerP,plot=False)
                if Ka2:
                    pos2 = pos+lamRatio*tand(pos/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    iBeg = np.searchsorted(xdata,pos2-fmin)
                    iFin = np.searchsorted(xdata,pos2+fmin)
                    if iBeg-iFin:
                        LaueFringePeakCalc(xdata,yc,lam2,pos2,intens*kRatio,sig,gam,shol,ncells,clat,dampM,dampP,fmin,fitPowerM,fitPowerP)
            except KeyError:        #no more peaks to process
                return yb+yc
    elif dataType[2] in ['A','B','C']:
        if 'B' not in dataType:
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
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                sig = max(sig,0.001)          #avoid neg sigma^2
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
                    gam = parmDict[gamName]
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                gam = max(gam,0.001)             #avoid neg gamma
                if 'C' not in dataType:
                    alpName = 'alp'+str(iPeak)
                    if alpName in varyList or not peakInstPrmMode:
                        alp = parmDict[alpName]
                    else:
                        alp = G2mth.getPinkAlpha(parmDict,tth)
                    alp = max(0.1,alp)
                    betName = 'bet'+str(iPeak)
                    if betName in varyList or not peakInstPrmMode:
                        bet = parmDict[betName]
                    else:
                        bet = G2mth.getPinkBeta(parmDict,tth)
                    bet = max(0.1,bet)
                if 'C' in dataType:
                    Wd,fmin,fmax = getWidthsCW(pos,sig,gam,shl)
                elif 'B' in dataType:
                    Wd,fmin,fmax = getWidthsCWB(pos,alp,bet,sig,gam)
                else: #'A'
                    Wd,fmin,fmax = getWidthsCWA(pos,alp,bet,sig,gam,shl)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    return yb+yc
                elif iFin-iBeg < 2:
                    iPeak += 1
                    continue
                if 'C' in dataType:
                    fp = getFCJVoigt3(pos,sig,gam,shl,xdata[iBeg:iFin])[0]
                elif 'B' in dataType:
                    fp = getEpsVoigt(pos,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])[0]/100.
                else: #'A'
                    fp = getExpFCJVoigt3(pos,alp,bet,sig,gam,shl,xdata[iBeg:iFin])[0]
                yc[iBeg:iFin] += intens*fp
                if Ka2:
                    pos2 = pos+lamRatio*tand(pos/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    iBeg = np.searchsorted(xdata,pos2-fmin)
                    iFin = np.searchsorted(xdata,pos2+fmin)
                    if iBeg-iFin:
                        if 'C' in dataType:
                            fp2 = getFCJVoigt3(pos2,sig,gam,shl,xdata[iBeg:iFin])[0]
                        elif 'B' in dataType:
                            fp2 = getEpsVoigt(pos2,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])[0]/100.
                        else: #'A'
                            fp2 = getExpFCJVoigt3(pos2,alp,bet,sig,gam,shl,xdata[iBeg:iFin])[0]
                        yc[iBeg:iFin] += intens*kRatio*fp2
                iPeak += 1
            except KeyError:        #no more peaks to process
                return yb+yc
    elif 'E' in dataType:
        iPeak = 0
        dsp = 1.0 #for now - fix later
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                intens = parmDict['int'+str(iPeak)]
                sigName = 'sig'+str(iPeak)
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getEDsig(parmDict,pos)
                sig = max(sig,0.001)          #avoid neg sigma^2
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
                    gam = parmDict[gamName]
                else:
                    gam = G2mth.getEDgam(parmDict,pos)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsED(pos,sig,gam)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = max(iBeg+3,np.searchsorted(xdata,pos+fmin))
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    return yb+yc
                yc[iBeg:iFin] += intens*getPsVoigt(pos,sig*10.**4,gam*100.,xdata[iBeg:iFin])[0]
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
                if alpName in varyList or not peakInstPrmMode:
                    alp = parmDict[alpName]
                else:
                    if len(Pdabc):
                        alp = np.interp(dsp,Pdabc[0],Pdabc[1])
                    else:
                        alp = G2mth.getTOFalpha(parmDict,dsp)
                alp = max(0.1,alp)
                betName = 'bet'+str(iPeak)
                if betName in varyList or not peakInstPrmMode:
                    bet = parmDict[betName]
                else:
                    if len(Pdabc):
                        bet = np.interp(dsp,Pdabc[0],Pdabc[2])
                    else:
                        bet = G2mth.getTOFbeta(parmDict,dsp)
                bet = max(0.0001,bet)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                else:
                    sig = G2mth.getTOFsig(parmDict,dsp)
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
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
                yc[iBeg:iFin] += intens*getEpsVoigt(pos,alp,bet,sig,gam,xdata[iBeg:iFin])[0]
                iPeak += 1
            except KeyError:        #no more peaks to process
                return yb+yc

def getPeakProfileDerv(dataType,parmDict,xdata,fixback,varyList,bakType):
    '''Computes the profile derivatives for a powder pattern for single peak fitting

    return: np.array([dMdx1,dMdx2,...]) in same order as varylist = backVary,insVary,peakVary order

    NB: not used for Rietveld refinement
    '''
    dMdv = np.zeros(shape=(len(varyList),len(xdata)))
    dMdb,dMddb,dMdpk,dMdfb = getBackgroundDerv('',parmDict,bakType,dataType,xdata,fixback)
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
    if 'LF' in dataType:
        for i,name in enumerate(varyList):
            if not np.all(dMdv[i] == 0): continue
            deltaParmDict = parmDict.copy()
            delta = max(parmDict[name]/1e5,0.001)
            deltaParmDict[name] += delta
            #print('num. deriv for',name,'val',deltaParmDict[name],'delta',delta)
            intArrP = getPeakProfile(dataType,deltaParmDict,xdata,fixback,varyList,bakType)
            deltaParmDict[name] -= 2*delta
            intArrM = getPeakProfile(dataType,deltaParmDict,xdata,fixback,varyList,bakType)
            dMdv[i] = 0.5 * (intArrP - intArrM) / delta
        return dMdv
    if dataType[2] in ['A','B','C']:
        if 'B' not in dataType:
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
                if 'C' not in dataType:
                    alpName = 'alp'+str(iPeak)
                    if alpName in varyList or not peakInstPrmMode:
                        alp = parmDict[alpName]
                        dada0 = dada1 = 0.0
                    else:
                        alp = G2mth.getPinkAlpha(parmDict,tth)
                        dada0,dada1 = G2mth.getPinkAlphaDeriv(tth)
                    alp = max(0.0001,alp)
                    betName = 'bet'+str(iPeak)
                    if betName in varyList or not peakInstPrmMode:
                        bet = parmDict[betName]
                        dbdb0 = dbdb1 = 0.0
                    else:
                        bet = G2mth.getPinkBeta(parmDict,tth)
                        dbdb0,dbdb1 = G2mth.getPinkBetaDeriv(tth)
                    bet = max(0.0001,bet)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                    dsdU = dsdV = dsdW = 0
                else:
                    sig = G2mth.getCWsig(parmDict,tth)
                    dsdU,dsdV,dsdW = G2mth.getCWsigDeriv(tth)
                sig = max(sig,0.001)          #avoid neg sigma
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
                    gam = parmDict[gamName]
                    dgdX = dgdY = dgdZ = 0
                else:
                    gam = G2mth.getCWgam(parmDict,tth)
                    dgdX,dgdY,dgdZ = G2mth.getCWgamDeriv(tth)
                gam = max(gam,0.001)             #avoid neg gamma
                if 'C' in dataType:
                    Wd,fmin,fmax = getWidthsCW(pos,sig,gam,shl)
                elif 'B' in dataType:
                    Wd,fmin,fmax = getWidthsCWB(pos,alp,bet,sig,gam)
                else: #'A'
                    Wd,fmin,fmax = getWidthsCWA(pos,alp,bet,sig,gam,shl)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = max(iBeg+3,np.searchsorted(xdata,pos+fmin))
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    break
                if 'C' in dataType:
                    dMdpk = np.zeros(shape=(6,len(xdata)))
                    dMdipk = getdFCJVoigt3(pos,sig,gam,shl,xdata[iBeg:iFin])
                    for i in range(1,5):
                        dMdpk[i][iBeg:iFin] += intens*dMdipk[i]
                    dMdpk[0][iBeg:iFin] += dMdipk[0]
                    dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4]}
                elif 'B' in dataType:
                    dMdpk = np.zeros(shape=(7,len(xdata)))
                    dMdipk = getdEpsVoigt(pos,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])
                    for i in range(1,6):
                        dMdpk[i][iBeg:iFin] += intens*dMdipk[i]/100.
                    dMdpk[0][iBeg:iFin] += dMdipk[0]/100.
                    dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4]/1.e4,'gam':dMdpk[5]/100.}
                else: #'A'
                    dMdpk = np.zeros(shape=(8,len(xdata)))
                    dMdipk = getdExpFCJVoigt3(pos,alp,bet,sig,gam,shl,xdata[iBeg:iFin])
                    for i in range(1,7):
                        dMdpk[i][iBeg:iFin] += intens*dMdipk[i]
                    dMdpk[0][iBeg:iFin] += dMdipk[0]
                    dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5],'shl':dMdpk[6]}
                if Ka2:
                    pos2 = pos+lamRatio*tand(pos/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    iBeg = np.searchsorted(xdata,pos2-fmin)
                    iFin = np.searchsorted(xdata,pos2+fmin)
                    if iBeg-iFin:
                        if 'C' in dataType:
                            dMdpk2 = np.zeros(shape=(6,len(xdata)))
                            dMdipk2 = getdFCJVoigt3(pos2,sig,gam,shl,xdata[iBeg:iFin])
                            for i in range(1,5):
                                dMdpk2[i][iBeg:iFin] += intens*kRatio*dMdipk2[i]
                            dMdpk[0][iBeg:iFin] += kRatio*dMdipk2[0]
                            dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':dMdpk[5]*intens}
                        elif 'B' in dataType:
                            dMdpk2 = np.zeros(shape=(7,len(xdata)))
                            dMdipk2 = getdEpsVoigt(pos,alp,bet,sig/1.e4,gam/100.,xdata[iBeg:iFin])
                            for i in range(1,6):
                                dMdpk2[i][iBeg:iFin] += intens*kRatio*dMdipk2[i]/100.
                            dMdpk[0][iBeg:iFin] += kRatio*dMdipk2[0]/100.
                            dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4]/1.e4,'gam':dMdpk[5]/100.,'L1/L2':dMdpk[6]*intens}
                        else: #'A'
                            dMdpk2 = np.zeros(shape=(8,len(xdata)))
                            dMdipk2 = getdExpFCJVoigt3(pos,alp,bet,sig,gam,shl,xdata[iBeg:iFin])
                            for i in range(1,7):
                                dMdpk2[i][iBeg:iFin] += intens*kRatio*dMdipk2[i]
                            dMdpk[0][iBeg:iFin] += kRatio*dMdipk2[0]
                            dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5],'shl':dMdpk[6],'L1/L2':dMdpk[7]*intens}
                for parmName in ['pos','int','alp','bet','sig','gam','shl','L1/L2']:
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
                    dMdv[varyList.index('SH/L')] += dervDict['shl']
                if 'alpha-0' in varyList:
                    dMdv[varyList.index('alpha-0')] += dada0*dervDict['alp']
                if 'alpha-1' in varyList:
                    dMdv[varyList.index('alpha-1')] += dada1*dervDict['alp']
                if 'beta-0' in varyList:
                    dMdv[varyList.index('beta-0')] += dbdb0*dervDict['bet']
                if 'beta-1' in varyList:
                    dMdv[varyList.index('beta-1')] += dbdb1*dervDict['bet']
                if 'I(L2)/I(L1)' in varyList:
                    dMdv[varyList.index('I(L2)/I(L1)')] += dervDict['L1/L2']
                iPeak += 1
            except KeyError:        #no more peaks to process
                break
    elif 'E' in dataType:
        iPeak = 0
        while True:
            try:
                pos = parmDict['pos'+str(iPeak)]
                intens = parmDict['int'+str(iPeak)]
                sigName = 'sig'+str(iPeak)
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                    dsdA = dsdB = dsdC = 0
                else:
                    sig = G2mth.getEDsig(parmDict,pos)
                    dsdA,dsdB,dsdC = G2mth.getEDsigDeriv(parmDict,pos)
                sig = max(sig,0.001)          #avoid neg sigma
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
                    gam = parmDict[gamName]
                    dgdX = dgdY = dgdZ = 0
                else:
                    gam = G2mth.getEDgam(parmDict,pos)
                    dgdX,dgdY,dgdZ = G2mth.getEDgamDeriv(parmDict,pos)
                gam = max(gam,0.001)             #avoid neg gamma
                Wd,fmin,fmax = getWidthsED(pos,sig,gam)
                iBeg = np.searchsorted(xdata,pos-fmin)
                iFin = np.searchsorted(xdata,pos+fmin)
                if not iBeg+iFin:       #peak below low limit
                    iPeak += 1
                    continue
                elif not iBeg-iFin:     #peak above high limit
                    break
                dMdpk = np.zeros(shape=(4,len(xdata)))
                dMdipk = getdPsVoigt(pos,sig*10.**4,gam*100.,xdata[iBeg:iFin])
                dMdpk[0][iBeg:iFin] += dMdipk[0]
                for i in range(1,4):
                    dMdpk[i][iBeg:iFin] += intens*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':-dMdpk[1],'sig':dMdpk[2]*10**4,'gam':dMdpk[3]*100.}
                for parmName in ['pos','int','sig','gam']:
                    try:
                        idx = varyList.index(parmName+str(iPeak))
                        dMdv[idx] = dervDict[parmName]
                    except ValueError:
                        pass
                if 'A' in varyList:
                    dMdv[varyList.index('A')] += dsdA*dervDict['sig']
                if 'B' in varyList:
                    dMdv[varyList.index('B')] += dsdB*dervDict['sig']
                if 'C' in varyList:
                    dMdv[varyList.index('C')] += dsdC*dervDict['sig']
                if 'X' in varyList:
                    dMdv[varyList.index('X')] += dgdX*dervDict['gam']
                if 'Y' in varyList:
                    dMdv[varyList.index('Y')] += dgdY*dervDict['gam']
                if 'Z' in varyList:
                    dMdv[varyList.index('Z')] += dgdZ*dervDict['gam']
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
                if alpName in varyList or not peakInstPrmMode:
                    alp = parmDict[alpName]
                else:
                    if len(Pdabc):
                        alp = np.interp(dsp,Pdabc[0],Pdabc[1])
                        dada0 = 0
                    else:
                        alp = G2mth.getTOFalpha(parmDict,dsp)
                        dada0 = G2mth.getTOFalphaDeriv(dsp)
                betName = 'bet'+str(iPeak)
                if betName in varyList or not peakInstPrmMode:
                    bet = parmDict[betName]
                else:
                    if len(Pdabc):
                        bet = np.interp(dsp,Pdabc[0],Pdabc[2])
                        dbdb0 = dbdb1 = dbdb2 = 0
                    else:
                        bet = G2mth.getTOFbeta(parmDict,dsp)
                        dbdb0,dbdb1,dbdb2 = G2mth.getTOFbetaDeriv(dsp)
                sigName = 'sig'+str(iPeak)
                if sigName in varyList or not peakInstPrmMode:
                    sig = parmDict[sigName]
                    dsds0 = dsds1 = dsds2 = dsds3 = 0
                else:
                    sig = G2mth.getTOFsig(parmDict,dsp)
                    dsds0,dsds1,dsds2,dsds3 = G2mth.getTOFsigDeriv(dsp)
                gamName = 'gam'+str(iPeak)
                if gamName in varyList or not peakInstPrmMode:
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
                    dMdpk[i][iBeg:iFin] += intens*dMdipk[i]
                dMdpk[0][iBeg:iFin] += dMdipk[0]
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
    if 'BF mult' in varyList:
        dMdv[varyList.index('BF mult')] = fixback

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
    if 'background PWDR' in Background[1]:
        backDict['Back File'] = Background[1]['background PWDR'][0]
        backDict['BF mult'] = Background[1]['background PWDR'][1]
        if len(Background[1]['background PWDR']) > 2:
            if Background[1]['background PWDR'][2]:
                backVary += ['BF mult',]
    return bakType,backDict,backVary

def autoBkgCalc(bkgdict,ydata):
    '''Compute the autobackground using the selected pybaselines function

    :param dict bkgdict: background parameters
    :param np.array ydata: array of Y values
    :returns: points for background intensity at each Y position
    '''
    try:
        import pybaselines.whittaker
    except:
        raise Exception("import of pybaselines failed. autobackground requires this")
    lamb = int(10**bkgdict['autoPrms']['logLam'])
    if bkgdict['autoPrms']['opt'] == 0:
        func = pybaselines.whittaker.arpls
    else:
        func = pybaselines.whittaker.iarpls
    return func(ydata, lam=lamb, max_iter=10)[0]

def DoCalibInst(IndexPeaks,Inst,Sample):

    def SetInstParms():
        dataType = Inst['Type'][0]
        insVary = []
        insNames = []
        insVals = []
        for parm in Inst:
            insNames.append(parm)
            insVals.append(Inst[parm][1])
            if parm in ['Lam','difC','difA','difB','Zero','2-theta','XE','YE','ZE']:
                if Inst[parm][2]:
                    insVary.append(parm)
        if 'C' in dataType or 'B' in dataType and 'Debye' in Sample['Type']:
            insNames.append('radius')
            insVals.append(Sample['Gonio. radius'])
            insNames.append('DisplaceX')
            insVals.append(Sample['DisplaceX'][0])
            if Sample['DisplaceX'][1]:
                insVary.append('DisplaceX')
            insNames.append('DisplaceY')
            insVals.append(Sample['DisplaceY'][0])
            if Sample['DisplaceY'][1]:
                insVary.append('DisplaceY')
        instDict = dict(zip(insNames,insVals))
        return dataType,instDict,insVary

    def GetInstParms(parmDict,varyList):
        for name in Inst:
            Inst[name][1] = parmDict[name]
        for name in Sample:
            if name in ['DisplaceX','DisplaceY']: # for CW only
                try:
                    Sample[name][0] = parmDict[name]
                except:
                    pass

    def InstPrint(sigDict):
        print ('Instrument/Sample Parameters:')
        if 'C' in Inst['Type'][0] or 'B' in Inst['Type'][0]:
            ptfmt = "%12.6f"
        else:
            ptfmt = "%12.3f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for parm in Inst:
            if parm in  ['Lam','difC','difA','difB','Zero','2-theta','XE','YE','ZE']:
                ptlbls += "%s" % (parm.center(12))
                ptstr += ptfmt % (Inst[parm][1])
                if parm in sigDict:
                    sigstr += ptfmt % (sigDict[parm])
                else:
                    sigstr += 12*' '
        for parm in Sample:
            if parm in  ['DisplaceX','DisplaceY']:
                ptlbls += "%s" % (parm.center(12))
                ptstr += ptfmt % (Sample[parm][0])
                if parm in sigDict:
                    sigstr += ptfmt % (sigDict[parm])
                else:
                    sigstr += 12*' '
        print (ptlbls)
        print (ptstr)
        print (sigstr)

    def errPeakPos(values,peakDsp,peakPos,peakWt,dataType,parmDict,varyList):
        parmDict.update(zip(varyList,values))
        calcPos = G2lat.getPeakPos(dataType,parmDict,peakDsp)
        if 'C' in dataType or 'B' in dataType:
            const = 18.e-2/(np.pi*parmDict['radius'])
            shft = -const*(parmDict['DisplaceX']*npcosd(calcPos)+parmDict['DisplaceY']*npsind(calcPos))+parmDict['Zero']
            # DThX = npasind(10**-3*parmDict['DisplaceX']*npcosd(calcPos)/parmDict['radius'])
            # DThY = -npasind(10**-3*parmDict['DisplaceY']*npsind(calcPos)/parmDict['radius'])
            # shft = DThX+DThY+parmDict['Zero']
            return np.sqrt(peakWt)*(calcPos+shft-peakPos)
        else:
            return np.sqrt(peakWt)*(calcPos-peakPos)

    peakPos = []
    peakDsp = []
    peakWt = []
    for peak,sig in zip(IndexPeaks[0],IndexPeaks[1]):
        if peak[2] and peak[3] and sig > 0.:
            peakPos.append(peak[0])
            peakDsp.append(peak[-1])    #d-calc
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
            return False

    sigDict = dict(zip(varyList,sig))
    GetInstParms(parmDict,varyList)
    InstPrint(sigDict)
    return True

def getHeaderInfo(dataType):
    '''Provide parameter name, label name and formatting information for the
    contents of the Peak Table and where used in DoPeakFit
    '''
    names = ['pos','int']
    lnames = ['position','intensity']
    if 'LF' in dataType:
        names = ['int','sig','gam','dampM','dampP','l','ttheta']
        lnames = ['intensity','sigma\u00b2','gamma','damping\nminus',
                      'damping\nplus','00l',
                      #'2theta    '
                      '2\u03B8'
                      ]
        fmt = ["%10.2f","%10.3f","%10.3f","%10.3f","%10.3f","%4.0f","%8.3f"]
    elif 'C' in dataType:
        names += ['sig','gam']
        lnames += ['sigma\u00b2','gamma']
        fmt = ["%10.5f","%10.1f","%10.3f","%10.3f"]
    elif 'T' in dataType:
        names += ['alp','bet','sig','gam']
        lnames += ['alpha','beta','sigma\u00b2','gamma']
        fmt = ["%10.2f","%10.4f","%8.3f","%8.5f","%10.3f","%10.3f"]
    elif 'E' in dataType:
        names += ['sig','gam']
        lnames += ['sigma\u00b2','gamma']
        fmt = ["%10.5f","%10.1f","%8.3f","%10.3f"]
    else: # 'B'
        names += ['alp','bet','sig','gam']
        lnames += ['alpha','beta','sigma\u00b2','gamma']
        fmt = ["%10.5f","%10.1f","%8.2f","%8.4f","%10.3f","%10.3f"]
    return names, fmt, lnames

def DoPeakFit(FitPgm,Peaks,Background,Limits,Inst,Inst2,data,fixback=None,prevVaryList=[],
                  oneCycle=False,controls=None,wtFactor=1.0,dlg=None,noFit=False):
    '''Called to perform a peak fit, refining the selected items in the peak
    table as well as selected items in the background.

    :param str FitPgm: type of fit to perform. At present this is ignored.
    :param list Peaks: a list of peaks. Each peak entry is a list with paired values:
      The number of pairs depends on the data type (see :func:`getHeaderInfo`).
      For CW data there are
      four values each followed by a refine flag where the values are: position, intensity,
      sigma (Gaussian width) and gamma (Lorentzian width). From the Histogram/"Peak List"
      tree entry, dict item "peaks". For some types of fits, overall parameters are placed
      in a dict entry.
    :param list Background: describes the background. List with two items.
      Item 0 specifies a background model and coefficients. Item 1 is a dict.
      From the Histogram/Background tree entry.
    :param list Limits: min and max x-value to use
    :param dict Inst: Instrument parameters
    :param dict Inst2: more Instrument parameters
    :param numpy.array data: a 5xn array. data[0] is the x-values,
      data[1] is the y-values, data[2] are weight values, data[3], [4] and [5] are
      calc, background and difference intensities, respectively.
    :param array fixback: fixed background array; same size as data[0-5]
    :param list prevVaryList: Used in sequential refinements to override the
      variable list. Defaults as an empty list.
    :param bool oneCycle: True if only one cycle of fitting should be performed
    :param dict controls: a dict specifying two values, Ftol = controls['min dM/M']
      and derivType = controls['deriv type']. If None default values are used.
    :param float wtFactor: weight multiplier; = 1.0 by default
    :param wx.Dialog dlg: A dialog box that is updated with progress from the fit.
      Defaults to None, which means no updates are done.
    :param bool noFit: When noFit is True, a refinement is not performed. Default
      is False.

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
        if 'BF mult' in parmList:
            Background[1]['background PWDR'][1] = parmList['BF mult']

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
        if 'BF mult' in sigDict:
            print('Background file mult: %.3f(%d)'%(Background[1]['background PWDR'][1],int(1000*sigDict['BF mult'])))

    def SetInstParms(Inst):
        dataType = Inst['Type'][0]
        insVary = []
        insNames = []
        insVals = []
        for parm in Inst:
            insNames.append(parm)
            insVals.append(Inst[parm][1])
            if parm in ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha','A','B','C',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','alpha-0','alpha-1'] and Inst[parm][2]:
                    insVary.append(parm)
        instDict = dict(zip(insNames,insVals))
        if 'SH/L' in instDict:
            instDict['SH/L'] = max(instDict['SH/L'],0.002)
        return dataType,instDict,insVary

    def GetPkInstParms(parmDict,Inst,varyList):
        for name in Inst:
            Inst[name][1] = parmDict[name]
        iPeak = 0
        while True:
            try:
                sigName = 'sig'+str(iPeak)
                pos = parmDict['pos'+str(iPeak)]
                if sigName not in varyList and peakInstPrmMode:
                    if 'T' in Inst['Type'][0]:
                        dsp = G2lat.Pos2dsp(Inst,pos)
                        parmDict[sigName] = G2mth.getTOFsig(parmDict,dsp)
                    if 'E' in Inst['Type'][0]:
                        parmDict[sigName] = G2mth.getEDsig(parmDict,pos)
                    else:
                        parmDict[sigName] = G2mth.getCWsig(parmDict,pos)
                gamName = 'gam'+str(iPeak)
                if gamName not in varyList and peakInstPrmMode:
                    if 'T' in Inst['Type'][0]:
                        dsp = G2lat.Pos2dsp(Inst,pos)
                        parmDict[gamName] = G2mth.getTOFgamma(parmDict,dsp)
                    if 'E' in Inst['Type'][0]:
                        parmDict[gamName] = G2mth.getEDgam(parmDict,pos)
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
            if parm in  ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha','A','B','C',
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
        '''Set the contents of peakDict from list Peaks
        '''
        peakDict = {}
        peakVary = []
        names,_,_ = getHeaderInfo(dataType)
        if 'LF' in dataType:
            off = 2
            names = names[:-2] # drop 00l & 2theta from header
        else:
            off = 0
        for i,peak in enumerate(Peaks):
            if type(peak) is dict:
                peakDict.update(peak)
                continue
            if 'LF' in dataType: peakDict['l'+str(i)] = peak[12]
            for j,name in enumerate(names):
                parName = name+str(i)
                peakDict[parName] = peak[off+2*j]
                if peak[off+2*j+1]:
                    peakVary.append(parName)
        return peakDict,peakVary

    def GetPeaksParms(Inst,parmDict,Peaks,varyList):
        '''Put values into the Peaks list from the refinement results from inside
        the parmDict array
        '''
        names,_,_ = getHeaderInfo(Inst['Type'][0])
        off = 0
        if 'LF' in Inst['Type'][0]:
            off = 2
            if 'clat' in varyList:
                Peaks[-1]['clat'] = parmDict['clat']
            names = names[:-1] # drop 2nd 2theta value
            for i,peak in enumerate(Peaks):
                if type(peak) is dict: continue
                parmDict['ttheta'+str(i)] = peak[-1]
        for i,peak in enumerate(Peaks):
            if type(peak) is dict:
                continue
            for j in range(len(names)):
                parName = names[j]+str(i)
                if parName in varyList or not peakInstPrmMode:
                    peak[2*j+off] = parmDict[parName]
            if 'pos'+str(i) not in parmDict: continue
            pos = parmDict['pos'+str(i)]
            if 'LF' in Inst['Type'][0]:
                peak[0] = pos
                peak[-1] = pos
            if 'difC' in Inst:
                dsp = pos/Inst['difC'][1]
            for j in range(len(names)):
                parName = names[j]+str(i)
                if peak[2*j+off + 1] or not peakInstPrmMode: continue
                if 'alp' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getTOFalpha(parmDict,dsp)
                    else: #'B'
                        peak[2*j+off] = G2mth.getPinkAlpha(parmDict,pos)
                elif 'bet' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getTOFbeta(parmDict,dsp)
                    else:   #'B'
                        peak[2*j+off] = G2mth.getPinkBeta(parmDict,pos)
                elif 'sig' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getTOFsig(parmDict,dsp)
                    elif 'E' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getEDsig(parmDict,pos)
                    else:   #'C' & 'B'
                        peak[2*j+off] = G2mth.getCWsig(parmDict,pos)
                elif 'gam' in parName:
                    if 'T' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getTOFgamma(parmDict,dsp)
                    elif 'E' in Inst['Type'][0]:
                        peak[2*j+off] = G2mth.getEDgam(parmDict,pos)
                    else:   #'C' & 'B'
                        peak[2*j+off] = G2mth.getCWgam(parmDict,pos)

    def PeaksPrint(dataType,parmDict,sigDict,varyList,ptsperFW):
        if 'clat' in varyList:
            print('c = {:.6f} esd {:.6f}'.format(parmDict['clat'],sigDict['clat']))
        print ('Peak coefficients:')
        names,fmt,_ = getHeaderInfo(dataType)
        head = 13*' '
        for name in names:
            if name == 'l':
                head += name
            elif name == 'ttheta':
                head += name.center(8)
            elif name in ['alp','bet']:
                head += name.center(8)+'esd'.center(8)
            else:
                head += name.center(10)+'esd'.center(10)
        head += 'bins'.center(12)
        print (head)
        ptfmt = dict(zip(names,fmt))
        for i,peak in enumerate(Peaks):
            if type(peak) is dict:
                continue
            ptstr =  ':'
            for j in range(len(names)):
                name = names[j]
                parName = name+str(i)
                if parName not in parmDict: continue
                ptstr += ptfmt[name] % (parmDict[parName])
                if name == 'l' or name == 'ttheta':
                    continue
                if parName in varyList:
                    ptstr += ptfmt[name] % (sigDict[parName])
                else:
                    if name in ['alp','bet']:
                        ptstr += 8*' '
                    else:
                        ptstr += 10*' '
            ptstr += '%8.1f'%(ptsperFW[i])
            print ('%s'%(('Peak'+str(i+1)).center(8)),ptstr)

    def devPeakProfile(values,xdata,ydata,fixback, weights,dataType,parmdict,varylist,bakType,dlg):
        '''Computes a matrix where each row is the derivative of the calc-obs
        values (see :func:`errPeakProfile`) with respect to each parameter
        in backVary,insVary,peakVary. Used for peak fitting.
        '''
        parmdict.update(zip(varylist,values))
        return np.sqrt(weights)*getPeakProfileDerv(dataType,parmdict,xdata,fixback,varylist,bakType)

    def errPeakProfile(values,xdata,ydata,fixback,weights,dataType,parmdict,varylist,bakType,dlg):
        '''Computes a vector with the weighted calc-obs values differences
        for peak fitting
        '''
        parmdict.update(zip(varylist,values))
        M = np.sqrt(weights)*(getPeakProfile(dataType,parmdict,xdata,fixback,varylist,bakType)-ydata)
        Rwp = min(100.,np.sqrt(np.sum(M**2)/np.sum(weights*ydata**2))*100.)
        if dlg:
            dlg.Raise()
            GoOn = dlg.Update(int(Rwp),newmsg='%s%8.3f%s'%('Peak fit Rwp =',Rwp,'%'))[0]
            if not GoOn:
                return -M           #abort!!
        return M

    #---- beginning of DoPeakFit ---------------------------------------------
    if controls:
        Ftol = controls['min dM/M']
    else:
        Ftol = 0.0001
    if oneCycle:
        Ftol = 1.0
    x,y,w,yc,yb,yd = data   #these are numpy arrays - remove masks!
    if fixback is None:
        fixback = np.zeros_like(y)
    yc.fill(0.)                            #set calcd ones to zero
    yb.fill(0.)
    yd.fill(0.)
    xBeg = np.searchsorted(x,Limits[0])
    xFin = np.searchsorted(x,Limits[1])+1
    # find out what is varied
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
        if 'LF' in Inst['Type'][0] and Peaks:
            if Peaks[-1].get('clat-ref'): varyList += ['clat']
    fullvaryList = varyList[:]
    remVary = []
    if not peakInstPrmMode:
        for v in ('U','V','W','X','Y','Z','alpha','alpha-0','alpha-1','A','B','C',
            'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',):
            if v in varyList:
                remVary.append(v)
        if remVary:
            print('Instrumental profile terms cannot be varied '+
                          'after setPeakInstPrmMode(False) is used'+
                          f'\nremoving vars: {" ,".join(remVary)}')
            varyList = [v for v in varyList if v not in remVary]
    if 'LF' in Inst['Type'][0]:
        warn = []
        for v in ('U','V','W','X','Y','Z','alpha','alpha-0','alpha-1',
            'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',):
            if v in varyList:
                warn.append(v)
                del varyList[varyList.index(v)]
        if warn:
            print('Instrumental profile terms cannot be varied '+
                                    'in Laue Fringe fits:',warn)

    while not noFit:
        begin = time.time()
        values =  np.array(Dict2Values(parmDict, varyList))
        Rvals = {}
        badVary = []
        try:
            result = so.leastsq(errPeakProfile,values,Dfun=devPeakProfile,full_output=True,ftol=Ftol,col_deriv=True,
               args=(x[xBeg:xFin],y[xBeg:xFin],fixback[xBeg:xFin],wtFactor*w[xBeg:xFin],dataType,parmDict,varyList,bakType,dlg))
        except Exception as msg:
            if GSASIIpath.GetConfigValue('debug'):
                print('peak fit failure\n',msg)
                import traceback
                print (traceback.format_exc())
            else:
                print('peak fit failure')
            return
        ncyc = int(result[2]['nfev']/2)
        runtime = time.time()-begin
        chisq = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
        Rvals['Rwp'] = np.sqrt(chisq/np.sum(wtFactor*w[xBeg:xFin]*y[xBeg:xFin]**2))*100.      #to %
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
    yb[xBeg:xFin] = getBackground('',parmDict,bakType,dataType,x[xBeg:xFin],fixback[xBeg:xFin])[0]
    yc[xBeg:xFin] = getPeakProfile(dataType,parmDict,x[xBeg:xFin],fixback[xBeg:xFin],varyList,bakType)
    yd[xBeg:xFin] = y[xBeg:xFin]-yc[xBeg:xFin]
    if noFit:
        GetPeaksParms(Inst,parmDict,Peaks,varyList)
        return
    sigDict = dict(zip(varyList,sig))
    GetBackgroundParms(parmDict,Background)
    if bakVary: BackgroundPrint(Background,sigDict)
    GetPkInstParms(parmDict,Inst,varyList)
    if insVary: InstPrint(Inst,sigDict)
    GetPeaksParms(Inst,parmDict,Peaks,varyList)
    binsperFWHM = []
    for peak in Peaks:
        if type(peak) is dict:
            continue
        FWHM = getFWHM(peak[0],Inst)
        try:
            xpk = x.searchsorted(peak[0])
            cw = x[xpk]-x[xpk-1]
            binsperFWHM.append(FWHM/cw)
        except IndexError:
            binsperFWHM.append(0.)
    if peakVary: PeaksPrint(dataType,parmDict,sigDict,varyList,binsperFWHM)
    if len(binsperFWHM):
        if min(binsperFWHM) < 1.:
            G2fil.G2Print ('*** Warning: calculated peak widths are too narrow to refine profile coefficients ***')
            G2fil.G2Print (' Manually make individually refined sigma/gamma terms positive')
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

#### RMCutilities ################################################################################
def MakeInst(PWDdata,Name,Size,Mustrain,useSamBrd):
    inst = PWDdata['Instrument Parameters'][0]
    sample = PWDdata['Sample Parameters']
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
                'difC','difA','Zero','2-theta','difB',
                'alpha','beta-0','beta-1','beta-q',
                'sig-0','sig-1','sig-2','sig-q',
                'Z','X','Y']
        fname = Name+'.inst'
        fl = open(fname,'w')
        fl.write('1\n')
        fl.write('%d\n'%int(inst[prms[0]][1]))
        fl.write('%19.11f%19.11f%19.11f%19.11f%19.11f\n'%(inst[prms[1]][1],inst[prms[2]][1],inst[prms[3]][1],inst[prms[4]][1],inst[prms[5]][1],))
        fl.write('%12.6e%14.6e%14.6e%14.6e\n'%(inst[prms[6]][1],inst[prms[7]][1],inst[prms[8]][1],inst[prms[9]][1]))
        fl.write('%12.6e%14.6e%14.6e%14.6e\n'%(inst[prms[10]][1],inst[prms[11]][1],inst[prms[12]][1],inst[prms[13]][1]))
        fl.write('%12.6e%14.6e%14.6e%14.6e%14.6e\n'%(inst[prms[14]][1],inst[prms[15]][1]+Ysb,inst[prms[16]][1]+Xsb,0.0,0.0))
        fl.write('%12.6e\n\n\n'%(sample['Absorption'][0]))
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
        fl.write('%10.5f%10.5f%10.4f%10d\n'%(inst[prms[1]][1],100.*inst[prms[2]][1],inst[prms[3]][1],0))
        fl.write('%10.3f%10.3f%10.3f\n'%(inst[prms[4]][1],inst[prms[5]][1],inst[prms[6]][1]))
        fl.write('%10.3f%10.3f%10.3f\n'%(inst[prms[7]][1]+Xsb,inst[prms[8]][1]+Ysb,0.0))
        fl.write('%10.3f%10.3f%10.3f\n'%(0.0,0.0,0.0))
        fl.write('%12.6e\n\n\n'%(sample['Absorption'][0]))
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
    Natoms = RMCPdict['NoAtoms']
    sumatms = np.sum(np.array([Natoms[iatm] for iatm in Natoms]))
    Isotope = RMCPdict['Isotope']
    Isotopes = RMCPdict['Isotopes']
    Atypes = RMCPdict['aTypes']
    if 'Va' in Atypes:
        Isotope['Va'] = 'Nat. Abund.'
        Isotopes['Va'] = {'Nat. Abund.':{'SL':[0.0,0.0]}}
    atPairs = RMCPdict['Pairs']
    Files = RMCPdict['files']
    BraggWt = RMCPdict['histogram'][1]
    inst = PWDdata['Instrument Parameters'][0]
    try:
        #pName = Phase['General']['Name']
        refList = PWDdata['Reflection Lists'][Name]['RefList']
    except TypeError:
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
    Ncoeff = []
    Nblen = [Isotopes[at][Isotope[at]]['SL'][0] for at in Atypes]
    for pair in [[' %s-%s'%(Atseq[i],Atseq[j]) for j in range(i,lenA)] for i in range(lenA)]:
        Pairs += pair
    for pair in Pairs:
        pair = pair.replace(' ','')
        at1,at2 = pair.split('-')
        if at1 == 'Va' or at2 == 'Va':
            ncoef = 0.0
        else:
            ncoef = Isotopes[at1][Isotope[at1]]['SL'][0]*Natoms[at1]/sumatms
            ncoef *= Isotopes[at2][Isotope[at2]]['SL'][0]*Natoms[at2]/sumatms
        if at1 != at2:
            ncoef *= 2.
        Ncoeff += [ncoef,]
    pairMin = [atPairs[pair] if pair in atPairs else [0.0,0.,0.] for pair in Pairs ]
    maxMoves = [Atypes[atm] if atm in Atypes else 0.0 for atm in Atseq ]
    fname = Name+'.dat'
    fl = open(fname,'w')
    fl.write(' %% Hand edit the following as needed\n')
    fl.write('TITLE :: '+Name+'\n')
    fl.write('MATERIAL :: '+Meta['material']+'\n')
    fl.write('PHASE :: '+Meta['phase']+'\n')
    fl.write('TEMPERATURE :: '+str(Meta['temperature'])+'\n')
    fl.write('INVESTIGATOR :: '+Meta['owner']+'\n')
    if RMCPdict.get('useGPU',False):
        fl.write('GPU_ACCELERATOR :: 0\n')
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
    if 'T' in inst['Type'][1]:
        fl.write('NEUTRON_COEFFICIENTS :: '+''.join(['%9.5f'%coeff for coeff in Ncoeff])+'\n')
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
    if len(RMCPdict['Swaps']):
        fl.write('\n')
        fl.write('SWAP_MULTI ::\n')
        for swap in RMCPdict['Swaps']:
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
                fqdata.close()
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
    if 'T' in inst['Type'][1]:
        fl.write('  > SCATTERING LENGTH :: '+''.join(['%8.4f'%blen for blen in Nblen])+'\n')
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

def findfullrmc():
    '''Find where fullrmc is installed. Tries the following:

         1. Returns the Config var 'fullrmc_exec', if defined. If an executable
            is found at that location it is assumed to run and supply
            fullrmc 5.0+
         2. The path is checked for a fullrmc image as named by Bachir

    :returns: the full path to a python executable that is assumed to
      have fullrmc installed or None, if it was not found.
    '''
    fullrmc_exe = GSASIIpath.GetConfigValue('fullrmc_exec')
    if fullrmc_exe is not None and is_exe(fullrmc_exe):
        return fullrmc_exe
    pathlist = os.environ["PATH"].split(os.pathsep)
    for p in (GSASIIpath.path2GSAS2,GSASIIpath.binaryPath,os.getcwd(),
                  os.path.split(sys.executable)[0]):
        if p not in pathlist: pathlist.append(p)
    import glob
    for p in pathlist:
        if sys.platform == "win32":
            lookfor = "fullrmc5*.exe"
        else:
            lookfor = "fullrmc5*64bit"
        fl = glob.glob(os.path.join(p,lookfor))
        if len(fl) > 0:
            fullrmc_exe = os.path.abspath(sorted(fl)[0])
            if GSASIIpath.GetConfigValue('debug'):
                print('fullrmc found as',fullrmc_exe)
            return fullrmc_exe

def fullrmcDownload():
    '''Downloads the fullrmc executable from Bachir's site to the current
    GSAS-II binary directory.

    Does some error checking.
    '''
    import os
    import requests
    import platform
    if platform.architecture()[0] != '64bit':
        return "fullrmc is only available for 64 bit machines. This is 32 bit"
    setXbit = True
    if sys.platform == "darwin":
        URL = "https://github.com/bachiraoun/fullrmc/raw/master/standalones/fullrmc500_3p8p6_macOS-10p16-x86_64-i386-64bit"
    elif sys.platform == "win32":
        setXbit = False
        URL = "https://github.com/bachiraoun/fullrmc/raw/master/standalones/fullrmc500_3p8p10_Windows-10-10p0p19041-SP0.exe"
    else:
        if 'aarch' in platform.machine() or 'arm' in platform.machine():
            return "Sorry, fullrmc is only available for Intel-compatible machines."
        URL = "https://github.com/bachiraoun/fullrmc/raw/master/standalones/fullrmc500_3p8p5_Linux-4p19p121-linuxkit-x86_64-with-glibc2p29"

    GSASIIpath.SetBinaryPath()
    fil = os.path.join(GSASIIpath.binaryPath,os.path.split(URL)[1])
    print('Starting installation of fullrmc\nDownloading from',
              'https://github.com/bachiraoun/fullrmc/tree/master/standalones',
              '\nCreating '+fil,
              '\nThis may take a while...')
    with open(fil, "wb") as fp:
        fp.write(requests.get(URL).content)
    print('...Download completed')
    if setXbit:
        import stat
        os.chmod(fil, os.stat(fil).st_mode | stat.S_IEXEC)
    return ''

def findPDFfit():
    '''Find if PDFfit2 is installed (may be local to GSAS-II). Does the following:
    :returns: the full path to a python executable or None, if it was not found.
    '''
    if GSASIIpath.GetConfigValue('pdffit2_exec') is not None and is_exe(
            GSASIIpath.GetConfigValue('pdffit2_exec')):
        return GSASIIpath.GetConfigValue('pdffit2_exec')
    try:
        from diffpy.pdffit2 import PdfFit
        import diffpy
        PdfFit
        diffpy
        return sys.executable
    except Exception as msg:
        print('Error importing PDFfit2:\n',msg)
        return None

def GetPDFfitAtomVar(Phase,RMCPdict):
    ''' Find dict of independent "@n" variables for PDFfit in atom constraints
    '''
    General = Phase['General']
    Atoms = Phase['Atoms']
    cx,ct,cs,cia = General['AtomPtrs']
    AtomVar = RMCPdict['AtomVar']
    varnames = []
    for iat,atom in enumerate(RMCPdict['AtomConstr']):
        for it,item in enumerate(atom):
            if it > 1 and item:
                itms = item.split('@')
                for itm in itms[1:]:
                    itnum = itm[:2]
                    varname = '@%s'%itnum
                    varnames.append(varname)
                    if it < 6:
                        if varname not in AtomVar:
                            AtomVar[varname] = 0.0      #put ISODISTORT mode displ here?
                    else:
                        for i in range(3):
                            if varname not in AtomVar:
                                AtomVar[varname] = Atoms[iat][cia+i+2]
    varnames = set(varnames)
    for name in list(AtomVar.keys()):       #clear out unused parameters
        if name not in varnames:
            del AtomVar[name]

def MakePDFfitAtomsFile(Phase,RMCPdict):
    '''Make the PDFfit atoms file
    '''
    General = Phase['General']
    if General['SGData']['SpGrp'] != 'P 1':
        return 'Space group symmetry must be lowered to P 1 for PDFfit'
    fName = General['Name']+'-PDFfit.stru'
    fName = fName.replace(' ','_')
    if 'sequential' in RMCPdict['refinement']:
        fName = 'Sequential_PDFfit.stru'
    fatm = open(fName,'w')
    fatm.write('title  structure of '+General['Name']+'\n')
    fatm.write('format pdffit\n')
    fatm.write('scale   1.000000\n')    #fixed
    sharp = '%10.6f,%10.6f,%10.6f,%10.6f\n'%(RMCPdict['delta2'][0],RMCPdict['delta1'][0],RMCPdict['sratio'][0],RMCPdict['rcut'])
    fatm.write('sharp '+sharp)
    shape = ''
    if RMCPdict['shape'] == 'sphere' and RMCPdict['spdiameter'][0] > 0.:
        shape = '   sphere, %10.6f\n'%RMCPdict['spdiameter'][0]
    elif RMCPdict['stepcut'] > 0.:
        shape = 'stepcut, %10.6f\n'%RMCPdict['stepcut']
    if shape:
        fatm.write('shape  '+shape)
    fatm.write('spcgr   %s\n'%RMCPdict['SGData']['SpGrp'].replace(' ',''))
    cell = General['Cell'][1:7]
    fatm.write('cell  %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f\n'%(
        cell[0],cell[1],cell[2],cell[3],cell[4],cell[5]))
    fatm.write('dcell '+5*'  0.000000,'+'  0.000000\n')
    Atoms = Phase['Atoms']
    fatm.write('ncell %8d,%8d,%8d,%10d\n'%(1,1,1,len(Atoms)))
    fatm.write('atoms\n')
    cx,ct,cs,cia = General['AtomPtrs']
    for atom in Atoms:
        fatm.write('%4s%18.8f%18.8f%18.8f%13.4f\n'%(atom[ct][:2].ljust(2),atom[cx],atom[cx+1],atom[cx+2],atom[cx+3]))
        fatm.write('    '+'%18.8f%18.8f%18.8f%13.4f\n'%(0.,0.,0.,0.))
        fatm.write('    '+'%18.8f%18.8f%18.8f\n'%(atom[cia+2],atom[cia+3],atom[cia+4]))
        fatm.write('    '+'%18.8f%18.8f%18.8f\n'%(0.,0.,0.,))
        fatm.write('    '+'%18.8f%18.8f%18.8f\n'%(atom[cia+5],atom[cia+6],atom[cia+7]))
        fatm.write('    '+'%18.8f%18.8f%18.8f\n'%(0.,0.,0.))
    fatm.close()

def MakePDFfitRunFile(Phase,RMCPdict):
    '''Make the PDFfit python run file
    '''

    def GetCellConstr(SGData):
        if SGData['SGLaue'] in ['m3', 'm3m']:
            return [1,1,1,0,0,0]
        elif SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm','4/m','4/mmm']:
            return [1,1,2,0,0,0]
        elif SGData['SGLaue'] in ['3R','3mR']:
            return [1,1,1,2,2,2]
        elif SGData['SGLaue'] == 'mmm':
            return [1,2,3,0,0,0]
        elif SGData['SGLaue'] == '2/m':
            if SGData['SGUniq'] == 'a':
                return [1,2,3,4,0,0]
            elif SGData['SGUniq'] == 'b':
                return [1,2,3,0,4,0]
            elif SGData['SGUniq'] == 'c':
                return [1,2,3,0,0,4]
        else:
            return [1,2,3,4,5,6]

    General = Phase['General']
    Cell = General['Cell'][1:7]
    rundata = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys,os
datadir = r'{:}'
pathWrap = lambda f: os.path.join(datadir,f)
'''.format(os.path.abspath(os.getcwd()))
    PDFfit_exe = findPDFfit()  # returns python loc
    if not PDFfit_exe:
        print('PDFfit2 is not found. Creating .sh file without paths.')
    rundata += 'from diffpy.pdffit2 import PdfFit\n'
    rundata += 'pf = PdfFit()\n'
    Nd = 0
    Np = 0
    parms = {}
    parmNames = {}
    if 'sequential' in RMCPdict['refinement']:
        Np = 3
        rundata += '#sequential data here\n'
    else:
        for fil in RMCPdict['files']:
            filNam = RMCPdict['files'][fil][0]
            if 'Select' in filNam:
                continue
            if 'Neutron' in fil:
                Nd += 1
                dType = 'Ndata'
            else:
                Nd += 1
                dType = 'Xdata'
            rundata += "pf.read_data(pathWrap(r'%s'), '%s', 30.0, %.4f)\n"%(filNam,dType[0],RMCPdict[dType]['qdamp'][0])
            rundata += 'pf.setdata(%d)\n'%Nd
            rundata += 'pf.pdfrange(%d, %6.2f, %6.2f)\n'%(Nd,RMCPdict[dType]['Fitrange'][0],RMCPdict[dType]['Fitrange'][1])
            for item in ['dscale','qdamp','qbroad']:
                if RMCPdict[dType][item][1]:
                    Np += 1
                    rundata += 'pf.constrain(pf.%s(),"@%d")\n'%(item,Np)
                    parms[Np] = RMCPdict[dType][item][0]
                    parmNames[Np] = item
    fName = General['Name']+'-PDFfit.stru'
    fName = fName.replace(' ','_')
    if 'sequential' in RMCPdict['refinement']:
        fName = 'Sequential_PDFfit.stru'
    Np = 9
    rundata += "pf.read_struct(pathWrap(r'{:}'))\n".format(fName)
    for item in ['delta1','delta2','sratio']:
        if RMCPdict[item][1]:
            Np += 1
            rundata += 'pf.constrain(pf.%s,"@%d")\n'%(item,Np)
            parms[Np] = RMCPdict[item][0]
            parmNames[Np] = item
    if 'sphere' in RMCPdict['shape'] and RMCPdict['spdiameter'][1]:
        Np += 1
        rundata += 'pf.constrain(pf.spdiameter,"@%d")\n'%Np
        parms[Np] = RMCPdict['spdiameter'][0]
        parmNames[Np] = 'spdiameter'

    if RMCPdict['cellref']:
        cellconst = GetCellConstr(RMCPdict['SGData'])
        used = []
        cellNames = ['a','b','c','alpha','beta','gamma']
        for ic in range(6):
            if cellconst[ic]:
                rundata += 'pf.constrain(pf.lat(%d), "@%d")\n'%(ic+1,Np+cellconst[ic])
                if cellconst[ic] not in used:
                    parms[Np+cellconst[ic]] = Cell[ic]
                    parmNames[Np+cellconst[ic]] = cellNames[ic]
                used.append(cellconst[ic])
#Atom constraints here -------------------------------------------------------
    AtomVar = RMCPdict['AtomVar']
    used = []
    for iat,atom in enumerate(RMCPdict['AtomConstr']):
        for it,item in enumerate(atom):
            names = ['pf.x(%d)'%(iat+1),'pf.y(%d)'%(iat+1),'pf.z(%d)'%(iat+1),'pf.occ(%d)'%(iat+1)]
            if it > 1 and item:
                itms = item.split('@')
                once = False
                for itm in itms[1:]:
                    try:
                        itnum = int(itm[:2])
                    except ValueError:
                        print(' *** ERROR - invalid string in atom constraint %s ***'%(item))
                        return None
                    if it < 6:
                        if not once:
                            rundata += 'pf.constrain(%s,"%s")\n'%(names[it-2],item)
                            once = True
                        if itnum not in used:
                            parms[itnum] = AtomVar['@%d'%itnum]
                            parmNames[itnum] = names[it-2].split('.')[1]
                            used.append(itnum)
                    else:
                        uijs = ['pf.u11(%d)'%(iat+1),'pf.u22(%d)'%(iat+1),'pf.u33(%d)'%(iat+1)]
                        for i in range(3):
                            rundata += 'pf.constrain(%s,"%s")\n'%(uijs[i],item)
                            if itnum not in used:
                                parms[itnum] = AtomVar['@%d'%itnum]
                                parmNames[itnum] = uijs[i].split('.')[1]
                                used.append(itnum)

    if 'sequential' in RMCPdict['refinement']:
        rundata += '#parameters here\n'
        RMCPdict['Parms'] = parms           #{'n':val,...}
        RMCPdict['ParmNames'] = parmNames   #{'n':name,...}
    else:
# set parameter values
        for iprm in parms:
            rundata += 'pf.setpar(%d,%.6f)\n'%(iprm,parms[iprm])

# Save results ---------------------------------------------------------------
    rundata += 'pf.refine()\n'
    if 'sequential' in RMCPdict['refinement']:
        fName = 'Sequential_PDFfit'
        rfile = open('Seq_PDFfit_template.py','w')
        rundata += 'pf.save_pdf(1, pathWrap("%s"))\n'%(fName+'.fgr')
    else:
        fName = General['Name'].replace(' ','_')+'-PDFfit'
        rfile = open(fName+'.py','w')
        Nd = 0
        for file in RMCPdict['files']:
            if 'Select' in RMCPdict['files'][file][0]:      #skip unselected
                continue
            Nd += 1
            rundata += 'pf.save_pdf(%d, pathWrap("%s"))\n'%(Nd,fName+file[0]+'.fgr')

    rundata += 'pf.save_struct(1, pathWrap("%s"))\n'%(fName+'.rstr')
    rundata += 'pf.save_res(pathWrap("%s"))\n'%(fName+'.res')

    rfile.writelines(rundata)
    rfile.close()

    return fName+'.py'

def GetSeqCell(SGData,parmDict):
    ''' For use in processing PDFfit sequential results
    '''
    try:
        if SGData['SGLaue'] in ['m3', 'm3m']:
            cell = [parmDict['11'][0],parmDict['11'][0],parmDict['11'][0],90.,90.,90.]
        elif SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm','4/m','4/mmm']:
            cell = [parmDict['11'][0],parmDict['11'][0],parmDict['12'][0],90.,90.,90.]
        elif SGData['SGLaue'] in ['3R','3mR']:
            cell = [parmDict['11'][0],parmDict['11'][0],parmDict['11'][0],
                parmDict['12'][0],parmDict['12'][0],parmDict['12'][0]]
        elif SGData['SGLaue'] == 'mmm':
            cell = [parmDict['11'][0],parmDict['12'][0],parmDict['13'][0],90.,90.,90.]
        elif SGData['SGLaue'] == '2/m':
            if SGData['SGUniq'] == 'a':
                cell = [parmDict['11'][0],parmDict['12'][0],parmDict['13'][0],parmDict['14'][0],90.,90.]
            elif SGData['SGUniq'] == 'b':
                cell = [parmDict['11'][0],parmDict['12'][0],parmDict['13'][0],90.,parmDict['14'][0],90.]
            elif SGData['SGUniq'] == 'c':
                cell = [parmDict['11'][0],parmDict['12'][0],parmDict['13'][0],90.,90.,parmDict['14'][0]]
        else:
            cell = [parmDict['11'][0],parmDict['12'][0],parmDict['13'][0],
                parmDict['14'][0],parmDict['15'][0],parmDict['16'][0]]
        return G2lat.cell2A(cell)
    except KeyError:
         return None

def UpdatePDFfit(Phase,RMCPdict):
    ''' Updates various PDFfit parameters held in GSAS-II
    '''

    General = Phase['General']
    if RMCPdict['refinement'] == 'normal':
        fName = General['Name']+'-PDFfit.rstr'
        try:
            rstr = open(fName.replace(' ','_'),'r')
        except FileNotFoundError:
            return [fName,'Not found - PDFfit failed']
        lines = rstr.readlines()
        rstr.close()
        header = [line[:-1].split(' ',1) for line in lines[:7]]
        resdict = dict(header)
        for item in ['sharp','cell']:
            resdict[item] = [float(val) for val in resdict[item].split(',')]
        General['Cell'][1:7] = resdict['cell']
        for inam,name in enumerate(['delta2','delta1','sratio']):
            RMCPdict[name][0] = float(resdict['sharp'][inam])
        if 'shape' in resdict:
            if 'sphere' in resdict['shape']:
                RMCPdict['spdiameter'][0] = float(resdict['shape'].split()[-1])
            else:
                RMCPdict['stepcut'][0] = float(resdict['shape'][-1])
        cx,ct,cs,ci = G2mth.getAtomPtrs(Phase)
        Atoms = Phase['Atoms']
        atmBeg = 0
        for line in lines:
            atmBeg += 1
            if 'atoms' in line:
                break
        for atom in Atoms:
            atstr = lines[atmBeg][:-1].split()
            Uiistr = lines[atmBeg+2][:-1].split()
            Uijstr = lines[atmBeg+4][:-1].split()
            atom[cx:cx+4] = [float(atstr[1]),float(atstr[2]),float(atstr[3]),float(atstr[4])]
            atom[ci] = 'A'
            atom[ci+2:ci+5] = [float(Uiistr[0]),float(Uiistr[1]),float(Uiistr[2])]
            atom[ci+5:ci+8] = [float(Uijstr[0]),float(Uijstr[1]),float(Uijstr[2])]
            atmBeg += 6
        fName = General['Name']+'-PDFfit.res'
    else:
        fName = 'Sequential_PDFfit.res'
    try:
        res = open(fName.replace(' ','_'),'r')
    except FileNotFoundError:
        return [fName,'Not found - PDFfit failed']
    lines = res.readlines()
    res.close()
    Ibeg = False
    resline = ''
    XNdata = {'Xdata':RMCPdict['Xdata'],'Ndata':RMCPdict['Ndata']}
    for line in lines:
        if 'Radiation' in line and 'X-Rays' in line:
            dkey = 'Xdata'
        if 'Radiation' in line and'Neutrons' in line:
            dkey = 'Ndata'
        if 'Qdamp' in line and '(' in line:
            XNdata[dkey]['qdamp'][0] = float(line.split()[4])
        if 'Qbroad' in line and '(' in line:
            XNdata[dkey]['qbroad'][0] = float(line.split()[4])
        if 'Scale' in line and '(' in line:
            XNdata[dkey]['dscale'][0] = float(line.split()[3])

    for iline,line in enumerate(lines):
        if 'Refinement parameters' in line:
            Ibeg = True
            continue
        if Ibeg:
            if '---------' in line:
                break
            resline += line[:-1]
    for iline,line in enumerate(lines):
        if 'Rw - ' in line:
            if 'nan' in line:
                Rwp = 100.0
            else:
                Rwp = float(line.split(':')[1])
    results = resline.replace('(','').split(')')[:-1]
    results = ['@'+result.lstrip() for result in results]
    results = [item.split() for item in results]
    RMCPdict['Parms'] = dict([[item[0][1:-1],float(item[1])] for item in results])      #{'n':val,...}
    if RMCPdict['refinement'] == 'normal':
        fName = General['Name']+'-PDFfit.py'
        py = open(fName.replace(' ','_'),'r')
        pylines = py.readlines()
        py.close()
        py = open(fName.replace(' ','_'),'w')
        newpy = []
        for pyline in pylines:
            if 'setpar' in pyline:
                parm = pyline.split('(')[1].split(',')[0]
                newpy.append('pf.setpar(%s,%.5f)\n'%(parm,RMCPdict['Parms'][parm]))
            else:
                newpy.append(pyline)
        py.writelines(newpy)
        py.close()
        RMCPdict.update(XNdata)
        results = dict([[item[0][:-1],float(item[1])] for item in results if item[0][:-1] in RMCPdict['AtomVar']])
        RMCPdict['AtomVar'].update(results)
        return None
    else:   #sequential
        newParms = dict([[item[0][1:-1],[float(item[1]),float(item[2])]] for item in results])  #{'n':[val,esd],...}
        return newParms,Rwp

def MakefullrmcSupercell(Phase,RMCPdict):
    '''Create a fullrmc supercell from GSAS-II

    :param dict Phase: phase information from data tree
    :param dict RMCPdict: fullrmc parameters from GUI
    :param list grpDict: a list of lists where the inner list
      contains the atom numbers contained in each group. e.g.
      [[0,1,2,3,4],[5,6],[4,6]] creates three groups with
      atoms 0-4 in the first
      atoms 5 & 6 in the second and
      atoms 4 & 6 in the third. Note that it is fine that
      atom 4 appears in two groups.
    '''
    #for i in (0,1): grpDict[i].append(1)    # debug: 1st & 2nd atoms in 2nd group
    cell = Phase['General']['Cell'][1:7]
    A,B = G2lat.cell2AB(cell)
    cx,ct,cs,cia = Phase['General']['AtomPtrs']
    SGData = Phase['General']['SGData']
    atomlist = []
    for i,atom in enumerate(Phase['Atoms']):
        el = ''.join([i for i in atom[ct] if i.isalpha()])
        grps = [j for j,g in enumerate(RMCPdict.get('Groups',[])) if i in g]
        atomlist.append((el, atom[ct-1], grps))
    # create a list of coordinates with symmetry & unit cell translation duplicates
    coordlist = []
    cellnum = -1
    for a in range(int(0.5-RMCPdict['SuperCell'][0]/2),int(1+RMCPdict['SuperCell'][0]/2)):
        for b in range(int(0.5-RMCPdict['SuperCell'][1]/2),int(1+RMCPdict['SuperCell'][1]/2)):
            for c in range(int(0.5-RMCPdict['SuperCell'][2]/2),int(1+RMCPdict['SuperCell'][2]/2)):
                cellnum += 1
                for i,atom in enumerate(Phase['Atoms']):
                    for item in G2spc.GenAtom(atom[cx:cx+3],SGData,Move=False):
#                        if i == 0: print(item[0]+[a,b,c])
                        xyzOrth = np.inner(A,item[0]+[a,b,c])
                        #coordlist.append((i,list(xyzOrth),cellnum,list(item[0]+[a,b,c])))
                        coordlist.append((item[1],cellnum,i,list(xyzOrth)))
    return atomlist,coordlist

def MakefullrmcRun(pName,Phase,RMCPdict):
    '''Creates a script to run fullrmc. Returns the name of the file that was
    created.
    '''
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
    # rmin = RMCPdict['min Contact']
    cell = Phase['General']['Cell'][1:7]
    SymOpList = G2spc.AllOps(Phase['General']['SGData'])[0]
    cx,ct,cs,cia = Phase['General']['AtomPtrs']
    atomsList = []
    for atom in Phase['Atoms']:
        el = ''.join([i for i in atom[ct] if i.isalpha()])
        atomsList.append([el] + atom[cx:cx+4])
    projDir,projName = os.path.split(os.path.abspath(pName))
    scrname = pName+'-fullrmc.py'
    #restart = '%s_restart.pdb'%pName
    Files = RMCPdict['files']
    rundata = ''
    rundata += '## fullrmc %s file ##\n## OK to edit this by hand ##\n'%scrname
    rundata += '# created in '+__file__+" v"+filversion
    rundata += dt.datetime.strftime(dt.datetime.now()," at %Y-%m-%dT%H:%M\n")
    rundata += '''
# fullrmc imports (all that are potentially useful)
import os,glob
import time
import pickle
import types
import copy
import numpy as np
import matplotlib as mpl
import fullrmc
from pdbparser import pdbparser
from pdbparser.Utilities.Database import __ATOM__
from fullrmc.Core import Collection
from fullrmc.Engine import Engine
import fullrmc.Constraints.PairDistributionConstraints as fPDF
from fullrmc.Constraints.StructureFactorConstraints import ReducedStructureFactorConstraint, StructureFactorConstraint
from fullrmc.Constraints.RadialDistributionConstraints import RadialDistributionConstraint
from fullrmc.Constraints.StructureFactorConstraints import NormalizedStructureFactorConstraint
from fullrmc.Constraints.DistanceConstraints import DistanceConstraint
from fullrmc.Constraints.BondConstraints import BondConstraint
from fullrmc.Constraints.AngleConstraints import BondsAngleConstraint
from fullrmc.Constraints.DihedralAngleConstraints import DihedralAngleConstraint
from fullrmc.Generators.Swaps import SwapPositionsGenerator
from fullrmc.Core.MoveGenerator import MoveGeneratorCollector
from fullrmc.Generators.Translations import TranslationGenerator
from fullrmc.Generators.Rotations import RotationGenerator

# utility routines
def writeHeader(ENGINE,statFP):
    """header for stats file"""
    statFP.write('generated-steps, total-error, ')
    for c in ENGINE.constraints:
        statFP.write(c.constraintName)
        statFP.write(', ')
    statFP.write('\\n')
    statFP.flush()

def writeCurrentStatus(ENGINE,statFP,plotF):
    """line in stats file & current constraint plots"""
    statFP.write(str(ENGINE.generated))
    statFP.write(', ')
    statFP.write(str(ENGINE.totalStandardError))
    statFP.write(', ')
    for c in ENGINE.constraints:
        statFP.write(str(c.standardError))
        statFP.write(', ')
    statFP.write('\\n')
    statFP.flush()
    mpl.use('agg')
    fp = open(plotF,'wb')
    for c in ENGINE.constraints:
        p = c.plot(show=False)
        p[0].canvas.draw()
        image = p[0].canvas.buffer_rgba()
        pickle.dump(c.constraintName,fp)
        pickle.dump(np.array(image),fp)
    fp.close()

def calcRmax(ENGINE):
    """from Bachir, works for non-orthorhombic cells"""
    a,b,c = ENGINE.basisVectors
    lens = []
    ts    = np.linalg.norm(np.cross(a,b))/2
    lens.extend( [ts/np.linalg.norm(a), ts/np.linalg.norm(b)] )
    ts = np.linalg.norm(np.cross(b,c))/2
    lens.extend( [ts/np.linalg.norm(b), ts/np.linalg.norm(c)] )
    ts = np.linalg.norm(np.cross(a,c))/2
    lens.extend( [ts/np.linalg.norm(a), ts/np.linalg.norm(c)] )
    return min(lens)
'''
    if RMCPdict.get('Groups',[]): rundata += '''
def makepdb(atoms, coords, bbox=None):
    """creates a supercell directly from atom info"""
    # used when ENGINE.build_crystal_set_pdb is not called
    prevcell = None
    rec = copy.copy(__ATOM__)
    rec['residue_name'] = 'MOL'
    records = []
    seqNum  = 0
    segId   = '0'
    groups = {}
    for symcell in set([(sym,cell) for sym,cell,atm,xyz in coords]):
        seqNum += 1
        if seqNum == 9999:
            seqNum = 1
            segId  = str(int(segId) + 1)
        for i,(sym,cell,atm,(x,y,z)) in enumerate(coords):
            if (sym,cell) != symcell: continue
            rec   = copy.copy(rec)
            for grp in atoms[atm][2]:
                if (sym,cell) not in groups:
                    groups[(sym,cell)] = {}
                if grp not in groups[(sym,cell)]:
                    groups[(sym,cell)][grp] = [len(records)]
                else:
                    groups[(sym,cell)][grp].append(len(records))
            rec['coordinates_x']      = x
            rec['coordinates_y']      = y
            rec['coordinates_z']      = z
            rec['element_symbol']     = atoms[atm][0]
            rec['atom_name']          = atoms[atm][1]
            rec['sequence_number']    = seqNum
            rec['segment_identifier'] = segId
            records.append(rec)
    # create pdb
    pdb = pdbparser()
    pdb.records = records
    if groups:
        return pdb,[groups[j][i] for j in groups for i in groups[j]]
    else:
        return pdb,[]
'''
    rundata += '''
### When True, erases an existing engine to provide a fresh start
FRESH_START = {:}
dirName = "{:}"
prefix = "{:}"
project = prefix + "-fullrmc"
time0 = time.time()
'''.format(RMCPdict['ReStart'][0],projDir,projName)

    rundata += '# setup structure\n'
    rundata += 'cell = ' + str(cell) + '\n'
    rundata += 'supercell = ' + str(RMCPdict['SuperCell']) + '\n'
    rundata += '\n# define structure info\n'
    if RMCPdict.get('Groups',[]):
        # compute bounding box coordinates
        bbox = []
        A,B = G2lat.cell2AB(cell)
        for i in range(3):
            for val in int(0.5-RMCPdict['SuperCell'][i]/2),int(1+RMCPdict['SuperCell'][0]/2):
                fpos = [0,0,0]
                fpos[i] = val
                bbox.append(np.inner(A,fpos))
        rundata += 'bboxlist = [     # orthogonal coordinate for supercell corners\n'
        for i in bbox:
            rundata += '  '+str(list(i))+',\n'
        rundata += ' ] # bboxlist\n\n'
        atomlist,coordlist = MakefullrmcSupercell(Phase,RMCPdict)
        rundata += 'atomlist = [  # [element, label, grouplist]\n'
        for i in atomlist:
            rundata += '  '+str(i)+',\n'
        rundata += ' ] # atomlist\n\n'
        rundata += 'coordlist = [     # (sym#, cell#, atom#, [ortho coords],)\n'
        for i in coordlist:
            rundata += '  '+str(i)+',\n'
        rundata += ' ] # coordlist\n'
    else:
        rundata += "SymOpList = "+str([i.lower() for i in SymOpList]) + '\n'
        rundata += 'atomList = ' + str(atomsList).replace('],','],\n  ') + '\n'

    rundata += '\n# initialize engine\n'
    rundata += '''
engineFileName = os.path.join(dirName, project + '.rmc')
projectStats = os.path.join(dirName, project + '.stats')
projectPlots = os.path.join(dirName, project + '.plots')
projectXYZ = os.path.join(dirName, project + '.atoms')
pdbFile = os.path.join(dirName, project + '_restart.pdb')
# check Engine exists if so (and not FRESH_START) load it otherwise build it
ENGINE = Engine(path=None)
if not ENGINE.is_engine(engineFileName) or FRESH_START:
    ENGINE = Engine(path=engineFileName, freshStart=True)
'''
    if RMCPdict.get('Groups',[]):
        rundata += '''
    # create structure from GSAS-II constructed supercell
    bbox = (np.array(bboxlist[1::2])-np.array(bboxlist[0::2])).flatten()
    pdb,grouplist = makepdb(atomlist,coordlist,bbox)
    ENGINE.set_pdb(pdb)
    ENGINE.set_boundary_conditions(bbox)
    if grouplist: ENGINE.set_groups(grouplist)
'''
        if RMCPdict.get('GroupMode',0) == 0:   # 'Rotate & Translate'
            rundata += '''
    for g in ENGINE.groups:
        TMG = TranslationGenerator(amplitude=0.2) # create translation generator
        if len(g) > 1:  # create rotation generator for groups with more than 1 atom
            RMG = RotationGenerator(amplitude=2)
            MG  = MoveGeneratorCollector(collection=[TMG,RMG],randomize=True)
        else:
            MG  = MoveGeneratorCollector(collection=[TMG],randomize=True)
        g.set_move_generator( MG )
'''
        elif RMCPdict.get('GroupMode',0) == 1: # 'Rotate only'
            rundata += '''
    for g in ENGINE.groups:
        if len(g) > 1:  # create rotation generator for groups with more than 1 atom
            RMG = RotationGenerator(amplitude=2)
            g.set_move_generator( RMG )
'''
        else:                                  # 'Translate only'
            rundata += '        # translate only set by default'
    else:
        rundata += '''
    # create structure, let fullrmc construct supercell
    ENGINE.build_crystal_set_pdb(symOps     = SymOpList,
                                 atoms      = atomList,
                                 unitcellBC = cell,
                                 supercell  = supercell)
    ENGINE.set_groups_as_atoms()
'''
    rundata += '    rho0 = len(ENGINE.allNames)/ENGINE.volume\n'
    rundata += '\n    # "Constraints" (includes experimental data) setup\n'
    # settings that require a new Engine
    for File in Files:
        filDat = RMCPdict['files'][File]
        if not os.path.exists(filDat[0]): continue
        sfwt = 'neutronCohb'
        if 'Xray' in File:
            sfwt = 'atomicNumber'
        if 'G(r)' in File:
            rundata += '    GR = np.loadtxt(os.path.join(dirName,"%s")).T\n'%filDat[0]
            if filDat[3] == 0:
                #rundata += '''    # read and xform G(r) as defined in RMCProfile
    # see eq. 44 in Keen, J. Appl. Cryst. (2001) 34 172-177\n'''
                #rundata += '    GR[1] *= 4 * np.pi * GR[0] * rho0 / sumCiBi2\n'
                #rundata += '    GofR = fPDF.PairDistributionConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
                rundata += '    # G(r) as defined in RMCProfile\n'
                rundata += '    GofR = RadialDistributionConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            elif filDat[3] == 1:
                rundata += '    # This is G(r) as defined in PDFFIT\n'
                rundata += '    GofR = fPDF.PairDistributionConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            elif filDat[3] == 2:
                rundata += '    # This is g(r)\n'
                rundata += '    GofR = fPDF.PairCorrelationConstraint(experimentalData=GR.T, weighting="%s")\n'%sfwt
            else:
                raise ValueError('Invalid G(r) type: '+str(filDat[3]))
            rundata += '    ENGINE.add_constraints([GofR])\n'
            rundata += '    GofR.set_limits((None, calcRmax(ENGINE)))\n'
            if RMCPdict['addThermalBroadening']:
                rundata += "    GofR.set_thermal_corrections({'defaultFactor': 0.001})\n"
                rundata += "    GofR.thermalCorrections['factors'] = {\n"
                RMCPdict['addThermalBroadening']
                for atm1 in RMCPdict['aTypes']:
                    for atm2 in RMCPdict['aTypes']:
                        rundata += "         ('{}', '{}'): {},\n".format(
                            atm1,atm2,
                            (RMCPdict['ThermalU'].get(atm1,0.005)+                                           RMCPdict['ThermalU'].get(atm2,0.005))/2)
                rundata += '    }\n'
        elif '(Q)' in File:
            rundata += '    SOQ = np.loadtxt(os.path.join(dirName,"%s")).T\n'%filDat[0]
            if filDat[3] == 0:
                rundata += '    # F(Q) as defined in RMCProfile\n'
                #rundata += '    SOQ[1] *= 1 / sumCiBi2\n'
                if filDat[4]:
                    rundata += '    SOQ[1] = Collection.sinc_convolution(q=SOQ[0],sq=SOQ[1],rmax=calcRmax(ENGINE))\n'
                rundata += '    SofQ = NormalizedStructureFactorConstraint(experimentalData=SOQ.T, weighting="%s")\n'%sfwt
            elif filDat[3] == 1:
                rundata += '    # S(Q) as defined in PDFFIT\n'
                rundata += '    SOQ[1] -= 1\n'
                if filDat[4]:
                    rundata += '    SOQ[1] = Collection.sinc_convolution(q=SOQ[0],sq=SOQ[1],rmax=calcRmax(ENGINE))\n'
                rundata += '    SofQ = ReducedStructureFactorConstraint(experimentalData=SOQ.T, weighting="%s")\n'%sfwt
            else:
                raise ValueError('Invalid S(Q) type: '+str(filDat[3]))
            rundata += '    ENGINE.add_constraints([SofQ])\n'
        else:
            print('What is this?')
    minDists = ''
    if BondList and RMCPdict.get('useBondConstraints',True):
        rundata += '''    B_CONSTRAINT   = BondConstraint()
    ENGINE.add_constraints(B_CONSTRAINT)
    B_CONSTRAINT.create_supercell_bonds(bondsDefinition=[
'''
        for pair in BondList:
            e1,e2 = pair.split('-')
            d1,d2 = BondList[pair]
            if d1 == 0: continue
            if d2 == 0:
                minDists += '("element","{}","{}",{}),'.format(e1.strip(),e2.strip(),d1)
            else:
                rundata += '            ("element","{}","{}",{},{}),\n'.format(
                                        e1.strip(),e2.strip(),d1,d2)
        rundata += '             ])\n'
    rundata += '    D_CONSTRAINT = DistanceConstraint(defaultLowerDistance={})\n'.format(RMCPdict['min Contact'])
    if minDists:
        rundata += "    D_CONSTRAINT.set_pairs_definition( {'inter':[" + minDists + "]})\n"
    rundata += '    ENGINE.add_constraints(D_CONSTRAINT)\n'

    if AngleList:
        rundata += '''    A_CONSTRAINT   = BondsAngleConstraint()
    ENGINE.add_constraints(A_CONSTRAINT)
    A_CONSTRAINT.create_supercell_angles(anglesDefinition=[
'''
        for item in AngleList:
            rundata += ('            '+
               '("element","{1}","{0}","{2}",{5},{6},{5},{6},{3},{4}),\n'.format(*item))
        rundata += '             ])\n'
    rundata += '''
    for f in glob.glob(os.path.join(dirName,prefix+"_*.log")): os.remove(f)
    ENGINE.save()
else:
    ENGINE = ENGINE.load(path=engineFileName)

ENGINE.set_log_file(os.path.join(dirName,prefix))
'''
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
    rundata += '''for c in ENGINE.constraints:
    if hasattr(c, '_ExperimentalConstraint__adjustScaleFactor'):
        def _constraint_copy_needs_lut(self, *args, **kwargs):
            result =  fPDF.PairDistributionConstraint._constraint_copy_needs_lut(self, *args, **kwargs)
            result['_ExperimentalConstraint__adjustScaleFactor'] = '_ExperimentalConstraint__adjustScaleFactor'
            return result
        c._constraint_copy_needs_lut = types.MethodType(_constraint_copy_needs_lut, c)
    if c.__class__.__name__ in ('ReducedStructureFactorConstraint', 'StructureFactorConstraint'):
        def _constraint_copy_needs_lut(self, *args, **kwargs):
            result =  fPDF.PairDistributionConstraint._constraint_copy_needs_lut(self, *args, **kwargs)
            result.pop('_PairDistributionConstraint__histogramSize', None)
            result.pop('_PairDistributionConstraint__shellVolumes', None)
            result.pop('_PairDistributionConstraint__shellCenters', None)
            result.pop('_PairDistributionConstraint__windowArray', None)
            result.pop('_PairDistributionConstraint__experimentalDistances', None)
            result.pop('_PairDistributionConstraint__experimentalPDF', None)
            result.pop('_PairDistributionConstraint__minimumDistance', None)
            result.pop('_PairDistributionConstraint__maximumDistance', None)
            result.pop('_PairDistributionConstraint__distanceBin', None)
            result.pop('_shapeFuncParams', None)
            result.pop('_shapeArray', None)
            result['_StructureFactorConstraint__Gr2SqMatrix'] = '_StructureFactorConstraint__Gr2SqMatrix'
            result['_StructureFactorConstraint__histogramSize'] = '_StructureFactorConstraint__histogramSize'
            result['_StructureFactorConstraint__shellVolumes'] = '_StructureFactorConstraint__shellVolumes'
            result['_StructureFactorConstraint__shellCenters'] = '_StructureFactorConstraint__shellCenters'
            result['_StructureFactorConstraint__windowArray'] = '_StructureFactorConstraint__windowArray'
            result['_StructureFactorConstraint__experimentalQValues'] = '_StructureFactorConstraint__experimentalQValues'
            result['_StructureFactorConstraint__experimentalSF'] = '_StructureFactorConstraint__experimentalSF'
            result['_StructureFactorConstraint__minimumDistance'] = '_StructureFactorConstraint__minimumDistance'
            result['_StructureFactorConstraint__maximumDistance'] = '_StructureFactorConstraint__maximumDistance'
            result['_StructureFactorConstraint__distanceBin'] = '_StructureFactorConstraint__distanceBin'
            return result
        c._constraint_copy_needs_lut = types.MethodType(_constraint_copy_needs_lut, c)

'''
    # rundata += '\n# set weights -- do this now so values can be changed without a restart\n'
    # rundata += 'wtDict = {}\n'
    # for File in Files:
    #     filDat = RMCPdict['files'][File]
    #     if not os.path.exists(filDat[0]): continue
    #     if 'Xray' in File:
    #         sfwt = 'atomicNumber'
    #     else:
    #         sfwt = 'neutronCohb'
    #     if 'G(r)' in File:
    #         typ = 'Pair'
    #     elif '(Q)' in File:
    #         typ = 'Struct'
    #     rundata += 'wtDict["{}-{}"] = {}\n'.format(typ,sfwt,filDat[1])
    rundata += '\n# set PDF fitting range\n'
    rundata += 'for c in ENGINE.constraints:  # loop over predefined constraints\n'
    rundata += '    if type(c) is fPDF.PairDistributionConstraint:\n'
    # rundata += '        c.set_variance_squared(1./wtDict["Pair-"+c.weighting])\n'
    rundata += '        c.set_limits((None,calcRmax(ENGINE)))\n'
    if RMCPdict['FitScale']:
        rundata += '        c.set_adjust_scale_factor((10, 0.01, 100.))\n'
    # rundata += '        c.set_variance_squared(1./wtDict["Struct-"+c.weighting])\n'
    if RMCPdict['FitScale']:
        rundata += '    elif type(c) is ReducedStructureFactorConstraint:\n'
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
    rundata += '''
if FRESH_START:
    # initialize engine with one step to get starting config energetics
    ENGINE.run(restartPdb=pdbFile,numberOfSteps=1, saveFrequency=1)
    statFP = open(projectStats,'w')
    writeHeader(ENGINE,statFP)
    writeCurrentStatus(ENGINE,statFP,projectPlots)
else:
    statFP = open(projectStats,'a')

# setup runs for fullrmc
'''
    rundata += 'steps = {}\n'.format(RMCPdict['Steps/cycle'])
    rundata += 'for _ in range({}):\n'.format(RMCPdict['Cycles'])
    rundata += '    expected = ENGINE.generated+steps\n'

    rundata += '    ENGINE.run(restartPdb=pdbFile,numberOfSteps=steps, saveFrequency=steps)\n'
    rundata += '    writeCurrentStatus(ENGINE,statFP,projectPlots)\n'
    rundata += '    if ENGINE.generated != expected: break # run was stopped'
    rundata += '''
statFP.close()
fp = open(projectXYZ,'w')  # save final atom positions
fp.write('cell: {} {} {} {} {} {}\\n')
fp.write('supercell: {} {} {}\\n')
'''.format(*cell,*RMCPdict['SuperCell'])
    rundata += '''# loop over atoms
for n,e,(x,y,z) in zip(ENGINE.allNames,
                       ENGINE.allElements,ENGINE.realCoordinates):
    fp.write('{} {} {:.5f} {:.5f} {:.5f}\\n'.format(n,e,x,y,z))
fp.close()
print("ENGINE run time %.2f s"%(time.time()-time0))
'''
    rfile = open(scrname,'w')
    rfile.writelines(rundata)
    rfile.close()
    return scrname

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

def ISO2PDFfit(Phase):
    ''' Creates new phase structure to be used for PDFfit from an ISODISTORT mode displacement phase.
    It builds the distortion mode parameters to be used as PDFfit variables for atom displacements from
    the original parent positions as transformed to the child cell wiht symmetry defined from ISODISTORT.

    :param Phase: dict GSAS-II Phase structure; must contain ISODISTORT dict. NB: not accessed otherwise

    :returns: dict: GSAS-II Phase structure; will contain ['RMC']['PDFfit'] dict
    '''

    Trans = np.eye(3)
    Uvec = np.zeros(3)
    Vvec = np.zeros(3)
    Phase = copy.deepcopy(Phase)
    Atoms = Phase['Atoms']
    parentXYZ = Phase['ISODISTORT']['G2parentCoords']           #starting point for mode displacements
    cx,ct,cs,cia = Phase['General']['AtomPtrs']
    for iat,atom in enumerate(Atoms):
        atom[cx:cx+3] = parentXYZ[iat]
    SGData = copy.deepcopy(Phase['General']['SGData'])
    SGOps = SGData['SGOps']
    newPhase = copy.deepcopy(Phase)
    newPhase['ranId'] = rand.randint(0,sys.maxsize)
    newPhase['General']['Name'] += '_PDFfit'
    newPhase['General']['SGData'] = G2spc.SpcGroup('P 1')[1]    #this is for filled unit cell
    newPhase,atCodes = G2lat.TransformPhase(Phase,newPhase,Trans,Uvec,Vvec,False)
    newPhase['Histograms'] = {}
    newPhase['Drawing'] = []
    Atoms = newPhase['Atoms']
    RMCPdict = newPhase['RMC']['PDFfit']
    ISOdict = newPhase['ISODISTORT']
    RMCPdict['AtomConstr'] = []
    RMCPdict['SGData'] = copy.deepcopy(SGData)      #this is from the ISODISTORT child; defines PDFfit constraints
    Norms = ISOdict['NormList']
    ModeMatrix = ISOdict['Mode2VarMatrix']
    RMCPdict['AtomVar'] = {'@%d'%(itm+21):val for itm,val in enumerate(ISOdict['modeDispl'])}
    for iatm,[atom,atcode] in enumerate(zip(Atoms,atCodes)):
        pid,opid = [int(item) for item in atcode.split(':')]
        atmConstr = [atom[ct-1],atom[ct],'','','','','',atcode]
        for ip,pname in enumerate(['%s_d%s'%(atom[ct-1],x) for x in ['x','y','z']]):
            try:
                conStr = ''
                if Atoms[iatm][cx+ip]:
                    conStr += '%.5f'%Atoms[iatm][cx+ip]
                pid = ISOdict['IsoVarList'].index(pname)
                consVec = ModeMatrix[pid]
                for ic,citm in enumerate(consVec):      #NB: this assumes orthorhombic or lower symmetry
                    if opid < 0:
                        citm *= -SGOps[100-opid%100-1][0][ip][ip]   #remove centering, if any
                    else:
                        citm *= SGOps[opid%100-1][0][ip][ip]
                    if citm > 0.:
                        conStr += '+%.5f*@%d'%(citm*Norms[ic],ic+21)
                    elif citm < 0.:
                        conStr += '%.5f*@%d'%(citm*Norms[ic],ic+21)
                atmConstr[ip+2] = conStr
            except ValueError:
                atmConstr[ip+2] = ''
        RMCPdict['AtomConstr'].append(atmConstr)
    return newPhase

def GetAtmDispList(ISOdata):
    atmDispList = []
    MT = ISOdata['Mode2VarMatrix'].T
    DispList = ISOdata['IsoVarList']
    N = len(DispList)
    for I in range(N):
        vary = []
        for i in range(N):
            if MT[I,i]:
                vary.append(DispList[i])
        atmDispList.append(vary)
    return atmDispList

#### Reflectometry calculations ################################################################################
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


#### Stacking fault simulation codes ################################################################################
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
    from . import atmdata
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
    from . import atmdata
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

#### Maximum Entropy Method - Dysnomia ###############################################################################
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
    fba.close()
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

#===Laue Fringe code ===================================================================
from . import NIST_profile as FP

class profileObj(FP.FP_profile):
    def conv_Lauefringe(self):
        """Compute the FT of the Laue Fringe function"""

        me=self.get_function_name() #the name of this convolver,as a string
        wave = self.param_dicts['conv_global']['dominant_wavelength']*1.e10 # in A
        pos = np.rad2deg(self.param_dicts["conv_global"]["twotheta0"])  # peak position as 2theta in deg
        posQ = np.pi * 4 * np.sin(self.param_dicts["conv_global"]["twotheta0"]/2) / wave # peak position as Q
        ttwid = self.twotheta_window_fullwidth_deg
        ncell = self.param_dicts[me]['Ncells']
        co2 = self.param_dicts[me]['clat'] / 2.
        dampM =  self.param_dicts[me]['dampM']
        dampP =  self.param_dicts[me]['dampP']
        fpowM =  self.param_dicts[me]['fitPowerM']
        fpowP =  self.param_dicts[me]['fitPowerP']
        ttlist = np.linspace(pos-ttwid/2,pos+ttwid/2,len(self._epsb2))
        Qs = np.pi * 4 * np.sin(np.deg2rad(ttlist/2)) / wave
        w =  np.exp(-1*10**((dampM) * np.abs(Qs - posQ)**fpowM))
        w2 = np.exp(-1*10**((dampP) * np.abs(Qs - posQ)**fpowP))
        w[len(w)//2:] = w2[len(w)//2:]
        weqdiv = w * np.sin(Qs * ncell * co2)**2 / (np.sin(Qs * co2)**2)
        weqdiv[:np.searchsorted(Qs,posQ - np.pi/self.param_dicts[me]['clat'])] = 0  # isolate central peak, if needed
        weqdiv[np.searchsorted(Qs,posQ + np.pi/self.param_dicts[me]['clat']):] = 0
        conv = FP.best_rfft(weqdiv)
        conv[1::2] *= -1 #flip center
        return conv

    def conv_Lorentzian(self):
        """Compute the FT of a Lorentz function where gamma is the FWHM"""
        ttwid = self.twotheta_window_fullwidth_deg
        me=self.get_function_name() #the name of this convolver,as a string
        g2gam = self.param_dicts[me]['g2gam'] # gsas-ii gamma in centidegrees
        gamma = g2gam/100 # deg
        ttlist = np.linspace(-ttwid/2,ttwid/2,len(self._epsb2))
        eqdiv = (0.5 * gamma / np.pi) / (gamma**2/4. + ttlist**2)
        conv = FP.best_rfft(eqdiv)
        conv[1::2] *= -1 #flip center
        return conv

    def conv_Gaussian(self):
        """Compute the FT of a Gaussian where sigma**2 is the variance"""
        ttwid = self.twotheta_window_fullwidth_deg
        me=self.get_function_name() #the name of this convolver,as a string
        g2sig2 = self.param_dicts[me]['g2sig2'] # gsas-ii sigma**2 in centidegr**2
        sigma = math.sqrt(g2sig2)/100.
        ttlist = np.linspace(-ttwid/2,ttwid/2,len(self._epsb2))
        eqdiv = np.exp(-0.5*ttlist**2/sigma**2) / math.sqrt(2*np.pi*sigma**2)
        conv = FP.best_rfft(eqdiv)
        conv[1::2] *= -1 #flip center
        return conv

def LaueFringePeakCalc(ttArr,intArr,lam,peakpos,intens,sigma2,gamma,shol,ncells,clat,dampM,dampP,calcwid,fitPowerM=2,fitPowerP=2,plot=False):
    '''Compute the peakshape for a Laue Fringe peak convoluted with a Gaussian, Lorentzian &
    an axial divergence asymmetry correction.

    :param np.array ttArr: Array of two-theta values (in degrees)
    :param np.array intArr: Array of intensity values (peaks are added to this)
    :param float lam: wavelength in Angstrom
    :param float peakpos: peak position in two-theta (deg.)
    :param float intens: intensity factor for peak
    :param float sigma2: Gaussian variance (in centidegrees**2) **
    :param float gamma:  Lorenzian FWHM (in centidegrees) **
    :param float shol: FCJ (S + H)/L where S=sample-half height, H=slit half-height, L=radius **
    :param float ncells: number of unit cells in specular direction **
    :param float clat: c lattice parameter **
    :param float dampM:
    :param float dampP:
    :param float calcwid: two-theta (deg.) width for cutoff of peak computation.
       Defaults to 5
    :param float fitPowerM: exponent used for damping fall-off on minus side of peak
    :param float fitPowerP: exponent used for damping fall-off on plus side of peak
    :param bool plot: for debugging, shows contributions to peak

    **  If term is <= zero, item is removed from convolution
    '''
    # def LaueFringePeakPlot(ttArr,intArr):
    #     import matplotlib.pyplot as plt
    #     refColors = ['xkcd:blue','xkcd:red','xkcd:green','xkcd:cyan','xkcd:magenta','xkcd:black',
    #                      'xkcd:pink','xkcd:brown','xkcd:teal','xkcd:orange','xkcd:grey','xkcd:violet',]
    #     fig, ax = plt.subplots()
    #     ax.set(title='Peak convolution functions @ 2theta={:.3f}'.format(peakpos),
    #                xlabel=r'$\Delta 2\theta, deg$',
    #                ylabel=r'Intensity (arbitrary)')
    #     ax.set_yscale("log",nonpositive='mask')
    #     ttmin = ttmax = 0
    #     for i,conv in enumerate(convList):
    #         f = NISTpk.convolver_funcs[conv]()
    #         if f is None: continue
    #         FFT = FP.best_irfft(f)
    #         if f[1].real > 0: FFT = np.roll(FFT,int(len(FFT)/2.))
    #         FFT /= FFT.max()
    #         if i == 0:
    #             tt = np.linspace(-NISTpk.twotheta_window_fullwidth_deg/2,
    #                                 NISTpk.twotheta_window_fullwidth_deg/2,len(FFT))
    #         ttmin = min(ttmin,tt[np.argmax(FFT>.005)])
    #         ttmax = max(ttmax,tt[::-1][np.argmax(FFT[::-1]>.005)])
    #         color = refColors[i%len(refColors)]
    #         ax.plot(tt,FFT,color,label=conv[5:])
    #     color = refColors[(i+1)%len(refColors)]
    #     ax.plot(ttArr-peakpos,intArr/max(intArr),color,label='Convolution')
    #     ax.set_xlim((ttmin,ttmax))
    #     ax.legend(loc='best')
    #     plt.show()
    # hardcoded constants
    diffRadius = 220 # diffractometer radius in mm; needed for axial divergence, etc, but should not matter
    axial_factor = 1.5  # fudge factor to bring sh/l broadening to ~ agree with FPA
    equatorial_divergence_deg = 0.5 # not sure exactly what this impacts
    NISTparms = {
        "": {
            'equatorial_divergence_deg' : equatorial_divergence_deg,
            'dominant_wavelength' : 1.e-10 * lam,
            'diffractometer_radius' : 1e-3* diffRadius, # diffractometer radius in m
            'oversampling' : 8,
            },
        "emission": {
            'emiss_wavelengths' : 1.e-10 * np.array([lam]),
            'emiss_intensities' : np.array([1.]),
            'emiss_gauss_widths' : 1.e-10 * 1.e-3 * np.array([0.001]),
            'emiss_lor_widths' : 1.e-10 * 1.e-3 * np.array([0.001]),
            'crystallite_size_gauss' : 1.e-9 * 1e6,
            'crystallite_size_lor' : 1.e-9 * 1e6
            },
        "axial": {
            'axDiv':"full",
            'slit_length_source' : 1e-3 * diffRadius * shol * axial_factor,
            'slit_length_target' : 1e-3 * diffRadius * shol * 1.00001 * axial_factor, # != 'slit_length_source'
            'length_sample' : 1e-3 * diffRadius * shol * axial_factor,
            'n_integral_points' : 10,
            'angI_deg' : 2.5,
            'angD_deg': 2.5,
            },
        'Gaussian': {'g2sig2': sigma2},
        'Lorentzian': {'g2gam': gamma},
        'Lauefringe': {'Ncells': ncells, 'clat':clat, 'dampM': dampM,
                        'dampP': dampP, 'fitPowerM':fitPowerM, 'fitPowerP':fitPowerP},
        }
    NISTpk=profileObj(anglemode="twotheta",
                    output_gaussian_smoother_bins_sigma=1.0,
                    oversampling=NISTparms.get('oversampling',10))
    NISTpk.debug_cache=False
    for key in NISTparms: #set parameters for each convolver
        if key:
            NISTpk.set_parameters(convolver=key,**NISTparms[key])
        else:
            NISTpk.set_parameters(**NISTparms[key])
    # find closest point to peak location (which may be outside limits of the array)
    center_bin_idx=min(ttArr.searchsorted(peakpos),len(ttArr)-1)
    step = (ttArr[-1]-ttArr[0])/(len(ttArr)-1)
    NISTpk.set_optimized_window(twotheta_exact_bin_spacing_deg=step,
                    twotheta_window_center_deg=ttArr[center_bin_idx],
                    twotheta_approx_window_fullwidth_deg=calcwid,
                    )
    NISTpk.set_parameters(twotheta0_deg=peakpos)
    convList = ['conv_emission']
    if ncells: convList += ['conv_Lauefringe']
    if sigma2 > 0: convList += ['conv_Gaussian']
    if gamma > 0: convList += ['conv_Lorentzian']
    if shol > 0: convList += ['conv_axial']

    # global deriv
    # if deriv:
    #     peakObj = NISTpk.compute_line_profile(convolver_names=convList,compute_derivative=True)
    # else:
    #     peakObj = NISTpk.compute_line_profile(convolver_names=convList)
    peakObj = NISTpk.compute_line_profile(convolver_names=convList)

    pkPts = len(peakObj.peak)
    pkMax = peakObj.peak.max()
    startInd = center_bin_idx-(pkPts//2)
    istart = None
    pstart = None
    iend = None
    pend = None
    # adjust data range if peak calc begins below data range or ends above data range
    # but range of peak calc should not extend past both ends of ttArr
    if startInd < 0:
        iend = startInd+pkPts
        pstart = -startInd
    elif startInd > len(intArr):
        return
#    elif startInd+pkPts >= len(intArr):
    elif startInd+pkPts > len(intArr):
        offset = pkPts - len( intArr[startInd:] )
        istart = startInd
        iend = startInd+pkPts-offset
        pend = -offset
    else:
        istart = startInd
        iend = startInd+pkPts
    intArr[istart:iend] += intens * peakObj.peak[pstart:pend]/pkMax
#    if plot:
#        LaueFringePeakPlot(ttArr[istart:iend], (intens * peakObj.peak[pstart:pend]/pkMax))

def LaueSatellite(peakpos,wave,c,ncell,j=[-4,-3,-2,-1,0,1,2,3,4]):
    '''Returns the locations of the Laue satellite positions relative
    to the peak position

    :param float peakpos: the peak position in degrees 2theta
    :param float ncell:   Laue fringe parameter, number of unit cells in layer
    :param list j: the satellite order, where j=-1 is the first satellite
      on the lower 2theta side and j=1 is the first satellite on the high
      2theta side. j=0 gives the peak position
    '''
    Qpos = 4 * np.pi * np.sin(peakpos * np.pi / 360) / wave
    dQvals = (2 * np.array(j) + np.sign(j)) * np.pi / (c * ncell)
    return np.arcsin((Qpos+dQvals)*wave/(4*np.pi)) * (360 / np.pi)

def SetDefaultSubstances():
    'Fills in default items for the SASD Substances dictionary'
    return {'Substances':{'vacuum':{'Elements':{},'Volume':1.0,'Density':0.0,'Scatt density':0.0,'XImag density':0.0},
        'unit scatter':{'Elements':None,'Volume':None,'Density':None,'Scatt density':1.0,'XImag density':1.0}}}

def SetDefaultSASDModel():
    'Fills in default items for the SASD Models dictionary'
    return {'Back':[0.0,False],
        'Size':{'MinDiam':50,'MaxDiam':10000,'Nbins':100,'logBins':True,'Method':'MaxEnt',
                'Distribution':[],'Shape':['Spheroid',1.0],
                'MaxEnt':{'Niter':100,'Precision':0.01,'Sky':-3},
                'IPG':{'Niter':100,'Approach':0.8,'Power':-1},'Reg':{},},
        'Pair':{'Method':'Moore','MaxRadius':100.,'NBins':100,'Errors':'User',
                'Percent error':2.5,'Background':[0,False],'Distribution':[],
                'Moore':10,'Dist G':100.,'Result':[],},
        'Particle':{'Matrix':{'Name':'vacuum','VolFrac':[0.0,False]},'Levels':[],},
        'Shapes':{'outName':'run','NumAA':100,'Niter':1,'AAscale':1.0,'Symm':1,'bias-z':0.0,
                 'inflateV':1.0,'AAglue':0.0,'pdbOut':False,'boxStep':4.0},
        'Current':'Size dist.','BackFile':'',
        }

def SetDefaultREFDModel():
    '''Fills in default items for the REFD Models dictionary which are
    defined as follows for each layer:

    * Name: name of substance
    * Thick: thickness of layer in Angstroms (not present for top & bottom layers)
    * Rough: upper surface roughness for layer (not present for toplayer)
    * Penetration: mixing of layer substance into layer above-is this needed?
    * DenMul: multiplier for layer scattering density (default = 1.0)

    Top layer defaults to vacuum (or air/any gas); can be substituted for some other substance.

    Bottom layer default: infinitely thisck Silicon; can be substituted for some other substance.
    '''
    return {'Layers':[{'Name':'vacuum','DenMul':[1.0,False],},                                  #top layer
        {'Name':'vacuum','Rough':[0.,False],'Penetration':[0.,False],'DenMul':[1.0,False],}],   #bottom layer
        'Scale':[1.0,False],'FltBack':[0.0,False],'Zero':'Top','dQ type':'None','Layer Seq':[],               #globals
        'Minimizer':'LMLS','Resolution':[0.,'Const dq/q'],'Recomb':0.5,'Toler':0.5,             #minimizer controls
        'DualFitFiles':['',],'DualFltBacks':[[0.0,False],],'DualScales':[[1.0,False],]}         #optional stuff for multidat fits?

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
    gplot.plot(xdata,getPeakProfile(parmDict0,xdata,np.zeros_like(xdata),varyList,bakType))
    fplot = plotter.add('FCJ-Voigt, Ka1+2').gca()
    fplot.plot(xdata,getBackground('',parmDict1,bakType,'PXC',xdata)[0])
    fplot.plot(xdata,getPeakProfile(parmDict1,xdata,np.zeros_like(xdata),varyList,bakType))

def test1():
    if NeedTestData: TestData()
    time0 = time.time()
    for i in range(100):
        getPeakProfile(parmDict1,xdata,np.zeros_like(xdata),varyList,bakType)
    G2fil.G2Print ('100+6*Ka1-2 peaks=1200 peaks %.2f'%time.time()-time0)

def test2(name,delt):
    if NeedTestData: TestData()
    varyList = [name,]
    xdata = np.linspace(5.6,5.8,400)
    hplot = plotter.add('derivatives test for '+name).gca()
    hplot.plot(xdata,getPeakProfileDerv(parmDict2,xdata,np.zeros_like(xdata),varyList,bakType)[0])
    y0 = getPeakProfile(parmDict2,xdata,np.zeros_like(xdata),varyList,bakType)
    parmDict2[name] += delt
    y1 = getPeakProfile(parmDict2,xdata,np.zeros_like(xdata),varyList,bakType)
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
    y0 = getFCJVoigt3(myDict['pos'],myDict['sig'],myDict['gam'],myDict['shl'],xdata)[0]
    myDict[name] += delt
    y1 = getFCJVoigt3(myDict['pos'],myDict['sig'],myDict['gam'],myDict['shl'],xdata)[0]
    hplot.plot(xdata,(y1-y0)/delt,'r+')

if __name__ == '__main__':
    from . import GSASIItestplot as plot
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

#    GSASIIpath.SetBinaryPath(True,False)
#    print('found',findfullrmc())
