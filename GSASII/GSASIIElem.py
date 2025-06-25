# -*- coding: utf-8 -*-
# Copyright: 2008, Robert B. Von Dreele & Brian H. Toby (Argonne National Laboratory)
"""
Routines used to define element settings follow.
"""

import math
import sys
import os.path
from . import GSASIIpath
import copy
import numpy as np
from . import atmdata
from . import GSASIImath as G2mth
from . import ElementTable as ET

nxs = np.newaxis
Bohr = 0.529177

# Wavelength data
#These from Intl. Tables C, Table 4.2.2.1, p. 177-179
waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
        'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'GaKa':[1.34003,1.34394],
        'MoKa':[0.70926,0.713543],'AgKa':[0.559363,0.563775],'InKa':[0.512094,0.516525]}
# meanwaves computed as (2*Ka1+Ka2)/3
meanwaves = {'CuKa':1.54178,'TiKa':2.74963,'CrKa':2.29092,'FeKa':1.93728,
        'CoKa':1.79021,'MoKa':0.71069,'AgKa':0.56083,'GaKa':1.34134,'Inka':0.51357}

getElSym = lambda sym: sym.split('+')[0].split('-')[0].capitalize()
def GetFormFactorCoeff(El):
    """Read X-ray form factor coefficients from `atomdata.py` file

    :param str El: element 1-2 character symbol, case irrevelant
    :return: `FormFactors`: list of form factor dictionaries

    Each X-ray form factor dictionary is:

    * `Symbol`: 4 character element symbol with valence (e.g. 'NI+2')
    * `Z`: atomic number
    * `fa`: 4 A coefficients
    * `fb`: 4 B coefficients
    * `fc`: C coefficient

    """

    Els = El.capitalize().strip()
    valences = [ky for ky in atmdata.XrayFF.keys() if Els == getElSym(ky)]
    FormFactors = [atmdata.XrayFF[val] for val in valences]
    for Sy,FF in zip(valences,FormFactors):
        FF.update({'Symbol':Sy.upper()})
    return FormFactors

def GetEFormFactorCoeff(El):
    """Read electron form factor coefficients from `atomdata.py` file

    :param str El: element 1-2 character symbol, case irrevelant
    :return: `FormFactors`: list of form factor dictionaries

    Each electrn form factor dictionary is:

    * `Symbol`: 4 character element symbol (no valence)
    * `Z`: atomic number
    * `fa`: 5 A coefficients
    * `fb`: 5 B coefficients

    """

    Els = El.capitalize().strip()
    valences = [ky for ky in atmdata.ElecFF.keys() if Els == getElSym(ky)] #will only be one
    FormFactors = [atmdata.ElecFF[val] for val in valences]
    for Sy,FF in zip(valences,FormFactors):
        FF.update({'Symbol':Sy.upper()})
    return FormFactors

def GetFFtable(atomTypes):
    ''' returns a dictionary of form factor data for atom types found in atomTypes

    :param list atomTypes: list of atom types
    :return: FFtable, dictionary of form factor data; key is atom type

    '''
    FFtable = {}
    for El in atomTypes:
        FFs = GetFormFactorCoeff(getElSym(El))
        for item in FFs:
            if item['Symbol'] == El.upper():
                FFtable[El] = item
    return FFtable

def GetEFFtable(atomTypes):
    ''' returns a dictionary of electron form factor data for atom types found in atomTypes
    might not be needed?

    :param list atomTypes: list of atom types
    :return: FFtable, dictionary of form factor data; key is atom type

    '''
    FFtable = {}
    for El in atomTypes:
        FFs = GetEFormFactorCoeff(getElSym(El))
        for item in FFs:
            if item['Symbol'] == El.upper():
                FFtable[El] = item
    return FFtable

def GetORBtable(atomTypes):
    ''' returns a dictionary of orbital form factor data for atom types found in atomTypes

    :param list atomTypes: list of atom types
    :return: ORBtable, dictionary of orbital form factor data; key is atom type

    '''
    ORBtable = {}
    for El in atomTypes:
        ORBtable[El] = copy.deepcopy(atmdata.OrbFF[El])
    return ORBtable

def GetMFtable(atomTypes,Landeg):
    ''' returns a dictionary of magnetic form factor data for atom types found in atomTypes

    :param list atomTypes: list of atom types
    :param list Landeg: Lande g factors for atomTypes
    :return: FFtable, dictionary of form factor data; key is atom type

    '''
    MFtable = {}
    for El,gfac in zip(atomTypes,Landeg):
        MFs = GetMagFormFacCoeff(getElSym(El))
        for item in MFs:
            if item['Symbol'] == El.upper():
                item['gfac'] = gfac
                MFtable[El] = item
    return MFtable

def GetBLtable(General):
    ''' returns a dictionary of neutron scattering length data for atom types & isotopes found in General

    :param dict General: dictionary of phase info.; includes AtomTypes & Isotopes
    :return: BLtable, dictionary of scattering length data; key is atom type
    '''
    atomTypes = General['AtomTypes']
    BLtable = {}
    isotope = General['Isotope']
    for El in atomTypes:
        ElS = getElSym(El)
        if 'Nat' in isotope[El]:
            BLtable[El] = [isotope[El],atmdata.AtmBlens[ElS+'_']]
        else:
            BLtable[El] = [isotope[El],atmdata.AtmBlens[ElS+'_'+isotope[El]]]
    return BLtable

def getFFvalues(FFtables,SQ,ifList=False):
    'Needs a doc string'
    if ifList:
        FFvals = []
        for El in FFtables:
            FFvals.append(ScatFac(FFtables[El],SQ)[0])
    else:
        FFvals = {}
        for El in FFtables:
            FFvals[El] = ScatFac(FFtables[El],SQ)[0]
    return FFvals

def getBLvalues(BLtables,ifList=False):
    'Needs a doc string'
    if ifList:
        BLvals = []
        for El in BLtables:
            if 'BW-LS' in El:
                BLvals.append(BLtables[El][1]['BW-LS'][0])
            else:
                BLvals.append(BLtables[El][1]['SL'][0])
    else:
        BLvals = {}
        for El in BLtables:
            if 'BW-LS' in El:
                BLvals[El] = BLtables[El][1]['BW-LS'][0]
            else:
                BLvals[El] = BLtables[El][1]['SL'][0]
    return BLvals

def getMFvalues(MFtables,SQ,ifList=False):
    'Needs a doc string'
    if ifList:
        MFvals = []
        for El in MFtables:
            MFvals.append(MagScatFac(MFtables[El],SQ)[0])
    else:
        MFvals = {}
        for El in MFtables:
            MFvals[El] = MagScatFac(MFtables[El],SQ)[0]
    return MFvals

def GetFFC5(ElSym):
    '''Get 5 term form factor and Compton scattering data

    :param ElSym: str(1-2 character element symbol with proper case);
    :return El: dictionary with 5 term form factor & compton coefficients
    '''
    from . import FormFactors as FF
    El = {}
    FF5 = FF.FFac5term[ElSym]
    El['fa'] = FF5[:5]
    El['fc'] = FF5[5]
    El['fb'] = FF5[6:]
    Cmp5 = FF.Compton[ElSym]
    El['cmpz'] = Cmp5[0]
    El['cmpa'] = Cmp5[1:6]
    El['cmpb'] = Cmp5[6:]
    return El

def GetBVS(Pair,atSeq,Valences):
    Els = Pair.strip().split('-')
    iAt = atSeq.index(Els[0])
    iVal = Valences[iAt][0]
    if Els[1] in ['O','F','Cl']:
        iEls = ['O','F','Cl'].index(Els[1])
        if iVal in atmdata.BVScoeff:
            return atmdata.BVScoeff[iVal][iEls]
        else:
            return 0.0
    elif Els[1] in ['Br','I','S','Se','Te','N','P','As','H','D']:
        iEls = ['Br','I','S','Se','Te','N','P','As','H','D'].index(Els[1])
        if Els[0] in atmdata.BVSnotOFCl:
            return atmdata.BVSnotOFCl[Els[0]][iEls]
        else:
            return 0.0
    else:
        return 0.0

def CheckElement(El):
    '''Check if element El is in the periodic table

    :param str El: One or two letter element symbol, capitaliztion ignored
    :returns: True if the element is found

    '''
    Elements = []
    for elem in ET.ElTable:
        Elements.append(elem[0][0])
    if El.capitalize() in Elements:
        return True
    else:
        return False

def FixValence(El):
    'Returns the element symbol, even when a valence is present'
    if '+' in El[-1]: #converts An+/- to A+/-n
        num = El[-2]
        El = El.split(num)[0]+'+'+num
    if '+0' in El:
        El = El.split('+0')[0]
    if '-' in El[-1]:
        num = El[-2]
        El = El.split(num)[0]+'-'+num
    if '-0' in El:
        El = El.split('-0')[0]
    return El

def SetAtomColor(El,RGB):
    'Overrides the default color in the atoms table; not saved'
    Elem = ET.ElTable
    Elements = [elem[0][0] for elem in Elem]
    if 'Q' in El: El = 'Q'      #patch - remove Qa, etc.
    ElS = getElSym(El)
    ET.ElTable[Elements.index(ElS)] = ET.ElTable[Elements.index(ElS)][0:6] + (
        tuple(RGB),)

def GetAtomInfo(El,ifMag=False):
    'reads element information from atmdata.py'
    Elem = ET.ElTable
    if ifMag:
        Elem = ET.MagElTable
    Elements = [elem[0][0] for elem in Elem]
    AtomInfo = {}
    if 'Q' in El: El = 'Q'      #patch - remove Qa, etc.
    ElS = getElSym(El)
    if El not in atmdata.XrayFF and El not in atmdata.MagFF:
        if ElS not in atmdata.XrayFF:
            if ElS.endswith('0') and ElS[:-1] in atmdata.XrayFF:
                ElS = ElS[:-1]
            else:
                ElS = 'H'
        print('Atom type '+El+' not found, using '+ElS)
        El = ElS
    AtomInfo.update(dict(zip(['Drad','Arad','Vdrad','Hbrad'],atmdata.AtmSize[ElS])))
    AtomInfo['Symbol'] = El
    AtomInfo['Color'] = ET.ElTable[Elements.index(ElS)][6]
    AtomInfo['Z'] = atmdata.XrayFF[ElS]['Z']
    isotopes = [ky for ky in atmdata.AtmBlens.keys() if ElS == ky.split('_')[0]]
    isotopes.sort()
    AtomInfo['Mass'] = atmdata.AtmBlens[isotopes[0]]['Mass']    #default to nat. abund.
    AtomInfo['Isotopes'] = {}
    for isotope in isotopes:
        data = atmdata.AtmBlens[isotope]
        if isotope == ElS+'_':
            AtomInfo['Isotopes']['Nat. Abund.'] = data
        else:
            AtomInfo['Isotopes'][isotope.split('_')[1]] = data
    AtomInfo['Lande g'] = 2.0
    return AtomInfo

def GetElInfo(El,inst):
    ElemSym = El.strip().capitalize()
    if 'X' in inst['Type'][0]:
        keV = 12.397639/G2mth.getWave(inst)
        FpMu = FPcalc(GetXsectionCoeff(ElemSym), keV)
        ElData = GetFormFactorCoeff(ElemSym)[0]
        ElData['FormulaNo'] = 0.0
        ElData.update(GetAtomInfo(ElemSym))
        ElData.update(dict(zip(['fp','fpp','mu'],FpMu)))
        ElData.update(GetFFC5(El))
    else: #'N'eutron
        ElData = {}
        ElData.update(GetAtomInfo(ElemSym))
        ElData['FormulaNo'] = 0.0
        ElData.update({'mu':0.0,'fp':0.0,'fpp':0.0})
    return ElData

def GetXsectionCoeff(El):
    """Read atom orbital scattering cross sections for fprime calculations via Cromer-Lieberman algorithm

    :param El: 2 character element symbol
    :return: Orbs: list of orbitals each a dictionary with detailed orbital information used by FPcalc

    each dictionary is:

    * 'OrbName': Orbital name read from file
    * 'IfBe' 0/2 depending on orbital
    * 'BindEn': binding energy
    * 'BB': BindEn/0.02721
    * 'XSectIP': 5 cross section inflection points
    * 'ElEterm': energy correction term
    * 'SEdge': absorption edge for orbital
    * 'Nval': 10/11 depending on IfBe
    * 'LEner': 10/11 values of log(energy)
    * 'LXSect': 10/11 values of log(cross section)

    """
    AU = 2.80022e+7
    C1 = 0.02721
    ElS = El.upper()
    ElS = ElS.ljust(2)
    filename = os.path.join(GSASIIpath.path2GSAS2,'inputs','Xsect.dat')
    if not os.path.exists(filename):  # patch 3/2024 for svn dir organization
        filename = os.path.join(GSASIIpath.path2GSAS2,'Xsect.dat')
    try:
        xsec = open(filename,'r')
    except:
        print (f'**** ERROR - File Xsect.dat not found in directory {os.path.dirname(filename)}')
        sys.exit()
    S = '1'
    Orbs = []
    while S:
        S = xsec.readline()
        if S[:2] == ElS:
            S = S[:-1]+xsec.readline()[:-1]+xsec.readline()
            OrbName = S[9:14]
            S = S[14:]
            IfBe = int(S[0])
            S = S[1:]
            val = S.split()
            BindEn = float(val[0])
            BB = BindEn/C1
            Orb = {'OrbName':OrbName,'IfBe':IfBe,'BindEn':BindEn,'BB':BB}
            Energy = []
            XSect = []
            for i in range(11):
                Energy.append(float(val[2*i+1]))
                XSect.append(float(val[2*i+2]))
            XSecIP = []
            for i in range(5): XSecIP.append(XSect[i+5]/AU)
            Orb['XSecIP'] = XSecIP
            if IfBe == 0:
                Orb['SEdge'] = XSect[10]/AU
                Nval = 11
            else:
                Orb['ElEterm'] = XSect[10]
                del Energy[10]
                del XSect[10]
                Nval = 10
                Orb['SEdge'] = 0.0
            Orb['Nval'] = Nval
            D = dict(zip(Energy,XSect))
            Energy.sort()
            X = []
            for key in Energy:
                X.append(D[key])
            XSect = X
            LEner = []
            LXSect = []
            for i in range(Nval):
                LEner.append(math.log(Energy[i]))
                if XSect[i] > 0.0:
                    LXSect.append(math.log(XSect[i]))
                else:
                    LXSect.append(0.0)
            Orb['LEner'] = LEner
            Orb['LXSect'] = LXSect
            Orbs.append(Orb)
    xsec.close()
    return Orbs

def GetMagFormFacCoeff(El):
    """Read magnetic form factor data from atmdata.py

    :param El: 2 character element symbol
    :return: MagFormFactors: list of all magnetic form factors dictionaries for element El.

    each dictionary contains:

    * 'Symbol':Symbol
    * 'Z':Z
    * 'mfa': 4 MA coefficients
    * 'nfa': 4 NA coefficients
    * 'mfb': 4 MB coefficients
    * 'nfb': 4 NB coefficients
    * 'mfc': MC coefficient
    * 'nfc': NC coefficient

    """
    Els = El.capitalize().strip()
    MagFormFactors = []
    mags = [ky for ky in atmdata.MagFF.keys() if Els == getElSym(ky)]
    for mag in mags:
        magData = {}
        data = atmdata.MagFF[mag]
        magData['Symbol'] = mag.upper()
        magData['Z'] = atmdata.XrayFF[getElSym(mag)]['Z']
        magData['mfa'] = [data['M'][i] for i in [0,2,4,6]]
        magData['mfb'] = [data['M'][i] for i in [1,3,5,7]]
        magData['mfc'] = data['M'][8]
        magData['nfa'] = [data['N'][i] for i in [0,2,4,6]]
        magData['nfb'] = [data['N'][i] for i in [1,3,5,7]]
        magData['nfc'] = data['N'][8]
        MagFormFactors.append(magData)
    return MagFormFactors

def ScatFac(El, SQ):
    """compute value of form factor

    :param El: element dictionary defined in GetFormFactorCoeff
    :param SQ: (sin-theta/lambda)**2
    :return: real part of form factor
    """
    fa = np.array(El['fa'])
    fb = np.array(El['fb'])
    t = -fb[:,np.newaxis]*SQ
    return np.sum(fa[:,np.newaxis]*np.exp(t)[:],axis=0)+El.get('fc',0.0)

def ScatFacDer(El, SQ):
    """compute derivative of form factor wrt SQ

    :param El: element dictionary defined in GetFormFactorCoeff
    :param SQ: (sin-theta/lambda)**2
    :return: real part of form factor
    """
    fa = np.array(El['fa'])
    fb = np.array(El['fb'])
    t = -fb[:,np.newaxis]*SQ
    return -np.sum(fa[:,np.newaxis]*fb[:,np.newaxis]*np.exp(t)[:],axis=0)

def MagScatFac(El, SQ):
    """compute value of form factor

    :param El: element dictionary defined in GetFormFactorCoeff
    :param SQ: (sin-theta/lambda)**2
    :param gfac: Lande g factor (normally = 2.0)
    :return: real part of form factor
    """
    mfa = np.array(El['mfa'])
    mfb = np.array(El['mfb'])
    nfa = np.array(El['nfa'])
    nfb = np.array(El['nfb'])
    mt = -mfb[:,np.newaxis]*SQ
    nt = -nfb[:,np.newaxis]*SQ
    MMF = np.sum(mfa[:,np.newaxis]*np.exp(mt)[:],axis=0)+El['mfc']
    MMF0 = np.sum(mfa)+El['mfc']
    NMF = np.sum(nfa[:,np.newaxis]*np.exp(nt)[:],axis=0)+El['nfc']
    NMF0 = np.sum(nfa)+El['nfc']
    MF0 = MMF0+(2.0/El['gfac']-1.0)*NMF0
    return (MMF+(2.0/El['gfac']-1.0)*NMF)/MF0
  
def ClosedFormFF(Z,SQ,k,N):
    """Closed form expressions for FT Slater fxns. IT B Table 1.2.7.4
    (not used at present - doesn't make sense yet)

    :param Z: element zeta factor
    :param SQ: (sin-theta/lambda)**2
    :param k: int principal Bessel fxn order as in <jk>
    :param N: int power

    return: form factor
    """
    Z2 = Z**2
    K2 = 16.0*SQ*np.pi**2
    K2pZ2 = K2+Z2
    K = np.sqrt(K2)
    if k == 0:
        if N == 1:
            return 1.0/K2pZ2
        elif N == 2:
            return 2.0*Z/K2pZ2**2
        elif N == 3:
            return 2.0*(3.0*Z2-K2)/K2pZ2**3
        elif N == 4:
            return 24.0*Z*(Z2-K2)/K2pZ2**4
        elif N == 5:
            return 24.0*(5.0*Z2-10.0*K2*Z2+K2**2)/K2pZ2**5
        elif N == 6:
            return 240.0*(K2-3.0*Z2)*(3.0*K2-Z2)/K2pZ2**6
        elif N == 7:
            return 720.0*(7.0*Z2**3-35.0*K2*Z2**2+21.0*Z2*K2**2-K2**3)/K2pZ2**7
        elif N == 8:
            return 40320.0*(Z*Z2**3-7.0*K2*Z*Z2**2+7.0*K2**2*Z*Z2-Z*K2**3)/K2pZ2**8
    elif k == 1:
        if N == 2:
            return 2.0*K/K2pZ2**2
        elif N == 3:
            return 8.0*K*Z/K2pZ2**3
        elif N == 4:
            return 8.0*K*(5.0*Z2-K2)/K2pZ2**4
        elif N == 5:
            return 48.0*K*Z*(5.0*Z2-3.0*K2)/K2pZ2**5
        elif N == 6:
            return 48.0*K*(35.0*Z2**2-42.0*K2*Z2+3.0*K2**2)/K2pZ2**6
        elif N == 7:
            return 1920.0*K*Z*(7.0*Z2**2-14.0*K2*Z2+3.0*K2**2)/K2pZ2**7
        elif N == 8:
            return 5760.0*K*(21.0*Z2**3-63.0*K2*Z2**2+27.0*K2**2*Z2-K2**3)/K2pZ2**8
    elif k == 2:
        if N == 3:
            return 8.0*K2/K2pZ2**3
        elif N == 4:
            return 48.0*K2*Z/K2pZ2**4
        elif N == 5:
            return 48.0*K2*(7.0*Z2-K2)/K2pZ2**5
        elif N == 6:
            return 384.0*K2*Z*(7.0*Z2-3.0*K2)/K2pZ2**6
        elif N == 7:
            return 1152.0*K2*(21.0*Z2**2-18.0*K2*Z2+K2**2)/K2pZ2**7
        elif N == 8:
            return 11520.0*K2*Z*(21.0*Z2**2-30.0*K2*Z2+5.0*K2**2)/K2pZ2**8
    elif k == 3:
        if N == 4:
            return 48.0*K**3/K2pZ2**4
        elif N == 5:
            return 384.0*K**3*Z/K2pZ2**5
        elif N == 6:
            return 384.0*K**3*(9.0*Z2-K2)/K2pZ2**6
        elif N == 7:
            return 11520.0*K**3*Z*(3.0*Z2-K2)/K2pZ2**7
        elif N == 8:
            return 11520.0*K**3*(33.0*Z2**2-22.0*K2*Z2+K2**2)/K2pZ2**8
    elif k == 4:
        if N == 5:
            return 384.0*K2**2/K2pZ2**5
        elif N == 6:
            return 3840.0*K2**2*Z/K2pZ2**6
        elif N == 7:
            return 3840.0*K2**2*(11.0*Z2-K2)/K2pZ2**7
        elif N == 8:
            return 46080.0*K**5*(13.0*Z2-K2)/K2pZ2**8
    elif k == 5:
        if N == 6:
            return 3840.0*K**5/K2pZ2**6
        elif N == 7:
            return 46080.0*Z*K**5/K2pZ2**7
        elif N == 8:
            return 46080.0*K**5*(13.0*Z2-K2)/K2pZ2**8
    elif k == 6:
        if N == 7:
            return 46080.0*K**6/K2pZ2**7
        elif N == 8:
            return 645120.0*Z*K2**3/K2pZ2**8
    elif k == 7:
        if N == 8:
            return 645120.0*K**7/K2pZ2**8

def BlenResCW(Els,BLtables,wave):
    ''' Computes resonant scattering lengths - single wavelength version (CW)
    returns bo+b' and b"'
    '''
    FP = np.zeros(len(Els))
    FPP = np.zeros(len(Els))
    for i,El in enumerate(Els):
        BL = BLtables[El][1]
        if 'BW-LS' in BL:
            Re,Im,E0,gam,A,E1,B,E2 = BL['BW-LS'][1:]
            Emev = 81.80703/wave**2
            T0 = Emev-E0
            T1 = Emev-E1
            T2 = Emev-E2
            D0 = T0**2+gam**2/4.
            D1 = T1**2+gam**2/4.
            D2 = T2**2+gam**2/4.
            FP[i] = Re*(T0/D0+A*T1/D1+B*T2/D2)+BL['BW-LS'][0]
            FPP[i] = Im*(1/D0+A/D1+B/D2)
        else:
            FPP[i] = BL['SL'][1]    #for Li, B, etc.
    return FP,FPP

def BlenResTOF(Els,BLtables,wave):
    ''' Computes resonant scattering lengths - multiple wavelength version (TOF)
    returns bo+b' and b"'
    '''
    FP = np.zeros((len(Els),len(wave)))
    FPP = np.zeros((len(Els),len(wave)))
    BL = [BLtables[el][1] for el in Els]
    for i,El in enumerate(Els):
        if 'BW-LS' in BL[i]:
            Re,Im,E0,gam,A,E1,B,E2 = BL[i]['BW-LS'][1:]
            Emev = 81.80703/wave**2
            T0 = Emev-E0
            T1 = Emev-E1
            T2 = Emev-E2
            D0 = T0**2+gam**2/4.
            D1 = T1**2+gam**2/4.
            D2 = T2**2+gam**2/4.
            FP[i] = Re*(T0/D0+A*T1/D1+B*T2/D2)+BL[i]['BW-LS'][0]
            FPP[i] = Im*(1/D0+A/D1+B/D2)
        else:
            FPP[i] = np.ones(len(wave))*BL[i]['SL'][1]    #for Li, B, etc.
    return FP,FPP

def ComptonFac(El,SQ):
    """compute Compton scattering factor

    :param El: element dictionary
    :param SQ: (sin-theta/lambda)**2
    :return: compton scattering factor
    """
    ca = np.array(El['cmpa'])
    cb = np.array(El['cmpb'])
    t = -cb[:,np.newaxis]*SQ
    return El['cmpz']-np.sum(ca[:,np.newaxis]*np.exp(t),axis=0)

def FPcalc(Orbs, KEv):
    """Compute real & imaginary resonant X-ray scattering factors

    :param Orbs: list of orbital dictionaries as defined in GetXsectionCoeff
    :param KEv: x-ray energy in keV
    :return: C: (f',f",mu): real, imaginary parts of resonant scattering & atomic absorption coeff.
    """
    def Aitken(Orb, LKev):
        Nval = Orb['Nval']
        j = Nval-1
        LEner = Orb['LEner']
        for i in range(Nval):
            if LEner[i] <= LKev: j = i
        if j > Nval-3: j= Nval-3
        T = [0,0,0,0,0,0]
        LXSect = Orb['LXSect']
        for i in range(3):
           T[i] = LXSect[i+j]
           T[i+3] = LEner[i+j]-LKev
        T[1] = (T[0]*T[4]-T[1]*T[3])/(LEner[j+1]-LEner[j])
        T[2] = (T[0]*T[5]-T[2]*T[3])/(LEner[j+2]-LEner[j])
        T[2] = (T[1]*T[5]-T[2]*T[4])/(LEner[j+2]-LEner[j+1])
        C = T[2]
        return C

    def DGauss(Orb,CX,RX,ISig):
        ALG = (0.11846344252810,0.23931433524968,0.284444444444,
        0.23931433524968,0.11846344252810)
        XLG = (0.04691007703067,0.23076534494716,0.5,
        0.76923465505284,0.95308992296933)

        D = 0.0
        B2 = Orb['BB']**2
        R2 = RX**2
        XSecIP = Orb['XSecIP']
        for i in range(5):
            X = XLG[i]
            X2 = X**2
            XS = XSecIP[i]
            if ISig == 0:
                S = BB*(XS*(B2/X2)-CX*R2)/(R2*X2-B2)
            elif ISig == 1:
                S = 0.5*BB*B2*XS/(math.sqrt(X)*(R2*X2-X*B2))
            elif ISig == 2:
                T = X*X2*R2-B2/X
                S = 2.0*BB*(XS*B2/(T*X2**2)-(CX*R2/T))
            else:
                S = BB*B2*(XS-Orb['SEdge']*X2)/(R2*X2**2-X2*B2)
            A = ALG[i]
            D += A*S
        return D

    AU = 2.80022e+7
    C1 = 0.02721
    C = 137.0367
    FP = 0.0
    FPP = 0.0
    Mu = 0.0
    LKev = math.log(KEv)
    RX = KEv/C1
    if Orbs:
        for Orb in Orbs:
            CX = 0.0
            BB = Orb['BB']
            BindEn = Orb['BindEn']
            if Orb['IfBe'] != 0: ElEterm = Orb['ElEterm']
            if BindEn <= KEv:
                CX = math.exp(Aitken(Orb,LKev))
                Mu += CX
                CX /= AU
            Corr = 0.0
            if Orb['IfBe'] == 0 and BindEn >= KEv:
                CX = 0.0
                FPI = DGauss(Orb,CX,RX,3)
                Corr = 0.5*Orb['SEdge']*BB**2*math.log((RX-BB)/(-RX-BB))/RX
            else:
                FPI = DGauss(Orb,CX,RX,Orb['IfBe'])
                if CX != 0.0: Corr = -0.5*CX*RX*math.log((RX+BB)/(RX-BB))
            FPI = (FPI+Corr)*C/(2.0*math.pi**2)
            FPPI = C*CX*RX/(4.0*math.pi)
            FP += FPI
            FPP += FPPI
        FP -= ElEterm

    return (FP, FPP, Mu)

mapDefault = {'MapType':'','RefList':'','GridStep':0.25,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}

def SetupGeneral(data, dirname):
    '''Initialize the General sections of the Phase tree contents. Should
    be done after changes to the Atoms array.

    Called by routine SetupGeneral (in :func:`GSASIIphsGUI.UpdatePhaseData`),
    :func:`GSASIIphsGUI.makeIsoNewPhase`, :func:`GSASIImiscGUI.saveNewPhase`,
    and in :func:`GSASIIscriptable.SetupGeneral`.
    '''
    generalData = data['General']
    atomData = data['Atoms']
    generalData['AtomTypes'] = []
    generalData['Isotopes'] = {}
    RBModels = data.get('RBModels',{})
# various patches
    if 'Isotope' not in generalData:
        generalData['Isotope'] = {}
    if 'Data plot type' not in generalData:
        generalData['Data plot type'] = 'Mustrain'
    if 'POhkl' not in generalData:
        generalData['POhkl'] = [0,0,1]
    if 'Map' not in generalData:
        generalData['Map'] = mapDefault.copy()
    if 'Flip' not in generalData:
        generalData['Flip'] = {'RefList':'','GridStep':0.25,'Norm element':'None',
            'k-factor':0.1,'k-Max':20.,}
    if 'testHKL' not in generalData['Flip']:
        generalData['Flip']['testHKL'] = [[0,0,2],[2,0,0],[1,1,1],[0,2,0],[1,2,3]]
    if 'doPawley' not in generalData:
        generalData['doPawley'] = False     #ToDo: change to ''
    if 'Pawley dmin' not in generalData:
        generalData['Pawley dmin'] = 1.0
    if 'Pawley dmax' not in generalData:
        generalData['Pawley dmax'] = 100.0
    if 'Pawley neg wt' not in generalData:
        generalData['Pawley neg wt'] = 0.0
    if '3Dproj' not in generalData:
        generalData['3Dproj'] = ''
    if 'doDysnomia' not in generalData:
        generalData['doDysnomia'] = False
    if 'Algolrithm' in generalData.get('MCSA controls',{}) or \
        'MCSA controls' not in generalData:
        generalData['MCSA controls'] = {'Data source':'','Annealing':[0.7,0.1,250],
        'dmin':2.8,'Algorithm':'log','fast parms':[0.8,0.6],'log slope':0.9,
        'Cycles':1,'Results':[],'newDmin':True}
    if 'AtomPtrs' not in generalData:
        generalData['AtomPtrs'] = [3,1,7,9]
        if generalData['Type'] == 'macromolecular':
            generalData['AtomPtrs'] = [6,4,10,12]
        elif generalData['Type'] == 'magnetic':
            generalData['AtomPtrs'] = [3,1,10,12]
    if generalData['Modulated']:
        if 'Super' not in generalData:
            generalData['Super'] = 1
            generalData['SuperVec'] = [[0.,0.,0.],False,1]
            generalData['SSGData'] = {}
        if '4DmapData' not in generalData:
            generalData['4DmapData'] = mapDefault.copy()
            generalData['4DmapData'].update({'MapType':'Fobs'})
        atomData = data['Atoms']
        for atom in atomData:
#                if 'SS1' not in atom:
#                    atom += [[],[],{'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}]
            if isinstance(atom[-1],dict) and 'waveType' in atom[-1]['SS1']:
                waveType = atom[-1]['SS1']['waveType']
                for parm in ['Sfrac','Spos','Sadp','Smag']:
                    if len(atom[-1]['SS1'][parm]):
                        wType = 'Fourier'
                        if parm == 'Sfrac':
                            if 'Crenel' in waveType:
                                wType = 'Crenel'
                        elif parm == 'Spos':
                            if not 'Crenel' in waveType:
                                wType = waveType
                        atom[-1]['SS1'][parm] = [wType,]+list(atom[-1]['SS1'][parm])
                del atom[-1]['SS1']['waveType']
    else:
        generalData['Super'] = 0
    if 'Modulated' not in generalData:
        generalData['Modulated'] = False
    if 'HydIds' not in generalData:
        generalData['HydIds'] = {}
    if generalData['Type'] == 'magnetic':
        if 'SGGray' not in generalData['SGData']:
            generalData['SGData']['SGGray'] = False
    if 'Resolution' in generalData['Map']:
        generalData['Map']['GridStep'] = generalData['Map']['Resolution']/2.
        generalData['Flip']['GridStep'] = generalData['Flip']['Resolution']/2.
        del generalData['Map']['Resolution']
        del generalData['Flip']['Resolution']
    if 'Compare' not in generalData:
        generalData['Compare'] = {'Oatoms':'','Tatoms':'',
                'Tilts':{'Otilts':[],'Ttilts':[]},
            'Bonds':{'Obonds':[],'Tbonds':[]},'Vects':{'Ovec':[],'Tvec':[]},
            'dVects':{'Ovec':[],'Tvec':[]},'Sampling':1.0}
    if 'Sampling' not in generalData['Compare']:
        generalData['Compare']['Sampling'] = 1.0
    generalData['SpnIds'] = generalData.get('SpnIds',{})

# end of patches
    cx,ct,cs,cia = generalData['AtomPtrs']
    generalData['NoAtoms'] = {}
    generalData['BondRadii'] = []
    generalData['AngleRadii'] = []
    generalData['vdWRadii'] = []
    generalData['AtomMass'] = []
    generalData['Color'] = []
    if generalData['Type'] == 'magnetic':
        generalData['MagDmin'] = generalData.get('MagDmin',1.0)
        landeg = generalData.get('Lande g',[])
    generalData['Mydir'] = dirname
    badList = {}
    for iat,atom in enumerate(atomData):
        atom[ct] = atom[ct].lower().capitalize()              #force to standard form
        if generalData['AtomTypes'].count(atom[ct]):
            generalData['NoAtoms'][atom[ct]] += atom[cx+3]*float(atom[cs+1])
        elif atom[ct] != 'UNK':
            Info = GetAtomInfo(atom[ct])
            if not Info:
                if atom[ct] not in badList:
                    badList[atom[ct]] = 0
                badList[atom[ct]] += 1
                atom[ct] = 'UNK'
                continue
            atom[ct] = Info['Symbol'] # N.B. symbol might be changed by GetAtomInfo
            generalData['AtomTypes'].append(atom[ct])
            generalData['Z'] = Info['Z']
            generalData['Isotopes'][atom[ct]] = Info['Isotopes']
            generalData['BondRadii'].append(Info['Drad'])
            generalData['AngleRadii'].append(Info['Arad'])
            generalData['vdWRadii'].append(Info['Vdrad'])
            if atom[ct] in generalData['Isotope']:
                if generalData['Isotope'][atom[ct]] not in generalData['Isotopes'][atom[ct]]:
                    isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]]['Mass'])
            else:
                generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                if 'Nat. Abund.' not in generalData['Isotopes'][atom[ct]]:
                    isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Mass'])
            generalData['NoAtoms'][atom[ct]] = atom[cx+3]*float(atom[cs+1])
            generalData['Color'].append(Info['Color'])
            if generalData['Type'] == 'magnetic':
                if len(landeg) < len(generalData['AtomTypes']):
                    landeg.append(2.0)
        if 'Q' in atom[ct]:
            atom[ct] = 'Q'  #patch - remove 'QA', etc.
            for Srb in RBModels.get('Spin',[]):
                if Srb['Ids'][0] != atom[cia+8]:
                    continue
                nSh = len(Srb['RBId'])
                for iSh in range(nSh):
                    Info = GetAtomInfo(Srb['atType'][iSh])
                    if Info['Symbol'] not in generalData['AtomTypes']:
                        generalData['AtomTypes'].append(Info['Symbol'])
                        generalData['Z'] = Info['Z']
                        generalData['Isotopes'][Info['Symbol']] = Info['Isotopes']
                        generalData['BondRadii'].append(Info['Drad'])
                        generalData['AngleRadii'].append(Info['Arad'])
                        generalData['vdWRadii'].append(Info['Vdrad'])
                        if Info['Symbol'] in generalData['Isotope']:
                            if generalData['Isotope'][Info['Symbol']] not in generalData['Isotopes'][Info['Symbol']]:
                                isotope = list(generalData['Isotopes'][Info['Symbol']].keys())[-1]
                                generalData['Isotope'][Info['Symbol']] = isotope
                            generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][Info['Symbol']]]['Mass'])
                        else:
                            generalData['Isotope'][Info['Symbol']] = 'Nat. Abund.'
                            if 'Nat. Abund.' not in generalData['Isotopes'][Info['Symbol']]:
                                isotope = list(generalData['Isotopes'][Info['Symbol']].keys())[-1]
                                generalData['Isotope'][Info['Symbol']] = isotope
                            generalData['AtomMass'].append(Info['Mass'])
                        generalData['NoAtoms'][Info['Symbol']] = atom[cx+3]*atom[cs+1]*Srb['Natoms'][iSh]
                        generalData['Color'].append(Info['Color'])
                    else:
                        generalData['NoAtoms'][Info['Symbol']] += atom[cx+3]*atom[cs+1]*Srb['Natoms'][iSh]

    if generalData['Type'] == 'magnetic':
        generalData['Lande g'] = landeg[:len(generalData['AtomTypes'])]

    F000X = 0.
    F000N = 0.
    for i,elem in enumerate(generalData['AtomTypes']):
        F000X += generalData['NoAtoms'][elem]*generalData['Z']
        isotope = generalData['Isotope'][elem]
        F000N += generalData['NoAtoms'][elem]*generalData['Isotopes'][elem][isotope]['SL'][0]
    generalData['F000X'] = F000X
    generalData['F000N'] = F000N
    generalData['Mass'] = G2mth.getMass(generalData)

    if badList:
        msg = 'Warning: element symbol(s) not found:'
        for key in badList:
            msg += '\n\t' + key
            if badList[key] > 1:
                msg += ' (' + str(badList[key]) + ' times)'
        #wx.MessageBox(msg,caption='Element symbol error')
        raise ValueError("Phase error:\n" + msg)
