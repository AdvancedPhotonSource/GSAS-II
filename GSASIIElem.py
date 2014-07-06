# -*- coding: utf-8 -*-
"""
*GSASIIElem: functions for element types*
-----------------------------------------

"""
# Copyright: 2008, Robert B. Von Dreele & Brian H. Toby (Argonne National Laboratory)
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################

import math
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import numpy as np
import atmdata

getElSym = lambda sym: sym.split('+')[0].split('-')[0]
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
        
def GetFFC5(ElSym):
    '''Get 5 term form factor and Compton scattering data

    :param ElSym: str(1-2 character element symbol with proper case);
    :return El: dictionary with 5 term form factor & compton coefficients
    '''
    import FormFactors as FF
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
    
def CheckElement(El):
    '''Check if element El is in the periodic table

    :param str El: One or two letter element symbol, capitaliztion ignored
    :returns: True if the element is found

    '''
    import ElementTable as ET
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
    
def GetAtomInfo(El):
    'reads element information from atmdata.py'
    import ElementTable as ET
    Elements = [elem[0][0] for elem in ET.ElTable]
    AtomInfo = {}
    ElS = getElSym(El)
    AtomInfo.update(dict(zip(['Drad','Arad','Vdrad','Hbrad'],atmdata.AtmSize[ElS])))
    AtomInfo['Symbol'] = El
    AtomInfo['Color'] = ET.ElTable[Elements.index(ElS)][6]
    AtomInfo['Z'] = atmdata.XrayFF[El]['Z']
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
    return AtomInfo
        
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
    filename = os.path.join(os.path.split(__file__)[0],'Xsect.dat')
    try:
        xsec = open(filename,'Ur')
    except:
        print '**** ERROR - File Xsect.dat not found in directory %s' % os.path.split(filename)[0]
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
    """Read magnetic form factor data from atomdata.asc file

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
    MagFormFactors = []
    mags = [ky for ky in atmdata.MagFF.keys() if El == getElSym(ky)]
    for mag in mags:
        magData = {}
        data = atmdata.MagFF[mag]
        magData['Symbol'] = mag
        magData['Z'] = atmdata.XrayFF[getElSym(mags)]['Z']
        magData['mfa'] = [data['M'][i] for i in [0,2,4,6]]
        magdata['mfb'] = [data['M'][i] for i in [1,3,5,7]]
        magdata['mfc'] = data['M'][8]
        magData['nfa'] = [data['N'][i] for i in [0,2,4,6]]
        magdata['nfb'] = [data['N'][i] for i in [1,3,5,7]]
        magdata['nfc'] = data['N'][8]
        magdata['g-fac'] = data['N'][9]
        MagFormFactors.append(magdata)
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
    return np.sum(fa[:,np.newaxis]*np.exp(t)[:],axis=0)+El['fc']
        
def BlenResCW(Els,BLtables,wave):
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
            D0 = T0**2+gam**2
            D1 = T1**2+gam**2
            D2 = T2**2+gam**2
            FP[i] = Re*(T0/D0+A*T1/D1+B*T2/D2)
            FPP[i] = Im*(1/D0+A/D1+B/D2)
        else:
            FPP[i] = BL['SL'][1]    #for Li, B, etc.
    return FP,FPP
    
def BlenResTOF(El,BLtables,wave):
    FP = np.zeros(len(wave))
    FPP = np.zeros(len(wave))
    BL = BLtables[El][1]
    print BL
    if 'BW-LS' in BL:
        Re,Im,E0,gam,A,E1,B,E2 = BL['BW-LS'][1:]
        Emev = 81.80703/wave**2
        T0 = Emev-E0
        T1 = Emev-E1
        T2 = Emev-E2
        D0 = T0**2+gam**2
        D1 = T1**2+gam**2
        D2 = T2**2+gam**2
        FP = Re*(T0/D0+A*T1/D1+B*T2/D2)
        FPP = Im*(1/D0+A/D1+B/D2)
    else:
        FPP = np.ones(len(wave))*BL['SL'][1]    #for Li, B, etc.
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
    

