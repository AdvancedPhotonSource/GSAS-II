# -*- coding: utf-8 -*-
"""Element: functions for element types
   Copyright: 2008, Robert B. Von Dreele & Brian H. Toby (Argonne National Laboratory)
"""
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
import numpy as np

def GetFormFactorCoeff(El):
    """Read X-ray form factor coefficients from `atomdata.asc` file

    :param El: element 1-2 character symbol case irrevelant
    :return: `FormFactors`: list of form factor dictionaries
    
    Each X-ray form factor dictionary is:
    
    * `Symbol`: 4 character element symbol with valence (e.g. 'NI+2')
    * `Z`: atomic number
    * `fa`: 4 A coefficients
    * `fb`: 4 B coefficients
    * `fc`: C coefficient
    
    """
    ElS = El.upper()
    ElS = ElS.rjust(2)
    filename = os.path.join(os.path.split(__file__)[0],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        print "**** ERROR - File atmdata.dat not found in directory %s" % sys.path[0]
        sys.exit()
    S = '1'
    FormFactors = []
    while S:
        S = FFdata.readline()
        if S[3:5] == ElS:
            if S[5:6] != '_' and S[8] not in ['N','M']:
                Z=int(S[:2])
                Symbol = S[3:7].strip()
                S = S[12:]
                fa = (float(S[:7]),float(S[14:21]),float(S[28:35]),float(S[42:49]))
                fb = (float(S[7:14]),float(S[21:28]),float(S[35:42]),float(S[49:56]))
                FormFac = {'Symbol':Symbol,'Z':Z,'fa':fa,'fb':fb,'fc':float(S[56:63])}
                FormFactors.append(FormFac)               
    FFdata.close()
    return FormFactors
    
def GetFFC5(ElSym):
    '''Get 5 term form factor and Compton scattering data
    @param ElSym: str(1-2 character element symbol with proper case);
    @return El: dictionary with 5 term form factor & compton coefficients
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
    
def GetAtomInfo(El):
    
    import ElementTable as ET
    Elements = []
    for elem in ET.ElTable:
        Elements.append(elem[0][0])
    if len(El) in [2,4]:
        ElS = El.upper()[:2].rjust(2)
    else:
        ElS = El.upper()[:1].rjust(2)
    filename = os.path.join(os.path.split(__file__)[0],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        print '**** ERROR - File atmdata.dat not found in directory %s' % sys.path[0]
        sys.exit()
    S = '1'
    AtomInfo = {}
    Isotopes = {}
    Mass = []
    while S:
        S = FFdata.readline()
        if S[3:5] == ElS:
            if S[5:6] == '_':
                if not Mass:                                 #picks 1st one; natural abundance or 1st isotope
                    Mass = float(S[10:19])
                if S[6] in [' ','1','2','3','4','5','6','7','8','9']:                        
                    isoName = S[6:9]
                    if isoName == '   ':
                        isoName = 'Nat. Abund.'              #natural abundance
                    if S[76:78] in ['LS','BW']:     #special anomalous scattering length info
                        St = [S[10:19],S[19:25],S[25:31],S[31:38],S[38:44],S[44:50],
                            S[50:56],S[56:62],S[62:68],S[68:74],]
                        Vals = []
                        for item in St:
                            if item.strip():
                                Vals.append(float(item.strip()))
                        Isotopes[isoName.rstrip()] = Vals                        
                    else:
                        Isotopes[isoName.rstrip()] = [float(S[10:19]),float(S[19:25])]
                elif S[5:9] == '_SIZ':
                    Z=int(S[:2])
                    Symbol = S[3:5].strip().lower().capitalize()
                    Drad = float(S[12:22])
                    Arad = float(S[22:32])
                    Vdrad = float(S[32:38])
                    Color = ET.ElTable[Elements.index(Symbol)][6]
    FFdata.close()
    AtomInfo={'Symbol':Symbol,'Isotopes':Isotopes,'Mass':Mass,'Z':Z,'Drad':Drad,'Arad':Arad,'Vdrad':Vdrad,'Color':Color}    
    return AtomInfo
      
def GetXsectionCoeff(El):
    """Read atom orbital scattering cross sections for fprime calculations via Cromer-Lieberman algorithm
    @param El: 2 character element symbol
    @return: Orbs: list of orbitals each a dictionary with detailed orbital information used by FPcalc
    each dictionary is:
    'OrbName': Orbital name read from file
    'IfBe' 0/2 depending on orbital
    'BindEn': binding energy
    'BB': BindEn/0.02721
    'XSectIP': 5 cross section inflection points
    'ElEterm': energy correction term
    'SEdge': absorption edge for orbital
    'Nval': 10/11 depending on IfBe
    'LEner': 10/11 values of log(energy)
    'LXSect': 10/11 values of log(cross section)
    """
    AU = 2.80022e+7
    C1 = 0.02721
    ElS = El.upper()
    ElS = ElS.ljust(2)
    filename = os.path.join(os.path.split(__file__)[0],'Xsect.dat')
    try:
        xsec = open(filename,'Ur')
    except:
        print '**** ERROR - File Xsect.dat not found in directory %s' % sys.path[0]
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
    @param El: 2 character element symbol
    @return: MagFormFactors: list of all magnetic form factors dictionaries for element El.
    each dictionary contains:
    'Symbol':Symbol
    'Z':Z
    'mfa': 4 MA coefficients
    'nfa': 4 NA coefficients
    'mfb': 4 MB coefficients
    'nfb': 4 NB coefficients
    'mfc': MC coefficient
    'nfc': NC coefficient
    """
    ElS = El.upper()
    ElS = ElS.rjust(2)
    filename = os.path.join(os.path.split(__file__)[0],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        print '**** ERROR - File atmdata.dat not found in directory %s' % sys.path[0]
        sys.exit()
    S = '1'
    MagFormFactors = []
    while S:
        S = FFdata.readline()
        if S[3:5] == ElS:
            if S[8:9] == 'M':
                SN = FFdata.readline()               #'N' is assumed to follow 'M' in Atomdata.asc
                Z=int(S[:2])
                Symbol = S[3:7]
                S = S[12:]
                SN = SN[12:]
                mfa = (float(S[:7]),float(S[14:21]),float(S[28:35]),float(S[42:49]))
                mfb = (float(S[7:14]),float(S[21:28]),float(S[35:42]),float(S[49:56]))
                nfa = (float(SN[:7]),float(SN[14:21]),float(SN[28:35]),float(SN[42:49]))
                nfb = (float(SN[7:14]),float(SN[21:28]),float(SN[35:42]),float(SN[49:56]))
                FormFac = {'Symbol':Symbol,'Z':Z,'mfa':mfa,'nfa':nfa,'mfb':mfb,'nfb':nfb,
                    'mfc':float(S[56:63]),'nfc':float(SN[56:63])}
                MagFormFactors.append(FormFac)
    FFdata.close()
    return MagFormFactors

def ScatFac(El, SQ):
    """compute value of form factor
    @param El: element dictionary defined in GetFormFactorCoeff 
    @param SQ: (sin-theta/lambda)**2
    @return: real part of form factor
    """
    fa = np.array(El['fa'])
    fb = np.array(El['fb'])
    t = -fb[:,np.newaxis]*SQ
    return np.sum(fa[:,np.newaxis]*np.exp(t)[:],axis=0)+El['fc']
        
def BlenRes(Elist,BLtables,wave):
    FP = np.zeros(len(Elist))
    FPP = np.zeros(len(Elist))
    Emev = 81.80703/wave**2
    for i,El in enumerate(Elist):
        BL = BLtables[El]
        if len(BL) >= 6:
            Emev = 81.80703/wave**2
            G2 = BL[5]**2
            T = [Emev-BL[4],0,0]
            D = [T**2+G2,0,0]
            fp = T/D
            fpp = 1.0/D
            if len(BL) == 8:
                T = Emev-BL[7]
                D = T**2+G2
                fp += BL[6]*T/D
                fpp += BL[6]/D
            if len(BL) == 10:
                T = Emev-BL[9]
                D = T**2+G2
                fp += BL[8]*T/D
                fpp += BL[8]/D
            FP[i] = (BL[2]*fp)
            FPP[i] = (-BL[3]*fpp)
        else:
            FP[i] = 0.0
            FPP[i] = 0.0
    return FP,FPP
    
#def BlenRes(BLdata,wave):
#    FP = []
#    FPP = []
#    Emev = 81.80703/wave**2
#    for BL in BLdata:
#        if len(BL) >= 6:
#            Emev = 81.80703/wave**2
#            G2 = BL[5]**2
#            T = [Emev-BL[4],0,0]
#            D = [T**2+G2,0,0]
#            fp = T/D
#            fpp = 1.0/D
#            if len(BL) == 8:
#                T = Emev-BL[7]
#                D = T**2+G2
#                fp += BL[6]*T/D
#                fpp += BL[6]/D
#            if len(BL) == 10:
#                T = Emev-BL[9]
#                D = T**2+G2
#                fp += BL[8]*T/D
#                fpp += BL[8]/D
#            FP.append(BL[2]*fp)
#            FPP.append(-BL[3]*fpp)
#        else:
#            FP.append(0.0)
#            FPP.append(0.0)
#    return np.array(FP),np.array(FPP)
    
def ComptonFac(El,SQ):
    """compute Compton scattering factor
    @param El: element dictionary 
    @param SQ: (sin-theta/lambda)**2
    @return: compton scattering factor
    """    
    ca = np.array(El['cmpa'])
    cb = np.array(El['cmpb'])
    t = -cb[:,np.newaxis]*SQ       
    return El['cmpz']-np.sum(ca[:,np.newaxis]*np.exp(t),axis=0) 
            
def FPcalc(Orbs, KEv):
    """Compute real & imaginary resonant X-ray scattering factors
    @param Orbs: list of orbital dictionaries as defined in GetXsectionCoeff
    @param KEv: x-ray energy in keV
    @return: C: (f',f",mu): real, imaginary parts of resonant scattering & atomic absorption coeff.
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
    

