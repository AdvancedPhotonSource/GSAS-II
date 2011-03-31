"""Element: functions for element types
   Copyright: 2008, Robert B. Von Dreele (Argonne National Laboratory)
"""
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################

import wx
import math
import sys
import os.path
import wx.lib.colourselect as wscs
import GSASIIpath
import numpy as np

def GetFormFactorCoeff(El):
    """Read form factor coefficients from atomdata.asc file
    @param El: element 1-2 character symbol case irrevelant
    @return: FormFactors: list of form factor dictionaries
    each dictionary is:
    'Symbol':4 character element symbol with valence (e.g. 'NI+2')
    'Z': atomic number
    'fa': 4 A coefficients
    'fb':4 B coefficients
    'fc': C coefficient 
    """
    ElS = El.upper()
    ElS = ElS.rjust(2)
    filename = os.path.join(sys.path[1],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        wx.MessageBox(message="File atmdata.dat not found in directory %s" % sys.path[0],
            caption="No atmdata.dat file",style=wx.OK | wx.ICON_EXCLAMATION | wx.STAY_ON_TOP)
        sys.exit()
    S = '1'
    FormFactors = []
    while S:
        S = FFdata.readline()
        if S[3:5] == ElS:
            if S[5:6] != '_':
                Z=int(S[:2])
                Symbol = S[3:7].strip()
                S = S[12:]
                fa = (float(S[:7]),float(S[14:21]),float(S[28:35]),float(S[42:49]))
                fb = (float(S[7:14]),float(S[21:28]),float(S[35:42]),float(S[49:56]))
                FormFac = {'Symbol':Symbol,'Z':Z,'fa':fa,'fb':fb,'fc':float(S[56:63])}
                FormFactors.append(FormFac)               
    FFdata.close()
    return FormFactors
    
def GetAtomInfo(El):
    
    import ElementTable as ET
    Elements = []
    for elem in ET.ElTable:
        Elements.append(elem[0][0])
    if len(El) in [2,4]:
        ElS = El.upper()[:2].rjust(2)
    else:
        ElS = El.upper()[:1].rjust(2)
    filename = os.path.join(sys.path[1],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        wx.MessageBox(message="File atmdata.dat not found in directory %s" % sys.path[0],
            caption="No atmdata.dat file",style=wx.OK | wx.ICON_EXCLAMATION | wx.STAY_ON_TOP)
        sys.exit()
    S = '1'
    AtomInfo = {}
    Mass = []
    while S:
        S = FFdata.readline()
        if S[3:5] == ElS:
            if S[5:6] == '_':
                if not Mass:                                 #picks 1st one; natural abundance or 1st isotope
                    Mass = float(S[10:19])
                if S[5:9] == '_SIZ':
                    Z=int(S[:2])
                    Symbol = S[3:5].strip().lower().capitalize()
                    Drad = float(S[12:22])
                    Arad = float(S[22:32])
                    Vdrad = float(S[32:38])
                    Color = ET.ElTable[Elements.index(Symbol)][6]
    FFdata.close()
    AtomInfo={'Symbol':Symbol,'Mass':Mass,'Z':Z,'Drad':Drad,'Arad':Arad,'Vdrad':Vdrad,'Color':Color}    
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
    filename = os.path.join(sys.path[1],'Xsect.dat')
    try:
        xsec = open(filename,'Ur')
    except:
        wx.MessageBox(message="File Xsect.dat not found in directory %s" % sys.path[0],
            caption="No Xsect.dat file",style=wx.OK | wx.ICON_EXCLAMATION |wx.STAY_ON_TOP)
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
    filename = os.path.join(sys.path[1],'atmdata.dat')
    try:
        FFdata = open(filename,'Ur')
    except:
        wx.MessageBox(message="File atmdata.dat not found in directory %s" % sys.path[0],
            caption="No atmdata.dat file",style=wx.OK | wx.ICON_EXCLAMATION |wx.STAY_ON_TOP)
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

def ScatFac(FormFac, SQ):
    """compute value of form factor
    @param FormFac: dictionary  defined in GetFormFactorCoeff 
    @param SQ: (sin-theta/lambda)**2
    @return: f: real part of form factor
    """
    fa = np.array(FormFac['fa'])
    fb = np.array(FormFac['fb'])
    t = -fb*SQ
    return np.sum(fa*np.exp(t))+FormFac['fc']
    
def ComptonFac(ComptonCoeff,SQ):
    """compute Compton scattering factor
    @param ComptonCoeff: list [Z, a1:a5, b1:b5]
    @param SQ: (sin-theta/lambda)**2
    @return: comp: compton scattering factor
    """    
    ca = np.array(ComptonCoeff[1:6])
    cb = np.array(ComptonCoeff[6:11])
    t = -cb*SQ
    return ComptonCoeff[0]-np.sum(ca*np.exp(t)) 
            
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
    

class PickElement(wx.Dialog):
    "Makes periodic table widget for picking element - caller maintains element list"
    Elem=None
    def _init_ctrls(self, prnt,oneOnly):
        wx.Dialog.__init__(self, id=-1, name='PickElement',
              parent=prnt, pos=wx.DefaultPosition, 
              style=wx.DEFAULT_DIALOG_STYLE, title='Pick Element')
        import ElementTable as ET
        self.SetClientSize(wx.Size(770, 250))
        
        i=0
        for E in ET.ElTable:
            if oneOnly:
                color=E[4]
            else:
                color=E[6]
            PickElement.ElButton(self,name=E[0],
               pos=wx.Point(E[1]*40+25,E[2]*24+24),tip=E[3],color=color,oneOnly=oneOnly)
            i+=1

    def __init__(self, parent,oneOnly=False):
        self._init_ctrls(parent,oneOnly)
        
    def ElButton(self, name, pos, tip, color, oneOnly):
        Black = wx.Colour(0,0,0)
        if oneOnly:
            El = wscs.ColourSelect(label=name[0], parent=self,colour=color,
                pos=pos, size=wx.Size(40,23), style=wx.RAISED_BORDER)
#            El.SetLabel(name)
            El.Bind(wx.EVT_BUTTON, self.OnElButton)
        else:
            El = wx.ComboBox(choices=name, parent=self, pos=pos, size=wx.Size(40,23),
                style=wx.CB_READONLY, value=name[0])
            El.Bind(wx.EVT_COMBOBOX,self.OnElButton)
        
        El.SetBackgroundColour(color)
        El.SetToolTipString(tip)

    def OnElButton(self, event):
        El = event.GetEventObject().GetLabel()
        self.Elem = El
        self.EndModal(wx.ID_OK)        
        
class DeleteElement(wx.Dialog):
    "Delete element from selected set widget"
    def _init_ctrls(self, parent):
        l = len(DeleteElement.Elems)-1
        wx.Dialog.__init__(self, id=-1, name='Delete', parent=parent,
              pos=wx.DefaultPosition, size=wx.Size(max(128,64+l*24), 87),
              style=wx.DEFAULT_DIALOG_STYLE, title='Delete Element')
        self.Show(True)
        self.SetAutoLayout(True)
        self.SetHelpText('Select element to delete')
        self.SetWindowVariant(wx.WINDOW_VARIANT_SMALL)

        i = 0
        Elem = []
        for Elem in DeleteElement.Elems:
            name = Elem[0].lower().capitalize()
            self.ElButton(id=-1,name=name,pos=wx.Point(16+i*24, 16))
            i+=1
              
    def __init__(self, parent):
        DeleteElement.Elems = parent.Elems
        DeleteElement.El = ' '
        self._init_ctrls(parent)

    def ElButton(self, id, name, pos):
        White = wx.Colour(255, 255, 255)
        El = wscs.ColourSelect(label=name, parent=self, colour = White,
            pos=pos, size=wx.Size(24, 23), style=wx.RAISED_BORDER)
        El.Bind(wx.EVT_BUTTON, self.OnDeleteButton)
    
    def OnDeleteButton(self, event):
        DeleteElement.El=event.GetEventObject().GetLabel()
        self.EndModal(wx.ID_OK)
        
    def GetDeleteElement(self):
        return DeleteElement.El
        

