# -*- coding: utf-8 -*-
"""
*GSASIIspc: Space group module*
-------------------------------

Space group interpretation routines. Note that space group information is
stored in a :ref:`Space Group (SGData)<SGData_table>` object.

"""
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import scipy.optimize as so
import math
import sys
import os.path as ospath

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import pyspg

npsind = lambda x: np.sin(x*np.pi/180.)
npcosd = lambda x: np.cos(x*np.pi/180.)
    
################################################################################
#### Space group codes
################################################################################

def SpcGroup(SGSymbol):
    """
    Determines cell and symmetry information from a short H-M space group name

    :param SGSymbol: space group symbol (string) with spaces between axial fields
    :returns: (SGError,SGData)
       * SGError = 0 for no errors; >0 for errors (see SGErrors below for details)
       * SGData - is a dict (see :ref:`Space Group object<SGData_table>`) with entries:
       
             * 'SpGrp': space group symbol, slightly cleaned up
             * 'SGLaue':  one of '-1', '2/m', 'mmm', '4/m', '4/mmm', '3R',
               '3mR', '3', '3m1', '31m', '6/m', '6/mmm', 'm3', 'm3m'
             * 'SGInv': boolean; True if centrosymmetric, False if not
             * 'SGLatt': one of 'P', 'A', 'B', 'C', 'I', 'F', 'R'
             * 'SGUniq': one of 'a', 'b', 'c' if monoclinic, '' otherwise
             * 'SGCen': cell centering vectors [0,0,0] at least
             * 'SGOps': symmetry operations as [M,T] so that M*x+T = x'
             * 'SGSys': one of 'triclinic', 'monoclinic', 'orthorhombic',
               'tetragonal', 'rhombohedral', 'trigonal', 'hexagonal', 'cubic'
             * 'SGPolax': one of '', 'x', 'y', 'x y', 'z', 'x z', 'y z',
               'xyz', '111' for arbitrary axes
             * 'SGPtGrp': one of 32 point group symbols (with some permutations)
                - filled by SGPtGroup - is external (KE) part of supersymmetry point group
             * 'SSGKl': default internal (Kl) part of supersymmetry point group; modified 
             in supersymmetry stuff depending on chosen modulation vector for Mono & Ortho

    """
    LaueSym = ('-1','2/m','mmm','4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm','m3','m3m')
    LattSym = ('P','A','B','C','I','F','R')
    UniqSym = ('','','a','b','c','',)
    SysSym = ('triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','trigonal','hexagonal','cubic')
    SGData = {}
    SGInfo = pyspg.sgforpy(SGSymbol)
    SGData['SpGrp'] = SGSymbol.strip().lower().capitalize()
    SGData['SGLaue'] = LaueSym[SGInfo[0]-1]
    SGData['SGInv'] = bool(SGInfo[1])
    SGData['SGLatt'] = LattSym[SGInfo[2]-1]
    SGData['SGUniq'] = UniqSym[SGInfo[3]+1]
    if SGData['SGLatt'] == 'P':
        SGData['SGCen'] = np.array(([0,0,0],))
    elif SGData['SGLatt'] == 'A':
        SGData['SGCen'] = np.array(([0,0,0],[0,.5,.5]))
    elif SGData['SGLatt'] == 'B':
        SGData['SGCen'] = np.array(([0,0,0],[.5,0,.5]))
    elif SGData['SGLatt'] == 'C':
        SGData['SGCen'] = np.array(([0,0,0],[.5,.5,0,]))
    elif SGData['SGLatt'] == 'I':
        SGData['SGCen'] = np.array(([0,0,0],[.5,.5,.5]))
    elif SGData['SGLatt'] == 'F':
        SGData['SGCen'] = np.array(([0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0,]))
    elif SGData['SGLatt'] == 'R':
        SGData['SGCen'] = np.array(([0,0,0],[1./3.,2./3.,2./3.],[2./3.,1./3.,1./3.]))
    SGData['SGOps'] = []
    for i in range(SGInfo[5]):
        Mat = np.array(SGInfo[6][i])
        Trns = np.array(SGInfo[7][i])
        SGData['SGOps'].append([Mat,Trns])
    if SGData['SGLaue'] in '-1':
        SGData['SGSys'] = SysSym[0]
    elif SGData['SGLaue'] in '2/m':
        SGData['SGSys'] = SysSym[1]
    elif SGData['SGLaue'] in 'mmm':
        SGData['SGSys'] = SysSym[2]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        SGData['SGSys'] = SysSym[3]
    elif SGData['SGLaue'] in ['3R','3mR']:
        SGData['SGSys'] = SysSym[4]
    elif SGData['SGLaue'] in ['3','3m1','31m']:
        SGData['SGSys'] = SysSym[5]
    elif SGData['SGLaue'] in ['6/m','6/mmm']:
        SGData['SGSys'] = SysSym[6]
    elif SGData['SGLaue'] in ['m3','m3m']:
        SGData['SGSys'] = SysSym[7]
    SGData['SGPolax'] = SGpolar(SGData)
    SGData['SGPtGrp'],SGData['SSGKl'] = SGPtGroup(SGData)
    return SGInfo[8],SGData

def SGErrors(IErr):
    '''
    Interprets the error message code from SpcGroup. Used in SpaceGroup.
    
    :param IErr: see SGError in :func:`SpcGroup`
    :returns:
        ErrString - a string with the error message or "Unknown error"
    '''

    ErrString = [' ',
        'Less than 2 operator fields were found',
        'Illegal Lattice type, not P, A, B, C, I, F or R',
        'Rhombohedral lattice requires a 3-axis',
        'Minus sign does not preceed 1, 2, 3, 4 or 6',
        'Either a 5-axis anywhere or a 3-axis in field not allowed',
        ' ',
        'I for COMPUTED GO TO out of range.',
        'An a-glide mirror normal to A not allowed',
        'A b-glide mirror normal to B not allowed',
        'A c-glide mirror normal to C not allowed',
        'D-glide in a primitive lattice not allowed',
        'A 4-axis not allowed in the 2nd operator field',
        'A 6-axis not allowed in the 2nd operator field',
        'More than 24 matrices needed to define group',
        ' ',
        'Improper construction of a rotation operator',
        'Mirror following a / not allowed',
        'A translation conflict between operators',
        'The 2bar operator is not allowed',
        '3 fields are legal only in R & m3 cubic groups',
        'Syntax error. Expected I -4 3 d at this point',
        ' ',
        'A or B centered tetragonal not allowed',
        ' ','unknown error in sgroup',' ',' ',' ',
        'Illegal character in the space group symbol',
        ]
    try:
        return ErrString[IErr]
    except:
        return "Unknown error"

def SGpolar(SGData):
    '''
    Determine identity of polar axes if any
    '''
    POL = ('','x','y','x y','z','x z','y z','xyz','111')
    NP = [1,2,4]
    NPZ = [0,1]
    for M,T in SGData['SGOps']:
        for i in range(3):
            if M[i][i] <= 0.: NP[i] = 0
        if M[0][2] > 0: NPZ[0] = 8
        if M[1][2] > 0: NPZ[1] = 0
    NPol = (NP[0]+NP[1]+NP[2]+NPZ[0]*NPZ[1])*(1-int(SGData['SGInv']))
    return POL[NPol]
    
def SGPtGroup(SGData):
    '''
    Determine point group of the space group - done after space group symbol has
    been evaluated by SpcGroup. Only short symbols are allowed
    
    :param SGData: from :func SpcGroup
    returns SSGPtGrp & SSGKl (only defaults for Mono & Ortho)
    '''
    Flds = SGData['SpGrp'].split(' ')
    if SGData['SGLaue'] == '-1':    #triclinic
        if '-' in Flds[1]:
            return '-1',[-1,]
        else:
            return '1',[1,]
    elif SGData['SGLaue'] == '2/m': #monoclinic - default for 2D modulation vector
        if '/' in SGData['SpGrp']:
            return '2/m',[-1,1]
        elif '2' in SGData['SpGrp']:
            return '2',[-1,]
        else:
            return 'm',[1,]
    elif SGData['SGLaue'] == 'mmm': #orthorhombic
        if SGData['SpGrp'].count('2') == 3:
            return '222',[1,1,1]
        elif SGData['SpGrp'].count('2') == 1:
            if SGData['SGPolax'] == 'x':
                return '2mm',[1,1,1]
            elif SGData['SGPolax'] == 'y':
                return 'm2m',[1,1,1]
            elif SGData['SGPolax'] == 'z':
                return 'mm2',[1,1,1]
        else:
            return 'mmm',[1,1,-1]
    elif SGData['SGLaue'] == '4/m': #tetragonal
        if '/' in SGData['SpGrp']:
            return '4/m',[1,-1]
        elif '-' in Flds[1]:
            return '-4',[-1,]
        else:
            return '4',[1,]
    elif SGData['SGLaue'] == '4/mmm':
        if '/' in SGData['SpGrp']:
            return '4/mmm',[1,-1,1,1]
        elif '-' in Flds[1]:
            if '2' in Flds[2]:
                return '-42m',[-1,-1,1]
            else:
                return '-4m2',[-1,1,-1]              
        elif '2' in Flds[2:]:
            return '422',[1,-1,-1]
        else:
            return '4mm',[1,1,1]
    elif SGData['SGLaue'] in ['3','3R']:  #trigonal/rhombohedral
        if '-' in Flds[1]:
            return '-3',[-1,]
        else:
            return '3',[1,]
    elif SGData['SGLaue'] == '3mR' or 'R' in Flds[0]:
        if '2' in Flds[2]:
            return '32',[1,-1]
        elif '-' in Flds[1]:
            return '-3m',[-1,1]
        else:
            return '3m',[1,1]
    elif SGData['SGLaue'] == '3m1':
        if '2' in Flds[2]:
            return '321',[1,-1,1]
        elif '-' in Flds[1]:
            return '-3m1',[-1,1,1]
        else:
            return '3m1',[1,1,1]
    elif SGData['SGLaue'] == '31m':
        if '2' in Flds[3]:
            return '312',[1,1,-1]
        elif '-' in Flds[1]:
            return '-31m',[-1,1,1]
        else:
            return '31m',[1,1,1]
    elif SGData['SGLaue'] == '6/m': #hexagonal
        if '/' in SGData['SpGrp']:
            return '6/m',[1,-1]
        elif '-' in SGData['SpGrp']:
            return '-6',[-1,]
        else:
            return '6',[1,]
    elif SGData['SGLaue'] == '6/mmm':
        if '/' in SGData['SpGrp']:
            return '6/mmm',[1,-1,1,1]
        elif '-' in Flds[1]:
            if '2' in Flds[2]:
                return '-62m',[-1,-1,1]
            else:
                return '-6m2',[-1,1,-1]                 
        elif '2' in Flds[2:]:
            return '622',[1,-1,-1]
        else:
            return '6mm',[1,1,1]   
    elif SGData['SGLaue'] == 'm3':      #cubic - no (3+1) supersymmetry
        if '2' in Flds[1]:
            return '23',[]
        else:  
            return 'm3',[]
    elif SGData['SGLaue'] == 'm3m':
        if '4' in Flds[1]:
            if '-' in Flds[1]:
                return '-43m',[]
            else:
                return '432',[]
        else:
            return 'm-3m',[]
    
def SGPrint(SGData):
    '''
    Print the output of SpcGroup in a nicely formatted way. Used in SpaceGroup

    :param SGData: from :func:`SpcGroup`
    :returns:
        SGText - list of strings with the space group details
        SGTable - list of strings for each of the operations
    '''
    Mult = len(SGData['SGCen'])*len(SGData['SGOps'])*(int(SGData['SGInv'])+1)
    SGText = []
    SGText.append(' Space Group: '+SGData['SpGrp'])
    CentStr = 'centrosymmetric'
    if not SGData['SGInv']:
        CentStr = 'non'+CentStr
    if SGData['SGLatt'] in 'ABCIFR':
        SGText.append(' The lattice is '+CentStr+' '+SGData['SGLatt']+'-centered '+SGData['SGSys'].lower())
    else:
        SGText.append(' The lattice is '+CentStr+' '+'primitive '+SGData['SGSys'].lower()) 
    SGText.append(' The Laue symmetry is '+SGData['SGLaue'])
    if 'SGPtGrp' in SGData:         #patch
        SGText.append(' The lattice point group is '+SGData['SGPtGrp'])
    SGText.append(' Multiplicity of a general site is '+str(Mult))
    if SGData['SGUniq'] in ['a','b','c']:
        SGText.append(' The unique monoclinic axis is '+SGData['SGUniq'])
    if SGData['SGInv']:
        SGText.append(' The inversion center is located at 0,0,0')
    if SGData['SGPolax']:
        SGText.append(' The location of the origin is arbitrary in '+SGData['SGPolax'])
    SGText.append(' ')
    if SGData['SGLatt'] == 'P':
        SGText.append(' The equivalent positions are:\n')
    else:    
        SGText.append(' The equivalent positions are:')
        SGText.append(' ('+Latt2text(SGData['SGLatt'])+')+\n')
    SGTable = []
    for i,Opr in enumerate(SGData['SGOps']):
        SGTable.append('(%2d) %s'%(i+1,MT2text(Opr)))
    return SGText,SGTable

def AllOps(SGData):
    '''
    Returns a list of all operators for a space group, including those for
    centering and a center of symmetry
    
    :param SGData: from :func:`SpcGroup`
    :returns: (SGTextList,offsetList,symOpList,G2oprList) where

      * SGTextList: a list of strings with formatted and normalized
        symmetry operators.
      * offsetList: a tuple of (dx,dy,dz) offsets that relate the GSAS-II
        symmetry operation to the operator in SGTextList and symOpList.
        these dx (etc.) values are added to the GSAS-II generated
        positions to provide the positions that are generated
        by the normalized symmetry operators.        
      * symOpList: a list of tuples with the normalized symmetry
        operations as (M,T) values
        (see ``SGOps`` in the :ref:`Space Group object<SGData_table>`)
      * G2oprList: The GSAS-II operations for each symmetry operation as
        a tuple with (center,mult,opnum), where center is (0,0,0), (0.5,0,0),
        (0.5,0.5,0.5),...; where mult is 1 or -1 for the center of symmetry
        and opnum is the number for the symmetry operation, in ``SGOps``
        (starting with 0).
    '''
    SGTextList = []
    offsetList = []
    symOpList = []
    G2oprList = []
    onebar = (1,)
    if SGData['SGInv']:
        onebar += (-1,)
    for cen in SGData['SGCen']:
        for mult in onebar:
            for j,(M,T) in enumerate(SGData['SGOps']):
                offset = [0,0,0]
                Tprime = (mult*T)+cen
                for i in range(3):
                    while Tprime[i] < 0:
                        Tprime[i] += 1
                        offset[i] += 1
                    while Tprime[i] >= 1:
                        Tprime[i] += -1
                        offset[i] += -1
                Opr = [mult*M,Tprime]
                OPtxt = MT2text(Opr)
                SGTextList.append(OPtxt.replace(' ',''))
                offsetList.append(tuple(offset))
                symOpList.append((mult*M,Tprime))
                G2oprList.append((cen,mult,j))
    return SGTextList,offsetList,symOpList,G2oprList
    
def MT2text(Opr):
    "From space group matrix/translation operator returns text version"
    XYZ = ('-Z','-Y','-X','X-Y','ERR','Y-X','X','Y','Z')
    TRA = ('   ','ERR','1/6','1/4','1/3','ERR','1/2','ERR','2/3','3/4','5/6','ERR')
    Fld = ''
    M,T = Opr
    for j in range(3):
        IJ = int(round(2*M[j][0]+3*M[j][1]+4*M[j][2]+4))%12
        IK = int(round(T[j]*12))%12
        if IK:
            if IJ < 3:
                Fld += (TRA[IK]+XYZ[IJ]).rjust(5)
            else:
                Fld += (TRA[IK]+'+'+XYZ[IJ]).rjust(5)
        else:
            Fld += XYZ[IJ].rjust(5)
        if j != 2: Fld += ', '
    return Fld
    
def Latt2text(Latt):
    "From lattice type ('P',A', etc.) returns ';' delimited cell centering vectors"
    lattTxt = {'A':'0,0,0; 0,1/2,1/2','B':'0,0,0; 1/2,0,1/2',
        'C':'0,0,0; 1/2,1/2,0','I':'0,0,0; 1/2,1/2,1/2',
        'F':'0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0',
        'R':'0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3','P':'0,0,0'}
    return lattTxt[Latt]    
        
def SpaceGroup(SGSymbol):
    '''
    Print the output of SpcGroup in a nicely formatted way. 

    :param SGSymbol: space group symbol (string) with spaces between axial fields
    :returns: nothing
    '''
    E,A = SpcGroup(SGSymbol)
    if E > 0:
        print SGErrors(E)
        return
    for l in SGPrint(A):
        print l
        
################################################################################
#### Superspace group codes
################################################################################
        
def SSpcGroup(SGData,SSymbol):
    """
    Determines supersymmetry information from superspace group name; currently only for (3+1) superlattices

    :param SGData: space group data structure as defined in SpcGroup above.
    :param SSymbol: superspace group symbol extension (string) defining modulation direction & generator info.
    :returns: (SSGError,SSGData)
       * SGError = 0 for no errors; >0 for errors (see SGErrors below for details)
       * SSGData - is a dict (see :ref:`Superspace Group object<SSGData_table>`) with entries:
       
             * 'SSpGrp': superspace group symbol extension to space group symbol, accidental spaces removed
             * 'SSGCen': 4D cell centering vectors [0,0,0,0] at least
             * 'SSGOps': 4D symmetry operations as [M,T] so that M*x+T = x'

    """
    
    def splitSSsym(SSymbol):
        '''
        Splits supersymmetry symbol into two lists of strings
        '''
        modsym,gensym = SSymbol.replace(' ','').split(')')
        nfrac = modsym.count('/')
        modsym = modsym.lstrip('(')
        if nfrac == 0:
            modsym = list(modsym)
        elif nfrac == 1:
            pos = modsym.find('/')
            if pos == 1:
                modsym = [modsym[:3],modsym[3],modsym[4]]
            elif pos == 2:
                modsym = [modsym[0],modsym[1:4],modsym[4]]
            else:
                modsym = [modsym[0],modsym[1],modsym[2:]]
        else:
            lpos = modsym.find('/')
            rpos = modsym.rfind('/')
            if lpos == 1 and rpos == 4:
                modsym = [modsym[:3],modsym[3:6],modsym[6]]
            elif lpos == 1 and rpos == 5:
                modsym = [modsym[:3],modsym[3],modsym[4:]]
            else:
                modsym = [modsym[0],modsym[1:4],modsym[4:]]
        gensym = list(gensym)
        return modsym,gensym
        
    def checkModSym():
        ''' 
        Checks to see if proposed modulation form is allowed for Laue group
        '''
        if LaueId in [0,] and LaueModId in [0,]:
            return True
        elif LaueId in [1,]:
            try:
                if modsym.index('1/2') != ['A','B','C'].index(SGData['SGLatt']):
                    return False
                if 'I'.index(SGData['SGLatt']) and modsym.count('1/2') not in [0,2]:
                    return False
            except ValueError:
                pass
            if SGData['SGUniq'] == 'a' and LaueModId in [5,6,7,8,9,]:
                return True
            elif SGData['SGUniq'] == 'b' and LaueModId in [3,4,13,14,15,]:
                return True
            elif SGData['SGUniq'] == 'c' and LaueModId in [1,2,19,20,21,]:
                return True
        elif LaueId in [2,] and LaueModId in [i+7 for i in range(18)]:
            try:
                if modsym.index('1/2') != ['A','B','C'].index(SGData['SGLatt']):
                    return False
                if SGData['SGLatt'] in ['I','F',] and modsym.index('1/2'):
                    return False
            except ValueError:
                pass
            return True
        elif LaueId in [3,4,] and LaueModId in [19,22,]:
            try:
                if SGData['SGLatt'] == 'I' and modsym.count('1/2'):
                    return False
            except ValueError:
                pass
            return True
        elif LaueId in [7,8,9,] and LaueModId in [19,25,]:
            if (SGData['SGLatt'] == 'R' or SGData['SGPtGrp'] in ['3m1','-3m1']) and modsym.count('1/3'):
                return False
            return True
        elif LaueId in [10,11,] and LaueModId in [19,]:
            return True
        return False
        
    def fixMonoOrtho():
        mod = ''.join(modsym).replace('1/2','0').replace('1','0')
        if SGData['SGPtGrp'] in ['2','m']:  #OK
            if mod in ['a00','0b0','00g']:
                result = [i*-1 for i in SGData['SSGKl']]
            else:
                result = SGData['SSGKl'][:]
            if '/' in mod:
                return [i*-1 for i in result]
            else:
                return result
        elif SGData['SGPtGrp'] == '2/m':    #OK
            if mod in ['a00','0b0','00g']:
                result =  SGData['SSGKl'][:]
            else:
                result = [i*-1 for i in SGData['SSGKl']]
            if '/' in mod:
                return [i*-1 for i in result]
            else:
                return result
        else:   #orthorhombic
            if SGData['SGPtGrp'] == '222':
                return [1 if i in ['a','b','g'] else -1 for i in mod]
            elif SGData['SGPtGrp'] == 'mm2':
                if 'g' in mod:
                    return [1,1,1]
                elif 'b' in mod:
                    return [1,-1,-1]
                else:
                    return [-1,1,-1]
            elif SGData['SGPtGrp'] == 'm2m':
                if 'b' in mod:
                    return [1,1,1]
                elif 'g' in mod:
                    return [1,-1,-1]
                else:
                    return [-1,-1,1]                
            elif SGData['SGPtGrp'] == '2mm':
                if 'a' in mod:
                    return [1,1,1]
                elif 'b' in mod:
                    return [-1,-1,1]
                else:
                    return [-1,1,-1]
            else:
                return [-1 if i in ['a','b','g'] else 1 for i in mod]
                
    def extendSSGOps(SSGOps):
        nOps = len(SSGOps)
        for i in range(nOps):
            if np.allclose(SSGOps[i][0][3],np.zeros(4)):
                continue
            for j in range(nOps):
                if np.allclose(SSGOps[j][0][3],np.zeros(4)):
                    continue
                OpC = list(SGProd(SSGOps[j],SSGOps[i]))
                OpC[1] %= 1.
                for k in range(nOps):
                    OpD = SSGOps[k]
                    if SSMT2text(OpC) == SSMT2text(OpD):
                        continue
                    elif np.allclose(OpC[0][:3,:3],OpD[0][:3,:3]):
                        if np.allclose(OpD[0][3],np.zeros(4)):
                            SSGOps[k] = OpC
                        elif np.any([np.allclose(OpC[0][3][:3],cen) for cen in SGData['SGCen']]):   #?
                            continue
                        else:
                            OpCtxt = SSMT2text(OpC).replace(' ','')
                            OpDtxt = SSMT2text(OpD).replace(' ','')
                            print 'OpC',OpCtxt,'OpD',OpDtxt
                            return False,OpCtxt+' conflict with '+OpDtxt
        return True,SSGOps
                
    def genSSGOps():
        SSGOps = SSGData['SSGOps'][:]
        iFrac = {}
        for i,frac in enumerate(SSGData['modSymb']):
            if frac in ['1/2','1/3','1/4','1/6','1']:
                iFrac[i] = frac
        print SGData['SpGrp']+SSymbol
        print 'SSGKl',SSGKl,'genQ',genQ,'iFrac',iFrac
# set identity & 1,-1; triclinic
        SSGOps[0][0][3,3] = 1.
# expand if centrosymmetric
        if SGData['SGInv']:
            SSGOps += [[-1*M,V] for M,V in SSGOps[:]]
# monoclinic - all done
        if SGData['SGPtGrp'] in ['2','m']:  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
            for i in iFrac:
                SSGOps[1][0][3,i] = -SSGKl[0]
        elif SGData['SGPtGrp'] == '2/m':    #OK
            for i,j in enumerate([1,3]):
                SSGOps[j][0][3,3] = -SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                for k in iFrac:
                    SSGOps[j][0][3,k] = SSGKl[i]
                E,SSGOps = extendSSGOps(SSGOps)
            
# orthorhombic
        elif SGData['SGPtGrp'] in ['222','mm2','m2m','2mm','mmm']:
            for i in [0,1,2]:
                SSGOps[i+1][0][3,3] = SSGKl[i]
                SSGOps[i+1][1][3] = genQ[i]
            for i in iFrac:
                SSGOps[1][0][3,i] = -1
            print SSMT2text(SSGOps[1]).replace(' ',''),SSMT2text(SSGOps[2]).replace(' ',''), \
                SSMT2text(SSGOps[3]).replace(' ','')
                
# tetragonal
        elif SGData['SGPtGrp'] == '4':  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
        elif SGData['SGPtGrp'] == '-4': #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = 1
        elif SGData['SGPtGrp'] in ['4/m',]:
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            for i,j in enumerate([1,6]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['422','4mm','-42m','-4m2',]:
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            for i,j in enumerate([1,4,5]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['4/mmm',]:
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            for i,j in enumerate([1,10,6,7]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
#                for k in iFrac:
#                    SSGOps[j][0][3,k] = SSGKl[i]
                E,SSGOps = extendSSGOps(SSGOps)
                
# trigonal - all done
        elif SGData['SGPtGrp'] == '3':  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            if '1/3' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            SSGOps[1][1][3] = genQ[0]
        elif SGData['SGPtGrp'] == '-3': #OK
            SSGOps[1][0][3,3] = -SSGKl[0]
            if '1/3' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            SSGOps[1][1][3] = genQ[0]
        elif SGData['SGPtGrp'] in ['312','3m','-3m','-3m1','3m1']:   #OK
            if '1/3' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
            for i,j in enumerate([1,5]):
                if SGData['SGPtGrp'] in ['3m','-3m']:
                    SSGOps[j][0][3,3] = SSGKl[i]
                else:                    
                    SSGOps[j][0][3,3] = SSGKl[i+1]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
        elif SGData['SGPtGrp'] in ['321','32']:   #OK
            for i,j in enumerate([1,4]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
        elif SGData['SGPtGrp'] in ['31m','-31m']:   #OK
            ids = [1,3]
            if SGData['SGPtGrp'] == '-31m':
                ids = [7,3]
            if '1/3' in SSGData['modSymb']:
                SSGOps[ids[0]][0][3,1] = -1
            for i,j in enumerate(ids):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i+1]:
                    SSGOps[j][1][3] = genQ[i+1]
                     
# hexagonal - all done
        elif SGData['SGPtGrp'] == '6':  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
        elif SGData['SGPtGrp'] == '-6': #OK
            SSGOps[1][0][3,3] = SSGKl[0]
        elif SGData['SGPtGrp'] in ['6/m',]: #OK
            SSGOps[1][0][3,3] = -SSGKl[1]
            SSGOps[1][1][3] = genQ[0]
            SSGOps[2][1][3] = genQ[1]
        elif SGData['SGPtGrp'] in ['622','6mm','-62m','-6m2',]: #OK
            for i,j in enumerate([1,10,11]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['6/mmm',]: #OK
            for i,j in enumerate([1,15,10,11]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['1','-1']: #triclinic - done
            return True,SSGOps
        E,SSGOps = extendSSGOps(SSGOps)
        return E,SSGOps
        
    def specialGen(gensym):
        sym = ''.join(gensym)
        if SGData['SGPtGrp'] in ['2/m',] and 'n' in SGData['SpGrp']:
            if 's' in sym:
                gensym = 'ss'
        if SGData['SGPtGrp'] in ['-62m',] and sym == '00s':
            gensym = '0ss'
        elif SGData['SGPtGrp'] in ['222',]:
            if sym == '00s':
                gensym = '0ss'
            elif sym == '0s0':
                gensym = 'ss0'
            elif sym == 's00':
                gensym = 's0s'
        return gensym
                    
    def checkGen(gensym):
        sym = ''.join(gensym)
        print str(SSGKl),sym
# monoclinic - all done
        if str(SSGKl) == '[-1]' and sym == 's':
            return False
        elif SGData['SGPtGrp'] in ['2/m',]:
            if str(SSGKl) == '[-1, 1]' and sym == '0s':
                return False
            elif str(SSGKl) == '[1, -1]' and sym == 's0':
                return False
#orthorhombic - all 
        elif SGData['SGPtGrp'] in ['222',] and sym not in ['','s00','0s0','00s']:
            return False 
        elif SGData['SGPtGrp'] in ['2mm','m2m','mm2','mmm'] and sym not in GenSymList[4:15]:
            return False 
#tetragonal - all done
        elif SGData['SGPtGrp'] in ['4',] and sym not in ['','s','q']:
            return False 
        elif SGData['SGPtGrp'] in ['-4',] and sym not in ['',]:
            return False             
        elif SGData['SGPtGrp'] in ['4/m',] and sym not in ['','s0','q0']:
            return False
        elif SGData['SGPtGrp'] in ['422',] and sym not in ['','q00','s00']:
            return False         
        elif SGData['SGPtGrp'] in ['4mm',] and sym not in ['','ss0','s0s','0ss','qq0','qqs']:
            return False
        elif SGData['SGPtGrp'] in ['-4m2',] and sym not in ['','00s','00q']:
            return False
        elif SGData['SGPtGrp'] in ['-42m',] and sym not in ['','0s0','0q0']:
            return False
        elif SGData['SGPtGrp'] in ['4/mmm',] and sym not in ['','s00s','s0s0','00ss','q0q0','q0qs']:
            return False
#trigonal/rhombohedral - all done
        elif SGData['SGPtGrp'] in ['3',] and sym not in ['','t']:
            return False 
        elif SGData['SGPtGrp'] in ['-3',] and sym not in ['',]:
            return False 
        elif SGData['SGPtGrp'] in ['32',] and sym not in ['','t0']:
            return False 
        elif SGData['SGPtGrp'] in ['321','312'] and sym not in ['','t00']:
            return False 
        elif SGData['SGPtGrp'] in ['3m','-3m'] and sym not in ['','0s']:
            return False 
        elif SGData['SGPtGrp'] in ['3m1','-3m1'] and sym not in ['','0s0']:
            return False 
        elif SGData['SGPtGrp'] in ['31m','-31m'] and sym not in ['','00s']:
            return False 
#hexagonal - all done
        elif SGData['SGPtGrp'] in ['6',] and sym not in ['','s','h','t']:
            return False 
        elif SGData['SGPtGrp'] in ['-6',] and sym not in ['',]:
            return False
        elif SGData['SGPtGrp'] in ['6/m',] and sym not in ['','s0']:
            return False
        elif SGData['SGPtGrp'] in ['622',] and sym not in ['','h00','t00','s00']:
            return False         
        elif SGData['SGPtGrp'] in ['6mm',] and sym not in ['','ss0','s0s','0ss']:
            return False
        elif SGData['SGPtGrp'] in ['-6m2',] and sym not in ['','0s0']:
            return False
        elif SGData['SGPtGrp'] in ['-62m',] and sym not in ['','0ss']:
            return False
        elif SGData['SGPtGrp'] in ['6/mmm',] and sym not in ['','s00s','s0s0','00ss']:
            return False
        return True
        
    LaueModList = ['abg', 'ab0', 'ab1/2', 'a0g', 'a1/2g','0bg', '1/2bg',
               'a00', 'a01/2', 'a1/20', 'a1/21/2', 'a01', 'a10', 
               '0b0', '0b1/2', '1/2b0', '1/2b1/2', '0b1', '1b0',
               '00g', '01/2g', '1/20g', '1/21/2g', '01g', '10g','1/31/3g']
    LaueList = ['-1','2/m','mmm','4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm','m3','m3m']
    GenSymList = ['','s','0s','s0','00s','0s0','s00','s0s','ss0','0ss','q00','0q0','00q','qq0','q0q','0qq',
        'q','q0','0q','qqs','s0s0','00ss','s00s','q0q0','q0qs','t','t00','t0','h','h00']
    Fracs = {'1/2':0.5,'1/3':1./3,'1':1.0,'0':0.,'s':.5,'t':1./3,'q':.25,'h':1./6,'a':0.,'b':0.,'g':0.}
    LaueId = LaueList.index(SGData['SGLaue'])
    if SGData['SGLaue'] in ['m3','m3m']:
        return '(3+1) superlattices not defined for cubic space groups',None
    elif SGData['SGLaue'] in ['3R','3mR']:
        return '(3+1) superlattices not defined for rhombohedral settings - use hexagonal setting',None
    try:
        modsym,gensym = splitSSsym(SSymbol)
    except ValueError:
        return 'Error in superspace symbol '+SSymbol,None
    if ''.join(gensym) not in GenSymList:
        return 'unknown generator symbol '+''.join(gensym),None
    try:
        print modsym,''.join(modsym)
        LaueModId = LaueModList.index(''.join(modsym))
    except ValueError:
        return 'Unknown modulation symbol '+''.join(modsym),None
    if not checkModSym():
        return 'Modulation '+''.join(modsym)+' not consistent with space group '+SGData['SpGrp'],None
    modQ = [Fracs[mod] for mod in modsym]
    SSGKl = SGData['SSGKl'][:]
    if SGData['SGLaue'] in ['2/m','mmm']:
        SSGKl = fixMonoOrtho()
    if len(gensym) and len(gensym) != len(SSGKl):
        return 'Wrong number of items in generator symbol '+''.join(gensym),None
    if not checkGen(gensym):
        return 'Generator '+''.join(gensym)+' not consistent with space group '+SGData['SpGrp'],None
    gensym = specialGen(gensym)
    genQ = [Fracs[mod] for mod in gensym]
    if not genQ:
        genQ = [0,0,0,0]
    SSGData = {'SSpGrp':SGData['SpGrp']+SSymbol,'modQ':modQ,'modSymb':modsym}
    SSCen = np.zeros((len(SGData['SGCen']),4))
    for icen,cen in enumerate(SGData['SGCen']):
        SSCen[icen,0:3] = cen
    SSCen[0] = np.zeros(4)
    SSGData['SSGCen'] = SSCen
    SSGData['SSGOps'] = []
    for iop,op in enumerate(SGData['SGOps']):
        T = np.zeros(4)
        ssop = np.zeros((4,4))
        ssop[:3,:3] = op[0]
        T[:3] = op[1]
        SSGData['SSGOps'].append([ssop,T])
    E,Result = genSSGOps()
    if E:
        SSGData['SSGOps'] = Result                     
        return None,SSGData
    else:
        return Result+'\nOperator conflict - incorrect superspace symbol',None

def SSGPrint(SGData,SSGData):
    '''
    Print the output of SSpcGroup in a nicely formatted way. Used in SSpaceGroup

    :param SGData: space group data structure as defined in SpcGroup above.
    :param SSGData: from :func:`SSpcGroup`
    :returns:
        SSGText - list of strings with the superspace group details
        SGTable - list of strings for each of the operations
    '''
    Mult = len(SSGData['SSGCen'])*len(SSGData['SSGOps'])
    SSGText = []
    SSGText.append(' Superspace Group: '+SSGData['SSpGrp'])
    CentStr = 'centrosymmetric'
    if not SGData['SGInv']:
        CentStr = 'non'+CentStr
    if SGData['SGLatt'] in 'ABCIFR':
        SSGText.append(' The lattice is '+CentStr+' '+SGData['SGLatt']+'-centered '+SGData['SGSys'].lower())
    else:
        SSGText.append(' The superlattice is '+CentStr+' '+'primitive '+SGData['SGSys'].lower())        
    SSGText.append(' The Laue symmetry is '+SGData['SGLaue'])
    SSGText.append(' The superlattice point group is '+SGData['SGPtGrp']+','+''.join([str(i) for i in SGData['SSGKl']]))
    SSGText.append(' The number of superspace group generators is '+str(len(SGData['SSGKl'])))
    SSGText.append(' Multiplicity of a general site is '+str(Mult))
    if SGData['SGUniq'] in ['a','b','c']:
        SSGText.append(' The unique monoclinic axis is '+SGData['SGUniq'])
    if SGData['SGInv']:
        SSGText.append(' The inversion center is located at 0,0,0')
    if SGData['SGPolax']:
        SSGText.append(' The location of the origin is arbitrary in '+SGData['SGPolax'])
    SSGText.append(' ')
    if len(SSGData['SSGCen']) > 1:
        SSGText.append(' The equivalent positions are:')
        SSGText.append(' ('+SSLatt2text(SSGData['SSGCen'])+')+\n')
    else:
        SSGText.append(' The equivalent positions are:\n')
    SSGTable = []
    for i,Opr in enumerate(SSGData['SSGOps']):
        SSGTable.append('(%2d) %s'%(i+1,SSMT2text(Opr)))
    return SSGText,SSGTable
    
def SSGModCheck(Vec,SSGData):
    ''' Checks modulation vector compatibility with supersymmetry space group symbol. 
    Superspace group symbol takes precidence & the vector will be modified accordingly
    '''
    modQ = SSGData['modQ']
    modSymb = SSGData['modSymb']
    Vec = [0.1 if (vec == 0.0 and mod in ['a','b','g']) else vec for [vec,mod] in zip(Vec,modSymb)]
    return [Q if mod not in ['a','b','g'] and vec != Q else vec for [vec,mod,Q] in zip(Vec,modSymb,modQ)]

def SSMT2text(Opr):
    "From superspace group matrix/translation operator returns text version"
    XYZS = ('x','y','z','t')    #Stokes, Campbell & van Smaalen notation
    TRA = ('   ','ERR','1/6','1/4','1/3','ERR','1/2','ERR','2/3','3/4','5/6','ERR')
    Fld = ''
    M,T = Opr
    for j in range(4):
        IJ = ''
        for k in range(4):
            txt = str(int(round(M[j][k])))
            txt = txt.replace('1',XYZS[k]).replace('0','')
            if '2' in txt:
                txt += XYZS[k]
            if IJ and M[j][k] > 0:
                IJ += '+'+txt
            else:
                IJ += txt
        IK = int(round(T[j]*12))%12
        if IK:
            if not IJ:
                break
            if IJ[0] == '-':
                Fld += (TRA[IK]+IJ).rjust(8)
            else:
                Fld += (TRA[IK]+'+'+IJ).rjust(8)
        else:
            Fld += IJ.rjust(8)
        if j != 3: Fld += ', '
    return Fld
    
def SSLatt2text(SSGCen):
    "Lattice centering vectors to text"
    lattTxt = ''
    for vec in SSGCen:
        lattTxt += ' '
        for item in vec:
            if int(item*12.):
                lattTxt += '1/%d,'%(12/int(item*12))
            else:
                lattTxt += '0,'
        lattTxt = lattTxt.rstrip(',')
        lattTxt += ';'
    lattTxt = lattTxt.rstrip(';').lstrip(' ')
    return lattTxt
        
def SSpaceGroup(SGSymbol,SSymbol):
    '''
    Print the output of SSpcGroup in a nicely formatted way. 

    :param SGSymbol: space group symbol with spaces between axial fields.
    :param SSymbol: superspace group symbol extension (string).
    :returns: nothing
    '''

    E,A = SpcGroup(SGSymbol)
    if E > 0:
        print SGErrors(E)
        return
    E,B = SSpcGroup(A,SSymbol)    
    if E > 0:
        print E
        return
    for l in SSGPrint(B):
        print l
        
def SGProd(OpA,OpB):
    '''
    Form space group operator product. OpA & OpB are [M,V] pairs; 
        both must be of same dimension (3 or 4). Returns [M,V] pair
    '''
    A,U = OpA
    B,V = OpB
    M = np.inner(B.T,A)
    W = np.inner(B,U)+V
    return M.T,W
        
def MoveToUnitCell(xyz):
    '''
    Translates a set of coordinates so that all values are >=0 and < 1 

    :param xyz: a list or numpy array of fractional coordinates
    :returns: XYZ - numpy array of new coordinates now 0 or greater and less than 1
    '''
    XYZ = np.zeros(3)
    for i,x in enumerate(xyz):
        XYZ[i] = (x-int(x))%1.0
    return XYZ
        
def Opposite(XYZ,toler=0.0002):
    '''
    Gives opposite corner, edge or face of unit cell for position within tolerance. 
        Result may be just outside the cell within tolerance 

    :param XYZ: 0 >= np.array[x,y,z] > 1 as by MoveToUnitCell
    :param toler: unit cell fraction tolerance making opposite
    :returns:
        XYZ: array of opposite positions; always contains XYZ
    '''
    perm3 = [[1,1,1],[0,1,1],[1,0,1],[1,1,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0]]
    TB = np.where(abs(XYZ-1)<toler,-1,0)+np.where(abs(XYZ)<toler,1,0)
    perm = TB*perm3
    cperm = ['%d%d%d'%(i,j,k) for i,j,k in perm]
    D = dict(zip(cperm,perm))
    new = []
    for key in D:
        new.append(np.array(D[key])+np.array(XYZ))
    return new
        
def GenAtom(XYZ,SGData,All=False,Uij=[],Move=True):
    '''
    Generates the equivalent positions for a specified coordinate and space group

    :param XYZ: an array, tuple or list containing 3 elements: x, y & z
    :param SGData: from :func:`SpcGroup`
    :param All: True return all equivalent positions including duplicates;
      False return only unique positions
    :param Uij: [U11,U22,U33,U12,U13,U23] or [] if no Uij
    :param Move: True move generated atom positions to be inside cell
      False do not move atoms       
    :return: [[XYZEquiv],Idup,[UijEquiv]]

      *  [XYZEquiv] is list of equivalent positions (XYZ is first entry)
      *  Idup = [-][C]SS where SS is the symmetry operator number (1-24), C (if not 0,0,0)
      * is centering operator number (1-4) and - is for inversion
        Cell = unit cell translations needed to put new positions inside cell
        [UijEquiv] - equivalent Uij; absent if no Uij given
        
    '''
    XYZEquiv = []
    UijEquiv = []
    Idup = []
    Cell = []
    X = np.array(XYZ)
    if Move:
        X = MoveToUnitCell(X)
    for ic,cen in enumerate(SGData['SGCen']):
        C = np.array(cen)
        for invers in range(int(SGData['SGInv']+1)):
            for io,[M,T] in enumerate(SGData['SGOps']):
                idup = ((io+1)+100*ic)*(1-2*invers)
                XT = np.inner(M,X)+T
                if len(Uij):
                    U = Uij2U(Uij)
                    U = np.inner(M,np.inner(U,M).T)
                    newUij = U2Uij(U)
                if invers:
                    XT = -XT
                XT += C
                if Move:
                    newX = MoveToUnitCell(XT)
                else:
                    newX = XT
                cell = np.asarray(np.rint(newX-XT),dtype=np.int32)
                if All:
                    if np.allclose(newX,X,atol=0.0002):
                        idup = False
                else:
                    if True in [np.allclose(newX,oldX,atol=0.0002) for oldX in XYZEquiv]:
                        idup = False
                if All or idup:
                    XYZEquiv.append(newX)
                    Idup.append(idup)
                    Cell.append(cell)
                    if len(Uij):
                        UijEquiv.append(newUij)                    
    if len(Uij):
        return zip(XYZEquiv,UijEquiv,Idup,Cell)
    else:
        return zip(XYZEquiv,Idup,Cell)

def GenHKLf(HKL,SGData):
    '''
    Uses old GSAS Fortran routine genhkl.for

    :param HKL:  [h,k,l]
    :param SGData: space group data obtained from SpcGroup
    :returns: iabsnt,mulp,Uniq,phi

     *   iabsnt = True if reflection is forbidden by symmetry
     *   mulp = reflection multiplicity including Friedel pairs
     *   Uniq = numpy array of equivalent hkl in descending order of h,k,l

    '''
    hklf = HKL+[0,]
    Ops = SGData['SGOps']
    OpM = np.array([op[0] for op in Ops])
    OpT = np.array([op[1] for op in Ops])
    Inv = SGData['SGInv']
    Cen = np.array([cen for cen in SGData['SGCen']])
    
    Nuniq,Uniq,iabsnt,mulp = pyspg.genhklpy(hklf,len(Ops),OpM,OpT,SGData['SGInv'],len(Cen),Cen)
    h,k,l,f = Uniq
    Uniq=np.array(zip(h[:Nuniq],k[:Nuniq],l[:Nuniq]))
    phi = f[:Nuniq]
    
    return iabsnt,mulp,Uniq,phi
                                  
def GetOprPtrName(key):
    'Needs a doc string'
    OprPtrName = {
        '-6643':[   2,' 1bar ', 1],'6479' :[  10,'  2z  ', 2],'-6479':[   9,'  mz  ', 3],
        '6481' :[   7,'  my  ', 4],'-6481':[   6,'  2y  ', 5],'6641' :[   4,'  mx  ', 6],
        '-6641':[   3,'  2x  ', 7],'6591' :[  28,' m+-0 ', 8],'-6591':[  27,' 2+-0 ', 9],
        '6531' :[  25,' m110 ',10],'-6531':[  24,' 2110 ',11],'6537' :[  61,'  4z  ',12],
        '-6537':[  62,' -4z  ',13],'975'  :[  68,' 3+++1',14],'6456' :[ 114,'  3z1 ',15],
        '-489' :[  73,' 3+-- ',16],'483'  :[  78,' 3-+- ',17],'-969' :[  83,' 3--+ ',18],
        '819'  :[  22,' m+0- ',19],'-819' :[  21,' 2+0- ',20],'2431' :[  16,' m0+- ',21],
        '-2431':[  15,' 20+- ',22],'-657' :[  19,' m101 ',23],'657'  :[  18,' 2101 ',24],
        '1943' :[  48,' -4x  ',25],'-1943':[  47,'  4x  ',26],'-2429':[  13,' m011 ',27],
        '2429' :[  12,' 2011 ',28],'639'  :[  55,' -4y  ',29],'-639' :[  54,'  4y  ',30],
        '-6484':[ 146,' 2010 ', 4],'6484' :[ 139,' m010 ', 5],'-6668':[ 145,' 2100 ', 6],
        '6668' :[ 138,' m100 ', 7],'-6454':[ 148,' 2120 ',18],'6454' :[ 141,' m120 ',19],
        '-6638':[ 149,' 2210 ',20],'6638' :[ 142,' m210 ',21],              #search ends here
        '2223' :[  68,' 3+++2',39],
        '6538' :[ 106,'  6z1 ',40],'-2169':[  83,' 3--+2',41],'2151' :[  73,' 3+--2',42],
        '2205' :[  79,'-3-+-2',43],'-2205':[  78,' 3-+-2',44],'489'  :[  74,'-3+--1',45],
        '801'  :[  53,'  4y1 ',46],'1945' :[  47,'  4x3 ',47],'-6585':[  62,' -4z3 ',48],
        '6585' :[  61,'  4z3 ',49],'6584' :[ 114,'  3z2 ',50],'6666' :[ 106,'  6z5 ',51],
        '6643' :[   1,' Iden ',52],'-801' :[  55,' -4y1 ',53],'-1945':[  48,' -4x3 ',54],
        '-6666':[ 105,' -6z5 ',55],'-6538':[ 105,' -6z1 ',56],'-2223':[  69,'-3+++2',57],
        '-975' :[  69,'-3+++1',58],'-6456':[ 113,' -3z1 ',59],'-483' :[  79,'-3-+-1',60],
        '969'  :[  84,'-3--+1',61],'-6584':[ 113,' -3z2 ',62],'2169' :[  84,'-3--+2',63],
        '-2151':[  74,'-3+--2',64],'0':[0,' ????',0]
        }
    return OprPtrName[key]

def GetKNsym(key):
    'Needs a doc string'
    KNsym = {
        '0'         :'    1   ','1'         :'   -1   ','64'        :'    2(x)','32'        :'    m(x)',
        '97'        :'  2/m(x)','16'        :'    2(y)','8'         :'    m(y)','25'        :'  2/m(y)',
        '2'         :'    2(z)','4'         :'    m(z)','7'         :'  2/m(z)','134217728' :'   2(yz)',
        '67108864'  :'   m(yz)','201326593' :' 2/m(yz)','2097152'   :'  2(0+-)','1048576'   :'  m(0+-)',
        '3145729'   :'2/m(0+-)','8388608'   :'   2(xz)','4194304'   :'   m(xz)','12582913'  :' 2/m(xz)',
        '524288'    :'  2(+0-)','262144'    :'  m(+0-)','796433'    :'2/m(+0-)','1024'      :'   2(xy)',
        '512'       :'   m(xy)','1537'      :' 2/m(xy)','256'       :'  2(+-0)','128'       :'  m(+-0)',
        '385'       :'2/m(+-0)','76'        :'  mm2(x)','52'        :'  mm2(y)','42'        :'  mm2(z)',
        '135266336' :' mm2(yz)','69206048'  :'mm2(0+-)','8650760'   :' mm2(xz)','4718600'   :'mm2(+0-)',
        '1156'      :' mm2(xy)','772'       :'mm2(+-0)','82'        :'  222   ','136314944' :'  222(x)',
        '8912912'   :'  222(y)','1282'      :'  222(z)','127'       :'  mmm   ','204472417' :'  mmm(x)',
        '13369369'  :'  mmm(y)','1927'      :'  mmm(z)','33554496'  :'  4(100)','16777280'  :' -4(100)',
        '50331745'  :'4/m(100)','169869394' :'422(100)','84934738'  :'-42m 100','101711948' :'4mm(100)',
        '254804095' :'4/mmm100','536870928 ':'  4(010)','268435472' :' -4(010)','805306393' :'4/m (10)',
        '545783890' :'422(010)','272891986' :'-42m 010','541327412' :'4mm(010)','818675839' :'4/mmm010',
        '2050'      :'  4(001)','4098'      :' -4(001)','6151'      :'4/m(001)','3410'      :'422(001)',
        '4818'      :'-42m 001','2730'      :'4mm(001)','8191'      :'4/mmm001','8192'      :'  3(111)',
        '8193'      :' -3(111)','2629888'   :' 32(111)','1319040'   :' 3m(111)','3940737'   :'-3m(111)',
        '32768'     :'  3(+--)','32769'     :' -3(+--)','10519552'  :' 32(+--)','5276160'   :' 3m(+--)',
        '15762945'  :'-3m(+--)','65536'     :'  3(-+-)','65537'     :' -3(-+-)','134808576' :' 32(-+-)',
        '67437056'  :' 3m(-+-)','202180097' :'-3m(-+-)','131072'    :'  3(--+)','131073'    :' -3(--+)',
        '142737664' :' 32(--+)','71434368'  :' 3m(--+)','214040961' :'-3m(--+)','237650'    :'   23   ',
        '237695'    :'   m3   ','715894098' :'   432  ','358068946' :'  -43m  ','1073725439':'   m3m  ',
        '68157504'  :' mm2d100','4456464'   :' mm2d010','642'       :' mm2d001','153092172' :'-4m2 100',
        '277348404' :'-4m2 010','5418'      :'-4m2 001','1075726335':'  6/mmm ','1074414420':'-6m2 100',
        '1075070124':'-6m2 120','1075069650':'   6mm  ','1074414890':'   622  ','1073758215':'   6/m  ',
        '1073758212':'   -6   ','1073758210':'    6   ','1073759865':'-3m(100)','1075724673':'-3m(120)',
        '1073758800':' 3m(100)','1075069056':' 3m(120)','1073759272':' 32(100)','1074413824':' 32(120)',
        '1073758209':'   -3   ','1073758208':'    3   ','1074135143':'mmm(100)','1075314719':'mmm(010)',
        '1073743751':'mmm(110)','1074004034':' mm2z100','1074790418':' mm2z010','1073742466':' mm2z110',
        '1074004004':'mm2(100)','1074790412':'mm2(010)','1073742980':'mm2(110)','1073872964':'mm2(120)',
        '1074266132':'mm2(210)','1073742596':'mm2(+-0)','1073872930':'222(100)','1074266122':'222(010)',
        '1073743106':'222(110)','1073741831':'2/m(001)','1073741921':'2/m(100)','1073741849':'2/m(010)',
        '1073743361':'2/m(110)','1074135041':'2/m(120)','1075314689':'2/m(210)','1073742209':'2/m(+-0)',
        '1073741828':' m(001) ','1073741888':' m(100) ','1073741840':' m(010) ','1073742336':' m(110) ',
        '1074003968':' m(120) ','1074790400':' m(210) ','1073741952':' m(+-0) ','1073741826':' 2(001) ',
        '1073741856':' 2(100) ','1073741832':' 2(010) ','1073742848':' 2(110) ','1073872896':' 2(120) ',
        '1074266112':' 2(210) ','1073742080':' 2(+-0) ','1073741825':'   -1   '
        }
    return KNsym[key]        

def GetNXUPQsym(siteSym):        
    'Needs a doc string'
    NXUPQsym = {
        '    1   ':(28,29,28,28),'   -1   ':( 1,29,28, 0),'    2(x)':(12,18,12,25),'    m(x)':(25,18,12,25),
        '  2/m(x)':( 1,18, 0,-1),'    2(y)':(13,17,13,24),'    m(y)':(24,17,13,24),'  2/m(y)':( 1,17, 0,-1),
        '    2(z)':(14,16,14,23),'    m(z)':(23,16,14,23),'  2/m(z)':( 1,16, 0,-1),'   2(yz)':(10,23,10,22),
        '   m(yz)':(22,23,10,22),' 2/m(yz)':( 1,23, 0,-1),'  2(0+-)':(11,24,11,21),'  m(0+-)':(21,24,11,21),
        '2/m(0+-)':( 1,24, 0,-1),'   2(xz)':( 8,21, 8,20),'   m(xz)':(20,21, 8,20),' 2/m(xz)':( 1,21, 0,-1),
        '  2(+0-)':( 9,22, 9,19),'  m(+0-)':(19,22, 9,19),'2/m(+0-)':( 1,22, 0,-1),'   2(xy)':( 6,19, 6,18),
        '   m(xy)':(18,19, 6,18),' 2/m(xy)':( 1,19, 0,-1),'  2(+-0)':( 7,20, 7,17),'  m(+-0)':(17,20, 7,17),
        '2/m(+-0)':( 1,20, 0,-1),'  mm2(x)':(12,10, 0,-1),'  mm2(y)':(13,10, 0,-1),'  mm2(z)':(14,10, 0,-1),
        ' mm2(yz)':(10,13, 0,-1),'mm2(0+-)':(11,13, 0,-1),' mm2(xz)':( 8,12, 0,-1),'mm2(+0-)':( 9,12, 0,-1),
        ' mm2(xy)':( 6,11, 0,-1),'mm2(+-0)':( 7,11, 0,-1),'  222   ':( 1,10, 0,-1),'  222(x)':( 1,13, 0,-1),
        '  222(y)':( 1,12, 0,-1),'  222(z)':( 1,11, 0,-1),'  mmm   ':( 1,10, 0,-1),'  mmm(x)':( 1,13, 0,-1),
        '  mmm(y)':( 1,12, 0,-1),'  mmm(z)':( 1,11, 0,-1),'  4(100)':(12, 4,12, 0),' -4(100)':( 1, 4,12, 0),
        '4/m(100)':( 1, 4,12,-1),'422(100)':( 1, 4, 0,-1),'-42m 100':( 1, 4, 0,-1),'4mm(100)':(12, 4, 0,-1),
        '4/mmm100':( 1, 4, 0,-1),'  4(010)':(13, 3,13, 0),' -4(010)':( 1, 3,13, 0),'4/m (10)':( 1, 3,13,-1),
        '422(010)':( 1, 3, 0,-1),'-42m 010':( 1, 3, 0,-1),'4mm(010)':(13, 3, 0,-1),'4/mmm010':(1, 3, 0,-1,),
        '  4(001)':(14, 2,14, 0),' -4(001)':( 1, 2,14, 0),'4/m(001)':( 1, 2,14,-1),'422(001)':( 1, 2, 0,-1),
        '-42m 001':( 1, 2, 0,-1),'4mm(001)':(14, 2, 0,-1),'4/mmm001':( 1, 2, 0,-1),'  3(111)':( 2, 5, 2, 0),
        ' -3(111)':( 1, 5, 2, 0),' 32(111)':( 1, 5, 0, 2),' 3m(111)':( 2, 5, 0, 2),'-3m(111)':( 1, 5, 0,-1),
        '  3(+--)':( 5, 8, 5, 0),' -3(+--)':( 1, 8, 5, 0),' 32(+--)':( 1, 8, 0, 5),' 3m(+--)':( 5, 8, 0, 5),
        '-3m(+--)':( 1, 8, 0,-1),'  3(-+-)':( 4, 7, 4, 0),' -3(-+-)':( 1, 7, 4, 0),' 32(-+-)':( 1, 7, 0, 4),
        ' 3m(-+-)':( 4, 7, 0, 4),'-3m(-+-)':( 1, 7, 0,-1),'  3(--+)':( 3, 6, 3, 0),' -3(--+)':( 1, 6, 3, 0),
        ' 32(--+)':( 1, 6, 0, 3),' 3m(--+)':( 3, 6, 0, 3),'-3m(--+)':( 1, 6, 0,-1),'   23   ':( 1, 1, 0, 0),
        '   m3   ':( 1, 1, 0, 0),'   432  ':( 1, 1, 0, 0),'  -43m  ':( 1, 1, 0, 0),'   m3m  ':( 1, 1, 0, 0),
        ' mm2d100':(12,13, 0,-1),' mm2d010':(13,12, 0,-1),' mm2d001':(14,11, 0,-1),'-4m2 100':( 1, 4, 0,-1),
        '-4m2 010':( 1, 3, 0,-1),'-4m2 001':( 1, 2, 0,-1),'  6/mmm ':( 1, 9, 0,-1),'-6m2 100':( 1, 9, 0,-1),
        '-6m2 120':( 1, 9, 0,-1),'   6mm  ':(14, 9, 0,-1),'   622  ':( 1, 9, 0,-1),'   6/m  ':( 1, 9,14,-1),
        '   -6   ':( 1, 9,14, 0),'    6   ':(14, 9,14, 0),'-3m(100)':( 1, 9, 0,-1),'-3m(120)':( 1, 9, 0,-1),
        ' 3m(100)':(14, 9, 0,14),' 3m(120)':(14, 9, 0,14),' 32(100)':( 1, 9, 0,14),' 32(120)':( 1, 9, 0,14),
        '   -3   ':( 1, 9,14, 0),'    3   ':(14, 9,14, 0),'mmm(100)':( 1,14, 0,-1),'mmm(010)':( 1,15, 0,-1),
        'mmm(110)':( 1,11, 0,-1),' mm2z100':(14,14, 0,-1),' mm2z010':(14,15, 0,-1),' mm2z110':(14,11, 0,-1),
        'mm2(100)':(12,14, 0,-1),'mm2(010)':(13,15, 0,-1),'mm2(110)':( 6,11, 0,-1),'mm2(120)':(15,14, 0,-1),
        'mm2(210)':(16,15, 0,-1),'mm2(+-0)':( 7,11, 0,-1),'222(100)':( 1,14, 0,-1),'222(010)':( 1,15, 0,-1),
        '222(110)':( 1,11, 0,-1),'2/m(001)':( 1,16,14,-1),'2/m(100)':( 1,25,12,-1),'2/m(010)':( 1,28,13,-1),
        '2/m(110)':( 1,19, 6,-1),'2/m(120)':( 1,27,15,-1),'2/m(210)':( 1,26,16,-1),'2/m(+-0)':( 1,20,17,-1),
        ' m(001) ':(23,16,14,23),' m(100) ':(26,25,12,26),' m(010) ':(27,28,13,27),' m(110) ':(18,19, 6,18),
        ' m(120) ':(24,27,15,24),' m(210) ':(25,26,16,25),' m(+-0) ':(17,20, 7,17),' 2(001) ':(14,16,14,23),
        ' 2(100) ':(12,25,12,26),' 2(010) ':(13,28,13,27),' 2(110) ':( 6,19, 6,18),' 2(120) ':(15,27,15,24),
        ' 2(210) ':(16,26,16,25),' 2(+-0) ':( 7,20, 7,17),'   -1   ':( 1,29,28, 0)
        }
    return NXUPQsym[siteSym]

def GetCSxinel(siteSym):  
    'Needs a doc string'
    CSxinel = [[],                         # 0th empty - indices are Fortran style
        [[0,0,0],[ 0.0, 0.0, 0.0]],      #1  0  0  0
        [[1,1,1],[ 1.0, 1.0, 1.0]],      #2  X  X  X
        [[1,1,1],[ 1.0, 1.0,-1.0]],      #3  X  X -X
        [[1,1,1],[ 1.0,-1.0, 1.0]],      #4  X -X  X
        [[1,1,1],[ 1.0,-1.0,-1.0]],      #5 -X  X  X
        [[1,1,0],[ 1.0, 1.0, 0.0]],      #6  X  X  0
        [[1,1,0],[ 1.0,-1.0, 0.0]],      #7  X -X  0
        [[1,0,1],[ 1.0, 0.0, 1.0]],      #8  X  0  X
        [[1,0,1],[ 1.0, 0.0,-1.0]],      #9  X  0 -X
        [[0,1,1],[ 0.0, 1.0, 1.0]],      #10  0  Y  Y
        [[0,1,1],[ 0.0, 1.0,-1.0]],      #11 0  Y -Y
        [[1,0,0],[ 1.0, 0.0, 0.0]],      #12  X  0  0
        [[0,1,0],[ 0.0, 1.0, 0.0]],      #13  0  Y  0
        [[0,0,1],[ 0.0, 0.0, 1.0]],      #14  0  0  Z
        [[1,1,0],[ 1.0, 2.0, 0.0]],      #15  X 2X  0
        [[1,1,0],[ 2.0, 1.0, 0.0]],      #16 2X  X  0
        [[1,1,2],[ 1.0, 1.0, 1.0]],      #17  X  X  Z
        [[1,1,2],[ 1.0,-1.0, 1.0]],      #18  X -X  Z
        [[1,2,1],[ 1.0, 1.0, 1.0]],      #19  X  Y  X
        [[1,2,1],[ 1.0, 1.0,-1.0]],      #20  X  Y -X
        [[1,2,2],[ 1.0, 1.0, 1.0]],      #21  X  Y  Y
        [[1,2,2],[ 1.0, 1.0,-1.0]],      #22  X  Y -Y
        [[1,2,0],[ 1.0, 1.0, 0.0]],      #23  X  Y  0
        [[1,0,2],[ 1.0, 0.0, 1.0]],      #24  X  0  Z
        [[0,1,2],[ 0.0, 1.0, 1.0]],      #25  0  Y  Z
        [[1,1,2],[ 1.0, 2.0, 1.0]],      #26  X 2X  Z
        [[1,1,2],[ 2.0, 1.0, 1.0]],      #27 2X  X  Z
        [[1,2,3],[ 1.0, 1.0, 1.0]],      #28  X  Y  Z
        ]
    indx = GetNXUPQsym(siteSym)
    return CSxinel[indx[0]]
    
def GetCSuinel(siteSym):
    "returns Uij terms, multipliers, GUI flags & Uiso2Uij multipliers"
    CSuinel = [[],                                             # 0th empty - indices are Fortran style
        [[1,1,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,0,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #1  A  A  A  0  0  0
        [[1,1,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,0,1,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #2  A  A  C  0  0  0
        [[1,2,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #3  A  B  A  0  0  0
        [[1,2,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #4  A  B  B  0  0  0
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #5  A  A  A  D  D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #6  A  A  A  D -D -D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #7  A  A  A  D -D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #8  A  A  A  D  D -D
        [[1,1,2,1,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,0,1,0,0,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #9  A  A  C A/2 0  0
        [[1,2,3,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #10  A  B  C  0  0  0
        [[1,1,2,3,0,0],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0],[1,0,1,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #11  A  A  C  D  0  0
        [[1,2,1,0,3,0],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],[1,1,0,0,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #12  A  B  A  0  E  0
        [[1,2,2,0,0,3],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],[1,1,0,0,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #13  A  B  B  0  0  F
        [[1,2,3,2,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.5,0.0]],    #14  A  B  C B/2 0  0
        [[1,2,3,1,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.5,0.0]],    #15  A  B  C A/2 0  0
        [[1,2,3,4,0,0],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0],[1,1,1,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #16  A  B  C  D  0  0
        [[1,2,3,0,4,0],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],[1,1,1,0,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #17  A  B  C  0  E  0
        [[1,2,3,0,0,4],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #18  A  B  C  0  0  F
        [[1,1,2,3,4,4],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,0,1,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #19  A  A  C  D  E -E
        [[1,1,2,3,4,4],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,0,1,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #20  A  A  C  D  E  E
        [[1,2,1,3,4,3],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,1,0,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #21  A  B  A  D  E -D
        [[1,2,1,3,4,3],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,0,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #22  A  B  A  D  E  D
        [[1,2,2,3,3,4],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0],[1,1,0,1,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #23  A  B  B  D -D  F
        [[1,2,2,3,3,4],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,0,1,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #24  A  B  B  D  D  F
        [[1,2,3,2,4,4],[ 1.0, 1.0, 1.0, 0.5, 0.5, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #25  A  B  C B/2 F/2 F
        [[1,2,3,1,0,4],[ 1.0, 1.0, 1.0, 0.5, 0.0, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #26  A  B  C A/2  0  F
        [[1,2,3,2,4,0],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.0],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #27  A  B  C B/2  E  0
        [[1,2,3,1,4,4],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.5],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #28  A  B  C A/2  E E/2
        [[1,2,3,4,5,6],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,1,1,1,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #29  A  B  C  D  E   F
        ]
    indx = GetNXUPQsym(siteSym)
    return CSuinel[indx[1]]
    
def MustrainNames(SGData):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    if laue in ['m3','m3m']:
        return ['S400','S220']
    elif laue in ['6/m','6/mmm','3m1']:
        return ['S400','S004','S202']
    elif laue in ['31m','3']:
        return ['S400','S004','S202','S211']
    elif laue in ['3R','3mR']:
        return ['S400','S220','S310','S211']
    elif laue in ['4/m','4/mmm']:
        return ['S400','S004','S220','S022']
    elif laue in ['mmm']:
        return ['S400','S040','S004','S220','S202','S022']
    elif laue in ['2/m']:
        SHKL = ['S400','S040','S004','S220','S202','S022']
        if uniq == 'a':
            SHKL += ['S013','S031','S211']
        elif uniq == 'b':
            SHKL += ['S301','S103','S121']
        elif uniq == 'c':
            SHKL += ['S130','S310','S112']
        return SHKL
    else:
        SHKL = ['S400','S040','S004','S220','S202','S022']
        SHKL += ['S310','S103','S031','S130','S301','S013']
        SHKL += ['S211','S121','S112']
        return SHKL

def HStrainNames(SGData):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    if laue in ['m3','m3m']:
        return ['D11','eA']         #add cubic strain term
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        return ['D11','D33']
    elif laue in ['3R','3mR']:
        return ['D11','D12']
    elif laue in ['4/m','4/mmm']:
        return ['D11','D33']
    elif laue in ['mmm']:
        return ['D11','D22','D33']
    elif laue in ['2/m']:
        Dij = ['D11','D22','D33']
        if uniq == 'a':
            Dij += ['D23']
        elif uniq == 'b':
            Dij += ['D13']
        elif uniq == 'c':
            Dij += ['D12']
        return Dij
    else:
        Dij = ['D11','D22','D33','D12','D13','D23']
        return Dij
    
def MustrainCoeff(HKL,SGData):
    'Needs a doc string'
    #NB: order of terms is the same as returned by MustrainNames
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = HKL
    Strm = []
    if laue in ['m3','m3m']:
        Strm.append(h**4+k**4+l**4)
        Strm.append(3.0*((h*k)**2+(h*l)**2+(k*l)**2))
    elif laue in ['6/m','6/mmm','3m1']:
        Strm.append(h**4+k**4+2.0*k*h**3+2.0*h*k**3+3.0*(h*k)**2)
        Strm.append(l**4)
        Strm.append(3.0*((h*l)**2+(k*l)**2+h*k*l**2))
    elif laue in ['31m','3']:
        Strm.append(h**4+k**4+2.0*k*h**3+2.0*h*k**3+3.0*(h*k)**2)
        Strm.append(l**4)
        Strm.append(3.0*((h*l)**2+(k*l)**2+h*k*l**2))
        Strm.append(4.0*h*k*l*(h+k))
    elif laue in ['3R','3mR']:
        Strm.append(h**4+k**4+l**4)
        Strm.append(3.0*((h*k)**2+(h*l)**2+(k*l)**2))
        Strm.append(2.0*(h*l**3+l*k**3+k*h**3)+2.0*(l*h**3+k*l**3+l*k**3))
        Strm.append(4.0*(k*l*h**2+h*l*k**2+h*k*l**2))
    elif laue in ['4/m','4/mmm']:
        Strm.append(h**4+k**4)
        Strm.append(l**4)
        Strm.append(3.0*(h*k)**2)
        Strm.append(3.0*((h*l)**2+(k*l)**2))
    elif laue in ['mmm']:
        Strm.append(h**4)
        Strm.append(k**4)
        Strm.append(l**4)
        Strm.append(3.0*(h*k)**2)
        Strm.append(3.0*(h*l)**2)
        Strm.append(3.0*(k*l)**2)
    elif laue in ['2/m']:
        Strm.append(h**4)
        Strm.append(k**4)
        Strm.append(l**4)
        Strm.append(3.0*(h*k)**2)
        Strm.append(3.0*(h*l)**2)
        Strm.append(3.0*(k*l)**2)
        if uniq == 'a':
            Strm.append(2.0*k*l**3)
            Strm.append(2.0*l*k**3)
            Strm.append(4.0*k*l*h**2)
        elif uniq == 'b':
            Strm.append(2.0*l*h**3)
            Strm.append(2.0*h*l**3)
            Strm.append(4.0*h*l*k**2)
        elif uniq == 'c':
            Strm.append(2.0*h*k**3)
            Strm.append(2.0*k*h**3)
            Strm.append(4.0*h*k*l**2)
    else:
        Strm.append(h**4)
        Strm.append(k**4)
        Strm.append(l**4)
        Strm.append(3.0*(h*k)**2)
        Strm.append(3.0*(h*l)**2)
        Strm.append(3.0*(k*l)**2)
        Strm.append(2.0*k*h**3)
        Strm.append(2.0*h*l**3)
        Strm.append(2.0*l*k**3)
        Strm.append(2.0*h*k**3)
        Strm.append(2.0*l*h**3)
        Strm.append(2.0*k*l**3)
        Strm.append(4.0*k*l*h**2)
        Strm.append(4.0*h*l*k**2)
        Strm.append(4.0*k*h*l**2)
    return Strm
    
def Muiso2Shkl(muiso,SGData,cell):
    "this is to convert isotropic mustrain to generalized Shkls"
    import GSASIIlattice as G2lat
    A = G2lat.cell2AB(cell)[0]
    
    def minMus(Shkl,muiso,H,SGData,A):
        U = np.inner(A.T,H)
        S = np.array(MustrainCoeff(U,SGData))
        Sum = np.sqrt(np.sum(np.multiply(S,Shkl[:,np.newaxis]),axis=0))
        rad = np.sqrt(np.sum((Sum[:,np.newaxis]*H)**2,axis=1))
        return (muiso-rad)**2
        
    laue = SGData['SGLaue']
    PHI = np.linspace(0.,360.,60,True)
    PSI = np.linspace(0.,180.,60,True)
    X = np.outer(npsind(PHI),npsind(PSI))
    Y = np.outer(npcosd(PHI),npsind(PSI))
    Z = np.outer(np.ones(np.size(PHI)),npcosd(PSI))
    HKL = np.dstack((X,Y,Z))
    if laue in ['m3','m3m']:
        S0 = [1000.,1000.]
    elif laue in ['6/m','6/mmm','3m1']:
        S0 = [1000.,1000.,1000.]
    elif laue in ['31m','3']:
        S0 = [1000.,1000.,1000.,1000.]
    elif laue in ['3R','3mR']:
        S0 = [1000.,1000.,1000.,1000.]
    elif laue in ['4/m','4/mmm']:
        S0 = [1000.,1000.,1000.,1000.]
    elif laue in ['mmm']:
        S0 = [1000.,1000.,1000.,1000.,1000.,1000.]
    elif laue in ['2/m']:
        S0 = [1000.,1000.,1000.,0.,0.,0.,0.,0.,0.]
    else:
        S0 = [1000.,1000.,1000.,1000.,1000., 1000.,1000.,1000.,1000.,1000., 
            1000.,1000.,0.,0.,0.]
    S0 = np.array(S0)
    HKL = np.reshape(HKL,(-1,3))
    result = so.leastsq(minMus,S0,(np.ones(HKL.shape[0])*muiso,HKL,SGData,A))
    return result[0]
       
def SytSym(XYZ,SGData):
    '''
    Generates the number of equivalent positions and a site symmetry code for a specified coordinate and space group

    :param XYZ: an array, tuple or list containing 3 elements: x, y & z
    :param SGData: from SpcGroup
    :Returns: a two element tuple:

     * The 1st element is a code for the site symmetry (see GetKNsym)
     * The 2nd element is the site multiplicity

    '''
    def PackRot(SGOps):
        IRT = []
        for ops in SGOps:
            M = ops[0]
            irt = 0
            for j in range(2,-1,-1):
                for k in range(2,-1,-1):
                    irt *= 3
                    irt += M[k][j]
            IRT.append(int(irt))
        return IRT 
        
    SymName = ''
    Mult = 1
    Isym = 0
    if SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
        Isym = 1073741824
    Jdup = 0
    Xeqv = GenAtom(XYZ,SGData,True)
    IRT = PackRot(SGData['SGOps'])
    L = -1
    for ic,cen in enumerate(SGData['SGCen']):
        for invers in range(int(SGData['SGInv']+1)):
            for io,ops in enumerate(SGData['SGOps']):
                irtx = (1-2*invers)*IRT[io]
                L += 1
                if not Xeqv[L][1]:
                    Jdup += 1
                    jx = GetOprPtrName(str(irtx))
                    if jx[2] < 39:
                        Isym += 2**(jx[2]-1)
    if Isym == 1073741824: Isym = 0
    Mult = len(SGData['SGOps'])*len(SGData['SGCen'])*(int(SGData['SGInv'])+1)/Jdup
          
    return GetKNsym(str(Isym)),Mult
    
def ElemPosition(SGData):
    ''' Under development. 
    Object here is to return a list of symmetry element types and locations suitable
    for say drawing them.
    So far I have the element type... getting all possible locations without lookup may be impossible!
    '''
    SymElements = []
    Inv = SGData['SGInv']
    Cen = SGData['SGCen']
    eleSym = {-3:['','-1'],-2:['',-6],-1:['2','-4'],0:['3','-3'],1:['4','m'],2:['6',''],3:['1','']}
    # get operators & expand if centrosymmetric
    Ops = SGData['SGOps']
    opM = np.array([op[0].T for op in Ops])
    opT = np.array([op[1] for op in Ops])
    if Inv:
        opM = np.concatenate((opM,-opM))
        opT = np.concatenate((opT,-opT))
    opMT = zip(opM,opT)
    for M,T in opMT[1:]:        #skip I
        Dt = int(nl.det(M))
        Tr = int(np.trace(M))
        Dt = -(Dt-1)/2
        Es = eleSym[Tr][Dt]
        if Dt:              #rotation-inversion
            I = np.eye(3)
            if Tr == 1:     #mirrors/glides
                if np.any(T):       #glide
                    M2 = np.inner(M,M)
                    MT = np.inner(M,T)+T
                    print 'glide',Es,MT
                    print M2
                else:               #mirror
                    print 'mirror',Es,T
                    print I-M
                X = [-1,-1,-1]
            elif Tr == -3:  # pure inversion
                X = np.inner(nl.inv(I-M),T)
                print 'inversion',Es,X
            else:           #other rotation-inversion
                M2 = np.inner(M,M)
                MT = np.inner(M,T)+T
                print 'rot-inv',Es,MT
                print M2
                X = [-1,-1,-1]
        else:               #rotations
            print 'rotation',Es
            X = [-1,-1,-1]
        #SymElements.append([Es,X])
        
    return #SymElements
    
def ApplyStringOps(A,SGData,X,Uij=[]):
    'Needs a doc string'
    SGOps = SGData['SGOps']
    SGCen = SGData['SGCen']
    Ax = A.split('+')
    Ax[0] = int(Ax[0])
    iC = 0
    if Ax[0] < 0:
        iC = 1
    Ax[0] = abs(Ax[0])
    nA = Ax[0]%100-1
    cA = Ax[0]/100
    Cen = SGCen[cA]
    M,T = SGOps[nA]
    if len(Ax)>1:
        cellA = Ax[1].split(',')
        cellA = np.array([int(a) for a in cellA])
    else:
        cellA = np.zeros(3)
    newX = (1-2*iC)*(Cen+np.inner(M,X)+T)+cellA
    if len(Uij):
        U = Uij2U(Uij)
        U = np.inner(M,np.inner(U,M).T)
        newUij = U2Uij(U)
        return [newX,newUij]
    else:
        return newX
        
def StringOpsProd(A,B,SGData):
    """
    Find A*B where A & B are in strings '-' + '100*c+n' + '+ijk'
    where '-' indicates inversion, c(>0) is the cell centering operator, 
    n is operator number from SgOps and ijk are unit cell translations (each may be <0).
    Should return resultant string - C. SGData - dictionary using entries:

       *  'SGCen': cell centering vectors [0,0,0] at least
       *  'SGOps': symmetry operations as [M,T] so that M*x+T = x'

    """
    SGOps = SGData['SGOps']
    SGCen = SGData['SGCen']
    #1st split out the cell translation part & work on the operator parts
    Ax = A.split('+'); Bx = B.split('+')
    Ax[0] = int(Ax[0]); Bx[0] = int(Bx[0])
    iC = 0
    if Ax[0]*Bx[0] < 0:
        iC = 1
    Ax[0] = abs(Ax[0]); Bx[0] = abs(Bx[0])
    nA = Ax[0]%100-1;  nB = Bx[0]%100-1
    cA = Ax[0]/100;  cB = Bx[0]/100
    Cen = (SGCen[cA]+SGCen[cB])%1.0
    cC = np.nonzero([np.allclose(C,Cen) for C in SGCen])[0][0]
    Ma,Ta = SGOps[nA]; Mb,Tb = SGOps[nB]
    Mc = np.inner(Ma,Mb.T)
#    print Ma,Mb,Mc
    Tc = (np.add(np.inner(Mb,Ta)+1.,Tb))%1.0
#    print Ta,Tb,Tc
#    print [np.allclose(M,Mc)&np.allclose(T,Tc) for M,T in SGOps]
    nC = np.nonzero([np.allclose(M,Mc)&np.allclose(T,Tc) for M,T in SGOps])[0][0]
    #now the cell translation part
    if len(Ax)>1:
        cellA = Ax[1].split(',')
        cellA = [int(a) for a in cellA]
    else:
        cellA = [0,0,0]
    if len(Bx)>1:
        cellB = Bx[1].split(',')
        cellB = [int(b) for b in cellB]
    else:
        cellB = [0,0,0]
    cellC = np.add(cellA,cellB)
    C = str(((nC+1)+(100*cC))*(1-2*iC))+'+'+ \
        str(int(cellC[0]))+','+str(int(cellC[1]))+','+str(int(cellC[2]))
    return C
            
def U2Uij(U):
    #returns the UIJ vector U11,U22,U33,U12,U13,U23 from tensor U
    return [U[0][0],U[1][1],U[2][2],2.*U[0][1],2.*U[0][2],2.*U[1][2]]
    
def Uij2U(Uij):
    #returns the thermal motion tensor U from Uij as numpy array
    return np.array([[Uij[0],Uij[3]/2.,Uij[4]/2.],[Uij[3]/2.,Uij[1],Uij[5]/2.],[Uij[4]/2.,Uij[5]/2.,Uij[2]]])

def StandardizeSpcName(spcgroup):
    '''Accept a spacegroup name where spaces may have not been used
    in the names according to the GSAS convention (spaces between symmetry
    for each axis) and return the space group name as used in GSAS
    '''
    rspc = spcgroup.replace(' ','').upper()
    # deal with rhombohedral and hexagonal setting designations
    rhomb = ''
    if rspc[-1:] == 'R':
        rspc = rspc[:-1]
        rhomb = ' R'
    if rspc[-1:] == 'H': # hexagonal is assumed and thus can be ignored
        rspc = rspc[:-1]
    # look for a match in the spacegroup lists
    for i in spglist.values():
        for spc in i:
            if rspc == spc.replace(' ','').upper():
                return spc + rhomb
    # how about the post-2002 orthorhombic names?
    for i,spc in sgequiv_2002_orthorhombic:
        if rspc == i.replace(' ','').upper():
            return spc
    # not found
    return ''

    
spglist = {}
'''A dictionary of space groups as ordered and named in the pre-2002 International 
Tables Volume A, except that spaces are used following the GSAS convention to 
separate the different crystallographic directions.
Note that the symmetry codes here will recognize many non-standard space group 
symbols with different settings. They are ordered by Laue group
'''
spglist = {
    'P1' : ('P 1','P -1',), # 1-2
    'P2/m': ('P 2','P 21','P m','P a','P c','P n',
        'P 2/m','P 21/m','P 2/c','P 2/a','P 2/n','P 21/c','P 21/a','P 21/n',), #3-15
    'C2/m':('C 2','C m','C c','C n',
        'C 2/m','C 2/c','C 2/n',),
    'Pmmm':('P 2 2 2',
        'P 2 2 21','P 2 21 2','P 21 2 2',
        'P 21 21 2','P 21 2 21','P 2 21 21',
        'P 21 21 21',
        'P m m 2','P m 2 m','P 2 m m',
        'P m c 21','P c m 21','P 21 m a','P 21 a m','P b 21 m','P m 21 b',
        'P c c 2','P 2 a a','P b 2 b',
        'P m a 2','P b m 2','P 2 m b','P 2 c m','P c 2 m','P m 2 a',
        'P c a 21','P b c 21','P 21 a b','P 21 c a','P c 21 b','P b 21 a',
        'P n c 2','P c n 2','P 2 n a','P 2 a n','P b 2 n','P n 2 b',
        'P m n 21','P n m 21','P 21 m n','P 21 n m','P n 21 m','P m 21 n',
        'P b a 2','P 2 c b','P c 2 a',
        'P n a 21','P b n 21','P 21 n b','P 21 c n','P c 21 n','P n 21 a',
        'P n n 2','P 2 n n','P n 2 n',
        'P m m m','P n n n',
        'P c c m','P m a a','P b m b',
        'P b a n','P n c b','P c n a',
        'P m m a','P m m b','P b m m','P c m m','P m c m','P m a m',
        'P n n a','P n n b','P b n n','P c n n','P n c n','P n a n',
        'P m n a','P n m b','P b m n','P c n m','P n c m','P m a n',
        'P c c a','P c c b','P b a a','P c a a','P b c b','P b a b',
        'P b a m','P m c b','P c m a',
        'P c c n','P n a a','P b n b',
        'P b c m','P c a m','P m c a','P m a b','P b m a','P c m b',
        'P n n m','P m n n','P n m n',
        'P m m n','P n m m','P m n m',
        'P b c n','P c a n','P n c a','P n a b','P b n a','P c n b',
        'P b c a','P c a b',
        'P n m a','P m n b','P b n m','P c m n','P m c n','P n a m',
        ),
    'Cmmm':('C 2 2 21','C 2 2 2','C m m 2','C m c 21','C c c 2','C m 2 m','C 2 m m',
        'C m 2 a','C 2 m b','C 2 c m','C c 2 m','C 2 c m',
        'C m c a','C m m m','C c c m','C m m a','C c c a','C m c m',),
    'Immm':('I 2 2 2','I 21 21 21','I m m m',
        'I m m 2','I m 2 m','I 2 m m',
        'I b a 2','I 2 c b','I c 2 a',
        'I m a 2','I b m 2','I 2 m b','I 2 c m','I c 2 m','I m 2 a',
        'I b a m','I m c b','I c m a',
        'I b c a','I c a b',
        'I m m a','I m m b','I b m m ','I c m m','I m c m','I m a m',),
    'Fmmm':('F 2 2 2','F m m m', 'F d d d',
        'F m m 2','F m 2 m','F 2 m m',
        'F d d 2','F d 2 d','F 2 d d',),
    'P4/mmm':('P 4','P 41','P 42','P 43','P -4','P 4/m','P 42/m','P 4/n','P 42/n',
        'P 4 2 2','P 4 21 2','P 41 2 2','P 41 21 2','P 42 2 2',
        'P 42 21 2','P 43 2 2','P 43 21 2','P 4 m m','P 4 b m','P 42 c m',
        'P 42 n m','P 4 c c','P 4 n c','P 42 m c','P 42 b c','P -4 2 m',
        'P -4 2 c','P -4 21 m','P -4 21 c','P -4 m 2','P -4 c 2','P -4 b 2',
        'P -4 n 2','P 4/m m m','P 4/m c c','P 4/n b m','P 4/n n c','P 4/m b m',
        'P 4/m n c','P 4/n m m','P 4/n c c','P 42/m m c','P 42/m c m',
        'P 42/n b c','P 42/n n m','P 42/m b c','P 42/m n m','P 42/n m c',
        'P 42/n c m',),
    'I4/mmm':('I 4','I 41','I -4','I 4/m','I 41/a','I 4 2 2','I 41 2 2','I 4 m m',
        'I 4 c m','I 41 m d','I 41 c d',
        'I -4 m 2','I -4 c 2','I -4 2 m','I -4 2 d','I 4/m m m','I 4/m c m',
        'I 41/a m d','I 41/a c d'),
    'R3-H':('R 3','R -3','R 3 2','R 3 m','R 3 c','R -3 m','R -3 c',),
    'P6/mmm': ('P 3','P 31','P 32','P -3','P 3 1 2','P 3 2 1','P 31 1 2',
        'P 31 2 1','P 32 1 2','P 32 2 1', 'P 3 m 1','P 3 1 m','P 3 c 1',
        'P 3 1 c','P -3 1 m','P -3 1 c','P -3 m 1','P -3 c 1','P 6','P 61',
        'P 65','P 62','P 64','P 63','P -6','P 6/m','P 63/m','P 6 2 2',
        'P 61 2 2','P 65 2 2','P 62 2 2','P 64 2 2','P 63 2 2','P 6 m m',
        'P 6 c c','P 63 c m','P 63 m c','P -6 m 2','P -6 c 2','P -6 2 m',
        'P -6 2 c','P 6/m m m','P 6/m c c','P 63/m c m','P 63/m m c',),
    'Pm3m': ('P 2 3','P 21 3','P m 3','P n 3','P a 3','P 4 3 2','P 42 3 2',
        'P 43 3 2','P 41 3 2','P -4 3 m','P -4 3 n','P m 3 m','P n 3 n',
        'P m 3 n','P n 3 m',),
    'Im3m':('I 2 3','I 21 3','I m -3','I a -3', 'I 4 3 2','I 41 3 2',
        'I -4 3 m', 'I -4 3 d','I m -3 m','I m 3 m','I a -3 d',),
    'Fm3m':('F 2 3','F m -3','F d -3','F 4 3 2','F 41 3 2','F -4 3 m',
        'F -4 3 c','F m -3 m','F m 3 m','F m -3 c','F d -3 m','F d -3 c',),
}

ssdict = {}
'''A dictionary of superspace group symbols allowed for each entry in spglist
(except cubics). Monoclinics are all b-unique setting.
'''
ssdict = {
    'P 1':['(abg)',],'P -1':['(abg)',],
    #monoclinic - done
    'P 2':['(a0g)','(a1/2g)','(0b0)','(0b0)s','(1/2b0)','(0b1/2)',],
    'P 21':['(a0g)','(0b0)','(0b0)s','(1/2b0)','(0b1/2)','(1/2b0)s','(0b1/2)s',],
    'P m':['(a0g)','(a0g)s','(a1/2g)','(a1/2g)s','(0b0)','(1/2b0)','(0b1/2)',],
    'P a':['(a0g)','(a1/2g)','(a0g)s','(a1/2g)s','(0b0)','(0b1/2)',],
    'P c':['(a0g)','(a1/2g)','(a0g)s','(a1/2g)s','(0b0)','(1/2b0)',],
    'P n':['(a0g)','(a1/2g)','(a0g)s','(a1/2g)s','(0b0)','(1/2b1/2)',],
    'P 2/m':['(a0g)','(a1/2g)','(a0g)0s','(a1/2g)0s',
        '(0b0)','(0b0)s0','(1/2b0)','(0b1/2)','(1/2b0)s0','(0b1/2)s0',],
    'P 21/m':['(a0g)','(a0g)0s','(0b0)','(0b0)s',
        '(1/2b0)','(0b1/2)','(1/2b0)s0','(0b1/2)s0'],
    'P 2/c':['(a0g)','(a1/2g)','(a0g)0s','(a1/2g)0s',
        '(0b0)','(0b0)s0','(1/2b0)','(1/2b0)s0',],
    'P 2/a':['(a0g)','(a1/2g)','(a0g)0s','(a1/2g)0s',
        '(0b0)','(0b0)s0','(0b1/2)','(0b1/2)s0',],
    'P 2/n':['(a0g)','(a1/2g)','(a0g)0s','(a1/2g)0s',
        '(0b0)','(0b0)s0','(1/2b1/2)','(1/2b1/2)s0',],
    'P 21/c':['(a0g)','(a0g)0s','(0b0)','(0b0)s0','(1/2b0)','(1/2b0)s0',],
    'P 21/a':['(a0g)','(a0g)0s','(0b0)','(0b0)s0','(0b1/2)','(0b1/2)s0',],
    'P 21/n':['(a0g)','(a0g)0s','(0b0)','(0b0)s0','(1/2b1/2)','(1/2b1/2)s0',],
    'C 2':['(a0g)','(0b0)','(0b0)s','(0b1/2)','(0b1/2)s',],
    'C m':['(a0g)','(a0g)s','(0b0)','(0b1/2)',],
    'C c':['(a0g)','(a0g)s','(0b0)',],
    'C n':['(a0g)','(a0g)s','(0b0)',],
    'C 2/m':['(a0g)','(a0g)0s','(0b0)','(0b0)s0',
        '(1/2b0)','(0b1/2)','(1/2b0)s0','(0b1/2)s0',],
    'C 2/c':['(a0g)','(a0g)0s','(0b0)','(0b0)s0',],
    'C 2/n':['(a0g)','(a0g)0s','(0b0)','(0b0)s0',],
    #orthorhombic    
    'P 2 2 2':['(a00)','(0b0)','(00g)','(a00)s00','(0b0)0s0','(00g)00s',
        '(a1/20)','(a01/2)','(0b1/2)','(1/2b0)','(01/2g)','(1/20g)',
        '(a1/21/2)','(1/2b1/2)','(1/21/2g)',],
        
    'P 2 2 21':['(a00)','(0b0)','(00g)','(a00)s00','(0b0)0s0',
        '(a1/20)','(1/2b0)','(01/2g)','(1/20g)','(1/21/2g)',],
    'P 2 21 2':['(a00)','(0b0)','(00g)','(a00)s00','(00g)00s',
        '(a01/2)','(1/20g)','(0b1/2)','(1/2b0)','(1/2b1/2)',],
    'P 21 2 2':['(a00)','(0b0)','(00g)','(0b0)0s0','(00g)00s',
        '(01/2g)','(0b1/2)','(a01/2)','(a1/20)','(a1/21/2)',],
        
    'P 21 21 2':['(a00)','(0b0)','(00g)','(00g)00s','(a01/2)','(0b1/2)',],
    'P 21 2 21':['(a00)','(0b0)','(00g)','(0b0)0s0','(a1/20)','(01/2g)',],
    'P 2 21 21':['(a00)','(0b0)','(00g)','(a00)s00','(1/2b0)','(1/20g)',],
        
    'P 21 21 21':['(a00)','(0b0)','(00g)',],
        
    'P m m 2':['(a00)','(0b0)','(00g)',
        '(a00)s00','(0b0)0s0','(00g)s0s','(00g)ss0','(00g)0ss',
        '(a1/20)','(a01/2)','(0b1/2)','(1/2b0)','(01/2g)','(1/20g)',
        '(a1/21/2)','(1/2b1/2)','(1/21/2g)',],
    'P m 2 m':[],
    'P 2 m m':[],
        
    'P m c 21':[],
    'P c m 21':[],
    'P 21 m a':[],
    'P 21 a m':[],
    'P b 21 m':[],
    'P m 21 b':[],
        
    'P c c 2':[],
    'P 2 a a':[],
    'P b 2 b':[],
        
    'P m a 2':[],
    'P b m 2':[],
    'P 2 m b':[],
    'P 2 c m':[],
    'P c 2 m':[],
    'P m 2 a':[],
        
    'P c a 21':[],
    'P b c 21':[],
    'P 21 a b':[],
    'P 21 c a':[],
    'P c 21 b':[],
    'P b 21 a':[],
        
    'P n c 2':[],
    'P c n 2':[],
    'P 2 n a':[],
    'P 2 a n':[],
    'P b 2 n':[],
    'P n 2 b':[],
        
    'P m n 21':[],
    'P n m 21':[],
    'P 21 m n':[],
    'P 21 n m':[],
    'P n 21 m':[],
    'P m 21 n':[],
        
    'P b a 2':[],
    'P 2 c b':[],
    'P c 2 a':[],
        
    'P n a 21':[],
    'P b n 21':[],
    'P 21 n b':[],
    'P 21 c n':[],
    'P c 21 n':[],
    'P n 21 a':[],
        
    'P n n 2':[],
    'P 2 n n':[],
    'P n 2 n':[],
        
    'P m m m':[],
        
    'P n n n':[],
        
    'P c c m':[],
    'P m a a':[],
    'P b m b':[],
        
    'P b a n':[],
    'P n c b':[],
    'P c n a':[],
        
    'P m m a':[],
    'P m m b':[],
    'P b m m':[],
    'P c m m':[],
    'P m c m':[],
    'P m a m':[],
        
    'P n n a':[],
    'P n n b':[],
    'P b n n':[],
    'P c n n':[],
    'P n c n':[],
    'P n a n':[],
        
    'P m n a':[],
    'P n m b':[],
    'P b m n':[],
    'P c n m':[],
    'P n c m':[],
    'P m a n':[],
        
    'P c c a':[],
    'P c c b':[],
    'P b a a':[],
    'P c a a':[],
    'P b c b':[],
    'P b a b':[],
        
    'P b a m':[],
    'P m c b':[],
    'P c m a':[],
        
    'P c c n':[],
    'P n a a':[],
    'P b n b':[],
        
    'P b c m':[],
    'P c a m':[],
    'P m c a':[],
    'P m a b':[],
    'P b m a':[],
    'P c m b':[],
        
    'P n n m':[],
    'P m n n':[],
    'P n m n':[],
        
    'P m m n':[],
    'P n m m':[],
    'P m n m':[],
        
    'P b c n':[],
    'P c a n':[],
    'P n c a':[],
    'P n a b':[],
    'P b n a':[],
    'P c n b':[],
        
    'P b c a':[],
    'P c a b':[],
        
    'P n m a':[],
    'P m n b':[],
    'P b n m':[],
    'P c m n':[],
    'P m c n':[],
    'P n a m':[],
        
    'C 2 2 21':['(a00)','(0b0)','(00g)','(10g)','(01g)',],
    'C 2 2 2':[],
    'C m m 2':[],
    'C m c 21':[],
    'C c c 2':[],
        
    'C m 2 m':[],
    'C 2 m m':[],
        
    'C m 2 a':[],
    'C 2 m b':[],
        
    'C 2 c m':[],
    'C c 2 m':[],
    'C 2 c m':[],
        
    'C m c a':[],
    'C m m m':[],
    'C c c m':[],
    'C m m a':[],
    'C c c a':[],
    'C m c m':[],
    'I 2 2 2':[],
    'I 21 21 21':[],
    'I m m m':[],
        
    'I m m 2':[],
    'I m 2 m':[],
    'I 2 m m':[],
        
    'I b a 2':[],
    'I 2 c b':[],
    'I c 2 a':[],
        
    'I m a 2':[],
    'I b m 2':[],
    'I 2 m b':[],
        
    'I 2 c m':[],
    'I c 2 m':[],
    'I m 2 a':[],
        
    'I b a m':[],
    'I m c b':[],
    'I c m a':[],
        
    'I b c a':[],
    'I c a b':[],
        
    'I m m a':[],
    'I m m b':[],
    'I b m m ':[],
    'I c m m':[],
    'I m c m':[],
    'I m a m':[],
        
    'F 2 2 2':[],
    'F m m m':[],
    'F d d d':[],
        
    'F m m 2':[],
    'F m 2 m':[],
    'F 2 m m':[],
        
    'F d d 2':[],
    'F d 2 d':[],
    'F 2 d d':[],        
    #tetragonal - done
    'P 4':['(00g)','(00g)q','(00g)s','(1/21/2g)','(1/21/2g)q',],
    'P 41':['(00g)','(00g)q','(00g)s','(1/21/2g)','(1/21/2g)q',],
    'P 42':['(00g)','(00g)q','(00g)s','(1/21/2g)','(1/21/2g)q',],
    'P 43':['(00g)','(00g)q','(00g)s','(1/21/2g)','(1/21/2g)q',],
    'P -4':['(00g)','(1/21/2g)',],
    'P 4/m':['(00g)','(00g)s0','(1/21/2g)','(1/21/2g)s0',],
    'P 42/m':['(00g)','(00g)s0','(1/21/2g)','(1/21/2g)s0',],
    'P 4/n':['(00g)','(00g)s0','(1/21/2g)','(1/21/2g)q0',],
    'P 42/n':['(00g)','(00g)s0','(1/21/2g)','(1/21/2g)q0',],
    'P 4 2 2':['(00g)','(00g)q00','(00g)s00','(1/21/2g)','(1/21/2g)q00','(1/21/2g)s00',],
    'P 4 21 2':['(00g)','(00g)q00','(00g)s00',],
    'P 41 2 2':['(00g)','(00g)q00','(00g)s00','(1/21/2g)','(1/21/2g)q00','(1/21/2g)s00',],
    'P 41 21 2':['(00g)','(00g)q00','(00g)s00',],
    'P 42 2 2':['(00g)','(00g)q00','(00g)s00','(1/21/2g)','(1/21/2g)q00','(1/21/2g)s00',],
    'P 42 21 2':['(00g)','(00g)q00','(00g)s00',],
    'P 43 2 2':['(00g)','(00g)q00','(00g)s00','(1/21/2g)','(1/21/2g)q00','(1/21/2g)s00',],
    'P 43 21 2':['(00g)','(00g)q00','(00g)s00',],
    'P 4 m m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',
        '(1/21/2g)','(1/21/2g)ss0','(1/21/2g)0ss','(1/21/2g)s0s',],
    'P 4 b m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s','(1/21/2g)','(1/21/2g)qq0','(1/21/2g)qqs',],
    'P 42 c m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',
        '(1/21/2g)','(1/21/2g)ss0','(1/21/2g)0ss','(1/21/2g)s0s',],
    'P 42 n m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s','(1/21/2g)','(1/21/2g)qq0','(1/21/2g)qqs',],
    'P 4 c c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',
        '(1/21/2g)','(1/21/2g)ss0','(1/21/2g)0ss','(1/21/2g)s0s',],
    'P 4 n c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s','(1/21/2g)','(1/21/2g)qq0','(1/21/2g)qqs'],
    'P 42 m c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',
        '(1/21/2g)','(1/21/2g)ss0','(1/21/2g)0ss','(1/21/2g)s0s',],
    'P 42 b c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s','(1/21/2g)','(1/21/2g)qq0','(1/21/2g)qqs'],
    'P -4 2 m':['(00g)','(00g)00s','(1/21/2g)','(1/21/2g)00s',],
    'P -4 2 c':['(00g)','(00g)00s',],
    'P -4 21 m':['(00g)','(00g)00s',],
    'P -4 21 c':['(00g)','(00g)00s',],
    'P -4 m 2':['(00g)','(00g)0s0','(1/21/2g)','(1/21/2g)0s0',],
    'P -4 c 2':['(00g)','(00g)0s0','(1/21/2g)','(1/21/2g)0s0',],
    'P -4 b 2':['(00g)','(00g)0s0','(1/21/2g)','(1/21/2g)0q0',],
    'P -4 n 2':['(00g)','(00g)0s0','(1/21/2g)','(1/21/2g)0q0',],
    'P 4/m m m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)s0s0','(1/21/2g)00ss','(1/21/2g)s00s',],
    'P 4/m c c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)s0s0','(1/21/2g)00ss','(1/21/2g)s00s',],
    'P 4/n b m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)q0q0','(1/21/2g)q0qs',],
    'P 4/n n c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)q0q0','(1/21/2g)q0qs',],
    'P 4/m b m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 4/m n c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 4/n m m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 4/n c c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 42/m m c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)s0s0','(1/21/2g)00ss','(1/21/2g)s00s',],
    'P 42/m c m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)s0s0','(1/21/2g)00ss','(1/21/2g)s00s',],
    'P 42/n b c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)q0q0','(1/21/2g)q0qs',],
    'P 42/n n m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',
        '(1/21/2g)','(1/21/2g)q0q0','(1/21/2g)q0qs',],
    'P 42/m b c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 42/m n m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 42/n m c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 42/n c m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'I 4':['(00g)','(00g)q','(00g)s',],
    'I 41':['(00g)','(00g)q','(00g)s',],
    'I -4':['(00g)',],
    'I 4/m':['(00g)','(00g)s0',],
    'I 41/a':['(00g)','(00g)s0',],  #s0?
    'I 4 2 2':['(00g)','(00g)q00','(00g)s00',],
    'I 41 2 2':['(00g)','(00g)q00','(00g)s00',],
    'I 4 m m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'I 4 c m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'I 41 m d':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'I 41 c d':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'I -4 m 2':['(00g)','(00g)0s0',],
    'I -4 c 2':['(00g)','(00g)0s0',],
    'I -4 2 m':['(00g)','(00g)00s',],
    'I -4 2 d':['(00g)','(00g)00s',],
    'I 4/m m m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'I 4/m c m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'I 41/a m d':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'I 41/a c d':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    #trigonal/rhombahedral - done & checked
    'R 3':['(00g)','(00g)t',],
    'R -3':['(00g)',],
    'R 3 2':['(00g)','(00g)t0',],
    'R 3 m':['(00g)','(00g)0s',],
    'R 3 c':['(00g)','(00g)0s',],   #not 0s0
    'R -3 m':['(00g)','(00g)0s',],
    'R -3 c':['(00g)','(00g)0s',],  #not 0s0
    'P 3':['(00g)','(00g)t','(1/31/3g)','(1/31/3g)t',],
    'P 31':['(00g)','(00g)t','(1/31/3g)','(1/31/3g)t',],
    'P 32':['(00g)','(00g)t','(1/31/3g)','(1/31/3g)t',],
    'P -3':['(00g)','(1/31/3g)',],
    'P 3 1 2':['(00g)','(00g)t00','(1/31/3g)','(1/31/3g)t00',],
    'P 3 2 1':['(00g)','(00g)t00',],
    'P 31 1 2':['(00g)','(00g)t00','(1/31/3g)','(1/31/3g)t00',],
    'P 31 2 1':['(00g)','(00g)t00',],
    'P 32 1 2':['(00g)','(00g)t00','(1/31/3g)','(1/31/3g)t00',],
    'P 32 2 1':['(00g)','(00g)t00',],
    'P 3 m 1':['(00g)','(00g)0s0',],
    'P 3 1 m':['(00g)','(00g)00s','(1/31/3g)','(1/31/3g)00s',],
    'P 3 c 1':['(00g)','(00g)0s0',],
    'P 3 1 c':['(00g)','(00g)00s','(1/31/3g)','(1/31/3g)00s',],
    'P -3 1 m':['(00g)','(00g)00s','(1/31/3g)','(1/31/3g)00s',],
    'P -3 1 c':['(00g)','(00g)00s','(1/31/3g)','(1/31/3g)00s',],
    'P -3 m 1':['(00g)','(00g)0s0',],
    'P -3 c 1':['(00g)','(00g)0s0',],
    #hexagonal - done & checked
    'P 6':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P 61':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P 65':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P 62':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P 64':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P 63':['(00g)','(00g)h','(00g)t','(00g)s',],
    'P -6':['(00g)',],
    'P 6/m':['(00g)','(00g)s0',],
    'P 63/m':['(00g)','(00g)s0'],
    'P 6 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 61 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 65 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 62 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 64 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 63 2 2':['(00g)','(00g)h00','(00g)t00','(00g)s00',],
    'P 6 m m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'P 6 c c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'P 63 c m':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'P 63 m c':['(00g)','(00g)ss0','(00g)0ss','(00g)s0s',],
    'P -6 m 2':['(00g)','(00g)0s0',],
    'P -6 c 2':['(00g)','(00g)0s0',],
    'P -6 2 m':['(00g)','(00g)00s',],
    'P -6 2 c':['(00g)','(00g)00s',],
    'P 6/m m m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 6/m c c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 63/m c m':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    'P 63/m m c':['(00g)','(00g)s0s0','(00g)00ss','(00g)s00s',],
    }

#'A few non-standard space groups for test use'
nonstandard_sglist = ('P 21 1 1','P 1 21 1','P 1 1 21','R 3 r','R 3 2 h', 
                      'R -3 r', 'R 3 2 r','R 3 m h', 'R 3 m r',
                      'R 3 c r','R -3 c r','R -3 m r',),

#A list of orthorhombic space groups that were renamed in the 2002 Volume A,
# along with the pre-2002 name. The e designates a double glide-plane'''
sgequiv_2002_orthorhombic= (('A e m 2', 'A b m 2',),
                            ('A e a 2', 'A b a 2',),
                            ('C m c e', 'C m c a',),
                            ('C m m e', 'C m m a',),
                            ('C c c e', 'C c c a'),)
#Use the space groups types in this order to list the symbols in the 
#order they are listed in the International Tables, vol. A'''
symtypelist = ('triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 
               'trigonal', 'hexagonal', 'cubic')

# self-test materials follow. Requires files in directory testinp
selftestlist = []
'''Defines a list of self-tests'''
selftestquiet = True
def _ReportTest():
    'Report name and doc string of current routine when ``selftestquiet`` is False'
    if not selftestquiet:
        import inspect
        caller = inspect.stack()[1][3]
        doc = eval(caller).__doc__
        if doc is not None:
            print('testing '+__file__+' with '+caller+' ('+doc+')')
        else:
            print('testing '+__file__()+" with "+caller)
def test0():
    '''self-test #0: exercise MoveToUnitCell'''
    _ReportTest()
    msg = "MoveToUnitCell failed"
    assert (MoveToUnitCell([1,2,3]) == [0,0,0]).all, msg
    assert (MoveToUnitCell([2,-1,-2]) == [0,0,0]).all, msg
    assert abs(MoveToUnitCell(np.array([-.1]))[0]-0.9) < 1e-6, msg
    assert abs(MoveToUnitCell(np.array([.1]))[0]-0.1) < 1e-6, msg
selftestlist.append(test0)

def test1():
    '''self-test #1: SpcGroup against previous results'''
    #'''self-test #1: SpcGroup and SGPrint against previous results'''
    _ReportTest()
    testdir = ospath.join(ospath.split(ospath.abspath( __file__ ))[0],'testinp')
    if ospath.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    import spctestinp
    def CompareSpcGroup(spc, referr, refdict, reflist): 
        'Compare output from GSASIIspc.SpcGroup with results from a previous run'
        # if an error is reported, the dictionary can be ignored
        msg0 = "CompareSpcGroup failed on space group %s" % spc
        result = SpcGroup(spc)
        if result[0] == referr and referr > 0: return True
        keys = result[1].keys()
        #print result[1]['SpGrp']
        #msg = msg0 + " in list lengths"
        #assert len(keys) == len(refdict.keys()), msg
        for key in refdict.keys():
            if key == 'SGOps' or  key == 'SGCen':
                msg = msg0 + (" in key %s length" % key)
                assert len(refdict[key]) == len(result[1][key]), msg
                for i in range(len(refdict[key])):
                    msg = msg0 + (" in key %s level 0" % key)
                    assert np.allclose(result[1][key][i][0],refdict[key][i][0]), msg
                    msg = msg0 + (" in key %s level 1" % key)
                    assert np.allclose(result[1][key][i][1],refdict[key][i][1]), msg
            else:
                msg = msg0 + (" in key %s" % key)
                assert result[1][key] == refdict[key], msg
        msg = msg0 + (" in key %s reflist" % key)
        #for (l1,l2) in zip(reflist, SGPrint(result[1])):
        #    assert l2.replace('\t','').replace(' ','') == l1.replace(' ',''), 'SGPrint ' +msg
        # for now disable SGPrint testing, output has changed
        #assert reflist == SGPrint(result[1]), 'SGPrint ' +msg
    for spc in spctestinp.SGdat:
        CompareSpcGroup(spc, 0, spctestinp.SGdat[spc], spctestinp.SGlist[spc] )
selftestlist.append(test1)

def test2():
    '''self-test #2: SpcGroup against cctbx (sgtbx) computations'''
    _ReportTest()
    testdir = ospath.join(ospath.split(ospath.abspath( __file__ ))[0],'testinp')
    if ospath.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    import sgtbxtestinp
    def CompareWcctbx(spcname, cctbx_in, debug=0):
        'Compare output from GSASIIspc.SpcGroup with results from cctbx.sgtbx'
        cctbx = cctbx_in[:] # make copy so we don't delete from the original
        spc = (SpcGroup(spcname))[1]
        if debug: print spc['SpGrp']
        if debug: print spc['SGCen']
        latticetype = spcname.strip().upper()[0]
        # lattice type of R implies Hexagonal centering", fix the rhombohedral settings
        if latticetype == "R" and len(spc['SGCen']) == 1: latticetype = 'P'
        assert latticetype == spc['SGLatt'], "Failed: %s does not match Lattice: %s" % (spcname, spc['SGLatt'])
        onebar = [1]
        if spc['SGInv']: onebar.append(-1)
        for (op,off) in spc['SGOps']:
            for inv in onebar:
                for cen in spc['SGCen']:
                    noff = off + cen
                    noff = MoveToUnitCell(noff)
                    mult = tuple((op*inv).ravel().tolist())
                    if debug: print "\n%s: %s + %s" % (spcname,mult,noff)
                    for refop in cctbx:
                        if debug: print refop
                        # check the transform
                        if refop[:9] != mult: continue
                        if debug: print "mult match"
                        # check the translation
                        reftrans = list(refop[-3:])
                        reftrans = MoveToUnitCell(reftrans)
                        if all(abs(noff - reftrans) < 1.e-5):
                            cctbx.remove(refop)
                            break
                    else:
                        assert False, "failed on %s:\n\t %s + %s" % (spcname,mult,noff)
    for key in sgtbxtestinp.sgtbx:
        CompareWcctbx(key, sgtbxtestinp.sgtbx[key])
selftestlist.append(test2)

def test3(): 
    '''self-test #3: exercise SytSym (includes GetOprPtrName, GenAtom, GetKNsym)
     for selected space groups against info in IT Volume A '''
    _ReportTest()
    def ExerciseSiteSym (spc, crdlist):
        'compare site symmetries and multiplicities for a specified space group'
        msg = "failed on site sym test for %s" % spc
        (E,S) = SpcGroup(spc)
        assert not E, msg
        for t in crdlist:
            symb, m = SytSym(t[0],S)
            if symb.strip() != t[2].strip() or m != t[1]:
                print spc,t[0],m,symb,t[2]
            assert m == t[1]
            #assert symb.strip() == t[2].strip()

    ExerciseSiteSym('p 1',[
            ((0.13,0.22,0.31),1,'1'),
            ((0,0,0),1,'1'),
            ])
    ExerciseSiteSym('p -1',[
            ((0.13,0.22,0.31),2,'1'),
            ((0,0.5,0),1,'-1'),
            ])
    ExerciseSiteSym('C 2/c',[
            ((0.13,0.22,0.31),8,'1'),
            ((0.0,.31,0.25),4,'2(y)'),
            ((0.25,.25,0.5),4,'-1'),
            ((0,0.5,0),4,'-1'),
            ])
    ExerciseSiteSym('p 2 2 2',[
            ((0.13,0.22,0.31),4,'1'),
            ((0,0.5,.31),2,'2(z)'),
            ((0.5,.31,0.5),2,'2(y)'),
            ((.11,0,0),2,'2(x)'),
            ((0,0.5,0),1,'222'),
            ])
    ExerciseSiteSym('p 4/n',[
            ((0.13,0.22,0.31),8,'1'),
            ((0.25,0.75,.31),4,'2(z)'),
            ((0.5,0.5,0.5),4,'-1'),
            ((0,0.5,0),4,'-1'),
            ((0.25,0.25,.31),2,'4(001)'),
            ((0.25,.75,0.5),2,'-4(001)'),
            ((0.25,.75,0.0),2,'-4(001)'),
            ])
    ExerciseSiteSym('p 31 2 1',[
            ((0.13,0.22,0.31),6,'1'),
            ((0.13,0.0,0.833333333),3,'2(100)'),
            ((0.13,0.13,0.),3,'2(110)'),
            ])
    ExerciseSiteSym('R 3 c',[
            ((0.13,0.22,0.31),18,'1'),
            ((0.0,0.0,0.31),6,'3'),
            ])
    ExerciseSiteSym('R 3 c R',[
            ((0.13,0.22,0.31),6,'1'),
            ((0.31,0.31,0.31),2,'3(111)'),
            ])
    ExerciseSiteSym('P 63 m c',[
            ((0.13,0.22,0.31),12,'1'),
            ((0.11,0.22,0.31),6,'m(100)'),
            ((0.333333,0.6666667,0.31),2,'3m(100)'),
            ((0,0,0.31),2,'3m(100)'),
            ])
    ExerciseSiteSym('I a -3',[
            ((0.13,0.22,0.31),48,'1'),
            ((0.11,0,0.25),24,'2(x)'),
            ((0.11,0.11,0.11),16,'3(111)'),
            ((0,0,0),8,'-3(111)'),
            ])
selftestlist.append(test3)

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    for test in selftestlist:
        test()
    print "OK"
