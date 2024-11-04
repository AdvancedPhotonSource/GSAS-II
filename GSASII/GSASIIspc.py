# -*- coding: utf-8 -*-
''':mod:`GSASIIspc` Classes & routines follow
'''
from __future__ import division, print_function
import numpy as np
import numpy.linalg as nl
import scipy.optimize as so
import sys
import copy
import os.path as ospath

import GSASIIpath

npsind = lambda x: np.sin(x*np.pi/180.)
npcosd = lambda x: np.cos(x*np.pi/180.)
nxs = np.newaxis
DEBUG = False
    
################################################################################
#### Space group codes
################################################################################

def SpcGroup(SGSymbol):
    """
    Determines cell and symmetry information from a short H-M space group name

    :param SGSymbol: space group symbol (string) with spaces between axial fields
    :returns: (SGError,SGData)
    
       * SGError = 0 for no errors; >0 for errors (see :func:`SGErrors` for details)
       * SGData - is a dict (see :ref:`Space Group object<SGData_table>`) with entries:
       
             * 'SpGrp': space group symbol, slightly cleaned up
             * 'SGFixed': True if space group data can not be changed, e.g. from magnetic cif; otherwise False
             * 'SGGray': True if 1' in symbol - gray group for mag. incommensurate phases
             * 'SGLaue':  one of '-1', '2/m', 'mmm', '4/m', '4/mmm', '3R',
               '3mR', '3', '3m1', '31m', '6/m', '6/mmm', 'm3', 'm3m'
             * 'SGInv': boolean; True if centrosymmetric, False if not
             * 'SGLatt': one of 'P', 'A', 'B', 'C', 'I', 'F', 'R'
             * 'SGUniq': one of 'a', 'b', 'c' if monoclinic, '' otherwise
             * 'SGCen': cell centering vectors [0,0,0] at least
             * 'SGOps': symmetry operations as [M,T] so that M*x+T = x'
             * 'SGSys': one of 'triclinic', 'monoclinic', 'orthorhombic',
               'tetragonal', 'rhombohedral', 'trigonal', 'hexagonal', 'cubic'
             * 'SGPolax': one of ' ', 'x', 'y', 'x y', 'z', 'x z', 'y z',
               'xyz', '111' for arbitrary axes
             * 'SGPtGrp': one of 32 point group symbols (with some permutations), which
               is filled by SGPtGroup, is external (KE) part of supersymmetry point group
             * 'SSGKl': default internal (Kl) part of supersymmetry point group; modified 
               in supersymmetry stuff depending on chosen modulation vector for Mono & Ortho
             * 'BNSlattsym': BNS lattice symbol & cenering op - used for magnetic structures

    """
    LaueSym = ('-1','2/m','mmm','4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm','m3','m3m')
    LattSym = ('P','A','B','C','I','F','R')
    UniqSym = ('','','a','b','c','',)
    SysSym = ('triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','trigonal','hexagonal','cubic')
    SGData = {}
    if len(SGSymbol.split()) < 2:
        return SGErrors(0),SGData
    if ':R' in SGSymbol:
        SGSymbol = SGSymbol.replace(':',' ')    #get rid of ':' in R space group symbols from some cif files
    SGData['SGGray'] = False
    if "1'" in SGSymbol:        #set for incommensurate magnetic
        SGData['SGGray'] = True
        SGSymbol = SGSymbol.replace("1'",'')
    SGSymbol = SGSymbol.split(':')[0]   #remove :1/2 setting symbol from some cif files
    if '-2' in SGSymbol:    #replace bad but legal symbols with correct equivalents
        SGSymbol = SGSymbol.replace('-2','m')
    if SGSymbol.split()[1] =='3/m':
        SGSymbol = SGSymbol.replace('3/m','-6')
    import pyspg
    SGInfo = pyspg.sgforpy(SGSymbol)
    SGData['SpGrp'] = SGSymbol.strip().lower().capitalize()
    SGData['SGLaue'] = LaueSym[SGInfo[0]-1]
    SGData['SGInv'] = bool(SGInfo[1])
    SGData['SGLatt'] = LattSym[SGInfo[2]-1]
    SGData['SGUniq'] = UniqSym[SGInfo[3]+1]
    SGData['SGFixed'] = False
    SGData['SGOps'] = []
    SGData['SGGen'] = []
    for i in range(SGInfo[5]):
        Mat = np.array(SGInfo[6][i])
        Trns = np.array(SGInfo[7][i])
        SGData['SGOps'].append([Mat,Trns])
        if 'array' in str(type(SGInfo[8])):        #patch for old fortran bin?
            SGData['SGGen'].append(int(SGInfo[8][i]))
    SGData['BNSlattsym'] = [LattSym[SGInfo[2]-1],[0,0,0]]
    lattSpin = []
    if SGData['SGLatt'] == 'P':
        SGData['SGCen'] = np.array(([0,0,0],))
    elif SGData['SGLatt'] == 'A':
        SGData['SGCen'] = np.array(([0,0,0],[0,.5,.5]))
        lattSpin += [1,]
    elif SGData['SGLatt'] == 'B':
        SGData['SGCen'] = np.array(([0,0,0],[.5,0,.5]))
        lattSpin += [1,]
    elif SGData['SGLatt'] == 'C':
        SGData['SGCen'] = np.array(([0,0,0],[.5,.5,0,]))
        lattSpin += [1,]
    elif SGData['SGLatt'] == 'I':
        SGData['SGCen'] = np.array(([0,0,0],[.5,.5,.5]))
        lattSpin += [1,]
    elif SGData['SGLatt'] == 'F':
        SGData['SGCen'] = np.array(([0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0,]))
        lattSpin += [1,1,1,1]
    elif SGData['SGLatt'] == 'R':
        SGData['SGCen'] = np.array(([0,0,0],[2./3,1./3,1./3],[1./3,2./3,2./3]))

    if SGData['SGInv']:
        if SGData['SGLaue'] in ['-1','2/m','mmm']:
            Ibar = 7
        elif SGData['SGLaue'] in ['4/m','4/mmm']:
            Ibar = 1
        elif SGData['SGLaue'] in ['3R','3mR','3','3m1','31m','6/m','6/mmm']:
            Ibar = 15 #8+4+2+1
        else:
            Ibar = 4
        Ibarx = Ibar&14
    else:
        Ibarx = 8
        if SGData['SGLaue'] in ['-1','2/m','mmm','m3','m3m']:
            Ibarx = 0
    moregen = []
    for i,gen in enumerate(SGData['SGGen']):
        if SGData['SGLaue'] in ['m3','m3m']:
            if gen in [1,2,4]:
                SGData['SGGen'][i] = 4
            elif gen < 7:
                SGData['SGGen'][i] = 0
        elif SGData['SGLaue'] in ['4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm']:
            if gen == 2:
                SGData['SGGen'][i] = 4
            elif gen in [3,5]:
                SGData['SGGen'][i] = 3
            elif gen == 6:
                if SGData['SGLaue'] in ['4/m','4/mmm']:
                    SGData['SGGen'][i] = 128
                else:
                    SGData['SGGen'][i] = 16
            elif not SGData['SGInv'] and gen == 12:
                SGData['SGGen'][i] = 8
            elif (not SGData['SGInv']) and (SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']) and (gen == 1):
                SGData['SGGen'][i] = 24
        gen1 = SGData['SGGen'][i]
        if gen1 == 99:
            gen1 = 8
            if SGData['SGLaue'] in ['3m1','31m','6/m','6/mmm']:
                gen1 = 3
            elif SGData['SGLaue'] == 'm3m':
                gen1 = 12
            SGData['SGGen'][i] = gen1
        elif gen1 == 98:
            gen1 = 8
            if SGData['SGLaue'] in ['3m1','31m','6/m','6/mmm']:
                gen1 = 4
            SGData['SGGen'][i] = gen1
        elif not SGData['SGInv'] and gen1 in [23,] and SGData['SGLaue'] in ['m3','m3m']:
            SGData['SGGen'][i] = 24
        elif gen1 >= 16 and gen1 != 128:
            if not SGData['SGInv']:
                gen1 = 31
            else:
                gen1 ^= Ibarx 
            SGData['SGGen'][i] = gen1
        if SGData['SGInv']:
            if gen1 < 128:
                moregen.append(SGData['SGGen'][i]^Ibar)
            else:
                moregen.append(1)
    SGData['SGGen'] += moregen
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

    if SGData['SGLatt'] == 'R':
        if SGData['SGPtGrp'] in ['3',]:
            SGData['SGSpin'] = 3*[1,]
        elif SGData['SGPtGrp'] in  ['-3','32','3m']:
            SGData['SGSpin'] = 4*[1,]
        elif SGData['SGPtGrp'] in  ['-3m',]:
            SGData['SGSpin'] = 5*[1,]
        
    else:
        if SGData['SGPtGrp'] in ['1','3','23',]:
            SGData['SGSpin'] = lattSpin+[1,]
        elif SGData['SGPtGrp'] in ['-1','2','m','4','-4','-3','312','321','3m1','31m','6','-6','432','-43m']:
            SGData['SGSpin'] = lattSpin+[1,1,]
        elif SGData['SGPtGrp'] in ['2/m','4/m','422','4mm','-42m','-4m2','-3m1','-31m',
            '6/m','622','6mm','-6m2','-62m','m3','m3m']:
            SGData['SGSpin'] = lattSpin+[1,1,1,]
        else: #'222'-'mmm','4/mmm','6/mmm'
            SGData['SGSpin'] = lattSpin+[1,1,1,1,]
    return SGInfo[-1],SGData

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
    
    :param SGData: from :func:`SpcGroup`
    :returns: SSGPtGrp & SSGKl (only defaults for Mono & Ortho)
    '''
    Flds = SGData['SpGrp'].split()
    if len(Flds) < 2:
        return '',[]
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
            return '222',[-1,-1,-1]
        elif SGData['SpGrp'].count('2') == 1:
            if SGData['SGPolax'] == 'x':
                return '2mm',[-1,1,1]
            elif SGData['SGPolax'] == 'y':
                return 'm2m',[1,-1,1]
            elif SGData['SGPolax'] == 'z':
                return 'mm2',[1,1,-1]
        else:
            return 'mmm',[1,1,1]
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
            return 'm3m',[]
    
def SGPrint(SGData,AddInv=False):
    '''
    Print the output of SpcGroup in a nicely formatted way. Used in SpaceGroup

    :param SGData: from :func:`SpcGroup`
    :returns:
        SGText - list of strings with the space group details
        SGTable - list of strings for each of the operations
    '''
    if SGData.get('SGFixed',False):       #inverses included in ops for cif fixed
        Mult = len(SGData['SGCen'])*len(SGData['SGOps'])
    else:
        Mult = len(SGData['SGCen'])*len(SGData['SGOps'])*(int(SGData['SGInv'])+1)
    SGText = []
    SGText.append(' Space Group: '+SGData['SpGrp'])
    SGCen = list(SGData['SGCen'])
    if SGData.get('SGGray',False):
        SGText[-1] += " 1'"
        if SGData.get('SGFixed',False): 
            Mult //= 2
        else:
            SGCen += list(SGData['SGCen']+[0,0,0])
            SGCen =  np.array(SGCen)%1.
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
    if len(SGData['SGCen']) == 1:
        SGText.append(' The equivalent positions are:\n')
    else:    
        SGText.append(' The equivalent positions are:\n')
        SGText.append(' ('+Latt2text(SGCen)+')+\n')
    SGTable = []
    for i,Opr in enumerate(SGData['SGOps']):
        SGTable.append('(%2d) %s'%(i+1,MT2text(Opr)))
    if AddInv and SGData['SGInv']:
        for i,Opr in enumerate(SGData['SGOps']):
            IOpr = [-Opr[0],-Opr[1]]
            SGTable.append('(%2d) %s'%(i+1,MT2text(IOpr)))        
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
      * G2oprList: a list with the GSAS-II operations for each symmetry operation as
        a tuple with (center,mult,opnum,opcode), where center is (0,0,0), (0.5,0,0),
        (0.5,0.5,0.5),...; where mult is 1 or -1 for the center of symmetry
        where opnum is the number for the symmetry operation, in ``SGOps``
        (starting with 0) and opcode is mult*(100*icen+j+1).
      * G2opcodes: a list with the name that GSAS-II uses for each symmetry 
        operation (same as opcode, above)
    '''
    SGTextList = []
    offsetList = []
    symOpList = []
    G2oprList = []
    G2opcodes = []
    onebar = (1,)
    if SGData['SGInv']:
        onebar += (-1,)
    for icen,cen in enumerate(SGData['SGCen']):
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
                G2opcodes.append(mult*(100*icen+j+1))
    return SGTextList,offsetList,symOpList,G2oprList,G2opcodes

def TextOps(text,table,reverse=False):
    ''' Makes formatted operator list
        :param text,table: arrays of text made by SGPrint
        :param reverse: True for x+1/2 form; False for 1/2+x form
        :returns: OpText: full list of symmetry operators; one operation per line
        generally printed to console for use via cut/paste in other programs, but
        could be used for direct input
    '''
    OpText = []
    Inv = True
    Super = False
    if 'noncentro' in text[1]:
        Inv = False
    if 'Super' in text[0]:
        Super = True
    Cent = [[0,0,0],]
    if Super:
        Cent = [[0,0,0,0],]
    if '0,0,0' in text[-1]:
        Cent = np.array(eval(text[-1].split('+')[0].replace(';','),(')))
    OpsM = []
    OpsT = []
    for item in table:
        if 'for' in item: continue
        if Super:
            M,T = MagSSText2MTS(item.split(')')[1].replace(' ',''),G2=True)[:2]
        else:
            M,T = Text2MT(item.split(')')[1].replace(' ',''),CIF=True)
        OpsM.append(M)
        OpsT.append(T)
    OpsM = np.array(OpsM)
    OpsT = np.array(OpsT)
    if Inv and not Super:   #inversion ops altready listed in supersymmetries
        OpsM = np.concatenate((OpsM,-OpsM))
        OpsT = np.concatenate((OpsT,-OpsT%1.))
    for cent in Cent:
        for iop,opM in enumerate(list(OpsM)):
            if Super:
                txt = SSMT2text([opM,(OpsT[iop]+cent)%1.])
            else:
                txt = MT2text([opM,(OpsT[iop]+cent)%1.],reverse)
            OpText.append(txt.replace(' ','').lower())
    return OpText

def TextGen(SGData,reverse=False):      #does not always work correctly - not used anyway
    GenSym,GenFlg,BNSsym = GetGenSym(SGData)
    SGData['GenSym'] = GenSym
    SGData['GenFlg'] = GenFlg
    text,table = SGPrint(SGData)
    GenText = []
    OprNames = GetOprNames(SGData)
    OpText = TextOps(text,table,reverse)
    for name in SGData['GenSym']:
        gid = OprNames.index(name.replace(' ',''))
        GenText.append(OpText[gid])
    if len(SGData['SGCen']) > 1:
        GenText.append(OpText[-1])
    return GenText

def GetOprNames(SGData):
    OprNames = [GetOprPtrName(str(irtx)) for irtx in PackRot(SGData['SGOps'])]
    if SGData['SGInv']:
        OprNames += [GetOprPtrName(str(-irtx)) for irtx in PackRot(SGData['SGOps'])]
    return OprNames
    
def MT2text(Opr,reverse=False):
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
                if reverse:
                    Fld += (XYZ[IJ]+'+'+TRA[IK]).rjust(5)
                else:
                    Fld += (TRA[IK]+XYZ[IJ]).rjust(5)
            else:
                if reverse:
                    Fld += (XYZ[IJ]+'+'+TRA[IK]).rjust(5)
                else:
                    Fld += (TRA[IK]+'+'+XYZ[IJ]).rjust(5)
        else:
            Fld += XYZ[IJ].rjust(5)
        if j != 2: Fld += ', '
    return Fld
    
def Latt2text(Cen):
    "From lattice centering vectors returns ';' delimited cell centering vectors"
    lattTxt = ''
    fracList = ['1/2','1/3','2/3','1/4','3/4','1/5','2/5','3/5','4/5','1/6','5/6',
        '1/7','2/7','3/7','4/7','5/7','6/7','1/8','3/8','5/8','7/8','1/9','2/9','4/9','5/9','7/9','8/9']
    mulList = [2,3,3,4,4,5,5,5,5,6,6,7,7,7,7,7,7,8,8,8,8,9,9,9,9,9,9]
    prodList = [1.,1.,2.,1.,3.,1.,2.,3.,4.,1.,5.,1.,2.,3.,4.,5.,6.,1.,3.,5.,7.,1.,2.,4.,5.,7.,8.]
    nCen = len(Cen)
    for i,cen in enumerate(Cen):
        txt = ''
        for icen in cen:
            if icen == 1:
                txt += '1,'
                continue
            elif icen == -1:
                txt += '-1,'
                continue
            if not icen:
                txt += '0,'
                continue
            if icen < 0:
                txt += '-'
                icen *= -1
            for mul,prod,frac in zip(mulList,prodList,fracList):
                if abs(icen*mul-prod) < 1.e-5:
                    txt += frac+','
                    break
        lattTxt += txt[:-1]+'; '
        if i and not i%8 and i < nCen-1:    #not for the last cen!
            lattTxt += '\n     '
    return lattTxt[:-2]

def SpaceGroup(SGSymbol):
    '''
    Print the output of SpcGroup in a nicely formatted way. 

    :param SGSymbol: space group symbol (string) with spaces between axial fields
    :returns: nothing
    '''
    E,A = SpcGroup(SGSymbol)
    if E > 0:
        print (SGErrors(E))
        return
    for l in SGPrint(A):
        print (l)
################################################################################
#### Magnetic space group stuff
################################################################################

def SplitMagSpSG(MSpSg):
    SpSG = StandardizeSpcName(MSpSg.replace("'",''))    #removes all mag. op. marks & makes axial fields
    Sflds = SpSG.split()                #2-4 axial fields
    Axf = []
    MSpSg += ' '
    Psym = MSpSg[0]
    Ib = 0
    if '_' in MSpSg:
        Psym = MSpSg[:3]
        Ib = 2
    Ib += 1
    if len(Sflds) == 2:                 #simple: one axial field
        return Psym+' '+MSpSg[Ib:]
    if '/' in MSpSg[Ib:]:               #more fields - do 1st axial field
        Is = MSpSg.index('/')
        If = Is+1
        if "'" in MSpSg[Is+2]:
            If += 1
        Axf.append(MSpSg[Ib:If+1])
        Ib = If+1
    else:
        If = Ib+len(Sflds[1])
        if "'" == MSpSg[If]:
            If += 1
        Axf.append(MSpSg[Ib:If])
        Ib = If
    for i in range(len(Sflds)-2):   #do rest
        If = Ib+len(Sflds[i+2])
        if "'" == MSpSg[If]:
            If += 1
        Axf.append(MSpSg[Ib:If])
        Ib = If
    MSpcGp = Psym+' '+' '.join(Axf)
    return MSpcGp
        
def SetMagnetic(SGData):
    GenSym,GenFlg,BNSsym = GetGenSym(SGData)
    SGData['GenSym'] = GenSym
    SGData['GenFlg'] = GenFlg
    OprNames,SpnFlp = GenMagOps(SGData)
    SGData['SpnFlp'] = SpnFlp
    SGData['MagSpGrp'] = MagSGSym(SGData)

def GetGenSym(SGData):
    '''
    Get the space group generator symbols
    :param SGData: from :func:`SpcGroup`
    LaueSym = ('-1','2/m','mmm','4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm','m3','m3m')
    LattSym = ('P','A','B','C','I','F','R')
    
    '''
    BNSsym = {}
    OprNames = [GetOprPtrName(str(irtx)) for irtx in PackRot(SGData['SGOps'])]
    if SGData['SGInv']:
        OprNames += [GetOprPtrName(str(-irtx)) for irtx in PackRot(SGData['SGOps'])]
    # for oprname in OprNames: 
    #     print(oprname)
    Nsyms = len(SGData['SGOps'])
    if SGData['SGInv'] and not SGData['SGFixed']: Nsyms *= 2
    UsymOp = ['1',]
    OprFlg = [0,]  
    if Nsyms == 2:                    #Centric triclinic or acentric monoclinic
        UsymOp.append(OprNames[1])
        OprFlg.append(SGData['SGGen'][1])
    elif Nsyms == 4:                    #Point symmetry 2/m, 222, 22m, or 4
        if '4z' in OprNames[1]:          #Point symmetry 4 or -4
            UsymOp.append(OprNames[1])
            OprFlg.append(SGData['SGGen'][1])
        elif not SGData['SGInv']:       #Acentric Orthorhombic
            if 'm' in OprNames[1:4]:    #22m, 2m2 or m22
                if '2' in OprNames[1]:      #Acentric orthorhombic, 2mm
                    UsymOp.append(OprNames[2])
                    OprFlg.append(SGData['SGGen'][2])
                    UsymOp.append(OprNames[3])
                    OprFlg.append(SGData['SGGen'][3])
                elif '2' in OprNames[2]:    #Acentric orthorhombic, m2m
                    UsymOp.append(OprNames[1])
                    OprFlg.append(SGData['SGGen'][1])
                    UsymOp.append(OprNames[3])
                    OprFlg.append(SGData['SGGen'][3])
                else:                       #Acentric orthorhombic, mm2
                    UsymOp.append(OprNames[1])
                    OprFlg.append(SGData['SGGen'][1])
                    UsymOp.append(OprNames[2])
                    OprFlg.append(SGData['SGGen'][2])
            else:                           #Acentric orthorhombic, 222
                SGData['SGGen'][1:] = [4,2,1]
                UsymOp.append(OprNames[1])
                OprFlg.append(SGData['SGGen'][1])
                UsymOp.append(OprNames[2])
                OprFlg.append(SGData['SGGen'][2])
                UsymOp.append(OprNames[3])
                OprFlg.append(SGData['SGGen'][3])
        else:                               #Centric Monoclinic
            UsymOp.append(OprNames[1])
            OprFlg.append(SGData['SGGen'][1])
            UsymOp.append(OprNames[3])
            OprFlg.append(SGData['SGGen'][3])
    elif Nsyms == 6:                    #Point symmetry 32, 3m or 6
            if '6' in OprNames[1]:      #Hexagonal 6 Laue symmetry
                UsymOp.append(OprNames[1])
                OprFlg.append(SGData['SGGen'][1])
            else:                       #Trigonal
                UsymOp.append(OprNames[4])
                OprFlg.append(SGData['SGGen'][3])
                if '2110' in OprNames[1]: UsymOp[-1] = ' 2100 '
    elif Nsyms == 8:                    #Point symmetry mmm, 4/m, or 422, etc
        if '4' in OprNames[1]:           #Tetragonal
            if SGData['SGInv']:         #4/m
                UsymOp.append(OprNames[1])
                OprFlg.append(SGData['SGGen'][1])
                UsymOp.append(OprNames[6])
                OprFlg.append(SGData['SGGen'][6])
            else:
                if 'x' in OprNames[4]:      #4mm type group
                    UsymOp.append(OprNames[4])
                    OprFlg.append(6)
                    UsymOp.append(OprNames[7])
                    OprFlg.append(8)
                else:                       #-42m, -4m2, and 422 type groups
                    UsymOp.append(OprNames[6])
                    OprFlg.append(19)
                    UsymOp.append(OprNames[5])
                    OprFlg.append(8)
        else:                               #Orthorhombic, mmm
            UsymOp.append(OprNames[1])
            OprFlg.append(SGData['SGGen'][1])
            UsymOp.append(OprNames[2])
            OprFlg.append(SGData['SGGen'][2])
            UsymOp.append(OprNames[7])
            OprFlg.append(SGData['SGGen'][7])
    elif Nsyms == 12 and '3' in OprNames[1] and SGData['SGSys'] != 'cubic':        #Trigonal
        UsymOp.append(OprNames[3])
        OprFlg.append(SGData['SGGen'][3])
        UsymOp.append(OprNames[9])
        OprFlg.append(SGData['SGGen'][9])
    elif Nsyms == 12 and '6' in OprNames[1]:        #Hexagonal
        if 'mz' in OprNames[9]:                     #6/m
            UsymOp.append(OprNames[1])
            OprFlg.append(SGData['SGGen'][1])
            UsymOp.append(OprNames[9])
            OprFlg.append(SGData['SGGen'][9])
        else:                                       #6mm, -62m, -6m2 or 622
            UsymOp.append(OprNames[6])
            OprFlg.append(18)
            if 'm' in UsymOp[-1]: OprFlg[-1] = 20
            UsymOp.append(OprNames[7])
            OprFlg.append(24)
    elif Nsyms in [16,24]:
        if '3' in OprNames[1]:
            UsymOp.append('')
            OprFlg.append(SGData['SGGen'][3])
            for i in range(Nsyms):
                if 'mx' in OprNames[i]:
                    UsymOp[-1] = OprNames[i]
                elif 'm11' in OprNames[i]:
                    UsymOp[-1] = OprNames[i]
                elif '211' in OprNames[i]:
                    UsymOp[-1] = OprNames[i]
                    OprFlg[-1] = 24
        else:                                     #4/mmm or 6/mmm
            UsymOp.append('  mz  ')
            OprFlg.append(1)
            if '4' in OprNames[1]:                  #4/mmm
                UsymOp.append('  mx  ')
                OprFlg.append(20)
                UsymOp.append(' m110 ')
                OprFlg.append(24)
            else:                                   #6/mmm
                UsymOp.append(' m110 ')
                OprFlg.append(4)
                UsymOp.append(' m+-0 ')
                OprFlg.append(8)
    else:                                           #System is cubic
        if Nsyms == 48:
            UsymOp.append('  mx  ')
            OprFlg.append(4)
            UsymOp.append(' m110 ')
            OprFlg.append(24)
            
    if 'P' in SGData['SGLatt']:
        if SGData['SGSys'] == 'triclinic':
            BNSsym = {'P_a':[.5,0,0],'P_b':[0,.5,0],'P_c':[0,0,.5]}            
        elif SGData['SGSys'] == 'monoclinic':
            BNSsym = {'P_a':[.5,0,0],'P_b':[0,.5,0],'P_c':[0,0,.5],'P_I':[.5,.5,.5]}
            if SGData['SGUniq'] == 'a':
                BNSsym.update({'P_B':[.5,0,.5],'P_C':[.5,.5,0]})
            elif SGData['SGUniq'] == 'b':
                BNSsym.update({'P_A':[.5,.5,0],'P_C':[0,.5,.5]})
            elif SGData['SGUniq'] == 'c':
                BNSsym.update({'P_A':[0,.5,.5],'P_B':[.5,0,.5]})
        elif SGData['SGSys'] == 'orthorhombic':
            BNSsym = {'P_a':[.5,0,0],'P_b':[0,.5,0],'P_c':[0,0,.5],
                'P_A':[0,.5,.5],'P_B':[.5,0,.5],'P_C':[.5,.5,0],'P_I':[.5,.5,.5]}
        elif SGData['SGSys'] == 'tetragonal':
            BNSsym = {'P_c':[0,0,.5],'P_C':[.5,.5,0],'P_I':[.5,.5,.5]}            
        elif SGData['SGSys'] in ['trigonal','hexagonal']:
            BNSsym = {'P_c':[0,0,.5]}            
        elif SGData['SGSys'] == 'cubic':
            BNSsym = {'P_I':[.5,.5,.5]}            
            
    elif 'A' in SGData['SGLatt']:
        if SGData['SGSys'] == 'monoclinic':
            BNSsym = {}
            if SGData['SGUniq'] == 'b':
                BNSsym.update({'A_a':[.5,0,0],'A_c':[0,0,.5]})
            elif SGData['SGUniq'] == 'c':
                BNSsym.update({'A_a':[.5,0,0],'A_b':[0,.5,0]})
        elif SGData['SGSys'] == 'orthorhombic':
            BNSsym = {'A_a':[.5,0,0],'A_b':[0,.5,0],'A_c':[0,0,.5],
               'A_B':[.5,0,.5],'A_C':[.5,.5,0]}   
        elif SGData['SGSys'] == 'triclinic':
            BNSsym = {'A_a':[.5,0,0],'A_b':[0,.5,0],'A_c':[0,0,.5]}   
            
    elif 'B' in SGData['SGLatt']:
        if SGData['SGSys'] == 'monoclinic':
            BNSsym = {}
            if SGData['SGUniq'] == 'a':
                BNSsym.update({'B_b':[0,.5,0],'B_c':[0,0,.5]})
            elif SGData['SGUniq'] == 'c':
                BNSsym.update({'B_a':[.5,0,0],'B_b':[0,.5,0]})
        elif SGData['SGSys'] == 'orthorhombic':
            BNSsym = {'B_a':[.5,0,0],'B_b':[0,.5,0],'B_c':[0,0,.5],
                'B_A':[0,.5,.5],'B_C':[.5,.5,0]}     
        elif SGData['SGSys'] == 'triclinic':
            BNSsym = {'B_a':[.5,0,0],'B_b':[0,.5,0],'B_c':[0,0,.5]}     
            
    elif 'C' in SGData['SGLatt']:
        if SGData['SGSys'] == 'monoclinic':
            BNSsym = {}
            if SGData['SGUniq'] == 'a':
                BNSsym.update({'C_b':[0,.5,.0],'C_c':[0,0,.5]})
            elif SGData['SGUniq'] == 'b':
                BNSsym.update({'C_a':[.5,0,0],'C_c':[0,0,.5],'C_B':[.5,0.,.5]})
        elif SGData['SGSys'] == 'orthorhombic':
            BNSsym = {'C_a':[.5,0,0],'C_b':[0,.5,0],'C_c':[0,0,.5],
                'C_A':[0,.5,.5],'C_B':[.5,0,.5]}      
        elif SGData['SGSys'] == 'triclinic':
            BNSsym = {'C_a':[.5,0,0],'C_b':[0,.5,0],'C_c':[0,0,.5]}      
            
    elif 'I' in SGData['SGLatt']:
        if SGData['SGSys'] in ['monoclinic','orthorhombic','triclinic']:
            BNSsym = {'I_a':[.5,0,0],'I_b':[0,.5,0],'I_c':[0,0,.5]}
        elif SGData['SGSys'] == 'tetragonal':
            BNSsym = {'I_c':[0,0,.5]}
        elif SGData['SGSys'] == 'cubic':
            BNSsym = {} 
            
    elif 'F' in SGData['SGLatt']:
        if SGData['SGSys'] in ['monoclinic','orthorhombic','cubic','triclinic']:
            BNSsym = {'F_S':[.5,.5,.5]}
            
    elif 'R' in SGData['SGLatt']:
        BNSsym = {'R_I':[0,0,.5]}
        
    if SGData['SGGray']:
        for bns in BNSsym:
            BNSsym[bns].append(0.5)
            
    return UsymOp,OprFlg,BNSsym

def ApplyBNSlatt(SGData,BNSlatt):
    Tmat = np.eye(3)
    BNS = BNSlatt[0]
    A = np.array(BNSlatt[1])
    Laue = SGData['SGLaue']
    SGCen = SGData['SGCen']
    if '_a' in BNS:
        Tmat[0,0] = 2.0
    elif '_b' in BNS:
        Tmat[1,1] = 2.0
    elif '_c' in BNS:
        Tmat[2,2] = 2.0
    elif '_A' in BNS:
        Tmat[0,0] = 2.0
    elif '_B' in BNS:
        Tmat[1,1] = 2.0
    elif '_C' in BNS:
        Tmat[2,2] = 2.0
    elif '_I' in BNS:
        Tmat *= 2.0
        if 'R' in Laue:
            SGData['SGSpin'][-1] = -1
        else:
            SGData['SGSpin'].append(-1)
    elif '_S' in BNS:
        SGData['SGSpin'][-1] = -1
        SGData['SGSpin'] += [-1,-1,-1,]
        Tmat *= 2.0
    else:
        return Tmat
    if SGData.get('SGGray',False):
        SGData['SGSpin'].append(1)     #BNS centering are spin invrsion
    else:
        SGData['SGSpin'].append(-1)     #BNS centering in grey are not spin invrsion
    C = SGCen+A[:3]
    SGData['SGCen'] = np.vstack((SGCen,C))%1.
    return Tmat
        
def CheckSpin(isym,SGData):
    ''' Check for exceptions in spin rules
    '''
    if SGData['SGPtGrp'] in ['222','mm2','2mm','m2m']:      #only 2/3 can be red; not 1/3 or 3/3
        if SGData['SGSpin'][1]*SGData['SGSpin'][2]*SGData['SGSpin'][3] < 0:
            SGData['SGSpin'][(isym+1)%3+1] *= -1
        if SGData['SpGrp'][0] == 'F' and isym > 2:
            SGData['SGSpin'][(isym+1)%3+3] == 1
    elif SGData['SGPtGrp'] == 'mmm':
        if SGData['SpGrp'][0] == 'F' and isym > 2:
            SGData['SGSpin'][(isym+1)%3+3] == 1

def MagSGSym(SGData):       #needs to use SGPtGrp not SGLaue!
    SGLaue = SGData['SGLaue']
    if '1' not in SGData['GenSym']:        #patch for old gpx files
        SGData['GenSym'] = ['1',]+SGData['GenSym']
        SGData['SGSpin'] = [1,]+list(SGData['SGSpin'])
    if len(SGData['SGSpin'])<len(SGData['GenSym']):
        SGData['SGSpin'] = [1,]+list(SGData['SGSpin'])      #end patch
    GenSym = SGData['GenSym'][1:]       #skip identity
    SpnFlp = SGData['SGSpin']
    # print('SpnFlp',SpnFlp,' GenSym',GenSym)
    SGPtGrp = SGData['SGPtGrp']
    if len(SpnFlp) == 1:
        SGData['MagPtGp'] = SGPtGrp
        return SGData['SpGrp']
    magSym = SGData['SpGrp'].split()
    if SGLaue in ['-1',]:
        SGData['MagPtGp'] = SGPtGrp
        if SpnFlp[1] == -1:
            magSym[1] += "'"
            SGData['MagPtGp'] += "'"
    elif SGLaue in ['2/m','4/m','6/m']: #all ok
        Uniq = {'a':1,'b':2,'c':3,'':1}
        Id = [0,1]
        if len(magSym) > 2:
            Id = [0,Uniq[SGData['SGUniq']]]
        sym = magSym[Id[1]].split('/')
        Ptsym = SGLaue.split('/')
        if len(GenSym) == 3:
            for i in [0,1,2]:
                if SpnFlp[i+1] < 0:
                    sym[i] += "'"
                    Ptsym[i] += "'"
        else:
            for i in range(len(GenSym)):
                if SpnFlp[i+1] < 0:                      
                    sym[i] += "'"
                    Ptsym[i] += "'"
        SGData['MagPtGp'] = '/'.join(Ptsym)
        magSym[Id[1]] = '/'.join(sym)
    elif SGPtGrp in ['mmm','mm2','m2m','2mm','222']:
        SGData['MagPtGp'] = ''
        for i in [0,1,2]:
            SGData['MagPtGp'] += SGPtGrp[i]
            if SpnFlp[i+1] < 0:
                magSym[i+1] += "'"
                SGData['MagPtGp'] += "'"
    elif SGLaue == '6/mmm': #ok
        if '-' in SGPtGrp:
            magPtGp = list(SGPtGrp[1:])
            magPtGp[0] = '-'+magPtGp[0]
        else:
            magPtGp = list(SGPtGrp)
        if len(GenSym) == 2:
            for i in [0,1]:
                if SpnFlp[i+1] < 0:
                    magSym[i+2] += "'"
                    magPtGp[i+1] += "'"
            if SpnFlp[1]*SpnFlp[2] < 0:
                magSym[1] += "'"
                magPtGp[0] += "'"
        else:
            sym = magSym[1].split('/')
            Ptsym = ['6','m']
            magPtGp = ['','m','m']
            for i in [0,1,2]:
                if SpnFlp[i+1] < 0:
                    if i:
                        magSym[i+1] += "'"
                        magPtGp[i] += "'"
                    else:
                        sym[1] += "'"
                        Ptsym[0] += "'"
            if SpnFlp[2]*SpnFlp[3] < 0:
                sym[0] += "'"                    
                Ptsym[0] += "'"                    
            magSym[1] = '/'.join(sym)
            magPtGp[0] = '/'.join(Ptsym)
        SGData['MagPtGp'] = ''.join(magPtGp)
    elif SGLaue == '4/mmm':
        if '-' in SGPtGrp:
            magPtGp = list(SGPtGrp[1:])
            magPtGp[0] = '-'+magPtGp[0]
        else:
            magPtGp = list(SGPtGrp)
        if len(GenSym) == 2:
            for i in [0,1]:
                if SpnFlp[i+1] < 0:
                    magSym[i+2] += "'"
                    magPtGp[i+1] += "'"
            if SpnFlp[1]*SpnFlp[2] < 0:
                magSym[1] += "'"
                magPtGp[0] += "'"
        else:
            if '/' in magSym[1]:    #P 4/m m m, etc.
                sym = magSym[1].split('/')
                Ptsym = ['4','m']
                magPtGp = ['','m','m']
                for i in [0,1,2]:
                    if SpnFlp[i+1] < 0:
                        if i:
                            magSym[i+1] += "'"
                            magPtGp[i] += "'"
                        else:
                            sym[1] += "'"
                            Ptsym[1] += "'"
                if SpnFlp[2]*SpnFlp[3] < 0:
                    sym[0] += "'"                    
                    Ptsym[0] += "'"                    
                magSym[1] = '/'.join(sym)
                magPtGp[0] = '/'.join(Ptsym)
            else:
                for i in [0,1]:
                    if SpnFlp[i+1] < 0:
                        magSym[i+2] += "'"
                if SpnFlp[1]*SpnFlp[2] < 0:
                    magSym[1] += "'"
        SGData['MagPtGp'] = ''.join(magPtGp)
    elif SGLaue in ['3','3m1','31m']:   #ok 
        if '-' in SGPtGrp:
            Ptsym = list(SGPtGrp[1:])
            Ptsym[0] = '-'+Ptsym[0]
        else:
            Ptsym = list(SGPtGrp)
        if len(GenSym) == 1:    #all ok
            Id = 2
            if (len(magSym) == 4) and (magSym[2] == '1'):
                Id = 3
            if '3' in GenSym[0]:
                Id = 1
            magSym[Id].strip("'")
            if SpnFlp[1] < 0:
                magSym[Id] += "'"
                Ptsym[Id-1] += "'"
        elif len(GenSym) == 2:
            if 'R' in GenSym[1]:
                magSym[-1].strip("'")
                if SpnFlp[1] < 0:
                    magSym[-1] += "'"
                    Ptsym[-1] += "'"
            else:
                i,j = [1,2]
                if magSym[2] == '1':
                    i,j = [1,3]
                magSym[i].strip("'")
                Ptsym[i-1].strip("'")
                magSym[j].strip("'")
                Ptsym[j-1].strip("'")
                if SpnFlp[1:3] == [1,-1]:
                    magSym[i] += "'"
                    Ptsym[i-1] += "'"
                elif SpnFlp[1:3] == [-1,-1]:
                    magSym[j] += "'"
                    Ptsym[j-1] += "'"
                elif SpnFlp[1:3] == [-1,1]:
                    magSym[i] += "'"
                    Ptsym[i-1] += "'"
                    magSym[j] += "'"
                    Ptsym[j-1] += "'"
        elif len(GenSym):
            if 'c' not in magSym[2]:
                i,j = [1,2]
                magSym[i].strip("'")
                Ptsym[i-1].strip("'")
                magSym[j].strip("'")
                Ptsym[j-1].strip("'")
                if SpnFlp[1:3] == [1,-1]:
                    magSym[i] += "'"
                    Ptsym[i-1] += "'"
                elif SpnFlp[1:3] == [-1,-1]:
                    magSym[j] += "'"
                    Ptsym[j-1] += "'"
                elif SpnFlp[2:4] == [-1,1]:
                    magSym[i] += "'"
                    Ptsym[i-1] += "'"
                    magSym[j] += "'"
                    Ptsym[j-1] += "'"
        SGData['MagPtGp'] = ''.join(Ptsym)
    elif SGData['SGPtGrp'] == '23' and len(magSym):
        SGData['MagPtGp'] = '23'
    elif SGData['SGPtGrp'] == 'm3':
        SGData['MagPtGp'] = "m3"
        if SpnFlp[1] < 0:
            magSym[1] += "'"
            magSym[2] += "'"
            SGData['MagPtGp'] = "m'3'"
        if SpnFlp[1] < 0:
            if not 'm' in magSym[1]:    #only Ia3
                magSym[1].strip("'")
                SGData['MagPtGp'] = "m3'"
    elif SGData['SGPtGrp'] in ['432','-43m']:
        Ptsym = SGData['SGPtGrp'].split('3')
        if SpnFlp[1] < 0:
            magSym[1] += "'"
            Ptsym[0] += "'"
            magSym[3] += "'"
            Ptsym[1] += "'"
        SGData['MagPtGp'] = '3'.join(Ptsym)
    elif SGData['SGPtGrp'] == 'm3m':
        Ptsym = ['m','3','m']
        if SpnFlp[1:3] == [-1,1]:
            magSym[1] += "'"
            Ptsym[0] += "'"
            magSym[2] += "'"
            Ptsym[1] += "'"
        elif SpnFlp[1:3] == [1,-1]:
            magSym[3] += "'"
            Ptsym[2] += "'"
        elif SpnFlp[1:3] == [-1,-1]:
            magSym[1] += "'"
            Ptsym[0] += "'"
            magSym[2] += "'"
            Ptsym[1] += "'"
            magSym[3] += "'"
            Ptsym[2] += "'"
        SGData['MagPtGp'] = ''.join(Ptsym)
#    print SpnFlp
    magSym[0] = SGData.get('BNSlattsym',[SGData['SGLatt'],[0,0,0]])[0]
    return ' '.join(magSym)

def fixMono(SpGrp):
    'fixes b-unique monoclinics in e.g. P 1 2/1c 1 --> P 21/c '
    Flds = SpGrp.split()
    if len(Flds) == 4:
        if Flds[2] != '1':
            return '%s %s'%(Flds[0],Flds[2])
        else:
            return None
    else:
        return SpGrp

def Trans2Text(Trans):
    "from transformation matrix to text"
    cells = ['a','b','c']
    Text = ''
    for row in Trans:
        Fld = ''
        for i in [0,1,2]:
            if row[i]:
                if Fld and row[i] > 0.:
                    Fld += '+'
                Fld += '%3.1f'%(row[i])+cells[i]
        Text += Fld
        Text += ','
        Text = Text.replace('1.0','').replace('.0','').replace('0.5','1/2')
    return Text[:-1]

def getlattSym(Trans):
    Fives = {'ababc':'abc','bcbca':'cba','acacb':'acb','cabab':'cab','abcab':'acb'}
    transText = Trans2Text(Trans)
    lattSym = ''
    for fld in transText.split(','):
        if 'a' in fld: lattSym += 'a'
        if 'b' in fld: lattSym += 'b'
        if 'c' in fld: lattSym += 'c'
    if len(lattSym) != 3:
        lattSym = 'abc'
#        lattSym = Fives[lattSym]
    return lattSym

def Text2MT(mcifOpr,CIF=True):
    "From space group cif text returns matrix/translation"
    XYZ = {'x':[1,0,0],'+x':[1,0,0],'-x':[-1,0,0],'y':[0,1,0],'+y':[0,1,0],'-y':[0,-1,0],
           'z':[0,0,1],'+z':[0,0,1],'-z':[0,0,-1],'x-y':[1,-1,0],'-x+y':[-1,1,0],'y-x':[-1,1,0],
           '+x-y':[1,-1,0],'+y-x':[-1,1,0]}
    ops = mcifOpr.split(",")
    M = []
    T = []
    for op in ops[:3]:
        ip = len(op)
        if '/' in op:
            try:    #mcif format
                nP = op.count('+')
                opMT = op.split('+')
                T.append(eval(opMT[nP]))
                if nP == 2:
                    opMT[0] = '+'.join(opMT[0:2])
            except NameError:   #normal cif format
                ip = op.index('/')
                T.append(eval(op[:ip+2]))
                opMT = [op[ip+2:],'']
        else:
            opMT = [op,'']
            T.append(0.)
        M.append(XYZ[opMT[0].lower()])
    return np.array(M),np.array(T)
            
def MagText2MTS(mcifOpr,CIF=True):
    "From magnetic space group cif text returns matrix/translation + spin flip"
    XYZ = {'x':[1,0,0],'+x':[1,0,0],'-x':[-1,0,0],'y':[0,1,0],'+y':[0,1,0],'-y':[0,-1,0],
           'z':[0,0,1],'+z':[0,0,1],'-z':[0,0,-1],'x-y':[1,-1,0],'-x+y':[-1,1,0],'y-x':[-1,1,0],
           '+x-y':[1,-1,0],'+y-x':[-1,1,0]}
    ops = mcifOpr.split(",")
    M = []
    T = []
    for op in ops[:3]:
        ip = len(op)
        if '/' in op:
            try:    #mcif format
                nP = op.count('+')
                opMT = op.split('+')
                T.append(eval(opMT[nP]))
                if nP == 2:
                    opMT[0] = '+'.join(opMT[0:2])
            except NameError:   #normal cif format
                ip = op.index('/')
                T.append(eval(op[:ip+2]))
                opMT = [op[ip+2:],'']
        else:
            opMT = [op,'']
            T.append(0.)
        M.append(XYZ[opMT[0].lower()])
    spnflp = 1
    if '-1' in ops[3]:
        spnflp = -1
    return np.array(M),np.array(T),spnflp
            
def MagSSText2MTS(Opr,G2=False):
    "From magnetic super space group cif text returns matrix/translation + spin flip"
    XYZ = {'x1':[1,0,0,0],'-x1':[-1,0,0,0],
           'x2':[0,1,0,0],'-x2':[0,-1,0,0],
           'x3':[0,0,1,0],'-x3':[0,0,-1,0],
           'x4':[0,0,0,1],'-x4':[0,0,0,-1],
           'x1-x2':[1,-1,0,0],'-x1+x2':[-1,1,0,0],
           'x1-x4':[1,0,0,-1],'-x1+x4':[-1,0,0,1],
           'x2-x4':[0,1,0,-1],'-x2+x4':[0,-1,0,1],
           '-x1-x2+x4':[-1,-1,0,1],'x1+x2-x4':[1,1,0,-1]}
    if G2:
        XYZ = {'x':[1,0,0,0],'+x':[1,0,0,0],'-x':[-1,0,0,0],
               'y':[0,1,0,0],'+y':[0,1,0,0],'-y':[0,-1,0,0],
               'z':[0,0,1,0],'+z':[0,0,1,0],'-z':[0,0,-1,0],
               't':[0,0,0,1],'+t':[0,0,0,1],'-t':[0,0,0,-1],
               'x-y':[1,-1,0,0],'+x-y':[1,-1,0,0],'-x+y':[-1,1,0,0],
               'x-t':[1,0,0,-1],'+x-t':[1,0,0,-1],'-x+t':[-1,0,0,1],
               'y-t':[0,1,0,-1],'+y-t':[0,1,0,-1],'-y+t':[0,-1,0,1],
               'x+y-t':[1,1,0,-1],'+x+y-t':[1,1,0,-1],'-x-y+t':[-1,-1,0,1]}
       
    ops = Opr.split(",")
    M = []
    T = []
    for op in ops[:4]:
        ip = len(op)
        if '/' in op:
            ip = op.index('/')
            if G2:
                T.append(eval(op[:ip+2]))                
                M.append(XYZ[op[ip+2:]])
            else:
                T.append(eval(op[ip-2:]))
                M.append(XYZ[op[:ip-2]])
        else:
            T.append(0.)
            M.append(XYZ[op])
    spnflp = 1
    if len(ops) == 5:
        if '-1' in ops[4]:
            spnflp = -1
    return np.array(M),np.array(T),spnflp

def GetSGSpin(SGData,MSgSym):
    'get spin generators from magnetic space group symbol'

    mSgSym = MSgSym+' '
    mFlds = mSgSym.split()
    mSpn = SGData['SGSpin']
    SGPtGrp = SGData['SGPtGrp']
    SGLaue = SGData['SGLaue']
    GenSym = GetGenSym(SGData)[0][1:]       #skip identity
    # print('MSgSym',MSgSym,' GenSym',GenSym)
    if SGLaue in ['-1',]:
        if "'" in mFlds[1]:
            mSpn[1] = -1
    elif SGLaue in ['2/m','4/m','6/m']: #all ok
        Uniq = {'a':1,'b':2,'c':3,'':1}
        Id = [0,1]
        if len(mFlds) > 2:
            Id = [0,Uniq[SGData['SGUniq']]]
        sym = mFlds[Id[1]].split('/')
        for i in range(len(GenSym)):
            if "'" in sym[i]:
                mSpn[i+1] = -1
    elif SGPtGrp in ['mmm','mm2','m2m','2mm','222']:
        for i in [0,1,2]:
            if "'" in mFlds[i+1]:
                mSpn[i+1] = -1
    elif SGLaue in ['4/mmm','6/mmm']: #ok
        if len(GenSym) == 2: #422, -4m2, -42m, 4mm, 622, -6m2, -62m, 6mm
            for i in [0,1]:
                if "'" in mFlds[i+1]:
                    mSpn[i] = -1
            mSpn[2] = mSpn[0]*mSpn[1]
        else: #4/mmm, 6/mmm
            sym = mFlds[1].split('/')
            if "'" in sym[1]:
                mSpn[1] = -1
            for i in [1,2]:
                if "'" in mFlds[i+1]:
                    mSpn[i+1] = -1
    elif SGLaue in ['3','3m1','31m']:   #ok 
        if len(GenSym) == 1:    #all ok
            Id = 2
            if (len(mFlds) == 4) and (mFlds[2] == '1'):
                Id = 3
            if '3' in GenSym[0]:
                Id = 1
            if "'" in mFlds[Id]:
                mSpn[1] = -1
        elif len(GenSym) == 2:
            if 'R' in GenSym[1]:
                if "'" in mFlds[-1]:
                    mSpn[1] = -1
            else:
                i,j = [1,2]
                if mFlds[2] == '1':
                    i,j = [1,3]
                if "'" in mFlds[i]:
                    mSpn[1:3] = [1,-1] 
                elif "'" in mFlds[i]:
                    mSpn[1:3] = [-1,-1]
                elif "'" in mFlds[i] and "'" in mFlds[i]:
                    mSpn[1:3] = [-1,1]
        elif len(GenSym):
            if 'c' not in mFlds[2]:
                i,j = [1,2]
                if "'" in mFlds[i]:
                    mSpn[1:3] = [1,-1] 
                elif "'" in mFlds[i]:
                    mSpn[1:3] = [-1,-1]
                elif "'" in mFlds[i] and "'" in mFlds[i]:
                    mSpn[2:4] = [-1,1]
    elif SGData['SGPtGrp'] == 'm3':
        if "'" in mFlds[1] and "'" in mFlds[2]:
            mSpn[1]  = -1
    elif SGData['SGPtGrp'] in ['432','-43m']:
        if "'" in mFlds[1] and "'" in mFlds[3]:
            mSpn[1] = -1
    elif SGData['SGPtGrp'] == 'm3m':
        if "'" in mFlds[1] and "'" in mFlds[2]:
            mSpn[1:3] = [-1,1]
        elif "'" in mFlds[3]:
            mSpn[1:3] = [1,-1]
        elif "'" in mFlds[1] and "'" in mFlds[2] and "'" in mFlds[3]:
            mSpn[1:3] = [-1,-1]
    return mSpn

def GenMagOps(SGData):
    FlpSpn = SGData['SGSpin']
    Nsym = len(SGData['SGOps'])
    Ncv = len(SGData['SGCen'])
    sgOp = [M for M,T in SGData['SGOps']]
    detM = [nl.det(M) for M in sgOp]
    oprName = [GetOprPtrName(str(irtx)) for irtx in PackRot(SGData['SGOps'])]
    if SGData['SGInv'] and not SGData['SGFixed']:
        Nsym *= 2
        detM += [nl.det(-M) for M in sgOp]
        sgOp += [-M for M,T in SGData['SGOps']]
        oprName += [GetOprPtrName(str(-irtx)) for irtx in PackRot(SGData['SGOps'])]
    Nsyms = 0
    sgOps = []
    OprNames = []
    detMs = []
    for incv in range(Ncv):
        Nsyms += Nsym
        sgOps += sgOp
        detMs += detM
        OprNames += oprName
    if SGData['SGFixed']:
        SpnFlp = SGData['SpnFlp']
    else:
        SpnFlp = np.ones(Nsym,dtype=np.int32)
        GenFlg = SGData.get('GenFlg',[0])
        Ngen = len(SGData['SGGen'])
        Nfl = len(GenFlg)
        for ieqv in range(Nsym):
            for iunq in range(Nfl):
                if SGData['SGGen'][ieqv%Ngen] & GenFlg[iunq]:
                    try:
                        SpnFlp[ieqv] *= FlpSpn[iunq]
                    except IndexError:
                        print('index error: ',Nsym,ieqv,Nfl,iunq)
                        FlpSpn = FlpSpn+[1,]
                        SpnFlp[ieqv] *= FlpSpn[iunq]
        if SGData['SGLaue'] == '6/m': # treat special as algorithm above fails 
            if SGData['SGSpin'] == [1,-1,1]:
                SpnFlp = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
            elif SGData['SGSpin'] == [1,1,-1]:
                SpnFlp = [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1]
            elif SGData['SGSpin'] == [1,-1,-1]:
                SpnFlp = [1,-1,1,-1,1,-1,-1,1,-1,1,-1,1]
        for incv in range(Ncv):
            if incv:
                try:
                    SpnFlp = np.concatenate((SpnFlp,SpnFlp[:Nsym]*FlpSpn[Nfl+incv-1]))
                except IndexError:
                    FlpSpn = [1,]+FlpSpn
                    SpnFlp = np.concatenate((SpnFlp,SpnFlp[:Nsym]*FlpSpn[Nfl+incv-1]))
        if SGData['SGGray']:
           SpnFlp = np.concatenate((SpnFlp,-SpnFlp))
           detMs = 2*detMs                   
    MagMom = SpnFlp*np.array(detMs)      #duplicate for no. centerings
    SGData['MagMom'] = MagMom
    return OprNames,SpnFlp
    
def GetOpNum(Opr,SGData):
    Nops = len(SGData['SGOps'])
    opNum = abs(Opr)%100
    cent = abs(Opr)//100
    if Opr < 0 and not SGData['SGFixed']:
        opNum += Nops
    if SGData['SGInv'] and not SGData['SGFixed']:
        Nops *= 2
    opNum += cent*Nops
    return opNum
        
################################################################################
#### Superspace group codes
################################################################################
        
def SSpcGroup(SGData,SSymbol):
    """
    Determines supersymmetry information from superspace group name; currently only for (3+1) superlattices

    :param SGData: space group data structure as defined in SpcGroup above (see :ref:`SGData<SGData_table>`).
    :param SSymbol: superspace group symbol extension (string) defining modulation direction & generator info.
    :returns: (SSGError,SSGData)
    
       * SGError = 0 for no errors; >0 for errors (see SGErrors below for details)
       * SSGData - is a dict (see :ref:`Superspace Group object<SSGData_table>`) with entries:
       
             * 'SSpGrp': full superspace group symbol, accidental spaces removed; for display only
             * 'SSGCen': 4D cell centering vectors [0,0,0,0] at least
             * 'SSGOps': 4D symmetry operations as [M,T] so that M*x+T = x'

    """
            
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
            if mod:
                return [-SSGKl[i] if mod[i] in ['a','b','g'] else SSGKl[i] for i in range(3)]
            else:
                return [SSGKl[i] for i in range(3)]
        
    def extendSSGOps(SSGOps):
        for OpA in SSGOps:
            OpAtxt = SSMT2text(OpA)
            if 't' not in OpAtxt:
                continue
            for OpB in SSGOps:
                OpBtxt = SSMT2text(OpB)
                if 't' not in OpBtxt:
                    continue
                OpC = list(SGProd(OpB,OpA))
                OpC[1] %= 1.
                OpCtxt = SSMT2text(OpC)
#                print OpAtxt.replace(' ','')+' * '+OpBtxt.replace(' ','')+' = '+OpCtxt.replace(' ','')
                for k,OpD in enumerate(SSGOps):
                    OpDtxt = SSMT2text(OpD)
                    OpDtxt2 = ''
                    if SGData['SGGray']:                        
                        OpDtxt2 = SSMT2text([OpD[0],OpD[1]+np.array([0.,0.,0.,.5])])
#                    print '    ('+OpCtxt.replace(' ','')+' = ? '+OpDtxt.replace(' ','')+')'
                    if OpCtxt == OpDtxt:
                        continue
                    elif OpCtxt == OpDtxt2:
                        continue
                    elif OpCtxt.split(',')[:3] == OpDtxt.split(',')[:3]:
                        if 't' not in OpDtxt:
                            SSGOps[k] = OpC
#                            print k,'   new:',OpCtxt.replace(' ','')
                            break
                        else:
                            OpCtxt = OpCtxt.replace(' ','')
                            OpDtxt = OpDtxt.replace(' ','')
                            Txt = OpCtxt+' conflicts with '+OpDtxt
#                            print (Txt)
                            return False,Txt
        return True,SSGOps
        
    def findMod(modSym):
        for a in ['a','b','g']:
            if a in modSym:
                return a
                
    def genSSGOps():
        SSGOps = SSGData['SSGOps'][:]
        iFrac = {}
        for i,frac in enumerate(SSGData['modSymb']):
            if frac in ['1/2','1/3','1/4','1/6','1']:
                iFrac[i] = frac+'.'
#        print SGData['SpGrp']+SSymbol
#        print 'SSGKl',SSGKl,'genQ',genQ,'iFrac',iFrac,'modSymb',SSGData['modSymb']
# set identity & 1,-1; triclinic
        SSGOps[0][0][3,3] = 1.
## expand if centrosymmetric
#        if SGData['SGInv']:
#            SSGOps += [[-1*M,V] for M,V in SSGOps[:]]
# monoclinic - all done & all checked
        if SGData['SGPtGrp'] in ['2','m']:  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
            for i in iFrac:
                SSGOps[1][0][3,i] = -SSGKl[0]
        elif SGData['SGPtGrp'] == '2/m':    #OK
            SSGOps[1][0][3,3] = SSGKl[1]
            if 's' in gensym:
                SSGOps[1][1][3] = 0.5
            for i in iFrac:
                SSGOps[1][0][3,i] = SSGKl[0]
            
# orthorhombic - all OK not fully checked
        elif SGData['SGPtGrp'] in ['222','mm2','m2m','2mm']:    #OK
            if SGData['SGPtGrp'] == '222':
                OrOps = {'g':{0:[1,3],1:[2,3]},'a':{1:[1,2],2:[1,3]},'b':{2:[3,2],0:[1,2]}} #OK
            elif SGData['SGPtGrp'] == 'mm2':
                OrOps = {'g':{0:[1,3],1:[2,3]},'a':{1:[2,1],2:[3,1]},'b':{0:[1,2],2:[3,2]}} #OK
            elif SGData['SGPtGrp'] == 'm2m':
                OrOps = {'b':{0:[1,2],2:[3,2]},'g':{0:[1,3],1:[2,3]},'a':{1:[2,1],2:[3,1]}} #OK
            elif SGData['SGPtGrp'] == '2mm':
                OrOps = {'a':{1:[2,1],2:[3,1]},'b':{0:[1,2],2:[3,2]},'g':{0:[1,3],1:[2,3]}} #OK
            a = findMod(SSGData['modSymb'])
            OrFrac = OrOps[a]
            for j in iFrac:
                for i in OrFrac[j]:
                    SSGOps[i][0][3,j] = -2.*eval(iFrac[j])*SSGKl[i-1]
            for i in [0,1,2]:
                SSGOps[i+1][0][3,3] = SSGKl[i]
                SSGOps[i+1][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
                if not E:
                    return E,SSGOps
        elif SGData['SGPtGrp'] == 'mmm':    #OK
            OrOps = {'g':{0:[1,3],1:[2,3]},'a':{1:[2,1],2:[3,1]},'b':{0:[1,2],2:[3,2]}} 
            a = findMod(SSGData['modSymb'])
            if a == 'g':
                SSkl = [1,1,1]
            elif a == 'a':
                SSkl = [-1,1,-1]
            else:
                SSkl = [1,-1,-1]
            OrFrac = OrOps[a]
            for j in iFrac:
                for i in OrFrac[j]:
                    SSGOps[i][0][3,j] = -2.*eval(iFrac[j])*SSkl[i-1]
            for i in [0,1,2]:
                SSGOps[i+1][0][3,3] = SSkl[i]
                SSGOps[i+1][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
                if not E:
                    return E,SSGOps                
# tetragonal - all done & checked
        elif SGData['SGPtGrp'] == '4':  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -1
        elif SGData['SGPtGrp'] == '-4': #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = 1
        elif SGData['SGPtGrp'] in ['4/m',]: #OK
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -SSGKl[0]
            for i,j in enumerate([1,3]):
                SSGOps[j][0][3,3] = 1
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
                if not E:
                    return E,SSGOps
        elif SGData['SGPtGrp'] in ['422','4mm','-42m','-4m2',]: #OK
            iGens = [1,4,5]
            if SGData['SGPtGrp'] in ['4mm','-4m2',]:
                iGens = [1,6,7]
            for i,j in enumerate(iGens):
                if '1/2' in SSGData['modSymb'] and i < 2:
                    SSGOps[j][0][3,1] = SSGKl[i]
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    if 's' in gensym and j == 6:
                        SSGOps[j][1][3] = -genQ[i]
                    else:
                        SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
                if not E:
                    return E,SSGOps
        elif SGData['SGPtGrp'] in ['4/mmm',]:#OK
            if '1/2' in SSGData['modSymb']:
                SSGOps[1][0][3,1] = -SSGKl[0]
                SSGOps[6][0][3,1] = SSGKl[1]
                if modsym:
                   SSGOps[1][1][3]  = -genQ[3]
            for i,j in enumerate([1,2,6,7]):
                SSGOps[j][0][3,3] = 1
                SSGOps[j][1][3] = genQ[i]
                E,Result = extendSSGOps(SSGOps)
                if not E:
                    return E,Result
                else:
                    SSGOps = Result
                
# trigonal - all done & checked
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
                    SSGOps[j][0][3,3] = 1
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
                ids = [1,3]
            if '1/3' in SSGData['modSymb']:
                SSGOps[ids[0]][0][3,1] = -SSGKl[0]
            for i,j in enumerate(ids):
                SSGOps[j][0][3,3] = 1
                if genQ[i+1]:
                    SSGOps[j][1][3] = genQ[i+1]
                     
# hexagonal all done & checked
        elif SGData['SGPtGrp'] == '6':  #OK
            SSGOps[1][0][3,3] = SSGKl[0]
            SSGOps[1][1][3] = genQ[0]
        elif SGData['SGPtGrp'] == '-6': #OK
            SSGOps[1][0][3,3] = SSGKl[0]
        elif SGData['SGPtGrp'] in ['6/m',]: #OK
            SSGOps[1][0][3,3] = -SSGKl[1]
            SSGOps[1][1][3] = genQ[0]
            SSGOps[2][1][3] = genQ[1]
        elif SGData['SGPtGrp'] in ['622',]: #OK
            for i,j in enumerate([1,9,8]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = -genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
            
        elif SGData['SGPtGrp'] in ['6mm','-62m','-6m2',]: #OK
            for i,j in enumerate([1,6,7]):
                SSGOps[j][0][3,3] = SSGKl[i]
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['6/mmm',]: # OK
            for i,j in enumerate([1,2,10,11]):
                SSGOps[j][0][3,3] = 1
                if genQ[i]:
                    SSGOps[j][1][3] = genQ[i]
                E,SSGOps = extendSSGOps(SSGOps)
        elif SGData['SGPtGrp'] in ['1','-1']: #triclinic - done
            return True,SSGOps
        E,SSGOps = extendSSGOps(SSGOps)
        return E,SSGOps
        
    def specialGen(gensym,modsym):
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
        elif SGData['SGPtGrp'] in ['mmm',]:
            if 'g' in modsym:
                if sym == 's00':
                    gensym = 's0s'
                elif sym == '0s0':
                    gensym = '0ss'
            elif 'a' in modsym:
                if sym == '0s0':
                    gensym = 'ss0'
                elif sym == '00s':
                    gensym = 's0s'
            elif 'b' in modsym:
                if sym == '00s':
                    gensym = '0ss'
                elif sym == 's00':
                    gensym = 'ss0'
        return gensym
                            
    Fracs = {'1/2':0.5,'1/3':1./3,'1':1.0,'0':0.,'s':.5,'t':1./3,'q':.25,'h':-1./6,'a':0.,'b':0.,'g':0.}
    if SGData['SGLaue'] in ['m3','m3m']:
        return '(3+1) superlattices not defined for cubic space groups',None
    elif SGData['SGLaue'] in ['3R','3mR']:
        return '(3+1) superlattices not defined for rhombohedral settings - use hexagonal setting',None
    try:
        modsym,gensym = splitSSsym(SSymbol)
    except ValueError:
        return 'Error in superspace symbol '+SSymbol,None
    modQ = [Fracs[mod] for mod in modsym]
    SSGKl = SGData['SSGKl'][:]
    if SGData['SGLaue'] in ['2/m','mmm']:
        SSGKl = fixMonoOrtho()
    Ngen = len(gensym)
    if SGData.get('SGGray',False):
        Ngen -= 1
    if len(gensym) and Ngen != len(SSGKl):
        return 'Wrong number of items in generator symbol '+''.join(gensym),None
    gensym = specialGen(gensym[:Ngen],modsym)
    genQ = [Fracs[mod] for mod in gensym[:Ngen]]
    if not genQ:
        genQ = [0,0,0,0]
    SSgSpc = SGData['SpGrp']+SSymbol
    if SGData['SGGray']:
        SSgSpc = SSgSpc.replace('('," 1'(")
    SSGData = {'SSpGrp':SSgSpc,'modQ':modQ,'modSymb':modsym,'SSGKl':SSGKl}
    SSCen = np.zeros((len(SGData['SGCen']),4))
    for icen,cen in enumerate(SGData['SGCen']):
        SSCen[icen,0:3] = cen
    if 'BNSlattsym' in SGData and '_' in SGData['BNSlattsym'][0]:
        Ncen = len(SGData['SGCen'])
        for icen in range(Ncen//2,Ncen):
            SSCen[icen,3] = 0.5
    SSGData['SSGCen'] = SSCen%1.
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
        if DEBUG:
            print ('Super spacegroup operators for '+SSGData['SSpGrp'])
            for Op in Result:
                print (SSMT2text(Op).replace(' ',''))
            if SGData['SGInv']:                                 
                for Op in Result:
                    Op = [-Op[0],-Op[1]%1.]
                    print (SSMT2text(Op).replace(' ',''))                                 
        return None,SSGData
    else:
        return Result+'\nOperator conflict - incorrect superspace symbol',None
    
def SSChoice(SGData):
    '''
    Gets the unique set of possible super space groups for a given space group
    '''
    ptgpSS = {'1':['(abg)',],'-1':['(abg)',],
                   
        '2':['(a0g)','(a1/2g)','(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)'],
        'm':['(a0g)','(a1/2g)','(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)'],
        '2/m':['(a0g)','(a1/2g)','(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)'],
        
        '222':['(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)',
               '(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',],
        'mm2':['(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)',
               '(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',],
        'm2m':['(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)',
               '(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',],
        '2mm':['(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)',
               '(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',],
        'mmm':['(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)',
               '(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',],
               
        '4':['(00g)','(1/21/2g)'],'4mm':['(00g)','(1/21/2g)'],
        '4/m':['(00g)','(1/21/2g)'],
        '422':['(00g)','(1/21/2g)'],'-4m2':['(00g)','(1/21/2g)'],'-42m':['(00g)','(1/21/2g)'],
        '4/mmm':['(00g)','(1/21/2g)'],
        
        '3':['(00g)','(1/31/3g)'],'-3':['(00g)','(1/31/3g)'],
        '32':['(00g)'],'3m':['(00g)'],'-3m':['(00g)'],
        '321':['(00g)'],'3m1':['(00g)'],'-3m1':['(00g)'],
        '312':['(00g)','(1/31/3g)'],'31m':['(00g)','(1/31/3g)'],'-31m':['(00g)','(1/31/3g)'],
        
        '6':['(00g)',],'6/m':['(00g)',],'-62m':['(00g)',],'-6m2':['(00g)',],
        '622':['(00g)',],'6/mmm':['(00g)',],'6mm':['(00g)',],
        
        '23':['',],'m3':['',],'432':['',],'-43m':['',],'m3m':['',]}
            
    ptgpTS = {'1':['0',],'-1':['0',],
              
        '2':['0','s'],'m':['0','s'],
        '2/m':['00','0s','ss','s0'],
        
        '222':['000','s00','0s0','00s',],
        'mm2':['000','s00','0s0','00s','ss0','s0s','0ss','q00','0q0','00q','0qq','q0q','qq0'],
        'm2m':['000','s00','0s0','00s','ss0','s0s','0ss','q00','0q0','00q','0qq','q0q','qq0'],
        '2mm':['000','s00','0s0','00s','ss0','s0s','0ss','q00','0q0','00q','0qq','q0q','qq0'],
        'mmm':['000','s00','0s0','00s','ss0','s0s','0ss','q00','0q0','00q','0qq','q0q','qq0'],
        
        '4':['0','q','s'],'4mm':['000','q00','s00','s0s','ss0','0ss','qq0','qqs'],
        '4/m':['00','s0'],'-4m2':['000','0s0','0q0'],'-42m':['000','00s'],
        '422':['000','q00','s00','s0s','ss0','0ss','qq0','qqs','0q0'],
        '4/mmm':['0000','s0s0','00ss','s00s','ss00','0ss0','0s0s'],
        
        '3':['0','t'],'-3':['0','t'],
        '32':['00','t0'],'3m':['00','0s'],'-3m':['00','0s'],
        '321':['000','t00'],'3m1':['000','0s0'],'-3m1':['000','0s0'],
        '312':['000','t00'],'31m':['000','00s'],'-31m':['000','00s'],
        
        '6':['0','h','t','s'],
        '6/m':['00','s0'],'-62m':['000','00s'],'-6m2':['000','0s0'],
        '622':['000','h00','t00','s00',],'6mm':['000','ss0','s0s','0ss',],
        '6/mmm':['0000','s0s0','00ss','s00s','ss00','0ss0','0s0s'],
        
        '23':['',],'m3':['',],'432':['',],'-43m':['',],'m3m':['',]}
    
    ptgp = SGData['SGPtGrp']
    SSChoice = []
    for ax in ptgpSS[ptgp]:
        for sx in ptgpTS[ptgp]:
            SSChoice.append(ax+sx)
            if SGData['SGGray']: SSChoice[-1] += 's'
    ssChoice = []
    ssHash = []
    for item in SSChoice:
        E,SSG = SSpcGroup(SGData,item)
        if SSG:
            sshash = hash(str(SSGPrint(SGData,SSG)[1]))
            if sshash not in ssHash:
                ssHash.append(sshash)
                ssChoice.append(item)
    return ssChoice
           
def splitSSsym(SSymbol):
    '''
    Splits supersymmetry symbol into two lists of strings
    '''
    mssym = SSymbol.replace(' ','').split(')')
    if len(mssym) > 1:
        modsym,gensym = mssym
    else:
        modsym = mssym[0]
        gensym = ''
    modsym = modsym.replace(',','')
    if "1'" in modsym:
        gensym = gensym[:-1]
    modsym = modsym.replace("1'",'')
    if gensym in ['0','00','000','0000']:       #get rid of extraneous symbols
        gensym = ''
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
        
def SSGPrint(SGData,SSGData,AddInv=False):
    '''
    Print the output of SSpcGroup in a nicely formatted way. Used in SSpaceGroup

    :param SGData: space group data structure as defined in SpcGroup above.
    :param SSGData: from :func:`SSpcGroup`
    :returns:
        SSGText - list of strings with the superspace group details
        SGTable - list of strings for each of the operations
    '''
    nCen = len(SSGData['SSGCen'])
    Mult = nCen*len(SSGData['SSGOps'])*(int(SGData['SGInv'])+1)
    if SGData.get('SGFixed',False):
        Mult = len(SSGData['SSGCen'])*len(SSGData['SSGOps'])
    SSsymb = SSGData['SSpGrp']
    if 'BNSlattsym' in SGData and '_' in SGData['BNSlattsym'][0]:
        SSsymb = SGData['BNSlattsym'][0]+SSsymb[1:]
    SSGCen = list(SSGData['SSGCen'])
    if SGData.get('SGGray',False):
        if SGData.get('SGFixed',False): 
            Mult //= 2
        else:
            SSGCen += list(SSGData['SSGCen']+[0,0,0,0.5])
            SSGCen =  np.array(SSGCen)%1.
    else:
        if "1'" in SSsymb:  #leftover in nonmag phase in mcif file
            nCen //= 2
            Mult //= 2
            SSsymb = SSsymb.replace("1'",'')[:-1]
    SSGText = []
    SSGText.append(' Superspace Group: '+SSsymb)
    CentStr = 'centrosymmetric'
    if not SGData['SGInv']:
        CentStr = 'non'+CentStr
    if SGData['SGLatt'] in 'ABCIFR':
        SSGText.append(' The lattice is '+CentStr+' '+SGData['SGLatt']+'-centered '+SGData['SGSys'].lower())
    else:
        SSGText.append(' The superlattice is '+CentStr+' '+'primitive '+SGData['SGSys'].lower())
    SSGText.append(' The Laue symmetry is '+SGData['SGLaue'])
    SGptGp = SGData['SGPtGrp']
    if SGData['SGGray']:
        SGptGp += "1'"
    SSGText.append(' The superlattice point group is '+SGptGp+', '+''.join([str(i) for i in SSGData['SSGKl']]))
    SSGText.append(' The number of superspace group generators is '+str(len(SGData['SSGKl'])))
    SSGText.append(' Multiplicity of a general site is '+str(Mult))
    if SGData['SGUniq'] in ['a','b','c']:
        SSGText.append(' The unique monoclinic axis is '+SGData['SGUniq'])
    if SGData['SGInv']:
        SSGText.append(' The inversion center is located at 0,0,0')
    if SGData['SGPolax']:
        SSGText.append(' The location of the origin is arbitrary in '+SGData['SGPolax'])
    SSGText.append(' ')
    if len(SSGCen) > 1:
        SSGText.append(' The equivalent positions are:')
        SSGText.append(' ('+SSLatt2text(SSGCen)+')+\n')
    else:
        SSGText.append(' The equivalent positions are:\n')
    SSGTable = []
    for i,Opr in enumerate(SSGData['SSGOps']):
        SSGTable.append('(%2d) %s'%(i+1,SSMT2text(Opr)))
    if AddInv and SGData['SGInv']:
        for i,Opr in enumerate(SSGData['SSGOps']):
            IOpr = [-Opr[0],-Opr[1]]
            SSGTable.append('(%2d) %s'%(i+1+len(SSGData['SSGOps']),SSMT2text(IOpr)))        
    return SSGText,SSGTable
    
def SSGModCheck(Vec,modSymb,newMod=True):
    ''' Checks modulation vector compatibility with supersymmetry space group symbol. 
    if newMod: Superspace group symbol takes precidence & the vector will be modified accordingly
    '''
    Fracs = {'1/2':0.5,'1/3':1./3,'1':1.0,'0':0.,'a':0.,'b':0.,'g':0.}
    modQ = [Fracs[mod] for mod in modSymb]
    if newMod:
        newVec = Vec
        if not np.any(Vec):
            newVec = [0.1 if (vec == 0.0 and mod in ['a','b','g']) else vec for [vec,mod] in zip(Vec,modSymb)]
        return [Q if mod not in ['a','b','g'] and vec != Q else vec for [vec,mod,Q] in zip(newVec,modSymb,modQ)],  \
            [True if mod in ['a','b','g'] else False for mod in modSymb]
    else:
        return Vec,[True if mod in ['a','b','g'] else False for mod in modSymb]

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
    lattDir = {4:'1/3',6:'1/2',8:'2/3',0:'0'}
    for vec in SSGCen:
        lattTxt += ' '
        for item in vec:
            lattTxt += '%s,'%(lattDir[int(item*12)])
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
        print (SGErrors(E))
        return
    E,B = SSpcGroup(A,SSymbol)    
    if E > 0:
        print (E)
        return
    for l in SSGPrint(B):
        print (l)
        
def SGProd(OpA,OpB):
    '''
    Form space group operator product. OpA & OpB are [M,V] pairs; 
        both must be of same dimension (3 or 4). Returns [M,V] pair
    '''
    A,U = OpA
    B,V = OpB
    M = np.inner(B,A.T)
    W = np.inner(B,U)+V
    return M,W
        
def GetLittleGrpOps(SGData,vec):
    ''' Find rotation part of operators that leave vec unchanged
    
    :param SGData: space group data structure as defined in SpcGroup above.
    :param vec: a numpy array of fractional vector coordinates
    :returns: Little - list of operators [M,T] that form the little gropu
    '''
    Little = []
    Ops = SGData['SGOps'][:]
    if SGData['SGInv']:
        Ops += [[-M,-T] for [M,T] in Ops]
    for [M,T] in Ops:
        tvec = np.inner(M,vec)%1.
        if np.allclose(tvec,vec%1.):
            Little.append([M,T])
    return Little
        
def MoveToUnitCell(xyz):
    '''
    Translates a set of coordinates so that all values are >=0 and < 1 

    :param xyz: a list or numpy array of fractional coordinates
    :returns: XYZ - numpy array of new coordinates now 0 or greater and less than 1
    '''
    XYZ = np.array(xyz)%1.
    cell = np.asarray(np.rint(XYZ-xyz),dtype=np.int32)
    return XYZ,cell
        
def Opposite(XYZ,toler=0.0002):
    '''
    Gives opposite corner, edge or face of unit cell for position within tolerance. 
        Result may be just outside the cell within tolerance 

    :param XYZ: 0 >= np.array[x,y,z] > 1 as by MoveToUnitCell
    :param toler: unit cell fraction tolerance making opposite
    :returns:
        XYZ: dict of opposite positions; key=unit cell & always contains XYZ
    '''
    perm3 = [[1,1,1],[0,1,1],[1,0,1],[1,1,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0]]
    TB = np.where(abs(XYZ-1)<toler,-1,0)+np.where(abs(XYZ)<toler,1,0)
    perm = TB*perm3
    cperm = ['%d,%d,%d'%(i,j,k) for i,j,k in perm]
    D = dict(zip(cperm,perm))
    new = {}
    for key in D:
        new[key] = np.array(D[key])+np.array(XYZ)
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
    :return: [[XYZEquiv],Idup,[UijEquiv],spnflp]

      *  [XYZEquiv] is list of equivalent positions (XYZ is first entry)
      *  Idup = [-][C]SS where SS is the symmetry operator number (1-24), C (if not 0,0,0)
      * is centering operator number (1-4) and - is for inversion
        Cell = unit cell translations needed to put new positions inside cell
        [UijEquiv] - equivalent Uij; absent if no Uij given
      * +1/-1 for spin inversion of operator - empty if not magnetic
        
    '''
    XYZEquiv = []
    UijEquiv = []
    Idup = []
    Cell = []
    inv = int(SGData['SGInv']+1)
    icen = SGData['SGCen']
    if SGData.get('SGFixed',False):
        inv = 1
    SpnFlp = SGData.get('SpnFlp',[])
    spnflp = []
    X = np.array(XYZ)
    mj = 0
    for ic,cen in enumerate(icen):
        C = np.array(cen)
        for invers in range(inv):
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
                cell = np.zeros(3,dtype=np.int32)
                if Move:
                    newX,cell = MoveToUnitCell(XT)
                else:
                    newX = XT
                if All:
                    if np.allclose(newX,X,atol=0.0002):     #do we want %1. here?
                        idup = False
                else:
                    if True in [np.allclose(newX%1.,oldX%1.,atol=0.0002) for oldX in XYZEquiv]:
                        idup = False
                if All or idup:
                    XYZEquiv.append(newX)
                    Idup.append(idup)
                    Cell.append(cell)
                    if len(Uij):
                        UijEquiv.append(newUij)
                    if len(SpnFlp):
                        spnflp.append(SpnFlp[mj])
                    else:
                        spnflp.append(1)
                mj += 1
    if len(Uij):
        return zip(XYZEquiv,UijEquiv,Idup,Cell,spnflp)
    else:
        return zip(XYZEquiv,Idup,Cell,spnflp)
        
def GenHKL(HKL,SGData):
    ''' Generates all equivlent reflections including Friedel pairs
    :param HKL:  [h,k,l] must be integral values
    :param SGData: space group data obtained from SpcGroup
    :returns: array Uniq: equivalent reflections
    '''
    
    Ops = SGData['SGOps']
    OpM = np.array([op[0] for op in Ops])
    Uniq = np.inner(OpM,HKL)
    Uniq = list(Uniq)+list(-1*Uniq)
    return np.array(Uniq)

def GenHKLf(HKL,SGData):
    '''
    Uses old GSAS Fortran routine genhkl.for

    :param HKL:  [h,k,l] must be integral values for genhkl.for to work
    :param SGData: space group data obtained from SpcGroup
    :returns: iabsnt,mulp,Uniq,phi

     *   iabsnt = True if reflection is forbidden by symmetry
     *   mulp = reflection multiplicity including Friedel pairs
     *   Uniq = numpy array of equivalent hkl in descending order of h,k,l
     *   phi = phase offset for each equivalent h,k,l

    '''
    hklf = list(HKL)+[0,]       #could be numpy array!
    Ops = SGData['SGOps']
    OpM = np.array([op[0] for op in Ops],order='F')
    OpT = np.array([op[1] for op in Ops])
    Cen = np.array([cen for cen in SGData['SGCen']],order='F')
    
    import pyspg
    Nuniq,Uniq,iabsnt,mulp = pyspg.genhklpy(hklf,len(Ops),OpM,OpT,SGData['SGInv'],len(Cen),Cen)
    h,k,l,f = Uniq
    Uniq=np.array(list(zip(h[:Nuniq],k[:Nuniq],l[:Nuniq])))
    phi = f[:Nuniq]
    return iabsnt,mulp,Uniq,phi
    
def checkSSLaue(HKL,SGData,SSGData):
    #Laue check here - Toss HKL if outside unique Laue part
    h,k,l,m = HKL
    if SGData['SGLaue'] == '2/m':
        if SGData['SGUniq'] == 'a':
            if 'a' in SSGData['modSymb'] and h == 0 and m < 0:
                return False
            elif 'b' in SSGData['modSymb'] and k == 0 and l ==0 and m < 0:
                return False
            else:
                return True
        elif SGData['SGUniq'] == 'b':
            if 'b' in SSGData['modSymb'] and k == 0 and m < 0:
                return False
            elif 'a' in SSGData['modSymb'] and h == 0 and l ==0 and m < 0:
                return False
            else:
                return True
        elif SGData['SGUniq'] == 'c':
            if 'g' in SSGData['modSymb'] and l == 0 and m < 0:
                return False
            elif 'a' in SSGData['modSymb'] and h == 0 and k ==0 and m < 0:
                return False
            else:
                return True
    elif SGData['SGLaue'] == 'mmm':
        if 'a' in SSGData['modSymb']:
            if h == 0 and m < 0:
                return False
            else:
                return True
        elif 'b' in SSGData['modSymb']:
            if k == 0 and m < 0:
                return False
            else:
                return True
        elif 'g' in SSGData['modSymb']:
            if l == 0 and m < 0:
                return False
            else:
                return True
    else:   #tetragonal, trigonal, hexagonal (& triclinic?)
        if l == 0 and m < 0:
            return False
        else:
            return True
        
def checkHKLextc(HKL,SGData):
    '''
    Checks if reflection extinct - does not check centering

    :param HKL:  [h,k,l] 
    :param SGData: space group data obtained from SpcGroup
    :returns: True if extinct; False if allowed

    '''
    Ops = SGData['SGOps']
    OpM = np.array([op[0] for op in Ops])
    OpT = np.array([op[1] for op in Ops])
    HKLS = np.array([HKL,-HKL])     #Freidel's Law
    DHKL = np.reshape(np.inner(HKLS,OpM)-HKL,(-1,3))
    PHKL = np.reshape(np.inner(HKLS,OpT),(-1,))
    for dhkl,phkl in zip(DHKL,PHKL)[1:]:    #skip identity
        if dhkl.any():
            continue
        else:
            if phkl%1.:
                return True
    return False

def checkMagextc(HKL,SGData):
    '''
    Checks if reflection magnetically extinct; does fullcheck (centering, too)
    uses algorthm from Gallego, et al., J. Appl. Cryst. 45, 1236-1247 (2012)

    :param HKL:  [h,k,l] 
    :param SGData: space group data obtained from SpcGroup; must have magnetic symmetry SpnFlp data
    :returns: True if magnetically extinct; False if allowed (to match GenHKLf)

    '''
    Ops = SGData['SGOps']
    Ncen = len(SGData['SGCen'])
    OpM = np.array([op[0] for op in Ops])
    OpT = np.array([op[1] for op in Ops])
    if SGData['SGInv'] and not SGData['SGFixed']:
        OpM = np.vstack((OpM,-OpM))
        OpT = np.vstack((OpT,-OpT))%1.
    OpM = np.reshape(np.array(list(OpM)*Ncen),(-1,3,3))
    OpT = np.reshape(np.array([OpT+cen for cen in SGData['SGCen']]),(-1,3))
    Spn = SGData['SpnFlp'][:len(OpM)]
    Mag = np.array([nl.det(opm) for opm in OpM])*Spn
    DHKL = np.reshape(np.inner(HKL,OpM),(-1,3))
    PHKL = np.reshape(np.cos(2.0*np.pi*np.inner(HKL,OpT))*Mag,(-1,))[:,nxs,nxs]*OpM     #compute T(R,theta) eq(7)
    Ftest = np.random.rand(3)       #random magnetic moment
    Psum = np.zeros(3)
    nsum = 0.
    nA = 0
    for dhkl,phkl in zip(DHKL,PHKL):
        if not np.allclose(dhkl,HKL,atol=1.e-3):           #test for eq(5)
            continue
        else:
            nA += 1
            nsum += np.trace(phkl)          #eq(8)
            pterm = np.inner(Ftest,phkl)    #eq(9)
            Psum += pterm
    if nsum/nA > 1.:        #only need to look at nA=1 from eq(8)
        return False
    if np.allclose(Psum,np.zeros(3),atol=1.e-3):
        return True
    else:
        if np.abs(np.inner(HKL,Psum)) > 1.e-3:
            return True
        return False
    
def checkSSextc(HKL,SSGData):
    Ops = SSGData['SSGOps']
    OpM = np.array([op[0] for op in Ops])
    OpT = np.array([op[1] for op in Ops])
    HKLS = np.array([HKL,-HKL])     #Freidel's Law
    DHKL = np.reshape(np.inner(HKLS,OpM)-HKL,(-1,4))
    PHKL = np.reshape(np.inner(HKLS,OpT),(-1,))
    for dhkl,phkl in list(zip(DHKL,PHKL))[1:]:    #skip identity
        if dhkl.any():
            continue
        else:
            if phkl%1.:
                return False
    return True
    
################################################################################
#### Site symmetry tables
################################################################################
      
OprName = {
    '-6643':       ['-1',1],'6479' :    ['2(z)',2],'-6479':     ['m(z)',3],
    '6481' :     ['m(y)',4],'-6481':    ['2(y)',5],'6641' :     ['m(x)',6],
    '-6641':     ['2(x)',7],'6591' :  ['m(+-0)',8],'-6591':   ['2(+-0)',9],
    '6531' :  ['m(110)',10],'-6531': ['2(110)',11],'6537' :    ['4(z)',12],
    '-6537':   ['-4(z)',13],'975'  : ['3(111)',14],'6456' :       ['3',15],
    '-489' :  ['3(+--)',16],'483'  : ['3(-+-)',17],'-969' :  ['3(--+)',18],
    '819'  :  ['m(+0-)',19],'-819' : ['2(+0-)',20],'2431' :  ['m(0+-)',21],
    '-2431':  ['2(0+-)',22],'-657' :  ['m(xz)',23],'657'  :   ['2(xz)',24],
    '1943' :   ['-4(x)',25],'-1943':   ['4(x)',26],'-2429':   ['m(yz)',27],
    '2429' :   ['2(yz)',28],'639'  :  ['-4(y)',29],'-639' :    ['4(y)',30],
    '-6484':   ['2(010)',4],'6484' :  ['m(010)',5],'-6668':   ['2(100)',6],
    '6668' :   ['m(100)',7],'-6454': ['2(120)',18],'6454' :  ['m(120)',19],
    '-6638':  ['2(210)',20],'6638' : ['m(210)',21],   #search in SytSym ends at m(210)
    '2223' : ['3(+++)2',39],
    '6538' :   ['6(z)1',40],'-2169':['3(--+)2',41],'2151' : ['3(+--)2',42],
    '2205' :['-3(-+-)2',43],'-2205':[' (-+-)2',44],'489'  :['-3(+--)1',45],
    '801'  :   ['4(y)1',46],'1945' :  ['4(x)3',47],'-6585': ['-4(z)3 ',48],
    '6585' :   ['4(z)3',49],'6584' :  ['3(z)2',50],'6666' :  ['6(z)5 ',51],
    '6643' :       ['1',52],'-801' : ['-4(y)1',53],'-1945': ['-4(x)3 ',54],
    '-6666':  ['-6(z)5',55],'-6538': ['-6(z)1',56],'-2223':['-3(+++)2',57],
    '-975' :['-3(+++)1',58],'-6456': ['-3(z)1',59],'-483' :['-3(-+-)1',60],
    '969'  :['-3(--+)1',61],'-6584': ['-3(z)2',62],'2169' :['-3(--+)2',63],
    '-2151':['-3(+--)2',64],   }                               

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
    '13369369'  :'  mmm(y)','1927'      :'  mmm(z)','33554496'  :'  4(x)','16777280'  :' -4(x)',
    '50331745'  :'4/m(x)'  ,'169869394' :'422(x)','84934738'  :'-42m(x)','101711948' :'4mm(x)',
    '254804095' :'4/mmm(x)','536870928 ':'  4(y)','268435472' :' -4(y)','805306393' :'4/m(y)',
    '545783890' :'422(y)','272891986' :'-42m(y)','541327412' :'4mm(y)','818675839' :'4/mmm(y)',
    '2050'      :'  4(z)','4098'      :' -4(z)','6151'      :'4/m(z)','3410'      :'422(z)',
    '4818'      :'-42m(z)','2730'      :'4mm(z)','8191'      :'4/mmm(z)','8192'      :'  3(111)',
    '8193'      :' -3(111)','2629888'   :' 32(111)','1319040'   :' 3m(111)','3940737'   :'-3m(111)',
    '32768'     :'  3(+--)','32769'     :' -3(+--)','10519552'  :' 32(+--)','5276160'   :' 3m(+--)',
    '15762945'  :'-3m(+--)','65536'     :'  3(-+-)','65537'     :' -3(-+-)','134808576' :' 32(-+-)',
    '67437056'  :' 3m(-+-)','202180097' :'-3m(-+-)','131072'    :'  3(--+)','131073'    :' -3(--+)',
    '142737664' :' 32(--+)','71434368'  :' 3m(--+)','214040961' :'-3m(--+)','237650'    :'   23   ',
    '237695'    :'   m3   ','715894098' :'   432  ','358068946' :'  -43m  ','1073725439':'   m3m  ',
    '68157504'  :' mm2(d100)','4456464'   :' mm2(d010)','642'       :' mm2(d001)','153092172' :'-4m2(x)',
    '277348404' :'-4m2(y)','5418'      :'-4m2(z)','1075726335':'  6/mmm ','1074414420':'-6m2(100)',
    '1075070124':'-6m2(120)','1075069650':'   6mm  ','1074414890':'   622  ','1073758215':'   6/m  ',
    '1073758212':'   -6   ','1073758210':'    6   ','1073759865':'-3m(100)','1075724673':'-3m(120)',
    '1073758800':' 3m(100)','1075069056':' 3m(120)','1073759272':' 32(100)','1074413824':' 32(120)',
    '1073758209':'   -3   ','1073758208':'    3   ','1074135143':'mmm(100)','1075314719':'mmm(010)',
    '1073743751':'mmm(110)','1074004034':' mm2(z100)','1074790418':' mm2(z010)','1073742466':' mm2(z110)',
    '1074004004':'mm2(100)','1074790412':'mm2(010)','1073742980':'mm2(110)','1073872964':'mm2(120)',
    '1074266132':'mm2(210)','1073742596':'mm2(+-0)','1073872930':'222(100)','1074266122':'222(010)',
    '1073743106':'222(110)','1073741831':'2/m(001)','1073741921':'2/m(100)','1073741849':'2/m(010)',
    '1073743361':'2/m(110)','1074135041':'2/m(120)','1075314689':'2/m(210)','1073742209':'2/m(+-0)',
    '1073741828':' m(001) ','1073741888':' m(100) ','1073741840':' m(010) ','1073742336':' m(110) ',
    '1074003968':' m(120) ','1074790400':' m(210) ','1073741952':' m(+-0) ','1073741826':' 2(001) ',
    '1073741856':' 2(100) ','1073741832':' 2(010) ','1073742848':' 2(110) ','1073872896':' 2(120) ',
    '1074266112':' 2(210) ','1073742080':' 2(+-0) ','1073741825':'   -1   ',
    }

NXUPQsym = {
    '1'        :(28,29,28,28),'-1'       :( 1,29,28, 0),'2(x)'     :(12,18,12,25),'m(x)'     :(25,18,12,25),
    '2/m(x)'   :( 1,18, 0,-1),'2(y)'     :(13,17,13,24),'m(y)'     :(24,17,13,24),'2/m(y)'   :( 1,17, 0,-1),
    '2(z)'     :(14,16,14,23),'m(z)'     :(23,16,14,23),'2/m(z)'   :( 1,16, 0,-1),'2(yz)'    :(10,23,10,22),
    'm(yz)'    :(22,23,10,22),' 2/m(yz)' :( 1,23, 0,-1),'2(0+-)'   :(11,24,11,21),'m(0+-)'   :(21,24,11,21),
    '2/m(0+-)' :( 1,24, 0,-1),'2(xz)'    :( 8,21, 8,20),'m(xz)'    :(20,21, 8,20),'2/m(xz)'  :( 1,21, 0,-1),
    '2(+0-)'   :( 9,22, 9,19),'m(+0-)'   :(19,22, 9,19),'2/m(+0-)' :( 1,22, 0,-1),'2(xy)'    :( 6,19, 6,18),
    'm(xy)'    :(18,19, 6,18),' 2/m(xy)' :( 1,19, 0,-1),'2(+-0)'   :( 7,20, 7,17),'m(+-0)'   :(17,20, 7,17),
    'mm2(x)'  :(12,10, 0,-1),'mm2(y)'   :(13,10, 0,-1),'mm2(z)'   :(14,10, 0,-1),
    'mm2(yz)'  :(10,13, 0,-1),'mm2(0+-)' :(11,13, 0,-1),'mm2(xz)'  :( 8,12, 0,-1),'mm2(+0-)' :( 9,12, 0,-1),
    'mm2(xy)'  :( 6,11, 0,-1),'222'      :( 1,10, 0,-1),'222(x)'   :( 1,13, 0,-1),
    '222(y)'   :( 1,12, 0,-1),'222(z)'   :( 1,11, 0,-1),'mmm'      :( 1,10, 0,-1),'mmm(x)'   :( 1,13, 0,-1),
    'mmm(y)'   :( 1,12, 0,-1),'mmm(z)'   :( 1,11, 0,-1),'4(x)'     :(12, 4,12, 0),'-4(x)'    :( 1, 4,12, 0),
    '4/m(x)'   :( 1, 4,12,-1),'422(x)'   :( 1, 4, 0,-1),'-42m(x)'  :( 1, 4, 0,-1),'4mm(x)'   :(12, 4, 0,-1),
    '4/mmm(x)' :( 1, 4, 0,-1),'4(y)'     :(13, 3,13, 0),'-4(y)'    :( 1, 3,13, 0),'4/m(y)'   :( 1, 3,13,-1),
    '422(y)'   :( 1, 3, 0,-1),'-42m(y)'  :( 1, 3, 0,-1),'4mm(y)'   :(13, 3, 0,-1),'4/mmm(y)' :(1, 3, 0,-1,),
    '4(z)'     :(14, 2,14, 0),'-4(z)'    :( 1, 2,14, 0),'4/m(z)'   :( 1, 2,14,-1),'422(z)'   :( 1, 2, 0,-1),
    '-42m(z)'  :( 1, 2, 0,-1),'4mm(z)'   :(14, 2, 0,-1),'4/mmm(z)' :( 1, 2, 0,-1),'3(111)'   :( 2, 5, 2, 0),
    '-3(111)'  :( 1, 5, 2, 0),'32(111)'  :( 1, 5, 0, 2),'3m(111)'  :( 2, 5, 0, 2),'-3m(111)' :( 1, 5, 0,-1),
    '3(+--)'   :( 5, 8, 5, 0),'-3(+--)'  :( 1, 8, 5, 0),'32(+--)'  :( 1, 8, 0, 5),'3m(+--)'  :( 5, 8, 0, 5),
    '-3m(+--)' :( 1, 8, 0,-1),'3(-+-)'   :( 4, 7, 4, 0),'-3(-+-)'  :( 1, 7, 4, 0),'32(-+-)'  :( 1, 7, 0, 4),
    '3m(-+-)'  :( 4, 7, 0, 4),'-3m(-+-)' :( 1, 7, 0,-1),'3(--+)'   :( 3, 6, 3, 0),'-3(--+)'  :( 1, 6, 3, 0),
    '32(--+)'  :( 1, 6, 0, 3),'3m(--+)'  :( 3, 6, 0, 3),'-3m(--+)' :( 1, 6, 0,-1),'23'       :( 1, 1, 0, 0),
    'm3'       :( 1, 1, 0, 0),'432'      :( 1, 1, 0, 0),'-43m'     :( 1, 1, 0, 0),'m3m'      :( 1, 1, 0, 0),
    'mm2(d100)':(12,13, 0,-1),'mm2(d010)':(13,12, 0,-1),'mm2(d001)':(14,11, 0,-1),'-4m2(x)'  :( 1, 4, 0,-1),
    '-4m2(y)'  :( 1, 3, 0,-1),'-4m2(z)'  :( 1, 2, 0,-1),'6/mmm'    :( 1, 9, 0,-1),'-6m2(100)':( 1, 9, 0,-1),
    '-6m2(120)':( 1, 9, 0,-1),'6mm'      :(14, 9, 0,-1),'622'      :( 1, 9, 0,-1),'6/m'      :( 1, 9,14,-1),
    '-6'       :( 1, 9,14, 0),'6'        :(14, 9,14, 0),'-3m(100)' :( 1, 9, 0,-1),'-3m(120)' :( 1, 9, 0,-1),
    '3m(100)'  :(14, 9, 0,14),'3m(120)'  :(14, 9, 0,14),'32(100)'  :( 1, 9, 0,14),'32(120)'  :( 1, 9, 0,14),
    '-3'       :( 1, 9,14, 0),'3'        :(14, 9,14, 0),'mmm(100)' :( 1,14, 0,-1),'mmm(010)' :( 1,15, 0,-1),
    'mmm(110)' :( 1,11, 0,-1),'mm2(z100)':(14,14, 0,-1),'mm2(z010)':(14,15, 0,-1),'mm2(z110)':(14,11, 0,-1),
    'mm2(100)' :(12,14, 0,-1),'mm2(010)' :(13,15, 0,-1),'mm2(110)' :( 6,11, 0,-1),'mm2(120)' :(15,14, 0,-1),
    'mm2(210)' :(16,15, 0,-1),'mm2(+-0)' :( 7,11, 0,-1),'222(100)' :( 1,14, 0,-1),'222(010)' :( 1,15, 0,-1),
    '222(110)' :( 1,11, 0,-1),'2/m(001)' :( 1,16,14,-1),'2/m(100)' :( 1,25,12,-1),'2/m(010)' :( 1,28,13,-1),
    '2/m(110)' :( 1,19, 6,-1),'2/m(120)' :( 1,27,15,-1),'2/m(210)' :( 1,26,16,-1),'2/m(+-0)' :( 1,20,17,-1),
    'm(001)'   :(23,16,14,23),'m(100)'   :(26,25,12,26),'m(010)'   :(27,28,13,27),'m(110)'   :(18,19, 6,18),
    'm(120)'   :(24,27,15,24),'m(210)'   :(25,26,16,25),'m(+-0)'   :(17,20, 7,17),'2(001)'   :(14,16,14,23),
    '2(100)'   :(12,25,12,26),'2(010)'   :(13,28,13,27),'2(110)'   :( 6,19, 6,18),'2(120)'   :(15,27,15,24),
    '2(210)'   :(16,26,16,25),'2(+-0)'   :( 7,20, 7,17),'-1'       :( 1,29,28, 0)
    }
        
CSxinel = [[],      # 0th empty - indices are Fortran style
    [[0,0,0],[ 0.0, 0.0, 0.0],[0,0,0]],      #1  0  0  0
    [[1,1,1],[ 1.0, 1.0, 1.0],[1,1,1]],      #2  X  X  X
    [[1,1,1],[ 1.0, 1.0,-1.0],[1,1,1]],      #3  X  X -X
    [[1,1,1],[ 1.0,-1.0, 1.0],[1,1,1]],      #4  X -X  X
    [[1,1,1],[ 1.0,-1.0,-1.0],[1,1,1]],      #5 -X  X  X
    [[1,1,0],[ 1.0, 1.0, 0.0],[1,1,0]],      #6  X  X  0
    [[1,1,0],[ 1.0,-1.0, 0.0],[1,1,0]],      #7  X -X  0
    [[1,0,1],[ 1.0, 0.0, 1.0],[1,0,1]],      #8  X  0  X
    [[1,0,1],[ 1.0, 0.0,-1.0],[1,0,1]],      #9  X  0 -X
    [[0,1,1],[ 0.0, 1.0, 1.0],[0,2,2]],      #10  0  Y  Y
    [[0,1,1],[ 0.0, 1.0,-1.0],[0,2,2]],      #11 0  Y -Y
    [[1,0,0],[ 1.0, 0.0, 0.0],[1,0,0]],      #12  X  0  0
    [[0,1,0],[ 0.0, 1.0, 0.0],[0,2,0]],      #13  0  Y  0
    [[0,0,1],[ 0.0, 0.0, 1.0],[0,0,3]],      #14  0  0  Z
    [[1,1,0],[ 1.0, 2.0, 0.0],[1,1,0]],      #15  X 2X  0
    [[1,1,0],[ 2.0, 1.0, 0.0],[1,1,0]],      #16 2X  X  0
    [[1,1,2],[ 1.0, 1.0, 1.0],[1,1,3]],      #17  X  X  Z
    [[1,1,2],[ 1.0,-1.0, 1.0],[1,1,3]],      #18  X -X  Z
    [[1,2,1],[ 1.0, 1.0, 1.0],[1,2,1]],      #19  X  Y  X
    [[1,2,1],[ 1.0, 1.0,-1.0],[1,2,1]],      #20  X  Y -X
    [[1,2,2],[ 1.0, 1.0, 1.0],[1,2,2]],      #21  X  Y  Y
    [[1,2,2],[ 1.0, 1.0,-1.0],[1,2,2]],      #22  X  Y -Y
    [[1,2,0],[ 1.0, 1.0, 0.0],[1,2,0]],      #23  X  Y  0
    [[1,0,2],[ 1.0, 0.0, 1.0],[1,0,3]],      #24  X  0  Z
    [[0,1,2],[ 0.0, 1.0, 1.0],[0,2,3]],      #25  0  Y  Z
    [[1,1,2],[ 1.0, 2.0, 1.0],[1,1,3]],      #26  X 2X  Z
    [[1,1,2],[ 2.0, 1.0, 1.0],[1,1,3]],      #27 2X  X  Z
    [[1,2,3],[ 1.0, 1.0, 1.0],[1,2,3]],      #28  X  Y  Z
    ]

CSuinel = [[],      # 0th empty - indices are Fortran style
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
    [[1,2,3,2,4,4],[ 1.0, 1.0, 1.0, 0.5, 1.0, 2.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #25  A  B  C B/2 F/2 F
    [[1,2,3,1,0,4],[ 1.0, 1.0, 1.0, 0.5, 0.0, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #26  A  B  C A/2  0  F
    [[1,2,3,2,4,0],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.0],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #27  A  B  C B/2  E  0
    [[1,2,3,1,4,4],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.5],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #28  A  B  C A/2  E E/2
    [[1,2,3,4,5,6],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,1,1,1,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #29  A  B  C  D  E   F
    ]
    
################################################################################
#### Site symmetry routines
################################################################################
    
def GetOprPtrName(key):
    'Needs a doc string'
    try:
        oprName = OprName[key][0]
    except KeyError:
        return key
    return oprName.replace('(','').replace(')','')

def GetOprPtrNumber(key):
    'Needs a doc string'
    try:
        return OprName[key][1]
    except KeyError:
        return key

def GetOprName(key):
    'Needs a doc string'
    return OprName[key][0]

def GetKNsym(key):
    'Needs a doc string'
    try:
        return KNsym[key].strip()
    except KeyError:
        return 'sp'

def GetNXUPQsym(siteSym):
    '''        
    The codes XUPQ are for lookup of symmetry constraints for position(X), thermal parm(U) & magnetic moments (P & Q) 
    '''
    return NXUPQsym[siteSym]

def GetCSxinel(siteSym):  
    "returns Xyz terms, multipliers, GUI flags"
    indx = GetNXUPQsym(siteSym.strip())
    return CSxinel[indx[0]]
    
def GetCSuinel(siteSym):
    "returns Uij terms, multipliers, GUI flags & Uiso2Uij multipliers"
    indx = GetNXUPQsym(siteSym.strip())
    return CSuinel[indx[1]]
    
def GetCSpqinel(SpnFlp,dupDir):  
    "returns Mxyz terms, multipliers, GUI flags"
    CSI = [[1,2,3],[1.0,1.0,1.0]]
    for sopr in dupDir:
#        print (sopr,dupDir[sopr])
        opr = sopr.replace("'",'')
        indx = GetNXUPQsym(opr)
        if SpnFlp[dupDir[sopr]] > 0:
            csi = CSxinel[indx[2]]  #P
        else:
            csi = CSxinel[indx[3]]  #Q
#        print(opr,indx,csi,CSI)
        if not len(csi):
            return [[0,0,0],[0.,0.,0.]]
        for kcs in [0,1,2]:
            if csi[0][kcs] == 0 and CSI[0][kcs] != 0:
                jcs = CSI[0][kcs]
                for ics in [0,1,2]:
                    if CSI[0][ics] == jcs:
                        CSI[0][ics] = 0
                        CSI[1][ics] = 0.
                    elif CSI[0][ics] > jcs:
                        CSI[0][ics] = CSI[0][ics]-1
            elif (CSI[0][kcs] == csi[0][kcs]) and (CSI[1][kcs] != csi[1][kcs]):
                CSI[1][kcs] = csi[1][kcs]
            elif CSI[0][kcs] >= csi[0][kcs]:
                CSI[0][kcs] = min(CSI[0][kcs],csi[0][kcs])
                if CSI[1][kcs] != csi[1][kcs]:
                    if CSI[1][kcs] == 1.:
                        CSI[1][kcs] = csi[1][kcs]
#        print(CSI)
    return CSI
    
def getTauT(tau,sop,ssop,XYZ,wave=np.zeros(3)):
    phase = np.sum(XYZ*wave)
    ssopinv = nl.inv(ssop[0])
    mst = ssopinv[3][:3]
    epsinv = ssopinv[3][3]
    sdet = nl.det(sop[0])
    ssdet = nl.det(ssop[0])
    dtau = mst*(XYZ-sop[1])-epsinv*ssop[1][3]
    dT = 1.0
    if np.any(dtau%.5):
        sumdtau = np.sum(dtau%.5)
        dT = 0.
        if np.abs(sumdtau-.5) > 1.e-4:
            dT = np.tan(np.pi*sumdtau)
    tauT = np.inner(mst,XYZ-sop[1])+epsinv*(tau-ssop[1][3]+phase)
    return sdet,ssdet,dtau,dT,tauT
    
def OpsfromStringOps(A,SGData,SSGData):
    SGOps = SGData['SGOps']
    SSGOps = SSGData['SSGOps']
    Ax = A.split('+')
    Ax[0] = int(Ax[0])
    iC = 1
    if Ax[0] < 0:
        iC = -1
    iAx = abs(Ax[0])
    nA = iAx%100-1
    nC = iAx//100
    unit = [0,0,0]
    if len(Ax) > 1:
        unit = eval('['+Ax[1]+']')
    return SGOps[nA],SSGOps[nA],iC,SGData['SGCen'][nC],unit
    
def GetSSfxuinel(waveType,Stype,nH,XYZ,SGData,SSGData,debug=False):
    
    def orderParms(CSI):
        parms = [0,]
        for csi in CSI:
            for i in [0,1,2]:
                if csi[i] not in parms:
                    parms.append(csi[i])
        for csi in CSI:
            for i in [0,1,2]:
                csi[i] = parms.index(csi[i])
        return CSI
        
    def fracCrenel(tau,Toff,Twid):
        Tau = (tau-Toff[:,nxs])%1.
        A = np.where(Tau<Twid[:,nxs],1.,0.)
        return A
        
    def fracFourier(tau,nH,fsin,fcos):
        SA = np.sin(2.*nH*np.pi*tau)
        CB = np.cos(2.*nH*np.pi*tau)
        A = SA[nxs,nxs,:]*fsin[:,:,nxs]
        B = CB[nxs,nxs,:]*fcos[:,:,nxs]
        return A+B
        
    def posFourier(tau,nH,psin,pcos):
        SA = np.sin(2*nH*np.pi*tau)
        CB = np.cos(2*nH*np.pi*tau)
        A = SA[nxs,nxs,:]*psin[:,:,nxs]
        B = CB[nxs,nxs,:]*pcos[:,:,nxs]
        return A+B    

    def posZigZag(tau,Tmm,XYZmax):
        DT = Tmm[1]-Tmm[0]
        slopeUp = 2.*XYZmax/DT
        slopeDn = 2.*XYZmax/(1.-DT)
        A = np.array([np.where(0. < t-(Tmm[0])%1. <= DT,-XYZmax+slopeUp*((t-Tmm[0])%1.),XYZmax-slopeDn*((t-Tmm[1])%1.)) for t in tau])
        return A

    def posBlock(tau,Tmm,XYZmax):
        A = np.array([np.where(Tmm[0] < t <= Tmm[1],XYZmax,-XYZmax) for t in tau])
        return A
        
    def DoFrac():
        delt2 = np.eye(2)*0.001
        dF = fracFourier(tau,nH,delt2[:1],delt2[1:]).squeeze()
        dFTP = []
        if siteSym == '1':
            CSI = [[1,0],[2,0]],2*[[1.,0.],]
        elif siteSym == '-1':
            CSI = [[1,0],[0,0]],2*[[1.,0.],]
        else:
            FSC = np.ones(2,dtype='i')
            CSI = [np.zeros((2),dtype='i'),np.zeros(2)]
            if 'Crenel' in waveType:
                dF = np.zeros_like(tau)
            else:
                dF = fracFourier(tau,nH,delt2[:1],delt2[1:]).squeeze()
            dFT = np.zeros_like(dF)
            dFTP = []
            for i in SdIndx:
                sop = Sop[i]
                ssop = SSop[i]            
                sdet,ssdet,dtau,dT,tauT = getTauT(tau,sop,ssop,XYZ)
                fsc = np.ones(2,dtype='i')
                if 'Crenel' in waveType:
                    dFT = np.zeros_like(tau)
                    fsc = [1,1]
                else:   #Fourier
                    dFT = fracFourier(tauT,nH,delt2[:1],delt2[1:]).squeeze()
                    dFT = nl.det(sop[0])*dFT
                    dFT = dFT[:,np.argsort(tauT)]
                    dFT[0] *= ssdet
                    dFT[1] *= sdet
                    dFTP.append(dFT)
                
                    if np.any(dtau%.5) and ('1/2' in SSGData['modSymb'] or '1' in SSGData['modSymb']):
                        fsc = [1,1]
                        if dT:
                            CSI = [[[1,0],[1,0]],[[1.,0.],[1/dT,0.]]]
                        else:
                            CSI = [[[1,0],[0,0]],[[1.,0.],[0.,0.]]]
                        FSC = np.zeros(2,dtype='i')
                        return CSI,dF,dFTP
                    else:
                        for i in range(2):
                            if np.allclose(dF[i,:],dFT[i,:],atol=1.e-6):
                                fsc[i] = 1
                            else:
                                fsc[i] = 0
                        FSC &= fsc
                        if debug: print (SSMT2text(ssop).replace(' ',''),sdet,ssdet,epsinv,fsc)
            n = -1
            for i,F in enumerate(FSC):
                if F:
                    n += 1
                    CSI[0][i] = n+1
                    CSI[1][i] = 1.0
            
        return CSI,dF,dFTP
        
    def DoXYZ():
        delt5 = np.ones(5)*0.001
        delt6 = np.eye(6)*0.001
        if 'Fourier' in waveType:
            dX = posFourier(tau,nH,delt6[:3],delt6[3:]) #+np.array(XYZ)[:,nxs,nxs]
              #3x6x12 modulated position array (X,Spos,tau)& force positive
        elif waveType in ['ZigZag','Block']:
            if waveType == 'ZigZag':
                dX = posZigZag(tau,delt5[:2],delt5[2:])
            else:
                dX = posBlock(tau,delt5[:2],delt5[2:])
        dXTP = []
        if siteSym == '1':
            CSI = [[1,0,0],[2,0,0],[3,0,0], [4,0,0],[5,0,0],[6,0,0]],6*[[1.,0.,0.],]
        elif siteSym == '-1':
            CSI = [[1,0,0],[2,0,0],[3,0,0], [0,0,0],[0,0,0],[0,0,0]],3*[[1.,0.,0.],]+3*[[0.,0.,0.],]
        else:
            if 'Fourier' in waveType:
                CSI = [np.zeros((6,3),dtype='i'),np.zeros((6,3))]
            elif waveType in ['ZigZag','Block']:
                CSI = [np.array([[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0]]),
                    np.array([[1.0,.0,.0],[1.0,.0,.0],[1.0,.0,.0],[1.0,.0,.0],[1.0,.0,.0]])]
            XSC = np.ones(6,dtype='i')
            dXTP = []
            for i in SdIndx:
                sop = Sop[i]
                ssop = SSop[i]
                sdet,ssdet,dtau,dT,tauT = getTauT(tau,sop,ssop,XYZ)
                xsc = np.ones(6,dtype='i')
                if 'Fourier' in waveType:
                    dXT = posFourier(np.sort(tauT),nH,delt6[:3],delt6[3:])   #+np.array(XYZ)[:,nxs,nxs]
                elif waveType == 'ZigZag':
                    dXT = posZigZag(tauT,delt5[:2],delt5[2:])+np.array(XYZ)[:,nxs,nxs]
                elif waveType == 'Block':
                    dXT = posBlock(tauT,delt5[:2],delt5[2:])+np.array(XYZ)[:,nxs,nxs]
                dXT = np.inner(sop[0],dXT.T)    # X modulations array(3x6x49) -> array(3x49x6)
                dXT = np.swapaxes(dXT,1,2)      # back to array(3x6x49)
                dXT[:,:3,:] *= (ssdet*sdet)            # modify the sin component
                dXTP.append(dXT)
                if waveType == 'Fourier':
                    for i in range(3):
                        if not np.allclose(dX[i,i,:],dXT[i,i,:]):
                            xsc[i] = 0
                        if not np.allclose(dX[i,i+3,:],dXT[i,i+3,:]):
                            xsc[i+3] = 0
                    if np.any(dtau%.5) and ('1/2' in SSGData['modSymb'] or '1' in SSGData['modSymb']):
                        xsc[3:6] = 0
                        CSI = [[[1,0,0],[2,0,0],[3,0,0], [1,0,0],[2,0,0],[3,0,0]],
                            [[1.,0.,0.],[1.,0.,0.],[1.,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]]
                        if dT:
                            if '(x)' in siteSym:
                                CSI[1][3:] = [1./dT,0.,0.],[-dT,0.,0.],[-dT,0.,0.]
                                if 'm' in siteSym and len(SdIndx) == 1:
                                    CSI[1][3:] = [-dT,0.,0.],[1./dT,0.,0.],[1./dT,0.,0.]
                            elif '(y)' in siteSym:
                                CSI[1][3:] = [-dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]
                                if 'm' in siteSym and len(SdIndx) == 1:
                                    CSI[1][3:] = [1./dT,0.,0.],[-dT,0.,0.],[1./dT,0.,0.]
                            elif '(z)' in siteSym:
                                CSI[1][3:] = [-dT,0.,0.],[-dT,0.,0.],[1./dT,0.,0.]
                                if 'm' in siteSym and len(SdIndx) == 1:
                                    CSI[1][3:] = [1./dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]
                        else:
                            CSI[1][3:] = [0.,0.,0.],[0.,0.,0.],[0.,0.,0.]
                    if '4/mmm' in laue:
                        if np.any(dtau%.5) and '1/2' in SSGData['modSymb']:
                            if '(xy)' in siteSym:
                                CSI[0] = [[1,0,0],[1,0,0],[2,0,0], [1,0,0],[1,0,0],[2,0,0]]
                                if dT:
                                    CSI[1][3:] = [[1./dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]]
                                else:
                                    CSI[1][3:] = [0.,0.,0.],[0.,0.,0.],[0.,0.,0.]
                        if '(xy)' in siteSym or '(+-0)' in siteSym:
                            mul = 1
                            if '(+-0)' in siteSym:
                                mul = -1
                            if np.allclose(dX[0,0,:],dXT[1,0,:]):
                                CSI[0][3:5] = [[11,0,0],[11,0,0]]
                                CSI[1][3:5] = [[1.,0,0],[mul,0,0]]
                                xsc[3:5] = 0
                            if np.allclose(dX[0,3,:],dXT[0,4,:]):
                                CSI[0][:2] = [[12,0,0],[12,0,0]]
                                CSI[1][:2] = [[1.,0,0],[mul,0,0]]
                                xsc[:2] = 0
                else:
                    for i in range(3):
                        if not np.allclose(dX[:,i],dXT[i,:,i]):
                            xsc[i] = 0
                XSC &= xsc
                if debug: print (SSMT2text(ssop).replace(' ',''),sdet,ssdet,epsinv,xsc)
            if waveType == 'Fourier':
                n = -1
                if debug: print (XSC)
                for i,X in enumerate(XSC):
                    if X:
                        n += 1
                        CSI[0][i][0] = n+1
                        CSI[1][i][0] = 1.0
            
        return list(CSI),dX,dXTP
        
    def DoUij():
        delt12 = np.eye(12)*0.0001
        dU = posFourier(tau,nH,delt12[:6],delt12[6:])                  #Uij modulations - 6x12x12 array
        dUTP = []
        if siteSym == '1':
            CSI = [[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0], 
                [7,0,0],[8,0,0],[9,0,0],[10,0,0],[11,0,0],[12,0,0]],12*[[1.,0.,0.],]
        elif siteSym == '-1':
            CSI = 6*[[0,0,0],]+[[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]],   \
                6*[[0.,0.,0.],]+[[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]
        else:
            CSI = [np.zeros((12,3),dtype='i'),np.zeros((12,3))]
            USC = np.ones(12,dtype='i')
            dUTP = []
            dtau = 0.
            for i in SdIndx:
                sop = Sop[i]
                ssop = SSop[i]
                sdet,ssdet,dtau,dT,tauT = getTauT(tau,sop,ssop,XYZ)
                usc = np.ones(12,dtype='i')
                dUT = posFourier(tauT,nH,delt12[:6],delt12[6:])                  #Uij modulations - 6x12x49 array
                dUijT = np.rollaxis(np.rollaxis(np.array(Uij2U(dUT)),3),3)    #convert dUT to 12x49x3x3 
                dUijT = np.rollaxis(np.inner(np.inner(sop[0],dUijT),sop[0].T),3) #transform by sop - 3x3x12x49
                dUT = np.array(U2Uij(dUijT))    #convert to 6x12x49
                dUT = dUT[:,:,np.argsort(tauT)]
                dUT[:,:6,:] *=(ssdet*sdet)
                dUTP.append(dUT)
                if np.any(dtau%.5) and ('1/2' in SSGData['modSymb'] or '1' in SSGData['modSymb']):
                    if dT:
                        CSI = [[[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0], 
                        [1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]],
                        [[1.,0.,0.],[1.,0.,0.],[1.,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.],
                        [1./dT,0.,0.],[1./dT,0.,0.],[1./dT,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]]
                    else:
                        CSI = [[[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0], 
                        [1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]],
                        [[1.,0.,0.],[1.,0.,0.],[1.,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.],
                        [0.,0.,0.],[0.,0.,0.],[0.,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]]
                    if 'mm2(x)' in siteSym and dT:
                        CSI[1][9:] = [0.,0.,0.],[-dT,0.,0.],[0.,0.,0.]
                        USC = [1,1,1,0,1,0,1,1,1,0,1,0]
                    elif '(xy)' in siteSym and dT:
                        CSI[0] = [[1,0,0],[1,0,0],[2,0,0],[3,0,0],[4,0,0],[4,0,0],
                            [1,0,0],[1,0,0],[2,0,0],[3,0,0],[4,0,0],[4,0,0]]
                        CSI[1][9:] = [[1./dT,0.,0.],[-dT,0.,0.],[-dT,0.,0.]]
                        USC = [1,1,1,1,1,1,1,1,1,1,1,1]                              
                    elif '(x)' in siteSym and dT:
                        CSI[1][9:] = [-dT,0.,0.],[-dT,0.,0.],[1./dT,0.,0.]
                    elif '(y)' in siteSym and dT:
                        CSI[1][9:] = [-dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]
                    elif '(z)' in siteSym and dT:
                        CSI[1][9:] = [1./dT,0.,0.],[-dT,0.,0.],[-dT,0.,0.]
                    for i in range(6):
                        if not USC[i]:
                            CSI[0][i] = [0,0,0]
                            CSI[1][i] = [0.,0.,0.]
                            CSI[0][i+6] = [0,0,0]
                            CSI[1][i+6] = [0.,0.,0.]
                else:                        
                    for i in range(6):
                        if not np.allclose(dU[i,i,:],dUT[i,i,:]):  #sin part
                            usc[i] = 0
                        if not np.allclose(dU[i,i+6,:],dUT[i,i+6,:]):   #cos part
                            usc[i+6] = 0
                    if np.any(dUT[1,0,:]):
                        if '4/m' in siteSym:
                            CSI[0][6:8] = [[12,0,0],[12,0,0]]
                            if ssop[1][3]:
                                CSI[1][6:8] = [[1.,0.,0.],[-1.,0.,0.]]
                                usc[9] = 1
                            else:
                                CSI[1][6:8] = [[1.,0.,0.],[1.,0.,0.]]
                                usc[9] = 0
                        elif '4' in siteSym:
                            CSI[0][6:8] = [[12,0,0],[12,0,0]]
                            CSI[0][:2] = [[11,0,0],[11,0,0]]
                            if ssop[1][3]:
                                CSI[1][:2] = [[1.,0.,0.],[-1.,0.,0.]]
                                CSI[1][6:8] = [[1.,0.,0.],[-1.,0.,0.]]
                                usc[2] = 0
                                usc[8] = 0
                                usc[3] = 1
                                usc[9] = 1
                            else:
                                CSI[1][:2] = [[1.,0.,0.],[1.,0.,0.]]
                                CSI[1][6:8] = [[1.,0.,0.],[1.,0.,0.]]
                                usc[2] = 1
                                usc[8] = 1
                                usc[3] = 0                
                                usc[9] = 0
                        elif 'xy' in siteSym or '+-0' in siteSym:
                            if np.allclose(dU[0,0,:],dUT[0,1,:]*sdet):
                                CSI[0][4:6] = [[12,0,0],[12,0,0]]
                                CSI[0][6:8] = [[11,0,0],[11,0,0]]
                                CSI[1][4:6] = [[1.,0.,0.],[sdet,0.,0.]]
                                CSI[1][6:8] = [[1.,0.,0.],[sdet,0.,0.]]
                                usc[4:6] = 0
                                usc[6:8] = 0
                            
                    if debug: print (SSMT2text(ssop).replace(' ',''),sdet,ssdet,epsinv,usc)
                USC &= usc
            if debug: print (USC)
            if not np.any(dtau%.5):
                n = -1
                for i,U in enumerate(USC):
                    if U:
                        n += 1
                        CSI[0][i][0] = n+1
                        CSI[1][i][0] = 1.0
    
        return list(CSI),dU,dUTP
    
    def DoMag():
        delt6 = np.eye(6)*0.001
        dM = posFourier(tau,nH,delt6[:3],delt6[3:]) #+np.array(Mxyz)[:,nxs,nxs]
        dMTP = []
        CSI = [np.zeros((6,3),dtype='i'),np.zeros((6,3))]
        if siteSym == '1':
            CSI = [[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]],6*[[1.,0.,0.],]
        elif siteSym in ['-1','mmm',]:
            CSI = 3*[[0,0,0],]+[[1,0,0],[2,0,0],[3,0,0]],3*[[0.,0.,0.],]+3*[[1.,0.,0.],]
        elif siteSym in ['4(z)','422(z)']:
            CSI[0][0][0] = CSI[0][4][1] = 1
            CSI[1][0][0] = 1.0
            CSI[1][4][1] = -1.0
        elif siteSym in ['-4m2(z)','422(z)',]:
            CSI[0][5][0] = 1
            CSI[1][5][0] = 1.0
        elif siteSym in ['-32(100)','-3',]:
            CSI[0][2][0] = 1
            CSI[1][2][0] = 1.0
        elif siteSym in ['3',]:
            CSI[0][0][0] = CSI[0][3][0] = CSI[0][4][0] = 1
            CSI[1][0][0] = -np.sqrt(3.0)
            CSI[1][3][0] = 2.0
            CSI[1][4][0] = 1.0
        elif siteSym in ['622','2(100)','32(100)',]:
            CSI[0][0][0] = CSI[0][1][0] = CSI[0][3][0] = 1
            CSI[1][0][0] = 1.0
            CSI[1][1][0] = 2.0
            CSI[1][3][0] = np.sqrt(3.0)
        else:
              #3x6x12 modulated moment array (M,Spos,tau)& force positive
            CSI = [np.zeros((6,3),dtype='i'),np.zeros((6,3))]
            MSC = np.ones(6,dtype='i')
            dMTP = []
            for i in SdIndx:
                sop = Sop[i]
                ssop = SSop[i]
                sdet,ssdet,dtau,dT,tauT = getTauT(tau,sop,ssop,XYZ)
                msc = np.ones(6,dtype='i')
                dMT = posFourier(np.sort(tauT),nH,delt6[:3],delt6[3:])   #+np.array(XYZ)[:,nxs,nxs]
                dMT = np.inner(sop[0],dMT.T)    # X modulations array(3x6x49) -> array(3x49x6)
                dMT = np.swapaxes(dMT,1,2)      # back to array(3x6x49)
                dMT[:,:3,:] *= (ssdet*sdet)            # modify the sin component
                dMTP.append(dMT)
                for i in range(3):
                    if not np.allclose(dM[i,i,:],sdet*dMT[i,i,:]):
                        msc[i] = 0
                    if not np.allclose(dM[i,i+3,:],sdet*dMT[i,i+3,:]):
                        msc[i+3] = 0
                if np.any(dtau%.5) and ('1/2' in SSGData['modSymb'] or '1' in SSGData['modSymb']):
                    msc[3:6] = 0
                    CSI = [[[1,0,0],[2,0,0],[3,0,0], [1,0,0],[2,0,0],[3,0,0]],
                        [[1.,0.,0.],[1.,0.,0.],[1.,0.,0.], [1.,0.,0.],[1.,0.,0.],[1.,0.,0.]]]
                    if dT:
                        if '(x)' in siteSym:
                            CSI[1][3:] = [1./dT,0.,0.],[-dT,0.,0.],[-dT,0.,0.]
                            if 'm' in siteSym and len(SdIndx) == 1:
                                CSI[1][3:] = [1./dT,0.,0.],[-dT,0.,0.],[-dT,0.,0.]
                        elif '(y)' in siteSym:
                            CSI[1][3:] = [-dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]
                            if 'm' in siteSym and len(SdIndx) == 1:
                                CSI[1][3:] = [-dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]
                        elif '(z)' in siteSym:
                            CSI[1][3:] = [-dT,0.,0.],[-dT,0.,0.],[1./dT,0.,0.]
                            if 'm' in siteSym and len(SdIndx) == 1:
                                CSI[1][3:] = [-dT,0.,0.],[-dT,0.,0.],[1./dT,0.,0.]
                    else:
                        CSI[1][3:] = [0.,0.,0.],[0.,0.,0.],[0.,0.,0.]
                if '4/mmm' in laue:
                    if siteSym in ['4/mmm(z)',]:
                        CSI = 3*[[0,0,0],]+[[0,0,0],[0,0,0],[1,0,0]],3*[[0.,0.,0.],]+3*[[1.,0.,0.],]
                    if np.any(dtau%.5) and '1/2' in SSGData['modSymb']:
                        if '(xy)' in siteSym:
                            CSI[0] = [[1,0,0],[1,0,0],[2,0,0], [1,0,0],[1,0,0],[2,0,0]]
                            if dT:
                                CSI[1][3:] = [[1./dT,0.,0.],[1./dT,0.,0.],[-dT,0.,0.]]
                            else:
                                CSI[1][3:] = [0.,0.,0.],[0.,0.,0.],[0.,0.,0.]
                    if '(xy)' in siteSym or '(+-0)' in siteSym:
                        mul = 1
                        if '(+-0)' in siteSym:
                            mul = -1
                        if np.allclose(dM[0,0,:],dMT[1,0,:]):
                            CSI[0][3:5] = [[11,0,0],[11,0,0]]
                            CSI[1][3:5] = [[1.,0,0],[mul,0,0]]
                            msc[3:5] = 0
                        if np.allclose(dM[0,3,:],dMT[0,4,:]):
                            CSI[0][:2] = [[12,0,0],[12,0,0]]
                            CSI[1][:2] = [[1.,0,0],[mul,0,0]]
                            msc[:2] = 0
                MSC &= msc
                if debug: print (SSMT2text(ssop).replace(' ',''),sdet,ssdet,epsinv,msc)
            n = -1
            if debug: print (MSC)
            for i,M in enumerate(MSC):
                if M:
                    n += 1
                    CSI[0][i][0] = n+1
                    CSI[1][i][0] = 1.0

        return list(CSI),dM,dMTP
        
    if debug: print ('super space group: '+SSGData['SSpGrp'])
    xyz = np.array(XYZ)%1.
    SGOps = copy.deepcopy(SGData['SGOps'])
    laue = SGData['SGLaue']
    siteSym = SytSym(XYZ,SGData)[0].strip()
    if debug: print ('siteSym: '+siteSym)
    SSGOps = copy.deepcopy(SSGData['SSGOps'])
    #expand ops to include inversions if any
    if SGData['SGInv'] and not SGData['SGFixed']:
        for op,sop in zip(SGData['SGOps'],SSGData['SSGOps']):
            SGOps.append([-op[0],-op[1]%1.])
            SSGOps.append([-sop[0],-sop[1]%1.])
    #build set of sym ops around special position        
    SSop = []
    Sop = []
    Sdtau = []
    for iop,Op in enumerate(SGOps):         
        nxyz = (np.inner(Op[0],xyz)+Op[1])%1.
        if np.allclose(xyz,nxyz,1.e-4) and iop and MT2text(Op).replace(' ','') != '-X,-Y,-Z':
            SSop.append(SSGOps[iop])
            Sop.append(SGOps[iop])
            ssopinv = nl.inv(SSGOps[iop][0])
            mst = ssopinv[3][:3]
            epsinv = ssopinv[3][3]
            Sdtau.append(np.sum(mst*(XYZ-SGOps[iop][1])-epsinv*SSGOps[iop][1][3]))
    SdIndx = np.argsort(np.array(Sdtau))     # just to do in sensible order
    if debug: print ('special pos super operators: ',[SSMT2text(ss).replace(' ','') for ss in SSop])
    #setup displacement arrays
    tau = np.linspace(-1,1,49,True)
    #make modulation arrays - one parameter at a time
    if Stype == 'Sfrac':
        CSI,dF,dFTP = DoFrac()
    elif Stype == 'Spos':
        CSI,dF,dFTP = DoXYZ()
        CSI[0] = orderParms(CSI[0])
    elif Stype == 'Sadp':
        CSI,dF,dFTP = DoUij()
        CSI[0] = orderParms(CSI[0]) 
    elif Stype == 'Smag':
        CSI,dF,dFTP = DoMag()

    if debug:
        return CSI,dF,dFTP
    else:
        return CSI,[],[]
    
def MustrainNames(SGData):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    if laue in ['m3','m3m']:
        return ['S400','S220']
    elif laue in ['6/m','6/mmm','3m1','3']:
        return ['S400','S004','S202']
    elif laue in ['31m',]:
        return ['S400','S004','S202','S301']
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
        
def HStrainVals(HSvals,SGData):
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    DIJ = np.zeros(6)
    if laue in ['m3','m3m']:
        DIJ[:3] = [HSvals[0],HSvals[0],HSvals[0]]
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        DIJ[:4] = [HSvals[0],HSvals[0],HSvals[1],HSvals[0]]
    elif laue in ['3R','3mR']:
        DIJ = [HSvals[0],HSvals[0],HSvals[0],HSvals[1],HSvals[1],HSvals[1]]
    elif laue in ['4/m','4/mmm']:
        DIJ[:3] = [HSvals[0],HSvals[0],HSvals[1]]
    elif laue in ['mmm']:
        DIJ[:3] = [HSvals[0],HSvals[1],HSvals[2]]
    elif laue in ['2/m']:
        DIJ[:3] = [HSvals[0],HSvals[1],HSvals[2]]
        if uniq == 'a':
            DIJ[5] = HSvals[3]
        elif uniq == 'b':
            DIJ[4] = HSvals[3]
        elif uniq == 'c':
            DIJ[3] = HSvals[3]
    else:
        DIJ = [HSvals[0],HSvals[1],HSvals[2],HSvals[3],HSvals[4],HSvals[5]]
    return DIJ

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
    elif laue in ['6/m','6/mmm','3m1','3']:
        Strm.append(h**4+k**4+2.0*k*h**3+2.0*h*k**3+3.0*(h*k)**2)
        Strm.append(l**4)
        Strm.append(3.0*((h*l)**2+(k*l)**2+h*k*l**2))
    elif laue in ['31m',]:
        Strm.append(h**4+k**4+2.0*k*h**3+2.0*h*k**3+3.0*(h*k)**2)
        Strm.append(l**4)
        Strm.append(3.0*((h*l)**2+(k*l)**2+h*k*l**2))
        Strm.append(4.0*l*h**3)
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

def MuShklMean(SGData,Amat,Shkl):
    
    def genMustrain(xyz,Shkl):
        uvw = np.inner(Amat.T,xyz)
        Strm = np.array(MustrainCoeff(uvw,SGData))
        Sum = np.sum(np.multiply(Shkl,Strm))
        Sum = np.where(Sum > 0.01,Sum,0.01)
        Sum = np.sqrt(Sum)
        return Sum*xyz
        
    PHI = np.linspace(0.,360.,30,True)
    PSI = np.linspace(0.,180.,30,True)
    X = np.outer(npcosd(PHI),npsind(PSI))
    Y = np.outer(npsind(PHI),npsind(PSI))
    Z = np.outer(np.ones(np.size(PHI)),npcosd(PSI))
    XYZ = np.dstack((X,Y,Z))
    XYZ = np.nan_to_num(np.apply_along_axis(genMustrain,2,XYZ,Shkl))
    return np.sqrt(np.sum(XYZ**2)/900.)
    
def Muiso2Shkl(muiso,SGData,cell):
    "this is to convert isotropic mustrain to generalized Shkls"
    import GSASIIlattice as G2lat
    A = G2lat.cell2AB(cell)[0]
    
    def minMus(Shkl,muiso,H,SGData,A):
        U = np.inner(A.T,H)
        S = np.array(MustrainCoeff(U,SGData))
        nS = S.shape[0]
        Sum = np.sqrt(np.sum(np.multiply(S,Shkl[:nS,nxs]),axis=0))
        rad = np.sqrt(np.sum((Sum[:,nxs]*H)**2,axis=1))
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
    elif laue in ['6/m','6/mmm']:
        S0 = [1000.,1000.,1000.]
    elif laue in ['31m','3','3m1']:
        S0 = [1000.,1000.,1000.]
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
        
def SytSym(XYZ,SGData):
    '''
    Generates the number of equivalent positions and a site symmetry code for a specified coordinate and space group

    :param XYZ: an array, tuple or list containing 3 elements: x, y & z
    :param SGData: from SpcGroup
    :Returns: a four element tuple:

     * The 1st element is a code for the site symmetry (see GetKNsym)
     * The 2nd element is the site multiplicity
     * Ndup number of overlapping operators
     * dupDir Dict - dictionary of overlapping operators

    '''
    if SGData['SpGrp'] == 'P 1':
        return '1',1,1,{}
    Mult = 1
    Isym = 0
    if SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
        Isym = 1073741824
    Jdup = 0
    Ndup = 0
    dupDir = {}
    inv = SGData['SGInv']+1
    icen = SGData['SGCen']
    Ncen = len(icen)
    if SGData['SGFixed']:       #already in list of operators
        inv = 1
    if SGData['SGGray'] and Ncen > 1: Ncen //= 2
    Xeqv = list(GenAtom(np.array(XYZ)%1.,SGData,True))
#    for xeqv in Xeqv:   print(xeqv)
    IRT = PackRot(SGData['SGOps'])
    L = -1
    for ic,cen in enumerate(icen[:Ncen]):
        for invers in range(int(inv)):
            for io,ops in enumerate(SGData['SGOps']):
                irtx = (1-2*invers)*IRT[io]
                L += 1
                if not Xeqv[L][1]:
                    Ndup = io
                    Jdup += 1
                    jx = GetOprPtrNumber(str(irtx))   #[KN table no,op name,KNsym ptr]
                    if jx < 39:
                        px = GetOprName(str(irtx))
                        if Xeqv[L][-1] < 0:
                            if '(' in px:
                                px = px.split('(')
                                px[0] += "'"
                                px = '('.join(px)
                            else:    
                                px += "'"
                        dupDir[px] = L
                        Isym += 2**(jx-1)
    if Isym == 1073741824: Isym = 0
    try:
        Mult = len(SGData['SGOps'])*Ncen*inv//Jdup
    except: # patch because Jdup is not getting incremented for most atoms!
        Mult = 0
        
    return GetKNsym(str(Isym)),Mult,Ndup,dupDir

def AtomDxSymFix(Dx,SytSym,CSIX):
    ''' Applies site symmetry restrictions to atom position shifts. 1st parameter value 
    of each kind encountered is assumed to be the independent one. Needed for ISODISTORT mode shifts.
    

    '''
    if SytSym == '1':
        return Dx    
    newDx = []
    for ix in [0,1,2]:
        cx = CSIX[2][ix]
        if cx > 0:
            newDx.append(CSIX[1][ix]*Dx[cx-1])
        else:
            newDx.append(0.0)
    return np.array(newDx)
    
    
   
def MagSytSym(SytSym,dupDir,SGData):
    '''
    site sym operations: 1,-1,2,3,-3,4,-4,6,-6,m need to be marked if spin inversion
    '''
    SGData['GenSym'],SGData['GenFlg'] = GetGenSym(SGData)[:2]
#    print('SGPtGrp',SGData['SGPtGrp'],'SytSym',SytSym,'MagSpGrp',SGData['MagSpGrp'])
#    print('dupDir',dupDir)
    SplitSytSym = SytSym.split('(')
    if SGData['SGGray']:
        return SytSym+"1'"
    if SytSym in ['1','-1','2','3','6','m']:       #axial position
        return SytSym
    if SplitSytSym[0] == SGData['SGPtGrp']:     #simple cases
        try:
            MagSytSym = SGData['MagSpGrp'].split()[1]
        except IndexError:
            MagSytSym = SGData['MagSpGrp'][1:].strip("1'")
        if len(SplitSytSym) > 1:
            MagSytSym += '('+SplitSytSym[1]
        return MagSytSym
    if len(dupDir) == 1:
        return list(dupDir.keys())[0]
    
    
    if '2/m' in SytSym:         #done I think; last 2wo might be not needed
        ops = {'(x)':['2(x)','m(x)'],'(y)':['2(y)','m(y)'],'(z)':['2(z)','m(z)'],
               '(100)':['2(100)','m(100)'],'(010)':['2(010)','m(010)'],'(001)':['2(001)','m(001)'],
               '(120)':['2(120)','m(120)'],'(210)':['2(210)','m(210)'],'(+-0)':['2(+-0)','m(+-0)'],
               '(110)':['2(110)','m(110)']}
    
    elif '4/mmm' in SytSym:
        ops = {'(x)':['4(x)','m(x)','m(y)','m(0+-)'],   #m(0+-) for cubic m3m?
               '(y)':['4(y)','m(y)','m(z)','m(+0-)'],   #m(+0-)
               '(z)':['4(z)','m(z)','m(x)','m(+-0)']}   #m(+-0)
    elif '4mm' in SytSym:
        ops = {'(x)':['4(x)','m(y)','m(yz)'],'(y)':['4(y)','m(z)','m(xz)'],'(z)':['4(z)','m(x)','m(110)']}
    elif '422' in SytSym:
        ops = {'(x)':['4(x)','2(y)','2(yz)'],'(y)':['4(y)','2(z)','2(xz)'],'(z)':['4(z)','2(x)','2(110)']}
    elif '-4m2' in SytSym:
        ops = {'(x)':['-4(x)','m(x)','2(yz)'],'(y)':['-4(y)','m(y)','2(xz)'],'(z)':['-4(z)','m(z)','2(110)']}
    elif '-42m' in SytSym:
        ops = {'(x)':['-4(x)','2(y)','m(yz)'],'(y)':['-4(y)','2(z)','m(xz)'],'(z)':['-4(z)','2(x)','m(110)']}
    elif '-4' in SytSym:
        ops = {'(x)':['-4(x)',],'(y)':['-4(y)',],'(z)':['-4(z)',],}
    elif '4' in SytSym:
        ops = {'(x)':['4(x)',],'(y)':['4(y)',],'(z)':['4(z)',],}

    elif '222' in SytSym:
        ops = {'':['2(x)','2(y)','2(z)'],
                   '(x)':['2(y)','2(z)','2(x)'],'(y)':['2(x)','2(z)','2(y)'],'(z)':['2(x)','2(y)','2(z)'],
                   '(100)':['2(z)','2(100)','2(120)',],'(010)':['2(z)','2(010)','2(210)',],
                   '(110)':['2(z)','2(110)','2(+-0)',],}
    elif 'mm2' in SytSym:
        ops = {'(x)':['m(y)','m(z)','2(x)'],'(y)':['m(x)','m(z)','2(y)'],'(z)':['m(x)','m(y)','2(z)'],
               '(xy)':['m(+-0)','m(z)','2(110)'],'(yz)':['m(0+-)','m(xz)','2(yz)'],     #not 2(xy)!
               '(xz)':['m(+0-)','m(y)','2(xz)'],'(z100)':['m(100)','m(120)','2(z)'],
               '(z010)':['m(010)','m(210)','2(z)'],'(z110)':['m(110)','m(+-0)','2(z)'],
               '(+-0)':[ 'm(110)','m(z)','2(+-0)'],'(d100)':['m(yz)','m(0+-)','2(xz)'],
               '(d010)':['m(xz)','m(+0-)','2(y)'],'(d001)':['m(110)','m(+-0)','2(z)'],
               '(210)':['m(z)','m(010)','2(210)'],'(120)':['m(z)','m(100)','2(120)'],
               '(100)':['m(z)','m(120)','2(100)',],'(010)':['m(z)','m(210)','2(010)',],
               '(110)':['m(z)','m(+-0)','2(110)',],}
    elif 'mmm' in SytSym:
        ops = {'':['m(x)','m(y)','m(z)'],
                   '(100)':['m(z)','m(100)','m(120)',],'(010)':['m(z)','m(010)','m(210)',],
                   '(110)':['m(z)','m(110)','m(+-0)',],
                   '(x)':['m(x)','m(y)','m(z)'],'(y)':['m(x)','m(y)','m(z)'],'(z)':['m(x)','m(y)','m(z)'],}
        
    elif '32' in SytSym:
        ops = {'(120)':['3','2(120)',],'(100)':['3','2(100)'],'(111)':['3(111)','2(x)']}
    elif '23' in SytSym:
        ops = {'':['2(x)','3(111)']}
    elif 'm3' in SytSym:
        ops = {'':['m(x)','3(111)',]}
    elif '3m' in SytSym:
        ops = {'(111)':['3(111)','m(+-0)',],'(+--)':['3(+--)','m(0+-)',],
               '(-+-)':['3(-+-)','m(+0-)',],'(--+)':['3(--+)','m(+-0)',],
               '(100)':['3','m(100)'],'(120)':['3','m(210)',]}
    
    if SytSym.split('(')[0] in ['6/m','6mm','-6m2','622','-6','-3','-3m','-43m',]:     #not simple cases
        MagSytSym = SytSym
        if "-1'" in dupDir:
            if '-6' in SytSym:
                MagSytSym = MagSytSym.replace('-6',"-6'")
            elif '-3m' in SytSym:
                MagSytSym = MagSytSym.replace('-3m',"-3'm'")
            elif '-3' in SytSym:
                MagSytSym = MagSytSym.replace('-3',"-3'")
        elif '-6m2' in SytSym:
            if "m'(110)" in dupDir:
                MagSytSym = "-6m'2'("+SytSym.split('(')[1]
        elif '6/m' in SytSym:
            if "m'(z)" in dupDir:
                MagSytSym = "6'/m'"
        elif '6mm' in SytSym:
            if "m'(110)" in dupDir:
                MagSytSym = "6'm'm"
        elif '-43m' in SytSym:
            if "m'(110)" in dupDir:
                MagSytSym = "-43m'"
        return MagSytSym
    try:
        axis = '('+SytSym.split('(')[1]
    except IndexError:
        axis = ''
    MagSytSym = ''
    for m in ops[axis]:
        if m in dupDir:
            MagSytSym += m.split('(')[0]
        else:
            MagSytSym += m.split('(')[0]+"'"
        if '2/m' in SytSym and '2' in m:
            MagSytSym += '/'
        if '-3/m' in SytSym:
            MagSytSym = '-'+MagSytSym
        
    MagSytSym += axis
# some exceptions & special rules          
    if MagSytSym == "4'/m'm'm'": MagSytSym = "4/m'm'm'"
    return MagSytSym
    
#    if len(GenSym) == 3:
#        if SGSpin[1] < 0:
#            if 'mm2' in SytSym:
#                MagSytSym = "m'm'2"+'('+SplitSytSym[1]
#            else:   #bad rule for I41/a
#                MagSytSym = SplitSytSym[0]+"'"
#                if len(SplitSytSym) > 1:
#                    MagSytSym += '('+SplitSytSym[1]
#        else:
#            MagSytSym = SytSym
#        if len(SplitSytSym) >1:
#            if "-4'"+'('+SplitSytSym[1] in dupDir:
#                MagSytSym = MagSytSym.replace('-4',"-4'")
#            if "-6'"+'('+SplitSytSym[1] in dupDir:
#                MagSytSym = MagSytSym.replace('-6',"-6'")
#        return MagSytSym
#            
    return SytSym

def UpdateSytSym(Phase):
    ''' Update site symmetry/site multiplicity after space group/BNS lattice change
    '''
    generalData = Phase['General']
    SGData = generalData['SGData']
    Atoms = Phase['Atoms']
    cx,ct,cs,cia = generalData['AtomPtrs']
    for atom in Atoms:
        XYZ = atom[cx:cx+3]
        sytsym,Mult = SytSym(XYZ,SGData)[:2]
        sytSym,Mul,Nop,dupDir = SytSym(XYZ,SGData)
        atom[cs] = sytsym
        if generalData['Type'] == 'magnetic':
            magSytSym = MagSytSym(sytSym,dupDir,SGData)
            atom[cs] = magSytSym
        atom[cs+1] = Mult
    return
    
def ElemPosition(SGData):
    ''' Under development. 
    Object here is to return a list of symmetry element types and locations suitable
    for say drawing them.
    So far I have the element type... getting all possible locations without lookup may be impossible!
    '''
    Inv = SGData['SGInv']
    eleSym = {-3:['','-1'],-2:['',-6],-1:['2','-4'],0:['3','-3'],1:['4','m'],2:['6',''],3:['1','']}
    # get operators & expand if centrosymmetric
    SymElements = []
    Ops = SGData['SGOps']
    opM = np.array([op[0].T for op in Ops])
    opT = np.array([op[1] for op in Ops])
    if Inv:
        opM = np.concatenate((opM,-opM))
        opT = np.concatenate((opT,-opT))
    opMT = list(zip(opM,opT))
    for M,T in opMT[1:]:        #skip I
        Dt = int(nl.det(M))
        Tr = int(np.trace(M))
        Dt = -(Dt-1)//2
        Es = eleSym[Tr][Dt]
        if Dt:              #rotation-inversion
            I = np.eye(3)
            if Tr == 1:     #mirrors/glides
                if np.any(T):       #glide
                    M2 = np.inner(M,M)
                    MT = np.inner(M,T)+T
                    print ('glide',Es,MT)
                    print (M2)
                else:               #mirror
                    print ('mirror',Es,T)
                    print (I-M)
                X = [-1,-1,-1]
            elif Tr == -3:  # pure inversion
                X = np.inner(nl.inv(I-M),T)
                print ('inversion',Es,X)
            else:           #other rotation-inversion
                M2 = np.inner(M,M)
                MT = np.inner(M,T)+T
                print ('rot-inv',Es,MT)
                print (M2)
                X = [-1,-1,-1]
        else:               #rotations
            print ('rotation',Es)
            X = [-1,-1,-1]
        SymElements.append([Es,X])
        
    return SymElements
    
def ApplyStringOps(A,SGData,X,Uij=[]):
    'Needs a doc string'
    SGOps = SGData['SGOps']
    SGCen = SGData['SGCen']
    Ax = A.split('+')
    Ax[0] = int(Ax[0])
    iC = 1
    if Ax[0] < 0:
        iC = -1
    Ax[0] = abs(Ax[0])
    nA = Ax[0]%100-1
    cA = Ax[0]//100
    Cen = SGCen[cA]
    M,T = SGOps[nA]
    if len(Ax)>1:
        cellA = Ax[1].split(',')
        cellA = np.array([int(a) for a in cellA])
    else:
        cellA = np.zeros(3)
    newX = Cen+iC*(np.inner(M,X).T+T)+cellA
    if len(Uij):
        U = Uij2U(Uij)
        U = np.inner(M,np.inner(U,M).T)
        newUij = U2Uij(U)
        return [newX,newUij]
    else:
        return newX
        
def ApplyStringOpsMom(A,SGData,SSGData,Mom):
    '''Applies string operations to modulated magnetic moment components used in drawing
    Drawing matches Bilbao MVISUALIZE
    '''
    SGOps = SGData['SGOps']
    Ax = A.split('+')
    Ax[0] = int(Ax[0])
    iAx = abs(Ax[0])
    nA = iAx%100-1
    if SGData['SGInv'] and not SGData['SGFixed']:
        nC = 2*len(SGOps)*(iAx//100)
    else:
        nC = len(SGOps)*(iAx//100)
    NA = nA
    if Ax[0] < 0:
        NA += len(SGOps)
    M,T = SGOps[nA]
    newMom = np.inner(Mom,M).T*SGData['SpnFlp'][NA+nC]*nl.det(M)
    if SSGData is not None:
        if SSGData['SSGCen'][iAx//100][3]:     #flip spin for BNS centered atoms
            newMom *= -1.
    return newMom
        
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
    cA = Ax[0]//100;  cB = Bx[0]//100
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
    return [U[0][0],U[1][1],U[2][2],U[0][1],U[0][2],U[1][2]]
    
def Uij2U(Uij):
    #returns the thermal motion tensor U from Uij as numpy array
    return np.array([[Uij[0],Uij[3],Uij[4]],[Uij[3],Uij[1],Uij[5]],[Uij[4],Uij[5],Uij[2]]])

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
    gray = ''
    if "1'" in rspc:
        gray = " 1'"
        rspc = rspc.replace("1'",'')
    rspc = rspc.replace("'",'')
    if rspc[-1:] == 'H': # hexagonal is assumed and thus can be ignored
        rspc = rspc[:-1]
    if rspc[1:3] in ['M3','N3','A3','D3']:      #fix cubic old style
        rspc.replace('3','-3')
    bns = -1
    try:
        bns = rspc.index('_')
        rspc = rspc.replace(rspc[bns:bns+2],'')
    except ValueError:
        pass
    # look for a match in the spacegroup lists
    for i in spglist.values():
        for spc in i:
            if rspc == spc.replace(' ','').upper():
                return spc+gray+rhomb
    # how about the post-2002 orthorhombic names?
    if rspc in sgequiv_2002_orthorhombic:
        return sgequiv_2002_orthorhombic[rspc]+gray
    else:
    # not found
        return ''
    
def fullHM2shortHM(SpcGp):
    ''' Accepts a full H-M space group symbol and returns a short H-M symbol that the space group 
    interpreter can translate
    '''
    fields = SpcGp.split()
    if len(fields) == 4 and fields[1] == '1' and fields[3] == '1': #b-unique monoclinics
        return fields[0]+' '+fields[2]
    if '/' not in SpcGp or len(fields) < 3:
        return SpcGp
    if '-3' in fields[1]:   #trigonals
        fields[2] = fields.split('/')[1]
        return ' '.join(fields)
    if '-3' in fields[2]:       #cubics
        fields[1] = fields.split('/')[1]
        if len(fields) == 4:
            fields[3].split('/')[1]
        return ' '.join(fields)
    if '6' in fields[1] and '/' in fields[2]:    #hexagonals
        fields[2] = fields[2].split('/')[1]
        fields[3] = fields[3].split('/')[1]
        return ' '.join(fields)
    if '4' in fields[1] and '/' in fields[2] and len(fields) == 4:   #tetragonals
        fields[2] = fields[2].split('/')[1]
        fields[3] = fields[3].split('/')[1]
        return ' '.join(fields)
    if len(fields) == 4 and '/' in fields[1] and '/' in fields[2] and '/' in fields[3]:    #orthorhombics
        if fields[2] == '1':     #skip a-unique monoclinics
            return SpcGp
        fields[1] = fields[1].split('/')[1]
        fields[2] = fields[2].split('/')[1]
        fields[3] = fields[3].split('/')[1]
        return ' '.join(fields)
    return SpcGp
    
def offsetNorm(x):
    '''Translate a coordinate (or vector of symmetry offsets) into 
    the range -1/6 < x <= 5/6, matching GSAS-II use.
    '''
    negoff = 1/6  # negative offset
    return np.round(((x + negoff - 1e-9) % 1) - negoff + 1e-9, 8)
        
def ParseXYZ(sym):
    '''Parse a set of space group operations, such as '-x,-y,-z' 
    or 'x-y,z,-x' etc. into an algebraic form returning a 
    3x3 matrix, M, and a length 3 vector, O, where the symmetry 
    operation can be applied as 

       M . (x,y,z) + O 

    where M is a 3x3 matrix and O is a displacement vector. 
    This is the same form generated by :func:`SpcGroup`. 

    Note that this code will process offsets either before the coordinate 
    multipliers or after, so that either of these

       1/2-y,1/2+x,z or 
       -y+1/2,x-0.5,z

    will both be parsed properly. Returns None if a symbol cannot be parsed.
    '''
    xyz = sym.split(',')
    if len(xyz) != 3:
        print(f'error with: {sym}')
        return None
    multMat = []
    offsetVec = []
    for o in xyz:
        o = o.lower()
        op = [0,0,0]
        offset = copy.copy(o)
        for i,v in enumerate(('x','y','z')):
            if '-'+v in o:
                op[i] = -1
                offset = offset.replace('-'+v,'')
            elif v in o:
                op[i] = 1
                offset = offset.replace('+'+v,'')
                offset = offset.replace(v,'')
        if not offset.strip():
            offset = 0.0
        else:
            offset = eval(offset)
        offsetVec.append(offset)
#        offsetVec.append(offsetNorm(offset))
        multMat.append(op)
    return np.array(multMat),np.array(offsetVec)

def CompareSym(symList,sgName=None,SGData=None):
    '''Compare symmetry generated by GSAS-II from a space group name
    with a list of operators from some other source.

    :param list symList: a list of symmetry operations (such as 'x,y,z'
      or '-Y,X,Z'). The code is fairly forgiving in that 1/2-X or X-1/2 or 
      even -X+0.5 can all be accepted. This is typically read from a CIF. 
    :param str sgName: a space group name. Need not follow GSAS-II 
      convention for spaces, etc (see :func:`GSASIIspc.StandardizeSpcName`)
    :param SGData: a GSAS-II symmetry objectspace group name. Need not follow GSAS-II 
      convention for spaces, etc (see :func:`GSASIIspc.StandardizeSpcName`)
    '''
    if sgName:
        sgName = StandardizeSpcName(sgName) # Deal with non-standard spaces in name
        # use GSAS-II sym code to get a list of symmetry matrices from name
        matList = AllOps(SpcGroup(sgName)[1])[2]
    elif SGData:
        matList = AllOps(SGData)[2]
    else:
        raise Exception('CompareSym called without symmetry input')

    for sym in symList:
        syMat = ParseXYZ(sym)

        for i,mat in enumerate(matList):
            if np.allclose(mat[0],syMat[0]) and np.allclose(offsetNorm(mat[1]),offsetNorm(syMat[1])):
                break
        else:
            print(f'Symmetry element from input, {sym}, not matched in GSAS-II setting')
            return False
        del matList[i] # matched to this, remove it
    return True

def SpaceGroupNumber(spcgroup):
    '''Determine the space group number for a space group from H-M name.
    Will work for non-standard groups insofar as we can normalize -- far 
    from perfect.
    '''
    SGNo = -1
    SpcGp = StandardizeSpcName(spcgroup)
    if not SpcGp:
        return SGNo
    try:
        SGNo = spgbyNum.index(SpcGp)
    except ValueError:
        pass
    return SGNo

def GetHallSpaceGroup(SGData):
    '''Determine the Hall space group symbol for a GSAS-II space group object, 
    if it exists (it will not if non-standard centering is used, perhaps also 
    for other cases). Will return None if not found.
    '''
    try:
        key = SGData['SpGrp'].lower().replace(' ','')
        return spgHall[key][0]
    except:
        pass
    # H is also allowed as suffix for Hexagonal & trigonal SGs (remove it)
    if key[0] == 'p' and key[-1] == 'h':
        try:
            return spgHall[key[:-1]][0]
        except:
            pass
    if key[1:3] in ['m3','n3','a3','d3']: # fix old-style cubic names
        keyB3 = key.replace('3','-3')
        if keyB3 in spgHall: return spgHall[keyB3][0]
    # give up
    return None

spgbyNum = []
'''Space groups indexed by number'''
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    spgbyNum = [None,
        'P 1','P -1',                                                   #1-2
        'P 2','P 21','C 2','P m','P c','C m','C c','P 2/m','P 21/m',
        'C 2/m','P 2/c','P 21/c','C 2/c',                               #3-15
        'P 2 2 2','P 2 2 21','P 21 21 2','P 21 21 21',
        'C 2 2 21','C 2 2 2','F 2 2 2','I 2 2 2','I 21 21 21',
        'P m m 2','P m c 21','P c c 2','P m a 2','P c a 21',
        'P n c 2','P m n 21','P b a 2','P n a 21','P n n 2',
        'C m m 2','C m c 21','C c c 2',
        'A m m 2','A b m 2','A m a 2','A b a 2',
        'F m m 2','F d d 2','I m m 2','I b a 2','I m a 2',
        'P m m m','P n n n','P c c m','P b a n',
        'P m m a','P n n a','P m n a','P c c a','P b a m','P c c n',
        'P b c m','P n n m','P m m n','P b c n','P b c a','P n m a',
        'C m c m','C m c a','C m m m','C c c m','C m m a','C c c a',
        'F m m m', 'F d d d',
        'I m m m','I b a m','I b c a','I m m a',                        #16-74
        'P 4','P 41','P 42','P 43',
        'I 4','I 41',
        'P -4','I -4','P 4/m','P 42/m','P 4/n','P 42/n',
        'I 4/m','I 41/a',
        'P 4 2 2','P 4 21 2','P 41 2 2','P 41 21 2','P 42 2 2',
        'P 42 21 2','P 43 2 2','P 43 21 2',
        'I 4 2 2','I 41 2 2',
        'P 4 m m','P 4 b m','P 42 c m','P 42 n m','P 4 c c','P 4 n c',
        'P 42 m c','P 42 b c',
        'I 4 m m','I 4 c m','I 41 m d','I 41 c d',
        'P -4 2 m','P -4 2 c','P -4 21 m','P -4 21 c','P -4 m 2',
        'P -4 c 2','P -4 b 2','P -4 n 2',
        'I -4 m 2','I -4 c 2','I -4 2 m','I -4 2 d',
        'P 4/m m m','P 4/m c c','P 4/n b m','P 4/n n c','P 4/m b m',
        'P 4/m n c','P 4/n m m','P 4/n c c','P 42/m m c','P 42/m c m',
        'P 42/n b c','P 42/n n m','P 42/m b c','P 42/m n m','P 42/n m c',
        'P 42/n c m',
        'I 4/m m m','I 4/m c m','I 41/a m d','I 41/a c d',
        'P 3','P 31','P 32','R 3','P -3','R -3',
        'P 3 1 2','P 3 2 1','P 31 1 2','P 31 2 1','P 32 1 2','P 32 2 1',
        'R 3 2',
        'P 3 m 1','P 3 1 m','P 3 c 1','P 3 1 c',
        'R 3 m','R 3 c',
        'P -3 1 m','P -3 1 c','P -3 m 1','P -3 c 1',
        'R -3 m','R -3 c',                                               #75-167
        'P 6','P 61',
        'P 65','P 62','P 64','P 63','P -6','P 6/m','P 63/m','P 6 2 2',
        'P 61 2 2','P 65 2 2','P 62 2 2','P 64 2 2','P 63 2 2','P 6 m m',
        'P 6 c c','P 63 c m','P 63 m c','P -6 m 2','P -6 c 2','P -6 2 m',
        'P -6 2 c','P 6/m m m','P 6/m c c','P 63/m c m','P 63/m m c',   #168-194
        'P 2 3','F 2 3','I 2 3','P 21 3','I 21 3','P m 3','P n 3',
        'F m -3','F d -3','I m -3',
        'P a -3','I a -3','P 4 3 2','P 42 3 2','F 4 3 2','F 41 3 2',
        'I 4 3 2','P 43 3 2','P 41 3 2','I 41 3 2','P -4 3 m',
        'F -4 3 m','I -4 3 m','P -4 3 n','F -4 3 c','I -4 3 d',
        'P m -3 m','P n -3 n','P m -3 n','P n -3 m',
        'F m -3 m','F m -3 c','F d -3 m','F d -3 c',
        'I m -3 m','I a -3 d',]                                       #195-230
altSettingOrtho = {}
''' A dictionary of alternate settings for orthorhombic unit cells
'''
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    altSettingOrtho = {
        'P 2 2 2' :{'abc':'P 2 2 2','cab':'P 2 2 2','bca':'P 2 2 2','acb':'P 2 2 2','bac':'P 2 2 2','cba':'P 2 2 2'},
        'P 2 2 21' :{'abc':'P 2 2 21','cab':'P 21 2 2','bca':'P 2 21 2','acb':'P 2 21 2','bac':'P 2 2 21','cba':'P 21 2 2'},
        'P 21 21 2':{'abc':'P 21 21 2','cab':'P 2 21 21','bca':'P 21 2 21','acb':'P 21 2 21','bac':'P 21 21 2','cba':'P 2 21 21'},
        'P 21 21 21':{'abc':'P 21 21 21','cab':'P 21 21 21','bca':'P 21 21 21','acb':'P 21 21 21','bac':'P 21 21 21','cba':'P 21 21 21'},
        'C 2 2 21':{'abc':'C 2 2 21','cab':'A 21 2 2','bca':'B 2 21 2','acb':'B 2 21 2','bac':'C 2 2 21','cba':'A 21 2 2'},
        'C 2 2 2':{'abc':'C 2 2 2','cab':'A 2 2 2','bca':'B 2 2 2','acb':'B 2 2 2','bac':'C 2 2 2','cba':'A 2 2 2'},
        'F 2 2 2':{'abc':'F 2 2 2','cab':'F 2 2 2','bca':'F 2 2 2','acb':'F 2 2 2','bac':'F 2 2 2','cba':'F 2 2 2'},
        'I 2 2 2':{'abc':'I 2 2 2','cab':'I 2 2 2','bca':'I 2 2 2','acb':'I 2 2 2','bac':'I 2 2 2','cba':'I 2 2 2'},
        'I 21 21 21':{'abc':'I 21 21 21','cab':'I 21 21 21','bca':'I 21 21 21','acb':'I 21 21 21','bac':'I 21 21 21','cba':'I 21 21 21'},
        'P m m 2':{'abc':'P m m 2','cab':'P 2 m m','bca':'P m 2 m','acb':'P m 2 m','bac':'P m m 2','cba':'P 2 m m'},
        'P m c 21':{'abc':'P m c 21','cab':'P 21 m a','bca':'P b 21 m','acb':'P m 21 b','bac':'P c m 21','cba':'P 21 a m'},
        'P c c 2':{'abc':'P c c 2','cab':'P 2 a a','bca':'P b 2 b','acb':'P b 2 b','bac':'P c c 2','cba':'P 2 a a'},
        'P m a 2':{'abc':'P m a 2','cab':'P 2 m b','bca':'P c 2 m','acb':'P m 2 a','bac':'P b m 2','cba':'P 2 c m'},
        'P c a 21':{'abc':'P c a 21','cab':'P 21 a b','bca':'P c 21 b','acb':'P b 21 a','bac':'P b c 21','cba':'P 21 c a'},
        'P n c 2':{'abc':'P n c 2','cab':'P 2 n a','bca':'P b 2 n','acb':'P n 2 b','bac':'P c n 2','cba':'P 2 a n'},
        'P m n 21':{'abc':'P m n 21','cab':'P 21 m n','bca':'P n 21 m','acb':'P m 21 n','bac':'P n m 21','cba':'P 21 n m'},
        'P b a 2':{'abc':'P b a 2','cab':'P 2 c b','bca':'P c 2 a','acb':'P c 2 a','bac':'P b a 2','cba':'P 2 c b'},
        'P n a 21':{'abc':'P n a 21','cab':'P 21 n b','bca':'P c 21 n','acb':'P n 21 a','bac':'P b n 21','cba':'P 21 c n'},
        'P n n 2':{'abc':'P n n 2','cab':'P 2 n n','bca':'P n 2 n','acb':'P n 2 n','bac':'P n n 2','cba':'P 2 n n'},
        'C m m 2':{'abc':'C m m 2','cab':'A 2 m m','bca':'B m 2 m','acb':'B m 2 m','bac':'C m m 2','cba':'A 2 m m'},
        'C m c 21':{'abc':'C m c 21','cab':'A 21 m a','bca':'B b 21 m','acb':'B m 21 b','bac':'C c m 21','cba':'A 21 a m'},
        'C c c 2':{'abc':'C c c 2','cab':'A 2 a a','bca':'B b 2 b','acb':'B b 2 b','bac':'C c c 2','cba':'A 2 a a'},
        'A m m 2':{'abc':'A m m 2','cab':'B 2 m m','bca':'C m 2 m','acb':'A m 2 m','bac':'B m m 2','cba':'C 2 m m'},
        'A b m 2':{'abc':'A b m 2','cab':'B 2 c m','bca':'C m 2 a','acb':'A c 2 m','bac':'B m a 2','cba':'C 2 m b'},
        'A m a 2':{'abc':'A m a 2','cab':'B 2 m b','bca':'C c 2 m','acb':'A m 2 a','bac':'B b m 2','cba':'C 2 c m'},
        'A b a 2':{'abc':'A b a 2','cab':'B 2 c b','bca':'C c 2 a','acb':'A c 2 a','bac':'B b a 2','cba':'C 2 c b'},
        'F m m 2':{'abc':'F m m 2','cab':'F 2 m m','bca':'F m 2 m','acb':'F m 2 m','bac':'F m m 2','cba':'F 2 m m'},
        'F d d 2':{'abc':'F d d 2','cab':'F 2 d d','bca':'F d 2 d','acb':'F d 2 d','bac':'F d d 2','cba':'F 2 d d'},
        'I m m 2':{'abc':'I m m 2','cab':'I 2 m m','bca':'I m 2 m','acb':'I m 2 m','bac':'I m m 2','cba':'I 2 m m'},
        'I b a 2':{'abc':'I b a 2','cab':'I 2 c b','bca':'I c 2 a','acb':'I c 2 a','bac':'I b a 2','cba':'I 2 c b'},
        'I m a 2':{'abc':'I m a 2','cab':'I 2 m b','bca':'I c 2 m','acb':'I m 2 a','bac':'I b m 2','cba':'I 2 c m'},
        'P m m m':{'abc':'P m m m','cab':'P m m m','bca':'P m m m','acb':'P m m m','bac':'P m m m','cba':'P m m m'},
        'P n n n':{'abc':'P n n n','cab':'P n n n','bca':'P n n n','acb':'P n n n','bac':'P n n n','cba':'P n n n'},
        'P c c m':{'abc':'P c c m','cab':'P m a a','bca':'P b m b','acb':'P b m b','bac':'P c c m','cba':'P m a a'},
        'P b a n':{'abc':'P b a n','cab':'P n c b','bca':'P c n a','acb':'P c n a','bac':'P b a n','cba':'P n c b'},
        'P m m a':{'abc':'P m m a','cab':'P b m m','bca':'P m c m','acb':'P m a m','bac':'P m m b','cba':'P c m m'},
        'P n n a':{'abc':'P n n a','cab':'P b n n','bca':'P n c n','acb':'P n a n','bac':'P n n b','cba':'P c n n'},
        'P m n a':{'abc':'P m n a','cab':'P b m n','bca':'P n c m','acb':'P m a n','bac':'P n m b','cba':'P c n m'},
        'P c c a':{'abc':'P c c a','cab':'P b a a','bca':'P b c b','acb':'P b a b','bac':'P c c b','cba':'P c a a'},
        'P b a m':{'abc':'P b a m','cab':'P m c b','bca':'P c m a','acb':'P c m a','bac':'P b a m','cba':'P m c b'},
        'P c c n':{'abc':'P c c n','cab':'P n a a','bca':'P b n b','acb':'P b n b','bac':'P c c n','cba':'P n a a'},
        'P b c m':{'abc':'P b c m','cab':'P m c a','bca':'P b m a','acb':'P c m b','bac':'P c a m','cba':'P m a b'},
        'P n n m':{'abc':'P n n m','cab':'P m n n','bca':'P n m n','acb':'P n m n','bac':'P n n m','cba':'P m n n'},
        'P m m n':{'abc':'P m m n','cab':'P n m m','bca':'P m n m','acb':'P m n m','bac':'P m m n','cba':'P n m m'},
        'P b c n':{'abc':'P b c n','cab':'P n c a','bca':'P b n a','acb':'P c n b','bac':'P c a n','cba':'P n a b'},
        'P b c a':{'abc':'P b c a','cab':'P b c a','bca':'P b c a','acb':'P c a b','bac':'P c a b','cba':'P c a b'},
        'P n m a':{'abc':'P n m a','cab':'P b n m','bca':'P m c n','acb':'P n a m','bac':'P m n b','cba':'P c m n'},
        'C m c m':{'abc':'C m c m','cab':'A m m a','bca':'B b m m','acb':'B m m b','bac':'C c m m','cba':'A m a m'},
        'C m c a':{'abc':'C m c a','cab':'A b m a','bca':'B b c m','acb':'B m a b','bac':'C c m b','cba':'A c a m'},
        'C m m m':{'abc':'C m m m','cab':'A m m m','bca':'B m m m','acb':'B m m m','bac':'C m m m','cba':'A m m m'},
        'C c c m':{'abc':'C c c m','cab':'A m a a','bca':'B b m b','acb':'B b m b','bac':'C c c m','cba':'A m a a'},
        'C m m a':{'abc':'C m m a','cab':'A b m m','bca':'B m c m','acb':'B m a m','bac':'C m m b','cba':'A c m m'},
        'C c c a':{'abc':'C c a a','cab':'A b a a','bca':'B b c b','acb':'B b a b','bac':'C c c b','cba':'A c a a'},
        'F m m m':{'abc':'F m m m','cab':'F m m m','bca':'F m m m','acb':'F m m m','bac':'F m m m','cba':'F m m m'},
        'F d d d':{'abc':'F d d d','cab':'F d d d','bca':'F d d d','acb':'F d d d','bac':'F d d d','cba':'F d d d'},
        'I m m m':{'abc':'I m m m','cab':'I m m m','bca':'I m m m','acb':'I m m m','bac':'I m m m','cba':'I m m m'},
        'I b a m':{'abc':'I b a m','cab':'I m c b','bca':'I c m a','acb':'I c m a','bac':'I b a m','cba':'I m c b'},
        'I b c a':{'abc':'I b c a','cab':'I b c a','bca':'I b c a','acb':'I c a b','bac':'I c a b','cba':'I c a b'},
        'I m m a':{'abc':'I m m a','cab':'I b m m','bca':'I m c m','acb':'I m a m','bac':'I m m  b','cba':'I c m m'},
        }
spg2origins = {}
''' A dictionary of all spacegroups that have 2nd settings; the value is the 
1st --> 2nd setting transformation vector as X(2nd) = X(1st)-V, nonstandard ones are included.
'''
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    spg2origins = {
        'P n n n':[-.25,-.25,-.25],
        'P b a n':[-.25,-.25,0],'P n c b':[0,-.25,-.25],'P c n a':[-.25,0,-.25],
        'P m m n':[-.25,-.25,0],'P n m m':[0,-.25,-.25],'P m n m':[-.25,0,-.25],
        'C c c a':[0,-.25,-.25],'C c c b':[-.25,0,-.25],'A b a a':[-.25,0,-.25],
        'A c a a':[-.25,-.25,0],'B b c b':[-.25,-.25,0],'B b a b':[0,-.25,-.25],
        'F d d d':[-.125,-.125,-.125],
        'P 4/n':[-.25,-.25,0],'P 42/n':[-.25,-.25,-.25],'I 41/a':[0,-.25,-.125],
        'P 4/n b m':[-.25,-.25,0],'P 4/n n c':[-.25,-.25,-.25],'P 4/n m m':[-.25,-.25,0],'P 4/n c c':[-.25,-.25,0],
        'P 42/n b c':[-.25,-.25,-.25],'P 42/n n m':[-.25,.25,-.25],'P 42/n m c':[-.25,.25,-.25],'P 42/n c m':[-.25,.25,-.25],
        'I 41/a m d':[0,.25,-.125],'I 41/a c d':[0,.25,-.125],
        'p n -3':[-.25,-.25,-.25],'F d -3':[-.125,-.125,-.125],'P n -3 n':[-.25,-.25,-.25],
        'P n -3 m':[-.25,-.25,-.25],'F d -3 m':[-.125,-.125,-.125],'F d -3 c':[-.375,-.375,-.375],
        'p n 3':[-.25,-.25,-.25],'F d 3':[-.125,-.125,-.125],'P n 3 n':[-.25,-.25,-.25],
        'P n 3 m':[-.25,-.25,-.25],'F d 3 m':[-.125,-.125,-.125],'F d - c':[-.375,-.375,-.375]}
spglist = {}
'''A dictionary of space groups as ordered and named in the pre-2002 International 
Tables Volume A, except that spaces are used following the GSAS convention to 
separate the different crystallographic directions.
Note that the symmetry codes here will recognize many non-standard space group 
symbols with different settings. They are ordered by Laue group
'''
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    spglist = {
    'P1' : ('P 1','P -1',), # 1-2
    'C1' : ('C 1','C -1',),
    'P2/m': ('P 2','P 21','P m','P a','P c','P n',
        'P 2/m','P 21/m','P 2/c','P 2/a','P 2/n','P 21/c','P 21/a','P 21/n',), #3-15
    'C2/m':('C 2','C m','C c','C n',
        'C 2/m','C 2/c','C 2/n',),
    'A2/m':('A 2','A m','A a','A n',
        'A 2/m','A 2/a','A 2/n',),
    'I2/m':('I 2','I m','I a','I n','I c',
        'I 2/m','I 2/a','I 2/c','I 2/n',),
   'Pmmm':('P 2 2 2',
        'P 2 2 21','P 21 2 2','P 2 21 2',
        'P 21 21 2','P 2 21 21','P 21 2 21',
        'P 21 21 21',
        'P m m 2','P 2 m m','P m 2 m',
        'P m c 21','P 21 m a','P b 21 m','P m 21 b','P c m 21','P 21 a m',
        'P c c 2','P 2 a a','P b 2 b',
        'P m a 2','P 2 m b','P c 2 m','P m 2 a','P b m 2','P 2 c m',
        'P c a 21','P 21 a b','P c 21 b','P b 21 a','P b c 21','P 21 c a',
        'P n c 2','P 2 n a','P b 2 n','P n 2 b','P c n 2','P 2 a n',
        'P m n 21','P 21 m n','P n 21 m','P m 21 n','P n m 21','P 21 n m',
        'P b a 2','P 2 c b','P c 2 a',
        'P n a 21','P 21 n b','P c 21 n','P n 21 a','P b n 21','P 21 c n',
        'P n n 2','P 2 n n','P n 2 n',
        'P m m m','P n n n',
        'P c c m','P m a a','P b m b',
        'P b a n','P n c b','P c n a',
        'P m m a','P b m m','P m c m','P m a m','P m m b','P c m m',
        'P n n a','P b n n','P n c n','P n a n','P n n b','P c n n',
        'P m n a','P b m n','P n c m','P m a n','P n m b','P c n m',
        'P c c a','P b a a','P b c b','P b a b','P c c b','P c a a',
        'P b a m','P m c b','P c m a',
        'P c c n','P n a a','P b n b',
        'P b c m','P m c a','P b m a','P c m b','P c a m','P m a b',
        'P n n m','P m n n','P n m n',
        'P m m n','P n m m','P m n m',
        'P b c n','P n c a','P b n a','P c n b','P c a n','P n a b',
        'P b c a','P c a b',
        'P n m a','P b n m','P m c n','P n a m','P m n b','P c m n',
        ),
    'Cmmm':('C 2 2 21','C 2 2 2','C m m 2',
        'C m c 21','C c m 21','C c c 2','C m 2 m','C 2 m m',
        'C m 2 a','C 2 m b','C c 2 m','C 2 c m','C c 2 a','C 2 c b',
        'C m c m','C c m m','C m c a','C c m b',
        'C m m m','C c c m','C m m a','C m m b','C c c a','C c c b',),
    'Ammm':('A 21 2 2','A 2 2 2','A 2 m m',
        'A 21 m a','A 21 a m','A 2 a a','A m 2 m','A m m 2',
        'A b m 2','A c 2 m','A m a 2','A m 2 a','A b a 2','A c 2 a',
        'A m m a','A m a m','A b m a','A c a m',
        'A m m m','A m a a','A b m m','A c m m','A c a a','A b a a',),
    'Bmmm':('B 2 21 2','B 2 2 2','B m 2 m',
        'B m 21 b','B b 21 m','B b 2 b','B m m 2','B 2 m m',
        'B 2 c m','B m a 2','B 2 m b','B b m 2','B 2 c b','B b a 2',
        'B b m m','B m m b','B b c m','B m a b',
        'B m m m','B b m b','B m a m','B m c m','B b a b','B b c b',),
    'Immm':('I 2 2 2','I 21 21 21',
        'I m m 2','I m 2 m','I 2 m m',
        'I b a 2','I 2 c b','I c 2 a',
        'I m a 2','I 2 m b','I c 2 m','I m 2 a','I b m 2','I 2 c m',
        'I m m m','I b a m','I m c b','I c m a',
        'I b c a','I c a b',
        'I m m a','I b m m ','I m c m','I m a m','I m m b','I c m m',),
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
    'Pm3m': ('P 2 3','P 21 3','P m 3','P m -3','P n 3','P n -3','P a 3','P a -3',
        'P 4 3 2','P 42 3 2','P 43 3 2','P 41 3 2','P -4 3 m','P -4 3 n',
        'P m 3 m','P m -3 m','P n 3 n','P n -3 n',
        'P m 3 n','P m -3 n','P n 3 m','P n -3 m',),
    'Im3m':('I 2 3','I 21 3','I m 3','I m -3','I a 3','I a -3', 'I 4 3 2','I 41 3 2',
        'I -4 3 m', 'I -4 3 d','I m -3 m','I m 3 m','I a 3 d','I a -3 d','I n 3 n','I n -3 n'),
    'Fm3m':('F 2 3','F m 3','F m -3','F d 3','F d -3',
        'F 4 3 2','F 41 3 2','F -4 3 m','F -4 3 c',
        'F m 3 m','F m -3 m','F m 3 c','F m -3 c',
        'F d 3 m','F d -3 m','F d 3 c','F d -3 c',),
}
sgequiv_2002_orthorhombic = {}
''' A dictionary of orthorhombic space groups that were renamed in the 2002 Volume A,
 along with the pre-2002 name. The e designates a double glide-plane
'''
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    sgequiv_2002_orthorhombic = {
        'AEM2':'A b m 2','B2EM':'B 2 c m','CM2E':'C m 2 a',
        'AE2M':'A c 2 m','BME2':'B m a 2','C2ME':'C 2 m b',
        'AEA2':'A b a 2','B2EB':'B 2 c b','CC2E':'C c 2 a',
        'AE2A':'A c 2 a','BBE2':'B b a 2','C2CE':'C 2 c b',
        'CMCE':'C m c a','AEMA':'A b m a','BBEM':'B b c m',
        'BMEB':'B m a b','CCME':'C c m b','AEAM':'A c a m',
        'CMME':'C m m a','AEMM':'A b m m','BMEM':'B m c m',
        'CCCE':'C c c a','AEAA':'A b a a','BBEB':'B b c b'}

spgHall = []
'''Hall space group symbols indexed by GSAS-II space group name. This is indexed 
by a key which is generated from the name used in GSAS-II with spaces removed 
and all letters in lowercase. The value associated with the key consists of a 
list of two values: the first is the Hall name and the second is a
Hermann-Mauguin name. Note that there may be several names used by 
GSAS-II that map to the same symmetry operators and the same Hall symbol
(for example "P 21" and "P 1 21 1" and "P 63/m" and "P 63/m H"). There
are also space group names that GSAS-II will accept that do not have Hall
symbols (such as "I -1" or "F 21 21 21"). 
'''
# generated from https://cci.lbl.gov/sginfo/hall_symbols.html by Sydney R. Hall & Ralf W. Grosse-Kunstleve
if 'sphinx' not in sys.modules: # skip if generating sphinx docs to unclutter 
    spgHall = {
        'p1':('P 1', 'P 1', 1),   # P 1
        'p-1':('-P 1', 'P -1', 2),   # P -1
        'p121':('P 2y', 'P 1 2 1', 3),   # P 1 2 1
        'p2':('P 2y', 'P 1 2 1', 3),   # P 2
        'p112':('P 2', 'P 1 1 2', 3),   # P 1 1 2
        'p211':('P 2x', 'P 2 1 1', 3),   # P 2 1 1
        'p1211':('P 2yb', 'P 1 21 1', 4),   # P 1 21 1
        'p21':('P 2yb', 'P 1 21 1', 4),   # P 21
        'p1121':('P 2c', 'P 1 1 21', 4),   # P 1 1 21
        'p2111':('P 2xa', 'P 21 1 1', 4),   # P 21 1 1
        'c121':('C 2y', 'C 1 2 1', 5),   # C 1 2 1
        'c2':('C 2y', 'C 1 2 1', 5),   # C 2
        'a121':('A 2y', 'A 1 2 1', 5),   # A 1 2 1
        'a2':('A 2y', 'A 1 2 1', 5),   # A 2
        'i121':('I 2y', 'I 1 2 1', 5),   # I 1 2 1
        'i2':('I 2y', 'I 1 2 1', 5),   # I 2
        'a112':('A 2', 'A 1 1 2', 5),   # A 1 1 2
        'b112':('B 2', 'B 1 1 2', 5),   # B 1 1 2
        'i112':('I 2', 'I 1 1 2', 5),   # I 1 1 2
        'b211':('B 2x', 'B 2 1 1', 5),   # B 2 1 1
        'c211':('C 2x', 'C 2 1 1', 5),   # C 2 1 1
        'i211':('I 2x', 'I 2 1 1', 5),   # I 2 1 1
        'p1m1':('P -2y', 'P 1 m 1', 6),   # P 1 m 1
        'pm':('P -2y', 'P 1 m 1', 6),   # P m
        'p11m':('P -2', 'P 1 1 m', 6),   # P 1 1 m
        'pm11':('P -2x', 'P m 1 1', 6),   # P m 1 1
        'p1c1':('P -2yc', 'P 1 c 1', 7),   # P 1 c 1
        'pc':('P -2yc', 'P 1 c 1', 7),   # P c
        'p1n1':('P -2yac', 'P 1 n 1', 7),   # P 1 n 1
        'pn':('P -2yac', 'P 1 n 1', 7),   # P n
        'p1a1':('P -2ya', 'P 1 a 1', 7),   # P 1 a 1
        'pa':('P -2ya', 'P 1 a 1', 7),   # P a
        'p11a':('P -2a', 'P 1 1 a', 7),   # P 1 1 a
        'p11n':('P -2ab', 'P 1 1 n', 7),   # P 1 1 n
        'p11b':('P -2b', 'P 1 1 b', 7),   # P 1 1 b
        'pb11':('P -2xb', 'P b 1 1', 7),   # P b 1 1
        'pn11':('P -2xbc', 'P n 1 1', 7),   # P n 1 1
        'pc11':('P -2xc', 'P c 1 1', 7),   # P c 1 1
        'c1m1':('C -2y', 'C 1 m 1', 8),   # C 1 m 1
        'cm':('C -2y', 'C 1 m 1', 8),   # C m
        'a1m1':('A -2y', 'A 1 m 1', 8),   # A 1 m 1
        'am':('A -2y', 'A 1 m 1', 8),   # A m
        'i1m1':('I -2y', 'I 1 m 1', 8),   # I 1 m 1
        'im':('I -2y', 'I 1 m 1', 8),   # I m
        'a11m':('A -2', 'A 1 1 m', 8),   # A 1 1 m
        'b11m':('B -2', 'B 1 1 m', 8),   # B 1 1 m
        'i11m':('I -2', 'I 1 1 m', 8),   # I 1 1 m
        'bm11':('B -2x', 'B m 1 1', 8),   # B m 1 1
        'cm11':('C -2x', 'C m 1 1', 8),   # C m 1 1
        'im11':('I -2x', 'I m 1 1', 8),   # I m 1 1
        'c1c1':('C -2yc', 'C 1 c 1', 9),   # C 1 c 1
        'cc':('C -2yc', 'C 1 c 1', 9),   # C c
        'a1n1':('A -2yac', 'A 1 n 1', 9),   # A 1 n 1
        'an':('A -2yac', 'A 1 n 1', 9),   # A n
        'i1a1':('I -2ya', 'I 1 a 1', 9),   # I 1 a 1
        'ia':('I -2ya', 'I 1 a 1', 9),   # I a
        'a1a1':('A -2ya', 'A 1 a 1', 9),   # A 1 a 1
        'aa':('A -2ya', 'A 1 a 1', 9),   # A a
        'c1n1':('C -2ybc', 'C 1 n 1', 9),   # C 1 n 1
        'cn':('C -2ybc', 'C 1 n 1', 9),   # C n
        'i1c1':('I -2yc', 'I 1 c 1', 9),   # I 1 c 1
        'ic':('I -2yc', 'I 1 c 1', 9),   # I c
        'a11a':('A -2a', 'A 1 1 a', 9),   # A 1 1 a
        'b11n':('B -2bc', 'B 1 1 n', 9),   # B 1 1 n
        'i11b':('I -2b', 'I 1 1 b', 9),   # I 1 1 b
        'b11b':('B -2b', 'B 1 1 b', 9),   # B 1 1 b
        'a11n':('A -2ac', 'A 1 1 n', 9),   # A 1 1 n
        'i11a':('I -2a', 'I 1 1 a', 9),   # I 1 1 a
        'bb11':('B -2xb', 'B b 1 1', 9),   # B b 1 1
        'cn11':('C -2xbc', 'C n 1 1', 9),   # C n 1 1
        'ic11':('I -2xc', 'I c 1 1', 9),   # I c 1 1
        'cc11':('C -2xc', 'C c 1 1', 9),   # C c 1 1
        'bn11':('B -2xbc', 'B n 1 1', 9),   # B n 1 1
        'ib11':('I -2xb', 'I b 1 1', 9),   # I b 1 1
        'p12/m1':('-P 2y', 'P 1 2/m 1', 10),   # P 1 2/m 1
        'p2/m':('-P 2y', 'P 1 2/m 1', 10),   # P 2/m
        'p112/m':('-P 2', 'P 1 1 2/m', 10),   # P 1 1 2/m
        'p2/m11':('-P 2x', 'P 2/m 1 1', 10),   # P 2/m 1 1
        'p121/m1':('-P 2yb', 'P 1 21/m 1', 11),   # P 1 21/m 1
        'p21/m':('-P 2yb', 'P 1 21/m 1', 11),   # P 21/m
        'p1121/m':('-P 2c', 'P 1 1 21/m', 11),   # P 1 1 21/m
        'p21/m11':('-P 2xa', 'P 21/m 1 1', 11),   # P 21/m 1 1
        'c12/m1':('-C 2y', 'C 1 2/m 1', 12),   # C 1 2/m 1
        'c2/m':('-C 2y', 'C 1 2/m 1', 12),   # C 2/m
        'a12/m1':('-A 2y', 'A 1 2/m 1', 12),   # A 1 2/m 1
        'a2/m':('-A 2y', 'A 1 2/m 1', 12),   # A 2/m
        'i12/m1':('-I 2y', 'I 1 2/m 1', 12),   # I 1 2/m 1
        'i2/m':('-I 2y', 'I 1 2/m 1', 12),   # I 2/m
        'a112/m':('-A 2', 'A 1 1 2/m', 12),   # A 1 1 2/m
        'b112/m':('-B 2', 'B 1 1 2/m', 12),   # B 1 1 2/m
        'i112/m':('-I 2', 'I 1 1 2/m', 12),   # I 1 1 2/m
        'b2/m11':('-B 2x', 'B 2/m 1 1', 12),   # B 2/m 1 1
        'c2/m11':('-C 2x', 'C 2/m 1 1', 12),   # C 2/m 1 1
        'i2/m11':('-I 2x', 'I 2/m 1 1', 12),   # I 2/m 1 1
        'p12/c1':('-P 2yc', 'P 1 2/c 1', 13),   # P 1 2/c 1
        'p2/c':('-P 2yc', 'P 1 2/c 1', 13),   # P 2/c
        'p12/n1':('-P 2yac', 'P 1 2/n 1', 13),   # P 1 2/n 1
        'p2/n':('-P 2yac', 'P 1 2/n 1', 13),   # P 2/n
        'p12/a1':('-P 2ya', 'P 1 2/a 1', 13),   # P 1 2/a 1
        'p2/a':('-P 2ya', 'P 1 2/a 1', 13),   # P 2/a
        'p112/a':('-P 2a', 'P 1 1 2/a', 13),   # P 1 1 2/a
        'p112/n':('-P 2ab', 'P 1 1 2/n', 13),   # P 1 1 2/n
        'p112/b':('-P 2b', 'P 1 1 2/b', 13),   # P 1 1 2/b
        'p2/b11':('-P 2xb', 'P 2/b 1 1', 13),   # P 2/b 1 1
        'p2/n11':('-P 2xbc', 'P 2/n 1 1', 13),   # P 2/n 1 1
        'p2/c11':('-P 2xc', 'P 2/c 1 1', 13),   # P 2/c 1 1
        'p121/c1':('-P 2ybc', 'P 1 21/c 1', 14),   # P 1 21/c 1
        'p21/c':('-P 2ybc', 'P 1 21/c 1', 14),   # P 21/c
        'p121/n1':('-P 2yn', 'P 1 21/n 1', 14),   # P 1 21/n 1
        'p21/n':('-P 2yn', 'P 1 21/n 1', 14),   # P 21/n
        'p121/a1':('-P 2yab', 'P 1 21/a 1', 14),   # P 1 21/a 1
        'p21/a':('-P 2yab', 'P 1 21/a 1', 14),   # P 21/a
        'p1121/a':('-P 2ac', 'P 1 1 21/a', 14),   # P 1 1 21/a
        'p1121/n':('-P 2n', 'P 1 1 21/n', 14),   # P 1 1 21/n
        'p1121/b':('-P 2bc', 'P 1 1 21/b', 14),   # P 1 1 21/b
        'p21/b11':('-P 2xab', 'P 21/b 1 1', 14),   # P 21/b 1 1
        'p21/n11':('-P 2xn', 'P 21/n 1 1', 14),   # P 21/n 1 1
        'p21/c11':('-P 2xac', 'P 21/c 1 1', 14),   # P 21/c 1 1
        'c12/c1':('-C 2yc', 'C 1 2/c 1', 15),   # C 1 2/c 1
        'c2/c':('-C 2yc', 'C 1 2/c 1', 15),   # C 2/c
        'a12/n1':('-A 2yac', 'A 1 2/n 1', 15),   # A 1 2/n 1
        'a2/n':('-A 2yac', 'A 1 2/n 1', 15),   # A 2/n
        'i12/a1':('-I 2ya', 'I 1 2/a 1', 15),   # I 1 2/a 1
        'i2/a':('-I 2ya', 'I 1 2/a 1', 15),   # I 2/a
        'a12/a1':('-A 2ya', 'A 1 2/a 1', 15),   # A 1 2/a 1
        'a2/a':('-A 2ya', 'A 1 2/a 1', 15),   # A 2/a
        'c12/n1':('-C 2ybc', 'C 1 2/n 1', 15),   # C 1 2/n 1
        'c2/n':('-C 2ybc', 'C 1 2/n 1', 15),   # C 2/n
        'i12/c1':('-I 2yc', 'I 1 2/c 1', 15),   # I 1 2/c 1
        'i2/c':('-I 2yc', 'I 1 2/c 1', 15),   # I 2/c
        'a112/a':('-A 2a', 'A 1 1 2/a', 15),   # A 1 1 2/a
        'b112/n':('-B 2bc', 'B 1 1 2/n', 15),   # B 1 1 2/n
        'i112/b':('-I 2b', 'I 1 1 2/b', 15),   # I 1 1 2/b
        'b112/b':('-B 2b', 'B 1 1 2/b', 15),   # B 1 1 2/b
        'a112/n':('-A 2ac', 'A 1 1 2/n', 15),   # A 1 1 2/n
        'i112/a':('-I 2a', 'I 1 1 2/a', 15),   # I 1 1 2/a
        'b2/b11':('-B 2xb', 'B 2/b 1 1', 15),   # B 2/b 1 1
        'c2/n11':('-C 2xbc', 'C 2/n 1 1', 15),   # C 2/n 1 1
        'i2/c11':('-I 2xc', 'I 2/c 1 1', 15),   # I 2/c 1 1
        'c2/c11':('-C 2xc', 'C 2/c 1 1', 15),   # C 2/c 1 1
        'b2/n11':('-B 2xbc', 'B 2/n 1 1', 15),   # B 2/n 1 1
        'i2/b11':('-I 2xb', 'I 2/b 1 1', 15),   # I 2/b 1 1
        'p222':('P 2 2', 'P 2 2 2', 16),   # P 2 2 2
        'p2221':('P 2c 2', 'P 2 2 21', 17),   # P 2 2 21
        'p2122':('P 2a 2a', 'P 21 2 2', 17),   # P 21 2 2
        'p2212':('P 2 2b', 'P 2 21 2', 17),   # P 2 21 2
        'p21212':('P 2 2ab', 'P 21 21 2', 18),   # P 21 21 2
        'p22121':('P 2bc 2', 'P 2 21 21', 18),   # P 2 21 21
        'p21221':('P 2ac 2ac', 'P 21 2 21', 18),   # P 21 2 21
        'p212121':('P 2ac 2ab', 'P 21 21 21', 19),   # P 21 21 21
        'c2221':('C 2c 2', 'C 2 2 21', 20),   # C 2 2 21
        'a2122':('A 2a 2a', 'A 21 2 2', 20),   # A 21 2 2
        'b2212':('B 2 2b', 'B 2 21 2', 20),   # B 2 21 2
        'c222':('C 2 2', 'C 2 2 2', 21),   # C 2 2 2
        'a222':('A 2 2', 'A 2 2 2', 21),   # A 2 2 2
        'b222':('B 2 2', 'B 2 2 2', 21),   # B 2 2 2
        'f222':('F 2 2', 'F 2 2 2', 22),   # F 2 2 2
        'i222':('I 2 2', 'I 2 2 2', 23),   # I 2 2 2
        'i212121':('I 2b 2c', 'I 21 21 21', 24),   # I 21 21 21
        'pmm2':('P 2 -2', 'P m m 2', 25),   # P m m 2
        'p2mm':('P -2 2', 'P 2 m m', 25),   # P 2 m m
        'pm2m':('P -2 -2', 'P m 2 m', 25),   # P m 2 m
        'pmc21':('P 2c -2', 'P m c 21', 26),   # P m c 21
        'pcm21':('P 2c -2c', 'P c m 21', 26),   # P c m 21
        'p21ma':('P -2a 2a', 'P 21 m a', 26),   # P 21 m a
        'p21am':('P -2 2a', 'P 21 a m', 26),   # P 21 a m
        'pb21m':('P -2 -2b', 'P b 21 m', 26),   # P b 21 m
        'pm21b':('P -2b -2', 'P m 21 b', 26),   # P m 21 b
        'pcc2':('P 2 -2c', 'P c c 2', 27),   # P c c 2
        'p2aa':('P -2a 2', 'P 2 a a', 27),   # P 2 a a
        'pb2b':('P -2b -2b', 'P b 2 b', 27),   # P b 2 b
        'pma2':('P 2 -2a', 'P m a 2', 28),   # P m a 2
        'pbm2':('P 2 -2b', 'P b m 2', 28),   # P b m 2
        'p2mb':('P -2b 2', 'P 2 m b', 28),   # P 2 m b
        'p2cm':('P -2c 2', 'P 2 c m', 28),   # P 2 c m
        'pc2m':('P -2c -2c', 'P c 2 m', 28),   # P c 2 m
        'pm2a':('P -2a -2a', 'P m 2 a', 28),   # P m 2 a
        'pca21':('P 2c -2ac', 'P c a 21', 29),   # P c a 21
        'pbc21':('P 2c -2b', 'P b c 21', 29),   # P b c 21
        'p21ab':('P -2b 2a', 'P 21 a b', 29),   # P 21 a b
        'p21ca':('P -2ac 2a', 'P 21 c a', 29),   # P 21 c a
        'pc21b':('P -2bc -2c', 'P c 21 b', 29),   # P c 21 b
        'pb21a':('P -2a -2ab', 'P b 21 a', 29),   # P b 21 a
        'pnc2':('P 2 -2bc', 'P n c 2', 30),   # P n c 2
        'pcn2':('P 2 -2ac', 'P c n 2', 30),   # P c n 2
        'p2na':('P -2ac 2', 'P 2 n a', 30),   # P 2 n a
        'p2an':('P -2ab 2', 'P 2 a n', 30),   # P 2 a n
        'pb2n':('P -2ab -2ab', 'P b 2 n', 30),   # P b 2 n
        'pn2b':('P -2bc -2bc', 'P n 2 b', 30),   # P n 2 b
        'pmn21':('P 2ac -2', 'P m n 21', 31),   # P m n 21
        'pnm21':('P 2bc -2bc', 'P n m 21', 31),   # P n m 21
        'p21mn':('P -2ab 2ab', 'P 21 m n', 31),   # P 21 m n
        'p21nm':('P -2 2ac', 'P 21 n m', 31),   # P 21 n m
        'pn21m':('P -2 -2bc', 'P n 21 m', 31),   # P n 21 m
        'pm21n':('P -2ab -2', 'P m 21 n', 31),   # P m 21 n
        'pba2':('P 2 -2ab', 'P b a 2', 32),   # P b a 2
        'p2cb':('P -2bc 2', 'P 2 c b', 32),   # P 2 c b
        'pc2a':('P -2ac -2ac', 'P c 2 a', 32),   # P c 2 a
        'pna21':('P 2c -2n', 'P n a 21', 33),   # P n a 21
        'pbn21':('P 2c -2ab', 'P b n 21', 33),   # P b n 21
        'p21nb':('P -2bc 2a', 'P 21 n b', 33),   # P 21 n b
        'p21cn':('P -2n 2a', 'P 21 c n', 33),   # P 21 c n
        'pc21n':('P -2n -2ac', 'P c 21 n', 33),   # P c 21 n
        'pn21a':('P -2ac -2n', 'P n 21 a', 33),   # P n 21 a
        'pnn2':('P 2 -2n', 'P n n 2', 34),   # P n n 2
        'p2nn':('P -2n 2', 'P 2 n n', 34),   # P 2 n n
        'pn2n':('P -2n -2n', 'P n 2 n', 34),   # P n 2 n
        'cmm2':('C 2 -2', 'C m m 2', 35),   # C m m 2
        'a2mm':('A -2 2', 'A 2 m m', 35),   # A 2 m m
        'bm2m':('B -2 -2', 'B m 2 m', 35),   # B m 2 m
        'cmc21':('C 2c -2', 'C m c 21', 36),   # C m c 21
        'ccm21':('C 2c -2c', 'C c m 21', 36),   # C c m 21
        'a21ma':('A -2a 2a', 'A 21 m a', 36),   # A 21 m a
        'a21am':('A -2 2a', 'A 21 a m', 36),   # A 21 a m
        'bb21m':('B -2 -2b', 'B b 21 m', 36),   # B b 21 m
        'bm21b':('B -2b -2', 'B m 21 b', 36),   # B m 21 b
        'ccc2':('C 2 -2c', 'C c c 2', 37),   # C c c 2
        'a2aa':('A -2a 2', 'A 2 a a', 37),   # A 2 a a
        'bb2b':('B -2b -2b', 'B b 2 b', 37),   # B b 2 b
        'amm2':('A 2 -2', 'A m m 2', 38),   # A m m 2
        'bmm2':('B 2 -2', 'B m m 2', 38),   # B m m 2
        'b2mm':('B -2 2', 'B 2 m m', 38),   # B 2 m m
        'c2mm':('C -2 2', 'C 2 m m', 38),   # C 2 m m
        'cm2m':('C -2 -2', 'C m 2 m', 38),   # C m 2 m
        'am2m':('A -2 -2', 'A m 2 m', 38),   # A m 2 m
        'abm2':('A 2 -2c', 'A b m 2', 39),   # A b m 2
        'bma2':('B 2 -2c', 'B m a 2', 39),   # B m a 2
        'b2cm':('B -2c 2', 'B 2 c m', 39),   # B 2 c m
        'c2mb':('C -2b 2', 'C 2 m b', 39),   # C 2 m b
        'cm2a':('C -2b -2b', 'C m 2 a', 39),   # C m 2 a
        'ac2m':('A -2c -2c', 'A c 2 m', 39),   # A c 2 m
        'ama2':('A 2 -2a', 'A m a 2', 40),   # A m a 2
        'bbm2':('B 2 -2b', 'B b m 2', 40),   # B b m 2
        'b2mb':('B -2b 2', 'B 2 m b', 40),   # B 2 m b
        'c2cm':('C -2c 2', 'C 2 c m', 40),   # C 2 c m
        'cc2m':('C -2c -2c', 'C c 2 m', 40),   # C c 2 m
        'am2a':('A -2a -2a', 'A m 2 a', 40),   # A m 2 a
        'aba2':('A 2 -2ac', 'A b a 2', 41),   # A b a 2
        'bba2':('B 2 -2bc', 'B b a 2', 41),   # B b a 2
        'b2cb':('B -2bc 2', 'B 2 c b', 41),   # B 2 c b
        'c2cb':('C -2bc 2', 'C 2 c b', 41),   # C 2 c b
        'cc2a':('C -2bc -2bc', 'C c 2 a', 41),   # C c 2 a
        'ac2a':('A -2ac -2ac', 'A c 2 a', 41),   # A c 2 a
        'fmm2':('F 2 -2', 'F m m 2', 42),   # F m m 2
        'f2mm':('F -2 2', 'F 2 m m', 42),   # F 2 m m
        'fm2m':('F -2 -2', 'F m 2 m', 42),   # F m 2 m
        'fdd2':('F 2 -2d', 'F d d 2', 43),   # F d d 2
        'f2dd':('F -2d 2', 'F 2 d d', 43),   # F 2 d d
        'fd2d':('F -2d -2d', 'F d 2 d', 43),   # F d 2 d
        'imm2':('I 2 -2', 'I m m 2', 44),   # I m m 2
        'i2mm':('I -2 2', 'I 2 m m', 44),   # I 2 m m
        'im2m':('I -2 -2', 'I m 2 m', 44),   # I m 2 m
        'iba2':('I 2 -2c', 'I b a 2', 45),   # I b a 2
        'i2cb':('I -2a 2', 'I 2 c b', 45),   # I 2 c b
        'ic2a':('I -2b -2b', 'I c 2 a', 45),   # I c 2 a
        'ima2':('I 2 -2a', 'I m a 2', 46),   # I m a 2
        'ibm2':('I 2 -2b', 'I b m 2', 46),   # I b m 2
        'i2mb':('I -2b 2', 'I 2 m b', 46),   # I 2 m b
        'i2cm':('I -2c 2', 'I 2 c m', 46),   # I 2 c m
        'ic2m':('I -2c -2c', 'I c 2 m', 46),   # I c 2 m
        'im2a':('I -2a -2a', 'I m 2 a', 46),   # I m 2 a
        'pmmm':('-P 2 2', 'P m m m', 47),   # P m m m
        'pnnn':('-P 2ab 2bc', 'P n n n', 48),   # P n n n
        'pccm':('-P 2 2c', 'P c c m', 49),   # P c c m
        'pmaa':('-P 2a 2', 'P m a a', 49),   # P m a a
        'pbmb':('-P 2b 2b', 'P b m b', 49),   # P b m b
        'pban':('-P 2ab 2b', 'P b a n', 50),   # P b a n
        'pncb':('-P 2b 2bc', 'P n c b', 50),   # P n c b
        'pcna':('-P 2a 2c', 'P c n a', 50),   # P c n a
        'pmma':('-P 2a 2a', 'P m m a', 51),   # P m m a
        'pmmb':('-P 2b 2', 'P m m b', 51),   # P m m b
        'pbmm':('-P 2 2b', 'P b m m', 51),   # P b m m
        'pcmm':('-P 2c 2c', 'P c m m', 51),   # P c m m
        'pmcm':('-P 2c 2', 'P m c m', 51),   # P m c m
        'pmam':('-P 2 2a', 'P m a m', 51),   # P m a m
        'pnna':('-P 2a 2bc', 'P n n a', 52),   # P n n a
        'pnnb':('-P 2b 2n', 'P n n b', 52),   # P n n b
        'pbnn':('-P 2n 2b', 'P b n n', 52),   # P b n n
        'pcnn':('-P 2ab 2c', 'P c n n', 52),   # P c n n
        'pncn':('-P 2ab 2n', 'P n c n', 52),   # P n c n
        'pnan':('-P 2n 2bc', 'P n a n', 52),   # P n a n
        'pmna':('-P 2ac 2', 'P m n a', 53),   # P m n a
        'pnmb':('-P 2bc 2bc', 'P n m b', 53),   # P n m b
        'pbmn':('-P 2ab 2ab', 'P b m n', 53),   # P b m n
        'pcnm':('-P 2 2ac', 'P c n m', 53),   # P c n m
        'pncm':('-P 2 2bc', 'P n c m', 53),   # P n c m
        'pman':('-P 2ab 2', 'P m a n', 53),   # P m a n
        'pcca':('-P 2a 2ac', 'P c c a', 54),   # P c c a
        'pccb':('-P 2b 2c', 'P c c b', 54),   # P c c b
        'pbaa':('-P 2a 2b', 'P b a a', 54),   # P b a a
        'pcaa':('-P 2ac 2c', 'P c a a', 54),   # P c a a
        'pbcb':('-P 2bc 2b', 'P b c b', 54),   # P b c b
        'pbab':('-P 2b 2ab', 'P b a b', 54),   # P b a b
        'pbam':('-P 2 2ab', 'P b a m', 55),   # P b a m
        'pmcb':('-P 2bc 2', 'P m c b', 55),   # P m c b
        'pcma':('-P 2ac 2ac', 'P c m a', 55),   # P c m a
        'pccn':('-P 2ab 2ac', 'P c c n', 56),   # P c c n
        'pnaa':('-P 2ac 2bc', 'P n a a', 56),   # P n a a
        'pbnb':('-P 2bc 2ab', 'P b n b', 56),   # P b n b
        'pbcm':('-P 2c 2b', 'P b c m', 57),   # P b c m
        'pcam':('-P 2c 2ac', 'P c a m', 57),   # P c a m
        'pmca':('-P 2ac 2a', 'P m c a', 57),   # P m c a
        'pmab':('-P 2b 2a', 'P m a b', 57),   # P m a b
        'pbma':('-P 2a 2ab', 'P b m a', 57),   # P b m a
        'pcmb':('-P 2bc 2c', 'P c m b', 57),   # P c m b
        'pnnm':('-P 2 2n', 'P n n m', 58),   # P n n m
        'pmnn':('-P 2n 2', 'P m n n', 58),   # P m n n
        'pnmn':('-P 2n 2n', 'P n m n', 58),   # P n m n
        'pmmn':('-P 2ab 2a', 'P m m n', 59),   # P m m n
        'pnmm':('-P 2c 2bc', 'P n m m', 59),   # P n m m
        'pmnm':('-P 2c 2a', 'P m n m', 59),   # P m n m
        'pbcn':('-P 2n 2ab', 'P b c n', 60),   # P b c n
        'pcan':('-P 2n 2c', 'P c a n', 60),   # P c a n
        'pnca':('-P 2a 2n', 'P n c a', 60),   # P n c a
        'pnab':('-P 2bc 2n', 'P n a b', 60),   # P n a b
        'pbna':('-P 2ac 2b', 'P b n a', 60),   # P b n a
        'pcnb':('-P 2b 2ac', 'P c n b', 60),   # P c n b
        'pbca':('-P 2ac 2ab', 'P b c a', 61),   # P b c a
        'pcab':('-P 2bc 2ac', 'P c a b', 61),   # P c a b
        'pnma':('-P 2ac 2n', 'P n m a', 62),   # P n m a
        'pmnb':('-P 2bc 2a', 'P m n b', 62),   # P m n b
        'pbnm':('-P 2c 2ab', 'P b n m', 62),   # P b n m
        'pcmn':('-P 2n 2ac', 'P c m n', 62),   # P c m n
        'pmcn':('-P 2n 2a', 'P m c n', 62),   # P m c n
        'pnam':('-P 2c 2n', 'P n a m', 62),   # P n a m
        'cmcm':('-C 2c 2', 'C m c m', 63),   # C m c m
        'ccmm':('-C 2c 2c', 'C c m m', 63),   # C c m m
        'amma':('-A 2a 2a', 'A m m a', 63),   # A m m a
        'amam':('-A 2 2a', 'A m a m', 63),   # A m a m
        'bbmm':('-B 2 2b', 'B b m m', 63),   # B b m m
        'bmmb':('-B 2b 2', 'B m m b', 63),   # B m m b
        'cmca':('-C 2bc 2', 'C m c a', 64),   # C m c a
        'ccmb':('-C 2bc 2bc', 'C c m b', 64),   # C c m b
        'abma':('-A 2ac 2ac', 'A b m a', 64),   # A b m a
        'acam':('-A 2 2ac', 'A c a m', 64),   # A c a m
        'bbcm':('-B 2 2bc', 'B b c m', 64),   # B b c m
        'bmab':('-B 2bc 2', 'B m a b', 64),   # B m a b
        'cmmm':('-C 2 2', 'C m m m', 65),   # C m m m
        'ammm':('-A 2 2', 'A m m m', 65),   # A m m m
        'bmmm':('-B 2 2', 'B m m m', 65),   # B m m m
        'cccm':('-C 2 2c', 'C c c m', 66),   # C c c m
        'amaa':('-A 2a 2', 'A m a a', 66),   # A m a a
        'bbmb':('-B 2b 2b', 'B b m b', 66),   # B b m b
        'cmma':('-C 2b 2', 'C m m a', 67),   # C m m a
        'cmmb':('-C 2b 2b', 'C m m b', 67),   # C m m b
        'abmm':('-A 2c 2c', 'A b m m', 67),   # A b m m
        'acmm':('-A 2 2c', 'A c m m', 67),   # A c m m
        'bmcm':('-B 2 2c', 'B m c m', 67),   # B m c m
        'bmam':('-B 2c 2', 'B m a m', 67),   # B m a m
        'ccca':('-C 2b 2bc', 'C c c a', 68),   # C c c a
        'cccb':('-C 2b 2c', 'C c c b', 68),   # C c c b
        'abaa':('-A 2a 2c', 'A b a a', 68),   # A b a a
        'acaa':('-A 2ac 2c', 'A c a a', 68),   # A c a a
        'bbcb':('-B 2bc 2b', 'B b c b', 68),   # B b c b
        'bbab':('-B 2b 2bc', 'B b a b', 68),   # B b a b
        'fmmm':('-F 2 2', 'F m m m', 69),   # F m m m
        'fddd':('-F 2uv 2vw', 'F d d d', 70),   # F d d d
        'immm':('-I 2 2', 'I m m m', 71),   # I m m m
        'ibam':('-I 2 2c', 'I b a m', 72),   # I b a m
        'imcb':('-I 2a 2', 'I m c b', 72),   # I m c b
        'icma':('-I 2b 2b', 'I c m a', 72),   # I c m a
        'ibca':('-I 2b 2c', 'I b c a', 73),   # I b c a
        'icab':('-I 2a 2b', 'I c a b', 73),   # I c a b
        'imma':('-I 2b 2', 'I m m a', 74),   # I m m a
        'immb':('-I 2a 2a', 'I m m b', 74),   # I m m b
        'ibmm':('-I 2c 2c', 'I b m m', 74),   # I b m m
        'icmm':('-I 2 2b', 'I c m m', 74),   # I c m m
        'imcm':('-I 2 2a', 'I m c m', 74),   # I m c m
        'imam':('-I 2c 2', 'I m a m', 74),   # I m a m
        'p4':('P 4', 'P 4', 75),   # P 4
        'p41':('P 4w', 'P 41', 76),   # P 41
        'p42':('P 4c', 'P 42', 77),   # P 42
        'p43':('P 4cw', 'P 43', 78),   # P 43
        'i4':('I 4', 'I 4', 79),   # I 4
        'i41':('I 4bw', 'I 41', 80),   # I 41
        'p-4':('P -4', 'P -4', 81),   # P -4
        'i-4':('I -4', 'I -4', 82),   # I -4
        'p4/m':('-P 4', 'P 4/m', 83),   # P 4/m
        'p42/m':('-P 4c', 'P 42/m', 84),   # P 42/m
        'p4/n':('-P 4a', 'P 4/n', 85),   # P 4/n
        'p42/n':('-P 4bc', 'P 42/n', 86),   # P 42/n
        'i4/m':('-I 4', 'I 4/m', 87),   # I 4/m
        'i41/a':('-I 4ad', 'I 41/a', 88),   # I 41/a
        'p422':('P 4 2', 'P 4 2 2', 89),   # P 4 2 2
        'p4212':('P 4ab 2ab', 'P 42 1 2', 90),   # P 42 1 2
        'p4122':('P 4w 2c', 'P 41 2 2', 91),   # P 41 2 2
        'p41212':('P 4abw 2nw', 'P 41 21 2', 92),   # P 41 21 2
        'p4222':('P 4c 2', 'P 42 2 2', 93),   # P 42 2 2
        'p42212':('P 4n 2n', 'P 42 21 2', 94),   # P 42 21 2
        'p4322':('P 4cw 2c', 'P 43 2 2', 95),   # P 43 2 2
        'p43212':('P 4nw 2abw', 'P 43 21 2', 96),   # P 43 21 2
        'i422':('I 4 2', 'I 4 2 2', 97),   # I 4 2 2
        'i4122':('I 4bw 2bw', 'I 41 2 2', 98),   # I 41 2 2
        'p4mm':('P 4 -2', 'P 4 m m', 99),   # P 4 m m
        'p4bm':('P 4 -2ab', 'P 4 b m', 100),   # P 4 b m
        'p42cm':('P 4c -2c', 'P 42 c m', 101),   # P 42 c m
        'p42nm':('P 4n -2n', 'P 42 n m', 102),   # P 42 n m
        'p4cc':('P 4 -2c', 'P 4 c c', 103),   # P 4 c c
        'p4nc':('P 4 -2n', 'P 4 n c', 104),   # P 4 n c
        'p42mc':('P 4c -2', 'P 42 m c', 105),   # P 42 m c
        'p42bc':('P 4c -2ab', 'P 42 b c', 106),   # P 42 b c
        'i4mm':('I 4 -2', 'I 4 m m', 107),   # I 4 m m
        'i4cm':('I 4 -2c', 'I 4 c m', 108),   # I 4 c m
        'i41md':('I 4bw -2', 'I 41 m d', 109),   # I 41 m d
        'i41cd':('I 4bw -2c', 'I 41 c d', 110),   # I 41 c d
        'p-42m':('P -4 2', 'P -4 2 m', 111),   # P -4 2 m
        'p-42c':('P -4 2c', 'P -4 2 c', 112),   # P -4 2 c
        'p-421m':('P -4 2ab', 'P -4 21 m', 113),   # P -4 21 m
        'p-421c':('P -4 2n', 'P -4 21 c', 114),   # P -4 21 c
        'p-4m2':('P -4 -2', 'P -4 m 2', 115),   # P -4 m 2
        'p-4c2':('P -4 -2c', 'P -4 c 2', 116),   # P -4 c 2
        'p-4b2':('P -4 -2ab', 'P -4 b 2', 117),   # P -4 b 2
        'p-4n2':('P -4 -2n', 'P -4 n 2', 118),   # P -4 n 2
        'i-4m2':('I -4 -2', 'I -4 m 2', 119),   # I -4 m 2
        'i-4c2':('I -4 -2c', 'I -4 c 2', 120),   # I -4 c 2
        'i-42m':('I -4 2', 'I -4 2 m', 121),   # I -4 2 m
        'i-42d':('I -4 2bw', 'I -4 2 d', 122),   # I -4 2 d
        'p4/mmm':('-P 4 2', 'P 4/m m m', 123),   # P 4/m m m
        'p4/mcc':('-P 4 2c', 'P 4/m c c', 124),   # P 4/m c c
        'p4/nbm':('-P 4a 2b', 'P 4/n b m', 125),   # P 4/n b m
        'p4/nnc':('-P 4a 2bc', 'P 4/n n c', 126),   # P 4/n n c
        'p4/mbm':('-P 4 2ab', 'P 4/m b m', 127),   # P 4/m b m
        'p4/mnc':('-P 4 2n', 'P 4/m n c', 128),   # P 4/m n c
        'p4/nmm':('-P 4a 2a', 'P 4/n m m', 129),   # P 4/n m m
        'p4/ncc':('-P 4a 2ac', 'P 4/n c c', 130),   # P 4/n c c
        'p42/mmc':('-P 4c 2', 'P 42/m m c', 131),   # P 42/m m c
        'p42/mcm':('-P 4c 2c', 'P 42/m c m', 132),   # P 42/m c m
        'p42/nbc':('-P 4ac 2b', 'P 42/n b c', 133),   # P 42/n b c
        'p42/nnm':('-P 4ac 2bc', 'P 42/n n m', 134),   # P 42/n n m
        'p42/mbc':('-P 4c 2ab', 'P 42/m b c', 135),   # P 42/m b c
        'p42/mnm':('-P 4n 2n', 'P 42/m n m', 136),   # P 42/m n m
        'p42/nmc':('-P 4ac 2a', 'P 42/n m c', 137),   # P 42/n m c
        'p42/ncm':('-P 4ac 2ac', 'P 42/n c m', 138),   # P 42/n c m
        'i4/mmm':('-I 4 2', 'I 4/m m m', 139),   # I 4/m m m
        'i4/mcm':('-I 4 2c', 'I 4/m c m', 140),   # I 4/m c m
        'i41/amd':('-I 4bd 2', 'I 41/a m d', 141),   # I 41/a m d
        'i41/acd':('-I 4bd 2c', 'I 41/a c d', 142),   # I 41/a c d
        'p3':('P 3', 'P 3', 143),   # P 3
        'p31':('P 31', 'P 31', 144),   # P 31
        'p32':('P 32', 'P 32', 145),   # P 32
        'r3h':('R 3', 'R 3 H', 146),   # R 3 H
        'r3':('R 3', 'R 3 H', 146),   # R 3
        'r3r':('P 3*', 'R 3 R', 146),   # R 3 R
        'p-3':('-P 3', 'P -3', 147),   # P -3
        'r-3h':('-R 3', 'R -3 H', 148),   # R -3 H
        'r-3':('-R 3', 'R -3 H', 148),   # R -3
        'r-3r':('-P 3*', 'R -3 R', 148),   # R -3 R
        'p312':('P 3 2', 'P 3 1 2', 149),   # P 3 1 2
        'p321':('P 3 2"', 'P 3 2 1', 150),   # P 3 2 1
        'p3112':('P 31 2c (0 0 1)', 'P 31 1 2', 151),   # P 31 1 2
        'p3121':('P 31 2"', 'P 31 2 1', 152),   # P 31 2 1
        'p3212':('P 32 2c (0 0 -1)', 'P 32 1 2', 153),   # P 32 1 2
        'p3221':('P 32 2"', 'P 32 2 1', 154),   # P 32 2 1
        'r32h':('R 3 2"', 'R 32 H', 155),   # R 32 H
        'r32':('R 3 2"', 'R 32 H', 155),   # R 32
        'r32r':('P 3* 2', 'R 32 R', 155),   # R 32 R
        'p3m1':('P 3 -2"', 'P 3 m 1', 156),   # P 3 m 1
        'p31m':('P 3 -2', 'P 3 1 m', 157),   # P 3 1 m
        'p3c1':('P 3 -2"c', 'P 3 c 1', 158),   # P 3 c 1
        'p31c':('P 3 -2c', 'P 3 1 c', 159),   # P 3 1 c
        'r3mh':('R 3 -2"', 'R 3 m H', 160),   # R 3 m H
        'r3m':('R 3 -2"', 'R 3 m H', 160),   # R 3 m
        'r3mr':('P 3* -2', 'R 3 m R', 160),   # R 3 m R
        'r3ch':('R 3 -2"c', 'R 3 c H', 161),   # R 3 c H
        'r3c':('R 3 -2"c', 'R 3 c H', 161),   # R 3 c
        'r3cr':('P 3* -2n', 'R 3 c R', 161),   # R 3 c R
        'p-31m':('-P 3 2', 'P -3 1 m', 162),   # P -3 1 m
        'p-31c':('-P 3 2c', 'P -3 1 c', 163),   # P -3 1 c
        'p-3m1':('-P 3 2"', 'P -3 m 1', 164),   # P -3 m 1
        'p-3c1':('-P 3 2"c', 'P -3 c 1', 165),   # P -3 c 1
        'r-3mh':('-R 3 2"', 'R -3 m H', 166),   # R -3 m H
        'r-3m':('-R 3 2"', 'R -3 m H', 166),   # R -3 m
        'r-3mr':('-P 3* 2', 'R -3 m R', 166),   # R -3 m R
        'r-3ch':('-R 3 2"c', 'R -3 c H', 167),   # R -3 c H
        'r-3c':('-R 3 2"c', 'R -3 c H', 167),   # R -3 c
        'r-3cr':('-P 3* 2n', 'R -3 c R', 167),   # R -3 c R
        'p6':('P 6', 'P 6', 168),   # P 6
        'p61':('P 61', 'P 61', 169),   # P 61
        'p65':('P 65', 'P 65', 170),   # P 65
        'p62':('P 62', 'P 62', 171),   # P 62
        'p64':('P 64', 'P 64', 172),   # P 64
        'p63':('P 6c', 'P 63', 173),   # P 63
        'p-6':('P -6', 'P -6', 174),   # P -6
        'p6/m':('-P 6', 'P 6/m', 175),   # P 6/m
        'p63/m':('-P 6c', 'P 63/m', 176),   # P 63/m
        'p622':('P 6 2', 'P 6 2 2', 177),   # P 6 2 2
        'p6122':('P 61 2 (0 0 -1)', 'P 61 2 2', 178),   # P 61 2 2
        'p6522':('P 65 2 (0 0 1)', 'P 65 2 2', 179),   # P 65 2 2
        'p6222':('P 62 2c (0 0 1)', 'P 62 2 2', 180),   # P 62 2 2
        'p6422':('P 64 2c (0 0 -1)', 'P 64 2 2', 181),   # P 64 2 2
        'p6322':('P 6c 2c', 'P 63 2 2', 182),   # P 63 2 2
        'p6mm':('P 6 -2', 'P 6 m m', 183),   # P 6 m m
        'p6cc':('P 6 -2c', 'P 6 c c', 184),   # P 6 c c
        'p63cm':('P 6c -2', 'P 63 c m', 185),   # P 63 c m
        'p63mc':('P 6c -2c', 'P 63 m c', 186),   # P 63 m c
        'p-6m2':('P -6 2', 'P -6 m 2', 187),   # P -6 m 2
        'p-6c2':('P -6c 2', 'P -6 c 2', 188),   # P -6 c 2
        'p-62m':('P -6 -2', 'P -6 2 m', 189),   # P -6 2 m
        'p-62c':('P -6c -2c', 'P -6 2 c', 190),   # P -6 2 c
        'p6/mmm':('-P 6 2', 'P 6/m m m', 191),   # P 6/m m m
        'p6/mcc':('-P 6 2c', 'P 6/m c c', 192),   # P 6/m c c
        'p63/mcm':('-P 6c 2', 'P 63/m c m', 193),   # P 63/m c m
        'p63/mmc':('-P 6c 2c', 'P 63/m m c', 194),   # P 63/m m c
        'p23':('P 2 2 3', 'P 2 3', 195),   # P 2 3
        'f23':('F 2 2 3', 'F 2 3', 196),   # F 2 3
        'i23':('I 2 2 3', 'I 2 3', 197),   # I 2 3
        'p213':('P 2ac 2ab 3', 'P 21 3', 198),   # P 21 3
        'i213':('I 2b 2c 3', 'I 21 3', 199),   # I 21 3
        'pm-3':('-P 2 2 3', 'P m -3', 200),   # P m -3
        'pn-3':('-P 2ab 2bc 3', 'P n -3', 201),   # P n -3
        'fm-3':('-F 2 2 3', 'F m -3', 202),   # F m -3
        'fd-3':('-F 2uv 2vw 3', 'F d -3', 203),   # F d -3
        'im-3':('-I 2 2 3', 'I m -3', 204),   # I m -3
        'pa-3':('-P 2ac 2ab 3', 'P a -3', 205),   # P a -3
        'ia-3':('-I 2b 2c 3', 'I a -3', 206),   # I a -3
        'p432':('P 4 2 3', 'P 4 3 2', 207),   # P 4 3 2
        'p4232':('P 4n 2 3', 'P 42 3 2', 208),   # P 42 3 2
        'f432':('F 4 2 3', 'F 4 3 2', 209),   # F 4 3 2
        'f4132':('F 4d 2 3', 'F 41 3 2', 210),   # F 41 3 2
        'i432':('I 4 2 3', 'I 4 3 2', 211),   # I 4 3 2
        'p4332':('P 4acd 2ab 3', 'P 43 3 2', 212),   # P 43 3 2
        'p4132':('P 4bd 2ab 3', 'P 41 3 2', 213),   # P 41 3 2
        'i4132':('I 4bd 2c 3', 'I 41 3 2', 214),   # I 41 3 2
        'p-43m':('P -4 2 3', 'P -4 3 m', 215),   # P -4 3 m
        'f-43m':('F -4 2 3', 'F -4 3 m', 216),   # F -4 3 m
        'i-43m':('I -4 2 3', 'I -4 3 m', 217),   # I -4 3 m
        'p-43n':('P -4n 2 3', 'P -4 3 n', 218),   # P -4 3 n
        'f-43c':('F -4c 2 3', 'F -4 3 c', 219),   # F -4 3 c
        'i-43d':('I -4bd 2c 3', 'I -4 3 d', 220),   # I -4 3 d
        'pm-3m':('-P 4 2 3', 'P m -3 m', 221),   # P m -3 m
        'pn-3n':('-P 4a 2bc 3', 'P n -3 n', 222),   # P n -3 n
        'pm-3n':('-P 4n 2 3', 'P m -3 n', 223),   # P m -3 n
        'pn-3m':('-P 4bc 2bc 3', 'P n -3 m', 224),   # P n -3 m
        'fm-3m':('-F 4 2 3', 'F m -3 m', 225),   # F m -3 m
        'fm-3c':('-F 4c 2 3', 'F m -3 c', 226),   # F m -3 c
        'fd-3m':('-F 4vw 2vw 3', 'F d -3 m', 227),   # F d -3 m
        'fd-3c':('-F 4cvw 2vw 3', 'F d -3 c', 228),   # F d -3 c
        'im-3m':('-I 4 2 3', 'I m -3 m', 229),   # I m -3 m
        'ia-3d':('-I 4bd 2c 3', 'I a -3 d', 230),   # I a -3 d
        }
        
#'A few non-standard space groups for test use'
nonstandard_sglist = ('P 21 1 1','P 1 21 1','P 1 1 21','R 3 r','R 3 2 h', 
                      'R -3 r', 'R 3 2 r','R 3 m h', 'R 3 m r',
                      'R 3 c r','R -3 c r','R -3 m r',),

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
    assert (MoveToUnitCell([1,2,3])[0] == [0,0,0]).all, msg
    assert (MoveToUnitCell([2,-1,-2])[0] == [0,0,0]).all, msg
    assert abs(MoveToUnitCell(np.array([-.1]))[0]-0.9)[0] < 1e-6, msg
    assert abs(MoveToUnitCell(np.array([.1]))[0]-0.1)[0] < 1e-6, msg
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
#        #print result[1]['SpGrp']
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
        if debug: print (spc['SpGrp'])
        if debug: print (spc['SGCen'])
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
                    noff = MoveToUnitCell(noff)[0]
                    mult = tuple((op*inv).ravel().tolist())
                    if debug: print ("\n%s: %s + %s" % (spcname,mult,noff))
                    for refop in cctbx:
                        if debug: print (refop)
                        # check the transform
                        if refop[:9] != mult: continue
                        if debug: print ("mult match")
                        # check the translation
                        reftrans = list(refop[-3:])
                        reftrans = MoveToUnitCell(reftrans)[0]
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
            symb, m, n, od = SytSym(t[0],S)
            if symb.strip() != t[2].strip() or m != t[1]:
                print (spc,t[0],m,n,symb,t[2],od)
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
    print ("OK")
