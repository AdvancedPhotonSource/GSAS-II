# -*- coding: utf-8 -*-
"GSASII - Space group interpretion routines"
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
import math
import sys
import os.path as ospath

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import pyspg

def SpcGroup(SGSymbol):
    '''
    Determines cell and symmetry information from a short H-M space group name
    input:
        SGSymbol - space group symbol (string) with spaces between axial fields
    returns:
        SGError = 0 for no errors; >0 for errors (see SGErrors below for details)
        SGData - dictionary with entries:
             'SpGrp': space group symbol slightly cleaned up
             'Laue':  one of '-1','2/m','mmm','4/m','4/mmm','3R','3mR','3',
                      '3m1','31m','6/m','6/mmm','m3','m3m'
             'SGInv': boolean; True if centrosymmetric, False if not
             'SGLatt': one of 'P','A','B','C','I','F','R'
             'SGUniq': one of 'a','b','c' if monoclinic, '' otherwise
             'SGCen': cell centering vectors [0,0,0] at least
             'SGOps': symmetry operations as [M,T] so that M*x+T = x'
             'SGSys': one of 'triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','trigonal','hexagonal','cubic'
             'SGPolax': one of '','x','y','x y','z','x z','y z','xyz','111' for arbitrary axes
       '''
    LaueSym = ('-1','2/m','mmm','4/m','4/mmm','3R','3mR','3','3m1','31m','6/m','6/mmm','m3','m3m')
    LattSym = ('P','A','B','C','I','F','R')
    UniqSym = ('','','a','b','c','',)
    SysSym = ('triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','trigonal','hexagonal','cubic')
    SGData = {}
    SGData['SpGrp'] = SGSymbol.strip().lower().capitalize()
    SGInfo = pyspg.sgforpy(SGSymbol)
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
    return SGInfo[8],SGData

def SGErrors(IErr):
    '''
    Interprets the error message code from SpcGroup. Used in SpaceGroup.
    input:
        SGError - from SpcGroup
    returns:
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
    
def SGPrint(SGData):
    '''
    Print the output of SpcGroup in a nicely formatted way. Used in SpaceGroup
    input:
        SGData - from SpcGroup
    returns:
        SGText - list of strings with the space group details
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
    SGText.append(' Multiplicity of a general site is '+str(Mult))
    SGText.append(' The Laue symmetry is '+SGData['SGLaue'])
    if SGData['SGUniq'] in ['a','b','c']:
        SGText.append(' The unique monoclinic axis is '+SGData['SGUniq'])
    if SGData['SGInv']:
        SGText.append(' The inversion center is located at 0,0,0')
    if SGData['SGPolax']:
        SGText.append(' The location of the origin is arbitrary in '+SGData['SGPolax'])
    SGText.append('\n'+' The equivalent positions are:')
    if SGData['SGLatt'] != 'P':
        SGText.append('\n ('+Latt2text(SGData['SGLatt'])+')+\n')
    Ncol = 2
    line = ' '
    col = 0
    for iop,[M,T] in enumerate(SGData['SGOps']):
        OPtxt = MT2text(M,T)
        Fld = '(%2i) '%(iop+1)+OPtxt+'\t'
        line += Fld
        if '/' not in Fld:
            line += '\t'
        col += 1
        if col == Ncol:
            SGText.append(line)        
            line = ' '
            col = 0
    SGText.append(line)        
    return SGText
    
def MT2text(M,T):
    #From space group matrix/translation operator returns text version
    XYZ = ('-Z ','-Y ','-X ','X-Y','ERR','Y-X',' X ',' Y ',' Z ','+X ','+Y ','+Z ')
    TRA = ('   ','ERR','1/6','1/4','1/3','ERR','1/2','ERR','2/3','3/4','5/6','ERR')
    Fld = ''
    for j in range(3):
        IJ = int(round(2*M[j][0]+3*M[j][1]+4*M[j][2]+4))%12
        IK = int(round(T[j]*12))%12
        if IK > 0 and IJ > 4: IJ += 3
        Fld += TRA[IK]+XYZ[IJ]
        if j != 2: Fld += ','
    return Fld
    
def Latt2text(Latt):
    #From lattice type ('P',A', etc.) returns ';' delimited cell centering vectors
    lattTxt = {'A':'0,0,0; 0,1/2,1/2','B':'0,0,0; 1/2,0,1/2',
        'C':'0,0,0; 1/2,1/2,0','I':'0,0,0; 1/2,1/2,1/2',
        'F':'0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0',
        'R':'0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3','P':'0,0,0'}
    return lattTxt[Latt]    
        
def SpaceGroup(SGSymbol):
    '''
    Print the output of SpcGroup in a nicely formatted way. 
    input: 
        SGSymbol - space group symbol (string) with spaces between axial fields
    returns:
        nothing
    '''
    E,A = SpcGroup(SGSymbol)
    if E > 0:
        print SGErrors(E)
        return
    for l in SGPrint(A):
        print l

def MoveToUnitCell(xyz):
    '''
    Translates a set of coordinates so that all values are >=0 and < 1 
    input:
        xyz - a list or numpy array of fractional coordinates
    returns: 
        XYZ - numpy array of new coordinates inside 0-1
    '''
    XYZ = np.zeros(3)
    for i,x in enumerate(xyz):
        XYZ[i] = (x-int(x))%1.0
    return XYZ
        
def Opposite(XYZ,toler=0.0002):
    '''
    Gives opposite corner, edge or face of unit cell for position within tolerance. 
        Result may be just outside the cell within tolerance 
    input:
        XYZ: 0 >= np.array[x,y,z] > 1 as by MoveToUnitCell
        toler: unit cell fraction tolerance making opposite
    returns:
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
    input:  
        XYZ an array, tuple or list containing 3 elements: x, y & z
        SGData, from SpcGroup
        All  = True return all equivalent positions including duplicates
             = False return only unique positions
        Uij  = [U11,U22,U33,U12,U13,U23] or [] if no Uij
        Move = True move generated atom positions to be inside cell
             = False do not move atoms       
    return: [[XYZEquiv],Idup,[UijEquiv]]
        [XYZEquiv] is list of equivalent positions (XYZ is first entry)
        Idup = [-][C]SS where SS is the symmetry operator number (1-24), C (if not 0,0,0)
        is centering operator number (1-4) and - is for inversion
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
                newX = MoveToUnitCell(XT)
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
    input:
        HKL - [h,k,l]
        SGData - space group data obtained from SpcGroup
    returns:
        iabsnt = True is reflection is forbidden by symmetry
        mulp = reflection multiplicity including Friedel pairs
        Uniq = numpy array of equivalent hkl in descending order of h,k,l
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
    KNsym = {
        '0'         :'    1   ','1'         :'   -1   ','64'        :'  2(100)','32'        :'  m(100)',
        '97'        :'2/m(100)','16'        :'  2(010)','8'         :'  m(010)','25'        :'2/m(010)',
        '2'         :'  2(001)','4'         :'  m(001)','7'         :'2/m(001)','134217728' :'  2(011)',
        '67108864'  :'  m(011)','201326593' :'2/m(011)','2097152'   :'  2(0+-)','1048576'   :'  m(0+-)',
        '3145729'   :'2/m(0+-)','8388608'   :'  2(101)','4194304'   :'  m(101)','12582913'  :'2/m(101)',
        '524288'    :'  2(+0-)','262144'    :'  m(+0-)','796433'    :'2/m(+0-)','1024'      :'  2(110)',
        '512'       :'  m(110)','1537'      :'2/m(110)','256'       :'  2(+-0)','128'       :'  m(+-0)',
        '385'       :'2/m(+-0)','76'        :'mm2(100)','52'        :'mm2(010)','42'        :'mm2(001)',
        '135266336' :'mm2(011)','69206048'  :'mm2(0+-)','8650760'   :'mm2(101)','4718600'   :'mm2(+0-)',
        '1156'      :'mm2(110)','772'       :'mm2(+-0)','82'        :'  222   ','136314944' :'222(100)',
        '8912912'   :'222(010)','1282'      :'222(001)','127'       :'  mmm   ','204472417' :'mmm(100)',
        '13369369'  :'mmm(010)','1927'      :'mmm(001)','33554496'  :'  4(100)','16777280'  :' -4(100)',
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
    NXUPQsym = {
        '    1   ':(28,29,28,28),'   -1   ':( 1,29,28, 0),'  2(100)':(12,18,12,25),'  m(100)':(25,18,12,25),
        '2/m(100)':( 1,18, 0,-1),'  2(010)':(13,17,13,24),'  m(010)':(24,17,13,24),'2/m(010)':( 1,17, 0,-1),
        '  2(001)':(14,16,14,23),'  m(001)':(23,16,14,23),'2/m(001)':( 1,16, 0,-1),'  2(011)':(10,23,10,22),
        '  m(011)':(22,23,10,22),'2/m(011)':( 1,23, 0,-1),'  2(0+-)':(11,24,11,21),'  m(0+-)':(21,24,11,21),
        '2/m(0+-)':( 1,24, 0,-1),'  2(101)':( 8,21, 8,20),'  m(101)':(20,21, 8,20),'2/m(101)':( 1,21, 0,-1),
        '  2(+0-)':( 9,22, 9,19),'  m(+0-)':(19,22, 9,19),'2/m(+0-)':( 1,22, 0,-1),'  2(110)':( 6,19, 6,18),
        '  m(110)':(18,19, 6,18),'2/m(110)':( 1,19, 0,-1),'  2(+-0)':( 7,20, 7,17),'  m(+-0)':(17,20, 7,17),
        '2/m(+-0)':( 1,20, 0,-1),'mm2(100)':(12,10, 0,-1),'mm2(010)':(13,10, 0,-1),'mm2(001)':(14,10, 0,-1),
        'mm2(011)':(10,13, 0,-1),'mm2(0+-)':(11,13, 0,-1),'mm2(101)':( 8,12, 0,-1),'mm2(+0-)':( 9,12, 0,-1),
        'mm2(110)':( 6,11, 0,-1),'mm2(+-0)':( 7,11, 0,-1),'  222   ':( 1,10, 0,-1),'222(100)':( 1,13, 0,-1),
        '222(010)':( 1,12, 0,-1),'222(001)':( 1,11, 0,-1),'  mmm   ':( 1,10, 0,-1),'mmm(100)':( 1,13, 0,-1),
        'mmm(010)':( 1,12, 0,-1),'mmm(001)':( 1,11, 0,-1),'  4(100)':(12, 4,12, 0),' -4(100)':( 1, 4,12, 0),
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
    CSxinel = [[],                         # 0th empty - indices are Fortran style
        [[0,0,0],[ 0.0, 0.0, 0.0]],      #  0  0  0
        [[1,1,1],[ 1.0, 1.0, 1.0]],      #  X  X  X
        [[1,1,1],[ 1.0, 1.0,-1.0]],      #  X  X -X
        [[1,1,1],[ 1.0,-1.0, 1.0]],      #  X -X  X
        [[1,1,1],[ 1.0,-1.0,-1.0]],      # -X  X  X
        [[1,1,0],[ 1.0, 1.0, 0.0]],      #  X  X  0
        [[1,1,0],[ 1.0,-1.0, 0.0]],      #  X -X  0
        [[1,0,1],[ 1.0, 0.0, 1.0]],      #  X  0  X
        [[1,0,1],[ 1.0, 0.0,-1.0]],      #  X  0 -X
        [[0,1,1],[ 0.0, 1.0, 1.0]],      #  0  Y  Y
        [[0,1,1],[ 0.0, 1.0,-1.0]],      #  0  Y -Y
        [[1,0,0],[ 1.0, 0.0, 0.0]],      #  X  0  0
        [[0,1,0],[ 0.0, 1.0, 0.0]],      #  0  Y  0
        [[0,0,1],[ 0.0, 0.0, 1.0]],      #  0  0  Z
        [[1,1,0],[ 1.0, 2.0, 0.0]],      #  X 2X  0
        [[1,1,0],[ 2.0, 1.0, 0.0]],      # 2X  X  0
        [[1,1,2],[ 1.0, 1.0, 1.0]],      #  X  X  Z
        [[1,1,2],[ 1.0,-1.0, 1.0]],      #  X -X  Z
        [[1,2,1],[ 1.0, 1.0, 1.0]],      #  X  Y  X
        [[1,2,1],[ 1.0, 1.0,-1.0]],      #  X  Y -X
        [[1,2,2],[ 1.0, 1.0, 1.0]],      #  X  Y  Y
        [[1,2,2],[ 1.0, 1.0,-1.0]],      #  X  Y -Y
        [[1,2,0],[ 1.0, 1.0, 0.0]],      #  X  Y  0
        [[1,0,2],[ 1.0, 0.0, 1.0]],      #  X  0  Z
        [[0,1,2],[ 0.0, 1.0, 1.0]],      #  0  Y  Z
        [[1,1,2],[ 1.0, 2.0, 1.0]],      #  X 2X  Z
        [[1,1,2],[ 2.0, 1.0, 1.0]],      # 2X  X  Z
        [[1,2,3],[ 1.0, 1.0, 1.0]],      #  X  Y  Z
        ]
    indx = GetNXUPQsym(siteSym)
    return CSxinel[indx[0]]
    
def GetCSuinel(siteSym):
    # returns Uij terms, multipliers, GUI flags & Uiso2Uij multipliers
    CSuinel = [[],                                             # 0th empty - indices are Fortran style
        [[1,1,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,0,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  A  0  0  0
        [[1,1,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,0,1,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  C  0  0  0
        [[1,2,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  A  0  0  0
        [[1,2,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,0,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  B  0  0  0
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  A  D  D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  A  D -D -D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  A  D -D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,0,0,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  A  D  D -D
        [[1,1,2,1,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,0,1,0,0,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #  A  A  C A/2 0  0
        [[1,2,3,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  C  0  0  0
        [[1,1,2,3,0,0],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0],[1,0,1,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  C  D  0  0
        [[1,2,1,0,3,0],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],[1,1,0,0,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  A  0  E  0
        [[1,2,2,0,0,3],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],[1,1,0,0,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  B  0  0  F
        [[1,2,3,2,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.5,0.0]],    #  A  B  C B/2 0  0
        [[1,2,3,1,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],[1,1,1,0,0,0],[1.0,1.0,1.0,0.0,0.5,0.0]],    #  A  B  C A/2 0  0
        [[1,2,3,4,0,0],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0],[1,1,1,1,0,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  C  D  0  0
        [[1,2,3,0,4,0],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],[1,1,1,0,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  C  0  E  0
        [[1,2,3,0,0,4],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  C  0  0  F
        [[1,1,2,3,4,4],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,0,1,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  C  D  E -E
        [[1,1,2,3,4,4],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,0,1,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  A  C  D  E  E
        [[1,2,1,3,4,3],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0],[1,1,0,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  A  D  E -D
        [[1,2,1,3,4,3],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,0,1,1,0],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  A  D  E  D
        [[1,2,2,3,3,4],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0],[1,1,0,1,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  B  D -D  F
        [[1,2,2,3,3,4],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,0,1,0,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  B  D  D  F
        [[1,2,3,2,4,4],[ 1.0, 1.0, 1.0, 0.5, 0.5, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #  A  B  C B/2 F/2 F
        [[1,2,3,1,0,4],[ 1.0, 1.0, 1.0, 0.5, 0.0, 1.0],[1,1,1,0,0,1],[1.0,1.0,1.0,0.5,0.0,0.0]],    #  A  B  C A/2  0  F
        [[1,2,3,2,4,0],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.0],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #  A  B  C B/2  E  0
        [[1,2,3,1,4,4],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.5],[1,1,1,0,1,0],[1.0,1.0,1.0,0.5,0.0,0.0]],    #  A  B  C A/2  E E/2
        [[1,2,3,4,5,6],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],[1,1,1,1,1,1],[1.0,1.0,1.0,0.0,0.0,0.0]],    #  A  B  C  D  E   F
        ]
    indx = GetNXUPQsym(siteSym)
    return CSuinel[indx[1]]
    
def MustrainNames(SGData):
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
    #this is to convert isotropic mustrain to generalized Shkls - doesn't work just now
    import GSASIIlattice as G2lat
    from scipy.optimize import fmin
    A = G2lat.cell2AB(cell)[0]
    def minMus(Shkl,H,muiso,SGData,A):
        U = np.inner(A.T,H)
        S = np.array(MustrainCoeff(H.T,SGData))
        sum = np.sqrt(np.sum(np.multiply(S,Shkl)))
        return abs(muiso-sum*H)
    laue = SGData['SGLaue']
    if laue in ['m3','m3m']:
        H = [[1,0,0],[1,1,0]]
        S0 = [0.01,0.01]
    elif laue in ['6/m','6/mmm','3m1']:
        H = [[1,0,0],[0,0,1],[1,0,1]]
        S0 = [0.01,0.01,0.01]
    elif laue in ['31m','3']:
        H = [[1,0,0],[0,0,1],[1,0,1],[1,1,1]]
        S0 = [0.01,0.01,0.01,0.01]
    elif laue in ['3R','3mR']:
        H = [[1,0,0],[1,1,0],[1,0,1],[1,1,1]]
        S0 = [0.01,0.01,0.01,0.01]
    elif laue in ['4/m','4/mmm']:
        H = [[1,0,0],[0,0,1],[1,1,0],[1,0,1]]
        S0 = [0.01,0.01,0.01,0.01]
    elif laue in ['mmm']:
        H = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]]
        S0 = [0.01,0.01,0.01,0.01,0.01,0.01]
    elif laue in ['2/m']:
        H = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]]
        if uniq == 'a':
            H.append([0,1,-1])
            H.append([0,-2,1])
        elif uniq == 'b':
            H.append([1,0,-1])
            H.append([-2,0,1])
        elif uniq == 'c':
            H.append([1,-1,0])
            H.append([-2,1,0])
        H.append([1,1,1])
        S0 = [9*[0.01,]]
    else:
        H = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],
            [-1,1,0],[1,0,-1],[0,-1,1],[1,-2,0],[-2,0,1],[0,1,-2],
            [1,-1,1],[-1, 1, 1],[1,-1,1]]
        S0 = [15*[0.01,]]
    H = np.array(H)
    S0 = np.array(S0)
    return fmin(minMus,S0,(H,muiso,SGData,A))
       
def SytSym(XYZ,SGData):
    '''
    Generates the number of equivalent positions and a site symmetry code for a specified coordinate and space group
    input:  
       XYZ: an array, tuple or list containing 3 elements: x, y & z
       SGData: from SpcGroup
    Returns a two element tuple:
       The 1st element is a code for the site symmetry (see GetKNsym)
       The 2nd element is the site multiplicity
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
    ''' Under development
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
    ''' Find A*B where A & B are in strings '-' + '100*c+n' + '+ijk'
    where '-' indicates inversion, c(>0) is the cell centering operator, 
    n is operator number from SgOps and ijk are unit cell translations (each may be <0).
    Should return resultant string - C.
        SGData - dictionary using entries:
             'SGCen': cell centering vectors [0,0,0] at least
             'SGOps': symmetry operations as [M,T] so that M*x+T = x'
    '''
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
    return [U[0][0],U[1][1],U[2][2],U[0][1],U[0][2],U[1][2]]
    
def Uij2U(Uij):
    #returns the thermal motion tensor U from Uij as numpy array
    return np.array([[Uij[0],Uij[3],Uij[4]],[Uij[3],Uij[1],Uij[5]],[Uij[4],Uij[5],Uij[2]]])
    
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
        'C m 2 a','C 2 m b','C 2 c m','C c 2 m','C 2 c m','C c 2 m',
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
        'F d d 2','F d 2 d','F d 2 d',),
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
    'Pm3m': ('P 2 3','P 21 3','P m -3','P n -3','P a -3','P 4 3 2','P 42 3 2',
        'P 43 3 2','P 41 3 2','P -4 3 m','P -4 3 n','P m -3 m','P n -3 n',
        'P m -3 n','P n -3 m',),
    'Im3m':('I 2 3','I 21 3','I m -3','I a -3', 'I 4 3 2','I 41 3 2',
        'I -4 3 m', 'I -4 3 d','I m -3 m','I a -3 d',),
    'Fm3m':('F 2 3','F m -3','F d -3','F 4 3 2','F 41 3 2','F -4 3 m',
        'F -4 3 c','F m -3 m','F m -3 c','F d -3 m','F d -3 c',),
}
'A few non-standard space groups for test use'
nonstandard_sglist = ('P 21 1 1','P 1 21 1','P 1 1 21','R 3 r','R 3 2 h', 
                      'R -3 r', 'R 3 2 r','R 3 m h', 'R 3 m r',
                      'R 3 c r','R -3 c r','R -3 m r',),
'''A list of orthorhombic space groups that were renamed in the 2002 Volume A,
along with the pre-2002 name. The e designates a double glide-plane'''
sgequiv_2002_orthorhombic= (('A e m 2', 'A b m 2',),
                            ('A e a 2', 'A b a 2',),
                            ('C m c e', 'C m c a',),
                            ('C m m e', 'C m m a',),
                            ('C c c e', 'C c c a'),)
'''Use the space groups types in this order to list the symbols in the 
order they are listed in the International Tables, vol. A'''
symtypelist = ('triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 
               'trigonal', 'hexagonal', 'cubic')

# self-test materials follow. Requires files in directory testinp
def test0():
    '''test #0: exercise MoveToUnitCell'''
    msg = "MoveToUnitCell failed"
    assert (MoveToUnitCell([1,2,3]) == [0,0,0]).all, msg
    assert (MoveToUnitCell([2,-1,-2]) == [0,0,0]).all, msg
    assert abs(MoveToUnitCell(np.array([-.1]))[0]-0.9) < 1e-6, msg
    assert abs(MoveToUnitCell(np.array([.1]))[0]-0.1) < 1e-6, msg

def test1():
    ''' test #1: SpcGroup and SGPrint against previous results'''
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
        msg = msg0 + " in list lengths"
        assert len(keys) == len(refdict.keys()), msg
        for key in keys:
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
        assert reflist == SGPrint(result[1]), 'SGPrint ' +msg
    for spc in spctestinp.SGdat:
        CompareSpcGroup(spc, 0, spctestinp.SGdat[spc], spctestinp.SGlist[spc] )

def test2():
    ''' test #2: SpcGroup against cctbx (sgtbx) computations'''
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

def test3(): 
    ''' test #3: exercise SytSym (includes GetOprPtrName, GenAtom, GetKNsym)
     for selected space groups against info in IT Volume A '''
    def ExerciseSiteSym (spc, crdlist):
        'compare site symmetries and multiplicities for a specified space group'
        msg = "failed on site sym test for %s" % spc
        (E,S) = SpcGroup(spc)
        assert not E, msg
        for t in crdlist:
            symb, m = SytSym(t[0],S)
            if symb.strip() != t[2].strip() or m != t[1]:
                print spc,t[0],m,symb
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
            ((0.0,.31,0.25),4,'2(010)'),
            ((0.25,.25,0.5),4,'-1'),
            ((0,0.5,0),4,'-1'),
            ])
    ExerciseSiteSym('p 2 2 2',[
            ((0.13,0.22,0.31),4,'1'),
            ((0,0.5,.31),2,'2(001)'),
            ((0.5,.31,0.5),2,'2(010)'),
            ((.11,0,0),2,'2(100)'),
            ((0,0.5,0),1,'222'),
            ])
    ExerciseSiteSym('p 4/n',[
            ((0.13,0.22,0.31),8,'1'),
            ((0.25,0.75,.31),4,'2(001)'),
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
            ((0.11,0,0.25),24,'2(100)'),
            ((0.11,0.11,0.11),16,'3(111)'),
            ((0,0,0),8,'-3(111)'),
            ])

if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()
    print "OK"
