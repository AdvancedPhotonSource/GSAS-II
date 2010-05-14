"GSASII - Space group interpretion routines"

import numpy as np
import sys
import os.path as ospath

import GSASIIpath
import pyspg

def SpcGroup(SGSymbol):
    '''
   Determines cell and symmetry information from a short H-M space group name
   input: space group symbol (string) with spaces between axial fields
   returns [SGError,SGData]
       SGError = 0 for no errors; >0 for errors (see SGErrors below for details)
       returns dictionary SGData with entries:
         'SpGrp': space group symbol slightly cleaned up
         'Laue':  one of '-1','2/m','mmm','4/m','4/mmm','3R','3mR','3',
                  '3m1','31m','6/m','6/mmm','m3','m3m'
         'SGInv': boolean; True if centrosymmetric, False if not
         'SGLatt': one of 'P','A','B','C','I','F','R'
         'SGUniq': one of 'a','b','c' if monoclinic, '' otherwise
         'SGCen': cell centering vectors [0,0,0] at least
         'SGOps': symmetry operations as [M,T] so that M*x+T = x'
         'SGSys': one of 'triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','trigonal','hexagonal','cubic'
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
    return SGInfo[8],SGData

def SGErrors(IErr):
    '''Interprets the error message code from SpcGroup. Used in SpaceGroup.
    input:  SGError, from SpcGroup
    returns a string with the error message or "Unknown error"
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
    
def SGPrint(SGData):
    '''
    Print the output of SpcGroup in a nicely formatted way. Used in SpaceGroup
    input:  SGData, from SpcGroup
    returns a list of strings with the space group details
    '''
    XYZ = ('-Z ','-Y ','-X ','X-Y','ERR','Y-X',' X ',' Y ',' Z ','+X ','+Y ','+Z ')
    TRA = ('   ','ERR','1/6','1/4','1/3','ERR','1/2','ERR','2/3','3/4','5/6','ERR')
    POL = (' ','x','y','x y','z','x z','y z','xyz','111')
    Mult = len(SGData['SGCen'])*len(SGData['SGOps'])*(int(SGData['SGInv'])+1)
    NP = [1,2,4]
    NPZ = [0,1]
    for M,T in SGData['SGOps']:
        for i in range(3):
            if M[i][i] <= 0.: NP[i] = 0
        if M[0][2] > 0: NPZ[0] = 8
        if M[1][2] > 0: NPZ[1] = 0
    NPol = (NP[0]+NP[1]+NP[2]+NPZ[0]*NPZ[1])*(1-int(SGData['SGInv']))
    SGText = []
    SGText.append('Space Group '+SGData['SpGrp'])
    CentStr = 'centrosymmetric'
    if not SGData['SGInv']:
        CentStr = 'non'+CentStr
    if SGData['SGLatt'] in 'ABCIFR':
        SGText.append('The lattice is '+CentStr+' '+SGData['SGLatt']+'-centered '+SGData['SGSys'].lower())
    else:
        SGText.append('The lattice is '+CentStr+' '+'primitive '+SGData['SGSys'].lower())        
    SGText.append('Multiplicity of a general site is '+str(Mult))
    SGText.append('The Laue symmetry is '+SGData['SGLaue'])
    if SGData['SGUniq'] in ['a','b','c']:
        SGText.append('The unique monoclinic axis is '+SGData['SGUniq'])
    if SGData['SGInv']:
        SGText.append('The inversion center is located at 0,0,0')
    if NPol:
        SGText.append('The location of the origin is arbitrary in '+POL[NPol])
    SGText.append('\n'+'The equivalent positions are:')
    if SGData['SGLatt'] in 'A':
        SGText.append('\n'+'    (0,0,0; 0,1/2,1/2)+')
    elif SGData['SGLatt'] in 'B':
        SGText.append('\n'+'    (0,0,0; 1/2,0,1/2)+')
    elif SGData['SGLatt'] in 'C':
        SGText.append('\n'+'    (0,0,0; 1/2,1/2,0)+')
    elif SGData['SGLatt'] in 'I':
        SGText.append('\n'+'    (0,0,0; 1/2,1/2,1/2)+')
    elif SGData['SGLatt'] in 'F':
        SGText.append('\n'+'    (0,0,0; 0,1/2,1/2; 1/2,0,1/2; 1/2,1/2,0)+')
    elif SGData['SGLatt'] in 'R':
        SGText.append('\n'+'    (0,0,0; 1/3,2/3,2/3; 2/3,1/3,1/3)+')
    if SGData['SGLaue'] in ['-1','2/m','mmm','4/m','4/mmm']:
        Ncol = 2
    else:
        Ncol = 3
    line = ''
    for iop,[M,T] in enumerate(SGData['SGOps']):
        if iop % Ncol == 0:
            SGText.append(line)        
            line = ''
        Fld = '(%2i) ' % (iop+1)
        for j in range(3):
            IJ = int(round(2*M[j][0]+3*M[j][1]+4*M[j][2]+4)) % 12
            IK = int(round(T[j]*12)) % 12
            if IK > 0 and IJ > 4: IJ += 3
            Fld += TRA[IK]+XYZ[IJ]
            if j != 2: Fld += ','
        line += Fld
    SGText.append(line)        
    return SGText
    
def SpaceGroup(SgSym):
    '''
    Print the output of SpcGroup in a nicely formatted way. 
      input: space group symbol (string) with spaces between axial fields
      returns nothing
    '''
    E,A = SpcGroup(SgSym)
    if E > 0:
        print SGErrors(E)
        return
    for l in SGPrint(A):
        print l

def MoveToUnitCell(XYZ):
    '''
    Translates a set of coordinates so that all values are >=0 and < 1 
      input: a list or numpy array of any length. Note that the object is modified  in place.
      output: none
    '''
    for i,x in enumerate(XYZ):
        x = ((x % 1.0)+1.0) % 1.0
        if x > 0.9999: x = 0.0
        XYZ[i] = x
        
def GenAtom(XYZ,SGData,ifAll=False):
    '''
    Generates the equivalent positions for a specified coordinate and space group
    input:  
       XYZ an array, tuple or list containing 3 elements: x, y & z
       SGData, from SpcGroup
       ifAll=True causes the return to provide the unique set of 
                  equivalent positions
            =False causes the input position to be repeated. This is the default,
                   but why someone would want this, I am not sure.
    Returns a list of two element tuples: 
       The first element is the coordinate as a three-element array and 
       the second describes the symmetry used to generate the site, of form [-][C]SS
          C indicates a centering operation was used (omitted if the 1st, [0,0,0])
          SS is the symmetry operator number (1-24)
          - indicates the center of symmetry was used (omitted otherwise)      
    '''
    XYZEquiv = []
    Idup = []
    X = np.array(XYZ)
    MoveToUnitCell(X)
    XYZEquiv.append(np.array(X))
    Idup.append(1)
    for ic,cen in enumerate(SGData['SGCen']):
        C = np.array(cen)
        for invers in range(int(SGData['SGInv']+1)):
            for io,ops in enumerate(SGData['SGOps']):
                idup = ((io+1)+100*ic)*(1-2*invers)
                T = np.array(ops[1])
                M =  np.array(ops[0])
                newX = np.sum(M*X,axis=1)+T
                if invers:
                    newX = -newX
                newX += C
                MoveToUnitCell(newX)
                New = True
                if ifAll:
                    if np.allclose(newX,X,atol=0.0002):
                        New = False
                        idup = 0                    
                    XYZEquiv.append(newX)
                else:
                    for oldX in XYZEquiv[:-1]:
                        if np.allclose(newX,oldX,atol=0.0002):
                            New = False
                            idup = 0
                    if New or ifAll:
                        XYZEquiv.append(newX)
                if ifAll and len(XYZEquiv) == 2:
                    Idup.append(1)
                else:
                    Idup.append(idup)
    if ifAll:
        return zip(XYZEquiv[1:],Idup[1:])                  #eliminate duplicate initial entry 
    else:
        return zip(XYZEquiv,Idup)
                                   
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
    CSuinel = [[],                                             # 0th empty - indices are Fortran style
        [[1,1,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],    #  A  A  A  0  0  0
        [[1,1,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],    #  A  A  C  0  0  0
        [[1,2,1,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],    #  A  B  A  0  0  0
        [[1,2,2,0,0,0],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],    #  A  B  B  0  0  0
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],    #  A  A  A  D  D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0,-1.0]],    #  A  A  A  D -D -D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0]],    #  A  A  A  D -D  D
        [[1,1,1,2,2,2],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0]],    #  A  A  A  D  D -D
        [[1,1,2,1,0,0],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0]],    #  A  A  C A/2 0  0
        [[1,2,3,0,0,],[ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]],    #  A  B  C  0  0  0
        [[1,1,2,3,0,],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]],    #  A  A  C  D  0  0
        [[1,2,1,0,3,],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]],    #  A  B  A  0  E  0
        [[1,2,2,0,0,],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]],    #  A  B  B  0  0  F
        [[1,2,3,2,0,],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0]],    #  A  B  C B/2 0  0
        [[1,2,3,1,0,],[ 1.0, 1.0, 1.0, 0.5, 0.0, 0.0]],    #  A  B  C A/2 0  0
        [[1,2,3,4,0,],[ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]],    #  A  B  C  D  0  0
        [[1,2,3,0,4,],[ 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]],    #  A  B  C  0  E  0
        [[1,2,3,0,0,],[ 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]],    #  A  B  C  0  0  F
        [[1,1,2,3,4,],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0]],    #  A  A  C  D  E -E
        [[1,1,2,3,4,],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],    #  A  A  C  D  E  E
        [[1,2,1,3,4,],[ 1.0, 1.0, 1.0, 1.0, 1.0,-1.0]],    #  A  B  A  D  E -D
        [[1,2,1,3,4,],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],    #  A  B  A  D  E  D
        [[1,2,2,3,3,],[ 1.0, 1.0, 1.0, 1.0,-1.0, 1.0]],    #  A  B  B  D -D  F
        [[1,2,2,3,3,],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],    #  A  B  B  D  D  F
        [[1,2,3,2,4,],[ 1.0, 1.0, 1.0, 0.5, 0.5, 1.0]],    #  A  B  C B/2 F/2 F
        [[1,2,3,1,0,],[ 1.0, 1.0, 1.0, 0.5, 0.0, 1.0]],    #  A  B  C A/2  0  F
        [[1,2,3,2,4,],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.0]],    #  A  B  C B/2  E  0
        [[1,2,3,1,4,],[ 1.0, 1.0, 1.0, 0.5, 1.0, 0.5]],    #  A  B  C A/2  E E/2
        [[1,2,3,4,5,],[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],    #  A  B  C  D  E   F
        ]
    indx = GetNXUPQsym(siteSym)
    return CSuinel[indx[1]]
        
def SytSym(XYZ,SGData):
    '''
    Generates the number of equivalent positions and a site symmetry code for a specified coordinate and space group
    input:  
       XYZ: an array, tuple or list containing 3 elements: x, y & z
       SGData: from SpcGroup
    Returns a two element tuple:
       The 1st element is a code for the site symmetry (see GetOprPtrName)
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
    Jdup = 1
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
    
# self-test materials follow. Requires files in directory testinp
def test0():
    '''test #0: exercise MoveToUnitCell'''
    msg = "MoveToUnitCell failed"
    v = [0,1,2,-1,-2]; MoveToUnitCell(v); assert v==[0,0,0,0,0], msg
    v = np.array([-.1]); MoveToUnitCell(v); assert abs(v-0.9) < 1e-6, msg
    v = np.array([.1]); MoveToUnitCell(v); assert abs(v-0.1) < 1e-6, msg

def test1():
    ''' test #1: SpcGroup and SGPrint against previous results'''
    testdir = ospath.join(mypath,'testinp')
    if ospath.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    import spctestinp
    def CompareSpcGroup(spc, referr, refdict, reflist): 
        'Compare output from GSASIIspc.SpcGroup with results from a previous run'
        # if an error is reported, the dictionary can be ignored
        msg = "failed on space group %s" % spc
        result = SpcGroup(spc)
        if result[0] == referr and referr > 0: return True
        keys = result[1].keys()
        #print result[1]['SpGrp']
        assert len(keys) == len(refdict.keys()), msg
        for key in keys:
        #print key, type(refdict[key])
            if key == 'SGOps' or  key == 'SGCen':
                assert len(refdict[key]) == len(result[1][key]), msg
                for i in range(len(refdict[key])):
                    assert np.allclose(result[1][key][i][0],refdict[key][i][0]), msg
                    assert np.allclose(result[1][key][i][1],refdict[key][i][1]), msg
            else:
                assert result[1][key] == refdict[key], msg
        assert reflist == SGPrint(result[1]), 'SGPrint ' +msg
    for spc in spctestinp.SGdat:
        CompareSpcGroup(spc, 0, spctestinp.SGdat[spc], spctestinp.SGlist[spc] )

def test2():
    ''' test #2: SpcGroup against cctbx (sgtbx) computations'''
    testdir = ospath.join(mypath,'testinp')
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
                    MoveToUnitCell(noff)
                    mult = tuple((op*inv).ravel().tolist())
                    if debug: print "\n%s: %s + %s" % (spcname,mult,noff)
                    for refop in cctbx:
                        if debug: print refop
                        # check the transform
                        if refop[:9] != mult: continue
                        if debug: print "mult match"
                        # check the translation
                        reftrans = list(refop[-3:])
                        MoveToUnitCell(reftrans)
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
