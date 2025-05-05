# -*- coding: utf-8 -*-
'''
:mod:`GSASIIlattice` Classes & routines follow
'''
from __future__ import division, print_function
import math
import copy
import sys
import random as ran
import numpy as np
import numpy.linalg as nl
import scipy.special as spsp
from . import GSASIImath as G2mth
from . import GSASIIspc as G2spc
from . import GSASIIElem as G2elem
from . import GSASIIpath as GSASIIpath
try:
    if GSASIIpath.binaryPath:
        import pytexture as ptx
    else:
        from . import pytexture as ptx
        #import GSASII.pytexture as ptx
except ImportError: # ignore; will report this as an error in GSASIIplot import
    pass

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
tand = lambda x: np.tan(x*np.pi/180.)
atand = lambda x: 180.*np.arctan(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
cosd = lambda x: np.cos(x*np.pi/180.)
acosd = lambda x: 180.*np.arccos(x)/np.pi
rdsq2d = lambda x,p: round(1.0/math.sqrt(x),p)
nprdsq2d = lambda x,p: np.round(1.0/np.sqrt(x),p)
try:  # fails on doc build
    rpd = np.pi/180.
    RSQ2PI = 1./np.sqrt(2.*np.pi)
    SQ2 = np.sqrt(2.)
    RSQPI = 1./np.sqrt(np.pi)
    R2pisq = 1./(2.*np.pi**2)
    Forpi = 4.0*np.pi
except TypeError:
    pass
nxs = np.newaxis

def sec2HMS(sec):
    """Convert time in sec to H:M:S string

    :param sec: time in seconds
    :return: H:M:S string (to nearest 100th second)

    """
    H = int(sec//3600)
    M = int(sec//60-H*60)
    S = sec-3600*H-60*M
    return '%d:%2d:%.2f'%(H,M,S)

def rotdMat(angle,axis=0):
    """Prepare rotation matrix for angle in degrees about axis(=0,1,2)

    :param angle: angle in degrees
    :param axis:  axis (0,1,2 = x,y,z) about which for the rotation
    :return: rotation matrix - 3x3 numpy array

    """
    if axis == 2:
        return np.array([[cosd(angle),-sind(angle),0],[sind(angle),cosd(angle),0],[0,0,1]])
    elif axis == 1:
        return np.array([[cosd(angle),0,-sind(angle)],[0,1,0],[sind(angle),0,cosd(angle)]])
    else:
        return np.array([[1,0,0],[0,cosd(angle),-sind(angle)],[0,sind(angle),cosd(angle)]])

def rotdMat4(angle,axis=0):
    """Prepare rotation matrix for angle in degrees about axis(=0,1,2) with scaling for OpenGL

    :param angle: angle in degrees
    :param axis:  axis (0,1,2 = x,y,z) about which for the rotation
    :return: rotation matrix - 4x4 numpy array (last row/column for openGL scaling)

    """
    Mat = rotdMat(angle,axis)
    return np.concatenate((np.concatenate((Mat,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)

def fillgmat(cell):
    """Compute lattice metric tensor from unit cell constants

    :param cell: tuple with a,b,c,alpha, beta, gamma (degrees)
    :return: 3x3 numpy array

    """
    a,b,c,alp,bet,gam = cell
    g = np.array([
        [a*a,  a*b*cosd(gam),  a*c*cosd(bet)],
        [a*b*cosd(gam),  b*b,  b*c*cosd(alp)],
        [a*c*cosd(bet) ,b*c*cosd(alp),   c*c]])
    return g

def cell2Gmat(cell):
    """Compute real and reciprocal lattice metric tensor from unit cell constants

    :param cell: tuple with a,b,c,alpha, beta, gamma (degrees)
    :return: reciprocal (G) & real (g) metric tensors (list of two numpy 3x3 arrays)

    """
    g = fillgmat(cell)
    G = nl.inv(g)
    return G,g

def A2Gmat(A,inverse=True):
    """Fill real & reciprocal metric tensor (G) from A.

    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :param bool inverse: if True return both G and g; else just G
    :return: reciprocal (G) & real (g) metric tensors (list of two numpy 3x3 arrays)

    """
    G = np.array([
        [A[0],  A[3]/2.,  A[4]/2.],
        [A[3]/2.,A[1],    A[5]/2.],
        [A[4]/2.,A[5]/2.,    A[2]]])
    if inverse:
        g = nl.inv(G)
        return G,g
    else:
        return G

def Gmat2A(G):
    """Extract A from reciprocal metric tensor (G)

    :param G: reciprocal metric tensor (3x3 numpy array)
    :return: A = [G11,G22,G33,2*G12,2*G13,2*G23]

    """
    return [G[0][0],G[1][1],G[2][2],2.*G[0][1],2.*G[0][2],2.*G[1][2]]

def cell2A(cell):
    """Obtain A = [G11,G22,G33,2*G12,2*G13,2*G23] from lattice parameters

    :param cell: [a,b,c,alpha,beta,gamma] (degrees)
    :return: G reciprocal metric tensor as 3x3 numpy array

    """
    G,g = cell2Gmat(cell)
    return Gmat2A(G)

def A2cell(A):
    """Compute unit cell constants from A

    :param A: [G11,G22,G33,2*G12,2*G13,2*G23] G - reciprocal metric tensor
    :return: a,b,c,alpha, beta, gamma (degrees) - lattice parameters

    """
    G,g = A2Gmat(A)
    return Gmat2cell(g)

def Gmat2cell(g):
    """Compute real/reciprocal lattice parameters from real/reciprocal metric tensor (g/G)
    The math works the same either way.

    :param g (or G): real (or reciprocal) metric tensor 3x3 array
    :return: a,b,c,alpha, beta, gamma (degrees) (or a*,b*,c*,alpha*,beta*,gamma* degrees)

    """
    oldset = np.seterr('raise')
    a = np.sqrt(max(0,g[0][0]))
    b = np.sqrt(max(0,g[1][1]))
    c = np.sqrt(max(0,g[2][2]))
    alp = acosd(g[2][1]/(b*c))
    bet = acosd(g[2][0]/(a*c))
    gam = acosd(g[0][1]/(a*b))
    np.seterr(**oldset)
    return a,b,c,alp,bet,gam

def invcell2Gmat(invcell):
    """Compute real and reciprocal lattice metric tensor from reciprocal
       unit cell constants

    :param invcell: [a*,b*,c*,alpha*, beta*, gamma*] (degrees)
    :return: reciprocal (G) & real (g) metric tensors (list of two 3x3 arrays)

    """
    G = fillgmat(invcell)
    g = nl.inv(G)
    return G,g

def cellDijFill(pfx,phfx,SGData,parmDict):
    '''Returns the filled-out reciprocal cell (A) terms
    from the parameter dictionaries corrected for Dij.

    :param str pfx: parameter prefix ("n::", where n is a phase index)
    :param str phfx: parameter prefix ("n:h:", where n is a phase index 
       and h is a histogram index)
    :param dict SGdata: a symmetry object
    :param dict parmDict: a dictionary of parameters

    :returns: A, a list of six terms
    '''
    if phfx+'D11' not in parmDict:
        return None
    if SGData['SGLaue'] in ['-1',]:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A1']+parmDict[phfx+'D22'],
             parmDict[pfx+'A2']+parmDict[phfx+'D33'],
             parmDict[pfx+'A3']+parmDict[phfx+'D12'],parmDict[pfx+'A4']+parmDict[phfx+'D13'],
             parmDict[pfx+'A5']+parmDict[phfx+'D23']]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A1']+parmDict[phfx+'D22'],
                 parmDict[pfx+'A2']+parmDict[phfx+'D33'],0,0,parmDict[pfx+'A5']+parmDict[phfx+'D23']]
        elif SGData['SGUniq'] == 'b':
            A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A1']+parmDict[phfx+'D22'],
                 parmDict[pfx+'A2']+parmDict[phfx+'D33'],0,parmDict[pfx+'A4']+parmDict[phfx+'D13'],0]
        else:
            A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A1']+parmDict[phfx+'D22'],
                 parmDict[pfx+'A2']+parmDict[phfx+'D33'],parmDict[pfx+'A3']+parmDict[phfx+'D12'],0,0]
    elif SGData['SGLaue'] in ['mmm',]:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A1']+parmDict[phfx+'D22'],
             parmDict[pfx+'A2']+parmDict[phfx+'D33'],0,0,0]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A0']+parmDict[phfx+'D11'],
             parmDict[pfx+'A2']+parmDict[phfx+'D33'],0,0,0]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A0']+parmDict[phfx+'D11'],
             parmDict[pfx+'A2']+parmDict[phfx+'D33'],parmDict[pfx+'A0']+parmDict[phfx+'D11'],0,0]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A0']+parmDict[phfx+'D11'],
            parmDict[pfx+'A0']+parmDict[phfx+'D11'],
            parmDict[pfx+'A3']+parmDict[phfx+'D23'],parmDict[pfx+'A3']+parmDict[phfx+'D23'],
            parmDict[pfx+'A3']+parmDict[phfx+'D23']]
    elif SGData['SGLaue'] in ['m3m','m3']:
        A = [parmDict[pfx+'A0']+parmDict[phfx+'D11'],parmDict[pfx+'A0']+parmDict[phfx+'D11'],
             parmDict[pfx+'A0']+parmDict[phfx+'D11'],0,0,0]
    return A

def getCellSU(pId,hId,SGData,parmDict,covData):
    '''Compute the unit cell parameters and standard uncertainties
    where lattice parameters and Hstrain (Dij) may be refined. 

    :param pId: phase index
    :param hId: histogram index
    :param SGdata: space group information for current phase
    :param dict parmDict: parameter dict, must have all non-zero Dij and Ai terms
    :param dict covData: covariance tree item 
    :returns: cellList,cellSig where each term is a list of 7 items, 
       a, b,... Vol and sig(a), sig(b),... sig(Vol), respectively
    '''

    Dnames = [f'{pId}:{hId}:D{i}' for i in ['11','22','33','12','13','23']]
    Anames = [f'{pId}::A{i}' for i in range(6)]
    pfx = f'{pId}::'
    phfx = f'{pId}:{hId}:'
    A = cellDijFill(pfx,phfx,SGData,parmDict)
    cell = list(A2cell(A)) + [calc_V(A)]
    
    rVsq = calc_rVsq(A)
    G,g = A2Gmat(A)       #get recip. & real metric tensors
    varyList = covData['varyList']
    covMatrix = covData['covMatrix']
    if len(covMatrix):
        vcov = G2mth.getVCov(Anames+Dnames,varyList,covMatrix)
        for i in [0,6]:
            for j in [0,6]:
                if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm']:
                    vcov[1+i,1+j] = vcov[3+i,3+j] = vcov[i,1+j] = vcov[1+i,j] = vcov[i,j]
                    vcov[1+i,3+j] = vcov[3+i,1+j] = vcov[i,3+j] = vcov[3+i,j] = vcov[i,j]
                    vcov[1+i,2+j] = vcov[2+i,1+j] = vcov[2+i,3+j] = vcov[3+i,2+j] = vcov[i,2+j]
                elif SGData['SGLaue'] in ['m3','m3m']:
                    vcov[i:3+i,j:3+j] = vcov[i,j]
                elif SGData['SGLaue'] in ['4/m', '4/mmm']:
                    vcov[i:2+i,j:2+j] = vcov[i,j]
                    vcov[1+i,2+j] = vcov[2+i,1+j] = vcov[i,2+j]
                elif SGData['SGLaue'] in ['3R','3mR']:
                    vcov[i:3+j,i:3+j] = vcov[i,j]
                    #        vcov[4,4] = vcov[5,5] = vcov[3,3]
                    vcov[3+i:6+i,3+j:6+j] = vcov[3,3+j]
                    vcov[i:3+i,3+j:6+j] = vcov[i,3+j]
                    vcov[3+i:6+i,j:3+j] = vcov[3+i,j]
    else:
        vcov = np.eye(12)
    delt = 1.e-9
    drVdA = np.zeros(12)
    for i in range(12):
        A[i%6] += delt
        drVdA[i] = calc_rVsq(A)
        A[i%6] -= 2*delt
        drVdA[i] -= calc_rVsq(A)
        A[i%6] += delt
    drVdA /= 2.*delt
    srcvlsq = np.inner(drVdA,np.inner(drVdA,vcov))
    Vol = 1/np.sqrt(rVsq)
    sigVol = Vol**3*np.sqrt(srcvlsq)/2.         #ok - checks with GSAS
    
    dcdA = np.zeros((12,12))
    for i in range(12):
        pdcdA =np.zeros(12)
        A[i%6] += delt
        pdcdA += A2cell(A) + A2cell(A)
        A[i%6] -= 2*delt
        pdcdA -= A2cell(A) + A2cell(A)
        A[i%6] += delt
        dcdA[i] = pdcdA/(2.*delt)
    sigMat = np.inner(dcdA,np.inner(dcdA,vcov))
    var = np.diag(sigMat)
    CS = np.where(var>0.,np.sqrt(var),0.)
    if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm','m3','m3m','4/m','4/mmm']:
        CS[3:6] = 0.0
    # show s.u. values only for the unique values
    if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm','4/m', '4/mmm']:
        CS[1] = -CS[1]
    elif SGData['SGLaue'] in ['m3','m3m']:
        CS[1] = -CS[1]
        CS[2] = -CS[2]
    elif SGData['SGLaue'] in ['3R','3mR']:
        CS[1] = -CS[1]
        CS[2] = -CS[2]
        CS[4] = -CS[4]
        CS[3] = -CS[3]
    return cell,[CS[0],CS[1],CS[2],CS[5],CS[4],CS[3],sigVol]

def showCellSU(cellList,cellSig,SGData,cellNames=None):
    '''Produce the cell parameters from :func:`getCellSU` as a nicely 
    formatted string

    :param list cellList: a list of 7 items, a, b,... Vol, from getCellSU
    :param list cellSig: a list of 7 items, sig(a), sig(b),... sig(Vol), 
      from getCellSU
    :param SGdata: space group information for current phase
    :param list cellNames: if specified, should be the labels to be 
      used for a, b,... volume. Defaults to cellAlbl & 'vol', but for 
      on-screen use, cellUlbl might be better than cellAlbl
    '''
    defsigL = 3*[-0.00001] + 3*[-0.001] + [-0.01] # significance to use when no sigma
    if cellNames is None:
        cellNames = list(cellAlbl) + ['vol']
    laue = SGData['SGLaue']
    if laue == '2/m':
        laue += SGData['SGUniq']
    for symlist,celllist in UniqueCellByLaue:
        if laue in symlist:
            uniqCellIndx = celllist
            break
    else: # should not happen
        uniqCellIndx = list(range(6))
    prevsig = 0
    s = ''
    for i,(lbl,defsig,val,sig) in enumerate(zip(cellNames,defsigL,cellList,cellSig)):
        if i != 6 and i not in uniqCellIndx: continue
        if sig:
            txt = G2mth.ValEsd(val,sig)
            prevsig = -sig # use this as the significance for next value
        else:
            txt = G2mth.ValEsd(val,min(defsig,prevsig),True)
        s += f' {lbl}={txt}'
    return s

def getCellEsd(pfx,SGData,A,covData,unique=False):
    '''Compute the standard uncertainty on cell parameters
    
    :param str pfx: prefix of form "p::"
    :param SGdata: space group information
    :param list A: Reciprocal cell Ai terms
    :param dict covData: covariance tree item 
    :param bool unique: when True, only directly refined parameters
      (a in cubic, a & alpha in rhombohedral cells) are assigned 
      positive s.u. values. Used for CIF generation.
    '''
    rVsq = calc_rVsq(A)
    G,g = A2Gmat(A)       #get recip. & real metric tensors
    RMnames = [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    varyList = covData['varyList']
    covMatrix = covData['covMatrix']
    if len(covMatrix):
        vcov = G2mth.getVCov(RMnames,varyList,covMatrix)
        if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm']:
            vcov[1,1] = vcov[3,3] = vcov[0,1] = vcov[1,0] = vcov[0,0]
            vcov[1,3] = vcov[3,1] = vcov[0,3] = vcov[3,0] = vcov[0,0]
            vcov[1,2] = vcov[2,1] = vcov[2,3] = vcov[3,2] = vcov[0,2]
        elif SGData['SGLaue'] in ['m3','m3m']:
            vcov[0:3,0:3] = vcov[0,0]
        elif SGData['SGLaue'] in ['4/m', '4/mmm']:
            vcov[0:2,0:2] = vcov[0,0]
            vcov[1,2] = vcov[2,1] = vcov[0,2]
        elif SGData['SGLaue'] in ['3R','3mR']:
            vcov[0:3,0:3] = vcov[0,0]
    #        vcov[4,4] = vcov[5,5] = vcov[3,3]
            vcov[3:6,3:6] = vcov[3,3]
            vcov[0:3,3:6] = vcov[0,3]
            vcov[3:6,0:3] = vcov[3,0]
    else:
        vcov = np.eye(6)
    delt = 1.e-9
    drVdA = np.zeros(6)
    for i in range(6):
        A[i] += delt
        drVdA[i] = calc_rVsq(A)
        A[i] -= 2*delt
        drVdA[i] -= calc_rVsq(A)
        A[i] += delt
    drVdA /= 2.*delt    
    srcvlsq = np.inner(drVdA,np.inner(drVdA,vcov))
    Vol = 1/np.sqrt(rVsq)
    sigVol = Vol**3*np.sqrt(srcvlsq)/2.         #ok - checks with GSAS
    
    dcdA = np.zeros((6,6))
    for i in range(6):
        pdcdA =np.zeros(6)
        A[i] += delt
        pdcdA += A2cell(A)
        A[i] -= 2*delt
        pdcdA -= A2cell(A)
        A[i] += delt
        dcdA[i] = pdcdA/(2.*delt)
    
    sigMat = np.inner(dcdA,np.inner(dcdA,vcov))
    var = np.diag(sigMat)
    CS = np.where(var>0.,np.sqrt(var),0.)
    if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm','m3','m3m','4/m','4/mmm']:
        CS[3:6] = 0.0
    # show s.u. values only for the unique values
    if not unique:
        pass
    elif SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm','4/m', '4/mmm']:
        CS[1] = -CS[1]
    elif SGData['SGLaue'] in ['m3','m3m']:
        CS[1] = -CS[1]
        CS[2] = -CS[2]
    elif SGData['SGLaue'] in ['3R','3mR']:
        CS[1] = -CS[1]
        CS[2] = -CS[2]
        CS[4] = -CS[4]
        CS[3] = -CS[3]
    return [CS[0],CS[1],CS[2],CS[5],CS[4],CS[3],sigVol]

def CellDijCorr(Cell,SGData,Data,hist):
    '''Returns the cell corrected for Dij values.

    :param list Cell: lattice parameters
    :param dict SGdata: a symmetry object
    :param dict Data: phase data structure; contains set of Dij values
    :param str hist: histogram name

    :returns: cell corrected for Dij values
    '''
    A = cell2A(Cell)
    Dij = Data[hist]['HStrain'][0]
    newA = AplusDij(A,Dij,SGData)
    return A2cell(newA)

def AplusDij(A,Dij,SGData):
    ''' returns the A corrected by Dij

    :param list A: reciprocal metric terms A0-A5
    :param array Dij: unique Dij values as stored in Hstrain
    :param dict SGdata: a symmetry object

    :returns list newA: A corrected by Dij
    '''
    if SGData['SGLaue'] in ['-1',]:
        newA = [A[0]+Dij[0],A[1]+Dij[1],A[2]+Dij[2],A[3]+Dij[3],A[4]+Dij[4],A[5]+Dij[5]]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            newA = [A[0]+Dij[0],A[1]+Dij[1],A[2]+Dij[2],0,0,A[5]+Dij[3]]
        elif SGData['SGUniq'] == 'b':
            newA = [A[0]+Dij[0],A[1]+Dij[1],A[2]+Dij[2],0,A[4]+Dij[3],0]
        else:
            newA = [A[0]+Dij[0],A[1]+Dij[1],A[2]+Dij[2],A[3]+Dij[3],0,0]
    elif SGData['SGLaue'] in ['mmm',]:
        newA = [A[0]+Dij[0],A[1]+Dij[1],A[2]+Dij[2],0,0,0]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        newA = [A[0]+Dij[0],A[0]+Dij[0],A[2]+Dij[1],0,0,0]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        newA = [A[0]+Dij[0],A[0]+Dij[0],A[2]+Dij[1],A[0]+Dij[0],0,0]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        newA = [A[0]+Dij[0],A[0]+Dij[0],A[0]+Dij[0],A[3]+Dij[1],A[3]+Dij[1],A[3]+Dij[1]]
    elif SGData['SGLaue'] in ['m3m','m3']:
        newA = [A[0]+Dij[0],A[0]+Dij[0],A[0]+Dij[0],0,0,0]

    return newA

def prodMGMT(G,Mat):
    '''Transform metric tensor by matrix

    :param G: array metric tensor
    :param Mat: array transformation matrix
    :return: array new metric tensor

    '''
    return np.inner(np.inner(Mat,G),Mat)        #right
#    return np.inner(Mat,np.inner(Mat,G))       #right
#    return np.inner(np.inner(G,Mat).T,Mat)      #right
#    return np.inner(Mat,np.inner(G,Mat).T)     #right

def TransformCell(cell,Trans):
    '''Transform lattice parameters by matrix

    :param cell: list a,b,c,alpha,beta,gamma,(volume)
    :param Trans: array transformation matrix
    :return: array transformed a,b,c,alpha,beta,gamma,volume

    '''
    newCell = np.zeros(7)
    g = cell2Gmat(cell)[1]
    newg = prodMGMT(g,Trans)
    newCell[:6] = Gmat2cell(newg)
    newCell[6] = calc_V(cell2A(newCell[:6]))
    return newCell

# code to generate lattice constraint relationships between two phases
# (chemical & magnetic) related by a transformation matrix.

def symInner(M1,M2):
    '''Compute inner product of two square matrices with symbolic processing
    Use dot product because sympy does not define an inner product primitive

    This requires that M1 & M2 be two sympy objects, as created in
    GenerateCellConstraints().

    Note that this is only used to do the symbolic math needed to generate
    cell relationships. It is not used normally in GSAS-II.
    '''
    import sympy as sym
    prodOuter = []
    for i in range(3):
        prod = []
        for j in range(3):
            prod.append(M1[i,:].dot(M2[j,:]))
        prodOuter.append(prod)
    return sym.Matrix(prodOuter)

def GenerateCellConstraints():
    '''Generate unit cell constraints for transforming one set of A tensor
    values to another using symbolic math (requires the sympy package)

    Note that this is only used to do the symbolic math needed to generate
    cell relationships. It is not used normally in GSAS-II.
    '''
    import sympy as sym
    # define A tensor for starting cell
    A0, A1, A2, A3, A4, A5 = sym.symbols('A0, A1, A2, A3, A4, A5')
    G = sym.Matrix([ [A0,    A3/2.,  A4/2.],
                     [A3/2.,  A1,    A5/2.],
                     [A4/2.,  A5/2.,   A2 ]] )
    # define transformation matrix
    T00, T10, T20, T01, T11, T21, T02, T12, T22 = sym.symbols(
        'T00, T10, T20, T01, T11, T21, T02, T12, T22')
    Tr = sym.Matrix([ [T00, T10, T20], [T01, T11, T21], [T02, T12, T22],])
    # apply transform
    newG = symInner(symInner(Tr,G),Tr).expand()
    # define A tensor for converted cell
    return [newG[0,0],newG[1,1],newG[2,2],2.*newG[0,1],2.*newG[0,2],2.*newG[1,2]]

def subVals(expr,A,T):
    '''Evaluate the symbolic expressions by substituting for A0-A5 & Tij

    This can be used on the cell relationships created in
    :func:`GenerateCellConstraints` like this::

       Trans = np.array([ [2/3, 4/3, 1/3], [-1, 0, 0], [-1/3, -2/3, 1/3] ])
       T = np.linalg.inv(Trans).T
       print([subVals(i,Aold,T) for i in GenerateCellConstraints()])

    :param list expr: a list of sympy expressions.
    :param list A: This is the A* tensor as defined above.
    :param np.array T: a 3x3 transformation matrix where,
       Trans = np.array([ [2/3, 4/3, 1/3], [-1, 0, 0], [-1/3, -2/3, 1/3] ])
       (for a' = 2/3a + 4/3b + 1/3c; b' = -a; c' = -1/3, -2/3, 1/3)
       then T = np.linalg.inv(Trans).T

    Note that this is only used to do the symbolic math needed to generate
    cell relationships. It is not used normally in GSAS-II.
    '''
    import sympy as sym
    A0, A1, A2, A3, A4, A5 = sym.symbols('A0, A1, A2, A3, A4, A5')
    # transformation matrix
    T00, T10, T20, T01, T11, T21, T02, T12, T22 = sym.symbols(
        'T00, T10, T20, T01, T11, T21, T02, T12, T22')
    vals = dict(zip([T00, T10, T20, T01, T11, T21, T02, T12, T22],T.ravel()))
    vals.update(dict(zip([A0, A1, A2, A3, A4, A5],A)))
    return float(expr.subs(vals))

# some sample test code using the routines above follows::
# Trans = np.array([ [2/3, 4/3, 1/3], [-1, 0, 0], [-1/3, -2/3, 1/3] ])
# Mat = np.linalg.inv(Trans).T
# Aold = [0.05259986634758891, 0.05259986634758891, 0.005290771904550856, 0.052599866347588925, 0, 0]
# Anew = [0.018440738491448085, 0.03944989976069168, 0.034313054205100654, 0, -0.00513684555559103, 0]
# cellConstr = GenerateCellConstraints()
# calcA = [subVals(i,Aold,Mat) for i in cellConstr]
# print('original   xform A',Anew)
# print('calculated xfrom A',calcA)
# print('input')
# print('  old cell',A2cell(Aold))
# print('  new cell',A2cell(Anew))
# print('derived results')
# print('  from eq.',A2cell(calcA))
# print('  diffs   ',np.array(A2cell(calcA)) - A2cell(Anew))

def fmtCellConstraints(cellConstr):
    '''Format the cell relationships created in :func:`GenerateCellConstraints`
    in a format that can be used to generate constraints.

    Use::

      cXforms = G2lat.fmtCellConstraints(G2lat.GenerateCellConstraints())

    Note that this is only used to do the symbolic math needed to generate
    cell relationships. It is not used normally in GSAS-II.
    '''
    import re
    import sympy as sym
    A3, A4, A5 = sym.symbols('A3, A4, A5')
    consDict = {}
    for num,cons in enumerate(cellConstr):
        cons = str(cons.factor(A3,A4,A5,deep=True).simplify())
        cons = re.sub('T([0-2]?)([0-2]?)',r'T[\2,\1]',cons) # Tij to T[j,i]
        l = []
        for i in str(cons).split('+'):
            if ')' in i:
                l[-1] += ' + ' + i.strip()
            else:
                l.append(i.strip())
        print("\nA'{} = ".format(num),str(cons))
        consDict[num] = l
    return consDict

cellXformRelations = {0: ['1.0*A0*T[0,0]**2',
                              '1.0*A1*T[0,1]**2',
                              '1.0*A2*T[0,2]**2',
                              '1.0*A3*T[0,0]*T[0,1]',
                              '1.0*A4*T[0,0]*T[0,2]',
                              '1.0*A5*T[0,1]*T[0,2]'],
                    1: ['1.0*A0*T[1,0]**2',
                            '1.0*A1*T[1,1]**2',
                            '1.0*A2*T[1,2]**2',
                            '1.0*A3*T[1,0]*T[1,1]',
                            '1.0*A4*T[1,0]*T[1,2]',
                            '1.0*A5*T[1,1]*T[1,2]'],
                    2: ['1.0*A0*T[2,0]**2',
                            '1.0*A1*T[2,1]**2',
                            '1.0*A2*T[2,2]**2',
                            '1.0*A3*T[2,0]*T[2,1]',
                            '1.0*A4*T[2,0]*T[2,2]',
                            '1.0*A5*T[2,1]*T[2,2]'],
                    3: ['2.0*A0*T[0,0]*T[1,0]',
                            '2.0*A1*T[0,1]*T[1,1]',
                            '2.0*A2*T[0,2]*T[1,2]',
                            '1.0*A3*(T[0,0]*T[1,1] + T[1,0]*T[0,1])',
                            '1.0*A4*(T[0,0]*T[1,2] + T[1,0]*T[0,2])',
                            '1.0*A5*(T[0,1]*T[1,2] + T[1,1]*T[0,2])'],
                    4: ['2.0*A0*T[0,0]*T[2,0]',
                            '2.0*A1*T[0,1]*T[2,1]',
                            '2.0*A2*T[0,2]*T[2,2]',
                            '1.0*A3*(T[0,0]*T[2,1] + T[2,0]*T[0,1])',
                            '1.0*A4*(T[0,0]*T[2,2] + T[2,0]*T[0,2])',
                            '1.0*A5*(T[0,1]*T[2,2] + T[2,1]*T[0,2])'],
                    5: ['2.0*A0*T[1,0]*T[2,0]',
                            '2.0*A1*T[1,1]*T[2,1]',
                            '2.0*A2*T[1,2]*T[2,2]',
                            '1.0*A3*(T[1,0]*T[2,1] + T[2,0]*T[1,1])',
                            '1.0*A4*(T[1,0]*T[2,2] + T[2,0]*T[1,2])',
                            '1.0*A5*(T[1,1]*T[2,2] + T[2,1]*T[1,2])']}

'''cellXformRelations provide the constraints on newA[i] values for a new
cell generated from oldA[i] values.
'''
# cellXformRelations values were generated using::
#   from GSASIIlattice import fmtCellConstraints,GenerateCellConstraints
#   cellXformRelations = fmtCellConstraints(GenerateCellConstraints())

def GenCellConstraints(Trans,origPhase,newPhase,origA,oSGLaue,nSGLaue,debug=False):
    '''Generate the constraints between two unit cells constants for a phase transformed
    by matrix Trans.

    :param np.array Trans: a 3x3 direct cell transformation matrix where,
       Trans = np.array([ [2/3, 4/3, 1/3], [-1, 0, 0], [-1/3, -2/3, 1/3] ])
       (for a' = 2/3a + 4/3b + 1/3c; b' = -a; c' = -1/3, -2/3, 1/3)
    :param int origPhase: phase id (pId) for original phase
    :param int newPhase: phase id for the transformed phase to be constrained from
      original phase
    :param list origA: reciprocal cell ("A*") tensor (used for debug only)
    :param dict oSGLaue: space group info for original phase
    :param dict nSGLaue: space group info for transformed phase
    :param bool debug: If true, the constraint input is used to compute and print A*
      and from that the direct cell for the transformed phase.
    '''
    from . import GSASIIobj as G2obj
    Anew = []
    constrList = []
    uniqueAnew = cellUnique(nSGLaue)
    zeroAorig = cellZeros(oSGLaue)
    for i in range(6):
        constr = [[-1.0,G2obj.G2VarObj('{}::A{}'.format(newPhase,i))]]
        mult = []
        for j,item in enumerate(cellXformRelations[i]):
            const, aTerm, tTerm = item.split('*',2)
            const = float(const) * eval(tTerm)
            mult.append(const)
            # skip over A terms that are required to be zero
            if zeroAorig[int(aTerm[1])]: continue   # only add non-zero terms
            # ignore terms where either the Transform contribution is zero [= abs() < 1e-8]
            # If the multiplier term is zero I don't think this accidental
            # but since it will not change there is no reason to include that
            # term in any case
            if abs(const) < 1e-8: continue
            constr.append([const,G2obj.G2VarObj('{}::{}'.format(origPhase,aTerm))])
        if i in uniqueAnew:
            constrList.append(constr + [0.0,None,'c'])
        if debug: Anew.append(np.dot(origA,mult))
    if debug:
        print('xformed A*  ',Anew)
        print('xformed cell',A2cell(Anew))
    return constrList

def cellUnique(SGData):
    '''Returns the indices for the unique A tensor terms
    based on the Laue class.
    Any terms that are determined from others or are zero are not included.

    :param dict SGdata: a symmetry object
    :returns: a list of 0 to 6 terms with indices of the unique A terms
    '''
    if SGData['SGLaue'] in ['-1',]:
        return [0,1,2,3,4,5]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [0,1,2,5]
        elif SGData['SGUniq'] == 'b':
            return [0,1,2,4]
        else:
            return [0,1,2,3]
    elif SGData['SGLaue'] in ['mmm',]:
        return [0,1,2]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        return [0,2]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        return [0,2]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        return [0,3]
    elif SGData['SGLaue'] in ['m3m','m3']:
        return [0,]

def cellZeros(SGData):
    '''Returns a list with the A terms required to be zero based on Laue symmetry

    :param dict SGdata: a symmetry object
    :returns: A list of six terms where the values are True if the
      A term must be zero, False otherwise.
    '''
    if SGData['SGLaue'] in ['-1',]:
        return 6*[False]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [False,False,False,True,True,False]
        elif SGData['SGUniq'] == 'b':
            return [False,False,False,True,False,True]
        else:
            return [False,False,False,False,True,True]
    elif SGData['SGLaue'] in ['mmm',]:
        return [False,False,False,True,True,True]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        return [False,False,False,True,True,True]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        return [False,False,False,False,True,True]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        return 6*[False]
    elif SGData['SGLaue'] in ['m3m','m3']:
        return [False,False,False,True,True,True]

def TransformXYZ(XYZ,Trans,Vec):
    return np.inner(XYZ,Trans)+Vec

def TransformU6(U6,Trans):
    Uij = np.inner(Trans,np.inner(U6toUij(U6),Trans).T)/nl.det(Trans)
    return UijtoU6(Uij)

def ExpandCell(Atoms,atCodes,cx,Trans):
    Unit = [int(max(abs(np.array(unit)))-1) for unit in Trans.T]
    nUnit = (Unit[0]+1)*(Unit[1]+1)*(Unit[2]+1)
    Ugrid = np.mgrid[0:Unit[0]+1,0:Unit[1]+1,0:Unit[2]+1]
    Ugrid = np.reshape(Ugrid,(3,nUnit)).T
    Codes = copy.deepcopy(atCodes)
    newAtoms = copy.deepcopy(Atoms)
    for unit in Ugrid[1:]:
        moreAtoms = copy.deepcopy(Atoms)
        for atom in moreAtoms:
            atom[cx:cx+3] += unit
        newAtoms += moreAtoms
        codes = copy.deepcopy(atCodes)
        moreCodes = [code+'+%d,%d,%d'%(unit[0],unit[1],unit[2]) for code in codes]
        Codes += moreCodes
    return newAtoms,Codes

def TransformPhase(oldPhase,newPhase,Trans,Uvec,Vvec,ifMag,Force=True):
    '''Transform atoms from oldPhase to newPhase
    M' is inv(M)
    does X' = M(X-U)+V transformation for coordinates and U' = MUM/det(M)
    for anisotropic thermal parameters

    :param oldPhase: dict G2 phase info for old phase
    :param newPhase: dict G2 phase info for new phase; with new cell & space group
            atoms are from oldPhase & will be transformed
    :param Trans: lattice transformation matrix M
    :param Uvec: array parent coordinates transformation vector U
    :param Vvec: array child coordinate transformation vector V
    :param ifMag: bool True if convert to magnetic phase;
        if True all nonmagnetic atoms will be removed

    :return: newPhase dict modified G2 phase info
    :return: atCodes list atom transformation codes

    '''
    cx,ct,cs,cia = oldPhase['General']['AtomPtrs']
    cm = 0
    if oldPhase['General']['Type'] == 'magnetic':
        cm = cx+4
    oAmat,oBmat = cell2AB(oldPhase['General']['Cell'][1:7])
    nAmat,nBmat = cell2AB(newPhase['General']['Cell'][1:7])
    SGData = newPhase['General']['SGData']
    invTrans = nl.inv(Trans)
    newAtoms,atCodes = FillUnitCell(oldPhase,Force)
    newAtoms,atCodes = ExpandCell(newAtoms,atCodes,cx,Trans)
    if ifMag:
        cia += 3
        cs += 3
        newPhase['General']['Type'] = 'magnetic'
        newPhase['General']['AtomPtrs'] = [cx,ct,cs,cia]
        magAtoms = []
        magatCodes = []
        Landeg = 2.0
        for iat,atom in enumerate(newAtoms):
            if len(G2elem.GetMFtable([atom[ct],],[Landeg,])):
                magAtoms.append(atom[:cx+4]+[0.,0.,0.]+atom[cx+4:])
                magatCodes.append(atCodes[iat])
        newAtoms = magAtoms
        atCodes = magatCodes
        newPhase['Draw Atoms'] = []
    for atom in newAtoms:
        xyz = TransformXYZ(atom[cx:cx+3]+Uvec,invTrans.T,Vvec)
        if Force:
            xyz = np.around(xyz,6)%1.
        atom[cx:cx+3] = xyz
        if atom[cia] == 'A':
            atom[cia+2:cia+8] = TransformU6(atom[cia+2:cia+8],Trans)
        atom[cs:cs+2] = G2spc.SytSym(atom[cx:cx+3],SGData)[:2]
        atom[cia+8] = ran.randint(0,sys.maxsize)
        if cm:
            mag = np.sqrt(np.sum(np.array(atom[cm:cm+3])**2))
            if mag:
                mom = np.inner(np.array(atom[cm:cm+3]),oBmat)
                mom = np.inner(mom,invTrans)
                mom = np.inner(mom,nAmat)
                mom /= np.sqrt(np.sum(mom**2))
                atom[cm:cm+3] = mom*mag
    newPhase['Atoms'] = newAtoms
    if SGData['SpGrp'] != 'P 1':
        newPhase['Atoms'],atCodes = GetUnique(newPhase,atCodes)
    newPhase['Drawing'] = []
    newPhase['ranId'] = ran.randint(0,sys.maxsize)
    return newPhase,atCodes

def FindNonstandard(controls,Phase):
    '''
    Find nonstandard setting of magnetic cell that aligns with parent nuclear cell

    :param controls: list unit cell indexing controls
    :param Phase: dict new magnetic phase data (NB:not G2 phase construction); modified here
    :return: None

    '''
    abc = np.eye(3)
    cba = np.rot90(np.eye(3))
    cba[1,1] *= -1      #makes c-ba
    Mats = {'abc':abc,'cab':np.roll(abc,2,1),'bca':np.roll(abc,1,1),
            'acb':np.roll(cba,1,1),'bac':np.roll(cba,2,1),'cba':cba}        #ok
    BNS = {'A':{'abc':'A','cab':'C','bca':'B','acb':'A','bac':'B','cba':'C'},
           'B':{'abc':'B','cab':'A','bca':'C','acb':'C','bac':'A','cba':'B'},
           'C':{'abc':'C','cab':'B','bca':'A','acb':'B','bac':'C','cba':'A'},
           'a':{'abc':'a','cab':'c','bca':'b','acb':'a','bac':'b','cba':'c'},   #Ok
           'b':{'abc':'b','cab':'a','bca':'c','acb':'c','bac':'a','cba':'b'},
           'c':{'abc':'c','cab':'b','bca':'a','acb':'b','bac':'c','cba':'a'},
           'S':{'abc':'S','cab':'S','bca':'S','acb':'S','bac':'S','cba':'S'},
           'I':{'abc':'I','cab':'I','bca':'I','acb':'I','bac':'I','cba':'I'},
           }
    Trans = Phase['Trans']
    Uvec = Phase['Uvec']
    SGData = Phase['SGData']
    MSG = SGData.get('MagSpGrp',SGData['SpGrp']).split(' ',1)
    MSG[0] += ' '
    bns = ''
    if '_' in MSG[0]:
        bns = MSG[0][2]
    spn = SGData.get('SGSpin',[])
    if 'ortho' in SGData['SGSys']:
        lattSym = G2spc.getlattSym(Trans)
        SpGrp = SGData['SpGrp']
        NTrans = np.inner(Mats[lattSym].T,Trans.T)        #ok
        if len(spn): spn[1:4] = np.inner(np.abs(nl.inv(Mats[lattSym])),spn[1:4])         #ok
        SGsym = G2spc.getlattSym(nl.inv(Mats[lattSym]))

        if lattSym != 'abc':
            NSG = G2spc.altSettingOrtho[SpGrp][SGsym].replace("'",'').split(' ')
            if ' '.join(NSG) in ['P 2 21 2',]:
                Uvec[1] += .25
            elif ' '.join(NSG) in ['P 21 2 2',]:
                Uvec[0] += .25
            elif ' '.join(NSG) in ['P 2 2 21',]:
                Uvec[2] += .25
            Bns = ''
            if bns:
                Bns = BNS[bns][lattSym]
                NSG[0] += '_'+Bns+' '
            elif len(spn):
                for ifld in [1,2,3]:
                    if spn[ifld] < 0:
                        NSG[ifld] += "'"
            Nresult = [''.join(NSG)+'  ',Bns]
            return Nresult,Uvec,NTrans
        else:
            return None
    elif 'mono' in SGData['SGSys']: # and not 'P_A' in Phase['Name']:  #skip the one that doesn't work
        newcell = TransformCell(controls[6:12],Trans)
        MatsA = np.array([[1.,0.,0.],[0.,1.,0.],[1.,0,1.]])
        MatsB = np.array([[1.,0.,0.],[0.,1.,0.],[-1.,0,1.]])
        if not 70. < newcell[4] < 110.:
            MSG[1] = MSG[1].replace('c','n')
            MSG[0] = MSG[0].replace('C_c','C_B').replace('P_A','P ')
            if '_' in MSG[0]:
                bns = MSG[0][2]
            if newcell[4] > 110.:
                if newcell[2] > newcell[0]:
                    Mats = MatsA
                else:
                    MSG[1] = MSG[1].replace('n','c')
                    MSG[0] = MSG[0].replace('C ','I ')
                    Mats = MatsA.T
            elif newcell[4] < 70.:
                if newcell[2] > newcell[0]:
                    Mats = MatsB
                else:
                    MSG[1] = MSG[1].replace('n','c')
                    MSG[0] = MSG[0].replace('C ','I ')
                    Mats = MatsB.T
            Nresult = [' '.join(MSG)+' ',bns]
            NTrans = np.inner(Mats,Trans.T)
            return Nresult,Uvec,NTrans
    return None

def makeBilbaoPhase(result,uvec,trans,ifMag=False):
    phase = {}
    phase['Name'] = result[0].strip()
    phase['Uvec'] = uvec
    phase['Trans'] = trans
    phase['Keep'] = False
    phase['Use'] = False
    phase['aType'] = ''
    SpGp = result[0].replace("'",'')
    SpGrp = G2spc.StandardizeSpcName(SpGp)
    phase['SGData'] = G2spc.SpcGroup(SpGrp)[1]
    if ifMag:
        BNSlatt = phase['SGData']['SGLatt']
        if not result[1]:
            MSpGrp = G2spc.SplitMagSpSG(result[0])
            phase['SGData']['SGSpin'] = G2spc.GetSGSpin(phase['SGData'],MSpGrp)
        phase['SGData']['GenSym'],phase['SGData']['GenFlg'],BNSsym = G2spc.GetGenSym(phase['SGData'])
        if result[1]:
            BNSlatt += '_'+result[1]
            if 'P_S' in BNSlatt: BNSlatt = 'P_c'    #triclinic fix
            phase['SGData']['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
            G2spc.ApplyBNSlatt(phase['SGData'],phase['SGData']['BNSlattsym'])
        phase['SGData']['SpnFlp'] = G2spc.GenMagOps(phase['SGData'])[1]
        phase['SGData']['MagSpGrp'] = G2spc.MagSGSym(phase['SGData'])
    return phase

def FillUnitCell(Phase,Force=True):
    Atoms = copy.deepcopy(Phase['Atoms'])
    atomData = []
    atCodes = []
    SGData = Phase['General']['SGData']
    SpnFlp = SGData.get('SpnFlp',[])
    Amat,Bmat = cell2AB(Phase['General']['Cell'][1:7])
    cx,ct,cs,cia = Phase['General']['AtomPtrs']
    cm = 0
    if Phase['General']['Type'] == 'magnetic':
        cm = cx+4
    for iat,atom in enumerate(Atoms):
        XYZ = np.array(atom[cx:cx+3])
        xyz = XYZ
        cellj = np.zeros(3,dtype=np.int32)
        if Force:
            xyz,cellj = G2spc.MoveToUnitCell(xyz)
        if atom[cia] == 'A':
            Uij = atom[cia+2:cia+8]
            result = G2spc.GenAtom(xyz,SGData,False,Uij,Force)
            for item in result:
                item = list(item)
                item[2] += cellj
#                if item[0][2] >= .95: item[0][2] -= 1.
                atom[cx:cx+3] = item[0]
                atom[cia+2:cia+8] = item[1]
                if cm:
                    Opr = abs(item[2])%100
                    M = SGData['SGOps'][Opr-1][0]
                    opNum = G2spc.GetOpNum(item[2],SGData)
                    mom = np.inner(np.array(atom[cm:cm+3]),Bmat)
                    atom[cm:cm+3] = np.inner(np.inner(mom,M),Amat)*nl.det(M)*SpnFlp[opNum-1]
                atCodes.append('%d:%s'%(iat,str(item[2])))
                atomData.append(atom[:cia+9])  #not SS stuff
        else:
            result = G2spc.GenAtom(xyz,SGData,False,Move=Force)
            for item in result:
                item = list(item)
                item[2] += cellj
#                if item[0][2] >= .95: item[0][2] -= 1.
                atom[cx:cx+3] = item[0]
                if cm:
                    Opr = abs(item[1])%100
                    M = SGData['SGOps'][Opr-1][0]
                    opNum = G2spc.GetOpNum(item[1],SGData)
                    mom = np.inner(np.array(atom[cm:cm+3]),Bmat)
                    atom[cm:cm+3] = np.inner(np.inner(mom,M),Amat)*nl.det(M)*SpnFlp[opNum-1]
                atCodes.append('%d:%s'%(iat,str(item[1])))
                atomData.append(atom[:cia+9])  #not SS stuff

    return atomData,atCodes

def GetUnique(Phase,atCodes):

    def noDuplicate(xyzA,XYZ):
        if True in [np.allclose(xyzA%1.,xyzB%1.,atol=0.0002) for xyzB in XYZ]:
            return False
        return True

    cx,ct = Phase['General']['AtomPtrs'][:2]
    SGData = Phase['General']['SGData']
    Atoms = Phase['Atoms']
    Ind = len(Atoms)
    newAtoms = []
    newAtCodes = []
    Indx = {}
    XYZ = {}
    for ind in range(Ind):
        XYZ[ind] = np.array(Atoms[ind][cx:cx+3])%1.
        Indx[ind] = True
    for ind in range(Ind):
        if Indx[ind]:
            xyz = XYZ[ind]
            for jnd in range(Ind):
                if Atoms[ind][ct-1] == Atoms[jnd][ct-1]:
                    if ind != jnd and Indx[jnd]:
                        Equiv = G2spc.GenAtom(XYZ[jnd],SGData,Move=True)
                        xyzs = np.array([equiv[0] for equiv in Equiv])
                        Indx[jnd] = noDuplicate(xyz,xyzs)
    Ind = []
    for ind in Indx:
        if Indx[ind]:
            newAtoms.append(Atoms[ind])
            newAtCodes.append(atCodes[ind])
    return newAtoms,newAtCodes

def calc_rVsq(A):
    """Compute the square of the reciprocal lattice volume (1/V**2) from A'

    """
    G,g = A2Gmat(A)
    rVsq = nl.det(G)
    if rVsq < 0:
        return 1
    return rVsq

def calc_rV(A):
    """Compute the reciprocal lattice volume (V*) from A
    """
    return np.sqrt(calc_rVsq(A))

def calc_V(A):
    """Compute the real lattice volume (V) from A
    """
    return 1./calc_rV(A)

def A2invcell(A):
    """Compute reciprocal unit cell constants from A
    returns tuple with a*,b*,c*,alpha*, beta*, gamma* (degrees)
    """
    G,g = A2Gmat(A)
    return Gmat2cell(G)

def Gmat2AB(G):
    """Computes orthogonalization matrix from reciprocal metric tensor G

    :returns: tuple of two 3x3 numpy arrays (A,B)

       * A for crystal to Cartesian transformations (A*x = np.inner(A,x) = X)
       * B (= inverse of A) for Cartesian to crystal transformation (B*X = np.inner(B,X) = x)

    """
#    cellstar = Gmat2cell(G)
    g = nl.inv(G)
    cell = Gmat2cell(g)
#    A = np.zeros(shape=(3,3))
    return cell2AB(cell)
#    # from Giacovazzo (Fundamentals 2nd Ed.) p.75
#    A[0][0] = cell[0]                # a
#    A[0][1] = cell[1]*cosd(cell[5])  # b cos(gamma)
#    A[0][2] = cell[2]*cosd(cell[4])  # c cos(beta)
#    A[1][1] = cell[1]*sind(cell[5])  # b sin(gamma)
#    A[1][2] = -cell[2]*cosd(cellstar[3])*sind(cell[4]) # - c cos(alpha*) sin(beta)
#    A[2][2] = 1./cellstar[2]         # 1/c*
#    B = nl.inv(A)
#    return A,B

def cell2AB(cell,alt=False):
    """Computes orthogonalization matrix from unit cell constants

    :param tuple cell: a,b,c, alpha, beta, gamma (degrees)
    :returns: tuple of two 3x3 numpy arrays (A,B)
       A for crystal to Cartesian transformations A*x = np.inner(A,x) = X
       B (= inverse of A) for Cartesian to crystal transformation B*X = np.inner(B,X) = x
       both rounded to 12 places (typically zero terms = +/-10e-6 otherwise)
    """
    G,g = cell2Gmat(cell)
    cellstar = Gmat2cell(G)
    A = np.zeros(shape=(3,3))
    if alt: #as used in RMCProfile!!
        A[0][0] = 1./cellstar[0]
        A[0][1] = cell[0]*cosd(cell[5])*sind(cell[3])
        A[0][2] = cell[0]*cosd(cell[4])
        A[1][1] = cell[1]*sind(cell[3])
        A[1][2] = cell[1]*cosd(cell[3])
        A[2][2] = cell[2]
        A = np.around(A,12)
        B = nl.inv(A)
        return A,B
    # from Giacovazzo (Fundamentals 2nd Ed.) p.75
    A[0][0] = cell[0]                # a
    A[0][1] = cell[1]*cosd(cell[5])  # b cos(gamma)
    A[0][2] = cell[2]*cosd(cell[4])  # c cos(beta)
    A[1][1] = cell[1]*sind(cell[5])  # b sin(gamma)
    A[1][2] = -cell[2]*cosd(cellstar[3])*sind(cell[4]) # - c cos(alpha*) sin(beta)
    A[2][2] = 1./cellstar[2]         # 1/c*
    A = np.around(A,12)
    B = nl.inv(A)
    return A,B

def HKL2SpAng(H,cell,SGData):
    """Computes spherical coords for hkls; view along 001

    :param array H: arrays of hkl
    :param tuple cell: a,b,c, alpha, beta, gamma (degrees)
    :param dict SGData: space group dictionary
    :returns: arrays of r,phi,psi (radius,inclination,azimuth) about 001
    """
    A,B = cell2AB(cell)
    xH = np.inner(B.T,H)
    r = np.sqrt(np.sum(xH**2,axis=0))
    phi = acosd(xH[2]/r)
    psi = atan2d(xH[1],xH[0])
    phi = np.where(phi>90.,180.-phi,phi)
    return r,phi,psi

def U6toUij(U6):
    """Fill matrix (Uij) from U6 = [U11,U22,U33,U12,U13,U23]
    NB: there is a non numpy version in GSASIIspc: U2Uij

    :param list U6: 6 terms of u11,u22,...
    :returns:
        Uij - numpy [3][3] array of uij
    """
    U = np.array([
        [U6[0],  U6[3],  U6[4]],
        [U6[3],  U6[1],  U6[5]],
        [U6[4],  U6[5],  U6[2]]])
    return U

def UijtoU6(U):
    """Fill vector [U11,U22,U33,U12,U13,U23] from Uij
    NB: there is a non numpy version in GSASIIspc: Uij2U
    """
    U6 = np.array([U[0][0],U[1][1],U[2][2],U[0][1],U[0][2],U[1][2]])
    return U6

def betaij2Uij(betaij,G):
    """
    Convert beta-ij to Uij tensors

    :param beta-ij - numpy array [beta-ij]
    :param G: reciprocal metric tensor
    :returns: Uij: numpy array [Uij]
    """
    ast = np.sqrt(np.diag(G))   #a*, b*, c*
    Mast = np.multiply.outer(ast,ast)
    return R2pisq*UijtoU6(U6toUij(betaij)/Mast)

def Uij2betaij(Uij,G):
    """
    Convert Uij to beta-ij tensors -- stub for eventual completion

    :param Uij: numpy array [Uij]
    :param G: reciprocal metric tensor
    :returns: beta-ij - numpy array [beta-ij]
    """
    pass

def cell2GS(cell):
    ''' returns Uij to betaij conversion matrix'''
    G,g = cell2Gmat(cell)
    GS = G
    GS[0][1] = GS[1][0] = math.sqrt(GS[0][0]*GS[1][1])
    GS[0][2] = GS[2][0] = math.sqrt(GS[0][0]*GS[2][2])
    GS[1][2] = GS[2][1] = math.sqrt(GS[1][1]*GS[2][2])
    return GS

def Uij2Ueqv(Uij,GS,Amat):
    ''' returns 1/3 trace of diagonalized U matrix
    :param Uij: numpy array [Uij]
    :param GS: Uij too betaij conversion matrix
    :param Amat: crystal to Cartesian transformation matrix
    :returns: 1/3 trace of diagonalized U matrix
    :returns: True if nonpositive-definite; False otherwise
    '''
    U = np.multiply(U6toUij(Uij),GS)
    U = np.inner(Amat,np.inner(U,Amat).T)
    E,R = nl.eigh(U)
    return np.sum(E)/3.,E[0] < 0.

def CosAngle(U,V,G):
    """ calculate cos of angle between U & V in generalized coordinates
    defined by metric tensor G

    :param U: 3-vectors assume numpy arrays, can be multiple reflections as (N,3) array
    :param V: 3-vectors assume numpy arrays, only as (3) vector
    :param G: metric tensor for U & V defined space assume numpy array
    :returns:
        cos(phi)
    """
    u = (U.T/np.sqrt(np.sum(np.inner(U,G)*U,axis=1))).T
    v = V/np.sqrt(np.inner(V,np.inner(G,V)))
    cosP = np.inner(u,np.inner(G,v))
    return cosP

def CosSinAngle(U,V,G):
    """ calculate sin & cos of angle between U & V in generalized coordinates
    defined by metric tensor G

    :param U: 3-vectors assume numpy arrays
    :param V: 3-vectors assume numpy arrays
    :param G: metric tensor for U & V defined space assume numpy array
    :returns:
        cos(phi) & sin(phi)
    """
    u = U/np.sqrt(np.inner(U,np.inner(G,U)))
    v = V/np.sqrt(np.inner(V,np.inner(G,V)))
    cosP = np.inner(u,np.inner(G,v))
    sinP = np.sqrt(max(0.0,1.0-cosP**2))
    return cosP,sinP

def criticalEllipse(prob):
    """
    Calculate critical values for probability ellipsoids from probability
    """
    if not ( 0.01 <= prob < 1.0):
        return 1.54
    coeff = np.array([6.44988E-09,4.16479E-07,1.11172E-05,1.58767E-04,0.00130554,
        0.00604091,0.0114921,-0.040301,-0.6337203,1.311582])
    llpr = math.log(-math.log(prob))
    return np.polyval(coeff,llpr)

def CellBlock(nCells):
    """
    Generate block of unit cells n*n*n on a side; [0,0,0] centered, n = 2*nCells+1
    currently only works for nCells = 0 or 1 (not >1)
    """
    if nCells:
        N = 2*nCells+1
        N2 = N*N
        N3 = N*N*N
        cellArray = []
        A = np.array(range(N3))
        cellGen = np.array([A//N2-1,A//N%N-1,A%N-1]).T
        for cell in cellGen:
            cellArray.append(cell)
        return cellArray
    else:
        return [0,0,0]

def CellAbsorption(ElList,Volume):
    '''Compute unit cell absorption

    :param dict ElList: dictionary of element contents including mu and
      number of atoms be cell
    :param float Volume: unit cell volume
    :returns: mu-total/Volume
    '''
    muT = 0
    for El in ElList:
        muT += ElList[El]['mu']*ElList[El]['FormulaNo']
    return muT/Volume

#Permutations and Combinations
# Four routines: combinations,uniqueCombinations, selections & permutations
#These taken from Python Cookbook, 2nd Edition. 19.15 p724-726
#
def _combinators(_handle, items, n):
    """ factored-out common structure of all following combinators """
    if n==0:
        yield [ ]
        return
    for i, item in enumerate(items):
        this_one = [ item ]
        for cc in _combinators(_handle, _handle(items, i), n-1):
            yield this_one + cc
def combinations(items, n):
    """ take n distinct items, order matters """
    def skipIthItem(items, i):
        return items[:i] + items[i+1:]
    return _combinators(skipIthItem, items, n)
def uniqueCombinations(items, n):
    """ take n distinct items, order is irrelevant """
    def afterIthItem(items, i):
        return items[i+1:]
    return _combinators(afterIthItem, items, n)
def selections(items, n):
    """ take n (not necessarily distinct) items, order matters """
    def keepAllItems(items, i):
        return items
    return _combinators(keepAllItems, items, n)
def permutations(items):
    """ take all items, order matters """
    return combinations(items, len(items))

#reflection generation routines
#for these: H = [h,k,l]; A is as used in calc_rDsq; G - inv metric tensor, g - metric tensor;
#           cell - a,b,c,alp,bet,gam in A & deg

def Pos2dsp(Inst,pos):
    ''' convert powder pattern position (2-theta or TOF, musec) to d-spacing
    is currently only approximate for EDX data; accurate for others.
    '''
    if 'T' in Inst['Type'][0]:
        return TOF2dsp(Inst,pos)
    elif 'E' in Inst['Type'][0]:
        return 12.398/(2.0*pos*sind(Inst['2-theta'][1]/2.0))
    else:   #'PKS', 'C' or 'B'
        wave = G2mth.getWave(Inst)
        return wave/(2.0*sind((pos-Inst.get('Zero',[0,0])[1])/2.0))

def TOF2dsp(Inst,Pos):
    ''' convert powder pattern TOF, musec to d-spacing by successive approximation
    Pos can be numpy array
    '''
    def func(d,pos,Inst):
        return (pos-Inst['difA'][1]*d**2-Inst['Zero'][1]-Inst['difB'][1]/d)/Inst['difC'][1]
    dsp0 = Pos/Inst['difC'][1]
    N = 0
    while True:      #successive approximations
        dsp = func(dsp0,Pos,Inst)
        if np.allclose(dsp,dsp0,atol=0.000001):
            return dsp
        dsp0 = dsp
        N += 1
        if N > 10:
            return dsp

def Dsp2pos(Inst,dsp):
    ''' convert d-spacing to powder pattern position (2-theta or TOF, musec)
    '''
    if 'T' in Inst['Type'][0]:
        pos = Inst['difC'][1]*dsp+Inst['difA'][1]*dsp**2+Inst.get('difB',[0,0,False])[1]/dsp+Inst['Zero'][1]
    elif 'E' in Inst['Type'][0]:
        return 12.398/(2.0*dsp*sind(Inst['2-theta'][1]/2.0))+Inst['ZE'][1]+Inst['YE'][1]*dsp+Inst['XE'][1]*dsp**2
    else:   #'A', 'B', 'C' or 'PKS'
        wave = G2mth.getWave(Inst)
        val = min(0.995,wave/(2.*dsp))  #set max at 168deg
        pos = 2.0*asind(val)+Inst.get('Zero',[0,0])[1]
    return pos

def getPeakPos(dataType,parmdict,dsp):
    ''' convert d-spacing to powder pattern position (2-theta, E or TOF, musec)
    '''
    if 'T' in dataType:
        pos = parmdict['difC']*dsp+parmdict['difA']*dsp**2+parmdict['difB']/dsp+parmdict['Zero']
    elif 'E'in dataType:
        pos = 12.398/(2.0*dsp*sind(parmdict['2-theta']/2.0)+parmdict['ZE']+parmdict['YE']*dsp+parmdict['XE']*dsp**2)
    else:   #'A', 'B' or 'C'
        pos = 2.0*asind(parmdict['Lam']/(2.*dsp))+parmdict['Zero']
    return pos

def calc_rDsq(H,A):
    'calc 1/d^2 from individual hkl and A-terms'
    B = np.array([H[0]**2,H[1]**2,H[2]**2,H[0]*H[1],H[0]*H[2],H[1]*H[2]])
    return np.sum(A*B)

def calc_rDsqA(H,A):
    'calc array of 1/d^2 from array of hkl & A-terms'
    # B = np.array([H[0]**2,H[1]**2,H[2]**2,H[0]*H[1],H[0]*H[2],H[1]*H[2]])
    # return np.sum(np.array(A)[:,nxs]*B,axis=0)
    h,k,l = H
    return A[0]*h*h+A[1]*k*k+A[2]*l*l+A[3]*h*k+A[4]*h*l+A[5]*k*l #quicker

def calc_rDsq2(H,G):
    'computes 1/d^2 from one hkl & reciprocal metric tensor G'
    return np.inner(H,np.inner(G,H))

def calc_rDsqSS(H,A,vec):
    'computes 1/d^2 from one hklm, reciprocal metric tensor A & k-vector'
    rdsq = calc_rDsq(H[:3]+(H[3]*vec).T,A)
    return rdsq

def calc_rDsqZ(H,A,Z,tth,lam):
    'computes 1/d^2 from hkl array & reciprocal metric tensor A with CW ZERO shift'
    rdsq = calc_rDsqA(H,A)+Z*sind(tth)*2.0*rpd/lam**2
    return rdsq

def calc_rDsqZSS(H,A,vec,Z,tth,lam):
    'computes 1/d^2 from hklm array, reciprocal metric tensor A & k-vector with CW Z shift'
    rdsq = calc_rDsqA(H[:3]+(H[3][:,np.newaxis]*vec).T,A)+Z*sind(tth)*2.0*rpd/lam**2
    return rdsq

def calc_rDsqT(H,A,Z,tof,difC):
    'computes 1/d^2 from hkl array & reciprocal metric tensor A with TOF ZERO shift'
    rdsq = calc_rDsqA(H,A)+Z/difC
    return rdsq

def calc_rDsqTSS(H,A,vec,Z,tof,difC):
    'computes 1/d^2 from hklm array, reciprocal metric tensor A & k-vector with TOF Z shift'
    rdsq = calc_rDsqA(H[:3]+(H[3][:,np.newaxis]*vec).T,A)+Z/difC
    return rdsq

def PlaneIntercepts(Amat,H,phase,stack):
    ''' find unit cell intercepts for a stack of hkl planes
    '''
    Steps = range(-1,2,2)
    if stack:
        Steps = range(-10,10,1)
    Stack = []
    Ux = np.array([[0,0],[1,0],[1,1],[0,1]])
    for step in Steps:
        HX = []
        for i in [0,1,2]:
            if H[i]:
               h,k,l = [(i+1)%3,(i+2)%3,(i+3)%3]
               for j in [0,1,2,3]:
                    hx = [0,0,0]
                    intcpt = ((phase)/360.+step-H[h]*Ux[j,0]-H[k]*Ux[j,1])/H[l]
                    if 0. <= intcpt <= 1.:
                        hx[h] = Ux[j,0]
                        hx[k] = Ux[j,1]
                        hx[l] = intcpt
                        HX.append(hx)
        if len(HX)> 2:
            HX = np.array(HX)
            DX = np.inner(HX-HX[0],Amat)
            D = np.sqrt(np.sum(DX**2,axis=1))
            Dsort = np.argsort(D)
            HX = HX[Dsort]
            DX = DX[Dsort]
            D = D[Dsort]
            DX[1:,:] = DX[1:,:]/D[1:,nxs]
            A = 2.*np.ones(HX.shape[0])
            A[1:] = [np.dot(DX[1],dx) for dx in DX[1:]]
            HX = HX[np.argsort(A)]
            Stack.append(HX)
    return Stack

def MaxIndex(dmin,A):
    'needs doc string'
    Hmax = [0,0,0]
    try:
        cell = A2cell(A)
    except:
        cell = [1.,1.,1.,90.,90.,90.]
    for i in range(3):
        Hmax[i] = int(np.round(cell[i]/dmin))
    return Hmax

def transposeHKLF(transMat,Super,refList):
    ''' Apply transformation matrix to hkl(m)
    param: transmat: 3x3 or 4x4 array
    param: Super: 0 or 1 for extra index
    param: refList list of h,k,l,....
    return: newRefs transformed list of h',k',l',,,
    return: badRefs list of noninteger h',k',l'...
    '''
    newRefs = np.copy(refList)
    badRefs = []
    for H in newRefs:
        newH = np.inner(transMat,H[:3+Super])
        H[:3+Super] = np.rint(newH)
        if not np.allclose(newH,H[:3+Super],atol=0.01):
            badRefs.append(newH)
    return newRefs,badRefs

def sortHKLd(HKLd,ifreverse,ifdup,ifSS=False):
    '''sort reflection list on d-spacing; can sort in either order

    :param HKLd: a list of [h,k,l,d,...];
    :param ifreverse: True for largest d first
    :param ifdup: True if duplicate d-spacings allowed
    :return: sorted reflection list
    '''
    T = []
    N = 3
    if ifSS:
        N = 4
    for i,H in enumerate(HKLd):
        if ifdup:
            T.append((H[N],i))
        else:
            T.append(H[N])
    D = dict(zip(T,HKLd))
    T.sort()
    if ifreverse:
        T.reverse()
    X = []
    okey = ''
    for key in T:
        if key != okey: X.append(D[key])    #remove duplicate d-spacings
        okey = key
    return X

def SwapIndx(Axis,H):
    'needs doc string'
    if Axis in [1,-1]:
        return H
    elif Axis in [2,-3]:
        return [H[1],H[2],H[0]]
    else:
        return [H[2],H[0],H[1]]

def SwapItems(Alist,pos1,pos2):
    'exchange 2 items in a list'
    try:
        get = Alist[pos1],Alist[pos2]
        Alist[pos2],Alist[pos1] = get
    except IndexError:
        pass
    return Alist

def Rh2Hx(Rh):
    'needs doc string'
    Hx = [0,0,0]
    Hx[0] = Rh[0]-Rh[1]
    Hx[1] = Rh[1]-Rh[2]
    Hx[2] = np.sum(Rh)
    return Hx

def Hx2Rh(Hx):
    'needs doc string'
    Rh = [0,0,0]
    itk = -Hx[0]+Hx[1]+Hx[2]
    if itk%3 != 0:
        return 0        #error - not rhombohedral reflection
    else:
        Rh[1] = itk//3
        Rh[0] = Rh[1]+Hx[0]
        Rh[2] = Rh[1]-Hx[1]
        if Rh[0] < 0:
            for i in range(3):
                Rh[i] = -Rh[i]
        return Rh

def CentCheck(Cent,H):
    'checks individual hkl for centering extinction; returns True for allowed, False otherwise - slow'
    h,k,l = H
    if Cent == 'A' and (k+l)%2:
        return False
    elif Cent == 'B' and (h+l)%2:
        return False
    elif Cent == 'C' and (h+k)%2:
        return False
    elif Cent == 'I' and (h+k+l)%2:
        return False
    elif Cent == 'F' and ((h+k)%2 or (h+l)%2 or (k+l)%2):
        return False
    elif Cent == 'R' and (-h+k+l)%3:
        return False
    else:
        return True

def newCentCheck(Cent,H):
    'checks np.array of HKls for centering extinction; returns allowed HKLs - fast'
    K = H.T
    if Cent == 'A':
        return H[np.where((K[1]+K[2])%2 == 0)]
    elif Cent == 'B':
        return H[np.where((K[0]+K[2])%2 == 0)]
    elif Cent == 'C':
        return H[np.where((K[0]+K[1])%2 == 0)]
    elif Cent == 'I':
        return H[np.where((K[0]+K[1]+K[2])%2 == 0)]
    elif Cent == 'F':
        H = H[np.where((K[0]+K[1])%2 == 0)]
        K = H.T
        H = H[np.where((K[0]+K[2])%2 == 0)]
        K = H.T
        return H[np.where((K[1]+K[2])%2 == 0)]
    elif Cent == 'R':
        return H[np.where((-K[0]+K[1]+K[2])%3 == 0)]
    else:
        return H

def RBsymCheck(Atoms,ct,cx,cs,AtLookUp,Amat,RBObjIds,SGData):
    """ Checks members of a rigid body to see if one is a symmetry equivalent of another.
    If so the atom site frac is set to zero.

    :param Atoms: atom array as defined in GSAS-II; modified here
    :param ct: int location of atom type in Atoms item
    :param cx: int location of x,y,z,frac in Atoms item
    :param dict AtLookUp: atom lookup by Id table
    :param np.array Amat: crystal-to-Cartesian transformation matrix
    :param list RBObjIds: atom Id belonging to rigid body being tested
    :param dict SGData: GSAS-II space group info.
    :returns: Atoms with modified atom frac entries

    """
    for i,Id in enumerate(RBObjIds):
        XYZo = np.array(Atoms[AtLookUp[Id]][cx:cx+3])%1.
        typo = Atoms[AtLookUp[Id]][ct]
        for Jd in RBObjIds[i+1:]:
            if Atoms[AtLookUp[Jd]][ct] == typo:
                XYZt = Atoms[AtLookUp[Jd]][cx:cx+3]
                Xeqv = list(G2spc.GenAtom(np.array(XYZt)%1.,SGData,True))
                close = [np.allclose(np.inner(Amat,XYZo),np.inner(Amat,eqv[0]),atol=0.1) for eqv in Xeqv]
                if True in close:
                    Atoms[AtLookUp[Jd]][cx+3] = 0.0
        Sytsym,Mult = G2spc.SytSym(Atoms[AtLookUp[Id]][cx:cx+3],SGData)[:2]
        Atoms[AtLookUp[Id]][cs] = Sytsym
        Atoms[AtLookUp[Id]][cs+1] = Mult
    return Atoms

def GetBraviasNum(center,system):
    """Determine the Bravais lattice number, as used in GenHBravais

    :param center: one of: 'P', 'C', 'I', 'F', 'R' (see SGLatt from GSASIIspc.SpcGroup)
    :param system: one of 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'trigonal' (for R)
      'monoclinic', 'triclinic' (see SGSys from GSASIIspc.SpcGroup)
    :return: a number between 0 and 13
      or throws a ValueError exception if the combination of center, system is not found (i.e. non-standard)

    """
    if center.upper() == 'F' and system.lower() == 'cubic':
        return 0
    elif center.upper() == 'I' and system.lower() == 'cubic':
        return 1
    elif center.upper() == 'P' and system.lower() == 'cubic':
        return 2
    elif center.upper() == 'R' and system.lower() == 'trigonal':
        return 3
    elif center.upper() == 'P' and system.lower() == 'hexagonal':
        return 4
    elif center.upper() == 'I' and system.lower() == 'tetragonal':
        return 5
    elif center.upper() == 'P' and system.lower() == 'tetragonal':
        return 6
    elif center.upper() == 'F' and system.lower() == 'orthorhombic':
        return 7
    elif center.upper() == 'I' and system.lower() == 'orthorhombic':
        return 8
    elif center.upper() == 'A' and system.lower() == 'orthorhombic':
        return 9
    elif center.upper() == 'B' and system.lower() == 'orthorhombic':
        return 10
    elif center.upper() == 'C' and system.lower() == 'orthorhombic':
        return 11
    elif center.upper() == 'P' and system.lower() == 'orthorhombic':
        return 12
    elif center.upper() == 'I' and system.lower() == 'monoclinic':
        return 13
    elif center.upper() == 'A' and system.lower() == 'monoclinic':
        return 14
    elif center.upper() == 'C' and system.lower() == 'monoclinic':
        return 15
    elif center.upper() == 'P' and system.lower() == 'monoclinic':
        return 16
    elif center.upper() == 'P' and system.lower() == 'triclinic':
        return 17
    raise ValueError('non-standard Bravais lattice center=%s, cell=%s' % (center,system))

def _GenHBravais_cctbx(dmin, Bravais, A, sg_type, uctbx_unit_cell, miller_index_generator):
    '''Alternate form of :func:`GenHBravais` that uses CCTBX internals
    '''
    g_inv = np.array([[A[0],   A[3]/2, A[4]/2],
                      [A[3]/2, A[1],   A[5]/2],
                      [A[4]/2, A[5]/2, A[2]]])
    g = np.linalg.inv(g_inv)
    g_elems = (g[0][0], g[1][1], g[2][2], g[0][1], g[0][2], g[1][2])
    try:
        uc = uctbx_unit_cell(metrical_matrix=g_elems)
    except ValueError: # this function sometimes receives an A matrix that gives
                           # numbers <0 in the diagonal elems of g. Not sure why.
        return []
    #if sg_type is None:
    #    sg_type = make_sgtype(Bravais)
    mig = miller_index_generator(uc, sg_type, 0, dmin)
    result = []
    for h,k,l in mig:
        d = uc.d((h,k,l))
        result.append([h, k, l, d, -1])
    result.sort(key=lambda l: l[3], reverse=True)
    return result

def GenHBravais(dmin, Bravais, A, cctbx_args=None,ifList=False):
    """Generate the positionally unique powder diffraction reflections

    :param dmin: minimum d-spacing in A
    :param Bravais: lattice type (see GetBraviasNum). Bravais is one of:

            * 0 F cubic
            * 1 I cubic
            * 2 P cubic
            * 3 R hexagonal (trigonal not rhombohedral)
            * 4 P hexagonal
            * 5 I tetragonal
            * 6 P tetragonal
            * 7 F orthorhombic
            * 8 I orthorhombic
            * 9 A orthorhombic
            * 10 B orthorhombic
            * 11 C orthorhombic
            * 12 P orthorhombic
            * 13 I monoclinic
            * 14 A monoclinic
            * 15 C monoclinic
            * 16 P monoclinic
            * 17 P triclinic

    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :param dict cctbx_args: items defined in CCTBX:

         * 'sg_type': value from cctbx.sgtbx.space_group_type(symmorphic_sgs[ibrav])
         * 'uctbx_unit_cell': pointer to :meth:`cctbx.uctbx.unit_cell`
         * 'miller_index_generator':  pointer to :meth:`cctbx.miller.index_generator`
    :param ifList: if True output is 2D list of HKL; if False (default) output is 2D nd.array of HKLs
         * ifList=False is fast & suitable for indexing routines; ifList=True is needed for graphics.
         * Currently applies only to Triclinic, Monoclinic & Orthorhombic; otherwise output is 2D list.

    :returns: HKL unique d list of [h,k,l,d,-1] sorted with largest d first

    """
    if cctbx_args:
        return _GenHBravais_cctbx(dmin, Bravais, A,
                    cctbx_args['sg_type'], cctbx_args['uctbx_unit_cell'], cctbx_args['miller_index_generator'])
    if Bravais in [9,14]:
        Cent = 'A'
    elif Bravais in [10,]:
        Cent = 'B'
    elif Bravais in [11,15]:
        Cent = 'C'
    elif Bravais in [1,5,8,13]:
        Cent = 'I'
    elif Bravais in [0,7]:
        Cent = 'F'
    elif Bravais in [3]:
        Cent = 'R'
    else:
        Cent = 'P'
    Hmax = MaxIndex(dmin,A)
    HKL = []
    #generate HKL block first - no extinction conditions; obey Laue symmetry

    if Bravais == 17:                       #triclinic

        H0 = np.reshape(np.mgrid[0:1,1:Hmax[1]+1,0:1],(3,-1)).T     #0k0, k>0
        H1 = np.reshape(np.mgrid[0:1,1:Hmax[1]+1,-Hmax[2]:0],(3,-1)).T     #0kl, k>=0,l>0
        H2 = np.reshape(np.mgrid[0:1,0:Hmax[1]+1,1:Hmax[2]+1],(3,-1)).T     #0kl, k>=0,l>0
        H = np.reshape(np.mgrid[1:Hmax[0]+1,-Hmax[1]:Hmax[1]+1,-Hmax[2]:Hmax[2]+1],(3,-1)).T
        H = np.vstack((H0,H1,H2,H))

    elif Bravais in [13,14,15,16]:          #monoclinic - b unique

        H0 = np.reshape(np.mgrid[0:1,1:Hmax[1]+1,0:1],(3,-1)).T     #0k0, k>0
        H1 = np.reshape(np.mgrid[0:1,0:Hmax[1]+1,1:Hmax[2]+1],(3,-1)).T     #0kl, k>=0,l>0
        H = np.reshape(np.mgrid[1:Hmax[0]+1,0:Hmax[1]+1,-Hmax[2]:Hmax[2]+1],(3,-1)).T
        H = np.vstack((H0,H1,H))

    elif Bravais in [7,8,9,10,11,12]:       #orthorhombic

        H0 = np.reshape(np.mgrid[0:1,1:Hmax[1]+1,0:1],(3,-1)).T     #0k0, k>0
        H1 = np.reshape(np.mgrid[0:1,0:Hmax[1]+1,1:Hmax[2]+1],(3,-1)).T     #0kl, k>=0,l>0
        H = np.reshape(np.mgrid[1:Hmax[0]+1,0:Hmax[1]+1,0:Hmax[2]+1],(3,-1)).T
        H = np.vstack((H0,H1,H))

    elif Bravais in [5,6]:                  #tetragonal
        H = []
        for l in range(Hmax[2]+1):
            for k in range(Hmax[1]+1):
                for h in range(k,Hmax[0]+1):
                    if [h,k,l] == [0,0,0]:
                        continue
                    H.append([h,k,l])
        H = np.array(H)


    elif Bravais in [3,4]:                  #hexagonal/trigonal
        H = []
        lmin = 0
        if Bravais == 3: lmin = -Hmax[2]
        for l in range(lmin,Hmax[2]+1):
            for k in range(Hmax[1]+1):
                hmin = k
                if l < 0: hmin += 1
                for h in range(hmin,Hmax[0]+1):
                    if [h,k,l] == [0,0,0]:
                        continue
                    H.append([h,k,l])
        H = np.array(H)


    else:                                   #cubic
        H = []
        for l in range(Hmax[2]+1):
            for k in range(l,Hmax[1]+1):
                for h in range(k,Hmax[0]+1):
                    if [h,k,l] == [0,0,0]:
                        continue
                    H.append([h,k,l])
        H = np.array(H)

    #then enforce centering rules & d >= dmin limit; make array or list as needed

    H = newCentCheck(Cent,H)
    rdsq = nprdsq2d(calc_rDsqA(H.T,A),6)
    if ifList:  #structure needed for plotting
       HKL = np.vstack((H.T,rdsq)).T
       HKL = [[int(h[0]),int(h[1]),int(h[2]),h[3],-1] for h in HKL if h[3] >= dmin]
       return sortHKLd(HKL,True,False)
    else:   #fast for indexing
        HKL = np.vstack((H.T,rdsq,[-1]*len(H))).T[rdsq >= dmin]
        if HKL.shape[0] > 1:
            D = HKL[:,3]
            Indx = np.argsort(D)
            return np.flipud(HKL[Indx])
        else:
            return HKL

def getHKLmax(dmin,SGData,A):
    'finds maximum allowed hkl for given A within dmin'
    SGLaue = SGData['SGLaue']
    if SGLaue in ['3R','3mR']:        #Rhombohedral axes
        Hmax = [0,0,0]
        cell = A2cell(A)
        aHx = cell[0]*math.sqrt(2.0*(1.0-cosd(cell[3])))
        cHx = cell[0]*math.sqrt(3.0*(1.0+2.0*cosd(cell[3])))
        Hmax[0] = Hmax[1] = int(round(aHx/dmin))
        Hmax[2] = int(round(cHx/dmin))
        #print Hmax,aHx,cHx
    else:                           # all others
        Hmax = MaxIndex(dmin,A)
    return Hmax

def GenHLaue(dmin,SGData,A):
    """Generate the crystallographically unique powder diffraction reflections
    for a lattice and Bravais type

    :param dmin: minimum d-spacing
    :param SGData: space group dictionary with at least

        * 'SGLaue': Laue group symbol: one of '-1','2/m','mmm','4/m','6/m','4/mmm','6/mmm', '3m1', '31m', '3', '3R', '3mR', 'm3', 'm3m'
        * 'SGLatt': lattice centering: one of 'P','A','B','C','I','F'
        * 'SGUniq': code for unique monoclinic axis one of 'a','b','c' (only if 'SGLaue' is '2/m') otherwise an empty string

    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :return: HKL = list of [h,k,l,d] sorted with largest d first and is unique
            part of reciprocal space ignoring anomalous dispersion

    """
    import math
    SGLaue = SGData['SGLaue']
    SGLatt = SGData['SGLatt']
    SGUniq = SGData['SGUniq']
    #finds maximum allowed hkl for given A within dmin
    Hmax = getHKLmax(dmin,SGData,A)

    dminsq = 1./(dmin**2)
    HKL = []
    if SGLaue == '-1':                       #triclinic
        for l in range(-Hmax[2],Hmax[2]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                hmin = 0
                if (k < 0) or (k ==0 and l < 0): hmin = 1
                for h in range(hmin,Hmax[0]+1):
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1./math.sqrt(rdsq)])
    elif SGLaue == '2/m':                #monoclinic
        axisnum = 1 + ['a','b','c'].index(SGUniq)
        Hmax = SwapIndx(axisnum,Hmax)
        for h in range(Hmax[0]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                lmin = 0
                if k < 0:lmin = 1
                for l in range(lmin,Hmax[2]+1):
                    [h,k,l] = SwapIndx(-axisnum,[h,k,l])
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1./math.sqrt(rdsq)])
                    [h,k,l] = SwapIndx(axisnum,[h,k,l])
    elif SGLaue in ['mmm','4/m','6/m']:            #orthorhombic
        for l in range(Hmax[2]+1):
            for h in range(Hmax[0]+1):
                kmin = 1
                if SGLaue == 'mmm' or h ==0: kmin = 0
                for k in range(kmin,Hmax[1]+1):
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1./math.sqrt(rdsq)])
    elif SGLaue in ['4/mmm','6/mmm']:                  #tetragonal & hexagonal
        for l in range(Hmax[2]+1):
            for h in range(Hmax[0]+1):
                for k in range(h+1):
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1./math.sqrt(rdsq)])
    elif SGLaue in ['3m1','31m','3','3R','3mR']:                  #trigonals
        for l in range(-Hmax[2],Hmax[2]+1):
            hmin = 0
            if l < 0: hmin = 1
            for h in range(hmin,Hmax[0]+1):
                if SGLaue in ['3R','3']:
                    kmax = h
                    kmin = -int((h-1.)/2.)
                else:
                    kmin = 0
                    kmax = h
                    if SGLaue in ['3m1','3mR'] and l < 0: kmax = h-1
                    if SGLaue == '31m' and l < 0: kmin = 1
                for k in range(kmin,kmax+1):
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if SGLaue in ['3R','3mR']:
                        H = Hx2Rh(H)
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([H[0],H[1],H[2],1./math.sqrt(rdsq)])
    else:                                   #cubic
        for h in range(Hmax[0]+1):
            for k in range(h+1):
                lmin = 0
                lmax = k
                if SGLaue =='m3':
                    lmax = h-1
                    if h == k: lmax += 1
                for l in range(lmin,lmax+1):
                    H = []
                    if CentCheck(SGLatt,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1./math.sqrt(rdsq)])
    return sortHKLd(HKL,True,True)

def GenPfHKLs(nMax,SGData,A):
    """Generate the unique pole figure reflections for a lattice and Bravais type.
    Min d-spacing=1.0A & no more than nMax returned

    :param nMax: maximum number of hkls returned
    :param SGData: space group dictionary with at least

        * 'SGLaue': Laue group symbol: one of '-1','2/m','mmm','4/m','6/m','4/mmm','6/mmm', '3m1', '31m', '3', '3R', '3mR', 'm3', 'm3m'
        * 'SGLatt': lattice centering: one of 'P','A','B','C','I','F'
        * 'SGUniq': code for unique monoclinic axis one of 'a','b','c' (only if 'SGLaue' is '2/m') otherwise an empty string

    :param A: reciprocal metric tensor elements as [G11,G22,G33,2*G12,2*G13,2*G23]
    :return: HKL = list of 'h k l' strings sorted with largest d first; no duplicate zones

    """
    HKL = np.array(GenHLaue(1.0,SGData,A)).T[:3].T     #strip d-spacings
    N = min(nMax,len(HKL))
    return ['%d %d %d'%(h[0],h[1],h[2]) for h in HKL[:N]]

def GenSSHLaue(dmin,SGData,SSGData,Vec,maxH,A):
    'needs a doc string'
    ifMag = False
    if 'MagSpGrp' in SGData:
        ifMag = True
    HKLs = []
    vec = np.array(Vec)
    vstar = np.sqrt(calc_rDsq(vec,A))     #find extra needed for -n SS reflections
    dvec = 1./(maxH*vstar+1./dmin)
    HKL = GenHLaue(dvec,SGData,A)
    SSdH = [vec*h for h in range(-maxH,maxH+1)]
    SSdH = dict(zip(range(-maxH,maxH+1),SSdH))
    for h,k,l,d in HKL:
        ext = G2spc.GenHKLf([h,k,l],SGData)[0]  #h,k,l must be integral values here
        if not ext and d >= dmin:
            HKLs.append([h,k,l,0,d])
        for dH in SSdH:
            if dH:
                DH = SSdH[dH]
                H = [h+DH[0],k+DH[1],l+DH[2]]
                d = 1./np.sqrt(calc_rDsq(H,A))
                if d >= dmin:
                    HKLM = np.array([h,k,l,dH])
                    if (G2spc.checkSSLaue([h,k,l,dH],SGData,SSGData) and G2spc.checkSSextc(HKLM,SSGData)) or ifMag:
                        HKLs.append([h,k,l,dH,d])
    return HKLs

def LaueUnique2(SGData,refList):
    ''' Impose Laue symmetry on hkl

    :param SGData: space group data from 'P '+Laue
    :param HKLF: np.array([[h,k,l,...]]) reflection set to be converted

    :return: HKLF new reflection array with imposed Laue symmetry
    '''
    for ref in refList:
        H = ref[:3]
        Uniq = G2spc.GenHKLf(H,SGData)[2]
        Uniq = G2mth.sortArray(G2mth.sortArray(G2mth.sortArray(Uniq,2),1),0)
        ref[:3] = Uniq[-1]
    return refList

def LaueUnique(Laue,HKLF):
    ''' Impose Laue symmetry on hkl

    :param str Laue: Laue symbol, as below

      centrosymmetric Laue groups::

            ['-1','2/m','112/m','2/m11','mmm','-42m','-4m2','4/mmm','-3','-3m',
            '-31m','-3m1','6/m','6/mmm','m3','m3m']

      noncentrosymmetric Laue groups::

           ['1','2','211','112','m','m11','11m','222','mm2','m2m','2mm',
           '4','-4','422','4mm','3','312','321','3m','31m','3m1','6','-6',
           '622','6mm','-62m','-6m2','23','432','-43m']

    :param HKLF: np.array([[h,k,l,...]]) reflection set to be converted

    :returns: HKLF new reflection array with imposed Laue symmetry
    '''

    HKLFT = HKLF.T
    mat41 = np.array([[0,1,0],[-1,0,0],[0,0,1]])    #hkl -> k,-h,l
    mat43 = np.array([[0,-1,0],[1,0,0],[0,0,1]])    #hkl -> -k,h,l
    mat4bar = np.array([[0,-1,0],[1,0,0],[0,0,-1]]) #hkl -> k,-h,-l
    mat31 = np.array([[-1,-1,0],[1,0,0],[0,0,1]])   #hkl -> ihl = -h-k,h,l
    mat32 = np.array([[0,1,0],[-1,-1,0],[0,0,1]])   #hkl -> kil = k,-h-k,l
    matd3 = np.array([[0,1,0],[0,0,1],[1,0,0]])     #hkl -> k,l,h
    matd3q = np.array([[0,0,-1],[-1,0,0],[0,1,0]])  #hkl -> -l,-h,k
    matd3t = np.array([[0,0,-1],[1,0,0],[0,-1,0]])  #hkl -> -l,h,-k
    mat6 = np.array([[1,1,0],[-1,0,0],[0,0,1]])     #hkl -> h+k,-h,l really 65
    matdm = np.array([[0,1,0],[1,0,0],[0,0,1]])     #hkl -> k,h,l
    matdmp = np.array([[-1,-1,0],[0,1,0],[0,0,1]])  #hkl -> -h-k,k,l
    matkm = np.array([[-1,0,0],[1,1,0],[0,0,1]])    #hkl -> -h,h+k,l
    matd2 = np.array([[0,1,0],[1,0,0],[0,0,-1]])    #hkl -> k,h,-l
    matdm3 = np.array([[1,0,0],[0,0,1],[0,1,0]])    #hkl -> h,l,k
    mat2d43 = np.array([[0,1,0],[1,0,0],[0,0,1]])   #hkl -> k,-h,l
    matk2 = np.array([[-1,0,0],[1,1,0],[0,0,-1]])   #hkl -> -h,-i,-l
    #triclinic
    if Laue == '1': #ok
        pass
    elif Laue == '-1':  #ok
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]<0),HKLFT[:3]*np.array([-1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([-1,-1,-1])[:,nxs],HKLFT[:3])
    #monoclinic
    #noncentrosymmetric - all ok
    elif Laue == '2':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([-1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '1 1 2':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]<0),HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
    elif Laue == '2 1 1':
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
    elif Laue == 'm':
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
    elif Laue == 'm 1 1':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
    elif Laue == '1 1 m':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    #centrosymmetric - all ok
    elif Laue == '2/m 1 1':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]*HKLFT[0]==0)&(HKLFT[1]<0),HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
    elif Laue == '2/m':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]*HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '1 1 2/m':
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]*HKLFT[2]==0)&(HKLFT[0]<0),HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
    #orthorhombic
    #noncentrosymmetric - all OK
    elif Laue == '2 2 2':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == 'm m 2':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
    elif Laue == '2 m m':
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == 'm 2 m':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    #centrosymmetric - all ok
    elif Laue == 'm m m':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    #tetragonal
    #noncentrosymmetric - all ok
    elif Laue == '4':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]>0),np.squeeze(np.inner(HKLF[:,:3],mat41[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-4':
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '4 2 2':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[1]<HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]==0,np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])   #in lieu od 2-fold
    elif Laue == '4 m m':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-4 2 m':
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '-4 m 2':
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[1]<=0),np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]<0),HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[1]==0),np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[0]>HKLFT[1]),np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    #centrosymmetric - all ok
    elif Laue == '4/m':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]>0),np.squeeze(np.inner(HKLF[:,:3],mat41[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '4/m m m':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],mat41[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
    #trigonal - all hex cell
    #noncentrosymmetric - all ok
    elif Laue == '3':
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '3 1 2':
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matk2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],matk2[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '3 2 1':
        HKLFT[:3] = np.where(HKLFT[0]<=-2*HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<-2*HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>0)&(HKLFT[1]==HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T
        HKLFT[:3] = np.where((HKLFT[0]!=0)&(HKLFT[2]>0)&(HKLFT[0]==-2*HKLFT[1]),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '3 1 m':
        HKLFT[:3] = np.where(HKLFT[0]>=HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(2*HKLFT[1]<-HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]>-2*HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matdmp[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T
    elif (Laue == '3 m 1' or Laue == '3 m'):
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]+HKLFT[0])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],matkm[nxs,:,:])).T,HKLFT[:3])
    #centrosymmetric
    elif Laue == '-3':  #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([-1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[0]<0),-np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],-mat31[nxs,:,:])).T,HKLFT[:3])
    elif (Laue == '-3 m 1' or Laue == '-3 m'):  #ok
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]+HKLFT[0])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],matkm[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[1]<HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-3 1 m':  #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([-1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,np.squeeze(np.inner(HKLF[:,:3],-mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    #hexagonal
    #noncentrosymmetric
    elif Laue == '6':   #ok
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]==0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-6':  #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '6 2 2':   #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[0]>HKLFT[1]),np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '6 m m':   #ok
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]==0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]>HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-6 m 2':  #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matk2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat31[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],matk2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '-6 2 m':  #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=-2*HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<-2*HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>0)&(HKLFT[1]==HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T
        HKLFT[:3] = np.where(HKLFT[2]<0,np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]>HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    #centrosymmetric
    elif Laue == '6/m': #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]==0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '6/m m m': #ok
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]+HKLFT[1])<0,np.squeeze(np.inner(HKLF[:,:3],mat32[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,np.squeeze(np.inner(HKLF[:,:3],mat6[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]>HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],matdm.T[nxs,:,:])).T,HKLFT[:3])
    #cubic - all ok
    #noncentrosymmetric -
    elif Laue == '2 3':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]<0)&((HKLFT[0]>-HKLFT[2])|(HKLFT[1]>-HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3t[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]<0)&((HKLFT[0]>-HKLFT[2])|(HKLFT[1]>=-HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3t[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([-1,1,-1])[:,nxs],HKLFT[:3])
    elif Laue == '4 3 2':
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,-1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,np.squeeze(np.inner(HKLF[:,:3],mat43[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]==0)&(HKLFT[1]<HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matd2[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]==0,np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])   #in lieu od 2-fold
        HKLFT[:3] = np.where((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2]),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2]),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]==0,np.squeeze(np.inner(HKLF[:,:3],mat2d43[nxs,:,:])).T,HKLFT[:3])
    elif Laue == '-4 3 m':
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]<=0,HKLFT[:3]*np.array([-1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<=0,np.squeeze(np.inner(HKLF[:,:3],mat4bar[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[1]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<HKLFT[0],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]==0)&(HKLFT[2]<0),HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&(HKLFT[1]<HKLFT[0]),np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([-1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]<0)&(HKLFT[2]<-HKLFT[0])&(HKLFT[1]>HKLFT[2]),np.squeeze(np.inner(HKLF[:,:3],matd3q[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[0]<0)&(HKLFT[2]>=-HKLFT[0])&(HKLFT[1]>HKLFT[2]),np.squeeze(np.inner(HKLF[:,:3],matdm3[nxs,:,:])).T,HKLFT[:3])
    #centrosymmetric
    elif Laue == 'm 3':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
    elif Laue == 'm 3 m':
        HKLFT[:3] = np.where(HKLFT[0]<0,HKLFT[:3]*np.array([-1,1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[1]<0,HKLFT[:3]*np.array([1,-1,1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[2]<0,HKLFT[:3]*np.array([1,1,-1])[:,nxs],HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where((HKLFT[2]>=0)&((HKLFT[0]>=HKLFT[2])|(HKLFT[1]>HKLFT[2])),np.squeeze(np.inner(HKLF[:,:3],matd3[nxs,:,:])).T,HKLFT[:3])
        HKLFT[:3] = np.where(HKLFT[0]>HKLFT[1],np.squeeze(np.inner(HKLF[:,:3],matdm[nxs,:,:])).T,HKLFT[:3])
    return HKLFT.T

#Spherical harmonics routines
def RBChk(sytsym,L,M):
    '''finds symmetry rules for spherical harmonic coefficients for site symmetries
    :param str sytsym: atom site symmetry symbol
    :param int L: principal harmonic term L>0
    :param int M: second harmonic term; can be -L <= M <= L
    :returns True if allowed and sign for term
    NB: not complete for all possible site symmetries! Many are missing
    Based on Tables 2 & 4 of M. Kara & K. Kurki-Suonio, Acta Cryst. A37, 201-210 (1981).
    '''
    if M <= L:
        if sytsym == '53m':
            if not L%2 and M > 0:
                if L in [6,10,12,16,18]:
                    if L%12 == 2:
                        if M <= L//12: return True,1.0
                    else:
                        if M <= L//12+1: return True,1.0
        elif sytsym == '23':   #cubics use different Fourier expansion than those below
            if 2 < L < 11 and [L,M] in [[3,1],[4,1],[6,1],[6,2],[7,1],[8,1],[9,1],[9,2],[10,1],[10,2]]:
                return True,1.0
        elif  sytsym == 'm3':
            if 2 < L < 11 and [L,M] in [[4,1],[6,1],[6,2],[8,1],[10,1],[10,2]]:
                return True,1.0
        elif sytsym == '432':
            if 2 < L < 11 and [L,M] in [[4,1],[6,1],[8,1],[9,2],[10,1]]:
                return True,1.0
        elif sytsym == '-43m':
            if 2 < L < 11 and [L,M] in [[3,1],[4,1],[6,1],[7,1],[8,1],[9,1],[10,1]]:
                return True,1.0
        elif sytsym == 'm3m':   #correct for L < 21 by generator
            if not L%2 and M > 0:
                if L%12 == 2:
                    if M <= L//12: return True,1.0
                else:
                    if M <= L//12+1: return True,1.0
        elif sytsym == '6':
            if not M%6: return True,1.0     #P?
        elif sytsym == '-6':    #L=M+2J
            if L != 1 and not M%3:          #P?
                if not L%2 and not M%6: return True,1.0
                elif L%2 and (M//3)%2: return True,1.0
        elif sytsym == '6/m':
            if not L%2 and not M%6: return True,1.0   #P?
        elif sytsym == '622':
            if not M%6: return True,-1.**M
        elif sytsym == '6mm':
            if not M%6: return True,1.0
        elif sytsym in ['-6m2(100)','-6m2']:   #L=M+2J
            if L != 1 and not M%3:
                if not L%2 and not M%6: return True,1.0
                elif L%2 and (M//3)%2: return True,1.0
        elif sytsym == '-6m2(120)':   #L=M+2J
            if L != 1 and not M%3:
                if not L%2 and not M%6: return True,1.0
                elif L%2 and (M//3)%2: return True,-1.**M
        elif sytsym == '6/mmm':
            if not L%2 and not M%6: return True,1.0
        elif sytsym == '4(z)':
            if not M%4: return True,1.0     #P?
        elif sytsym == '-4(z)':   #m=2l-4j
            if L%2 and (M//2)%2: return True,1.0    #P?
            if not L%2 and not (M//2)%2: return True,1.0
        elif sytsym == '4/m(z)':
            if not M%4: return True,1.0   #P?
        elif sytsym == '422(z)':
            if not M%4: return True,-1.0**L
        elif sytsym == '4mm(z)':
            if not M%4: return True,1.0
        elif sytsym in ['-42m(z)','-42m']:   #m=2l-4j
            if L%2 and (M//2)%2: return True,1.0
            if not L%2 and not (M//2)%2: return True,-1.0**L
        elif sytsym == '-4m2(z)':   #m=2l-4j
            if L%2 and (M//2)%2: return True,1.0
            if not L%2 and not (M//2)%2: return True,1.0
        elif sytsym == '4/mmm(z)':
            if not L%2 and not M%4: return True,1.0
        elif sytsym in ['3','3(111)']:
            if not M%3: return True,1.0     #P?
        elif sytsym in ['-3','-3(111)']:
            if not L%2 and not M%3: return True,1.0    #P?
        elif sytsym in ['32','32(100)','32(111)']:
            if not M%3: return True,-1.0**(L-M)
        elif sytsym == '32(120)':
            if not M%3: return True,-1.0**(L-M)
        elif sytsym in ['3m','3m(100)','3m(111)']:
            if not M%3: return True,-1.0**M
        elif sytsym == '3m(120)':
            if not M%3: return True,1.0
        elif sytsym in ['-3m(100)','-3m(111)','-3m']:
            if not L%2 and not M%3: return True,-1.0**M
        elif sytsym == '-3m(120)':
            if not L%2 and not M%3: return True,1.0
        elif '222' in sytsym:
            if M%2: return True,-1.0**L
        elif 'mm2(x)' in sytsym:  #m=l-2j
            if L%2 and M%2: return True,1.0  #both odd
            if not L%2 and not M%2: return True,1.0    #both even
        elif 'mm2(y)' in sytsym:  #m=l-2j
            if L%2 and M%2: return True,-1.0**L  #both odd
            if not L%2 and not M%2: return True,-1.0**L     #both even
        elif 'mm2(z)' in sytsym:
            if M%2: return True,1.0
        elif 'mmm' in sytsym :
            if not L%2 and not M%2: return True,1.0
        elif sytsym in ['2(x)','2(100)']:
            return True,-1.0**(L-M)
        elif sytsym == '2(y)':
            return True,-1.0**L
        elif sytsym == '2(z)':
            if not M%2: return True,1.0     #P?
        elif sytsym == 'm(x)':
            if not L%2 : return True,-1.0**M
        elif sytsym == 'm(y)':
            return True,1.0
        elif sytsym == 'm(z)':  #m=l-2j
            if L%2 and M%2: return True,1.0       #P?
            if not L%2 and not M%2: return True,1.0         #P?
        elif sytsym == '2/m(x)':
            if not L%2 : return True,-1.0**M
        elif sytsym in ['2/m(y)','2/m']:
            if not L%2: return True,1.0
        elif sytsym == '2/m(z)':
            if not L%2 and not M%2: return True,1.0
        elif sytsym == '1':     #P?
            return True,1.0
        elif sytsym == '-1':    #P?
            if not L%2: return True,1.0
    return False,0.

def RBsymChk(RBsym,cubic,coefNames,L=18):
    '''imposes rigid body symmetry on spherical harmonics terms
    Key problem is noncubic RB symmetries in cubic site symmetries & vice versa.
    :param str RBsym:  molecular point symmetry symbol
    :param bool cubic: True if atom site symmetry is cubic
    :param list coefNames: sp. harm coefficient names to be checked/converted
    :param int L:  maximum spherical harmonic order no. for cubic generation if needed
    '''
    # cubicsigns = {'C(3,1)c':[-1,],'C(4,1)c':[1,1,],'C(6,1)c':[1,-1],'C(6,2)c':[1,-1],'C(7,1)c':[1,-1],'C(8,1)c':[1,1,1],
    #     'C(9,1)c':[1,-1],'C(9,2)c':[1,-1],'C(10,1)c':[1,-1,-1],'C(10,2)c':[1,1,-1]}
    # cubicnames = {'C(3,1)c':['C(3,2)',],'C(4,1)c':['C(4,0)','C(4,4)'],'C(6,1)c':['C(6,0)','C(6,4)'],
    #     'C(6,2)c':['C(6,2)','C(6,6)'],'C(7,1)c':['C(7,2)','C(7,6)'],'C(8,1)c':['C(8,0)','C(8,4)','C(8,8)'],
    #     'C(9,1)c':['C(9,2)','C((9,6)'],'C(9,2)c':['C(9,4)','C(9,8)'],
    #     'C(10,1)c':['C(10,0)','C(10,4)','C(10,8)'],'C(10,2)c':['C(10,2)','C(10,6)','C(10,10)']}
    newNames = []
    newSgns = []
    if cubic:       #sytsym is a cubic site
        if RBsym in ['53m','532']:
            for name in coefNames:
                LM = eval(name[1:-1])
                if LM[0] in [6,10,12,16,18]:
                    newNames.append(name)
                    newSgns.append(1.0)
        elif RBsym in ['m3m','-43m']:   #take all terms?
            for name in coefNames:
                LM = eval(name[1:-1])
                rbChk,sgn = RBChk(RBsym,LM[0],LM[1])
                if rbChk:
                    newNames.append(name)
                    newSgns.append(1.0)
        else:   #RBsym not cubic or icosahedral
            for name in coefNames:  #these are cubic names
                LM = eval(name[1:-1])
                if (LM[0]+LM[1])%2:     #even L odd M or vv
                    if LM[0]%2:
                        M = [4*m for m in range(LM[0]//2)[1:] if 4*m <= LM[0]]
                    else:
                        M = [4*m for m in range(LM[0]//2) if 4*m <= LM[0]]
                else:       #both even or both odd
                    M = [4*m+2 for m in range(LM[0]//2) if 4*m+2 <= LM[0]]
                for m in M:
                    rbChk,sgn = RBChk(RBsym,LM[0],m)
                    if rbChk:
                        newNames.append('C(%d,%d)'%(LM[0],m))
                        newSgns.append(sgn)
    else:
        if RBsym in ['m3m','-43m','53m']:   #force mol. sym. here
            for L in range(L+1):
                cubNames,cubSgns = GenShCoeff(RBsym,L)
                newNames += cubNames
                newSgns += cubSgns
        else:
            for name in coefNames:
                LM = eval(name[1:])
                rbChk,sgn = RBChk(RBsym,LM[0],LM[1])
                if rbChk:
                    newNames.append(name)
                    newSgns.append(sgn)
    return newNames,newSgns

def GenRBCoeff(sytsym,RBsym,L):
    '''imposes rigid body symmetry on spherical harmonics terms
    Key problem is noncubic RB symmetries in cubic site symmetries & vice versa.
    :param str sytsym: atom position site symmetry symbol
    :param str RBsym: molecular point symmetry symbol
    :param int L: spherical harmonic order no.
    :returns list newNames: spherical harmonic term of order L as either C(L,M) or C(L,M)c for cubic terms
    :returns list newSgns: matching coefficient signs as +/- 1.0
    '''
    coefNames = []
    coefSgns = []
    cubic = False
    if sytsym in ['23','m3','432','-43m','m3m']:
        cubic = True
    for iord in range(1,L+1):
        for n in range(-iord,iord+1):
            rbChk,sgn = RBChk(sytsym,iord,n)
            if rbChk:
                if cubic:
                    coefNames.append('C(%d,%d)c'%(iord,n))
                else:
                    coefNames.append('C(%d,%d)'%(iord,n))
                coefSgns.append(sgn)
    if RBsym == '1':
        return coefNames,coefSgns
    newNames,newSgns = RBsymChk(RBsym,cubic,coefNames,L)
    return newNames,newSgns

def GenShCoeff(sytsym,L):
    '''Generate spherical harmonic coefficient names for atom site symmetry
    :param str sytsym: site symmetry or perhaps molecular symmetry
    :param int L:spherical harmonic order no.
    :returns list newNames: spherical harmonic term of order L as either C(L,M) or C(L,M)c for cubic terms
    :returns list newSgns: matching coefficient signs as +/- 1.0
    '''
    coefNames = []
    coefSgns = []
    cubic = False
    if sytsym in ['23','m3','432','-43m','m3m','53m']:
        cubic = True
    for n in range(-L,L+1):
        rbChk,sgn = RBChk(sytsym,L,n)
        if rbChk:
            if cubic:
                coefNames.append('C(%d,%d)c'%(L,n))
            else:
                coefNames.append('C(%d,%d)'%(L,n))
            coefSgns.append(sgn)
    newNames,newSgns = RBsymChk(sytsym,cubic,coefNames,L)
    return newNames,newSgns

def OdfChk(SGLaue,L,M):
    '''finds symmetry rules for spherical harmonic coefficients for Laue groups
    :param str SGLaue: Laue symbol
    :param int L: principal harmonic term; only evens are used
    :param int M: second harmonic term; can be -L <= M <= L
    :returns True if allowed
    '''
    if not L%2 and abs(M) <= L:
        if SGLaue == '0':                      #cylindrical symmetry
            if M == 0: return True
        elif SGLaue == '-1':
            return True
        elif SGLaue == '2/m':
            if not abs(M)%2: return True
        elif SGLaue == 'mmm':
            if not abs(M)%2 and M >= 0: return True
        elif SGLaue == '4/m':
            if not abs(M)%4: return True
        elif SGLaue == '4/mmm':
            if not abs(M)%4 and M >= 0: return True
        elif SGLaue in ['3R','3']:
            if not abs(M)%3: return True
        elif SGLaue in ['3mR','3m1','31m']:
            if not abs(M)%3 and M >= 0: return True
        elif SGLaue == '6/m':
            if not abs(M)%6: return True
        elif SGLaue == '6/mmm':
            if not abs(M)%6 and M >= 0: return True
        elif SGLaue in ['m3']:   #cubics use different Fourier expansion than those above
            if M > 0:
                if L%12 == 2:
                    if M <= L//12: return True
                else:
                    if M <= L//12+1: return True
        elif SGLaue in ['m3m']:
            if M > 0:
                if L%12 == 2:
                    if M <= L//12: return True
                else:
                    if M <= L//12+1: return True
    return False

def GenSHCoeff(SGLaue,SamSym,L,IfLMN=True):
    '''Generate spherical harmonics coefficient names for texture
    :param str SGLaue: Laue symbol
    :param str SamSym: sample symmetry symbol
    :param int L: spherical harmonic order no.
    :param bool IfLMN: if TRUE return sp.harm. name as C(L,M,N); else return C(L,N)
    :returns coefficient name as C(L,M,N) or C(L,N)
    '''
    coeffNames = []
    for iord in [2*i+2 for i in range(L//2)]:
        for m in [i-iord for i in range(2*iord+1)]:
            if OdfChk(SamSym,iord,m):
                for n in [i-iord for i in range(2*iord+1)]:
                    if OdfChk(SGLaue,iord,n):
                        if IfLMN:
                            coeffNames.append('C(%d,%d,%d)'%(iord,m,n))
                        else:
                            coeffNames.append('C(%d,%d)'%(iord,n))
    return coeffNames

def CrsAng(H,cell,SGData):
    '''Convert HKL to polar coordinates with proper orientation WRT space group point group
    :param array H: hkls
    :param list cell: lattice parameters
    :param dict SGData: space group data
    :returns arrays phi,beta: polar, azimuthal angles for HKL
    '''
    a,b,c,al,be,ga = cell
    SQ3 = 1.732050807569
    H1 = np.array([1,0,0])
    H2 = np.array([0,1,0])
    H3 = np.array([0,0,1])
    H4 = np.array([1,1,1])
    G,g = cell2Gmat(cell)
    Laue = SGData['SGLaue']
    Naxis = SGData['SGUniq']
    if len(H.shape) == 1:
        DH = np.inner(H,np.inner(G,H))
    else:
        DH = np.array([np.inner(h,np.inner(G,h)) for h in H])
    if Laue == '2/m':
        if Naxis == 'a':
            DR = np.inner(H1,np.inner(G,H1))
            DHR = np.inner(H,np.inner(G,H1))
        elif Naxis == 'b':
            DR = np.inner(H2,np.inner(G,H2))
            DHR = np.inner(H,np.inner(G,H2))
        else:
            DR = np.inner(H3,np.inner(G,H3))
            DHR = np.inner(H,np.inner(G,H3))
    elif Laue in ['R3','R3m']:
        DR = np.inner(H4,np.inner(G,H4))
        DHR = np.inner(H,np.inner(G,H4))
    else:
        DR = np.inner(H3,np.inner(G,H3))
        DHR = np.inner(H,np.inner(G,H3))
    DHR /= np.sqrt(DR*DH)
    phi = acosd(max(-1,min(DHR,1.)))
    if Laue == '-1':
        BA = H.T[1]*a/(b-H.T[0]*cosd(ga))
        BB = H.T[0]*sind(ga)**2
    elif Laue == '2/m':
        if Naxis == 'a':
            BA = H.T[2]*b/(c-H.T[1]*cosd(al))
            BB = H.T[1]*sind(al)**2
        elif Naxis == 'b':
            BA = H.T[0]*c/(a-H.T[2]*cosd(be))
            BB = H.T[2]*sind(be)**2
        else:
            BA = H.T[1]*a/(b-H.T[0]*cosd(ga))
            BB = H.T[0]*sind(ga)**2
    elif Laue in ['mmm','4/m','4/mmm']:
        BA = H.T[1]*a
        BB = H.T[0]*b
    elif Laue in ['3R','3mR']:
        BA = H.T[0]+H.T[1]-2.0*H.T[2]
        BB = SQ3*(H.T[0]-H.T[1])
    elif Laue in ['m3','m3m']:
        BA = H.T[1]
        BB = H.T[0]
    else:
        BA = H.T[0]+2.0*H.T[1]
        BB = SQ3*H.T[0]
    beta = atan2d(BA,BB)
    return phi,beta

def SamAng(Tth,Gangls,Sangl,IFCoup):
    """Compute sample orientation angles vs laboratory coord. system

    :param Tth:        Signed theta
    :param Gangls:     Sample goniometer angles phi,chi,omega,azmuth
    :param Sangl:      Sample angle zeros om-0, chi-0, phi-0
    :param IFCoup:     True if omega & 2-theta coupled in CW scan
    :returns:
        psi,gam:    Sample odf angles
        dPSdA,dGMdA:    Angle zero derivatives
    """

    if IFCoup:
        GSomeg = sind(Gangls[2]+Tth)
        GComeg = cosd(Gangls[2]+Tth)
    else:
        GSomeg = sind(Gangls[2])
        GComeg = cosd(Gangls[2])
    GSTth = sind(Tth)
    GCTth = cosd(Tth)
    GSazm = sind(Gangls[3])
    GCazm = cosd(Gangls[3])
    GSchi = sind(Gangls[1])
    GCchi = cosd(Gangls[1])
    GSphi = sind(Gangls[0]+Sangl[2])
    GCphi = cosd(Gangls[0]+Sangl[2])
    SSomeg = sind(Sangl[0])
    SComeg = cosd(Sangl[0])
    SSchi = sind(Sangl[1])
    SCchi = cosd(Sangl[1])
    AT = -GSTth*GComeg+GCTth*GCazm*GSomeg
    BT = GSTth*GSomeg+GCTth*GCazm*GComeg
    CT = -GCTth*GSazm*GSchi
    DT = -GCTth*GSazm*GCchi

    BC1 = -AT*GSphi+(CT+BT*GCchi)*GCphi
    BC2 = DT-BT*GSchi
    BC3 = AT*GCphi+(CT+BT*GCchi)*GSphi

    BC = BC1*SComeg*SCchi+BC2*SComeg*SSchi-BC3*SSomeg
    psi = acosd(BC)

    BD = 1.0-BC**2
    C = np.where(BD>1.e-6,rpd/np.sqrt(BD),0.)
    dPSdA = [-C*(-BC1*SSomeg*SCchi-BC2*SSomeg*SSchi-BC3*SComeg),
        -C*(-BC1*SComeg*SSchi+BC2*SComeg*SCchi),
        -C*(-BC1*SSomeg-BC3*SComeg*SCchi)]

    BA = -BC1*SSchi+BC2*SCchi
    BB = BC1*SSomeg*SCchi+BC2*SSomeg*SSchi+BC3*SComeg
    gam = atan2d(BB,BA)

    BD = (BA**2+BB**2)/rpd

    dBAdO = 0
    dBAdC = -BC1*SCchi-BC2*SSchi
    dBAdF = BC3*SSchi

    dBBdO = BC1*SComeg*SCchi+BC2*SComeg*SSchi-BC3*SSomeg
    dBBdC = -BC1*SSomeg*SSchi+BC2*SSomeg*SCchi
    dBBdF = BC1*SComeg-BC3*SSomeg*SCchi

    dGMdA = np.where(BD > 1.e-6,[(BA*dBBdO-BB*dBAdO)/BD,(BA*dBBdC-BB*dBAdC)/BD, \
        (BA*dBBdF-BB*dBAdF)/BD],[np.zeros_like(BD),np.zeros_like(BD),np.zeros_like(BD)])

    return psi,gam,dPSdA,dGMdA

BOH = {
'L=2':[[],[],[]],
'L=4':[[0.30469720,0.36418281],[],[]],
'L=6':[[-0.14104740,0.52775103],[],[]],
'L=8':[[0.28646862,0.21545346,0.32826995],[],[]],
'L=10':[[-0.16413497,0.33078546,0.39371345],[],[]],
'L=12':[[0.26141975,0.27266871,0.03277460,0.32589402],
    [0.09298802,-0.23773812,0.49446631,0.0],[]],
'L=14':[[-0.17557309,0.25821932,0.27709173,0.33645360],[],[]],
'L=16':[[0.24370673,0.29873515,0.06447688,0.00377,0.32574495],
    [0.12039646,-0.25330128,0.23950998,0.40962508,0.0],[]],
'L=18':[[-0.16914245,0.17017340,0.34598142,0.07433932,0.32696037],
    [-0.06901768,0.16006562,-0.24743528,0.47110273,0.0],[]],
'L=20':[[0.23067026,0.31151832,0.09287682,0.01089683,0.00037564,0.32573563],
    [0.13615420,-0.25048007,0.12882081,0.28642879,0.34620433,0.0],[]],
'L=22':[[-0.16109560,0.10244188,0.36285175,0.13377513,0.01314399,0.32585583],
    [-0.09620055,0.20244115,-0.22389483,0.17928946,0.42017231,0.0],[]],
'L=24':[[0.22050742,0.31770654,0.11661736,0.02049853,0.00150861,0.00003426,0.32573505],
    [0.13651722,-0.21386648,0.00522051,0.33939435,0.10837396,0.32914497,0.0],
    [0.05378596,-0.11945819,0.16272298,-0.26449730,0.44923956,0.0,0.0]],
'L=26':[[-0.15435003,0.05261630,0.35524646,0.18578869,0.03259103,0.00186197,0.32574594],
    [-0.11306511,0.22072681,-0.18706142,0.05439948,0.28122966,0.35634355,0.0],[]],
'L=28':[[0.21225019,0.32031716,0.13604702,0.03132468,0.00362703,0.00018294,0.00000294,0.32573501],
    [0.13219496,-0.17206256,-0.08742608,0.32671661,0.17973107,0.02567515,0.32619598,0.0],
    [0.07989184,-0.16735346,0.18839770,-0.20705337,0.12926808,0.42715602,0.0,0.0]],
'L=30':[[-0.14878368,0.01524973,0.33628434,0.22632587,0.05790047,0.00609812,0.00022898,0.32573594],
    [-0.11721726,0.20915005,-0.11723436,-0.07815329,0.31318947,0.13655742,0.33241385,0.0],
    [-0.04297703,0.09317876,-0.11831248,0.17355132,-0.28164031,0.42719361,0.0,0.0]],
'L=32':[[0.20533892,0.32087437,0.15187897,0.04249238,0.00670516,0.00054977,0.00002018,0.00000024,0.32573501],
    [0.12775091,-0.13523423,-0.14935701,0.28227378,0.23670434,0.05661270,0.00469819,0.32578978,0.0],
    [0.09703829,-0.19373733,0.18610682,-0.14407046,0.00220535,0.26897090,0.36633402,0.0,0.0]],
'L=34':[[-0.14409234,-0.01343681,0.31248977,0.25557722,0.08571889,0.01351208,0.00095792,0.00002550,0.32573508],
    [-0.11527834,0.18472133,-0.04403280,-0.16908618,0.27227021,0.21086614,0.04041752,0.32688152,0.0],
    [-0.06773139,0.14120811,-0.15835721,0.18357456,-0.19364673,0.08377174,0.43116318,0.0,0.0]]
}

Lnorm = lambda L: 4.*np.pi/(2.0*L+1.)

def GetKcl(L,N,SGLaue,phi,beta):
    'needs doc string'
#    from . import pytexture as ptx
    if SGLaue in ['m3','m3m']:
        if 'array' in str(type(phi)) and np.any(phi.shape):
            Kcl = np.zeros_like(phi)
        else:
            Kcl = 0.
        for j in range(0,L+1,4):
            im = j//4
            if 'array' in str(type(phi)) and np.any(phi.shape):
                pcrs = ptx.pyplmpsi(L,j,len(phi),phi)[0]
            else:
                pcrs = ptx.pyplmpsi(L,j,1,phi)[0]
            Kcl += BOH['L=%d'%(L)][N-1][im]*pcrs*cosd(j*beta)
    else:
        if 'array' in str(type(phi)) and np.any(phi.shape):
            pcrs = ptx.pyplmpsi(L,N,len(phi),phi)[0]
        else:
            pcrs = ptx.pyplmpsi(L,N,1,phi)[0]
        pcrs *= RSQ2PI
        if N:
            pcrs *= SQ2
        if SGLaue in ['mmm','4/mmm','6/mmm','R3mR','3m1','31m']:
            if SGLaue in ['3mR','3m1','31m']:
                if N%6 == 3:
                    Kcl = pcrs*sind(N*beta)
                else:
                    Kcl = pcrs*cosd(N*beta)
            else:
                Kcl = pcrs*cosd(N*beta)
        else:
            Kcl = pcrs*(cosd(N*beta)+sind(N*beta))
    return Kcl

def GetKsl(L,M,SamSym,psi,gam):
    'needs doc string'
#    from . import pytexture as ptx
    if 'array' in str(type(psi)) and np.any(psi.shape):
        psrs,dpdps = ptx.pyplmpsi(L,M,len(psi),psi)
    else:
        psrs,dpdps = ptx.pyplmpsi(L,M,1,psi)
    psrs *= RSQ2PI
    dpdps *= RSQ2PI
    if M:
        psrs *= SQ2
        dpdps *= SQ2
    if SamSym in ['mmm',]:
        dum = cosd(M*gam)
        Ksl = psrs*dum
        dKsdp = dpdps*dum
        dKsdg = -psrs*M*sind(M*gam)
    else:
        dum = cosd(M*gam)+sind(M*gam)
        Ksl = psrs*dum
        dKsdp = dpdps*dum
        dKsdg = psrs*M*(-sind(M*gam)+cosd(M*gam))
    return Ksl,dKsdp,dKsdg

def GetKclKsl(L,N,SGLaue,psi,phi,beta):
    """
    This is used for spherical harmonics description of preferred orientation;
        cylindrical symmetry only (M=0) and no sample angle derivatives returned
    """
#    from . import pytexture as ptx
    Ksl,x = ptx.pyplmpsi(L,0,1,psi)
    Ksl *= RSQ2PI
    if SGLaue in ['m3','m3m']:
        Kcl = 0.0
        for j in range(0,L+1,4):
            im = j//4
            pcrs,dum = ptx.pyplmpsi(L,j,1,phi)
            Kcl += BOH['L=%d'%(L)][N-1][im]*pcrs*cosd(j*beta)
    else:
        pcrs,dum = ptx.pyplmpsi(L,N,1,phi)
        pcrs *= RSQ2PI
        if N:
            pcrs *= SQ2
        if SGLaue in ['mmm','4/mmm','6/mmm','R3mR','3m1','31m']:
            if SGLaue in ['3mR','3m1','31m']:
                if N%6 == 3:
                    Kcl = pcrs*sind(N*beta)
                else:
                    Kcl = pcrs*cosd(N*beta)
            else:
                Kcl = pcrs*cosd(N*beta)
        else:
            Kcl = pcrs*(cosd(N*beta)+sind(N*beta))
    return Kcl*Ksl,Lnorm(L)

def H2ThPh2(H,Bmat):
    '''Convert HKL to spherical polar & azimuth angles

    :param array H: array of hkl as [n,3]
    :param [3,3] array Bmat: inv crystal to Cartesian transformation
    :returns array Th: HKL azimuth angles
    :returns array Ph: HKL polar angles
    '''
    Hcart = np.inner(H,Bmat)
    R = np.sqrt(np.sum(np.square(Hcart),axis=1))
    Hcart /= R[:,nxs]
    Pl = acosd(Hcart[:,2])
    Az = atan2d(Hcart[:,1],Hcart[:,0])
    return R,Az,Pl

# def H2ThPh(H,Bmat,Q):
#     '''Convert HKL to spherical polar & azimuth angles - wrong!

#     :param array H: array of hkl as [n,3]
#     :param [3,3] array Bmat: inv crystal to Cartesian transformation
#     :param array Q: quaternion for rotation of HKL to new polar axis
#     :returns array Th: HKL azimuth angles
#     :returns array Ph: HKL polar angles
#     '''
#     # A,V = G2mth.Q2AVdeg(Q)
#     # QR,R = G2mth.make2Quat(V,np.array([0.,0.,1.0]))
#     # QA = G2mth.AVdeg2Q(A,np.array([0.,0.,1.0]))
#     # Q2 = G2mth.prodQQ(QR,QA)
#     Qmat = G2mth.Q2Mat(Q)
#     CH1 = np.inner(H,Bmat.T)
#     CH = np.inner(CH1,Qmat.T)
#     N = nl.norm(CH,axis=1)
#     CH /= N[:,nxs]
#     H3 = np.array([0,0,1.])
#     DHR = np.inner(CH,H3)
#     Ph = np.where(DHR <= 1.0,acosd(DHR),0.0)    #polar angle 0<=Ph<=180.; correct
#     TH = CH*np.array([1.,1.,0.])[nxs,:]     #projection of CH onto xy plane
#     N = nl.norm(TH,axis=1)
#     N = np.where(N > 1.e-5,N,1.)
#     TH /= N[:,nxs]
#     Th = atan2d(TH[:,1],TH[:,0])                #azimuth angle 0<=Th<360<
#     Th = np.where(Th<0.,Th+360.,Th)
#     return Th,Ph        #azimuth,polar angles

def SetUVvec(Neigh):
    ''' Set deformation coordinate choices from neighbors; called in G2phsGUI/UpdateDeformation
    
    param: list neigh: list of neighboring atoms; each with name, dist & cartesian vector

    :returns list UVvec: list of normalized vectors
    :returns list UVchoice: list of names for each
    '''
    UVvec = []
    UVchoice = []    
    choice = ['A','B','C','D','E','F','G','H']
    for N,neigh in enumerate(Neigh):
        UVvec += [neigh[3]/neigh[2],]
        UVchoice += choice[N]
    Nneigh = len(Neigh)
    if Nneigh >= 2:
        UVvec += [(Neigh[0][3]+Neigh[1][3])/np.sqrt(Neigh[0][2]**2+Neigh[1][2]**2),]
        UVchoice += ['A+B',]
    if Nneigh >= 3:
        UVvec += [(Neigh[0][3]+Neigh[1][3]+Neigh[2][3])/np.sqrt(Neigh[0][2]**2+Neigh[1][2]**2+Neigh[2][2]**2),]
        UVchoice += ['A+B+C',]
    return UVvec,UVchoice

def SHarmcal(SytSym,SHFln,psi,gam):
    '''Perform a surface spherical harmonics computation.
    Presently only used for plotting
    Note that the the number of gam values must either be 1 or must match psi

    :param str SytSym: sit symmetry - only looking for cubics - remove this
    :param dict SHFln: spherical harmonics coefficients; key has L & M
    :param float/array psi: Azimuthal coordinate 0 <= Th <= 360
    :param float/array gam: Polar coordinate 0<= Ph <= 180

    :returns array SHVal: spherical harmonics array for psi,gam values
    '''
    SHVal = np.ones_like(psi)/(4.*np.pi)
    for term in SHFln:
        trm = term.strip('+').strip('-')    #patch
        if 'C(' in term[:3]:
            l,m = eval(trm.strip('C').strip('c'))
            if SytSym in ['m3m','m3','43m','432','23'] or 'c' in trm:
                Ksl = CubicSHarm(l,m,psi,gam)
            else:
                p = SHFln[term][2]
                Ksl = SphHarmAng(l,m,p,psi,gam)
            SHVal += SHFln[term][0]*Ksl
    return SHVal

def KslCalc(trm,psi,gam):
    '''Compute one angular part term in spherical harmonics
    :param str trm:sp. harm term name in the form of 'C(l,m)' or 'C(l,m)c' for cubic
    :param float/array psi: Azimuthal coordinate 0 <= Th <= 360
    :param float/array gam: Polar coordinate 0<= Ph <= 180

    :returns array Ksl: spherical harmonics angular part for psi,gam pairs
    '''
    l,m = eval(trm.strip('C').strip('c'))
    if 'c' in trm:
        return CubicSHarm(l,m,psi,gam)
    else:
        return SphHarmAng(l,m,1.0,psi,gam)

def SphHarmAng(L,M,P,Th,Ph):
    ''' Compute spherical harmonics values using scipy.special.sph_harm

    :param int L: degree of the harmonic (L >= 0)
    :param int M: order number (\\|M\\| <= L)
    :param int P: sign flag = -1 or 1
    :param float/array Th: Azimuthal coordinate 0 <= Th <= 360
    :param float/array Ph: Polar coordinate 0<= Ph <= 180

    :returns ylmp value/array: as reals
    '''
    ylmp = spsp.sph_harm(M,L,rpd*Th,rpd*Ph)   #wants radians; order then degree

    if M > 0:
        return (-1)**M*P*np.real(ylmp)*SQ2
    elif M == 0:
        return P*np.real(ylmp)
    else:
        return (-1)**M*P*np.imag(ylmp)*SQ2

def CubicSHarm(L,M,Th,Ph):
    '''Calculation of the cubic harmonics given in Table 3 in M.Kara & K. Kurki-Suonio,
    Acta Cryt. A37, 201 (1981). For L = 14,20 only for m3m from F.M. Mueller and M.G. Priestley,
    Phys Rev 148, 638 (1966)

    :param int L: degree of the harmonic (L >= 0)
    :param int M: order number [\\|M\\| <= L]
    :param float/array Th: Azimuthal coordinate 0 <= Th <= 360
    :param float/array Ph: Polar coordinate 0<= Ph <= 180

    :returns klm value/array: cubic harmonics

    '''
    if L == 0:
        return SphHarmAng(L,M,1,Th,Ph)
    elif L == 3:
        return SphHarmAng(3,2,-1,Th,Ph)
    elif L == 4:
        klm = 0.5*np.sqrt(7.0/3.0)*SphHarmAng(4,0,1,Th,Ph)
        klm += 0.5*np.sqrt(5.0/3.0)*SphHarmAng(4,4,1,Th,Ph)
    elif L == 6:
        if M == 1:
            klm = 0.5*np.sqrt(0.5)*SphHarmAng(6,0,1,Th,Ph)
            klm -= 0.5*np.sqrt(7.0/2.0)*SphHarmAng(6,4,1,Th,Ph)
        else:
            klm = 0.25*np.sqrt(11.0)*SphHarmAng(6,2,1,Th,Ph)
            klm -= 0.25*np.sqrt(5.0)*SphHarmAng(6,6,1,Th,Ph)
    elif L == 7:
        klm = 0.5*np.sqrt(13./6.)*SphHarmAng(7,2,-1,Th,Ph)
        klm += 0.5*np.sqrt(11./6.)*SphHarmAng(7,6,-1,Th,Ph)
    elif L == 8:
        klm = 0.125*np.sqrt(33.)*SphHarmAng(8,0,1,Th,Ph)
        klm += 0.25*np.sqrt(7./3.)*SphHarmAng(8,4,1,Th,Ph)
        klm += 0.125*np.sqrt(65./3.)*SphHarmAng(8,8,1,Th,Ph)
    elif L == 9:
        if M == 1:
            klm = 0.25*np.sqrt(3.)*SphHarmAng(9,2,-1,Th,Ph)
            klm -= 0.25*np.sqrt(13.)*SphHarmAng(9,6,-1,Th,Ph)
        else:
            klm = 0.5*np.sqrt(17./6.)*SphHarmAng(9,4,-1,Th,Ph)
            klm -= 0.5*np.sqrt(7./6.)*SphHarmAng(9,8,-1,Th,Ph)
    elif L == 10:
        if M == 1:
            klm = 0.125*np.sqrt(65./6.)*SphHarmAng(10,0,1,Th,Ph)
            klm -= 0.25*np.sqrt(11.0/2.0)*SphHarmAng(10,4,1,Th,Ph)
            klm -= 0.125*np.sqrt(187.0/6.0)*SphHarmAng(10,8,1,Th,Ph)
        else:
            klm = 0.125*np.sqrt(247./6.)*SphHarmAng(10,2,1,Th,Ph)
            klm += (1./16.)*np.sqrt(19./3.)*SphHarmAng(10,6,1,Th,Ph)
            klm -= (1./16.)*np.sqrt(85.)*SphHarmAng(10,10,1,Th,Ph)
#only m3m cubics from here down; from F.M. Mueller and M.G. Priestley, Phys Rev 148, 638 (1966)
    elif L == 12:
        if M == 1:
            klm = 0.69550265*SphHarmAng(12,0,1,Th,Ph)
            klm += 0.31412555*SphHarmAng(12,4,1,Th,Ph)
            klm += 0.34844954*SphHarmAng(12,8,1,Th,Ph)
            klm += 0.54422797*SphHarmAng(12,12,1,Th,Ph)
        else:
            klm = 0.55897937*SphHarmAng(12,4,1,Th,Ph)
            klm -= 0.80626751*SphHarmAng(12,8,1,Th,Ph)
            klm += 0.19358400*SphHarmAng(12,12,1,Th,Ph)
    elif L == 14:
        klm = 0.44009645*SphHarmAng(14,0,1,Th,Ph)
        klm -= 0.45768183*SphHarmAng(14,4,1,Th,Ph)
        klm -= 0.49113230*SphHarmAng(14,8,1,Th,Ph)
        klm -= 0.59634848*SphHarmAng(14,12,1,Th,Ph)
    elif L == 16:
        if M == 1:
            klm = 0.68136168*SphHarmAng(16,0,1,Th,Ph)
            klm += 0.27586801*SphHarmAng(16,4,1,Th,Ph)
            klm += 0.29048987*SphHarmAng(16,8,1,Th,Ph)
            klm += 0.32756975*SphHarmAng(16,12,1,Th,Ph)
            klm += 0.51764542*SphHarmAng(16,16,1,Th,Ph)
        else:
            klm = 0.63704821*SphHarmAng(16,4,1,Th,Ph)
            klm -= 0.32999033*SphHarmAng(16,8,1,Th,Ph)
            klm -= 0.64798073*SphHarmAng(16,12,1,Th,Ph)
            klm += 0.25572816*SphHarmAng(16,16,1,Th,Ph)
    elif L == 18:
        if M == 1:
            klm = 0.45791513*SphHarmAng(18,0,1,Th,Ph)
            klm -= 0.38645598*SphHarmAng(18,4,1,Th,Ph)
            klm -= 0.40209462*SphHarmAng(18,8,1,Th,Ph)
            klm -= 0.43746593*SphHarmAng(18,12,1,Th,Ph)
            klm -= 0.53657149*SphHarmAng(18,16,1,Th,Ph)
        else:
            klm = 0.14872751*SphHarmAng(18,4,1,Th,Ph)
            klm -= 0.63774601*SphHarmAng(18,8,1,Th,Ph)
            klm += 0.72334167*SphHarmAng(18,12,1,Th,Ph)
            klm -= 0.21894515*SphHarmAng(18,16,1,Th,Ph)
    elif L == 20:
        if M == 1:
            klm = 0.67141495*SphHarmAng(20,0,1,Th,Ph)
            klm += 0.24982619*SphHarmAng(20,4,1,Th,Ph)
            klm += 0.25782846*SphHarmAng(20,8,1,Th,Ph)
            klm += 0.27469333*SphHarmAng(20,12,1,Th,Ph)
            klm += 0.31248919*SphHarmAng(20,16,1,Th,Ph)
            klm += 0.49719956*SphHarmAng(20,20,1,Th,Ph)
        else:
            klm = 0.66299538*SphHarmAng(20,4,1,Th,Ph)
            klm -= 0.11295259*SphHarmAng(20,8,1,Th,Ph)
            klm -= 0.42738441*SphHarmAng(20,12,1,Th,Ph)
            klm -= 0.52810433*SphHarmAng(20,16,1,Th,Ph)
            klm += 0.29347435*SphHarmAng(20,20,1,Th,Ph)
    else:   #shouldn't happen
        return 0.0
    return klm

def Glnh(SHCoef,psi,gam,SamSym):
    'needs doc string'
#    from . import pytexture as ptx

    Fln = np.zeros(len(SHCoef))
    for i,term in enumerate(SHCoef):
        l,m,n = eval(term.strip('C'))
        pcrs,dum = ptx.pyplmpsi(l,m,1,psi)
        pcrs *= RSQPI
        if m == 0:
            pcrs /= SQ2
        if SamSym in ['mmm',]:
            Ksl = pcrs*cosd(m*gam)
        else:
            Ksl = pcrs*(cosd(m*gam)+sind(m*gam))
        Fln[i] = SHCoef[term]*Ksl*Lnorm(l)
    ODFln = dict(zip(SHCoef.keys(),list(zip(SHCoef.values(),Fln))))
    return ODFln

def Flnh(SHCoef,phi,beta,SGData):
    'needs doc string'
#    from . import pytexture as ptx

    Fln = np.zeros(len(SHCoef))
    for i,term in enumerate(SHCoef):
        l,m,n = eval(term.strip('C'))
        if SGData['SGLaue'] in ['m3','m3m']:
            Kcl = 0.0
            for j in range(0,l+1,4):
                im = j//4
                pcrs,dum = ptx.pyplmpsi(l,j,1,phi)
                Kcl += BOH['L='+str(l)][n-1][im]*pcrs*cosd(j*beta)
        else:                #all but cubic
            pcrs,dum = ptx.pyplmpsi(l,n,1,phi)
            pcrs *= RSQPI
            if n == 0:
                pcrs /= SQ2
            if SGData['SGLaue'] in ['mmm','4/mmm','6/mmm','R3mR','3m1','31m']:
               if SGData['SGLaue'] in ['3mR','3m1','31m']:
                   if n%6 == 3:
                       Kcl = pcrs*sind(n*beta)
                   else:
                       Kcl = pcrs*cosd(n*beta)
               else:
                   Kcl = pcrs*cosd(n*beta)
            else:
                Kcl = pcrs*(cosd(n*beta)+sind(n*beta))
        Fln[i] = SHCoef[term]*Kcl*Lnorm(l)
    ODFln = dict(zip(SHCoef.keys(),list(zip(SHCoef.values(),Fln))))
    return ODFln

def polfcal(ODFln,SamSym,psi,gam):
    '''Perform a pole figure computation.
    Note that the the number of gam values must either be 1 or must
    match psi. Updated for numpy 1.8.0
    '''
#    from . import pytexture as ptx
    PolVal = np.ones_like(psi)
    for term in ODFln:
        if abs(ODFln[term][1]) > 1.e-3:
            l,m,n = eval(term.strip('C'))
            psrs,dum = ptx.pyplmpsi(l,m,len(psi),psi)
            if SamSym in ['-1','2/m']:
                if m:
                    Ksl = RSQPI*psrs*(cosd(m*gam)+sind(m*gam))
                else:
                    Ksl = RSQPI*psrs/SQ2
            else:
                if m:
                    Ksl = RSQPI*psrs*cosd(m*gam)
                else:
                    Ksl = RSQPI*psrs/SQ2
            PolVal += ODFln[term][1]*Ksl
    return PolVal

def invpolfcal(ODFln,SGData,phi,beta):
    'needs doc string'
#    from . import pytexture as ptx

    invPolVal = np.ones_like(beta)
    for term in ODFln:
        if abs(ODFln[term][1]) > 1.e-3:
            l,m,n = eval(term.strip('C'))
            if SGData['SGLaue'] in ['m3','m3m']:
                Kcl = 0.0
                for j in range(0,l+1,4):
                    im = j//4
                    pcrs,dum = ptx.pyplmpsi(l,j,len(beta),phi)
                    Kcl += BOH['L=%d'%(l)][n-1][im]*pcrs*cosd(j*beta)
            else:                #all but cubic
                pcrs,dum = ptx.pyplmpsi(l,n,len(beta),phi)
                pcrs *= RSQPI
                if n == 0:
                    pcrs /= SQ2
                if SGData['SGLaue'] in ['mmm','4/mmm','6/mmm','R3mR','3m1','31m']:
                   if SGData['SGLaue'] in ['3mR','3m1','31m']:
                       if n%6 == 3:
                           Kcl = pcrs*sind(n*beta)
                       else:
                           Kcl = pcrs*cosd(n*beta)
                   else:
                       Kcl = pcrs*cosd(n*beta)
                else:
                    Kcl = pcrs*(cosd(n*beta)+sind(n*beta))
            invPolVal += ODFln[term][1]*Kcl
    return invPolVal


def textureIndex(SHCoef):
    'needs doc string'
    Tindx = 1.0
    for term in SHCoef:
        l = eval(term.strip('C'))[0]
        Tindx += SHCoef[term]**2/(2.0*l+1.)
    return Tindx

UniqueCellByLaue = [
        [['m3','m3m'],(0,)],
        [['3R','3mR'],(0,3)],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],(0,2)],
        [['mmm'],(0,1,2)],
        [['2/m'+'a'],(0,1,2,3)],
        [['2/m'+'b'],(0,1,2,4)],
        [['2/m'+'c'],(0,1,2,5)],
        [['-1'],(0,1,2,3,4,5)],
    ]
'''List the unique cell terms by index for each Laue class'''

cellAlbl = ('a','b','c', 'alpha', 'beta', 'gamma')
'ASCII labels for a, b, c, alpha, beta, gamma'

cellUlbl = ('a','b','c',u'\u03B1',u'\u03B2',u'\u03B3')
'unicode labels for a, b, c, alpha, beta, gamma'
