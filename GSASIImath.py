# -*- coding: utf-8 -*-
#GSASIImath - major mathematics routines
########### SVN repository information ###################
# $Date: 2012-01-13 11:48:53 -0600 (Fri, 13 Jan 2012) $
# $Author: vondreele & toby $
# $Revision: 451 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIImath.py $
# $Id: GSASIImath.py 451 2012-01-13 17:48:53Z vondreele $
########### SVN repository information ###################
import sys
import os
import os.path as ospath
import random as rn
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import cPickle
import time
import math
import copy
import GSASIIpath
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import scipy.optimize as so
import scipy.linalg as sl
import numpy.fft as fft

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

def HessianLSQ(func,x0,Hess,args=(),ftol=1.49012e-8,xtol=1.49012e-8, maxcyc=0):
    
    """
    Minimize the sum of squares of a set of equations.

    ::
    
                    Nobs
        x = arg min(sum(func(y)**2,axis=0))
                    y=0

    Parameters
    ----------
    func : callable
        should take at least one (possibly length N vector) argument and
        returns M floating point numbers.
    x0 : ndarray
        The starting estimate for the minimization of length N
    Hess : callable
        A required function or method to compute the weighted vector and Hessian for func.
        It must be a symmetric NxN array
    args : tuple
        Any extra arguments to func are placed in this tuple.
    ftol : float 
        Relative error desired in the sum of squares.
    xtol : float
        Relative error desired in the approximate solution.
    maxcyc : int
        The maximum number of cycles of refinement to execute, if -1 refine 
        until other limits are met (ftol, xtol)

    Returns
    -------
    x : ndarray
        The solution (or the result of the last iteration for an unsuccessful
        call).
    cov_x : ndarray
        Uses the fjac and ipvt optional outputs to construct an
        estimate of the jacobian around the solution.  ``None`` if a
        singular matrix encountered (indicates very flat curvature in
        some direction).  This matrix must be multiplied by the
        residual standard deviation to get the covariance of the
        parameter estimates -- see curve_fit.
    infodict : dict
        a dictionary of optional outputs with the key s::

            - 'fvec' : the function evaluated at the output


    Notes
    -----

    """
                
    x0 = np.array(x0, ndmin=1)      #might be redundant?
    n = len(x0)
    if type(args) != type(()):
        args = (args,)
        
    icycle = 0
    One = np.ones((n,n))
    lam = 0.001
    lamMax = lam
    nfev = 0
    while icycle < maxcyc:
        lamMax = max(lamMax,lam)
        M = func(x0,*args)
        nfev += 1
        chisq0 = np.sum(M**2)
        Yvec,Amat = Hess(x0,*args)
        Adiag = np.sqrt(np.diag(Amat))
        psing = np.where(np.abs(Adiag) < 1.e-14,True,False)
        if np.any(psing):                #hard singularity in matrix
            return [x0,None,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':psing}]
        Anorm = np.outer(Adiag,Adiag)
        Yvec /= Adiag
        Amat /= Anorm        
        while True:
            Lam = np.eye(Amat.shape[0])*lam
            Amatlam = Amat*(One+Lam)
            try:
                Xvec = nl.solve(Amatlam,Yvec)
            except nl.LinAlgError:
                print 'ouch #1'
                psing = list(np.where(np.diag(nl.qr(Amatlam)[1]) < 1.e-14)[0])
                return [x0,None,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':psing}]
            Xvec /= Adiag
            M2 = func(x0+Xvec,*args)
            nfev += 1
            chisq1 = np.sum(M2**2)
            if chisq1 > chisq0:
                lam *= 10.
            else:
                x0 += Xvec
                lam /= 10.
                break
        if (chisq0-chisq1)/chisq0 < ftol:
            break
        icycle += 1
    M = func(x0,*args)
    nfev += 1
    Yvec,Amat = Hess(x0,*args)
    try:
        Bmat = nl.inv(Amat)
        return [x0,Bmat,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':[]}]
    except nl.LinAlgError:
        print 'ouch #2'
        psing = []
        if maxcyc:
            psing = list(np.where(np.diag(nl.qr(Amat)[1]) < 1.e-14)[0])
        return [x0,None,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':psing}]
   
def getVCov(varyNames,varyList,covMatrix):
    vcov = np.zeros((len(varyNames),len(varyNames)))
    for i1,name1 in enumerate(varyNames):
        for i2,name2 in enumerate(varyNames):
            try:
                vcov[i1][i2] = covMatrix[varyList.index(name1)][varyList.index(name2)]
            except ValueError:
                vcov[i1][i2] = 0.0
    return vcov
    
def getDistDerv(Oxyz,Txyz,Amat,Tunit,Top,SGData):
    
    def calcDist(Ox,Tx,U,inv,C,M,T,Amat):
        TxT = inv*(np.inner(M,Tx)+T)+C+U
        return np.sqrt(np.sum(np.inner(Amat,(TxT-Ox))**2))
        
    inv = Top/abs(Top)
    cent = abs(Top)/100
    op = abs(Top)%100-1
    M,T = SGData['SGOps'][op]
    C = SGData['SGCen'][cent]
    dx = .00001
    deriv = np.zeros(6)
    for i in [0,1,2]:
        Oxyz[i] += dx
        d0 = calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)
        Oxyz[i] -= 2*dx
        deriv[i] = (calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)-d0)/(2.*dx)
        Oxyz[i] += dx
        Txyz[i] += dx
        d0 = calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)
        Txyz[i] -= 2*dx
        deriv[i+3] = (calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)-d0)/(2.*dx)
        Txyz[i] += dx
    return deriv
    
def getAngSig(VA,VB,Amat,SGData,covData={}):
    
    def calcVec(Ox,Tx,U,inv,C,M,T,Amat):
        TxT = inv*(np.inner(M,Tx)+T)+C
        TxT = G2spc.MoveToUnitCell(TxT)+U
        return np.inner(Amat,(TxT-Ox))
        
    def calcAngle(Ox,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat):
        VecA = calcVec(Ox,TxA,unitA,invA,CA,MA,TA,Amat)
        VecA /= np.sqrt(np.sum(VecA**2))
        VecB = calcVec(Ox,TxB,unitB,invB,CB,MB,TB,Amat)
        VecB /= np.sqrt(np.sum(VecB**2))
        edge = VecB-VecA
        edge = np.sum(edge**2)
        angle = (2.-edge)/2.
        angle = max(angle,-1.)
        return acosd(angle)
        
    OxAN,OxA,TxAN,TxA,unitA,TopA = VA
    OxBN,OxB,TxBN,TxB,unitB,TopB = VB
    invA = invB = 1
    invA = TopA/abs(TopA)
    invB = TopB/abs(TopB)
    centA = abs(TopA)/100
    centB = abs(TopB)/100
    opA = abs(TopA)%100-1
    opB = abs(TopB)%100-1
    MA,TA = SGData['SGOps'][opA]
    MB,TB = SGData['SGOps'][opB]
    CA = SGData['SGCen'][centA]
    CB = SGData['SGCen'][centB]
    if 'covMatrix' in covData:
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        AngVcov = getVCov(OxAN+TxAN+TxBN,varyList,covMatrix)
        dx = .00001
        dadx = np.zeros(9)
        Ang = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
        for i in [0,1,2]:
            OxA[i] += dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            OxA[i] -= 2*dx
            dadx[i] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/dx
            OxA[i] += dx
            
            TxA[i] += dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            TxA[i] -= 2*dx
            dadx[i+3] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/dx
            TxA[i] += dx
            
            TxB[i] += dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            TxB[i] -= 2*dx
            dadx[i+6] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/dx
            TxB[i] += dx
            
        sigAng = np.sqrt(np.inner(dadx,np.inner(AngVcov,dadx)))
        if sigAng < 0.01:
            sigAng = 0.0
        return Ang,sigAng
    else:
        return calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat),0.0

def GetDistSig(Oatoms,Atoms,Amat,SGData,covData={}):

    def calcDist(Atoms,SyOps,Amat):
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        return np.sqrt(np.sum(V1**2))
        
    Inv = []
    SyOps = []
    names = []
    for i,atom in enumerate(Oatoms):
        names += atom[-1]
        Op,unit = Atoms[i][-1]
        inv = Op/abs(Op)
        m,t = SGData['SGOps'][abs(Op)%100-1]
        c = SGData['SGCen'][abs(Op)/100]
        SyOps.append([inv,m,t,c,unit])
    Dist = calcDist(Oatoms,SyOps,Amat)
    
    sig = -0.001
    if 'covMatrix' in covData:
        parmNames = []
        dx = .00001
        dadx = np.zeros(6)
        for i in range(6):
            ia = i/3
            ix = i%3
            Oatoms[ia][ix+1] += dx
            a0 = calcDist(Oatoms,SyOps,Amat)
            Oatoms[ia][ix+1] -= 2*dx
            dadx[i] = (calcDist(Oatoms,SyOps,Amat)-a0)/(2.*dx)
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        DistVcov = getVCov(names,varyList,covMatrix)
        sig = np.sqrt(np.inner(dadx,np.inner(DistVcov,dadx)))
        if sig < 0.001:
            sig = -0.001
    
    return Dist,sig

def GetAngleSig(Oatoms,Atoms,Amat,SGData,covData={}):
        
    def calcAngle(Atoms,SyOps,Amat):
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        V1 /= np.sqrt(np.sum(V1**2))
        V2 = XYZ[1]-XYZ[2]
        V2 /= np.sqrt(np.sum(V2**2))
        V3 = V2-V1
        cang = min(1.,max((2.-np.sum(V3**2))/2.,-1.))
        return acosd(cang)

    Inv = []
    SyOps = []
    names = []
    for i,atom in enumerate(Oatoms):
        names += atom[-1]
        Op,unit = Atoms[i][-1]
        inv = Op/abs(Op)
        m,t = SGData['SGOps'][abs(Op)%100-1]
        c = SGData['SGCen'][abs(Op)/100]
        SyOps.append([inv,m,t,c,unit])
    Angle = calcAngle(Oatoms,SyOps,Amat)
    
    sig = -0.01
    if 'covMatrix' in covData:
        parmNames = []
        dx = .00001
        dadx = np.zeros(9)
        for i in range(9):
            ia = i/3
            ix = i%3
            Oatoms[ia][ix+1] += dx
            a0 = calcAngle(Oatoms,SyOps,Amat)
            Oatoms[ia][ix+1] -= 2*dx
            dadx[i] = (calcAngle(Oatoms,SyOps,Amat)-a0)/(2.*dx)
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        AngVcov = getVCov(names,varyList,covMatrix)
        sig = np.sqrt(np.inner(dadx,np.inner(AngVcov,dadx)))
        if sig < 0.01:
            sig = -0.01
    
    return Angle,sig

def GetTorsionSig(Oatoms,Atoms,Amat,SGData,covData={}):
    
    def calcTorsion(Atoms,SyOps,Amat):
        
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        V2 = XYZ[2]-XYZ[1]
        V3 = XYZ[3]-XYZ[2]
        V1 /= np.sqrt(np.sum(V1**2))
        V2 /= np.sqrt(np.sum(V2**2))
        V3 /= np.sqrt(np.sum(V3**2))
        M = np.array([V1,V2,V3])
        D = nl.det(M)
        Ang = 1.0
        P12 = np.dot(V1,V2)
        P13 = np.dot(V1,V3)
        P23 = np.dot(V2,V3)
        Tors = acosd((P12*P23-P13)/(np.sqrt(1.-P12**2)*np.sqrt(1.-P23**2)))*D/abs(D)
        return Tors
            
    Inv = []
    SyOps = []
    names = []
    for i,atom in enumerate(Oatoms):
        names += atom[-1]
        Op,unit = Atoms[i][-1]
        inv = Op/abs(Op)
        m,t = SGData['SGOps'][abs(Op)%100-1]
        c = SGData['SGCen'][abs(Op)/100]
        SyOps.append([inv,m,t,c,unit])
    Tors = calcTorsion(Oatoms,SyOps,Amat)
    
    sig = -0.01
    if 'covMatrix' in covData:
        parmNames = []
        dx = .00001
        dadx = np.zeros(12)
        for i in range(12):
            ia = i/3
            ix = i%3
            Oatoms[ia][ix+1] += dx
            a0 = calcTorsion(Oatoms,SyOps,Amat)
            Oatoms[ia][ix+1] -= 2*dx
            dadx[i] = (calcTorsion(Oatoms,SyOps,Amat)-a0)/(2.*dx)
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        TorVcov = getVCov(names,varyList,covMatrix)
        sig = np.sqrt(np.inner(dadx,np.inner(TorVcov,dadx)))
        if sig < 0.01:
            sig = -0.01
    
    return Tors,sig
        
def GetDATSig(Oatoms,Atoms,Amat,SGData,covData={}):
    
    def calcDist(Atoms,SyOps,Amat):
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        return np.sqrt(np.sum(V1**2))
        
    def calcAngle(Atoms,SyOps,Amat):
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        V1 /= np.sqrt(np.sum(V1**2))
        V2 = XYZ[1]-XYZ[2]
        V2 /= np.sqrt(np.sum(V2**2))
        V3 = V2-V1
        cang = min(1.,max((2.-np.sum(V3**2))/2.,-1.))
        return acosd(cang)

    def calcTorsion(Atoms,SyOps,Amat):
        
        XYZ = []
        for i,atom in enumerate(Atoms):
            Inv,M,T,C,U = SyOps[i]
            XYZ.append(np.array(atom[1:4]))
            XYZ[-1] = Inv*(np.inner(M,np.array(XYZ[-1]))+T)+C+U
            XYZ[-1] = np.inner(Amat,XYZ[-1]).T
        V1 = XYZ[1]-XYZ[0]
        V2 = XYZ[2]-XYZ[1]
        V3 = XYZ[3]-XYZ[2]
        V1 /= np.sqrt(np.sum(V1**2))
        V2 /= np.sqrt(np.sum(V2**2))
        V3 /= np.sqrt(np.sum(V3**2))
        M = np.array([V1,V2,V3])
        D = nl.det(M)
        Ang = 1.0
        P12 = np.dot(V1,V2)
        P13 = np.dot(V1,V3)
        P23 = np.dot(V2,V3)
        Tors = acosd((P12*P23-P13)/(np.sqrt(1.-P12**2)*np.sqrt(1.-P23**2)))*D/abs(D)
        return Tors
            
    Inv = []
    SyOps = []
    names = []
    for i,atom in enumerate(Oatoms):
        names += atom[-1]
        Op,unit = Atoms[i][-1]
        inv = Op/abs(Op)
        m,t = SGData['SGOps'][abs(Op)%100-1]
        c = SGData['SGCen'][abs(Op)/100]
        SyOps.append([inv,m,t,c,unit])
    M = len(Oatoms)
    if M == 2:
        Val = calcDist(Oatoms,SyOps,Amat)
    elif M == 3:
        Val = calcAngle(Oatoms,SyOps,Amat)
    else:
        Val = calcTorsion(Oatoms,SyOps,Amat)
    
    sigVals = [-0.001,-0.01,-0.01]
    sig = sigVals[M-3]
    if 'covMatrix' in covData:
        parmNames = []
        dx = .00001
        N = M*3
        dadx = np.zeros(N)
        for i in range(N):
            ia = i/3
            ix = i%3
            Oatoms[ia][ix+1] += dx
            if M == 2:
                a0 = calcDist(Oatoms,SyOps,Amat)
            elif M == 3:
                a0 = calcAngle(Oatoms,SyOps,Amat)
            else:
                a0 = calcTorsion(Oatoms,SyOps,Amat)
            Oatoms[ia][ix+1] -= 2*dx
            if M == 2:
                dadx[i] = (calcDist(Oatoms,SyOps,Amat)-a0)/(2.*dx)                
            elif M == 3:
                dadx[i] = (calcAngle(Oatoms,SyOps,Amat)-a0)/(2.*dx)                
            else:
                dadx[i] = (calcTorsion(Oatoms,SyOps,Amat)-a0)/(2.*dx)
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        Vcov = getVCov(names,varyList,covMatrix)
        sig = np.sqrt(np.inner(dadx,np.inner(Vcov,dadx)))
        if sig < sigVals[M-3]:
            sig = sigVals[M-3]
    
    return Val,sig
        
    
def ValEsd(value,esd=0,nTZ=False):                  #NOT complete - don't use
    # returns value(esd) string; nTZ=True for no trailing zeros
    # use esd < 0 for level of precision shown e.g. esd=-0.01 gives 2 places beyond decimal
    #get the 2 significant digits in the esd 
    edig = lambda esd: int(round(10**(math.log10(esd) % 1+1)))
    #get the number of digits to represent them 
    epl = lambda esd: 2+int(1.545-math.log10(10*edig(esd)))
    
    mdec = lambda esd: -int(round(math.log10(abs(esd))))+1
    ndec = lambda esd: int(1.545-math.log10(abs(esd)))
    if esd > 0:
        fmt = '"%.'+str(ndec(esd))+'f(%d)"'
        return str(fmt%(value,int(round(esd*10**(mdec(esd)))))).strip('"')
    elif esd < 0:
         return str(round(value,mdec(esd)-1))
    else:
        text = str("%f"%(value))
        if nTZ:
            return text.rstrip('0')
        else:
            return text

def FourierMap(data,reflData):
    
    generalData = data['General']
    if not generalData['Map']['MapType']:
        print '**** ERROR - Fourier map not defined'
        return
    mapData = generalData['Map']
    dmin = mapData['Resolution']
    SGData = generalData['SGData']
    cell = generalData['Cell'][1:8]        
    A = G2lat.cell2A(cell[:6])
    Hmax = np.asarray(G2lat.getHKLmax(dmin,SGData,A),dtype='i')+1
    Fhkl = np.zeros(shape=2*Hmax,dtype='c16')
#    Fhkl[0,0,0] = generalData['F000X']
    time0 = time.time()
    for ref in reflData:
        if ref[4] >= dmin:
            Fosq,Fcsq,ph = ref[8:11]
            for i,hkl in enumerate(ref[11]):
                hkl = np.asarray(hkl,dtype='i')
                dp = 360.*ref[12][i]
                a = cosd(ph+dp)
                b = sind(ph+dp)
                phasep = complex(a,b)
                phasem = complex(a,-b)
                if 'Fobs' in mapData['MapType']:
                    F = np.sqrt(Fosq)
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = F*phasep
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = F*phasem
                elif 'Fcalc' in mapData['MapType']:
                    F = np.sqrt(Fcsq)
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = F*phasep
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = F*phasem
                elif 'delt-F' in mapData['MapType']:
                    dF = np.sqrt(Fosq)-np.sqrt(Fcsq)
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = dF*phasep
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = dF*phasem
                elif '2*Fo-Fc' in mapData['MapType']:
                    F = 2.*np.sqrt(Fosq)-np.sqrt(Fcsq)
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = F*phasep
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = F*phasem
                elif 'Patterson' in mapData['MapType']:
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = complex(Fosq,0.)
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = complex(Fosq,0.)
    rho = fft.fftn(fft.fftshift(Fhkl))/cell[6]
    print 'Fourier map time: %.4f'%(time.time()-time0),'no. elements: %d'%(Fhkl.size)
    mapData['rho'] = np.real(rho)
    mapData['rhoMax'] = max(np.max(mapData['rho']),-np.min(mapData['rho']))
    return mapData
    
# map printing for testing purposes
def printRho(SGLaue,rho,rhoMax):                          
    dim = len(rho.shape)
    if dim == 2:
        ix,jy = rho.shape
        for j in range(jy):
            line = ''
            if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                line += (jy-j)*'  '
            for i in range(ix):
                r = int(100*rho[i,j]/rhoMax)
                line += '%4d'%(r)
            print line+'\n'
    else:
        ix,jy,kz = rho.shape
        for k in range(kz):
            print 'k = ',k
            for j in range(jy):
                line = ''
                if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                    line += (jy-j)*'  '
                for i in range(ix):
                    r = int(100*rho[i,j,k]/rhoMax)
                    line += '%4d'%(r)
                print line+'\n'
## keep this
                
def findOffset(SGData,rho,Fhkl):
    
    def calcPhase(DX,DH,Dphi):
        H,K,L = DH.T
        Phi = (H*DX[0]+K*DX[1]+L*DX[2]+0.5) % 1.
        return Dphi-Phi
    
    if SGData['SpGrp'] == 'P 1':
        return [0,0,0]    
# will need to consider 'SGPolax': one of '','x','y','x y','z','x z','y z','xyz','111'
    mapShape = rho.shape
    steps = np.array(mapShape)
    hklShape = Fhkl.shape
    mapHalf = np.array(mapShape)/2
    Fmax = np.max(np.absolute(Fhkl))
    hklHalf = np.array(hklShape)/2
    sortHKL = np.argsort(Fhkl.flatten())
    Fdict = {}
    for hkl in sortHKL:
        HKL = np.unravel_index(hkl,hklShape)
        F = Fhkl[HKL[0]][HKL[1]][HKL[2]]
        if F == 0.:
            break
        Fdict['%.6f'%(np.absolute(F))] = hkl
    Flist = np.flipud(np.sort(Fdict.keys()))
    F = str(1.e6)
    i = 0
    DH = []
    Dphi = []
    while i < 20:
        F = Flist[i]
        hkl = np.unravel_index(Fdict[F],hklShape)
        iabsnt,mulp,Uniq,Phi = G2spc.GenHKLf(list(hkl-hklHalf),SGData)
        Uniq = np.array(Uniq,dtype='i')
        Phi = np.array(Phi)
        Uniq = np.concatenate((Uniq,-Uniq))+hklHalf         # put in Friedel pairs & make as index to Farray
        Phi = np.concatenate((Phi,-Phi))                      # and their phase shifts
        Fh0 = Fhkl[hkl[0],hkl[1],hkl[2]]
        ang0 = np.angle(Fh0,deg=True)/360.
        for j,H in enumerate(Uniq[1:]):
            ang = (np.angle(Fhkl[H[0],H[1],H[2]],deg=True)/360.-Phi[j+1])
            dH = H-hkl
            dang = ang-ang0
            if np.any(np.abs(dH)-8 > 0):    #keep low order DHs
                continue
            DH.append(dH)
            Dphi.append((dang+0.5) % 1.0)
        i += 1
    DH = np.array(DH)
    Dphi = np.array(Dphi)
#    for item in zip(DH,Dphi):
#        print item[0],'%.4f'%(item[1])
    DX = np.zeros(3)
    X,Y,Z = np.mgrid[0:1:10j,0:1:10j,0:1:10j]
    XYZ = np.array(zip(X.flatten(),Y.flatten(),Z.flatten()))
    Mmin = 1.e10
    
    for xyz in XYZ:             #do a global search for best roll
        M = np.sum(calcPhase(xyz,DH,Dphi)**2)
        if M < Mmin:
            DX = xyz
            Mmin = M
    
    result = so.leastsq(calcPhase,DX,full_output=True,args=(DH,Dphi))
    for item in zip(DH,Dphi,result[2]['fvec']):
        print item[0],'%.4f %.4f'%(item[1],item[2])
    chisq = np.sum(result[2]['fvec']**2)
    DX = np.array(np.fix(-result[0]*steps),dtype='i')
    print ' map offset chi**2: %.3f, map offset: %d %d %d'%(chisq,DX[0],DX[1],DX[2])
    return DX
    
def findOffset2(SGData,rho,Fhkl):
    if SGData['SpGrp'] == 'P 1':
        return [0,0,0]
    mapShape = rho.shape
    mapHalf = np.array(mapShape)/2
    steps = np.array(mapShape)
    hklShape = Fhkl.shape
    hklHalf = np.array(hklShape)/2
    for M,T in SGData['SGOps'][1:]:
        FThkl = np.zeros(shape=hklShape,dtype='c16')
        for ind,F in enumerate(Fhkl.flatten()):
            if F:
                HKL = np.unravel_index(ind,hklShape)-hklHalf
                HKLT = np.inner(HKL,M.T)
                phi = (np.angle(F,deg=True)/360.+np.vdot(T,HKL) % 1.)*360.
                phasep = complex(cosd(phi),-sind(phi))      #want F*
                h,k,l = HKLT+hklHalf
                FThkl[h,k,l] = np.absolute(F)*phasep
        Cd = np.real(fft.fftn(fft.fftshift(Fhkl*FThkl)))
        CdMax = np.array(np.unravel_index(np.argmax(Cd),mapShape))
        print CdMax,CdMax/2
    DX = CdMax/2            
    
    print ' Test map offset : %d %d %d'%(DX[0],DX[1],DX[2])
    return DX


def ChargeFlip(data,reflData,pgbar):
    generalData = data['General']
    mapData = generalData['Map']
    flipData = generalData['Flip']
    FFtable = {}
    if 'None' not in flipData['Norm element']:
        normElem = flipData['Norm element'].upper()
        FFs = G2el.GetFormFactorCoeff(normElem.split('+')[0].split('-')[0])
        for ff in FFs:
            if ff['Symbol'] == normElem:
                FFtable.update(ff)
    dmin = flipData['Resolution']
    SGData = generalData['SGData']
    cell = generalData['Cell'][1:8]        
    A = G2lat.cell2A(cell[:6])
    Vol = cell[6]
    Hmax = np.asarray(G2lat.getHKLmax(dmin,SGData,A),dtype='i')+1
    Ehkl = np.zeros(shape=2*Hmax,dtype='c16')       #2X64bits per complex no.
    time0 = time.time()
    for ref in reflData:
        dsp = ref[4]
        if dsp >= dmin:
            ff = 0.1*Vol    #est. no. atoms for ~10A**3/atom
            if FFtable:
                SQ = 0.25/dsp**2
                ff *= G2el.ScatFac(FFtable,SQ)[0]
            if ref[8] > 0.:
                E = np.sqrt(ref[8])/ff
            else:
                E = 0.
            ph = ref[10]
            ph = rn.uniform(0.,360.)
            for i,hkl in enumerate(ref[11]):
                hkl = np.asarray(hkl,dtype='i')
                dp = 360.*ref[12][i]
                a = cosd(ph+dp)
                b = sind(ph+dp)
                phasep = complex(a,b)
                phasem = complex(a,-b)
                h,k,l = hkl+Hmax
                Ehkl[h,k,l] = E*phasep
                h,k,l = -hkl+Hmax       #Friedel pair refl.
                Ehkl[h,k,l] = E*phasem
#    Ehkl[Hmax] = 0.00001           #this to preserve F[0,0,0]
    CEhkl = copy.copy(Ehkl)
    MEhkl = ma.array(Ehkl,mask=(Ehkl==0.0))
    Emask = ma.getmask(MEhkl)
    sumE = np.sum(ma.array(np.absolute(CEhkl),mask=Emask))
    Ncyc = 0
    old = np.seterr(all='raise')
    while True:        
        try:
            CErho = np.real(fft.fftn(fft.fftshift(CEhkl)))*(1.+0j)
            CEsig = np.std(CErho)
            CFrho = np.where(np.real(CErho) >= flipData['k-factor']*CEsig,CErho,-CErho)
            CFhkl = fft.ifftshift(fft.ifftn(CFrho))
            phase = CFhkl/np.absolute(CFhkl)
            CEhkl = np.absolute(Ehkl)*phase
            Ncyc += 1
            sumCF = np.sum(ma.array(np.absolute(CFhkl),mask=Emask))
            DEhkl = np.absolute(np.absolute(Ehkl)/sumE-np.absolute(CFhkl)/sumCF)
            Rcf = min(100.,np.sum(ma.array(DEhkl,mask=Emask)*100.))
            if Rcf < 5.:
                break
            GoOn = pgbar.Update(Rcf,newmsg='%s%8.3f%s\n%s %d'%('Residual Rcf =',Rcf,'%','No.cycles = ',Ncyc))[0]
            if not GoOn:
                break
        except FloatingPointError:
            Rcf = 100.
            break
    np.seterr(**old)
    print ' Charge flip time: %.4f'%(time.time()-time0),'no. elements: %d'%(Ehkl.size)
    print ' No.cycles = ',Ncyc,'Residual Rcf =%8.3f%s'%(Rcf,'%')+' Map size:',CErho.shape
    CErho = np.real(fft.fftn(fft.fftshift(CEhkl)))
#    roll = findOffset2(SGData,CErho,CEhkl) #doesn't work
    roll = findOffset(SGData,CErho,CEhkl)
        
    mapData['Rcf'] = Rcf
    mapData['rho'] = np.roll(np.roll(np.roll(CErho,roll[0],axis=0),roll[1],axis=1),roll[2],axis=2)
    mapData['rhoMax'] = max(np.max(mapData['rho']),-np.min(mapData['rho']))
    mapData['rollMap'] = [0,0,0]
    return mapData
    
def SearchMap(data,keepDup=False):
    
    norm = 1./(np.sqrt(3.)*np.sqrt(2.*np.pi)**3)
    
    def noDuplicate(xyz,peaks,SGData,incre):                  #be careful where this is used - it's slow
        equivs = G2spc.GenAtom(xyz,SGData)
        xyzs = [equiv[0] for equiv in equivs]
        for x in xyzs:
            if True in [np.allclose(x,peak,atol=0.02) for peak in peaks]:
                return False
        return True
            
    def findRoll(rhoMask,mapHalf):
        indx = np.array(ma.nonzero(rhoMask)).T
        rhoList = np.array([rho[i,j,k] for i,j,k in indx])
        rhoMax = np.max(rhoList)
        return mapHalf-indx[np.argmax(rhoList)]
        
    def rhoCalc(parms,rX,rY,rZ,res,SGLaue):
        Mag,x0,y0,z0,sig = parms
        if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
            return norm*Mag*np.exp(-((x0-rX)**2+(y0-rY)**2+(x0-rX)*(y0-rY)+(z0-rZ)**2)/(2.*sig**2))/(sig*res**3)
        else:
            return norm*Mag*np.exp(-((x0-rX)**2+(y0-rY)**2+(z0-rZ)**2)/(2.*sig**2))/(sig*res**3)
        
    def peakFunc(parms,rX,rY,rZ,rho,res,SGLaue):
        Mag,x0,y0,z0,sig = parms
        M = rho-rhoCalc(parms,rX,rY,rZ,res,SGLaue)
        return M
        
    def peakHess(parms,rX,rY,rZ,rho,res,SGLaue):
        Mag,x0,y0,z0,sig = parms
        dMdv = np.zeros(([5,]+list(rX.shape)))
        delt = .01
        for i in range(5):
            parms[i] -= delt
            rhoCm = rhoCalc(parms,rX,rY,rZ,res,SGLaue)
            parms[i] += 2.*delt
            rhoCp = rhoCalc(parms,rX,rY,rZ,res,SGLaue)
            parms[i] -= delt
            dMdv[i] = (rhoCp-rhoCm)/(2.*delt)
        rhoC = rhoCalc(parms,rX,rY,rZ,res,SGLaue)
        Vec = np.sum(np.sum(np.sum(dMdv*(rho-rhoC),axis=3),axis=2),axis=1)
        dMdv = np.reshape(dMdv,(5,rX.size))
        Hess = np.inner(dMdv,dMdv)
        
        return Vec,Hess
        
    generalData = data['General']
    phaseName = generalData['Name']
    SGData = generalData['SGData']
    cell = generalData['Cell'][1:8]        
    A = G2lat.cell2A(cell[:6])
    drawingData = data['Drawing']
    peaks = []
    mags = []
    try:
        mapData = generalData['Map']
        contLevel = mapData['cutOff']*mapData['rhoMax']/100.
        rho = copy.copy(mapData['rho'])     #don't mess up original
        mapHalf = np.array(rho.shape)/2
        res = mapData['Resolution']
        incre = 1./np.array(rho.shape)
        step = max(1.0,1./res)+1
        steps = np.array(3*[step,])
    except KeyError:
        print '**** ERROR - Fourier map not defined'
        return peaks,mags
    while True:
        rhoMask = ma.array(rho,mask=(rho<contLevel))
        if not ma.count(rhoMask):
            break
        rMI = findRoll(rhoMask,mapHalf)
        rho = np.roll(np.roll(np.roll(rho,rMI[0],axis=0),rMI[1],axis=1),rMI[2],axis=2)
        rMM = mapHalf-steps
        rMP = mapHalf+steps+1
        rhoPeak = rho[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]]
        peakInt = np.sum(rhoPeak)*res**3
        rX,rY,rZ = np.mgrid[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]]
        x0 = [peakInt,mapHalf[0],mapHalf[1],mapHalf[2],2.0]          #magnitude, position & width(sig)
        result = HessianLSQ(peakFunc,x0,Hess=peakHess,
            args=(rX,rY,rZ,rhoPeak,res,SGData['SGLaue']),ftol=.01,maxcyc=10)
        x1 = result[0]
        if np.any(x1 < 0):
            break
        peak = (np.array(x1[1:4])-rMI)*incre
        if not len(peaks):
            peaks.append(peak)
            mags.append(x1[0])
        else:
            if keepDup:
                peaks.append(peak)
                mags.append(x1[0])
            elif noDuplicate(peak,peaks,SGData,incre) and x1[0] > 0.:
                peaks.append(peak)
                mags.append(x1[0])
            if len(peaks) > 300:
                break
        rho[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]] = peakFunc(x1,rX,rY,rZ,rhoPeak,res,SGData['SGLaue'])
        rho = np.roll(np.roll(np.roll(rho,-rMI[2],axis=2),-rMI[1],axis=1),-rMI[0],axis=0)
    return np.array(peaks),np.array([mags,]).T