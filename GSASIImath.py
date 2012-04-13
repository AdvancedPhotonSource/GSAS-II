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
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import cPickle
import time
import math
import copy
import GSASIIpath
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import scipy.optimize as so
import scipy.linalg as sl

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
        if 0.0 in Adiag:                #hard singularity in matrix
            psing = list(np.where(Adiag == 0.)[0])
            return [x0,None,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':psing}]
        Anorm = np.outer(Adiag,Adiag)
        Yvec /= Adiag
        Amat /= Anorm        
        while True:
            Lam = np.eye(Amat.shape[0])*lam
            Amatlam = Amat*(One+Lam)
            try:
                Xvec = nl.solve(Amatlam,Yvec)
            except LinAlgError:
                psing = list(np.where(np.diag(nl.gr(Amatlam)[1]) < 1.e-14)[0])
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
    
    import scipy.fftpack as fft
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
            for i,hkl in enumerate(ref[11]):
                hkl = np.asarray(hkl,dtype='i')
                Fosq,Fcsq,ph = ref[8:11]
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
    Fhkl = fft.fftshift(Fhkl)
    rho = fft.fftn(Fhkl)/cell[6]
    print 'Fourier map time: %.4f'%(time.time()-time0),'no. elements: %d'%(Fhkl.size)
    mapData['rho'] = np.real(rho)
    mapData['rhoMax'] = max(np.max(mapData['rho']),-np.min(mapData['rho']))
    return mapData
    
def SearchMap(data,keepDup=False):
    
    norm = 1./(np.sqrt(3.)*np.sqrt(2.*np.pi)**3)
    
    def noDuplicate(xyz,peaks,SGData):                  #be careful where this is used - it's slow
        equivs = G2spc.GenAtom(xyz,SGData)
        xyzs = [equiv[0] for equiv in equivs]
        for x in xyzs:
            if True in [np.allclose(x,peak,atol=0.002) for peak in peaks]:
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
            args=(rX,rY,rZ,rhoPeak,res,SGData['SGLaue']),ftol=.0001,maxcyc=100)
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
            elif noDuplicate(peak,peaks,SGData) and x1[0] > 0.:
                peaks.append(peak)
                mags.append(x1[0])
        rho[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]] = peakFunc(result[0],rX,rY,rZ,rhoPeak,res,SGData['SGLaue'])
        rho = np.roll(np.roll(np.roll(rho,-rMI[2],axis=2),-rMI[1],axis=1),-rMI[0],axis=0)
    return np.array(peaks),np.array([mags,]).T