# -*- coding: utf-8 -*-
#GSASIImath - major mathematics routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
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
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import numpy.fft as fft

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atand = lambda x: 180.*np.arctan(x)/np.pi
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
        print 'ouch #2 linear algebra error in LS'
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

def FindAtomIndexByIDs(atomData,IDs,Draw=True):
    indx = []
    for i,atom in enumerate(atomData):
        if Draw and atom[-3] in IDs:
            indx.append(i)
        elif atom[-1] in IDs:
            indx.append(i)
    return indx

def FillAtomLookUp(atomData):
    atomLookUp = {}
    for iatm,atom in enumerate(atomData):
        atomLookUp[atom[-1]] = iatm
    return atomLookUp

def GetAtomsById(atomData,atomLookUp,IdList):
    atoms = []
    for id in IdList:
        atoms.append(atomData[atomLookUp[id]])
    return atoms
    
def GetAtomItemsById(atomData,atomLookUp,IdList,itemLoc,numItems=1):
    Items = []
    if not isinstance(IdList,list):
        IdList = [IdList,]
    for id in IdList:
        if numItems == 1:
            Items.append(atomData[atomLookUp[id]][itemLoc])
        else:
            Items.append(atomData[atomLookUp[id]][itemLoc:itemLoc+numItems])
    return Items
    
def GetAtomCoordsByID(pId,parmDict,AtLookup,indx):
    pfx = [str(pId)+'::A'+i+':' for i in ['x','y','z']]
    dpfx = [str(pId)+'::dA'+i+':' for i in ['x','y','z']]
    XYZ = []
    for ind in indx:
        names = [pfx[i]+str(AtLookup[ind]) for i in range(3)]
        dnames = [dpfx[i]+str(AtLookup[ind]) for i in range(3)]
        XYZ.append([parmDict[name]+parmDict[dname] for name,dname in zip(names,dnames)])
    return XYZ

def AtomUij2TLS(atomData,atPtrs,Amat,Bmat,rbObj):   #unfinished & not used
    for atom in atomData:
        XYZ = np.inner(Amat,atom[cx:cx+3])
        if atom[cia] == 'A':
            UIJ = atom[cia+2:cia+8]
                
def TLS2Uij(xyz,g,Amat,rbObj):
    TLStype,TLS = rbObj['ThermalMotion'][:2]
    Tmat = np.zeros((3,3))
    Lmat = np.zeros((3,3))
    Smat = np.zeros((3,3))
    gvec = np.sqrt(np.array([g[0][0]**2,g[1][1]**2,g[2][2]**2,
        g[0][0]*g[1][1],g[0][0]*g[2][2],g[1][1]*g[2][2]]))
    if 'T' in TLStype:
        Tmat = G2lat.U6toUij(TLS[:6])
    if 'L' in TLStype:
        Lmat = G2lat.U6toUij(TLS[6:12])
    if 'S' in TLStype:
        Smat = np.array([[TLS[18],TLS[12],TLS[13]],[TLS[14],TLS[19],TLS[15]],[TLS[16],TLS[17],0] ])
    XYZ = np.inner(Amat,xyz)
    Axyz = np.array([[ 0,XYZ[2],-XYZ[1]], [-XYZ[2],0,XYZ[0]], [XYZ[1],-XYZ[0],0]] )
    Umat = Tmat+np.inner(Axyz,Smat)+np.inner(Smat.T,Axyz.T)+np.inner(np.inner(Axyz,Lmat),Axyz.T)
    beta = np.inner(np.inner(g,Umat),g)
    return G2lat.UijtoU6(beta)*gvec    
        
def AtomTLS2UIJ(atomData,atPtrs,Amat,rbObj):
    cx,ct,cs,cia = atPtrs
    TLStype,TLS = rbObj['ThermalMotion'][:2]
    Tmat = np.zeros((3,3))
    Lmat = np.zeros((3,3))
    Smat = np.zeros((3,3))
    G,g = G2lat.A2Gmat(Amat)
    gvec = 1./np.sqrt(np.array([g[0][0],g[1][1],g[2][2],g[0][1],g[0][2],g[1][2]]))
    if 'T' in TLStype:
        Tmat = G2lat.U6toUij(TLS[:6])
    if 'L' in TLStype:
        Lmat = G2lat.U6toUij(TLS[6:12])
    if 'S' in TLStype:
        Smat = np.array([ [TLS[18],TLS[12],TLS[13]], [TLS[14],TLS[19],TLS[15]], [TLS[16],TLS[17],0] ])
    for atom in atomData:
        XYZ = np.inner(Amat,atom[cx:cx+3])
        Axyz = np.array([ 0,XYZ[2],-XYZ[1], -XYZ[2],0,XYZ[0], XYZ[1],-XYZ[0],0],ndmin=2 )
        if 'U' in TSLtype:
            atom[cia+1] = TLS[0]
            atom[cia] = 'I'
        else:
            atom[cia] = 'A'
            Umat = Tmat+np.inner(Axyz,Smat)+np.inner(Smat.T,Axyz.T)+np.inner(np.inner(Axyz,Lmat),Axyz.T)
            beta = np.inner(np.inner(g,Umat),g)
            atom[cia+2:cia+8] = G2spc.U2Uij(beta/gvec)

def GetXYZDist(xyz,XYZ,Amat):
    ''' gets distance from position xyz to all XYZ, xyz & XYZ are np.array 
        and are in crystal coordinates; Amat is crystal to Cart matrix
    '''
    return np.sqrt(np.sum(np.inner(Amat,XYZ-xyz)**2,axis=0))

def getAtomXYZ(atoms,cx):
    XYZ = []
    for atom in atoms:
        XYZ.append(atom[cx:cx+3])
    return np.array(XYZ)

def UpdateRBXYZ(Bmat,RBObj,RBData,RBType):
    ''' Returns crystal coordinates for atoms described by RBObj
    '''
    RBRes = RBData[RBType][RBObj['RBId']]
    if RBType == 'Vector':
        vecs = RBRes['rbVect']
        mags = RBRes['VectMag']
        Cart = np.zeros_like(vecs[0])
        for vec,mag in zip(vecs,mags):
            Cart += vec*mag
    elif RBType == 'Residue':
        Cart = np.array(RBRes['rbXYZ'])
        for tor,seq in zip(RBObj['Torsions'],RBRes['rbSeq']):
            QuatA = AVdeg2Q(tor[0],Cart[seq[0]]-Cart[seq[1]])
            for ride in seq[3]:
                Cart[ride] = prodQVQ(QuatA,Cart[ride]-Cart[seq[1]])+Cart[seq[1]]
    XYZ = np.zeros_like(Cart)
    for i,xyz in enumerate(Cart):
        X = prodQVQ(RBObj['Orient'][0],xyz)
        XYZ[i] = np.inner(Bmat,X)+RBObj['Orig'][0]
    return XYZ,Cart
    
def UpdateRBUIJ(Bmat,Cart,RBObj):
    ''' Returns atom I/A, Uiso or UIJ for atoms at XYZ as described by RBObj
    '''
    TLStype,TLS = RBObj['ThermalMotion'][:2]
    T = np.zeros(6)
    L = np.zeros(6)
    S = np.zeros(8)
    if 'T' in TLStype:
        T = TLS[:6]
    if 'L' in TLStype:
        L = np.array(TLS[6:12])*(np.pi/180.)**2
    if 'S' in TLStype:
        S = np.array(TLS[12:])*(np.pi/180.)
    g = np.inner(Bmat,Bmat.T)
    gvec = 1./np.sqrt(np.array([g[0][0]**2,g[1][1]**2,g[2][2]**2,
        g[0][0]*g[1][1],g[0][0]*g[2][2],g[1][1]*g[2][2]]))
    Uout = []
    Q = RBObj['Orient'][0]
    for X in Cart:
        X = prodQVQ(Q,X)
        if 'U' in TLStype:
            Uout.append(['I',TLS[0],0,0,0,0,0,0])
        elif not 'N' in TLStype:
            U = [0,0,0,0,0,0]
            U[0] = T[0]+L[1]*X[2]**2+L[2]*X[1]**2-2.0*L[5]*X[1]*X[2]+2.0*(S[2]*X[2]-S[4]*X[1])
            U[1] = T[1]+L[0]*X[2]**2+L[2]*X[0]**2-2.0*L[4]*X[0]*X[2]+2.0*(S[5]*X[0]-S[0]*X[2])
            U[2] = T[2]+L[1]*X[0]**2+L[0]*X[1]**2-2.0*L[3]*X[1]*X[0]+2.0*(S[1]*X[1]-S[3]*X[0])
            U[3] = T[3]+L[4]*X[1]*X[2]+L[5]*X[0]*X[2]-L[3]*X[2]**2-L[2]*X[0]*X[1]+  \
                S[4]*X[0]-S[5]*X[1]-(S[6]+S[7])*X[2]
            U[4] = T[4]+L[3]*X[1]*X[2]+L[5]*X[0]*X[1]-L[4]*X[1]**2-L[1]*X[0]*X[2]+  \
                S[3]*X[2]-S[2]*X[0]+S[6]*X[1]
            U[5] = T[5]+L[3]*X[0]*X[2]+L[4]*X[0]*X[1]-L[5]*X[0]**2-L[0]*X[2]*X[1]+  \
                S[0]*X[1]-S[1]*X[2]+S[7]*X[0]
            Umat = G2lat.U6toUij(U)
            beta = np.inner(np.inner(Bmat,Umat),Bmat.T)
            Uout.append(['A',0.0,]+list(G2lat.UijtoU6(beta)*gvec))
        else:
            Uout.append(['N',])
    return Uout
    
def GetSHCoeff(pId,parmDict,SHkeys):
    SHCoeff = {}
    for shkey in SHkeys:
        shname = str(pId)+'::'+shkey
        SHCoeff[shkey] = parmDict[shname]
    return SHCoeff
        
def getMass(generalData):
    mass = 0.
    for i,elem in enumerate(generalData['AtomTypes']):
        mass += generalData['NoAtoms'][elem]*generalData['AtomMass'][i]
    return mass    

def getDensity(generalData):
    
    mass = getMass(generalData)
    Volume = generalData['Cell'][7]
    density = mass/(0.6022137*Volume)
    return density,Volume/mass
    
def getSyXYZ(XYZ,ops,SGData):
    XYZout = np.zeros_like(XYZ)
    for i,[xyz,op] in enumerate(zip(XYZ,ops)):
        if op == '1':
            XYZout[i] = xyz
        else:
            oprs = op.split('+')
            unit = [0,0,0]
            if oprs[1]:
                unit = np.array(list(eval(oprs[1])))
            syop =int(oprs[0])
            inv = syop/abs(syop)
            syop *= inv
            cent = syop/100
            syop %= 100
            syop -= 1
            M,T = SGData['SGOps'][syop]
            XYZout[i] = (np.inner(M,xyz)+T)*inv+SGData['SGCen'][cent]+unit
    return XYZout
    
def getRestDist(XYZ,Amat):
    return np.sqrt(np.sum(np.inner(Amat,(XYZ[1]-XYZ[0]))**2))
    
def getRestDeriv(Func,XYZ,Amat,ops,SGData):
    deriv = np.zeros((len(XYZ),3))
    dx = 0.00001
    for j,xyz in enumerate(XYZ):
        for i,x in enumerate(np.array([[dx,0,0],[0,dx,0],[0,0,dx]])):
            XYZ[j] -= x
            d1 = Func(getSyXYZ(XYZ,ops,SGData),Amat)
            XYZ[j] += 2*x
            d2 = Func(getSyXYZ(XYZ,ops,SGData),Amat)
            XYZ[j] -= x
            deriv[j][i] = (d1-d2)/(2*dx)
    return deriv.flatten()

def getRestAngle(XYZ,Amat):
    
    def calcVec(Ox,Tx,Amat):
        return np.inner(Amat,(Tx-Ox))

    VecA = calcVec(XYZ[1],XYZ[0],Amat)
    VecA /= np.sqrt(np.sum(VecA**2))
    VecB = calcVec(XYZ[1],XYZ[2],Amat)
    VecB /= np.sqrt(np.sum(VecB**2))
    edge = VecB-VecA
    edge = np.sum(edge**2)
    angle = (2.-edge)/2.
    angle = max(angle,-1.)
    return acosd(angle)
    
def getRestPlane(XYZ,Amat):
    sumXYZ = np.zeros(3)
    for xyz in XYZ:
        sumXYZ += xyz
    sumXYZ /= len(XYZ)
    XYZ = np.array(XYZ)-sumXYZ
    XYZ = np.inner(Amat,XYZ).T
    Zmat = np.zeros((3,3))
    for i,xyz in enumerate(XYZ):
        Zmat += np.outer(xyz.T,xyz)
    Evec,Emat = nl.eig(Zmat)
    Evec = np.sqrt(Evec)/(len(XYZ)-3)
    Order = np.argsort(Evec)
    return Evec[Order[0]]
    
def getRestChiral(XYZ,Amat):    
    VecA = np.empty((3,3))    
    VecA[0] = np.inner(XYZ[1]-XYZ[0],Amat)
    VecA[1] = np.inner(XYZ[2]-XYZ[0],Amat)
    VecA[2] = np.inner(XYZ[3]-XYZ[0],Amat)
    return nl.det(VecA)
    
def getRestTorsion(XYZ,Amat):
    VecA = np.empty((3,3))
    VecA[0] = np.inner(XYZ[1]-XYZ[0],Amat)
    VecA[1] = np.inner(XYZ[2]-XYZ[1],Amat)
    VecA[2] = np.inner(XYZ[3]-XYZ[2],Amat)
    D = nl.det(VecA)
    Mag = np.sqrt(np.sum(VecA*VecA,axis=1))
    P12 = np.sum(VecA[0]*VecA[1])/(Mag[0]*Mag[1])
    P13 = np.sum(VecA[0]*VecA[2])/(Mag[0]*Mag[2])
    P23 = np.sum(VecA[1]*VecA[2])/(Mag[1]*Mag[2])
    Ang = 1.0
    if abs(P12) < 1.0 and abs(P23) < 1.0:
        Ang = (P12*P23-P13)/(np.sqrt(1.-P12**2)*np.sqrt(1.-P23**2))
    TOR = (acosd(Ang)*D/abs(D)+720.)%360.
    return TOR
    
def calcTorsionEnergy(TOR,Coeff=[]):
    sum = 0.
    if len(Coeff):
        cof = np.reshape(Coeff,(3,3)).T
        delt = TOR-cof[1]
        delt = np.where(delt<-180.,delt+360.,delt)
        delt = np.where(delt>180.,delt-360.,delt)
        term = -cof[2]*delt**2
        val = cof[0]*np.exp(term/1000.0)
        pMax = cof[0][np.argmin(val)]
        Eval = np.sum(val)
        sum = Eval-pMax
    return sum,Eval

def getTorsionDeriv(XYZ,Amat,Coeff):
    deriv = np.zeros((len(XYZ),3))
    dx = 0.00001
    for j,xyz in enumerate(XYZ):
        for i,x in enumerate(np.array([[dx,0,0],[0,dx,0],[0,0,dx]])):
            XYZ[j] -= x
            tor = getRestTorsion(XYZ,Amat)
            p,d1 = calcTorsionEnergy(tor,Coeff)
            XYZ[j] += 2*x
            tor = getRestTorsion(XYZ,Amat)
            p,d2 = calcTorsionEnergy(tor,Coeff)            
            XYZ[j] -= x
            deriv[j][i] = (d2-d1)/(2*dx)
    return deriv.flatten()

def getRestRama(XYZ,Amat):
    phi = getRestTorsion(XYZ[:5],Amat)
    psi = getRestTorsion(XYZ[1:],Amat)
    return phi,psi
    
def calcRamaEnergy(phi,psi,Coeff=[]):
    sum = 0.
    if len(Coeff):
        cof = Coeff.T
        dPhi = phi-cof[1]
        dPhi = np.where(dPhi<-180.,dPhi+360.,dPhi)
        dPhi = np.where(dPhi>180.,dPhi-360.,dPhi)
        dPsi = psi-cof[2]
        dPsi = np.where(dPsi<-180.,dPsi+360.,dPsi)
        dPsi = np.where(dPsi>180.,dPsi-360.,dPsi)
        val = -cof[3]*dPhi**2-cof[4]*dPsi**2-2.0*cof[5]*dPhi*dPsi
        val = cof[0]*np.exp(val/1000.)
        pMax = cof[0][np.argmin(val)]
        Eval = np.sum(val)
        sum = Eval-pMax
    return sum,Eval

def getRamaDeriv(XYZ,Amat,Coeff):
    deriv = np.zeros((len(XYZ),3))
    dx = 0.00001
    for j,xyz in enumerate(XYZ):
        for i,x in enumerate(np.array([[dx,0,0],[0,dx,0],[0,0,dx]])):
            XYZ[j] -= x
            phi,psi = getRestRama(XYZ,Amat)
            p,d1 = calcRamaEnergy(phi,psi,Coeff)
            XYZ[j] += 2*x
            phi,psi = getRestRama(XYZ,Amat)
            p,d2 = calcRamaEnergy(phi,psi,Coeff)
            XYZ[j] -= x
            deriv[j][i] = (d2-d1)/(2*dx)
    return deriv.flatten()

def getRestPolefig(ODFln,SamSym,Grid):
    X,Y = np.meshgrid(np.linspace(1.,-1.,Grid),np.linspace(-1.,1.,Grid))
    R,P = np.sqrt(X**2+Y**2).flatten(),atan2d(Y,X).flatten()
    R = np.where(R <= 1.,2.*atand(R),0.0)
    Z = np.zeros_like(R)
    Z = G2lat.polfcal(ODFln,SamSym,R,P)
    Z = np.reshape(Z,(Grid,Grid))
    return np.reshape(R,(Grid,Grid)),np.reshape(P,(Grid,Grid)),Z

def getRestPolefigDerv(HKL,Grid,SHCoeff):
    pass
        
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
        Oxyz[i] -= dx
        d0 = calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)
        Oxyz[i] += 2*dx
        deriv[i] = (calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)-d0)/(2.*dx)
        Oxyz[i] -= dx
        Txyz[i] -= dx
        d0 = calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)
        Txyz[i] += 2*dx
        deriv[i+3] = (calcDist(Oxyz,Txyz,Tunit,inv,C,M,T,Amat)-d0)/(2.*dx)
        Txyz[i] -= dx
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
            OxA[i] -= dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            OxA[i] += 2*dx
            dadx[i] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/(2*dx)
            OxA[i] -= dx
            
            TxA[i] -= dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            TxA[i] += 2*dx
            dadx[i+3] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/(2*dx)
            TxA[i] -= dx
            
            TxB[i] -= dx
            a0 = calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)
            TxB[i] += 2*dx
            dadx[i+6] = (calcAngle(OxA,TxA,TxB,unitA,unitB,invA,CA,MA,TA,invB,CB,MB,TB,Amat)-a0)/(2*dx)
            TxB[i] -= dx
            
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
            Oatoms[ia][ix+1] -= dx
            a0 = calcTorsion(Oatoms,SyOps,Amat)
            Oatoms[ia][ix+1] += 2*dx
            dadx[i] = (calcTorsion(Oatoms,SyOps,Amat)-a0)/(2.*dx)
            Oatoms[ia][ix+1] -= dx            
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
            
def adjHKLmax(SGData,Hmax):
    if SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
        Hmax[0] = ((Hmax[0]+3)/6)*6
        Hmax[1] = ((Hmax[1]+3)/6)*6
        Hmax[2] = ((Hmax[2]+1)/4)*4
    else:
        Hmax[0] = ((Hmax[0]+2)/4)*4
        Hmax[1] = ((Hmax[1]+2)/4)*4
        Hmax[2] = ((Hmax[2]+1)/4)*4

def OmitMap(data,reflData):
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
    adjHKLmax(SGData,Hmax)
    Fhkl = np.zeros(shape=2*Hmax,dtype='c16')
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
                F = np.sqrt(Fosq)
                h,k,l = hkl+Hmax
                Fhkl[h,k,l] = F*phasep
                h,k,l = -hkl+Hmax
                Fhkl[h,k,l] = F*phasem
    rho = fft.fftn(fft.fftshift(Fhkl))/cell[6]
    print 'NB: this is just an Fobs map for now - under development'
    print 'Omit map time: %.4f'%(time.time()-time0),'no. elements: %d'%(Fhkl.size)
    print rho.shape
    mapData['rho'] = np.real(rho)
    mapData['rhoMax'] = max(np.max(mapData['rho']),-np.min(mapData['rho']))
    return mapData

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
    adjHKLmax(SGData,Hmax)
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
                
def findOffset(SGData,A,Fhkl):    
    if SGData['SpGrp'] == 'P 1':
        return [0,0,0]    
    hklShape = Fhkl.shape
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
    Hmax = 2*np.asarray(G2lat.getHKLmax(3.5,SGData,A),dtype='i')
    for F in Flist:
        hkl = np.unravel_index(Fdict[F],hklShape)
        iabsnt,mulp,Uniq,Phi = G2spc.GenHKLf(list(hkl-hklHalf),SGData)
        Uniq = np.array(Uniq,dtype='i')
        Phi = np.array(Phi)
        Uniq = np.concatenate((Uniq,-Uniq))+hklHalf         # put in Friedel pairs & make as index to Farray
        Phi = np.concatenate((Phi,-Phi))                      # and their phase shifts
        Fh0 = Fhkl[hkl[0],hkl[1],hkl[2]]
        ang0 = np.angle(Fh0,deg=True)/360.
        for H,phi in zip(Uniq,Phi)[1:]:
            ang = (np.angle(Fhkl[H[0],H[1],H[2]],deg=True)/360.-phi)
            dH = H-hkl
            dang = ang-ang0
            if np.any(np.abs(dH)-Hmax > 0):    #keep low order DHs
                continue
            DH.append(dH)
            Dphi.append((dang+.5) % 1.0)
        if i > 20 or len(DH) > 30:
            break
        i += 1
    DH = np.array(DH)
    print ' map offset no.of terms: %d from %d reflections'%(len(DH),len(Flist))
    Dphi = np.array(Dphi)
    steps = np.array(hklShape)
    X,Y,Z = np.mgrid[0:1:1./steps[0],0:1:1./steps[1],0:1:1./steps[2]]
    XYZ = np.array(zip(X.flatten(),Y.flatten(),Z.flatten()))
    Dang = (np.dot(XYZ,DH.T)+.5)%1.-Dphi
    Mmap = np.reshape(np.sum((Dang)**2,axis=1),newshape=steps)/len(DH)
    hist,bins = np.histogram(Mmap,bins=1000)
#    for i,item in enumerate(hist[:10]):
#        print item,bins[i]
    chisq = np.min(Mmap)
    DX = -np.array(np.unravel_index(np.argmin(Mmap),Mmap.shape))
    print ' map offset chi**2: %.3f, map offset: %d %d %d'%(chisq,DX[0],DX[1],DX[2])
#    print (np.dot(DX,DH.T)+.5)%1.-Dphi
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
    adjHKLmax(SGData,Hmax)
    Ehkl = np.zeros(shape=2*Hmax,dtype='c16')       #2X64bits per complex no.
    time0 = time.time()
    for ref in reflData:
        dsp = ref[4]
        if dsp >= dmin:
            ff = 0.1*Vol    #est. no. atoms for ~10A**3/atom
            if FFtable:
                SQ = 0.25/dsp**2
                ff *= G2el.ScatFac(FFtable,SQ)[0]
            if ref[8] > 0.:         #use only +ve Fobs**2
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
        CErho = np.real(fft.fftn(fft.fftshift(CEhkl)))*(1.+0j)
        CEsig = np.std(CErho)
        CFrho = np.where(np.real(CErho) >= flipData['k-factor']*CEsig,CErho,-CErho)
        CFrho = np.where(np.real(CErho) <= flipData['k-Max']*CEsig,CFrho,-CFrho)      #solves U atom problem! make 20. adjustible
        CFhkl = fft.ifftshift(fft.ifftn(CFrho))
        CFhkl = np.where(CFhkl,CFhkl,1.0)           #avoid divide by zero
        phase = CFhkl/np.absolute(CFhkl)
        CEhkl = np.absolute(Ehkl)*phase
        Ncyc += 1
        sumCF = np.sum(ma.array(np.absolute(CFhkl),mask=Emask))
        DEhkl = np.absolute(np.absolute(Ehkl)/sumE-np.absolute(CFhkl)/sumCF)
        Rcf = min(100.,np.sum(ma.array(DEhkl,mask=Emask)*100.))
        if Rcf < 5.:
            break
        GoOn = pgbar.Update(Rcf,newmsg='%s%8.3f%s\n%s %d'%('Residual Rcf =',Rcf,'%','No.cycles = ',Ncyc))[0]
        if not GoOn or Ncyc > 10000:
            break
    np.seterr(**old)
    print ' Charge flip time: %.4f'%(time.time()-time0),'no. elements: %d'%(Ehkl.size)
    CErho = np.real(fft.fftn(fft.fftshift(CEhkl)))
    print ' No.cycles = ',Ncyc,'Residual Rcf =%8.3f%s'%(Rcf,'%')+' Map size:',CErho.shape
    roll = findOffset(SGData,A,CEhkl)
        
    mapData['Rcf'] = Rcf
    mapData['rho'] = np.roll(np.roll(np.roll(CErho,roll[0],axis=0),roll[1],axis=1),roll[2],axis=2)
    mapData['rhoMax'] = max(np.max(mapData['rho']),-np.min(mapData['rho']))
    mapData['rollMap'] = [0,0,0]
    return mapData
    
def SearchMap(data):
    rollMap = lambda rho,roll: np.roll(np.roll(np.roll(rho,roll[0],axis=0),roll[1],axis=1),roll[2],axis=2)
    
    norm = 1./(np.sqrt(3.)*np.sqrt(2.*np.pi)**3)
    
    def noDuplicate(xyz,peaks,Amat):
        XYZ = np.inner(Amat,xyz)
        if True in [np.allclose(XYZ,np.inner(Amat,peak),atol=0.5) for peak in peaks]:
            print ' Peak',xyz,' <0.5A from another peak'
            return False
        return True
                            
    def fixSpecialPos(xyz,SGData,Amat):
        equivs = G2spc.GenAtom(xyz,SGData,Move=True)
        X = []
        xyzs = [equiv[0] for equiv in equivs]
        for x in xyzs:
            if np.sqrt(np.sum(np.inner(Amat,xyz-x)**2,axis=0))<0.5:
                X.append(x)
        if len(X) > 1:
            return np.average(X,axis=0)
        else:
            return xyz
        
    def rhoCalc(parms,rX,rY,rZ,res,SGLaue):
        Mag,x0,y0,z0,sig = parms
        z = -((x0-rX)**2+(y0-rY)**2+(z0-rZ)**2)/(2.*sig**2)
#        return norm*Mag*np.exp(z)/(sig*res**3)     #not slower but some faults in LS
        return norm*Mag*(1.+z+z**2/2.)/(sig*res**3)
        
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
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    drawingData = data['Drawing']
    peaks = []
    mags = []
    dzeros = []
    try:
        mapData = generalData['Map']
        contLevel = mapData['cutOff']*mapData['rhoMax']/100.
        rho = copy.copy(mapData['rho'])     #don't mess up original
        mapHalf = np.array(rho.shape)/2
        res = mapData['Resolution']
        incre = np.array(rho.shape,dtype=np.float)
        step = max(1.0,1./res)+1
        steps = np.array(3*[step,])
    except KeyError:
        print '**** ERROR - Fourier map not defined'
        return peaks,mags
    rhoMask = ma.array(rho,mask=(rho<contLevel))
    indices = (-1,0,1)
    rolls = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    for roll in rolls:
        if np.any(roll):
            rhoMask = ma.array(rhoMask,mask=(rhoMask-rollMap(rho,roll)<=0.))
    indx = np.transpose(rhoMask.nonzero())
    peaks = indx/incre
    mags = rhoMask[rhoMask.nonzero()]
    for i,[ind,peak,mag] in enumerate(zip(indx,peaks,mags)):
        rho = rollMap(rho,ind)
        rMM = mapHalf-steps
        rMP = mapHalf+steps+1
        rhoPeak = rho[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]]
        peakInt = np.sum(rhoPeak)*res**3
        rX,rY,rZ = np.mgrid[rMM[0]:rMP[0],rMM[1]:rMP[1],rMM[2]:rMP[2]]
        x0 = [peakInt,mapHalf[0],mapHalf[1],mapHalf[2],2.0]          #magnitude, position & width(sig)
        result = HessianLSQ(peakFunc,x0,Hess=peakHess,
            args=(rX,rY,rZ,rhoPeak,res,SGData['SGLaue']),ftol=.01,maxcyc=10)
        x1 = result[0]
        if not np.any(x1 < 0):
            mag = x1[0]
            peak = (np.array(x1[1:4])-ind)/incre
        peak = fixSpecialPos(peak,SGData,Amat)
        rho = rollMap(rho,-ind)        
    dzeros = np.sqrt(np.sum(np.inner(Amat,peaks)**2,axis=0))
    return np.array(peaks),np.array([mags,]).T,np.array([dzeros,]).T
    
def sortArray(data,pos,reverse=False):
    #data is a list of items
    #sort by pos in list; reverse if True
    T = []
    for i,M in enumerate(data):
        T.append((M[pos],i))
    D = dict(zip(T,data))
    T.sort()
    if reverse:
        T.reverse()
    X = []
    for key in T:
        X.append(D[key])
    return X

def PeaksEquiv(data,Ind):

    def Duplicate(xyz,peaks,Amat):
        if True in [np.allclose(np.inner(Amat,xyz),np.inner(Amat,peak),atol=0.5) for peak in peaks]:
            return True
        return False
                            
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    A = G2lat.cell2A(cell)
    SGData = generalData['SGData']
    mapPeaks = data['Map Peaks']
    XYZ = np.array([xyz[1:4] for xyz in mapPeaks])
    Indx = {}
    for ind in Ind:
        xyz = np.array(mapPeaks[ind][1:4])
        xyzs = np.array([equiv[0] for equiv in G2spc.GenAtom(xyz,SGData,Move=True)])
#        for x in xyzs: print x 
        for jnd,xyz in enumerate(XYZ):       
            Indx[jnd] = Duplicate(xyz,xyzs,Amat)
    Ind = []
    for ind in Indx:
        if Indx[ind]:
            Ind.append(ind)
    return Ind
                
def PeaksUnique(data,Ind):
#    XYZE = np.array([[equiv[0] for equiv in G2spc.GenAtom(xyz[1:4],SGData,Move=True)] for xyz in mapPeaks]) #keep this!!

    def noDuplicate(xyz,peaks,Amat):
        if True in [np.allclose(np.inner(Amat,xyz),np.inner(Amat,peak),atol=0.5) for peak in peaks]:
            return False
        return True
                            
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    A = G2lat.cell2A(cell)
    SGData = generalData['SGData']
    mapPeaks = data['Map Peaks']
    Indx = {}
    XYZ = {}
    for ind in Ind:
        XYZ[ind] = np.array(mapPeaks[ind][1:4])
        Indx[ind] = True
    for ind in Ind:
        if Indx[ind]:
            xyz = XYZ[ind]
            for jnd in Ind:
                if ind != jnd and Indx[jnd]:                        
                    Equiv = G2spc.GenAtom(XYZ[jnd],SGData,Move=True)
                    xyzs = np.array([equiv[0] for equiv in Equiv])
                    Indx[jnd] = noDuplicate(xyz,xyzs,Amat)
    Ind = []
    for ind in Indx:
        if Indx[ind]:
            Ind.append(ind)
    return Ind
    
def setPeakparms(Parms,Parms2,pos,mag,ifQ=False,useFit=False):
    ind = 0
    if useFit:
        ind = 1
    ins = {}
    if 'C' in Parms['Type'][0]:                            #CW data - TOF later in an elif
        for x in ['U','V','W','X','Y']:
            ins[x] = Parms[x][ind]
        if ifQ:                              #qplot - convert back to 2-theta
            pos = 2.0*asind(pos*wave/(4*math.pi))
        sig = ins['U']*tand(pos/2.0)**2+ins['V']*tand(pos/2.0)+ins['W']
        gam = ins['X']/cosd(pos/2.0)+ins['Y']*tand(pos/2.0)           
        XY = [pos,0, mag,1, sig,0, gam,0]       #default refine intensity 1st
    else:
        if ifQ:
            dsp = 2.*np.pi/pos
            pos = Parms['difC']*dsp
        else:
            dsp = pos/Parms['difC'][1]
        if 'Pdabc' in Parms2:
            for x in ['sig-0','sig-1','X','Y']:
                ins[x] = Parms[x][ind]
            Pdabc = Parms2['Pdabc'].T
            alp = np.interp(dsp,Pdabc[0],Pdabc[1])
            bet = np.interp(dsp,Pdabc[0],Pdabc[2])
        else:
            for x in ['alpha','beta-0','beta-1','sig-0','sig-1','X','Y']:
                ins[x] = Parms[x][ind]
            alp = ins['alpha']/dsp
            bet = ins['beta-0']+ins['beta-1']/dsp**4
        sig = ins['sig-0']+ins['sig-1']*dsp**2
        gam = ins['X']*dsp+ins['Y']*dsp**2
        XY = [pos,0,mag,1,alp,0,bet,0,sig,0,gam,0]
    return XY
        
def getWave(Parms):
    try:
        return Parms['Lam'][1]
    except KeyError:
        return Parms['Lam1'][1]
    
def prodQQ(QA,QB):
    ''' Grassman quaternion product
        QA,QB quaternions; q=r+ai+bj+ck
    '''
    D = np.zeros(4)
    D[0] = QA[0]*QB[0]-QA[1]*QB[1]-QA[2]*QB[2]-QA[3]*QB[3]
    D[1] = QA[0]*QB[1]+QA[1]*QB[0]+QA[2]*QB[3]-QA[3]*QB[2]
    D[2] = QA[0]*QB[2]-QA[1]*QB[3]+QA[2]*QB[0]+QA[3]*QB[1]
    D[3] = QA[0]*QB[3]+QA[1]*QB[2]-QA[2]*QB[1]+QA[3]*QB[0]
    return D
    
def normQ(QA):
    ''' get length of quaternion & normalize it
        q=r+ai+bj+ck
    '''
    n = np.sqrt(np.sum(np.array(QA)**2))
    return QA/n
    
def invQ(Q):
    '''
        get inverse of quaternion
        q=r+ai+bj+ck; q* = r-ai-bj-ck
    '''
    return Q*np.array([1,-1,-1,-1])
    
def prodQVQ(Q,V):
    ''' compute the quaternion vector rotation qvq-1 = v'
        q=r+ai+bj+ck
    '''
    VP = np.zeros(3)
    T2 = Q[0]*Q[1]
    T3 = Q[0]*Q[2]
    T4 = Q[0]*Q[3]
    T5 = -Q[1]*Q[1]
    T6 = Q[1]*Q[2]
    T7 = Q[1]*Q[3]
    T8 = -Q[2]*Q[2]
    T9 = Q[2]*Q[3]
    T10 = -Q[3]*Q[3]
    VP[0] = 2.*((T8+T10)*V[0]+(T6-T4)*V[1]+(T3+T7)*V[2])+V[0]
    VP[1] = 2.*((T4+T6)*V[0]+(T5+T10)*V[1]+(T9-T2)*V[2])+V[1]
    VP[2] = 2.*((T7-T3)*V[0]+(T2+T9)*V[1]+(T5+T8)*V[2])+V[2] 
    return VP   
    
def Q2Mat(Q):
    ''' make rotation matrix from quaternion
        q=r+ai+bj+ck
    '''
    aa = Q[0]**2
    ab = Q[0]*Q[1]
    ac = Q[0]*Q[2]
    ad = Q[0]*Q[3]
    bb = Q[1]**2
    bc = Q[1]*Q[2]
    bd = Q[1]*Q[3]
    cc = Q[2]**2
    cd = Q[2]*Q[3]
    dd = Q[3]**2
    M = [[aa+bb-cc-dd, 2.*(bc-ad),  2.*(ac+bd)],
        [2*(ad+bc),   aa-bb+cc-dd,  2.*(cd-ab)],
        [2*(bd-ac),    2.*(ab+cd), aa-bb-cc+dd]]
    return np.array(M)
    
def AV2Q(A,V):
    ''' convert angle (radians) & vector to quaternion
        q=r+ai+bj+ck
    '''
    Q = np.zeros(4)
    d = np.sqrt(np.sum(np.array(V)**2))
    if d:
        V /= d
    else:
        return [0.,0.,0.,1.]    #identity
    p = A/2.
    Q[0] = np.cos(p)
    Q[1:4] = V*np.sin(p)
    return Q
    
def AVdeg2Q(A,V):
    ''' convert angle (degrees) & vector to quaternion
        q=r+ai+bj+ck
    '''
    Q = np.zeros(4)
    d = np.sqrt(np.sum(np.array(V)**2))
    if d:
        V /= d
    else:
        return [0.,0.,0.,1.]    #identity
    p = A/2.
    Q[0] = cosd(p)
    Q[1:4] = V*sind(p)
    return Q
    
def Q2AVdeg(Q):
    ''' convert quaternion to angle (degrees 0-360) & normalized vector
        q=r+ai+bj+ck
    '''
    A = 2.*acosd(Q[0])
    V = np.array(Q[1:])
    if nl.norm(Q[1:]):
        V = Q[1:]/nl.norm(Q[1:])
    else:
        A = 0.
        V = np.array([0.,0.,1.])
    return A,V
    
def Q2AV(Q):
    ''' convert quaternion to angle (radians 0-2pi) & normalized vector
        q=r+ai+bj+ck
    '''
    A = 2.*np.arccos(Q[0])
    V = np.array(Q[1:])
    if nl.norm(Q[1:]):
        V = Q[1:]/nl.norm(Q[1:])
    else:
        A = 0.
        V = np.array([0.,0.,1.])
    return A,V
    
def makeQuat(A,B,C):
    ''' Make quaternion from rotation of A vector to B vector about C axis
        A,B,C are np.array Cartesian 3-vectors
    Returns quaternion & rotation angle in radians
        q=r+ai+bj+ck
    '''

    V1 = np.cross(A,C)
    V2 = np.cross(B,C)
    if nl.norm(V1)*nl.norm(V2):
        V1 /= nl.norm(V1)
        V2 /= nl.norm(V2)
        V3 = np.cross(V1,V2)
    else:
        V3 = np.zeros(3)
    Q = np.array([0.,0.,0.,1.])
    D = 0.
    if nl.norm(V3):
        V3 /= nl.norm(V3)
        D1 = min(1.0,max(-1.0,np.vdot(V1,V2)))
        D = np.arccos(D1)/2.0
        V1 = C-V3
        V2 = C+V3
        DM = nl.norm(V1)
        DP = nl.norm(V2)
        S = np.sin(D)
        Q[0] = np.cos(D)
        Q[1:] = V3*S
        D *= 2.
        if DM > DP:
            D *= -1.
    return Q,D
    
