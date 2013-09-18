# -*- coding: utf-8 -*-
#GSASIImath - major mathematics routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIImath: computation module*
================================

Routines for least-squares minimization and other stuff

'''
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
import GSASIIpwd as G2pwd
import numpy.fft as fft
import pypowder as pyd

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atand = lambda x: 180.*np.arctan(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

def HessianLSQ(func,x0,Hess,args=(),ftol=1.49012e-8,xtol=1.49012e-8, maxcyc=0,Print=False):
    
    """
    Minimize the sum of squares of a function (:math:`f`) evaluated on a series of
    values (y): :math:`\sum_{y=0}^{N_{obs}} f(y,{args})`
    
    ::
    
                    Nobs
        x = arg min(sum(func(y)**2,axis=0))
                    y=0

    :param function func: callable method or function
        should take at least one (possibly length N vector) argument and
        returns M floating point numbers.
    :param np.ndarray x0: The starting estimate for the minimization of length N
    :param function Hess: callable method or function
        A required function or method to compute the weighted vector and Hessian for func.
        It must be a symmetric NxN array
    :param tuple args: Any extra arguments to func are placed in this tuple.
    :param float ftol: Relative error desired in the sum of squares.
    :param float xtol: Relative error desired in the approximate solution.
    :param int maxcyc: The maximum number of cycles of refinement to execute, if -1 refine 
        until other limits are met (ftol, xtol)
    :param bool Print: True for printing results (residuals & times) by cycle

    :returns: (x,cov_x,infodict) where

      * x : ndarray
        The solution (or the result of the last iteration for an unsuccessful
        call).
      * cov_x : ndarray
        Uses the fjac and ipvt optional outputs to construct an
        estimate of the jacobian around the solution.  ``None`` if a
        singular matrix encountered (indicates very flat curvature in
        some direction).  This matrix must be multiplied by the
        residual standard deviation to get the covariance of the
        parameter estimates -- see curve_fit.
      * infodict : dict
        a dictionary of optional outputs with the keys:

         * 'fvec' : the function evaluated at the output
         * 'num cyc':
         * 'nfev':
         * 'lamMax':
         * 'psing':
            
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
    if Print:
        print ' Hessian refinement on %d variables:'%(n)
    while icycle < maxcyc:
        time0 = time.time()
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
            if lam > 10.e5:
                print 'ouch #3 chisq1 ',chisq1,' stuck > chisq0 ',chisq0
                break
        lamMax = max(lamMax,lam)
        if Print:
            print ' Cycle: %d, Time: %.2fs, Chi**2: %.3g, Lambda: %.3g'%(icycle,time.time()-time0,chisq1,lam)
        if (chisq0-chisq1)/chisq0 < ftol:
            break
        icycle += 1
    M = func(x0,*args)
    nfev += 1
    Yvec,Amat = Hess(x0,*args)
    Amatlam = Amat          #*(One+Lam)/Anorm
    try:
        Bmat = nl.inv(Amatlam)          #*(One+Lam)/Anorm
        return [x0,Bmat,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':[]}]
    except nl.LinAlgError:
        print 'ouch #2 linear algebra error in LS'
        psing = []
        if maxcyc:
            psing = list(np.where(np.diag(nl.qr(Amat)[1]) < 1.e-14)[0])
        return [x0,None,{'num cyc':icycle,'fvec':M,'nfev':nfev,'lamMax':lamMax,'psing':psing}]
   
def getVCov(varyNames,varyList,covMatrix):
    '''obtain variance-covariance terms for a set of variables. NB: the varyList 
    and covMatrix were saved by the last least squares refinement so they must match.
    
    :param list varyNames: variable names to find v-cov matric for
    :param list varyList: full list of all variables in v-cov matrix
    :param nparray covMatrix: full variance-covariance matrix from the last 
     least squares refinement
    
    :returns: nparray vcov: variance-covariance matrix for the variables given
     in varyNames
    
    '''
    vcov = np.zeros((len(varyNames),len(varyNames)))
    for i1,name1 in enumerate(varyNames):
        for i2,name2 in enumerate(varyNames):
            try:
                vcov[i1][i2] = covMatrix[varyList.index(name1)][varyList.index(name2)]
            except ValueError:
                vcov[i1][i2] = 0.0
    return vcov

def FindAtomIndexByIDs(atomData,IDs,Draw=True):
    '''finds the set of atom array indices for a list of atom IDs. Will search 
    either the Atom table or the drawAtom table.
    
    :param list atomData: Atom or drawAtom table containting coordinates, etc.
    :param list IDs: atom IDs to be found
    :param bool Draw: True if drawAtom table to be searched; False if Atom table
      is searched
    
    :returns: list indx: atom (or drawAtom) indices
    
    '''
    indx = []
    for i,atom in enumerate(atomData):
        if Draw and atom[-3] in IDs:
            indx.append(i)
        elif atom[-1] in IDs:
            indx.append(i)
    return indx

def FillAtomLookUp(atomData):
    '''create a dictionary of atom indexes with atom IDs as keys
    
    :param list atomData: Atom table to be used
    
    :returns: dict atomLookUp: dictionary of atom indexes with atom IDs as keys
    
    '''
    atomLookUp = {}
    for iatm,atom in enumerate(atomData):
        atomLookUp[atom[-1]] = iatm
    return atomLookUp

def GetAtomsById(atomData,atomLookUp,IdList):
    '''gets a list of atoms from Atom table that match a set of atom IDs
    
    :param list atomData: Atom table to be used
    :param dict atomLookUp: dictionary of atom indexes with atom IDs as keys
    :param list IdList: atom IDs to be found
    
    :returns: list atoms: list of atoms found
    
    '''
    atoms = []
    for id in IdList:
        atoms.append(atomData[atomLookUp[id]])
    return atoms
    
def GetAtomItemsById(atomData,atomLookUp,IdList,itemLoc,numItems=1):
    '''gets atom parameters for atoms using atom IDs
    
    :param list atomData: Atom table to be used
    :param dict atomLookUp: dictionary of atom indexes with atom IDs as keys
    :param list IdList: atom IDs to be found
    :param int itemLoc: pointer to desired 1st item in an atom table entry
    :param int numItems: number of items to be retrieved
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    pfx = [str(pId)+'::A'+i+':' for i in ['x','y','z']]
    dpfx = [str(pId)+'::dA'+i+':' for i in ['x','y','z']]
    XYZ = []
    for ind in indx:
        names = [pfx[i]+str(AtLookup[ind]) for i in range(3)]
        dnames = [dpfx[i]+str(AtLookup[ind]) for i in range(3)]
        XYZ.append([parmDict[name]+parmDict[dname] for name,dname in zip(names,dnames)])
    return XYZ

def AtomUij2TLS(atomData,atPtrs,Amat,Bmat,rbObj):   #unfinished & not used
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    for atom in atomData:
        XYZ = np.inner(Amat,atom[cx:cx+3])
        if atom[cia] == 'A':
            UIJ = atom[cia+2:cia+8]
                
def TLS2Uij(xyz,g,Amat,rbObj):    #not used anywhere, but could be?
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
        
def AtomTLS2UIJ(atomData,atPtrs,Amat,rbObj):    #not used anywhere, but could be?
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''gets distance from position xyz to all XYZ, xyz & XYZ are np.array 
        and are in crystal coordinates; Amat is crystal to Cart matrix
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return np.sqrt(np.sum(np.inner(Amat,XYZ-xyz)**2,axis=0))

def getAtomXYZ(atoms,cx):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    XYZ = []
    for atom in atoms:
        XYZ.append(atom[cx:cx+3])
    return np.array(XYZ)

def RotateRBXYZ(Bmat,Cart,oriQ):
    '''rotate & transform cartesian coordinates to crystallographic ones
    no translation applied. To be used for numerical derivatives 
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    ''' returns crystal coordinates for atoms described by RBObj
    '''
    XYZ = np.zeros_like(Cart)
    for i,xyz in enumerate(Cart):
        X = prodQVQ(oriQ,xyz)
        XYZ[i] = np.inner(Bmat,X)
    return XYZ

def UpdateRBXYZ(Bmat,RBObj,RBData,RBType):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    ''' returns crystal coordinates for atoms described by RBObj
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

def UpdateMCSAxyz(Bmat,MCSA):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    xyz = []
    atTypes = []
    iatm = 0
    for model in MCSA['Models'][1:]:        #skip the MD model
        if model['Type'] == 'Atom':
            xyz.append(model['Pos'][0])
            atTypes.append(model['atType'])
            iatm += 1
        else:
            RBRes = MCSA['rbData'][model['Type']][model['RBId']]
            Pos = np.array(model['Pos'][0])
            Ori = np.array(model['Ori'][0])
            Qori = AVdeg2Q(Ori[0],Ori[1:])
            if model['Type'] == 'Vector':
                vecs = RBRes['rbVect']
                mags = RBRes['VectMag']
                Cart = np.zeros_like(vecs[0])
                for vec,mag in zip(vecs,mags):
                    Cart += vec*mag
            elif model['Type'] == 'Residue':
                Cart = np.array(RBRes['rbXYZ'])
                for itor,seq in enumerate(RBRes['rbSeq']):
                    QuatA = AVdeg2Q(model['Tor'][0][itor],Cart[seq[0]]-Cart[seq[1]])
                    for ride in seq[3]:
                        Cart[ride] = prodQVQ(QuatA,Cart[ride]-Cart[seq[1]])+Cart[seq[1]]
            if model['MolCent'][1]:
                Cart -= model['MolCent'][0]
            for i,x in enumerate(Cart):
                xyz.append(np.inner(Bmat,prodQVQ(Qori,x))+Pos)
                atType = RBRes['rbTypes'][i]
                atTypes.append(atType)
                iatm += 1
    return np.array(xyz),atTypes
    
def SetMolCent(model,RBData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    rideList = []
    RBRes = RBData[model['Type']][model['RBId']]
    if model['Type'] == 'Vector':
        vecs = RBRes['rbVect']
        mags = RBRes['VectMag']
        Cart = np.zeros_like(vecs[0])
        for vec,mag in zip(vecs,mags):
            Cart += vec*mag
    elif model['Type'] == 'Residue':
        Cart = np.array(RBRes['rbXYZ'])
        for seq in RBRes['rbSeq']:
            rideList += seq[3]
    centList = set(range(len(Cart)))-set(rideList)
    cent = np.zeros(3)
    for i in centList:
        cent += Cart[i]
    model['MolCent'][0] = cent/len(centList)    
    
def UpdateRBUIJ(Bmat,Cart,RBObj):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    ''' returns atom I/A, Uiso or UIJ for atoms at XYZ as described by RBObj
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
    g = nl.inv(np.inner(Bmat,Bmat))
    gvec = np.sqrt(np.array([g[0][0]**2,g[1][1]**2,g[2][2]**2,
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
            beta = np.inner(np.inner(Bmat.T,Umat),Bmat)
            Uout.append(['A',0.0,]+list(G2lat.UijtoU6(beta)*gvec))
        else:
            Uout.append(['N',])
    return Uout
    
def GetSHCoeff(pId,parmDict,SHkeys):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    SHCoeff = {}
    for shkey in SHkeys:
        shname = str(pId)+'::'+shkey
        SHCoeff[shkey] = parmDict[shname]
    return SHCoeff
        
def getMass(generalData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    'Computes mass of unit cell contents'
    mass = 0.
    for i,elem in enumerate(generalData['AtomTypes']):
        mass += generalData['NoAtoms'][elem]*generalData['AtomMass'][i]
    return mass    

def getDensity(generalData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    mass = getMass(generalData)
    Volume = generalData['Cell'][7]
    density = mass/(0.6022137*Volume)
    return density,Volume/mass
    
def getWave(Parms):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    try:
        return Parms['Lam'][1]
    except KeyError:
        return Parms['Lam1'][1]
    
################################################################################
##### distance, angle, planes, torsion stuff stuff
################################################################################

def getSyXYZ(XYZ,ops,SGData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return np.sqrt(np.sum(np.inner(Amat,(XYZ[1]-XYZ[0]))**2))
    
def getRestDeriv(Func,XYZ,Amat,ops,SGData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    VecA = np.empty((3,3))    
    VecA[0] = np.inner(XYZ[1]-XYZ[0],Amat)
    VecA[1] = np.inner(XYZ[2]-XYZ[0],Amat)
    VecA[2] = np.inner(XYZ[3]-XYZ[0],Amat)
    return nl.det(VecA)
    
def getRestTorsion(XYZ,Amat):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''Computes a pair of torsion angles in a 5 atom string
    
    :param nparray XYZ: crystallographic coordinates of 5 atoms
    :param nparray Amat: crystal to cartesian transformation matrix
    
    :returns: list (phi,psi) two torsion angles in degrees
    
    '''
    phi = getRestTorsion(XYZ[:5],Amat)
    psi = getRestTorsion(XYZ[1:],Amat)
    return phi,psi
    
def calcRamaEnergy(phi,psi,Coeff=[]):
    '''Computes pseudo potential energy from a pair of torsion angles and a
    numerical description of the potential energy surface. Used to create 
    penalty function in LS refinement:     
    :math:`Eval(\\phi,\\psi) = C[0]*exp(-V/1000)`

    where :math:`V = -C[3] * (\\phi-C[1])^2 - C[4]*(\\psi-C[2])^2 - 2*(\\phi-C[1])*(\\psi-C[2])`
    
    :param float phi: first torsion angle (:math:`\\phi`)
    :param float psi: second torsion angle (:math:`\\psi`)
    :param list Coeff: pseudo potential coefficients
    
    :returns: list (sum,Eval): pseudo-potential difference from minimum & value;
      sum is used for penalty function.
    
    '''
    sum = 0.
    Eval = 0.
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
    '''Computes numerical derivatives of torsion angle pair pseudo potential
    with respect of crystallographic atom coordinates of the 5 atom sequence 
    
    :param nparray XYZ: crystallographic coordinates of 5 atoms
    :param nparray Amat: crystal to cartesian transformation matrix
    :param list Coeff: pseudo potential coefficients
    
    :returns: list (deriv) derivatives of pseudopotential with respect to 5 atom
     crystallographic xyz coordinates.
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    X,Y = np.meshgrid(np.linspace(1.,-1.,Grid),np.linspace(-1.,1.,Grid))
    R,P = np.sqrt(X**2+Y**2).flatten(),atan2d(Y,X).flatten()
    R = np.where(R <= 1.,2.*atand(R),0.0)
    Z = np.zeros_like(R)
    Z = G2lat.polfcal(ODFln,SamSym,R,P)
    Z = np.reshape(Z,(Grid,Grid))
    return np.reshape(R,(Grid,Grid)),np.reshape(P,(Grid,Grid)),Z

def getRestPolefigDerv(HKL,Grid,SHCoeff):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    pass
        
def getDistDerv(Oxyz,Txyz,Amat,Tunit,Top,SGData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''

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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''

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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''

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
        
def ValEsd(value,esd=0,nTZ=False):
    '''Format a floating point number with a given level of precision or
    with in crystallographic format with a "esd", as value(esd). If esd is
    negative the number is formatted with the level of significant figures
    appropriate if abs(esd) were the esd, but the esd is not included.
    if the esd is zero, approximately 6 significant figures are printed.
    nTZ=True causes "extra" zeros to be removed after the decimal place.
    for example:

      * "1.235(3)" for value=1.2346 & esd=0.003
      * "1.235(3)e4" for value=12346. & esd=30
      * "1.235(3)e6" for value=0.12346e7 & esd=3000
      * "1.235" for value=1.2346 & esd=-0.003
      * "1.240" for value=1.2395 & esd=-0.003
      * "1.24" for value=1.2395 & esd=-0.003 with nTZ=True
      * "1.23460" for value=1.2346 & esd=0.0

    :param float value: number to be formatted
    :param float esd: uncertainty or if esd < 0, specifies level of
      precision to be shown e.g. esd=-0.01 gives 2 places beyond decimal
    :param bool nTZ: True to remove trailing zeros (default is False)
    :returns: value(esd) or value as a string

    '''
    # Note: this routine is Python 3 compatible -- I think
    if math.isnan(value): # invalid value, bail out
        return '?'
    if math.isnan(esd): # invalid esd, treat as zero
        esd = 0
        esdoff = 5
    elif esd != 0:
        # transform the esd to a one or two digit integer
        l = math.log10(abs(esd)) % 1
        # cut off of 19 1.9==>(19) but 1.95==>(2) N.B. log10(1.95) = 0.2900...
        if l < 0.290034611362518: l+= 1.        
        intesd = int(round(10**l)) # esd as integer
        # determine the number of digits offset for the esd
        esdoff = int(round(math.log10(intesd*1./abs(esd))))
    else:
        esdoff = 5
    valoff = 0
    if abs(value) < abs(esdoff): # value is effectively zero
        pass
    elif esdoff < 0 or abs(value) > 1.0e6 or abs(value) < 1.0e-4: # use scientific notation
        # where the digit offset is to the left of the decimal place or where too many
        # digits are needed
        if abs(value) > 1:
            valoff = int(math.log10(abs(value)))
        elif abs(value) > 0:
            valoff = int(math.log10(abs(value))-0.9999999)
        else:
            valoff = 0
    if esd != 0:
        if valoff+esdoff < 0:
            valoff = esdoff = 0
        out = ("{:."+str(valoff+esdoff)+"f}").format(value/10**valoff) # format the value
    elif valoff != 0: # esd = 0; exponential notation ==> esdoff decimal places
        out = ("{:."+str(esdoff)+"f}").format(value/10**valoff) # format the value
    else: # esd = 0; non-exponential notation ==> esdoff+1 significant digits
        if abs(value) > 0:            
            extra = -math.log10(abs(value))
        else:
            extra = 0
        if extra > 0: extra += 1
        out = ("{:."+str(max(0,esdoff+int(extra)))+"f}").format(value) # format the value
    if esd > 0:
        out += ("({:d})").format(intesd)  # add the esd
    elif nTZ and '.' in out:
        out = out.rstrip('0')  # strip zeros to right of decimal
        out = out.rstrip('.')  # and decimal place when not needed
    if valoff != 0:
        out += ("e{:d}").format(valoff) # add an exponent, when needed
    return out

################################################################################
##### Fourier & charge flip stuff
################################################################################

def adjHKLmax(SGData,Hmax):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    if SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
        Hmax[0] = ((Hmax[0]+3)/6)*6
        Hmax[1] = ((Hmax[1]+3)/6)*6
        Hmax[2] = ((Hmax[2]+1)/4)*4
    else:
        Hmax[0] = ((Hmax[0]+2)/4)*4
        Hmax[1] = ((Hmax[1]+2)/4)*4
        Hmax[2] = ((Hmax[2]+1)/4)*4

def OmitMap(data,reflData):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
                    F = np.where(Fosq>0.,np.sqrt(Fosq),0.)
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
                    dF = np.where(Fosq>0.,np.sqrt(Fosq),0.)-np.sqrt(Fcsq)
                    h,k,l = hkl+Hmax
                    Fhkl[h,k,l] = dF*phasep
                    h,k,l = -hkl+Hmax
                    Fhkl[h,k,l] = dF*phasem
                elif '2*Fo-Fc' in mapData['MapType']:
                    F = 2.*np.where(Fosq>0.,np.sqrt(Fosq),0.)-np.sqrt(Fcsq)
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
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
    '''Does a search of a density map for peaks meeting the criterion of peak
    height is greater than mapData['cutOff']/100 of mapData['rhoMax'] where 
    mapData is data['General']['mapData']; the map is also in mapData.

    :param data: the phase data structure

    :returns: (peaks,mags,dzeros) where

        * peaks : ndarray
          x,y,z positions of the peaks found in the map
        * mags : ndarray
          the magnitudes of the peaks
        * dzeros : ndarray
          the distance of the peaks from  the unit cell origin

    '''        
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
    '''data is a list of items
    sort by pos in list; reverse if True
    '''
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
    '''Find the equivalent map peaks for those selected. Works on the 
    contents of data['Map Peaks'].

    :param data: the phase data structure
    :param list Ind: list of selected peak indices
    :returns: augmented list of peaks including those related by symmetry to the
      ones in Ind

    '''        
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
    '''Finds the symmetry unique set of peaks from those selected. Works on the 
    contents of data['Map Peaks']. 

    :param data: the phase data structure
    :param list Ind: list of selected peak indices
    :returns: the list of symmetry unique peaks from among those given in Ind

    '''        
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
    
################################################################################
##### single peak fitting profile fxn stuff
################################################################################

def getCWsig(ins,pos):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    tp = tand(pos/2.0)
    return ins['U']*tp**2+ins['V']*tp+ins['W']
    
def getCWsigDeriv(pos):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    tp = tand(pos/2.0)
    return tp**2,tp,1.0
    
def getCWgam(ins,pos):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return ins['X']/cosd(pos/2.0)+ins['Y']*tand(pos/2.0)
    
def getCWgamDeriv(pos):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return 1./cosd(pos/2.0),tand(pos/2.0)
    
def getTOFsig(ins,dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return ins['sig-0']+ins['sig-1']*dsp**2+ins['sig-q']*dsp
    
def getTOFsigDeriv(dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return 1.0,dsp**2,dsp
    
def getTOFgamma(ins,dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return ins['X']*dsp+ins['Y']*dsp**2
    
def getTOFgammaDeriv(dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return dsp,dsp**2
    
def getTOFbeta(ins,dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return ins['beta-0']+ins['beta-1']/dsp**4+ins['beta-q']/dsp
    
def getTOFbetaDeriv(dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return 1.0,1./dsp**4,1./dsp
    
def getTOFalpha(ins,dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return ins['alpha']/dsp
    
def getTOFalphaDeriv(dsp):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    return 1./dsp
    
def setPeakparms(Parms,Parms2,pos,mag,ifQ=False,useFit=False):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    ind = 0
    if useFit:
        ind = 1
    ins = {}
    if 'C' in Parms['Type'][0]:                            #CW data - TOF later in an elif
        for x in ['U','V','W','X','Y']:
            ins[x] = Parms[x][ind]
        if ifQ:                              #qplot - convert back to 2-theta
            pos = 2.0*asind(pos*wave/(4*math.pi))
        sig = getCWsig(ins,pos)
        gam = getCWgam(ins,pos)           
        XY = [pos,0, mag,1, sig,0, gam,0]       #default refine intensity 1st
    else:
        if ifQ:
            dsp = 2.*np.pi/pos
            pos = Parms['difC']*dsp
        else:
            dsp = pos/Parms['difC'][1]
        if 'Pdabc' in Parms2:
            for x in ['sig-0','sig-1','sig-q','X','Y']:
                ins[x] = Parms[x][ind]
            Pdabc = Parms2['Pdabc'].T
            alp = np.interp(dsp,Pdabc[0],Pdabc[1])
            bet = np.interp(dsp,Pdabc[0],Pdabc[2])
        else:
            for x in ['alpha','beta-0','beta-1','beta-q','sig-0','sig-1','sig-q','X','Y']:
                ins[x] = Parms[x][ind]
            alp = getTOFalpha(ins,dsp)
            bet = getTOFbeta(ins,dsp)
        sig = getTOFsig(ins,dsp)
        gam = getTOFgamma(ins,dsp)
        XY = [pos,0,mag,1,alp,0,bet,0,sig,0,gam,0]
    return XY
    
################################################################################
##### MC/SA stuff
################################################################################

#scipy/optimize/anneal.py code modified by R. Von Dreele 2013
# Original Author: Travis Oliphant 2002
# Bug-fixes in 2006 by Tim Leslie


import numpy
from numpy import asarray, tan, exp, ones, squeeze, sign, \
     all, log, sqrt, pi, shape, array, minimum, where
from numpy import random

#__all__ = ['anneal']

_double_min = numpy.finfo(float).min
_double_max = numpy.finfo(float).max
class base_schedule(object):
    def __init__(self):
        self.dwell = 20
        self.learn_rate = 0.5
        self.lower = -10
        self.upper = 10
        self.Ninit = 50
        self.accepted = 0
        self.tests = 0
        self.feval = 0
        self.k = 0
        self.T = None

    def init(self, **options):
        self.__dict__.update(options)
        self.lower = asarray(self.lower)
        self.lower = where(self.lower == numpy.NINF, -_double_max, self.lower)
        self.upper = asarray(self.upper)
        self.upper = where(self.upper == numpy.PINF, _double_max, self.upper)
        self.k = 0
        self.accepted = 0
        self.feval = 0
        self.tests = 0

    def getstart_temp(self, best_state):
        """ Find a matching starting temperature and starting parameters vector
        i.e. find x0 such that func(x0) = T0.

        Parameters
        ----------
        best_state : _state
            A _state object to store the function value and x0 found.

        returns
        -------
        x0 : array
            The starting parameters vector.
        """

        assert(not self.dims is None)
        lrange = self.lower
        urange = self.upper
        fmax = _double_min
        fmin = _double_max
        for _ in range(self.Ninit):
            x0 = random.uniform(size=self.dims)*(urange-lrange) + lrange
            fval = self.func(x0, *self.args)
            self.feval += 1
            if fval > fmax:
                fmax = fval
            if fval < fmin:
                fmin = fval
                best_state.cost = fval
                best_state.x = array(x0)

        self.T0 = (fmax-fmin)*1.5
        return best_state.x

    def accept_test(self, dE):
        T = self.T
        self.tests += 1
        if dE < 0:
            self.accepted += 1
            return 1
        p = exp(-dE*1.0/self.boltzmann/T)
        if (p > random.uniform(0.0, 1.0)):
            self.accepted += 1
            return 1
        return 0

    def update_guess(self, x0):
        return np.squeeze(np.random.uniform(0.,1.,size=self.dims))*(self.upper-self.lower)+self.lower

    def update_temp(self, x0):
        pass


#  A schedule due to Lester Ingber modified to use bounds - OK
class fast_sa(base_schedule):
    def init(self, **options):
        self.__dict__.update(options)
        if self.m is None:
            self.m = 1.0
        if self.n is None:
            self.n = 1.0
        self.c = self.m * exp(-self.n * self.quench)

    def update_guess(self, x0):
        x0 = asarray(x0)
        u = squeeze(random.uniform(0.0, 1.0, size=self.dims))
        T = self.T
        xc = (sign(u-0.5)*T*((1+1.0/T)**abs(2*u-1)-1.0)+1.0)/2.0
        xnew = xc*(self.upper - self.lower)+self.lower
        return xnew
#        y = sign(u-0.5)*T*((1+1.0/T)**abs(2*u-1)-1.0)
#        xc = y*(self.upper - self.lower)
#        xnew = x0 + xc
#        return xnew

    def update_temp(self):
        self.T = self.T0*exp(-self.c * self.k**(self.quench))
        self.k += 1
        return

class cauchy_sa(base_schedule):     #modified to use bounds - not good
    def update_guess(self, x0):
        x0 = asarray(x0)
        numbers = squeeze(random.uniform(-pi/4, pi/4, size=self.dims))
        xc = (1.+(self.learn_rate * self.T * tan(numbers))%1.)
        xnew = xc*(self.upper - self.lower)+self.lower
        return xnew
#        numbers = squeeze(random.uniform(-pi/2, pi/2, size=self.dims))
#        xc = self.learn_rate * self.T * tan(numbers)
#        xnew = x0 + xc
#        return xnew

    def update_temp(self):
        self.T = self.T0/(1+self.k)
        self.k += 1
        return

class boltzmann_sa(base_schedule):
#    def update_guess(self, x0):
#        std = minimum(sqrt(self.T)*ones(self.dims), (self.upper-self.lower)/3.0/self.learn_rate)
#        x0 = asarray(x0)
#        xc = squeeze(random.normal(0, 1.0, size=self.dims))
#
#        xnew = x0 + xc*std*self.learn_rate
#        return xnew

    def update_temp(self):
        self.k += 1
        self.T = self.T0 / log(self.k+1.0)
        return

class log_sa(base_schedule):        #OK

    def init(self,**options):
        self.__dict__.update(options)
        
#    def update_guess(self,x0):
#        return np.squeeze(np.random.uniform(0.,1.,size=self.dims))*(self.upper-self.lower)+self.lower
        
    def update_temp(self):
        self.k += 1
        self.T = self.T0*self.slope**self.k
        
class _state(object):
    def __init__(self):
        self.x = None
        self.cost = None

# TODO:
#     allow for general annealing temperature profile
#     in that case use update given by alpha and omega and
#     variation of all previous updates and temperature?

# Simulated annealing

def anneal(func, x0, args=(), schedule='fast', full_output=0,
           T0=None, Tf=1e-12, maxeval=None, maxaccept=None, maxiter=400,
           boltzmann=1.0, learn_rate=0.5, feps=1e-6, quench=1.0, m=1.0, n=1.0,
           lower=-100, upper=100, dwell=50, slope=0.9,ranStart=True,dlg=None):
    """Minimize a function using simulated annealing.

    Schedule is a schedule class implementing the annealing schedule.
    Available ones are 'fast', 'cauchy', 'boltzmann'

    :param callable func: f(x, \*args)
        Function to be optimized.
    :param ndarray x0:
        Initial guess.
    :param tuple args: 
        Extra parameters to `func`.
    :param base_schedule schedule: 
        Annealing schedule to use (a class).
    :param bool full_output:
        Whether to return optional outputs.
    :param float T0: 
        Initial Temperature (estimated as 1.2 times the largest
        cost-function deviation over random points in the range).
    :param float Tf: 
        Final goal temperature.
    :param int maxeval: 
        Maximum function evaluations.
    :param int maxaccept:
        Maximum changes to accept.
    :param int maxiter: 
        Maximum cooling iterations.
    :param float learn_rate:
        Scale constant for adjusting guesses.
    :param float boltzmann: 
        Boltzmann constant in acceptance test
        (increase for less stringent test at each temperature).
    :param float feps:
        Stopping relative error tolerance for the function value in
        last four coolings.
    :param float quench,m,n:
        Parameters to alter fast_sa schedule.
    :param float/ndarray lower,upper: 
        Lower and upper bounds on `x`.
    :param int dwell:
        The number of times to search the space at each temperature.
    :param float slope: 
        Parameter for log schedule
    :param bool ranStart=True:
        False for fixed point start

    :returns: (xmin, Jmin, T, feval, iters, accept, retval) where

     * xmin (ndarray): Point giving smallest value found.
     * Jmin (float): Minimum value of function found.
     * T (float): Final temperature.
     * feval (int): Number of function evaluations.
     * iters (int): Number of cooling iterations.
     * accept (int): Number of tests accepted.
     * retval (int): Flag indicating stopping condition:

              *  0: Points no longer changing
              *  1: Cooled to final temperature
              *  2: Maximum function evaluations
              *  3: Maximum cooling iterations reached
              *  4: Maximum accepted query locations reached
              *  5: Final point not the minimum amongst encountered points

    *Notes*:
    Simulated annealing is a random algorithm which uses no derivative
    information from the function being optimized. In practice it has
    been more useful in discrete optimization than continuous
    optimization, as there are usually better algorithms for continuous
    optimization problems.

    Some experimentation by trying the difference temperature
    schedules and altering their parameters is likely required to
    obtain good performance.

    The randomness in the algorithm comes from random sampling in numpy.
    To obtain the same results you can call numpy.random.seed with the
    same seed immediately before calling scipy.optimize.anneal.

    We give a brief description of how the three temperature schedules
    generate new points and vary their temperature. Temperatures are
    only updated with iterations in the outer loop. The inner loop is
    over xrange(dwell), and new points are generated for
    every iteration in the inner loop. (Though whether the proposed
    new points are accepted is probabilistic.)

    For readability, let d denote the dimension of the inputs to func.
    Also, let x_old denote the previous state, and k denote the
    iteration number of the outer loop. All other variables not
    defined below are input variables to scipy.optimize.anneal itself.

    In the 'fast' schedule the updates are ::

        u ~ Uniform(0, 1, size=d)
        y = sgn(u - 0.5) * T * ((1+ 1/T)**abs(2u-1) -1.0)
        xc = y * (upper - lower)
        x_new = x_old + xc

        c = n * exp(-n * quench)
        T_new = T0 * exp(-c * k**quench)


    In the 'cauchy' schedule the updates are ::

        u ~ Uniform(-pi/2, pi/2, size=d)
        xc = learn_rate * T * tan(u)
        x_new = x_old + xc

        T_new = T0 / (1+k)

    In the 'boltzmann' schedule the updates are ::

        std = minimum( sqrt(T) * ones(d), (upper-lower) / (3*learn_rate) )
        y ~ Normal(0, std, size=d)
        x_new = x_old + learn_rate * y

        T_new = T0 / log(1+k)

    """
    x0 = asarray(x0)
    lower = asarray(lower)
    upper = asarray(upper)

    schedule = eval(schedule+'_sa()')
    #   initialize the schedule
    schedule.init(dims=shape(x0),func=func,args=args,boltzmann=boltzmann,T0=T0,
                  learn_rate=learn_rate, lower=lower, upper=upper,
                  m=m, n=n, quench=quench, dwell=dwell, slope=slope)

    current_state, last_state, best_state = _state(), _state(), _state()
    if T0 is None:
        x0 = schedule.getstart_temp(best_state)
    else:
        if ranStart:
            x0 = random.uniform(size=len(x0))*(upper-lower) + lower #comment to avoid random start
        best_state.x = None
        best_state.cost = numpy.Inf

    last_state.x = asarray(x0).copy()
    fval = func(x0,*args)
    schedule.feval += 1
    last_state.cost = fval
    if last_state.cost < best_state.cost:
        best_state.cost = fval
        best_state.x = asarray(x0).copy()
    schedule.T = schedule.T0
    fqueue = [100, 300, 500, 700]
    iters = 1
    keepGoing = True
    bestn = 0
    while keepGoing:
        retval = 0
        for n in xrange(dwell):
            current_state.x = schedule.update_guess(last_state.x)
            current_state.cost = func(current_state.x,*args)
            schedule.feval += 1

            dE = current_state.cost - last_state.cost
            if schedule.accept_test(dE):
                last_state.x = current_state.x.copy()
                last_state.cost = current_state.cost
                if last_state.cost < best_state.cost:
                    best_state.x = last_state.x.copy()
                    best_state.cost = last_state.cost
                    bestn = n
        if dlg:
            GoOn = dlg.Update(min(100.,best_state.cost*100),
                newmsg='%s%8.5f, %s%d\n%s%8.4f%s'%('Temperature =',schedule.T, \
                    'Best trial:',bestn,  \
                    'MC/SA Residual:',best_state.cost*100,'%', \
                    ))[0]
            if not GoOn:
                best_state.x = last_state.x.copy()
                best_state.cost = last_state.cost
                retval = 5
        schedule.update_temp()
        iters += 1
        # Stopping conditions
        # 0) last saved values of f from each cooling step
        #     are all very similar (effectively cooled)
        # 1) Tf is set and we are below it
        # 2) maxeval is set and we are past it
        # 3) maxiter is set and we are past it
        # 4) maxaccept is set and we are past it
        # 5) user canceled run via progress bar

        fqueue.append(squeeze(last_state.cost))
        fqueue.pop(0)
        af = asarray(fqueue)*1.0
        if retval == 5:
            print ' User terminated run; incomplete MC/SA'
            keepGoing = False
            break
        if all(abs((af-af[0])/af[0]) < feps):
            retval = 0
            if abs(af[-1]-best_state.cost) > feps*10:
                retval = 5
#                print "Warning: Cooled to %f at %s but this is not" \
#                      % (squeeze(last_state.cost), str(squeeze(last_state.x))) \
#                      + " the smallest point found."
            break
        if (Tf is not None) and (schedule.T < Tf):
            retval = 1
            break
        if (maxeval is not None) and (schedule.feval > maxeval):
            retval = 2
            break
        if (iters > maxiter):
            print "Warning: Maximum number of iterations exceeded."
            retval = 3
            break
        if (maxaccept is not None) and (schedule.accepted > maxaccept):
            retval = 4
            break

    if full_output:
        return best_state.x, best_state.cost, schedule.T, \
               schedule.feval, iters, schedule.accepted, retval
    else:
        return best_state.x, retval

def worker(iCyc,data,RBdata,reflType,reflData,covData,out_q):
    outlist = []
    for n in range(iCyc):
        result = mcsaSearch(data,RBdata,reflType,reflData,covData,None)
        outlist.append(result[0])
        print ' MC/SA residual: %.3f%% structure factor time: %.3f'%(100*result[0][2],result[1])
    out_q.put(outlist)

def MPmcsaSearch(nCyc,data,RBdata,reflType,reflData,covData):
    import multiprocessing as mp
    
    nprocs = mp.cpu_count()
    out_q = mp.Queue()
    procs = []
    iCyc = np.zeros(nprocs)
    for i in range(nCyc):
        iCyc[i%nprocs] += 1
    for i in range(nprocs):
        p = mp.Process(target=worker,args=(int(iCyc[i]),data,RBdata,reflType,reflData,covData,out_q))
        procs.append(p)
        p.start()
    resultlist = []
    for i in range(nprocs):
        resultlist += out_q.get()
    for p in procs:
        p.join()
    return resultlist

def mcsaSearch(data,RBdata,reflType,reflData,covData,pgbar):
    '''default doc string
    
    :param type name: description
    
    :returns: type name: description
    
    '''
    
    twopi = 2.0*np.pi
    global tsum
    tsum = 0.
    
    def getMDparms(item,pfx,parmDict,varyList):
        parmDict[pfx+'MDaxis'] = item['axis']
        parmDict[pfx+'MDval'] = item['Coef'][0]
        if item['Coef'][1]:
            varyList += [pfx+'MDval',]
            limits = item['Coef'][2]
            lower.append(limits[0])
            upper.append(limits[1])
                        
    def getAtomparms(item,pfx,aTypes,SGData,parmDict,varyList):
        parmDict[pfx+'Atype'] = item['atType']
        aTypes |= set([item['atType'],]) 
        pstr = ['Ax','Ay','Az']
        XYZ = [0,0,0]
        for i in range(3):
            name = pfx+pstr[i]
            parmDict[name] = item['Pos'][0][i]
            XYZ[i] = parmDict[name]
            if item['Pos'][1][i]:
                varyList += [name,]
                limits = item['Pos'][2][i]
                lower.append(limits[0])
                upper.append(limits[1])
        parmDict[pfx+'Amul'] = len(G2spc.GenAtom(XYZ,SGData))
            
    def getRBparms(item,mfx,aTypes,RBdata,SGData,atNo,parmDict,varyList):
        parmDict[mfx+'MolCent'] = item['MolCent']
        parmDict[mfx+'RBId'] = item['RBId']
        pstr = ['Px','Py','Pz']
        ostr = ['Qa','Qi','Qj','Qk']    #angle,vector not quaternion
        for i in range(3):
            name = mfx+pstr[i]
            parmDict[name] = item['Pos'][0][i]
            if item['Pos'][1][i]:
                varyList += [name,]
                limits = item['Pos'][2][i]
                lower.append(limits[0])
                upper.append(limits[1])
        AV = item['Ori'][0]
        A = AV[0]
        V = AV[1:]
        for i in range(4):
            name = mfx+ostr[i]
            if i:
                parmDict[name] = V[i-1]
            else:
                parmDict[name] = A
            if item['Ovar'] == 'AV':
                varyList += [name,]
                limits = item['Ori'][2][i]
                lower.append(limits[0])
                upper.append(limits[1])
            elif item['Ovar'] == 'A' and not i:
                varyList += [name,]
                limits = item['Ori'][2][i]
                lower.append(limits[0])
                upper.append(limits[1])
        if 'Tor' in item:      #'Tor' not there for 'Vector' RBs
            for i in range(len(item['Tor'][0])):
                name = mfx+'Tor'+str(i)
                parmDict[name] = item['Tor'][0][i]
                if item['Tor'][1][i]:
                    varyList += [name,]
                    limits = item['Tor'][2][i]
                    lower.append(limits[0])
                    upper.append(limits[1])
        atypes = RBdata[item['Type']][item['RBId']]['rbTypes']
        aTypes |= set(atypes)
        atNo += len(atypes)
        return atNo
                
    def GetAtomM(Xdata,SGData):
        Mdata = []
        for xyz in Xdata:
            Mdata.append(float(len(G2spc.GenAtom(xyz,SGData))))
        return np.array(Mdata)
        
    def GetAtomTX(RBdata,parmDict):
        'Needs a doc string'
        Bmat = parmDict['Bmat']
        atNo = parmDict['atNo']
        nfixAt = parmDict['nfixAt']
        Tdata = atNo*[' ',]
        Xdata = np.zeros((3,atNo))
        keys = {':Atype':Tdata,':Ax':Xdata[0],':Ay':Xdata[1],':Az':Xdata[2]}
        for iatm in range(nfixAt):
            for key in keys:
                parm = ':'+str(iatm)+key
                if parm in parmDict:
                    if key == ':Atype':
                        keys[key][iatm] = aTypes.index(parmDict[parm])
                    else:
                        keys[key][iatm] = parmDict[parm]
        iatm = nfixAt
        for iObj in range(parmDict['nObj']):
            pfx = str(iObj)+':'
            if parmDict[pfx+'Type'] in ['Vector','Residue']:
                if parmDict[pfx+'Type'] == 'Vector':
                    RBRes = RBdata['Vector'][parmDict[pfx+'RBId']]
                    vecs = RBRes['rbVect']
                    mags = RBRes['VectMag']
                    Cart = np.zeros_like(vecs[0])
                    for vec,mag in zip(vecs,mags):
                        Cart += vec*mag
                elif parmDict[pfx+'Type'] == 'Residue':
                    RBRes = RBdata['Residue'][parmDict[pfx+'RBId']]
                    Cart = np.array(RBRes['rbXYZ'])
                    for itor,seq in enumerate(RBRes['rbSeq']):
                        QuatA = AVdeg2Q(parmDict[pfx+'Tor'+str(itor)],Cart[seq[0]]-Cart[seq[1]])
                        for ride in seq[3]:
                            Cart[ride] = prodQVQ(QuatA,Cart[ride]-Cart[seq[1]])+Cart[seq[1]]
                if parmDict[pfx+'MolCent'][1]:
                    Cart -= parmDict[pfx+'MolCent'][0]
                Qori = AVdeg2Q(parmDict[pfx+'Qa'],[parmDict[pfx+'Qi'],parmDict[pfx+'Qj'],parmDict[pfx+'Qk']])
                Pos = np.array([parmDict[pfx+'Px'],parmDict[pfx+'Py'],parmDict[pfx+'Pz']])
                for i,x in enumerate(Cart):
                    X = np.inner(Bmat,prodQVQ(Qori,x))+Pos
                    for j in range(3):
                        Xdata[j][iatm] = X[j]
                    Tdata[iatm] = aTypes.index(RBRes['rbTypes'][i])
                    iatm += 1
            elif parmDict[pfx+'Type'] == 'Atom':
                atNo = parmDict[pfx+'atNo']
                for key in keys:
                    parm = pfx+key[1:]              #remove extra ':'
                    if parm in parmDict:
                        if key == ':Atype':
                            keys[key][atNo] = aTypes.index(parmDict[parm])
                        else:
                            keys[key][atNo] = parmDict[parm]
                iatm += 1
            else:
                continue        #skips March Dollase
        return Tdata,Xdata.T
        
    def getAllTX(Tdata,Mdata,Xdata,SGM,SGT):
        allX = np.inner(Xdata,SGM)+SGT
        allT = np.repeat(Tdata,allX.shape[1])
        allM = np.repeat(Mdata,allX.shape[1])
        allX = np.reshape(allX,(-1,3))
        return allT,allM,allX

    def getAllX(Xdata,SGM,SGT):
        allX = np.inner(Xdata,SGM)+SGT
        allX = np.reshape(allX,(-1,3))
        return allX
        
    def normQuaternions(RBdata,parmDict,varyList,values):
        for iObj in range(parmDict['nObj']):
            pfx = str(iObj)+':'
            if parmDict[pfx+'Type'] in ['Vector','Residue']:
                Qori = AVdeg2Q(parmDict[pfx+'Qa'],[parmDict[pfx+'Qi'],parmDict[pfx+'Qj'],parmDict[pfx+'Qk']])
                A,V = Q2AVdeg(Qori)
                for i,name in enumerate(['Qa','Qi','Qj','Qk']):
                    if i:
                        parmDict[pfx+name] = V[i-1]
                    else:
                        parmDict[pfx+name] = A
        
    def mcsaCalc(values,refList,rcov,ifInv,allFF,RBdata,varyList,parmDict):
        ''' Compute structure factors for all h,k,l for phase
        input:
            refList: [ref] where each ref = h,k,l,m,d,...
            rcov:   array[nref,nref] covariance terms between Fo^2 values
            ifInv:  bool True if centrosymmetric
            allFF: array[nref,natoms] each value is mult*FF(H)/max(mult)
            RBdata: [dict] rigid body dictionary
            varyList: [list] names of varied parameters in MC/SA (not used here)           
            ParmDict: [dict] problem parameters
        puts result F^2 in each ref[5] in refList
        returns:
            delt-F*rcov*delt-F/sum(Fo^2)^2
        '''       
        global tsum
        parmDict.update(dict(zip(varyList,values)))             #update parameter tables
        t0 = time.time()
        Xdata = GetAtomTX(RBdata,parmDict)[1]                   #get new atom coords from RB
        tsum += (time.time()-t0)
        allX = getAllX(Xdata,SGM,SGT)                           #fill unit cell - dups. OK
        MDval = parmDict['0:MDval']                             #get March-Dollase coeff
        MDaxis = parmDict['0:MDaxis']
        Gmat = parmDict['Gmat']
        HX2pi = 2.*np.pi*np.inner(allX,refList[:3].T)           #form 2piHX for every H,X pair
        Aterm = refList[3]*np.sum(allFF*np.cos(HX2pi),axis=0)**2    #compute real part for all H
        if ifInv:
            refList[5] = Aterm
        else:
            Bterm = refList[3]*np.sum(allFF*np.sin(HX2pi),axis=0)**2    #imaginary part for all H
            refList[5] = Aterm+Bterm
        sumFcsq = np.sum(refList[5])
        scale = parmDict['sumFosq']/sumFcsq
        refList[5] *= scale
        refList[6] = refList[4]-refList[5]
        M = np.inner(refList[6],np.inner(rcov,refList[6]))
#        print M,parmDict['sumFosq'],np.sum(refList[6]**2),np.sum(refList[4]**2)
#        print np.sum(refList[6]**2)/np.sum(refList[4]**2)
        return M/np.sum(refList[4]**2)

    sq8ln2 = np.sqrt(8*np.log(2))
    sq2pi = np.sqrt(2*np.pi)
    sq4pi = np.sqrt(4*np.pi)
    generalData = data['General']
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    Gmat = G2lat.cell2Gmat(generalData['Cell'][1:7])[0]
    SGData = generalData['SGData']
    SGM = [SGData['SGOps'][i][0] for i in range(len(SGData['SGOps']))]
    SGT = [SGData['SGOps'][i][1] for i in range(len(SGData['SGOps']))]
    fixAtoms = data['Atoms']                       #if any
    cx,ct,cs = generalData['AtomPtrs'][:3]
    aTypes = set([])
    parmDict = {'Bmat':Bmat,'Gmat':Gmat}
    varyList = []
    atNo = 0
    for atm in fixAtoms:
        pfx = ':'+str(atNo)+':'
        parmDict[pfx+'Atype'] = atm[ct]
        aTypes |= set([atm[ct],])
        pstr = ['Ax','Ay','Az']
        parmDict[pfx+'Amul'] = atm[cs+1]
        for i in range(3):
            name = pfx+pstr[i]
            parmDict[name] = atm[cx+i]
        atNo += 1
    parmDict['nfixAt'] = len(fixAtoms)        
    MCSA = generalData['MCSA controls']
    reflName = MCSA['Data source']
    phaseName = generalData['Name']
    MCSAObjs = data['MCSA']['Models']               #list of MCSA models
    upper = []
    lower = []
    for i,item in enumerate(MCSAObjs):
        mfx = str(i)+':'
        parmDict[mfx+'Type'] = item['Type']
        if item['Type'] == 'MD':
            getMDparms(item,mfx,parmDict,varyList)
        elif item['Type'] == 'Atom':
            getAtomparms(item,mfx,aTypes,SGData,parmDict,varyList)
            parmDict[mfx+'atNo'] = atNo
            atNo += 1
        elif item['Type'] in ['Residue','Vector']:
            atNo = getRBparms(item,mfx,aTypes,RBdata,SGData,atNo,parmDict,varyList)
    parmDict['atNo'] = atNo                 #total no. of atoms
    parmDict['nObj'] = len(MCSAObjs)
    aTypes = list(aTypes)
    Tdata,Xdata = GetAtomTX(RBdata,parmDict)
    Mdata = GetAtomM(Xdata,SGData)
    allT,allM = getAllTX(Tdata,Mdata,Xdata,SGM,SGT)[:2]
    FFtables = G2el.GetFFtable(aTypes)
    refs = []
    allFF = []
    sumFosq = 0
    if 'PWDR' in reflName:
        for ref in reflData:
            h,k,l,m,d,pos,sig,gam,f = ref[:9]
            if d >= MCSA['dmin']:
                sig = G2pwd.getgamFW(sig,gam)/sq8ln2        #--> sig from FWHM
                SQ = 0.25/d**2
                allFF.append(allM*[G2el.getFFvalues(FFtables,SQ,True)[i] for i in allT]/np.max(allM))
                refs.append([h,k,l,m,f*m,pos,sig])
                sumFosq += f*m
        nRef = len(refs)
        rcov = np.zeros((nRef,nRef))
        for iref,refI in enumerate(refs):
            rcov[iref][iref] = 1./(sq4pi*refI[6])
            for jref,refJ in enumerate(refs[:iref]):
                t1 = refI[6]**2+refJ[6]**2
                t2 = (refJ[5]-refI[5])**2/(2.*t1)
                if t2 > 10.:
                    rcov[iref][jref] = 0.
                else:
                    rcov[iref][jref] = 1./(sq2pi*np.sqrt(t1)*np.exp(t2))
        rcov += (rcov.T-np.diagflat(np.diagonal(rcov)))
        Rdiag = np.sqrt(np.diag(rcov))
        Rnorm = np.outer(Rdiag,Rdiag)
        rcov /= Rnorm
    elif 'Pawley' in reflName:  #need a bail out if Pawley cov matrix doesn't exist.
        vList = covData['varyList']
        for iref,refI in enumerate(reflData):
            h,k,l,m,d,v,f,s = refI
            if d >= MCSA['dmin'] and v:       #skip unrefined ones
                SQ = 0.25/d**2
                allFF.append(allM*[G2el.getFFvalues(FFtables,SQ,True)[i] for i in allT]/np.max(allM))
                refs.append([h,k,l,m,f*m,iref,0.])
                sumFosq += f*m
        nRef = len(refs)
        pfx = str(data['pId'])+'::PWLref:'
        if covData['freshCOV'] and generalData['doPawley']:
            covMatrix = covData['covMatrix']
            rcov = np.zeros((nRef,nRef))        
            for iref,refI in enumerate(refs):
                I = refI[5]
                nameI = pfx+str(I)
                if nameI in vList:
                    Iindx = vList.index(nameI)
                    rcov[iref][iref] = covMatrix[Iindx][Iindx]
                    for jref,refJ in enumerate(refs[:iref]):
                        J = refJ[5]
                        nameJ = pfx+str(J)
                        try:
                            rcov[iref][jref] = covMatrix[vList.index(nameI)][vList.index(nameJ)]
                        except ValueError:
                            rcov[iref][jref] = rcov[iref][jref-1]
                else:
                    rcov[iref] = rcov[iref-1]
                    rcov[iref][iref] = rcov[iref-1][iref-1]
            rcov += (rcov.T-np.diagflat(np.diagonal(rcov)))
            Rdiag = np.sqrt(np.diag(rcov))
            Rnorm = np.outer(Rdiag,Rdiag)
            rcov /= Rnorm
            MCSA['rcov'] = rcov
            covData['freshCOV'] = False
        else:
            rcov = MCSA['rcov']
    elif 'HKLF' in reflName:
        for ref in reflData:
            [h,k,l,m,d],f = ref[:5],ref[6]
            if d >= MCSA['dmin']:
                SQ = 0.25/d**2
                allFF.append(allM*[G2el.getFFvalues(FFtables,SQ,True)[i] for i in allT]/np.max(allM))
                refs.append([h,k,l,m,f*m,0.,0.])
                sumFosq += f*m
        nRef = len(refs)
        rcov = np.identity(len(refs))
    allFF = np.array(allFF).T
    refs = np.array(refs).T
    print ' Minimum d-spacing used: %.2f No. reflections used: %d'%(MCSA['dmin'],nRef)
    print ' Number of parameters varied: %d'%(len(varyList))
    parmDict['sumFosq'] = sumFosq
    x0 = [parmDict[val] for val in varyList]
    ifInv = SGData['SGInv']
    results = anneal(mcsaCalc,x0,args=(refs,rcov,ifInv,allFF,RBdata,varyList,parmDict),
        schedule=MCSA['Algorithm'], full_output=True,
        T0=MCSA['Annealing'][0], Tf=MCSA['Annealing'][1],dwell=MCSA['Annealing'][2],
        boltzmann=MCSA['boltzmann'], learn_rate=0.5,  
        quench=MCSA['fast parms'][0], m=MCSA['fast parms'][1], n=MCSA['fast parms'][2],
        lower=lower, upper=upper, slope=MCSA['log slope'],ranStart=MCSA.get('ranStart',True),dlg=pgbar)
    M = mcsaCalc(results[0],refs,rcov,ifInv,allFF,RBdata,varyList,parmDict)
#    for ref in refs.T:
#        print ' %4d %4d %4d %10.3f %10.3f %10.3f'%(int(ref[0]),int(ref[1]),int(ref[2]),ref[4],ref[5],ref[6])
#    print np.sqrt((np.sum(refs[6]**2)/np.sum(refs[4]**2)))
    Result = [False,False,results[1],results[2],]+list(results[0])
    Result.append(varyList)
    return Result,tsum

        
################################################################################
##### Quaternion stuff
################################################################################

def prodQQ(QA,QB):
    ''' Grassman quaternion product
        QA,QB quaternions; q=r+ai+bj+ck
    '''
    D = np.zeros(4)
    D[0] = QA[0]*QB[0]-QA[1]*QB[1]-QA[2]*QB[2]-QA[3]*QB[3]
    D[1] = QA[0]*QB[1]+QA[1]*QB[0]+QA[2]*QB[3]-QA[3]*QB[2]
    D[2] = QA[0]*QB[2]-QA[1]*QB[3]+QA[2]*QB[0]+QA[3]*QB[1]
    D[3] = QA[0]*QB[3]+QA[1]*QB[2]-QA[2]*QB[1]+QA[3]*QB[0]
    
#    D[0] = QA[0]*QB[0]-np.dot(QA[1:],QB[1:])
#    D[1:] = QA[0]*QB[1:]+QB[0]*QA[1:]+np.cross(QA[1:],QB[1:])
    
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
    """
    compute the quaternion vector rotation qvq-1 = v'
    q=r+ai+bj+ck
    """
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
    QN = normQ(Q)
    aa = QN[0]**2
    ab = QN[0]*QN[1]
    ac = QN[0]*QN[2]
    ad = QN[0]*QN[3]
    bb = QN[1]**2
    bc = QN[1]*QN[2]
    bd = QN[1]*QN[3]
    cc = QN[2]**2
    cd = QN[2]*QN[3]
    dd = QN[3]**2
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
        p = A/2.
        Q[0] = np.cos(p)
        Q[1:4] = V*np.sin(p)
    else:
        Q[3] = 1.
    return Q
    
def AVdeg2Q(A,V):
    ''' convert angle (degrees) & vector to quaternion
        q=r+ai+bj+ck
    '''
    Q = np.zeros(4)
    d = np.sqrt(np.sum(np.array(V)**2))
    if d:
        V /= d
        p = A/2.
        Q[0] = cosd(p)
        Q[1:4] = V*sind(p)
    else:
        Q[3] = 1.
    return Q
    
def Q2AVdeg(Q):
    ''' convert quaternion to angle (degrees 0-360) & normalized vector
        q=r+ai+bj+ck
    '''
    A = 2.*acosd(Q[0])
    V = np.array(Q[1:])
    V /= sind(A/2.)
    return A,V
    
def Q2AV(Q):
    ''' convert quaternion to angle (radians 0-2pi) & normalized vector
        q=r+ai+bj+ck
    '''
    A = 2.*np.arccos(Q[0])
    V = np.array(Q[1:])
    V /= np.sin(A/2.)
    return A,V
    
def makeQuat(A,B,C):
    ''' Make quaternion from rotation of A vector to B vector about C axis

        :param np.array A,B,C: Cartesian 3-vectors
        :returns: quaternion & rotation angle in radians q=r+ai+bj+ck
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
    
if __name__ == "__main__":
    from numpy import cos
    # minimum expected at ~-0.195
    func = lambda x: cos(14.5*x-0.3) + (x+0.2)*x
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='cauchy')
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='fast')
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='boltzmann')

    # minimum expected at ~[-0.195, -0.1]
    func = lambda x: cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='cauchy')
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='fast')
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='boltzmann')
