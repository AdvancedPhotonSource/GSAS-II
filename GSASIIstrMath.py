# -*- coding: utf-8 -*-
'''
*GSASIIstrMath - structure math routines*
-----------------------------------------
'''
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
from __future__ import division, print_function
import time
import copy
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import scipy.stats as st
import multiprocessing as mp
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIpwd as G2pwd
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIobj as G2obj
import GSASIImpsubs as G2mp
#G2mp.InitMP(False)  # This disables multiprocessing 

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
try:  # fails on doc build
    ateln2 = 8.0*np.log(2.0)
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
except TypeError:
    pass
nxs = np.newaxis

################################################################################
##### Rigid Body Models
################################################################################
        
def ApplyRBModels(parmDict,Phases,rigidbodyDict,Update=False):
    ''' Takes RB info from RBModels in Phase and RB data in rigidbodyDict along with
    current RB values in parmDict & modifies atom contents (xyz & Uij) of parmDict
    '''
    atxIds = ['Ax:','Ay:','Az:']
    atuIds = ['AU11:','AU22:','AU33:','AU12:','AU13:','AU23:']
    RBIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})  #these are lists of rbIds
    if not RBIds['Vector'] and not RBIds['Residue']:
        return
    VRBIds = RBIds['Vector']
    RRBIds = RBIds['Residue']
    if Update:
        RBData = rigidbodyDict
    else:
        RBData = copy.deepcopy(rigidbodyDict)     # don't mess with original!
    if RBIds['Vector']:                       # first update the vector magnitudes
        VRBData = RBData['Vector']
        for i,rbId in enumerate(VRBIds):
            if VRBData[rbId]['useCount']:
                for j in range(len(VRBData[rbId]['VectMag'])):
                    name = '::RBV;'+str(j)+':'+str(i)
                    VRBData[rbId]['VectMag'][j] = parmDict[name]
    for phase in Phases:
        Phase = Phases[phase]
        General = Phase['General']
        cx,ct,cs,cia = General['AtomPtrs']
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        AtLookup = G2mth.FillAtomLookUp(Phase['Atoms'],cia+8)
        pfx = str(Phase['pId'])+'::'
        if Update:
            RBModels = Phase['RBModels']
        else:
            RBModels =  copy.deepcopy(Phase['RBModels']) # again don't mess with original!
        for irb,RBObj in enumerate(RBModels.get('Vector',[])):
            jrb = VRBIds.index(RBObj['RBId'])
            rbsx = str(irb)+':'+str(jrb)
            for i,px in enumerate(['RBVPx:','RBVPy:','RBVPz:']):
                RBObj['Orig'][0][i] = parmDict[pfx+px+rbsx]
            for i,po in enumerate(['RBVOa:','RBVOi:','RBVOj:','RBVOk:']):
                RBObj['Orient'][0][i] = parmDict[pfx+po+rbsx]
            RBObj['Orient'][0] = G2mth.normQ(RBObj['Orient'][0])
            TLS = RBObj['ThermalMotion']
            if 'T' in TLS[0]:
                for i,pt in enumerate(['RBVT11:','RBVT22:','RBVT33:','RBVT12:','RBVT13:','RBVT23:']):
                    TLS[1][i] = parmDict[pfx+pt+rbsx]
            if 'L' in TLS[0]:
                for i,pt in enumerate(['RBVL11:','RBVL22:','RBVL33:','RBVL12:','RBVL13:','RBVL23:']):
                    TLS[1][i+6] = parmDict[pfx+pt+rbsx]
            if 'S' in TLS[0]:
                for i,pt in enumerate(['RBVS12:','RBVS13:','RBVS21:','RBVS23:','RBVS31:','RBVS32:','RBVSAA:','RBVSBB:']):
                    TLS[1][i+12] = parmDict[pfx+pt+rbsx]
            if 'U' in TLS[0]:
                TLS[1][0] = parmDict[pfx+'RBVU:'+rbsx]
            XYZ,Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Vector')
            UIJ = G2mth.UpdateRBUIJ(Bmat,Cart,RBObj)
            for i,x in enumerate(XYZ):
                atId = RBObj['Ids'][i]
                for j in [0,1,2]:
                    parmDict[pfx+atxIds[j]+str(AtLookup[atId])] = x[j]
                if UIJ[i][0] == 'A':
                    for j in range(6):
                        parmDict[pfx+atuIds[j]+str(AtLookup[atId])] = UIJ[i][j+2]
                elif UIJ[i][0] == 'I':
                    parmDict[pfx+'AUiso:'+str(AtLookup[atId])] = UIJ[i][1]
            
        for irb,RBObj in enumerate(RBModels.get('Residue',[])):
            jrb = RRBIds.index(RBObj['RBId'])
            rbsx = str(irb)+':'+str(jrb)
            for i,px in enumerate(['RBRPx:','RBRPy:','RBRPz:']):
                RBObj['Orig'][0][i] = parmDict[pfx+px+rbsx]
            for i,po in enumerate(['RBROa:','RBROi:','RBROj:','RBROk:']):
                RBObj['Orient'][0][i] = parmDict[pfx+po+rbsx]                
            RBObj['Orient'][0] = G2mth.normQ(RBObj['Orient'][0])
            TLS = RBObj['ThermalMotion']
            if 'T' in TLS[0]:
                for i,pt in enumerate(['RBRT11:','RBRT22:','RBRT33:','RBRT12:','RBRT13:','RBRT23:']):
                    RBObj['ThermalMotion'][1][i] = parmDict[pfx+pt+rbsx]
            if 'L' in TLS[0]:
                for i,pt in enumerate(['RBRL11:','RBRL22:','RBRL33:','RBRL12:','RBRL13:','RBRL23:']):
                    RBObj['ThermalMotion'][1][i+6] = parmDict[pfx+pt+rbsx]
            if 'S' in TLS[0]:
                for i,pt in enumerate(['RBRS12:','RBRS13:','RBRS21:','RBRS23:','RBRS31:','RBRS32:','RBRSAA:','RBRSBB:']):
                    RBObj['ThermalMotion'][1][i+12] = parmDict[pfx+pt+rbsx]
            if 'U' in TLS[0]:
                RBObj['ThermalMotion'][1][0] = parmDict[pfx+'RBRU:'+rbsx]
            for itors,tors in enumerate(RBObj['Torsions']):
                tors[0] = parmDict[pfx+'RBRTr;'+str(itors)+':'+rbsx]
            XYZ,Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')
            UIJ = G2mth.UpdateRBUIJ(Bmat,Cart,RBObj)
            for i,x in enumerate(XYZ):
                atId = RBObj['Ids'][i]
                for j in [0,1,2]:
                    parmDict[pfx+atxIds[j]+str(AtLookup[atId])] = x[j]
                if UIJ[i][0] == 'A':
                    for j in range(6):
                        parmDict[pfx+atuIds[j]+str(AtLookup[atId])] = UIJ[i][j+2]
                elif UIJ[i][0] == 'I':
                    parmDict[pfx+'AUiso:'+str(AtLookup[atId])] = UIJ[i][1]
                    
def ApplyRBModelDervs(dFdvDict,parmDict,rigidbodyDict,Phase):
    'Computes rigid body derivatives'
    atxIds = ['dAx:','dAy:','dAz:']
    atuIds = ['AU11:','AU22:','AU33:','AU12:','AU13:','AU23:']
    OIds = ['Oa:','Oi:','Oj:','Ok:']
    RBIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})  #these are lists of rbIds
    if not RBIds['Vector'] and not RBIds['Residue']:
        return
    VRBIds = RBIds['Vector']
    RRBIds = RBIds['Residue']
    RBData = rigidbodyDict
    for item in parmDict:
        if 'RB' in item:
            dFdvDict[item] = 0.        #NB: this is a vector which is no. refl. long & must be filled!
    General = Phase['General']
    cx,ct,cs,cia = General['AtomPtrs']
    cell = General['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(cell)
    rpd = np.pi/180.
    rpd2 = rpd**2
    g = nl.inv(np.inner(Bmat,Bmat))
    gvec = np.sqrt(np.array([g[0][0]**2,g[1][1]**2,g[2][2]**2,
        g[0][0]*g[1][1],g[0][0]*g[2][2],g[1][1]*g[2][2]]))
    AtLookup = G2mth.FillAtomLookUp(Phase['Atoms'],cia+8)
    pfx = str(Phase['pId'])+'::'
    RBModels =  Phase['RBModels']
    for irb,RBObj in enumerate(RBModels.get('Vector',[])):
        symAxis = RBObj.get('symAxis')
        VModel = RBData['Vector'][RBObj['RBId']]
        Q = RBObj['Orient'][0]
        jrb = VRBIds.index(RBObj['RBId'])
        rbsx = str(irb)+':'+str(jrb)
        dXdv = []
        for iv in range(len(VModel['VectMag'])):
            dCdv = []
            for vec in VModel['rbVect'][iv]:
                dCdv.append(G2mth.prodQVQ(Q,vec))
            dXdv.append(np.inner(Bmat,np.array(dCdv)).T)
        XYZ,Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Vector')
        for ia,atId in enumerate(RBObj['Ids']):
            atNum = AtLookup[atId]
            dx = 0.00001
            for iv in range(len(VModel['VectMag'])):
                for ix in [0,1,2]:
                    dFdvDict['::RBV;'+str(iv)+':'+str(jrb)] += dXdv[iv][ia][ix]*dFdvDict[pfx+atxIds[ix]+str(atNum)]
            for i,name in enumerate(['RBVPx:','RBVPy:','RBVPz:']):
                dFdvDict[pfx+name+rbsx] += dFdvDict[pfx+atxIds[i]+str(atNum)]
            for iv in range(4):
                Q[iv] -= dx
                XYZ1 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q),symAxis)
                Q[iv] += 2.*dx
                XYZ2 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q),symAxis)
                Q[iv] -= dx
                dXdO = (XYZ2[ia]-XYZ1[ia])/(2.*dx)
                for ix in [0,1,2]:
                    dFdvDict[pfx+'RBV'+OIds[iv]+rbsx] += dXdO[ix]*dFdvDict[pfx+atxIds[ix]+str(atNum)]
            X = G2mth.prodQVQ(Q,Cart[ia])
            dFdu = np.array([dFdvDict[pfx+Uid+str(AtLookup[atId])] for Uid in atuIds]).T/gvec
            dFdu = G2lat.U6toUij(dFdu.T)
            dFdu = np.tensordot(Amat,np.tensordot(Amat,dFdu,([1,0])),([0,1]))            
            dFdu = G2lat.UijtoU6(dFdu)
            atNum = AtLookup[atId]
            if 'T' in RBObj['ThermalMotion'][0]:
                for i,name in enumerate(['RBVT11:','RBVT22:','RBVT33:','RBVT12:','RBVT13:','RBVT23:']):
                    dFdvDict[pfx+name+rbsx] += dFdu[i]
            if 'L' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBVL11:'+rbsx] += rpd2*(dFdu[1]*X[2]**2+dFdu[2]*X[1]**2-dFdu[5]*X[1]*X[2])
                dFdvDict[pfx+'RBVL22:'+rbsx] += rpd2*(dFdu[0]*X[2]**2+dFdu[2]*X[0]**2-dFdu[4]*X[0]*X[2])
                dFdvDict[pfx+'RBVL33:'+rbsx] += rpd2*(dFdu[0]*X[1]**2+dFdu[1]*X[0]**2-dFdu[3]*X[0]*X[1])
                dFdvDict[pfx+'RBVL12:'+rbsx] += rpd2*(-dFdu[3]*X[2]**2-2.*dFdu[2]*X[0]*X[1]+
                    dFdu[4]*X[1]*X[2]+dFdu[5]*X[0]*X[2])
                dFdvDict[pfx+'RBVL13:'+rbsx] += rpd2*(-dFdu[4]*X[1]**2-2.*dFdu[1]*X[0]*X[2]+
                    dFdu[3]*X[1]*X[2]+dFdu[5]*X[0]*X[1])
                dFdvDict[pfx+'RBVL23:'+rbsx] += rpd2*(-dFdu[5]*X[0]**2-2.*dFdu[0]*X[1]*X[2]+
                    dFdu[3]*X[0]*X[2]+dFdu[4]*X[0]*X[1])
            if 'S' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBVS12:'+rbsx] += rpd*(dFdu[5]*X[1]-2.*dFdu[1]*X[2])
                dFdvDict[pfx+'RBVS13:'+rbsx] += rpd*(-dFdu[5]*X[2]+2.*dFdu[2]*X[1])
                dFdvDict[pfx+'RBVS21:'+rbsx] += rpd*(-dFdu[4]*X[0]+2.*dFdu[0]*X[2])
                dFdvDict[pfx+'RBVS23:'+rbsx] += rpd*(dFdu[4]*X[2]-2.*dFdu[2]*X[0])
                dFdvDict[pfx+'RBVS31:'+rbsx] += rpd*(dFdu[3]*X[0]-2.*dFdu[0]*X[1])
                dFdvDict[pfx+'RBVS32:'+rbsx] += rpd*(-dFdu[3]*X[1]+2.*dFdu[1]*X[0])
                dFdvDict[pfx+'RBVSAA:'+rbsx] += rpd*(dFdu[4]*X[1]-dFdu[3]*X[2])
                dFdvDict[pfx+'RBVSBB:'+rbsx] += rpd*(dFdu[5]*X[0]-dFdu[3]*X[2])
            if 'U' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBVU:'+rbsx] += dFdvDict[pfx+'AUiso:'+str(AtLookup[atId])]


    for irb,RBObj in enumerate(RBModels.get('Residue',[])):
        symAxis = RBObj.get('symAxis')
        Q = RBObj['Orient'][0]
        jrb = RRBIds.index(RBObj['RBId'])
        torData = RBData['Residue'][RBObj['RBId']]['rbSeq']
        rbsx = str(irb)+':'+str(jrb)
        XYZ,Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')
        for itors,tors in enumerate(RBObj['Torsions']):     #derivative error?
            tname = pfx+'RBRTr;'+str(itors)+':'+rbsx           
            orId,pvId = torData[itors][:2]
            pivotVec = Cart[orId]-Cart[pvId]
            QA = G2mth.AVdeg2Q(-0.001,pivotVec)
            QB = G2mth.AVdeg2Q(0.001,pivotVec)
            for ir in torData[itors][3]:
                atNum = AtLookup[RBObj['Ids'][ir]]
                rVec = Cart[ir]-Cart[pvId]
                dR = G2mth.prodQVQ(QB,rVec)-G2mth.prodQVQ(QA,rVec)
                dRdT = np.inner(Bmat,G2mth.prodQVQ(Q,dR))/.002
                for ix in [0,1,2]:
                    dFdvDict[tname] += dRdT[ix]*dFdvDict[pfx+atxIds[ix]+str(atNum)]
        for ia,atId in enumerate(RBObj['Ids']):
            atNum = AtLookup[atId]
            dx = 0.00001
            for i,name in enumerate(['RBRPx:','RBRPy:','RBRPz:']):
                dFdvDict[pfx+name+rbsx] += dFdvDict[pfx+atxIds[i]+str(atNum)]
            for iv in range(4):
                Q[iv] -= dx
                XYZ1 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q),symAxis)
                Q[iv] += 2.*dx
                XYZ2 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q),symAxis)
                Q[iv] -= dx
                dXdO = (XYZ2[ia]-XYZ1[ia])/(2.*dx)
                for ix in [0,1,2]:
                    dFdvDict[pfx+'RBR'+OIds[iv]+rbsx] += dXdO[ix]*dFdvDict[pfx+atxIds[ix]+str(atNum)]
            X = G2mth.prodQVQ(Q,Cart[ia])
            dFdu = np.array([dFdvDict[pfx+Uid+str(AtLookup[atId])] for Uid in atuIds]).T/gvec
            dFdu = G2lat.U6toUij(dFdu.T)
            dFdu = np.tensordot(Amat.T,np.tensordot(Amat,dFdu,([1,0])),([0,1]))
            dFdu = G2lat.UijtoU6(dFdu)
            atNum = AtLookup[atId]
            if 'T' in RBObj['ThermalMotion'][0]:
                for i,name in enumerate(['RBRT11:','RBRT22:','RBRT33:','RBRT12:','RBRT13:','RBRT23:']):
                    dFdvDict[pfx+name+rbsx] += dFdu[i]
            if 'L' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBRL11:'+rbsx] += rpd2*(dFdu[1]*X[2]**2+dFdu[2]*X[1]**2-dFdu[5]*X[1]*X[2])
                dFdvDict[pfx+'RBRL22:'+rbsx] += rpd2*(dFdu[0]*X[2]**2+dFdu[2]*X[0]**2-dFdu[4]*X[0]*X[2])
                dFdvDict[pfx+'RBRL33:'+rbsx] += rpd2*(dFdu[0]*X[1]**2+dFdu[1]*X[0]**2-dFdu[3]*X[0]*X[1])
                dFdvDict[pfx+'RBRL12:'+rbsx] += rpd2*(-dFdu[3]*X[2]**2-2.*dFdu[2]*X[0]*X[1]+
                    dFdu[4]*X[1]*X[2]+dFdu[5]*X[0]*X[2])
                dFdvDict[pfx+'RBRL13:'+rbsx] += rpd2*(dFdu[4]*X[1]**2-2.*dFdu[1]*X[0]*X[2]+
                    dFdu[3]*X[1]*X[2]+dFdu[5]*X[0]*X[1])
                dFdvDict[pfx+'RBRL23:'+rbsx] += rpd2*(dFdu[5]*X[0]**2-2.*dFdu[0]*X[1]*X[2]+
                    dFdu[3]*X[0]*X[2]+dFdu[4]*X[0]*X[1])
            if 'S' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBRS12:'+rbsx] += rpd*(dFdu[5]*X[1]-2.*dFdu[1]*X[2])
                dFdvDict[pfx+'RBRS13:'+rbsx] += rpd*(-dFdu[5]*X[2]+2.*dFdu[2]*X[1])
                dFdvDict[pfx+'RBRS21:'+rbsx] += rpd*(-dFdu[4]*X[0]+2.*dFdu[0]*X[2])
                dFdvDict[pfx+'RBRS23:'+rbsx] += rpd*(dFdu[4]*X[2]-2.*dFdu[2]*X[0])
                dFdvDict[pfx+'RBRS31:'+rbsx] += rpd*(dFdu[3]*X[0]-2.*dFdu[0]*X[1])
                dFdvDict[pfx+'RBRS32:'+rbsx] += rpd*(-dFdu[3]*X[1]+2.*dFdu[1]*X[0])
                dFdvDict[pfx+'RBRSAA:'+rbsx] += rpd*(dFdu[4]*X[1]-dFdu[3]*X[2])
                dFdvDict[pfx+'RBRSBB:'+rbsx] += rpd*(dFdu[5]*X[0]-dFdu[3]*X[2])
            if 'U' in RBObj['ThermalMotion'][0]:
                dFdvDict[pfx+'RBRU:'+rbsx] += dFdvDict[pfx+'AUiso:'+str(AtLookup[atId])]
    
################################################################################
##### Penalty & restraint functions 
################################################################################

def penaltyFxn(HistoPhases,calcControls,parmDict,varyList):
    'Compute user-supplied and built-in restraint functions'
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    pNames = []
    pVals = []
    pWt = []
    negWt = {}
    pWsum = {}
    pWnum = {}
    for phase in Phases:
        pId = Phases[phase]['pId']
        negWt[pId] = Phases[phase]['General']['Pawley neg wt']
        General = Phases[phase]['General']
        cx,ct,cs,cia = General['AtomPtrs']
        textureData = General['SH Texture']
        SGData = General['SGData']
        Atoms = Phases[phase]['Atoms']
        AtLookup = G2mth.FillAtomLookUp(Phases[phase]['Atoms'],cia+8)
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        if phase not in restraintDict:
            continue
        phaseRest = restraintDict[phase]
        names = [['Bond','Bonds'],['Angle','Angles'],['Plane','Planes'],
            ['Chiral','Volumes'],['Torsion','Torsions'],['Rama','Ramas'],
            ['ChemComp','Sites'],['Texture','HKLs'],['General','General'],]
        for name,rest in names:
            pWsum[name] = 0.
            pWnum[name] = 0
            if name not in phaseRest:
                continue
            itemRest = phaseRest[name]
            if itemRest[rest] and itemRest['Use']:
                wt = itemRest.get('wtFactor',1.)
                if name in ['Bond','Angle','Plane','Chiral']:
                    for i,[indx,ops,obs,esd] in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        if name == 'Bond':
                            calc = G2mth.getRestDist(XYZ,Amat)
                        elif name == 'Angle':
                            calc = G2mth.getRestAngle(XYZ,Amat)
                        elif name == 'Plane':
                            calc = G2mth.getRestPlane(XYZ,Amat)
                        elif name == 'Chiral':
                            calc = G2mth.getRestChiral(XYZ,Amat)
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)
                        pWsum[name] += wt*((obs-calc)/esd)**2
                        pWnum[name] += 1
                elif name in ['Torsion','Rama']:
                    coeffDict = itemRest['Coeff']
                    for i,[indx,ops,cofName,esd] in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        if name == 'Torsion':
                            tor = G2mth.getRestTorsion(XYZ,Amat)
                            restr,calc = G2mth.calcTorsionEnergy(tor,coeffDict[cofName])
                        else:
                            phi,psi = G2mth.getRestRama(XYZ,Amat)
                            restr,calc = G2mth.calcRamaEnergy(phi,psi,coeffDict[cofName])                               
                        pVals.append(restr)
                        pWt.append(wt/esd**2)
                        pWsum[name] += wt*(restr/esd)**2
                        pWnum[name] += 1
                elif name == 'ChemComp':
                    for i,[indx,factors,obs,esd] in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs+1))
                        frac = np.array(G2mth.GetAtomFracByID(pId,parmDict,AtLookup,indx))
                        calc = np.sum(mul*frac*factors)
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)                    
                        pWsum[name] += wt*((obs-calc)/esd)**2
                        pWnum[name] += 1
                elif name == 'Texture':
                    SHkeys = list(textureData['SH Coeff'][1].keys())
                    SHCoef = G2mth.GetSHCoeff(pId,parmDict,SHkeys)
                    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
                    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
                    for i,[hkl,grid,esd1,ifesd2,esd2] in enumerate(itemRest[rest]):
                        PH = np.array(hkl)
                        phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                        ODFln = G2lat.Flnh(False,SHCoef,phi,beta,SGData)
                        R,P,Z = G2mth.getRestPolefig(ODFln,SamSym[textureData['Model']],grid)
                        Z1 = ma.masked_greater(Z,0.0)           #is this + or -?
                        IndZ1 = np.array(ma.nonzero(Z1))
                        for ind in IndZ1.T:
                            pNames.append('%d:%s:%d:%.2f:%.2f'%(pId,name,i,R[ind[0],ind[1]],P[ind[0],ind[1]]))
                            pVals.append(Z1[ind[0]][ind[1]])
                            pWt.append(wt/esd1**2)
                            pWsum[name] += wt*(-Z1[ind[0]][ind[1]]/esd1)**2
                            pWnum[name] += 1
                        if ifesd2:
                            Z2 = 1.-Z
                            for ind in np.ndindex(grid,grid):
                                pNames.append('%d:%s:%d:%.2f:%.2f'%(pId,name+'-unit',i,R[ind[0],ind[1]],P[ind[0],ind[1]]))
                                pVals.append(Z2[ind[0]][ind[1]])
                                pWt.append(wt/esd2**2)
                                pWsum[name] += wt*(Z2/esd2)**2
                                pWnum[name] += 1
                elif name == 'General':
                    for i,(eq,obs,esd) in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        calcobj = G2obj.ExpressionCalcObj(eq)
                        calcobj.SetupCalc(parmDict)
                        calc = calcobj.EvalExpression()
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)                    
                        pWsum[name] += wt*((obs-calc)/esd)**2
                        pWnum[name] += 1
        
    for phase in Phases:
        name = 'SH-Pref.Ori.'
        pId = Phases[phase]['pId']
        General = Phases[phase]['General']
        SGData = General['SGData']
        cell = General['Cell'][1:7]
        pWsum[name] = 0.0
        pWnum[name] = 0
        for hist in Phases[phase]['Histograms']:
            if not Phases[phase]['Histograms'][hist]['Use']:
                continue
            if hist in Histograms and 'PWDR' in hist:
                hId = Histograms[hist]['hId']
                phfx = '%d:%d:'%(pId,hId)
                if calcControls[phfx+'poType'] == 'SH':
                    toler = calcControls[phfx+'SHtoler']
                    wt = 1./toler**2
                    HKLs = np.array(calcControls[phfx+'SHhkl'])
                    SHnames = calcControls[phfx+'SHnames']
                    SHcof = dict(zip(SHnames,[parmDict[phfx+cof] for cof in SHnames]))
                    for i,PH in enumerate(HKLs):
                        phi,beta = G2lat.CrsAng(PH,cell,SGData)
                        SH3Coef = {}
                        for item in SHcof:
                            L,N = eval(item.strip('C'))
                            SH3Coef['C%d,0,%d'%(L,N)] = SHcof[item]                        
                        ODFln = G2lat.Flnh(False,SH3Coef,phi,beta,SGData)
                        X = np.linspace(0,90.0,26)
                        Y = ma.masked_greater(G2lat.polfcal(ODFln,'0',X,0.0),0.0)       #+ or -?
                        IndY = ma.nonzero(Y)
                        for ind in IndY[0]:
                            pNames.append('%d:%d:%s:%d:%.2f'%(pId,hId,name,i,X[ind]))
                            pVals.append(Y[ind])
                            pWt.append(wt)
                            pWsum[name] += wt*(Y[ind])**2
                            pWnum[name] += 1
    pWsum['PWLref'] = 0.
    pWnum['PWLref'] = 0
    for item in varyList:
        if 'PWLref' in item and parmDict[item] < 0.:
            pId = int(item.split(':')[0])
            if negWt[pId]:
                pNames.append(item)
                pVals.append(parmDict[item])
                pWt.append(negWt[pId])
                pWsum['PWLref'] += negWt[pId]*(parmDict[item])**2
                pWnum['PWLref'] += 1
    pVals = np.array(pVals)
    pWt = np.array(pWt)         #should this be np.sqrt?
    return pNames,pVals,pWt,pWsum,pWnum
    
def penaltyDeriv(pNames,pVal,HistoPhases,calcControls,parmDict,varyList):
    '''Compute derivatives on user-supplied and built-in restraint 
    (penalty) functions

    where pNames is list of restraint labels

    returns pDerv with partial derivatives by variable# in varList and 
       restraint# in pNames (pDerv[variable#][restraint#])
    '''
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    pDerv = np.zeros((len(varyList),len(pVal)))
    for pName in pNames: # loop over restraints
        if 'General' == pName.split(':')[1]:
            # initialize for General restraint(s) here
            GeneralInit = True
            parmDict0 = parmDict.copy()
            # setup steps for each parameter
            stepDict = {}
            for parm in varyList:
                stepDict[parm] = G2obj.getVarStep(parm,parmDict)
            break
    for phase in Phases:
#        if phase not in restraintDict:
#            continue
        pId = Phases[phase]['pId']
        General = Phases[phase]['General']
        cx,ct,cs,cia = General['AtomPtrs']
        SGData = General['SGData']
        Atoms = Phases[phase]['Atoms']
        AtLookup = G2mth.FillAtomLookUp(Phases[phase]['Atoms'],cia+8)
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        textureData = General['SH Texture']

        SHkeys = list(textureData['SH Coeff'][1].keys())
        SHCoef = G2mth.GetSHCoeff(pId,parmDict,SHkeys)
        shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
        SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
        sam = SamSym[textureData['Model']]
        phaseRest = restraintDict.get(phase,{})
        names = {'Bond':'Bonds','Angle':'Angles','Plane':'Planes',
            'Chiral':'Volumes','Torsion':'Torsions','Rama':'Ramas',
            'ChemComp':'Sites','Texture':'HKLs'}
        lasthkl = np.array([0,0,0])
        for ip,pName in enumerate(pNames): # loop over restraints
            pnames = pName.split(':')
            if pId == int(pnames[0]):
                name = pnames[1]
                if 'PWL' in pName:
                    pDerv[varyList.index(pName)][ip] += 1.
                    continue
                elif 'SH-' in pName:
                    continue
                Id = int(pnames[2]) 
                itemRest = phaseRest[name]
                if name in ['Bond','Angle','Plane','Chiral']:
                    indx,ops,obs,esd = itemRest[names[name]][Id]
                    dNames = []
                    for ind in indx:
                        dNames += [str(pId)+'::dA'+Xname+':'+str(AtLookup[ind]) for Xname in ['x','y','z']]
                    XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                    if name == 'Bond':
                        deriv = G2mth.getRestDeriv(G2mth.getRestDist,XYZ,Amat,ops,SGData)
                    elif name == 'Angle':
                        deriv = G2mth.getRestDeriv(G2mth.getRestAngle,XYZ,Amat,ops,SGData)
                    elif name == 'Plane':
                        deriv = G2mth.getRestDeriv(G2mth.getRestPlane,XYZ,Amat,ops,SGData)
                    elif name == 'Chiral':
                        deriv = G2mth.getRestDeriv(G2mth.getRestChiral,XYZ,Amat,ops,SGData)
                elif name in ['Torsion','Rama']:
                    coffDict = itemRest['Coeff']
                    indx,ops,cofName,esd = itemRest[names[name]][Id]
                    dNames = []
                    for ind in indx:
                        dNames += [str(pId)+'::dA'+Xname+':'+str(AtLookup[ind]) for Xname in ['x','y','z']]
                    XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                    if name == 'Torsion':
                        deriv = G2mth.getTorsionDeriv(XYZ,Amat,coffDict[cofName])
                    else:
                        deriv = G2mth.getRamaDeriv(XYZ,Amat,coffDict[cofName])
                elif name == 'ChemComp':
                    indx,factors,obs,esd = itemRest[names[name]][Id]
                    dNames = []
                    for ind in indx:
                        dNames += [str(pId)+'::Afrac:'+str(AtLookup[ind])]
                        mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs+1))
                        deriv = mul*factors
                elif 'Texture' in name:
                    deriv = []
                    dNames = []
                    hkl,grid,esd1,ifesd2,esd2 = itemRest[names[name]][Id]
                    hkl = np.array(hkl)
                    if np.any(lasthkl-hkl):
                        phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                        ODFln = G2lat.Flnh(False,SHCoef,phi,beta,SGData)
                        lasthkl = copy.copy(hkl)                        
                    if 'unit' in name:
                        pass
                    else:
                        gam = float(pnames[3])
                        psi = float(pnames[4])
                        for SHname in ODFln:
                            l,m,n = eval(SHname[1:])
                            Ksl = G2lat.GetKsl(l,m,sam,psi,gam)[0]
                            dNames += [str(pId)+'::'+SHname]
                            deriv.append(-ODFln[SHname][0]*Ksl/SHCoef[SHname])
                elif name == 'General':
                    deriv = []
                    dNames = []
                    eq,obs,esd = itemRest[name][Id]
                    calcobj = G2obj.ExpressionCalcObj(eq)
                    parmlist = list(eq.assgnVars.values()) # parameters used in this expression
                    for parm in parmlist: # expand list if any parms are determined by constraints
                        if parm in G2mv.dependentVars:
                            parmlist += G2mv.independentVars
                            break
                    for ind,var in enumerate(varyList):
                        drv = 0
                        if var in parmlist:
                            step = stepDict.get(var,1e-5)
                            calc = []
                            # apply step to parameter
                            oneparm = True
                            for s in -step,2*step:
                                parmDict[var] += s
                                # extend shift if needed to other parameters
                                if var in G2mv.independentVars:
                                    G2mv.Dict2Map(parmDict,[])
                                    oneparm = False
                                elif var in G2mv.dependentVars:
                                    G2mv.Map2Dict(parmDict,[])
                                    oneparm = False
                                if 'RB' in var:
                                    ApplyRBModels(parmDict,Phases,rigidbodyDict)
# test
                                    oneparm = False
                                calcobj.SetupCalc(parmDict)
                                calc.append(calcobj.EvalExpression())
                            drv = (calc[1]-calc[0])*.5/step
                            # restore the dict
                            if oneparm:
                                parmDict[var] = parmDict0[var]
                            else:
                                parmDict = parmDict0.copy()
                        else:
                            drv = 0
                        pDerv[ind][ip] = drv
                # Add derivatives into matrix, if needed
                for dName,drv in zip(dNames,deriv):
                    try:
                        ind = varyList.index(dName)
                        pDerv[ind][ip] += drv
                    except ValueError:
                        pass
        
        lasthkl = np.array([0,0,0])
        for ip,pName in enumerate(pNames):
            deriv = []
            dNames = []
            pnames = pName.split(':')
            if 'SH-' in pName and pId == int(pnames[0]):
                hId = int(pnames[1])
                phfx = '%d:%d:'%(pId,hId)
                psi = float(pnames[4])
                HKLs = calcControls[phfx+'SHhkl']
                SHnames = calcControls[phfx+'SHnames']
                SHcof = dict(zip(SHnames,[parmDict[phfx+cof] for cof in SHnames]))
                hkl = np.array(HKLs[int(pnames[3])])     
                if np.any(lasthkl-hkl):
                    phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                    SH3Coef = {}
                    for item in SHcof:
                        L,N = eval(item.strip('C'))
                        SH3Coef['C%d,0,%d'%(L,N)] = SHcof[item]                        
                    ODFln = G2lat.Flnh(False,SH3Coef,phi,beta,SGData)
                    lasthkl = copy.copy(hkl)                        
                for SHname in SHnames:
                    l,n = eval(SHname[1:])
                    SH3name = 'C%d,0,%d'%(l,n)
                    Ksl = G2lat.GetKsl(l,0,'0',psi,0.0)[0]
                    dNames += [phfx+SHname]
                    deriv.append(ODFln[SH3name][0]*Ksl/SHcof[SHname])
            for dName,drv in zip(dNames,deriv):
                try:
                    ind = varyList.index(dName)
                    pDerv[ind][ip] += drv
                except ValueError:
                    pass
    return pDerv

################################################################################
##### Function & derivative calculations
################################################################################        
                    
def GetAtomFXU(pfx,calcControls,parmDict):
    'Needs a doc string'
    Natoms = calcControls['Natoms'][pfx]
    Tdata = Natoms*[' ',]
    Mdata = np.zeros(Natoms)
    IAdata = Natoms*[' ',]
    Fdata = np.zeros(Natoms)
    Xdata = np.zeros((3,Natoms))
    dXdata = np.zeros((3,Natoms))
    Uisodata = np.zeros(Natoms)
    Uijdata = np.zeros((6,Natoms))
    Gdata = np.zeros((3,Natoms))
    keys = {'Atype:':Tdata,'Amul:':Mdata,'Afrac:':Fdata,'AI/A:':IAdata,
        'dAx:':dXdata[0],'dAy:':dXdata[1],'dAz:':dXdata[2],
        'Ax:':Xdata[0],'Ay:':Xdata[1],'Az:':Xdata[2],'AUiso:':Uisodata,
        'AU11:':Uijdata[0],'AU22:':Uijdata[1],'AU33:':Uijdata[2],
        'AU12:':Uijdata[3],'AU13:':Uijdata[4],'AU23:':Uijdata[5],
        'AMx:':Gdata[0],'AMy:':Gdata[1],'AMz:':Gdata[2],}
    for iatm in range(Natoms):
        for key in keys:
            parm = pfx+key+str(iatm)
            if parm in parmDict:
                keys[key][iatm] = parmDict[parm]
    Fdata = np.where(Fdata,Fdata,1.e-8)         #avoid divide by zero in derivative calc.
    Gdata = np.where(Gdata,Gdata,1.e-8)         #avoid divide by zero in derivative calc.
    
    return Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata
    
def GetAtomSSFXU(pfx,calcControls,parmDict):
    'Needs a doc string'
    Natoms = calcControls['Natoms'][pfx]
    maxSSwave = calcControls['maxSSwave'][pfx]
    Nwave = {'F':maxSSwave['Sfrac'],'X':maxSSwave['Spos'],'Y':maxSSwave['Spos'],'Z':maxSSwave['Spos'],
        'U':maxSSwave['Sadp'],'M':maxSSwave['Smag'],'T':maxSSwave['Spos']}
    XSSdata = np.zeros((6,maxSSwave['Spos'],Natoms))
    FSSdata = np.zeros((2,maxSSwave['Sfrac'],Natoms))
    USSdata = np.zeros((12,maxSSwave['Sadp'],Natoms))
    MSSdata = np.zeros((6,maxSSwave['Smag'],Natoms))
    waveTypes = []
    keys = {'Fsin:':FSSdata[0],'Fcos:':FSSdata[1],'Fzero:':FSSdata[0],'Fwid:':FSSdata[1],
        'Tmin:':XSSdata[0],'Tmax:':XSSdata[1],'Xmax:':XSSdata[2],'Ymax:':XSSdata[3],'Zmax:':XSSdata[4],
        'Xsin:':XSSdata[0],'Ysin:':XSSdata[1],'Zsin:':XSSdata[2],'Xcos:':XSSdata[3],'Ycos:':XSSdata[4],'Zcos:':XSSdata[5],
        'U11sin:':USSdata[0],'U22sin:':USSdata[1],'U33sin:':USSdata[2],'U12sin:':USSdata[3],'U13sin:':USSdata[4],'U23sin:':USSdata[5],
        'U11cos:':USSdata[6],'U22cos:':USSdata[7],'U33cos:':USSdata[8],'U12cos:':USSdata[9],'U13cos:':USSdata[10],'U23cos:':USSdata[11],
        'MXsin:':MSSdata[0],'MYsin:':MSSdata[1],'MZsin:':MSSdata[2],'MXcos:':MSSdata[3],'MYcos:':MSSdata[4],'MZcos:':MSSdata[5]}
    for iatm in range(Natoms):
        wavetype = [parmDict.get(pfx+kind+'waveType:'+str(iatm),'') for kind in ['F','P','A','M']]
        waveTypes.append(wavetype)
        for key in keys:
            for m in range(Nwave[key[0]]):
                parm = pfx+key+str(iatm)+':%d'%(m)
                if parm in parmDict:
                    keys[key][m][iatm] = parmDict[parm]
    return waveTypes,FSSdata,XSSdata,USSdata,MSSdata
    
def StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Compute structure factors for all h,k,l for phase
    puts the result, F^2, in each ref[8] in refList
    operates on blocks of 100 reflections for speed
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:

    '''        
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Amat,Bmat = G2lat.Gmat2AB(G)
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
    TwinLaw = np.array([[[1,0,0],[0,1,0],[0,0,1]],])
    TwDict = refDict.get('TwDict',{})           
    if 'S' in calcControls[hfx+'histType']:
        NTL = calcControls[phfx+'NTL']
        NM = calcControls[phfx+'TwinNMN']+1
        TwinLaw = calcControls[phfx+'TwinLaw']
        TwinFr = np.array([parmDict[phfx+'TwinFr:'+str(i)] for i in range(len(TwinLaw))])
        TwinInv = list(np.where(calcControls[phfx+'TwinInv'],-1,1))
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    blkSize = 100       #no. of reflections in a block - size seems optimal
    nRef = refDict['RefList'].shape[0]
    SQ = 1./(2.*refDict['RefList'].T[4])**2
    if 'N' in calcControls[hfx+'histType']:
        dat = G2el.getBLvalues(BLtables)
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.ones((nRef,len(dat)))*list(dat.values())
    else:       #'X'
        dat = G2el.getFFvalues(FFtables,0.)
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.zeros((nRef,len(dat)))
        for iel,El in enumerate(refDict['FF']['El']):
            refDict['FF']['FF'].T[iel] = G2el.ScatFac(FFtables[El],SQ)
#reflection processing begins here - big arrays!
    iBeg = 0
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:3]                          #array(blkSize,3)
        H = np.squeeze(np.inner(H.T,TwinLaw))   #maybe array(blkSize,nTwins,3) or (blkSize,3)
        TwMask = np.any(H,axis=-1)
        if TwinLaw.shape[0] > 1 and TwDict: #need np.inner(TwinLaw[?],TwDict[iref][i])*TwinInv[i]
            for ir in range(blkSize):
                iref = ir+iBeg
                if iref in TwDict:
                    for i in TwDict[iref]:
                        for n in range(NTL):
                            H[ir][i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
            TwMask = np.any(H,axis=-1)
        SQ = 1./(2.*refl.T[4])**2               #array(blkSize)
        SQfactor = 4.0*SQ*twopisq               #ditto prev.
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,len(SGT)*len(TwinLaw),axis=0)
            FPP = np.repeat(FPP.T,len(SGT)*len(TwinLaw),axis=0)
        Uniq = np.inner(H,SGMT)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T).T+Phi.T).T
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT)*len(TwinLaw),axis=1).T
        HbH = -np.sum(Uniq.T*np.swapaxes(np.inner(bij,Uniq),2,-1),axis=1)
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0).T
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/len(SGMT)
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,len(SGT)*len(TwinLaw),axis=0)
        Bab = np.repeat(parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor),len(SGT)*len(TwinLaw))
        if 'T' in calcControls[hfx+'histType']: #fa,fb are 2 X blkSize X nTwin X nOps x nAtoms
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-np.reshape(Flack*FPP,sinp.shape)*sinp*Tcorr])
            fb = np.array([np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr,np.reshape(Flack*FPP,cosp.shape)*cosp*Tcorr])
        else:
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-Flack*FPP*sinp*Tcorr])
            fb = np.array([np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr,Flack*FPP*cosp*Tcorr])
        fas = np.sum(np.sum(fa,axis=-1),axis=-1)  #real 2 x blkSize x nTwin; sum over atoms & uniq hkl
        fbs = np.sum(np.sum(fb,axis=-1),axis=-1)  #imag 
        if SGData['SGInv']: #centrosymmetric; B=0
            fbs[0] *= 0.
            fas[1] *= 0.
        if 'P' in calcControls[hfx+'histType']:     #PXC, PNC & PNT: F^2 = A[0]^2 + A[1]^2 + B[0]^2 + B[1]^2
            refl.T[9] = np.sum(fas**2,axis=0)+np.sum(fbs**2,axis=0) #add fam**2 & fbm**2 here    
            refl.T[10] = atan2d(fbs[0],fas[0])  #ignore f' & f"
        else:                                       #HKLF: F^2 = (A[0]+A[1])^2 + (B[0]+B[1])^2
            if len(TwinLaw) > 1:
                refl.T[9] = np.sum(fas[:,:,0],axis=0)**2+np.sum(fbs[:,:,0],axis=0)**2   #FcT from primary twin element
                refl.T[7] = np.sum(TwinFr*TwMask*np.sum(fas,axis=0)**2,axis=-1)+   \
                    np.sum(TwinFr*TwMask*np.sum(fbs,axis=0)**2,axis=-1)                        #Fc sum over twins
                refl.T[10] = atan2d(fbs[0].T[0],fas[0].T[0])  #ignore f' & f" & use primary twin
            else:   # checked correct!!
                refl.T[9] = np.sum(fas,axis=0)**2+np.sum(fbs,axis=0)**2
                refl.T[7] = np.copy(refl.T[9])                
                refl.T[10] = atan2d(fbs[0],fas[0])  #ignore f' & f"
#                refl.T[10] = atan2d(np.sum(fbs,axis=0),np.sum(fas,axis=0)) #include f' & f"
        iBeg += blkSize
#    print 'sf time %.4f, nref %d, blkSize %d'%(time.time()-time0,nRef,blkSize)
    
def StructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    '''Compute structure factor derivatives on blocks of reflections - for powders/nontwins only
    faster than StructureFactorDerv - correct for powders/nontwins!!
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filled in below
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict parmDict:
    
    :returns: dict dFdvDict: dictionary of derivatives
    '''
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Amat,Bmat = G2lat.Gmat2AB(G)
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return {}
    mSize = len(Mdata)
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    dFdvDict = {}
    dFdfr = np.zeros((nRef,mSize))
    dFdx = np.zeros((nRef,mSize,3))
    dFdui = np.zeros((nRef,mSize))
    dFdua = np.zeros((nRef,mSize,6))
    dFdbab = np.zeros((nRef,2))
    dFdfl = np.zeros((nRef))
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
#    time0 = time.time()
#reflection processing begins here - big arrays!
    iBeg = 0
    blkSize = 32       #no. of reflections in a block - optimized for speed
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:3].T
        SQ = 1./(2.*refl.T[4])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,len(SGT),axis=0)
            FPP = np.repeat(FPP.T,len(SGT),axis=0)
        dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
        Bab = np.repeat(parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor),len(SGT))
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,len(SGT),axis=0)
        Uniq = np.inner(H,SGMT)             # array(nSGOp,3)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T).T+Phi.T).T
        sinp = np.sin(phase)        #refBlk x nOps x nAtoms
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(SGT)
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT),axis=1).T
        HbH = -np.sum(Uniq.T*np.swapaxes(np.inner(bij,Uniq),2,-1),axis=1)
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0).T
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/len(SGMT)
        Hij = np.array([Mast*np.multiply.outer(U,U) for U in np.reshape(Uniq,(-1,3))])      #Nref*Nops,3,3
        Hij = np.reshape(np.array([G2lat.UijtoU6(uij) for uij in Hij]),(-1,len(SGT),6))     #Nref,Nops,6
        fot = np.reshape(((FF+FP).T-Bab).T,cosp.shape)*Tcorr
        if len(FPP.shape) > 1:
            fotp = np.reshape(FPP,cosp.shape)*Tcorr
        else:
            fotp = FPP*Tcorr     
        if 'T' in calcControls[hfx+'histType']:
            fa = np.array([fot*cosp,-np.reshape(Flack*FPP,sinp.shape)*sinp*Tcorr])
            fb = np.array([fot*sinp,np.reshape(Flack*FPP,cosp.shape)*cosp*Tcorr])
        else:
            fa = np.array([fot*cosp,-Flack*FPP*sinp*Tcorr])
            fb = np.array([fot*sinp,Flack*FPP*cosp*Tcorr])
        fas = np.sum(np.sum(fa,axis=-1),axis=-1)      #real sum over atoms & unique hkl array(2,refBlk,nTwins)
        fbs = np.sum(np.sum(fb,axis=-1),axis=-1)      #imag sum over atoms & uniq hkl
        fax = np.array([-fot*sinp,-fotp*cosp])   #positions array(2,refBlk,nEqv,nAtoms)
        fbx = np.array([fot*cosp,-fotp*sinp])
        #sum below is over Uniq
        dfadfr = np.sum(fa/occ,axis=-2)        #array(2,refBlk,nAtom) Fdata != 0 avoids /0. problem 
        dfadba = np.sum(-cosp*Tcorr,axis=-2)  #array(refBlk,nAtom)
        dfadx = np.sum(twopi*Uniq[nxs,:,nxs,:,:]*np.swapaxes(fax,-2,-1)[:,:,:,:,nxs],axis=-2)
        dfadui = np.sum(-SQfactor[nxs,:,nxs,nxs]*fa,axis=-2) #array(Ops,refBlk,nAtoms)
        dfadua = np.sum(-Hij[nxs,:,nxs,:,:]*np.swapaxes(fa,-2,-1)[:,:,:,:,nxs],axis=-2)
        # array(2,refBlk,nAtom,3) & array(2,refBlk,nAtom,6)
        if not SGData['SGInv']:
            dfbdfr = np.sum(fb/occ,axis=-2)        #Fdata != 0 avoids /0. problem
            dfbdba = np.sum(-sinp*Tcorr,axis=-2)
            dfadfl = np.sum(np.sum(-fotp*sinp,axis=-1),axis=-1)
            dfbdfl = np.sum(np.sum(fotp*cosp,axis=-1),axis=-1)
            dfbdx = np.sum(twopi*Uniq[nxs,:,nxs,:,:]*np.swapaxes(fbx,-2,-1)[:,:,:,:,nxs],axis=-2)           
            dfbdui = np.sum(-SQfactor[nxs,:,nxs,nxs]*fb,axis=-2)
            dfbdua = np.sum(-Hij[nxs,:,nxs,:,:]*np.swapaxes(fb,-2,-1)[:,:,:,:,nxs],axis=-2)
        else:
            dfbdfr = np.zeros_like(dfadfr)
            dfbdx = np.zeros_like(dfadx)
            dfbdui = np.zeros_like(dfadui)
            dfbdua = np.zeros_like(dfadua)
            dfbdba = np.zeros_like(dfadba)
            dfadfl = 0.0
            dfbdfl = 0.0
        #NB: the above have been checked against PA(1:10,1:2) in strfctr.for for Al2O3!    
        SA = fas[0]+fas[1]
        SB = fbs[0]+fbs[1]
        if 'P' in calcControls[hfx+'histType']: #checked perfect for centro & noncentro
            dFdfr[iBeg:iFin] = 2.*np.sum(fas[:,:,nxs]*dfadfr+fbs[:,:,nxs]*dfbdfr,axis=0)*Mdata/len(SGMT)
            dFdx[iBeg:iFin] = 2.*np.sum(fas[:,:,nxs,nxs]*dfadx+fbs[:,:,nxs,nxs]*dfbdx,axis=0)
            dFdui[iBeg:iFin] = 2.*np.sum(fas[:,:,nxs]*dfadui+fbs[:,:,nxs]*dfbdui,axis=0)
            dFdua[iBeg:iFin] = 2.*np.sum(fas[:,:,nxs,nxs]*dfadua+fbs[:,:,nxs,nxs]*dfbdua,axis=0)
        else:
            dFdfr[iBeg:iFin] = (2.*SA[:,nxs]*(dfadfr[0]+dfadfr[1])+2.*SB[:,nxs]*(dfbdfr[0]+dfbdfr[1]))*Mdata/len(SGMT)
            dFdx[iBeg:iFin] = 2.*SA[:,nxs,nxs]*(dfadx[0]+dfadx[1])+2.*SB[:,nxs,nxs]*(dfbdx[0]+dfbdx[1])
            dFdui[iBeg:iFin] = 2.*SA[:,nxs]*(dfadui[0]+dfadui[1])+2.*SB[:,nxs]*(dfbdui[0]+dfbdui[1])
            dFdua[iBeg:iFin] = 2.*SA[:,nxs,nxs]*(dfadua[0]+dfadua[1])+2.*SB[:,nxs,nxs]*(dfbdua[0]+dfbdua[1])
            dFdfl[iBeg:iFin] = -SA*dfadfl-SB*dfbdfl  #array(nRef,)
        dFdbab[iBeg:iFin] = 2.*(fas[0,nxs]*np.array([np.sum(dfadba.T*dBabdA,axis=0),np.sum(-dfadba.T*parmDict[phfx+'BabA']*SQfactor*dBabdA,axis=0)])+ \
                            fbs[0,nxs]*np.array([np.sum(dfbdba.T*dBabdA,axis=0),np.sum(-dfbdba.T*parmDict[phfx+'BabA']*SQfactor*dBabdA,axis=0)])).T
        iBeg += blkSize
#    print 'derv time %.4f, nref %d, blkSize %d'%(time.time()-time0,nRef,blkSize)
        #loop over atoms - each dict entry is list of derivatives for all the reflections
    for i in range(len(Mdata)):
        dFdvDict[pfx+'Afrac:'+str(i)] = dFdfr.T[i]
        dFdvDict[pfx+'dAx:'+str(i)] = dFdx.T[0][i]
        dFdvDict[pfx+'dAy:'+str(i)] = dFdx.T[1][i]
        dFdvDict[pfx+'dAz:'+str(i)] = dFdx.T[2][i]
        dFdvDict[pfx+'AUiso:'+str(i)] = dFdui.T[i]
        dFdvDict[pfx+'AU11:'+str(i)] = dFdua.T[0][i]
        dFdvDict[pfx+'AU22:'+str(i)] = dFdua.T[1][i]
        dFdvDict[pfx+'AU33:'+str(i)] = dFdua.T[2][i]
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*dFdua.T[3][i]
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*dFdua.T[4][i]
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*dFdua.T[5][i]
    dFdvDict[phfx+'Flack'] = 4.*dFdfl.T
    dFdvDict[phfx+'BabA'] = dFdbab.T[0]
    dFdvDict[phfx+'BabU'] = dFdbab.T[1]
    return dFdvDict
    
def MagStructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Compute neutron magnetic structure factors for all h,k,l for phase
    puts the result, F^2, in each ref[8] in refList
    operates on blocks of 100 reflections for speed
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:
        
    :returns: copy of new refList - used in calculating numerical derivatives

    '''        
    g = nl.inv(G)
    ast = np.sqrt(np.diag(G))
    ainv = np.sqrt(np.diag(g))
    GS = G/np.outer(ast,ast)
    Ginv = g/np.outer(ainv,ainv)
    uAmat = G2lat.Gmat2AB(GS)[0]
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    Ncen = len(SGData['SGCen'])
    Nops = len(SGMT)*Ncen
    if not SGData['SGFixed']:
        Nops *= (1+SGData['SGInv'])
    MFtables = calcControls['MFtables']
    Bmat = G2lat.Gmat2AB(G)[1]
    TwinLaw = np.ones(1)
#    TwinLaw = np.array([[[1,0,0],[0,1,0],[0,0,1]],])
#    TwDict = refDict.get('TwDict',{})           
#    if 'S' in calcControls[hfx+'histType']:
#        NTL = calcControls[phfx+'NTL']
#        NM = calcControls[phfx+'TwinNMN']+1
#        TwinLaw = calcControls[phfx+'TwinLaw']
#        TwinFr = np.array([parmDict[phfx+'TwinFr:'+str(i)] for i in range(len(TwinLaw))])
#        TwinInv = list(np.where(calcControls[phfx+'TwinInv'],-1,1))
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return
    Mag = np.array([np.sqrt(np.inner(mag,np.inner(mag,Ginv))) for mag in Gdata.T])
    Gdata = np.inner(Gdata.T,np.swapaxes(SGMT,1,2)).T            #apply sym. ops.
    if SGData['SGInv'] and not SGData['SGFixed']:
        Gdata = np.hstack((Gdata,-Gdata))       #inversion if any
    Gdata = np.hstack([Gdata for icen in range(Ncen)])        #dup over cell centering--> [Mxyz,nops,natms]
    Gdata = SGData['MagMom'][nxs,:,nxs]*Gdata   #flip vectors according to spin flip * det(opM)
    Mag = np.tile(Mag[:,nxs],Nops).T  #make Mag same length as Gdata
    Kdata = np.inner(Gdata.T,uAmat).T
    Kmean = np.mean(np.sqrt(np.sum(Kdata**2,axis=0)),axis=0)
    Kdata /= Kmean     #Cartesian unit vectors
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    blkSize = 100       #no. of reflections in a block - size seems optimal
    nRef = refDict['RefList'].shape[0]
    SQ = 1./(2.*refDict['RefList'].T[4])**2
    refDict['FF']['El'] = list(MFtables.keys())
    refDict['FF']['MF'] = np.zeros((nRef,len(MFtables)))
    for iel,El in enumerate(refDict['FF']['El']):
        refDict['FF']['MF'].T[iel] = G2el.MagScatFac(MFtables[El],SQ)
#reflection processing begins here - big arrays!
    iBeg = 0
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:3].T                          #array(blkSize,3)
#        H = np.squeeze(np.inner(H.T,TwinLaw))   #maybe array(blkSize,nTwins,3) or (blkSize,3)
#        TwMask = np.any(H,axis=-1)
#        if TwinLaw.shape[0] > 1 and TwDict: #need np.inner(TwinLaw[?],TwDict[iref][i])*TwinInv[i]
#            for ir in range(blkSize):
#                iref = ir+iBeg
#                if iref in TwDict:
#                    for i in TwDict[iref]:
#                        for n in range(NTL):
#                            H[ir][i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
#            TwMask = np.any(H,axis=-1)
        SQ = 1./(2.*refl.T[4])**2               #array(blkSize)
        SQfactor = 4.0*SQ*twopisq               #ditto prev.
        Uniq = np.inner(H,SGMT)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T).T+Phi.T).T
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT)*len(TwinLaw),axis=1).T
        HbH = -np.sum(Uniq.T*np.swapaxes(np.inner(bij,Uniq),2,-1),axis=1)
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0).T
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        MF = refDict['FF']['MF'][iBeg:iFin].T[Tindx].T   #Nref,Natm
        TMcorr = 0.539*(np.reshape(Tiso,Tuij.shape)*Tuij)[:,0,:]*Fdata*Mdata*MF/(2*Nops)     #Nref,Natm
        if SGData['SGInv']:
            if not SGData['SGFixed']:
                mphase = np.hstack((phase,-phase))  #OK
            else:
                mphase = phase
        else:
            mphase = phase                    #
        mphase = np.array([mphase+twopi*np.inner(cen,H)[:,nxs,nxs] for cen in SGData['SGCen']])
        mphase = np.concatenate(mphase,axis=1)              #Nref,full Nop,Natm
        sinm = np.sin(mphase)                               #ditto - match magstrfc.for
        cosm = np.cos(mphase)                               #ditto
        HM = np.inner(Bmat,H)                             #put into cartesian space
        HM = HM/np.sqrt(np.sum(HM**2,axis=0))               #Kdata = MAGS & HM = UVEC in magstrfc.for both OK
        eDotK = np.sum(HM[:,:,nxs,nxs]*Kdata[:,nxs,:,:],axis=0)
        Q = HM[:,:,nxs,nxs]*eDotK[nxs,:,:,:]-Kdata[:,nxs,:,:] #xyz,Nref,Nop,Natm = BPM in magstrfc.for OK
        fam = Q*TMcorr[nxs,:,nxs,:]*cosm[nxs,:,:,:]*Mag[nxs,nxs,:,:]    #ditto
        fbm = Q*TMcorr[nxs,:,nxs,:]*sinm[nxs,:,:,:]*Mag[nxs,nxs,:,:]    #ditto
        fams = np.sum(np.sum(fam,axis=-1),axis=-1)     #Mxyz,Nref  Sum(sum(fam,atoms),ops)
        fbms = np.sum(np.sum(fbm,axis=-1),axis=-1)     #ditto
        refl.T[9] = np.sum(fams**2,axis=0)+np.sum(fbms**2,axis=0)   #Sum(fams**2,Mxyz) Re + Im
        refl.T[7] = np.copy(refl.T[9])                
        refl.T[10] = atan2d(fbms[0],fams[0]) #- what is phase for mag refl?
#        if 'P' in calcControls[hfx+'histType']:     #PXC, PNC & PNT: F^2 = A[0]^2 + A[1]^2 + B[0]^2 + B[1]^2
#            refl.T[9] = np.sum(fas**2,axis=0)+np.sum(fbs**2,axis=0) #add fam**2 & fbm**2 here    
#            refl.T[10] = atan2d(fbs[0],fas[0])  #ignore f' & f"
#        else:                                       #HKLF: F^2 = (A[0]+A[1])^2 + (B[0]+B[1])^2
#            if len(TwinLaw) > 1:
#                refl.T[9] = np.sum(fas[:,:,0],axis=0)**2+np.sum(fbs[:,:,0],axis=0)**2   #FcT from primary twin element
#                refl.T[7] = np.sum(TwinFr*TwMask*np.sum(fas,axis=0)**2,axis=-1)+   \
#                    np.sum(TwinFr*TwMask*np.sum(fbs,axis=0)**2,axis=-1)                        #Fc sum over twins
#                refl.T[10] = atan2d(fbs[0].T[0],fas[0].T[0])  #ignore f' & f" & use primary twin
#            else:   # checked correct!!
#                refl.T[9] = np.sum(fas,axis=0)**2+np.sum(fbs,axis=0)**2
#                refl.T[7] = np.copy(refl.T[9])                
#                refl.T[10] = atan2d(fbs[0],fas[0])  #ignore f' & f"
##                refl.T[10] = atan2d(np.sum(fbs,axis=0),np.sum(fas,axis=0)) #include f' & f"
        iBeg += blkSize
#    print 'sf time %.4f, nref %d, blkSize %d'%(time.time()-time0,nRef,blkSize)
    return copy.deepcopy(refDict['RefList'])

def MagStructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    '''Compute magnetic structure factor derivatives numerically - for powders/nontwins only
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filled in below
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict parmDict:
    
    :returns: dict dFdvDict: dictionary of magnetic derivatives
    '''
    
    trefDict = copy.deepcopy(refDict)
    dM = 1.e-6
    dFdvDict = {}
    for parm in parmDict:
        if 'AM' in parm:
            parmDict[parm] += dM
            prefList = MagStructureFactor2(trefDict,G,hfx,pfx,SGData,calcControls,parmDict)
            parmDict[parm] -= 2*dM
            mrefList = MagStructureFactor2(trefDict,G,hfx,pfx,SGData,calcControls,parmDict)
            parmDict[parm] += dM
            dFdvDict[parm] = (prefList[:,9]-mrefList[:,9])/(2.*dM)
    return dFdvDict
            
def MagStructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    '''Compute nonmagnetic structure factor derivatives on blocks of reflections in magnetic structures - for powders/nontwins only
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filled in below
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict parmDict:
    
    :returns: dict dFdvDict: dictionary of derivatives
    '''
    
    g = nl.inv(G)
    ast = np.sqrt(np.diag(G))
    ainv = np.sqrt(np.diag(g))
    GS = G/np.outer(ast,ast)
    Ginv = g/np.outer(ainv,ainv)
    uAmat = G2lat.Gmat2AB(GS)[0]
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    Ncen = len(SGData['SGCen'])
    Nops = len(SGMT)*Ncen
    if not SGData['SGFixed']:
        Nops *= (1+SGData['SGInv'])
    Bmat = G2lat.Gmat2AB(G)[1]
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return {}
    mSize = len(Mdata)
    Mag = np.array([np.sqrt(np.inner(mag,np.inner(mag,Ginv))) for mag in Gdata.T])
    Gones = np.ones_like(Gdata)
    Gdata = np.inner(Gdata.T,np.swapaxes(SGMT,1,2)).T            #apply sym. ops.
    Gones = np.inner(Gones.T,SGMT).T
    if SGData['SGInv'] and not SGData['SGFixed']:
        Gdata = np.hstack((Gdata,-Gdata))       #inversion if any
        Gones = np.hstack((Gones,-Gones))       #inversion if any
    Gdata = np.hstack([Gdata for icen in range(Ncen)])        #dup over cell centering
    Gones = np.hstack([Gones for icen in range(Ncen)])        #dup over cell centering
    Gdata = SGData['MagMom'][nxs,:,nxs]*Gdata   #flip vectors according to spin flip
    Gones = SGData['MagMom'][nxs,:,nxs]*Gones   #flip vectors according to spin flip
    Mag = np.tile(Mag[:,nxs],Nops).T  #make Mag same length as Gdata
    Kdata = np.inner(Gdata.T,uAmat).T     #Cartesian unit vectors
    Kmean = np.mean(np.sqrt(np.sum(Kdata**2,axis=0)),axis=0)
    Kdata /= Kmean
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    dFdvDict = {}
    dFdfr = np.zeros((nRef,mSize))
    dFdx = np.zeros((nRef,mSize,3))
    dFdui = np.zeros((nRef,mSize))
    dFdua = np.zeros((nRef,mSize,6))
    time0 = time.time()
#reflection processing begins here - big arrays!
    iBeg = 0
    blkSize = 5       #no. of reflections in a block - optimized for speed
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:3].T
        SQ = 1./(2.*refl.T[4])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        Uniq = np.inner(H,SGMT)             # array(nSGOp,3)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T).T+Phi.T).T
        occ = Mdata*Fdata/Nops
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT),axis=1).T
        HbH = -np.sum(Uniq.T*np.swapaxes(np.inner(bij,Uniq),2,-1),axis=1)
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0).T
        Hij = np.array([Mast*np.multiply.outer(U,U) for U in np.reshape(Uniq,(-1,3))])
        Hij = np.reshape(np.array([G2lat.UijtoU6(uij) for uij in Hij]),(-1,len(SGT),6))
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        MF = refDict['FF']['MF'][iBeg:iFin].T[Tindx].T   #Nref,Natm
        TMcorr = 0.539*(np.reshape(Tiso,Tuij.shape)*Tuij)[:,0,:]*Fdata*Mdata*MF/(2*Nops)     #Nref,Natm
        if SGData['SGInv']:
            if not SGData['SGFixed']:
                mphase = np.hstack((phase,-phase))  #OK
                Uniq = np.hstack((Uniq,-Uniq))      #Nref,Nops,hkl
                Hij = np.hstack((Hij,Hij))
            else:
                mphase = phase
        else:
            mphase = phase                    #
        Hij = np.concatenate(np.array([Hij for cen in SGData['SGCen']]),axis=1)
        Uniq = np.hstack([Uniq for cen in SGData['SGCen']])
        mphase = np.array([mphase+twopi*np.inner(cen,H)[:,nxs,nxs] for cen in SGData['SGCen']])
        mphase = np.concatenate(mphase,axis=1)              #Nref,Nop,Natm
        sinm = np.sin(mphase)                               #ditto - match magstrfc.for
        cosm = np.cos(mphase)                               #ditto
        HM = np.inner(Bmat.T,H)                             #put into cartesian space
        HM = HM/np.sqrt(np.sum(HM**2,axis=0))               #unit cartesian vector for H
        eDotK = np.sum(HM[:,:,nxs,nxs]*Kdata[:,nxs,:,:],axis=0)
        Q = HM[:,:,nxs,nxs]*eDotK[nxs,:,:,:]-Kdata[:,nxs,:,:] #Mxyz,Nref,Nop,Natm = BPM in magstrfc.for OK
        
        fam = Q*TMcorr[nxs,:,nxs,:]*cosm[nxs,:,:,:]*Mag[nxs,nxs,:,:]    #Mxyz,Nref,Nop,Natm
        fbm = Q*TMcorr[nxs,:,nxs,:]*sinm[nxs,:,:,:]*Mag[nxs,nxs,:,:]
        fams = np.sum(np.sum(fam,axis=-1),axis=-1)                      #Mxyz,Nref
        fbms = np.sum(np.sum(fbm,axis=-1),axis=-1)
        famx = -Q*TMcorr[nxs,:,nxs,:]*Mag[nxs,nxs,:,:]*sinm[nxs,:,:,:]   #Mxyz,Nref,Nops,Natom
        fbmx = Q*TMcorr[nxs,:,nxs,:]*Mag[nxs,nxs,:,:]*cosm[nxs,:,:,:]
        #sums below are over Nops - real part
        dfadfr = np.sum(fam/occ,axis=2)        #array(Mxyz,refBlk,nAtom) Fdata != 0 avoids /0. problem deriv OK
        dfadx = np.sum(twopi*Uniq[nxs,:,:,nxs,:]*famx[:,:,:,:,nxs],axis=2)          #deriv OK
        dfadui = np.sum(-SQfactor[:,nxs,nxs]*fam,axis=2) #array(Ops,refBlk,nAtoms)  deriv OK
        dfadua = np.sum(-Hij[nxs,:,:,nxs,:]*fam[:,:,:,:,nxs],axis=2)            #deriv OK
        # imaginary part; array(3,refBlk,nAtom,3) & array(3,refBlk,nAtom,6)
        dfbdfr = np.sum(fbm/occ,axis=2)        #array(mxyz,refBlk,nAtom) Fdata != 0 avoids /0. problem 
        dfbdx = np.sum(twopi*Uniq[nxs,:,:,nxs,:]*fbmx[:,:,:,:,nxs],axis=2)
        dfbdui = np.sum(-SQfactor[:,nxs,nxs]*fbm,axis=2) #array(Ops,refBlk,nAtoms)
        dfbdua = np.sum(-Hij[nxs,:,:,nxs,:]*fbm[:,:,:,:,nxs],axis=2)
        #accumulate derivatives    
        dFdfr[iBeg:iFin] = 2.*np.sum((fams[:,:,nxs]*dfadfr+fbms[:,:,nxs]*dfbdfr)*Mdata/Nops,axis=0) #ok
        dFdx[iBeg:iFin] = 2.*np.sum(fams[:,:,nxs,nxs]*dfadx+fbms[:,:,nxs,nxs]*dfbdx,axis=0)         #ok
        dFdui[iBeg:iFin] = 2.*np.sum(fams[:,:,nxs]*dfadui+fbms[:,:,nxs]*dfbdui,axis=0)              #ok
        dFdua[iBeg:iFin] = 2.*np.sum(fams[:,:,nxs,nxs]*dfadua+fbms[:,:,nxs,nxs]*dfbdua,axis=0)      #ok
        iBeg += blkSize
    print (' %d derivative time %.4f\r'%(nRef,time.time()-time0))
        #loop over atoms - each dict entry is list of derivatives for all the reflections
    for i in range(len(Mdata)):
        dFdvDict[pfx+'Afrac:'+str(i)] = dFdfr.T[i]
        dFdvDict[pfx+'dAx:'+str(i)] = dFdx.T[0][i]
        dFdvDict[pfx+'dAy:'+str(i)] = dFdx.T[1][i]
        dFdvDict[pfx+'dAz:'+str(i)] = dFdx.T[2][i]
        dFdvDict[pfx+'AUiso:'+str(i)] = dFdui.T[i]
        dFdvDict[pfx+'AU11:'+str(i)] = dFdua.T[0][i]
        dFdvDict[pfx+'AU22:'+str(i)] = dFdua.T[1][i]
        dFdvDict[pfx+'AU33:'+str(i)] = dFdua.T[2][i]
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*dFdua.T[3][i]
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*dFdua.T[4][i]
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*dFdua.T[5][i]
    return dFdvDict
        
def StructureFactorDervTw2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    '''Compute structure factor derivatives on blocks of reflections - for twins only
    faster than StructureFactorDervTw
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,it,d,...
        'FF' dict of form factors - filled in below
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict parmDict:
    
    :returns: dict dFdvDict: dictionary of derivatives
    '''
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    TwDict = refDict.get('TwDict',{})           
    NTL = calcControls[phfx+'NTL']
    NM = calcControls[phfx+'TwinNMN']+1
    TwinLaw = calcControls[phfx+'TwinLaw']
    TwinFr = np.array([parmDict[phfx+'TwinFr:'+str(i)] for i in range(len(TwinLaw))])
    TwinInv = list(np.where(calcControls[phfx+'TwinInv'],-1,1))
    nTwin = len(TwinLaw)        
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return {}
    mSize = len(Mdata)
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    dFdvDict = {}
    dFdfr = np.zeros((nRef,nTwin,mSize))
    dFdx = np.zeros((nRef,nTwin,mSize,3))
    dFdui = np.zeros((nRef,nTwin,mSize))
    dFdua = np.zeros((nRef,nTwin,mSize,6))
    dFdbab = np.zeros((nRef,nTwin,2))
    dFdtw = np.zeros((nRef,nTwin))
    time0 = time.time()
#reflection processing begins here - big arrays!
    iBeg = 0
    blkSize = 16       #no. of reflections in a block - optimized for speed
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:3]
        H = np.inner(H.T,TwinLaw)   #array(3,nTwins)
        TwMask = np.any(H,axis=-1)
        for ir in range(blkSize):
            iref = ir+iBeg
            if iref in TwDict:
                for i in TwDict[iref]:
                    for n in range(NTL):
                        H[ir][i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
        TwMask = np.any(H,axis=-1)
        SQ = 1./(2.*refl.T[4])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,len(SGT)*len(TwinLaw),axis=0)
            FPP = np.repeat(FPP.T,len(SGT)*len(TwinLaw),axis=0)
        dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
        Bab = np.repeat(parmDict[phfx+'BabA']*dBabdA,len(SGT)*nTwin)
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,len(SGT)*len(TwinLaw),axis=0)
        Uniq = np.inner(H,SGMT)             # (nTwin,nSGOp,3)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T).T+Phi.T).T
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(SGT)
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT)*nTwin,axis=1)
        HbH = -np.sum(Uniq.T*np.swapaxes(np.inner(bij,Uniq),2,-1),axis=1)
        Hij = np.array([Mast*np.multiply.outer(U,U) for U in np.reshape(Uniq,(-1,3))])
        Hij = np.reshape(np.array([G2lat.UijtoU6(uij) for uij in Hij]),(-1,nTwin,len(SGT),6))
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = (np.reshape(Tiso,Tuij.shape)*Tuij).T*Mdata*Fdata/len(SGMT)
        fot = np.reshape(((FF+FP).T-Bab).T,cosp.shape)*Tcorr
        fotp = FPP*Tcorr        
        if 'T' in calcControls[hfx+'histType']: #fa,fb are 2 X blkSize X nTwin X nOps x nAtoms
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-np.reshape(FPP,sinp.shape)*sinp*Tcorr])
            fb = np.array([np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr,np.reshape(FPP,cosp.shape)*cosp*Tcorr])
        else:
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-FPP*sinp*Tcorr])
            fb = np.array([np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr,FPP*cosp*Tcorr])
        fas = np.sum(np.sum(fa,axis=-1),axis=-1)      #real sum over atoms & unique hkl array(2,nTwins)
        fbs = np.sum(np.sum(fb,axis=-1),axis=-1)      #imag sum over atoms & uniq hkl
        if SGData['SGInv']: #centrosymmetric; B=0
            fbs[0] *= 0.
            fas[1] *= 0.
        fax = np.array([-fot*sinp,-fotp*cosp])   #positions array(2,nRef,ntwi,nEqv,nAtoms)
        fbx = np.array([fot*cosp,-fotp*sinp])
        #sum below is over Uniq 
        dfadfr = np.sum(np.sum(fa/occ,axis=-2),axis=0)        #array(2,nRef,ntwin,nAtom) Fdata != 0 avoids /0. problem 
        dfadba = np.sum(-cosp*Tcorr[:,nxs],axis=1)
        dfadui = np.sum(np.sum(-SQfactor[nxs,:,nxs,nxs,nxs]*fa,axis=-2),axis=0)            
        dfadx = np.sum(np.sum(twopi*Uniq[nxs,:,:,:,nxs,:]*fax[:,:,:,:,:,nxs],axis=-3),axis=0) # nRef x nTwin x nAtoms x xyz; sum on ops & A,A'
        dfadua = np.sum(np.sum(-Hij[nxs,:,:,:,nxs,:]*fa[:,:,:,:,:,nxs],axis=-3),axis=0) 
        if not SGData['SGInv']:
            dfbdfr = np.sum(np.sum(fb/occ,axis=-2),axis=0)        #Fdata != 0 avoids /0. problem
            dfadba /= 2.
#            dfbdba = np.sum(-sinp*Tcorr[:,nxs],axis=1)/2.
            dfbdui = np.sum(np.sum(-SQfactor[nxs,:,nxs,nxs,nxs]*fb,axis=-2),axis=0)
            dfbdx = np.sum(np.sum(twopi*Uniq[nxs,:,:,:,nxs,:]*fbx[:,:,:,:,:,nxs],axis=-3),axis=0)
            dfbdua = np.sum(np.sum(-Hij[nxs,:,:,:,nxs,:]*fb[:,:,:,:,:,nxs],axis=-3),axis=0)
        else:
            dfbdfr = np.zeros_like(dfadfr)
            dfbdx = np.zeros_like(dfadx)
            dfbdui = np.zeros_like(dfadui)
            dfbdua = np.zeros_like(dfadua)
#            dfbdba = np.zeros_like(dfadba)
        SA = fas[0]+fas[1]
        SB = fbs[0]+fbs[1]
        dFdfr[iBeg:iFin] = ((2.*TwMask*SA)[:,:,nxs]*dfadfr+(2.*TwMask*SB)[:,:,nxs]*dfbdfr)*Mdata[nxs,nxs,:]/len(SGMT)
        dFdx[iBeg:iFin] = (2.*TwMask*SA)[:,:,nxs,nxs]*dfadx+(2.*TwMask*SB)[:,:,nxs,nxs]*dfbdx
        dFdui[iBeg:iFin] = (2.*TwMask*SA)[:,:,nxs]*dfadui+(2.*TwMask*SB)[:,:,nxs]*dfbdui
        dFdua[iBeg:iFin] = (2.*TwMask*SA)[:,:,nxs,nxs]*dfadua+(2.*TwMask*SB)[:,:,nxs,nxs]*dfbdua
        if SGData['SGInv']: #centrosymmetric; B=0
            dFdtw[iBeg:iFin] = np.sum(TwMask[nxs,:]*fas,axis=0)**2
        else:                
            dFdtw[iBeg:iFin] = np.sum(TwMask[nxs,:]*fas,axis=0)**2+np.sum(TwMask[nxs,:]*fbs,axis=0)**2
#        dFdbab[iBeg:iFin] = fas[0,:,nxs]*np.array([np.sum(dfadba*dBabdA),np.sum(-dfadba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T+ \
#            fbs[0,:,nxs]*np.array([np.sum(dfbdba*dBabdA),np.sum(-dfbdba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
        iBeg += blkSize
    print (' %d derivative time %.4f\r'%(len(refDict['RefList']),time.time()-time0))
    #loop over atoms - each dict entry is list of derivatives for all the reflections
    for i in range(len(Mdata)):     #these all OK
        dFdvDict[pfx+'Afrac:'+str(i)] = np.sum(dFdfr.T[i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'dAx:'+str(i)] = np.sum(dFdx.T[0][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'dAy:'+str(i)] = np.sum(dFdx.T[1][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'dAz:'+str(i)] = np.sum(dFdx.T[2][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AUiso:'+str(i)] = np.sum(dFdui.T[i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU11:'+str(i)] = np.sum(dFdua.T[0][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU22:'+str(i)] = np.sum(dFdua.T[1][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU33:'+str(i)] = np.sum(dFdua.T[2][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*np.sum(dFdua.T[3][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*np.sum(dFdua.T[4][i]*TwinFr[:,nxs],axis=0)
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*np.sum(dFdua.T[5][i]*TwinFr[:,nxs],axis=0)
    dFdvDict[phfx+'BabA'] = dFdbab.T[0]
    dFdvDict[phfx+'BabU'] = dFdbab.T[1]
    for i in range(nTwin):
        dFdvDict[phfx+'TwinFr:'+str(i)] = dFdtw.T[i]
    return dFdvDict
    
def SStructureFactor(refDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict):
    ''' 
    Compute super structure factors for all h,k,l,m for phase - no twins
    puts the result, F^2, in each ref[9] in refList
    works on blocks of 32 reflections for speed
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,it,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:

    '''
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    GS = G/np.outer(ast,ast)
    uAmat,uBmat = G2lat.Gmat2AB(GS)
    Mast = twopisq*np.multiply.outer(ast,ast)    
    SGInv = SGData['SGInv']
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    Nops = len(SGMT)*(1+SGData['SGInv'])
    if SGData['SGGray']:
        Nops *= 2
    SSGMT = np.array([ops[0].T for ops in SSGData['SSGOps']])
    SSGT = np.array([ops[1] for ops in SSGData['SSGOps']])
    SSCen = SSGData['SSGCen']
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    MFtables = calcControls['MFtables']
    Amat,Bmat = G2lat.Gmat2AB(G)
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return
    waveTypes,FSSdata,XSSdata,USSdata,MSSdata = GetAtomSSFXU(pfx,calcControls,parmDict)
    ngl,nWaves,Fmod,Xmod,Umod,Mmod,glTau,glWt = G2mth.makeWaves(waveTypes,FSSdata,XSSdata,USSdata,MSSdata,Mast)      #NB: Mmod is ReIm,Mxyz,Ntau,Natm
    modQ = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])

    if parmDict[pfx+'isMag']:       #This part correct for making modulated mag moments on equiv atoms - Mmod matched drawing & Bilbao drawings
    
        mXYZ = np.array([[XYZ[0] for XYZ in list(G2spc.GenAtom(xyz,SGData,All=True,Move=True))] for xyz in (Xdata+dXdata).T]) #Natn,Nop,xyz
        if SGData['SGGray']:
            mXYZ = np.hstack((mXYZ,mXYZ))

        MmodA,MmodB = G2mth.MagMod(glTau,mXYZ,modQ,MSSdata,SGData,SSGData)  #Ntau,Nops,Natm,Mxyz cos,sim parts sum matches drawing
        
        if not SGData['SGGray']:    #for fixed Mx,My,Mz
            GSdata = np.inner(Gdata.T,np.swapaxes(SGMT,1,2))  #apply sym. ops.--> Natm,Nops,Nxyz
            if SGData['SGInv'] and not SGData['SGFixed']:   #inversion if any
                GSdata = np.hstack((GSdata,-GSdata))      
            GSdata = np.hstack([GSdata for cen in SSCen])        #dup over cell centering - Natm,Nops,Mxyz
            GSdata = SGData['MagMom'][nxs,:,nxs]*GSdata   #flip vectors according to spin flip * det(opM)
            GSdata = np.swapaxes(GSdata,0,1)    #Nop,Natm,Mxyz
        
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata)).T
    bij = Mast*Uij
    blkSize = 48       #no. of reflections in a block
    nRef = refDict['RefList'].shape[0]
    SQ = 1./(2.*refDict['RefList'].T[5])**2
    if 'N' in calcControls[hfx+'histType']:
        dat = G2el.getBLvalues(BLtables)
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.ones((nRef,len(dat)))*list(dat.values())
        refDict['FF']['MF'] = np.zeros((nRef,len(dat)))
        for iel,El in enumerate(refDict['FF']['El']):
            if El in MFtables:
                refDict['FF']['MF'].T[iel] = G2el.MagScatFac(MFtables[El],SQ)
    else:
        dat = G2el.getFFvalues(FFtables,0.)
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.zeros((nRef,len(dat)))
        for iel,El in enumerate(refDict['FF']['El']):
            refDict['FF']['FF'].T[iel] = G2el.ScatFac(FFtables[El],SQ)
#reflection processing begins here - big arrays!
    iBeg = 0
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        mRef= iFin-iBeg
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl.T[:4]                          #array(blkSize,4)
        HP = H[:3]+modQ[:,nxs]*H[3:]            #projected hklm to hkl
        SQ = 1./(2.*refl.T[5])**2               #array(blkSize)
        SQfactor = 4.0*SQ*twopisq               #ditto prev.
        Uniq = np.inner(H.T,SSGMT)
        UniqP = np.inner(HP.T,SGMT)
        Phi = np.inner(H.T,SSGT)
        if SGInv and not SGData['SGFixed']:   #if centro - expand HKL sets
            Uniq = np.hstack((Uniq,-Uniq))
            Phi = np.hstack((Phi,-Phi))
            UniqP = np.hstack((UniqP,-UniqP))
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,Uniq.shape[1],axis=0)
            FPP = np.repeat(FPP.T,Uniq.shape[1],axis=0)
        Bab = 0.
        if phfx+'BabA' in parmDict:
            Bab = np.repeat(parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor),Uniq.shape[1])
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,Uniq.shape[1],axis=0)
        phase = twopi*(np.inner(Uniq[:,:,:3],(dXdata.T+Xdata.T))-Phi[:,:,nxs])
        phase = np.hstack([phase for cen in SSCen])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),Uniq.shape[1],axis=1).T
        HbH = -np.sum(UniqP[:,:,nxs,:]*np.inner(UniqP[:,:,:],bij),axis=-1)  #use hklt proj to hkl
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/Uniq.shape[1]  #refBlk x ops x atoms

        if 'N' in calcControls[hfx+'histType'] and parmDict[pfx+'isMag']:            
            
            phasem = twopi*np.inner(mXYZ,HP.T).T    #2pi(Q.r)
            cosm = np.cos(phasem)                   #Nref,nops,natm
            sinm = np.sin(phasem)
            MF = refDict['FF']['MF'][iBeg:iFin].T[Tindx].T   #Nref,Natm
            TMcorr = 0.539*(np.reshape(Tiso,Tuij.shape)*Tuij)[:,0,:]*Mdata*Fdata*MF/(2*Nops)     #Nref,Natm
            HM = np.inner(uAmat,HP.T)                    #put into cartesian space X||H,Z||H*L; uAmat better than uBmat
            eM = (HM/np.sqrt(np.sum(HM**2,axis=0))).T    # normalize  HP  Nref,hkl=Unit vectors || Q

            if not SGData['SGGray']:     #correct -fixed Mx,My,Mz contribution              
                fam0 = TMcorr[:,nxs,:,nxs]*GSdata[nxs,:,:,:]*cosm[:,:,:,nxs]    #Nref,Nops,Natm,Mxyz
                fbm0 = TMcorr[:,nxs,:,nxs]*GSdata[nxs,:,:,:]*sinm[:,:,:,nxs]
#  calc mag. structure factors; Nref,Ntau,Nops,Natm,Mxyz                            
            fams = TMcorr[:,nxs,nxs,:,nxs]*np.array([np.where(H[3,i]!=0,(MmodA*cosm[i,nxs,:,:,nxs]+    
                H[3,i]*MmodB*sinm[i,nxs,:,:,nxs]),0.) for i in range(mRef)])          #Nref,Ntau,Nops,Natm,Mxyz
                        
            fbms = TMcorr[:,nxs,nxs,:,nxs]*np.array([np.where(H[3,i]!=0,(MmodA*sinm[i,nxs,:,:,nxs]+    
                H[3,i]*MmodB*cosm[i,nxs,:,:,nxs]),0.) for i in range(mRef)])          #Nref,Ntau,Nops,Natm,Mxyz

            if not SGData['SGGray']:            
                fams += fam0[:,nxs,:,:,:]
                fbms += fbm0[:,nxs,:,:,:]
                                
#sum ops & atms                                
            fasm = np.sum(np.sum(fams,axis=-2),axis=-2)    #Nref,Mxyz; sum ops & atoms
            fbsm = np.sum(np.sum(fbms,axis=-2),axis=-2)
#put into cartesian space
            facm = np.inner(fasm,uBmat.T)
            fbcm = np.inner(fbsm,uBmat.T)
#form e.F dot product
            eDotFa = np.sum(eM[:,nxs,:]*facm,axis=-1)    #Nref,Ntau        
            eDotFb = np.sum(eM[:,nxs,:]*fbcm,axis=-1)
#intensity Halpern & Johnson
            fass = np.sum((facm-eM[:,nxs,:]*eDotFa[:,:,nxs])**2,axis=-1)
            fbss = np.sum((fbcm-eM[:,nxs,:]*eDotFb[:,:,nxs])**2,axis=-1)
## #do integration            
            fas = np.sum(fass*glWt[nxs,:],axis=1)
            fbs = np.sum(fbss*glWt[nxs,:],axis=1)
            
            refl.T[10] = fas+fbs   #Sum(fams**2,Mxyz) Re + Im
            refl.T[11] = atan2d(fbs,fas)

        else:
            GfpuA = G2mth.Modulation(Uniq,UniqP,nWaves,Fmod,Xmod,Umod,glTau,glWt) #2 x refBlk x sym X atoms
            if 'T' in calcControls[hfx+'histType']:
                fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-np.reshape(Flack*FPP,sinp.shape)*sinp*Tcorr])
                fb = np.array([np.reshape(Flack*FPP,cosp.shape)*cosp*Tcorr,np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr])
            else:
                fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-Flack*FPP*sinp*Tcorr])
                fb = np.array([Flack*FPP*cosp*Tcorr,np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr])
            fag = fa*GfpuA[0]-fb*GfpuA[1]   #real; 2 x refBlk x sym x atoms
            fbg = fb*GfpuA[0]+fa*GfpuA[1]
            fas = np.sum(np.sum(fag,axis=-1),axis=-1)   #2 x refBlk; sum sym & atoms
            fbs = np.sum(np.sum(fbg,axis=-1),axis=-1)
            
            refl.T[10] = np.sum(fas,axis=0)**2+np.sum(fbs,axis=0)**2    #square of sums
            refl.T[11] = atan2d(fbs[0],fas[0])  #use only tau=0; ignore f' & f"
        if 'P' not in calcControls[hfx+'histType']:
            refl.T[8] = np.copy(refl.T[10])                
        iBeg += blkSize
    return copy.deepcopy(refDict['RefList'])

def SStructureFactorTw(refDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict):
    ''' 
    Compute super structure factors for all h,k,l,m for phase - twins only
    puts the result, F^2, in each ref[8+im] in refList
    works on blocks of 32 reflections for speed
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,it,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:

    '''
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)    
    SGInv = SGData['SGInv']
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SSGMT = np.array([ops[0].T for ops in SSGData['SSGOps']])
    SSGT = np.array([ops[1] for ops in SSGData['SSGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    MFtables = calcControls['MFtables']
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
    TwinLaw = np.array([[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],])    #4D?
    TwDict = refDict.get('TwDict',{})           
    if 'S' in calcControls[hfx+'histType']:
        NTL = calcControls[phfx+'NTL']
        NM = calcControls[phfx+'TwinNMN']+1
        TwinLaw = calcControls[phfx+'TwinLaw']  #this'll have to be 4D also...
        TwinFr = np.array([parmDict[phfx+'TwinFr:'+str(i)] for i in range(len(TwinLaw))])
        TwinInv = list(np.where(calcControls[phfx+'TwinInv'],-1,1))
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return
    waveTypes,FSSdata,XSSdata,USSdata,MSSdata = GetAtomSSFXU(pfx,calcControls,parmDict)
    ngl,nWaves,Fmod,Xmod,Umod,Mmod,glTau,glWt = G2mth.makeWaves(waveTypes,FSSdata,XSSdata,USSdata,Mast)
    modQ = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata)).T
    bij = Mast*Uij
    blkSize = 32       #no. of reflections in a block
    nRef = refDict['RefList'].shape[0]
    if not len(refDict['FF']):                #no form factors - 1st time thru StructureFactor
        SQ = 1./(2.*refDict['RefList'].T[5])**2
        if 'N' in calcControls[hfx+'histType']:
            dat = G2el.getBLvalues(BLtables)
            refDict['FF']['El'] = list(dat.keys())
            refDict['FF']['FF'] = np.ones((nRef,len(dat)))*list(dat.values())
            refDict['FF']['MF'] = np.zeros((nRef,len(dat)))
            for iel,El in enumerate(refDict['FF']['El']):
                if El in MFtables:
                    refDict['FF']['MF'].T[iel] = G2el.MagScatFac(MFtables[El],SQ)
        else:
            dat = G2el.getFFvalues(FFtables,0.)
            refDict['FF']['El'] = list(dat.keys())
            refDict['FF']['FF'] = np.zeros((nRef,len(dat)))
            for iel,El in enumerate(refDict['FF']['El']):
                refDict['FF']['FF'].T[iel] = G2el.ScatFac(FFtables[El],SQ)
#    time0 = time.time()
#reflection processing begins here - big arrays!
    iBeg = 0
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]    #array(blkSize,nItems)
        H = refl[:,:4]                          #array(blkSize,4)
        H3 = refl[:,:3]
        HP = H[:,:3]+modQ[nxs,:]*H[:,3:]        #projected hklm to hkl
        HP = np.inner(HP,TwinLaw)             #array(blkSize,nTwins,4)
        H3 = np.inner(H3,TwinLaw)        
        TwMask = np.any(HP,axis=-1)
        if TwinLaw.shape[0] > 1 and TwDict: #need np.inner(TwinLaw[?],TwDict[iref][i])*TwinInv[i]
            for ir in range(blkSize):
                iref = ir+iBeg
                if iref in TwDict:
                    for i in TwDict[iref]:
                        for n in range(NTL):
                            HP[ir][i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
                            H3[ir][i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
            TwMask = np.any(HP,axis=-1)
        SQ = 1./(2.*refl.T[5])**2               #array(blkSize)
        SQfactor = 4.0*SQ*twopisq               #ditto prev.
        Uniq = np.inner(H,SSGMT)
        Uniq3 = np.inner(H3,SGMT)
        UniqP = np.inner(HP,SGMT)
        Phi = np.inner(H,SSGT)
        if SGInv:   #if centro - expand HKL sets
            Uniq = np.hstack((Uniq,-Uniq))
            Uniq3 = np.hstack((Uniq3,-Uniq3))
            Phi = np.hstack((Phi,-Phi))
            UniqP = np.hstack((UniqP,-UniqP))
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,Uniq.shape[1]*len(TwinLaw),axis=0)
            FPP = np.repeat(FPP.T,Uniq.shape[1]*len(TwinLaw),axis=0)
        Bab = np.repeat(parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor),Uniq.shape[1]*len(TwinLaw))
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,Uniq.shape[1]*len(TwinLaw),axis=0)
        phase = twopi*(np.inner(Uniq3,(dXdata.T+Xdata.T))-Phi[:,nxs,:,nxs])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),Uniq.shape[1]*len(TwinLaw),axis=1).T
        HbH = -np.sum(UniqP[:,:,:,nxs]*np.inner(UniqP[:,:,:],bij),axis=-1)  #use hklt proj to hkl
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/Uniq.shape[1]  #refBlk x ops x atoms
        if 'T' in calcControls[hfx+'histType']:
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-np.reshape(Flack*FPP,sinp.shape)*sinp*Tcorr])
            fb = np.array([np.reshape(Flack*FPP,cosp.shape)*cosp*Tcorr,np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr])
        else:
            fa = np.array([np.reshape(((FF+FP).T-Bab).T,cosp.shape)*cosp*Tcorr,-Flack*FPP*sinp*Tcorr])
            fb = np.array([Flack*FPP*cosp*Tcorr,np.reshape(((FF+FP).T-Bab).T,sinp.shape)*sinp*Tcorr])
        GfpuA = G2mth.ModulationTw(Uniq,UniqP,nWaves,Fmod,Xmod,Umod,glTau,glWt) #2 x refBlk x sym X atoms
        fag = fa*GfpuA[0]-fb*GfpuA[1]   #real; 2 x refBlk x sym x atoms
        fbg = fb*GfpuA[0]+fa*GfpuA[1]
        fas = np.sum(np.sum(fag,axis=-1),axis=-1)   #2 x refBlk; sum sym & atoms
        fbs = np.sum(np.sum(fbg,axis=-1),axis=-1)
        refl.T[10] = np.sum(fas[:,:,0],axis=0)**2+np.sum(fbs[:,:,0],axis=0)**2                  #FcT from primary twin element
        refl.T[8] = np.sum(TwinFr*np.sum(TwMask[nxs,:,:]*fas,axis=0)**2,axis=-1)+   \
            np.sum(TwinFr*np.sum(TwMask[nxs,:,:]*fbs,axis=0)**2,axis=-1)                 #Fc sum over twins
        refl.T[11] = atan2d(fbs[0].T[0],fas[0].T[0])  #ignore f' & f"
        iBeg += blkSize
#    print ('nRef %d time %.4f\r'%(nRef,time.time()-time0))

def SStructureFactorDerv(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict):
    ''' 
    Compute super structure factor derivatives for all h,k,l,m for phase - no twins
    Only Fourier component are done analytically here
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,it,d,...
        'FF' dict of form factors - filled in below
    :param int im: = 1 (could be eliminated)
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict SSGData: super space group info.
    :param dict calcControls:
    :param dict ParmDict:
    
    :returns: dict dFdvDict: dictionary of derivatives
    '''
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGInv = SGData['SGInv']
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SSGMT = np.array([ops[0].T for ops in SSGData['SSGOps']])
    SSGT = np.array([ops[1] for ops in SSGData['SSGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return {}
    mSize = len(Mdata)  #no. atoms
    waveTypes,FSSdata,XSSdata,USSdata,MSSdata = GetAtomSSFXU(pfx,calcControls,parmDict)
    ngl,nWaves,Fmod,Xmod,Umod,Mmod,glTau,glWt = G2mth.makeWaves(waveTypes,FSSdata,XSSdata,USSdata,MSSdata,Mast)
    waveShapes,SCtauF,SCtauX,SCtauU,UmodAB = G2mth.makeWavesDerv(ngl,waveTypes,FSSdata,XSSdata,USSdata,Mast)
    modQ = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata)).T
    bij = Mast*Uij
    if not len(refDict['FF']):
        if 'N' in calcControls[hfx+'histType']:
            dat = G2el.getBLvalues(BLtables)        #will need wave here for anom. neutron b's
        else:
            dat = G2el.getFFvalues(FFtables,0.)        
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.zeros((len(refDict['RefList']),len(dat)))
    dFdvDict = {}
    dFdfr = np.zeros((nRef,mSize))
    dFdx = np.zeros((nRef,mSize,3))
    dFdui = np.zeros((nRef,mSize))
    dFdua = np.zeros((nRef,mSize,6))
    dFdbab = np.zeros((nRef,2))
    dFdfl = np.zeros((nRef))
    dFdGf = np.zeros((nRef,mSize,FSSdata.shape[1],2))
    dFdGx = np.zeros((nRef,mSize,XSSdata.shape[1],6))
    dFdGu = np.zeros((nRef,mSize,USSdata.shape[1],12))
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
    time0 = time.time()
    nRef = len(refDict['RefList'])/100
    for iref,refl in enumerate(refDict['RefList']):
        if 'T' in calcControls[hfx+'histType']:
            FP,FPP = G2el.BlenResCW(Tdata,BLtables,refl.T[12+im])
        H = np.array(refl[:4])
        HP = H[:3]+modQ*H[3:]            #projected hklm to hkl
        SQ = 1./(2.*refl[4+im])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        Bab = 0.0
        if phfx+'BabA' in parmDict:
            dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
            Bab = parmDict[phfx+'BabA']*dBabdA
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = refDict['FF']['FF'][iref].T[Tindx]
        Uniq = np.inner(H,SSGMT)
        Phi = np.inner(H,SSGT)
        UniqP = np.inner(HP,SGMT)
        if SGInv:   #if centro - expand HKL sets
            Uniq = np.vstack((Uniq,-Uniq))
            Phi = np.hstack((Phi,-Phi))
            UniqP = np.vstack((UniqP,-UniqP))
        phase = twopi*(np.inner(Uniq[:,:3],(dXdata+Xdata).T)+Phi[:,nxs])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/Uniq.shape[0]
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),Uniq.shape[0],axis=1).T    #ops x atoms
        HbH = -np.sum(UniqP[:,nxs,:3]*np.inner(UniqP[:,:3],bij),axis=-1)  #ops x atoms
        Hij = np.array([Mast*np.multiply.outer(U[:3],U[:3]) for U in UniqP]) #atoms x 3x3
        Hij = np.array([G2lat.UijtoU6(uij) for uij in Hij])                     #atoms x 6
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)     #ops x atoms
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/Uniq.shape[0]  #ops x atoms
        fot = (FF+FP-Bab)*Tcorr     #ops x atoms
        fotp = FPP*Tcorr            #ops x atoms
        GfpuA = G2mth.Modulation(Uniq,UniqP,nWaves,Fmod,Xmod,Umod,glTau,glWt) #2 x sym X atoms
        dGdf,dGdx,dGdu = G2mth.ModulationDerv(Uniq,UniqP,Hij,nWaves,waveShapes,Fmod,Xmod,UmodAB,SCtauF,SCtauX,SCtauU,glTau,glWt)
        # GfpuA is 2 x ops x atoms
        # derivs are: ops x atoms x waves x 2,6,12, or 5 parms as [real,imag] parts
        fa = np.array([((FF+FP).T-Bab).T*cosp*Tcorr,-Flack*FPP*sinp*Tcorr]) # array(2,nEqv,nAtoms)
        fb = np.array([((FF+FP).T-Bab).T*sinp*Tcorr,Flack*FPP*cosp*Tcorr])  #or array(2,nEqv,nAtoms)
        fag = fa*GfpuA[0]-fb*GfpuA[1]
        fbg = fb*GfpuA[0]+fa*GfpuA[1]
        
        fas = np.sum(np.sum(fag,axis=1),axis=1)     # 2 x twin
        fbs = np.sum(np.sum(fbg,axis=1),axis=1)
        fax = np.array([-fot*sinp,-fotp*cosp])   #positions; 2 x ops x atoms
        fbx = np.array([fot*cosp,-fotp*sinp])
        fax = fax*GfpuA[0]-fbx*GfpuA[1]
        fbx = fbx*GfpuA[0]+fax*GfpuA[1]
        #sum below is over Uniq
        dfadfr = np.sum(fag/occ,axis=1)        #Fdata != 0 ever avoids /0. problem
        dfbdfr = np.sum(fbg/occ,axis=1)        #Fdata != 0 avoids /0. problem
        dfadba = np.sum(-cosp*Tcorr[:,nxs],axis=1)
        dfbdba = np.sum(-sinp*Tcorr[:,nxs],axis=1)
        dfadui = np.sum(-SQfactor*fag,axis=1)
        dfbdui = np.sum(-SQfactor*fbg,axis=1)
        dfadx = np.sum(twopi*Uniq[:,:3]*np.swapaxes(fax,-2,-1)[:,:,:,nxs],axis=-2)  #2 x nAtom x 3xyz; sum nOps
        dfbdx = np.sum(twopi*Uniq[:,:3]*np.swapaxes(fbx,-2,-1)[:,:,:,nxs],axis=-2)           
        dfadua = np.sum(-Hij*np.swapaxes(fag,-2,-1)[:,:,:,nxs],axis=-2)             #2 x nAtom x 6Uij; sum nOps
        dfbdua = np.sum(-Hij*np.swapaxes(fbg,-2,-1)[:,:,:,nxs],axis=-2)         #these are correct also for twins above
        # array(2,nAtom,nWave,2) & array(2,nAtom,nWave,6) & array(2,nAtom,nWave,12); sum on nOps
        dfadGf = np.sum(fa[:,:,:,nxs,nxs]*dGdf[0][nxs,:,:,:,:]-fb[:,:,:,nxs,nxs]*dGdf[1][nxs,:,:,:,:],axis=1)
        dfbdGf = np.sum(fb[:,:,:,nxs,nxs]*dGdf[0][nxs,:,:,:,:]+fa[:,:,:,nxs,nxs]*dGdf[1][nxs,:,:,:,:],axis=1)
        dfadGx = np.sum(fa[:,:,:,nxs,nxs]*dGdx[0][nxs,:,:,:,:]-fb[:,:,:,nxs,nxs]*dGdx[1][nxs,:,:,:,:],axis=1)
        dfbdGx = np.sum(fb[:,:,:,nxs,nxs]*dGdx[0][nxs,:,:,:,:]+fa[:,:,:,nxs,nxs]*dGdx[1][nxs,:,:,:,:],axis=1)
        dfadGu = np.sum(fa[:,:,:,nxs,nxs]*dGdu[0][nxs,:,:,:,:]-fb[:,:,:,nxs,nxs]*dGdu[1][nxs,:,:,:,:],axis=1)
        dfbdGu = np.sum(fb[:,:,:,nxs,nxs]*dGdu[0][nxs,:,:,:,:]+fa[:,:,:,nxs,nxs]*dGdu[1][nxs,:,:,:,:],axis=1)   
        if not SGData['SGInv']:   #Flack derivative
            dfadfl = np.sum(-FPP*Tcorr*sinp)
            dfbdfl = np.sum(FPP*Tcorr*cosp)
        else:
            dfadfl = 1.0
            dfbdfl = 1.0
        SA = fas[0]+fas[1]      #float = A+A' 
        SB = fbs[0]+fbs[1]      #float = B+B' 
        if 'P' in calcControls[hfx+'histType']: #checked perfect for centro & noncentro?
            dFdfl[iref] = -SA*dfadfl-SB*dfbdfl                  #array(nRef,) 
            dFdfr[iref] = 2.*(fas[0]*dfadfr[0]+fas[1]*dfadfr[1])*Mdata/len(Uniq)+   \
                2.*(fbs[0]*dfbdfr[0]-fbs[1]*dfbdfr[1])*Mdata/len(Uniq)
            dFdx[iref] = 2.*(fas[0]*dfadx[0]+fas[1]*dfadx[1])+  \
                2.*(fbs[0]*dfbdx[0]+fbs[1]*dfbdx[1])
            dFdui[iref] = 2.*(fas[0]*dfadui[0]+fas[1]*dfadui[1])+   \
                2.*(fbs[0]*dfbdui[0]-fbs[1]*dfbdui[1])
            dFdua[iref] = 2.*(fas[0]*dfadua[0]+fas[1]*dfadua[1])+   \
                2.*(fbs[0]*dfbdua[0]+fbs[1]*dfbdua[1])
            dFdGf[iref] = 2.*(fas[0]*dfadGf[0]+fas[1]*dfadGf[1])+  \
                2.*(fbs[0]*dfbdGf[0]+fbs[1]*dfbdGf[1])
            dFdGx[iref] = 2.*(fas[0]*dfadGx[0]+fas[1]*dfadGx[1])+  \
                2.*(fbs[0]*dfbdGx[0]-fbs[1]*dfbdGx[1])
            dFdGu[iref] = 2.*(fas[0]*dfadGu[0]+fas[1]*dfadGu[1])+  \
                2.*(fbs[0]*dfbdGu[0]+fbs[1]*dfbdGu[1])
        else:                       #OK, I think
            dFdfr[iref] = 2.*(SA*dfadfr[0]+SA*dfadfr[1]+SB*dfbdfr[0]+SB*dfbdfr[1])*Mdata/len(Uniq) #array(nRef,nAtom)
            dFdx[iref] = 2.*(SA*dfadx[0]+SA*dfadx[1]+SB*dfbdx[0]+SB*dfbdx[1])    #array(nRef,nAtom,3)
            dFdui[iref] = 2.*(SA*dfadui[0]+SA*dfadui[1]+SB*dfbdui[0]+SB*dfbdui[1])   #array(nRef,nAtom)
            dFdua[iref] = 2.*(SA*dfadua[0]+SA*dfadua[1]+SB*dfbdua[0]+SB*dfbdua[1])    #array(nRef,nAtom,6)
            dFdfl[iref] = -SA*dfadfl-SB*dfbdfl                  #array(nRef,) 
                           
            dFdGf[iref] = 2.*(SA*dfadGf[0]+SB*dfbdGf[1])      #array(nRef,natom,nwave,2)
            dFdGx[iref] = 2.*(SA*dfadGx[0]+SB*dfbdGx[1])      #array(nRef,natom,nwave,6)
            dFdGu[iref] = 2.*(SA*dfadGu[0]+SB*dfbdGu[1])      #array(nRef,natom,nwave,12)
        if phfx+'BabA' in parmDict:
            dFdbab[iref] = 2.*fas[0]*np.array([np.sum(dfadba*dBabdA),np.sum(-dfadba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T+ \
                2.*fbs[0]*np.array([np.sum(dfbdba*dBabdA),np.sum(-dfbdba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
        #loop over atoms - each dict entry is list of derivatives for all the reflections
        if not iref%100 :
            print (' %d derivative time %.4f\r'%(iref,time.time()-time0),end='')
    for i in range(len(Mdata)):     #loop over atoms
        dFdvDict[pfx+'Afrac:'+str(i)] = dFdfr.T[i]
        dFdvDict[pfx+'dAx:'+str(i)] = dFdx.T[0][i]
        dFdvDict[pfx+'dAy:'+str(i)] = dFdx.T[1][i]
        dFdvDict[pfx+'dAz:'+str(i)] = dFdx.T[2][i]
        dFdvDict[pfx+'AUiso:'+str(i)] = dFdui.T[i]
        dFdvDict[pfx+'AU11:'+str(i)] = dFdua.T[0][i]
        dFdvDict[pfx+'AU22:'+str(i)] = dFdua.T[1][i]
        dFdvDict[pfx+'AU33:'+str(i)] = dFdua.T[2][i]
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*dFdua.T[3][i]
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*dFdua.T[4][i]
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*dFdua.T[5][i]
        for j in range(FSSdata.shape[1]):        #loop over waves Fzero & Fwid?
            dFdvDict[pfx+'Fsin:'+str(i)+':'+str(j)] = dFdGf.T[0][j][i]
            dFdvDict[pfx+'Fcos:'+str(i)+':'+str(j)] = dFdGf.T[1][j][i]
        nx = 0
        if waveTypes[i] in ['Block','ZigZag']:
            nx = 1 
        for j in range(XSSdata.shape[1]-nx):       #loop over waves 
            dFdvDict[pfx+'Xsin:'+str(i)+':'+str(j+nx)] = dFdGx.T[0][j][i]
            dFdvDict[pfx+'Ysin:'+str(i)+':'+str(j+nx)] = dFdGx.T[1][j][i]
            dFdvDict[pfx+'Zsin:'+str(i)+':'+str(j+nx)] = dFdGx.T[2][j][i]
            dFdvDict[pfx+'Xcos:'+str(i)+':'+str(j+nx)] = dFdGx.T[3][j][i]
            dFdvDict[pfx+'Ycos:'+str(i)+':'+str(j+nx)] = dFdGx.T[4][j][i]
            dFdvDict[pfx+'Zcos:'+str(i)+':'+str(j+nx)] = dFdGx.T[5][j][i]
        for j in range(USSdata.shape[1]):       #loop over waves
            dFdvDict[pfx+'U11sin:'+str(i)+':'+str(j)] = dFdGu.T[0][j][i]
            dFdvDict[pfx+'U22sin:'+str(i)+':'+str(j)] = dFdGu.T[1][j][i]
            dFdvDict[pfx+'U33sin:'+str(i)+':'+str(j)] = dFdGu.T[2][j][i]
            dFdvDict[pfx+'U12sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[3][j][i]
            dFdvDict[pfx+'U13sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[4][j][i]
            dFdvDict[pfx+'U23sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[5][j][i]
            dFdvDict[pfx+'U11cos:'+str(i)+':'+str(j)] = dFdGu.T[6][j][i]
            dFdvDict[pfx+'U22cos:'+str(i)+':'+str(j)] = dFdGu.T[7][j][i]
            dFdvDict[pfx+'U33cos:'+str(i)+':'+str(j)] = dFdGu.T[8][j][i]
            dFdvDict[pfx+'U12cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[9][j][i]
            dFdvDict[pfx+'U13cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[10][j][i]
            dFdvDict[pfx+'U23cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[11][j][i]
            
    dFdvDict[phfx+'Flack'] = 4.*dFdfl.T
    dFdvDict[phfx+'BabA'] = dFdbab.T[0]
    dFdvDict[phfx+'BabU'] = dFdbab.T[1]
    return dFdvDict

def SStructureFactorDerv2(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict):
    '''
    Compute super structure factor derivatives for all h,k,l,m for phase - no twins
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,it,d,...
        'FF' dict of form factors - filled in below
    :param int im: = 1 (could be eliminated)
    :param np.array G:      reciprocal metric tensor
    :param str hfx:    histogram id string
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict SSGData: super space group info.
    :param dict calcControls:
    :param dict ParmDict:
    
    :returns: dict dFdvDict: dictionary of derivatives
    '''

    trefDict = copy.deepcopy(refDict)
    dM = 1.e-4
    dFdvDict = {}
    for parm in parmDict:
        if ':' not in parm:
            continue
        if parm.split(':')[2] in ['Tmin','Tmax','Xmax','Ymax','Zmax','Fzero','Fwid',
            'MXsin','MXcos','MYsin','MYcos','MZsin','MZcos','AMx','AMy','AMz',]:
            parmDict[parm] += dM
            prefList = SStructureFactor(trefDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
            parmDict[parm] -= 2*dM
            mrefList = SStructureFactor(trefDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
            parmDict[parm] += dM
            dFdvDict[parm] = (prefList[:,9+im]-mrefList[:,9+im])/(2.*dM)
    return dFdvDict
    
def SStructureFactorDervTw(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict):
    'Needs a doc string'
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGInv = SGData['SGInv']
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SSGMT = np.array([ops[0].T for ops in SSGData['SSGOps']])
    SSGT = np.array([ops[1] for ops in SSGData['SSGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    TwinLaw = np.array([[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],])
    TwDict = refDict.get('TwDict',{})           
    if 'S' in calcControls[hfx+'histType']:
        NTL = calcControls[phfx+'NTL']
        NM = calcControls[phfx+'TwinNMN']+1
        TwinLaw = calcControls[phfx+'TwinLaw']
        TwinInv = list(np.where(calcControls[phfx+'TwinInv'],-1,1))
    nTwin = len(TwinLaw)        
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata,Gdata = \
        GetAtomFXU(pfx,calcControls,parmDict)
    if not Xdata.size:          #no atoms in phase!
        return {}
    mSize = len(Mdata)  #no. atoms
    waveTypes,FSSdata,XSSdata,USSdata,MSSdata = GetAtomSSFXU(pfx,calcControls,parmDict)
    ngl,nWaves,Fmod,Xmod,Umod,Mmod,glTau,glWt = G2mth.makeWaves(waveTypes,FSSdata,XSSdata,USSdata,MSSdata,Mast)     #NB: Mmod is ReIm,Mxyz,Ntau,Natm
    waveShapes,SCtauF,SCtauX,SCtauU,UmodAB = G2mth.makeWavesDerv(ngl,waveTypes,FSSdata,XSSdata,USSdata,Mast)
    modQ = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType'] or 'NB' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata)).T
    bij = Mast*Uij
    if not len(refDict['FF']):
        if 'N' in calcControls[hfx+'histType']:
            dat = G2el.getBLvalues(BLtables)        #will need wave here for anom. neutron b's
        else:
            dat = G2el.getFFvalues(FFtables,0.)        
        refDict['FF']['El'] = list(dat.keys())
        refDict['FF']['FF'] = np.zeros((len(refDict['RefList']),len(dat)))
    dFdvDict = {}
    dFdfr = np.zeros((nRef,nTwin,mSize))
    dFdx = np.zeros((nRef,nTwin,mSize,3))
    dFdui = np.zeros((nRef,nTwin,mSize))
    dFdua = np.zeros((nRef,nTwin,mSize,6))
    dFdbab = np.zeros((nRef,nTwin,2))
    dFdtw = np.zeros((nRef,nTwin))
    dFdGf = np.zeros((nRef,nTwin,mSize,FSSdata.shape[1]))
    dFdGx = np.zeros((nRef,nTwin,mSize,XSSdata.shape[1],3))
    dFdGz = np.zeros((nRef,nTwin,mSize,5))
    dFdGu = np.zeros((nRef,nTwin,mSize,USSdata.shape[1],6))
    Flack = 1.0
    if not SGData['SGInv'] and 'S' in calcControls[hfx+'histType'] and phfx+'Flack' in parmDict:
        Flack = 1.-2.*parmDict[phfx+'Flack']
    time0 = time.time()
    nRef = len(refDict['RefList'])/100
    for iref,refl in enumerate(refDict['RefList']):
        if 'T' in calcControls[hfx+'histType']:
            FP,FPP = G2el.BlenResCW(Tdata,BLtables,refl.T[12+im])
        H = np.array(refl[:4])
        HP = H[:3]+modQ*H[3:]            #projected hklm to hkl
        H = np.inner(H.T,TwinLaw)   #maybe array(4,nTwins) or (4)
        TwMask = np.any(H,axis=-1)
        if TwinLaw.shape[0] > 1 and TwDict:
            if iref in TwDict:
                for i in TwDict[iref]:
                    for n in range(NTL):
                        H[i+n*NM] = np.inner(TwinLaw[n*NM],np.array(TwDict[iref][i])*TwinInv[i+n*NM])
            TwMask = np.any(H,axis=-1)
        SQ = 1./(2.*refl[4+im])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
        Bab = parmDict[phfx+'BabA']*dBabdA
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = refDict['FF']['FF'][iref].T[Tindx]
        Uniq = np.inner(H,SSGMT)
        Phi = np.inner(H,SSGT)
        UniqP = np.inner(HP,SGMT)
        if SGInv:   #if centro - expand HKL sets
            Uniq = np.vstack((Uniq,-Uniq))
            Phi = np.hstack((Phi,-Phi))
            UniqP = np.vstack((UniqP,-UniqP))
        phase = twopi*(np.inner(Uniq[:,:3],(dXdata+Xdata).T)+Phi[:,nxs])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/Uniq.shape[0]
        biso = -SQfactor*Uisodata[:,nxs]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),Uniq.shape[0]*len(TwinLaw),axis=1).T    #ops x atoms
        HbH = -np.sum(UniqP[:,nxs,:3]*np.inner(UniqP[:,:3],bij),axis=-1)  #ops x atoms
        Hij = np.array([Mast*np.multiply.outer(U[:3],U[:3]) for U in UniqP]) #atoms x 3x3
        Hij = np.squeeze(np.reshape(np.array([G2lat.UijtoU6(uij) for uij in Hij]),(nTwin,-1,6)))
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)     #ops x atoms
        Tcorr = np.reshape(Tiso,Tuij.shape)*Tuij*Mdata*Fdata/Uniq.shape[0]  #ops x atoms
        fot = (FF+FP-Bab)*Tcorr     #ops x atoms
        fotp = FPP*Tcorr            #ops x atoms
        GfpuA = G2mth.Modulation(Uniq,UniqP,nWaves,Fmod,Xmod,Umod,glTau,glWt) #2 x sym X atoms
        dGdf,dGdx,dGdu,dGdz = G2mth.ModulationDerv(Uniq,UniqP,Hij,nWaves,waveShapes,Fmod,Xmod,UmodAB,SCtauF,SCtauX,SCtauU,glTau,glWt)
        # GfpuA is 2 x ops x atoms
        # derivs are: ops x atoms x waves x 2,6,12, or 5 parms as [real,imag] parts
        fa = np.array([((FF+FP).T-Bab).T*cosp*Tcorr,-Flack*FPP*sinp*Tcorr]) # array(2,nTwin,nEqv,nAtoms)
        fb = np.array([((FF+FP).T-Bab).T*sinp*Tcorr,Flack*FPP*cosp*Tcorr])  #or array(2,nEqv,nAtoms)
        fag = fa*GfpuA[0]-fb*GfpuA[1]
        fbg = fb*GfpuA[0]+fa*GfpuA[1]
        
        fas = np.sum(np.sum(fag,axis=1),axis=1)     # 2 x twin
        fbs = np.sum(np.sum(fbg,axis=1),axis=1)
        fax = np.array([-fot*sinp,-fotp*cosp])   #positions; 2 x twin x ops x atoms
        fbx = np.array([fot*cosp,-fotp*sinp])
        fax = fax*GfpuA[0]-fbx*GfpuA[1]
        fbx = fbx*GfpuA[0]+fax*GfpuA[1]
        #sum below is over Uniq
        dfadfr = np.sum(fag/occ,axis=1)        #Fdata != 0 ever avoids /0. problem
        dfbdfr = np.sum(fbg/occ,axis=1)        #Fdata != 0 avoids /0. problem
        dfadba = np.sum(-cosp*Tcorr[:,nxs],axis=1)
        dfbdba = np.sum(-sinp*Tcorr[:,nxs],axis=1)
        dfadui = np.sum(-SQfactor*fag,axis=1)
        dfbdui = np.sum(-SQfactor*fbg,axis=1)
        dfadx = np.array([np.sum(twopi*Uniq[it,:,:3]*np.swapaxes(fax,-2,-1)[:,it,:,:,nxs],axis=-2) for it in range(nTwin)])
        dfbdx = np.array([np.sum(twopi*Uniq[it,:,:3]*np.swapaxes(fbx,-2,-1)[:,it,:,:,nxs],axis=-2) for it in range(nTwin)])           
        dfadua = np.array([np.sum(-Hij[it]*np.swapaxes(fag,-2,-1)[:,it,:,:,nxs],axis=-2) for it in range(nTwin)])
        dfbdua = np.array([np.sum(-Hij[it]*np.swapaxes(fbg,-2,-1)[:,it,:,:,nxs],axis=-2) for it in range(nTwin)])
        # array(2,nTwin,nAtom,3) & array(2,nTwin,nAtom,6) & array(2,nTwin,nAtom,12)
        dfadGf = np.sum(fa[:,it,:,:,nxs,nxs]*dGdf[0][nxs,nxs,:,:,:,:]-fb[:,it,:,:,nxs,nxs]*dGdf[1][nxs,nxs,:,:,:,:],axis=1)
        dfbdGf = np.sum(fb[:,it,:,:,nxs,nxs]*dGdf[0][nxs,nxs,:,:,:,:]+fa[:,it,:,:,nxs,nxs]*dGdf[1][nxs,nxs,:,:,:,:],axis=1)
        dfadGx = np.sum(fa[:,it,:,:,nxs,nxs]*dGdx[0][nxs,nxs,:,:,:,:]-fb[:,it,:,:,nxs,nxs]*dGdx[1][nxs,nxs,:,:,:,:],axis=1)
        dfbdGx = np.sum(fb[:,it,:,:,nxs,nxs]*dGdx[0][nxs,nxs,:,:,:,:]+fa[:,it,:,:,nxs,nxs]*dGdx[1][nxs,nxs,:,:,:,:],axis=1)
        dfadGz = np.sum(fa[:,it,:,0,nxs,nxs]*dGdz[0][nxs,nxs,:,:,:]-fb[:,it,:,0,nxs,nxs]*dGdz[1][nxs,nxs,:,:,:],axis=1)
        dfbdGz = np.sum(fb[:,it,:,0,nxs,nxs]*dGdz[0][nxs,nxs,:,:,:]+fa[:,it,:,0,nxs,nxs]*dGdz[1][nxs,nxs,:,:,:],axis=1)
        dfadGu = np.sum(fa[:,it,:,:,nxs,nxs]*dGdu[0][nxs,nxs,:,:,:,:]-fb[:,it,:,:,nxs,nxs]*dGdu[1][nxs,nxs,:,:,:,:],axis=1)
        dfbdGu = np.sum(fb[:,it,:,:,nxs,nxs]*dGdu[0][nxs,nxs,:,:,:,:]+fa[:,it,:,:,nxs,nxs]*dGdu[1][nxs,nxs,:,:,:,:],axis=1)
#        GSASIIpath.IPyBreak()
        #NB: the above have been checked against PA(1:10,1:2) in strfctr.for for Al2O3!    
        SA = fas[0]+fas[1]      #float = A+A' (might be array[nTwin])
        SB = fbs[0]+fbs[1]      #float = B+B' (might be array[nTwin])
        dFdfr[iref] = [2.*TwMask[it]*(SA[it]*dfadfr[0][it]+SA[it]*dfadfr[1][it]+SB[it]*dfbdfr[0][it]+SB[it]*dfbdfr[1][it])*Mdata/len(Uniq[it]) for it in range(nTwin)]
        dFdx[iref] = [2.*TwMask[it]*(SA[it]*dfadx[it][0]+SA[it]*dfadx[it][1]+SB[it]*dfbdx[it][0]+SB[it]*dfbdx[it][1]) for it in range(nTwin)]
        dFdui[iref] = [2.*TwMask[it]*(SA[it]*dfadui[it][0]+SA[it]*dfadui[it][1]+SB[it]*dfbdui[it][0]+SB[it]*dfbdui[it][1]) for it in range(nTwin)]
        dFdua[iref] = [2.*TwMask[it]*(SA[it]*dfadua[it][0]+SA[it]*dfadua[it][1]+SB[it]*dfbdua[it][0]+SB[it]*dfbdua[it][1]) for it in range(nTwin)]
        dFdtw[iref] = np.sum(TwMask*fas,axis=0)**2+np.sum(TwMask*fbs,axis=0)**2

        dFdGf[iref] = [2.*TwMask[it]*(SA[it]*dfadGf[1]+SB[it]*dfbdGf[1]) for it in range(nTwin)]
        dFdGx[iref] = [2.*TwMask[it]*(SA[it]*dfadGx[1]+SB[it]*dfbdGx[1]) for it in range(nTwin)]
        dFdGz[iref] = [2.*TwMask[it]*(SA[it]*dfadGz[1]+SB[it]*dfbdGz[1]) for it in range(nTwin)]
        dFdGu[iref] = [2.*TwMask[it]*(SA[it]*dfadGu[1]+SB[it]*dfbdGu[1]) for it in range(nTwin)]                
#            GSASIIpath.IPyBreak()
        dFdbab[iref] = 2.*fas[0]*np.array([np.sum(dfadba*dBabdA),np.sum(-dfadba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T+ \
            2.*fbs[0]*np.array([np.sum(dfbdba*dBabdA),np.sum(-dfbdba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
        #loop over atoms - each dict entry is list of derivatives for all the reflections
        if not iref%100 :
            print (' %d derivative time %.4f\r'%(iref,time.time()-time0),end='')
    for i in range(len(Mdata)):     #loop over atoms
        dFdvDict[pfx+'Afrac:'+str(i)] = dFdfr.T[i]
        dFdvDict[pfx+'dAx:'+str(i)] = dFdx.T[0][i]
        dFdvDict[pfx+'dAy:'+str(i)] = dFdx.T[1][i]
        dFdvDict[pfx+'dAz:'+str(i)] = dFdx.T[2][i]
        dFdvDict[pfx+'AUiso:'+str(i)] = dFdui.T[i]
        dFdvDict[pfx+'AU11:'+str(i)] = dFdua.T[0][i]
        dFdvDict[pfx+'AU22:'+str(i)] = dFdua.T[1][i]
        dFdvDict[pfx+'AU33:'+str(i)] = dFdua.T[2][i]
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*dFdua.T[3][i]
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*dFdua.T[4][i]
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*dFdua.T[5][i]
        for j in range(FSSdata.shape[1]):        #loop over waves Fzero & Fwid?
            dFdvDict[pfx+'Fsin:'+str(i)+':'+str(j)] = dFdGf.T[0][j][i]
            dFdvDict[pfx+'Fcos:'+str(i)+':'+str(j)] = dFdGf.T[1][j][i]
        nx = 0
        if waveTypes[i] in ['Block','ZigZag']:
            nx = 1 
            dFdvDict[pfx+'Tmin:'+str(i)+':0'] = dFdGz.T[0][i]   #ZigZag/Block waves (if any)
            dFdvDict[pfx+'Tmax:'+str(i)+':0'] = dFdGz.T[1][i]
            dFdvDict[pfx+'Xmax:'+str(i)+':0'] = dFdGz.T[2][i]
            dFdvDict[pfx+'Ymax:'+str(i)+':0'] = dFdGz.T[3][i]
            dFdvDict[pfx+'Zmax:'+str(i)+':0'] = dFdGz.T[4][i]
        for j in range(XSSdata.shape[1]-nx):       #loop over waves 
            dFdvDict[pfx+'Xsin:'+str(i)+':'+str(j+nx)] = dFdGx.T[0][j][i]
            dFdvDict[pfx+'Ysin:'+str(i)+':'+str(j+nx)] = dFdGx.T[1][j][i]
            dFdvDict[pfx+'Zsin:'+str(i)+':'+str(j+nx)] = dFdGx.T[2][j][i]
            dFdvDict[pfx+'Xcos:'+str(i)+':'+str(j+nx)] = dFdGx.T[3][j][i]
            dFdvDict[pfx+'Ycos:'+str(i)+':'+str(j+nx)] = dFdGx.T[4][j][i]
            dFdvDict[pfx+'Zcos:'+str(i)+':'+str(j+nx)] = dFdGx.T[5][j][i]
        for j in range(USSdata.shape[1]):       #loop over waves
            dFdvDict[pfx+'U11sin:'+str(i)+':'+str(j)] = dFdGu.T[0][j][i]
            dFdvDict[pfx+'U22sin:'+str(i)+':'+str(j)] = dFdGu.T[1][j][i]
            dFdvDict[pfx+'U33sin:'+str(i)+':'+str(j)] = dFdGu.T[2][j][i]
            dFdvDict[pfx+'U12sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[3][j][i]
            dFdvDict[pfx+'U13sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[4][j][i]
            dFdvDict[pfx+'U23sin:'+str(i)+':'+str(j)] = 2.*dFdGu.T[5][j][i]
            dFdvDict[pfx+'U11cos:'+str(i)+':'+str(j)] = dFdGu.T[6][j][i]
            dFdvDict[pfx+'U22cos:'+str(i)+':'+str(j)] = dFdGu.T[7][j][i]
            dFdvDict[pfx+'U33cos:'+str(i)+':'+str(j)] = dFdGu.T[8][j][i]
            dFdvDict[pfx+'U12cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[9][j][i]
            dFdvDict[pfx+'U13cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[10][j][i]
            dFdvDict[pfx+'U23cos:'+str(i)+':'+str(j)] = 2.*dFdGu.T[11][j][i]
            
#        GSASIIpath.IPyBreak()
    dFdvDict[phfx+'BabA'] = dFdbab.T[0]
    dFdvDict[phfx+'BabU'] = dFdbab.T[1]
    return dFdvDict
    
def SCExtinction(ref,im,phfx,hfx,pfx,calcControls,parmDict,varyList):
    ''' Single crystal extinction function; returns extinction & derivative
    '''
    extCor = 1.0
    dervDict = {}
    dervCor = 1.0
    if calcControls[phfx+'EType'] != 'None':
        SQ = 1/(4.*ref[4+im]**2)
        if 'C' in parmDict[hfx+'Type']:            
            cos2T = 1.0-2.*SQ*parmDict[hfx+'Lam']**2           #cos(2theta)
        else:   #'T'
            cos2T = 1.0-2.*SQ*ref[12+im]**2                       #cos(2theta)            
        if 'SXC' in parmDict[hfx+'Type']:
            AV = 7.9406e5/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            P12 = (calcControls[phfx+'Cos2TM']+cos2T**4)/(calcControls[phfx+'Cos2TM']+cos2T**2)
            PLZ = AV*P12*ref[9+im]*parmDict[hfx+'Lam']**2
        elif 'SNT' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = SQ
            PLZ = AV*ref[9+im]*ref[12+im]**2
        elif 'SNC' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            PLZ = AV*ref[9+im]*parmDict[hfx+'Lam']**2
            
        if 'Primary' in calcControls[phfx+'EType']:
            PLZ *= 1.5
        else:
            if 'C' in parmDict[hfx+'Type']:
                PLZ *= calcControls[phfx+'Tbar']
            else: #'T'
                PLZ *= ref[13+im]      #t-bar
        if 'Primary' in calcControls[phfx+'EType']:
            PLZ *= 1.5
            PSIG = parmDict[phfx+'Ep']
        elif 'I & II' in calcControls[phfx+'EType']:
            PSIG = parmDict[phfx+'Eg']/np.sqrt(1.+(parmDict[phfx+'Es']*PL/parmDict[phfx+'Eg'])**2)
        elif 'Type II' in calcControls[phfx+'EType']:
            PSIG = parmDict[phfx+'Es']
        else:       # 'Secondary Type I'
            PSIG = parmDict[phfx+'Eg']/PL
            
        AG = 0.58+0.48*cos2T+0.24*cos2T**2
        AL = 0.025+0.285*cos2T
        BG = 0.02-0.025*cos2T
        BL = 0.15-0.2*(0.75-cos2T)**2
        if cos2T < 0.:
            BL = -0.45*cos2T
        CG = 2.
        CL = 2.
        PF = PLZ*PSIG
        
        if 'Gaussian' in calcControls[phfx+'EApprox']:
            PF4 = 1.+CG*PF+AG*PF**2/(1.+BG*PF)
            extCor = np.sqrt(PF4)
            PF3 = 0.5*(CG+2.*AG*PF/(1.+BG*PF)-AG*PF**2*BG/(1.+BG*PF)**2)/(PF4*extCor)
        else:
            PF4 = 1.+CL*PF+AL*PF**2/(1.+BL*PF)
            extCor = np.sqrt(PF4)
            PF3 = 0.5*(CL+2.*AL*PF/(1.+BL*PF)-AL*PF**2*BL/(1.+BL*PF)**2)/(PF4*extCor)

        dervCor = (1.+PF)*PF3   #extinction corr for other derivatives
        if 'Primary' in calcControls[phfx+'EType'] and phfx+'Ep' in varyList:
            dervDict[phfx+'Ep'] = -ref[7+im]*PLZ*PF3
        if 'II' in calcControls[phfx+'EType'] and phfx+'Es' in varyList:
            dervDict[phfx+'Es'] = -ref[7+im]*PLZ*PF3*(PSIG/parmDict[phfx+'Es'])**3
        if 'I' in calcControls[phfx+'EType'] and phfx+'Eg' in varyList:
            dervDict[phfx+'Eg'] = -ref[7+im]*PLZ*PF3*(PSIG/parmDict[phfx+'Eg'])**3*PL**2
               
    return 1./extCor,dervDict,dervCor
    
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def GetNewCellParms(parmDict,varyList):
    '''Compute unit cell tensor terms from varied Aij and Dij values.
    Terms are included in the dict only if Aij or Dij is varied.
    '''
    newCellDict = {}
    Anames = ['A'+str(i) for i in range(6)]
    Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],Anames))
    for item in varyList:
        keys = item.split(':')
        if keys[2] in Ddict:
            key = keys[0]+'::'+Ddict[keys[2]]       #key is e.g. '0::A0'
            parm = keys[0]+'::'+keys[2]             #parm is e.g. '0::D11'
            newCellDict[parm] = [key,parmDict[key]+parmDict[item]]
    return newCellDict          # is e.g. {'0::D11':A0-D11}
    
def ApplyXYZshifts(parmDict,varyList):
    '''
    takes atom x,y,z shift and applies it to corresponding atom x,y,z value
    
    :param dict parmDict: parameter dictionary
    :param list varyList: list of variables (not used!)
    :returns: newAtomDict - dictionary of new atomic coordinate names & values; key is parameter shift name

    '''
    newAtomDict = {}
    for item in parmDict:
        if 'dA' in item:
            parm = ''.join(item.split('d'))
            parmDict[parm] += parmDict[item]
            newAtomDict[item] = [parm,parmDict[parm]]
    return newAtomDict
    
def SHTXcal(refl,im,g,pfx,hfx,SGData,calcControls,parmDict):
    'Spherical harmonics texture'
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5+im]
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Phi'],parmDict[hfx+'Chi'],parmDict[hfx+'Omega'],parmDict[hfx+'Azimuth']]
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(tth/2.,Gangls,Sangls,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],parmDict[pfx+'SHmodel'],parmDict[pfx+'SHorder'])
    for item in SHnames:
        L,M,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,x,x = G2lat.GetKsl(L,M,parmDict[pfx+'SHmodel'],psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[pfx+item]*Lnorm*Kcl*Ksl
    return odfCor
    
def SHTXcalDerv(refl,im,g,pfx,hfx,SGData,calcControls,parmDict):
    'Spherical harmonics texture derivatives'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5+im]
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    odfCor = 1.0
    dFdODF = {}
    dFdSA = [0,0,0]
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Phi'],parmDict[hfx+'Chi'],parmDict[hfx+'Omega'],parmDict[hfx+'Azimuth']]
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,dPSdA,dGMdA = G2lat.SamAng(tth/2.,Gangls,Sangls,IFCoup)
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],parmDict[pfx+'SHmodel'],parmDict[pfx+'SHorder'])
    for item in SHnames:
        L,M,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,dKsdp,dKsdg = G2lat.GetKsl(L,M,parmDict[pfx+'SHmodel'],psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[pfx+item]*Lnorm*Kcl*Ksl
        dFdODF[pfx+item] = Lnorm*Kcl*Ksl
        for i in range(3):
            dFdSA[i] += parmDict[pfx+item]*Lnorm*Kcl*(dKsdp*dPSdA[i]+dKsdg*dGMdA[i])
    return odfCor,dFdODF,dFdSA
    
def SHPOcal(refl,im,g,phfx,hfx,SGData,calcControls,parmDict):
    'spherical harmonics preferred orientation (cylindrical symmetry only)'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5+im]
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [parmDict[hfx+'Phi'],parmDict[hfx+'Chi'],parmDict[hfx+'Omega'],parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(tth/2.,Gangls,Sangls,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = calcControls[phfx+'SHnames']
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,x,x = G2lat.GetKsl(L,0,'0',psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[phfx+item]*Lnorm*Kcl*Ksl
    return np.squeeze(odfCor)
    
def SHPOcalDerv(refl,im,g,phfx,hfx,SGData,calcControls,parmDict):
    'spherical harmonics preferred orientation derivatives (cylindrical symmetry only)'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5+im]
    odfCor = 1.0
    dFdODF = {}
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [parmDict[hfx+'Phi'],parmDict[hfx+'Chi'],parmDict[hfx+'Omega'],parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(tth/2.,Gangls,Sangls,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = calcControls[phfx+'SHnames']
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,x,x = G2lat.GetKsl(L,0,'0',psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[phfx+item]*Lnorm*Kcl*Ksl
        dFdODF[phfx+item] = Kcl*Ksl*Lnorm
    return odfCor,dFdODF
    
def GetPrefOri(uniq,G,g,phfx,hfx,SGData,calcControls,parmDict):
    'March-Dollase preferred orientation correction'
    POcorr = 1.0
    MD = parmDict[phfx+'MD']
    if MD != 1.0:
        MDAxis = calcControls[phfx+'MDAxis']
        sumMD = 0
        for H in uniq:            
            cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
            A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
            sumMD += A**3
        POcorr = sumMD/len(uniq)
    return POcorr
    
def GetPrefOriDerv(refl,im,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict):
    'Needs a doc string'
    POcorr = 1.0
    POderv = {}
    if calcControls[phfx+'poType'] == 'MD':
        MD = parmDict[phfx+'MD']
        MDAxis = calcControls[phfx+'MDAxis']
        sumMD = 0
        sumdMD = 0
        for H in uniq:            
            cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
            A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
            sumMD += A**3
            sumdMD -= (1.5*A**5)*(2.0*MD*cosP**2-(sinP/MD)**2)
        POcorr = sumMD/len(uniq)
        POderv[phfx+'MD'] = sumdMD/len(uniq)
    else:   #spherical harmonics
        if calcControls[phfx+'SHord']:
            POcorr,POderv = SHPOcalDerv(refl,im,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr,POderv
    
def GetAbsorb(refl,im,hfx,calcControls,parmDict):
    'Needs a doc string'
    if 'Debye' in calcControls[hfx+'instType']:
        if 'T' in calcControls[hfx+'histType']:
            return G2pwd.Absorb('Cylinder',parmDict[hfx+'Absorption']*refl[14+im],abs(parmDict[hfx+'2-theta']),0,0)
        else:
            return G2pwd.Absorb('Cylinder',parmDict[hfx+'Absorption'],refl[5+im],0,0)
    else:
        return G2pwd.SurfaceRough(parmDict[hfx+'SurfRoughA'],parmDict[hfx+'SurfRoughB'],refl[5+im])
    
def GetAbsorbDerv(refl,im,hfx,calcControls,parmDict):
    'Needs a doc string'
    if 'Debye' in calcControls[hfx+'instType']:
        if 'T' in calcControls[hfx+'histType']:
            return G2pwd.AbsorbDerv('Cylinder',parmDict[hfx+'Absorption']*refl[14+im],abs(parmDict[hfx+'2-theta']),0,0)
        else:
            return G2pwd.AbsorbDerv('Cylinder',parmDict[hfx+'Absorption'],refl[5+im],0,0)
    else:
        return np.array(G2pwd.SurfaceRoughDerv(parmDict[hfx+'SurfRoughA'],parmDict[hfx+'SurfRoughB'],refl[5+im]))
        
def GetPwdrExt(refl,im,pfx,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    coef = np.array([-0.5,0.25,-0.10416667,0.036458333,-0.0109375,2.8497409E-3])
    pi2 = np.sqrt(2./np.pi)
    if 'T' in calcControls[hfx+'histType']:
        sth2 = sind(abs(parmDict[hfx+'2-theta'])/2.)**2
        wave = refl[14+im]
    else:   #'C'W
        sth2 = sind(refl[5+im]/2.)**2
        wave = parmDict.get(hfx+'Lam',parmDict.get(hfx+'Lam1',1.0))
    c2th = 1.-2.0*sth2
    flv2 = refl[9+im]*(wave/parmDict[pfx+'Vol'])**2
    if 'X' in calcControls[hfx+'histType']:
        flv2 *= 0.079411*(1.0+c2th**2)/2.0
    xfac = flv2*parmDict[phfx+'Extinction']
    exb = 1.0
    if xfac > -1.:
        exb = 1./np.sqrt(1.+xfac)
    exl = 1.0
    if 0 < xfac <= 1.:
        xn = np.array([xfac**(i+1) for i in range(6)])
        exl += np.sum(xn*coef)
    elif xfac > 1.:
        xfac2 = 1./np.sqrt(xfac)
        exl = pi2*(1.-0.125/xfac)*xfac2
    return exb*sth2+exl*(1.-sth2)
    
def GetPwdrExtDerv(refl,im,pfx,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    coef = np.array([-0.5,0.25,-0.10416667,0.036458333,-0.0109375,2.8497409E-3])
    pi2 = np.sqrt(2./np.pi)
    if 'T' in calcControls[hfx+'histType']:
        sth2 = sind(abs(parmDict[hfx+'2-theta'])/2.)**2
        wave = refl[14+im]
    else:   #'C'W
        sth2 = sind(refl[5+im]/2.)**2
        wave = parmDict.get(hfx+'Lam',parmDict.get(hfx+'Lam1',1.0))
    c2th = 1.-2.0*sth2
    flv2 = refl[9+im]*(wave/parmDict[pfx+'Vol'])**2
    if 'X' in calcControls[hfx+'histType']:
        flv2 *= 0.079411*(1.0+c2th**2)/2.0
    xfac = flv2*parmDict[phfx+'Extinction']
    dbde = -500.*flv2
    if xfac > -1.:
        dbde = -0.5*flv2/np.sqrt(1.+xfac)**3
    dlde = 0.
    if 0 < xfac <= 1.:
        xn = np.array([i*flv2*xfac**i for i in [1,2,3,4,5,6]])
        dlde = np.sum(xn*coef)/xfac
    elif xfac > 1.:
        xfac2 = 1./np.sqrt(xfac)
        dlde = 0.5*flv2*pi2*xfac2*(-1./xfac+0.375/xfac**2)
        
    return dbde*sth2+dlde*(1.-sth2)
    
def GetIntensityCorr(refl,im,uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    'Needs a doc string'    #need powder extinction!
    Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3+im]               #scale*multiplicity
    if 'X' in parmDict[hfx+'Type']:
        Icorr *= G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5+im],parmDict[hfx+'Azimuth'])[0]
    POcorr = 1.0
    if pfx+'SHorder' in parmDict:                 #generalized spherical harmonics texture - takes precidence
        POcorr = SHTXcal(refl,im,g,pfx,hfx,SGData,calcControls,parmDict)
    elif calcControls[phfx+'poType'] == 'MD':         #March-Dollase
        POcorr = GetPrefOri(uniq,G,g,phfx,hfx,SGData,calcControls,parmDict)
    elif calcControls[phfx+'SHord']:                #cylindrical spherical harmonics
        POcorr = SHPOcal(refl,im,g,phfx,hfx,SGData,calcControls,parmDict)
    Icorr *= POcorr
    AbsCorr = 1.0
    AbsCorr = GetAbsorb(refl,im,hfx,calcControls,parmDict)
    Icorr *= AbsCorr
    ExtCorr = GetPwdrExt(refl,im,pfx,phfx,hfx,calcControls,parmDict)
    Icorr *= ExtCorr
    return Icorr,POcorr,AbsCorr,ExtCorr
    
def GetIntensityDerv(refl,im,wave,uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    'Needs a doc string'    #need powder extinction derivs!
    dIdsh = 1./parmDict[hfx+'Scale']
    dIdsp = 1./parmDict[phfx+'Scale']
    if 'X' in parmDict[hfx+'Type']:
        pola,dIdPola = G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5+im],parmDict[hfx+'Azimuth'])
        dIdPola /= pola
    else:       #'N'
        dIdPola = 0.0
    dFdODF = {}
    dFdSA = [0,0,0]
    dIdPO = {}
    if pfx+'SHorder' in parmDict:
        odfCor,dFdODF,dFdSA = SHTXcalDerv(refl,im,g,pfx,hfx,SGData,calcControls,parmDict)
        for iSH in dFdODF:
            dFdODF[iSH] /= odfCor
        for i in range(3):
            dFdSA[i] /= odfCor
    elif calcControls[phfx+'poType'] == 'MD' or calcControls[phfx+'SHord']:
        POcorr,dIdPO = GetPrefOriDerv(refl,im,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict)        
        for iPO in dIdPO:
            dIdPO[iPO] /= POcorr
    if 'T' in parmDict[hfx+'Type']:
        dFdAb = GetAbsorbDerv(refl,im,hfx,calcControls,parmDict)*wave/refl[16+im] #wave/abs corr
        dFdEx = GetPwdrExtDerv(refl,im,pfx,phfx,hfx,calcControls,parmDict)/refl[17+im]    #/ext corr
    else:
        dFdAb = GetAbsorbDerv(refl,im,hfx,calcControls,parmDict)*wave/refl[13+im] #wave/abs corr
        dFdEx = GetPwdrExtDerv(refl,im,pfx,phfx,hfx,calcControls,parmDict)/refl[14+im]    #/ext corr        
    return dIdsh,dIdsp,dIdPola,dIdPO,dFdODF,dFdSA,dFdAb,dFdEx
        
def GetSampleSigGam(refl,im,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:     #All checked & OK
        costh = cosd(refl[5+im]/2.)
        #crystallite size
        if calcControls[phfx+'SizeType'] == 'isotropic':
            Sgam = 1.8*wave/(np.pi*parmDict[phfx+'Size;i']*costh)
        elif calcControls[phfx+'SizeType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Sgam = (1.8*wave/np.pi)/(parmDict[phfx+'Size;i']*parmDict[phfx+'Size;a']*costh)
            Sgam *= np.sqrt((sinP*parmDict[phfx+'Size;a'])**2+(cosP*parmDict[phfx+'Size;i'])**2)
        else:           #ellipsoidal crystallites
            Sij =[parmDict[phfx+'Size;%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR = G2pwd.ellipseSize(H,Sij,GB)
            Sgam = 1.8*wave/(np.pi*costh*lenR)
        #microstrain                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 0.018*parmDict[phfx+'Mustrain;i']*tand(refl[5+im]/2.)/np.pi
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            Mgam = 0.018*Si*Sa*tand(refl[5+im]/2.)/(np.pi*np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
        else:       #generalized - P.W. Stephens model
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain;'+str(i)]*strm
            Mgam = 0.018*refl[4+im]**2*tand(refl[5+im]/2.)*np.sqrt(Sum)/np.pi
    elif 'T' in calcControls[hfx+'histType']:       #All checked & OK
        #crystallite size
        if calcControls[phfx+'SizeType'] == 'isotropic':    #OK
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2/parmDict[phfx+'Size;i']
        elif calcControls[phfx+'SizeType'] == 'uniaxial':   #OK
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2/(parmDict[phfx+'Size;i']*parmDict[phfx+'Size;a'])
            Sgam *= np.sqrt((sinP*parmDict[phfx+'Size;a'])**2+(cosP*parmDict[phfx+'Size;i'])**2)
        else:           #ellipsoidal crystallites   #OK
            Sij =[parmDict[phfx+'Size;%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR = G2pwd.ellipseSize(H,Sij,GB)
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2/lenR
        #microstrain                
        if calcControls[phfx+'MustrainType'] == 'isotropic':    #OK
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4+im]*parmDict[phfx+'Mustrain;i']
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':   #OK
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4+im]*Si*Sa/np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
        else:       #generalized - P.W. Stephens model  OK
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain;'+str(i)]*strm
            Mgam = 1.e-6*parmDict[hfx+'difC']*np.sqrt(Sum)*refl[4+im]**3
            
    gam = Sgam*parmDict[phfx+'Size;mx']+Mgam*parmDict[phfx+'Mustrain;mx']
    sig = (Sgam*(1.-parmDict[phfx+'Size;mx']))**2+(Mgam*(1.-parmDict[phfx+'Mustrain;mx']))**2
    sig /= ateln2
    return sig,gam
        
def GetSampleSigGamDerv(refl,im,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    gamDict = {}
    sigDict = {}
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:         #All checked & OK
        costh = cosd(refl[5+im]/2.)
        tanth = tand(refl[5+im]/2.)
        #crystallite size derivatives
        if calcControls[phfx+'SizeType'] == 'isotropic':
            Sgam = 1.8*wave/(np.pi*costh*parmDict[phfx+'Size;i'])
            gamDict[phfx+'Size;i'] = -1.8*wave*parmDict[phfx+'Size;mx']/(np.pi*costh*parmDict[phfx+'Size;i']**2)
            sigDict[phfx+'Size;i'] = -3.6*Sgam*wave*(1.-parmDict[phfx+'Size;mx'])**2/(np.pi*costh*ateln2)
        elif calcControls[phfx+'SizeType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Size;i']
            Sa = parmDict[phfx+'Size;a']
            gami = 1.8*wave/(costh*np.pi*Si*Sa)
            sqtrm = np.sqrt((sinP*Sa)**2+(cosP*Si)**2)
            Sgam = gami*sqtrm
            dsi = gami*Si*cosP**2/sqtrm-Sgam/Si
            dsa = gami*Sa*sinP**2/sqtrm-Sgam/Sa
            gamDict[phfx+'Size;i'] = dsi*parmDict[phfx+'Size;mx']
            gamDict[phfx+'Size;a'] = dsa*parmDict[phfx+'Size;mx']
            sigDict[phfx+'Size;i'] = 2.*dsi*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
            sigDict[phfx+'Size;a'] = 2.*dsa*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        else:           #ellipsoidal crystallites
            const = 1.8*wave/(np.pi*costh)
            Sij =[parmDict[phfx+'Size;%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
            Sgam = const/lenR
            for i,item in enumerate([phfx+'Size;%d'%(j) for j in range(6)]):
                gamDict[item] = -(const/lenR**2)*dRdS[i]*parmDict[phfx+'Size;mx']
                sigDict[item] = -2.*Sgam*(const/lenR**2)*dRdS[i]*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        gamDict[phfx+'Size;mx'] = Sgam
        sigDict[phfx+'Size;mx'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])/ateln2
                
        #microstrain derivatives                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 0.018*parmDict[phfx+'Mustrain;i']*tand(refl[5+im]/2.)/np.pi
            gamDict[phfx+'Mustrain;i'] =  0.018*tanth*parmDict[phfx+'Mustrain;mx']/np.pi
            sigDict[phfx+'Mustrain;i'] =  0.036*Mgam*tanth*(1.-parmDict[phfx+'Mustrain;mx'])**2/(np.pi*ateln2)        
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            gami = 0.018*Si*Sa*tanth/np.pi
            sqtrm = np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
            Mgam = gami/sqtrm
            dsi = -gami*Si*cosP**2/sqtrm**3
            dsa = -gami*Sa*sinP**2/sqtrm**3
            gamDict[phfx+'Mustrain;i'] = (Mgam/Si+dsi)*parmDict[phfx+'Mustrain;mx']
            gamDict[phfx+'Mustrain;a'] = (Mgam/Sa+dsa)*parmDict[phfx+'Mustrain;mx']
            sigDict[phfx+'Mustrain;i'] = 2*(Mgam/Si+dsi)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2
            sigDict[phfx+'Mustrain;a'] = 2*(Mgam/Sa+dsa)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2       
        else:       #generalized - P.W. Stephens model
            const = 0.018*refl[4+im]**2*tanth/np.pi
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain;'+str(i)]*strm
                gamDict[phfx+'Mustrain;'+str(i)] = strm*parmDict[phfx+'Mustrain;mx']/2.
                sigDict[phfx+'Mustrain;'+str(i)] = strm*(1.-parmDict[phfx+'Mustrain;mx'])**2
            Mgam = const*np.sqrt(Sum)
            for i in range(len(Strms)):
                gamDict[phfx+'Mustrain;'+str(i)] *= Mgam/Sum
                sigDict[phfx+'Mustrain;'+str(i)] *= const**2/ateln2
        gamDict[phfx+'Mustrain;mx'] = Mgam
        sigDict[phfx+'Mustrain;mx'] = -2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])/ateln2
    else:   #'T'OF - All checked & OK
        if calcControls[phfx+'SizeType'] == 'isotropic':    #OK
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2/parmDict[phfx+'Size;i']
            gamDict[phfx+'Size;i'] = -Sgam*parmDict[phfx+'Size;mx']/parmDict[phfx+'Size;i']
            sigDict[phfx+'Size;i'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])**2/(ateln2*parmDict[phfx+'Size;i'])
        elif calcControls[phfx+'SizeType'] == 'uniaxial':   #OK
            const = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Size;i']
            Sa = parmDict[phfx+'Size;a']
            gami = const/(Si*Sa)
            sqtrm = np.sqrt((sinP*Sa)**2+(cosP*Si)**2)
            Sgam = gami*sqtrm
            dsi = gami*Si*cosP**2/sqtrm-Sgam/Si
            dsa = gami*Sa*sinP**2/sqtrm-Sgam/Sa
            gamDict[phfx+'Size;i'] = dsi*parmDict[phfx+'Size;mx']
            gamDict[phfx+'Size;a'] = dsa*parmDict[phfx+'Size;mx']
            sigDict[phfx+'Size;i'] = 2.*dsi*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
            sigDict[phfx+'Size;a'] = 2.*dsa*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        else:           #OK  ellipsoidal crystallites 
            const = 1.e-4*parmDict[hfx+'difC']*refl[4+im]**2
            Sij =[parmDict[phfx+'Size;%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
            Sgam = const/lenR
            for i,item in enumerate([phfx+'Size;%d'%(j) for j in range(6)]):
                gamDict[item] = -(const/lenR**2)*dRdS[i]*parmDict[phfx+'Size;mx']
                sigDict[item] = -2.*Sgam*(const/lenR**2)*dRdS[i]*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        gamDict[phfx+'Size;mx'] = Sgam  #OK
        sigDict[phfx+'Size;mx'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])/ateln2  #OK
                
        #microstrain derivatives                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4+im]*parmDict[phfx+'Mustrain;i']
            gamDict[phfx+'Mustrain;i'] =  1.e-6*refl[4+im]*parmDict[hfx+'difC']*parmDict[phfx+'Mustrain;mx']   #OK
            sigDict[phfx+'Mustrain;i'] =  2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])**2/(ateln2*parmDict[phfx+'Mustrain;i'])        
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            gami = 1.e-6*parmDict[hfx+'difC']*refl[4+im]*Si*Sa
            sqtrm = np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
            Mgam = gami/sqtrm
            dsi = -gami*Si*cosP**2/sqtrm**3
            dsa = -gami*Sa*sinP**2/sqtrm**3
            gamDict[phfx+'Mustrain;i'] = (Mgam/Si+dsi)*parmDict[phfx+'Mustrain;mx']
            gamDict[phfx+'Mustrain;a'] = (Mgam/Sa+dsa)*parmDict[phfx+'Mustrain;mx']
            sigDict[phfx+'Mustrain;i'] = 2*(Mgam/Si+dsi)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2
            sigDict[phfx+'Mustrain;a'] = 2*(Mgam/Sa+dsa)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2       
        else:       #generalized - P.W. Stephens model OK
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            const = 1.e-6*parmDict[hfx+'difC']*refl[4+im]**3
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain;'+str(i)]*strm
                gamDict[phfx+'Mustrain;'+str(i)] = strm*parmDict[phfx+'Mustrain;mx']/2.
                sigDict[phfx+'Mustrain;'+str(i)] = strm*(1.-parmDict[phfx+'Mustrain;mx'])**2
            Mgam = const*np.sqrt(Sum)
            for i in range(len(Strms)):
                gamDict[phfx+'Mustrain;'+str(i)] *= Mgam/Sum
                sigDict[phfx+'Mustrain;'+str(i)] *= const**2/ateln2        
        gamDict[phfx+'Mustrain;mx'] = Mgam
        sigDict[phfx+'Mustrain;mx'] = -2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])/ateln2
        
    return sigDict,gamDict
        
def GetReflPos(refl,im,wave,A,pfx,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    if im:
        h,k,l,m = refl[:4]
        vec = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])
        d = 1./np.sqrt(G2lat.calc_rDsqSS(np.array([h,k,l,m]),A,vec))
    else:
        h,k,l = refl[:3]
        d = 1./np.sqrt(G2lat.calc_rDsq(np.array([h,k,l]),A))
    refl[4+im] = d
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:
        pos = 2.0*asind(wave/(2.0*d))+parmDict[hfx+'Zero']
        const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            pos -= const*(4.*parmDict[hfx+'Shift']*cosd(pos/2.0)+ \
                parmDict[hfx+'Transparency']*sind(pos)*100.0)            #trans(=1/mueff) in cm
        else:               #Debye-Scherrer - simple but maybe not right
            pos -= const*(parmDict[hfx+'DisplaceX']*cosd(pos)+(parmDict[hfx+'DisplaceY']+parmDict[phfx+'LayerDisp'])*sind(pos))
    elif 'T' in calcControls[hfx+'histType']:
        pos = parmDict[hfx+'difC']*d+parmDict[hfx+'difA']*d**2+parmDict[hfx+'difB']/d+parmDict[hfx+'Zero']
        #do I need sample position effects - maybe?
    return pos

def GetReflPosDerv(refl,im,wave,A,pfx,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    dpr = 180./np.pi
    if im:
        h,k,l,m = refl[:4]
        vec = np.array([parmDict[pfx+'mV0'],parmDict[pfx+'mV1'],parmDict[pfx+'mV2']])
        dstsq = G2lat.calc_rDsqSS(np.array([h,k,l,m]),A,vec)
        h,k,l = [h+m*vec[0],k+m*vec[1],l+m*vec[2]]          #do proj of hklm to hkl so dPdA & dPdV come out right
    else:
        m = 0
        h,k,l = refl[:3]        
        dstsq = G2lat.calc_rDsq(np.array([h,k,l]),A)
    dst = np.sqrt(dstsq)
    dsp = 1./dst
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:
        pos = refl[5+im]-parmDict[hfx+'Zero']
        const = dpr/np.sqrt(1.0-wave**2*dstsq/4.0)
        dpdw = const*dst
        dpdA = np.array([h**2,k**2,l**2,h*k,h*l,k*l])*const*wave/(2.0*dst)
        dpdZ = 1.0
        dpdV = np.array([2.*h*A[0]+k*A[3]+l*A[4],2*k*A[1]+h*A[3]+l*A[5],
            2*l*A[2]+h*A[4]+k*A[5]])*m*const*wave/(2.0*dst)
        shft = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            dpdSh = -4.*shft*cosd(pos/2.0)
            dpdTr = -shft*sind(pos)*100.0
            return dpdA,dpdw,dpdZ,dpdSh,dpdTr,0.,0.,dpdV
        else:               #Debye-Scherrer - simple but maybe not right
            dpdXd = -shft*cosd(pos)
            dpdYd = -shft*sind(pos)
            return dpdA,dpdw,dpdZ,0.,0.,dpdXd,dpdYd,dpdV
    elif 'T' in calcControls[hfx+'histType']:
        dpdA = -np.array([h**2,k**2,l**2,h*k,h*l,k*l])*parmDict[hfx+'difC']*dsp**3/2.
        dpdZ = 1.0
        dpdDC = dsp
        dpdDA = dsp**2
        dpdDB = 1./dsp
        dpdV = np.array([2.*h*A[0]+k*A[3]+l*A[4],2*k*A[1]+h*A[3]+l*A[5],
            2*l*A[2]+h*A[4]+k*A[5]])*m*parmDict[hfx+'difC']*dsp**3/2.
        return dpdA,dpdZ,dpdDC,dpdDA,dpdDB,dpdV
            
def GetHStrainShift(refl,im,SGData,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+l**2)+ \
            refl[4+im]**2*parmDict[phfx+'eA']*((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+h*k)+parmDict[phfx+'D33']*l**2
    elif laue in ['3R','3mR']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+l**2)+parmDict[phfx+'D12']*(h*k+h*l+k*l)
    elif laue in ['4/m','4/mmm']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2)+parmDict[phfx+'D33']*l**2
    elif laue in ['mmm']:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2
    elif laue in ['2/m']:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2
        if uniq == 'a':
            Dij += parmDict[phfx+'D23']*k*l
        elif uniq == 'b':
            Dij += parmDict[phfx+'D13']*h*l
        elif uniq == 'c':
            Dij += parmDict[phfx+'D12']*h*k
    else:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2+ \
            parmDict[phfx+'D12']*h*k+parmDict[phfx+'D13']*h*l+parmDict[phfx+'D23']*k*l
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:
        return -180.*Dij*refl[4+im]**2*tand(refl[5+im]/2.0)/np.pi
    else:
        return -Dij*parmDict[hfx+'difC']*refl[4+im]**2/2.
            
def GetHStrainShiftDerv(refl,im,SGData,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        dDijDict = {phfx+'D11':h**2+k**2+l**2,
            phfx+'eA':refl[4+im]**2*((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2}
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        dDijDict = {phfx+'D11':h**2+k**2+h*k,phfx+'D33':l**2}
    elif laue in ['3R','3mR']:
        dDijDict = {phfx+'D11':h**2+k**2+l**2,phfx+'D12':h*k+h*l+k*l}
    elif laue in ['4/m','4/mmm']:
        dDijDict = {phfx+'D11':h**2+k**2,phfx+'D33':l**2}
    elif laue in ['mmm']:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2}
    elif laue in ['2/m']:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2}
        if uniq == 'a':
            dDijDict[phfx+'D23'] = k*l
        elif uniq == 'b':
            dDijDict[phfx+'D13'] = h*l
        elif uniq == 'c':
            dDijDict[phfx+'D12'] = h*k
    else:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2,
            phfx+'D12':h*k,phfx+'D13':h*l,phfx+'D23':k*l}
    if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:
        for item in dDijDict:
            dDijDict[item] *= 180.0*refl[4+im]**2*tand(refl[5+im]/2.0)/np.pi
    else:
        for item in dDijDict:
            dDijDict[item] *= -parmDict[hfx+'difC']*refl[4+im]**3/2.
    return dDijDict
    
def GetDij(phfx,SGData,parmDict):
    HSvals = [parmDict[phfx+name] for name in G2spc.HStrainNames(SGData)]
    return G2spc.HStrainVals(HSvals,SGData)
                
def GetFobsSq(Histograms,Phases,parmDict,calcControls):
    '''Compute the observed structure factors for Powder histograms and store in reflection array
    Multiprocessing support added
    '''
    if GSASIIpath.GetConfigValue('Show_timing',False):
        starttime = time.time() #; print 'start GetFobsSq'
    histoList = list(Histograms.keys())
    histoList.sort()
    Ka2 = shl = lamRatio = kRatio = None
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            if 'C' in calcControls[hfx+'histType']:
                shl = max(parmDict[hfx+'SH/L'],0.0005)
                Ka2 = False
                kRatio = 0.0
                if hfx+'Lam1' in list(parmDict.keys()):
                    Ka2 = True
                    lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
                    kRatio = parmDict[hfx+'I(L2)/I(L1)']
            x,y,w,yc,yb,yd = Histogram['Data']
            xMask = ma.getmaskarray(x)
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            ymb = np.array(y-yb)
            ymb = np.where(ymb,ymb,1.0)
            ycmb = np.array(yc-yb)
            ratio = 1./np.where(ycmb,ycmb/ymb,1.e10)          
            refLists = Histogram['Reflection Lists']
            for phase in refLists:
                if phase not in Phases:     #skips deleted or renamed phases silently!
                    continue
                Phase = Phases[phase]
                if histogram not in Phase['Histograms']:
                    continue
                im = 0
                if Phase['General'].get('Modulated',False):
                    im = 1
                pId = Phase['pId']
                phfx = '%d:%d:'%(pId,hId)
                refDict = refLists[phase]
                sumFo = 0.0
                sumdF = 0.0
                sumFosq = 0.0
                sumdFsq = 0.0
                sumInt = 0.0
                nExcl = 0
                # test to see if we are using multiprocessing below
                useMP,ncores = G2mp.InitMP()
                if len(refDict['RefList']) < 100: useMP = False        
                if useMP: # multiprocessing: create a set of initialized Python processes
                    MPpool = mp.Pool(G2mp.ncores,G2mp.InitFobsSqGlobals,
                                    [x,ratio,shl,xB,xF,im,lamRatio,kRatio,xMask,Ka2])
                    profArgs = [[] for i in range(G2mp.ncores)]
                else:
                    G2mp.InitFobsSqGlobals(x,ratio,shl,xB,xF,im,lamRatio,kRatio,xMask,Ka2)
                if 'C' in calcControls[hfx+'histType']:
                    # are we multiprocessing?
                    for iref,refl in enumerate(refDict['RefList']):
                        if useMP: 
                            profArgs[iref%G2mp.ncores].append((refl,iref))
                        else:
                            icod= G2mp.ComputeFobsSqCW(refl,iref)
                            if type(icod) is tuple:
                                refl[8+im] = icod[0]
                                sumInt += icod[1]
                                if parmDict[phfx+'LeBail']: refl[9+im] = refl[8+im]
                            elif icod == -1:
                                refl[3+im] *= -1
                                nExcl += 1
                            elif icod == -2:
                                break
                    if useMP:
                        for sInt,resList in MPpool.imap_unordered(G2mp.ComputeFobsSqCWbatch,profArgs):
                            sumInt += sInt
                            for refl8im,irefl in resList:
                                if refl8im is None:
                                    refDict['RefList'][irefl][3+im] *= -1
                                    nExcl += 1
                                else:
                                    refDict['RefList'][irefl][8+im] = refl8im
                                    if parmDict[phfx+'LeBail']:
                                        refDict['RefList'][irefl][9+im] = refDict['RefList'][irefl][8+im]
                elif 'T' in calcControls[hfx+'histType']:
                    for iref,refl in enumerate(refDict['RefList']):
                        if useMP: 
                            profArgs[iref%G2mp.ncores].append((refl,iref))
                        else:
                            icod= G2mp.ComputeFobsSqTOF(refl,iref)
                            if type(icod) is tuple:
                                refl[8+im] = icod[0]
                                sumInt += icod[1]
                                if parmDict[phfx+'LeBail']: refl[9+im] = refl[8+im]
                            elif icod == -1:
                                refl[3+im] *= -1
                                nExcl += 1
                            elif icod == -2:
                                break
                    if useMP:
                        for sInt,resList in MPpool.imap_unordered(G2mp.ComputeFobsSqTOFbatch,profArgs):
                            sumInt += sInt
                            for refl8im,irefl in resList:
                                if refl8im is None:
                                    refDict['RefList'][irefl][3+im] *= -1
                                    nExcl += 1
                                else:
                                    refDict['RefList'][irefl][8+im] = refl8im
                                    if parmDict[phfx+'LeBail']:
                                        refDict['RefList'][irefl][9+im] = refDict['RefList'][irefl][8+im]
                elif 'B' in calcControls[hfx+'histType']:
                    for iref,refl in enumerate(refDict['RefList']):
                        if useMP: 
                            profArgs[iref%G2mp.ncores].append((refl,iref))
                        else:
                            icod= G2mp.ComputeFobsSqPink(refl,iref)
                            if type(icod) is tuple:
                                refl[8+im] = icod[0]
                                sumInt += icod[1]
                                if parmDict[phfx+'LeBail']: refl[9+im] = refl[8+im]
                            elif icod == -1:
                                refl[3+im] *= -1
                                nExcl += 1
                            elif icod == -2:
                                break
                    if useMP:
                        for sInt,resList in MPpool.imap_unordered(G2mp.ComputeFobsSqPinkbatch,profArgs):
                            sumInt += sInt
                            for refl8im,irefl in resList:
                                if refl8im is None:
                                    refDict['RefList'][irefl][3+im] *= -1
                                    nExcl += 1
                                else:
                                    refDict['RefList'][irefl][8+im] = refl8im
                                    if parmDict[phfx+'LeBail']:
                                        refDict['RefList'][irefl][9+im] = refDict['RefList'][irefl][8+im]
                if useMP: MPpool.terminate()
                sumFo = 0.0
                sumdF = 0.0
                sumFosq = 0.0
                sumdFsq = 0.0
                for iref,refl in enumerate(refDict['RefList']):
                    Fo = np.sqrt(np.abs(refl[8+im]))
                    Fc = np.sqrt(np.abs(refl[9]+im))
                    sumFo += Fo
                    sumFosq += refl[8+im]**2
                    sumdF += np.abs(Fo-Fc)
                    sumdFsq += (refl[8+im]-refl[9+im])**2
                if sumFo:
                    Histogram['Residuals'][phfx+'Rf'] = min(100.,(sumdF/sumFo)*100.)
                    Histogram['Residuals'][phfx+'Rf^2'] = min(100.,np.sqrt(sumdFsq/sumFosq)*100.)
                else:
                    Histogram['Residuals'][phfx+'Rf'] = 100.
                    Histogram['Residuals'][phfx+'Rf^2'] = 100.
                Histogram['Residuals'][phfx+'sumInt'] = sumInt
                Histogram['Residuals'][phfx+'Nref'] = len(refDict['RefList'])-nExcl
                Histogram['Residuals']['hId'] = hId
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            Histogram['Residuals']['hId'] = Histograms[histogram]['hId']
    if GSASIIpath.GetConfigValue('Show_timing',False):
        print ('GetFobsSq t=',time.time()-starttime)
                
def getPowderProfile(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup,histogram=None):
    'Computes the powder pattern for a histogram based on contributions from all used phases'
    if GSASIIpath.GetConfigValue('Show_timing',False): starttime = time.time()
    
    def GetReflSigGamCW(refl,im,wave,G,GB,phfx,calcControls,parmDict):
        U = parmDict[hfx+'U']
        V = parmDict[hfx+'V']
        W = parmDict[hfx+'W']
        X = parmDict[hfx+'X']
        Y = parmDict[hfx+'Y']
        Z = parmDict[hfx+'Z']
        tanPos = tand(refl[5+im]/2.0)
        Ssig,Sgam = GetSampleSigGam(refl,im,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict)
        sig = U*tanPos**2+V*tanPos+W+Ssig     #save peak sigma
        sig = max(0.001,sig)
        gam = X/cosd(refl[5+im]/2.0)+Y*tanPos+Sgam+Z     #save peak gamma
        gam = max(0.001,gam)
        return sig,gam
                
    def GetReflSigGamTOF(refl,im,G,GB,phfx,calcControls,parmDict):
        sig = parmDict[hfx+'sig-0']+parmDict[hfx+'sig-1']*refl[4+im]**2+   \
            parmDict[hfx+'sig-2']*refl[4+im]**4+parmDict[hfx+'sig-q']*refl[4+im]
        gam = parmDict[hfx+'X']*refl[4+im]+parmDict[hfx+'Y']*refl[4+im]**2+parmDict[hfx+'Z']
        Ssig,Sgam = GetSampleSigGam(refl,im,0.0,G,GB,SGData,hfx,phfx,calcControls,parmDict)
        sig += Ssig
        gam += Sgam
        return sig,gam
        
    def GetReflAlpBet(refl,im,hfx,parmDict):
        alp = parmDict[hfx+'alpha']/refl[4+im]
        bet = parmDict[hfx+'beta-0']+parmDict[hfx+'beta-1']/refl[4+im]**4+parmDict[hfx+'beta-q']/refl[4+im]**2
        return alp,bet
        
    def GetPinkReflAlpBet(refl,im,hfx,parmDict):
        tanPos = tand(refl[5+im]/2.0)
        alp = parmDict[hfx+'alpha-0']+parmDict[hfx+'alpha-0']*tanPos
        bet = parmDict[hfx+'beta-0']+parmDict[hfx+'beta-1']*tanPos
        return alp,bet

    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    fixedBkg = {i:Histogram['Background'][1][i] for i in Histogram['Background'][1] if i.startswith("_")} 
    yb,Histogram['sumBk'] = G2pwd.getBackground(hfx,parmDict,bakType,calcControls[hfx+'histType'],x,fixedBkg)
    yc = np.zeros_like(yb)
    cw = np.diff(ma.getdata(x))
    cw = np.append(cw,cw[-1])
        
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.002)
        Ka2 = False
        if hfx+'Lam1' in (parmDict.keys()):
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    elif 'B' in calcControls[hfx+'histType']:
            wave = parmDict[hfx+'Lam']
    else:
        shl = 0.
    for phase in Histogram['Reflection Lists']:
        refDict = Histogram['Reflection Lists'][phase]
        if phase not in Phases:     #skips deleted or renamed phases silently!
            continue
        Phase = Phases[phase]
        if histogram and not histogram in Phase['Histograms']:
            continue
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        hfx = ':%d:'%(hId)
        SGData = Phase['General']['SGData']
        SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
        im = 0
        if Phase['General'].get('Modulated',False):
            SSGData = Phase['General']['SSGData']
            im = 1  #offset in SS reflection list
            #??
        Dij = GetDij(phfx,SGData,parmDict)
        A = [parmDict[pfx+'A%d'%(i)]+Dij[i] for i in range(6)]  #TODO: need to do someting if Dij << 0. 
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        if np.any(np.diag(G)<0.) or np.any(np.isnan(A)):
            raise G2obj.G2Exception('invalid metric tensor \n cell/Dij refinement not advised')
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        Vst = np.sqrt(nl.det(G))    #V*
        if not Phase['General'].get('doPawley') and not parmDict[phfx+'LeBail']:
            if im:
                SStructureFactor(refDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
            elif parmDict[pfx+'isMag'] and 'N' in calcControls[hfx+'histType']:
                MagStructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
            else:
                StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
        badPeak = False
        # test to see if we are using multiprocessing here
        useMP,ncores = G2mp.InitMP()
        if len(refDict['RefList']) < 100: useMP = False        
        if useMP: # multiprocessing: create a set of initialized Python processes
            MPpool = mp.Pool(ncores,G2mp.InitPwdrProfGlobals,[im,shl,x])
            profArgs = [[] for i in range(ncores)]
        if 'C' in calcControls[hfx+'histType']:
            for iref,refl in enumerate(refDict['RefList']):
                if im:
                    h,k,l,m = refl[:4]
                else:
                    h,k,l = refl[:3]
                Uniq = np.inner(refl[:3],SGMT)
                refl[5+im] = GetReflPos(refl,im,wave,A,pfx,hfx,phfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = 1./(2.*sind(refl[5+im]/2.)**2*cosd(refl[5+im]/2.))           #Lorentz correction
                refl[6+im:8+im] = GetReflSigGamCW(refl,im,wave,G,GB,phfx,calcControls,parmDict)    #peak sig & gam
                refl[11+im:15+im] = GetIntensityCorr(refl,im,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                refl[11+im] *= Vst*Lorenz
                 
                if Phase['General'].get('doPawley'):
                    try:
                        if im:
                            pInd = pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d,%d'%(h,k,l,m)])
                        else:
                            pInd = pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        refl[9+im] = parmDict[pInd]
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5+im],refl[6+im],refl[7+im],shl)
                iBeg = np.searchsorted(x,refl[5+im]-fmin)
                iFin = np.searchsorted(x,refl[5+im]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                elif iBeg > iFin:   #bad peak coeff - skip
                    badPeak = True
                    continue
                if useMP:
                    profArgs[iref%ncores].append((refl[5+im],refl,iBeg,iFin,1.))
                else:
                    yc[iBeg:iFin] += refl[11+im]*refl[9+im]*G2pwd.getFCJVoigt3(refl[5+im],refl[6+im],refl[7+im],shl,ma.getdata(x[iBeg:iFin]))    #>90% of time spent here
                if Ka2:
                    pos2 = refl[5+im]+lamRatio*tand(refl[5+im]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6+im],refl[7+im],shl)
                    iBeg = np.searchsorted(x,pos2-fmin)
                    iFin = np.searchsorted(x,pos2+fmax)
                    if not iBeg+iFin:       #peak below low limit - skip peak
                        continue
                    elif not iBeg-iFin:     #peak above high limit - done
                        return yc,yb
                    elif iBeg > iFin:   #bad peak coeff - skip
                        continue
                    if useMP:
                        profArgs[iref%ncores].append((pos2,refl,iBeg,iFin,kRatio))
                    else:
                        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6+im],refl[7+im],shl,ma.getdata(x[iBeg:iFin]))        #and here
        elif 'B' in calcControls[hfx+'histType']:
            for iref,refl in enumerate(refDict['RefList']):
                if im:
                    h,k,l,m = refl[:4]
                else:
                    h,k,l = refl[:3]
                Uniq = np.inner(refl[:3],SGMT)
                refl[5+im] = GetReflPos(refl,im,wave,A,pfx,hfx,phfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = 1./(2.*sind(refl[5+im]/2.)**2*cosd(refl[5+im]/2.))           #Lorentz correction
                refl[6+im:8+im] = GetReflSigGamCW(refl,im,wave,G,GB,phfx,calcControls,parmDict)    #peak sig & gam
                refl[12+im:14+im] = GetPinkReflAlpBet(refl,im,hfx,parmDict)             #TODO - skip if alp, bet tabulated?
                refl[11+im],refl[14+im],refl[15+im],refl[16+im] = GetIntensityCorr(refl,im,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                refl[11+im] *= Vst*Lorenz
                 
                if Phase['General'].get('doPawley'):
                    try:
                        if im:
                            pInd = pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d,%d'%(h,k,l,m)])
                        else:
                            pInd = pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        refl[9+im] = parmDict[pInd]
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.)
                iBeg = np.searchsorted(x,refl[5+im]-fmin)
                iFin = np.searchsorted(x,refl[5+im]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                elif iBeg > iFin:   #bad peak coeff - skip
                    badPeak = True
                    continue
                if useMP:
                    profArgs[iref%ncores].append((refl[5+im],refl,iBeg,iFin,1.))
                else:
                    yc[iBeg:iFin] += refl[11+im]*refl[9+im]*G2pwd.getEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.,ma.getdata(x[iBeg:iFin]))
        elif 'T' in calcControls[hfx+'histType']:
            for iref,refl in enumerate(refDict['RefList']):
                if im:
                    h,k,l,m = refl[:4]
                else:
                    h,k,l = refl[:3]
                Uniq = np.inner(refl[:3],SGMT)
                refl[5+im] = GetReflPos(refl,im,0.0,A,pfx,hfx,phfx,calcControls,parmDict)         #corrected reflection position - #TODO - what about tabluated offset?
                Lorenz = sind(abs(parmDict[hfx+'2-theta'])/2)*refl[4+im]**4                                                #TOF Lorentz correction
#                refl[5+im] += GetHStrainShift(refl,im,SGData,phfx,hfx,calcControls,parmDict)               #apply hydrostatic strain shift
                refl[6+im:8+im] = GetReflSigGamTOF(refl,im,G,GB,phfx,calcControls,parmDict)    #peak sig & gam
                refl[12+im:14+im] = GetReflAlpBet(refl,im,hfx,parmDict)             #TODO - skip if alp, bet tabulated?
                refl[11+im],refl[15+im],refl[16+im],refl[17+im] = GetIntensityCorr(refl,im,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                refl[11+im] *= Vst*Lorenz
                if Phase['General'].get('doPawley'):
                    try:
                        if im:
                            pInd =pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d,%d'%(h,k,l,m)])
                        else:
                            pInd =pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        refl[9+im] = parmDict[pInd]
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im])
                iBeg = np.searchsorted(x,refl[5+im]-fmin)
                iFin = np.searchsorted(x,refl[5+im]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                elif iBeg > iFin:   #bad peak coeff - skip
                    badPeak = True
                    continue
                if useMP:
                    profArgs[iref%ncores].append((refl[5+im],refl,iBeg,iFin))
                else:
                    yc[iBeg:iFin] += refl[11+im]*refl[9+im]*G2pwd.getEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],ma.getdata(x[iBeg:iFin]))/cw[iBeg:iFin]
#        print 'profile calc time: %.3fs'%(time.time()-time0)
        if useMP and 'C' in calcControls[hfx+'histType']:
            for y in MPpool.imap_unordered(G2mp.ComputePwdrProfCW,profArgs):
                yc += y
            MPpool.terminate()
        elif useMP and 'T' in calcControls[hfx+'histType']:
            for y in MPpool.imap_unordered(G2mp.ComputePwdrProfTOF,profArgs):
                yc += y
            MPpool.terminate()
        elif useMP and 'B' in calcControls[hfx+'histType']:
            for y in MPpool.imap_unordered(G2mp.ComputePwdrProfPink,profArgs):
                yc += y
            MPpool.terminate()
    if badPeak:
        print ('ouch #4 bad profile coefficients yield negative peak width; some reflections skipped')
    if GSASIIpath.GetConfigValue('Show_timing',False):
        print ('getPowderProfile t=%.3f'%(time.time()-starttime))
    return yc,yb
    
def getPowderProfileDervMP(args):
    '''Computes the derivatives of the computed powder pattern with respect to all
    refined parameters.
    Multiprocessing version.
    '''
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics for each processor
    parmDict,x,varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars = args[:9]
    prc,tprc,histogram = 0,1,None
    if len(args) >= 10: prc=args[9]
    if len(args) >= 11: tprc=args[10]
    if len(args) >= 12: histogram=args[11]
    def cellVaryDerv(pfx,SGData,dpdA): 
        if SGData['SGLaue'] in ['-1',]:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],
                [pfx+'A3',dpdA[3]],[pfx+'A4',dpdA[4]],[pfx+'A5',dpdA[5]]]
        elif SGData['SGLaue'] in ['2/m',]:
            if SGData['SGUniq'] == 'a':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A5',dpdA[5]]]
            elif SGData['SGUniq'] == 'b':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A4',dpdA[4]]]
            else:
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A3',dpdA[3]]]
        elif SGData['SGLaue'] in ['mmm',]:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['4/m','4/mmm']:
            return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
            return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['3R', '3mR']:
            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]],[pfx+'A3',dpdA[3]+dpdA[4]+dpdA[5]]]                       
        elif SGData['SGLaue'] in ['m3m','m3']:
            return [[pfx+'A0',dpdA[0]]]
            
    # create a list of dependent variables and set up a dictionary to hold their derivatives
#    dependentVars = G2mv.GetDependentVars()
    depDerivDict = {}
    for j in dependentVars:
        depDerivDict[j] = np.zeros(shape=(len(x)))
#    print 'dependent vars',dependentVars
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    dMdv = np.zeros(shape=(len(varylist),len(x)))
    dMdb,dMddb,dMdpk = G2pwd.getBackgroundDerv(hfx,parmDict,bakType,calcControls[hfx+'histType'],x)
    if prc == 0 and hfx+'Back;0' in varylist: # for now assume that Back;x vars to not appear in constraints
        bBpos = varylist.index(hfx+'Back;0')
        dMdv[bBpos:bBpos+len(dMdb)] += dMdb     #TODO crash if bck parms tossed
    names = [hfx+'DebyeA',hfx+'DebyeR',hfx+'DebyeU']
    for name in varylist:
        if prc == 0 and 'Debye' in name:
            Id = int(name.split(';')[-1])
            parm = name[:int(name.rindex(';'))]
            ip = names.index(parm)
            dMdv[varylist.index(name)] += dMddb[3*Id+ip]
    names = [hfx+'BkPkpos',hfx+'BkPkint',hfx+'BkPksig',hfx+'BkPkgam']
    for name in varylist:
        if prc == 0 and 'BkPk' in name:
            parm,Id = name.split(';')
            Id = int(Id)
            if parm in names:
                ip = names.index(parm)
                dMdv[varylist.index(name)] += dMdpk[4*Id+ip]
    cw = np.diff(ma.getdata(x))
    cw = np.append(cw,cw[-1])
    Ka2 = False #also for TOF!
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.002)
        if hfx+'Lam1' in (parmDict.keys()):
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    elif 'B' in calcControls[hfx+'histType']:
        wave = parmDict[hfx+'Lam']
    for phase in Histogram['Reflection Lists']:
        refDict = Histogram['Reflection Lists'][phase]
        if phase not in Phases:     #skips deleted or renamed phases silently!
            continue
        Phase = Phases[phase]
        if histogram and histogram not in Phase['Histograms']:
            continue
        SGData = Phase['General']['SGData']
        SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
        im = 0
        if Phase['General'].get('Modulated',False):
            SSGData = Phase['General']['SSGData']
            im = 1  #offset in SS reflection list
            #??
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        Dij = GetDij(phfx,SGData,parmDict)
        A = [parmDict[pfx+'A%d'%(i)]+Dij[i] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        if not Phase['General'].get('doPawley') and not parmDict[phfx+'LeBail']:
            if im:
                dFdvDict = SStructureFactorDerv(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
                dFdvDict.update(SStructureFactorDerv2(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict))
            else:
                if Phase['General']['Type'] == 'magnetic': 
                    dFdvDict = MagStructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
                    dFdvDict.update(MagStructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict))
                else:
                    dFdvDict = StructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
            ApplyRBModelDervs(dFdvDict,parmDict,rigidbodyDict,Phase)
        # determine the parameters that will have derivatives computed only at end
        nonatomvarylist = []
        for name in varylist:
            if '::RBV;' not in name:
                try:
                    aname = name.split(pfx)[1][:2]
                    if aname not in ['Af','dA','AU','RB','AM','Xs','Xc','Ys','Yc','Zs','Zc',    \
                        'Tm','Xm','Ym','Zm','U1','U2','U3','MX','MY','MZ']: continue # skip anything not an atom or rigid body param
                except IndexError:
                    continue
            nonatomvarylist.append(name)
        nonatomdependentVars = []
        for name in dependentVars:
            if '::RBV;' not in name:
                try:
                    aname = name.split(pfx)[1][:2]
                    if aname not in ['Af','dA','AU','RB','AM','Xs','Xc','Ys','Yc','Zs','Zc',    \
                        'Tm','Xm','Ym','Zm','U1','U2','U3','MX','MY','MZ']: continue # skip anything not an atom or rigid body param
                except IndexError:
                    continue
            nonatomdependentVars.append(name)
        #==========================================================================================
        #==========================================================================================
        for iref in range(prc,len(refDict['RefList']),tprc):
            refl = refDict['RefList'][iref]
            if im:
                h,k,l,m = refl[:4]
            else:
                h,k,l = refl[:3]
            Uniq = np.inner(refl[:3],SGMT)
            if 'T' in calcControls[hfx+'histType']:
                wave = refl[14+im]
            dIdsh,dIdsp,dIdpola,dIdPO,dFdODF,dFdSA,dFdAb,dFdEx = GetIntensityDerv(refl,im,wave,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
            if 'C' in calcControls[hfx+'histType']:        #CW powder
                Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5+im],refl[6+im],refl[7+im],shl)
            elif 'T' in calcControls[hfx+'histType']:
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im])
            elif 'B' in calcControls[hfx+'histType']:
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.)
            iBeg = np.searchsorted(x,refl[5+im]-fmin)
            iFin = np.searchsorted(x,refl[5+im]+fmax)
            if not iBeg+iFin:       #peak below low limit - skip peak
                continue
            elif not iBeg-iFin:     #peak above high limit - done
                break
            pos = refl[5+im]
            if 'C' in calcControls[hfx+'histType']:
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                lenBF = iFin-iBeg
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdFCJVoigt3(refl[5+im],refl[6+im],refl[7+im],shl,ma.getdata(x[iBeg:iFin]))
                for i in range(5):
                    dMdpk[i] += 100.*cw[iBeg:iFin]*refl[11+im]*refl[9+im]*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':np.zeros_like(dMdpk[0])}
                if Ka2:
                    pos2 = refl[5+im]+lamRatio*tanth       # + 360/pi * Dlam/lam * tan(th)
                    iBeg2 = np.searchsorted(x,pos2-fmin)
                    iFin2 = np.searchsorted(x,pos2+fmax)
                    if iBeg2-iFin2:
                        lenBF2 = iFin2-iBeg2
                        dMdpk2 = np.zeros(shape=(6,lenBF2))
                        dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6+im],refl[7+im],shl,ma.getdata(x[iBeg2:iFin2]))
                        for i in range(5):
                            dMdpk2[i] = 100.*cw[iBeg2:iFin2]*refl[11+im]*refl[9+im]*kRatio*dMdipk2[i]
                        dMdpk2[5] = 100.*cw[iBeg2:iFin2]*refl[11+im]*dMdipk2[0]
                        dervDict2 = {'int':dMdpk2[0],'pos':dMdpk2[1],'sig':dMdpk2[2],'gam':dMdpk2[3],'shl':dMdpk2[4],'L1/L2':dMdpk2[5]*refl[9]}
            elif 'T' in calcControls[hfx+'histType']:
                lenBF = iFin-iBeg
                if lenBF < 0:   #bad peak coeff
                    break
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],ma.getdata(x[iBeg:iFin]))
                for i in range(6):
                    dMdpk[i] += refl[11+im]*refl[9+im]*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5]}            
            elif 'B' in calcControls[hfx+'histType']:
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                lenBF = iFin-iBeg
                if lenBF < 0:   #bad peak coeff
                    break
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im]/1.e4,refl[7+im]/100.,ma.getdata(x[iBeg:iFin]))
                for i in range(6):
                    dMdpk[i] = cw[iBeg:iFin]*refl[11+im]*refl[9+im]*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4]/1.e4,'gam':dMdpk[5]/100.}            
            if Phase['General'].get('doPawley'):
                dMdpw = np.zeros(len(x))
                try:
                    if im:
                        pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d,%d'%(h,k,l,m)])
                    else:
                        pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                    idx = varylist.index(pIdx)
                    dMdpw[iBeg:iFin] = dervDict['int']/refl[9+im]
                    if Ka2: #not for TOF either
                        dMdpw[iBeg2:iFin2] += dervDict2['int']/refl[9+im]
                    dMdv[idx] = dMdpw
                except: # ValueError:
                    pass
            if 'C' in calcControls[hfx+'histType']:
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY,dpdV = GetReflPosDerv(refl,im,wave,A,pfx,hfx,phfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],
                    hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'Z':[1.0,'gam'],hfx+'SH/L':[1.0,'shl'],
                    hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
                    hfx+'DisplaceY':[dpdY,'pos'],phfx+'Extinction':[dFdEx,'int'],phfx+'LayerDisp':[dpdY,'pos']}
                if 'Bragg' in calcControls[hfx+'instType']:
                    names.update({hfx+'SurfRoughA':[dFdAb[0],'int'],
                        hfx+'SurfRoughB':[dFdAb[1],'int'],})
                else:
                    names.update({hfx+'Absorption':[dFdAb,'int'],})
            elif 'T' in calcControls[hfx+'histType']:   #'T'OF
                dpdA,dpdZ,dpdDC,dpdDA,dpdDB,dpdV = GetReflPosDerv(refl,im,0.0,A,pfx,hfx,phfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'difC':[dpdDC,'pos'],hfx+'difA':[dpdDA,'pos'],hfx+'difB':[dpdDB,'pos'],
                    hfx+'Zero':[dpdZ,'pos'],hfx+'X':[refl[4+im],'gam'],hfx+'Y':[refl[4+im]**2,'gam'],hfx+'Z':[1.0,'gam'],
                    hfx+'alpha':[1./refl[4+im],'alp'],hfx+'beta-0':[1.0,'bet'],hfx+'beta-1':[1./refl[4+im]**4,'bet'],
                    hfx+'beta-q':[1./refl[4+im]**2,'bet'],hfx+'sig-0':[1.0,'sig'],hfx+'sig-1':[refl[4+im]**2,'sig'],
                    hfx+'sig-2':[refl[4+im]**4,'sig'],hfx+'sig-q':[refl[4+im],'sig'],
                    hfx+'Absorption':[dFdAb,'int'],phfx+'Extinction':[dFdEx,'int'],}
            elif 'B' in calcControls[hfx+'histType']:   #'T'OF
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY,dpdV = GetReflPosDerv(refl,im,0.0,A,pfx,hfx,phfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],phfx+'Scale':[dIdsp,'int'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Zero':[dpdZ,'pos'],hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'Z':[1.0,'gam'],
                    hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],hfx+'Polariz.':[dIdpola,'int'],
                    hfx+'alpha-0':[1.0,'alp'],hfx+'alpha-1':[tanth,'alp'],hfx+'beta-0':[1.0,'bet'],hfx+'beta-1':[tanth,'bet'],
                    hfx+'Absorption':[dFdAb,'int'],phfx+'Extinction':[dFdEx,'int'],}
            for name in names:
                item = names[name]
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += item[0]*dervDict[item[1]]
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += item[0]*dervDict[item[1]]
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
            for iPO in dIdPO:
                if iPO in varylist:
                    dMdv[varylist.index(iPO)][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(iPO)][iBeg2:iFin2] += dIdPO[iPO]*dervDict2['int']
                elif iPO in dependentVars:
                    depDerivDict[iPO][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[iPO][iBeg2:iFin2] += dIdPO[iPO]*dervDict2['int']
            for i,name in enumerate(['omega','chi','phi']):
                aname = pfx+'SH '+name
                if aname in varylist:
                    dMdv[varylist.index(aname)][iBeg:iFin] += dFdSA[i]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(aname)][iBeg2:iFin2] += dFdSA[i]*dervDict2['int']
                elif aname in dependentVars:
                    depDerivDict[aname][iBeg:iFin] += dFdSA[i]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[aname][iBeg2:iFin2] += dFdSA[i]*dervDict2['int']
            for iSH in dFdODF:
                if iSH in varylist:
                    dMdv[varylist.index(iSH)][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(iSH)][iBeg2:iFin2] += dFdODF[iSH]*dervDict2['int']
                elif iSH in dependentVars:
                    depDerivDict[iSH][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[iSH][iBeg2:iFin2] += dFdODF[iSH]*dervDict2['int']
            cellDervNames = cellVaryDerv(pfx,SGData,dpdA)
            for name,dpdA in cellDervNames:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dpdA*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += dpdA*dervDict2['pos']
                elif name in dependentVars: #need to scale for mixed phase constraints?
                    depDerivDict[name][iBeg:iFin] += dpdA*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += dpdA*dervDict2['pos']
            dDijDict = GetHStrainShiftDerv(refl,im,SGData,phfx,hfx,calcControls,parmDict)
            for name in dDijDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
#TODO: need Layer Disp deriv here
            for i,name in enumerate([pfx+'mV0',pfx+'mV1',pfx+'mV2']):
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dpdV[i]*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += dpdV[i]*dervDict2['pos']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += dpdV[i]*dervDict['pos']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += dpdV[i]*dervDict2['pos']
            if 'C' in calcControls[hfx+'histType'] or 'B' in calcControls[hfx+'histType']:
                sigDict,gamDict = GetSampleSigGamDerv(refl,im,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict)
            else:   #'T'OF
                sigDict,gamDict = GetSampleSigGamDerv(refl,im,0.0,G,GB,SGData,hfx,phfx,calcControls,parmDict)
            for name in gamDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += gamDict[name]*dervDict['gam']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += gamDict[name]*dervDict2['gam']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += gamDict[name]*dervDict['gam']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += gamDict[name]*dervDict2['gam']
            for name in sigDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += sigDict[name]*dervDict['sig']
                    if Ka2 and iFin2-iBeg2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += sigDict[name]*dervDict2['sig']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += sigDict[name]*dervDict['sig']
                    if Ka2 and iFin2-iBeg2:
                        depDerivDict[name][iBeg2:iFin2] += sigDict[name]*dervDict2['sig']
            for name in ['BabA','BabU']:
                if refl[9+im]:
                    if phfx+name in varylist:
                        dMdv[varylist.index(phfx+name)][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict['int']/refl[9+im]
                        if Ka2 and iFin2-iBeg2:
                            dMdv[varylist.index(phfx+name)][iBeg2:iFin2] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict2['int']/refl[9+im]
                    elif phfx+name in dependentVars:                    
                        depDerivDict[phfx+name][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict['int']/refl[9+im]
                        if Ka2 and iFin2-iBeg2:
                            depDerivDict[phfx+name][iBeg2:iFin2] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict2['int']/refl[9+im]                  
            if not Phase['General'].get('doPawley') and not parmDict[phfx+'LeBail']:
                #do atom derivatives -  for RB,F,X & U so far - how do I scale mixed phase constraints?
                corr = 0.
                corr2 = 0.
                if refl[9+im]:             
                    corr = dervDict['int']/refl[9+im]
                    if Ka2 and iFin2-iBeg2:
                        corr2 = dervDict2['int']/refl[9+im]
                for name in nonatomvarylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dFdvDict[name][iref]*corr
                    if Ka2 and iFin2-iBeg2:
                       dMdv[varylist.index(name)][iBeg2:iFin2] += dFdvDict[name][iref]*corr2
                for name in nonatomdependentVars:
                   depDerivDict[name][iBeg:iFin] += dFdvDict[name][iref]*corr
                   if Ka2 and iFin2-iBeg2:
                       depDerivDict[name][iBeg2:iFin2] += dFdvDict[name][iref]*corr2
    # now process derivatives in constraints
    dMdv[:,ma.getmaskarray(x)] = 0.  # instead of masking, zero out masked values
    #G2mv.Dict2Deriv(varylist,depDerivDict,dMdv)
    return dMdv,depDerivDict
    
def UserRejectHKL(ref,im,userReject):
    if ref[5+im]/ref[6+im] < userReject['minF/sig']:
        return False
    elif userReject['MaxD'] < ref[4+im] > userReject['MinD']:
        return False
    elif ref[11+im] < userReject['MinExt']:
        return False
    elif abs(ref[5+im]-ref[7+im])/ref[6+im] > userReject['MaxDF/F']:
        return False
    return True
    
def dervHKLF(Histogram,Phase,calcControls,varylist,parmDict,rigidbodyDict):
    '''Loop over reflections in a HKLF histogram and compute derivatives of the fitting
    model (M) with respect to all parameters.  Independent and dependant dM/dp arrays 
    are returned to either dervRefine or HessRefine.

    :returns: 
    '''
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    pfx = '%d::'%(Phase['pId'])
    phfx = '%d:%d:'%(Phase['pId'],hId)
    SGData = Phase['General']['SGData']
    im = 0
    if Phase['General'].get('Modulated',False):
        SSGData = Phase['General']['SSGData']
        im = 1  #offset in SS reflection list
    A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
    G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
    TwinLaw = calcControls[phfx+'TwinLaw']
    refDict = Histogram['Data']
    if parmDict[phfx+'Scale'] < 0.:
        parmDict[phfx+'Scale'] = .001
    if im: # split to nontwin/twin versions
        if len(TwinLaw) > 1:
            dFdvDict = SStructureFactorDervTw(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)    #not ok
        else:
            dFdvDict = SStructureFactorDerv(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)  #OK
            dFdvDict.update(SStructureFactorDerv2(refDict,im,G,hfx,pfx,SGData,SSGData,calcControls,parmDict))
    else:
        if len(TwinLaw) > 1:
            dFdvDict = StructureFactorDervTw2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
        else:   #correct!!
            if Phase['General']['Type'] == 'magnetic': 
                dFdvDict = MagStructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
                dFdvDict.update(MagStructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict))
            else:
                dFdvDict = StructureFactorDerv2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
    ApplyRBModelDervs(dFdvDict,parmDict,rigidbodyDict,Phase)
    dMdvh = np.zeros((len(varylist),len(refDict['RefList'])))
    dependentVars = G2mv.GetDependentVars()
    depDerivDict = {}
    for j in dependentVars:
        depDerivDict[j] = np.zeros(shape=(len(refDict['RefList'])))
    wdf = np.zeros(len(refDict['RefList']))
    if calcControls['F**2']:
        for iref,ref in enumerate(refDict['RefList']):
            if ref[6+im] > 0:
                dervDict,dervCor = SCExtinction(ref,im,phfx,hfx,pfx,calcControls,parmDict,varylist+dependentVars)[1:]
                w = 1.0/ref[6+im]
                if ref[3+im] > 0:
                    wdf[iref] = w*(ref[5+im]-ref[7+im])
                    for j,var in enumerate(varylist):
                        if var in dFdvDict:
                            dMdvh[j][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11+im]
                    for var in dependentVars:
                        if var in dFdvDict:
                            depDerivDict[var][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11+im]
                    if phfx+'Scale' in varylist:
                        dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[7+im]*ref[11+im]/parmDict[phfx+'Scale']  #OK
                    elif phfx+'Scale' in dependentVars:
                        depDerivDict[phfx+'Scale'][iref] = w*ref[7+im]*ref[11+im]/parmDict[phfx+'Scale']   #OK
                    for item in ['Ep','Es','Eg']:
                        if phfx+item in varylist and phfx+item in dervDict:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]/ref[11+im]  #OK
                        elif phfx+item in dependentVars and phfx+item in dervDict:
                            depDerivDict[phfx+item][iref] = w*dervDict[phfx+item]/ref[11+im]  #OK
                    for item in ['BabA','BabU']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dFdvDict[phfx+item][iref]*parmDict[phfx+'Scale']*ref[11+im]
                        elif phfx+item in dependentVars:
                            depDerivDict[phfx+item][iref] = w*dFdvDict[phfx+item][iref]*parmDict[phfx+'Scale']*ref[11+im]
    else:   #F refinement
        for iref,ref in enumerate(refDict['RefList']):
            if ref[5+im] > 0.:
                dervDict,dervCor = SCExtinction(ref,im,phfx,hfx,pfx,calcControls,parmDict,varylist+dependentVars)[1:]
                Fo = np.sqrt(ref[5+im])
                Fc = np.sqrt(ref[7+im])
                w = 1.0/ref[6+im]
                if ref[3+im] > 0:
                    wdf[iref] = 2.0*Fc*w*(Fo-Fc)
                    for j,var in enumerate(varylist):
                        if var in dFdvDict:
                            dMdvh[j][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11+im]
                    for var in dependentVars:
                        if var in dFdvDict:
                            depDerivDict[var][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11+im]
                    if phfx+'Scale' in varylist:
                        dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[7+im]*ref[11+im]/parmDict[phfx+'Scale']  #OK
                    elif phfx+'Scale' in dependentVars:
                        depDerivDict[phfx+'Scale'][iref] = w*ref[7+im]*ref[11+im]/parmDict[phfx+'Scale']   #OK                    
                    for item in ['Ep','Es','Eg']:   #OK!
                        if phfx+item in varylist and phfx+item in dervDict:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]/ref[11+im] 
                        elif phfx+item in dependentVars and phfx+item in dervDict:
                            depDerivDict[phfx+item][iref] = w*dervDict[phfx+item]/ref[11+im]
                    for item in ['BabA','BabU']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dFdvDict[phfx+item][iref]*parmDict[phfx+'Scale']*ref[11+im]
                        elif phfx+item in dependentVars:
                            depDerivDict[phfx+item][iref] = w*dFdvDict[phfx+item][iref]*parmDict[phfx+'Scale']*ref[11+im]
    return dMdvh,depDerivDict,wdf
    

def dervRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,dlg):
    '''Loop over histograms and compute derivatives of the fitting
    model (M) with respect to all parameters.  Results are returned in
    a Jacobian matrix (aka design matrix) of dimensions (n by m) where
    n is the number of parameters and m is the number of data
    points. This can exceed memory when m gets large. This routine is
    used when refinement derivatives are selected as "analtytic
    Jacobian" in Controls.

    :returns: Jacobian numpy.array dMdv for all histograms concatinated
    '''
    parmDict.update(zip(varylist,values))
    G2mv.Dict2Map(parmDict,varylist)
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    dependentVars = G2mv.GetDependentVars()
    histoList = list(Histograms.keys())
    histoList.sort()
    First = True
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])+1
            dMdv,depDerivDict = getPowderProfileDervMP([parmDict,x[xB:xF],
                varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars])
            G2mv.Dict2Deriv(varylist,depDerivDict,dMdv)
            dMdvh = np.sqrt(w[xB:xF])*dMdv
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            dMdvh,depDerivDict,wdf = dervHKLF(Histogram,Phase,calcControls,varylist,parmDict,rigidbodyDict)
            hfx = ':%d:'%(Histogram['hId'])
            wtFactor = calcControls[hfx+'wtFactor']
            # now process derivatives in constraints
            G2mv.Dict2Deriv(varylist,depDerivDict,dMdvh)
        else:
            continue        #skip non-histogram entries
        if First:
            dMdV = np.sqrt(wtFactor)*dMdvh
            First = False
        else:
            dMdV = np.concatenate((dMdV.T,np.sqrt(wtFactor)*dMdvh.T)).T
            
    GetFobsSq(Histograms,Phases,parmDict,calcControls)
    pNames,pVals,pWt,pWsum,pWnum = penaltyFxn(HistoPhases,calcControls,parmDict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,calcControls,parmDict,varylist)
        dMdV = np.concatenate((dMdV.T,(np.sqrt(pWt)*dpdv).T)).T
        
    return dMdV

def HessRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,dlg):
    '''Loop over histograms and compute derivatives of the fitting
    model (M) with respect to all parameters.  For each histogram, the
    Jacobian matrix, dMdv, with dimensions (n by m) where n is the
    number of parameters and m is the number of data points *in the
    histogram*. The (n by n) Hessian is computed from each Jacobian
    and it is returned.  This routine is used when refinement
    derivatives are selected as "analtytic Hessian" in Controls.

    :returns: Vec,Hess where Vec is the least-squares vector and Hess is the Hessian
    '''
    parmDict.update(zip(varylist,values))
    G2mv.Dict2Map(parmDict,varylist)
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    dependentVars = G2mv.GetDependentVars()
    #fixup H atom positions here?
    ApplyRBModels(parmDict,Phases,rigidbodyDict)        #,Update=True??
    Hess = np.empty(0)
    Vec = np.empty(0)
    histoList = list(Histograms.keys())
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            W = wtFactor*w
            dy = y-yc
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])+1
            useMP,ncores = G2mp.InitMP()
            if GSASIIpath.GetConfigValue('Show_timing',False): starttime = time.time()
            if useMP:
                MPpool = mp.Pool(ncores)
                dMdvh = None
                depDerivDict = None
                # old approach, create all args prior to use
#                profArgs = [
#                    (parmDict,x[xB:xF],varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars,
#                     i,ncores,histogram) for i in range(ncores)]
#                for dmdv,depDerivs in MPpool.imap_unordered(getPowderProfileDervMP,profArgs):
                # better, use a generator so arg is created as used
                profGenArgs = (
                    (parmDict,x[xB:xF],varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars,
                     i,ncores,histogram) for i in range(ncores))
                for dmdv,depDerivs in MPpool.imap_unordered(getPowderProfileDervMP,profGenArgs):
                    if dMdvh is None:
                       dMdvh = dmdv
                       depDerivDict = depDerivs
                    else: 
                       dMdvh += dmdv
                       for key in depDerivs.keys(): depDerivDict[key] += depDerivs[key]                        
                MPpool.terminate()
            else:
                dMdvh,depDerivDict = getPowderProfileDervMP([parmDict,x[xB:xF],
                    varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars,0,1,histogram])
                #dMdvh = getPowderProfileDerv(parmDict,x[xB:xF],
                #    varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup,dependentVars)
            G2mv.Dict2Deriv(varylist,depDerivDict,dMdvh)
            if GSASIIpath.GetConfigValue('Show_timing',False): print ('getPowderProfileDerv t=%.3f'%(time.time()-starttime))
            Wt = ma.sqrt(W[xB:xF])[nxs,:]
            Dy = dy[xB:xF][nxs,:]
            dMdvh *= Wt
            if dlg:
                GoOn = dlg.Update(Histogram['Residuals']['wR'],newmsg='Hessian for histogram %d\nAll data Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))
                if type(GoOn) is tuple:
                    if not GoOn[0]:
                        raise G2obj.G2RefineCancel('Cancel pressed')
                elif not GoOn:
                    raise G2obj.G2RefineCancel('Cancel pressed')
                dlg.Raise()
            if len(Hess):
                Hess += np.inner(dMdvh,dMdvh)
                dMdvh *= Wt*Dy
                Vec += np.sum(dMdvh,axis=1)
            else:
                Hess = np.inner(dMdvh,dMdvh)
                dMdvh *= Wt*Dy
                Vec = np.sum(dMdvh,axis=1)
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            dMdvh,depDerivDict,wdf = dervHKLF(Histogram,Phase,calcControls,varylist,parmDict,rigidbodyDict)
            hId = Histogram['hId']
            hfx = ':%d:'%(Histogram['hId'])
            wtFactor = calcControls[hfx+'wtFactor']
            # now process derivatives in constraints
            G2mv.Dict2Deriv(varylist,depDerivDict,dMdvh)
#            print 'matrix build time: %.3f'%(time.time()-time0)

            if dlg:
                GoOn = dlg.Update(Histogram['Residuals']['wR'],newmsg='Hessian for histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))
                if type(GoOn) is tuple:
                    if not GoOn[0]:
                        raise G2obj.G2RefineCancel('Cancel pressed')
                elif not GoOn:
                    raise G2obj.G2RefineCancel('Cancel pressed')
                dlg.Raise()
            if len(Hess):
                Vec += wtFactor*np.sum(dMdvh*wdf,axis=1)
                Hess += wtFactor*np.inner(dMdvh,dMdvh)
            else:
                Vec = wtFactor*np.sum(dMdvh*wdf,axis=1)
                Hess = wtFactor*np.inner(dMdvh,dMdvh)
        else:
            continue        #skip non-histogram entries
    GetFobsSq(Histograms,Phases,parmDict,calcControls)
    pNames,pVals,pWt,pWsum,pWnum = penaltyFxn(HistoPhases,calcControls,parmDict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,calcControls,parmDict,varylist)
        Vec -= np.sum(dpdv*pWt*pVals,axis=1)
        Hess += np.inner(dpdv*pWt,dpdv)
    return Vec,Hess

def errRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,dlg=None):        
    '''Computes the point-by-point discrepancies between every data point in every histogram
    and the observed value. Used in the Jacobian, Hessian & numeric least-squares to compute function
    
    :returns: an np array of differences between observed and computed diffraction values.
    '''
    Values2Dict(parmDict, varylist, values)
    G2mv.Dict2Map(parmDict,varylist)
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    M = np.empty(0)
    SumwYo = 0
    Nobs = 0
    Nrej = 0
    Next = 0
    ApplyRBModels(parmDict,Phases,rigidbodyDict)
    #fixup Hatom positions here....
    histoList = list(Histograms.keys())
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            yc *= 0.0                           #zero full calcd profiles
            yb *= 0.0
            yd *= 0.0
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])+1
            yc[xB:xF],yb[xB:xF] = getPowderProfile(parmDict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup,histogram)
            yc[xB:xF] += yb[xB:xF]
            if not np.any(y):                   #fill dummy data
                try:
                    rv = st.poisson(yc[xB:xF])
                    y[xB:xF] = rv.rvs()
                except ValueError:
                    y[xB:xF] = yc[xB:xF]
                Z = np.ones_like(yc[xB:xF])
                Z[1::2] *= -1
                y[xB:xF] = yc[xB:xF]+np.abs(y[xB:xF]-yc[xB:xF])*Z
                w[xB:xF] = np.where(y[xB:xF]>0.,1./y[xB:xF],1.0)
            yd[xB:xF] = y[xB:xF]-yc[xB:xF]
            W = wtFactor*w
            wdy = -ma.sqrt(w[xB:xF])*(yd[xB:xF])
            Histogram['Residuals']['Durbin-Watson'] = ma.sum(ma.diff(wdy)**2)/ma.sum(wdy**2)
            wdy *= np.sqrt(wtFactor)
            Histogram['Residuals']['Nobs'] = ma.count(x[xB:xF])
            Nobs += Histogram['Residuals']['Nobs']
            Histogram['Residuals']['sumwYo'] = ma.sum(W[xB:xF]*y[xB:xF]**2)
            SumwYo += Histogram['Residuals']['sumwYo']
            Histogram['Residuals']['R'] = min(100.,ma.sum(ma.abs(yd[xB:xF]))/ma.sum(y[xB:xF])*100.)
            Histogram['Residuals']['wR'] = min(100.,ma.sqrt(ma.sum(wdy**2)/Histogram['Residuals']['sumwYo'])*100.)
            sumYmB = ma.sum(ma.where(yc[xB:xF]!=yb[xB:xF],ma.abs(y[xB:xF]-yb[xB:xF]),0.))
            sumwYmB2 = ma.sum(ma.where(yc[xB:xF]!=yb[xB:xF],W[xB:xF]*(y[xB:xF]-yb[xB:xF])**2,0.))
            sumYB = ma.sum(ma.where(yc[xB:xF]!=yb[xB:xF],ma.abs(y[xB:xF]-yc[xB:xF])*ma.abs(y[xB:xF]-yb[xB:xF])/y[xB:xF],0.))
            sumwYB2 = ma.sum(ma.where(yc[xB:xF]!=yb[xB:xF],W[xB:xF]*(ma.abs(y[xB:xF]-yc[xB:xF])*ma.abs(y[xB:xF]-yb[xB:xF])/y[xB:xF])**2,0.))
            Histogram['Residuals']['Rb'] = min(100.,100.*sumYB/sumYmB)
            Histogram['Residuals']['wRb'] = min(100.,100.*ma.sqrt(sumwYB2/sumwYmB2))
            Histogram['Residuals']['wRmin'] = min(100.,100.*ma.sqrt(Histogram['Residuals']['Nobs']/Histogram['Residuals']['sumwYo']))
            if dlg:
                GoOn = dlg.Update(Histogram['Residuals']['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))
                if type(GoOn) is tuple:
                    if not GoOn[0]:
                        raise G2obj.G2RefineCancel('Cancel pressed')
                elif not GoOn:
                    raise G2obj.G2RefineCancel('Cancel pressed')
                dlg.Raise()
            M = np.concatenate((M,wdy))
#end of PWDR processing
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            Histogram['Residuals'] = {}
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            pfx = '%d::'%(Phase['pId'])
            phfx = '%d:%d:'%(Phase['pId'],hId)
            SGData = Phase['General']['SGData']
            TwinLaw = calcControls[phfx+'TwinLaw']
            im = 0
            if parmDict[phfx+'Scale'] < 0.:
                parmDict[phfx+'Scale'] = .001                
            if Phase['General'].get('Modulated',False):
                SSGData = Phase['General']['SSGData']
                im = 1  #offset in SS reflection list
            A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
            G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
            refDict = Histogram['Data']
            if im:
                if len(TwinLaw) > 1:
                    SStructureFactorTw(refDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
                else:
                    SStructureFactor(refDict,G,hfx,pfx,SGData,SSGData,calcControls,parmDict)
            else:
                StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
#            print 'sf-calc time: %.3f'%(time.time()-time0)
            df = np.zeros(len(refDict['RefList']))
            sumwYo = 0
            sumFo = 0
            sumFo2 = 0
            sumFc2 = 0
            sumdF = 0
            sumdF2 = 0
            if im:
                sumSSFo = np.zeros(10)
                sumSSFo2 = np.zeros(10)
                sumSSdF = np.zeros(10)
                sumSSdF2 = np.zeros(10)
                sumSSwYo = np.zeros(10)
                sumSSwdf2 = np.zeros(10)
                SSnobs = np.zeros(10)
            nobs = 0
            nrej = 0
            next = 0
            maxH = 0
            if calcControls['F**2']:
                for i,ref in enumerate(refDict['RefList']):
                    if ref[6+im] > 0:
                        ref[11+im] = SCExtinction(ref,im,phfx,hfx,pfx,calcControls,parmDict,varylist)[0]
                        w = 1.0/ref[6+im]   # 1/sig(F^2)
                        ref[7+im] *= parmDict[phfx+'Scale']*ref[11+im]  #correct Fc^2 for extinction
                        ref[8+im] = ref[5+im]/(parmDict[phfx+'Scale']*ref[11+im])
                        if UserRejectHKL(ref,im,calcControls['UsrReject']) and ref[3+im]:    #skip sp.gp. absences (mul=0)
                            ref[3+im] = abs(ref[3+im])      #mark as allowed
                            Fo = np.sqrt(ref[5+im])
                            sumFo += Fo
                            sumFo2 += ref[5+im]
                            sumFc2 += ref[7+im]
                            sumdF += abs(Fo-np.sqrt(ref[7+im]))
                            sumdF2 += abs(ref[5+im]-ref[7+im])
                            nobs += 1
                            df[i] = -w*(ref[5+im]-ref[7+im])
                            sumwYo += (w*ref[5+im])**2      #w*Fo^2
                            if im:  #accumulate super lattice sums
                                ind = int(abs(ref[3]))
                                sumSSFo[ind] += Fo
                                sumSSFo2[ind] += ref[5+im]
                                sumSSdF[ind] += abs(Fo-np.sqrt(ref[7+im]))
                                sumSSdF2[ind] += abs(ref[5+im]-ref[7+im])
                                sumSSwYo[ind] += (w*ref[5+im])**2      #w*Fo^2
                                sumSSwdf2[ind] +=  df[i]**2
                                SSnobs[ind] += 1
                                maxH = max(maxH,ind)                           
                        else:
                            if ref[3+im]:
                                ref[3+im] = -abs(ref[3+im])      #mark as rejected
                                nrej += 1
                            else:   #sp.gp.extinct
                                next += 1
            else:
                for i,ref in enumerate(refDict['RefList']):
                    if ref[5+im] > 0.:
                        ref[11+im] = SCExtinction(ref,im,phfx,hfx,pfx,calcControls,parmDict,varylist)[0]
                        ref[7+im] *= parmDict[phfx+'Scale']*ref[11+im]    #correct Fc^2 for extinction
                        ref[8+im] = ref[5+im]/(parmDict[phfx+'Scale']*ref[11+im])
                        Fo = np.sqrt(ref[5+im])
                        Fc = np.sqrt(ref[7+im])
                        w = 2.0*Fo/ref[6+im]    # 1/sig(F)?
                        if UserRejectHKL(ref,im,calcControls['UsrReject']) and ref[3+im]:    #skip sp.gp. absences (mul=0)
                            ref[3+im] = abs(ref[3+im])      #mark as allowed
                            sumFo += Fo
                            sumFo2 += ref[5+im]
                            sumFc2 += ref[7+im]
                            sumdF += abs(Fo-Fc)
                            sumdF2 += abs(ref[5+im]-ref[7+im])
                            nobs += 1
                            df[i] = -w*(Fo-Fc)
                            sumwYo += (w*Fo)**2
                            if im:
                                ind = int(abs(ref[3]))
                                sumSSFo[ind] += Fo
                                sumSSFo2[ind] += ref[5+im]
                                sumSSdF[ind] += abs(Fo-Fc)
                                sumSSdF2[ind] += abs(ref[5+im]-ref[7+im])
                                sumSSwYo[ind] += (w*Fo)**2                                                            
                                sumSSwdf2[ind] +=  df[i]**2
                                SSnobs[ind] += 1                           
                                maxH = max(maxH,ind)                           
                        else:
                            if ref[3+im]:
                                ref[3+im] = -abs(ref[3+im])      #mark as rejected
                                nrej += 1
                            else:   #sp.gp.extinct
                                next += 1
            Scale = sumFo2/sumFc2
            if (Scale < 0.8 or Scale > 1.2) and phfx+'Scale' in varylist:
                print ('New scale: %.4f'%(Scale*parmDict[phfx+'Scale']))
                indx = varylist.index(phfx+'Scale')
                values[indx] = Scale*parmDict[phfx+'Scale']              
            Histogram['Residuals']['Nobs'] = nobs
            Histogram['Residuals']['sumwYo'] = sumwYo
            SumwYo += sumwYo
            Histogram['Residuals']['wR'] = min(100.,np.sqrt(np.sum(df**2)/sumwYo)*100.)
            Histogram['Residuals'][phfx+'Rf'] = 100.*sumdF/sumFo
            Histogram['Residuals'][phfx+'Rf^2'] = 100.*sumdF2/sumFo2
            Histogram['Residuals'][phfx+'Nref'] = nobs
            Histogram['Residuals'][phfx+'Nrej'] = nrej
            Histogram['Residuals'][phfx+'Next'] = next
            if im:
                Histogram['Residuals'][phfx+'SSRf'] = 100.*sumSSdF[:maxH+1]/sumSSFo[:maxH+1]
                Histogram['Residuals'][phfx+'SSRf^2'] = 100.*sumSSdF2[:maxH+1]/sumSSFo2[:maxH+1]
                Histogram['Residuals'][phfx+'SSNref'] = SSnobs[:maxH+1]
                Histogram['Residuals']['SSwR'] = np.sqrt(sumSSwdf2[:maxH+1]/sumSSwYo[:maxH+1])*100.                
            Nobs += nobs
            Nrej += nrej
            Next += next
            if dlg:
                GoOn = dlg.Update(Histogram['Residuals']['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))
                if type(GoOn) is tuple:
                    if not GoOn[0]:
                        raise G2obj.G2RefineCancel('Cancel pressed')
                elif not GoOn:
                    raise G2obj.G2RefineCancel('Cancel pressed')
                dlg.Raise()
            M = np.concatenate((M,wtFactor*df))
# end of HKLF processing
#    GetFobsSq(Histograms,Phases,parmDict,calcControls)
    Histograms['sumwYo'] = SumwYo
    Histograms['Nobs'] = Nobs
    Histograms['Nrej'] = Nrej
    Histograms['Next'] = Next
    Rw = min(100.,np.sqrt(np.sum(M**2)/SumwYo)*100.)
    if dlg:
        GoOn = dlg.Update(Rw,newmsg='%s%8.3f%s'%('All data Rw =',Rw,'%'))
        if type(GoOn) is tuple:
            if not GoOn[0]:
                parmDict['saved values'] = values
                raise G2obj.G2RefineCancel('Cancel pressed')
        elif not GoOn:
            parmDict['saved values'] = values
            raise G2obj.G2RefineCancel('Cancel pressed')
        dlg.Raise()
    pDict,pVals,pWt,pWsum,pWnum = penaltyFxn(HistoPhases,calcControls,parmDict,varylist)
    if len(pVals):
        pSum = np.sum(pWt*pVals**2)
        for name in pWsum:
            if pWsum[name]:
                print ('  Penalty function for %5d %8ss = %12.5g'%(pWnum[name],name,pWsum[name]))
        print ('Total penalty function: %12.5g on %d terms'%(pSum,len(pVals)))
        Nobs += len(pVals)
        M = np.concatenate((M,np.sqrt(pWt)*pVals))
    return M
