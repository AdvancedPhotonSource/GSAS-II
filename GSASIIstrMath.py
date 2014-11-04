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
import time
import math
import copy
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import scipy.optimize as so
import scipy.stats as st
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIpwd as G2pwd
import GSASIImapvars as G2mv
import GSASIImath as G2mth

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ateln2 = 8.0*math.log(2.0)

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
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        AtLookup = G2mth.FillAtomLookUp(Phase['Atoms'])
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
    'Needs a doc string'
    atxIds = ['dAx:','dAy:','dAz:']
    atuIds = ['AU11:','AU22:','AU33:','AU12:','AU13:','AU23:']
    TIds = ['T11:','T22:','T33:','T12:','T13:','T23:']
    LIds = ['L11:','L22:','L33:','L12:','L13:','L23:']
    SIds = ['S12:','S13:','S21:','S23:','S31:','S32:','SAA:','SBB:']
    PIds = ['Px:','Py:','Pz:']
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
    cell = General['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(cell)
    rpd = np.pi/180.
    rpd2 = rpd**2
    g = nl.inv(np.inner(Bmat,Bmat))
    gvec = np.sqrt(np.array([g[0][0]**2,g[1][1]**2,g[2][2]**2,
        g[0][0]*g[1][1],g[0][0]*g[2][2],g[1][1]*g[2][2]]))
    AtLookup = G2mth.FillAtomLookUp(Phase['Atoms'])
    pfx = str(Phase['pId'])+'::'
    RBModels =  Phase['RBModels']
    for irb,RBObj in enumerate(RBModels.get('Vector',[])):
        VModel = RBData['Vector'][RBObj['RBId']]
        Q = RBObj['Orient'][0]
        Pos = RBObj['Orig'][0]
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
                XYZ1 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q))
                Q[iv] += 2.*dx
                XYZ2 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q))
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
        Q = RBObj['Orient'][0]
        Pos = RBObj['Orig'][0]
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
                XYZ1 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q))
                Q[iv] += 2.*dx
                XYZ2 = G2mth.RotateRBXYZ(Bmat,Cart,G2mth.normQ(Q))
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

def penaltyFxn(HistoPhases,parmDict,varyList):
    'Needs a doc string'
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    pNames = []
    pVals = []
    pWt = []
    negWt = {}
    pWsum = {}
    for phase in Phases:
        pId = Phases[phase]['pId']
        negWt[pId] = Phases[phase]['General']['Pawley neg wt']
        General = Phases[phase]['General']
        textureData = General['SH Texture']
        SGData = General['SGData']
        AtLookup = G2mth.FillAtomLookUp(Phases[phase]['Atoms'])
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        if phase not in restraintDict:
            continue
        phaseRest = restraintDict[phase]
        names = [['Bond','Bonds'],['Angle','Angles'],['Plane','Planes'],
            ['Chiral','Volumes'],['Torsion','Torsions'],['Rama','Ramas'],
            ['ChemComp','Sites'],['Texture','HKLs']]
        for name,rest in names:
            pWsum[name] = 0.
            itemRest = phaseRest[name]
            if itemRest[rest] and itemRest['Use']:
                wt = itemRest['wtFactor']
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
                elif name == 'ChemComp':
                    for i,[indx,factors,obs,esd] in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs+1))
                        frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs-1))
                        calc = np.sum(mul*frac*factors)
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)                    
                        pWsum[name] += wt*((obs-calc)/esd)**2
                elif name == 'Texture':
                    SHkeys = textureData['SH Coeff'][1].keys()
                    SHCoef = G2mth.GetSHCoeff(pId,parmDict,SHkeys)
                    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
                    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
                    for i,[hkl,grid,esd1,ifesd2,esd2] in enumerate(itemRest[rest]):
                        PH = np.array(hkl)
                        phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                        ODFln = G2lat.Flnh(False,SHCoef,phi,beta,SGData)
                        R,P,Z = G2mth.getRestPolefig(ODFln,SamSym[textureData['Model']],grid)
                        Z1 = -ma.masked_greater(Z,0.0)
                        IndZ1 = np.array(ma.nonzero(Z1))
                        for ind in IndZ1.T:
                            pNames.append('%d:%s:%d:%.2f:%.2f'%(pId,name,i,R[ind[0],ind[1]],P[ind[0],ind[1]]))
                            pVals.append(Z1[ind[0]][ind[1]])
                            pWt.append(wt/esd1**2)
                            pWsum[name] += wt*((obs-calc)/esd)**2
                        if ifesd2:
                            Z2 = 1.-Z
                            for ind in np.ndindex(grid,grid):
                                pNames.append('%d:%s:%d:%.2f:%.2f'%(pId,name+'-unit',i,R[ind[0],ind[1]],P[ind[0],ind[1]]))
                                pVals.append(Z1[ind[0]][ind[1]])
                                pWt.append(wt/esd2**2)
                                pWsum[name] += wt*((obs-calc)/esd)**2
         
    pWsum['PWLref'] = 0.
    for item in varyList:
        if 'PWLref' in item and parmDict[item] < 0.:
            pId = int(item.split(':')[0])
            if negWt[pId]:
                pNames.append(item)
                pVals.append(-parmDict[item])
                pWt.append(negWt[pId])
                pWsum['PWLref'] += negWt[pId]*(-parmDict[item])**2
    pVals = np.array(pVals)
    pWt = np.array(pWt)         #should this be np.sqrt?
    return pNames,pVals,pWt,pWsum
    
def penaltyDeriv(pNames,pVal,HistoPhases,parmDict,varyList):
    'Needs a doc string'
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    pDerv = np.zeros((len(varyList),len(pVal)))
    for phase in Phases:
#        if phase not in restraintDict:
#            continue
        pId = Phases[phase]['pId']
        General = Phases[phase]['General']
        SGData = General['SGData']
        AtLookup = G2mth.FillAtomLookUp(Phases[phase]['Atoms'])
        cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(cell)
        textureData = General['SH Texture']

        SHkeys = textureData['SH Coeff'][1].keys()
        SHCoef = G2mth.GetSHCoeff(pId,parmDict,SHkeys)
        shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
        SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
        sam = SamSym[textureData['Model']]
        phaseRest = restraintDict.get(phase,{})
        names = {'Bond':'Bonds','Angle':'Angles','Plane':'Planes',
            'Chiral':'Volumes','Torsion':'Torsions','Rama':'Ramas',
            'ChemComp':'Sites','Texture':'HKLs'}
        lasthkl = np.array([0,0,0])
        for ip,pName in enumerate(pNames):
            pnames = pName.split(':')
            if pId == int(pnames[0]):
                name = pnames[1]
                if 'PWL' in pName:
                    pDerv[varyList.index(pName)][ip] += 1.
                    continue
                id = int(pnames[2]) 
                itemRest = phaseRest[name]
                if name in ['Bond','Angle','Plane','Chiral']:
                    indx,ops,obs,esd = itemRest[names[name]][id]
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
                    indx,ops,cofName,esd = itemRest[names[name]][id]
                    dNames = []
                    for ind in indx:
                        dNames += [str(pId)+'::dA'+Xname+':'+str(AtLookup[ind]) for Xname in ['x','y','z']]
                    XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                    if name == 'Torsion':
                        deriv = G2mth.getTorsionDeriv(XYZ,Amat,coffDict[cofName])
                    else:
                        deriv = G2mth.getRamaDeriv(XYZ,Amat,coffDict[cofName])
                elif name == 'ChemComp':
                    indx,factors,obs,esd = itemRest[names[name]][id]
                    dNames = []
                    for ind in indx:
                        dNames += [str(pId)+'::Afrac:'+str(AtLookup[ind])]
                        mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs+1))
                        deriv = mul*factors
                elif 'Texture' in name:
                    deriv = []
                    dNames = []
                    hkl,grid,esd1,ifesd2,esd2 = itemRest[names[name]][id]
                    hkl = np.array(hkl)
                    if np.any(lasthkl-hkl):
                        PH = np.array(hkl)
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
    FFdata = []
    BLdata = []
    Xdata = np.zeros((3,Natoms))
    dXdata = np.zeros((3,Natoms))
    Uisodata = np.zeros(Natoms)
    Uijdata = np.zeros((6,Natoms))
    keys = {'Atype:':Tdata,'Amul:':Mdata,'Afrac:':Fdata,'AI/A:':IAdata,
        'dAx:':dXdata[0],'dAy:':dXdata[1],'dAz:':dXdata[2],
        'Ax:':Xdata[0],'Ay:':Xdata[1],'Az:':Xdata[2],'AUiso:':Uisodata,
        'AU11:':Uijdata[0],'AU22:':Uijdata[1],'AU33:':Uijdata[2],
        'AU12:':Uijdata[3],'AU13:':Uijdata[4],'AU23:':Uijdata[5]}
    for iatm in range(Natoms):
        for key in keys:
            parm = pfx+key+str(iatm)
            if parm in parmDict:
                keys[key][iatm] = parmDict[parm]
    Fdata = np.where(Fdata,Fdata,1.e-8)         #avoid divide by zero in derivative calc.?
    return Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata
    
def StructureFactor(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Not Used: here only for comparison the StructureFactor2 - faster version
    Compute structure factors for all h,k,l for phase
    puts the result, F^2, in each ref[8] in refList
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:

    '''        
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,calcControls,parmDict)
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    else:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    if not len(refDict['FF']):
        if 'N' in calcControls[hfx+'histType']:
            dat = G2el.getBLvalues(BLtables)        #will need wave here for anom. neutron b's
        else:
            dat = G2el.getFFvalues(FFtables,0.)        
        refDict['FF']['El'] = dat.keys()
        refDict['FF']['FF'] = np.zeros((len(refDict['RefList']),len(dat)))   
    for iref,refl in enumerate(refDict['RefList']):
        if 'NT' in calcControls[hfx+'histType']:
            FP,FPP = G2el.BlenResCW(Tdata,BLtables,refl[14])
        fbs = np.array([0,0])
        H = refl[:3]
        SQ = 1./(2.*refl[4])**2
        SQfactor = 4.0*SQ*twopisq
        Bab = parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor)
        if not np.any(refDict['FF']['FF'][iref]):                #no form factors - 1st time thru StructureFactor
            if 'N' in calcControls[hfx+'histType']:
                dat = G2el.getBLvalues(BLtables)
                refDict['FF']['FF'][iref] = dat.values()
            else:       #'X'
                dat = G2el.getFFvalues(FFtables,SQ)
                refDict['FF']['FF'][iref] = dat.values()
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = refDict['FF']['FF'][iref][Tindx]
        Uniq = np.inner(H,SGMT)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner(Uniq,(dXdata.T+Xdata.T))+Phi[:,np.newaxis])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        biso = -SQfactor*Uisodata
        Tiso = np.where(biso<1.,np.exp(biso),1.0)
        HbH = np.array([-np.inner(h,np.inner(bij,h)) for h in Uniq])
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = Tiso*Tuij*Mdata*Fdata/len(Uniq)
        fa = np.array([(FF+FP-Bab)*cosp*Tcorr,-FPP*sinp*Tcorr])
        fas = np.sum(np.sum(fa,axis=1),axis=1)        #real
        if not SGData['SGInv']:
            fb = np.array([(FF+FP-Bab)*sinp*Tcorr,FPP*cosp*Tcorr])
            fbs = np.sum(np.sum(fb,axis=1),axis=1)
        fasq = fas**2
        fbsq = fbs**2        #imaginary
        refl[9] = np.sum(fasq)+np.sum(fbsq)
        refl[10] = atan2d(fbs[0],fas[0])
    
def StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Compute structure factors for all h,k,l for phase
    puts the result, F^2, in each ref[8] in refList
    input:
    
    :param dict refDict: where
        'RefList' list where each ref = h,k,l,m,d,...
        'FF' dict of form factors - filed in below
    :param np.array G:      reciprocal metric tensor
    :param str pfx:    phase id string
    :param dict SGData: space group info. dictionary output from SpcGroup
    :param dict calcControls:
    :param dict ParmDict:

    '''        
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,calcControls,parmDict)
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType']:
        FP,FPP = G2el.BlenResCW(Tdata,BLtables,parmDict[hfx+'Lam'])
    elif 'X' in calcControls[hfx+'histType']:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    blkSize = 100       #no. of reflections in a block
    nRef = refDict['RefList'].shape[0]
    if not len(refDict['FF']):                #no form factors - 1st time thru StructureFactor
        if 'N' in calcControls[hfx+'histType']:
            dat = G2el.getBLvalues(BLtables)
            refDict['FF']['El'] = dat.keys()
            refDict['FF']['FF'] = np.ones((nRef,len(dat)))*dat.values()            
        else:       #'X'
            dat = G2el.getFFvalues(FFtables,0.)
            refDict['FF']['El'] = dat.keys()
            refDict['FF']['FF'] = np.ones((nRef,len(dat)))
            for iref,ref in enumerate(refDict['RefList']):
                SQ = 1./(2.*ref[4])**2
                dat = G2el.getFFvalues(FFtables,SQ)
                refDict['FF']['FF'][iref] *= dat.values()
#reflection processing begins here - big arrays!
    iBeg = 0            
    while iBeg < nRef:
        iFin = min(iBeg+blkSize,nRef)
        refl = refDict['RefList'][iBeg:iFin]
        H = refl.T[:3]
        SQ = 1./(2.*refl.T[4])**2
        SQfactor = 4.0*SQ*twopisq
        if 'T' in calcControls[hfx+'histType']:
            if 'P' in calcControls[hfx+'histType']:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[14])
            else:
                FP,FPP = G2el.BlenResTOF(Tdata,BLtables,refl.T[12])
            FP = np.repeat(FP.T,len(SGT),axis=0)
            FPP = np.repeat(FPP.T,len(SGT),axis=0)
        Bab = np.repeat(parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor),len(SGT))
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = np.repeat(refDict['FF']['FF'][iBeg:iFin].T[Tindx].T,len(SGT),axis=0)
        Uniq = np.reshape(np.inner(H.T,SGMT),(-1,3))
        Phi = np.inner(H.T,SGT).flatten()
        phase = twopi*(np.inner(Uniq,(dXdata+Xdata).T)+Phi[:,np.newaxis])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        biso = -SQfactor*Uisodata[:,np.newaxis]
        Tiso = np.repeat(np.where(biso<1.,np.exp(biso),1.0),len(SGT),axis=1).T
        HbH = -np.sum(Uniq.T*np.inner(bij,Uniq),axis=1)
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0).T
        Tcorr = Tiso*Tuij*Mdata*Fdata/len(SGMT)
        fa = np.array([((FF+FP).T-Bab).T*cosp*Tcorr,-FPP*sinp*Tcorr])
        fa = np.reshape(fa,(2,len(refl),len(SGT),len(Mdata)))
        fas = np.sum(np.sum(fa,axis=2),axis=2)        #real
        fbs = np.zeros_like(fas)
        if not SGData['SGInv']:
            fb = np.array([((FF+FP).T-Bab).T*sinp*Tcorr,FPP*cosp*Tcorr])
            fb = np.reshape(fb,(2,len(refl),len(SGT),len(Mdata)))
            fbs = np.sum(np.sum(fb,axis=2),axis=2)
        fasq = fas**2
        fbsq = fbs**2        #imaginary
        refl.T[9] = np.sum(fasq,axis=0)+np.sum(fbsq,axis=0)
        refl.T[10] = atan2d(fbs[0],fas[0])
        iBeg += blkSize
    
def StructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict):
    'Needs a doc string'
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
    SGT = np.array([ops[1] for ops in SGData['SGOps']])
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    nRef = len(refDict['RefList'])
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,calcControls,parmDict)
    mSize = len(Mdata)
    FF = np.zeros(len(Tdata))
    if 'NC' in calcControls[hfx+'histType']:
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
    for iref,refl in enumerate(refDict['RefList']):
        if 'T' in calcControls[hfx+'histType']:
            FP,FPP = G2el.BlenResCW(Tdata,BLtables,refl.T[12])
        H = np.array(refl[:3])
        SQ = 1./(2.*refl[4])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
        Bab = parmDict[phfx+'BabA']*dBabdA
        Tindx = np.array([refDict['FF']['El'].index(El) for El in Tdata])
        FF = refDict['FF']['FF'][iref].T[Tindx]
        Uniq = np.inner(H,SGMT)
        Phi = np.inner(H,SGT)
        phase = twopi*(np.inner((dXdata.T+Xdata.T),Uniq)+Phi[np.newaxis,:])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(Uniq)
        biso = -SQfactor*Uisodata
        Tiso = np.where(biso<1.,np.exp(biso),1.0)
        HbH = -np.inner(H,np.inner(bij,H))
        Hij = np.array([Mast*np.multiply.outer(U,U) for U in Uniq])
        Hij = np.array([G2lat.UijtoU6(Uij) for Uij in Hij])
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = Tiso*Tuij
        fot = (FF+FP-Bab)*occ*Tcorr
        fotp = FPP*occ*Tcorr
        fa = np.array([fot[:,np.newaxis]*cosp,fotp[:,np.newaxis]*cosp])       #non positions
        fb = np.array([fot[:,np.newaxis]*sinp,-fotp[:,np.newaxis]*sinp])
        
        fas = np.sum(np.sum(fa,axis=1),axis=1)
        fbs = np.sum(np.sum(fb,axis=1),axis=1)
        fax = np.array([-fot[:,np.newaxis]*sinp,-fotp[:,np.newaxis]*sinp])   #positions
        fbx = np.array([fot[:,np.newaxis]*cosp,-fot[:,np.newaxis]*cosp])
        #sum below is over Uniq
        dfadfr = np.sum(fa/occ[:,np.newaxis],axis=2)        #Fdata != 0 ever avoids /0. problem
        dfadx = np.sum(twopi*Uniq*fax[:,:,:,np.newaxis],axis=2)
        dfadui = np.sum(-SQfactor*fa,axis=2)
        dfadua = np.sum(-Hij*fa[:,:,:,np.newaxis],axis=2)
        dfadba = np.sum(-cosp*(occ*Tcorr)[:,np.newaxis],axis=1)
        #NB: the above have been checked against PA(1:10,1:2) in strfctr.for for al2O3!    
        dFdfr[iref] = 2.*(fas[0]*dfadfr[0]+fas[1]*dfadfr[1])*Mdata/len(Uniq)
        dFdx[iref] = 2.*(fas[0]*dfadx[0]+fas[1]*dfadx[1])
        dFdui[iref] = 2.*(fas[0]*dfadui[0]+fas[1]*dfadui[1])
        dFdua[iref] = 2.*(fas[0]*dfadua[0]+fas[1]*dfadua[1])
        dFdbab[iref] = 2.*fas[0]*np.array([np.sum(dfadba*dBabdA),np.sum(-dfadba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
        if not SGData['SGInv']:
            dfbdfr = np.sum(fb/occ[:,np.newaxis],axis=2)        #Fdata != 0 ever avoids /0. problem
            dfbdx = np.sum(twopi*Uniq*fbx[:,:,:,np.newaxis],axis=2)           
            dfbdui = np.sum(-SQfactor*fb,axis=2)
            dfbdua = np.sum(-Hij*fb[:,:,:,np.newaxis],axis=2)
            dfbdba = np.sum(-sinp*(occ*Tcorr)[:,np.newaxis],axis=1)
            dFdfr[iref] += 2.*(fbs[0]*dfbdfr[0]-fbs[1]*dfbdfr[1])*Mdata/len(Uniq)
            dFdx[iref] += 2.*(fbs[0]*dfbdx[0]+fbs[1]*dfbdx[1])
            dFdui[iref] += 2.*(fbs[0]*dfbdui[0]-fbs[1]*dfbdui[1])
            dFdua[iref] += 2.*(fbs[0]*dfbdua[0]+fbs[1]*dfbdua[1])
            dFdbab[iref] += 2.*fbs[0]*np.array([np.sum(dfbdba*dBabdA),np.sum(-dfbdba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
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
    dFdvDict[pfx+'BabA'] = dFdbab.T[0]
    dFdvDict[pfx+'BabU'] = dFdbab.T[1]
    return dFdvDict
    
def SCExtinction(ref,phfx,hfx,pfx,calcControls,parmDict,varyList):
    ''' Single crystal extinction function; returns extinction & derivative
    '''
    extCor = 1.0
    dervDict = {}
    if calcControls[phfx+'EType'] != 'None':
        SQ = 1/(4.*ref[4]**2)
        if 'C' in parmDict[hfx+'Type']:            
            cos2T = 1.0-2.*SQ*parmDict[hfx+'Lam']**2           #cos(2theta)
        else:   #'T'
            cos2T = 1.0-2.*SQ*ref[12]**2                       #cos(2theta)            
        if 'SXC' in parmDict[hfx+'Type']:
            AV = 7.9406e5/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            P12 = (calcControls[phfx+'Cos2TM']+cos2T**4)/(calcControls[phfx+'Cos2TM']+cos2T**2)
            PLZ = AV*P12*ref[7]*parmDict[hfx+'Lam']**2
        elif 'SNT' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = SQ
            PLZ = AV*ref[7]*ref[12]**2
        elif 'SNC' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            PLZ = AV*ref[9]*parmDict[hfx+'Lam']**2      #Fcsq as per GSAS, why not FcTsq (ref[9])?
            
        if 'Primary' in calcControls[phfx+'EType']:
            PLZ *= 1.5
        else:
            if 'C' in parmDict[hfx+'Type']:
                PLZ *= calcControls[phfx+'Tbar']
            else: #'T'
                PLZ *= ref[13]      #t-bar
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

        if 'Primary' in calcControls[phfx+'EType'] and phfx+'Ep' in varyList:
            dervDict[phfx+'Ep'] = -ref[7]*PLZ*PF3
        if 'II' in calcControls[phfx+'EType'] and phfx+'Es' in varyList:
            dervDict[phfx+'Es'] = -ref[7]*PLZ*PF3*(PSIG/parmDict[phfx+'Es'])**3
        if 'I' in calcControls[phfx+'EType'] and phfx+'Eg' in varyList:
            dervDict[phfx+'Eg'] = -ref[7]*PLZ*PF3*(PSIG/parmDict[phfx+'Eg'])**3*PL**2
               
    return 1./extCor,dervDict
    
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def GetNewCellParms(parmDict,varyList):
    'Needs a doc string'
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
    
def SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    'Spherical harmonics texture'
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5]
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Omega'],parmDict[hfx+'Chi'],parmDict[hfx+'Phi'],parmDict[hfx+'Azimuth']]
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
    
def SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    'Spherical harmonics texture derivatives'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5]
    FORPI = 4.0*np.pi
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    odfCor = 1.0
    dFdODF = {}
    dFdSA = [0,0,0]
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Omega'],parmDict[hfx+'Chi'],parmDict[hfx+'Phi'],parmDict[hfx+'Azimuth']]
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
    
def SHPOcal(refl,g,phfx,hfx,SGData,calcControls,parmDict):
    'spherical harmonics preferred orientation (cylindrical symmetry only)'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5]
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangl = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [0.,0.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(tth/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta)
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
    return np.squeeze(odfCor)
    
def SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict):
    'spherical harmonics preferred orientation derivatives (cylindrical symmetry only)'
    if 'T' in calcControls[hfx+'histType']:
        tth = parmDict[hfx+'2-theta']
    else:
        tth = refl[5]
    FORPI = 12.5663706143592
    odfCor = 1.0
    dFdODF = {}
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangl = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [0.,0.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(tth/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta) 
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
        dFdODF[phfx+item] = Kcsl*Lnorm
    return odfCor,dFdODF
    
def GetPrefOri(refl,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict):
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
    
def GetPrefOriDerv(refl,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict):
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
            POcorr,POderv = SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr,POderv
    
def GetAbsorb(refl,hfx,calcControls,parmDict):
    'Needs a doc string'
    if 'Debye' in calcControls[hfx+'instType']:
        if 'T' in calcControls[hfx+'histType']:
            return G2pwd.Absorb('Cylinder',parmDict[hfx+'Absorption']*refl[14],parmDict[hfx+'2-theta'],0,0)
        else:
            return G2pwd.Absorb('Cylinder',parmDict[hfx+'Absorption'],refl[5],0,0)
    else:
        return G2pwd.SurfaceRough(parmDict[hfx+'SurfRoughA'],parmDict[hfx+'SurfRoughB'],refl[5])
    
def GetAbsorbDerv(refl,hfx,calcControls,parmDict):
    'Needs a doc string'
    if 'Debye' in calcControls[hfx+'instType']:
        if 'T' in calcControls[hfx+'histType']:
            return G2pwd.AbsorbDerv('Cylinder',parmDict[hfx+'Absorption']*refl[14],parmDict[hfx+'2-theta'],0,0)
        else:
            return G2pwd.AbsorbDerv('Cylinder',parmDict[hfx+'Absorption'],refl[5],0,0)
    else:
        return np.array(G2pwd.SurfaceRoughDerv(parmDict[hfx+'SurfRoughA'],parmDict[hfx+'SurfRoughB'],refl[5]))
        
def GetPwdrExt(refl,pfx,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    coef = np.array([-0.5,0.25,-0.10416667,0.036458333,-0.0109375,2.8497409E-3])
    pi2 = np.sqrt(2./np.pi)
    if 'T' in calcControls[hfx+'histType']:
        sth2 = sind(parmDict[hfx+'2-theta']/2.)**2
        wave = refl[14]
    else:   #'C'W
        sth2 = sind(refl[5]/2.)**2
        wave = parmDict.get(hfx+'Lam',parmDict.get(hfx+'Lam1',1.0))
    c2th = 1.-2.0*sth2
    flv2 = refl[9]*(wave/parmDict[pfx+'Vol'])**2
    if 'X' in calcControls[hfx+'histType']:
        flv2 *= 0.079411*(1.0+c2th**2)/2.0
    xfac = flv2*parmDict[phfx+'Extinction']
    exb = 1.0
    if xfac > -1.:
        exb = 1./(1.+xfac)
    exl = 1.0
    if 0 < xfac <= 1.:
        xn = np.array([xfac**(i+1) for i in range(6)])
        exl = np.sum(xn*coef)
    elif xfac > 1.:
        xfac2 = 1./np.sqrt(xfac)
        exl = pi2*(1.-0.125/xfac)*xfac2
    return exb*sth2+exl*(1.-sth2)
    
def GetPwdrExtDerv(refl,pfx,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    delt = 0.001
    parmDict[phfx+'Extinction'] += delt
    plus = GetPwdrExt(refl,pfx,phfx,hfx,calcControls,parmDict)
    parmDict[phfx+'Extinction'] -= 2.*delt
    minus = GetPwdrExt(refl,pfx,phfx,hfx,calcControls,parmDict)
    parmDict[phfx+'Extinction'] += delt
    return (plus-minus)/(2.*delt)    
    
def GetIntensityCorr(refl,uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    'Needs a doc string'    #need powder extinction!
    Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3]               #scale*multiplicity
    if 'X' in parmDict[hfx+'Type']:
        Icorr *= G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])[0]
    POcorr = 1.0
    if pfx+'SHorder' in parmDict:                 #generalized spherical harmonics texture
        POcorr = SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict)
    elif calcControls[phfx+'poType'] == 'MD':         #March-Dollase
        POcorr = GetPrefOri(refl,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict)
    elif calcControls[phfx+'SHord']:                #cylindrical spherical harmonics
        POcorr = SHPOcal(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    Icorr *= POcorr
    AbsCorr = 1.0
    AbsCorr = GetAbsorb(refl,hfx,calcControls,parmDict)
    Icorr *= AbsCorr
    ExtCorr = GetPwdrExt(refl,pfx,phfx,hfx,calcControls,parmDict)
    Icorr *= ExtCorr
    return Icorr,POcorr,AbsCorr,ExtCorr
    
def GetIntensityDerv(refl,wave,uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    'Needs a doc string'    #need powder extinction derivs!
    dIdsh = 1./parmDict[hfx+'Scale']
    dIdsp = 1./parmDict[phfx+'Scale']
    if 'X' in parmDict[hfx+'Type']:
        pola,dIdPola = G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])
        dIdPola /= pola
    else:       #'N'
        dIdPola = 0.0
    dFdODF = {}
    dFdSA = [0,0,0]
    dIdPO = {}
    if pfx+'SHorder' in parmDict:
        odfCor,dFdODF,dFdSA = SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict)
        for iSH in dFdODF:
            dFdODF[iSH] /= odfCor
        for i in range(3):
            dFdSA[i] /= odfCor
    elif calcControls[phfx+'poType'] == 'MD' or calcControls[phfx+'SHord']:
        POcorr,dIdPO = GetPrefOriDerv(refl,uniq,G,g,phfx,hfx,SGData,calcControls,parmDict)        
        for iPO in dIdPO:
            dIdPO[iPO] /= POcorr
    if 'T' in parmDict[hfx+'Type']:
        dFdAb = GetAbsorbDerv(refl,hfx,calcControls,parmDict)*wave/refl[16] #wave/abs corr
        dFdEx = GetPwdrExtDerv(refl,pfx,phfx,hfx,calcControls,parmDict)/refl[17]    #/ext corr
    else:
        dFdAb = GetAbsorbDerv(refl,hfx,calcControls,parmDict)*wave/refl[13] #wave/abs corr
        dFdEx = GetPwdrExtDerv(refl,pfx,phfx,hfx,calcControls,parmDict)/refl[14]    #/ext corr        
    return dIdsh,dIdsp,dIdPola,dIdPO,dFdODF,dFdSA,dFdAb,dFdEx
        
def GetSampleSigGam(refl,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    if 'C' in calcControls[hfx+'histType']:     #All checked & OK
        costh = cosd(refl[5]/2.)
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
            Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR = G2pwd.ellipseSize(H,Sij,GB)
            Sgam = 1.8*wave/(np.pi*costh*lenR)
        #microstrain                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 0.018*parmDict[phfx+'Mustrain;i']*tand(refl[5]/2.)/np.pi
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            Mgam = 0.018*Si*Sa*tand(refl[5]/2.)/(np.pi*np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
        else:       #generalized - P.W. Stephens model
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain:'+str(i)]*strm
            Mgam = 0.018*refl[4]**2*tand(refl[5]/2.)*np.sqrt(Sum)/np.pi
    elif 'T' in calcControls[hfx+'histType']:       #All checked & OK
        #crystallite size
        if calcControls[phfx+'SizeType'] == 'isotropic':    #OK
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4]**2/parmDict[phfx+'Size;i']
        elif calcControls[phfx+'SizeType'] == 'uniaxial':   #OK
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4]**2/(parmDict[phfx+'Size;i']*parmDict[phfx+'Size;a'])
            Sgam *= np.sqrt((sinP*parmDict[phfx+'Size;a'])**2+(cosP*parmDict[phfx+'Size;i'])**2)
        else:           #ellipsoidal crystallites   #OK
            Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR = G2pwd.ellipseSize(H,Sij,GB)
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4]**2/lenR
        #microstrain                
        if calcControls[phfx+'MustrainType'] == 'isotropic':    #OK
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4]*parmDict[phfx+'Mustrain;i']
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':   #OK
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4]*Si*Sa/np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
        else:       #generalized - P.W. Stephens model  OK
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain:'+str(i)]*strm
            Mgam = 1.e-6*parmDict[hfx+'difC']*np.sqrt(Sum)*refl[4]**3
            
    gam = Sgam*parmDict[phfx+'Size;mx']+Mgam*parmDict[phfx+'Mustrain;mx']
    sig = (Sgam*(1.-parmDict[phfx+'Size;mx']))**2+(Mgam*(1.-parmDict[phfx+'Mustrain;mx']))**2
    sig /= ateln2
    return sig,gam
        
def GetSampleSigGamDerv(refl,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict):
    'Needs a doc string'
    gamDict = {}
    sigDict = {}
    if 'C' in calcControls[hfx+'histType']:         #All checked & OK
        costh = cosd(refl[5]/2.)
        tanth = tand(refl[5]/2.)
        #crystallite size derivatives
        if calcControls[phfx+'SizeType'] == 'isotropic':
            Sgam = 1.8*wave/(np.pi*parmDict[phfx+'Size;i']*costh)
            gamDict[phfx+'Size;i'] = -1.8*wave*parmDict[phfx+'Size;mx']/(np.pi*costh)
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
            Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
            Sgam = const/lenR
            for i,item in enumerate([phfx+'Size:%d'%(j) for j in range(6)]):
                gamDict[item] = -(const/lenR**2)*dRdS[i]*parmDict[phfx+'Size;mx']
                sigDict[item] = -2.*Sgam*(const/lenR**2)*dRdS[i]*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        gamDict[phfx+'Size;mx'] = Sgam
        sigDict[phfx+'Size;mx'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])/ateln2
                
        #microstrain derivatives                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 0.018*parmDict[phfx+'Mustrain;i']*tand(refl[5]/2.)/np.pi
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
            const = 0.018*refl[4]**2*tanth/np.pi
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain:'+str(i)]*strm
                gamDict[phfx+'Mustrain:'+str(i)] = strm*parmDict[phfx+'Mustrain;mx']/2.
                sigDict[phfx+'Mustrain:'+str(i)] = strm*(1.-parmDict[phfx+'Mustrain;mx'])**2
            Mgam = const*np.sqrt(Sum)
            for i in range(len(Strms)):
                gamDict[phfx+'Mustrain:'+str(i)] *= Mgam/Sum
                sigDict[phfx+'Mustrain:'+str(i)] *= const**2/ateln2
        gamDict[phfx+'Mustrain;mx'] = Mgam
        sigDict[phfx+'Mustrain;mx'] = -2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])/ateln2
    else:   #'T'OF - All checked & OK
        if calcControls[phfx+'SizeType'] == 'isotropic':    #OK
            Sgam = 1.e-4*parmDict[hfx+'difC']*refl[4]**2/parmDict[phfx+'Size;i']
            gamDict[phfx+'Size;i'] = -Sgam*parmDict[phfx+'Size;mx']/parmDict[phfx+'Size;i']
            sigDict[phfx+'Size;i'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])**2/(ateln2*parmDict[phfx+'Size;i'])
        elif calcControls[phfx+'SizeType'] == 'uniaxial':   #OK
            const = 1.e-4*parmDict[hfx+'difC']*refl[4]**2
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
            const = 1.e-4*parmDict[hfx+'difC']*refl[4]**2
            Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
            H = np.array(refl[:3])
            lenR,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
            Sgam = const/lenR
            for i,item in enumerate([phfx+'Size:%d'%(j) for j in range(6)]):
                gamDict[item] = -(const/lenR**2)*dRdS[i]*parmDict[phfx+'Size;mx']
                sigDict[item] = -2.*Sgam*(const/lenR**2)*dRdS[i]*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        gamDict[phfx+'Size;mx'] = Sgam  #OK
        sigDict[phfx+'Size;mx'] = -2.*Sgam**2*(1.-parmDict[phfx+'Size;mx'])/ateln2  #OK
                
        #microstrain derivatives                
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            Mgam = 1.e-6*parmDict[hfx+'difC']*refl[4]*parmDict[phfx+'Mustrain;i']
            gamDict[phfx+'Mustrain;i'] =  1.e-6*refl[4]*parmDict[hfx+'difC']*parmDict[phfx+'Mustrain;mx']   #OK
            sigDict[phfx+'Mustrain;i'] =  2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])**2/(ateln2*parmDict[phfx+'Mustrain;i'])        
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain;i']
            Sa = parmDict[phfx+'Mustrain;a']
            gami = 1.e-6*parmDict[hfx+'difC']*refl[4]*Si*Sa
            sqtrm = np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
            Mgam = gami/sqtrm
            dsi = -gami*Si*cosP**2/sqtrm**3
            dsa = -gami*Sa*sinP**2/sqtrm**3
            gamDict[phfx+'Mustrain;i'] = (Mgam/Si+dsi)*parmDict[phfx+'Mustrain;mx']
            gamDict[phfx+'Mustrain;a'] = (Mgam/Sa+dsa)*parmDict[phfx+'Mustrain;mx']
            sigDict[phfx+'Mustrain;i'] = 2*(Mgam/Si+dsi)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2
            sigDict[phfx+'Mustrain;a'] = 2*(Mgam/Sa+dsa)*Mgam*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2       
        else:       #generalized - P.W. Stephens model OK
            pwrs = calcControls[phfx+'MuPwrs']
            Strms = G2spc.MustrainCoeff(refl[:3],SGData)
            const = 1.e-6*parmDict[hfx+'difC']*refl[4]**3
            Sum = 0
            for i,strm in enumerate(Strms):
                Sum += parmDict[phfx+'Mustrain:'+str(i)]*strm
                gamDict[phfx+'Mustrain:'+str(i)] = strm*parmDict[phfx+'Mustrain;mx']/2.
                sigDict[phfx+'Mustrain:'+str(i)] = strm*(1.-parmDict[phfx+'Mustrain;mx'])**2
            Mgam = const*np.sqrt(Sum)
            for i in range(len(Strms)):
                gamDict[phfx+'Mustrain:'+str(i)] *= Mgam/Sum
                sigDict[phfx+'Mustrain:'+str(i)] *= const**2/ateln2        
        gamDict[phfx+'Mustrain;mx'] = Mgam
        sigDict[phfx+'Mustrain;mx'] = -2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])/ateln2
        
    return sigDict,gamDict
        
def GetReflPos(refl,wave,A,hfx,calcControls,parmDict):
    'Needs a doc string'
    h,k,l = refl[:3]
    d = 1./np.sqrt(G2lat.calc_rDsq(np.array([h,k,l]),A))

    refl[4] = d
    if 'C' in calcControls[hfx+'histType']:
        pos = 2.0*asind(wave/(2.0*d))+parmDict[hfx+'Zero']
        const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            pos -= const*(4.*parmDict[hfx+'Shift']*cosd(pos/2.0)+ \
                parmDict[hfx+'Transparency']*sind(pos)*100.0)            #trans(=1/mueff) in cm
        else:               #Debye-Scherrer - simple but maybe not right
            pos -= const*(parmDict[hfx+'DisplaceX']*cosd(pos)+parmDict[hfx+'DisplaceY']*sind(pos))
    elif 'T' in calcControls[hfx+'histType']:
        pos = parmDict[hfx+'difC']*d+parmDict[hfx+'difA']*d**2+parmDict[hfx+'difB']/d+parmDict[hfx+'Zero']
        #do I need sample position effects - maybe?
    return pos

def GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict):
    'Needs a doc string'
    dpr = 180./np.pi
    h,k,l = refl[:3]
    dstsq = G2lat.calc_rDsq(np.array([h,k,l]),A)
    dst = np.sqrt(dstsq)
    dsp = 1./dst
    if 'C' in calcControls[hfx+'histType']:
        pos = refl[5]-parmDict[hfx+'Zero']
        const = dpr/np.sqrt(1.0-wave**2*dstsq/4.0)
        dpdw = const*dst
        dpdA = np.array([h**2,k**2,l**2,h*k,h*l,k*l])*const*wave/(2.0*dst)
        dpdZ = 1.0
        const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            dpdSh = -4.*const*cosd(pos/2.0)
            dpdTr = -const*sind(pos)*100.0
            return dpdA,dpdw,dpdZ,dpdSh,dpdTr,0.,0.
        else:               #Debye-Scherrer - simple but maybe not right
            dpdXd = -const*cosd(pos)
            dpdYd = -const*sind(pos)
            return dpdA,dpdw,dpdZ,0.,0.,dpdXd,dpdYd
    elif 'T' in calcControls[hfx+'histType']:
        dpdA = -np.array([h**2,k**2,l**2,h*k,h*l,k*l])*parmDict[hfx+'difC']*dsp**3/2.
        dpdZ = 1.0
        dpdDC = dsp
        dpdDA = dsp**2
        dpdDB = 1./dsp
        return dpdA,dpdZ,dpdDC,dpdDA,dpdDB
            
def GetHStrainShift(refl,SGData,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+l**2)+ \
            refl[4]**2*parmDict[phfx+'eA']*((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2
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
    if 'C' in calcControls[hfx+'histType']:
        return -180.*Dij*refl[4]**2*tand(refl[5]/2.0)/np.pi
    else:
        return -Dij*parmDict[hfx+'difC']*refl[4]**2/2.
            
def GetHStrainShiftDerv(refl,SGData,phfx,hfx,calcControls,parmDict):
    'Needs a doc string'
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        dDijDict = {phfx+'D11':h**2+k**2+l**2,
            phfx+'eA':refl[4]**2*((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2}
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
            names.append()
    else:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2,
            phfx+'D12':h*k,phfx+'D13':h*l,phfx+'D23':k*l}
    if 'C' in calcControls[hfx+'histType']:
        for item in dDijDict:
            dDijDict[item] *= 180.0*refl[4]**2*tand(refl[5]/2.0)/np.pi
    else:
        for item in dDijDict:
            dDijDict[item] *= -parmDict[hfx+'difC']*refl[4]**3/2.
    return dDijDict
    
def GetDij(phfx,SGData,parmDict):
    HSvals = [parmDict[phfx+name] for name in G2spc.HStrainNames(SGData)]
    return G2spc.HStrainVals(HSvals,SGData)
                
def GetFobsSq(Histograms,Phases,parmDict,calcControls):
    'Needs a doc string'
    histoList = Histograms.keys()
    histoList.sort()
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
                if hfx+'Lam1' in parmDict.keys():
                    Ka2 = True
                    lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
                    kRatio = parmDict[hfx+'I(L2)/I(L1)']
            x,y,w,yc,yb,yd = Histogram['Data']
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            ymb = np.array(y-yb)
            ymb = np.where(ymb,ymb,1.0)
            ycmb = np.array(yc-yb)
            ratio = 1./np.where(ycmb,ycmb/ymb,1.e10)          
            refLists = Histogram['Reflection Lists']
            for phase in refLists:
                Phase = Phases[phase]
                pId = Phase['pId']
                phfx = '%d:%d:'%(pId,hId)
                refDict = refLists[phase]
                sumFo = 0.0
                sumdF = 0.0
                sumFosq = 0.0
                sumdFsq = 0.0
                for refl in refDict['RefList']:
                    if 'C' in calcControls[hfx+'histType']:
                        yp = np.zeros_like(yb)
                        Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5],refl[6],refl[7],shl)
                        iBeg = max(xB,np.searchsorted(x,refl[5]-fmin))
                        iFin = max(xB,min(np.searchsorted(x,refl[5]+fmax),xF))
                        iFin2 = iFin
                        if not iBeg+iFin:       #peak below low limit - skip peak
                            continue
                        elif not iBeg-iFin:     #peak above high limit - done
                            break
                        elif iBeg < iFin:
                            yp[iBeg:iFin] = refl[11]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,ma.getdata(x[iBeg:iFin]))    #>90% of time spent here
                            if Ka2:
                                pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                                Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6],refl[7],shl)
                                iBeg2 = max(xB,np.searchsorted(x,pos2-fmin))
                                iFin2 = min(np.searchsorted(x,pos2+fmax),xF)
                                yp[iBeg2:iFin2] += refl[11]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,ma.getdata(x[iBeg2:iFin2]))        #and here
                            refl[8] = np.sum(np.where(ratio[iBeg:iFin2]>0.,yp[iBeg:iFin2]*ratio[iBeg:iFin2]/(refl[11]*(1.+kRatio)),0.0))
                    elif 'T' in calcControls[hfx+'histType']:
                        yp = np.zeros_like(yb)
                        Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5],refl[12],refl[13],refl[6],refl[7])
                        iBeg = max(xB,np.searchsorted(x,refl[5]-fmin))
                        iFin = max(xB,min(np.searchsorted(x,refl[5]+fmax),xF))
                        if iBeg < iFin:
                            yp[iBeg:iFin] = refl[11]*refl[9]*G2pwd.getEpsVoigt(refl[5],refl[12],refl[13],refl[6],refl[7],ma.getdata(x[iBeg:iFin]))  #>90% of time spent here
                            refl[8] = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin]/refl[11],0.0))
                    Fo = np.sqrt(np.abs(refl[8]))
                    Fc = np.sqrt(np.abs(refl[9]))
                    sumFo += Fo
                    sumFosq += refl[8]**2
                    sumdF += np.abs(Fo-Fc)
                    sumdFsq += (refl[8]-refl[9])**2
                Histogram['Residuals'][phfx+'Rf'] = min(100.,(sumdF/sumFo)*100.)
                Histogram['Residuals'][phfx+'Rf^2'] = min(100.,np.sqrt(sumdFsq/sumFosq)*100.)
                Histogram['Residuals'][phfx+'Nref'] = len(refDict['RefList'])
                Histogram['Residuals']['hId'] = hId
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            Histogram['Residuals']['hId'] = Histograms[histogram]['hId']
                
def getPowderProfile(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    'Needs a doc string'
    
    def GetReflSigGamCW(refl,wave,G,GB,phfx,calcControls,parmDict):
        U = parmDict[hfx+'U']
        V = parmDict[hfx+'V']
        W = parmDict[hfx+'W']
        X = parmDict[hfx+'X']
        Y = parmDict[hfx+'Y']
        tanPos = tand(refl[5]/2.0)
        Ssig,Sgam = GetSampleSigGam(refl,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict)
        sig = U*tanPos**2+V*tanPos+W+Ssig     #save peak sigma
        sig = max(0.001,sig)
        gam = X/cosd(refl[5]/2.0)+Y*tanPos+Sgam     #save peak gamma
        gam = max(0.001,gam)
        return sig,gam
                
    def GetReflSigGamTOF(refl,G,GB,phfx,calcControls,parmDict):
        sig = parmDict[hfx+'sig-0']+parmDict[hfx+'sig-1']*refl[4]**2+   \
            parmDict[hfx+'sig-2']*refl[4]**4+parmDict[hfx+'sig-q']/refl[4]**2
        gam = parmDict[hfx+'X']*refl[4]+parmDict[hfx+'Y']*refl[4]**2
        Ssig,Sgam = GetSampleSigGam(refl,0.0,G,GB,SGData,hfx,phfx,calcControls,parmDict)
        sig += Ssig
        gam += Sgam
        return sig,gam
        
    def GetReflAlpBet(refl,hfx,parmDict):
        alp = parmDict[hfx+'alpha']/refl[4]
        bet = parmDict[hfx+'beta-0']+parmDict[hfx+'beta-1']/refl[4]**4+parmDict[hfx+'beta-q']/refl[4]**2
        return alp,bet
        
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    yb = G2pwd.getBackground(hfx,parmDict,bakType,calcControls[hfx+'histType'],x)
    yc = np.zeros_like(yb)
    cw = np.diff(x)
    cw = np.append(cw,cw[-1])
        
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.002)
        Ka2 = False
        if hfx+'Lam1' in parmDict.keys():
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    for phase in Histogram['Reflection Lists']:
        refDict = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        hfx = ':%d:'%(hId)
        SGData = Phase['General']['SGData']
        SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
        Dij = GetDij(phfx,SGData,parmDict)
        A = [parmDict[pfx+'A%d'%(i)]+Dij[i] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        Vst = np.sqrt(nl.det(G))    #V*
        if not Phase['General'].get('doPawley'):
            time0 = time.time()
            StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
#            print 'sf calc time: %.3fs'%(time.time()-time0)
        time0 = time.time()
        badPeak = False
        for iref,refl in enumerate(refDict['RefList']):
            if 'C' in calcControls[hfx+'histType']:
                h,k,l = refl[:3]
                Uniq = np.inner(refl[:3],SGMT)
                refl[5] = GetReflPos(refl,wave,A,hfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = 1./(2.*sind(refl[5]/2.)**2*cosd(refl[5]/2.))           #Lorentz correction
#                refl[5] += GetHStrainShift(refl,SGData,phfx,hfx,calcControls,parmDict)               #apply hydrostatic strain shift
                refl[6:8] = GetReflSigGamCW(refl,wave,G,GB,phfx,calcControls,parmDict)    #peak sig & gam
                refl[11:15] = GetIntensityCorr(refl,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                refl[11] *= Vst*Lorenz
                 
                if Phase['General'].get('doPawley'):
                    try:
                        pInd =pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        refl[9] = parmDict[pInd]
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                elif iBeg > iFin:   #bad peak coeff - skip
                    badPeak = True
                    continue
                yc[iBeg:iFin] += refl[11]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,ma.getdata(x[iBeg:iFin]))    #>90% of time spent here
                if Ka2:
                    pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6],refl[7],shl)
                    iBeg = np.searchsorted(x,pos2-fmin)
                    iFin = np.searchsorted(x,pos2+fmax)
                    if not iBeg+iFin:       #peak below low limit - skip peak
                        continue
                    elif not iBeg-iFin:     #peak above high limit - done
                        return yc,yb
                    elif iBeg > iFin:   #bad peak coeff - skip
                        continue
                    yc[iBeg:iFin] += refl[11]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,ma.getdata(x[iBeg:iFin]))        #and here
            elif 'T' in calcControls[hfx+'histType']:
                h,k,l = refl[:3]
                Uniq = np.inner(refl[:3],SGMT)
                refl[5] = GetReflPos(refl,0.0,A,hfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = sind(parmDict[hfx+'2-theta']/2)*refl[4]**4                                                #TOF Lorentz correction
#                refl[5] += GetHStrainShift(refl,SGData,phfx,hfx,calcControls,parmDict)               #apply hydrostatic strain shift
                refl[6:8] = GetReflSigGamTOF(refl,G,GB,phfx,calcControls,parmDict)    #peak sig & gam
                refl[12:14] = GetReflAlpBet(refl,hfx,parmDict)
                refl[11],refl[15],refl[16],refl[17] = GetIntensityCorr(refl,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                refl[11] *= Vst*Lorenz
                if Phase['General'].get('doPawley'):
                    try:
                        pInd =pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        refl[9] = parmDict[pInd]
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5],refl[12],refl[13],refl[6],refl[7])
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                elif iBeg > iFin:   #bad peak coeff - skip
                    badPeak = True
                    continue
                yc[iBeg:iFin] += refl[11]*refl[9]*G2pwd.getEpsVoigt(refl[5],refl[12],refl[13],refl[6],refl[7],ma.getdata(x[iBeg:iFin]))/cw[iBeg:iFin]
#        print 'profile calc time: %.3fs'%(time.time()-time0)
    if badPeak:
        print 'ouch #4 bad profile coefficients yield negative peak width; some reflections skipped' 
    return yc,yb
    
def getPowderProfileDerv(parmDict,x,varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup):
    'Needs a doc string'
    
    def cellVaryDerv(pfx,SGData,dpdA): 
        if SGData['SGLaue'] in ['-1',]:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],
                [pfx+'A3',dpdA[3]],[pfx+'A4',dpdA[4]],[pfx+'A5',dpdA[5]]]
        elif SGData['SGLaue'] in ['2/m',]:
            if SGData['SGUniq'] == 'a':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A3',dpdA[3]]]
            elif SGData['SGUniq'] == 'b':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A4',dpdA[4]]]
            else:
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A5',dpdA[5]]]
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
    dependentVars = G2mv.GetDependentVars()
    depDerivDict = {}
    for j in dependentVars:
        depDerivDict[j] = np.zeros(shape=(len(x)))
    #print 'dependent vars',dependentVars
    lenX = len(x)                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    dMdv = np.zeros(shape=(len(varylist),len(x)))
    dMdb,dMddb,dMdpk = G2pwd.getBackgroundDerv(hfx,parmDict,bakType,calcControls[hfx+'histType'],x)
    if hfx+'Back;0' in varylist: # for now assume that Back;x vars to not appear in constraints
        bBpos =varylist.index(hfx+'Back;0')
        dMdv[bBpos:bBpos+len(dMdb)] = dMdb
    names = [hfx+'DebyeA',hfx+'DebyeR',hfx+'DebyeU']
    for name in varylist:
        if 'Debye' in name:
            id = int(name.split(';')[-1])
            parm = name[:int(name.rindex(';'))]
            ip = names.index(parm)
            dMdv[varylist.index(name)] = dMddb[3*id+ip]
    names = [hfx+'BkPkpos',hfx+'BkPkint',hfx+'BkPksig',hfx+'BkPkgam']
    for name in varylist:
        if 'BkPk' in name:
            parm,id = name.split(';')
            id = int(id)
            if parm in names:
                ip = names.index(parm)
                dMdv[varylist.index(name)] = dMdpk[4*id+ip]
    cw = np.diff(x)
    cw = np.append(cw,cw[-1])
    Ka2 = False #also for TOF!
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.002)
        if hfx+'Lam1' in parmDict.keys():
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    for phase in Histogram['Reflection Lists']:
        refDict = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        SGData = Phase['General']['SGData']
        SGMT = np.array([ops[0].T for ops in SGData['SGOps']])
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        Dij = GetDij(phfx,SGData,parmDict)
        A = [parmDict[pfx+'A%d'%(i)]+Dij[i] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        if not Phase['General'].get('doPawley'):
            time0 = time.time()
            dFdvDict = StructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
#            print 'sf-derv time %.3fs'%(time.time()-time0)
            ApplyRBModelDervs(dFdvDict,parmDict,rigidbodyDict,Phase)
        time0 = time.time()
        for iref,refl in enumerate(refDict['RefList']):
            h,k,l = refl[:3]
            Uniq = np.inner(refl[:3],SGMT)
            if 'T' in calcControls[hfx+'histType']:
                wave = refl[14]
            dIdsh,dIdsp,dIdpola,dIdPO,dFdODF,dFdSA,dFdAb,dFdEx = GetIntensityDerv(refl,wave,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
            if 'C' in calcControls[hfx+'histType']:        #CW powder
                Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5],refl[6],refl[7],shl)
            else: #'T'OF
                Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5],refl[12],refl[13],refl[6],refl[7])
            iBeg = np.searchsorted(x,refl[5]-fmin)
            iFin = np.searchsorted(x,refl[5]+fmax)
            if not iBeg+iFin:       #peak below low limit - skip peak
                continue
            elif not iBeg-iFin:     #peak above high limit - done
                break
            pos = refl[5]
            if 'C' in calcControls[hfx+'histType']:
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                lenBF = iFin-iBeg
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdFCJVoigt3(refl[5],refl[6],refl[7],shl,ma.getdata(x[iBeg:iFin]))
                for i in range(5):
                    dMdpk[i] += 100.*cw[iBeg:iFin]*refl[11]*refl[9]*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':np.zeros_like(dMdpk[0])}
                if Ka2:
                    pos2 = refl[5]+lamRatio*tanth       # + 360/pi * Dlam/lam * tan(th)
                    iBeg2 = np.searchsorted(x,pos2-fmin)
                    iFin2 = np.searchsorted(x,pos2+fmax)
                    if iBeg2-iFin2:
                        lenBF2 = iFin2-iBeg2
                        dMdpk2 = np.zeros(shape=(6,lenBF2))
                        dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6],refl[7],shl,ma.getdata(x[iBeg2:iFin2]))
                        for i in range(5):
                            dMdpk2[i] = 100.*cw[iBeg2:iFin2]*refl[11]*refl[9]*kRatio*dMdipk2[i]
                        dMdpk2[5] = 100.*cw[iBeg2:iFin2]*refl[11]*dMdipk2[0]
                        dervDict2 = {'int':dMdpk2[0],'pos':dMdpk2[1],'sig':dMdpk2[2],'gam':dMdpk2[3],'shl':dMdpk2[4],'L1/L2':dMdpk2[5]*refl[9]}
            else:   #'T'OF
                lenBF = iFin-iBeg
                if lenBF < 0:   #bad peak coeff
                    break
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdEpsVoigt(refl[5],refl[12],refl[13],refl[6],refl[7],ma.getdata(x[iBeg:iFin]))
                for i in range(6):
                    dMdpk[i] += refl[11]*refl[9]*dMdipk[i]      #cw[iBeg:iFin]*
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5]}            
            if Phase['General'].get('doPawley'):
                dMdpw = np.zeros(len(x))
                try:
                    pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                    idx = varylist.index(pIdx)
                    dMdpw[iBeg:iFin] = dervDict['int']/refl[9]
                    if Ka2: #not for TOF either
                        dMdpw[iBeg2:iFin2] += dervDict2['int']/refl[9]
                    dMdv[idx] = dMdpw
                except: # ValueError:
                    pass
            if 'C' in calcControls[hfx+'histType']:
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY = GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],
                    hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'SH/L':[1.0,'shl'],
                    hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
                    hfx+'DisplaceY':[dpdY,'pos'],}
                if 'Bragg' in calcControls[hfx+'instType']:
                    names.update({hfx+'SurfRoughA':[dFdAb[0],'int'],
                        hfx+'SurfRoughB':[dFdAb[1],'int'],})
                else:
                    names.update({hfx+'Absorption':[dFdAb,'int'],})
            else:   #'T'OF
                dpdA,dpdZ,dpdDC,dpdDA,dpdDB = GetReflPosDerv(refl,0.0,A,hfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'difC':[dpdDC,'pos'],hfx+'difA':[dpdDA,'pos'],hfx+'difB':[dpdDB,'pos'],
                    hfx+'Zero':[dpdZ,'pos'],hfx+'X':[refl[4],'gam'],hfx+'Y':[refl[4]**2,'gam'],
                    hfx+'alpha':[1./refl[4],'alp'],hfx+'beta-0':[1.0,'bet'],hfx+'beta-1':[1./refl[4]**4,'bet'],
                    hfx+'beta-q':[1./refl[4]**2,'bet'],hfx+'sig-0':[1.0,'sig'],hfx+'sig-1':[refl[4]**2,'sig'],
                    hfx+'sig-2':[refl[4]**4,'sig'],hfx+'sig-q':[1./refl[4]**2,'sig'],
                    hfx+'Absorption':[dFdAb,'int'],phfx+'Extinction':[dFdEx,'int'],}
            for name in names:
                item = names[name]
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += item[0]*dervDict[item[1]]
                    if Ka2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += item[0]*dervDict[item[1]]
                    if Ka2:
                        depDerivDict[name][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
            for iPO in dIdPO:
                if iPO in varylist:
                    dMdv[varylist.index(iPO)][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
                    if Ka2:
                        dMdv[varylist.index(iPO)][iBeg2:iFin2] += dIdPO[iPO]*dervDict2['int']
                elif iPO in dependentVars:
                    depDerivDict[iPO][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
                    if Ka2:
                        depDerivDict[iPO][iBeg2:iFin2] += dIdPO[iPO]*dervDict2['int']
            for i,name in enumerate(['omega','chi','phi']):
                aname = pfx+'SH '+name
                if aname in varylist:
                    dMdv[varylist.index(aname)][iBeg:iFin] += dFdSA[i]*dervDict['int']
                    if Ka2:
                        dMdv[varylist.index(aname)][iBeg2:iFin2] += dFdSA[i]*dervDict2['int']
                elif aname in dependentVars:
                    depDerivDict[aname][iBeg:iFin] += dFdSA[i]*dervDict['int']
                    if Ka2:
                        depDerivDict[aname][iBeg2:iFin2] += dFdSA[i]*dervDict2['int']
            for iSH in dFdODF:
                if iSH in varylist:
                    dMdv[varylist.index(iSH)][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
                    if Ka2:
                        dMdv[varylist.index(iSH)][iBeg2:iFin2] += dFdODF[iSH]*dervDict2['int']
                elif iSH in dependentVars:
                    depDerivDict[iSH][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
                    if Ka2:
                        depDerivDict[iSH][iBeg2:iFin2] += dFdODF[iSH]*dervDict2['int']
            cellDervNames = cellVaryDerv(pfx,SGData,dpdA)
            for name,dpdA in cellDervNames:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dpdA*dervDict['pos']
                    if Ka2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += dpdA*dervDict2['pos']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += dpdA*dervDict['pos']
                    if Ka2:
                        depDerivDict[name][iBeg2:iFin2] += dpdA*dervDict2['pos']
            dDijDict = GetHStrainShiftDerv(refl,SGData,phfx,hfx,calcControls,parmDict)
            for name in dDijDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                    if Ka2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                    if Ka2:
                        depDerivDict[name][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
            if 'C' in calcControls[hfx+'histType']:
                sigDict,gamDict = GetSampleSigGamDerv(refl,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict)
            else:   #'T'OF
                sigDict,gamDict = GetSampleSigGamDerv(refl,0.0,G,GB,SGData,hfx,phfx,calcControls,parmDict)
            for name in gamDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += gamDict[name]*dervDict['gam']
                    if Ka2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += gamDict[name]*dervDict2['gam']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += gamDict[name]*dervDict['gam']
                    if Ka2:
                        depDerivDict[name][iBeg2:iFin2] += gamDict[name]*dervDict2['gam']
            for name in sigDict:
                if name in varylist:
                    dMdv[varylist.index(name)][iBeg:iFin] += sigDict[name]*dervDict['sig']
                    if Ka2:
                        dMdv[varylist.index(name)][iBeg2:iFin2] += sigDict[name]*dervDict2['sig']
                elif name in dependentVars:
                    depDerivDict[name][iBeg:iFin] += sigDict[name]*dervDict['sig']
                    if Ka2:
                        depDerivDict[name][iBeg2:iFin2] += sigDict[name]*dervDict2['sig']
            for name in ['BabA','BabU']:
                if refl[9]:
                    if phfx+name in varylist:
                        dMdv[varylist.index(phfx+name)][iBeg:iFin] += dFdvDict[pfx+name][iref]*dervDict['int']/refl[9]
                        if Ka2:
                            dMdv[varylist.index(phfx+name)][iBeg2:iFin2] += dFdvDict[pfx+name][iref]*dervDict2['int']/refl[9]
                    elif phfx+name in dependentVars:                    
                        depDerivDict[phfx+name][iBeg:iFin] += dFdvDict[pfx+name][iref]*dervDict['int']/refl[9]
                        if Ka2:
                            depDerivDict[phfx+name][iBeg2:iFin2] += dFdvDict[pfx+name][iref]*dervDict2['int']/refl[9]                  
            if not Phase['General'].get('doPawley'):
                #do atom derivatives -  for RB,F,X & U so far              
                corr = dervDict['int']/refl[9]
                if Ka2:
                    corr2 = dervDict2['int']/refl[9]
                for name in varylist+dependentVars:
                    if '::RBV;' in name:
                        pass
                    else:
                        try:
                            aname = name.split(pfx)[1][:2]
                            if aname not in ['Af','dA','AU','RB']: continue # skip anything not an atom or rigid body param
                        except IndexError:
                            continue
                    if name in varylist:
                        dMdv[varylist.index(name)][iBeg:iFin] += dFdvDict[name][iref]*corr
                        if Ka2:
                            dMdv[varylist.index(name)][iBeg2:iFin2] += dFdvDict[name][iref]*corr2
                    elif name in dependentVars:
                        depDerivDict[name][iBeg:iFin] += dFdvDict[name][iref]*corr
                        if Ka2:
                            depDerivDict[name][iBeg2:iFin2] += dFdvDict[name][iref]*corr2
    #        print 'profile derv time: %.3fs'%(time.time()-time0)
    # now process derivatives in constraints
    G2mv.Dict2Deriv(varylist,depDerivDict,dMdv)
    return dMdv
    
def dervHKLF(Histogram,Phase,calcControls,varylist,parmDict,rigidbodyDict):
    '''Loop over reflections in a HKLF histogram and compute derivatives of the fitting
    model (M) with respect to all parameters.  Independent and dependant dM/dp arrays 
    are returned to either dervRefine or HessRefine.

    :returns: 
    '''
    nobs = Histogram['Residuals']['Nobs']
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    pfx = '%d::'%(Phase['pId'])
    phfx = '%d:%d:'%(Phase['pId'],hId)
    SGData = Phase['General']['SGData']
    A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
    G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
    refDict = Histogram['Data']
    dFdvDict = StructureFactorDerv(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
    ApplyRBModelDervs(dFdvDict,parmDict,rigidbodyDict,Phase)
    dMdvh = np.zeros((len(varylist),len(refDict['RefList'])))
    dependentVars = G2mv.GetDependentVars()
    depDerivDict = {}
    for j in dependentVars:
        depDerivDict[j] = np.zeros(shape=(len(refDict['RefList'])))
    wdf = np.zeros(len(refDict['RefList']))
    if calcControls['F**2']:
        for iref,ref in enumerate(refDict['RefList']):
            if ref[6] > 0:
                dervDict = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmDict,varylist+dependentVars)[1] 
                w = 1.0/ref[6]
                if w*ref[5] >= calcControls['minF/sig']:
                    wdf[iref] = w*(ref[5]-ref[7])
                    for j,var in enumerate(varylist):
                        if var in dFdvDict:
                            dMdvh[j][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11]
                    for var in dependentVars:
                        if var in dFdvDict:
                            depDerivDict[var][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11]
                    if phfx+'Scale' in varylist:
                        dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*ref[11]
                    elif phfx+'Scale' in dependentVars:
                        depDerivDict[phfx+'Scale'][iref] = w*ref[9]*ref[11]
                    for item in ['Ep','Es','Eg']:
                        if phfx+item in varylist and phfx+item in dervDict:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]/ref[11]  #OK
                        elif phfx+item in dependentVars and phfx+item in dervDict:
                            depDerivDict[phfx+item][iref] = w*dervDict[phfx+item]/ref[11]  #OK
                    for item in ['BabA','BabU']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dFdvDict[pfx+item][iref]*parmDict[phfx+'Scale']*ref[11]
                        elif phfx+item in dependentVars:
                            depDerivDict[phfx+item][iref] = w*dFdvDict[pfx+item][iref]*parmDict[phfx+'Scale']*ref[11]
    else:
        for iref,ref in enumerate(refDict['RefList']):
            if ref[5] > 0.:
                dervDict = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmDict,varylist+dependentVars)[1]
                Fo = np.sqrt(ref[5])
                Fc = np.sqrt(ref[7])
                w = 1.0/ref[6]
                if 2.0*Fo*w*Fo >= calcControls['minF/sig']:
                    wdf[iref] = 2.0*Fo*w*(Fo-Fc)
                    for j,var in enumerate(varylist):
                        if var in dFdvDict:
                            dMdvh[j][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11]
                    for var in dependentVars:
                        if var in dFdvDict:
                            depDerivDict[var][iref] = w*dFdvDict[var][iref]*parmDict[phfx+'Scale']*ref[11]
                    if phfx+'Scale' in varylist:
                        dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*ref[11]
                    elif phfx+'Scale' in dependentVars:
                        depDerivDict[phfx+'Scale'][iref] = w*ref[9]*ref[11]                           
                    for item in ['Ep','Es','Eg']:
                        if phfx+item in varylist and phfx+item in dervDict:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]/ref[11]  #correct
                        elif phfx+item in dependentVars and phfx+item in dervDict:
                            depDerivDict[phfx+item][iref] = w*dervDict[phfx+item]/ref[11]
                    for item in ['BabA','BabU']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dFdvDict[pfx+item][iref]*parmDict[phfx+'Scale']*ref[11]
                        elif phfx+item in dependentVars:
                            depDerivDict[phfx+item][iref] = w*dFdvDict[pfx+item][iref]*parmDict[phfx+'Scale']*ref[11]
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
    nvar = len(varylist)
    dMdv = np.empty(0)
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            dMdvh = np.sqrt(w[xB:xF])*getPowderProfileDerv(parmDict,x[xB:xF],
                varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup)
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
        if len(dMdv):
            dMdv = np.concatenate((dMdv.T,np.sqrt(wtFactor)*dMdvh.T)).T
        else:
            dMdv = np.sqrt(wtFactor)*dMdvh
            
    pNames,pVals,pWt,pWsum = penaltyFxn(HistoPhases,parmDict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,parmDict,varylist)
        dMdv = np.concatenate((dMdv.T,(np.sqrt(pWt)*dpdv).T)).T
        
    return dMdv

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
    ApplyRBModels(parmDict,Phases,rigidbodyDict)        #,Update=True??
    nvar = len(varylist)
    Hess = np.empty(0)
    histoList = Histograms.keys()
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
            xF = np.searchsorted(x,Limits[1])
            dMdvh = getPowderProfileDerv(parmDict,x[xB:xF],
                varylist,Histogram,Phases,rigidbodyDict,calcControls,pawleyLookup)
            Wt = ma.sqrt(W[xB:xF])[np.newaxis,:]
            Dy = dy[xB:xF][np.newaxis,:]
            dMdvh *= Wt
            if dlg:
                dlg.Update(Histogram['Residuals']['wR'],newmsg='Hessian for histogram %d\nAll data Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))[0]
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
                dlg.Update(Histogram['Residuals']['wR'],newmsg='Hessian for histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))[0]
            if len(Hess):
                Vec += wtFactor*np.sum(dMdvh*wdf,axis=1)
                Hess += wtFactor*np.inner(dMdvh,dMdvh)
            else:
                Vec = wtFactor*np.sum(dMdvh*wdf,axis=1)
                Hess = wtFactor*np.inner(dMdvh,dMdvh)
        else:
            continue        #skip non-histogram entries
    pNames,pVals,pWt,pWsum = penaltyFxn(HistoPhases,parmDict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,parmDict,varylist)
        Vec += np.sum(dpdv*pWt*pVals,axis=1)
        Hess += np.inner(dpdv*pWt,dpdv)
    return Vec,Hess

def errRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,dlg):        
    'Needs a doc string'
    Values2Dict(parmDict, varylist, values)
    G2mv.Dict2Map(parmDict,varylist)
    Histograms,Phases,restraintDict,rigidbodyDict = HistoPhases
    M = np.empty(0)
    SumwYo = 0
    Nobs = 0
    ApplyRBModels(parmDict,Phases,rigidbodyDict)
    histoList = Histograms.keys()
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
            xF = np.searchsorted(x,Limits[1])
            yc[xB:xF],yb[xB:xF] = getPowderProfile(parmDict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            yc[xB:xF] += yb[xB:xF]
            if not np.any(y):                   #fill dummy data
                rv = st.poisson(yc[xB:xF])
                y[xB:xF] = rv.rvs()
                Z = np.ones_like(yc[xB:xF])
                Z[1::2] *= -1
                y[xB:xF] = yc[xB:xF]+np.abs(y[xB:xF]-yc[xB:xF])*Z
                w[xB:xF] = np.where(y[xB:xF]>0.,1./y[xB:xF],1.0)
            yd[xB:xF] = y[xB:xF]-yc[xB:xF]
            W = wtFactor*w
            wdy = -ma.sqrt(W[xB:xF])*(yd[xB:xF])
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
                dlg.Update(Histogram['Residuals']['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))[0]
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
            A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
            G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
            refDict = Histogram['Data']
            time0 = time.time()
            StructureFactor2(refDict,G,hfx,pfx,SGData,calcControls,parmDict)
#            print 'sf-calc time: %.3f'%(time.time()-time0)
            df = np.zeros(len(refDict['RefList']))
            sumwYo = 0
            sumFo = 0
            sumFo2 = 0
            sumdF = 0
            sumdF2 = 0
            nobs = 0
            if calcControls['F**2']:
                for i,ref in enumerate(refDict['RefList']):
                    if ref[6] > 0:
                        ref[11] = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmDict,varylist)[0]
                        w = 1.0/ref[6]
                        ref[7] = parmDict[phfx+'Scale']*ref[9]*ref[11]  #correct Fc^2 for extinction
                        ref[8] = ref[5]/(parmDict[phfx+'Scale']*ref[11])
                        if w*ref[5] >= calcControls['minF/sig']:
                            Fo = np.sqrt(ref[5])
                            sumFo += Fo
                            sumFo2 += ref[5]
                            sumdF += abs(Fo-np.sqrt(ref[7]))
                            sumdF2 += abs(ref[5]-ref[7])
                            nobs += 1
                            df[i] = -w*(ref[5]-ref[7])
                            sumwYo += (w*ref[5])**2
            else:
                for i,ref in enumerate(refDict['RefList']):
                    if ref[5] > 0.:
                        ref[11] = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmDict,varylist)[0]
                        ref[7] = parmDict[phfx+'Scale']*ref[9]*ref[11]    #correct Fc^2 for extinction
                        ref[8] = ref[5]/(parmDict[phfx+'Scale']*ref[11])
                        Fo = np.sqrt(ref[5])
                        Fc = np.sqrt(ref[7])
                        w = 2.0*Fo/ref[6]
                        if w*Fo >= calcControls['minF/sig']:
                            sumFo += Fo
                            sumFo2 += ref[5]
                            sumdF += abs(Fo-Fc)
                            sumdF2 += abs(ref[5]-ref[7])
                            nobs += 1
                            df[i] = -w*(Fo-Fc)
                            sumwYo += (w*Fo)**2
            Histogram['Residuals']['Nobs'] = nobs
            Histogram['Residuals']['sumwYo'] = sumwYo
            SumwYo += sumwYo
            Histogram['Residuals']['wR'] = min(100.,np.sqrt(np.sum(df**2)/Histogram['Residuals']['sumwYo'])*100.)
            Histogram['Residuals'][phfx+'Rf'] = 100.*sumdF/sumFo
            Histogram['Residuals'][phfx+'Rf^2'] = 100.*sumdF2/sumFo2
            Histogram['Residuals'][phfx+'Nref'] = nobs
            Nobs += nobs
            if dlg:
                dlg.Update(Histogram['Residuals']['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['Residuals']['wR'],'%'))[0]
            M = np.concatenate((M,wtFactor*df))
# end of HKLF processing
    Histograms['sumwYo'] = SumwYo
    Histograms['Nobs'] = Nobs
    Rw = min(100.,np.sqrt(np.sum(M**2)/SumwYo)*100.)
    if dlg:
        GoOn = dlg.Update(Rw,newmsg='%s%8.3f%s'%('All data Rw =',Rw,'%'))[0]
        if not GoOn:
            parmDict['saved values'] = values
            dlg.Destroy()
            raise Exception         #Abort!!
    pDict,pVals,pWt,pWsum = penaltyFxn(HistoPhases,parmDict,varylist)
    if len(pVals):
        pSum = np.sum(pWt*pVals**2)
        for name in pWsum:
            if pWsum:
                print '  Penalty function for %8s = %12.5g'%(name,pWsum[name])
        print 'Total penalty function: %12.5g on %d terms'%(pSum,len(pVals))
        Nobs += len(pVals)
        M = np.concatenate((M,np.sqrt(pWt)*pVals))
    return M
                        
