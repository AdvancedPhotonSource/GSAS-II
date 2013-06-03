# -*- coding: utf-8 -*-
'''
*GSASIIstrMain: main structure routine*
---------------------------------------

'''
########### SVN repository information ###################
# $Date: 2013-05-17 12:13:15 -0500 (Fri, 17 May 2013) $
# $Author: vondreele $
# $Revision: 920 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIIstrMain.py $
# $Id: GSASIIstrMain.py 920 2013-05-17 17:13:15Z vondreele $
########### SVN repository information ###################
import sys
import os
import os.path as ospath
import time
import math
import copy
import random
import cPickle
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import scipy.optimize as so
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 920 $")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIstrIO as G2stIO
import GSASIIstrMath as G2stMth

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ateln2 = 8.0*math.log(2.0)
DEBUG = False


def Refine(GPXfile,dlg):
    'Needs a doc string'
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    G2stIO.ShowBanner(printFile)
    varyList = []
    parmDict = {}
    G2mv.InitVars()    
    Controls = G2stIO.GetControls(GPXfile)
    G2stIO.ShowControls(Controls,printFile)
    calcControls = {}
    calcControls.update(Controls)            
    constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no phases! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception        
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,pFile=printFile)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,pFile=printFile)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,pFile=printFile)
    calcControls.update(controlDict)
    varyList = rbVary+phaseVary+hapVary+histVary
    parmDict.update(rbDict)
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    G2stIO.GetFprime(calcControls,Histograms)
    # do constraint processing
    try:
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList)
    except:
        print ' *** ERROR - your constraints are internally inconsistent ***'
        errmsg, warnmsg = G2stIO.CheckConstraints(varylist,constrDict,fixedList)
        print 'Errors',errmsg
        if warnmsg: print 'Warnings',warnmsg
        raise Exception(' *** Refine aborted ***')
    # # check to see which generated parameters are fully varied
    # msg = G2mv.SetVaryFlags(varyList)
    # if msg:
    #     print ' *** ERROR - you have not set the refine flags for constraints consistently! ***'
    #     print msg
    #     raise Exception(' *** Refine aborted ***')
    #print G2mv.VarRemapShow(varyList)
    G2mv.Map2Dict(parmDict,varyList)
    Rvals = {}
    while True:
        begin = time.time()
        values =  np.array(G2stMth.Dict2Values(parmDict, varyList))
        Ftol = Controls['min dM/M']
        Factor = Controls['shift factor']
        maxCyc = Controls['max cyc']
        if 'Jacobian' in Controls['deriv type']:            
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
        elif 'Hessian' in Controls['deriv type']:
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,maxcyc=maxCyc,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
        else:           #'numeric'
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/len(varyList))
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        G2stMth.Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict,varyList)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = Rvals['chisq']/(Histograms['Nobs']-len(varyList))
        print >>printFile,'\n Refinement results:'
        print >>printFile,135*'-'
        print >>printFile,' Number of function calls:',result[2]['nfev'],' Number of observations: ',Histograms['Nobs'],' Number of parameters: ',len(varyList)
        print >>printFile,' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles'%(runtime,runtime/ncyc,ncyc)
        print >>printFile,' wR = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
        try:
            covMatrix = result[1]*Rvals['GOF']
            sig = np.sqrt(np.diag(covMatrix))
            if np.any(np.isnan(sig)):
                print '*** Least squares aborted - some invalid esds possible ***'
#            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except TypeError:          #result[1] is None on singular matrix
            print '**** Refinement failed - singular matrix ****'
            if 'Hessian' in Controls['deriv type']:
                num = len(varyList)-1
                for i,val in enumerate(np.flipud(result[2]['psing'])):
                    if val:
                        print 'Removing parameter: ',varyList[num-i]
                        del(varyList[num-i])                    
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        print 'Removing parameter: ',varyList[ipvt-1]
                        del(varyList[ipvt-1])
                        break

#    print 'dependentParmList: ',G2mv.dependentParmList
#    print 'arrayList: ',G2mv.arrayList
#    print 'invarrayList: ',G2mv.invarrayList
#    print 'indParmList: ',G2mv.indParmList
#    print 'fixedDict: ',G2mv.fixedDict
#    print 'test1'
    G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
#    print 'test2'
    sigDict = dict(zip(varyList,sig))
    newCellDict = G2stMth.GetNewCellParms(parmDict,varyList)
    newAtomDict = G2stMth.ApplyXYZshifts(parmDict,varyList)
    covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
        'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
    # add the uncertainties into the esd dictionary (sigDict)
    sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
    G2mv.PrintIndependentVars(parmDict,varyList,sigDict,pFile=printFile)
    G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
    G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
    G2stIO.SetPhaseData(parmDict,sigDict,Phases,rbIds,covData,restraintDict,printFile)
    G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,pFile=printFile)
    G2stIO.SetHistogramData(parmDict,sigDict,Histograms,pFile=printFile)
    G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,rigidbodyDict,covData)
    printFile.close()
    print ' Refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst'
    print ' ***** Refinement successful *****'
    
#for testing purposes!!!
    if DEBUG:
#needs: values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup
        import cPickle
        fl = open('testDeriv.dat','wb')
        cPickle.dump(result[0],fl,1)
        cPickle.dump([Histograms,Phases,restraintDict,rigidbodyDict],fl,1)
        cPickle.dump(parmDict,fl,1)
        cPickle.dump(varyList,fl,1)
        cPickle.dump(calcControls,fl,1)
        cPickle.dump(pawleyLookup,fl,1)
        fl.close()

    if dlg:
        return Rvals['Rwp']

def SeqRefine(GPXfile,dlg):
    'Needs a doc string'
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    print ' Sequential Refinement'
    G2stIO.ShowBanner(printFile)
    G2mv.InitVars()    
    Controls = G2stIO.GetControls(GPXfile)
    G2stIO.ShowControls(Controls,printFile)            
    constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no histograms to refine! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,False,printFile)
    for item in phaseVary:
        if '::A0' in item:
            print '**** WARNING - lattice parameters should not be refined in a sequential refinement ****'
            print '****           instead use the Dij parameters for each powder histogram            ****'
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else:
        histNames = G2stIO.GetHistogramNames(GPXfile,['PWDR',])
    if 'Reverse Seq' in Controls:
        if Controls['Reverse Seq']:
            histNames.reverse()
    SeqResult = {'histNames':histNames}
    makeBack = True
    for ihst,histogram in enumerate(histNames):
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['atomIndx'] = atomIndx
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['BLtables'] = BLtables
        varyList = []
        parmDict = {}
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histo,False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = rbVary+phaseVary+hapVary+histVary
        if not ihst:
            saveVaryList = varyList[:]
            for i,item in enumerate(saveVaryList):
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                saveVaryList[i] = item
            SeqResult['varyList'] = saveVaryList
        else:
            newVaryList = varyList[:]
            for i,item in enumerate(newVaryList):
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                newVaryList[i] = item
            if newVaryList != SeqResult['varyList']:
                print newVaryList
                print SeqResult['varyList']
                print '**** ERROR - variable list for this histogram does not match previous'
                raise Exception
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        G2stIO.GetFprime(calcControls,Histo)
        # do constraint processing
        try:
            groups,parmlist = G2mv.GroupConstraints(constrDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList)
        except:
            print ' *** ERROR - your constraints are internally inconsistent ***'
            raise Exception(' *** Refine aborted ***')
        # check to see which generated parameters are fully varied
        # msg = G2mv.SetVaryFlags(varyList)
        # if msg:
        #     print ' *** ERROR - you have not set the refine flags for constraints consistently! ***'
        #     print msg
        #     raise Exception(' *** Refine aborted ***')
        #print G2mv.VarRemapShow(varyList)
        G2mv.Map2Dict(parmDict,varyList)
        Rvals = {}
        while True:
            begin = time.time()
            values =  np.array(G2stMth.Dict2Values(parmDict,varyList))
            Ftol = Controls['min dM/M']
            Factor = Controls['shift factor']
            maxCyc = Controls['max cyc']

            if 'Jacobian' in Controls['deriv type']:            
                result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                    ftol=Ftol,col_deriv=True,factor=Factor,
                    args=([Histo,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/2)
            elif 'Hessian' in Controls['deriv type']:
                result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,maxcyc=maxCyc,
                    args=([Histo,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = result[2]['num cyc']+1                           
            else:           #'numeric'
                result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                    args=([Histo,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/len(varyList))

            runtime = time.time()-begin
            Rvals['chisq'] = np.sum(result[2]['fvec']**2)
            G2stMth.Values2Dict(parmDict, varyList, result[0])
            G2mv.Dict2Map(parmDict,varyList)
            
            Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histo['sumwYo'])*100.      #to %
            Rvals['GOF'] = Rvals['Rwp']/(Histo['Nobs']-len(varyList))
            Rvals['Nobs'] = Histo['Nobs']
            print >>printFile,'\n Refinement results for histogram: v'+histogram
            print >>printFile,135*'-'
            print >>printFile,' Number of function calls:',result[2]['nfev'],' Number of observations: ',Histo['Nobs'],' Number of parameters: ',len(varyList)
            print >>printFile,' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles'%(runtime,runtime/ncyc,ncyc)
            print >>printFile,' wRp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
            try:
                covMatrix = result[1]*Rvals['GOF']
                sig = np.sqrt(np.diag(covMatrix))
                if np.any(np.isnan(sig)):
                    print '*** Least squares aborted - some invalid esds possible ***'
                    ifPrint = True
                break                   #refinement succeeded - finish up!
            except TypeError:          #result[1] is None on singular matrix
                print '**** Refinement failed - singular matrix ****'
                if 'Hessian' in Controls['deriv type']:
                    num = len(varyList)-1
                    for i,val in enumerate(np.flipud(result[2]['psing'])):
                        if val:
                            print 'Removing parameter: ',varyList[num-i]
                            del(varyList[num-i])                    
                else:
                    Ipvt = result[2]['ipvt']
                    for i,ipvt in enumerate(Ipvt):
                        if not np.sum(result[2]['fjac'],axis=1)[i]:
                            print 'Removing parameter: ',varyList[ipvt-1]
                            del(varyList[ipvt-1])
                            break
    
        G2stMth.GetFobsSq(Histo,Phases,parmDict,calcControls)
        sigDict = dict(zip(varyList,sig))
        newCellDict = copy.deepcopy(G2stMth.GetNewCellParms(parmDict,varyList))
        newAtomDict = copy.deepcopy(G2stMth.ApplyXYZshifts(parmDict,varyList))
        covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
            'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
        # add the uncertainties into the esd dictionary (sigDict)
        G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
#??        SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
        G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,ifPrint,printFile)
        G2stIO.SetHistogramData(parmDict,sigDict,Histo,ifPrint,printFile)
        SeqResult[histogram] = covData
        G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histo,Phases,rigidbodyDict,covData,makeBack)
        makeBack = False
    G2stIO.SetSeqResult(GPXfile,Histograms,SeqResult)
    printFile.close()
    print ' Sequential refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst'
    print ' ***** Sequential refinement successful *****'

def DistAngle(DisAglCtls,DisAglData):
    'Needs a doc string'
    import numpy.ma as ma
    
    def ShowBanner(name):
        print 80*'*'
        print '   Interatomic Distances and Angles for phase '+name
        print 80*'*','\n'

    ShowBanner(DisAglCtls['Name'])
    SGData = DisAglData['SGData']
    SGtext = G2spc.SGPrint(SGData)
    for line in SGtext: print line
    Cell = DisAglData['Cell']
    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if 'covData' in DisAglData:   
        covData = DisAglData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(DisAglData['pId'])+'::'
        A = G2lat.cell2A(Cell[:6])
        cellSig = G2stIO.getCellEsd(pfx,SGData,A,covData)
        names = [' a = ',' b = ',' c = ',' alpha = ',' beta = ',' gamma = ',' Volume = ']
        valEsd = [G2mth.ValEsd(Cell[i],cellSig[i],True) for i in range(7)]
        line = '\n Unit cell:'
        for name,vals in zip(names,valEsd):
            line += name+vals  
        print line
    else: 
        print '\n Unit cell: a = ','%.5f'%(Cell[0]),' b = ','%.5f'%(Cell[1]),' c = ','%.5f'%(Cell[2]), \
            ' alpha = ','%.3f'%(Cell[3]),' beta = ','%.3f'%(Cell[4]),' gamma = ', \
            '%.3f'%(Cell[5]),' volume = ','%.3f'%(Cell[6])
    Factor = DisAglCtls['Factors']
    Radii = dict(zip(DisAglCtls['AtomTypes'],zip(DisAglCtls['BondRadii'],DisAglCtls['AngleRadii'])))
    indices = (-1,0,1)
    Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    origAtoms = DisAglData['OrigAtoms']
    targAtoms = DisAglData['TargAtoms']
    for Oatom in origAtoms:
        OxyzNames = ''
        IndBlist = []
        Dist = []
        Vect = []
        VectA = []
        angles = []
        for Tatom in targAtoms:
            Xvcov = []
            TxyzNames = ''
            if 'covData' in DisAglData:
                OxyzNames = [pfx+'dAx:%d'%(Oatom[0]),pfx+'dAy:%d'%(Oatom[0]),pfx+'dAz:%d'%(Oatom[0])]
                TxyzNames = [pfx+'dAx:%d'%(Tatom[0]),pfx+'dAy:%d'%(Tatom[0]),pfx+'dAz:%d'%(Tatom[0])]
                Xvcov = G2mth.getVCov(OxyzNames+TxyzNames,varyList,covMatrix)
            result = G2spc.GenAtom(Tatom[3:6],SGData,False,Move=False)
            BsumR = (Radii[Oatom[2]][0]+Radii[Tatom[2]][0])*Factor[0]
            AsumR = (Radii[Oatom[2]][1]+Radii[Tatom[2]][1])*Factor[1]
            for Txyz,Top,Tunit in result:
                Dx = (Txyz-np.array(Oatom[3:6]))+Units
                dx = np.inner(Amat,Dx)
                dist = ma.masked_less(np.sqrt(np.sum(dx**2,axis=0)),0.5)
                IndB = ma.nonzero(ma.masked_greater(dist-BsumR,0.))
                if np.any(IndB):
                    for indb in IndB:
                        for i in range(len(indb)):
                            if str(dx.T[indb][i]) not in IndBlist:
                                IndBlist.append(str(dx.T[indb][i]))
                                unit = Units[indb][i]
                                tunit = '[%2d%2d%2d]'%(unit[0]+Tunit[0],unit[1]+Tunit[1],unit[2]+Tunit[2])
                                pdpx = G2mth.getDistDerv(Oatom[3:6],Tatom[3:6],Amat,unit,Top,SGData)
                                sig = 0.0
                                if len(Xvcov):
                                    sig = np.sqrt(np.inner(pdpx,np.inner(Xvcov,pdpx)))
                                Dist.append([Oatom[1],Tatom[1],tunit,Top,ma.getdata(dist[indb])[i],sig])
                                if (Dist[-1][-2]-AsumR) <= 0.:
                                    Vect.append(dx.T[indb][i]/Dist[-1][-2])
                                    VectA.append([OxyzNames,np.array(Oatom[3:6]),TxyzNames,np.array(Tatom[3:6]),unit,Top])
                                else:
                                    Vect.append([0.,0.,0.])
                                    VectA.append([])
        Vect = np.array(Vect)
        angles = np.zeros((len(Vect),len(Vect)))
        angsig = np.zeros((len(Vect),len(Vect)))
        for i,veca in enumerate(Vect):
            if np.any(veca):
                for j,vecb in enumerate(Vect):
                    if np.any(vecb):
                        angles[i][j],angsig[i][j] = G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)
        line = ''
        for i,x in enumerate(Oatom[3:6]):
            if len(Xvcov):
                line += '%12s'%(G2mth.ValEsd(x,np.sqrt(Xvcov[i][i]),True))
            else:
                line += '%12.5f'%(x)
        print '\n Distances & angles for ',Oatom[1],' at ',line
        print 80*'*'
        line = ''
        for dist in Dist[:-1]:
            line += '%12s'%(dist[1].center(12))
        print '  To       cell +(sym. op.)      dist.  ',line
        for i,dist in enumerate(Dist):
            line = ''
            for j,angle in enumerate(angles[i][0:i]):
                sig = angsig[i][j]
                if angle:
                    if sig:
                        line += '%12s'%(G2mth.ValEsd(angle,sig,True).center(12))
                    else:
                        val = '%.3f'%(angle)
                        line += '%12s'%(val.center(12))
                else:
                    line += 12*' '
            if dist[5]:            #sig exists!
                val = G2mth.ValEsd(dist[4],dist[5])
            else:
                val = '%8.4f'%(dist[4])
            print '  %8s%10s+(%4d) %12s'%(dist[1].ljust(8),dist[2].ljust(10),dist[3],val.center(12)),line

def DisAglTor(DATData):
    'Needs a doc string'
    SGData = DATData['SGData']
    Cell = DATData['Cell']
    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    pfx = ''
    if 'covData' in DATData:   
        covData = DATData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(DATData['pId'])+'::'
    Datoms = []
    Oatoms = []
    for i,atom in enumerate(DATData['Datoms']):
        symop = atom[-1].split('+')
        if len(symop) == 1:
            symop.append('0,0,0')        
        symop[0] = int(symop[0])
        symop[1] = eval(symop[1])
        atom.append(symop)
        Datoms.append(atom)
        oatom = DATData['Oatoms'][i]
        names = ['','','']
        if pfx:
            names = [pfx+'dAx:'+str(oatom[0]),pfx+'dAy:'+str(oatom[0]),pfx+'dAz:'+str(oatom[0])]
        oatom += [names,]
        Oatoms.append(oatom)
    atmSeq = [atom[1]+'('+atom[-2]+')' for atom in Datoms]
    if DATData['Natoms'] == 4:  #torsion
        Tors,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print ' Torsion angle for '+DATData['Name']+' atom sequence: ',atmSeq,'=',G2mth.ValEsd(Tors,sig)
        print ' NB: Atom sequence determined by selection order'
        return      # done with torsion
    elif DATData['Natoms'] == 3:  #angle
        Ang,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print ' Angle in '+DATData['Name']+' for atom sequence: ',atmSeq,'=',G2mth.ValEsd(Ang,sig)
        print ' NB: Atom sequence determined by selection order'
    else:   #2 atoms - distance
        Dist,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print ' Distance in '+DATData['Name']+' for atom sequence: ',atmSeq,'=',G2mth.ValEsd(Dist,sig)
                
def BestPlane(PlaneData):
    'Needs a doc string'

    def ShowBanner(name):
        print 80*'*'
        print '   Best plane result for phase '+name
        print 80*'*','\n'

    ShowBanner(PlaneData['Name'])

    Cell = PlaneData['Cell']    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])        
    Atoms = PlaneData['Atoms']
    sumXYZ = np.zeros(3)
    XYZ = []
    Natoms = len(Atoms)
    for atom in Atoms:
        xyz = np.array(atom[3:6])
        XYZ.append(xyz)
        sumXYZ += xyz
    sumXYZ /= Natoms
    XYZ = np.array(XYZ)-sumXYZ
    XYZ = np.inner(Amat,XYZ).T
    Zmat = np.zeros((3,3))
    for i,xyz in enumerate(XYZ):
        Zmat += np.outer(xyz.T,xyz)
    print ' Selected atoms centered at %10.5f %10.5f %10.5f'%(sumXYZ[0],sumXYZ[1],sumXYZ[2])
    Evec,Emat = nl.eig(Zmat)
    Evec = np.sqrt(Evec)/(Natoms-3)
    Order = np.argsort(Evec)
    XYZ = np.inner(XYZ,Emat.T).T
    XYZ = np.array([XYZ[Order[2]],XYZ[Order[1]],XYZ[Order[0]]]).T
    print ' Atoms in Cartesian best plane coordinates:'
    print ' Name         X         Y         Z'
    for i,xyz in enumerate(XYZ):
        print ' %6s%10.3f%10.3f%10.3f'%(Atoms[i][1].ljust(6),xyz[0],xyz[1],xyz[2])
    print '\n Best plane RMS X =%8.3f, Y =%8.3f, Z =%8.3f'%(Evec[Order[2]],Evec[Order[1]],Evec[Order[0]])   

            
def main():
    'Needs a doc string'
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print 'ERROR - ',GPXfile," doesn't exist!"
            exit()
        GPXpath = ospath.dirname(arg[1])
        Refine(GPXfile,None)
    else:
        print 'ERROR - missing filename'
        exit()
    print "Done"
         
if __name__ == '__main__':
    main()
