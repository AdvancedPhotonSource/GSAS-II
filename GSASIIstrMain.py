# -*- coding: utf-8 -*-
'''
*GSASIIstrMain: main structure routine*
---------------------------------------

'''
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
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIstrIO as G2stIO
import GSASIIstrMath as G2stMth
import GSASIIobj as G2obj

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ateln2 = 8.0*math.log(2.0)
DEBUG = True

def RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
    calcControls,pawleyLookup,ifPrint,printFile,dlg):
    'Core optimization routines, shared between SeqRefine and Refine'
#    print 'current',varyList
#    for item in parmDict: print item,parmDict[item] ######### show dict just before refinement
    G2mv.Map2Dict(parmDict,varyList)
    Rvals = {}
    while True:
        begin = time.time()
        values =  np.array(G2stMth.Dict2Values(parmDict, varyList))
        # test code to compute GOF and save for external repeat
        #args = ([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg)
        #print '*** before fit chi**2',np.sum(G2stMth.errRefine(values,*args)**2)            
        #fl = open('beforeFit.cpickle','wb')
        #import cPickle
        #cPickle.dump(values,fl,1)
        #cPickle.dump(args[:-1],fl,1)
        #fl.close()
        Ftol = Controls['min dM/M']
        Factor = Controls['shift factor']
        if 'Jacobian' in Controls['deriv type']:            
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
        elif 'Hessian' in Controls['deriv type']:
            maxCyc = Controls['max cyc']
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,maxcyc=maxCyc,Print=ifPrint,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
        else:           #'numeric'
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = 1
            if len(varyList):
                ncyc = int(result[2]['nfev']/len(varyList))
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['converged'] = result[2].get('Converged')
        Rvals['DelChi2'] = result[2].get('DelChi2',-1.)
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        G2stMth.Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict,varyList)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs']-len(varyList)))
        print >>printFile,' Number of function calls:',result[2]['nfev'],   \
            ' No. of observations: ',Histograms['Nobs'],' No. of parameters: ',len(varyList),   \
            ' User rejected: ',Histograms['Nrej'],' Sp. gp. extinct: ',Histograms['Next']
        print >>printFile,' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles'%(runtime,runtime/ncyc,ncyc)
        print >>printFile,' wR = %7.2f%%, chi**2 = %12.6g, GOF = %6.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
        IfOK = True
        try:
            covMatrix = result[1]*Rvals['GOF']**2
            sig = np.sqrt(np.diag(covMatrix))
            if np.any(np.isnan(sig)):
                print '*** Least squares aborted - some invalid esds possible ***'
#            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except TypeError,FloatingPointError:          #result[1] is None on singular matrix
            IfOK = False
            if not len(varyList):
                covMatrix = []
                sig = []
                break
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
    G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
    return IfOK,Rvals,result,covMatrix,sig

def Refine(GPXfile,dlg):
    'Global refinement -- refines to minimize against all histograms'
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
        print ' *** ERROR - you have no phases to refine! ***'
        print ' *** Refine aborted ***'
        return False,'No phases'
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,maxSSwave = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,pFile=printFile)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    calcControls['maxSSwave'] = maxSSwave
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
    varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed
    try:
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict)
    except:
        print ' *** ERROR - your constraints are internally inconsistent ***'
        #errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
        #print 'Errors',errmsg
        #if warnmsg: print 'Warnings',warnmsg
        return False,' Constraint error'
#    print G2mv.VarRemapShow(varyList)
    
    ifPrint = True
    print >>printFile,'\n Refinement results:'
    print >>printFile,135*'-'
    try:
        IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histograms,Phases,restraintDict,
            rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifPrint,printFile,dlg)
        sigDict = dict(zip(varyList,sig))
        newCellDict = G2stMth.GetNewCellParms(parmDict,varyList)
        newAtomDict = G2stMth.ApplyXYZshifts(parmDict,varyList)
        covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                   'varyListStart':varyListStart,
                   'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,
                   'newCellDict':newCellDict,'freshCOV':True}
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
    except G2stMth.UserAbort:
        printFile.close()
        return False,'Refinement aborted by user'
    
#for testing purposes!!!
    if DEBUG:
#needs: values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup
        import cPickle
        fl = open('testDeriv.dat','wb')
        cPickle.dump(result[0],fl,1)
        cPickle.dump([Histograms,Phases,restraintDict,rigidbodyDict],fl,1)
        cPickle.dump([G2mv.dependentParmList,G2mv.arrayList,G2mv.invarrayList,
            G2mv.indParmList,G2mv.invarrayList],fl,1)
        cPickle.dump(parmDict,fl,1)
        cPickle.dump(varyList,fl,1)
        cPickle.dump(calcControls,fl,1)
        cPickle.dump(pawleyLookup,fl,1)
        fl.close()

    if dlg:
        return True,Rvals['Rwp']

def SeqRefine(GPXfile,dlg):
    '''Perform a sequential refinement -- cycles through all selected histgrams,
    one at a time
    '''
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    print 'Starting Sequential Refinement'
    G2stIO.ShowBanner(printFile)
    Controls = G2stIO.GetControls(GPXfile)
    G2stIO.ShowControls(Controls,printFile,SeqRef=True)            
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no phases to refine! ***'
        print ' *** Refine aborted ***'
        return False,'No phases'
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,maxSSwave = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,False,printFile)
    for item in phaseVary:
        if '::A0' in item:
            print '**** WARNING - lattice parameters should not be refined in a sequential refinement ****'
            print '****           instead use the Dij parameters for each powder histogram            ****'
            return False,'Lattice parameter refinement error - see console message'
        if '::C(' in item:
            print '**** WARNING - phase texture parameters should not be refined in a sequential refinement ****'
            print '****           instead use the C(L,N) parameters for each powder histogram               ****'
            return False,'Phase texture refinement error - see console message'
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else:
        histNames = G2stIO.GetHistogramNames(GPXfile,['PWDR',])
    if 'Reverse Seq' in Controls:
        if Controls['Reverse Seq']:
            histNames.reverse()
    SeqResult = {'histNames':histNames}
    makeBack = True
    Histo = {}
    NewparmDict = {}
    for ihst,histogram in enumerate(histNames):
        print('Refining with '+str(histogram))
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['atomIndx'] = atomIndx
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['BLtables'] = BLtables
        calcControls['maxSSwave'] = maxSSwave
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histo,Print=False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = rbVary+phaseVary+hapVary+histVary
        if not ihst:
            # save the initial vary list, but without histogram numbers on parameters
            saveVaryList = varyList[:]
            for i,item in enumerate(saveVaryList):
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                saveVaryList[i] = item
            SeqResult['varyList'] = saveVaryList
        origvaryList = varyList[:]
        parmDict = {}
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        if Controls['Copy2Next']:
            parmDict.update(NewparmDict)
        G2stIO.GetFprime(calcControls,Histo)
        # do constraint processing
        #reload(G2mv) # debug
        G2mv.InitVars()    
        constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
        varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed
        try:
            groups,parmlist = G2mv.GroupConstraints(constrDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict,SeqHist=ihst)
            constraintInfo = (groups,parmlist,constrDict,fixedList,ihst)
        except:
            print ' *** ERROR - your constraints are internally inconsistent ***'
            #errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
            #print 'Errors',errmsg
            #if warnmsg: print 'Warnings',warnmsg
            return False,' Constraint error'
        #print G2mv.VarRemapShow(varyList)
        if not ihst:
            # first histogram to refine against
            firstVaryList = []
            for item in varyList:
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                firstVaryList.append(item)
            newVaryList = firstVaryList
        else:
            newVaryList = []
            for item in varyList:
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                newVaryList.append(item)
        if newVaryList != firstVaryList and Controls['Copy2Next']:
            # variable lists are expected to match between sequential refinements when Copy2Next is on
            print '**** ERROR - variable list for this histogram does not match previous'
            print '     Copy of variables is not possible'
            print '\ncurrent histogram',histogram,'has',len(newVaryList),'variables'
            combined = list(set(firstVaryList+newVaryList))
            c = [var for var in combined if var not in newVaryList]
            p = [var for var in combined if var not in firstVaryList]
            line = 'Variables in previous but not in current: '
            if c:
                for var in c:
                    if len(line) > 100:
                        print line
                        line = '    '
                    line += var + ', '
            else:
                line += 'none'
            print line
            print '\nPrevious refinement has',len(firstVaryList),'variables'
            line = 'Variables in current but not in previous: '
            if p:
                for var in p:
                    if len(line) > 100:
                        print line
                        line = '    '
                    line += var + ', '
            else:
                line += 'none'
            print line
            return False,line
        
        ifPrint = False
        print >>printFile,'\n Refinement results for histogram: v'+histogram
        print >>printFile,135*'-'
        try:
            IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histo,Phases,restraintDict,
                rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifPrint,printFile,dlg)
    
            print '  wR = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f, last delta chi = %.4f'%(
                Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']**2,Rvals['DelChi2'])
            # add the uncertainties into the esd dictionary (sigDict)
            sigDict = dict(zip(varyList,sig))
            # the uncertainties for dependent constrained parms into the esd dict
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
    
            # a dict with values & esds for dependent (constrained) parameters
            depParmDict = {i:(parmDict[i],sigDict[i]) for i in varyListStart
                           if i not in varyList}
            newCellDict = copy.deepcopy(G2stMth.GetNewCellParms(parmDict,varyList))
            newAtomDict = copy.deepcopy(G2stMth.ApplyXYZshifts(parmDict,varyList))
            histRefData = {
                'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                'varyListStart':varyListStart,
                'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,
                'newCellDict':newCellDict,'depParmDict':depParmDict,
                'constraintInfo':constraintInfo,
                'parmDict':parmDict}
            SeqResult[histogram] = histRefData
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
    #        G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,ifPrint,printFile)
            G2stIO.SetHistogramData(parmDict,sigDict,Histo,ifPrint,printFile)
            G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histo,Phases,rigidbodyDict,histRefData,makeBack)
            makeBack = False
            NewparmDict = {}
            # make dict of varied parameters in current histogram, renamed to
            # next histogram, for use in next refinement. 
            if Controls['Copy2Next'] and ihst < len(histNames)-1:
                hId = Histo[histogram]['hId'] # current histogram
                nexthId = Histograms[histNames[ihst+1]]['hId']
                for parm in set(list(varyList)+list(varyListStart)):
                    items = parm.split(':')
                    if len(items) < 3: continue
                    if str(hId) in items[1]:
                        items[1] = str(nexthId)
                        newparm = ':'.join(items)
                        NewparmDict[newparm] = parmDict[parm]
        except G2stMth.UserAbort:
            printFile.close()
            print ' ***** Refinement aborted *****'
            return False,' Refinement aborted by user'
    G2stIO.SetSeqResult(GPXfile,Histograms,SeqResult)
    printFile.close()
    print ' Sequential refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst'
    print ' ***** Sequential refinement successful *****'
    return True,'Success'

def RetDistAngle(DisAglCtls,DisAglData):
    '''Compute and return distances and angles

    :param dict DisAglCtls: contains distance/angle radii usually defined using
       :func:`GSASIIgrid.DisAglDialog`
    :param dict DisAglData: contains phase data: 
       Items 'OrigAtoms' and 'TargAtoms' contain the atoms to be used
       for distance/angle origins and atoms to be used as targets.
       Item 'SGData' has the space group information (see :ref:`Space Group object<SGData_table>`)

    :returns: AtomLabels,DistArray,AngArray where:

      **AtomLabels** is a dict of atom labels, keys are the atom number

      **DistArray** is a dict keyed by the origin atom number where the value is a list
      of distance entries. The value for each distance is a list containing:
      
        0) the target atom number (int);
        1) the unit cell offsets added to x,y & z (tuple of int values)
        2) the symmetry operator number (which may be modified to indicate centering and center of symmetry)
        3) an interatomic distance in A (float)
        4) an uncertainty on the distance in A or 0.0 (float)

      **AngArray** is a dict keyed by the origin (central) atom number where
      the value is a list of
      angle entries. The value for each angle entry consists of three values:

        0) a distance item reference for one neighbor (int)
        1) a distance item reference for a second neighbor (int)
        2) a angle, uncertainty pair; the s.u. may be zero (tuple of two floats)

      The AngArray distance reference items refer directly to the index of the items in the
      DistArray item for the list of distances for the central atom. 
    '''
    import numpy.ma as ma
    
    SGData = DisAglData['SGData']
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

    Factor = DisAglCtls['Factors']
    Radii = dict(zip(DisAglCtls['AtomTypes'],zip(DisAglCtls['BondRadii'],DisAglCtls['AngleRadii'])))
    indices = (-1,0,1)
    Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    origAtoms = DisAglData['OrigAtoms']
    targAtoms = DisAglData['TargAtoms']
    AtomLabels = {}
    for Oatom in origAtoms:
        AtomLabels[Oatom[0]] = Oatom[1]
    for Oatom in targAtoms:
        AtomLabels[Oatom[0]] = Oatom[1]
    DistArray = {}
    AngArray = {}
    for Oatom in origAtoms:
        DistArray[Oatom[0]] = []
        AngArray[Oatom[0]] = []
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
                                tunit = (unit[0]+Tunit[0],unit[1]+Tunit[1],unit[2]+Tunit[2])
                                pdpx = G2mth.getDistDerv(Oatom[3:6],Tatom[3:6],Amat,unit,Top,SGData)
                                sig = 0.0
                                if len(Xvcov):
                                    sig = np.sqrt(np.inner(pdpx,np.inner(Xvcov,pdpx)))
                                Dist.append([Oatom[0],Tatom[0],tunit,Top,ma.getdata(dist[indb])[i],sig])
                                if (Dist[-1][-2]-AsumR) <= 0.:
                                    Vect.append(dx.T[indb][i]/Dist[-1][-2])
                                    VectA.append([OxyzNames,np.array(Oatom[3:6]),TxyzNames,np.array(Tatom[3:6]),unit,Top])
                                else:
                                    Vect.append([0.,0.,0.])
                                    VectA.append([])
        for D in Dist:
            DistArray[Oatom[0]].append(D[1:])
        Vect = np.array(Vect)
        angles = np.zeros((len(Vect),len(Vect)))
        angsig = np.zeros((len(Vect),len(Vect)))
        for i,veca in enumerate(Vect):
            if np.any(veca):
                for j,vecb in enumerate(Vect):
                    if np.any(vecb):
                        angles[i][j],angsig[i][j] = G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)
                        if i <= j: continue
                        AngArray[Oatom[0]].append((i,j,
                            G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)))
    return AtomLabels,DistArray,AngArray

def PrintDistAngle(DisAglCtls,DisAglData,out=sys.stdout):
    '''Print distances and angles

    :param dict DisAglCtls: contains distance/angle radii usually defined using
       :func:`GSASIIgrid.DisAglDialog`
    :param dict DisAglData: contains phase data: 
       Items 'OrigAtoms' and 'TargAtoms' contain the atoms to be used
       for distance/angle origins and atoms to be used as targets.
       Item 'SGData' has the space group information (see :ref:`Space Group object<SGData_table>`)
    :param file out: file object for output. Defaults to sys.stdout.    
    '''
    import numpy.ma as ma
    def MyPrint(s):
        out.write(s+'\n')
        # print(s,file=out) # use in Python 3
    
    def ShowBanner(name):
        MyPrint(80*'*')
        MyPrint('   Interatomic Distances and Angles for phase '+name)
        MyPrint((80*'*')+'\n')

    ShowBanner(DisAglCtls['Name'])
    SGData = DisAglData['SGData']
    SGtext,SGtable = G2spc.SGPrint(SGData)
    for line in SGtext: MyPrint(line)
    if len(SGtable):
        for i,item in enumerate(SGtable[::2]):
            line = ' %s %s'%(item.ljust(30),SGtable[2*i+1].ljust(30))
            MyPrint(line)   
    else:
        MyPrint(' ( 1)    %s'%(SGtable[0])) 
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
        MyPrint(line)
    else: 
        MyPrint('\n Unit cell: a = '+('%.5f'%Cell[0])+' b = '+('%.5f'%Cell[1])+' c = '+('%.5f'%Cell[2])+
            ' alpha = '+('%.3f'%Cell[3])+' beta = '+('%.3f'%Cell[4])+' gamma = '+
            ('%.3f'%Cell[5])+' volume = '+('%.3f'%Cell[6]))

    AtomLabels,DistArray,AngArray = RetDistAngle(DisAglCtls,DisAglData)
    origAtoms = DisAglData['OrigAtoms']
    targAtoms = DisAglData['TargAtoms']
    for Oatom in origAtoms:
        i = Oatom[0]
        Dist = DistArray[i]
        nDist = len(Dist)
        angles = np.zeros((nDist,nDist))
        angsig = np.zeros((nDist,nDist))
        for k,j,tup in AngArray[i]:
            angles[k][j],angsig[k][j] = angles[j][k],angsig[j][k] = tup
        line = ''
        for i,x in enumerate(Oatom[3:6]):
            line += ('%12.5f'%x).rstrip('0')
        MyPrint('\n Distances & angles for '+Oatom[1]+' at '+line.rstrip())
        MyPrint(80*'*')
        line = ''
        for dist in Dist[:-1]:
            line += '%12s'%(AtomLabels[dist[0]].center(12))
        MyPrint('  To       cell +(sym. op.)      dist.  '+line.rstrip())
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
            if dist[4]:            #sig exists!
                val = G2mth.ValEsd(dist[3],dist[4])
            else:
                val = '%8.4f'%(dist[3])
            tunit = '[%2d%2d%2d]'% dist[1]
            MyPrint(('  %8s%10s+(%4d) %12s'%(AtomLabels[dist[0]].ljust(8),tunit.ljust(10),dist[2],val.center(12)))+line.rstrip())

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
