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
from __future__ import division, print_function
import platform
import sys
import os.path as ospath
import time
import math
import copy
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import numpy as np
import numpy.linalg as nl
import scipy.optimize as so
import GSASIIpath
GSASIIpath.SetBinaryPath()
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIstrIO as G2stIO
import GSASIIstrMath as G2stMth
import GSASIIobj as G2obj
import GSASIIfiles as G2fil

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

ateln2 = 8.0*math.log(2.0)
DEBUG = True

def RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
    calcControls,pawleyLookup,ifSeq,printFile,dlg,refPlotUpdate=None):
    '''Core optimization routines, shared between SeqRefine and Refine

    :returns: 5-tuple of ifOk (bool), Rvals (dict), result, covMatrix, sig
    '''
    #patch (added Oct 2020) convert variable names for parm limits to G2VarObj
    import GSASIIscriptable as G2sc
    G2sc.patchControls(Controls)
    # end patch
#    print 'current',varyList
#    for item in parmDict: print item,parmDict[item] ######### show dict just before refinement
    G2mv.Map2Dict(parmDict,varyList)
    ifPrint = True
    if ifSeq:
        ifPrint = False
    Rvals = {}
    chisq0 = None
    while True:
        begin = time.time()
        values =  np.array(G2stMth.Dict2Values(parmDict, varyList))
        if np.any(np.isnan(values)):
            raise G2obj.G2Exception('ERROR - nan found in LS parameters - use Calculate/View LS parms to locate')
        # test code to compute GOF and save for external repeat
        #args = ([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg)
        #print '*** before fit chi**2',np.sum(G2stMth.errRefine(values,*args)**2)
        #fl = open('beforeFit.cpickle','wb')
        #cPickle.dump(values,fl,1)
        #cPickle.dump(args[:-1],fl,1)
        #fl.close()
        Ftol = Controls['min dM/M']
        Xtol = Controls['SVDtol']
        Factor = Controls['shift factor']
        if 'Jacobian' in Controls['deriv type']:
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
            if refPlotUpdate is not None: refPlotUpdate(Histograms)   # update plot after completion
        elif 'analytic Hessian' in Controls['deriv type']:
            Lamda = Controls.get('Marquardt',-3)
            maxCyc = Controls['max cyc']
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,xtol=Xtol,maxcyc=maxCyc,Print=ifPrint,lamda=Lamda,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg),
                refPlotUpdate=refPlotUpdate)
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
            Controls['Marquardt'] = -3  #reset to default
            if 'chisq0' in result[2] and chisq0 is None:
                chisq0 = result[2]['chisq0']
        elif 'Hessian SVD' in Controls['deriv type']:
            maxCyc = Controls['max cyc']
            result = G2mth.HessianSVD(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,xtol=Xtol,maxcyc=maxCyc,Print=ifPrint,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg),
                refPlotUpdate=refPlotUpdate)
            if result[1] is None:
                IfOK = False
                covMatrix = []
                sig = len(varyList)*[None,]
                break
            ncyc = result[2]['num cyc']+1
            if 'chisq0' in result[2] and chisq0 is None:
                chisq0 = result[2]['chisq0']
        else:           #'numeric'
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = 1
            if len(varyList):
                ncyc = int(result[2]['nfev']/len(varyList))
            if refPlotUpdate is not None: refPlotUpdate(Histograms)   # update plot
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['SVD0'] = result[2].get('SVD0',0)
        Rvals['converged'] = result[2].get('Converged')
        Rvals['DelChi2'] = result[2].get('DelChi2',-1.)
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        G2stMth.Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict,varyList)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs']-len(varyList)))
        printFile.write(' Number of function calls: %d No. of observations: %d No. of parameters: %d User rejected: %d Sp. gp. extinct: %d\n'%  \
            (result[2]['nfev'],Histograms['Nobs'],len(varyList),Histograms['Nrej'],Histograms['Next']))
        if ncyc:
            printFile.write(' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles\n'%(runtime,runtime/ncyc,ncyc))
        printFile.write(' wR = %7.2f%%, chi**2 = %12.6g, GOF = %6.2f\n'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']))
        sig = len(varyList)*[None,]
        if 'None' in str(type(result[1])) and ifSeq:    #this bails out of a sequential refinement on singular matrix
            IfOK = False
            covMatrix = []
            G2fil.G2Print ('Warning: **** Refinement failed - singular matrix ****')
            if 'Hessian' in Controls['deriv type']:
                num = len(varyList)-1
                for i,val in enumerate(np.flipud(result[2]['psing'])):
                    if val:
                        G2fil.G2Print('Bad parameter: '+varyList[num-i],mode='warn')
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        G2fil.G2Print('Bad parameter: '+varyList[ipvt-1],mode='warn')
            break
        IfOK = True
        try:
            covMatrix = result[1]*Rvals['GOF']**2
            sig = np.sqrt(np.diag(covMatrix))
            Lastshft = result[2].get('Xvec',None)
            if Lastshft is None:
                Rvals['Max shft/sig'] = 0.0
            else:
                Rvals['Max shft/sig'] = np.max(Lastshft/sig)
            if np.any(np.isnan(sig)) or not sig.shape:
                G2fil.G2Print ('*** Least squares aborted - some invalid esds possible ***',mode='error')
#            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except TypeError:          #result[1] is None on singular matrix or LinAlgError
            IfOK = False
            if not len(varyList):
                covMatrix = []
                break
            G2fil.G2Print ('**** Refinement failed - singular matrix ****',mode='error')
            if 'Hessian' in Controls['deriv type']:
                if result[1] is None:
                    IfOK = False
                    covMatrix = []
                    sig = len(varyList)*[None,]
                    break
                num = len(varyList)-1
                for i,val in enumerate(np.flipud(result[2]['psing'])):
                    if val:
                        G2fil.G2Print ('Removing parameter: '+varyList[num-i])
                        del(varyList[num-i])
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        G2fil.G2Print ('Removing parameter: '+varyList[ipvt-1])
                        del(varyList[ipvt-1])
                        break
    if IfOK:
        G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
    if chisq0 is not None:
        Rvals['GOF0'] = np.sqrt(chisq0/(Histograms['Nobs']-len(varyList)))
    return IfOK,Rvals,result,covMatrix,sig

def Refine(GPXfile,dlg=None,makeBack=True,refPlotUpdate=None):
    'Global refinement -- refines to minimize against all histograms'
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
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
        G2fil.G2Print (' *** ERROR - you have no phases to refine! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,'No phases'
    if not Histograms:
        G2fil.G2Print (' *** ERROR - you have no data to refine with! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,
         maxSSwave) = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,pFile=printFile)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    calcControls['MFtables'] = MFtables
    calcControls['maxSSwave'] = maxSSwave
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,pFile=printFile)
    TwConstr,TwFixed = G2stIO.makeTwinFrConstr(Phases,Histograms,hapVary)
    constrDict += TwConstr
    fixedList += TwFixed
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
    msg = G2mv.EvaluateMultipliers(constrDict,parmDict)
    if msg:
        return False,'Unable to interpret multiplier(s): '+msg
    try:
        G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict)
        #print(G2mv.VarRemapShow(varyList))
        #print('DependentVars',G2mv.GetDependentVars())
        #print('IndependentVars',G2mv.GetIndependentVars())
    except G2mv.ConstraintException:
        G2fil.G2Print (' *** ERROR - your constraints are internally inconsistent ***')
        #errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
        #print 'Errors',errmsg
        #if warnmsg: print 'Warnings',warnmsg
        return False,' Constraint error'
#    print G2mv.VarRemapShow(varyList)

    # remove frozen vars from refinement
    if 'parmFrozen' not in Controls:
        Controls['parmFrozen'] = {}
    if 'FrozenList' not in Controls['parmFrozen']: 
        Controls['parmFrozen']['FrozenList'] = []
    parmFrozenList = Controls['parmFrozen']['FrozenList']
    frozenList = [i for i in varyList if i in parmFrozenList]
    if len(frozenList) != 0: 
        varyList = [i for i in varyList if i not in parmFrozenList]
        G2fil.G2Print(
            'Frozen refined variables (due to exceeding limits)\n\t:{}'
            .format(frozenList))
        
    ifSeq = False
    printFile.write('\n Refinement results:\n')
    printFile.write(135*'-'+'\n')
    try:
        covData = {}
        IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histograms,Phases,restraintDict,
            rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg,
            refPlotUpdate=refPlotUpdate)
        if IfOK:
            sigDict = dict(zip(varyList,sig))
            newCellDict = G2stMth.GetNewCellParms(parmDict,varyList)
            newAtomDict = G2stMth.ApplyXYZshifts(parmDict,varyList)
            covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                       'varyListStart':varyListStart,
                       'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,
                       'newCellDict':newCellDict,'freshCOV':True}
            # add the uncertainties into the esd dictionary (sigDict)
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
            # check for variables outside their allowed range, reset and freeze them
            frozen = dropOOBvars(varyList,parmDict,sigDict,Controls,parmFrozenList)
            G2mv.PrintIndependentVars(parmDict,varyList,sigDict,pFile=printFile)
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
            G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
            G2stIO.SetPhaseData(parmDict,sigDict,Phases,rbIds,covData,restraintDict,printFile)
            G2stIO.PrintISOmodes(printFile,Phases,parmDict,sigDict)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,calcControls['FFtables'],pFile=printFile)
            G2stIO.SetHistogramData(parmDict,sigDict,Histograms,calcControls['FFtables'],pFile=printFile)
            if len(frozen):
                msg = ('Warning: {} variable(s) refined outside limits and were frozen ({} total frozen)'
                    .format(len(frozen),len(parmFrozenList))
                    )
                G2fil.G2Print(msg)
                Rvals['msg'] = msg
            elif len(parmFrozenList):
                msg = ('Note: a total of {} variable(s) are frozen due to refining outside limits'
                    .format(len(parmFrozenList))
                    )
                G2fil.G2Print('Note: ',msg)
                Rvals['msg'] = msg
            G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,rigidbodyDict,covData,parmFrozenList,makeBack)
            printFile.close()
            G2fil.G2Print (' Refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst')
            G2fil.G2Print (' ***** Refinement successful *****')
        else:
            G2fil.G2Print ('****ERROR - Refinement failed')
            raise G2obj.G2Exception('****ERROR - Refinement failed')
    except G2obj.G2RefineCancel as Msg:
        printFile.close()
        G2fil.G2Print (' ***** Refinement stopped *****')
        return False,Msg.msg
    except G2obj.G2Exception as Msg:  # cell metric error, others?
        printFile.close()
        G2fil.G2Print (' ***** Refinement error *****')
        return False,Msg.msg

#for testing purposes, create a file for testderiv
    if GSASIIpath.GetConfigValue('debug'):   # and IfOK:
#needs: values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup
        fl = open(ospath.splitext(GPXfile)[0]+'.testDeriv','wb')
        cPickle.dump(result[0],fl,1)
        cPickle.dump([Histograms,Phases,restraintDict,rigidbodyDict],fl,1)
        cPickle.dump([constrDict,fixedList,G2mv.GetDependentVars()],fl,1)
        cPickle.dump(parmDict,fl,1)
        cPickle.dump(varyList,fl,1)
        cPickle.dump(calcControls,fl,1)
        cPickle.dump(pawleyLookup,fl,1)
        fl.close()
    if dlg:
        return True,Rvals

def phaseCheck(phaseVary,Phases,histogram):
    '''
    Removes unused parameters from phase varylist if phase not in histogram
    #TODO - implement "Fix FXU" for seq refinement here - done?
    '''
    NewVary = []
    for phase in Phases:
        if histogram not in Phases[phase]['Histograms']: continue
        if Phases[phase]['Histograms'][histogram]['Use']:
            pId = Phases[phase]['pId']
            newVary = [item for item in phaseVary if item.split(':')[0] == str(pId)]
            FixVals = Phases[phase]['Histograms'][histogram].get('Fix FXU',' ')
            if 'F' in FixVals:
                newVary = [item for item in newVary if not 'Afrac' in item]
            if 'X' in FixVals:
                newVary = [item for item in newVary if not 'dA' in item]
            if 'U' in FixVals:
                newVary = [item for item in newVary if not 'AU' in item]
            if 'M' in FixVals:
                newVary = [item for item in newVary if not 'AM' in item]
            NewVary += newVary
    return NewVary

def SeqRefine(GPXfile,dlg,refPlotUpdate=None):
    '''Perform a sequential refinement -- cycles through all selected histgrams,
    one at a time
    '''
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    msgs = {}
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    G2fil.G2Print ('Starting Sequential Refinement')
    G2stIO.ShowBanner(printFile)
    Controls = G2stIO.GetControls(GPXfile)
    preFrozenCount = 0
    for h in Controls['parmFrozen']:
        if h == 'FrozenList': continue
        preFrozenCount += len(Controls['parmFrozen'][h])    
    G2stIO.ShowControls(Controls,printFile,SeqRef=True,
                            preFrozenCount=preFrozenCount)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        G2fil.G2Print (' *** ERROR - you have no phases to refine! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,'No phases'
    if not Histograms:
        G2fil.G2Print (' *** ERROR - you have no data to refine with! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    G2mv.InitVars()
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,
         maxSSwave) = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,
                                    Print=False,pFile=printFile,seqRef=True)
    for item in phaseVary:
        if '::A0' in item:
            G2fil.G2Print ('**** WARNING - lattice parameters should not be refined in a sequential refinement ****')
            G2fil.G2Print ('****           instead use the Dij parameters for each powder histogram            ****')
            return False,'Lattice parameter refinement error - see console message'
        if '::C(' in item:
            G2fil.G2Print ('**** WARNING - phase texture parameters should not be refined in a sequential refinement ****')
            G2fil.G2Print ('****           instead use the C(L,N) parameters for each powder histogram               ****')
            return False,'Phase texture refinement error - see console message'
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else: # patch from before Controls['Seq Data'] was implemented? 
        histNames = G2stIO.GetHistogramNames(GPXfile,['PWDR',])
    if Controls.get('Reverse Seq'):
        histNames.reverse()
    SeqResult = G2stIO.GetSeqResult(GPXfile)
#    SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
    Histo = {}
    NewparmDict = {}
    G2stIO.SetupSeqSavePhases(GPXfile)
    msgs['steepestNum'] = 0
    msgs['maxshift/sigma'] = []
    for ihst,histogram in enumerate(histNames):
        if GSASIIpath.GetConfigValue('Show_timing'): t1 = time.time()
        G2fil.G2Print('\nRefining with '+str(histogram))
        G2mv.InitVars()
        (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,
             FFtables,BLtables,MFtables,maxSSwave) = G2stIO.GetPhaseData(
                 Phases,restraintDict,rbIds,
                 Print=False,pFile=printFile,seqRef=True)
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['atomIndx'] = atomIndx
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['BLtables'] = BLtables
        calcControls['MFtables'] = MFtables
        calcControls['maxSSwave'] = maxSSwave
        if histogram not in Histograms:
            G2fil.G2Print("Error: not found!")
            continue
    #TODO - implement "Fix FXU" for seq refinement here - done?
        hId = Histograms[histogram]['hId']
        redphaseVary = phaseCheck(phaseVary,Phases,histogram)
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histo,Print=False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = rbVary+redphaseVary+hapVary+histVary
#        if not ihst:
            # save the initial vary list, but without histogram numbers on parameters
        saveVaryList = varyList[:]
        for i,item in enumerate(saveVaryList):
            items = item.split(':')
            if items[1]:
                items[1] = ''
            item = ':'.join(items)
            saveVaryList[i] = item
        if not ihst:
            SeqResult['varyList'] = saveVaryList
        else:
            SeqResult['varyList'] = list(set(SeqResult['varyList']+saveVaryList))
        parmDict = {}
        parmDict.update(rbDict)
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        if Controls['Copy2Next']:   # update with parms from last histogram
            #parmDict.update(NewparmDict) # don't use in case extra entries would cause a problem
            for parm in NewparmDict:
                if parm in parmDict:
                    parmDict[parm] = NewparmDict[parm]
        elif histogram in SeqResult:  # update phase from last seq ref
            NewparmDict = SeqResult[histogram].get('parmDict',{})
            for parm in NewparmDict:
                if '::' in parm and parm in parmDict:
                    parmDict[parm] = NewparmDict[parm]
            
        G2stIO.GetFprime(calcControls,Histo)
        # do constraint processing
        #reload(G2mv) # debug
        constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
        varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed
        msg = G2mv.EvaluateMultipliers(constrDict,parmDict)
        if msg:
            return False,'Unable to interpret multiplier(s): '+msg
        try:
            groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict,SeqHist=hId)
#            if GSASIIpath.GetConfigValue('debug'): print("DBG_"+
#                G2mv.VarRemapShow(varyList,True))
#            print('DependentVars',G2mv.GetDependentVars())
#            print('IndependentVars',G2mv.GetIndependentVars())
            constraintInfo = (groups,parmlist,constrDict,fixedList,ihst)
        except G2mv.ConstraintException:
            G2fil.G2Print (' *** ERROR - your constraints are internally inconsistent ***')
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
            #print '**** ERROR - variable list for this histogram does not match previous'
            #print '     Copy of variables is not possible'
            #print '\ncurrent histogram',histogram,'has',len(newVaryList),'variables'
            combined = list(set(firstVaryList+newVaryList))
            c = [var for var in combined if var not in newVaryList]
            p = [var for var in combined if var not in firstVaryList]
            G2fil.G2Print('*** Variables change ***')
            for typ,vars in [('Removed',c),('Added',p)]:
                line = '  '+typ+': '
                if vars:
                    for var in vars:
                        if len(line) > 70:
                            G2fil.G2Print(line)
                            line = '    '
                        line += var + ', '
                else:
                        line += 'none, '
                G2fil.G2Print(line[:-2])
            firstVaryList = newVaryList

        ifSeq = True
        printFile.write('\n Refinement results for histogram id {}: {}\n'
                            .format(hId,histogram))
        printFile.write(135*'-'+'\n')
        # remove frozen vars
        if 'parmFrozen' not in Controls:
            Controls['parmFrozen'] = {}
        if histogram not in Controls['parmFrozen']: 
            Controls['parmFrozen'][histogram] = []
        parmFrozenList = Controls['parmFrozen'][histogram]
        frozenList = [i for i in varyList if i in parmFrozenList]
        if len(frozenList) != 0: 
           varyList = [i for i in varyList if i not in parmFrozenList]
           s = ''
           for a in frozenList:
               if s:
                   s+= ', '
               s += a
           printFile.write(
               ' The following refined variables have previously been frozen due to exceeding limits:\n\t{}\n'
               .format(s))
        try:
            IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histo,Phases,restraintDict,
                rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg,
                refPlotUpdate=refPlotUpdate)
            G2fil.G2Print ('  wR = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f, last delta chi = %.4f, last shft/sig = %.4f'%(
                Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']**2,Rvals['DelChi2'],Rvals.get('Max shft/sig',np.nan)))
            if Rvals.get('lamMax',0) >= 10.:
                msgs['steepestNum'] += 1
            if Rvals.get('Max shft/sig'):
                msgs['maxshift/sigma'].append(Rvals['Max shft/sig'])
            # add the uncertainties into the esd dictionary (sigDict)
            if not IfOK:
                G2fil.G2Print('***** Sequential refinement failed at histogram '+histogram,mode='warn')
                break
            sigDict = dict(zip(varyList,sig))
            # the uncertainties for dependent constrained parms into the esd dict
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
            # check for variables outside their allowed range, reset and freeze them
            frozen = dropOOBvars(varyList,parmDict,sigDict,Controls,parmFrozenList)
            msg = None
            if len(frozen) > 0:
               msg = ('Hist {}: {} variables were outside limits and were frozen (now {} frozen total)'
                   .format(ihst,len(frozen),len(parmFrozenList)))
               G2fil.G2Print(msg)
               msg = (' {} variables were outside limits and were frozen (now {} frozen total)'
                   .format(len(frozen),len(parmFrozenList)))
               for p in frozen:
                   if p not in varyList:
                       print('Frozen Warning: {} not in varyList. This should not happen!'.format(p))
                       continue
                   i = varyList.index(p)
                   result[0][i] = parmDict[p]
                   sig[i] = -0.1
            # a dict with values & esds for dependent (constrained) parameters - avoid extraneous holds
            depParmDict = {i:(parmDict[i],sigDict[i]) for i in varyListStart if i in sigDict and i not in varyList}
            newCellDict = copy.deepcopy(G2stMth.GetNewCellParms(parmDict,varyList))
            newAtomDict = copy.deepcopy(G2stMth.ApplyXYZshifts(parmDict,varyList))
            histRefData = {
                'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                'varyListStart':varyListStart,
                'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,
                'newCellDict':newCellDict,'depParmDict':depParmDict,
                'constraintInfo':constraintInfo,
                'parmDict':parmDict,
                }
            SeqResult[histogram] = histRefData
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
#            G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,None,ifPrint,printFile)
            G2stIO.SetHistogramData(parmDict,sigDict,Histo,None,ifPrint,printFile,seq=True)
            G2stIO.SaveUpdatedHistogramsAndPhases(GPXfile,Histo,Phases,rigidbodyDict,histRefData,Controls['parmFrozen'])
            if msg: 
                printFile.write(msg+'\n')
            NewparmDict = {}
            # make dict of varied parameters in current histogram, renamed to
            # next histogram, for use in next refinement.
            if Controls['Copy2Next'] and ihst < len(histNames)-1:
                hId = Histo[histogram]['hId'] # current histogram
                nexthId = Histograms[histNames[ihst+1]]['hId']
                for parm in set(list(varyList)+list(varyListStart)):
                    items = parm.split(':')
                    if len(items) < 3: 
                        continue
                    if str(hId) in items[1]:
                        items[1] = str(nexthId)
                        newparm = ':'.join(items)
                        NewparmDict[newparm] = parmDict[parm]
                    else:
                        if items[2].startswith('dA'): parm = parm.replace(':dA',':A') 
                        NewparmDict[parm] = parmDict[parm]
                    
        except G2obj.G2RefineCancel as Msg:
            printFile.close()
            G2fil.G2Print (' ***** Refinement stopped *****')
            return False,Msg.msg
        except G2obj.G2Exception as Msg:  # cell metric error, others?
            printFile.close()
            G2fil.G2Print (' ***** Refinement error *****')
            return False,Msg.msg
        if GSASIIpath.GetConfigValue('Show_timing'):
            t2 = time.time()
            G2fil.G2Print("Fit step time {:.2f} sec.".format(t2-t1))
            t1 = t2
    SeqResult['histNames'] = [itm for itm in G2stIO.GetHistogramNames(GPXfile,['PWDR',]) if itm in SeqResult.keys()]
    try:
        G2stIO.SetSeqResult(GPXfile,Histograms,SeqResult)
    except Exception as msg:
        print('Error reading Sequential results')
        if GSASIIpath.GetConfigValue('debug'):
            import traceback
            print(traceback.format_exc())        
    postFrozenCount = 0
    for h in Controls['parmFrozen']:
        if h == 'FrozenList': continue
        postFrozenCount += len(Controls['parmFrozen'][h])
    if postFrozenCount:
        msgs['Frozen'] = 'Ending refinement with {} Frozen variables ({} added now)\n'.format(postFrozenCount,postFrozenCount-preFrozenCount)
        printFile.write('\n'+msgs['Frozen'])
    printFile.close()
    G2fil.G2Print (' Sequential refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst')
    G2fil.G2Print (' ***** Sequential refinement successful *****')
    return True,msgs

def dropOOBvars(varyList,parmDict,sigDict,Controls,parmFrozenList):
    '''Find variables in the parameters dict that are outside the ranges 
    (in parmMinDict and parmMaxDict) and set them to the limits values. 
    Add any such variables into the list of frozen variable 
    (parmFrozenList). Returns a list of newly frozen variables, if any.
    '''
    parmMinDict = Controls.get('parmMinDict',{})
    parmMaxDict = Controls.get('parmMaxDict',{})
    freeze = []
    if parmMinDict or parmMaxDict:
        for name in varyList:
            if name not in parmDict: continue
            n,val = G2obj.prmLookup(name,parmMinDict)
            if n is not None:
                if parmDict[name] < parmMinDict[n]:
                    parmDict[name] = parmMinDict[n]
                    sigDict[name] = 0.0
                    freeze.append(name)
                    continue
            n,val = G2obj.prmLookup(name,parmMaxDict)
            if n is not None:
                if parmDict[name] > parmMaxDict[n]:
                    parmDict[name] = parmMaxDict[n]
                    sigDict[name] = 0.0
                    freeze.append(name)
                    continue
        for v in freeze:
            if v not in parmFrozenList:
                parmFrozenList.append(v)
    return freeze

def RetDistAngle(DisAglCtls,DisAglData,dlg=None):
    '''Compute and return distances and angles

    :param dict DisAglCtls: contains distance/angle radii usually defined using
       :func:`GSASIIctrlGUI.DisAglDialog`
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
    if len(DisAglData.get('covData',{})):
        covData = DisAglData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(DisAglData['pId'])+'::'

    Factor = DisAglCtls['Factors']
    Radii = dict(zip(DisAglCtls['AtomTypes'],zip(DisAglCtls['BondRadii'],DisAglCtls['AngleRadii'])))
    indices = (-2,-1,0,1,2)
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
    for iO,Oatom in enumerate(origAtoms):
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
            if len(DisAglData.get('covData',{})):
                OxyzNames = [pfx+'dAx:%d'%(Oatom[0]),pfx+'dAy:%d'%(Oatom[0]),pfx+'dAz:%d'%(Oatom[0])]
                TxyzNames = [pfx+'dAx:%d'%(Tatom[0]),pfx+'dAy:%d'%(Tatom[0]),pfx+'dAz:%d'%(Tatom[0])]
                Xvcov = G2mth.getVCov(OxyzNames+TxyzNames,varyList,covMatrix)
            result = G2spc.GenAtom(Tatom[3:6],SGData,False,Move=False)
            BsumR = (Radii[Oatom[2]][0]+Radii[Tatom[2]][0])*Factor[0]
            AsumR = (Radii[Oatom[2]][1]+Radii[Tatom[2]][1])*Factor[1]
            for [Txyz,Top,Tunit,Spn] in result:
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
                                sig = 0.0
                                if len(Xvcov):
                                    pdpx = G2mth.getDistDerv(Oatom[3:6],Tatom[3:6],Amat,unit,Top,SGData)
                                    sig = np.sqrt(np.inner(pdpx,np.inner(pdpx,Xvcov)))
                                Dist.append([Oatom[0],Tatom[0],tunit,Top,ma.getdata(dist[indb])[i],sig])
                                if (Dist[-1][-2]-AsumR) <= 0.:
                                    Vect.append(dx.T[indb][i]/Dist[-1][-2])
                                    VectA.append([OxyzNames,np.array(Oatom[3:6]),TxyzNames,np.array(Tatom[3:6]),unit,Top])
                                else:
                                    Vect.append([0.,0.,0.])
                                    VectA.append([])
        if dlg is not None:
            dlg.Update(iO,newmsg='Atoms done=%d'%(iO))
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
       :func:`GSASIIctrlGUI.DisAglDialog`
    :param dict DisAglData: contains phase data:
       Items 'OrigAtoms' and 'TargAtoms' contain the atoms to be used
       for distance/angle origins and atoms to be used as targets.
       Item 'SGData' has the space group information (see :ref:`Space Group object<SGData_table>`)
    :param file out: file object for output. Defaults to sys.stdout.
    '''
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
    if len(SGtable) > 1:
        for i,item in enumerate(SGtable[::2]):
            if 2*i+1 == len(SGtable):
                line = ' %s'%(item.ljust(30))
            else:
                line = ' %s %s'%(item.ljust(30),SGtable[2*i+1].ljust(30))
            MyPrint(line)
    else:
        MyPrint(' ( 1)    %s'%(SGtable[0])) #triclinic case
    Cell = DisAglData['Cell']

    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if len(DisAglData.get('covData',{})):
        covData = DisAglData['covData']
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
            ('%.3f'%Cell[5])+' Volume = '+('%.3f'%Cell[6]))

    AtomLabels,DistArray,AngArray = RetDistAngle(DisAglCtls,DisAglData)
    origAtoms = DisAglData['OrigAtoms']
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
        G2fil.G2Print (' Torsion angle for %s atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Tors,sig)))
        G2fil.G2Print (' NB: Atom sequence determined by selection order')
        return      # done with torsion
    elif DATData['Natoms'] == 3:  #angle
        Ang,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        G2fil.G2Print (' Angle in %s for atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Ang,sig)))
        G2fil.G2Print (' NB: Atom sequence determined by selection order')
    else:   #2 atoms - distance
        Dist,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        G2fil.G2Print (' Distance in %s for atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Dist,sig)))

def BestPlane(PlaneData):
    'Needs a doc string'

    def ShowBanner(name):
        G2fil.G2Print (80*'*')
        G2fil.G2Print ('   Best plane result for phase '+name)
        G2fil.G2Print (80*'*','\n')

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
    G2fil.G2Print (' Selected atoms centered at %10.5f %10.5f %10.5f'%(sumXYZ[0],sumXYZ[1],sumXYZ[2]))
    Evec,Emat = nl.eig(Zmat)
    Evec = np.sqrt(Evec)/(Natoms-3)
    Order = np.argsort(Evec)
    XYZ = np.inner(XYZ,Emat.T).T
    XYZ = np.array([XYZ[Order[2]],XYZ[Order[1]],XYZ[Order[0]]]).T
    G2fil.G2Print (' Atoms in Cartesian best plane coordinates:')
    G2fil.G2Print (' Name         X         Y         Z')
    for i,xyz in enumerate(XYZ):
        G2fil.G2Print (' %6s%10.3f%10.3f%10.3f'%(Atoms[i][1].ljust(6),xyz[0],xyz[1],xyz[2]))
    G2fil.G2Print ('\n Best plane RMS X =%8.3f, Y =%8.3f, Z =%8.3f'%(Evec[Order[2]],Evec[Order[1]],Evec[Order[0]]))

def main():
    'Called to run a refinement when this module is executed '
    starttime = time.time()
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            G2fil.G2Print ('ERROR - '+GPXfile+" doesn't exist!")
            exit()
    else:
        G2fil.G2Print ('ERROR - missing filename')
        exit()
    # TODO: figure out if this is a sequential refinement and call SeqRefine(GPXfile,None)
    Refine(GPXfile,None)
    G2fil.G2Print("Done. Execution time {:.2f} sec.".format(time.time()-starttime))

if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    main()
