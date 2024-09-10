# -*- coding: utf-8 -*-
'''
:mod:`GSASIIstrMain` routines, used for refinement are found below.
'''
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
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIstrIO as G2stIO
import GSASIIstrMath as G2stMth
import GSASIIobj as G2obj
import GSASIIfiles as G2fil
import GSASIIElem as G2elem
import GSASIIscriptable as G2sc
import atmdata

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

ateln2 = 8.0*math.log(2.0)
DEBUG = True
#PhFrExtPOSig = None

def ReportProblems(result,Rvals,varyList):
    '''Create a message based results from the refinement
    '''
    #report on SVD 0's and highly correlated variables
    msg = ''
    # process singular variables; all vars go to console, first 10 to
    # dialog window
    psing = result[2].get('psing',[])
    if len(psing):
        if msg: msg += '\n'
        m = 'Error: {} Parameter(s) dropped:'.format(len(psing))
        msg += m
        G2fil.G2Print(m, mode='warn')
        m = ''
        for i,val in enumerate(psing):
            if i == 0:
                msg += '\n{}'.format(varyList[val])
                m = '  {}'.format(varyList[val])
            else:
                if len(m) > 70:
                    G2fil.G2Print(m, mode='warn')
                    m = '  '
                else:
                    m += ', '
                m += '{}'.format(varyList[val])
                if i == 10:
                    msg += ', {}... see console for full list'.format(varyList[val])
                elif i > 10:
                    pass
                else:
                    msg += ', {}'.format(varyList[val])
        if m: G2fil.G2Print(m, mode='warn')
    SVD0 = result[2].get('SVD0',0)
    if SVD0 == 1: 
        msg += 'Warning: Soft (SVD) singularity in the Hessian'
    elif SVD0 > 0: 
        msg += 'Warning: {} soft (SVD) Hessian singularities'.format(SVD0)
    SVDsing = result[2].get('SVDsing',[])
    if len(SVDsing):
        if msg: msg += '\n'
        m = 'SVD problem(s) likely from:'
        msg += m
        G2fil.G2Print(m, mode='warn')
        m = ''
        for i,val in enumerate(SVDsing):
            if i == 0:
                msg += '\n{}'.format(varyList[val])
                m = '  {}'.format(varyList[val])
            else:
                if len(m) > 70:
                    G2fil.G2Print(m, mode='warn')
                    m = '  '
                else:
                    m += ', '
                m += '{}'.format(varyList[val])
                if i == 10:
                    msg += ', {}... see console for full list'.format(varyList[val])
                elif i > 10:
                    pass
                else:
                    msg += ', {}'.format(varyList[val])
        if m: G2fil.G2Print(m, mode='warn')
    #report on highly correlated variables
    Hcorr = result[2].get('Hcorr',[])
    if len(Hcorr) > 0: 
        if msg: msg += '\n'
        m = 'Note highly correlated parameters:'
        G2fil.G2Print(m, mode='warn')
        msg += m
    elif SVD0 > 0:
        if msg: msg += '\n'
        m = 'Check covariance matrix for parameter correlation'
        G2fil.G2Print(m, mode='warn')
        msg += m
    for i,(v1,v2,corr) in enumerate(Hcorr):
        if corr > .95:
            stars = '**'
        else:
            stars = '   '
        m = ' {} {} and {} (@{:.2f}%)'.format(
            stars,varyList[v1],varyList[v2],100.*corr)
        G2fil.G2Print(m, mode='warn')
        if i == 5:
            msg += '\n' + m
            msg += '\n   ... check console for more'
        elif i < 5:
            msg += '\n' + m
    if msg:
        if 'msg' not in Rvals: Rvals['msg'] = ''
        Rvals['msg'] += msg

def IgnoredLatticePrms(Phases):
    ignore = []
    copydict = {}
    for p in Phases:
        pfx = str(Phases[p]['pId']) + '::'
        laue = Phases[p]['General']['SGData']['SGLaue']
        axis = Phases[p]['General']['SGData']['SGUniq']
        if laue in ['-1',]:
            pass
        elif laue in ['2/m',]:
            if axis == 'a':
                ignore += [pfx+'A4',pfx+'A5']
            elif axis == 'b':
                ignore += [pfx+'A3',pfx+'A5']
            else:
                ignore += [pfx+'A3',pfx+'A4']
        elif laue in ['mmm',]:
            ignore += [pfx+'A3',pfx+'A4',pfx+'A5']
        elif laue in ['4/m','4/mmm']:
            ignore += [pfx+'A1',pfx+'A3',pfx+'A4',pfx+'A5']
            copydict[pfx+'A0':[pfx+'A1']]
        elif laue in ['6/m','6/mmm','3m1', '31m', '3']:
            ignore += [pfx+'A1',pfx+'A3',pfx+'A4',pfx+'A5']
            copydict[pfx+'A0'] = [pfx+'A1',pfx+'A3']
        elif laue in ['3R', '3mR']:
            ignore += [pfx+'A1',pfx+'A2',pfx+'A4',pfx+'A5']
            copydict[pfx+'A0'] = [pfx+'A1',pfx+'A2']
            copydict[pfx+'A3'] = [pfx+'A4',pfx+'A5']
        elif laue in ['m3m','m3']:
            ignore += [pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
            copydict[pfx+'A0'] = [pfx+'A1',pfx+'A2']
    return ignore,copydict

def AllPrmDerivs(Controls,Histograms,Phases,restraintDict,rigidbodyDict,
    parmDict,varyList,calcControls,pawleyLookup,symHold,dlg=None):
    '''Computes the derivative of the fitting function (total Chi**2) with 
    respect to every parameter in the parameter dictionary (parmDict)
    by applying shift below the parameter value as well as above. 
    
    :returns: a dict with the derivatives keyed by variable number. 
      Derivatives are a list with three values: evaluated over
      v-d to v; v-d to v+d; v to v+d where v is the current value for the
      variable and d is a small delta value chosen for that variable type.
    '''
    import re
    rms = lambda y: np.sqrt(np.mean(y**2))
    G2mv.Map2Dict(parmDict,varyList)
    begin = time.time()
    seqList = Controls.get('Seq Data',[])
    hId = '*'
    if seqList:
        hId = str(Histograms[seqList[0]]['hId'])
#    values =  np.array(G2stMth.Dict2Values(parmDict, varyList))
#    if np.any(np.isnan(values)):
#        raise G2obj.G2Exception('ERROR - nan found in LS parameters - use Calculate/View LS parms to locate')
    latIgnoreLst,latCopyDict = IgnoredLatticePrms(Phases)
    HistoPhases = [Histograms,Phases,restraintDict,rigidbodyDict]
    origDiffs = G2stMth.errRefine([],HistoPhases,parmDict,[],calcControls,pawleyLookup,None)
    chiStart = rms(origDiffs)
    origParms = copy.deepcopy(parmDict)
    #print('after 1st calc',time.time()-begin)
    derivCalcs = {}
    if dlg: dlg.SetRange(len(origParms))
    for i,prm in enumerate(origParms):
        if dlg:
            if not dlg.Update(i)[0]:
                return None
        parmDict = copy.deepcopy(origParms)
        p,h,nam = prm.split(':')[:3]
        if 'UVmat' in nam:
            continue
        if hId != '*' and h != '' and h != hId: continue
        if (type(parmDict[prm]) is bool or type(parmDict[prm]) is str or
                 type(parmDict[prm]) is int): continue
        if type(parmDict[prm]) is not float and type(parmDict[prm]) is not np.float64: 
            print('*** unexpected type for ',prm,parmDict[prm],type(parmDict[prm]))
            continue
        if prm in latIgnoreLst: continue # remove unvaried lattice params
        if re.match(r'\d:\d:D[012][012]',prm): continue   # don't need Dij terms
        if nam in ['Vol','Gonio. radius']: continue
        if nam.startswith('dA') and nam[2] in ['x','y','z']: continue
        delta = max(abs(parmDict[prm])*0.0001,1e-6)
        if nam in ['Shift','DisplaceX','DisplaceY',]:
            delta = 0.1
        elif nam.startswith('AUiso'):
            delta = 1e-5
        if nam[0] == 'A' and nam[1] in ['x','y','z']: 
            dprm = prm.replace('::A','::dA')
            if dprm in symHold: continue # held by symmetry
            delta = 1e-6
        if nam in ['A0','A1','A2','A3','A4','A5'] and 'PWDR' not in Histograms.keys():
            continue
        else:
            dprm = prm
        #print('***',prm,type(parmDict[prm]))
        #origVal = parmDict[dprm]
        parmDict[dprm] -= delta
        G2mv.Dict2Map(parmDict)
        if dprm in latCopyDict:         # apply contraints on lattice parameters 
            for i in latCopyDict:
                parmDict[i] = parmDict[dprm]
        #for i in parmDict:
        #    if origParms[i] != parmDict[i]: print('changed',i,origParms[i],parmDict[i])
        chiLow = rms(G2stMth.errRefine([],HistoPhases,parmDict,[],calcControls,pawleyLookup,None))
        parmDict[dprm] += 2*delta
        G2mv.Dict2Map(parmDict)
        if dprm in latCopyDict:         # apply contraints on lattice parameters 
            for i in latCopyDict:
                parmDict[i] = parmDict[dprm]
        #for i in parmDict:
        #    if origParms[i] != parmDict[i]: print('changed',i,origParms[i],parmDict[i])
        chiHigh = rms(G2stMth.errRefine([],HistoPhases,parmDict,[],calcControls,pawleyLookup,None))
        #print('===>',prm,parmDict[dprm],delta)
        #print(chiLow,chiStart,chiHigh)
        #print((chiLow-chiStart)/delta,0.5*(chiLow-chiHigh)/delta,(chiStart-chiHigh)/delta)
        derivCalcs[prm] = ((chiLow-chiStart)/delta,0.5*(chiLow-chiHigh)/delta,(chiStart-chiHigh)/delta)
    print('derivative computation time',time.time()-begin)
    return derivCalcs

def RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
    calcControls,pawleyLookup,ifSeq,printFile,dlg,refPlotUpdate=None):
    '''Core optimization routines, shared between SeqRefine and Refine

    :returns: 5-tuple of ifOk (bool), Rvals (dict), result, covMatrix, sig
    '''
    #patch (added Oct 2020) convert variable names for parm limits to G2VarObj
    G2sc.patchControls(Controls)
    # end patch
#    print 'current',varyList
#    for item in parmDict: print item,parmDict[item] ######### show dict just before refinement
    ifPrint = True
    if ifSeq:
        ifPrint = False
    Rvals = {}
    chisq0 = None
    Lastshft = None
    IfOK = True
    while True:
        G2mv.Map2Dict(parmDict,varyList)
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
            maxCyc = Controls.get('max cyc',1)
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
            result[2]['num cyc'] = ncyc
            if refPlotUpdate is not None: refPlotUpdate(Histograms)   # update plot after completion
        elif 'analytic Hessian' in Controls['deriv type']:
            Lamda = Controls.get('Marquardt',-3)
            maxCyc = Controls['max cyc']
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,xtol=Xtol,maxcyc=maxCyc,Print=ifPrint,lamda=Lamda,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg),
                refPlotUpdate=refPlotUpdate)
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
            if 'lastShifts' in result[2]:
                Rvals['lastShifts'] = dict(zip(varyList,result[2]['lastShifts']))
            if 'Ouch#4' in  result[2]:
                Rvals['Aborted'] = True
            if 'msg' in result[2]:
                Rvals['msg'] = result[2]['msg']
            Controls['Marquardt'] = -3  #reset to default
            if 'chisq0' in result[2] and chisq0 is None:
                chisq0 = result[2]['chisq0']
            if maxCyc == 0:
                covMatrix = []
                sig = len(varyList)*[None,]
            elif result[1] is None:
                IfOK = False
                covMatrix = []
                sig = len(varyList)*[None,]
                break
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
            maxCyc = Controls.get('max cyc',1)
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = 1
            result[2]['num cyc'] = ncyc
            if len(varyList):
                ncyc = int(result[2]['nfev']/len(varyList))
            if refPlotUpdate is not None: refPlotUpdate(Histograms)   # update plot
        #table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
        #for item in table: print(item,table[item])               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['SVD0'] = result[2].get('SVD0',0)
        Rvals['converged'] = result[2].get('Converged')
        Rvals['DelChi2'] = result[2].get('DelChi2',-1.)
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        G2stMth.Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Nvars'] = len(varyList)
        Rvals['RestraintSum'] = Histograms.get('RestraintSum',0.)
        Rvals['RestraintTerms'] = Histograms.get('RestraintTerms',0)
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs']+Rvals['RestraintTerms']-len(varyList)))
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
                # BHT -- I am not sure if this works correctly:
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
        if not len(varyList) or not maxCyc:
            covMatrix = []
            break
        try:
            covMatrix = result[1]*Rvals['GOF']**2
            sig = np.sqrt(np.diag(covMatrix))
            Lastshft = result[0]-values     #NOT last shift since values is starting set before current refinement
            #table = dict(zip(varyList,zip(values,result[0],Lastshft,Lastshft/sig)))
            #for item in table: print(item,table[item])               #useful debug
            Rvals['Max shft/sig'] = np.max(np.nan_to_num(Lastshft/sig))
            if np.any(np.isnan(sig)) or not sig.shape:
                G2fil.G2Print ('*** Least squares aborted - some invalid esds possible ***',mode='error')
            else:
                print('Maximum shift/esd = {:.3f} for all cycles'.format(Rvals['Max shft/sig']))
            # report on refinement issues. Result in Rvals['msg']
            ReportProblems(result,Rvals,varyList)
            break                   #refinement succeeded - finish up!
        except TypeError:
            # if we get here, no result[1] (covar matrix) was returned or other calc error: refinement failed
            IfOK = False
            if 'Hessian' in Controls['deriv type']:
                SVD0 = result[2].get('SVD0')
                if SVD0 == -1:
                    G2fil.G2Print ('**** Refinement failed - singular matrix ****',mode='error')
                elif SVD0 == -2:
                    G2fil.G2Print ('**** Refinement failed - other problem ****',mode='error')
                elif SVD0 > 0:
                    G2fil.G2Print ('**** Refinement failed with {} SVD singularities ****'.format(SVD0),mode='error')
                else:
                    G2fil.G2Print ('**** Refinement failed ****',mode='error')
                if result[1] is None:
                    IfOK = False
                    covMatrix = []
                    sig = len(varyList)*[None,]
                # report on highly correlated variables
                ReportProblems(result,Rvals,varyList)
                # process singular variables
                if dlg: break # refining interactively
            else:
                G2fil.G2Print ('**** Refinement failed - singular matrix ****',mode='error')
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        G2fil.G2Print ('Removing parameter: '+varyList[ipvt-1])
                        del(varyList[ipvt-1])
                        break
    if IfOK:
        if CheckLeBail(Phases):   # only needed for LeBail extraction
            G2stMth.errRefine([],[Histograms,Phases,restraintDict,rigidbodyDict],
                parmDict,[],calcControls,pawleyLookup,dlg)
        G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
    if chisq0 is not None:
        Rvals['GOF0'] = np.sqrt(chisq0/(Histograms['Nobs']-len(varyList)))
    return IfOK,Rvals,result,covMatrix,sig,Lastshft

def Refine(GPXfile,dlg=None,makeBack=True,refPlotUpdate=None,newLeBail=False,allDerivs=False):
    '''Global refinement -- refines to minimize against all histograms. 
    This can be called in one of three ways, from :meth:`GSASIIdataGUI.GSASII.OnRefine` in an 
    interactive refinement, where dlg will be a wx.ProgressDialog, or non-interactively from 
    :meth:`GSASIIscriptable.G2Project.refine` or from :func:`do_refine`, where dlg will be None.
    '''
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

    if allDerivs:
        printFile = open(ospath.splitext(GPXfile)[0]+'.junk','w')
    else:
        printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    G2stIO.ShowBanner(printFile)
    varyList = []
    parmDict = {}
    G2mv.InitVars()
    Controls = G2stIO.GetControls(GPXfile)
    Controls['newLeBail'] = newLeBail
    G2stIO.ShowControls(Controls,printFile)
    calcControls = {}
    calcControls.update(Controls)
    constrDict,fixedList = G2stIO.ReadConstraints(GPXfile)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        G2fil.G2Print (' *** ERROR - you have no phases to refine! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,{'msg':'No phases'}
    if not Histograms:
        G2fil.G2Print (' *** ERROR - you have no data to refine with! ***')
        G2fil.G2Print (' *** Refine aborted ***')
        return False,{'msg':'No data'}
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    symHold = None
    if allDerivs: #=============  develop partial derivative map
        symHold = []
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,
         maxSSwave) = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,pFile=printFile,symHold=symHold)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['EFtables'] = EFtables
    calcControls['ORBtables'] = ORBtables
    calcControls['BLtables'] = BLtables
    calcControls['MFtables'] = MFtables
    calcControls['maxSSwave'] = maxSSwave
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,Controls=calcControls,pFile=printFile)
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
    if allDerivs: #=============  develop partial derivative map
        varyListStart = varyList
        varyList = None
    if msg:
        return False,{'msg':'Unable to interpret multiplier(s): '+msg}
    try:
        errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict)
        G2mv.normParms(parmDict)
        G2mv.Map2Dict(parmDict,varyList)   # changes varyList
    except G2mv.ConstraintException:
        G2fil.G2Print (' *** ERROR - your constraints are internally inconsistent ***')
        return False,{'msg':' Constraint error'}

    # remove frozen vars from refinement
    if 'parmFrozen' not in Controls:
        Controls['parmFrozen'] = {}
    if 'FrozenList' not in Controls['parmFrozen']: 
        Controls['parmFrozen']['FrozenList'] = []
    if varyList is not None:
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
    Rvals = {}
    G2mv.Dict2Map(parmDict)  # impose constraints initially
    if allDerivs: #=============  develop partial derivative map
        derivDict = AllPrmDerivs(Controls, Histograms, Phases, restraintDict,
            rigidbodyDict, parmDict, varyList, calcControls,pawleyLookup,symHold,dlg)
        printFile.close() #closes the .junk file
        return derivDict,varyListStart
    try:
        covData = {}
        IfOK,Rvals,result,covMatrix,sig,Lastshft = RefineCore(Controls,Histograms,Phases,restraintDict,
            rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg,
            refPlotUpdate=refPlotUpdate)
        if IfOK:
            if len(covMatrix):      #empty for zero cycle refinement
                sigDict = dict(zip(varyList,sig))
                newCellDict = G2stMth.GetNewCellParms(parmDict,varyList)
                newAtomDict = G2stMth.ApplyXYZshifts(parmDict,varyList)
                covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                           'varyListStart':varyListStart,'Lastshft':Lastshft,
                           'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,
                           'newCellDict':newCellDict,'freshCOV':True}
                # add indirectly computed uncertainties into the esd dict
                sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList))
                G2stIO.PrintIndependentVars(parmDict,varyList,sigDict,pFile=printFile)
                G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
                G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
                G2stIO.SetPhaseData(parmDict,sigDict,Phases,rbIds,covData,restraintDict,printFile)
                G2stIO.SetISOmodes(parmDict,sigDict,Phases,printFile)
                G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,calcControls,
                    pFile=printFile,covMatrix=covMatrix,varyList=varyList)
                G2stIO.SetHistogramData(parmDict,sigDict,Histograms,calcControls,pFile=printFile)
                # check for variables outside their allowed range, reset and freeze them
                frozen = dropOOBvars(varyList,parmDict,sigDict,Controls,parmFrozenList)
                # covData['depSig'] = G2stIO.PhFrExtPOSig  # created in G2stIO.SetHistogramData, no longer used?
                covData['depSigDict'] = {i:(parmDict[i],sigDict[i]) for i in parmDict if i in sigDict}
                if len(frozen):
                    if 'msg' in Rvals:
                        Rvals['msg'] += '\n'
                    else:
                        Rvals['msg'] = ''
                    msg = ('Warning: {} variable(s) refined outside limits and were frozen ({} total frozen)'
                        .format(len(frozen),len(parmFrozenList))
                        )
                    G2fil.G2Print(msg)
                    Rvals['msg'] += msg
                elif len(parmFrozenList):
                    if 'msg' in Rvals:
                        Rvals['msg'] += '\n'
                    else:
                        Rvals['msg'] = ''
                    msg = ('Note: a total of {} variable(s) are frozen due to refining outside limits'
                        .format(len(parmFrozenList))
                        )
                    G2fil.G2Print('Note: ',msg)
                    Rvals['msg'] += msg
            G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,rigidbodyDict,covData,parmFrozenList,makeBack)
            printFile.close()
            G2fil.G2Print (' Refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst')
            G2fil.G2Print (' ***** Refinement successful *****')
        else:
            G2fil.G2Print ('****ERROR - Refinement failed',mode='error')
            if 'msg' in Rvals:
                G2fil.G2Print ('Note refinement problem:',mode='warn')
                G2fil.G2Print (Rvals['msg'],mode='warn')
            raise G2obj.G2Exception('**** ERROR: Refinement failed ****')
    except G2obj.G2RefineCancel as Msg:
        printFile.close()
        G2fil.G2Print (' ***** Refinement stopped *****')
        if not hasattr(Msg,'msg'): Msg.msg = str(Msg)
        if 'msg' in Rvals:
            Rvals['msg'] += '\n'
            Rvals['msg'] += Msg.msg
            if not dlg:
                G2fil.G2Print ('Note refinement problem:',mode='warn')
                G2fil.G2Print (Rvals['msg'],mode='warn')
        else:
            Rvals['msg'] = Msg.msg
        return False,Rvals
#    except G2obj.G2Exception as Msg:  # cell metric error, others?
    except Exception as Msg:  # cell metric error, others?
        if GSASIIpath.GetConfigValue('debug'):
            import traceback
            print(traceback.format_exc())        
        if not hasattr(Msg,'msg'): Msg.msg = str(Msg)
        printFile.close()
        G2fil.G2Print (' ***** Refinement error *****')
        if 'msg' in Rvals:
            Rvals['msg'] += '\n\n'
            Rvals['msg'] += Msg.msg
            if not dlg:
                G2fil.G2Print ('Note refinement problem:',mode='warn')
                G2fil.G2Print (Rvals['msg'],mode='warn')
        else:
            Rvals['msg'] = Msg.msg
        return False,Rvals

    # document the refinement further: RB, constraints, restraints, what's varied
    Rvals['varyList'] = 'Varied: ' + ', '.join(varyList)
    s = G2mv.VarRemapSumm()
    if s: Rvals['contrSumm'] = f'Constraints: {s}' 
    Rvals['restrSumm'] = G2stIO.SummRestraints(restraintDict)
    Rvals['RBsumm'] = ''
    for ph in Phases:
        s = ''
        for i in 'Vector','Residue':
            try: 
                l = len(Phases[ph]['RBModels'][i])
                if s: s += '; '
                s += f'{l} {i} bodies'
            except:
                pass
        if s:
            if not Rvals['RBsumm']: Rvals['RBsumm'] += 'Rigid Bodies: '
            Rvals['RBsumm'] += f'{ph}: {s}'
    
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
    elif 'msg' in Rvals:
        G2fil.G2Print ('Reported from refinement:',mode='warn')
        G2fil.G2Print (Rvals['msg'],mode='warn')

def CheckLeBail(Phases):
    '''Check if there is a LeBail extraction in any histogram

    :returns: True if there is at least one LeBail flag turned on, False otherwise
    '''
    for key in Phases:
        phase = Phases[key]
        for h in phase['Histograms']:
            #phase['Histograms'][h]
            if not phase['Histograms'][h]['Use']: continue
            try:
                if phase['Histograms'][h]['LeBail']:
                     return True
            except KeyError:    #HKLF & old gpx files
                pass
    return False

def DoNoFit(GPXfile,key):
    '''Compute the diffraction pattern with no refinement of parameters.

    TODO: At present, this will compute intensities all diffraction patterns
    in the project, but this likely can be made faster by dropping
    all the histograms except key from Histograms.

    :param str GPXfile: G2 .gpx file name
    :param str key: name of histogram to be computed
    :returns: the computed diffraction pattern for the selected histogram
    '''
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

    parmDict = {}
    Controls = G2stIO.GetControls(GPXfile)
    calcControls = {}
    calcControls.update(Controls)
    constrDict,fixedList = G2stIO.ReadConstraints(GPXfile)
    restraintDict = {}
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        G2fil.G2Print (' *** ERROR - you have no phases to refine! ***')
        return False,{'msg':'No phases'}
    if not Histograms:
        G2fil.G2Print (' *** ERROR - you have no data to refine with! ***')
        return False,{'msg':'No data'}
    if key not in Histograms:
        print(f"Error: no histogram by name {key}")
        return
    #TODO: Histograms = {key:Histograms[key]}
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,
         maxSSwave) = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,Print=False)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['EFtables'] = EFtables
    calcControls['ORBtables'] = ORBtables
    calcControls['BLtables'] = BLtables
    calcControls['MFtables'] = MFtables
    calcControls['maxSSwave'] = maxSSwave
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,Controls=calcControls,Print=False)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,Print=False)
    calcControls.update(controlDict)
    parmDict.update(rbDict)
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    G2stIO.GetFprime(calcControls,Histograms)
    
    M = G2stMth.errRefine([],[Histograms,Phases,restraintDict,rigidbodyDict],parmDict,[],calcControls,pawleyLookup,None)
    return Histograms[key]['Data'][3]

def DoLeBail(GPXfile,dlg=None,cycles=10,refPlotUpdate=None,seqList=None):
    '''Fit LeBail intensities without changes to any other refined parameters.
    This is a stripped-down version of :func:`Refine` that does not perform 
    any refinement cycles

    :param str GPXfile: G2 .gpx file name
    :param wx.ProgressDialog dlg: optional progress window to update. 
      Default is None, which means no calls are made to this. 
    :param int cycles: Number of LeBail cycles to perform
    :param function refPlotUpdate: Optional routine used to plot results. 
      Default is None, which means no calls are made to this. 
    :param list seqList: List of histograms to be processed. Default 
      is None which means that all used histograms in .gpx file are processed.
    '''
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

    #varyList = []
    parmDict = {}
    Controls = G2stIO.GetControls(GPXfile)
    calcControls = {}
    calcControls.update(Controls)
    constrDict,fixedList = G2stIO.ReadConstraints(GPXfile)
    restraintDict = {}
    Histograms_All,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if seqList:
        Histograms = {i:Histograms_All[i] for i in seqList}
    else:
        Histograms = Histograms_All
    if not Phases:
        G2fil.G2Print (' *** ERROR - you have no phases to refine! ***')
        return False,{'msg':'No phases'}
    if not Histograms:
        G2fil.G2Print (' *** ERROR - you have no data to refine with! ***')
        return False,{'msg':'No data'}
    if not CheckLeBail(Phases):
        msg = 'Warning: There are no histograms with LeBail extraction enabled'
        G2fil.G2Print ('*** '+msg+' ***')
        return False,{'msg': msg}
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,
         maxSSwave) = G2stIO.GetPhaseData(Phases,restraintDict,rbIds,Print=False)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['EFtables'] = EFtables
    calcControls['ORBtables'] = ORBtables
    calcControls['BLtables'] = BLtables
    calcControls['MFtables'] = MFtables
    calcControls['maxSSwave'] = maxSSwave
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,Controls=calcControls,Print=False)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,Print=False)
    calcControls.update(controlDict)
    parmDict.update(rbDict)
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    G2stIO.GetFprime(calcControls,Histograms)
    try:
        for i in range(cycles):
            M = G2stMth.errRefine([],[Histograms,Phases,restraintDict,rigidbodyDict],parmDict,[],calcControls,pawleyLookup,dlg)
            G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
            if refPlotUpdate is not None: refPlotUpdate(Histograms,i)
        Rvals = {}
        Rvals['chisq'] = np.sum(M**2)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs'])) # no variables

        covData = {'variables':0,'varyList':[],'sig':[],'Rvals':Rvals,'varyListStart':[],
            'covMatrix':None,'title':GPXfile,'freshCOV':True}   #np.zeros([0,0])?
          # ??  'newAtomDict':newAtomDict,'newCellDict':newCellDict,
        
        G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,rigidbodyDict,covData,[],True)
        G2fil.G2Print (' ***** LeBail fit completed *****')
        return True,Rvals
    except Exception as Msg:
        G2fil.G2Print (' ***** LeBail fit error *****')
        if not hasattr(Msg,'msg'): Msg.msg = str(Msg)
        if GSASIIpath.GetConfigValue('debug'):
            import traceback
            print(traceback.format_exc())        
        return False,{'msg':Msg.msg}

def phaseCheck(phaseVary,Phases,histogram):
    '''
    Removes unused parameters from phase varylist if phase not in histogram
    for seq refinement removes vars in "Fix FXU" and "FixedSeqVars" here
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
            removeVars = Phases[phase]['Histograms'][histogram].get('FixedSeqVars',[])
            newVary = [item for item in newVary if item not in removeVars]
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
        if h == 'FrozenList':
            continue
        preFrozenCount += len(Controls['parmFrozen'][h])    
    G2stIO.ShowControls(Controls,printFile,SeqRef=True,preFrozenCount=preFrozenCount)
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
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,BLtables,ORBtables,MFtables,maxSSwave) = \
        G2stIO.GetPhaseData(Phases,restraintDict,rbIds,Print=False,pFile=printFile,seqHistName='All')
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
    lasthist = ''
    for ihst,histogram in enumerate(histNames):
        if GSASIIpath.GetConfigValue('Show_timing'): t1 = time.time()
        G2fil.G2Print('\nRefining with '+str(histogram))
        G2mv.InitVars()
        (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave) = \
            G2stIO.GetPhaseData(Phases,restraintDict,rbIds,Print=False,pFile=printFile,seqHistName=histogram)
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['atomIndx'] = atomIndx
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['EFtables'] = EFtables
        calcControls['ORBtables'] = ORBtables
        calcControls['BLtables'] = BLtables
        calcControls['MFtables'] = MFtables
        calcControls['maxSSwave'] = maxSSwave
        if histogram not in Histograms:
            G2fil.G2Print("Error: not found!")
            raise G2obj.G2Exception("refining with invalid histogram {}".format(histogram))
        hId = Histograms[histogram]['hId']
        redphaseVary = phaseCheck(phaseVary,Phases,histogram)
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histo,Controls=calcControls,Print=False)
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
            for phase in Phases:
                if Phases[phase]['Histograms'][histogram].get('LeBail',False) and lasthist:
                    oldFsqs = Histograms[lasthist]['Reflection Lists'][phase]['RefList'].T[8:10]    #assume no superlattice!
                    newRefs = Histograms[histogram]['Reflection Lists'][phase]['RefList']
                    if len(newRefs) == len(oldFsqs.T):
                        newRefs.T[8:10] = copy.copy(oldFsqs)
                        # for i,ref in enumerate(newRefs):
                        #     ref[8:10] = oldFsqs.T[i]
                    else:
                        print('ERROR - mismatch in reflection list length bewteen %s and %s; no copy done'%(lasthist,histogram))
####TBD: if LeBail copy reflections here?
        elif histogram in SeqResult:  # update phase from last seq ref
            NewparmDict = SeqResult[histogram].get('parmDict',{})
            for parm in NewparmDict:
                if '::' in parm and parm in parmDict:
                    parmDict[parm] = NewparmDict[parm]
            
        G2stIO.GetFprime(calcControls,Histo)
        # do constraint processing (again, if called from GSASIIdataGUI.GSASII.OnSeqRefine)
        constrDict,fixedList = G2stIO.ReadConstraints(GPXfile,seqHist=hId)
        varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed

        msg = G2mv.EvaluateMultipliers(constrDict,phaseDict,hapDict,histDict)
        if msg:
            return False,'Unable to interpret multiplier(s): '+msg
      
        try:
            errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict,
                seqHistNum=hId,raiseException=True)
            constraintInfo = (groups,parmlist,constrDict,fixedList,ihst)
            G2mv.normParms(parmDict)
            G2mv.Map2Dict(parmDict,varyList)   # changes varyList
        except G2mv.ConstraintException:
            G2fil.G2Print (' *** ERROR - your constraints are internally inconsistent for histogram {}***'.format(hId))
            return False,' Constraint error'
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
        lasthist = histogram
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
        G2mv.Dict2Map(parmDict)  # impose constraints initially
        try:
            IfOK,Rvals,result,covMatrix,sig,Lastshft = RefineCore(Controls,Histo,Phases,restraintDict,
                rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg,
                refPlotUpdate=refPlotUpdate)
            try:
                shft = '%.4f'% Rvals['Max shft/sig']
            except:
                shft = '?'
            G2fil.G2Print ('  wR = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f, last delta chi = %.4f, last shft/sig = %s'%(
                Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']**2,Rvals['DelChi2'],shft))
            if Rvals.get('lamMax',0) >= 10.:
                msgs['steepestNum'] += 1
            if Rvals.get('Max shft/sig'):
                msgs['maxshift/sigma'].append(Rvals['Max shft/sig'])
            # add the uncertainties into the esd dictionary (sigDict)
            if not IfOK:
                G2fil.G2Print('***** Sequential refinement failed at histogram '+histogram,mode='warn')
                break
            sigDict = dict(zip(varyList,sig))
            # add indirectly computed uncertainties into the esd dict
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList))
        
            newCellDict = copy.deepcopy(G2stMth.GetNewCellParms(parmDict,varyList))
            newAtomDict = copy.deepcopy(G2stMth.ApplyXYZshifts(parmDict,varyList))
            SeqResult[histogram] = {
                'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                'varyListStart':varyListStart,'Lastshft':Lastshft,
                'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,
                'newCellDict':newCellDict,'depParmDict':{},
                'constraintInfo':constraintInfo,
                'parmDict':parmDict,
                }
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
            SeqResult[histogram]['RBsuDict'] = G2stMth.computeRBsu(parmDict,Phases,rigidbodyDict,
                            covMatrix,varyList,sig)
            G2stIO.SetISOmodes(parmDict,sigDict,Phases,None)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,None,ifPrint,
                                         pFile=printFile,covMatrix=covMatrix,varyList=varyList)
            G2stIO.SetHistogramData(parmDict,sigDict,Histo,None,ifPrint,printFile,seq=True)
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
            SeqResult[histogram]['depParmDict'] = {i:(parmDict[i],sigDict[i]) for i in sigDict if i not in varyList}
            
            
            G2stIO.SaveUpdatedHistogramsAndPhases(GPXfile,Histo,Phases,
                rigidbodyDict,SeqResult[histogram],Controls['parmFrozen'])
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
            if not hasattr(Msg,'msg'): Msg.msg = str(Msg)
            printFile.close()
            G2fil.G2Print (' ***** Refinement stopped *****')
            return False,Msg.msg
        except (G2obj.G2Exception,Exception) as Msg:  # cell metric error, others?
            if not hasattr(Msg,'msg'): Msg.msg = str(Msg)
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
        print('Error reading Sequential results\n',str(msg))
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
            BsumR = (Radii[Oatom[2]][0]+Radii[Tatom[2]][0])*Factor[0]
            AsumR = (Radii[Oatom[2]][1]+Radii[Tatom[2]][1])*Factor[1]
            for [Txyz,Top,Tunit,Spn] in G2spc.GenAtom(Tatom[3:6],SGData,False,Move=False):
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
        BVS = {}
        BVdat = {}
        Otyp = G2elem.FixValence(Oatom[2]).split('+')[0].split('-')[0]
        BVox = [BV for BV in atmdata.BVSoxid[Otyp] if '+' in BV]
        if len(BVox):
            BVS = {BV:0.0 for BV in BVox}
            BVdat = {BV:dict(zip(['O','F','Cl'],atmdata.BVScoeff[BV])) for BV in BVox}
            pvline = 'Bond Valence sums for: '
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
            if len(BVox):
                Tatm = G2elem.FixValence(DisAglData['TargAtoms'][dist[0]][2]).split('-')[0]
                if Tatm in ['O','F','Cl']:
                    for BV in BVox:
                        BVS[BV] += np.exp((BVdat[BV][Tatm]-dist[3])/0.37)                
            tunit = '[%2d%2d%2d]'% dist[1]
            MyPrint(('  %8s%10s+(%4d) %12s'%(AtomLabels[dist[0]].ljust(8),tunit.ljust(10),dist[2],val.center(12)))+line.rstrip())
        if len(BVox):
            MyPrint(80*'*')
            for BV in BVox:
                pvline += ' %s: %.2f  '%(BV,BVS[BV])
            MyPrint(pvline)

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

def do_refine(*args):
    'Called to run a refinement when this module is executed '
    starttime = time.time()
    #arg = sys.argv
    if len(args) >= 1:
        files = args
    elif len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        G2fil.G2Print ('ERROR GSASIIstrMain.do_refine error - missing filename')
        G2fil.G2Print ('Use "python GSASIIstrMain.py f1.gpx [f2.gpx f3.gpx...]" to run')
        G2fil.G2Print ('or call GSASIIstrMain.do_refine directly')
        sys.exit()
    for GPXfile in files:
        if not ospath.exists(GPXfile):
            G2fil.G2Print ('ERROR - '+GPXfile+" doesn't exist! Skipping.")
            continue
        # TODO: test below
        # figure out if this is a sequential refinement and call SeqRefine(GPXfile,None)
        #Controls = G2stIO.GetControls(GPXfile)
        #if Controls.get('Seq Data',[]): 
            Refine(GPXfile,None)
        #else:
        #    SeqRefine(GPXfile,None)
        G2fil.G2Print("Done with {}.\nExecution time {:.2f} sec.".format(GPXfile,time.time()-starttime))

if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    do_refine()
