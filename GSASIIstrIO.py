# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIstrIO: structure I/O routines*
-------------------------------------

'''
from __future__ import division, print_function
import platform
import os
import os.path as ospath
import time
import math
import copy
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import numpy as np
import numpy.ma as ma
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIobj as G2obj
import GSASIImapvars as G2mv
import GSASIImath as G2mth

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ateln2 = 8.0*math.log(2.0)

def cPickleLoad(fp):
    if '2' in platform.python_version_tuple()[0]:
        return cPickle.load(fp)
    else:
       return cPickle.load(fp,encoding='latin-1')

def GetControls(GPXfile):
    ''' Returns dictionary of control items found in GSASII gpx file

    :param str GPXfile: full .gpx file name
    :return: dictionary of control items
    '''
    Controls = copy.copy(G2obj.DefaultControls)
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Controls':
            Controls.update(datum[1])
    fl.close()
    return Controls
    
def GetConstraints(GPXfile):
    '''Read the constraints from the GPX file and interpret them

    called in :func:`ReadCheckConstraints`, :func:`GSASIIstrMain.Refine`
    and :func:`GSASIIstrMain.SeqRefine`. 
    '''
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Constraints':
            constList = []
            for item in datum[1]:
                if item.startswith('_'): continue
                constList += datum[1][item]
            fl.close()
            constDict,fixedList,ignored = ProcessConstraints(constList)
            if ignored:
                print ('%d Constraints were rejected. Was a constrained phase, histogram or atom deleted?'%ignored)
            return constDict,fixedList
    fl.close()
    raise Exception("No constraints in GPX file")
    
def ProcessConstraints(constList):
    """Interpret the constraints in the constList input into a dictionary, etc.
    All :class:`GSASIIobj.G2VarObj` objects are mapped to the appropriate
    phase/hist/atoms based on the object internals (random Ids). If this can't be
    done (if a phase has been deleted, etc.), the variable is ignored.
    If the constraint cannot be used due to too many dropped variables,
    it is counted as ignored.
    
    :param list constList: a list of lists where each item in the outer list
      specifies a constraint of some form, as described in the :mod:`GSASIIobj`
      :ref:`Constraint definition <Constraint_definitions_table>`.

    :returns:  a tuple of (constDict,fixedList,ignored) where:
      
      * constDict (list of dicts) contains the constraint relationships
      * fixedList (list) contains the fixed values for each type
        of constraint.
      * ignored (int) counts the number of invalid constraint items
        (should always be zero!)
    """
    constDict = []
    fixedList = []
    ignored = 0
    for item in constList:
        if item[-1] == 'h':
            # process a hold
            fixedList.append('0')
            var = str(item[0][1])
            if '?' not in var:
                constDict.append({var:0.0})
            else:
                ignored += 1
        elif item[-1] == 'f':
            # process a new variable
            fixedList.append(None)
            D = {}
            varyFlag = item[-2]
            varname = item[-3]
            for term in item[:-3]:
                var = str(term[1])
                if '?' not in var:
                    D[var] = term[0]
            if len(D) > 1:
                # add extra dict terms for input variable name and vary flag
                if varname is not None:                    
                    if varname.startswith('::'):
                        varname = varname[2:].replace(':',';')
                    else:
                        varname = varname.replace(':',';')
                    D['_name'] = '::' + varname
                D['_vary'] = varyFlag == True # force to bool
                constDict.append(D)
            else:
                ignored += 1
            #constFlag[-1] = ['Vary']
        elif item[-1] == 'c': 
            # process a contraint relationship
            D = {}
            for term in item[:-3]:
                var = str(term[1])
                if '?' not in var:
                    D[var] = term[0]
            if len(D) >= 1:
                fixedList.append(str(item[-3]))
                constDict.append(D)
            else:
                ignored += 1
        elif item[-1] == 'e':
            # process an equivalence
            firstmult = None
            eqlist = []
            for term in item[:-3]:
                if term[0] == 0: term[0] = 1.0
                var = str(term[1])
                if '?' in var: continue
                if firstmult is None:
                    firstmult = term[0]
                    firstvar = var
                else:
                    eqlist.append([var,firstmult/term[0]])
            if len(eqlist) > 0:
                G2mv.StoreEquivalence(firstvar,eqlist,False)
            else:
                ignored += 1
        else:
            ignored += 1
    return constDict,fixedList,ignored

def ReadCheckConstraints(GPXfile):
    '''Load constraints and related info and return any error or warning messages'''
    # init constraints
    G2mv.InitVars()    
    # get variables
    Histograms,Phases = GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        return 'Error: No phases or no histograms in phases!',''
    if not Histograms:
        return 'Error: no diffraction data',''
    constrDict,fixedList = GetConstraints(GPXfile) # load user constraints before internally generated ones
    rigidbodyDict = GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = GetRigidBodyModels(rigidbodyDict,Print=False)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,maxSSwave = \
        GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False) # generates atom symmetry constraints
    hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histograms,Print=False,resetRefList=True)
    histVary,histDict,controlDict = GetHistogramData(Histograms,Print=False)
    varyList = rbVary+phaseVary+hapVary+histVary
    errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
    if errmsg:
        # print some diagnostic info on the constraints
        print('Error in constraints:\n'+errmsg+
              '\nRefinement not possible due to conflict in constraints, see below:\n')
        print (G2mv.VarRemapShow(varyList,True))
    return errmsg, warnmsg
    
def makeTwinFrConstr(Phases,Histograms,hapVary):
    TwConstr = []
    TwFixed = []
    for Phase in Phases:
        pId = Phases[Phase]['pId']
        for Histogram in Phases[Phase]['Histograms']:
            try:
                hId = Histograms[Histogram]['hId']
                phfx = '%d:%d:'%(pId,hId)
                if phfx+'TwinFr:0' in hapVary:
                    TwFixed.append('1.0')     #constraint value
                    nTwin = len(Phases[Phase]['Histograms'][Histogram]['Twins'])
                    TwConstr.append({phfx+'TwinFr:'+str(i):'1.0' for i in range(nTwin)})
            except KeyError:    #unused histograms?
                pass
    return TwConstr,TwFixed   
    
def GetRestraints(GPXfile):
    '''Read the restraints from the GPX file.
    Throws an exception if not found in the .GPX file
    '''
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Restraints':
            restraintDict = datum[1]
    fl.close()
    return restraintDict
    
def GetRigidBodies(GPXfile):
    '''Read the rigid body models from the GPX file
    '''
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Rigid bodies':
            rigidbodyDict = datum[1]
    fl.close()
    return rigidbodyDict
        
def GetFprime(controlDict,Histograms):
    'Needs a doc string'
    FFtables = controlDict['FFtables']
    if not FFtables:
        return
    histoList = list(Histograms.keys())
    histoList.sort()
    for histogram in histoList:
        if histogram[:4] in ['PWDR','HKLF']:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            if 'X' in controlDict[hfx+'histType']:
                keV = controlDict[hfx+'keV']
                for El in FFtables:
                    Orbs = G2el.GetXsectionCoeff(El.split('+')[0].split('-')[0])
                    FP,FPP,Mu = G2el.FPcalc(Orbs, keV)
                    FFtables[El][hfx+'FP'] = FP
                    FFtables[El][hfx+'FPP'] = FPP
                    
def PrintFprime(FFtables,pfx,pFile):
    pFile.write('\n Resonant form factors:\n')
    Elstr = ' Element:'
    FPstr = " f'     :"
    FPPstr = ' f"     :'
    for El in FFtables:
        Elstr += ' %8s'%(El)
        FPstr += ' %8.3f'%(FFtables[El][pfx+'FP'])
        FPPstr += ' %8.3f'%(FFtables[El][pfx+'FPP'])
    pFile.write(Elstr+'\n')
    pFile.write(FPstr+'\n')
    pFile.write(FPPstr+'\n')
            
def GetPhaseNames(GPXfile):
    ''' Returns a list of phase names found under 'Phases' in GSASII gpx file

    :param str GPXfile: full .gpx file name
    :return: list of phase names
    '''
    fl = open(GPXfile,'rb')
    PhaseNames = []
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if 'Phases' == datum[0]:
            for datus in data[1:]:
                PhaseNames.append(datus[0])
    fl.close()
    return PhaseNames

def GetAllPhaseData(GPXfile,PhaseName):
    ''' Returns the entire dictionary for PhaseName from GSASII gpx file

    :param str GPXfile: full .gpx file name
    :param str PhaseName: phase name
    :return: phase dictionary
    '''        
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if 'Phases' == datum[0]:
            for datus in data[1:]:
                if datus[0] == PhaseName:
                    break
    fl.close()
    return datus[1]
    
def GetHistograms(GPXfile,hNames):
    """ Returns a dictionary of histograms found in GSASII gpx file

    :param str GPXfile: full .gpx file name
    :param str hNames: list of histogram names
    :return: dictionary of histograms (types = PWDR & HKLF)

    """
    fl = open(GPXfile,'rb')
    Histograms = {}
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        hist = datum[0]
        if hist in hNames:
            if 'PWDR' in hist[:4]:
                PWDRdata = {}
                PWDRdata.update(datum[1][0])        #weight factor
                PWDRdata['Data'] = ma.array(ma.getdata(datum[1][1]))          #masked powder data arrays/clear previous masks
                PWDRdata[data[2][0]] = data[2][1]       #Limits & excluded regions (if any)
                PWDRdata[data[3][0]] = data[3][1]       #Background
                PWDRdata[data[4][0]] = data[4][1]       #Instrument parameters
                PWDRdata[data[5][0]] = data[5][1]       #Sample parameters
                try:
                    PWDRdata[data[9][0]] = data[9][1]       #Reflection lists might be missing
                except IndexError:
                    PWDRdata['Reflection Lists'] = {}
                PWDRdata['Residuals'] = {}
    
                Histograms[hist] = PWDRdata
            elif 'HKLF' in hist[:4]:
                HKLFdata = {}
                HKLFdata.update(datum[1][0])        #weight factor
#patch
                if 'list' in str(type(datum[1][1])):
                #if isinstance(datum[1][1],list):
                    RefData = {'RefList':[],'FF':{}}
                    for ref in datum[1][1]:
                        RefData['RefList'].append(ref[:11]+[ref[13],])
                    RefData['RefList'] = np.array(RefData['RefList'])
                    datum[1][1] = RefData
#end patch
                datum[1][1]['FF'] = {}
                HKLFdata['Data'] = datum[1][1]
                HKLFdata[data[1][0]] = data[1][1]       #Instrument parameters
                HKLFdata['Reflection Lists'] = None
                HKLFdata['Residuals'] = {}
                Histograms[hist] = HKLFdata           
    fl.close()
    return Histograms
    
def GetHistogramNames(GPXfile,hType):
    """ Returns a list of histogram names found in GSASII gpx file

    :param str GPXfile: full .gpx file name
    :param str hType: list of histogram types
    :return: list of histogram names (types = PWDR & HKLF)

    """
    fl = open(GPXfile,'rb')
    HistogramNames = []
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0][:4] in hType:
            HistogramNames.append(datum[0])
    fl.close()
    return HistogramNames
    
def GetUsedHistogramsAndPhases(GPXfile):
    ''' Returns all histograms that are found in any phase
    and any phase that uses a histogram. This also
    assigns numbers to used phases and histograms by the
    order they appear in the file. 

    :param str GPXfile: full .gpx file name
    :returns: (Histograms,Phases)

     * Histograms = dictionary of histograms as {name:data,...}
     * Phases = dictionary of phases that use histograms

    '''
    phaseNames = GetPhaseNames(GPXfile)
    histoList = GetHistogramNames(GPXfile,['PWDR','HKLF'])
    allHistograms = GetHistograms(GPXfile,histoList)
    phaseData = {}
    for name in phaseNames: 
        phaseData[name] =  GetAllPhaseData(GPXfile,name)
    Histograms = {}
    Phases = {}
    for phase in phaseData:
        Phase = phaseData[phase]
        if Phase['General']['Type'] == 'faulted': continue      #don't use faulted phases!
        if Phase['Histograms']:
            if phase not in Phases:
                pId = phaseNames.index(phase)
                Phase['pId'] = pId
                Phases[phase] = Phase
            for hist in Phase['Histograms']:
                if 'Use' not in Phase['Histograms'][hist]:      #patch
                    Phase['Histograms'][hist]['Use'] = True         
                if hist not in Histograms and Phase['Histograms'][hist]['Use']:
                    try:
                        Histograms[hist] = allHistograms[hist]
                        hId = histoList.index(hist)
                        Histograms[hist]['hId'] = hId
                    except KeyError: # would happen if a referenced histogram were
                        # renamed or deleted
                        print('For phase "'+phase+
                              '" unresolved reference to histogram "'+hist+'"')
    G2obj.IndexAllIds(Histograms=Histograms,Phases=Phases)
    return Histograms,Phases
    
def getBackupName(GPXfile,makeBack):
    '''
    Get the name for the backup .gpx file name
    
    :param str GPXfile: full .gpx file name
    :param bool makeBack: if True the name of a new file is returned, if
      False the name of the last file that exists is returned
    :returns: the name of a backup file
    
    '''
    GPXpath,GPXname = ospath.split(GPXfile)
    if GPXpath == '': GPXpath = '.'
    Name = ospath.splitext(GPXname)[0]
    files = os.listdir(GPXpath)
    last = 0
    for name in files:
        name = name.split('.')
        if len(name) == 3 and name[0] == Name and 'bak' in name[1]:
            if makeBack:
                last = max(last,int(name[1].strip('bak'))+1)
            else:
                last = max(last,int(name[1].strip('bak')))
    GPXback = ospath.join(GPXpath,ospath.splitext(GPXname)[0]+'.bak'+str(last)+'.gpx')
    return GPXback    
        
def GPXBackup(GPXfile,makeBack=True):
    '''
    makes a backup of the current .gpx file (?)
    
    :param str GPXfile: full .gpx file name
    :param bool makeBack: if True (default), the backup is written to
      a new file; if False, the last backup is overwritten
    :returns: the name of the backup file that was written
    '''
    import distutils.file_util as dfu
    GPXback = getBackupName(GPXfile,makeBack)
    tries = 0
    while True:
        try:
            dfu.copy_file(GPXfile,GPXback)
            break
        except:
            tries += 1
            if tries > 10:
                return GPXfile  #failed!
            time.sleep(1)           #just wait a second!         
    return GPXback

def SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,RigidBodies,CovData,makeBack=True):
    ''' Updates gpxfile from all histograms that are found in any phase
    and any phase that used a histogram. Also updates rigid body definitions.


    :param str GPXfile: full .gpx file name
    :param dict Histograms: dictionary of histograms as {name:data,...}
    :param dict Phases: dictionary of phases that use histograms
    :param dict RigidBodies: dictionary of rigid bodies
    :param dict CovData: dictionary of refined variables, varyList, & covariance matrix
    :param bool makeBack: True if new backup of .gpx file is to be made; else use the last one made

    '''
                        
    import distutils.file_util as dfu
    GPXback = GPXBackup(GPXfile,makeBack)
    print ('Read from file:'+GPXback)
    print ('Save to file  :'+GPXfile)
    infile = open(GPXback,'rb')
    outfile = open(GPXfile,'wb')
    while True:
        try:
            data = cPickleLoad(infile)
        except EOFError:
            break
        datum = data[0]
#        print 'read: ',datum[0]
        if datum[0] == 'Phases':
            for iphase in range(len(data)):
                if data[iphase][0] in Phases:
                    phaseName = data[iphase][0]
                    data[iphase][1].update(Phases[phaseName])
        elif datum[0] == 'Covariance':
            data[0][1] = CovData
        elif datum[0] == 'Rigid bodies':
            data[0][1] = RigidBodies
        try:
            histogram = Histograms[datum[0]]
#            print 'found ',datum[0]
            data[0][1][1] = histogram['Data']
            data[0][1][0].update(histogram['Residuals'])
            for datus in data[1:]:
#                print '    read: ',datus[0]
                if datus[0] in ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']:
                    datus[1] = histogram[datus[0]]
        except KeyError:
            pass
        try:                        
            cPickle.dump(data,outfile,1)
        except AttributeError:
            print ('ERROR - bad data in least squares result')
            infile.close()
            outfile.close()
            dfu.copy_file(GPXback,GPXfile)
            print ('GPX file save failed - old version retained')
            return
            
    print ('GPX file save successful')
    
def GetSeqResult(GPXfile):
    '''
    Needs doc string
    
    :param str GPXfile: full .gpx file name
    '''
    fl = open(GPXfile,'rb')
    SeqResult = {}
    while True:
        try:
            data = cPickleLoad(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Sequential results':
            SeqResult = datum[1]
    fl.close()
    return SeqResult
    
def SetSeqResult(GPXfile,Histograms,SeqResult):
    '''
    Needs doc string
    
    :param str GPXfile: full .gpx file name
    '''
    GPXback = GPXBackup(GPXfile)
    print ('Read from file:'+GPXback)
    print ('Save to file  :'+GPXfile)
    infile = open(GPXback,'rb')
    outfile = open(GPXfile,'wb')
    while True:
        try:
            data = cPickleLoad(infile)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Sequential results':
            data[0][1] = SeqResult
        # reset the Copy Next flag, since it should not be needed twice in a row
        if datum[0] == 'Controls':
            data[0][1]['Copy2Next'] = False
        try:
            histogram = Histograms[datum[0]]
            data[0][1][1] = list(histogram['Data'])
            for datus in data[1:]:
                if datus[0] in ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']:
                    datus[1] = histogram[datus[0]]
        except KeyError:
            pass
                                
        cPickle.dump(data,outfile,1)
    infile.close()
    outfile.close()
    print ('GPX file save successful')
                        
def ShowBanner(pFile=None):
    'Print authorship, copyright and citation notice'
    pFile.write(80*'*'+'\n')
    pFile.write('   General Structure Analysis System-II Crystal Structure Refinement\n')
    pFile.write('              by Robert B. Von Dreele & Brian H. Toby\n')
    pFile.write('                Argonne National Laboratory(C), 2010\n')
    pFile.write(' This product includes software developed by the UChicago Argonne, LLC,\n')
    pFile.write('            as Operator of Argonne National Laboratory.\n')
    pFile.write('                          Please cite:\n')
    pFile.write('   B.H. Toby & R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013)\n')

    pFile.write(80*'*'+'\n')

def ShowControls(Controls,pFile=None,SeqRef=False):
    'Print controls information'
    pFile.write(' Least squares controls:\n')
    pFile.write(' Refinement type: %s\n'%Controls['deriv type'])
    if 'Hessian' in Controls['deriv type']:
        pFile.write(' Maximum number of cycles: %d\n'%Controls['max cyc'])
    else:
        pFile.write(' Minimum delta-M/M for convergence: %.2g\n'%(Controls['min dM/M']))
    pFile.write(' Regularize hydrogens (if any): %s\n'%Controls.get('HatomFix',False))
    pFile.write(' Initial shift factor: %.3f\n'%(Controls['shift factor']))
    if SeqRef:
        pFile.write(' Sequential refinement controls:\n')
        pFile.write(' Copy of histogram results to next: %s\n'%(Controls['Copy2Next']))
        pFile.write(' Process histograms in reverse order: %s\n'%(Controls['Reverse Seq']))
    
def GetPawleyConstr(SGLaue,PawleyRef,im,pawleyVary):
    'needs a doc string'
#    if SGLaue in ['-1','2/m','mmm']:
#        return                      #no Pawley symmetry required constraints
    eqvDict = {}
    for i,varyI in enumerate(pawleyVary):
        eqvDict[varyI] = []
        refI = int(varyI.split(':')[-1])
        ih,ik,il = PawleyRef[refI][:3]
        dspI = PawleyRef[refI][4+im]
        for varyJ in pawleyVary[i+1:]:
            refJ = int(varyJ.split(':')[-1])
            jh,jk,jl = PawleyRef[refJ][:3]
            dspJ = PawleyRef[refJ][4+im]
            if SGLaue in ['4/m','4/mmm']:
                isum = ih**2+ik**2
                jsum = jh**2+jk**2
                if abs(il) == abs(jl) and isum == jsum:
                    eqvDict[varyI].append(varyJ) 
            elif SGLaue in ['3R','3mR']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2+jl**2
                isum2 = ih*ik+ih*il+ik*il
                jsum2 = jh*jk+jh*jl+jk*jl
                if isum == jsum and isum2 == jsum2:
                    eqvDict[varyI].append(varyJ) 
            elif SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                isum = ih**2+ik**2+ih*ik
                jsum = jh**2+jk**2+jh*jk
                if abs(il) == abs(jl) and isum == jsum:
                    eqvDict[varyI].append(varyJ) 
            elif SGLaue in ['m3','m3m']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2+jl**2
                if isum == jsum:
                    eqvDict[varyI].append(varyJ)
            elif abs(dspI-dspJ)/dspI < 1.e-4:
                eqvDict[varyI].append(varyJ) 
    for item in pawleyVary:
        if eqvDict[item]:
            for item2 in pawleyVary:
                if item2 in eqvDict[item]:
                    eqvDict[item2] = []
            G2mv.StoreEquivalence(item,eqvDict[item])
                    
def cellVary(pfx,SGData): 
    'needs a doc string'
    if SGData['SGLaue'] in ['-1',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A5']
        elif SGData['SGUniq'] == 'b':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A4']
        else:
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3']
    elif SGData['SGLaue'] in ['mmm',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2']
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        G2mv.StoreEquivalence(pfx+'A0',(pfx+'A1',))
        return [pfx+'A0',pfx+'A1',pfx+'A2']
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        G2mv.StoreEquivalence(pfx+'A0',(pfx+'A1',pfx+'A3',))
        return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3']
    elif SGData['SGLaue'] in ['3R', '3mR']:
        G2mv.StoreEquivalence(pfx+'A0',(pfx+'A1',pfx+'A2',))
        G2mv.StoreEquivalence(pfx+'A3',(pfx+'A4',pfx+'A5',))
        return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']                       
    elif SGData['SGLaue'] in ['m3m','m3']:
        G2mv.StoreEquivalence(pfx+'A0',(pfx+'A1',pfx+'A2',))
        return [pfx+'A0',pfx+'A1',pfx+'A2']
        
def modVary(pfx,SSGData):
    vary = []
    for i,item in enumerate(SSGData['modSymb']):
        if item in ['a','b','g']:
            vary.append(pfx+'mV%d'%(i))
    return vary
        
################################################################################
##### Rigid Body Models and not General.get('doPawley')
################################################################################
        
def GetRigidBodyModels(rigidbodyDict,Print=True,pFile=None):
    'needs a doc string'
    
    def PrintResRBModel(RBModel):
        atNames = RBModel['atNames']
        rbRef = RBModel['rbRef']
        rbSeq = RBModel['rbSeq']
        pFile.write('Residue RB name: %s no.atoms: %d, No. times used: %d\n'%
            (RBModel['RBname'],len(RBModel['rbTypes']),RBModel['useCount']))
        pFile.write('    At name       x          y          z\n')
        for name,xyz in zip(atNames,RBModel['rbXYZ']):
            pFile.write('  %8s %10.4f %10.4f %10.4f\n'%(name,xyz[0],xyz[1],xyz[2]))
        pFile.write('Orientation defined by: %s -> %s & %s -> %s\n'%
            (atNames[rbRef[0]],atNames[rbRef[1]],atNames[rbRef[0]],atNames[rbRef[2]]))
        if rbSeq:
            for i,rbseq in enumerate(rbSeq):
#                nameLst = [atNames[i] for i in rbseq[3]]
                pFile.write('Torsion sequence %d Bond: %s - %s riding: %s\n'%
                    (i,atNames[rbseq[0]],atNames[rbseq[1]],str([atNames[i] for i in rbseq[3]])))
        
    def PrintVecRBModel(RBModel):
        rbRef = RBModel['rbRef']
        atTypes = RBModel['rbTypes']
        pFile.write('Vector RB name: %s no.atoms: No. times used: %d\n'%
            (RBModel['RBname'],len(RBModel['rbTypes']),RBModel['useCount']))
        for i in range(len(RBModel['VectMag'])):
            pFile.write('Vector no.: %d Magnitude: %8.4f Refine? %s\n'%(i,RBModel['VectMag'][i],RBModel['VectRef'][i]))
            pFile.write('  No. Type     vx         vy         vz\n')
            for j,[name,xyz] in enumerate(zip(atTypes,RBModel['rbVect'][i])):
                pFile.write('  %d   %2s %10.4f %10.4f %10.4f\n'%(j,name,xyz[0],xyz[1],xyz[2]))
        pFile.write('  No. Type      x          y          z\n')
        for i,[name,xyz] in enumerate(zip(atTypes,RBModel['rbXYZ'])):
            pFile.write('  %d   %2s %10.4f %10.4f %10.4f\n'%(i,name,xyz[0],xyz[1],xyz[2]))
        pFile.write('Orientation defined by: atom %s -> atom %s & atom %s -> atom %s\n'%
            (rbRef[0],rbRef[1],rbRef[0],rbRef[2]))
    rbVary = []
    rbDict = {}
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    if len(rbIds['Vector']):
        for irb,item in enumerate(rbIds['Vector']):
            if rigidbodyDict['Vector'][item]['useCount']:
                RBmags = rigidbodyDict['Vector'][item]['VectMag']
                RBrefs = rigidbodyDict['Vector'][item]['VectRef']
                for i,[mag,ref] in enumerate(zip(RBmags,RBrefs)):
                    pid = '::RBV;'+str(i)+':'+str(irb)
                    rbDict[pid] = mag
                    if ref:
                        rbVary.append(pid)
                if Print:
                    pFile.write('\nVector rigid body model:\n')
                    PrintVecRBModel(rigidbodyDict['Vector'][item])
    if len(rbIds['Residue']):
        for item in rbIds['Residue']:
            if rigidbodyDict['Residue'][item]['useCount']:
                if Print:
                    pFile.write('\nResidue rigid body model:\n')
                    PrintResRBModel(rigidbodyDict['Residue'][item])
    return rbVary,rbDict
    
def SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,pFile=None):
    'needs a doc string'
    
    def PrintRBVectandSig(VectRB,VectSig):
        pFile.write('\n Rigid body vector magnitudes for %s:\n'%VectRB['RBname'])
        namstr = '  names :'
        valstr = '  values:'
        sigstr = '  esds  :'
        for i,[val,sig] in enumerate(zip(VectRB['VectMag'],VectSig)):
            namstr += '%12s'%('Vect '+str(i))
            valstr += '%12.4f'%(val)
            if sig:
                sigstr += '%12.4f'%(sig)
            else:
                sigstr += 12*' '
        pFile.write(namstr+'\n')
        pFile.write(valstr+'\n')
        pFile.write(sigstr+'\n')        
        
    RBIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})  #these are lists of rbIds
    if not RBIds['Vector']:
        return
    for irb,item in enumerate(RBIds['Vector']):
        if rigidbodyDict['Vector'][item]['useCount']:
            VectSig = []
            RBmags = rigidbodyDict['Vector'][item]['VectMag']
            for i,mag in enumerate(RBmags):
                name = '::RBV;'+str(i)+':'+str(irb)
                if name in sigDict:
                    VectSig.append(sigDict[name])
            PrintRBVectandSig(rigidbodyDict['Vector'][item],VectSig)    
        
################################################################################
##### Phase data
################################################################################        
                    
def GetPhaseData(PhaseData,RestraintDict={},rbIds={},Print=True,pFile=None,seqRef=False):
    'needs a doc string'
            
    def PrintFFtable(FFtable):
        pFile.write('\n X-ray scattering factors:\n')
        pFile.write('   Symbol     fa                                      fb                                      fc\n')
        pFile.write(99*'-'+'\n')
        for Ename in FFtable:
            ffdata = FFtable[Ename]
            fa = ffdata['fa']
            fb = ffdata['fb']
            pFile.write(' %8s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n'%
                (Ename.ljust(8),fa[0],fa[1],fa[2],fa[3],fb[0],fb[1],fb[2],fb[3],ffdata['fc']))
                
    def PrintMFtable(MFtable):
        pFile.write('\n <j0> Magnetic scattering factors:\n')
        pFile.write('   Symbol     mfa                                    mfb                                     mfc\n')
        pFile.write(99*'-'+'\n')
        for Ename in MFtable:
            mfdata = MFtable[Ename]
            fa = mfdata['mfa']
            fb = mfdata['mfb']
            pFile.write(' %8s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n'%
                (Ename.ljust(8),fa[0],fa[1],fa[2],fa[3],fb[0],fb[1],fb[2],fb[3],mfdata['mfc']))
        pFile.write('\n <j2> Magnetic scattering factors:\n')
        pFile.write('   Symbol     nfa                                    nfb                                     nfc\n')
        pFile.write(99*'-'+'\n')
        for Ename in MFtable:
            mfdata = MFtable[Ename]
            fa = mfdata['nfa']
            fb = mfdata['nfb']
            pFile.write(' %8s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n'%
                (Ename.ljust(8),fa[0],fa[1],fa[2],fa[3],fb[0],fb[1],fb[2],fb[3],mfdata['nfc']))
                
    def PrintBLtable(BLtable):
        pFile.write('\n Neutron scattering factors:\n')
        pFile.write('   Symbol   isotope       mass       b       resonant terms\n')
        pFile.write(99*'-'+'\n')
        for Ename in BLtable:
            bldata = BLtable[Ename]
            isotope = bldata[0]
            mass = bldata[1]['Mass']
            if 'BW-LS' in bldata[1]:
                bres = bldata[1]['BW-LS']
                blen = 0
            else:
                blen = bldata[1]['SL'][0]
                bres = []
            line = ' %8s%11s %10.3f %8.3f'%(Ename.ljust(8),isotope.center(11),mass,blen)
            for item in bres:
                line += '%10.5g'%(item)
            pFile.write(line+'\n')
            
    def PrintRBObjects(resRBData,vecRBData):
        
        def PrintRBThermals():
            tlstr = ['11','22','33','12','13','23']
            sstr = ['12','13','21','23','31','32','AA','BB']
            TLS = RB['ThermalMotion'][1]
            TLSvar = RB['ThermalMotion'][2]
            if 'T' in RB['ThermalMotion'][0]:
                pFile.write('TLS data\n')
                text = ''
                for i in range(6):
                    text += 'T'+tlstr[i]+' %8.4f %s '%(TLS[i],str(TLSvar[i])[0])
                pFile.write(text+'\n')
                if 'L' in RB['ThermalMotion'][0]: 
                    text = ''
                    for i in range(6,12):
                        text += 'L'+tlstr[i-6]+' %8.2f %s '%(TLS[i],str(TLSvar[i])[0])
                    pFile.write(text+'\n')
                if 'S' in RB['ThermalMotion'][0]:
                    text = ''
                    for i in range(12,20):
                        text += 'S'+sstr[i-12]+' %8.3f %s '%(TLS[i],str(TLSvar[i])[0])
                    pFile.write(text+'\n')
            if 'U' in RB['ThermalMotion'][0]:
                pFile.write('Uiso data\n')
                text = 'Uiso'+' %10.3f %s'%(TLS[0],str(TLSvar[0])[0])           
                pFile.write(text+'\n')
            
        if len(resRBData):
            for RB in resRBData:
                Oxyz = RB['Orig'][0]
                Qrijk = RB['Orient'][0]
                Angle = 2.0*acosd(Qrijk[0])
                pFile.write('\nRBObject %s at %10.4f %10.4f %10.4f Refine? %s\n'%
                    (RB['RBname'],Oxyz[0],Oxyz[1],Oxyz[2],RB['Orig'][1]))
                pFile.write('Orientation angle,vector: %10.3f %10.4f %10.4f %10.4f Refine? %s\n'%
                    (Angle,Qrijk[1],Qrijk[2],Qrijk[3],RB['Orient'][1]))
                Torsions = RB['Torsions']
                if len(Torsions):
                    text = 'Torsions: '
                    for torsion in Torsions:
                        text += '%10.4f Refine? %s'%(torsion[0],torsion[1])
                    pFile.write(text+'\n')
                PrintRBThermals()
        if len(vecRBData):
            for RB in vecRBData:
                Oxyz = RB['Orig'][0]
                Qrijk = RB['Orient'][0]
                Angle = 2.0*acosd(Qrijk[0])
                pFile.write('\nRBObject %s',' at %10.4f %10.4f %10.4f Refine? %s\n'%
                    (RB['RBname'],Oxyz[0],Oxyz[1],Oxyz[2],RB['Orig'][1]))
                pFile.write('Orientation angle,vector: %10.3f %10.4f %10.4f %10.4f Refine? %s\n'%
                    (Angle,Qrijk[1],Qrijk[2],Qrijk[3],RB['Orient'][1]))
                PrintRBThermals()
                
    def PrintAtoms(General,Atoms):
        cx,ct,cs,cia = General['AtomPtrs']
        pFile.write('\n Atoms:\n')
        line = '   name    type  refine?   x         y         z    '+ \
            '  frac site sym  mult I/A   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'macromolecular':
            line = ' res no residue chain'+line
        pFile.write(line+'\n')
        if General['Type'] in ['nuclear','magnetic','faulted',]:
            pFile.write(135*'-'+'\n')
            for i,at in enumerate(Atoms):
                line = '%7s'%(at[ct-1])+'%7s'%(at[ct])+'%7s'%(at[ct+1])+'%10.5f'%(at[cx])+'%10.5f'%(at[cx+1])+ \
                    '%10.5f'%(at[cx+2])+'%8.3f'%(at[cx+3])+'%7s'%(at[cs])+'%5d'%(at[cs+1])+'%5s'%(at[cia])
                if at[cia] == 'I':
                    line += '%8.5f'%(at[cia+1])+48*' '
                else:
                    line += 8*' '
                    for j in range(6):
                        line += '%8.5f'%(at[cia+2+j])
                pFile.write(line+'\n')
        elif General['Type'] == 'macromolecular':
            pFile.write(135*'-'+'\n')            
            for i,at in enumerate(Atoms):
                line = '%7s'%(at[0])+'%7s'%(at[1])+'%7s'%(at[2])+'%7s'%(at[ct-1])+'%7s'%(at[ct])+'%7s'%(at[ct+1])+'%10.5f'%(at[cx])+'%10.5f'%(at[cx+1])+ \
                    '%10.5f'%(at[cx+2])+'%8.3f'%(at[cx+3])+'%7s'%(at[cs])+'%5d'%(at[cs+1])+'%5s'%(at[cia])
                if at[cia] == 'I':
                    line += '%8.4f'%(at[cia+1])+48*' '
                else:
                    line += 8*' '
                    for j in range(6):
                        line += '%8.4f'%(at[cia+2+j])
                pFile.write(line+'\n')
                
    def PrintMoments(General,Atoms):
        cx,ct,cs,cia = General['AtomPtrs']
        cmx = cx+4
        AtInfo = dict(zip(General['AtomTypes'],General['Lande g']))
        pFile.write('\n Magnetic moments:\n')
        line = '   name    type  refine?  Mx        My        Mz    '
        pFile.write(line+'\n')
        pFile.write(135*'-'+'\n')
        for i,at in enumerate(Atoms):
            if AtInfo[at[ct]]:
                line = '%7s'%(at[ct-1])+'%7s'%(at[ct])+'%7s'%(at[ct+1])+'%10.5f'%(at[cmx])+'%10.5f'%(at[cmx+1])+ \
                    '%10.5f'%(at[cmx+2])
                pFile.write(line+'\n')
        
                
    def PrintWaves(General,Atoms):
        cx,ct,cs,cia = General['AtomPtrs']
        pFile.write('\n Modulation waves\n')
        names = {'Sfrac':['Fsin','Fcos','Fzero','Fwid'],'Spos':['Xsin','Ysin','Zsin','Xcos','Ycos','Zcos','Tmin','Tmax','Xmax','Ymax','Zmax'],
            'Sadp':['U11sin','U22sin','U33sin','U12sin','U13sin','U23sin','U11cos','U22cos',
            'U33cos','U12cos','U13cos','U23cos'],'Smag':['MXsin','MYsin','MZsin','MXcos','MYcos','MZcos']}
        pFile.write(135*'-'+'\n')
        for i,at in enumerate(Atoms):
            AtomSS = at[-1]['SS1']
            for Stype in ['Sfrac','Spos','Sadp','Smag']:
                Waves = AtomSS[Stype]
                if len(Waves):
                    pFile.write(' atom: %s, site sym: %s, type: %s wave type: %s:\n'%
                        (at[ct-1],at[cs],Stype,AtomSS['waveType']))
                for iw,wave in enumerate(Waves):                    
                    line = ''
                    if AtomSS['waveType'] in ['Block','ZigZag'] and Stype == 'Spos' and not iw:
                        for item in names[Stype][6:]:
                            line += '%8s '%(item)                        
                    else:
                        if Stype == 'Spos':
                            for item in names[Stype][:6]:
                                line += '%8s '%(item)
                        else:
                            for item in names[Stype]:
                                line += '%8s '%(item)
                    pFile.write(line+'\n')
                    line = ''
                    for item in wave[0]:
                        line += '%8.4f '%(item)
                    line += ' Refine? '+str(wave[1])
                    pFile.write(line+'\n')
        
    def PrintTexture(textureData):
        topstr = '\n Spherical harmonics texture: Order:' + \
            str(textureData['Order'])
        if textureData['Order']:
            pFile.write('%s Refine? %s\n'%(topstr,textureData['SH Coeff'][0]))
        else:
            pFile.write(topstr+'\n')
            return
        names = ['omega','chi','phi']
        line = '\n'
        for name in names:
            line += ' SH '+name+':'+'%12.4f'%(textureData['Sample '+name][1])+' Refine? '+str(textureData['Sample '+name][0])
        pFile.write(line+'\n')
        pFile.write('\n Texture coefficients:\n')
        SHcoeff = textureData['SH Coeff'][1]
        SHkeys = list(SHcoeff.keys())
        nCoeff = len(SHcoeff)
        nBlock = nCoeff//10+1
        iBeg = 0
        iFin = min(iBeg+10,nCoeff)
        for block in range(nBlock):
            ptlbls = ' names :'
            ptstr =  ' values:'
            for item in SHkeys[iBeg:iFin]:
                ptlbls += '%12s'%(item)
                ptstr += '%12.4f'%(SHcoeff[item]) 
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            iBeg += 10
            iFin = min(iBeg+10,nCoeff)
        
    def MakeRBParms(rbKey,phaseVary,phaseDict):
        rbid = str(rbids.index(RB['RBId']))
        pfxRB = pfx+'RB'+rbKey+'P'
        pstr = ['x','y','z']
        ostr = ['a','i','j','k']
        for i in range(3):
            name = pfxRB+pstr[i]+':'+str(iRB)+':'+rbid
            phaseDict[name] = RB['Orig'][0][i]
            if RB['Orig'][1]:
                phaseVary += [name,]
        pfxRB = pfx+'RB'+rbKey+'O'
        for i in range(4):
            name = pfxRB+ostr[i]+':'+str(iRB)+':'+rbid
            phaseDict[name] = RB['Orient'][0][i]
            if RB['Orient'][1] == 'AV' and i:
                phaseVary += [name,]
            elif RB['Orient'][1] == 'A' and not i:
                phaseVary += [name,]
            
    def MakeRBThermals(rbKey,phaseVary,phaseDict):
        rbid = str(rbids.index(RB['RBId']))
        tlstr = ['11','22','33','12','13','23']
        sstr = ['12','13','21','23','31','32','AA','BB']
        if 'T' in RB['ThermalMotion'][0]:
            pfxRB = pfx+'RB'+rbKey+'T'
            for i in range(6):
                name = pfxRB+tlstr[i]+':'+str(iRB)+':'+rbid
                phaseDict[name] = RB['ThermalMotion'][1][i]
                if RB['ThermalMotion'][2][i]:
                    phaseVary += [name,]
        if 'L' in RB['ThermalMotion'][0]:
            pfxRB = pfx+'RB'+rbKey+'L'
            for i in range(6):
                name = pfxRB+tlstr[i]+':'+str(iRB)+':'+rbid
                phaseDict[name] = RB['ThermalMotion'][1][i+6]
                if RB['ThermalMotion'][2][i+6]:
                    phaseVary += [name,]
        if 'S' in RB['ThermalMotion'][0]:
            pfxRB = pfx+'RB'+rbKey+'S'
            for i in range(8):
                name = pfxRB+sstr[i]+':'+str(iRB)+':'+rbid
                phaseDict[name] = RB['ThermalMotion'][1][i+12]
                if RB['ThermalMotion'][2][i+12]:
                    phaseVary += [name,]
        if 'U' in RB['ThermalMotion'][0]:
            name = pfx+'RB'+rbKey+'U:'+str(iRB)+':'+rbid
            phaseDict[name] = RB['ThermalMotion'][1][0]
            if RB['ThermalMotion'][2][0]:
                phaseVary += [name,]
                
    def MakeRBTorsions(rbKey,phaseVary,phaseDict):
        rbid = str(rbids.index(RB['RBId']))
        pfxRB = pfx+'RB'+rbKey+'Tr;'
        for i,torsion in enumerate(RB['Torsions']):
            name = pfxRB+str(i)+':'+str(iRB)+':'+rbid
            phaseDict[name] = torsion[0]
            if torsion[1]:
                phaseVary += [name,]
                    
    if Print:
        pFile.write('\n Phases:\n')
    phaseVary = []
    phaseDict = {}
    pawleyLookup = {}
    FFtables = {}                   #scattering factors - xrays
    MFtables = {}                   #Mag. form factors
    BLtables = {}                   # neutrons
    Natoms = {}
    maxSSwave = {}
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    atomIndx = {}
    for name in PhaseData:
        General = PhaseData[name]['General']
        pId = PhaseData[name]['pId']
        pfx = str(pId)+'::'
        FFtable = G2el.GetFFtable(General['AtomTypes'])
        BLtable = G2el.GetBLtable(General)
        FFtables.update(FFtable)
        BLtables.update(BLtable)
        phaseDict[pfx+'isMag'] = False
        SGData = General['SGData']
        SGtext,SGtable = G2spc.SGPrint(SGData)
        if General['Type'] == 'magnetic':
            MFtable = G2el.GetMFtable(General['AtomTypes'],General['Lande g'])
            MFtables.update(MFtable)
            phaseDict[pfx+'isMag'] = True
            SpnFlp = SGData['SpnFlp']
        Atoms = PhaseData[name]['Atoms']
        if Atoms and not General.get('doPawley'):
            cx,ct,cs,cia = General['AtomPtrs']
            AtLookup = G2mth.FillAtomLookUp(Atoms,cia+8)
        PawleyRef = PhaseData[name].get('Pawley ref',[])
        cell = General['Cell']
        A = G2lat.cell2A(cell[1:7])
        phaseDict.update({pfx+'A0':A[0],pfx+'A1':A[1],pfx+'A2':A[2],
            pfx+'A3':A[3],pfx+'A4':A[4],pfx+'A5':A[5],pfx+'Vol':G2lat.calc_V(A)})
        if cell[0]:
            phaseVary += cellVary(pfx,SGData)       #also fills in symmetry required constraints 
        SSGtext = []    #no superstructure
        im = 0
        if General.get('Modulated',False):
            im = 1      #refl offset
            Vec,vRef,maxH = General['SuperVec']
            phaseDict.update({pfx+'mV0':Vec[0],pfx+'mV1':Vec[1],pfx+'mV2':Vec[2]})
            SSGData = General['SSGData']
            SSGtext,SSGtable = G2spc.SSGPrint(SGData,SSGData)
            if vRef:
                phaseVary += modVary(pfx,SSGData)        
        resRBData = PhaseData[name]['RBModels'].get('Residue',[])
        if resRBData:
            rbids = rbIds['Residue']    #NB: used in the MakeRB routines
            for iRB,RB in enumerate(resRBData):
                MakeRBParms('R',phaseVary,phaseDict)
                MakeRBThermals('R',phaseVary,phaseDict)
                MakeRBTorsions('R',phaseVary,phaseDict)
        
        vecRBData = PhaseData[name]['RBModels'].get('Vector',[])
        if vecRBData:
            rbids = rbIds['Vector']    #NB: used in the MakeRB routines
            for iRB,RB in enumerate(vecRBData):
                MakeRBParms('V',phaseVary,phaseDict)
                MakeRBThermals('V',phaseVary,phaseDict)
                    
        Natoms[pfx] = 0
        maxSSwave[pfx] = {'Sfrac':0,'Spos':0,'Sadp':0,'Smag':0}
        if Atoms and not General.get('doPawley'):
            cx,ct,cs,cia = General['AtomPtrs']
            Natoms[pfx] = len(Atoms)
            for i,at in enumerate(Atoms):
                atomIndx[at[cia+8]] = [pfx,i]      #lookup table for restraints
                phaseDict.update({pfx+'Atype:'+str(i):at[ct],pfx+'Afrac:'+str(i):at[cx+3],pfx+'Amul:'+str(i):at[cs+1],
                    pfx+'Ax:'+str(i):at[cx],pfx+'Ay:'+str(i):at[cx+1],pfx+'Az:'+str(i):at[cx+2],
                    pfx+'dAx:'+str(i):0.,pfx+'dAy:'+str(i):0.,pfx+'dAz:'+str(i):0.,         #refined shifts for x,y,z
                    pfx+'AI/A:'+str(i):at[cia],})
                if at[cia] == 'I':
                    phaseDict[pfx+'AUiso:'+str(i)] = at[cia+1]
                else:
                    phaseDict.update({pfx+'AU11:'+str(i):at[cia+2],pfx+'AU22:'+str(i):at[cia+3],pfx+'AU33:'+str(i):at[cia+4],
                        pfx+'AU12:'+str(i):at[cia+5],pfx+'AU13:'+str(i):at[cia+6],pfx+'AU23:'+str(i):at[cia+7]})
                if General['Type'] == 'magnetic':
                    phaseDict.update({pfx+'AMx:'+str(i):at[cx+4],pfx+'AMy:'+str(i):at[cx+5],pfx+'AMz:'+str(i):at[cx+6]})
                if 'F' in at[ct+1]:
                    phaseVary.append(pfx+'Afrac:'+str(i))
                if 'X' in at[ct+1]:
                    try:    #patch for sytsym name changes
                        xId,xCoef = G2spc.GetCSxinel(at[cs])
                    except KeyError:
                        Sytsym = G2spc.SytSym(at[cx:cx+3],SGData)[0]
                        at[cs] = Sytsym
                        xId,xCoef = G2spc.GetCSxinel(at[cs])
                    names = [pfx+'dAx:'+str(i),pfx+'dAy:'+str(i),pfx+'dAz:'+str(i)]
                    equivs = [[],[],[]]
                    for j in range(3):
                        if xId[j] > 0:                               
                            phaseVary.append(names[j])
                            equivs[xId[j]-1].append([names[j],xCoef[j]])
                    for equiv in equivs:
                        if len(equiv) > 1:
                            name = equiv[0][0]
                            coef = equiv[0][1]
                            for eqv in equiv[1:]:
                                eqv[1] /= coef
                                G2mv.StoreEquivalence(name,(eqv,))
                if 'U' in at[ct+1]:
                    if at[cia] == 'I':
                        phaseVary.append(pfx+'AUiso:'+str(i))
                    else:
                        try:    #patch for sytsym name changes
                            uId,uCoef = G2spc.GetCSuinel(at[cs])[:2]
                        except KeyError:
                            Sytsym = G2spc.SytSym(at[cx:cx+3],SGData)[0]
                            at[cs] = Sytsym
                            uId,uCoef = G2spc.GetCSuinel(at[cs])[:2]
                        names = [pfx+'AU11:'+str(i),pfx+'AU22:'+str(i),pfx+'AU33:'+str(i),
                            pfx+'AU12:'+str(i),pfx+'AU13:'+str(i),pfx+'AU23:'+str(i)]
                        equivs = [[],[],[],[],[],[]]
                        for j in range(6):
                            if uId[j] > 0:                               
                                phaseVary.append(names[j])
                                equivs[uId[j]-1].append([names[j],uCoef[j]])
                        for equiv in equivs:
                            if len(equiv) > 1:
                                name = equiv[0][0]
                                coef = equiv[0][1]
                                for eqv in equiv[1:]:
                                    eqv[1] /= coef
                                    G2mv.StoreEquivalence(name,(eqv,))
                if 'M' in at[ct+1]:
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(at[cx:cx+3],SGData)
                    mId,mCoef = G2spc.GetCSpqinel(SytSym,SpnFlp,dupDir)
                    names = [pfx+'AMx:'+str(i),pfx+'AMy:'+str(i),pfx+'AMz:'+str(i)]
                    equivs = [[],[],[]]
                    for j in range(3):
                        if mId[j] > 0:
                            phaseVary.append(names[j])
                            equivs[mId[j]-1].append([names[j],mCoef[j]])
                    for equiv in equivs:
                        if len(equiv) > 1:
                            name = equiv[0][0]
                            coef = equiv[0][1]
                            for eqv in equiv[1:]:
                                eqv[1] /= coef
                                G2mv.StoreEquivalence(name,(eqv,))
                if General.get('Modulated',False):
                    AtomSS = at[-1]['SS1']
                    waveType = AtomSS['waveType']
                    phaseDict[pfx+'waveType:'+str(i)] = waveType
                    for Stype in ['Sfrac','Spos','Sadp','Smag']:
                        Waves = AtomSS[Stype]
                        nx = 0
                        for iw,wave in enumerate(Waves):
                            if not iw:
                                if waveType in ['ZigZag','Block']:
                                    nx = 1
                                CSI = G2spc.GetSSfxuinel(waveType,1,at[cx:cx+3],SGData,SSGData)
                            else:
                                CSI = G2spc.GetSSfxuinel('Fourier',iw+1-nx,at[cx:cx+3],SGData,SSGData)
                            uId,uCoef = CSI[Stype]
                            stiw = str(i)+':'+str(iw)
                            if Stype == 'Spos':
                                if waveType in ['ZigZag','Block',] and not iw:
                                    names = [pfx+'Tmin:'+stiw,pfx+'Tmax:'+stiw,pfx+'Xmax:'+stiw,pfx+'Ymax:'+stiw,pfx+'Zmax:'+stiw]
                                    equivs = [[],[], [],[],[]]
                                else:
                                    names = [pfx+'Xsin:'+stiw,pfx+'Ysin:'+stiw,pfx+'Zsin:'+stiw,
                                        pfx+'Xcos:'+stiw,pfx+'Ycos:'+stiw,pfx+'Zcos:'+stiw]
                                    equivs = [[],[],[], [],[],[]]
                            elif Stype == 'Sadp':
                                names = [pfx+'U11sin:'+stiw,pfx+'U22sin:'+stiw,pfx+'U33sin:'+stiw,
                                    pfx+'U12sin:'+stiw,pfx+'U13sin:'+stiw,pfx+'U23sin:'+stiw,
                                    pfx+'U11cos:'+stiw,pfx+'U22cos:'+stiw,pfx+'U33cos:'+stiw,
                                    pfx+'U12cos:'+stiw,pfx+'U13cos:'+stiw,pfx+'U23cos:'+stiw]
                                equivs = [[],[],[],[],[],[], [],[],[],[],[],[]]
                            elif Stype == 'Sfrac':
                                equivs = [[],[]]
                                if 'Crenel' in waveType and not iw:
                                    names = [pfx+'Fzero:'+stiw,pfx+'Fwid:'+stiw]
                                else:
                                    names = [pfx+'Fsin:'+stiw,pfx+'Fcos:'+stiw]
                            elif Stype == 'Smag':
                                equivs = [[],[],[], [],[],[]]
                                names = [pfx+'MXsin:'+stiw,pfx+'MYsin:'+stiw,pfx+'MZsin:'+stiw,
                                    pfx+'MXcos:'+stiw,pfx+'MYcos:'+stiw,pfx+'MZcos:'+stiw]
                            phaseDict.update(dict(zip(names,wave[0])))
                            if wave[1]: #what do we do here for multiple terms in modulation constraints?
                                for j in range(len(equivs)):
                                    if uId[j][0] > 0:                               
                                        phaseVary.append(names[j])
                                        equivs[uId[j][0]-1].append([names[j],uCoef[j][0]])
                                for equiv in equivs:
                                    if len(equiv) > 1:
                                        name = equiv[0][0]
                                        coef = equiv[0][1]
                                        for eqv in equiv[1:]:
                                            eqv[1] /= coef
                                            G2mv.StoreEquivalence(name,(eqv,))
                            maxSSwave[pfx][Stype] = max(maxSSwave[pfx][Stype],iw+1)
            textureData = General['SH Texture']
            if textureData['Order'] and not seqRef:
                phaseDict[pfx+'SHorder'] = textureData['Order']
                phaseDict[pfx+'SHmodel'] = SamSym[textureData['Model']]
                for item in ['omega','chi','phi']:
                    phaseDict[pfx+'SH '+item] = textureData['Sample '+item][1]
                    if textureData['Sample '+item][0]:
                        phaseVary.append(pfx+'SH '+item)
                for item in textureData['SH Coeff'][1]:
                    phaseDict[pfx+item] = textureData['SH Coeff'][1][item]
                    if textureData['SH Coeff'][0]:
                        phaseVary.append(pfx+item)
                
            if Print:
                pFile.write('\n Phase name: %s\n'%General['Name'])
                pFile.write(135*'='+'\n')
                PrintFFtable(FFtable)
                PrintBLtable(BLtable)
                if General['Type'] == 'magnetic':
                    PrintMFtable(MFtable)
                pFile.write('\n')
                #how do we print magnetic symmetry table? TBD
                if len(SSGtext):    #if superstructure
                    for line in SSGtext: pFile.write(line+'\n')
                    if len(SSGtable):
                        for item in SSGtable:
                            line = ' %s '%(item)
                            pFile.write(line+'\n') 
                    else:
                        pFile.write(' ( 1)    %s\n'%(SSGtable[0]))
                else:
                    for line in SGtext: pFile.write(line+'\n')
                    if len(SGtable):
                        for item in SGtable:
                            line = ' %s '%(item)
                            pFile.write(line+'\n') 
                    else:
                        pFile.write(' ( 1)    %s\n'%(SGtable[0]))
                PrintRBObjects(resRBData,vecRBData)
                PrintAtoms(General,Atoms)
                if General['Type'] == 'magnetic':
                    PrintMoments(General,Atoms)
                if General.get('Modulated',False):
                    PrintWaves(General,Atoms)
                pFile.write('\n Unit cell: a = %.5f b = %.5f c = %.5f alpha = %.3f beta = %.3f gamma = %.3f volume = %.3f Refine? %s\n'%
                    (cell[1],cell[2],cell[3],cell[4],cell[5],cell[6],cell[7],cell[0]))
                if len(SSGtext):    #if superstructure
                    pFile.write('\n Modulation vector: mV0 = %.4f mV1 = %.4f mV2 = %.4f max mod. index = %d Refine? %s\n'%
                        (Vec[0],Vec[1],Vec[2],maxH,vRef))
                if not seqRef:
                    PrintTexture(textureData)
                if name in RestraintDict:
                    PrintRestraints(cell[1:7],SGData,General['AtomPtrs'],Atoms,AtLookup,
                        textureData,RestraintDict[name],pFile)
                    
        elif PawleyRef:
            if Print:
                pFile.write('\n Phase name: %s\n'%General['Name'])
                pFile.write(135*'='+'\n')
                pFile.write('\n')
                if len(SSGtext):    #if superstructure
                    for line in SSGtext: pFile.write(line+'\n')
                    if len(SSGtable):
                        for item in SSGtable:
                            line = ' %s '%(item)
                            pFile.write(line+'\n')  
                    else:
                        pFile.write(' ( 1)    %s\n'%SSGtable[0])
                else:
                    for line in SGtext: pFile.write(line+'\n')
                    if len(SGtable):
                        for item in SGtable:
                            line = ' %s '%(item)
                            pFile.write(line+'\n')
                    else:
                        pFile.write(' ( 1)    %s\n'%(SGtable[0]))
                pFile.write('\n Unit cell: a = %.5f b = %.5f c = %.5f alpha = %.3f beta = %.3f gamma = %.3f volume = %.3f Refine? %s\n'%
                    (cell[1],cell[2],cell[3],cell[4],cell[5],cell[6],cell[7],cell[0]))
                if len(SSGtext):    #if superstructure
                    pFile.write('\n Modulation vector: mV0 = %.4f mV1 = %.4f mV2 = %.4f max mod. index = %d Refine? %s\n'%
                        (Vec[0],Vec[1],Vec[2],maxH,vRef))
            pawleyVary = []
            for i,refl in enumerate(PawleyRef):
                phaseDict[pfx+'PWLref:'+str(i)] = refl[6+im]
                if im:
                    pawleyLookup[pfx+'%d,%d,%d,%d'%(refl[0],refl[1],refl[2],refl[3])] = i
                else:
                    pawleyLookup[pfx+'%d,%d,%d'%(refl[0],refl[1],refl[2])] = i
                if refl[5+im]:
                    pawleyVary.append(pfx+'PWLref:'+str(i))
            GetPawleyConstr(SGData['SGLaue'],PawleyRef,im,pawleyVary)      #does G2mv.StoreEquivalence
            phaseVary += pawleyVary
                
    return Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,maxSSwave
    
def cellFill(pfx,SGData,parmDict,sigDict): 
    '''Returns the filled-out reciprocal cell (A) terms and their uncertainties
    from the parameter and sig dictionaries.

    :param str pfx: parameter prefix ("n::", where n is a phase number)
    :param dict SGdata: a symmetry object
    :param dict parmDict: a dictionary of parameters
    :param dict sigDict:  a dictionary of uncertainties on parameters

    :returns: A,sigA where each is a list of six terms with the A terms 
    '''
    if SGData['SGLaue'] in ['-1',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
            parmDict[pfx+'A3'],parmDict[pfx+'A4'],parmDict[pfx+'A5']]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,0,parmDict[pfx+'A5']]
        elif SGData['SGUniq'] == 'b':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,parmDict[pfx+'A4'],0]
        else:
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                parmDict[pfx+'A3'],0,0]
    elif SGData['SGLaue'] in ['mmm',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],
            parmDict[pfx+'A0'],0,0]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],
            parmDict[pfx+'A3'],parmDict[pfx+'A3'],parmDict[pfx+'A3']]
    elif SGData['SGLaue'] in ['m3m','m3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],0,0,0]

    try:
        if SGData['SGLaue'] in ['-1',]:
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                sigDict[pfx+'A3'],sigDict[pfx+'A4'],sigDict[pfx+'A5']]
        elif SGData['SGLaue'] in ['2/m',]:
            if SGData['SGUniq'] == 'a':
                sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                    0,0,sigDict[pfx+'A5']]
            elif SGData['SGUniq'] == 'b':
                sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                    0,sigDict[pfx+'A4'],0]
            else:
                sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                    sigDict[pfx+'A3'],0,0]
        elif SGData['SGLaue'] in ['mmm',]:
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],0,0,0]
        elif SGData['SGLaue'] in ['4/m','4/mmm']:
            sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
        elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
            sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
        elif SGData['SGLaue'] in ['3R', '3mR']:
            sigA = [sigDict[pfx+'A0'],0,0,sigDict[pfx+'A3'],0,0]
        elif SGData['SGLaue'] in ['m3m','m3']:
            sigA = [sigDict[pfx+'A0'],0,0,0,0,0]
    except KeyError:
        sigA = [0,0,0,0,0,0]
    return A,sigA
        
def PrintRestraints(cell,SGData,AtPtrs,Atoms,AtLookup,textureData,phaseRest,pFile):
    'needs a doc string'
    if phaseRest:
        Amat = G2lat.cell2AB(cell)[0]
        cx,ct,cs = AtPtrs[:3]
        names = [['Bond','Bonds'],['Angle','Angles'],['Plane','Planes'],
            ['Chiral','Volumes'],['Torsion','Torsions'],['Rama','Ramas'],
            ['ChemComp','Sites'],['Texture','HKLs']]
        for name,rest in names:
            itemRest = phaseRest[name]
            if rest in itemRest and itemRest[rest] and itemRest['Use']:
                pFile.write('\n  %s restraint weight factor %10.3f Use: %s\n'%(name,itemRest['wtFactor'],str(itemRest['Use'])))
                if name in ['Bond','Angle','Plane','Chiral']:
                    pFile.write('     calc       obs      sig   delt/sig  atoms(symOp): \n')
                    for indx,ops,obs,esd in itemRest[rest]:
                        try:
                            AtNames = G2mth.GetAtomItemsById(Atoms,AtLookup,indx,ct-1)
                            AtName = ''
                            for i,Aname in enumerate(AtNames):
                                AtName += Aname
                                if ops[i] == '1':
                                    AtName += '-'
                                else:
                                    AtName += '+('+ops[i]+')-'
                            XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cx,3))
                            XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                            if name == 'Bond':
                                calc = G2mth.getRestDist(XYZ,Amat)
                            elif name == 'Angle':
                                calc = G2mth.getRestAngle(XYZ,Amat)
                            elif name == 'Plane':
                                calc = G2mth.getRestPlane(XYZ,Amat)
                            elif name == 'Chiral':
                                calc = G2mth.getRestChiral(XYZ,Amat)
                            pFile.write(' %9.3f %9.3f %8.3f %8.3f   %s\n'%(calc,obs,esd,(obs-calc)/esd,AtName[:-1]))
                        except KeyError:
                            del itemRest[rest]
                elif name in ['Torsion','Rama']:
                    pFile.write('  atoms(symOp)  calc  obs  sig  delt/sig  torsions: \n')
                    coeffDict = itemRest['Coeff']
                    for indx,ops,cofName,esd in itemRest[rest]:
                        AtNames = G2mth.GetAtomItemsById(Atoms,AtLookup,indx,ct-1)
                        AtName = ''
                        for i,Aname in enumerate(AtNames):
                            AtName += Aname+'+('+ops[i]+')-'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        if name == 'Torsion':
                            tor = G2mth.getRestTorsion(XYZ,Amat)
                            restr,calc = G2mth.calcTorsionEnergy(tor,coeffDict[cofName])
                            pFile.write(' %8.3f %8.3f %.3f %8.3f %8.3f %s\n'%(calc,obs,esd,(obs-calc)/esd,tor,AtName[:-1]))
                        else:
                            phi,psi = G2mth.getRestRama(XYZ,Amat)
                            restr,calc = G2mth.calcRamaEnergy(phi,psi,coeffDict[cofName])                               
                            pFile.write(' %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %s\n'%(calc,obs,esd,(obs-calc)/esd,phi,psi,AtName[:-1]))
                elif name == 'ChemComp':
                    pFile.write('     atoms   mul*frac  factor     prod\n')
                    for indx,factors,obs,esd in itemRest[rest]:
                        try:
                            atoms = G2mth.GetAtomItemsById(Atoms,AtLookup,indx,ct-1)
                            mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs+1))
                            frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs-1))
                            mulfrac = mul*frac
                            calcs = mul*frac*factors
                            for iatm,[atom,mf,fr,clc] in enumerate(zip(atoms,mulfrac,factors,calcs)):
                                pFile.write(' %10s %8.3f %8.3f %8.3f\n'%(atom,mf,fr,clc))
                            pFile.write(' Sum:                   calc: %8.3f obs: %8.3f esd: %8.3f\n'%(np.sum(calcs),obs,esd))
                        except KeyError:
                            del itemRest[rest]
                elif name == 'Texture' and textureData['Order']:
                    Start = False
                    SHCoef = textureData['SH Coeff'][1]
                    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
                    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
                    pFile.write ('    HKL  grid  neg esd  sum neg  num neg use unit?  unit esd \n')
                    for hkl,grid,esd1,ifesd2,esd2 in itemRest[rest]:
                        phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                        ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
                        R,P,Z = G2mth.getRestPolefig(ODFln,SamSym[textureData['Model']],grid)
                        Z = ma.masked_greater(Z,0.0)
                        num = ma.count(Z)
                        sum = 0
                        if num:
                            sum = np.sum(Z)
                        pFile.write ('   %d %d %d  %d %8.3f %8.3f %8d   %s    %8.3f\n'%(hkl[0],hkl[1],hkl[2],grid,esd1,sum,num,str(ifesd2),esd2))
        
def getCellEsd(pfx,SGData,A,covData):
    'needs a doc string'
    rVsq = G2lat.calc_rVsq(A)
    G,g = G2lat.A2Gmat(A)       #get recip. & real metric tensors
    RMnames = [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    varyList = covData['varyList']
    covMatrix = covData['covMatrix']
    vcov = G2mth.getVCov(RMnames,varyList,covMatrix)
    if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm']:
        vcov[1,1] = vcov[3,3] = vcov[0,1] = vcov[1,0] = vcov[0,0]
        vcov[1,3] = vcov[3,1] = vcov[0,3] = vcov[3,0] = vcov[0,0]
        vcov[1,2] = vcov[2,1] = vcov[2,3] = vcov[3,2] = vcov[0,2]
    elif SGData['SGLaue'] in ['m3','m3m']:
        vcov[0:3,0:3] = vcov[0,0]
    elif SGData['SGLaue'] in ['4/m', '4/mmm']:
        vcov[0:2,0:2] = vcov[0,0]
        vcov[1,2] = vcov[2,1] = vcov[0,2]
    elif SGData['SGLaue'] in ['3R','3mR']:
        vcov[0:3,0:3] = vcov[0,0]
#        vcov[4,4] = vcov[5,5] = vcov[3,3]
        vcov[3:6,3:6] = vcov[3,3]
        vcov[0:3,3:6] = vcov[0,3]
        vcov[3:6,0:3] = vcov[3,0]
    delt = 1.e-9
    drVdA = np.zeros(6)
    for i in range(6):
        A[i] += delt
        drVdA[i] = G2lat.calc_rVsq(A)
        A[i] -= 2*delt
        drVdA[i] -= G2lat.calc_rVsq(A)
        A[i] += delt
    drVdA /= 2.*delt    
    srcvlsq = np.inner(drVdA,np.inner(drVdA,vcov))
    Vol = 1/np.sqrt(rVsq)
    sigVol = Vol**3*np.sqrt(srcvlsq)/2.         #ok - checks with GSAS
    
    dcdA = np.zeros((6,6))
    for i in range(6):
        pdcdA =np.zeros(6)
        A[i] += delt
        pdcdA += G2lat.A2cell(A)
        A[i] -= 2*delt
        pdcdA -= G2lat.A2cell(A)
        A[i] += delt
        dcdA[i] = pdcdA/(2.*delt)
    
    sigMat = np.inner(dcdA,np.inner(dcdA,vcov))
    var = np.diag(sigMat)
    CS = np.where(var>0.,np.sqrt(var),0.)
    if SGData['SGLaue'] in ['3', '3m1', '31m', '6/m', '6/mmm','m3','m3m','4/m','4/mmm']:
        CS[3:6] = 0.0
    return [CS[0],CS[1],CS[2],CS[5],CS[4],CS[3],sigVol]
    
def SetPhaseData(parmDict,sigDict,Phases,RBIds,covData,RestraintDict=None,pFile=None):
    'needs a doc string'
    
    def PrintAtomsAndSig(General,Atoms,atomsSig):
        pFile.write('\n Atoms:\n')
        line = '   name      x         y         z      frac   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'macromolecular':
            line = ' res no residue chain '+line
        cx,ct,cs,cia = General['AtomPtrs']
        pFile.write(line+'\n')
        pFile.write(135*'-'+'\n')
        fmt = {0:'%7s',ct:'%7s',cx:'%10.5f',cx+1:'%10.5f',cx+2:'%10.5f',cx+3:'%8.3f',cia+1:'%8.5f',
            cia+2:'%8.5f',cia+3:'%8.5f',cia+4:'%8.5f',cia+5:'%8.5f',cia+6:'%8.5f',cia+7:'%8.5f'}
        noFXsig = {cx:[10*' ','%10s'],cx+1:[10*' ','%10s'],cx+2:[10*' ','%10s'],cx+3:[8*' ','%8s']}
        for atyp in General['AtomTypes']:       #zero composition
            General['NoAtoms'][atyp] = 0.
        for i,at in enumerate(Atoms):
            General['NoAtoms'][at[ct]] += at[cx+3]*float(at[cx+5])     #new composition
            if General['Type'] == 'macromolecular':
                name = ' %s %s %s %s:'%(at[0],at[1],at[2],at[3])
                valstr = ' values:          '
                sigstr = ' sig   :          '
            else:
                name = fmt[0]%(at[ct-1])+fmt[1]%(at[ct])+':'
                valstr = ' values:'
                sigstr = ' sig   :'
            for ind in range(cx,cx+4):
                sigind = str(i)+':'+str(ind)
                valstr += fmt[ind]%(at[ind])                    
                if sigind in atomsSig:
                    sigstr += fmt[ind]%(atomsSig[sigind])
                else:
                    sigstr += noFXsig[ind][1]%(noFXsig[ind][0])
            if at[cia] == 'I':
                valstr += fmt[cia+1]%(at[cia+1])
                if '%d:%d'%(i,cia+1) in atomsSig:
                    sigstr += fmt[cia+1]%(atomsSig['%d:%d'%(i,cia+1)])
                else:
                    sigstr += 8*' '
            else:
                valstr += 8*' '
                sigstr += 8*' '
                for ind in range(cia+2,cia+8):
                    sigind = str(i)+':'+str(ind)
                    valstr += fmt[ind]%(at[ind])
                    if sigind in atomsSig:                        
                        sigstr += fmt[ind]%(atomsSig[sigind])
                    else:
                        sigstr += 8*' '
            pFile.write(name+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
            
    def PrintMomentsAndSig(General,Atoms,atomsSig):
        pFile.write('\n Magnetic Moments:\n')    #add magnitude & angle, etc.? TBD
        line = '   name      Mx        My        Mz'
        cx,ct,cs,cia = General['AtomPtrs']
        cmx = cx+4
        AtInfo = dict(zip(General['AtomTypes'],General['Lande g']))
        pFile.write(line+'\n')
        pFile.write(135*'-'+'\n')
        fmt = {0:'%7s',ct:'%7s',cmx:'%10.5f',cmx+1:'%10.5f',cmx+2:'%10.5f'}
        noFXsig = {cmx:[10*' ','%10s'],cmx+1:[10*' ','%10s'],cmx+2:[10*' ','%10s']}
        for i,at in enumerate(Atoms):
            if AtInfo[at[ct]]:
                name = fmt[0]%(at[ct-1])+fmt[1]%(at[ct])+':'
                valstr = ' values:'
                sigstr = ' sig   :'
                for ind in range(cmx,cmx+3):
                    sigind = str(i)+':'+str(ind)
                    valstr += fmt[ind]%(at[ind])                    
                    if sigind in atomsSig:
                        sigstr += fmt[ind]%(atomsSig[sigind])
                    else:
                        sigstr += noFXsig[ind][1]%(noFXsig[ind][0])
                pFile.write(name+'\n')
                pFile.write(valstr+'\n')
                pFile.write(sigstr+'\n')
            
    def PrintWavesAndSig(General,Atoms,wavesSig):
        cx,ct,cs,cia = General['AtomPtrs']
        pFile.write('\n Modulation waves\n')
        names = {'Sfrac':['Fsin','Fcos','Fzero','Fwid'],'Spos':['Xsin','Ysin','Zsin','Xcos','Ycos','Zcos','Tmin','Tmax','Xmax','Ymax','Zmax'],
            'Sadp':['U11sin','U22sin','U33sin','U12sin','U13sin','U23sin','U11cos','U22cos',
            'U33cos','U12cos','U13cos','U23cos'],'Smag':['MXsin','MYsin','MZsin','MXcos','MYcos','MZcos']}
        pFile.write(135*'-'+'\n')
        for i,at in enumerate(Atoms):
            AtomSS = at[-1]['SS1']
            waveType = AtomSS['waveType']
            for Stype in ['Sfrac','Spos','Sadp','Smag']:
                Waves = AtomSS[Stype]
                if len(Waves):
                    pFile.write(' atom: %s, site sym: %s, type: %s wave type: %s:\n'%
                        (at[ct-1],at[cs],Stype,waveType))
                    for iw,wave in enumerate(Waves):
                        stiw = ':'+str(i)+':'+str(iw)
                        namstr = '  names :'
                        valstr = '  values:'
                        sigstr = '  esds  :'
                        if Stype == 'Spos':
                            nt = 6
                            ot = 0
                            if waveType in ['ZigZag','Block',] and not iw:
                                nt = 5
                                ot = 6
                            for j in range(nt):
                                name = names['Spos'][j+ot]
                                namstr += '%12s'%(name)
                                valstr += '%12.4f'%(wave[0][j])
                                if name+stiw in wavesSig:
                                    sigstr += '%12.4f'%(wavesSig[name+stiw])
                                else:
                                    sigstr += 12*' '
                        elif Stype == 'Sfrac':
                            ot = 0
                            if 'Crenel' in waveType and not iw:
                                ot = 2
                            for j in range(2):
                                name = names['Sfrac'][j+ot]
                                namstr += '%12s'%(names['Sfrac'][j+ot])
                                valstr += '%12.4f'%(wave[0][j])
                                if name+stiw in wavesSig:
                                    sigstr += '%12.4f'%(wavesSig[name+stiw])
                                else:
                                    sigstr += 12*' '
                        elif Stype == 'Sadp':
                            for j in range(12):
                                name = names['Sadp'][j]
                                namstr += '%10s'%(names['Sadp'][j])
                                valstr += '%10.6f'%(wave[0][j])
                                if name+stiw in wavesSig:
                                    sigstr += '%10.6f'%(wavesSig[name+stiw])
                                else:
                                    sigstr += 10*' '
                        elif Stype == 'Smag':
                            for j in range(6):
                                name = names['Smag'][j]
                                namstr += '%12s'%(names['Smag'][j])
                                valstr += '%12.4f'%(wave[0][j])
                                if name+stiw in wavesSig:
                                    sigstr += '%12.4f'%(wavesSig[name+stiw])
                                else:
                                    sigstr += 12*' '
                                
                    pFile.write(namstr+'\n')
                    pFile.write(valstr+'\n')
                    pFile.write(sigstr+'\n')
        
                
    def PrintRBObjPOAndSig(rbfx,rbsx):
        namstr = '  names :'
        valstr = '  values:'
        sigstr = '  esds  :'
        for i,px in enumerate(['Px:','Py:','Pz:']):
            name = pfx+rbfx+px+rbsx
            namstr += '%12s'%('Pos '+px[1])
            valstr += '%12.5f'%(parmDict[name])
            if name in sigDict:
                sigstr += '%12.5f'%(sigDict[name])
            else:
                sigstr += 12*' '
        for i,po in enumerate(['Oa:','Oi:','Oj:','Ok:']):
            name = pfx+rbfx+po+rbsx
            namstr += '%12s'%('Ori '+po[1])
            valstr += '%12.5f'%(parmDict[name])
            if name in sigDict:
                sigstr += '%12.5f'%(sigDict[name])
            else:
                sigstr += 12*' '
        pFile.write(namstr+'\n')
        pFile.write(valstr+'\n')
        pFile.write(sigstr+'\n')
        
    def PrintRBObjTLSAndSig(rbfx,rbsx,TLS):
        namstr = '  names :'
        valstr = '  values:'
        sigstr = '  esds  :'
        if 'N' not in TLS:
            pFile.write(' Thermal motion:\n')
        if 'T' in TLS:
            for i,pt in enumerate(['T11:','T22:','T33:','T12:','T13:','T23:']):
                name = pfx+rbfx+pt+rbsx
                namstr += '%12s'%(pt[:3])
                valstr += '%12.5f'%(parmDict[name])
                if name in sigDict:
                    sigstr += '%12.5f'%(sigDict[name])
                else:
                    sigstr += 12*' '
            pFile.write(namstr+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
        if 'L' in TLS:
            namstr = '  names :'
            valstr = '  values:'
            sigstr = '  esds  :'
            for i,pt in enumerate(['L11:','L22:','L33:','L12:','L13:','L23:']):
                name = pfx+rbfx+pt+rbsx
                namstr += '%12s'%(pt[:3])
                valstr += '%12.3f'%(parmDict[name])
                if name in sigDict:
                    sigstr += '%12.3f'%(sigDict[name])
                else:
                    sigstr += 12*' '
            pFile.write(namstr+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
        if 'S' in TLS:
            namstr = '  names :'
            valstr = '  values:'
            sigstr = '  esds  :'
            for i,pt in enumerate(['S12:','S13:','S21:','S23:','S31:','S32:','SAA:','SBB:']):
                name = pfx+rbfx+pt+rbsx
                namstr += '%12s'%(pt[:3])
                valstr += '%12.4f'%(parmDict[name])
                if name in sigDict:
                    sigstr += '%12.4f'%(sigDict[name])
                else:
                    sigstr += 12*' '
            pFile.write(namstr+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
        if 'U' in TLS:
            name = pfx+rbfx+'U:'+rbsx
            namstr = '  names :'+'%12s'%('Uiso')
            valstr = '  values:'+'%12.5f'%(parmDict[name])
            if name in sigDict:
                sigstr = '  esds  :'+'%12.5f'%(sigDict[name])
            else:
                sigstr = '  esds  :'+12*' '
            pFile.write(namstr+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
        
    def PrintRBObjTorAndSig(rbsx):
        namstr = '  names :'
        valstr = '  values:'
        sigstr = '  esds  :'
        nTors = len(RBObj['Torsions'])
        if nTors:
            pFile.write(' Torsions:\n')
            for it in range(nTors):
                name = pfx+'RBRTr;'+str(it)+':'+rbsx
                namstr += '%12s'%('Tor'+str(it))
                valstr += '%12.4f'%(parmDict[name])
                if name in sigDict:
                    sigstr += '%12.4f'%(sigDict[name])
            pFile.write(namstr+'\n')
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
                
    def PrintSHtextureAndSig(textureData,SHtextureSig):
        pFile.write('\n Spherical harmonics texture: Order: %d\n'%textureData['Order'])
        names = ['omega','chi','phi']
        namstr = '  names :'
        ptstr =  '  values:'
        sigstr = '  esds  :'
        for name in names:
            namstr += '%12s'%(name)
            ptstr += '%12.3f'%(textureData['Sample '+name][1])
            if 'Sample '+name in SHtextureSig:
                sigstr += '%12.3f'%(SHtextureSig['Sample '+name])
            else:
                sigstr += 12*' '
        pFile.write(namstr+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(sigstr+'\n')
        pFile.write('\n Texture coefficients:\n')
        SHcoeff = textureData['SH Coeff'][1]
        SHkeys = list(SHcoeff.keys())
        nCoeff = len(SHcoeff)
        nBlock = nCoeff//10+1
        iBeg = 0
        iFin = min(iBeg+10,nCoeff)
        for block in range(nBlock):
            namstr = '  names :'
            ptstr =  '  values:'
            sigstr = '  esds  :'
            for name in SHkeys[iBeg:iFin]:
                namstr += '%12s'%(name)
                ptstr += '%12.3f'%(SHcoeff[name])
                if name in SHtextureSig:
                    sigstr += '%12.3f'%(SHtextureSig[name])
                else:
                    sigstr += 12*' '
            pFile.write(namstr+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(sigstr+'\n')
            iBeg += 10
            iFin = min(iBeg+10,nCoeff)
            
    pFile.write('\n Phases:\n')
    for phase in Phases:
        pFile.write(' Result for phase: %s\n'%phase)
        pFile.write(135*'='+'\n')
        Phase = Phases[phase]
        General = Phase['General']
        SGData = General['SGData']
        Atoms = Phase['Atoms']
        AtLookup = []
        if Atoms and not General.get('doPawley'):
            cx,ct,cs,cia = General['AtomPtrs']
            AtLookup = G2mth.FillAtomLookUp(Atoms,cia+8)
        cell = General['Cell']
        pId = Phase['pId']
        pfx = str(pId)+'::'
        if cell[0]:
            A,sigA = cellFill(pfx,SGData,parmDict,sigDict)
            cellSig = getCellEsd(pfx,SGData,A,covData)  #includes sigVol
            pFile.write(' Reciprocal metric tensor: \n')
            ptfmt = "%15.9f"
            names = ['A11','A22','A33','A12','A13','A23']
            namstr = '  names :'
            ptstr =  '  values:'
            sigstr = '  esds  :'
            for name,a,siga in zip(names,A,sigA):
                namstr += '%15s'%(name)
                ptstr += ptfmt%(a)
                if siga:
                    sigstr += ptfmt%(siga)
                else:
                    sigstr += 15*' '
            pFile.write(namstr+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(sigstr+'\n')
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            pFile.write(' New unit cell:\n')
            ptfmt = ["%12.6f","%12.6f","%12.6f","%12.4f","%12.4f","%12.4f","%12.3f"]
            names = ['a','b','c','alpha','beta','gamma','Volume']
            namstr = '  names :'
            ptstr =  '  values:'
            sigstr = '  esds  :'
            for name,fmt,a,siga in zip(names,ptfmt,cell[1:8],cellSig):
                namstr += '%12s'%(name)
                ptstr += fmt%(a)
                if siga:
                    sigstr += fmt%(siga)
                else:
                    sigstr += 12*' '
            pFile.write(namstr+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(sigstr+'\n')
        ik = 6  #for Pawley stuff below
        if General.get('Modulated',False):
            ik = 7
            Vec,vRef,maxH = General['SuperVec']
            if vRef:
                pFile.write(' New modulation vector:\n')
                namstr = '  names :'
                ptstr =  '  values:'
                sigstr = '  esds  :'
                for var in ['mV0','mV1','mV2']:
                    namstr += '%12s'%(pfx+var)
                    ptstr += '%12.6f'%(parmDict[pfx+var])
                    if pfx+var in sigDict:
                        sigstr += '%12.6f'%(sigDict[pfx+var])
                    else:
                        sigstr += 12*' '
                pFile.write(namstr+'\n')
                pFile.write(ptstr+'\n')
                pFile.write(sigstr+'\n')
            
        General['Mass'] = 0.
        if Phase['General'].get('doPawley'):
            pawleyRef = Phase['Pawley ref']
            for i,refl in enumerate(pawleyRef):
                key = pfx+'PWLref:'+str(i)
                refl[ik] = parmDict[key]
                if key in sigDict:
                    refl[ik+1] = sigDict[key]
                else:
                    refl[ik+1] = 0
        else:
            VRBIds = RBIds['Vector']
            RRBIds = RBIds['Residue']
            RBModels = Phase['RBModels']
            for irb,RBObj in enumerate(RBModels.get('Vector',[])):
                jrb = VRBIds.index(RBObj['RBId'])
                rbsx = str(irb)+':'+str(jrb)
                pFile.write(' Vector rigid body parameters:\n')
                PrintRBObjPOAndSig('RBV',rbsx)
                PrintRBObjTLSAndSig('RBV',rbsx,RBObj['ThermalMotion'][0])
            for irb,RBObj in enumerate(RBModels.get('Residue',[])):
                jrb = RRBIds.index(RBObj['RBId'])
                rbsx = str(irb)+':'+str(jrb)
                pFile.write(' Residue rigid body parameters:\n')
                PrintRBObjPOAndSig('RBR',rbsx)
                PrintRBObjTLSAndSig('RBR',rbsx,RBObj['ThermalMotion'][0])
                PrintRBObjTorAndSig(rbsx)
            atomsSig = {}
            wavesSig = {}
            cx,ct,cs,cia = General['AtomPtrs']
            for i,at in enumerate(Atoms):
                names = {cx:pfx+'Ax:'+str(i),cx+1:pfx+'Ay:'+str(i),cx+2:pfx+'Az:'+str(i),cx+3:pfx+'Afrac:'+str(i),
                    cia+1:pfx+'AUiso:'+str(i),cia+2:pfx+'AU11:'+str(i),cia+3:pfx+'AU22:'+str(i),cia+4:pfx+'AU33:'+str(i),
                    cia+5:pfx+'AU12:'+str(i),cia+6:pfx+'AU13:'+str(i),cia+7:pfx+'AU23:'+str(i),
                    cx+4:pfx+'AMx:'+str(i),cx+5:pfx+'AMy:'+str(i),cx+6:pfx+'AMz:'+str(i)}
                for ind in range(cx,cx+4):
                    at[ind] = parmDict[names[ind]]
                    if ind in range(cx,cx+3):
                        name = names[ind].replace('A','dA')
                    else:
                        name = names[ind]
                    if name in sigDict:
                        atomsSig[str(i)+':'+str(ind)] = sigDict[name]
                if at[cia] == 'I':
                    at[cia+1] = parmDict[names[cia+1]]
                    if names[cia+1] in sigDict:
                        atomsSig['%d:%d'%(i,cia+1)] = sigDict[names[cia+1]]
                else:
                    for ind in range(cia+2,cia+8):
                        at[ind] = parmDict[names[ind]]
                        if names[ind] in sigDict:
                            atomsSig[str(i)+':'+str(ind)] = sigDict[names[ind]]
                if General['Type'] == 'magnetic':
                    for ind in range(cx+4,cx+7):
                        at[ind] = parmDict[names[ind]]
                        if names[ind] in sigDict:
                            atomsSig[str(i)+':'+str(ind)] = sigDict[names[ind]]
                ind = General['AtomTypes'].index(at[ct])
                General['Mass'] += General['AtomMass'][ind]*at[cx+3]*at[cx+5]
                if General.get('Modulated',False):
                    AtomSS = at[-1]['SS1']
                    waveType = AtomSS['waveType']
                    for Stype in ['Sfrac','Spos','Sadp','Smag']:
                        Waves = AtomSS[Stype]
                        for iw,wave in enumerate(Waves):
                            stiw = str(i)+':'+str(iw)
                            if Stype == 'Spos':
                                if waveType in ['ZigZag','Block',] and not iw:
                                    names = ['Tmin:'+stiw,'Tmax:'+stiw,'Xmax:'+stiw,'Ymax:'+stiw,'Zmax:'+stiw]
                                else:
                                    names = ['Xsin:'+stiw,'Ysin:'+stiw,'Zsin:'+stiw,
                                        'Xcos:'+stiw,'Ycos:'+stiw,'Zcos:'+stiw]
                            elif Stype == 'Sadp':
                                names = ['U11sin:'+stiw,'U22sin:'+stiw,'U33sin:'+stiw,
                                    'U12sin:'+stiw,'U13sin:'+stiw,'U23sin:'+stiw,
                                    'U11cos:'+stiw,'U22cos:'+stiw,'U33cos:'+stiw,
                                    'U12cos:'+stiw,'U13cos:'+stiw,'U23cos:'+stiw]
                            elif Stype == 'Sfrac':
                                if 'Crenel' in waveType and not iw:
                                    names = ['Fzero:'+stiw,'Fwid:'+stiw]
                                else:
                                    names = ['Fsin:'+stiw,'Fcos:'+stiw]
                            elif Stype == 'Smag':
                                names = ['MXsin:'+stiw,'MYsin:'+stiw,'MZsin:'+stiw,
                                    'MXcos:'+stiw,'MYcos:'+stiw,'MZcos:'+stiw]
                            for iname,name in enumerate(names):
                                AtomSS[Stype][iw][0][iname] = parmDict[pfx+name]
                                if pfx+name in sigDict:
                                    wavesSig[name] = sigDict[pfx+name]
                    
            PrintAtomsAndSig(General,Atoms,atomsSig)
            if General['Type'] == 'magnetic':
                PrintMomentsAndSig(General,Atoms,atomsSig)
            if General.get('Modulated',False):
                PrintWavesAndSig(General,Atoms,wavesSig)
            
        
        textureData = General['SH Texture']    
        if textureData['Order']:
            SHtextureSig = {}
            for name in ['omega','chi','phi']:
                aname = pfx+'SH '+name
                textureData['Sample '+name][1] = parmDict[aname]
                if aname in sigDict:
                    SHtextureSig['Sample '+name] = sigDict[aname]
            for name in textureData['SH Coeff'][1]:
                aname = pfx+name
                textureData['SH Coeff'][1][name] = parmDict[aname]
                if aname in sigDict:
                    SHtextureSig[name] = sigDict[aname]
            PrintSHtextureAndSig(textureData,SHtextureSig)
        if phase in RestraintDict and not Phase['General'].get('doPawley'):
            PrintRestraints(cell[1:7],SGData,General['AtomPtrs'],Atoms,AtLookup,
                textureData,RestraintDict[phase],pFile)
                    
################################################################################
##### Histogram & Phase data
################################################################################        
                    
def GetHistogramPhaseData(Phases,Histograms,Print=True,pFile=None,resetRefList=True):
    '''Loads the HAP histogram/phase information into dicts

    :param dict Phases: phase information
    :param dict Histograms: Histogram information
    :param bool Print: prints information as it is read
    :param file pFile: file object to print to (the default, None causes printing to the console)
    :param bool resetRefList: Should the contents of the Reflection List be initialized
      on loading. The default, True, initializes the Reflection List as it is loaded.

    :returns: (hapVary,hapDict,controlDict)
      * hapVary: list of refined variables
      * hapDict: dict with refined variables and their values
      * controlDict: dict with fixed parameters
    '''
    
    def PrintSize(hapData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Size model    : %9s'%(hapData[0])
            line += ' equatorial:'+'%12.3f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.3f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            line += '\n\t LG mixing coeff.: %12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            pFile.write(line+'\n')
        else:
            pFile.write('\n Size model    : %s\n\t LG mixing coeff.:%12.4f Refine? %s\n'%
                (hapData[0],hapData[1][2],hapData[2][2]))
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.3f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(varstr+'\n')
        
    def PrintMuStrain(hapData,SGData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Mustrain model: %9s'%(hapData[0])
            line += ' equatorial:'+'%12.1f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.1f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            line +='\n\t LG mixing coeff.:%12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            pFile.write(line+'\n')
        else:
            pFile.write('\n Mustrain model: %s\n\t LG mixing coeff.:%12.4f Refine? %s\n'%
                (hapData[0],hapData[1][2],hapData[2][2]))
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.1f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(varstr+'\n')

    def PrintHStrain(hapData,SGData):
        pFile.write('\n Hydrostatic/elastic strain:\n')
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' names :'
        ptstr =  ' values:'
        varstr = ' refine:'
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.4g' % (hapData[0][i])
            varstr += '%12s' % (str(hapData[1][i]))
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(varstr+'\n')

    def PrintSHPO(hapData):
        pFile.write('\n Spherical harmonics preferred orientation: Order: %d Refine? %s\n'%(hapData[4],hapData[2]))
        ptlbls = ' names :'
        ptstr =  ' values:'
        for item in hapData[5]:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[5][item]) 
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
    
    def PrintBabinet(hapData):
        pFile.write('\n Babinet form factor modification:\n')
        ptlbls = ' names :'
        ptstr =  ' values:'
        varstr = ' refine:'
        for name in ['BabA','BabU']:
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (hapData[name][0])
            varstr += '%12s' % (str(hapData[name][1]))
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(varstr+'\n')
        
    hapDict = {}
    hapVary = []
    controlDict = {}
    
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        SGData = Phases[phase]['General']['SGData']
        cell = Phases[phase]['General']['Cell'][1:7]
        A = G2lat.cell2A(cell)
        if Phases[phase]['General'].get('Modulated',False):
            SSGData = Phases[phase]['General']['SSGData']
            Vec,x,maxH = Phases[phase]['General']['SuperVec']
        pId = Phases[phase]['pId']
        histoList = list(HistoPhase.keys())
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            if not HistoPhase[histogram]['Use']:        #remove previously created  & now unused phase reflection list
                if phase in Histograms[histogram]['Reflection Lists']:
                    del Histograms[histogram]['Reflection Lists'][phase]
                continue
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            if 'PWDR' in histogram:
                limits = Histogram['Limits'][1]
                inst = Histogram['Instrument Parameters'][0]    #TODO - grab table here if present
                if 'C' in inst['Type'][1]:
                    try:
                        wave = inst['Lam'][1]
                    except KeyError:
                        wave = inst['Lam1'][1]
                    dmin = wave/(2.0*sind(limits[1]/2.0))
                elif 'T' in inst['Type'][0]:
                    dmin = limits[0]/inst['difC'][1]
                pfx = str(pId)+':'+str(hId)+':'
                if Phases[phase]['General']['doPawley']:
                    hapDict[pfx+'LeBail'] = False           #Pawley supercedes LeBail
                    hapDict[pfx+'newLeBail'] = True
                    Tmin = G2lat.Dsp2pos(inst,dmin)
                    if 'C' in inst['Type'][1]:
                        limits[1] = min(limits[1],Tmin)
                    else:
                        limits[0] = max(limits[0],Tmin)
                else:
                    hapDict[pfx+'LeBail'] = hapData.get('LeBail',False)
                    hapDict[pfx+'newLeBail'] = hapData.get('newLeBail',True)
                if Phases[phase]['General']['Type'] == 'magnetic':
                    dmin = max(dmin,Phases[phase]['General']['MagDmin'])
                for item in ['Scale','Extinction']:
                    hapDict[pfx+item] = hapData[item][0]
                    if hapData[item][1] and not hapDict[pfx+'LeBail']:
                        hapVary.append(pfx+item)
                names = G2spc.HStrainNames(SGData)
                HSvals = []
                for i,name in enumerate(names):
                    hapDict[pfx+name] = hapData['HStrain'][0][i]
                    HSvals.append(hapDict[pfx+name])
                    if hapData['HStrain'][1][i]:
#                    if hapData['HStrain'][1][i] and not hapDict[pfx+'LeBail']:
                        hapVary.append(pfx+name)
                controlDict[pfx+'poType'] = hapData['Pref.Ori.'][0]
                if hapData['Pref.Ori.'][0] == 'MD':
                    hapDict[pfx+'MD'] = hapData['Pref.Ori.'][1]
                    controlDict[pfx+'MDAxis'] = hapData['Pref.Ori.'][3]
                    if hapData['Pref.Ori.'][2] and not hapDict[pfx+'LeBail']:
                        hapVary.append(pfx+'MD')
                else:                           #'SH' spherical harmonics
                    controlDict[pfx+'SHord'] = hapData['Pref.Ori.'][4]
                    controlDict[pfx+'SHncof'] = len(hapData['Pref.Ori.'][5])
                    controlDict[pfx+'SHnames'] = G2lat.GenSHCoeff(SGData['SGLaue'],'0',controlDict[pfx+'SHord'],False)
                    controlDict[pfx+'SHhkl'] = []
                    try: #patch for old Pref.Ori. items
                        controlDict[pfx+'SHtoler'] = 0.1
                        if hapData['Pref.Ori.'][6][0] != '':
                            controlDict[pfx+'SHhkl'] = [eval(a.replace(' ',',')) for a in hapData['Pref.Ori.'][6]]
                        controlDict[pfx+'SHtoler'] = hapData['Pref.Ori.'][7]
                    except IndexError:
                        pass
                    for item in hapData['Pref.Ori.'][5]:
                        hapDict[pfx+item] = hapData['Pref.Ori.'][5][item]
                        if hapData['Pref.Ori.'][2] and not hapDict[pfx+'LeBail']:
                            hapVary.append(pfx+item)
                for item in ['Mustrain','Size']:
                    controlDict[pfx+item+'Type'] = hapData[item][0]
                    hapDict[pfx+item+';mx'] = hapData[item][1][2]
                    if hapData[item][2][2]:
                        hapVary.append(pfx+item+';mx')
                    if hapData[item][0] in ['isotropic','uniaxial']:
                        hapDict[pfx+item+';i'] = hapData[item][1][0]
                        if hapData[item][2][0]:
                            hapVary.append(pfx+item+';i')
                        if hapData[item][0] == 'uniaxial':
                            controlDict[pfx+item+'Axis'] = hapData[item][3]
                            hapDict[pfx+item+';a'] = hapData[item][1][1]
                            if hapData[item][2][1]:
                                hapVary.append(pfx+item+';a')
                    else:       #generalized for mustrain or ellipsoidal for size
                        Nterms = len(hapData[item][4])
                        if item == 'Mustrain':
                            names = G2spc.MustrainNames(SGData)
                            pwrs = []
                            for name in names:
                                h,k,l = name[1:]
                                pwrs.append([int(h),int(k),int(l)])
                            controlDict[pfx+'MuPwrs'] = pwrs
                        for i in range(Nterms):
                            sfx = ';'+str(i)
                            hapDict[pfx+item+sfx] = hapData[item][4][i]
                            if hapData[item][5][i]:
                                hapVary.append(pfx+item+sfx)
                if Phases[phase]['General']['Type'] != 'magnetic':
                    for bab in ['BabA','BabU']:
                        hapDict[pfx+bab] = hapData['Babinet'][bab][0]
                        if hapData['Babinet'][bab][1] and not hapDict[pfx+'LeBail']:
                            hapVary.append(pfx+bab)
                                
                if Print: 
                    pFile.write('\n Phase: %s in histogram: %s\n'%(phase,histogram))
                    pFile.write(135*'='+'\n')
                    if hapDict[pfx+'LeBail']:
                        pFile.write(' Perform LeBail extraction\n')                     
                    else:
                        pFile.write(' Phase fraction  : %10.4f Refine? %s\n'%(hapData['Scale'][0],hapData['Scale'][1]))
                        pFile.write(' Extinction coeff: %10.4f Refine? %s\n'%(hapData['Extinction'][0],hapData['Extinction'][1]))
                        if hapData['Pref.Ori.'][0] == 'MD':
                            Ax = hapData['Pref.Ori.'][3]
                            pFile.write(' March-Dollase PO: %10.4f Refine? %s Axis: %d %d %d\n'%
                                (hapData['Pref.Ori.'][1],hapData['Pref.Ori.'][2],Ax[0],Ax[1],Ax[2]))
                        else: #'SH' for spherical harmonics
                            PrintSHPO(hapData['Pref.Ori.'])
                            pFile.write(' Penalty hkl list: %s tolerance: %.2f\n'%(controlDict[pfx+'SHhkl'],controlDict[pfx+'SHtoler']))
                    PrintSize(hapData['Size'])
                    PrintMuStrain(hapData['Mustrain'],SGData)
                    PrintHStrain(hapData['HStrain'],SGData)
                    if Phases[phase]['General']['Type'] != 'magnetic':
                        if hapData['Babinet']['BabA'][0]:
                            PrintBabinet(hapData['Babinet'])                        
                if resetRefList and (not hapDict[pfx+'LeBail'] or (hapData['LeBail'] and hapData['newLeBail'])):
                    if hapData.get('LeBail',True):         #stop regeneating reflections for LeBail
                        hapData['newLeBail'] = False
                    refList = []
                    Uniq = []
                    Phi = []
                    useExt = 'magnetic' in Phases[phase]['General']['Type'] and 'N' in inst['Type'][0]
                    if Phases[phase]['General'].get('Modulated',False):
                        ifSuper = True
                        HKLd = np.array(G2lat.GenSSHLaue(dmin,SGData,SSGData,Vec,maxH,A))
                        HKLd = G2mth.sortArray(HKLd,4,reverse=True)
                        for h,k,l,m,d in HKLd:
                            ext,mul,uniq,phi = G2spc.GenHKLf([h,k,l],SGData)
                            mul *= 2      # for powder overlap of Friedel pairs
                            if m or not ext or useExt:
                                if 'C' in inst['Type'][0]:
                                    pos = G2lat.Dsp2pos(inst,d)
                                    if limits[0] < pos < limits[1]:
                                        refList.append([h,k,l,m,mul,d, pos,0.0,0.0,0.0,100., 0.0,0.0,1.0,1.0,1.0])
                                        #... sig,gam,fotsq,fctsq, phase,icorr,prfo,abs,ext
                                        Uniq.append(uniq)
                                        Phi.append(phi)
                                elif 'T' in inst['Type'][0]:
                                    pos = G2lat.Dsp2pos(inst,d)
                                    if limits[0] < pos < limits[1]:
                                        wave = inst['difC'][1]*d/(252.816*inst['fltPath'][0])
                                        refList.append([h,k,l,m,mul,d, pos,0.0,0.0,0.0,100., 0.0,0.0,0.0,0.0,wave, 1.0,1.0,1.0])
                                        # ... sig,gam,fotsq,fctsq, phase,icorr,alp,bet,wave, prfo,abs,ext
                                        #TODO - if tabulated put alp & bet in here
                                        Uniq.append(uniq)
                                        Phi.append(phi)
                    else:
                        ifSuper = False
                        HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
                        HKLd = G2mth.sortArray(HKLd,3,reverse=True)
                        for h,k,l,d in HKLd:
                            ext,mul,uniq,phi = G2spc.GenHKLf([h,k,l],SGData)
                            mul *= 2      # for powder overlap of Friedel pairs
                            if ext and not useExt:
                                continue
                            if 'C' in inst['Type'][0]:
                                pos = G2lat.Dsp2pos(inst,d)
                                if limits[0] < pos < limits[1]:
                                    refList.append([h,k,l,mul,d, pos,0.0,0.0,0.0,100., 0.0,0.0,1.0,1.0,1.0])
                                    #... sig,gam,fotsq,fctsq, phase,icorr,prfo,abs,ext
                                    Uniq.append(uniq)
                                    Phi.append(phi)
                            elif 'T' in inst['Type'][0]:
                                pos = G2lat.Dsp2pos(inst,d)
                                if limits[0] < pos < limits[1]:
                                    wave = inst['difC'][1]*d/(252.816*inst['fltPath'][0])
                                    refList.append([h,k,l,mul,d, pos,0.0,0.0,0.0,100., 0.0,0.0,0.0,0.0,wave, 1.0,1.0,1.0])
                                    # ... sig,gam,fotsq,fctsq, phase,icorr,alp,bet,wave, prfo,abs,ext
                                    Uniq.append(uniq)
                                    Phi.append(phi)
                    Histogram['Reflection Lists'][phase] = {'RefList':np.array(refList),'FF':{},'Type':inst['Type'][0],'Super':ifSuper}
            elif 'HKLF' in histogram:
                inst = Histogram['Instrument Parameters'][0]
                hId = Histogram['hId']
                hfx = ':%d:'%(hId)
                for item in inst:
                    if type(inst) is not list and item != 'Type': continue # skip over non-refined items (such as InstName)
                    hapDict[hfx+item] = inst[item][1]
                pfx = str(pId)+':'+str(hId)+':'
                hapDict[pfx+'Scale'] = hapData['Scale'][0]
                if hapData['Scale'][1]:
                    hapVary.append(pfx+'Scale')
                                
                extApprox,extType,extParms = hapData['Extinction']
                controlDict[pfx+'EType'] = extType
                controlDict[pfx+'EApprox'] = extApprox
                if 'C' in inst['Type'][0]:
                    controlDict[pfx+'Tbar'] = extParms['Tbar']
                    controlDict[pfx+'Cos2TM'] = extParms['Cos2TM']
                if 'Primary' in extType:
                    Ekey = ['Ep',]
                elif 'I & II' in extType:
                    Ekey = ['Eg','Es']
                elif 'Secondary Type II' == extType:
                    Ekey = ['Es',]
                elif 'Secondary Type I' == extType:
                    Ekey = ['Eg',]
                else:   #'None'
                    Ekey = []
                for eKey in Ekey:
                    hapDict[pfx+eKey] = extParms[eKey][0]
                    if extParms[eKey][1]:
                        hapVary.append(pfx+eKey)
                for bab in ['BabA','BabU']:
                    hapDict[pfx+bab] = hapData['Babinet'][bab][0]
                    if hapData['Babinet'][bab][1]:
                        hapVary.append(pfx+bab)
                Twins = hapData.get('Twins',[[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,0]],])
                if len(Twins) == 1:
                    hapDict[pfx+'Flack'] = hapData.get('Flack',[0.,False])[0]
                    if hapData.get('Flack',[0,False])[1]:
                        hapVary.append(pfx+'Flack')
                sumTwFr = 0.
                controlDict[pfx+'TwinLaw'] = []
                controlDict[pfx+'TwinInv'] = []
                NTL = 0            
                for it,twin in enumerate(Twins):
                    if 'bool' in str(type(twin[0])):
                        controlDict[pfx+'TwinInv'].append(twin[0])
                        controlDict[pfx+'TwinLaw'].append(np.zeros((3,3)))
                    else:
                        NTL += 1
                        controlDict[pfx+'TwinInv'].append(False)
                        controlDict[pfx+'TwinLaw'].append(twin[0])
                    if it:
                        hapDict[pfx+'TwinFr:'+str(it)] = twin[1]
                        sumTwFr += twin[1]
                    else:
                        hapDict[pfx+'TwinFr:0'] = twin[1][0]
                        controlDict[pfx+'TwinNMN'] = twin[1][1]
                    if Twins[0][1][1]:
                        hapVary.append(pfx+'TwinFr:'+str(it))
                controlDict[pfx+'NTL'] = NTL
                #need to make constraint on TwinFr
                controlDict[pfx+'TwinLaw'] = np.array(controlDict[pfx+'TwinLaw'])
                if len(Twins) > 1:    #force sum to unity
                    hapDict[pfx+'TwinFr:0'] = 1.-sumTwFr
                if Print: 
                    pFile.write('\n Phase: %s in histogram: %s\n'%(phase,histogram))
                    pFile.write(135*'='+'\n')
                    pFile.write(' Scale factor     : %10.4f Refine? %s\n'%(hapData['Scale'][0],hapData['Scale'][1]))
                    if extType != 'None':
                        pFile.write(' Extinction  Type: %15s approx: %10s\n'%(extType,extApprox))
                        text = ' Parameters       :'
                        for eKey in Ekey:
                            text += ' %4s : %10.3e Refine? '%(eKey,extParms[eKey][0])+str(extParms[eKey][1])
                        pFile.write(text+'\n')
                    if hapData['Babinet']['BabA'][0]:
                        PrintBabinet(hapData['Babinet'])
                    if not SGData['SGInv'] and len(Twins) == 1:
                        pFile.write(' Flack parameter: %10.3f Refine? %s\n'%(hapData['Flack'][0],hapData['Flack'][1]))
                    if len(Twins) > 1:
                        for it,twin in enumerate(Twins):
                            if 'bool' in str(type(twin[0])):
                                pFile.write(' Nonmerohedral twin fr.: %5.3f Inv? %s Refine?\n'%
                                    (hapDict[pfx+'TwinFr:'+str(it)],str(controlDict[pfx+'TwinInv'][it])),Twins[0][1][1]) 
                            else:
                                pFile.write(' Twin law: %s Twin fr.: %5.3f Refine? %s\n'%
                                    (str(twin[0]).replace('\n',','),hapDict[pfx+'TwinFr:'+str(it)],Twins[0][1][1]))
                        
                Histogram['Reflection Lists'] = phase       
                
    return hapVary,hapDict,controlDict
    
def SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,FFtables,Print=True,pFile=None):
    'needs a doc string'
    
    def PrintSizeAndSig(hapData,sizeSig):
        line = '\n Size model:     %9s'%(hapData[0])
        refine = False
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.5f'%(hapData[1][0])
            if sizeSig[0][0]:
                line += ', sig:%8.4f'%(sizeSig[0][0])
                refine = True
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.4f'%(hapData[1][1])
                if sizeSig[0][1]:
                    refine = True
                    line += ', sig:%8.4f'%(sizeSig[0][1])
            line += ' LG mix coeff.:%12.4f'%(hapData[1][2])
            if sizeSig[0][2]:
                refine = True
                line += ', sig:%8.4f'%(sizeSig[0][2])
            if refine:
                pFile.write(line+'\n')
        else:
            line += ' LG mix coeff.:%12.4f'%(hapData[1][2])
            if sizeSig[0][2]:
                refine = True
                line += ', sig:%8.4f'%(sizeSig[0][2])
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                if sizeSig[1][i]:
                    refine = True
                    sigstr += '%12.6f' % (sizeSig[1][i])
                else:
                    sigstr += 12*' '
            if refine:
                pFile.write(line+'\n')
                pFile.write(ptlbls+'\n')
                pFile.write(ptstr+'\n')
                pFile.write(sigstr+'\n')
        
    def PrintMuStrainAndSig(hapData,mustrainSig,SGData):
        line = '\n Mustrain model: %9s\n'%(hapData[0])
        refine = False
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.1f'%(hapData[1][0])
            if mustrainSig[0][0]:
                line += ', sig:%8.1f'%(mustrainSig[0][0])
                refine = True
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.1f'%(hapData[1][1])
                if mustrainSig[0][1]:
                     line += ', sig:%8.1f'%(mustrainSig[0][1])
            line += ' LG mix coeff.:%12.4f'%(hapData[1][2])
            if mustrainSig[0][2]:
                refine = True
                line += ', sig:%8.3f'%(mustrainSig[0][2])
            if refine:
                pFile.write(line+'\n')
        else:
            line += ' LG mix coeff.:%12.4f'%(hapData[1][2])
            if mustrainSig[0][2]:
                refine = True
                line += ', sig:%8.3f'%(mustrainSig[0][2])
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.1f' % (hapData[4][i])
                if mustrainSig[1][i]:
                    refine = True
                    sigstr += '%12.1f' % (mustrainSig[1][i])
                else:
                    sigstr += 12*' '
            if refine:
                pFile.write(line+'\n')
                pFile.write(ptlbls+'\n')
                pFile.write(ptstr+'\n')
                pFile.write(sigstr+'\n')
            
    def PrintHStrainAndSig(hapData,strainSig,SGData):
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' name  :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        refine = False
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.4g' % (hapData[0][i])
            if name in strainSig:
                refine = True
                sigstr += '%12.4g' % (strainSig[name])
            else:
                sigstr += 12*' '
        if refine:
            pFile.write('\n Hydrostatic/elastic strain:\n')
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(sigstr+'\n')
        
    def PrintSHPOAndSig(pfx,hapData,POsig):
        pFile.write('\n Spherical harmonics preferred orientation: Order: %d\n'%hapData[4])
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        for item in hapData[5]:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[5][item])
            if pfx+item in POsig:
                sigstr += '%12.3f'%(POsig[pfx+item])
            else:
                sigstr += 12*' ' 
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(sigstr+'\n')
        
    def PrintExtAndSig(pfx,hapData,ScalExtSig):
        pFile.write('\n Single crystal extinction: Type: %s Approx: %s\n'%(hapData[0],hapData[1]))
        text = ''
        for item in hapData[2]:
            if pfx+item in ScalExtSig:
                text += '       %s: '%(item)
                text += '%12.2e'%(hapData[2][item][0])
                if pfx+item in ScalExtSig:
                    text += ' sig: %12.2e'%(ScalExtSig[pfx+item])
        pFile.write(text+'\n')    
        
    def PrintBabinetAndSig(pfx,hapData,BabSig):
        pFile.write('\n Babinet form factor modification:\n')
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        for item in hapData:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[item][0])
            if pfx+item in BabSig:
                sigstr += '%12.3f'%(BabSig[pfx+item])
            else:
                sigstr += 12*' ' 
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(sigstr+'\n')
        
    def PrintTwinsAndSig(pfx,twinData,TwinSig):
        pFile.write('\n Twin Law fractions :\n')
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        for it,item in enumerate(twinData):
            ptlbls += '     twin #%d'%(it)
            if it:
                ptstr += '%12.3f'%(item[1])
            else:
                ptstr += '%12.3f'%(item[1][0])
            if pfx+'TwinFr:'+str(it) in TwinSig:
                sigstr += '%12.3f'%(TwinSig[pfx+'TwinFr:'+str(it)])
            else:
                sigstr += 12*' ' 
        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(sigstr+'\n')
        
    
    PhFrExtPOSig = {}
    SizeMuStrSig = {}
    ScalExtSig = {}
    BabSig = {}
    TwinFrSig = {}
    wtFrSum = {}
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        General = Phases[phase]['General']
        SGData = General['SGData']
        pId = Phases[phase]['pId']
        histoList = list(HistoPhase.keys())
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            if not Phases[phase]['Histograms'][histogram]['Use']:
                #skip if phase absent from this histogram
                continue
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            pfx = str(pId)+':'+str(hId)+':'
            if hId not in wtFrSum:
                wtFrSum[hId] = 0.
            if 'PWDR' in histogram:
                for item in ['Scale','Extinction']:
                    hapData[item][0] = parmDict[pfx+item]
                    if pfx+item in sigDict and not parmDict[pfx+'LeBail']:
                        PhFrExtPOSig.update({pfx+item:sigDict[pfx+item],})
                wtFrSum[hId] += hapData['Scale'][0]*General['Mass']
                if hapData['Pref.Ori.'][0] == 'MD':
                    hapData['Pref.Ori.'][1] = parmDict[pfx+'MD']
                    if pfx+'MD' in sigDict and not parmDict[pfx+'LeBail']:
                        PhFrExtPOSig.update({pfx+'MD':sigDict[pfx+'MD'],})
                else:                           #'SH' spherical harmonics
                    for item in hapData['Pref.Ori.'][5]:
                        hapData['Pref.Ori.'][5][item] = parmDict[pfx+item]
                        if pfx+item in sigDict and not parmDict[pfx+'LeBail']:
                            PhFrExtPOSig.update({pfx+item:sigDict[pfx+item],})
                SizeMuStrSig.update({pfx+'Mustrain':[[0,0,0],[0 for i in range(len(hapData['Mustrain'][4]))]],
                    pfx+'Size':[[0,0,0],[0 for i in range(len(hapData['Size'][4]))]],
                    pfx+'HStrain':{}})                  
                for item in ['Mustrain','Size']:
                    hapData[item][1][2] = parmDict[pfx+item+';mx']
#                    hapData[item][1][2] = min(1.,max(0.,hapData[item][1][2]))
                    if pfx+item+';mx' in sigDict:
                        SizeMuStrSig[pfx+item][0][2] = sigDict[pfx+item+';mx']
                    if hapData[item][0] in ['isotropic','uniaxial']:                    
                        hapData[item][1][0] = parmDict[pfx+item+';i']
                        if item == 'Size':
                            hapData[item][1][0] = min(10.,max(0.001,hapData[item][1][0]))
                        if pfx+item+';i' in sigDict: 
                            SizeMuStrSig[pfx+item][0][0] = sigDict[pfx+item+';i']
                        if hapData[item][0] == 'uniaxial':
                            hapData[item][1][1] = parmDict[pfx+item+';a']
                            if item == 'Size':
                                hapData[item][1][1] = min(10.,max(0.001,hapData[item][1][1]))                        
                            if pfx+item+';a' in sigDict:
                                SizeMuStrSig[pfx+item][0][1] = sigDict[pfx+item+';a']
                    else:       #generalized for mustrain or ellipsoidal for size
                        Nterms = len(hapData[item][4])
                        for i in range(Nterms):
                            sfx = ';'+str(i)
                            hapData[item][4][i] = parmDict[pfx+item+sfx]
                            if pfx+item+sfx in sigDict:
                                SizeMuStrSig[pfx+item][1][i] = sigDict[pfx+item+sfx]
                names = G2spc.HStrainNames(SGData)
                for i,name in enumerate(names):
                    hapData['HStrain'][0][i] = parmDict[pfx+name]
                    if pfx+name in sigDict:
                        SizeMuStrSig[pfx+'HStrain'][name] = sigDict[pfx+name]
                if Phases[phase]['General']['Type'] != 'magnetic':
                    for name in ['BabA','BabU']:
                        hapData['Babinet'][name][0] = parmDict[pfx+name]
                        if pfx+name in sigDict and not parmDict[pfx+'LeBail']:
                            BabSig[pfx+name] = sigDict[pfx+name]                
                
            elif 'HKLF' in histogram:
                for item in ['Scale','Flack']:
                    if parmDict.get(pfx+item):
                        hapData[item][0] = parmDict[pfx+item]
                        if pfx+item in sigDict:
                            ScalExtSig[pfx+item] = sigDict[pfx+item]
                for item in ['Ep','Eg','Es']:
                    if parmDict.get(pfx+item):
                        hapData['Extinction'][2][item][0] = parmDict[pfx+item]
                        if pfx+item in sigDict:
                            ScalExtSig[pfx+item] = sigDict[pfx+item]
                for item in ['BabA','BabU']:
                    hapData['Babinet'][item][0] = parmDict[pfx+item]
                    if pfx+item in sigDict:
                        BabSig[pfx+item] = sigDict[pfx+item]
                if 'Twins' in hapData:
                    it = 1
                    sumTwFr = 0.
                    while True:
                        try:
                            hapData['Twins'][it][1] = parmDict[pfx+'TwinFr:'+str(it)]
                            if pfx+'TwinFr:'+str(it) in sigDict:
                                TwinFrSig[pfx+'TwinFr:'+str(it)] = sigDict[pfx+'TwinFr:'+str(it)]
                            if it:
                                sumTwFr += hapData['Twins'][it][1]
                            it += 1
                        except KeyError:
                            break
                    hapData['Twins'][0][1][0] = 1.-sumTwFr

    if Print:
        for phase in Phases:
            HistoPhase = Phases[phase]['Histograms']
            General = Phases[phase]['General']
            SGData = General['SGData']
            pId = Phases[phase]['pId']
            histoList = list(HistoPhase.keys())
            histoList.sort()
            for histogram in histoList:
                try:
                    Histogram = Histograms[histogram]
                except KeyError:                        
                    #skip if histogram not included e.g. in a sequential refinement
                    continue
                hapData = HistoPhase[histogram]
                hId = Histogram['hId']
                Histogram['Residuals'][str(pId)+'::Name'] = phase
                pfx = str(pId)+':'+str(hId)+':'
                hfx = ':%s:'%(hId)
                if pfx+'Nref' not in Histogram['Residuals']:    #skip not used phase in histogram
                    continue
                pFile.write('\n Phase: %s in histogram: %s\n'%(phase,histogram))
                pFile.write(135*'='+'\n')
                if 'PWDR' in histogram:
                    pFile.write(' Final refinement RF, RF^2 = %.2f%%, %.2f%% on %d reflections\n'%
                        (Histogram['Residuals'][pfx+'Rf'],Histogram['Residuals'][pfx+'Rf^2'],Histogram['Residuals'][pfx+'Nref']))
                    pFile.write(' Durbin-Watson statistic = %.3f\n'%(Histogram['Residuals']['Durbin-Watson']))
                    pFile.write(' Bragg intensity sum = %.3g\n'%(Histogram['Residuals'][pfx+'sumInt']))
                    
                    if parmDict[pfx+'LeBail']:
                        pFile.write(' Performed LeBail extraction for phase %s in histogram %s\n'%(phase,histogram))
                    else:
                        if pfx+'Scale' in PhFrExtPOSig:
                            wtFr = hapData['Scale'][0]*General['Mass']/wtFrSum[hId]
                            sigwtFr = PhFrExtPOSig[pfx+'Scale']*wtFr/hapData['Scale'][0]
                            pFile.write(' Phase fraction  : %10.5f, sig %10.5f Weight fraction  : %8.5f, sig %10.5f\n'%
                                (hapData['Scale'][0],PhFrExtPOSig[pfx+'Scale'],wtFr,sigwtFr))
                        if pfx+'Extinction' in PhFrExtPOSig:
                            pFile.write(' Extinction coeff: %10.4f, sig %10.4f\n'%(hapData['Extinction'][0],PhFrExtPOSig[pfx+'Extinction']))
                        if hapData['Pref.Ori.'][0] == 'MD':
                            if pfx+'MD' in PhFrExtPOSig:
                                pFile.write(' March-Dollase PO: %10.4f, sig %10.4f\n'%(hapData['Pref.Ori.'][1],PhFrExtPOSig[pfx+'MD']))
                        else:
                            PrintSHPOAndSig(pfx,hapData['Pref.Ori.'],PhFrExtPOSig)
                    PrintSizeAndSig(hapData['Size'],SizeMuStrSig[pfx+'Size'])
                    PrintMuStrainAndSig(hapData['Mustrain'],SizeMuStrSig[pfx+'Mustrain'],SGData)
                    PrintHStrainAndSig(hapData['HStrain'],SizeMuStrSig[pfx+'HStrain'],SGData)
                    if Phases[phase]['General']['Type'] != 'magnetic' and not parmDict[pfx+'LeBail']:
                        if len(BabSig):
                            PrintBabinetAndSig(pfx,hapData['Babinet'],BabSig)
                    
                elif 'HKLF' in histogram:
                    Inst = Histogram['Instrument Parameters'][0]
                    pFile.write(' Final refinement RF, RF^2 = %.2f%%, %.2f%% on %d reflections (%d user rejected, %d sp.gp.extinct)\n'%
                        (Histogram['Residuals'][pfx+'Rf'],Histogram['Residuals'][pfx+'Rf^2'],Histogram['Residuals'][pfx+'Nref'],
                        Histogram['Residuals'][pfx+'Nrej'],Histogram['Residuals'][pfx+'Next']))
                    if FFtables != None and 'N' not in Inst['Type'][0]:
                        PrintFprime(FFtables,hfx,pFile)
                    pFile.write(' HKLF histogram weight factor = %.3f\n'%(Histogram['wtFactor']))
                    if pfx+'Scale' in ScalExtSig:
                        pFile.write(' Scale factor : %10.4f, sig %10.4f\n'%(hapData['Scale'][0],ScalExtSig[pfx+'Scale']))
                    if hapData['Extinction'][0] != 'None':
                        PrintExtAndSig(pfx,hapData['Extinction'],ScalExtSig)
                    if len(BabSig):
                        PrintBabinetAndSig(pfx,hapData['Babinet'],BabSig)
                    if pfx+'Flack' in ScalExtSig:
                        pFile.write(' Flack parameter : %10.3f, sig %10.3f\n'%(hapData['Flack'][0],ScalExtSig[pfx+'Flack']))
                    if len(TwinFrSig):
                        PrintTwinsAndSig(pfx,hapData['Twins'],TwinFrSig)

################################################################################
##### Histogram data
################################################################################        
                    
def GetHistogramData(Histograms,Print=True,pFile=None):
    'needs a doc string'
    
    def GetBackgroundParms(hId,Background):
        Back = Background[0]
        DebyePeaks = Background[1]
        bakType,bakFlag = Back[:2]
        backVals = Back[3:]
        backNames = [':'+str(hId)+':Back;'+str(i) for i in range(len(backVals))]
        backDict = dict(zip(backNames,backVals))
        backVary = []
        if bakFlag:
            backVary = backNames
        backDict[':'+str(hId)+':nDebye'] = DebyePeaks['nDebye']
        backDict[':'+str(hId)+':nPeaks'] = DebyePeaks['nPeaks']
        debyeDict = {}
        debyeList = []
        for i in range(DebyePeaks['nDebye']):
            debyeNames = [':'+str(hId)+':DebyeA;'+str(i),':'+str(hId)+':DebyeR;'+str(i),':'+str(hId)+':DebyeU;'+str(i)]
            debyeDict.update(dict(zip(debyeNames,DebyePeaks['debyeTerms'][i][::2])))
            debyeList += zip(debyeNames,DebyePeaks['debyeTerms'][i][1::2])
        debyeVary = []
        for item in debyeList:
            if item[1]:
                debyeVary.append(item[0])
        backDict.update(debyeDict)
        backVary += debyeVary
        peakDict = {}
        peakList = []
        for i in range(DebyePeaks['nPeaks']):
            peakNames = [':'+str(hId)+':BkPkpos;'+str(i),':'+str(hId)+ \
                ':BkPkint;'+str(i),':'+str(hId)+':BkPksig;'+str(i),':'+str(hId)+':BkPkgam;'+str(i)]
            peakDict.update(dict(zip(peakNames,DebyePeaks['peaksList'][i][::2])))
            peakList += zip(peakNames,DebyePeaks['peaksList'][i][1::2])
        peakVary = []
        for item in peakList:
            if item[1]:
                peakVary.append(item[0])
        backDict.update(peakDict)
        backVary += peakVary
        return bakType,backDict,backVary            
        
    def GetInstParms(hId,Inst):
        #patch
        if 'Z' not in Inst:
            Inst['Z'] = [0.0,0.0,False]
        dataType = Inst['Type'][0]
        instDict = {}
        insVary = []
        pfx = ':'+str(hId)+':'
        insKeys = list(Inst.keys())
        insKeys.sort()
        for item in insKeys:
            insName = pfx+item
            instDict[insName] = Inst[item][1]
            if len(Inst[item]) > 2 and Inst[item][2]:
                insVary.append(insName)
        if 'C' in dataType:
            instDict[pfx+'SH/L'] = max(instDict[pfx+'SH/L'],0.0005)
        elif 'T' in dataType:   #trap zero alp, bet coeff.
            if not instDict[pfx+'alpha']:
                instDict[pfx+'alpha'] = 1.0
            if not instDict[pfx+'beta-0'] and not instDict[pfx+'beta-1']:
                instDict[pfx+'beta-1'] = 1.0
        return dataType,instDict,insVary
        
    def GetSampleParms(hId,Sample):
        sampVary = []
        hfx = ':'+str(hId)+':'        
        sampDict = {hfx+'Gonio. radius':Sample['Gonio. radius'],hfx+'Omega':Sample['Omega'],
            hfx+'Chi':Sample['Chi'],hfx+'Phi':Sample['Phi']}
        for key in ('Temperature','Pressure','FreePrm1','FreePrm2','FreePrm3'):
            if key in Sample:
                sampDict[hfx+key] = Sample[key]
        Type = Sample['Type']
        if 'Bragg' in Type:             #Bragg-Brentano
            for item in ['Scale','Shift','Transparency','SurfRoughA','SurfRoughB']:
                sampDict[hfx+item] = Sample[item][0]
                if Sample[item][1]:
                    sampVary.append(hfx+item)
        elif 'Debye' in Type:        #Debye-Scherrer
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                sampDict[hfx+item] = Sample[item][0]
                if Sample[item][1]:
                    sampVary.append(hfx+item)
        return Type,sampDict,sampVary
        
    def PrintBackground(Background):
        Back = Background[0]
        DebyePeaks = Background[1]
        pFile.write('\n Background function: %s Refine? %s\n'%(Back[0],Back[1]))
        line = ' Coefficients: '
        for i,back in enumerate(Back[3:]):
            line += '%10.3f'%(back)
            if i and not i%10:
                line += '\n'+15*' '
        pFile.write(line+'\n')
        if DebyePeaks['nDebye']:
            pFile.write('\n Debye diffuse scattering coefficients\n')
            parms = ['DebyeA','DebyeR','DebyeU']
            line = ' names :  '
            for parm in parms:
                line += '%8s refine?'%(parm)
            pFile.write(line+'\n')
            for j,term in enumerate(DebyePeaks['debyeTerms']):
                line = ' term'+'%2d'%(j)+':'
                for i in range(3):
                    line += '%10.3f %5s'%(term[2*i],bool(term[2*i+1]))                    
                pFile.write(line+'\n')
        if DebyePeaks['nPeaks']:
            pFile.write('\n Single peak coefficients\n')
            parms =    ['BkPkpos','BkPkint','BkPksig','BkPkgam']
            line = ' names :  '
            for parm in parms:
                line += '%8s refine?'%(parm)
            pFile.write(line+'\n')
            for j,term in enumerate(DebyePeaks['peaksList']):
                line = ' peak'+'%2d'%(j)+':'
                for i in range(4):
                    line += '%12.3f %5s'%(term[2*i],bool(term[2*i+1]))                    
                pFile.write(line+'\n')
        
    def PrintInstParms(Inst):
        pFile.write('\n Instrument Parameters:\n')
        insKeys = list(Inst.keys())
        insKeys.sort()
        iBeg = 0
        Ok = True
        while Ok:
            ptlbls = ' name  :'
            ptstr =  ' value :'
            varstr = ' refine:'
            iFin = min(iBeg+9,len(insKeys))
            for item in insKeys[iBeg:iFin]:
                if item not in ['Type','Source']:
                    ptlbls += '%12s' % (item)
                    ptstr += '%12.6f' % (Inst[item][1])
                    if item in ['Lam1','Lam2','Azimuth','fltPath','2-theta',]:
                        varstr += 12*' '
                    else:
                        varstr += '%12s' % (str(bool(Inst[item][2])))
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(varstr+'\n')
            iBeg = iFin
            if iBeg == len(insKeys):
                Ok = False
            else:
                pFile.write('\n')
        
    def PrintSampleParms(Sample):
        pFile.write('\n Sample Parameters:\n')
        pFile.write(' Goniometer omega = %.2f, chi = %.2f, phi = %.2f\n'%
            (Sample['Omega'],Sample['Chi'],Sample['Phi']))
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        if 'Bragg' in Sample['Type']:
            for item in ['Scale','Shift','Transparency','SurfRoughA','SurfRoughB']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))
            
        elif 'Debye' in Type:        #Debye-Scherrer
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))

        pFile.write(ptlbls+'\n')
        pFile.write(ptstr+'\n')
        pFile.write(varstr+'\n')
        
    histDict = {}
    histVary = []
    controlDict = {}
    histoList = list(Histograms.keys())
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            controlDict[pfx+'wtFactor'] = Histogram['wtFactor']
            controlDict[pfx+'Limits'] = Histogram['Limits'][1]
            controlDict[pfx+'Exclude'] = Histogram['Limits'][2:]
            for excl in controlDict[pfx+'Exclude']:
                Histogram['Data'][0] = ma.masked_inside(Histogram['Data'][0],excl[0],excl[1])
            if controlDict[pfx+'Exclude']:
                ma.mask_rows(Histogram['Data'])
            Background = Histogram['Background']
            Type,bakDict,bakVary = GetBackgroundParms(hId,Background)
            controlDict[pfx+'bakType'] = Type
            histDict.update(bakDict)
            histVary += bakVary
            
            Inst = Histogram['Instrument Parameters']        #TODO ? ignores tabulated alp,bet & delt for TOF
            if 'T' in Type and len(Inst[1]):    #patch -  back-to-back exponential contribution to TOF line shape is removed
                print ('Warning: tabulated profile coefficients are ignored')
            Type,instDict,insVary = GetInstParms(hId,Inst[0])
            controlDict[pfx+'histType'] = Type
            if 'XC' in Type:
                if pfx+'Lam1' in instDict:
                    controlDict[pfx+'keV'] = 12.397639/instDict[pfx+'Lam1']
                else:
                    controlDict[pfx+'keV'] = 12.397639/instDict[pfx+'Lam']            
            histDict.update(instDict)
            histVary += insVary
            
            Sample = Histogram['Sample Parameters']
            Type,sampDict,sampVary = GetSampleParms(hId,Sample)
            controlDict[pfx+'instType'] = Type
            histDict.update(sampDict)
            histVary += sampVary
            
    
            if Print: 
                pFile.write('\n Histogram: %s histogram Id: %d\n'%(histogram,hId))
                pFile.write(135*'='+'\n')
                Units = {'C':' deg','T':' msec'}
                units = Units[controlDict[pfx+'histType'][2]]
                Limits = controlDict[pfx+'Limits']
                pFile.write(' Instrument type: %s\n'%Sample['Type'])
                pFile.write(' Histogram limits: %8.2f%s to %8.2f%s\n'%(Limits[0],units,Limits[1],units))
                if len(controlDict[pfx+'Exclude']):
                    excls = controlDict[pfx+'Exclude']
                    for excl in excls:
                        pFile.write(' Excluded region:  %8.2f%s to %8.2f%s\n'%(excl[0],units,excl[1],units)) 
                PrintSampleParms(Sample)
                PrintInstParms(Inst[0])
                PrintBackground(Background)
        elif 'HKLF' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            controlDict[pfx+'wtFactor'] = Histogram['wtFactor']
            Inst = Histogram['Instrument Parameters'][0]
            controlDict[pfx+'histType'] = Inst['Type'][0]
            if 'X' in Inst['Type'][0]:
                histDict[pfx+'Lam'] = Inst['Lam'][1]
                controlDict[pfx+'keV'] = 12.397639/histDict[pfx+'Lam']                    
    return histVary,histDict,controlDict
    
def SetHistogramData(parmDict,sigDict,Histograms,FFtables,Print=True,pFile=None):
    'needs a doc string'
    
    def SetBackgroundParms(pfx,Background,parmDict,sigDict):
        Back = Background[0]
        DebyePeaks = Background[1]
        lenBack = len(Back[3:])
        backSig = [0 for i in range(lenBack+3*DebyePeaks['nDebye']+4*DebyePeaks['nPeaks'])]
        for i in range(lenBack):
            Back[3+i] = parmDict[pfx+'Back;'+str(i)]
            if pfx+'Back;'+str(i) in sigDict:
                backSig[i] = sigDict[pfx+'Back;'+str(i)]
        if DebyePeaks['nDebye']:
            for i in range(DebyePeaks['nDebye']):
                names = [pfx+'DebyeA;'+str(i),pfx+'DebyeR;'+str(i),pfx+'DebyeU;'+str(i)]
                for j,name in enumerate(names):
                    DebyePeaks['debyeTerms'][i][2*j] = parmDict[name]
                    if name in sigDict:
                        backSig[lenBack+3*i+j] = sigDict[name]            
        if DebyePeaks['nPeaks']:
            for i in range(DebyePeaks['nPeaks']):
                names = [pfx+'BkPkpos;'+str(i),pfx+'BkPkint;'+str(i),
                    pfx+'BkPksig;'+str(i),pfx+'BkPkgam;'+str(i)]
                for j,name in enumerate(names):
                    DebyePeaks['peaksList'][i][2*j] = parmDict[name]
                    if name in sigDict:
                        backSig[lenBack+3*DebyePeaks['nDebye']+4*i+j] = sigDict[name]
        return backSig
        
    def SetInstParms(pfx,Inst,parmDict,sigDict):
        instSig = {}
        insKeys = list(Inst.keys())
        insKeys.sort()
        for item in insKeys:
            insName = pfx+item
            Inst[item][1] = parmDict[insName]
            if insName in sigDict:
                instSig[item] = sigDict[insName]
            else:
                instSig[item] = 0
        return instSig
        
    def SetSampleParms(pfx,Sample,parmDict,sigDict):
        if 'Bragg' in Sample['Type']:             #Bragg-Brentano
            sampSig = [0 for i in range(5)]
            for i,item in enumerate(['Scale','Shift','Transparency','SurfRoughA','SurfRoughB']):
                Sample[item][0] = parmDict[pfx+item]
                if pfx+item in sigDict:
                    sampSig[i] = sigDict[pfx+item]
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            sampSig = [0 for i in range(4)]
            for i,item in enumerate(['Scale','Absorption','DisplaceX','DisplaceY']):
                Sample[item][0] = parmDict[pfx+item]
                if pfx+item in sigDict:
                    sampSig[i] = sigDict[pfx+item]
        return sampSig
        
    def PrintBackgroundSig(Background,backSig):
        Back = Background[0]
        DebyePeaks = Background[1]
        valstr = ' value : '
        sigstr = ' sig   : '
        refine = False
        for i,back in enumerate(Back[3:]):
            valstr += '%10.4g'%(back)
            if Back[1]:
                refine = True
                sigstr += '%10.4g'%(backSig[i])
            else:
                sigstr += 10*' '
        if refine:
            pFile.write('\n Background function: %s\n'%Back[0])
            pFile.write(valstr+'\n')
            pFile.write(sigstr+'\n')
        if DebyePeaks['nDebye']:
            ifAny = False
            ptfmt = "%12.3f"
            names =  ' names :'
            ptstr =  ' values:'
            sigstr = ' esds  :'
            for item in sigDict:
                if 'Debye' in item:
                    ifAny = True
                    names += '%12s'%(item)
                    ptstr += ptfmt%(parmDict[item])
                    sigstr += ptfmt%(sigDict[item])
            if ifAny:
                pFile.write('\n Debye diffuse scattering coefficients\n')
                pFile.write(names+'\n')
                pFile.write(ptstr+'\n')
                pFile.write(sigstr+'\n')
        if DebyePeaks['nPeaks']:
            pFile.write('\n Single peak coefficients:\n')
            parms =    ['BkPkpos','BkPkint','BkPksig','BkPkgam']
            line = ' peak no. '
            for parm in parms:
                line += '%14s%12s'%(parm.center(14),'esd'.center(12))
            pFile.write(line+'\n')
            for ip in range(DebyePeaks['nPeaks']):
                ptstr = ' %4d '%(ip)
                for parm in parms:
                    fmt = '%14.3f'
                    efmt = '%12.3f'
                    if parm == 'BkPkpos':
                        fmt = '%14.4f'
                        efmt = '%12.4f'
                    name = pfx+parm+';%d'%(ip)
                    ptstr += fmt%(parmDict[name])
                    if name in sigDict:
                        ptstr += efmt%(sigDict[name])
                    else:
                        ptstr += 12*' '
                pFile.write(ptstr+'\n')
        sumBk = np.array(Histogram['sumBk'])
        pFile.write(' Background sums: empirical %.3g, Debye %.3g, peaks %.3g, Total %.3g\n'%
            (sumBk[0],sumBk[1],sumBk[2],np.sum(sumBk)))
        
    def PrintInstParmsSig(Inst,instSig):
        refine = False
        insKeys = list(instSig.keys())
        insKeys.sort()
        iBeg = 0
        Ok = True
        while Ok:
            ptlbls = ' names :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            iFin = min(iBeg+9,len(insKeys))
            for name in insKeys[iBeg:iFin]:
                if name not in  ['Type','Lam1','Lam2','Azimuth','Source','fltPath']:
                    ptlbls += '%12s' % (name)
                    ptstr += '%12.6f' % (Inst[name][1])
                    if instSig[name]:
                        refine = True
                        sigstr += '%12.6f' % (instSig[name])
                    else:
                        sigstr += 12*' '
            if refine:
                pFile.write('\n Instrument Parameters:\n')
                pFile.write(ptlbls+'\n')
                pFile.write(ptstr+'\n')
                pFile.write(sigstr+'\n')
            iBeg = iFin
            if iBeg == len(insKeys):
                Ok = False
        
    def PrintSampleParmsSig(Sample,sampleSig):
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        refine = False
        if 'Bragg' in Sample['Type']:
            for i,item in enumerate(['Scale','Shift','Transparency','SurfRoughA','SurfRoughB']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                if sampleSig[i]:
                    refine = True
                    sigstr += '%14.4f'%(sampleSig[i])
                else:
                    sigstr += 14*' '
            
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            for i,item in enumerate(['Scale','Absorption','DisplaceX','DisplaceY']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                if sampleSig[i]:
                    refine = True
                    sigstr += '%14.4f'%(sampleSig[i])
                else:
                    sigstr += 14*' '

        if refine:
            pFile.write('\n Sample Parameters:\n')
            pFile.write(ptlbls+'\n')
            pFile.write(ptstr+'\n')
            pFile.write(sigstr+'\n')
        
    histoList = list(Histograms.keys())
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            Background = Histogram['Background']
            backSig = SetBackgroundParms(pfx,Background,parmDict,sigDict)
            
            Inst = Histogram['Instrument Parameters'][0]
            instSig = SetInstParms(pfx,Inst,parmDict,sigDict)
        
            Sample = Histogram['Sample Parameters']
            sampSig = SetSampleParms(pfx,Sample,parmDict,sigDict)

            pFile.write('\n Histogram: %s histogram Id: %d\n'%(histogram,hId))
            pFile.write(135*'='+'\n')
            pFile.write(' PWDR histogram weight factor = '+'%.3f\n'%(Histogram['wtFactor']))
            pFile.write(' Final refinement wR = %.2f%% on %d observations in this histogram\n'%
                (Histogram['Residuals']['wR'],Histogram['Residuals']['Nobs']))
            pFile.write(' Other residuals: R = %.2f%%, R-bkg = %.2f%%, wR-bkg = %.2f%% wRmin = %.2f%%\n'%
                (Histogram['Residuals']['R'],Histogram['Residuals']['Rb'],Histogram['Residuals']['wR'],Histogram['Residuals']['wRmin']))
            if Print:
                pFile.write(' Instrument type: %s\n'%Sample['Type'])
                if FFtables != None and 'N' not in Inst['Type'][0]:
                    PrintFprime(FFtables,pfx,pFile)
                PrintSampleParmsSig(Sample,sampSig)
                PrintInstParmsSig(Inst,instSig)
                PrintBackgroundSig(Background,backSig)
                
