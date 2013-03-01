# -*- coding: utf-8 -*-
#GSASIIstructure - structure computation routines
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
DEBUG = False

def GetControls(GPXfile):
    ''' Returns dictionary of control items found in GSASII gpx file
    input:
        GPXfile = .gpx full file name
    return:
        Controls = dictionary of control items
    '''
    Controls = {'deriv type':'analytic Hessian','max cyc':3,'max Hprocess':1,
        'max Rprocess':1,'min dM/M':0.0001,'shift factor':1.}
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Controls':
            Controls.update(datum[1])
    fl.close()
    return Controls
    
def GetConstraints(GPXfile):
    '''Read the constraints from the GPX file and interpret them
    '''
    constList = []
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Constraints':
            constDict = datum[1]
            for item in constDict:
                constList += constDict[item]
    fl.close()
    constDict,fixedList,ignored = ProcessConstraints(constList)
    if ignored:
        print ignored,'old-style Constraints were rejected'
    return constDict,fixedList
    
def ProcessConstraints(constList):
    "interpret constraints"
    constDict = []
    fixedList = []
    ignored = 0
    for item in constList:
        if item[-1] == 'h':
            # process a hold
            fixedList.append('0')
            constDict.append({item[0][1]:0.0})
        elif item[-1] == 'f':
            # process a new variable
            fixedList.append(None)
            constDict.append({})
            for term in item[:-3]:
                constDict[-1][term[1]] = term[0]
            #constFlag[-1] = ['Vary']
        elif item[-1] == 'c': 
            # process a contraint relationship
            fixedList.append(str(item[-3]))
            constDict.append({})
            for term in item[:-3]:
                constDict[-1][term[1]] = term[0]
            #constFlag[-1] = ['VaryFree']
        elif item[-1] == 'e':
            # process an equivalence
            firstmult = None
            eqlist = []
            for term in item[:-3]:
                if term[0] == 0: term[0] = 1.0
                if firstmult is None:
                    firstmult,firstvar = term
                else:
                    eqlist.append([term[1],firstmult/term[0]])
            G2mv.StoreEquivalence(firstvar,eqlist)
        else:
            ignored += 1
    return constDict,fixedList,ignored

def CheckConstraints(GPXfile):
    '''Load constraints and related info and return any error or warning messages'''
    # init constraints
    G2mv.InitVars()    
    # get variables
    Histograms,Phases = GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        return 'Error: No Phases!',''
    if not Histograms:
        return 'Error: no diffraction data',''
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = GetPhaseData(Phases,RestraintDict=None,Print=False)
    hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histograms,Print=False)
    histVary,histDict,controlDict = GetHistogramData(Histograms,Print=False)
    varyList = phaseVary+hapVary+histVary
    constrDict,fixedList = GetConstraints(GPXfile)
    return G2mv.CheckConstraints(varyList,constrDict,fixedList)
    
def GetRestraints(GPXfile):
    '''Read the restraints from the GPX file
    '''
    fl = open(GPXfile,'rb')
    while True:
        try:
            data = cPickle.load(fl)
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
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Rigid bodies':
            rigidbodyDict = datum[1]
    fl.close()
    return rigidbodyDict
        
def GetPhaseNames(GPXfile):
    ''' Returns a list of phase names found under 'Phases' in GSASII gpx file
    input: 
        GPXfile = gpx full file name
    return: 
        PhaseNames = list of phase names
    '''
    fl = open(GPXfile,'rb')
    PhaseNames = []
    while True:
        try:
            data = cPickle.load(fl)
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
    input:
        GPXfile = gpx full file name
        PhaseName = phase name
    return:
        phase dictionary
    '''        
    fl = open(GPXfile,'rb')
    General = {}
    Atoms = []
    while True:
        try:
            data = cPickle.load(fl)
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
    input: 
        GPXfile = .gpx full file name
        hNames = list of histogram names 
    return: 
        Histograms = dictionary of histograms (types = PWDR & HKLF)
    """
    fl = open(GPXfile,'rb')
    Histograms = {}
    while True:
        try:
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        hist = datum[0]
        if hist in hNames:
            if 'PWDR' in hist[:4]:
                PWDRdata = {}
                try:
                    PWDRdata.update(datum[1][0])        #weight factor
                except ValueError:
                    PWDRdata['wtFactor'] = 1.0          #patch
                PWDRdata['Data'] = datum[1][1]          #powder data arrays
                PWDRdata[data[2][0]] = data[2][1]       #Limits
                PWDRdata[data[3][0]] = data[3][1]       #Background
                PWDRdata[data[4][0]] = data[4][1]       #Instrument parameters
                PWDRdata[data[5][0]] = data[5][1]       #Sample parameters
                try:
                    PWDRdata[data[9][0]] = data[9][1]       #Reflection lists might be missing
                except IndexError:
                    PWDRdata['Reflection Lists'] = {}
    
                Histograms[hist] = PWDRdata
            elif 'HKLF' in hist[:4]:
                HKLFdata = {}
                try:
                    HKLFdata.update(datum[1][0])        #weight factor
                except ValueError:
                    HKLFdata['wtFactor'] = 1.0          #patch
                HKLFdata['Data'] = datum[1][1]
                HKLFdata[data[1][0]] = data[1][1]       #Instrument parameters
                HKLFdata['Reflection Lists'] = None
                Histograms[hist] = HKLFdata           
    fl.close()
    return Histograms
    
def GetHistogramNames(GPXfile,hType):
    """ Returns a list of histogram names found in GSASII gpx file
    input: 
        GPXfile = .gpx full file name
        hType = list ['PWDR','HKLF'] 
    return: 
        HistogramNames = list of histogram names (types = PWDR & HKLF)
    """
    fl = open(GPXfile,'rb')
    HistogramNames = []
    while True:
        try:
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        if datum[0][:4] in hType:
            HistogramNames.append(datum[0])
    fl.close()
    return HistogramNames
    
def GetUsedHistogramsAndPhases(GPXfile):
    ''' Returns all histograms that are found in any phase
    and any phase that uses a histogram
    input:
        GPXfile = .gpx full file name
    return:
        Histograms = dictionary of histograms as {name:data,...}
        Phases = dictionary of phases that use histograms
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
        if Phase['Histograms']:
            if phase not in Phases:
                pId = phaseNames.index(phase)
                Phase['pId'] = pId
                Phases[phase] = Phase
            for hist in Phase['Histograms']:
                if 'Use' not in Phase['Histograms'][hist]:      #patch
                    Phase['Histograms'][hist]['Use'] = True         
                if hist not in Histograms and Phase['Histograms'][hist]['Use']:
                    Histograms[hist] = allHistograms[hist]
                    hId = histoList.index(hist)
                    Histograms[hist]['hId'] = hId
    return Histograms,Phases
    
def getBackupName(GPXfile,makeBack):
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
    import distutils.file_util as dfu
    GPXback = getBackupName(GPXfile,makeBack)
    dfu.copy_file(GPXfile,GPXback)
    return GPXback

def SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,CovData,makeBack=True):
    ''' Updates gpxfile from all histograms that are found in any phase
    and any phase that used a histogram
    input:
        GPXfile = .gpx full file name
        Histograms = dictionary of histograms as {name:data,...}
        Phases = dictionary of phases that use histograms
        CovData = dictionary of refined variables, varyList, & covariance matrix
        makeBack = True if new backup of .gpx file is to be made; else use the last one made
    '''
                        
    GPXback = GPXBackup(GPXfile,makeBack)
    print 'Read from file:',GPXback
    print 'Save to file  :',GPXfile
    infile = open(GPXback,'rb')
    outfile = open(GPXfile,'wb')
    while True:
        try:
            data = cPickle.load(infile)
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
        try:
            histogram = Histograms[datum[0]]
#            print 'found ',datum[0]
            data[0][1][1] = histogram['Data']
            for datus in data[1:]:
#                print '    read: ',datus[0]
                if datus[0] in ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']:
                    datus[1] = histogram[datus[0]]
        except KeyError:
            pass
                                
        cPickle.dump(data,outfile,1)
    infile.close()
    outfile.close()
    print 'GPX file save successful'
    
def SetSeqResult(GPXfile,Histograms,SeqResult):
    GPXback = GPXBackup(GPXfile)
    print 'Read from file:',GPXback
    print 'Save to file  :',GPXfile
    infile = open(GPXback,'rb')
    outfile = open(GPXfile,'wb')
    while True:
        try:
            data = cPickle.load(infile)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Sequental results':
            data[0][1] = SeqResult
        try:
            histogram = Histograms[datum[0]]
            data[0][1][1] = histogram['Data']
            for datus in data[1:]:
                if datus[0] in ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']:
                    datus[1] = histogram[datus[0]]
        except KeyError:
            pass
                                
        cPickle.dump(data,outfile,1)
    infile.close()
    outfile.close()
    print 'GPX file save successful'
                        
def ShowBanner(pFile=None):
    print >>pFile,80*'*'
    print >>pFile,'   General Structure Analysis System-II Crystal Structure Refinement'
    print >>pFile,'              by Robert B. Von Dreele & Brian H. Toby'
    print >>pFile,'                Argonne National Laboratory(C), 2010'
    print >>pFile,' This product includes software developed by the UChicago Argonne, LLC,' 
    print >>pFile,'            as Operator of Argonne National Laboratory.'
    print >>pFile,80*'*','\n'

def ShowControls(Controls,pFile=None):
    print >>pFile,' Least squares controls:'
    print >>pFile,' Refinement type: ',Controls['deriv type']
    if 'Hessian' in Controls['deriv type']:
        print >>pFile,' Maximum number of cycles:',Controls['max cyc']
    else:
        print >>pFile,' Minimum delta-M/M for convergence: ','%.2g'%(Controls['min dM/M'])
    print >>pFile,' Initial shift factor: ','%.3f'%(Controls['shift factor'])
    
def GetFFtable(General):
    ''' returns a dictionary of form factor data for atom types found in General
    input:
        General = dictionary of phase info.; includes AtomTypes
    return:
        FFtable = dictionary of form factor data; key is atom type
    '''
    atomTypes = General['AtomTypes']
    FFtable = {}
    for El in atomTypes:
        FFs = G2el.GetFormFactorCoeff(El.split('+')[0].split('-')[0])
        for item in FFs:
            if item['Symbol'] == El.upper():
                FFtable[El] = item
    return FFtable
    
def GetBLtable(General):
    ''' returns a dictionary of neutron scattering length data for atom types & isotopes found in General
    input:
        General = dictionary of phase info.; includes AtomTypes & Isotopes
    return:
        BLtable = dictionary of scattering length data; key is atom type
    '''
    atomTypes = General['AtomTypes']
    BLtable = {}
    isotopes = General['Isotopes']
    isotope = General['Isotope']
    for El in atomTypes:
        BLtable[El] = [isotope[El],isotopes[El][isotope[El]]]
    return BLtable
        
def GetPawleyConstr(SGLaue,PawleyRef,pawleyVary):
#    if SGLaue in ['-1','2/m','mmm']:
#        return                      #no Pawley symmetry required constraints
    eqvDict = {}
    for i,varyI in enumerate(pawleyVary):
        eqvDict[varyI] = []
        refI = int(varyI.split(':')[-1])
        ih,ik,il = PawleyRef[refI][:3]
        dspI = PawleyRef[refI][4]
        for varyJ in pawleyVary[i+1:]:
            refJ = int(varyJ.split(':')[-1])
            jh,jk,jl = PawleyRef[refJ][:3]
            dspJ = PawleyRef[refJ][4]
            if SGLaue in ['4/m','4/mmm']:
                isum = ih**2+ik**2
                jsum = jh**2+jk**2
                if abs(il) == abs(jl) and isum == jsum:
                    eqvDict[varyI].append(varyJ) 
            elif SGLaue in ['3R','3mR']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2*jl**2
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
    if SGData['SGLaue'] in ['-1',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3']
        elif SGData['SGUniq'] == 'b':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A4']
        else:
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A5']
    elif SGData['SGLaue'] in ['mmm',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2']
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        return [pfx+'A0',pfx+'A2']
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        return [pfx+'A0',pfx+'A2']
    elif SGData['SGLaue'] in ['3R', '3mR']:
        return [pfx+'A0',pfx+'A3']                       
    elif SGData['SGLaue'] in ['m3m','m3']:
        return [pfx+'A0',]
        
################################################################################
##### Phase data
################################################################################        
                    
def GetPhaseData(PhaseData,RestraintDict={},Print=True,pFile=None):
            
    def PrintFFtable(FFtable):
        print >>pFile,'\n X-ray scattering factors:'
        print >>pFile,'   Symbol     fa                                      fb                                      fc'
        print >>pFile,99*'-'
        for Ename in FFtable:
            ffdata = FFtable[Ename]
            fa = ffdata['fa']
            fb = ffdata['fb']
            print >>pFile,' %8s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f' %  \
                (Ename.ljust(8),fa[0],fa[1],fa[2],fa[3],fb[0],fb[1],fb[2],fb[3],ffdata['fc'])
                
    def PrintBLtable(BLtable):
        print >>pFile,'\n Neutron scattering factors:'
        print >>pFile,'   Symbol   isotope       mass       b       resonant terms'
        print >>pFile,99*'-'
        for Ename in BLtable:
            bldata = BLtable[Ename]
            isotope = bldata[0]
            mass = bldata[1][0]
            blen = bldata[1][1]
            bres = []
            if len(bldata[1]) > 2:
                bres = bldata[1][2:]
            line = ' %8s%11s %10.3f %8.3f'%(Ename.ljust(8),isotope.center(11),mass,blen)
            for item in bres:
                line += '%10.5g'%(item)
            print >>pFile,line
            
    #def PrintRBObjects()
                
    def PrintAtoms(General,Atoms):
        cx,ct,cs,cia = General['AtomPtrs']
        print >>pFile,'\n Atoms:'
        line = '   name    type  refine?   x         y         z    '+ \
            '  frac site sym  mult I/A   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'magnetic':
            line += '   Mx     My     Mz'
        elif General['Type'] == 'macromolecular':
            line = ' res no residue chain'+line
        print >>pFile,line
        if General['Type'] == 'nuclear':
            print >>pFile,135*'-'
            for i,at in enumerate(Atoms):
                line = '%7s'%(at[ct-1])+'%7s'%(at[ct])+'%7s'%(at[ct+1])+'%10.5f'%(at[cx])+'%10.5f'%(at[cx+1])+ \
                    '%10.5f'%(at[cx+2])+'%8.3f'%(at[cx+3])+'%7s'%(at[cs])+'%5d'%(at[cs+1])+'%5s'%(at[cia])
                if at[cia] == 'I':
                    line += '%8.4f'%(at[cia+1])+48*' '
                else:
                    line += 8*' '
                    for j in range(6):
                        line += '%8.4f'%(at[cia+1+j])
                print >>pFile,line
        elif General['Type'] == 'macromolecular':
            print >>pFile,135*'-'            
            for i,at in enumerate(Atoms):
                line = '%7s'%(at[0])+'%7s'%(at[1])+'%7s'%(at[2])+'%7s'%(at[ct-1])+'%7s'%(at[ct])+'%7s'%(at[ct+1])+'%10.5f'%(at[cx])+'%10.5f'%(at[cx+1])+ \
                    '%10.5f'%(at[cx+2])+'%8.3f'%(at[cx+3])+'%7s'%(at[cs])+'%5d'%(at[cs+1])+'%5s'%(at[cia])
                if at[cia] == 'I':
                    line += '%8.4f'%(at[cia+1])+48*' '
                else:
                    line += 8*' '
                    for j in range(6):
                        line += '%8.4f'%(at[cia+1+j])
                print >>pFile,line
        
    def PrintTexture(textureData):
        topstr = '\n Spherical harmonics texture: Order:' + \
            str(textureData['Order'])
        if textureData['Order']:
            print >>pFile,topstr+' Refine? '+str(textureData['SH Coeff'][0])
        else:
            print >>pFile,topstr
            return
        names = ['omega','chi','phi']
        line = '\n'
        for name in names:
            line += ' SH '+name+':'+'%12.4f'%(textureData['Sample '+name][1])+' Refine? '+str(textureData['Sample '+name][0])
        print >>pFile,line
        print >>pFile,'\n Texture coefficients:'
        ptlbls = ' names :'
        ptstr =  ' values:'
        SHcoeff = textureData['SH Coeff'][1]
        for item in SHcoeff:
            ptlbls += '%12s'%(item)
            ptstr += '%12.4f'%(SHcoeff[item]) 
        print >>pFile,ptlbls
        print >>pFile,ptstr
        
    if Print:print  >>pFile,' Phases:'
    phaseVary = []
    phaseDict = {}
    phaseConstr = {}
    pawleyLookup = {}
    FFtables = {}                   #scattering factors - xrays
    BLtables = {}                   # neutrons
    Natoms = {}
    AtMults = {}
    AtIA = {}
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    atomIndx = {}
    for name in PhaseData:
        General = PhaseData[name]['General']
        pId = PhaseData[name]['pId']
        pfx = str(pId)+'::'
        FFtable = GetFFtable(General)
        BLtable = GetBLtable(General)
        FFtables.update(FFtable)
        BLtables.update(BLtable)
        Atoms = PhaseData[name]['Atoms']
        AtLookup = G2mth.FillAtomLookUp(Atoms)
        try:
            PawleyRef = PhaseData[name]['Pawley ref']
        except KeyError:
            PawleyRef = []
        SGData = General['SGData']
        SGtext = G2spc.SGPrint(SGData)
        cell = General['Cell']
        A = G2lat.cell2A(cell[1:7])
        phaseDict.update({pfx+'A0':A[0],pfx+'A1':A[1],pfx+'A2':A[2],
            pfx+'A3':A[3],pfx+'A4':A[4],pfx+'A5':A[5],pfx+'Vol':G2lat.calc_V(A)})
        if cell[0]:
            phaseVary += cellVary(pfx,SGData)
        #rigid body model input here 
        Natoms[pfx] = 0
        if Atoms and not General.get('doPawley'):
            cx,ct,cs,cia = General['AtomPtrs']
            if General['Type'] in ['nuclear','macromolecular']:
                Natoms[pfx] = len(Atoms)
                for i,at in enumerate(Atoms):
                    atomIndx[at[-1]] = [pfx,i]      #lookup table for restraints
                    phaseDict.update({pfx+'Atype:'+str(i):at[ct],pfx+'Afrac:'+str(i):at[cx+3],pfx+'Amul:'+str(i):at[cs+1],
                        pfx+'Ax:'+str(i):at[cx],pfx+'Ay:'+str(i):at[cx+1],pfx+'Az:'+str(i):at[cx+2],
                        pfx+'dAx:'+str(i):0.,pfx+'dAy:'+str(i):0.,pfx+'dAz:'+str(i):0.,         #refined shifts for x,y,z
                        pfx+'AI/A:'+str(i):at[cia],})
                    if at[cia] == 'I':
                        phaseDict[pfx+'AUiso:'+str(i)] = at[cia+1]
                    else:
                        phaseDict.update({pfx+'AU11:'+str(i):at[cia+2],pfx+'AU22:'+str(i):at[cia+3],pfx+'AU33:'+str(i):at[cia+4],
                            pfx+'AU12:'+str(i):at[cia+5],pfx+'AU13:'+str(i):at[cia+6],pfx+'AU23:'+str(i):at[cia+7]})
                    if 'F' in at[ct+1]:
                        phaseVary.append(pfx+'Afrac:'+str(i))
                    if 'X' in at[ct+1]:
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
                                for eqv in equiv[1:]:
                                    G2mv.StoreEquivalence(name,(eqv,))
                    if 'U' in at[ct+1]:
                        if at[9] == 'I':
                            phaseVary.append(pfx+'AUiso:'+str(i))
                        else:
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
                                    for eqv in equiv[1:]:
                                        G2mv.StoreEquivalence(name,(eqv,))
#            elif General['Type'] == 'magnetic':
#            elif General['Type'] == 'macromolecular':
            textureData = General['SH Texture']
            if textureData['Order']:
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
                print >>pFile,'\n Phase name: ',General['Name']
                print >>pFile,135*'-'
                PrintFFtable(FFtable)
                PrintBLtable(BLtable)
                print >>pFile,''
                for line in SGtext: print >>pFile,line
                #PrintRBObjects(whatever is needed here)
                PrintAtoms(General,Atoms)
                print >>pFile,'\n Unit cell: a =','%.5f'%(cell[1]),' b =','%.5f'%(cell[2]),' c =','%.5f'%(cell[3]), \
                    ' alpha =','%.3f'%(cell[4]),' beta =','%.3f'%(cell[5]),' gamma =', \
                    '%.3f'%(cell[6]),' volume =','%.3f'%(cell[7]),' Refine?',cell[0]
                PrintTexture(textureData)
                if name in RestraintDict:
                    PrintRestraints(cell[1:7],SGData,General['AtomPtrs'],Atoms,AtLookup,
                        textureData,RestraintDict[name],pFile)
                    
        elif PawleyRef:
            pawleyVary = []
            for i,refl in enumerate(PawleyRef):
                phaseDict[pfx+'PWLref:'+str(i)] = refl[6]
                pawleyLookup[pfx+'%d,%d,%d'%(refl[0],refl[1],refl[2])] = i
                if refl[5]:
                    pawleyVary.append(pfx+'PWLref:'+str(i))
            GetPawleyConstr(SGData['SGLaue'],PawleyRef,pawleyVary)      #does G2mv.StoreEquivalence
            phaseVary += pawleyVary
                
    return Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables
    
def cellFill(pfx,SGData,parmDict,sigDict): 
    if SGData['SGLaue'] in ['-1',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
            parmDict[pfx+'A3'],parmDict[pfx+'A4'],parmDict[pfx+'A5']]
        sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
            sigDict[pfx+'A3'],sigDict[pfx+'A4'],sigDict[pfx+'A5']]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                parmDict[pfx+'A3'],0,0]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                sigDict[pfx+'A3'],0,0]
        elif SGData['SGUniq'] == 'b':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,parmDict[pfx+'A4'],0]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                0,sigDict[pfx+'A4'],0]
        else:
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,0,parmDict[pfx+'A5']]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                0,0,sigDict[pfx+'A5']]
    elif SGData['SGLaue'] in ['mmm',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],0,0,0]
        sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],0,0,0]
        sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],
            parmDict[pfx+'A0'],0,0]
        sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],
            parmDict[pfx+'A3'],parmDict[pfx+'A3'],parmDict[pfx+'A3']]
        sigA = [sigDict[pfx+'A0'],0,0,sigDict[pfx+'A3'],0,0]
    elif SGData['SGLaue'] in ['m3m','m3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],0,0,0]
        sigA = [sigDict[pfx+'A0'],0,0,0,0,0]
    return A,sigA
        
def PrintRestraints(cell,SGData,AtPtrs,Atoms,AtLookup,textureData,phaseRest,pFile):
    if phaseRest:
        Amat = G2lat.cell2AB(cell)[0]
        cx,ct,cs = AtPtrs[:3]
        names = [['Bond','Bonds'],['Angle','Angles'],['Plane','Planes'],
            ['Chiral','Volumes'],['Torsion','Torsions'],['Rama','Ramas'],
            ['ChemComp','Sites'],['Texture','HKLs']]
        for name,rest in names:
            itemRest = phaseRest[name]
            if itemRest[rest] and itemRest['Use']:
                print >>pFile,'\n  %s %10.3f Use: %s'%(name+' restraint weight factor',itemRest['wtFactor'],str(itemRest['Use']))
                if name in ['Bond','Angle','Plane','Chiral']:
                    print >>pFile,'     calc       obs      sig   delt/sig  atoms(symOp): '
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
                            print >>pFile,' %9.3f %9.3f %8.3f %8.3f   %s'%(calc,obs,esd,(obs-calc)/esd,AtName[:-1])
                        except KeyError:
                            del itemRest[rest]
                elif name in ['Torsion','Rama']:
                    print >>pFile,'  atoms(symOp)  calc  obs  sig  delt/sig  torsions: '
                    coeffDict = itemRest['Coeff']
                    for indx,ops,cofName,esd in enumerate(itemRest[rest]):
                        AtNames = G2mth.GetAtomItemsById(Atoms,AtLookup,indx,ct-1)
                        AtName = ''
                        for i,Aname in enumerate(AtNames):
                            AtName += Aname+'+('+ops[i]+')-'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        if name == 'Torsion':
                            tor = G2mth.getRestTorsion(XYZ,Amat)
                            restr,calc = G2mth.calcTorsionEnergy(tor,coeffDict[cofName])
                            print >>pFile,' %8.3f %8.3f %.3f %8.3f %8.3f %s'%(AtName[:-1],calc,obs,esd,(obs-calc)/esd,tor)
                        else:
                            phi,psi = G2mth.getRestRama(XYZ,Amat)
                            restr,calc = G2mth.calcRamaEnergy(phi,psi,coeffDict[cofName])                               
                            print >>pFile,' %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %s'%(AtName[:-1],calc,obs,esd,(obs-calc)/esd,phi,psi)
                elif name == 'ChemComp':
                    print >>pFile,'     atoms   mul*frac  factor     prod'
                    for indx,factors,obs,esd in itemRest[rest]:
                        try:
                            atoms = G2mth.GetAtomItemsById(Atoms,AtLookup,indx,ct-1)
                            mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs+1))
                            frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookup,indx,cs-1))
                            mulfrac = mul*frac
                            calcs = mul*frac*factors
                            for iatm,[atom,mf,fr,clc] in enumerate(zip(atoms,mulfrac,factors,calcs)):
                                print >>pFile,' %10s %8.3f %8.3f %8.3f'%(atom,mf,fr,clc)
                            print >>pFile,' Sum:                   calc: %8.3f obs: %8.3f esd: %8.3f'%(np.sum(calcs),obs,esd)
                        except KeyError:
                            del itemRest[rest]
                elif name == 'Texture' and textureData['Order']:
                    Start = False
                    SHCoef = textureData['SH Coeff'][1]
                    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
                    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
                    print '    HKL  grid  neg esd  sum neg  num neg use unit?  unit esd  '
                    for hkl,grid,esd1,ifesd2,esd2 in itemRest[rest]:
                        PH = np.array(hkl)
                        phi,beta = G2lat.CrsAng(np.array(hkl),cell,SGData)
                        ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
                        R,P,Z = G2mth.getRestPolefig(ODFln,SamSym[textureData['Model']],grid)
                        Z = ma.masked_greater(Z,0.0)
                        num = ma.count(Z)
                        sum = 0
                        if num:
                            sum = np.sum(Z)
                        print '   %d %d %d  %d %8.3f %8.3f %8d   %s    %8.3f'%(hkl[0],hkl[1],hkl[2],grid,esd1,sum,num,str(ifesd2),esd2)
        
def getCellEsd(pfx,SGData,A,covData):
    dpr = 180./np.pi
    rVsq = G2lat.calc_rVsq(A)
    G,g = G2lat.A2Gmat(A)       #get recip. & real metric tensors
    cell = np.array(G2lat.Gmat2cell(g))   #real cell
    cellst = np.array(G2lat.Gmat2cell(G)) #recip. cell
    scos = cosd(cellst[3:6])
    ssin = sind(cellst[3:6])
    scot = scos/ssin
    rcos = cosd(cell[3:6])
    rsin = sind(cell[3:6])
    rcot = rcos/rsin
    RMnames = [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    varyList = covData['varyList']
    covMatrix = covData['covMatrix']
    vcov = G2mth.getVCov(RMnames,varyList,covMatrix)
    Ax = np.array(A)
    Ax[3:] /= 2.
    drVdA = np.array([Ax[1]*Ax[2]-Ax[5]**2,Ax[0]*Ax[2]-Ax[4]**2,Ax[0]*Ax[1]-Ax[3]**2,
        Ax[4]*Ax[5]-Ax[2]*Ax[3],Ax[3]*Ax[5]-Ax[1]*Ax[4],Ax[3]*Ax[4]-Ax[0]*Ax[5]])
    srcvlsq = np.inner(drVdA,np.inner(vcov,drVdA.T))
    Vol = 1/np.sqrt(rVsq)
    sigVol = Vol**3*np.sqrt(srcvlsq)/2.
    R123 = Ax[0]*Ax[1]*Ax[2]
    dsasdg = np.zeros((3,6))
    dadg = np.zeros((6,6))
    for i0 in range(3):         #0  1   2
        i1 = (i0+1)%3           #1  2   0
        i2 = (i1+1)%3           #2  0   1
        i3 = 5-i2               #3  5   4
        i4 = 5-i1               #4  3   5
        i5 = 5-i0               #5  4   3
        dsasdg[i0][i1] = 0.5*scot[i0]*scos[i0]/Ax[i1]
        dsasdg[i0][i2] = 0.5*scot[i0]*scos[i0]/Ax[i2]
        dsasdg[i0][i5] = -scot[i0]/np.sqrt(Ax[i1]*Ax[i2])
        denmsq = Ax[i0]*(R123-Ax[i1]*Ax[i4]**2-Ax[i2]*Ax[i3]**2+(Ax[i4]*Ax[i3])**2)
        denom = np.sqrt(denmsq)
        dadg[i5][i0] = -Ax[i5]/denom-rcos[i0]/denmsq*(R123-0.5*Ax[i1]*Ax[i4]**2-0.5*Ax[i2]*Ax[i3]**2)
        dadg[i5][i1] = -0.5*rcos[i0]/denmsq*(Ax[i0]**2*Ax[i2]-Ax[i0]*Ax[i4]**2)
        dadg[i5][i2] = -0.5*rcos[i0]/denmsq*(Ax[i0]**2*Ax[i1]-Ax[i0]*Ax[i3]**2)
        dadg[i5][i3] = Ax[i4]/denom+rcos[i0]/denmsq*(Ax[i0]*Ax[i2]*Ax[i3]-Ax[i3]*Ax[i4]**2)
        dadg[i5][i4] = Ax[i3]/denom+rcos[i0]/denmsq*(Ax[i0]*Ax[i1]*Ax[i4]-Ax[i3]**2*Ax[i4])
        dadg[i5][i5] = -Ax[i0]/denom
    for i0 in range(3):
        i1 = (i0+1)%3
        i2 = (i1+1)%3
        i3 = 5-i2
        for ij in range(6):
            dadg[i0][ij] = cell[i0]*(rcot[i2]*dadg[i3][ij]/rsin[i2]-dsasdg[i1][ij]/ssin[i1])
            if ij == i0:
                dadg[i0][ij] = dadg[i0][ij]-0.5*cell[i0]/Ax[i0]
            dadg[i3][ij] = -dadg[i3][ij]*rsin[2-i0]*dpr
    sigMat = np.inner(dadg,np.inner(vcov,dadg.T))
    var = np.diag(sigMat)
    CS = np.where(var>0.,np.sqrt(var),0.)
    cellSig = [CS[0],CS[1],CS[2],CS[5],CS[4],CS[3],sigVol]  #exchange sig(alp) & sig(gam) to get in right order
    return cellSig            
    
def SetPhaseData(parmDict,sigDict,Phases,covData,RestraintDict=None,pFile=None):
    
    def PrintAtomsAndSig(General,Atoms,atomsSig):
        print >>pFile,'\n Atoms:'
        line = '   name      x         y         z      frac   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'magnetic':
            line += '   Mx     My     Mz'
        elif General['Type'] == 'macromolecular':
            line = ' res no  residue  chain '+line
        print >>pFile,line
        if General['Type'] == 'nuclear':
            print >>pFile,135*'-'
            fmt = {0:'%7s',1:'%7s',3:'%10.5f',4:'%10.5f',5:'%10.5f',6:'%8.3f',10:'%8.5f',
                11:'%8.5f',12:'%8.5f',13:'%8.5f',14:'%8.5f',15:'%8.5f',16:'%8.5f'}
            noFXsig = {3:[10*' ','%10s'],4:[10*' ','%10s'],5:[10*' ','%10s'],6:[8*' ','%8s']}
            for atyp in General['AtomTypes']:       #zero composition
                General['NoAtoms'][atyp] = 0.
            for i,at in enumerate(Atoms):
                General['NoAtoms'][at[1]] += at[6]*float(at[8])     #new composition
                name = fmt[0]%(at[0])+fmt[1]%(at[1])+':'
                valstr = ' values:'
                sigstr = ' sig   :'
                for ind in [3,4,5,6]:
                    sigind = str(i)+':'+str(ind)
                    valstr += fmt[ind]%(at[ind])                    
                    if sigind in atomsSig:
                        sigstr += fmt[ind]%(atomsSig[sigind])
                    else:
                        sigstr += noFXsig[ind][1]%(noFXsig[ind][0])
                if at[9] == 'I':
                    valstr += fmt[10]%(at[10])
                    if str(i)+':10' in atomsSig:
                        sigstr += fmt[10]%(atomsSig[str(i)+':10'])
                    else:
                        sigstr += 8*' '
                else:
                    valstr += 8*' '
                    sigstr += 8*' '
                    for ind in [11,12,13,14,15,16]:
                        sigind = str(i)+':'+str(ind)
                        valstr += fmt[ind]%(at[ind])
                        if sigind in atomsSig:                        
                            sigstr += fmt[ind]%(atomsSig[sigind])
                        else:
                            sigstr += 8*' '
                print >>pFile,name
                print >>pFile,valstr
                print >>pFile,sigstr
                
    #def PrintRBObjectsAndSig()
                
    def PrintSHtextureAndSig(textureData,SHtextureSig):
        print >>pFile,'\n Spherical harmonics texture: Order:' + str(textureData['Order'])
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
        print >>pFile,namstr
        print >>pFile,ptstr
        print >>pFile,sigstr
        print >>pFile,'\n Texture coefficients:'
        namstr = '  names :'
        ptstr =  '  values:'
        sigstr = '  esds  :'
        SHcoeff = textureData['SH Coeff'][1]
        for name in SHcoeff:
            namstr += '%12s'%(name)
            ptstr += '%12.3f'%(SHcoeff[name])
            if name in SHtextureSig:
                sigstr += '%12.3f'%(SHtextureSig[name])
            else:
                sigstr += 12*' '
        print >>pFile,namstr
        print >>pFile,ptstr
        print >>pFile,sigstr
        
            
    print >>pFile,'\n Phases:'
    for phase in Phases:
        print >>pFile,' Result for phase: ',phase
        Phase = Phases[phase]
        General = Phase['General']
        SGData = General['SGData']
        Atoms = Phase['Atoms']
        AtLookup = G2mth.FillAtomLookUp(Atoms)
        cell = General['Cell']
        pId = Phase['pId']
        pfx = str(pId)+'::'
        if cell[0]:
            A,sigA = cellFill(pfx,SGData,parmDict,sigDict)
            cellSig = getCellEsd(pfx,SGData,A,covData)  #includes sigVol
            print >>pFile,' Reciprocal metric tensor: '
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
            print >>pFile,namstr
            print >>pFile,ptstr
            print >>pFile,sigstr
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            print >>pFile,' New unit cell:'
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
            print >>pFile,namstr
            print >>pFile,ptstr
            print >>pFile,sigstr
            
        General['Mass'] = 0.
        if Phase['General'].get('doPawley'):
            pawleyRef = Phase['Pawley ref']
            for i,refl in enumerate(pawleyRef):
                key = pfx+'PWLref:'+str(i)
                refl[6] = parmDict[key]
                if key in sigDict:
                    refl[7] = sigDict[key]
                else:
                    refl[7] = 0
        else:
            #new RBObj parms here, mod atoms, & PrintRBObjectsAndSig
            atomsSig = {}
            if General['Type'] == 'nuclear':        #this needs macromolecular variant!
                for i,at in enumerate(Atoms):
                    names = {3:pfx+'Ax:'+str(i),4:pfx+'Ay:'+str(i),5:pfx+'Az:'+str(i),6:pfx+'Afrac:'+str(i),
                        10:pfx+'AUiso:'+str(i),11:pfx+'AU11:'+str(i),12:pfx+'AU22:'+str(i),13:pfx+'AU33:'+str(i),
                        14:pfx+'AU12:'+str(i),15:pfx+'AU13:'+str(i),16:pfx+'AU23:'+str(i)}
                    for ind in [3,4,5,6]:
                        at[ind] = parmDict[names[ind]]
                        if ind in [3,4,5]:
                            name = names[ind].replace('A','dA')
                        else:
                            name = names[ind]
                        if name in sigDict:
                            atomsSig[str(i)+':'+str(ind)] = sigDict[name]
                    if at[9] == 'I':
                        at[10] = parmDict[names[10]]
                        if names[10] in sigDict:
                            atomsSig[str(i)+':10'] = sigDict[names[10]]
                    else:
                        for ind in [11,12,13,14,15,16]:
                            at[ind] = parmDict[names[ind]]
                            if names[ind] in sigDict:
                                atomsSig[str(i)+':'+str(ind)] = sigDict[names[ind]]
                    ind = General['AtomTypes'].index(at[1])
                    General['Mass'] += General['AtomMass'][ind]*at[6]*at[8]
            PrintAtomsAndSig(General,Atoms,atomsSig)
        
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
        if phase in RestraintDict:
            PrintRestraints(cell[1:7],SGData,General['AtomPtrs'],Atoms,AtLookup,
                textureData,RestraintDict[phase],pFile)
                    
################################################################################
##### Histogram & Phase data
################################################################################        
                    
def GetHistogramPhaseData(Phases,Histograms,Print=True,pFile=None):
    
    def PrintSize(hapData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Size model    : %9s'%(hapData[0])
            line += ' equatorial:'+'%12.3f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.3f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            line += '\n\t LG mixing coeff.: %12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            print >>pFile,line
        else:
            print >>pFile,'\n Size model    : %s'%(hapData[0])+ \
                '\n\t LG mixing coeff.:%12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print >>pFile,ptlbls
            print >>pFile,ptstr
            print >>pFile,varstr
        
    def PrintMuStrain(hapData,SGData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Mustrain model: %9s'%(hapData[0])
            line += ' equatorial:'+'%12.1f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.1f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            line +='\n\t LG mixing coeff.:%12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            print >>pFile,line
        else:
            print >>pFile,'\n Mustrain model: %s'%(hapData[0])+ \
                '\n\t LG mixing coeff.:%12.4f'%(hapData[1][2])+' Refine? '+str(hapData[2][2])
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print >>pFile,ptlbls
            print >>pFile,ptstr
            print >>pFile,varstr

    def PrintHStrain(hapData,SGData):
        print >>pFile,'\n Hydrostatic/elastic strain: '
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' names :'
        ptstr =  ' values:'
        varstr = ' refine:'
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (hapData[0][i])
            varstr += '%12s' % (str(hapData[1][i]))
        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,varstr

    def PrintSHPO(hapData):
        print >>pFile,'\n Spherical harmonics preferred orientation: Order:' + \
            str(hapData[4])+' Refine? '+str(hapData[2])
        ptlbls = ' names :'
        ptstr =  ' values:'
        for item in hapData[5]:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[5][item]) 
        print >>pFile,ptlbls
        print >>pFile,ptstr
    
    def PrintBabinet(hapData):
        print >>pFile,'\n Babinet form factor modification: '
        ptlbls = ' names :'
        ptstr =  ' values:'
        varstr = ' refine:'
        for name in ['BabA','BabU']:
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (hapData[name][0])
            varstr += '%12s' % (str(hapData[name][1]))
        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,varstr
        
    
    hapDict = {}
    hapVary = []
    controlDict = {}
    poType = {}
    poAxes = {}
    spAxes = {}
    spType = {}
    
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        SGData = Phases[phase]['General']['SGData']
        cell = Phases[phase]['General']['Cell'][1:7]
        A = G2lat.cell2A(cell)
        pId = Phases[phase]['pId']
        histoList = HistoPhase.keys()
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            if 'PWDR' in histogram:
                limits = Histogram['Limits'][1]
                inst = Histogram['Instrument Parameters'][0]
                Zero = inst['Zero'][1]
                if 'C' in inst['Type'][1]:
                    try:
                        wave = inst['Lam'][1]
                    except KeyError:
                        wave = inst['Lam1'][1]
                    dmin = wave/(2.0*sind(limits[1]/2.0))
                pfx = str(pId)+':'+str(hId)+':'
                for item in ['Scale','Extinction']:
                    hapDict[pfx+item] = hapData[item][0]
                    if hapData[item][1]:
                        hapVary.append(pfx+item)
                names = G2spc.HStrainNames(SGData)
                for i,name in enumerate(names):
                    hapDict[pfx+name] = hapData['HStrain'][0][i]
                    if hapData['HStrain'][1][i]:
                        hapVary.append(pfx+name)
                controlDict[pfx+'poType'] = hapData['Pref.Ori.'][0]
                if hapData['Pref.Ori.'][0] == 'MD':
                    hapDict[pfx+'MD'] = hapData['Pref.Ori.'][1]
                    controlDict[pfx+'MDAxis'] = hapData['Pref.Ori.'][3]
                    if hapData['Pref.Ori.'][2]:
                        hapVary.append(pfx+'MD')
                else:                           #'SH' spherical harmonics
                    controlDict[pfx+'SHord'] = hapData['Pref.Ori.'][4]
                    controlDict[pfx+'SHncof'] = len(hapData['Pref.Ori.'][5])
                    for item in hapData['Pref.Ori.'][5]:
                        hapDict[pfx+item] = hapData['Pref.Ori.'][5][item]
                        if hapData['Pref.Ori.'][2]:
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
                            sfx = ':'+str(i)
                            hapDict[pfx+item+sfx] = hapData[item][4][i]
                            if hapData[item][5][i]:
                                hapVary.append(pfx+item+sfx)
                for bab in ['BabA','BabU']:
                    hapDict[pfx+bab] = hapData['Babinet'][bab][0]
                    if hapData['Babinet'][bab][1]:
                        hapVary.append(pfx+bab)
                                
                if Print: 
                    print >>pFile,'\n Phase: ',phase,' in histogram: ',histogram
                    print >>pFile,135*'-'
                    print >>pFile,' Phase fraction  : %10.4f'%(hapData['Scale'][0]),' Refine?',hapData['Scale'][1]
                    print >>pFile,' Extinction coeff: %10.4f'%(hapData['Extinction'][0]),' Refine?',hapData['Extinction'][1]
                    if hapData['Pref.Ori.'][0] == 'MD':
                        Ax = hapData['Pref.Ori.'][3]
                        print >>pFile,' March-Dollase PO: %10.4f'%(hapData['Pref.Ori.'][1]),' Refine?',hapData['Pref.Ori.'][2], \
                            ' Axis: %d %d %d'%(Ax[0],Ax[1],Ax[2])
                    else: #'SH' for spherical harmonics
                        PrintSHPO(hapData['Pref.Ori.'])
                    PrintSize(hapData['Size'])
                    PrintMuStrain(hapData['Mustrain'],SGData)
                    PrintHStrain(hapData['HStrain'],SGData)
                    PrintBabinet(hapData['Babinet'])
                HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
                refList = []
                for h,k,l,d in HKLd:
                    ext,mul,Uniq,phi = G2spc.GenHKLf([h,k,l],SGData)
                    mul *= 2      # for powder overlap of Friedel pairs
                    if ext:
                        continue
                    if 'C' in inst['Type'][0]:
                        pos = 2.0*asind(wave/(2.0*d))+Zero
                        if limits[0] < pos < limits[1]:
                            refList.append([h,k,l,mul,d,pos,0.0,0.0,0.0,0.0,0.0,Uniq,phi,0.0,{}])
                            #last item should contain form factor values by atom type
                    else:
                        raise ValueError 
                Histogram['Reflection Lists'][phase] = refList
            elif 'HKLF' in histogram:
                inst = Histogram['Instrument Parameters'][0]
                hId = Histogram['hId']
                hfx = ':%d:'%(hId)
                for item in inst:
                    hapDict[hfx+item] = inst[item][1]
                pfx = str(pId)+':'+str(hId)+':'
                hapDict[pfx+'Scale'] = hapData['Scale'][0]
                if hapData['Scale'][1]:
                    hapVary.append(pfx+'Scale')
                                
                extApprox,extType,extParms = hapData['Extinction']
                controlDict[pfx+'EType'] = extType
                controlDict[pfx+'EApprox'] = extApprox
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
                if Print: 
                    print >>pFile,'\n Phase: ',phase,' in histogram: ',histogram
                    print >>pFile,135*'-'
                    print >>pFile,' Scale factor     : %10.4f'%(hapData['Scale'][0]),' Refine?',hapData['Scale'][1]
                    if extType != 'None':
                        print >>pFile,' Extinction  Type: %15s'%(extType),' approx: %10s'%(extApprox),' tbar: %6.3f'%(extParms['Tbar'])
                        text = ' Parameters       :'
                        for eKey in Ekey:
                            text += ' %4s : %10.3e Refine? '%(eKey,extParms[eKey][0])+str(extParms[eKey][1])
                        print >>pFile,text
                    PrintBabinet(hapData['Babinet'])
                Histogram['Reflection Lists'] = phase       
                
    return hapVary,hapDict,controlDict
    
def SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,Print=True,pFile=None):
    
    def PrintSizeAndSig(hapData,sizeSig):
        line = '\n Size model:     %9s'%(hapData[0])
        refine = False
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.4f'%(hapData[1][0])
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
                print >>pFile,line
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
                print >>pFile,line
                print >>pFile,ptlbls
                print >>pFile,ptstr
                print >>pFile,sigstr
        
    def PrintMuStrainAndSig(hapData,mustrainSig,SGData):
        line = '\n Mustrain model: %9s'%(hapData[0])
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
                print >>pFile,line
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
                ptstr += '%12.6f' % (hapData[4][i])
                if mustrainSig[1][i]:
                    refine = True
                    sigstr += '%12.6f' % (mustrainSig[1][i])
                else:
                    sigstr += 12*' '
            if refine:
                print >>pFile,line
                print >>pFile,ptlbls
                print >>pFile,ptstr
                print >>pFile,sigstr
            
    def PrintHStrainAndSig(hapData,strainSig,SGData):
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' name  :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        refine = False
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6g' % (hapData[0][i])
            if name in strainSig:
                refine = True
                sigstr += '%12.6g' % (strainSig[name])
            else:
                sigstr += 12*' '
        if refine:
            print >>pFile,'\n Hydrostatic/elastic strain: '
            print >>pFile,ptlbls
            print >>pFile,ptstr
            print >>pFile,sigstr
        
    def PrintSHPOAndSig(pfx,hapData,POsig):
        print >>pFile,'\n Spherical harmonics preferred orientation: Order:'+str(hapData[4])
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
        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,sigstr
        
    def PrintExtAndSig(pfx,hapData,ScalExtSig):
        print >>pFile,'\n Single crystal extinction: Type: ',hapData[0],' Approx: ',hapData[1]
        text = ''
        for item in hapData[2]:
            if pfx+item in ScalExtSig:
                text += '       %s: '%(item)
                text += '%12.2e'%(hapData[2][item][0])
                if pfx+item in ScalExtSig:
                    text += ' sig: %12.2e'%(ScalExtSig[pfx+item])
        print >>pFile,text        
        
    def PrintBabinetAndSig(pfx,hapData,BabSig):
        print >>pFile,'\n Babinet form factor modification: '
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
        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,sigstr
    
    PhFrExtPOSig = {}
    SizeMuStrSig = {}
    ScalExtSig = {}
    BabSig = {}
    wtFrSum = {}
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        General = Phases[phase]['General']
        SGData = General['SGData']
        pId = Phases[phase]['pId']
        histoList = HistoPhase.keys()
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            pfx = str(pId)+':'+str(hId)+':'
            if hId not in wtFrSum:
                wtFrSum[hId] = 0.
            if 'PWDR' in histogram:
                for item in ['Scale','Extinction']:
                    hapData[item][0] = parmDict[pfx+item]
                    if pfx+item in sigDict:
                        PhFrExtPOSig.update({pfx+item:sigDict[pfx+item],})
                wtFrSum[hId] += hapData['Scale'][0]*General['Mass']
                if hapData['Pref.Ori.'][0] == 'MD':
                    hapData['Pref.Ori.'][1] = parmDict[pfx+'MD']
                    if pfx+'MD' in sigDict:
                        PhFrExtPOSig.update({pfx+'MD':sigDict[pfx+'MD'],})
                else:                           #'SH' spherical harmonics
                    for item in hapData['Pref.Ori.'][5]:
                        hapData['Pref.Ori.'][5][item] = parmDict[pfx+item]
                        if pfx+item in sigDict:
                            PhFrExtPOSig.update({pfx+item:sigDict[pfx+item],})
                SizeMuStrSig.update({pfx+'Mustrain':[[0,0,0],[0 for i in range(len(hapData['Mustrain'][4]))]],
                    pfx+'Size':[[0,0,0],[0 for i in range(len(hapData['Size'][4]))]],
                    pfx+'HStrain':{}})                  
                for item in ['Mustrain','Size']:
                    hapData[item][1][2] = parmDict[pfx+item+';mx']
                    hapData[item][1][2] = min(1.,max(0.1,hapData[item][1][2]))
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
                            sfx = ':'+str(i)
                            hapData[item][4][i] = parmDict[pfx+item+sfx]
                            if pfx+item+sfx in sigDict:
                                SizeMuStrSig[pfx+item][1][i] = sigDict[pfx+item+sfx]
                names = G2spc.HStrainNames(SGData)
                for i,name in enumerate(names):
                    hapData['HStrain'][0][i] = parmDict[pfx+name]
                    if pfx+name in sigDict:
                        SizeMuStrSig[pfx+'HStrain'][name] = sigDict[pfx+name]
                for name in ['BabA','BabU']:
                    hapData['Babinet'][name][0] = parmDict[pfx+name]
                    if pfx+name in sigDict:
                        BabSig[pfx+name] = sigDict[pfx+name]                
                
            elif 'HKLF' in histogram:
                for item in ['Scale',]:
                    if parmDict.get(pfx+item):
                        hapData[item][0] = parmDict[pfx+item]
                        if pfx+item in sigDict:
                            ScalExtSig[pfx+item] = sigDict[pfx+item]
                for item in ['Ep','Eg','Es']:
                    if parmDict.get(pfx+item):
                        hapData['Extinction'][2][item][0] = parmDict[pfx+item]
                        if pfx+item in sigDict:
                            ScalExtSig[pfx+item] = sigDict[pfx+item]
                for name in ['BabA','BabU']:
                    hapData['Babinet'][name][0] = parmDict[pfx+name]
                    if pfx+name in sigDict:
                        BabSig[pfx+name] = sigDict[pfx+name]                

    if Print:
        for phase in Phases:
            HistoPhase = Phases[phase]['Histograms']
            General = Phases[phase]['General']
            SGData = General['SGData']
            pId = Phases[phase]['pId']
            histoList = HistoPhase.keys()
            histoList.sort()
            for histogram in histoList:
                try:
                    Histogram = Histograms[histogram]
                except KeyError:                        
                    #skip if histogram not included e.g. in a sequential refinement
                    continue
                print >>pFile,'\n Phase: ',phase,' in histogram: ',histogram
                print >>pFile,130*'-'
                hapData = HistoPhase[histogram]
                hId = Histogram['hId']
                pfx = str(pId)+':'+str(hId)+':'
                if 'PWDR' in histogram:
                    print >>pFile,' Final refinement RF, RF^2 = %.2f%%, %.2f%% on %d reflections'   \
                        %(Histogram[pfx+'Rf'],Histogram[pfx+'Rf^2'],Histogram[pfx+'Nref'])
                
                    if pfx+'Scale' in PhFrExtPOSig:
                        wtFr = hapData['Scale'][0]*General['Mass']/wtFrSum[hId]
                        sigwtFr = PhFrExtPOSig[pfx+'Scale']*wtFr/hapData['Scale'][0]
                        print >>pFile,' Phase fraction  : %10.5f, sig %10.5f Weight fraction  : %8.5f, sig %10.5f' \
                            %(hapData['Scale'][0],PhFrExtPOSig[pfx+'Scale'],wtFr,sigwtFr)
                    if pfx+'Extinction' in PhFrExtPOSig:
                        print >>pFile,' Extinction coeff: %10.4f, sig %10.4f'%(hapData['Extinction'][0],PhFrExtPOSig[pfx+'Extinction'])
                    if hapData['Pref.Ori.'][0] == 'MD':
                        if pfx+'MD' in PhFrExtPOSig:
                            print >>pFile,' March-Dollase PO: %10.4f, sig %10.4f'%(hapData['Pref.Ori.'][1],PhFrExtPOSig[pfx+'MD'])
                    else:
                        PrintSHPOAndSig(pfx,hapData['Pref.Ori.'],PhFrExtPOSig)
                    PrintSizeAndSig(hapData['Size'],SizeMuStrSig[pfx+'Size'])
                    PrintMuStrainAndSig(hapData['Mustrain'],SizeMuStrSig[pfx+'Mustrain'],SGData)
                    PrintHStrainAndSig(hapData['HStrain'],SizeMuStrSig[pfx+'HStrain'],SGData)
                    PrintBabinetAndSig(pfx,hapData['Babinet'],BabSig)
                    
                elif 'HKLF' in histogram:
                    print >>pFile,' Final refinement RF, RF^2 = %.2f%%, %.2f%% on %d reflections'   \
                        %(Histogram[pfx+'Rf'],Histogram[pfx+'Rf^2'],Histogram[pfx+'Nref'])
                    print >>pFile,' HKLF histogram weight factor = ','%.3f'%(Histogram['wtFactor'])
                    if pfx+'Scale' in ScalExtSig:
                        print >>pFile,' Scale factor : %10.4f, sig %10.4f'%(hapData['Scale'][0],ScalExtSig[pfx+'Scale'])
                    if hapData['Extinction'][0] != 'None':
                        PrintExtAndSig(pfx,hapData['Extinction'],ScalExtSig)
                    if len(BabSig):
                        PrintBabinetAndSig(pfx,hapData['Babinet'],BabSig)

################################################################################
##### Histogram data
################################################################################        
                    
def GetHistogramData(Histograms,Print=True,pFile=None):
    
    def GetBackgroundParms(hId,Background):
        Back = Background[0]
        DebyePeaks = Background[1]
        bakType,bakFlag = Back[:2]
        backVals = Back[3:]
        backNames = [':'+str(hId)+':Back:'+str(i) for i in range(len(backVals))]
        backDict = dict(zip(backNames,backVals))
        backVary = []
        if bakFlag:
            backVary = backNames
        backDict[':'+str(hId)+':nDebye'] = DebyePeaks['nDebye']
        backDict[':'+str(hId)+':nPeaks'] = DebyePeaks['nPeaks']
        debyeDict = {}
        debyeList = []
        for i in range(DebyePeaks['nDebye']):
            debyeNames = [':'+str(hId)+':DebyeA:'+str(i),':'+str(hId)+':DebyeR:'+str(i),':'+str(hId)+':DebyeU:'+str(i)]
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
            peakNames = [':'+str(hId)+':BkPkpos:'+str(i),':'+str(hId)+ \
                ':BkPkint:'+str(i),':'+str(hId)+':BkPksig:'+str(i),':'+str(hId)+':BkPkgam:'+str(i)]
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
        dataType = Inst['Type'][0]
        instDict = {}
        insVary = []
        pfx = ':'+str(hId)+':'
        insKeys = Inst.keys()
        insKeys.sort()
        for item in insKeys:
            insName = pfx+item
            instDict[insName] = Inst[item][1]
            if Inst[item][2]:
                insVary.append(insName)
        instDict[pfx+'X'] = max(instDict[pfx+'X'],0.001)
        instDict[pfx+'Y'] = max(instDict[pfx+'Y'],0.001)
        instDict[pfx+'SH/L'] = max(instDict[pfx+'SH/L'],0.0005)
        return dataType,instDict,insVary
        
    def GetSampleParms(hId,Sample):
        sampVary = []
        hfx = ':'+str(hId)+':'        
        sampDict = {hfx+'Gonio. radius':Sample['Gonio. radius'],hfx+'Omega':Sample['Omega'],
            hfx+'Chi':Sample['Chi'],hfx+'Phi':Sample['Phi']}
        Type = Sample['Type']
        if 'Bragg' in Type:             #Bragg-Brentano
            for item in ['Scale','Shift','Transparency']:       #surface roughness?, diffuse scattering?
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
        print >>pFile,'\n Background function: ',Back[0],' Refine?',bool(Back[1])
        line = ' Coefficients: '
        for i,back in enumerate(Back[3:]):
            line += '%10.3f'%(back)
            if i and not i%10:
                line += '\n'+15*' '
        print >>pFile,line
        if DebyePeaks['nDebye']:
            print >>pFile,'\n Debye diffuse scattering coefficients'
            parms = ['DebyeA','DebyeR','DebyeU']
            line = ' names :  '
            for parm in parms:
                line += '%8s refine?'%(parm)
            print >>pFile,line
            for j,term in enumerate(DebyePeaks['debyeTerms']):
                line = ' term'+'%2d'%(j)+':'
                for i in range(3):
                    line += '%10.3f %5s'%(term[2*i],bool(term[2*i+1]))                    
                print >>pFile,line
        if DebyePeaks['nPeaks']:
            print >>pFile,'\n Single peak coefficients'
            parms =    ['BkPkpos','BkPkint','BkPksig','BkPkgam']
            line = ' names :  '
            for parm in parms:
                line += '%8s refine?'%(parm)
            print >>pFile,line
            for j,term in enumerate(DebyePeaks['peaksList']):
                line = ' peak'+'%2d'%(j)+':'
                for i in range(4):
                    line += '%10.3f %5s'%(term[2*i],bool(term[2*i+1]))                    
                print >>pFile,line
        
    def PrintInstParms(Inst):
        print >>pFile,'\n Instrument Parameters:'
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        insKeys = Inst.keys()
        insKeys.sort()
        for item in insKeys:
            if item != 'Type':
                ptlbls += '%12s' % (item)
                ptstr += '%12.6f' % (Inst[item][1])
                if item in ['Lam1','Lam2','Azimuth']:
                    varstr += 12*' '
                else:
                    varstr += '%12s' % (str(bool(Inst[item][2])))
        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,varstr
        
    def PrintSampleParms(Sample):
        print >>pFile,'\n Sample Parameters:'
        print >>pFile,' Goniometer omega = %.2f, chi = %.2f, phi = %.2f'% \
            (Sample['Omega'],Sample['Chi'],Sample['Phi'])
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        if 'Bragg' in Sample['Type']:
            for item in ['Scale','Shift','Transparency']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))
            
        elif 'Debye' in Type:        #Debye-Scherrer
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))

        print >>pFile,ptlbls
        print >>pFile,ptstr
        print >>pFile,varstr
        
    histDict = {}
    histVary = []
    controlDict = {}
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            controlDict[pfx+'wtFactor'] = Histogram['wtFactor']
            controlDict[pfx+'Limits'] = Histogram['Limits'][1]
            
            Background = Histogram['Background']
            Type,bakDict,bakVary = GetBackgroundParms(hId,Background)
            controlDict[pfx+'bakType'] = Type
            histDict.update(bakDict)
            histVary += bakVary
            
            Inst = Histogram['Instrument Parameters'][0]
            Type,instDict,insVary = GetInstParms(hId,Inst)
            controlDict[pfx+'histType'] = Type
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
                print >>pFile,'\n Histogram: ',histogram,' histogram Id: ',hId
                print >>pFile,135*'-'
                Units = {'C':' deg','T':' msec'}
                units = Units[controlDict[pfx+'histType'][2]]
                Limits = controlDict[pfx+'Limits']
                print >>pFile,' Instrument type: ',Sample['Type']
                print >>pFile,' Histogram limits: %8.2f%s to %8.2f%s'%(Limits[0],units,Limits[1],units)     
                PrintSampleParms(Sample)
                PrintInstParms(Inst)
                PrintBackground(Background)
        elif 'HKLF' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            controlDict[pfx+'wtFactor'] =Histogram['wtFactor']
            Inst = Histogram['Instrument Parameters'][0]
            controlDict[pfx+'histType'] = Inst['Type'][0]
            histDict[pfx+'Lam'] = Inst['Lam'][1]
            controlDict[pfx+'keV'] = 12.397639/histDict[pfx+'Lam']                    
    return histVary,histDict,controlDict
    
def SetHistogramData(parmDict,sigDict,Histograms,Print=True,pFile=None):
    
    def SetBackgroundParms(pfx,Background,parmDict,sigDict):
        Back = Background[0]
        DebyePeaks = Background[1]
        lenBack = len(Back[3:])
        backSig = [0 for i in range(lenBack+3*DebyePeaks['nDebye']+4*DebyePeaks['nPeaks'])]
        for i in range(lenBack):
            Back[3+i] = parmDict[pfx+'Back:'+str(i)]
            if pfx+'Back:'+str(i) in sigDict:
                backSig[i] = sigDict[pfx+'Back:'+str(i)]
        if DebyePeaks['nDebye']:
            for i in range(DebyePeaks['nDebye']):
                names = [pfx+'DebyeA:'+str(i),pfx+'DebyeR:'+str(i),pfx+'DebyeU:'+str(i)]
                for j,name in enumerate(names):
                    DebyePeaks['debyeTerms'][i][2*j] = parmDict[name]
                    if name in sigDict:
                        backSig[lenBack+3*i+j] = sigDict[name]            
        if DebyePeaks['nPeaks']:
            for i in range(DebyePeaks['nPeaks']):
                names = [pfx+'BkPkpos:'+str(i),pfx+'BkPkint:'+str(i),
                    pfx+'BkPksig:'+str(i),pfx+'BkPkgam:'+str(i)]
                for j,name in enumerate(names):
                    DebyePeaks['peaksList'][i][2*j] = parmDict[name]
                    if name in sigDict:
                        backSig[lenBack+3*DebyePeaks['nDebye']+4*i+j] = sigDict[name]
        return backSig
        
    def SetInstParms(pfx,Inst,parmDict,sigDict):
        instSig = {}
        insKeys = Inst.keys()
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
            sampSig = [0 for i in range(3)]
            for i,item in enumerate(['Scale','Shift','Transparency']):       #surface roughness?, diffuse scattering?
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
        lenBack = len(Back[3:])
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
            print >>pFile,'\n Background function: ',Back[0]
            print >>pFile,valstr
            print >>pFile,sigstr 
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
                print >>pFile,'\n Debye diffuse scattering coefficients'
                print >>pFile,names
                print >>pFile,ptstr
                print >>pFile,sigstr
        if DebyePeaks['nPeaks']:
            ifAny = False
            ptfmt = "%14.3f"
            names =  ' names :'
            ptstr =  ' values:'
            sigstr = ' esds  :'
            for item in sigDict:
                if 'BkPk' in item:
                    ifAny = True
                    names += '%14s'%(item)
                    ptstr += ptfmt%(parmDict[item])
                    sigstr += ptfmt%(sigDict[item])
            if ifAny:
                print >>pFile,'\n Single peak coefficients'
                print >>pFile,names
                print >>pFile,ptstr
                print >>pFile,sigstr
        
    def PrintInstParmsSig(Inst,instSig):
        ptlbls = ' names :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        refine = False
        insKeys = instSig.keys()
        insKeys.sort()
        for name in insKeys:
            if name not in  ['Type','Lam1','Lam2','Azimuth']:
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (Inst[name][1])
                if instSig[name]:
                    refine = True
                    sigstr += '%12.6f' % (instSig[name])
                else:
                    sigstr += 12*' '
        if refine:
            print >>pFile,'\n Instrument Parameters:'
            print >>pFile,ptlbls
            print >>pFile,ptstr
            print >>pFile,sigstr
        
    def PrintSampleParmsSig(Sample,sampleSig):
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        refine = False
        if 'Bragg' in Sample['Type']:
            for i,item in enumerate(['Scale','Shift','Transparency']):
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
            print >>pFile,'\n Sample Parameters:'
            print >>pFile,ptlbls
            print >>pFile,ptstr
            print >>pFile,sigstr
        
    histoList = Histograms.keys()
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

            print >>pFile,'\n Histogram: ',histogram,' histogram Id: ',hId
            print >>pFile,135*'-'
            print >>pFile,' Final refinement wR = %.2f%% on %d observations in this histogram'%(Histogram['wR'],Histogram['Nobs'])
            print >>pFile,' PWDR histogram weight factor = '+'%.3f'%(Histogram['wtFactor'])
            if Print:
                print >>pFile,' Instrument type: ',Sample['Type']
                PrintSampleParmsSig(Sample,sampSig)
                PrintInstParmsSig(Inst,instSig)
                PrintBackgroundSig(Background,backSig)
                
################################################################################
##### Penalty & restraint functions 
################################################################################

def penaltyFxn(HistoPhases,parmDict,varyList):
    Histograms,Phases,restraintDict = HistoPhases
    pNames = []
    pVals = []
    pWt = []
    negWt = {}
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
                elif name in ['Torsion','Rama']:
                    coeffDict = itemRest['Coeff']
                    for i,[indx,ops,cofName,esd] in enumerate(torsionList):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        XYZ = np.array(G2mth.GetAtomCoordsByID(pId,parmDict,AtLookup,indx))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        if name == 'Torsion':
                            tor = G2mth.getRestTorsion(XYZ,Amat)
                            restr,calc = G2mth.calcTorsionEnergy(tor,coeffDict[cofName])
                        else:
                            phi,psi = G2mth.getRestRama(XYZ,Amat)
                            restr,calc = G2mth.calcRamaEnergy(phi,psi,coeffDict[cofName])                               
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)
                elif name == 'ChemComp':
                    for i,[indx,factors,obs,esd] in enumerate(itemRest[rest]):
                        pNames.append(str(pId)+':'+name+':'+str(i))
                        mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs+1))
                        frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs-1))
                        calc = np.sum(mul*frac*factors)
                        pVals.append(obs-calc)
                        pWt.append(wt/esd**2)                    
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
                        if ifesd2:
                            Z2 = 1.-Z
                            for ind in np.ndindex(grid,grid):
                                pNames.append('%d:%s:%d:%.2f:%.2f'%(pId,name+'-unit',i,R[ind[0],ind[1]],P[ind[0],ind[1]]))
                                pVals.append(Z1[ind[0]][ind[1]])
                                pWt.append(wt/esd2**2)
         
    for item in varyList:
        if 'PWLref' in item and parmDict[item] < 0.:
            pId = int(item.split(':')[0])
            if negWt[pId]:
                pNames.append(item)
                pVals.append(-parmDict[item])
                pWt.append(negWt[pId])
    pVals = np.array(pVals)
    pWt = np.array(pWt)
    return pNames,pVals,pWt
    
def penaltyDeriv(pNames,pVal,HistoPhases,parmDict,varyList):
    Histograms,Phases,restraintDict = HistoPhases
    pDerv = np.zeros((len(varyList),len(pVal)))
    for phase in Phases:
        if phase not in restraintDict:
            continue
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
        phaseRest = restraintDict[phase]
        names = {'Bond':'Bonds','Angle':'Angles','Plane':'Planes',
            'Chiral':'Volumes','Torsion':'Torsions','Rama':'Ramas',
            'ChemComp':'Sites','Texture':'HKLs'}
        lasthkl = np.array([0,0,0])
        for ip,pName in enumerate(pNames):
            pnames = pName.split(':')
            if pId == int(pnames[0]):
                name = pnames[1]
                if not name:        #empty for Pawley restraints; pName has '::' in it
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
    print np.sum(pDerv),pDerv.shape
    return pDerv

################################################################################
##### Function & derivative calculations
################################################################################        
                    
def GetAtomFXU(pfx,calcControls,parmDict):
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
    return Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata
    
def getFFvalues(FFtables,SQ):
    FFvals = {}
    for El in FFtables:
        FFvals[El] = G2el.ScatFac(FFtables[El],SQ)[0]
    return FFvals
    
def getBLvalues(BLtables):
    BLvals = {}
    for El in BLtables:
        BLvals[El] = BLtables[El][1][1]
    return BLvals
        
def StructureFactor(refList,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Compute structure factors for all h,k,l for phase
    input:
        refList: [ref] where each ref = h,k,l,m,d,...,[equiv h,k,l],phase[equiv] 
        G:      reciprocal metric tensor
        pfx:    phase id string
        SGData: space group info. dictionary output from SpcGroup
        calcControls:
        ParmDict:
    puts result F^2 in each ref[8] in refList
    '''        
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,calcControls,parmDict)
    FF = np.zeros(len(Tdata))
    if 'N' in parmDict[hfx+'Type']:
        FP,FPP = G2el.BlenRes(Tdata,BLtables,parmDict[hfx+'Lam'])
    else:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    maxPos = len(SGData['SGOps'])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    for refl in refList:
        fbs = np.array([0,0])
        H = refl[:3]
        SQ = 1./(2.*refl[4])**2
        SQfactor = 4.0*SQ*twopisq
        Bab = parmDict[phfx+'BabA']*np.exp(-parmDict[phfx+'BabU']*SQfactor)
        if not len(refl[-1]):                #no form factors
            if 'N' in parmDict[hfx+'Type']:
                refl[-1] = getBLvalues(BLtables)
            else:       #'X'
                refl[-1] = getFFvalues(FFtables,SQ)
        for i,El in enumerate(Tdata):
            FF[i] = refl[-1][El]           
        Uniq = refl[11]
        phi = refl[12]
        phase = twopi*(np.inner(Uniq,(dXdata.T+Xdata.T))+phi[:,np.newaxis])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(Uniq)
        biso = -SQfactor*Uisodata
        Tiso = np.where(biso<1.,np.exp(biso),1.0)
        HbH = np.array([-np.inner(h,np.inner(bij,h)) for h in Uniq])
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = Tiso*Tuij
        fa = np.array([(FF+FP-Bab)*occ*cosp*Tcorr,-FPP*occ*sinp*Tcorr])
        fas = np.sum(np.sum(fa,axis=1),axis=1)        #real
        if not SGData['SGInv']:
            fb = np.array([(FF+FP-Bab)*occ*sinp*Tcorr,FPP*occ*cosp*Tcorr])
            fbs = np.sum(np.sum(fb,axis=1),axis=1)
        fasq = fas**2
        fbsq = fbs**2        #imaginary
        refl[9] = np.sum(fasq)+np.sum(fbsq)
        refl[10] = atan2d(fbs[0],fas[0])
    return refList
    
def StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmDict):
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    phfx = pfx.split(':')[0]+hfx
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    Tdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,calcControls,parmDict)
    FF = np.zeros(len(Tdata))
    if 'N' in parmDict[hfx+'Type']:
        FP = 0.
        FPP = 0.
    else:
        FP = np.array([FFtables[El][hfx+'FP'] for El in Tdata])
        FPP = np.array([FFtables[El][hfx+'FPP'] for El in Tdata])
    maxPos = len(SGData['SGOps'])       
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    dFdvDict = {}
    dFdfr = np.zeros((len(refList),len(Mdata)))
    dFdx = np.zeros((len(refList),len(Mdata),3))
    dFdui = np.zeros((len(refList),len(Mdata)))
    dFdua = np.zeros((len(refList),len(Mdata),6))
    dFdbab = np.zeros((len(refList),2))
    for iref,refl in enumerate(refList):
        H = np.array(refl[:3])
        SQ = 1./(2.*refl[4])**2             # or (sin(theta)/lambda)**2
        SQfactor = 8.0*SQ*np.pi**2
        dBabdA = np.exp(-parmDict[phfx+'BabU']*SQfactor)
        Bab = parmDict[phfx+'BabA']*dBabdA
        for i,El in enumerate(Tdata):            
            FF[i] = refl[-1][El]           
        Uniq = refl[11]
        phi = refl[12]
        phase = twopi*(np.inner((dXdata.T+Xdata.T),Uniq)+phi[np.newaxis,:])
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
        dfadfr = np.sum(fa/occ[:,np.newaxis],axis=2)
        dfadx = np.sum(twopi*Uniq*fax[:,:,:,np.newaxis],axis=2)
        dfadui = np.sum(-SQfactor*fa,axis=2)
        dfadua = np.sum(-Hij*fa[:,:,:,np.newaxis],axis=2)
        dfadba = np.sum(-cosp*(occ*Tcorr)[:,np.newaxis],axis=1)
        #NB: the above have been checked against PA(1:10,1:2) in strfctr.for      
        dFdfr[iref] = 2.*(fas[0]*dfadfr[0]+fas[1]*dfadfr[1])*Mdata/len(Uniq)
        dFdx[iref] = 2.*(fas[0]*dfadx[0]+fas[1]*dfadx[1])
        dFdui[iref] = 2.*(fas[0]*dfadui[0]+fas[1]*dfadui[1])
        dFdua[iref] = 2.*(fas[0]*dfadua[0]+fas[1]*dfadua[1])
        dFdbab[iref] = np.array([np.sum(dfadba*dBabdA),np.sum(-dfadba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
        if not SGData['SGInv']:
            dfbdfr = np.sum(fb/occ[:,np.newaxis],axis=2)        #problem here if occ=0 for some atom
            dfbdx = np.sum(twopi*Uniq*fbx[:,:,:,np.newaxis],axis=2)          
            dfbdui = np.sum(-SQfactor*fb,axis=2)
            dfbdua = np.sum(-Hij*fb[:,:,:,np.newaxis],axis=2)
            dfbdba = np.sum(-sinp*(occ*Tcorr)[:,np.newaxis],axis=1)
            dFdfr[iref] += 2.*(fbs[0]*dfbdfr[0]-fbs[1]*dfbdfr[1])*Mdata/len(Uniq)
            dFdx[iref] += 2.*(fbs[0]*dfbdx[0]+fbs[1]*dfbdx[1])
            dFdui[iref] += 2.*(fbs[0]*dfbdui[0]-fbs[1]*dfbdui[1])
            dFdua[iref] += 2.*(fbs[0]*dfbdua[0]+fbs[1]*dfbdua[1])
            dFdbab[iref] += np.array([np.sum(dfbdba*dBabdA),np.sum(-dfbdba*parmDict[phfx+'BabA']*SQfactor*dBabdA)]).T
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
    ''' Single crystal extinction function; puts correction in ref[13] and returns
    corrections needed for derivatives
    '''
    ref[13] = 1.0
    dervCor = 1.0
    dervDict = {}
    if calcControls[phfx+'EType'] != 'None':
        cos2T = 1.0-0.5*(parmDict[hfx+'Lam']/ref[4])**2         #cos(2theta)
        if 'SXC' in parmDict[hfx+'Type']:
            AV = 7.9406e5/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            P12 = (calcControls[phfx+'Cos2TM']+cos2T**4)/(calcControls[phfx+'Cos2TM']+cos2T**2)
        elif 'SNT' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = 1./(4.*refl[4]**2)
            P12 = 1.0
        elif 'SNC' in parmDict[hfx+'Type']:
            AV = 1.e7/parmDict[pfx+'Vol']**2
            PL = np.sqrt(1.0-cos2T**2)/parmDict[hfx+'Lam']
            P12 = 1.0
            
        PLZ = AV*P12*parmDict[hfx+'Lam']**2*ref[7]
        if 'Primary' in calcControls[phfx+'EType']:
            PLZ *= 1.5
        else:
            PLZ *= calcControls[phfx+'Tbar']
                        
        if 'Primary' in calcControls[phfx+'EType']:
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

        dervCor = (1.+PF)*PF3
        if 'Primary' in calcControls[phfx+'EType'] and phfx+'Ep' in varyList:
            dervDict[phfx+'Ep'] = -ref[7]*PLZ*PF3
        if 'II' in calcControls[phfx+'EType'] and phfx+'Es' in varyList:
            dervDict[phfx+'Es'] = -ref[7]*PLZ*PF3*(PSIG/parmDict[phfx+'Es'])**3
        if 'I' in calcControls[phfx+'EType'] and phfx+'Eg' in varyList:
            dervDict[phfx+'Eg'] = -ref[7]*PLZ*PF3*(PSIG/parmDict[phfx+'Eg'])**3*PL**2
               
        ref[13] = 1./extCor
    return dervCor,dervDict
        
    
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def GetNewCellParms(parmDict,varyList):
    newCellDict = {}
    Anames = ['A'+str(i) for i in range(6)]
    Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],Anames))
    for item in varyList:
        keys = item.split(':')
        if keys[2] in Ddict:
            key = keys[0]+'::'+Ddict[keys[2]]       #key is e.g. '0::A0'
            parm = keys[0]+'::'+keys[2]             #parm is e.g. '0::D11'
            newCellDict[parm] = [key,parmDict[key]+parmDict[item]]
    return newCellDict          # is e.g. {'0::D11':A0+D11}
    
def ApplyXYZshifts(parmDict,varyList):
    ''' takes atom x,y,z shift and applies it to corresponding atom x,y,z value
        input:
            parmDict - parameter dictionary
            varyList - list of variables
        returns:
            newAtomDict - dictionary of new atomic coordinate names & values; 
                key is parameter shift name
    '''
    newAtomDict = {}
    for item in parmDict:
        if 'dA' in item:
            parm = ''.join(item.split('d'))
            parmDict[parm] += parmDict[item]
            newAtomDict[item] = [parm,parmDict[parm]]
    return newAtomDict
    
def SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    #Spherical harmonics texture
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Omega'],parmDict[hfx+'Chi'],parmDict[hfx+'Phi'],parmDict[hfx+'Azimuth']]
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangls,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],parmDict[pfx+'SHmodel'],parmDict[pfx+'SHorder'])
    for item in SHnames:
        L,M,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,x,x = G2lat.GetKsl(L,M,parmDict[pfx+'SHmodel'],psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[pfx+item]*Lnorm*Kcl*Ksl
    return odfCor
    
def SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    #Spherical harmonics texture derivatives
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
    psi,gam,dPSdA,dGMdA = G2lat.SamAng(refl[5]/2.,Gangls,Sangls,IFCoup)
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
    #sphericaql harmonics preferred orientation (cylindrical symmetry only)
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
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta) 
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
    return odfCor
    
def SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict):
    #spherical harmonics preferred orientation derivatives (cylindrical symmetry only)
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
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta) 
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
        dFdODF[phfx+item] = Kcsl*Lnorm
    return odfCor,dFdODF
    
def GetPrefOri(refl,G,g,phfx,hfx,SGData,calcControls,parmDict):
    POcorr = 1.0
    if calcControls[phfx+'poType'] == 'MD':
        MD = parmDict[phfx+'MD']
        if MD != 1.0:
            MDAxis = calcControls[phfx+'MDAxis']
            sumMD = 0
            for H in refl[11]:            
                cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
                A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
                sumMD += A**3
            POcorr = sumMD/len(refl[11])
    else:   #spherical harmonics
        if calcControls[phfx+'SHord']:
            POcorr = SHPOcal(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr
    
def GetPrefOriDerv(refl,G,g,phfx,hfx,SGData,calcControls,parmDict):
    POcorr = 1.0
    POderv = {}
    if calcControls[phfx+'poType'] == 'MD':
        MD = parmDict[phfx+'MD']
        MDAxis = calcControls[phfx+'MDAxis']
        sumMD = 0
        sumdMD = 0
        for H in refl[11]:            
            cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
            A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
            sumMD += A**3
            sumdMD -= (1.5*A**5)*(2.0*MD*cosP**2-(sinP/MD)**2)
        POcorr = sumMD/len(refl[11])
        POderv[phfx+'MD'] = sumdMD/len(refl[11])
    else:   #spherical harmonics
        if calcControls[phfx+'SHord']:
            POcorr,POderv = SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr,POderv
    
def GetAbsorb(refl,hfx,calcControls,parmDict):
    if 'Debye' in calcControls[hfx+'instType']:
        return G2pwd.Absorb('Cylinder',parmDict[hfx+'Absorption'],refl[5],0,0)
    else:
        return 1.0
    
def GetAbsorbDerv(refl,hfx,calcControls,parmDict):
    if 'Debye' in calcControls[hfx+'instType']:
        return G2pwd.AbsorbDerv('Cylinder',parmDict[hfx+'Absorption'],refl[5],0,0)
    else:
        return 0.0
    
def GetIntensityCorr(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3]               #scale*multiplicity
    if 'X' in parmDict[hfx+'Type']:
        Icorr *= G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])[0]
    Icorr *= GetPrefOri(refl,G,g,phfx,hfx,SGData,calcControls,parmDict)
    if pfx+'SHorder' in parmDict:
        Icorr *= SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict)
    Icorr *= GetAbsorb(refl,hfx,calcControls,parmDict)
    refl[13] = Icorr        
    
def GetIntensityDerv(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    dIdsh = 1./parmDict[hfx+'Scale']
    dIdsp = 1./parmDict[phfx+'Scale']
    if 'X' in parmDict[hfx+'Type']:
        pola,dIdPola = G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])
        dIdPola /= pola
    else:       #'N'
        dIdPola = 0.0
    POcorr,dIdPO = GetPrefOriDerv(refl,G,g,phfx,hfx,SGData,calcControls,parmDict)
    for iPO in dIdPO:
        dIdPO[iPO] /= POcorr
    dFdODF = {}
    dFdSA = [0,0,0]
    if pfx+'SHorder' in parmDict:
        odfCor,dFdODF,dFdSA = SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict)
        for iSH in dFdODF:
            dFdODF[iSH] /= odfCor
        for i in range(3):
            dFdSA[i] /= odfCor
    dFdAb = GetAbsorbDerv(refl,hfx,calcControls,parmDict)
    return dIdsh,dIdsp,dIdPola,dIdPO,dFdODF,dFdSA,dFdAb
        
def GetSampleSigGam(refl,wave,G,GB,phfx,calcControls,parmDict):
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
        pwrs = calcControls[phfx+'MuPwrs']
        sum = 0
        for i,pwr in enumerate(pwrs):
            sum += parmDict[phfx+'Mustrain:'+str(i)]*refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
        Mgam = 0.018*refl[4]**2*tand(refl[5]/2.)*sum
    gam = Sgam*parmDict[phfx+'Size;mx']+Mgam*parmDict[phfx+'Mustrain;mx']
    sig = (Sgam*(1.-parmDict[phfx+'Size;mx']))**2+(Mgam*(1.-parmDict[phfx+'Mustrain;mx']))**2
    sig /= ateln2
    return sig,gam
        
def GetSampleSigGamDerv(refl,wave,G,GB,phfx,calcControls,parmDict):
    gamDict = {}
    sigDict = {}
    costh = cosd(refl[5]/2.)
    tanth = tand(refl[5]/2.)
    #crystallite size derivatives
    if calcControls[phfx+'SizeType'] == 'isotropic':
        Sgam = 1.8*wave/(np.pi*parmDict[phfx+'Size;i']*costh)
        gamDict[phfx+'Size;i'] = -900.*wave*parmDict[phfx+'Size;mx']/(np.pi*costh)
        sigDict[phfx+'Size;i'] = -1800.*Sgam*wave*(1.-parmDict[phfx+'Size;mx'])**2/(np.pi*costh*ateln2)
    elif calcControls[phfx+'SizeType'] == 'uniaxial':
        H = np.array(refl[:3])
        P = np.array(calcControls[phfx+'SizeAxis'])
        cosP,sinP = G2lat.CosSinAngle(H,P,G)
        Si = parmDict[phfx+'Size;i']
        Sa = parmDict[phfx+'Size;a']
        gami = (1.8*wave/np.pi)/(Si*Sa)
        sqtrm = np.sqrt((sinP*Sa)**2+(cosP*Si)**2)
        Sgam = gami*sqtrm
        gam = Sgam/costh
        dsi = (gami*Si*cosP**2/(sqtrm*costh)-gam/Si)
        dsa = (gami*Sa*sinP**2/(sqtrm*costh)-gam/Sa)
        gamDict[phfx+'Size;i'] = dsi*parmDict[phfx+'Size;mx']
        gamDict[phfx+'Size;a'] = dsa*parmDict[phfx+'Size;mx']
        sigDict[phfx+'Size;i'] = 2.*dsi*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
        sigDict[phfx+'Size;a'] = 2.*dsa*Sgam*(1.-parmDict[phfx+'Size;mx'])**2/ateln2
    else:           #ellipsoidal crystallites
        const = 1.8*wave/(np.pi*costh)
        Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
        H = np.array(refl[:3])
        lenR,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
        Sgam = 1.8*wave/(np.pi*costh*lenR)
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
        pwrs = calcControls[phfx+'MuPwrs']
        const = 0.018*refl[4]**2*tanth
        sum = 0
        for i,pwr in enumerate(pwrs):
            term = refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
            sum += parmDict[phfx+'Mustrain:'+str(i)]*term
            gamDict[phfx+'Mustrain:'+str(i)] = const*term*parmDict[phfx+'Mustrain;mx']
            sigDict[phfx+'Mustrain:'+str(i)] = \
                2.*const*term*(1.-parmDict[phfx+'Mustrain;mx'])**2/ateln2
        Mgam = 0.018*refl[4]**2*tand(refl[5]/2.)*sum
        for i in range(len(pwrs)):
            sigDict[phfx+'Mustrain:'+str(i)] *= Mgam
    gamDict[phfx+'Mustrain;mx'] = Mgam
    sigDict[phfx+'Mustrain;mx'] = -2.*Mgam**2*(1.-parmDict[phfx+'Mustrain;mx'])/ateln2
    return sigDict,gamDict
        
def GetReflPos(refl,wave,G,hfx,calcControls,parmDict):
    h,k,l = refl[:3]
    dsq = 1./G2lat.calc_rDsq2(np.array([h,k,l]),G)
    d = np.sqrt(dsq)

    refl[4] = d
    pos = 2.0*asind(wave/(2.0*d))+parmDict[hfx+'Zero']
    const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
    if 'Bragg' in calcControls[hfx+'instType']:
        pos -= const*(4.*parmDict[hfx+'Shift']*cosd(pos/2.0)+ \
            parmDict[hfx+'Transparency']*sind(pos)*100.0)            #trans(=1/mueff) in cm
    else:               #Debye-Scherrer - simple but maybe not right
        pos -= const*(parmDict[hfx+'DisplaceX']*cosd(pos)+parmDict[hfx+'DisplaceY']*sind(pos))
    return pos

def GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict):
    dpr = 180./np.pi
    h,k,l = refl[:3]
    dstsq = G2lat.calc_rDsq(np.array([h,k,l]),A)
    dst = np.sqrt(dstsq)
    pos = refl[5]-parmDict[hfx+'Zero']
    const = dpr/np.sqrt(1.0-wave**2*dstsq/4.0)
    dpdw = const*dst
    dpdA = np.array([h**2,k**2,l**2,h*k,h*l,k*l])
    dpdA *= const*wave/(2.0*dst)
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
            
def GetHStrainShift(refl,SGData,phfx,parmDict):
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
    return -Dij*refl[4]**2*tand(refl[5]/2.0)
            
def GetHStrainShiftDerv(refl,SGData,phfx):
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
    for item in dDijDict:
        dDijDict[item] *= -refl[4]**2*tand(refl[5]/2.0)
    return dDijDict
    
def GetFprime(controlDict,Histograms):
    FFtables = controlDict['FFtables']
    if not FFtables:
        return
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if histogram[:4] in ['PWDR','HKLF']:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            keV = controlDict[hfx+'keV']
            for El in FFtables:
                Orbs = G2el.GetXsectionCoeff(El.split('+')[0].split('-')[0])
                FP,FPP,Mu = G2el.FPcalc(Orbs, keV)
                FFtables[El][hfx+'FP'] = FP
                FFtables[El][hfx+'FPP'] = FPP                
            
def GetFobsSq(Histograms,Phases,parmDict,calcControls):
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            shl = max(parmDict[hfx+'SH/L'],0.002)
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
                refList = refLists[phase]
                sumFo = 0.0
                sumdF = 0.0
                sumFosq = 0.0
                sumdFsq = 0.0
                for refl in refList:
                    if 'C' in calcControls[hfx+'histType']:
                        yp = np.zeros_like(yb)
                        Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5],refl[6],refl[7],shl)
                        iBeg = max(xB,np.searchsorted(x,refl[5]-fmin))
                        iFin = max(xB,min(np.searchsorted(x,refl[5]+fmax),xF))
                        iFin2 = iFin
                        yp[iBeg:iFin] = refl[13]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])    #>90% of time spent here
                        if Ka2:
                            pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                            Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6],refl[7],shl)
                            iBeg2 = max(xB,np.searchsorted(x,pos2-fmin))
                            iFin2 = min(np.searchsorted(x,pos2+fmax),xF)
                            if not iBeg2+iFin2:       #peak below low limit - skip peak
                                continue
                            elif not iBeg2-iFin2:     #peak above high limit - done
                                break
                            yp[iBeg2:iFin2] += refl[13]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg2:iFin2])        #and here
                        refl[8] = np.sum(np.where(ratio[iBeg:iFin2]>0.,yp[iBeg:iFin2]*ratio[iBeg:iFin2]/(refl[13]*(1.+kRatio)),0.0))
                    elif 'T' in calcControls[hfx+'histType']:
                        print 'TOF Undefined at present'
                        raise Exception    #no TOF yet
                    Fo = np.sqrt(np.abs(refl[8]))
                    Fc = np.sqrt(np.abs(refl[9]))
                    sumFo += Fo
                    sumFosq += refl[8]**2
                    sumdF += np.abs(Fo-Fc)
                    sumdFsq += (refl[8]-refl[9])**2
                Histogram[phfx+'Rf'] = min(100.,(sumdF/sumFo)*100.)
                Histogram[phfx+'Rf^2'] = min(100.,np.sqrt(sumdFsq/sumFosq)*100.)
                Histogram[phfx+'Nref'] = len(refList)
                
def getPowderProfile(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
    def GetReflSigGam(refl,wave,G,GB,hfx,phfx,calcControls,parmDict):
        U = parmDict[hfx+'U']
        V = parmDict[hfx+'V']
        W = parmDict[hfx+'W']
        X = parmDict[hfx+'X']
        Y = parmDict[hfx+'Y']
        tanPos = tand(refl[5]/2.0)
        Ssig,Sgam = GetSampleSigGam(refl,wave,G,GB,phfx,calcControls,parmDict)
        sig = U*tanPos**2+V*tanPos+W+Ssig     #save peak sigma
        sig = max(0.001,sig)
        gam = X/cosd(refl[5]/2.0)+Y*tanPos+Sgam     #save peak gamma
        gam = max(0.001,gam)
        return sig,gam
                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    yb = G2pwd.getBackground(hfx,parmDict,bakType,x)
    yc = np.zeros_like(yb)
        
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
    else:
        print 'TOF Undefined at present'
        raise ValueError
    for phase in Histogram['Reflection Lists']:
        refList = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        hfx = ':%d:'%(hId)
        SGData = Phase['General']['SGData']
        A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        Vst = np.sqrt(nl.det(G))    #V*
        if not Phase['General'].get('doPawley'):
            refList = StructureFactor(refList,G,hfx,pfx,SGData,calcControls,parmDict)
        for refl in refList:
            if 'C' in calcControls[hfx+'histType']:
                h,k,l = refl[:3]
                refl[5] = GetReflPos(refl,wave,G,hfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = 1./(2.*sind(refl[5]/2.)**2*cosd(refl[5]/2.))           #Lorentz correction
                refl[5] += GetHStrainShift(refl,SGData,phfx,parmDict)               #apply hydrostatic strain shift
                refl[6:8] = GetReflSigGam(refl,wave,G,GB,hfx,phfx,calcControls,parmDict)    #peak sig & gam
                GetIntensityCorr(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)    #puts corrections in refl[13]
                refl[13] *= Vst*Lorenz
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
                yc[iBeg:iFin] += refl[13]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])    #>90% of time spent here
                if Ka2:
                    pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6],refl[7],shl)
                    iBeg = np.searchsorted(x,pos2-fmin)
                    iFin = np.searchsorted(x,pos2+fmax)
                    if not iBeg+iFin:       #peak below low limit - skip peak
                        continue
                    elif not iBeg-iFin:     #peak above high limit - done
                        return yc,yb
                    yc[iBeg:iFin] += refl[13]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg:iFin])        #and here
            elif 'T' in calcControls[hfx+'histType']:
                print 'TOF Undefined at present'
                raise Exception    #no TOF yet
    return yc,yb
    
def getPowderProfileDerv(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
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
    dMdb,dMddb,dMdpk = G2pwd.getBackgroundDerv(hfx,parmDict,bakType,x)
    if hfx+'Back:0' in varylist: # for now assume that Back:x vars to not appear in constraints
        bBpos =varylist.index(hfx+'Back:0')
        dMdv[bBpos:bBpos+len(dMdb)] = dMdb
    names = [hfx+'DebyeA',hfx+'DebyeR',hfx+'DebyeU']
    for name in varylist:
        if 'Debye' in name:
            id = int(name.split(':')[-1])
            parm = name[:int(name.rindex(':'))]
            ip = names.index(parm)
            dMdv[varylist.index(name)] = dMddb[3*id+ip]
    names = [hfx+'BkPkpos',hfx+'BkPkint',hfx+'BkPksig',hfx+'BkPkgam']
    for name in varylist:
        if 'BkPk' in name:
            id = int(name.split(':')[-1])
            parm = name[:int(name.rindex(':'))]
            ip = names.index(parm)
            dMdv[varylist.index(name)] = dMdpk[4*id+ip]
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
    else:
        print 'TOF Undefined at present'
        raise ValueError
    for phase in Histogram['Reflection Lists']:
        refList = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        SGData = Phase['General']['SGData']
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        if not Phase['General'].get('doPawley'):
            dFdvDict = StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmDict)
        for iref,refl in enumerate(refList):
            if 'C' in calcControls[hfx+'histType']:        #CW powder
                h,k,l = refl[:3]
                dIdsh,dIdsp,dIdpola,dIdPO,dFdODF,dFdSA,dFdAb = GetIntensityDerv(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                pos = refl[5]
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                lenBF = iFin-iBeg
                dMdpk = np.zeros(shape=(6,lenBF))
                dMdipk = G2pwd.getdFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])
                for i in range(5):
                    dMdpk[i] += 100.*cw[iBeg:iFin]*refl[13]*refl[9]*dMdipk[i]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':np.zeros_like(dMdpk[0])}
                if Ka2:
                    pos2 = refl[5]+lamRatio*tanth       # + 360/pi * Dlam/lam * tan(th)
                    iBeg2 = np.searchsorted(x,pos2-fmin)
                    iFin2 = np.searchsorted(x,pos2+fmax)
                    if iBeg2-iFin2:
                        lenBF2 = iFin2-iBeg2
                        dMdpk2 = np.zeros(shape=(6,lenBF2))
                        dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg2:iFin2])
                        for i in range(5):
                            dMdpk2[i] = 100.*cw[iBeg2:iFin2]*refl[13]*refl[9]*kRatio*dMdipk2[i]
                        dMdpk2[5] = 100.*cw[iBeg2:iFin2]*refl[13]*dMdipk2[0]
                        dervDict2 = {'int':dMdpk2[0],'pos':dMdpk2[1],'sig':dMdpk2[2],'gam':dMdpk2[3],'shl':dMdpk2[4],'L1/L2':dMdpk2[5]*refl[9]}
                if Phase['General'].get('doPawley'):
                    dMdpw = np.zeros(len(x))
                    try:
                        pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
                        idx = varylist.index(pIdx)
                        dMdpw[iBeg:iFin] = dervDict['int']/refl[9]
                        if Ka2:
                            dMdpw[iBeg2:iFin2] += dervDict2['int']/refl[9]
                        dMdv[idx] = dMdpw
                    except: # ValueError:
                        pass
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY = GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],
                    hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'SH/L':[1.0,'shl'],
                    hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
                    hfx+'DisplaceY':[dpdY,'pos'],hfx+'Absorption':[dFdAb,'int'],}
                for name in names:
                    item = names[name]
                    if name in varylist:
                        dMdv[varylist.index(name)][iBeg:iFin] += item[0]*dervDict[item[1]]
                        if Ka2:
                            dMdv[varylist.index(name)][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
                    elif name in dependentVars:
                        if Ka2:
                            depDerivDict[name][iBeg2:iFin2] += item[0]*dervDict2[item[1]]
                        depDerivDict[name][iBeg:iFin] += item[0]*dervDict[item[1]]
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
                dDijDict = GetHStrainShiftDerv(refl,SGData,phfx)
                for name in dDijDict:
                    if name in varylist:
                        dMdv[varylist.index(name)][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                        if Ka2:
                            dMdv[varylist.index(name)][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
                    elif name in dependentVars:
                        depDerivDict[name][iBeg:iFin] += dDijDict[name]*dervDict['pos']
                        if Ka2:
                            depDerivDict[name][iBeg2:iFin2] += dDijDict[name]*dervDict2['pos']
                sigDict,gamDict = GetSampleSigGamDerv(refl,wave,G,GB,phfx,calcControls,parmDict)
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
                    if phfx+name in varylist:
                        dMdv[varylist.index(phfx+name)][iBeg:iFin] += dFdvDict[pfx+name][iref]*dervDict['int']*cw[iBeg:iFin]
                        if Ka2:
                            dMdv[varylist.index(phfx+name)][iBeg2:iFin2] += dFdvDict[pfx+name][iref]*dervDict2['int']*cw[iBeg2:iFin2]
                    elif phfx+name in dependentVars:                    
                        depDerivDict[phfx+name][iBeg:iFin] += dFdvDict[pfx+name][iref]*dervDict['int']*cw[iBeg:iFin]
                        if Ka2:
                            depDerivDict[phfx+name][iBeg2:iFin2] += dFdvDict[pfx+name][iref]*dervDict2['int']*cw[iBeg2:iFin2]                  
            elif 'T' in calcControls[hfx+'histType']:
                print 'TOF Undefined at present'
                raise Exception    #no TOF yet
            #do atom derivatives -  for F,X & U so far              
            corr = dervDict['int']/refl[9]
            if Ka2:
                corr2 = dervDict2['int']/refl[9]
            for name in varylist+dependentVars:
                try:
                    aname = name.split(pfx)[1][:2]
                    if aname not in ['Af','dA','AU']: continue # skip anything not an atom param
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
    # now process derivatives in constraints
    G2mv.Dict2Deriv(varylist,depDerivDict,dMdv)
    return dMdv

def dervRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):
    parmdict.update(zip(varylist,values))
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases,restraintDict = HistoPhases
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
            W = wtFactor*w
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            dMdvh = np.sqrt(W[xB:xF])*getPowderProfileDerv(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            nobs = Histogram['Nobs']
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            pfx = '%d::'%(Phase['pId'])
            phfx = '%d:%d:'%(Phase['pId'],hId)
            SGData = Phase['General']['SGData']
            A = [parmdict[pfx+'A%d'%(i)] for i in range(6)]
            G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
            refList = Histogram['Data']
            dFdvDict = StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmdict)
            dMdvh = np.zeros((len(varylist),len(refList)))
            for iref,ref in enumerate(refList):
                if ref[6] > 0:
                    dervCor,dervDict = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmdict,varylist) #puts correction in refl[13]
                    if calcControls['F**2']:
                        if ref[5]/ref[6] >= calcControls['minF/sig']:
                            w = wtFactor/ref[6]
                            for j,var in enumerate(varylist):
                                if var in dFdvDict:
                                    dMdvh[j][iref] = w*dFdvDict[var][iref]*dervCor
                            if phfx+'Scale' in varylist:
                                dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*dervCor
                    else:
                        Fo = np.sqrt(ref[5])
                        Fc = np.sqrt(ref[7])
                        sig = ref[6]/(2.0*Fo)
                        if Fo/sig >= calcControls['minF/sig']:
                            w = wtFactor/sig
                            for j,var in enumerate(varylist):
                                if var in dFdvDict:
                                    dMdvh[j][iref] = w*dFdvDict[var][iref]*np.sqrt(dervCor)
                            if phfx+'Scale' in varylist:
                                dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*np.sqrt(dervCor)                            
                    for item in ['Ep','Es','Eg']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]
        else:
            continue        #skip non-histogram entries
        if len(dMdv):
            dMdv = np.concatenate((dMdv.T,dMdvh.T)).T
        else:
            dMdv = dMdvh
            
    pNames,pVals,pWt = penaltyFxn(HistoPhases,parmdict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,parmdict,varylist)
        dMdv = np.concatenate((dMdv.T,dpdv.T)).T
        
    return dMdv

def HessRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):
    parmdict.update(zip(varylist,values))
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases,restraintDict = HistoPhases
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
            dMdvh = getPowderProfileDerv(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            Wt = np.sqrt(W[xB:xF])[np.newaxis,:]
            Dy = dy[xB:xF][np.newaxis,:]
            dMdvh *= Wt
            if dlg:
                dlg.Update(Histogram['wR'],newmsg='Hessian for histogram %d\nAll data Rw=%8.3f%s'%(hId,Histogram['wR'],'%'))[0]
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
            nobs = Histogram['Nobs']
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            pfx = '%d::'%(Phase['pId'])
            phfx = '%d:%d:'%(Phase['pId'],hId)
            SGData = Phase['General']['SGData']
            A = [parmdict[pfx+'A%d'%(i)] for i in range(6)]
            G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
            refList = Histogram['Data']
            dFdvDict = StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmdict)
            dMdvh = np.zeros((len(varylist),len(refList)))
            wdf = np.zeros(len(refList))
            for iref,ref in enumerate(refList):
                if ref[6] > 0:
                    dervCor,dervDict = SCExtinction(ref,phfx,hfx,pfx,calcControls,parmdict,varylist) #puts correction in refl[13]
                    if calcControls['F**2']:
                        if ref[5]/ref[6] >= calcControls['minF/sig']:
                            w =  wtFactor/ref[6]
                            wdf[iref] = w*(ref[5]-ref[7])
                            for j,var in enumerate(varylist):
                                if var in dFdvDict:
                                    dMdvh[j][iref] = w*dFdvDict[var][iref]*dervCor
                            if phfx+'Scale' in varylist:
                                dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*dervCor
                    else:
                        if ref[5] > 0.:
                            Fo = np.sqrt(ref[5])
                            Fc = np.sqrt(ref[7])
                            sig = ref[6]/(2.0*Fo)
                            w = wtFactor/sig
                            wdf[iref] = w*(Fo-Fc)
                            if Fo/sig >= calcControls['minF/sig']:
                                for j,var in enumerate(varylist):
                                    if var in dFdvDict:
                                        dMdvh[j][iref] = w*dFdvDict[var][iref]*np.sqrt(dervCor)
                                if phfx+'Scale' in varylist:
                                    dMdvh[varylist.index(phfx+'Scale')][iref] = w*ref[9]*np.sqrt(dervCor)                           
                    for item in ['Ep','Es','Eg']:
                        if phfx+item in varylist:
                            dMdvh[varylist.index(phfx+item)][iref] = w*dervDict[phfx+item]
            if dlg:
                dlg.Update(Histogram['wR'],newmsg='Hessian for histogram %d Rw=%8.3f%s'%(hId,Histogram['wR'],'%'))[0]
            if len(Hess):
                Vec += np.sum(dMdvh*wdf,axis=1)
                Hess += np.inner(dMdvh,dMdvh)
            else:
                Vec = np.sum(dMdvh*wdf,axis=1)
                Hess = np.inner(dMdvh,dMdvh)
        else:
            continue        #skip non-histogram entries
    pNames,pVals,pWt = penaltyFxn(HistoPhases,parmdict,varylist)
    if np.any(pVals):
        dpdv = penaltyDeriv(pNames,pVals,HistoPhases,parmdict,varylist)
        Vec += np.sum(dpdv*pWt*pVals,axis=1)
        Hess += np.inner(dpdv*pWt,dpdv)
    return Vec,Hess

def errRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):        
    parmdict.update(zip(varylist,values))
    Values2Dict(parmdict, varylist, values)
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases,restraintDict = HistoPhases
    M = np.empty(0)
    SumwYo = 0
    Nobs = 0
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
            yc *= 0.0                           #zero full calcd profiles
            yb *= 0.0
            yd *= 0.0
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            Histogram['Nobs'] = xF-xB
            Nobs += Histogram['Nobs']
            Histogram['sumwYo'] = np.sum(W[xB:xF]*y[xB:xF]**2)
            SumwYo += Histogram['sumwYo']
            yc[xB:xF],yb[xB:xF] = getPowderProfile(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            yc[xB:xF] += yb[xB:xF]
            yd[xB:xF] = y[xB:xF]-yc[xB:xF]
            Histogram['sumwYd'] = np.sum(np.sqrt(W[xB:xF])*(yd[xB:xF]))
            wdy = -np.sqrt(W[xB:xF])*(yd[xB:xF])
            Histogram['wR'] = min(100.,np.sqrt(np.sum(wdy**2)/Histogram['sumwYo'])*100.)
            if dlg:
                dlg.Update(Histogram['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['wR'],'%'))[0]
            M = np.concatenate((M,wdy))
#end of PWDR processing
        elif 'HKLF' in histogram[:4]:
            Histogram = Histograms[histogram]
            phase = Histogram['Reflection Lists']
            Phase = Phases[phase]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            wtFactor = calcControls[hfx+'wtFactor']
            pfx = '%d::'%(Phase['pId'])
            phfx = '%d:%d:'%(Phase['pId'],hId)
            SGData = Phase['General']['SGData']
            A = [parmdict[pfx+'A%d'%(i)] for i in range(6)]
            G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
            refList = Histogram['Data']
            refList = StructureFactor(refList,G,hfx,pfx,SGData,calcControls,parmdict)
            df = np.zeros(len(refList))
            sumwYo = 0
            sumFo = 0
            sumFo2 = 0
            sumdF = 0
            sumdF2 = 0
            nobs = 0
            for i,ref in enumerate(refList):
                if ref[6] > 0:
                    SCExtinction(ref,phfx,hfx,pfx,calcControls,parmdict,varylist) #puts correction in refl[13]
                    ref[7] = parmdict[phfx+'Scale']*ref[9]
                    ref[7] *= ref[13]
                    ref[8] = ref[5]/parmdict[phfx+'Scale']
                    if calcControls['F**2']:
                        if ref[5]/ref[6] >= calcControls['minF/sig']:
                            sumFo2 += ref[5]
                            Fo = np.sqrt(ref[5])
                            sumFo += Fo
                            sumFo2 += ref[5]
                            sumdF += abs(Fo-np.sqrt(ref[7]))
                            sumdF2 += abs(ref[5]-ref[7])
                            nobs += 1
                            df[i] = -np.sqrt(wtFactor)*(ref[5]-ref[7])/ref[6]
                            sumwYo += wtFactor*(ref[5]/ref[6])**2
                    else:
                        Fo = np.sqrt(ref[5])
                        Fc = np.sqrt(ref[7])
                        sig = ref[6]/(2.0*Fo)
                        if Fo/sig >= calcControls['minF/sig']:
                            sumFo += Fo
                            sumFo2 += ref[5]
                            sumdF += abs(Fo-Fc)
                            sumdF2 += abs(ref[5]-ref[7])
                            nobs += 1
                            df[i] = -np.sqrt(wtFactor)*(Fo-Fc)/sig
                            sumwYo += wtFactor*(Fo/sig)**2
            Histogram['Nobs'] = nobs
            Histogram['sumwYo'] = sumwYo
            SumwYo += sumwYo
            Histogram['wR'] = min(100.,np.sqrt(np.sum(df**2)/Histogram['sumwYo'])*100.)
            Histogram[phfx+'Rf'] = 100.*sumdF/sumFo
            Histogram[phfx+'Rf^2'] = 100.*sumdF2/sumFo2
            Histogram[phfx+'Nref'] = nobs
            Nobs += nobs
            if dlg:
                dlg.Update(Histogram['wR'],newmsg='For histogram %d Rw=%8.3f%s'%(hId,Histogram['wR'],'%'))[0]
            M = np.concatenate((M,df))
# end of HKLF processing
    Histograms['sumwYo'] = SumwYo
    Histograms['Nobs'] = Nobs
    Rw = min(100.,np.sqrt(np.sum(M**2)/SumwYo)*100.)
    if dlg:
        GoOn = dlg.Update(Rw,newmsg='%s%8.3f%s'%('All data Rw =',Rw,'%'))[0]
        if not GoOn:
            parmdict['saved values'] = values
            dlg.Destroy()
            raise Exception         #Abort!!
    pDict,pVals,pWt = penaltyFxn(HistoPhases,parmdict,varylist)
    if np.any(pVals):
        print 'Penalty function: %.3f on %d terms'%(np.sum(pWt*pVals**2),len(pVals))
        Nobs += len(pVals)
        M = np.concatenate((M,pWt*pVals))
    return M
                        
def Refine(GPXfile,dlg):
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    ShowBanner(printFile)
    varyList = []
    parmDict = {}
    G2mv.InitVars()    
    Controls = GetControls(GPXfile)
    ShowControls(Controls,printFile)
    calcControls = {}
    calcControls.update(Controls)            
    constrDict,fixedList = GetConstraints(GPXfile)
    restraintDict = GetRestraints(GPXfile)
    rigidbodyDict = GetRigidBodies(GPXfile)
    Histograms,Phases = GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no phases! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception        
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = GetPhaseData(Phases,restraintDict,pFile=printFile)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histograms,pFile=printFile)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = GetHistogramData(Histograms,pFile=printFile)
    calcControls.update(controlDict)
    varyList = phaseVary+hapVary+histVary
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    GetFprime(calcControls,Histograms)
    # do constraint processing
    try:
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList)
    except:
        print ' *** ERROR - your constraints are internally inconsistent ***'
        # traceback for debug
        #print 'varyList',varyList
        #print 'constrDict',constrDict
        #print 'fixedList',fixedList
        #import traceback
        #print traceback.format_exc()
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
        values =  np.array(Dict2Values(parmDict, varyList))
        Ftol = Controls['min dM/M']
        Factor = Controls['shift factor']
        maxCyc = Controls['max cyc']
        if 'Jacobian' in Controls['deriv type']:            
            result = so.leastsq(errRefine,values,Dfun=dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
        elif 'Hessian' in Controls['deriv type']:
            result = G2mth.HessianLSQ(errRefine,values,Hess=HessRefine,ftol=Ftol,maxcyc=maxCyc,
                args=([Histograms,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
        else:           #'numeric'
            result = so.leastsq(errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/len(varyList))
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
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
    GetFobsSq(Histograms,Phases,parmDict,calcControls)
#    print 'test2'
    sigDict = dict(zip(varyList,sig))
    newCellDict = GetNewCellParms(parmDict,varyList)
    newAtomDict = ApplyXYZshifts(parmDict,varyList)
    covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
        'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
    # add the uncertainties into the esd dictionary (sigDict)
    sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
    G2mv.PrintIndependentVars(parmDict,varyList,sigDict,pFile=printFile)
    SetPhaseData(parmDict,sigDict,Phases,covData,restraintDict,printFile)
    SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,pFile=printFile)
    SetHistogramData(parmDict,sigDict,Histograms,pFile=printFile)
    SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,covData)
    printFile.close()
    print ' Refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst'
    print ' ***** Refinement successful *****'
    
#for testing purposes!!!
    if DEBUG:
        import cPickle
        fl = open('structTestdata.dat','wb')
        cPickle.dump(parmDict,fl,1)
        cPickle.dump(varyList,fl,1)
        for histogram in Histograms:
            if 'PWDR' in histogram[:4]:
                Histogram = Histograms[histogram]
        cPickle.dump(Histogram,fl,1)
        cPickle.dump(Phases,fl,1)
        cPickle.dump(calcControls,fl,1)
        cPickle.dump(pawleyLookup,fl,1)
        fl.close()

    if dlg:
        return Rvals['Rwp']

def SeqRefine(GPXfile,dlg):
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    print ' Sequential Refinement'
    ShowBanner(printFile)
    G2mv.InitVars()    
    Controls = GetControls(GPXfile)
    ShowControls(Controls,printFile)            
    constrDict,fixedList = GetConstraints(GPXfile)
    restraintDict = GetRestraints(GPXfile)
    Histograms,Phases = GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no histograms to refine! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = GetPhaseData(Phases,restraintDict,False,printFile)
    for item in phaseVary:
        if '::A0' in item:
            print '**** WARNING - lattice parameters should not be refined in a sequential refinement ****'
            print '****           instead use the Dij parameters for each powder histogram            ****'
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else:
        histNames = GetHistogramNames(GPXfile,['PWDR',])
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
        hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histo,False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = phaseVary+hapVary+histVary
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
        GetFprime(calcControls,Histo)
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
            values =  np.array(Dict2Values(parmDict,varyList))
            Ftol = Controls['min dM/M']
            Factor = Controls['shift factor']
            maxCyc = Controls['max cyc']

            if 'Jacobian' in Controls['deriv type']:            
                result = so.leastsq(errRefine,values,Dfun=dervRefine,full_output=True,
                    ftol=Ftol,col_deriv=True,factor=Factor,
                    args=([Histo,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/2)
            elif 'Hessian' in Controls['deriv type']:
                result = G2mth.HessianLSQ(errRefine,values,Hess=HessRefine,ftol=Ftol,maxcyc=maxCyc,
                    args=([Histo,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = result[2]['num cyc']+1                           
            else:           #'numeric'
                result = so.leastsq(errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                    args=([Histo,Phases,restraintDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/len(varyList))

            runtime = time.time()-begin
            Rvals['chisq'] = np.sum(result[2]['fvec']**2)
            Values2Dict(parmDict, varyList, result[0])
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
    
        GetFobsSq(Histo,Phases,parmDict,calcControls)
        sigDict = dict(zip(varyList,sig))
        newCellDict = GetNewCellParms(parmDict,varyList)
        newAtomDict = ApplyXYZshifts(parmDict,varyList)
        covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
            'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
        # add the uncertainties into the esd dictionary (sigDict)
        SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,ifPrint,printFile)
        SetHistogramData(parmDict,sigDict,Histo,ifPrint,printFile)
        SeqResult[histogram] = covData
        SetUsedHistogramsAndPhases(GPXfile,Histo,Phases,covData,makeBack)
        makeBack = False
    SetSeqResult(GPXfile,Histograms,SeqResult)
    printFile.close()
    print ' Sequential refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst'
    print ' ***** Sequential refinement successful *****'

def DistAngle(DisAglCtls,DisAglData):
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
        cellSig = getCellEsd(pfx,SGData,A,covData)
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
