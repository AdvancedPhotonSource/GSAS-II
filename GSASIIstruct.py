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
import numpy as np
import cPickle
import time
import math
import wx
import GSASIIpath
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIpwd as G2pwd
import GSASIImapvars as G2mv
import scipy.optimize as so

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi


def ShowBanner():
    print 80*'*'
    print '   General Structure Analysis System-II Crystal Structure Refinement'
    print '     by Robert B. Von Dreele, Argonne National Laboratory(C), 2010'
    print ' This product includes software developed by the UChicago Argonne, LLC,' 
    print '            as Operator of Argonne National Laboratory.'
    print 80*'*','\n'

def GetControls(GPXfile):
    ''' Returns dictionary of control items found in GSASII gpx file
    input:
        GPXfile = .gpx full file name
    return:
        Controls = dictionary of control items
    '''
    file = open(GPXfile,'rb')
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == 'Controls':
            Controls = datum[1]
    file.close()
    return Controls
    
def ShowControls(Controls):
    print ' Least squares controls:'
    print ' Derivative type: ',Controls['deriv type']
    print ' Minimum delta-M/M for convergence: ','%.2g'%(Controls['min dM/M'])
    
def GetPhaseNames(GPXfile):
    ''' Returns a list of phase names found under 'Phases' in GSASII gpx file
    input: 
        GPXfile = gpx full file name
    return: 
        PhaseNames = list of phase names
    '''
    file = open(GPXfile,'rb')
    PhaseNames = []
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if 'Phases' == datum[0]:
            for datus in data[1:]:
                PhaseNames.append(datus[0])
    file.close()
    return PhaseNames

def GetAllPhaseData(GPXfile,PhaseName):
    ''' Returns the entire dictionary for PhaseName from GSASII gpx file
    input:
        GPXfile = gpx full file name
        PhaseName = phase name
    return:
        phase dictionary
    '''        
    file = open(GPXfile,'rb')
    General = {}
    Atoms = []
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if 'Phases' == datum[0]:
            for datus in data[1:]:
                if datus[0] == PhaseName:
                    break
    file.close()
    return datus[1]
    
def GetHistogramNames(GPXfile):
    ''' Returns a list of histogram names found in GSASII gpx file
    input: 
        GPXfile = .gpx full file name
    return: 
        HistogramNames = list of histogram names (types = PWDR & HKLF)
    '''
    file = open(GPXfile,'rb')
    HistogramNames = []
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if datum[0][:4] in ['PWDR','HKLF']:
            HistogramNames.append(datum[0])
    file.close()
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
    phaseData = {}
    for name in phaseNames: 
        phaseData[name] =  GetAllPhaseData(GPXfile,name)
    Histograms = {}
    Phases = {}
    pId = 0
    hId = 0
    for phase in phaseData:
        Phase = phaseData[phase]
        if Phase['Histograms']:
            if phase not in Phases:
                Phase['pId'] = pId
                pId += 1
                Phases[phase] = Phase
            for hist in Phase['Histograms']:
                if hist not in Histograms:
                    if 'PWDR' in hist[:4]: 
                        Histograms[hist] = GetPWDRdata(GPXfile,hist)
                    elif 'HKLF' in hist[:4]:
                        Histograms[hist] = GetHKLFdata(GPXfile,hist)
                    #future restraint, etc. histograms here            
                    Histograms[hist]['hId'] = hId
                    hId += 1
    return Histograms,Phases
    
def SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases):
    ''' Updates gpxfile from all histograms that are found in any phase
    and any phase that used a histogram
    input:
        GPXfile = .gpx full file name
        Histograms = dictionary of histograms as {name:data,...}
        Phases = dictionary of phases that use histograms
    '''
                        
    def GPXBackup(GPXfile):
        import distutils.file_util as dfu
        GPXpath,GPXname = ospath.split(GPXfile)
        Name = ospath.splitext(GPXname)[0]
        files = os.listdir(GPXpath)
        last = 0
        for name in files:
            name = name.split('.')
            if len(name) == 3 and name[0] == Name and 'bak' in name[1]:
                last = max(last,int(name[1].strip('bak'))+1)
        GPXback = ospath.join(GPXpath,ospath.splitext(GPXname)[0]+'.bak'+str(last)+'.gpx')
        dfu.copy_file(GPXfile,GPXback)
        return GPXback
        
    GPXback = GPXBackup(GPXfile)
    print '\n',130*'-'
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
        print 'read: ',datum[0]
        if datum[0] == 'Phases':
            for iphase in range(len(data)):
                if data[iphase][0] in Phases:
                    phaseName = data[iphase][0]
                    data[iphase][1] = Phases[phaseName]
        try:
            histogram = Histograms[datum[0]]
            print 'found ',datum[0]
            data[0][1][1] = histogram['Data']
            for datus in data[1:]:
                print '    read: ',datus[0]
                if datus[0] in ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']:
                    datus[1] = histogram[datus[0]]
        except KeyError:
            pass
                                
        cPickle.dump(data,outfile,1)
    infile.close()
    outfile.close()
    print 'refinement save successful'
                    
def GetPWDRdata(GPXfile,PWDRname):
    ''' Returns powder data from GSASII gpx file
    input: 
        GPXfile = .gpx full file name
        PWDRname = powder histogram name as obtained from GetHistogramNames
    return: 
        PWDRdata = powder data dictionary with:
            Data - powder data arrays, Limits, Instrument Parameters, Sample Parameters
        
    '''
    file = open(GPXfile,'rb')
    PWDRdata = {}
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == PWDRname:
            PWDRdata['Data'] = datum[1][1]          #powder data arrays
            PWDRdata[data[2][0]] = data[2][1]       #Limits
            PWDRdata[data[3][0]] = data[3][1]       #Background
            PWDRdata[data[4][0]] = data[4][1]       #Instrument parameters
            PWDRdata[data[5][0]] = data[5][1]       #Sample parameters
            try:
                PWDRdata[data[9][0]] = data[9][1]       #Reflection lists might be missing
            except IndexError:
                PWDRdata['Reflection lists'] = {}
    file.close()
    return PWDRdata
    
def GetHKLFdata(GPXfile,HKLFname):
    ''' Returns single crystal data from GSASII gpx file
    input: 
        GPXfile = .gpx full file name
        HKLFname = single crystal histogram name as obtained from GetHistogramNames
    return: 
        HKLFdata = single crystal data list of reflections: for each reflection:
            HKLF = [np.array([h,k,l]),FoSq,sigFoSq,FcSq,Fcp,Fcpp,phase]
    '''
    file = open(GPXfile,'rb')
    HKLFdata = []
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == HKLFname:
            HKLFdata = datum[1:][0]
    file.close()
    return HKLFdata
    
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
            if item['Symbol'] == El:
                FFtable[El] = item
    return FFtable
    
def GetPawleyConstr(SGLaue,PawleyRef,pawleyVary):
    if SGLaue in ['-1','2/m','mmm']:
        return                      #no Pawley symmetry required constraints
    for i,varyI in enumerate(pawleyVary):
        refI = int(varyI.split(':')[-1])
        ih,ik,il = PawleyRef[refI][:3]
        for varyJ in pawleyVary[0:i]:
            refJ = int(varyJ.split(':')[-1])
            jh,jk,jl = PawleyRef[refJ][:3]
            if SGLaue in ['4/m','4/mmm']:
                isum = ih**2+ik**2
                jsum = jh**2+jk**2
                if abs(il) == abs(jl) and isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['3R','3mR']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2*jl**2
                isum2 = ih*ik+ih*il+ik*il
                jsum2 = jh*jk+jh*jl+jk*jl
                if isum == jsum and isum2 == jsum2:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                isum = ih**2+ik**2+ih*ik
                jsum = jh**2+jk**2+jh*jk
                if abs(il) == abs(jl) and isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['m3','m3m']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2+jl**2
                if isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
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
        return [pfx+'A0']
    
        
def GetPhaseData(PhaseData):
    
        
    print ' Phases:'
    phaseVary = []
    phaseDict = {}
    phaseConstr = {}
    pawleyLookup = {}
    for name in PhaseData:
        General = PhaseData[name]['General']
        pId = PhaseData[name]['pId']
        pfx = str(pId)+'::'
        Atoms = PhaseData[name]['Atoms']
        try:
            PawleyRef = PhaseData[name]['Pawley ref']
        except KeyError:
            PawleyRef = []
        print '\n Phase name: ',General['Name']
        SGData = General['SGData']
        SGtext = G2spc.SGPrint(SGData)
        for line in SGtext: print line
        cell = General['Cell']
        A = G2lat.cell2A(cell[1:7])
        phaseDict.update({pfx+'A0':A[0],pfx+'A1':A[1],pfx+'A2':A[2],pfx+'A3':A[3],pfx+'A4':A[4],pfx+'A5':A[5]})
        if cell[0]:
            phaseVary = cellVary(pfx,SGData)
        print '\n Unit cell: a =','%.5f'%(cell[1]),' b =','%.5f'%(cell[2]),' c =','%.5f'%(cell[3]), \
            ' alpha =','%.3f'%(cell[4]),' beta =','%.3f'%(cell[5]),' gamma =', \
            '%.3f'%(cell[6]),' volume =','%.3f'%(cell[7]),' Refine?',cell[0]
        if Atoms:
            print '\n Atoms:'
            line = '   name    type  refine?   x         y         z    '+ \
                '  frac site sym  mult I/A   Uiso     U11     U22     U33     U12     U13     U23'
            if General['Type'] == 'magnetic':
                line += '   Mx     My     Mz'
            elif General['Type'] == 'macromolecular':
                line = ' res no  residue  chain '+line
            print line
            if General['Type'] == 'nuclear':
                print 135*'-'
                for at in Atoms:
                    line = '%7s'%(at[0])+'%7s'%(at[1])+'%7s'%(at[2])+'%10.5f'%(at[3])+'%10.5f'%(at[4])+ \
                        '%10.5f'%(at[5])+'%8.3f'%(at[6])+'%7s'%(at[7])+'%5d'%(at[8])+'%5s'%(at[9])
                    if at[9] == 'I':
                        line += '%8.4f'%(at[10])+48*' '
                    else:
                        line += 8*' '
                        for i in range(6):
                            line += '%8.4f'%(at[11+i])
                    print line
                    if 'X' in at[2]:
                        xId,xCoef = G2spc.GetCSxinel(at[7])
                    if 'U' in at[2]:
                        uId,uCoef = G2spc.GetCSuinel(at[7])[:2]
#        elif General['Type'] == 'magnetic':
#        elif General['Type'] == 'macromolecular':
#       PWDR: harmonics texture, unit cell, atom parameters
        elif PawleyRef:
            pawleyVary = []
            for i,refl in enumerate(PawleyRef):
                phaseDict[pfx+'PWLref:'+str(i)] = refl[6]
                pawleyLookup[pfx+'%d,%d,%d'%(refl[0],refl[1],refl[2])] = i
                if refl[5]:
                    pawleyVary.append(pfx+'PWLref:'+str(i))
            GetPawleyConstr(SGData['SGLaue'],PawleyRef,pawleyVary)      #does G2mv.StoreEquivalence
            phaseVary += pawleyVary    
                
    return phaseVary,phaseDict,pawleyLookup
    
def SetPhaseData(parmDict,sigDict,Phases):
    
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
            
    print '\n Phases:'
    for phase in Phases:
        print ' Result for phase: ',phase
        Phase = Phases[phase]
        General = Phase['General']
        SGData = General['SGData']
        cell = General['Cell']
        pId = Phase['pId']
        pfx = str(pId)+'::'
        if cell[0]:
            A,sigA = cellFill(pfx,SGData,parmDict,sigDict)
            print ' Reciprocal metric tensor: '
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
            print namstr
            print ptstr
            print sigstr
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            print ' New unit cell: a = %.5f'%(cell[1]),' b = %.5f'%(cell[2]), \
                ' c = %.5f'%(cell[3]),' alpha = %.3f'%(cell[4]),' beta = %.3f'%(cell[5]), \
                ' gamma = %.3f'%(cell[6]),' volume = %.3f'%(cell[7])
        if 'Pawley' in Phase['General']['Type']:
            pawleyRef = Phase['Pawley ref']
            for i,refl in enumerate(pawleyRef):
                key = pfx+'PWLref:'+str(i)
                refl[6] = abs(parmDict[key])        #suppress negative Fsq
                if key in sigDict:
                    refl[7] = sigDict[key]
                else:
                    refl[7] = 0

def GetHistogramPhaseData(Phases,Histograms):
    
    def PrintSize(hapData):
        line = '\n Size model    : '+hapData[0]
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:'+'%12.2f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.2f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            print line
        else:
            print line
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print ptlbls
            print ptstr
            print varstr
        
    def PrintMuStrain(hapData,SGData):
        line = '\n Mustrain model: '+hapData[0]
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:'+'%12.4f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.4f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            print line
        else:
            print line
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print ptlbls
            print ptstr
            print varstr
        
    
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
        for histogram in HistoPhase:
            print '\n Phase: ',phase,' in histogram: ',histogram
            hapData = HistoPhase[histogram]
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            limits = Histogram['Limits'][1]
            inst = Histogram['Instrument Parameters']
            inst = dict(zip(inst[3],inst[1]))
            Zero = inst['Zero']
            if 'C' in inst['Type']:
                try:
                    wave = inst['Lam']
                except KeyError:
                    wave = inst['Lam1']
                dmin = wave/(2.0*sind(limits[1]/2.0))
            pfx = str(pId)+':'+str(hId)+':'
            for item in ['Scale','Extinction']:
                hapDict[pfx+item] = hapData[item][0]
                if hapData[item][1]:
                    hapVary.append(pfx+item)
            controlDict[pfx+'poType'] = hapData['Pref.Ori.'][0]
            if hapData['Pref.Ori.'][0] == 'MD':
                hapDict[pfx+'MD'] = hapData['Pref.Ori.'][1]
                controlDict[pfx+'MDAxis'] = hapData['Pref.Ori.'][3]
                if hapData['Pref.Ori.'][2]:
                    hapVary.append(pfx+'MD')
            else:                           #'SH' spherical harmonics
                for item in hapData['Pref.Ori.'][5]:
                    hapDict[pfx+item] = hapData['Pref.Ori.'][5][item]
                    if hapData['Pref.Ori.'][2]:
                        hapVary.append(pfx+item)                    
            for item in ['Mustrain','Size']:
                controlDict[pfx+item+'Type'] = hapData[item][0]
                if hapData[item][0] in ['isotropic','uniaxial']:
                    hapDict[pfx+item+':0'] = hapData[item][1][0]
                    if hapData[item][2][0]:
                        hapVary.append(pfx+item+':0')
                    if hapData[item][0] == 'uniaxial':
                        controlDict[pfx+item+'Axis'] = hapData[item][3]
                        hapDict[pfx+item+':1'] = hapData[item][1][1]
                        if hapData[item][2][1]:
                            hapVary.append(pfx+item+':1')
                else:       #generalized for mustrain or ellipsoidal for size
                    if item == 'Mustrain':
                        names = G2spc.MustrainNames(SGData)
                        pwrs = []
                        for name in names:
                            h,k,l = name[1:]
                            pwrs.append([int(h),int(k),int(l)])
                        controlDict[pfx+'MuPwrs'] = pwrs
                    for i in range(len(hapData[item][4])):
                        sfx = ':'+str(i)
                        hapDict[pfx+item+sfx] = hapData[item][4][i]
                        if hapData[item][5][i]:
                            hapVary.append(pfx+item+sfx)
                            
            print '\n Phase fraction  : %10.4f'%(hapData['Scale'][0]),' Refine?',hapData['Scale'][1]
            print ' Extinction coeff: %10.4f'%(hapData['Extinction'][0]),' Refine?',hapData['Extinction'][1]
            if hapData['Pref.Ori.'][0] == 'MD':
                Ax = hapData['Pref.Ori.'][3]
                print ' March-Dollase PO: %10.4f'%(hapData['Pref.Ori.'][1]),' Refine?',hapData['Pref.Ori.'][2], \
                    ' Axis: %d %d %d'%(Ax[0],Ax[1],Ax[2])                
            PrintSize(hapData['Size'])
            PrintMuStrain(hapData['Mustrain'],SGData)
            HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
            refList = []
            for h,k,l,d in HKLd:
                ext,mul,Uniq = G2spc.GenHKLf([h,k,l],SGData)
                if ext:
                    continue
                if 'C' in inst['Type']:
                    pos = 2.0*asind(wave/(2.0*d))
                    if limits[0] < pos < limits[1]:
                        refList.append([h,k,l,mul,d,pos,0.0,0.0,0.0,0.0,Uniq])
                else:
                    raise ValueError 
            Histogram['Reflection Lists'][phase] = refList
    return hapVary,hapDict,controlDict
    
def SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms):
    
    def PrintSizeAndSig(hapData,sizeSig):
        line = '\n Size model: '+hapData[0]
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.2f'%(hapData[1][0])
            if sizeSig[0][0]:
                line += ', sig: %8.2f'%(sizeSig[0][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.2f'%(hapData[1][1])
                if sizeSig[0][1]:
                    line += ', sig: %8.2f'%(sizeSig[0][1])
            print line
        else:
            print line
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                if sizeSig[1][i]:
                    sigstr += '%12.6f' % (sizeSig[1][i])
                else:
                    sigstr += 12*' '
            print ptlbls
            print ptstr
            print sigstr
        
    def PrintMuStrainAndSig(hapData,mustrainSig,SGData):
        line = '\n Mustrain model: '+hapData[0]
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.4f'%(hapData[1][0])
            if mustrainSig[0][0]:
                line += ',sig: %12.4f'%(mustrainSig[0][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.4f'%(hapData[1][1])
                if mustrainSig[0][1]:
                     line += ',sig: %12.4f'%(mustrainSig[0][1])
            print line
        else:
            print line
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                if mustrainSig[1][i]:
                    sigstr += '%12.6f' % (mustrainSig[1][i])
                else:
                    sigstr += 12*' '
            print ptlbls
            print ptstr
            print sigstr
        
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        SGData = Phases[phase]['General']['SGData']
        pId = Phases[phase]['pId']
        for histogram in HistoPhase:
            print '\n Phase: ',phase,' in histogram: ',histogram
            hapData = HistoPhase[histogram]
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = str(pId)+':'+str(hId)+':'
            
            PhFrExtPOSig = {}
            for item in ['Scale','Extinction']:
                hapData[item][0] = parmDict[pfx+item]
                if hapData[item][1]:
                    PhFrExtPOSig[item] = sigDict[pfx+item]            
            if hapData['Pref.Ori.'][0] == 'MD':
                hapData['Pref.Ori.'][1] = parmDict[pfx+'MD']
                if hapData['Pref.Ori.'][2]:
                    PhFrExtPOSig[item] = sigDict[pfx+item]
            else:                           #'SH' spherical harmonics
                for item in hapData['Pref.Ori.'][5]:
                    hapData['Pref.Ori.'][5][item] = parmDict[pfx+item]
                    if hapData['Pref.Ori.'][2]:
                        PhFrExtPOSig[item] = sigDict[pfx+item]
#            print '\n Phase fraction  : %10.4f, sig %10.4f'%(hapData['Scale'][0],PhFrExtPOSig['Scale'])
#            print ' Extinction coeff: %10.4f, sig %10.4f'%(hapData['Extinction'][0],PhFrExtPOSig['Extinction'])
#            if hapData['Pref.Ori.'][0] == 'MD':
#                Ax = hapData['Pref.Ori.'][3]
#                print ' March-Dollase PO: %10.4f'%(hapData['Pref.Ori.'][1]),' Refine?',hapData['Pref.Ori.'][2], \
#                    ' Axis: %d %d %d'%(Ax[0],Ax[1],Ax[2]) 
               
            SizeMuStrSig = {'Mustrain':[[0,0],[0 for i in range(len(hapData['Mustrain'][4]))]],
                'Size':[[0,0],[0 for i in range(len(hapData['Size'][4]))]]}                  
            for item in ['Mustrain','Size']:
                if hapData[item][0] in ['isotropic','uniaxial']:
                    hapData[item][1][0] = parmDict[pfx+item+':0']
                    if hapData[item][2][0]: 
                        SizeMuStrSig[item][0][0] = sigDict[pfx+item+':0']
                    if hapData[item][0] == 'uniaxial':
                        hapData[item][1][1] = parmDict[pfx+item+':1']
                        if hapData[item][2][1]:
                            SizeMuStrSig[item][0][1] = sigDict[pfx+item+':1']
                else:       #generalized for mustrain or ellipsoidal for size
                    for i in range(len(hapData[item][4])):
                        sfx = ':'+str(i)
                        hapData[item][4][i] = parmDict[pfx+item+sfx]
                        if hapData[item][5][i]:
                            SizeMuStrSig[item][1][i] = sigDict[pfx+item+sfx]
                            
            PrintSizeAndSig(hapData['Size'],SizeMuStrSig['Size'])
            PrintMuStrainAndSig(hapData['Mustrain'],SizeMuStrSig['Mustrain'],SGData)
    
def GetHistogramData(Histograms):
    
    def GetBackgroundParms(hId,Background):
        bakType,bakFlag = Background[:2]
        backVals = Background[3:]
        backNames = [':'+str(hId)+':Back:'+str(i) for i in range(len(backVals))]
        if bakFlag:                                 #returns backNames as varyList = backNames
            return bakType,dict(zip(backNames,backVals)),backNames
        else:                                       #no background varied; varyList = []
            return bakType,dict(zip(backNames,backVals)),[]
        
    def GetInstParms(hId,Inst):
        insVals,insFlags,insNames = Inst[1:4]
        dataType = insVals[0]
        instDict = {}
        insVary = []
        pfx = ':'+str(hId)+':'
        for i,flag in enumerate(insFlags):
            insName = pfx+insNames[i]
            instDict[insName] = insVals[i]
            if flag:
                insVary.append(insName)
        instDict[pfx+'X'] = max(instDict[pfx+'X'],0.01)
        instDict[pfx+'Y'] = max(instDict[pfx+'Y'],0.01)
        instDict[pfx+'SH/L'] = max(instDict[pfx+'SH/L'],0.0005)
        return dataType,instDict,insVary
        
    def GetSampleParms(hId,Sample):
        sampVary = []
        hfx = ':'+str(hId)+':'        
        sampDict = {hfx+'Gonio. radius':Sample['Gonio. radius']}
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
        print '\n Background function: ',Background[0],' Refine?',bool(Background[1])
        line = ' Coefficients: '
        for back in Background[3:]:
            line += '%10.3f'%(back)
        print line 
        
    def PrintInstParms(Inst):
        print '\n Instrument Parameters:'
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        instNames = Inst[3][1:]
        for i,name in enumerate(instNames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (Inst[1][i+1])
            if name in ['Lam1','Lam2','Azimuth']:
                varstr += 12*' '
            else:
                varstr += '%12s' % (str(bool(Inst[2][i+1])))
        print ptlbls
        print ptstr
        print varstr
        
    def PrintSampleParms(Sample):
        print '\n Sample Parameters:'
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

        print ptlbls
        print ptstr
        print varstr
        

    histDict = {}
    histVary = []
    controlDict = {}
    for histogram in Histograms:
        Histogram = Histograms[histogram]
        hId = Histogram['hId']
        pfx = ':'+str(hId)+':'
        controlDict[pfx+'Limits'] = Histogram['Limits'][1]
        
        Background = Histogram['Background'][0]
        Type,bakDict,bakVary = GetBackgroundParms(hId,Background)
        controlDict[pfx+'bakType'] = Type
        histDict.update(bakDict)
        histVary += bakVary
        
        Inst = Histogram['Instrument Parameters']
        Type,instDict,insVary = GetInstParms(hId,Inst)
        controlDict[pfx+'histType'] = Type
        histDict.update(instDict)
        histVary += insVary
        
        Sample = Histogram['Sample Parameters']
        Type,sampDict,sampVary = GetSampleParms(hId,Sample)
        controlDict[pfx+'instType'] = Type
        histDict.update(sampDict)
        histVary += sampVary

        print '\n Histogram: ',histogram,' histogram Id: ',hId
        print 135*'-'
        Units = {'C':' deg','T':' msec'}
        units = Units[controlDict[pfx+'histType'][2]]
        Limits = controlDict[pfx+'Limits']
        print ' Instrument type: ',Sample['Type']
        print ' Histogram limits: %8.2f%s to %8.2f%s'%(Limits[0],units,Limits[1],units)     
        PrintSampleParms(Sample)
        PrintInstParms(Inst)
        PrintBackground(Background)
        
    return histVary,histDict,controlDict
    
def SetHistogramData(parmDict,sigDict,Histograms):
    
    def SetBackgroundParms(pfx,Background,parmDict,sigDict):
        lenBack = len(Background[3:])
        backSig = [0 for i in range(lenBack)]
        for i in range(lenBack):
            Background[3+i] = parmDict[pfx+'Back:'+str(i)]
            if Background[1]:
                backSig[i] = sigDict[pfx+'Back:'+str(i)]
        return backSig
        
    def SetInstParms(pfx,Inst,parmDict,sigDict):
        insVals,insFlags,insNames = Inst[1:4]
        instSig = [0 for i in range(len(insVals))]
        for i,flag in enumerate(insFlags):
            insName = pfx+insNames[i]
            insVals[i] = parmDict[insName]
            if flag:
                instSig[i] = sigDict[insName]
        return instSig
        
    def SetSampleParms(pfx,Sample,parmDict,sigDict):
        if 'Bragg' in Sample['Type']:             #Bragg-Brentano
            sampSig = [0 for i in range(3)]
            for i,item in enumerate(['Scale','Shift','Transparency']):       #surface roughness?, diffuse scattering?
                Sample[item][0] = parmDict[pfx+item]
                if Sample[item][1]:
                    sampSig[i] = sigDict[pfx+item]
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            sampSig = [0 for i in range(4)]
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                Sample[item][0] = parmDict[pfx+item]
                if Sample[item][1]:
                    sampSig[i] = sigDict[pfx+item]
        return sampSig
        
    def PrintBackgroundSig(Background,backSig):
        print '\n Background function: ',Background[0]
        valstr = ' value : '
        sigstr = ' sig   : '
        for i,back in enumerate(Background[3:]):
            valstr += '%10.4f'%(back)
            if Background[1]:
                sigstr += '%10.4f'%(backSig[i])
            else:
                sigstr += 10*' '
        print valstr
        print sigstr 
        
    def PrintInstParmsSig(Inst,instSig):
        print '\n Instrument Parameters:'
        ptlbls = ' names :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        instNames = Inst[3][1:]
        for i,name in enumerate(instNames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (Inst[1][i+1])
            if instSig[i+1]:
                sigstr += '%12.6f' % (instSig[i+1])
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr
        
    def PrintSampleParmsSig(Sample,sampleSig):
        print '\n Sample Parameters:'
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        if 'Bragg' in Sample['Type']:
            for i,item in enumerate(['Scale','Shift','Transparency']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                sigstr += '%14.4f'%(sampleSig[i])
            
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            for i,item in enumerate(['Scale','Absorption','DisplaceX','DisplaceY']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                sigstr += '%14.4f'%(sampleSig[i])

        print ptlbls
        print ptstr
        print sigstr
        
    for histogram in Histograms:
        if 'PWDR' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            Background = Histogram['Background'][0]
            backSig = SetBackgroundParms(pfx,Background,parmDict,sigDict)
            
            Inst = Histogram['Instrument Parameters']
            instSig = SetInstParms(pfx,Inst,parmDict,sigDict)
        
            Sample = Histogram['Sample Parameters']
            sampSig = SetSampleParms(pfx,Sample,parmDict,sigDict)

            print '\n Histogram: ',histogram,' histogram Id: ',hId
            print 135*'-'
            print ' Instrument type: ',Sample['Type']
            PrintSampleParmsSig(Sample,sampSig)
            PrintInstParmsSig(Inst,instSig)
            PrintBackgroundSig(Background,backSig)


def SetupSFcalc(General,Atoms):
    ''' setup data for use in StructureFactor; mostly rearranging arrays
    input:
        General = dictionary of general phase info.
        Atoms = list of atom parameters
    returns:
        G = reciprocal metric tensor
        StrData = list of arrays; one entry per atom:
            T = atom types
            R = refinement flags, e.g. 'FXU'
            F = site fractions
            X = atom coordinates as numpy array
            M = atom multiplicities
            IA = isotropic/anisotropic thermal flags
            Uiso = isotropic thermal parameters if IA = 'I'; else = 0
            Uij = numpy array of 6 anisotropic thermal parameters if IA='A'; else = zeros
    '''            
    SGData = General['SGData']
    Cell = General['Cell']
    G,g = G2lat.cell2Gmat(Cell[1:7])        #skip refine & volume; get recip & real metric tensors
    cx,ct,cs,cia = General['AtomPtrs']
    X = [];F = [];T = [];IA = [];Uiso = [];Uij = [];R = [];M = []
    for atom in Atoms:
        T.append(atom[ct])
        R.append(atom[ct+1])
        F.append(atom[cx+3])
        X.append(np.array(atom[cx:cx+3]))
        M.append(atom[cia-1])
        IA.append(atom[cia])
        Uiso.append(atom[cia+1])
        Uij.append(np.array(atom[cia+2:cia+8]))
    return G,[T,R,np.array(F),np.array(X),np.array(M),IA,np.array(Uiso),np.array(Uij)]
            
def StructureFactor(H,G,SGData,StrData,FFtable):
    ''' Compute structure factor for a single h,k,l
    H = np.array(h,k,l)
    G = reciprocal metric tensor
    SGData = space group info. dictionary output from SpcGroup
    StrData = [
        [atom types], 
        [refinement flags], 
        [site fractions],
        np.array(coordinates), 
        [multiplicities], 
        [I/A flag], 
        [Uiso], 
        np.array(Uij)]
    FFtable = dictionary of form factor coeff. for atom types used in StrData
    '''
    twopi = 2.0*math.pi
    twopisq = 2.0*math.pi**2
    SQ = G2lat.calc_rDsq2(H,G)/4.0          # SQ = (sin(theta)/lambda)**2
    SQfactor = 8.0*SQ*math.pi**2
    print 'SQ',SQfactor
    FF = {}
    for type in FFtable.keys():
        FF[type] = G2el.ScatFac(FFtable[type],SQ)
    print 'FF',FF 
    iabsnt,mulp,Uniq,Phs = G2spc.GenHKL(H,SGData,False)
    fa = [0,0,0]        #real
    fb = [0,0,0]        #imaginary
    if not iabsnt:
        phase = twopi*np.inner(Uniq,StrData[3])       #2pi[hx+ky+lz] for each atom in each equiv. hkl
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = StrData[2]*StrData[4]/mulp
        Tiso = np.multiply(StrData[6],SQfactor)
        Tuij = np.multiply(-SQfactor,np.inner(H,np.inner(G2spc.Uij2U(StrData[7]),H)))
        print 'sinp',sinp
        print 'cosp',cosp
        print 'occ',occ
        print 'Tiso',Tiso
        print 'Tuij',Tuij
    else:
        print 'Absent'
        
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def getPowderProfile(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
    def GetSampleGam(refl,wave,G,phfx,calcControls,parmDict,sizeEllipse):
        costh = cosd(refl[5]/2.)
        if calcControls[phfx+'SizeType'] == 'isotropic':
            gam = 1.8*wave/(np.pi*parmDict[phfx+'Size:0']*costh)
        elif calcControls[phfx+'SizeType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            gam = (1.8*wave/np.pi)/parmDict[phfx+'Size:0']*parmDict[phfx+'Size:1']
            gam *= np.sqrt((cosP*parmDict[phfx+'Size:1'])**2+(sinP*parmDict[phfx+'Size:0'])**2)/costh
        else:           #ellipsoidal crystallites
            H = np.array(refl[:3])
            gam += 1.8*wave/(np.pi*costh*np.inner(H,np.inner(sizeEllipse,H)))            
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            gam += 0.018*parmDict[phfx+'Mustrain:0']*tand(refl[5]/2.)/np.pi
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain:0']
            Sa = parmDict[phfx+'Mustrain:1']
            gam += 0.018*Si*Sa*tand(refl[5]/2.)/(np.pi*np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
        else:       #generalized - P.W. Stephens model
            pwrs = calcControls[phfx+'MuPwrs']
            sum = 0
            for i,pwr in enumerate(pwrs):
                sum += parmDict[phfx+'Mustrain:'+str(i)]*refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
            gam += 0.018*refl[4]**2*tand(refl[5]/2.)*sum            
        return gam
        
    def GetIntensityCorr(refl,phfx,hfx,calcControls,parmDict):
        Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3]               #scale*multiplicity
        Icorr *= G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])[0]
        
        return Icorr
        
    def GetReflPos(refl,wave,G,hfx,calcControls,parmDict):
        h,k,l = refl[:3]
        dsq = 1./G2lat.calc_rDsq2(np.array([h,k,l]),G)
        d = np.sqrt(dsq)
        pos = 2.0*asind(wave/(2.0*d))+parmDict[hfx+'Zero']
        const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            pos -= const*(4.*parmDict[hfx+'Shift']*cosd(pos/2.0)+ \
                1.e-7*parmDict[hfx+'Transparency']*sind(pos))            #trans(=1/mueff) in Angstroms
        else:               #Debye-Scherrer - simple but maybe not right
            pos -= const*(parmDict[hfx+'DisplaceX']*cosd(pos)+parmDict[hfx+'DisplaceY']*sind(pos))
        return pos
    
    def GetReflSIgGam(refl,wave,G,hfx,phfx,calcControls,parmDict,sizeEllipse):
        U = parmDict[hfx+'U']
        V = parmDict[hfx+'V']
        W = parmDict[hfx+'W']
        X = parmDict[hfx+'X']
        Y = parmDict[hfx+'Y']
        tanPos = tand(refl[5]/2.0)
        sig = U*tanPos**2+V*tanPos+W        #save peak sigma
        gam = X/cosd(refl[5]/2.0)+Y*tanPos+GetSampleGam(refl,wave,G,phfx,calcControls,parmDict,sizeEllipse) #save peak gamma
        return sig,gam
                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    yb = G2pwd.getBackground(hfx,parmDict,bakType,x)
    yc = np.zeros_like(yb)
        
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.0005)
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
        A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        sizeEllipse = []
        if calcControls[phfx+'SizeType'] == 'ellipsoidal':
            sizeEllipse = G2lat.U6toUij([parmDIct[phfx+'Size:%d'%(i)] for i in range(6)])
        for refl in refList:
            if 'C' in calcControls[hfx+'histType']:
                h,k,l = refl[:3]
                refl[5] = GetReflPos(refl,wave,G,hfx,calcControls,parmDict)         #corrected reflection position
                refl[6:8] = GetReflSIgGam(refl,wave,G,hfx,phfx,calcControls,parmDict,sizeEllipse)    #peak sig & gam
                Icorr = GetIntensityCorr(refl,phfx,hfx,calcControls,parmDict)
                if 'Pawley' in Phase['General']['Type']:
                    try:
                        refl[8] = abs(parmDict[pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])])
                    except KeyError:
                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                else:
                    raise ValueError       #wants strctrfacr calc here
                Wd,fmin,fmax = G2pwd.getWidths(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    return yc,yb
                yc[iBeg:iFin] += Icorr*refl[8]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])    #>90% of time spent here
                if Ka2:
                    pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    Wd,fmin,fmax = G2pwd.getWidths(pos2,refl[6],refl[7],shl)
                    iBeg = np.searchsorted(x,pos2-fmin)
                    iFin = np.searchsorted(x,pos2+fmax)
                    yc[iBeg:iFin] += Icorr*refl[8]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg:iFin])        #and here
            else:
                raise ValueError
    return yc,yb    
            
def getPowderProfileDerv(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
    def GetSampleGamDerv(refl,wave,G,phfx,calcControls,parmDict,sizeEllipse):
        gamDict = {}
        costh = cosd(refl[5]/2.)
        tanth = tand(refl[5]/2.)
        if calcControls[phfx+'SizeType'] == 'isotropic':
            gam = 1.8*wave/(np.pi*parmDict[phfx+'Size:0']*costh)
            gamDict[phfx+'Size:0'] = -gam/parmDict[phfx+'Size:0']
        elif calcControls[phfx+'SizeType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'SizeAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Size:0']
            Sa = parmDict[phfx+'Size:1']
            gami = (1.8*wave/np.pi)/(Si*Sa)
            sqtrm = np.sqrt((cosP*Sa)**2+(sinP*Si)**2)
            gam = gami*sqtrm/costh            
            gamDict[phfx+'Size:0'] = gami*Si*sinP**2/(sqtrm*costh)-gam/Si
            gamDict[phfx+'Size:1'] = gami*Sa*cosP**2/(sqtrm*costh)-gam/Sa          
        else:           #ellipsoidal crystallites - do numerically?
            H = np.array(refl[:3])
            gam = 1.8*wave/(np.pi*costh*np.inner(H,np.inner(sizeEllipse,H)))
                        
        if calcControls[phfx+'MustrainType'] == 'isotropic':
            gamDict[phfx+'Mustrain:0'] =  0.018*tanth/np.pi            
        elif calcControls[phfx+'MustrainType'] == 'uniaxial':
            H = np.array(refl[:3])
            P = np.array(calcControls[phfx+'MustrainAxis'])
            cosP,sinP = G2lat.CosSinAngle(H,P,G)
            Si = parmDict[phfx+'Mustrain:0']
            Sa = parmDict[phfx+'Mustrain:1']
            gami = 0.018*Si*Sa*tanth/np.pi
            sqtrm = np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
            gam = gami/sqtrm
            gamDict[phfx+'Mustrain:0'] = gam/Si-gami*Si*cosP**2/sqtrm**3
            gamDict[phfx+'Mustrain:1'] = gam/Sa-gami*Sa*sinP**2/sqtrm**3
        else:       #generalized - P.W. Stephens model
            pwrs = calcControls[phfx+'MuPwrs']
            const = 0.018*refl[4]**2*tanth
            for i,pwr in enumerate(pwrs):
                gamDict[phfx+'Mustrain:'+str(i)] = const*refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
        return gamDict
        
    def GetIntensityDerv(refl,phfx,hfx,calcControls,parmDict):
        Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3]               #scale*multiplicity
        pola,dpdPola = G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])
        dIdpola = Icorr*dpdPola
        Icorr *= pola
        dIdsh = Icorr/parmDict[hfx+'Scale']
        dIdsp = Icorr/parmDict[phfx+'Scale']
        
        return Icorr,dIdsh,dIdsp,dIdpola
        
    def GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict):
        dpr = 180./np.pi
        h,k,l = refl[:3]
        dstsq = G2lat.calc_rDsq(np.array([h,k,l]),A)
        dst = np.sqrt(dstsq)
        pos = refl[5]
        const = dpr/np.sqrt(1.0-wave*dst/4.0)
        dpdw = const*dst
        dpdA = np.array([h**2,k**2,l**2,h*k,h*l,k*l])
        dpdA *= const*wave/(2.0*dst)
        dpdZ = 1.0
        const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
        if 'Bragg' in calcControls[hfx+'instType']:
            dpdSh = -4.*const*cosd(pos/2.0)
            dpdTr = -1.e-7*const*sind(pos)
            return dpdA,dpdw,dpdZ,dpdSh,dpdTr,0.,0.
        else:               #Debye-Scherrer - simple but maybe not right
            dpdXd = -const*cosd(pos)
            dpdYd = -const*sind(pos)
            return dpdA,dpdw,dpdZ,0.,0.,dpdXd,dpdYd
            
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
            return [[pfx+'A0',dpdA[0]+dpdA[1]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[3]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['3R', '3mR']:
            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]],[pfx+'A3',dpdA[3]+dpdA[4]+dpdA[5]]]                       
        elif SGData['SGLaue'] in ['m3m','m3']:
            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]]]
    
    lenX = len(x)                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    dMdv = np.zeros(shape=(len(varylist),len(x)))
    if hfx+'Back:0' in varylist:
        dMdb = G2pwd.getBackgroundDerv(hfx,parmDict,bakType,x)
        bBpos =varylist.index(hfx+'Back:0')
        dMdv[bBpos:bBpos+len(dMdb)] = dMdb
        
    if 'C' in calcControls[hfx+'histType']:    
        dx = x[1]-x[0]
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
        sizeEllipse = []
        if calcControls[phfx+'SizeType'] == 'ellipsoidal':
            sizeEllipse = G2lat.U6toUij([parmDIct[phfx+'Size:%d'%(i)] for i in range(6)])
        for iref,refl in enumerate(refList):
            if 'C' in calcControls[hfx+'histType']:
                Icorr,dIdsh,dIdsp,dIdpola = GetIntensityDerv(refl,phfx,hfx,calcControls,parmDict)
                hkl = refl[:3]
                pos = refl[5]
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                dsdU = tanth**2
                dsdV = tanth
                dsdW = 1.0
                dgdX = 1.0/costh
                dgdY = tanth
                Wd,fmin,fmax = G2pwd.getWidths(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                dMdpk = np.zeros(shape=(6,len(x)))
                dMdipk = G2pwd.getdFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])
                for i in range(1,5):
                    dMdpk[i][iBeg:iFin] += 100.*dx*Icorr*refl[8]*dMdipk[i]
                dMdpk[0][iBeg:iFin] += 100.*dx*Icorr*dMdipk[0]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4]}
                if Ka2:
                    pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    kdelt = int((pos2-refl[5])/dx)               
                    iBeg = min(lenX,iBeg+kdelt)
                    iFin = min(lenX,iFin+kdelt)
                    if iBeg-iFin:
                        dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6],refl[7],shl,xdata[iBeg:iFin])
                        for i in range(1,5):
                            dMdpk[i][iBeg:iFin] += 100.*dx*Icorr*refl[8]*kRatio*dMdipk2[i]
                        dMdpk[0][iBeg:iFin] += 100.*dx*kRatio*dMdipk2[0]
                        dMdpk[5][iBeg:iFin] += 100.*dx*dMdipk2[0]
                        dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':dMdpk[5]*Icorr*refl[8]}
                try:
                    idx = varylist.index(pfx+'PWLref:'+str(iref))
                    dMdv[idx] = dervDict['int']
                except ValueError:
                    pass
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY = GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'U':[dsdU,'sig'],hfx+'V':[dsdV,'sig'],hfx+'W':[dsdW,'sig'],
                    hfx+'X':[dgdX,'gam'],hfx+'Y':[dgdY,'gam'],hfx+'SH/L':[1.0,'shl'],
                    hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
                    hfx+'DisplaceY':[dpdY,'pos'],}
                for name in names:
                    if name in varylist:
                        item = names[name]
                        dMdv[varylist.index(name)] += item[0]*dervDict[item[1]]
                cellDervNames = cellVaryDerv(pfx,SGData,dpdA)
                for name,dpdA in cellDervNames:
                    if name in varylist:
                        dMdv[varylist.index(name)] += dpdA*dervDict['pos']
                gamDict = GetSampleGamDerv(refl,wave,G,phfx,calcControls,parmDict,sizeEllipse)
                for name in gamDict:
                    if name in varylist:
                        dMdv[varylist.index(name)] += gamDict[name]*dervDict['gam']                       
            else:
                raise ValueError
            
    return dMdv    
                    
def Refine(GPXfile):
    
    def dervRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):
        parmdict.update(zip(varylist,values))
        G2mv.Dict2Map(parmDict)
        Histograms,Phases = HistoPhases
        dMdv = np.empty(0)
        for histogram in Histograms:
            if 'PWDR' in histogram[:4]:
                Histogram = Histograms[histogram]
                hId = Histogram['hId']
                hfx = ':%d:'%(hId)
                Limits = calcControls[hfx+'Limits']
                x,y,w,yc,yb,yd = Histogram['Data']
                xB = np.searchsorted(x,Limits[0])
                xF = np.searchsorted(x,Limits[1])
                dMdvh = np.sqrt(w[xB:xF])*getPowderProfileDerv(parmdict,x[xB:xF],
                    varylist,Histogram,Phases,calcControls,pawleyLookup)
                if len(dMdv):
                    dMdv = np.concatenate((dMdv,dMdvh))
                else:
                    dMdv = dMdvh
        return dMdv
    
    def errRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):        
        parmdict.update(zip(varylist,values))
        G2mv.Dict2Map(parmDict)
        Histograms,Phases = HistoPhases
        M = np.empty(0)
        sumwYo = 0
        Nobs = 0
        for histogram in Histograms:
            if 'PWDR' in histogram[:4]:
                Histogram = Histograms[histogram]
                hId = Histogram['hId']
                hfx = ':%d:'%(hId)
                Limits = calcControls[hfx+'Limits']
                x,y,w,yc,yb,yd = Histogram['Data']
                yc *= 0.0                           #zero full calcd profiles
                yb *= 0.0
                yd *= 0.0
                xB = np.searchsorted(x,Limits[0])
                xF = np.searchsorted(x,Limits[1])
                Nobs += xF-xB
                Histogram['sumwYo'] = np.sum(w[xB:xF]*y[xB:xF]**2)
                sumwYo += Histogram['sumwYo']
                yc[xB:xF],yb[xB:xF] = getPowderProfile(parmdict,x[xB:xF],
                    varylist,Histogram,Phases,calcControls,pawleyLookup)
                yc[xB:xF] *= parmDict[hfx+'Scale']
                yc[xB:xF] += yb[xB:xF]
                yd[xB:xF] = yc[xB:xF]-y[xB:xF]          #yc-yo then all dydv have no '-' needed
                Histogram['sumwYd'] = np.sum(np.sqrt(w[xB:xF])*(yd[xB:xF]))
                M = np.concatenate((M,np.sqrt(w[xB:xF])*(yd[xB:xF])))
        Histograms['sumwYo'] = sumwYo
        Histograms['Nobs'] = Nobs
        Rwp = min(100.,np.sqrt(np.sum(M**2)/sumwYo)*100.)
        if dlg:
            GoOn = dlg.Update(Rwp,newmsg='%s%8.3f%s'%('Powder profile Rwp =',Rwp,'%'))[0]
            if not GoOn:
                return -M           #abort!!
        return M
    
    ShowBanner()
    varyList = []
    parmDict = {}
    calcControls = {}    
    Controls = GetControls(GPXfile)
    ShowControls(Controls)            
    Histograms,Phases = GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no histograms to refine! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception
    phaseVary,phaseDict,pawleyLookup = GetPhaseData(Phases)
    hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histograms)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = GetHistogramData(Histograms)
    calcControls.update(controlDict)
    varyList = phaseVary+histVary+hapVary
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    constrDict,constrFlag,fixedList = G2mv.InputParse([])        #constraints go here?
    groups,parmlist = G2mv.GroupConstraints(constrDict)
    G2mv.GenerateConstraints(groups,parmlist,constrDict,constrFlag,fixedList)
    G2mv.Map2Dict(parmDict,varyList)

    while True:
        begin = time.time()
        values =  np.array(Dict2Values(parmDict, varyList))
        dlg = wx.ProgressDialog('Residual','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        Ftol = Controls['min dM/M']
        try:
            if Controls['deriv type'] == 'analytic':
                result = so.leastsq(errRefine,values,Dfun=dervRefine,full_output=True,ftol=Ftol,col_deriv=True,
                    args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/2)                
            else:           #'numeric'
                result = so.leastsq(errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,
                    args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/len(varyList))
        finally:
            dlg.Destroy()
        runtime = time.time()-begin
        chisq = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict)
        Rwp = np.sqrt(chisq/Histograms['sumwYo'])*100.      #to %
        GOF = chisq/(Histograms['Nobs']-len(varyList))
        print '\n Refinement results:'
        print 135*'-'
        print 'Number of function calls:',result[2]['nfev'],' Number of observations: ',Histograms['Nobs'],' Number of parameters: ',len(varyList)
        print 'Refinement time = %8.3fs, %8.3fs/cycle'%(runtime,runtime/ncyc)
        print 'Rwp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rwp,chisq,GOF)
        try:
            sig = np.sqrt(np.diag(result[1])*GOF)
            if np.any(np.isnan(sig)):
                print '*** Least squares aborted - some invalid esds possible ***'
            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except ValueError:          #result[1] is None on singular matrix
            print '**** Refinement failed - singular matrix ****'
            Ipvt = result[2]['ipvt']
            for i,ipvt in enumerate(Ipvt):
                if not np.sum(result[2]['fjac'],axis=1)[i]:
                    print 'Removing parameter: ',varyList[ipvt-1]
                    del(varyList[ipvt-1])
                    break

    sigDict = dict(zip(varyList,sig))
    SetPhaseData(parmDict,sigDict,Phases)
    SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms)
    SetHistogramData(parmDict,sigDict,Histograms)
    SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases)
#for testing purposes!!!
#    import cPickle
#    file = open('structTestdata.dat','wb')
#    cPickle.dump(parmDict,file,1)
#    cPickle.dump(varyList,file,1)
#    for histogram in Histograms:
#        if 'PWDR' in histogram[:4]:
#            Histogram = Histograms[histogram]
#    cPickle.dump(Histogram,file,1)
#    cPickle.dump(Phases,file,1)
#    cPickle.dump(calcControls,file,1)
#    cPickle.dump(pawleyLookup,file,1)
#    file.close()

def main():
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print 'ERROR - ',GPXfile," doesn't exist!"
            exit()
        GPXpath = ospath.dirname(arg[1])
        Refine(GPXfile)
    else:
        print 'ERROR - missing filename'
        exit()
    print "Done"
         
if __name__ == '__main__':
    main()
