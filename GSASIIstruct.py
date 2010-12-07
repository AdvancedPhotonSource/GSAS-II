#GSASIIstructure - structure computation routines
import sys
import os.path as ospath
import numpy as np
import cPickle
import time
import math
import GSASIIpath
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc

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
    HistogramNames = []
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
    print ' Controls:'
    if Controls['bandWidth']:
        print ' Band approximated least squares matrix refinement, width: ',Controls['bandWidth']
    else:
        print ' Full matrix least squares refinement'
    print ' Maximum number of refinement cycles: ',Controls['Ncycles']
    print ' Minimum sum(shift/esd)^2 for convergence: ','%.2f'%(Controls['minSumShftESD'])
    print ' The Marquardt damping factor: ','%.2f'%(Controls['Marquardt'])
    print ' Maximum allowed atom shift: ','%.2f'%(Controls['maxShift']),'A'
    if Controls['restraintWeight']:
        print ' The weights of the restraint histograms will be modified to normalize their contribution','\n'
    else:
        print ' The restraint histogram weights are fixed','\n'
    
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

def GetPhaseData(GPXfile,PhaseName):
    ''' Returns the 'General' and 'Atoms' objects for PhaseName from GSASII gpx file
    input:
        GPXfile = gpx full file name
        PhaseName = phase name
    return:
        General = dictionary of general phase info.
        Atoms = list of atom parameters
        these are returned empty if PhaseName not found
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
                    General = datus[1]['General']
                    Atoms = datus[1]['Atoms']                                            
    file.close()
    return {'General':General,'Atoms':Atoms}
    
def ShowPhaseData(phaseNames,PhaseData):
    print ' Phases:'
    for name in phaseNames:
        General = PhaseData[name]['General']
        Atoms = PhaseData[name]['Atoms']
        print '\n Phase name: ',General['Name']
        SGtext = G2spc.SGPrint(General['SGData'])
        for line in SGtext: print line
        cell = General['Cell']
        print '\n Unit cell: a =','%.5f'%(cell[1]),' b =','%.5f'%(cell[2]),' c =','%.5f'%(cell[3]), \
            ' alpha =','%.3f'%(cell[4]),' beta =','%.3f'%(cell[5]),' gamma =','%.3f'%(cell[6]),' volume =','%.3f'%(cell[7])
        print ' Refine?',cell[0]
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
#        elif General['Type'] == 'magnetic':
#        elif General['Type'] == 'macromolecular':
        

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
    
def GetPWDRdata(GPXfile,PWDRname):
    ''' Returns powder data from GSASII gpx file
    input: 
        GPXfile = .gpx full file name
        PWDRname = powder histogram name as obtained from GetHistogramNames
    return: 
        PWDRdata = powder data list:
        
    '''
    file = open(GPXfile,'rb')
    PWDRdata = []
    while True:
        try:
            data = cPickle.load(file)
        except EOFError:
            break
        datum = data[0]
        if datum[0] == PWDRname:
            PWDRdata = datum[1:][0]
                                                
            
    file.close()
    return PWDRdata#,Limits,InstrumentParms
    
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
        
def Refine(GPXfile):
    ShowBanner()
    Controls = GetControls(GPXfile)
    ShowControls(Controls)
    
    phaseNames = GetPhaseNames(GPXfile)
    phaseData = {}
    for name in phaseNames: 
        phaseData[name] =  GetPhaseData(GPXfile,name)
    ShowPhaseData(phaseNames,phaseData)
       
    histograms = GetHistogramNames(GPXfile)
    for hist in histograms:
        if 'PWDR' in hist[:4]: 
            print hist
            PWDRdata = GetPWDRdata(GPXfile,hist)
            print PWDRdata

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
