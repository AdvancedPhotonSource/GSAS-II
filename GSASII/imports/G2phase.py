# -*- coding: utf-8 -*-
#
'''There are several classes in :mod:`~GSASII.imports.G2phase`. 
The documentation for them follows.
'''

from __future__ import division, print_function
import sys
import os.path
import math
import random as ran
import numpy as np
try:
    import wx
except ImportError:
    wx = None
from .. import GSASIIobj as G2obj
from .. import GSASIIspc as G2spc
from .. import GSASIIlattice as G2lat
try:  # fails on doc build
    R2pisq = 1./(2.*np.pi**2)
except TypeError:
    pass

class PDB_ReaderClass(G2obj.ImportPhase):
    'Routine to import Phase information from a PDB file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.pdb','.ent','.PDB','.ENT'),
            strictExtension=True,
            formatName = 'PDB',
            longFormatName = 'Original Protein Data Bank (.pdb file) import'
            )
    def ContentsValidator(self, filename):
        '''Taking a stab a validating a PDB file
        (look for cell & at least one atom)
        '''
        fp = open(filename,'r')
#        for i,l in enumerate(fp):
#            if l.startswith('CRYST1'):
#                break
#        else:
#            self.errors = 'no CRYST1 record found'
#            fp.close()
#            return False
        for i,l in enumerate(fp):
            if l.startswith('ATOM') or l.startswith('HETATM'):
                fp.close()
                return True
        self.errors = 'no ATOM records found after CRYST1 record'
        fp.close()
        return False

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read a PDB file using :meth:`ReadPDBPhase`'
        self.Phase = self.ReadPDBPhase(filename, ParentFrame)
        return True
        
    def ReadPDBPhase(self,filename,parent=None):
        '''Read a phase from a PDB file.
        '''
        EightPiSq = 8.*math.pi**2
        self.errors = 'Error opening file'
        file = open(filename, 'r')
        Phase = {}
        Title = os.path.basename(filename)
        RES = Title[:3]
        
        Compnd = ''
        Atoms = []
        A = np.zeros(shape=(3,3))
        S = file.readline()
        line = 1
        SGData = None
        cell = None
        Dummy = True
        Anum = 0
        while S:
            self.errors = 'Error reading at line '+str(line)
            Atom = []
            if 'TITLE' in S[:5]:
                Title = S[10:72].strip()
            elif 'COMPND    ' in S[:10]:
                Compnd = S[10:72].strip()
            elif 'CRYST' in S[:5]:
                Dummy = False
                abc = S[7:34].split()
                angles = S[34:55].split()
                cell=[float(abc[0]),float(abc[1]),float(abc[2]),
                    float(angles[0]),float(angles[1]),float(angles[2])]
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                AA,AB = G2lat.cell2AB(cell)
                SpGrp = S[55:65]
                E,SGData = G2spc.SpcGroup(SpGrp)
                # space group processing failed, try to look up name in table
                if E:
                    SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                    if SpGrpNorm:
                        E,SGData = G2spc.SpcGroup(SpGrpNorm)
                while E:
                    dlg = wx.TextEntryDialog(parent,
                        SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                        'ERROR in space group symbol','',style=wx.OK)
                    if dlg.ShowModal() == wx.ID_OK:
                        SpGrp = dlg.GetValue()
                        E,SGData = G2spc.SpcGroup(SpGrp)
                    else:
                        SGData = G2obj.P1SGData # P 1
                        self.warnings += '\nThe space group was not interpreted and has been set to "P 1".'
                        self.warnings += "Change this in phase's General tab."            
                    dlg.Destroy()
#                SGlines = G2spc.SGPrint(SGData)
#                for l in SGlines: print (l)
            elif 'SCALE' in S[:5]:
                V = S[10:41].split()
                A[int(S[5])-1] = [float(V[0]),float(V[1]),float(V[2])]
            elif 'ATOM' in S[:4] or 'HETATM' in S[:6]:
                if not SGData:
                    self.warnings += '\nThe space group was not read before atoms and has been set to "P 1". '
                    self.warnings += "Change this in phase's General tab."
                    SGData = G2obj.P1SGData # P 1
                    cell = [20.0,20.0,20.0,90.,90.,90.]
                    Volume = G2lat.calc_V(G2lat.cell2A(cell))
                    AA,AB = G2lat.cell2AB(cell)
                    Anum = 1                    
                XYZ = [float(S[31:39]),float(S[39:47]),float(S[47:55])]
                XYZ = np.inner(AB,XYZ)
                XYZ = np.where(abs(XYZ)<0.00001,0,XYZ)
                SytSym,Mult = G2spc.SytSym(XYZ,SGData)[:2]
                Uiso = float(S[61:67])/EightPiSq
                Type = S[76:78].lower()
                if Dummy and S[12:17].strip() == 'CA':
                    Type = 'C'
                Aname = S[12:17].strip()
                if Anum:
                    Aname += '%d'%Anum
                if S[17:20].upper() != 'UNL':
                    RES = S[17:20].upper() 
                Atom = [S[22:27].strip(),RES,S[20:22],
                    Aname,Type.strip().capitalize(),'',XYZ[0],XYZ[1],XYZ[2],
                    float(S[55:61]),SytSym,Mult,'I',Uiso,0,0,0,0,0,0]
                if S[16] in [' ','A','B']:
                    Atom[3] = Atom[3][:3]
                    Atom.append(ran.randint(0,sys.maxsize))
                    Atoms.append(Atom)
                if Anum:
                    Anum += 1
            elif 'ANISOU' in S[:6]:
                Uij = S[30:72].split()
                Uij = [float(Uij[0])/10000.,float(Uij[1])/10000.,float(Uij[2])/10000.,
                    float(Uij[3])/10000.,float(Uij[4])/10000.,float(Uij[5])/10000.]
                Atoms[-1] = Atoms[-1][:14]+Uij
                Atoms[-1][12] = 'A'
                Atoms[-1].append(ran.randint(0,sys.maxsize))
            S = file.readline()
            line += 1
        file.close()
        self.errors = 'Error after read complete'
        if Title:
            PhaseName = Title
        elif Compnd:
            PhaseName = Compnd
        else:
            PhaseName = 'None'
        if not SGData:
            raise self.ImportException("No space group (CRYST entry) found")
        if not cell:
            raise self.ImportException("No cell (CRYST entry) found")
        Phase = G2obj.SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Type'] = 'macromolecular'
        Phase['General']['AtomPtrs'] = [6,4,10,12]
        Phase['Atoms'] = Atoms
        return Phase

class EXP_ReaderClass(G2obj.ImportPhase):
    'Routine to import Phase information from GSAS .EXP files'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.EXP','.exp'),
            strictExtension=True,
            formatName = 'GSAS .EXP',
            longFormatName = 'GSAS Experiment (.EXP file) import'
            )
        
    def ContentsValidator(self, filename):
        'Look for a VERSION tag in 1st line' 
        fp = open(filename,'r')
        if fp.read(13) == '     VERSION ':
            fp.close()
            return True
        self.errors = 'File does not begin with VERSION tag'
        fp.close()
        return False

    def Reader(self,filename,ParentFrame=None,usedRanIdList=[],**unused):
        '''Read a phase from a GSAS .EXP file using 
        :meth:`~EXP_ReaderClass.ReadEXPPhase`
        '''
        self.Phase = G2obj.SetNewPhase(Name='new phase') # create a new empty phase dict
        while self.Phase['ranId'] in usedRanIdList:
            self.Phase['ranId'] = ran.randint(0,sys.maxsize)
        # make sure the ranId is really unique!
        self.MPhase = G2obj.SetNewPhase(Name='new phase') # create a new empty phase dict
        while self.MPhase['ranId'] in usedRanIdList:
            self.MPhase['ranId'] = ran.randint(0,sys.maxsize)
        fp = open(filename,'r')
        self.ReadEXPPhase(ParentFrame, fp)
        fp.close()
        return True

    def ReadEXPPhase(self, G2frame,filepointer):
        '''Read a phase from a GSAS .EXP file.
        '''
        shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
        textureData = {'Order':0,'Model':'cylindrical','Sample omega':[False,0.0],
            'Sample chi':[False,0.0],'Sample phi':[False,0.0],'SH Coeff':[False,{}],
            'SHShow':False,'PFhkl':[0,0,1],'PFxyz':[0,0,1],'PlotType':'Pole figure'}
        shNcof = 0
        S = 1
        NPhas = []
        Expr = [{},{},{},{},{},{},{},{},{}] # GSAS can have at most 9 phases
        for line,S in enumerate(filepointer):
            self.errors = 'Error reading at line '+str(line+1)
            if 'EXPR NPHAS' in S[:12]:
                NPhas = S[12:-1].split()
            if 'CRS' in S[:3]:
                N = int(S[3:4])-1
                Expr[N][S[:12]] = S[12:-1]
        PNames = []
        if not NPhas:
            raise self.ImportException("No EXPR NPHAS record found")
        self.errors = 'Error interpreting file'
        for n,N in enumerate(NPhas):
            if N != '0':
                result = n
                key = 'CRS'+str(n+1)+'    PNAM'
                PNames.append(Expr[n][key])
        if len(PNames) == 0:
            raise self.ImportException("No phases found")            
        elif len(PNames) > 1:
            dlg = wx.SingleChoiceDialog(G2frame, 'Which phase to read?', 'Read phase data', PNames, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelection() # I think this breaks is there are skipped phases. Cant this happen?
            finally:
                dlg.Destroy()        
        EXPphase = Expr[result]
        keyList = list(EXPphase.keys())
        keyList.sort()
        SGData = {}
        MPtype = ''
        if NPhas[result] == '1':
            Ptype = 'nuclear'
        elif NPhas[result] =='2':
            Ptype = 'nuclear'
            MPtype = 'magnetic'
            MagDmin = 1.0
        elif NPhas[result] =='3':
            Ptype = 'magnetic'
            MagDmin = 1.0
        elif NPhas[result] == '4':
            Ptype = 'macromolecular'
        elif NPhas[result] == '10':
            Ptype = 'Pawley'
        else:
            raise self.ImportException("Phase type not recognized")  
            
        for key in keyList:
            if 'PNAM' in key:
               PhaseName = EXPphase[key].strip()
            elif 'ABC   ' in key:
                abc = [float(EXPphase[key][:10]),float(EXPphase[key][10:20]),float(EXPphase[key][20:30])]                        
            elif 'ANGLES' in key:
                angles = [float(EXPphase[key][:10]),float(EXPphase[key][10:20]),float(EXPphase[key][20:30])]                                                
            elif 'SG SYM' in key:
                SpGrp = EXPphase[key][:15].strip()
                E,SGData = G2spc.SpcGroup(SpGrp)
                if E:
                    SGData = G2obj.P1SGData # P 1 -- unlikely to need this!
                    self.warnings += '\nThe GSAS space group was not interpreted(!) and has been set to "P 1".'
                    self.warnings += "Change this in phase's General tab."                       
            elif 'SPNFLP' in key:
                SpnFlp = np.array([int(float(s)) for s in EXPphase[key].split()])
                SpnFlp = np.where(SpnFlp==0,1,SpnFlp)
                SpnFlp = [1,]+list(SpnFlp)
                if SGData['SpGrp'][0] in ['A','B','C','I','R','F']:
                    SpnFlp = list(SpnFlp)+[1,1,1,1]
            elif 'MXDSTR' in key:
                MagDmin = float(EXPphase[key][:10])               
            elif 'OD    ' in key:
                SHdata = EXPphase[key].split() # may not have all 9 values
                SHvals = 9*[0]
                for i in range(9):
                    try:
                        float(SHdata[i])
                        SHvals[i] = SHdata[i]
                    except:
                        pass
                textureData['Order'] = int(SHvals[0])
                textureData['Model'] = shModels[int(SHvals[2])]
                textureData['Sample omega'] = [False,float(SHvals[6])]
                textureData['Sample chi'] = [False,float(SHvals[7])]
                textureData['Sample phi'] = [False,float(SHvals[8])]
                shNcof = int(SHvals[1])
        Volume = G2lat.calc_V(G2lat.cell2A(abc+angles))

        Atoms = []
        MAtoms = []
        Bmat = G2lat.cell2AB(abc+angles)[1]
        if Ptype == 'macromolecular':
            for key in keyList:
                if 'AT' in key[6:8]:
                    S = EXPphase[key]
                    Atom = [S[56:60].strip(),S[50:54].strip().upper(),S[54:56],
                        S[46:51].strip(),S[:8].strip().capitalize(),'',
                        float(S[16:24]),float(S[24:32]),float(S[32:40]),
                        float(S[8:16]),'1',1,'I',float(S[40:46]),0,0,0,0,0,0]
                    XYZ = Atom[6:9]
                    Atom[10],Atom[11] = G2spc.SytSym(XYZ,SGData)[:2]
                    Atom.append(ran.randint(0,sys.maxsize))
                    Atoms.append(Atom)
        else:
            for key in keyList:
                if 'AT' in key:
                    if key[11:] == 'A':
                        S = EXPphase[key]
                    elif key[11:] == 'B':
                        S1 = EXPphase[key]
                        Atom = [S[50:58].strip(),S[:10].strip().capitalize(),'',
                            float(S[10:20]),float(S[20:30]),float(S[30:40]),
                            float(S[40:50]),'',int(S[60:62]),S1[62:63]]
                            #float(S[40:50]),'',int(S[60:62]),S1[130:131]]
                        if Atom[9] == 'I':
                            Atom += [float(S1[0:10]),0.,0.,0.,0.,0.,0.]
                        elif Atom[9] == 'A':
                            Atom += [0.0,
                                float(S1[ 0:10]),float(S1[10:20]),
                                float(S1[20:30]),float(S1[30:40]),
                                float(S1[40:50]),float(S1[50:60])]
                        else:
                            print('Error in line with key: '+key)
                            Atom += [0.,0.,0.,0.,0.,0.,0.]
                        XYZ = Atom[3:6]
                        Atom[7],Atom[8] = G2spc.SytSym(XYZ,SGData)[:2]
                        Atom.append(ran.randint(0,sys.maxsize))
                        Atoms.append(Atom)
                    elif key[11:] == 'M' and key[6:8] == 'AT':
                        S = EXPphase[key]
                        mom = np.array([float(S[:10]),float(S[10:20]),float(S[20:30])])
                        mag = np.sqrt(np.sum(mom**2))
                        mom = np.inner(Bmat,mom)*mag
                        MAtoms.append(Atom)
                        MAtoms[-1] = Atom[:7]+list(mom)+Atom[7:]
                        
        if shNcof:
            shCoef = {}
            nRec = [i+1 for i in range((shNcof-1)//6+1)]
            for irec in nRec:
                ODkey = keyList[0][:6]+'OD'+'%3dA'%(irec)
                indx = EXPphase[ODkey].split()
                ODkey = ODkey[:-1]+'B'
                vals = EXPphase[ODkey].split()
                for i,val in enumerate(vals):
                    key = 'C(%s,%s,%s)'%(indx[3*i],indx[3*i+1],indx[3*i+2])
                    shCoef[key] = float(val)
            textureData['SH Coeff'] = [False,shCoef]
            
        if not SGData:
            raise self.ImportException("No space group found in phase")
        if not abc:
            raise self.ImportException("No cell lengths found in phase")
        if not angles:
            raise self.ImportException("No cell angles found in phase")
        if not Atoms:
            raise self.ImportException("No atoms found in phase")
            
        self.Phase['General'].update({'Type':Ptype,'Name':PhaseName,'Cell':[False,]+abc+angles+[Volume,],'SGData':SGData})
        if MPtype == 'magnetic':
            self.MPhase['General'].update({'Type':'magnetic','Name':PhaseName+' mag','Cell':[False,]+abc+angles+[Volume,],'SGData':SGData})
        else:
            self.MPhase = None
            
        if Ptype =='macromolecular':
            self.Phase['General']['AtomPtrs'] = [6,4,10,12]
            self.Phase['Atoms'] = Atoms
        elif Ptype == 'magnetic':
            self.Phase['General']['AtomPtrs'] = [3,1,10,12]
            self.Phase['General']['SGData']['SGSpin'] = SpnFlp
            self.Phase['General']['MagDmin'] = MagDmin
            self.Phase['Atoms'] = MAtoms            
        else:   #nuclear
            self.Phase['General']['AtomPtrs'] = [3,1,7,9]    
            self.Phase['General']['SH Texture'] = textureData
            self.Phase['Atoms'] = Atoms
        if MPtype =='magnetic':
            self.MPhase['General']['AtomPtrs'] = [3,1,10,12]
            self.MPhase['General']['SGData']['SGSpin'] = SpnFlp
            self.MPhase['General']['MagDmin'] = MagDmin
            self.MPhase['Atoms'] = MAtoms

class JANA_ReaderClass(G2obj.ImportPhase):
    'Routine to import Phase information from a JANA2006 file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.m50','.M50'),
            strictExtension=True,
            formatName = 'JANA m50',
            longFormatName = 'JANA2006 phase import'
            )
    def ContentsValidator(self, filename):
        '''Taking a stab a validating a .m50 file
        (look for cell & at least one atom)
        '''
        fp = open(filename,'r')
        for i,l in enumerate(fp):
            if l.startswith('cell'):
                break
        else:
            self.errors = 'no cell record found'
            fp.close()
            return False
        for i,l in enumerate(fp):
            if l.startswith('spgroup'):
                fp.close()
                return True
        self.errors = 'no spgroup record found after cell record'
        fp.close()
        return False
        
    def Reader(self,filename, ParentFrame=None, **unused):
        'Read a m50 file using :meth:`ReadJANAPhase`'
        self.Phase = self.ReadJANAPhase(filename, ParentFrame)
        return True
        
    def ReadJANAPhase(self,filename,parent=None):
        '''Read a phase from a JANA2006 m50 & m40 files.
        '''
        self.errors = 'Error opening file'
        fp = open(filename, 'r') #contains only cell & spcgroup
        Phase = {}
        Title = os.path.basename(filename)
        Type = 'nuclear'
        Atoms = []
        Atypes = []
        SuperVec = [[0,0,.1],False,4]
        S = fp.readline()
        line = 1
        SGData = None
        SuperSg = ''
        cell = None
        nqi = 0
        version = '2000'
        while S:
            self.errors = 'Error reading at line '+str(line)
            if 'title' in S and S != 'title\n':
                Title = S.split()[1]
            elif 'Jana2006' in S:
                self.warnings += '\nJana2006 file detected'
                version = '2006'
            elif 'cell' in S[:4]:
                cell = S[5:].split()
                cell=[float(cell[0]),float(cell[1]),float(cell[2]),
                    float(cell[3]),float(cell[4]),float(cell[5])]
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                G,g = G2lat.cell2Gmat(cell)
                ast = np.sqrt(np.diag(G))
                Mast = np.multiply.outer(ast,ast)    
                
            elif 'spgroup' in S:
                if 'X' in S:
                    raise self.ImportException("Ad hoc Supersymmetry centering "+S+" not allowed in GSAS-II")            
                SpGrp = S.split()[1]
                SuperSg = ''
                if '(' in SpGrp:    #supercell symmetry - split in 2
                    SuperStr = SpGrp.split('(')
                    SpGrp = SuperStr[0]
                    SuperSg = '('+SuperStr[1]
                SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                E,SGData = G2spc.SpcGroup(SpGrpNorm)
                # space group processing failed, try to look up name in table
                while E:
                    print (G2spc.SGErrors(E))
                    dlg = wx.TextEntryDialog(parent,
                        SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                        'ERROR in space group symbol','',style=wx.OK)
                    if dlg.ShowModal() == wx.ID_OK:
                        SpGrp = dlg.GetValue()
                        E,SGData = G2spc.SpcGroup(SpGrp)
                    else:
                        SGData = G2obj.P1SGData # P 1
                        self.warnings += '\nThe space group was not interpreted and has been set to "P 1".'
                        self.warnings += "Change this in phase's General tab."            
                    dlg.Destroy()
                G2spc.SGPrint(SGData) #silent check of space group symbol
            elif 'qi' in S[:2]:
                if nqi:
                    raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")            
                vec = S.split()[1:]
                SuperVec = [[float(vec[i]) for i in range(3)],False,4]
                nqi += 1
            elif 'atom' in S[:4]:
                Atypes.append(S.split()[1])
            S = fp.readline()
            line += 1
        fp.close()
        #read atoms from m40 file
        if not SGData:
            self.warnings += '\nThe space group was not read before atoms and has been set to "P 1". '
            self.warnings += "Change this in phase's General tab."
            SGData = G2obj.P1SGData # P 1
        waveTypes = ['Fourier','Sawtooth','ZigZag',]
        filename2 = os.path.splitext(filename)[0]+'.m40'
        file2 = open(filename2,'r')
        S = file2.readline()
        line = 1
        self.errors = 'Error reading at line '+str(line)
        nAtoms = int(S.split()[0])
        for i in range(4):
            S = file2.readline()
        for i in range(nAtoms):
            S1N = [0,0,0]            
            Spos = []
            Sadp = []
            Sfrac = []
            Smag = []
            S1 = file2.readline().strip()
            if len(S1) > 55:
                S1N = S1.split()[-3:]   # no. occ, no. pos waves, no. ADP waves
                S1N = [int(i) for i in S1N]
                S1T = list(S1[60:63])
                waveType = waveTypes[int(S1T[1])]
            XYZ = [float(S1[27:36]),float(S1[36:45]),float(S1[45:54])]
            SytSym,Mult = G2spc.SytSym(XYZ,SGData)[:2]
            aType = Atypes[int(S1[9:11])-1]
            Name = S1[:8].strip()
            if S1[11:15].strip() == '1':
                S2 = file2.readline()
                Uiso = float(S2[:9])
                if version == '2000':
                    Uiso = R2pisq*float(Uiso)/4.      #Biso -> Uiso
                Uij = [0,0,0,0,0,0]
                IA = 'I'
            elif S1[11:15].strip() == '2':
                S2 = file2.readline()
                IA = 'A'
                Uiso = 0.
                Uij = [float(S2[:9]),float(S2[9:18]),float(S2[18:27]),
                    float(S2[27:36]),float(S2[36:45]),float(S2[45:54])] #Uij in Jana2006!
                if version == '2000':
                    Uij = R2pisq*G2lat.UijtoU6(G2lat.U6toUij(Uij)/Mast) #these things are betaij in Jana2000! need to convert to Uij
            for i in range(S1N[0]):
                if not i:
                    FS = file2.readline()
                    Sfrac.append(FS[:9])    #'O' or 'delta' = 'length' for crenel
                    if int(S1T[0]):  #"", "Legendre" or "Xharm" in 18:27 for "crenel"!
                        waveType = 'Crenel/Fourier' #all waves 'Fourier' no other choice
                Sfrac.append(file2.readline()[:18]) #if not crenel = Osin & Ocos
                # else Osin & Ocos except last one is X40 = 'Center'
            for i in range(S1N[1]):  
                Spos.append(file2.readline()[:54])
            for i in range(S1N[2]):
                Sadp.append(file2.readline()[:54]+file2.readline())
            if sum(S1N):    #if any waves: skip mystery line?
                file2.readline()
            for i,it in enumerate(Sfrac):
                if not i:
                    if 'Crenel' in waveType:
                        vals = [float(it),float(Sfrac[-1][:9])]
                    else:
                        vals = [float(it),]
                else:
                    vals = [float(it[:9]),float(it[9:18])]
                if 'Crenel' in waveType and i == len(Sfrac)-1:
                    del Sfrac[-1]
                    break                
                Sfrac[i] = [vals,False]
            for i,it in enumerate(Spos):
                if waveType in ['Sawtooth',] and not i:
                    vals = [float(it[:9]),float(it[9:18]),float(it[18:27]),float(it[27:36])]
                else:
                    vals = [float(it[:9]),float(it[9:18]),float(it[18:27]),float(it[27:36]),float(it[36:45]),float(it[45:54])]
                Spos[i] = [vals,False]
            for i,it in enumerate(Sadp):
                vals = [float(it[:9]),float(it[9:18]),float(it[18:27]),float(it[27:36]),float(it[36:45]),float(it[45:54]),
                    float(it[54:63]),float(it[63:72]),float(it[72:81]),float(it[81:90]),float(it[90:99]),float(it[99:108])]
                #these are betaij modulations in Jana2000! need to convert to Uij modulations
                if version == '2000':               
                    vals[:6] = R2pisq*G2lat.UijtoU6(G2lat.U6toUij(vals[:6])/Mast)    #convert sin bij to Uij
                    vals[6:] = R2pisq*G2lat.UijtoU6(G2lat.U6toUij(vals[6:])/Mast)    #convert cos bij to Uij
                Sadp[i] = [vals,False]
            Atom = [Name,aType,'',XYZ[0],XYZ[1],XYZ[2],1.0,SytSym,Mult,IA,Uiso]
            Atom += Uij
            Atom.append(ran.randint(0,sys.maxsize))
            if len(S1) > 55:
                Atom.append({'SS1':{'Sfrac':[waveType,]+Sfrac,'Spos':[waveType,]+Spos,'Sadp':['Fourier',]+Sadp,'Smag':['Fourier',]+Smag}})    #SS2 is for (3+2), etc.
            Atoms.append(Atom)
        file2.close()
        self.errors = 'Error after read complete'
        if not SGData:
            raise self.ImportException("No space group (spcgroup entry) found")
        if not cell:
            raise self.ImportException("No cell found")
        Phase = G2obj.SetNewPhase(Name=Title,SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Type'] = Type
        Phase['General']['Modulated'] = False
        Phase['General']['AtomPtrs'] = [3,1,7,9]
        if len(S1) > 55:
            Phase['General']['Modulated'] = True
            Phase['General']['Super'] = nqi
            Phase['General']['SuperVec'] = SuperVec
            Phase['General']['SuperSg'] = SuperSg
            if SuperSg:
                Phase['General']['SSGData'] = G2spc.SSpcGroup(SGData,SuperSg)[1]
        Phase['Atoms'] = Atoms
        return Phase
    
class PDF_ReaderClass(G2obj.ImportPhase):
    '''Routine to import Phase information from ICDD Powder Diffraction 
    File(r) Card, exported by their software.'''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.str',),
            strictExtension=True,
            formatName = 'ICDD .str',
            longFormatName = 'ICDD PDF Card (.str file) import'
            )
        
    def ContentsValidator(self, filename):
        'Look for a str tag in 1st line' 
        fp = open(filename,'r')
        if fp.read(3) == 'str':
            fp.close()
            return True
        self.errors = 'File does not begin with str tag'
        fp.close()
        return False

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read phase from a ICDD .str file using :meth:`ReadPDFPhase`'
        fp = open(filename,'r')
        self.Phase = self.ReadPDFPhase(ParentFrame, fp)
        fp.close()
        return True

    def ReadPDFPhase(self, G2frame,fp):
        '''Read a phase from a ICDD .str file.
        '''
        EightPiSq = 8.*math.pi**2
        self.errors = 'Error opening file'
        Phase = {}
        Atoms = []
        S = fp.readline()
        line = 1
        SGData = G2obj.P1SGData
        SGset = False
        cell = []
        cellkey = []
        while S:
            if 'space_group' in S:
                break
            S = fp.readline()
            line += 1
        while S:
            self.errors = 'Error reading at line '+str(line)
            if 'phase_name' in S:
                Title = S.split('"')[1]
            elif 'Space group (HMS)' in S or 'space_group' in S:
                SpGrp = S.split()[-1]
                SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                E,SGData = G2spc.SpcGroup(SpGrpNorm)
                if not E: SGset = True
                # space group processing failed, try to look up name in table
                while E:
                    print (G2spc.SGErrors(E))
                    dlg = wx.TextEntryDialog(G2frame,
                        SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                        'ERROR in space group symbol','',style=wx.OK)
                    if dlg.ShowModal() == wx.ID_OK:
                        SpGrp = dlg.GetValue()
                        E,SGData = G2spc.SpcGroup(SpGrp)
                        if not E:
                            SGset = True
                            break
                    else:
                        SGData = G2obj.P1SGData # P 1
                        SGset = False

                    dlg.Destroy()
                G2spc.SGPrint(SGData) #silent check of space group symbol
            elif 'a a_' in S[:7]:
                data = S.split()
                cell.append(float(data[2]))
                cellkey.append(data[1])
            elif 'b b_' in S[:7]:
                data = S.split()
                cell.append(float(data[2]))
                cellkey.append(data[1])
            elif 'b =' in S[:6]:
                data = S.split('=')
                indx = cellkey.index(data[1].split(';')[0])
                cell.append(cell[indx])
            elif 'c c_' in S[:7]:
                data = S.split()
                cell.append(float(data[2]))
            elif 'c =' in S[:6]:
                data = S.split('=')
                indx = cellkey.index(data[1].split(';')[0])
                cell.append(cell[indx])
            elif 'al' in S[:5]:
                cell.append(float(S.split()[1]))
            elif 'be' in S[:5]:
                cell.append(float(S.split()[1]))
            elif 'ga' in S[:5]:
                cell.append(float(S.split()[1]))
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                break
            S = fp.readline()
        S = fp.readline()
        while S:
            if '/*' in S[:5]:
                break
            if 'site' in S[:7]:
                atom = []
                xyzkey = []
                data = S.split()
                atom.append(data[1])    #name
                pos = data.index('occ')+1
                atom.append(data[pos])  #type
                atom.append('')         #refine
                for xid in ['x =','y =','z =']:
                    if xid in S:
                        xpos = S.index(xid)+3
                        xend = xpos+S[xpos:].index(';')
                        if S[xpos:xend] in xyzkey:
                            atom.append(atom[3+xyzkey.index(S[xpos:xend])])
                        else:
                            atom.append(eval(S[xpos:xend]+'.'))
                    else:
                        xpos = data.index(xid[0])+2
                        xyzkey.append(data[xpos-1][1:])
                        atom.append(float(data[xpos]))
                atom.append(float(data[pos+2]))
                SytSym,Mult = G2spc.SytSym(np.array(atom[3:6]),SGData)[:2]
                atom.append(SytSym)
                atom.append(Mult)
                if 'beq' in S:
                    atom.append('I')
                    upos = data.index('beq')
                    atom.append(float(data[upos+2])/EightPiSq)
                    atom += [0.,0.,0.,0.,0.,0.,]
                elif 'ADPs' in S:
                    upos = data.index('ADPs')
                    atom.append('A')
                    atom.append(0.0)
                    for uid in ['Bani11','Bani22','Bani33','Bani12','Bani13','Bani23']:
                        upos = data.index(uid)+1
                        atom.append(float(data[upos])/EightPiSq)
                else:
                    atom.append('I')
                    atom += [0.02,0.,0.,0.,0.,0.,0.,]                    
                atom.append(ran.randint(0,sys.maxsize))
                Atoms.append(atom)
            S = fp.readline()                
        fp.close()
        if not SGset:
            self.warnings += '\nThe space group was not interpreted and has been set to "P 1".'
            self.warnings += " Change this in phase's General tab."
        self.errors = 'Error after read complete'
        if not SGData:
            raise self.ImportException("No space group (spcgroup entry) found")
        if not cell:
            raise self.ImportException("No cell found")
        Phase = G2obj.SetNewPhase(Name=Title,SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Type'] = 'nuclear'
        Phase['General']['AtomPtrs'] = [3,1,7,9]
        Phase['Atoms'] = Atoms
        return Phase
