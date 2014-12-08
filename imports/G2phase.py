# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
#
'''
*Module G2phase: PDB, .EXP & JANA m40,m50*
-------------------------------------------

A set of short routines to read in phases using routines that were
previously implemented in GSAS-II: PDB, GSAS .EXP and JANA m40-m50 file formats

'''

import sys
import os.path
import math
import random as ran
import traceback
import numpy as np
import wx
import GSASIIIO as G2IO
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

class PDB_ReaderClass(G2IO.ImportPhase):
    'Routine to import Phase information from a PDB file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.pdb','.ent','.PDB','.ENT'),
            strictExtension=True,
            formatName = 'PDB',
            longFormatName = 'Original Protein Data Bank (.pdb file) import'
            )
    def ContentsValidator(self, filepointer):
        '''Taking a stab a validating a PDB file
        (look for cell & at least one atom)
        '''
        for i,l in enumerate(filepointer):
            if l.startswith('CRYST1'):
                break
        else:
            self.errors = 'no CRYST1 record found'
            return False
        for i,l in enumerate(filepointer):
            if l.startswith('ATOM'):
                return True
        self.errors = 'no ATOM records found after CRYST1 record'
        return False

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a PDF file using :meth:`ReadPDBPhase`'
        try:
            self.Phase = self.ReadPDBPhase(filename, ParentFrame)
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print 'PDB read error:',detail # for testing
            traceback.print_exc(file=sys.stdout)
            return False
        
    def ReadPDBPhase(self,filename,parent=None):
        '''Read a phase from a PDB file.
        '''
        EightPiSq = 8.*math.pi**2
        self.errors = 'Error opening file'
        file = open(filename, 'Ur')
        Phase = {}
        Title = ''
        Compnd = ''
        Atoms = []
        A = np.zeros(shape=(3,3))
        S = file.readline()
        line = 1
        SGData = None
        cell = None
        while S:
            self.errors = 'Error reading at line '+str(line)
            Atom = []
            if 'TITLE' in S[:5]:
                Title = S[10:72].strip()
            elif 'COMPND    ' in S[:10]:
                Compnd = S[10:72].strip()
            elif 'CRYST' in S[:5]:
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
                    print G2spc.SGErrors(E)
                    dlg = wx.TextEntryDialog(parent,
                        SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                        'ERROR in space group symbol','',style=wx.OK)
                    if dlg.ShowModal() == wx.ID_OK:
                        SpGrp = dlg.GetValue()
                        E,SGData = G2spc.SpcGroup(SpGrp)
                    else:
                        SGData = G2IO.SGData # P 1
                        self.warnings += '\nThe space group was not interpreted and has been set to "P 1".'
                        self.warnings += "Change this in phase's General tab."            
                    dlg.Destroy()
                SGlines = G2spc.SGPrint(SGData)
                for l in SGlines: print l
            elif 'SCALE' in S[:5]:
                V = S[10:41].split()
                A[int(S[5])-1] = [float(V[0]),float(V[1]),float(V[2])]
            elif 'ATOM' in S[:4] or 'HETATM' in S[:6]:
                if not SGData:
                    self.warnings += '\nThe space group was not read before atoms and has been set to "P 1". '
                    self.warnings += "Change this in phase's General tab."
                    SGData = G2IO.SGData # P 1
                XYZ = [float(S[31:39]),float(S[39:47]),float(S[47:55])]
                XYZ = np.inner(AB,XYZ)
                XYZ = np.where(abs(XYZ)<0.00001,0,XYZ)
                SytSym,Mult = G2spc.SytSym(XYZ,SGData)
                Uiso = float(S[61:67])/EightPiSq
                Type = S[12:14].lower()
                if Type[0] in '123456789':
                    Type = Type[1:]
                Atom = [S[22:27].strip(),S[17:20].upper(),S[20:22],
                    S[12:17].strip(),Type.strip().capitalize(),'',XYZ[0],XYZ[1],XYZ[2],
                    float(S[55:61]),SytSym,Mult,'I',Uiso,0,0,0,0,0,0]
                S = file.readline()
                line += 1
                if 'ANISOU' in S[:6]:
                    Uij = S[30:72].split()
                    Uij = [float(Uij[0])/10000.,float(Uij[1])/10000.,float(Uij[2])/10000.,
                        float(Uij[3])/10000.,float(Uij[4])/10000.,float(Uij[5])/10000.]
                    Atom = Atom[:14]+Uij
                    Atom[12] = 'A'
                Atom.append(ran.randint(0,sys.maxint))
                Atoms.append(Atom)
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
        Phase = G2IO.SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Type'] = 'macromolecular'
        Phase['General']['AtomPtrs'] = [6,4,10,12]
        Phase['Atoms'] = Atoms
        return Phase

class EXP_ReaderClass(G2IO.ImportPhase):
    'Routine to import Phase information from GSAS .EXP files'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.EXP','.exp'),
            strictExtension=True,
            formatName = 'GSAS .EXP',
            longFormatName = 'GSAS Experiment (.EXP file) import'
            )
        
    def ContentsValidator(self, filepointer):
        'Look for a VERSION tag in 1st line' 
        if filepointer.read(13) == '     VERSION ':
            return True
        self.errors = 'File does not begin with VERSION tag'
        return False

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a phase from a GSAS .EXP file using :meth:`ReadEXPPhase`'
        try:
            self.Phase = self.ReadEXPPhase(ParentFrame, filepointer)
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print 'GSAS .EXP read error:',detail # for testing
            traceback.print_exc(file=sys.stdout)
            return False

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
                Num = S[12:-1].count('0')
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
        keyList = EXPphase.keys()
        keyList.sort()
        SGData = {}
        if NPhas[result] == '1':
            Ptype = 'nuclear'
        elif NPhas[result] in ['2','3']:
            Ptype = 'magnetic'
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
                    SGData = G2IO.SGData # P 1 -- unlikely to need this!
                    self.warnings += '\nThe GSAS space group was not interpreted(!) and has been set to "P 1".'
                    self.warnings += "Change this in phase's General tab."                       
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
        Atoms = []
        if Ptype == 'nuclear':
            for key in keyList:
                if 'AT' in key:
                    if key[11:] == 'A':
                        S = EXPphase[key]
                    elif key[11:] == 'B':
                        S += EXPphase[key]
                        Atom = [S[50:58].strip(),S[:10].strip().capitalize(),'',
                            float(S[10:20]),float(S[20:30]),float(S[30:40]),
                            float(S[40:50]),'',int(S[60:62]),S[130:131]]
                        if Atom[9] == 'I':
                            Atom += [float(S[68:78]),0.,0.,0.,0.,0.,0.]
                        elif Atom[9] == 'A':
                            Atom += [0.0,float(S[68:78]),float(S[78:88]),
                                float(S[88:98]),float(S[98:108]),
                                float(S[108:118]),float(S[118:128])]
                        XYZ = Atom[3:6]
                        Atom[7],Atom[8] = G2spc.SytSym(XYZ,SGData)
                        Atom.append(ran.randint(0,sys.maxint))
                        Atoms.append(Atom)
        elif Ptype == 'macromolecular':
            for key in keyList:
                if 'AT' in key[6:8]:
                    S = EXPphase[key]
                    Atom = [S[56:60],S[50:54].strip().upper(),S[54:56],
                        S[46:51].strip(),S[:8].strip().capitalize(),'',
                        float(S[16:24]),float(S[24:32]),float(S[32:40]),
                        float(S[8:16]),'1',1,'I',float(S[40:46]),0,0,0,0,0,0]
                    XYZ = Atom[6:9]
                    Atom[10],Atom[11] = G2spc.SytSym(XYZ,SGData)
                    Atom.append(ran.randint(0,sys.maxint))
                    Atoms.append(Atom)
        Volume = G2lat.calc_V(G2lat.cell2A(abc+angles))
        if shNcof:
            shCoef = {}
            nRec = [i+1 for i in range((shNcof-1)/6+1)]
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
        Phase = G2IO.SetNewPhase(Name=PhaseName,SGData=SGData,cell=abc+angles+[Volume,])
        general = Phase['General']
        general['Type'] = Ptype
        if general['Type'] =='macromolecular':
            general['AtomPtrs'] = [6,4,10,12]
        else:
            general['AtomPtrs'] = [3,1,7,9]    
        general['SH Texture'] = textureData
        Phase['Atoms'] = Atoms
        return Phase

class JANA_ReaderClass(G2IO.ImportPhase):
    'Routine to import Phase information from a JANA2006 file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.m50','.M50'),
            strictExtension=True,
            formatName = 'JANA m50',
            longFormatName = 'JANA2006 phase import'
            )
    def ContentsValidator(self, filepointer):
        '''Taking a stab a validating a .m50 file
        (look for cell & at least one atom)
        '''
        for i,l in enumerate(filepointer):
            if l.startswith('cell'):
                break
        else:
            self.errors = 'no cell record found'
            return False
        for i,l in enumerate(filepointer):
            if l.startswith('spgroup'):
                return True
        self.errors = 'no spgroup record found after cell record'
        return False
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a m50 file using :meth:`ReadJANAPhase`'
        try:
            self.Phase = self.ReadJANAPhase(filename, ParentFrame)
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print 'JANA read error:',detail # for testing
            traceback.print_exc(file=sys.stdout)
            return False
        
    def ReadJANAPhase(self,filename,parent=None):
        '''Read a phase from a JANA2006 m50 & m40 files.
        '''
        self.errors = 'Error opening file'
        file = open(filename, 'Ur') #contains only cell & spcgroup
        Phase = {}
        Title = os.path.basename(filename)
        Compnd = ''
        Type = 'nuclear'
        Atoms = []
        Atypes = []
        SuperVec = [[[0,0,.1],False,4],[[0,0,.1],False,4],[[0,0,.1],False,4]]
        S = file.readline()
        line = 1
        SGData = None
        SuperSg = ''
        cell = None
        nqi = 0
        while S:
            self.errors = 'Error reading at line '+str(line)
            if 'title' in S and S != 'title\n':
                Title = S.split()[1]
            elif 'cell' in S[:4]:
                cell = S[5:].split()
                cell=[float(cell[0]),float(cell[1]),float(cell[2]),
                    float(cell[3]),float(cell[4]),float(cell[5])]
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
            elif 'spgroup' in S:
                if 'X' in S:
                    raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")            
                SpGrp = S.split()[1]
                if '(' in SpGrp:    #supercell symmetry - split in 2
                    SuperStr = SpGrp.split('(')
                    SpGrp = SuperStr[0]
                    SuperSg = '('+SuperStr[1]
                SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                E,SGData = G2spc.SpcGroup(SpGrp)
                # space group processing failed, try to look up name in table
                while E:
                    print G2spc.SGErrors(E)
                    dlg = wx.TextEntryDialog(parent,
                        SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                        'ERROR in space group symbol','',style=wx.OK)
                    if dlg.ShowModal() == wx.ID_OK:
                        SpGrp = dlg.GetValue()
                        E,SGData = G2spc.SpcGroup(SpGrp)
                    else:
                        SGData = G2IO.SGData # P 1
                        self.warnings += '\nThe space group was not interpreted and has been set to "P 1".'
                        self.warnings += "Change this in phase's General tab."            
                    dlg.Destroy()
                SGlines = G2spc.SGPrint(SGData)
                for l in SGlines: print l
            elif 'qi' in S[:2]:
                if nqi:
                    raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")            
                Type = 'modulated'
                vec = S.split()[1:]
                SuperVec = [[float(vec[i]) for i in range(3)],False,4]
                nqi += 1
            elif 'atom' in S[:4]:
                Atypes.append(S.split()[1])
            S = file.readline()
            line += 1
        file.close()
        #read atoms from m40 file
        if not SGData:
            self.warnings += '\nThe space group was not read before atoms and has been set to "P 1". '
            self.warnings += "Change this in phase's General tab."
            SGData = G2IO.SGData # P 1
        waveTypes = ['Fourier','Sawtooth','ZigZag',]
        filename2 = os.path.splitext(filename)[0]+'.m40'
        file2 = open(filename2,'Ur')
        S = file2.readline()
        line = 1
        self.errors = 'Error reading at line '+str(line)
        nAtoms = int(S.split()[0])
        for i in range(4):
            S = file2.readline()            
        for i in range(nAtoms):
            S1 = file2.readline()
            S1N = S1.split()[-3:]   # no. occ, no. pos waves, no. ADP waves
            S1N = [int(i) for i in S1N]
            S1T = list(S1[60:63])
            waveType = waveTypes[int(S1T[1])]
            crenelType = ''
            Spos = []
            Sadp = []
            Sfrac = []
            Smag = []
            XYZ = [float(S1[27:36]),float(S1[36:45]),float(S1[45:54])]
            SytSym,Mult = G2spc.SytSym(XYZ,SGData)
            aType = Atypes[int(S1[9:11])-1]
            Name = S1[:8].strip()
            if S1[11:15].strip() == '1':
                S2 = file2.readline()
                Uiso = float(S2[:9])
                Uij = [0,0,0,0,0,0]
                IA = 'I'
            elif S1[11:15].strip() == '2':
                S2 = file2.readline()
                IA = 'A'
                Uiso = 0.
                Uij = [float(S2[:9]),float(S2[9:18]),float(S2[18:27]),
                    float(S2[27:36]),float(S2[36:45]),float(S2[45:54])]
            for i in range(S1N[0]):
                if not i:
                    FS = file2.readline()
                    Sfrac.append(FS[:9])    #'O' or 'delta' = 'length' for crenel
                    if int(S1T[0]):  #"", "Legendre" or "Xharm" in 18:27 for "crenel"!
                        waveType = 'Crenel/Fourier' #all waves 'Fourier' no other choice
                        crenelType = FS[18:27]
                Sfrac.append(file2.readline()[:18]) #if not crenel = Osin & Ocos
                # else Osin & Ocos except last one is X40 = 'Center'
            for i in range(S1N[1]):  
                Spos.append(file2.readline()[:54])
            for i in range(S1N[2]):
                Sadp.append(file2.readline()[:54]+file2.readline())
            if sum(S1N):    #if any waves: skip mystery line?
                file2.readline()
            for i,it in enumerate(Sfrac):
                print i,it
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
                print Sfrac[i]
            for i,it in enumerate(Spos):
                vals = [float(it[:9]),float(it[9:18]),float(it[18:27]),float(it[27:36]),float(it[36:45]),float(it[45:54])]
                Spos[i] = [vals,False]
            for i,it in enumerate(Sadp):
                vals = [float(it[:9]),float(it[9:18]),float(it[18:27]),float(it[27:36]),float(it[36:45]),float(it[45:54]),
                    float(it[54:63]),float(it[63:72]),float(it[72:81]),float(it[81:90]),float(it[90:99]),float(it[99:108])]                
                Sadp[i] = [vals,False]
            Atom = [Name,aType,'',XYZ[0],XYZ[1],XYZ[2],1.0,SytSym,Mult,IA,Uiso]
            Atom += Uij
            Atom.append({'SS1':{'waveType':waveType,'crenelType':crenelType,'Sfrac':Sfrac,'Spos':Spos,'Sadp':Sadp,'Smag':Smag}})    #SS2 is for (3+2), etc.
            Atom.append(ran.randint(0,sys.maxint))
            Atoms.append(Atom)
        file2.close()
        self.errors = 'Error after read complete'
        if not SGData:
            raise self.ImportException("No space group (spcgroup entry) found")
        if not cell:
            raise self.ImportException("No cell found")
        Phase = G2IO.SetNewPhase(Name=Title,SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Type'] = Type
        Phase['General']['Super'] = nqi
        Phase['General']['SuperVec'] = SuperVec
        Phase['General']['SuperSg'] = SuperSg
        Phase['General']['AtomPtrs'] = [3,1,7,9]
        Phase['Atoms'] = Atoms
        return Phase
