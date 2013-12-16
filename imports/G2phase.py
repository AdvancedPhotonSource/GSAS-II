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
*Module G2phase: PDB and .EXP*
------------------------------------

A set of short routines to read in phases using routines that were
previously implemented in GSAS-II: PDB and GSAS .EXP file formats

'''

import sys
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
