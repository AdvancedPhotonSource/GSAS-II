# -*- coding: utf-8 -*-
'''Class to read a phase from a RMCprofile output file
'''
from __future__ import division, print_function
import sys
import os.path
import numpy as np
import random as ran
from .. import GSASIIobj as G2obj
from .. import GSASIIlattice as G2lat

class PhaseReaderClass(G2obj.ImportPhase):
    'Opens a .rmc6f file and pulls out the phase'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.rmc6f',),
            strictExtension=True,
            formatName = 'RMCProfile .rmc6f',
            longFormatName = 'RMCProfile structure (*.rmc6f) file import'
            )
        
    def ContentsValidator(self, filename):
        "Test if the rmc6f file has a CELL record"
        fp = open(filename,'r')
        if fp.readline()[:-1] != '(Version 6f format configuration file)':
            self.errors = 'This is not a valid .rmc6f file.'
            fp.close()
            return False
        fp.close()
        return True

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read a rmc6f file using 
        :meth:`~GSASII.imports.G2phase_rmc6f.PhaseReaderClass.Readrmc6fPhase`
        '''
        self.Phase = self.Readrmc6fPhase(filename, ParentFrame)
        return True

    def Readrmc6fPhase(self,filename,parent=None):
        '''Read a phase from a rmc6f file.
        '''
        self.errors = 'Error opening file'
        fp = open(filename, 'r')
        Phase = {}
        Title = os.path.split(filename)
        import GSASIIctrlGUI as G2G
        G2G.SaveGPXdirectory(Title[0])
        Title = os.path.splitext(Title[1])[0]
        Atoms = []
        S = fp.readline()
        line = 1
        SGData = None
        cell = None
        IA = 'I'
        Uiso = 0.01
        Uij = [0. for i in range(6)]
        while S:
            self.errors = 'Error reading at line '+str(line)
            Atom = []
            if 'Cell' in S[:4]:
                cellRec = S.split(':')[1].split()
                abc = cellRec[:3]
                angles = cellRec[3:]
                cell=[float(abc[0]),float(abc[1]),float(abc[2]),
                    float(angles[0]),float(angles[1]),float(angles[2])]
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                AA,AB = G2lat.cell2AB(cell)
                SGData = G2obj.P1SGData # P 1
            elif 'Atoms' in S[:5]:
                S = fp.readline()[:-1]
                AtRec = S.split()
                for ix,s in enumerate(AtRec):
                    if '.' in s:
                        break       #is points at x
                while S:
                    AtRec = S.split()
                    Atype = AtRec[1]
                    Aname = Atype+AtRec[0]
                    Afrac = 1.0
                    x,y,z = AtRec[ix:ix+3]
                    XYZ = np.array([float(x),float(y),float(z)])
                    SytSym,Mult = '1',1
                    Atom = [Aname,Atype,'',XYZ[0],XYZ[1],XYZ[2],Afrac,SytSym,Mult,IA,Uiso]
                    Atom += Uij
                    Atom.append(ran.randint(0,sys.maxsize))
                    Atoms.append(Atom)
                    S = fp.readline()[:-1]
            S = fp.readline()
            line += 1
        fp.close()
        self.errors = 'Error after read complete'
        Phase = G2obj.SetNewPhase(Name='RMCProfile phase',SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Name'] = Title
        Phase['General']['Type'] = 'nuclear'
        Phase['General']['AtomPtrs'] = [3,1,7,9]
        Phase['Atoms'] = Atoms
        return Phase
