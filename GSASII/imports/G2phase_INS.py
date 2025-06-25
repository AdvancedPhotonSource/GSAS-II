# -*- coding: utf-8 -*-
'''
'''
from __future__ import division, print_function
import sys
import numpy as np
import random as ran
from .. import GSASIIobj as G2obj
from .. import GSASIIspc as G2spc
from .. import GSASIIlattice as G2lat

class PhaseReaderClass(G2obj.ImportPhase):
    'Opens a .INS file and pulls out a selected phase'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.ins','.INS','.res','.RES'),
            strictExtension=True,
            formatName = 'SHELX ins, res',
            longFormatName = 'SHELX input (*.ins, *.res) file import'
            )
        
    def ContentsValidator(self, filename):
        "Test if the ins file has a CELL record"
        fp = open(filename,'r')
        for i,l in enumerate(fp):
            if l.startswith('CELL'):
                break
        else:
            self.errors = 'no CELL record found'
            self.errors = 'This is not a valid .ins file.'
            fp.close()
            return False
        fp.close()
        return True

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a ins file using :meth:`ReadINSPhase`'
        self.Phase = self.ReadINSPhase(filename, ParentFrame)
        return True

    def ReadINSPhase(self,filename,parent=None):
        '''Read a phase from a INS file.
        '''
        Shelx = ['TITL','CELL','ZERR','LATT','SYMM','SFAC','DISP','UNIT','LAUE','EADP',
            'MORE','TIME','HKLF','OMIT','SHEL','BASF','TWIN','EXTI','SWAT',
            'HOPE','MERG','SPEC','RESI','RTAB','MPLA','HFIX','MOVE','ANIS','AFIX',
            'FRAG','FEND','EXYZ','EDAP','EQIV','CONN','PART','BIND','FREE','DFIX','DANG',
            'BUMP','SAME','SADI','CHIV','FLAT','DELU','SIMU','DEFS','ISOR','NCSY','RIGU',
            'SUMP','L.S.','CGLS','BLOC','DAMP','STIR','WGHT','FVAR','BOND','CONF','MPLA',
            'HTAB','LIST','ACTA','SIZE','TEMP','WPDB','FMAP','GRID','PLAN','MOLE']
        self.errors = 'Error opening file'
        fp = open(filename, 'r')
        Phase = {}
        Title = ''
        Atoms = []
        aTypes = []
        S = fp.readline()
        line = 1
        SGData = None
        cell = None
        SFACcontinue = True
        while S:
            if '!' in S:
                S = S.split('!')[0]
            self.errors = 'Error reading at line '+str(line)
            Atom = []
            if 'TITL' in S[:4].upper():
                Title = S[4:72].strip()
                ifTITL = True
            elif not S.strip():
                pass
            elif 'CELL' in S[:4].upper():
                cellRec = S.split()
                abc = cellRec[2:5]
                angles = cellRec[5:8]
                cell=[float(abc[0]),float(abc[1]),float(abc[2]),
                    float(angles[0]),float(angles[1]),float(angles[2])]
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                AA,AB = G2lat.cell2AB(cell)
                SGData = G2obj.P1SGData # P 1
                self.warnings += '\nThe space group is not given in an ins file and has been set to "P 1".'
                self.warnings += "\nChange this in phase's General tab; NB: it might be in the Phase name."
            elif S[:4].upper() in 'SFAC' or SFACcontinue:
                aTypes = S[4:].split()
                if 'H' in aTypes:
                    self.warnings += '\n\nHydrogen atoms found; consider replacing them with stereochemically tied ones'
                    self.warnings += '\nas Shelx constraints & HFIX commands are ignored.'
                    self.warnings += "\nDo 'Edit/Insert H atoms' in this phase's Atoms tab after deleting the old ones."
                if "=" in aTypes:
                    SFACcontinue = True
                else:
                    SFACcontinue = False
            elif S[0] == 'Q':
                pass
            elif '\x1a' in S[:4]:
                pass
            elif S[:3].upper() == 'REM':
                pass
            elif S[:3].upper() == 'END':
                pass
            elif S[:4].strip().upper() not in Shelx:   #this will find an atom record or an extended title!
                if ifTITL:
                    Title += S.strip()
                else:
                    AtRec = S.split()
                    Atype = aTypes[int(AtRec[1])-1]
                    Aname = AtRec[0]
                    Afrac = abs(float(AtRec[5]))%10.
                    x,y,z = AtRec[2:5]
                    XYZ = np.array([float(x),float(y),float(z)])
                    XYZ = np.where(np.abs(XYZ)<0.00001,0,XYZ)
                    SytSym,Mult = G2spc.SytSym(XYZ,SGData)[:2]
                    if '=' not in S:
                        IA = 'I'
                        Uiso = float(AtRec[6])
                        if Uiso < 0. or Uiso > 1.0:
                            Uiso = 0.025
                        Uij = [0. for i in range(6)]
                    else:
                        IA = 'A'
                        Uiso = 0.
                        Ustr = AtRec[6:8]
                        S = fp.readline()
                        if '!' in S:
                            S = S.split('!')[0]
                        AtRec = S.split()
                        line += 1
                        Ustr += AtRec
                        Uij = [float(Ustr[i]) for i in range(6)]
                        Uij = Uij[0:3]+[Uij[5],Uij[4],Uij[3]]
                    Atom = [Aname,Atype,'',XYZ[0],XYZ[1],XYZ[2],Afrac,SytSym,Mult,IA,Uiso]
                    Atom += Uij
                    Atom.append(ran.randint(0,sys.maxsize))
                    Atoms.append(Atom)
            S = fp.readline()
            if S[:4].strip().upper() in Shelx:
                ifTITL = False
            line += 1
        fp.close()
        self.errors = 'Error after read complete'
        Phase = G2obj.SetNewPhase(Name='ShelX phase',SGData=SGData,cell=cell+[Volume,])
        Phase['General']['Name'] = Title
        Phase['General']['Type'] = 'nuclear'
        Phase['General']['AtomPtrs'] = [3,1,7,9]
        Phase['Atoms'] = Atoms
        return Phase
