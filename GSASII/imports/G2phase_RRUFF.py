'''Reads structure files downloaded from the RRUFF database either as ASCII
text files or .rtf files that somehow are generated
'''
# import sys
# import numpy as np
# import random as ran
# #from .. import GSASIIobj as G2obj

# class RRUFFReader(G2obj.ImportPhase):
#     '''A fairly quickly-written importer to pull out the phase info from a 
#     RRUFF database (https://rruff.info) text file. 
#     '''
#     pass
    # def __init__(self):
    #     super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
    #         extensionlist=('.txt'),
    #         strictExtension=True,
    #         formatName = 'RRUFF DIF file',
    #         longFormatName = 'RRUFF database DIF file download (*.txt)'
    #         )
        
    # def ContentsValidator(self, filename):
    #     "Test if the file has CELL PARAMETERS: and SPACE GROUP: entries"
    #     with open(filename,'r') as fp:
    #         txt = fp.read()
    #     status = True
    #     for s in ('CELL PARAMETERS:','SPACE GROUP:'):
    #         if not s in txt:
    #             print(s,'not found')
    #             if self.errors is None: self.errors = ''
    #             self.errors += f'no {s} record found; not valid. '
    #             status = False
    #     if not ' ATOM ' in txt:
    #         print('ATOM line not found')
    #         if self.errors is None: self.errors = ''
    #         self.errors += 'no ATOM record found; no structure. '
    #         status = False           
    #     fp.close()
    #     return status

    # def Reader(self,filename,ParentFrame=None, **unused):
    #     'Read a DIF file from RRUFF'
    #     #self.errors = 'Error opening file'
    #     from .. import GSASIIspc as G2spc
    #     from .. import GSASIIlattice as G2lat
    #     Title = ''
    #     atomsmode = False
    #     fp = open(filename, 'r')
    #     SGData = None
    #     Atoms = []
    #     cell = None
    #     for line in fp.readlines():
    #         if len(line) == 0: continue
    #         line = line.strip()
    #         if atomsmode:
    #             # all DIF files appear to have atoms as element, X, Y, Z, OCCUPANCY, ISO(B)
    #             if r'\par' in line: # some files seem to have some formatting
    #                 line = line.split('\\par')[0]
    #             try:
    #                 Atype = line.split()[0]
    #                 x,y,z,Afrac,B = [float(i) for i in line.split()[1:]]
    #                 Uiso = B/(8*np.pi**2)
    #                 XYZ = np.array([float(x),float(y),float(z)])
    #                 XYZ = np.where(np.abs(XYZ)<0.00001,0,XYZ)
    #                 SytSym,Mult = G2spc.SytSym(XYZ,SGData)[:2]
    #                 IA = 'I'
    #                 i = 1
    #                 while f"{Atype}{i}" in [i[0] for i in Atoms]:
    #                     i += 1
    #                     if i > 999:
    #                         Aname = f"{Atype}?"
    #                         break
    #                 else:
    #                     Aname = f"{Atype}{i}"
    #                 if Atype.upper() in ['WA','OH','OW','OA','OB','OC','OD','OE','OL','OP','OO']:
    #                     Atype = 'O'
    #                 Atom = [Aname,Atype,'',x,y,z,Afrac,SytSym,Mult,IA,Uiso]
    #                 Atom += 6*[0]
    #                 Atom.append(ran.randint(0,sys.maxsize))
    #                 Atoms.append(Atom)
    #             except:
    #                 atomsmode = False
    #                 break
    #         if Title == '':
    #             if '\\rtf' in line: continue
    #             Title = line
    #             if r'\par' in line: # some files seem to have some formatting
    #                 Title = line.split()[-1].split('\\')[0]
    #         if 'CELL PARAMETERS:' in line:
    #             L = line.split(':')[1]
    #             if r'\par' in line: # some files seem to have some formatting
    #                 L = L.split('\\par')[0]
    #             cellRec = L.split()
    #             abc = cellRec[0:3]
    #             angles = cellRec[3:]
    #             cell=[float(abc[0]),float(abc[1]),float(abc[2]),
    #                 float(angles[0]),float(angles[1]),float(angles[2])]
    #             Volume = float(G2lat.calc_V(G2lat.cell2A(cell)))
    #         elif 'SPACE GROUP:' in line:
    #             SGData = G2obj.P1SGData # P 1
    #             S = line.split(':')[1].strip()
    #             if r'\par' in S: # some files seem to have some formatting
    #                 S = S.split()[0]
    #             Sinit = S
    #             E,SGData = G2spc.SpcGroup(S)
    #             if E:
    #                 SpGrpNorm = G2spc.StandardizeSpcName(S)
    #                 if SpGrpNorm:
    #                     E,SGData = G2spc.SpcGroup(SpGrpNorm)
    #             if E and '_' in S:
    #                 S = S.replace('_','')
    #             if E:
    #                 E,SGData = G2spc.SpcGroup(S[0]+' '+S[1:])
    #             # some space group names that show up in RRUFF
    #             SGfixDict = {'B2/b': 'B 2/b 1 1', 'Bb21m':'B b 21 m',
    #                              'P21/b': 'P 21/b 1 1', 'P21212': 'P 21 21 2',
    #                              'P21ca': 'P 21 c a', 'P21cn': 'P 21 c n',
    #                              'P21nb': 'P21 n b', 'P21nm': 'P21 n m',
    #                              'P63cm': 'P63 c m', 'P63mc': 'P63 m c',
    #                              'Pn21a': 'P n 21 a', 'Pn21m': 'P n 21 m'}

    #             if E and S in SGfixDict:
    #                 E,SGData = G2spc.SpcGroup(SGfixDict[S])
    #             if E:
    #                 self.warnings += f'ERROR in space group symbol {Sinit!r}'
    #                 self.warnings += '\nThe space group has been set to "P 1". '
    #                 self.warnings += "Change this in phase's General tab."
    #                 self.warnings += " Error msg="+G2spc.SGErrors(E)
    #             elif SGData['SpGrp'] in G2spc.spg2origins:
    #                 self.warnings += f"WARNING space group {SGData['SpGrp']} has two Origins"
    #         elif 'ATOM ' in line:
    #             atomsmode = True
    #     if SGData is None:
    #         if self.errors is None: self.errors = ''
    #         self.errors += 'no Space Group record found; not valid. '
    #         return False
    #     if cell is None:
    #         if self.errors is None: self.errors = ''
    #         self.errors += 'no CELL record found; not valid. '
    #         return False
    #     if len(Atoms) == 0:
    #         if self.errors is None: self.errors = ''
    #         self.errors += 'no Atoms found; not valid. '
    #         return False
    #     if cell[3] == cell[4] == cell[5] and cell[3] != 90:
    #         self.warnings += f'Note: {filename!r} has rhombohedral cell, changing space group'
    #         E,SGData = G2spc.SpcGroup(SGData['SpGrp'] + ' R')
    #     self.MPhase = None
    #     self.Phase = G2obj.SetNewPhase(Name=Title,SGData=SGData,cell=cell+[Volume,])
    #     self.Phase['General']['Name'] = Title
    #     self.Phase['General']['Type'] = 'nuclear'
    #     self.Phase['General']['AtomPtrs'] = [3,1,7,9]
    #     self.Phase['Atoms'] = Atoms
    #     fp.close()
    #     return True
