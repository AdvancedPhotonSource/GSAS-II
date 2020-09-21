# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
#
'''
*Module G2phase_xyz: read coordinates from an xyz file*
-------------------------------------------------------

A short routine to read in a phase from an xyz Cartesian coordinate file

'''

from __future__ import division, print_function
import sys
import os.path
import math
import random as ran
import numpy as np
import GSASIIobj as G2obj
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")

class XYZ_ReaderClass(G2obj.ImportPhase):
    'Routine to import Phase information from a XYZ file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.xyz','.XYZ'),
            strictExtension=True,
            formatName = 'XYZ',
            longFormatName = 'XYZ Cartesian coordinate file import'
            )
    def ContentsValidator(self, filename):
        '''Taking a stab a validating: 1st line should be a number
        '''
        fp = open(filename,'r')
        try:
            int(fp.readline().strip())
        except:
            fp.close()
            return False
        fp.close()
        return True

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read a PDF file using :meth:`ReadPDBPhase`'
        self.errors = 'Error opening file'
        fp = open(filename, 'Ur')
        self.Phase = {}
        natom = int(fp.readline().strip())
        Title = os.path.basename(filename)
        skip = fp.readline()
        line = 2
        SGData = G2obj.P1SGData # P 1
        self.warnings += '\nNo space group in file, set to "P 1".'
        self.warnings += "Change this in phase's General tab."
        cell = [10.,10.,10.,90.,90.,90.]
        Volume = G2lat.calc_V(G2lat.cell2A(cell))

        counts = {}
        Atoms = []
        for i in range(natom):
            line += 1
            self.errors = 'Error reading at line '+str(line)
            l = fp.readline()
            Type = l.split()[0]
            XYZ = [float(x)/10. for x in l.split()[1:4]]
            if Type not in counts:
                counts[Type] = 0
            counts[Type] += 1
            Aname = Type + str(counts[Type])
            Atoms.append([Aname,Type.strip().capitalize(),'',XYZ[0],XYZ[1],XYZ[2],
                    1.0,'1',1,'I',0.0,0,0,0,0,0,0,ran.randint(0,sys.maxsize)])
        fp.close()
        self.errors = 'Error after read complete'
        self.Phase = G2obj.SetNewPhase(Name=Title,SGData=SGData,cell=cell+[Volume,])
        self.Phase['General']['Type'] = 'nuclear'
        self.Phase['General']['AtomPtrs'] = [3,1,7,9]    
        self.Phase['Atoms'] = Atoms
        return True
