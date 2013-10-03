# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
# a set of short routines to read in phases from various file formats
# 
import sys
import traceback
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

class PDB_ReaderClass(G2IO.ImportPhase):
    'Routines to import Phase information from a PDB file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.pdb','.ent','.PDB','.ENT'),
            strictExtension=True,
            formatName = 'PDB',
            longFormatName = 'Original Protein Data Bank (.pdb file) import'
            )
    # I don't know enough to validate the contents
    #def ContentsValidator(self, filepointer):
    #    filepointer.seek(0) # rewind the file pointer
    #    return True

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        try:
            self.Phase = G2IO.ReadPDBPhase(filename)
            return True
        except Exception as detail:
            print 'PDB read error:',detail # for testing
            traceback.print_exc(file=sys.stdout)
            return False

class EXP_ReaderClass(G2IO.ImportPhase):
    ' Routines to import Phase information from GSAS .EXP files'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.EXP','.exp'),
            strictExtension=True,
            formatName = 'GSAS .EXP',
            longFormatName = 'GSAS Experiment (.EXP file) import'
            )
    def ContentsValidator(self, filepointer):
        filepointer.seek(0) # rewind the file pointer
        # first 13 characters should be VERSION tag -- I think
        try:
            if filepointer.read(13) == '     VERSION ':
                return True
        except: pass
        return False

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        try:
            self.Phase = G2IO.ReadEXPPhase(ParentFrame, filename)
            return True
        except Exception as detail:
            print 'GSAS .EXP read error:',detail # for testing
            traceback.print_exc(file=sys.stdout)
            return False
