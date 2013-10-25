# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2phase_GPX: Import phase from GSAS-II project*
--------------------------------------------------------

Copies a phase from another GSAS-II project file into the
current project.

'''
# Routines to import Phase information from GSAS-II .gpx files
import cPickle
import GSASIIIO as G2IO
import GSASIIstrIO as G2stIO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

class PhaseReaderClass(G2IO.ImportPhase):
    'Opens a .GPX file and pulls out a selected phase'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.gpx',),
            strictExtension=True,
            formatName = 'GSAS-II gpx',
            longFormatName = 'GSAS-II project (.gpx file) import'
            )
    def ContentsValidator(self, filepointer):
        # if the 1st section can't be read as a cPickle file, it can't be!
        try: 
            cPickle.load(filepointer)
        except:
            return False
        return True
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        try:
            phasenames = G2stIO.GetPhaseNames(filename)
        except:
            return False
        if not phasenames:
            return False            # no blocks with coordinates
        elif len(phasenames) == 1: # no choices
            selblk = 0
        else:                       # choose from options                
            selblk = self.PhaseSelector(
                phasenames,
                ParentFrame=ParentFrame,
                title= 'Select a phase from the list below',
                )
            if selblk is None: return False # User pressed cancel
        try:
            self.Phase = G2stIO.GetAllPhaseData(filename,phasenames[selblk])
            self.Phase['Histograms'] = {}       #remove any histograms
            self.Phase['Pawley ref'] = []       # & any Pawley refl.
            return True
        except Exception as detail:
            import sys
            print self.formatName+' error:',detail # for testing
            print sys.exc_info()[0] # for testing
            import traceback
            print traceback.format_exc()
