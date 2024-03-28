# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2phase_GPX.py $
# $Id: G2phase_GPX.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''
from __future__ import division, print_function
import platform
import sys
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import random as ran
import GSASIIobj as G2obj
import GSASIIstrIO as G2stIO
import GSASIIpath
try:
    import GSASIIctrlGUI as G2G
except ImportError:
    pass
GSASIIpath.SetVersionNumber("$Revision: 5577 $")

class PhaseReaderClass(G2obj.ImportPhase):
    'Opens a .GPX file and pulls out a selected phase'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.gpx',),
            strictExtension=True,
            formatName = 'GSAS-II gpx',
            longFormatName = 'GSAS-II project (.gpx file) import'
            )
        
    def ContentsValidator(self, filename):
        "Test if the 1st section can be read as a cPickle block, if not it can't be .GPX!"
        if True:
            fp = open(filename,'rb')
        try: 
            if '2' in platform.python_version_tuple()[0]:
                data = cPickle.load(fp)
            else:
                data = cPickle.load(fp,encoding='latin-1')
        except:
            self.errors = 'This is not a valid .GPX file. Not recognized by cPickle'
            fp.close()
            return False
        fp.close()
        return True

    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read a phase from a .GPX file. Does not (yet?) support selecting and reading
        more than one phase at a time.'''
        try:
            phasenames = G2stIO.GetPhaseNames(filename)
        except:
            self.errors = 'Reading of phase names failed'
            return False
        if not phasenames:
            self.errors = 'No phases found in '+str(filename)
            return False            # no blocks with coordinates
        elif len(phasenames) == 1: # one block, no choices
            selblk = 0
        else:                       # choose from options                
            selblk = G2G.PhaseSelector(phasenames,ParentFrame=ParentFrame,
                title= 'Select a phase from the list below',)
            if selblk is None:
                self.errors = 'No phase selected'
                return False # User pressed cancel
        self.Phase = G2stIO.GetAllPhaseData(filename,phasenames[selblk])
        self.Phase['Histograms'] = {}       #remove any histograms
        self.Phase['Pawley ref'] = []       # & any Pawley refl.
        self.Phase['RBModels'] = {}
        self.Phase['Drawing'] = {}
        if 'MCSA' in self.Phase:
            del self.Phase['MCSA']
        if 'Map Peaks' in self.Phase:
            del self.Phase['Map Peaks']
        if 'Map' in self.Phase['General']:
            del self.Phase['General']['Map']
        self.Phase['ranId'] = ran.randint(0,sys.maxsize)
        return True
