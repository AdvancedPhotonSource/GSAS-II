#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: $
# $Revision: -1 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/exports/G2export_CIF.py $
# $Id: G2export_CIF.py -1   $
########### SVN repository information ###################
'''Code to demonstrate how export routines are created
'''
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: -1 $")
import GSASIIIO as G2IO
#import GSASIIgrid as G2gd
import GSASIIstrIO as G2stIO
#import GSASIImath as G2mth
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
#import GSASIIphsGUI as G2pg
#import GSASIIstrMain as G2stMn

class ExportPhaseShelx(G2IO.ExportBaseclass):
    '''Used to create a SHELX .ins file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'SHELX',
            extension='.ins',
            longFormatName = 'Export phase as SHELX .ins file'
            )
        self.exporttype = ['phase']
        self.multiple = True

    def Exporter(self,event=None):
        '''Export as a SHELX .ins file
        '''
        # the export process starts here
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.SetupExport(event,                         # set export parameters
                            AskFile=True
                            ): return 
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            i = self.Phases[phasenam]['pId']
            if len(self.phasenam) > 1: # if more than one filename is included, add a phase #
                nam,ext = os.path.splitext(self.filename)
                fil = nam+"_"+str(i)+ext
            else:
                fil = self.filename
            fp = self.OpenFile(fil)
            self.Write("TITL from "+str(self.G2frame.GSASprojectfile)+" phase "+str(phasenam))
            # get cell parameters 
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            #cellSig = G2stIO.getCellEsd(pfx,
            #                           phasedict['General']['SGData'],A,
            #                           self.OverallParms['Covariance'])  # returns 7 vals, includes sigVol
            #cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            # write out cell parameters with dummy (0.5A) wavelength
            self.Write("CELL 0.5 {:.5f} {:.5f} {:.5f} {:.3f} {:.3f} {:.3f}".format(*G2lat.A2cell(A)))
            

            self.CloseFile()
            print('Phase '+str(phasenam)+' written to file '+str(fil))
