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

class ExportPhasePDB(G2IO.ExportBaseclass):
    '''Used to create a PDB file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'PDB',
            extension='.PDB',
            longFormatName = 'Export phase as .PDB file'
            )
        self.exporttype = ['phase']
        self.multiple = True

    def Exporter(self,event=None):
        '''Export as a PDB file
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
            General = phasedict['General']
            if General['Type'] != 'macromolecular':
                print 'not macromolecular phase'
                return
            print phasedict.keys()
            print General['SGData'].keys()           
            i = self.Phases[phasenam]['pId']
            if len(self.phasenam) > 1: # if more than one filename is included, add a phase #
                nam,ext = os.path.splitext(self.filename)
                fil = nam+"_"+str(i)+ext
            else:
                fil = self.filename
            fp = self.OpenFile(fil)
            self.Write("HEADER phase "+str(phasenam)+" from "+str(self.G2frame.GSASprojectfile))
            # get cell parameters 
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            line = "CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} ".format(*G2lat.A2cell(A))
            line += General['SGData']['SpGrp'].ljust(13)
            line += '%2d'%(len(General['SGData']['SGOps'])*len(General['SGData']['SGCen']))
            self.Write(line)
            self.Write('ORIGX1      1.000000  0.000000  0.000000        0.00000')
            self.Write('ORIGX2      0.000000  1.000000  0.000000        0.00000')
            self.Write('ORIGX3      0.000000  0.000000  1.000000        0.00000')
            G = G2lat.A2Gmat(A)[0]
            A,B = G2lat.Gmat2AB(G)
            self.Write('SCALE1     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[0]))
            self.Write('SCALE2     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[1]))
            self.Write('SCALE3     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[2]))
            Atoms = phasedict['Atoms']
           

            self.CloseFile()
            print('Phase '+str(phasenam)+' written to file '+str(fil))
