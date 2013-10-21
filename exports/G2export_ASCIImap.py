#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''Code to demonstrate how export routines are created: Export a Fourier or
Charge-flip map. 
'''
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
#import GSASIIgrid as G2gd
#import GSASIIstrIO as G2stIO
#import GSASIImath as G2mth
#import GSASIIlattice as G2lat
#import GSASIIspc as G2spc
#import GSASIIphsGUI as G2pg
#import GSASIIstrMain as G2stMn

class ExportMapASCII(G2IO.ExportBaseclass):
    '''Used to create a text file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'FOX/DrawXTL file',
            extension='.grd',
            longFormatName = 'Export map as text (.grd) file'
            )
        self.exporttype = ['map']
        self.multiple = False 

    def Exporter(self,event=None):
        '''Export a map as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile=True     # prompt the user for a file name
            ): return 
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            rho = phasedict['General']['Map'].get('rho',[])
            if not len(rho):
                return
            self.OpenFile(self.filename)
            self.Write("Map of Phase "+str(phasenam)+" from "+str(self.G2frame.GSASprojectfile))
            # get cell parameters & print them
            cellList,cellSig = self.GetCell(phasenam)
            fmt = 3*" {:9.5f}" + 3*" {:9.3f}"
            self.Write(fmt.format(*cellList[:6]))
            nx,ny,nz = rho.shape
            self.Write(" {:3d} {:3d} {:3d}".format(nx,ny,nz))
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        self.Write(str(rho[i,j,k]))
            print('map from Phase '+str(phasenam)+' written to file '+str(self.filename))
            self.CloseFile()
            return
