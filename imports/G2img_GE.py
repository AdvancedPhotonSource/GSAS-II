# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_GE: summed GE image file*
---------------------------------------

Routine to read a summed GE image from APS Sector 1

'''

import sys
import os
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class GEsum_ReaderClass(G2IO.ImportImage):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.sum','.cor','.avg','.ge','.ge1','.ge2','.ge3','.ge4'),
            strictExtension=True,
            formatName = 'GE image',
            longFormatName = 'Summed GE image file'
            )

    def ContentsValidator(self, filepointer):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using Bob's routine
        '''
        self.Comments,self.Data,self.Npix,self.Image = G2IO.GetGEsumData(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        return True
