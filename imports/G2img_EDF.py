# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_EDF: .edf image file*
--------------------------------------

'''

import sys
import os
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class EDF_ReaderClass(G2IO.ImportImage):
    '''Routine to read a Read European detector data .edf file.
    This is a particularly nice standard. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.edf',),
            strictExtension=True,
            formatName = 'EDF image',
            longFormatName = 'European Data Format image file'
            )

    def ContentsValidator(self, filepointer):        
        '''no test used at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using Bob's routine :func:`GSASIIIO.GetEdfData`
        (to be moved to this file, eventually)
        '''
        self.Comments,self.Data,self.Npix,self.Image = G2IO.GetEdfData(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
