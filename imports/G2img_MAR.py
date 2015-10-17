# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_MAR: MAR image files*
--------------------------------------

Routine to read several MAR formats, .mar3450,.mar2300,.mar2560

'''

import sys
import os
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class MAR_ReaderClass(G2IO.ImportImage):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.mar3450','.mar2300','.mar2560'),
            strictExtension=True,
            formatName = 'MAR image',
            longFormatName = 'MAR Research 345, 230 and 256 image files'
            )

    def ContentsValidator(self, filepointer):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using Bob's routine
        '''
        self.Comments,self.Data,self.Npix,self.Image = G2IO.GetMAR345Data(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
