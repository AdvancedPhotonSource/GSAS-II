# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_Rigaku: .stl image file*
--------------------------------------

Routine to read a Rigaku R-Axis IV image file

'''

import sys
import os
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class Rigaku_ReaderClass(G2IO.ImportImage):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.stl',),
            strictExtension=True,
            formatName = 'Rigaku image',
            longFormatName = 'Read Rigaku R-Axis IV image file'
            )

    def ContentsValidator(self, filepointer):
        
        '''Does the file size make sense?
        '''
        fileSize = os.stat(fp.name).st_size
        Npix = (fileSize-6000)/2
        if Npix == 9000000 or Npix == 2250000 or Npix == 36000000:
            return True
        return False # not valid size
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using Bob's routine
        '''

        self.Comments,self.Data,self.Npix,self.Image = G2IO.GetRigaku(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        return True
