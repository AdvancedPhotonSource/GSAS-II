# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_G2: Python pickled image*
---------------------------------------

Routine to read an image that has been pickled in Python. Images
in this format are created by the "Sum image data" command. 

'''

import sys
import os
import cPickle
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class G2_ReaderClass(G2IO.ImportImage):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.G2img',),
            strictExtension=True,
            formatName = 'Summed GSAS-II image',
            longFormatName = 'Pickled image from GSAS-II "Sum image data" command'
            )

    def ContentsValidator(self, filepointer):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using scipy PNG reader
        '''
        import scipy.misc
        Fp = open(filename,'rb')
        self.Comments,self.Data,self.Npix,self.Image = cPickle.load(Fp)
        Fp.close()
        self.LoadImage(ParentFrame,filename)
        return True
