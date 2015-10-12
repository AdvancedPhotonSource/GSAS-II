# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_png: png image file*
---------------------------------------

Routine to read an image in .png (Portable Network Graphics) format.
For now, the only known use of this is with converted CheMin tif files
so default parameters are that machine.

'''

import sys
import os
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class png_ReaderClass(G2IO.ImportImage):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.png',),
            strictExtension=True,
            formatName = 'PNG image',
            longFormatName = 'PNG image from CheMin'
            )

    def ContentsValidator(self, filepointer):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using scipy PNG reader
        '''
        import scipy.misc
        self.Image = scipy.misc.imread(filename,flatten=True)
        self.Npix = self.Image.size
        self.Comments = ['no metadata']
        pixy = list(self.Image.shape)
        sizexy = [40,40]
        self.Data = {'wavelength': 1.78892, 'pixelSize': sizexy, 'distance': 18.0,'size':pixy}
        self.Data['center'] = [pixy[0]*sizexy[0]/1000,pixy[1]*sizexy[1]/2000]
        if self.Npix == 0 or not self.Comments:
            return False
        if ParentFrame:
            G2IO.EditImageParms(ParentFrame,self.Data,self.Comments,self.Image,filename)
        return True
# N.B. This replaces G2IO.GetPNGData
