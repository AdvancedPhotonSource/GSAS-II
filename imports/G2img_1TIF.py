# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_TIF: .tif image file*
--------------------------------------

Routine to read an image in Tagged-image file (TIF) format as well as a variety
of slightly incorrect pseudo-TIF formats used at instruments around the world. 

'''

import sys
import struct as st
import os.path as ospath
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class TIF_ReaderClass(G2IO.ImportImage):
    '''Routine to read an image in Tagged-image file (TIF) format as well as a variety
    of slightly incorrect pseudo-TIF formats
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.tif','.tiff'),
            strictExtension=False,
            formatName = 'TIF image',
            longFormatName = 'Various .tif and pseudo-TIF formats'
            )

    def ContentsValidator(self, filepointer):
        '''Does the header match the required TIF header?
        '''
        tag = filepointer.read(2)
        if tag == 'II' and int(st.unpack('<h',filepointer.read(2))[0]) == 42: #little endian
            pass
        elif tag == 'MM' and int(st.unpack('>h',filepointer.read(2))[0]) == 42: #big endian
            pass
        else:
            return False # header not found; not valid TIF
        return True
    
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read the TIF file using Bob's routine
        '''
        
        self.Comments,self.Data,self.Npix,self.Image = G2IO.GetTifData(filename)
        if self.Npix == 0:
            print("GetTifData failed to read "+str(filename)+" Trying SciPy")
            import scipy.misc
            self.Image = scipy.misc.imread(filename,flatten=True)
            self.Npix = self.Image.size
            if ParentFrame:
                self.Comments = ['no metadata']
                self.Data = {'wavelength': 0.1, 'pixelSize': [200, 200], 'distance': 100.0}
                self.Data['size'] = list(self.Image.shape)
                self.Data['center'] = [int(i/2) for i in self.Image.shape]
                G2IO.EditImageParms(ParentFrame,self.Data,self.Comments,self.Image,filename)
        if self.Npix == 0:
            return False
        self.LoadImage(ParentFrame,filename)
        return True

