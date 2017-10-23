# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2img_png: png image file*
---------------------------------------

Routine to read an image in .png (Portable Network Graphics) format.
For now, the only known use of this is with converted Mars Rover (CheMin)
tif files, so default parameters are for that.

'''

from __future__ import division, print_function
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
class png_ReaderClass(G2obj.ImportImage):
    '''Reads standard PNG images; parameters are set to those of the
    Mars Rover (CheMin) diffractometer.
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.png',),
            strictExtension=True,
            formatName = 'PNG image',
            longFormatName = 'PNG image from CheMin'
            )

    def ContentsValidator(self, filename):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Reads using standard scipy PNG reader
        '''
        import scipy.misc
        self.Image = scipy.misc.imread(filename,flatten=True)
        self.Npix = self.Image.size
        if self.Npix == 0:
            return False
        if ParentFrame:
            self.SciPy = True
            self.Comments = ['no metadata']
            pixy = list(self.Image.shape)
            sizexy = [40.,40.]
            self.Data = {'wavelength': 1.78892, 'pixelSize': sizexy, 'distance': 18.0,'size':pixy}
            self.Data['center'] = [pixy[0]*sizexy[0]/1000.,pixy[1]*sizexy[1]/2000.]
        self.LoadImage(ParentFrame,filename)
        return True
