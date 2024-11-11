# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
import GSASIIobj as G2obj
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
            self.Data = {'wavelength': 1.78892, 'pixelSize': sizexy, 'distance': 18.0,'size':pixy,'det2theta':0.0}
            self.Data['center'] = [pixy[0]*sizexy[0]/1000.,pixy[1]*sizexy[1]/2000.]
        self.LoadImage(ParentFrame,filename)
        return True
