# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
from .. import GSASIIobj as G2obj
from .. import GSASIIpath
from .. import GSASIIfiles as G2fil
try:
    import imageio
except ImportError:
    imageio = None
class png_ReaderClass(G2obj.ImportImage):
    '''Reads standard PNG images; parameters are set to those of the
    Mars Rover (CheMin) diffractometer.
    '''
    def __init__(self):
        if imageio is None:
            self.UseReader = False
            msg = 'CheMin Reader skipped because imageio library is not installed'
            if GSASIIpath.condaTest():
                msg += ' To fix this use command:\n\tconda install imageio'
            G2fil.ImportErrorMsg(msg,{'CheMin image importer':['imageio']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.png',),
            strictExtension=True,
            formatName = 'CheMin PNG image',
            longFormatName = 'PNG image from CheMin'
            )

    def ContentsValidator(self, filename):
        '''no test at this time
        '''
        return True

    def Reader(self,filename, ParentFrame=None, **unused):
        '''Reads using standard scipy PNG reader
        '''
        self.Image = imageio.imread(filename,flatten=True)
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
