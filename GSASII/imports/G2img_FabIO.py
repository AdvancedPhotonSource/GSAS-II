# -*- coding: utf-8 -*-
'''Importer that interfaces with the FabIO image handling library.
'''

from __future__ import division, print_function
try:
    import fabio
except ImportError:
    fabio = None
import GSASIIobj as G2obj
import GSASIIfiles as G2fil
import GSASIIpath

class FabIO_Reader(G2obj.ImportImage):
    '''Read an image using the FabIO package
    '''
    def __init__(self):
        if fabio is None:
            self.UseReader = False
            msg = 'FabIO Reader skipped because FabIO module is not installed'
            if GSASIIpath.condaTest():
                msg += ' To fix this use command:\n\tconda install h5py hdf5'
            G2fil.ImportErrorMsg(msg,{'FabIO image importer':['fabio']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hdf5','.hd5','.h5','.hdf'),              # this is a list of likely file extensions
            strictExtension=False,                                    # False allows files to be opened that have other extensions
            formatName = 'FabIO image',
            longFormatName = 'Use FabIO to import many image file types'
            )

    def ContentsValidator(self, filename):
        '''This method is optional. Use if there is a quick test that can be made to see
        if the selected file is appropriate for the FabIO library. 
        '''
        # try:
        #     fp = h5py.File(filename, 'r')
        #     fp.close()
        #     return True
        # except IOError:
        #     return False               
        return True

    def Reader(self, filename, ParentFrame=None, **kwarg):
        '''Read file 
        then read one image using :meth:`readDataset`. Save map of file structure in
        buffer arg, if used. 
        '''
        imagenum = kwarg.get('blocknum')   # use this for files that have more than one image
        if imagenum is None: imagenum = 1
        self.buffer = kwarg.get('buffer',{})   # use this to save information that is costly to
        # reproduce when reading multiple images from a file

        try:
            self.Comments = ''   # a list of lines of text info about image (more is better)
            self.Data = {'pixelSize':[200.,200.],    # size of pixel in x & y direction
                          'wavelength':0.123,      # in A
                          'distance':100.,         # actual sample to detector distance in cm (might be fit)
                          'center':[204.8, 204.8], # measured from top left corner of the detector (mm)
                          'size':[2048, 2048],     # pixels on x & y
                          'setdist':100.,          # nominal sample to detector distance, as reported by instrument
                          'PolaVal':[0.99,False]}  # False is refinement flag
            self.Image = # should be np.array dtype=int32
            self.Npix = self.Image.size # number of pixels
            self.repeatcount = imagenum
            self.repeat = False # set to True if there are more images in current file to be read
            return True   # indicates success on import
        except IOError:
            print (f'cannot open/read file {filename}')
            return False   # indicates failure to import
        finally:
            fp.close()
