# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*Module G2img_HDF5: summed HDF5 image file*
---------------------------------------

Reads all images found in a HDF5 file.

'''

try:
    import h5py
except ImportError:
    h5py = None
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")

class HDF5_Reader(G2IO.ImportImage):
    '''Routine to read a HD5 image, typically from APS Sector 6.
    B. Frosik/SDM. Initial version.
                   Refactored.        
    '''
    dsetlist = []
    buffer = {}
    init = False

    def __init__(self):
        if h5py is None:
            self.UseReader = False
            print('HDF5 Reader skipped because h5py library is not installed')
        super(self.__class__,self).__init__( # fancy way to self-reference
                                             extensionlist=('.hdf5','.hd5','.h5','.hdf'),
                                             strictExtension=True,
                                             formatName = 'HDF5 image',
                                             longFormatName = 'HDF5 image file'
                                             )

    def ContentsValidator(self, filepointer):
        '''Test if valid by seeing if the HDF5 library recognizes the file.
        '''
        try:
            f = h5py.File(filepointer.name, 'r')
            f.close()
            return True
        except IOError:
            return False               

    def Reader(self, filename, filepointer, ParentFrame=None, **kwarg):
        '''Scan file structure using :meth:`visit` and map out locations of image(s)
        then read one image using :meth:`readDataset`. Save map of file structure in
        buffer arg, if used. 
        '''
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        self.buffer = kwarg.get('buffer',{})
        try:
            fp = h5py.File(filename, 'r')
            if not self.buffer.get('init'):
                self.buffer['init'] = True
                self.visit(fp)
            self.Comments,self.Data,self.Npix,self.Image = self.readDataset(fp,imagenum)
            if self.Npix == 0 or not self.Comments:
                return False
            self.LoadImage(ParentFrame,filename,imagenum)
            self.repeatcount = imagenum 
            self.repeat = imagenum < len(self.buffer['imagemap'])
            if GSASIIpath.GetConfigValue('debug'): print('Read image #'+str(imagenum)+' from file '+filename)
            return True
        except IOError:
            print 'cannot open file ', filename
            return False
        finally:
            fp.close()

    def visit(self, fp):
        '''Recursively visit each node in an HDF5 file. For nodes
        ending in 'data' look at dimensions of contents. If the shape is
        length 2 or 4 assume an image and index in self.buffer['imagemap']
        ''' 
        datakeyword = 'data'
        def func(name, dset):
            if not hasattr(dset,'shape'): return # not array, can't be image
            if isinstance(dset, h5py.Dataset):
                if dset.name.endswith(datakeyword):
                    dims = dset.shape
                    if len(dims) == 4:
                        self.buffer['imagemap'] += [
                            (dset.name,i) for i in range(dims[1])]
                    elif len(dims) == 2:
                        self.buffer['imagemap'] += [(dset.name,None)]
        if GSASIIpath.GetConfigValue('debug'): print 'visit'
        self.buffer['imagemap'] = []
        fp.visititems(func)
        
    def readDataset(self,fp,imagenum=1):
        '''Read a specified image number from a file
        '''
        name,num = self.buffer['imagemap'][imagenum-1] # look up in map
        dset = fp[name]
        if num is None:
            image = dset[()]
        else:
            image = dset[0,num,...]
        head = ['raw data']
        sizexy = list(image.shape) 
        Npix = sizexy[0]*sizexy[1]             
        data = {'pixelSize':[200,200],'wavelength':0.15,'distance':250.0,
                'center':[204.8,204.8],'size':sizexy}
        return head,data,Npix,image
