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

This shows an example of an importer that will handle files with
more than a single image. 

'''

import numpy as np
import h5py
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")

class Hdf5_Reader(G2IO.ImportImage):
    '''Routine to read a HD5 image, typically from APS Sector 6.
    B. Frosik/SDM. Initial version.
        
    '''
    dsetlist = []

    def __init__(self):
        print 'start'
        super(self.__class__,self).__init__( # fancy way to self-reference
                                             extensionlist=('.hdf5','.hd5','.h5','.hdf'),
                                             strictExtension=True,
                                             formatName = 'HDF5 image',
                                             longFormatName = 'HDF5 image file'
                                             )

    def ContentsValidator(self, filepointer):
        '''no test at this time
        '''
        try:
            # the following does not work, filepointer is not a filename
            f = h5py.File(filepointer.name, 'r')
            f.close()
            return True
        except IOError:
            return False               

    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        '''Read using HDF5 file reader, :func:`ReadData`
        '''
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        f = None
        try:
            f = h5py.File(filename, 'r')
            #print f.keys()
            self.visit(f)
            self.Comments,self.Data,self.Npix,self.Image,more = self.read_set(f,imagenum=imagenum)
            self.LoadImage(ParentFrame,filename,imagenum)
            self.repeatcount = imagenum
            self.repeat = more
            print('Read image #'+str(imagenum)+' from file '+filename)
            return True
        except IOError:
            print 'cannot open file ', filename
            return False
        finally:
            if f is not None:
               f.close()
            
    def visit(self, f):         
        def func(name, dset):
            datakeyword = 'data'
            if isinstance(dset, h5py.Dataset):
                self.dsetlist.append(dset.name)
                lastindex = dset.name.rfind(datakeyword)
                if lastindex is not -1:
                    index = len(dset.name)-len(datakeyword)
                    if index == lastindex:
                        self.image = np.array(f[dset.name],dtype=np.int32)
        f.visititems(func)
        
    def read_set(self,f,imagenum=1):
         more = False
         head = ['raw data']
         #if multiple frames the array size is 1 x number_frames x 2048 x 2048
         num_dim = len(self.image.shape)
         #GSASIIpath.IPyBreak()
         if num_dim > 2:
             image = self.image[0,imagenum-1]
             more = self.image.shape[1] > imagenum
         else:
             image = self.image
         x_dim = image.shape[0]
         y_dim = image.shape[1]
         sizexy = [x_dim,y_dim]
         Npix = sizexy[0]*sizexy[1]             
         data = {'pixelSize':[200,200],'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy}
         return head,data,Npix,image,more
