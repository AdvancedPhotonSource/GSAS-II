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

import h5py
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")

class Hdf5_Reader(G2IO.ImportImage):
    '''Routine to read a HD5 image, typically from APS Sector 6.
    B. Frosik/SDM. Initial version.
                   Refactored.        
    '''
    dsetlist = []
    buffer = {}
    init = False

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
        #self.buffer = kwarg.get('buffer')
        
        if self.init is False:
            self.init = True
            self.readDataset(filename)

        len = self.buffer['len']
        if len is 1:
            image = self.buffer['image']
        else:
            images = self.buffer['image']
            image = images[imagenum-1]

        self.Comments,self.Data,self.Npix,self.Image = self.read_set(image,imagenum=imagenum)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        self.repeatcount = imagenum 
        self.repeat = imagenum < len
        #print('Read image #'+str(imagenum)+' from file '+filename)
        return True
        
    def readDataset(self,filename):
            imagenum = 1
            try:
                f = h5py.File(filename, 'r')
                self.visit(f, imagenum)
            except IOError:
                print 'cannot open file ', filename
                return False
            finally:
                if f is not None:
                   f.close()
        
         
    def visit(self, f, imagenum):
        #GSASIIpath.IPyBreak()   
        print 'visiting'      
        def func(name, dset):
            datakeyword = 'data'
            if isinstance(dset, h5py.Dataset):
                self.dsetlist.append(dset.name)
                lastindex = dset.name.rfind(datakeyword)
                if lastindex is not -1:
                    index = len(dset.name)-len(datakeyword)
                    if index == lastindex:
                        if self.buffer is not None: 
                            if len(dset.shape) > 2:
                                    image = dset[0,...]
                                    self.buffer['image'] = dset[0,...]
                                    self.buffer['len'] = image.shape[0]
                            else:
                                self.buffer['image'] = dset[()]
                                self.buffer['len'] = 1
                            self.buffer['name'] = dset.name
                        #GSASIIpath.IPyBreak()
        f.visititems(func)
        
    def read_set(self,image,imagenum=1):
         #more = False
         head = ['raw data']
         x_dim = image.shape[0]
         y_dim = image.shape[1]
         sizexy = [x_dim,y_dim]
         Npix = sizexy[0]*sizexy[1]             
         data = {'pixelSize':[200,200],'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy}
         return head,data,Npix,image
