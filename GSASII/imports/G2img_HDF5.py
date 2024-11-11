# -*- coding: utf-8 -*-
'''A reader for HDF-5 files. This should be as generic as possible, but 
at present this is pretty much customized for XSD-MPE (APS) uses. 
'''

from __future__ import division, print_function
try:
    import h5py
except ImportError:
    h5py = None
import GSASIIobj as G2obj
import GSASIIfiles as G2fil
import GSASIIpath

class HDF5_Reader(G2obj.ImportImage):
    '''Routine to read a HDF-5 image, typically from APS Sector 6.
    initial version from B. Frosik/SDM. Updated to also handle Varex images
    from APS 1-ID-C.
    '''
    dsetlist = []
    buffer = {}
    init = False

    def __init__(self):
        if h5py is None:
            self.UseReader = False
            msg = 'HDF5 Reader skipped because h5py library is not installed'
            if GSASIIpath.condaTest():
                msg += ' To fix this use command:\n\tconda install h5py hdf5'
            G2fil.ImportErrorMsg(msg,{'HDF5 image importer':['h5py','hdf5']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hdf5','.hd5','.h5','.hdf'),strictExtension=True,
            formatName = 'HDF5 image',longFormatName = 'HDF5 image file')

    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the HDF5 library recognizes the file.
        '''
        try:
            fp = h5py.File(filename, 'r')
            fp.close()
            return True
        except IOError:
            return False               

    def Reader(self, filename, ParentFrame=None, **kwarg):
        '''Scan file structure using :meth:`visit` and map out locations of image(s)
        then read one image using :meth:`readDataset`. Save map of file structure in
        buffer arg, if used. 
        '''
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        self.buffer = kwarg.get('buffer',{})
        try:
            fp = h5py.File(filename, 'r')
            # set up a buffer saving multiple images (saves lots of time)
            if not self.buffer.get('init'):
                self.buffer['init'] = True
                self.Comments = self.visit(fp)
                if len(self.buffer['imagemap']) == 0:
                    self.errors = 'No valid images found in file'
                    return False
                if imagenum > len(self.buffer['imagemap']):
                    self.errors = f"Only {len(self.buffer['imagemap'])} images found in file. {imagenum} cannot be read."
                    return False
                self.buffer['misc'] = {}
                # get 1ID pixel size info. Currently an array, but this may change
                try:
                    for key in 'DetSizeX','DetSizeY':
                        self.buffer['misc'][key] = [i for i in fp['misc'][key]]
                    for key in 'DetPixelSizeX','DetPixelSizeY':
                        self.buffer['misc'][key] = [float(i) for i in fp['misc'][key]]
                except:
                    pass
            self.Data,self.Npix,self.Image = self.readDataset(fp,imagenum)
            if self.Npix == 0:
                self.errors = 'No valid images found in file'
                return False
            self.LoadImage(ParentFrame,filename,imagenum)
            self.repeatcount = imagenum 
            self.repeat = imagenum < len(self.buffer['imagemap'])
            tag = self.buffer['imagemap'][imagenum-1][0]
            self.Data['ImageSection'] = tag # save section of file here
            if GSASIIpath.GetConfigValue('debug'): print(f'Read image #{imagenum} ({tag}) from file {filename}')
            return True
        except IOError:
            print ('cannot open file '+ filename)
            return False
        finally:
            fp.close()

    def visit(self, fp):
        '''Recursively visit each node in an HDF5 file. For nodes
        ending in 'data' look at dimensions of contents. If the shape is
        length 2, 3, or 4 assume an image and index in self.buffer['imagemap']
        '''
        # "root" names where images may be found
        #    legacy: images
        #    EPICS Area Detector: data, dark, bright
        datakeywords = ['images','data','dark','bright']
        head = []
        def func(name, dset):
            if not hasattr(dset,'shape'): return # not array, can't be image
            if isinstance(dset, h5py.Dataset):
                if len(dset.shape) < 2:
                    head.append('%s: %s'%(dset.name,str(dset[()][0])))
                for datakeyword in datakeywords:
                    if dset.name.endswith(datakeyword):
                        dims = dset.shape
                        if len(dims) == 4:
                            self.buffer['imagemap'] += [(dset.name,i) for i in range(dims[1])]
                        elif len(dims) == 3:
                            self.buffer['imagemap'] += [(dset.name,i) for i in range(dims[0])]
                        elif len(dims) == 2:
                            self.buffer['imagemap'] += [(dset.name,None)]
                        else:
                            print('Skipping entry '+str(dset.name)+'. Shape is '+str(dims))
                        break
        self.buffer['imagemap'] = []
        fp.visititems(func)
        return head
        
    def readDataset(self,fp,imagenum=1):
        '''Read a specified image number from a file
        '''
        name,num = self.buffer['imagemap'][imagenum-1] # look up image in map
        dset = fp[name]
        if num is None:
            image = dset[()]
        elif len(dset.shape) == 4:
            image = dset[0,num,...]
        elif len(dset.shape) == 3:
            image = dset[num,...]
        else:
            msg = 'Unexpected image dimensions '+name
            print(msg)
            raise Exception(msg)
        sizexy = list(image.shape) 
        Npix = sizexy[0]*sizexy[1]
        # default pixel size (for APS sector 6?)
        pixelsize = [74.8,74.8]
        j = 0   # use 1st size/bin entry for all images
        try:
            misc = self.buffer['misc']
            if 'DetSizeX' in misc and 'DetSizeY' in misc:
                pixelsize = [misc[f'DetSize{i}'][j]*misc[f'DetPixelSize{i}'][j] for i in ('X','Y')]
                print(f'Using DetSize* & DetPixelSize* for Pixel size: {pixelsize}.')
            else:
                pixelsize = [misc[f'DetPixelSize{i}'][j] for i in ('X','Y')]
                print(f'Using DetPixelSize* for Pixel size: {pixelsize}.')
        except:
            print(f'No DetSize* or DetPixelSize* Pixel size defaulting to {pixelsize}.')
        data = {'pixelSize':pixelsize,'wavelength':0.15,'distance':1000.,
                'center':[sizexy[0]*0.1,sizexy[1]*0.1],'size':sizexy,'det2theta':0.0}
        for item in self.Comments:
            name,val = item.split(':',1)
            if 'wavelength' in name and 'spread' not in name:
                try:
                    data['wavelength'] = float(val)
                except ValueError:
                    pass
            elif 'distance' in name:
                data['distance'] = float(val)
            elif 'x_pixel_size' in name:
                data['pixelSize'][0] = float(val)*1000.
            elif 'y_pixel_size' in name:
                data['pixelSize'][1] = float(val)*1000.
            elif 'beam_center_x' in name: 
                data['center'][0] = float(val)
            elif 'beam_center_y' in name: 
                data['center'][1] = float(val)                
        return data,Npix,image.T
