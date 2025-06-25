# -*- coding: utf-8 -*-
'''A reader for HDF-5 files. This should be as generic as possible, but
at present this is pretty much customized for XSD-MPE (APS) uses.
'''
try:
    import h5py
except ImportError:
    h5py = None
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil
from .. import GSASIIpath

class HDF5_Reader(G2obj.ImportImage):
    '''Routine to read a HDF-5 image, typically from APS Sector 1 or 6.
    Initial version from B. Frosik/SDM.
    Updated to also handle images from APS 1-ID-C.
    '''
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
        '''Read an image from a HDF5 file. Note that images are using :meth:`readDataset`.

        When called the first time on a file, the file structure is scanned
        using :meth:`visit` to map out locations of image(s). On subsequent calls,
        if more than one image is in the file, the map of file structure (in
        buffer arg) is reused. Depending on the Config setting for HDF5selection,
        a window may be opened to allow selection of which images will be read.

        When an image is reread, the blocknum will be a list item with the location
        to be read, so the file scan can be skipped.
        '''
        try:
            fp = h5py.File(filename, 'r')
        except IOError:
            return False
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        # do we have a image number or a map to the section with the image?
        try:
            int(imagenum)
            # set up an index as to where images are found
            self.buffer = kwarg.get('buffer',{})
            if not self.buffer.get('imagemap'):
                if GSASIIpath.GetConfigValue('debug'): print('Scanning for image map')
                self.buffer['imagemap'] = []
                self.Comments = self.visit(fp)
                if len(self.buffer['imagemap']) == 0:
                    self.errors = 'No valid images found in file'
                    fp.close()
                    return False
                if imagenum > len(self.buffer['imagemap']):
                    self.errors = f"Only {len(self.buffer['imagemap'])} images found in file. {imagenum} cannot be read."
                    fp.close()
                    return False
                nsel = GSASIIpath.GetConfigValue('HDF5selection',getDefault=True)
                self.buffer['selectedImages'] = list(range(len(self.buffer['imagemap'])))
                if ParentFrame and len(self.buffer['imagemap']) > nsel and nsel >= 0:
                    import wx
                    from .. import GSASIIctrlGUI as G2G
                    choices = []
                    for loc,num,siz in self.buffer['imagemap']:
                        if num is None:
                            choices.append(f'image in {loc} size={siz}')
                        else:
                            choices.append(f'image in {loc} sec {num} size={siz}')
                    dlg = G2G.G2MultiChoiceDialog(ParentFrame,'Select images to read',
                                                      'Choose images',choices)
                    dlg.Layout()
                    dlg.SendSizeEvent()
                    if dlg.ShowModal() == wx.ID_OK:
                        self.buffer['selectedImages'] = dlg.GetSelections()
                    dlg.Destroy()
                if len(self.buffer['selectedImages']) == 0:
                    self.errors = 'No images selected from file'
                    fp.close()
                    return False
            self.buffer['selectedImages'] = self.buffer.get('selectedImages',
                                                list(range(len(self.buffer['imagemap']))))
            # get the first selected image
            while imagenum <= len(self.buffer['imagemap']):
                if imagenum-1 in self.buffer['selectedImages']:
                    del self.buffer['selectedImages'][self.buffer['selectedImages'].index(imagenum-1)]
                    break
                else:
                    imagenum += 1
            else:  # unexpected!
                self.errors = 'No images selected from file'
                fp.close()
                return False
            kwargs = {'imagenum':imagenum}
            quick = False
        except:
            kwargs = {'name':imagenum[0],'num':imagenum[1]}
            quick = True
        # we have been passed a map to images
        self.Data,self.Npix,self.Image = self.readDataset(fp,**kwargs)
        if quick:
            fp.close()
            return True
        if self.Npix == 0:
            self.errors = 'No valid images found in file'
            fp.close()
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        tag = self.buffer['imagemap'][imagenum-1][0]
        sec = self.buffer['imagemap'][imagenum-1][1]
        self.imageEntry = imagenum-1
        if sec is None:
            sec = "(none)"
            self.repeatcount = 0
        else:
            self.repeatcount = self.buffer['imagemap'][imagenum-1][1]+1
        if GSASIIpath.GetConfigValue('debug'): print(f'Read image #{imagenum} ({tag} section {sec}) from file {filename}')
        self.Data['ImageSection'] = tag # save section of file here
        # look for next image to read
        while imagenum <= len(self.buffer['imagemap']):
            if imagenum in self.buffer['selectedImages']:
                self.repeat = True
                break
            else:
                imagenum += 1
        else:
            self.repeat = False
        fp.close()
        return True

    def visit(self, fp):
        '''Recursively visit every node in an HDF5 file & look at dimensions
        of contents. If the shape is length 2, 3, or 4 assume an image
        and index in self.buffer['imagemap']
        '''
        head = []
        def func(name, dset):
            if not hasattr(dset,'shape'): return # not array, can't be image
            if isinstance(dset, h5py.Dataset):
                dims = dset.shape
                if len(dims) < 2:
                    head.append('%s: %s'%(dset.name,str(dset[()][0])))
                elif len(dims) == 4:
                    size = dims[2:]
                    self.buffer['imagemap'] += [(dset.name,i,size) for i in range(dims[1])]
                elif len(dims) == 3:
                    size = dims[1:]
                    self.buffer['imagemap'] += [(dset.name,i,size) for i in range(dims[0])]
                elif len(dims) == 2:
                    size = dims
                    self.buffer['imagemap'] += [(dset.name,None,size)]
                else:
                    print('Skipping entry '+str(dset.name)+'. Shape is '+str(dims))
        fp.visititems(func)
        return head

    def readDataset(self,fp,imagenum=1,name=None,num=None):
        '''Read a specified image number from a file
        '''
        if name is None:
            name,num,size = self.buffer['imagemap'][imagenum-1] # look up image in map
            quick = False
        else:
            quick = True
        dset = fp[name]
        if num == None:
            image = dset[()]
        elif len(dset.shape) == 4:
            image = dset[0,num,...]
        elif len(dset.shape) == 3:
            image = dset[num,...]
        else:
            msg = 'Unexpected image dimensions '+name
            print(msg)
            raise Exception(msg)
        if quick:
            return {},None,image.T
        sizexy = list(image.shape)
        Npix = sizexy[0]*sizexy[1]
        # default pixel size (for APS sector 6?)
        pixelsize = [74.8,74.8]
        j = 0   # use 1st size/bin entry for all images
        # get 1ID pixel size info. Currently an array, but this may change
        try:
            misc = {}
            for key in 'DetSizeX','DetSizeY':
                misc[key] = [i for i in fp['misc'][key]]
            for key in 'DetPixelSizeX','DetPixelSizeY':
                misc[key] = [float(i) for i in fp['misc'][key]]
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
