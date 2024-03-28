# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_PILTIF.py $
# $Id: G2img_PILTIF.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
The metadata can be specified in a file with the same name and path as 
the TIFF file except that the the extension is .metadata.

The contents of that file are a series of lines of form::

     keyword = value

Note that capitalization of keywords is ignored. Defined keywords are in table below. Any line 
without one of these keywords will be ignored. 

.. Next command allows \\AA to be used in HTML

.. tabularcolumns:: |l|p{4.5in}|

==============================  ====================================================
  keyword                        explanation
==============================  ====================================================
wavelength                       Wavelength in :math:`\\AA`
distance                         Distance to sample in mm 
polarization                     Percentage polarized in horizontal plane 
sampleChangerCoordinate          Used for sample changers to track sample
pixelSizeX                       Pixel size in X direction (microns)
pixelSizeY                       Pixel size in Y direction (microns)
CenterPixelX                     Location of beam center as a pixel number (in X)
CenterPixelY                     Location of beam center as a pixel number (in X)
==============================  ====================================================

'''

from __future__ import division, print_function
import struct as st
import numpy as np
import time
import GSASIIobj as G2obj
import GSASIIpath
import GSASIIfiles as G2fil
import G2img_1TIF
DEBUG = False
GSASIIpath.SetVersionNumber("$Revision: 5577 $")

class TIF_LibraryReader(G2obj.ImportImage):
    '''Reads TIF files using a standard library routine. Metadata (such as pixel 
    size) must be specified by user, either in GUI or via a metadata file. 
    The library TIF reader can handle compression and other things that are not 
    commonly used at beamlines. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.tif','.tiff'),
            strictExtension=True,
            formatName = 'Standard TIF image; metadata req.',
            longFormatName = 'TIFF images read with standard library (metadata must be supplied)'
            )
        self.scriptable = True

    def ContentsValidator(self, filename):
        '''Does the header match the required TIF header?
        '''
        return G2img_1TIF.TIFValidator(filename)
    
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read the TIF file using the PIL/Pillow reader and give the 
        user a chance to edit the likely wrong default image parameters. 
        '''
        import PIL.Image as PI
        self.Image = PI.open(filename,mode='r')
        self.Npix = self.Image.size
        if ParentFrame:
            self.SciPy = True
            self.Comments = ['no metadata']
            self.Data = {'wavelength': 0.1, 'pixelSize': [200., 200.], 'distance': 100.0}
            self.Data['size'] = list(self.Image.size)
            self.Data['center'] = [int(i/2) for i in self.Image.size]
            try:
                Meta = open(filename+'.metadata','r')
                head = Meta.readlines()
                for line in head:
                    line = line.strip()
                    try:
                        if '=' not in line: continue
                        keyword = line.split('=')[0].strip()
                        if 'wavelength' == keyword.lower():
                            self.Data['wavelength'] = float(line.split('=')[1])
                        elif 'distance' == keyword.lower():
                            self.Data['distance'] = float(line.split('=')[1])
                        elif 'polarization' == keyword.lower():
                            polarization = float(line.split('=')[1])
                            self.Data['PolaVal'] = [polarization,False]
                        elif 'samplechangercoordinate' == keyword.lower():
                            self.Data['samplechangerpos'] = float(line.split('=')[1])
                        elif 'pixelsizex' == keyword.lower():
                            self.Data['pixelSize'][0] = float(line.split('=')[1])
                        elif 'pixelsizey' == keyword.lower():
                            self.Data['pixelSize'][1] = float(line.split('=')[1])
                        elif 'centerpixelx' == keyword.lower():
                            self.Data['center'][0] = float(line.split('=')[1])
                        elif 'centerpixely' == keyword.lower():
                            self.Data['center'][1] = float(line.split('=')[1])
                    except:
                        G2fil.G2Print('error reading metadata: '+line)
                Meta.close()
                self.SciPy = False
            except IOError:
                G2fil.G2Print ('no metadata file found - image params must be set manually')
                head = ['no metadata file found',]
        if self.Npix == 0:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
