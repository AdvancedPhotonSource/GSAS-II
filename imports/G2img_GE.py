# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2img_GE: summed GE image file*
---------------------------------------

This shows an example of an importer that will handle files with
more than a single image. 

'''

import sys
import os
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
class GE_ReaderClass(G2IO.ImportImage):
    '''Routine to read a GE image, typically from APS Sector 1.
        
        The image files may be of form .geX (where X is ' ', 1, 2, 3 or 4),
        which is a raw image from the detector. These files may contain more
        than one image and have a rudimentary header. 
        Files with extension .sum or .cor are 4 byte integers/pixel, one image/file.
        Files with extension .avg are 2 byte integers/pixel, one image/file.
    '''

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.sum','.cor','.avg','.ge','.ge1','.ge2','.ge3','.ge4'),
            strictExtension=True,
            formatName = 'GE image',
            longFormatName = 'Summed GE image file'
            )

    def ContentsValidator(self, filepointer):
        '''just a test on file size
        '''
        if '.sum' not in str(filepointer):
            try:
                statinfo = os.stat(str(filepointer).split("'")[1])
                fsize = statinfo.st_size
                self.nimages = (fsize-8192)/(2*2048**2)
            except:
                return False    #bad file size
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        '''Read using GE file reader, :func:`GetGEsumData`
        '''
        #rdbuffer = kwarg.get('buffer')
        imagenum = kwarg.get('blocknum')
        #sum = kwarg.get('sum')
        if imagenum is None: imagenum = 1
        self.Comments,self.Data,self.Npix,self.Image,more = \
            GetGEsumData(filename,imagenum=imagenum)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        self.repeatcount = imagenum
        self.repeat = more
        return True

class GEsum_ReaderClass(G2IO.ImportImage):
    '''Routine to read multiple GE images & sum them, typically from APS Sector 1.
        
        The image files may be of form .geX (where X is ' ', 1, 2, 3 or 4),
        which is a raw image from the detector. These files may contain more
        than one image and have a rudimentary header. 
        Files with extension .sum or .cor are 4 byte integers/pixel, one image/file.
        Files with extension .avg are 2 byte integers/pixel, one image/file.
    '''

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.ge1','.ge2','.ge3','.ge4'),
            strictExtension=True,
            formatName = 'GE multi-image summed',
            longFormatName = 'sum of GE multi-image file'
            )

    def ContentsValidator(self, filepointer):
        '''just a test on file size
        '''
        try:
            statinfo = os.stat(str(filepointer).split("'")[1])
            fsize = statinfo.st_size
            nimages = (fsize-8192)/(2*2048**2)
        except:
            return False    #bad file size
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        '''Read using GE file reader, :func:`GetGEsumData`
        '''
        #rdbuffer = kwarg.get('buffer')
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        self.Comments,self.Data,self.Npix,self.Image,more = \
            GetGEsumData(filename,imagenum=imagenum,sum=True)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        self.repeatcount = imagenum
        self.repeat = more
        return True

def GetGEsumData(filename,imagenum=1,sum=False):
    '''Read G.E. detector images from various files as produced at 1-ID and
    with Detector Pool detector. Also sums multiple image files if desired
    '''
    import struct as st
    import array as ar
    more = False
    File = open(filename,'rb')
    if '.sum' in filename or '.cor' in filename:
        head = ['GE detector sum or cor data from APS 1-ID',]
        sizexy = [2048,2048]
        Npix = sizexy[0]*sizexy[1]
        image = np.array(ar.array('f',File.read(4*Npix)),dtype=np.int32)
    elif '.avg' in filename:
        head = ['GE detector avg or ge* data from APS 1-ID',]
        sizexy = [2048,2048]
        Npix = sizexy[0]*sizexy[1]
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    else:
        head = ['GE detector raw data',]
        File.seek(18)
        size,nframes = st.unpack('<ih',File.read(6))
        # number of frames seems to be 3 for single-image files
        if size != 2048:
            print('Warning GE image size unexpected: '+str(size))
            return 0,0,0,0,False # probably should quit now
        if imagenum > nframes:
            print('Error: attempt to read image #'+str(imagenum)+
                  ' from file with '+str(nframes)+' images.')
            return 0,0,0,0,False
        elif imagenum < nframes:
            more = True
        sizexy = [2048,2048]
        Npix = sizexy[0]*sizexy[1]
        pos = 8192 + (imagenum-1)*2*Npix
        File.seek(pos)
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
        if len(image) != sizexy[1]*sizexy[0]:
            print('not enough images while reading GE file: '+filename+'image #'+str(imagenum))
            return 0,0,0,0,False            
        head += ['file: '+filename+' image #'+str(imagenum),]
        if sum:    #will ignore imagenum
            while nframes > 1: #OK, this will sum the frames.
                image += np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
                nframes -= 1
            more = False
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':[200,200],'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy}
    File.close()
    if GSASIIpath.GetConfigValue('debug'):
        print 'Read GE file: '+filename+' image #'+str(imagenum)
    return head,data,Npix,image,more
