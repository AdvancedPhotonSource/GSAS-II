# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_GE.py $
# $Id: G2img_GE.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import os
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
class GE_ReaderClass(G2obj.ImportImage):
    '''Routine to read a GE image, typically from APS Sector 1.
        
        The image files may be of form .geX (where X is ' ', 1, 2, 3, 4 or 5),
        which is a raw image from the detector. These files may contain more
        than one image and have a rudimentary header. 
        Files with extension .sum or .cor are 4 byte integers/pixel, one image/file.
        Files with extension .avg are 2 byte integers/pixel, one image/file.
    '''

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.sum','.cor','.cor32','.avg','.ge','.ge1','.ge2','.ge3','.ge4','.ge5'),
            strictExtension=True,
            formatName = 'GE image',
            longFormatName = 'Summed GE image file'
            )

    def ContentsValidator(self, filename):
        '''just a test on file size
        '''
        if '.sum' not in str(filename):
            try:
                fp = open(filename,'rb')
                statinfo = os.stat(str(fp).split("'")[1])
                fsize = statinfo.st_size
                self.nimages = (fsize-8192)/(2*2048**2)
                fp.close()
            except:
                return False    #bad file size
        return True
        
    def Reader(self,filename, ParentFrame=None, **kwarg):
        '''Read using GE file reader, :func:`GetGEsumData`
        '''
        #rdbuffer = kwarg.get('buffer')
        imagenum = kwarg.get('blocknum')
        #sum = kwarg.get('sum')
        if imagenum is None: imagenum = 1
        self.Comments,self.Data,self.Npix,self.Image,more = \
            GetGEsumData(self,filename,imagenum=imagenum)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        self.repeatcount = imagenum
        self.repeat = more
        return True

class GEsum_ReaderClass(G2obj.ImportImage):
    '''Routine to read multiple GE images & sum them, typically from APS Sector 1.
        
        The image files may be of form .geX (where X is ' ', 1, 2, 3, 4 or 5),
        which is a raw image from the detector. These files may contain more
        than one image and have a rudimentary header. 
        Files with extension .sum or .cor are 4 byte integers/pixel, one image/file.
        Files with extension .avg are 2 byte integers/pixel, one image/file.
    '''

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.ge1','.ge2','.ge3','.ge4','.ge5'),
            strictExtension=True,
            formatName = 'sum GE multi-image',
            longFormatName = 'sum of GE multi-image file'
            )

    def ContentsValidator(self, filename):
        '''just a test on file size
        '''
        try:
            fp = open(filename,'rb')
            statinfo = os.stat(str(fp).split("'")[1])
            fsize = statinfo.st_size
            nimages = (fsize-8192)/(2*2048**2)
            fp.close()
        except:
            return False    #bad file size
        return True
        
    def Reader(self,filename, ParentFrame=None, **kwarg):
        '''Read using GE file reader, :func:`GetGEsumData`
        '''
        #rdbuffer = kwarg.get('buffer')
        imagenum = kwarg.get('blocknum')
        if imagenum is None: imagenum = 1
        self.Comments,self.Data,self.Npix,self.Image,more = GetGEsumData(
            self,filename,imagenum=imagenum,sum=True)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename,imagenum)
        self.repeatcount = imagenum
        self.repeat = more
        return True

def GetGEsumData(self,filename,imagenum=1,sum=False):
    '''Read G.E. detector images from various files as produced at 1-ID and
    with Detector Pool detector. Also sums multiple image files if desired
    '''
    import struct as st
    import platform
    if '2' in platform.python_version_tuple()[0]:
        import cPickle
    else:
        import pickle as cPickle
    import time
    more = False
    time0 = time.time()
    File = open(filename,'rb')
    if filename.split('.')[-1] in ['sum','cor32']:
        head = ['GE detector sum/corrected data from APS 1-ID',]
        sizexy = [2048,2048]
        Npix = sizexy[0]*sizexy[1]
        image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.float32),dtype=np.int32)
    elif filename.split('.')[-1] in ['avg','cor']:
        File.seek(0,2)
        last = File.tell()
        pos = last-2*(2048**2)
        File.seek(pos)
        head = ['GE detector avg or cor data from APS 1-ID',]
        sizexy = [2048,2048]
        Npix = sizexy[0]*sizexy[1]
        image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.int16),dtype=np.int32)
    else:
        head = ['GE detector raw data',]
        File.seek(18)
        size,nframes = st.unpack('<ih',File.read(6))
        # number of frames seems to be 3 for single-image files
        if size != 2048:
            print('Warning GE image size unexpected: '+str(size))
            print('Assumed 2048x2048')
            size = 2048
            statinfo = os.stat(str(File).split("'")[1])
            fsize = statinfo.st_size
            nframes = (fsize-8192)/(2*2048**2)
#            return 0,0,0,0,False # probably should quit now
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
        image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.int16),dtype=np.int32)
        if len(image) != sizexy[1]*sizexy[0]:
            print('not enough images while reading GE file: '+filename+'image #'+str(imagenum))
            return 0,0,0,0,False            
        head += ['file: '+filename+' image #'+str(imagenum),]
        if sum:    #will ignore imagenum
            print ('Frames to read %d,'%(nframes),end='')
            while nframes > 1: #OK, this will sum the frames.
                try:
                    image += np.array(np.frombuffer(File.read(2*Npix),dtype=np.int16),dtype=np.int32)
                except ValueError:
                    break
                nframes -= 1
                print ('%d,'%(nframes),end='')
            print ('') 
            more = False
            filename = os.path.splitext(filename)[0]+'.G2img'
            File = open(filename,'wb')
            Data = {'pixelSize':[200.,200.],'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy}
            image = np.reshape(image,(sizexy[1],sizexy[0]))
            cPickle.dump([head,Data,Npix,image],File,1)
            File.close()
            self.sumfile = filename
            self.formatName = 'GSAS-II image'
            sum = False
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':[200.,200.],'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy,'det2theta':0.0}
    File.close()
    if GSASIIpath.GetConfigValue('debug'):
        print ('Image read time %.2fs'%(time.time()-time0))
        print ('Read GE file: '+filename+' image #'+'%04d'%(imagenum))
    return head,data,Npix,image,more
