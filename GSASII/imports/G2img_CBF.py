# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_CBF.py $
# $Id: G2img_CBF.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import time
import GSASIIobj as G2obj
import GSASIIpath
import struct as st
import numpy as np
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
class CBF_ReaderClass(G2obj.ImportImage):
    '''Routine to read a Read cif image data .cbf file.
    This is used by Pilatus. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.cbf',),
            strictExtension=True,
            formatName = 'CBF image',
            longFormatName = 'CIF Binary Data Format image file (NB: Slow!)'
            )

    def ContentsValidator(self, filename):        
        '''no test used at this time
        '''
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read using Bob's routine :func:`GetCbfData`
        '''
        self.Comments,self.Data,self.Npix,self.Image = GetCbfData(self,filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
        
def GetCbfData(self,filename):    
    'Read cif binarydetector data cbf file'
    
    import unpack_cbf as cbf
    if GSASIIpath.GetConfigValue('debug'):
        print ('Read cif binary detector data cbf file: '+filename)
    File = open(filename,'rb')
    sizexy = [0,0]
    pixSize = [172,172]     #Pixium4700?
    cent = [0,0]
    wave = 1.54187  #default <CuKa>
    det2theta = 0.0
    dist = 1000.
    byteOrd = '<'
    stream = File.read()
    if 'bytes' in str(type(stream)):
        stream = stream.decode('latin-1')
    starter = '\x0c\x1a\x04\xd5'
    imageBeg = stream.find(starter)+4
    head = stream[:imageBeg]
    term = '\r\n'       #CR-LF
    if term in head:
        pass
    else:
        term = '\n' #LF only
    head = head.split(term)
    for line in head:
        fields = line.split()
        if 'Wavelength' in line:
            wave = float(fields[2])
        elif 'Detector_distance' in line:
            dist = float(fields[2])*1000.
        elif 'Pixel_size' in line:
            pixSize = [float(fields[2])*1.e6,float(fields[5])*1.e6]
        elif 'Beam_xy' in line:
            cent = [float(fields[2].strip('(').strip(',')),float(fields[3].strip(')'))]
        elif 'X-Binary-Size:' in line:
            compImageSize = int(fields[1])
        elif 'BIG_ENDIAN' in line:
            byteOrd = '>'
        elif 'Fastest-Dimension' in line:
            sizexy[0] = int(fields[1])
        elif 'Second-Dimension' in line:
            sizexy[1] = int(fields[1])
        elif 'Number-of-Elements' in line:
            Npix = int(fields[1])
        elif 'Detector_2theta' in line:
            det2theta = float(fields[2])
    nxy = sizexy[0]*sizexy[1]
    cent = [cent[0]*pixSize[0]/1000.,cent[1]*pixSize[1]/1000.]
    File.seek(0)
    img = File.read()[imageBeg:imageBeg+compImageSize]
    File.close()
    nimg = len(img)
    image = np.zeros(nxy,dtype=np.int32)
    time0 = time.time()
    if 'bytes' in str(type(img)):
        img = np.frombuffer(img,dtype=np.uint8)
        image = cbf.unpack_cbf3(nimg,img,nxy,image)
    else:
        image = cbf.unpack_cbf(nimg,img,nxy,image)
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    print ('import time: %.3f'%(time.time()-time0))
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist,'center':cent,'size':sizexy,'det2theta':det2theta}
    Npix = sizexy[0]*sizexy[1]
    
    return head,data,Npix,image
        
