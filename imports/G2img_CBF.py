# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2016-01-22 13:05:12 -0600 (Fri, 22 Jan 2016) $
# $Author: toby $
# $Revision: 2133 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_CBF.py $
# $Id: G2img_CBF.py 2133 2016-01-22 19:05:12Z toby $
########### SVN repository information ###################
'''
*Module G2img_CBF: .cbf cif image file*
--------------------------------------

'''

import sys
import os
import time
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2133 $")
class CBF_ReaderClass(G2IO.ImportImage):
    '''Routine to read a Read cif image data .cbf file.
    This is used by Pilatus. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.cbf',),
            strictExtension=True,
            formatName = 'CBF image (Slow!)',
            longFormatName = 'CIF Binary Data Format image file (NB: Slow!)'
            )

    def ContentsValidator(self, filepointer):        
        '''no test used at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        '''Read using Bob's routine :func:`GetCbfData`
        '''
        self.Comments,self.Data,self.Npix,self.Image = GetCbfData(self,filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
        
def GetCbfData(self,filename):    
    'Read cif binarydetector data cbf file'
    
    import numpy as np
    import cPickle
    
    def analyse(stream):
        """
        Analyze a stream of char with any length of exception (2,4, or 8 bytes integers)

        @return list of NParrays
        """
        listnpa = []
        key16 = "\x80"
        key32 = "\x00\x80"
        key64 = "\x00\x00\x00\x80"
        idx = 0
        shift = 1
        position = 0
        while True:
            lns = len(stream)
            idx = stream.find(key16)
            if idx == -1:
                listnpa.append(np.fromstring(stream, dtype="int8"))
                break
            listnpa.append(np.fromstring(stream[:idx], dtype="int8"))
            position += listnpa[-1].size

            if stream[idx + 1:idx + 3] == key32:
                if stream[idx + 3:idx + 7] == key64:
                    listnpa.append(np.fromstring(stream[idx + 7:idx + 15], dtype="int64"))
                    position += 1
#                    print "loop64 x=%4i y=%4i in idx %4i lns %4i value=%s" % ((position % 2463), (position // 2463), idx, lns, listnpa[-1])
                    shift = 15
                else: #32 bit int
                    listnpa.append(np.fromstring(stream[idx + 3:idx + 7], dtype="int32"))
                    position += 1
#                    print "loop32 x=%4i y=%4i in idx %4i lns %4i value=%s" % ((position % 2463), (position // 2463), idx, lns, listnpa[-1])
                    shift = 7
            else: #int16 
                listnpa.append(np.fromstring(stream[idx + 1:idx + 3], dtype="int16"))
                position += 1
#                print "loop16 x=%4i y=%4i in idx %4i lns %4i value=%s" % ((position % 2463), (position // 2463), idx, lns, listnpa[-1])
                shift = 3
            stream = stream[idx + shift:]
        return  listnpa
    
    
    if GSASIIpath.GetConfigValue('debug'):
        print 'Read cif binary detector data cbf file: ',filename
    File = open(filename,'rb')
    sizexy = [0,0]
    pixSize = [154,154]     #Pixium4700?
    cent = [0,0]
    wave = 1.54187  #default <CuKa>
    dist = 1000.
    byteOrd = '<'
    stream = File.read()
    File.close()
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
    cent = [cent[0]*pixSize[0]/1000.,cent[1]*pixSize[1]/1000.]
    compImage = stream[imageBeg:imageBeg+compImageSize]
#    GSASIIpath.IPyBreak()
    time0 = time.time()
    image = np.hstack(analyse(compImage)).cumsum()
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    print 'import time:',time.time()-time0
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist,'center':cent,'size':sizexy}
    Npix = sizexy[0]*sizexy[1]
    
    filename = os.path.splitext(filename)[0]+'.G2img'
    File = open(filename,'wb')
    cPickle.dump([head,data,Npix,image],File,1)
    File.close()
    self.sumfile = filename
    self.formatName = 'GSAS-II image'
    return head,data,Npix,image
        
