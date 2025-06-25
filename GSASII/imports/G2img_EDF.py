# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
import os
import numpy as np
from .. import GSASIIobj as G2obj
class EDF_ReaderClass(G2obj.ImportImage):
    '''Routine to read a Read European detector data .edf file.
    This is a particularly nice standard. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.edf',),
            strictExtension=True,
            formatName = 'EDF image',
            longFormatName = 'European Data Format image file'
            )

    def ContentsValidator(self, filename):        
        '''no test used at this time
        '''
        return True
        
    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        self.Comments,self.Data,self.Npix,self.Image = GetEdfData(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True

def GetEdfData(filename,imageOnly=False):    
    'Read European detector data edf file'
    if not imageOnly:
        print ('Read European detector data edf file: '+filename)
    File = open(filename,'rb')
    fileSize = os.stat(filename).st_size
    head = File.read(3072).decode(encoding='latin-1')
    lines = head.split('\n')
    sizexy = [0,0]
    pixSize = [154,154]     #Pixium4700?
    cent = [0,0]
    pcent = [0,0]
    wave = 1.54187  #default <CuKa>
    dist = 1000.
    Temperature = 300.
    head = ['European detector data',]
    for line in lines:
        line = line.replace(';',' ').strip()
        fields = line.split()
        if 'Dim_1' in line:
            sizexy[0] = int(fields[2])
        elif 'Dim_2' in line:
            sizexy[1] = int(fields[2])
        elif 'DataType' in line:
            dType = fields[2]
        elif 'wavelength' in line:
            wave = float(fields[2])
        elif 'Wavelength' in line:
            wave = float(fields[2])*1.e10
        elif 'Size' == fields[0]:
            imSize = int(fields[2])
#        elif 'DataType' in lines:
#            dType = fields[2]
        elif 'pixel_size_x' in line:
            pixSize[0] = float(fields[2])
        elif 'PSize_1' == fields[0]:
            pixSize[0] = float(fields[2])*1.e6
        elif 'pixel_size_y' in line:
            pixSize[1] = float(fields[2])
        elif 'PSize_2' == fields[0]:
            pixSize[1] = float(fields[2])*1.e6
        elif 'beam_center_x' in line: 
            cent[0] = float(fields[2])
        elif 'Center_1' in line:
            pcent[0] = float(fields[2])
        elif 'beam_center_y' in line:
            cent[1] = float(fields[2])
        elif 'Center_2' in line:
            pcent[1] = float(fields[2])
        elif 'refined_distance' in line:
            dist = float(fields[2])
        elif 'Temperature' in line:
            Temperature = float(fields[2])
        elif 'SampleDistance' in line:
            dist = float(fields[2])
        if line:
            head.append(line)
        else:   #blank line at end of header
            break
    if any(pcent):
        cent[0] = pcent[0]*pixSize[0]/1000.
        cent[1] = pcent[1]*pixSize[1]/1000.
    File.seek(fileSize-imSize)
    if dType == 'UnsignedShort':        
        image = np.array(np.frombuffer(File.read(imSize),dtype=np.int16),dtype=np.int32)
    else:
        image = np.array(np.frombuffer(File.read(imSize),dtype=np.int32),dtype=np.int32)
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist*1000.,'center':cent,
            'size':sizexy,'Temperature':Temperature+273.}
    Npix = sizexy[0]*sizexy[1]
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
        
