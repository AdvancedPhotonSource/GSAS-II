# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_EDF.py $
# $Id: G2img_EDF.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import os
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
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
    wave = 1.54187  #default <CuKa>
    dist = 1000.
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
        elif 'Size' in line:
            imSize = int(fields[2])
#        elif 'DataType' in lines:
#            dType = fields[2]
        elif 'pixel_size_x' in line:
            pixSize[0] = float(fields[2])
        elif 'pixel_size_y' in line:
            pixSize[1] = float(fields[2])
        elif 'beam_center_x' in line:
            cent[0] = float(fields[2])
        elif 'beam_center_y' in line:
            cent[1] = float(fields[2])
        elif 'refined_distance' in line:
            dist = float(fields[2])
        if line:
            head.append(line)
        else:   #blank line at end of header
            break  
    File.seek(fileSize-imSize)
    if dType == 'UnsignedShort':        
        image = np.array(np.frombuffer(File.read(imSize),dtype=np.int16),dtype=np.int32)
    else:
        image = np.array(np.frombuffer(File.read(imSize),dtype=np.int32),dtype=np.int32)
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist,'center':cent,'size':sizexy}
    Npix = sizexy[0]*sizexy[1]
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
        
