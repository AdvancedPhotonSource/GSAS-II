# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_ADSC.py $
# $Id: G2img_ADSC.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import GSASIIobj as G2obj
import GSASIIpath
import numpy as np
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
class ADSC_ReaderClass(G2obj.ImportImage):
    '''Reads an ADSC .img file
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.img',),
            strictExtension=True,
            formatName = 'ADSC image',
            longFormatName = 'ADSC image file'
            )

    def ContentsValidator(self, filename):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        self.Comments,self.Data,self.Npix,self.Image = GetImgData(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True

def GetImgData(filename,imageOnly=False):
    'Read an ADSC image file'
    import array as ar
    if not imageOnly:
        print ('Read ADSC img file: '+filename)
    File = open(filename,'rb')
    head = File.read(511)
    lines = head.split('\n')
    head = []
    center = [0,0]
    for line in lines[1:-2]:
        line = line.strip()[:-1]
        if line:
            if 'SIZE1' in line:
                size = int(line.split('=')[1])
                Npix = size*size
            elif 'WAVELENGTH' in line:
                wave = float(line.split('=')[1])
            elif 'BIN' in line:
                if line.split('=')[1] == '2x2':
                    pixel=(102,102)
                else:
                    pixel = (51,51)
            elif 'DISTANCE' in line:
                distance = float(line.split('=')[1])
            elif 'CENTER_X' in line:
                center[0] = float(line.split('=')[1])
            elif 'CENTER_Y' in line:
                center[1] = float(line.split('=')[1])
            head.append(line)
    data = {'pixelSize':pixel,'wavelength':wave,'distance':distance,'center':center,'size':[size,size],'det2theta':0.0}
    image = []
    pos = 512
    File.seek(pos)
    image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    image = np.reshape(image,(size,size))
    File.close()
    if imageOnly:
        return image
    else:
        return lines[1:-2],data,Npix,image
       
