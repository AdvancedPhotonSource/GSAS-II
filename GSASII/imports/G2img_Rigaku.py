# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
import os
from .. import GSASIIobj as G2obj
import numpy as np
class Rigaku_ReaderClass(G2obj.ImportImage):
    '''Routine to read a Rigaku R-Axis IV image file.
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.stl',),
            strictExtension=True,
            formatName = 'Rigaku image',
            longFormatName = 'Read Rigaku R-Axis IV image file'
            )

    def ContentsValidator(self, filename):        
        '''Test by checking if the file size makes sense.
        '''
        fileSize = os.stat(filename).st_size
        Npix = (fileSize-6000)/2
        if Npix == 9000000 or Npix == 2250000 or Npix == 36000000:
            return True
        return False # not valid size
        
    def Reader(self,filename, ParentFrame=None, **unused):
        self.Comments,self.Data,self.Npix,self.Image = GetRigaku(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True

def GetRigaku(filename,imageOnly=False):
    'Read Rigaku R-Axis IV image file'
    import array as ar
    if not imageOnly:
        print ('Read Rigaku R-Axis IV file: '+filename)   
    File = open(filename,'rb')
    fileSize = os.stat(filename).st_size
    Npix = (fileSize-6000)/2
    File.read(6000)
    head = ['Rigaku R-Axis IV detector data',]
    image = np.array(ar.array('H',File.read(fileSize-6000)),dtype=np.int32)
    print ('%s %s'%(fileSize,str(image.shape)))
    print (head)
    if Npix == 9000000:
        sizexy = [3000,3000]
        pixSize = [100.,100.]        
    elif Npix == 2250000:
        sizexy = [1500,1500]
        pixSize = [200.,200.]
    else:
        sizexy = [6000,6000]
        pixSize = [50.,50.] 
    image = np.reshape(image,(sizexy[1],sizexy[0]))        
    data = {'pixelSize':pixSize,'wavelength':1.5428,'distance':250.0,'center':[150.,150.],'size':sizexy,'det2theta':0.0}  
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
    
