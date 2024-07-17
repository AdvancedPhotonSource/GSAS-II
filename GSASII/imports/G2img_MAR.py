# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_MAR.py $
# $Id: G2img_MAR.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import platform
import time
import GSASIIobj as G2obj
import GSASIIpath
import numpy as np
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
class MAR_ReaderClass(G2obj.ImportImage):
    '''Routine to read several MAR formats, .mar3450,.mar2300,.mar2560
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.mar3450','.mar2300','.mar2560'),
            strictExtension=True,
            formatName = 'MAR image',
            longFormatName = 'MAR Research 345, 230 and 256 image files'
            )

    def ContentsValidator(self, filename):
        '''no test at this time
        '''
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        self.Comments,self.Data,self.Npix,self.Image = GetMAR345Data(filename)
        if self.Npix == 0 or not self.Comments:
            return False
        self.LoadImage(ParentFrame,filename)
        return True

def GetMAR345Data(filename,imageOnly=False):
    'Read a MAR-345 image plate image'
    try:
        import pack_f as pf
    except:
        print ('**** ERROR - Unable to load the GSAS-II MAR image decompression, pack_f')
        return None,None,None,None

    if not imageOnly:
        print ('Read Mar345 file: '+filename)
    File = open(filename,'rb')
    head = File.read(4095).decode(encoding='latin-1')
    lines = head[128:].split('\n')
    head = []
    for line in lines:
        line = line.strip()
        if 'PIXEL' in line:
            values = line.split()
            pixel = (int(values[2]),int(values[4]))     #in microns
        elif 'WAVELENGTH' in line:
            wave = float(line.split()[1])
        elif 'DISTANCE' in line:
            distance = float(line.split()[1])           #in mm
            if not distance:
                distance = 500.
        elif 'CENTER' in line:
            values = line.split()
            center = [float(values[2])/10.,float(values[4])/10.]    #make in mm from pixels
        if line: 
            head.append(line)
    data = {'pixelSize':pixel,'wavelength':wave,'distance':distance,'center':center,'det2theta':0.0}
    for line in head:
        if 'FORMAT' in line[0:6]:
            items = line.split()
            sizex = int(items[1])
            Npix = int(items[3])
            sizey = int(Npix/sizex)
    pos = 4096
    data['size'] = [sizex,sizey]
    File.seek(pos)
    line = File.read(8).decode(encoding='latin-1')
    while 'CCP4' not in line:       #get past overflow list for now
        line = File.read(8).decode(encoding='latin-1')
        pos += 8
    pos += 37
    File.seek(pos)
    image = np.zeros(shape=(sizex,sizey),dtype=np.int32)    
    time0 = time.time()
    if '2' in platform.python_version_tuple()[0]:
        raw = File.read()
        image = np.flipud(pf.pack_f(len(raw),raw,sizex,sizey,image).T)  #transpose to get it right way around & flip
    else:
        raw = np.frombuffer(File.read(),dtype=np.uint8)
        image = np.flipud(pf.pack_f3(len(raw),raw,sizex,sizey,image).T)  #transpose to get it right way around & flip
    if GSASIIpath.GetConfigValue('debug'):
        print ('image read time: %.3f'%(time.time()-time0))
    File.close()
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
        
