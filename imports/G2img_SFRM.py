# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2017-12-26 20:18:10 -0600 (Tue, 26 Dec 2017) $
# $Author: toby $
# $Revision: 3207 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_CBF.py $
# $Id: G2img_CBF.py 3207 2017-12-27 02:18:10Z toby $
########### SVN repository information ###################
'''
*Module G2img_SFRM: .sfrm image file*
---------------------------------------

'''

from __future__ import division, print_function
import time
import GSASIIobj as G2obj
import GSASIIpath
import numpy as np
GSASIIpath.SetVersionNumber("$Revision: 3207 $")
class SFRM_ReaderClass(G2obj.ImportImage):
    '''Routine to read a Read Bruker Advance image data .sfrm file.
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.sfrm',),
            strictExtension=True,
            formatName = 'SFRM image',
            longFormatName = 'Bruker SFRM Binary Data Format image file'
            )

    def ContentsValidator(self, filename):        
        '''no test used at this time
        '''
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read using Bob's routine :func:`GetSFRMData`
        '''
        self.Comments,self.Data,self.Npix,self.Image = GetSFRMData(self,filename)
        if self.Npix == 0 or not len(self.Comments):
            return False
        self.LoadImage(ParentFrame,filename)
        return True
        
def GetSFRMData(self,filename):    
    'Read cbf compressed binarydetector data sfrm file'
    
    if GSASIIpath.GetConfigValue('debug'):
        print ('Read cbf compressed binary detector data sfrm file: '+filename)
    File = open(filename,'rb')
    sizexy = [0,0]
    pixSize = [135.3,135.3]     #Pixium4700?
    cent = [0,0]
    wave = 1.54187  #default <CuKa>
    dist = 250.
    stream = File.read()
    if 'bytes' in str(type(stream)):
        stream = stream.decode('latin-1')
    starter = 'IMG: '
    meanwaves = {'Cu':1.54051,'Ti':2.74841,'Cr':2.28962,'Fe':1.93597,
        'Co':1.78892,'Mo':0.70926,'Ag':0.559363}
    imageBeg = stream.find(starter)+4
    head = np.array(list(stream[:imageBeg].split('CFR:')[0]))
    head = head.reshape(-1,80)
    lines = []
    for line in head:
        line = ''.join(line)
        lines.append(line)
        fields = line.split(':')[1].split()
        if 'TARGET' in line:
            wave = meanwaves[fields[0]]
            target = fields[0].capitalize()
        elif 'DISTANC' in line:
            dist = float(fields[1])*10.
        elif 'ANGLES' in line:
            twoth = float(fields[0])
        elif 'CENTER' in line:
            cent = [float(fields[0]),float(fields[1])]
        elif 'NROWS' in line:
            sizexy[1] = int(fields[0])
        elif 'NCOLS' in line:
            sizexy[0] = int(fields[0])
        elif 'FORMAT' in line:
            frmt = int(fields[0])
        elif 'HDRBLKS' in line:
            imageBeg = 512*int(fields[0])
        elif 'NOVERFL' in line:
            Nunder = int(fields[0])
            N2byte = 2*int(fields[1])
            if N2byte%16:
                N2byte = (N2byte//16+1)*16
            N4byte = 4*int(fields[2])
            if N4byte%16:
                N4byte = (N4byte//16+1)*16
    if frmt == 86:
        lines = ['FORMAT 86 Bruker files currently not readible by GSAS-II',]
        return lines,0,0,0
    nxy = sizexy[0]*sizexy[1]
    cent = [cent[0]*pixSize[0]/1000.,cent[1]*pixSize[1]/1000.]
    cent[0] += dist*np.tan(np.pi*twoth/180.)
    File.seek(imageBeg)
    img = File.read(nxy)
    img2byte = File.read(N2byte)
    img4byte = File.read(N4byte)
    time0 = time.time()
    img = np.array(np.frombuffer(img,dtype='u1'),dtype=np.int32)
    img2byte = np.array(np.frombuffer(img2byte,dtype='u2'),dtype=np.int32)
    img4byte = np.array(np.frombuffer(img4byte,dtype='u4'),dtype=np.int32)
    ins2byte = np.argwhere(img==255)
    for j,i in enumerate(list(ins2byte)):
        img[i] = img2byte[j]
    ins4byte = np.argwhere(img==65535)
    for j,i in enumerate(list(ins4byte)):
        img[i] = img4byte[j]
    image = np.reshape(img,(sizexy[1],sizexy[0]))
    print ('import time: %.3f'%(time.time()-time0))
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist,'center':cent,
            'size':sizexy,'target':target,'tilt':-twoth,'rotation':90.,'twoth':str(round(twoth,1))}
    data['pixLimit'] = 5
    data['calibdmin'] = 1.0
    data['cutoff'] = .5
    Npix = sizexy[0]*sizexy[1]
    
    return lines,data,Npix,image
        
