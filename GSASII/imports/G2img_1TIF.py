# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-10-28 16:43:36 -0500 (Sat, 28 Oct 2023) $
# $Author: vondreele $
# $Revision: 5687 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_1TIF.py $
# $Id: G2img_1TIF.py 5687 2023-10-28 21:43:36Z vondreele $
########### SVN repository information ###################
'''
Note that the name ``G2img_1TIF`` is used so that this file will
sort to the top of the image formats and thus show up first in the menu.
(It is the most common, alas).
'''

from __future__ import division, print_function
import struct as st
import GSASIIobj as G2obj
import GSASIIpath
import GSASIIfiles as G2fil
import numpy as np
import time
DEBUG = False
GSASIIpath.SetVersionNumber("$Revision: 5687 $")
class TIF_ReaderClass(G2obj.ImportImage):
    '''Reads TIF files using a routine (:func:`GetTifData`) that looks
    for files that can be identified from known instruments and will
    correct for slightly incorrect TIF usage. 
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.tif','.tiff'),
            strictExtension=False,
            formatName = 'GSAS-II known TIF image',
            longFormatName = 'Various .tif and pseudo-TIF formats using GSAS-II reader'
            )
        self.scriptable = True

    def ContentsValidator(self, filename):
        '''Does the header match the required TIF header?
        '''
        return TIFValidator(filename)
    
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read the TIF file using :func:`GetTifData` which attempts to 
        recognize the detector type and set various parameters
        '''
        self.Npix = 0
        self.Comments,self.Data,self.Npix,self.Image = GetTifData(filename)
        if self.Npix == 0:
            return False
        self.LoadImage(ParentFrame,filename)
        return True
    
def GetTifData(filename):
    '''Read an image in a pseudo-tif format,
    as produced by a wide variety of software, almost always
    incorrectly in some way. 
    '''
    import struct as st
    import array as ar
    import ReadMarCCDFrame as rmf
    image = None
    File = open(filename,'rb')
    dataType = 5
    center = [None,None]
    wavelength = None
    distance = None
    polarization = None
    samplechangerpos = None
    xpixelsize = None
    ypixelsize = None
    try:
        Meta = open(filename+'.metadata','r')
        head = Meta.readlines()
        for line in head:
            line = line.strip()
            try:
                if '=' not in line: continue
                keyword = line.split('=')[0].strip()
                if 'dataType' == keyword:
                    dataType = int(line.split('=')[1])
                elif 'wavelength' == keyword.lower():
                    wavelength = float(line.split('=')[1])
                elif 'distance' == keyword.lower():
                    distance = float(line.split('=')[1])
                elif 'polarization' == keyword.lower():
                    polarization = float(line.split('=')[1])
                elif 'samplechangercoordinate' == keyword.lower():
                    samplechangerpos = float(line.split('=')[1])
                elif 'detectorxpixelsize' == keyword.lower():
                    xpixelsize = float(line.split('=')[1])
                elif 'detectorypixelsize' == keyword.lower():
                    ypixelsize = float(line.split('=')[1])
            except:
                G2fil.G2Print('error reading metadata: '+line)
        Meta.close()
    except IOError:
        if DEBUG:
            G2fil.G2Print ('no metadata file found - will try to read file anyway')
        head = ['no metadata file found',]
        
    tag = File.read(2)
    if 'bytes' in str(type(tag)):
        tag = tag.decode('latin-1')
    byteOrd = '<'
    if tag == 'II' and int(st.unpack('<h',File.read(2))[0]) == 42:     #little endian
        IFD = int(st.unpack(byteOrd+'i',File.read(4))[0])
    elif tag == 'MM' and int(st.unpack('>h',File.read(2))[0]) == 42:   #big endian
        byteOrd = '>'
        IFD = int(st.unpack(byteOrd+'i',File.read(4))[0])        
    else:
        lines = ['not a detector tiff file',]
        return lines,0,0,0
    File.seek(IFD)                                                  #get number of directory entries
    NED = int(st.unpack(byteOrd+'h',File.read(2))[0])
    IFD = {}
    nSlice = 1
    if DEBUG: print('byteorder:',byteOrd)
    for ied in range(NED):
        Tag,Type = st.unpack(byteOrd+'Hh',File.read(4))
        nVal = st.unpack(byteOrd+'i',File.read(4))[0]
        if DEBUG: print ('Try:',Tag,Type,nVal)
        if Type == 1:
            Value = st.unpack(byteOrd+nVal*'b',File.read(nVal))
        elif Type == 2:
            Value = st.unpack(byteOrd+'i',File.read(4))
        elif Type == 3:
            Value = st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
            st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
        elif Type == 4:
            if Tag in [273,279]:
                nSlice = nVal
                nVal = 1
            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
        elif Type == 5:
            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
        elif Type == 11:
            Value = st.unpack(byteOrd+nVal*'f',File.read(nVal*4))
        IFD[Tag] = [Type,nVal,Value]
        if DEBUG: print (Tag,IFD[Tag])
    sizexy = [IFD[256][2][0],IFD[257][2][0]]
    [nx,ny] = sizexy
    Npix = nx*ny
    time0 = time.time()
    if 34710 in IFD:
        G2fil.G2Print ('Read MAR CCD tiff file: '+filename)
        marFrame = rmf.marFrame(File,byteOrd,IFD)
        image = np.array(np.asarray(marFrame.image),dtype=np.int32)
        image = np.reshape(image,sizexy)
        if marFrame.origin:     #not upper left?
            image = np.flipud(image)
        if marFrame.viewDirection:  #view through sample to detector instead of TOWARD_SOURCE
            image = np.fliplr(image)
        tifType = marFrame.filetitle
        pixy = [marFrame.pixelsizeX/1000.0,marFrame.pixelsizeY/1000.0]
        head = marFrame.outputHead()
# extract resonable wavelength from header
        wavelength = marFrame.sourceWavelength*1e-5
        wavelength = (marFrame.opticsWavelength > 0) and marFrame.opticsWavelength*1e-5 or wavelength
        wavelength = (wavelength <= 0) and None or wavelength
# extract resonable distance from header
        distance = (marFrame.startXtalToDetector+marFrame.endXtalToDetector)*5e-4
        distance = (distance <= marFrame.startXtalToDetector*5e-4) and marFrame.xtalToDetector*1e-3 or distance
        distance = (distance <= 0) and None or distance
# extract resonable center from header
        center = [marFrame.beamX*marFrame.pixelsizeX*1e-9,marFrame.beamY*marFrame.pixelsizeY*1e-9]
        center = (center[0] != 0 and center[1] != 0) and center or [None,None]
#print head,tifType,pixy
    elif nSlice > 1:    #CheMin multislice tif file!
        try:
            import Image as Im
        except ImportError:
            try:
                from PIL import Image as Im
            except ImportError:
                G2fil.G2Print ("PIL/pillow Image module not present. This TIF cannot be read without this")
                #raise Exception("PIL/pillow Image module not found")
                lines = ['not a detector tiff file',]
                return lines,0,0,0
        tifType = 'CheMin'
        pixy = [40.,40.]
        image = np.flipud(np.array(Im.open(filename)))*10.
        distance = 18.0
        center = [pixy[0]*sizexy[0]/2000,0]     #the CheMin beam stop is here
        wavelength = 1.78892
    elif 272 in IFD:
        ifd = IFD[272]
        File.seek(ifd[2][0])
        S = File.read(ifd[1])
        if b'PILATUS' in S:
            tifType = 'Pilatus'
            dataType = 0
            pixy = [172.,172.]
            File.seek(4096)
            G2fil.G2Print ('Read Pilatus tiff file: '+filename)
            image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.int32),dtype=np.int32)
        else:
            if IFD[258][2][0] == 16:
                if sizexy == [3888,3072] or sizexy == [3072,3888]:
                    tifType = 'Dexela'
                    pixy = [74.8,74.8]
                    G2fil.G2Print ('Read Dexela detector tiff file: '+filename)
                else:
                    tifType = 'GE'
                    pixy = [200.,200.]
                    G2fil.G2Print ('Read GE-detector tiff file: '+filename)
                File.seek(8)
                image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
            elif IFD[258][2][0] == 32:
                # includes CHESS & Pilatus files from Area Detector
                tifType = 'CHESS'
                pixy = [200.,200.]
                File.seek(8)
                G2fil.G2Print ('Read as 32-bit unsigned (CHESS) tiff file: '+filename)
                image = np.array(ar.array('I',File.read(4*Npix)),dtype=np.uint32)
    elif 270 in IFD:
        File.seek(IFD[270][2][0])
        S = File.read(IFD[273][2][0]-IFD[270][2][0])
        if b'Pilatus3' in S:
            tifType = 'Pilatus3'
            dataType = 0
            pixy = [172.,172.]
            File.seek(IFD[273][2][0])
            G2fil.G2Print ('Read Pilatus3 tiff file: '+filename)
            image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.int32),dtype=np.int32)
        elif b'ImageJ' in S:
            tifType = 'ImageJ'
            dataType = 0
            pixy = [200.,200.]*IFD[277][2][0]
            File.seek(IFD[273][2][0])
            G2fil.G2Print ('Read ImageJ tiff file: '+filename)
            if IFD[258][2][0] == 32:
                image = File.read(4*Npix)
                image = np.array(np.frombuffer(image,dtype=byteOrd+'i4'),dtype=np.int32)
            elif IFD[258][2][0] == 16:
                image = File.read(2*Npix)
                if sizexy == [400,250]:
                    pixy = [175.2,175.2]        #for Sect36 ImageJ files
                else:
                    pixy = [109.92,109.92]      #for LCLS ImageJ tif files
                image = np.array(np.frombuffer(image,dtype=byteOrd+'u2'),dtype=np.int32)
        else:   #gain map from  11-ID-C?
            pixy = [200.,200.]
            tifType = 'Gain map'
            image = File.read(4*Npix)
            image = np.array(np.frombuffer(image,dtype=byteOrd+'f4')*1000,dtype=np.int32)
            
    elif 262 in IFD and IFD[262][2][0] > 4:
        tifType = 'DND'
        pixy = [158.,158.]
        File.seek(512)
        G2fil.G2Print ('Read DND SAX/WAX-detector tiff file: '+filename)
        image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
    elif sizexy == [1536,1536]:
        tifType = 'APS Gold'
        pixy = [150.,150.]
        File.seek(64)
        G2fil.G2Print ('Read Gold tiff file:'+filename)
        image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
    elif sizexy == [2048,2048] or sizexy == [1024,1024] or sizexy == [3072,3072]:
        if IFD[273][2][0] == 8:
            if IFD[258][2][0] == 32:
                tifType = 'PE'
                pixy = [200.,200.]
                File.seek(8)
                G2fil.G2Print ('Read APS PE-detector tiff file: '+filename)
                if dataType == 5:
                    image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.float32),dtype=np.int32)  #fastest
                else:
                    image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.int32),dtype=np.int32)
            elif IFD[258][2][0] == 16: 
                tifType = 'MedOptics D1'
                pixy = [46.9,46.9]
                File.seek(8)
                G2fil.G2Print ('Read MedOptics D1 tiff file: '+filename)
                image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
                  
        elif IFD[273][2][0] == 4096:
            if sizexy[0] == 3072:
                pixy =  [73.,73.]
                tifType = 'MAR225'            
            else:
                pixy = [158.,158.]
                tifType = 'MAR325'            
            File.seek(4096)
            G2fil.G2Print ('Read MAR CCD tiff file: '+filename)
            image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
        elif IFD[273][2][0] == 512:
            tifType = '11-ID-C'
            pixy = [200.,200.]
            File.seek(512)
            G2fil.G2Print ('Read 11-ID-C tiff file: '+filename)
            image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
                    
    elif sizexy == [4096,4096]:
        if IFD[273][2][0] == 8:
            if IFD[258][2][0] == 16:
                tifType = 'scanCCD'
                pixy = [9.,9.]
                File.seek(8)
                G2fil.G2Print ('Read APS scanCCD tiff file: '+filename)
                image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
            elif IFD[258][2][0] == 32:
                tifType = 'PE4k'
                pixy = [100.,100.]
                File.seek(8)
                G2fil.G2Print ('Read PE 4Kx4K tiff file: '+filename)
                image = np.array(np.frombuffer(File.read(4*Npix),dtype=np.float32)/2.**4,dtype=np.int32)
        elif IFD[273][2][0] == 4096:
            tifType = 'Rayonix'
            pixy = [73.242,73.242]
            File.seek(4096)
            G2fil.G2Print ('Read Rayonix MX300HE tiff file: '+filename)
            image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.uint16),dtype=np.int32)
    elif sizexy == [391,380]:
        pixy = [109.92,109.92]
        File.seek(8)
        image = np.array(np.frombuffer(File.read(2*Npix),dtype=np.int16),dtype=np.int32)
    elif sizexy == [380,391]:
        File.seek(110)
        pixy = [109.92,109.92]
        image = np.array(np.frombuffer(File.read(Npix),dtype=np.uint8),dtype=np.int32)
    elif sizexy ==  [825,830]:
        pixy = [109.92,109.92]
        File.seek(8)
        image = np.array(np.frombuffer(File.read(Npix),dtype=np.uint8),dtype=np.int32)
    elif sizexy ==  [1800,1800]:
        pixy = [109.92,109.92]
        File.seek(110)
        image = np.array(np.frombuffer(File.read(Npix),dtype=np.uint8),dtype=np.int32)
    elif sizexy == [2880,2880]:
        pixy = [150.,150.]
        File.seek(8)
        dt = np.dtype(np.float32)
        dt = dt.newbyteorder(byteOrd)
        image = np.array(np.frombuffer(File.read(Npix*4),dtype=dt),dtype=np.int32)
    elif sizexy == [3070,1102]:
        G2fil.G2Print ('Read Dectris Eiger 1M tiff file: '+filename)
        pixy = [75.,75.]
        File.seek(8)
        dt = np.dtype(np.float32)
        dt = dt.newbyteorder(byteOrd)
        image = np.array(np.frombuffer(File.read(Npix*4),dtype=np.uint32),dtype=np.int32)
    elif sizexy == [1024,402]:
        pixy = [56.,56.]
        File.seek(8)
        dt = np.dtype(np.float32)
        dt = dt.newbyteorder(byteOrd)
        image = np.array(np.frombuffer(File.read(Npix*2),dtype=np.uint16),dtype=np.int32)
        
#    elif sizexy == [960,960]:
#        tiftype = 'PE-BE'
#        pixy = (200,200)
#        File.seek(8)
#        if not imageOnly:
#            print 'Read Gold tiff file:',filename
#        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
           
    if image is None:
        lines = ['not a known detector tiff file',]
        File.close()    
        return lines,0,0,0
        
    if sizexy[1]*sizexy[0] != image.size: # test is resize is allowed
        lines = ['not a known detector tiff file',]
        File.close()    
        return lines,0,0,0
#    if GSASIIpath.GetConfigValue('debug'):
    if DEBUG:
        G2fil.G2Print ('image read time: %.3f'%(time.time()-time0))
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    center = (not center[0]) and [pixy[0]*sizexy[0]/2000,pixy[1]*sizexy[1]/2000] or center
    wavelength = (not wavelength) and 0.10 or wavelength
    distance = (not distance) and 100.0 or distance
    polarization = (not polarization) and 0.99 or polarization
    samplechangerpos = (not samplechangerpos) and 0.0 or samplechangerpos
    if xpixelsize is not None and ypixelsize is not None:
        pixy = [xpixelsize,ypixelsize]
        if GSASIIpath.GetConfigValue('debug'):
            G2fil.G2Print ('pixel size from metadata: '+str(pixy))
    data = {'pixelSize':pixy,'wavelength':wavelength,'distance':distance,'center':center,'size':sizexy,
            'setdist':distance,'PolaVal':[polarization,False],'samplechangerpos':samplechangerpos,'det2theta':0.0}
    File.close()    
    return head,data,Npix,image

def TIFValidator(filename):
    '''Does the header match the required TIF header?
    '''
    fp = open(filename,'rb')
    tag = fp.read(2)
    if 'bytes' in str(type(tag)):
        tag = tag.decode('latin-1')
    if tag == 'II' and int(st.unpack('<h',fp.read(2))[0]) == 42: #little endian
        pass
    elif tag == 'MM' and int(st.unpack('>h',fp.read(2))[0]) == 42: #big endian
        pass
    else:
        return False # header not found; not valid TIF
        fp.close()
    fp.close()
    return True
