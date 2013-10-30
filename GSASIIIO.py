# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIIO: Misc I/O routines*
=============================

Module with miscellaneous routines for input and output. Many
are GUI routines to interact with user.

Includes support for image reading.

Also includes base classes for data import routines.

'''
"""GSASIIIO: functions for IO of data
   Copyright: 2008, Robert B. Von Dreele (Argonne National Laboratory)
"""
import wx
import math
import numpy as np
import cPickle
import sys
import re
import random as ran
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIgrid as G2gd
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import GSASIIpwdGUI as G2pdG
import GSASIIElem as G2el
import GSASIIstrIO as G2stIO
import GSASIImapvars as G2mv
import os
import os.path as ospath

def sfloat(S):
    'Convert a string to float. An empty field is treated as zero'
    if S.strip():
        return float(S)
    else:
        return 0.0

def sint(S):
    'Convert a string to int. An empty field is treated as zero'
    if S.strip():
        return int(S)
    else:
        return 0

def trim(val):
    '''Simplify a string containing leading and trailing spaces
    as well as newlines, tabs, repeated spaces etc. into a shorter and
    more simple string, by replacing all ranges of whitespace
    characters with a single space. 

    :param str val: the string to be simplified

    :returns: the (usually) shortened version of the string
    '''
    return re.sub('\s+', ' ', val).strip()

def makeInstDict(names,data,codes):
    inst = dict(zip(names,zip(data,data,codes)))
    for item in inst:
        inst[item] = list(inst[item])
    return inst

def FileDlgFixExt(dlg,file):
    'this is needed to fix a problem in linux wx.FileDialog'
    ext = dlg.GetWildcard().split('|')[2*dlg.GetFilterIndex()+1].strip('*')
    if ext not in file:
        file += ext
    return file
        
def GetPowderPeaks(fileName):
    'Read powder peaks from a file'
    sind = lambda x: math.sin(x*math.pi/180.)
    asind = lambda x: 180.*math.asin(x)/math.pi
    Cuka = 1.54052
    File = open(fileName,'Ur')
    Comments = []
    peaks = []
    S = File.readline()
    while S:
        if S[:1] == '#':
            Comments.append(S[:-1])
        else:
            item = S.split()
            if len(item) == 1:
                peaks.append([float(item[0]),1.0])
            elif len(item) > 1:
                peaks.append([float(item[0]),float(item[0])])
        S = File.readline()
    File.close()
    if Comments:
       print 'Comments on file:'
       for Comment in Comments: print Comment
    Peaks = []
    if peaks[0][0] > peaks[-1][0]:          # d-spacings - assume CuKa
        for peak in peaks:
            dsp = peak[0]
            sth = Cuka/(2.0*dsp)
            if sth < 1.0:
                tth = 2.0*asind(sth)
            else:
                break
            Peaks.append([tth,peak[1],True,False,0,0,0,dsp,0.0])
    else:                                   #2-thetas - assume Cuka (for now)
        for peak in peaks:
            tth = peak[0]
            dsp = Cuka/(2.0*sind(tth/2.0))
            Peaks.append([tth,peak[1],True,False,0,0,0,dsp,0.0])
    return Comments,Peaks

def CheckImageFile(G2frame,imagefile):
    '''Get an new image file name if the specified one does not
    exist

    :param wx.Frame G2frame: main GSAS-II Frame and data object

    :param str imagefile: name of image file

    :returns: imagefile, if it exists, or the name of a file
      that does exist or False if the user presses Cancel

    '''
    if not ospath.exists(imagefile):
        dlg = wx.FileDialog(G2frame, 'Bad image file name; choose name', '.', '',\
        'Any image file (*.edf;*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img)\
        |*.edf;*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img|\
        European detector file (*.edf)|*.edf|\
        Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|\
        MAR file (*.mar*)|*.mar*|\
        GE Image (*.avg;*.sum)|*.avg;*.sum|\
        ADSC Image (*.img)|*.img|\
        All files (*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            dlg.SetFilename(''+ospath.split(imagefile)[1])
            if dlg.ShowModal() == wx.ID_OK:
                imagefile = dlg.GetPath()
            else:
                imagefile = False
        finally:
            dlg.Destroy()
    return imagefile
        
def GetImageData(G2frame,imagefile,imageOnly=False):
    '''Read an image with the file reader keyed by the
    file extension

    :param wx.Frame G2frame: main GSAS-II Frame and data object. Note: not used!

    :param str imagefile: name of image file

    :param bool imageOnly: If True return only the image,
      otherwise  (default) return more (see below)

    :returns: an image as a numpy array or a list of four items:
      Comments, Data, Npix and the Image, as selected by imageOnly

    '''
    ext = ospath.splitext(imagefile)[1]
    Comments = []
    if ext == '.tif' or ext == '.tiff':
        Comments,Data,Npix,Image = GetTifData(imagefile)
    elif ext == '.edf':
        Comments,Data,Npix,Image = GetEdfData(imagefile)
    elif ext == '.img':
        Comments,Data,Npix,Image = GetImgData(imagefile)
        Image[0][0] = 0
    elif ext == '.mar3450' or ext == '.mar2300':
        Comments,Data,Npix,Image = GetMAR345Data(imagefile)
    elif ext in ['.sum','.avg','']:
        Comments,Data,Npix,Image = GetGEsumData(imagefile)
    elif ext == '.G2img':
        Comments,Data,Npix,Image = GetG2Image(imagefile)
    if imageOnly:
        return Image
    else:
        return Comments,Data,Npix,Image
        
def PutG2Image(filename,Comments,Data,Npix,image):
    'Write an image as a python pickle'
    File = open(filename,'wb')
    cPickle.dump([Comments,Data,Npix,image],File,1)
    File.close()
    return
    
def GetG2Image(filename):
    'Read an image as a python pickle'
    File = open(filename,'rb')
    Comments,Data,Npix,image = cPickle.load(File)
    File.close()
    return Comments,Data,Npix,image
    
def GetEdfData(filename,imageOnly=False):    
    'Read European detector data edf file'
    import struct as st
    import array as ar
    if not imageOnly:
        print 'Read European detector data edf file: ',filename
    File = open(filename,'rb')
    fileSize = os.stat(filename).st_size
    head = File.read(3072)
    lines = head.split('\n')
    sizexy = [0,0]
    pixSize = [0,0]
    cent = [0,0]
    head = ['European detector data',]
    for line in lines:
        fields = line.split()
        if 'Dim_1' in line:
            sizexy[0] = int(fields[2])
        elif 'Dim_2' in line:
            sizexy[1] = int(fields[2])
        elif 'DataType' in line:
            dType = fields[2]
        elif 'refined_wavelength' in line:
            wave = float(fields[2])
        elif 'Size' in line:
            imSize = int(fields[2])
        elif 'DataType' in lines:
            dType = fields[2]
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
    File.seek(fileSize-imSize)
    if dType == 'UnsignedShort':        
        image = np.array(ar.array('H',File.read(imSize)),dtype=np.int32)
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':pixSize,'wavelength':wave,'distance':dist,'center':cent,'size':sizexy}
    Npix = sizexy[0]*sizexy[1]
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
        
def GetGEsumData(filename,imageOnly=False):
    'Read SUM file as produced at 1-ID from G.E. images'
    import struct as st
    import array as ar
    if not imageOnly:
        print 'Read GE sum file: ',filename    
    File = open(filename,'rb')
    if '.sum' in filename:
        head = ['GE detector sum data from APS 1-ID',]
        sizexy = [2048,2048]
    elif '.avg' in filename:
        head = ['GE detector avg data from APS 1-ID',]
        sizexy = [2048,2048]
    else:
        head = ['GE detector raw data from APS 1-ID',]
        File.seek(18)
        size,nframes = st.unpack('<ih',File.read(6))
        sizexy = [2048,2048]
        pos = 8192
        File.seek(pos)
    Npix = sizexy[0]*sizexy[1]
    if '.sum' in filename:
        image = np.array(ar.array('f',File.read(4*Npix)),dtype=np.int32)
    elif '.avg' in filename:
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    else:
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
        while nframes > 1:
            image += np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
            nframes -= 1
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    data = {'pixelSize':(200,200),'wavelength':0.15,'distance':250.0,'center':[204.8,204.8],'size':sizexy}  
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
        
def GetImgData(filename,imageOnly=False):
    'Read an ADSC image file'
    import struct as st
    import array as ar
    if not imageOnly:
        print 'Read ADSC img file: ',filename
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
    data = {'pixelSize':pixel,'wavelength':wave,'distance':distance,'center':center,'size':[size,size]}
    image = []
    row = 0
    pos = 512
    File.seek(pos)
    image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    image = np.reshape(image,(sizexy[1],sizexy[0]))
#    image = np.zeros(shape=(size,size),dtype=np.int32)    
#    while row < size:
#        File.seek(pos)
#        line = ar.array('H',File.read(2*size))
#        image[row] = np.asarray(line)
#        row += 1
#        pos += 2*size
    File.close()
    if imageOnly:
        return image
    else:
        return lines[1:-2],data,Npix,image
       
def GetMAR345Data(filename,imageOnly=False):
    'Read a MAR-345 image plate image'
    import array as ar
    import struct as st
    try:
        import pack_f as pf
    except:
        msg = wx.MessageDialog(None, message="Unable to load the GSAS MAR image decompression, pack_f",
                               caption="Import Error",
                               style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
        msg.ShowModal()
        return None,None,None,None

    if not imageOnly:
        print 'Read Mar345 file: ',filename
    File = open(filename,'rb')
    head = File.read(4095)
    numbers = st.unpack('<iiiiiiiiii',head[:40])
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
        elif 'CENTER' in line:
            values = line.split()
            center = [float(values[2])/10.,float(values[4])/10.]    #make in mm from pixels
        if line: 
            head.append(line)
    data = {'pixelSize':pixel,'wavelength':wave,'distance':distance,'center':center}
    for line in head:
        if 'FORMAT' in line[0:6]:
            items = line.split()
            size = int(items[1])
            Npix = size*size
    pos = 4096
    data['size'] = [size,size]
    File.seek(pos)
    line = File.read(8)
    while 'CCP4' not in line:       #get past overflow list for now
        line = File.read(8)
        pos += 8
    pos += 37
    File.seek(pos)
    raw = File.read()
    File.close()
    image = np.zeros(shape=(size,size),dtype=np.int32)
    image = np.flipud(pf.pack_f(len(raw),raw,size,image).T)  #transpose to get it right way around & flip
    if imageOnly:
        return image
    else:
        return head,data,Npix,image

def GetTifData(filename,imageOnly=False):
    '''Read an image in a pseudo-tif format,
    as produced by a wide variety of software, almost always
    incorrectly in some way. 
    '''
    import struct as st
    import array as ar
    import ReadMarCCDFrame as rmf
    File = open(filename,'rb')
    dataType = 5
    try:
        Meta = open(filename+'.metadata','Ur')
        head = Meta.readlines()
        for line in head:
            line = line.strip()
            if 'dataType=' in line:
                dataType = int(line.split('=')[1])
        Meta.close()
    except IOError:
        print 'no metadata file found - will try to read file anyway'
        head = ['no metadata file found',]
        
    tag = File.read(2)
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
    for ied in range(NED):
        Tag,Type = st.unpack(byteOrd+'Hh',File.read(4))
        nVal = st.unpack(byteOrd+'i',File.read(4))[0]
        if Type == 1:
            Value = st.unpack(byteOrd+nVal*'b',File.read(nVal))
        elif Type == 2:
            Value = st.unpack(byteOrd+'i',File.read(4))
        elif Type == 3:
            Value = st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
            x = st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
        elif Type == 4:
            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
        elif Type == 5:
            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
        elif Type == 11:
            Value = st.unpack(byteOrd+nVal*'f',File.read(nVal*4))
        IFD[Tag] = [Type,nVal,Value]
        #print Tag,IFD[Tag]
    sizexy = [IFD[256][2][0],IFD[257][2][0]]
    [nx,ny] = sizexy
    Npix = nx*ny
    if 34710 in IFD:
        if not imageOnly:
            print 'Read MAR CCD tiff file: ',filename
        marFrame = rmf.marFrame(File,byteOrd,IFD)
        image = np.flipud(np.array(np.asarray(marFrame.image),dtype=np.int32))
        tifType = marFrame.filetitle
        pixy = (marFrame.pixelsizeX/1000.0,marFrame.pixelsizeY/1000.0)
        head = marFrame.outputHead()
    elif 272 in IFD:
        ifd = IFD[272]
        File.seek(ifd[2][0])
        S = File.read(ifd[1])
        if 'PILATUS' in S:
            tifType = 'Pilatus'
            dataType = 0
            pixy = (172,172)
            File.seek(4096)
            if not imageOnly:
                print 'Read Pilatus tiff file: ',filename
            image = ar.array('L',File.read(4*Npix))
            image = np.array(np.asarray(image),dtype=np.int32)
        else:
            if IFD[258][2][0] == 16:
                print sizexy
                tifType = 'GE'
                pixy = (200,200)
                File.seek(8)
                if not imageOnly:
                    print 'Read GE-detector tiff file: ',filename
                image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
            
    elif 262 in IFD and IFD[262][2][0] > 4:
        tifType = 'DND'
        pixy = (158,158)
        File.seek(512)
        if not imageOnly:
            print 'Read DND SAX/WAX-detector tiff file: ',filename
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    elif sizexy == [1536,1536]:
        tifType = 'APS Gold'
        pixy = (150,150)
        File.seek(64)
        if not imageOnly:
            print 'Read Gold tiff file:',filename
        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
    elif sizexy == [2048,2048] or sizexy == [1024,1024] or sizexy == [3072,3072]:
        if IFD[273][2][0] == 8:
            if IFD[258][2][0] == 32:
                tifType = 'PE'
                pixy = (200,200)
                File.seek(8)
                if not imageOnly:
                    print 'Read APS PE-detector tiff file: ',filename
                if dataType == 5:
                    image = np.array(ar.array('f',File.read(4*Npix)),dtype=np.float32)
                else:
                    image = np.array(ar.array('I',File.read(4*Npix)),dtype=np.int32)
                
        elif IFD[273][2][0] == 4096:
            if sizexy[0] == 3072:
                pixy =  (73,73)
                tifType = 'MAR225'            
            else:
                pixy = (158,158)
                tifType = 'MAR325'            
            File.seek(4096)
            if not imageOnly:
                print 'Read MAR CCD tiff file: ',filename
            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
        elif IFD[273][2][0] == 512:
            tiftype = '11-ID-C'
            pixy = [200,200]
            File.seek(512)
            if not imageOnly:
                print 'Read 11-ID-C tiff file: ',filename
            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)            
    elif sizexy == [4096,4096]:
        if IFD[273][2][0] == 8:
            if IFD[258][2][0] == 16:
                tifType = 'scanCCD'
                pixy = (9,9)
                File.seek(8)
                if not imageOnly:
                    print 'Read APS scanCCD tiff file: ',filename
                image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
        elif IFD[273][2][0] == 4096:
            tifType = 'Rayonix'
            pixy = (73.242,73.242)
            File.seek(4096)
            if not imageOnly:
                print 'Read Rayonix MX300HE tiff file: ',filename
            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#    elif sizexy == [960,960]:
#        tiftype = 'PE-BE'
#        pixy = (200,200)
#        File.seek(8)
#        if not imageOnly:
#            print 'Read Gold tiff file:',filename
#        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
           
    else:
        lines = ['not a known detector tiff file',]
        return lines,0,0,0
        
    image = np.reshape(image,(sizexy[1],sizexy[0]))
    center = [pixy[0]*sizexy[0]/2000,pixy[1]*sizexy[1]/2000]
    data = {'pixelSize':pixy,'wavelength':0.10,'distance':100.0,'center':center,'size':sizexy}
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,Npix,image
    
#def GetTifData(filename,imageOnly=False):
#    import struct as st
#    import array as ar
#    File = open(filename,'rb')
#    dataType = 5
#    try:
#        Meta = open(filename+'.metadata','Ur')
#        head = Meta.readlines()
#        for line in head:
#            line = line.strip()
#            if 'dataType=' in line:
#                dataType = int(line.split('=')[1])
#        Meta.close()
#    except IOError:
#        print 'no metadata file found - will try to read file anyway'
#        head = ['no metadata file found',]
#        
#    tag = File.read(2)
#    byteOrd = '<'
#    if tag == 'II' and int(st.unpack('<h',File.read(2))[0]) == 42:     #little endian
#        IFD = int(st.unpack(byteOrd+'i',File.read(4))[0])
#    elif tag == 'MM' and int(st.unpack('>h',File.read(2))[0]) == 42:   #big endian
#        byteOrd = '>'
#        IFD = int(st.unpack(byteOrd+'i',File.read(4))[0])        
#    else:
#        lines = ['not a detector tiff file',]
#        return lines,0,0,0
#    File.seek(IFD)                                                  #get number of directory entries
#    NED = int(st.unpack(byteOrd+'h',File.read(2))[0])
#    IFD = {}
#    for ied in range(NED):
#        Tag,Type = st.unpack(byteOrd+'Hh',File.read(4))
#        nVal = st.unpack(byteOrd+'i',File.read(4))[0]
#        if Type == 1:
#            Value = st.unpack(byteOrd+nVal*'b',File.read(nVal))
#        elif Type == 2:
#            Value = st.unpack(byteOrd+'i',File.read(4))
#        elif Type == 3:
#            Value = st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
#            x = st.unpack(byteOrd+nVal*'h',File.read(nVal*2))
#        elif Type == 4:
#            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
#        elif Type == 5:
#            Value = st.unpack(byteOrd+nVal*'i',File.read(nVal*4))
#        elif Type == 11:
#            Value = st.unpack(byteOrd+nVal*'f',File.read(nVal*4))
#        IFD[Tag] = [Type,nVal,Value]
##        print Tag,IFD[Tag]
#    sizexy = [IFD[256][2][0],IFD[257][2][0]]
#    [nx,ny] = sizexy
#    Npix = nx*ny
#    if 272 in IFD:
#        ifd = IFD[272]
#        File.seek(ifd[2][0])
#        S = File.read(ifd[1])
#        if 'PILATUS' in S:
#            tifType = 'Pilatus'
#            dataType = 0
#            pixy = (172,172)
#            File.seek(4096)
#            if not imageOnly:
#                print 'Read Pilatus tiff file: ',filename
#            image = ar.array('L',File.read(4*Npix))
#            image = np.array(np.asarray(image),dtype=np.int32)
#    elif 262 in IFD and IFD[262][2][0] > 4:
#        tifType = 'DND'
#        pixy = (158,158)
#        File.seek(512)
#        if not imageOnly:
#            print 'Read DND SAX/WAX-detector tiff file: ',filename
#        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#    elif sizexy == [1536,1536]:
#        tifType = 'APS Gold'
#        pixy = (150,150)
#        File.seek(64)
#        if not imageOnly:
#            print 'Read Gold tiff file:',filename
#        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#    elif sizexy == [2048,2048] or sizexy == [1024,1024] or sizexy == [3072,3072]:
#        if IFD[273][2][0] == 8:
#            if IFD[258][2][0] == 32:
#                tifType = 'PE'
#                pixy = (200,200)
#                File.seek(8)
#                if not imageOnly:
#                    print 'Read APS PE-detector tiff file: ',filename
#                if dataType == 5:
#                    image = np.array(ar.array('f',File.read(4*Npix)),dtype=np.float32)
#                else:
#                    image = np.array(ar.array('I',File.read(4*Npix)),dtype=np.int32)
#        elif IFD[273][2][0] == 4096:
#            if sizexy[0] == 3072:
#                pixy =  (73,73)
#                tifType = 'MAR225'            
#            else:
#                pixy = (158,158)
#                tifType = 'MAR325'            
#            File.seek(4096)
#            if not imageOnly:
#                print 'Read MAR CCD tiff file: ',filename
#            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#        elif IFD[273][2][0] == 512:
#            tiftype = '11-ID-C'
#            pixy = [200,200]
#            File.seek(512)
#            if not imageOnly:
#                print 'Read 11-ID-C tiff file: ',filename
#            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)            
#    elif sizexy == [4096,4096]:
#        if IFD[273][2][0] == 8:
#            if IFD[258][2][0] == 16:
#                tifType = 'scanCCD'
#                pixy = (9,9)
#                File.seek(8)
#                if not imageOnly:
#                    print 'Read APS scanCCD tiff file: ',filename
#                image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#        elif IFD[273][2][0] == 4096:
#            tifType = 'Rayonix'
#            pixy = (73.242,73.242)
#            File.seek(4096)
#            if not imageOnly:
#                print 'Read Rayonix MX300HE tiff file: ',filename
#            image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
##    elif sizexy == [960,960]:
##        tiftype = 'PE-BE'
##        pixy = (200,200)
##        File.seek(8)
##        if not imageOnly:
##            print 'Read Gold tiff file:',filename
##        image = np.array(ar.array('H',File.read(2*Npix)),dtype=np.int32)
#           
#    else:
#        lines = ['not a known detector tiff file',]
#        return lines,0,0,0
#        
#    image = np.reshape(image,(sizexy[1],sizexy[0]))
#    center = [pixy[0]*sizexy[0]/2000,pixy[1]*sizexy[1]/2000]
#    data = {'pixelSize':pixy,'wavelength':0.10,'distance':100.0,'center':center,'size':sizexy}
#    File.close()    
#    if imageOnly:
#        return image
#    else:
#        return head,data,Npix,image
#    
def ProjFileOpen(G2frame):
    'Read a GSAS-II project file'
    file = open(G2frame.GSASprojectfile,'rb')
    print 'load from file: ',G2frame.GSASprojectfile
    G2frame.SetTitle("GSAS-II data tree: "+
                     os.path.split(G2frame.GSASprojectfile)[1])
    wx.BeginBusyCursor()
    try:
        while True:
            try:
                data = cPickle.load(file)
            except EOFError:
                break
            datum = data[0]
            
            Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text=datum[0])
            if 'PWDR' in datum[0]:                
                G2frame.PatternTree.SetItemPyData(Id,datum[1][:3])     #temp. trim off junk
            else:
                G2frame.PatternTree.SetItemPyData(Id,datum[1])
            for datus in data[1:]:
                sub = G2frame.PatternTree.AppendItem(Id,datus[0])
#patch
                if datus[0] == 'Instrument Parameters' and len(datus[1]) == 1:
                    if 'PWDR' in datum[0]:
                        datus[1] = [dict(zip(datus[1][3],zip(datus[1][0],datus[1][1],datus[1][2]))),{}]
                    else:
                        datus[1] = [dict(zip(datus[1][2],zip(datus[1][0],datus[1][1]))),{}]
                    for item in datus[1][0]:               #zip makes tuples - now make lists!
                        datus[1][0][item] = list(datus[1][0][item])
#end patch
                G2frame.PatternTree.SetItemPyData(sub,datus[1])
            if 'IMG' in datum[0]:                   #retrieve image default flag & data if set
                Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Image Controls'))
                if Data['setDefault']:
                    G2frame.imageDefault = Data                
        file.close()
        print 'project load successful'
        G2frame.NewPlot = True
    except:
        msg = wx.MessageDialog(G2frame,message="Error reading file "+
            str(G2frame.GSASprojectfile)+". This is not a GSAS-II .gpx file",
            caption="Load Error",style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
        msg.ShowModal()
    finally:
        wx.EndBusyCursor()
    
def ProjFileSave(G2frame):
    'Save a GSAS-II project file'
    if not G2frame.PatternTree.IsEmpty():
        file = open(G2frame.GSASprojectfile,'wb')
        print 'save to file: ',G2frame.GSASprojectfile
        wx.BeginBusyCursor()
        try:
            item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while item:
                data = []
                name = G2frame.PatternTree.GetItemText(item)
                data.append([name,G2frame.PatternTree.GetItemPyData(item)])
                item2, cookie2 = G2frame.PatternTree.GetFirstChild(item)
                while item2:
                    name = G2frame.PatternTree.GetItemText(item2)
                    data.append([name,G2frame.PatternTree.GetItemPyData(item2)])
                    item2, cookie2 = G2frame.PatternTree.GetNextChild(item, cookie2)                            
                item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)                            
                cPickle.dump(data,file,1)
            file.close()
        finally:
            wx.EndBusyCursor()
        print 'project save successful'

def SaveIntegration(G2frame,PickId,data):
    'Save image integration results as powder pattern(s)'
    azms = G2frame.Integrate[1]
    X = G2frame.Integrate[2][:-1]
    Xminmax = [X[0],X[-1]]
    N = len(X)
    Id = G2frame.PatternTree.GetItemParent(PickId)
    name = G2frame.PatternTree.GetItemText(Id)
    name = name.replace('IMG ','PWDR ')
    Comments = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Comments'))
    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
    codes = [0 for i in range(12)]
    LRazm = data['LRazimuth']
    Azms = []
    if data['fullIntegrate'] and data['outAzimuths'] == 1:
        Azms = [45.0,]                              #a poor man's average?
    else:
        for i,azm in enumerate(azms[:-1]):
            Azms.append((azms[i+1]+azm)/2.)
    for i,azm in enumerate(azms[:-1]):
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        Id = 0
        while item:
            Name = G2frame.PatternTree.GetItemText(item)
            if name == Name:
                Id = item
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        parms = ['PXC',data['wavelength'],0.0,0.99,1.0,-0.10,0.4,0.30,1.0,0.0001,Azms[i]]    #set polarization for synchrotron radiation!
        Y = G2frame.Integrate[0][i]
        W = 1./Y                    #probably not true
        Sample = G2pdG.SetDefaultSample()
        Sample['Gonio. radius'] = data['distance']
        Sample['Omega'] = data['GonioAngles'][0]
        Sample['Chi'] = data['GonioAngles'][1]
        Sample['Phi'] = data['GonioAngles'][2]
        if Id:
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Comments'),Comments)                    
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Limits'),[tuple(Xminmax),Xminmax])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Background'),[['chebyschev',1,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
            inst = [dict(zip(names,zip(parms,parms,codes))),{}]
            for item in inst[0]:
                inst[0][item] = list(inst[0][item])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'),inst)
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Peak List'),[])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Index Peak List'),[])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Unit Cells List'),[])             
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Reflection Lists'),{})             
        else:
            Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text=name+" Azm= %.2f"%(Azms[i]))
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',1,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
            inst = [dict(zip(names,zip(parms,parms,codes))),{}]
            for item in inst[0]:
                inst[0][item] = list(inst[0][item])
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Instrument Parameters'),inst)
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Peak List'),[])
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Index Peak List'),[])
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Unit Cells List'),[])
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
        G2frame.PatternTree.SetItemPyData(Id,[[''],[np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]])
    G2frame.PatternTree.SelectItem(Id)
    G2frame.PatternTree.Expand(Id)
    G2frame.PatternId = Id
            
def powderFxyeSave(G2frame,exports,powderfile):
    'Save a powder histogram as a GSAS FXYE file'
    head,tail = ospath.split(powderfile)
    name,ext = tail.split('.')
    for i,export in enumerate(exports):
        filename = ospath.join(head,name+'-%03d.'%(i)+ext)
        prmname = filename.strip(ext)+'prm'
        prm = open(prmname,'w')      #old style GSAS parm file
        PickId = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, export)
        Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame, \
            PickId, 'Instrument Parameters'))[0]
        prm.write( '            123456789012345678901234567890123456789012345678901234567890        '+'\n')
        prm.write( 'INS   BANK      1                                                               '+'\n')
        prm.write(('INS   HTYPE   %sR                                                              '+'\n')%(Inst['Type'][0]))
        if 'Lam1' in Inst:              #Ka1 & Ka2
            prm.write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   '+'\n')%(Inst['Lam1'][0],Inst['Lam2'][0]))
        elif 'Lam' in Inst:             #single wavelength
            prm.write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   '+'\n')%(Inst['Lam'][1],0.0))
        prm.write( 'INS  1 IRAD     0                                                               '+'\n')
        prm.write( 'INS  1I HEAD                                                                    '+'\n')
        prm.write( 'INS  1I ITYP    0    0.0000  180.0000         1                                 '+'\n')
        prm.write(('INS  1DETAZM%10.3f                                                          '+'\n')%(Inst['Azimuth'][0]))
        prm.write( 'INS  1PRCF1     3    8   0.00100                                                '+'\n')
        prm.write(('INS  1PRCF11     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(Inst['U'][1],Inst['V'][1],Inst['W'][1],0.0))
        prm.write(('INS  1PRCF12     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(Inst['X'][1],Inst['Y'][1],Inst['SH/L'][1]/2.,Inst['SH/L'][1]/2.))
        prm.close()
        file = open(filename,'w')
        print 'save powder pattern to file: ',filename
        x,y,w,yc,yb,yd = G2frame.PatternTree.GetItemPyData(PickId)[1]
        file.write(powderfile+'\n')
        file.write('Instrument parameter file:'+ospath.split(prmname)[1]+'\n')
        file.write('BANK 1 %d %d CONS %.2f %.2f 0 0 FXYE\n'%(len(x),len(x),\
            100.*x[0],100.*(x[1]-x[0])))
        s = list(np.sqrt(1./np.array(w)))        
        XYW = zip(x,y,s)
        for X,Y,S in XYW:
            file.write("%15.6g %15.6g %15.6g\n" % (100.*X,Y,max(S,1.0)))
        file.close()
        print 'powder pattern file '+filename+' written'
        
def powderXyeSave(G2frame,exports,powderfile):
    'Save a powder histogram as a Topas XYE file'
    head,tail = ospath.split(powderfile)
    name,ext = tail.split('.')
    for i,export in enumerate(exports):
        filename = ospath.join(head,name+'-%03d.'%(i)+ext)
        PickId = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, export)
        file = open(filename,'w')
        file.write('#%s\n'%(export))
        print 'save powder pattern to file: ',filename
        x,y,w,yc,yb,yd = G2frame.PatternTree.GetItemPyData(PickId)[1]
        s = list(np.sqrt(1./np.array(w)))        
        XYW = zip(x,y,s)
        for X,Y,W in XYW:
            file.write("%15.6g %15.6g %15.6g\n" % (X,Y,W))
        file.close()
        print 'powder pattern file '+filename+' written'
        
def PDFSave(G2frame,exports):
    'Save a PDF G(r) and S(Q) in column formats'
    for export in exports:
        PickId = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, export)
        SQname = 'S(Q)'+export[4:]
        GRname = 'G(R)'+export[4:]
        sqfilename = ospath.join(G2frame.dirname,export.replace(' ','_')[5:]+'.sq')
        grfilename = ospath.join(G2frame.dirname,export.replace(' ','_')[5:]+'.gr')
        sqId = G2gd.GetPatternTreeItemId(G2frame, PickId, SQname)
        grId = G2gd.GetPatternTreeItemId(G2frame, PickId, GRname)
        sqdata = np.array(G2frame.PatternTree.GetItemPyData(sqId)[1][:2]).T
        grdata = np.array(G2frame.PatternTree.GetItemPyData(grId)[1][:2]).T
        sqfile = open(sqfilename,'w')
        grfile = open(grfilename,'w')
        sqfile.write('#T S(Q) %s\n'%(export))
        grfile.write('#T G(R) %s\n'%(export))
        sqfile.write('#L Q     S(Q)\n')
        grfile.write('#L R     G(R)\n')
        for q,sq in sqdata:
            sqfile.write("%15.6g %15.6g\n" % (q,sq))
        sqfile.close()
        for r,gr in grdata:
            grfile.write("%15.6g %15.6g\n" % (r,gr))
        grfile.close()
    
def PeakListSave(G2frame,file,peaks):
    'Save powder peaks to a data file'
    print 'save peak list to file: ',G2frame.peaklistfile
    if not peaks:
        dlg = wx.MessageDialog(G2frame, 'No peaks!', 'Nothing to save!', wx.OK)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        return
    for peak in peaks:
        file.write("%10.4f %12.2f %10.3f %10.3f \n" % \
            (peak[0],peak[2],peak[4],peak[6]))
    print 'peak list saved'
              
def IndexPeakListSave(G2frame,peaks):
    'Save powder peaks from the indexing list'
    file = open(G2frame.peaklistfile,'wa')
    print 'save index peak list to file: ',G2frame.peaklistfile
    wx.BeginBusyCursor()
    try:
        if not peaks:
            dlg = wx.MessageDialog(G2frame, 'No peaks!', 'Nothing to save!', wx.OK)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            return
        for peak in peaks:
            file.write("%12.6f\n" % (peak[7]))
        file.close()
    finally:
        wx.EndBusyCursor()
    print 'index peak list saved'
    
def SetNewPhase(Name='New Phase',SGData=None,cell=None):
    '''Create a new phase with default values for various parameters

    :param str Name: Name for new Phase

    :param dict SGData: space group data from :func:`GSASIIspc:SpcGroup`;
      defaults to data for P 1

    :param list cell: unit cell parameter list; defaults to
      [1.0,1.0,1.0,90.,90,90.,1.]

    '''
    if SGData is None: SGData = G2spc.SpcGroup('P 1')[1]
    if cell is None: cell=[1.0,1.0,1.0,90.,90,90.,1.]
    phaseData = {
        'ranId':ran.randint(0,sys.maxint),
        'General':{
            'Name':Name,
            'Type':'nuclear',
            'SGData':SGData,
            'Cell':[False,]+cell,
            'Pawley dmin':1.0,
            'Data plot type':'None',
            'SH Texture':{
                'Order':0,
                'Model':'cylindrical',
                'Sample omega':[False,0.0],
                'Sample chi':[False,0.0],
                'Sample phi':[False,0.0],
                'SH Coeff':[False,{}],
                'SHShow':False,
                'PFhkl':[0,0,1],
                'PFxyz':[0,0,1],
                'PlotType':'Pole figure'}},
        'Atoms':[],
        'Drawing':{},
        'Histograms':{},
        'Pawley ref':[],
        'RBModels':{},
        }
    return phaseData
    
def ReadEXPPhase(G2frame,filename):
    '''Read a phase from a GSAS .EXP file.
    Called in the GSAS phase import routine (see imports/G2phase.py)
    '''
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    textureData = {'Order':0,'Model':'cylindrical','Sample omega':[False,0.0],
        'Sample chi':[False,0.0],'Sample phi':[False,0.0],'SH Coeff':[False,{}],
        'SHShow':False,'PFhkl':[0,0,1],'PFxyz':[0,0,1],'PlotType':'Pole figure'}
    shNcof = 0
    file = open(filename, 'Ur')
    S = 1
    Expr = [{},{},{},{},{},{},{},{},{}]
    while S:
        S = file.readline()
        if 'EXPR NPHAS' in S[:12]:
            Num = S[12:-1].count('0')
            NPhas = S[12:-1].split()
        if 'CRS' in S[:3]:
            N = int(S[3:4])-1
            Expr[N][S[:12]] = S[12:-1]
    file.close()
    PNames = []
    for n,N in enumerate(NPhas):
        if N != '0':
            result = n
            key = 'CRS'+str(n+1)+'    PNAM'
            PNames.append(Expr[n][key])
    if Num < 8:
        dlg = wx.SingleChoiceDialog(G2frame, 'Which phase to read?', 'Read phase data', PNames, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelection()
        finally:
            dlg.Destroy()        
    EXPphase = Expr[result]
    keyList = EXPphase.keys()
    keyList.sort()
    SGData = {}
    if NPhas[result] == '1':
        Ptype = 'nuclear'
    elif NPhas[result] in ['2','3']:
        Ptype = 'magnetic'
    elif NPhas[result] == '4':
        Ptype = 'macromolecular'
    elif NPhas[result] == '10':
        Ptype = 'Pawley'
    for key in keyList:
        if 'PNAM' in key:
           PhaseName = EXPphase[key].strip()
        elif 'ABC   ' in key:
            abc = [float(EXPphase[key][:10]),float(EXPphase[key][10:20]),float(EXPphase[key][20:30])]                        
        elif 'ANGLES' in key:
            angles = [float(EXPphase[key][:10]),float(EXPphase[key][10:20]),float(EXPphase[key][20:30])]                                                
        elif 'SG SYM' in key:
            SpGrp = EXPphase[key][:15].strip()
            E,SGData = G2spc.SpcGroup(SpGrp)
            if E:
                E,SGData = G2spc.SpcGroup('P 1') # unlikely to need this
        elif 'OD    ' in key:
            SHdata = EXPphase[key].split() # may not have all 9 values
            SHvals = 9*[0]
            for i in range(9):
                try:
                    float(SHdata[i])
                    SHvals[i] = SHdata[i]
                except:
                    pass
            textureData['Order'] = int(SHvals[0])
            textureData['Model'] = shModels[int(SHvals[2])]
            textureData['Sample omega'] = [False,float(SHvals[6])]
            textureData['Sample chi'] = [False,float(SHvals[7])]
            textureData['Sample phi'] = [False,float(SHvals[8])]
            shNcof = int(SHvals[1])
    Atoms = []
    if Ptype == 'nuclear':
        for key in keyList:
            if 'AT' in key:
                if key[11:] == 'A':
                    S = EXPphase[key]
                elif key[11:] == 'B':
                    S += EXPphase[key]
                    Atom = [S[50:58].strip(),S[:10].strip().capitalize(),'',
                        float(S[10:20]),float(S[20:30]),float(S[30:40]),
                        float(S[40:50]),'',int(S[60:62]),S[130:131]]
                    if Atom[9] == 'I':
                        Atom += [float(S[68:78]),0.,0.,0.,0.,0.,0.]
                    elif Atom[9] == 'A':
                        Atom += [0.0,float(S[68:78]),float(S[78:88]),
                            float(S[88:98]),float(S[98:108]),
                            float(S[108:118]),float(S[118:128])]
                    XYZ = Atom[3:6]
                    Atom[7],Atom[8] = G2spc.SytSym(XYZ,SGData)
                    Atom.append(ran.randint(0,sys.maxint))
                    Atoms.append(Atom)
    elif Ptype == 'macromolecular':
        for key in keyList:
            if 'AT' in key[6:8]:
                S = EXPphase[key]
                Atom = [S[56:60],S[50:54].strip().upper(),S[54:56],
                    S[46:51].strip(),S[:8].strip().capitalize(),'',
                    float(S[16:24]),float(S[24:32]),float(S[32:40]),
                    float(S[8:16]),'1',1,'I',float(S[40:46]),0,0,0,0,0,0]
                XYZ = Atom[6:9]
                Atom[10],Atom[11] = G2spc.SytSym(XYZ,SGData)
                Atom.append(ran.randint(0,sys.maxint))
                Atoms.append(Atom)
    Volume = G2lat.calc_V(G2lat.cell2A(abc+angles))
    if shNcof:
        shCoef = {}
        nRec = [i+1 for i in range((shNcof-1)/6+1)]
        for irec in nRec:
            ODkey = keyList[0][:6]+'OD'+'%3dA'%(irec)
            indx = EXPphase[ODkey].split()
            ODkey = ODkey[:-1]+'B'
            vals = EXPphase[ODkey].split()
            for i,val in enumerate(vals):
                key = 'C(%s,%s,%s)'%(indx[3*i],indx[3*i+1],indx[3*i+2])
                shCoef[key] = float(val)
        textureData['SH Coeff'] = [False,shCoef]
        
    Phase = SetNewPhase(Name=PhaseName,SGData=SGData,cell=abc+angles+[Volume,])
    general = Phase['General']
    general['Type'] = Ptype
    general['SH Texture'] = textureData
    Phase['Atoms'] = Atoms
    return Phase
       
def ReadPDBPhase(filename):
    '''Read a phase from a PDB file.
    Called in the PDB phase import routine (see imports/G2phase.py)
    '''
    EightPiSq = 8.*math.pi**2
    file = open(filename, 'Ur')
    Phase = {}
    Title = ''
    Compnd = ''
    Atoms = []
    A = np.zeros(shape=(3,3))
    S = file.readline()
    while S:
        Atom = []
        if 'TITLE' in S[:5]:
            Title = S[10:72].strip()
            S = file.readline()
        elif 'COMPND    ' in S[:10]:
            Compnd = S[10:72].strip()
            S = file.readline()
        elif 'CRYST' in S[:5]:
            abc = S[7:34].split()
            angles = S[34:55].split()
            cell=[float(abc[0]),float(abc[1]),float(abc[2]),
                float(angles[0]),float(angles[1]),float(angles[2])]
            Volume = G2lat.calc_V(G2lat.cell2A(cell))
            AA,AB = G2lat.cell2AB(cell)
            SpGrp = S[55:65]
            E,SGData = G2spc.SpcGroup(SpGrp)
            # space group processing failed, try to look up name in table
            if E:
                SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                if SpGrpNorm:
                    E,SGData = G2spc.SpcGroup(SpGrpNorm)
            while E:
                print G2spc.SGErrors(E)
                dlg = wx.TextEntryDialog(None,
                    SpGrp[:-1]+' is invalid \nN.B.: make sure spaces separate axial fields in symbol',
                    'ERROR in space group symbol','',style=wx.OK)
                if dlg.ShowModal() == wx.ID_OK:
                    SpGrp = dlg.GetValue()
                    E,SGData = G2spc.SpcGroup(SpGrp)
                else:
                    return None
                dlg.Destroy()                
            SGlines = G2spc.SGPrint(SGData)
            for line in SGlines: print line
            S = file.readline()
        elif 'SCALE' in S[:5]:
            V = (S[10:41].split())
            A[int(S[5])-1] = [float(V[0]),float(V[1]),float(V[2])]
            S = file.readline()
        elif 'ATOM' in S[:4] or 'HETATM' in S[:6]:
            XYZ = [float(S[31:39]),float(S[39:47]),float(S[47:55])]
            XYZ = np.inner(AB,XYZ)
            XYZ = np.where(abs(XYZ)<0.00001,0,XYZ)
            SytSym,Mult = G2spc.SytSym(XYZ,SGData)
            Uiso = float(S[61:67])/EightPiSq
            Type = S[12:14].lower()
            if Type[0] in '123456789':
                Type = Type[1:]
            Atom = [S[22:27].strip(),S[17:20].upper(),S[20:22],
                S[12:17].strip(),Type.strip().capitalize(),'',XYZ[0],XYZ[1],XYZ[2],
                float(S[55:61]),SytSym,Mult,'I',Uiso,0,0,0,0,0,0]
            S = file.readline()
            if 'ANISOU' in S[:6]:
                Uij = S[30:72].split()
                Uij = [float(Uij[0])/10000.,float(Uij[1])/10000.,float(Uij[2])/10000.,
                    float(Uij[3])/10000.,float(Uij[4])/10000.,float(Uij[5])/10000.]
                Atom = Atom[:14]+Uij
                Atom[12] = 'A'
                S = file.readline()
            Atom.append(ran.randint(0,sys.maxint))
            Atoms.append(Atom)
        else:           
            S = file.readline()
    file.close()
    if Title:
        PhaseName = Title
    elif Compnd:
        PhaseName = Compnd
    else:
        PhaseName = 'None'
    Phase = SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell+[Volume,])
    Phase['General']['Type'] = 'macromolecular'
    Phase['Atoms'] = Atoms
    
    return Phase

class MultipleChoicesDialog(wx.Dialog):
    '''A dialog that offers a series of choices, each with a
    title and a wx.Choice widget. Intended to be used Modally. 
    typical input:

        *  choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
        *  headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
        
    selections are placed in self.chosen when OK is pressed
    '''
    def __init__(self,choicelist,headinglist,
                 head='Select options',
                 title='Please select from options below',
                 parent=None):
        self.chosen = []
        wx.Dialog.__init__(
            self,parent,wx.ID_ANY,head, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((10,10),1)
        topLabl = wx.StaticText(panel,wx.ID_ANY,title)
        mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.CENTER,10)
        self.ChItems = []
        for choice,lbl in zip(choicelist,headinglist):
            mainSizer.Add((10,10),1)
            self.chosen.append(0)
            topLabl = wx.StaticText(panel,wx.ID_ANY,' '+lbl)
            mainSizer.Add(topLabl,0,wx.ALIGN_LEFT,10)
            self.ChItems.append(wx.Choice(self, wx.ID_ANY, (100, 50), choices = choice))
            mainSizer.Add(self.ChItems[-1],0,wx.ALIGN_CENTER,10)

        OkBtn = wx.Button(panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()
        self.Fit()
        
    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        # save the results from the choice widgets
        self.chosen = []
        for w in self.ChItems:
            self.chosen.append(w.GetSelection())
        self.EndModal(wx.ID_OK)              
            
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.chosen = []
        self.EndModal(wx.ID_CANCEL)              
            
def ExtractFileFromZip(filename, selection=None, confirmread=True,
                       confirmoverwrite=True, parent=None,
                       multipleselect=False):
    '''If the filename is a zip file, extract a file from that
    archive.

    :param list Selection: used to predefine the name of the file
      to be extracted. Filename case and zip directory name are
      ignored in selection; the first matching file is used.

    :param bool confirmread: if True asks the user to confirm before expanding
      the only file in a zip

    :param bool confirmoverwrite: if True asks the user to confirm
      before overwriting if the extracted file already exists

    :param bool multipleselect: if True allows more than one zip
      file to be extracted, a list of file(s) is returned.
      If only one file is present, do not ask which one, otherwise
      offer a list of choices (unless selection is used).
    
    :returns: the name of the file that has been created or a
      list of files (see multipleselect)

    If the file is not a zipfile, return the name of the input file.
    If the zipfile is empty or no file has been selected, return None
    '''
    import zipfile # do this now, since we can save startup time by doing this only on need
    import shutil
    zloc = os.path.split(filename)[0]
    if not zipfile.is_zipfile(filename):
        #print("not zip")
        return filename

    z = zipfile.ZipFile(filename,'r')
    zinfo = z.infolist()

    if len(zinfo) == 0:
        #print('Zip has no files!')
        zlist = [-1]
    if selection:
        choices = [os.path.split(i.filename)[1].lower() for i in zinfo]
        if selection.lower() in choices:
            zlist = [choices.index(selection.lower())]
        else:
            print('debug: file '+str(selection)+' was not found in '+str(filename))
            zlist = [-1]
    elif len(zinfo) == 1 and confirmread:
        result = wx.ID_NO
        dlg = wx.MessageDialog(
            parent,
            'Is file '+str(zinfo[0].filename)+
            ' what you want to extract from '+
            str(os.path.split(filename)[1])+'?',
            'Confirm file', 
            wx.YES_NO | wx.ICON_QUESTION)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_NO:
            zlist = [-1]
        else:
            zlist = [0]
    elif len(zinfo) == 1:
        zlist = [0]
    elif multipleselect:
        # select one or more from a from list
        choices = [i.filename for i in zinfo]
        dlg = wx.MultiChoiceDialog(parent,'Select file(s) to extract from zip file'+str(filename),
            'Choose file(s)',choices,wx.CHOICEDLG_STYLE,)
        if dlg.ShowModal() == wx.ID_OK:
            zlist = dlg.GetSelections()
        else:
            zlist = []
        dlg.Destroy()
    else:
        # select one from a from list
        choices = [i.filename for i in zinfo]
        dlg = wx.SingleChoiceDialog(parent,
            'Select file to extract from zip file'+str(filename),'Choose file',
            choices,)
        if dlg.ShowModal() == wx.ID_OK:
            zlist = [dlg.GetSelection()]
        else:
            zlist = [-1]
        dlg.Destroy()
        
    outlist = []
    for zindex in zlist:
        if zindex >= 0:
            efil = os.path.join(zloc, os.path.split(zinfo[zindex].filename)[1])
            if os.path.exists(efil) and confirmoverwrite:
                result = wx.ID_NO
                dlg = wx.MessageDialog(parent,
                    'File '+str(efil)+' already exists. OK to overwrite it?',
                    'Confirm overwrite',wx.YES_NO | wx.ICON_QUESTION)
                try:
                    result = dlg.ShowModal()
                finally:
                    dlg.Destroy()
                if result == wx.ID_NO:
                    zindex = -1
        if zindex >= 0:
            # extract the file to the current directory, regardless of it's original path
            #z.extract(zinfo[zindex],zloc)
            eloc,efil = os.path.split(zinfo[zindex].filename)
            outfile = os.path.join(zloc, efil)
            fpin = z.open(zinfo[zindex])
            fpout = file(outfile, "wb")
            shutil.copyfileobj(fpin, fpout)
            fpin.close()
            fpout.close()
            outlist.append(outfile)
    z.close()
    if multipleselect and len(outlist) >= 1:
        return outlist
    elif len(outlist) == 1:
        return outlist[0]
    else:
        return None

######################################################################
# base classes for reading various types of data files
#   not used directly, only by subclassing
######################################################################
E,SGData = G2spc.SpcGroup('P 1') # data structure for default space group
class ImportBaseclass(object):
    '''Defines a base class for the reading of input files (diffraction
    data, coordinates,...
    '''
    def __init__(self,
                 formatName,
                 longFormatName=None,
                 extensionlist=[],
                 strictExtension=False,
                 ):
        self.formatName = formatName # short string naming file type
        if longFormatName: # longer string naming file type
            self.longFormatName = longFormatName
        else:
            self.longFormatName = formatName
        # define extensions that are allowed for the file type
        # for windows, remove any extensions that are duplicate, as case is ignored
        if sys.platform == 'windows' and extensionlist:
            extensionlist = list(set([s.lower() for s in extensionlist]))
        self.extensionlist = extensionlist
        # If strictExtension is True, the file will not be read, unless
        # the extension matches one in the extensionlist
        self.strictExtension = strictExtension
        self.warnings = ''
        # used for readers that will use multiple passes to read
        # more than one data block
        self.repeat = False
        self.repeatcount = 0
        #print 'created',self.__class__

    def BlockSelector(self, ChoiceList, ParentFrame=None,
                      title='Select a block',
                      size=None, header='Block Selector',
                      useCancel=True):
        ''' Provide a wx dialog to select a block if the file contains more
        than one set of data and one must be selected
        '''
        if useCancel:
            dlg = wx.SingleChoiceDialog(
                ParentFrame,title, header,ChoiceList)
        else:
            dlg = wx.SingleChoiceDialog(
                ParentFrame,title, header,ChoiceList,
                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.OK|wx.CENTRE)
        if size: dlg.SetSize(size)
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            return sel
        else:
            return None
        dlg.Destroy()

    def MultipleBlockSelector(self, ChoiceList, ParentFrame=None,
        title='Select a block',size=None, header='Block Selector'):
        '''Provide a wx dialog to select a block of data if the
        file contains more than one set of data and one must be
        selected.

        :returns: a list of the selected blocks
        '''
        dlg = wx.MultiChoiceDialog(ParentFrame,title, header,ChoiceList+['Select all'],
            wx.CHOICEDLG_STYLE)
        dlg.CenterOnParent()
        if size: dlg.SetSize(size)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
        else:
            return []
        dlg.Destroy()
        selected = []
        if len(ChoiceList) in sel:
            return range(len(ChoiceList))
        else:
            return sel
        return selected

    def MultipleChoicesDialog(self, choicelist, headinglist, ParentFrame=None, **kwargs):
        '''A modal dialog that offers a series of choices, each with a title and a wx.Choice
        widget. Typical input:
        
           * choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
           
           * headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
           
        optional keyword parameters are: head (window title) and title
        returns a list of selected indicies for each choice (or None)
        '''
        result = None
        dlg = MultipleChoicesDialog(choicelist,headinglist,
            parent=ParentFrame, **kwargs)          
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            result = dlg.chosen
        dlg.Destroy()
        return result

    def ShowBusy(self):
        wx.BeginBusyCursor()

    def DoneBusy(self):
        wx.EndBusyCursor()
        
#    def Reader(self, filename, filepointer, ParentFrame=None, **unused):
#        '''This method must be supplied in the child class
#        it will read the file
#        '''
#        return True # if read OK
#        return False # if an error occurs

    def ExtensionValidator(self, filename):
        '''This methods checks if the file has the correct extension
        Return False if this filename will not be supported by this reader
        Return True if the extension matches the list supplied by the reader
        Return None if the reader allows un-registered extensions
        '''
        if filename:
            ext = os.path.splitext(filename)[1]
            if sys.platform == 'windows': ext = ext.lower()
            if ext in self.extensionlist: return True
            if self.strictExtension: return False
        return None

    def ContentsValidator(self, filepointer):
        '''This routine will attempt to determine if the file can be read
        with the current format.
        This will typically be overridden with a method that 
        takes a quick scan of [some of]
        the file contents to do a "sanity" check if the file
        appears to match the selected format. 
        Expected to be called via self.Validator()
        '''
        #filepointer.seek(0) # rewind the file pointer
        return True

class ImportPhase(ImportBaseclass):
    '''Defines a base class for the reading of files with coordinates
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):
        # call parent __init__
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)
        # define a default Phase structure
        self.Phase = SetNewPhase(Name='new phase',SGData=SGData)

    def PhaseSelector(self, ChoiceList, ParentFrame=None,
        title='Select a phase', size=None,header='Phase Selector'):
        ''' Provide a wx dialog to select a phase if the file contains more
        than one phase
        '''
        return self.BlockSelector(ChoiceList,ParentFrame,title,
            size,header)

class ImportStructFactor(ImportBaseclass):
    '''Defines a base class for the reading of files with tables
    of structure factors

    Note that the default controls are stored in self.Controls and the
    default instrument parameters are stored in self.Parameters.
    These can be changed, but any changes will be the defaults for all
    subsequent uses of the :class:`ImportStructFactor` derived classes
    until :meth:`InitControls` and :meth:`InitParameters` are
    called. Probably better to use :meth:`UpdateControls` and
    :meth:`UpdateControls` (adding new args if needed) to change
    values.
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)

        # define contents of Structure Factor entry
        self.InitParameters()
        self.InitControls()
        self.RefDict = {'RefList':[],'FF':[]}
        
    def InitControls(self):
        'initialize the controls structure'
        self.Controls = { # dictionary with plotting controls
            'Type' : 'Fosq',
            'ifFc' : False,    # 
            'HKLmax' : [None,None,None],
            'HKLmin' : [None,None,None],
            'FoMax' : None,   # maximum observed structure factor as Fo
            'Zone' : '001',
            'Layer' : 0,
            'Scale' : 1.0,
            'log-lin' : 'lin',
            }

    def InitParameters(self):
        'initialize the instrument parameters structure'
        Lambda = 0.70926
        HistType = 'SXC'
        self.Parameters = [{'Type':[HistType,HistType], # create the structure
                            'Lam':[Lambda,Lambda]
                            }, {}]

    def UpdateParameters(self,Type=None,Wave=None):
        'Revise the instrument parameters'
        if Type is not None:
            self.Parameters[0]['Type'] = [Type,Type]
        if Wave is not None:
            self.Parameters[0]['Lam'] = [Wave,Wave]
            
    def UpdateControls(self,Type='Fosq',FcalcPresent=False):
        '''Scan through the reflections to update the Controls dictionary
        '''
        self.Controls['Type'] = Type
        self.Controls['ifFc'] = FcalcPresent
        HKLmax = [None,None,None]
        HKLmin = [None,None,None]
        Fo2max = None
        for refl in self.RefDict['RefList']:
            HKL = refl[:3]
            if Fo2max is None:
                Fo2max = refl[8]
            else:
                Fo2max = max(Fo2max,refl[8])
            for i,hkl in enumerate(HKL):
                if HKLmax[i] is None:
                    HKLmax[i] = hkl
                    HKLmin[i] = hkl
                else:
                    HKLmax[i] = max(HKLmax[i],hkl)
                    HKLmin[i] = min(HKLmin[i],hkl)
        self.Controls['HKLmax'] = HKLmax
        self.Controls['HKLmin'] = HKLmin
        if Type ==  'Fosq':
            self.Controls['FoMax'] = np.sqrt(Fo2max)
        elif Type ==  'Fo':
            self.Controls['FoMax'] = Fo2max
        else:
            print "Unsupported Struct Fact type in ImportStructFactor.UpdateControls"
            raise Exception,"Unsupported Struct Fact type in ImportStructFactor.UpdateControls"

######################################################################
class ImportPowderData(ImportBaseclass):
    '''Defines a base class for the reading of files with powder data
    '''
    # define some default instrument parameter files
    # just like GSAS, sigh
    defaultIparm_lbl = []
    defaultIparms = []
    defaultIparm_lbl.append('CuKa lab data')
    defaultIparms.append({
        'INS   HTYPE ':'PXC ',
        'INS  1 ICONS':'  1.540500  1.544300       0.0         0       0.7    0       0.5   ',
        'INS  1PRCF1 ':'    3    8      0.01                                                ',
        'INS  1PRCF11':'   2.000000E+00  -2.000000E+00   5.000000E+00   0.000000E+00        ',
        'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.150000E-01   0.150000E-01        ',
        })
    defaultIparm_lbl.append('0.6A synch')
    defaultIparms.append({
        'INS   HTYPE ':'PXC ',
        'INS  1 ICONS':'  0.600000  0.000000       0.0         0      0.99    0       0.5   ',
        'INS  1PRCF1 ':'    3    8      0.01                                                ',
        'INS  1PRCF11':'   1.000000E+00  -1.000000E+00   0.300000E+00   0.000000E+00        ',
        'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.100000E-01   0.100000E-01        ',
        })
    defaultIparm_lbl.append('1.5A CW neutron data')
    defaultIparms.append({
        'INS   HTYPE ':'PNC',
        'INS  1 ICONS':'   1.54020   0.00000   0.04000         0',
        'INS  1PRCF1 ':'    3    8     0.005',
        'INS  1PRCF11':'   0.239700E+03  -0.298200E+03   0.180800E+03   0.000000E+00',
        'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.400000E-01   0.300000E-01',
        })
    defaultIparm_lbl.append('10m TOF backscattering bank')
    defaultIparms.append({
        'INS   HTYPE ':'PNT',
        'INS  1 ICONS':'   5000.00      0.00      0.00',
        'INS  1BNKPAR':'    1.0000   150.000',       
        'INS  1PRCF1 ':'    1    8   0.01000',
        'INS  1PRCF11':'   0.000000E+00   5.000000E+00   3.000000E-02   1.000000E-03',
        'INS  1PRCF12':'   0.000000E+00   4.000000E+01   0.000000E+00   0.000000E+00',        
        })
    defaultIparm_lbl.append('10m TOF 90deg bank')
    defaultIparms.append({
        'INS   HTYPE ':'PNT',
        'INS  1 ICONS':'   3500.00      0.00      0.00',
        'INS  1BNKPAR':'    1.0000    90.000',       
        'INS  1PRCF1 ':'    1    8   0.01000',
        'INS  1PRCF11':'   0.000000E+00   5.000000E+00   3.000000E-02   4.000000E-03',
        'INS  1PRCF12':'   0.000000E+00   8.000000E+01   0.000000E+00   0.000000E+00',        
        })
    defaultIparm_lbl.append('63m POWGEN 90deg bank')
    defaultIparms.append({
        'INS   HTYPE ':'PNT',
        'INS  1 ICONS':'  22585.80      0.00      0.00',
        'INS  1BNKPAR':'    1.0000    90.000',       
        'INS  1PRCF1 ':'    1    8   0.01000',
        'INS  1PRCF11':'   0.000000E+00   1.000000E+00   3.000000E-02   4.000000E-03',
        'INS  1PRCF12':'   0.000000E+00   8.000000E+01   0.000000E+00   0.000000E+00',        
        })
    def __init__(self,
                 formatName,
                 longFormatName=None,
                 extensionlist=[],
                 strictExtension=False,
                 ):
        ImportBaseclass.__init__(self,formatName,
                                            longFormatName,
                                            extensionlist,
                                            strictExtension)
        self.powderentry = ['',None,None] #  (filename,Pos,Bank)
        self.powderdata = [] # Powder dataset
        '''A powder data set is a list with items [x,y,w,yc,yb,yd]:
                np.array(x), # x-axis values
                np.array(y), # powder pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.array(yc), # calc. intensities (zero)
                np.array(yb), # calc. background (zero)
                np.array(yd), # obs-calc profiles
        '''                            
        self.comments = []
        self.idstring = ''
        self.Sample = G2pdG.SetDefaultSample()
        self.GSAS = None     # used in TOF
        self.clockWd = None  # used in TOF
        self.repeat_instparm = True # Should a parm file be
        #                             used for multiple histograms? 
        self.instparm = None # name hint 
        self.instfile = '' # full path name to instrument parameter file
        self.instbank = '' # inst parm bank number
        self.instmsg = ''  # a label that gets printed to show
                           # where instrument parameters are from
        self.numbanks = 1
        self.instdict = {} # place items here that will be transferred to the instrument parameters
######################################################################
class ExportBaseclass(object):
    '''Defines a base class for the exporting of GSAS-II results
    '''
    def __init__(self,
                 G2frame,
                 formatName,
                 extension,
                 longFormatName=None,
                 ):
        self.G2frame = G2frame
        self.formatName = formatName # short string naming file type
        self.extension = extension
        if longFormatName: # longer string naming file type
            self.longFormatName = longFormatName
        else:
            self.longFormatName = formatName
        self.OverallParms = {}
        self.Phases = {}
        self.Histograms = {}
        self.powderDict = {}
        self.xtalDict = {}
        self.parmDict = {}
        self.sigDict = {}
        # updated in InitExport:
        self.currentExportType = None # type of export that has been requested
        # updated in ExportSelect (when used):
        self.phasenam = None # a list of selected phases
        self.histnam = None # a list of selected histograms
        self.filename = None # name of file to be written
        
        # items that should be defined in a subclass of this class
        self.exporttype = []  # defines the type(s) of exports that the class can handle.
        # The following types are defined: 'project', "phase", "powder", "single"
        self.multiple = False # set as True if the class can export multiple phases or histograms
        # self.multiple is ignored for "project" exports

    def InitExport(self,event):
        '''Determines the type of menu that called the Exporter.
        '''
        if event:
            self.currentExportType = self.G2frame.ExportLookup.get(event.Id)

    def ExportSelect(self,AskFile=True):
        '''Selects histograms or phases when needed. Sets a default file name.

        :param bool AskFile: if AskFile is True (default) get the name of the file
          in a dialog
        :returns: True in case of an error
        '''
        
        if self.currentExportType == 'phase':
            if len(self.Phases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any phases.')
                return True
            elif len(self.Phases) == 1:
                self.phasenam = self.Phases.keys()
            elif self.multiple: 
                choices = sorted(self.Phases.keys())
                phasenum = G2gd.ItemSelector(choices,self.G2frame,multiple=True)
                if phasenum is None: return True
                self.phasenam = [choices[i] for i in phasenum]
                if not self.phasenam: return True
            else:
                choices = sorted(self.Phases.keys())
                phasenum = G2gd.ItemSelector(choices,self.G2frame)
                if phasenum is None: return True
                self.phasenam = [choices[phasenum]]
        elif self.currentExportType == 'single':
            if len(self.xtalDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any single crystal data.')
                return True
            elif len(self.xtalDict) == 1:
                self.histnam = self.xtalDict.values()
            elif self.multiple:
                choices = sorted(self.xtalDict.values())
                hnum = G2gd.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
            else:
                choices = sorted(self.xtalDict.values())
                hnum = G2gd.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
        elif self.currentExportType == 'powder':
            if len(self.powderDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any powder data.')
                return True
            elif len(self.powderDict) == 1:
                self.histnam = self.powderDict.values()
            elif self.multiple:
                choices = sorted(self.powderDict.values())
                hnum = G2gd.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
            else:
                choices = sorted(self.powderDict.values())
                hnum = G2gd.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
        elif self.currentExportType == 'image':
            if len(self.Histograms) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any images.')
                return True
            elif len(self.Histograms) == 1:
                self.histnam = self.Histograms.keys()
            else:
                choices = sorted(self.Histograms.keys())
                hnum = G2gd.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not hnum: return True
                    self.histnam = [choices[i] for i in hnum]
                else:
                    if hnum is None: return True
                    self.histnam = [choices[hnum]]
        if self.currentExportType == 'map':
            # search for phases with maps
            mapPhases = []
            choices = []
            for phasenam in sorted(self.Phases):
                phasedict = self.Phases[phasenam] # pointer to current phase info            
                if len(phasedict['General']['Map'].get('rho',[])):
                    mapPhases.append(phasenam)
                    if phasedict['General']['Map'].get('Flip'):
                        choices.append('Charge flip map: '+str(phasenam))
                    elif phasedict['General']['Map'].get('MapType'):
                        choices.append(
                            str(phasedict['General']['Map'].get('MapType'))
                            + ' map: ' + str(phasenam))
                    else:
                        choices.append('unknown map: '+str(phasenam))
            # select a map if needed
            if len(mapPhases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any maps.')
                return True
            elif len(mapPhases) == 1:
                self.phasenam = mapPhases
            else: 
                phasenum = G2gd.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not phasenum: return True
                    self.phasenam = [mapPhases[i] for i in phasenum]
                else:
                    if phasenum is None: return True
                    self.phasenam = [mapPhases[phasenum]]

        if AskFile:
            self.filename = self.askSaveFile()
        else:
            self.filename = self.defaultSaveFile()
        if not self.filename: return True
        
    # def SetupExport(self,event,AskFile=True):
    #     '''Determines the type of menu that called the Exporter. Selects histograms
    #     or phases when needed. Better to replace with individual calls to 
    #     self.InitExport() and self.ExportSelect() so that the former can be called prior
    #     to self.LoadTree()

    #     :param bool AskFile: if AskFile is True (default) get the name of the file
    #       in a dialog
    #     :returns: True in case of an error
    #     '''
    #     self.ExportInit(event)
    #     return self.ExportSelect(AskFile)
                    
    def loadParmDict(self):
        '''Load the GSAS-II refinable parameters from the tree into a dict (self.parmDict). Update
        refined values to those from the last cycle and set the uncertainties for the
        refined parameters in another dict (self.sigDict).

        Expands the parm & sig dicts to include values derived from constraints.
        '''
        self.parmDict = {}
        self.sigDict = {}
        rigidbodyDict = {}
        covDict = {}
        consDict = {}
        Histograms,Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        if self.G2frame.PatternTree.IsEmpty(): return # nothing to do
        item, cookie = self.G2frame.PatternTree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.PatternTree.GetItemText(item)
            if name == 'Rigid bodies':
                 rigidbodyDict = self.G2frame.PatternTree.GetItemPyData(item)
            elif name == 'Covariance':
                 covDict = self.G2frame.PatternTree.GetItemPyData(item)
            elif name == 'Constraints':
                 consDict = self.G2frame.PatternTree.GetItemPyData(item)
            item, cookie = self.G2frame.PatternTree.GetNextChild(self.G2frame.root, cookie)
        rbVary,rbDict =  G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        self.parmDict.update(rbDict)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables =  G2stIO.GetPhaseData(
            Phases,RestraintDict=None,rbIds=rbIds,Print=False)
        self.parmDict.update(phaseDict)
        hapVary,hapDict,controlDict =  G2stIO.GetHistogramPhaseData(
            Phases,Histograms,Print=False,resetRefList=False)
        self.parmDict.update(hapDict)
        histVary,histDict,controlDict =  G2stIO.GetHistogramData(Histograms,Print=False)
        self.parmDict.update(histDict)
        self.parmDict.update(zip(
            covDict.get('varyList',[]),
            covDict.get('variables',[])))
        self.sigDict = dict(zip(
            covDict.get('varyList',[]),
            covDict.get('sig',[])))
        # expand to include constraints: first compile a list of constraints
        constList = []
        for item in consDict:
            constList += consDict[item]
        # now process the constraints
        G2mv.InitVars()
        constDict,fixedList,ignored = G2stIO.ProcessConstraints(constList)
        varyList = covDict.get('varyListStart')
        if varyList is None and len(constDict) == 0:
            # no constraints can use varyList
            varyList = covDict.get('varyList')
        elif varyList is None:
            # old GPX file from before pre-constraint varyList is saved
            print ' *** Old refinement: Please use Calculate/Refine to redo  ***'
            raise Exception(' *** CIF creation aborted ***')
        else:
            varyList = list(varyList)
        try:
            groups,parmlist = G2mv.GroupConstraints(constDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constDict,fixedList)
        except:
            # this really should not happen
            print ' *** ERROR - constraints are internally inconsistent ***'
            errmsg, warnmsg = G2mv.CheckConstraints(varyList,constDict,fixedList)
            print 'Errors',errmsg
            if warnmsg: print 'Warnings',warnmsg
            raise Exception(' *** CIF creation aborted ***')
        # add the constrained values to the parameter dictionary
        G2mv.Dict2Map(self.parmDict,varyList)
        # and add their uncertainties into the esd dictionary (sigDict)
        if covDict.get('covMatrix') is not None:
            self.sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],covDict['varyList'],self.parmDict))

    def loadTree(self):
        '''Load the contents of the data tree into a set of dicts
        (self.OverallParms, self.Phases and self.Histogram as well as self.powderDict
        & self.xtalDict)
        
        * The childrenless data tree items are overall parameters/controls for the
          entire project and are placed in self.OverallParms
        * Phase items are placed in self.Phases
        * Data items are placed in self.Histogram. The key for these data items
          begin with a keyword, such as PWDR, IMG, HKLF,... that identifies the data type.
        '''
        self.OverallParms = {}
        self.powderDict = {}
        self.xtalDict = {}
        if self.G2frame.PatternTree.IsEmpty(): return # nothing to do
        histType = None        
        if self.currentExportType == 'phase':
            # if exporting phases load them here
            sub = G2gd.GetPatternTreeItemId(self.G2frame,self.G2frame.root,'Phases')
            if not sub:
                print 'no phases found'
                return True
            item, cookie = self.G2frame.PatternTree.GetFirstChild(sub)
            while item:
                phaseName = self.G2frame.PatternTree.GetItemText(item)
                self.Phases[phaseName] =  self.G2frame.PatternTree.GetItemPyData(item)
                item, cookie = self.G2frame.PatternTree.GetNextChild(sub, cookie)
            return
        elif self.currentExportType == 'single':
            histType = 'HKLF'
        elif self.currentExportType == 'powder':
            histType = 'PWDR'
        elif self.currentExportType == 'image':
            histType = 'IMG'

        if histType: # Loading just one kind of tree entry
            item, cookie = self.G2frame.PatternTree.GetFirstChild(self.G2frame.root)
            while item:
                name = self.G2frame.PatternTree.GetItemText(item)
                if name.startswith(histType):
                    if self.Histograms.get(name): # there is already an item with this name
                        if name[-1] == '9':
                            name = name[:-1] + '10'
                        elif name[-1] in '012345678':
                            name = name[:-1] + str(int(name[:-1])+1)
                        else:                            
                            name += '-1'
                    self.Histograms[name] = {}
                    # the main info goes into Data, but the 0th
                    # element contains refinement results, carry
                    # that over too now. 
                    self.Histograms[name]['Data'] = self.G2frame.PatternTree.GetItemPyData(item)[1]
                    self.Histograms[name][0] = self.G2frame.PatternTree.GetItemPyData(item)[0]
                    item2, cookie2 = self.G2frame.PatternTree.GetFirstChild(item)
                    while item2: 
                        child = self.G2frame.PatternTree.GetItemText(item2)
                        self.Histograms[name][child] = self.G2frame.PatternTree.GetItemPyData(item2)
                        item2, cookie2 = self.G2frame.PatternTree.GetNextChild(item, cookie2)
                item, cookie = self.G2frame.PatternTree.GetNextChild(self.G2frame.root, cookie)
            # index powder and single crystal histograms by number
            for hist in self.Histograms:
                if hist.startswith("PWDR"): 
                    d = self.powderDict
                elif hist.startswith("HKLF"): 
                    d = self.xtalDict
                else:
                    return                    
                i = self.Histograms[hist].get('hId')
                if i is None and not d.keys():
                    i = 0
                elif i is None or i in d.keys():
                    i = max(d.keys())+1
                d[i] = hist
            return
        # else standard load: using all interlinked phases and histograms
        self.Histograms,self.Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        item, cookie = self.G2frame.PatternTree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.PatternTree.GetItemText(item)
            item2, cookie2 = self.G2frame.PatternTree.GetFirstChild(item)
            if not item2: 
                self.OverallParms[name] = self.G2frame.PatternTree.GetItemPyData(item)
            item, cookie = self.G2frame.PatternTree.GetNextChild(self.G2frame.root, cookie)
        # index powder and single crystal histograms
        for hist in self.Histograms:
            i = self.Histograms[hist]['hId']
            if hist.startswith("PWDR"): 
                self.powderDict[i] = hist
            elif hist.startswith("HKLF"): 
                self.xtalDict[i] = hist

    def dumpTree(self,mode='type'):
        '''Print out information on the data tree dicts loaded in loadTree
        '''
        print '\nOverall'
        if mode == 'type':
            def Show(arg): return type(arg)
        else:
            def Show(arg): return arg
        for key in self.OverallParms:
            print '  ',key,Show(self.OverallParms[key])
        print 'Phases'
        for key1 in self.Phases:
            print '    ',key1,Show(self.Phases[key1])
        print 'Histogram'
        for key1 in self.Histograms:
            print '    ',key1,Show(self.Histograms[key1])
            for key2 in self.Histograms[key1]:
                print '      ',key2,Show(self.Histograms[key1][key2])

    def defaultSaveFile(self):
        return os.path.abspath(
            os.path.splitext(self.G2frame.GSASprojectfile
                             )[0]+self.extension)
        
    def askSaveFile(self):
        '''Ask the user to supply a file name

        :returns: a file name (str)
        '''
        defnam = os.path.splitext(
            os.path.split(self.G2frame.GSASprojectfile)[1]
            )[0]+self.extension
        dlg = wx.FileDialog(
            self.G2frame, 'Input name for file to write', '.', defnam,
            self.longFormatName+' (*'+self.extension+')|*'+self.extension,
            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is correct
                filename = os.path.splitext(filename)[0]+self.extension
            else:
                filename = None
        finally:
            dlg.Destroy()
        return filename

    # Tools for file writing. 
    def OpenFile(self,fil=None):
        '''Open the output file

        :param str fil: The name of the file to open. If None (default)
          the name defaults to self.filename.
        :returns: the file object opened by the routine which is also
          saved as self.fp
        '''
        if not fil:
            fil = self.filename
        self.fp = open(fil,'w')
        return self.fp
    def Write(self,line):
        '''write a line of output, attaching a line-end character

        :param str line: the text to be written. 
        '''
        self.fp.write(line+'\n')
    def CloseFile(self,fp=None):
        '''Close a file opened in OpenFile

        :param file fp: the file object to be closed. If None (default)
          file object self.fp is closed. 
        '''
        if fp is None:
            fp = self.fp
            self.fp = None
        fp.close()
    # Tools to pull information out of the data arrays
    def GetCell(self,phasenam):
        """Gets the unit cell parameters and their s.u.'s for a selected phase

        :param str phasenam: the name for the selected phase
        :returns: `cellList,cellSig` where each is a 7 element list corresponding
          to a, b, c, alpha, beta, gamma, volume where `cellList` has the
          cell values and `cellSig` has their uncertainties.
        """
        phasedict = self.Phases[phasenam] # pointer to current phase info
        try:
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            cellSig = G2stIO.getCellEsd(pfx,
                                        phasedict['General']['SGData'],A,
                                        self.OverallParms['Covariance'])  # returns 7 vals, includes sigVol
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            return cellList,cellSig
        except KeyError:
            cell = phasedict['General']['Cell'][1:]
            return cell,7*[0]
    
    def GetAtoms(self,phasenam):
        """Gets the atoms associated with a phase. Can be used with standard
        or macromolecular phases

        :param str phasenam: the name for the selected phase
        :returns: a list of items for eac atom where each item is a list containing:
          label, typ, mult, xyz, and td, where

          * label and typ are the atom label and the scattering factor type (str)
          * mult is the site multiplicity (int)
          * xyz is contains a list with four pairs of numbers:
            x, y, z and fractional occupancy and
            their standard uncertainty (or a negative value)
          * td is contains a list with either one or six pairs of numbers:
            if one number it is U\ :sub:`iso` and with six numbers it is
            U\ :sub:`11`, U\ :sub:`22`, U\ :sub:`33`, U\ :sub:`12`, U\ :sub:`13` & U\ :sub:`23`
            paired with their standard uncertainty (or a negative value)
        """
        phasedict = self.Phases[phasenam] # pointer to current phase info            
        cx,ct,cs,cia = phasedict['General']['AtomPtrs']
        cfrac = cx+3
        fpfx = str(phasedict['pId'])+'::Afrac:'        
        atomslist = []
        for i,at in enumerate(phasedict['Atoms']):
            if phasedict['General']['Type'] == 'macromolecular':
                label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
            else:
                label = at[ct-1]
            fval = self.parmDict.get(fpfx+str(i),at[cfrac])
            fsig = self.sigDict.get(fpfx+str(i),-0.009)
            mult = at[cs+1]
            typ = at[ct]
            xyz = []
            for j,v in enumerate(('x','y','z')):
                val = at[cx+j]
                pfx = str(phasedict['pId'])+'::dA'+v+':'+str(i)
                sig = self.sigDict.get(pfx,-0.000009)
                xyz.append((val,sig))
            xyz.append((fval,fsig))
            td = []
            if at[cia] == 'I':
                pfx = str(phasedict['pId'])+'::AUiso:'+str(i)
                val = self.parmDict.get(pfx,at[cia+1])
                sig = self.sigDict.get(pfx,-0.0009)
                td.append((val,sig))
            else:
                for i,var in enumerate(('AU11','AU22','AU33','AU12','AU13','AU23')):
                    pfx = str(phasedict['pId'])+'::'+var+':'+str(i)
                    val = self.parmDict.get(pfx,at[cia+2+i])
                    sig = self.sigDict.get(pfx,-0.0009)
                    td.append((val,sig))
            atomslist.append((label,typ,mult,xyz,td))
        return atomslist
######################################################################

def ReadCIF(URLorFile):
    '''Open a CIF, which may be specified as a file name or as a URL using PyCifRW
    (from James Hester).
    The open routine gets confused with DOS names that begin with a letter and colon
    "C:\dir\" so this routine will try to open the passed name as a file and if that
    fails, try it as a URL

    :param str URLorFile: string containing a URL or a file name. Code will try first
      to open it as a file and then as a URL.

    :returns: a PyCifRW CIF object.
    '''
    import CifFile as cif # PyCifRW from James Hester

    # alternate approach:
    #import urllib
    #ciffile = 'file:'+urllib.pathname2url(filename)
   
    try:
        fp = open(URLorFile,'r')
        cf = cif.ReadCif(fp)
        fp.close()
        return cf
    except IOError:
        return cif.ReadCif(URLorFile)

if __name__ == '__main__':
    app = wx.PySimpleApp()
    frm = wx.Frame(None) # create a frame
    frm.Show(True)
    filename = '/tmp/notzip.zip'
    filename = '/tmp/all.zip'
    #filename = '/tmp/11bmb_7652.zip'
    
    #selection=None, confirmoverwrite=True, parent=None
    #print ExtractFileFromZip(filename, selection='11bmb_7652.fxye',parent=frm)
    print ExtractFileFromZip(filename,multipleselect=True)
                             #confirmread=False, confirmoverwrite=False)

    # choicelist=[ ('a','b','c'),
    #              ('test1','test2'),('no choice',)]
    # titles = [ 'a, b or c', 'tests', 'No option here']
    # dlg = MultipleChoicesDialog(
    #     choicelist,titles,
    #     parent=frm)
    # if dlg.ShowModal() == wx.ID_OK:
    #     print 'Got OK'
