# -*- coding: utf-8 -*-
"""GSASIIIO: functions for IO of data
   Copyright: 2008, Robert B. Von Dreele (Argonne National Laboratory)
"""
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import wx
import math
import numpy as np
import cPickle
import sys
import random as ran
import GSASIIpath
import GSASIIgrid as G2gd
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import GSASIIpwdGUI as G2pdG
import GSASIIElem as G2el
import os
import os.path as ospath

def sfloat(S):
    if S.strip():
        return float(S)
    else:
        return 0.0

def sint(S):
    if S.strip():
        return int(S)
    else:
        return 0

def FileDlgFixExt(dlg,file):
    #this is needed to fix a problem in linux wx.FileDialog
    ext = dlg.GetWildcard().split('|')[2*dlg.GetFilterIndex()+1].strip('*')
    if ext not in file:
        file += ext
    return file
    
# to be removed
def SelectPowderData(G2frame, filename):
    """Selects banks of data from a filename of any GSAS powder data format
    Input - filename: any GSAS powder data formatted file (currently STD, FXYE, FXY & ESD)
    Returns - a list of banks to be read; each entry in list is a tuple containing:
    filename: same as input filename
    Pos: position for start of data; record just after BANK record
    Bank: the BANK record
    """
    File = open(filename,'Ur')
    Title = '''
First line of this file:
'''+File.readline()
    dlg = wx.MessageDialog(G2frame, Title, 'Is this the file you want?', 
        wx.YES_NO | wx.ICON_QUESTION)
    try:
        result = dlg.ShowModal()
    finally:
        dlg.Destroy()
    if result == wx.ID_NO: return (0,0)
    Temperature = 300
    
    if '.xye' in filename:      #Topas style xye file (e.g. 2-th, I, sig) - no iparm file/no BANK record
        dlg = wx.MessageDialog(G2frame,'''Is this laboratory Cu Ka1/Ka2 data? 
(No = 0.6A wavelength synchrotron data)
Change wavelength in Instrument Parameters if needed''','Data type?',
            wx.YES_NO | wx.ICON_QUESTION)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        print result
        if result == wx.ID_YES:
            Iparm = {}                                               #Assume CuKa lab data
            Iparm['INS   HTYPE '] = 'PXC '
            Iparm['INS  1 ICONS'] = '  1.540500  1.544300       0.0         0       0.7    0       0.5   '
            Iparm['INS  1PRCF1 '] = '    3    8      0.01                                                '
            Iparm['INS  1PRCF11'] = '   2.000000E+00  -2.000000E+00   5.000000E+00   0.000000E+00        '
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.100000E-01   0.100000E-01        '
        else:
            Iparm = {}                                               #Assume 0.6A synchrotron data
            Iparm['INS   HTYPE '] = 'PXC '
            Iparm['INS  1 ICONS'] = '  0.600000  0.000000       0.0         0      0.99    0       0.5   '
            Iparm['INS  1PRCF1 '] = '    3    8      0.01                                                '
            Iparm['INS  1PRCF11'] = '   1.000000E+00  -1.000000E+00   0.300000E+00   0.000000E+00        '
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.100000E-02   0.100000E-02        '
                        
        
    else:                       #GSAS style fxye or fxy file (e.g. 100*2-th, I, sig)
        G2frame.IparmName = GetInstrumentFile(G2frame,filename)
        if G2frame.IparmName:
            Iparm = GetInstrumentData(G2frame.IparmName)
        else:
            Iparm = {}                                               #Assume CuKa lab data if no iparm file
            Iparm['INS   HTYPE '] = 'PXC '
            Iparm['INS  1 ICONS'] = '  1.540500  1.544300       0.0         0       0.7    0       0.5   '
            Iparm['INS  1PRCF1 '] = '    3    8      0.01                                                '
            Iparm['INS  1PRCF11'] = '   2.000000E+00  -2.000000E+00   5.000000E+00   0.000000E+00        '
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.100000E-01   0.100000E-01        '
    S = 1
    Banks = []
    Pos = []
    FoundData = []
    Comments = []
    wx.BeginBusyCursor()
    try:
        while S:
            S = File.readline()
            if S[:1] != '#':
                if S[:4] == 'BANK':
                    Banks.append(S)
                    Pos.append(File.tell())
                elif '.xye' in filename:    #No BANK in a xye file
                    Banks.append('BANK 1 XYE')
                    Pos.append(File.tell())
                    break
            else:
                Comments.append(S[:-1])
                if 'Temp' in S.split('=')[0]:
                    Temperature = float(S.split('=')[1])
        File.close()
    finally:
        wx.EndBusyCursor()
    if Comments:
       print 'Comments on file:'
       for Comment in Comments: print Comment
    if Banks:
        result = [0]
        if len(Banks) >= 2:
            dlg = wx.MultiChoiceDialog(G2frame, 'Which scans do you want?', 'Select scans', Banks, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                else:
                    result = []
            finally:
                dlg.Destroy()
        for i in result:
            FoundData.append((filename,Pos[i],Banks[i]))
    else:
        dlg = wx.MessageDialog(G2frame, 'ERROR - this is not a GSAS powder data file', 'No BANK records', wx.OK | wx.ICON_ERROR)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
    return FoundData,Iparm,Comments,Temperature

# to be removed
def GetInstrumentFile(G2frame,filename):
    import os.path as op
    dlg = wx.FileDialog(G2frame,'Choose an instrument file','.', '', 'GSAS iparm file (*.prm)|*.prm|All files(*.*)|*.*', 
        wx.OPEN|wx.CHANGE_DIR)
    Tname = filename[:filename.index('.')]+'.prm'
    if op.exists(Tname):
        G2frame.IparmName = Tname        
    if G2frame.IparmName: dlg.SetFilename(G2frame.IparmName)
    filename = ''
    try:
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
    finally:
        dlg.Destroy()
    return filename

# to be removed
def GetInstrumentData(IparmName):
    file = open(IparmName, 'Ur')
    S = 1
    Iparm = {}
    while S:
        S = file.readline()
        Iparm[S[:12]] = S[12:-1]
    return Iparm
    
def GetPowderPeaks(fileName):
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

def GetPawleyPeaks(filename):
    rt2ln2x2 = 2.35482
    File = open(filename,'Ur')
    PawleyPeaks = []
    S = File.readline()         #skip header
    S = File.readline()
    item = S.split()
    while S:
        h,k,l = int(item[0]),int(item[1]),int(item[2])
        mult = int(item[3])
        tth = float(item[5])
        sig = float(item[6])/rt2ln2x2
        Iobs = float(item[7])*mult
        PawleyPeaks.append([h,k,l,mult,tth,False,Iobs,0.0])
        S = File.readline()
        item = S.split()
        if item[3] == '-100.0000':       #find trailer
            break
    File.close()
    return PawleyPeaks
    
# this will be removed eventually
def GetHKLData(filename):
    print 'Reading: '+filename
    File = open(filename,'Ur')
    HKLref = []
    HKLmin = [1000,1000,1000]
    HKLmax = [0,0,0]
    FoMax = 0
    ifFc = False
    S = File.readline()
    while '#' in S[0]:        #get past comments if any
        S = File.readline()        
    if '_' in S:         #cif style .hkl file
        while 'loop_' not in S:         #skip preliminaries if any - can't have 'loop_' in them!
            S = File.readline()        
        S = File.readline()             #get past 'loop_' line
        pos = 0
        hpos = kpos = lpos = Fosqpos = Fcsqpos = sigpos = -1
        while S:
            if '_' in S:
                if 'index_h' in S:
                    hpos = pos
                elif 'index_k' in S:
                    kpos = pos
                elif 'index_l' in S:
                    lpos = pos
                elif 'F_squared_meas' in S:
                    Fosqpos = pos
                elif 'F_squared_calc' in S:
                    Fcsqpos = pos
                elif 'F_squared_sigma' in S:
                    sigpos = pos
                pos += 1
            else:
                data = S.split()
                if data:                    #avoid blank lines
                    HKL = np.array([int(data[hpos]),int(data[kpos]),int(data[lpos])])
                    h,k,l = HKL
                    Fosq = float(data[Fosqpos])
                    if sigpos != -1:
                        sigFosq = float(data[sigpos])
                    else:
                        sigFosq = 1.
                    if Fcsqpos != -1:
                        Fcsq = float(data[Fcsqpos])
                        if Fcsq:
                            ifFc = True
                    else:
                        Fcsq = 0.
                        
                    HKLmin = [min(h,HKLmin[0]),min(k,HKLmin[1]),min(l,HKLmin[2])]
                    HKLmax = [max(h,HKLmax[0]),max(k,HKLmax[1]),max(l,HKLmax[2])]
                    FoMax = max(FoMax,Fosq)
                    HKLref.append([HKL,Fosq,sigFosq,Fcsq,0,0,0])                 #room for Fcp, Fcpp & phase
            S = File.readline()
    else:                   #dumb h,k,l,Fo,sigFo .hkl file
        while S:
            h,k,l,Fo,sigFo = S.split()
            HKL = np.array([int(h),int(k),int(l)])
            h,k,l = HKL
            Fo = float(Fo)
            sigFo = float(sigFo)
            HKLmin = [min(h,HKLmin[0]),min(k,HKLmin[1]),min(l,HKLmin[2])]
            HKLmax = [max(h,HKLmax[0]),max(k,HKLmax[1]),max(l,HKLmax[2])]
            FoMax = max(FoMax,Fo)
            HKLref.append([HKL,Fo**2,2.*Fo*sigFo,0,0,0,0])                 #room for Fc, Fcp, Fcpp & phase
            S = File.readline()
    File.close()
    return HKLref,HKLmin,HKLmax,FoMax,ifFc

# to be removed
def GetPowderData(filename,Pos,Bank,DataType):
    '''Reads one BANK of data from GSAS raw powder data file
    input:
    filename: GSAS raw powder file dataname
    Pos: start of data in file just after BANK record
    Bank: the BANK record
    DataType: powder data type, e.g. "PXC" for Powder X-ray CW data
    returns: list [x,y,e,yc,yb]
    x: np.array of x-axis values
    y: np.array of powder pattern intensities
    w: np.array of w=sig(intensity)^2 values
    yc: np.array of calc. intensities (zero)
    yb: np.array of calc. background (zero)
    yd: np.array of obs-calc profiles
    '''
    print 'Reading: '+filename
    print 'Bank:    '+Bank[:-1]
    if 'FXYE' in Bank:
        return GetFXYEdata(filename,Pos,Bank,DataType)
    elif ' XYE' in Bank:
        return GetXYEdata(filename,Pos,Bank,DataType)
    elif 'FXY' in Bank:
        return GetFXYdata(filename,Pos,Bank,DataType)
    elif 'ESD' in Bank:
        return GetESDdata(filename,Pos,Bank,DataType)
    elif 'STD' in Bank:
        return GetSTDdata(filename,Pos,Bank,DataType)
    else:
        return GetSTDdata(filename,Pos,Bank,DataType)
    return []

# to be removed
def GetFXYEdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    File.seek(Pos)
    x = []
    y = []
    w = []
    S = File.readline()
    while S and S[:4] != 'BANK':
        vals = S.split()
        if DataType[2] == 'C':
            x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
        elif DataType[2] == 'T':
            x.append(float(vals[0])/1000.0)             #TOF: from musec to millisec
        f = float(vals[1])
        if f <= 0.0:
            y.append(0.0)
            w.append(1.0)
        else:
            y.append(float(vals[1]))
            w.append(1.0/float(vals[2])**2)
        S = File.readline()
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
    
# to be removed
def GetXYEdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    File.seek(Pos)
    x = []
    y = []
    w = []
    S = File.readline()
    while S:
        vals = S.split()
        try:
            x.append(float(vals[0]))
            f = float(vals[1])
            if f <= 0.0:
                y.append(0.0)
                w.append(1.0)
            else:
                y.append(float(vals[1]))
                w.append(1.0/float(vals[2])**2)
            S = File.readline()
        except ValueError:
            break
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
    
    
# to be removed
def GetFXYdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    File.seek(Pos)
    x = []
    y = []
    w = []
    S = File.readline()
    while S and S[:4] != 'BANK':
        vals = S.split()
        if DataType[2] == 'C':
            x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
        elif DataType[2] == 'T':
            x.append(float(vals[0])/1000.0)             #TOF: from musec to millisec
        f = float(vals[1])
        if f > 0.0:
            y.append(f)
            w.append(1.0/f)
        else:              
            y.append(0.0)
            w.append(1.0)
        S = File.readline()
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
    
# to be removed
def GetESDdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    cons = Bank.split()
    if DataType[2] == 'C':
        start = float(cons[5])/100.0               #CW: from centidegrees to degrees
        step = float(cons[6])/100.0
    elif DataType[2] == 'T':
        start = float(cons[5])/1000.0              #TOF: from musec to millisec
        step = float(cons[6])/1000.0
    File.seek(Pos)
    x = []
    y = []
    w = []
    S = File.readline()
    j = 0
    while S and S[:4] != 'BANK':
        for i in range(0,80,16):
            xi = start+step*j
            yi = sfloat(S[i:i+8])
            ei = sfloat(S[i+8:i+16])
            x.append(xi)
            if yi > 0.0:
                y.append(yi)
                w.append(1.0/ei**2)
            else:              
                y.append(0.0)
                w.append(1.0)
            j += 1
        S = File.readline()
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]

# to be removed
def GetSTDdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    cons = Bank.split()
    Nch = int(cons[2])
    if DataType[2] == 'C':
        start = float(cons[5])/100.0               #CW: from centidegrees to degrees
        step = float(cons[6])/100.0
    elif DataType[2] == 'T':
        start = float(cons[5])/1000.0              #TOF: from musec to millisec - not likely!
        step = float(cons[6])/1000.0
    File.seek(Pos)
    x = []
    y = []
    w = []
    S = File.readline()
    j = 0
    while S and S[:4] != 'BANK':
        for i in range(0,80,8):
            xi = start+step*j
            ni = max(sint(S[i:i+2]),1)
            yi = max(sfloat(S[i+2:i+8]),0.0)
            if yi:
                vi = yi/ni
            else:
                yi = 0.0
                vi = 1.0
            j += 1
            if j < Nch:
                x.append(xi)
                y.append(yi)
                w.append(1.0/vi)
        S = File.readline()
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
    
def CheckImageFile(G2frame,imagefile):
    if not ospath.exists(imagefile):
        dlg = wx.FileDialog(G2frame, 'Bad image file name; choose name', '.', '',\
        'Any image file (*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img)\
        |*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img|\
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
    ext = ospath.splitext(imagefile)[1]
    Comments = []
    if ext == '.tif' or ext == '.tiff':
        Comments,Data,Npix,Image = GetTifData(imagefile)
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
    File = open(filename,'wb')
    cPickle.dump([Comments,Data,Npix,image],File,1)
    File.close()
    return
    
def GetG2Image(filename):
    File = open(filename,'rb')
    Comments,Data,Npix,image = cPickle.load(File)
    File.close()
    return Comments,Data,Npix,image
    
def GetGEsumData(filename,imageOnly=False):
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
    image = pf.pack_f(len(raw),raw,size,image)
    if imageOnly:
        return image.T              #transpose to get it right way around
    else:
        return head,data,Npix,image.T
        
def GetTifData(filename,imageOnly=False):
    import struct as st
    import array as ar
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
#    for key in IFD:
#        print key,IFD[key]
    sizexy = [IFD[256][2][0],IFD[257][2][0]]
    [nx,ny] = sizexy
    Npix = nx*ny
    if 272 in IFD:
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
    elif sizexy == [2048,2048] or sizexy == [1024,1024]:
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
            tifType = 'MAR'
            pixy = (158,158)
            File.seek(4096)
            if not imageOnly:
                print 'Read MAR CCD tiff file: ',filename
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
    
def ProjFileOpen(G2frame):
    file = open(G2frame.GSASprojectfile,'rb')
    print 'load from file: ',G2frame.GSASprojectfile
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
                G2frame.PatternTree.SetItemPyData(sub,datus[1])
            if 'IMG' in datum[0]:                   #retreive image default flag & data if set
                Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Image Controls'))
                if Data['setDefault']:
                    G2frame.imageDefault = Data                
        file.close()
        
    finally:
        wx.EndBusyCursor()
    print 'project load successful'
    G2frame.NewPlot = True
    
def ProjFileSave(G2frame):
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
    azms = G2frame.Integrate[1]
    X = G2frame.Integrate[2][:-1]
    Xminmax = [X[0],X[-1]]
    N = len(X)
    Id = G2frame.PatternTree.GetItemParent(PickId)
    name = G2frame.PatternTree.GetItemText(Id)
    name = name.replace('IMG ','PWDR ')
    Comments = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Comments'))
    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
    codes = [0 for i in range(11)]
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
        if Id:
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Comments'),Comments)                    
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Limits'),[tuple(Xminmax),Xminmax])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Background'),[['chebyschev',1,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'),[tuple(parms),parms[:],codes,names])
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
            G2frame.PatternTree.SetItemPyData(G2frame.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(parms),parms[:],codes,names])
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
    head,tail = ospath.split(powderfile)
    name,ext = tail.split('.')
    for i,export in enumerate(exports):
        filename = ospath.join(head,name+'-%03d.'%(i)+ext)
        prmname = filename.strip(ext)+'prm'
        prm = open(prmname,'w')      #old style GSAS parm file
        PickId = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, export)
        Values,Names = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame, \
            PickId, 'Instrument Parameters'))[1::2]     #get values & names
        Inst = dict(zip(Names,Values))
        print Inst['Type']
        prm.write( '            123456789012345678901234567890123456789012345678901234567890        '+'\n')
        prm.write( 'INS   BANK      1                                                               '+'\n')
        prm.write(('INS   HTYPE   %sR                                                              '+'\n')%(Inst['Type']))
        if 'Lam1' in Inst:              #Ka1 & Ka2
            prm.write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   '+'\n')%(Inst['Lam1'],Inst['Lam2']))
        elif 'Lam' in Inst:             #single wavelength
            prm.write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   '+'\n')%(Inst['Lam'],0.0))
        prm.write( 'INS  1 IRAD     0                                                               '+'\n')
        prm.write( 'INS  1I HEAD                                                                    '+'\n')
        prm.write( 'INS  1I ITYP    0    0.0000  180.0000         1                                 '+'\n')
        prm.write(('INS  1DETAZM%10.3f                                                          '+'\n')%(Inst['Azimuth']))
        prm.write( 'INS  1PRCF1     3    8   0.00100                                                '+'\n')
        prm.write(('INS  1PRCF11     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(Inst['U'],Inst['V'],Inst['W'],0.0))
        prm.write(('INS  1PRCF12     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(Inst['X'],Inst['Y'],Inst['SH/L']/2.,Inst['SH/L']/2.))
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
    
def SetNewPhase(Name='New Phase',SGData=G2spc.SpcGroup('P 1')[1],cell=[1.0,1.0,1.0,90.,90,90.,1.]):
    phaseData = {
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
        'Models':{},
        }
    return phaseData
    
def ReadEXPPhase(G2frame,filename):
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
                    Atom = [S[50:58].strip(),S[:10].strip(),'',
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
                    S[46:51].strip(),S[:8].strip(),'',
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
            Type = S[12:14].upper()
            if Type[0] in '123456789':
                Type = Type[1:]
            Atom = [S[22:27].strip(),S[17:20].upper(),S[20:22],
                S[12:17].strip(),Type.strip(),'',XYZ[0],XYZ[1],XYZ[2],
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
    '''A dialog that offers a series of choices, each with a title and a wx.Choice
    widget. Intended to be used Modally. 
    typical input:
          choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
          headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
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
            
######################################################################
# base classes for reading various types of data files
#   not used directly, only by subclassing
######################################################################
E,SGData = G2spc.SpcGroup('P 1') # data structure for default space group
class ImportBaseclass(object):
    '''Defines a base class for the importing of data files (diffraction
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
        self.errors = ''
        # used for readers that will use multiple passes to read
        # more than one data block
        self.repeat = False
        self.repeatcount = 0
        #print 'created',self.__class__

    def BlockSelector(self, ChoiceList, ParentFrame=None,
                      title='Select a block',
                      size=None, header='Block Selector'):
        ''' Provide a wx dialog to select a block if the file contains more
        than one set of data and one must be selected
        '''
        dlg = wx.SingleChoiceDialog(
            ParentFrame,
            title, header,
            ChoiceList,
            )
        if size: dlg.SetSize(size)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            dlg.Destroy()
            return sel
        else:
            dlg.Destroy()
            return None

    def MultipleBlockSelector(self, ChoiceList, ParentFrame=None,
                      title='Select a block',
                      size=None, header='Block Selector'):
        ''' Provide a wx dialog to select a block of data if the file contains more
        than one set of data and one must be selected.
        Returns a list of the selected blocks
        '''
        dlg = wx.MultiChoiceDialog(
            ParentFrame,
            title, header,
            ChoiceList+['Select all'],
            wx.CHOICEDLG_STYLE
            )
        if size: dlg.SetSize(size)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return []
        selected = []
        if len(ChoiceList) in sel:
            return range(len(ChoiceList))
        else:
            return sel
        return selected

    def MultipleChoicesDialog(self, choicelist, headinglist, ParentFrame=None, **kwargs):
        '''A modal dialog that offers a series of choices, each with a title and a wx.Choice
        widget. 
        typical input:
           choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
           headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
        optional keyword parameters are: head (window title) and title
        returns a list of selected indicies for each choice (or None)
        '''
        result = None
        dlg = MultipleChoicesDialog(choicelist,headinglist,
                                    parent=ParentFrame, **kwargs)          
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
    def __init__(self,
                 formatName,
                 longFormatName=None,
                 extensionlist=[],
                 strictExtension=False,
                 ):
        # call parent __init__
        ImportBaseclass.__init__(self,formatName,
                                            longFormatName,
                                            extensionlist,
                                            strictExtension)
        # define a default Phase structure
        self.Phase = SetNewPhase(Name='new phase',SGData=SGData)

    def PhaseSelector(self, ChoiceList, ParentFrame=None,
                      title='Select a phase', size=None,
                      header='Phase Selector'):
        ''' Provide a wx dialog to select a phase if the file contains more
        than one phase
        '''
        return self.BlockSelector(ChoiceList,
                                  ParentFrame,
                                  title,
                                  size,
                                  header)

######################################################################
class ImportStructFactor(ImportBaseclass):
    '''Defines a base class for the reading of files with tables
    of structure factors
    '''
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

        # define contents of Structure Factor entry
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
        self.Parameters = [ # list with data collection parameters
            ('SXC',1.5428),
            ['SXC',1.5428],
            ['Type','Lam']
            ]
        self.RefList = []

    def UpdateParameters(self,Type=None,Wave=None):
        HistType = self.Parameters[0][0]
        HistWave = self.Parameters[0][1]
        if Type is not None:
            HistType = Type
        if Wave is not None:
            HistWave = Wave
        self.Parameters = [ # overwrite entire list 
            (HistType,HistWave),
            [HistType,HistWave],
            ['Type','Lam']
            ]
            
    def UpdateControls(self,Type='Fosq',FcalcPresent=False):
        '''Scan through the reflections to update the Controls dictionary
        '''
        self.Controls['Type'] = Type
        self.Controls['iffc'] = FcalcPresent
        HKLmax = [None,None,None]
        HKLmin = [None,None,None]
        Fo2max = None
        for HKL,Fo2,SFo2,Fc,Fcp,Fcpp,phase in self.RefList:
            if Fo2max is None:
                Fo2max = Fo2
            else:
                Fo2max = max(Fo2max,Fo2)
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
            print "Unsupported Stract Fact type in ImportStructFactor.UpdateControls"
            raise Exception,"Unsupported Stract Fact type in ImportStructFactor.UpdateControls"

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
        'INS  1PRCF1 ':'    3    8      0.01                                                ',
        'INS  1PRCF1 ':'    3    8     0.005',
        'INS  1PRCF11':'   0.239700E+03  -0.298200E+03   0.180800E+03   0.000000E+00',
        'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.400000E-01   0.300000E-01',
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
        self.instparm = None # name hint 
        self.instfile = '' # full path name to instrument parameter file
        self.instbank = '' # inst parm bank number
        self.instmsg = ''  # a label that gets printed to show
                           # where instrument parameters are from
        self.numbanks = 1
        self.instdict = {} # place items here that will be transferred to the instrument parameters

if __name__ == '__main__':
    app = wx.PySimpleApp()
    frm = wx.Frame(None) # create a frame
    choicelist=[ ('a','b','c'),
                 ('test1','test2'),('no choice',)]
    titles = [ 'a, b or c', 'tests', 'No option here']
    dlg = MultipleChoicesDialog(
        choicelist,titles,
        parent=frm)
    if dlg.ShowModal() == wx.ID_OK:
        print 'Got OK'
