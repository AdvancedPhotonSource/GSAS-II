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

def SelectPowderData(self, filename):
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
    dlg = wx.MessageDialog(self, Title, 'Is this the file you want?', 
        wx.YES_NO | wx.ICON_QUESTION)
    try:
        result = dlg.ShowModal()
    finally:
        dlg.Destroy()
    if result == wx.ID_NO: return (0,0)
    Temperature = 300
    
    if '.xye' in filename:      #Topas style xye file (e.g. 2-th, I, sig) - no iparm file/no BANK record
        dlg = wx.MessageDialog(self,'''Is this laboratory Cu Ka1/Ka2 data? 
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
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.150000E-01   0.150000E-01        '
        else:
            Iparm = {}                                               #Assume 0.6A synchrotron data
            Iparm['INS   HTYPE '] = 'PXC '
            Iparm['INS  1 ICONS'] = '  0.600000  0.000000       0.0         0      0.99    0       0.5   '
            Iparm['INS  1PRCF1 '] = '    3    8      0.01                                                '
            Iparm['INS  1PRCF11'] = '   1.000000E+00  -1.000000E+00   0.300000E+00   0.000000E+00        '
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.100000E-01   0.100000E-01        '
                        
        
    else:                       #GSAS style fxye or fxy file (e.g. 100*2-th, I, sig)
        self.IparmName = GetInstrumentFile(self,filename)
        if self.IparmName:
            Iparm = GetInstrumentData(self.IparmName)
        else:
            Iparm = {}                                               #Assume CuKa lab data if no iparm file
            Iparm['INS   HTYPE '] = 'PXC '
            Iparm['INS  1 ICONS'] = '  1.540500  1.544300       0.0         0       0.7    0       0.5   '
            Iparm['INS  1PRCF1 '] = '    3    8      0.01                                                '
            Iparm['INS  1PRCF11'] = '   2.000000E+00  -2.000000E+00   5.000000E+00   0.000000E+00        '
            Iparm['INS  1PRCF12'] = '   0.000000E+00   0.000000E+00   0.150000E-01   0.150000E-01        '
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
            dlg = wx.MultiChoiceDialog(self, 'Which scans do you want?', 'Select scans', Banks, wx.CHOICEDLG_STYLE)
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
        dlg = wx.MessageDialog(self, 'ERROR - this is not a GSAS powder data file', 'No BANK records', wx.OK | wx.ICON_ERROR)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
    return FoundData,Iparm,Comments,Temperature

def GetInstrumentFile(self,filename):
    import os.path as op
    dlg = wx.FileDialog(self,'Choose an instrument file','.', '', 'GSAS iparm file (*.prm)|*.prm|All files(*.*)|*.*', wx.OPEN)
    if self.dirname: 
        dlg.SetDirectory(self.dirname)
        Tname = filename[:filename.index('.')]+'.prm'
        if op.exists(Tname):
            self.IparmName = Tname        
    if self.IparmName: dlg.SetFilename(self.IparmName)
    filename = ''
    try:
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
    finally:
        dlg.Destroy()
    return filename

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

def GetSTDdata(filename,Pos,Bank,DataType):
    File = open(filename,'Ur')
    cons = Bank.split()
    Nch = cons[2]
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
                ei = math.sqrt(yi*ni)
            else:
                yi = 0.0
                ei = 1.0
            j += 1
            if j < Nch:
                x.append(xi)
                y.append(yi)
                w.append(1.0/ei**2)
        S = File.readline()
    File.close()
    N = len(x)
    return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
    
def CheckImageFile(self,imagefile):
    if not ospath.exists(imagefile):
        dlg = wx.FileDialog(self, 'Bad image file name; choose name', '.', '',\
        'Any image file (*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img)\
        |*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img|\
        Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|\
        MAR file (*.mar*)|*.mar*|\
        GE Image (*.avg;*.sum)|*.avg;*.sum|\
        ADSC Image (*.img)|*.img|\
        All files (*.*)|*.*',wx.OPEN)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            dlg.SetFilename(ospath.split(imagefile)[1])
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetDirectory()
                imagefile = dlg.GetPath()
            else:
                imagefile = False
        finally:
            dlg.Destroy()
    return imagefile
        
def GetImageData(self,imagefile,imageOnly=False):        
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
            if 'dataType' in line:
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
    
def ProjFileOpen(self):
    file = open(self.GSASprojectfile,'rb')
    print 'load from file: ',self.GSASprojectfile
    wx.BeginBusyCursor()
    try:
        while True:
            try:
                data = cPickle.load(file)
            except EOFError:
                break
            datum = data[0]
            
            Id = self.PatternTree.AppendItem(parent=self.root,text=datum[0])
            if 'PWDR' in datum[0]:                
                self.PatternTree.SetItemPyData(Id,datum[1][:3])     #temp. trim off junk
            else:
                self.PatternTree.SetItemPyData(Id,datum[1])
            for datus in data[1:]:
                sub = self.PatternTree.AppendItem(Id,datus[0])
                self.PatternTree.SetItemPyData(sub,datus[1])
            if 'IMG' in datum[0]:                   #retreive image default flag & data if set
                Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Image Controls'))
                if Data['setDefault']:
                    self.imageDefault = Data                
        file.close()
        
    finally:
        wx.EndBusyCursor()
    print 'project load successful'
    self.NewPlot = True
    
def ProjFileSave(self):
    if not self.PatternTree.IsEmpty():
        file = open(self.GSASprojectfile,'wb')
        print 'save to file: ',self.GSASprojectfile
        wx.BeginBusyCursor()
        try:
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                data = []
                name = self.PatternTree.GetItemText(item)
                data.append([name,self.PatternTree.GetItemPyData(item)])
                item2, cookie2 = self.PatternTree.GetFirstChild(item)
                while item2:
                    name = self.PatternTree.GetItemText(item2)
                    data.append([name,self.PatternTree.GetItemPyData(item2)])
                    item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                cPickle.dump(data,file,1)
            file.close()
        finally:
            wx.EndBusyCursor()
        print 'project save successful'
        
def SaveIntegration(self,PickId,data):
    azms = self.Integrate[1]
    X = self.Integrate[2][:-1]
    Xminmax = [X[0],X[-1]]
    N = len(X)
    Id = self.PatternTree.GetItemParent(PickId)
    name = self.PatternTree.GetItemText(Id)
    name = name.replace('IMG ','PWDR ')
    Comments = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,Id, 'Comments'))
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
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        Id = 0
        while item:
            Name = self.PatternTree.GetItemText(item)
            if name == Name:
                Id = item
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        parms = ['PXC',data['wavelength'],0.0,0.99,1.0,-0.10,0.4,0.30,1.0,0.0001,Azms[i]]    #set polarization for synchrotron radiation!
        Y = self.Integrate[0][i]
        W = 1./Y                    #probably not true
        Sample = G2pdG.SetDefaultSample()
        if Id:
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id, 'Comments'),Comments)                    
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Limits'),[tuple(Xminmax),Xminmax])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Instrument Parameters'),[tuple(parms),parms[:],codes,names])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Peak List'),[])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Index Peak List'),[])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Unit Cells List'),[])             
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'),{})             
        else:
            Id = self.PatternTree.AppendItem(parent=self.root,text=name+" Azm= %.2f"%(Azms[i]))
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(parms),parms[:],codes,names])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
        self.PatternTree.SetItemPyData(Id,[[''],[np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]])
    self.PatternTree.SelectItem(Id)
    self.PatternTree.Expand(Id)
    self.PatternId = Id
            
def powderFxyeSave(self,exports,powderfile):
    head,tail = ospath.split(powderfile)
    name,ext = tail.split('.')
    wx.BeginBusyCursor()
    for i,export in enumerate(exports):
        filename = ospath.join(head,name+'-%03d.'%(i)+ext)
        prmname = filename.strip(ext)+'prm'
        prm = open(prmname,'w')      #old style GSAS parm file
        PickId = G2gd.GetPatternTreeItemId(self, self.root, export)
        Values,Names = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self, \
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
        try:
            x,y,w,yc,yb,yd = self.PatternTree.GetItemPyData(PickId)[1]
            file.write(powderfile+'\n')
            file.write('BANK 1 %d %d CONS %.2f %.2f 0 0 FXYE\n'%(len(x),len(x),\
                100.*x[0],100.*(x[1]-x[0])))
            s = list(np.sqrt(1./np.array(w)))        
            XYW = zip(x,y,s)
            for X,Y,S in XYW:
                file.write("%15.6g %15.6g %15.6g\n" % (100.*X,Y,max(S,1.0)))
            file.close()
        finally:
            wx.EndBusyCursor()
        print 'powder pattern file written'
        
def powderXyeSave(self,exports,powderfile):
    head,tail = ospath.split(powderfile)
    name,ext = tail.split('.')
    for i,export in enumerate(exports):
        filename = ospath.join(head,name+'-%03d.'%(i)+ext)
        PickId = G2gd.GetPatternTreeItemId(self, self.root, export)
        file = open(filename,'w')
        file.write('#%s\n'%(export))
        print 'save powder pattern to file: ',filename
        wx.BeginBusyCursor()
        try:
            x,y,w,yc,yb,yd = self.PatternTree.GetItemPyData(PickId)[1]
            s = list(np.sqrt(1./np.array(w)))        
            XYW = zip(x,y,s)
            for X,Y,W in XYW:
                file.write("%15.6g %15.6g %15.6g\n" % (X,Y,W))
            file.close()
        finally:
            wx.EndBusyCursor()
        print 'powder pattern file written'
        
def PDFSave(self,exports):    
    for export in exports:
        PickId = G2gd.GetPatternTreeItemId(self, self.root, export)
        SQname = 'S(Q)'+export[4:]
        GRname = 'G(R)'+export[4:]
        sqfilename = ospath.join(self.dirname,export.replace(' ','_')[5:]+'.sq')
        grfilename = ospath.join(self.dirname,export.replace(' ','_')[5:]+'.gr')
        sqId = G2gd.GetPatternTreeItemId(self, PickId, SQname)
        grId = G2gd.GetPatternTreeItemId(self, PickId, GRname)
        sqdata = np.array(self.PatternTree.GetItemPyData(sqId)[1][:2]).T
        grdata = np.array(self.PatternTree.GetItemPyData(grId)[1][:2]).T
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
    
def PeakListSave(self,file,peaks):
    print 'save peak list to file: ',self.peaklistfile
    if not peaks:
        dlg = wx.MessageDialog(self, 'No peaks!', 'Nothing to save!', wx.OK)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        return
    for peak in peaks:
        file.write("%10.4f %12.2f %10.3f %10.3f \n" % \
            (peak[0],peak[2],peak[4],peak[6]))
    print 'peak list saved'
              
def IndexPeakListSave(self,peaks):
    file = open(self.peaklistfile,'wa')
    print 'save index peak list to file: ',self.peaklistfile
    wx.BeginBusyCursor()
    try:
        if not peaks:
            dlg = wx.MessageDialog(self, 'No peaks!', 'Nothing to save!', wx.OK)
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
    
def ReadEXPPhase(self,filename):
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
        dlg = wx.SingleChoiceDialog(self, 'Which phase to read?', 'Read phase data', PNames, wx.CHOICEDLG_STYLE)
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
            SHdata = EXPphase[key].split()
            textureData['Order'] = int(SHdata[0])
            textureData['Model'] = shModels[int(SHdata[2])]
            textureData['Sample omega'] = [False,float(SHdata[6])]
            textureData['Sample chi'] = [False,float(SHdata[7])]
            textureData['Sample phi'] = [False,float(SHdata[8])]
            shNcof = int(SHdata[1])
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
                Atom = [int(S[56:60]),S[50:54].strip().upper(),S[54:56],
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
        
    Phase = {
            'General':{
                'Name':PhaseName,
                'Type':Ptype,
                'SGData':SGData,
                'Cell':[False,]+abc+angles+[Volume,],
                'Pawley dmin':1.0},
            'Atoms':Atoms,
            'Drawing':{},
            'Histograms':{},
            'Pawley ref':[],
            'Models':{},
            'SH Texture':textureData
            }
            
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
            if E: 
                print ' ERROR in space group symbol ',SpGrp,' in file ',filename
                print ' N.B.: make sure spaces separate axial fields in symbol' 
                print G2spc.SGErrors(E)
                return None
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
    Phase['General'] = {'Name':PhaseName,'Type':'macromolecular','SGData':SGData,
        'Cell':[False,]+cell+[Volume,]}
    Phase['Atoms'] = Atoms
    Phase['Drawing'] = {}
    Phase['Histograms'] = {}
    
    return Phase
    
def ReadCIFPhase(filename):
    anisoNames = ['aniso_u_11','aniso_u_22','aniso_u_33','aniso_u_12','aniso_u_13','aniso_u_23']
    file = open(filename, 'Ur')
    Phase = {}
    Title = ospath.split(filename)[-1]
    print '\n Reading cif file: ',Title
    Compnd = ''
    Atoms = []
    A = np.zeros(shape=(3,3))
    S = file.readline()
    while S:
        if '_symmetry_space_group_name_H-M' in S:
            SpGrp = S.split("_symmetry_space_group_name_H-M")[1].strip().strip('"').strip("'")
            E,SGData = G2spc.SpcGroup(SpGrp)
            if E:
                print ' ERROR in space group symbol ',SpGrp,' in file ',filename
                print ' N.B.: make sure spaces separate axial fields in symbol' 
                print G2spc.SGErrors(E)
                return None
            S = file.readline()
        elif '_cell' in S:
            if '_cell_length_a' in S:
                a = S.split('_cell_length_a')[1].strip().strip('"').strip("'").split('(')[0]
            elif '_cell_length_b' in S:
                b = S.split('_cell_length_b')[1].strip().strip('"').strip("'").split('(')[0]
            elif '_cell_length_c' in S:
                c = S.split('_cell_length_c')[1].strip().strip('"').strip("'").split('(')[0]
            elif '_cell_angle_alpha' in S:
                alp = S.split('_cell_angle_alpha')[1].strip().strip('"').strip("'").split('(')[0]
            elif '_cell_angle_beta' in S:
                bet = S.split('_cell_angle_beta')[1].strip().strip('"').strip("'").split('(')[0]
            elif '_cell_angle_gamma' in S:
                gam = S.split('_cell_angle_gamma')[1].strip().strip('"').strip("'").split('(')[0]
            S = file.readline()
        elif 'loop_' in S:
            labels = {}
            i = 0
            while S:
                S = file.readline()
                if '_atom_site' in S.strip()[:10]:
                    labels[S.strip().split('_atom_site_')[1].lower()] = i
                    i += 1
                else:
                    break
            if labels:
                if 'aniso_label' not in labels:
                    while S:
                        atom = ['','','',0,0,0,1.0,'','','I',0.01,0,0,0,0,0,0]
                        S.strip()
                        if len(S.split()) != len(labels):
                            if 'loop_' in S:
                                break
                            S += file.readline().strip()
                        data = S.split()
                        if len(data) != len(labels):
                            break
                        for key in labels:
                            if key == 'type_symbol':
                                atom[1] = data[labels[key]]
                            elif key == 'label':
                                atom[0] = data[labels[key]]
                            elif key == 'fract_x':
                                atom[3] = float(data[labels[key]].split('(')[0])
                            elif key == 'fract_y':
                                atom[4] = float(data[labels[key]].split('(')[0])
                            elif key == 'fract_z':
                                atom[5] = float(data[labels[key]].split('(')[0])
                            elif key == 'occupancy':
                                atom[6] = float(data[labels[key]].split('(')[0])
                            elif key == 'thermal_displace_type':
                                if data[labels[key]].lower() == 'uiso':
                                    atom[9] = 'I'
                                    atom[10] = float(data[labels['u_iso_or_equiv']].split('(')[0])
                                else:
                                    atom[9] = 'A'
                                    atom[10] = 0.0
                                    
                        atom[7],atom[8] = G2spc.SytSym(atom[3:6],SGData)
                        atom.append(ran.randint(0,sys.maxint))
                        Atoms.append(atom)
                        S = file.readline()
                else:
                    while S:
                        S.strip()
                        data = S.split()
                        if len(data) != len(labels):
                            break
                        name = data[labels['aniso_label']]
                        for atom in Atoms:
                            if name == atom[0]:
                                for i,uname in enumerate(anisoNames):
                                    atom[i+11] = float(data[labels[uname]].split('(')[0])
                        S = file.readline()
                                                                        
        else:           
            S = file.readline()
    file.close()
    if Title:
        PhaseName = Title
    else:
        PhaseName = 'None'
    SGlines = G2spc.SGPrint(SGData)
    for line in SGlines: print line
    cell = [float(a),float(b),float(c),float(alp),float(bet),float(gam)]
    Volume = G2lat.calc_V(G2lat.cell2A(cell))
    Phase['General'] = {'Name':PhaseName,'Type':'nuclear','SGData':SGData,
        'Cell':[False,]+cell+[Volume,]}
    Phase['Atoms'] = Atoms
    Phase['Drawing'] = {}
    Phase['Histograms'] = {}
    
    return Phase

def ValEsd(value,esd=0,nTZ=False):                  #NOT complete - don't use
    # returns value(esd) string; nTZ=True for no trailing zeros
    # use esd < 0 for level of precision shown e.g. esd=-0.01 gives 2 places beyond decimal
    #get the 2 significant digits in the esd 
    edig = lambda esd: int(round(10**(math.log10(esd) % 1+1)))
    #get the number of digits to represent them 
    epl = lambda esd: 2+int(1.545-math.log10(10*edig(esd)))
    
    mdec = lambda esd: -int(math.log10(abs(esd)))
    ndec = lambda esd: int(1.545-math.log10(abs(esd)))
    if esd > 0:
        fmt = '"%.'+str(ndec(esd))+'f(%d)"'
        print fmt,ndec(esd),esd*10**(mdec(esd)+1)
        return fmt%(value,int(esd*10**(mdec(esd)+1)))
    elif esd < 0:
         return str(round(value,mdec(esd)))
    else:
        text = "%F"%(value)
        if nTZ:
            return text.rstrip('0')
        else:
            return text

def Fesd(value,esd=0,nTZ=False):
#pythonized version of fortran routine in GSAS cifsubs directory - doesn't work correctly
    nint = lambda x: int(round(x))
    iExp = 0
    if value == 0. and esd == 0.:
        iDec = 1
        iFld = 5
    elif value == 0.:
        iDec = max(0.,1.545-math.log10(abs(esd)))
        iFld = 4+iDec
    elif esd == 0.:
        iDec = 5
        iFld = max(1.,math.log10(abs(value)))+3+iDec
    else:
        iFld = math.log10(max(abs(esd),abs(value)))
        if iFld < -4:
            iExp = 1-iFld
            iFld -= iExp
        elif iFld > 8:
            iExp = -iFld
            iFld += iExp
        if iExp:
            value *= 10.**iExp
            esd *= 10.**iExp
        iDec = min(7,int(max(0.,1.545-math.log10(max(0.000001*abs(value),abs(esd))))))
        iFld = max(1,iFld)+3+iDec
    if esd <= 0.:
        iSigw = 0
    else:
        iSig = nint(esd*(10.**iDec))
        iSigw = 1
        if iSig > 0:
            iSigw = int(1.+math.log10(1.*iSig))
    if iSigw > 2:
        xmult = 10.**(iSigw-2)
        value = xmult*nint(value/xmult)
        iSig = xmult*nint(iSig/xmult)            
        iSigw = int(1.+math.log10(1.*iSig))
    if iSigw == 0:
        fmt = '%.'+str(iDec)+'f'
        string = fmt%(value) 
    elif iDec > 0:
        fmt = '%.'+str(iDec)+'f(%d)'
        string = fmt%(value,iSig)
    else:
        fmt = '%'+str(iFld)+'d(%d)'
        string = fmt%(nint(value),iSig)
    if iExp:
        iDig = 1+math.log10(abs(1.*iExp))
        if iExp > 0:
            iDig += 1
            string += str(-iExp)
    return string
    

