"""GSASIIIO: functions for IO of data
   Copyright: 2008, Robert B. Von Dreele (Argonne National Laboratory)
"""

import wx
import math
import numpy as np
import cPickle
import sys
import GSASIIpath
import GSASIIgrid as G2gd
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
            else:
                Comments.append(S[:-1])
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
    return FoundData,Iparm,Comments

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
                    HKLref.append([HKL,Fosq,sigFosq,Fcsq,0,0,0])                 #room for Fc, Fcp, Fcpp & phase
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
    
def GetImageData(imagefile,imageOnly=False):
    ext = ospath.splitext(imagefile)[1]
    Comments = []
    if ext == '.tif':
        Comments,Data,Size,Image = GetTifData(imagefile)
    elif ext == '.img':
        Comments,Data,Size,Image = GetImgData(imagefile)
        Image[0][0] = 0
    elif ext == '.mar3450' or ext == '.mar2300':
        Comments,Data,Size,Image = GetMAR345Data(imagefile)
    elif ext in ['.sum','.avg']:
        Comments,Data,Size,Image = GetGEsumData(imagefile)
    elif ext == '.G2img':
        return GetG2Image(imagefile)
    if imageOnly:
        return Image
    else:
        return Comments,Data,Size,Image
        
def PutG2Image(filename,image):
    File = open(filename,'wb')
    cPickle.dump(image,File,1)
    File.close()
    return
    
def GetG2Image(filename):
    File = open(filename,'rb')
    image = cPickle.load(File)
    File.close()
    return image
    
def GetGEsumData(filename,imageOnly=False):
    import array as ar
    if not imageOnly:
        print 'Read GE sum file: ',filename    
    File = open(filename,'rb')
    size = 2048
    if '.sum' in filename:
        head = ['GE detector sum data from APS 1-ID',]
    if '.avg' in filename:
        head = ['GE detector avg data from APS 1-ID',]
    image = np.zeros(shape=(size,size),dtype=np.int32)
    row = 0
    pos = 0
    while row < size:
        File.seek(pos)
        if '.sum' in filename:
            line = ar.array('f',File.read(4*size))
            pos += 4*size
        elif '.avg' in filename:
            line = ar.array('H',File.read(2*size))
            pos += 2*size
        image[row] = np.asarray(line)
        row += 1
    data = {'pixelSize':(200,200),'wavelength':0.15,'distance':250.0,'center':[204.8,204.8]}  
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,size,image
        
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
    data = {'pixelSize':pixel,'wavelength':wave,'distance':distance,'center':center}
    image = []
    row = 0
    pos = 512
    image = np.zeros(shape=(size,size),dtype=np.int32)    
    while row < size:
        File.seek(pos)
        line = ar.array('H',File.read(2*size))
        image[row] = np.asarray(line)
        row += 1
        pos += 2*size
    File.close()
    if imageOnly:
        return image
    else:
        return lines[1:-2],data,size,image
       
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
    pos = 4096
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
        return head,data,size,image.T
    
def GetTifData(filename,imageOnly=False):
    # only works for APS Perkin-Elmer detector data files in "TIFF" format that are readable by Fit2D
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
        head = 'no metadata file found'
    tag = File.read(3)
    if tag != 'II*':
        lines = ['not a detector tiff file',]
        return lines,0,0
    size,Ityp = st.unpack('<ii',File.read(8))
    if Ityp == 0:
        tifType = 'Pilatus'
        pixy = (172,172)
        pos = 4096
        if not imageOnly:
            print 'Read Pilatus tiff file: ',filename
    elif Ityp == 1:
        tifType = 'PE'
        pixy = (200,200)
        pos = 8
        if not imageOnly:
            print 'Read APS PE-detector tiff file: ',filename
    elif Ityp == 3328:
        tifType = 'MAR'
        pixy = (79,79)
        pos = 4096
        if not imageOnly:
            print 'Read MAR CCD tiff file: ',filename
    else:
        lines = 'unknown tif type'
        return lines,0,0
    image = np.zeros(shape=(size,size),dtype=np.int32)
    row = 0
    while row < size:
        File.seek(pos)
        if 'PE' in tifType: 
            if dataType == 5:
                line = ar.array('f',File.read(4*size))
            else:
                line = ar.array('l',File.read(4*size))
            pos += 4*size
        elif 'MAR' in tifType:
            line = ar.array('H',File.read(2*size))
            pos += 2*size
        image[row] = np.asarray(line)
        row += 1
    data = {'pixelSize':pixy,'wavelength':0.10,'distance':100.0,'center':[204.8,204.8]}
    File.close()    
    if imageOnly:
        return image
    else:
        return head,data,size,image
    
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
            print 'load: ',datum[0]
            
            #temporary fixes to old project files
            #fix to convert old style list arrays to numpy arrays
            if 'PWDR' in datum[0] and 'list' in str(type(datum[1][1][0])):      
                X = datum[1][1]
                X = [np.array(X[0]),np.array(X[1]),np.array(X[2]),np.array(X[3]),np.array(X[4]),np.array(X[5])]
                datum[1] = [datum[1][0],X]
                print 'powder data converted to numpy arrays'
            #temporary fix to insert 'PWDR' in front of powder names
            if 'PKS' not in datum[0] and 'IMG' not in datum[0] and 'SNGL' not in datum[0]:
                if datum[0] not in ['Notebook','Controls','Phases'] and 'PWDR' not in datum[0]:
                    datum[0] = 'PWDR '+datum[0]
                    print 'add PWDR to powder names'
            #end of temporary fixes
            
            Id = self.PatternTree.AppendItem(parent=self.root,text=datum[0])
            self.PatternTree.SetItemPyData(Id,datum[1])
            for datus in data[1:]:
                print '    load: ',datus[0]
                
                #temporary fix to add azimuthal angle to instrument parameters
                if 'PWDR' in datum[0] and 'Instrument Parameters' in datus[0]:
                    if len(datus[1][0]) == 10 or len(datus[1][0]) == 12:
                        datus[1][0] += (0.0,)                   #add missing azimuthal angle
                        datus[1][1].append(0.0)
                        datus[1][2].append(0.0)
                        datus[1][3].append('Azimuth')
                        print 'add azimuth to instrument parameters'
                #end of temporary fix        
                
                sub = self.PatternTree.AppendItem(Id,datus[0])
                self.PatternTree.SetItemPyData(sub,datus[1])
                
            #temporary fix to add Comments to powder data sets
            if 'PWDR' in datum[0] and not G2gd.GetPatternTreeItemId(self,Id, 'Comments'):
                print 'no comments - add comments'
                sub = self.PatternTree.AppendItem(Id,'Comments')
                self.PatternTree.SetItemPyData(sub,['no comments'])
            #end of temporary fix
                                
            if 'IMG' in datum[0]:                   #retreive image default flag & data if set
                Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Image Controls'))
                if Data['setDefault']:
                    self.imageDefault = Data
                #temporary fix to add masks
                if not G2gd.GetPatternTreeItemId(self,Id, 'Masks'):
                    Imin,Imax = Data['range'][0]
                    Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                    self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                #end of temporary fix
                
        file.close()
        
        #temporary fix to add Notebook & Controls to project
        if not G2gd.GetPatternTreeItemId(self,self.root,'Notebook'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Notebook')
            self.PatternTree.SetItemPyData(sub,[''])
            sub = self.PatternTree.AppendItem(parent=self.root,text='Controls')
            self.PatternTree.SetItemPyData(sub,[0])
            print 'add Notebook and Controls to project'
            
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
                print 'save: ',name
                data.append([name,self.PatternTree.GetItemPyData(item)])
                item2, cookie2 = self.PatternTree.GetFirstChild(item)
                while item2:
                    name = self.PatternTree.GetItemText(item2)
                    print '    save: ',name
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
    X = self.Integrate[2].flatten()[:-1]
    Xminmax = [X[0],X[-1]]
    N = len(X)
    Id = self.PatternTree.GetItemParent(PickId)
    name = self.PatternTree.GetItemText(Id)
    name = name.replace('IMG ','PWDR ')
    Comments = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,Id, 'Comments'))
    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
    codes = [0 for i in range(11)]
    parms = ['PXC',data['wavelength'],0.0,0.0,1.0,-1.0,0.3,0.0,1.0,0.0,0.0]
    Azms = [(azms[i+1]+azms[i])/2. for i in range(len(azms)-1)]
    for i,azm in enumerate(Azms):
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        Id = 0
        while item:
            Name = self.PatternTree.GetItemText(item)
            if name == Name:
                Id = item
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        parms[10] = azm
        Y = self.Integrate[0][i].flatten()
        W = np.sqrt(Y)
        if Id:
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id, 'Comments'),Comments)                    
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Limits'),[tuple(Xminmax),Xminmax])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Instrument Parameters'),[tuple(parms),parms,codes,names])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Peak List'),[])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Index Peak List'),[])
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,Id,'Unit Cells List'),[])             
        else:
            Id = self.PatternTree.AppendItem(parent=self.root,text=name+" Azm= %.2f"%(azm))
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(parms),parms,codes,names])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
        self.PatternTree.SetItemPyData(Id,[[''],[np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]])
    self.PatternTree.SelectItem(Id)
    self.PatternTree.Expand(Id)
    self.PatternId = Id
            
def powderFxyeSave(self,powderfile):
    file = open(powderfile,'w')
    prm = open(powderfile.strip('fxye')+'prm','w')      #old style GSAS parm file
    print 'save powder pattern to file: ',powderfile
    wx.BeginBusyCursor()
    Inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self, \
                    self.PickId, 'Instrument Parameters'))[1]
    if len(Inst) == 11:             #single wavelength
        lam1 = Inst[1]
        lam2 = 0.0
        GU,GV,GW = Inst[4:7]
        LX,LY = Inst[7:9]
        SL = HL = Inst[9]/2.0   
    else:                           #Ka1 & Ka2
        lam1 = Inst[1]
        lam2 = Inst[2]
        GU,GV,GW = Inst[6:9]
        LX,LY = Inst[9:11]
        SL = HL = Inst[11]/2.0   
    prm.write( '            123456789012345678901234567890123456789012345678901234567890        '+'\n')
    prm.write( 'INS   BANK      1                                                               '+'\n')
    prm.write( 'INS   HTYPE   PXCR                                                              '+'\n')
    prm.write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   '+'\n')%(lam1,lam2))
    prm.write( 'INS  1 IRAD     0                                                               '+'\n')
    prm.write( 'INS  1I HEAD                                                                    '+'\n')
    prm.write( 'INS  1I ITYP    0    0.0000  180.0000         1                                 '+'\n')
    prm.write( 'INS  1PRCF1     3    8   0.00100                                                '+'\n')
    prm.write(('INS  1PRCF11     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(GU,GV,GW,0.0))
    prm.write(('INS  1PRCF12     %15.6g%15.6g%15.6g%15.6g   '+'\n')%(LX,LY,SL,HL))
    prm.close()
    try:
        x,y,w,yc,yb,yd = self.PatternTree.GetItemPyData(self.PickId)[1]
        file.write(powderfile+'\n')
        file.write('BANK 1 %d %d CONS %.2f %.2f 0 0 FXYE\n'%(len(x),len(x),\
            100.*x[0],100.*(x[1]-x[0])))        
        XYW = zip(x,y,w)
        for X,Y,W in XYW:
            file.write("%15.6g %15.6g %15.6g\n" % (100.*X,Y,W))
        file.close()
    finally:
        wx.EndBusyCursor()
    print 'powder pattern file written'
        
def powderXyeSave(self,powderfile):
    file = open(powderfile,'w')
    print 'save powder pattern to file: ',powderfile
    wx.BeginBusyCursor()
    try:
        x,y,w,yc,yb,yd = self.PatternTree.GetItemPyData(self.PickId)[1]
        XYW = zip(x,y,w)
        for X,Y,W in XYW:
            file.write("%15.6g %15.6g %15.6g\n" % (X,Y,W))
        file.close()
    finally:
        wx.EndBusyCursor()
    print 'powder pattern file written'
    
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
    
def ReadEXPPhase(filename):
    import GSASIIspc as G2spc
    import GSASIIlattice as G2lat
    import GSASIIElem as G2el
    file = open(filename, 'Ur')
    Phase = {}
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
                Atoms.append(Atom)
    Volume = G2lat.calc_V(G2lat.cell2A(abc+angles))
    Phase['General'] = [PhaseName,Ptype,SGData,[False,]+abc+angles+[Volume,],[False,1.0],
        0,0,0,0,0]
    Phase['Atoms'] = Atoms
    return Phase
       
def ReadPDBPhase(filename):
    import GSASIIspc as G2spc
    import GSASIIlattice as G2lat
    import GSASIIElem as G2el
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
            S = file.readline()
        elif 'SCALE' in S[:5]:
            V = (S[10:41].split())
            A[int(S[5])-1] = [float(V[0]),float(V[1]),float(V[2])]
            S = file.readline()
        elif 'ATOM' in S[:4] or 'HETATM' in S[:6]:
            XYZ = [float(S[31:39]),float(S[39:47]),float(S[47:55])]
            XYZ = np.sum(AA*XYZ,axis=1)
            SytSym,Mult = G2spc.SytSym(XYZ,SGData)
            Uiso = float(S[61:67])/EightPiSq
            Type = S[12:14].upper()
            if Type[0] in '123456789':
                Type = Type[1:]
            Atom = [int(S[22:27]),S[17:20].upper(),S[20:22],
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
    Phase['General'] = [PhaseName,'macromolecular',SGData,[False,]+cell+[Volume,],[False,1.0],
        0,0,0,0,0]
    Phase['Atoms'] = Atoms
    
    return Phase
    
def ReadCIFAtoms(self,data):
    print data
