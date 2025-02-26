# -*- coding: utf-8 -*-
'''
'''
import os
import struct as st
import numpy as np
from .. import GSASIIobj as G2obj
class raw_ReaderClass(G2obj.ImportPowderData):
    'Routines to import powder data from a binary Bruker .RAW file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.RAW',),
            strictExtension=False,
            formatName = 'Bruker RAW',
            longFormatName = 'Bruker .RAW powder data file'
            )
        self.scriptable = True

    # Validate the contents -- make sure we only have valid lines
    def Read(self,fp,nbytes):
        data = fp.read(nbytes)
        if 'bytes' in str(type(data)):
            data = data.decode('latin-1')
        return data
    
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid Bruker RAW file'
        fp = open(filename,'rb')
        head = self.Read(fp,7)
        if 'bytes' in str(type(head)):
            head = head.decode('latin-1')
        if head[:4] == 'RAW ':
            self.formatName = 'Bruker RAW ver. 1'
        elif head[:4] == 'RAW2':
            self.formatName = 'Bruker RAW ver. 2'
        elif head == 'RAW1.01':
            self.formatName = 'Bruker RAW ver. 3'
        elif head == 'RAW4.00':
            self.formatName = 'Bruker RAW ver. 4'
            pwdrscan = fp.read()
            nBanks = pwdrscan.count(b'2Theta')
            if not len(self.selections):
                self.selections = list(range(nBanks))
                self.numbanks = nBanks
            for i in range(nBanks):
                self.dnames.append(os.path.basename(filename)+' scan '+(str(i)))
        else:
            self.errors = 'Unexpected information in header: '
            if all([ord(c) < 128 and ord(c) != 0 for c in str(head)]): # show only if ASCII
                self.errors += '  '+str(head)
            else: 
                self.errors += '  (binary)'
            fp.close()
            return False
        fp.close()
        return True
            
    def Reader(self,filename, ParentFrame=None, **kwarg):
        'Read a Bruker RAW file'
        self.comments = []
        self.powderentry[0] = filename
        fp = open(filename,'rb')
        if 'ver. 1' in self.formatName:
            raise Exception('Read of Bruker "RAW " (pre-version #) file not supported')    #for now
        elif 'ver. 2' in self.formatName:
            fp.seek(4)
            nBlock = int(st.unpack('<i',fp.read(4))[0])
            fp.seek(168)
            self.comments.append('Date/Time='+self.Read(fp,20))
            self.comments.append('Anode='+self.Read(fp,2))
            self.comments.append('Ka1=%.5f'%(st.unpack('<f',fp.read(4))[0]))
            self.comments.append('Ka2=%.5f'%(st.unpack('<f',fp.read(4))[0]))
            self.comments.append('Ka2/Ka1=%.5f'%(st.unpack('<f',fp.read(4))[0]))
            fp.seek(206)
            self.comments.append('Kb=%.5f'%(st.unpack('<f',fp.read(4))[0]))
            pos = 256
            fp.seek(pos)
            blockNum = kwarg.get('blocknum',0)
            self.idstring = os.path.basename(filename) + ' Scan '+str(blockNum)
            if blockNum <= nBlock:
                for iBlock in range(blockNum):
                    headLen = int(st.unpack('<H',fp.read(2))[0])
                    nSteps = int(st.unpack('<H',fp.read(2))[0])
                    if iBlock+1 == blockNum:
                        fp.seek(pos+12)
                        step = st.unpack('<f',fp.read(4))[0]
                        start2Th = st.unpack('<f',fp.read(4))[0]
                        pos += headLen      #position at start of data block
                        fp.seek(pos)                                    
                        x = np.array([start2Th+i*step for i in range(nSteps)])
                        y = np.array([max(1.,st.unpack('<f',fp.read(4))[0]) for i in range(nSteps)])
                        y = np.where(y<0.,1.,y)
                        w = 1./y
                        self.powderdata = [x,y,w,np.zeros(nSteps),np.zeros(nSteps),np.zeros(nSteps)]
                        break
                    pos += headLen+4*nSteps
                    fp.seek(pos)
                if blockNum == nBlock:
                    self.repeat = False                                   
                else:
                    self.repeat = True
            fp.close()
        elif 'ver. 3' in self.formatName:
            fp.seek(12)
            nBlock = int(st.unpack('<i',fp.read(4))[0])
            self.comments.append('Date='+self.Read(fp,10))
            self.comments.append('Time='+self.Read(fp,10))
            fp.seek(326)
            self.comments.append('Sample='+self.Read(fp,60))
            fp.seek(564)
            radius = st.unpack('<f',fp.read(4))[0]
            self.comments.append('Gonio. radius=%.2f'%(radius))
            self.Sample['Gonio. radius'] = radius
            fp.seek(608)
            self.comments.append('Anode='+self.Read(fp,4))
            fp.seek(616)
            self.comments.append('Ka mean=%.5f'%(st.unpack('<d',fp.read(8))[0]))
            self.comments.append('Ka1=%.5f'%(st.unpack('<d',fp.read(8))[0]))
            self.comments.append('Ka2=%.5f'%(st.unpack('<d',fp.read(8))[0]))
            self.comments.append('Kb=%.5f'%(st.unpack('<d',fp.read(8))[0]))
            self.comments.append('Ka2/Ka1=%.5f'%(st.unpack('<d',fp.read(8))[0]))
            pos = 712
            fp.seek(pos)      #position at 1st block header
            blockNum = kwarg.get('blocknum',0)
            self.idstring = os.path.basename(filename) + ' Scan '+str(blockNum)
            if blockNum <= nBlock:
                for iBlock in range(blockNum):
                    headLen = int(st.unpack('<i',fp.read(4))[0])
                    nSteps = int(st.unpack('<i',fp.read(4))[0])
                    if not nSteps: break
                    if nBlock > 1:
                        fp.seek(pos+256)
                        headLen += st.unpack('<i',fp.read(4))[0]
                    else:
                        headLen += 40
                    if iBlock+1 == blockNum:
                        fp.seek(pos+8)
                        st.unpack('<d',fp.read(8))[0]
                        start2Th = st.unpack('<d',fp.read(8))[0]
                        fp.seek(pos+212)
                        temp = st.unpack('<f',fp.read(4))[0]
                        if temp > 0.:
                            self.Sample['Temperature'] = temp                                                        
                        fp.seek(pos+176)
                        step = st.unpack('<d',fp.read(8))[0]
                        pos += headLen      #position at start of data block
                        fp.seek(pos)
                        x = np.array([start2Th+i*step for i in range(nSteps)])
                        try:
                            y = np.array([max(1.,st.unpack('<f',fp.read(4))[0]) for i in range(nSteps)])
                        except: #this is absurd
                            fp.seek(pos-40)
                            y = np.array([max(1.,st.unpack('<f',fp.read(4))[0]) for i in range(nSteps)])
                        w = 1./y
                        self.powderdata = [x,y,w,np.zeros(nSteps),np.zeros(nSteps),np.zeros(nSteps)]
                        break
                    pos += headLen+4*nSteps
                    fp.seek(pos)
                if blockNum == nBlock:
                    self.repeat = False
                else:
                    self.repeat = True
            fp.close()
            
        elif 'ver. 4' in self.formatName: 
            driveNo = 0
            fp.seek(12)   #ok
            self.comments.append('Date='+self.Read(fp,12).strip('\x00'))
            self.comments.append('Time='+self.Read(fp,10).strip('\x00'))
            fp.seek(61)     #start of header segments
            nBank = 0
            blockNum = kwarg.get('blocknum',0)
            while nBank < self.numbanks:
                while True:     #read block header
                    segtype = st.unpack('<I',fp.read(4))[0]
                    if not segtype or segtype == 160:    
                        break           # done with header
                    seglen = max(st.unpack('<I',fp.read(4))[0],8)
                    if segtype == 10:
                        fp.read(4)    #skip these
                        self.comments.append('%s=%s'%(self.Read(fp,24).strip('\x00'),self.Read(fp,seglen-36).strip('\x00')))
                    elif segtype == 30: #x-ray source info
                        fp.read(64)
                        self.comments.append('Ka mean=%.5f'%(st.unpack('<d',fp.read(8))[0]))
                        self.comments.append('Ka1=%.5f'%(st.unpack('<d',fp.read(8))[0]))
                        self.comments.append('Ka2=%.5f'%(st.unpack('<d',fp.read(8))[0]))
                        self.comments.append('Kb=%.5f'%(st.unpack('<d',fp.read(8))[0]))
                        self.comments.append('Ka2/Ka1=%.5f'%(st.unpack('<d',fp.read(8))[0]))
                        fp.read(4)
                        self.comments.append('Anode='+self.Read(fp,4).strip('\x00'))
                        fp.read(seglen-120)
                    elif segtype == 60:
                        alignFlag = st.unpack('<I',fp.read(4))[0]
                        driveName = self.Read(fp,24).strip('\x00')
                        fp.read(32)
                        Delt = st.unpack('<d',fp.read(8))[0]
                        fp.read(seglen-76)
                        self.comments.append('Drive %s: align flag %d'%(driveName,alignFlag))
                        self.comments.append('Drive %s: delta %f'%(driveName,Delt))
                        driveNo += 1
                    else:
                        fp.read(seglen-8)
                if (segtype == 0 or segtype == 160):    #read data block
                    self.idstring = self.dnames[nBank]
                    meta = {}
                    fp.read(28)
                    meta['ScanType'] = self.Read(fp,24).strip('\x00')
                    if meta['ScanType'] not in ['Locked Coupled','Unlocked Coupled','Detector Scan']:
                        return False
                    fp.read(16)
                    startAngle = st.unpack('<d',fp.read(8))[0]
                    meta['startAngle'] = '%.4f'%startAngle
                    stepSize = st.unpack('<d',fp.read(8))[0]
                    meta['stepSize'] = '%.4f'%stepSize
                    Nsteps = st.unpack('<I',fp.read(4))[0]
                    meta['Nsteps'] = '%d'%Nsteps
                    meta['stepTime(ms)'] = st.unpack('<f',fp.read(4))[0]
                    fp.read(4)
                    meta['generatorVoltage(kV)'] = st.unpack('<f',fp.read(4))[0]
                    meta['generatorCurrent(mA)'] = st.unpack('<f',fp.read(4))[0]
                    fp.read(4)
                    meta['usedWave'] = st.unpack('<d',fp.read(8))[0]
                    fp.read(16)
                    datumSize = st.unpack('<I',fp.read(4))[0]
                    hdrSize = st.unpack('<I',fp.read(4))[0]
                    fp.read(16)
                    if meta['ScanType'] in ['Locked Coupled','Unlocked Coupled','Detector Scan']:
                        while hdrSize > 0:
                            segtype = st.unpack('<I',fp.read(4))[0]
                            seglen = max(st.unpack('<I',fp.read(4))[0],8)
                            if segtype == 50:
                                fp.read(4)
                                segName = self.Read(fp,24).strip('\x00')
                                if segName in ['Theta','2Theta','Chi','Phi','BeamTranslation','Z-Drive','Divergence Slit']:
                                    fp.read(20)
                                    meta['start %s'%segName] = '%.4f'%(st.unpack('<d',fp.read(8))[0])
                                    fp.read(seglen-64)
                                else:
                                    fp.read(seglen-36)
                            else:
                                fp.read(seglen-8)
                            hdrSize -= seglen
                        #end of reading scan header
                        pos = fp.tell()
                        fp.seek(pos-16)
                        meta['Temperature'] = st.unpack('<f',fp.read(4))[0]
                        if meta['Temperature'] > 7.:  #one raw4 file had int4='9999' in this place & <7K unlikely for lab data
                            self.Sample['Temperature'] = meta['Temperature']
                        self.Sample['Omega'] = meta['start Theta']
                        fp.read(12)
                        x = np.array([startAngle+i*stepSize for i in range(Nsteps)])
                        y = np.array([max(1.,st.unpack('<f',fp.read(4))[0]) for i in range(Nsteps)])
                        w = 1./y
                        if nBank == blockNum-1:
                            self.powderdata = [x,y,w,np.zeros(Nsteps),np.zeros(Nsteps),np.zeros(Nsteps)]
                            for item in meta:
                                self.comments.append('%s = %s'%(item,str(meta[item])))
                            fp.close()
                            self.repeat = True
                            if nBank == self.numbanks-1:
                                self.repeat = False
                            return True
                    else:
                        meta['Unknown range/scan type'] = True
                        fp.read(hdrSize)
                        fp.read(datumSize*Nsteps)
                nBank += 1
        else:
            return False
        self.repeat = False    
        return True

