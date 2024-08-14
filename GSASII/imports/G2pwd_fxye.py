# -*- coding: utf-8 -*-
'''
'''
from __future__ import division, print_function
import os.path as ospath
import platform
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath

class GSAS_ReaderClass(G2obj.ImportPowderData):
    'Routines to import powder data from a GSAS files'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.fxye','.raw','.gsas','.gda','.gsa','.gss','.RAW','.GSAS','.GDA','.GSA','.dat'),
            strictExtension=False,
            formatName = 'GSAS powder data',
            longFormatName = 'GSAS powder data files (.fxye, .raw, .gsas...)'
            )
        self.clockWd = {}
        self.TimeMap = {}
        self.dnames = []
        self.scriptable = True
        self.Iparm = {} #only filled for EDS data

    # Validate the contents -- look for a bank line
    def ContentsValidator(self, filename):
        'Validate by checking to see if the file has BANK lines & count them'
        #print 'ContentsValidator: '+self.formatName
        nBanks= 0
        if '2' in platform.python_version_tuple()[0]:
            fp = open(filename,'Ur')
        else:
            fp = open(filename,'r',encoding='latin-1')
        fname = ospath.basename(fp.name)
        for i,line in enumerate(fp):
            self.GSAS = True
            if i==0: # first line is always a comment
                continue
            if i==1 and line[:4].lower() == 'inst':
                # 2nd line is optional instrument parameter file
                continue
            if line[0] == '#': continue
            if line[:4] == 'BANK':
                self.dnames.append(fname+' '+' '.join(line.split()[:2]))
                nBanks += 1
                continue
            elif line[:7] == 'Monitor': continue
            elif line [:8] == 'TIME_MAP':          #LANSCE TOF data
                continue
#            else:
#                if not all(c in ' 0123456789.eE+-' for c in line[:-1]):
#                    self.errors = 'Unexpected information in line: '+str(i+1)
#                    self.errors += '  '+str(line)
#                    return False
        if nBanks:
            if not len(self.selections):
                self.selections = range(nBanks)
            fp.close()
            return True
        self.errors = 'No BANK records found'
        fp.close()
        return False # no bank records

    def Reader(self,filename, ParentFrame=None, **kwarg):
        '''Read a GSAS (old formats) file of type FXY, FXYE, ESD or STD types.
        If multiple datasets are requested, use self.repeat and buffer caching.
        
        EDS data is only in the STD format (10 values per line separated by spaces); 
        the 1st line contains at col 60 the word "Two-Theta " followed by the appropriate value. 
        The BANK record contains the 3 values (4th not used) after 'EDS' for converting MCA 
        channel number (c) to keV via E = A + Bc + Cc^2; these coefficients are 
        generally predetermined by calibration of the MCA. They & 2-theta are transferred to 
        the Instrument parameters data.
        '''
        def GetFXYEdata(File,Pos,Bank):
            File.seek(Pos)
            x = []
            y = []
            w = []
            S = File.readline()
            while S and S[:4] != 'BANK' and S[0] != '#':
                vals = S.split()
                x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
                f = float(vals[1])
                s = float(vals[2])
                if f <= 0.0 or s <= 0.0:
                    y.append(0.0)
                    w.append(0.0)
                else:
                    y.append(float(vals[1]))
                    w.append(1.0/float(vals[2])**2)
                S = File.readline()
            N = len(x)
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]    
            
        def GetFXYdata(File,Pos,Bank):
            File.seek(Pos)
            x = []
            y = []
            w = []
            S = File.readline()
            while S and S[:4] != 'BANK' and S[0] != '#':
                vals = S.split()
                x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
                f = float(vals[1])
                if f > 0.0:
                    y.append(f)
                    w.append(1.0/f)
                else:              
                    y.append(0.0)
                    w.append(0.0)
                S = File.readline()
            N = len(x)
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            
        def GetESDdata(File,Pos,Bank):
            File.seek(Pos)
            cons = Bank.split()
            if 'TIME_MAP' == cons[4]:
                start = 0
                step = 1
            else:
                start = float(cons[5])/100.0               #CW: from centidegrees to degrees
                step = float(cons[6])/100.0
            x = []
            y = []
            w = []
            S = File.readline()
            j = 0
            while S and S[:4] != 'BANK' and S[0] != '#':
                if 'TIME_MAP' in S:
                    break
                for i in range(0,80,16):
                    if S[i:i+8] == 8*' ':
                        break
                    xi = start+step*j
                    yi = sfloat(S[i:i+8])
                    ei = sfloat(S[i+8:i+16])
                    x.append(xi)
                    if yi > 0.0:
                        y.append(yi)
                        w.append(1.0/ei**2)
                    else:              
                        y.append(0.0)
                        w.append(0.0)
                    j += 1
                S = File.readline()
            N = len(x)
            if self.clockWd:
                if len(self.clockWd) == 1:
                    x = Tmap2TOF(self.TimeMap['1'],self.clockWd['1'])
                else:
                    x = Tmap2TOF(self.TimeMap[cons[5]],self.clockWd[cons[5]])
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
        
        def GetSTDdata(File,Pos,Bank):
            File.seek(Pos)
            cons = Bank.split()
            Nch = int(cons[2])
            if self.clockWd:
                start = 0
                step = 1
            else:
                start = float(cons[5])/100.0               #CW: from centidegrees to degrees
                step = float(cons[6])/100.0                 #NB TOF 0.1*ms!
            x = []
            y = []
            w = []
            S = File.readline()
            j = 0
            Ecoef = np.zeros(4)
            if 'EDS' in cons[4]:
               Ecoef = np.fromstring(' '.join(cons[5:9]),sep=' ')
               self.Inst = {'XE':[Ecoef[0],Ecoef[0],False],'YE':[Ecoef[1],Ecoef[1],False],
                'ZE':[Ecoef[2],Ecoef[2],False],'WE':[Ecoef[3],Ecoef[3],False],'2-theta':[5.,5.,False]}
            while S and S[:4] != 'BANK' and S[0] != '#':
                if 'TIME_MAP' in S or '\x1a' in S:
                    break
                for i in range(0,80,8):
                    if S[i:i+10] == 10*' ':
                        break
                    xi = start+step*j
                    if 'EDS' in Bank:   #Energy dispersive - convert CN to EkeV
                        xi = Ecoef[0]+(Ecoef[1]+Ecoef[2]*j)*j+Ecoef[3]*j**3
                    ni = max(sint(S[i:i+2]),1)
                    yi = max(sfloat(S[i+2:i+8]),0.0)
                    if yi:
                        vi = yi/ni
                    else:
                        yi = 0.0
                        vi = 0.0
                    j += 1
                    if j < Nch:
                        x.append(xi)
                        if vi <= 0.:
                            y.append(0.)
                            w.append(0.)
                        else:
                            y.append(yi)
                            w.append(1.0/vi)
                S = File.readline()
            N = len(x)
            if self.clockWd:
                if len(self.clockWd) == 1:
                    x = Tmap2TOF(self.TimeMap['1'],self.clockWd['1'])
                else:
                    x = Tmap2TOF(self.TimeMap[cons[5]],self.clockWd[cons[5]])
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
           
        def GetALTdata(File,Pos,Bank):
            File.seek(Pos)
            cons = Bank.split()
            x = []
            y = []
            w = []
            S = File.readline()
            j = 0
            while S and S[:4] != 'BANK' and S[0] != '#':
                for i in range(0,80,20):
                    if S[i:i+8] == 8*' ':
                        break
                    xi = sfloat(S[i:i+8])/3200.
                    if not xi:
                        break
                    yi = sfloat(S[i+8:i+15])/1000.
                    ei = sfloat(S[i+15:i+21])/1000.
                    x.append(xi)
                    if yi > 0.0:
                        y.append(yi)
                        w.append(1.0/ei**2)
                    else:              
                        y.append(0.0)
                        w.append(0.0)
                    j += 1
                S = File.readline()
            N = len(x)
            if self.clockWd:
                x = Tmap2TOF(self.TimeMap[cons[5]],self.clockWd[cons[5]])
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            
        def GetTimeMap(File,Pos,TimeMap):
            File.seek(Pos)
            cons = TimeMap[8:].split()
            mapNo = cons[0]
            if 'TIME_MAP10' in TimeMap: mapNo = '1'   #HIPD cluge!
            Nch = int(cons[1])
            Nrec = int(cons[2])
            clockWd = float(cons[4])/1000.          #in mus
            TMap = np.zeros(Nch+2,dtype=int)
            ind = 0
            for i in range(Nrec):
                S = File.readline().rstrip('\n')
                vals = S.split()
                for val in vals:
                    TMap[ind] = int(val)
                    ind += 1
            TMap = np.reshape(TMap,(-1,3))
            TMax = TMap[-1][0]
            Nch = TMap[-2][0]+(TMax-TMap[-2][1]+TMap[-2][2]-1)/TMap[-2][2]
            TMap[-1] = [Nch+1,TMax,0]
            TMap = TMap.T
            TMap[0] -= 1
            return TMap.T,clockWd,mapNo
            
        def Tmap2TOF(TMap,clockWd):
            TOF = []
            Tch,T,Step = TMap[0]
            for tmap in TMap[1:]:
                tch,t,step = tmap
                TOF += [T+Step*(i-Tch) for i in range(Tch,tch)]
                Tch,T,Step = tmap
            TOF = np.array(TOF)*clockWd
            return TOF

        Banks = []
        Pos = []
        rdbuffer = kwarg.get('buffer')
        title = ''
        comments = None
        if '2' in platform.python_version_tuple()[0]:
            fp = open(filename,'Ur')
        else:
            fp = open(filename,'r',encoding='latin-1')
        # reload previously saved values - used for multibank reads
        if self.repeat and rdbuffer is not None:
            Banks = rdbuffer.get('Banks')
            Pos = rdbuffer.get('Pos')
            self.selections = rdbuffer.get('selections')
            comments = rdbuffer.get('comments')

        # read through the file and find the beginning of each bank
        # Save the offset (Pos), BANK line (Banks), comments for each bank
        #
        # This is going to need a fair amount of work to track line numbers
        # in the input file. 
        if len(Banks) != len(Pos) or len(Banks) == 0:
            i = -1
            while True:
                i += 1
                S = fp.readline()
                if len(S) == 0: break
                    
                if i==0: # first line is always a comment
                    self.errors = 'Error reading title'
                    title = S[:-1]
                    comments = [[title,]]
                    continue
                if i==1 and S[:4].lower() == 'inst' and ':' in S:
                    # 2nd line is instrument parameter file (optional)
                    self.errors = 'Error reading instrument parameter filename'
                    self.instparm = S.split(':')[1].strip('[]').strip()
                    continue
                if S[0] == '#': # allow comments anywhere in the file
                    # comments in fact should only preceed BANK lines
                    comments[-1].append(S[:-1])
                    continue
                if S[:4] == 'BANK':
                    self.errors = 'Error reading bank:'
                    self.errors += '  '+str(S)
                    comments.append([title,])
                    Banks.append(S)
                    Pos.append(fp.tell())
                if S[:8] == 'TIME_MAP':     #assumes one time map; HIPPO has multiple time maps
                    if len(Banks) == 0:
                        self.errors = 'Error reading time map before any bank lines'
                    else:
                        self.errors = 'Error reading time map after bank:\n  '+str(Banks[-1])
                    timemap,clockwd,mapNo = GetTimeMap(fp,fp.tell(),S)
                    self.TimeMap[mapNo] = timemap
                    self.clockWd[mapNo] = clockwd 
                    

        # Now select the bank to read
        if not Banks: # use of ContentsValidator should prevent this error
            print (self.formatName+' scan error: no BANK records')
            selblk = None # no block to choose
            self.errors = 'No BANK records found (strange!)'
            return False
        elif len(Banks) == 1: # only one Bank, don't ask
            selblk = 0
        elif self.repeat and self.selections is not None:
            # we were called to repeat the read
            #print 'debug: repeat #',self.repeatcount,'selection',self.selections[self.repeatcount]
            selblk = self.selections[self.repeatcount]
            self.repeatcount += 1
            if self.repeatcount >= len(self.selections): self.repeat = False
        else:                       # choose from options
            selblk = self.selections[0] # select first in list
            if len(self.selections) > 1: # prepare to loop through again
                self.repeat = True
                self.repeatcount = 1
                if rdbuffer is not None:
                    rdbuffer['Banks'] = Banks
                    rdbuffer['Pos'] = Pos
                    rdbuffer['selections'] = self.selections
                    rdbuffer['comments'] = comments

        # got a selection, now read it
        Bank = Banks[selblk]
        try:
            bnkNo = int(Bank.split()[1])
        except ValueError:
            bnkNo = 1
        if 'FXYE' in Bank:
            self.errors = 'Error reading FXYE data in Bank\n  '+Banks[selblk]
            self.powderdata = GetFXYEdata(fp,Pos[selblk],Banks[selblk])
        elif 'FXY' in Bank:
            self.errors = 'Error reading FXY data in Bank\n  '+Banks[selblk]
            self.powderdata = GetFXYdata(fp,Pos[selblk],Banks[selblk])
        elif 'ESD' in Bank:
            self.errors = 'Error reading ESD data in Bank\n  '+Banks[selblk]
            self.powderdata = GetESDdata(fp,Pos[selblk],Banks[selblk])
        elif 'STD' in Bank:
            self.errors = 'Error reading STD data in Bank\n  '+Banks[selblk]
            self.powderdata = GetSTDdata(fp,Pos[selblk],Banks[selblk])
        elif 'ALT' in Bank:
            self.errors = 'Error reading ALT data in Bank\n  '+Banks[selblk]
            self.powderdata = GetALTdata(fp,Pos[selblk],Banks[selblk])
        else:
            self.errors = 'Error reading STD data in Bank\n  '+Banks[selblk]
            self.powderdata = GetSTDdata(fp,Pos[selblk],Banks[selblk])

        self.errors = 'Error processing information after read complete'
        if comments is not None:
            self.comments = comments[selblk]
        self.powderentry[0] = filename
        self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = bnkNo #selblk+1 # bank number
        self.idstring = ospath.basename(filename) + ' Bank '+str(bnkNo) #selblk+1)
        self.numbanks=len(Banks)
        # scan comments for temperature & radius
        Temperature = 300.
        for S in self.comments:
            if 'Temp' in S.split('=')[0]:
                try:
                    Temperature = float(S.split('=')[1])
                except:
                    pass
            elif 'Gonio' in S.split('=')[0]:
                try:
                    self.Sample['Gonio. radius'] = float(S.split('=')[1])
                except:
                    pass
            elif 'Omega' in S.split('=')[0] or 'Theta' in S.split('=')[0]:  #HIPD weirdness
                try:
                    self.Sample['Omega'] = float(S.split('=')[1])
                except:
                    pass
            elif 'Chi' in S.split('=')[0]:
                try:
                    self.Sample['Chi'] = float(S.split('=')[1])
                except:
                    pass                    
            elif 'Phi' in S.split('=')[0]:
                try:
                    self.Sample['Phi'] = float(S.split('=')[1])
                except:
                    pass
        if 'EDS' in Bank:
            S = self.comments[0].lower().split('theta')
            if len(S) > 1:
                self.Inst['2-theta'] = [float(S[1]),float(S[1]),False]                    
        self.Sample['Temperature'] = Temperature
        fp.close()
        return True        

def sfloat(S):
    'convert a string to a float, treating an all-blank string as zero'
    if S.strip():
        return float(S)
    else:
        return 0.0

def sint(S):
    'convert a string to an integer, treating an all-blank string as zero'
    if S.strip():
        return int(S)
    else:
        return 0

