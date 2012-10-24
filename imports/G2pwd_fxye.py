# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2012-02-13 11:33:35 -0600 (Mon, 13 Feb 2012) $
# $Author: vondreele & toby $
# $Revision: 482 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/G2importphase.py $
# $Id: G2importphase.py 482 2012-02-13 17:33:35Z vondreele $
########### SVN repository information ###################
# a routine to read in powder data from a GSAS-compatible files
# 
import sys
import os.path as ospath
import numpy as np
import GSASIIIO as G2IO

class GSAS_ReaderClass(G2IO.ImportPowderData):
    'Routines to import powder data from a GSAS files'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.fxye','.raw','.gsas','.gsa','.RAW','.GSAS','.GSA'),
            strictExtension=False,
            formatName = 'GSAS',
            longFormatName = 'GSAS powder data files'
            )

    # Validate the contents -- look for a bank line
    def ContentsValidator(self, filepointer):
        #print 'ContentsValidator: '+self.formatName
        for i,line in enumerate(filepointer):
            self.TimeMap = False
            if i==0: # first line is always a comment
                continue
            if i==1 and line[:4].lower() == 'inst' and ':' in line:
                # 2nd line is optional instrument parameter file
                continue
            if line[0] == '#': continue
            if line[:4] == 'BANK':
                return True
            if line [:8] == 'TIME_MAP':          #LANSCE TOF data
                self.TimeMap = True
                return True
            else:
                print 'ContentsValidator: '+self.formatName
                print 'Unexpected information in line:',i+1 # debug info
                print line
                return False
        return False # no bank records

    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        clockWd = None
        
        def GetFXYEdata(File,Pos,Bank):
            File.seek(Pos)
            x = []
            y = []
            w = []
            S = File.readline()
            while S and S[:4] != 'BANK':
                vals = S.split()
                x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
                f = float(vals[1])
                if f <= 0.0:
                    y.append(0.0)
                    w.append(1.0)
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
            while S and S[:4] != 'BANK':
                vals = S.split()
                x.append(float(vals[0])/100.)               #CW: from centidegrees to degrees
                f = float(vals[1])
                if f > 0.0:
                    y.append(f)
                    w.append(1.0/f)
                else:              
                    y.append(0.0)
                    w.append(1.0)
                S = File.readline()
            N = len(x)
            return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            
        def GetESDdata(File,Pos,Bank):
            File.seek(Pos)
            cons = Bank.split()
            if clockWd:
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
            N = len(x)
            if clockWd:
                x,cw = Tmap2TOF(self.TimeMap,clockWd)
                x = x[:-2]+cw[:-1]/2.
                return [x,np.array(y)/cw[:-1],np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            else:
                return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
        
        def GetSTDdata(File,Pos,Bank):
            File.seek(Pos)
            cons = Bank.split()
            Nch = int(cons[2])
            if clockWd:
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
            N = len(x)
            if clockWd:
                x,cw = Tmap2TOF(self.TimeMap,clockWd)
                x = x[:-2]+cw[:-1]/2.
                return [x,np.array(y)/cw[:-1],np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            else:
                return [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
            
        def GetTimeMap(File,Pos,TimeMap):
            File.seek(Pos)
            cons = TimeMap.split()
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
            return TMap.T,clockWd
            
        def Tmap2TOF(TMap,clockWd):
            TOF = []
            chWdt = []
            Tch,T,Step = TMap[0]
            for tmap in TMap[1:]:
                tch,t,step = tmap
                TOF += [T+Step*(i-Tch) for i in range(Tch,tch)]
                Tch,T,Step = tmap
            TOF = np.array(TOF)*clockWd
            chWdt = np.diff(TOF)
            return TOF,chWdt

        x = []
        y = []
        w = []
        Banks = []
        Pos = []
        rdbuffer = kwarg.get('buffer')
        title = ''
        comments = None
        selections = None

        # reload previously saved values
        if self.repeat and rdbuffer is not None:
            Banks = rdbuffer.get('Banks')
            Pos = rdbuffer.get('Pos')
            selections = rdbuffer.get('selections')
            comments = rdbuffer.get('comments')

        # read through the file and find the beginning of each bank
        # Save the offset (Pos), BANK line (Banks), comments for each
        # bank
        if len(Banks) != len(Pos) or len(Banks) == 0:
            try:
                i = -1
                while True:
                    i += 1
                    S = filepointer.readline()
                    if len(S) == 0: break
                        
                    if i==0: # first line is always a comment
                        title = S[:-1]
                        comments = [[title,]]
                        continue
                    if i==1 and S[:4].lower() == 'inst' and ':' in S:
                        # 2nd line is instrument parameter file (optional)
                        self.instparm = S.split(':')[1].strip()
                        continue
                    if S[0] == '#': # allow comments anywhere in the file
                        # comments in fact should only preceed BANK lines
                        comments[-1].append(S[:-1])
                        continue
                    if S[:4] == 'BANK':
                        comments.append([title,])
                        Banks.append(S)
                        Pos.append(filepointer.tell())
                    if S[:8] == 'TIME_MAP':
                        self.TimeMap,clockWd = GetTimeMap(filepointer,filepointer.tell(),S)
            except Exception as detail:
                print self.formatName+' scan error:'+str(detail) # for testing
                import traceback
                traceback.print_exc(file=sys.stdout)
                return False

        # Now select the bank to read
        if not Banks: # use of ContentsValidator should prevent this error
            print self.formatName+' scan error: no BANK records'
            selblk = None # no block to choose
            return False
        elif len(Banks) == 1: # only one Bank, don't ask
            selblk = 0
        elif self.repeat and selections is not None:
            # we were called to repeat the read
            print 'debug: repeat #',self.repeatcount,'selection',selections[self.repeatcount]
            selblk = selections[self.repeatcount]
            self.repeatcount += 1
            if self.repeatcount >= len(selections): self.repeat = False
        else:                       # choose from options
            selections = self.MultipleBlockSelector(
                Banks,
                ParentFrame=ParentFrame,
                title='Select Bank(s) to read from the list below',
                size=(600,100),
                header='Dataset Selector')
            if len(selections) == 0: return False
            selblk = selections[0] # select first in list
            if len(selections) > 1: # prepare to loop through again
                self.repeat = True
                self.repeatcount = 1
                if rdbuffer is not None:
                    rdbuffer['Banks'] = Banks
                    rdbuffer['Pos'] = Pos
                    rdbuffer['selections'] = selections
                    rdbuffer['comments'] = comments

        # got a selection, now read it
        Bank = Banks[selblk]
        try:
            if 'FXYE' in Bank:
                self.powderdata = GetFXYEdata(filepointer,
                                              Pos[selblk],Banks[selblk])
            elif 'FXY' in Bank:
                self.powderdata = GetFXYdata(filepointer,
                                             Pos[selblk],Banks[selblk])
            elif 'ESD' in Bank:
                self.powderdata = GetESDdata(filepointer,
                                             Pos[selblk],Banks[selblk])
            elif 'STD' in Bank:
                self.powderdata = GetSTDdata(filepointer,
                                             Pos[selblk],Banks[selblk])
            else:
                self.powderdata = GetSTDdata(filepointer,
                                             Pos[selblk],Banks[selblk])
        except Exception as detail:
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False
        if comments is not None:
            self.comments = comments[selblk]
        self.powderentry[0] = filename
        self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = selblk+1 # bank number
        self.idstring = ospath.basename(filename) + ' Bank '+str(selblk+1)
        self.numbanks=len(Banks)
        # scan comments for temperature
        Temperature = 300
        for S in self.comments:
            if 'Temp' in S.split('=')[0]:
                try:
                    Temperature = float(S.split('=')[1])
                except:
                    pass
        self.Sample['Temperature'] = Temperature
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
