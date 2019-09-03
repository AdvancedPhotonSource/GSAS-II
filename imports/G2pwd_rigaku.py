# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: toby $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################

from __future__ import division, print_function
import os
import platform
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class Rigaku_txtReaderClass(G2obj.ImportPowderData):
    '''Routines to import powder data from a Rigaku .txt file with an angle and
    then 1 or 11(!) intensity values on the line. The example file is proceeded
    with 10 of blank lines, but I have assumed they could be any sort of text.
    This code should work with an angle and any number of intensity values/line
    as long as the number is the same on each line. The step size may not change. The
    number of comment lines can also change, but should not appear to be intensity
    values (numbers only).
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.txt','.TXT'),
            strictExtension=True,
            formatName = 'Rigaku .txt exported',
            longFormatName = 'Rigaku powder data exported as .txt'
            )
        self.scriptable = True
        self.vals = None
        self.stepsize = None
        self.skip = 0

    # Validate the contents -- make sure we only have valid lines and set
    # values we will need for later read.
    def ContentsValidator(self, filename):
        self.vals = None
        self.stepsize = None
        j = 0
        prevAngle = None
        header = True
        self.skip = -1
        fp = open(filename,'r')
        for i,line in enumerate(fp):
            sline = line.split()
            vals = len(sline)
            if header:
                self.skip += 1
                if not line.strip(): continue # ignore blank lines
                err = False
                for item in sline:
                    try:
                        float(item)
                    except:
                        err = True
                        break
                if err: continue
                if vals < 1: continue
                header = False # found first non-header line
            if vals < 2:
                print('Too few values for Rigaku .txt file')
                fp.close()
                return False
            if self.vals is None:
                self.vals = vals
            elif self.vals != vals:
                print('Inconsistent numbers values for Rigaku .txt file on line '+str(i+1))
                fp.close()
                return False
            else:
                j += 1
            try: 
                angle = float(sline[0])
            except:
                print('Unable to read angle on line '+str(i+1))
                fp.close()
                return False
            if prevAngle is None:
                prevAngle = angle
                continue
            stepsize = (angle-prevAngle)/(vals-1)
            prevAngle = angle
            if self.stepsize is None:
                self.stepsize = stepsize
            elif abs(self.stepsize - stepsize) > max(abs(stepsize),abs(self.stepsize))/10000. :
                print('Inconsistent step size for Rigaku .txt file on line '+
                        str(i+1) + ' here '+ repr(stepsize) + ' prev '+ repr(self.stepsize))
                fp.close()
                return False
            if j > 30:
                fp.close()
                return True
        fp.close()
        return False
            
    def Reader(self,filename, ParentFrame=None, **kwarg):
        'Read a Rigaku .txt file'
        x = []
        y = []
        w = []
        fp = open(filename,'r')
        for i,line in enumerate(fp):
            if i < self.skip: continue
            sline = line.split()
            try: 
                angle = float(sline[0])
            except:
                print('Unable to read angle on line '+str(i+1))
                self.errors = 'Error reading line: '+str(i+1)
                return False
            for j in sline[1:]:
                x.append(angle)
                angle += self.stepsize
                try: 
                    y.append(float(j))
                except:
                    print('Unable to read intensity on line '+str(i+1))
                    self.errors = 'Error reading line: '+str(i+1)
                    return False
                w.append(1.0/max(1.,float(j)))
        N = len(x)
        self.powderdata = [
            np.array(x), # x-axis values
            np.array(y), # powder pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # calc. background (zero)
            np.zeros(N), # obs-calc profiles
            ]
        self.powderentry[0] = filename
        #self.powderentry[2] = 1 # xye file only has one bank
        self.idstring = os.path.basename(filename)
        return True
class Rigaku_rasReaderClass(G2obj.ImportPowderData):
    '''Routines to import powder data from a Rigaku .ras file with multiple scans. 
    All scans will be imported as individual PWDR entries
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.ras','.RAS'),
            strictExtension=True,
            formatName = 'Rigaku .ras file',
            longFormatName = 'Rigaku .ras raw multipattern powder data'
            )
        self.scriptable = True
        self.vals = None
        self.stepsize = None
        self.skip = 0

    # Validate the contents -- make sure we only have valid lines and set
    # values we will need for later read.
    def ContentsValidator(self, filename):
        if '2' in platform.python_version_tuple()[0]:
            fp = open(filename,'Ur')
        else:
            fp = open(filename,'r',encoding='latin-1')
        self.vals = None
        self.stepsize = None
        fp.seek(0)
        if fp.readline()[:-1] != '*RAS_DATA_START':
            self.errors = 'Bad ras file'
            fp.close()
            return False
        nBanks= 0
        for i,line in enumerate(fp):
            if line[:-1] == '*RAS_HEADER_START':
                nBanks += 1
                self.dnames.append(os.path.basename(filename)+' sample '+(str(nBanks)))
        if nBanks:
            if not len(self.selections):
                self.selections = list(range(nBanks))
                self.numbanks = nBanks
        fp.close()
        return True

    def Reader(self,filename, ParentFrame=None, **kwarg):
        'Read a Rigaku .ras file'
        if '2' in platform.python_version_tuple()[0]:
            fp = open(filename,'Ur')
        else:
            fp = open(filename,'r',encoding='latin-1')
        blockNum = self.selections[0]
        x = []
        y = []
        w = []
        block = 0
        while True:
            line = fp.readline()[:-1]
            if line != '*RAS_INT_START':
                continue
            if block == blockNum:
                line = fp.readline()[:-1]
                while True:
                    if line == '*RAS_INT_END':
                        break
                    sline = line.split()
                    x.append(float(sline[0]))
                    y.append(float(sline[1]))
                    w.append(1.0/max(1.,float(y[-1])))
                    line = fp.readline()[:-1]
                break
            block += 1            
        N = len(x)
        self.powderdata = [
            np.array(x), # x-axis values
            np.array(y), # powder pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # calc. background (zero)
            np.zeros(N), # obs-calc profiles
            ]
        self.powderentry[0] = self.dnames[blockNum]
        self.idstring = self.dnames[blockNum]
        self.selections.remove(blockNum)
        self.repeat = False
        if len(self.selections):
            self.repeat = True
        return True
