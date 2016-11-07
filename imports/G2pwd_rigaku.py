# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: toby $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
import sys
import os
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: $")
class Rigaku_ReaderClass(G2IO.ImportPowderData):
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
        self.vals = None
        self.stepsize = None
        self.skip = 0

    # Validate the contents -- make sure we only have valid lines and set
    # values we will need for later read.
    def ContentsValidator(self, filepointer):
        self.vals = None
        self.stepsize = None
        j = 0
        prevAngle = None
        header = True
        self.skip = -1
        for i,line in enumerate(filepointer):
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
                return False
            if self.vals is None:
                self.vals = vals
            elif self.vals != vals:
                print('Inconsistent numbers values for Rigaku .txt file on line '+
                      str(i+1))
                return False
            else:
                j += 1
            try: 
                angle = float(sline[0])
            except:
                print('Unable to read angle on line '+str(i+1))
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
                return False
            if j > 30: return True
        return False
            
    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        'Read a Rigaku .txt file'
        x = []
        y = []
        w = []
        for i,line in enumerate(filepointer):
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
