# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: vondreele $
# $Revision: 1620 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2pwd_xye.py $
# $Id: G2pwd_xye.py 1620 2014-12-27 17:14:59Z vondreele $
########### SVN repository information ###################
'''
*Module G2pwd_FP: FullProf .dat data*
------------------------------------

Routine to read in powder data from a FullProf .dat file

'''

import sys
import os.path as ospath
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1620 $")
class xye_ReaderClass(G2IO.ImportPowderData):
    'Routines to import powder data from a FullProf 1-10 column .dat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.dat',),
            strictExtension=False,
            formatName = 'FullProf .dat',
            longFormatName = 'FullProf 1-10 column .dat powder data file'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filepointer):
        'Look through the file for expected types of lines in a valid FullProf file'
        gotCcomment = False
        begin = True
        steps = False
        self.GSAS = False
        for i,S in enumerate(filepointer):
            if i > 1000: break
            if begin:
                if gotCcomment and S.find('*/') > -1:
                    begin = False
                    continue
                if S.strip().startswith('/*'):
                    gotCcomment = True
                    continue   
                if S.lstrip()[0] in ["'",'#','!',]:
                    continue       #ignore comments, if any
                else:
                    begin = False
                # valid line to read? 
            vals = S.split()
            try:    #look for start,step,stop card
                for j,val in enumerate(vals):
                    num = float(val)
                    if j == 2:
                        break
            except ValueError:
                self.errors = 'Unexpected information in line: '+str(i+1)
                if all([ord(c) < 128 and ord(c) != 0 for c in str(S)]): # show only if ASCII
                    self.errors += '  '+str(S)
                else: 
                    self.errors += '  (binary)'
                return False
        return True # no errors encountered

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a FullProf file'
        x = []
        y = []
        w = []
        try:
            gotCcomment = False
            begin = True
            steps = False
            Stop = False
            N = 0
            for i,S in enumerate(filepointer):
                self.errors = 'Error reading line: '+str(i+1)
                # or a block of comments delimited by /* and */
                # or (GSAS style) each line can begin with '#' or '!'
                if begin:
                    if gotCcomment and S.find('*/') > -1:
                        self.comments.append(S[:-1])
                        begin = False
                        continue
                    if S.strip().startswith('/*'):
                        self.comments.append(S[:-1])
                        gotCcomment = True
                        continue   
                    if S.lstrip()[0] in ["'",'#','!',]:
                        self.comments.append(S[:-1])
                        continue       #ignore comments, if any
                    else:
                        begin = False
                # valid line to read
                if not steps:
                    vals = S.split(None,4)
                    if len(vals) >= 3:
                        steps = True
                        start = float(vals[0])
                        step = float(vals[1])
                        stop = float(vals[2])
                        if len(vals) > 3:
                            self.comments.append(vals[3][:-1])
                        continue
                vals = S.split()    #data strings
                try:
                    for j in range(len(vals)):
                        x.append(start+N*step)
                        f = float(vals[j])
                        if f <= 0.0:
                            y.append(0.0)
                            w.append(0.0)
                        else:
                            y.append(float(vals[j]))
                            w.append(1.0/float(vals[j]))
                        if x[-1] >= stop:
                            Stop = True
                            break
                        N += 1
                    if Stop:
                        break
                except ValueError:
                    msg = 'Error parsing number in line '+str(i+1)
                    print msg
                    print S
                    break
                except:
                    msg = 'Error in line '+str(i+1)
                    print msg
                    print S
                    break
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
            #self.powderentry[1] = pos # bank offset (N/A here)
            self.powderentry[2] = 1 # xye file only has one bank
            self.idstring = ospath.basename(filename)
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
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

