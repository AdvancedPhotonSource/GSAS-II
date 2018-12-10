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
-------------------------------------

Routine to read in powder data from a FullProf .dat file

'''

from __future__ import division, print_function
import os.path as ospath
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1620 $")
class fp_ReaderClass(G2obj.ImportPowderData):
    'Routines to import powder data from a FullProf 1-10 column .dat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.dat',),
            strictExtension=False,
            formatName = 'FullProf .dat',
            longFormatName = 'FullProf 1-10 column .dat powder data file'
            )
        self.scriptable = True

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid FullProf file'
        gotCcomment = False
        begin = True
        self.GSAS = False
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if i > 50: break
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
                    float(val)
                    if j == 2:
                        break
            except ValueError:
                self.errors = 'Unexpected information in line: '+str(i+1)
                if all([ord(c) < 128 and ord(c) != 0 for c in str(S)]): # show only if ASCII
                    self.errors += '  '+str(S)
                else: 
                    self.errors += '  (binary)'
                if i > 2:
                    fp.close()
                    return False
        fp.close()
        return True # no errors encountered

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read a FullProf file'
        x = []
        y = []
        w = []
        gotCcomment = False
        begin = True
        steps = False
        Stop = False
        N = 0
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            self.errors = 'Error reading line: '+str(i+1)
            # Allow a block of comments delimited by /* and */
            # or (GSAS style) each comment line can begin with '#' or '!'
            if S.lstrip()[0] in ["'",'#','!',]:
                self.comments.append(S[:-1])
                continue       # store comments, if any
            if begin:
                if gotCcomment and S.find('*/') > -1:
                    self.comments.append(S[:-1])
                    begin = False
                    continue
                if S.strip().startswith('/*'):
                    self.comments.append(S[:-1])
                    gotCcomment = True
                    continue   
            # look for a line with start, steps etc. in 1st 4 lines
            if not steps:
                vals = S.replace(',',' ').split(None,4)
                if 'lambda' in S:
                    self.instdict['wave'] = float(vals[1])
                    continue
                elif len(vals) >= 3:
                    try:
                        start = float(vals[0])
                        step = float(vals[1])
                        stop = float(vals[2])
                        steps = True
                        begin = False
                        if len(vals) > 3:
                            self.comments.append(vals[3][:-1])
                    except:
                        print('Skipping line ',S)
                    continue
                elif i<3:
                    print('Skipping header line ',S)
                    continue
            # should be a valid line to read
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
                if GSASIIpath.GetConfigValue('debug'):
                    print (msg)
                    print (S.strip())
                break
            except:
                msg = 'Error in line '+str(i+1)
                if GSASIIpath.GetConfigValue('debug'):
                    print (msg)
                    print (S.strip())
                break
        fp.close()
        N = len(x)
        if N <= 1: return False
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
        #self.powderentry[2] = 1 # xye file only has one bank
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

