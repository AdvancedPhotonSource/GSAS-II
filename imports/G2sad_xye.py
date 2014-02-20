# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2013-12-16 10:43:01 -0600 (Mon, 16 Dec 2013) $
# $Author: toby $
# $Revision: 1168 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2sad_xye.py $
# $Id: G2sad_xye.py 1168 2013-12-16 16:43:01Z toby $
########### SVN repository information ###################
'''
*Module G2sad_xye: small angle q step .xye data*
------------------------------------

Routine to read in small angle data from an .xye file

'''

import sys
import os.path as ospath
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1168 $")
npasind = lambda x: 180.*np.arcsin(x)/np.pi

class xye_ReaderClass(G2IO.ImportSmallAngleData):
    'Routines to import q SAXD data from a .xye file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.xye','.txt','.dat'),
            strictExtension=False,
            formatName = 'q step xye',
            longFormatName = 'q stepped data file'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filepointer):
        'Look through the file for expected types of lines in a valid q-step file'
        gotCcomment = False
        begin = True
        for i,S in enumerate(filepointer):
            if i > 1000: break
            if begin:
                if gotCcomment and S.find('*/') > -1:
                    begin = False
                    continue
                if S.strip().startswith('/*'):
                    gotCcomment = True
                    continue   
                if S[0] == '#':
                    continue       #ignore comments, if any
                elif len(S) == 1:
                    continue        #ignore blank lines
                elif not gotCcomment:
                    begin = False
                # valid line to read?
            if begin:
                continue    
            vals = S.split()
            if 2 <= len(vals) <= 3:
                continue
            else:
                self.errors = 'Unexpected information in line: '+str(i+1)
                self.errors += '  '+str(S)
                return False
        return True # no errors encountered

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a q-step file'
        x = []
        y = []
        w = []
        try:
            wave = 1.5428   #Cuka default
            Temperature = 300
            gotCcomment = False
            begin = True
            for i,S in enumerate(filepointer):
                self.errors = 'Error reading line: '+str(i+1)
                # or a block of comments delimited by /* and */
                # or (GSAS style) each line can begin with '#'
                if begin:
                    if 'Wave' in S.split('=')[0]:
                        try:
                            wave = float(S.split('=')[1])
                        except:
                            pass
                    if gotCcomment:
                        if S.find('*/') == -1:
                            self.comments.append(S[:-1])
                        else:
                            self.comments.append(S[:-1])
                            begin = False
                        continue
                    if S.strip().startswith('/*'):
                        self.comments.append(S[:-1])
                        gotCcomment = True
                        continue   
                    if S[0] == '#':
                        self.comments.append(S[:-1])
                        continue       #ignore comments, if any
                    elif len(S) == 1:      #blank line only CR/LF
                        continue
                    elif  not gotCcomment:
                        begin = False
                # valid line to read
                if begin:
                    continue
                vals = S.split()
                try:
                    x.append(float(vals[0]))
#                    x.append(2.*npasind(wave*float(vals[0])/(4.*np.pi)))
                    f = float(vals[1])
                    if f <= 0.0:
                        y.append(0.0)
                        w.append(1.0)
                    elif len(vals) == 3:
                        y.append(float(vals[1]))
                        w.append(1.0/float(vals[2])**2)
                    else:
                        y.append(float(vals[1]))
                        w.append(1.0/float(vals[1]))
                except ValueError:
                    msg = 'Error in line '+str(i+1)
                    print msg
                    break
            N = len(x)
            for S in self.comments:
                if 'Temp' in S.split('=')[0]:
                    try:
                        Temperature = float(S.split('=')[1])
                    except:
                        pass
            self.instdict['wave'] = wave
            self.instdict['type'] = 'LXC'
            x = np.array(x)
            if np.any(x > 2.):      #nanometers-1?
                x /= 10.
            self.smallangledata = [
                x, # x-axis values - q
                np.array(y), # small angle pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.zeros(N), # calc. intensities (zero)
                np.zeros(N), # obs-calc profiles
                ]
            self.smallangleentry[0] = filename
            self.smallangleentry[2] = 1 # xye file only has one bank
            self.idstring = ospath.basename(filename)
            # scan comments for temperature
            self.Sample['Temperature'] = Temperature

            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

class txt_ReaderClass(G2IO.ImportSmallAngleData):
    'Routines to import q SAXD data from a .txt file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.txt','.dat'),
            strictExtension=False,
            formatName = 'q step txt',
            longFormatName = 'q stepped text data file'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filepointer):
        'Look through the file for expected types of lines in a valid txt q-step file'
        Ndata = 0
        for i,S in enumerate(filepointer):
            vals = S.split()
            if 2 <= len(vals):
                try:
                    data = [float(val) for val in vals]
                    Ndata += 1
                except ValueError:
                    pass
        if not Ndata:     
            self.errors = 'No 2 or more column numeric data found'
            return False
        return True # no errors encountered

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a q-step file'
        x = []
        y = []
        w = []
        try:
            wave = 1.5428   #Cuka default
            Temperature = 300
            Ndata = 0
            for i,S in enumerate(filepointer):
                if len(S) == 1:     #skip blank line
                    continue
                if '=' in S:
                    self.comments.append(S[:-1])
                    if 'wave' in S.split('=')[0].lower():
                        try:
                            wave = float(S.split('=')[1])
                        except:
                            pass
                    continue
                vals = S.split()
                if 2 <= len(vals):
                    try:
                        data = [float(val) for val in vals]
                        x.append(float(data[0]))
                        f = float(data[1])
                        if f <= 0.0:
                            y.append(0.0)
                            w.append(1.0)
                        elif len(vals) > 2:
                            y.append(float(data[1]))
                            w.append(1.0/float(data[2])**2)
                        else:
                            y.append(float(data[1]))
                            w.append(1.0/float(data[1]))
                    except ValueError:
                        msg = 'Error in line '+str(i+1)
                        print msg
                        continue
            N = len(x)
            for S in self.comments:
                if 'Temp' in S.split('=')[0]:
                    try:
                        Temperature = float(S.split('=')[1])
                    except:
                        pass
            self.instdict['wave'] = wave
            self.instdict['type'] = 'LXC'
            x = np.array(x)
            if np.any(x > 2.):         #q must be nm-1
                x /= 10.
#            x = 2.*npasind(wave*x)/(4.*np.pi)   #convert to 2-theta
            self.smallangledata = [
                x, # x-axis values q
                np.array(y), # small angle pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.zeros(N), # calc. intensities (zero)
                np.zeros(N), # obs-calc profiles
                ]
            self.smallangleentry[0] = filename
            self.smallangleentry[2] = 1 # xye file only has one bank
            self.idstring = ospath.basename(filename)
            # scan comments for temperature
            self.Sample['Temperature'] = Temperature
    
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False
