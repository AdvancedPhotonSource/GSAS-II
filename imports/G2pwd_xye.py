# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2pwd_xye: Topas .xye data*
------------------------------------

Routine to read in powder data from a Topas-compatible .xye file

'''

import sys
import os.path as ospath
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
class xye_ReaderClass(G2IO.ImportPowderData):
    'Routines to import powder data from a .xye file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.xye',),
            strictExtension=False,
            formatName = 'Topas xye',
            longFormatName = 'Topas .xye powder data file'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filepointer):
        'Look through the file for expected types of lines in a valid Topas file'
        gotCcomment = False
        begin = True
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
                if S[0] == '#':
                    continue       #ignore comments, if any
                else:
                    begin = False
                # valid line to read? 
            vals = S.split()
            if len(vals) == 2 or len(vals) == 3:
                continue
            else:
                self.errors = 'Unexpected information in line: '+str(i+1)
                self.errors += '  '+str(S)
                return False
        return True # no errors encountered

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a Topas file'
        x = []
        y = []
        w = []
        try:
            gotCcomment = False
            begin = True
            for i,S in enumerate(filepointer):
                self.errors = 'Error reading line: '+str(i+1)
                # or a block of comments delimited by /* and */
                # or (GSAS style) each line can begin with '#'
                if begin:
                    if gotCcomment and S.find('*/') > -1:
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
                    else:
                        begin = False
                # valid line to read
                vals = S.split()
                if len(vals) < 2:
                    print 'Line '+str(i+1)+' cannot be read:\n\t'+S
                    continue
                try:
                    x.append(float(vals[0]))
                    f = float(vals[1])
                    if f <= 0.0:
                        y.append(0.0)
                        w.append(0.0)
                    elif len(vals) == 3:
                        y.append(float(vals[1]))
                        w.append(1.0/float(vals[2])**2)
                    else:
                        y.append(float(vals[1]))
                        w.append(1.0/float(vals[1]))
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

