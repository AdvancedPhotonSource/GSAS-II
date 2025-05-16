# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
import os.path as ospath
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIIpath
class csv_ReaderClass(G2obj.ImportPowderData):
    'Routines to import powder data from a .xye file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.csv','.xy','.XY',),
            strictExtension=True,
            formatName = 'comma/tab/semicolon separated',
            longFormatName = 'Worksheet-type .csv powder data file'
            )
        self.scriptable = True

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        good = 0
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if S.strip().startswith('#'): continue
            if i > 1000: break
            vals = S.replace(',',' ').replace(';',' ').split()
            if len(vals) >= 2:
                for j,v in enumerate(vals):
                    if len(v) == 0: continue
                    if j == 3: break
                    try:
                        float(v)
                    except ValueError:
                        if good > 1: 
                            fp.close()
                            return False
                        continue
                good += 1
                continue
            elif good > 1:
                fp.close()
                return False
        fp.close()
        return True # no errors encountered

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read a csv file'
        x = []
        y = []
        w = []
        positions = [0,1,2]
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if i <= 2 and 'x=' in S: # header entry specifying columns
                for v in S.strip().replace(',',';').split(';'):
                    if 'x=' in v:
                        j = 0
                    elif 'y=' in v:
                        j = 1
                    elif 'e=' in v:
                        j = 2
                    else:
                        continue
                    try:
                        positions[j] = int(v.strip().split('=')[1])
                        if j == 1: positions[2] = -1
                    except:
                        print('Error parsing: "'+S+'"')
                print('positions=',positions)
            if S.strip().startswith('#'): continue
            vals = S.replace(',',' ').replace(';',' ').split()
            if len(vals) < 2 and i > 0:
                print ('Line '+str(i+1)+' cannot be read:\n\t'+S)
                continue
            try:
                x.append(float(vals[positions[0]]))
                f = float(vals[positions[1]])
                if f <= 0.0:
                    y.append(0.0)
                    w.append(0.0)
                elif len(vals) > positions[2] and positions[2] >= 0:
                    y.append(f)
                    w.append(1.0/float(vals[positions[2]])**2)
                else:
                    y.append(f)
                    w.append(1.0/f)
                err = False
            except ValueError:
                err = True
                msg = 'Error parsing number in line '+str(i+1)
            except:
                err = True
                msg = 'Error in line '+str(i+1)
            if err and i > 0:
                if GSASIIpath.GetConfigValue('debug'):
                    print (msg)
                    print (S.strip())
                break
        fp.close()
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
