# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2pwd_GPX.py $
# $Id: G2pwd_GPX.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''
from __future__ import division, print_function
import platform
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5577 $")

def cPickleLoad(fp):
    if '2' in platform.python_version_tuple()[0]:
        return cPickle.load(fp)
    else:
        return cPickle.load(fp,encoding='latin-1')

class GSAS2_ReaderClass(G2obj.ImportPowderData):
    """Routines to import powder data from a GSAS-II file
    This should work to pull data out from a out of date .GPX file
    as long as the details of the histogram data itself don't change
    """
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.gpx',),
            strictExtension=True,
            formatName = 'GSAS-II gpx',
            longFormatName = 'GSAS-II project (.gpx file) import'
            )
        
    def ContentsValidator(self, filename):
        "Test if the 1st section can be read as a cPickle block, if not it can't be .GPX!"
        fp = open(filename,'rb')
        try: 
            data = cPickleLoad(fp)
        except:
            self.errors = 'This is not a valid .GPX file. Not recognized by cPickle'
            fp.close()
            return False
        fp.seek(0)
        nhist = 0
        while True:
            try:
                data = cPickleLoad(fp)
            except EOFError:
                break
            if data[0][0][:4] == 'PWDR':
                nhist += 1
                self.dnames.append(data[0][0])
        if nhist:
            if not len(self.selections):    #no PWDR entries
                self.selections = range(nhist)
                fp.close()
                return True
        self.errors = 'No powder entries found'
        fp.close()
        return False

    def Reader(self,filename, ParentFrame=None, **kwarg):
        '''Read a dataset from a .GPX file.
        If multiple datasets are requested, use self.repeat and buffer caching.
        '''
        histnames = []
        poslist = []
        rdbuffer = kwarg.get('buffer')
        # reload previously saved values
        if self.repeat and rdbuffer is not None:
            self.selections = rdbuffer.get('selections')
            poslist = rdbuffer.get('poslist')
            histnames = rdbuffer.get('histnames')
        else:
            try:
                fp = open(filename,'rb')
                while True:
                    pos = fp.tell()
                    try:
                        data = cPickleLoad(fp)
                    except EOFError:
                        break
                    if data[0][0][:4] == 'PWDR':
                        histnames.append(data[0][0])
                        poslist.append(pos)
            except:
                self.errors = 'Reading of histogram names failed'
                return False
            finally:
                fp.close()
        if not histnames:
            return False            # no blocks with powder data
        elif len(histnames) == 1: # one block, no choices
            selblk = 0
        elif self.repeat and self.selections is not None:
            # we were called to repeat the read
            #print 'debug: repeat #',self.repeatcount,'selection',selections[self.repeatcount]
            selblk = self.selections[self.repeatcount]
            self.repeatcount += 1
            if self.repeatcount >= len(self.selections): self.repeat = False
        else:                       # choose from options                
            selblk = self.selections[0] # select first in list
            if len(self.selections) > 1: # prepare to loop through again
                self.repeat = True
                self.repeatcount = 1
                if rdbuffer is not None:
                    rdbuffer['poslist'] = poslist
                    rdbuffer['histnames'] = histnames
                    rdbuffer['selections'] = self.selections

        fp = open(filename,'rb')
        fp.seek(poslist[selblk])
        data = cPickleLoad(fp)
        N = len(data[0][1][1][0])
        #self.powderdata = data[0][1][1]
        self.powderdata = [
            data[0][1][1][0], # x-axis values
            data[0][1][1][1], # powder pattern intensities
            data[0][1][1][2], # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # calc. background (zero)
            np.zeros(N), # obs-calc profiles
        ]
        self.powderentry[0] = filename
        self.powderentry[2] = selblk+1
        # pull some sections from the PWDR children
        for i in range(1,len(data)):
            if data[i][0] == 'Comments':
                self.comments = data[i][1]
                continue
            elif data[i][0] == 'Sample Parameters':
                self.Sample = data[i][1]
                continue
            for keepitem in ('Limits','Background','Instrument Parameters'): 
                if data[i][0] == keepitem:
                    self.pwdparms[data[i][0]] = data[i][1]
                    break
        self.idstring = data[0][0][5:]
        self.repeat_instparm = False # prevent reuse of iparm when several hists are read
        fp.close()
        return True
