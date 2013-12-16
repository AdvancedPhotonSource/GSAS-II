# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2pwd_GPX: GSAS-II projects*
------------------------------------
Routine to import powder data from GSAS-II .gpx files

'''
import sys
import cPickle
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

class GSAS2_ReaderClass(G2IO.ImportPowderData):
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
        
    def ContentsValidator(self, filepointer):
        "Test if the 1st section can be read as a cPickle block, if not it can't be .GPX!"
        try: 
            cPickle.load(filepointer)
        except:
            self.errors = 'This is not a valid .GPX file. Not recognized by cPickle'
            return False
        return True

    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        '''Read a dataset from a .GPX file.
        If multiple datasets are requested, use self.repeat and buffer caching.
        '''
        histnames = []
        poslist = []
        rdbuffer = kwarg.get('buffer')
        # reload previously saved values
        if self.repeat and rdbuffer is not None:
            selections = rdbuffer.get('selections')
            poslist = rdbuffer.get('poslist')
            histnames = rdbuffer.get('histnames')
        else:
            try:
                fl = open(filename,'rb')
                while True:
                    pos = fl.tell()
                    try:
                        data = cPickle.load(fl)
                    except EOFError:
                        break
                    if data[0][0][:4] == 'PWDR':
                        histnames.append(data[0][0])
                        poslist.append(pos)
            except:
                self.errors = 'Reading of histogram names failed'
                return False
            finally:
                fl.close()
        if not histnames:
            return False            # no blocks with powder data
        elif len(histnames) == 1: # one block, no choices
            selblk = 0
        elif self.repeat and selections is not None:
            # we were called to repeat the read
            #print 'debug: repeat #',self.repeatcount,'selection',selections[self.repeatcount]
            selblk = selections[self.repeatcount]
            self.repeatcount += 1
            if self.repeatcount >= len(selections): self.repeat = False
        else:                       # choose from options                
            selections = self.MultipleBlockSelector(
                histnames,
                ParentFrame=ParentFrame,
                title='Select histogram(s) to read from the list below',
                size=(600,100),
                header='Dataset Selector')
            if len(selections) == 0:
                self.errors = 'No histogram selected'
                return False
            selblk = selections[0] # select first in list
            if len(selections) > 1: # prepare to loop through again
                self.repeat = True
                self.repeatcount = 1
                if rdbuffer is not None:
                    rdbuffer['poslist'] = poslist
                    rdbuffer['histnames'] = histnames
                    rdbuffer['selections'] = selections

        try:
            fl = open(filename,'rb')
            fl.seek(poslist[selblk])
            data = cPickle.load(fl)
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
            self.idstring = data[0][0][5:]
            # pull out wavelength 
            try:
                if len(data[4][1]) == 2: # current GPX file
                    if data[4][1][0].get('Lam'):
                        self.instdict['wave'] = [data[4][1][0].get('Lam')[1]]
                    elif data[4][1][0].get('Lam1') and data[4][1][0].get('Lam2'):
                        self.instdict['wave'] = [
                            data[4][1][0].get('Lam1')[1],
                            data[4][1][0].get('Lam2')[1]
                            ]
                elif len(data[4][1]) == 4: # original GPX file
                    pos = data[4][1][3].index('Lam')
                    self.instdict['wave'] = [data[4][1][1][pos],]
            except:
                self.warnings += "Failed to read wavelength"
                self.warnings += '\n  '+str(detail)
            # pull out temperature
            try:
                if data[5][1].get('Temperature'):
                    self.Sample['Temperature'] = data[5][1]['Temperature']
            except:
                self.warnings += "Failed to read temperature"
                self.warnings += '\n  '+str(detail)
            self.repeat_instparm = False # prevent reuse of iparm when several hists are read
            return True
        except Exception as detail:
            self.errors = 'Error reading selected histogram\n  '+str(detail)
            print self.formatName+' error:',detail # for testing
            print sys.exc_info()[0] # for testing
            import traceback
            print traceback.format_exc()
            return False
        finally:
            fl.close()
