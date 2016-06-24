# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2014-12-27 11:14:59 -0600 (Sat, 27 Dec 2014) $
# $Author: vondreele $
# $Revision: 1620 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2pwd_xye.py $
# $Id: G2pwd_xye.py 1620 2014-12-27 17:14:59Z vondreele $
########### SVN repository information ###################
'''
*Module G2pwd_Pnn: GSAS .Pnn data*
------------------------------------

Routine to read in powder data from a GSAS .Pnn binary blocked file

'''

import sys
import os.path
import numpy as np
import array as ar
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1620 $")
class Pnn_ReaderClass(G2IO.ImportPowderData):
    'Routines to import powder data from a .Pnn file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.P01','.P02','.P03','.P04','.P05','.P06',
                           '.P07','.P08','.P09','.P*',),
            strictExtension=False,
            formatName = 'GSAS .Pnn',
            longFormatName = 'GSAS .Pnn powder data file'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filepointer):
        'Look through the file for expected types of lines in a valid GSAS Pnn file'
        gotCcomment = False
        self.GSAS = False
        # file extension must be .Pxx or .pxx
        ext = os.path.splitext(filepointer.name)[1]
        if ext[1].upper() != 'P' or len(ext) != 4:
            return False
        if ext[2].isdigit() and ext[3].isdigit():
            return True # no errors encountered
        return False
        

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read a GSAS Pnn file'
        x = []
        y = []
        w = []
        buf = filepointer.read()
        print len(buf)
        remain = len(buf)%512
        if remain:
            buf += (512-remain)*' '
        print len(buf)
        buf = np.array(ar.array('f',buf),dtype=np.float32)
        print buf.shape
        buf = np.reshape(buf,(-1,128))
        print buf.shape
        for block in buf:
            block = np.reshape(block[:126],(14,9)).T
            x += list(block[1]-block[7]/2.)
            y += list(block[2]*block[4])
            w += list(block[6])
        N = len(x)
        print N
        self.powderdata = [np.array(x),np.array(y),np.array(w),np.zeros(N),np.zeros(N),np.zeros(N)]
        self.powderentry[0] = filename
        self.powderentry[2] = 1 # Pnn file only has one bank
        self.idstring = os.path.basename(filename)
        return True
     
            
            

