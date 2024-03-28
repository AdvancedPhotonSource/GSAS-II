#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 14:22:54 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5576 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/exports/G2export_pwdr.py $
# $Id: G2export_pwdr.py 5576 2023-05-11 19:22:54Z toby $
########### SVN repository information ###################
'''Classes in :mod:`G2export_pwdr` follow:
'''
from __future__ import division, print_function
import os.path
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5576 $")
import GSASIIIO as G2IO
import GSASIIobj as G2obj
import GSASIIfiles as G2fil

class ExportPowderFXYE(G2IO.ExportBaseclass):
    '''Used to create a FXYE file for a powder data set

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'GSAS FXYE file',
            extension='.fxye',
            longFormatName = 'Export powder data as GSAS FXYE (column) file'
            )
        self.exporttype = ['powder']
        self.multiple = True

    def WriteInstFile(self,hist,Inst):
        '''Write an instrument parameter file
        '''
        if 'T' in Inst['Type'][0] or 'B' in Inst['Type'][0]:      #can't do TOF or pink iparm files
            return None
        prmname = os.path.splitext(self.filename)[0] + '.prm'
        prmname = os.path.join(self.dirname,prmname)
        self.OpenFile(prmname)
        self.Write( '            123456789012345678901234567890123456789012345678901234567890        ')
        self.Write( 'INS   BANK      1                                                               ')
        self.Write(('INS   HTYPE   %sR                                                              ')%(Inst['Type'][0]))
        if 'Lam1' in Inst:              #Ka1 & Ka2
            self.Write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   ')%(Inst['Lam1'][0],Inst['Lam2'][0]))
        elif 'Lam' in Inst:             #single wavelength
            self.Write(('INS  1 ICONS%10.7f%10.7f    0.0000               0.990    0     0.500   ')%(Inst['Lam'][1],0.0))
        self.Write( 'INS  1 IRAD     0                                                               ')
        self.Write( 'INS  1I HEAD                                                                    ')
        self.Write( 'INS  1I ITYP    0    0.0000  180.0000         1                                 ')
        self.Write(('INS  1DETAZM%10.3f                                                          ')%(Inst['Azimuth'][0]))
        self.Write( 'INS  1PRCF1     3    8   0.00100                                                ')
        self.Write(('INS  1PRCF11%15.6e%15.6e%15.6e%15.6e   ')%(Inst['U'][1],Inst['V'][1],Inst['W'][1],0.0))
        self.Write(('INS  1PRCF12%15.6e%15.6e%15.6e%15.6e   ')%(Inst['X'][1],Inst['Y'][1],Inst['SH/L'][1]/2.,Inst['SH/L'][1]/2.))
        self.CloseFile()
        print('Parameters from '+hist+' written to '+prmname)
        return prmname

    def Writer(self,TreeName,filename=None,prmname=''):
        '''Write a single PWDR entry to a FXYE file
        '''
        histblk = self.Histograms[TreeName]
        self.OpenFile(filename)
        self.Write(TreeName[5:])
        if prmname: self.Write('Instrument parameter file:'+os.path.split(prmname)[1])
        x = np.array(histblk['Data'][0])
        if 'T' in histblk['Instrument Parameters'][0]['Type'][0]:
            cw = np.diff(x)
            x[:-1] += cw
        else:
            x *= 100.
        # convert weights to sigmas; use largest weight as minimum esd
        s = np.sqrt(np.maximum(0.,np.array(histblk['Data'][2])))
        s[s==0] = np.max(s)
        s = 1./s
        self.Write('BANK 1 %d %d CONS %.2f %.2f 0 0 FXYE' % (
            len(x),len(x),x[0],(x[1]-x[0])
            ))
#            for X,Y,S in zip(x,histblk['Data'][1],s):
#                self.Write("{:15.6g} {:15.6g} {:15.6g}".format(X,Y,S))
        for XYS in zip(x,histblk['Data'][1],s):
            line = ''
            for val in XYS:
                line += G2fil.FormatPadValue(val,(15,6))
            self.Write(line)
        self.CloseFile()
        
    def Exporter(self,event=None):
        '''Export one or more sets of powder data as FXYE file(s)
        '''
        # the export process starts here
        self.InitExport(event)
        self.loadTree() # load all of the tree into a set of dicts
        if self.ExportSelect( # set export parameters
            AskFile='single' # get a file name/directory to save in
            ): return
        filenamelist = []
        for hist in self.histnam:
            if len(self.histnam) == 1:
                name = self.filename
            else:    # multiple files: create a unique name from the histogram
                name = self.MakePWDRfilename(hist)
            fileroot = os.path.splitext(G2obj.MakeUniqueLabel(name,filenamelist))[0]
            self.filename = os.path.join(self.dirname,fileroot + self.extension)
            histblk = self.Histograms[hist]
            # create an instrument parameter file
            prmname = self.WriteInstFile(hist,histblk['Instrument Parameters'][0])
            self.Writer(hist,prmname=prmname)
            print('Histogram '+hist+' written to '+self.fullpath)

class ExportPowderXYE(G2IO.ExportBaseclass):
    '''Used to create a Topas XYE file for a powder data set

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Topas XYE file',
            extension='.xye',
            longFormatName = 'Export powder data as Topas XYE (column) file'
            )
        self.exporttype = ['powder']
        self.multiple = True
        
    def Writer(self,TreeName,filename=None):
        self.OpenFile(filename)
        histblk = self.Histograms[TreeName]
        self.Write('/*')    #The ugly c comment delimiter used in topas!
        self.Write('# '+TreeName[5:])  #evidently this by itself fails in topas
        self.Write('*/')
        x = np.array(histblk['Data'][0])
        # convert weights to sigmas; use largest weight as minimum esd
        s = np.sqrt(np.maximum(0.,np.array(histblk['Data'][2])))
        s[s==0] = np.max(s)
        s = 1./s
        for XYS in zip(x,histblk['Data'][1],s):
            line = ''
            for val in XYS:
                line += G2fil.FormatPadValue(val,(15,6))
            self.Write(line)
        self.CloseFile()

    def Exporter(self,event=None):
        '''Export one or more sets of powder data as XYE file(s)
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile='single' # get a file name/directory to save in
            ): return
        filenamelist = []
        for hist in self.histnam:
            if len(self.histnam) == 1:
                name = self.filename
            else:    # multiple files: create a unique name from the histogram
                name = self.MakePWDRfilename(hist)
            fileroot = os.path.splitext(G2obj.MakeUniqueLabel(name,filenamelist))[0]
            self.filename = os.path.join(self.dirname,fileroot + self.extension)
            self.Writer(hist)
            print('Histogram '+hist+' written to '+self.fullpath)
