#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2export_pwdr: Export powder input files*
-------------------------------------------------

Creates files used by GSAS (FXYE) & TOPAS (XYE) as input

'''
import os.path
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
import GSASIIpy3 as G2py3
import GSASIIobj as G2obj

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
        self.Write(('INS  1PRCF11     %15.6g%15.6g%15.6g%15.6g   ')%(Inst['U'][1],Inst['V'][1],Inst['W'][1],0.0))
        self.Write(('INS  1PRCF12     %15.6g%15.6g%15.6g%15.6g   ')%(Inst['X'][1],Inst['Y'][1],Inst['SH/L'][1]/2.,Inst['SH/L'][1]/2.))
        self.CloseFile()
        print('Parameters from '+str(hist)+' written to file '+str(prmname))
        return prmname

    def Exporter(self,event=None):
        '''Export one or more sets of powder data as FXYE file(s)
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
            if len(self.histnam) > 1:
                # multiple files: create a unique name from the histogram
                fileroot = G2obj.MakeUniqueLabel(self.MakePWDRfilename(hist),filenamelist)
                # create an instrument parameter file
                self.filename = os.path.join(self.dirname,fileroot + self.extension)
            else:
                # use the supplied name, but force the extension
                self.filename= os.path.splitext(self.filename)[0] + self.extension
                
            histblk = self.Histograms[hist]
            prmname = self.WriteInstFile(hist,histblk['Instrument Parameters'][0])
            
            self.OpenFile()
            self.Write(hist[5:])
            self.Write('Instrument parameter file:'+os.path.split(prmname)[1])
            x = 100*np.array(histblk['Data'][0])
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
                    line += G2py3.FormatPadValue(val,(15,6))
                self.Write(line)
            self.CloseFile()
            print('Histogram '+str(hist)+' written to file '+str(self.fullpath))


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
            if len(self.histnam) > 1:
                # multiple files: create a unique name from the histogram
                fileroot = G2obj.MakeUniqueLabel(self.MakePWDRfilename(hist),filenamelist)
                # create an instrument parameter file
                self.filename = os.path.join(self.dirname,fileroot + self.extension)
            else:
                # use the supplied name, but force the extension
                self.filename= os.path.splitext(self.filename)[0] + self.extension

            self.OpenFile()
            histblk = self.Histograms[hist]
            self.Write('# '+hist[5:])
            x = np.array(histblk['Data'][0])
            # convert weights to sigmas; use largest weight as minimum esd
            s = np.sqrt(np.maximum(0.,np.array(histblk['Data'][2])))
            s[s==0] = np.max(s)
            s = 1./s
            for XYS in zip(x,histblk['Data'][1],s):
                line = ''
                for val in XYS:
                    line += G2py3.FormatPadValue(val,(15,6))
                self.Write(line)
            self.CloseFile()
            print('Histogram '+str(hist)+' written to file '+str(self.fullpath))
