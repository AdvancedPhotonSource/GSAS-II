# -*- coding: utf-8 -*-
'''Classes in :mod:`~GSASII.exports.G2export_map` follow:
'''
from __future__ import division, print_function
import platform
import os
import numpy as np
from .. import GSASIIfiles as G2fil

class ExportMapASCII(G2fil.ExportBaseclass):
    '''Used to create a text file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'FOX/DrawXTL file',
            extension='.grd',
            longFormatName = 'Export map as text (.grd) file'
            )
        self.exporttype = ['map']
        self.multiple = False 

    def Exporter(self,event=None):
        '''Export a map as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect(): return  # set export parameters, get file name
        filename = self.filename
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            mapData = phasedict['General']['Map']
            rho = phasedict['General']['Map'].get('rho',[])
            if not len(rho):
                print ("There is no map for phase "+phasenam)
                continue
            if len(self.phasenam) > 1: # if more than one filename is written, add a phase # -- not in use yet
                i = self.Phases[phasenam]['pId']
                self.filename = os.path.splitext(filename)[1] + "_" + mapData['MapType'] + str(i) + self.extension
            self.OpenFile()
            self.Write(u"Map of Phase "+phasenam+u" from "+self.G2frame.GSASprojectfile)
            # get cell parameters & print them
            cellList,cellSig = self.GetCell(phasenam)
            fmt = 3*" {:9.5f}" + 3*" {:9.3f}"
            self.Write(fmt.format(*cellList[:6]))
            nx,ny,nz = rho.shape
            self.Write(" {:3d} {:3d} {:3d}".format(nx,ny,nz))
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        self.Write(str(rho[i,j,k]))
            self.CloseFile()
            print(u'map from Phase '+phasenam+u' written to file '+self.fullpath)

class ExportMapCCP4(G2fil.ExportBaseclass):
    '''Used to create a text file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'CCP4 map file',
            extension='.map',
            longFormatName = 'Export CCP4 .map file'
            )
        self.exporttype = ['map']
        self.multiple = False 

    # Tools for file writing. 
    def Write(self,data,dtype):
        import struct
        '''write a block of output

        :param data: the data to be written. 
        '''
        self.fp.write(struct.pack(dtype,data))
        
    def Exporter(self,event=None):
        '''Export a map as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect(): return  # set export parameters, get file name
        filename = self.filename
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info 
            mapData = phasedict['General']['Map']
            rho = mapData.get('rho',[])
            
            if not len(rho):
                print ("There is no map for phase "+phasenam)
                continue
            if len(self.phasenam) > 1: # if more than one filename is written, add a phase # -- not in use yet
                i = self.Phases[phasenam]['pId']
                self.filename = os.path.splitext(filename)[1] + "_" + mapData['MapType'] + str(i) + self.extension
            cell = phasedict['General']['Cell'][1:7]
            nx,ny,nz = rho.shape
            self.OpenFile(mode='wb')
            for n in rho.shape: self.Write(n,'i')  #nX,nY,nZ
            self.Write(2,'i')           #mode=2 float map
            for i in [0,0,0]: self.Write(i,'i')    #1st position on x,y,z
            for n in rho.shape: self.Write(n,'i')  #nX,nY,nZ
            for c in cell: self.Write(c,'f')
            for i in [1,2,3]: self.Write(i,'i')    #axes order = x,y,z
            self.Write(np.min(rho),'f')
            self.Write(np.max(rho),'f')
            self.Write(np.mean(rho),'f')
            self.Write(0,'i')
            for i in range(24,53):
                self.Write(0,'i')
            if '2' in platform.python_version_tuple()[0]:
                for s in ['M','A','P',' ']: self.Write(s,'s')
            else:
                for s in [b'M',b'A',b'P',b' ']: self.Write(s,'c')
            self.Write(0x44410000,'i')
            self.Write(np.std(rho),'f')
            for i in range(56,257):
                self.Write(0,'i')
            for x in rho.flatten('F'):
                self.Write(x,'f')
            self.CloseFile()
            print(u'map from Phase '+phasenam+u' written to file '+self.fullpath)
