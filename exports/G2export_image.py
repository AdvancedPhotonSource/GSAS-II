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
*Module G2export_image: 2D Image data export*
------------------------------------------------------

Demonstrates how an image is retrieved and written. Uses
a SciPy routine to write a PNG format file. 
'''
import os.path
import scipy.misc
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
import GSASIImath as G2mth

class ExportImagePNG(G2IO.ExportBaseclass):
    '''Used to create a PNG file for a GSAS-II image

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'PNG',
            extension='.png',
            longFormatName = 'Export image in PNG format'
            )
        self.exporttype = ['image']
        #self.multiple = True
    def Exporter(self,event=None):
        '''Export an image
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect(False): return
        # process the selected image(s) (at present only one image)
        for i in sorted(self.histnam): 
            imgFile = self.Histograms[i].get('Data',(None,None))
            Comments,Data,Npix,Image = G2IO.GetImageData(self.G2frame,imgFile)
            scipy.misc.imsave(self.filename,Image)
            print('Image '+str(imgFile)+' written to file '+str(self.filename))                   
