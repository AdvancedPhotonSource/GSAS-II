# -*- coding: utf-8 -*-
'''Classes in :mod:`G2export_image` follow:
'''
from __future__ import division, print_function
import os.path
import scipy.misc
import GSASIIfiles as G2fil

class ExportImagePNG(G2fil.ExportBaseclass):
    '''Used to create a PNG file for a GSAS-II image

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'PNG image file',
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
        if self.ExportSelect(): return # select one image; ask for a file name
        # process the selected image(s) (at present only one image)
        for i in sorted(self.histnam): 
            filename = os.path.join(
                self.dirname,
                os.path.splitext(self.filename)[0] + self.extension
                )
            imgFile = self.Histograms[i].get('Data',(None,None))
            Image = G2fil.GetImageData(self.G2frame,imgFile,imageOnly=True)
            scipy.misc.imsave(filename,Image)
            print('Image '+imgFile+' written to file '+filename)
            
