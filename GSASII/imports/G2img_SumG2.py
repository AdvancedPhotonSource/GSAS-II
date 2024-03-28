# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2img_SumG2.py $
# $Id: G2img_SumG2.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import platform
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5577 $")
class G2_ReaderClass(G2obj.ImportImage):
    '''Routine to read an image that has been pickled in Python. Images
    in this format are created by the "Sum image data" command. At least for
    now, only one image is permitted per file.
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.G2img',),
            strictExtension=True,
            formatName = 'GSAS-II image',
            longFormatName = 'cPickled image from GSAS-II'
            )

    def ContentsValidator(self, filename):
        '''test by trying to unpickle (should be quick)
        '''
        try:
            fp = open(filename,'rb')
            cPickle.load(fp)
            fp.close()
        except:
            return False
        return True
        
    def Reader(self,filename, ParentFrame=None, **unused):
        '''Read using cPickle
        '''
        Fp = open(filename,'rb')
        self.Comments,self.Data,self.Npix,self.Image = cPickle.load(Fp)
        Fp.close()
        self.LoadImage(ParentFrame,filename)
        return True
