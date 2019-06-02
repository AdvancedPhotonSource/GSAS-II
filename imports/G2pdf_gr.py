# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2017-03-04 08:34:05 -0600 (Sat, 04 Mar 2017) $
# $Author: vondreele $
# $Revision: 2738 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2sad_xye.py $
# $Id: G2sad_xye.py 2738 2017-03-04 14:34:05Z vondreele $
########### SVN repository information ###################
'''
*Module G2pdf_gr: read PDF G(R) data*
------------------------------------------------

Routines to read in G(R) data from an .gr type file, with
Angstrom steps. 

'''

from __future__ import division, print_function
import os.path as ospath
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2738 $")

class txt_PDFReaderClass(G2obj.ImportPDFData):
    'Routines to import PDF G(R) data from a .gr file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.gr',),
            strictExtension=False,
            formatName = 'r (A) step G(r) data',
            longFormatName = 'r (A) stepped G(r) PDF data from pdfGet or GSAS-II'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid r-step file'
        filepointer = open(filename,'r')
        Ndata = 0
        for i,S in enumerate(filepointer):
            if '#L' in S[:2]:
                break
        for i,S in enumerate(filepointer):            
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    Ndata += 1
                except ValueError:
                    pass
        if not Ndata:     
            self.errors = 'No 2 or more column numeric data found'
            filepointer.close()
            return False
        filepointer.close()
        return True # no errors encountered

    def Reader(self,filename,ParentFrame=None, **unused):
        print ('Read a q-step text file')
        x = []
        y = []
        ifData = False
        filepointer = open(filename,'r')
        for i,S in enumerate(filepointer):
            if not ifData:
                if len(S) == 1:     #skip blank line
                    continue
                self.comments.append(S[:-1])
                if '#L' in S[:2]:
                    ifData = True
            else:
                vals = S.split()
                if len(vals) >= 2:
                    try:
                        data = [float(val) for val in vals]
                        x.append(float(data[0]))
                        y.append(float(data[1]))
                    except ValueError:
                        msg = 'Error in line '+str(i+1)
                        print (msg)
                        continue
        self.pdfdata = np.array([
            np.array(x), # x-axis values r
            np.array(y), # pdf g(r)
            ])
        self.pdfentry[0] = filename
        self.pdfentry[2] = 1 # xy file only has one bank
        self.idstring = ospath.basename(filename)

        return True

