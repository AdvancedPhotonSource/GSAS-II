# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-05-11 18:08:12 -0500 (Thu, 11 May 2023) $
# $Author: toby $
# $Revision: 5577 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2pdf_gr.py $
# $Id: G2pdf_gr.py 5577 2023-05-11 23:08:12Z toby $
########### SVN repository information ###################
'''
'''

from __future__ import division, print_function
import os.path as ospath
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5577 $")

class txt_FSQReaderClass(G2obj.ImportPDFData):
    'Routines to import S(Q) data from a .fq file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.fq','.sq'),
            strictExtension=False,
            formatName = 'q (A-1) step S(Q) data',
            longFormatName = 'q (A-1) stepped S(Q) PDF data from pdfGet or GSAS-II'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid r-step file'
        filepointer = open(filename,'r')
        Ndata = 0
        for i,S in enumerate(filepointer):
            if '#L' in S[:2]:
                break
            if len(S.split()) != 2:
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
                if '#L' in S[:2]:
                    ifData = True
                    continue
                if len(S) == 1:     #skip blank line
                    continue
                if len(S.split()) != 2:
                    continue
                self.comments.append(S[:-1])
            else:
                vals = S.split()
                if len(vals) >= 2:
                    try:
                        data = [float(val) for val in vals]
                        x.append(float(data[0]))
                        y.append(float(data[1]))
                    except ValueError:
                        msg = 'Error in line '+str(i+1)
                        print (msg,S)
                        continue
        self.Type = 'X f(q)'             #assume X-ray PDF
        self.pdfdata = np.array([
            np.array(x), # x-axis values q
            np.array(y), # pdf f(q))
            ])
        self.pdfentry[0] = filename
        self.pdfentry[2] = 1 # xy file only has one bank
        self.idstring = ospath.basename(filename)
        filepointer.close()
        return True

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
            if '#L r' in S[:4]:
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
        print ('Read a r-step text file')
        x = []
        y = []
        ifData = False
        filepointer = open(filename,'r')
        self.Type = 'N g(r)'
        for i,S in enumerate(filepointer):
            if not ifData:
                if len(S) == 1:     #skip blank line
                    continue
                if '#L r' in S[:4]:
                    ifData = True
                if 'X-Ray' in S:
                    self.Type = 'X g(r)'
                self.comments.append(S[:-1])
            else:
                vals = S.split()
                if len(vals) >= 2:
                    try:
                        data = [float(val) for val in vals]
                        x.append(float(data[0]))
                        y.append(float(data[1]))
                    except ValueError:
                        msg = 'Error in line '+str(i+1)
                        print (msg,S[:-1])
                        continue
        self.pdfdata = np.array([
            np.array(x), # x-axis values r
            np.array(y), # pdf g(r)
            ])
        self.pdfentry[0] = filename
        self.pdfentry[2] = 1 # xy file only has one bank
        self.idstring = ospath.basename(filename)
        filepointer.close()
        return True

class txt_PDFReaderClassG(G2obj.ImportPDFData):
    'Routines to import PDF G(R) data from a .dat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.dat',),
            strictExtension=False,
            formatName = 'gudrun r (A) step G(r) data',
            longFormatName = 'r (A) stepped G(r) PDF data from gudrun'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid r-step file'
        filepointer = open(filename,'r')
        Ndata = 0
        for i,S in enumerate(filepointer):
            if i < 2:
                continue
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
        print ('Read a r-step text file')
        x = []
        y = []
        filepointer = open(filename,'r')
        for i,S in enumerate(filepointer):
            if i < 2:
                continue
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
            else:
                break
        self.pdfdata = np.array([
            np.array(x), # x-axis values r
            np.array(y), # pdf g(r)
            ])
        self.Type = 'X g(r)'         #assume X-ray PDF
        self.pdfentry[0] = filename
        self.pdfentry[2] = 1 # xy file only has one bank
        self.idstring = ospath.basename(filename)
        filepointer.close()
        return True

