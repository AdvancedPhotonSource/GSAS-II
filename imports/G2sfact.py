# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2sfact: simple HKL import*
-----------------------------------
Read structure factors from a simple hkl file. Two routines are
provided to read from files containing F or F\ :sup:`2` values.

'''
import sys
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

def ColumnValidator(parent, filepointer):
    'Validate a file to check that it contains columns of numbers'
    l = S = filepointer.readline()
    line = 1
    while '#' in S[0]:        #get past comments, if any
        l = filepointer.readline()        
        line += 1
    for i in range(10): # scan a few lines
        S = l.split()
        if len(S) != 5:
            parent.errors = 'line '+str(line)+': invalid input\n'+l
            return False
        for v in S:
            try:
                float(v)
            except ValueError:
                parent.errors = 'line '+str(line)+': string found where a number is expected\n'+l
                return False            
        l = filepointer.readline()
        line += 1
    return True


class HKLF_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F, sig(F) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = 'F containing HKL',
            longFormatName = 'Simple [hkl, Fo, sig(Fo)] Structure factor text file'
            )

    def ContentsValidator(self, filepointer):
        'Make sure file contains the expected columns on numbers'
        return ColumnValidator(self, filepointer)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                h,k,l,Fo,sigFo = S.split()
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,0,0,Fo**2,2.*Fo*sigFo,0,Fo**2,0,0,0])
                #self.RefDict['FF'].append({}) # now done in OnImportSfact
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.UpdateControls(Type='Fosq',FcalcPresent=False) # set Fobs type & if Fcalc values are loaded
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

class HKLF2_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = u'F\u00b2 containing HKL',
            longFormatName = u'Simple [hkl, Fo\u00b2, sig(Fo\u00b2)] Structure factor text file'
            )

    def ContentsValidator(self, filepointer):
        'Make sure file contains the expected columns on numbers'
        return ColumnValidator(self, filepointer)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                h,k,l,Fo,sigFo = S.split()
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0])
                #self.RefDict['FF'].append({}) # now done in OnImportSfact
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.UpdateControls(Type='Fosq',FcalcPresent=False) # set Fobs type & if Fcalc values are loaded
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

