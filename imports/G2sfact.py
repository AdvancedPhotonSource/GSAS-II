# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
# short routines to read in structure factors from simple file formats
# 
import sys
import numpy as np
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

class HKLF_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F, sig(F) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = 'F HKL',
            longFormatName = 'Simple (hkl Fo sig(Fo)) Structure factor text file'
            )
    # Validate the contents
    def ContentsValidator(self, filepointer):
        S = filepointer.readline() 
        while '#' in S[0]:        #get past comments, if any
            S = filepointer.readline()        
        for i in range(3): # scan a few lines
            S = S.split()
            if len(S) != 5: return False
            for v in S:
                try:
                    float(v)
                except ValueError:
                    return False            
            S = filepointer.readline()
        return True

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        try:
            for S in filepointer:
                if S[0] == '#': continue       #ignore comments, if any
                h,k,l,Fo,sigFo = S.split()
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefList.append([h,k,l,0,0,Fo**2,2.*Fo*sigFo,0,Fo**2,0,0,[],[],0,{}])
            self.UpdateControls(Type='Fosq',FcalcPresent=False) # set Fobs type & if Fcalc values are loaded
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

class HKLF2_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = 'F**2 HKL',
            longFormatName = 'Simple (hkl Fo**2 sig(Fo**2)) Structure factor text file'
            )
    # Validate the contents
    def ContentsValidator(self, filepointer):
        S = filepointer.readline() 
        while '#' in S[0]:        #get past comments, if any
            S = filepointer.readline()        
        for i in range(3): # scan a few lines
            S = S.split()
            if len(S) != 5: return False
            for v in S:
                try:
                    float(v)
                except ValueError:
                    return False            
            S = filepointer.readline()
        return True

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        try:
            for S in filepointer:
                if S[0] == '#': continue       #ignore comments, if any
                h,k,l,Fo,sigFo = S.split()
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefList.append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,[],[],0,{}])
            self.UpdateControls(Type='Fosq',FcalcPresent=False) # set Fobs type & if Fcalc values are loaded
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False
