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
import copy
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

def ColumnValidator(parent, filepointer,nCol=5):
    'Validate a file to check that it contains columns of numbers'
    l = filepointer.readline()
    line = 1
    while l[0] in ['#','(']:        #get past comments & fortran formats, if any
        l = filepointer.readline()        
        line += 1
    for i in range(10): # scan a few lines
        S = l.split()
        if len(S) < nCol:
            parent.errors = 'line '+str(line)+': invalid input\n'+l
            return False
        for v in S[:nCol]:
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
            formatName = 'HKL F',
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
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = 0
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
            formatName = u'HKL F\u00b2',
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
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = 0
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False
            
class M90_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a JANA M90 file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.m90','.M90'),
            strictExtension=False,
            formatName = u'JANA M90',
            longFormatName = u'Simple [hkl, Fo\u00b2, sig(Fo\u00b2)] Structure factor text file'
            )
        self.Super = 0

    def ContentsValidator(self, filepointer):
        'Discover how many columns are in the m90 file - could be 9-12 depending on satellites'
        numCols = 0
        for i,line in enumerate(filepointer):
            if 'Data' in line:
                startData = i
                break
        for i,line in enumerate(filepointer):
            if i > startData:
                numCols = max(numCols,len(line.split()))
            if i > startData+20:
                break
        self.Super = numCols-9     #= 0,1,2,or 3
        return True #ColumnValidator(self, filepointer)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                try:
                    if self.Super == 0:
                        h,k,l,Fo,sigFo = S.split()[:5]
                        h,k,l = [int(h),int(k),int(l)]
                    elif self.Super == 1:
                        h,k,l,m1,Fo,sigFo = S.split()[:6]
                        h,k,l,m1 = [int(h),int(k),int(l),int(m1)]
                    elif self.Super == 2:
                        h,k,l,m1,m2,Fo,sigFo = S.split()[:7]
                        h,k,l,m1,m2 = [int(h),int(k),int(l),int(m1),int(m2)]
                    elif self.Super == 3:
                        h,k,l,m1,m2,m3,Fo,sigFo = S.split()[:8]
                        h,k,l,m1,m2,m3 = [int(h),int(k),int(l),int(m1),int(m2),int(m3)]                        
                except ValueError:  #skipping text at front
                    text = S.split()
                    if text[0] == 'lambda':
                        wave = float(text[1])
                    continue
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                if self.Super == 0:
                    self.RefDict['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 1:
                    self.RefDict['RefList'].append([h,k,l,m1,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 2:
                    self.RefDict['RefList'].append([h,k,l,m1,m2,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 3:
                    self.RefDict['RefList'].append([h,k,l,m1,m2,m3,0,0,Fo,sigFo,0,Fo,0,0,0])
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = self.Super
            self.UpdateParameters(Type='SXC',Wave=wave) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False
            
class SHELX5_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a fixed format SHELX HKLF5 file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = u'SHELX HKL F\u00b2',
            longFormatName = u'SHELX HKLF5 [hkl, Fo\u00b2, sig(Fo\u00b2)] Structure factor text file'
            )
        self.Super = 0

    def ContentsValidator(self, filepointer):
        'Discover how many characters are in the SHELX file - could be 32-44 depending on satellites'
        numCols = 0
        for i,line in enumerate(filepointer):
            numCols = max(numCols,len(line))
            if i > 20:
                break
        self.Super = (numCols-33)/4     #= 0,1,2,or 3
        print numCols,self.Super
        return True #ColumnValidator(self, filepointer)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if self.Super == 0:
                    h,k,l,Fo,sigFo = S[:4],S[4:8],S[8:12],S[12:20],S[20:28]
                    h,k,l = [int(h),int(k),int(l)]
                elif self.Super == 1:
                    h,k,l,m1,Fo,sigFo = S[:4],S[4:8],S[8:12],S[12:16],S[16:24],S[24:32]
                    h,k,l,m1 = [int(h),int(k),int(l),int(m1)]
                elif self.Super == 2:
                    h,k,l,m1,m2,Fo,sigFo = S[:4],S[4:8],S[8:12],S[12:16],S[16:20],S[20:28],S[28:36]
                    h,k,l,m1,m2 = [int(h),int(k),int(l),int(m1),int(m2)]
                elif self.Super == 3:
                    h,k,l,m1,m2,m3,Fo,sigFo = S[:4],S[4:8],S[8:12],S[12:16],S[16:20],S[20:24],S[24:32],S[32:40]
                    h,k,l,m1,m2,m3 = [int(h),int(k),int(l),int(m1),int(m2),int(m3)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                if self.Super == 0:
                    self.RefDict['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 1:
                    self.RefDict['RefList'].append([h,k,l,m1,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 2:
                    self.RefDict['RefList'].append([h,k,l,m1,m2,0,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 3:
                    self.RefDict['RefList'].append([h,k,l,m1,m2,m3,0,0,Fo,sigFo,0,Fo,0,0,0])
                #self.RefDict['FF'].append({}) # now done in OnImportSfact
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = self.Super
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

class NT_HKLF2_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import neutron TOF F**2, sig(F**2) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = u'Neutron TOF HKL F\u00b2',
            longFormatName = u'Neutron TOF [hkl, Fo\u00b2, sig(Fo\u00b2),...] Structure factor text file'
            )

    def ContentsValidator(self, filepointer):
        'Make sure file contains the expected columns on numbers & count number of data blocks - "Banks"'
        oldNo = -1
        for line,S in enumerate(filepointer):
            if not S:   #empty line terminates read
                break
            if S[0] == '#': continue       #ignore comments, if any
            if S.split()[:3] == ['0','0','0']:
                break
            bankNo = S.split()[5]
            if bankNo != oldNo:
                self.Banks.append({'RefDict':{'RefList':[],}})
                oldNo = bankNo
        filepointer.seek(0)
        return ColumnValidator(self, filepointer,nCol=8)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        filepointer.seek(0)
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                data = S.split()
                h,k,l,Fo,sigFo,bN,wave,tbar = data[:8]  #bN = 1..., 6 dir cos next                    
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                wave = float(wave)
                tbar = float(tbar)
                if len(self.Banks):
                    self.Banks[int(bN)-1]['RefDict']['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
                else:
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    self.RefDict['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
            if len(self.Banks):
                self.UpdateParameters(Type='SNT',Wave=None) # histogram type
                for Bank in self.Banks:
                    Bank['RefDict']['RefList'] = np.array(Bank['RefDict']['RefList'])
                    Bank['RefDict']['Type'] = 'SNT'                    
                    Bank['RefDict']['Super'] = 0
            else:
                self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
                self.RefDict['Type'] = 'SNT'
                self.RefDict['Super'] = 0
                self.errors = 'Error after reading reflections (unexpected!)'
                self.UpdateParameters(Type='SNT',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

class NT_JANA2K_ReaderClass(G2IO.ImportStructFactor):
    'Routines to import neutron TOF F**2, sig(F**2) reflections from a JANA2000 file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.int','.INT'),
            strictExtension=False,
            formatName = u'Neutron TOF JANA2000 F\u00b2',
            longFormatName = u'Neutron TOF [hkl, Fo\u00b2, sig(Fo\u00b2),...] Structure factor text file'
            )

    def ContentsValidator(self, filepointer):
        'Make sure file contains the expected columns on numbers & count number of data blocks - "Banks"'
        oldNo = -1
        for line,S in enumerate(filepointer):
            if not S:   #empty line terminates read
                break
            if S[0] in ['#','(']: continue       #ignore comments & fortran format line
            bankNo = S.split()[5]
            if bankNo != oldNo:
                self.Banks.append({'RefDict':{'RefList':[],}})
                oldNo = bankNo
        filepointer.seek(0)
        return ColumnValidator(self, filepointer,nCol=10)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        filepointer.seek(0)
        try:
            for line,S in enumerate(filepointer):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] in ['#','(']: continue       #ignore comments & fortran format line
                data = S.split()
                h,k,l,Fo,sigFo,bN,wave,x,x,tbar = data[:10]  #bN = 1..., 6 dir cos next                    
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                wave = float(wave)
                tbar = float(tbar)
                if len(self.Banks):
                    self.Banks[int(bN)-1]['RefDict']['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
                else:
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    self.RefDict['RefList'].append([h,k,l,0,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
            if len(self.Banks):
                self.UpdateParameters(Type='SNT',Wave=None) # histogram type
                for Bank in self.Banks:
                    Bank['RefDict']['RefList'] = np.array(Bank['RefDict']['RefList'])
                    Bank['RefDict']['Type'] = 'SNT'                    
                    Bank['RefDict']['Super'] = 0        #for now                    
            else:
                self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
                self.RefDict['Type'] = 'SNT'
                self.RefDict['Super'] = 0   #for now
                self.errors = 'Error after reading reflections (unexpected!)'
                self.UpdateParameters(Type='SNT',Wave=None) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print '\n\n'+self.formatName+' read error: '+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
            return False

