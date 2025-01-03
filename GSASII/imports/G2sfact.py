# -*- coding: utf-8 -*-
'''Classes to read single crystal reflection files in formats used by:
Shelx, Jana, REMOS, TOPAS (SNS), HB-3A (HIFR)
'''
from __future__ import division, print_function
import sys
import numpy as np
from .. import GSASIIobj as G2obj

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


class HKLF_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F, sig(F) reflections from a HKLF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName = 'HKL F',
            longFormatName = 'Simple [hkl, Fo, sig(Fo)] Structure factor text file'
            )

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers'
        fp = open(filename,'r')
        valid = ColumnValidator(self, fp)
        fp.close()
        return valid

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                items = S.split()
                h,k,l,Fo,sigFo = items[:5]
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,1,0,Fo**2,2.*Fo*sigFo,0,Fo**2,0,0,0])
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = 0
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except:
            return False

class HKLMF_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F, reflections from a REMOS HKLMF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.fo','.FO'),
            strictExtension=False,
            formatName = 'HKLM F',
            longFormatName = 'REMOS [hklm, Fo] Structure factor text file'
            )

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers'
        fp = open(filename,'r')
        valid = ColumnValidator(self, fp)
        fp.close()
        return valid

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                h,k,l,m,Fo= S.split()
                h,k,l,m = [int(h),int(k),int(l),int(m)]
                if h == 999 or not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo2 = Fo
                if Fo < 1.0:
                    sigFo2 = 1.0
               # h,k,l,m,tw,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,m,1,0,Fo**2,sigFo2,0,Fo**2,0,0,0])
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = 1
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except:
            return False

class SHELX4_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a Shelx HKLF 4 file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = 'HKLF 4'
            longFormatName = 'Shelx HKLF 4 [hkl, Fo2, sig(Fo2)] Structure factor text file'
        else:
            formatName = u'Shelx HKLF 4 F\u00b2'
            longFormatName = u'Shelx HKLF 4 [hkl, Fo\u00b2, sig(Fo\u00b2)] Structure factor text file'
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers'
        return True
#        return ColumnValidator(self, filepointer)

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                try:   # use a simple split if possible
                    items = S.split()
                    h,k,l,Fo,sigFo = items[:5]
                    h,k,l = [int(h),int(k),int(l)]
                except: # but sometimes if no space between k & F
                        # need to use some fixed formatting
                    h,k,l = S[:12].split()
                    h,k,l = [int(h),int(k),int(l)]
                    Fo,sigFo = S[12:].split()[:2]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0])
                #self.RefDict['FF'].append({}) # now done in OnImportSfact
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = 0
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except:
            return False
            
class SHELX5_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F**2, sig(F**2) twin/incommensurate reflections from a fixed format SHELX HKLF5 file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = 'Shelx HKLF 5 F2 Tw/Incom'
            longFormatName = 'Shelx HKLF 5 [hklm, Fo2, sig(Fo2), Tind] Twin/incommensurate structure factor text file'
        else:
            formatName = u'Shelx HKLF 5 F\u00b2 Tw/Incom'
            longFormatName = u'Shelx HKLF 5 [hklm, Fo\u00b2, sig(Fo\u00b2), Tind] Twin/incommensurate structure factor text file'        
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)
        self.Super = 0

    def ContentsValidator(self, filename):
        '''Discover how many columns before F^2 are in the SHELX HKL5 file 
        - could be 3-6 depending on satellites'''
        numCols = 0
        fp = open(filename,'r')
        for i,line in enumerate(fp):
            for j,item in enumerate(line.split()):  #find 1st col with '.'; has F^2
                if '.' in item:
                    numCols = max(numCols,j)
                    break
            if i > 20:
                break
        fp.close()
        self.Super = numCols-3     #= 0,1,2,or 3
        if self.Super > 1:
            raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")            
        return True

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        TwDict = {}
        TwSet = {}
        TwMax = [-1,[]]
        first = True
        try:
            fp = open(filename,'r')
            m1 = 0
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if self.Super == 0:
                    h,k,l,Fo,sigFo,Tw = S.split()
    #                    h,k,l,Fo,sigFo,Tw = S[:4],S[4:8],S[8:12],S[12:20],S[20:28],S[28:32]
                    h,k,l = [int(h),int(k),int(l)]
                elif self.Super == 1:
                    h,k,l,m1,Fo,sigFo,Tw = S.split()
    #                    h,k,l,m1,Fo,sigFo,Tw = S[:4],S[4:8],S[8:12],S[12:16],S[16:24],S[24:32],S[32:36]
                    h,k,l,m1 = [int(h),int(k),int(l),int(m1)]
                Tw = Tw.strip()
                if Tw in ['','0']:
                    Tw = '1'
                if not any([h,k,l,m1]):
                    break
                if '-' in Tw:
                    if Tw == '-1':  #fix reversed twin ids
                        Tw = '-2'
                        if first:
                            self.warnings += '\nPrimary twin id changed to 1'
                            first = False
                    TwId = -int(Tw)-1
                    TwSet[TwId] = np.array([h,k,l])
                    if TwId not in TwMax[1]:
                        TwMax[1].append(TwId)
                else:
                    if Tw != '1':  #fix reversed twin ids
                        if first:
                            self.warnings += '\nPrimary twin id changed to 1\nNB: multiple primary twins not working'
                            first = False
                        Tw = '1'
                    TwId = int(Tw)-1
                    if TwSet:
                        TwDict[len(self.RefDict['RefList'])] = TwSet
                        TwSet = {}    
                    Fo = float(Fo)
                    sigFo = float(sigFo)
                    # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    if self.Super == 0:
                        self.RefDict['RefList'].append([h,k,l,int(Tw),0,Fo,sigFo,0,Fo,0,0,0])
                    elif self.Super == 1:
                        self.RefDict['RefList'].append([h,k,l,m1,int(Tw),0,Fo,sigFo,0,Fo,0,0,0])
                TwMax[0] = max(TwMax[0],TwId)
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = self.Super
            self.RefDict['TwDict'] = TwDict
            self.RefDict['TwMax'] = TwMax
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except:
            return False

class SHELX6_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F**2, sig(F**2) twin/incommensurate reflections from a fixed format SHELX HKLF6 file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = 'Shelx HKLF 6 F2 Tw/Incom'
            longFormatName = 'Shelx HKLF 6 [hklm, Fo2, sig(Fo2), Tind] Twin/incommensurate structure factor text file'
        else:
            formatName = u'Shelx HKLF 6 F\u00b2 Tw/Incom'
            longFormatName = u'Shelx HKLF 6 [hklm, Fo\u00b2, sig(Fo\u00b2), Tind] Twin/incommensurate structure factor text file'        
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hk6','.HK6'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)
        self.Super = 0

    def ContentsValidator(self, filename):
        '''Discover how many columns before F^2 are in the SHELX HKL6 file 
        - could be 3-6 depending on satellites'''
        numCols = 0
        fp = open(filename,'r')
        for i,line in enumerate(fp):
            for j,item in enumerate(line.split()):  #find 1st col with '.'; has F^2
                if '.' in item:
                    numCols = max(numCols,j)
                    break
            if i > 20:
                break
        fp.close()
        if numCols != 6:
            self.warnings += '\nInvalid hk6 file; wrong number of columns'
            #raise self.ImportException('Invalid hk6 file; wrong number of columns')
            return False
        self.Super = 1
        return True

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        TwDict = {}
        TwSet = {}
        TwMax = [-1,[]]
        first = True
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                h,k,l,m1,m2,m3,Fo,sigFo,Tw = S[:4],S[4:8],S[8:12],S[12:16],S[16:20],S[20:24],S[24:32],S[32:40],S[40:44]
                h,k,l,m1,m2,m3 = [int(h),int(k),int(l),int(m1),int(m2),int(m3)]
                if m2 or m3:
                    self.warnings += '\n(3+2) & (3+3) Supersymmetry not allowed in GSAS-II. Reformulate as twinned (3+1) supersymmetry'
                    raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")                                
                Tw = Tw.strip()
                if Tw in ['','0']:
                    Tw = '1'
                if not any([h,k,l]):    #look for 0 0 0 or blank line
                    break
                if '-' in Tw:
                    if Tw == '-1':  #fix reversed twin ids
                        Tw = '-2'
                        if first:
                            self.warnings += '\nPrimary twin id changed to 1'
                            first = False
                    TwId = -int(Tw)-1
                    TwSet[TwId] = np.array([h,k,l])
                    if TwId not in TwMax[1]:
                        TwMax[1].append(TwId)
                else:
                    if Tw != '1':  #fix reversed twin ids
                        if first:
                            self.warnings += '\nPrimary twin id changed to 1\nNB: multiple primary twins not working'
                            first = False
                        Tw = '1'
                    TwId = int(Tw)-1
                    if TwSet:
                        TwDict[len(self.RefDict['RefList'])] = TwSet
                        TwSet = {}    
                    Fo = float(Fo)
                    sigFo = float(sigFo)
                    # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    self.RefDict['RefList'].append([h,k,l,m1,int(Tw),0,Fo,sigFo,0,Fo,0,0,0])
                TwMax[0] = max(TwMax[0],TwId)
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = self.Super
            self.RefDict['TwDict'] = TwDict
            self.RefDict['TwMax'] = TwMax
            self.UpdateParameters(Type='SXC',Wave=None) # histogram type
            return True
        except:
            return False

class M90_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import F**2, sig(F**2) reflections from a JANA M90 file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            longFormatName = 'JANA [hkl, Fo2, sig(Fo2)] Structure factor text file'
        else:
            longFormatName = u'JANA [hkl, Fo\u00b2, sig(Fo\u00b2)] Structure factor text file'
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.m90','.M90'),
            strictExtension=False,
            formatName = u'JANA M90',
            longFormatName = longFormatName
            )
        self.Super = 0

    def ContentsValidator(self, filename):
        'Discover how many columns are in the m90 file - could be 9-12 depending on satellites'
        numCols = 0
        fp = open(filename,'r')
        for i,line in enumerate(fp):
            if 'Data' in line:
                startData = i
                break
        for i,line in enumerate(fp):
            if i > startData:
                numCols = max(numCols,len(line.split()))
            if i > startData+20:
                break
        fp.close()
        self.Super = numCols-9     #= 0,1,2,or 3
        if self.Super > 1:
            raise self.ImportException("Supersymmetry too high; GSAS-II limited to (3+1) supersymmetry")            
        return True #ColumnValidator(self, filepointer)

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                try:
                    if self.Super == 0:
                        h,k,l,Fo,sigFo = S.split()[:5]
                        h,k,l = [int(h),int(k),int(l)]
                    elif self.Super == 1:
                        h,k,l,m1,Fo,sigFo = S.split()[:6]
                        h,k,l,m1 = [int(h),int(k),int(l),int(m1)]
                except ValueError:  #skipping text at front
                    if not S:
                        break
                    text = S.split()
                    if not text:
                        break
                    print(line,text)
                    if text[0] == 'lambda':
                        wave = float(text[1])
                    continue
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                if self.Super == 0:
                    self.RefDict['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0])
                elif self.Super == 1:
                    self.RefDict['RefList'].append([h,k,l,m1,1,0,Fo,sigFo,0,Fo,0,0,0])
            fp.close()
            self.errors = 'Error after reading reflections (unexpected!)'
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SXC'
            self.RefDict['Super'] = self.Super
            self.UpdateParameters(Type='SXC',Wave=wave) # histogram type
            return True
        except:
            return False
            
class NT_HKLF2_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import neutron TOF F**2, sig(F**2) reflections from a HKLF file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = 'Neutron SNS TOF HKL F2'
            longFormatName = 'Neutron SNS TOF [hkl, Fo2, sig(Fo2),...] Structure factor text file'
        else:
            formatName = u'Neutron SNS TOF HKL F\u00b2'
            longFormatName = u'Neutron SNS TOF [hkl, Fo\u00b2, sig(Fo\u00b2),...] Structure factor text file'
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hkl','.HKL'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers & count number of data blocks - "Banks"'
        oldNo = -1
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                if not S:   #empty line terminates read
                    break
                if S[0] == '#': continue       #ignore comments, if any
                if S.split()[:3] == ['0','0','0']:
                    break
                try:
                    bankNo = S.split()[5]
                except:
                    return False
                if bankNo != oldNo:
                    self.Banks.append({'RefDict':{'RefList':[],}})
                    oldNo = bankNo
            fp.seek(0)
            valid = ColumnValidator(self, fp,nCol=8)
            fp.close()
            return valid
        except:
            return False

    def Reader(self,filename, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                if not S:
                    break
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
                    self.Banks[int(bN)-1]['RefDict']['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
                else:
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    self.RefDict['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
            fp.close()
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
        except:
            return False

class NT_JANA2K_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import neutron TOF F**2, sig(F**2) reflections from a JANA2000 file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = 'Neutron TOF JANA2000 F2'
            longFormatName = 'Neutron TOF [hkl, Fo2, sig(Fo2),...] Structure factor text file'
        else:
            formatName = u'Neutron TOF JANA2000 F\u00b2'
            longFormatName = u'Neutron TOF [hkl, Fo\u00b2, sig(Fo\u00b2),...] Structure factor text file'
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.int','.INT'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers & count number of data blocks - "Banks"'
        oldNo = -1
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                if not S:   #empty line terminates read
                    break
                if S[0] in ['#','(']: continue       #ignore comments & fortran format line
                bankNo = S.split()[5]
                if bankNo != oldNo:
                    self.Banks.append({'RefDict':{'RefList':[],}})
                    oldNo = bankNo
            fp.seek(0)
            valid = ColumnValidator(self, fp,nCol=10)
            fp.close()
            return valid
        except:
            return False

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
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
                    self.Banks[int(bN)-1]['RefDict']['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
                else:
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                    self.RefDict['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0,wave,tbar])
            fp.close()
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
        except:
            return False

class hb3a_INT_ReaderClass(G2obj.ImportStructFactor):
    'Routines to import neutron CW F**2, sig(F**2) reflections from a NIST hb3a int file'
    def __init__(self):
        if 'linux' in sys.platform:  # wx 3.0.0.0 on gtk does not like Unicode in menus
            formatName = u'Neutron HFIR HB-3A  CW HKL F2'
            longFormatName = u'Neutron HIFR HB-3A  CW HKL [hkl, Fo2, sig(Fo2),...] 5 column Structure factor text file'
        else:
            formatName = u'Neutron NIST HIFR  CW HKL F\u00b2'
            longFormatName = u'Neutron HFIR HB-3A  CW HKL [hkl, Fo\u00b2, sig(Fo\u00b2),...] 5 column Structure factor text file'
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.int','.INT'),
            strictExtension=False,
            formatName=formatName,
            longFormatName=longFormatName)

    def ContentsValidator(self, filename):
        'Make sure file contains the expected columns on numbers & count number of data blocks - "Banks"'
        fp = open(filename,'r')
        for line,S in enumerate(fp):
            if not S:   #empty line terminates read
                break
            if S[0] == '#': continue       #ignore comments, if any
        fp.seek(0)
        valid = ColumnValidator(self, fp,nCol=5)
        fp.close()
        return valid

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        'Read the file'
        try:
            fp = open(filename,'r')
            for line,S in enumerate(fp):
                self.errors = '  Error reading line '+str(line+1)
                if S[0] == '#': continue       #ignore comments, if any
                data = S.split()
                h,k,l,Fo,sigFo = data                   
                h,k,l = [int(h),int(k),int(l)]
                if not any([h,k,l]):
                    break
                Fo = float(Fo)
                sigFo = float(sigFo)
                # h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                self.RefDict['RefList'].append([h,k,l,1,0,Fo,sigFo,0,Fo,0,0,0])
            fp.close()
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.RefDict['Type'] = 'SNC'
            self.RefDict['Super'] = 0
            self.errors = 'Error after reading reflections (unexpected!)'
            self.UpdateParameters(Type='SNC',Wave=None) # histogram type
            return True
        except:
            return False
