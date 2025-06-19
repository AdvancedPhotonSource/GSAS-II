# -*- coding: utf-8 -*-
'''Classes in :mod:`~GSASII.exports.G2export_examples` follow:
'''
# note documentation in docs/source/exports.rst
#
from __future__ import division, print_function
import os
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIImath as G2mth
from .. import GSASIIpwd as G2pwd
from .. import GSASIIfiles as G2fil

class ExportPhaseText(G2fil.ExportBaseclass):
    '''Used to create a text file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Text file',
            extension='.txt',
            longFormatName = 'Export phase as text file'
            )
        self.exporttype = ['phase']
        self.multiple = True # allow multiple phases to be selected

    def Exporter(self,event=None):
        '''Export a phase as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.ExportSelect(): return # set export parameters; prompt for file name
        self.OpenFile()
        # if more than one format is selected, put them into a single file
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            i = self.Phases[phasenam]['pId']
            self.Write('\n'+80*'=')
            self.Write("Phase "+str(phasenam)+" from "+str(self.G2frame.GSASprojectfile))
            self.Write("\nSpace group = "+str(phasedict['General']['SGData']['SpGrp'].strip()))
            # get cell parameters & print them
            cellList,cellSig = self.GetCell(phasenam)
            prevsig = 0
            for lbl,defsig,val,sig in zip(
                ['a','b','c','alpha','beta ','gamma','volume'],
                3*[-0.00001] + 3*[-0.001] + [-0.01], # sign values to use when no sigma
                cellList,cellSig
                ):
                if sig:
                    txt = G2mth.ValEsd(val,sig)
                    prevsig = -sig # use this as the significance for next value
                else:
                    txt = G2mth.ValEsd(val,min(defsig,prevsig),True)
                self.Write(lbl+' = '+txt)
            # get atoms and print them in nice columns
            AtomsList = self.GetAtoms(phasenam)
            fmt = "{:8s} {:4s} {:4s} {:12s} {:12s} {:12s} {:10s} {:10s}"
            self.Write('\nAtoms\n'+80*'-')
            self.Write(fmt.format("label","elem","mult","x","y","z","frac","Uiso"))
            self.Write(80*'-')
            aniso = False
            for lbl,typ,mult,xyz,td in AtomsList:
                vals = [lbl,typ,str(mult)]
                if xyz[3][0] == 0: continue
                for val,sig in xyz:
                    vals.append(G2mth.ValEsd(val,sig))
                if len(td) == 1:
                    vals.append(G2mth.ValEsd(td[0][0],td[0][1]))
                else:
                    vals.append("aniso")
                    aniso = True
                self.Write(fmt.format(*vals))
            # print anisotropic values, if any
            if aniso:
                self.Write('\nAnisotropic parameters')
                self.Write(80*'-')
                fmt = "{:8s} {:4s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s}"
                self.Write(fmt.format("label","elem",'U11','U22','U33','U12','U13','U23'))
                self.Write(80*'-')
                for lbl,typ,mult,xyz,td in AtomsList:
                    if len(td) == 1: continue
                    if xyz[3][0] == 0: continue
                    vals = [lbl,typ]
                    for val,sig in td:
                        vals.append(G2mth.ValEsd(val,sig))
                    self.Write(fmt.format(*vals))
            print('Phase '+phasenam+' written to file '+self.fullpath)
        self.CloseFile()

class ExportPowderText(G2fil.ExportBaseclass):
    '''Used to create a text file for a powder data set

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Text file',
            extension='.txt',
            longFormatName = 'Export powder data as text file'
            )
        self.exporttype = ['powder']
        self.multiple = True # allow one or more histogram(s) to be selected

    def Writer(self,TreeName,filename=None):
        self.OpenFile(filename)
        histblk = self.Histograms[TreeName]
        hfmt = 5*"{:12s} "
        digitList = 2*((13,3),) + ((13,5),) + 2*((13,3),)
        
        self.Write(hfmt.format("x","y_obs","weight","y_calc","y_bkg"))
        for vallist in zip(histblk['Data'][0],
                           histblk['Data'][1],
                           histblk['Data'][2],
                           histblk['Data'][3],
                           histblk['Data'][4],
                           #histblk['Data'][5],
                           ):
            strg = ''
            for val,digits in zip(vallist,digitList):
                strg += G2fil.FormatPadValue(val,digits)
            self.Write(strg)
        self.CloseFile()
        
    def Exporter(self,event=None):
        '''Export a set of powder data as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile='single' # selected one or more histograms; get file name (1 hist) or a directory (>1)
            ): return 
        filenamelist = []
        for hist in self.histnam:
            if len(self.histnam) == 1:
                name = self.filename
            else:    # multiple files: create a unique name from the histogram
                name = self.MakePWDRfilename(hist)
            fileroot = os.path.splitext(G2obj.MakeUniqueLabel(name,filenamelist))[0]
            self.filename = os.path.join(self.dirname,fileroot + self.extension)
            self.Writer(hist)
            print(hist+' written to file '+self.fullpath)
        
class ExportPowderReflText(G2fil.ExportBaseclass):
    '''Used to create a text file of reflections from a powder data set

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'reflection list as text',
            extension='.txt',
            longFormatName = 'Export powder reflection list as a text file'
            )
        self.exporttype = ['powder']
        self.multiple = False # only allow one histogram to be selected

    def Exporter(self,event=None):
        '''Export a set of powder reflections as a text file
        '''
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile='default' # base name on the GPX file name
            ): return 
        self.OpenFile()
        hist = list(self.histnam)[0] # there should only be one histogram, in any case take the 1st
        self.Write('\nHistogram '+hist)
        histblk = self.Histograms[hist]
        for phasenam in histblk['Reflection Lists']:
            phasDict = histblk['Reflection Lists'][phasenam]
            tname = {'T':'TOF','C':'2-theta','B':'2-theta'}[phasDict['Type'][2]]
            self.Write('\nPhase '+str(phasenam))
            if phasDict.get('Super',False):
                self.Write(96*'=')
                hklfmt = "{:.0f},{:.0f},{:.0f},{:.0f}"
                hfmt = "{:>10s} {:>8s} {:>12s} {:>12s} {:>7s} {:>6s} {:>8s} {:>8s} {:>8s} {:>8s}"
                if 'T' in phasDict['Type']:
                    fmt = "{:>10s} {:8.3f} {:12.3f} {:12.3f} {:7.2f} {:6.0f} {:8.3f} {:8.3f} {:8.3f} {:8.4f}"
                else:
                    fmt = "{:>10s} {:8.3f} {:12.3f} {:12.3f} {:7.2f} {:6.0f} {:8.5f} {:8.5f} {:8.5f} {:8.4f}"
                self.Write(hfmt.format("h,k,l,m",tname,"F_obs","F_calc","phase","mult","sig","gam","FWHM","Prfo"))
                self.Write(96*'=')
                refList = phasDict['RefList']
                for refItem in refList:
                    if 'T' in phasDict['Type']:
                        h,k,l,m,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,x,x,x,prfo = refItem[:17]
                        FWHM = G2pwd.getgamFW(gam,sig)
                        self.Write(fmt.format(hklfmt.format(h,k,l,m),pos,Fobs,Fcalc,phase,mult,sig,gam,FWHM,prfo))
                    elif 'C' in phasDict['Type']:
                        h,k,l,m,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,prfo = refItem[:14]
                        g = gam/100.    #centideg -> deg
                        s = np.sqrt(max(sig,0.0001))/100.   #var -> sig in deg
                        FWHM = G2pwd.getgamFW(g,s)  #FWHM
                        self.Write(fmt.format(hklfmt.format(h,k,l,m),pos,Fobs,Fcalc,phase,mult,s,g,FWHM,prfo))
                    elif 'B' in phasDict['Type']:
                        h,k,l,m,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,x,x,x,prfo = refItem[:17]
                        g = gam/100.    #centideg -> deg
                        s = np.sqrt(max(sig,0.0001))/100.   #var -> sig in deg
                        FWHM = G2pwd.getgamFW(g,s)  #FWHM
                        self.Write(fmt.format(hklfmt.format(h,k,l,m),pos,Fobs,Fcalc,phase,mult,s,g,FWHM,prfo))
            else:
                self.Write(94*'=')
                hklfmt = "{:.0f},{:.0f},{:.0f}"
                hfmt = "{:>8s} {:>8s} {:>12s} {:>12s} {:>7s} {:>6s} {:>8s} {:>8s} {:>8s} {:>8s}"
                if 'T' in phasDict['Type']:
                    fmt = "{:>8s} {:8.3f} {:12.3f} {:12.3f} {:7.2f} {:6.0f} {:8.3f} {:8.3f} {:8.3f} {:8.4f}"
                else:
                    fmt = "{:>8s} {:8.3f} {:12.3f} {:12.3f} {:7.2f} {:6.0f} {:8.5f} {:8.5f} {:8.5f} {:8.4f}"
                self.Write(hfmt.format("h,k,l",tname,"F_obs","F_calc","phase","mult","sig","gam","FWHM","Prfo"))
                self.Write(94*'=')
                refList = phasDict['RefList']
                for refItem in refList:
                    if 'T' in phasDict['Type']:
                        h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,x,x,x,prfo = refItem[:16]
                        FWHM = G2pwd.getgamFW(gam,sig)
                        self.Write(fmt.format(hklfmt.format(h,k,l),pos,Fobs,Fcalc,phase,mult,np.sqrt(max(sig,0.0001)),gam,FWHM,prfo))
                    elif 'C' in phasDict['Type']:
                        h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,prfo = refItem[:13]
                        g = gam/100.    #centideg -> deg
                        s = np.sqrt(max(sig,0.0001))/100.   #var -> sig in deg
                        FWHM = G2pwd.getgamFW(g,s)
                        self.Write(fmt.format(hklfmt.format(h,k,l),pos,Fobs,Fcalc,phase,mult,   \
                            s,g,FWHM,prfo))
                    elif 'B' in phasDict['Type']:
                        h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,x,x,x,x,prfo = refItem[:16]
                        g = gam/100.    #centideg -> deg
                        s = np.sqrt(max(sig,0.0001))/100.   #var -> sig in deg
                        FWHM = G2pwd.getgamFW(g,s)
                        self.Write(fmt.format(hklfmt.format(h,k,l),pos,Fobs,Fcalc,phase,mult,   \
                            s,g,FWHM,prfo))
        self.CloseFile()
        print(hist+'reflections written to file '+self.fullpath)                        

class ExportSingleText(G2fil.ExportBaseclass):
    '''Used to create a text file with single crystal reflection data
    skips user rejected & space group extinct reflections

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Text file',
            extension='.txt',
            longFormatName = 'Export reflection list as a text file'
            )
        self.exporttype = ['single']
        self.multiple = False # only allow one histogram to be selected

    def Writer(self,hist,filename=None):
        self.OpenFile(filename)
        hklfmt = "{:.0f},{:.0f},{:.0f}"
        hfmt = "{:>10s} {:>8s} {:>12s} {:>12s} {:>12s} {:>7s} {:>6s}"
        fmt = "{:>10s} {:8.3f} {:12.2f} {:12.4f} {:12.2f} {:7.2f} {:6.0f}"
        self.Write(80*'=')
        self.Write(hfmt.format("h,k,l","d-space","F_obs","sig(Fobs)","F_calc","phase","twin"))
        self.Write(80*'=')
        for (
            h,k,l,twin,dsp,Fobs,sigFobs,Fcalc,FobsT,FcalcT,phase,Icorr
            ) in hist.data['data'][1]['RefList']:
            if twin > 0:
                self.Write(fmt.format(hklfmt.format(h,k,l),dsp,Fobs,sigFobs,Fcalc,phase,twin))
        self.CloseFile()
        print(hist.name+' written to file '+self.fullpath)
        
    def Exporter(self,event=None):
        '''Export a set of single crystal data as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile='default' # base name on the GPX file name
            ): return 
        self.OpenFile()
        hist = list(self.histnam)[0] # there should only be one histogram, in any case take the 1st
        histblk = self.Histograms[hist]
        hklfmt = "{:.0f},{:.0f},{:.0f}"
        hfmt = "{:>10s} {:>8s} {:>12s} {:>12s} {:>12s} {:>7s} {:>6s}"
        fmt = "{:>10s} {:8.3f} {:12.2f} {:12.4f} {:12.2f} {:7.2f} {:6.0f}"
        self.Write(80*'=')
        self.Write(hfmt.format("h,k,l","d-space","F_obs","sig(Fobs)","F_calc","phase","twin"))
        self.Write(80*'=')
        for (
            h,k,l,twin,dsp,Fobs,sigFobs,Fcalc,FobsT,FcalcT,phase,Icorr
            ) in histblk['Data']['RefList']:
            if twin > 0:
                self.Write(fmt.format(hklfmt.format(h,k,l),dsp,Fobs,sigFobs,Fcalc,phase,twin))
        self.CloseFile()
        print(hist+' written to file '+self.fullpath)                        

