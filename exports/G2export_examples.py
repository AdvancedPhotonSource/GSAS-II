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
*Module G2export_examples: Examples*
-------------------------------------------

Code to demonstrate how GSAS-II data export routines are created. The
classes defined here, :class:`ExportPhaseText`, 
:class:`ExportSingleText`, :class:`ExportPowderReflText`, 
and :class:`ExportPowderText` each demonstrate a different type
of export. Also see :class:`G2export_map.ExportMapASCII` for an
example of a map export.

'''
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
#import GSASIIgrid as G2gd
#import GSASIIstrIO as G2stIO
import GSASIImath as G2mth
#import GSASIIlattice as G2lat
#import GSASIIspc as G2spc
#import GSASIIphsGUI as G2pg
#import GSASIIstrMain as G2stMn

class ExportPhaseText(G2IO.ExportBaseclass):
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
        if self.ExportSelect( # set export parameters
            AskFile=True     # prompt the user for a file name
            ): return 
        self.OpenFile(self.filename)
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
            print('Phase '+str(phasenam)+' written to file '+str(self.filename))                        
        self.CloseFile()

class ExportPowderText(G2IO.ExportBaseclass):
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
        self.multiple = False # only allow one histogram to be selected

    def Exporter(self,event=None):
        '''Export a set of powder data as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile=False # use the default file name
            #AskFile=True
            ): return 
        self.OpenFile()
        hist = self.histnam[0] # there should only be one histogram, in any case take the 1st
        histblk = self.Histograms[hist]
        fmt = 2*"{:12.3f} " + "{:12.5f} " + 2*"{:12.3f} "
        hfmt = 5*"{:>12s} "
        self.Write(hfmt.format("x","y_obs","weight","y_calc","y_bkg"))
        for x,yobs,yw,ycalc,ybkg,obsmcalc in zip(histblk['Data'][0],
                                                 histblk['Data'][1],
                                                 histblk['Data'][2],
                                                 histblk['Data'][3],
                                                 histblk['Data'][4],
                                                 histblk['Data'][5],
                                                 ):
            self.Write(fmt.format(x,yobs,yw,ycalc,ybkg))
        self.CloseFile()
        print(str(hist)+' written to file '+str(self.filename))                        
class ExportPowderReflText(G2IO.ExportBaseclass):
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
            AskFile=False # use the default file name
            #AskFile=True
            ): return 
        self.OpenFile()
        hist = self.histnam[0] # there should only be one histogram, in any case take the 1st
        histblk = self.Histograms[hist]
        hklfmt = "{:.0f},{:.0f},{:.0f}"
        hfmt = "{:>8s} {:>8s} {:>12s} {:>12s} {:>7s} {:>6s}"
        fmt = "{:>8s} {:8.3f} {:12.3f} {:12.3f} {:7.2f} {:6.0f}"
        for phasenam in histblk['Reflection Lists']:
            self.Write('\nPhase '+str(phasenam))
            self.Write(80*'=')
            self.Write(hfmt.format("h,k,l","2-theta","F_obs","F_calc","phase","mult"))
            self.Write(80*'=')
            for (
                h,k,l,mult,dsp,pos,sig,gam,Fobs,Fcalc,phase,eqlist,phaselist,Icorr,FFdict
                ) in histblk['Reflection Lists'][phasenam]:
                self.Write(fmt.format(hklfmt.format(h,k,l),pos,Fobs,Fcalc,phase,mult))
        self.CloseFile()
        print(str(hist)+'reflections written to file '+str(self.filename))                        

class ExportSingleText(G2IO.ExportBaseclass):
    '''Used to create a text file with single crystal reflection data

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

    def Exporter(self,event=None):
        '''Export a set of single crystal data as a text file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        if self.ExportSelect( # set export parameters
            AskFile=False # use the default file name
            #AskFile=True
            ): return 
        self.OpenFile()
        hist = self.histnam[0] # there should only be one histogram, in any case take the 1st
        histblk = self.Histograms[hist]
        hklfmt = "{:.0f},{:.0f},{:.0f}"
        hfmt = "{:>10s} {:>8s} {:>12s} {:>12s} {:>12s} {:>7s} {:>6s}"
        fmt = "{:>10s} {:8.3f} {:12.2f} {:12.4f} {:12.2f} {:7.2f} {:6.0f}"
        self.Write(80*'=')
        self.Write(hfmt.format("h,k,l","d-space","F_obs","sig(Fobs)","F_calc","phase","mult"))
        self.Write(80*'=')
        for (
            h,k,l,mult,dsp,Fobs,sigFobs,Fcalc,FobsT,FcalcT,phase,eqlist,phaselist,Icorr,FFdict
            ) in histblk['Data']:
            self.Write(fmt.format(hklfmt.format(h,k,l),dsp,Fobs,sigFobs,Fcalc,phase,mult))
        self.CloseFile()
        print(str(hist)+' written to file '+str(self.filename))                        

