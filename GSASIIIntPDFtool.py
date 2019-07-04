'''Independent-running GSAS-II based auto-integration program with minimal
GUI, no visualization but intended to implement significant levels of 
parallelization.
'''
# Autointegration from 
# $Id: GSASIIimgGUI.py 3926 2019-04-23 18:11:07Z toby $
# hacked for stand-alone use
#
# idea: select image file type & set filter from that
#
from __future__ import division, print_function
import os
import copy
import glob
import time
import re
import math
import sys
import wx
import wx.lib.mixins.listctrl  as  listmix
import wx.grid as wg
import numpy as np
import GSASIIpath
GSASIIpath.SetBinaryPath(True)
GSASIIpath.SetVersionNumber("$Revision: $")
import GSASIIIO as G2IO
import GSASIIctrlGUI as G2G
import GSASIIobj as G2obj
import GSASIIpy3 as G2py3
import GSASIIimgGUI as G2imG
import GSASIIfiles as G2fil
import GSASIIscriptable as G2sc
import multiprocessing as mp

class AutoIntFrame(wx.Frame):
    '''Creates a wx.Frame window for the Image AutoIntegration.
    The intent is that this will be used as a non-modal dialog window.
    
    Implements a Start button that morphs into a pause and resume button.
    This button starts a processing loop that is repeated every
    :meth:`PollTime` seconds.

    :param wx.Frame G2frame: main GSAS-II frame
    :param float PollTime: frequency in seconds to repeat calling the
      processing loop. (Default is 30.0 seconds.)
    '''
    def __init__(self,G2frame,PollTime=30.0):
        def OnStart(event):
            '''Called when the start button is pressed. Changes button label 
            to Pause. When Pause is pressed the label changes to Resume.
            When either Start or Resume is pressed, the processing loop
            is started. When Pause is pressed, the loop is stopped. 
            '''
            self.Pause = False
            # change button label
            if self.btnstart.GetLabel() != 'Pause':
                self.btnstart.SetLabel('Pause')
                self.Status.SetStatusText('Press Pause to delay integration or Reset to prepare to reintegrate all images')
                if self.timer.IsRunning(): self.timer.Stop()
                self.PreventTimerReEntry = False
                if self.StartLoop():
                    G2G.G2MessageBox(self,'Error in setting up integration. See console')
                    return
                # delete pool of slaved interpreters in case input has changed
                if self.MPpool:
                    self.MPpool.terminate()
                    self.MPpool = None
#                if self.ncores >= 1: # debug
                if self.ncores > 1:
#                    print('Creating pool of ',self.ncores,'cores') # debug
                    self.MPpool = mp.Pool(self.ncores)
                self.OnTimerLoop(None) # run once immediately
                if not self.Pause:
                    # no pause, so start timer to check for new files
                    self.timer.Start(int(1000*PollTime),oneShot=False)
                    return
            # we will get to this point if Paused
            self.OnPause()
            
        def OnReset(event):
            '''Called when Reset button is pressed. This stops the
            processing loop and resets the list of integrated files so
            all images can be reintegrated. 
            '''
            self.btnstart.SetLabel('Restart')
            self.Status.SetStatusText('Press Restart to reload and re-integrate images matching filter')
            if self.timer.IsRunning(): self.timer.Stop()
            self.Reset = True
            self.Pause = True
            self.ProcessedList = []
            self.ShowMatchingFiles(None)
            
        def OnQuit(event):
            '''Stop the processing loop and close the Frame
            '''
            if self.timer.IsRunning(): self.timer.Stop() # make sure we stop first
            wx.CallAfter(self.Destroy)
            
        def OnBrowse(event):
            '''Responds when the Browse button is pressed to load a file.
            The routine determines which button was pressed and gets the
            appropriate file type and loads it into the appropriate place
            in the dict.
            '''
            if btn3 == event.GetEventObject():
                dlg = wx.DirDialog(
                    self, 'Select directory for output files',
                    self.params['outdir'],wx.DD_DEFAULT_STYLE)
                dlg.CenterOnParent()
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        self.params['outdir'] = dlg.GetPath()
                        fInp3.SetValue(self.params['outdir'])
                finally:
                    dlg.Destroy()
                return
                                
        def showPDFctrls(event):
            '''Called to show or hide AutoPDF widgets. Note that TextCtrl's 
            must be included in the sizer layout with .Show(True) before 
            .Show(False) will work properly.
            '''
            TestInput()
            if len(self.pdfList) == 0 or self.params['TableMode']:
                lbl4b.SetValue(False)
                self.params['ComputePDF'] = False
                lbl4b.Enable(False)
            else:
                lbl4b.Enable(True)

            for w in [self.pbkg[i][4] for i in (0,1,2)]:
                w.Enable(self.params['ComputePDF'])
                w.Show(self.params['ComputePDF'])
            #for o in pdfwidgets+[self.pbkg[i][2] for i in (0,1,2)]:
            for o in pdfwidgets+[self.pdfSel]+[self.pbkg[i][1] for i in range(3)]:
                o.Enable(self.params['ComputePDF'])
            if self.params['ComputePDF']:
                c = "black"
            else:
                c = "gray"
            for l in ([lbl4,lbl4a,lbl5,lbl5a,lbl5b]
                          +[self.pbkg[i][0] for i in (0,1,2)]
                          +[self.pbkg[i][5] for i in (0,1,2)]):
                l.SetForegroundColour(c)
            checkPDFselection()
            self.SendSizeEvent()
                                    
        def GetGPXInputFile(event):
            'Get and read input from input GPX file'
            pth = self.params['readdir']
            lbl = 'Select a project input file'
            filtyp = 'project file (.gpx)|*.gpx'
            if event is not None:
                dlg = wx.FileDialog(self,lbl, pth,
                        style=wx.FD_OPEN,wildcard=filtyp)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    self.gpxin[3] = dlg.GetPath()
                dlg.Destroy()
            if not os.path.exists(self.gpxin[3]):
                G2G.G2MessageBox(self,'Error: file {} not found.'.format(self.gpxin[3]))
                return
            SetGPXInputFile()
        def SetGPXInputFile():
            gpx = G2sc.G2Project(self.gpxin[3])
            self.imgList = gpx.images()
            self.histList = gpx.histograms()
            self.pdfList = gpx.pdfs()
            if not self.imgList:
                G2G.G2MessageBox(self,'Error: no images in {}.'.format(self.gpxin[3]))
                return
            self.gpxin[1].SetValue(self.gpxin[3])            
            self.imprm[1].Clear()
            self.imprm[1].AppendItems([i.name for i in self.imgList])
            if len(self.imgList) == 1:
                self.imprm[1].SetSelection(0)
            self.maskfl[1].Clear()
            self.maskfl[1].AppendItems(['']+[i.name for i in self.imgList])
            for i in range(3):
                self.pbkg[i][1].Clear()
                self.pbkg[i][1].AppendItems(['']+[i.name for i in self.histList])
            self.pdfSel.Clear()
            self.pdfSel.AppendItems([i.name for i in self.pdfList])
            showPDFctrls(None)

        def TestInput(*args,**kwargs):
            'Determine if the start button should be enabled'
            writingSomething = False
            writingPDF = False
            # is at least one output file selected?
            for dfmt in self.fmtlist:
                fmt = dfmt[1:]
                if self.params['outsel'][fmt]:
                    writingSomething = True
                    break
            if self.params['ComputePDF']:
                for fmt in self.PDFformats:
                    if self.params['outsel'][fmt]:
                        writingPDF = writingSomething = True
                        break
            if not writingSomething: 
                self.EnableIntButtons(False)
                return
            # do we have integration input?
            if not self.params['TableMode']:
                if not self.gpxin[3]:
                    self.EnableIntButtons(False)
                    return
                elif not os.path.exists(self.gpxin[3]):
                    self.EnableIntButtons(False)
                    return
                elif len(self.imgList) == 0:
                    self.EnableIntButtons(False)
                    return
                else:
                    self.EnableIntButtons(True)
            else:
                if self.ImgTblParms:
                    self.EnableIntButtons(True)
                else:
                    self.EnableIntButtons(False)
                    return
            # do we have PDF input, if requested
            if self.params['ComputePDF']:
                if len(self.pdfList) == 0 or not writingPDF:
                    self.EnableIntButtons(False)
                elif 'Error' in self.formula:
                    self.EnableIntButtons(False)
                
        def checkPDFselection():
            'Read PDF entry from input GPX file & show in GUI'
            pdfEntry = self.pdfSel.GetStringSelection()
            if not self.pdfList:
                self.params['ComputePDF'] = False
                lbl4b.Enable(False)
                return
            else:
                lbl4b.Enable(True)
            if pdfEntry not in [i.name for i in self.pdfList]:
                # strange something -- not selected
                self.pdfSel.SetSelection(0)
                pdfEntry = self.pdfSel.GetStringSelection()
            if self.gpxInp is None or self.gpxInp.filename != self.gpxin[3]:
                self.gpxInp = G2sc.G2Project(self.gpxin[3])
            try: 
                PDFobj = self.gpxInp.pdf(pdfEntry)
            except KeyError:
                print("PDF entry not found: {}".format(pdfEntry))
                return
            histNames = [i.name for i in self.histList]
            for i,lbl in enumerate(('Sample Bkg.','Container',
                                    'Container Bkg.')):
                self.pbkg[i][4].SetValue(str(PDFobj.data['PDF Controls'][lbl]['Mult']))
                self.pbkg[i][6] = PDFobj.data['PDF Controls'][lbl]['Mult']
                try:
                    i = 1 + histNames.index(PDFobj.data['PDF Controls'][lbl]['Name'])
                    self.pbkg[i][1].SetSelection(i)
                except ValueError:
                    i = 0
                    self.pbkg[i][1].SetSelection(0)
                    if PDFobj.data['PDF Controls'][lbl]['Name']:
                        print('PDF {} hist entry {} not found'.format(
                            lbl,PDFobj.data['PDF Controls'][lbl]['Name']))
                        PDFobj.data['PDF Controls'][lbl]['Name'] = ''
            self.formula = ''
            for el in PDFobj.data['PDF Controls']['ElList']:
                i = PDFobj.data['PDF Controls']['ElList'][el]['FormulaNo']
                if i <= 0:
                    continue
                elif i == 1:
                    if self.formula: self.formula += ' '
                    self.formula += '{}'.format(el)
                else:
                    if self.formula: self.formula += ' '
                    self.formula += '{}({:.1f})'.format(el,i)
            if not self.formula:
                self.formula = 'Error: no chemical formula'
            lbl5b.SetLabel(self.formula)
            TestInput()
                    
        def ShowbyMode():
            'create table or non-table integration section of GUI'
            intPrmSizer.Clear(True)
            if not self.params['TableMode']:
                sizer = wx.BoxSizer(wx.HORIZONTAL)
                self.gpxin[0] = wx.StaticText(mnpnl, wx.ID_ANY,'Project (gpx) file:')
                sizer.Add(self.gpxin[0])
                self.gpxin[1] = G2G.ValidatedTxtCtrl(mnpnl,self.gpxin,3,
                                OKcontrol=TestInput,OnLeave=TestInput)
                sizer.Add(self.gpxin[1],1,wx.EXPAND,1)
                self.gpxin[2] = wx.Button(mnpnl, wx.ID_ANY, "Browse")
                self.gpxin[2].Bind(wx.EVT_BUTTON, GetGPXInputFile)
                sizer.Add(self.gpxin[2],0,wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
                intPrmSizer.Add(sizer,0,wx.EXPAND,0)

                sizer = wx.BoxSizer(wx.HORIZONTAL)
                self.imprm[0] = wx.StaticText(mnpnl, wx.ID_ANY,'Image parms from:')
                sizer.Add(self.imprm[0])
                self.imprm[3] = 0
                self.imprm[1] = G2G.G2ChoiceButton(mnpnl,[''],self.imprm,3,onChoice=TestInput)
                sizer.Add(self.imprm[1],1,wx.EXPAND,1)
                intPrmSizer.Add(sizer,0,wx.EXPAND,0)

                sizer = wx.BoxSizer(wx.HORIZONTAL)
                self.maskfl[0] = wx.StaticText(mnpnl, wx.ID_ANY,'Mask parms from:')
                sizer.Add(self.maskfl[0])
                self.maskfl[3] = 0
                self.maskfl[1] = G2G.G2ChoiceButton(mnpnl,[''],self.maskfl,3,onChoice=TestInput)
                sizer.Add(self.maskfl[1],1,wx.EXPAND,1)
                intPrmSizer.Add(sizer,0,wx.EXPAND,0)
            else:
                sizer = wx.BoxSizer(wx.HORIZONTAL)
                self.table = [None,None,None,None]
                self.table[0] = wx.Button(mnpnl,  wx.ID_ANY, "Create table")
                sizer.Add(self.table[0],0,wx.ALIGN_LEFT|wx.ALL,5)
                self.table[0].Bind(wx.EVT_BUTTON, OnTableButton)
                self.table[1] = wx.Button(mnpnl,  wx.ID_ANY, "Read table")
                sizer.Add(self.table[1],0,wx.ALIGN_LEFT|wx.ALL,5)
                self.table[1].Bind(wx.EVT_BUTTON, OnTableButton)
                self.table[2] = wx.Button(mnpnl,  wx.ID_ANY, "Edit table")
                sizer.Add(self.table[2],0,wx.ALIGN_LEFT|wx.ALL,5)
                self.table[2].Bind(wx.EVT_BUTTON, OnTableButton)
                #self.table[3] = wx.Button(mnpnl,  wx.ID_ANY, "Save table")
                #sizer.Add(self.table[3],0,wx.ALIGN_LEFT|wx.ALL,5)
                #self.table[3].Bind(wx.EVT_BUTTON, OnTableButton)
                intPrmSizer.Add(sizer,0,wx.EXPAND,0)
            # enable/disable based on status of files/table
            TestInput()
            mnsizer.Fit(self)
            
        def OnTableButton(event):
            '''Called to edit/create the distance-dependent parameter look-up table.
            '''
            pth = self.params['readdir']
            readFileList = []
            parms,fileList = [], []
            if event.GetEventObject() == self.table[0]:
                dlg = wx.FileDialog(self, 'Build new table by selecting image control files', pth,
                    style=wx.FD_OPEN| wx.FD_MULTIPLE,
                    wildcard='Image control files (.imctrl)|*.imctrl')
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    readFileList = dlg.GetPaths()
                dlg.Destroy()
                if len(readFileList) <= 0: return
            elif event.GetEventObject() == self.table[1]:
                dlg = wx.FileDialog(self, 'Reload table by selecting saved file', pth,
                    style=wx.FD_OPEN,
                    wildcard='Integration table (*.imtbl)|*.imtbl')
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    readFileList = [dlg.GetPath()]
                dlg.Destroy()
                if len(readFileList) <= 0: return
            elif event.GetEventObject() == self.table[2]:
                parms = copy.deepcopy(self.ImgTblParms)
                fileList = copy.copy(self.IMfileList)
                if not parms:
                    G2G.G2MessageBox(self,'Create or Read table first')
                    return
            dlg = None
            try:
                dlg = G2imG.IntegParmTable(self,parms,fileList,readFileList)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    self.params['InterVals'] = SetupInterpolation(dlg)
                    self.ImgTblParms = dlg.ReadImageParmTable()
                    self.IMfileList = dlg.IMfileList
                    self.params['TableMode'] = True
                    self.params['ControlsTable'] = {}
                    self.params['MaskTable'] = {}
                    for f,m in zip(self.IMfileList,self.ImgTblParms[-1]):
                        n = os.path.split(f)[1]
                        if n in self.params['ControlsTable']:
                            print('Warning overwriting entry {}'.format(n))
                        self.params['ControlsTable'][n] = G2imG.ReadControls(f)
                        if m and os.path.exists(m):
                            self.params['MaskTable'][n] = G2imG.ReadMask(m)
                        elif m != "(none)":
                            print("Error: Mask file {} not found".format(m))
                else:
                    self.params['TableMode'] = False
                    self.params['ControlsTable'] = {}
                    self.params['MaskTable'] = {}
                    self.imageBase = G2frame.Image
            finally:
                if dlg: dlg.Destroy()
            TestInput()
            
        ##################################################
        # beginning of __init__ processing
        ##################################################
        self.G2frame = G2frame
        self.ImgTblParms = None
        self.IMfileList = None
        self.params = {}
        self.Reset = False
        self.Pause = False
        self.PreventReEntryShowMatch = False
        self.PreventTimerReEntry = False
        self.params['ControlsTable'] = {}
        self.params['MaskTable'] = {}
        G2sc.LoadG2fil()
        self.fmtlist = G2sc.exportersByExtension.get('powder',{})
        self.fmtlist['.gpx'] = None
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnTimerLoop)
        self.imageBase = G2frame.Image
        self.params['ComputePDF'] = False
        self.params['pdfDmax'] = 0.0
        self.params['pdfprm'] = ''
        self.params['optPDF'] = True
        self.params['TableMode'] = False
        self.params['outsel'] = {}
        self.formula = 'Error'
        self.pdfControls = {}
        self.imgList = []
        self.histList = []
        self.pdfList = []
        self.Scale = [1.0,]
        self.ProcessedList = [] # files that have been integrated
        self.currImageList = [] # files that will be integrated
        self.gpxInp = None
        self.gpxin = [None,None,None,'']
        self.imprm = [None,None,None,'']
        self.maskfl = [None,None,None,'']
        self.ncores = mp.cpu_count()//2
        self.MPpool = None
        self.params['readdir'] = os.getcwd()
        self.params['filter'] = '*.tif'
        self.params['outdir'] = os.getcwd()
        self.PDFformats = ('I(Q)', 'S(Q)', 'F(Q)', 'G(r)', 'PDFgui')
        #GSASIIpath.IPyBreak_base()
        
        wx.Frame.__init__(self, None, title='Automatic Integration',
                          style=wx.DEFAULT_FRAME_STYLE)
        self.Status = self.CreateStatusBar()
        self.Status.SetStatusText('Press Start to load and integrate images matching filter')
        mnpnl = wx.Panel(self)
        mnsizer = wx.BoxSizer(wx.VERTICAL)
        # box for integration controls & masks input
        intSizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((-1,-1),1,wx.EXPAND,1)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Calibration/Integration parameters'))
        sizer.Add((-1,-1),1,wx.EXPAND,1)
        intSizer.Add(sizer,1,wx.EXPAND,1)
        def ontblModeBtn(event):
            if tblModeBtn.GetValue():
                self.params['TableMode'] = True
            else:
                self.params['TableMode'] = False
            ShowbyMode()
        tblModeBtn = wx.CheckBox(mnpnl,label='Use distance lookup table')
        tblModeBtn.SetValue(False)
        tblModeBtn.Bind(wx.EVT_CHECKBOX, ontblModeBtn)
        intSizer.Add(tblModeBtn)
        mnsizer.Add(intSizer,0,wx.EXPAND,0)
        intPrmSizer = wx.BoxSizer(wx.VERTICAL)
        
        mnsizer.Add(intPrmSizer,0,wx.EXPAND,0)
        # file filter stuff
        mnsizer.Add((-1,15))
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Read images from '))
        self.readDir = G2G.ValidatedTxtCtrl(mnpnl,self.params,'readdir',
                            OnLeave=self.ShowMatchingFiles,size=(200,-1))
        sizer.Add(self.readDir,1,wx.EXPAND,1)
        btn3 = wx.Button(mnpnl, wx.ID_ANY, "Browse")
        btn3.Bind(wx.EVT_BUTTON, self.SetSourceDir)
        sizer.Add(btn3,0,wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        mnsizer.Add(sizer,0,wx.EXPAND,0)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((-1,-1),1,wx.EXPAND,1)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'  Image filter'))
        flterInp = G2G.ValidatedTxtCtrl(mnpnl,self.params,'filter',
                                        OnLeave=self.ShowMatchingFiles)
        sizer.Add(flterInp)
        mnsizer.Add(sizer,0,wx.EXPAND,0)
        
        self.ListBox = wx.ListBox(mnpnl,size=(-1,100))
        mnsizer.Add(self.ListBox,1,wx.EXPAND,1)
        self.ShowMatchingFiles(None)

        # box for output selections
        mnsizer.Add((-1,15))
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Integration output")
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Write to: '),0,wx.ALIGN_CENTER_VERTICAL)
        fInp3 = G2G.ValidatedTxtCtrl(mnpnl,self.params,'outdir',notBlank=False,size=(300,-1))
        sizer.Add(fInp3,1,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        btn3 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
        btn3.Bind(wx.EVT_BUTTON, OnBrowse)
        sizer.Add(btn3,0,wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        lblsizr.Add(sizer,0,wx.EXPAND)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Select format(s):'))
        for dfmt in self.fmtlist:
            sizer.Add((6,2)) # add a bit of extra space
            fmt = dfmt[1:]
            if fmt not in self.params['outsel']: self.params['outsel'][fmt] = False
            btn = G2G.G2CheckBox(mnpnl,dfmt,self.params['outsel'],fmt,
                                     OnChange=TestInput)
            sizer.Add(btn)
        lblsizr.Add(sizer)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Separate dir for each format: '))
        self.params['SeparateDir'] = False
        sizer.Add(G2G.G2CheckBox(mnpnl,'',self.params,'SeparateDir'))
        lblsizr.Add(sizer)
        mnsizer.Add(lblsizr,0,wx.ALIGN_CENTER|wx.EXPAND,1)
        
        mnsizer.Add((-1,15))
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "PDF settings")
        pdfwidgets = []
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Autocompute PDF:'),0,wx.ALIGN_CENTER_VERTICAL)
        lbl4b = G2G.G2CheckBox(mnpnl,'',self.params,'ComputePDF',
                                     OnChange=showPDFctrls)
        sizer.Add(lbl4b)
        lbl4a = wx.StaticText(mnpnl, wx.ID_ANY,'Max detector distance: ')
        sizer.Add(lbl4a,0,wx.ALIGN_CENTER_VERTICAL)
        fInp4a = G2G.ValidatedTxtCtrl(mnpnl,self.params,'pdfDmax',min=0.0)
        pdfwidgets.append(fInp4a)
        sizer.Add(fInp4a,0,wx.ALIGN_CENTER_VERTICAL)
        cOpt = G2G.G2CheckBox(mnpnl,'Optimize',self.params,'optPDF')
        pdfwidgets.append(cOpt)
        sizer.Add(cOpt)
        lblsizr.Add(sizer,0)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        lbl4 = wx.StaticText(mnpnl, wx.ID_ANY,'PDF control: ')
        sizer.Add(lbl4,0,wx.ALIGN_CENTER_VERTICAL)
        self.pdfSel = G2G.G2ChoiceButton(mnpnl,[''],self.params,'pdfprm',
                                       onChoice=checkPDFselection)
        sizer.Add(self.pdfSel,1,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND,1)
        lblsizr.Add(sizer,0,wx.EXPAND)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        lbl5a = wx.StaticText(mnpnl, wx.ID_ANY,'Chemical formula: ')
        sizer.Add(lbl5a)
        lbl5b = wx.StaticText(mnpnl, wx.ID_ANY,'(formula)')
        sizer.Add(lbl5b)
        lblsizr.Add(sizer,0,wx.EXPAND)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        lbl5 = wx.StaticText(mnpnl, wx.ID_ANY,'Select format(s):')
        sizer.Add(lbl5)
        for fmt in self.PDFformats:
            sizer.Add((6,2)) # add a bit of extra space
            if fmt not in self.params['outsel']: self.params['outsel'][fmt] = False
            btn = G2G.G2CheckBox(mnpnl,fmt,self.params['outsel'],fmt,
                                     OnChange=TestInput)
            sizer.Add(btn)
            pdfwidgets.append(btn)
        lblsizr.Add(sizer,0,wx.EXPAND)

        self.pbkg = 3*[None]
        for i,lbl in enumerate((' Sample bkg:',' Container:',
                                   'Container bkg:')):
            self.pbkg[i] = [None,None,None,'',None,None,-1.0]
            sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.pbkg[i][0] = wx.StaticText(mnpnl, wx.ID_ANY,lbl)
            sizer.Add(self.pbkg[i][0])
            self.pbkg[i][1] = G2G.G2ChoiceButton(mnpnl,[''],self.pbkg[i],3)
            sizer.Add(self.pbkg[i][1],1,wx.EXPAND,1)
            self.pbkg[i][5] = wx.StaticText(mnpnl, wx.ID_ANY,' mult:')
            sizer.Add(self.pbkg[i][5])
            self.pbkg[i][4] = G2G.ValidatedTxtCtrl(mnpnl,self.pbkg[i],6,
                            (6,3),typeHint=float,size=(50,-1),
                            OnLeave=TestInput,notBlank=False)
            sizer.Add(self.pbkg[i][4])
            lblsizr.Add(sizer,0,wx.EXPAND,0)
        
        mnsizer.Add(lblsizr,0,wx.ALIGN_CENTER|wx.EXPAND,1)

        # buttons on bottom
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((20,-1))
        self.btnstart = wx.Button(mnpnl,  wx.ID_ANY, "Start")
        self.btnstart.Bind(wx.EVT_BUTTON, OnStart)
        sizer.Add(self.btnstart)
        self.btnreset = wx.Button(mnpnl,  wx.ID_ANY, "Reset")
        self.btnreset.Bind(wx.EVT_BUTTON, OnReset)
        sizer.Add(self.btnreset)
        sizer.Add((20,-1),wx.EXPAND,1)
        
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,' cores:'))
        spin = wx.SpinCtrl(mnpnl, wx.ID_ANY,size=(50,-1))
        spin.SetRange(1, mp.cpu_count())
#        spin.SetRange(0, mp.cpu_count()) # debug
        spin.SetValue(self.ncores)
        def _onSpin(event):
            if event: event.Skip()
            self.ncores = spin.GetValue()
        spin.Bind(wx.EVT_SPINCTRL, _onSpin)
        spin.Bind(wx.EVT_KILL_FOCUS, _onSpin)
        sizer.Add(spin)

        sizer.Add((20,-1),wx.EXPAND,1)
        self.btnclose = wx.Button(mnpnl,  wx.ID_ANY, "Exit")
        self.btnclose.Bind(wx.EVT_BUTTON, OnQuit)
        self.EnableIntButtons(False)
        sizer.Add(self.btnclose)
        sizer.Add((20,-1))
        mnsizer.Add(sizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP,5)
        # finish up window
        mnpnl.SetSizer(mnsizer)
        #OnRadioSelect(None) # disable widgets
        mnsizer.Fit(self)
        ShowbyMode()
        if len(sys.argv) > 1:
            fil = os.path.splitext(sys.argv[1])[0]+'.gpx'
            if os.path.exists(fil):
                self.gpxin[3] = fil
                SetGPXInputFile()
        showPDFctrls(None)
    
    def SetSourceDir(self,event):
        '''Use a dialog to get a directory for image files
        '''
        dlg = wx.DirDialog(self, 'Select directory for image files',
                        self.params['readdir'],wx.DD_DEFAULT_STYLE)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.params['readdir'] = dlg.GetPath()
            self.readDir.SetValue(self.params['readdir'])
            self.ShowMatchingFiles(None)
        finally:
            dlg.Destroy()
        return
        
    def ShowMatchingFiles(self,value,invalid=False,**kwargs):
        '''Find and image files matching the image
        file directory (self.params['readdir']) and the image file filter
        (self.params['filter']) and add this information to the GUI list box
        '''
        if invalid: return
        if self.PreventReEntryShowMatch: return
        self.PreventReEntryShowMatch = True
        filmsg = ""
        self.currImageList = []
        if os.path.exists(self.params['readdir']): 
            imageList = sorted(
                glob.glob(os.path.join(self.params['readdir'],self.params['filter'])))
            if not imageList:
                msg = 'Warning: No files match search string '+os.path.join(
                        self.params['readdir'],self.params['filter'])
            else:
                for fil in imageList:
                    if fil not in self.ProcessedList:
                        filmsg += '\n  '+fil
                        self.currImageList.append(fil)
                if filmsg:
                    msg = 'Files to integrate from '+os.path.join(
                            self.params['readdir'],self.params['filter'])+filmsg
                else:
                    msg = 'No files found to process in '+self.params['readdir']
        else:
            msg = 'Warning, does not exist: '+self.params['readdir']
        if self.ProcessedList:
            msg += '\nIntegrated files:'
            for fil in self.ProcessedList:
                msg += '\n  '+fil
        self.ListBox.Clear()
        self.ListBox.AppendItems(msg.split('\n'))
        self.PreventReEntryShowMatch = False
        return
        
    def OnPause(self):
        '''Respond to Pause, changes text on button/Status line, if needed
        Stops timer
        self.Pause should already be True
        '''
        if self.timer.IsRunning(): self.timer.Stop()
        if self.btnstart.GetLabel() == 'Restart':
            return
        if self.btnstart.GetLabel() != 'Resume':
            print('\nPausing autointegration\n')
            self.btnstart.SetLabel('Resume')
            self.Status.SetStatusText(
                    'Press Resume to continue integration or Reset to prepare to reintegrate all images')
        self.Pause = True
            
    def EnableIntButtons(self,flag):
        for item in (self.btnstart,self.btnreset): item.Enable(flag)
            
    def StartLoop(self):
        '''Prepare to start autointegration timer loop. 
        Save current Image params for use in future integrations
        also label the window so users understand what is being used
        '''
        print('\nStarting new autointegration\n')
        # make sure all output directories exist
        if self.params['SeparateDir']:
            for dfmt in self.fmtlist:
                if not self.params['outsel'][dfmt[1:]]: continue
                dir = os.path.join(self.params['outdir'],dfmt[1:])
                if not os.path.exists(dir): os.makedirs(dir)
        else:
            if not os.path.exists(self.params['outdir']):
                os.makedirs(self.params['outdir'])
        if self.Reset: # special things to do after Reset has been pressed
            self.G2frame.IntegratedList = []
            wx.Yield()
            self.Reset = False
        if self.params['ComputePDF'] and self.params['SeparateDir']:
            for fmt in self.PDFformats:
                if not self.params['outsel'][fmt]: continue
                dir = os.path.join(self.params['outdir'],
                                   fmt.replace("(","_").replace(")",""))
                if not os.path.exists(dir): os.makedirs(dir)            
        return False
                
    def ArgGen(self,PDFobj,imgprms,mskprms,xydata):
        '''generator for arguments for integration/PDF calc
        '''
        for newImage in self.currImageList:
            self.Pause |= self.G2frame.PauseIntegration
            if self.Pause:
                self.OnPause()
                self.PreventTimerReEntry = False
                self.Raise()
                return
            TableMode = self.params['TableMode']
            ComputePDF = self.params['ComputePDF']
            SeparateDir = self.params['SeparateDir']
            optPDF = self.params['optPDF']
            outdir = self.params['outdir']
            calcModes = (TableMode,ComputePDF,SeparateDir,optPDF)
            InterpVals = self.params.get('InterVals')
            outputSelect = self.params['outsel']
            PDFformats = self.PDFformats
            outputModes = (outputSelect,PDFformats,self.fmtlist,outdir)
            if PDFobj:
                PDFdict = PDFobj.data
            else:
                PDFdict = None
            yield (newImage,imgprms,mskprms,xydata,PDFdict,InterpVals,calcModes,outputModes)
    def OnTimerLoop(self,event):
        '''A method that is called every :meth:`PollTime` seconds that is
        used to check for new files and process them. Integrates new images.
        Also optionally sets up and computes PDF. 
        This is called only after the "Start" button is pressed (then its label reads "Pause").
        '''
            
        if GSASIIpath.GetConfigValue('debug'):
            import datetime
            print ("DBG_Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))
        if self.PreventTimerReEntry: return
        self.PreventTimerReEntry = True
        self.ShowMatchingFiles(None)
        if not self.currImageList:
            self.PreventTimerReEntry = False
            return
        updateList = False

        # get input for integration
        imgprms = mskprms = None
        if not self.params['TableMode']:
            # read in image controls/masks, used below in loop. In Table mode
            # we will get this image-by image. 
            gpxinp = G2sc.G2Project(self.gpxin[3])
            print('reading template project',gpxinp.filename)
            img = gpxinp.image(self.imprm[1].GetStringSelection())
            imgprms = img.getControls(True)
            if self.maskfl[1].GetStringSelection().strip():
                img = gpxinp.image(self.maskfl[1].GetStringSelection())
                mskprms = img.getMasks()
        # setup shared input for PDF computation (for now will not be table mode)
        xydata = {}
        if self.params['ComputePDF']:
            pdfEntry = self.pdfSel.GetStringSelection()
            try: 
                PDFobj = gpxinp.pdf(pdfEntry)
            except KeyError:
                print("PDF entry not found: {}".format(pdfEntry))
            # update with GUI input
            for i,lbl in enumerate(('Sample Bkg.','Container',
                                    'Container Bkg.')):
                name = self.pbkg[i][1].GetStringSelection()
                try:
                    xydata[lbl] = gpxinp.histogram(name).data['data']
                except AttributeError:
                    pass
                PDFobj.data['PDF Controls'][lbl]['Mult'] = self.pbkg[i][6]
                PDFobj.data['PDF Controls'][lbl]['Name'] = name
        else:
            PDFobj = None
        if self.MPpool:
            self.MPpool.imap_unordered(ProcessImageMP,
                            self.ArgGen(PDFobj,imgprms,mskprms,xydata))
        else:
            for intArgs in self.ArgGen(PDFobj,imgprms,mskprms,xydata):
                newImage = intArgs[0]
                print('processing ',newImage)
                ProcessImage(*intArgs)
                updateList = True
        for newImage in self.currImageList:
            self.ProcessedList.append(newImage)
        if updateList: self.ShowMatchingFiles(None)
        self.PreventTimerReEntry = False
        self.Raise()
        
def ProcessImageMP(onearg):
    newImage = onearg[0]
    print('processing ',newImage)
    t0 = time.time() # debug
    ProcessImage(*onearg)
    print('ProcessImageMP time',time.time()-t0)

MapCache = {'maskMap':{}, 'ThetaAzimMap':{}, 'distanceList':[]}
'caches for TA and Mask maps'

def ProcessImage(newImage,imgprms,mskprms,xydata,PDFdict,InterpVals,calcModes,outputModes):
    '''Process one image that is read from file newImage and is integrated into 
    one or more diffraction patterns and optionally each diffraction pattern can 
    be transformed into a pair distribution function.

    :param str newImage: file name (full path) for input image
    :param dict imgprms: dict with some nested lists & dicts describing the image
      settings and integration parameters
    :param dict mskprms: dict with areas of image to be masked
    :param dict xydata: contains histogram information with about background 
      contributions, used for PDF computation (used if ComputePDF is True)
    :param PDFdict: contains PDF parameters (used if ComputePDF is True)
    :param InterpVals: contains interpolation table (used if TableMode is True)
    :param tuple calcModes: set of values for which computations are 
      performed and how
    :param tuple outputModes: determines which files are written and where
    '''
    (TableMode,ComputePDF,SeparateDir,optPDF) = calcModes
    (outputSelect,PDFformats,fmtlist,outdir) = outputModes            
    if SeparateDir:
        savedir = os.path.join(outdir,'gpx')
        if not os.path.exists(savedir): os.makedirs(savedir)
    else:
        savedir = outdir
    outgpx = os.path.join(savedir,os.path.split(os.path.splitext(newImage)[0]+'.gpx')[1])
    gpxout = G2sc.G2Project(filename=outgpx)
    print('creating',gpxout.filename)
    # looped because a file can contain multiple images
    if TableMode: # look up parameter values from table
        imgprms,mskprms = LookupFromTable(im.data['Image Controls'].get('setdist'),
                                                  InterpVals)    
    for im in gpxout.add_image(newImage):
        # apply image parameters
        im.setControls(imgprms)
        setdist = '{:.2f}'.format(im.getControls()['setdist']) # ignore differences in position less than 0.01 mm
        if setdist not in MapCache['distanceList']:
            if mskprms:
                im.setMasks(mskprms)
            else:
                im.initMasks()
            MapCache['distanceList'].append(setdist)
            MapCache['maskMap'][setdist] = G2sc.calcMaskMap(im.getControls(),
                                                            im.getMasks())
            MapCache['ThetaAzimMap'][setdist] = G2sc.calcThetaAzimMap(im.getControls())
#        else: # debug
#            print('*** reusing',setdist)
        #if mskprms:
        #    im.setMasks(mskprms)
        #else:
        #    im.initMasks()
        hists = im.Integrate(MaskMap=MapCache['maskMap'][setdist],
                             ThetaAzimMap=MapCache['ThetaAzimMap'][setdist])
        # write requested files
        for dfmt in fmtlist:
            fmt = dfmt[1:]
            if not outputSelect[fmt]: continue
            if fmtlist[dfmt] is None: continue
            if SeparateDir:
                savedir = os.path.join(outdir,fmt)
            else:
                savedir = outdir
            if not os.path.exists(savedir): os.makedirs(savedir)
            # loop over created histgrams (multiple if caked), writing them as requested
            for i,h in enumerate(hists):
                fname = h.name[5:].replace(' ','_')
                try:
                    fil = os.path.join(savedir,fname)
                    print('Wrote',h.Export(fil,dfmt))
                except Exception as msg:
                    print('Failed to write {} as {}. Error msg\n{}'
                              .format(fname,dfmt,msg))
        if ComputePDF:  # compute PDF
            for h in hists:
                pdf = gpxout.copy_PDF(PDFdict,h)
                pdf.data['PDF Controls']['Sample']['Name'] = h.name
                xydata['Sample'] = h.data['data']
                fname = h.name[5:].replace(' ','_')
                limits = h.data['Limits'][1]
                inst = h.data['Instrument Parameters'][0]
                pdf.calculate(copy.deepcopy(xydata),limits,inst)
                if optPDF:
                    for i in range(5):
                        if pdf.optimize(True,5,copy.deepcopy(xydata),limits,inst):
                            break
                    pdf.calculate(copy.deepcopy(xydata),limits,inst)
                for fmt in PDFformats:
                    if not outputSelect[fmt]: continue
                    if SeparateDir:
                        savedir = os.path.join(outdir,fmt.replace("(","_").replace(")",""))
                    else:
                        savedir = outdir
                    pdf.export(os.path.join(savedir,fname),fmt)
    if outputSelect.get('gpx'):
        gpxout.save()
    else:
        del gpxout
# Autointegration end

def SetupInterpolation(dlg):
    '''Creates an object for interpolating image parameters at a given distance value
    '''
    parms = dlg.ReadImageParmTable()
    IMfileList = dlg.IMfileList
    cols = dlg.list.GetColumnCount()
    ParmList = dlg.ParmList
    nonInterpVars = dlg.nonInterpVars
    ControlsTable = {}
    MaskTable = {}
    for f,m in zip(IMfileList,parms[-1]):
        n = os.path.split(f)[1]
        if n in ControlsTable:
            print('Warning overwriting entry {}'.format(n))
        ControlsTable[n] = G2imG.ReadControls(f)
        if m and os.path.exists(m):
            MaskTable[n] = G2imG.ReadMask(m)
        elif m != "(none)":
            print("Error: Mask file {} not found".format(m))
    return copy.deepcopy([cols, parms, IMfileList, ParmList, nonInterpVars,ControlsTable,MaskTable])

def LookupFromTable(dist,parmList):
    '''Interpolate image parameters for a supplied distance value

    :param float dist: distance to use for interpolation
    :returns: a list with 2 items:
      * a dict with interpolated parameter values,
      * the closest imctrl
    '''
    cols, parms, IMfileList, ParmList, nonInterpVars,ControlsTable,MaskTable = parmList
    x = np.array([float(i) for i in parms[0]])
    closest = abs(x-dist).argmin()
    D = {'setdist':dist}
    imctfile = IMfileList[closest]
    for c in range(1,cols-1):
        lbl = ParmList[c]
        if lbl in nonInterpVars:
            if lbl in ['outChannels',]:
                D[lbl] = int(float(parms[c][closest]))
            else:
                D[lbl] = float(parms[c][closest])
        else:
            y = np.array([float(i) for i in parms[c]])
            D[lbl] = np.interp(dist,x,y)
    # full integration when angular range is 0
    D['fullIntegrate'] = (D['LRazimuth_min'] == D['LRazimuth_max'])
    # conversion for paired values
    for a,b in ('center_x','center_y'),('LRazimuth_min','LRazimuth_max'),('IOtth_min','IOtth_max'):
        r = a.split('_')[0]
        D[r] = [D[a],D[b]]
        if r in ['LRazimuth',]:
            D[r] = [int(D[a]),int(D[b])]
        del D[a]
        del D[b]
    interpDict,imgctrl = D,imctfile
    if GSASIIpath.GetConfigValue('debug'):
        print ('DBG_interpolated values: ',interpDict)
    f = os.path.split(imgctrl)[1]
    ImageControls = ControlsTable[f]
    ImageControls.update(interpDict)
    ImageControls['showLines'] = True
    ImageControls['ring'] = []
    ImageControls['rings'] = []
    ImageControls['ellipses'] = []
    ImageControls['setDefault'] = False
    for i in 'range','size','GonioAngles':
        if i in ImageControls: del ImageControls[i]
    ImageMasks = MaskTable.get(f)
    return ImageControls,ImageMasks
    
###########################################################################
if __name__ == "__main__":
    GSASIIpath.InvokeDebugOpts()
    App = wx.App()
    class dummyClass(object):
        '''An empty class where a few values needed from parent are placed
        '''
        def __init__(self): 
            self.Image = None
            self.PauseIntegration = False
            self.TutorialImportDir = None
            self.GSASprojectfile = ''
            self.LastExportDir = ''
            self.LastGPXdir = ''
            
    G2frame = dummyClass()
    frm = AutoIntFrame(G2frame,5)
    App.GetTopWindow().Show(True)
    App.MainLoop()
