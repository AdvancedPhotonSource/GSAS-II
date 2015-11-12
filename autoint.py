import os
import wx
import copy
import glob
import re
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIctrls as G2G
import GSASIIgrid as G2gd
import GSASIIimgGUI as G2imG
'''
Define a class to be used for Andrey's AutoIntegration process
'''

class AutoIntFrame(wx.Frame):
    '''Creates a wx.Frame window for the Image AutoIntegration.
    The intent is that this will be used as a non-modal dialog window.
    
    Implements a Start button that morphs into a pause and resume button.
    This button starts a processing loop that is repeated every
    :meth:`PollTime` seconds.

    :param wx.Frame G2frame: main GSAS-II frame
    :param float PollTime: frequency in seconds to repeat calling the
      processing loop. (Default is 3.0 seconds.)
    '''
    def OnTimerLoop(self,event):
        '''A method that is called every :meth:`PollTime` seconds that is
        used to check for new files and process them. This is called only
        after the "Start" button is pressed (when its label reads "Pause").
        '''
        G2frame = self.G2frame
        try:
            self.currImageList = sorted(
                glob.glob(os.path.join(self.imagedir,self.params['filter'])))
        except IndexError:
            self.currImageList = []
            return

        # index of image tree items by file name:
        imageDict = {G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img))[1]:
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
                     for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG '])}
        createdImageIdList = []
        # loop over files that are found, reading in new ones
        for newImage in self.currImageList:
            if newImage in self.G2frame.IntegratedList: continue # already integrated
            # has this image already been loaded?
            if newImage not in imageDict:
                Comments,Data,Npix,Image = G2IO.GetImageData(G2frame,newImage)
                if not Npix:
                    print('problem reading '+newImage)
                    continue
                G2IO.LoadImage2Tree(newImage,G2frame,Comments,Data,Npix,Image)
            else:
                G2frame.Image = imageDict[newImage]
            # update controls from master
            controlsDict = G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
            controlsDict.update(self.ImageControls)
            # update masks from master
            ImageMasks = G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
            createdImageIdList.append(G2frame.Image) # save IMG Id
            self.G2frame.IntegratedList.append(newImage) # save name of image so we don't process it again
            #print('debug: read '+newImage)

        # now integrate the images we have read
        for newImagId in createdImageIdList:
            G2frame.Image = newImagId
            G2frame.PickId = G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls')
            #  integrate in this entry
            size,imagefile = G2frame.PatternTree.GetItemPyData(G2frame.Image)
            masks = G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
            data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
            G2frame.ImageZ = G2IO.GetImageData(G2frame,imagefile,True)
            self.oldImagefile = '' # mark image as changed; reread as needed
            # simulate a Image Controls press, since that is where the
            # integration is hidden
            G2imG.UpdateImageControls(G2frame,data,masks,IntegrateOnly=True)
            # split name and control number
            s = re.split(r'(\d+)\Z',os.path.split(os.path.splitext(imagefile)[0])[1])
            namepre = s[0]
            if len(s) > 1:
                namenum = s[1]
            else:
                namenum = ''
            # write out the images in the selected formats and save the names,
            # reset will delete them
            for Id in G2frame.IntgOutList:
                treename = G2frame.PatternTree.GetItemText(Id)
                self.G2frame.AutointPWDRnames.append(treename)
                Sdata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Sample Parameters'))
                # determine the name for the current file
                fileroot = namepre
                if len(G2frame.IntgOutList) > 1:
                    fileroot += "_AZM"
                    if 'Azimuth' in Sdata:
                        fileroot += str(int(10*Sdata['Azimuth']))
                    fileroot += "_" 
                fileroot += namenum
                # loop over selected formats
                for dfmt in self.fmtlist:
                    if not self.params[dfmt[1:]]: continue
                    if self.params['SeparateDir']:
                        subdir = dfmt[1:]
                    else:
                        subdir = ''
                    fil = os.path.join(self.params['outdir'],subdir,fileroot)
                    print('writing file '+fil+dfmt)
                    G2IO.ExportPowder(G2frame,treename,fil,dfmt)
        
        if GSASIIpath.GetConfigValue('debug'):
            import datetime
            print ("Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))

    def StartLoop(self):
        '''Save current Image params for use in future integrations
        also label the window so users understand whatis being used
        '''
        print '\nStarting new autointegration\n'
        # show current IMG base
        self.ControlBaseLbl.SetLabel(self.G2frame.PatternTree.GetItemText(self.G2frame.Image))
        if self.params['Mode'] == 'file':
            'get file info'
            GSASIIpath.IPyBreak()
        else:
            # load copy of Image Controls from current image and clean up
            # items that should not be copied
            self.ImageControls = copy.deepcopy(
                self.G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(
                    self.G2frame,self.G2frame.Image, 'Image Controls')))
            self.ImageControls['showLines'] = True
            self.ImageControls['ring'] = []
            self.ImageControls['rings'] = []
            self.ImageControls['ellipses'] = []
            self.ImageControls['setDefault'] = False
            del self.ImageControls['range']
            del self.ImageControls['size']
            del self.ImageControls['GonioAngles']
            # load copy of Image Masks, keep thresholds
            self.ImageMasks = copy.deepcopy(
                self.G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(self.G2frame,self.G2frame.Image, 'Masks')))
            self.Thresholds = self.ImageMasks['Thresholds'][:]
        # make sure all output directories exist
        if self.params['SeparateDir']:
            for dfmt in self.fmtlist:
                if not self.params[dfmt[1:]]: continue
                dir = os.path.join(self.params['outdir'],dfmt[1:])
                if not os.path.exists(dir): os.makedirs(dir)
        else:
            if not os.path.exists(self.params['outdir']):
                os.makedirs(self.params['outdir'])
        if self.Reset: # after Reset has been pressed, delete all PWDR items
            # created after last Start was pressed 
            G2frame = self.G2frame
            idlist = []
            item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while item:
                itemName = G2frame.PatternTree.GetItemText(item)
                if itemName in self.G2frame.AutointPWDRnames:
                    idlist.append(item)
                item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            for item in idlist:
                G2frame.PatternTree.Delete(item)
        self.Reset = False
        self.G2frame.AutointPWDRnames = [] # list of created PWDR tree item names

    def __init__(self,G2frame,PollTime=60.0):
        def OnStart(event):
            '''Called when the start button is pressed. Changes button label 
            to Pause. When Pause is pressed the label changes to Resume.
            When either Start or Resume is pressed, the processing loop
            is started. When Pause is pressed, the loop is stopped. 
            '''
            #print self.params # for debug

            # check inputs before starting
            err = ''
            #if not any([self.params[fmt] for fmt in self.fmtlist]):
            #    err += '\nPlease select at least one output format\n'
            if (self.params['Mode'] == 'file' and not
                    os.path.exists(self.params['IMGfile'])):
                err += '\nThe image controls file could not be found\n'
            if (self.params['Mode'] == 'file' and
                not self.params['IgnoreMask']
                ) and not os.path.exists(self.params['MaskFile']): 
                err += '\nThe mask file could not be found\n'
            if err:
                G2G.G2MessageBox(self,err)
                return
            # change button label
            if btnstart.GetLabel() != 'Pause':
                btnstart.SetLabel('Pause')
                if self.timer.IsRunning(): self.timer.Stop()
                self.StartLoop()
                self.OnTimerLoop(None) # run once immediately and again after delay
                self.timer.Start(int(1000*PollTime),oneShot=False)
                self.Status.SetStatusText('Press Pause to delay integration or Reset to prepare to reintegrate all images')
            else:
                btnstart.SetLabel('Resume')
                if self.timer.IsRunning(): self.timer.Stop()
                print('\nPausing autointegration\n')
                self.Status.SetStatusText('Press Resume to continue integration or Reset to prepare to reintegrate all images')

        def OnReset(event):
            '''Called when Reset button is pressed. This stops the
            processing loop and resets the list of integrated files so
            all images can be reintegrated. 
            '''
            btnstart.SetLabel('Restart')
            self.Status.SetStatusText('Press Restart to reload and re-integrate images matching filter')
            if self.timer.IsRunning(): self.timer.Stop()
            self.Reset = True
            self.G2frame.IntegratedList = []
            
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
            if btn1 == event.GetEventObject():
                ext = '.imctrl'
                title = 'Image control'
            else:
                ext = '.immask'
                title = 'Image masks'                
            dlg = wx.FileDialog(
                self, 'Select name for '+title+' file to read',
                '.', '',
                title+'file (*'+ext+')|*'+ext,
                wx.OPEN|wx.CHANGE_DIR)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    # make sure extension is correct
                    #filename = os.path.splitext(filename)[0]+ext
                    if btn1 == event.GetEventObject():
                        fInp1.SetValue(filename)
                    else:
                        fInp2.SetValue(filename)
                else:
                    filename = None
            finally:
                dlg.Destroy()
                
        def OnRadioSelect(event):
            '''Respond to a radiobutton selection and enable or
            disable widgets accordingly. Also gets called when the
            "Don't Use" flag for Mask use is called. 
            '''
            lbl1.Disable()
            fInp1.Disable()
            btn1.Disable()
            lbl2.Disable()
            fInp2.Disable()
            ign2.Disable()
            btn2.Disable()
            if r2.GetValue():
                self.params['Mode'] = 'file'
                fInp1.Enable()
                btn1.Enable()
                lbl1.Enable()
                ign2.Enable()
                if not self.params['IgnoreMask']:
                    fInp2.Enable()
                    btn2.Enable()
                    lbl2.Enable()
            else:
                self.params['Mode'] = 'active'
        ##################################################
        # beginning of __init__ processing
        ##################################################
        self.G2frame = G2frame
        self.params = {}
        self.Reset = False
        self.params['IMGfile'] = ''
        self.params['MaskFile'] = ''
        self.params['IgnoreMask'] = True
        self.fmtlist = G2IO.ExportPowderList(G2frame)
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnTimerLoop)

        controlsId = G2frame.PatternTree.GetSelection()
        size,imagefile = G2frame.PatternTree.GetItemPyData(G2frame.Image)
        self.imagedir,fileroot = os.path.split(imagefile)
        self.params['filter'] = '*'+os.path.splitext(fileroot)[1]
        self.params['outdir'] = os.path.abspath(self.imagedir)
        wx.Frame.__init__(self, G2frame,title='Automatic Integration')
        self.Status = self.CreateStatusBar()
        self.Status.SetStatusText('Press Start to load and integrate images matching filter')
        mnpnl = wx.Panel(self)
        mnsizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Integration based on: '))
        self.ControlBaseLbl = wx.StaticText(mnpnl, wx.ID_ANY,'?')
        self.ControlBaseLbl.SetLabel(G2frame.PatternTree.GetItemText(G2frame.Image))
        sizer.Add(self.ControlBaseLbl)
        mnsizer.Add(sizer,0,wx.ALIGN_LEFT,1)
        # file filter stuff
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Image filter'))
        flterInp = G2G.ValidatedTxtCtrl(mnpnl,self.params,'filter')
        sizer.Add(flterInp)
        mnsizer.Add(sizer,0,wx.ALIGN_RIGHT,1)
        # box for integration controls & masks input
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Integration Controls/Masks source")
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        r1 = wx.RadioButton(mnpnl, wx.ID_ANY, "Use Active Image",
                            style = wx.RB_GROUP)
        r1.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
        lblsizr.Add(r1)
        r1.SetValue(True)
        r2 = wx.RadioButton(mnpnl, wx.ID_ANY, "Use from file(s)")
        lblsizr.Add(r2)
        r2.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
        r2.Disable()         # deactivate this until implemented
        # Image controls file
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((20,-1))
        lbl1 = wx.StaticText(mnpnl, wx.ID_ANY,'IMG control file: ')
        sizer.Add(lbl1)
        fInp1 = G2G.ValidatedTxtCtrl(mnpnl,self.params,'IMGfile',
                                       notBlank=False,size=(300,-1))
        sizer.Add(fInp1)
        btn1 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
        btn1.Bind(wx.EVT_BUTTON, OnBrowse)
        sizer.Add(btn1)
        lblsizr.Add(sizer)
        # Masks input file
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((20,-1))
        lbl2 = wx.StaticText(mnpnl, wx.ID_ANY,'Mask file: ')
        sizer.Add(lbl2)
        fInp2 = G2G.ValidatedTxtCtrl(mnpnl,self.params,'MaskFile',
                                       notBlank=False,size=(300,-1))
        sizer.Add(fInp2)
        ign2 = G2G.G2CheckBox(mnpnl,"Don't use",self.params,'IgnoreMask',
                              OnChange=OnRadioSelect)
        sizer.Add(ign2)
        btn2 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
        btn2.Bind(wx.EVT_BUTTON, OnBrowse)
        sizer.Add(btn2)
        lblsizr.Add(sizer)
        mnsizer.Add(lblsizr)

        # box for output selections
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Output settings")
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Write to: '))
        fInp3 = G2G.ValidatedTxtCtrl(mnpnl,self.params,'outdir',
                                       notBlank=False,size=(300,-1))
        sizer.Add(fInp3)
        btn3 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
        btn3.Bind(wx.EVT_BUTTON, OnBrowse)
        sizer.Add(btn3)
        lblsizr.Add(sizer)
        #lblsizr.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Select format(s): '))
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Select format(s): '))
        for dfmt in self.fmtlist:
            fmt = dfmt[1:]
            self.params[fmt] = False
            btn = G2G.G2CheckBox(mnpnl,dfmt,self.params,fmt)
            sizer.Add(btn)
        lblsizr.Add(sizer)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Separate dir for each format: '))
        self.params['SeparateDir'] = False
        sizer.Add(G2G.G2CheckBox(mnpnl,'',self.params,'SeparateDir'))
        lblsizr.Add(sizer)
        mnsizer.Add(lblsizr,0,wx.ALIGN_CENTER,1)

        # buttons on bottom
        mnsizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'AutoIntegration controls'),0,wx.TOP,5)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((20,-1))
        btnstart = wx.Button(mnpnl,  wx.ID_ANY, "Start")
        btnstart.Bind(wx.EVT_BUTTON, OnStart)
        sizer.Add(btnstart)
        btnstop = wx.Button(mnpnl,  wx.ID_ANY, "Reset")
        btnstop.Bind(wx.EVT_BUTTON, OnReset)
        sizer.Add(btnstop)
        sizer.Add((20,-1),wx.EXPAND,1)
        btnquit = wx.Button(mnpnl,  wx.ID_ANY, "Close")
        btnquit.Bind(wx.EVT_BUTTON, OnQuit)
        sizer.Add(btnquit)
        sizer.Add((20,-1))
        mnsizer.Add(sizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP,5)
        
        # finish up window
        mnpnl.SetSizer(mnsizer)
        OnRadioSelect(None) # disable widgets
        mnsizer.Fit(self)
        self.CenterOnParent()
        self.Show()  

if __name__ == '__main__':
    app = wx.PySimpleApp()
    G2frame = wx.Frame(None) # create a top-level frame as a stand-in for the GSAS-II data tree
    G2frame.Show() 
    frm = AutoIntFrame(G2frame) # test the one above
    app.MainLoop()
 
