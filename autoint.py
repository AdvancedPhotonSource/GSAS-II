import os
import wx
import copy
import glob
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIctrls as G2G
import GSASIIgrid as G2gd
#('xye','fxye','xy','chi')
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
        try:
            self.currImageList = sorted(
                glob.glob(os.path.join(self.imagedir,self.params['filter'])))
        except IndexError:
            self.currImageList = []
            return

        createdImageIdList = []
        for newImage in self.currImageList:
            if newImage in self.IntegratedList: continue
            # need to read in this file
            Comments,Data,Npix,Image = G2IO.GetImageData(self.G2frame,newImage)
            if not Npix:
                print('problem reading '+newImage)
                continue
            G2IO.LoadImage(newImage,self.G2frame,Comments,Data,Npix,Image)
            controlsDict = self.G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self.G2frame,self.G2frame.Image, 'Image Controls'))
            controlsDict.update(self.ImageControls)
            ImageMasks = self.G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self.G2frame,self.G2frame.Image, 'Masks'))
            createdImageIdList.append(self.G2frame.Image)
            self.IntegratedList.append(newImage)
            print('debug: read '+newImage)

        for newImagId in createdImageIdList:
            print('debug: process '+str(newImagId))
            pass
            # need to integrate in this entry
        import datetime
        print ("Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))
        #GSASIIpath.IPyBreak_base()

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
        
    def __init__(self,G2frame,PollTime=3.0):
        def OnStart(event):
            '''Called when the start button is pressed. Changes button label 
            to Pause. When Pause is pressed the label changes to Resume.
            When either Start or Resume is pressed, the processing loop
            is started. When Pause is pressed, the loop is stopped. 
            '''
            print self.params # for debug

            # check inputs before starting
            err = ''
            #if not any([self.params[fmt] for fmt in fmtlist]):
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
            else:
                btnstart.SetLabel('Resume')
                if self.timer.IsRunning(): self.timer.Stop()
                print('\nPausing autointegration\n')

        def OnStop(event):
            '''Called when Stop button is pressed. At present this only
            stops the processing loop (same as pressing Pause except the
            label is reset to "Start".) 
            '''
            btnstart.SetLabel('Start')
            if self.timer.IsRunning(): self.timer.Stop()
            
        def OnQuit(event):
            '''Stop the processing loop and close the Frame
            '''
            # make sure we stop first
            OnStop(event)
            self.Destroy()
            
        def OnBrowse(event):
            '''Responds when the Browse button is pressed to load a file.
            The routine determines which button was pressed and gets the
            appropriate file type and loads it into the appropriate place
            in the dict.
            '''
            if btn1 == event.GetEventObject():
                ext = '.imctrl'
                title = 'Image control'
                print 'button 1'
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
        self.params['IMGfile'] = ''
        self.params['MaskFile'] = ''
        self.params['IgnoreMask'] = True
        fmtlist = G2IO.ExportPowderList(G2frame)
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnTimerLoop)

        controlsId = G2frame.PatternTree.GetSelection()
        size,imagefile = G2frame.PatternTree.GetItemPyData(G2frame.Image)
        self.imagedir,fileroot = os.path.split(imagefile)
        self.params['filter'] = '*'+os.path.splitext(fileroot)[1]
        os.chdir(self.imagedir)
        # get image names that have already been read
        self.IntegratedList = []
        for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
            self.IntegratedList.append(G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
                )[1])
            
        GSASIIpath.IPyBreak()
        
        wx.Frame.__init__(self, G2frame)
        mnpnl = wx.Panel(self)
        mnsizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Integration based on: '))
        self.ControlBaseLbl = wx.StaticText(mnpnl, wx.ID_ANY,'?')
        self.ControlBaseLbl.SetLabel(G2frame.PatternTree.GetItemText(G2frame.Image))
        sizer.Add(self.ControlBaseLbl)
        mnsizer.Add(sizer,1,wx.ALIGN_LEFT,1)
        mnpnl.SetSizer(mnsizer)
        #GSASIIpath.IPyBreak()
        # file filter stuff
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Image filter'))
        flterInp = G2G.ValidatedTxtCtrl(mnpnl,self.params,'filter')
        sizer.Add(flterInp)
        mnsizer.Add(sizer,1,wx.ALIGN_RIGHT,1)
        # box for integration controls & masks input
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Integration Controls/Masks from")
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        r1 = wx.RadioButton(mnpnl, wx.ID_ANY, "Use Active Image",
                            style = wx.RB_GROUP)
        r1.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
        lblsizr.Add(r1)
        r1.SetValue(True)
        r2 = wx.RadioButton(mnpnl, wx.ID_ANY, "Use from file(s)")
        lblsizr.Add(r2)
        r2.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
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
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Select output format(s)")
        lblsizr = wx.StaticBoxSizer(lbl, wx.HORIZONTAL)
        for dfmt in fmtlist:
            fmt = dfmt[1:]
            self.params[fmt] = False
            btn = G2G.G2CheckBox(mnpnl,dfmt,self.params,fmt)
            lblsizr.Add(btn)
        mnsizer.Add(lblsizr,1,wx.ALIGN_CENTER,1)

        # buttons on bottom
        mnsizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'AutoIntegration controls'))
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((20,-1))
        btnstart = wx.Button(mnpnl,  wx.ID_ANY, "Start")
        btnstart.Bind(wx.EVT_BUTTON, OnStart)
        sizer.Add(btnstart)
        btnstop = wx.Button(mnpnl,  wx.ID_ANY, "Stop")
        btnstop.Bind(wx.EVT_BUTTON, OnStop)
        sizer.Add(btnstop)
        sizer.Add((20,-1),wx.EXPAND,1)
        btnquit = wx.Button(mnpnl,  wx.ID_ANY, "Close")
        btnquit.Bind(wx.EVT_BUTTON, OnQuit)
        sizer.Add(btnquit)
        sizer.Add((20,-1))
        mnsizer.Add(sizer,1,wx.EXPAND|wx.BOTTOM,1)
        
        # finish up window
        OnRadioSelect(None) # disable widgets
        mnsizer.Fit(self)
        self.Show()
    

if __name__ == '__main__':
    app = wx.PySimpleApp()
    G2frame = wx.Frame(None) # create a top-level frame as a stand-in for the GSAS-II data tree
    G2frame.Show() 
    frm = AutoIntFrame(G2frame) # test the one above
    app.MainLoop()
 
