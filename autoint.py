import os
import sys
import copy
import glob
import re
import bisect
import numpy as np
import wx
import wx.lib.mixins.listctrl  as  listmix
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIctrls as G2G
import GSASIIgrid as G2gd
import GSASIIimgGUI as G2imG
import GSASIIpy3 as G2py3

print 'loading autoint'

def ReadMask(filename):
    'Read a mask (.immask) file'
    File = open(filename,'r')
    save = {}
    S = File.readline()
    while S:
        if S[0] == '#':
            S = File.readline()
            continue
        [key,val] = S[:-1].split(':')
        if key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds']:
            save[key] = eval(val)
        S = File.readline()
    File.close()
    G2imG.CleanupMasks(save)
    return save

def ReadControls(filename):
    'read an image controls (.imctrl) file'
    cntlList = ['wavelength','distance','tilt','invert_x','invert_y','type',
            'fullIntegrate','outChannels','outAzimuths','LRazimuth','IOtth','azmthOff','DetDepth',
            'calibskip','pixLimit','cutoff','calibdmin','chisq','Flat Bkg',
            'PolaVal','SampleAbs','dark image','background image']
    File = open(filename,'r')
    save = {}
    S = File.readline()
    while S:
        if S[0] == '#':
            S = File.readline()
            continue
        [key,val] = S[:-1].split(':')
        if key in ['type','calibrant','binType','SampleShape',]:    #strings
            save[key] = val
        elif key in ['rotation']:
            save[key] = float(val)
        elif key in ['center',]:
            if ',' in val:
                save[key] = eval(val)
            else:
                vals = val.strip('[] ').split()
                save[key] = [float(vals[0]),float(vals[1])] 
        elif key in cntlList:
            save[key] = eval(val)
        S = File.readline()
    File.close()
    return save

def Read_imctrl(imctrl_file):
    '''Read an image control file and record control parms into a dict, with some simple
    type conversions
    '''
    file_opt = options = {}
    save = {'filename':imctrl_file}
    immask_file = os.path.splitext(imctrl_file)[0]+'.immask'
    if os.path.exists(immask_file):
        save['maskfile'] = immask_file
    else:
        save['maskfile'] = '(none)'
    cntlList = ['wavelength','distance','tilt','invert_x','invert_y','type',
                        'fullIntegrate','outChannels','outAzimuths','LRazimuth','IOtth','azmthOff','DetDepth',
                        'calibskip','pixLimit','cutoff','calibdmin','chisq','Flat Bkg',
                        'PolaVal','SampleAbs','dark image','background image']
    File = open(imctrl_file,'r')
    fullIntegrate = False
    try:
        S = File.readline()
        while S:
            if S[0] == '#':
                S = File.readline()
                continue
            [key,val] = S[:-1].split(':')
            if key in ['type','calibrant','binType','SampleShape',]:    #strings
                save[key] = val
            elif key == 'rotation':
                save[key] = float(val)
            elif key == 'fullIntegrate':
                fullIntegrate = eval(val)
            elif key == 'LRazimuth':
                save['LRazimuth_min'],save['LRazimuth_max'] = eval(val)[0:2]
            elif key == 'IOtth':
                save['IOtth_min'],save['IOtth_max'] = eval(val)[0:2]
            elif key == 'center':
                if ',' in val:
                    vals = eval(val)
                else:
                    vals = val.strip('[] ').split()
                    vals = [float(vals[0]),float(vals[1])] 
                save['center_x'],save['center_y'] = vals[0:2]
            elif key in cntlList:
                save[key] = eval(val)
            S = File.readline()
    finally:
        File.close()
        if fullIntegrate: save['LRazimuth_min'],save['LRazimuth_max'] = 0,0
    return save
    
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

        # Create a list of image files that have already been read
        imageFileList = []
        for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
            imgId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
            size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(imgId)
            if imagefile not in imageFileList: imageFileList.append(imagefile)
        # loop over image files matching glob, reading in new ones
        for newImage in self.currImageList:
            if newImage in imageFileList: continue # already read
            for imgId in G2IO.ReadImages(G2frame,newImage):
                controlsDict = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Image Controls'))
                ImageMasks = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Masks'))
                if self.params['Mode'] == 'table':
                    dist = controlsDict['distance']
                    interpDict,imgctrl,immask = self.Evaluator(dist) # interpolated calibration values
                    self.ImageControls = ReadControls(imgctrl)
                    self.ImageControls.update(interpDict)
                    self.ImageControls['showLines'] = True
                    self.ImageControls['ring'] = []
                    self.ImageControls['rings'] = []
                    self.ImageControls['ellipses'] = []
                    self.ImageControls['setDefault'] = False
                    for i in 'range','size','GonioAngles':
                        if i in self.ImageControls:
                            del self.ImageControls[i]
                    # load copy of Image Masks
                    if immask:
                        self.ImageMasks = ReadMask(immask)
                        del self.Thresholds['Thresholds']
                    else:
                        self.ImageMasks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[]}
                # update controls from master
                controlsDict.update(self.ImageControls)
                # update masks from master w/o Thresholds
                ImageMasks.update(self.ImageMasks)
        # now integrate the images that have not already been processed before
        for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
            if img in G2frame.IntegratedList: continue
            G2frame.IntegratedList.append(img)
            imgId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
            G2frame.Image = imgId
            G2frame.PickId = G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls')
            #  integrate in this entry
            size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(imgId)
            G2frame.ImageZ = G2IO.GetImageData(G2frame,imagefile,True,imagetag)
            masks = G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
            data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
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
                G2frame.AutointPWDRnames.append(treename)
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
        G2frame = self.G2frame
        # show current IMG base
        self.ControlBaseLbl.SetLabel(G2frame.PatternTree.GetItemText(G2frame.Image))
        if self.params['Mode'] != 'table':
            # load copy of Image Controls from current image and clean up
            # items that should not be copied
            self.ImageControls = copy.deepcopy(
                G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(
                    G2frame,G2frame.Image, 'Image Controls')))
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
                G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks')))
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
        if self.Reset: # special things to do after Reset has been pressed
            # reset controls and masks for all IMG items in tree to master
            for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
                # update controls from master
                controlsDict = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
                controlsDict.update(self.ImageControls)
                # update masks from master
                ImageMasks = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
                ImageMasks.update(self.ImageMasks)
            # delete all PWDR items created after last Start was pressed 
            idlist = []
            item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while item:
                itemName = G2frame.PatternTree.GetItemText(item)
                if itemName in G2frame.AutointPWDRnames:
                    idlist.append(item)
                item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            for item in idlist:
                G2frame.PatternTree.Delete(item)
            self.Reset = False
        G2frame.AutointPWDRnames = [] # list of created PWDR tree item names

    def __init__(self,G2frame,PollTime=60.0):
        def OnStart(event):
            '''Called when the start button is pressed. Changes button label 
            to Pause. When Pause is pressed the label changes to Resume.
            When either Start or Resume is pressed, the processing loop
            is started. When Pause is pressed, the loop is stopped. 
            '''
            # check inputs for errors before starting
            #err = ''
            #if not any([self.params[fmt] for fmt in self.fmtlist]):
            #    err += '\nPlease select at least one output format\n'
            #if err:
            #    G2G.G2MessageBox(self,err)
            #    return
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
                
        def OnRadioSelect(event):
            '''Respond to a radiobutton selection and when in table
            mode, get parameters from user. 
            '''
            self.Evaluator = None
            if r2.GetValue():
                self.params['Mode'] = 'table'
                try:
                    dlg = IntegParmTable(self.G2frame) # create the dialog
                    if dlg.ShowModal() == wx.ID_OK:
                        self.Evaluator = DefineEvaluator(dlg)
                    else:
                        r1.SetValue(True)
                finally:
                    dlg.Destroy()
            else:
                self.params['Mode'] = 'active'
        ##################################################
        # beginning of __init__ processing
        ##################################################
        self.G2frame = G2frame
        self.Evaluator = None
        self.params = {}
        self.Reset = False
        self.params['IMGfile'] = ''
        self.params['MaskFile'] = ''
        self.params['IgnoreMask'] = True
        self.fmtlist = G2IO.ExportPowderList(G2frame)
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnTimerLoop)

        controlsId = G2frame.PatternTree.GetSelection()
        size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(G2frame.Image)        
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
        r2 = wx.RadioButton(mnpnl, wx.ID_ANY, "Use from table")
        lblsizr.Add(r2)
        r2.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
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

def DefineEvaluator(dlg):
    '''Creates a function that provides interpolated values for a given distance value
    '''
    def Evaluator(dist):
        '''Interpolate image parameters for a supplied distance value

        :param float dist: distance to use for interpolation
        :returns: a list with 3 items:

          * a dict with parameter values,
          * the closest imctrl and
          * the closest maskfile (or None)
        '''            
        x = np.array([float(i) for i in parms[0]])
        closest = abs(x-dist).argmin()
        closeX = x[closest]
        D = {'distance':dist}
        imctfile = IMfileList[closest]
        if parms[-1][closest].lower() != '(none)':
            maskfile = parms[-1][closest]
        else:
            maskfile = None
        for c in range(1,cols-1):
            lbl = ParmList[c]
            if lbl in nonInterpVars:
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
            del D[a]
            del D[b]
        return D,imctfile,maskfile
    # save local copies of values needed in Evaluator
    parms = dlg.ReadImageParmTable()
    IMfileList = dlg.IMfileList
    cols = dlg.list.GetColumnCount()
    ParmList = dlg.ParmList
    nonInterpVars = dlg.nonInterpVars
    return Evaluator

class IntegParmTable(wx.Dialog):
    '''Creates a dialog window with a table of integration parameters.
    :meth:`ShowModal` will return wx.ID_OK if the process has been successful.
    In this case, :func:`DefineEvaluator` should be called to obtain a function that
    creates a dictionary with interpolated parameter values.
    '''
    ParmList = ('distance','center_x','center_y','wavelength','tilt','rotation','DetDepth',
            'LRazimuth_min','LRazimuth_max','IOtth_min','IOtth_max','outChannels',
            'maskfile',
            )
    nonInterpVars = ('tilt','rotation','LRazimuth_min','LRazimuth_max','IOtth_min','IOtth_max',
                     'outChannels')  # values in this list are taken from nearest rather than interpolated
    HeaderList = ('Det Dist','X cntr','Y cntr','wavelength','tilt','rotation','DetDepth',
            'Azimuth min','Azimuth max','2Th min','2Th max','Int. pts',
            'Mask File',
            )
    def __init__(self,G2frame):
        self.G2frame = G2frame
        self.parms = [] # list of values by column
        self.IMfileList = [] # list of .imctrl file names for each entry in table
        wx.Dialog.__init__(self,G2frame,style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE)
        files = []
        try:
            dlg = wx.FileDialog(self, 'Select image control files or previous table', 
                                style=wx.OPEN| wx.MULTIPLE,
                                wildcard='image control files (.imctrl)|*.imctrl|Integration table (*.imtbl)|*.imtbl')
            if dlg.ShowModal() == wx.ID_OK:
                files = dlg.GetPaths()
                self.parms,self.IMfileList = self.ReadFiles(files)
        finally:
            dlg.Destroy()
        if not files:
            wx.CallAfter(self.EndModal,wx.ID_CANCEL)
            return
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.list = ImgIntLstCtrl(self, wx.ID_ANY,
                      style=wx.LC_REPORT 
                          | wx.BORDER_SUNKEN
                         #| wx.BORDER_NONE
                         )
        mainSizer.Add(self.list,1,wx.EXPAND,1)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(self, wx.ID_OK)
        btnsizer.Add(btn)
        btn = wx.Button(self, wx.ID_ANY,'Save')
        btn.Bind(wx.EVT_BUTTON,self._onSave)
        btnsizer.Add(btn)
        btn = wx.Button(self, wx.ID_CLOSE,'Quit')
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)    
        self.SetSizer(mainSizer)
        self.list.FillList(self.parms)
        mainSizer.Layout()
        mainSizer.Fit(self)
        
    def ReadFiles(self,files):
        '''Reads a list of .imctrl files or a single .imtbl file
        '''
        tmpDict = {}
        if not files: return
        # option 1, a dump from a previous save
        if os.path.splitext(files[0])[1] == '.imtbl':
            fp = open(files[0],'r')
            S = fp.readline()
            while S:
                if S[0] != '#':
                    [key,val] = S[:-1].split(':')
                    tmpDict[key] = eval(val)
                S = fp.readline()
            fp.close()
            # delete entries
            m1 = [i for i,f in enumerate(tmpDict['filenames']) if not os.path.exists(f)]
            if m1:
                print('\nimctrl file not found:')
                for i in m1: print('\t#'+str(i)+': '+tmpDict['filenames'][i])
            m2 = [i for i,f in enumerate(tmpDict['maskfile']) if not (os.path.exists(f) or f.startswith('('))]
            if m2:
                print('\nmask file not found')
                for i in m2: print('\t#'+str(i)+': '+tmpDict['maskfile'][i])
            m3 = [i for i,d in enumerate(tmpDict['distance']) if d < 0]
            if m3:
                print('\nDropping entries due to negative distance: '+str(m3))
            m = sorted(set(m1 + m2 + m3))
            m.reverse()
            for c in m:
                for key in tmpDict:
                    del tmpDict[key][c]
            fileList = tmpDict.get('filenames','[]')
            parms = []
            for key in self.ParmList:
                try:
                    float(tmpDict[key][0])
                    parms.append([str(G2py3.FormatSigFigs(val,sigfigs=5)) for val in tmpDict[key]])
                except ValueError:
                    parms.append(tmpDict[key])
            return parms,fileList
        # option 2, read in a list of files
        for file in files: # read all files; place in dict by distance
            imgDict = Read_imctrl(file)
            tmpDict[imgDict.get('distance')] = imgDict
        parms = [[] for key in self.ParmList]
        fileList = []
        for d in sorted(tmpDict):
            fileList.append(tmpDict[d].get('filename'))
            if d is None: continue
            if d < 0: continue
            for i,key in enumerate(self.ParmList):
                val = tmpDict[d].get(key)
                try:
                    val = str(G2py3.FormatSigFigs(val,sigfigs=5))
                except:
                    val = str(val)
                parms[i].append(val)
        return parms,fileList
    
    def ReadImageParmTable(self):
        '''Reads possibly edited values from the ListCtrl table and returns a list
        of values for each column.
        '''
        rows = self.list.GetItemCount()
        cols = self.list.GetColumnCount()
        parms = []
        for c in range(cols):
            lbl = self.ParmList[c]
            parms.append([])
            for r in range(rows):
                parms[c].append(self.list.GetItem(r,c).GetText())
        return parms

    def _onClose(self,event):
        'Called when Cancel button is pressed'
        self.EndModal(wx.ID_CANCEL)
        
    def _onSave(self,event):
        'Called when save button is pressed; creates a .imtbl file'
        fil = ''
        if self.G2frame.GSASprojectfile:
            fil = os.path.splitext(self.G2frame.GSASprojectfile)[0]+'.imtbl'
        dir,f = os.path.split(fil)
        try:
            dlg = wx.FileDialog(self, 'Save table data as',
                        defaultDir=dir, defaultFile=f, style=wx.SAVE)
            if dlg.ShowModal() != wx.ID_OK: return
            fil = dlg.GetPath()
            fil = os.path.splitext(fil)[0]+'.imtbl'
        finally:
            dlg.Destroy()        
        parms = self.ReadImageParmTable()
        print('Writing image parameter table as '+fil)
        fp = open(fil,'w')
        for c in range(len(parms)-1):
            lbl = self.ParmList[c]
            fp.write(lbl+': '+str([eval(i) for i in parms[c]])+'\n')
        lbl = self.ParmList[c+1]
        fp.write(lbl+': '+str(parms[c+1])+'\n')
        lbl = 'filenames'
        fp.write(lbl+': '+str(self.IMfileList)+'\n')
        fp.close()
    
class ImgIntLstCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin,listmix.TextEditMixin):
    '''Creates a custom ListCtrl for editing Image Integration parameters
    '''
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=(1000,200), style=0):
        self.parent=parent
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.Bind(wx.EVT_LEFT_DCLICK, self.OnDouble)
        #self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick)
    def FillList(self,parms):
        'Places the current parms into the table'
        self.ClearAll()
        self.rowlen = len(self.parent.ParmList)
        for i,lbl in enumerate(self.parent.HeaderList):
            self.InsertColumn(i, lbl)
        for r,d in enumerate(parms[0]):
            if float(d) < 0: continue
            index = self.InsertStringItem(sys.maxint, d)
            for j in range(1,len(parms)):
                self.SetStringItem(index, j, parms[j][r])
        for i,lbl in enumerate(self.parent.ParmList):
            self.SetColumnWidth(i, wx.LIST_AUTOSIZE)

    def OnDouble(self,evt):
        'respond to a double-click'
        self.CloseEditor()
        fil = '(none)'
        try:
            dlg = wx.FileDialog(G2frame, 'Select mask or control file to add (Press cancel if none)', 
                                style=wx.OPEN,
                                wildcard='Add GSAS-II mask file (.immask)|*.immask|add image control file (.imctrl)|*.imctrl')
            if dlg.ShowModal() == wx.ID_OK:
                fil = dlg.GetPath()
        finally:
            dlg.Destroy()
        if os.path.splitext(fil)[1] != '.imctrl':
            self.SetStringItem(self.curRow, self.rowlen-1, fil)
            self.SetColumnWidth(self.rowlen-1, wx.LIST_AUTOSIZE)
        else:
            # insert or overwrite an instrument parameter set
            if not os.path.exists(fil):
                print('Does not exist: '+fil)
                return
            imgDict = Read_imctrl(fil)
            dist = imgDict['distance']
            parms = self.parent.ReadImageParmTable()
            x = np.array([float(i) for i in parms[0]])
            closest = abs(x-dist).argmin()
            closeX = x[closest]
            # fix IMfileList
            for c,lbl in enumerate(self.parent.ParmList):
                try:
                    vali = G2py3.FormatSigFigs(float(imgDict[lbl]),sigfigs=5)
                except ValueError:
                    vali = imgDict[lbl]
                if abs(closeX-dist) < 1.: # distance is within 1 mm, replace
                    parms[c][closest] = vali
                elif dist > closeX: # insert after
                    parms[c].insert(closest+1,vali)
                else:
                    parms[c].insert(closest,vali)
            if abs(closeX-dist) < 1.: # distance is within 1 mm, replace
                self.parent.IMfileList[closest] = fil
            elif dist > closeX: # insert after
                self.parent.IMfileList.insert(closest+1,fil)
            else:
                self.parent.IMfileList.insert(closest,fil)
            self.FillList(parms)

if __name__ == '__main__':
    app = wx.PySimpleApp()
    G2frame = wx.Frame(None) # create a top-level frame as a stand-in for the GSAS-II data tree
    G2frame.Show() 
    frm = AutoIntFrame(G2frame) # test the one above
    app.MainLoop()
 
