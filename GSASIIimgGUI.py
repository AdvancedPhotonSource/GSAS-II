#GSASII - image data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import wx
import wx.grid as wg
import matplotlib as mpl
import math
import time
import cPickle
import GSASIIpath
import GSASIIimage as G2img
import GSASIIplot as G2plt
import GSASIIIO as G2IO
import GSASIIgrid as G2gd

VERY_LIGHT_GREY = wx.Colour(235,235,235)

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi

                    
def UpdateImageControls(self,data,masks):
    import ImageCalibrants as calFile
    
    def OnDataType(event):
        data['type'] = typeSel.GetValue()[:4]
    
    def OnNewColorBar(event):
        data['color'] = colSel.GetValue()
        G2plt.PlotExposedImage(self,event=event)
        
    def OnNewCalibrant(event):
        data['calibrant'] = calSel.GetValue()
        data['calibskip'] = calFile.Calibrants[data['calibrant']][2]
        limits = calFile.Calibrants[data['calibrant']][3]
        data['calibdmin'],data['pixLimit'],data['cutoff'] = limits
        pixLimit.SetValue(str(limits[1]))
        cutOff.SetValue('%.1f'%(limits[2]))
        calibSkip.SetValue(str(data['calibskip']))
        calibDmin.SetValue('%.1f'%(limits[0]))
        
    def OnPixLimit(event):
        data['pixLimit'] = int(pixLimit.GetValue())
        
    def OnCalibSkip(event):
        data['calibskip'] = int(calibSkip.GetValue())
        
    def OnSetRings(event):
        if data['setRings']:
            data['setRings'] = False
        else:
            data['setRings'] = True
        G2plt.PlotExposedImage(self,event=event)

    def OnCalibDmin(event):
        try:
            dmin = float(calibDmin.GetValue())
            if dmin < 0.5:
                raise ValueError
            data['calibdmin'] = dmin
        except ValueError:
            pass
        calibDmin.SetValue("%.2f"%(data['calibdmin']))          #reset in case of error  
                
    def OnAzmthOff(event):
        try:
            azmthoff = float(azmthOff.GetValue())
            data['azmthOff'] = azmthoff
        except ValueError:
            pass
        azmthOff.SetValue("%.2f"%(data['azmthOff']))          #reset in case of error  
        
    def OnCutOff(event):
        try:
            cutoff = float(cutOff.GetValue())
            if cutoff < 0.1:
                raise ValueError
            data['cutoff'] = cutoff
        except ValueError:
            pass
        cutOff.SetValue("%.1f"%(data['cutoff']))          #reset in case of error  
        
    def OnMaxVal(event):
        try:
            value = min(data['range'][0][1],int(maxVal.GetValue()))
            if value < data['range'][1][0]+1:
                raise ValueError
            data['range'][1][1] = value
        except ValueError:
            pass
        maxVal.SetValue('%.0f'%(data['range'][1][1]))
        DeltOne = data['range'][1][1]-max(0.0,data['range'][0][0])
        sqrtDeltOne = math.sqrt(DeltOne)
        maxSel.SetValue(int(100*sqrtDeltOne/sqrtDeltZero))
        minSel.SetValue(int(100*(data['range'][1][0]/DeltOne)))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnMinVal(event):
        try:
            value = int(minVal.GetValue())
            if value > data['range'][1][1]-1:
                raise ValueError
            data['range'][1][0] = value
        except ValueError:
            pass
        minVal.SetValue('%.0f'%(data['range'][1][0]))
        minSel.SetValue(int(100*(data['range'][1][0]-max(0.0,data['range'][0][0]))/DeltOne))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnMaxSlider(event):
        sqrtDeltZero = math.sqrt(data['range'][0][1])
        imax = int(maxSel.GetValue())*sqrtDeltZero/100.
        data['range'][1][1] = imax**2
        data['range'][1][0] = max(0.0,min(data['range'][1][1]-1,data['range'][1][0]))
        DeltOne = max(1.0,data['range'][1][1]-data['range'][1][0])
        minSel.SetValue(int(100*(data['range'][1][0]/DeltOne)))
        maxVal.SetValue('%.0f'%(data['range'][1][1]))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnMinSlider(event):
        DeltOne = data['range'][1][1]-data['range'][1][0]
        imin = int(minSel.GetValue())*DeltOne/100.
        data['range'][1][0] = max(0.0,min(data['range'][1][1]-1,imin))
        minVal.SetValue('%.0f'%(data['range'][1][0]))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnNumOutChans(event):
        try:
            numChans = int(outChan.GetValue())
            if numChans < 1:
                raise ValueError
            data['outChannels'] = numChans
        except ValueError:
            pass
        outChan.SetValue(str(data['outChannels']))          #reset in case of error        
        
    def OnNumOutAzms(event):
        try:
            numAzms = int(outAzim.GetValue())
            if numAzms < 1:
                raise ValueError
            data['outAzimuths'] = numAzms            
        except ValueError:
            pass
        outAzim.SetValue(str(data['outAzimuths']))          #reset in case of error        
        G2plt.PlotExposedImage(self,event=event)
        
    def OnWavelength(event):
        try:
            wave = float(waveSel.GetValue())
            if wave < .01:
                raise ValueError
            data['wavelength'] = wave
        except ValueError:
            pass
        waveSel.SetValue("%6.5f" % (data['wavelength']))          #reset in case of error          
        
    def OnShowLines(event):
        if data['showLines']:
            data['showLines'] = False
        else:
            data['showLines'] = True
        G2plt.PlotExposedImage(self,event=event)
        
    def OnFullIntegrate(event):
        if data['fullIntegrate']:
            data['fullIntegrate'] = False
            self.Lazim.SetEditable(True)            
            self.Razim.SetEditable(True)            
        else:
            data['fullIntegrate'] = True
            self.Lazim.SetEditable(False)            
            self.Razim.SetEditable(False)
        UpdateImageControls(self,data,masks)            
        G2plt.PlotExposedImage(self,event=event)
        
    def OnSetDefault(event):
        import copy
        if data['setDefault']:
            self.imageDefault = {}
            data['setDefault'] = False
        else:
            self.imageDefault = copy.copy(data)
            data['setDefault'] = True
            
    def OnIOtth(event):
        Ltth = max(float(self.InnerTth.GetValue()),0.001)
        Utth = float(self.OuterTth.GetValue())
        if Ltth > Utth:
            Ltth,Utth = Utth,Ltth
        data['IOtth'] = [Ltth,Utth]
        self.InnerTth.SetValue("%8.3f" % (Ltth))
        self.OuterTth.SetValue("%8.2f" % (Utth))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnLRazim(event):
        Lazm = min(180,int(self.Lazim.GetValue()))
        Razm = max(-180,int(self.Razim.GetValue()))
        data['LRazimuth'] = [Lazm,Razm]
        G2plt.PlotExposedImage(self,event=event)
            
    def OnClearCalib(event):
        data['ring'] = []
        data['rings'] = []
        data['ellipses'] = []
        self.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMCLEARCALIB,enable=False)    
        G2plt.PlotExposedImage(self,event=event)
            
    def OnCalibrate(event):        
        self.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMCLEARCALIB,enable=True)    
        self.dataFrame.GetStatusBar().SetStatusText('Select > 4 points on 1st used ring; LB to pick, RB on point to delete else RB to finish')
        self.ifGetRing = True
        
    def OnIntegrate(event):
        self.Integrate = G2img.ImageIntegrate(self.ImageZ,data,masks)
        G2plt.PlotIntegration(self,newPlot=True)
        G2IO.SaveIntegration(self,self.PickId,data)
        self.MakePDF.Enable(True)
        
    def OnIntegrateAll(event):
        print 'integrate all'
        TextList = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append([False,name,id])
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if not len(TextList):
                self.ErrorDialog('Nothing to integrate','There must some "IMG" patterns')
                return
            dlg = self.CopyDialog(self,'Image integration controls','Select images to integrate:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    for item in result:
                        ifintegrate,name,id = item
                        if ifintegrate:
                            id = G2gd.GetPatternTreeItemId(self, self.root, name)
                            Npix,imagefile = self.PatternTree.GetItemPyData(id)
                            image = G2IO.GetImageData(self,imagefile,True)
                            Id = G2gd.GetPatternTreeItemId(self,id, 'Image Controls')
                            Data = self.PatternTree.GetItemPyData(Id)
                            try:
                                Masks = self.PatternTree.GetItemPyData(
                                    G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))
                            except TypeError:       #missing Masks
                                Imin,Imax = Data['Range']
                                Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                                self.PatternTree.SetItemPyData(
                                    G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'),Masks)                                
                            self.Integrate = G2img.ImageIntegrate(image,Data,Masks)
#                            G2plt.PlotIntegration(self,newPlot=True,event=event)
                            G2IO.SaveIntegration(self,Id,Data)
            finally:
                dlg.Destroy()
        
    def OnCopyControls(event):
        import copy
        TextList = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    if id == self.Image:
                        Source = name
                        Data = copy.deepcopy(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Image Controls')))
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['ellipses'] = []
                        Data['setDefault'] = False
                    else:
                        TextList.append([False,name,id])
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if not len(TextList):
                self.ErrorDialog('Nothing to copy controls to','There must be more than one "IMG" pattern')
                return
            dlg = self.CopyDialog(self,'Copy image controls','Copy controls from '+Source+' to:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    for i,item in enumerate(result):
                        ifcopy,name,id = item
                        if ifcopy:
                            oldData = copy.deepcopy(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Image Controls')))
                            Data['range'] = oldData['range']
                            Data['size'] = oldData['size']                                
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['ellipses'] = []
                            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Image Controls'),copy.deepcopy(Data))
            finally:
                dlg.Destroy()
                
    def OnSaveControls(event):
        dlg = wx.FileDialog(self, 'Choose image controls file', '.', '', 
            'image control files (*.imctrl)|*.imctrl',wx.OPEN)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'w')
                save = {}
                keys = ['type','wavelength','calibrant','distance','center','tilt','rotation','azmthOff']
                for key in keys:
                    File.write(key+':'+str(data[key])+'\n')
                File.close()
        finally:
            dlg.Destroy()
        
    def OnLoadControls(event):
        dlg = wx.FileDialog(self, 'Choose image controls file', '.', '', 
            'image control files (*.imctrl)|*.imctrl',wx.OPEN)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                save = {}
                S = File.readline()
                while S:
                    if S[0] == '#':
                        S = File.readline()
                        continue
                    [key,val] = S[:-1].split(':')
                    if key in ['type','calibrant']:
                        save[key] = val
                    elif key in ['wavelength','distance','tilt','rotation']:
                        save[key] = float(val)
                    elif key in ['center',]:
                        vals = val.strip('[] ').split()
                        save[key] = [float(vals[0]),float(vals[1])]                    
                    S = File.readline()
                data.update(save)
                calSel.SetValue(data['calibrant']) 
                waveSel.SetValue("%6.5f" % (data['wavelength']))
                cent = data['center']
                centText.SetValue(("%8.3f,%8.3f" % (cent[0],cent[1])))
                distSel.SetValue("%8.3f"%(data['distance']))
                tiltSel.SetValue("%9.3f"%(data['tilt']))            
                rotSel.SetValue("%9.3f"%(data['rotation']))
                File.close()
        finally:
            dlg.Destroy()
                                        
    #fix for old files:
    if 'azmthOff' not in data:
        data['azmthOff'] = 0.0
    #end fix
    
    colorList = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
    calList = [m for m in calFile.Calibrants.keys()]
    typeList = ['PWDR - powder diffraction data','SASD - small angle scattering data',
        'REFL - reflectometry data']
    if not data.get('type'):                        #patch for old project files
        data['type'] = 'PWDR'
    typeDict = {'PWDR':typeList[0],'SASD':typeList[1],'REFL':typeList[2]}
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.ImageMenu)
    if not self.dataFrame.GetStatusBar():
        self.dataFrame.CreateStatusBar()
    self.dataFrame.Bind(wx.EVT_MENU, OnCalibrate, id=G2gd.wxID_IMCALIBRATE)
    self.dataFrame.Bind(wx.EVT_MENU, OnClearCalib, id=G2gd.wxID_IMCLEARCALIB)
    if not data['rings']:
        self.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMCLEARCALIB,enable=False)    
    self.dataFrame.Bind(wx.EVT_MENU, OnIntegrate, id=G2gd.wxID_IMINTEGRATE)
    self.dataFrame.Bind(wx.EVT_MENU, OnIntegrateAll, id=G2gd.wxID_INTEGRATEALL)
    self.dataFrame.Bind(wx.EVT_MENU, OnCopyControls, id=G2gd.wxID_IMCOPYCONTROLS)
    self.dataFrame.Bind(wx.EVT_MENU, OnSaveControls, id=G2gd.wxID_IMSAVECONTROLS)
    self.dataFrame.Bind(wx.EVT_MENU, OnLoadControls, id=G2gd.wxID_IMLOADCONTROLS)
    self.dataDisplay = wx.Panel(self.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Type of image data: '),0,
        wx.ALIGN_CENTER_VERTICAL)
    typeSel = wx.ComboBox(parent=self.dataDisplay,value=typeDict[data['type']],choices=typeList,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    typeSel.SetValue(data['type'])
    typeSel.Bind(wx.EVT_COMBOBOX, OnDataType)
    comboSizer.Add(typeSel,0,wx.ALIGN_CENTER_VERTICAL)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Color bar '),0,
        wx.ALIGN_CENTER_VERTICAL)
    colSel = wx.ComboBox(parent=self.dataDisplay,value=data['color'],choices=colorList,
        style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
    colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
    comboSizer.Add(colSel,0,wx.ALIGN_CENTER_VERTICAL)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Azimuth offset '),0,
        wx.ALIGN_CENTER_VERTICAL)
    azmthOff = wx.TextCtrl(parent=self.dataDisplay,value=("%.2f" % (data['azmthOff'])),
        style=wx.TE_PROCESS_ENTER)
    azmthOff.Bind(wx.EVT_TEXT_ENTER,OnAzmthOff)
    azmthOff.Bind(wx.EVT_KILL_FOCUS,OnAzmthOff)
    comboSizer.Add(azmthOff,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(comboSizer,0,wx.ALIGN_LEFT)
    mainSizer.Add((5,5),0)
            
    maxSizer = wx.FlexGridSizer(2,3,0,5)
    maxSizer.AddGrowableCol(1,1)
    maxSizer.SetFlexibleDirection(wx.HORIZONTAL)
    sqrtDeltZero = math.sqrt(data['range'][0][1]-max(0.0,data['range'][0][0]))
    DeltOne = data['range'][1][1]-max(0.0,data['range'][0][0])
    sqrtDeltOne = math.sqrt(DeltOne)
    maxSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Max intensity'),0,
        wx.ALIGN_CENTER_VERTICAL)
    maxSel = wx.Slider(parent=self.dataDisplay,style=wx.SL_HORIZONTAL,
        value=int(100*sqrtDeltOne/sqrtDeltZero))
    maxSizer.Add(maxSel,1,wx.EXPAND)
    maxSel.Bind(wx.EVT_SLIDER, OnMaxSlider)
    maxVal = wx.TextCtrl(parent=self.dataDisplay,value='%.0f'%(data['range'][1][1]))
    maxVal.Bind(wx.EVT_TEXT_ENTER,OnMaxVal)    
    maxVal.Bind(wx.EVT_KILL_FOCUS,OnMaxVal)
    maxSizer.Add(maxVal,0,wx.ALIGN_CENTER_VERTICAL)    
    maxSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Min intensity'),0,
        wx.ALIGN_CENTER_VERTICAL)
    minSel = wx.Slider(parent=self.dataDisplay,style=wx.SL_HORIZONTAL,
        value=int(100*(data['range'][1][0]-max(0.0,data['range'][0][0]))/DeltOne))
    maxSizer.Add(minSel,1,wx.EXPAND)
    minSel.Bind(wx.EVT_SLIDER, OnMinSlider)
    minVal = wx.TextCtrl(parent=self.dataDisplay,value='%.0f'%(data['range'][1][0]))
    minVal.Bind(wx.EVT_TEXT_ENTER,OnMinVal)    
    minVal.Bind(wx.EVT_KILL_FOCUS,OnMinVal)
    maxSizer.Add(minVal,0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(maxSizer,0,wx.ALIGN_LEFT|wx.EXPAND)
    
    dataSizer = wx.FlexGridSizer(6,4,5,5)
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calibration coefficients'),0,
        wx.ALIGN_CENTER_VERTICAL)    
    dataSizer.Add((5,0),0)
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Integration coefficients'),0,
        wx.ALIGN_CENTER_VERTICAL)    
    dataSizer.Add((5,0),0)
    
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Beam center X,Y'),0,
        wx.ALIGN_CENTER_VERTICAL)
    cent = data['center']
    centText = wx.TextCtrl(parent=self.dataDisplay,value=("%8.3f,%8.3f" % (cent[0],cent[1])),style=wx.TE_READONLY)
    centText.SetBackgroundColour(VERY_LIGHT_GREY)
    dataSizer.Add(centText,0,wx.ALIGN_CENTER_VERTICAL)
    
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Inner/Outer 2-theta'),0,
        wx.ALIGN_CENTER_VERTICAL)
        
    IOtth = data['IOtth']
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    self.InnerTth = wx.TextCtrl(parent=self.dataDisplay,
        value=("%8.3f" % (IOtth[0])),style=wx.TE_PROCESS_ENTER)
    self.InnerTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
    self.InnerTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
    littleSizer.Add(self.InnerTth,0,wx.ALIGN_CENTER_VERTICAL)
    self.OuterTth = wx.TextCtrl(parent=self.dataDisplay,
        value=("%8.2f" % (IOtth[1])),style=wx.TE_PROCESS_ENTER)
    self.OuterTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
    self.OuterTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
    littleSizer.Add(self.OuterTth,0,wx.ALIGN_CENTER_VERTICAL)
    dataSizer.Add(littleSizer,0,)
       
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Wavelength'),0,
        wx.ALIGN_CENTER_VERTICAL)
    waveSel = wx.TextCtrl(parent=self.dataDisplay,value=("%6.5f" % (data['wavelength'])),
        style=wx.TE_PROCESS_ENTER)
    waveSel.Bind(wx.EVT_TEXT_ENTER,OnWavelength)
    waveSel.Bind(wx.EVT_KILL_FOCUS,OnWavelength)
    dataSizer.Add(waveSel,0,wx.ALIGN_CENTER_VERTICAL)
         
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Start/End azimuth'),0,
        wx.ALIGN_CENTER_VERTICAL)
    LRazim = data['LRazimuth']
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    self.Lazim = wx.TextCtrl(parent=self.dataDisplay,
        value=("%6d" % (LRazim[0])),style=wx.TE_PROCESS_ENTER)
    self.Lazim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
    self.Lazim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
    littleSizer.Add(self.Lazim,0,wx.ALIGN_CENTER_VERTICAL)
    self.Razim = wx.TextCtrl(parent=self.dataDisplay,
        value=("%6d" % (LRazim[1])),style=wx.TE_PROCESS_ENTER)
    self.Razim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
    self.Razim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
    littleSizer.Add(self.Razim,0,wx.ALIGN_CENTER_VERTICAL)
    dataSizer.Add(littleSizer,0,)
       
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Distance'),0,
        wx.ALIGN_CENTER_VERTICAL)
    distSel = wx.TextCtrl(parent=self.dataDisplay,value=("%8.3f"%(data['distance'])),style=wx.TE_READONLY)
    distSel.SetBackgroundColour(VERY_LIGHT_GREY)
    dataSizer.Add(distSel,0,wx.ALIGN_CENTER_VERTICAL)

    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' No. 2-theta/azimuth bins'),0,
        wx.ALIGN_CENTER_VERTICAL)
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    outChan = wx.TextCtrl(parent=self.dataDisplay,value=str(data['outChannels']),style=wx.TE_PROCESS_ENTER)
    outChan.Bind(wx.EVT_TEXT_ENTER,OnNumOutChans)
    outChan.Bind(wx.EVT_KILL_FOCUS,OnNumOutChans)
    littleSizer.Add(outChan,0,wx.ALIGN_CENTER_VERTICAL)
    outAzim = wx.TextCtrl(parent=self.dataDisplay,value=str(data['outAzimuths']),style=wx.TE_PROCESS_ENTER)
    outAzim.Bind(wx.EVT_TEXT_ENTER,OnNumOutAzms)
    outAzim.Bind(wx.EVT_KILL_FOCUS,OnNumOutAzms)
    littleSizer.Add(outAzim,0,wx.ALIGN_CENTER_VERTICAL)
    dataSizer.Add(littleSizer,0,)

    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Tilt angle'),0,
        wx.ALIGN_CENTER_VERTICAL)
    tiltSel = wx.TextCtrl(parent=self.dataDisplay,value=("%9.3f"%(data['tilt'])),style=wx.TE_READONLY)
    tiltSel.SetBackgroundColour(VERY_LIGHT_GREY)
    dataSizer.Add(tiltSel,0,wx.ALIGN_CENTER_VERTICAL)
    showLines = wx.CheckBox(parent=self.dataDisplay,label='Show integration limits?')
    dataSizer.Add(showLines,0,wx.ALIGN_CENTER_VERTICAL)
    showLines.Bind(wx.EVT_CHECKBOX, OnShowLines)
    showLines.SetValue(data['showLines'])
    fullIntegrate = wx.CheckBox(parent=self.dataDisplay,label='Do full integration?')
    dataSizer.Add(fullIntegrate,0,wx.ALIGN_CENTER_VERTICAL)
    fullIntegrate.Bind(wx.EVT_CHECKBOX, OnFullIntegrate)
    fullIntegrate.SetValue(data['fullIntegrate'])
    
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Tilt rotation'),0,
        wx.ALIGN_CENTER_VERTICAL)
    rotSel = wx.TextCtrl(parent=self.dataDisplay,value=("%9.3f"%(data['rotation']-90.)),style=wx.TE_READONLY)
    rotSel.SetBackgroundColour(VERY_LIGHT_GREY)
    dataSizer.Add(rotSel,0,wx.ALIGN_CENTER_VERTICAL)
    setDefault = wx.CheckBox(parent=self.dataDisplay,label='Use as default for all images?')
    dataSizer.Add(setDefault,0,wx.ALIGN_CENTER_VERTICAL)
    setDefault.Bind(wx.EVT_CHECKBOX, OnSetDefault)
    setDefault.SetValue(data['setDefault'])

    mainSizer.Add((5,5),0)        
    mainSizer.Add(dataSizer,0)
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calibration controls:'),0,
        wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    calibSizer = wx.FlexGridSizer(2,3,5,5)
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calibrant '),0,
        wx.ALIGN_CENTER_VERTICAL)
    calSel = wx.ComboBox(parent=self.dataDisplay,value=data['calibrant'],choices=calList,
        style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
    calSel.Bind(wx.EVT_COMBOBOX, OnNewCalibrant)
    comboSizer.Add(calSel,0,wx.ALIGN_CENTER_VERTICAL)
    calibSizer.Add(comboSizer,0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calib lines to skip   '),0,
        wx.ALIGN_CENTER_VERTICAL)
    calibSkip  = wx.ComboBox(parent=self.dataDisplay,value=str(data['calibskip']),choices=['0','1','2','3','4','5','6','7','8','9','10'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    calibSkip.Bind(wx.EVT_COMBOBOX, OnCalibSkip)
    comboSizer.Add(calibSkip,0,wx.ALIGN_CENTER_VERTICAL)
    calibSizer.Add(comboSizer,0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)        
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Min calib d-spacing '),0,
        wx.ALIGN_CENTER_VERTICAL)
    calibDmin = wx.TextCtrl(parent=self.dataDisplay,value=("%.2f" % (data['calibdmin'])),
        style=wx.TE_PROCESS_ENTER)
    calibDmin.Bind(wx.EVT_TEXT_ENTER,OnCalibDmin)
    calibDmin.Bind(wx.EVT_KILL_FOCUS,OnCalibDmin)
    comboSizer.Add(calibDmin,0,wx.ALIGN_CENTER_VERTICAL)
    calibSizer.Add(comboSizer,0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Min ring I/Ib '),0,
        wx.ALIGN_CENTER_VERTICAL)
    cutOff = wx.TextCtrl(parent=self.dataDisplay,value=("%.1f" % (data['cutoff'])),
        style=wx.TE_PROCESS_ENTER)
    cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)
    cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
    comboSizer.Add(cutOff,0,wx.ALIGN_CENTER_VERTICAL)
    calibSizer.Add(comboSizer,0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Pixel search range '),0,
        wx.ALIGN_CENTER_VERTICAL)
    pixLimit = wx.ComboBox(parent=self.dataDisplay,value=str(data['pixLimit']),choices=['1','2','5','10','15','20'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    pixLimit.Bind(wx.EVT_COMBOBOX, OnPixLimit)
    comboSizer.Add(pixLimit,0,wx.ALIGN_CENTER_VERTICAL)
    calibSizer.Add(comboSizer,0)
    
    comboSizer = wx.BoxSizer(wx.HORIZONTAL)
    setRings = wx.CheckBox(parent=self.dataDisplay,label='Show ring picks?')
    comboSizer.Add(setRings,0)
    setRings.Bind(wx.EVT_CHECKBOX, OnSetRings)
    setRings.SetValue(data['setRings'])
    calibSizer.Add(comboSizer,0)
    
    mainSizer.Add(calibSizer,0,wx.ALIGN_CENTER_VERTICAL)
        
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    fitSize = mainSizer.Fit(self.dataFrame)
    self.dataFrame.setSizePosLeft(fitSize)
    self.dataDisplay.SetSize(fitSize)
    
def UpdateMasks(self,data):
    
    def OnTextMsg(event):
        Obj = event.GetEventObject()
        Obj.SetToolTipString('Drag this mask on 2D Powder Image with mouse to change ')
        
    def OnThreshold(event):
        try:
            lower = max(int(lowerThreshold.GetValue()),thresh[0][0])
        except ValueError:
            lower = thresh[0][0]
        try:
            upper = min(int(upperThreshold.GetValue()),thresh[0][1])
        except ValueError:
            upper = thresh[0][1]
        data['Thresholds'][1] = [lower,upper]
        lowerThreshold.SetValue("%8d" % (lower))
        upperThreshold.SetValue("%8d" % (upper))
        G2plt.PlotExposedImage(self,event=event)
        
    def OnSpotDiameter(event):
        Obj = event.GetEventObject()
        try:
            diameter = min(100.,max(0.1,float(Obj.GetValue())))
        except ValueError:
            diameter = 1.0
        Obj.SetValue("%.2f"%(diameter))
        data['Points'][spotIds.index(Obj.GetId())][2] = diameter
        G2plt.PlotExposedImage(self,event=event)
        
    def OnDeleteSpot(event):
        Obj = event.GetEventObject()
        del(data['Points'][delSpotId.index(Obj)])
        UpdateMasks(self,data)           
        G2plt.PlotExposedImage(self,event=event)
        
    def OnRingThickness(event):
        Obj = event.GetEventObject()
        try:
            thick = min(1.0,max(0.001,float(Obj.GetValue())))
        except ValueError:
            thick = 0.1
        Obj.SetValue("%.3f"%(thick))
        data['Rings'][ringIds.index(Obj.GetId())][1] = thick
        G2plt.PlotExposedImage(self,event=event)
        
    def OnDeleteRing(event):
        Obj = event.GetEventObject()
        del(data['Rings'][delRingId.index(Obj)])
        UpdateMasks(self,data)           
        G2plt.PlotExposedImage(self,event=event)

    def OnArcThickness(event):
        Obj = event.GetEventObject()
        try:
            thick = min(20.0,max(0.001,float(Obj.GetValue())))
        except ValueError:
            thick = 0.1
        Obj.SetValue("%.3f"%(thick))
        data['Arcs'][arcIds.index(Obj.GetId())][2] = thick
        G2plt.PlotExposedImage(self,event=event)
        
    def OnDeleteArc(event):
        Obj = event.GetEventObject()
        del(data['Arcs'][delArcId.index(Obj)])
        UpdateMasks(self,data)           
        G2plt.PlotExposedImage(self,event=event)

    def OnDeletePoly(event):
        Obj = event.GetEventObject()
        del(data['Polygons'][delPolyId.index(Obj)])
        UpdateMasks(self,data)           
        G2plt.PlotExposedImage(self,event=event)

    def OnCopyMask(event):
        import copy
        TextList = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    if id == self.Image:
                        Source = name
                        Mask = copy.copy(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Masks')))
                        del Mask['Thresholds']
                    else:
                        TextList.append([False,name,id])
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if not len(TextList):
                self.ErrorDialog('Nothing to copy mask to','There must be more than one "IMG" pattern')
                return
            dlg = self.CopyDialog(self,'Copy mask information','Copy mask from '+Source+' to:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    for i,item in enumerate(result):
                        ifcopy,name,id = item
                        if ifcopy:
                            mask = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Masks'))
                            mask.update(Mask)                                
                            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'Masks'),mask)
            finally:
                dlg.Destroy()
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.MaskMenu)
    self.dataFrame.Bind(wx.EVT_MENU, OnCopyMask, id=G2gd.wxID_MASKCOPY)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
        Status.SetStatusText("To add mask: On 2D Powder Image, key a:arc, r:ring, s:spot, p:polygon")
    self.dataDisplay = wx.Panel(self.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)

    thresh = data['Thresholds']         #min/max intensity range
    spots = data['Points']               #x,y,radius in mm
    rings = data['Rings']               #radius, thickness
    polygons = data['Polygons']         #3+ x,y pairs
    arcs = data['Arcs']                 #radius, start/end azimuth, thickness
    
    littleSizer = wx.FlexGridSizer(2,3,0,5)
    littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Lower/Upper limits '),0,
        wx.ALIGN_CENTER_VERTICAL)
    Text = wx.TextCtrl(self.dataDisplay,value=("%8d" % (thresh[0][0])),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,wx.ALIGN_CENTER_VERTICAL)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    Text = wx.TextCtrl(self.dataDisplay,value=("%8d" % (thresh[0][1])),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,wx.ALIGN_CENTER_VERTICAL)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Lower/Upper thresholds '),
        0,wx.ALIGN_CENTER_VERTICAL)
    lowerThreshold = wx.TextCtrl(parent=self.dataDisplay,
        value=("%8d" % (thresh[1][0])),style=wx.TE_PROCESS_ENTER)
    lowerThreshold.Bind(wx.EVT_TEXT_ENTER,OnThreshold)
    lowerThreshold.Bind(wx.EVT_KILL_FOCUS,OnThreshold)
    littleSizer.Add(lowerThreshold,0,wx.ALIGN_CENTER_VERTICAL)
    upperThreshold = wx.TextCtrl(parent=self.dataDisplay,
        value=("%8d" % (thresh[1][1])),style=wx.TE_PROCESS_ENTER)
    upperThreshold.Bind(wx.EVT_TEXT_ENTER,OnThreshold)
    upperThreshold.Bind(wx.EVT_KILL_FOCUS,OnThreshold)
    littleSizer.Add(upperThreshold,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0,)
    spotIds = []
    delSpotId = []
    if spots:
        littleSizer = wx.FlexGridSizer(len(spots)+2,3,0,5)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Spot masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' position, mm'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' diameter, mm'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for x,y,d in spots:
            spotText = wx.TextCtrl(parent=self.dataDisplay,value=("%.2f,%.2f" % (x,y)),
                style=wx.TE_READONLY)
            spotText.SetBackgroundColour(VERY_LIGHT_GREY)
            littleSizer.Add(spotText,0,wx.ALIGN_CENTER_VERTICAL)
            spotText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
            spotDiameter = wx.TextCtrl(parent=self.dataDisplay,value=("%.2f" % (d)),
                style=wx.TE_PROCESS_ENTER)
            littleSizer.Add(spotDiameter,0,wx.ALIGN_CENTER_VERTICAL)
            spotDiameter.Bind(wx.EVT_TEXT_ENTER,OnSpotDiameter)
            spotDiameter.Bind(wx.EVT_KILL_FOCUS,OnSpotDiameter)
            spotIds.append(spotDiameter.GetId())
            spotDelete = wx.CheckBox(parent=self.dataDisplay,label='delete?')
            spotDelete.Bind(wx.EVT_CHECKBOX,OnDeleteSpot)
            delSpotId.append(spotDelete)
            littleSizer.Add(spotDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    ringIds = []
    delRingId = []
    if rings:
        littleSizer = wx.FlexGridSizer(len(rings)+2,3,0,5)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Ring masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' 2-theta,deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' thickness, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for tth,thick in rings:
            ringText = wx.TextCtrl(parent=self.dataDisplay,value=("%.3f" % (tth)),
                style=wx.TE_READONLY)
            ringText.SetBackgroundColour(VERY_LIGHT_GREY)
            ringText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
            littleSizer.Add(ringText,0,wx.ALIGN_CENTER_VERTICAL)
            ringThick = wx.TextCtrl(parent=self.dataDisplay,value=("%.3f" % (thick)),
                style=wx.TE_PROCESS_ENTER)
            littleSizer.Add(ringThick,0,wx.ALIGN_CENTER_VERTICAL)
            ringThick.Bind(wx.EVT_TEXT_ENTER,OnRingThickness)
            ringThick.Bind(wx.EVT_KILL_FOCUS,OnRingThickness)
            ringIds.append(ringThick.GetId())
            ringDelete = wx.CheckBox(parent=self.dataDisplay,label='delete?')
            ringDelete.Bind(wx.EVT_CHECKBOX,OnDeleteRing)
            delRingId.append(ringDelete)
            littleSizer.Add(ringDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    arcIds = []
    delArcId = []
    if arcs:
        littleSizer = wx.FlexGridSizer(len(rings)+2,4,0,5)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Arc masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' 2-theta,deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' azimuth, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' thickness, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for tth,azimuth,thick in arcs:
            arcText = wx.TextCtrl(parent=self.dataDisplay,value=("%.3f" % (tth)),
                style=wx.TE_READONLY)
            arcText.SetBackgroundColour(VERY_LIGHT_GREY)
            arcText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
            littleSizer.Add(arcText,0,wx.ALIGN_CENTER_VERTICAL)
            azmText = wx.TextCtrl(parent=self.dataDisplay,value=("%d,%d" % (azimuth[0],azimuth[1])),
                style=wx.TE_READONLY)
            azmText.SetBackgroundColour(VERY_LIGHT_GREY)
            azmText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
            littleSizer.Add(azmText,0,wx.ALIGN_CENTER_VERTICAL)
            arcThick = wx.TextCtrl(parent=self.dataDisplay,value=("%.3f" % (thick)),
                style=wx.TE_PROCESS_ENTER)
            littleSizer.Add(arcThick,0,wx.ALIGN_CENTER_VERTICAL)
            arcThick.Bind(wx.EVT_TEXT_ENTER,OnArcThickness)
            arcThick.Bind(wx.EVT_KILL_FOCUS,OnArcThickness)
            arcIds.append(arcThick.GetId())
            arcDelete = wx.CheckBox(parent=self.dataDisplay,label='delete?')
            arcDelete.Bind(wx.EVT_CHECKBOX,OnDeleteArc)
            delArcId.append(arcDelete)
            littleSizer.Add(arcDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    polyIds = []
    delPolyId = []
    if polygons:
        littleSizer = wx.FlexGridSizer(len(polygons)+2,2,0,5)
        littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Polygon masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for polygon in polygons:
            if polygon:
                polyList = []
                for x,y in polygon:
                    polyList.append("%.2f, %.2f"%(x,y))
                polyText = wx.ComboBox(self.dataDisplay,value=polyList[0],choices=polyList,style=wx.CB_READONLY)
                littleSizer.Add(polyText,0,wx.ALIGN_CENTER_VERTICAL)
                polyDelete = wx.CheckBox(parent=self.dataDisplay,label='delete?')
                polyDelete.Bind(wx.EVT_CHECKBOX,OnDeletePoly)
                delPolyId.append(polyDelete)
                littleSizer.Add(polyDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataDisplay.SetSize(mainSizer.Fit(self.dataFrame))
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))    
