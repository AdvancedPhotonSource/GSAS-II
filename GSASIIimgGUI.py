# -*- coding: utf-8 -*-
#GSASII - image data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIimgGUI: Image GUI*
-------------------------

Control image display and processing

'''
import wx
import matplotlib as mpl
import math
import time
import copy
import cPickle
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIimage as G2img
import GSASIIplot as G2plt
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import numpy as np

VERY_LIGHT_GREY = wx.Colour(235,235,235)

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi

################################################################################
##### Image Controls
################################################################################                    
def UpdateImageControls(G2frame,data,masks):
    '''Shows and handles the controls on the "Image Controls"
    data tree entry
    '''
    import ImageCalibrants as calFile
#patch
    if 'GonioAngles' not in data:
        data['GonioAngles'] = [0.,0.,0.]
    if 'DetDepth' not in data:
        data['DetDepth'] = 0.
        data['DetDepthRef'] = False
#end patch

    
# Menu items
            
    def OnCalibrate(event):        
        G2frame.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMRECALIBRATE,enable=True)    
        G2frame.dataFrame.GetStatusBar().SetStatusText('Select > 4 points on 1st used ring; LB to pick, RB on point to delete else RB to finish')
        G2frame.ifGetRing = True
        
    def OnRecalibrate(event):
        G2img.ImageRecalibrate(G2frame,data)
        UpdateImageControls(G2frame,data,masks)
        
    def OnClearCalib(event):
        data['ring'] = []
        data['rings'] = []
        data['ellipses'] = []
#        G2frame.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMRECALIBRATE,enable=False)    
        G2plt.PlotExposedImage(G2frame,event=event)
            
    def OnIntegrate(event):
        dlg = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        try:
            if data['background image'][0]:
                maskCopy = copy.deepcopy(masks)
                backImg = data['background image'][0]
                backScale = data['background image'][1]
                id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, backImg)
                Npix,imagefile = G2frame.PatternTree.GetItemPyData(id)
                backImage = G2IO.GetImageData(G2frame,imagefile,True)*backScale
                sumImage = G2frame.ImageZ+backImage
                sumMin = np.min(sumImage)
                sumMax = np.max(sumImage)
                maskCopy['Thresholds'] = [(sumMin,sumMax),[sumMin,sumMax]]
                G2frame.Integrate = G2img.ImageIntegrate(sumImage,data,maskCopy,dlg)
            else:
                G2frame.Integrate = G2img.ImageIntegrate(G2frame.ImageZ,data,masks,dlg)
    #        G2plt.PlotIntegration(G2frame,newPlot=True)
            G2IO.SaveIntegration(G2frame,G2frame.PickId,data)
        finally:
            dlg.Destroy()
        for item in G2frame.MakePDF: item.Enable(True)
        
    def OnIntegrateAll(event):
        print 'integrate all'
        TextList = [[False,'All IMG',0]]
        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append([False,name,id])
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            if len(TextList) == 1:
                G2frame.ErrorDialog('Nothing to integrate','There must some "IMG" patterns')
                return
            dlg = G2frame.CopyDialog(G2frame,'Image integration controls','Select images to integrate:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    if result[0][0]:                    #the 'All IMG' is True
                        result = TextList[1:]
                        for item in result: item[0] = True
                    for item in result:
                        ifintegrate,name,id = item
                        if ifintegrate:
                            dlgp = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
                                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
                            try:
                                id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, name)
                                Npix,imagefile = G2frame.PatternTree.GetItemPyData(id)
                                image = G2IO.GetImageData(G2frame,imagefile,True)
                                Id = G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls')
                                Data = G2frame.PatternTree.GetItemPyData(Id)
                                backImage = []
                                if Data['background image'][0]:
                                    backImg = Data['background image'][0]
                                    backScale = Data['background image'][1]
                                    id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, backImg)
                                    Npix,imagefile = G2frame.PatternTree.GetItemPyData(id)
                                    backImage = G2IO.GetImageData(G2frame,imagefile,True)*backScale
                                try:
                                    Masks = G2frame.PatternTree.GetItemPyData(
                                        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
                                except TypeError:       #missing Masks
                                    Imin,Imax = Data['Range']
                                    Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                                    G2frame.PatternTree.SetItemPyData(
                                        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'),Masks)
                                if len(backImage):                                
                                    G2frame.Integrate = G2img.ImageIntegrate(image+backImage,Data,Masks,dlgp)
                                else:
                                    G2frame.Integrate = G2img.ImageIntegrate(image,Data,Masks,dlgp)
#                               G2plt.PlotIntegration(G2frame,newPlot=True,event=event)
                                G2IO.SaveIntegration(G2frame,Id,Data)
                            finally:
                                dlgp.Destroy()
            finally:
                dlg.Destroy()
        
    def OnCopyControls(event):
        import copy
        TextList = [[False,'All IMG',0]]
        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    if id == G2frame.Image:
                        Source = name
                        Data = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls')))
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['ellipses'] = []
                        Data['setDefault'] = False
                    else:
                        TextList.append([False,name,id])
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            if len(TextList) == 1:
                G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "IMG" pattern')
                return
            dlg = G2frame.CopyDialog(G2frame,'Copy image controls','Copy controls from '+Source+' to:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    if result[0][0]:
                        result = TextList[1:]
                        for item in result: item[0] = True
                    for i,item in enumerate(result):
                        ifcopy,name,id = item
                        if ifcopy:
                            oldData = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls')))
                            Data['range'] = oldData['range']
                            Data['size'] = oldData['size']
                            Data['GonioAngles'] = oldData.get('GonioAngles', [0.,0.,0.])
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['ellipses'] = []
                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls'),copy.deepcopy(Data))
            finally:
                dlg.Destroy()
                
    def OnSaveControls(event):
        dlg = wx.FileDialog(G2frame, 'Choose image controls file', '.', '', 
            'image control files (*.imctrl)|*.imctrl',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'w')
                save = {}
                keys = ['type','wavelength','calibrant','distance','center',
                    'tilt','rotation','azmthOff','fullIntegrate','LRazimuth',
                    'IOtth','outAzimuths','invert_x','invert_y']
                for key in keys:
                    if key in ['rotation']:
                        File.write(key+':'+str(data[key])+'\n')                        
                    else:
                        File.write(key+':'+str(data[key])+'\n')
                File.close()
        finally:
            dlg.Destroy()
        
    def OnLoadControls(event):
        cntlList = ['wavelength','distance','tilt','invert_x','invert_y',
            'fullIntegrate','outAzimuths','LRazimuth','IOtth','azmthOff']
        dlg = wx.FileDialog(G2frame, 'Choose image controls file', '.', '', 
            'image control files (*.imctrl)|*.imctrl',wx.OPEN|wx.CHANGE_DIR)
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
                    if key in ['type','calibrant',]:
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
                data.update(save)
                UpdateImageControls(G2frame,data,masks)
                G2plt.PlotExposedImage(G2frame,event=event)
                
                File.close()
        finally:
            dlg.Destroy()
            
# Sizers
                                        
    def ComboSizer():

        def OnDataType(event):
            data['type'] = typeSel.GetValue()[:4]
            wx.CallAfter(UpdateImageControls,G2frame,data,masks)
    
        def OnNewColorBar(event):
            data['color'] = colSel.GetValue()
            G2plt.PlotExposedImage(G2frame,event=event)
        
        def OnAzmthOff(event):
            try:
                azmthoff = float(azmthOff.GetValue())
                data['azmthOff'] = azmthoff
            except ValueError:
                pass
            azmthOff.SetValue("%.2f"%(data['azmthOff']))          #reset in case of error  
            G2plt.PlotExposedImage(G2frame,event=event)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Type of image data: '),0,
            wx.ALIGN_CENTER_VERTICAL)
        typeSel = wx.ComboBox(parent=G2frame.dataDisplay,value=typeDict[data['type']],choices=typeList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        typeSel.SetValue(data['type'])
        typeSel.Bind(wx.EVT_COMBOBOX, OnDataType)
        comboSizer.Add(typeSel,0,wx.ALIGN_CENTER_VERTICAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Color bar '),0,
            wx.ALIGN_CENTER_VERTICAL)
        colSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['color'],choices=colorList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
        colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
        comboSizer.Add(colSel,0,wx.ALIGN_CENTER_VERTICAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Azimuth offset '),0,
            wx.ALIGN_CENTER_VERTICAL)
        azmthOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f" % (data['azmthOff'])),
            style=wx.TE_PROCESS_ENTER)
        azmthOff.Bind(wx.EVT_TEXT_ENTER,OnAzmthOff)
        azmthOff.Bind(wx.EVT_KILL_FOCUS,OnAzmthOff)
        comboSizer.Add(azmthOff,0,wx.ALIGN_CENTER_VERTICAL)
        return comboSizer
        
    def MaxSizer():
                
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
            G2plt.PlotExposedImage(G2frame,event=event)
            
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
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnMaxSlider(event):
            sqrtDeltZero = math.sqrt(data['range'][0][1])
            imax = int(maxSel.GetValue())*sqrtDeltZero/100.
            data['range'][1][1] = imax**2
            data['range'][1][0] = max(0.0,min(data['range'][1][1]-1,data['range'][1][0]))
            DeltOne = max(1.0,data['range'][1][1]-data['range'][1][0])
            minSel.SetValue(int(100*(data['range'][1][0]/DeltOne)))
            maxVal.SetValue('%.0f'%(data['range'][1][1]))
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnMinSlider(event):
            DeltOne = data['range'][1][1]-data['range'][1][0]
            imin = int(minSel.GetValue())*DeltOne/100.
            data['range'][1][0] = max(0.0,min(data['range'][1][1]-1,imin))
            minVal.SetValue('%.0f'%(data['range'][1][0]))
            G2plt.PlotExposedImage(G2frame,event=event)
            
        maxSizer = wx.FlexGridSizer(2,3,0,5)
        maxSizer.AddGrowableCol(1,1)
        maxSizer.SetFlexibleDirection(wx.HORIZONTAL)
        sqrtDeltZero = math.sqrt(data['range'][0][1]-max(0.0,data['range'][0][0]))
        DeltOne = data['range'][1][1]-max(0.0,data['range'][0][0])
        sqrtDeltOne = math.sqrt(DeltOne)
        maxSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Max intensity'),0,
            wx.ALIGN_CENTER_VERTICAL)
        maxSel = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(100*sqrtDeltOne/sqrtDeltZero))
        maxSizer.Add(maxSel,1,wx.EXPAND)
        maxSel.Bind(wx.EVT_SLIDER, OnMaxSlider)
        maxVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.0f'%(data['range'][1][1]))
        maxVal.Bind(wx.EVT_TEXT_ENTER,OnMaxVal)    
        maxVal.Bind(wx.EVT_KILL_FOCUS,OnMaxVal)
        maxSizer.Add(maxVal,0,wx.ALIGN_CENTER_VERTICAL)    
        maxSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min intensity'),0,
            wx.ALIGN_CENTER_VERTICAL)
        minSel = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(100*(data['range'][1][0]-max(0.0,data['range'][0][0]))/DeltOne))
        maxSizer.Add(minSel,1,wx.EXPAND)
        minSel.Bind(wx.EVT_SLIDER, OnMinSlider)
        minVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.0f'%(data['range'][1][0]))
        minVal.Bind(wx.EVT_TEXT_ENTER,OnMinVal)    
        minVal.Bind(wx.EVT_KILL_FOCUS,OnMinVal)
        maxSizer.Add(minVal,0,wx.ALIGN_CENTER_VERTICAL)
        return maxSizer
        
    def CalibCoeffSizer():
        
        def OnWavelength(event):
            try:
                wave = float(waveSel.GetValue())
                if wave < .01:
                    raise ValueError
                data['wavelength'] = wave
            except ValueError:
                pass
            waveSel.SetValue("%6.5f" % (data['wavelength']))          #reset in case of error
            
        def OnDetDepthRef(event):
            data['DetDepthRef'] = penSel.GetValue()
            
        def OnDetDepth(event):
            try:
                data['DetDepth'] = float(penVal.GetValue())
            except ValueError:
                pass
            penVal.SetValue("%6.3f" % (data['DetDepth']))          #reset in case of error                      
            
        calibSizer = wx.FlexGridSizer(5,2,5,5)
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibration coefficients'),0,
            wx.ALIGN_CENTER_VERTICAL)    
        calibSizer.Add((5,0),0)        
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Beam center X,Y'),0,
            wx.ALIGN_CENTER_VERTICAL)
        cent = data['center']
        centText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%8.3f,%8.3f" % (cent[0],cent[1])),style=wx.TE_READONLY)
        centText.SetBackgroundColour(VERY_LIGHT_GREY)
        calibSizer.Add(centText,0,wx.ALIGN_CENTER_VERTICAL)        
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Wavelength'),0,
            wx.ALIGN_CENTER_VERTICAL)
        waveSel = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%6.5f" % (data['wavelength'])),
            style=wx.TE_PROCESS_ENTER)
        waveSel.Bind(wx.EVT_TEXT_ENTER,OnWavelength)
        waveSel.Bind(wx.EVT_KILL_FOCUS,OnWavelength)
        calibSizer.Add(waveSel,0,wx.ALIGN_CENTER_VERTICAL)             
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Distance'),0,
            wx.ALIGN_CENTER_VERTICAL)
        distSel = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%8.3f"%(data['distance'])),style=wx.TE_READONLY)
        distSel.SetBackgroundColour(VERY_LIGHT_GREY)
        calibSizer.Add(distSel,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Tilt angle'),0,
            wx.ALIGN_CENTER_VERTICAL)
        tiltSel = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%9.3f"%(data['tilt'])),style=wx.TE_READONLY)
        tiltSel.SetBackgroundColour(VERY_LIGHT_GREY)
        calibSizer.Add(tiltSel,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Tilt rotation'),0,
            wx.ALIGN_CENTER_VERTICAL)
        rotSel = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%9.3f"%(data['rotation'])),style=wx.TE_READONLY)
        rotSel.SetBackgroundColour(VERY_LIGHT_GREY)
        calibSizer.Add(rotSel,0,wx.ALIGN_CENTER_VERTICAL)
        penSel = wx.CheckBox(parent=G2frame.dataDisplay,label='Penetration?')
        calibSizer.Add(penSel,0,wx.ALIGN_CENTER_VERTICAL)
        penSel.Bind(wx.EVT_CHECKBOX, OnDetDepthRef)
        penSel.SetValue(data['DetDepthRef'])
        penVal = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%6.5f" % (data['DetDepth'])),
            style=wx.TE_PROCESS_ENTER)
        penVal.Bind(wx.EVT_TEXT_ENTER,OnDetDepth)
        penVal.Bind(wx.EVT_KILL_FOCUS,OnDetDepth)
        calibSizer.Add(penVal,0,wx.ALIGN_CENTER_VERTICAL)             
        
        return calibSizer
    
    def IntegrateSizer():
        
        def OnIOtth(event):
            Ltth = max(float(G2frame.InnerTth.GetValue()),0.001)
            Utth = float(G2frame.OuterTth.GetValue())
            if Ltth > Utth:
                Ltth,Utth = Utth,Ltth
            data['IOtth'] = [Ltth,Utth]
            G2frame.InnerTth.SetValue("%8.3f" % (Ltth))
            G2frame.OuterTth.SetValue("%8.2f" % (Utth))
            G2plt.PlotExposedImage(G2frame,event=event)
        
        def OnLRazim(event):
            Lazm =int(G2frame.Lazim.GetValue())
            if data['fullIntegrate']:
               G2frame.Razim.SetValue("%6d" % (Lazm+360))
            Razm = int(G2frame.Razim.GetValue())
            if Lazm > Razm:
                Lazm -= 360
            data['LRazimuth'] = [Lazm,Razm]
            G2plt.PlotExposedImage(G2frame,event=event)
        
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
            G2plt.PlotExposedImage(G2frame,event=event)
        
        def OnOblique(event):
            if data['Oblique'][1]:
                data['Oblique'][1] = False
            else:
                data['Oblique'][1] = True
                
        def OnObliqVal(event):
            try:
                value = float(obliqVal.GetValue())
                if 0.01 <= value <= 0.99:
                    data['Oblique'][0] = value
                else:
                    raise ValueError
            except ValueError:
                pass
            obliqVal.SetValue('%.3f'%(data['Oblique'][0]))
                           
        def OnShowLines(event):
            if data['showLines']:
                data['showLines'] = False
            else:
                data['showLines'] = True
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnFullIntegrate(event):
            Lazm =int(G2frame.Lazim.GetValue())
            if data['fullIntegrate']:
                data['fullIntegrate'] = False
                data['LRazimuth'] = [Lazm,Lazm+20]
            else:
                data['fullIntegrate'] = True
                data['LRazimuth'] = [Lazm,Lazm+360]
            UpdateImageControls(G2frame,data,masks)
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnSetDefault(event):
            import copy
            if data['setDefault']:
                G2frame.imageDefault = {}
                data['setDefault'] = False
            else:
                G2frame.imageDefault = copy.copy(data)
                data['setDefault'] = True
                
        def OnCenterAzm(event):
            if data['centerAzm']:
                data['centerAzm'] = False
            else:
                data['centerAzm'] = True
            G2plt.PlotExposedImage(G2frame,event=event)
                
        def OnApplyPola(event):
            if data['PolaVal'][1]:
                data['PolaVal'][1] = False
            else:
                data['PolaVal'][1] = True
                
        def OnPolaVal(event):
            try:
                value = float(polaVal.GetValue())
                if 0.001 <= value <= 0.999:
                    data['PolaVal'][0] = value
                else:
                    raise ValueError
            except ValueError:
                pass
            polaVal.SetValue('%.3f'%(data['PolaVal'][0]))
                           
        dataSizer = wx.FlexGridSizer(5,2,5,3)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Integration coefficients'),0,
            wx.ALIGN_CENTER_VERTICAL)    
        dataSizer.Add((5,0),0)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Inner/Outer 2-theta'),0,
            wx.ALIGN_CENTER_VERTICAL)
            
        IOtth = data['IOtth']
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        G2frame.InnerTth = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%8.3f" % (IOtth[0])),style=wx.TE_PROCESS_ENTER)
        G2frame.InnerTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
        G2frame.InnerTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
        littleSizer.Add(G2frame.InnerTth,0,wx.ALIGN_CENTER_VERTICAL)
        G2frame.OuterTth = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%8.2f" % (IOtth[1])),style=wx.TE_PROCESS_ENTER)
        G2frame.OuterTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
        G2frame.OuterTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
        littleSizer.Add(G2frame.OuterTth,0,wx.ALIGN_CENTER_VERTICAL)
        dataSizer.Add(littleSizer,0,)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Start/End azimuth'),0,
            wx.ALIGN_CENTER_VERTICAL)
        LRazim = data['LRazimuth']
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        G2frame.Lazim = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%6d" % (LRazim[0])),style=wx.TE_PROCESS_ENTER)
        G2frame.Lazim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
        G2frame.Lazim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
        littleSizer.Add(G2frame.Lazim,0,wx.ALIGN_CENTER_VERTICAL)
        G2frame.Razim = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%6d" % (LRazim[1])),style=wx.TE_PROCESS_ENTER)
        G2frame.Razim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
        G2frame.Razim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
        if data['fullIntegrate']:
            G2frame.Razim.Enable(False)
            G2frame.Razim.SetBackgroundColour(VERY_LIGHT_GREY)
            G2frame.Razim.SetValue("%6d" % (LRazim[0]+360))
        littleSizer.Add(G2frame.Razim,0,wx.ALIGN_CENTER_VERTICAL)
        dataSizer.Add(littleSizer,0,)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' No. 2-theta/azimuth bins'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        outChan = wx.TextCtrl(parent=G2frame.dataDisplay,value=str(data['outChannels']),style=wx.TE_PROCESS_ENTER)
        outChan.Bind(wx.EVT_TEXT_ENTER,OnNumOutChans)
        outChan.Bind(wx.EVT_KILL_FOCUS,OnNumOutChans)
        littleSizer.Add(outChan,0,wx.ALIGN_CENTER_VERTICAL)
        outAzim = wx.TextCtrl(parent=G2frame.dataDisplay,value=str(data['outAzimuths']),style=wx.TE_PROCESS_ENTER)
        outAzim.Bind(wx.EVT_TEXT_ENTER,OnNumOutAzms)
        outAzim.Bind(wx.EVT_KILL_FOCUS,OnNumOutAzms)
        littleSizer.Add(outAzim,0,wx.ALIGN_CENTER_VERTICAL)
        dataSizer.Add(littleSizer,0,)
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        oblique = wx.CheckBox(parent=G2frame.dataDisplay,label='Apply detector absorption?')
        dataSizer.Add(oblique,0,wx.ALIGN_CENTER_VERTICAL)
        oblique.Bind(wx.EVT_CHECKBOX, OnOblique)
        oblique.SetValue(data['Oblique'][1])
        littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Value (0.01-0.99)  '),0,
            wx.ALIGN_CENTER_VERTICAL)
        obliqVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.3f'%(data['Oblique'][0]),style=wx.TE_PROCESS_ENTER)
        obliqVal.Bind(wx.EVT_TEXT_ENTER,OnObliqVal)
        obliqVal.Bind(wx.EVT_KILL_FOCUS,OnObliqVal)
        littleSizer.Add(obliqVal,0,wx.ALIGN_CENTER_VERTICAL)
        dataSizer.Add(littleSizer,0,)
        if 'SASD' in data['type']:
            littleSizer = wx.BoxSizer(wx.HORIZONTAL)
            setPolariz = wx.CheckBox(parent=G2frame.dataDisplay,label='Apply polarization?')
            dataSizer.Add(setPolariz,0,wx.ALIGN_CENTER_VERTICAL)
            setPolariz.Bind(wx.EVT_CHECKBOX, OnApplyPola)
            setPolariz.SetValue(data['PolaVal'][1])
            littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Value (0.001-0.999)  '),0,
                wx.ALIGN_CENTER_VERTICAL)
            polaVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.3f'%(data['PolaVal'][0]),
                style=wx.TE_PROCESS_ENTER)
            polaVal.Bind(wx.EVT_TEXT_ENTER,OnPolaVal)
            polaVal.Bind(wx.EVT_KILL_FOCUS,OnPolaVal)
            littleSizer.Add(polaVal,0,wx.ALIGN_CENTER_VERTICAL)
            dataSizer.Add(littleSizer,0,)
        
        showLines = wx.CheckBox(parent=G2frame.dataDisplay,label='Show integration limits?')
        dataSizer.Add(showLines,0,wx.ALIGN_CENTER_VERTICAL)
        showLines.Bind(wx.EVT_CHECKBOX, OnShowLines)
        showLines.SetValue(data['showLines'])
        fullIntegrate = wx.CheckBox(parent=G2frame.dataDisplay,label='Do full integration?')
        dataSizer.Add(fullIntegrate,0,wx.ALIGN_CENTER_VERTICAL)
        fullIntegrate.Bind(wx.EVT_CHECKBOX, OnFullIntegrate)
        fullIntegrate.SetValue(data['fullIntegrate'])
        setDefault = wx.CheckBox(parent=G2frame.dataDisplay,label='Use as default for all images?')
        dataSizer.Add(setDefault,0,wx.ALIGN_CENTER_VERTICAL)
        setDefault.Bind(wx.EVT_CHECKBOX, OnSetDefault)
        setDefault.SetValue(data['setDefault'])
        centerAzm = wx.CheckBox(parent=G2frame.dataDisplay,label='Azimuth at bin center?')
        dataSizer.Add(centerAzm,0,wx.ALIGN_CENTER_VERTICAL)
        centerAzm.Bind(wx.EVT_CHECKBOX, OnCenterAzm)
        centerAzm.SetValue(data['centerAzm'])
        return dataSizer
        
    def BackSizer():
        
        def OnBackImage(event):
            data['background image'][0] = backImage.GetValue()
            
        def OnBackMult(event):
            try:
                mult = float(backMult.GetValue())
                data['background image'][1] = mult
            except ValueError:
                pass
            backMult.SetValue("%.3f" % (data['background image'][1]))          #reset in case of error 
        
        backSizer = wx.FlexGridSizer(1,4,5,5)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background image'),0,wx.ALIGN_CENTER_VERTICAL)
        Choices = ['',]+G2gd.GetPatternTreeDataNames(G2frame,['IMG ',])
        backImage = wx.ComboBox(parent=G2frame.dataDisplay,value=data['background image'][0],choices=Choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        backImage.Bind(wx.EVT_COMBOBOX,OnBackImage)
        backSizer.Add(backImage)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' multiplier'),0,wx.ALIGN_CENTER_VERTICAL)
        backMult =  wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (data['background image'][1])),
            style=wx.TE_PROCESS_ENTER)
        backMult.Bind(wx.EVT_TEXT_ENTER,OnBackMult)
        backMult.Bind(wx.EVT_KILL_FOCUS,OnBackMult)
        backSizer.Add(backMult,0,wx.ALIGN_CENTER_VERTICAL)
        return backSizer
                        
    def CalibSizer():
                
        def OnNewCalibrant(event):
            data['calibrant'] = calSel.GetValue()
            data['calibskip'] = calFile.Calibrants[data['calibrant']][2]
            limits = calFile.Calibrants[data['calibrant']][3]
            data['calibdmin'],data['pixLimit'],data['cutoff'] = limits
            pixLimit.SetValue(str(limits[1]))
            cutOff.SetValue('%.1f'%(limits[2]))
            calibSkip.SetValue(str(data['calibskip']))
            calibDmin.SetValue('%.1f'%(limits[0]))
            
        def OnCalibSkip(event):
            data['calibskip'] = int(calibSkip.GetValue())
            
        def OnCalibDmin(event):
            try:
                dmin = float(calibDmin.GetValue())
                if dmin < 0.25:
                    raise ValueError
                data['calibdmin'] = dmin
            except ValueError:
                pass
            calibDmin.SetValue("%.2f"%(data['calibdmin']))          #reset in case of error  
                    
        def OnCutOff(event):
            try:
                cutoff = float(cutOff.GetValue())
                if cutoff < 0.1:
                    raise ValueError
                data['cutoff'] = cutoff
            except ValueError:
                pass
            cutOff.SetValue("%.1f"%(data['cutoff']))          #reset in case of error  
        
        def OnPixLimit(event):
            data['pixLimit'] = int(pixLimit.GetValue())
            
        def OnSetRings(event):
            if data['setRings']:
                data['setRings'] = False
            else:
                data['setRings'] = True
            G2plt.PlotExposedImage(G2frame,event=event)
    
        calibSizer = wx.FlexGridSizer(2,3,5,5)
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibrant '),0,
            wx.ALIGN_CENTER_VERTICAL)
        calSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['calibrant'],choices=calList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
        calSel.Bind(wx.EVT_COMBOBOX, OnNewCalibrant)
        comboSizer.Add(calSel,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calib lines to skip   '),0,
            wx.ALIGN_CENTER_VERTICAL)
        calibSkip  = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['calibskip']),choices=[str(i) for i in range(25)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        calibSkip.Bind(wx.EVT_COMBOBOX, OnCalibSkip)
        comboSizer.Add(calibSkip,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)        
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min calib d-spacing '),0,
            wx.ALIGN_CENTER_VERTICAL)
        calibDmin = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f" % (data['calibdmin'])),
            style=wx.TE_PROCESS_ENTER)
        calibDmin.Bind(wx.EVT_TEXT_ENTER,OnCalibDmin)
        calibDmin.Bind(wx.EVT_KILL_FOCUS,OnCalibDmin)
        comboSizer.Add(calibDmin,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min ring I/Ib '),0,
            wx.ALIGN_CENTER_VERTICAL)
        cutOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.1f" % (data['cutoff'])),
            style=wx.TE_PROCESS_ENTER)
        cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)
        cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
        comboSizer.Add(cutOff,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Pixel search range '),0,
            wx.ALIGN_CENTER_VERTICAL)
        pixLimit = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['pixLimit']),choices=['1','2','5','10','15','20'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        pixLimit.Bind(wx.EVT_COMBOBOX, OnPixLimit)
        comboSizer.Add(pixLimit,0,wx.ALIGN_CENTER_VERTICAL)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        setRings = wx.CheckBox(parent=G2frame.dataDisplay,label='Show ring picks?')
        comboSizer.Add(setRings,0)
        setRings.Bind(wx.EVT_CHECKBOX, OnSetRings)
        setRings.SetValue(data['setRings'])
        calibSizer.Add(comboSizer,0)
        return calibSizer
        
    def GonioSizer():
        
        ValObj = {}
        
        def OnGonioAngle(event):
            Obj = event.GetEventObject()
            item = ValObj[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = data['GonioAngles'][item]
            data['GonioAngles'][item] = value
            Obj.SetValue('%8.2f'%(value))
        
        gonioSizer = wx.BoxSizer(wx.HORIZONTAL)
        names = ['Omega','Chi','Phi']
        gonioSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'Sample goniometer angles: '),0,wx.ALIGN_CENTER_VERTICAL)
        for i,name in enumerate(names):
            gonioSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,name),0,wx.ALIGN_CENTER_VERTICAL)
            angle = wx.TextCtrl(G2frame.dataDisplay,-1,value='%8.2f'%(data['GonioAngles'][i]),
                style=wx.TE_PROCESS_ENTER)
            angle.Bind(wx.EVT_TEXT_ENTER,OnGonioAngle)
            angle.Bind(wx.EVT_KILL_FOCUS,OnGonioAngle)
            ValObj[angle.GetId()] = i
            gonioSizer.Add(angle,0,wx.ALIGN_CENTER_VERTICAL)
        return gonioSizer
        
# Image Controls main code             
                            
    #fix for old files:
    if 'azmthOff' not in data:
        data['azmthOff'] = 0.0
    if 'background image' not in data:
        data['background image'] = ['',1.0]
    if 'centerAzm' not in data:
        data['centerAzm'] = False
    if 'Oblique' not in data:
        data['Oblique'] = [0.5,False]
    if 'PolaVal' not in data:
        data['PolaVal'] = [0.99,False]
    #end fix
    
    colorList = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
    calList = [m for m in calFile.Calibrants.keys()]
    typeList = ['PWDR - powder diffraction data','SASD - small angle scattering data',
        'REFL - reflectometry data']
    if not data.get('type'):                        #patch for old project files
        data['type'] = 'PWDR'
    typeDict = {'PWDR':typeList[0],'SASD':typeList[1],'REFL':typeList[2]}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ImageMenu)
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCalibrate, id=G2gd.wxID_IMCALIBRATE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnRecalibrate, id=G2gd.wxID_IMRECALIBRATE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnClearCalib, id=G2gd.wxID_IMCLEARCALIB)
    if not data['rings']:
        G2frame.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMRECALIBRATE,enable=False)    
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnIntegrate, id=G2gd.wxID_IMINTEGRATE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnIntegrateAll, id=G2gd.wxID_INTEGRATEALL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyControls, id=G2gd.wxID_IMCOPYCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveControls, id=G2gd.wxID_IMSAVECONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadControls, id=G2gd.wxID_IMLOADCONTROLS)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)

    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)    
    mainSizer.Add(ComboSizer(),0,wx.ALIGN_LEFT)
    mainSizer.Add((5,5),0)            
    mainSizer.Add(MaxSizer(),0,wx.ALIGN_LEFT|wx.EXPAND)
    
    mainSizer.Add((5,5),0)
    DataSizer = wx.FlexGridSizer(1,2,5,5)
    DataSizer.Add(CalibCoeffSizer(),0)
    DataSizer.Add(IntegrateSizer(),0)        
    mainSizer.Add(DataSizer,0)
    mainSizer.Add((5,5),0)            
    mainSizer.Add(BackSizer(),0)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibration controls:'),0,
        wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(CalibSizer(),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(GonioSizer(),0,wx.ALIGN_CENTER_VERTICAL)   
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    fitSize = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(fitSize)
    G2frame.dataDisplay.SetSize(fitSize)
    
################################################################################
##### Masks
################################################################################
    
def UpdateMasks(G2frame,data):
    '''Shows and handles the controls on the "Masks"
    data tree entry
    '''
    
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
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def OnSpotDiameter(event):
        Obj = event.GetEventObject()
        try:
            diameter = min(100.,max(0.1,float(Obj.GetValue())))
        except ValueError:
            diameter = 1.0
        Obj.SetValue("%.2f"%(diameter))
        data['Points'][spotIds.index(Obj.GetId())][2] = diameter
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def OnDeleteSpot(event):
        Obj = event.GetEventObject()
        del(data['Points'][delSpotId.index(Obj)])
        UpdateMasks(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def OnRingThickness(event):
        Obj = event.GetEventObject()
        try:
            thick = min(1.0,max(0.001,float(Obj.GetValue())))
        except ValueError:
            thick = 0.1
        Obj.SetValue("%.3f"%(thick))
        data['Rings'][ringIds.index(Obj.GetId())][1] = thick
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def OnDeleteRing(event):
        Obj = event.GetEventObject()
        del(data['Rings'][delRingId.index(Obj)])
        UpdateMasks(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def OnArcThickness(event):
        Obj = event.GetEventObject()
        try:
            thick = min(20.0,max(0.001,float(Obj.GetValue())))
        except ValueError:
            thick = 0.1
        Obj.SetValue("%.3f"%(thick))
        data['Arcs'][arcIds.index(Obj.GetId())][2] = thick
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def OnDeleteArc(event):
        Obj = event.GetEventObject()
        del(data['Arcs'][delArcId.index(Obj)])
        UpdateMasks(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def OnDeletePoly(event):
        Obj = event.GetEventObject()
        del(data['Polygons'][delPolyId.index(Obj)])
        UpdateMasks(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def OnDeleteFrame(event):
        data['Frames'] = []
        UpdateMasks(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def OnCopyMask(event):
        import copy
        TextList = [[False,'All IMG',0]]
        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    if id == G2frame.Image:
                        Source = name
                        Mask = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Masks')))
                    else:
                        TextList.append([False,name,id])
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            if len(TextList) == 1:
                G2frame.ErrorDialog('Nothing to copy mask to','There must be more than one "IMG" pattern')
                return
            dlg = G2frame.CopyDialog(G2frame,'Copy mask information','Copy mask from '+Source+' to:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    if result[0][0]:
                        result = TextList[1:]
                        for item in result: item[0] = True
                    for i,item in enumerate(result):
                        ifcopy,name,id = item
                        if ifcopy:
                            mask = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Masks'))
                            Mask['Thresholds'][0] = mask['Thresholds'][0]
                            Mask['Thresholds'][1][1] = min(mask['Thresholds'][1][1],Mask['Thresholds'][1][1])
                            mask.update(Mask)                                
                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Masks'),copy.deepcopy(mask))
            finally:
                dlg.Destroy()
                
    def OnSaveMask(event):
        dlg = wx.FileDialog(G2frame, 'Choose image mask file', '.', '', 
            'image mask files (*.immask)|*.immask',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'w')
                save = {}
                keys = ['Points','Rings','Arcs','Polygons','Frames','Thresholds']
                for key in keys:
                    File.write(key+':'+str(data[key])+'\n')
                File.close()
        finally:
            dlg.Destroy()
        
    def OnLoadMask(event):
        dlg = wx.FileDialog(G2frame, 'Choose image mask file', '.', '', 
            'image mask files (*.immask)|*.immask',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                save = {}
                oldThreshold = data['Thresholds'][0]
                S = File.readline()
                while S:
                    if S[0] == '#':
                        S = File.readline()
                        continue
                    [key,val] = S[:-1].split(':')
                    if key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds']:
                        save[key] = eval(val)
                        if key == 'Thresholds':
                            save[key][0] = oldThreshold
                            save[key][1][1] = min(oldThreshold[1],save[key][1][1])
                    S = File.readline()
                data.update(save)
                UpdateMasks(G2frame,data)
                G2plt.PlotExposedImage(G2frame,event=event)
                
                File.close()
        finally:
            dlg.Destroy()
            
    def OnNewSpotMask(event):
        'Start a new spot mask'
        G2frame.MaskKey = 's'
        G2plt.OnStartMask(G2frame)
        
    def OnNewArcMask(event):
        'Start a new arc mask'
        G2frame.MaskKey = 'a'
        G2plt.OnStartMask(G2frame)
        
    def OnNewRingMask(event):
        'Start a new ring mask'
        G2frame.MaskKey = 'r'
        G2plt.OnStartMask(G2frame)
        
    def OnNewPolyMask(event):
        'Start a new polygon mask'
        G2frame.MaskKey = 'p'
        G2plt.OnStartMask(G2frame)
        
    def OnNewFrameMask(event):
        'Start a new Frame mask'
        G2frame.MaskKey = 'f'
        G2plt.OnStartMask(G2frame)
                
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.MaskMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyMask, id=G2gd.wxID_MASKCOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadMask, id=G2gd.wxID_MASKLOAD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveMask, id=G2gd.wxID_MASKSAVE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnNewSpotMask, id=G2gd.wxID_NEWMASKSPOT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnNewArcMask, id=G2gd.wxID_NEWMASKARC)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnNewRingMask, id=G2gd.wxID_NEWMASKRING)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnNewPolyMask, id=G2gd.wxID_NEWMASKPOLY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnNewFrameMask, id=G2gd.wxID_NEWMASKFRAME)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    if G2frame.MaskKey == 'f':
        G2frame.dataFrame.GetStatusBar().SetStatusText('Frame mask active - LB pick next point, RB close polygon')
    elif G2frame.MaskKey == 'p':
        G2frame.dataFrame.GetStatusBar().SetStatusText('Polygon mask active - LB pick next point, RB close polygon')
    elif G2frame.MaskKey == 's':
        G2frame.dataFrame.GetStatusBar().SetStatusText('Spot mask active - LB pick spot location')
    elif G2frame.MaskKey == 'a':
        G2frame.dataFrame.GetStatusBar().SetStatusText('Arc mask active - LB pick arc location')
    elif G2frame.MaskKey == 'r':
        G2frame.dataFrame.GetStatusBar().SetStatusText('Ring mask active - LB pick ring location')
    else:
        G2frame.dataFrame.GetStatusBar().SetStatusText("To add mask: On 2D Powder Image, key a:arc, r:ring, s:spot, p:polygon, f:frame")
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)

    thresh = data['Thresholds']         #min/max intensity range
    spots = data['Points']               #x,y,radius in mm
    rings = data['Rings']               #radius, thickness
    polygons = data['Polygons']         #3+ x,y pairs
    if 'Frames' not in data:
        data['Frames'] = []
    frame = data['Frames']             #3+ x,y pairs
    arcs = data['Arcs']                 #radius, start/end azimuth, thickness
    
    littleSizer = wx.FlexGridSizer(2,3,0,5)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Lower/Upper limits '),0,
        wx.ALIGN_CENTER_VERTICAL)
    Text = wx.TextCtrl(G2frame.dataDisplay,value=("%8d" % (thresh[0][0])),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,wx.ALIGN_CENTER_VERTICAL)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    Text = wx.TextCtrl(G2frame.dataDisplay,value=("%8d" % (thresh[0][1])),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,wx.ALIGN_CENTER_VERTICAL)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Lower/Upper thresholds '),
        0,wx.ALIGN_CENTER_VERTICAL)
    lowerThreshold = wx.TextCtrl(parent=G2frame.dataDisplay,
        value=("%8d" % (thresh[1][0])),style=wx.TE_PROCESS_ENTER)
    lowerThreshold.Bind(wx.EVT_TEXT_ENTER,OnThreshold)
    lowerThreshold.Bind(wx.EVT_KILL_FOCUS,OnThreshold)
    littleSizer.Add(lowerThreshold,0,wx.ALIGN_CENTER_VERTICAL)
    upperThreshold = wx.TextCtrl(parent=G2frame.dataDisplay,
        value=("%8d" % (thresh[1][1])),style=wx.TE_PROCESS_ENTER)
    upperThreshold.Bind(wx.EVT_TEXT_ENTER,OnThreshold)
    upperThreshold.Bind(wx.EVT_KILL_FOCUS,OnThreshold)
    littleSizer.Add(upperThreshold,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0,)
    spotIds = []
    delSpotId = []
    if spots:
        littleSizer = wx.FlexGridSizer(len(spots)+2,3,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Spot masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' position, mm'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' diameter, mm'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for spot in spots:
            if spot:
                x,y,d = spot
                spotText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f,%.2f" % (x,y)),
                    style=wx.TE_READONLY)
                spotText.SetBackgroundColour(VERY_LIGHT_GREY)
                littleSizer.Add(spotText,0,wx.ALIGN_CENTER_VERTICAL)
                spotText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                spotDiameter = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f" % (d)),
                    style=wx.TE_PROCESS_ENTER)
                littleSizer.Add(spotDiameter,0,wx.ALIGN_CENTER_VERTICAL)
                spotDiameter.Bind(wx.EVT_TEXT_ENTER,OnSpotDiameter)
                spotDiameter.Bind(wx.EVT_KILL_FOCUS,OnSpotDiameter)
                spotIds.append(spotDiameter.GetId())
                spotDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
                spotDelete.Bind(wx.EVT_CHECKBOX,OnDeleteSpot)
                delSpotId.append(spotDelete)
                littleSizer.Add(spotDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    ringIds = []
    delRingId = []
    if rings:
        littleSizer = wx.FlexGridSizer(len(rings)+2,3,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Ring masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' 2-theta,deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' thickness, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for ring in rings:
            if ring:
                tth,thick = ring
                ringText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (tth)),
                    style=wx.TE_READONLY)
                ringText.SetBackgroundColour(VERY_LIGHT_GREY)
                ringText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(ringText,0,wx.ALIGN_CENTER_VERTICAL)
                ringThick = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (thick)),
                    style=wx.TE_PROCESS_ENTER)
                littleSizer.Add(ringThick,0,wx.ALIGN_CENTER_VERTICAL)
                ringThick.Bind(wx.EVT_TEXT_ENTER,OnRingThickness)
                ringThick.Bind(wx.EVT_KILL_FOCUS,OnRingThickness)
                ringIds.append(ringThick.GetId())
                ringDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
                ringDelete.Bind(wx.EVT_CHECKBOX,OnDeleteRing)
                delRingId.append(ringDelete)
                littleSizer.Add(ringDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    arcIds = []
    delArcId = []
    if arcs:
        littleSizer = wx.FlexGridSizer(len(rings)+2,4,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Arc masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add((5,0),0)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' 2-theta,deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' azimuth, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' thickness, deg'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for arc in arcs:
            if arc:
                tth,azimuth,thick = arc
                arcText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (tth)),
                    style=wx.TE_READONLY)
                arcText.SetBackgroundColour(VERY_LIGHT_GREY)
                arcText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(arcText,0,wx.ALIGN_CENTER_VERTICAL)
                azmText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%d,%d" % (azimuth[0],azimuth[1])),
                    style=wx.TE_READONLY)
                azmText.SetBackgroundColour(VERY_LIGHT_GREY)
                azmText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(azmText,0,wx.ALIGN_CENTER_VERTICAL)
                arcThick = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (thick)),
                    style=wx.TE_PROCESS_ENTER)
                littleSizer.Add(arcThick,0,wx.ALIGN_CENTER_VERTICAL)
                arcThick.Bind(wx.EVT_TEXT_ENTER,OnArcThickness)
                arcThick.Bind(wx.EVT_KILL_FOCUS,OnArcThickness)
                arcIds.append(arcThick.GetId())
                arcDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
                arcDelete.Bind(wx.EVT_CHECKBOX,OnDeleteArc)
                delArcId.append(arcDelete)
                littleSizer.Add(arcDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    polyIds = []
    delPolyId = []
    delFrameId = []
    if polygons:
        littleSizer = wx.FlexGridSizer(len(polygons)+2,2,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Polygon masks:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        for polygon in polygons:
            if polygon:
                polyList = []
                for x,y in polygon:
                    polyList.append("%.2f, %.2f"%(x,y))
                polyText = wx.ComboBox(G2frame.dataDisplay,value=polyList[0],choices=polyList,style=wx.CB_READONLY)
                littleSizer.Add(polyText,0,wx.ALIGN_CENTER_VERTICAL)
                polyDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
                polyDelete.Bind(wx.EVT_CHECKBOX,OnDeletePoly)
                delPolyId.append(polyDelete)
                littleSizer.Add(polyDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    if frame:
        littleSizer = wx.FlexGridSizer(3,2,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Frame mask:'),0,
            wx.ALIGN_CENTER_VERTICAL)
        littleSizer.Add((5,0),0)
        frameList = []
        for x,y in frame:
            frameList.append("%.2f, %.2f"%(x,y))
        frameText = wx.ComboBox(G2frame.dataDisplay,value=frameList[0],choices=frameList,style=wx.CB_READONLY)
        littleSizer.Add(frameText,0,wx.ALIGN_CENTER_VERTICAL)
        frameDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
        frameDelete.Bind(wx.EVT_CHECKBOX,OnDeleteFrame)
        delFrameId.append(frameDelete)
        littleSizer.Add(frameDelete,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(littleSizer,0,)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] = 450
    G2frame.dataFrame.setSizePosLeft(Size)    

################################################################################
##### Stress/Strain
################################################################################

def UpdateStressStrain(G2frame,data):
    '''Shows and handles the controls on the "Stress/Strain"
    data tree entry
    '''
    
    def OnAppendDzero(event):
        data['d-zero'].append({'Dset':1.0,'Dcalc':0.0,'pixLimit':10,'cutoff':10.0,'ImxyObs':[[],[]],'ImxyCalc':[[],[]]})
        UpdateStressStrain(G2frame,data)
            
    def OnCopyStrSta(event):
        import copy
        TextList = [[False,'All IMG',0]]
        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
                Names.append(name)
                if 'IMG' in name:
                    if id == G2frame.Image:
                        Source = name
                        Data = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Stress/Strain')))
                    else:
                        TextList.append([False,name,id])
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            if len(TextList) == 1:
                G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "IMG" pattern')
                return
            dlg = G2frame.CopyDialog(G2frame,'Copy stress/strain controls','Copy controls from '+Source+' to:',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    if result[0][0]:
                        result = TextList[1:]
                        for item in result: item[0] = True
                    for i,item in enumerate(result):
                        ifcopy,name,id = item
                        if ifcopy:
                            oldData = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Stress/Strain')))
                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Stress/Strain'),copy.deepcopy(Data))
            finally:
                dlg.Destroy()

    def OnLoadStrSta(event):
        print 'Load stress/strain data - does nothing yet'
        event.Skip()

    def OnSaveStrSta(event):
        dlg = wx.FileDialog(G2frame, 'Choose stress/strain file', '.', '', 
            'image control files (*.strsta)|*.strsta',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'w')
                save = {}
                keys = ['Type','Sample phi','Sample z','strain']
                keys2 = ['Dset','Dcalc','pixLimit','cutoff']
                File.write('{\n\t')
                for key in keys:
                    if key in 'strain':
                        File.write("'"+key+"':["+str(data[key][0])+','+str(data[key][1])+','+str(data[key][2])+'],')
                    else:
                        File.write("'"+key+"':"+str(data[key])+',')
                File.write('\n\t'+"'d-zero':[\n")
                for data2 in data['d-zero']:
                    File.write('\t\t{')
                    for key in keys2:
                        File.write("'"+key+"':"+':'+str(data2[key])+',')
                    File.write("'ImxyObs':[[],[]],'Imxycalc':[[],[]]},\n")
                File.write('\t]\n}')
                File.close()
        finally:
            dlg.Destroy()

    def OnFitStrSta(event):
        Masks = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
        Controls = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
        G2img.FitStrSta(G2frame.ImageZ,data,Controls,Masks)
        G2plt.PlotExposedImage(G2frame,event=event)
        
    def SamSizer():
        
        def OnStrainType(event):
            data['Type'] = strType.GetValue()
        
        def OnSamPhi(event):
            try:
                value = float(samPhi.GetValue())
            except ValueError:
                value = data['Sample phi']
            data['Sample phi'] = value
            samPhi.SetValue("%.3f" % (data['Sample phi']))
                
        def OnSamZ(event):
            try:
                value = float(samZ.GetValue())
            except ValueError:
                value = data['Sample z']
            data['Sample z'] = value
            samZ.SetValue("%.3f" % (data['Sample z']))
                
        samSizer = wx.BoxSizer(wx.HORIZONTAL)
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Strain type: '),0,wx.ALIGN_CENTER_VERTICAL)
        strType = wx.ComboBox(G2frame.dataDisplay,value=data['Type'],choices=['True','Conventional'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        strType.SetValue(data['Type'])
        strType.Bind(wx.EVT_COMBOBOX, OnStrainType)
        samSizer.Add(strType,0,wx.ALIGN_CENTER_VERTICAL)
        
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Sample phi: '),0,wx.ALIGN_CENTER_VERTICAL)
        samPhi = wx.TextCtrl(G2frame.dataDisplay,-1,value=("%.3f" % (data['Sample phi'])),
            style=wx.TE_PROCESS_ENTER)
        samSizer.Add(samPhi,0,wx.ALIGN_CENTER_VERTICAL)
        samPhi.Bind(wx.EVT_TEXT_ENTER,OnSamPhi)
        samPhi.Bind(wx.EVT_KILL_FOCUS,OnSamPhi)
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Sample delta-z(mm): '),0,wx.ALIGN_CENTER_VERTICAL)
        samZ = wx.TextCtrl(G2frame.dataDisplay,-1,value=("%.3f" % (data['Sample z'])),
            style=wx.TE_PROCESS_ENTER)
        samSizer.Add(samZ,0,wx.ALIGN_CENTER_VERTICAL)
        samZ.Bind(wx.EVT_TEXT_ENTER,OnSamZ)
        samZ.Bind(wx.EVT_KILL_FOCUS,OnSamZ)
        return samSizer
        
    def DzeroSizer():
    
        def OnDzero(event):
            Obj = event.GetEventObject()
            try:
                value = min(10.0,max(1.0,float(Obj.GetValue())))
            except ValueError:
                value = 1.0
            Obj.SetValue("%.5f"%(value))
            data['d-zero'][Indx[Obj.GetId()]]['Dset'] = value
            
        def OnDeleteDzero(event):
            Obj = event.GetEventObject()
            del(data['d-zero'][delIndx.index(Obj)])
            UpdateStressStrain(G2frame,data)
        
        def OnCutOff(event):
            Obj = event.GetEventObject()
            try:
                value = min(10.0,max(0.5,float(Obj.GetValue())))
            except ValueError:
                value = 10.0
            Obj.SetValue("%.1f"%(value))
            data['d-zero'][Indx[Obj.GetId()]]['cutoff'] = value 
        
        def OnPixLimit(event):
            Obj = event.GetEventObject()
            data['d-zero'][Indx[Obj.GetId()]]['pixLimit'] = int(Obj.GetValue())
            
        Indx = {}
        delIndx = []    
        dzeroSizer = wx.FlexGridSizer(1,8,5,5)
        for id,dzero in enumerate(data['d-zero']):
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' d-zero #%d: '%(id))),
                0,wx.ALIGN_CENTER_VERTICAL)
            dZero = wx.TextCtrl(G2frame.dataDisplay,-1,value=('%.5f'%(dzero['Dset'])),
                style=wx.TE_PROCESS_ENTER)
            dzeroSizer.Add(dZero,0,wx.ALIGN_CENTER_VERTICAL)
            dZero.Bind(wx.EVT_TEXT_ENTER,OnDzero)
            dZero.Bind(wx.EVT_KILL_FOCUS,OnDzero)
            Indx[dZero.GetId()] = id
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' d-zero calc: %.5f'%(dzero['Dcalc']))),
                0,wx.ALIGN_CENTER_VERTICAL)
                
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min ring I/Ib '),0,
                wx.ALIGN_CENTER_VERTICAL)
            cutOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.1f" % (dzero['cutoff'])),
                style=wx.TE_PROCESS_ENTER)
            cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)
            cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
            Indx[cutOff.GetId()] = id
            dzeroSizer.Add(cutOff,0,wx.ALIGN_CENTER_VERTICAL)
        
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Pixel search range '),0,
                wx.ALIGN_CENTER_VERTICAL)
            pixLimit = wx.ComboBox(parent=G2frame.dataDisplay,value=str(dzero['pixLimit']),choices=['1','2','5','10','15','20'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            pixLimit.Bind(wx.EVT_COMBOBOX, OnPixLimit)
            Indx[pixLimit.GetId()] = id
            dzeroSizer.Add(pixLimit,0,wx.ALIGN_CENTER_VERTICAL)                
                
            dzeroDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
            dzeroDelete.Bind(wx.EVT_CHECKBOX,OnDeleteDzero)
            delIndx.append(dzeroDelete)
            dzeroSizer.Add(dzeroDelete,0,wx.ALIGN_CENTER_VERTICAL)
        return dzeroSizer
        
    def StrainSizer():
        
        strainSizer = wx.BoxSizer(wx.VERTICAL)
        strainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' Strain tensor:')),
            0,wx.ALIGN_CENTER_VERTICAL)
        tensorSizer = wx.FlexGridSizer(3,6,5,5)
        names = [[' e11','e12','e13'],[' e21','e22','e23'],[' e31','e32','e33']]
        for i in range(3):
            for j in range(3):
                tensorSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=names[i][j]),0,wx.ALIGN_CENTER_VERTICAL)
                tensorElem = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2f'%(data['strain'][i][j]),style=wx.TE_READONLY)
                tensorElem.SetBackgroundColour(VERY_LIGHT_GREY)
                tensorSizer.Add(tensorElem,0,wx.ALIGN_CENTER_VERTICAL)
        strainSizer.Add(tensorSizer)
        return strainSizer

    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.StrStaMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAppendDzero, id=G2gd.wxID_APPENDDZERO)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitStrSta, id=G2gd.wxID_STRSTAFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyStrSta, id=G2gd.wxID_STRSTACOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadStrSta, id=G2gd.wxID_STRSTALOAD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveStrSta, id=G2gd.wxID_STRSTASAVE)    
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText(" test  ")
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    mainSizer.Add(SamSizer())
    mainSizer.Add((5,10),0)
    mainSizer.Add(DzeroSizer())
    mainSizer.Add((5,10),0)
    mainSizer.Add(StrainSizer())
    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))    
