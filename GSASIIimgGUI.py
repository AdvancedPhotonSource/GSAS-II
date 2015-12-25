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
import os
import copy
import glob
import re
import bisect
import math
import time
import copy
import sys
import wx
import wx.lib.scrolledpanel as wxscroll
import wx.lib.mixins.listctrl  as  listmix
import matplotlib as mpl
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIimage as G2img
import GSASIImath as G2mth
import GSASIIplot as G2plt
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import GSASIIctrls as G2G
import GSASIIpy3 as G2py3

VERY_LIGHT_GREY = wx.Colour(235,235,235)
WACV = wx.ALIGN_CENTER_VERTICAL

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
    
################################################################################
##### Image Data
################################################################################

def UpdateImageData(G2frame,data):
    
    def OnPixVal(event):
        Obj = event.GetEventObject()
        id = Indx[Obj.GetId()]
        try:
            data['pixelSize'][id] = min(500,max(10,float(Obj.GetValue())))
        except ValueError:
            pass
        Obj.SetValue('%.3f'%(data['pixelSize'][id]))
        G2plt.PlotExposedImage(G2frame,newPlot=True,event=event)
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,
        label='Do not change anything here unless you are absolutely sure!'),0,WACV)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Image size: %d by %d'%(data['size'][0],data['size'][1])),0,WACV)
    pixSize = wx.FlexGridSizer(0,4,5,5)
    pixLabels = [u' Pixel X-dimension (\xb5m)',u' Pixel Y-dimension (\xb5m)']
    Indx = {}
    for i,[pixLabel,pix] in enumerate(zip(pixLabels,data['pixelSize'])):
        pixSize.Add(wx.StaticText(G2frame.dataDisplay,label=pixLabel),0,WACV)
        pixVal = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(pix),style=wx.TE_PROCESS_ENTER)
        Indx[pixVal.GetId()] = i
        pixVal.Bind(wx.EVT_TEXT_ENTER,OnPixVal)
        pixVal.Bind(wx.EVT_KILL_FOCUS,OnPixVal)
        pixSize.Add(pixVal,0,WACV)
    mainSizer.Add(pixSize,0)
    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    fitSize = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(fitSize)
    G2frame.dataDisplay.SetSize(fitSize)

################################################################################
##### Image Controls
################################################################################                    
def UpdateImageControls(G2frame,data,masks,IntegrateOnly=False):
    '''Shows and handles the controls on the "Image Controls"
    data tree entry
    '''
    import ImageCalibrants as calFile
#patch
    if 'Flat Bkg' not in data:
        data['Flat Bkg'] = 0.0
    if 'GonioAngles' not in data:
        data['GonioAngles'] = [0.,0.,0.]
    if 'DetDepth' not in data:
        data['DetDepth'] = 0.
    if 'SampleAbs' not in data:
        data['SampleShape'] = 'Cylinder'
        data['SampleAbs'] = [0.0,False]
    if 'binType' not in data:
        if 'PWDR' in data['type']:
            data['binType'] = '2-theta'
        elif 'SASD' in data['type']:
            data['binType'] = 'log(q)'
    if 'varyList' not in data:
        data['varyList'] = {'dist':True,'det-X':True,'det-Y':True,'tilt':True,'phi':True,'dep':False,'wave':False}
#end patch

# Menu items
            
    def OnCalibrate(event):
        if not data.get('calibrant'):
            G2G.G2MessageBox(G2frame,'A calibrant must first be selected')
        G2frame.dataFrame.GetStatusBar().SetStatusText('Select > 4 points on 1st used ring; LB to pick, RB on point to delete else RB to finish')
        G2frame.ifGetRing = True
                
    def OnRecalibrate(event):
        G2img.ImageRecalibrate(G2frame,data,masks)
        wx.CallLater(100,UpdateImageControls,G2frame,data,masks)
        
    def OnClearCalib(event):
        data['ring'] = []
        data['rings'] = []
        data['ellipses'] = []
        G2plt.PlotExposedImage(G2frame,event=event)
            
    def OnIntegrate(event):
        '''Integrate image in response to a menu event or from the AutoIntegrate
        dialog. In the latter case, event=None. 
        '''
        CleanupMasks(masks)
        blkSize = 128   #this seems to be optimal; will break in polymask if >1024
        Nx,Ny = data['size']
        nXBlks = (Nx-1)/blkSize+1
        nYBlks = (Ny-1)/blkSize+1
        Nup = nXBlks*nYBlks*3+3
        dlg = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        try:
            sumImg = G2frame.ImageZ
            darkImg,darkScale = data['dark image']
            if darkImg:
                Did = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, darkImg)
                Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(Did)
                darkImage = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)
                sumImg += darkImage*darkScale
            backImg,backScale = data['background image']            
            if backImg:     #ignores any transmission effect in the background image
                Bid = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, backImg)
                Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(Bid)
                backImage = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)
                Bdata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Bid,'Image Controls'))
                BdarkImg,BdarkScale = Bdata['dark image']
                if BdarkImg:
                    BDid = G2gd.GetPatternTreeItemId(G2frame, G2frame.root,BdarkImg)
                    Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(BDid)
                    BdarkImage = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)
                    backImage += BdarkImage*BdarkScale                
                sumImg += backImage*backScale
            G2frame.Integrate = G2img.ImageIntegrate(sumImg-data['Flat Bkg'],data,masks,blkSize,dlg)
#            G2plt.PlotIntegration(G2frame,newPlot=True)
            Id = G2IO.SaveIntegration(G2frame,G2frame.PickId,data,(event is None))
            G2frame.PatternId = Id
            G2frame.PatternTree.SelectItem(Id)
            G2frame.PatternTree.Expand(Id)
        finally:
            dlg.Destroy()
        for item in G2frame.MakePDF: item.Enable(True)
        
    def OnIntegrateAll(event):
        TextList = [[False,'All IMG',0]]
#        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
#                Names.append(name)
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
                    G2frame.EnablePlot = False
                    for item in result:
                        ifintegrate,name,id = item
                        if ifintegrate:
                            Id = G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls')
                            Data = G2frame.PatternTree.GetItemPyData(Id)
                            blkSize = 128   #this seems to be optimal; will break in polymask if >1024
                            Nx,Ny = Data['size']
                            nXBlks = (Nx-1)/blkSize+1
                            nYBlks = (Ny-1)/blkSize+1
                            Nup = nXBlks*nYBlks*3+3
                            dlgp = wx.ProgressDialog("Elapsed time","2D image integration",Nup,
                                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
                            try:
                                id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, name)
                                Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(id)
                                image = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)
                                backImage = []
                                if Data['background image'][0]:
                                    backImg = Data['background image'][0]
                                    backScale = Data['background image'][1]
                                    id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, backImg)
                                    Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(id)
                                    backImage = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)*backScale
                                FlatBkg = Data.get('Flat Bkg',0.0)
                                try:
                                    Masks = G2frame.PatternTree.GetItemPyData(
                                        G2gd.GetPatternTreeItemId(G2frame,id, 'Masks'))
                                except TypeError:       #missing Masks
                                    # I think the next line should be 'range'! (BHT)
                                    Imin,Imax = Data['Range']
                                    Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                                    G2frame.PatternTree.SetItemPyData(
                                        G2gd.GetPatternTreeItemId(G2frame,id, 'Masks'),Masks)
                                CleanupMasks(Masks)
                                if len(backImage):                                
                                    G2frame.Integrate = G2img.ImageIntegrate(image+backImage-FlatBkg,Data,Masks,blkSize,dlgp)
                                else:
                                    G2frame.Integrate = G2img.ImageIntegrate(image-FlatBkg,Data,Masks,blkSize,dlgp)
                                pId = G2IO.SaveIntegration(G2frame,Id,Data)
                            finally:
                                dlgp.Destroy()
                    else:
                        G2frame.EnablePlot = True
                        G2frame.PatternTree.SelectItem(pId)
                        G2frame.PatternTree.Expand(pId)
                        G2frame.PatternId = pId
                        
            finally:
                dlg.Destroy()
        
    def OnCopyControls(event):
        TextList = [[False,'All IMG',0]]
#        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
#                Names.append(name)
                if 'IMG' in name:
                    if id == G2frame.Image:
                        Source = name
                        Data = copy.deepcopy(data)
#                        Data = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls')))
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
                G2frame.PatternTree.SelectItem(G2frame.PickId)
                
    def OnSaveControls(event):
        dlg = wx.FileDialog(G2frame, 'Choose image controls file', '.', '', 
            'image control files (*.imctrl)|*.imctrl',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .imctrl
                filename = os.path.splitext(filename)[0]+'.imctrl'
                File = open(filename,'w')
                keys = ['type','wavelength','calibrant','distance','center',
                    'tilt','rotation','azmthOff','fullIntegrate','LRazimuth',
                    'IOtth','outChannels','outAzimuths','invert_x','invert_y','DetDepth',
                    'calibskip','pixLimit','cutoff','calibdmin','chisq','Flat Bkg',
                    'binType','SampleShape','PolaVal','SampleAbs','dark image','background image']
                for key in keys:
                    if key not in data:     #uncalibrated!
                        continue
                    File.write(key+':'+str(data[key])+'\n')
                File.close()
        finally:
            dlg.Destroy()
        
    def OnLoadControls(event):
        cntlList = ['wavelength','distance','tilt','invert_x','invert_y','type',
            'fullIntegrate','outChannels','outAzimuths','LRazimuth','IOtth','azmthOff','DetDepth',
            'calibskip','pixLimit','cutoff','calibdmin','chisq','Flat Bkg',
            'PolaVal','SampleAbs','dark image','background image']
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
                data.update(save)
                # next line removed. Previous updates tree contents. The next
                # makes a copy of data, puts it into tree and "disconnects" data
                # from tree contents (later changes to data are lost!)
                #G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'),copy.deepcopy(data))
                File.close()
        finally:
            dlg.Destroy()
        G2plt.PlotExposedImage(G2frame,event=event)
        wx.CallLater(100,UpdateImageControls,G2frame,data,masks)
            
# Sizers
                                        
    def ComboSizer():

        def OnDataType(event):
            data['type'] = typeSel.GetValue()[:4]
            if 'SASD' in data['type']:
                data['SampleAbs'][0] = np.exp(-data['SampleAbs'][0]) #switch from muT to trans!
                if data['binType'] == '2-theta': data['binType'] = 'log(q)'  #switch default bin type
            elif 'PWDR' in data['type']:
                data['SampleAbs'][0] = -np.log(data['SampleAbs'][0])  #switch from trans to muT!
                if data['binType'] == 'log(q)': data['binType'] = '2-theta'  #switch default bin type                 
            wx.CallLater(100,UpdateImageControls,G2frame,data,masks)
    
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
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Type of image data: '),0,WACV)
        typeSel = wx.ComboBox(parent=G2frame.dataDisplay,value=typeDict[data['type']],choices=typeList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        typeSel.SetValue(data['type'])
        typeSel.Bind(wx.EVT_COMBOBOX, OnDataType)
        comboSizer.Add(typeSel,0,WACV)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Color bar '),0,WACV)
        colSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['color'],choices=colorList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
        comboSizer.Add(colSel,0,WACV)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Azimuth offset '),0,WACV)
        azmthOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f" % (data['azmthOff'])),
            style=wx.TE_PROCESS_ENTER)
        azmthOff.Bind(wx.EVT_TEXT_ENTER,OnAzmthOff)
        azmthOff.Bind(wx.EVT_KILL_FOCUS,OnAzmthOff)
        comboSizer.Add(azmthOff,0,WACV)
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
            
        maxSizer = wx.FlexGridSizer(0,3,0,5)
        maxSizer.AddGrowableCol(1,1)
        maxSizer.SetFlexibleDirection(wx.HORIZONTAL)
        sqrtDeltZero = math.sqrt(data['range'][0][1]-max(0.0,data['range'][0][0]))
        DeltOne = data['range'][1][1]-max(0.0,data['range'][0][0])
        sqrtDeltOne = math.sqrt(DeltOne)
        maxSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Max intensity'),0,WACV)
        maxSel = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(100*sqrtDeltOne/sqrtDeltZero))
        maxSizer.Add(maxSel,1,wx.EXPAND)
        maxSel.Bind(wx.EVT_SLIDER, OnMaxSlider)
        maxVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.0f'%(data['range'][1][1]))
        maxVal.Bind(wx.EVT_TEXT_ENTER,OnMaxVal)    
        maxVal.Bind(wx.EVT_KILL_FOCUS,OnMaxVal)
        maxSizer.Add(maxVal,0,WACV)    
        maxSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min intensity'),0,WACV)
        minSel = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(100*(data['range'][1][0]-max(0.0,data['range'][0][0]))/DeltOne))
        maxSizer.Add(minSel,1,wx.EXPAND)
        minSel.Bind(wx.EVT_SLIDER, OnMinSlider)
        minVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.0f'%(data['range'][1][0]))
        minVal.Bind(wx.EVT_TEXT_ENTER,OnMinVal)    
        minVal.Bind(wx.EVT_KILL_FOCUS,OnMinVal)
        maxSizer.Add(minVal,0,WACV)
        return maxSizer
        
    def CalibCoeffSizer():
        
        def OnCalRef(event):
            Obj = event.GetEventObject()
            name = Indx[Obj]
            data['varyList'][name] = Obj.GetValue()
            
        def OnCalVal(event):
            Obj = event.GetEventObject()
            name = Indx[Obj]
            try:
                value = float(Obj.GetValue())
                if name == 'wave' and value < 0.01:
                    raise ValueError
            except ValueError:
                value = Parms[name][2]
            if name == 'dist':
                data['distance'] = value
            elif name == 'det-X':
                data['center'][0] = value
            elif name == 'det-Y':
                data['center'][1] = value
            elif name == 'tilt':
                data['tilt'] = value
            elif name == 'phi':
                data['rotation'] = value
            elif name == 'wave':
                data['wavelength'] = value
            elif name == 'dep':
                data['DetDepth'] = value                                
            Parms[name][2] = value
            Obj.SetValue(Parms[name][1]%(value))
            
        calibSizer = wx.FlexGridSizer(0,2,5,5)
        calibSizer.SetFlexibleDirection(wx.HORIZONTAL)
        calibSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibration coefficients'),0,WACV)    
        calibSizer.Add((5,0),0)
        cent = data['center']
        Names = ['det-X','det-Y','wave','dist','tilt','phi']
        if 'PWDR' in data['type']:
            Names.append('dep') 
        Parms = {'dist':['Distance','%.3f',data['distance']],'det-X':['Beam center X','%.3f',data['center'][0]],
            'det-Y':['Beam center Y','%.3f',data['center'][1]],'tilt':['Tilt angle','%.3f',data['tilt']],
            'phi':['Tilt rotation','%.2f',data['rotation']],'dep':['Penetration','%.2f',data['DetDepth']],
            'wave':['Wavelength','%.6f',data['wavelength']]}
        Indx = {}
        for name in Names:
            calSel = wx.CheckBox(parent=G2frame.dataDisplay,label=Parms[name][0])
            calibSizer.Add(calSel,0,WACV)
            calSel.Bind(wx.EVT_CHECKBOX, OnCalRef)
            calSel.SetValue(data['varyList'][name])
            Indx[calSel] = name
            calVal = wx.TextCtrl(G2frame.dataDisplay,value=(Parms[name][1]%(Parms[name][2])),style=wx.TE_PROCESS_ENTER)
            calVal.Bind(wx.EVT_TEXT_ENTER,OnCalVal)
            calVal.Bind(wx.EVT_KILL_FOCUS,OnCalVal)
            Indx[calVal] = name
            calibSizer.Add(calVal,0,WACV)
        return calibSizer
    
    def IntegrateSizer():
        
        def OnNewBinType(event):
            data['binType'] = binSel.GetValue()
            wx.CallLater(100,UpdateImageControls,G2frame,data,masks)
        
        def OnIOtth(event):
            Ltth = max(float(G2frame.InnerTth.GetValue()),0.001)
            Utth = float(G2frame.OuterTth.GetValue())
            if Ltth > Utth:
                Ltth,Utth = Utth,Ltth
            if 'Q' in data['binType']:
                data['IOtth'] = [2.*asind(Ltth*wave/(4.*math.pi)),2.*asind(Utth*wave/(4.*math.pi))]
            else:
                data['IOtth'] = [Ltth,Utth]
            G2frame.InnerTth.SetValue("%8.3f" % (Ltth))
            G2frame.OuterTth.SetValue("%8.3f" % (Utth))
            G2plt.PlotExposedImage(G2frame,event=event)
        
        def OnLRazim(event):
            Lazm = int(G2frame.Lazim.GetValue())%360
            Razm = int(G2frame.Razim.GetValue())%360
            if Lazm > Razm:
                Razm += 360
            if data['fullIntegrate']:
                Razm = Lazm+360
            G2frame.Lazim.SetValue("%6d" % (Lazm))
            G2frame.Razim.SetValue("%6d" % (Razm))
            data['LRazimuth'] = [Lazm,Razm]
            G2plt.PlotExposedImage(G2frame,event=event)
        
        def OnNumOutChans(event):
            try:
                numChans = int(outChan.GetValue())
                if numChans < 10:
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
                           
        def OnSamAbs(event):
            if data['SampleAbs'][1]:
                data['SampleAbs'][1] = False
            else:
                data['SampleAbs'][1] = True
                
        def OnSamAbsVal(event):
            try:
                value = float(samabsVal.GetValue())
                minmax = [0.,2.]
                if 'SASD' in data['type']:
                    minmax = [.05,1.0]
                if minmax[0] <= value <= minmax[1]:
                    data['SampleAbs'][0] = value
                else:
                    raise ValueError
            except ValueError:
                pass
            samabsVal.SetValue('%.3f'%(data['SampleAbs'][0]))
                           
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
            wx.CallLater(100,UpdateImageControls,G2frame,data,masks)
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnSetDefault(event):
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
                           
        dataSizer = wx.FlexGridSizer(0,2,5,3)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Integration coefficients'),0,WACV)    
        dataSizer.Add((5,0),0)
        if 'PWDR' in data['type']:
            binChoice = ['2-theta','Q']
        elif 'SASD' in data['type']:
            binChoice = ['2-theta','Q','log(q)']
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Bin style: Constant step bins in'),0,WACV)            
        binSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['binType'],choices=binChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        binSel.Bind(wx.EVT_COMBOBOX, OnNewBinType)
        dataSizer.Add(binSel,0,WACV)
        binType = '2-theta'
        if 'q' in data['binType'].lower():
            binType = 'Q'
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Inner/Outer '+binType),0,WACV)            
        IOtth = data['IOtth'][:]
        if 'Q' in data['binType']:
            wave = data['wavelength']
            IOtth = [4.*math.pi*sind(IOtth[0]/2.)/wave,4.*math.pi*sind(IOtth[1]/2.)/wave]
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        G2frame.InnerTth = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%8.3f" % (IOtth[0])),style=wx.TE_PROCESS_ENTER)
        G2frame.InnerTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
        G2frame.InnerTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
        littleSizer.Add(G2frame.InnerTth,0,WACV)
        G2frame.OuterTth = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%8.2f" % (IOtth[1])),style=wx.TE_PROCESS_ENTER)
        G2frame.OuterTth.Bind(wx.EVT_TEXT_ENTER,OnIOtth)
        G2frame.OuterTth.Bind(wx.EVT_KILL_FOCUS,OnIOtth)
        littleSizer.Add(G2frame.OuterTth,0,WACV)
        dataSizer.Add(littleSizer,0,)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Start/End azimuth'),0,WACV)
        LRazim = data['LRazimuth']
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        G2frame.Lazim = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%6d" % (LRazim[0])),style=wx.TE_PROCESS_ENTER)
        G2frame.Lazim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
        G2frame.Lazim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
        littleSizer.Add(G2frame.Lazim,0,WACV)
        G2frame.Razim = wx.TextCtrl(parent=G2frame.dataDisplay,
            value=("%6d" % (LRazim[1])),style=wx.TE_PROCESS_ENTER)
        G2frame.Razim.Bind(wx.EVT_TEXT_ENTER,OnLRazim)
        G2frame.Razim.Bind(wx.EVT_KILL_FOCUS,OnLRazim)
        if data['fullIntegrate']:
            G2frame.Razim.Enable(False)
            G2frame.Razim.SetBackgroundColour(VERY_LIGHT_GREY)
            G2frame.Razim.SetValue("%6d" % (LRazim[0]+360))
        littleSizer.Add(G2frame.Razim,0,WACV)
        dataSizer.Add(littleSizer,0,)
        dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' No. 2-theta/azimuth bins'),0,WACV)
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        outChan = wx.TextCtrl(parent=G2frame.dataDisplay,value=str(data['outChannels']),style=wx.TE_PROCESS_ENTER)
        outChan.Bind(wx.EVT_TEXT_ENTER,OnNumOutChans)
        outChan.Bind(wx.EVT_KILL_FOCUS,OnNumOutChans)
        littleSizer.Add(outChan,0,WACV)
        outAzim = wx.TextCtrl(parent=G2frame.dataDisplay,value=str(data['outAzimuths']),style=wx.TE_PROCESS_ENTER)
        outAzim.Bind(wx.EVT_TEXT_ENTER,OnNumOutAzms)
        outAzim.Bind(wx.EVT_KILL_FOCUS,OnNumOutAzms)
        littleSizer.Add(outAzim,0,WACV)
        dataSizer.Add(littleSizer,0,)
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        samabs = wx.CheckBox(parent=G2frame.dataDisplay,label='Apply sample absorption?')
        dataSizer.Add(samabs,0,WACV)
        samabs.Bind(wx.EVT_CHECKBOX, OnSamAbs)
        samabs.SetValue(data['SampleAbs'][1])
        if 'PWDR' in data['type']:
            littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='mu/R (0.00-2.0) '),0,WACV)
        elif 'SASD' in data['type']:
            littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='transmission '),0,WACV)
        samabsVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.3f'%(data['SampleAbs'][0]),style=wx.TE_PROCESS_ENTER)            
        samabsVal.Bind(wx.EVT_TEXT_ENTER,OnSamAbsVal)
        samabsVal.Bind(wx.EVT_KILL_FOCUS,OnSamAbsVal)
        littleSizer.Add(samabsVal,0,WACV)
        dataSizer.Add(littleSizer,0,)
        if 'PWDR' in data['type']:
            littleSizer = wx.BoxSizer(wx.HORIZONTAL)
            oblique = wx.CheckBox(parent=G2frame.dataDisplay,label='Apply detector absorption?')
            dataSizer.Add(oblique,0,WACV)
            oblique.Bind(wx.EVT_CHECKBOX, OnOblique)
            oblique.SetValue(data['Oblique'][1])
            littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Value (0.01-0.99)  '),0,WACV)
            obliqVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.3f'%(data['Oblique'][0]),style=wx.TE_PROCESS_ENTER)
            obliqVal.Bind(wx.EVT_TEXT_ENTER,OnObliqVal)
            obliqVal.Bind(wx.EVT_KILL_FOCUS,OnObliqVal)
            littleSizer.Add(obliqVal,0,WACV)
            dataSizer.Add(littleSizer,0,)
        if 'SASD' in data['type']:
            littleSizer = wx.BoxSizer(wx.HORIZONTAL)
            setPolariz = wx.CheckBox(parent=G2frame.dataDisplay,label='Apply polarization?')
            dataSizer.Add(setPolariz,0,WACV)
            setPolariz.Bind(wx.EVT_CHECKBOX, OnApplyPola)
            setPolariz.SetValue(data['PolaVal'][1])
            littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Value (0.001-0.999)  '),0,WACV)
            polaVal = wx.TextCtrl(parent=G2frame.dataDisplay,value='%.3f'%(data['PolaVal'][0]),
                style=wx.TE_PROCESS_ENTER)
            polaVal.Bind(wx.EVT_TEXT_ENTER,OnPolaVal)
            polaVal.Bind(wx.EVT_KILL_FOCUS,OnPolaVal)
            littleSizer.Add(polaVal,0,WACV)
            dataSizer.Add(littleSizer,0,)
        
        showLines = wx.CheckBox(parent=G2frame.dataDisplay,label='Show integration limits?')
        dataSizer.Add(showLines,0,WACV)
        showLines.Bind(wx.EVT_CHECKBOX, OnShowLines)
        showLines.SetValue(data['showLines'])
        fullIntegrate = wx.CheckBox(parent=G2frame.dataDisplay,label='Do full integration?')
        dataSizer.Add(fullIntegrate,0,WACV)
        fullIntegrate.Bind(wx.EVT_CHECKBOX, OnFullIntegrate)
        fullIntegrate.SetValue(data['fullIntegrate'])
        setDefault = wx.CheckBox(parent=G2frame.dataDisplay,label='Use as default for all images?')
        dataSizer.Add(setDefault,0,WACV)
        setDefault.Bind(wx.EVT_CHECKBOX, OnSetDefault)
        setDefault.SetValue(data['setDefault'])
        centerAzm = wx.CheckBox(parent=G2frame.dataDisplay,label='Azimuth at bin center?')
        dataSizer.Add(centerAzm,0,WACV)
        centerAzm.Bind(wx.EVT_CHECKBOX, OnCenterAzm)
        centerAzm.SetValue(data['centerAzm'])
        return dataSizer
        
    def BackSizer():
        
        def OnBackImage(event):
            data['background image'][0] = backImage.GetValue()
            
        def OnDarkImage(event):
            data['dark image'][0] = darkImage.GetValue()
            G2plt.PlotExposedImage(G2frame,event=event)
            
        def OnFlatBkg(event):
            try:
                value = float(flatbkg.GetValue())
                data['Flat Bkg'] = value
            except ValueError:
                pass
            flatbkg.SetValue("%.0f"%(data['Flat Bkg']))    
            G2plt.PlotExposedImage(G2frame,event=event)

        def OnBackMult(event):
            try:
                mult = float(backMult.GetValue())
                data['background image'][1] = mult
            except ValueError:
                pass
            backMult.SetValue("%.3f" % (data['background image'][1]))          #reset in case of error 
        
        def OnDarkMult(event):
            try:
                mult = float(darkMult.GetValue())
                data['dark image'][1] = mult
            except ValueError:
                pass
            darkMult.SetValue("%.3f" % (data['dark image'][1]))          #reset in case of error 
            G2plt.PlotExposedImage(G2frame,event=event)
        
        backSizer = wx.FlexGridSizer(0,6,5,5)

        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Dark image'),0,WACV)
        Choices = ['',]+G2gd.GetPatternTreeDataNames(G2frame,['IMG ',])
        darkImage = wx.ComboBox(parent=G2frame.dataDisplay,value=data['dark image'][0],choices=Choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        darkImage.Bind(wx.EVT_COMBOBOX,OnDarkImage)
        backSizer.Add(darkImage)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' multiplier'),0,WACV)
        darkMult =  wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (data['dark image'][1])),
            style=wx.TE_PROCESS_ENTER)
        darkMult.Bind(wx.EVT_TEXT_ENTER,OnDarkMult)
        darkMult.Bind(wx.EVT_KILL_FOCUS,OnDarkMult)
        backSizer.Add(darkMult,0,WACV)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Flat Bkg: '),0,WACV)
        flatbkg = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.0f" % (data['Flat Bkg'])),
            style=wx.TE_PROCESS_ENTER)
        flatbkg.Bind(wx.EVT_TEXT_ENTER,OnFlatBkg)
        flatbkg.Bind(wx.EVT_KILL_FOCUS,OnFlatBkg)
        backSizer.Add(flatbkg,0,WACV)

        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background image'),0,WACV)
        Choices = ['',]+G2gd.GetPatternTreeDataNames(G2frame,['IMG ',])
        backImage = wx.ComboBox(parent=G2frame.dataDisplay,value=data['background image'][0],choices=Choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        backImage.Bind(wx.EVT_COMBOBOX,OnBackImage)
        backSizer.Add(backImage)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' multiplier'),0,WACV)
        backMult =  wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (data['background image'][1])),
            style=wx.TE_PROCESS_ENTER)
        backMult.Bind(wx.EVT_TEXT_ENTER,OnBackMult)
        backMult.Bind(wx.EVT_KILL_FOCUS,OnBackMult)
        backSizer.Add(backMult,0,WACV)
        return backSizer
                        
    def CalibSizer():
                
        def OnNewCalibrant(event):
            data['calibrant'] = calSel.GetValue()
            data['calibskip'] = calFile.Calibrants[data['calibrant']][3]
            limits = calFile.Calibrants[data['calibrant']][4]
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
    
        calibSizer = wx.FlexGridSizer(0,3,5,5)
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibrant '),0,WACV)
        calSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['calibrant'],choices=calList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        calSel.Bind(wx.EVT_COMBOBOX, OnNewCalibrant)
        comboSizer.Add(calSel,0,WACV)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)    
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calib lines to skip   '),0,WACV)
        calibSkip  = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['calibskip']),choices=[str(i) for i in range(25)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        calibSkip.Bind(wx.EVT_COMBOBOX, OnCalibSkip)
        comboSizer.Add(calibSkip,0,WACV)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)        
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min calib d-spacing '),0,WACV)
        calibDmin = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f" % (data['calibdmin'])),
            style=wx.TE_PROCESS_ENTER)
        calibDmin.Bind(wx.EVT_TEXT_ENTER,OnCalibDmin)
        calibDmin.Bind(wx.EVT_KILL_FOCUS,OnCalibDmin)
        comboSizer.Add(calibDmin,0,WACV)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Min ring I/Ib '),0,WACV)
        cutOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.1f" % (data['cutoff'])),
            style=wx.TE_PROCESS_ENTER)
        cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)
        cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
        comboSizer.Add(cutOff,0,WACV)
        calibSizer.Add(comboSizer,0)
        
        comboSizer = wx.BoxSizer(wx.HORIZONTAL)
        comboSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Pixel search range '),0,WACV)
        pixLimit = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['pixLimit']),choices=['1','2','5','10','15','20'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        pixLimit.Bind(wx.EVT_COMBOBOX, OnPixLimit)
        comboSizer.Add(pixLimit,0,WACV)
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
            
        def OnGlobalEdit(event):
            Names = []
            Items = []
            if G2frame.PatternTree.GetCount():
                id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
                while id:
                    name = G2frame.PatternTree.GetItemText(id)
                    if 'IMG' in name:
                        ctrls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id,'Image Controls'))
                        Names.append(name)
                        Items.append(ctrls['GonioAngles'])
                    id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
                if len(Names) == 1:
                    G2frame.ErrorDialog('Nothing for global editing','There must be more than one "IMG" pattern')
                    return
                dlg = G2G.G2HistoDataDialog(G2frame,' Edit sample goniometer data:',
                    'Edit data',['Omega','Chi','Phi'],['%.2f','%.2f','%.2f'],Names,Items)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetData()
                    id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
                    while id:
                        name = G2frame.PatternTree.GetItemText(id)
                        if 'IMG' in name:
                            ctrls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id,'Image Controls'))
                            vals = Items[Names.index(name)]
                            ctrls['GonioAngles'] = vals
#                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Image Controls'),ctrls)
                        id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            finally:
                dlg.Destroy()
                G2frame.PatternTree.SelectItem(G2frame.PickId)
        
        gonioSizer = wx.BoxSizer(wx.HORIZONTAL)
        names = ['Omega','Chi','Phi']
        gonioSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'Sample goniometer angles: '),0,WACV)
        for i,name in enumerate(names):
            gonioSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,name),0,WACV)
            angle = wx.TextCtrl(G2frame.dataDisplay,-1,value='%8.2f'%(data['GonioAngles'][i]),
                style=wx.TE_PROCESS_ENTER)
            angle.Bind(wx.EVT_TEXT_ENTER,OnGonioAngle)
            angle.Bind(wx.EVT_KILL_FOCUS,OnGonioAngle)
            ValObj[angle.GetId()] = i
            gonioSizer.Add(angle,0,WACV)
        globEdit = wx.Button(G2frame.dataDisplay,-1,'Global edit')
        globEdit.Bind(wx.EVT_BUTTON,OnGlobalEdit)
        gonioSizer.Add(globEdit,0,WACV)
        return gonioSizer
        
# Image Controls main code             
                            
    #fix for old files:
    if 'azmthOff' not in data:
        data['azmthOff'] = 0.0
    if 'background image' not in data:
        data['background image'] = ['',-1.0]
    if 'dark image' not in data:
        data['dark image'] = ['',-1.0]
    if 'centerAzm' not in data:
        data['centerAzm'] = False
    if 'Oblique' not in data:
        data['Oblique'] = [0.5,False]
    if 'PolaVal' not in data:
        data['PolaVal'] = [0.99,False]
    #end fix
    
    if IntegrateOnly:
        OnIntegrate(None)
        return
    
    colorList = sorted([m for m in mpl.cm.datad.keys() if not m.endswith("_r")],key=lambda s: s.lower())
    calList = sorted([m for m in calFile.Calibrants.keys()],key=lambda s: s.lower())
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
    if data.get('calibrant'):
        G2frame.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMRECALIBRATE,enable=True)
    else:
        G2frame.dataFrame.ImageEdit.Enable(id=G2gd.wxID_IMRECALIBRATE,enable=False)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnIntegrate, id=G2gd.wxID_IMINTEGRATE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnIntegrateAll, id=G2gd.wxID_INTEGRATEALL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyControls, id=G2gd.wxID_IMCOPYCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveControls, id=G2gd.wxID_IMSAVECONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadControls, id=G2gd.wxID_IMLOADCONTROLS)
    def OnDestroy(event):
        G2frame.autoIntFrame = None
    def OnAutoInt(event):
        if G2frame.autoIntFrame: # ensure only one open at a time
            G2frame.autoIntFrame.Raise()
            return
        G2frame.autoIntFrame = AutoIntFrame(G2frame,PollTime=10.0)
        G2frame.autoIntFrame.Bind(wx.EVT_WINDOW_DESTROY,OnDestroy) # clean up name on window close 
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAutoInt, id=G2gd.wxID_IMAUTOINTEG)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)

    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)    
    mainSizer.Add(ComboSizer(),0,wx.ALIGN_LEFT)
    mainSizer.Add((5,5),0)            
    mainSizer.Add(MaxSizer(),0,wx.ALIGN_LEFT|wx.EXPAND)
    
    mainSizer.Add((5,5),0)
    DataSizer = wx.FlexGridSizer(0,2,5,5)
    DataSizer.Add(CalibCoeffSizer(),0)
    DataSizer.Add(IntegrateSizer(),0)        
    mainSizer.Add(DataSizer,0)
    mainSizer.Add((5,5),0)            
    mainSizer.Add(BackSizer(),0)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Calibration controls:'),0,WACV)
    mainSizer.Add((5,5),0)
    mainSizer.Add(CalibSizer(),0,WACV)
    mainSizer.Add((5,5),0)
    mainSizer.Add(GonioSizer(),0,WACV)   
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    fitSize = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(fitSize)
    G2frame.dataDisplay.SetSize(fitSize)
    
################################################################################
##### Masks
################################################################################
def CleanupMasks(data):
    '''If a mask creation is not completed, an empty mask entry is created in the
    masks array. This cleans them out. It is called when the masks page is first loaded
    and before saving them or after reading them in. This should also probably be done
    before they are used for integration.
    '''
    for key in ['Points','Rings','Arcs','Polygons']:
        data[key] = data.get(key,[])
        l1 = len(data[key])
        data[key] = [i for i in data[key] if i]
        l2 = len(data[key])
        if GSASIIpath.GetConfigValue('debug') and l1 != l2:
            print 'Mask Cleanup:',key,'was',l1,'entries','now',l2
    
def UpdateMasks(G2frame,data):
    '''Shows and handles the controls on the "Masks" data tree entry
    '''
    
    def OnTextMsg(event):
        Obj = event.GetEventObject()
        Obj.SetToolTipString('Drag this mask on 2D Powder Image with mouse to change ')

    def Replot(*args,**kwargs):
        #G2plt.PlotExposedImage(G2frame,newPlot=True)        
        G2plt.PlotExposedImage(G2frame) # remove newPlot=True to prevent zoom-out

    def onDeleteMask(event):
        Obj = event.GetEventObject()
        typ = Obj.locationcode.split('+')[1]
        num = int(Obj.locationcode.split('+')[2])
        del(data[typ][num])
        wx.CallAfter(UpdateMasks,G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def onDeleteFrame(event):
        data['Frames'] = []
        wx.CallAfter(UpdateMasks,G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)

    def OnCopyMask(event):
        TextList = [[False,'All IMG',0]]
#        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
#                Names.append(name)
                if 'IMG' in name:
                    if id == G2frame.Image:
                        Source = name
                        Mask = copy.deepcopy(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Masks')))
                        Thresh = Mask.pop('Thresholds')  #remove Thresholds from source mask & save it for later
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
#                            Mask['Thresholds'][0] = mask['Thresholds'][0]
#                            Mask['Thresholds'][1][1] = min(mask['Thresholds'][1][1],Mask['Thresholds'][1][1])
                            mask.update(Mask)
                            mask['Thresholds'][1][0] = Thresh[1][0]  #copy only lower threshold                             
                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Masks'),copy.deepcopy(mask))
            finally:
                dlg.Destroy()
                
    def OnSaveMask(event):
        CleanupMasks(data)
        dlg = wx.FileDialog(G2frame, 'Choose image mask file', '.', '', 
            'image mask files (*.immask)|*.immask',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.immask'
                File = open(filename,'w')
                keys = ['Points','Rings','Arcs','Polygons','Frames','Thresholds']
                for key in keys:
                    File.write(key+':'+str(data[key])+'\n')
                File.close()
        finally:
            dlg.Destroy()
        
    def OnLoadMask(event):
        if event.Id == G2gd.wxID_MASKLOADNOT:
            ignoreThreshold = True
        else:
            ignoreThreshold = False
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
                        if ignoreThreshold and key == 'Thresholds':
                            S = File.readline() 
                            continue
                        save[key] = eval(val)
                        if key == 'Thresholds':
                            save[key][0] = oldThreshold
                            save[key][1][1] = min(oldThreshold[1],save[key][1][1])
                    S = File.readline()
                File.close()
                data.update(save)
                CleanupMasks(data)
                wx.CallAfter(UpdateMasks,G2frame,data)
                G2plt.PlotExposedImage(G2frame,event=event)                
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

    startScroll = None
    if G2frame.dataDisplay:
        startScroll = G2frame.dataDisplay.GetScrollPos(wx.VERTICAL) # save scroll position
        G2frame.dataDisplay.Destroy()
    else:
        CleanupMasks(data) # posting page for 1st time; clean out anything unfinished
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.MaskMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyMask, id=G2gd.wxID_MASKCOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadMask, id=G2gd.wxID_MASKLOAD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadMask, id=G2gd.wxID_MASKLOADNOT)
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
        G2frame.dataFrame.GetStatusBar().SetStatusText("To add mask: press a,r,s,p or f on 2D image for arc/ring/spot/polygon/frame")
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)

    thresh = data['Thresholds']         #min/max intensity range
    Spots = data['Points']               #x,y,radius in mm
    Rings = data['Rings']               #radius, thickness
    Polygons = data['Polygons']         #3+ x,y pairs
    if 'Frames' not in data:
        data['Frames'] = []
    frame = data['Frames']             #3+ x,y pairs
    Arcs = data['Arcs']                 #radius, start/end azimuth, thickness
    
    littleSizer = wx.FlexGridSizer(0,3,0,5)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Lower/Upper limits '),0,WACV)
    Text = wx.TextCtrl(G2frame.dataDisplay,value=str(thresh[0][0]),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,WACV)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    Text = wx.TextCtrl(G2frame.dataDisplay,value=str(thresh[0][1]),style=wx.TE_READONLY)
    littleSizer.Add(Text,0,WACV)
    Text.SetBackgroundColour(VERY_LIGHT_GREY)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Lower/Upper thresholds '),0,WACV)
    lowerThreshold = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,loc=thresh[1],key=0,
                                           min=thresh[0][0],OnLeave=Replot,typeHint=int)
    littleSizer.Add(lowerThreshold,0,WACV)
    upperThreshold = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,loc=thresh[1],key=1,
                                           max=thresh[0][1],OnLeave=Replot,typeHint=int)
    littleSizer.Add(upperThreshold,0,WACV)
    mainSizer.Add(littleSizer,0,)
    if Spots:
        lbl = wx.StaticText(parent=G2frame.dataDisplay,label=' Spot masks')
        lbl.SetBackgroundColour(wx.Colour(200,200,210))
        mainSizer.Add(lbl,0,wx.EXPAND|wx.ALIGN_CENTER,0)
        littleSizer = wx.FlexGridSizer(0,3,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' position, mm'),0,WACV)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' diameter, mm'),0,WACV)
        littleSizer.Add((5,0),0)
        for i in range(len(Spots)):
            if Spots[i]:
                x,y,d = Spots[i]
                spotText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.2f,%.2f" % (x,y)),
                    style=wx.TE_READONLY)
                spotText.SetBackgroundColour(VERY_LIGHT_GREY)
                littleSizer.Add(spotText,0,WACV)
                spotText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                spotDiameter = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,loc=Spots[i],key=2,
                                           max=100.,OnLeave=Replot,nDig=[8,2])
                littleSizer.Add(spotDiameter,0,WACV)
                spotDelete = G2G.G2LoggedButton(G2frame.dataDisplay,label='delete?',
                                            locationcode='Delete+Points+'+str(i),
                                            handler=onDeleteMask)
                littleSizer.Add(spotDelete,0,WACV)
        mainSizer.Add(littleSizer,0,)
    if Rings:
        lbl = wx.StaticText(parent=G2frame.dataDisplay,label=' Ring masks')
        lbl.SetBackgroundColour(wx.Colour(200,200,210))
        mainSizer.Add(lbl,0,wx.EXPAND|wx.ALIGN_CENTER,0)
        littleSizer = wx.FlexGridSizer(0,3,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' 2-theta,deg'),0,WACV)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' thickness, deg'),0,WACV)
        littleSizer.Add((5,0),0)
        for i in range(len(Rings)):
            if Rings[i]:
                ringText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (Rings[i][0])),
                    style=wx.TE_READONLY)
                ringText.SetBackgroundColour(VERY_LIGHT_GREY)
                ringText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(ringText,0,WACV)
                ringThick = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,loc=Rings[i],key=1,
                                           min=0.001,max=1.,OnLeave=Replot,nDig=[8,3])
                littleSizer.Add(ringThick,0,WACV)
                ringDelete = G2G.G2LoggedButton(G2frame.dataDisplay,label='delete?',
                                            locationcode='Delete+Rings+'+str(i),
                                            handler=onDeleteMask)
                littleSizer.Add(ringDelete,0,WACV)
        mainSizer.Add(littleSizer,0,)
    if Arcs:
        lbl = wx.StaticText(parent=G2frame.dataDisplay,label=' Arc masks')
        lbl.SetBackgroundColour(wx.Colour(200,200,210))
        mainSizer.Add(lbl,0,wx.EXPAND|wx.ALIGN_CENTER,0)
        littleSizer = wx.FlexGridSizer(0,4,0,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' 2-theta,deg'),0,WACV)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' azimuth, deg'),0,WACV)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' thickness, deg'),0,WACV)
        littleSizer.Add((5,0),0)
        for i in range(len(Arcs)):
            if Arcs[i]:
                tth,azimuth,thick = Arcs[i]
                arcText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.3f" % (tth)),
                    style=wx.TE_READONLY)
                arcText.SetBackgroundColour(VERY_LIGHT_GREY)
                arcText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(arcText,0,WACV)
                azmText = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%d,%d" % (azimuth[0],azimuth[1])),
                    style=wx.TE_READONLY)
                azmText.SetBackgroundColour(VERY_LIGHT_GREY)
                azmText.Bind(wx.EVT_ENTER_WINDOW,OnTextMsg)
                littleSizer.Add(azmText,0,WACV)
                arcThick = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,loc=Arcs[i],key=2,
                                           min=0.001,max=20.,OnLeave=Replot,nDig=[8,3])
                littleSizer.Add(arcThick,0,WACV)
                arcDelete = G2G.G2LoggedButton(G2frame.dataDisplay,label='delete?',
                                            locationcode='Delete+Arcs+'+str(i),
                                            handler=onDeleteMask)
                littleSizer.Add(arcDelete,0,WACV)
        mainSizer.Add(littleSizer,0,)
    if Polygons:
        lbl = wx.StaticText(parent=G2frame.dataDisplay,
            label=' Polygon masks (on plot RB vertex drag to move,\nLB vertex drag to insert)')
        lbl.SetBackgroundColour(wx.Colour(200,200,210))
        mainSizer.Add(lbl,0,wx.EXPAND|wx.ALIGN_CENTER,0)
        littleSizer = wx.FlexGridSizer(0,2,0,5)
        for i in range(len(Polygons)):
            if Polygons[i]:
                polyList = []
                for x,y in Polygons[i]:
                    polyList.append("%.2f, %.2f"%(x,y))
                polyText = wx.ComboBox(G2frame.dataDisplay,value=polyList[0],choices=polyList,style=wx.CB_READONLY)
                littleSizer.Add(polyText,0,WACV)
                polyDelete = G2G.G2LoggedButton(G2frame.dataDisplay,label='delete?',
                                            locationcode='Delete+Polygons+'+str(i),
                                            handler=onDeleteMask)
                littleSizer.Add(polyDelete,0,WACV)
        mainSizer.Add(littleSizer,0,)
    if frame:
        lbl = wx.StaticText(parent=G2frame.dataDisplay,
            label=' Frame mask (on plot RB vertex drag to move,LB vertex drag to insert)')
        lbl.SetBackgroundColour(wx.Colour(200,200,210))
        mainSizer.Add(lbl,0,wx.EXPAND|wx.ALIGN_CENTER,0)
        littleSizer = wx.FlexGridSizer(0,2,0,5)
        frameList = []
        for x,y in frame:
            frameList.append("%.2f, %.2f"%(x,y))
        frameText = wx.ComboBox(G2frame.dataDisplay,value=frameList[0],choices=frameList,style=wx.CB_READONLY)
        littleSizer.Add(frameText,0,WACV)
        frameDelete = G2G.G2LoggedButton(G2frame.dataDisplay,label='delete?',
                                            locationcode='Delete+Frame',
                                            handler=onDeleteFrame)
        littleSizer.Add(frameDelete,0,WACV)
        mainSizer.Add(littleSizer,0,)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 50 # room for scrollbar & status msg
    Size[1] = min(Size[1],500)
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)    
    wx.Yield()
    if startScroll: # reset scroll to saved position
        G2frame.dataDisplay.Scroll(0,startScroll) # set to saved scroll position
        wx.Yield()

################################################################################
##### Stress/Strain
################################################################################

def UpdateStressStrain(G2frame,data):
    '''Shows and handles the controls on the "Stress/Strain"
    data tree entry
    '''
    
    def OnAppendDzero(event):
        data['d-zero'].append({'Dset':1.0,'Dcalc':0.0,'pixLimit':10,'cutoff':1.0,
            'ImxyObs':[[],[]],'ImtaObs':[[],[]],'ImtaCalc':[[],[]],'Emat':[1.0,1.0,1.0]})
        UpdateStressStrain(G2frame,data)
        
    def OnUpdateDzero(event):
        for item in data['d-zero']:
            if item['Dcalc']:   #skip unrefined ones
                item['Dset'] = item['Dcalc']
        UpdateStressStrain(G2frame,data)
            
    def OnCopyStrSta(event):
        TextList = [[False,'All IMG',0,0]]
#        Names = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                name = G2frame.PatternTree.GetItemText(id)
#                Names.append(name)
                if 'IMG' in name:
                    Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Stress/Strain'))
                    if id == G2frame.Image:
                        Source = name
                    else:
                        TextList.append([False,name,id,Data.get('Sample load',0.0)])
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
                        ifcopy,name,id,load = item
                        if ifcopy:
                            Data = copy.deepcopy(data)
                            Data['Sample load'] = load
                            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'Stress/Strain'),Data)
            finally:
                dlg.Destroy()

    def OnLoadStrSta(event):
        dlg = wx.FileDialog(G2frame, 'Choose stress/strain file', '.', '', 
            'image control files (*.strsta)|*.strsta',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.read()
                data = eval(S)
                Controls = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
                G2img.FitStrSta(G2frame.ImageZ,data,Controls)
                UpdateStressStrain(G2frame,data)
                G2plt.PlotExposedImage(G2frame,event=event)
                G2plt.PlotStrain(G2frame,data,newPlot=True)
                File.close()
        finally:
            dlg.Destroy()

    def OnSaveStrSta(event):
        dlg = wx.FileDialog(G2frame, 'Choose stress/strain file', '.', '', 
            'image control files (*.strsta)|*.strsta',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'w')
                save = {}
                keys = ['Type','Sample phi','Sample z','Sample load']
                keys2 = ['Dset','Dcalc','pixLimit','cutoff','Emat']
                File.write('{\n\t')
                for key in keys:
                    if key in 'Type':
                        File.write("'"+key+"':'"+data[key]+"',")
                    else:
                        File.write("'"+key+"':"+str(data[key])+',')
                File.write('\n\t'+"'d-zero':[\n")
                for data2 in data['d-zero']:
                    File.write('\t\t{')
                    for key in keys2:
                        File.write("'"+key+"':"+str(data2[key])+',')
                    File.write("'ImxyObs':[[],[]],'ImtaObs':[[],[]],'ImtaCalc':[[],[]]},\n")
                File.write('\t]\n}')
                File.close()
        finally:
            dlg.Destroy()
            
    def OnStrStaSample(event):
        filename = ''
        dlg = wx.FileDialog(G2frame, 'Choose multihistogram metadata text file', '.', '', 
            'metadata file (*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                newItems = []
                itemNames = []
                Comments = []
                while S:
                    if S[0] == '#':
                        Comments.append(S)
                        S = File.readline()
                        continue
                    S = S.replace(',',' ').replace('\t',' ')
                    Stuff = S[:-1].split()
                    itemNames.append(Stuff[0])
                    newItems.append(Stuff[1:])
                    S = File.readline()                
                File.close()
        finally:
            dlg.Destroy()
        if not filename:
            G2frame.ErrorDialog('Nothing to do','No file selected')
            return
        dataDict = dict(zip(itemNames,newItems))
        ifany = False
        Names = [' ','Sample phi','Sample z','Sample load']
        dlg = G2G.G2ColumnIDDialog( G2frame,' Choose multihistogram metadata columns:',
            'Select columns',Comments,Names,np.array(newItems).T)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                colNames,newData = dlg.GetSelection()
                dataDict = dict(zip(itemNames,newData.T))
                for item in colNames:
                    if item != ' ':
                        ifany = True
        finally:
            dlg.Destroy()
        if not ifany:
            G2frame.ErrorDialog('Nothing to do','No columns identified')
            return
        histList = []
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)        
        while item:
            name = G2frame.PatternTree.GetItemText(item)
            if name.startswith('IMG'):
                histList.append(name)
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        colIds = {}
        for i,name in enumerate(colNames):
            if name != ' ':
                colIds[name] = i
        for hist in histList:
            name = hist.split()[1]  #this is file name
            if name in dataDict:
                newItems = {}
                for item in colIds:
                    newItems[item] = float(dataDict[name][colIds[item]])
                Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,hist)
                stsrData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Stress/Strain'))
                stsrData.update(newItems)        
        UpdateStressStrain(G2frame,data)        
    
    def OnFitStrSta(event):
        Controls = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
        G2img.FitStrSta(G2frame.ImageZ,data,Controls)
        print 'Strain fitting finished'
        UpdateStressStrain(G2frame,data)
        G2plt.PlotExposedImage(G2frame,event=event)
        G2plt.PlotStrain(G2frame,data,newPlot=True)
        
    def OnFitAllStrSta(event):
        TextList = [[False,'All IMG',0]]
#        Names = []
        if G2frame.PatternTree.GetCount():
            choices = G2gd.GetPatternTreeDataNames(G2frame,['IMG ',])
            if len(choices) == 1:
                G2frame.ErrorDialog('Nothing to fit','There must some "IMG" patterns')
                return
            sel = []
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Stress/Strain fitting','Select images to fit:',choices)
            dlg.SetSelections(sel)
            names = []
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    names.append(choices[sel])
            dlg.Destroy()
            SeqResult = {}
            dlg = wx.ProgressDialog('Sequential IMG Strain fit','Data set name = '+names[0],len(names), 
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)          
            wx.BeginBusyCursor()
            goodnames = []
            try:
                for i,name in enumerate(names):
                    print ' Sequential strain fit for ',name
                    GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                    if not GoOn:
                        break
                    Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name)
                    Controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Image Controls'))
                    StaCtrls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Stress/Strain'))
                    if not len(StaCtrls['d-zero']):
                        continue
                    goodnames.append(name)
                    id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root, name)
                    Npix,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(Id)
                    image = G2IO.GetImageData(G2frame,imagefile,True,ImageTag=imagetag)
                    dark = Controls['dark image']
                    if dark[0]:
                        id = G2gd.GetPatternTreeItemId(G2frame, G2frame.root,dark[0])
                        Npix,darkfile,imagetag = G2frame.PatternTree.GetImageLoc(id)
                        darkImg = G2IO.GetImageData(G2frame,darkfile,True,ImageTag=imagetag)
                        image += dark[1]*darkImg
                    G2img.FitStrSta(image,StaCtrls,Controls)
                    G2plt.PlotStrain(G2frame,StaCtrls,newPlot=True)
                    parmDict = {'Sample load':StaCtrls['Sample load'],}
                    varyNames = ['e11','e12','e22']
                    sig = []
                    varyList = []
                    variables = []
                    for i,item in enumerate(StaCtrls['d-zero']):
                        variables += item['Emat']
                        sig += item['Esig']
                        varylist = ['%d%s%s'%(i,';',Name) for Name in varyNames]
                        varyList += varylist
                        parmDict.update(dict(zip(varylist,item['Emat'])))
                        parmDict['%d:Dcalc'%(i)] = item['Dcalc']
                    SeqResult[name] = {'variables':variables,'varyList':varyList,'sig':sig,'Rvals':[],
                        'covMatrix':np.eye(len(variables)),'title':name,'parmDict':parmDict}
                else:
                    SeqResult['histNames'] = goodnames
                    dlg.Destroy()
                    print ' ***** Sequential strain refinement successful *****'
            finally:
                wx.EndBusyCursor()    
            Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequential results')
            if Id:
                G2frame.PatternTree.SetItemPyData(Id,SeqResult)
            else:
                Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Sequential results')
                G2frame.PatternTree.SetItemPyData(Id,SeqResult)
            G2frame.PatternTree.SelectItem(Id)
            print 'All images fitted'
        
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
                
        def OnSamLoad(event):
            try:
                value = float(samLoad.GetValue())
            except ValueError:
                value = data['Sample load']
            data['Sample load'] = value
            samLoad.SetValue("%.3f" % (data['Sample load']))
                
        samSizer = wx.BoxSizer(wx.HORIZONTAL)
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Strain type: '),0,WACV)
        strType = wx.ComboBox(G2frame.dataDisplay,value=data['Type'],choices=['True','Conventional'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        strType.SetValue(data['Type'])
        strType.Bind(wx.EVT_COMBOBOX, OnStrainType)
        samSizer.Add(strType,0,WACV)
        
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Sample phi: '),0,WACV)
        samPhi = wx.TextCtrl(G2frame.dataDisplay,-1,value=("%.3f" % (data['Sample phi'])),
            style=wx.TE_PROCESS_ENTER)
        samSizer.Add(samPhi,0,WACV)
        samPhi.Bind(wx.EVT_TEXT_ENTER,OnSamPhi)
        samPhi.Bind(wx.EVT_KILL_FOCUS,OnSamPhi)
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Sample delta-z(mm): '),0,WACV)
        samZ = wx.TextCtrl(G2frame.dataDisplay,-1,value=("%.3f" % (data['Sample z'])),
            style=wx.TE_PROCESS_ENTER)
        samSizer.Add(samZ,0,WACV)
        samZ.Bind(wx.EVT_TEXT_ENTER,OnSamZ)
        samZ.Bind(wx.EVT_KILL_FOCUS,OnSamZ)
        samSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=' Sample load(MPa): '),0,WACV)
        samLoad = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Sample load',
                nDig=[8,3],typeHint=float,)
        samSizer.Add(samLoad,0,WACV)

        return samSizer
        
    def DzeroSizer():
    
        def OnDzero(event):
            Obj = event.GetEventObject()
            try:
                value = min(20.0,max(0.25,float(Obj.GetValue())))
            except ValueError:
                value = 1.0
            Obj.SetValue("%.5f"%(value))
            data['d-zero'][Indx[Obj.GetId()]]['Dset'] = value
            data['d-zero'] = G2mth.sortArray(data['d-zero'],'Dset',reverse=True)
            Ring,R = G2img.MakeStrStaRing(data['d-zero'][Indx[Obj.GetId()]],G2frame.ImageZ,Controls)
            if len(Ring):
                data['d-zero'][Indx[Obj.GetId()]].update(R)
            else:
                G2frame.ErrorDialog('Strain peak selection','WARNING - No points found for this ring selection')
                
            UpdateStressStrain(G2frame,data)
            G2plt.PlotExposedImage(G2frame,event=event,newPlot=False)
            G2plt.PlotStrain(G2frame,data,newPlot=True)
            
        def OnDeleteDzero(event):
            Obj = event.GetEventObject()
            del(data['d-zero'][delIndx.index(Obj)])
            UpdateStressStrain(G2frame,data)
            G2plt.PlotExposedImage(G2frame,event=event,newPlot=True)
            G2plt.PlotStrain(G2frame,data,newPlot=True)
        
        def OnCutOff(event):
            Obj = event.GetEventObject()
            try:
                value = min(10.0,max(0.5,float(Obj.GetValue())))
            except ValueError:
                value = 1.0
            Obj.SetValue("%.1f"%(value))
            data['d-zero'][Indx[Obj.GetId()]]['cutoff'] = value 
            Ring,R = G2img.MakeStrStaRing(data['d-zero'][Indx[Obj.GetId()]],G2frame.ImageZ,Controls)
            G2plt.PlotExposedImage(G2frame,event=event)
            G2plt.PlotStrain(G2frame,data,newPlot=True)
        
        def OnPixLimit(event):
            Obj = event.GetEventObject()
            data['d-zero'][Indx[Obj.GetId()]]['pixLimit'] = int(Obj.GetValue())
            Ring,R = G2img.MakeStrStaRing(data['d-zero'][Indx[Obj.GetId()]],G2frame.ImageZ,Controls)
            G2plt.PlotExposedImage(G2frame,event=event)
            G2plt.PlotStrain(G2frame,data,newPlot=True)
            
        Indx = {}
        delIndx = []    
        dzeroSizer = wx.FlexGridSizer(0,8,5,5)
        for id,dzero in enumerate(data['d-zero']):
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' d-zero #%d: '%(id))),0,WACV)
            dZero = wx.TextCtrl(G2frame.dataDisplay,-1,value=('%.5f'%(dzero['Dset'])),
                style=wx.TE_PROCESS_ENTER)
            dzeroSizer.Add(dZero,0,WACV)
            dZero.Bind(wx.EVT_TEXT_ENTER,OnDzero)
            dZero.Bind(wx.EVT_KILL_FOCUS,OnDzero)
            Indx[dZero.GetId()] = id
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' d-zero ave: %.5f'%(dzero['Dcalc']))),0,WACV)
                
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min ring I/Ib '),0,WACV)
            cutOff = wx.TextCtrl(parent=G2frame.dataDisplay,value=("%.1f" % (dzero['cutoff'])),
                style=wx.TE_PROCESS_ENTER)
            cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)
            cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
            Indx[cutOff.GetId()] = id
            dzeroSizer.Add(cutOff,0,WACV)
        
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Pixel search range '),0,WACV)
            pixLimit = wx.ComboBox(parent=G2frame.dataDisplay,value=str(dzero['pixLimit']),choices=['1','2','5','10','15','20'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            pixLimit.Bind(wx.EVT_COMBOBOX, OnPixLimit)
            Indx[pixLimit.GetId()] = id
            dzeroSizer.Add(pixLimit,0,WACV)                
                
            dzeroDelete = wx.CheckBox(parent=G2frame.dataDisplay,label='delete?')
            dzeroDelete.Bind(wx.EVT_CHECKBOX,OnDeleteDzero)
            delIndx.append(dzeroDelete)
            dzeroSizer.Add(dzeroDelete,0,WACV)
            
            dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=(' Strain tensor:')),WACV)
            names = ['e11','e12','e22']
            for i in range(3):
                dzeroSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label=names[i]),0,WACV)
                tensorElem = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2f'%(dzero['Emat'][i]),style=wx.TE_READONLY)
                tensorElem.SetBackgroundColour(VERY_LIGHT_GREY)
                dzeroSizer.Add(tensorElem,0,WACV)
            dzeroSizer.Add((5,5),0)             
        return dzeroSizer
        
# patches
    if 'Sample load' not in data:
        data['Sample load'] = 0.0
# end patches
    
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    Controls = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))        
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.StrStaMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAppendDzero, id=G2gd.wxID_APPENDDZERO)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnUpdateDzero, id=G2gd.wxID_UPDATEDZERO)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitStrSta, id=G2gd.wxID_STRSTAFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitAllStrSta, id=G2gd.wxID_STRSTAALLFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyStrSta, id=G2gd.wxID_STRSTACOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadStrSta, id=G2gd.wxID_STRSTALOAD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveStrSta, id=G2gd.wxID_STRSTASAVE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnStrStaSample, id=G2gd.wxID_STRSTSAMPLE)        
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    if G2frame.StrainKey == 'a':    #probably doesn't happen
        G2frame.dataFrame.GetStatusBar().SetStatusText('Add strain ring active - LB pick d-zero value')
    else:
        G2frame.dataFrame.GetStatusBar().SetStatusText("To add strain data: On 2D Powder Image, key a:add ring")
        
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    mainSizer.Add(SamSizer())
    mainSizer.Add((5,10),0)
    mainSizer.Add(DzeroSizer())
    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)    

###########################################################################
# Autointegration follows 
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
    CleanupMasks(save)
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
            self.Pause = False
            # change button label
            if btnstart.GetLabel() != 'Pause':
                btnstart.SetLabel('Pause')
                if self.timer.IsRunning(): self.timer.Stop()
                self.PreventReEntryTimer = False
                self.StartLoop()
                self.OnTimerLoop(None) # run once immediately and again after delay
                self.timer.Start(int(1000*PollTime),oneShot=False)
                self.Status.SetStatusText('Press Pause to delay integration or Reset to prepare to reintegrate all images')
            else:
                btnstart.SetLabel('Resume')
                if self.timer.IsRunning(): self.timer.Stop()
                print('\nPausing autointegration\n')
                self.Status.SetStatusText('Press Resume to continue integration or Reset to prepare to reintegrate all images')
                self.Pause = True

        def OnReset(event):
            '''Called when Reset button is pressed. This stops the
            processing loop and resets the list of integrated files so
            all images can be reintegrated. 
            '''
            btnstart.SetLabel('Restart')
            self.Status.SetStatusText('Press Restart to reload and re-integrate images matching filter')
            if self.timer.IsRunning(): self.timer.Stop()
            self.Reset = True
            self.Pause = True
            
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
            mode, get distance-dependent parameters from user. 
            '''
            self.Evaluator = None
            if r2.GetValue():
                try:
                    dlg = IntegParmTable(self.G2frame) # create the dialog
                    dlg.CenterOnParent()
                    if dlg.ShowModal() == wx.ID_OK:
                        self.ImgTblParms = dlg.parms
                        self.IMfileList = dlg.IMfileList
                        self.Evaluator = DefineEvaluator(dlg)
                        self.params['Mode'] = 'table'
                        r2E.Enable(True)
                    else:
                        self.useActive.SetValue(True)
                finally:
                    dlg.Destroy()
            if self.useActive.GetValue():
                self.params['Mode'] = 'active'
                self.imageBase = G2frame.Image
                self.useActive.SetLabel("Active Image: "+
                        G2frame.PatternTree.GetItemText(self.imageBase))
                r2E.Enable(False)

        def OnEditTable(event):
            '''Called to edit the distance-dependent parameter look-up table.
            Should be called only when table is defined and active.
            '''
            try:
                dlg = IntegParmTable(self.G2frame,self.ImgTblParms,self.IMfileList)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    self.ImgTblParms = dlg.parms
                    self.IMfileList = dlg.IMfileList
                    self.Evaluator = DefineEvaluator(dlg)
                    self.params['Mode'] = 'table'
                    r2E.Enable(True)
                else:
                    self.useActive.SetValue(True)
                    self.params['Mode'] = 'active'
                    self.imageBase = G2frame.Image
                    self.useActive.SetLabel("Active Image: "+
                            G2frame.PatternTree.GetItemText(self.imageBase))
                    r2E.Enable(False)
            finally:
                dlg.Destroy()

        ##################################################
        # beginning of __init__ processing
        ##################################################
        self.G2frame = G2frame
        self.ImgTblParms = None
        self.IMfileList = None
        self.Evaluator = None
        self.params = {}
        self.Reset = False
        self.Pause = False
        self.PreventReEntryShowMatch = False
        self.PreventReEntryTimer = False
        self.params['IMGfile'] = ''
        self.params['MaskFile'] = ''
        self.params['IgnoreMask'] = True
        self.fmtlist = G2IO.ExportPowderList(G2frame)
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnTimerLoop)
        self.imageBase = G2frame.Image

        controlsId = G2frame.PatternTree.GetSelection()
        size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(self.imageBase)        
        self.imagedir,fileroot = os.path.split(imagefile)
        self.params['filter'] = '*'+os.path.splitext(fileroot)[1]
        self.params['outdir'] = os.path.abspath(self.imagedir)
        wx.Frame.__init__(self, G2frame, title='Automatic Integration',
                          style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.Status = self.CreateStatusBar()
        self.Status.SetStatusText('Press Start to load and integrate images matching filter')
        mnpnl = wx.Panel(self)
        mnsizer = wx.BoxSizer(wx.VERTICAL)
        # box for integration controls & masks input
        lbl = wx.StaticBox(mnpnl, wx.ID_ANY, "Integration Control")
        lblsizr = wx.StaticBoxSizer(lbl, wx.VERTICAL)
        lblsizr.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Use integration parameters from:'))
        self.useActive = wx.RadioButton(mnpnl, wx.ID_ANY, 
                            style = wx.RB_GROUP)
        self.useActive.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
        self.useActive.SetLabel("Active Image: "+
                                G2frame.PatternTree.GetItemText(self.imageBase))
        lblsizr.Add(self.useActive,1,wx.EXPAND,1)
        self.useActive.SetValue(True)
        minisizer = wx.BoxSizer(wx.HORIZONTAL)
        r2 = wx.RadioButton(mnpnl, wx.ID_ANY, "From distance look-up table")
        minisizer.Add(r2,0,wx.ALIGN_LEFT|wx.ALL,1)
        r2.Bind(wx.EVT_RADIOBUTTON, OnRadioSelect)
        r2E = wx.Button(mnpnl,  wx.ID_ANY, "Edit table")
        minisizer.Add(r2E,0,wx.ALIGN_LEFT,10)
        r2E.Enable(False)
        r2E.Bind(wx.EVT_BUTTON, OnEditTable)
        # bind button and deactivate be default
        lblsizr.Add(minisizer)
        mnsizer.Add(lblsizr,1,wx.EXPAND,1)

        # file filter stuff
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Image filter'))
        flterInp = G2G.ValidatedTxtCtrl(mnpnl,self.params,'filter',
                                        OnLeave=self.ShowMatchingFiles)
        sizer.Add(flterInp)
        mnsizer.Add(sizer,0,wx.ALIGN_RIGHT,1)
        self.ListBox = wx.ListBox(mnpnl,size=(-1,100))
        mnsizer.Add(self.ListBox,0,wx.EXPAND,1)
        self.ShowMatchingFiles(self.params['filter'])
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

    def ShowMatchingFiles(self,value,invalid=False,**kwargs):
        G2frame = self.G2frame
        if invalid: return
        msg = ''
        if self.PreventReEntryShowMatch: return
        self.PreventReEntryShowMatch = True
        imageFileList = []
        for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
            imgId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
            size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(imgId)
            if imagefile not in imageFileList: imageFileList.append(imagefile)
            if img not in G2frame.IntegratedList:
                if msg: msg += '\n'
                msg += '  ' + img
        if msg: msg = "Loaded images to integrate:\n" + msg + "\n"
        msg1 = ""
        try:
            imageList = sorted(
                glob.glob(os.path.join(self.imagedir,value)))
            if not imageList:
                msg1 = 'Warning: No files match search string '+os.path.join(self.imagedir,value)
            else:
                for fil in imageList:
                    if fil not in imageFileList: msg1 += '\n  '+fil
                if msg1:
                    msg += 'Files to integrate from '+os.path.join(self.imagedir,value)+msg1
                else:
                    msg += 'All files integrated'
        except IndexError:
            msg += 'Error searching for files named '+os.path.join(self.imagedir,value)
        self.ListBox.Clear()
        self.ListBox.AppendItems(msg.split('\n'))
        self.PreventReEntryShowMatch = False
        return
        
    def IntegrateImage(self,img):
        '''Integrates a single image'''
        G2frame = self.G2frame
        imgId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
        G2frame.Image = imgId
        G2frame.PickId = G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls')
        # do integration
        size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(imgId)
        G2frame.ImageZ = G2IO.GetImageData(G2frame,imagefile,True,imagetag)
        masks = G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
        data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
        # simulate a Image Controls press, since that is where the
        # integration is hidden
        UpdateImageControls(G2frame,data,masks,IntegrateOnly=True)
        G2frame.IntegratedList.append(img) # note this as integrated
        # split name and control number
        s = re.split(r'(\d+)\Z',os.path.split(os.path.splitext(imagefile)[0])[1])
        namepre = s[0]
        if len(s) > 1:
            namenum = s[1]
        else:
            namenum = ''
        for Id in G2frame.IntgOutList: # loop over newly created PDWR entry(ies)
            # save the created PWDR tree names so that a reset can delete them
            treename = G2frame.PatternTree.GetItemText(Id)
            G2frame.AutointPWDRnames.append(treename)
            # write out the images in the selected formats
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
                
    def ResetFromTable(self,dist):
        '''Sets integration parameters based on values from
        the lookup table
        '''
        #dist = self.controlsDict['distance']
        interpDict,imgctrl,immask = self.Evaluator(dist) # interpolated calibration values
        #if GSASIIpath.GetConfigValue('debug'):
        if GSASIIpath.GetConfigValue('debug'):
            print 'interpolated values: ',interpDict
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
        
    def StartLoop(self):
        '''Save current Image params for use in future integrations
        also label the window so users understand what is being used
        '''
        print '\nStarting new autointegration\n'
        G2frame = self.G2frame
        # show current IMG base
        if self.params['Mode'] != 'table':
            self.useActive.SetLabel("Active Image: "+
                                    G2frame.PatternTree.GetItemText(self.imageBase))
            # load copy of Image Controls from current image and clean up
            # items that should not be copied
            self.ImageControls = copy.deepcopy(
                G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(
                    G2frame,self.imageBase, 'Image Controls')))
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
                    G2gd.GetPatternTreeItemId(G2frame,self.imageBase, 'Masks')))
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
            self.G2frame.IntegratedList = []
            
            if self.params['Mode'] != 'table': # reset controls and masks for all IMG items in tree to master
                for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
                    # update controls from master
                    controlsDict = G2frame.PatternTree.GetItemPyData(
                        G2gd.GetPatternTreeItemId(G2frame,self.imageBase, 'Image Controls'))
                    controlsDict.update(self.ImageControls)
                    # update masks from master
                    ImageMasks = G2frame.PatternTree.GetItemPyData(
                        G2gd.GetPatternTreeItemId(G2frame,self.imageBase, 'Masks'))
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
            wx.Yield()
            self.Reset = False
        G2frame.AutointPWDRnames = [] # list of created PWDR tree item names

    def OnTimerLoop(self,event):
        '''A method that is called every :meth:`PollTime` seconds that is
        used to check for new files and process them. This is called only
        after the "Start" button is pressed (when its label reads "Pause").
        '''
        G2frame = self.G2frame
        try:
            self.currImageList = sorted(
                glob.glob(os.path.join(self.imagedir,self.params['filter'])))
            self.ShowMatchingFiles(self.params['filter'])
        except IndexError:
            self.currImageList = []
            return

        if self.PreventReEntryTimer: return
        self.PreventReEntryTimer = True
        imageFileList = []
        # integrate the images that have already been read in, but
        # have not yet been processed
        for img in G2gd.GetPatternTreeDataNames(G2frame,['IMG ']):
            imgId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,img)
            size,imagefile,imagetag = G2frame.PatternTree.GetImageLoc(imgId)
            # Create a list of image files that have been read in
            if imagefile not in imageFileList: imageFileList.append(imagefile)
            # skip if already integrated
            if img in G2frame.IntegratedList: continue
            if self.params['Mode'] == 'table': # look up parameter values from table
                controlsDict = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Image Controls'))
                ImageMasks = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Masks'))
                self.ResetFromTable(controlsDict['distance'])
                # update controls from master
                controlsDict.update(self.ImageControls)
                # update masks from master w/o Thresholds
                ImageMasks.update(self.ImageMasks)
            self.IntegrateImage(img)
            self.G2frame.oldImagefile = '' # mark image as changed; reread as needed
            wx.Yield()
            self.ShowMatchingFiles(self.params['filter'])
            wx.Yield()
            if self.Pause: return

        # loop over image files matching glob, reading in any new ones
        for newImage in self.currImageList:
            if newImage in imageFileList: continue # already read?
            for imgId in G2IO.ReadImages(G2frame,newImage):
                controlsDict = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Image Controls'))
                ImageMasks = G2frame.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(G2frame,imgId, 'Masks'))
                if self.params['Mode'] == 'table': # look up parameter values from table
                    self.ResetFromTable(controlsDict['distance'])
                # update controls from master
                controlsDict.update(self.ImageControls)
                # update masks from master w/o Thresholds
                ImageMasks.update(self.ImageMasks)
                # now integrate the image
                img = G2frame.PatternTree.GetItemText(imgId)
                self.IntegrateImage(img)
                self.G2frame.oldImagefile = '' # mark image as changed; reread as needed
                wx.Yield()
                self.ShowMatchingFiles(self.params['filter'])
                wx.Yield()
            if self.Pause: return
        
        if GSASIIpath.GetConfigValue('debug'):
            import datetime
            print ("Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))
        self.PreventReEntryTimer = False

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
    def __init__(self,G2frame,parms=None,IMfileList=None):
        self.G2frame = G2frame
        wx.Dialog.__init__(self,G2frame,style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE)
        if parms:
            self.parms = parms # list of values by column
            self.IMfileList = IMfileList # list of .imctrl file names for each entry in table
        else:
            self.parms = [] # list of values by column
            self.IMfileList = [] # list of .imctrl file names for each entry in table
            files = []
            try:
                dlg = wx.FileDialog(self, 'Select image control files or previous table', 
                                    style=wx.OPEN| wx.MULTIPLE,
                                    wildcard='image control files (.imctrl)|*.imctrl|Integration table (*.imtbl)|*.imtbl')
                dlg.CenterOnParent()
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
        btn = wx.Button(self, wx.ID_ANY,'Save as file')
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
                    [key,val] = S[:-1].split(':',1)
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
                        defaultDir=dir, defaultFile=f, style=wx.SAVE,
                        wildcard='G2 Image Param Table file (*.imtbl)|*.imtbl')
            dlg.CenterOnParent()
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
            dlg = wx.FileDialog(self, 'Select mask or control file to add (Press cancel if none)', 
                                style=wx.OPEN,
                                wildcard='Add GSAS-II mask file (.immask)|*.immask|add image control file (.imctrl)|*.imctrl')
            dlg.CenterOnParent()
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
# Autointegration end
###########################################################################
