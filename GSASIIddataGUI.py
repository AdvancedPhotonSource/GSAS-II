# -*- coding: utf-8 -*-
#GSASII - phase data display routines
########### SVN repository information ###################
# $Date: 2013-02-01 15:23:56 -0600 (Fri, 01 Feb 2013) $
# $Author: vondreele $
# $Revision: 844 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIIddataGUI.py $
# $Id: GSASIIddataGUI.py 844 2013-02-01 21:23:56Z vondreele $
########### SVN repository information ###################
'''
*GSASIIddataGUI: Phase Diffraction Data GUI*
--------------------------------------------

Module to create the GUI for display of diffraction data * phase
information that is shown in the data display window
(when a phase is selected.)

'''
import wx
import wx.grid as wg
import wx.lib.gridmovers as wgmove
import matplotlib as mpl
import math
import copy
import time
import sys
import random as ran
import cPickle
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 844 $")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIIO as G2IO
import GSASIIstruct as G2str
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma

VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)
BLACK = wx.Colour(0,0,0)
mapDefault = {'MapType':'','RefList':'','Resolution':0.5,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}
# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi

################################################################################
##### DData routines
################################################################################        
def UpdateDData(G2frame,DData,data):
    '''Display the Diffraction Data associated with a phase
    (items where there is a value for each histogram and phase)

    :param wx.frame G2frame: the main GSAS-II frame object

    :param wx.ScrolledWindow DData: notebook page to be used for the display

    :param dict data: all the information on the phase in a dictionary

    '''
    G2frame.dataFrame.SetStatusText('')
    UseList = data['Histograms']
    if UseList:
        G2frame.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,True)
        for item in G2frame.Refine: item.Enable(True)
    else:
        G2frame.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,False)
        for item in G2frame.Refine: item.Enable(False)
    generalData = data['General']
    PhaseName = generalData['Name']       
    SGData = generalData['SGData']
    keyList = UseList.keys()
    keyList.sort()
    PWDR = any(['PWDR' in item for item in keyList])
    Indx = {}
    
    def PlotSizer():

        def OnPlotSel(event):
            Obj = event.GetEventObject()
            generalData['Data plot type'] = Obj.GetStringSelection()
            wx.CallAfter(UpdateDData,G2frame,DData,data)
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnPOhkl(event):
            Obj = event.GetEventObject()
            Saxis = Obj.GetValue().split()
            try:
                hkl = [int(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                hkl = generalData['POhkl']
            if not np.any(np.array(hkl)):
                hkl = generalData['POhkl']
            generalData['POhkl'] = hkl
            h,k,l = hkl
            Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
            G2plt.PlotSizeStrainPO(G2frame,data)
        
        plotSizer = wx.BoxSizer(wx.VERTICAL)
        choice = ['None','Mustrain','Size','Preferred orientation']
        plotSel = wx.RadioBox(DData,-1,'Select plot type:',choices=choice,
            majorDimension=2,style=wx.RA_SPECIFY_COLS)
        plotSel.SetStringSelection(generalData['Data plot type'])
        plotSel.Bind(wx.EVT_RADIOBOX,OnPlotSel)    
        plotSizer.Add(plotSel)
        if generalData['Data plot type'] == 'Preferred orientation':
            POhklSizer = wx.BoxSizer(wx.HORIZONTAL)
            POhklSizer.Add(wx.StaticText(DData,-1,' Plot preferred orientation for H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
            h,k,l = generalData['POhkl']
            poAxis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
            poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOhkl)
            poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOhkl)
            POhklSizer.Add(poAxis,0,wx.ALIGN_CENTER_VERTICAL)
            plotSizer.Add(POhklSizer)            
        return plotSizer
       
    def ScaleSizer():
        
        def OnScaleRef(event):
            Obj = event.GetEventObject()
            UseList[Indx[Obj.GetId()]]['Scale'][1] = Obj.GetValue()
            
        def OnScaleVal(event):
            Obj = event.GetEventObject()
            try:
                scale = float(Obj.GetValue())
                if scale > 0:
                    UseList[Indx[Obj.GetId()]]['Scale'][0] = scale
            except ValueError:
                pass
            Obj.SetValue("%.4f"%(UseList[Indx[Obj.GetId()]]['Scale'][0]))          #reset in case of error
                        
        scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
        if 'PWDR' in item:
            scaleRef = wx.CheckBox(DData,-1,label=' Phase fraction: ')
        elif 'HKLF' in item:
            scaleRef = wx.CheckBox(DData,-1,label=' Scale factor: ')                
        scaleRef.SetValue(UseList[item]['Scale'][1])
        Indx[scaleRef.GetId()] = item
        scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        scaleSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
        scaleVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.4f'%(UseList[item]['Scale'][0]),style=wx.TE_PROCESS_ENTER)
        Indx[scaleVal.GetId()] = item
        scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
        scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
        scaleSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
        return scaleSizer
        
    def OnUseData(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        UseList[hist]['Use'] = Obj.GetValue()
        
    def OnShowData(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        UseList[hist]['Show'] = Obj.GetValue()
        wx.CallAfter(UpdateDData,G2frame,DData,data)
        G2plt.PlotSizeStrainPO(G2frame,data)
        
    def OnCopyData(event):
        #how about HKLF data? This is only for PWDR data
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        sourceDict = UseList[hist]
        copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet']
        copyDict = {}
        for name in copyNames: 
            copyDict[name] = copy.deepcopy(sourceDict[name])        #force copy
        keyList = ['All',]+UseList.keys()
        if UseList:
            copyList = []
            dlg = wx.MultiChoiceDialog(G2frame, 
                'Copy parameters to which histograms?', 'Copy parameters', 
                keyList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: 
                        copyList.append(keyList[i])
                    if 'All' in copyList: 
                        copyList = keyList[1:]
                    for item in copyList:
                        UseList[item].update(copy.deepcopy(copyDict))
                    wx.CallAfter(UpdateDData,G2frame,DData,data)
            finally:
                dlg.Destroy()
                
    def OnCopyFlags(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        sourceDict = UseList[hist]
        copyDict = {}
        copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet']
        babNames = ['BabA','BabU']
        for name in copyNames:
            if name in ['Scale','Extinction','HStrain']:
                copyDict[name] = sourceDict[name][1]
            elif name in ['Size','Mustrain']:
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2],sourceDict[name][4]]
            elif name == 'Pref.Ori.':
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2]]
                if sourceDict[name][0] == 'SH':
                    SHterms = sourceDict[name][5]
                    SHflags = {}
                    for item in SHterms:
                        SHflags[item] = SHterms[item][1]
                    copyDict[name].append(SHflags)
            elif name == 'Babinet':
                copyDict[name] = {}
                for bab in babNames:
                    copyDict[name][bab] = sourceDict[name][bab][1]                       
        keyList = ['All',]+UseList.keys()
        if UseList:
            copyList = []
            dlg = wx.MultiChoiceDialog(G2frame, 
                'Copy parameters to which histograms?', 'Copy parameters', 
                keyList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: 
                        copyList.append(keyList[i])
                    if 'All' in copyList: 
                        copyList = keyList[1:]
                    for item in copyList:
                        UseList[item]                            
                        for name in copyNames:
                            if name in ['Scale','Extinction','HStrain']:
                                UseList[item][name][1] = copy.copy(copyDict[name])
                            elif name in ['Size','Mustrain']:
                                UseList[item][name][0] = copy.copy(copyDict[name][0])
                                UseList[item][name][2] = copy.copy(copyDict[name][1])
                                UseList[item][name][4] = copy.copy(copyDict[name][2])
                            elif name == 'Pref.Ori.':
                                UseList[item][name][0] = copy.copy(copyDict[name][0])
                                UseList[item][name][2] = copy.copy(copyDict[name][1])
                                if sourceDict[name][0] == 'SH':
                                    SHflags = copy.copy(copyDict[name][2])
                                    SHterms = copy.copy(sourceDict[name][5])
                                    for item in SHflags:
                                        SHterms[item][1] = copy.copy(SHflags[item])
                            elif name == 'Babinet':
                                for bab in babNames:
                                    UseList[item][name][bab][1] = copy.copy(copyDict[name][bab])                                              
                    wx.CallAfter(UpdateDData,G2frame,DData,data)
            finally:
                dlg.Destroy()
        
        
    def OnLGmixRef(event):
        Obj = event.GetEventObject()
        hist,name = Indx[Obj.GetId()]
        UseList[hist][name][2][2] = Obj.GetValue()
        
    def OnLGmixVal(event):
        Obj = event.GetEventObject()
        hist,name = Indx[Obj.GetId()]
        try:
            value = float(Obj.GetValue())
            if 0.1 <= value <= 1:
                UseList[hist][name][1][2] = value
            else:
                raise ValueError
        except ValueError:
            pass
        Obj.SetValue("%.4f"%(UseList[hist][name][1][2]))          #reset in case of error

    def OnSizeType(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        UseList[hist]['Size'][0] = Obj.GetValue()
        G2plt.PlotSizeStrainPO(G2frame,data)
        wx.CallAfter(UpdateDData,G2frame,DData,data)
        
    def OnSizeRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[hist]['Size'][0] == 'ellipsoidal':
            UseList[hist]['Size'][5][pid] = Obj.GetValue()                
        else:
            UseList[hist]['Size'][2][pid] = Obj.GetValue()
        
    def OnSizeVal(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[hist]['Size'][0] == 'ellipsoidal':
            try:
                size = float(Obj.GetValue())
                if pid < 3 and size <= 0.001:            #10A lower limit!
                    raise ValueError
                UseList[hist]['Size'][4][pid] = size                    
            except ValueError:
                pass
            Obj.SetValue("%.3f"%(UseList[hist]['Size'][4][pid]))          #reset in case of error
        else:
            try:
                size = float(Obj.GetValue())
                if size <= 0.001:            #10A lower limit!
                    raise ValueError
                UseList[hist]['Size'][1][pid] = size
            except ValueError:
                pass
            Obj.SetValue("%.3f"%(UseList[hist]['Size'][1][pid]))          #reset in case of error
        G2plt.PlotSizeStrainPO(G2frame,data)
        
    def OnSizeAxis(event):            
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[hist]['Size'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[hist]['Size'][3]
        UseList[hist]['Size'][3] = hkl
        h,k,l = hkl
        Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
                    
    def OnResetSize(event):
        Obj = event.GetEventObject()
        Obj.SetValue(False)
        item,name = Indx[Obj.GetId()]
        if name == 'isotropic':
            UseList[item]['Size'][1][0] = 1.0
        elif name == 'uniaxial':
            UseList[item]['Size'][1][0] = 1.0
            UseList[item]['Size'][1][1] = 1.0
        elif name == 'ellipsoidal':
            for i in range(3):
                UseList[item]['Size'][4][i] = 1.0
                UseList[item]['Size'][4][i+3] = 0.0
        G2plt.PlotSizeStrainPO(G2frame,data)
        wx.CallAfter(UpdateDData,G2frame,DData,data)
            
    def OnStrainType(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        UseList[hist]['Mustrain'][0] = Obj.GetValue()
        wx.CallAfter(UpdateDData,G2frame,DData,data)
        G2plt.PlotSizeStrainPO(G2frame,data)
        
    def OnStrainRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[hist]['Mustrain'][0] == 'generalized':
            UseList[hist]['Mustrain'][5][pid] = Obj.GetValue()
        else:
            UseList[hist]['Mustrain'][2][pid] = Obj.GetValue()
        
    def OnStrainVal(event):
        Snames = G2spc.MustrainNames(SGData)
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        try:
            strain = float(Obj.GetValue())
            if UseList[hist]['Mustrain'][0] == 'generalized':
                if '4' in Snames[pid] and strain < 0:
                    raise ValueError
                UseList[hist]['Mustrain'][4][pid] = strain
            else:
                if strain <= 0:
                    raise ValueError
                UseList[hist]['Mustrain'][1][pid] = strain
        except ValueError:
            pass
        if UseList[hist]['Mustrain'][0] == 'generalized':
            Obj.SetValue("%.3f"%(UseList[hist]['Mustrain'][4][pid]))          #reset in case of error
        else:
            Obj.SetValue("%.1f"%(UseList[hist]['Mustrain'][1][pid]))          #reset in case of error
        G2plt.PlotSizeStrainPO(G2frame,data)
        
    def OnStrainAxis(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[hist]['Mustrain'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[hist]['Mustrain'][3]
        UseList[hist]['Mustrain'][3] = hkl
        h,k,l = hkl
        Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
        G2plt.PlotSizeStrainPO(G2frame,data)
        
    def OnResetStrain(event):
        Obj = event.GetEventObject()
        Obj.SetValue(False)
        item,name = Indx[Obj.GetId()]
        if name == 'isotropic':
            UseList[item]['Mustrain'][1][0] = 1000.0
        elif name == 'uniaxial':
            UseList[item]['Mustrain'][1][0] = 1000.0
            UseList[item]['Mustrain'][1][1] = 1000.0
        elif name == 'generalized':
            nTerm = len(UseList[item]['Mustrain'][4])
            for i in range(nTerm):
                UseList[item]['Mustrain'][4][i] = 0.01
        wx.CallAfter(UpdateDData,G2frame,DData,data)
        G2plt.PlotSizeStrainPO(G2frame,data)
            
    def OnHstrainRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        UseList[hist]['HStrain'][1][pid] = Obj.GetValue()
        
    def OnHstrainVal(event):
        Snames = G2spc.HStrainNames(SGData)
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        try:
            strain = float(Obj.GetValue())
            UseList[hist]['HStrain'][0][pid] = strain
        except ValueError:
            pass
        Obj.SetValue("%.5f"%(UseList[hist]['HStrain'][0][pid]))          #reset in case of error

    def OnPOVal(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        try:
            mdVal = float(Obj.GetValue())
            if mdVal > 0:
                UseList[hist]['Pref.Ori.'][1] = mdVal
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[hist]['Pref.Ori.'][1]))          #reset in case of error
        
    def OnPOAxis(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[hist]['Pref.Ori.'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[hist]['Pref.Ori.'][3]
        UseList[hist]['Pref.Ori.'][3] = hkl
        h,k,l = hkl
        Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
        
    def OnPOOrder(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        Order = int(Obj.GetValue())
        UseList[hist]['Pref.Ori.'][4] = Order
        UseList[hist]['Pref.Ori.'][5] = SetPOCoef(Order,hist)
        wx.CallAfter(UpdateDData,G2frame,DData,data)

    def OnPOType(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        if 'March' in Obj.GetValue():
            UseList[hist]['Pref.Ori.'][0] = 'MD'
        else:
            UseList[hist]['Pref.Ori.'][0] = 'SH'
        wx.CallAfter(UpdateDData,G2frame,DData,data)            

    def OnPORef(event):
        Obj = event.GetEventObject()
        hist = Indx[Obj.GetId()]
        UseList[hist]['Pref.Ori.'][2] = Obj.GetValue()
            
    def SetPOCoef(Order,hist):
        cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',Order,False)     #cylindrical & no M
        newPOCoef = dict(zip(cofNames,np.zeros(len(cofNames))))
        POCoeff = UseList[hist]['Pref.Ori.'][5]
        for cofName in POCoeff:
            if cofName in  cofNames:
                newPOCoef[cofName] = POCoeff[cofName]
        return newPOCoef
    
    def OnExtRef(event):
        Obj = event.GetEventObject()
        UseList[Indx[Obj.GetId()]]['Extinction'][1] = Obj.GetValue()
        
    def OnExtVal(event):
        Obj = event.GetEventObject()
        try:
            ext = float(Obj.GetValue())
            if ext >= 0:
                UseList[Indx[Obj.GetId()]]['Extinction'][0] = ext
        except ValueError:
            pass
        Obj.SetValue("%.2f"%(UseList[Indx[Obj.GetId()]]['Extinction'][0]))

    def OnBabRef(event):
        Obj = event.GetEventObject()
        item,bab = Indx[Obj.GetId()]
        UseList[item]['Babinet']['Bab'+bab][1] = Obj.GetValue()
        
    def OnBabVal(event):
        Obj = event.GetEventObject()
        item,bab = Indx[Obj.GetId()]
        try:
            val = float(Obj.GetValue())
            if val >= 0:
                UseList[item]['Babinet']['Bab'+bab][0] = val
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[item]['Babinet']['Bab'+bab][0]))

    def OnTbarVal(event):
        Obj = event.GetEventObject()
        try:
            tbar = float(Obj.GetValue())
            if tbar > 0:
                UseList[Indx[Obj.GetId()]]['Extinction'][2]['Tbar'] = tbar
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[Indx[Obj.GetId()]]['Extinction'][2]['Tbar']))

    def OnCos2TM(event):
        Obj = event.GetEventObject()
        try:
            val = float(Obj.GetValue())
            if 0. < val <= 1.:
                UseList[Indx[Obj.GetId()]]['Extinction'][2]['Cos2TM'] = val
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[Indx[Obj.GetId()]]['Extinction'][2]['Cos2TM']))
        
    def OnEval(event):
        Obj = event.GetEventObject()
        item = Indx[Obj.GetId()]
        try:
            val = float(Obj.GetValue())
            if val > 0:
                UseList[item[0]]['Extinction'][2][item[1]][0] = val
        except ValueError:
            pass
        Obj.SetValue("%10.3e"%(UseList[item[0]]['Extinction'][2][item[1]][0]))
        
    def OnEref(event):
        Obj = event.GetEventObject()
        item = Indx[Obj.GetId()]
        UseList[item[0]]['Extinction'][2][item[1]][1] = Obj.GetValue()

    def OnSCExtType(event):
        Obj = event.GetEventObject()
        item = Indx[Obj.GetId()]
        UseList[item[0]]['Extinction'][item[1]] = Obj.GetValue()
        wx.CallAfter(UpdateDData,G2frame,DData,data)
            
    def checkAxis(axis):
        if not np.any(np.array(axis)):
            return False
        return axis
        
    def TopSizer(name,choices,parm,OnType):
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(DData,-1,name),0,wx.ALIGN_CENTER_VERTICAL)
        sizeType = wx.ComboBox(DData,wx.ID_ANY,value=UseList[item][parm][0],choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        sizeType.Bind(wx.EVT_COMBOBOX, OnType)
        Indx[sizeType.GetId()] = item
        topSizer.Add(sizeType)
        topSizer.Add((5,0),0)
        return topSizer
        
    def LGmixSizer(name,OnVal,OnRef):
        lgmixSizer = wx.BoxSizer(wx.HORIZONTAL)
        lgmixRef = wx.CheckBox(DData,-1,label='LGmix')
        lgmixRef.thisown = False
        lgmixRef.SetValue(UseList[item][name][2][2])
        Indx[lgmixRef.GetId()] = [item,name]
        lgmixRef.Bind(wx.EVT_CHECKBOX, OnRef)
        lgmixSizer.Add(lgmixRef,0,wx.ALIGN_CENTER_VERTICAL)
        lgmixVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.4f'%(UseList[item][name][1][2]),style=wx.TE_PROCESS_ENTER)
        Indx[lgmixVal.GetId()] = [item,name]
        lgmixVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
        lgmixVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
        lgmixSizer.Add(lgmixVal,0,wx.ALIGN_CENTER_VERTICAL)
        return lgmixSizer
                    
    def ResetSizer(name,OnReset):
        resetSizer = wx.BoxSizer(wx.HORIZONTAL)
        resetSizer.Add((5,0),0)
        reset = wx.CheckBox(DData,-1,label='Reset?')
        reset.thisown = False
        reset.SetValue(False)
        Indx[reset.GetId()] = [item,name]
        reset.Bind(wx.EVT_CHECKBOX,OnReset)
        resetSizer.Add(reset,0,wx.ALIGN_CENTER_VERTICAL)
        return resetSizer
        
    def IsoSizer(name,parm,fmt,OnVal,OnRef):
        isoSizer = wx.BoxSizer(wx.HORIZONTAL)
        sizeRef = wx.CheckBox(DData,-1,label=name)
        sizeRef.thisown = False
        sizeRef.SetValue(UseList[item][parm][2][0])
        Indx[sizeRef.GetId()] = [item,0]
        sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
        isoSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
        sizeVal = wx.TextCtrl(DData,wx.ID_ANY,
            fmt%(UseList[item][parm][1][0]),style=wx.TE_PROCESS_ENTER)
        Indx[sizeVal.GetId()] = [item,0]
        sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
        sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
        isoSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
        return isoSizer
        
    def UniSizer(parm,OnAxis):
        uniSizer = wx.BoxSizer(wx.HORIZONTAL)
        uniSizer.Add(wx.StaticText(DData,-1,' Unique axis, H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
        h,k,l = UseList[item][parm][3]
        Axis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        Indx[Axis.GetId()] = item
        Axis.Bind(wx.EVT_TEXT_ENTER,OnAxis)
        Axis.Bind(wx.EVT_KILL_FOCUS,OnAxis)
        uniSizer.Add(Axis,0,wx.ALIGN_CENTER_VERTICAL)
        return uniSizer
        
    def UniDataSizer(parmName,parm,fmt,OnVal,OnRef):
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        parms = zip([' Equatorial '+parmName,' Axial '+parmName],
            UseList[item][parm][1],UseList[item][parm][2],range(2))
        for Pa,val,ref,id in parms:
            sizeRef = wx.CheckBox(DData,-1,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [item,id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
            dataSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
            sizeVal = wx.TextCtrl(DData,wx.ID_ANY,fmt%(val),style=wx.TE_PROCESS_ENTER)
            Indx[sizeVal.GetId()] = [item,id]
            sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
            sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
            dataSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
            dataSizer.Add((5,0),0)
        return dataSizer
        
    def EllSizeDataSizer():
        parms = zip(['S11','S22','S33','S12','S13','S23'],UseList[item]['Size'][4],
            UseList[item]['Size'][5],range(6))
        dataSizer = wx.FlexGridSizer(1,6,5,5)
        for Pa,val,ref,id in parms:
            sizeRef = wx.CheckBox(DData,-1,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [item,id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
            dataSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
            sizeVal = wx.TextCtrl(DData,wx.ID_ANY,'%.3f'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[sizeVal.GetId()] = [item,id]
            sizeVal.Bind(wx.EVT_TEXT_ENTER,OnSizeVal)
            sizeVal.Bind(wx.EVT_KILL_FOCUS,OnSizeVal)
            dataSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
        return dataSizer
        
    def GenStrainDataSizer():
        Snames = G2spc.MustrainNames(SGData)
        numb = len(Snames)
        if len(UseList[item]['Mustrain'][4]) < numb:
            UseList[item]['Mustrain'][4] = numb*[0.0,]
            UseList[item]['Mustrain'][5] = numb*[False,]
        parms = zip(Snames,UseList[item]['Mustrain'][4],UseList[item]['Mustrain'][5],range(numb))
        dataSizer = wx.FlexGridSizer(1,6,5,5)
        for Pa,val,ref,id in parms:
            strainRef = wx.CheckBox(DData,-1,label=Pa)
            strainRef.thisown = False
            strainRef.SetValue(ref)
            Indx[strainRef.GetId()] = [item,id]
            strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
            dataSizer.Add(strainRef,0,wx.ALIGN_CENTER_VERTICAL)
            strainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.5f'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[strainVal.GetId()] = [item,id]
            strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
            strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
            dataSizer.Add(strainVal,0,wx.ALIGN_CENTER_VERTICAL)
        return dataSizer

    def HstrainSizer():
        hstrainSizer = wx.FlexGridSizer(1,6,5,5)
        Hsnames = G2spc.HStrainNames(SGData)
        parms = zip(Hsnames,UseList[item]['HStrain'][0],UseList[item]['HStrain'][1],range(len(Hsnames)))
        for Pa,val,ref,id in parms:
            hstrainRef = wx.CheckBox(DData,-1,label=Pa)
            hstrainRef.thisown = False
            hstrainRef.SetValue(ref)
            Indx[hstrainRef.GetId()] = [item,id]
            hstrainRef.Bind(wx.EVT_CHECKBOX, OnHstrainRef)
            hstrainSizer.Add(hstrainRef,0,wx.ALIGN_CENTER_VERTICAL)
            hstrainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.5f'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[hstrainVal.GetId()] = [item,id]
            hstrainVal.Bind(wx.EVT_TEXT_ENTER,OnHstrainVal)
            hstrainVal.Bind(wx.EVT_KILL_FOCUS,OnHstrainVal)
            hstrainSizer.Add(hstrainVal,0,wx.ALIGN_CENTER_VERTICAL)
        return hstrainSizer
        
    def PoTopSizer(POData):
        poSizer = wx.FlexGridSizer(1,6,5,5)
        choice = ['March-Dollase','Spherical harmonics']
        POtype = choice[['MD','SH'].index(POData[0])]
        poSizer.Add(wx.StaticText(DData,-1,' Preferred orientation model '),0,wx.ALIGN_CENTER_VERTICAL)
        POType = wx.ComboBox(DData,wx.ID_ANY,value=POtype,choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[POType.GetId()] = item
        POType.Bind(wx.EVT_COMBOBOX, OnPOType)
        poSizer.Add(POType)
        if POData[0] == 'SH':
            poSizer.Add(wx.StaticText(DData,-1,' Harmonic order: '),0,wx.ALIGN_CENTER_VERTICAL)
            poOrder = wx.ComboBox(DData,wx.ID_ANY,value=str(POData[4]),choices=[str(2*i) for i in range(18)],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[poOrder.GetId()] = item
            poOrder.Bind(wx.EVT_COMBOBOX,OnPOOrder)
            poSizer.Add(poOrder,0,wx.ALIGN_CENTER_VERTICAL)
            poRef = wx.CheckBox(DData,-1,label=' Refine? ')
            poRef.SetValue(POData[2])
            Indx[poRef.GetId()] = item
            poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
            poSizer.Add(poRef,0,wx.ALIGN_CENTER_VERTICAL)
        return poSizer
       
    def MDDataSizer(POData):
        poSizer = wx.BoxSizer(wx.HORIZONTAL)
        poRef = wx.CheckBox(DData,-1,label=' March-Dollase ratio: ')
        poRef.SetValue(POData[2])
        Indx[poRef.GetId()] = item
        poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
        poSizer.Add(poRef,0,wx.ALIGN_CENTER_VERTICAL)
        poVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.3f'%(POData[1]),style=wx.TE_PROCESS_ENTER)
        Indx[poVal.GetId()] = item
        poVal.Bind(wx.EVT_TEXT_ENTER,OnPOVal)
        poVal.Bind(wx.EVT_KILL_FOCUS,OnPOVal)
        poSizer.Add(poVal,0,wx.ALIGN_CENTER_VERTICAL)
        poSizer.Add(wx.StaticText(DData,-1,' Unique axis, H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
        h,k,l =POData[3]
        poAxis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        Indx[poAxis.GetId()] = item
        poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOAxis)
        poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOAxis)
        poSizer.Add(poAxis,0,wx.ALIGN_CENTER_VERTICAL)
        return poSizer
        
    def SHDataSizer(POData):
        textJ = G2lat.textureIndex(POData[5])
        mainSizer.Add(wx.StaticText(DData,-1,' Spherical harmonic coefficients: '+'Texture index: %.3f'%(textJ)),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((0,5),0)
        ODFSizer = wx.FlexGridSizer(2,8,2,2)
        ODFIndx = {}
        ODFkeys = POData[5].keys()
        ODFkeys.sort()
        for odf in ODFkeys:
            ODFSizer.Add(wx.StaticText(DData,-1,odf),0,wx.ALIGN_CENTER_VERTICAL)
            ODFval = wx.TextCtrl(DData,wx.ID_ANY,'%8.3f'%(POData[5][odf]),style=wx.TE_PROCESS_ENTER)
            ODFIndx[ODFval.GetId()] = odf
#             ODFval.Bind(wx.EVT_TEXT_ENTER,OnODFValue)
#             ODFval.Bind(wx.EVT_KILL_FOCUS,OnODFValue)
            ODFSizer.Add(ODFval,0,wx.ALIGN_CENTER_VERTICAL)
        return ODFSizer
        
    def ExtSizer():            
        extSizer = wx.BoxSizer(wx.HORIZONTAL)
        extRef = wx.CheckBox(DData,-1,label=' Extinction: ')
        extRef.SetValue(UseList[item]['Extinction'][1])
        Indx[extRef.GetId()] = item
        extRef.Bind(wx.EVT_CHECKBOX, OnExtRef)
        extSizer.Add(extRef,0,wx.ALIGN_CENTER_VERTICAL)
        extVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.2f'%(UseList[item]['Extinction'][0]),style=wx.TE_PROCESS_ENTER)
        Indx[extVal.GetId()] = item
        extVal.Bind(wx.EVT_TEXT_ENTER,OnExtVal)
        extVal.Bind(wx.EVT_KILL_FOCUS,OnExtVal)
        extSizer.Add(extVal,0,wx.ALIGN_CENTER_VERTICAL)
        return extSizer
    
    def SCExtSizer():
        extSizer = wx.BoxSizer(wx.VERTICAL)
        typeSizer = wx.BoxSizer(wx.HORIZONTAL)            
        typeSizer.Add(wx.StaticText(DData,-1,' Extinction type: '),0,wx.ALIGN_CENTER_VERTICAL)
        Choices = ['None','Primary','Secondary Type I','Secondary Type II','Secondary Type I & II']
        typeTxt = wx.ComboBox(DData,-1,choices=Choices,value=UseList[item]['Extinction'][1],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[typeTxt.GetId()] = [item,1]
        typeTxt.Bind(wx.EVT_COMBOBOX,OnSCExtType)
        typeSizer.Add(typeTxt)
        typeSizer.Add(wx.StaticText(DData,-1,' Approx: '),0,wx.ALIGN_CENTER_VERTICAL)
        Choices=['Lorentzian','Gaussian']
        approxTxT = wx.ComboBox(DData,-1,choices=Choices,value=UseList[item]['Extinction'][0],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[approxTxT.GetId()] = [item,0]
        approxTxT.Bind(wx.EVT_COMBOBOX,OnSCExtType)
        typeSizer.Add(approxTxT)
        extSizer.Add(typeSizer,0,wx.ALIGN_CENTER_VERTICAL)
        if UseList[item]['Extinction'][1] != 'None':
            extSizer.Add((0,5),)
            valSizer =wx.BoxSizer(wx.HORIZONTAL)
            valSizer.Add(wx.StaticText(DData,-1,' Tbar(mm):'),0,wx.ALIGN_CENTER_VERTICAL)
            tbarVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.3f'%(UseList[item]['Extinction'][2]['Tbar']),style=wx.TE_PROCESS_ENTER)
            Indx[tbarVal.GetId()] = item
            tbarVal.Bind(wx.EVT_TEXT_ENTER,OnTbarVal)
            tbarVal.Bind(wx.EVT_KILL_FOCUS,OnTbarVal)
            valSizer.Add(tbarVal,0,wx.ALIGN_CENTER_VERTICAL)
            valSizer.Add(wx.StaticText(DData,-1,' cos(2ThM):'),0,wx.ALIGN_CENTER_VERTICAL)
            cos2tm = wx.TextCtrl(DData,wx.ID_ANY,
                '%.3f'%(UseList[item]['Extinction'][2]['Cos2TM']),style=wx.TE_PROCESS_ENTER)
            Indx[cos2tm.GetId()] = item
            cos2tm.Bind(wx.EVT_TEXT_ENTER,OnCos2TM)
            cos2tm.Bind(wx.EVT_KILL_FOCUS,OnCos2TM)
            valSizer.Add(cos2tm,0,wx.ALIGN_CENTER_VERTICAL)
            extSizer.Add(valSizer,0,wx.ALIGN_CENTER_VERTICAL)
            val2Sizer =wx.BoxSizer(wx.HORIZONTAL)
            if 'Primary' in UseList[item]['Extinction'][1]:
                Ekey = ['Ep',]
            elif 'Secondary Type II' == UseList[item]['Extinction'][1]:
                Ekey = ['Es',]
            elif 'Secondary Type I' == UseList[item]['Extinction'][1]:
                Ekey = ['Eg',]
            else:
                Ekey = ['Eg','Es']
            for ekey in Ekey:
                Eref = wx.CheckBox(DData,-1,label=ekey+' : ')
                Eref.SetValue(UseList[item]['Extinction'][2][ekey][1])
                Indx[Eref.GetId()] = [item,ekey]
                Eref.Bind(wx.EVT_CHECKBOX, OnEref)
                val2Sizer.Add(Eref,0,wx.ALIGN_CENTER_VERTICAL)
                Eval = wx.TextCtrl(DData,wx.ID_ANY,
                    '%10.3e'%(UseList[item]['Extinction'][2][ekey][0]),style=wx.TE_PROCESS_ENTER)
                Indx[Eval.GetId()] = [item,ekey]
                Eval.Bind(wx.EVT_TEXT_ENTER,OnEval)
                Eval.Bind(wx.EVT_KILL_FOCUS,OnEval)
                val2Sizer.Add(Eval,0,wx.ALIGN_CENTER_VERTICAL)

            extSizer.Add(val2Sizer,0,wx.ALIGN_CENTER_VERTICAL)
        return extSizer
        
    def BabSizer():
        babSizer = wx.BoxSizer(wx.HORIZONTAL)
        for bab in ['A','U']:
            babRef = wx.CheckBox(DData,-1,label=' Babinet '+bab+': ')
            babRef.SetValue(UseList[item]['Babinet']['Bab'+bab][1])
            Indx[babRef.GetId()] = [item,bab]
            babRef.Bind(wx.EVT_CHECKBOX, OnBabRef)
            babSizer.Add(babRef,0,wx.ALIGN_CENTER_VERTICAL)
            babVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.3f'%(UseList[item]['Babinet']['Bab'+bab][0]),style=wx.TE_PROCESS_ENTER)
            Indx[babVal.GetId()] = [item,bab]
            babVal.Bind(wx.EVT_TEXT_ENTER,OnBabVal)
            babVal.Bind(wx.EVT_KILL_FOCUS,OnBabVal)
            babSizer.Add(babVal,0,wx.ALIGN_CENTER_VERTICAL)
            babSizer.Add((5,5),0)
        return babSizer
        
    if DData.GetSizer():
        DData.GetSizer().Clear(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(DData,-1,'Histogram data for '+PhaseName+':'),0,wx.ALIGN_CENTER_VERTICAL)
    if PWDR:
        mainSizer.Add(PlotSizer())            
        
    for item in keyList:
        histData = UseList[item]
        if 'Use' not in UseList[item]:      #patch
            UseList[item]['Use'] = True
        if 'Babinet' not in UseList[item]:
            UseList[item]['Babinet'] = {'BabA':[0.0,False],'BabU':[0.0,False]}
        showSizer = wx.BoxSizer(wx.HORIZONTAL)
        showData = wx.CheckBox(DData,-1,label=' Show '+item)
        showData.SetValue(UseList[item]['Show'])
        Indx[showData.GetId()] = item
        showData.Bind(wx.EVT_CHECKBOX, OnShowData)
        showSizer.Add(showData,0,wx.ALIGN_CENTER_VERTICAL)
        useData = wx.CheckBox(DData,-1,label='Use?')
        Indx[useData.GetId()] = item
        showSizer.Add(useData,0,wx.ALIGN_CENTER_VERTICAL)
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(UseList[item]['Use'])
        copyData = wx.Button(DData,-1,label=' Copy?')
        Indx[copyData.GetId()] = item
        copyData.Bind(wx.EVT_BUTTON,OnCopyData)
        showSizer.Add(copyData,wx.ALIGN_CENTER_VERTICAL)
        copyFlags = wx.Button(DData,-1,label=' Copy flags?')
        Indx[copyFlags.GetId()] = item
        copyFlags.Bind(wx.EVT_BUTTON,OnCopyFlags)
        showSizer.Add(copyFlags,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(showSizer,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((0,5),0)
        
        if UseList[item]['Show']:
            mainSizer.Add(ScaleSizer())
            mainSizer.Add((0,5),0)
            
        if item[:4] == 'PWDR' and UseList[item]['Show']:
            if UseList[item]['Size'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                isoSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetSize),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(isoSizer)
                mainSizer.Add(IsoSizer(u' Cryst. size(\xb5m): ','Size','%.3f',
                    OnSizeVal,OnSizeRef),0,wx.ALIGN_CENTER_VERTICAL)
            elif UseList[item]['Size'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                uniSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetSize),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(UniSizer('Size',OnSizeAxis),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(uniSizer)
                mainSizer.Add(UniDataSizer(u'size(\xb5m): ','Size','%.3f',OnSizeVal,OnSizeRef))
            elif UseList[item]['Size'][0] == 'ellipsoidal':
                ellSizer = wx.BoxSizer(wx.HORIZONTAL)
                ellSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                ellSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                ellSizer.Add(ResetSizer('ellipsoidal',OnResetSize),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(ellSizer)
                mainSizer.Add(EllSizeDataSizer())
            mainSizer.Add((0,5),0)                    
            
            if UseList[item]['Mustrain'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                isoSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetStrain),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(isoSizer)
                mainSizer.Add(IsoSizer(' microstrain: ','Mustrain','%.1f',
                    OnStrainVal,OnStrainRef),0,wx.ALIGN_CENTER_VERTICAL)                   
                mainSizer.Add((0,5),0)
            elif UseList[item]['Mustrain'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                uniSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetStrain),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(uniSizer)
                mainSizer.Add(UniSizer('Mustrain',OnStrainAxis),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(UniDataSizer('mustrain: ','Mustrain','%.1f',OnStrainVal,OnStrainRef))
            elif UseList[item]['Mustrain'][0] == 'generalized':
                genSizer = wx.BoxSizer(wx.HORIZONTAL)
                genSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                genSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                genSizer.Add(ResetSizer('generalized',OnResetStrain),0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(genSizer)
                mainSizer.Add(GenStrainDataSizer())                        
            mainSizer.Add((0,5),0)
            
            mainSizer.Add(wx.StaticText(DData,-1,' Hydrostatic/elastic strain:'))
            mainSizer.Add(HstrainSizer())
                
            #texture  'Pref. Ori.':['MD',1.0,False,[0,0,1],0,[]] last two for 'SH' are SHorder & coeff
            poSizer = wx.BoxSizer(wx.VERTICAL)
            POData = UseList[item]['Pref.Ori.']
            poSizer.Add(PoTopSizer(POData))
            if POData[0] == 'MD':
                poSizer.Add(MDDataSizer(POData))
            else:           #'SH'
                if POData[4]:       #SH order > 0
                    poSizer.Add(SHDataSizer(POData))
                    
            mainSizer.Add(poSizer)
            mainSizer.Add((0,5),0)                
            mainSizer.Add(ExtSizer())
            mainSizer.Add((0,5),0)
            mainSizer.Add(BabSizer())
            mainSizer.Add((0,5),0)
        elif item[:4] == 'HKLF' and UseList[item]['Show']:
            mainSizer.Add((0,5),0)                
            mainSizer.Add(SCExtSizer())
            mainSizer.Add((0,5),0)
            mainSizer.Add(BabSizer())
            mainSizer.Add((0,5),0)
            pass
        #G2gd.HorizontalLine(mainSizer,DData)
    mainSizer.Add((5,5),0)

    DData.SetSizer(mainSizer,True)
    if G2frame.dataFrame.PhaseUserSize is None:
        Size = mainSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],290) + 35
        DData.SetSize(Size)
        DData.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],450)
        G2frame.dataFrame.setSizePosLeft(Size)
    else:
        Size = G2frame.dataFrame.PhaseUserSize
        DData.SetSize(G2frame.dataFrame.GetClientSize())
        DData.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.Update()
