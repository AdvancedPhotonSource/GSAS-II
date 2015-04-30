# -*- coding: utf-8 -*-
#GSASII - phase data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIddataGUI: Phase Diffraction Data GUI*
--------------------------------------------

Module to create the GUI for display of diffraction data * phase
information that is shown in the data display window
(when a phase is selected.)

'''
import wx
import math
import copy
import time
import sys
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIIO as G2IO
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import GSASIIphsGUI as G2phsGUI
import GSASIIctrls as G2G
import numpy as np

WACV = wx.ALIGN_CENTER_VERTICAL
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)
BLACK = wx.Colour(0,0,0)
mapDefault = {'MapType':'','RefList':'','Resolution':0.5,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}

################################################################################
##### DData routines
################################################################################        
def UpdateDData(G2frame,DData,data,hist=''):
    '''Display the Diffraction Data associated with a phase
    (items where there is a value for each histogram and phase)

    :param wx.frame G2frame: the main GSAS-II frame object

    :param wx.ScrolledWindow DData: notebook page to be used for the display

    :param dict data: all the information on the phase in a dictionary

    '''
    G2frame.dataFrame.SetStatusText('')
    keyList = G2frame.GetHistogramNames(['PWDR','HKLF'])
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
    if len(UseList) == 0: # no associated, don't display anything
        G2frame.hist = '' 
    elif hist: # something was input as a selection
        G2frame.hist = hist
    elif not G2frame.hist or G2frame.hist not in UseList: # no or bad selection but have data, take the first
        for key in keyList:
            if key in UseList:
                    G2frame.hist = key
                    break
    PWDR = any(['PWDR' in item for item in keyList])
    Indx = {}
    
    def PlotSizer():

        def OnPlotSel(event):
            Obj = event.GetEventObject()
            generalData['Data plot type'] = Obj.GetStringSelection()
            G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
            wx.CallLater(100,UpdateDData,G2frame,DData,data,G2frame.hist)
            
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
            G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        
        plotSizer = wx.BoxSizer(wx.VERTICAL)
        choice = ['None','Mustrain','Size','Preferred orientation']
        plotSel = wx.RadioBox(DData,wx.ID_ANY,'Select plot type:',choices=choice,
            majorDimension=1,style=wx.RA_SPECIFY_COLS)
        plotSel.SetStringSelection(generalData['Data plot type'])
        plotSel.Bind(wx.EVT_RADIOBOX,OnPlotSel)    
        plotSizer.Add(plotSel)
        if generalData['Data plot type'] == 'Preferred orientation':
            POhklSizer = wx.BoxSizer(wx.HORIZONTAL)
            POhklSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Plot preferred orientation for H K L: '),0,WACV)
            h,k,l = generalData['POhkl']
            poAxis = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
            poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOhkl)
            poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOhkl)
            POhklSizer.Add(poAxis,0,WACV)
            plotSizer.Add(POhklSizer)            
        return plotSizer
       
    def ScaleSizer():
        
        def OnScaleRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Scale'][1] = Obj.GetValue()
            
        def OnScaleVal(event):
            Obj = event.GetEventObject()
            try:
                scale = float(Obj.GetValue())
                if scale > 0:
                    UseList[G2frame.hist]['Scale'][0] = scale
            except ValueError:
                pass
            Obj.SetValue("%.4f"%(UseList[G2frame.hist]['Scale'][0]))          #reset in case of error
                        
        scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
        if 'PWDR' in G2frame.hist:
            scaleRef = wx.CheckBox(DData,wx.ID_ANY,label=' Phase fraction: ')
        elif 'HKLF' in G2frame.hist:
            scaleRef = wx.CheckBox(DData,wx.ID_ANY,label=' Scale factor: ')                
        scaleRef.SetValue(UseList[G2frame.hist]['Scale'][1])
        scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        scaleSizer.Add(scaleRef,0,WACV)
        scaleVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.4f'%(UseList[G2frame.hist]['Scale'][0]),style=wx.TE_PROCESS_ENTER)
        scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
        scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
        scaleSizer.Add(scaleVal,0,WACV)
        return scaleSizer
        
    def OnUseData(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Use'] = Obj.GetValue()
        
    def OnLGmixRef(event):
        Obj = event.GetEventObject()
        hist,name = Indx[Obj.GetId()]
        UseList[G2frame.hist][name][2][2] = Obj.GetValue()
        
    def OnLGmixVal(event):
        Obj = event.GetEventObject()
        hist,name = Indx[Obj.GetId()]
        try:
            value = float(Obj.GetValue())
            if 0 <= value <= 1:
                UseList[G2frame.hist][name][1][2] = value
            else:
                raise ValueError
        except ValueError:
            pass
        Obj.SetValue("%.4f"%(UseList[G2frame.hist][name][1][2]))          #reset in case of error

    def OnSizeType(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Size'][0] = Obj.GetValue()
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)
        
    def OnSizeRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[G2frame.hist]['Size'][0] == 'ellipsoidal':
            UseList[G2frame.hist]['Size'][5][pid] = Obj.GetValue()                
        else:
            UseList[G2frame.hist]['Size'][2][pid] = Obj.GetValue()
        
    def OnSizeVal(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[G2frame.hist]['Size'][0] == 'ellipsoidal':
            try:
                size = float(Obj.GetValue())
                if pid < 3 and size <= 0.001:            #10A lower limit!
                    raise ValueError
                UseList[G2frame.hist]['Size'][4][pid] = size                    
            except ValueError:
                pass
            Obj.SetValue("%.5f"%(UseList[G2frame.hist]['Size'][4][pid]))          #reset in case of error
        else:
            try:
                size = float(Obj.GetValue())
                if size <= 0.001:            #10A lower limit!
                    raise ValueError
                UseList[G2frame.hist]['Size'][1][pid] = size
            except ValueError:
                pass
            Obj.SetValue("%.5f"%(UseList[G2frame.hist]['Size'][1][pid]))          #reset in case of error
        G2plt.PlotSizeStrainPO(G2frame,data,hist)
        
    def OnSizeAxis(event):            
        Obj = event.GetEventObject()
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[G2frame.hist]['Size'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[G2frame.hist]['Size'][3]
        UseList[G2frame.hist]['Size'][3] = hkl
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
        G2plt.PlotSizeStrainPO(G2frame,data,item)
        wx.CallLater(100,RepaintHistogramInfo)
            
    def OnStrainType(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Mustrain'][0] = Obj.GetValue()
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)
        
    def OnStrainRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
            UseList[G2frame.hist]['Mustrain'][5][pid] = Obj.GetValue()
        else:
            UseList[G2frame.hist]['Mustrain'][2][pid] = Obj.GetValue()
        
    def OnStrainVal(event):
        Snames = G2spc.MustrainNames(SGData)
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        try:
            strain = float(Obj.GetValue())
            if UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
                if '4' in Snames[pid] and strain < 0:
                    raise ValueError
                UseList[G2frame.hist]['Mustrain'][4][pid] = strain
            else:
                if strain <= 0:
                    raise ValueError
                UseList[G2frame.hist]['Mustrain'][1][pid] = strain
        except ValueError:
            pass
        if UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
            Obj.SetValue("%.3f"%(UseList[G2frame.hist]['Mustrain'][4][pid]))          #reset in case of error
        else:
            Obj.SetValue("%.1f"%(UseList[G2frame.hist]['Mustrain'][1][pid]))          #reset in case of error
        G2plt.PlotSizeStrainPO(G2frame,data,hist)
        
    def OnStrainAxis(event):
        Obj = event.GetEventObject()
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[G2frame.hist]['Mustrain'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[G2frame.hist]['Mustrain'][3]
        UseList[G2frame.hist]['Mustrain'][3] = hkl
        h,k,l = hkl
        Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        
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
            muiso = 1000.
            cell = generalData['Cell'][1:7]
            vals = G2spc.Muiso2Shkl(muiso,SGData,cell)
            nTerm = len(UseList[item]['Mustrain'][4])
            for i in range(nTerm):
                UseList[item]['Mustrain'][4][i] = vals[i]
        G2plt.PlotSizeStrainPO(G2frame,data,item)
        wx.CallLater(100,RepaintHistogramInfo)
            
    def OnHstrainRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        UseList[G2frame.hist]['HStrain'][1][pid] = Obj.GetValue()
        
    def OnHstrainVal(event):
        Snames = G2spc.HStrainNames(SGData)
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        try:
            strain = float(Obj.GetValue())
            UseList[G2frame.hist]['HStrain'][0][pid] = strain
        except ValueError:
            pass
        Obj.SetValue("%.3g"%(UseList[G2frame.hist]['HStrain'][0][pid]))          #reset in case of error

    def OnPOVal(event):
        Obj = event.GetEventObject()
        try:
            mdVal = float(Obj.GetValue())
            if mdVal > 0:
                UseList[G2frame.hist]['Pref.Ori.'][1] = mdVal
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[G2frame.hist]['Pref.Ori.'][1]))          #reset in case of error
        
    def OnPOAxis(event):
        Obj = event.GetEventObject()
        Saxis = Obj.GetValue().split()
        try:
            hkl = [int(Saxis[i]) for i in range(3)]
        except (ValueError,IndexError):
            hkl = UseList[G2frame.hist]['Pref.Ori.'][3]
        if not np.any(np.array(hkl)):
            hkl = UseList[G2frame.hist]['Pref.Ori.'][3]
        UseList[G2frame.hist]['Pref.Ori.'][3] = hkl
        h,k,l = hkl
        Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
        
    def OnPOOrder(event):
        Obj = event.GetEventObject()
        Order = int(Obj.GetValue())
        UseList[G2frame.hist]['Pref.Ori.'][4] = Order
        UseList[G2frame.hist]['Pref.Ori.'][5] = SetPOCoef(Order,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)

    def OnPOType(event):
        Obj = event.GetEventObject()
        if 'March' in Obj.GetValue():
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'MD'
        else:
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'SH'
        wx.CallLater(100,RepaintHistogramInfo)

    def OnPORef(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Pref.Ori.'][2] = Obj.GetValue()
            
    def SetPOCoef(Order,hist):
        cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',Order,False)     #cylindrical & no M
        newPOCoef = dict(zip(cofNames,np.zeros(len(cofNames))))
        POCoeff = UseList[G2frame.hist]['Pref.Ori.'][5]
        for cofName in POCoeff:
            if cofName in  cofNames:
                newPOCoef[cofName] = POCoeff[cofName]
        return newPOCoef
        
    def OnExtRef(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Extinction'][1] = Obj.GetValue()
        
    def OnExtVal(event):
        Obj = event.GetEventObject()
        try:
            ext = float(Obj.GetValue())
            if ext >= 0:
                UseList[G2frame.hist]['Extinction'][0] = ext
        except ValueError:
            pass
        Obj.SetValue("%.2f"%(UseList[G2frame.hist]['Extinction'][0]))

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
                UseList[G2frame.hist]['Extinction'][2]['Tbar'] = tbar
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[G2frame.hist]['Extinction'][2]['Tbar']))

    def OnCos2TM(event):
        Obj = event.GetEventObject()
        try:
            val = float(Obj.GetValue())
            if 0. < val <= 1.:
                UseList[G2frame.hist]['Extinction'][2]['Cos2TM'] = val
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(UseList[G2frame.hist]['Extinction'][2]['Cos2TM']))
        
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
        wx.CallLater(100,RepaintHistogramInfo)
            
    def checkAxis(axis):
        if not np.any(np.array(axis)):
            return False
        return axis
        
    def TopSizer(name,choices,parm,OnType):
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(DData,wx.ID_ANY,name),0,WACV)
        sizeType = wx.ComboBox(DData,wx.ID_ANY,value=UseList[G2frame.hist][parm][0],choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        sizeType.Bind(wx.EVT_COMBOBOX, OnType)
        topSizer.Add(sizeType,0,WACV|wx.BOTTOM,5)
        return topSizer
        
    def LGmixSizer(name,OnVal,OnRef):
        lgmixSizer = wx.BoxSizer(wx.HORIZONTAL)
        lgmixRef = wx.CheckBox(DData,wx.ID_ANY,label='LGmix')
        lgmixRef.thisown = False
        lgmixRef.SetValue(UseList[G2frame.hist][name][2][2])
        Indx[lgmixRef.GetId()] = [G2frame.hist,name]
        lgmixRef.Bind(wx.EVT_CHECKBOX, OnRef)
        lgmixSizer.Add(lgmixRef,0,WACV)
        lgmixVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.4f'%(UseList[G2frame.hist][name][1][2]),style=wx.TE_PROCESS_ENTER)
        Indx[lgmixVal.GetId()] = [G2frame.hist,name]
        lgmixVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
        lgmixVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
        lgmixSizer.Add(lgmixVal,0,WACV)
        return lgmixSizer
                    
    def ResetSizer(name,OnReset):
        resetSizer = wx.BoxSizer(wx.HORIZONTAL)
        reset = wx.CheckBox(DData,wx.ID_ANY,label='Reset?')
        reset.thisown = False
        reset.SetValue(False)
        Indx[reset.GetId()] = [G2frame.hist,name]
        reset.Bind(wx.EVT_CHECKBOX,OnReset)
        resetSizer.Add(reset,0,WACV|wx.TOP,5)
        return resetSizer
        
    def IsoSizer(name,parm,fmt,OnVal,OnRef):
        isoSizer = wx.BoxSizer(wx.HORIZONTAL)
        sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=name)
        sizeRef.thisown = False
        sizeRef.SetValue(UseList[G2frame.hist][parm][2][0])
        Indx[sizeRef.GetId()] = [G2frame.hist,0]
        sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
        isoSizer.Add(sizeRef,0,WACV)
        sizeVal = wx.TextCtrl(DData,wx.ID_ANY,
            fmt%(UseList[G2frame.hist][parm][1][0]),style=wx.TE_PROCESS_ENTER)
        Indx[sizeVal.GetId()] = [G2frame.hist,0]
        sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
        sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
        isoSizer.Add(sizeVal,0,WACV)
        return isoSizer
        
    def UniSizer(parm,OnAxis):
        uniSizer = wx.BoxSizer(wx.HORIZONTAL)
        uniSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Unique axis, H K L: '),0,WACV)
        h,k,l = UseList[G2frame.hist][parm][3]
        Axis = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        Axis.Bind(wx.EVT_TEXT_ENTER,OnAxis)
        Axis.Bind(wx.EVT_KILL_FOCUS,OnAxis)
        uniSizer.Add(Axis,0,WACV)
        return uniSizer
        
    def UniDataSizer(parmName,parm,fmt,OnVal,OnRef):
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        parms = zip([' Equatorial '+parmName,' Axial '+parmName],
            UseList[G2frame.hist][parm][1],UseList[G2frame.hist][parm][2],range(2))
        for Pa,val,ref,id in parms:
            sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [G2frame.hist,id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
            dataSizer.Add(sizeRef,0,WACV)
            sizeVal = wx.TextCtrl(DData,wx.ID_ANY,fmt%(val),style=wx.TE_PROCESS_ENTER)
            Indx[sizeVal.GetId()] = [G2frame.hist,id]
            sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
            sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
            dataSizer.Add(sizeVal,0,WACV|wx.BOTTOM,5)
        return dataSizer
        
    def EllSizeDataSizer():
        parms = zip(['S11','S22','S33','S12','S13','S23'],UseList[G2frame.hist]['Size'][4],
            UseList[G2frame.hist]['Size'][5],range(6))
        dataSizer = wx.FlexGridSizer(0,6,5,5)
        for Pa,val,ref,id in parms:
            sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [G2frame.hist,id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
            dataSizer.Add(sizeRef,0,WACV)
            sizeVal = wx.TextCtrl(DData,wx.ID_ANY,'%.3f'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[sizeVal.GetId()] = [G2frame.hist,id]
            sizeVal.Bind(wx.EVT_TEXT_ENTER,OnSizeVal)
            sizeVal.Bind(wx.EVT_KILL_FOCUS,OnSizeVal)
            dataSizer.Add(sizeVal,0,WACV)
        return dataSizer
        
    def GenStrainDataSizer():
        Snames = G2spc.MustrainNames(SGData)
        numb = len(Snames)
        if len(UseList[G2frame.hist]['Mustrain'][4]) < numb:
            UseList[G2frame.hist]['Mustrain'][4] = numb*[0.0,]
            UseList[G2frame.hist]['Mustrain'][5] = numb*[False,]
        parms = zip(Snames,UseList[G2frame.hist]['Mustrain'][4],UseList[G2frame.hist]['Mustrain'][5],range(numb))
        dataSizer = wx.FlexGridSizer(0,6,5,5)
        for Pa,val,ref,id in parms:
            strainRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            strainRef.thisown = False
            strainRef.SetValue(ref)
            Indx[strainRef.GetId()] = [G2frame.hist,id]
            strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
            dataSizer.Add(strainRef,0,WACV)
            strainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.5f'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[strainVal.GetId()] = [G2frame.hist,id]
            strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
            strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
            dataSizer.Add(strainVal,0,WACV)
        return dataSizer

    def HstrainSizer():
        hstrainSizer = wx.FlexGridSizer(0,6,5,5)
        Hsnames = G2spc.HStrainNames(SGData)
        parms = zip(Hsnames,UseList[G2frame.hist]['HStrain'][0],UseList[G2frame.hist]['HStrain'][1],range(len(Hsnames)))
        for Pa,val,ref,id in parms:
            hstrainRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            hstrainRef.thisown = False
            hstrainRef.SetValue(ref)
            Indx[hstrainRef.GetId()] = [G2frame.hist,id]
            hstrainRef.Bind(wx.EVT_CHECKBOX, OnHstrainRef)
            hstrainSizer.Add(hstrainRef,0,WACV)
            hstrainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.3g'%(val),style=wx.TE_PROCESS_ENTER)
            Indx[hstrainVal.GetId()] = [G2frame.hist,id]
            hstrainVal.Bind(wx.EVT_TEXT_ENTER,OnHstrainVal)
            hstrainVal.Bind(wx.EVT_KILL_FOCUS,OnHstrainVal)
            hstrainSizer.Add(hstrainVal,0,WACV)
        return hstrainSizer
        
    def PoTopSizer(POData):
        poSizer = wx.FlexGridSizer(0,6,5,5)
        choice = ['March-Dollase','Spherical harmonics']
        POtype = choice[['MD','SH'].index(POData[0])]
        poSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Preferred orientation model '),0,WACV)
        POType = wx.ComboBox(DData,wx.ID_ANY,value=POtype,choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        POType.Bind(wx.EVT_COMBOBOX, OnPOType)
        poSizer.Add(POType)
        if POData[0] == 'SH':
            poSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Harmonic order: '),0,WACV)
            poOrder = wx.ComboBox(DData,wx.ID_ANY,value=str(POData[4]),choices=[str(2*i) for i in range(18)],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            poOrder.Bind(wx.EVT_COMBOBOX,OnPOOrder)
            poSizer.Add(poOrder,0,WACV)
            poRef = wx.CheckBox(DData,wx.ID_ANY,label=' Refine? ')
            poRef.SetValue(POData[2])
            poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
            poSizer.Add(poRef,0,WACV)
        return poSizer
       
    def MDDataSizer(POData):
        poSizer = wx.BoxSizer(wx.HORIZONTAL)
        poRef = wx.CheckBox(DData,wx.ID_ANY,label=' March-Dollase ratio: ')
        poRef.SetValue(POData[2])
        poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
        poSizer.Add(poRef,0,WACV)
        poVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.3f'%(POData[1]),style=wx.TE_PROCESS_ENTER)
        poVal.Bind(wx.EVT_TEXT_ENTER,OnPOVal)
        poVal.Bind(wx.EVT_KILL_FOCUS,OnPOVal)
        poSizer.Add(poVal,0,WACV)
        poSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Unique axis, H K L: '),0,WACV)
        h,k,l =POData[3]
        poAxis = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOAxis)
        poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOAxis)
        poSizer.Add(poAxis,0,WACV)
        return poSizer
        
    def SHDataSizer(POData):
        
        def OnODFValue(event):
            Obj = event.GetEventObject()
            odf = ODFIndx[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
                POData[5][odf] = value
            except ValueError:
                pass
            Obj.SetValue('%8.3f'%(POData[5][odf]))
            G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
    
        ODFSizer = wx.FlexGridSizer(0,8,2,2)
        ODFIndx = {}
        ODFkeys = POData[5].keys()
        ODFkeys.sort()
        for odf in ODFkeys:
            ODFSizer.Add(wx.StaticText(DData,wx.ID_ANY,odf),0,WACV)
            ODFval = wx.TextCtrl(DData,wx.ID_ANY,'%8.3f'%(POData[5][odf]),style=wx.TE_PROCESS_ENTER)
            ODFIndx[ODFval.GetId()] = odf
            ODFval.Bind(wx.EVT_TEXT_ENTER,OnODFValue)
            ODFval.Bind(wx.EVT_KILL_FOCUS,OnODFValue)
            ODFSizer.Add(ODFval,0,WACV)
        return ODFSizer
        
    def SHPenalty(POData):
        
        def OnHKLList(event):
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select penalty hkls',
                'Penalty hkls',hkls,filterBox=False)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    POData[6] = [hkls[i] for i in dlg.GetSelections()]
                    if not POData[6]:
                        POData[6] = ['',]
                else:
                    return
            finally:
                dlg.Destroy()
            wx.CallLater(100,RepaintHistogramInfo)
            
        def OnshToler(event):
            try:
                value = float(shToler.GetValue())
                POData[7] = value
            except ValueError:
                pass
            shToler.SetValue('%.2f'%(POData[7]))
        
        A = G2lat.cell2A(generalData['Cell'][1:7])
        hkls = G2lat.GenPfHKLs(10,SGData,A)    
        shPenalty = wx.BoxSizer(wx.HORIZONTAL)
        shPenalty.Add(wx.StaticText(DData,wx.ID_ANY,' Negative MRD penalty list: '),0,WACV)
        shPenalty.Add(wx.ComboBox(DData,value=POData[6][0],choices=POData[6],
            style=wx.CB_DROPDOWN),0,WACV)
        hklList = wx.Button(DData,label='Select penalty hkls')
        hklList.Bind(wx.EVT_BUTTON,OnHKLList)
        shPenalty.Add(hklList,0,WACV)
        shPenalty.Add(wx.StaticText(DData,wx.ID_ANY,' Zero MRD tolerance: '),0,WACV)
        shToler = wx.TextCtrl(DData,wx.ID_ANY,'%.2f'%(POData[7]),style=wx.TE_PROCESS_ENTER)
        shToler.Bind(wx.EVT_TEXT_ENTER,OnshToler)
        shToler.Bind(wx.EVT_KILL_FOCUS,OnshToler)
        shPenalty.Add(shToler,0,WACV)
        return shPenalty    
        
    def ExtSizer():            
        extSizer = wx.BoxSizer(wx.HORIZONTAL)
        extRef = wx.CheckBox(DData,wx.ID_ANY,label=' Extinction: ')
        extRef.SetValue(UseList[G2frame.hist]['Extinction'][1])
        extRef.Bind(wx.EVT_CHECKBOX, OnExtRef)
        extSizer.Add(extRef,0,WACV)
        extVal = wx.TextCtrl(DData,wx.ID_ANY,
            '%.2f'%(UseList[G2frame.hist]['Extinction'][0]),style=wx.TE_PROCESS_ENTER)
        extVal.Bind(wx.EVT_TEXT_ENTER,OnExtVal)
        extVal.Bind(wx.EVT_KILL_FOCUS,OnExtVal)
        extSizer.Add(extVal,0,WACV)
        return extSizer
    
    def SCExtSizer():
        extSizer = wx.BoxSizer(wx.VERTICAL)
        typeSizer = wx.BoxSizer(wx.HORIZONTAL)            
        typeSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Extinction type: '),0,WACV)
        Choices = ['None','Primary','Secondary Type I','Secondary Type II',]    # remove 'Secondary Type I & II'
        typeTxt = wx.ComboBox(DData,wx.ID_ANY,choices=Choices,value=UseList[G2frame.hist]['Extinction'][1],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[typeTxt.GetId()] = [G2frame.hist,1]
        typeTxt.Bind(wx.EVT_COMBOBOX,OnSCExtType)
        typeSizer.Add(typeTxt)
        typeSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Approx: '),0,WACV)
        Choices=['Lorentzian','Gaussian']
        approxTxT = wx.ComboBox(DData,wx.ID_ANY,choices=Choices,value=UseList[G2frame.hist]['Extinction'][0],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[approxTxT.GetId()] = [G2frame.hist,0]
        approxTxT.Bind(wx.EVT_COMBOBOX,OnSCExtType)
        typeSizer.Add(approxTxT)
        if UseList[G2frame.hist]['Extinction'][1] == 'None':
            extSizer.Add(typeSizer,0,WACV)
        else:
            extSizer.Add(typeSizer,0,WACV|wx.BOTTOM,5)        
            if 'Tbar' in UseList[G2frame.hist]['Extinction'][2]:       #skipped for TOF   
                valSizer =wx.BoxSizer(wx.HORIZONTAL)
                valSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Tbar(mm):'),0,WACV)
                tbarVal = wx.TextCtrl(DData,wx.ID_ANY,
                    '%.3f'%(UseList[G2frame.hist]['Extinction'][2]['Tbar']),style=wx.TE_PROCESS_ENTER)
                tbarVal.Bind(wx.EVT_TEXT_ENTER,OnTbarVal)
                tbarVal.Bind(wx.EVT_KILL_FOCUS,OnTbarVal)
                valSizer.Add(tbarVal,0,WACV)
                valSizer.Add(wx.StaticText(DData,wx.ID_ANY,' cos(2ThM):'),0,WACV)
                cos2tm = wx.TextCtrl(DData,wx.ID_ANY,
                    '%.3f'%(UseList[G2frame.hist]['Extinction'][2]['Cos2TM']),style=wx.TE_PROCESS_ENTER)
                cos2tm.Bind(wx.EVT_TEXT_ENTER,OnCos2TM)
                cos2tm.Bind(wx.EVT_KILL_FOCUS,OnCos2TM)
                valSizer.Add(cos2tm,0,WACV)
                extSizer.Add(valSizer,0,WACV)
            val2Sizer =wx.BoxSizer(wx.HORIZONTAL)
            if 'Primary' in UseList[G2frame.hist]['Extinction'][1]:
                Ekey = ['Ep',]
            elif 'Secondary Type II' == UseList[G2frame.hist]['Extinction'][1]:
                Ekey = ['Es',]
            elif 'Secondary Type I' == UseList[G2frame.hist]['Extinction'][1]:
                Ekey = ['Eg',]
            else:
                Ekey = ['Eg','Es']
            for ekey in Ekey:
                Eref = wx.CheckBox(DData,wx.ID_ANY,label=ekey+' : ')
                Eref.SetValue(UseList[G2frame.hist]['Extinction'][2][ekey][1])
                Indx[Eref.GetId()] = [G2frame.hist,ekey]
                Eref.Bind(wx.EVT_CHECKBOX, OnEref)
                val2Sizer.Add(Eref,0,WACV)
                Eval = wx.TextCtrl(DData,wx.ID_ANY,
                    '%10.3e'%(UseList[G2frame.hist]['Extinction'][2][ekey][0]),style=wx.TE_PROCESS_ENTER)
                Indx[Eval.GetId()] = [G2frame.hist,ekey]
                Eval.Bind(wx.EVT_TEXT_ENTER,OnEval)
                Eval.Bind(wx.EVT_KILL_FOCUS,OnEval)
                val2Sizer.Add(Eval,0,WACV)

            extSizer.Add(val2Sizer,0,WACV)
        return extSizer
        
    def BabSizer():
        babSizer = wx.BoxSizer(wx.HORIZONTAL)
        for bab in ['A','U']:
            babRef = wx.CheckBox(DData,wx.ID_ANY,label=' Babinet '+bab+': ')
            babRef.SetValue(UseList[G2frame.hist]['Babinet']['Bab'+bab][1])
            Indx[babRef.GetId()] = [G2frame.hist,bab]
            babRef.Bind(wx.EVT_CHECKBOX, OnBabRef)
            babSizer.Add(babRef,0,WACV)
            babVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.3f'%(UseList[G2frame.hist]['Babinet']['Bab'+bab][0]),style=wx.TE_PROCESS_ENTER)
            Indx[babVal.GetId()] = [G2frame.hist,bab]
            babVal.Bind(wx.EVT_TEXT_ENTER,OnBabVal)
            babVal.Bind(wx.EVT_KILL_FOCUS,OnBabVal)
            babSizer.Add(babVal,0,WACV|wx.BOTTOM,5)
        return babSizer
        
    def OnSelect(event):
        G2frame.hist = keyList[select.GetSelection()]
        selSpin.SetValue(keyList.index(G2frame.hist))
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)
       
    def OnSelSpin(event):
        G2frame.hist = keyList[selSpin.GetValue()]
        selSpin.SetValue(keyList.index(G2frame.hist))
        select.SetSelection(keyList.index(G2frame.hist))
        select.SetFirstItem(keyList.index(G2frame.hist))
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)
        
    def RepaintHistogramInfo():
        G2frame.bottomSizer.DeleteWindows()
        Indx.clear()
        G2frame.bottomSizer = ShowHistogramInfo()
        mainSizer.Add(G2frame.bottomSizer)
        mainSizer.Layout()
        G2frame.dataFrame.Refresh()
        DData.SetVirtualSize(mainSizer.GetMinSize())
        G2frame.dataFrame.SendSizeEvent()
        
    def ShowHistogramInfo():
        try:
            histData = UseList[G2frame.hist]
        except KeyError:
            G2frame.ErrorDialog('Missing data error',
                    G2frame.hist+' not in GSAS-II data tree')
            return
        if 'Use' not in UseList[G2frame.hist]:      #patch
            UseList[G2frame.hist]['Use'] = True
        if 'Babinet' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Babinet'] = {'BabA':[0.0,False],'BabU':[0.0,False]}
        bottomSizer = wx.BoxSizer(wx.VERTICAL)
        showSizer = wx.BoxSizer(wx.HORIZONTAL)
        useData = wx.CheckBox(DData,wx.ID_ANY,label='Use Histogram: '+G2frame.hist+' ?')
        showSizer.Add(useData,0,WACV|wx.TOP|wx.BOTTOM,5)
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(UseList[G2frame.hist]['Use'])
        bottomSizer.Add(showSizer,0,WACV|wx.TOP|wx.BOTTOM,5)
        
        bottomSizer.Add(ScaleSizer(),0,WACV|wx.BOTTOM,5)
            
        if G2frame.hist[:4] == 'PWDR':
            if UseList[G2frame.hist]['Size'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                isoSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetSize),0,WACV)
                bottomSizer.Add(isoSizer)
                bottomSizer.Add(IsoSizer(u'size(\xb5m): ','Size','%.5f',
                    OnSizeVal,OnSizeRef),0,WACV|wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Size'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                uniSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetSize),0,WACV)
                bottomSizer.Add(UniSizer('Size',OnSizeAxis),0,WACV)
                bottomSizer.Add(uniSizer)
                bottomSizer.Add(UniDataSizer(u'size(\xb5m): ','Size','%.5f',OnSizeVal,OnSizeRef)
                    ,0,WACV|wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Size'][0] == 'ellipsoidal':
                ellSizer = wx.BoxSizer(wx.HORIZONTAL)
                ellSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                ellSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                ellSizer.Add(ResetSizer('ellipsoidal',OnResetSize),0,WACV)
                bottomSizer.Add(ellSizer)
                bottomSizer.Add(EllSizeDataSizer(),0,WACV|wx.BOTTOM,5)
            
            if UseList[G2frame.hist]['Mustrain'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                isoSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetStrain),0,WACV)
                bottomSizer.Add(isoSizer)
                bottomSizer.Add(IsoSizer(' microstrain: ','Mustrain','%.1f',
                    OnStrainVal,OnStrainRef),0,WACV|wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Mustrain'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                uniSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetStrain),0,WACV)
                bottomSizer.Add(uniSizer)
                bottomSizer.Add(UniSizer('Mustrain',OnStrainAxis),0,WACV)
                bottomSizer.Add(UniDataSizer('mustrain: ','Mustrain','%.1f',OnStrainVal,OnStrainRef)
                                ,0,WACV|wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
                genSizer = wx.BoxSizer(wx.HORIZONTAL)
                genSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                genSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                genSizer.Add(ResetSizer('generalized',OnResetStrain),0,WACV)
                bottomSizer.Add(genSizer)
                bottomSizer.Add(GenStrainDataSizer(),0,WACV|wx.BOTTOM,5)
            
            bottomSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Hydrostatic/elastic strain:'))
            bottomSizer.Add(HstrainSizer())
                
            poSizer = wx.BoxSizer(wx.VERTICAL)
            POData = UseList[G2frame.hist]['Pref.Ori.']
# patch - add penalty items
            if len(POData) < 7:
                POData.append(['',])
                POData.append(0.1)
            if not POData[6]:
                POData[6] = ['',]
# end patch
            poSizer.Add(PoTopSizer(POData))
            if POData[0] == 'MD':
                poSizer.Add(MDDataSizer(POData))
            else:           #'SH'
                if POData[4]:       #SH order > 0
                    textJ = G2lat.textureIndex(POData[5])
                    poSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Spherical harmonic coefficients: '+'Texture index: %.3f'%(textJ))
                        ,0,WACV|wx.TOP|wx.BOTTOM,5)
                    poSizer.Add(SHDataSizer(POData),0,WACV|wx.TOP|wx.BOTTOM,5)
                    poSizer.Add(SHPenalty(POData),0,WACV|wx.TOP|wx.BOTTOM,5)
                    
            bottomSizer.Add(poSizer,0,WACV|wx.TOP|wx.BOTTOM,5)
            bottomSizer.Add(ExtSizer(),0,WACV|wx.TOP|wx.BOTTOM,5)
            bottomSizer.Add(BabSizer(),0,WACV|wx.BOTTOM,5)
        elif G2frame.hist[:4] == 'HKLF':
            bottomSizer.Add(SCExtSizer(),0,WACV|wx.TOP|wx.BOTTOM,5)
            bottomSizer.Add(BabSizer(),0,WACV|wx.BOTTOM,5)
    
        return bottomSizer
                
    #DData.DestroyChildren() # bad, deletes scrollbars on Mac!
    if DData.GetSizer():
        DData.GetSizer().Clear(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Histogram data for '+PhaseName+':'),0,WACV)
    if G2frame.hist != '':
        topSizer = wx.FlexGridSizer(1,2,5,5)
        selSizer = wx.BoxSizer(wx.HORIZONTAL)    
        selSpin = wx.SpinButton(DData,size=(20,120),style=wx.SP_VERTICAL|wx.SP_WRAP)
        selSpin.SetValue(keyList.index(G2frame.hist))
        selSpin.SetRange(0,len(keyList)-1)
        selSpin.Bind(wx.EVT_SPIN,OnSelSpin)
        selSizer.Add(selSpin)
        select = wx.ListBox(DData,choices=keyList,style=wx.LB_SINGLE,size=(-1,120))
        select.SetSelection(keyList.index(G2frame.hist))
        select.SetFirstItem(keyList.index(G2frame.hist))
        select.Bind(wx.EVT_LISTBOX,OnSelect)
        selSizer.Add(select,0,WACV)
        topSizer.Add(selSizer)
        if PWDR:
            topSizer.Add(PlotSizer())
        mainSizer.Add(topSizer)       
        G2frame.bottomSizer = ShowHistogramInfo()
        mainSizer.Add(G2frame.bottomSizer)
    elif not keyList:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,'  (This project has no data; use Import to read it)'),
                      0,WACV|wx.TOP,10)
    elif not UseList:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,'  (This phase has no associated data; use appropriate Edit/Add... menu item)'),
                      0,WACV|wx.TOP,10)
    else:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,'  (Strange, how did we get here?)'),
                      0,WACV|wx.TOP,10)
        
    G2phsGUI.SetPhaseWindow(G2frame.dataFrame,DData,mainSizer)
