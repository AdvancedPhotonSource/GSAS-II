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
from __future__ import division, print_function
import wx
import numpy as np
import numpy.linalg as nl
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIplot as G2plt
import GSASIIpwd as G2pwd
import GSASIIphsGUI as G2phG
import GSASIIctrlGUI as G2G
import GSASIIpy3 as G2py3

WACV = wx.ALIGN_CENTER_VERTICAL
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)
BLACK = wx.Colour(0,0,0)
mapDefault = {'MapType':'','RefList':'','GridStep':0.25,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}

################################################################################
##### DData routines
################################################################################        
def UpdateDData(G2frame,DData,data,hist='',Scroll=0):
    '''Display the Diffraction Data associated with a phase
    (items where there is a value for each histogram and phase)

    :param wx.frame G2frame: the main GSAS-II frame object
    :param wx.ScrolledWindow DData: notebook page to be used for the display
    :param dict data: all the information on the phase in a dictionary
    :param str hist: histogram name
    :param int Scroll: previous scroll position

    '''
    def PlotSizer():

        def OnPlotSel(event):
            Obj = event.GetEventObject()
            generalData['Data plot type'] = Obj.GetStringSelection()
            G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
            wx.CallLater(100,UpdateDData,G2frame,DData,data,G2frame.hist)
            
        def OnPOhkl(event):
            event.Skip()
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
            
        def OnProj(event):
            Obj = event.GetEventObject()
            generalData['3Dproj'] = Obj.GetValue()
            G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        
        plotSizer = wx.BoxSizer(wx.VERTICAL)
        choice = ['None','Mustrain','Size','Preferred orientation','St. proj. Inv. pole figure','Eq. area Inv. pole figure']
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
        elif generalData['Data plot type'] in ['Mustrain','Size']:
            projSizer = wx.BoxSizer(wx.HORIZONTAL)
            projSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Show projections for: '),0,WACV)
            proj = ['','x','y','z','xy','xz','yz','xyz']
            projType = wx.ComboBox(DData,wx.ID_ANY,value=generalData['3Dproj'],choices=proj,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            projType.Bind(wx.EVT_COMBOBOX, OnProj)
            projSizer.Add(projType,0,WACV)
            plotSizer.Add(projSizer)            
        return plotSizer
       
    def ScaleSizer():
        
        def OnScaleRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Scale'][1] = Obj.GetValue()
        def onChangeFraction(invalid,value,tc):
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
        scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
        if 'PWDR' in G2frame.hist:
            scaleRef = wx.CheckBox(DData,wx.ID_ANY,label=' Phase fraction: ')
        elif 'HKLF' in G2frame.hist:
            scaleRef = wx.CheckBox(DData,wx.ID_ANY,label=' Scale factor: ')                
        scaleRef.SetValue(UseList[G2frame.hist]['Scale'][1])
        scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        scaleSizer.Add(scaleRef,0,WACV|wx.LEFT,5)
        scaleVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Scale'],0,
            xmin=0.,nDig=(10,4),typeHint=float,OnLeave=onChangeFraction)
        scaleSizer.Add(scaleVal,0,WACV)
        if 'PWDR' in G2frame.hist and generalData['Type'] != 'magnetic':
            wtSum = G2pwd.PhaseWtSum(G2frame,G2frame.hist)
            if wtSum and UseList[G2frame.hist]['Use']:
                weightFr = UseList[G2frame.hist]['Scale'][0]*generalData['Mass']/wtSum
                scaleSizer.Add(wx.StaticText(DData,label=' Wt. fraction: %.3f'%(weightFr)),0,WACV)
        return scaleSizer
        
    def OnLGmixRef(event):
        Obj = event.GetEventObject()
        hist,name = Indx[Obj.GetId()]
        UseList[G2frame.hist][name][2][2] = Obj.GetValue()
        
    def OnSizeType(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Size'][0] = Obj.GetValue()
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
        
    def OnSizeRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[G2frame.hist]['Size'][0] == 'ellipsoidal':
            UseList[G2frame.hist]['Size'][5][pid] = Obj.GetValue()                
        else:
            UseList[G2frame.hist]['Size'][2][pid] = Obj.GetValue()
        
    def OnStrainType(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Mustrain'][0] = Obj.GetValue()
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
        
    def OnStrainRef(event):
        Obj = event.GetEventObject()
        hist,pid = Indx[Obj.GetId()]
        if UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
            UseList[G2frame.hist]['Mustrain'][5][pid] = Obj.GetValue()
        else:
            UseList[G2frame.hist]['Mustrain'][2][pid] = Obj.GetValue()
        
    def OnStrainAxis(event):
        event.Skip()
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
        wx.CallAfter(G2plt.PlotSizeStrainPO,G2frame,data,hist)
        
    def OnResetStrain(event):
        Obj = event.GetEventObject()
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
            UseList[item]['Mustrain'][4] = vals
        G2plt.PlotSizeStrainPO(G2frame,data,item)
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
    def OnPOAxis(event):
        event.Skip()
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
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

    def OnPOType(event):
        Obj = event.GetEventObject()
        if 'March' in Obj.GetValue():
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'MD'
        else:
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'SH'
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

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
        
    def checkAxis(axis):
        if not np.any(np.array(axis)):
            return False
        return axis
        
    def OnNewValue(invalid,value,tc):
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
            
    def OnNewValueReDraw(invalid,value,tc):
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
    def TopSizer(name,choices,parm,OnType):
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(DData,wx.ID_ANY,name),0,WACV)
        sizeType = wx.ComboBox(DData,wx.ID_ANY,value=UseList[G2frame.hist][parm][0],choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        sizeType.Bind(wx.EVT_COMBOBOX, OnType)
        topSizer.Add(sizeType,0,WACV|wx.BOTTOM,5)
        return topSizer
        
    def LGmixSizer(name,Limits,OnRef):
        lgmixSizer = wx.BoxSizer(wx.HORIZONTAL)
        lgmixRef = wx.CheckBox(DData,wx.ID_ANY,label='LGmix')
        lgmixRef.thisown = False
        lgmixRef.SetValue(UseList[G2frame.hist][name][2][2])
        Indx[lgmixRef.GetId()] = [G2frame.hist,name]
        lgmixRef.Bind(wx.EVT_CHECKBOX, OnRef)
        lgmixSizer.Add(lgmixRef,0,WACV|wx.LEFT,5)
        lgmixVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist][name][1],2,
            nDig=(10,3),xmin=Limits[0],xmax=Limits[1])
        lgmixSizer.Add(lgmixVal,0,WACV|wx.LEFT,5)
        return lgmixSizer
                    
    def ResetSizer(name,OnReset):
        resetSizer = wx.BoxSizer(wx.HORIZONTAL)
        reset = wx.Button(DData,wx.ID_ANY,label='Reset?')
        reset.thisown = False
        Indx[reset.GetId()] = [G2frame.hist,name]
        reset.Bind(wx.EVT_BUTTON,OnReset)
        resetSizer.Add(reset,0,WACV)
        return resetSizer
        
    def IsoSizer(name,parm,fmt,Limits,OnRef):
        isoSizer = wx.BoxSizer(wx.HORIZONTAL)
        sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=name)
        sizeRef.thisown = False
        sizeRef.SetValue(UseList[G2frame.hist][parm][2][0])
        Indx[sizeRef.GetId()] = [G2frame.hist,0]
        sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
        isoSizer.Add(sizeRef,0,WACV|wx.LEFT,5)
        sizeVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist][parm][1],0,
            nDig=fmt,xmin=Limits[0],xmax=Limits[1],OnLeave=OnNewValue)
        isoSizer.Add(sizeVal,0,WACV)
        return isoSizer
        
    def UniSizer(parm,OnAxis):
        uniSizer = wx.BoxSizer(wx.HORIZONTAL)
        uniSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Unique axis, H K L: '),0,WACV)
        h,k,l = UseList[G2frame.hist][parm][3]
        Axis = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        Axis.Bind(wx.EVT_TEXT_ENTER,OnAxis)
        Axis.Bind(wx.EVT_KILL_FOCUS,OnAxis)
        uniSizer.Add(Axis,0,WACV|wx.LEFT,5)
        return uniSizer
        
    def UniDataSizer(parmName,parm,fmt,Limits,OnRef):
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        parms = zip([' Equatorial '+parmName,' Axial '+parmName],
            UseList[G2frame.hist][parm][2],range(2))
        for Pa,ref,Id in parms:
            sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [G2frame.hist,Id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
            dataSizer.Add(sizeRef,0,WACV|wx.LEFT,5)
            sizeVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist][parm][1],
                Id,fmt,xmin=Limits[0],xmax=Limits[1],OnLeave=OnNewValue)
            dataSizer.Add(sizeVal,0,WACV)
        return dataSizer

    def EllSizeDataSizer():
        parms = zip(['S11','S22','S33','S12','S13','S23'],UseList[G2frame.hist]['Size'][4],
            UseList[G2frame.hist]['Size'][5],range(6))
        dataSizer = wx.BoxSizer(wx.VERTICAL)
        matrixSizer = wx.FlexGridSizer(0,6,5,5)
        Sij = []
        for Pa,val,ref,Id in parms:
            sizeRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            sizeRef.thisown = False
            sizeRef.SetValue(ref)
            Indx[sizeRef.GetId()] = [G2frame.hist,Id]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
            matrixSizer.Add(sizeRef,0,WACV)
            if Id < 3:
                sizeVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Size'][4],
                    Id,nDig=(10,3),xmin=0.,xmax=4.,OnLeave=OnNewValueReDraw)
            else:
                sizeVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Size'][4],
                    Id,nDig=(10,3),OnLeave=OnNewValueReDraw)
            # Create Sij matrix
            Sij += [val]
            matrixSizer.Add(sizeVal,0,WACV)
        dataSizer.Add(matrixSizer, 0)
        Esize,Rsize = nl.eigh(G2lat.U6toUij(np.asarray(Sij)))
        lengths = Esize
        G,g = G2lat.cell2Gmat(data['General']['Cell'][1:7])       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        hkls = [x/(sum(x**2)**0.5) for x in np.dot(Rsize, GA)]
        Ids = np.argsort(lengths)
        dataSizer.Add(wx.StaticText(DData,label=' Principal ellipsoid components:'),0)
        compSizer = wx.FlexGridSizer(3,3,5,5)
        Axes = [' Short Axis:',' Middle Axis:',' Long Axis:']
        for Id in Ids:
            compSizer.Add(wx.StaticText(DData,label=Axes[Id]),0,WACV)
            compSizer.Add(wx.StaticText(DData,label='(%.3f, %.3f, %.3f) '%(hkls[Id][0], hkls[Id][1], hkls[Id][2])),0,WACV)
            compSizer.Add(wx.StaticText(DData,label='Length: %.3f'%lengths[Id]),0,WACV)
        dataSizer.Add(compSizer)
        return dataSizer
        
    def GenStrainDataSizer():
        Snames = G2spc.MustrainNames(SGData)
        numb = len(Snames)
        onumb = len(UseList[G2frame.hist]['Mustrain'][4])
        while onumb < numb:
            UseList[G2frame.hist]['Mustrain'][4].append(0.0)
            UseList[G2frame.hist]['Mustrain'][5].append(False)
            onumb += 1
        muMean = G2spc.MuShklMean(SGData,Amat,UseList[G2frame.hist]['Mustrain'][4][:numb])
        parms = zip(Snames,UseList[G2frame.hist]['Mustrain'][5],range(numb))
        dataSizer = wx.FlexGridSizer(0,6,5,5)
        for Pa,ref,Id in parms:
            strainRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            strainRef.thisown = False
            strainRef.SetValue(ref)
            Indx[strainRef.GetId()] = [G2frame.hist,Id]
            strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
            dataSizer.Add(strainRef,0,WACV)
            strainVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Mustrain'][4],
                Id,nDig=(10,2),OnLeave=OnNewValueReDraw)
            dataSizer.Add(strainVal,0,WACV)
        dataSizer.Add(wx.StaticText(DData,label=' Mean mustrain %.1f'%muMean),0,WACV)
        return dataSizer

    def HstrainSizer():
        
        def OnHstrainRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            UseList[G2frame.hist]['HStrain'][1][pid] = Obj.GetValue()
            
        hSizer = wx.BoxSizer(wx.VERTICAL)
        hstrainSizer = wx.FlexGridSizer(0,6,5,5)
        Hsnames = G2spc.HStrainNames(SGData)
        parms = zip(Hsnames,UseList[G2frame.hist]['HStrain'][1],range(len(Hsnames)))
        allzero = True
        for Pa,ref,Id in parms:
            hstrainRef = wx.CheckBox(DData,wx.ID_ANY,label=Pa)
            hstrainRef.thisown = False
            hstrainRef.SetValue(ref)
            Indx[hstrainRef.GetId()] = [G2frame.hist,Id]
            hstrainRef.Bind(wx.EVT_CHECKBOX, OnHstrainRef)
            hstrainSizer.Add(hstrainRef,0,WACV|wx.LEFT,5)
            hstrainVal = G2G.ValidatedTxtCtrl(DData,
                        UseList[G2frame.hist]['HStrain'][0],Id,nDig=(10,3,'g'),
                        OnLeave=OnNewValueReDraw)
            if abs(UseList[G2frame.hist]['HStrain'][0][Id]) > 1e-8:
                allzero = False
            hstrainSizer.Add(hstrainVal,0,WACV)
        hSizer.Add(hstrainSizer,0)
        if not allzero:   # show Dij shifted unit cell
            DijVals = UseList[G2frame.hist]['HStrain'][0][:]
            # apply the Dij values to the reciprocal cell
            newA = []
            Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
            for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),
                            ['D11','D22','D33','D12','D13','D23']):
                newA.append(Aij + Dijdict.get(lbl,0.0))
            cell = G2lat.A2cell(newA)   # convert back to direct cell
            laue = generalData['SGData']['SGLaue']
            if laue == '2/m':
                laue += generalData['SGData']['SGUniq']
            for cellGUI in G2py3.cellGUIlist:
                if laue in cellGUI[0]:
                    useGUI = cellGUI
                    break
            else:
                return hSizer
            cellstr = ''
            for txt,fmt,ifEdit,Id in zip(*useGUI[2:]):
                if cellstr: cellstr += ", "
                cellstr += txt+fmt.format(cell[Id])
            cellstr += ', Vol = {:.3f}'.format(G2lat.calc_V(newA))
            hSizer.Add(wx.StaticText(DData,wx.ID_ANY,'     '+cellstr),0)
        return hSizer
        
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
        poSizer.Add(poRef,0,WACV|wx.LEFT,5)
        poVal = G2G.ValidatedTxtCtrl(DData,POData,1,nDig=(10,3),typeHint=float,xmin=0.)
        poSizer.Add(poVal,0,WACV)
        poSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Unique axis, H K L: '),0,WACV)
        h,k,l =POData[3]
        poAxis = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
        poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOAxis)
        poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOAxis)
        poSizer.Add(poAxis,0,WACV)
        return poSizer
        
    def SHDataSizer(POData):
        
        ODFSizer = wx.FlexGridSizer(0,8,2,2)
        ODFkeys = list(POData[5].keys())
        ODFkeys.sort()
        for odf in ODFkeys:
            ODFSizer.Add(wx.StaticText(DData,wx.ID_ANY,odf),0,WACV)
            ODFval = G2G.ValidatedTxtCtrl(DData,POData[5],odf,nDig=(8,3),typeHint=float,OnLeave=OnNewValue)
            ODFSizer.Add(ODFval,0,WACV|wx.LEFT,5)
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
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
        A = G2lat.cell2A(generalData['Cell'][1:7])
        hkls = G2lat.GenPfHKLs(10,SGData,A)    
        shPenalty = wx.BoxSizer(wx.HORIZONTAL)
        shPenalty.Add(wx.StaticText(DData,wx.ID_ANY,' Negative MRD penalty list: '),0,WACV)
        shPenalty.Add(wx.ComboBox(DData,value=POData[6][0],choices=POData[6],
            style=wx.CB_DROPDOWN),0,WACV|wx.LEFT,5)
        hklList = wx.Button(DData,label='Select penalty hkls')
        hklList.Bind(wx.EVT_BUTTON,OnHKLList)
        shPenalty.Add(hklList,0,WACV)
        shPenalty.Add(wx.StaticText(DData,wx.ID_ANY,' Zero MRD tolerance: '),0,WACV)
        shToler = G2G.ValidatedTxtCtrl(DData,POData,7,nDig=(10,2),typeHint=float)
        shPenalty.Add(shToler,0,WACV)
        return shPenalty    
        
    def ExtSizer(Type):
        
        def OnSCExtType(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            UseList[item[0]]['Extinction'][item[1]] = Obj.GetValue()
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
                
        def OnEref(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            UseList[item[0]]['Extinction'][2][item[1]][1] = Obj.GetValue()
    
        def OnExtRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Extinction'][1] = Obj.GetValue()
            
        if Type == 'HKLF':
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
                extSizer.Add(typeSizer,0)
            else:
                extSizer.Add(typeSizer,0,wx.BOTTOM,5)        
                if 'Tbar' in UseList[G2frame.hist]['Extinction'][2]:       #skipped for TOF   
                    valSizer =wx.BoxSizer(wx.HORIZONTAL)
                    valSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Tbar(mm):'),0,WACV)
                    tbarVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Extinction'][2],'Tbar',
                        xmin=0.,nDig=(10,3),typeHint=float)
                    valSizer.Add(tbarVal,0,WACV)
                    valSizer.Add(wx.StaticText(DData,wx.ID_ANY,' cos(2ThM):'),0,WACV)
                    cos2tm = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Extinction'][2],'Cos2TM',
                        xmin=0.,xmax=1.,nDig=(10,3),typeHint=float)
                    valSizer.Add(cos2tm,0,WACV)
                    extSizer.Add(valSizer,0)
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
                    val2Sizer.Add(Eref,0,WACV|wx.LEFT,5)
                    Eval = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Extinction'][2][ekey],0,
                        xmin=0.,nDig=(10,3,'g'),typeHint=float)
                    val2Sizer.Add(Eval,0,WACV)
                extSizer.Add(val2Sizer,0)
        else:   #PWDR
            extSizer = wx.BoxSizer(wx.HORIZONTAL)
            extRef = wx.CheckBox(DData,wx.ID_ANY,label=' Extinction: ')
            extRef.SetValue(UseList[G2frame.hist]['Extinction'][1])
            extRef.Bind(wx.EVT_CHECKBOX, OnExtRef)
            extSizer.Add(extRef,0,WACV|wx.LEFT,5)
            extVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Extinction'],0,
                xmin=0.,nDig=(10,2),typeHint=float)
            extSizer.Add(extVal,0,WACV)

        return extSizer
        
    def BabSizer():
        
        def OnBabRef(event):
            Obj = event.GetEventObject()
            item,bab = Indx[Obj.GetId()]
            UseList[item]['Babinet']['Bab'+bab][1] = Obj.GetValue()
        
        babSizer = wx.BoxSizer(wx.HORIZONTAL)
        for bab in ['A','U']:
            babRef = wx.CheckBox(DData,wx.ID_ANY,label=' Babinet '+bab+': ')
            babRef.SetValue(UseList[G2frame.hist]['Babinet']['Bab'+bab][1])
            Indx[babRef.GetId()] = [G2frame.hist,bab]
            babRef.Bind(wx.EVT_CHECKBOX, OnBabRef)
            babSizer.Add(babRef,0,WACV|wx.LEFT,5)
            babVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Babinet']['Bab'+bab],0,
                nDig=(10,3),xmin=0.,typeHint=float)
            babSizer.Add(babVal,0,WACV)
        return babSizer
        
    def FlackSizer():
        
        def OnFlackRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Flack'][1] = Obj.GetValue()
                
        flackSizer = wx.BoxSizer(wx.HORIZONTAL)
        flackRef = wx.CheckBox(DData,wx.ID_ANY,label=' Flack parameter: ')
        flackRef.SetValue(UseList[G2frame.hist]['Flack'][1])
        flackRef.Bind(wx.EVT_CHECKBOX, OnFlackRef)
        flackSizer.Add(flackRef,0,WACV|wx.LEFT,5)
        flackVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Flack'],0,nDig=(10,3),typeHint=float)
        flackSizer.Add(flackVal,0,WACV)
        return flackSizer
        
    def DispSizer():
        
        def OnDispRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Layer Disp'][1] = Obj.GetValue()
                
        dispSizer = wx.BoxSizer(wx.HORIZONTAL)
        dispRef = wx.CheckBox(DData,wx.ID_ANY,label=u' Layer displacement (\xb5m): ')
        dispRef.SetValue(UseList[G2frame.hist]['Layer Disp'][1])
        dispRef.Bind(wx.EVT_CHECKBOX, OnDispRef)
        dispSizer.Add(dispRef,0,WACV|wx.LEFT,5)
        dispSizer.Add(G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Layer Disp'],0,nDig=(10,2),typeHint=float),0,WACV)
        return dispSizer
        
    def twinSizer():
        
        def OnAddTwin(event):
            twinMat = np.array([[-1,0,0],[0,-1,0],[0,0,-1]])    #inversion by default
            twinVal = 0.0
            UseList[G2frame.hist]['Twins'].append([twinMat,twinVal])
            nNonM = UseList[G2frame.hist]['Twins'][0][1][2]
            for i in range(nNonM):
                UseList[G2frame.hist]['Twins'].append([False,0.0])
            addtwin.SetValue(False)
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
        def OnMat(event):
            event.Skip()
            Obj = event.GetEventObject()
            it,im = Indx[Obj.GetId()]
            newMat = Obj.GetValue().split()
            try:
                uvw = [int(newMat[i]) for i in range(3)]
            except ValueError:
                uvw = UseList[G2frame.hist]['Twins'][it][0][im]
            UseList[G2frame.hist]['Twins'][it][0][im] = uvw
            Obj.SetValue('%3d %3d %3d'%(uvw[0],uvw[1],uvw[2]))
            
        def OnTwinVal(invalid,value,tc):
            it = Indx[tc.GetId()]
            sumTw = 0.
            for it,twin in enumerate(UseList[G2frame.hist]['Twins']):
                if it:
                    sumTw += twin[1]
            UseList[G2frame.hist]['Twins'][0][1][0] = 1.-sumTw
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
        def OnTwinRef(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Twins'][0][1][1] = Obj.GetValue()
            
        def OnTwinInv(event):
            Obj = event.GetEventObject()
            it = Indx[Obj.GetId()]
            UseList[G2frame.hist]['Twins'][it][0] = Obj.GetValue()
                        
        def OnTwinDel(event):
            Obj = event.GetEventObject()
            it = Indx[Obj.GetId()]
            nNonM = UseList[G2frame.hist]['Twins'][0][1][2]
            for i in range(nNonM):
                del UseList[G2frame.hist]['Twins'][1+i+it]
            del UseList[G2frame.hist]['Twins'][it]
            sumTw = 0.
            for it,twin in enumerate(UseList[G2frame.hist]['Twins']):
                if it:
                    sumTw += twin[1]
            UseList[G2frame.hist]['Twins'][0][1][0] = 1.-sumTw
            if len(UseList[G2frame.hist]['Twins']) == 1:
                UseList[G2frame.hist]['Twins'][0][1][1] = False
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))           
            
        nTwin = len(UseList[G2frame.hist]['Twins'])
        twinsizer = wx.BoxSizer(wx.VERTICAL)
        topsizer = wx.BoxSizer(wx.HORIZONTAL)          
        topsizer.Add(wx.StaticText(DData,wx.ID_ANY,' Merohedral twins: '),0,WACV)
        #temporary - add twin not allowed if nonmerohedral twins present
#        if nTwin == 1 or 'bool' not in str(type(UseList[G2frame.hist]['Twins'][1][0])):
        addtwin = wx.CheckBox(DData,wx.ID_ANY,label=' Add Twin Law')
        addtwin.Bind(wx.EVT_CHECKBOX, OnAddTwin)
        topsizer.Add(addtwin,0,WACV|wx.LEFT,5)
        twinsizer.Add(topsizer)
        Indx = {}
        if nTwin > 1:
            for it,Twin in enumerate(UseList[G2frame.hist]['Twins']):
                twinMat,twinVal = Twin
                matSizer = wx.BoxSizer(wx.HORIZONTAL)
                if it:
                    Style = wx.TE_PROCESS_ENTER
                    TwVal = Twin[1]
                else:
                    Style = wx.TE_READONLY
                    TwVal = Twin[1][0]
                if 'bool' not in str(type(Twin[0])):
                    matSizer.Add(wx.StaticText(DData,-1,' Twin Law: '),0,WACV|wx.LEFT,5)
                    for im,Mat in enumerate(twinMat):
                        mat = wx.TextCtrl(DData,wx.ID_ANY,'%3d %3d %3d'%(Mat[0],Mat[1],Mat[2]),
                            style=Style)
                        if it:
                            Indx[mat.GetId()] = [it,im]
                            mat.Bind(wx.EVT_TEXT_ENTER,OnMat)
                            mat.Bind(wx.EVT_KILL_FOCUS,OnMat)
                        else:
                            mat.SetBackgroundColour(VERY_LIGHT_GREY)
                        matSizer.Add(mat,0,WACV|wx.LEFT,5)
                else:
                    matSizer.Add(wx.StaticText(DData,-1,' Nonmerohedral twin component %d: '%(it)),0,WACV|wx.LEFT,5)
                    if not SGData['SGInv']:
                        twinv = wx.CheckBox(DData,wx.ID_ANY,label=' Use enantiomorph?')
                        twinv.SetValue(Twin[0])
                        Indx[twinv.GetId()] = it
                        twinv.Bind(wx.EVT_CHECKBOX, OnTwinInv)
                        matSizer.Add(twinv,0,WACV)
                twinsizer.Add(matSizer,0,wx.LEFT,5)
                valSizer = wx.BoxSizer(wx.HORIZONTAL)
                valSizer.Add(wx.StaticText(DData,-1,label=' Twin element fraction:'),0,WACV)
                if it:
                    twinval = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Twins'][it],1,nDig=(10,3),
                        xmin=0.,xmax=1.,typeHint=float,OnLeave=OnTwinVal)
                    Indx[twinval.GetId()] = it
                else:
                    twinval = wx.TextCtrl(DData,-1,'%.3f'%(TwVal),style=Style)
                    twinval.SetBackgroundColour(VERY_LIGHT_GREY)
                valSizer.Add(twinval,0,WACV)
                if it and 'bool' not in str(type(Twin[0])):
                    twindel = wx.CheckBox(DData,wx.ID_ANY,label=' Delete?')
                    Indx[twindel.GetId()] = it
                    twindel.Bind(wx.EVT_CHECKBOX, OnTwinDel)
                    valSizer.Add(twindel,0,WACV)
                elif not it:
                    twinref = wx.CheckBox(DData,wx.ID_ANY,label=' Refine?')
                    twinref.SetValue(Twin[1][1])
                    twinref.Bind(wx.EVT_CHECKBOX, OnTwinRef)
                    valSizer.Add(twinref,0,WACV)
                twinsizer.Add(valSizer,0,wx.LEFT,5)
        return twinsizer
        
    def OnSelect(event):
        G2frame.hist = G2frame.dataWindow.HistsInPhase[DData.select.GetSelection()]
        oldFocus = wx.Window.FindFocus()
        G2plt.PlotSizeStrainPO(G2frame,data,G2frame.hist)
        wx.CallLater(100,RepaintHistogramInfo)
        if oldFocus: wx.CallAfter(oldFocus.SetFocus)
       
    def RepaintHistogramInfo(Scroll=0):
        if 'phoenix' in wx.version():
            G2frame.bottomSizer.Clear(True)
            # deal with case where this is called after another tree item has been selected
            try:
                DData.Shown
            except RuntimeError:
                if GSASIIpath.GetConfigValue('debug'):
                    print('DBG: DData window deleted. Ignoring RepaintHistogramInfo, forcing redraw')
                # Repaint called while DData window deleted, force redraw of entire window
                import GSASIIdataGUI
                G2frame.PickIdText = ''
                wx.CallLater(100,GSASIIdataGUI.SelectDataTreeItem,G2frame,G2frame.GPXtree.Selection)
                return
        else:
            # deal with case where this is called after another tree item has been selected
            if DData.__class__ is  not wx._windows.ScrolledWindow:
                # fix bug where this is called after the Window is deleted
                return
            G2frame.bottomSizer.DeleteWindows()
        Indx.clear()
        G2frame.bottomSizer = ShowHistogramInfo()
        mainSizer.Add(G2frame.bottomSizer)
        mainSizer.Layout()
        G2frame.dataWindow.Refresh()
        DData.SetVirtualSize(mainSizer.GetMinSize())
        DData.Scroll(0,Scroll)
        G2frame.dataWindow.SendSizeEvent()
        
    def ShowHistogramInfo():
        '''This creates a sizer with all the information pulled out from the Phase/data dict
        '''
        
        def OnUseData(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Use'] = Obj.GetValue()
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

        def OnLeBail(event):
            Obj = event.GetEventObject()
            if not UseList[G2frame.hist]['LeBail']:
                UseList[G2frame.hist]['newLeBail'] = True
                Obj.SetLabel('Do new LeBail extraction?')
            UseList[G2frame.hist]['LeBail'] = Obj.GetValue()

        def OnResetSize(event):
            Obj = event.GetEventObject()
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
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))
            
        def OnSizeAxis(event):            
            event.Skip()
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
            
        def OnFixVals(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Fix FXU'] = Obj.GetValue()

        if G2frame.hist not in UseList:                
            G2frame.ErrorDialog('Missing data error',
                    G2frame.hist+' not in GSAS-II data tree')
            return
#patch
        if 'Use' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Use'] = True
        if 'LeBail' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['LeBail'] = False
        if 'newLeBail' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['newLeBail'] = True
        if 'Babinet' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Babinet'] = {'BabA':[0.0,False],'BabU':[0.0,False]}
        if 'Fix FXU' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Fix FXU'] = ' '
        if 'Flack' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Flack'] = [0.0,False]
        if 'Twins' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Twins'] = [[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False]],]
        if 'Layer Disp' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Layer Disp'] = [0.0,False]
#end patch
        bottomSizer = wx.BoxSizer(wx.VERTICAL)
        useBox = wx.BoxSizer(wx.HORIZONTAL)
        useData = wx.CheckBox(DData,wx.ID_ANY,label='Use Histogram: '+G2frame.hist+' ?')
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(UseList[G2frame.hist]['Use'])
        useBox.Add(useData,0,WACV)
        if not generalData['doPawley'] and 'PWDR' in G2frame.hist[:4]:
            lbLabel = 'Redo LeBail extraction?   '
            if UseList[G2frame.hist]['newLeBail']:
                lbLabel = 'Do new LeBail extraction?'
            lebail = wx.CheckBox(DData,wx.ID_ANY,label=lbLabel)
            lebail.Bind(wx.EVT_CHECKBOX, OnLeBail)
            lebail.SetValue(UseList[G2frame.hist]['LeBail'])
            useBox.Add(lebail,0,WACV)
            if UseList[G2frame.hist]['LeBail']:
                G2frame.SetStatusText('To reset LeBail, cycle LeBail check box.',1)
        bottomSizer.Add(useBox,0,wx.TOP|wx.BOTTOM|wx.LEFT,5)
        fixBox = wx.BoxSizer(wx.HORIZONTAL)
        parmChoice = [' ','X','XU','U','F','FX','FXU','FU']
        if generalData['Type'] == 'magnetic':
            parmChoice += ['M','MX','MXU','MU','MF','MFX','MFXU','MFU']
        fixBox.Add(wx.StaticText(DData,label=' In sequential refinement, fix these in '+generalData['Name']+' for this histogram: '),0,WACV)
        fixVals = wx.ComboBox(DData,value=UseList[G2frame.hist]['Fix FXU'],choices=parmChoice,
            style=wx.CB_DROPDOWN)
        fixVals.Bind(wx.EVT_COMBOBOX,OnFixVals)
        fixBox.Add(fixVals,0,WACV)
        bottomSizer.Add(fixBox)
        #TODO - put Sequential refinement fix F? fix X? fix U? CheckBox here
        
        bottomSizer.Add(ScaleSizer(),0,wx.BOTTOM,5)
            
        if G2frame.hist[:4] == 'PWDR':
            if UseList[G2frame.hist]['Size'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                isoSizer.Add(LGmixSizer('Size',[0.,1.],OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetSize),0,WACV)
                bottomSizer.Add(isoSizer)
                bottomSizer.Add(IsoSizer(u'size(\xb5m): ','Size',(10,4),[0.,4.],OnSizeRef),0,wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Size'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                uniSizer.Add(LGmixSizer('Size',[0.,1.],OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetSize),0,WACV)
                bottomSizer.Add(UniSizer('Size',OnSizeAxis),0)
                bottomSizer.Add(uniSizer)
                bottomSizer.Add(UniDataSizer(u'size(\xb5m): ','Size',(10,3),[0.,4.],OnSizeRef),0,wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Size'][0] == 'ellipsoidal':
                ellSizer = wx.BoxSizer(wx.HORIZONTAL)
                ellSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                ellSizer.Add(LGmixSizer('Size',[0.,1.],OnLGmixRef))
                ellSizer.Add(ResetSizer('ellipsoidal',OnResetSize),0,WACV)
                bottomSizer.Add(ellSizer)
                bottomSizer.Add(EllSizeDataSizer(),0,wx.BOTTOM,5)
            
            if UseList[G2frame.hist]['Mustrain'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                isoSizer.Add(LGmixSizer('Mustrain',[0.,1.],OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetStrain),0,WACV)
                bottomSizer.Add(isoSizer)
                bottomSizer.Add(IsoSizer(' microstrain: ','Mustrain',(10,1),[0.,1.e5],OnStrainRef),0,wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Mustrain'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                uniSizer.Add(LGmixSizer('Mustrain',[0.,1.],OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetStrain),0,WACV)
                bottomSizer.Add(uniSizer)
                bottomSizer.Add(UniSizer('Mustrain',OnStrainAxis),0)
                bottomSizer.Add(UniDataSizer('mustrain: ','Mustrain',(10,1),[0.,1.e5],OnStrainRef),0,wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Mustrain'][0] == 'generalized':
                genSizer = wx.BoxSizer(wx.HORIZONTAL)
                genSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                    'Mustrain',OnStrainType),0,WACV)
                genSizer.Add(LGmixSizer('Mustrain',[0.,1.],OnLGmixRef))
                genSizer.Add(ResetSizer('generalized',OnResetStrain),0,WACV)
                bottomSizer.Add(genSizer)
                bottomSizer.Add(GenStrainDataSizer(),0,WACV|wx.BOTTOM,5)
            
            bottomSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Hydrostatic/elastic strain:'))
            bottomSizer.Add(HstrainSizer())
            bottomSizer.Add(DispSizer())
                
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
                    poSizer.Add(SHDataSizer(POData),0,wx.TOP|wx.BOTTOM,5)
                    poSizer.Add(SHPenalty(POData),0,wx.TOP|wx.BOTTOM,5)
                    
            bottomSizer.Add(poSizer,0,wx.TOP|wx.BOTTOM,5)
            bottomSizer.Add(ExtSizer('PWDR'),0,wx.TOP|wx.BOTTOM,5)
            if generalData['Type'] != 'magnetic': 
                bottomSizer.Add(BabSizer(),0,wx.BOTTOM,5)
        elif G2frame.hist[:4] == 'HKLF':
            bottomSizer.Add(ExtSizer('HKLF'),0,wx.BOTTOM,5)
            bottomSizer.Add(BabSizer(),0,wx.BOTTOM,5)
            if not SGData['SGInv'] and len(UseList[G2frame.hist]['Twins']) < 2:
                bottomSizer.Add(FlackSizer(),0,wx.BOTTOM,5)
            bottomSizer.Add(twinSizer(),0,wx.BOTTOM,5)
        return bottomSizer

    ######################################################################
    ### Beginning of UpdateDData execution here
    ######################################################################
    G2frame.SetStatusText('',1)
    keyList = G2frame.GetHistogramNames(['PWDR','HKLF'])
    UseList = data['Histograms']
    if UseList:
        G2frame.dataWindow.DataMenu.Enable(G2G.wxID_DATADELETE,True)
        for item in G2frame.Refine: item.Enable(True)
    else:
        G2frame.dataWindow.DataMenu.Enable(G2G.wxID_DATADELETE,False)
        for item in G2frame.Refine: item.Enable(False)
    # make a list of histograms (any type) used in this phase, ordered as in tree
    G2frame.dataWindow.HistsInPhase = [name for name in keyList if name in UseList]
    generalData = data['General']
    cell = generalData['Cell'][1:]
    Amat,Bmat = G2lat.cell2AB(cell[:6])
    PhaseName = generalData['Name']       
    SGData = generalData['SGData']
    if len(G2frame.dataWindow.HistsInPhase) == 0: # no associated histograms, nothing to display here
        G2frame.hist = ''
    elif hist and hist in G2frame.dataWindow.HistsInPhase: # something was input as a selection as an argument
        G2frame.hist = hist
    elif (not G2frame.hist) or (G2frame.hist not in G2frame.dataWindow.HistsInPhase): # no or bad selection but have data, take the first
        G2frame.hist = G2frame.dataWindow.HistsInPhase[0]
    Indx = {}
    
    if DData.GetSizer():
        try:
            if hasattr(DData,'select'): 
                DData.select.Unbind(wx.EVT_LISTBOX)  # remove binding to avoid event on Linux
        except:
            pass
        DData.GetSizer().Clear(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Histogram data for '+PhaseName+':'),0,wx.LEFT,5)
    if G2frame.hist:
        topSizer = wx.FlexGridSizer(1,2,5,5)
        DData.select = wx.ListBox(DData,choices=G2frame.dataWindow.HistsInPhase,
                            style=wx.LB_SINGLE,size=(-1,160))
        DData.select.SetSelection(G2frame.dataWindow.HistsInPhase.index(G2frame.hist))
        DData.select.SetFirstItem(G2frame.dataWindow.HistsInPhase.index(G2frame.hist))
        DData.select.Bind(wx.EVT_LISTBOX,OnSelect)
        topSizer.Add(DData.select,0,WACV|wx.LEFT,5)
        if any(['PWDR' in item for item in keyList]):
            topSizer.Add(PlotSizer())
        mainSizer.Add(topSizer)       
        G2frame.bottomSizer = ShowHistogramInfo()
        mainSizer.Add(G2frame.bottomSizer)
    elif not keyList:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,
            '  (This project has no data; use Import to read it)'),0,wx.TOP,10)
    elif not UseList in G2frame.dataWindow.HistsInPhase:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,
            '  (This phase has no associated data; use appropriate Edit/Add... menu item)'),0,wx.TOP,10)
    else:
        mainSizer.Add(wx.StaticText(DData,wx.ID_ANY,'  (Strange, how did we get here?)'),0,wx.TOP,10)
        
    G2phG.SetPhaseWindow(DData,mainSizer,Scroll=Scroll)
