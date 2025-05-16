# -*- coding: utf-8 -*-
#GSASII - phase data display routines
'''Routines for Data tab in Phase dataframe follows.
'''
from __future__ import division, print_function
import copy
import numpy as np
import numpy.linalg as nl
from . import GSASIIpath
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIplot as G2plt
from . import GSASIIpwd as G2pwd
from . import GSASIIphsGUI as G2phG
from . import GSASIIctrlGUI as G2G
from . import GSASIIdataGUI as G2gd
from . import GSASIIfiles as G2fil
from . import GSASIImath as G2mth

try:
    import wx
    WACV = wx.ALIGN_CENTER_VERTICAL
    #VERY_LIGHT_GREY = wx.Colour(235,235,235)
    #WHITE = wx.Colour(255,255,255)
    #BLACK = wx.Colour(0,0,0)
    VERY_LIGHT_GREY = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
    WHITE = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
    BLACK = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT)
except:
    pass
mapDefault = {'MapType':'','RefList':'','GridStep':0.25,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}

def UpdateDData(G2frame,DData,data,hist='',Scroll=0):
    '''Display the Diffraction Data associated with a phase
    (items where there is a value for each histogram and phase)
    Used in the Phase/Data tab or the Hist/Phase tree entry

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

        ifkeV = 'E' in UseList[G2frame.hist].get('Type','')
        plotSizer = wx.BoxSizer(wx.VERTICAL)
        if ifkeV:
            choice = ['None','Mustrain','Size']
        else:
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
        scaleRef.SetValue(bool(UseList[G2frame.hist]['Scale'][1]))
        scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        scaleSizer.Add(scaleRef,0,WACV|wx.LEFT,5)
        scaleVal = G2G.ValidatedTxtCtrl(DData,UseList[G2frame.hist]['Scale'],0,
            xmin=0.,nDig=(10,4,'g'),typeHint=float,OnLeave=onChangeFraction)
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
        '''Respond when the SH order is changed in the GUI.
        Updates the dict with the cylindrical Spherical harmonics 
        coefficients for the specified order, retaining values from 
        the previous dict, if values were already present
        '''
        Obj = event.GetEventObject()
        Order = int(Obj.GetValue())
        UseList[G2frame.hist]['Pref.Ori.'][4] = Order
        UseList[G2frame.hist]['Pref.Ori.'][5] = SetPOCoef(Order)
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

    def OnPOType(event):
        '''Respond when the SH type is changed in the GUI.
        Note that values set that are not used are also not changed.
        '''
        Obj = event.GetEventObject()
        if 'March' in Obj.GetValue():
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'MD'
        else:
            UseList[G2frame.hist]['Pref.Ori.'][0] = 'SH'
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

    def OnPORef(event):
        Obj = event.GetEventObject()
        UseList[G2frame.hist]['Pref.Ori.'][2] = Obj.GetValue()

    def OnAddFixed(event):
        fixedVars = UseList[G2frame.hist].get('FixedSeqVars',[])
        SeqId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, 'Sequential results')
        seqData = G2frame.GPXtree.GetItemPyData(SeqId)
        if G2frame.hist not in seqData:
            print('Strange: '+G2frame.hist+' not found')
            return
        parmDict = seqData[G2frame.hist].get('parmDict',[])
        # narrow down to items w/o a histogram & having float values
        phaseKeys = [i for i in parmDict if ':' in i and i.split(':')[1] == '']
        phaseKeys = [i for i in phaseKeys if type(parmDict[i]) not in (int,str,bool)]
        if len(phaseKeys) == 0: return
        selected = []
        for i,key in enumerate(phaseKeys):
            if key in fixedVars: selected.append(i)
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Choose phase vars to fix for this histogram only',
            'Choose items to edit', phaseKeys,selected=selected)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        UseList[G2frame.hist]['FixedSeqVars'] = [phaseKeys[i] for i in sel]
        wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

    def SetPOCoef(Order):
        '''Sets up a dict with the cylindrical Spherical harmonics 
        coefficients for the specified order. 
        Retains values from the previous dict, if values were already present
        '''
        cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',Order,False)     #cylindrical & no M
        newPOCoef = dict(zip(cofNames,len(cofNames)*[0.]))
        POCoeff = UseList[G2frame.hist]['Pref.Ori.'][5]
        for cofName in POCoeff:
            if cofName in cofNames:
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
        mainSizer = wx.BoxSizer(wx.VERTICAL)
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
        mainSizer.Add(dataSizer)
        mainSizer.Add(wx.StaticText(DData,label=' Mean mustrain %.1f'%muMean)
                          ,0,wx.ALIGN_CENTER_HORIZONTAL)
        return mainSizer

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
                UseList[G2frame.hist]['HStrain'][0],Id,nDig=(10,3,'g'),OnLeave=OnNewValueReDraw)
            if abs(UseList[G2frame.hist]['HStrain'][0][Id]) > 1e-8:
                allzero = False
            hstrainSizer.Add(hstrainVal,0,WACV)
        hSizer.Add(hstrainSizer,0)
        if not allzero:   # show Dij shifted unit cell
            DijVals = UseList[G2frame.hist]['HStrain'][0][:]
            A = G2lat.cell2A(data['General']['Cell'][1:7])
            # apply the Dij values to the reciprocal cell
            newA =  G2lat.AplusDij(A,DijVals,SGData)
            cell = G2lat.A2cell(newA)   # convert back to direct cell
            laue = generalData['SGData']['SGLaue']
            if laue == '2/m':
                laue += generalData['SGData']['SGUniq']
            for cellGUI in G2fil.cellGUIlist:
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
        h,k,l = POData[3]
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
                    POData[6] = [hkls[i].split(':')[0] for i in dlg.GetSelections()]
                    if not POData[6]:
                        POData[6] = ['',]
                else:
                    return
            finally:
                dlg.Destroy()
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,G2frame.hist)
        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Instrument Parameters'))[0]
        reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Reflection Lists'))
        hkls = []
        if generalData['Name'] in reflSets:
            reflData = reflSets[generalData['Name']]
            if 'RefList' in reflData:
                Super = generalData.get('Super',0)
                nprfo = 12
                if 'T' in Inst['Type'][0] or 'B' in Inst['Type'][0]:
                    nprfo = 14
                for ref in reflData['RefList']:
                    if ref[nprfo+Super] < 0.:
                        hkls += ['%d %d %d: %.3f'%(ref[0],ref[1],ref[2],ref[nprfo+Super]),]
                    if len(hkls) > 10:
                        break
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
                        xmin=0.,nDig=(10,6),typeHint=float)
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
                from . import GSASIIdataGUI
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
        G2frame.bottomSizer,LeBailMsg = ShowHistogramInfo()
        mainSizer.Add(G2frame.bottomSizer)
        mainSizer.Layout()
        G2frame.dataWindow.Refresh()
        DData.SetVirtualSize(mainSizer.GetMinSize())
        DData.Scroll(0,Scroll)
        G2frame.dataWindow.SendSizeEvent()
        if LeBailMsg:
            G2G.G2MessageBox(G2frame,LeBailMsg,title='LeBail refinement changes')

    def ShowHistogramInfo():
        '''This creates a sizer with all the information pulled out from the Phase/data dict
        '''

        def OnUseData(event):
            Obj = event.GetEventObject()
            UseList[G2frame.hist]['Use'] = Obj.GetValue()
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

        def OnLeBail(event):
            '''Toggle the LeBail flag (Phase['Histograms'][hist]['LeBail'] and
            when turned on, set the Controls['newLeBail'] flag to indicate that
            the a new Le Bail extraction will be performed at the next
            refinement
            '''
            UseList[G2frame.hist]['LeBail'] = not UseList[G2frame.hist]['LeBail']
            if UseList[G2frame.hist]['LeBail']:
                Controls['newLeBail'] = True
            else:
                Controls['newLeBail'] = False
            wx.CallLater(100,RepaintHistogramInfo,DData.GetScrollPos(wx.VERTICAL))

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
        if 'Babinet' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Babinet'] = {'BabA':[0.0,False],'BabU':[0.0,False]}
        if 'Fix FXU' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Fix FXU'] = ' '
        if 'FixedSeqVars' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['FixedSeqVars'] = []
        if 'Flack' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Flack'] = [0.0,False]
        if 'Twins' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Twins'] = [[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False]],]
        if 'Layer Disp' not in UseList[G2frame.hist]:
            UseList[G2frame.hist]['Layer Disp'] = [0.0,False]
#end patch
        ifkeV = 'E' in UseList[G2frame.hist].get('Type','')
        offMsg = ''
        bottomSizer = wx.BoxSizer(wx.VERTICAL)
        useBox = wx.BoxSizer(wx.HORIZONTAL)
        useData = wx.CheckBox(DData,wx.ID_ANY,label='Use Histogram: '+G2frame.hist+' ?')
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(UseList[G2frame.hist]['Use'])
        useBox.Add(useData,0,WACV)
        if not generalData['doPawley'] and 'PWDR' in G2frame.hist[:4]:
            lbLabel = 'Start Le Bail extraction?   '
            if UseList[G2frame.hist]['LeBail']:
                lbLabel = 'Stop Le Bail extraction?'
                G2frame.SetStatusText('To reset Le Bail extracted intensities, cycle Le Bail button.',1)
            lebail = wx.Button(DData,wx.ID_ANY,label=lbLabel)
            lebail.Bind(wx.EVT_BUTTON, OnLeBail)
            useBox.Add(lebail,0,WACV)
        bottomSizer.Add(useBox)
        G2G.HorizontalLine(bottomSizer,DData)
        if G2frame.testSeqRefineMode() and not UseList[G2frame.hist]['LeBail']:
            bottomSizer.Add(wx.StaticText(DData,label='     Sequential Refinement Options'))
            parmChoice = [' ','X','XU','U','F','FX','FXU','FU']
            if generalData['Type'] == 'magnetic':
                parmChoice += ['M','MX','MXU','MU','MF','MFX','MFXU','MFU']
            fixBox = wx.BoxSizer(wx.HORIZONTAL)
            fixBox.Add(wx.StaticText(DData,label='     Fix these var types: '),0,WACV)
            fixVals = wx.ComboBox(DData,value=UseList[G2frame.hist]['Fix FXU'],choices=parmChoice,
                style=wx.CB_DROPDOWN)
            fixVals.Bind(wx.EVT_COMBOBOX,OnFixVals)
            fixBox.Add(fixVals,0,WACV)
            fixBox.Add(wx.StaticText(DData,label=' in phase '+generalData['Name']+' for this histogram'),0,WACV)
            bottomSizer.Add(fixBox)
            SeqId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, 'Sequential results')
            if SeqId:
                fixBox = wx.BoxSizer(wx.HORIZONTAL)
                fixBox.Add(wx.StaticText(DData,label='     Specific phase variables to fix for this histogram: '),0,WACV)
                addFixed = wx.Button(DData,wx.ID_ANY,label='Select Vars')
                fixBox.Add(addFixed,0,WACV)
                addFixed.Bind(wx.EVT_BUTTON,OnAddFixed)
                fixedVars = UseList[G2frame.hist].get('FixedSeqVars',[])
                if len(fixedVars):
                    fixBox.Add(wx.StaticText(DData,label=' (currently {} fixed)'.format(len(fixedVars))),0,WACV)
                bottomSizer.Add(fixBox)

        if not UseList[G2frame.hist]['LeBail'] or 'HKLF' in G2frame.hist[:4]:
            bottomSizer.Add(ScaleSizer(),0,wx.BOTTOM,5)
        else:
            # if phase fraction/scale is hidden, turn off flag
            UseList[G2frame.hist]['Scale'][1] = False

        if G2frame.hist[:4] == 'PWDR':
            if UseList[G2frame.hist]['Size'][0] == 'isotropic':
                isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                isoSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                isoSizer.Add(LGmixSizer('Size',[0.,1.],OnLGmixRef))
                isoSizer.Add(ResetSizer('isotropic',OnResetSize),0,WACV)
                bottomSizer.Add(isoSizer)
                bottomSizer.Add(IsoSizer(u'size(\xb5m): ','Size',(10,4),[0.,10.],OnSizeRef),0,wx.BOTTOM,5)
            elif UseList[G2frame.hist]['Size'][0] == 'uniaxial':
                uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                uniSizer.Add(TopSizer(' Domain size model: ',['isotropic','uniaxial','ellipsoidal'],
                    'Size',OnSizeType),0,WACV)
                uniSizer.Add(LGmixSizer('Size',[0.,1.],OnLGmixRef))
                uniSizer.Add(ResetSizer('uniaxial',OnResetSize),0,WACV)
                bottomSizer.Add(UniSizer('Size',OnSizeAxis),0)
                bottomSizer.Add(uniSizer)
                bottomSizer.Add(UniDataSizer(u'size(\xb5m): ','Size',(10,3),[0.,10.],OnSizeRef),0,wx.BOTTOM,5)
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
                bottomSizer.Add(GenStrainDataSizer(),0,wx.BOTTOM,5)

            bottomSizer.Add(wx.StaticText(DData,wx.ID_ANY,' Hydrostatic/elastic strain:'))
            bottomSizer.Add(HstrainSizer())
            bottomSizer.Add(DispSizer())

            if not UseList[G2frame.hist]['LeBail'] and not ifkeV:
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
                            ,0,wx.TOP|wx.BOTTOM,5)
                        poSizer.Add(SHDataSizer(POData))  #,0,wx.TOP|wx.BOTTOM,5)
                        try:
                            poSizer.Add(SHPenalty(POData))  #,0,wx.TOP|wx.BOTTOM,5)
                        except:
                            print('SHPenalty error occurred')
                bottomSizer.Add(poSizer)    #,0,wx.TOP|wx.BOTTOM,5)
                bottomSizer.Add(ExtSizer('PWDR'),0,wx.TOP|wx.BOTTOM,5)
                if generalData['Type'] != 'magnetic':
                    bottomSizer.Add(BabSizer(),0,wx.BOTTOM,5)
            else:
                # turn off preferred orientation fitting in LeBail mode (hidden)
                UseList[G2frame.hist]['Pref.Ori.'][2] = False
                for bab in ['A','U']:
                    UseList[G2frame.hist]['Babinet']['Bab'+bab][1] = False
                # TODO: should turn off all Extinction refinement flags
                UseList[G2frame.hist]['Extinction'][1] = False
        elif G2frame.hist[:4] == 'HKLF':
            bottomSizer.Add(ExtSizer('HKLF'))  #,0,wx.BOTTOM,5)
            bottomSizer.Add(BabSizer())  #,0,wx.BOTTOM,5)
            if not SGData['SGInv'] and len(UseList[G2frame.hist]['Twins']) < 2:
                bottomSizer.Add(FlackSizer())  #,0,wx.BOTTOM,5)
            bottomSizer.Add(twinSizer())  #,0,wx.BOTTOM,5)
        return bottomSizer,offMsg

    ######################################################################
    #### Beginning of UpdateDData execution here
    ######################################################################
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.GPXtree.root, 'Controls'))
    G2frame.SetStatusText('',1)
    keyList = G2frame.GetHistogramNames(['PWDR','HKLF'])
    UseList = data['Histograms']
    # look for histgrams that are no longer in the project (can happen when histogram deleted from tree)
    broken = [i for i in UseList if i not in keyList]
    if broken:
        msg = 'Removing histogram(s) referenced in this phase that are no longer in project:\n\n'
        for i,j in enumerate(broken):
            if i > 0: msg += ', '
            msg += j
            del data['Histograms'][j]
        G2G.G2MessageBox(G2frame,msg,'Dereferencing Removed histograms')
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
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl= f"Histogram data for Phase {data['General']['Name']!r}"[:60]
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    LeBailMsg = None
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
        G2frame.bottomSizer,LeBailMsg = ShowHistogramInfo()
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
    if LeBailMsg:
        G2G.G2MessageBox(G2frame,LeBailMsg,title='LeBail refinement changes')

def MakeHistPhaseWin(G2frame):
    '''Display Phase/Data (HAP) info from a Hist/Phase tree item. Used when
    the HAP info is shown as if a top-level entry in the data tree (see the
    SeparateHistPhaseTreeItem config option). This code is not used when
    the HAP info is shown in its original location, the Data tab in the
    Phase info.
    '''
    TabSelectionIdDict = {}
    def OnSelectPage(event):
        'Called when an item is selected from the Select page menu'
        tabname = TabSelectionIdDict.get(event.GetId()) # lookup phase
        if not tabname:
            print ('Warning: menu item not in dict! id= %d'%event.GetId())
            return
        # find the tab matching the phase
        for i,page in enumerate(phaseList):
            if tabname == phaseList[i]:
                HAPBook.SetSelection(i)
                wx.CallAfter(FillDDataWindow,i) # may result in a double paint on some OSs
                return
        else:
            print ("Warning: tab "+tabname+" was not found")

    def OnPageChanged(event):
        'respond to a notebook tab'
        page = event.GetSelection()
        wx.CallAfter(FillDDataWindow,page)

    def getDDataWindow():
        'Get the current scrollwindow for selected phase'
        return DData[HAPBook.GetSelection()]

    def getDDataPhaseinfo():
        'Get the data tree entry for the currently selected phase'
        page = HAPBook.GetSelection()
        return G2frame.GPXtree.GetItemPyData(phaseIds[page])

    def FillDDataWindow(page):
        'display the DData info'
        G2frame.HistPhaseLastSel = phaseList[page]
        data = G2frame.GPXtree.GetItemPyData(phaseIds[page])
        G2plt.PlotSizeStrainPO(G2frame,data,hist='')
        UpdateDData(G2frame,DData[page],data)

    #### G2frame.dataWindow.DataMenu/"Edit Phase" menu routines follow
       # where these duplicate similar routines in GSASIIphsGUI.py.
    def OnHklfAdd(event):
        '''Called to link a Single Xtal (HKLF) dataset to a selected phase.
        Most commonly, the histogram and phase are linked when the latter
        item is read in (routines OnImportPhase or OnImportSfact in
        func:`GSASIIdataGUI.GSASIImain`, but one can defer this or change
        the linking later using this routine.

        Note that this nearly identical to routine OnHklfAdd
        inside :func:`GSASIIphsGUI.UpdatePhaseData`.
        '''
        data = getDDataPhaseinfo()
        result = G2phG.CheckAddHKLF(G2frame,data)
        if result is None: return
        wx.CallAfter(UpdateDData,G2frame,getDDataWindow(),data)

    def OnDataUse(event):
        data = getDDataPhaseinfo()
#        hist = G2frame.hist
        if data['Histograms']:
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Use histograms',
                'Use which histograms?',G2frame.dataWindow.HistsInPhase)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelections()
                    for Id,item in enumerate(G2frame.dataWindow.HistsInPhase):
                        if Id in sel:
                            data['Histograms'][item]['Use'] = True
                        else:
                            data['Histograms'][item]['Use'] = False
            finally:
                dlg.Destroy()
        wx.CallAfter(UpdateDData,G2frame,getDDataWindow(),data)

    def OnDataCopy(event):
        data = getDDataPhaseinfo()
        hist = G2frame.hist
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        sourceDict = copy.deepcopy(data['Histograms'][hist])
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        copyDict = {}
        for name in copyNames:
            if name not in sourceDict: continue
            copyDict[name] = copy.deepcopy(sourceDict[name])        #force copy
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy phase/histogram parameters\nfrom '+hist[5:][:35],
                'Copy phase/hist parameters', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    data['Histograms'][keyList[sel]].update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnDataCopyFlags(event):
        data = getDDataPhaseinfo()
        hist = G2frame.hist
        sourceDict = copy.deepcopy(data['Histograms'][hist])
        copyDict = {}
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        babNames = ['BabA','BabU']
        for name in copyNames:
            if name not in sourceDict: continue
            if name in ['Scale','Extinction','HStrain','Flack','Twins','Layer Disp']:
                if name == 'Extinction' and 'HKLF' in sourceDict['Histogram']:
                    copyDict[name] = {name:[sourceDict[name][:2]]}
                    for item in ['Eg','Es','Ep']:
                        copyDict[name][item] = sourceDict[name][2][item][1]
                elif name == 'Twins':
                    copyDict[name] = sourceDict[name][0][1][1]
                else:
                    copyDict[name] = sourceDict[name][1]
            elif name in ['Size','Mustrain']:
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2],sourceDict[name][5]]
            elif name == 'Pref.Ori.':
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2]]
                if sourceDict[name][0] == 'SH':
                    SHterms = sourceDict[name][5]
                    SHflags = {}
                    for item in SHterms:
                        SHflags[item] = SHterms[item]
                    copyDict[name].append(SHflags)
            elif name == 'Babinet':
                copyDict[name] = {}
                for bab in babNames:
                    copyDict[name][bab] = sourceDict[name][bab][1]
            elif name in ['Fix FXU','FixedSeqVars']:
                copyDict[name] = copy.deepcopy(sourceDict[name])
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy phase/histogram flags\nfrom '+hist[5:][:35],
                'Copy phase/hist flags', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    item = keyList[sel]
                    for name in copyNames:
                        if name not in sourceDict: continue
                        if name in ['Scale','Extinction','HStrain','Flack','Twins','Layer Disp']:
                            if name == 'Extinction' and 'HKLF' in sourceDict['Histogram']:
                                data['Histograms'][item][name][:2] = copy.deepcopy(sourceDict[name][:2])
                                for itm in ['Eg','Es','Ep']:
                                    data['Histograms'][item][name][2][itm][1] = copy.deepcopy(copyDict[name][itm])
                            elif name == 'Twins':
                                data['Histograms'][item]['Twins'][0][1][1] = copy.deepcopy(copyDict['Twins'])
                            else:
                                try:
                                    data['Histograms'][item][name][1] = copy.deepcopy(copyDict[name])
                                except KeyError:
                                    continue
                        elif name in ['Size','Mustrain']:
                            data['Histograms'][item][name][0] = copy.deepcopy(copyDict[name][0])
                            data['Histograms'][item][name][2] = copy.deepcopy(copyDict[name][1])
                            data['Histograms'][item][name][5] = copy.deepcopy(copyDict[name][2])
                        elif name == 'Pref.Ori.':
                            data['Histograms'][item][name][0] = copy.deepcopy(copyDict[name][0])
                            data['Histograms'][item][name][2] = copy.deepcopy(copyDict[name][1])
                            if sourceDict[name][0] == 'SH':
                               SHflags = copy.deepcopy(copyDict[name][2])
                               SHterms = copy.deepcopy(sourceDict[name][5])
                               data['Histograms'][item][name][6] = copy.deepcopy(sourceDict[name][6])
                               data['Histograms'][item][name][7] = copy.deepcopy(sourceDict[name][7])
                        elif name == 'Babinet':
                            for bab in babNames:
                                data['Histograms'][item][name][bab][1] = copy.deepcopy(copyDict[name][bab])
                        elif name in ['Fix FXU','FixedSeqVars']:
                            data['Histograms'][item][name] = copy.deepcopy(sourceDict[name])
        finally:
            dlg.Destroy()

    def OnSelDataCopy(event):
        data = getDDataPhaseinfo()
        hist = G2frame.hist
        sourceDict = data['Histograms'][hist]
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select which parameters to copy',
            'Select phase data parameters', copyNames)
        selectedItems = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                selectedItems = [copyNames[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if not selectedItems: return # nothing to copy
        copyDict = {}
        for parm in selectedItems:
            if parm not in sourceDict: continue
            copyDict[parm] = copy.deepcopy(sourceDict[parm])
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy selected phase/histogram parameters\nfrom '+hist[5:][:35],
            'Copy selected phase/hist parameters', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    data['Histograms'][keyList[sel]].update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnPwdrAdd(event):
        data = getDDataPhaseinfo()
        generalData = data['General']
        SGData = generalData['SGData']
        newList = []
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        keyList = data['Histograms'].keys()
        TextList = []
        if G2frame.GPXtree.GetCount():
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.GPXtree.GetItemText(item)
                if name not in keyList and 'PWDR' in name:
                    TextList.append(name)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
            if not TextList:
                G2G.G2MessageBox(G2frame,'No histograms')
                return
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select powder histograms to use',
                'Use data',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: newList.append(TextList[i])
                    if 'All PWDR' in newList:
                        newList = TextList[1:]
                    for histoName in newList:
                        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
                        data['Histograms'][histoName] = {'Histogram':histoName,'Show':False,'LeBail':False,
                            'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},['',],0.1],
                            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                                [1.,1.,1.,0.,0.,0.],6*[False,]],
                            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                            'HStrain':[NDij*[0.0,],NDij*[False,]],
                            'Layer Disp':[0.0,False],
                            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},'Fix FXU':' ','FixedSeqVars':[]}
                        refList = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Reflection Lists'))
                        refList[generalData['Name']] = {}
                    wx.CallAfter(UpdateDData,G2frame,getDDataWindow(),data)
            finally:
                dlg.Destroy()

    def OnDataDelete(event):
        data = getDDataPhaseinfo()
        if G2frame.dataWindow.HistsInPhase:
            DelList = []
            extraOpts= {'label_0':'Remove from all phases','value_0':False}
            h,pd = G2frame.GetUsedHistogramsAndPhasesfromTree()
            if len(pd) > 1:
                opts = extraOpts
            else:
                opts = {}
            dlg = G2G.G2MultiChoiceDialog(G2frame,
                'Select histogram(s) to remove   \nfrom this phase:',
                'Remove histograms', G2frame.dataWindow.HistsInPhase,
                extraOpts=opts)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    DelList = [G2frame.dataWindow.HistsInPhase[i] for i in dlg.GetSelections()]
            finally:
                dlg.Destroy()
            if extraOpts['value_0']:
                for p in pd:
                    for i in DelList:
                        if i in pd[p]['Histograms']: del pd[p]['Histograms'][i]
            else:
                for i in DelList:
                    del data['Histograms'][i]
            wx.CallLater(100,UpdateDData,G2frame,getDDataWindow(),data)

    def OnDataApplyStrain(event):
        data = getDDataPhaseinfo()
        SGData = data['General']['SGData']
        DijVals = data['Histograms'][G2frame.hist]['HStrain'][0][:]
        # apply the Dij values to the reciprocal cell
        newA = []
        Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
        for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),
                            ['D11','D22','D33','D12','D13','D23']):
            newA.append(Aij + Dijdict.get(lbl,0.0))
        # convert back to direct cell
        data['General']['Cell'][1:7] = G2lat.A2cell(newA)
        data['General']['Cell'][7] = G2lat.calc_V(newA)
        # subtract the selected histograms Dij values from all for this phase
        for hist in data['Histograms']:
            for i,val in enumerate(DijVals):
                data['Histograms'][hist]['HStrain'][0][i] -= val
        # for hist in sorted(data['Histograms']): # list effective lattice constants applying Dij values
        #     DijVals = data['Histograms'][hist]['HStrain'][0]
        #     newA = []
        #     Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
        #     for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),
        #                     ['D11','D22','D33','D12','D13','D23']):
        #         newA.append(Aij + Dijdict.get(lbl,0.0))
        #     print(hist, G2lat.A2cell(newA)[:3], G2lat.calc_V(newA))
        wx.CallAfter(UpdateDData,G2frame,getDDataWindow(),data)

    #### start of MakeHistPhaseWin.
    # binds the menu items and shows the Data Window info
    G2frame.dataWindow.ClearData()
    HAPBook = G2G.GSNoteBook(parent=G2frame.dataWindow)
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(HAPBook,1,wx.ALL|wx.EXPAND,1)
    phaseList = []
    phaseIds = []
    DData = []
    sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
    item, cookie = G2frame.GPXtree.GetFirstChild(sub)
    while item: # loop over phases
        phaseName = G2frame.GPXtree.GetItemText(item)
        phaseIds.append(item)
        phaseList.append(phaseName)
        item, cookie = G2frame.GPXtree.GetNextChild(sub, cookie)
        HAPtab = wx.ScrolledWindow(HAPBook)
        HAPBook.AddPage(HAPtab,phaseName)
        DData.append(HAPtab)
    HAPBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    # set up "Select tab" menu contents
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataMenu)
    mid = G2frame.dataWindow.DataMenu.FindMenu('Select tab')
    menu = G2frame.dataWindow.DataMenu.GetMenu(mid)
    items = menu.GetMenuItems()
    for item in items:
         menu.Remove(item)
    if len(phaseList) == 0: return
    for i,page in enumerate(phaseList):
        Id = wx.NewId()
        if menu.FindItem(page) >= 0: continue # is tab already in menu?
        menu.Append(Id,page,'')
        TabSelectionIdDict[Id] = page
        G2frame.Bind(wx.EVT_MENU, OnSelectPage, id=Id)
    # define commands in G2frame.dataWindow.DataMenu/"Edit Phase" menu
    G2frame.Bind(wx.EVT_MENU, OnDataUse, id=G2G.wxID_DATAUSE)
    G2frame.Bind(wx.EVT_MENU, OnDataCopy, id=G2G.wxID_DATACOPY)
    G2frame.Bind(wx.EVT_MENU, OnDataCopyFlags, id=G2G.wxID_DATACOPYFLAGS)
    G2frame.Bind(wx.EVT_MENU, OnSelDataCopy, id=G2G.wxID_DATASELCOPY)
    G2frame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2G.wxID_PWDRADD)
    G2frame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2G.wxID_HKLFADD)
    G2frame.Bind(wx.EVT_MENU, OnDataDelete, id=G2G.wxID_DATADELETE)
    G2frame.Bind(wx.EVT_MENU, OnDataApplyStrain, id=G2G.wxID_DATADIJ)
    # display the last-selected phase or the 1st
    try:
        G2frame.HistPhaseLastSel
    except:
        G2frame.HistPhaseLastSel = phaseList[0]
    if G2frame.HistPhaseLastSel in phaseList:
        page = phaseList.index(G2frame.HistPhaseLastSel)
    else:
        page = 0
    HAPBook.SetSelection(page)
    wx.CallAfter(FillDDataWindow,page)
