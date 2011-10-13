#GSASII - data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import wx
import wx.grid as wg
import numpy as np
import math
import time
import cPickle
import GSASIIpath
import GSASIIpwd as G2pwd
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIindex as G2indx
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIElem as G2elem

VERY_LIGHT_GREY = wx.Colour(235,235,235)

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
    
def IsHistogramInAnyPhase(self,histoName):
    phases = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
    if phases:
        item, cookie = self.PatternTree.GetFirstChild(phases)
        while item:
            data = self.PatternTree.GetItemPyData(item)
            histoList = data['Histograms'].keys()
            if histoName in histoList:
                return True
            item, cookie = self.PatternTree.GetNextChild(phases, cookie)
        return False
    else:
        return False
    
       
def UpdatePeakGrid(self, data):
    if self.dataDisplay:
        self.dataFrame.Clear()
    
    def OnUnDo(event):
        DoUnDo()
        self.dataFrame.UnDo.Enable(False)
        
    def DoUnDo():
        print 'Undo last refinement'
        file = open(self.undofile,'rb')
        PatternId = self.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item),cPickle.load(file))
            if self.dataDisplay.GetName() == item:
                if item == 'Background':
                    UpdateBackgroundGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
                elif item == 'Instrument Parameters':
                    UpdateInstrumentGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
                elif item == 'Peak List':
                    UpdatePeakGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
            print item,' recovered'
        file.close()
        
    def SaveState():
        self.undofile = self.dirname+'\\GSASII.save'
        file = open(self.undofile,'wb')
        PatternId = self.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            cPickle.dump(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId,item)),file,1)
        file.close()
        self.dataFrame.UnDo.Enable(True)
        
    def OnLSQPeakFit(event):
        OnPeakFit('LSQ')
        
    def OnOneCycle(event):
        OnPeakFit('LSQ',oneCycle=True)
        
    def OnClearPeaks(event):
        dlg = wx.MessageDialog(self,'Delete all peaks?','Clear peak list',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = []
        finally:
            dlg.Destroy()
        UpdatePeakGrid(self,peaks)
        G2plt.PlotPatterns(self)
        
    def OnPeakFit(FitPgm,oneCycle=False):
        SaveState()
        controls = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.root, 'Controls'))
        if not controls:
            controls = {'deriv type':'analytic','min dM/M':0.0001,}     #fill in defaults if needed
        print 'Peak Fitting with '+controls['deriv type']+' derivatives:'
        PatternId = self.PatternId
        PickId = self.PickId
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
        if not peaks:
            self.ErrorDialog('No peaks!','Nothing to fit!')
            return
        background = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Background'))[0]
        limits = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Limits'))[1]
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
        data = self.PatternTree.GetItemPyData(PatternId)[1]
        wx.BeginBusyCursor()
        try:
            G2pwd.DoPeakFit(FitPgm,peaks,background,limits,inst,data,oneCycle,controls)
        finally:
            wx.EndBusyCursor()    
        UpdatePeakGrid(self,peaks)
        G2plt.PlotPatterns(self)
        print 'finished'
        return
        
    def OnResetSigGam(event):
        PatternId = self.PatternId
        PickId = self.PickId
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
        if not peaks:
            self.ErrorDialog('No peaks!','Nothing to do!')
            return
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
        Inst = dict(zip(inst[3],inst[1]))
        print len(Inst['Type'])
        for peak in peaks:
            if Inst['Type'] in ['PXC','PNC']:
                peak[4] = Inst['U']*tand(peak[0]/2.0)**2+Inst['V']*tand(peak[0]/2.0)+Inst['W']
                peak[6] = Inst['X']/cosd(peak[0]/2.0)+Inst['Y']*tand(peak[0]/2.0)
        UpdatePeakGrid(self,peaks)
                
    def RefreshPeakGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        
        event.StopPropagation()
        data = self.PeakTable.GetData()
        T = []
        for peak in data:T.append(peak[0])
        D = dict(zip(T,data))
        T.sort()
        X = []
        for key in T: X.append(D[key])
        data = X        
        
    def setBackgroundColors():
       for r in range(self.dataDisplay.GetNumberRows()):
           for c in range(self.dataDisplay.GetNumberCols()):
               if self.dataDisplay.GetColLabelValue(c) in ['position','intensity','sigma','gamma']:
                   if float(self.dataDisplay.GetCellValue(r,c)) < 0.:
                       self.dataDisplay.SetCellBackgroundColour(r,c,wx.RED)
                   else:
                       self.dataDisplay.SetCellBackgroundColour(r,c,wx.WHITE)
                                                  
    def KeyEditPeakGrid(event):
        rowList = self.dataDisplay.GetSelectedRows()
        colList = self.dataDisplay.GetSelectedCols()
        selectList = self.dataDisplay.GetSelectedCells()
        data = self.PatternTree.GetItemPyData(self.PickId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif rowList:
            self.dataDisplay.ClearSelection()
            if event.GetKeyCode() == wx.WXK_DELETE:
                self.dataDisplay.ClearGrid()
                rowList.reverse()
                nDel = 0
                for row in rowList:
                    self.PeakTable.DeleteRow(row)
                    nDel += 1
                if nDel:
                    msg = wg.GridTableMessage(self.PeakTable, 
                        wg.GRIDTABLE_NOTIFY_ROWS_DELETED,0,nDel)
                    self.dataDisplay.ProcessTableMessage(msg)
                data = self.PeakTable.GetData()
                self.PatternTree.SetItemPyData(self.PickId,data[:-nDel])
                self.dataDisplay.ForceRefresh()
                setBackgroundColors()
                if not len(self.PatternTree.GetItemPyData(self.PickId)): 
                    self.dataFrame.PeakFit.Enable(False)
                        
        elif colList:
            self.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if self.PeakTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(self.PeakTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(self.PeakTable.GetNumberRows()): data[row][col]=False
        elif selectList:
            self.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for row,col in selectList:
                if self.PeakTable.GetTypeName(row,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        data[row][col]=True
                    elif key == 78:  #'N'
                        data[row][col]=False
        G2plt.PlotPatterns(self)
            
    self.dataFrame.SetMenuBar(self.dataFrame.PeakMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    self.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_UNDO)
    self.Bind(wx.EVT_MENU, OnLSQPeakFit, id=G2gd.wxID_LSQPEAKFIT)
    self.Bind(wx.EVT_MENU, OnOneCycle, id=G2gd.wxID_LSQONECYCLE)
    self.Bind(wx.EVT_MENU, OnClearPeaks, id=G2gd.wxID_CLEARPEAKS)
    self.Bind(wx.EVT_MENU, OnResetSigGam, id=G2gd.wxID_RESETSIGGAM)
    self.dataFrame.PeakFit.Enable(False)
    if data:
        self.dataFrame.PeakFit.Enable(True)
        self.dataFrame.PFOneCycle.Enable(True)
    self.PickTable = []
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','refine','intensity','refine','sigma','refine','gamma','refine']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
    T = []
    for peak in data:
        T.append(peak[0])
    D = dict(zip(T,data))
    T.sort()
    X = []
    for key in T: X.append(D[key])
    data = X        
    self.PatternTree.SetItemPyData(self.PickId,data)
    self.PeakTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Peak List')
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.PeakTable, True)
    setBackgroundColors()                         
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshPeakGrid)
    self.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([535,350])
        
def UpdateBackgroundGrid(self,data):
    ValObj = {}
    
    def OnNewType(event):
        data[0][0] = bakType.GetValue()
        
    def OnBakRef(event):
        data[0][1] = bakRef.GetValue()
        
    def OnBakTerms(event):
        data[0][2] = int(bakTerms.GetValue())
        M = len(data[0])
        N = data[0][2]+3
        item = data[0]
        if N > M:       #add terms
            for i in range(M,N): 
                item.append(0.0)
        elif N < M:     #delete terms
            for i in range(N,M):
                del(item[-1])
        self.PatternTree.SetItemPyData(BackId,data)
        UpdateBackgroundGrid(self,data)
        
    def OnBakVal(event):
        Obj = event.GetEventObject()
        item = ValObj[Obj.GetId()][0]
        try:
            value = float(Obj.GetValue())
        except ValueError:
            value = data[0][item]
        data[0][item] = value
        Obj.SetValue('%10.4f'%(value))
        
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataDisplay = wx.Panel(self.dataFrame)
    BackId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Background')
    Choices = ['chebyschev','cosine','lin interpolate','inv interpolate','log interpolate']
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    topSizer.Add(wx.StaticText(self.dataDisplay,-1,' Background function: '),0,wx.ALIGN_CENTER_VERTICAL)
    bakType = wx.ComboBox(self.dataDisplay,value=data[0][0],
            choices=Choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
    bakType.Bind(wx.EVT_COMBOBOX, OnNewType)
    topSizer.Add(bakType)
    topSizer.Add((5,0),0)
    bakRef = wx.CheckBox(self.dataDisplay,label=' Refine?')
    bakRef.SetValue(bool(data[0][1]))
    bakRef.Bind(wx.EVT_CHECKBOX, OnBakRef)
    topSizer.Add(bakRef,0,wx.ALIGN_CENTER_VERTICAL)
    topSizer.Add(wx.StaticText(self.dataDisplay,-1,' No. coeff.: '),0,wx.ALIGN_CENTER_VERTICAL)
    bakTerms = wx.ComboBox(self.dataDisplay,-1,value=str(data[0][2]),choices=[str(i+1) for i in range(36)],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    bakTerms.Bind(wx.EVT_COMBOBOX,OnBakTerms)
    topSizer.Add(bakTerms,0,wx.ALIGN_CENTER_VERTICAL)
    topSizer.Add((5,0),0)
    mainSizer.Add(topSizer)
    mainSizer.Add(wx.StaticText(self.dataDisplay,-1,' Background coefficients:'),0,wx.ALIGN_CENTER_VERTICAL)
    bakSizer = wx.FlexGridSizer(1,5,5,5)
    for i,value in enumerate(data[0][3:]):
        bakVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,'%10.4f'%(value),style=wx.TE_PROCESS_ENTER)
        bakSizer.Add(bakVal,0,wx.ALIGN_CENTER_VERTICAL)
        ValObj[bakVal.GetId()] = [i+3]
        bakVal.Bind(wx.EVT_TEXT_ENTER,OnBakVal)
        bakVal.Bind(wx.EVT_KILL_FOCUS,OnBakVal)        
    mainSizer.Add(bakSizer)
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
        
def UpdateLimitsGrid(self, data):
    if self.dataDisplay:
        self.dataFrame.Clear()
        
    def RefreshLimitsGrid(event):
        event.StopPropagation()
        data = self.LimitsTable.GetData()
        old = data[0]
        new = data[1]
        new[0] = max(old[0],new[0])
        new[1] = max(new[0],min(old[1],new[1]))
        data = [old,new]
        G2plt.PlotPatterns(self)
        
    self.LimitsTable = []
    colLabels = ['Tmin','Tmax']
    rowLabels = ['original','changed']
    Types = 2*[wg.GRID_VALUE_FLOAT+':10,3',]
    self.LimitsTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Limits')
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.LimitsTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshLimitsGrid)                
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([230,120])
    
def UpdateInstrumentGrid(self,data):
    if len(data) > 3:                   #powder data
        insVal = dict(zip(data[3],data[1]))
        insDef = dict(zip(data[3],data[0]))
        insRef = dict(zip(data[3],data[2]))
        if 'N' in insDef['Type']:
            del(insDef['Polariz.'])
            del(insVal['Polariz.'])
            del(insRef['Polariz.'])
    else:                               #single crystal data
        insVal = dict(zip(data[2],data[1]))
        insDef = dict(zip(data[2],data[0]))
        insRef = {}
    ValObj = {}
    RefObj = {}
    waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
        'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'MoKa':[0.70926,0.713543],
        'AgKa':[0.559363,0.563775]}
        
    def inst2data(inst,ref,data):
        if len(data) > 3:
            for i,item in enumerate(data[3]):
                try:
                    data[1][i] = inst[item]
                    data[2][i] = ref[item]
                except KeyError:
                    data[1][i] = 0
                    data[2][i] = 0                    
        else:
            for i,item in enumerate(data[2]):
                data[1][i] = inst[item]            
        return data
        
    def updateData(inst,ref):
        return inst2data(inst,ref,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,
            self.PatternId,'Instrument Parameters')))        
    
    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List'))
            if 'P' in insVal['Type']:                                       #update powder peak parameters
                for peak in peaks:
                    peak[4] = insVal['U']*tand(peak[0]/2.0)**2+insVal['V']*tand(peak[0]/2.0)+insVal['W']
                    peak[6] = insVal['X']/cosd(peak[0]/2.0)+insVal['Y']*tand(peak[0]/2.0)
                                                
    def OnReset(event):
        insVal.update(insDef)
        data = updateData(insVal,insRef)
        RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
        UpdateInstrumentGrid(self,data)
        
    def OnWaveChange(event):
        if 'Lam' in insVal:            
            data[0] = data[0][:1]+tuple(waves['CuKa'])+(.5,)+data[0][2:]
            data[1] = data[1][:1]+waves['CuKa']+[.5,]+data[1][2:]
            data[2] = data[2][:1]+[0,0,0,]+data[2][2:]
            data[3] = data[3][:1]+['Lam1','Lam2','I(L2)/I(L1)',]+data[3][2:]            
        else:
            data[0] = data[0][:2]+data[0][4:]
            data[1] = data[1][:2]+data[1][4:]
            data[2] = data[2][:2]+data[2][4:]
            data[3] = data[3][:1]+['Lam',]+data[3][4:]            
        UpdateInstrumentGrid(self,data)
                
    def OnNewType(event):
        insVal['Type'] = typePick.GetValue()
        data = updateData(insVal,insRef)
        UpdateInstrumentGrid(self,data)
        
    def OnLamPick(event):
        lamType = lamPick.GetValue()
        insVal['Lam1'] = waves[lamType][0]
        insVal['Lam2'] = waves[lamType][1]
        data = updateData(insVal,insRef)
        UpdateInstrumentGrid(self,data)
                 
    def OnRatValue(event):
        try:
            value = float(ratVal.GetValue())
            if value < 0:
                raise ValueError
        except ValueError:
            value = insVal['I(L2)/I(L1)']
        insVal['I(L2)/I(L1)'] = value
        ratVal.SetValue('%10.4f'%(value))
        data = updateData(insVal,insRef)
        
    def OnRatRef(event):
        insRef['I(L2)/I(L1)'] = ratRef.GetValue()
        data = updateData(insVal,insRef)
        
    def OnWaveValue(event):
        try:
            value = float(waveVal.GetValue())
            if value < 0:
                raise ValueError
        except ValueError:
            value = insVal['Lam']
        insVal['Lam'] = value
        waveVal.SetValue('%10.6f'%(value))
        data = updateData(insVal,insRef)
        
    def OnWaveRef(event):
        insRef['Lam'] = waveRef.GetValue()
        data = updateData(insVal,insRef)
        
    def OnItemValue(event):
        Obj = event.GetEventObject()
        item,fmt = ValObj[Obj.GetId()]
        try:
            value = float(Obj.GetValue())
        except ValueError:
            value = insVal[item]
        insVal[item] = value
        Obj.SetValue(fmt%(value))
        data = updateData(insVal,insRef)
        
    def OnItemRef(event):
        Obj = event.GetEventObject()
        item = RefObj[Obj.GetId()]
        insRef[item] = Obj.GetValue()
        data = updateData(insVal,insRef)
                
    if self.dataDisplay:
        self.dataFrame.Clear()
    histoName = self.PatternTree.GetItemPyData(self.PatternId)[-1]
    ifHisto = IsHistogramInAnyPhase(self,histoName)
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataDisplay = wx.Panel(self.dataFrame)
    instSizer = wx.FlexGridSizer(2,6,5,5)
    instSizer.Add(wx.StaticText(self.dataDisplay,-1,' Histogram Type:'),0,wx.ALIGN_CENTER_VERTICAL)
    if 'P' in insVal['Type']:                   #powder data
        self.dataFrame.SetMenuBar(self.dataFrame.InstMenu)
        if not self.dataFrame.GetStatusBar():
            Status = self.dataFrame.CreateStatusBar()
        self.Bind(wx.EVT_MENU, OnReset,id=G2gd.wxID_INSTPRMRESET)
        self.Bind(wx.EVT_MENU,OnWaveChange,id=G2gd.wxID_CHANGEWAVETYPE)
        typePick = wx.ComboBox(self.dataDisplay,value=insVal['Type'],
            choices=['PXC','PNC','PNT'],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        typePick.Bind(wx.EVT_COMBOBOX, OnNewType)
        instSizer.Add(typePick,0,wx.ALIGN_CENTER_VERTICAL)
        if 'C' in insVal['Type']:               #constant wavelength
            #patch
            if 'Azimuth' not in insVal:
                insVal['Azimuth'] = 0.0
                insDef['Azimuth'] = 0.0
                insRef['Azimuth'] = False
            #end of patch
            instSizer.Add(wx.StaticText(self.dataDisplay,-1,' Azimuth: %7.2f'%(insVal['Azimuth'])),0,wx.ALIGN_CENTER_VERTICAL)
            if 'Lam1' in insVal:
                instSizer.Add((5,5),0)
                instSizer.Add((5,5),0)
                instSizer.Add((5,5),0)
                instSizer.Add(wx.StaticText(self.dataDisplay,-1,' Ka1/Ka2:'),
                        0,wx.ALIGN_CENTER_VERTICAL)
                instSizer.Add(wx.StaticText(self.dataDisplay,-1,'%8.6f/%8.6f'%(insVal['Lam1'],insVal['Lam2'])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(self.dataDisplay,-1,'Select:'),0,wx.ALIGN_CENTER_VERTICAL)
                choice = ['TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa']
                lamPick = wx.ComboBox(self.dataDisplay,value=' ',choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0)
                instSizer.Add(waveSizer,0)
                instSizer.Add(wx.StaticText(self.dataDisplay,-1,' I(L2)/I(L1): (%10.4f)'%(insDef['I(L2)/I(L1)'])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                ratVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,'%10.4f'%(insVal['I(L2)/I(L1)']),style=wx.TE_PROCESS_ENTER)
                ratVal.Bind(wx.EVT_TEXT_ENTER,OnRatValue)
                ratVal.Bind(wx.EVT_KILL_FOCUS,OnRatValue)
                instSizer.Add(ratVal,0)
                ratRef = wx.CheckBox(self.dataDisplay,label=' Refine?')
                ratRef.SetValue(bool(insRef['I(L2)/I(L1)']))
                ratRef.Bind(wx.EVT_CHECKBOX, OnRatRef)
                instSizer.Add(ratRef,0,wx.ALIGN_CENTER_VERTICAL)
                
            else:
                instSizer.Add(wx.StaticText(self.dataDisplay,-1,' Lam: (%10.6f)'%(insDef['Lam'])),
                    0,wx.ALIGN_CENTER_VERTICAL)
                waveVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,'%10.6f'%(insVal['Lam']),style=wx.TE_PROCESS_ENTER)
                waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveValue)
                waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveValue)
                instSizer.Add(waveVal,0,wx.ALIGN_CENTER_VERTICAL)
                if ifHisto:
                    waveRef = wx.CheckBox(self.dataDisplay,label=' Refine?')
                    waveRef.SetValue(bool(insRef['Lam']))
                    waveRef.Bind(wx.EVT_CHECKBOX, OnWaveRef)
                    instSizer.Add(waveRef,0,wx.ALIGN_CENTER_VERTICAL)
                else:
                    instSizer.Add((5,5),0)
            for item in ['Zero','Polariz.']:
                fmt = '%10.3f'
                Fmt = ' %s: ('+fmt+')'
                if item in insDef:
                    instSizer.Add(wx.StaticText(self.dataDisplay,-1,Fmt%(item,insDef[item])),
                            0,wx.ALIGN_CENTER_VERTICAL)
                    itemVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                    ValObj[itemVal.GetId()] = [item,fmt]
                    itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                    itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                    instSizer.Add(itemVal,0,wx.ALIGN_CENTER_VERTICAL)
                    if ifHisto:
                        itemRef = wx.CheckBox(self.dataDisplay,wx.ID_ANY,label=' Refine?')
                        itemRef.SetValue(bool(insRef[item]))
                        RefObj[itemRef.GetId()] = item
                        itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                        instSizer.Add(itemRef,0,wx.ALIGN_CENTER_VERTICAL)
                    else:
                        instSizer.Add((5,5),0)
                else:                           #skip Polariz. for neutrons
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)
            for item in ['U','V','W','X','Y','SH/L']:
                fmt = '%10.3f'
                if item == 'SH/L':
                    fmt = '%10.5f'
                Fmt = ' %s: ('+fmt+')'
                instSizer.Add(wx.StaticText(self.dataDisplay,-1,Fmt%(item,insDef[item])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                itemVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                ValObj[itemVal.GetId()] = [item,fmt]
                itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                instSizer.Add(itemVal,0,wx.ALIGN_CENTER_VERTICAL)
                itemRef = wx.CheckBox(self.dataDisplay,wx.ID_ANY,label=' Refine?')
                itemRef.SetValue(bool(insRef[item]))
                RefObj[itemRef.GetId()] = item
                itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                instSizer.Add(itemRef,0,wx.ALIGN_CENTER_VERTICAL)
        else:                                   #time of flight (neutrons)
            pass                                #for now
        
        

    else:                       #single crystal data
        typePick = wx.ComboBox(self.dataDisplay,value=insVal['Type'],
            choices=['SXC','SNC','SNT'],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        typePick.Bind(wx.EVT_COMBOBOX, OnNewType)
        instSizer.Add(typePick,0,wx.ALIGN_CENTER_VERTICAL)
        if 'C' in insVal['Type']:               #constant wavelength
            instSizer.Add(wx.StaticText(self.dataDisplay,-1,' Lam: %10.6f'%(insDef['Lam'])),
                    0,wx.ALIGN_CENTER_VERTICAL)
        else:                                   #time of flight (neutrons)
            pass                                #for now
        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(instSizer,0)
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
    
def UpdateSampleGrid(self,data):
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()    
    self.dataDisplay = wx.Panel(self.dataFrame)

#patch
    if not 'Gonio. radius' in data:
        data['Gonio. radius'] = 200.0
#patch end
    
    parms = [['Gonio. radius',' Goniometer radius(mm): ','%.2f',]]
    if data['Type'] == 'Debye-Scherrer':
        parms += [['DisplaceX',' Sample X displacement(\xb5m): ','%.2f',],
            ['DisplaceY',' Sample Y displacement(\xb5m): ','%.2f',],
            ['Absorption',' Sample absorption(\xb5r): ','%.4f',],]
    elif data['Type'] == 'Bragg-Brentano':
        parms += [['Shift',' Sample displacement(\xb5m): ','%.2f',],
            ['Transparency',' Sample transparency(1/\xb5eff,cm): ','%.4f'],]
    parms.append(['Temperature',' Sample temperature(K): ','%.2f'])
    parms.append(['Pressure',' Sample pressure(MPa): ','%.3f'])
    parms.append(['Humidity',' Sample humidity(%): ','%.1f'])
    parms.append(['Voltage',' Sample voltage(V): ','%.3f'])
    parms.append(['Force',' Applied load(MN): ','%.3f'])
    objList = {}

    def OnScaleRef(event):
        Obj = event.GetEventObject()
        data['Scale'][1] = Obj.GetValue()
        
    def OnScaleVal(event):
        Obj = event.GetEventObject()
        try:
            scale = float(Obj.GetValue())
            if scale > 0:
                data['Scale'][0] = scale
        except ValueError:
            pass
        Obj.SetValue("%.4f"%(data['Scale'][0]))          #reset in case of error
        
    def OnHistoType(event):
        Obj = event.GetEventObject()
        data['Type'] = Obj.GetValue()
        if data['Type'] == 'Bragg-Brentano' and 'Shift' not in data:    #set up defaults for new type(s)
            data['Shift'] = [0.0,False]
            data['Transparency'] = [0.0,False]
        self.dataDisplay.Destroy()
        UpdateSampleGrid(self,data)
        
    def OnParmRef(event):
        Obj = event.GetEventObject()
        parm = objList[Obj.GetId()]
        data[parm][1] = Obj.GetValue()
        
    def OnParmVal(event):
        Obj = event.GetEventObject()
        parm = objList[Obj.GetId()]
        try:
            if 'list' in str(type(data[parm[0]])): 
                data[parm[0]][0] = float(Obj.GetValue())
            else:
                data[parm[0]] = float(Obj.GetValue())
        except ValueError:
            pass
        if 'list' in str(type(data[parm[0]])): 
            Obj.SetValue(parm[2]%(data[parm[0]][0]))          #reset in case of error
        else:
            Obj.SetValue(parm[2]%(data[parm[0]]))          #reset in case of error
                
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Sample parameters: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    parmSizer = wx.FlexGridSizer(9,2,5,0)
    scaleRef = wx.CheckBox(self.dataDisplay,label=' Histogram scale factor: ')
    scaleRef.SetValue(data['Scale'][1])
    scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
    parmSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
    scaleVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,
        '%.4f'%(data['Scale'][0]),style=wx.TE_PROCESS_ENTER)
    scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
    scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
    parmSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
    typeSizer = wx.BoxSizer(wx.HORIZONTAL)
    choices = ['Debye-Scherrer','Bragg-Brentano',]
    histoType = wx.ComboBox(self.dataDisplay,wx.ID_ANY,value=data['Type'],choices=choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    histoType.Bind(wx.EVT_COMBOBOX, OnHistoType)
    parmSizer.Add(histoType)
    parmSizer.Add((5,5),0)
    
    for parm in parms:
        if 'list' in str(type(data[parm[0]])):
            parmRef = wx.CheckBox(self.dataDisplay,label=parm[1])
            objList[parmRef.GetId()] = parm[0]
            parmRef.SetValue(data[parm[0]][1])
            parmRef.Bind(wx.EVT_CHECKBOX, OnParmRef)
            parmSizer.Add(parmRef,0,wx.ALIGN_CENTER_VERTICAL)
            parmVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,
                parm[2]%(data[parm[0]][0]),style=wx.TE_PROCESS_ENTER)
        else:
            parmSizer.Add(wx.StaticText(self.dataDisplay,label=parm[1]),
                0,wx.ALIGN_CENTER_VERTICAL)
            parmVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,
                parm[2]%(data[parm[0]]),style=wx.TE_PROCESS_ENTER)        
        objList[parmVal.GetId()] = parm
        parmVal.Bind(wx.EVT_TEXT_ENTER,OnParmVal)
        parmVal.Bind(wx.EVT_KILL_FOCUS,OnParmVal)
        parmSizer.Add(parmVal,1,wx.EXPAND)
    mainSizer.Add(parmSizer)
    mainSizer.Add((0,5),0)    
    
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(self.dataFrame)
    self.dataDisplay.SetSize(Size)
    self.dataFrame.setSizePosLeft(Size)
                
def UpdateIndexPeaksGrid(self, data):
    IndexId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Index Peak List')
    
    def RefreshIndexPeaksGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        data = self.IndexPeaksTable.GetData()
        if c == 2:
            if data[r][c]:
                data[r][c] = False
            else:
                data[r][c] = True
            self.IndexPeaksTable.SetData(data)
            self.PatternTree.SetItemPyData(IndexId,data)
            self.dataDisplay.ForceRefresh()
            
    def OnReload(event):
        data = []
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List'))
        for peak in peaks:
            dsp = inst[1]/(2.0*sind(peak[0]/2.0))
            data.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])
        self.PatternTree.SetItemPyData(IndexId,data)
        UpdateIndexPeaksGrid(self,data)
        
    def KeyEditPickGrid(event):
        colList = self.dataDisplay.GetSelectedCols()
        rowList = self.dataDisplay.GetSelectedRows()
        data = self.PatternTree.GetItemPyData(IndexId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif colList:
            self.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if self.IndexPeaksTable.GetColLabelValue(col) in ['use','refine']:
                    if key == 89: #'Y'
                        for row in range(self.IndexPeaksTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(self.IndexPeaksTable.GetNumberRows()): data[row][col]=False
            
    if self.dataDisplay:
        self.dataFrame.Clear()
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    if 'PWD' in self.PatternTree.GetItemText(self.PatternId):
        self.dataFrame.SetMenuBar(self.dataFrame.IndPeaksMenu)
        self.Bind(wx.EVT_MENU, OnReload, id=G2gd.wxID_INDXRELOAD)
    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
    self.dataFrame.IndexPeaks.Enable(False)
    self.IndexPeaksTable = []
    if data:
        self.dataFrame.IndexPeaks.Enable(True)
        cells = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Unit Cells List'))
        if cells:
            cellist = cells[2]
            dmin = cells[3]
            self.HKL = []
            for i,cell in enumerate(cellist):
                if cell[-1]:
                    ibrav = cell[2]
                    A = G2lat.cell2A(cell[3:9])
                    self.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                    G2indx.IndexPeaks(data,self.HKL)
                    for hkl in self.HKL:
                        hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','intensity','use','indexed','h','k','l','d-obs','d-calc']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_BOOL,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5']
    self.PatternTree.SetItemPyData(IndexId,data)
    self.IndexPeaksTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Index Peak List')
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)                
    self.dataDisplay.SetTable(self.IndexPeaksTable, True)
    for r in range(self.dataDisplay.GetNumberRows()):
        for c in range(self.dataDisplay.GetNumberCols()):
            if c == 2:
                self.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                self.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, RefreshIndexPeaksGrid)
    self.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)                 
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([490,300])
  
def UpdateUnitCellsGrid(self, data):
    UnitCellsId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Unit Cells List')
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    spaceGroups = ['F m 3 m','I m 3 m','P m 3 m','R -3 H','P 6/m m m','I 4/m m m',
        'P 4/m m m','F m m m','I m m m','C m m m','P m m m','C 2/m','P 2/m','P -1']
        
    def SetLattice(controls):
        ibrav = bravaisSymb.index(controls[5])
        if ibrav in [0,1,2]:
            controls[7] = controls[8] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
        elif ibrav in [3,4,5,6]:
            controls[7] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
            if ibrav in [3,4]:
                controls[11] = 120.
        elif ibrav in [7,8,9,10]:
            controls[9] = controls[10] = controls[11] = 90.
        elif ibrav in [11,12]:
            controls[9] = controls[11] = 90.  # b unique
        if len(controls) < 13: controls.append(0)
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        return ibrav
        
    def OnNcNo(event):
        controls[2] = NcNo.GetValue()
        
    def OnStartVol(event):
        try:
            stVol = int(float(startVol.GetValue()))
            if stVol < 25:
                raise ValueError
        except ValueError:
            stVol = 25
        controls[3] = stVol
        startVol.SetValue("%d"%(stVol))
        
    def OnBravais(event):
        Obj = event.GetEventObject()
        bravais[bravList.index(Obj.GetId())] = Obj.GetValue()
        
    def OnZero(event):
        try:
            Zero = min(0.1,max(-0.1,float(zero.GetValue())))
        except ValueError:
            Zero = 0.0
        controls[1] = Zero
        zero.SetValue("%.4f"%(Zero))
        
    def OnZeroVar(event):
        controls[0] = zeroVar.GetValue()
        
    def OnBravSel(event):
        controls[5] = bravSel.GetString(bravSel.GetSelection())       
        UpdateUnitCellsGrid(self,data)
        
    def OnCellChange(event):
        ibrav = bravaisSymb.index(controls[5])
        Obj = event.GetEventObject()
        ObjId = cellList.index(Obj.GetId())
        try:
            value = max(1.0,float(Obj.GetValue()))
        except ValueError:
            if ObjId < 3:               #bad cell edge - reset
                value = controls[6+ObjId]
            else:                       #bad angle
                value = 90.
        if ibrav in [0,1,2]:
            controls[6] = controls[7] = controls[8] = value
            controls[9] = controls[10] = controls[11] = 90.0
            Obj.SetValue("%.5f"%(controls[6]))
        elif ibrav in [3,4,5,6]:
            if ObjId == 0:
                controls[6] = controls[7] = value
                Obj.SetValue("%.5f"%(controls[6]))
            else:
                controls[8] = value
                Obj.SetValue("%.5f"%(controls[8]))
            controls[9] = controls[10] = controls[11] = 90.0
            if ibrav in [3,4]:
                controls[11] = 120.
        elif ibrav in [7,8,9,10]:
            controls[6+ObjId] = value
            Obj.SetValue("%.5f"%(controls[6+ObjId]))
            controls[9] = controls[10] = controls[11] = 90.0
        elif ibrav in [11,12]:
            controls[9] = controls[11] = 90.0
            if ObjId != 3:
                controls[6+ObjId] = value
                Obj.SetValue("%.5f"%(controls[6+ObjId]))
            else:
                controls[10] = value
                Obj.SetValue("%.3f"%(controls[10]))
        else:
            controls[6+ObjId] = value
            if ObjId < 3:
                Obj.SetValue("%.5f"%(controls[6+ObjId]))
            else:
                Obj.SetValue("%.3f"%(controls[6+ObjId]))
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        volVal.SetValue("%.3f"%(controls[12]))
        
    def OnHklShow(event):
        hklShow.SetValue(False)
        PatternId = self.PatternId
        PickId = self.PickId    
        limits = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Limits'))[1]
        controls,bravais,cells,dmin = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
        inst = dict(zip(inst[3],inst[1]))
        if 'Lam' in inst:
            wave = inst['Lam']
        else:
            wave = inst['Lam1']
        dmin = wave/(2.0*sind(limits[1]/2.0))
        self.HKL = G2lat.GenHBravais(dmin,ibrav,A)
        for hkl in self.HKL:
            hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1])             
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
        
    def CopyUnitCell(event):
        controls,bravais,cells,dmin = self.PatternTree.GetItemPyData(UnitCellsId)
        for Cell in cells:
            if Cell[-1]:
                break
        cell = Cell[2:9]
        controls[4] = 1
        controls[5] = bravaisSymb[cell[0]]
        controls[6:12] = cell[1:8]
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        self.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin])
        UpdateUnitCellsGrid(self,data)
        
        self.dataFrame.RefineCell.Enable(True)
                
    def RefineCell(event):
        def cellPrint(ibrav,A):
            cell = G2lat.A2cell(A)
            Vol = G2lat.calc_V(A)
            if ibrav in [0,1,2]:
                print "%s%10.6f" % ('a =',cell[0])
            elif ibrav in [3,4,5,6]:
                print "%s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],' c =',cell[2],' volume =',Vol)
            elif ibrav in [7,8,9,10]:
                print "%s%10.6f %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],' volume =',Vol)
            elif ibrav in [11,12]:
                print "%s%10.6f %s%10.6f %s%10.6f %s%8.3f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],'beta =',cell[4],' volume =',Vol)
            else:
                print "%s%10.6f %s%10.6f %s%10.6f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2])
                print "%s%8.3f %s%8.3f %s%8.3f %s%12.3f" % ('alpha =',cell[3],'beta =',cell[4],'gamma =',cell[5],' volume =',Vol)
             
        PatternId = self.PatternId
        PickId = self.PickId    
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
        if not peaks:
            self.ErrorDialog('No peaks!', 'Nothing to refine!')
            return        
        print 'Refine cell'
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))[1]
        controls,bravais,cells,dmin = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        dmin = G2indx.getDmin(peaks)-0.005
        self.HKL = G2lat.GenHBravais(dmin,ibrav,A)
        G2indx.IndexPeaks(peaks,self.HKL)
        if controls[0]:
            Lhkl,M20,X20,Aref,Zero = G2indx.refinePeaksZ(peaks,inst[1],ibrav,A,controls[1])            
            controls[1] = Zero
        else:
            Lhkl,M20,X20,Aref = G2indx.refinePeaks(peaks,ibrav,A)
        controls[6:12] = G2lat.A2cell(Aref)
        controls[12] = G2lat.calc_V(Aref)
        data = [controls,bravais,cells,dmin]
        cells = self.PatternTree.GetItemPyData(UnitCellsId)[2]
        for cell in cells:
            cell[-1] = False
        cells.insert(0,[M20,X20,ibrav]+controls[6:13]+[True,])
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'),data)
        self.HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
        UpdateUnitCellsGrid(self,data)
        print "%s%10.3f" % ('refinement M20 = ',M20)
        print 'unindexed lines = ',X20
        cellPrint(ibrav,Aref)
        for hkl in self.HKL:
            hkl.append(2.0*asind(inst[1]/(2.*hkl[3]))+controls[1])             
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
        
    def IndexPeaks(event):
        PatternId = self.PatternId    
#        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
#        if not peaks:
#            self.ErrorDialog('No peaks!', 'Nothing to index!')
#            return
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))[1]
        print 'Peak Indexing'
        try:
            controls,bravais,cells,dmin = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'))
            cells = []
        except ValueError:
            self.ErrorDialog('Error','Need to set controls in Unit Cell List first')
            return
        if True not in bravais:
            self.ErrorDialog('Error','No Bravais lattices selected')
            return
        self.dataFrame.CopyCell.Enable(False)
        self.dataFrame.RefineCell.Enable(False)
        OK,dmin,cells = G2indx.DoIndexPeaks(peaks,inst,controls,bravais)
        if OK:
            data = [controls,bravais,cells,dmin]
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'),data)
            UpdateUnitCellsGrid(self,data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                self.HKL = G2lat.GenHBravais(dmin,bestCell[2],G2lat.cell2A(bestCell[3:9]))
                for hkl in self.HKL:
                    hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
                if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
                    G2plt.PlotPowderLines(self)
                else:
                    G2plt.PlotPatterns(self)
            self.dataFrame.CopyCell.Enable(True)
            self.dataFrame.IndexPeaks.Enable(True)
            self.dataFrame.MakeNewPhase.Enable(True)
            UpdateUnitCellsGrid(self,data)
                
    def RefreshUnitCellsGrid(event):
        cells,dmin = self.PatternTree.GetItemPyData(UnitCellsId)[2:]
        r,c =  event.GetRow(),event.GetCol()
        if cells:
            if c == 2:
                for i in range(len(cells)):
                    cells[i][-1] = False
                    UnitCellsTable.SetValue(i,c,False)
                UnitCellsTable.SetValue(r,c,True)
                gridDisplay.ForceRefresh()
                cells[r][-1] = True
                ibrav = cells[r][2]
                A = G2lat.cell2A(cells[r][3:9])
                self.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                for hkl in self.HKL:
                    hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))
                if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
                    G2plt.PlotPowderLines(self)
                else:
                    G2plt.PlotPatterns(self)
        
    def MakeNewPhase(event):
        if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PhaseName = dlg.GetValue()
                cells = self.PatternTree.GetItemPyData(UnitCellsId)[2]
                for Cell in cells:
                    if Cell[-1]:
                        break
                cell = Cell[2:10]        
                sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
                E,SGData = G2spc.SpcGroup(spaceGroups[cell[0]])
                self.PatternTree.SetItemPyData(sub, \
                    {'General':{'Name':PhaseName,'Type':'nuclear','SGData':SGData,
                    'Cell':[False,]+cell[1:],
                    'Pawley dmin':1.00},'Atoms':[],'Drawing':{},'Histograms':{}})
                Status.SetStatusText('Change space group if needed')
        finally:
            dlg.Destroy()
            
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataFrame.SetMenuBar(self.dataFrame.IndexMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    self.Bind(wx.EVT_MENU, IndexPeaks, id=G2gd.wxID_INDEXPEAKS)
    self.Bind(wx.EVT_MENU, CopyUnitCell, id=G2gd.wxID_COPYCELL)
    self.Bind(wx.EVT_MENU, RefineCell, id=G2gd.wxID_REFINECELL)
    self.Bind(wx.EVT_MENU, MakeNewPhase, id=G2gd.wxID_MAKENEWPHASE)
    
    controls,bravais,cells,dmin = data
    if len(controls) < 13:              #add cell volume if missing
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    self.PatternTree.SetItemPyData(UnitCellsId,data)            #update with volume
    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
    bravaisNames = ['Cubic-F','Cubic-I','Cubic-P','Trigonal-R','Trigonal/Hexagonal-P',
        'Tetragonal-I','Tetragonal-P','Orthorhombic-F','Orthorhombic-I','Orthorhombic-C',
        'Orthorhombic-P','Monoclinic-C','Monoclinic-P','Triclinic']
    cellGUIlist = [[[0,1,2],4,zip([" Unit cell: a = "," Vol = "],["%.5f","%.3f"],[True,False],[0,0])],
    [[3,4,5,6],6,zip([" Unit cell: a = "," c = "," Vol = "],["%.5f","%.5f","%.3f"],[True,True,False],[0,2,0])],
    [[7,8,9,10],8,zip([" Unit cell: a = "," b = "," c = "," Vol = "],["%.5f","%.5f","%.5f","%.3f"],
        [True,True,True,False],[0,1,2,0])],
    [[11,12],10,zip([" Unit cell: a = "," b = "," c = "," beta = "," Vol = "],
        ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
    [[13,],8,zip([" Unit cell: a = "," b = "," c = "," Vol = "," alpha = "," beta = "," gamma = "],
        ["%.5f","%.5f","%.5f","%.3f","%.3f","%.3f","%.3f"],
        [True,True,True,False,True,True,True],[0,1,2,0,3,4,5])]]
    
    self.dataFrame.SetLabel('Unit Cells List')
    self.sp = wx.SplitterWindow(self.dataFrame)
    self.dataDisplay = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
    self.dataFrame.IndexPeaks.Enable(False)
    peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Index Peak List'))
    if peaks:
        self.dataFrame.IndexPeaks.Enable(True)
    self.dataFrame.RefineCell.Enable(False)
    if controls[12] > 1.0:                               #if a "real" volume (i.e. not default)
        self.dataFrame.RefineCell.Enable(True)    
    self.dataFrame.CopyCell.Enable(False)
    self.dataFrame.MakeNewPhase.Enable(False)        
    if cells:
        self.bottom = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
        self.sp.SplitHorizontally(self.dataDisplay,self.bottom,0)
        self.dataFrame.CopyCell.Enable(True)
        self.dataFrame.MakeNewPhase.Enable(True)        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Indexing controls: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,5,5,5)
    littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Max Nc/Nobs '),0,wx.ALIGN_CENTER_VERTICAL)
    NcNo = wx.SpinCtrl(self.dataDisplay)
    NcNo.SetRange(1,6)
    NcNo.SetValue(controls[2])
    NcNo.Bind(wx.EVT_SPINCTRL,OnNcNo)
    littleSizer.Add(NcNo,0,wx.ALIGN_CENTER_VERTICAL)
    littleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Start Volume '),0,wx.ALIGN_CENTER_VERTICAL)
    startVol = wx.TextCtrl(self.dataDisplay,value=str('%d'%(controls[3])),style=wx.TE_PROCESS_ENTER)
    startVol.Bind(wx.EVT_TEXT_ENTER,OnStartVol)
    startVol.Bind(wx.EVT_KILL_FOCUS,OnStartVol)
    littleSizer.Add(startVol,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Select Bravais Lattices for indexing: '),
        0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,7,5,5)
    bravList = []
    bravs = zip(bravais,bravaisNames)
    for brav,bravName in bravs:
        bravCk = wx.CheckBox(self.dataDisplay,label=bravName)
        bravList.append(bravCk.GetId())
        bravCk.SetValue(brav)
        bravCk.Bind(wx.EVT_CHECKBOX,OnBravais)
        littleSizer.Add(bravCk,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Cell Refinement: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    littleSizer.Add(wx.StaticText(self.dataDisplay,label=" Bravais lattice "),0,wx.ALIGN_CENTER_VERTICAL)
    bravSel = wx.Choice(self.dataDisplay,choices=bravaisSymb)
    bravSel.SetSelection(bravaisSymb.index(controls[5]))
    bravSel.Bind(wx.EVT_CHOICE,OnBravSel)
    littleSizer.Add(bravSel,0,wx.ALIGN_CENTER_VERTICAL)
    littleSizer.Add(wx.StaticText(self.dataDisplay,label=" Zero offset"),0,wx.ALIGN_CENTER_VERTICAL)
    zero = wx.TextCtrl(self.dataDisplay,value="%.4f"%(controls[1]),style=wx.TE_PROCESS_ENTER)
    zero.Bind(wx.EVT_TEXT_ENTER,OnZero)
    zero.Bind(wx.EVT_KILL_FOCUS,OnZero)
    littleSizer.Add(zero,0,wx.ALIGN_CENTER_VERTICAL)
    zeroVar = wx.CheckBox(self.dataDisplay,label="Refine?")
    zeroVar.SetValue(controls[0])
    zeroVar.Bind(wx.EVT_CHECKBOX,OnZeroVar)
    littleSizer.Add(zeroVar,0,wx.ALIGN_CENTER_VERTICAL)
    hklShow = wx.CheckBox(self.dataDisplay,label="  Show hkl positions")
    hklShow.Bind(wx.EVT_CHECKBOX,OnHklShow)
    littleSizer.Add(hklShow,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    
    mainSizer.Add((5,5),0)
    ibrav = SetLattice(controls)
    for cellGUI in cellGUIlist:
        if ibrav in cellGUI[0]:
            useGUI = cellGUI
    cellList = []
    littleSizer = wx.FlexGridSizer(2,useGUI[1],5,5)
    for txt,fmt,ifEdit,Id in useGUI[2]:
        littleSizer.Add(wx.StaticText(self.dataDisplay,label=txt),0,wx.ALIGN_CENTER_VERTICAL)
        if ifEdit:          #a,b,c,etc.
            cellVal = wx.TextCtrl(self.dataDisplay,value=(fmt%(controls[6+Id])),style=wx.TE_PROCESS_ENTER)
            cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
            cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
            littleSizer.Add(cellVal,0,wx.ALIGN_CENTER_VERTICAL)
            cellList.append(cellVal.GetId())
        else:               #volume
            volVal = wx.TextCtrl(self.dataDisplay,value=(fmt%(controls[12])),style=wx.TE_READONLY)
            volVal.SetBackgroundColour(VERY_LIGHT_GREY)
            littleSizer.Add(volVal,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
        
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    topSize = mainSizer.Fit(self.dataFrame)
    self.dataDisplay.SetSize(topSize)
    if cells:
        if ibrav == 13:
            topSize[1] += 230
        else:
            topSize[1] += 200
    self.dataFrame.setSizePosLeft(topSize)
    
    
    if cells:
        bottomSize = self.bottom.GetSize()
        if ibrav == 13:
            bottomSize[1] -= 240
        else:
            bottomSize[1] -= 210
        wx.StaticText(parent=self.bottom,label=' Indexing Result ')
        rowLabels = []
        colLabels = ['M20','X20','use','Bravais','a','b','c','alpha','beta','gamma','Volume']
        Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+ \
            [wg.GRID_VALUE_FLOAT+':10,2']
        numRows = len(cells)
        table = []
        for cell in cells:
            rowLabels.append('')
            row = cell[0:2]+[cell[-1]]+[bravaisSymb[cell[2]]]+cell[3:10]
            if cell[-1]:
                A = G2lat.cell2A(cell[3:9])
                self.HKL = G2lat.GenHBravais(dmin,cell[2],A)
                for hkl in self.HKL:
                    hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))
            table.append(row)
        UnitCellsTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        gridDisplay = G2gd.GSGrid(self.bottom)
        gridDisplay.SetPosition(wx.Point(0,20))                
        gridDisplay.SetTable(UnitCellsTable, True)
        self.dataFrame.CopyCell.Enable(True)
        gridDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,RefreshUnitCellsGrid)
        gridDisplay.SetMargins(0,0)
        gridDisplay.SetRowLabelSize(0)
        gridDisplay.AutoSizeColumns(False)
        for r in range(gridDisplay.GetNumberRows()):
            for c in range(gridDisplay.GetNumberCols()):
                if c == 2:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=True)
        gridDisplay.SetSize(bottomSize)

def UpdateReflectionGrid(self,data):
    if not data:
        print 'No phases, no reflections'
        return
    phases = data.keys()
    
    def OnSelectPhase(event):
        dlg = wx.SingleChoiceDialog(self,'Select','Phase',phases)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                self.RefList = phases[sel]
                UpdateReflectionGrid(self,data)
        finally:
            dlg.Destroy()
        G2plt.PlotPatterns(self)
        
        
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataFrame.SetMenuBar(self.dataFrame.ReflMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()    
    self.Bind(wx.EVT_MENU, OnSelectPhase, id=G2gd.wxID_SELECTPHASE)
    self.dataFrame.SelectPhase.Enable(False)
    if len(data) > 1:
        self.dataFrame.SelectPhase.Enable(True)
    rowLabels = []
    refList = [refl[:11] for refl in data[self.RefList]]
    for i in range(len(refList)): rowLabels.append(str(i))
    colLabels = ['H','K','L','mul','d','pos','sig','gam','Fosq','Fcsq','phase',]
    Types = 4*[wg.GRID_VALUE_LONG,]+4*[wg.GRID_VALUE_FLOAT+':10,4',]+2*[wg.GRID_VALUE_FLOAT+':10,2',]+[wg.GRID_VALUE_FLOAT+':10,3',]
    self.PeakTable = G2gd.Table(refList,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Reflection List for '+self.RefList)
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.PeakTable, True)
    self.dataDisplay.EnableEditing(False)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([555,350])

def UpdatePDFGrid(self,data):
    global inst
    tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
    dataFile = self.PatternTree.GetItemText(self.PatternId)
    powName = 'PWDR'+dataFile[4:]
    powId = G2gd.GetPatternTreeItemId(self,self.root, powName)
    fullLimits,limits = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,powId, 'Limits'))
    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,powId, 'Instrument Parameters'))
    inst = dict(zip(inst[3],inst[1]))
    if 'Lam' in inst:
        keV = 12.397639/inst['Lam']
    else:
        keV = 12.397639/inst['Lam1']
    wave = 12.397639/keV
    qLimits = [tth2q(fullLimits[0],wave),tth2q(fullLimits[1],wave)]
    data['QScaleLim'][1] = min(qLimits[1],data['QScaleLim'][1])
    if data['QScaleLim'][0]:
        data['QScaleLim'][0] = max(qLimits[0],data['QScaleLim'][0])
    else:                                #initial setting at 90% of max Q
        data['QScaleLim'][0] = 0.90*data['QScaleLim'][1]
    polariz = inst['Polariz.']
    azimuth = inst['Azimuth']
    itemDict = {}
    
    def FillFileSizer(fileSizer,key):
        #fileSizer is a FlexGridSizer(3,6)
        
        def OnSelectFile(event):
            Obj = event.GetEventObject()
            fileKey,itemKey,fmt = itemDict[Obj.GetId()]
            if itemKey == 'Name':
                value = Obj.GetValue()
            Obj.SetValue(fmt%(value))
            data[fileKey][itemKey] = value
            UpdatePDFGrid(self,data)
        
        def OnValueChange(event):
            Obj = event.GetEventObject()
            fileKey,itemKey,fmt = itemDict[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = -1.0
            Obj.SetValue(fmt%(value))
            data[fileKey][itemKey] = value
            auxPlot = ComputePDF(data)
            G2plt.PlotISFG(self,newPlot=True)
                        
        item = data[key]
        fileList = np.array(GetFileList('PWDR')).T[1]
        fileSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' '+key+' file:'),0,wx.ALIGN_CENTER_VERTICAL)
        fileName = wx.ComboBox(self.dataDisplay,value=item['Name'],choices=fileList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        itemDict[fileName.GetId()] = [key,'Name','%s']
        fileName.Bind(wx.EVT_COMBOBOX,OnSelectFile)        
        fileSizer.Add(fileName,0,)
        fileSizer.Add(wx.StaticText(parent=self.dataDisplay,label='Multiplier:'),0,wx.ALIGN_CENTER_VERTICAL)
        mult = wx.TextCtrl(self.dataDisplay,value='%.3f'%(item['Mult']),style=wx.TE_PROCESS_ENTER)
        itemDict[mult.GetId()] = [key,'Mult','%.3f']
        mult.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
        mult.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
        fileSizer.Add(mult,0,)
        fileSizer.Add(wx.StaticText(parent=self.dataDisplay,label='Add:'),0,wx.ALIGN_CENTER_VERTICAL)
        add = wx.TextCtrl(self.dataDisplay,value='%.0f'%(item['Add']),style=wx.TE_PROCESS_ENTER)
        itemDict[add.GetId()] = [key,'Add','%.0f']
        add.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
        add.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
        fileSizer.Add(add,0,)
        
    def SumElementVolumes():
        sumVol = 0.
        ElList = data['ElList']
        for El in ElList:
            Avol = (4.*math.pi/3.)*ElList[El]['Drad']**3
            sumVol += Avol*ElList[El]['FormulaNo']
        return sumVol
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)        
        
    def FillElemSizer(elemSizer,ElData):
        
        def OnFractionChange(event):
            try:
                value = max(0.0,float(num.GetValue()))
            except ValueError:
                value = 0.0
            num.SetValue('%.3f'%(value))
            ElData['FormulaNo'] = value
            data['Form Vol'] = max(10.0,SumElementVolumes())
            formVol.SetValue('%.2f'%(data['Form Vol']))
            UpdatePDFGrid(self,data)
            auxPlot = ComputePDF(data)
            G2plt.PlotISFG(self,newPlot=True)        
        
        elemSizer.Add(wx.StaticText(parent=self.dataDisplay,
            label=' Element: '+'%2s'%(ElData['Symbol'])+' * '),0,wx.ALIGN_CENTER_VERTICAL)
        num = wx.TextCtrl(self.dataDisplay,value='%.3f'%(ElData['FormulaNo']),style=wx.TE_PROCESS_ENTER)
        num.Bind(wx.EVT_TEXT_ENTER,OnFractionChange)        
        num.Bind(wx.EVT_KILL_FOCUS,OnFractionChange)
        elemSizer.Add(num,0,wx.ALIGN_CENTER_VERTICAL)
        elemSizer.Add(wx.StaticText(parent=self.dataDisplay,
            label="f': %.3f"%(ElData['fp'])+' f": %.3f'%(ElData['fpp'])+' mu: %.2f barns'%(ElData['mu']) ),
            0,wx.ALIGN_CENTER_VERTICAL)
            
    def OnGeometry(event):
        data['Geometry'] = geometry.GetValue()
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)        
        
    def OnDetType(event):
        data['DetType'] = detType.GetValue()
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)        
        
    def OnFormVol(event):
        try:
            value = float(formVol.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Form Vol']
        data['Form Vol'] = value
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)        
        
    def OnDiameter(event):
        try:
            value = float(diam.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Diam']
        data['Diam'] = value
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
        
    def OnPolaVal(event):
        try:
            value = float(polaVal.GetValue())
            if not (0.0 <= value <= 1.0):
                raise ValueError
        except ValueError:
            value = inst['Polariz.']
        inst['Polariz.'] = value
        polaVal.SetValue('%.2f'%(inst['Polariz.']))
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
                
    def OnAzimVal(event):
        try:
            value = float(azimVal.GetValue())
            if not (0. <= value <= 360.):
                raise ValueError
        except ValueError:
            value = inst['Azimuth']
        inst['Azimuth'] = value
        azimVal.SetValue('%.1f'%(inst['Azimuth']))
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
                        
    def OnObliqCoeff(event):
        try:
            value = float(obliqCoeff.GetValue())
            if value < 0.0:
                raise ValueError
            elif value > 1.0:
                value = 1.0
        except ValueError:
            value = data['ObliqCoeff']
        data['ObliqCoeff'] = value
        obliqCoeff.SetValue('%.3f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
        
    def OnRulandWdt(event):
        try:
            value = float(rulandWdt.GetValue())
            if value <= 0.001:
                raise ValueError
            elif value > 1.0:
                value = 1.0
        except ValueError:
            value = data['Ruland']
        data['Ruland'] = value
        rulandWdt.SetValue('%.3f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
        
    def OnRulSlider(event):
        value = int(rulandSldr.GetValue())/1000.
        data['Ruland'] = max(0.001,value)
        rulandWdt.SetValue('%.3f'%(data['Ruland']))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)
        
    def OnLorch(event):
        data['Lorch'] = lorch.GetValue()
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)        
                        
    def OnPacking(event):
        try:
            value = float(pack.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Pack']
        data['Pack'] = value
        UpdatePDFGrid(self,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=False)        
                
    def OnSQmin(event):
        try:
            value = float(SQmin.GetValue())
            if value < qLimits[0]:
                raise ValueError
        except ValueError:
            value = max(qLimits[0],data['QScaleLim'][0])
        data['QScaleLim'][0] = value
        SQmin.SetValue('%.1f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)        
        
    def OnSQmax(event):
        try:
            value = float(SQmax.GetValue())
            if value > qLimits[1]:
                raise ValueError
        except ValueError:
            value = min(qLimits[1],data['QScaleLim'][1])
        data['QScaleLim'][1] = value
        if value < data['QScaleLim'][0]:
            data['QScaleLim'][0] = 0.90*value
            SQmin.SetValue('%.1f'%(data['QScaleLim'][0]))
        SQmax.SetValue('%.1f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)
        
    def OnResetQ(event):
        resetQ.SetValue(False)
        data['QScaleLim'][1] = qLimits[1]
        SQmax.SetValue('%.1f'%(data['QScaleLim'][1]))
        data['QScaleLim'][0] = 0.9*qLimits[1]
        SQmin.SetValue('%.1f'%(data['QScaleLim'][0]))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(self,newPlot=True)        

    def GetFileList(fileType,skip=None):
        fileList = [[False,'',0]]
        Source = ''
        id, cookie = self.PatternTree.GetFirstChild(self.root)
        while id:
            name = self.PatternTree.GetItemText(id)
            if fileType in name:
                if id == skip:
                    Source = name
                else:
                    fileList.append([False,name,id])
            id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
        
    def OnCopyPDFControls(event):
        import copy
        TextList,Source = GetFileList('PDF',skip=self.PatternId)
        TextList[0] = [False,'All PDF',0]
        if len(TextList) == 1:
            self.ErrorDialog('Nothing to copy controls to','There must be more than one "PDF" pattern')
            return
        dlg = self.CopyDialog(self,'Copy PDF controls','Copy controls from '+Source+' to:',TextList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                if result[0][0]:
                    result = TextList[1:]
                    for item in result: item[0] = True
                for i,item in enumerate(result):
                    ifcopy,name,id = item
                    if ifcopy:
                        olddata = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'PDF Controls'))
                        sample = olddata['Sample']
                        olddata.update(data)
                        olddata['Sample'] = sample
                        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,id, 'PDF Controls'),olddata)
                Status.SetStatusText('PDF controls copied')
        finally:
            dlg.Destroy()
                
    def OnSavePDFControls(event):
        print 'save PDF controls?'
        
    def OnLoadPDFControls(event):
        print 'Load PDF controls?'
        
    def OnAddElement(event):
        ElList = data['ElList']
        PE = G2elem.PickElement(self,oneOnly=True)
        if PE.ShowModal() == wx.ID_OK:
            El = PE.Elem
            if El not in ElList:
                ElemSym = El.strip().upper()                
                FpMu = G2elem.FPcalc(G2elem.GetXsectionCoeff(ElemSym), keV)
                ElData = G2elem.GetFormFactorCoeff(ElemSym)[0]
                ElData['FormulaNo'] = 0.0
                ElData.update(G2elem.GetAtomInfo(ElemSym))
                ElData.update(dict(zip(['fp','fpp','mu'],FpMu)))
                ElData.update(G2elem.GetFFC5(El))
                data['ElList'][El] = ElData
            data['Form Vol'] = max(10.0,SumElementVolumes())
        PE.Destroy()
        UpdatePDFGrid(self,data)
        
    def OnDeleteElement(event):
        ElList = data['ElList']
        choice = ElList.keys()
        dlg = G2elem.DeleteElement(self,choice=choice)
        if dlg.ShowModal() == wx.ID_OK:
            del ElList[dlg.GetDeleteElement()]
        dlg.Destroy()
        UpdatePDFGrid(self,data)
                
    def ComputePDF(Data):
        xydata = {}
        for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
            name = Data[key]['Name']
            if name:
                xydata[key] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.root,name))
                PDFname = name
        powName = xydata['Sample'][2]
        powId = G2gd.GetPatternTreeItemId(self,self.root,powName)
        inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,powId,'Instrument Parameters'))
        inst = dict(zip(inst[3],inst[1]))
        auxPlot = G2pwd.CalcPDF(Data,inst,xydata)
        PDFId = G2gd.GetPatternTreeItemId(self,self.root,'PDF '+powName[4:])
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PDFId,'I(Q)'+powName[4:]),xydata['IofQ'])
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PDFId,'S(Q)'+powName[4:]),xydata['SofQ'])
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PDFId,'F(Q)'+powName[4:]),xydata['FofQ'])
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PDFId,'G(R)'+powName[4:]),xydata['GofR'])
        return auxPlot
        
    def OnComputePDF(event):
        print 'Calculating PDF:'
        auxPlot = ComputePDF(data)
        print 'Done calculating PDF:'
        Status.SetStatusText('PDF computed')
        for plot in auxPlot:
            G2plt.PlotXY(self,plot[:2],type=plot[2])
        
        G2plt.PlotISFG(self,newPlot=True,type='I(Q)')
        G2plt.PlotISFG(self,newPlot=True,type='S(Q)')
        G2plt.PlotISFG(self,newPlot=True,type='F(Q)')
        G2plt.PlotISFG(self,newPlot=True,type='G(R)')
        
    def OnComputeAllPDF(event):
        print 'Calculating PDFs:'
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                Name = self.PatternTree.GetItemText(id)
                if 'PDF' in Name:
                    Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,id,'PDF Controls'))
                    auxPlot = ComputePDF(Data)                    
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            Status.SetStatusText('All PDFs computed')
            G2plt.PlotISFG(self,newPlot=True,type='G(R)')
            print ' Done calculating PDFs:'
        
    def OnShowTip(self,tip):
        print tip    
                
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataFrame.SetMenuBar(self.dataFrame.PDFMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()    
    self.dataDisplay = wx.Panel(self.dataFrame)
    self.dataFrame.Bind(wx.EVT_MENU, OnCopyPDFControls, id=G2gd.wxID_PDFCOPYCONTROLS)
    self.dataFrame.Bind(wx.EVT_MENU, OnSavePDFControls, id=G2gd.wxID_PDFSAVECONTROLS)
    self.dataFrame.Bind(wx.EVT_MENU, OnLoadPDFControls, id=G2gd.wxID_PDFLOADCONTROLS)
    self.dataFrame.Bind(wx.EVT_MENU, OnAddElement, id=G2gd.wxID_PDFADDELEMENT)
    self.dataFrame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2gd.wxID_PDFDELELEMENT)
    self.dataFrame.Bind(wx.EVT_MENU, OnComputePDF, id=G2gd.wxID_PDFCOMPUTE)
    self.dataFrame.Bind(wx.EVT_MENU, OnComputeAllPDF, id=G2gd.wxID_PDFCOMPUTEALL)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' PDF data files: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    str = ' Sample file: PWDR %s   Wavelength, A: %.5f  Energy, keV: %.3f  Polariz.: %.2f '%(dataFile[3:],wave,keV,polariz)
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=str),0,wx.ALIGN_CENTER_VERTICAL)
#    dataSizer = wx.BoxSizer(wx.HORIZONTAL)
#    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label='Azimuth'),0,wx.ALIGN_CENTER_VERTICAL)
#    azimVal = wx.TextCtrl(self.dataDisplay,value='%.2f'%(inst['Azimuth']))
#    azimVal.Bind(wx.EVT_TEXT_ENTER,OnAzimVal)        
#    azimVal.Bind(wx.EVT_KILL_FOCUS,OnAzimVal)
#    dataSizer.Add(azimVal,0)    
#    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label='Polarization'),0,wx.ALIGN_CENTER_VERTICAL)
#    polaVal = wx.TextCtrl(self.dataDisplay,value='%.2f'%(inst['Polariz.']))
#    polaVal.Bind(wx.EVT_TEXT_ENTER,OnPolaVal)        
#    polaVal.Bind(wx.EVT_KILL_FOCUS,OnPolaVal)
#    dataSizer.Add(polaVal,0)    
#    mainSizer.Add(dataSizer,0)
    mainSizer.Add((5,5),0)
    fileSizer = wx.FlexGridSizer(3,6,5,1)
    select = ['Sample Bkg.','Container']
    if data['Container']['Name']:
        select.append('Container Bkg.')
    for key in select:
        FillFileSizer(fileSizer,key)
    mainSizer.Add(fileSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Sample information: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)    

    ElList = data['ElList']
    Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
    Trans = G2pwd.Transmission(data['Geometry'],Abs*data['Pack'],data['Diam'])
    elemSizer = wx.FlexGridSizer(3,3,5,1)
    for El in ElList:
        FillElemSizer(elemSizer,ElList[El])
    mainSizer.Add(elemSizer,0)
    mainSizer.Add((5,5),0)    
    midSizer = wx.BoxSizer(wx.HORIZONTAL)
    midSizer.Add(wx.StaticText(self.dataDisplay,label=' Formula volume: '),0,wx.ALIGN_CENTER_VERTICAL)
    formVol = wx.TextCtrl(self.dataDisplay,value='%.2f'%(data['Form Vol']))
    formVol.Bind(wx.EVT_TEXT_ENTER,OnFormVol)        
    formVol.Bind(wx.EVT_KILL_FOCUS,OnFormVol)
    midSizer.Add(formVol,0)
    midSizer.Add(wx.StaticText(self.dataDisplay,
        label=' Theoretical absorption: %.4f cm-1 Sample absorption: %.4f cm-1'%(Abs,Abs*data['Pack'])),
        0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(midSizer,0)
    mainSizer.Add((5,5),0)    

    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(self.dataDisplay,label=' Sample geometry: '),0,wx.ALIGN_CENTER_VERTICAL)
    choice = ['Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate']
    geometry = wx.ComboBox(self.dataDisplay,value=data['Geometry'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    geometry.Bind(wx.EVT_COMBOBOX, OnGeometry)
    geoBox.Add(geometry,0)
    geoBox.Add(wx.StaticText(self.dataDisplay,label=' Sample diameter/thickness, mm: '),0,wx.ALIGN_CENTER_VERTICAL)
    diam = wx.TextCtrl(self.dataDisplay,value='%.3f'%(data['Diam']))
    diam.Bind(wx.EVT_TEXT_ENTER,OnDiameter)        
    diam.Bind(wx.EVT_KILL_FOCUS,OnDiameter)
#    diam.Bind(wx.EVT_SET_FOCUS,OnShowTip(self,'tip')) #this doesn't work - what would????
    geoBox.Add(diam,0)
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(self.dataDisplay,label=' Packing: '),0,wx.ALIGN_CENTER_VERTICAL)
    pack = wx.TextCtrl(self.dataDisplay,value='%.2f'%(data['Pack']))
    pack.Bind(wx.EVT_TEXT_ENTER,OnPacking)        
    pack.Bind(wx.EVT_KILL_FOCUS,OnPacking)
    geoBox.Add(pack,0)
    geoBox.Add(wx.StaticText(self.dataDisplay,label=' Sample transmission: %.3f %%'%(Trans)),0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
        
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' S(Q)->F(Q)->G(R) controls: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(self.dataDisplay,label=' Detector type: '),0,wx.ALIGN_CENTER_VERTICAL)
    choice = ['Image plate','Point detector']
    detType = wx.ComboBox(self.dataDisplay,value=data['DetType'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    detType.Bind(wx.EVT_COMBOBOX, OnDetType)
    sqBox.Add(detType,0)
    if data['DetType'] == 'Image plate':
        sqBox.Add(wx.StaticText(self.dataDisplay,label=' IP transmission coeff.: '),0,wx.ALIGN_CENTER_VERTICAL)
        obliqCoeff = wx.TextCtrl(self.dataDisplay,value='%.3f'%(data['ObliqCoeff']))
        obliqCoeff.Bind(wx.EVT_TEXT_ENTER,OnObliqCoeff)        
        obliqCoeff.Bind(wx.EVT_KILL_FOCUS,OnObliqCoeff)
        sqBox.Add(obliqCoeff,0)
    mainSizer.Add(sqBox,0)
        
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(self.dataDisplay,label=' Ruland width: '),0,wx.ALIGN_CENTER_VERTICAL)    
    rulandSldr = wx.Slider(parent=self.dataDisplay,style=wx.SL_HORIZONTAL,
        value=int(1000*data['Ruland']))
    sqBox.Add(rulandSldr,1,wx.EXPAND)
    rulandSldr.Bind(wx.EVT_SLIDER, OnRulSlider)
    rulandWdt = wx.TextCtrl(self.dataDisplay,value='%.3f'%(data['Ruland']))
    rulandWdt.Bind(wx.EVT_TEXT_ENTER,OnRulandWdt)        
    rulandWdt.Bind(wx.EVT_KILL_FOCUS,OnRulandWdt)
    sqBox.Add(rulandWdt,0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(sqBox,0,wx.ALIGN_LEFT|wx.EXPAND)
    
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    lorch = wx.CheckBox(parent=self.dataDisplay,label='Lorch damping?')
    lorch.SetValue(data['Lorch'])
    lorch.Bind(wx.EVT_CHECKBOX, OnLorch)
    sqBox.Add(lorch,0,wx.ALIGN_CENTER_VERTICAL)
    sqBox.Add(wx.StaticText(self.dataDisplay,label=' Scaling q-range: '),0,wx.ALIGN_CENTER_VERTICAL)
    SQmin = wx.TextCtrl(self.dataDisplay,value='%.1f'%(data['QScaleLim'][0]))
    SQmin.Bind(wx.EVT_TEXT_ENTER,OnSQmin)        
    SQmin.Bind(wx.EVT_KILL_FOCUS,OnSQmin)    
    sqBox.Add(SQmin,0)
    sqBox.Add(wx.StaticText(self.dataDisplay,label=' to '),0,wx.ALIGN_CENTER_VERTICAL)
    SQmax = wx.TextCtrl(self.dataDisplay,value='%.1f'%(data['QScaleLim'][1]))
    SQmax.Bind(wx.EVT_TEXT_ENTER,OnSQmax)        
    SQmax.Bind(wx.EVT_KILL_FOCUS,OnSQmax)
    sqBox.Add(SQmax,0)
    resetQ = wx.CheckBox(parent=self.dataDisplay,label='Reset?')
    sqBox.Add(resetQ,0)
    resetQ.Bind(wx.EVT_CHECKBOX, OnResetQ)
    
    mainSizer.Add(sqBox,0)

    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(self.dataFrame)
    self.dataDisplay.SetSize(Size)
    self.dataFrame.setSizePosLeft(Size)
    