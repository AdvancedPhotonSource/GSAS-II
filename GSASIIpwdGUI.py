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
import copy
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
    
def IsHistogramInAnyPhase(G2frame,histoName):
    phases = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
    if phases:
        item, cookie = G2frame.PatternTree.GetFirstChild(phases)
        while item:
            data = G2frame.PatternTree.GetItemPyData(item)
            histoList = data['Histograms'].keys()
            if histoName in histoList:
                return True
            item, cookie = G2frame.PatternTree.GetNextChild(phases, cookie)
        return False
    else:
        return False

def SetDefaultSample():
    return {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],
        'DisplaceX':[0.0,False],'DisplaceY':[0.0,False],'Diffuse':[],
        'Temperature':300.,'Pressure':1.0,'Humidity':0.0,
        'Voltage':0.0,'Force':0.0,'Gonio. radius':200.0,
        'Omega':0.0,'Chi':0.0,'Phi':0.0}    
       
def UpdatePeakGrid(G2frame, data):
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    
    def OnUnDo(event):
        DoUnDo()
        G2frame.dataFrame.UnDo.Enable(False)
        
    def DoUnDo():
        print 'Undo last refinement'
        file = open(G2frame.undofile,'rb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item),cPickle.load(file))
            if G2frame.dataDisplay.GetName() == item:
                if item == 'Background':
                    UpdateBackground(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
                elif item == 'Instrument Parameters':
                    UpdateInstrumentGrid(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
                elif item == 'Peak List':
                    UpdatePeakGrid(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
            print item,' recovered'
        file.close()
        
    def SaveState():
        G2frame.undofile = G2frame.dirname+'\\GSASII.save'
        file = open(G2frame.undofile,'wb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            cPickle.dump(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,item)),file,1)
        file.close()
        G2frame.dataFrame.UnDo.Enable(True)
        
    def OnLSQPeakFit(event):
        if not G2frame.GSASprojectfile:            #force a save of the gpx file so SaveState can wirte in the same directory
            G2frame.OnFileSaveas(event)
        OnPeakFit('LSQ')
        
    def OnOneCycle(event):
        OnPeakFit('LSQ',oneCycle=True)
        
    def OnClearPeaks(event):
        dlg = wx.MessageDialog(G2frame,'Delete all peaks?','Clear peak list',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = []
        finally:
            dlg.Destroy()
        UpdatePeakGrid(G2frame,peaks)
        G2plt.PlotPatterns(G2frame)
        
    def OnPeakFit(FitPgm,oneCycle=False):
        SaveState()
        controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
        if not controls:
            controls = {'deriv type':'analytic','min dM/M':0.0001,}     #fill in defaults if needed
        print 'Peak Fitting with '+controls['deriv type']+' derivatives:'
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        data = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        wx.BeginBusyCursor()
        try:
            G2pwd.DoPeakFit(FitPgm,peaks,background,limits,inst,data,oneCycle,controls)
        finally:
            wx.EndBusyCursor()    
        UpdatePeakGrid(G2frame,peaks)
        G2plt.PlotPatterns(G2frame)
        print 'finished'
        return
        
    def OnResetSigGam(event):
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to do!')
            return
        inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        Inst = dict(zip(inst[3],inst[1]))
        print len(Inst['Type'])
        for peak in peaks:
            if Inst['Type'] in ['PXC','PNC']:
                peak[4] = Inst['U']*tand(peak[0]/2.0)**2+Inst['V']*tand(peak[0]/2.0)+Inst['W']
                peak[6] = Inst['X']/cosd(peak[0]/2.0)+Inst['Y']*tand(peak[0]/2.0)
        UpdatePeakGrid(G2frame,peaks)
                
    def RefreshPeakGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        
        event.StopPropagation()
        data = G2frame.PeakTable.GetData()
        T = []
        for peak in data:T.append(peak[0])
        D = dict(zip(T,data))
        T.sort()
        X = []
        for key in T: X.append(D[key])
        data = X        
        
    def setBackgroundColors():
       for r in range(G2frame.dataDisplay.GetNumberRows()):
           for c in range(G2frame.dataDisplay.GetNumberCols()):
               if G2frame.dataDisplay.GetColLabelValue(c) in ['position','intensity','sigma','gamma']:
                   if float(G2frame.dataDisplay.GetCellValue(r,c)) < 0.:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.RED)
                   else:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.WHITE)
                                                  
    def KeyEditPeakGrid(event):
        rowList = G2frame.dataDisplay.GetSelectedRows()
        colList = G2frame.dataDisplay.GetSelectedCols()
        selectList = G2frame.dataDisplay.GetSelectedCells()
        data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif rowList:
            G2frame.dataDisplay.ClearSelection()
            if event.GetKeyCode() == wx.WXK_DELETE:
                G2frame.dataDisplay.ClearGrid()
                rowList.reverse()
                nDel = 0
                for row in rowList:
                    G2frame.PeakTable.DeleteRow(row)
                    nDel += 1
                if nDel:
                    msg = wg.GridTableMessage(G2frame.PeakTable, 
                        wg.GRIDTABLE_NOTIFY_ROWS_DELETED,0,nDel)
                    G2frame.dataDisplay.ProcessTableMessage(msg)
                data = G2frame.PeakTable.GetData()
                G2frame.PatternTree.SetItemPyData(G2frame.PickId,data[:-nDel])
                G2frame.dataDisplay.ForceRefresh()
                setBackgroundColors()
                if not len(G2frame.PatternTree.GetItemPyData(G2frame.PickId)): 
                    G2frame.dataFrame.PeakFit.Enable(False)
                        
        elif colList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.PeakTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data[row][col]=False
        elif selectList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for row,col in selectList:
                if G2frame.PeakTable.GetTypeName(row,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        data[row][col]=True
                    elif key == 78:  #'N'
                        data[row][col]=False
        G2plt.PlotPatterns(G2frame)
            
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.PeakMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_UNDO)
    G2frame.Bind(wx.EVT_MENU, OnLSQPeakFit, id=G2gd.wxID_LSQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnOneCycle, id=G2gd.wxID_LSQONECYCLE)
    G2frame.Bind(wx.EVT_MENU, OnClearPeaks, id=G2gd.wxID_CLEARPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnResetSigGam, id=G2gd.wxID_RESETSIGGAM)
    G2frame.dataFrame.PeakFit.Enable(False)
    if data:
        G2frame.dataFrame.PeakFit.Enable(True)
        G2frame.dataFrame.PFOneCycle.Enable(True)
    G2frame.PickTable = []
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
    G2frame.PatternTree.SetItemPyData(G2frame.PickId,data)
    G2frame.PeakTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Peak List')
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.PeakTable, True)
    setBackgroundColors()                         
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshPeakGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([535,350])
        
def UpdateBackground(G2frame,data):
    if len(data) < 2:       #add Debye diffuse & peaks scattering here
        data.append({'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]})
    if 'nPeaks' not in data[1]:
        data[1].update({'nPeaks':0,'peaksList':[]})
    ValObj = {}
            
    def OnBackCopy(event):
        histList = ['All',]+G2gd.GetPatternTreeDataNames(G2frame,['PWDR',])
        copyList = []
        dlg = wx.MultiChoiceDialog(G2frame, 
            'Copy parameters to which histograms?', 'Copy parameters', 
            histList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    copyList.append(histList[i])
                if 'All' in copyList: 
                    copyList = histList[1:]
            for item in copyList:
                Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Background'),
                    copy.copy(data))
        finally:
            dlg.Destroy()
        
    def BackSizer():
        
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
            G2frame.PatternTree.SetItemPyData(BackId,data)
            UpdateBackground(G2frame,data)
            
        def OnBakVal(event):
            Obj = event.GetEventObject()
            item = ValObj[Obj.GetId()][0]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = data[0][item]
            data[0][item] = value
            Obj.SetValue('%10.4f'%(value))
        
        backSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background function: '),0,wx.ALIGN_CENTER_VERTICAL)
        bakType = wx.ComboBox(G2frame.dataDisplay,value=data[0][0],
                choices=Choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakType.Bind(wx.EVT_COMBOBOX, OnNewType)
        topSizer.Add(bakType)
        topSizer.Add((5,0),0)
        bakRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
        bakRef.SetValue(bool(data[0][1]))
        bakRef.Bind(wx.EVT_CHECKBOX, OnBakRef)
        topSizer.Add(bakRef,0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. coeff.: '),0,wx.ALIGN_CENTER_VERTICAL)
        bakTerms = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[0][2]),choices=[str(i+1) for i in range(36)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakTerms.Bind(wx.EVT_COMBOBOX,OnBakTerms)
        topSizer.Add(bakTerms,0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add((5,0),0)
        backSizer.Add(topSizer)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background coefficients:'),0,wx.ALIGN_CENTER_VERTICAL)
        bakSizer = wx.FlexGridSizer(1,5,5,5)
        for i,value in enumerate(data[0][3:]):
            bakVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.4f'%(value),style=wx.TE_PROCESS_ENTER)
            bakSizer.Add(bakVal,0,wx.ALIGN_CENTER_VERTICAL)
            ValObj[bakVal.GetId()] = [i+3]
            bakVal.Bind(wx.EVT_TEXT_ENTER,OnBakVal)
            bakVal.Bind(wx.EVT_KILL_FOCUS,OnBakVal)
        backSizer.Add(bakSizer)
        return backSizer
        
    def DebyeSizer():
        
        def OnDebTerms(event):
            data[1]['nDebye'] = int(debTerms.GetValue())
            M = len(data[1]['debyeTerms'])
            N = data[1]['nDebye']
            if N > M:       #add terms
                for i in range(M,N): 
                    data[1]['debyeTerms'].append([1.0,False,1.0,False,0.010,False])
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(data[1]['debyeTerms'][-1])
            UpdateBackground(G2frame,data)

        def KeyEditPeakGrid(event):
            colList = debyeGrid.GetSelectedCols()
            if event.GetKeyCode() == wx.WXK_RETURN:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_CONTROL:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_SHIFT:
                event.Skip(True)
            elif colList:
                debyeGrid.ClearSelection()
                key = event.GetKeyCode()
                for col in colList:
                    if debyeTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                        if key == 89: #'Y'
                            for row in range(debyeGrid.GetNumberRows()): data[1]['debyeTerms'][row][col]=True
                        elif key == 78:  #'N'
                            for row in range(debyeGrid.GetNumberRows()): data[1]['debyeTerms'][row][col]=False

        
        debSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Debye scattering: '),0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. coeff.: '),0,wx.ALIGN_CENTER_VERTICAL)
        debTerms = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[1]['nDebye']),choices=[str(i) for i in range(12)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        debTerms.Bind(wx.EVT_COMBOBOX,OnDebTerms)
        topSizer.Add(debTerms,0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add((5,0),0)
        debSizer.Add(topSizer)
        if data[1]['nDebye']:
            debSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Debye diffuse terms:'),0,wx.ALIGN_CENTER_VERTICAL)       
            rowLabels = []
            for i in range(len(data[1]['debyeTerms'])): rowLabels.append(str(i))
            colLabels = ['A','refine','R','refine','U','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            debyeTable = G2gd.Table(data[1]['debyeTerms'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            debyeGrid = G2gd.GSGrid(parent=G2frame.dataDisplay)
            debyeGrid.SetTable(debyeTable, True)
            debyeGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            debyeGrid.AutoSizeColumns(False)
            debSizer.Add(debyeGrid)        
        return debSizer
      
    def PeaksSizer():

        def OnPeaks(event):
            data[1]['nPeaks'] = int(peaks.GetValue())
            M = len(data[1]['peaksList'])
            N = data[1]['nPeaks']
            if N > M:       #add terms
                for i in range(M,N): 
                    data[1]['peaksList'].append([1.0,False,1.0,False,0.10,False,0.10,False])
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(data[1]['peaksList'][-1])
            UpdateBackground(G2frame,data)

        def KeyEditPeakGrid(event):
            colList = peaksGrid.GetSelectedCols()
            if event.GetKeyCode() == wx.WXK_RETURN:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_CONTROL:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_SHIFT:
                event.Skip(True)
            elif colList:
                peaksGrid.ClearSelection()
                key = event.GetKeyCode()
                for col in colList:
                    if peaksTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                        if key == 89: #'Y'
                            for row in range(peaksGrid.GetNumberRows()): data[1]['peaksList'][row][col]=True
                        elif key == 78:  #'N'
                            for row in range(peaksGrid.GetNumberRows()): data[1]['peaksList'][row][col]=False

        peaksSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Peaks in background: '),0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. peaks: '),0,wx.ALIGN_CENTER_VERTICAL)
        peaks = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[1]['nPeaks']),choices=[str(i) for i in range(12)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        peaks.Bind(wx.EVT_COMBOBOX,OnPeaks)
        topSizer.Add(peaks,0,wx.ALIGN_CENTER_VERTICAL)
        topSizer.Add((5,0),0)
        peaksSizer.Add(topSizer)
        if data[1]['nPeaks']:
            peaksSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Peak list:'),0,wx.ALIGN_CENTER_VERTICAL)       
            rowLabels = []
            for i in range(len(data[1]['peaksList'])): rowLabels.append(str(i))
            colLabels = ['pos','refine','int','refine','sig','refine','gam','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            peaksTable = G2gd.Table(data[1]['peaksList'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            peaksGrid = G2gd.GSGrid(parent=G2frame.dataDisplay)
            peaksGrid.SetTable(peaksTable, True)
            peaksGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            peaksGrid.AutoSizeColumns(False)
            peaksSizer.Add(peaksGrid)        
        return peaksSizer
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BackMenu)
    G2frame.dataFrame.SetLabel('Background')
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnBackCopy,id=G2gd.wxID_BACKCOPY)
    BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Background')
    Choices = ['chebyschev','cosine','lin interpolate','inv interpolate','log interpolate']
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(BackSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(DebyeSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(PeaksSizer())
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
        
def UpdateLimitsGrid(G2frame, data):
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
        
    def RefreshLimitsGrid(event):
        event.StopPropagation()
        data = G2frame.LimitsTable.GetData()
        old = data[0]
        new = data[1]
        new[0] = max(old[0],new[0])
        new[1] = max(new[0],min(old[1],new[1]))
        data = [old,new]
        G2plt.PlotPatterns(G2frame)
        
    def OnLimitCopy(event):
        histList = ['All',]+G2gd.GetPatternTreeDataNames(G2frame,['PWDR',])
        copyList = []
        dlg = wx.MultiChoiceDialog(G2frame, 
            'Copy limits to which histograms?', 'Copy limits', 
            histList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    copyList.append(histList[i])
                if 'All' in copyList: 
                    copyList = histList[1:]
            for item in copyList:
                Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Limits'),
                    copy.copy(data))
        finally:
            dlg.Destroy()
        
    G2frame.LimitsTable = []
    colLabels = ['Tmin','Tmax']
    rowLabels = ['original','changed']
    Types = 2*[wg.GRID_VALUE_FLOAT+':10,3',]
    G2frame.LimitsTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Limits')
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.LimitMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnLimitCopy,id=G2gd.wxID_LIMITCOPY)
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.LimitsTable, True)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshLimitsGrid)                
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([230,160])
    
def UpdateInstrumentGrid(G2frame,data):
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
        return inst2data(inst,ref,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters')))        
    
    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
            if 'P' in insVal['Type']:                                       #update powder peak parameters
                for peak in peaks:
                    peak[4] = insVal['U']*tand(peak[0]/2.0)**2+insVal['V']*tand(peak[0]/2.0)+insVal['W']
                    peak[6] = insVal['X']/cosd(peak[0]/2.0)+insVal['Y']*tand(peak[0]/2.0)
                                                
    def OnReset(event):
        insVal.update(insDef)
        data = updateData(insVal,insRef)
        RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
        UpdateInstrumentGrid(G2frame,data)
        
    def OnInstCopy(event):
        histList = ['All',]+G2gd.GetPatternTreeDataNames(G2frame,['PWDR',])
        copyList = []
        dlg = wx.MultiChoiceDialog(G2frame, 
            'Copy parameters to which histograms?', 'Copy parameters', 
            histList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    copyList.append(histList[i])
                if 'All' in copyList: 
                    copyList = histList[1:]
            for item in copyList:
                Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                instData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))
                if len(data) == len(instData):                          #don't mix lam & lam1/lam2 parms!
                    for i,item in enumerate(data[1:]):                  #skip default values in tuple
                        instData[i+1][:-1] = copy.copy(item[:-1])       #skip azimuth at end
                else:
                    print item+' not copied - instrument parameters not commensurate'
        finally:
            dlg.Destroy()
        
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
        UpdateInstrumentGrid(G2frame,data)
                
    def OnNewType(event):
        insVal['Type'] = typePick.GetValue()
        data = updateData(insVal,insRef)
        UpdateInstrumentGrid(G2frame,data)
        
    def OnLamPick(event):
        lamType = lamPick.GetValue()
        insVal['Lam1'] = waves[lamType][0]
        insVal['Lam2'] = waves[lamType][1]
        data = updateData(insVal,insRef)
        UpdateInstrumentGrid(G2frame,data)
                 
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
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    histoName = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[-1]
    ifHisto = IsHistogramInAnyPhase(G2frame,histoName)
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
    G2frame.dataFrame.SetLabel('Instrument Parameters')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    instSizer = wx.FlexGridSizer(2,6,5,5)
    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Histogram Type:'),0,wx.ALIGN_CENTER_VERTICAL)
    if 'P' in insVal['Type']:                   #powder data
        G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.InstMenu)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()
        G2frame.Bind(wx.EVT_MENU,OnReset,id=G2gd.wxID_INSTPRMRESET)
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2gd.wxID_INSTCOPY)
        G2frame.Bind(wx.EVT_MENU,OnWaveChange,id=G2gd.wxID_CHANGEWAVETYPE)        
        typePick = wx.ComboBox(G2frame.dataDisplay,value=insVal['Type'],
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
            instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Azimuth: %7.2f'%(insVal['Azimuth'])),0,wx.ALIGN_CENTER_VERTICAL)
            if 'Lam1' in insVal:
                instSizer.Add((5,5),0)
                instSizer.Add((5,5),0)
                instSizer.Add((5,5),0)
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Ka1/Ka2:'),
                        0,wx.ALIGN_CENTER_VERTICAL)
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'%8.6f/%8.6f'%(insVal['Lam1'],insVal['Lam2'])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'Select:'),0,wx.ALIGN_CENTER_VERTICAL)
                choice = ['TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa']
                lamPick = wx.ComboBox(G2frame.dataDisplay,value=' ',choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0)
                instSizer.Add(waveSizer,0)
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' I(L2)/I(L1): (%10.4f)'%(insDef['I(L2)/I(L1)'])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                ratVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.4f'%(insVal['I(L2)/I(L1)']),style=wx.TE_PROCESS_ENTER)
                ratVal.Bind(wx.EVT_TEXT_ENTER,OnRatValue)
                ratVal.Bind(wx.EVT_KILL_FOCUS,OnRatValue)
                instSizer.Add(ratVal,0)
                ratRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
                ratRef.SetValue(bool(insRef['I(L2)/I(L1)']))
                ratRef.Bind(wx.EVT_CHECKBOX, OnRatRef)
                instSizer.Add(ratRef,0,wx.ALIGN_CENTER_VERTICAL)
                
            else:
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Lam: (%10.6f)'%(insDef['Lam'])),
                    0,wx.ALIGN_CENTER_VERTICAL)
                waveVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.6f'%(insVal['Lam']),style=wx.TE_PROCESS_ENTER)
                waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveValue)
                waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveValue)
                instSizer.Add(waveVal,0,wx.ALIGN_CENTER_VERTICAL)
                if ifHisto:
                    waveRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
                    waveRef.SetValue(bool(insRef['Lam']))
                    waveRef.Bind(wx.EVT_CHECKBOX, OnWaveRef)
                    instSizer.Add(waveRef,0,wx.ALIGN_CENTER_VERTICAL)
                else:
                    instSizer.Add((5,5),0)
            for item in ['Zero','Polariz.']:
                fmt = '%10.4f'
                Fmt = ' %s: ('+fmt+')'
                if item in insDef:
                    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,Fmt%(item,insDef[item])),
                            0,wx.ALIGN_CENTER_VERTICAL)
                    itemVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                    ValObj[itemVal.GetId()] = [item,fmt]
                    itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                    itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                    instSizer.Add(itemVal,0,wx.ALIGN_CENTER_VERTICAL)
                    if ifHisto:
                        itemRef = wx.CheckBox(G2frame.dataDisplay,wx.ID_ANY,label=' Refine?')
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
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,Fmt%(item,insDef[item])),
                        0,wx.ALIGN_CENTER_VERTICAL)
                itemVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                ValObj[itemVal.GetId()] = [item,fmt]
                itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                instSizer.Add(itemVal,0,wx.ALIGN_CENTER_VERTICAL)
                itemRef = wx.CheckBox(G2frame.dataDisplay,wx.ID_ANY,label=' Refine?')
                itemRef.SetValue(bool(insRef[item]))
                RefObj[itemRef.GetId()] = item
                itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                instSizer.Add(itemRef,0,wx.ALIGN_CENTER_VERTICAL)
        else:                                   #time of flight (neutrons)
            pass                                #for now
        
        

    else:                       #single crystal data
        typePick = wx.ComboBox(G2frame.dataDisplay,value=insVal['Type'],
            choices=['SXC','SNC','SNT'],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        typePick.Bind(wx.EVT_COMBOBOX, OnNewType)
        instSizer.Add(typePick,0,wx.ALIGN_CENTER_VERTICAL)
        if 'C' in insVal['Type']:               #constant wavelength
            instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Lam: %10.6f'%(insDef['Lam'])),
                    0,wx.ALIGN_CENTER_VERTICAL)
        else:                                   #time of flight (neutrons)
            pass                                #for now
        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(instSizer,0)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
    
def UpdateSampleGrid(G2frame,data):
    
    def OnSampleCopy(event):
        histName = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        copyNames = ['Scale',]
        dataType = data['Type']
        histType = 'HKLF'
        if 'PWDR' in histName:          #else HKLF - only Scale
            histType = 'PWDR'
            if 'Debye' in dataType:
                copyNames += ['DisplaceX','DisplaceY','Absorption']
            else:       #Bragg-Brentano
                copyNames += ['Shift','Transparency']
        copyNames += ['Omega','Chi','Phi']
        copyDict = {}
        for parm in copyNames:
            copyDict[parm] = data[parm]
        histList = ['All '+histType,]
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while item:
            name = G2frame.PatternTree.GetItemText(item)
            if histType in name and name != histName:
                histList.append(name)
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        if len(histList) == 1:      #nothing to copy to!
            return
        copyList = []
        dlg = wx.MultiChoiceDialog(G2frame,'Copy parameters from\n'+histName,
            'Copy parameters to which histograms?',histList,wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    copyList.append(histList[i])
                if 'All '+histType in copyList: 
                    copyList = histList[1:]
            for item in copyList:
                Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                sampleData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
                sampleData.update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.SampleMenu)
    G2frame.dataFrame.SetLabel('Sample Parameters')
    G2frame.Bind(wx.EVT_MENU, OnSampleCopy, id=G2gd.wxID_SAMPLECOPY)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)

#patch
    if not 'Gonio. radius' in data:
        data['Gonio. radius'] = 200.0
    if not 'Omega' in data:
        data.update({'Omega':0.0,'Chi':0.0,'Phi':0.0})
#patch end
    
    parms = [['Gonio. radius',' Goniometer radius(mm): ','%.2f',]]
    if data['Type'] == 'Debye-Scherrer':
        parms += [['DisplaceX',u' Sample X displacement(\xb5m): ','%.2f',],
            ['DisplaceY',u' Sample Y displacement(\xb5m): ','%.2f',],
            ['Absorption',u' Sample absorption(\xb5r): ','%.4f',],]
    elif data['Type'] == 'Bragg-Brentano':
        parms += [['Shift',u' Sample displacement(\xb5m): ','%.2f',],
            ['Transparency',u' Sample transparency(1/\xb5eff,cm): ','%.4f'],]
    parms.append(['Omega','Goniometer omega:','%.2f'])
    parms.append(['Chi','Goniometer chi:','%.2f'])
    parms.append(['Phi','Goniometer phi:','%.2f'])
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
        G2frame.dataDisplay.Destroy()
        UpdateSampleGrid(G2frame,data)
        
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
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample parameters: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    parmSizer = wx.FlexGridSizer(9,2,5,0)
    scaleRef = wx.CheckBox(G2frame.dataDisplay,label=' Histogram scale factor: ')
    scaleRef.SetValue(data['Scale'][1])
    scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
    parmSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
    scaleVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,
        '%.4f'%(data['Scale'][0]),style=wx.TE_PROCESS_ENTER)
    scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
    scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
    parmSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
    typeSizer = wx.BoxSizer(wx.HORIZONTAL)
    choices = ['Debye-Scherrer','Bragg-Brentano',]
    histoType = wx.ComboBox(G2frame.dataDisplay,wx.ID_ANY,value=data['Type'],choices=choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    histoType.Bind(wx.EVT_COMBOBOX, OnHistoType)
    parmSizer.Add(histoType)
    parmSizer.Add((5,5),0)
    
    for parm in parms:
        if 'list' in str(type(data[parm[0]])):
            parmRef = wx.CheckBox(G2frame.dataDisplay,label=parm[1])
            objList[parmRef.GetId()] = parm[0]
            parmRef.SetValue(data[parm[0]][1])
            parmRef.Bind(wx.EVT_CHECKBOX, OnParmRef)
            parmSizer.Add(parmRef,0,wx.ALIGN_CENTER_VERTICAL)
            parmVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,
                parm[2]%(data[parm[0]][0]),style=wx.TE_PROCESS_ENTER)
        else:
            parmSizer.Add(wx.StaticText(G2frame.dataDisplay,label=parm[1]),
                0,wx.ALIGN_CENTER_VERTICAL)
            parmVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,
                parm[2]%(data[parm[0]]),style=wx.TE_PROCESS_ENTER)        
        objList[parmVal.GetId()] = parm
        parmVal.Bind(wx.EVT_TEXT_ENTER,OnParmVal)
        parmVal.Bind(wx.EVT_KILL_FOCUS,OnParmVal)
        parmSizer.Add(parmVal,1,wx.EXPAND)
    mainSizer.Add(parmSizer)
    mainSizer.Add((0,5),0)    
    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)
                
def UpdateIndexPeaksGrid(G2frame, data):
    IndexId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List')
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
    Inst = dict(zip(inst[3],inst[1]))
    try:
        wave = Inst['Lam']
    except KeyError:
        wave = Inst['Lam1']
    
    def RefreshIndexPeaksGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        data = G2frame.IndexPeaksTable.GetData()
        if c == 2:
            if data[r][c]:
                data[r][c] = False
            else:
                data[r][c] = True
            G2frame.IndexPeaksTable.SetData(data)
            G2frame.PatternTree.SetItemPyData(IndexId,data)
            G2frame.dataDisplay.ForceRefresh()
            
    def OnReload(event):
        data = []
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
        for peak in peaks:
            dsp = wave/(2.0*sind((peak[0]-Inst['Zero'])/2.0))
            data.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])
        G2frame.PatternTree.SetItemPyData(IndexId,data)
        UpdateIndexPeaksGrid(G2frame,data)
        
    def KeyEditPickGrid(event):
        colList = G2frame.dataDisplay.GetSelectedCols()
        rowList = G2frame.dataDisplay.GetSelectedRows()
        data = G2frame.PatternTree.GetItemPyData(IndexId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif colList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.IndexPeaksTable.GetColLabelValue(col) in ['use','refine']:
                    if key == 89: #'Y'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[row][col]=False
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    if 'PWD' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
        G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.IndPeaksMenu)
        G2frame.Bind(wx.EVT_MENU, OnReload, id=G2gd.wxID_INDXRELOAD)
    G2frame.dataFrame.IndexPeaks.Enable(False)
    G2frame.IndexPeaksTable = []
    if data:
        G2frame.dataFrame.IndexPeaks.Enable(True)
        cells = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List'))
        if cells:
            cellist = cells[2]
            dmin = cells[3]
            G2frame.HKL = []
            for i,cell in enumerate(cellist):
                if cell[-1]:
                    ibrav = cell[2]
                    A = G2lat.cell2A(cell[3:9])
                    G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                    G2indx.IndexPeaks(data,G2frame.HKL)
                    for hkl in G2frame.HKL:
                        hkl.append(2.0*asind(wave/(2.*hkl[3]))+Inst['Zero'])             
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','intensity','use','indexed','h','k','l','d-obs','d-calc']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_BOOL,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5']
    G2frame.PatternTree.SetItemPyData(IndexId,data)
    G2frame.IndexPeaksTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Index Peak List')
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)                
    G2frame.dataDisplay.SetTable(G2frame.IndexPeaksTable, True)
    for r in range(G2frame.dataDisplay.GetNumberRows()):
        for c in range(G2frame.dataDisplay.GetNumberCols()):
            if c == 2:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, RefreshIndexPeaksGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)                 
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([490,300])
  
def UpdateUnitCellsGrid(G2frame, data):
    UnitCellsId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    spaceGroups = ['F m 3 m','I m 3 m','P m 3 m','R -3 H','P 6/m m m','I 4/m m m',
        'P 4/m m m','F m m m','I m m m','C m m m','P m m m','C 2/m','P 2/m','P -1']
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
    Inst = dict(zip(inst[3],inst[1]))
    if 'Lam' in Inst:
        wave = Inst['Lam']
    else:
        wave = Inst['Lam1']
        
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
            Zero = min(5.0,max(-5.0,float(zero.GetValue())))
        except ValueError:
            Zero = 0.0
        controls[1] = Zero
        zero.SetValue("%.4f"%(Zero))
        
    def OnZeroVar(event):
        controls[0] = zeroVar.GetValue()
        
    def OnBravSel(event):
        controls[5] = bravSel.GetString(bravSel.GetSelection())       
        UpdateUnitCellsGrid(G2frame,data)
        
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
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId    
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        dmin = wave/(2.0*sind(limits[1]/2.0))
        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
        for hkl in G2frame.HKL:
            hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'])             
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
            
    def OnSortCells(event):
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        c =  event.GetCol()
        if colLabels[c] == 'M20':
            cells = G2indx.sortM20(cells)
        elif colLabels[c] == 'Volume':
            cells = G2indx.sortVolume(cells)
        else:
            return
        data = [controls,bravais,cells,dmin]
        G2frame.PatternTree.SetItemPyData(UnitCellsId,data)
        UpdateUnitCellsGrid(G2frame,data)
        
    def CopyUnitCell(event):
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        for Cell in cells:
            if Cell[-1]:
                break
        cell = Cell[2:9]
        controls[4] = 1
        controls[5] = bravaisSymb[cell[0]]
        controls[6:12] = cell[1:8]
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        G2frame.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin])
        UpdateUnitCellsGrid(G2frame,data)
        
        G2frame.dataFrame.RefineCell.Enable(True)
                
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
             
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId    
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!', 'Nothing to refine!')
            return        
        print 'Refine cell'
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        dmin = G2indx.getDmin(peaks)-0.005
        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
        G2indx.IndexPeaks(peaks,G2frame.HKL)
        Lhkl,M20,X20,Aref,Zero = G2indx.refinePeaksZ(peaks,wave,ibrav,A,controls[1],controls[0])            
        controls[1] = Zero
        controls[6:12] = G2lat.A2cell(Aref)
        controls[12] = G2lat.calc_V(Aref)
        data = [controls,bravais,cells,dmin]
        cells = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2]
        for cell in cells:
            cell[-1] = False
        cells.insert(0,[M20,X20,ibrav]+controls[6:13]+[True,])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
        UpdateUnitCellsGrid(G2frame,data)
        print "%s%10.3f" % ('refinement M20 = ',M20)
        print 'unindexed lines = ',X20
        cellPrint(ibrav,Aref)
        for hkl in G2frame.HKL:
            hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'])             
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
        
    def IndexPeaks(event):
        PatternId = G2frame.PatternId    
        print 'Peak Indexing'
        try:
            controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
            cells = []
        except ValueError:
            G2frame.ErrorDialog('Error','Need to set controls in Unit Cell List first')
            return
        if True not in bravais:
            G2frame.ErrorDialog('Error','No Bravais lattices selected')
            return
        G2frame.dataFrame.CopyCell.Enable(False)
        G2frame.dataFrame.RefineCell.Enable(False)
        OK,dmin,cells = G2indx.DoIndexPeaks(peaks,inst[1],controls,bravais)
        if OK:
            data = [controls,bravais,cells,dmin]
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
            UpdateUnitCellsGrid(G2frame,data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                G2frame.HKL = G2lat.GenHBravais(dmin,bestCell[2],G2lat.cell2A(bestCell[3:9]))
                for hkl in G2frame.HKL:
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'])             
                if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame)
                else:
                    G2plt.PlotPatterns(G2frame)
            G2frame.dataFrame.CopyCell.Enable(True)
            G2frame.dataFrame.IndexPeaks.Enable(True)
            G2frame.dataFrame.MakeNewPhase.Enable(True)
            UpdateUnitCellsGrid(G2frame,data)
                
    def RefreshUnitCellsGrid(event):
        cells,dmin = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2:]
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
                G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                for hkl in G2frame.HKL:
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'])             
                if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame)
                else:
                    G2plt.PlotPatterns(G2frame)
        
    def MakeNewPhase(event):
        if not G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases'):
            sub = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PhaseName = dlg.GetValue()
                cells = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2]
                for Cell in cells:
                    if Cell[-1]:
                        break
                cell = Cell[2:10]        
                sub = G2frame.PatternTree.AppendItem(parent=sub,text=PhaseName)
                E,SGData = G2spc.SpcGroup(spaceGroups[cell[0]])
                G2frame.PatternTree.SetItemPyData(sub, \
                    {'General':{'Name':PhaseName,'Type':'nuclear','SGData':SGData,
                    'Cell':[False,]+cell[1:],
                    'Pawley dmin':1.00},'Atoms':[],'Drawing':{},'Histograms':{}})
                Status.SetStatusText('Change space group if needed')
        finally:
            dlg.Destroy()
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.IndexMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU, IndexPeaks, id=G2gd.wxID_INDEXPEAKS)
    G2frame.Bind(wx.EVT_MENU, CopyUnitCell, id=G2gd.wxID_COPYCELL)
    G2frame.Bind(wx.EVT_MENU, RefineCell, id=G2gd.wxID_REFINECELL)
    G2frame.Bind(wx.EVT_MENU, MakeNewPhase, id=G2gd.wxID_MAKENEWPHASE)
    
    controls,bravais,cells,dmin = data
    if len(controls) < 13:              #add cell volume if missing
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    G2frame.PatternTree.SetItemPyData(UnitCellsId,data)            #update with volume
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
    
    G2frame.dataFrame.SetLabel('Unit Cells List')
    G2frame.sp = wx.SplitterWindow(G2frame.dataFrame)
    G2frame.dataDisplay = wx.Panel(G2frame.sp, style=wx.SUNKEN_BORDER)
    G2frame.dataFrame.IndexPeaks.Enable(False)
    peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
    if peaks:
        G2frame.dataFrame.IndexPeaks.Enable(True)
    G2frame.dataFrame.RefineCell.Enable(False)
    if controls[12] > 1.0:                               #if a "real" volume (i.e. not default)
        G2frame.dataFrame.RefineCell.Enable(True)    
    G2frame.dataFrame.CopyCell.Enable(False)
    G2frame.dataFrame.MakeNewPhase.Enable(False)        
    if cells:
        G2frame.bottom = wx.Panel(G2frame.sp, style=wx.SUNKEN_BORDER)
        G2frame.sp.SplitHorizontally(G2frame.dataDisplay,G2frame.bottom,0)
        G2frame.dataFrame.CopyCell.Enable(True)
        G2frame.dataFrame.MakeNewPhase.Enable(True)        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Indexing controls: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,5,5,5)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Max Nc/Nobs '),0,wx.ALIGN_CENTER_VERTICAL)
    NcNo = wx.SpinCtrl(G2frame.dataDisplay)
    NcNo.SetRange(1,6)
    NcNo.SetValue(controls[2])
    NcNo.Bind(wx.EVT_SPINCTRL,OnNcNo)
    littleSizer.Add(NcNo,0,wx.ALIGN_CENTER_VERTICAL)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Start Volume '),0,wx.ALIGN_CENTER_VERTICAL)
    startVol = wx.TextCtrl(G2frame.dataDisplay,value=str('%d'%(controls[3])),style=wx.TE_PROCESS_ENTER)
    startVol.Bind(wx.EVT_TEXT_ENTER,OnStartVol)
    startVol.Bind(wx.EVT_KILL_FOCUS,OnStartVol)
    littleSizer.Add(startVol,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Select Bravais Lattices for indexing: '),
        0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,7,5,5)
    bravList = []
    bravs = zip(bravais,bravaisNames)
    for brav,bravName in bravs:
        bravCk = wx.CheckBox(G2frame.dataDisplay,label=bravName)
        bravList.append(bravCk.GetId())
        bravCk.SetValue(brav)
        bravCk.Bind(wx.EVT_CHECKBOX,OnBravais)
        littleSizer.Add(bravCk,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Cell Refinement: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Bravais lattice "),0,wx.ALIGN_CENTER_VERTICAL)
    bravSel = wx.Choice(G2frame.dataDisplay,choices=bravaisSymb)
    bravSel.SetSelection(bravaisSymb.index(controls[5]))
    bravSel.Bind(wx.EVT_CHOICE,OnBravSel)
    littleSizer.Add(bravSel,0,wx.ALIGN_CENTER_VERTICAL)
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Zero offset"),0,wx.ALIGN_CENTER_VERTICAL)
    zero = wx.TextCtrl(G2frame.dataDisplay,value="%.4f"%(controls[1]),style=wx.TE_PROCESS_ENTER)
    zero.Bind(wx.EVT_TEXT_ENTER,OnZero)
    zero.Bind(wx.EVT_KILL_FOCUS,OnZero)
    littleSizer.Add(zero,0,wx.ALIGN_CENTER_VERTICAL)
    zeroVar = wx.CheckBox(G2frame.dataDisplay,label="Refine?")
    zeroVar.SetValue(controls[0])
    zeroVar.Bind(wx.EVT_CHECKBOX,OnZeroVar)
    littleSizer.Add(zeroVar,0,wx.ALIGN_CENTER_VERTICAL)
    hklShow = wx.Button(G2frame.dataDisplay,label="Show refined hkl positions + Zero offset")
    hklShow.Bind(wx.EVT_BUTTON,OnHklShow)
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
        littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=txt),0,wx.ALIGN_CENTER_VERTICAL)
        if ifEdit:          #a,b,c,etc.
            cellVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[6+Id])),style=wx.TE_PROCESS_ENTER)
            cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
            cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
            littleSizer.Add(cellVal,0,wx.ALIGN_CENTER_VERTICAL)
            cellList.append(cellVal.GetId())
        else:               #volume
            volVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[12])),style=wx.TE_READONLY)
            volVal.SetBackgroundColour(VERY_LIGHT_GREY)
            littleSizer.Add(volVal,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    topSize = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataDisplay.SetSize(topSize)
    if cells:
        if ibrav == 13:
            topSize[1] += 230
        else:
            topSize[1] += 200
    G2frame.dataFrame.setSizePosLeft(topSize)    
    
    if cells:
        bottomSize = G2frame.bottom.GetSize()
        if ibrav == 13:
            bottomSize[1] -= 240
        else:
            bottomSize[1] -= 210
        wx.StaticText(parent=G2frame.bottom,label=' Indexing Result ')
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
                G2frame.HKL = G2lat.GenHBravais(dmin,cell[2],A)
                for hkl in G2frame.HKL:
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'])             
            table.append(row)
        UnitCellsTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        gridDisplay = G2gd.GSGrid(G2frame.bottom)
        gridDisplay.SetPosition(wx.Point(0,20))                
        gridDisplay.SetTable(UnitCellsTable, True)
        G2frame.dataFrame.CopyCell.Enable(True)
        gridDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,RefreshUnitCellsGrid)
        gridDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnSortCells)
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

def UpdateReflectionGrid(G2frame,data):
    if not data:
        print 'No phases, no reflections'
        return
    phases = data.keys()
    
    def OnSelectPhase(event):
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Phase',phases)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.RefList = phases[sel]
                UpdateReflectionGrid(G2frame,data)
        finally:
            dlg.Destroy()
        G2plt.PlotPatterns(G2frame)
        
        
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.ReflMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    G2frame.Bind(wx.EVT_MENU, OnSelectPhase, id=G2gd.wxID_SELECTPHASE)
    G2frame.dataFrame.SelectPhase.Enable(False)
    if len(data) > 1:
        G2frame.dataFrame.SelectPhase.Enable(True)
    rowLabels = []
    refList = [refl[:11] for refl in data[G2frame.RefList]]
    for i in range(len(refList)): rowLabels.append(str(i))
    colLabels = ['H','K','L','mul','d','pos','sig','gam','Fosq','Fcsq','phase',]
    Types = 4*[wg.GRID_VALUE_LONG,]+4*[wg.GRID_VALUE_FLOAT+':10,4',]+2*[wg.GRID_VALUE_FLOAT+':10,2',]+[wg.GRID_VALUE_FLOAT+':10,3',]
    G2frame.PeakTable = G2gd.Table(refList,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Reflection List for '+G2frame.RefList)
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.PeakTable, True)
    G2frame.dataDisplay.EnableEditing(False)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([555,350])

def UpdatePDFGrid(G2frame,data):
    global inst
    tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
    dataFile = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    powName = 'PWDR'+dataFile[4:]
    powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, powName)
    fullLimits,limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Limits'))
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Instrument Parameters'))
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
            UpdatePDFGrid(G2frame,data)
        
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
            G2plt.PlotISFG(G2frame,newPlot=True)
                        
        item = data[key]
        fileList = np.array(GetFileList('PWDR')).T[1]
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' '+key+' file:'),0,wx.ALIGN_CENTER_VERTICAL)
        fileName = wx.ComboBox(G2frame.dataDisplay,value=item['Name'],choices=fileList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        itemDict[fileName.GetId()] = [key,'Name','%s']
        fileName.Bind(wx.EVT_COMBOBOX,OnSelectFile)        
        fileSizer.Add(fileName,0,)
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Multiplier:'),0,wx.ALIGN_CENTER_VERTICAL)
        mult = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(item['Mult']),style=wx.TE_PROCESS_ENTER)
        itemDict[mult.GetId()] = [key,'Mult','%.3f']
        mult.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
        mult.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
        fileSizer.Add(mult,0,)
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Add:'),0,wx.ALIGN_CENTER_VERTICAL)
        add = wx.TextCtrl(G2frame.dataDisplay,value='%.0f'%(item['Add']),style=wx.TE_PROCESS_ENTER)
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
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
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
            UpdatePDFGrid(G2frame,data)
            auxPlot = ComputePDF(data)
            G2plt.PlotISFG(G2frame,newPlot=True)        
        
        elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
            label=' Element: '+'%2s'%(ElData['Symbol'])+' * '),0,wx.ALIGN_CENTER_VERTICAL)
        num = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(ElData['FormulaNo']),style=wx.TE_PROCESS_ENTER)
        num.Bind(wx.EVT_TEXT_ENTER,OnFractionChange)        
        num.Bind(wx.EVT_KILL_FOCUS,OnFractionChange)
        elemSizer.Add(num,0,wx.ALIGN_CENTER_VERTICAL)
        elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
            label="f': %.3f"%(ElData['fp'])+' f": %.3f'%(ElData['fpp'])+' mu: %.2f barns'%(ElData['mu']) ),
            0,wx.ALIGN_CENTER_VERTICAL)
            
    def OnGeometry(event):
        data['Geometry'] = geometry.GetValue()
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def OnDetType(event):
        data['DetType'] = detType.GetValue()
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def OnFormVol(event):
        try:
            value = float(formVol.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Form Vol']
        data['Form Vol'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
        
    def OnDiameter(event):
        try:
            value = float(diam.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Diam']
        data['Diam'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnPolaVal(event):
        try:
            value = float(polaVal.GetValue())
            if not (0.0 <= value <= 1.0):
                raise ValueError
        except ValueError:
            value = inst['Polariz.']
        inst['Polariz.'] = value
        polaVal.SetValue('%.2f'%(inst['Polariz.']))
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
                
    def OnAzimVal(event):
        try:
            value = float(azimVal.GetValue())
            if not (0. <= value <= 360.):
                raise ValueError
        except ValueError:
            value = inst['Azimuth']
        inst['Azimuth'] = value
        azimVal.SetValue('%.1f'%(inst['Azimuth']))
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
                        
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
        G2plt.PlotISFG(G2frame,newPlot=False)
        
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
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnRulSlider(event):
        value = int(rulandSldr.GetValue())/1000.
        data['Ruland'] = max(0.001,value)
        rulandWdt.SetValue('%.3f'%(data['Ruland']))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnLorch(event):
        data['Lorch'] = lorch.GetValue()
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
                        
    def OnPacking(event):
        try:
            value = float(pack.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Pack']
        data['Pack'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
                
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
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
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
        G2plt.PlotISFG(G2frame,newPlot=True)
        
    def OnResetQ(event):
        resetQ.SetValue(False)
        data['QScaleLim'][1] = qLimits[1]
        SQmax.SetValue('%.1f'%(data['QScaleLim'][1]))
        data['QScaleLim'][0] = 0.9*qLimits[1]
        SQmin.SetValue('%.1f'%(data['QScaleLim'][0]))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        

    def GetFileList(fileType,skip=None):
        fileList = [[False,'',0]]
        Source = ''
        id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while id:
            name = G2frame.PatternTree.GetItemText(id)
            if fileType in name:
                if id == skip:
                    Source = name
                else:
                    fileList.append([False,name,id])
            id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
        
    def OnCopyPDFControls(event):
        import copy
        TextList,Source = GetFileList('PDF',skip=G2frame.PatternId)
        TextList[0] = [False,'All PDF',0]
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "PDF" pattern')
            return
        dlg = G2frame.CopyDialog(G2frame,'Copy PDF controls','Copy controls from '+Source+' to:',TextList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                if result[0][0]:
                    result = TextList[1:]
                    for item in result: item[0] = True
                for i,item in enumerate(result):
                    ifcopy,name,id = item
                    if ifcopy:
                        olddata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'))
                        sample = olddata['Sample']
                        olddata.update(copy.deepcopy(data))
                        olddata['Sample'] = sample
                        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'),olddata)
                Status.SetStatusText('PDF controls copied')
        finally:
            dlg.Destroy()
                
    def OnSavePDFControls(event):
        print 'save PDF controls?'
        
    def OnLoadPDFControls(event):
        print 'Load PDF controls?'
        
    def OnAddElement(event):
        ElList = data['ElList']
        PE = G2elem.PickElement(G2frame,oneOnly=True)
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
        UpdatePDFGrid(G2frame,data)
        
    def OnDeleteElement(event):
        ElList = data['ElList']
        choice = ElList.keys()
        dlg = G2elem.DeleteElement(G2frame,choice=choice)
        if dlg.ShowModal() == wx.ID_OK:
            del ElList[dlg.GetDeleteElement()]
        dlg.Destroy()
        UpdatePDFGrid(G2frame,data)
                
    def ComputePDF(Data):
        xydata = {}
        for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
            name = Data[key]['Name']
            if name:
                xydata[key] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name))
                PDFname = name
        powName = xydata['Sample'][2]
        powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,powName)
        inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Instrument Parameters'))
        inst = dict(zip(inst[3],inst[1]))
        auxPlot = G2pwd.CalcPDF(Data,inst,xydata)
        PDFId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'PDF '+powName[4:])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'I(Q)'+powName[4:]),xydata['IofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'S(Q)'+powName[4:]),xydata['SofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'F(Q)'+powName[4:]),xydata['FofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'G(R)'+powName[4:]),xydata['GofR'])
        return auxPlot
        
    def OnComputePDF(event):
        print 'Calculating PDF:'
        auxPlot = ComputePDF(data)
        print 'Done calculating PDF:'
        Status.SetStatusText('PDF computed')
        for plot in auxPlot:
            G2plt.PlotXY(G2frame,plot[:2],type=plot[2])
        
        G2plt.PlotISFG(G2frame,newPlot=True,type='I(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='S(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='F(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='G(R)')
        
    def OnComputeAllPDF(event):
        print 'Calculating PDFs:'
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                Name = G2frame.PatternTree.GetItemText(id)
                if 'PDF' in Name:
                    Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls'))
                    auxPlot = ComputePDF(Data)                    
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            Status.SetStatusText('All PDFs computed')
            G2plt.PlotISFG(G2frame,newPlot=True,type='G(R)')
            print ' Done calculating PDFs:'
        
    def OnShowTip(G2frame,tip):
        print tip    
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.PDFMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyPDFControls, id=G2gd.wxID_PDFCOPYCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSavePDFControls, id=G2gd.wxID_PDFSAVECONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadPDFControls, id=G2gd.wxID_PDFLOADCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddElement, id=G2gd.wxID_PDFADDELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2gd.wxID_PDFDELELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputePDF, id=G2gd.wxID_PDFCOMPUTE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputeAllPDF, id=G2gd.wxID_PDFCOMPUTEALL)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' PDF data files: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    str = ' Sample file: PWDR %s   Wavelength, A: %.5f  Energy, keV: %.3f  Polariz.: %.2f '%(dataFile[3:],wave,keV,polariz)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=str),0,wx.ALIGN_CENTER_VERTICAL)
#    dataSizer = wx.BoxSizer(wx.HORIZONTAL)
#    dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Azimuth'),0,wx.ALIGN_CENTER_VERTICAL)
#    azimVal = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(inst['Azimuth']))
#    azimVal.Bind(wx.EVT_TEXT_ENTER,OnAzimVal)        
#    azimVal.Bind(wx.EVT_KILL_FOCUS,OnAzimVal)
#    dataSizer.Add(azimVal,0)    
#    dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Polarization'),0,wx.ALIGN_CENTER_VERTICAL)
#    polaVal = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(inst['Polariz.']))
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
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample information: '),0,wx.ALIGN_CENTER_VERTICAL)
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
    midSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Formula volume: '),0,wx.ALIGN_CENTER_VERTICAL)
    formVol = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(data['Form Vol']))
    formVol.Bind(wx.EVT_TEXT_ENTER,OnFormVol)        
    formVol.Bind(wx.EVT_KILL_FOCUS,OnFormVol)
    midSizer.Add(formVol,0)
    midSizer.Add(wx.StaticText(G2frame.dataDisplay,
        label=' Theoretical absorption: %.4f cm-1 Sample absorption: %.4f cm-1'%(Abs,Abs*data['Pack'])),
        0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(midSizer,0)
    mainSizer.Add((5,5),0)    

    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample geometry: '),0,wx.ALIGN_CENTER_VERTICAL)
    choice = ['Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate']
    geometry = wx.ComboBox(G2frame.dataDisplay,value=data['Geometry'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    geometry.Bind(wx.EVT_COMBOBOX, OnGeometry)
    geoBox.Add(geometry,0)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample diameter/thickness, mm: '),0,wx.ALIGN_CENTER_VERTICAL)
    diam = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['Diam']))
    diam.Bind(wx.EVT_TEXT_ENTER,OnDiameter)        
    diam.Bind(wx.EVT_KILL_FOCUS,OnDiameter)
#    diam.Bind(wx.EVT_SET_FOCUS,OnShowTip(G2frame,'tip')) #this doesn't work - what would????
    geoBox.Add(diam,0)
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Packing: '),0,wx.ALIGN_CENTER_VERTICAL)
    pack = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(data['Pack']))
    pack.Bind(wx.EVT_TEXT_ENTER,OnPacking)        
    pack.Bind(wx.EVT_KILL_FOCUS,OnPacking)
    geoBox.Add(pack,0)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample transmission: %.3f %%'%(Trans)),0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
        
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' S(Q)->F(Q)->G(R) controls: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Detector type: '),0,wx.ALIGN_CENTER_VERTICAL)
    choice = ['Image plate','Point detector']
    detType = wx.ComboBox(G2frame.dataDisplay,value=data['DetType'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    detType.Bind(wx.EVT_COMBOBOX, OnDetType)
    sqBox.Add(detType,0)
    if data['DetType'] == 'Image plate':
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' IP transmission coeff.: '),0,wx.ALIGN_CENTER_VERTICAL)
        obliqCoeff = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['ObliqCoeff']))
        obliqCoeff.Bind(wx.EVT_TEXT_ENTER,OnObliqCoeff)        
        obliqCoeff.Bind(wx.EVT_KILL_FOCUS,OnObliqCoeff)
        sqBox.Add(obliqCoeff,0)
    mainSizer.Add(sqBox,0)
        
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Ruland width: '),0,wx.ALIGN_CENTER_VERTICAL)    
    rulandSldr = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
        value=int(1000*data['Ruland']))
    sqBox.Add(rulandSldr,1,wx.EXPAND)
    rulandSldr.Bind(wx.EVT_SLIDER, OnRulSlider)
    rulandWdt = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['Ruland']))
    rulandWdt.Bind(wx.EVT_TEXT_ENTER,OnRulandWdt)        
    rulandWdt.Bind(wx.EVT_KILL_FOCUS,OnRulandWdt)
    sqBox.Add(rulandWdt,0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(sqBox,0,wx.ALIGN_LEFT|wx.EXPAND)
    
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    lorch = wx.CheckBox(parent=G2frame.dataDisplay,label='Lorch damping?')
    lorch.SetValue(data['Lorch'])
    lorch.Bind(wx.EVT_CHECKBOX, OnLorch)
    sqBox.Add(lorch,0,wx.ALIGN_CENTER_VERTICAL)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Scaling q-range: '),0,wx.ALIGN_CENTER_VERTICAL)
    SQmin = wx.TextCtrl(G2frame.dataDisplay,value='%.1f'%(data['QScaleLim'][0]))
    SQmin.Bind(wx.EVT_TEXT_ENTER,OnSQmin)        
    SQmin.Bind(wx.EVT_KILL_FOCUS,OnSQmin)    
    sqBox.Add(SQmin,0)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' to '),0,wx.ALIGN_CENTER_VERTICAL)
    SQmax = wx.TextCtrl(G2frame.dataDisplay,value='%.1f'%(data['QScaleLim'][1]))
    SQmax.Bind(wx.EVT_TEXT_ENTER,OnSQmax)        
    SQmax.Bind(wx.EVT_KILL_FOCUS,OnSQmax)
    sqBox.Add(SQmax,0)
    resetQ = wx.CheckBox(parent=G2frame.dataDisplay,label='Reset?')
    sqBox.Add(resetQ,0)
    resetQ.Bind(wx.EVT_CHECKBOX, OnResetQ)
    
    mainSizer.Add(sqBox,0)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)
    