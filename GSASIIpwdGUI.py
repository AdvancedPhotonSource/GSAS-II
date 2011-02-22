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
import math
import time
import cPickle
import GSASIIpath
import GSASIIpeak as G2pk
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIindex as G2indx
import GSASIIplot as G2plt
import GSASIIgrid as G2gd

VERY_LIGHT_GREY = wx.Colour(235,235,235)

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
       
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
                
    def OnPeakFit(event):
        SaveState()
        print 'Peak Fitting - Do one cycle of peak fitting'
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
        OK,smin,Rwp,runtime,GoOn = G2pk.DoPeakFit(peaks,background,limits,inst,data)
        UpdatePeakGrid(self,peaks)
        G2plt.PlotPatterns(self)
        if not OK:
            print 'Refinement failed'
            dlg = wx.MessageDialog(self, 'Do you want to reload now?', 'Refinement failed',  wx.YES_NO)
            try:
                if dlg.ShowModal() == wx.ID_YES:
                    DoUnDo()
                    self.dataFrame.UnDo.Enable(False)
            finally:
                dlg.Destroy()
        else:
            self.dataFrame.UnDo.Enable(True)
            print "%s%7.2f%s%12.6g" % ('Rwp = ',Rwp,'%, Smin = ',smin)
            print "%s%8.3f%s " % ('fitpeak time =',runtime,'s')
            print 'finished'
        return
        
    def OnAutoPeakFit(event):
        SaveState()
        print 'AutoPeak Fitting - run until minimized'
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
        smin = 1.0e15
        GoOn = True
        while GoOn:
            osmin = smin
            OK,smin,Rwp,runtime,GoOn = G2pk.DoPeakFit(peaks,background,limits,inst,data)
            UpdatePeakGrid(self,peaks)
            if not OK:
                break
            G2plt.PlotPatterns(self)
            print "%s%7.2f%s%12.6g" % ('Rwp = ',Rwp,'%, Smin = ',smin)
            rat = (osmin-smin)/smin
            if rat < 1.0e-4: GoOn = False
        if not OK:
            print 'Refinement failed'
            dlg = wx.MessageDialog(self, 'Do you want to reload now?', 'Refinement failed',  wx.YES_NO)
            try:
                if dlg.ShowModal() == wx.ID_YES:
                    DoUnDo()
                    self.dataFrame.UnDo.Enable(False)
            finally:
                dlg.Destroy()
        else:
            self.dataFrame.UnDo.Enable(True)
            print "%s%8.3f%s " % ('fitpeak time =',runtime,'s per cycle')
            print 'finished'
        return        

    def RefreshPeakGrid(event):
        event.StopPropagation()
        data = self.PeakTable.GetData()
        T = []
        for peak in data:T.append(peak[0])
        D = dict(zip(T,data))
        T.sort()
        X = []
        for key in T: X.append(D[key])
        data = X        
        G2plt.PlotPatterns(self)
        
    def setBackgroundColors():
       for r in range(self.dataDisplay.GetNumberRows()):
           for c in range(self.dataDisplay.GetNumberCols()):
               if self.dataDisplay.GetColLabelValue(c) in ['position','intensity','sigma','gamma']:
                   if float(self.dataDisplay.GetCellValue(r,c)) < 0.:
                       self.dataDisplay.SetCellBackgroundColour(r,c,wx.RED)
                   else:
                       self.dataDisplay.SetCellBackgroundColour(r,c,wx.WHITE)
                       
    def RefineSelect(event):
        data = self.PatternTree.GetItemPyData(self.PickId)
        r,c =  event.GetRow(),event.GetCol()
        if r < 0 and self.dataDisplay.GetColLabelValue(c) == 'refine':
            self.dataDisplay.SelectCol(c,False)
        
                       
    def RowSelect(event):
        r,c =  event.GetRow(),event.GetCol()
        if r < 0 and c < 0:
            if self.dataDisplay.IsSelection():
                self.dataDisplay.ClearSelection()
        elif c < 0:                   #only row clicks
            if event.ControlDown():                    
                if r in self.dataDisplay.GetSelectedRows():
                    self.dataDisplay.DeselectRow(r)
                else:
                    self.dataDisplay.SelectRow(r,True)
            elif event.ShiftDown():
                for row in range(r+1):
                    self.dataDisplay.SelectRow(row,True)
            else:
                self.dataDisplay.ClearSelection()
                self.dataDisplay.SelectRow(r,True)                
        
                           
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
                    self.dataFrame.AutoPeakFit.Enable(False)
                        
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
    self.Bind(wx.EVT_MENU, OnPeakFit, id=G2gd.wxID_PEAKFIT)
    self.Bind(wx.EVT_MENU, OnAutoPeakFit, id=G2gd.wxID_AUTOPEAKFIT)
    self.dataFrame.PeakFit.Enable(False)
    self.dataFrame.AutoPeakFit.Enable(False)
    if data:
        self.dataFrame.PeakFit.Enable(True)
        self.dataFrame.AutoPeakFit.Enable(True)
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
    self.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)                 
    self.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefineSelect)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([535,350])
        
def UpdateBackgroundGrid(self,data):
    if self.dataDisplay:
        self.dataFrame.Clear()
    BackId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Background')
    maxTerm = 9
    Types = [wg.GRID_VALUE_CHOICE+':chebyschev,another,more',
        wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_NUMBER+':1,'+str(maxTerm)]
    for i in range(maxTerm):
        Types.append(wg.GRID_VALUE_FLOAT+':10,3')
    
    def RefreshBackgroundGrid(event):
        data = self.BackTable.GetData()
        M = len(data[0])
        N = data[0][2]+3
        item = data[0]
        if N > M:       #add terms
            for i in range(M,N): 
                item.append(0.0)
                self.BackTable.SetColLabelValue(i,str(i-2))
            data = [item]
            msg = wg.GridTableMessage(self.BackTable, 
                wg.GRIDTABLE_NOTIFY_COLS_APPENDED,0,N-M)
            self.dataDisplay.ProcessTableMessage(msg)                         
        elif N < M:     #delete terms
            new = []
            for i in range(N):
                new.append(item[i])
            data = [new]
            msg = wg.GridTableMessage(self.BackTable, 
                wg.GRIDTABLE_NOTIFY_COLS_DELETED,0,M-N)
            self.dataDisplay.ProcessTableMessage(msg)
        self.PatternTree.SetItemPyData(BackId,data)
        event.StopPropagation()
                  
    self.BackTable = []
    N = len(data[0])
    M = data[0][2]
    colLabels = ['function','refine','Nterms']
    rowLabels=['background']
    for i in range(M): colLabels.append(str(i+1))
    self.BackTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Background')
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    gridPanel = wx.Panel(self.dataFrame)
    self.dataDisplay = G2gd.GSGrid(gridPanel)                
    self.dataDisplay.SetTable(self.BackTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshBackgroundGrid)                
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(self.dataDisplay,0)
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
    Types = [wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3']
    self.LimitsTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Limits')
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.LimitsTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshLimitsGrid)                
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([230,120])
    
def UpdateInstrumentGrid(self, data):
    if self.dataDisplay:
        self.dataFrame.Clear()
    Ka2 = False
    if len(data[0]) == 13: 
        Ka2 = True
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    
    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List'))
            ins = data[1]
            if 'P' in ins[0]:                                       #update powder peak parameters
                for peak in peaks:
                    if Ka2:
                        peak[4] = ins[6]*tand(peak[0]/2.0)**2+ins[7]*tand(peak[0]/2.0)+ins[8]
                        peak[6] = ins[9]/cosd(peak[0]/2.0)+ins[10]*tand(peak[0]/2.0)
                    else:
                        peak[4] = ins[4]*tand(peak[0]/2.0)**2+ins[5]*tand(peak[0]/2.0)+ins[6]
                        peak[6] = ins[7]/cosd(peak[0]/2.0)+ins[8]*tand(peak[0]/2.0)
                        
    def OnReset(event):
        if Ka2:
            data[1][6:12] = data[0][6:12]
        else:
            data[1][4:10] = data[0][4:10]
        RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
        UpdateInstrumentGrid(self, data)
        
    self.InstrumentTable = []
    if 'P' in data[1][0]:                   #powder data
        self.dataFrame.SetMenuBar(self.dataFrame.InstMenu)
        if not self.dataFrame.GetStatusBar():
            Status = self.dataFrame.CreateStatusBar()
        self.Bind(wx.EVT_MENU, OnReset, id=G2gd.wxID_INSTPRMRESET)
        if Ka2:
            Types = [wg.GRID_VALUE_CHOICE+":PXC,PNC,PNT",wg.GRID_VALUE_FLOAT+':10,6',wg.GRID_VALUE_FLOAT+':10,6',               #type, lam-1 & lam-2
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #zero, ratio, pola
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #u,v,w
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,2']
        else:
            Types = [wg.GRID_VALUE_CHOICE+":PXC,PNC,PNT",wg.GRID_VALUE_FLOAT+':10,6',               #type & lam-1 
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #zero, pola
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #u,v,w
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,2']
        colLabels = data[3]
        rowLabels = ['default','changed','refine']
        self.InstrumentTable = G2gd.Table(data[:-1],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        self.dataFrame.SetLabel('Instrument Parameters')
        gridPanel = wx.Panel(self.dataFrame)
        self.dataDisplay = G2gd.GSGrid(gridPanel)                
        self.dataDisplay.SetTable(self.InstrumentTable, True)
        self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshInstrumentGrid)                
        self.dataDisplay.SetMargins(0,0)
        self.dataDisplay.AutoSizeColumns(False)
        beg = 4
        if Ka2: beg = 6
        for i in range(len(data[2])):
            if i < beg or i == beg+6:
                self.dataDisplay.SetCellRenderer(2,i,wg.GridCellStringRenderer())
                self.dataDisplay.SetCellValue(2,i,'')
                self.dataDisplay.SetReadOnly(2,i,isReadOnly=True)
            else:
                self.dataDisplay.SetCellRenderer(2,i,wg.GridCellBoolRenderer())
                self.dataDisplay.SetCellEditor(2,i,wg.GridCellBoolEditor())
    else:                       #single crystal data
        Types = [wg.GRID_VALUE_CHOICE+":SXC,SNC,SNT",wg.GRID_VALUE_FLOAT+':10,6']
        colLabels = data[2]
        rowLabels = ['original','changed']
        self.InstrumentTable = G2gd.Table(data[:-1],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        self.dataFrame.SetLabel('Instrument Parameters')
        gridPanel = wx.Panel(self.dataFrame)
        self.dataDisplay = G2gd.GSGrid(gridPanel)                
        self.dataDisplay.SetTable(self.InstrumentTable, True)
        self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshInstrumentGrid)                
        self.dataDisplay.SetMargins(0,0)
        self.dataDisplay.AutoSizeColumns(False)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(self.dataDisplay,0)
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
    
    if data['Type'] == 'Debye-Scherrer':
        parms = [['DisplaceX',' Sample X displacement: ','%.4f',],
            ['DisplaceY',' Sample Y displacement: ','%.4f',],
            ['Absorption',' Sample absorption: ','%.4f',],]
    elif data['Type'] == 'Bragg-Brentano':
        parms = [['Shift',' Sample displacement: ','%.4f',],
            ['Transparency',' Sample transparency: ','%.4f'],]
    parms.append(['Temperature',' Sample temperature: ','%.2f'])
    parms.append(['Pressure',' Sample pressure: ','%.3f'])
    parms.append(['Humidity',' Sample humidity: ','%.1f'])
    parms.append(['Voltage',' Sample voltage: ','%.3f'])
    parms.append(['Force',' Applied load: ','%.3f'])
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
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Sample parameters: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleRef = wx.CheckBox(self.dataDisplay,label=' Histogram scale factor: ')
    scaleRef.SetValue(data['Scale'][1])
    scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
    scaleSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
    scaleVal = wx.TextCtrl(self.dataDisplay,wx.ID_ANY,
        '%.4f'%(data['Scale'][0]),style=wx.TE_PROCESS_ENTER)
    scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
    scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
    scaleSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(scaleSizer)
    mainSizer.Add((0,5),0)
    typeSizer = wx.BoxSizer(wx.HORIZONTAL)
    choices = ['Debye-Scherrer','Bragg-Brentano',]
    histoType = wx.ComboBox(self.dataDisplay,wx.ID_ANY,value=data['Type'],choices=choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    histoType.Bind(wx.EVT_COMBOBOX, OnHistoType)
    typeSizer.Add(histoType)
    mainSizer.Add(typeSizer)
    mainSizer.Add((0,5),0)
    
    for parm in parms:
        parmSizer = wx.BoxSizer(wx.HORIZONTAL)
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
        parmSizer.Add(parmVal,0,wx.ALIGN_CENTER_VERTICAL)
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
        data = self.IndexPeaksTable.GetData()
        self.PatternTree.SetItemPyData(IndexId,data)
        
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
                if self.IndexPeaksTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(self.IndexPeaksTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(self.IndexPeaksTable.GetNumberRows()): data[row][col]=False
            
    if self.dataDisplay:
        self.dataFrame.Clear()
    self.dataFrame.SetMenuBar(self.dataFrame.IndPeaksMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
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
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshIndexPeaksGrid)
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
            stVol = int(startVol.GetValue())
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
            Zero = 0.1
        controls[1] = Zero
        zero.SetValue("%.2f"%(Zero))
        
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
        print controls[5]
        ibrav = bravaisSymb.index(controls[5])
        dmin = G2indx.getDmin(peaks)-0.005
        Lhkl,M20,X20 = G2indx.refinePeaks(peaks,ibrav,A)
        controls[6:12] = G2lat.A2cell(A)
        controls[12] = G2lat.calc_V(A)
        data = [controls,bravais,cells,dmin]
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'),data)
        self.HKL = G2lat.GenHBravais(dmin,ibrav,A)
        UpdateUnitCellsGrid(self,data)
        print "%s%10.3f" % ('refinement M20 = ',M20)
        print 'unindexed lines = ',X20
        cellPrint(ibrav,A)
        for hkl in self.HKL:
            hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
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
        
    def MakeNewPhase(event):
        if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            PhaseName = dlg.GetValue()
        dlg.Destroy()
        cells = self.PatternTree.GetItemPyData(UnitCellsId)[2]
        for Cell in cells:
            if Cell[-1]:
                break
        cell = Cell[2:10]        
        sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
        E,SGData = G2spc.SpcGroup(spaceGroups[cell[0]])
        self.PatternTree.SetItemPyData(sub, \
            {'General':{'Name':'phase name','Type':'nuclear','SGData':SGData,
            'Cell':[False,]+cell[1:],
            'Pawley dmin':0.25},'Atoms':[],'Drawing':{},'Histograms':{}})
        Status.SetStatusText('Change space group if needed')
            
    def RefreshUnitCellsGrid(event):
        cells,dmin = self.PatternTree.GetItemPyData(UnitCellsId)[2:]
        r,c =  event.GetRow(),event.GetCol()
        if cells:
            if c == 2:
                for i in range(len(cells)):
                    cells[i][-1] = False
                    UnitCellsTable.SetValue(i,c,False)
                UnitCellsTable.SetValue(r,c,True)
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
    startVol = wx.TextCtrl(self.dataDisplay,value=str(controls[3]),style=wx.TE_PROCESS_ENTER)
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
    littleSizer = wx.FlexGridSizer(1,3,5,5)
    littleSizer.Add(wx.StaticText(self.dataDisplay,label=" Zero offset"),0,wx.ALIGN_CENTER_VERTICAL)
    zero = wx.TextCtrl(self.dataDisplay,value=str(controls[1]),style=wx.TE_PROCESS_ENTER)
    zero.Bind(wx.EVT_TEXT_ENTER,OnZero)
    zero.Bind(wx.EVT_KILL_FOCUS,OnZero)
    littleSizer.Add(zero,0,wx.ALIGN_CENTER_VERTICAL)
    zeroVar = wx.CheckBox(self.dataDisplay,label="Vary? (not implemented)")
    zero.SetValue("%.2f"%(controls[1]))
    zeroVar.Bind(wx.EVT_CHECKBOX,OnZeroVar)
    littleSizer.Add(zeroVar,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Cell Refinement: '),0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(1,2,5,5)
    littleSizer.Add(wx.StaticText(self.dataDisplay,label=" Bravais lattice"),0,wx.ALIGN_CENTER_VERTICAL)
    bravSel = wx.Choice(self.dataDisplay,choices=bravaisSymb)
    bravSel.SetSelection(bravaisSymb.index(controls[5]))
    bravSel.Bind(wx.EVT_CHOICE,OnBravSel)
    littleSizer.Add(bravSel,0,wx.ALIGN_CENTER_VERTICAL)
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
        Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',
            wg.GRID_VALUE_FLOAT+':10,2']
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
        gridDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshUnitCellsGrid)
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