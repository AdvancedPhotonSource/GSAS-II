#GSASII - data display routines
import wx
import wx.grid as wg
import math
import time
import cPickle
import GSASIIpath
import GSASIIpeak as G2pk
import GSASIIlattice as G2lat
import GSASIIindex as G2indx
import GSASIIplot as G2plt
import GSASIIgrid as G2gd

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
       
def UpdatePeakGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    
    def OnUnDo(event):
        DoUnDo()
        self.dataFrame.UnDo.Enable(False)
        
    def DoUnDo():
        print 'Undo last refinement'
        file = open('GSASII.save','rb')
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
        
    def OnPeakFit(event):
        self.SaveState()
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
        self.SaveState()
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
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([550,350])
        
def UpdateBackgroundGrid(self,data):
    BackId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Background')
    
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
                  
    self.dataFrame.setSizePosLeft([700,150])
    maxTerm = 9
    self.BackTable = []
    N = len(data[0])
    M = data[0][2]
    colLabels = ['function','refine','Nterms']
    rowLabels=['background']
    for i in range(M): colLabels.append(str(i+1))
    Types = [wg.GRID_VALUE_CHOICE+':chebyschev,another,more',
        wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_NUMBER+':1,'+str(maxTerm)]
    for i in range(maxTerm):
        Types.append(wg.GRID_VALUE_FLOAT+':10,3')
    self.BackTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Background')
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshBackgroundGrid)                
    self.dataDisplay.SetTable(self.BackTable, True)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
        
def UpdateLimitsGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.setSizePosLeft([250,150])
    LimitId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Limits')
    def RefreshLimitsGrid(event):
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
    
def UpdateInstrumentGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    Ka2 = False
    Xwid = 700
    if len(data[0]) == 13: 
        Ka2 = True
        Xwid = 840        
    self.dataFrame.setSizePosLeft([Xwid,170])
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    InstId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters')
    
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
        self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)                
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
        self.InstrumentTable = Table(data[:-1],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        self.dataFrame.SetLabel('Instrument Parameters')
        self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)                
        self.dataDisplay.SetTable(self.InstrumentTable, True)
        self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshInstrumentGrid)                
        self.dataDisplay.SetMargins(0,0)
        self.dataDisplay.AutoSizeColumns(False)
                
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
        self.dataDisplay.Destroy()
    self.dataFrame.setSizePosLeft([500,300])
    self.dataFrame.SetMenuBar(self.dataFrame.IndPeaksMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    self.Bind(wx.EVT_MENU, OnReload, id=G2gd.wxID_INDXRELOAD)
    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
    self.IndexPeaksTable = []
    if data:
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

def UpdateUnitCellsGrid(self, data):
    UnitCellsId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Unit Cells List')
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
        
    def OnRefineCell(event):
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
             
        bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
            'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
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
        
    def OnIndexPeaks(event):
        PatternId = self.PatternId    
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
        if not peaks:
            self.ErrorDialog('No peaks!', 'Nothing to index!')
            return
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
        self.dataFrame.IndexPeaks.Enable(False)
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
        for i in range(4,13):
            self.UnitCellsTable.SetValue(i,1,controls[i])
        self.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin])
        self.dataDisplay.ForceRefresh()
        self.dataFrame.RefineCell.Enable(True)
            
    def RefreshUnitCellsGrid(event):
        cells,dmin = self.PatternTree.GetItemPyData(UnitCellsId)[2:]
        r,c =  event.GetRow(),event.GetCol()
        if cells:
            if c == 6:
                for i in range(min(self.UnitCellsTable.GetNumberRows(),len(cells))):
                    cells[i][-1] = False
                    self.UnitCellsTable.SetValue(i,c,0)
                self.UnitCellsTable.SetValue(r,c,1)
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
        controls = []
        bravais = [0,0,0,0,0,0,0, 0,0,0,0,0,0,0]
        table = self.UnitCellsTable.GetData()
        for i,row in enumerate(table):
            if i in [0,4]:
                if row[1]:
                    controls.append(1)
                else:
                    controls.append(0)
            elif i in [1,3]:
                controls.append(float(row[1]))
            elif i in [2]:
                controls.append(int(row[1]))
            elif i in [5]:
                controls.append(row[1])
            elif i in [6,7,8,9,10,11]:
                if controls[5] in bravaisSymb[:3]:              #cubic
                    if i in [6]:
                        controls.append(float(row[1]))
                        controls.append(float(row[1]))
                        controls.append(float(row[1]))
                        controls.append(90.)
                        controls.append(90.)
                        controls.append(90.)
                elif controls[5] in bravaisSymb[3:7]:           #hexagonal & tetragonal
                    if i in [6]:
                        controls.append(float(row[1]))
                        controls.append(float(row[1]))
                    elif i in [8]:
                        controls.append(float(row[1]))
                        controls.append(90.)
                        controls.append(90.)
                        if controls[5] in bravaisSymb[3:5]:     #hexagonal
                            controls.append(120.)
                        else:                                   #tetragonal
                            controls.append(90.)
                elif controls[5] in bravaisSymb[7:13]:          #orthorhombic & monoclinic
                    if i in [6,7,8]:
                        controls.append(float(row[1]))
                    if i in [9,10,11]:
                        if controls[5] in bravaisSymb[7:11]:
                            controls.append(90.)
                            controls.append(90.)
                            controls.append(90.)
                            break
                        else:
                            if i in [9,11]:
                                controls.append(90.)
                            else:
                                controls.append(float(row[1]))
                else:                                           #triclinic
                    controls.append(float(row[1]))
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))        #volume        
        for i,row in enumerate(table):
            if i < 14:
                bravais[i] = int(row[2])
            else:
                break
        if controls[4]:
            for i in range(6,13):
                self.UnitCellsTable.SetValue(i,1,controls[i])
        self.dataDisplay.ForceRefresh()
        if controls[4] and not False in controls[6:12]:
            self.dataFrame.RefineCell.Enable(True)
        else:
            self.dataFrame.RefineCell.Enable(False)
        data = [controls,bravais,cells,dmin]                    
        self.PatternTree.SetItemPyData(UnitCellsId,data)
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.IndexMenu)
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    self.Bind(wx.EVT_MENU, OnIndexPeaks, id=G2gd.wxID_INDEXPEAKS)
    self.Bind(wx.EVT_MENU, CopyUnitCell, id=G2gd.wxID_COPYCELL)
    self.Bind(wx.EVT_MENU, OnRefineCell, id=G2gd.wxID_REFINECELL)
    self.UnitCellsTable = []
    controls,bravais,cells,dmin = data
    if cells:
        self.dataFrame.setSizePosLeft([900,320])
    else:
        self.dataFrame.setSizePosLeft([280,320])
    if len(controls) < 13:
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    self.PatternTree.SetItemPyData(UnitCellsId,data)
    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
    if cells:
        colLabels = ['controls','value','try','Bravais cells',
            'M20','X20','use','Bravais','a','b','c','alpha','beta','gamma','Volume']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_FLOAT+":10,1",
            wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,wg.GRID_VALUE_FLOAT+':10,2',
            wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',
            wg.GRID_VALUE_FLOAT+':10,2']
    else:
        colLabels = ['controls','value','try','Bravais cells']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,
            wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING]
    controlNames = ['Vary zero?','Max zero error','Max Nc/Nobs','Start volume','refine cell?',
        'bravais','a=','b=','c=','alpha=','beta=','gamma=','Volume=']
    bravaisNames = ['Cubic-F','Cubic-I','Cubic-P','Trigonal-R','Trigonal/Hexagonal-P',
        'Tetragonal-I','Tetragonal-P','Orthorhombic-F','Orthorhombic-I','Orthorhombic-C',
        'Orthorhombic-P','Monoclinic-C','Monoclinic-P','Triclinic']
    rowLabels = []
    table = []
    numRows = max(len(bravais),len(cells))
    for i in range(numRows):
        rowLabels.append('')
        if i < 13:
            row = [controlNames[i],controls[i],bravais[i],bravaisSymb[i]]
        elif i < 14:
            row = ['','',bravais[i],bravaisSymb[i]]
        else:
            row = ['','','','']
        if cells:
            if i < len(cells):
                cell = cells[i]
                row += cell[0:2]+[cell[-1]]+[bravaisSymb[cell[2]]]+cell[3:10]
                if cell[-1]:
                    A = G2lat.cell2A(cell[3:9])
                    self.HKL = G2lat.GenHBravais(dmin,cell[2],A)
                    for hkl in self.HKL:
                        hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))
            else:
                row += 14*['',]
        table.append(row)
    self.UnitCellsTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Unit Cells List')
    self.dataDisplay = G2gd.GSGrid(parent=self.dataFrame)                
    self.dataDisplay.SetTable(self.UnitCellsTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshUnitCellsGrid)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.SetRowLabelSize(0)
    self.dataDisplay.SetCellRenderer(0,1,wg.GridCellBoolRenderer())
    self.dataDisplay.SetCellEditor(0,1,wg.GridCellBoolEditor())
    self.dataDisplay.SetCellRenderer(1,1,wg.GridCellFloatRenderer(5,2))
    self.dataDisplay.SetCellEditor(1,1,wg.GridCellFloatEditor(5,0))
    self.dataDisplay.SetCellRenderer(2,1,wg.GridCellNumberRenderer())
    self.dataDisplay.SetCellEditor(2,1,wg.GridCellNumberEditor(1,10))
    self.dataDisplay.SetCellRenderer(3,1,wg.GridCellFloatRenderer(5,0))
    self.dataDisplay.SetCellEditor(3,1,wg.GridCellFloatEditor(5,2))
    self.dataDisplay.SetCellRenderer(4,1,wg.GridCellBoolRenderer())
    self.dataDisplay.SetCellEditor(4,1,wg.GridCellBoolEditor())
    self.dataDisplay.SetCellRenderer(5,1,wg.GridCellStringRenderer())
    self.dataDisplay.SetCellEditor(5,1,wg.GridCellChoiceEditor(bravaisSymb,False))
    for i in range(6,9):
        self.dataDisplay.SetCellRenderer(i,1,wg.GridCellFloatRenderer(10,5))
        self.dataDisplay.SetCellEditor(i,1,wg.GridCellFloatEditor(10,5))
    for i in range(9,13):
        self.dataDisplay.SetCellRenderer(i,1,wg.GridCellFloatRenderer(10,3))
        self.dataDisplay.SetCellEditor(i,1,wg.GridCellFloatEditor(10,3))
    for i in range(14):
        self.dataDisplay.SetReadOnly(i,0,isReadOnly=True)
        self.dataDisplay.SetReadOnly(i,3,isReadOnly=True)
    if cells:
        self.dataFrame.CopyCell.Enable(True)
        for r in range(max(len(cells),14)):
            if r > 12:
                self.dataDisplay.SetCellRenderer(r,0,wg.GridCellStringRenderer())                    
                self.dataDisplay.SetCellRenderer(r,1,wg.GridCellStringRenderer())
            if r > 13:
                self.dataDisplay.SetCellRenderer(r,2,wg.GridCellStringRenderer())
            for c in range(4,15):
                if r >= len(cells):
                    self.dataDisplay.SetCellRenderer(r,c,wg.GridCellStringRenderer())
                if c != 6:
                    self.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    self.dataDisplay.AutoSizeColumns(False)
    if controls[4] and not False in controls[6:12]:
        self.dataFrame.RefineCell.Enable(True)
    else:
        self.dataFrame.RefineCell.Enable(False)
        
