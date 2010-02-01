#GSASII - data display routines
import wx
import wx.grid as wg
import matplotlib as mpl
import math
import time
import GSASIIcomp as G2cmp
import GSASIIspc as G2spc
import GSASIIElem as G2elem

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
       
[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, 
] = [wx.NewId() for _init_coll_Atom_Items in range(2)]

[ wxID_IMCALIBRATE, wxID_IMINTEGRATE, 
] = [wx.NewId() for _init_coll_IMAGE_Items in range(2)]

class DataFrame(wx.Frame):
    def _init_coll_BlankMenu(self,parent):
        parent.Append(menu=self.Blank,title='')
        
    def _init_coll_AtomsMenu(self,parent):
        parent.Append(menu=self.AtomEdit, title='Add atom')
                   
    def _init_coll_ImageMenu(self,parent):
        parent.Append(menu=self.ImageEdit, title='Image Operations')
        
    def _init_coll_Atom_Items(self,parent):
        self.Add = parent.Append(help='',id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append empty atom')
        self.Add = parent.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert empty atom',
            help='Double left click on atom row to Insert before')
            
    def _init_coll_Image_Items(self,parent):
        self.Add = parent.Append(help='',id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        self.Add = parent.Append(id=wxID_IMINTEGRATE, kind=wx.ITEM_NORMAL,text='Integrate',
            help='')
            
        
    def _init_utils(self):
        self.BlankMenu = wx.MenuBar()
        
        self.AtomsMenu = wx.MenuBar()
        self.ImageMenu = wx.MenuBar()
        self.AtomEdit = wx.Menu(title='')
        self.ImageEdit = wx.Menu(title='')
        self._init_coll_AtomsMenu(self.AtomsMenu)
        self._init_coll_Atom_Items(self.AtomEdit)
        self._init_coll_ImageMenu(self.ImageMenu)
        self._init_coll_Image_Items(self.ImageEdit)
        
    def _init_ctrls(self, parent,name=None,size=None,pos=None):
        wx.Frame.__init__(self,parent=parent,style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX,
            size=size,pos=pos,title='GSAS-II data display')
        self._init_utils()
        if name:
            self.SetLabel(name)
        self.Show()
        
    def __init__(self,parent,data=None,name=None, size=None,pos=None):
        self._init_ctrls(parent,name,size,pos)
        self.data = data
        self.screenSize = wx.DisplaySize()
        Size = self.GetSize()
        xPos = self.screenSize[0]-Size[0]
        self.SetPosition(wx.Point(xPos,250))
        self.dirname = ''
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        screenSize = wx.DisplaySize()
        self.SetSize(Width)
        self.SetPosition(wx.Point(screenSize[0]-Width[0],250))
        
    def Clear(self):
        self.ClearBackground()
        self.DestroyChildren()
                   
class GSGrid(wg.Grid):
    def __init__(self, parent, name=''):
        wg.Grid.__init__(self,parent,-1,name=name)                    
        self.SetSize(parent.GetClientSize())
            
    def Clear(self):
        wg.Grid.ClearGrid(self)
        
    def SetCellStyle(self,r,c,color="white",readonly=True):
        self.SetCellBackgroundColour(r,c,color)
        self.SetReadOnly(r,c,isReadOnly=readonly)
        
class GSNoteBook(wx.Notebook):
    def __init__(self, parent, name='',size = None):
        wx.Notebook.__init__(self, parent, -1, name=name, style= wx.BK_TOP)
        if size: self.SetSize(size)
                                                      
    def Clear(self):        
        GSNoteBook.DeleteAllPages(self)
        
class Table(wg.PyGridTableBase):
    def __init__(self, data=[], rowLabels=None, colLabels=None, types = None):
        wg.PyGridTableBase.__init__(self)
        self.colLabels = colLabels
        self.rowLabels = rowLabels
        self.dataTypes = types
        self.data = data
        
    def AppendRows(self, numRows=1):
        self.data.append([])
        return True
        
    def CanGetValueAs(self, row, col, typeName):
        if self.dataTypes:
            colType = self.dataTypes[col].split(':')[0]
            if typeName == colType:
                return True
            else:
                return False
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def DeleteRow(self,pos):
        data = self.GetData()
        self.SetData([])
        new = []
        for irow,row in enumerate(data):
            if irow <> pos:
                new.append(row)
        self.SetData(new)
        
    def GetColLabelValue(self, col):
        if self.colLabels:
            return self.colLabels[col]
            
    def GetData(self):
        data = []
        for row in range(self.GetNumberRows()):
            data.append(self.GetRowValues(row))
        return data
        
    def GetNumberCols(self):
        try:
            return len(self.colLabels)
        except TypeError:
            return None
        
    def GetNumberRows(self):
        return len(self.data)
        
    def GetRowLabelValue(self, row):
        if self.rowLabels:
            return self.rowLabels[row]
        
    def GetRowValues(self, row):
        data = []
        for col in range(self.GetNumberCols()):
            data.append(self.GetValue(row, col))
        return data
        
    def GetTypeName(self, row, col):
        try:
            return self.dataTypes[col]
        except TypeError:
            return None

    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return None
            
    def InsertRows(self, pos, rows):
        for row in range(rows):
            self.data.insert(pos,[])
            pos += 1
        
    def IsEmptyCell(self,row,col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
        
    def OnKeyPress(self, event):
        dellist = self.GetSelectedRows()
        if event.GetKeyCode() == wx.WXK_DELETE and dellist:
            grid = self.GetView()
            for i in dellist: grid.DeleteRow(i)
                
    def SetColLabelValue(self, col, label):
        numcols = self.GetNumberCols()
        if col > numcols-1:
            self.colLabels.append(label)
        else:
            self.colLabels[col]=label
        
    def SetData(self,data):
        for row in range(len(data)):
            self.SetRowValues(row,data[row])
                
    def SetRowLabelValue(self, row, label):
        self.rowLabels[row]=label
            
    def SetRowValues(self,row,data):
        self.data[row] = data
            
    def SetValue(self, row, col, value):
        def innerSetValue(row, col, value):
            try:
                self.data[row][col] = value
            except TypeError:
                return
            except IndexError:
                print row,col,value
                # add a new row
                if row > self.GetNumberRows():
                    self.data.append([''] * self.GetNumberCols())
                elif col > self.GetNumberCols():
                    for row in range(self.GetNumberRows):
                        self.data[row].append('')
                print self.data
                self.data[row][col] = value
        innerSetValue(row, col, value)
        
        
def UpdateNotebook(self,data):        
    if data:
        self.dataFrame.SetLabel('Notebook')
        self.dataDisplay = wx.TextCtrl(parent=self.dataFrame,size=self.dataFrame.GetClientSize(),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
        for line in data:
            self.dataDisplay.AppendText(line+"\n")
            self.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
            
def UpdateControls(self,data):
    if data:
        self.dataFrame.SetLabel('Controls')
        
     
def UpdateComments(self,data):                   
    if data:
        self.dataFrame.SetLabel('Comments')
        self.dataDisplay = wx.TextCtrl(parent=self.dataFrame,size=self.dataFrame.GetClientSize(),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
        for line in data:
            self.dataDisplay.AppendText(line+"\n")
     
def UpdatePeakGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    
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
        self.PlotPatterns()
        
    def setBackgroundColors():
       for r in range(self.dataDisplay.GetNumberRows()):
           for c in range(self.dataDisplay.GetNumberCols()):
               if self.dataDisplay.GetColLabelValue(c) in ['position','intensity','sigma','gamma','SH/L']:
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
        self.PlotPatterns()
            
    self.PickTable = []
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','refine','intensity','refine','sigma','refine','gamma','refine','SH/L','refine']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
    T = []
    for peak in data:T.append(peak[0])
    D = dict(zip(T,data))
    T.sort()
    X = []
    for key in T: X.append(D[key])
    data = X        
    self.PatternTree.SetItemPyData(self.PickId,data)
    self.PeakTable = Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Peak List')
    self.dataDisplay = GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.PeakTable, True)
    setBackgroundColors()                         
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshPeakGrid)
    self.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)                 
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    self.dataFrame.setSizePosLeft([650,350])
        
def UpdateBackgroundGrid(self,data):
    BackId = GetPatternTreeItemId(self,self.PatternId, 'Background')
    
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
    maxTerm = 7
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
    self.BackTable = Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Background')
    self.dataDisplay = GSGrid(parent=self.dataFrame)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshBackgroundGrid)                
    self.dataDisplay.SetTable(self.BackTable, True)
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
        
def UpdateLimitsGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.setSizePosLeft([250,150])
    LimitId = GetPatternTreeItemId(self,self.PatternId, 'Limits')
    def RefreshLimitsGrid(event):
        data = self.LimitsTable.GetData()
        old = data[0]
        new = data[1]
        new[0] = max(old[0],new[0])
        new[1] = max(new[0],min(old[1],new[1]))
        data = [old,new]
        self.PlotPatterns()
        
    self.LimitsTable = []
    colLabels = ['Tmin','Tmax']
    rowLabels = ['original','changed']
    Types = [wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3']
    self.LimitsTable = Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Limits')
    self.dataDisplay = GSGrid(parent=self.dataFrame)                
    self.dataDisplay.SetTable(self.LimitsTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshLimitsGrid)                
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)
    
def UpdateInstrumentGrid(self, data):
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    Ka2 = False
    Xwid = 720
    if len(data[0]) == 12: 
        Ka2 = True
        Xwid = 800        
    self.dataFrame.setSizePosLeft([Xwid,150])
    InstId = GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters')
    def RefreshInstrumentGrid(event):
        if event.GetRow() == 1:
            peaks = self.PatternTree.GetItemPyData(GetPatternTreeItemId(self,self.PatternId, 'Peak List'))
            ins = data[1]
            if 'P' in ins[0]:                                       #update powder peak parameters
                for peak in peaks:
                    if Ka2:
                        peak[4] = ins[6]*tand(peak[0]/2.0)**2+ins[7]*tand(peak[0]/2.0)+ins[8]
                        peak[6] = ins[9]/cosd(peak[0]/2.0)+ins[10]*tand(peak[0]/2.0)
                        peak[8] = ins[11]
                    else:
                        peak[4] = ins[4]*tand(peak[0]/2.0)**2+ins[5]*tand(peak[0]/2.0)+ins[6]
                        peak[6] = ins[7]/cosd(peak[0]/2.0)+ins[8]*tand(peak[0]/2.0)
                        peak[8] = ins[9]
        
    self.InstrumentTable = []
    if 'P' in data[1][0]:
        if Ka2:
            Types = [wg.GRID_VALUE_CHOICE+":PXC,PNC,PNT",wg.GRID_VALUE_FLOAT+':10,6',wg.GRID_VALUE_FLOAT+':10,6',               #type, lam-1 & lam-2
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #zero, ratio, pola
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #u,v,w
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,5']
        else:
            Types = [wg.GRID_VALUE_CHOICE+":PXC,PNC,PNT",wg.GRID_VALUE_FLOAT+':10,6',               #type & lam-1 
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #zero, pola
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3', #u,v,w
                wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_FLOAT+':10,5']
        colLabels = data[3]
        rowLabels = ['original','changed','refine']
        self.InstrumentTable = Table(data[:-1],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        self.dataFrame.SetLabel('Instrument Parameters')
        self.dataDisplay = GSGrid(parent=self.dataFrame)                
        self.dataDisplay.SetTable(self.InstrumentTable, True)
        self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshInstrumentGrid)                
        self.dataDisplay.SetMargins(0,0)
        self.dataDisplay.AutoSizeColumns(False)
        beg = 4
        if Ka2: beg = 6
        for i in range(len(data[2])):
            if i < beg:
                self.dataDisplay.SetCellRenderer(2,i,wg.GridCellStringRenderer())
                self.dataDisplay.SetCellValue(2,i,'')
                self.dataDisplay.SetReadOnly(2,i,isReadOnly=True)
            else:
                self.dataDisplay.SetCellRenderer(2,i,wg.GridCellBoolRenderer())
                self.dataDisplay.SetCellEditor(2,i,wg.GridCellBoolEditor())
    else:
        Types = [wg.GRID_VALUE_CHOICE+":SXC,SNC,SNT",wg.GRID_VALUE_FLOAT+':10,6']
        colLabels = data[2]
        rowLabels = ['original','changed']
        self.InstrumentTable = Table(data[:-1],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        self.dataFrame.SetLabel('Instrument Parameters')
        self.dataDisplay = GSGrid(parent=self.dataFrame)                
        self.dataDisplay.SetTable(self.InstrumentTable, True)
        self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshInstrumentGrid)                
        self.dataDisplay.SetMargins(0,0)
        self.dataDisplay.AutoSizeColumns(False)
                
def UpdateIndexPeaksGrid(self, data):
    IndexId = GetPatternTreeItemId(self,self.PatternId, 'Index Peak List')
    
    def RefreshIndexPeaksGrid(event):
        data = self.IndexPeaksTable.GetData()
        self.PatternTree.SetItemPyData(IndexId,data)
        
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
        elif event.GetKeyCode() == wx.WXK_DELETE:
            dlg = wx.MessageDialog(self, 'Delete Index Peak List?', ' ', wx.YES | wx.NO)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    oldlen = len(data)
                    data = []
                    self.PatternTree.SetItemPyData(IndexId,data)
                    self.dataDisplay.Clear() 
                    self.dataDisplay.Destroy()
                    self.IndexPeaksTable = []
            finally:
                dlg.Destroy()
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
    inst = self.PatternTree.GetItemPyData(GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
    self.IndexPeaksTable = []
    if not data:
        peaks = self.PatternTree.GetItemPyData(GetPatternTreeItemId(self,self.PatternId, 'Peak List'))
        for peak in peaks:
            dsp = inst[1]/(2.0*sind(peak[0]/2.0))
            data.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])
    else:
        cells = self.PatternTree.GetItemPyData(GetPatternTreeItemId(self,self.PatternId, 'Unit Cells List'))
        if cells:
            cellist = cells[2]
            dmin = cells[3]
            self.HKL = []
            for i,cell in enumerate(cellist):
                if cell[-1]:
                    ibrav = cell[2]
                    A = G2cmp.cell2A(cell[3:9])
                    self.HKL = G2cmp.GenHBravais(dmin,ibrav,A)
                    G2cmp.IndexPeaks(data,self.HKL)
                    for hkl in self.HKL:
                        hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','intensity','use','indexed','h','k','l','d-obs','d-calc']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_BOOL,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5']
    self.PatternTree.SetItemPyData(IndexId,data)
    self.IndexPeaksTable = Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Index Peak List')
    self.dataDisplay = GSGrid(parent=self.dataFrame)                
    self.dataDisplay.SetTable(self.IndexPeaksTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshIndexPeaksGrid)
    self.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)                 
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(False)

def UpdateUnitCellsGrid(self, data):
    UnitCellsId = GetPatternTreeItemId(self,self.PatternId, 'Unit Cells List')
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    def CopyUnitCell(event):
        if event.GetCol() == 6:
            row = event.GetRow()
            controls,bravais,cells,dmin = self.PatternTree.GetItemPyData(UnitCellsId)
            cell = cells[row][2:9]
            controls[4] = 1
            controls[5] = bravaisSymb[cell[0]]
            controls[6:12] = cell[1:8]
            controls[12] = G2cmp.calc_V(G2cmp.cell2A(controls[6:12]))
            for i in range(4,13):
                self.UnitCellsTable.SetValue(i,1,controls[i])
            self.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin])
            self.dataDisplay.ForceRefresh()
            self.RefineCell.Enable(True)
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
                A = G2cmp.cell2A(cells[r][3:9])
                self.HKL = G2cmp.GenHBravais(dmin,ibrav,A)
                for hkl in self.HKL:
                    hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))
                self.PlotPatterns()             
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
        controls.append(G2cmp.calc_V(G2cmp.cell2A(controls[6:12])))        #volume        
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
            self.RefineCell.Enable(True)
        else:
            self.RefineCell.Enable(False)
        data = [controls,bravais,cells,dmin]                    
        self.PatternTree.SetItemPyData(UnitCellsId,data)
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.UnitCellsTable = []
    controls,bravais,cells,dmin = data
    if cells:
        self.dataFrame.setSizePosLeft([900,320])
    else:
        self.dataFrame.setSizePosLeft([280,320])
    if len(controls) < 13:
        controls.append(G2cmp.calc_V(G2cmp.cell2A(controls[6:12])))
    self.PatternTree.SetItemPyData(UnitCellsId,data)
    inst = self.PatternTree.GetItemPyData(GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1]
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
                    A = G2cmp.cell2A(cell[3:9])
                    self.HKL = G2cmp.GenHBravais(dmin,cell[2],A)
                    for hkl in self.HKL:
                        hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))
            else:
                row += 14*['',]
        table.append(row)
    self.UnitCellsTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    self.dataFrame.SetLabel('Unit Cells List')
    self.dataDisplay = GSGrid(parent=self.dataFrame)                
    self.dataDisplay.SetTable(self.UnitCellsTable, True)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshUnitCellsGrid)
    self.dataDisplay.Bind(wg.EVT_GRID_CELL_RIGHT_DCLICK, CopyUnitCell)
    if cells:
        self.dataDisplay.GetGridWindow().SetToolTipString('Right d-click "use" box to copy cell to controls')
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
        self.RefineCell.Enable(True)
    else:
        self.RefineCell.Enable(False)
        
def UpdateHKLControls(self,data):
    
    def OnScaleSlider(event):
        scale = int(scaleSel.GetValue())/1000.
        scaleSel.SetValue(int(scale*1000.))
        data['Scale'] = scale*10.
        self.NewPlot = True
        self.PlotSngl()
        
    def OnLayerSlider(event):
        layer = layerSel.GetValue()
        data['Layer'] = layer
        self.NewPlot = True
        self.PlotSngl()
        
    def OnSelZone(event):
        data['Zone'] = zoneSel.GetValue()
        self.NewPlot = True
        self.PlotSngl()
        
    def OnSelType(event):
        data['Type'] = typeSel.GetValue()
        self.NewPlot = True
        self.PlotSngl()
        
    def SetStatusLine():
        Status.SetStatusText("look at me!!!")
                                      
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    Status = self.dataFrame.CreateStatusBar()
    SetStatusLine()
    zones = ['100','010','001']
    HKLmax = data['HKLmax']
    HKLmin = data['HKLmin']
    if data['ifFc']:
        typeChoices = ['Fosq','Fo','|DFsq|/sig','|DFsq|>sig','|DFsq|>3sig']
    else:
        typeChoices = ['Fosq','Fo']
    self.dataDisplay = wx.Panel(self.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Scale'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    scaleSel = wx.Slider(parent=self.dataDisplay,maxValue=1000,minValue=100,
        style=wx.SL_HORIZONTAL,value=int(data['Scale']*100))
    scaleSizer.Add(scaleSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    scaleSel.Bind(wx.EVT_SLIDER, OnScaleSlider)
    mainSizer.Add(scaleSizer,1,wx.EXPAND|wx.RIGHT)
    
    zoneSizer = wx.BoxSizer(wx.HORIZONTAL)
    zoneSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Zone  '),0,
        wx.ALIGN_CENTER_VERTICAL)
    zoneSel = wx.ComboBox(parent=self.dataDisplay,value=data['Zone'],choices=['100','010','001'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    zoneSel.Bind(wx.EVT_COMBOBOX, OnSelZone)
    zoneSizer.Add(zoneSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Plot type  '),0,
        wx.ALIGN_CENTER_VERTICAL)        
    typeSel = wx.ComboBox(parent=self.dataDisplay,value=data['Type'],choices=typeChoices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    typeSel.Bind(wx.EVT_COMBOBOX, OnSelType)
    zoneSizer.Add(typeSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add((10,0),0)    
    mainSizer.Add(zoneSizer,1,wx.EXPAND|wx.RIGHT)
        
    izone = zones.index(data['Zone'])
    layerSizer = wx.BoxSizer(wx.HORIZONTAL)
    layerSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Layer'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    layerSel = wx.Slider(parent=self.dataDisplay,maxValue=HKLmax[izone],minValue=HKLmin[izone],
        style=wx.SL_HORIZONTAL|wx.SL_AUTOTICKS|wx.SL_LABELS,value=0)
    layerSizer.Add(layerSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    layerSel.Bind(wx.EVT_SLIDER, OnLayerSlider)    
    layerSizer.Add((10,0),0)    
    mainSizer.Add(layerSizer,1,wx.EXPAND|wx.RIGHT)

        
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataDisplay.SetSize(mainSizer.Fit(self.dataFrame))
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
        
def UpdateImageControls(self,data):
    
    def OnNewColorBar(event):
        data['color'] = colSel.GetValue()
        self.PlotImage()
        
    def OnNewCalibrant(event):
        data['calibrant'] = calSel.GetValue()
        
    def OnMaxSlider(event):
        imax = max(data['range'][1][0],int(maxSel.GetValue()))
        maxSel.SetValue(imax)
        data['range'][1][1] = imax
        self.PlotImage()
        
    def OnMinSlider(event):
        imin = min(data['range'][1][1],int(minSel.GetValue()))
        minSel.SetValue(imin)
        data['range'][1][0] = imin
        self.PlotImage()
        
    def OnNumOutChans(event):
        try:
            numChans = int(outChan.GetValue())
            data['outChannels'] = numChans
        except ValueError:
            pass
        outChan.SetValue(str(data['outChannels']))          #reset in case of error        
        
    def OnWavelength(event):
        try:
            wave = float(waveSel.GetValue())
            data['wavelength'] = wave
        except ValueError:
            pass
        waveSel.SetValue("%6.5f" % (data['wavelength']))          #reset in case of error          
        
    def OnDistance(event):
        try:
            dist = float(distSel.GetValue())
            data['distance'] = dist
        except ValueError:
            pass
        distSel.SetValue("%8.3f"%(data['distance']))          #reset in case of error  
        
    def OnImRefine(event):
        ImageCalibRef[0] = centRef.GetValue()
        ImageCalibRef[1] = waveRef.GetValue()
        ImageCalibRef[2] = distRef.GetValue()
        ImageCalibRef[3] = tiltRef.GetValue()
        ImageCalibRef[4] = rotRef.GetValue()
        SetStatusLine()
                
    def OnShowLines(event):
        if data['showLines']:
            data['showLines'] = False
        else:
            data['showLines'] = True
        self.PlotImage()
        
    def OnFullIntegrate(event):
        if data['fullIntegrate']:
            data['fullIntegrate'] = False
            data['LRazimuth'] = [-45,45]
            self.LRazim.SetValue("%6d,%6d" % (-45,45))            
        else:
            data['fullIntegrate'] = True
            data['LRazimuth'] = [0,360]
            self.LRazim.SetValue("%6d,%6d" % (0,360))            
        self.PlotImage()
        
        
    def OnSetDefault(event):
        import copy
        if data['setDefault']:
            self.imageDefault = {}
            data['setDefault'] = False
        else:
            self.imageDefault = copy.copy(data)
            data['setDefault'] = True
        
    
    def OnCalibrate(event):
        G2cmp.ImageCalibrate(data,self.ImageZ,self.PlotImage())
        
    def OnIntegrate(event):
        print 'image integrate'
        
    def SetStatusLine():
        if data['refine'][0]:
            Status.SetStatusText("On Image: key 'c' to mark center or 'r' on inner ring for calibration")
        else:
            Status.SetStatusText("On Image: key 'r' on inner ring for calibration")
                              
    colorList = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
    import ImageCalibrants as calFile
    calList = [m for m in calFile.Calibrants.keys()]
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.ImageMenu)
    Status = self.dataFrame.CreateStatusBar()
    SetStatusLine()
    self.dataFrame.Bind(wx.EVT_MENU, OnCalibrate, id=wxID_IMCALIBRATE)
    self.dataFrame.Bind(wx.EVT_MENU, OnIntegrate, id=wxID_IMINTEGRATE)        
    self.dataDisplay = wx.Panel(self.dataFrame)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    maxSizer = wx.FlexGridSizer(2,2,0,5)
    maxSizer.AddGrowableCol(1,1)
    maxSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Max intensity'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    maxSel = wx.Slider(parent=self.dataDisplay,maxValue=data['range'][0][1],minValue=data['range'][0][0],
        style=wx.SL_HORIZONTAL,value=data['range'][1][1])
    maxSizer.Add(maxSel,1,wx.EXPAND|wx.RIGHT)
    maxSel.Bind(wx.EVT_SLIDER, OnMaxSlider)    
    maxSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Min intensity'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    minSel = wx.Slider(parent=self.dataDisplay,maxValue=data['range'][0][1],minValue=data['range'][0][0],
        style=wx.SL_HORIZONTAL,value=data['range'][1][0])
    maxSizer.Add(minSel,1,wx.EXPAND|wx.RIGHT)
    minSel.Bind(wx.EVT_SLIDER, OnMinSlider)
    mainSizer.Add(maxSizer,1,wx.EXPAND|wx.RIGHT)
    
    comboSizer = wx.FlexGridSizer(2,4,5,5)
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Color bar '),0,
        wx.ALIGN_CENTER_VERTICAL)
    colSel = wx.ComboBox(parent=self.dataDisplay,value=data['color'],choices=colorList,
        style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
    colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
    comboSizer.Add(colSel,0,wx.ALIGN_CENTER_VERTICAL)
    
    comboSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calibrant '),0,
        wx.ALIGN_CENTER_VERTICAL)
    calSel = wx.ComboBox(parent=self.dataDisplay,value=data['calibrant'],choices=calList,
        style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
    calSel.Bind(wx.EVT_COMBOBOX, OnNewCalibrant)
    comboSizer.Add(calSel,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(comboSizer,0,wx.ALIGN_CENTER_HORIZONTAL)
    mainSizer.Add((5,5),0)
         
    dataSizer = wx.FlexGridSizer(6,5,5,5)
    dataSizer.Add((5,0),0)
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Calibration coefficients'),0,
        wx.ALIGN_CENTER_VERTICAL)    
    dataSizer.Add((5,0),0)
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Integration coefficients'),0,
        wx.ALIGN_CENTER_VERTICAL)    
    dataSizer.Add((5,0),0)
    
    ImageCalibRef = data['refine']
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Beam center X,Y'),0,
        wx.ALIGN_CENTER_VERTICAL)
    cent = data['center']
    self.centText = wx.TextCtrl(parent=self.dataDisplay,value=("%8.3f,%8.3f" % (cent[0],cent[1])),style=wx.TE_READONLY)
    dataSizer.Add(self.centText,0,wx.ALIGN_CENTER_VERTICAL)
    centRef = wx.CheckBox(parent=self.dataDisplay,label='refine?')
    centRef.Bind(wx.EVT_CHECKBOX, OnImRefine)
    centRef.SetValue(ImageCalibRef[0])
    dataSizer.Add(centRef,0,wx.ALIGN_CENTER_VERTICAL)
    
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Inner/Outer radii'),0,
        wx.ALIGN_CENTER_VERTICAL)
    IOradii = data['IOradii']
    self.IOradText = wx.TextCtrl(parent=self.dataDisplay,
        value=("%8.3f,%8.3f" % (IOradii[0],IOradii[1])),style=wx.TE_READONLY)
    dataSizer.Add(self.IOradText,0,wx.ALIGN_CENTER_VERTICAL)
       
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Wavelength'),0,
        wx.ALIGN_CENTER_VERTICAL)
    waveSel = wx.TextCtrl(parent=self.dataDisplay,value=("%6.5f" % (data['wavelength'])),
        style=wx.TE_PROCESS_ENTER)
    waveSel.Bind(wx.EVT_TEXT_ENTER,OnWavelength)
    dataSizer.Add(waveSel,0,wx.ALIGN_CENTER_VERTICAL)
    waveRef = wx.CheckBox(parent=self.dataDisplay,label='refine?')
    waveRef.Bind(wx.EVT_CHECKBOX, OnImRefine)
    waveRef.SetValue(ImageCalibRef[1])
    dataSizer.Add(waveRef,0,wx.ALIGN_CENTER_VERTICAL)
         
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Start/End azimuth'),0,
        wx.ALIGN_CENTER_VERTICAL)
    LRazim = data['LRazimuth']
    self.LRazim = wx.TextCtrl(parent=self.dataDisplay,
        value=("%6d,%6d" % (LRazim[0],LRazim[1])),style=wx.TE_READONLY)
    dataSizer.Add(self.LRazim,0,wx.ALIGN_CENTER_VERTICAL)
       
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Distance'),0,
        wx.ALIGN_CENTER_VERTICAL)
    distSel = wx.TextCtrl(parent=self.dataDisplay,value=("%8.3f"%(data['distance'])),style=wx.TE_PROCESS_ENTER)
    distSel.Bind(wx.EVT_TEXT_ENTER,OnDistance)
    dataSizer.Add(distSel,0,wx.ALIGN_CENTER_VERTICAL)
    distRef = wx.CheckBox(parent=self.dataDisplay,label='refine?')
    distRef.Bind(wx.EVT_CHECKBOX, OnImRefine)
    distRef.SetValue(ImageCalibRef[2])
    dataSizer.Add(distRef,0,wx.ALIGN_CENTER_VERTICAL)

    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' No. bins'),0,
        wx.ALIGN_CENTER_VERTICAL)
    outChan = wx.TextCtrl(parent=self.dataDisplay,value=str(data['outChannels']),style=wx.TE_PROCESS_ENTER)
    outChan.Bind(wx.EVT_TEXT_ENTER,OnNumOutChans)
    dataSizer.Add(outChan,0,wx.ALIGN_CENTER_VERTICAL)

    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Tilt angle'),0,
        wx.ALIGN_CENTER_VERTICAL)
    self.tiltSel = wx.TextCtrl(parent=self.dataDisplay,value=("%9.3f"%(data['tilt'])),style=wx.TE_READONLY)
    dataSizer.Add(self.tiltSel,0,wx.ALIGN_CENTER_VERTICAL)
    tiltRef = wx.CheckBox(parent=self.dataDisplay,label='refine?')
    tiltRef.Bind(wx.EVT_CHECKBOX, OnImRefine)
    tiltRef.SetValue(ImageCalibRef[3])
    dataSizer.Add(tiltRef,0,wx.ALIGN_CENTER_VERTICAL)
    showLines = wx.CheckBox(parent=self.dataDisplay,label='Show integration limits?')
    dataSizer.Add(showLines,0)
    showLines.Bind(wx.EVT_CHECKBOX, OnShowLines)
    showLines.SetValue(data['showLines'])
    fullIntegrate = wx.CheckBox(parent=self.dataDisplay,label='Do full integration?')
    dataSizer.Add(fullIntegrate,0)
    fullIntegrate.Bind(wx.EVT_CHECKBOX, OnFullIntegrate)
    fullIntegrate.SetValue(data['fullIntegrate'])
    
    dataSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Tilt rotation'),0,
        wx.ALIGN_CENTER_VERTICAL)
    self.rotSel = wx.TextCtrl(parent=self.dataDisplay,value=("%9.3f"%(data['rotation'])),style=wx.TE_READONLY)
    dataSizer.Add(self.rotSel,0,wx.ALIGN_CENTER_VERTICAL)
    rotRef = wx.CheckBox(parent=self.dataDisplay,label='refine?')
    rotRef.Bind(wx.EVT_CHECKBOX, OnImRefine)
    rotRef.SetValue(ImageCalibRef[4])
    dataSizer.Add(rotRef,0,)
    setDefault = wx.CheckBox(parent=self.dataDisplay,label='Use as default for all images?')
    dataSizer.Add(setDefault,0)
    setDefault.Bind(wx.EVT_CHECKBOX, OnSetDefault)
    setDefault.SetValue(data['setDefault'])
    dataSizer.Add((5,5),0)
        
    mainSizer.Add(dataSizer,0)
    
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataDisplay.SetSize(mainSizer.Fit(self.dataFrame))
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
    
def UpdatePhaseData(self,item,data,oldPage):
    import GSASIIElem as G2el
    Atoms = []
    self.SelectedRow = 0
    
    def BookResize(event):
        w,h = self.GetSize()
        self.dataDisplay.SetSize(wx.Size(w,h))
        
    def FillGeneralGrid():
        def SetLatticeParametersStyle(SGData,table):
            if SGData['SGLaue'] in ['m3','m3m']:
                table[4][2] = table[4][3] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                General.SetCellStyle(4,3,"light grey",True)
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['3R','3mR']:
                table[4][2] = table[4][3] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                General.SetCellStyle(4,3,"light grey",True)
                table[4][5] = table[4][6] = table[4][4]
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
                table[4][2] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                table[4][4] = table[4][5] = 90.
                table[4][6] = 120.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['4/m','4/mmm']:
                table[4][2] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['mmm']:
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['2/m']:
                if SGData['SGUniq'] == 'a':
                    table[4][5]= table[4][6] = 90.
                    General.SetCellStyle(4,5,"light grey",True)
                    General.SetCellStyle(4,6,"light grey",True)
                if SGData['SGUniq'] == 'b':
                    table[4][4]= table[4][6] = 90.
                    General.SetCellStyle(4,4,"light grey",True)
                    General.SetCellStyle(4,6,"light grey",True)
                if SGData['SGUniq'] == 'c':
                    table[4][4]= table[4][5] = 90.
                    General.SetCellStyle(4,4,"light grey",True)
                    General.SetCellStyle(4,5,"light grey",True)
            
        def RefreshGeneralGrid(event):
                
            r,c =  event.GetRow(),event.GetCol()
            generalData[0] = table[0][0]
            self.PatternTree.SetItemText(item,generalData[0])
            generalData[1] = table[1][0]
            SpcGp = table[2][0]
            SGErr,SGData = G2spc.SpcGroup(SpcGp)
            if r == 2 and c == 0:
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    table[2][0] = generalData[2]['SpGrp']
                    msg = 'Space Group Error'
                    Style = wx.ICON_EXCLAMATION
                else:
                    text = G2spc.SGPrint(SGData)
                    generalData[2] = SGData
                    msg = 'Space Group Information'
                    Style = wx.ICON_INFORMATION
                Text = ''
                for line in text:
                    Text += line+'\n'
                wx.MessageBox(Text,caption=msg,style=Style)
            General.SetCellValue(4,0,str(generalData[3][0]))
            for c in range(1,7):
                General.SetCellStyle(4,c,"white",False)
                generalData[3][c] = float(General.GetCellValue(4,c))
            generalData[3][7] = G2cmp.calc_V(G2cmp.cell2A(generalData[3][1:7]))
            SetLatticeParametersStyle(SGData,table)
            generalData[4][1] = float(General.GetCellValue(5,1))
            print 'did general refresh'
            General.ForceRefresh()
                        
        rowLabels = ['Phase name','Phase type','Space group',
            'Lattice ',' parameters','Scale factor','Elements','No. per cell','Atom weight','','Bond radii','Angle radii']
        generalData = data['General']
        atomData = data['Atoms']
        AtomTypes = []
        NoAtoms = {}
        BondRadii = []
        AngleRadii = []
        AtomMass = []
        colType = 1
        colSS = 7
        self.dataFrame.setSizePosLeft([600,350])
        if generalData[1] =='macromolecular':
            colType = 4
            colSS = 10
        for atom in atomData:
            if AtomTypes.count(atom[colType]):
                NoAtoms[atom[colType]] += atom[colSS-1]*atom[colSS+1]
            else:
                Info = G2el.GetAtomInfo(atom[colType])
                AtomTypes.append(Info['Symbol'])
                BondRadii.append(Info['Drad'])
                AngleRadii.append(Info['Arad'])
                AtomMass.append(Info['Mass'])
                NoAtoms[atom[colType]] = atom[colSS-1]*atom[colSS+1]
        generalData[5:9] = [AtomTypes,NoAtoms,AtomMass,BondRadii,AngleRadii]
        colLabels = []
        colLabels += ['' for i in range(max(8,len(generalData[5])))]
        table = []
        table.append([generalData[0],'','','','','','','',''])      #phase name
        table.append([generalData[1],'','','','','','','',''])      #phase type
        E,SGData = G2spc.SpcGroup(generalData[2]['SpGrp'])
        table.append([SGData['SpGrp'],'','','','','','','',''])     #space group symbol
        table.append(['refine','a    ','b    ','c    ','alpha ','beta ','gamma','volume  '])
        table.append(generalData[3])                      #lattice parameters
        table.append([generalData[4][0],generalData[4][1],'','','','','',''])   #scale factor
        table.append(generalData[5]+['' for i in range(max(8,len(generalData[5])))]) #element list
        line = []
        mass = 0.
        for i,elem in enumerate(generalData[5]):
            mass += generalData[6][elem]*generalData[7][i]
            line.append(generalData[6][elem])
        Volume = generalData[3][7]
        table.append(line+['' for i in range(max(8,len(generalData[5])))]) #No. per cell
        table.append(generalData[7]+['' for i in range(max(8,len(generalData[5])))])  #At. wt.
        if generalData[1] == 'macromolecular' and mass > 0.0:
            table.append(['density',mass/(0.6022137*Volume),'Matthews coeff.',Volume/mass,'','','','',''])            
        else:
            table.append(['density',mass/(0.6022137*Volume),'','','','','','',''])
        table.append(generalData[8]+['' for i in range(max(8,len(generalData[5])))])
        table.append(generalData[9]+['' for i in range(max(8,len(generalData[5])))])
        Types = [wg.GRID_VALUE_STRING for i in range(max(8,len(generalData[5])))]
        generalTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        General.SetTable(generalTable, True)
        General.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshGeneralGrid)
        General.SetMargins(0,0)
        General.SetColSize(0,100)
        General.SetColLabelSize(0)
        for c in range(max(8,len(generalData[5]))):
            if c > 0:
                General.SetReadOnly(0,c,isReadOnly=True)
                General.SetReadOnly(1,c,isReadOnly=True)
                General.SetReadOnly(2,c,isReadOnly=True)
            General.SetReadOnly(3,c,isReadOnly=True)                         #unit cell labels
            General.SetCellAlignment(3,c,wx.ALIGN_RIGHT, wx.ALIGN_CENTRE)
            if c < 4:
                General.SetCellRenderer(4,c,wg.GridCellFloatRenderer(10,5))
                General.SetCellEditor(4,c,wg.GridCellFloatEditor(10,5))
                General.SetReadOnly(9,c,isReadOnly=True)
            else:
                General.SetCellRenderer(4,c,wg.GridCellFloatRenderer(10,3))
                General.SetCellEditor(4,c,wg.GridCellFloatEditor(10,3))
            for r in range(6,12):
                General.SetReadOnly(r,c,isReadOnly=True)
        General.SetReadOnly(4,7,isReadOnly=True)                            #cell volume - no edit
        General.SetCellEditor(1,0,wg.GridCellChoiceEditor(['nuclear','modulated',   #phase type 
            'magnetic','macromolecular','Pawley'],False))                           #- change only if no atoms
        if line:                                                    #no.of atoms not zero!
            General.SetReadOnly(1,0,isReadOnly=True)                #can't change phase type
        General.SetCellRenderer(4,0,wg.GridCellBoolRenderer())              #lattice parameters            
        General.SetCellEditor(4,0,wg.GridCellBoolEditor())
        SetLatticeParametersStyle(SGData,table)
        General.SetCellRenderer(5,1,wg.GridCellFloatRenderer(10,4))         #scale factor
        General.SetCellEditor(5,1,wg.GridCellFloatEditor(10,4))
        General.SetCellRenderer(5,0,wg.GridCellBoolRenderer())            
        General.SetCellEditor(5,0,wg.GridCellBoolEditor())
        General.SetCellRenderer(9,1,wg.GridCellFloatRenderer(8,3))
        General.SetCellRenderer(9,3,wg.GridCellFloatRenderer(8,3))
    
    def FillAtomsGrid():
        
        def RefreshAtomGrid(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0:                          #on col label!
                sel = -1
                if Atoms.GetColLabelValue(c) == 'refine':
                    choice = ['F - site fraction','X - coordinates','U - thermal parameters']
                    dlg = wx.MultiChoiceDialog(self,'Select','Refinement controls',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelections()
                        parms = ''
                        for x in sel:
                            parms += choice[x][0]                            
                elif Atoms.GetColLabelValue(c) == 'I/A':
                    choice = ['Isotropic','Anisotropic']
                    dlg = wx.SingleChoiceDialog(self,'Select','Thermal Motion',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel][0]
                if sel >= 0:
                    for r in range(Atoms.GetNumberRows()):
                        Atoms.SetCellValue(r,c,parms)
            elif c < 0:                    #picked atom row
                self.SelectedRow = r
            elif Atoms.GetColLabelValue(c) in ['x','y','z']:
                colLabel = Atoms.GetColLabelValue(c)
                if colLabel == 'x':
                    XYZ = [atomData[r][c],atomData[r][c+1],atomData[r][c+2]]
                elif colLabel == 'y':
                    XYZ = [atomData[r][c-1],atomData[r][c],atomData[r][c+1]]
                elif colLabel == 'z':
                    XYZ = [atomData[r][c-2],atomData[r][c-1],atomData[r][c]]
                if None in XYZ:
                    XYZ = [0,0,0]
                SScol = colLabels.index('site sym')
                Mulcol = colLabels.index('mult')
                E,SGData = G2spc.SpcGroup(generalData[2]['SpGrp'])
                Sytsym,Mult = G2spc.SytSym(XYZ,SGData)
                atomData[r][SScol] = Sytsym
                atomData[r][Mulcol] = Mult
                Atoms.ForceRefresh()
                    
        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if Atoms.GetColLabelValue(c) == 'Type':
                PE = G2elem.PickElement(self)
                if PE.ShowModal() == wx.ID_OK:
                    atomData[r][c] = PE.Elem.strip()
                PE.Destroy()
                Atoms.ForceRefresh()
            else:
                event.Skip()
        
        generalData = data['General']
        atomData = data['Atoms']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,X,XU,U,F,FX,FXU,FU",
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_STRING,wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_CHOICE+":I,A",
            wg.GRID_VALUE_FLOAT+':10,4',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
        colLabels = ['Name','Type','refine','x','y','z','frac','site sym','mult','I/A','Uiso','U11','U22','U33','U12','U13','U23']
        if generalData[1] == 'magnetic':
            colLabels += ['Mx','My','Mz']
            Types[2] = wg.GRID_VALUE_CHOICE+": ,X,XU,U,M,MX,MXU,MU,F,FX,FXU,FU,FM,FMX,FMU,"
            Types += [
                wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
        elif generalData[1] == 'macromolecular':
            colLabels = ['res no','residue','chain'] + colLabels
            Types = [wg.GRID_VALUE_NUMBER,
                wg.GRID_VALUE_CHOICE+": ,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE,HOH,UNK",
                wg.GRID_VALUE_STRING] + Types        
        table = []
        rowLabels = []
        for i,atom in enumerate(atomData):
            table.append(atom)
            rowLabels.append(str(i+1))
        atomTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        Atoms.SetTable(atomTable, True)
        Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_SELECT_CELL, AtomTypeSelect)
        Atoms.SetMargins(0,0)
        Atoms.AutoSizeColumns(True)
        colType = colLabels.index('Type')
        colSS = colLabels.index('site sym')
        colIA = colLabels.index('I/A')
        for row in range(Atoms.GetNumberRows()):
            Atoms.SetReadOnly(row,colSS,True)                         #site sym
            Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            if Atoms.GetCellValue(row,colIA) == 'I':
                for i in range(2,8):
                    Atoms.SetCellRenderer(row,colIA+i,wg.GridCellStringRenderer())
                    Atoms.SetReadOnly(row,colIA+i,isReadOnly=True)
                    Atoms.SetCellValue(row,colIA+i,'')
            elif Atoms.GetCellValue(row,colIA) == 'A':
                Atoms.SetCellRenderer(row,colIA+1,wg.GridCellStringRenderer())
                Atoms.SetReadOnly(row,colIA+1,isReadOnly=True)
                Atoms.SetCellValue(row,colIA+1,'')
        
    def AtomAdd(event):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        if generalData[1] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','UNK','',0,0,0,0,'',0,'I',0.10,0,0,0,0,0,0])
        elif generalData[1] == 'nuclear':
            atomData.append(['UNK','UNK','',0,0,0,0,'',0,'I',0.01,0,0,0,0,0,0])
        event.StopPropagation()
        FillAtomsGrid()
            
    def AtomInsert(event):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        if generalData[1][0] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','UNK','',0,0,0,0,'',0,'I',0.10,0,0,0,0,0,0])
        elif generalData[1][0] == 'nuclear':
            atomData.append(['UNK','UNK','',0,0,0,0,'',0,'I',0.01,0,0,0,0,0,0])
        event.StopPropagation()
        FillAtomsGrid()
        

    def UpdateDrawing():
        print 'Drawing'
        
    def FillPawleyReflectionsGrid():
        
        print 'Pawley reflections'
        
    def OnPageChanged(event):
        page = event.GetSelection()
        text = self.dataDisplay.GetPageText(page)
        if text == 'Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.AtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, AtomAdd, id=wxID_ATOMSEDITADD)
            self.dataFrame.Bind(wx.EVT_MENU, AtomInsert, id=wxID_ATOMSEDITINSERT)
            FillAtomsGrid()            
        else:
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        event.Skip()
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()                    
    PhaseName = self.PatternTree.GetItemText(item)
    self.dataFrame.SetLabel('Phase Data for '+PhaseName)
    self.dataDisplay = GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())
    
    General = GSGrid(parent=self.dataDisplay)
    FillGeneralGrid()
    self.dataDisplay.AddPage(General,'General')
     
    GeneralData = data['General']
    if GeneralData[3] == 'Pawley':
        PawleyRefl = GSGrid(parent=self.dataDisplay)
        self.dataDisplay.AddPage(PawleyRefl,'Pawley reflections')
        FillPawleyReflectionsGrid()
    else:
        Atoms = GSGrid(parent=self.dataDisplay)
        FillAtomsGrid()
        self.dataDisplay.AddPage(Atoms,'Atoms')

    Drawing = wx.Window(parent=self.dataDisplay)
    self.dataDisplay.AddPage(Drawing,'Drawing')
   
    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    
                          
def GetPatternTreeItemId(self, parentId, itemText):
    item, cookie = self.PatternTree.GetFirstChild(parentId)
    while item:
        if self.PatternTree.GetItemText(item) == itemText:
            return item
        item, cookie = self.PatternTree.GetNextChild(parentId, cookie)
    return 0                

def MovePatternTreeToGrid(self,item):
    
    oldPage = 0
    if self.dataFrame:
        if self.dataFrame.GetLabel() == 'Comments':
            data = [self.dataDisplay.GetValue()]
            self.dataDisplay.Clear() 
            Id = GetPatternTreeItemId(self,self.root, 'Comments')
            if Id: self.PatternTree.SetItemPyData(Id,data)
        if self.dataFrame.GetLabel() == 'Notebook':
            data = [self.dataDisplay.GetValue()]
            self.dataDisplay.Clear() 
            Id = GetPatternTreeItemId(self,self.root, 'Notebook')
            if Id: self.PatternTree.SetItemPyData(Id,data)
        if 'Phase Data for' in self.dataFrame.GetLabel():
            if self.dataDisplay: 
                oldPage = self.dataDisplay.GetSelection()
        self.dataFrame.Clear()
        self.dataFrame.SetLabel('')
    else:
       self.dataFrame = DataFrame(parent=self.mainPanel)
            
    self.PickId = 0
    self.PatternId = 0
    self.PeakFit.Enable(False)
    self.AutoPeakFit.Enable(False)
    self.IndexPeaks.Enable(False)
    self.RefineCell.Enable(False)
    parentID = self.root
    if item != self.root:
        parentID = self.PatternTree.GetItemParent(item)
    if self.PatternTree.GetItemParent(item) == self.root:
        self.PatternId = item
        self.ExportPattern.Enable(True)
        self.PickId = item
        if self.PatternTree.GetItemText(item) == 'Notebook':
            self.PatternId = 0
            self.ExportPattern.Enable(False)
            data = self.PatternTree.GetItemPyData(item)
            UpdateNotebook(self,data)
        elif self.PatternTree.GetItemText(item) == 'Controls':
            self.PatternId = 0
            self.ExportPattern.Enable(False)
            data = self.PatternTree.GetItemPyData(item)
            UpdateControls(self,data)
        elif 'IMG' in self.PatternTree.GetItemText(item):
            self.Image = item
            self.Img = 0
            self.PlotImage()
        elif 'PWDR' in self.PatternTree.GetItemText(item):
            self.PlotPatterns()
        elif 'SXTL ' in self.PatternTree.GetItemText(item):
            self.Sngl = item
            self.PlotSngl()
            
    elif self.PatternTree.GetItemText(parentID) == 'Phases':
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdatePhaseData(self,item,data,oldPage)
    elif self.PatternTree.GetItemText(item) == 'Comments':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateComments(self,data)
    elif self.PatternTree.GetItemText(item) == 'Image Controls':
        self.dataFrame.SetTitle('Image Controls')
        self.PickId = item
        self.Image = self.PatternTree.GetItemParent(item)
        data = self.PatternTree.GetItemPyData(item)
        UpdateImageControls(self,data)
        self.PlotImage()
    elif self.PatternTree.GetItemText(item) == 'HKL Plot Controls':
        self.PickId = item
        self.Sngl = self.PatternTree.GetItemParent(item)
        data = self.PatternTree.GetItemPyData(item)
        UpdateHKLControls(self,data)
        self.PlotSngl()               
    elif self.PatternTree.GetItemText(item) == 'Peak List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PeakFit.Enable(True)
        self.ExportPeakList.Enable(True)
        self.AutoPeakFit.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdatePeakGrid(self,data)
        self.PlotPatterns()
    elif self.PatternTree.GetItemText(item) == 'Background':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateBackgroundGrid(self,data)
        self.PlotPatterns()
    elif self.PatternTree.GetItemText(item) == 'Limits':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateLimitsGrid(self,data)
        self.PlotPatterns()
    elif self.PatternTree.GetItemText(item) == 'Instrument Parameters':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateInstrumentGrid(self,data)
        self.PlotPeakWidths()
    elif self.PatternTree.GetItemText(item) == 'Index Peak List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.ExportPeakList.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateIndexPeaksGrid(self,data)
        self.NewPlot = True
        self.PlotPatterns()
    elif self.PatternTree.GetItemText(item) == 'Unit Cells List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.IndexPeaks.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        if not data:
            data.append([0,0.1,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, max zero error, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            self.PatternTree.SetItemPyData(item,data)                             
        UpdateUnitCellsGrid(self,data)
        self.NewPlot = True
        self.PlotPatterns()
