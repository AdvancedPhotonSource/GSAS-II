#GSASII - data display routines
import wx
import wx.grid as wg
import time
import cPickle
import GSASIIpath
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, 
] = [wx.NewId() for _init_coll_Atom_Items in range(2)]

[ wxID_IMCALIBRATE, wxID_IMINTEGRATE, wxID_IMCLEARCALIB, wxID_SAVEINTG, 
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL,
] = [wx.NewId() for _init_coll_IMAGE_Items in range(6)]

[ wxID_MASKCOPY,
] = [wx.NewId() for _init_coll_MASK_Items in range(1)]

[ wxID_INSTPRMRESET,
] = [wx.NewId() for _init_coll_INST_Items in range(1)]

[ wxID_INDXRELOAD,
] = [wx.NewId() for _init_coll_IndPeaks_Items in range(1)]

[ wxID_UNDO,wxID_PEAKFIT,wxID_AUTOPEAKFIT,
] = [wx.NewId() for _init_coll_PEAK_Items in range(3)]

[  wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL,
] = [wx.NewId() for _init_coll_INDEX_Items in range(3)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

class DataFrame(wx.Frame):
    def _init_coll_BlankMenu(self,parent):
        parent.Append(menu=self.Blank,title='')
        
    def _init_coll_AtomsMenu(self,parent):
        parent.Append(menu=self.AtomEdit, title='Add atom')

    def _init_coll_IndPeaksMenu(self,parent):
        parent.Append(menu=self.IndPeaksEdit,title='Index Peaks Operations')
                   
    def _init_coll_ImageMenu(self,parent):
        parent.Append(menu=self.ImageEdit, title='Image Operations')
        
    def _init_coll_InstMenu(self,parent):
        parent.Append(menu=self.InstEdit, title='Inst. Parm. Operations')
        
    def _init_coll_MaskMenu(self,parent):
        parent.Append(menu=self.MaskEdit, title='Mask Operations')
        
    def _init_coll_PeakMenu(self,parent):
        parent.Append(menu=self.PeakEdit, title='Peak Fitting')

    def _init_coll_IndexMenu(self,parent):
        parent.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        
    def _init_coll_Atom_Items(self,parent):
        parent.Append(help='',id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append empty atom')
        parent.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert empty atom',
            help='Double left click on atom row to Insert before')
            
    def _init_coll_IndPeaks_Items(self,parent):
        parent.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
            
    def _init_coll_Image_Items(self,parent):
        parent.Append(help='Calibrate detector by fitting to calibrant lines', 
            id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        parent.Append(help='Clear calibration data points and rings',id=wxID_IMCLEARCALIB, 
            kind=wx.ITEM_NORMAL,text='Clear calibration')
        parent.Append(help='Integrate selected image',id=wxID_IMINTEGRATE, 
            kind=wx.ITEM_NORMAL,text='Integrate')
        parent.Append(help='Integrate all images selected from list',id=wxID_INTEGRATEALL,
            kind=wx.ITEM_NORMAL,text='Integrate all')
        parent.Append(help='Save integration results as a series of 1-D powder patterns', 
            id=wxID_SAVEINTG, kind=wx.ITEM_NORMAL,text='Save Integration')
        parent.Append(help='Copy image controls to other images', 
            id=wxID_IMCOPYCONTROLS, kind=wx.ITEM_NORMAL,text='Copy Controls')
                    
    def _init_coll_Mask_Items(self,parent):
        parent.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')

    def _init_coll_Inst_Items(self,parent):
        parent.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')

    def _init_coll_Peak_Items(self,parent):
        self.UnDo = parent.Append(help='Undo last least squares refinement', 
            id=wxID_UNDO, kind=wx.ITEM_NORMAL,text='UnDo')
        self.PeakFit = parent.Append(id=wxID_PEAKFIT, kind=wx.ITEM_NORMAL,text='PeakFit', 
            help='Do single cycle of peak fitting least-squares refinement' )
        self.AutoPeakFit = parent.Append(id=wxID_AUTOPEAKFIT, kind=wx.ITEM_NORMAL, 
            text='AutoPeakFit',help='Do peak fitting least-squares to convergence' )
            
    def _init_coll_Index_Items(self,parent):
        self.IndexPeaks = parent.Append(help='', id=wxID_INDEXPEAKS, kind=wx.ITEM_NORMAL,
            text='Index Cell')
        self.CopyCell = parent.Append( id=wxID_COPYCELL, kind=wx.ITEM_NORMAL,text='Copy Cell', 
            help='Copy selected unit cell from indexing to cell refinement fields')
        self.RefineCell = parent.Append( id=wxID_REFINECELL, kind=wx.ITEM_NORMAL, 
            text='Refine Cell',help='Refine unit cell parameters from indexed peaks')

    def _init_utils(self):
        self.BlankMenu = wx.MenuBar()
        
        self.AtomsMenu = wx.MenuBar()
        self.ImageMenu = wx.MenuBar()
        self.MaskMenu = wx.MenuBar()
        self.InstMenu = wx.MenuBar()
        self.PeakMenu = wx.MenuBar()
        self.IndPeaksMenu = wx.MenuBar()
        self.IndexMenu = wx.MenuBar()
        self.AtomEdit = wx.Menu(title='')
        self.ImageEdit = wx.Menu(title='')
        self.MaskEdit = wx.Menu(title='')
        self.InstEdit = wx.Menu(title='')
        self.PeakEdit = wx.Menu(title='')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndexEdit = wx.Menu(title='')
        self._init_coll_AtomsMenu(self.AtomsMenu)
        self._init_coll_Atom_Items(self.AtomEdit)
        self._init_coll_ImageMenu(self.ImageMenu)
        self._init_coll_Image_Items(self.ImageEdit)
        self._init_coll_MaskMenu(self.MaskMenu)
        self._init_coll_Mask_Items(self.MaskEdit)
        self._init_coll_InstMenu(self.InstMenu)
        self._init_coll_Inst_Items(self.InstEdit)
        self._init_coll_PeakMenu(self.PeakMenu)
        self._init_coll_Peak_Items(self.PeakEdit)
        self._init_coll_IndPeaksMenu(self.IndPeaksMenu)
        self._init_coll_IndPeaks_Items(self.IndPeaksEdit)
        self._init_coll_IndexMenu(self.IndexMenu)
        self._init_coll_Index_Items(self.IndexEdit)
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.AutoPeakFit.Enable(False)
        self.IndexPeaks.Enable(False)
        self.CopyCell.Enable(False)
        self.RefineCell.Enable(False)
               
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
             
def UpdateHKLControls(self,data):
    
    def OnScaleSlider(event):
        scale = int(scaleSel.GetValue())/1000.
        scaleSel.SetValue(int(scale*1000.))
        data['Scale'] = scale*10.
        G2plt.PlotSngl(self)
        
    def OnLayerSlider(event):
        layer = layerSel.GetValue()
        data['Layer'] = layer
        G2plt.PlotSngl(self)
        
    def OnSelZone(event):
        data['Zone'] = zoneSel.GetValue()
        G2plt.PlotSngl(self,newPlot=True)
        
    def OnSelType(event):
        data['Type'] = typeSel.GetValue()
        G2plt.PlotSngl(self)
        
    def SetStatusLine():
        Status.SetStatusText("look at me!!!")
                                      
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    if not self.dataFrame.GetStatusBar():
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
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleSizer.Add(wx.StaticText(parent=self.dataDisplay,label=' Scale'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    scaleSel = wx.Slider(parent=self.dataDisplay,maxValue=1000,minValue=100,
        style=wx.SL_HORIZONTAL,value=int(data['Scale']*100))
    scaleSizer.Add(scaleSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    scaleSel.SetLineSize(100)
    scaleSel.SetPageSize(900)
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
    layerSel.SetLineSize(1)
    layerSel.SetLineSize(5)
    layerSel.Bind(wx.EVT_SLIDER, OnLayerSlider)    
    layerSizer.Add(layerSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    layerSizer.Add((10,0),0)    
    mainSizer.Add(layerSizer,1,wx.EXPAND|wx.RIGHT)

        
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataDisplay.SetSize(mainSizer.Fit(self.dataFrame))
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
        
                          
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
        self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
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

    self.dataFrame.Raise()            
    self.PickId = 0
    parentID = self.root
    self.ExportPattern.Enable(False)
    if item != self.root:
        parentID = self.PatternTree.GetItemParent(item)
    if self.PatternTree.GetItemParent(item) == self.root:
        self.PatternId = item
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
            G2plt.PlotImage(self,newPlot=True)
        elif 'PKS' in self.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(self)
        elif 'PWDR' in self.PatternTree.GetItemText(item):            
            self.ExportPattern.Enable(True)
            G2plt.PlotPatterns(self,newPlot=True)
        elif 'SXTL' in self.PatternTree.GetItemText(item):
            self.Sngl = item
            G2plt.PlotSngl(self,newPlot=True)
            
    elif self.PatternTree.GetItemText(parentID) == 'Phases':
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2phG.UpdatePhaseData(self,item,data,oldPage)
    elif self.PatternTree.GetItemText(item) == 'Comments':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        UpdateComments(self,data)
    elif self.PatternTree.GetItemText(item) == 'Image Controls':
        self.dataFrame.SetTitle('Image Controls')
        self.PickId = item
        self.Image = self.PatternTree.GetItemParent(item)
        masks = self.PatternTree.GetItemPyData(
            GetPatternTreeItemId(self,self.Image, 'Masks'))
        data = self.PatternTree.GetItemPyData(item)
        G2imG.UpdateImageControls(self,data,masks)
        G2plt.PlotImage(self)
    elif self.PatternTree.GetItemText(item) == 'Masks':
        self.dataFrame.SetTitle('Masks')
        self.PickId = item
        self.Image = self.PatternTree.GetItemParent(item)
        data = self.PatternTree.GetItemPyData(item)
        G2imG.UpdateMasks(self,data)
        G2plt.PlotImage(self)
    elif self.PatternTree.GetItemText(item) == 'HKL Plot Controls':
        self.PickId = item
        self.Sngl = self.PatternTree.GetItemParent(item)
        data = self.PatternTree.GetItemPyData(item)
        UpdateHKLControls(self,data)
        G2plt.PlotSngl(self)               
    elif self.PatternTree.GetItemText(item) == 'Peak List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.ExportPeakList.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePeakGrid(self,data)
        G2plt.PlotPatterns(self)
    elif self.PatternTree.GetItemText(item) == 'Background':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdateBackgroundGrid(self,data)
        G2plt.PlotPatterns(self)
    elif self.PatternTree.GetItemText(item) == 'Limits':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdateLimitsGrid(self,data)
        G2plt.PlotPatterns(self)
    elif self.PatternTree.GetItemText(item) == 'Instrument Parameters':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdateInstrumentGrid(self,data)
        G2plt.PlotPeakWidths(self)
    elif self.PatternTree.GetItemText(item) == 'Index Peak List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.ExportPeakList.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdateIndexPeaksGrid(self,data)
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
    elif self.PatternTree.GetItemText(item) == 'Unit Cells List':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        if not data:
            data.append([0,0.1,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, max zero error, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            self.PatternTree.SetItemPyData(item,data)                             
        G2pdG.UpdateUnitCellsGrid(self,data)
        self.dataFrame.RefineCell.Enable(True)
        self.dataFrame.IndexPeaks.Enable(True)
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
