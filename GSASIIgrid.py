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
import time
import cPickle
import GSASIIpath
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSTESTADD, wxID_ATONTESTINSERT,
] = [wx.NewId() for _init_coll_Atom_Items in range(8)]

[ wxID_PWDRADD, wxID_HKLFADD, wxID_DATADELETE,
] = [wx.NewId() for _init_coll_Data_Items in range(3)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD,
] = [wx.NewId() for _init_coll_DrawAtom_Items in range(10)]

[ wxID_IMCALIBRATE, wxID_IMINTEGRATE, wxID_IMCLEARCALIB, wxID_SAVEINTG, 
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL, wxID_IMSAVECONTROLS, wxID_IMLOADCONTROLS,
] = [wx.NewId() for _init_coll_IMAGE_Items in range(8)]

[ wxID_MASKCOPY,
] = [wx.NewId() for _init_coll_MASK_Items in range(1)]

[ wxID_PAWLEYLOAD, wxID_PAWLEYIMPORT,
] = [wx.NewId() for _init_coll_PAWLEY_Items in range(2)]

[ wxID_INSTPRMRESET,
] = [wx.NewId() for _init_coll_INST_Items in range(1)]

[ wxID_INDXRELOAD,
] = [wx.NewId() for _init_coll_IndPeaks_Items in range(1)]

[ wxID_UNDO,wxID_PEAKFIT,wxID_AUTOPEAKFIT,
] = [wx.NewId() for _init_coll_PEAK_Items in range(3)]

[  wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL, wxID_MAKENEWPHASE,
] = [wx.NewId() for _init_coll_INDEX_Items in range(4)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

class DataFrame(wx.Frame):
    def _init_coll_BlankMenu(self,parent):
        parent.Append(menu=self.Blank,title='')
        
    def _init_coll_AtomsMenu(self,parent):
        parent.Append(menu=self.AtomEdit, title='Edit')
        
    def _init_coll_DataMenu(self,parent):
        parent.Append(menu=self.DataEdit, title='Edit')
        
    def _init_coll_DrawAtomsMenu(self,parent):
        parent.Append(menu=self.DrawAtomEdit, title='Edit')

    def _init_coll_PawleyMenu(self,parent):
        parent.Append(menu=self.PawleyEdit,title='Pawley Reflections Operations')
      
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
        parent.Append(id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append atom',
            help='Inserted as an H atom')
        parent.Append(id=wxID_ATOMSTESTADD, kind=wx.ITEM_NORMAL,text='Append test point',
            help='Inserted as an H atom')
        parent.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert atom',
            help='Select atom row to insert before; inserted as an H atom')
        parent.Append(id=wxID_ATONTESTINSERT, kind=wx.ITEM_NORMAL,text='Insert test point',
            help='Select atom row to insert before; inserted as an H atom')
        parent.Append(id=wxID_ATOMSEDITDELETE, kind=wx.ITEM_NORMAL,text='Delete atom',
            help='Select atoms to delete first')
        parent.Append(id=wxID_ATOMSREFINE, kind=wx.ITEM_NORMAL,text='Set atom refinement flags',
            help='Select atoms to refine first')
        parent.Append(id=wxID_ATOMSMODIFY, kind=wx.ITEM_NORMAL,text='Modify atom parameters',
            help='Select atoms to modify first')
        parent.Append(id=wxID_ATOMSTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Select atoms to transform first')
            
    def _init_coll_Data_Items(self,parent):
        parent.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        parent.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        parent.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Delete histograms',
            help='Delete histograms from use for this phase')
            
    def _init_coll_DrawAtom_Items(self,parent):
        parent.Append(id=wxID_DRAWATOMSTYLE, kind=wx.ITEM_NORMAL,text='Atom style',
            help='Select atoms first')
        parent.Append(id=wxID_DRAWATOMLABEL, kind=wx.ITEM_NORMAL,text='Atom label',
            help='Select atoms first')
        parent.Append(id=wxID_DRAWATOMCOLOR, kind=wx.ITEM_NORMAL,text='Atom color',
            help='Select atoms first')
        parent.Append(id=wxID_DRAWATOMRESETCOLOR, kind=wx.ITEM_NORMAL,text='Reset atom colors',
            help='Resets all atom colors to defaults')
        parent.Append(id=wxID_DRAWVIEWPOINT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st atom selected')
        parent.Append(id=wxID_DRAWADDEQUIV, kind=wx.ITEM_NORMAL,text='Add atoms',
            help='Add symmetry & cell equivalents to drawing set from selected atoms')
        parent.Append(id=wxID_DRAWTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Transform selected atoms by symmetry & cell translations')
        parent.Append(id=wxID_DRAWFILLCOORD, kind=wx.ITEM_NORMAL,text='Fill CN-sphere',
            help='Fill coordination sphere for selected atoms')            
        parent.Append(id=wxID_DRAWFILLCELL, kind=wx.ITEM_NORMAL,text='Fill unit cell',
            help='Fill unit cell with selected atoms')
        parent.Append(id=wxID_DRAWDELETE, kind=wx.ITEM_NORMAL,text='Delete atoms',
            help='Delete atoms from drawing set')

    def _init_coll_Pawley_Items(self,parent):
        parent.Append(id=wxID_PAWLEYLOAD, kind=wx.ITEM_NORMAL,text='Pawley create',
            help='Initialize Pawley reflection list')
        parent.Append(id=wxID_PAWLEYIMPORT, kind=wx.ITEM_NORMAL,text='Pawley import',
            help='Import Pawley reflection list')

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
        parent.Append(help='Save image controls to file', 
            id=wxID_IMSAVECONTROLS, kind=wx.ITEM_NORMAL,text='Save Controls')
        parent.Append(help='Load image controls from file', 
            id=wxID_IMLOADCONTROLS, kind=wx.ITEM_NORMAL,text='Load Controls')

                    
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
        self.MakeNewPhase = parent.Append( id=wxID_MAKENEWPHASE, kind=wx.ITEM_NORMAL,
            text='Make new phase',help='Make new phase from selected unit cell')

    def _init_utils(self):
        self.BlankMenu = wx.MenuBar()
        
        self.AtomsMenu = wx.MenuBar()
        self.DataMenu = wx.MenuBar()
        self.DrawAtomsMenu = wx.MenuBar()
        self.PawleyMenu = wx.MenuBar()
        self.ImageMenu = wx.MenuBar()
        self.MaskMenu = wx.MenuBar()
        self.InstMenu = wx.MenuBar()
        self.PeakMenu = wx.MenuBar()
        self.IndPeaksMenu = wx.MenuBar()
        self.IndexMenu = wx.MenuBar()
        self.AtomEdit = wx.Menu(title='')
        self.DataEdit = wx.Menu(title='')
        self.DrawAtomEdit = wx.Menu(title='')
        self.PawleyEdit = wx.Menu(title='')
        self.ImageEdit = wx.Menu(title='')
        self.MaskEdit = wx.Menu(title='')
        self.InstEdit = wx.Menu(title='')
        self.PeakEdit = wx.Menu(title='')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndexEdit = wx.Menu(title='')
        self._init_coll_AtomsMenu(self.AtomsMenu)
        self._init_coll_Atom_Items(self.AtomEdit)
        self._init_coll_DataMenu(self.DataMenu)
        self._init_coll_Data_Items(self.DataEdit)
        self._init_coll_DrawAtomsMenu(self.DrawAtomsMenu)
        self._init_coll_DrawAtom_Items(self.DrawAtomEdit)
        self._init_coll_PawleyMenu(self.PawleyMenu)
        self._init_coll_Pawley_Items(self.PawleyEdit)
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
        self.MakeNewPhase.Enable(False)
        
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
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]+250))
        self.dirname = ''
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        clientSize = wx.ClientDisplayRect()
        Width[1] = min(Width[1],clientSize[2]-300)
        self.SetSize(Width)
        self.SetPosition(wx.Point(clientSize[2]-Width[0],clientSize[1]+250))
        
    def Clear(self):
        self.ClearBackground()
        self.DestroyChildren()
                   
class GSNoteBook(wx.Notebook):
    def __init__(self, parent, name='',size = None):
        wx.Notebook.__init__(self, parent, -1, name=name, style= wx.BK_TOP)
        if size: self.SetSize(size)
                                                      
    def Clear(self):        
        GSNoteBook.DeleteAllPages(self)
        
    def FindPage(self,name):
        numPage = self.GetPageCount()
        for page in range(numPage):
            if self.GetPageText(page) == name:
                return page
        
class GSGrid(wg.Grid):
    def __init__(self, parent, name=''):
        wg.Grid.__init__(self,parent,-1,name=name)                    
        self.SetSize(parent.GetClientSize())
            
    def Clear(self):
        wg.Grid.ClearGrid(self)
        
    def SetCellStyle(self,r,c,color="white",readonly=True):
        self.SetCellBackgroundColour(r,c,color)
        self.SetReadOnly(r,c,isReadOnly=readonly)
        
    def GetSelection(self):
        #this is to satisfy structure drawing stuff in G2plt when focus changes
        return None
                        
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
    '''
    #Fourier controls
    'mapType':'Fobs','d-max':100.,'d-min':0.2,'histograms':[],
    'stepSize':[0.5,0.5,0.5],'minX':[0.,0.,0.],'maxX':[1.0,1.0,1.0],
    #distance/angle controls
    'distMax':0.0,'angleMax':0.0,'useMapPeaks':False}
    '''
        
    def SetStatusLine(text):
        Status.SetStatusText(text)
                                      
    def OnNumCycles(event):
        try:
            value = max(0,min(200,int(Ncyc.GetValue())))
        except ValueError:
            value = 3
        data['Ncycles'] = value
        Ncyc.SetValue('%d'%(value))
        
    def OnConvergence(event):
        try:
            value = max(0.01,min(100.,float(Cnvrg.GetValue())))
        except ValueError:
            value = 0.01
        data['minSumShftESD'] = value
        Cnvrg.SetValue('%.2f'%(value))
        
    def OnAtomShift(event):
        try:
            value = max(0.1,min(5.,float(AtShft.GetValue())))
        except ValueError:
            value = 2.0
        data['maxShift'] = value
        AtShft.SetValue('%.1f'%(value))
        
    def OnMarquardt(event):
        try:
            value = max(1.0,min(10.0,float(Marq.GetValue())))
        except ValueError:
            value = 1.0
        data['Marquardt'] = value
        Marq.SetValue('%.2f'%(value))
        
    def OnBandWidth(event):
        try:
            value = max(0,min(200,int(Band.GetValue())))
        except ValueError:
            value = 0
        data['bandWidth'] = value
        Band.SetValue('%d'%(value))
        
    def OnRestraint(event):
        data['restraintWeight'] = Restraint.GetValue()

    if self.dataDisplay:
        self.dataDisplay.Destroy()
    if not self.dataFrame.GetStatusBar():
        Status = self.dataFrame.CreateStatusBar()
    SetStatusLine('')
    self.dataFrame.SetLabel('Controls')
    self.dataDisplay = wx.Panel(self.dataFrame)
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Refinement Controls:'),0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Max cycles: '),0,wx.ALIGN_CENTER_VERTICAL)
    Ncyc = wx.TextCtrl(self.dataDisplay,-1,value='%d'%(data['Ncycles']),style=wx.TE_PROCESS_ENTER)
    Ncyc.Bind(wx.EVT_TEXT_ENTER,OnNumCycles)
    Ncyc.Bind(wx.EVT_KILL_FOCUS,OnNumCycles)
    LSSizer.Add(Ncyc,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Min sum(shift/esd)^2: '),0,wx.ALIGN_CENTER_VERTICAL)
    Cnvrg = wx.TextCtrl(self.dataDisplay,-1,value='%.2f'%(data['minSumShftESD']),style=wx.TE_PROCESS_ENTER)
    Cnvrg.Bind(wx.EVT_TEXT_ENTER,OnConvergence)
    Cnvrg.Bind(wx.EVT_KILL_FOCUS,OnConvergence)
    LSSizer.Add(Cnvrg,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Max atom shift: '),0,wx.ALIGN_CENTER_VERTICAL)
    AtShft = wx.TextCtrl(self.dataDisplay,-1,value='%.1f'%(data['maxShift']),style=wx.TE_PROCESS_ENTER)
    AtShft.Bind(wx.EVT_TEXT_ENTER,OnAtomShift)
    AtShft.Bind(wx.EVT_KILL_FOCUS,OnAtomShift)
    LSSizer.Add(AtShft,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Marquardt factor: '),0,wx.ALIGN_CENTER_VERTICAL)
    Marq = wx.TextCtrl(self.dataDisplay,-1,value='%.2f'%(data['Marquardt']),style=wx.TE_PROCESS_ENTER)
    Marq.Bind(wx.EVT_TEXT_ENTER,OnMarquardt)
    Marq.Bind(wx.EVT_KILL_FOCUS,OnMarquardt)
    LSSizer.Add(Marq,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Matrix band width: '),0,wx.ALIGN_CENTER_VERTICAL)
    Band = wx.TextCtrl(self.dataDisplay,-1,value='%d'%(data['bandWidth']),style=wx.TE_PROCESS_ENTER)
    Band.Bind(wx.EVT_TEXT_ENTER,OnBandWidth)
    Band.Bind(wx.EVT_KILL_FOCUS,OnBandWidth)
    LSSizer.Add(Band,0,wx.ALIGN_CENTER_VERTICAL)
    Restraint = wx.CheckBox(self.dataDisplay,-1,label='Modify restraint weights?')
    Restraint.Bind(wx.EVT_CHECKBOX, OnRestraint)
    Restraint.SetValue(data['restraintWeight'])
    LSSizer.Add(Restraint,0,wx.ALIGN_CENTER_VERTICAL)
    
    
    mainSizer.Add(LSSizer)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Density Map Controls:'),0,wx.ALIGN_CENTER_VERTICAL)

    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(self.dataDisplay,label=' Distance/angle Controls:'),0,wx.ALIGN_CENTER_VERTICAL)
        
    mainSizer.Layout()    
    self.dataDisplay.SetSizer(mainSizer)
    self.dataDisplay.SetSize(mainSizer.Fit(self.dataFrame))
    self.dataFrame.setSizePosLeft(mainSizer.Fit(self.dataFrame))
     
def UpdateComments(self,data):                   
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
            if data == [0] or data == {}:           #fill in defaults
                data = {
                    #least squares controls
                    'Ncycles':3,'maxShift':2.0,'bandWidth':0,'Marquardt':1.0,'restraintWeight':False,
                    'minSumShftESD':0.01,
                    #Fourier controls
                    'mapType':'Fobs','d-max':100.,'d-min':0.2,'histograms':[],
                    'stepSize':[0.5,0.5,0.5],'minX':[0.,0.,0.],'maxX':[1.0,1.0,1.0],
                    #distance/angle controls
                    'distMax':0.0,'angleMax':0.0,'useMapPeaks':False}
                self.PatternTree.SetItemPyData(item,data)                             
            self.Refine.Enable(True)
            UpdateControls(self,data)
        elif 'IMG' in self.PatternTree.GetItemText(item):
            self.Image = item
            G2plt.PlotImage(self,newPlot=True)
        elif 'PKS' in self.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(self)
        elif 'PWDR' in self.PatternTree.GetItemText(item):            
            self.ExportPattern.Enable(True)
            G2plt.PlotPatterns(self,newPlot=True)
        elif 'HKLF' in self.PatternTree.GetItemText(item):
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
    elif self.PatternTree.GetItemText(item) == 'Sample Parameters':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)

        if 'Temperature' not in data:           #temp fix for old gpx files
            data = {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],'DisplaceX':[0.0,False],
                'DisplaceY':[0.0,False],'Diffuse':[],'Temperature':300.,'Pressure':1.0,'Humidity':0.0,'Voltage':0.0,'Force':0.0}
            self.PatternTree.SetItemPyData(item,data)
    
        G2pdG.UpdateSampleGrid(self,data)
        G2plt.PlotPatterns(self)
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
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
