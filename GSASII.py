#GSASII

import os.path as ospath
import sys
import math
import cPickle
import time
import copy
import numpy as np
import wx
import matplotlib as mpl

# load the GSAS routines
import GSASIIIO as G2IO
import GSASIIcomp as G2cmp
import GSASIIgrid as G2gd
import GSASIIplot as G2plt

# print versions
print "Available python module versions for pyGSASII:"
print "python:     ",sys.version[:5]
print "wxpython:   ",wx.__version__
print "matplotlib: ",mpl.__version__
print "numpy:      ",np.__version__

__version__ = '0.1.3'

# useful degree trig functions
sind = lambda x: math.sin(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
acosd = lambda x: 180.*math.acos(x)/math.pi
atan2d = lambda x,y: 180.*math.atan2(y,x)/math.pi

def create(parent):
    return GSASII(parent)

[wxID_GSASII, wxID_GSASIIPATTERNTREE, wxID_GSASIIDATA, wxID_GSASIIPICKGRID,
] = [wx.NewId() for _init_ctrls in range(4)]

[wxID_GSASIIFILECLOSE, wxID_GSASIIFILEEXIT, wxID_GSASIIFILEOPEN, 
 wxID_GSASIIFILESAVE, wxID_GSASIIFILESAVEAS, wxID_GSASIIUNDO, 
] = [wx.NewId() for _init_coll_File_Items in range(6)]

[wxID_GSASIIPWDRREAD,wxID_GSASIISNGLREAD,wxID_GSASIIADDPHASE,wxID_GSASIIDELETEPHASE,
 wxID_GSASIIDATADELETE,wxID_GSASIIREADPEAKS,wxID_GSASIIPWDSUM,wxID_GSASIIIMGREAD,
 wxID_GSASIIIMSUM,
] = [wx.NewId() for _init_coll_Data_Items in range(9)]

[wxID_GSASIIIMPORT, wxID_GSASIIIMPORTPATTERN, wxID_GSASIIIMPORTHKL, wxID_GSASIIIMPORTPHASE,
wxID_GSASIIIMPORTCIF, wxID_GSASIIIMPORTPDB,  
] = [wx.NewId() for _init_coll_Import_Items in range(6)]

[wxID_GSAIIEXPORT, wxID_GSASIIEXPORTPATTERN, wxID_GSASIIEXPORTHKL, wxID_GSASIIEXPORTPHASE,
wxID_GSASIIEXPORTCIF, wxID_GSASIIEXPORTPEAKLIST
] = [wx.NewId() for _init_coll_Export_Items in range(6)]

[wxID_GSASIIHELPABOUT, wxID_GSASIIHELPHELP, 
] = [wx.NewId() for _init_coll_Help_Items in range(2)]

class GSASII(wx.Frame):
    
    def _init_coll_GSASIIMenu_Menus(self, parent):
        parent.Append(menu=self.File, title='File')
        parent.Append(menu=self.Data, title='Data')
        parent.Append(menu=self.Calculate, title='Calculate')
        parent.Append(menu=self.Import, title='Import')
        parent.Append(menu=self.Export, title='Export')
        parent.Append(menu=self.Help, title='Help')

    def _init_coll_File_Items(self, parent):
        parent.Append(help='Open a gsasii project file (*.gpx)', id=wxID_GSASIIFILEOPEN,
             kind=wx.ITEM_NORMAL,text='Open project')
        parent.Append(help='SAve project to old file', id=wxID_GSASIIFILESAVE, 
            kind=wx.ITEM_NORMAL,text='Save project')
        parent.Append(help='Save project to new file', id=wxID_GSASIIFILESAVEAS, 
            kind=wx.ITEM_NORMAL,text='SaveAs')
        parent.Append(help='Close project, saving is optional', id=wxID_GSASIIFILECLOSE, 
            kind=wx.ITEM_NORMAL,text='Close project')
        parent.Append(help='Exit from gsasii', id=wxID_GSASIIFILEEXIT, kind=wx.ITEM_NORMAL,
            text='Exit')
        self.Bind(wx.EVT_MENU, self.OnFileOpenMenu, id=wxID_GSASIIFILEOPEN)
        self.Bind(wx.EVT_MENU, self.OnFileSaveMenu, id=wxID_GSASIIFILESAVE)
        self.Bind(wx.EVT_MENU, self.OnFileSaveasMenu, id=wxID_GSASIIFILESAVEAS)
        self.Bind(wx.EVT_MENU, self.OnFileCloseMenu, id=wxID_GSASIIFILECLOSE)
        self.Bind(wx.EVT_MENU, self.OnFileExitMenu, id=wxID_GSASIIFILEEXIT)
        
    def _init_coll_Data_Items(self,parent):
        parent.Append(help='', id=wxID_GSASIIPWDRREAD, kind=wx.ITEM_NORMAL,
            text='Read powder data')
        parent.Append(help='',id=wxID_GSASIIIMGREAD, kind=wx.ITEM_NORMAL,
            text='Read image data')
        parent.Append(help='',id=wxID_GSASIIREADPEAKS, kind=wx.ITEM_NORMAL,
            text='Read Powder Pattern Peaks')
        parent.Append(help='', id=wxID_GSASIISNGLREAD, kind=wx.ITEM_NORMAL,
            text='Read single crystal data')
        parent.Append(help='', id=wxID_GSASIIPWDSUM, kind=wx.ITEM_NORMAL,
            text='Sum powder data')
        parent.Append(help='',id=wxID_GSASIIIMSUM, kind=wx.ITEM_NORMAL,
            text='Sum image data')
        parent.Append(help='', id=wxID_GSASIIADDPHASE, kind=wx.ITEM_NORMAL,
            text='Add phase')
        parent.Append(help='', id=wxID_GSASIIDELETEPHASE, kind=wx.ITEM_NORMAL,
            text='Delete phase')
        parent.Append(help='', id=wxID_GSASIIDATADELETE, kind=wx.ITEM_NORMAL,
            text='Delete data')
        self.Bind(wx.EVT_MENU, self.OnPwdrReadMenu, id=wxID_GSASIIPWDRREAD)
        self.Bind(wx.EVT_MENU, self.OnPwdrSumMenu, id=wxID_GSASIIPWDSUM)
        self.Bind(wx.EVT_MENU, self.OnReadPowderPeaks, id=wxID_GSASIIREADPEAKS)
        self.Bind(wx.EVT_MENU, self.OnImageRead, id=wxID_GSASIIIMGREAD)
        self.Bind(wx.EVT_MENU, self.OnImageSum, id=wxID_GSASIIIMSUM)
        self.Bind(wx.EVT_MENU, self.OnSnglReadMenu, id=wxID_GSASIISNGLREAD)
        self.Bind(wx.EVT_MENU, self.OnAddPhase, id=wxID_GSASIIADDPHASE)
        self.Bind(wx.EVT_MENU, self.OnDeletePhase, id=wxID_GSASIIDELETEPHASE)
        self.Bind(wx.EVT_MENU, self.OnDataDeleteMenu, id=wxID_GSASIIDATADELETE)
                
    def _init_coll_Calculate_Items(self,parent):
        self.UnDo = parent.Append(help='', id=wxID_GSASIIUNDO, kind=wx.ITEM_NORMAL,
            text='UnDo')
        self.UnDo.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnUnDo, id=wxID_GSASIIUNDO)
        
    def _init_coll_Import_Items(self,parent):
        self.ImportPhase = parent.Append(help='Import phase data from GSAS EXP file',
            id=wxID_GSASIIIMPORTPHASE, kind=wx.ITEM_NORMAL,text='Import GSAS EXP Phase')
        self.ImportPDB = parent.Append(help='Import phase data from PDB file',
            id=wxID_GSASIIIMPORTPDB, kind=wx.ITEM_NORMAL,text='Import PDB Phase')
        self.ImportPattern = parent.Append(help='',id=wxID_GSASIIIMPORTPATTERN, kind=wx.ITEM_NORMAL,
            text='Import Powder Pattern')
        self.ImportHKL = parent.Append(help='',id=wxID_GSASIIIMPORTHKL, kind=wx.ITEM_NORMAL,
            text='Import HKLs')
        self.ImportCIF = parent.Append(help='',id=wxID_GSASIIIMPORTCIF, kind=wx.ITEM_NORMAL,
            text='Import CIF')
        self.Bind(wx.EVT_MENU, self.OnImportPhase, id=wxID_GSASIIIMPORTPHASE)
        self.Bind(wx.EVT_MENU, self.OnImportPDB, id=wxID_GSASIIIMPORTPDB)
        self.Bind(wx.EVT_MENU, self.OnImportPattern, id=wxID_GSASIIIMPORTPATTERN)
        self.Bind(wx.EVT_MENU, self.OnImportHKL, id=wxID_GSASIIIMPORTHKL)
        self.Bind(wx.EVT_MENU, self.OnImportCIF, id=wxID_GSASIIIMPORTCIF)

    def _init_coll_Export_Items(self,parent):
        self.ExportPattern = parent.Append(help='',id=wxID_GSASIIEXPORTPATTERN, kind=wx.ITEM_NORMAL,
            text='Export Powder Pattern')
        self.ExportPeakList = parent.Append(help='',id=wxID_GSASIIEXPORTPEAKLIST, kind=wx.ITEM_NORMAL,
            text='Export All Peak Lists')
        self.ExportHKL = parent.Append(help='',id=wxID_GSASIIEXPORTHKL, kind=wx.ITEM_NORMAL,
            text='Export HKLs')
        self.ExportPhase = parent.Append(help='',id=wxID_GSASIIEXPORTPHASE, kind=wx.ITEM_NORMAL,
            text='Export Phase')
        self.ExportCIF = parent.Append(help='',id=wxID_GSASIIEXPORTCIF, kind=wx.ITEM_NORMAL,
            text='Export CIF')
        self.ExportPattern.Enable(False)
        self.ExportPeakList.Enable(True)
        self.ExportHKL.Enable(False)
        self.ExportPhase.Enable(False)
        self.ExportCIF.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPattern, id=wxID_GSASIIEXPORTPATTERN)
        self.Bind(wx.EVT_MENU, self.OnExportPeakList, id=wxID_GSASIIEXPORTPEAKLIST)
        self.Bind(wx.EVT_MENU, self.OnExportHKL, id=wxID_GSASIIEXPORTHKL)
        self.Bind(wx.EVT_MENU, self.OnExportPhase, id=wxID_GSASIIEXPORTPHASE)
        self.Bind(wx.EVT_MENU, self.OnExportCIF, id=wxID_GSASIIEXPORTCIF)
               
    def _init_coll_Help_Items(self, parent):
        parent.Append(help='', id=wxID_GSASIIHELPHELP, kind=wx.ITEM_NORMAL,
            text='Help')
        parent.Append(help='', id=wxID_GSASIIHELPABOUT, kind=wx.ITEM_NORMAL,
            text='About')
        self.Bind(wx.EVT_MENU, self.OnHelpHelpMenu, id=wxID_GSASIIHELPHELP)
        self.Bind(wx.EVT_MENU, self.OnHelpAboutMenu, id=wxID_GSASIIHELPABOUT)

    def _init_utils(self):
        self.GSASIIMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.Data = wx.Menu(title='')        
        self.Calculate = wx.Menu(title='')        
        self.Import = wx.Menu(title='')        
        self.Export = wx.Menu(title='')        
        self.Help = wx.Menu(title='')

        self._init_coll_GSASIIMenu_Menus(self.GSASIIMenu)
        self._init_coll_File_Items(self.File)
        self._init_coll_Data_Items(self.Data)
        self._init_coll_Calculate_Items(self.Calculate)
        self._init_coll_Import_Items(self.Import)
        self._init_coll_Export_Items(self.Export)
        self._init_coll_Help_Items(self.Help)
        
    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, id=wxID_GSASII, name='GSASII', parent=parent,
            size=wx.Size(300, 250),style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II')
        screenSize = wx.DisplaySize()
        Size = self.GetSize()
        xPos = screenSize[0]-Size[0]
        self.SetPosition(wx.Point(xPos,0))
        self._init_utils()
        self.SetMenuBar(self.GSASIIMenu)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.CreateStatusBar()
        self.mainPanel = wx.Panel(self,-1)
        
        self.PatternTree = wx.TreeCtrl(id=wxID_GSASIIPATTERNTREE,
            parent=self.mainPanel, pos=wx.Point(0, 0),style=wx.TR_DEFAULT_STYLE )
        self.PatternTree.Bind(wx.EVT_TREE_SEL_CHANGED,
            self.OnPatternTreeSelChanged, id=wxID_GSASIIPATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
            self.OnPatternTreeItemCollapsed, id=wxID_GSASIIPATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
            self.OnPatternTreeItemExpanded, id=wxID_GSASIIPATTERNTREE)
        self.root = self.PatternTree.AddRoot("Loaded Data")
        
        plotFrame = wx.Frame(None,-1,'GSASII Plots',size=wx.Size(700,600), \
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(plotFrame)
        plotFrame.Show()
        
        self.dataDisplay = None
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        self.dirname = ''
        self.GSASprojectfile = ''
        self.Offset = 0.0
        self.Weight = False
        self.IparmName = ''
        self.IfPlot = False
        self.PatternId = 0
        self.PickId = 0
        self.PeakTable = []
        self.LimitsTable = []
        self.HKL = []
        self.Lines = []
        self.itemPicked = None
        self.dataFrame = None
        self.Contour = False
        self.plotView = 0
        self.Image = 0
        self.Img = 0
        self.TA = 0
        self.Integrate = 0
        self.Pwdr = False
        self.imageDefault = {}
        self.Sngl = 0
        self.ifGetRing = False

    def OnSize(self,event):
        w,h = self.GetClientSizeTuple()
        self.mainPanel.SetSize(wx.Size(w,h))
        self.PatternTree.SetSize(wx.Size(w,h))
                        
    def OnPatternTreeSelChanged(self, event):
        pltNum = self.G2plotNB.nb.GetSelection()
        if pltNum >= 0:                         #to avoid the startup with no plot!
            pltPage = self.G2plotNB.nb.GetPage(pltNum)
            pltPlot = pltPage.figure.gca()
#            pltPage.xylim = [pltPlot.get_xlim(),pltPlot.get_ylim()]
#            print pltPage.xylim
#            pltPage.views = copy.deepcopy(pltPage.toolbar._views)
        item = event.GetItem()
        G2gd.MovePatternTreeToGrid(self,item)
        
    def OnPatternTreeItemCollapsed(self, event):
        event.Skip()

    def OnPatternTreeItemExpanded(self, event):
        self.PatternTree.ScrollTo(self.PatternTree.GetLastChild(event.GetItem()))
        
    def OnPatternTreeDeleteItem(self, event):
        event.Skip()

    def OnPatternTreeItemActivated(self, event):
        event.Skip()
        
    def OnPwdrReadMenu(self, event):
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose files', '.', '', 
            'GSAS fxye files (*.fxye)|*.fxye|GSAS fxy files (*.fxy)|*.fxy|All files (*.*)|*.*', 
            wx.OPEN | wx.MULTIPLE)
        if self.dirname: dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filenames = dlg.GetPaths()
                filenames.sort()
                self.dirname = dlg.GetDirectory()
                for filename in filenames:
                    Data,Iparm,Comments = G2IO.SelectPowderData(self, filename)              #Data: list of tuples (filename,Pos,Bank)
                    if not Data:                                                    #if Data rejected by user - go to next one
                        continue
                    DataType = Iparm['INS   HTYPE ']                                #expect only 4 char string
                    DataType = DataType.strip()[0:3]                                #just 1st 3 chars
                    wx.BeginBusyCursor()
                    try:
                        for Item in Data:
                            vals = Item[2].split()          #split up the BANK record
                            Id = self.PatternTree.AppendItem(parent=self.root,text='PWDR '+ospath.basename(Item[0])+': '+vals[0]+vals[1])
                            data = G2IO.GetPowderData(filename,Item[1],Item[2],DataType)
                            self.PatternTree.SetItemPyData(Id,[Item,data])
                            '''
                            Each tree item data is a list with:
                            Item: the (filename,Pos,Bank) tuple
                            data: (x,y,w,yc,yb,yd) list  of np.arrays from GetPowderData
                            '''
                            
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                           
                            Tmin = min(data[0])
                            Tmax = max(data[0])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[(Tmin,Tmax),[Tmin,Tmax]])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
        
                            data = [DataType,]
                            if 'C' in DataType:
                                s = Iparm['INS  1 ICONS']
                                v = (G2IO.sfloat(s[:10]),G2IO.sfloat(s[10:20]),G2IO.sfloat(s[20:30]),G2IO.sfloat(s[55:65]),G2IO.sfloat(s[40:50])) #get lam1, lam2, zero, pola & ratio
                                if not v[1]:
                                    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
                                    v = (v[0],v[2],v[4])
                                    codes = [0,0,0,0]
                                else:
                                    names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','SH/L','Azimuth']
                                    codes = [0,0,0,0,0,0]
                                data.extend(v)
                                v1 = Iparm['INS  1PRCF1 '].split()                                                  
                                v = Iparm['INS  1PRCF11'].split()
                                data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
                                v = Iparm['INS  1PRCF12'].split()
                                if v1[0] == 3:
                                    data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],0.0)])  #get LX, LY, S+H/L & azimuth
                                else:
                                    data.extend([0.0,0.0,0.002,0.0])                                      #OK defaults if fxn #3 not 1st in iprm file
                                codes.extend([0,0,0,0,0,0,0])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(data),data,codes,names])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                            self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
                    finally:
                        wx.EndBusyCursor()
                self.PatternTree.Expand(Id)
                self.PatternTree.SelectItem(Id)
    
        finally:
            dlg.Destroy()
        
    def OnReadPowderPeaks(self,event):
        Cuka = 1.54052
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose file with peak list', '.', '', 
            'peak files (*.txt)|*.txt|All files (*.*)|*.*',wx.OPEN)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.HKL = []
                self.powderfile = dlg.GetPath()
                self.dirname = dlg.GetDirectory()
                comments,peaks = G2IO.GetPowderPeaks(self.powderfile)
                Id = self.PatternTree.AppendItem(parent=self.root,text='PKS '+ospath.basename(self.powderfile))
                data = ['PKS',Cuka,0.0]
                names = ['Type','Lam','Zero'] 
                codes = [0,0]
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(data),data,codes,names])
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),comments)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),peaks)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                self.PatternTree.Expand(Id)
                self.PatternTree.SelectItem(Id)
        finally:
            dlg.Destroy()
            
    def OnImageRead(self,event):
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose image file', '.', '',\
        'MAR345 (*.mar3450;*.mar2300)|*.mar3450;*.mar2300|ADSC Image (*.img)\
        |*.img|Perkin-Elmer TIF (*.tif)|*.tif|GE Image sum (*.sum)\
        |*.sum|GE Image avg (*.avg)|*.avg|All files (*.*)|*.*',wx.OPEN)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.imagefile = dlg.GetPath()
                self.dirname = dlg.GetDirectory()
                ext = ospath.splitext(self.imagefile)[1]
                Comments = []
                if ext == '.tif':
                    Comments,Data,Size,Image = G2IO.GetTifData(self.imagefile)
                elif ext == '.img':
                    Comments,Data,Size,Image = G2IO.GetImgData(self.imagefile)
                    Image[0][0] = 0
                elif ext == '.mar3450' or ext == '.mar2300':
                    Comments,Data,Size,Image = G2IO.GetMAR345Data(self.imagefile)
                elif ext in ['.sum','.avg']:
                    Comments,Data,Size,Image = G2IO.GetGEsumData(self.imagefile)
                if Comments:
                    Id = self.PatternTree.AppendItem(parent=self.root,text='IMG '+ospath.basename(self.imagefile))
                    self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                    Imax = np.amax(Image)
                    Imin = np.amin(Image)
                    if self.imageDefault:
                        Data = copy.copy(self.imageDefault)
                        Data['refine'] = [False,False,False,False,False]
                        Data['showLines'] = True
                    else:
                        Data['color'] = 'binary'
                        Data['tilt'] = 0.0
                        Data['rotation'] = 0.0
                        Data['showLines'] = False
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['cutoff'] = 10
                        Data['pixLimit'] = 20
                        Data['ellipses'] = []
                        Data['masks'] = []
                        Data['calibrant'] = ''
                        Data['IOtth'] = [2.0,5.0]
                        Data['LRazimuth'] = [-135,-45]
                        Data['outChannels'] = 2500
                        Data['outAzimuths'] = 1
                        Data['fullIntegrate'] = False
                        Data['setRings'] = False
                    Data['setDefault'] = False
                    Data['range'] = [(Imin,Imax),[Imin,Imax]]
                    self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)
                    self.PatternTree.SetItemPyData(Id,[Size,Image])
                    self.PickId = Id
                    self.Image = Id
                    self.PatternTree.SelectItem(Id)
                    self.PatternTree.Expand(Id)
        finally:
            dlg.Destroy()
        
    def OnSnglReadMenu(self,event):
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose file', '.', '', 
            'hkl files (*.hkl)|*.hkl|All files (*.*)|*.*', 
            wx.OPEN)
        if self.dirname: dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                self.dirname = dlg.GetDirectory()
                wx.BeginBusyCursor()
                try:
                    Data = {}
                    names = ['Type','Lam']
                    HKLref,HKLmin,HKLmax,FoMax,ifFc = G2IO.GetHKLData(filename)
                    Id = self.PatternTree.AppendItem(parent=self.root,text='SXTL '+ospath.basename(filename))
                    self.PatternTree.SetItemPyData(Id,HKLref)
                    Sub = self.PatternTree.AppendItem(Id,text='Instrument Parameters')
                    data = ['SXC',1.5428,]
                    self.PatternTree.SetItemPyData(Sub,[tuple(data),data,names])
                    Data['Type'] = 'Fosq'
                    Data['ifFc'] = ifFc
                    Data['HKLmax'] = HKLmax
                    Data['HKLmin'] = HKLmin
                    Data['FoMax'] = FoMax
                    Data['Zone'] = '001'
                    Data['Layer'] = 0
                    Data['Scale'] = 1.0
                    Data['log-lin'] = 'lin'                    
                    self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='HKL Plot Controls'),Data)
                    self.PatternTree.SelectItem(Id)
                    self.PatternTree.Expand(Id)
                    self.Sngl = Id
                finally:
                    wx.EndBusyCursor()    
        finally:
            dlg.Destroy()
            
    def CheckNotebook(self):
        if not G2gd.GetPatternTreeItemId(self,self.root,'Notebook'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Notebook')
            self.PatternTree.SetItemPyData(sub,[''])
            sub = self.PatternTree.AppendItem(parent=self.root,text='Controls')
            self.PatternTree.SetItemPyData(sub,[0])
        
        
    class SumDialog(wx.Dialog):
        def __init__(self,parent,title,text,type,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            dataGridSizer = wx.FlexGridSizer(rows=len(data),cols=2,hgap=2,vgap=2)
            for id,item in enumerate(self.data[:-1]):
                name = wx.TextCtrl(panel,-1,item[1],size=wx.Size(200,20))
                name.SetEditable(False)
                scale = wx.TextCtrl(panel,id,str(item[0]),style=wx.TE_PROCESS_ENTER)
                scale.Bind(wx.EVT_TEXT,self.OnScaleChange)                    
                dataGridSizer.Add(scale,0,wx.LEFT,10)
                dataGridSizer.Add(name,0,wx.RIGHT,10)
            dataGridSizer.Add(wx.StaticText(panel,-1,'Sum result name: '+type),0, \
                wx.LEFT|wx.TOP|wx.ALIGN_CENTER_VERTICAL,10)
            self.name = wx.TextCtrl(panel,-1,self.data[-1],size=wx.Size(200,20),style=wx.TE_PROCESS_ENTER)
            self.name.Bind(wx.EVT_TEXT,self.OnNameChange)
            dataGridSizer.Add(self.name,0,wx.RIGHT|wx.TOP,10)
            mainSizer.Add(dataGridSizer,0,wx.EXPAND)
            OkBtn = wx.Button(panel,-1,"Ok")
            OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((20,20),1)
            
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            panel.SetSizer(mainSizer)
            panel.Fit()
            self.Fit()
            
        def OnNameChange(self,event):
            self.data[-1] = self.name.GetValue() 
            
        def OnScaleChange(self,event):
            id = event.GetId()
            value = self.FindWindowById(id).GetValue()
            try:
                self.data[id][0] = float(value)
            except ValueError:
                if value and '-' not in value[0]:
                    print 'bad input - numbers only'
                    self.FindWindowById(id).SetValue('0.0')
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.SetReturnCode(wx.ID_OK)
            self.MakeModal(False)              
            self.Destroy()
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.SetReturnCode(wx.ID_CANCEL)
            self.MakeModal(False)              
            self.Destroy()
            
        def GetData(self):
                return self.data
            
    def OnPwdrSumMenu(self,event):
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = []
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetItemPyData(item)[1])    # (x,y,w,yc,yb,yd)
                    if not Inst:
                        Inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item, 'Instrument Parameters'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "PWDR" pattern')
                return
            TextList.append('default sum name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each pattern in summation','PWDR',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    lenX = 0
                    Xminmax = [0,0]
                    Xsum = []
                    Ysum = []
                    Wsum = []
                    result = dlg.GetData()
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            x,y,w,yc,yb,yd = data
                            if lenX:
                                if lenX != len(x):
                                    self.ErrorDialog('Data length error','Data to be summed must have same number of points'+ \
                                        '\nExpected:'+str(lenX)+ \
                                        '\nFound:   '+str(len(x))+'\nfor '+name)
                                    return
                            else:
                                lenX = len(x)
                            if Xminmax[1]:
                                if Xminmax != [x[0],x[-1]]:
                                    self.ErrorDialog('Data range error','Data to be summed must span same range'+ \
                                        '\nExpected:'+str(Xminmax[0])+' '+str(Xminmax[1])+ \
                                        '\nFound:   '+str(x[0])+' '+str(x[-1])+'\nfor '+name)
                                    return
                                else:
                                    for j,yi in enumerate(y):
                                         Ysum[j] += scale*yi
                            else:
                                Xminmax = [x[0],x[-1]]
                                YCsum = YBsum = YDsum = [0.0 for i in range(lenX)]
                                for j,yi in enumerate(y):
                                    Xsum.append(x[j])
                                    Ysum.append(scale*yi)
                                    Wsum.append(w[j])
                    outname = 'PWDR '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                        finally:
                            dlg2.Destroy()
                    else:
                        Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        self.PatternTree.SetItemPyData(Id,[[''],[Xsum,Ysum,Wsum,YCsum,YBsum,YDsum]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',1,3,1.0,0.0,0.0]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),Inst)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
                    
            finally:
                dlg.Destroy()

    def OnImageSum(self,event):
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = []
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetItemPyData(item))        #Size,Image
                    Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item,'Image Controls'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "IMG" pattern')
                return
            TextList.append('default sum name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each image in summation','IMG',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    imSize = 0
                    newImage = []
                    result = dlg.GetData()
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            size,image = data
                            if imSize:
                                if imSize != size:
                                    self.ErrorDialog('Image size error','Images to be summed must be same size'+ \
                                        '\nExpected:'+str(imSize)+ \
                                        '\nFound:   '+str(size)+'\nfor '+name)
                                    return
                                newImage += scale*image
                            else:
                                imSize = size
                                newImage = scale*image
                    outname = 'IMG '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                        finally:
                            dlg2.Destroy()
                    else:
                        Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        self.PatternTree.SetItemPyData(Id,[imSize,newImage])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)                                            
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
                    self.PickId = Id
                    self.Image = Id
            finally:
                dlg.Destroy()
                      
    def OnAddPhase(self,event):
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
        sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
        SGData = {'SpGrp':'P 1'}
        self.PatternTree.SetItemPyData(sub, \
            {'General':[PhaseName,'nuclear',SGData,[False,10.,10.,10.,90.,90.,90,1000.],
            [False,1.0],[],{},[],[],[]],'Atoms':[]})
        
    def OnDeletePhase(self,event):
        if self.dataFrame:
            self.dataFrame.Clear() 
        TextList = []
        DelList = []
        DelItemList = []
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            return
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                TextList.append(self.PatternTree.GetItemText(item))
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)                
            dlg = wx.MultiChoiceDialog(self, 'Which phase to delete?', 'Delete phase', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append([i,TextList[i]])
                    item, cookie = self.PatternTree.GetFirstChild(sub)
                    i = 0
                    while item:
                        if [i,self.PatternTree.GetItemText(item)] in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(sub, cookie)
                        i += 1
                    for item in DelItemList:
                        self.PatternTree.Delete(item)
            finally:
                dlg.Destroy()       
        
    def OnDataDeleteMenu(self, event):
        TextList = []
        DelList = []
        DelItemList = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if 'PWDR' in name or 'SXTL' in name or 'IMG' in name:
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
            dlg = wx.MultiChoiceDialog(self, 'Which data to delete?', 'Delete data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append(TextList[i])
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        if self.PatternTree.GetItemText(item) in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
                    for item in DelItemList:
                        self.PatternTree.Delete(item)
                    G2plt.PlotPatterns(self,True)                        #so plot gets updated
            finally:
                dlg.Destroy()

    def OnFileOpenMenu(self, event):
        result = ''
        Id = 0
        if self.PatternTree.GetChildrenCount(self.root,False):
            if self.dataFrame:
                self.dataFrame.Clear() 
            dlg = wx.MessageDialog(self, 'Overwrite?','Project exists!',  wx.OK | wx.CANCEL)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_OK:
                    self.PatternTree.DeleteChildren(self.root)
                    print 'children deleted'
            finally:
                dlg.Destroy()
        if result != wx.ID_CANCEL:    
            if self.dataDisplay: self.dataDisplay.Destroy()
            dlg = wx.FileDialog(self, 'Choose GSAS-II project file', '.', '', 
                'GSAS-II project file (*.gpx)|*.gpx',wx.OPEN)
            if self.dirname: dlg.SetDirectory(self.dirname)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.GSASprojectfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
                    G2IO.ProjFileOpen(self)
                    self.PatternTree.Expand(self.root)
                    self.HKL = []
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item and not Id:
                        name = self.PatternTree.GetItemText(item)
                        if 'PWDR' in name or 'SXTL' in name or 'IMG' in name:
                            Id = item
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                    if Id:
                        self.PatternTree.SelectItem(Id)
            finally:
                dlg.Destroy()

    def OnFileCloseMenu(self, event):
        if self.dataFrame:
            self.dataFrame.Clear()
            self.dataFrame.SetLabel('GSAS-II data display') 
        dlg = wx.MessageDialog(self, 'Save current project?', ' ', wx.YES | wx.NO | wx.CANCEL)
        try:
            result = dlg.ShowModal()
            if result == wx.ID_OK:
                self.OnFileSaveMenu(event)
            if result != wx.ID_CANCEL:
                self.GSASprojectfile = ''
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
        finally:
            dlg.Destroy()

    def OnFileSaveMenu(self, event):
        if self.GSASprojectfile: 
            G2IO.ProjFileSave(self)
        else:
            self.OnFileSaveasMenu(event)

    def OnFileSaveasMenu(self, event):
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file name', '.', '', 
            'GSAS-II project file (*.gpx)|*.gpx',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.GSASprojectfile = dlg.GetPath()
                G2IO.ProjFileSave(self)
                self.dirname = dlg.GetDirectory()
        finally:
            dlg.Destroy()

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExitMenu(self, event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def OnImportPattern(self,event):
            dlg = wx.FileDialog(self, 'Choose nonGSAS powder file', '.', '', 
                '(*.*)|*.*',wx.OPEN)
            if self.dirname:
                dlg.SetDirectory(self.dirname)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.powderfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
            finally:
                dlg.Destroy()
                
    def OnImportHKL(self,event):
            dlg = wx.FileDialog(self, 'Choose structure factor file', '.', '', 
                '(*.*)|*.*',wx.OPEN)
            if self.dirname:
                dlg.SetDirectory(self.dirname)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.HKLfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
            finally:
                dlg.Destroy()
        
    def OnImportPhase(self,event):
            dlg = wx.FileDialog(self, 'Choose GSAS EXP file', '.', '', 
                'EXP file (*.EXP)|*.EXP',wx.OPEN)
            if self.dirname:
                dlg.SetDirectory(self.dirname)
            try:
                Phase = {}
                if dlg.ShowModal() == wx.ID_OK:
                    EXPfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
                    Phase = G2IO.ReadEXPPhase(EXPfile)
            finally:
                dlg.Destroy()
            if Phase:
                PhaseName = Phase['General'][0]
                if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
                    sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
                else:
                    sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
                sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
                self.PatternTree.SetItemPyData(sub,Phase)
                
    def OnImportPDB(self,event):
            dlg = wx.FileDialog(self, 'Choose PDB file', '.', '', 
                'PDB file (*.pdb,*.ent)|*.pdb;*.ent|All files (*.*)|*.*',wx.OPEN)
            if self.dirname:
                dlg.SetDirectory(self.dirname)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    PDBfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
                    Phase = G2IO.ReadPDBPhase(PDBfile)
            finally:
                dlg.Destroy()
            if Phase:
                PhaseName = Phase['General'][0]
                if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
                    sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
                else:
                    sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
                sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
                self.PatternTree.SetItemPyData(sub,Phase)        
        
    def OnImportCIF(self,event):
            dlg = wx.FileDialog(self, 'Choose CIF file', '.', '', 
                'CIF file (*.cif)|*.cif',wx.OPEN)
            if self.dirname:
                dlg.SetDirectory(self.dirname)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.CIFfile = dlg.GetPath()
                    self.dirname = dlg.GetDirectory()
            finally:
                dlg.Destroy()
        
    def OnExportPattern(self,event):
        dlg = wx.FileDialog(self, 'Choose output powder file name', '.', '', 
            'xye file (*.xye)|*.xye',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.powderfile = dlg.GetPath()
                G2IO.PowderxyeSave(self)
                self.dirname = dlg.GetDirectory()
        finally:
            dlg.Destroy()
        
    def OnExportPeakList(self,event):
        dlg = wx.FileDialog(self, 'Choose output peak list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                file = open(self.peaklistfile,'wa')                
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item:
                    name = self.PatternTree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.PatternTree.GetFirstChild(item)
                        while item2:
                            name2 = self.PatternTree.GetItemText(item2)
                            if name2 == 'Peak List':
                                peaks = self.PatternTree.GetItemPyData(item2)
                                file.write("%s \n" % (name+' Peak List'))                
                                for peak in peaks:
                                    file.write("%10.4f %12.2f %10.3f %10.3f \n" % \
                                        (peak[0],peak[2],peak[4],peak[6]))
                            item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                file.close()
                self.dirname = dlg.GetDirectory()
        finally:
            dlg.Destroy()
        
    def OnExportHKL(self,event):
        event.Skip()
        
    def OnExportPhase(self,event):
        event.Skip()
        
    def OnExportCIF(self,event):
        event.Skip()
        
    def OnUnDo(self,event):
        self.DoUnDo()
        self.UnDo.Enable(False)
        
    def DoUnDo(self):
        print 'Undo last refinement'
        file = open('GSASII.save','rb')
        PatternId = self.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item),cPickle.load(file))
            if self.dataDisplay.GetName() == item:
                if item == 'Background':
                    G2gd.UpdateBackgroundGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
                elif item == 'Instrument Parameters':
                    G2gd.UpdateInstrumentGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
                elif item == 'Peak List':
                    G2gd.UpdatePeakGrid(self,self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, item)))
            print item,' recovered'
        file.close()
        
    def SaveState(self):
        file = open('GSASII.save','wb')
        PatternId = self.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            cPickle.dump(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId,item)),file,1)
        file.close()
        self.UnDo.Enable(True)
                
    def ErrorDialog(self,title,message):
        dlg = wx.MessageDialog(self, message, title,  wx.OK)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()

    def OnHelpHelpMenu(self, event):
        event.Skip()
        
    def OnHelpAboutMenu(self, event):
        info = wx.AboutDialogInfo()
        info.Name = 'GSAS-II'
        info.Version = '0.0.1'
        info.Copyright = '''
Robert B. Von Dreele
Argonne National Laboratory(C)
This product includes software developed
by the UChicago Argonne, LLC, as 
Operator of Argonne National Laboratory.         '''
        info.Description = '''
General Structure Analysis System - II
        '''
        wx.AboutBox(info)
        
class GSASIImain(wx.App):
    def OnInit(self):
        self.main = GSASII(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    application = GSASIImain(0)
    application.MainLoop()

if __name__ == '__main__':
    main()
