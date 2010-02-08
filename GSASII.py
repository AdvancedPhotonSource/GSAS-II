#GSASII

import os.path as ospath
import sys
import math
import cPickle
import time
import numpy as np
import wx
import matplotlib as mpl
# use the newer wxmpl when needed
(main,sub) = mpl.__version__.split('.')[0:2]
if int(main) > 0 or int(sub) > 91: 
    import wxmpl131 as wxmpl
else:
    import wxmpl as wxmpl
import pylab

# determine a binary path pased on the host OS and the python version, path is relative to 
# location of this file
if sys.platform == "win32":
    bindir = 'binwin%d.%d' % sys.version_info[0:2]
elif sys.platform == "darwin":
    bindir = 'binmac%d.%d' % sys.version_info[0:2]
else:
    bindir = 'bin'
if ospath.exists(ospath.join(sys.path[0],bindir)): sys.path.insert(0,ospath.join(sys.path[0],bindir))
# load the GSAS routines
import GSASIIIO as G2IO
import GSASIIcomp as G2cmp
import GSASIIgrid as G2gd

# print versions
print "Available python module versions for pyGSASII:"
print "python:     ",sys.version[:5]
print "wxpython:   ",wx.__version__
print "matplotlib: ",mpl.__version__
print "numpy:      ",np.__version__
print "wxmpl:      ",wxmpl.__version__

__version__ = '0.1.2'

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
 wxID_GSASIIFILESAVE, wxID_GSASIIFILESAVEAS, wxID_GSASIIPEAKFIT, 
 wxID_GSASIIINDEX, wxID_GSASIIAUTOPEAKFIT, wxID_GSASIIUNDO, wxID_GSASIIREFINECELL,
] = [wx.NewId() for _init_coll_File_Items in range(10)]

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
        parent.Append(help='', id=wxID_GSASIIFILEOPEN, kind=wx.ITEM_NORMAL,
            text='Open project')
        parent.Append(help='', id=wxID_GSASIIFILESAVE, kind=wx.ITEM_NORMAL,
            text='Save project')
        parent.Append(help='', id=wxID_GSASIIFILESAVEAS, kind=wx.ITEM_NORMAL,
            text='SaveAs')
        parent.Append(help='', id=wxID_GSASIIFILECLOSE, kind=wx.ITEM_NORMAL,
            text='Close project')
        parent.Append(help='', id=wxID_GSASIIFILEEXIT, kind=wx.ITEM_NORMAL,
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
        self.PeakFit = parent.Append(help='', id=wxID_GSASIIPEAKFIT, kind=wx.ITEM_NORMAL,
            text='PeakFit')
        self.AutoPeakFit = parent.Append(help='', id=wxID_GSASIIAUTOPEAKFIT, kind=wx.ITEM_NORMAL,
            text='AutoPeakFit')
        self.RefineCell = parent.Append(help='', id=wxID_GSASIIREFINECELL, kind=wx.ITEM_NORMAL,
            text='RefineCell')
        self.IndexPeaks = parent.Append(help='', id=wxID_GSASIIINDEX, kind=wx.ITEM_NORMAL,
            text='IndexPeaks')
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.AutoPeakFit.Enable(False)
        self.IndexPeaks.Enable(False)
        self.RefineCell.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnUnDo, id=wxID_GSASIIUNDO)
        self.Bind(wx.EVT_MENU, self.OnPeakFit, id=wxID_GSASIIPEAKFIT)
        self.Bind(wx.EVT_MENU, self.OnAutoPeakFit, id=wxID_GSASIIAUTOPEAKFIT)
        self.Bind(wx.EVT_MENU, self.OnRefineCell, id=wxID_GSASIIREFINECELL)
        self.Bind(wx.EVT_MENU, self.OnIndexPeaks, id=wxID_GSASIIINDEX)
        
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
        
        self.dataDisplay = None
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        self.dirname = ''
        self.GSASprojectfile = ''
        self.Offset = 0.0
        self.Weight = False
        self.IparmName = ''
        self.NewPlot = True
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
        self.imageDefault = {}
        self.PDevent = []
        self.IMevent = []
        self.SCevent = []
        self.Sngl = 0

    def OnSize(self,event):
        w,h = self.GetClientSizeTuple()
        self.mainPanel.SetSize(wx.Size(w,h))
        self.PatternTree.SetSize(wx.Size(w,h))
                        
    def OnPatternTreeSelChanged(self, event):
        if self.PickId:
            if self.PatternTree.GetItemText(self.PickId) in ['Peak List','Limits','Peak Index List','Unit Cell List','Background']:
                ax = self.pdplot.gca()
                self.plotView = [ax.get_xlim(),ax.get_ylim()]
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
                            data: (x,y,w,yc,yb,yd) list from GetPowderData
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
                                    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L'] 
                                    v = (v[0],v[2],v[4])
                                    codes = [0,0,0,0]
                                else:
                                    names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','SH/L']
                                    codes = [0,0,0,0,0,0]
                                data.extend(v)
                                v1 = Iparm['INS  1PRCF1 '].split()                                                  
                                v = Iparm['INS  1PRCF11'].split()
                                data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
                                v = Iparm['INS  1PRCF12'].split()
                                if v1[0] == 3:
                                    data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3])])  #get LX, LY & S+H/L
                                else:
                                    data.extend([0.0,0.0,0.002])                                      #OK defaults if fxn #3 not 1st in iprm file
                                codes.extend([0,0,0,0,0,0])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(data),data,codes,names])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                            self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
                    finally:
                        wx.EndBusyCursor()
                self.PatternTree.Expand(Id)
                self.PatternTree.SelectItem(Id)
                self.NewPlot = True
                self.PlotPatterns()
    
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
                self.NewPlot = True
                self.PatternTree.Expand(Id)
        finally:
            dlg.Destroy()
            
    def OnImageRead(self,event):
        import copy
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose image file', '.', '', \
            'MAR345 (*.mar3450)|*.mar3450|ADSC Image (*.img)|*.img \
            |Perkin-Elmer TIF (*.tif)|*.tif|GE Image sum (*.sum)|*.sum|All files (*.*)|*.*',wx.OPEN)
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
                elif ext == '.mar3450':
                    Comments,Data,Size,Image = G2IO.GetMAR345Data(self.imagefile)
                elif ext == '.sum':
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
                        Data['refine'] = [True,False,True,True,True]
                        Data['showLines'] = False
                        Data['ring'] = []
                        Data['ellipses'] = []
                        Data['masks'] = []
                        Data['calibrant'] = ''
                        Data['IOradii'] = [10.,100.]
                        Data['LRazimuth'] = [-45,45]
                        Data['outChannels'] = 2500
                        Data['fullIntegrate'] = False
                    Data['setDefault'] = False
                    Data['range'] = [(Imin,Imax),[Imin,Imax]]
                    self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)
                    self.PatternTree.SetItemPyData(Id,[Size,Image])
                    self.PickId = Id
                    self.Image = Id
                    self.PlotImage()
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
                                YCsum = [0.0 for i in range(lenX)]
                                YBsum = [0.0 for i in range(lenX)]
                                YDsum = [0.0 for i in range(lenX)]
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
                    
                    self.PlotPatterns()
                    self.NewPlot = True
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
                    self.PickId = Id
                    self.Image = Id
                    self.PlotImage()
                    self.NewPlot = True
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
                    self.PlotPatterns()
                    self.NewPlot = True
            finally:
                dlg.Destroy()

    def OnFileOpenMenu(self, event):
        result = ''
        if self.PatternTree.GetChildrenCount(self.root,False):
            if self.dataFrame:
                self.dataFrame.Clear() 
            dlg = wx.MessageDialog(self, 'Overwrite?','Project exists!',  wx.OK | wx.CANCEL)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_OK:
                    self.PatternTree.DeleteChildren(self.root)
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
                    self.HKL = []
                    self.NewPlot = True
#                    self.PlotPatterns()
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
                self.NewPlot = True
                self.PlotPatterns()
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
        pylab.close('all')
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
        
    def OnPeakFit(self,event):
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
        OK,smin,Rwp,runtime,GoOn = G2cmp.DoPeakFit(peaks,background,limits,inst,data)
        G2gd.UpdatePeakGrid(self,peaks)
        self.PlotPatterns()
        if not OK:
            print 'Refinement failed'
            dlg = wx.MessageDialog(self, 'Do you want to reload now?', 'Refinement failed',  wx.YES_NO)
            try:
                if dlg.ShowModal() == wx.ID_YES:
                    self.DoUnDo()
                    self.UnDo.Enable(False)
            finally:
                dlg.Destroy()
        else:
            print "%s%7.2f%s%12.6g" % ('Rwp = ',Rwp,'%, Smin = ',smin)
            print "%s%8.3f%s " % ('fitpeak time =',runtime,'s')
            print 'finished'
        return
        
    def OnAutoPeakFit(self,event):
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
        smin = 1.0e10
        GoOn = True
        while GoOn:
            osmin = smin
            OK,smin,Rwp,runtime,GoOn = G2cmp.DoPeakFit(peaks,background,limits,inst,data)
            print GoOn
            G2gd.UpdatePeakGrid(self,peaks)
            if not OK:
                break
            self.PlotPatterns()
            print "%s%7.2f%s%12.6g" % ('Rwp = ',Rwp,'%, Smin = ',smin)
            rat = (osmin-smin)/smin
            if rat < 1.0e-4: GoOn = False
        if not OK:
            print 'Refinement failed'
            dlg = wx.MessageDialog(self, 'Do you want to reload now?', 'Refinement failed',  wx.YES_NO)
            try:
                if dlg.ShowModal() == wx.ID_YES:
                    self.DoUnDo()
                    self.UnDo.Enable(False)
            finally:
                dlg.Destroy()
        else:
            print "%s%8.3f%s " % ('fitpeak time =',runtime,'s per cycle')
            print 'finished'
        return
        
    def OnRefineCell(self,event):
        def cellPrint(ibrav,A):
            cell = G2cmp.A2cell(A)
            Vol = G2cmp.calc_V(A)
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
        A = G2cmp.cell2A(cell)
        print controls[5]
        ibrav = bravaisSymb.index(controls[5])
        dmin = G2cmp.getDmin(peaks)-0.05
        Lhkl,M20,X20 = G2cmp.refinePeaks(peaks,ibrav,A)
        controls[6:12] = G2cmp.A2cell(A)
        controls[12] = G2cmp.calc_V(A)
        data = [controls,bravais,cells,dmin]
        self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'),data)
        self.HKL = G2cmp.GenHBravais(dmin,ibrav,A)
        G2gd.UpdateUnitCellsGrid(self,data)
        print "%s%10.3f" % ('refinement M20 = ',M20)
        print 'unindexed lines = ',X20
        cellPrint(ibrav,A)
        for hkl in self.HKL:
            hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
        self.PlotPatterns()
        
    def SaveState(self):
        file = open('GSASII.save','wb')
        PatternId = self.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            cPickle.dump(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId,item)),file,1)
        file.close()
        self.UnDo.Enable(True)
         
    def OnIndexPeaks(self,event):
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
        self.IndexPeaks.Enable(False)
        self.RefineCell.Enable(False)
        OK,dmin,cells = G2cmp.DoIndexPeaks(peaks,inst,controls,bravais)
        if OK:
            data = [controls,bravais,cells,dmin]
            self.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Unit Cells List'),data)
            G2gd.UpdateUnitCellsGrid(self,data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                self.HKL = G2cmp.GenHBravais(dmin,bestCell[2],G2cmp.cell2A(bestCell[3:9]))
                for hkl in self.HKL:
                    hkl.append(2.0*asind(inst[1]/(2.*hkl[3])))             
                self.PlotPatterns()
        self.RefineCell.Enable(True)
        self.IndexPeaks.Enable(True)
        
    def ClearEventList(self,eventList):
        if eventList:
            for i in range(len(eventList)):
                self.pdplot.canvas.mpl_disconnect(eventList[i])
        return []

    def PlotPeakWidths(self):
        newPlot = False 
        PatternId = self.PatternId
        limitID = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
        if limitID:
            limits = self.PatternTree.GetItemPyData(limitID)
        else:
            return
        instParms = self.PatternTree.GetItemPyData( \
            G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
        if instParms[0][0] == 'PXC':
            lam = instParms[1][1]
            if len(instParms[1]) == 12:
                GU,GV,GW,LX,LY = instParms[0][6:11]
            else:
                GU,GV,GW,LX,LY = instParms[0][4:9]
        peakID = G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List')
        if peakID:
            peaks = self.PatternTree.GetItemPyData(peakID)
        else:
            peaks = []
        try:
            self.pdplot.clear()
            self.IMevent = self.ClearEventList(self.IMevent)
            self.PDevent = self.ClearEventList(self.PDevent)
            self.SCevent = self.ClearEventList(self.SCevent)
            self.pdplot.canvas.set_window_title('Peak Widths')
        except:
            self.pdplot.clear()
            self.pdplot = pylab.figure(facecolor='white')
            self.pdplot.canvas.set_window_title('Peak Widths')
            self.NewPlot = True
            newPlot = True
        self.pdplot.canvas.SetToolTipString('')
        colors=['b','g','r','c','m','k']
        Xmin,Xmax = limits[1]
        Xmin = min(0.5,max(Xmin,1))
        Xmin /= 2
        Xmax /= 2
        nPts = 100
        delt = (Xmax-Xmin)/nPts
        thetas = []
        for i in range(nPts):
            thetas.append(Xmin+i*delt)
        X = []
        Y = []
        Z = []
        W = []
        sig = lambda Th,U,V,W: 1.17741*math.sqrt(U*tand(Th)**2+V*tand(Th)+W)*math.pi/18000.
        gam = lambda Th,X,Y: (X/cosd(Th)+Y*tand(Th))*math.pi/18000.
        gamFW = lambda s,g: math.exp(math.log(g**5+2.69269*g**4*s+2.42843*g**3*s**2+4.47163*g**2*s**3+0.07842*g*s**4+s**5)/5.)
        for theta in thetas:
            X.append(4.0*math.pi*sind(theta)/lam)              #q
            s = sig(theta,GU,GV,GW)
            g = gam(theta,LX,LY)
            G = gamFW(g,s)
            Y.append(s/tand(theta))
            Z.append(g/tand(theta))
            W.append(G/tand(theta))
        ax = self.pdplot.add_subplot(111)
        ax.clear()
        ax.set_title('Instrument and sample peak widths')
        ax.set_ylabel(r'$\Delta q/q, \Delta d/d$',fontsize=14)
        ax.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
        ax.plot(X,Y,color='r',label='Gaussian')
        ax.plot(X,Z,color='g',label='Lorentzian')
        ax.plot(X,W,color='b',label='G+L')
        X = []
        Y = []
        Z = []
        W = []
        for peak in peaks:
            X.append(4.0*math.pi*sind(peak[0]/2.0)/lam)
            s = 1.17741*math.sqrt(peak[4])*math.pi/18000.
            g = peak[6]*math.pi/18000.
            G = gamFW(g,s)
            Y.append(s/tand(peak[0]/2.))
            Z.append(g/tand(peak[0]/2.))
            W.append(G/tand(peak[0]/2.))
        ax.plot(X,Y,'+',color='r',label='G peak')
        ax.plot(X,Z,'+',color='g',label='L peak')
        ax.plot(X,W,'+',color='b',label='G+L peak')
        ax.legend(loc='best')
        self.NewPlot = True
            
        if newPlot:
            pylab.show()
        else:                       #1st plot
            pylab.draw()
            
    def PlotImage(self):
        from matplotlib.patches import Ellipse

        def OnImMotion(event):
            self.pdplot.canvas.SetToolTipString('')
            if (xlim[0] < event.xdata < xlim[1]) & (ylim[0] > event.ydata > ylim[1]):
                item = self.itemPicked
                if item and self.PatternTree.GetItemText(self.PickId) == 'Image Controls':
                    if 'Text' in str(item) and Data['refine'][0]:
                        self.pdplot.canvas.SetToolTipString('%8.3f %8.3fmm'%(event.xdata/scalex,event.ydata/scaley))
                    else:
                        xcent,ycent = Data['center']
                        xpos = event.xdata-xcent*scalex
                        ypos = event.ydata-ycent*scaley
                        if 'line2' in  str(item) or 'line3' in str(item) and not Data['fullIntegrate']:
                            ang = int(atan2d(-ypos,xpos))
                            self.pdplot.canvas.SetToolTipString('%6d deg'%(ang))
                        elif 'line0' in  str(item) or 'line1' in str(item):
                            radius = math.sqrt(xpos**2+ypos**2)
                            self.pdplot.canvas.SetToolTipString('%8.3fmm'%(radius/scalex))                           
                else:
                    xpos = int(event.xdata)*self.imScale
                    ypos = int(event.ydata)*self.imScale
                    self.pdplot.canvas.SetToolTipString('%6d'%(self.ImageZ[ypos][xpos]))

        def OnImPlotKeyPress(event):
            if self.PatternTree.GetItemText(self.PickId) == 'Image Controls':
                Data = self.PatternTree.GetItemPyData(self.PickId)
                pixelSize = Data['pixelSize']
                size = len(self.ImageZ)
                Xpos = event.xdata
                if Xpos:
                    Xpos = int(Xpos)*self.imScale
                else:                   #got point out of frame
                    return
                Ypos = int(event.ydata)*self.imScale
                if event.key == 'c' and Data['refine'][0]:
                    cent = Data['center'] = [Xpos*pixelSize[0]/1000.,Ypos*pixelSize[1]/1000.] #convert to mm
                    self.centText.SetValue(("%8.3f,%8.3f" % (cent[0],cent[1])))
                elif event.key == 'r':
                    Xpos,Ypos,I,J = G2cmp.ImageLocalMax(self.ImageZ,20,Xpos,Ypos)
                    if I and J:
                        xpos = Xpos*pixelSize[0]/1000.
                        ypos = Ypos*pixelSize[1]/1000.
                        Data['ring'].append([xpos,ypos])
                elif event.key == 'd':
                    scale = self.imScale*pixelSize[0]/1000.
                    xypos = [event.xdata*scale,event.ydata*scale]
                    rings = Data['ring']
                    for ring in rings:
                        if np.allclose(ring,xypos,.01,0):
                            rings.remove(ring)                                               
                elif event.key == 'm':
                    xpos = Xpos*pixelSize[0]/1000.
                    ypos = Ypos*pixelSize[1]/1000.
                    print 'mask = ',xpos,ypos
                self.PlotImage()
                
        def OnImPick(event):
            if self.itemPicked is not None: return
            pick = event.artist
            self.itemPicked = pick
            
        def OnImRelease(event):
            if self.itemPicked is None or self.PatternTree.GetItemText(self.PickId) != 'Image Controls': return
            Data = self.PatternTree.GetItemPyData(self.PickId)
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                xcent,ycent = Data['center']
                xcent *= scalex
                ycent *= scaley
                ypos = event.ydata
                xpos -= xcent
                ypos -= ycent
                radius = math.sqrt(xpos**2+ypos**2)
                xpos /= radius
                ypos /= radius
                ang = int(atan2d(-ypos,xpos))
                if 'Line2D' in str(self.itemPicked):
                    if 'line2' in str(self.itemPicked) and not Data['fullIntegrate']:
                        Data['LRazimuth'][0] = ang
                    elif 'line3' in str(self.itemPicked) and not Data['fullIntegrate']:
                        Data['LRazimuth'][1] = ang
                    elif 'line0' in str(self.itemPicked):
                        Data['IOradii'][0] = radius/scalex
                    elif 'line1' in str(self.itemPicked):
                        Data['IOradii'][1] = radius/scalex
                    if Data['LRazimuth'][1] < Data['LRazimuth'][0]:
                        Data['LRazimuth'][1] += 360
                    if  Data['IOradii'][0] > Data['IOradii'][1]:
                        Data['IOradii'] = G2cmp.SwapXY(Data['IOradii'][0],Data['IOradii'][1])
                    self.IOradText.SetValue("%8.3f,%8.3f" % (Data['IOradii'][0],Data['IOradii'][1]))
                    self.LRazim.SetValue("%6d,%6d" % (Data['LRazimuth'][0],Data['LRazimuth'][1]))
                elif 'Text' in str(self.itemPicked) and Data['refine'][0]:
                    cent = Data['center'] = [event.xdata/scalex,event.ydata/scalex]
                    try:
                        self.centText.SetValue(("%8.3f,%8.3f" % (cent[0],cent[1])))
                    except AttributeError:
                        pass
                self.PlotImage()
            self.itemPicked = None
            
        newPlot = False
        self.NewPlot = True
        self.itemPicked = None 
        try:
            self.pdplot.clear()
            self.pdplot.canvas.toolbar.set_history_buttons()
            self.pdplot.canvas.set_window_title('2D Powder Image')
            self.PDevent = self.ClearEventList(self.PDevent)
            self.SCevent = self.ClearEventList(self.SCevent)
        except:
            self.pdplot = pylab.figure(facecolor='white')
            self.pdplot.clear()
            self.pdplot.canvas.set_window_title('2D Powder Image')
            self.NewPlot = True                     #to make sure subsequent 1-D plot will be OK
            newPlot = True
        if not self.IMevent:
            self.IMevent.append(self.pdplot.canvas.mpl_connect('key_press_event', OnImPlotKeyPress))
            self.IMevent.append(self.pdplot.canvas.mpl_connect('motion_notify_event', OnImMotion))
            self.IMevent.append(self.pdplot.canvas.mpl_connect('pick_event', OnImPick))
            self.IMevent.append(self.pdplot.canvas.mpl_connect('button_release_event', OnImRelease))            
        PickId = self.PickId
        ax = self.pdplot.add_subplot(111)
        self.PlotAX = ax
        ax.set_title(self.PatternTree.GetItemText(self.Image)[4:])
        size,self.ImageZ = self.PatternTree.GetItemPyData(self.Image)
        Data = self.PatternTree.GetItemPyData( \
            G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
        self.imScale = 1
        if len(self.ImageZ) > 1024:
            self.imScale = len(self.ImageZ)/1024
        xmax = len(self.ImageZ)/self.imScale
        xlim = (-0.5,xmax-.5)
        ylim = (xmax-.5,-0.5,)
        if self.Img:
            xlim = self.Img.axes.get_xlim()
            ylim = self.Img.axes.get_ylim()
        pixelSize = Data['pixelSize']
        Data['scalex'] = scalex = 1000./(pixelSize[0]*self.imScale)
        Data['scaley'] = scaley = 1000./(pixelSize[1]*self.imScale)
        Imin,Imax = Data['range'][1]
        acolor = mpl.cm.get_cmap(Data['color'])
        xcent,ycent = Data['center']
        xcent *= scalex
        ycent *= scaley
        ax.set_xlabel('Image x-axis/'+str(self.imScale),fontsize=12)
        ax.set_ylabel('Image y-axis/'+str(self.imScale),fontsize=12)
        self.Img = ax.imshow(self.ImageZ[::self.imScale,::self.imScale], \
            aspect='equal',origin='upper',cmap=acolor, \
            interpolation='nearest',vmin=Imin,vmax=Imax)
        ax.text(xcent,ycent,'+',ha='center',va='center',picker=3)
        if Data['showLines']:
            LRAzim = Data['LRazimuth']
            IOradii = Data['IOradii']
            arcxI = arcyI = np.array(range(LRAzim[0],LRAzim[1]+1))
            arcxI = np.sin(arcxI*math.pi/180.)*scalex*Data['IOradii'][0]+xcent
            arcyI = -np.cos(arcyI*math.pi/180.)*scaley*Data['IOradii'][0]+ycent
            ax.plot(arcxI,arcyI,picker=3)
            arcxO = arcyO = np.array(range(LRAzim[0],LRAzim[1]+1))
            arcxO = np.sin(arcxO*math.pi/180.)*scalex*Data['IOradii'][1]+xcent
            arcyO = -np.cos(arcyO*math.pi/180.)*scaley*Data['IOradii'][1]+ycent
            ax.plot(arcxO,arcyO,picker=3)
            if not Data['fullIntegrate']:
                xbeg = arcxI[0]
                ybeg = arcyI[0]
                ax.plot([xbeg,sind(LRAzim[0])*IOradii[1]*scalex+xcent],
                    [ybeg,-cosd(LRAzim[0])*IOradii[1]*scaley+ycent],picker=3)
                xbeg = arcxI[-1]
                ybeg = arcyI[-1]
                ax.plot([xbeg,sind(LRAzim[1])*IOradii[1]*scalex+xcent],
                    [ybeg,-cosd(LRAzim[1])*IOradii[1]*scaley+ycent],picker=3)
        for xring,yring in Data['ring']:
            xring *= scalex
            yring *= scaley
            ax.text(xring,yring,'+',ha='center',va='center',picker=3)
        for ellipse in Data['ellipses']:
            cent,phi,[width,height] = ellipse
            ax.add_artist(Ellipse([cent[0]*scalex,cent[1]*scaley],2*width*scalex,2*height*scalex,phi,fc=None))
        self.Img.axes.set_xlim(xlim)
        self.Img.axes.set_ylim(ylim)
        self.pdplot.colorbar(self.Img)
        if newPlot:
            pylab.show()
        else:                       #1st plot
            pylab.draw()        
                
    def PlotPatterns(self):
        
        def OnPick(event):
            if self.itemPicked is not None: return
            PatternId = self.PatternId
            PickId = self.PickId
            pick = event.artist
            mouse = event.mouseevent
            xpos = pick.get_xdata()
            ypos = pick.get_ydata()
            ind = event.ind
            xy = zip(xpos[ind],ypos[ind])
            if self.PatternTree.GetItemText(PickId) == 'Peak List':
                if ind.all() != [0]:                                    #picked a data point
                    inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
                    if len(inst[1]) == 10:
                        ins = inst[1][4:10]
                    else:
                        ins = inst[1][6:12]    
                    sig = ins[0]*tand(xy[0][0]/2.0)**2+ins[1]*tand(xy[0][0]/2.0)+ins[2]
                    gam = ins[3]/cosd(xy[0][0]/2.0)+ins[4]*tand(xy[0][0]/2.0)           
                    data = self.PatternTree.GetItemPyData(self.PickId)
                    XY = [xy[0][0],0, xy[0][1],1, sig,0, gam,0, ins[5],0]       #default refine intensity 1st   
                    data.append(XY)
                    G2gd.UpdatePeakGrid(self,data)
                else:                                                   #picked a peak list line
                    self.itemPicked = pick
            elif self.PatternTree.GetItemText(PickId) == 'Limits':
                if ind.all() != [0]:                                    #picked a data point
                    LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
                    data = self.PatternTree.GetItemPyData(LimitId)
                    if mouse.button==1:
                        data[1][0] = min(xy[0][0],data[1][1])
                    if mouse.button==3:
                        data[1][1] = max(xy[0][0],data[1][0])
                    self.PatternTree.SetItemPyData(LimitId,data)
                    G2gd.UpdateLimitsGrid(self,data)
                else:                                                   #picked a limit line
                    self.itemPicked = pick                
            self.PlotPatterns()
            
        def OnPlotKeyPress(event):
            if event.key == 'w':
                if self.Weight:
                    self.Weight = False
                else:
                    self.Weight = True
                self.PlotPatterns()
                print 'plot weighting:',self.Weight
            if self.PatternTree.GetChildrenCount(self.root,False) > 1:
                if event.key == 'u' and self.Offset < 100.:
                    self.Offset += 1.
                    self.PlotPatterns()
                elif event.key == 'd' and self.Offset > 0.:
                    self.Offset -= 1.
                    self.PlotPatterns()
                elif event.key == 'c':
                    print 'contouring'
                    if self.Contour:
                        self.Contour = False
                    else:
                        self.Contour = True
                    self.PlotPatterns()
    #            elif self.Contour and event.key == 'l':
    #                event.StopPropagation()
            else:
                event.Skip()
                            
        def OnMotion(event):
            if self.itemPicked:
                xpos = event.xdata
                if xpos:                                        #avoid out of frame mouse position
                    self.pdplot.canvas.SetToolTipString('%9.3f'%(xpos))
                       
        def OnRelease(event):
            if self.itemPicked is None: return
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                lines = []
                for line in self.Lines: lines.append(line.get_xdata()[0])
                lineNo = lines.index(self.itemPicked.get_xdata()[0])
                if  lineNo in [0,1]:
                    LimitId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Limits')
                    data = self.PatternTree.GetItemPyData(LimitId)
                    data[1][lineNo] = xpos
                    self.PatternTree.SetItemPyData(LimitId,data)
                    if self.PatternTree.GetItemText(self.PickId) == 'Limits':
                        G2gd.UpdateLimitsGrid(self,data)
                else:
                    PeakId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List')
                    data = self.PatternTree.GetItemPyData(PeakId)
                    data[lineNo-2][0] = xpos
                    self.PatternTree.SetItemPyData(PeakId,data)
                    G2gd.UpdatePeakGrid(self,data)
            self.PlotPatterns()
            self.itemPicked = None    
            
        try:
            if self.NewPlot:
                self.pdplot.clear()
            self.pdplot.canvas.toolbar.set_history_buttons()
            self.pdplot.canvas.set_window_title('Powder Patterns')
            self.IMevent = self.ClearEventList(self.IMevent)
            self.SCevent = self.ClearEventList(self.SCevent)
        except:
            self.pdplot = pylab.figure(facecolor='white')
            self.pdplot.clear()
            self.pdplot.canvas.set_window_title('Powder Patterns')
            self.NewPlot = True
        if not self.PDevent:
            self.PDevent.append(self.pdplot.canvas.mpl_connect('key_press_event', OnPlotKeyPress))
            self.PDevent.append(self.pdplot.canvas.mpl_connect('pick_event', OnPick))
            self.PDevent.append(self.pdplot.canvas.mpl_connect('button_release_event', OnRelease))
            self.PDevent.append(self.pdplot.canvas.mpl_connect('motion_notify_event', OnMotion))
        PickId = self.PickId
        PatternId = self.PatternId
        colors=['b','g','r','c','m','k']
        Ymax = 1.0
        PlotList = []
        Lines = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            if 'PWDR' in self.PatternTree.GetItemText(item):
                Pattern = self.PatternTree.GetItemPyData(item)
                Pattern.append(self.PatternTree.GetItemText(item))
                PlotList.append(Pattern)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
        for Pattern in PlotList:
            xye = Pattern[1]
            Ymax = max(Ymax,max(xye[1]))
        offset = self.Offset*Ymax/100.0
        ax = self.pdplot.add_subplot(111)
        ax.cla()
        if not self.NewPlot:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
        ax.clear()
        ax.set_title('Powder Patterns')
        ax.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
        ax.set_ylabel('Intensity',fontsize=12)
        if self.Contour:
            ContourZ = []
            ContourY = []
            Nseq = 0
        for Pattern in PlotList:
            ifpicked = False
            LimitId = 0
            if PickId:
                ifpicked = Pattern[2] in self.PatternTree.GetItemText(PatternId)
                LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
            xye = Pattern[1]
            N = PlotList.index(Pattern)
            X = np.array(xye[0])
            Y = np.array(xye[1])
            Y += offset*N
            if LimitId:
                limits = self.PatternTree.GetItemPyData(LimitId)
                Lines.append(ax.axvline(limits[1][0],color='g',dashes=(5,5),picker=3))    
                Lines.append(ax.axvline(limits[1][1],color='r',dashes=(5,5),picker=3))                    
            if self.Contour:
                ContourY.append(N)
                ContourZ.append(Y)
                ContourX = X
                Nseq += 1
                ax.set_ylabel('Data sequence',fontsize=12)
            else:
                if ifpicked:
                    Z = np.array(xye[3])+offset*N
                    W = np.array(xye[4])+offset*N
                    D = np.array(xye[5])+offset*N
                    if self.Weight:
                        W2 = np.sqrt(np.array(xye[2]))
                        D *= W2
                    ax.plot(X,Y,colors[N%6]+'+',picker=3)
                    ax.plot(X,Z,colors[(N+1)%6],picker=False)
                    ax.plot(X,W,colors[(N+2)%6],picker=False)
                    ax.plot(X,D,colors[(N+3)%6],picker=False)
                    ax.axhline(0.,color=wx.BLACK)
                    self.pdplot.canvas.SetToolTipString('')
                    if self.PatternTree.GetItemText(PickId) == 'Peak List':
                        tip = 'On data point: Pick peak - L or R MB.On line: MB down to move'
                        self.pdplot.canvas.SetToolTipString(tip)
                        data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
                        for item in data:
                            Lines.append(ax.axvline(item[0],color=colors[N%6],picker=2))
                    if self.PatternTree.GetItemText(PickId) == 'Limits':
                        tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                        self.pdplot.canvas.SetToolTipString(tip)
                        data = self.LimitsTable.GetData()
                else:
                    ax.plot(xye[0],Y,colors[N%6],picker=False)
        if PickId and self.PatternTree.GetItemText(PickId) in ['Index Peak List','Unit Cells List']:
            peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
            for peak in peaks:
                ax.axvline(peak[0],color='b')
            for hkl in self.HKL:
                ax.axvline(hkl[5],color='r',dashes=(5,5))
            if self.NewPlot and peaks:
                xmin = peaks[0][0]
                xmax = peaks[-1][0]
                delt = xmax-xmin
                xlim = [max(0,xmin-delt/20.),min(180.,xmax+delt/20.)]
                ax.set_xlim(xlim)                
        if self.Contour:
            ax.contourf(ContourX,ContourY,ContourZ)
        self.Lines = Lines
        if self.NewPlot:
            if self.plotView:
                ax.set_xlim(self.plotView[0])
                ax.set_ylim(self.plotView[1])
        else:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if self.Contour:
                ax.set_ylim(0,Nseq-1)
        if self.NewPlot:
            pylab.show()
            self.NewPlot = False
        else:                       #1st plot
            pylab.draw()
        
    def PlotSngl(self):
        from matplotlib.patches import Circle

        def OnSCMotion(event):
            xpos = event.xdata
            if xpos:
                xpos = round(xpos)                                        #avoid out of frame mouse position
                ypos = round(event.ydata)
                zpos = Data['Layer']
                if '100' in Data['Zone']:
                    self.pdplot.canvas.SetToolTipString('(%3d,%3d,%3d)'%(zpos,xpos,ypos))
                elif '010' in Data['Zone']:
                    self.pdplot.canvas.SetToolTipString('(%3d,%3d,%3d)'%(xpos,zpos,ypos))
                elif '001' in Data['Zone']:
                    self.pdplot.canvas.SetToolTipString('(%3d,%3d,%3d)'%(xpos,ypos,zpos))
                    
        def OnSCPick(event):
            zpos = Data['Layer']
            pos = event.artist.center
            if '100' in Data['Zone']:
                self.pdplot.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(zpos,pos[0],pos[1]))
            elif '010' in Data['Zone']:
                self.pdplot.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],zpos,pos[1]))
            elif '001' in Data['Zone']:
                self.pdplot.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],pos[1],zpos))                 
            
        def OnSCKeyPress(event):
            print event.key
                    
        try:
            if self.NewPlot:
                self.pdplot.clear()
            self.pdplot.canvas.toolbar.set_history_buttons()
            self.pdplot.canvas.set_window_title('Structure Factors')
            self.IMevent = self.ClearEventList(self.IMevent)
            self.PDevent = self.ClearEventList(self.PDevent)
        except:
            self.pdplot = pylab.figure(facecolor='white')
            self.pdplot.clear()
            self.pdplot.canvas.set_window_title('Structure Factors')
            self.NewPlot = True
        if not self.SCevent:
            self.SCevent.append(self.pdplot.canvas.mpl_connect('key_press_event', OnSCKeyPress))
            self.SCevent.append(self.pdplot.canvas.mpl_connect('pick_event', OnSCPick))
            self.SCevent.append(self.pdplot.canvas.mpl_connect('motion_notify_event', OnSCMotion))
        PickId = self.PickId
        ax = self.pdplot.add_subplot(111)
        ax.set_aspect(aspect='equal')
        HKLref = self.PatternTree.GetItemPyData(self.Sngl)
        Data = self.PatternTree.GetItemPyData( \
            G2gd.GetPatternTreeItemId(self,self.Sngl, 'HKL Plot Controls'))

        Type = Data['Type']            
        scale = Data['Scale']
        HKLmax = Data['HKLmax']
        HKLmin = Data['HKLmin']
        FosqMax = Data['FoMax']
        FoMax = math.sqrt(FosqMax)
        ifFc = Data['ifFc']
        xlabel = ['k, h=','h, k=','h, l=']
        ylabel = ['l','l','k']
        zones = ['100','010','001']
        pzone = [[1,2],[0,2],[0,1]]
        izone = zones.index(Data['Zone'])
        ax.set_title(self.PatternTree.GetItemText(self.Sngl)[5:])
        for h,k,l,Fosq,sig,Fcsq,x,x,x in HKLref:
            H = [h,k,l]
            if H[izone] == Data['Layer']:
                B = 0
                if Type == 'Fosq':
                    A = scale*Fosq/FosqMax
                    B = scale*Fcsq/FosqMax
                    C = abs(A-B)
                elif Type == 'Fo':
                    A = scale*math.sqrt(max(0,Fosq))/FoMax
                    B = scale*math.sqrt(max(0,Fcsq))/FoMax
                    C = abs(A-B)
                elif Type == '|DFsq|/sig':
                    A = abs(Fosq-Fcsq)/(scale*sig)
                elif Type == '|DFsq|>sig':
                    A = abs(Fosq-Fcsq)/(scale*sig)
                    if A < 1.0: A = 0                    
                elif Type == '|DFsq|>3sig':
                    A = abs(Fosq-Fcsq)/(scale*sig)
                    if A < 3.0: A = 0                    
                xy = (H[pzone[izone][0]],H[pzone[izone][1]])
                if A > 0.0:
                    ax.add_artist(Circle(xy,radius=A,ec='g',fc='w',picker=3))
                if B:
                    ax.add_artist(Circle(xy,radius=B,ec='b',fc='w'))
                    radius = C
                    if radius > 0:
                        if A > B:
                            ax.add_artist(Circle(xy,radius=radius,ec='r',fc='r'))
                        else:                    
                            ax.add_artist(Circle(xy,radius=radius,ec='g',fc='g'))
                    
        ax.set_xlabel(xlabel[izone]+str(Data['Layer']),fontsize=12)
        ax.set_ylabel(ylabel[izone],fontsize=12)
        ax.set_xlim((HKLmin[pzone[izone][0]],HKLmax[pzone[izone][0]]))
        ax.set_ylim((HKLmin[pzone[izone][1]],HKLmax[pzone[izone][1]]))

        if self.NewPlot:
            pylab.show()
            self.NewPlot = False
        else:                       #1st plot
            pylab.draw()
        
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
