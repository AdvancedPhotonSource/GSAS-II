#!/usr/bin/env python
# -*- coding: utf-8 -*-
#GSASII
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################

import os
import os.path as ospath
import sys
import math
#import cPickle
import time
import copy
import glob
import imp
import inspect
import numpy as np
import scipy as sp
import wx
import matplotlib as mpl
import wx.lib.inspection as wxeye
try:
    import OpenGL as ogl
except ImportError:
    print('*******************************************************')
    print('PyOpenGL is missing from your python installation')
    print('     - we will try to install it')
    print('*******************************************************')
    def install_with_easyinstall(package):
        easy_install.main(['-f',ospath.split(__file__)[0],package])
    from setuptools.command import easy_install
    install_with_easyinstall('PyOpenGl')
    print('*******************************************************')         
    print('OpenGL has been installed. Please restart GSAS-II')
    print('*******************************************************')         
    sys.exit()

# load the GSAS routines
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIspc as G2spc
import GSASIIstruct as G2str
import GSASIImapvars as G2mv
import GSASIIsolve as G2sol

#wx inspector - use as needed
wxInspector = False

# print versions
print "Available python module versions for GSASII:"
print "python:     ",sys.version[:5]
print "wxpython:   ",wx.__version__
print "matplotlib: ",mpl.__version__
print "numpy:      ",np.__version__
print "scipy:      ",sp.__version__
print "OpenGL:     ",ogl.__version__
try:
    import mkl
    print "Max threads ",mkl.get_max_threads()
except:
    print "MKL module not present"
__version__ = '0.1.5'
G2gd.__version__ = __version__
print "This is GSAS-II version:     ",__version__

# useful degree trig functions
sind = lambda x: math.sin(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
acosd = lambda x: 180.*math.acos(x)/math.pi
atan2d = lambda x,y: 180.*math.atan2(y,x)/math.pi

def create(parent):
    return GSASII(parent)

[wxID_PATTERNTREE, 
] = [wx.NewId() for _init_ctrls in range(1)]

[wxID_FILECLOSE, wxID_FILEEXIT, wxID_FILEOPEN,  wxID_FILESAVE, wxID_FILESAVEAS, 
wxID_REFINE,  wxID_MAKEPDFS, wxID_VIEWLSPARMS, wxID_SEQREFINE,
] = [wx.NewId() for _init_coll_File_Items in range(9)]

[wxID_PWDRREAD,wxID_SNGLREAD,wxID_ADDPHASE,wxID_DELETEPHASE,
 wxID_DATADELETE,wxID_READPEAKS,wxID_PWDSUM,wxID_IMGREAD,
 wxID_IMSUM, wxID_DATARENAME,
] = [wx.NewId() for _init_coll_Data_Items in range(10)]

[wxID_IMPORT, wxID_IMPORTPATTERN, wxID_IMPORTHKL
#, wxID_IMPORTPHASE,
#wxID_IMPORTCIF, wxID_IMPORTPDB,  
] = [wx.NewId() for _init_coll_Import_Items in range(3)]

[wxID_EXPORT, wxID_EXPORTPATTERN, wxID_EXPORTHKL, wxID_EXPORTPHASE,
wxID_EXPORTCIF, wxID_EXPORTPEAKLIST, wxID_EXPORTPDF,
] = [wx.NewId() for _init_coll_Export_Items in range(7)]

class GSASII(wx.Frame):
    
    def _init_coll_GSASIIMenu_Menus(self, parent):
        parent.Append(menu=self.File, title='File')
        parent.Append(menu=self.Data, title='Data')
        parent.Append(menu=self.Calculate, title='Calculate')
        parent.Append(menu=self.Import, title='Import')
        parent.Append(menu=self.Export, title='Export')
        parent.Append(menu=G2gd.MyHelp(self,helpType='Data tree'),title='&Help' )

    def _init_coll_File_Items(self, parent):
        parent.Append(help='Open a gsasii project file (*.gpx)', id=wxID_FILEOPEN,
             kind=wx.ITEM_NORMAL,text='Open project...')
        parent.Append(help='Save project to old file', id=wxID_FILESAVE, 
            kind=wx.ITEM_NORMAL,text='Save project')
        parent.Append(help='Save project to new file', id=wxID_FILESAVEAS, 
            kind=wx.ITEM_NORMAL,text='Save As...')
        parent.Append(help='Close project, saving is optional', id=wxID_FILECLOSE, 
            kind=wx.ITEM_NORMAL,text='Close project')
        parent.Append(help='Exit from gsasii', id=wxID_FILEEXIT, kind=wx.ITEM_NORMAL,
            text='Exit')
        self.Bind(wx.EVT_MENU, self.OnFileOpen, id=wxID_FILEOPEN)
        self.Bind(wx.EVT_MENU, self.OnFileSave, id=wxID_FILESAVE)
        self.Bind(wx.EVT_MENU, self.OnFileSaveas, id=wxID_FILESAVEAS)
        self.Bind(wx.EVT_MENU, self.OnFileClose, id=wxID_FILECLOSE)
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        
    def _init_coll_Data_Items(self,parent):
        parent.Append(help='', id=wxID_PWDRREAD, kind=wx.ITEM_NORMAL,
            text='Read powder data...')
        parent.Append(help='',id=wxID_IMGREAD, kind=wx.ITEM_NORMAL,
            text='Read image data...')
        parent.Append(help='',id=wxID_READPEAKS, kind=wx.ITEM_NORMAL,
            text='Read Powder Pattern Peaks...')
        parent.Append(help='', id=wxID_SNGLREAD, kind=wx.ITEM_NORMAL,
            text='Read single crystal data...')
        parent.Append(help='', id=wxID_PWDSUM, kind=wx.ITEM_NORMAL,
            text='Sum powder data')
        parent.Append(help='',id=wxID_IMSUM, kind=wx.ITEM_NORMAL,
            text='Sum image data')
        parent.Append(help='', id=wxID_ADDPHASE, kind=wx.ITEM_NORMAL,
            text='Add phase')
        parent.Append(help='', id=wxID_DELETEPHASE, kind=wx.ITEM_NORMAL,
            text='Delete phase')
        parent.Append(help='', id=wxID_DATARENAME, kind=wx.ITEM_NORMAL,
            text='Rename data') 
        parent.Append(help='', id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,
            text='Delete data')
        self.Bind(wx.EVT_MENU, self.OnPwdrRead, id=wxID_PWDRREAD)
        self.Bind(wx.EVT_MENU, self.OnPwdrSum, id=wxID_PWDSUM)
        self.Bind(wx.EVT_MENU, self.OnReadPowderPeaks, id=wxID_READPEAKS)
        self.Bind(wx.EVT_MENU, self.OnImageRead, id=wxID_IMGREAD)
        self.Bind(wx.EVT_MENU, self.OnImageSum, id=wxID_IMSUM)
        self.Bind(wx.EVT_MENU, self.OnSnglRead, id=wxID_SNGLREAD)
        self.Bind(wx.EVT_MENU, self.OnAddPhase, id=wxID_ADDPHASE)
        self.Bind(wx.EVT_MENU, self.OnDeletePhase, id=wxID_DELETEPHASE)
        self.Bind(wx.EVT_MENU, self.OnRenameData, id=wxID_DATARENAME)
        self.Bind(wx.EVT_MENU, self.OnDataDelete, id=wxID_DATADELETE)
                
    def _init_coll_Calculate_Items(self,parent):
        self.MakePDF = parent.Append(help='Make new PDFs from selected powder patterns', 
            id=wxID_MAKEPDFS, kind=wx.ITEM_NORMAL,text='Make new PDFs')
        self.Bind(wx.EVT_MENU, self.OnMakePDFs, id=wxID_MAKEPDFS)
        self.ViewLSParms = parent.Append(help='View least squares parameters', 
            id=wxID_VIEWLSPARMS, kind=wx.ITEM_NORMAL,text='View LS parms')
        self.Bind(wx.EVT_MENU, self.OnViewLSParms, id=wxID_VIEWLSPARMS)
        self.Refine = parent.Append(help='', id=wxID_REFINE, kind=wx.ITEM_NORMAL,
            text='Refine')
        self.Refine.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnRefine, id=wxID_REFINE)
        self.SeqRefine = parent.Append(help='', id=wxID_SEQREFINE, kind=wx.ITEM_NORMAL,
            text='Sequental refine')
        self.SeqRefine.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnSeqRefine, id=wxID_SEQREFINE)
        
    def _init_Import_Phase(self,parent):
        '''import all the G2importphase*.py files that are found in the
        path and configure the Import Phase menus accordingly
        '''

        path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # location of this file
        pathlist = sys.path[:]
        if path2GSAS2 not in pathlist: pathlist.append(path2GSAS2)
        filelist = []
        for path in pathlist:
            for filename in glob.iglob(os.path.join(path, "G2importphase*.py")):
                filelist.append(filename)    
                #print 'found',filename
        filelist = sorted(list(set(filelist))) # remove duplicates
        self.ImportPhaseReaderlist = []
        for filename in filelist:
            path,rootname = os.path.split(filename)
            pkg = os.path.splitext(rootname)[0]
            try:
                fp = None
                fp, fppath,desc = imp.find_module(pkg,[path,])
                pkg = imp.load_module(pkg,fp,fppath,desc)
                for clss in inspect.getmembers(pkg): # find classes defined in package
                    if clss[0].startswith('_'): continue
                    if inspect.isclass(clss[1]):
                        # check if we have the required methods
                        for m in 'Reader','ExtensionValidator','ContentsValidator':
                            if not hasattr(clss[1],m): break
                            if not callable(getattr(clss[1],m)): break
                        else:
                            reader = clss[1]() # create a phase import instance
                            self.ImportPhaseReaderlist.append(reader)
            except AttributeError:
                print 'Import_Phase: Attribute Error',filename
                pass
            except ImportError:
                print 'Import_Phase: Error importing file',filename
                pass
            finally:
                if fp: fp.close()
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Import Phase menu',
            submenu, help='Import phase data')
        self.PhaseImportMenuId = {}
        item = submenu.Append(wx.ID_ANY,
                              help='Import phase data based selected by extension',
                              kind=wx.ITEM_NORMAL,text='Import Phase by extension')
        self.Bind(wx.EVT_MENU, self.OnImportPhaseGeneric, id=item.GetId())
        for reader in self.ImportPhaseReaderlist:
            item = submenu.Append(wx.ID_ANY,
                help='Import specific format phase data',
                kind=wx.ITEM_NORMAL,text='Import Phase '+reader.formatName+'...')
            self.PhaseImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPhaseGeneric, id=item.GetId())

    def OnImportPhaseGeneric(self,event):
        # find out which format was requested
        reader = self.PhaseImportMenuId.get(event.GetId())
        if reader is None:
            #print "use all formats"
            readerlist = self.ImportPhaseReaderlist
            choices = "any file (*.*)|*.*"
            extdict = {}
            # compile a list of allowed extensions
            for rd in readerlist:
                fmt = rd.formatName
                for extn in rd.extensionlist:
                    if not extdict.get(extn): extdict[extn] = []
                    extdict[extn] += [fmt,]
            for extn in sorted(extdict.keys(),cmp=lambda x,y: cmp(x.lower(), y.lower())):
                fmt = ''
                for f in extdict[extn]:
                    if fmt != "": fmt += ', '
                    fmt += f
                choices += "|" + fmt + " file (*" + extn + ")|*" + extn
        else:
            readerlist = [reader,]
            # compile a list of allowed extensions
            choices = reader.formatName + " file ("
            w = ""
            for extn in reader.extensionlist:
                if w != "": w += ";"
                w += "*" + extn
            choices += w + ")|" + w
            if not reader.strictExtension:
                choices += "|any file (*.*)|*.*"
        # get the file
        dlg = wx.FileDialog(
            self, message="Choose phase input file",
            #defaultDir=os.getcwd(), 
            defaultFile="",
            wildcard=choices,
            style=wx.OPEN | wx.CHANGE_DIR
            )
        try:
            if dlg.ShowModal() == wx.ID_OK:
                file = dlg.GetPath()
            else: # cancel was pressed
                return
        finally:
            dlg.Destroy()
        # set what formats are compatible with this file
        primaryReaders = []
        secondaryReaders = []
        for reader in readerlist:
            flag = reader.ExtensionValidator(file)
            if flag is None: 
                secondaryReaders.append(reader)
            elif flag:
                primaryReaders.append(reader)
        if len(secondaryReaders) + len(primaryReaders) == 0:
            self.ErrorDialog('No matching format for file '+file,'No Format')
            return
        
        fp = None
        try:
            fp = open(file,'r')
            # try the file first with Readers that specify the
            # files extension and later with ones that allow it
            for rd in primaryReaders+secondaryReaders:
                if not rd.ContentsValidator(fp):
                    continue # rejected on cursory check
                #flag = rd.Reader(file,fp,self)
                try:
                    flag = rd.Reader(file,fp,self)
                except:
                    import traceback
                    print traceback.format_exc()
                    self.ErrorDialog('Error reading file '+file
                        +' with format '+ rd.formatName,'Read Error')
                    continue
                if not flag: continue
                dlg = wx.TextEntryDialog( # allow editing of phase name
                    self, 'Enter the name for the new phase',
                    'Edit phase name', rd.Phase['General']['Name'],
                    style=wx.OK)
                #dlg.SetValue("Python is the best!")
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    rd.Phase['General']['Name'] = dlg.GetValue()
                dlg.Destroy()
                PhaseName = rd.Phase['General']['Name']
                if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
                    sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
                else:
                    sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
                psub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
                self.PatternTree.SetItemPyData(psub,rd.Phase)
                self.PatternTree.Expand(self.root) # make sure phases are seen
                self.PatternTree.Expand(sub) 
                self.PatternTree.Expand(psub) 
                return # success
        except:
            import traceback
            print traceback.format_exc()
            self.ErrorDialog('Error on open of file '+file,'Open Error')
        finally:
            if fp: fp.close()

        return
        

    def _init_coll_Import_Items(self,parent):
        self.ImportPattern = parent.Append(help='',id=wxID_IMPORTPATTERN, kind=wx.ITEM_NORMAL,
            text='Import Powder Pattern...')
        self.ImportHKL = parent.Append(help='',id=wxID_IMPORTHKL, kind=wx.ITEM_NORMAL,
            text='Import HKLs...')
        self.Bind(wx.EVT_MENU, self.OnImportPattern, id=wxID_IMPORTPATTERN)
        self.Bind(wx.EVT_MENU, self.OnImportHKL, id=wxID_IMPORTHKL)

    def _init_coll_Export_Items(self,parent):
        self.ExportPattern = parent.Append(help='Select PWDR item to enable',id=wxID_EXPORTPATTERN, kind=wx.ITEM_NORMAL,
            text='Export Powder Patterns...')
        self.ExportPeakList = parent.Append(help='',id=wxID_EXPORTPEAKLIST, kind=wx.ITEM_NORMAL,
            text='Export All Peak Lists...')
        self.ExportHKL = parent.Append(help='',id=wxID_EXPORTHKL, kind=wx.ITEM_NORMAL,
            text='Export HKLs...')
        self.ExportPDF = parent.Append(help='Select PDF item to enable',id=wxID_EXPORTPDF, kind=wx.ITEM_NORMAL,
            text='Export PDF...')
        self.ExportPhase = parent.Append(help='',id=wxID_EXPORTPHASE, kind=wx.ITEM_NORMAL,
            text='Export Phase...')
        self.ExportCIF = parent.Append(help='',id=wxID_EXPORTCIF, kind=wx.ITEM_NORMAL,
            text='Export CIF...')
        self.ExportPattern.Enable(False)
        self.ExportPeakList.Enable(True)
        self.ExportHKL.Enable(False)
        self.ExportPDF.Enable(False)
        self.ExportPhase.Enable(False)
        self.ExportCIF.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPatterns, id=wxID_EXPORTPATTERN)
        self.Bind(wx.EVT_MENU, self.OnExportPeakList, id=wxID_EXPORTPEAKLIST)
        self.Bind(wx.EVT_MENU, self.OnExportHKL, id=wxID_EXPORTHKL)
        self.Bind(wx.EVT_MENU, self.OnExportPDF, id=wxID_EXPORTPDF)
        self.Bind(wx.EVT_MENU, self.OnExportPhase, id=wxID_EXPORTPHASE)
        self.Bind(wx.EVT_MENU, self.OnExportCIF, id=wxID_EXPORTCIF)
               
    def _init_utils(self):
        self.GSASIIMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.Data = wx.Menu(title='')        
        self.Calculate = wx.Menu(title='')        
        self.Import = wx.Menu(title='')        
        self.Export = wx.Menu(title='')        

        self._init_coll_GSASIIMenu_Menus(self.GSASIIMenu)
        self._init_coll_File_Items(self.File)
        self._init_coll_Data_Items(self.Data)
        self._init_coll_Calculate_Items(self.Calculate)
        self._init_Import_Phase(self.Import)
        self._init_coll_Import_Items(self.Import)
        self._init_coll_Export_Items(self.Export)
        
    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='GSASII', parent=parent,
            size=wx.Size(300, 250),style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II data tree')
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]))
        self._init_utils()
        self.SetMenuBar(self.GSASIIMenu)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.CreateStatusBar()
        self.mainPanel = wx.Panel(self,-1)
        
        self.PatternTree = wx.TreeCtrl(id=wxID_PATTERNTREE,
            parent=self.mainPanel, pos=wx.Point(0, 0),style=wx.TR_DEFAULT_STYLE )
        self.PatternTree.Bind(wx.EVT_TREE_SEL_CHANGED,
            self.OnPatternTreeSelChanged, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
            self.OnPatternTreeItemCollapsed, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
            self.OnPatternTreeItemExpanded, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_DELETE_ITEM,
            self.OnPatternTreeItemDelete, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_KEY_DOWN,
            self.OnPatternTreeKeyDown, id=wxID_PATTERNTREE)
        self.root = self.PatternTree.AddRoot('Loaded Data: ')
        
        plotFrame = wx.Frame(None,-1,'GSASII Plots',size=wx.Size(700,600), \
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(plotFrame)
        plotFrame.Show()
        
        self.dataDisplay = None
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        # various defaults
        self.oldFocus = None
        self.GSASprojectfile = ''
        self.dirname = ospath.expanduser('~')       #start in the users home directory by default; may be meaningless
        self.undofile = ''
        self.TreeItemDelete = False
        self.Offset = [0.0,0.0]
        self.delOffset = .02
        self.refOffset = -100.0
        self.refDelt = .01
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
        self.Interpolate = 'nearest'
        self.ContourColor = 'Paired'
        self.VcovColor = 'RdYlGn'
        self.Projection = 'equal area'
        self.logPlot = False
        self.qPlot = False
        self.Contour = False
        self.Legend = False
        self.SinglePlot = False
        self.plotView = 0
        self.Image = 0
        self.oldImagefile = ''
        self.ImageZ = []
        self.Integrate = 0
        self.Pwdr = False
        self.imageDefault = {}
        self.Sngl = 0
        self.ifGetRing = False
        self.setPoly = False
        arg = sys.argv
        if len(arg) > 1:
            self.GSASprojectfile = arg[1]
            self.dirname = ospath.dirname(arg[1])
            if self.dirname: os.chdir(self.dirname)
            G2IO.ProjFileOpen(self)
            self.PatternTree.Expand(self.root)
            self.Refine.Enable(True)
            self.SeqRefine.Enable(True)
            #self.Solve.Enable(True) # what is this? (BHT)

    def OnSize(self,event):
        w,h = self.GetClientSizeTuple()
        self.mainPanel.SetSize(wx.Size(w,h))
        self.PatternTree.SetSize(wx.Size(w,h))
                        
    def OnPatternTreeSelChanged(self, event):
        if self.TreeItemDelete:
            self.TreeItemDelete = False
        else:
            pltNum = self.G2plotNB.nb.GetSelection()
            if pltNum >= 0:                         #to avoid the startup with no plot!
                pltPage = self.G2plotNB.nb.GetPage(pltNum)
                pltPlot = pltPage.figure
            item = event.GetItem()
            G2gd.MovePatternTreeToGrid(self,item)
            if self.oldFocus:
                self.oldFocus.SetFocus()
        
    def OnPatternTreeItemCollapsed(self, event):
        event.Skip()

    def OnPatternTreeItemExpanded(self, event):
        event.Skip()
        
    def OnPatternTreeItemDelete(self, event):
        self.TreeItemDelete = True

    def OnPatternTreeItemActivated(self, event):
        event.Skip()
        
    def OnPatternTreeKeyDown(self,event):
        key = event.GetKeyCode()
        item = self.PickId
        if type(item) is int: return # is this the toplevel in tree?
        if key == wx.WXK_UP:
            self.oldFocus = wx.Window.FindFocus()
            self.PatternTree.GetPrevSibling(item)
        elif key == wx.WXK_DOWN:
            self.oldFocus = wx.Window.FindFocus()
            self.PatternTree.GetNextSibling(item)
                
    def OnPwdrRead(self, event):
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose files', '.', '', 
            'GSAS fxye files (*.fxye)|*.fxye|GSAS fxy files (*.fxy)|*.fxy|Topas xye files (*.xye)|*.xye|All files (*.*)|*.*', 
            wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filenames = dlg.GetPaths()
                filenames.sort()
                for filename in filenames:
                    Data,Iparm,Comments,Temperature = G2IO.SelectPowderData(self, filename)              #Data: list of tuples (filename,Pos,Bank)
                    if not Data:                                                    #if Data rejected by user - go to next one
                        continue
                    DataType = Iparm['INS   HTYPE ']                                #expect only 4 char string
                    DataType = DataType.strip()[0:3]                                #just 1st 3 chars
                    wx.BeginBusyCursor()
                    Sample = G2pdG.SetDefaultSample()
                    Sample['Temperature'] = Temperature
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
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',True,3,1.0,0.0,0.0],
                                {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
        
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
                                try:
                                    azm = float(Iparm['INS  1DETAZM'])
                                except KeyError:                                                #not in this Iparm file
                                    azm = 0.0
                                v = Iparm['INS  1PRCF12'].split()
                                if v1[0] == 3:
                                    data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],azm)])  #get LX, LY, S+H/L & azimuth
                                else:
                                    data.extend([0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
                                codes.extend([0,0,0,0,0,0,0])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),[tuple(data),data,codes,names])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
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
            'peak files (*.txt)|*.txt|All files (*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.HKL = []
                self.powderfile = dlg.GetPath()
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
        dlg = wx.FileDialog(self, 'Choose image files', '.', '',\
        'Any image file (*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img;*.G2img)\
        |*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img;*.G2img|\
        Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|\
        MAR file (*.mar*)|*.mar*|\
        GE Image (*.avg;*.sum)|*.avg;*.sum|\
        ADSC Image (*.img)|*.img|\
        GSAS-II Image (*.G2img)|*.G2img|\
        All files (*.*)|*.*',
        wx.OPEN | wx.MULTIPLE|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                imagefiles = dlg.GetPaths()
                imagefiles.sort()
                for imagefile in imagefiles:
                    Comments,Data,Npix,Image = G2IO.GetImageData(self,imagefile)
                    if Comments:
                        Id = self.PatternTree.AppendItem(parent=self.root,text='IMG '+ospath.basename(imagefile))
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        Imax = np.amax(Image)
                        Imin = max(0.0,np.amin(Image))          #force positive
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                            Data['showLines'] = True
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['cutoff'] = 10
                            Data['pixLimit'] = 20
                            Data['edgemin'] = 100000000
                            Data['calibdmin'] = 0.5
                            Data['calibskip'] = 0
                            Data['ellipses'] = []
                            Data['calibrant'] = ''
                        else:
                            Data['type'] = 'PWDR'
                            Data['color'] = 'binary'
                            Data['tilt'] = 0.0
                            Data['rotation'] = 0.0
                            Data['showLines'] = False
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['cutoff'] = 10
                            Data['pixLimit'] = 20
                            Data['calibdmin'] = 0.5
                            Data['calibskip'] = 0
                            Data['edgemin'] = 100000000
                            Data['ellipses'] = []
                            Data['calibrant'] = ''
                            Data['IOtth'] = [2.0,5.0]
                            Data['LRazimuth'] = [135,225]
                            Data['azmthOff'] = 0.0
                            Data['outChannels'] = 2500
                            Data['outAzimuths'] = 1
                            Data['centerAzm'] = False
                            Data['fullIntegrate'] = False
                            Data['setRings'] = False
                            Data['background image'] = ['',1.0]                            
                        Data['setDefault'] = False
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                        self.PatternTree.SetItemPyData(Id,[Npix,imagefile])
                        self.PickId = Id
                        self.Image = Id
                self.PatternTree.SelectItem(G2gd.GetPatternTreeItemId(self,Id,'Image Controls'))             #show last one
        finally:
            dlg.Destroy()
        
    def OnSnglRead(self,event):
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose file', '.', '', 
            'hkl files (*.hkl)|*.hkl|All files (*.*)|*.*', 
            wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                wx.BeginBusyCursor()
                try:
                    Data = {}
                    names = ['Type','Lam']
                    HKLref,HKLmin,HKLmax,FoMax,ifFc = G2IO.GetHKLData(filename)
                    Id = self.PatternTree.AppendItem(parent=self.root,text='HKLF '+ospath.basename(filename))
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
        if not G2gd.GetPatternTreeItemId(self,self.root,'Controls'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Controls')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Covariance'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Covariance')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Constraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Constraints')
            self.PatternTree.SetItemPyData(sub,{'Hist':[],'HAP':[],'Phase':[]})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Restraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Restraints')
            self.PatternTree.SetItemPyData(sub,{})
            
                
    class CopyDialog(wx.Dialog):
        def __init__(self,parent,title,text,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            ncols = len(data)/40+1
            dataGridSizer = wx.FlexGridSizer(rows=len(data),cols=ncols,hgap=2,vgap=2)
            for id,item in enumerate(self.data):
                ckbox = wx.CheckBox(panel,id,item[1])
                ckbox.Bind(wx.EVT_CHECKBOX,self.OnCopyChange)                    
                dataGridSizer.Add(ckbox,0,wx.LEFT,10)
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
        
        def OnCopyChange(self,event):
            id = event.GetId()
            self.data[id][0] = self.FindWindowById(id).GetValue()        
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            self.Destroy()
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            self.Destroy()
            
        def GetData(self):
            return self.data
        
    class SumDialog(wx.Dialog):
        def __init__(self,parent,title,text,dataType,data):
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
                scale = wx.TextCtrl(panel,id,'%.3f'%(item[0]),style=wx.TE_PROCESS_ENTER)
                scale.Bind(wx.EVT_TEXT_ENTER,self.OnScaleChange)
                scale.Bind(wx.EVT_KILL_FOCUS,self.OnScaleChange)
                dataGridSizer.Add(scale,0,wx.LEFT,10)
                dataGridSizer.Add(name,0,wx.RIGHT,10)
            if dataType:
                dataGridSizer.Add(wx.StaticText(panel,-1,'Sum result name: '+dataType),0, \
                    wx.LEFT|wx.TOP|wx.ALIGN_CENTER_VERTICAL,10)
                self.name = wx.TextCtrl(panel,-1,self.data[-1],size=wx.Size(200,20),style=wx.TE_PROCESS_ENTER)
                self.name.Bind(wx.EVT_TEXT_ENTER,self.OnNameChange)
                self.name.Bind(wx.EVT_KILL_FOCUS,self.OnNameChange)
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

    class ConstraintDialog(wx.Dialog):
        '''Window to edit Constraint values
        '''
        def __init__(self,parent,title,text,data):
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
                #name = wx.TextCtrl(panel,-1,item[1],size=wx.Size(200,20))
                name = wx.StaticText(panel,-1,item[1],size=wx.Size(200,20))
                #name.SetEditable(False)
                scale = wx.TextCtrl(panel,id,'%.3f'%(item[0]),style=wx.TE_PROCESS_ENTER)
                scale.Bind(wx.EVT_TEXT_ENTER,self.OnScaleChange)
                scale.Bind(wx.EVT_KILL_FOCUS,self.OnScaleChange)
                dataGridSizer.Add(scale,0,wx.LEFT,10)
                dataGridSizer.Add(name,0,wx.RIGHT,10)
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
                self.FindWindowById(id).SetValue('%.3f'%(self.data[id][0]))
            except ValueError:
                if value and '-' not in value[0]:
                    print 'bad input - numbers only'
                    self.FindWindowById(id).SetValue('0.000')
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            #self.Destroy() -- do this later, after using GetData
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            #self.Destroy()
            
        def GetData(self):
            return self.data
            
    def OnPwdrSum(self,event):
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
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each pattern in summation','PWDR',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    lenX = 0
                    Xminmax = [0,0]
                    Xsum = []
                    Ysum = []
                    Vsum = []
                    result = dlg.GetData()
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            x,y,w,yc,yb,yd = data   #numpy arrays!
                            v = 1./w
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
                                         Vsum[j] += abs(scale)*v[j]
                            else:
                                Xminmax = [x[0],x[-1]]
                                YCsum = YBsum = YDsum = [0.0 for i in range(lenX)]
                                for j,yi in enumerate(y):
                                    Xsum.append(x[j])
                                    Ysum.append(scale*yi)
                                    Vsum.append(abs(scale*v[j]))
                    Wsum = 1./np.array(Vsum)
                    outname = 'PWDR '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                                self.PatternTree.Delete(Id)
                        finally:
                            dlg2.Destroy()
                    Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        Sample = G2pdG.SetDefaultSample()
                        self.PatternTree.SetItemPyData(Id,[[''],[np.array(Xsum),np.array(Ysum),np.array(Wsum),
                            np.array(YCsum),np.array(YBsum),np.array(YDsum)]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',True,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),Inst)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
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
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each image in summation','IMG',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    imSize = 0
                    result = dlg.GetData()
                    First = True
                    Found = False
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Found = True                                
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            Npix,imagefile = data
                            image = G2IO.GetImageData(self,imagefile,imageOnly=True)
                            if First:
                                newImage = np.zeros_like(image)
                                First = False
                            if imSize:
                                if imSize != Npix:
                                    self.ErrorDialog('Image size error','Images to be summed must be same size'+ \
                                        '\nExpected:'+str(imSize)+ \
                                        '\nFound:   '+str(Npix)+'\nfor '+name)
                                    return
                                newImage = newImage+scale*image
                            else:
                                imSize = Npix
                                newImage = newImage+scale*image
                            del(image)
                    if not Found:
                        self.ErrorDialog('Image sum error','No nonzero image multipliers found')
                        return
                        
                    newImage = np.asfarray(newImage,dtype=np.float32)                        
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
                        dlg = wx.FileDialog(self, 'Choose sum image filename', '.', '', 
                            'G2img files (*.G2img)|*.G2img', 
                            wx.SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
                        if dlg.ShowModal() == wx.ID_OK:
                            newimagefile = dlg.GetPath()
                            newimagefile = G2IO.FileDlgFixExt(dlg,newimagefile)
                            G2IO.PutG2Image(newimagefile,Comments,Data,Npix,newImage)
                            Imax = np.amax(newImage)
                            Imin = np.amin(newImage)
                            newImage = []
                            self.PatternTree.SetItemPyData(Id,[imSize,newimagefile])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        del(newImage)
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['cutoff'] = 10
                        Data['pixLimit'] = 20
                        Data['ellipses'] = []
                        Data['calibrant'] = ''
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)                                            
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
                        self.PickId = G2gd.GetPatternTreeItemId(self,self.root,outname)
                        self.Image = self.PickId
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
        E,SGData = G2spc.SpcGroup('P 1')
        self.PatternTree.SetItemPyData(sub,G2IO.SetNewPhase(Name=PhaseName,SGData=SGData))
        
    def OnDeletePhase(self,event):
        #Hmm, also need to delete this phase from Reflection Lists for each PWDR histogram
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
                        name = self.PatternTree.GetItemText(item)
                        self.PatternTree.Delete(item)
                        self.G2plotNB.Delete(name)
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        name = self.PatternTree.GetItemText(item)
                        if 'PWDR' in name:
                            Id = G2gd.GetPatternTreeItemId(self,item, 'Reflection Lists')
                            refList = self.PatternTree.GetItemPyData(Id)
                            for i,item in DelList:
                                del(refList[item])
                            self.PatternTree.SetItemPyData(Id,refList)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            finally:
                dlg.Destroy()
                
    def OnRenameData(self,event):
        name = self.PatternTree.GetItemText(self.PickId)      
        if 'PWDR' in name or 'HKLF' in name or 'IMG' in name:
            dataType = name[:name.index(' ')+1]                 #includes the ' '
            dlg = wx.TextEntryDialog(self,'Data name: '+dataType,'Change data name',
                defaultValue=name[name.index(' ')+1:])
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.PatternTree.SetItemText(self.PickId,dataType+dlg.GetValue())
            finally:
                dlg.Destroy()
        
    def GetFileList(self,fileType,skip=None):        #potentially useful?
        fileList = []
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
            
    def OnDataDelete(self, event):
        TextList = ['All Data']
        DelList = []
        DelItemList = []
        ifPWDR = False
        ifIMG = False
        ifHKLF = False
        ifPDF = False
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name not in ['Notebook','Controls','Covariance','Constraints','Restraints','Phases']:
                    if 'PWDR' in name: ifPWDR = True
                    if 'IMG' in name: ifIMG = True
                    if 'HKLF' in name: ifHKLF = True
                    if 'PDF' in name: ifPDF = True
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if ifPWDR: TextList.insert(1,'All PWDR')
            if ifIMG: TextList.insert(1,'All IMG')
            if ifHKLF: TextList.insert(1,'All HKLF')
            if ifPDF: TextList.insert(1,'All PDF')                
            dlg = wx.MultiChoiceDialog(self, 'Which data to delete?', 'Delete data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append(TextList[i])
                    if 'All Data' in DelList:
                        DelList = [item for item in TextList if item[:3] != 'All']
                    elif 'All PWDR' in DelList:
                        DelList = [item for item in TextList if item[:4] == 'PWDR']
                    elif 'All IMG' in DelList:
                        DelList = [item for item in TextList if item[:3] == 'IMG']
                    elif 'All HKLF' in DelList:
                        DelList = [item for item in TextList if item[:4] == 'HKLF']
                    elif 'All PDF' in DelList:
                        DelList = [item for item in TextList if item[:3] == 'PDF']
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        if self.PatternTree.GetItemText(item) in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
                    for item in DelItemList:
                        self.PatternTree.Delete(item)
                    G2plt.PlotPatterns(self,True)                        #so plot gets updated
            finally:
                dlg.Destroy()

    def OnFileOpen(self, event):
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
                    self.GSASprojectfile = ''
#                    self.PatternTree.DeleteChildren(self.root)
                    if self.HKL: self.HKL = []
                    if self.G2plotNB.plotList:
                        self.G2plotNB.clear()
            finally:
                dlg.Destroy()
        if result != wx.ID_CANCEL:    
            if self.dataDisplay: self.dataDisplay.Destroy()
            dlg = wx.FileDialog(self, 'Choose GSAS-II project file', '.', '', 
                'GSAS-II project file (*.gpx)|*.gpx',wx.OPEN|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.GSASprojectfile = dlg.GetPath()
                    self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                    self.dirname = dlg.GetDirectory()
                    G2IO.ProjFileOpen(self)
                    self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
                    self.PatternTree.Expand(self.root)
                    self.HKL = []
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item and not Id:
                        name = self.PatternTree.GetItemText(item)
                        if name[:4] in ['PWDR','HKLF','IMG ','PDF ']:
                            Id = item
                        elif name == 'Controls':
                            data = self.PatternTree.GetItemPyData(item)
                            if data:
                                self.Refine.Enable(True)
                                self.SeqRefine.Enable(True)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                    if Id:
                        self.PatternTree.SelectItem(Id)
                    self.CheckNotebook()
            finally:
                dlg.Destroy()


    def OnFileClose(self, event):
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
                self.PatternTree.SetItemText(self.root,'Loaded Data: ')
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
        finally:
            dlg.Destroy()

    def OnFileSave(self, event):
        if self.GSASprojectfile: 
            self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
            G2IO.ProjFileSave(self)
        else:
            self.OnFileSaveas(event)

    def OnFileSaveas(self, event):
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file name', '.', '', 
            'GSAS-II project file (*.gpx)|*.gpx',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
                G2IO.ProjFileSave(self)
                self.dirname = dlg.GetDirectory()
        finally:
            dlg.Destroy()

    def ExitMain(self, event):
        if self.undofile:
            os.remove(self.undofile)
        sys.exit()
        
    def OnFileExit(self, event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def OnImportPattern(self,event):
        dlg = wx.FileDialog(self, 'Choose nonGSAS powder file', '.', '', 
            '(*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.powderfile = dlg.GetPath()
        finally:
            dlg.Destroy()
            
    def OnImportHKL(self,event):
        dlg = wx.FileDialog(self, 'Choose structure factor file', '.', '', 
            '(*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.HKLfile = dlg.GetPath()
        finally:
            dlg.Destroy()
        
    def OnExportPatterns(self,event):
        names = ['All']
        exports = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            name = self.PatternTree.GetItemText(item)
            if 'PWDR' in name:
                names.append(name)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if names:
            dlg = wx.MultiChoiceDialog(self,'Select','Powder patterns to export',names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                if sel[0] == 0:
                    exports = names[1:]
                else:
                    for x in sel:
                        exports.append(names[x])
            dlg.Destroy()
        if exports:
            dlg = wx.FileDialog(self, 'Choose output powder file name', '.', '', 
                'GSAS fxye file (*.fxye)|*.fxye|xye file (*.xye)|*.xye',
                wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    powderfile = dlg.GetPath()
                    powderfile = G2IO.FileDlgFixExt(dlg,powderfile)
                    if 'fxye' in powderfile:
                        G2IO.powderFxyeSave(self,exports,powderfile)
                    else:       #just xye
                        G2IO.powderXyeSave(self,exports,powderfile)
            finally:
                dlg.Destroy()
        
    def OnExportPeakList(self,event):
        dlg = wx.FileDialog(self, 'Choose output peak list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
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
        finally:
            dlg.Destroy()
        
    def OnExportHKL(self,event):
        event.Skip()
        
    def OnExportPDF(self,event):
        #need S(Q) and G(R) to be saved here - probably best from selection?
        names = ['All']
        exports = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            name = self.PatternTree.GetItemText(item)
            if 'PDF' in name:
                names.append(name)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if names:
            dlg = wx.MultiChoiceDialog(self,'Select','PDF patterns to export',names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                if sel[0] == 0:
                    exports = names[1:]
                else:
                    for x in sel:
                        exports.append(names[x])
            dlg.Destroy()
        if exports:
            G2IO.PDFSave(self,exports)
        
    def OnExportPhase(self,event):
        event.Skip()
        
    def OnExportCIF(self,event):
        event.Skip()

    def OnMakePDFs(self,event):
        tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
        TextList = ['All PWDR']
        PDFlist = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append(name)
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) == 1:
                self.ErrorDialog('Nothing to make PDFs for','There must be at least one "PWDR" pattern')
                return
            dlg = wx.MultiChoiceDialog(self,'Make PDF controls','Make PDF controls for:',TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: PDFlist.append(TextList[i])
                    if 0 in result:
                        PDFlist = [item for item in TextList if item[:4] == 'PWDR']                        
                    for item in PDFlist:
                        PWDRname = item[4:]
                        Id = self.PatternTree.AppendItem(parent=self.root,text='PDF '+PWDRname)
                        Data = {
                            'Sample':{'Name':item,'Mult':1.0,'Add':0.0},
                            'Sample Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},'ElList':{},
                            'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':10.0,
                            'DetType':'Image plate','ObliqCoeff':0.2,'Ruland':0.025,'QScaleLim':[0,100],
                            'Lorch':True,}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='PDF Controls'),Data)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='I(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='S(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='F(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='G(R)'+PWDRname),[])        
                self.ExportPDF.Enable(True)
            finally:
                dlg.Destroy()
                
    def GetPWDRdatafromTree(self,PWDRname):
        ''' Returns powder data from GSASII tree
        input: 
            PWDRname = powder histogram name as obtained from GetHistogramNames
        return: 
            PWDRdata = powder data dictionary with:
                Data - powder data arrays, Limits, Instrument Parameters, Sample Parameters            
        '''
        PWDRdata = {}
        PWDRdata['Data'] = self.PatternTree.GetItemPyData(PWDRname)[1]          #powder data arrays
        PWDRdata['Limits'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Limits'))
        PWDRdata['Background'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Background'))
        PWDRdata['Instrument Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Instrument Parameters'))
        PWDRdata['Sample Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Sample Parameters'))
        PWDRdata['Reflection Lists'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Reflection Lists'))
        return PWDRdata

    def GetHKLFdatafromTree(self,HKLFname):
        ''' Returns single crystal data from GSASII tree
        input: 
            HKLFname = single crystal histogram name as obtained from GetHistogramNames
        return: 
            HKLFdata = single crystal data list of reflections: for each reflection:
                HKLF = [np.array([h,k,l]),FoSq,sigFoSq,FcSq,Fcp,Fcpp,phase]
        '''
        HKLFdata = []
        while True:
            data = self.PatternTree.GetItemPyData(HKLFname)
            datum = data[0]
            if datum[0] == HKLFname:
                HKLFdata = datum[1:][0]
        return HKLFdata
                    
    def GetUsedHistogramsAndPhasesfromTree(self):
        ''' Returns all histograms that are found in any phase
        and any phase that uses a histogram
        return:
            Histograms = dictionary of histograms as {name:data,...}
            Phases = dictionary of phases that use histograms
        '''
        phaseData = {}
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            print 'no phases to be refined'
            return
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phaseData[self.PatternTree.GetItemText(item)] =  self.PatternTree.GetItemPyData(item)                
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)                
        Histograms = {}
        Phases = {}
        pId = 0
        hId = 0
        for phase in phaseData:
            Phase = phaseData[phase]
            if Phase['Histograms']:
                if phase not in Phases:
                    Phase['pId'] = pId
                    pId += 1
                    Phases[phase] = Phase
                for hist in Phase['Histograms']:
                    if hist not in Histograms:
                        item = G2gd.GetPatternTreeItemId(self,self.root,hist)
                        if 'PWDR' in hist[:4]: 
                            Histograms[hist] = self.GetPWDRdatafromTree(item)
                        elif 'HKLF' in hist[:4]:
                            Histograms[hist] = self.GetHKLFdatafromTree(item)
                        #future restraint, etc. histograms here            
                        Histograms[hist]['hId'] = hId
                        hId += 1
        return Histograms,Phases
        
    class ViewParmDialog(wx.Dialog):
        def __init__(self,parent,title,parmDict):
            wx.Dialog.__init__(self,parent,-1,title,size=(300,430),
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            panel = wx.Panel(self,size=(300,430))
            parmNames = parmDict.keys()
            parmNames.sort()
            parmText = ' p:h:Parameter       refine?              value\n'
            for name in parmNames:
                parmData = parmDict[name]
                try:
                    parmText += ' %s \t%12.4g \n'%(name.ljust(19)+'\t'+parmData[1],parmData[0])
                except TypeError:
                    pass
            parmTable = wx.TextCtrl(panel,-1,parmText,
                style=wx.TE_MULTILINE|wx.TE_READONLY,size=(290,400))
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add(parmTable)
            panel.SetSizer(mainSizer)
                            
    def OnViewLSParms(self,event):
        parmDict = {}
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        Natoms,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable = G2str.GetPhaseData(Phases,Print=False)        
        hapVary,hapDict,controlDict = G2str.GetHistogramPhaseData(Phases,Histograms,Print=False)
        histVary,histDict,controlDict = G2str.GetHistogramData(Histograms,Print=False)
        varyList = phaseVary+hapVary+histVary
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        for parm in parmDict:
            if parm.split(':')[-1] in ['Azimuth','Gonio. radius','Lam1','Lam2','Omega','Chi','Phi']:
                parmDict[parm] = [parmDict[parm],' ']
            elif parm.split(':')[-2] in ['Ax','Ay','Az','SHmodel','SHord']:
                parmDict[parm] = [parmDict[parm],' ']
            elif parm in varyList:
                parmDict[parm] = [parmDict[parm],'True']
            else:
                parmDict[parm] = [parmDict[parm],'False']
        parmDict[' Num refined'] = [len(varyList),'']
        dlg = self.ViewParmDialog(self,'Parameters for least squares',parmDict)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                print 'do something with changes?? - No!'
        finally:
            dlg.Destroy()
       
    def OnRefine(self,event):
        self.OnFileSave(event)
        #works - but it'd be better if it could restore plots
        dlg = wx.ProgressDialog('Residual','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        dlg.SetSize(Size)
        Rwp = 100.00
        try:
            Rwp = G2str.Refine(self.GSASprojectfile,dlg)
        finally:
            dlg.Destroy()        
        dlg = wx.MessageDialog(self,'Load new result?','Refinement results, Rwp =%.3f'%(Rwp),wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id = 0
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
                G2IO.ProjFileOpen(self)
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item and not Id:
                    name = self.PatternTree.GetItemText(item)
                    if name[:4] in ['PWDR','HKLF']:
                        Id = item
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                if Id:
                    self.PatternTree.SelectItem(Id)
        finally:
            dlg.Destroy()

    def OnSeqRefine(self,event):
        Id = G2gd.GetPatternTreeItemId(self,self.root,'Sequental results')
        if not Id:
            Id = self.PatternTree.AppendItem(self.root,text='Sequental results')
            self.PatternTree.SetItemPyData(Id,{})            
        self.OnFileSave(event)
        dlg = wx.ProgressDialog('Residual for histogram 0','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        dlg.SetSize(Size)
        try:
            G2str.SeqRefine(self.GSASprojectfile,dlg)
        finally:
            dlg.Destroy()        
        dlg = wx.MessageDialog(self,'Load new result?','Refinement results',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id = 0
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
                G2IO.ProjFileOpen(self)
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item and not Id:
                    name = self.PatternTree.GetItemText(item)
                    if name[:4] in ['PWDR','HKLF']:
                        Id = item
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                if Id:
                    self.PatternTree.SelectItem(Id)
        finally:
            dlg.Destroy()
        
    def ErrorDialog(self,title,message):
        dlg = wx.MessageDialog(self, message, title,  wx.OK)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()

class GSASIImain(wx.App):
    def OnInit(self):
        self.main = GSASII(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    application = GSASIImain(0)
    if wxInspector: wxeye.InspectionTool().Show()

    #application.main.OnRefine(None)
    application.MainLoop()
    
if __name__ == '__main__':
    main()
