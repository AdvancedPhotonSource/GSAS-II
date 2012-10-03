# -*- coding: utf-8 -*-
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
import sys
import numpy as np
import os.path
import wx.html        # could postpone this for quicker startup
import webbrowser     # could postpone this for quicker startup
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIImath as G2mth
import GSASIIIO as G2IO
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIstruct as G2str
import GSASIImapvars as G2mv

# globals we will use later
__version__ = None # gets overridden in GSASII.py
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # save location of this file
helpLocDict = {}
htmlPanel = None
htmlFrame = None
helpMode = 'browser'
if sys.platform.lower().startswith('win'): helpMode = 'internal' # need a global control to set this
    
htmlFirstUse = True

[ wxID_FOURCALC, wxID_FOURSEARCH, wxID_FOURCLEAR, wxID_PEAKSMOVE, wxID_PEAKSCLEAR, 
    wxID_CHARGEFLIP, wxID_PEAKSUNIQUE, wxID_PEAKSDELETE, wxID_PEAKSDA,
    wxID_PEAKSDISTVP, wxID_PEAKSVIEWPT, wxID_FINDEQVPEAKS,
] = [wx.NewId() for item in range(12)]

[ wxID_PWDRADD, wxID_HKLFADD, wxID_DATADELETE,
] = [wx.NewId() for item in range(3)]

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSVIEWADD, wxID_ATOMVIEWINSERT,
    wxID_RELOADDRAWATOMS,wxID_ATOMSDISAGL,
] = [wx.NewId() for item in range(10)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD, wxID_DRAWDISAGLTOR,  wxID_DRAWPLANE,
    wxID_DRAWDISTVP,
] = [wx.NewId() for item in range(13)]

[ wxID_DRAWRESTRBOND, wxID_DRAWRESTRANGLE, wxID_DRAWRESTRPLANE, wxID_DRAWRESTRCHIRAL,
] = [wx.NewId() for item in range(4)]

[ wxID_CLEARTEXTURE,wxID_REFINETEXTURE,
] = [wx.NewId() for item in range(2)]

[ wxID_PAWLEYLOAD, wxID_PAWLEYDELETE, wxID_PAWLEYESTIMATE,
    wxID_PAWLEYUPDATE,
] = [wx.NewId() for item in range(4)]

[ wxID_IMCALIBRATE,wxID_IMRECALIBRATE,wxID_IMINTEGRATE, wxID_IMCLEARCALIB,  
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL, wxID_IMSAVECONTROLS, wxID_IMLOADCONTROLS,
] = [wx.NewId() for item in range(8)]

[ wxID_MASKCOPY, wxID_MASKSAVE, wxID_MASKLOAD,
] = [wx.NewId() for item in range(3)]

[ wxID_STRSTACOPY, wxID_STRSTAFIT, wxID_STRSTASAVE, wxID_STRSTALOAD,wxID_APPENDDZERO,
] = [wx.NewId() for item in range(5)]

[ wxID_BACKCOPY,wxID_LIMITCOPY,wxID_SAMPLECOPY, wxID_BACKFLAGCOPY, wxID_SAMPLEFLAGCOPY,
    wxID_SAMPLESAVE, wxID_SAMPLELOAD,
] = [wx.NewId() for item in range(7)]

[ wxID_INSTPRMRESET,wxID_CHANGEWAVETYPE,wxID_INSTCOPY, wxID_INSTFLAGCOPY, wxID_INSTLOAD,
    wxID_INSTSAVE,
] = [wx.NewId() for item in range(6)]

[ wxID_UNDO,wxID_LSQPEAKFIT,wxID_LSQONECYCLE,wxID_RESETSIGGAM,wxID_CLEARPEAKS,
] = [wx.NewId() for item in range(5)]

[  wxID_INDXRELOAD, wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL, wxID_MAKENEWPHASE,
] = [wx.NewId() for item in range(5)]

[ wxID_CONSTRAINTADD,wxID_EQUIVADD,wxID_HOLDADD,wxID_FUNCTADD,
] = [wx.NewId() for item in range(4)]

[ wxID_RESTRAINTADD,wxID_PWDANALYSIS, wxID_RESTSELPHASE,wxID_RESTDELETE, wxID_RESRCHANGEVAL, 
    wxID_RESTCHANGEESD,
] = [wx.NewId() for item in range(6)]

[ wxID_SAVESEQSEL,
] = [wx.NewId() for item in range(1)]

[ wxID_SELECTPHASE,
] = [wx.NewId() for item in range(1)]

[ wxID_PDFCOPYCONTROLS, wxID_PDFSAVECONTROLS, wxID_PDFLOADCONTROLS, 
    wxID_PDFCOMPUTE, wxID_PDFCOMPUTEALL, wxID_PDFADDELEMENT, wxID_PDFDELELEMENT,
] = [wx.NewId() for item in range(7)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

def ShowHelp(helpType,frame):
    '''Called to bring up a web page for documentation.'''
    global htmlFirstUse
    # look up a definition for help info from dict
    helplink = helpLocDict.get(helpType)
    if helplink is None:
        # no defined link to use, create a default based on key
        helplink = 'gsasII.html#'+helpType.replace(' ','_')
    helplink = os.path.join(path2GSAS2,'help',helplink)
    if helpMode == 'internal':
        try:
            htmlPanel.LoadFile(helplink)
            htmlFrame.Raise()
        except:
            htmlFrame = wx.Frame(frame, -1, size=(610, 510))
            htmlFrame.Show(True)
            htmlFrame.SetTitle("HTML Window") # N.B. reset later in LoadFile
            htmlPanel = MyHtmlPanel(htmlFrame,-1)
            htmlPanel.LoadFile(helplink)
    else:
        if htmlFirstUse:
            webbrowser.open_new("file://"+helplink)
            htmlFirstUse = False
        else:
            webbrowser.open("file://"+helplink, new=0, autoraise=True)

class MyHelp(wx.Menu):
    '''This class creates the contents of a help menu.
    The menu will start with two entries:
      'Help on <helpType>': where helpType is a reference to an HTML page to
      be opened
      About: opens an About dialog using OnHelpAbout. N.B. on the Mac this
      gets moved to the App menu to be consistent with Apple style.
    NOTE: the title when appending this menu should be '&Help' so the wx handles
    it correctly. BHT
    '''
    def __init__(self,frame,helpType=None,helpLbl=None,morehelpitems=[],title=''):
        wx.Menu.__init__(self,title)
        self.HelpById = {}
        self.frame = frame
        # add a help item only when helpType is specified
        if helpType is not None:
            if helpLbl is None: helpLbl = helpType
            helpobj = self.Append(text='Help on '+helpLbl,
                                  id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
            frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
            self.HelpById[helpobj.GetId()] = helpType
        for lbl,indx in morehelpitems:
            helpobj = self.Append(text=lbl,
                id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
            frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
            self.HelpById[helpobj.GetId()] = indx
        self.Append(help='', id=wx.ID_ABOUT, kind=wx.ITEM_NORMAL,
            text='&About GSAS-II')
        frame.Bind(wx.EVT_MENU, self.OnHelpAbout, id=wx.ID_ABOUT)
        if GSASIIpath.whichsvn():
            helpobj = self.Append(
                help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
                text='&Check for updates')
            frame.Bind(wx.EVT_MENU, self.OnCheckUpdates, helpobj)
       
    def OnHelpById(self,event):
        '''Called when Help on... is pressed in a menu. Brings up
        a web page for documentation.
        '''
        helpType = self.HelpById.get(event.GetId())
        if helpType is None:
            print 'Error: help lookup failed!',event.GetEventObject()
            print 'id=',event.GetId()
        else:
            ShowHelp(helpType,self.frame)

    def OnHelpAbout(self, event):
        "Display an 'About GSAS-II' box"
        global __version__
        info = wx.AboutDialogInfo()
        info.Name = 'GSAS-II'
        info.Version = __version__ + ' Revision '+str(GSASIIpath.GetVersionNumber())
        info.Copyright = '''
Robert B. Von Dreele & Brian H. Toby
Argonne National Laboratory(C)
This product includes software developed
by the UChicago Argonne, LLC, as 
Operator of Argonne National Laboratory.         '''
        info.Description = '''
General Structure Analysis System - GSAS-II
'''
        wx.AboutBox(info)

    def OnCheckUpdates(self,event):
        '''Check if the GSAS-II repository has an update for the current source files
        and perform that update if requested.
        '''
        if not GSASIIpath.whichsvn():
            dlg = wx.MessageDialog(self,'No Subversion','Cannot update GSAS-II because subversion (svn) '+
                                   'was not found.'
                                   ,wx.OK)
            dlg.ShowModal()
            return
        wx.BeginBusyCursor()
        local = GSASIIpath.svnGetRev()
        if local is None: 
            wx.EndBusyCursor()
            dlg = wx.MessageDialog(self.frame,
                                   'Unable to run subversion on the GSAS-II current directory. Is GSAS-II installed correctly?',
                                   'Subversion error',
                                   wx.OK)
            dlg.ShowModal()
            return
        print 'Installed GSAS-II version: '+local
        repos = GSASIIpath.svnGetRev(local=False)
        wx.EndBusyCursor()
        if repos is None: 
            dlg = wx.MessageDialog(self.frame,
                                   'Unable to access the GSAS-II server. Is this computer on the internet?',
                                   'Server unavailable',
                                   wx.OK)
            dlg.ShowModal()
            return
        print 'GSAS-II version on server: '+repos
        if local == repos:
            dlg = wx.MessageDialog(self.frame,
                                   'GSAS-II is up-to-date. Version '+local+' is already loaded.',
                                   'GSAS-II Up-to-date',
                                   wx.OK)
            dlg.ShowModal()
            return
        mods = GSASIIpath.svnFindLocalChanges()
        if mods:
            dlg = wx.MessageDialog(self.frame,
                                   'You have version '+local+
                                   ' of GSAS-II installed, but the current version is '+repos+
                                   '. However, you have modified '+str(len(mods))+
                                   ' file(s) on your local computer have been modified.'
                                   ' Updating could wipe out your local changes. Press OK to start an update:',
                                   'Local GSAS-II Mods',
                                   wx.OK|wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK: return
        else:
            dlg = wx.MessageDialog(self.frame,
                                   'You have version '+local+
                                   ' of GSAS-II installed, but the current version is '+repos+
                                   '. Press OK to start an update:',
                                   'GSAS-II Updates',
                                   wx.OK|wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK: return
        print 'start updates'
        wx.BeginBusyCursor()
        moddict = GSASIIpath.svnUpdateDir()
        wx.EndBusyCursor()
        if moddict is None: 
            dlg = wx.MessageDialog(self.frame,
                                   'Error accessing the GSAS-II server or performing the update. '+
                                   'Try again later or perform a manual update',
                                   'Update Error',
                                   wx.OK)
            dlg.ShowModal()
            return
        modsbytype = {}
        for key in moddict:
            typ = moddict[key]
            if modsbytype.get(typ) is None:
                modsbytype[typ] = []
            modsbytype[typ].append(key)
        msg = 'Update was completed. Changes will take effect when GSAS-II is next updated. The following files were updated, ordered by status:'
        for key in modsbytype:
            msg += '\n' + key + ':\n\t'
            for fil in modsbytype:
                msg += fil + ', '
        dlg = wx.MessageDialog(self.frame,msg, 'Update Completed', wx.OK)
        dlg.ShowModal()
        return

class AddHelp(wx.Menu):
    '''This class a single entry for the help menu (used on the Mac only):
      'Help on <helpType>': where helpType is a reference to an HTML page to
      be opened
    NOTE: the title when appending this menu should be '&Help' so the wx handles
    it correctly. BHT
    '''
    def __init__(self,frame,helpType,helpLbl=None,title=''):
        wx.Menu.__init__(self,title)
        self.frame = frame
        if helpLbl is None: helpLbl = helpType
        # add a help item only when helpType is specified
        helpobj = self.Append(text='Help on '+helpLbl,
                              id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
        frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
        self.HelpById = helpType
       
    def OnHelpById(self,event):
        '''Called when Help on... is pressed in a menu. Brings up
        a web page for documentation.
        '''
        ShowHelp(self.HelpById,self.frame)

class MyHtmlPanel(wx.Panel):
    '''Defines a panel to display Help information'''
    def __init__(self, frame, id):
        self.frame = frame
        wx.Panel.__init__(self, frame, id)
        sizer = wx.BoxSizer(wx.VERTICAL)
        back = wx.Button(self, -1, "Back")
        back.Bind(wx.EVT_BUTTON, self.OnBack)
        self.htmlwin = G2HtmlWindow(self, id, size=(750,450))
        sizer.Add(self.htmlwin, 1,wx.EXPAND)
        sizer.Add(back, 0, wx.ALIGN_LEFT, 0)
        self.SetSizer(sizer)
        sizer.Fit(frame)        
        self.Bind(wx.EVT_SIZE,self.OnSize)
    def OnSize(self,event):         #does the job but weirdly!!
        anchor = self.htmlwin.GetOpenedAnchor()
        if anchor:            
            self.htmlwin.ScrollToAnchor(anchor)
            wx.CallAfter(self.htmlwin.ScrollToAnchor,anchor)
            event.Skip()
    def OnBack(self, event):
        self.htmlwin.HistoryBack()
    def LoadFile(self,file):
        pos = file.rfind('#')
        if pos != -1:
            helpfile = file[:pos]
            helpanchor = file[pos+1:]
        else:
            helpfile = file
            helpanchor = None
        self.htmlwin.LoadPage(helpfile)
        if helpanchor is not None:
            self.htmlwin.ScrollToAnchor(helpanchor)
            xs,ys = self.htmlwin.GetViewStart()
            self.htmlwin.Scroll(xs,ys-1)

class G2HtmlWindow(wx.html.HtmlWindow):
    '''Displays help information in a primitive HTML browser type window
    '''
    def __init__(self, parent, *args, **kwargs):
        self.parent = parent
        wx.html.HtmlWindow.__init__(self, parent, *args, **kwargs)
    def LoadPage(self, *args, **kwargs):
        wx.html.HtmlWindow.LoadPage(self, *args, **kwargs)
        self.TitlePage()
    def OnLinkClicked(self, *args, **kwargs):
        wx.html.HtmlWindow.OnLinkClicked(self, *args, **kwargs)
        xs,ys = self.GetViewStart()
        self.Scroll(xs,ys-1)
        self.TitlePage()
    def HistoryBack(self, *args, **kwargs):
        wx.html.HtmlWindow.HistoryBack(self, *args, **kwargs)
        self.TitlePage()
    def TitlePage(self):
        self.parent.frame.SetTitle(self.GetOpenedPage() + ' -- ' + 
            self.GetOpenedPageTitle())

class DataFrame(wx.Frame):
    '''Create the dataframe window and its menus
    '''
    def FillDataMenu(self,menu,helpType,helpLbl=None):
        '''Create the "standard" part of data frame menus. Note that on Linux and
        Windows, this is the standard help Menu. On Mac, this menu duplicates the
        tree menu, but adds an extra help command for the data item and a separator. 
        '''
        if sys.platform == "darwin": # mac                         
            menu.Append(AddHelp(self.G2frame,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')
            self.G2frame.FillMainMenu(menu) # add the data tree menu items
            menu.Append(wx.Menu(title=''),title='|') # add a separator
        else: # other
            menu.Append(menu=MyHelp(self,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')

    def _init_menus(self):
        
# define all GSAS-II data frame menus        

# for use where no menu or data frame help is provided
        self.BlankMenu = wx.MenuBar()
        
# Controls
        self.ControlsMenu = wx.MenuBar()
        self.FillDataMenu(self.ControlsMenu,helpType='Controls')
        
# Notebook
        self.DataNotebookMenu = wx.MenuBar()
        self.FillDataMenu(self.DataNotebookMenu,helpType='Notebook')
        
# Comments
        self.DataCommentsMenu = wx.MenuBar()
        self.FillDataMenu(self.DataCommentsMenu,helpType='Comments')
        
# Constraints
        self.ConstraintMenu = wx.MenuBar()
        self.FillDataMenu(self.ConstraintMenu,helpType='Constraints')
        self.ConstraintEdit = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintEdit, title='Edit')
        self.ConstraintEdit.Append(id=wxID_HOLDADD, kind=wx.ITEM_NORMAL,text='Add hold',
            help='Add hold on a parameter value')
        self.ConstraintEdit.Append(id=wxID_EQUIVADD, kind=wx.ITEM_NORMAL,text='Add equivalence',
            help='Add equivalence between parameter values')
        self.ConstraintEdit.Append(id=wxID_CONSTRAINTADD, kind=wx.ITEM_NORMAL,text='Add constraint',
            help='Add constraint on parameter values')
        self.ConstraintEdit.Append(id=wxID_FUNCTADD, kind=wx.ITEM_NORMAL,text='Add New Var',
            help='Add variable composed of existing parameter')
            
# Restraints
        self.RestraintMenu = wx.MenuBar()
        self.FillDataMenu(self.RestraintMenu,helpType='Restraints')
        self.RestraintEdit = wx.Menu(title='')
        self.RestraintMenu.Append(menu=self.RestraintEdit, title='Edit')
        self.RestraintEdit.Append(id=wxID_RESTSELPHASE, kind=wx.ITEM_NORMAL,text='Select phase',
            help='Select phase')
        self.RestraintEdit.Append(id=wxID_RESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add restraints',
            help='Add restraints')
        self.RestraintEdit.Append(id=wxID_RESRCHANGEVAL, kind=wx.ITEM_NORMAL,text='Change value',
            help='Change observed value')
        self.RestraintEdit.Append(id=wxID_RESTCHANGEESD, kind=wx.ITEM_NORMAL,text='Change esd',
            help='Change esd in observed value')
        self.RestraintEdit.Append(id=wxID_RESTDELETE, kind=wx.ITEM_NORMAL,text='Delete restraints',
            help='Delete selected restraints')
            
# Sequential results
        self.SequentialMenu = wx.MenuBar()
        self.FillDataMenu(self.SequentialMenu,helpType='Sequential',helpLbl='Sequential Refinement')
        self.SequentialFile = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialFile, title='File')
        self.SequentialFile.Append(id=wxID_SAVESEQSEL, kind=wx.ITEM_NORMAL,text='Save...',
            help='Save selected sequential refinement results')
            
# PDR
        self.ErrorMenu = wx.MenuBar()
        self.FillDataMenu(self.ErrorMenu,helpType='PWD Analysis',helpLbl='Powder Fit Error Analysis')
        self.ErrorAnal = wx.Menu(title='')
        self.ErrorMenu.Append(menu=self.ErrorAnal,title='Analysis')
        self.ErrorAnal.Append(id=wxID_PWDANALYSIS,kind=wx.ITEM_NORMAL,text='Analyze',
            help='Error analysis on powder pattern')
            
# PDR / Limits
        self.LimitMenu = wx.MenuBar()
        self.FillDataMenu(self.LimitMenu,helpType='Limits')
        self.LimitEdit = wx.Menu(title='')
        self.LimitMenu.Append(menu=self.LimitEdit, title='File')
        self.LimitEdit.Append(id=wxID_LIMITCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy limits to other histograms')
            
# PDR / Background
        self.BackMenu = wx.MenuBar()
        self.FillDataMenu(self.BackMenu,helpType='Background')
        self.BackEdit = wx.Menu(title='')
        self.BackMenu.Append(menu=self.BackEdit, title='File')
        self.BackEdit.Append(id=wxID_BACKCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy background parameters to other histograms')
        self.BackEdit.Append(id=wxID_BACKFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy background refinement flags to other histograms')
            
# PDR / Instrument Parameters
        self.InstMenu = wx.MenuBar()
        self.FillDataMenu(self.InstMenu,helpType='Instrument Parameters')
        self.InstEdit = wx.Menu(title='')
        self.InstMenu.Append(menu=self.InstEdit, title='Operations')
        self.InstEdit.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTLOAD, kind=wx.ITEM_NORMAL,text='Load profile...')
        self.InstEdit.Append(help='Load instrument profile parameters from file', 
            id=wxID_INSTSAVE, kind=wx.ITEM_NORMAL,text='Save profile...')
        self.InstEdit.Append(help='Save instrument profile parameters to file', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')
        self.InstEdit.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')
        self.InstEdit.Append(id=wxID_INSTFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy instrument parameter refinement flags to other histograms')
        self.InstEdit.Append(help='Change radiation type (Ka12 - synch)', 
            id=wxID_CHANGEWAVETYPE, kind=wx.ITEM_NORMAL,text='Change radiation')
        
# PDR / Sample Parameters
        self.SampleMenu = wx.MenuBar()
        self.FillDataMenu(self.SampleMenu,helpType='Sample Parameters')
        self.SampleEdit = wx.Menu(title='')
        self.SampleMenu.Append(menu=self.SampleEdit, title='File')
        self.SampleEdit.Append(id=wxID_SAMPLELOAD, kind=wx.ITEM_NORMAL,text='Load',
            help='Load sample parameters from file')
        self.SampleEdit.Append(id=wxID_SAMPLESAVE, kind=wx.ITEM_NORMAL,text='Save',
            help='Save sample parameters to file')
        self.SampleEdit.Append(id=wxID_SAMPLECOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy refinable sample parameters to other histograms')
        self.SampleEdit.Append(id=wxID_SAMPLEFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy sample parameter refinement flags to other histograms')

# PDR / Peak List
        self.PeakMenu = wx.MenuBar()
        self.FillDataMenu(self.PeakMenu,helpType='Peak List')
        self.PeakEdit = wx.Menu(title='')
        self.PeakMenu.Append(menu=self.PeakEdit, title='Peak Fitting')
        self.UnDo = self.PeakEdit.Append(help='Undo last least squares refinement', 
            id=wxID_UNDO, kind=wx.ITEM_NORMAL,text='UnDo')
        self.PeakFit = self.PeakEdit.Append(id=wxID_LSQPEAKFIT, kind=wx.ITEM_NORMAL,text='LSQ PeakFit', 
            help='Peak fitting via least-squares' )
        self.PFOneCycle = self.PeakEdit.Append(id=wxID_LSQONECYCLE, kind=wx.ITEM_NORMAL,text='LSQ one cycle', 
            help='One cycle of Peak fitting via least-squares' )
        self.PeakEdit.Append(id=wxID_RESETSIGGAM, kind=wx.ITEM_NORMAL, 
            text='Reset sig and gam',help='Reset sigma and gamma to global fit' )
        self.PeakEdit.Append(id=wxID_CLEARPEAKS, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the peak list' )
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.PFOneCycle.Enable(False)
        
# PDR / Index Peak List
        self.IndPeaksMenu = wx.MenuBar()
        self.FillDataMenu(self.IndPeaksMenu,helpType='Index Peak List')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndPeaksMenu.Append(menu=self.IndPeaksEdit,title='Operations')
        self.IndPeaksEdit.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
        
# PDR / Unit Cells List
        self.IndexMenu = wx.MenuBar()
        self.FillDataMenu(self.IndexMenu,helpType='Unit Cells List')
        self.IndexEdit = wx.Menu(title='')
        self.IndexMenu.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        self.IndexPeaks = self.IndexEdit.Append(help='', id=wxID_INDEXPEAKS, kind=wx.ITEM_NORMAL,
            text='Index Cell')
        self.CopyCell = self.IndexEdit.Append( id=wxID_COPYCELL, kind=wx.ITEM_NORMAL,text='Copy Cell', 
            help='Copy selected unit cell from indexing to cell refinement fields')
        self.RefineCell = self.IndexEdit.Append( id=wxID_REFINECELL, kind=wx.ITEM_NORMAL, 
            text='Refine Cell',help='Refine unit cell parameters from indexed peaks')
        self.MakeNewPhase = self.IndexEdit.Append( id=wxID_MAKENEWPHASE, kind=wx.ITEM_NORMAL,
            text='Make new phase',help='Make new phase from selected unit cell')
        self.IndexPeaks.Enable(False)
        self.CopyCell.Enable(False)
        self.RefineCell.Enable(False)
        self.MakeNewPhase.Enable(False)
        
# PDR / Reflection Lists
        self.ReflMenu = wx.MenuBar()
        self.FillDataMenu(self.ReflMenu,helpType='Reflection List')
        self.ReflEdit = wx.Menu(title='')
        self.ReflMenu.Append(menu=self.ReflEdit, title='Reflection List')
        self.SelectPhase = self.ReflEdit.Append(help='Select phase for reflection list',id=wxID_SELECTPHASE, 
            kind=wx.ITEM_NORMAL,text='Select phase')
        
# IMG / Image Controls
        self.ImageMenu = wx.MenuBar()
        self.FillDataMenu(self.ImageMenu,helpType='Image Controls')
        self.ImageEdit = wx.Menu(title='')
        self.ImageMenu.Append(menu=self.ImageEdit, title='Operations')
        self.ImageEdit.Append(help='Calibrate detector by fitting to calibrant lines', 
            id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        self.ImageEdit.Append(help='Recalibrate detector by fitting to calibrant lines', 
            id=wxID_IMRECALIBRATE, kind=wx.ITEM_NORMAL,text='Recalibrate')
        self.ImageEdit.Append(help='Clear calibration data points and rings',id=wxID_IMCLEARCALIB, 
            kind=wx.ITEM_NORMAL,text='Clear calibration')
        self.ImageEdit.Append(help='Integrate selected image',id=wxID_IMINTEGRATE, 
            kind=wx.ITEM_NORMAL,text='Integrate')
        self.ImageEdit.Append(help='Integrate all images selected from list',id=wxID_INTEGRATEALL,
            kind=wx.ITEM_NORMAL,text='Integrate all')
        self.ImageEdit.Append(help='Copy image controls to other images', 
            id=wxID_IMCOPYCONTROLS, kind=wx.ITEM_NORMAL,text='Copy Controls')
        self.ImageEdit.Append(help='Save image controls to file', 
            id=wxID_IMSAVECONTROLS, kind=wx.ITEM_NORMAL,text='Save Controls')
        self.ImageEdit.Append(help='Load image controls from file', 
            id=wxID_IMLOADCONTROLS, kind=wx.ITEM_NORMAL,text='Load Controls')
            
# IMG / Masks
        self.MaskMenu = wx.MenuBar()
        self.FillDataMenu(self.MaskMenu,helpType='Image Masks')
        self.MaskEdit = wx.Menu(title='')
        self.MaskMenu.Append(menu=self.MaskEdit, title='Operations')
        self.MaskEdit.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')
        self.MaskEdit.Append(help='Save mask to file', 
            id=wxID_MASKSAVE, kind=wx.ITEM_NORMAL,text='Save mask')
        self.MaskEdit.Append(help='Load mask from file', 
            id=wxID_MASKLOAD, kind=wx.ITEM_NORMAL,text='Load mask')
            
# IMG / Stress/Strain

        self.StrStaMenu = wx.MenuBar()
        self.FillDataMenu(self.StrStaMenu,helpType='Stress/Strain')
        self.StrStaEdit = wx.Menu(title='')
        self.StrStaMenu.Append(menu=self.StrStaEdit, title='Operations')
        self.StrStaEdit.Append(help='Append d-zero for one ring', 
            id=wxID_APPENDDZERO, kind=wx.ITEM_NORMAL,text='Append d-zero')
        self.StrStaEdit.Append(help='Fit stress/strain data', 
            id=wxID_STRSTAFIT, kind=wx.ITEM_NORMAL,text='Fit stress/strain')
        self.StrStaEdit.Append(help='Copy stress/strain data to other images', 
            id=wxID_STRSTACOPY, kind=wx.ITEM_NORMAL,text='Copy stress/strain')
        self.StrStaEdit.Append(help='Save stress/strain data to file', 
            id=wxID_STRSTASAVE, kind=wx.ITEM_NORMAL,text='Save stress/strain')
        self.StrStaEdit.Append(help='Load stress/strain data from file', 
            id=wxID_STRSTALOAD, kind=wx.ITEM_NORMAL,text='Load stress/strain')
            
# PDF / PDF Controls
        self.PDFMenu = wx.MenuBar()
        self.FillDataMenu(self.PDFMenu,helpType='PDF Controls')
        self.PDFEdit = wx.Menu(title='')
        self.PDFMenu.Append(menu=self.PDFEdit, title='PDF Controls')
        self.PDFEdit.Append(help='Add element to sample composition',id=wxID_PDFADDELEMENT, kind=wx.ITEM_NORMAL,
            text='Add element')
        self.PDFEdit.Append(help='Delete element from sample composition',id=wxID_PDFDELELEMENT, kind=wx.ITEM_NORMAL,
            text='Delete element')
        self.PDFEdit.Append(help='Copy PDF controls', id=wxID_PDFCOPYCONTROLS, kind=wx.ITEM_NORMAL,
            text='Copy controls')
#        self.PDFEdit.Append(help='Load PDF controls from file',id=wxID_PDFLOADCONTROLS, kind=wx.ITEM_NORMAL,
#            text='Load Controls')
#        self.PDFEdit.Append(help='Save PDF controls to file', id=wxID_PDFSAVECONTROLS, kind=wx.ITEM_NORMAL,
#            text='Save controls')
        self.PDFEdit.Append(help='Compute PDF', id=wxID_PDFCOMPUTE, kind=wx.ITEM_NORMAL,
            text='Compute PDF')
        self.PDFEdit.Append(help='Compute all PDFs', id=wxID_PDFCOMPUTEALL, kind=wx.ITEM_NORMAL,
            text='Compute all PDFs')
            
# Phase / General tab

        self.DataGeneral = wx.MenuBar()
        self.FillDataMenu(self.DataGeneral,helpType='General', helpLbl='Phase/General')
        self.GeneralCalc = wx.Menu(title='')
        self.DataGeneral.Append(menu=self.GeneralCalc,title='Compute')
        self.GeneralCalc.Append(help='Compute Fourier map',id=wxID_FOURCALC, kind=wx.ITEM_NORMAL,
            text='Fourier map')
        self.GeneralCalc.Append(help='Search Fourier map',id=wxID_FOURSEARCH, kind=wx.ITEM_NORMAL,
            text='Search map')
        self.GeneralCalc.Append(help='Run charge flipping',id=wxID_CHARGEFLIP, kind=wx.ITEM_NORMAL,
            text='Charge flipping')
        self.GeneralCalc.Append(help='Clear map',id=wxID_FOURCLEAR, kind=wx.ITEM_NORMAL,
            text='Clear map')
        
# Phase / Data tab
        self.DataMenu = wx.MenuBar()
        self.FillDataMenu(self.DataMenu,helpType='Data', helpLbl='Phase/Data')
        self.DataEdit = wx.Menu(title='')
        self.DataMenu.Append(menu=self.DataEdit, title='Edit')
        self.DataEdit.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Delete histograms',
            help='Delete histograms from use for this phase')
            
# Phase / Atoms tab
        self.AtomsMenu = wx.MenuBar()
        self.FillDataMenu(self.AtomsMenu,helpType='Atoms')
        self.AtomEdit = wx.Menu(title='')
        self.AtomCompute = wx.Menu(title='')
        self.AtomsMenu.Append(menu=self.AtomEdit, title='Edit')
        self.AtomsMenu.Append(menu=self.AtomCompute, title='Compute')
        self.AtomEdit.Append(id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append atom',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSVIEWADD, kind=wx.ITEM_NORMAL,text='Append view point',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert atom',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMVIEWINSERT, kind=wx.ITEM_NORMAL,text='Insert view point',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSEDITDELETE, kind=wx.ITEM_NORMAL,text='Delete atom',
            help='Select atoms to delete first')
        self.AtomEdit.Append(id=wxID_ATOMSREFINE, kind=wx.ITEM_NORMAL,text='Set atom refinement flags',
            help='Select atoms to refine first')
        self.AtomEdit.Append(id=wxID_ATOMSMODIFY, kind=wx.ITEM_NORMAL,text='Modify atom parameters',
            help='Select atoms to modify first')
        self.AtomEdit.Append(id=wxID_ATOMSTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Select atoms to transform first')
        self.AtomEdit.Append(id=wxID_RELOADDRAWATOMS, kind=wx.ITEM_NORMAL,text='Reload draw atoms',
            help='Reload atom drawing list')
        self.AtomCompute.Append(id=wxID_ATOMSDISAGL, kind=wx.ITEM_NORMAL,text='Distances && Angles',
            help='Compute distances & angles for selected atoms')
                 
# Phase / Draw Options tab
        self.DataDrawOptions = wx.MenuBar()
        self.FillDataMenu(self.DataDrawOptions,helpType='Draw Options', helpLbl='Phase/Draw Options')
        
# Phase / Draw Atoms tab
        self.DrawAtomsMenu = wx.MenuBar()
        self.FillDataMenu(self.DrawAtomsMenu,helpType='Draw Atoms')
        self.DrawAtomEdit = wx.Menu(title='')
        self.DrawAtomCompute = wx.Menu(title='')
        self.DrawAtomRestraint = wx.Menu(title='')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomEdit, title='Edit')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomCompute,title='Compute')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRestraint, title='Restraints')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMSTYLE, kind=wx.ITEM_NORMAL,text='Atom style',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMLABEL, kind=wx.ITEM_NORMAL,text='Atom label',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMCOLOR, kind=wx.ITEM_NORMAL,text='Atom color',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMRESETCOLOR, kind=wx.ITEM_NORMAL,text='Reset atom colors',
            help='Resets all atom colors to defaults')
        self.DrawAtomEdit.Append(id=wxID_DRAWVIEWPOINT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st atom selected')
        self.DrawAtomEdit.Append(id=wxID_DRAWADDEQUIV, kind=wx.ITEM_NORMAL,text='Add atoms',
            help='Add symmetry & cell equivalents to drawing set from selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Transform selected atoms by symmetry & cell translations')
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCOORD, kind=wx.ITEM_NORMAL,text='Fill CN-sphere',
            help='Fill coordination sphere for selected atoms')            
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCELL, kind=wx.ITEM_NORMAL,text='Fill unit cell',
            help='Fill unit cell with selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWDELETE, kind=wx.ITEM_NORMAL,text='Delete atoms',
            help='Delete atoms from drawing set')
        self.DrawAtomCompute.Append(id=wxID_DRAWDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected atoms from view point')   
        self.DrawAtomCompute.Append(id=wxID_DRAWDISAGLTOR, kind=wx.ITEM_NORMAL,text='Dist. Ang. Tors.',
            help='Compute distance, angle or torsion for 2-4 selected atoms')   
        self.DrawAtomCompute.Append(id=wxID_DRAWPLANE, kind=wx.ITEM_NORMAL,text='Best plane',
            help='Compute best plane for 4+ selected atoms')   
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRBOND, kind=wx.ITEM_NORMAL,text='Add bond restraint',
            help='Add bond restraint for selected atoms (2)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRANGLE, kind=wx.ITEM_NORMAL,text='Add angle restraint',
            help='Add angle restraint for selected atoms (3: one end 1st)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRPLANE, kind=wx.ITEM_NORMAL,text='Add plane restraint',
            help='Add plane restraint for selected atoms (4+)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRCHIRAL, kind=wx.ITEM_NORMAL,text='Add chiral restraint',
            help='Add chiral restraint for selected atoms (4: center atom 1st)')
            
# Phase / Texture tab
        self.TextureMenu = wx.MenuBar()
        self.FillDataMenu(self.TextureMenu,helpType='Texture')
        self.TextureEdit = wx.Menu(title='')
        self.TextureMenu.Append(menu=self.TextureEdit, title='Texture')
        self.TextureEdit.Append(id=wxID_REFINETEXTURE, kind=wx.ITEM_NORMAL,text='Refine texture', 
            help='Refine the texture coefficients from sequential Pawley results')
        self.TextureEdit.Append(id=wxID_CLEARTEXTURE, kind=wx.ITEM_NORMAL,text='Clear texture', 
            help='Clear the texture coefficients' )
            
# Phase / Pawley tab
        self.PawleyMenu = wx.MenuBar()
        self.FillDataMenu(self.PawleyMenu,helpType='Pawley')
        self.PawleyEdit = wx.Menu(title='')
        self.PawleyMenu.Append(menu=self.PawleyEdit,title='Operations')
        self.PawleyEdit.Append(id=wxID_PAWLEYLOAD, kind=wx.ITEM_NORMAL,text='Pawley create',
            help='Initialize Pawley reflection list')
        self.PawleyEdit.Append(id=wxID_PAWLEYESTIMATE, kind=wx.ITEM_NORMAL,text='Pawley estimate',
            help='Estimate initial Pawley intensities')
        self.PawleyEdit.Append(id=wxID_PAWLEYUPDATE, kind=wx.ITEM_NORMAL,text='Pawley update',
            help='Update Pawley intensities with abs(Fobs) from reflection list')
#        self.PawleyEdit.Append(id=wxID_PAWLEYDELETE, kind=wx.ITEM_NORMAL,text='Pawley delete',
#            help='Delete selected Pawley reflection')
            
# Phase / Map peaks tab
        self.MapPeaksMenu = wx.MenuBar()
        self.FillDataMenu(self.MapPeaksMenu,helpType='Map peaks')
        self.MapPeaksEdit = wx.Menu(title='')
        self.MapPeaksMenu.Append(menu=self.MapPeaksEdit, title='Map peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSMOVE, kind=wx.ITEM_NORMAL,text='Move peaks', 
            help='Move selected peaks to atom list')
        self.MapPeaksEdit.Append(id=wxID_PEAKSVIEWPT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st peak selected')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected peaks from view point')   
        self.MapPeaksEdit.Append(id=wxID_PEAKSDA, kind=wx.ITEM_NORMAL,text='Calc dist/ang', 
            help='Calculate distance or angle for selection')
        self.MapPeaksEdit.Append(id=wxID_FINDEQVPEAKS, kind=wx.ITEM_NORMAL,text='Equivalent peaks', 
            help='Find equivalent peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSUNIQUE, kind=wx.ITEM_NORMAL,text='Unique peaks', 
            help='Select unique set')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDELETE, kind=wx.ITEM_NORMAL,text='Delete peaks', 
            help='Delete selected peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSCLEAR, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the map peak list')
            
# end of GSAS-II menu definitions
        
    def _init_ctrls(self, parent,name=None,size=None,pos=None):
        wx.Frame.__init__(self,parent=parent,
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX | wx.FRAME_FLOAT_ON_PARENT ,
            size=size,pos=pos,title='GSAS-II data display')
        self._init_menus()
        if name:
            self.SetLabel(name)
        self.Show()
        
    def __init__(self,parent,frame,data=None,name=None, size=None,pos=None):
        self.G2frame = frame
        self._init_ctrls(parent,name,size,pos)
        self.data = data
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]+250))
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        clientSize = wx.ClientDisplayRect()
        Width[1] = min(Width[1],clientSize[2]-300)
        Width[0] = max(Width[0],300)
        self.SetSize(Width)
#        self.SetPosition(wx.Point(clientSize[2]-Width[0],clientSize[1]+250))
        
    def Clear(self):
        self.ClearBackground()
        self.DestroyChildren()
                   
################################################################################
#####  GSNotebook
################################################################################           
       
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
        
################################################################################
#####  GSGrid
################################################################################           
       
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
                        
################################################################################
#####  Table
################################################################################           
       
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
        
    def GetColValues(self, col):
        data = []
        for row in range(self.GetNumberRows()):
            data.append(self.GetValue(row, col))
        return data
        
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
                
################################################################################
#####  Notebook
################################################################################           
       
def UpdateNotebook(G2frame,data):        
    if data:
        G2frame.dataFrame.SetLabel('Notebook')
        G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
        for line in data:
            G2frame.dataDisplay.AppendText(line+"\n")
            G2frame.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
            
################################################################################
#####  Controls
################################################################################           
       
def UpdateControls(G2frame,data):
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytic Hessian'
        data['min dM/M'] = 0.0001
        data['shift factor'] = 1.
        data['max cyc'] = 3        
        data['F**2'] = True
        data['minF/sig'] = 0
    if 'shift factor' not in data:
        data['shift factor'] = 1.
    if 'max cyc' not in data:
        data['max cyc'] = 3
    if 'F**2' not in data:
        data['F**2'] = True
        data['minF/sig'] = 0
    #end patch

    def SeqSizer():
        
        def OnSelectData(event):
            choices = ['All',]+GetPatternTreeDataNames(G2frame,['PWDR',])
            sel = []
            if 'Seq Data' in data:
                for item in data['Seq Data']:
                    sel.append(choices.index(item))
            names = []
            dlg = wx.MultiChoiceDialog(G2frame,'Select data:','Sequential refinement',choices)
            dlg.SetSelections(sel)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for i in sel: names.append(choices[i])
                if 'All' in names:
                    names = choices[1:]
                data['Seq Data'] = names                
            dlg.Destroy()
            reverseSel.Enable(True)
            
        def OnReverse(event):
            data['Reverse Seq'] = reverseSel.GetValue()
                    
        seqSizer = wx.BoxSizer(wx.HORIZONTAL)
        seqSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sequential Refinement Powder Data: '),0,wx.ALIGN_CENTER_VERTICAL)
        selSeqData = wx.Button(G2frame.dataDisplay,-1,label=' Select data')
        selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
        seqSizer.Add(selSeqData,0,wx.ALIGN_CENTER_VERTICAL)
        seqSizer.Add((5,0),0)
        reverseSel = wx.CheckBox(G2frame.dataDisplay,-1,label=' Reverse order?')
        reverseSel.Bind(wx.EVT_CHECKBOX,OnReverse)
        if 'Seq Data' not in data:
            reverseSel.Enable(False)
        if 'Reverse Seq' in data:
            reverseSel.SetValue(data['Reverse Seq'])
        seqSizer.Add(reverseSel,0,wx.ALIGN_CENTER_VERTICAL)
        return seqSizer
        
    def LSSizer():        
        
        def OnDerivType(event):
            data['deriv type'] = derivSel.GetValue()
            derivSel.SetValue(data['deriv type'])
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnConvergence(event):
            try:
                value = max(1.e-9,min(1.0,float(Cnvrg.GetValue())))
            except ValueError:
                value = 0.0001
            data['min dM/M'] = value
            Cnvrg.SetValue('%.2g'%(value))
            
        def OnMaxCycles(event):
            data['max cyc'] = int(maxCyc.GetValue())
            maxCyc.SetValue(str(data['max cyc']))
                        
        def OnFactor(event):
            try:
                value = min(max(float(Factr.GetValue()),0.00001),100.)
            except ValueError:
                value = 1.0
            data['shift factor'] = value
            Factr.SetValue('%.5f'%(value))
            
        def OnFsqRef(event):
            data['F**2'] = fsqRef.GetValue()
        
        def OnMinSig(event):
            try:
                value = min(max(float(minSig.GetValue()),0.),5.)
            except ValueError:
                value = 1.0
            data['minF/sig'] = value
            minSig.SetValue('%.2f'%(value))

        LSSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement derivatives: '),0,wx.ALIGN_CENTER_VERTICAL)
        Choice=['analytic Jacobian','numeric','analytic Hessian']
        derivSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['deriv type'],choices=Choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        derivSel.SetValue(data['deriv type'])
        derivSel.Bind(wx.EVT_COMBOBOX, OnDerivType)
            
        LSSizer.Add(derivSel,0,wx.ALIGN_CENTER_VERTICAL)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min delta-M/M: '),0,wx.ALIGN_CENTER_VERTICAL)
        Cnvrg = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2g'%(data['min dM/M']),style=wx.TE_PROCESS_ENTER)
        Cnvrg.Bind(wx.EVT_TEXT_ENTER,OnConvergence)
        Cnvrg.Bind(wx.EVT_KILL_FOCUS,OnConvergence)
        LSSizer.Add(Cnvrg,0,wx.ALIGN_CENTER_VERTICAL)
        if 'Hessian' in data['deriv type']:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Max cycles: '),0,wx.ALIGN_CENTER_VERTICAL)
            Choice = ['0','1','2','3','5','10','15','20']
            maxCyc = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['max cyc']),choices=Choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            maxCyc.SetValue(str(data['max cyc']))
            maxCyc.Bind(wx.EVT_COMBOBOX, OnMaxCycles)
            LSSizer.Add(maxCyc,0,wx.ALIGN_CENTER_VERTICAL)
        else:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Initial shift factor: '),0,wx.ALIGN_CENTER_VERTICAL)
            Factr = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.5f'%(data['shift factor']),style=wx.TE_PROCESS_ENTER)
            Factr.Bind(wx.EVT_TEXT_ENTER,OnFactor)
            Factr.Bind(wx.EVT_KILL_FOCUS,OnFactor)
            LSSizer.Add(Factr,0,wx.ALIGN_CENTER_VERTICAL)
        if G2frame.Sngl:
            LSSizer.Add((1,0),)
            LSSizer.Add((1,0),)
            fsqRef = wx.CheckBox(G2frame.dataDisplay,-1,label='Refine HKLF as F^2? ')
            fsqRef.SetValue(data['F**2'])
            fsqRef.Bind(wx.EVT_CHECKBOX,OnFsqRef)
            LSSizer.Add(fsqRef,0,wx.ALIGN_CENTER_VERTICAL)
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label='Min obs/sig (0-5): '),0,wx.ALIGN_CENTER_VERTICAL)
            minSig = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2f'%(data['minF/sig']),style=wx.TE_PROCESS_ENTER)
            minSig.Bind(wx.EVT_TEXT_ENTER,OnMinSig)
            minSig.Bind(wx.EVT_KILL_FOCUS,OnMinSig)
            LSSizer.Add(minSig,0,wx.ALIGN_CENTER_VERTICAL)
        return LSSizer
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText('')
    G2frame.dataFrame.SetLabel('Controls')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    SetDataMenuBar(G2frame,G2frame.dataFrame.ControlsMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement Controls:'),0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(LSSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(SeqSizer())
    mainSizer.Add((5,5),0)
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
     
################################################################################
#####  Comments
################################################################################           
       
def UpdateComments(G2frame,data):                   
    G2frame.dataFrame.SetLabel('Comments')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    for line in data:
        if line[-1] == '\n':
            G2frame.dataDisplay.AppendText(line)
        else:
            G2frame.dataDisplay.AppendText(line+'\n')
            
################################################################################
#####  Sequential Results
################################################################################           
       
def UpdateSeqResults(G2frame,data):
    """ 
    input:
        data - dictionary
            'histNames' - list of histogram names in order as processed by Sequential Refinement
            'varyList' - list of variables - identical over all refinements insequence
            histName - dictionaries for all data sets processed:
                'variables'- result[0] from leastsq call
                'varyList' - list of variables; same as above
                'sig' - esds for variables
                'covMatrix' - covariance matrix from individual refinement
                'title' - histogram name; same as dict item name
                'newAtomDict' - new atom parameters after shifts applied
                'newCellDict' - new cell parameters after shifts to A0-A5 applied'
    """
    if not data:
        print 'No sequential refinement results'
        return
    histNames = data['histNames']
       
    def GetSampleParms():
        sampleParmDict = {'Temperature':[],'Pressure':[],'Humidity':[],'Voltage':[],'Force':[],}
        sampleParm = {}
        for name in histNames:
            Id = GetPatternTreeItemId(G2frame,G2frame.root,name)
            sampleData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
            for item in sampleParmDict:
                sampleParmDict[item].append(sampleData[item])
        for item in sampleParmDict:
            frstValue = sampleParmDict[item][0]
            if np.any(np.array(sampleParmDict[item])-frstValue):
                sampleParm[item] = sampleParmDict[item]            
        return sampleParm
            
    def GetRwps():
        Rwps = []
        for name in histNames:
            Rwps.append(data[name]['Rvals']['Rwp'])
        return Rwps
            
    def GetSigData(parm):
        sigData = []
        for name in histNames:
            sigList = data[name]['sig']
            if colLabels[parm] in atomList:
                sigData.append(sigList[colLabels.index(atomList[colLabels[parm]])])
            elif colLabels[parm] in cellList:
                sigData.append(sigList[colLabels.index(cellList[colLabels[parm]])])
            else:
                sigData.append(sigList[parm])
        return sigData
    
    def Select(event):
        cols = G2frame.dataDisplay.GetSelectedCols()
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            plotData = []
            plotSig = []
            plotNames = []
            for col in cols:
                plotData.append(G2frame.SeqTable.GetColValues(col))
                plotSig.append(GetSigData(col))
                plotNames.append(G2frame.SeqTable.GetColLabelValue(col))
            plotData = np.array(plotData)
            G2plt.PlotSeq(G2frame,plotData,plotSig,plotNames,sampleParms)
        elif rows:
            name = histNames[rows[0]]
            G2plt.PlotCovariance(G2frame,Data=data[name])
            
    def OnSaveSelSeq(event):        
        cols = G2frame.dataDisplay.GetSelectedCols()
        if cols:
            numRows = G2frame.SeqTable.GetNumberRows()
            dataNames = []
            saveNames = [G2frame.SeqTable.GetRowLabelValue(r) for r in range(numRows)]
            saveData = []
            for col in cols:
                dataNames.append(G2frame.SeqTable.GetColLabelValue(col))
                saveData.append(zip(G2frame.SeqTable.GetColValues(col),GetSigData(col)))
            lenName = len(saveNames[0])
            saveData = np.swapaxes(np.array(saveData),0,1)
            dlg = wx.FileDialog(G2frame, 'Choose text output file for your selection', '.', '', 
                'Text output file (*.txt)|*.txt',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    SeqTextFile = dlg.GetPath()
                    SeqTextFile = G2IO.FileDlgFixExt(dlg,SeqTextFile)
                    SeqFile = open(SeqTextFile,'w')
                    line = '  %s  '%('name'.center(lenName))
                    for item in dataNames:
                        line += ' %12s %12s '%(item.center(12),'esd'.center(12))
                    line += '\n'
                    SeqFile.write(line)
                    for i,item in enumerate(saveData):
                        line = " '%s' "%(saveNames[i])
                        for val,esd in item:
                            line += ' %12.6f %12.6f '%(val,esd)
                        line += '\n'
                        SeqFile.write(line)
                    SeqFile.close()
            finally:
                dlg.Destroy()
            
               
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    cellList = {}
    newCellDict = data[histNames[0]]['newCellDict']
    for item in newCellDict:
        if item in data['varyList']:
            cellList[newCellDict[item][0]] = item
    atomList = {}
    newAtomDict = data[histNames[0]]['newAtomDict']
    for item in newAtomDict:
        if item in data['varyList']:
            atomList[newAtomDict[item][0]] = item
    sampleParms = GetSampleParms()
    Rwps = GetRwps()
    SetDataMenuBar(G2frame,G2frame.dataFrame.SequentialMenu)
    G2frame.dataFrame.SetLabel('Sequental refinement results')
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=wxID_SAVESEQSEL)
    colLabels = ['Rwp',]+data['varyList']+atomList.keys()+cellList.keys()
    Types = (len(data['varyList']+atomList.keys()+cellList.keys())+1)*[wg.GRID_VALUE_FLOAT,]
    seqList = [[Rwps[i],]+list(data[name]['variables']) for i,name in enumerate(histNames)]    
    for i,item in enumerate(seqList):
        newAtomDict = data[histNames[i]]['newAtomDict']
        newCellDict = data[histNames[i]]['newCellDict']
        item += [newAtomDict[atomList[parm]][1] for parm in atomList.keys()]
        item += [newCellDict[cellList[parm]][1] for parm in cellList.keys()]
    G2frame.SeqTable = Table(seqList,colLabels=colLabels,rowLabels=histNames,types=Types)
    G2frame.dataDisplay = GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.SeqTable, True)
    G2frame.dataDisplay.EnableEditing(False)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, Select)
    G2frame.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(True)
    G2frame.dataFrame.setSizePosLeft([700,350])
    
################################################################################
#####  Constraints
################################################################################           
       
def UpdateConstraints(G2frame,data):
    '''Called when Constraints tree item is selected.
    Displays the constraints in the data window
    '''
    if not data:
        data.update({'Hist':[],'HAP':[],'Phase':[]})       #empty dict - fill it
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    AtomDict = dict([Phases[phase]['pId'],Phases[phase]['Atoms']] for phase in Phases)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable = G2str.GetPhaseData(Phases,Print=False)
    phaseList = []
    for item in phaseDict:
        if item.split(':')[2] not in ['Ax','Ay','Az','Amul','AI/A','Atype','SHorder']:
            phaseList.append(item)
    phaseList.sort()
    phaseAtNames = {}
    for item in phaseList:
        Split = item.split(':')
        if Split[2][:2] in ['AU','Af','dA']:
            phaseAtNames[item] = AtomDict[int(Split[0])][int(Split[3])][0]
        else:
            phaseAtNames[item] = ''
            
    hapVary,hapDict,controlDict = G2str.GetHistogramPhaseData(Phases,Histograms,Print=False)
    hapList = hapDict.keys()
    hapList.sort()
    histVary,histDict,controlDict = G2str.GetHistogramData(Histograms,Print=False)
    histList = []
    for item in histDict:
        if item.split(':')[2] not in ['Omega','Type','Chi','Phi','Azimuth','Gonio. radius','Lam1','Lam2','Back']:
            histList.append(item)
    histList.sort()
    Indx = {}
    scope = {}                          #filled out later
    G2frame.Page = [0,'phs']
    
    def GetPHlegends(Phases,Histograms):
        plegend = '\n In p::name'
        hlegend = '\n In :h:name'
        phlegend = '\n In p:h:name'
        for phase in Phases:
            plegend += '\n p:: = '+str(Phases[phase]['pId'])+':: for '+phase
            count = 0
            for histogram in Phases[phase]['Histograms']:
                if count < 3:
                    phlegend += '\n p:h: = '+str(Phases[phase]['pId'])+':'+str(Histograms[histogram]['hId'])+': for '+phase+' in '+histogram
                else:
                    phlegend += '\n ... etc.'
                    break
                count += 1
        count = 0
        for histogram in Histograms:
            if count < 3:
                hlegend += '\n :h: = :'+str(Histograms[histogram]['hId'])+': for '+histogram
            else:
                hlegend += '\n ... etc.'
                break
            count += 1
        return plegend,hlegend,phlegend
        
    def FindEquivVarb(name,nameList):
        outList = []
        namelist = [name.split(':')[2],]
        if 'dA' in name:
            namelist = ['dAx','dAy','dAz']
        elif 'AU' in name:
            namelist = ['AUiso','AU11','AU22','AU33','AU12','AU13','AU23']
        for item in nameList:
            key = item.split(':')[2]
            if key in namelist and item != name:
                outList.append(item)
        return outList
        
    def SelectVarbs(page,FrstVarb,varList,legend,constType):
        '''Select variables used in Constraints after one variable has
        been selected which determines the appropriate variables to be
        used here. Then creates the constraint and adds it to the
        constraints list.
        Called from OnAddEquivalence, OnAddFunction & OnAddConstraint
        '''
        #future -  add 'all:all:name', '0:all:name', etc. to the varList
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in varList]
            dlg = wx.MultiChoiceDialog(G2frame,'Select more variables:'+legend,
                'Constrain '+FrstVarb+' and...',atchoice)
        else:
            dlg = wx.MultiChoiceDialog(G2frame,'Select more variables:'+legend,
                'Constrain '+FrstVarb+' and...',varList)
        varbs = [FrstVarb,]
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                varbs.append(varList[x])
        dlg.Destroy()
        if len(varbs) > 1:
            if 'equivalence' in constType:
                constr = [[1.0,FrstVarb]]
                for item in varbs[1:]:
                    constr += [[1.0,item]]
                return [constr+[None,None,'e']]      # list of equivalent variables & mults
            elif 'function' in constType:
                constr = map(list,zip([1.0 for i in range(len(varbs))],varbs))
                return [constr+[None,False,'f']]         #just one constraint
            else:       #'constraint'
                constr = map(list,zip([1.0 for i in range(len(varbs))],varbs))
                return [constr+[1.0,None,'c']]          #just one constraint - default sum to one
        return []

    def CheckAddedConstraint(newcons):
        '''Check a new constraint that has just been input.
        If there is an error display a message and give the user a
        choice to keep or discard the last entry (why keep? -- they
        may want to delete something else or edit multipliers).
        Since the varylist is not available, no warning messages
        should be generated.
        Returns True if constraint should be added
        '''
        allcons = []
        for key in 'Hist','HAP','Phase':
            allcons += data[key]
        allcons += newcons
        if not len(allcons): return True
        G2mv.InitVars()    
        constDictList,fixedList,ignored = G2str.ProcessConstraints(allcons)
        errmsg, warnmsg = G2mv.CheckConstraints('',constDictList,fixedList)
        if errmsg:
            res = G2frame.ErrorDialog('Constraint Error',
                'Error with newly added constraint:\n'+errmsg+
                '\n\nDiscard newly added constraint?',parent=G2frame.dataFrame,
                wtype=wx.YES_NO)
            return res != wx.ID_YES
        elif warnmsg:
            print 'Unexpected contraint warning:\n',warnmsg
        return True

    def CheckChangedConstraint():
        '''Check all constraints after an edit has been made.
        If there is an error display a message and give the user a
        choice to keep or discard the last edit.
        Since the varylist is not available, no warning messages
        should be generated.
        Returns True if the edit should be retained
        '''
        allcons = []
        for key in 'Hist','HAP','Phase':
            allcons += data[key]
        if not len(allcons): return True
        G2mv.InitVars()    
        constDictList,fixedList,ignored = G2str.ProcessConstraints(allcons)
        errmsg, warnmsg = G2mv.CheckConstraints('',constDictList,fixedList)
        if errmsg:
            res = G2frame.ErrorDialog('Constraint Error',
                'Error after editing constraint:\n'+errmsg+
                '\n\nDiscard last constraint edit?',parent=G2frame.dataFrame,
                wtype=wx.YES_NO)
            return res != wx.ID_YES
        elif warnmsg:
            print 'Unexpected contraint warning:\n',warnmsg
        return True
             
    def OnAddHold(event):
        '''add a Hold constraint'''
        for phase in Phases:
            Phase = Phases[phase]
            Atoms = Phase['Atoms']
        constr = []
        page = G2frame.Page
        choice = scope[page[1]]
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in choice[2]]
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],atchoice)
        else:    
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],choice[2])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            FrstVarb = choice[2][sel]
            newcons = [[[0.0,FrstVarb],None,None,'h']]
            if CheckAddedConstraint(newcons):
                data[choice[3]] += newcons
        dlg.Destroy()
        choice[4]()
        
    def OnAddEquivalence(event):
        '''add an Equivalence constraint'''
        constr = []
        page = G2frame.Page
        choice = scope[page[1]]
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in choice[2]]
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],atchoice)
        else:    
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],choice[2])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            FrstVarb = choice[2][sel]
            moreVarb = FindEquivVarb(FrstVarb,choice[2])
            newcons = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'equivalence')
            if len(newcons) > 0:
                if CheckAddedConstraint(newcons):
                    data[choice[3]] += newcons
        dlg.Destroy()
        choice[4]()
   
    def OnAddFunction(event):
        '''add a Function (new variable) constraint'''
        constr = []
        page = G2frame.Page
        choice = scope[page[1]]
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in choice[2]]
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],atchoice)
        else:    
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],choice[2])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            FrstVarb = choice[2][sel]
            moreVarb = FindEquivVarb(FrstVarb,choice[2])
            newcons = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'function')
            if len(newcons) > 0:
                if CheckAddedConstraint(newcons):
                    data[choice[3]] += newcons
        dlg.Destroy()
        choice[4]()
                        
    def OnAddConstraint(event):
        '''add a constraint equation to the constraints list'''
        constr = []
        page = G2frame.Page
        choice = scope[page[1]]
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in choice[2]]
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],atchoice)
        else:    
            dlg = wx.SingleChoiceDialog(G2frame,'Select 1st variable:'+choice[1],choice[0],choice[2])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            FrstVarb = choice[2][sel]
            moreVarb = FindEquivVarb(FrstVarb,choice[2])
            newcons = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'constraint')
            if len(newcons) > 0:
                if CheckAddedConstraint(newcons):
                    data[choice[3]] += newcons
        dlg.Destroy()
        choice[4]()
                        
    def ConstSizer(name,pageDisplay):
        '''This creates a sizer displaying all of the constraints entered
        '''
        constSizer = wx.FlexGridSizer(1,4,0,0)
        maxlen = 70 # characters before wrapping a constraint
        for Id,item in enumerate(data[name]):
            eqString = ['',]
            if item[-1] == 'h':
                constSizer.Add((5,5),0)              # blank space for edit button
                typeString = ' FIXED   '
                eqString[-1] = item[0][1]+'   '
            elif isinstance(item[-1],str):
                constEdit = wx.Button(pageDisplay,-1,'Edit',style=wx.BU_EXACTFIT)
                constEdit.Bind(wx.EVT_BUTTON,OnConstEdit)
                constSizer.Add(constEdit)            # edit button
                Indx[constEdit.GetId()] = [Id,name]
                if item[-1] == 'f':
                    for term in item[:-3]:
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        m = term[0]
                        if eqString[-1] != '':
                            if m >= 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                                m = abs(m)
                        eqString[-1] += '%.3f*%s '%(m,term[1])
                    typeString = ' NEWVAR  '
                    eqString[-1] += ' = New Variable   '
                elif item[-1] == 'c':
                    for term in item[:-3]:
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        if eqString[-1] != '':
                            if term[0] > 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                        eqString[-1] += '%.3f*%s '%(abs(term[0]),term[1])
                    typeString = ' CONSTR  '
                    eqString[-1] += ' = %.3f'%(item[-3])+'  '
                elif item[-1] == 'e':
                    for term in item[:-3]:
                        if term[0] == 0: term[0] = 1.0
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        if eqString[-1] == '':
                            eqString[-1] += '%s '%(term[1])
                            first = term[0]
                        else:
                            eqString[-1] += ' = %.3f*%s '%(first/term[0],term[1])
                    typeString = ' EQUIV   '
                else:
                    print 'Unexpected constraint',item
            else:
                print 'Removing old-style constraints'
                data[name] = []
                return constSizer
            constDel = wx.Button(pageDisplay,-1,'Delete',style=wx.BU_EXACTFIT)
            constDel.Bind(wx.EVT_BUTTON,OnConstDel)
            Indx[constDel.GetId()] = [Id,name]
            constSizer.Add(constDel)             # delete button
            constSizer.Add(wx.StaticText(pageDisplay,-1,typeString),0,wx.ALIGN_CENTER_VERTICAL)
            EqSizer = wx.BoxSizer(wx.VERTICAL)
            for s in eqString:
                EqSizer.Add(wx.StaticText(pageDisplay,-1,s),0,wx.ALIGN_CENTER_VERTICAL)
            constSizer.Add(EqSizer,0,wx.ALIGN_CENTER_VERTICAL)
            # if item[-1] == 'f':
            #     constRef = wx.CheckBox(pageDisplay,-1,label=' Refine?') 
            #     constRef.SetValue(item[-2])
            #     constRef.Bind(wx.EVT_CHECKBOX,OnConstRef)
            #     Indx[constRef.GetId()] = item
            #     constSizer.Add(constRef)
            # else:
            #     constSizer.Add((5,5),0)
        return constSizer
                
    # def OnConstRef(event):
    #     Obj = event.GetEventObject()
    #     Indx[Obj.GetId()][-2] = Obj.GetValue()
        
    def OnConstDel(event):
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        del(data[name][Id])
        OnPageChanged(None)        
        
    def OnConstEdit(event):
        '''Called to edit an individual contraint by the Edit button'''
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        sep = '*'
        if data[name][Id][-1] == 'f':
            items = data[name][Id][:-3]+[[],]
            constType = 'New Variable'
            lbl = 'Enter value for each term in constraint; sum = new variable'
        elif data[name][Id][-1] == 'c':
            items = data[name][Id][:-3]+[
                [data[name][Id][-3],'fixed value ='],[]]
            constType = 'Constraint'
            lbl = 'Edit value for each term in constant constraint sum'
        elif data[name][Id][-1] == 'e':
            items = data[name][Id][:-3]+[[],]
            constType = 'Equivalence'
            lbl = 'The following terms are set to be equal:'
            sep = '/'
        else:
            return
        dlg = G2frame.ConstraintDialog(G2frame.dataFrame,constType,lbl,items,sep)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                prev = data[name][Id]
                result = dlg.GetData()
                if data[name][Id][-1] == 'c':
                    data[name][Id][:-3] = result[:-2]
                    data[name][Id][-3] = result[-2][0]
                else:
                    data[name][Id][:-3] = result[:-1]
                if not CheckChangedConstraint():
                    data[name][Id] = prev
        except:
            import traceback
            print traceback.format_exc()
        finally:
            dlg.Destroy()            
        OnPageChanged(None)                     
    
    def UpdateHAPConstr():
        '''Responds to press on Histogram/Phase Constraints tab,
        shows constraints in data window'''
        HAPConstr.DestroyChildren()
        HAPDisplay = wx.Panel(HAPConstr)
        HAPSizer = wx.BoxSizer(wx.VERTICAL)
        HAPSizer.Add((5,5),0)
        HAPSizer.Add(ConstSizer('HAP',HAPDisplay))
        HAPDisplay.SetSizer(HAPSizer,True)
        Size = HAPSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],250) + 20
        HAPDisplay.SetSize(Size)
        # scroll bar not working, at least not on Mac
        HAPConstr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],250)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def UpdateHistConstr():
        '''Responds to press on Histogram Constraints tab,
        shows constraints in data window'''
        HistConstr.DestroyChildren()
        HistDisplay = wx.Panel(HistConstr)
        HistSizer = wx.BoxSizer(wx.VERTICAL)
        HistSizer.Add((5,5),0)        
        HistSizer.Add(ConstSizer('Hist',HistDisplay))
        HistDisplay.SetSizer(HistSizer,True)
        Size = HistSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],250) + 20
        HistDisplay.SetSize(Size)
        HistConstr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],250)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def UpdatePhaseConstr():
        '''Responds to press on Phase Constraint tab,
        shows constraints in data window'''
        PhaseConstr.DestroyChildren()
        PhaseDisplay = wx.Panel(PhaseConstr)
        PhaseSizer = wx.BoxSizer(wx.VERTICAL)
        PhaseSizer.Add((5,5),0)        
        PhaseSizer.Add(ConstSizer('Phase',PhaseDisplay))
        PhaseDisplay.SetSizer(PhaseSizer,True)
        Size = PhaseSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],250) + 20
        PhaseDisplay.SetSize(Size)
        PhaseConstr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],250)
        G2frame.dataFrame.setSizePosLeft(Size)
    
    def OnPageChanged(event):
        if event:       #page change event!
            page = event.GetSelection()
        else:
            page = G2frame.dataDisplay.GetSelection()
        oldPage = G2frame.dataDisplay.ChangeSelection(page)
        text = G2frame.dataDisplay.GetPageText(page)
        if text == 'Histogram/Phase constraints':
            G2frame.Page = [page,'hap']
            UpdateHAPConstr()
        elif text == 'Histogram constraints':
            G2frame.Page = [page,'hst']
            UpdateHistConstr()
        elif text == 'Phase constraints':
            G2frame.Page = [page,'phs']
            UpdatePhaseConstr()

    def SetStatusLine(text):
        Status.SetStatusText(text)                                      
        
    plegend,hlegend,phlegend = GetPHlegends(Phases,Histograms)
    scope = {'hst':['Histogram contraints:',hlegend,histList,'Hist',UpdateHistConstr],
        'hap':['Histogram * Phase contraints:',phlegend,hapList,'HAP',UpdateHAPConstr],
        'phs':['Phase contraints:',plegend,phaseList,'Phase',UpdatePhaseConstr]}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    SetDataMenuBar(G2frame,G2frame.dataFrame.ConstraintMenu)
    G2frame.dataFrame.SetLabel('Constraints')
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    SetStatusLine('')
    
    SetDataMenuBar(G2frame,G2frame.dataFrame.ConstraintMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddConstraint, id=wxID_CONSTRAINTADD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddFunction, id=wxID_FUNCTADD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddEquivalence, id=wxID_EQUIVADD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddHold, id=wxID_HOLDADD)
    G2frame.dataDisplay = GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    
    PhaseConstr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(PhaseConstr,'Phase constraints')
    HAPConstr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(HAPConstr,'Histogram/Phase constraints')
    HistConstr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(HistConstr,'Histogram constraints')
    UpdatePhaseConstr()

    G2frame.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    # validate all the constrants -- should not see any errors here normally
    allcons = []
    for key in 'Hist','HAP','Phase':
        allcons += data[key]
    if not len(allcons): return
    G2mv.InitVars()    
    constDictList,fixedList,ignored = G2str.ProcessConstraints(allcons)
    errmsg, warnmsg = G2mv.CheckConstraints('',constDictList,fixedList)
    if errmsg:
        G2frame.ErrorDialog('Constraint Error','Error in constraints:\n'+errmsg,
            parent=G2frame.dataFrame)
    elif warnmsg:
        print 'Unexpected contraint warning:\n',warnmsg
        
################################################################################
#####  Restraints
################################################################################           
       
def UpdateRestraints(G2frame,data,Phases,phaseName):
    if not len(Phases):
        print 'There are no phases to form restraints'
        return
    phasedata = Phases[phaseName]
    if phaseName not in data:
        data[phaseName] = {}
    restrData = data[phaseName]
    if 'Bond' not in restrData:
        restrData['Bond'] = {'wtFactor':1.0,'Bonds':[],'Use':True}
    if 'Angle' not in restrData:
        restrData['Angle'] = {'wtFactor':1.0,'Angles':[],'Use':True}
    if 'Plane' not in restrData:
        restrData['Plane'] = {'wtFactor':1.0,'Planes':[],'Use':True}
    if 'Chiral' not in restrData:
        restrData['Chiral'] = {'wtFactor':1.0,'Volumes':[],'Use':True}
    
    def OnSelectPhase(event):
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Phase',Phases.keys())
        try:
            if dlg.ShowModal() == wx.ID_OK:
                phaseName = Phases.keys()[dlg.GetSelection()]
                UpdateRestraints(G2frame,data,Phases,phaseName)
        finally:
            dlg.Destroy()
    
    def OnAddRestraint(event):
        page = G2frame.dataDisplay.GetSelection()
        if 'Bond' in G2frame.dataDisplay.GetPageText(page):
            AddBondRestraint()
        elif 'Angle' in G2frame.dataDisplay.GetPageText(page):
            AddAngleRestraint()
        elif 'Plane' in G2frame.dataDisplay.GetPageText(page):
            AddPlaneRestraint()
        elif 'Chiral' in G2frame.dataDisplay.GetPageText(page):
            AddChiralRestraint()
            
    def AddBondRestraint():
        print 'Bond restraint'

    def AddAngleRestraint():
        print 'Angle restraint'

    def AddPlaneRestraint():
        print 'Plane restraint'

    def AddChiralRestraint():
        print 'Chiral restraint'
        
    def WtBox(wind,restData):
        
        def OnWtFactor(event):
            try:
                value = float(wtfactor.GetValue())
            except ValueError:
                value = 1.0
            restData['wtFactor'] = value
            wtfactor.SetValue('%.2f'%(value))
            
        def OnUseData(event):
            restData['Use'] = Obj.GetValue()

        wtBox = wx.BoxSizer(wx.HORIZONTAL)
        wtBox.Add(wx.StaticText(wind,-1,'Restraint weight factor:'),0,wx.ALIGN_CENTER_VERTICAL)
        wtfactor = wx.TextCtrl(wind,-1,value='%.2f'%(restData['wtFactor']),style=wx.TE_PROCESS_ENTER)
        wtfactor.Bind(wx.EVT_TEXT_ENTER,OnWtFactor)
        wtfactor.Bind(wx.EVT_KILL_FOCUS,OnWtFactor)
        wtBox.Add(wtfactor,0,wx.ALIGN_CENTER_VERTICAL)
        useData = wx.CheckBox(wind,-1,label=' Use?')
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(restData['Use'])        
        wtBox.Add(useData,0,wx.ALIGN_CENTER_VERTICAL)
        return wtBox
        
    def UpdateBondRestr(bondRestData):
        
        def OnColSort(event):
            r,c = event.GetRow(),event.GetCol()
            if r < 0 and c == 0:
                names = G2mth.sortArray(table,0)
                bonds = []
                for name in names:
                    idx = table.index(name)
                    bonds.append(bondList[idx])
                bondRestData['Bonds'] = bonds
                UpdateBondRestr(bondRestData)                
        
        def OnChangeValue(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][4]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new value for bond',val,[0.,5.])
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondList[r][4] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                

        def OnChangeEsd(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][5]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for bond',val,[0.,1.])
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondList[r][5] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                
                                
        def OnDeleteRestraint(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            rows.sort()
            rows.reverse()
            for row in rows:
                bondList.remove(bondList[row])
            UpdateBondRestr(bondRestData)                
            
        BondRestr.DestroyChildren()
        dataDisplay = wx.Panel(BondRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(BondRestr,bondRestData),0,wx.ALIGN_CENTER_VERTICAL)

        bondList = bondRestData['Bonds']
        if len(bondList):
            table = []
            rowLabels = []
            colLabels = ['A+SymOp  B+SymOp','d-calc','d-obs','esd']
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,3',]
            for i,[atoms,ops,indx,dcalc,dobs,esd] in enumerate(bondList):
                table.append([atoms[0]+'+ ('+ops[0]+')  '+atoms[1]+'+ ('+ops[1]+')',dcalc,dobs,esd])
                rowLabels.append(str(i))
            bondTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Bonds = GSGrid(BondRestr)
            Bonds.SetTable(bondTable, True)
            Bonds.AutoSizeColumns(False)
            for r in range(len(bondList)):
                for c in range(2):
                    Bonds.SetReadOnly(r,c,True)
                    Bonds.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Bonds.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnColSort)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=wxID_RESRCHANGEVAL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=wxID_RESTCHANGEESD)
            mainSizer.Add(Bonds,0,)
        else:
            mainSizer.Add(wx.StaticText(BondRestr,-1,'No bond distance restraints for this phase'),0,)

        BondRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] += 5
        Size[1] += 25       #make room for tab
        BondRestr.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def UpdateAngleRestr(angleRestData):
        
        def OnColSort(event):
            r,c = event.GetRow(),event.GetCol()
            if r < 0 and c == 0:
                names = G2mth.sortArray(table,0)
                angles = []
                for name in names:
                    idx = table.index(name)
                    angles.append(angleList[idx])
                angleRestData['Angles'] = angles
                UpdateAngleRestr(angleRestData)                
        
        def OnChangeValue(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][4]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new value for angle',val,[0.,360.])
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleList[r][4] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                

        def OnChangeEsd(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][5]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for angle',val,[0.,5.])
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleList[r][5] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                angleList.remove(angleList[row])
            UpdateAngleRestr(angleRestData)                
            
        AngleRestr.DestroyChildren()
        dataDisplay = wx.Panel(AngleRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(AngleRestr,angleRestData),0,wx.ALIGN_CENTER_VERTICAL)

        angleList = angleRestData['Angles']
        if len(angleList):
            table = []
            rowLabels = []
            colLabels = ['A+SymOp  B+SymOp  C+SymOp','calc','obs','esd']
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            for i,[atoms,ops,indx,dcalc,dobs,esd] in enumerate(angleList):
                table.append([atoms[0]+'+ ('+ops[0]+')  '+atoms[1]+'+ ('+ops[1]+')  '+atoms[2]+ \
                '+ ('+ops[2]+')',dcalc,dobs,esd])
                rowLabels.append(str(i))
            angleTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Angles = GSGrid(AngleRestr)
            Angles.SetTable(angleTable, True)
            Angles.AutoSizeColumns(False)
            for r in range(len(angleList)):
                for c in range(2):
                    Angles.SetReadOnly(r,c,True)
                    Angles.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Angles.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnColSort)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=wxID_RESRCHANGEVAL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=wxID_RESTCHANGEESD)
            mainSizer.Add(Angles,0,)
        else:
            mainSizer.Add(wx.StaticText(AngleRestr,-1,'No bond angle restraints for this phase'),0,)

        AngleRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] += 5
        Size[1] += 25       #make room for tab
        AngleRestr.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
    
    def UpdatePlaneRestr(planeRestData):
        
        items = G2frame.dataFrame.RestraintEdit.GetMenuItems()
        for item in items:
            if item.GetLabel() in ['Change value']:
                item.Enable(False)

        def OnDeleteRestraint(event):
            rows = Planes.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                planeList.remove(planeList[row])
            UpdatePlaneRestr(planeRestData)                
            
        PlaneRestr.DestroyChildren()
        dataDisplay = wx.Panel(PlaneRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(PlaneRestr,planeRestData),0,wx.ALIGN_CENTER_VERTICAL)

        planeList = planeRestData['Planes']
        if len(planeList):
            table = []
            rowLabels = []
            colLabels = ['atom+SymOp','calc','obs','esd']
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            for i,[atoms,ops,indx,dcalc,dobs,esd] in enumerate(planeList):
                atString = ''
                for a,atom in enumerate(atoms):
                    atString += atom+'+ ('+ops[a]+'),'
                    if (a+1)%3 == 0:
                        atString += '\n'
                table.append([atString[:-1],dcalc,dobs,esd])
                rowLabels.append(str(i))
            planeTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Planes = GSGrid(PlaneRestr)
            Planes.SetTable(planeTable, True)
            Planes.AutoSizeColumns(False)
            Planes.AutoSizeRows(False)
            for r in range(len(planeList)):
                for c in range(3):
                    Planes.SetReadOnly(r,c,True)
                    Planes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=wxID_RESTDELETE)
            mainSizer.Add(Planes,0,)
        else:
            mainSizer.Add(wx.StaticText(PlaneRestr,-1,'No plane restraints for this phase'),0,)

        PlaneRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] += 5
        Size[1] += 25       #make room for tab
        PlaneRestr.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
    
    def UpdateChiralRestr(chiralRestData):

        def OnDeleteRestraint(event):
            rows = Volumes.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                volumeList.remove(volumeList[row])
            UpdateChiralRestr(chiralRestData)                
            
        ChiralRestr.DestroyChildren()
        dataDisplay = wx.Panel(ChiralRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(ChiralRestr,chiralRestData),0,wx.ALIGN_CENTER_VERTICAL)

        volumeList = chiralRestData['Volumes']
        if len(volumeList):
            table = []
            rowLabels = []
            colLabels = ['O+SymOp  A+SymOp  B+SymOp  C+SymOp','calc','obs','esd']
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            for i,[atoms,ops,indx,dcalc,dobs,esd] in enumerate(volumeList):
                table.append([atoms[0]+'+ ('+ops[0]+') '+atoms[1]+'+ ('+ops[1]+') '+atoms[2]+ \
                '+ ('+ops[2]+') '+atoms[3]+'+ ('+ops[3]+')',dcalc,dobs,esd])
                rowLabels.append(str(i))
            volumeTable = Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Volumes = GSGrid(ChiralRestr)
            Volumes.SetTable(volumeTable, True)
            Volumes.AutoSizeColumns(False)
            for r in range(len(volumeList)):
                for c in range(2):
                    Volumes.SetReadOnly(r,c,True)
                    Volumes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=wxID_RESTDELETE)
            mainSizer.Add(Volumes,0,)
        else:
            mainSizer.Add(wx.StaticText(ChiralRestr,-1,'No chiral volume restraints for this phase'),0,)

        ChiralRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] += 5
        Size[1] += 25       #make room for tab
        ChiralRestr.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
    
    def OnPageChanged(event):
        page = event.GetSelection()
        text = G2frame.dataDisplay.GetPageText(page)
        if text == 'Bond restraints':
            SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            bondRestData = restrData['Bond']
            UpdateBondRestr(bondRestData)
        elif text == 'Angle restraints':
            SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            angleRestData = restrData['Angle']
            UpdateAngleRestr(angleRestData)
        elif text == 'Plane restraints':
            SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            planeRestData = restrData['Plane']
            UpdatePlaneRestr(planeRestData)
        elif text == 'Chiral restraints':
            SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            chiralRestData = restrData['Chiral']
            UpdateChiralRestr(chiralRestData)
        event.Skip()

    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
        
    SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
    G2frame.dataFrame.SetLabel('restraints for '+phaseName)
    G2frame.dataFrame.RestraintEdit.Enable(wxID_RESTSELPHASE,False)
    if len(Phases) > 1:
        G2frame.dataFrame.RestraintEdit.Enable(wxID_RESTSELPHASE,True)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnSelectPhase, id=wxID_RESTSELPHASE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddRestraint, id=wxID_RESTRAINTADD)
    G2frame.dataDisplay = GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    
    BondRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(BondRestr,'Bond restraints')
    AngleRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(AngleRestr,'Angle restraints')
    PlaneRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(PlaneRestr,'Plane restraints')
    ChiralRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(ChiralRestr,'Chiral restraints')
    UpdateBondRestr(restrData['Bond'])

    G2frame.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    
################################################################################
#####  Main PWDR panel
################################################################################           
       
def UpdatePWHKPlot(G2frame,kind,item):

    def OnErrorAnalysis(event):
        G2plt.PlotDeltSig(G2frame,kind)
        
    def OnWtFactor(event):
        try:
            val = float(wtval.GetValue())
        except ValueError:
            val = data[0]['wtFactor']
        data[0]['wtFactor'] = val
        wtval.SetValue('%.3f'%(val))
           
    data = G2frame.PatternTree.GetItemPyData(item)
    if 'wtFactor' not in data[0]:
        data[0] = {'wtFactor':1.0}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    SetDataMenuBar(G2frame,G2frame.dataFrame.ErrorMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU,OnErrorAnalysis, id=wxID_PWDANALYSIS)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),)
    wtSizer = wx.BoxSizer(wx.HORIZONTAL)
    wtSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Weight factor: '),0,wx.ALIGN_CENTER_VERTICAL)
    wtval = wx.TextCtrl(G2frame.dataDisplay,-1,'%.3f'%(data[0]['wtFactor']),style=wx.TE_PROCESS_ENTER)
    wtval.Bind(wx.EVT_TEXT_ENTER,OnWtFactor)
    wtval.Bind(wx.EVT_KILL_FOCUS,OnWtFactor)
    wtSizer.Add(wtval,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(wtSizer)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
    G2frame.PatternTree.SetItemPyData(item,data)
    if kind == 'PWDR':
        G2plt.PlotPatterns(G2frame,newPlot=True)
    elif kind == 'HKLF':
        G2plt.PlotSngl(G2frame,newPlot=True)
                 
################################################################################
#####  HKLF controls
################################################################################           
       
def UpdateHKLControls(G2frame,data):
    
    def OnScaleSlider(event):
        scale = int(scaleSel.GetValue())/1000.
        scaleSel.SetValue(int(scale*1000.))
        data['Scale'] = scale*10.
        G2plt.PlotSngl(G2frame)
        
    def OnLayerSlider(event):
        layer = layerSel.GetValue()
        data['Layer'] = layer
        G2plt.PlotSngl(G2frame)
        
    def OnSelZone(event):
        data['Zone'] = zoneSel.GetValue()
        izone = zones.index(data['Zone'])
        layerSel.SetRange(maxValue=HKLmax[izone],minValue=HKLmin[izone])
        G2plt.PlotSngl(G2frame,newPlot=True)
        
    def OnSelType(event):
        data['Type'] = typeSel.GetValue()
        G2plt.PlotSngl(G2frame)
        
    def SetStatusLine():
        Status.SetStatusText("")
                                      
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    SetStatusLine()
    zones = ['100','010','001']
    HKLmax = data['HKLmax']
    HKLmin = data['HKLmin']
    if data['ifFc']:
        typeChoices = ['Fosq','Fo','|DFsq|/sig','|DFsq|>sig','|DFsq|>3sig']
    else:
        typeChoices = ['Fosq','Fo']
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    SetDataMenuBar(G2frame)
    G2frame.dataFrame.SetTitle('HKL Plot Controls')
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Scale'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    scaleSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=1000,minValue=1,
        style=wx.SL_HORIZONTAL,value=int(data['Scale']*100))
    scaleSizer.Add(scaleSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    scaleSel.SetLineSize(10)
    scaleSel.SetPageSize(10)
    scaleSel.Bind(wx.EVT_SLIDER, OnScaleSlider)
    mainSizer.Add(scaleSizer,0,wx.EXPAND|wx.RIGHT)
    mainSizer.Add((0,10),0)    
    
    zoneSizer = wx.BoxSizer(wx.HORIZONTAL)
    zoneSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Zone  '),0,
        wx.ALIGN_CENTER_VERTICAL)
    zoneSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['Zone'],choices=['100','010','001'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    zoneSel.Bind(wx.EVT_COMBOBOX, OnSelZone)
    zoneSizer.Add(zoneSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Plot type  '),0,
        wx.ALIGN_CENTER_VERTICAL)        
    typeSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['Type'],choices=typeChoices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    typeSel.Bind(wx.EVT_COMBOBOX, OnSelType)
    zoneSizer.Add(typeSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add((10,0),0)    
    mainSizer.Add(zoneSizer,0,wx.EXPAND|wx.RIGHT)
    mainSizer.Add((0,10),0)    
        
    izone = zones.index(data['Zone'])
    layerSizer = wx.BoxSizer(wx.HORIZONTAL)
    layerSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Layer'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    layerSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=HKLmax[izone],minValue=HKLmin[izone],
        style=wx.SL_HORIZONTAL|wx.SL_AUTOTICKS|wx.SL_LABELS,value=0)
    layerSel.SetLineSize(1)
    layerSel.SetPageSize(1)
    layerSel.Bind(wx.EVT_SLIDER, OnLayerSlider)    
    layerSizer.Add(layerSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    layerSizer.Add((10,0),0)    
    mainSizer.Add(layerSizer,1,wx.EXPAND|wx.RIGHT)

        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))

################################################################################
#####  Pattern tree routines
################################################################################           
       
def GetPatternTreeDataNames(G2frame,dataTypes):
    names = []
    item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)        
    while item:
        name = G2frame.PatternTree.GetItemText(item)
        if name[:4] in dataTypes:
            names.append(name)
        item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    return names
                          
def GetPatternTreeItemId(G2frame, parentId, itemText):
    item, cookie = G2frame.PatternTree.GetFirstChild(parentId)
    while item:
        if G2frame.PatternTree.GetItemText(item) == itemText:
            return item
        item, cookie = G2frame.PatternTree.GetNextChild(parentId, cookie)
    return 0                

def MovePatternTreeToGrid(G2frame,item):
    
#    print G2frame.PatternTree.GetItemText(item)
    
    oldPage = 0
    if G2frame.dataFrame:
        SetDataMenuBar(G2frame)
        if G2frame.dataFrame.GetLabel() == 'Comments':
            data = [G2frame.dataDisplay.GetValue()]
            G2frame.dataDisplay.Clear() 
            Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Comments')
            if Id: G2frame.PatternTree.SetItemPyData(Id,data)
        elif G2frame.dataFrame.GetLabel() == 'Notebook':
            data = [G2frame.dataDisplay.GetValue()]
            G2frame.dataDisplay.Clear() 
            Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Notebook')
            if Id: G2frame.PatternTree.SetItemPyData(Id,data)
        elif 'Phase Data for' in G2frame.dataFrame.GetLabel():
            if G2frame.dataDisplay: 
                oldPage = G2frame.dataDisplay.GetSelection()
        G2frame.dataFrame.Clear()
        G2frame.dataFrame.SetLabel('')
    else:
        #create the frame for the data item window
        G2frame.dataFrame = DataFrame(parent=G2frame.mainPanel,frame=G2frame)

    G2frame.dataFrame.Raise()            
    G2frame.PickId = 0
    parentID = G2frame.root
    for i in G2frame.ExportPattern: i.Enable(False)
    defWid = [250,150]
    if item != G2frame.root:
        parentID = G2frame.PatternTree.GetItemParent(item)
    if G2frame.PatternTree.GetItemParent(item) == G2frame.root:
        G2frame.PatternId = item
        G2frame.PickId = item
        if G2frame.PatternTree.GetItemText(item) == 'Notebook':
            SetDataMenuBar(G2frame,G2frame.dataFrame.DataNotebookMenu)
            G2frame.PatternId = 0
            for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateNotebook(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Controls':
            G2frame.PatternId = 0
            for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            if not data:           #fill in defaults
                data = {
                    #least squares controls
                    'deriv type':'analytic Hessian','min dM/M':0.0001,'shift factor':1.0,'max cyc':3}
                G2frame.PatternTree.SetItemPyData(item,data)                             
            for i in G2frame.Refine: i.Enable(True)
            for i in G2frame.SeqRefine: i.Enable(True)
            UpdateControls(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Sequental results':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateSeqResults(G2frame,data)            
        elif G2frame.PatternTree.GetItemText(item) == 'Covariance':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2frame.dataFrame.setSizePosLeft(defWid)
            text = ''
            if 'Rvals' in data:
                Nvars = len(data['varyList'])
                Rvals = data['Rvals']
                text = '\nFinal residuals: \nRwp = %.3f%% \nchi**2 = %.1f \nGOF = %.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
                text += '\nNobs = %d \nNvals = %d'%(Rvals['Nobs'],Nvars)
                if 'lamMax' in Rvals:
                    text += '\nlog10 MaxLambda = %.1f'%(np.log10(Rvals['lamMax']))
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='See plot window for covariance display'+text,style=wx.TE_MULTILINE)
            G2plt.PlotCovariance(G2frame)
        elif G2frame.PatternTree.GetItemText(item) == 'Constraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateConstraints(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Restraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            Phases = G2frame.GetPhaseData()
            phase = ''
            phaseName = ''
            if Phases:
                phaseName = Phases.keys()[0]
            UpdateRestraints(G2frame,data,Phases,phaseName)
        elif 'IMG' in G2frame.PatternTree.GetItemText(item):
            G2frame.Image = item
            G2plt.PlotImage(G2frame,newPlot=True)
        elif 'PKS' in G2frame.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(G2frame)
        elif 'PWDR' in G2frame.PatternTree.GetItemText(item):
            for i in G2frame.ExportPattern: i.Enable(True)
            UpdatePWHKPlot(G2frame,'PWDR',item)
        elif 'HKLF' in G2frame.PatternTree.GetItemText(item):
            G2frame.Sngl = item
            UpdatePWHKPlot(G2frame,'HKLF',item)
        elif 'PDF' in G2frame.PatternTree.GetItemText(item):
            G2frame.PatternId = item
            for i in G2frame.ExportPDF: i.Enable(True)
            G2plt.PlotISFG(G2frame,type='S(Q)')
        elif G2frame.PatternTree.GetItemText(item) == 'Phases':
            G2frame.dataFrame.setSizePosLeft(defWid)
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='Select one phase to see its parameters')            
    elif 'I(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='I(Q)',newPlot=True)
    elif 'S(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='S(Q)',newPlot=True)
    elif 'F(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='F(Q)',newPlot=True)
    elif 'G(R)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='G(R)',newPlot=True)            
    elif G2frame.PatternTree.GetItemText(parentID) == 'Phases':
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2phG.UpdatePhaseData(G2frame,item,data,oldPage)
    elif G2frame.PatternTree.GetItemText(item) == 'Comments':
        SetDataMenuBar(G2frame,G2frame.dataFrame.DataCommentsMenu)
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateComments(G2frame,data)
    elif G2frame.PatternTree.GetItemText(item) == 'Image Controls':
        G2frame.dataFrame.SetTitle('Image Controls')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        masks = G2frame.PatternTree.GetItemPyData(
            GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateImageControls(G2frame,data,masks)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Masks':
        G2frame.dataFrame.SetTitle('Masks')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateMasks(G2frame,data)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Stress/Strain':
        G2frame.dataFrame.SetTitle('Stress/Strain')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateStressStrain(G2frame,data)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'HKL Plot Controls':
        G2frame.PickId = item
        G2frame.Sngl = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateHKLControls(G2frame,data)
        G2plt.PlotSngl(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'PDF Controls':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPDF: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='I(Q)')
        G2plt.PlotISFG(G2frame,type='S(Q)')
        G2plt.PlotISFG(G2frame,type='F(Q)')
        G2plt.PlotISFG(G2frame,type='G(R)')
    elif G2frame.PatternTree.GetItemText(item) == 'Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Background':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateBackground(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Limits':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateLimitsGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Instrument Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateInstrumentGrid(G2frame,data)
        G2plt.PlotPeakWidths(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Sample Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)

        if 'Temperature' not in data:           #temp fix for old gpx files
            data = {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],'DisplaceX':[0.0,False],
                'DisplaceY':[0.0,False],'Diffuse':[],'Temperature':300.,'Pressure':1.0,'Humidity':0.0,'Voltage':0.0,
                'Force':0.0,'Gonio. radius':200.0}
            G2frame.PatternTree.SetItemPyData(item,data)
    
        G2pdG.UpdateSampleGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Index Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateIndexPeaksGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Unit Cells List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        if not data:
            data.append([0,0.0,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, zero value, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            G2frame.PatternTree.SetItemPyData(item,data)                             
        G2pdG.UpdateUnitCellsGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection Lists':   #powder reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2frame.RefList = ''
        if len(data):
            G2frame.RefList = data.keys()[0]
        G2pdG.UpdateReflectionGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection List':    #HKLF reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        name = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        data = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)
        G2pdG.UpdateReflectionGrid(G2frame,data,HKLF=True,Name=name)

def SetDataMenuBar(G2frame,menu=None):
        '''Set the menu for the data frame. On the Mac put this
        menu for the data tree window instead.

        Note that data frame items do not have menus, for these (menu=None)
        display a blank menu or on the Mac display the standard menu for
        the data tree window.
        '''
        if sys.platform == "darwin":
            if menu is None:
                G2frame.SetMenuBar(G2frame.GSASIIMenu)
            else:
                G2frame.SetMenuBar(menu)
        else:
            if menu is None:
                G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
            else:
                G2frame.dataFrame.SetMenuBar(menu)
