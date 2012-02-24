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
import GSASIIIO as G2IO
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIstruct as G2str
#import GSASIImapvars as G2mv

# globals we will use later
__version__ = None # gets overridden in GSASII.py
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # save location of this file
helpLocDict = {}
htmlPanel = None
htmlFrame = None
if sys.platform.lower().startswith('win'): 
#    helpMode = 'browser'    # need a global control to set this
    helpMode = 'internal'    # need a global control to set this
else:
    helpMode = 'browser'    # need a global control to set this
htmlFirstUse = True

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSTESTADD, wxID_ATONTESTINSERT,
    wxID_RELOADDRAWATOMS,wxID_ATOMSDISAGL,
] = [wx.NewId() for _init_coll_Atom_Items in range(10)]

[ wxID_PWDRADD, wxID_HKLFADD, wxID_DATADELETE,
] = [wx.NewId() for _init_coll_Data_Items in range(3)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD, wxID_DRAWDISAGLTOR,  wxID_DRAWPLANE,
] = [wx.NewId() for _init_coll_DrawAtom_Items in range(12)]

[ wxID_IMCALIBRATE,wxID_IMRECALIBRATE,wxID_IMINTEGRATE, wxID_IMCLEARCALIB,  
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL, wxID_IMSAVECONTROLS, wxID_IMLOADCONTROLS,
] = [wx.NewId() for _init_coll_IMAGE_Items in range(8)]

[ wxID_MASKCOPY, wxID_MASKSAVE, wxID_MASKLOAD,
] = [wx.NewId() for _init_coll_MASK_Items in range(3)]

[ wxID_PAWLEYLOAD, wxID_PAWLEYIMPORT, wxID_PAWLEYDELETE, wxID_PAWLEYESTIMATE,
] = [wx.NewId() for _init_coll_PAWLEY_Items in range(4)]

[ wxID_INSTPRMRESET,wxID_CHANGEWAVETYPE,wxID_INSTCOPY,
] = [wx.NewId() for _init_coll_INST_Items in range(3)]

[ wxID_INDXRELOAD,
] = [wx.NewId() for _init_coll_IndPeaks_Items in range(1)]

[ wxID_UNDO,wxID_LSQPEAKFIT,wxID_LSQONECYCLE,wxID_RESETSIGGAM,wxID_CLEARPEAKS,
] = [wx.NewId() for _init_coll_PEAK_Items in range(5)]

[  wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL, wxID_MAKENEWPHASE,
] = [wx.NewId() for _init_coll_INDEX_Items in range(4)]

[ wxID_BACKCOPY,
] = [wx.NewId() for _init_coll_Back_Items in range(1)]

[ wxID_LIMITCOPY,
] = [wx.NewId() for _init_coll_Limit_Items in range(1)]

[ wxID_SAMPLECOPY,
] = [wx.NewId() for _init_coll_Sample_Items in range(1)]

[ wxID_CONSTRAINTADD,wxID_EQUIVADD,wxID_HOLDADD,wxID_FUNCTADD,
] = [wx.NewId() for _init_coll_Constraint_Items in range(4)]

[ wxID_RESTRAINTADD,
] = [wx.NewId() for _init_coll_Restraint_Items in range(1)]

[ wxID_SAVESEQSEL,
] = [wx.NewId() for _init_coll_Sequential_Items in range(1)]

[ wxID_SELECTPHASE,
] = [wx.NewId() for _init_coll_Refl_Items in range(1)]

[ wxID_FOURCALC,
] = [wx.NewId() for _init_coll_Genrl_Items in range(1)]

[ wxID_CLEARTEXTURE,wxID_REFINETEXTURE,
] = [wx.NewId() for _init_coll_Texture_Items in range(2)]

[ wxID_PDFCOPYCONTROLS, wxID_PDFSAVECONTROLS, wxID_PDFLOADCONTROLS, 
    wxID_PDFCOMPUTE, wxID_PDFCOMPUTEALL, wxID_PDFADDELEMENT, wxID_PDFDELELEMENT,
] = [wx.NewId() for _init_coll_PDF_Items in range(7)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

def ShowHelp(helpType,frame):
    '''Called to bring up a web page for documentation.'''
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
    NOTE: the title for this menu should be '&Help' so the wx handles
    it correctly. BHT
    '''
    def __init__(self,frame,title='',helpType=None):
        wx.Menu.__init__(self,title)
        self.helpType = helpType
        self.frame = frame
        # add a help item only when helpType is specified
        if helpType is not None:
            helpobj = self.Append(text='Help on '+helpType,
                id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
            frame.Bind(wx.EVT_MENU, self.OnHelp, helpobj)
        self.Append(help='', id=wx.ID_ABOUT, kind=wx.ITEM_NORMAL,
            text='&About GSAS-II')
        frame.Bind(wx.EVT_MENU, self.OnHelpAbout, id=wx.ID_ABOUT)

    def OnHelp(self,event):
        '''Called when Help on... is pressed in a menu. Brings up
        a web page for documentation.
        '''
        ShowHelp(self.helpType,self.frame)

    def OnHelpAbout(self, event):
        "Display an 'About GSAS-II' box"
        global __version__
        info = wx.AboutDialogInfo()
        info.Name = 'GSAS-II'
        info.Version = __version__
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

class MyHtmlPanel(wx.Panel):
    '''Defines a panel to display Help information'''
    def __init__(self, frame, id):
        self.frame = frame
        wx.Panel.__init__(self, frame, id)
        sizer = wx.BoxSizer(wx.VERTICAL)
        back = wx.Button(self, -1, "Back")
        back.Bind(wx.EVT_BUTTON, self.OnBack)
        self.htmlwin = G2HtmlWindow(self, id, size=(602,450))
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

    def _init_menus(self):
        
# define all GSAS-II menus        
        
        self.BlankMenu = wx.MenuBar()
        
# Controls
        self.ControlsMenu = wx.MenuBar()
        self.ControlsMenu.Append(menu=MyHelp(self,helpType='Controls'),title='&Help')
        
# Notebook
        self.DataNotebookMenu = wx.MenuBar()
        self.DataNotebookMenu.Append(menu=MyHelp(self,helpType='Notebook'),title='&Help')
        
# Comments
        self.DataCommentsMenu = wx.MenuBar()
        self.DataCommentsMenu.Append(menu=MyHelp(self,helpType='Comments'),title='&Help')
        
# Constraints
        self.ConstraintMenu = wx.MenuBar()
        self.ConstraintEdit = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintEdit, title='Edit')
        self.ConstraintMenu.Append(menu=MyHelp(self,helpType='Constraints'),title='&Help')
        self.ConstraintEdit.Append(id=wxID_HOLDADD, kind=wx.ITEM_NORMAL,text='Add hold',
            help='Add hold on a parameter value')
        self.ConstraintEdit.Append(id=wxID_EQUIVADD, kind=wx.ITEM_NORMAL,text='Add equivalence',
            help='Add equivalence between parameter values')
        self.ConstraintEdit.Append(id=wxID_CONSTRAINTADD, kind=wx.ITEM_NORMAL,text='Add constraint',
            help='Add constraint on parameter values')
        self.ConstraintEdit.Append(id=wxID_FUNCTADD, kind=wx.ITEM_NORMAL,text='Add function',
            help='Add function of parameter values')
            
# Restraints
        self.RestraintMenu = wx.MenuBar()
        self.RestraintEdit = wx.Menu(title='')
        self.RestraintMenu.Append(menu=self.RestraintEdit, title='Edit')
        self.RestraintMenu.Append(menu=MyHelp(self,helpType='Restraints'),title='&Help')
        self.RestraintEdit.Append(id=wxID_RESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add restraint',
            help='restraint dummy menu item')
            
# Sequential results
        self.SequentialMenu = wx.MenuBar()
        self.SequentialFile = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialFile, title='File')
        self.SequentialMenu.Append(menu=MyHelp(self,helpType='Sequential'),title='&Help')
        self.SequentialFile.Append(id=wxID_SAVESEQSEL, kind=wx.ITEM_NORMAL,text='Save...',
            help='Save selected sequential refinement results')
            
# PDR / Limits
        self.LimitMenu = wx.MenuBar()
        self.LimitEdit = wx.Menu(title='')
        self.LimitMenu.Append(menu=self.LimitEdit, title='File')
        self.LimitMenu.Append(menu=MyHelp(self,helpType='Limits'),title='&Help')
        self.LimitEdit.Append(id=wxID_LIMITCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy limits to other histograms')
            
# PDR / Background
        self.BackMenu = wx.MenuBar()
        self.BackEdit = wx.Menu(title='')
        self.BackMenu.Append(menu=self.BackEdit, title='File')
        self.BackMenu.Append(menu=MyHelp(self,helpType='Background'),title='&Help')
        self.BackEdit.Append(id=wxID_BACKCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy background parameters to other histograms')
            
# PDR / Instrument Parameters
        self.InstMenu = wx.MenuBar()
        self.InstEdit = wx.Menu(title='')
        self.InstMenu.Append(menu=self.InstEdit, title='Operations')
        self.InstMenu.Append(menu=MyHelp(self,helpType='Instrument Parameters'),title='&Help')
        self.InstEdit.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')
        self.InstEdit.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')
        self.InstEdit.Append(help='Change radiation type (Ka12 - synch)', 
            id=wxID_CHANGEWAVETYPE, kind=wx.ITEM_NORMAL,text='Change radiation')
        
# PDR / Sample Parameters
        self.SampleMenu = wx.MenuBar()
        self.SampleEdit = wx.Menu(title='')
        self.SampleMenu.Append(menu=self.SampleEdit, title='File')
        self.SampleMenu.Append(menu=MyHelp(self,helpType='Sample Parameters'),title='&Help')
        self.SampleEdit.Append(id=wxID_SAMPLECOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy refinable sample parameters to other histograms')

# PDR / Peak List
        self.PeakMenu = wx.MenuBar()
        self.PeakEdit = wx.Menu(title='')
        self.PeakMenu.Append(menu=self.PeakEdit, title='Peak Fitting')
        self.PeakMenu.Append(menu=MyHelp(self,helpType='Peak List'),title='&Help')
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
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndPeaksMenu.Append(menu=self.IndPeaksEdit,title='Operations')
        self.IndPeaksMenu.Append(menu=MyHelp(self,helpType='Index Peak List'),title='&Help')
        self.IndPeaksEdit.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
        
# PDR / Unit Cells List
        self.IndexMenu = wx.MenuBar()
        self.IndexEdit = wx.Menu(title='')
        self.IndexMenu.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        self.IndexMenu.Append(menu=MyHelp(self,helpType='Unit Cells List'),title='&Help')
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
        self.ReflEdit = wx.Menu(title='')
        self.ReflMenu.Append(menu=self.ReflEdit, title='Reflection List')
        self.ReflMenu.Append(menu=MyHelp(self,helpType='Reflection List'),title='&Help')
        self.SelectPhase = self.ReflEdit.Append(help='Select phase for reflection list',id=wxID_SELECTPHASE, 
            kind=wx.ITEM_NORMAL,text='Select phase')
        
# IMG / Image Controls
        self.ImageMenu = wx.MenuBar()
        self.ImageEdit = wx.Menu(title='')
        self.ImageMenu.Append(menu=self.ImageEdit, title='Operations')
        self.ImageMenu.Append(menu=MyHelp(self,helpType='Image Controls'),title='&Help')
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
        self.MaskEdit = wx.Menu(title='')
        self.MaskMenu.Append(menu=self.MaskEdit, title='Operations')
        self.MaskMenu.Append(menu=MyHelp(self,helpType='Image Masks'),title='&Help')
        self.MaskEdit.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')
        self.MaskEdit.Append(help='Save mask to file', 
            id=wxID_MASKSAVE, kind=wx.ITEM_NORMAL,text='Save mask')
        self.MaskEdit.Append(help='Load mask from file', 
            id=wxID_MASKLOAD, kind=wx.ITEM_NORMAL,text='Load mask')
            
# PDF / PDF Controls
        self.PDFMenu = wx.MenuBar()
        self.PDFEdit = wx.Menu(title='')
        self.PDFMenu.Append(menu=self.PDFEdit, title='PDF Controls')
        self.PDFMenu.Append(menu=MyHelp(self,helpType='PDF Controls'),title='&Help')
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
        self.GeneralCalc = wx.Menu(title='')
        self.DataGeneral.Append(menu=self.GeneralCalc,title='Compute')
        self.DataGeneral.Append(menu=MyHelp(self,helpType='General'),title='&Help')
        self.GeneralCalc.Append(help='Compute Fourier maps',id=wxID_FOURCALC, kind=wx.ITEM_NORMAL,
            text='Fourier maps')
        
# Phase / Atoms tab
        self.AtomsMenu = wx.MenuBar()
        self.AtomEdit = wx.Menu(title='')
        self.AtomCompute = wx.Menu(title='')
        self.AtomsMenu.Append(menu=self.AtomEdit, title='Edit')
        self.AtomsMenu.Append(menu=self.AtomCompute, title='Compute')
        self.AtomsMenu.Append(menu=MyHelp(self,helpType='Atoms'),title='&Help')
        self.AtomEdit.Append(id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append atom',
            help='Inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSTESTADD, kind=wx.ITEM_NORMAL,text='Append test point',
            help='Inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert atom',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATONTESTINSERT, kind=wx.ITEM_NORMAL,text='Insert test point',
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
        self.AtomCompute.Append(id=wxID_ATOMSDISAGL, kind=wx.ITEM_NORMAL,text='Distances & Angles',
            help='Compute distances & angles for selected atoms')   
                 
# Phase / Data tab
        self.DataMenu = wx.MenuBar()
        self.DataEdit = wx.Menu(title='')
        self.DataMenu.Append(menu=self.DataEdit, title='Edit')
        self.DataMenu.Append(menu=MyHelp(self,helpType='Data'),title='&Help')
        self.DataEdit.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Delete histograms',
            help='Delete histograms from use for this phase')
            
# Phase / Texture tab
        self.TextureMenu = wx.MenuBar()
        self.TextureEdit = wx.Menu(title='')
        self.TextureMenu.Append(menu=self.TextureEdit, title='Texture')
        self.TextureMenu.Append(menu=MyHelp(self,helpType='Texture'),title='&Help')
        self.TextureEdit.Append(id=wxID_REFINETEXTURE, kind=wx.ITEM_NORMAL,text='Refine texture', 
            help='Refine the texture coefficients from sequential Pawley results')
        self.TextureEdit.Append(id=wxID_CLEARTEXTURE, kind=wx.ITEM_NORMAL,text='Clear texture', 
            help='Clear the texture coefficients' )
            
# Phase / Pawley tab
        self.PawleyMenu = wx.MenuBar()
        self.PawleyEdit = wx.Menu(title='')
        self.PawleyMenu.Append(menu=self.PawleyEdit,title='Operations')
        self.PawleyMenu.Append(menu=MyHelp(self,helpType='Pawley'),title='&Help')
        self.PawleyEdit.Append(id=wxID_PAWLEYLOAD, kind=wx.ITEM_NORMAL,text='Pawley create',
            help='Initialize Pawley reflection list')
        self.PawleyEdit.Append(id=wxID_PAWLEYESTIMATE, kind=wx.ITEM_NORMAL,text='Pawley estimate',
            help='Estimate initial Pawley intensities')
        self.PawleyEdit.Append(id=wxID_PAWLEYDELETE, kind=wx.ITEM_NORMAL,text='Pawley delete',
            help='Delete Pawley reflection list')

# Phase / Draw Options tab
        self.DataDrawOptions = wx.MenuBar()
        self.DataDrawOptions.Append(menu=MyHelp(self,helpType='Draw Options'),title='&Help')
        
# Phase / Draw Atoms tab
        self.DrawAtomsMenu = wx.MenuBar()
        self.DrawAtomEdit = wx.Menu(title='')
        self.DrawAtomCompute = wx.Menu(title='')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomEdit, title='Edit')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomCompute,title='Compute')
        self.DrawAtomsMenu.Append(menu=MyHelp(self,helpType='Draw Atoms'),title='&Help')
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
        self.DrawAtomCompute.Append(id=wxID_DRAWDISAGLTOR, kind=wx.ITEM_NORMAL,text='Dist. Ang. Tors.',
            help='Compute distance, angle or torsion for 2-4 selected atoms')   
        self.DrawAtomCompute.Append(id=wxID_DRAWPLANE, kind=wx.ITEM_NORMAL,text='Best plane',
            help='Compute best plane for 4+ selected atoms')   
            
# end of GSAS-II menu definitions
        
    def _init_ctrls(self, parent,name=None,size=None,pos=None):
        wx.Frame.__init__(self,parent=parent,style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX,
            size=size,pos=pos,title='GSAS-II data display')
        self._init_menus()
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
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        clientSize = wx.ClientDisplayRect()
        Width[1] = min(Width[1],clientSize[2]-300)
        Width[0] = max(Width[0],300)
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
                
def UpdateNotebook(G2frame,data):        
    if data:
        G2frame.dataFrame.SetLabel('Notebook')
        G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
        for line in data:
            G2frame.dataDisplay.AppendText(line+"\n")
            G2frame.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
            
def UpdateControls(G2frame,data):
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytic Hessian'
        data['min dM/M'] = 0.0001
        data['shift factor'] = 1.
        data['max cyc'] = 3        
    if 'shift factor' not in data:
        data['shift factor'] = 1.
    if 'max cyc' not in data:
        data['max cyc'] = 3        
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
        
        LSSizer = wx.FlexGridSizer(cols=6,vgap=5,hgap=5)
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
        return LSSizer
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText('')
    G2frame.dataFrame.SetLabel('Controls')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.ControlsMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement Controls:'),0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(LSSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(SeqSizer())
    mainSizer.Add((5,5),0)
        
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Density Map Controls:'),0,wx.ALIGN_CENTER_VERTICAL)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
     
def UpdateComments(G2frame,data):                   
    G2frame.dataFrame.SetLabel('Comments')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    for line in data:
        if line[-1] == '\n':
            G2frame.dataDisplay.AppendText(line)
        else:
            G2frame.dataDisplay.AppendText(line+'\n')
            
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
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.SequentialMenu)
    G2frame.dataFrame.SetLabel('Sequental refinement results')
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=wxID_SAVESEQSEL)
    colLabels = data['varyList']+atomList.keys()+cellList.keys()
    Types = len(data['varyList']+atomList.keys()+cellList.keys())*[wg.GRID_VALUE_FLOAT,]
    seqList = [list(data[name]['variables']) for name in histNames]
    
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
    
def UpdateConstraints(G2frame,data):             
#    data.update({'Hist':[],'HAP':[],'Phase':[]})       #empty dict - fill it
    if not data:
        data.update({'Hist':[],'HAP':[],'Phase':[]})       #empty dict - fill it
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    AtomDict = dict([Phases[phase]['pId'],Phases[phase]['Atoms']] for phase in Phases)
    Natoms,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable = G2str.GetPhaseData(Phases,Print=False)
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
        #future -  add 'all:all:name', '0:all:name', etc. to the varList
        if page[1] == 'phs':
            atchoice = [item+' for '+phaseAtNames[item] for item in varList]
            dlg = wx.MultiChoiceDialog(G2frame,'Select more variables:'+legend,FrstVarb+' and:',atchoice)
        else:
            dlg = wx.MultiChoiceDialog(G2frame,'Select more variables:'+legend,FrstVarb+' and:',varList)
        varbs = [FrstVarb,]
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                varbs.append(varList[x])
        dlg.Destroy()
        if len(varbs) > 1:
            if 'equivalence' in constType:
                constr = []
                for item in varbs[1:]:
                    constr += [[[1.0,FrstVarb],[-1.0,item],None,None],]
                return constr           #multiple constraints
            elif 'function' in constType:
                constr = map(list,zip([1.0 for i in range(len(varbs))],varbs))
                return [constr+[0.0,False]]         #just one constraint
            else:       #'constraint'
                constr = map(list,zip([1.0 for i in range(len(varbs))],varbs))
                return [constr+[1.0,None]]          #just one constraint - default sum to one
        return []
             
    def OnAddHold(event):
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
            data[choice[3]] += [[[0.0,FrstVarb],None,None],]
        dlg.Destroy()
        choice[4]()
        
    def OnAddEquivalence(event):
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
            constr = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'equivalence')
            if len(constr) > 0:
                data[choice[3]] += constr
        dlg.Destroy()
        choice[4]()
   
    def OnAddFunction(event):
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
            constr = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'function')
            if len(constr) > 0:
                data[choice[3]] += constr
        dlg.Destroy()
        choice[4]()
                        
    def OnAddConstraint(event):
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
            constr = SelectVarbs(page,FrstVarb,moreVarb,choice[1],'constraint')
            if len(constr) > 0:
                data[choice[3]] += constr
        dlg.Destroy()
        choice[4]()
                        
    def ConstSizer(name,pageDisplay):
        constSizer = wx.FlexGridSizer(1,4,0,0)
        for Id,item in enumerate(data[name]):
            constDel = wx.Button(pageDisplay,-1,'Delete',style=wx.BU_EXACTFIT)
            constDel.Bind(wx.EVT_BUTTON,OnConstDel)
            Indx[constDel.GetId()] = [Id,name]
            if len(item) < 4:
                constSizer.Add((5,5),0)
                constSizer.Add(constDel)
                eqString = ' FIXED   '+item[0][1]+'   '
            else:
                constEdit = wx.Button(pageDisplay,-1,'Edit',style=wx.BU_EXACTFIT)
                constEdit.Bind(wx.EVT_BUTTON,OnConstEdit)
                Indx[constEdit.GetId()] = [Id,name]
                constSizer.Add(constEdit)            
                constSizer.Add(constDel)
                if isinstance(item[-1],bool):
                    eqString = ' FUNCT   '
                elif isinstance(item[-2],float):
                    eqString = ' CONSTR  '
                else:
                    eqString = ' EQUIV   '
                for term in item[:-2]:
                    eqString += '%+.3f*%s '%(term[0],term[1])
                if isinstance(item[-2],float):
                    eqString += ' = %.3f'%(item[-2])+'  '
                else:
                    eqString += ' = 0   '
            constSizer.Add(wx.StaticText(pageDisplay,-1,eqString),0,wx.ALIGN_CENTER_VERTICAL)
            if isinstance(item[-1],bool):
                constRef = wx.CheckBox(pageDisplay,-1,label=' Refine?')                    
                constRef.SetValue(item[-1])
                constRef.Bind(wx.EVT_CHECKBOX,OnConstRef)
                Indx[constRef.GetId()] = item
                constSizer.Add(constRef,0,wx.ALIGN_CENTER_VERTICAL)
            else:
                constSizer.Add((5,5),0)
        return constSizer
                
    def OnConstRef(event):
        Obj = event.GetEventObject()
        Indx[Obj.GetId()][-1] = Obj.GetValue()
        
    def OnConstDel(event):
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        del(data[name][Id])
        OnPageChanged(None)        
        
    def OnConstEdit(event):
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        const = data[name][Id][-2]        
        if isinstance(data[name][Id][-1],bool):
            items = data[name][Id][:-2]+[[],]
            constType = 'Function'
            extra = '; sum = new variable'
        elif isinstance(data[name][Id][-2],float):
            items = data[name][Id][:-2]+[[const,'= fixed value'],[]]
            constType = 'Constraint'
            extra = ' sum = constant'
        else:
            items = data[name][Id][:-2]+[[],]
            constType = 'Equivalence'
            extra = '; sum = 0'
        dlg = G2frame.SumDialog(G2frame,constType,'Enter value for each term in constraint'+extra,'',items)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                if isinstance(data[name][Id][-2],float):
                    data[name][Id][:-2] = result[:-2]
                    data[name][Id][-2] = result[-2][0]
                else:
                    data[name][Id][:-2] = result[:-1]
        finally:
            dlg.Destroy()            
        OnPageChanged(None)                     
    
    def UpdateHAPConstr():
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
        HAPConstr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],250)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def UpdateHistConstr():
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
    scope = {'hst':['Histogram variables:',hlegend,histList,'Hist',UpdateHistConstr],
        'hap':['HAP variables:',phlegend,hapList,'HAP',UpdateHAPConstr],
        'phs':['Phase variables:',plegend,phaseList,'Phase',UpdatePhaseConstr]}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.ConstraintMenu)
    G2frame.dataFrame.SetLabel('Constraints')
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    SetStatusLine('')
    
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.ConstraintMenu)
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
    
    
def UpdateRestraints(G2frame,data):

    def OnAddRestraint(event):
        page = G2frame.dataDisplay.GetSelection()
        print G2frame.dataDisplay.GetPageText(page)

    def UpdateAtomRestr():
        AtomRestr.DestroyChildren()
        dataDisplay = wx.Panel(AtomRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Atom restraint data:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)


        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def UpdatePhaseRestr():
        PhaseRestr.DestroyChildren()
        dataDisplay = wx.Panel(PhaseRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Phase restraint data:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)


        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)
    
    def OnPageChanged(event):
        page = event.GetSelection()
        text = G2frame.dataDisplay.GetPageText(page)
        if text == 'Atom restraints':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.RestraintMenu)
            UpdateAtomRestr()
        elif text == 'Phase restraints':
            UpdatePhaseRestr()
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.RestraintMenu)
        event.Skip()

    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.RestraintMenu)
    G2frame.dataFrame.SetLabel('restraints')
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddRestraint, id=wxID_RESTRAINTADD)
    G2frame.dataDisplay = GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    
    PhaseRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(PhaseRestr,'Phase restraints')
    AtomRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(AtomRestr,'Atom restraints')
    UpdatePhaseRestr()
#    AtomRestrData = data['AtomRestr']

    G2frame.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)        
             
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
        G2plt.PlotSngl(G2frame,newPlot=True)
        
    def OnSelType(event):
        data['Type'] = typeSel.GetValue()
        G2plt.PlotSngl(G2frame)
        
    def SetStatusLine():
        Status.SetStatusText("look at me!!!")
                                      
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
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Scale'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    scaleSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=1000,minValue=100,
        style=wx.SL_HORIZONTAL,value=int(data['Scale']*100))
    scaleSizer.Add(scaleSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    scaleSel.SetLineSize(100)
    scaleSel.SetPageSize(900)
    scaleSel.Bind(wx.EVT_SLIDER, OnScaleSlider)
    mainSizer.Add(scaleSizer,1,wx.EXPAND|wx.RIGHT)
    
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
    mainSizer.Add(zoneSizer,1,wx.EXPAND|wx.RIGHT)
        
    izone = zones.index(data['Zone'])
    layerSizer = wx.BoxSizer(wx.HORIZONTAL)
    layerSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Layer'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    layerSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=HKLmax[izone],minValue=HKLmin[izone],
        style=wx.SL_HORIZONTAL|wx.SL_AUTOTICKS|wx.SL_LABELS,value=0)
    layerSel.SetLineSize(1)
    layerSel.SetLineSize(5)
    layerSel.Bind(wx.EVT_SLIDER, OnLayerSlider)    
    layerSizer.Add(layerSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    layerSizer.Add((10,0),0)    
    mainSizer.Add(layerSizer,1,wx.EXPAND|wx.RIGHT)

        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))

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
        G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
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
        G2frame.dataFrame = DataFrame(parent=G2frame.mainPanel)

    G2frame.dataFrame.Raise()            
    G2frame.PickId = 0
    parentID = G2frame.root
    G2frame.ExportPattern.Enable(False)
    defWid = [250,150]
    if item != G2frame.root:
        parentID = G2frame.PatternTree.GetItemParent(item)
    if G2frame.PatternTree.GetItemParent(item) == G2frame.root:
        G2frame.PatternId = item
        G2frame.PickId = item
        if G2frame.PatternTree.GetItemText(item) == 'Notebook':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataNotebookMenu)
            G2frame.PatternId = 0
            G2frame.ExportPattern.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateNotebook(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Controls':
            G2frame.PatternId = 0
            G2frame.ExportPattern.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            if not data:           #fill in defaults
                data = {
                    #least squares controls
                    'deriv type':'analytic Jacobian','min dM/M':0.0001,'shift factor':1.0,'max cyc':3,
                    #Fourier controls
                    'mapType':'Fobs','d-max':100.,'d-min':0.2,'histograms':[],
                    'stepSize':[0.5,0.5,0.5],'minX':[0.,0.,0.],'maxX':[1.0,1.0,1.0],
                    #distance/angle controls
                    'distMax':0.0,'angleMax':0.0,'useMapPeaks':False}
                G2frame.PatternTree.SetItemPyData(item,data)                             
            G2frame.Refine.Enable(True)
            G2frame.SeqRefine.Enable(True)
            UpdateControls(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Sequental results':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateSeqResults(G2frame,data)            
        elif G2frame.PatternTree.GetItemText(item) == 'Covariance':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2frame.dataFrame.setSizePosLeft(defWid)
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                        value='See plot window for covariance display')
            G2plt.PlotCovariance(G2frame)
        elif G2frame.PatternTree.GetItemText(item) == 'Constraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateConstraints(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Restraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateRestraints(G2frame,data)
        elif 'IMG' in G2frame.PatternTree.GetItemText(item):
            G2frame.Image = item
            G2plt.PlotImage(G2frame,newPlot=True)
        elif 'PKS' in G2frame.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(G2frame)
        elif 'PWDR' in G2frame.PatternTree.GetItemText(item):
            G2frame.ExportPattern.Enable(True)
            G2frame.dataFrame.setSizePosLeft(defWid)
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                style=wx.TE_MULTILINE,
                value='See plot window for powder data display\nor select a data item in histogram')
            G2plt.PlotPatterns(G2frame,newPlot=True)
        elif 'HKLF' in G2frame.PatternTree.GetItemText(item):
            G2frame.Sngl = item
            G2plt.PlotSngl(G2frame,newPlot=True)
        elif 'PDF' in G2frame.PatternTree.GetItemText(item):
            G2frame.PatternId = item
            G2frame.ExportPDF.Enable(True)
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
        G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataCommentsMenu)
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
    elif G2frame.PatternTree.GetItemText(item) == 'HKL Plot Controls':
        G2frame.PickId = item
        G2frame.Sngl = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateHKLControls(G2frame,data)
        G2plt.PlotSngl(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'PDF Controls':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.ExportPDF.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='I(Q)')
        G2plt.PlotISFG(G2frame,type='S(Q)')
        G2plt.PlotISFG(G2frame,type='F(Q)')
        G2plt.PlotISFG(G2frame,type='G(R)')
    elif G2frame.PatternTree.GetItemText(item) == 'Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.ExportPeakList.Enable(True)
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
        G2frame.ExportPeakList.Enable(True)
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
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection Lists':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2frame.RefList = ''
        if len(data):
            G2frame.RefList = data.keys()[0]
        G2pdG.UpdateReflectionGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
