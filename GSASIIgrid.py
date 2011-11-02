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
import numpy as np
import GSASIIpath
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIstruct as G2str
import GSASIImapvars as G2mv

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSTESTADD, wxID_ATONTESTINSERT,
    wxID_RELOADDRAWATOMS,
] = [wx.NewId() for _init_coll_Atom_Items in range(9)]

[ wxID_PWDRADD, wxID_HKLFADD, wxID_DATADELETE,
] = [wx.NewId() for _init_coll_Data_Items in range(3)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD,
] = [wx.NewId() for _init_coll_DrawAtom_Items in range(10)]

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

[ wxID_SAMPLECOPY,
] = [wx.NewId() for _init_coll_Sample_Items in range(1)]

[ wxID_CONSTRAINTADD,
] = [wx.NewId() for _init_coll_Constraint_Items in range(1)]

[ wxID_RESTRAINTADD,
] = [wx.NewId() for _init_coll_Restraint_Items in range(1)]

[ wxID_SELECTPHASE,
] = [wx.NewId() for _init_coll_Refl_Items in range(1)]

[ wxID_CLEARTEXTURE,
] = [wx.NewId() for _init_coll_Texture_Items in range(1)]

[ wxID_PDFCOPYCONTROLS, wxID_PDFSAVECONTROLS, wxID_PDFLOADCONTROLS, 
    wxID_PDFCOMPUTE, wxID_PDFCOMPUTEALL, wxID_PDFADDELEMENT, wxID_PDFDELELEMENT,
] = [wx.NewId() for _init_coll_PDF_Items in range(7)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

class DataFrame(wx.Frame):
    def _init_coll_BlankMenu(self,parent):
        parent.Append(menu=self.Blank,title='')
        
    def _init_coll_AtomsMenu(self,parent):
        parent.Append(menu=self.AtomEdit, title='Edit')
        
    def _init_coll_ConstraintMenu(self,parent):
        parent.Append(menu=self.ConstraintEdit, title='Edit')
        
    def _init_coll_RestraintMenu(self,parent):
        parent.Append(menu=self.RestraintEdit, title='Edit')
        
    def _init_coll_DataMenu(self,parent):
        parent.Append(menu=self.DataEdit, title='Edit')
        
    def _init_coll_DrawAtomsMenu(self,parent):
        parent.Append(menu=self.DrawAtomEdit, title='Edit')

    def _init_coll_PawleyMenu(self,parent):
        parent.Append(menu=self.PawleyEdit,title='Operations')
      
    def _init_coll_IndPeaksMenu(self,parent):
        parent.Append(menu=self.IndPeaksEdit,title='Operations')
                   
    def _init_coll_ImageMenu(self,parent):
        parent.Append(menu=self.ImageEdit, title='Operations')
        
    def _init_coll_BackMenu(self,parent):
        parent.Append(menu=self.BackEdit, title='File')
        
    def _init_coll_InstMenu(self,parent):
        parent.Append(menu=self.InstEdit, title='Operations')
        
    def _init_coll_MaskMenu(self,parent):
        parent.Append(menu=self.MaskEdit, title='Operations')
        
    def _init_coll_SampleMenu(self,parent):
        parent.Append(menu=self.SampleEdit, title='File')
        
    def _init_coll_PeakMenu(self,parent):
        parent.Append(menu=self.PeakEdit, title='Peak Fitting')

    def _init_coll_IndexMenu(self,parent):
        parent.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        
    def _init_coll_ReflMenu(self,parent):
        parent.Append(menu=self.ReflEdit, title='Reflection List')

    def _init_coll_TextureMenu(self,parent):
        parent.Append(menu=self.TextureEdit, title='Texture')

    def _init_coll_PDFMenu(self,parent):
        parent.Append(menu=self.PDFEdit, title='PDF Controls')

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
        parent.Append(id=wxID_RELOADDRAWATOMS, kind=wx.ITEM_NORMAL,text='Reload draw atoms',
            help='Reload atom drawing list')
            
    def _init_coll_Constraint_Items(self,parent):
        parent.Append(id=wxID_CONSTRAINTADD, kind=wx.ITEM_NORMAL,text='Add constraint',
            help='constraint dummy menu item')
        
    def _init_coll_Restraint_Items(self,parent):
        parent.Append(id=wxID_RESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add restraint',
            help='restraint dummy menu item')
                    
    def _init_coll_Sample_Items(self,parent):
        parent.Append(id=wxID_SAMPLECOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy refinable sample parameters to other histograms')
                    
    def _init_coll_Back_Items(self,parent):
        parent.Append(id=wxID_BACKCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy background parameters to other histograms')
                    
    def _init_coll_Data_Items(self,parent):
        parent.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        parent.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        parent.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Delete histograms',
            help='Delete histograms from use for this phase')

    def _init_coll_Texture_Items(self,parent):
        self.ClearPeaks = parent.Append(id=wxID_CLEARTEXTURE, kind=wx.ITEM_NORMAL,text='Clear texture', 
            help='Clear the texture coefficients' )
            
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
#        parent.Append(id=wxID_PAWLEYIMPORT, kind=wx.ITEM_NORMAL,text='Pawley import',
#            help='Import Pawley reflection list')
        parent.Append(id=wxID_PAWLEYESTIMATE, kind=wx.ITEM_NORMAL,text='Pawley estimate',
            help='Estimate initial Pawley intensities')
        parent.Append(id=wxID_PAWLEYDELETE, kind=wx.ITEM_NORMAL,text='Pawley delete',
            help='Delete Pawley reflection list')

    def _init_coll_IndPeaks_Items(self,parent):
        parent.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
            
    def _init_coll_Refl_Items(self,parent):
        self.SelectPhase = parent.Append(help='Select phase for reflection list',id=wxID_SELECTPHASE, 
            kind=wx.ITEM_NORMAL,text='Select phase')
            
    def _init_coll_Image_Items(self,parent):
        parent.Append(help='Calibrate detector by fitting to calibrant lines', 
            id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        parent.Append(help='Recalibrate detector by fitting to calibrant lines', 
            id=wxID_IMRECALIBRATE, kind=wx.ITEM_NORMAL,text='Realibrate')
        parent.Append(help='Clear calibration data points and rings',id=wxID_IMCLEARCALIB, 
            kind=wx.ITEM_NORMAL,text='Clear calibration')
        parent.Append(help='Integrate selected image',id=wxID_IMINTEGRATE, 
            kind=wx.ITEM_NORMAL,text='Integrate')
        parent.Append(help='Integrate all images selected from list',id=wxID_INTEGRATEALL,
            kind=wx.ITEM_NORMAL,text='Integrate all')
        parent.Append(help='Copy image controls to other images', 
            id=wxID_IMCOPYCONTROLS, kind=wx.ITEM_NORMAL,text='Copy Controls')
        parent.Append(help='Save image controls to file', 
            id=wxID_IMSAVECONTROLS, kind=wx.ITEM_NORMAL,text='Save Controls')
        parent.Append(help='Load image controls from file', 
            id=wxID_IMLOADCONTROLS, kind=wx.ITEM_NORMAL,text='Load Controls')

                    
    def _init_coll_Mask_Items(self,parent):
        parent.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')
        parent.Append(help='Save mask to file', 
            id=wxID_MASKSAVE, kind=wx.ITEM_NORMAL,text='Save mask')
        parent.Append(help='Load mask from file', 
            id=wxID_MASKLOAD, kind=wx.ITEM_NORMAL,text='Load mask')

    def _init_coll_Inst_Items(self,parent):
        parent.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')
        parent.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')        
        parent.Append(help='Change radiation type (Ka12 - synch)', 
            id=wxID_CHANGEWAVETYPE, kind=wx.ITEM_NORMAL,text='Change radiation')

    def _init_coll_Peak_Items(self,parent):
        self.UnDo = parent.Append(help='Undo last least squares refinement', 
            id=wxID_UNDO, kind=wx.ITEM_NORMAL,text='UnDo')
        self.PeakFit = parent.Append(id=wxID_LSQPEAKFIT, kind=wx.ITEM_NORMAL,text='LSQ PeakFit', 
            help='Peak fitting via least-squares' )
        self.PFOneCycle = parent.Append(id=wxID_LSQONECYCLE, kind=wx.ITEM_NORMAL,text='LSQ one cycle', 
            help='One cycle of Peak fitting via least-squares' )
        self.ResetSigGam = parent.Append(id=wxID_RESETSIGGAM, kind=wx.ITEM_NORMAL, 
            text='Reset sig and gam',help='Reset sigma and gamma to global fit' )
        self.ClearPeaks = parent.Append(id=wxID_CLEARPEAKS, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the peak list' )
            
    def _init_coll_Index_Items(self,parent):
        self.IndexPeaks = parent.Append(help='', id=wxID_INDEXPEAKS, kind=wx.ITEM_NORMAL,
            text='Index Cell')
        self.CopyCell = parent.Append( id=wxID_COPYCELL, kind=wx.ITEM_NORMAL,text='Copy Cell', 
            help='Copy selected unit cell from indexing to cell refinement fields')
        self.RefineCell = parent.Append( id=wxID_REFINECELL, kind=wx.ITEM_NORMAL, 
            text='Refine Cell',help='Refine unit cell parameters from indexed peaks')
        self.MakeNewPhase = parent.Append( id=wxID_MAKENEWPHASE, kind=wx.ITEM_NORMAL,
            text='Make new phase',help='Make new phase from selected unit cell')
            
    def _init_coll_PDF_Items(self,parent):
        parent.Append(help='Add element to sample composition',id=wxID_PDFADDELEMENT, kind=wx.ITEM_NORMAL,
            text='Add element')
        parent.Append(help='Delete element from sample composition',id=wxID_PDFDELELEMENT, kind=wx.ITEM_NORMAL,
            text='Delete element')
        parent.Append(help='Copy PDF controls', id=wxID_PDFCOPYCONTROLS, kind=wx.ITEM_NORMAL,
            text='Copy controls')
#        parent.Append(help='Load PDF controls from file',id=wxID_PDFLOADCONTROLS, kind=wx.ITEM_NORMAL,
#            text='Load Controls')
#        parent.Append(help='Save PDF controls to file', id=wxID_PDFSAVECONTROLS, kind=wx.ITEM_NORMAL,
#            text='Save controls')
        self.PDFCompute = parent.Append(help='Compute PDF', id=wxID_PDFCOMPUTE, kind=wx.ITEM_NORMAL,
            text='Compute PDF')
        self.PDFCompute = parent.Append(help='Compute all PDFs', id=wxID_PDFCOMPUTEALL, kind=wx.ITEM_NORMAL,
            text='Compute all PDFs')
        

    def _init_utils(self):
        self.BlankMenu = wx.MenuBar()
        
        self.AtomsMenu = wx.MenuBar()
        self.ConstraintMenu = wx.MenuBar()
        self.RestraintMenu = wx.MenuBar()
        self.DataMenu = wx.MenuBar()
        self.TextureMenu = wx.MenuBar()
        self.DrawAtomsMenu = wx.MenuBar()
        self.PawleyMenu = wx.MenuBar()
        self.ImageMenu = wx.MenuBar()
        self.MaskMenu = wx.MenuBar()
        self.BackMenu = wx.MenuBar()
        self.InstMenu = wx.MenuBar()
        self.SampleMenu = wx.MenuBar()
        self.PeakMenu = wx.MenuBar()
        self.IndPeaksMenu = wx.MenuBar()
        self.IndexMenu = wx.MenuBar()
        self.ReflMenu = wx.MenuBar()
        self.PDFMenu = wx.MenuBar()
        self.AtomEdit = wx.Menu(title='')
        self.ConstraintEdit = wx.Menu(title='')
        self.RestraintEdit = wx.Menu(title='')
        self.DataEdit = wx.Menu(title='')
        self.TextureEdit = wx.Menu(title='')
        self.DrawAtomEdit = wx.Menu(title='')
        self.PawleyEdit = wx.Menu(title='')
        self.ImageEdit = wx.Menu(title='')
        self.MaskEdit = wx.Menu(title='')
        self.BackEdit = wx.Menu(title='')
        self.InstEdit = wx.Menu(title='')
        self.SampleEdit = wx.Menu(title='')
        self.PeakEdit = wx.Menu(title='')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndexEdit = wx.Menu(title='')
        self.ReflEdit = wx.Menu(title='')
        self.PDFEdit = wx.Menu(title='')
        self._init_coll_AtomsMenu(self.AtomsMenu)
        self._init_coll_Atom_Items(self.AtomEdit)
        self._init_coll_ConstraintMenu(self.ConstraintMenu)
        self._init_coll_Constraint_Items(self.ConstraintEdit)
        self._init_coll_RestraintMenu(self.RestraintMenu)
        self._init_coll_Restraint_Items(self.RestraintEdit)
        self._init_coll_DataMenu(self.DataMenu)
        self._init_coll_Data_Items(self.DataEdit)
        self._init_coll_TextureMenu(self.TextureMenu)
        self._init_coll_Texture_Items(self.TextureEdit)
        self._init_coll_DrawAtomsMenu(self.DrawAtomsMenu)
        self._init_coll_DrawAtom_Items(self.DrawAtomEdit)
        self._init_coll_PawleyMenu(self.PawleyMenu)
        self._init_coll_Pawley_Items(self.PawleyEdit)
        self._init_coll_ImageMenu(self.ImageMenu)
        self._init_coll_Image_Items(self.ImageEdit)
        self._init_coll_MaskMenu(self.MaskMenu)
        self._init_coll_Mask_Items(self.MaskEdit)
        self._init_coll_BackMenu(self.BackMenu)
        self._init_coll_Back_Items(self.BackEdit)
        self._init_coll_InstMenu(self.InstMenu)
        self._init_coll_Inst_Items(self.InstEdit)
        self._init_coll_SampleMenu(self.SampleMenu)
        self._init_coll_Sample_Items(self.SampleEdit)
        self._init_coll_PeakMenu(self.PeakMenu)
        self._init_coll_Peak_Items(self.PeakEdit)
        self._init_coll_IndPeaksMenu(self.IndPeaksMenu)
        self._init_coll_IndPeaks_Items(self.IndPeaksEdit)
        self._init_coll_IndexMenu(self.IndexMenu)
        self._init_coll_Index_Items(self.IndexEdit)
        self._init_coll_ReflMenu(self.ReflMenu)
        self._init_coll_Refl_Items(self.ReflEdit)
        self._init_coll_PDFMenu(self.PDFMenu)
        self._init_coll_PDF_Items(self.PDFEdit)
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.PFOneCycle.Enable(False)
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
                
def UpdateNotebook(self,data):        
    if data:
        self.dataFrame.SetLabel('Notebook')
        self.dataDisplay = wx.TextCtrl(parent=self.dataFrame,size=self.dataFrame.GetClientSize(),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
        for line in data:
            self.dataDisplay.AppendText(line+"\n")
            self.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
            
def UpdateControls(self,data):
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytical'
        data['min dM/M'] = 0.0001
        data['shift factor'] = 1.
    if 'shift factor' not in data:
        data['shift factor'] = 1.        
    #end patch
    '''
    #Fourier controls
    'mapType':'Fobs','d-max':100.,'d-min':0.2,'histograms':[],
    'stepSize':[0.5,0.5,0.5],'minX':[0.,0.,0.],'maxX':[1.0,1.0,1.0],
    #distance/angle controls
    'distMax':0.0,'angleMax':0.0,'useMapPeaks':False}
    '''
        
    def SetStatusLine(text):
        Status.SetStatusText(text)                                      
        
    def OnConvergence(event):
        try:
            value = max(1.e-9,min(1.0,float(Cnvrg.GetValue())))
        except ValueError:
            value = 0.0001
        data['min dM/M'] = value
        Cnvrg.SetValue('%.2g'%(value))
        
    def OnDerivType(event):
        data['deriv type'] = derivSel.GetValue()
        derivSel.SetValue(data['deriv type'])
        
    def OnFactor(event):
        try:
            value = min(max(float(Factr.GetValue()),0.00001),100.)
        except ValueError:
            value = 1.0
        data['shift factor'] = value
        Factr.SetValue('%.5f'%(value))
        
    def OnSelectData(event):
        choices = ['All',]+GetPatternTreeDataNames(self,['PWDR',])
        sel = []
        if 'Seq Data' in data:
            for item in data['Seq Data']:
                sel.append(choices.index(item))
        names = []
        dlg = wx.MultiChoiceDialog(self,'Select data:','Sequential refinement',choices)
        dlg.SetSelections(sel)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for i in sel: names.append(choices[i])
            if 'All' in names:
                names = choices[1:]
                
        dlg.Destroy()
        data['Seq Data'] = names
        reverseSel.Enable(True)
        
    def OnReverse(event):
        data['Reverse Seq'] = reverseSel.GetValue()
        
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
    LSSizer = wx.FlexGridSizer(cols=6,vgap=5,hgap=5)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Refinement derivatives: '),0,wx.ALIGN_CENTER_VERTICAL)
    Choice=['analytic','numeric']
    derivSel = wx.ComboBox(parent=self.dataDisplay,value=data['deriv type'],choices=Choice,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    derivSel.SetValue(data['deriv type'])
    derivSel.Bind(wx.EVT_COMBOBOX, OnDerivType)    
    LSSizer.Add(derivSel,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Min delta-M/M: '),0,wx.ALIGN_CENTER_VERTICAL)
    Cnvrg = wx.TextCtrl(self.dataDisplay,-1,value='%.2g'%(data['min dM/M']),style=wx.TE_PROCESS_ENTER)
    Cnvrg.Bind(wx.EVT_TEXT_ENTER,OnConvergence)
    Cnvrg.Bind(wx.EVT_KILL_FOCUS,OnConvergence)
    LSSizer.Add(Cnvrg,0,wx.ALIGN_CENTER_VERTICAL)
    LSSizer.Add(wx.StaticText(self.dataDisplay,label=' Initial shift factor: '),0,wx.ALIGN_CENTER_VERTICAL)
    Factr = wx.TextCtrl(self.dataDisplay,-1,value='%.5f'%(data['shift factor']),style=wx.TE_PROCESS_ENTER)
    Factr.Bind(wx.EVT_TEXT_ENTER,OnFactor)
    Factr.Bind(wx.EVT_KILL_FOCUS,OnFactor)
    LSSizer.Add(Factr,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(LSSizer)
    mainSizer.Add((5,5),0)
    
    SeqSizer = wx.BoxSizer(wx.HORIZONTAL)
    SeqSizer.Add(wx.StaticText(self.dataDisplay,label=' Sequential Refinement Powder Data: '),0,wx.ALIGN_CENTER_VERTICAL)
    selSeqData = wx.Button(self.dataDisplay,-1,label=' Select data')
    selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
    SeqSizer.Add(selSeqData,0,wx.ALIGN_CENTER_VERTICAL)
    SeqSizer.Add((5,0),0)
    reverseSel = wx.CheckBox(self.dataDisplay,-1,label=' Reverse order?')
    reverseSel.Bind(wx.EVT_CHECKBOX,OnReverse)
    if 'Seq Data' not in data:
        reverseSel.Enable(False)
    if 'Reverse Seq' in data:
        reverseSel.SetValue(data['Reverse Seq'])
    SeqSizer.Add(reverseSel,0,wx.ALIGN_CENTER_VERTICAL)
    
    
    mainSizer.Add(SeqSizer)
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
        if line[-1] == '\n':
            self.dataDisplay.AppendText(line)
        else:
            self.dataDisplay.AppendText(line+'\n')
            
def UpdateSeqResults(self,data):
    if not data:
        print 'No sequential refinement results'
        return
    histNames = data['histNames']
    
    def GetSigData(parm):
        sigData = []
        for name in histNames:
            sigList = data[name]['sig']
            sigData.append(sigList[parm])
        return sigData
    
    def Select(event):
        cols = self.dataDisplay.GetSelectedCols()
        rows = self.dataDisplay.GetSelectedRows()
        if cols:
            plotData = []
            plotSig = []
            plotNames = []
            for col in cols:
                plotData.append(self.SeqTable.GetColValues(col))
                plotSig.append(GetSigData(col))
                plotNames.append(self.SeqTable.GetColLabelValue(col))
            plotData = np.array(plotData)
            G2plt.PlotSeq(self,plotData,plotSig,plotNames)
        elif rows:
            name = histNames[rows[0]]
            G2plt.PlotCovariance(self,Data=data[name])
               
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataFrame.SetLabel('Sequental refinement results')
    self.dataFrame.CreateStatusBar()
    colLabels = data['varyList']
    Types = len(data['varyList'])*[wg.GRID_VALUE_FLOAT,]
    seqList = [list(data[name]['variables']) for name in histNames]
    self.SeqTable = Table(seqList,colLabels=colLabels,rowLabels=histNames,types=Types)
    self.dataDisplay = GSGrid(parent=self.dataFrame)
    self.dataDisplay.SetTable(self.SeqTable, True)
    self.dataDisplay.EnableEditing(False)
    self.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, Select)
    self.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    self.dataDisplay.SetMargins(0,0)
    self.dataDisplay.AutoSizeColumns(True)
    self.dataFrame.setSizePosLeft([700,350])
                
def UpdateConstraints(self,data):
    if not data:
        data.update({'Hist':[],'HAP':[],'Phase':[]})       #empty dict - fill it
    Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
    Natoms,phaseVary,phaseDict,pawleyLookup,FFtable = G2str.GetPhaseData(Phases,Print=False)        
    hapVary,hapDict,controlDict = G2str.GetHistogramPhaseData(Phases,Histograms,Print=False)
    histVary,histDict,controlDict = G2str.GetHistogramData(Histograms,Print=False)
    Indx = {}
    self.Page = 0
    
    def FindEquivVarb(name,nameList):
        outList = []
        for item in nameList:
            key = item.split(':')[2]
            if key in name and item != name:
                outList.append(item)
        return outList
        
    def SelectVarbs(FrstVarb,varList,legend):
        #future -  add 'all:all:name', '0:all:name', etc. to the varList
        dlg = wx.MultiChoiceDialog(self,'Select more variables:'+legend,FrstVarb+' and:',varList)
        varbs = [FrstVarb,]
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                varbs.append(varList[x])
        dlg.Destroy()
        if len(varbs) > 1:
            return map(list,zip([1.0 for i in range(len(varbs))],varbs))
        else:
            return [[0.0,FrstVarb],]        #setup for "fix" of variable
    
    def OnAddConstraint(event):
        constr = []
        plegend = '\n In p::name'
        hlegend = '\n In :h:name'
        phlegend = '\n In p:h:name'
        for phase in Phases:
            plegend += '\n p:: = '+str(Phases[phase]['pId'])+':: for '+phase
            for histogram in Phases[phase]['Histograms']:
                phlegend += '\n p:h: = '+str(Phases[phase]['pId'])+':'+str(Histograms[histogram]['hId'])+': for '+phase+' in '+histogram
        for histogram in Histograms:
            hlegend += '\n :h: = :'+str(Histograms[histogram]['hId'])+': for '+histogram
        page = self.Page
        if 'Histogram ' in self.dataDisplay.GetPageText(page):
            dlg = wx.SingleChoiceDialog(self,'Select 1st variable:'+hlegend,'Histogram variables:',histVary)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                FrstVarb = histVary[sel]
                moreVarb = FindEquivVarb(FrstVarb,histVary)
                constr = SelectVarbs(FrstVarb,moreVarb,hlegend)
                constr += [0.0,True]        #constant & refine flag
                data['Hist'].append(constr)
            dlg.Destroy()
            UpdateHistConstr()
        elif '/Phase' in self.dataDisplay.GetPageText(page):
            legend = 'Select 1st variable: \n '
            dlg = wx.SingleChoiceDialog(self,'Select 1st variable:'+phlegend,'HAP variables:',hapVary)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                FrstVarb = hapVary[sel]
                moreVarb = FindEquivVarb(FrstVarb,hapVary)
                constr = SelectVarbs(FrstVarb,moreVarb,phlegend)
                constr += [0.0,True]        #constant & refine flag
                data['HAP'].append(constr)
            dlg.Destroy()
            UpdateHAPConstr()
        elif 'Phase' in self.dataDisplay.GetPageText(page):
            #maybe get atom name in phaseVary items?
            for phase in Phases:
                Phase = Phases[phase]
                Atoms = Phase['Atoms']
            dlg = wx.SingleChoiceDialog(self,'Select 1st variable:'+plegend,'Phase variables:',phaseVary)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                FrstVarb = phaseVary[sel]
                moreVarb = FindEquivVarb(FrstVarb,phaseVary)
                constr = SelectVarbs(FrstVarb,moreVarb,plegend+'\nIf none selected then '+FrstVarb+' is held fixed')
                constr += [0.0,True]        #constant & refine flag
                data['Phase'].append(constr)
            dlg.Destroy()
            UpdatePhaseConstr()
            
    def ConstSizer(name,pageDisplay):
        constSizer = wx.FlexGridSizer(1,4,0,0)
        for Id,item in enumerate(data[name]):
            if len(item) < 4:
                constSizer.Add((5,0),0)
                eqString = ' FIXED   '+item[0][1]
                constSizer.Add((5,0),0)
            else:
                constEdit = wx.Button(pageDisplay,-1,'Edit',style=wx.BU_EXACTFIT)
                constEdit.Bind(wx.EVT_BUTTON,OnConstEdit)
                Indx[constEdit.GetId()] = [Id,name]
                constSizer.Add(constEdit)            
                constRef = wx.CheckBox(pageDisplay,-1,label=' Refine?')
                constRef.SetValue(item[-1])
                constRef.Bind(wx.EVT_CHECKBOX,OnConstRef)
                Indx[constRef.GetId()] = item
                constSizer.Add(constRef,0,wx.ALIGN_CENTER_VERTICAL)
                eqString = ''
                for term in item[:-2]:
                    eqString += '%+.3f*%s '%(term[0],term[1])
                eqString += ' = %.3f'%(item[-2])
            constSizer.Add(wx.StaticText(pageDisplay,-1,eqString),0,wx.ALIGN_CENTER_VERTICAL)
            constDel = wx.Button(pageDisplay,-1,'Delete',style=wx.BU_EXACTFIT)
            constDel.Bind(wx.EVT_BUTTON,OnConstDel)
            Indx[constDel.GetId()] = [Id,name]
            constSizer.Add(constDel)
        return constSizer
                
    def OnConstRef(event):
        Obj = event.GetEventObject()
        Indx[Obj.GetId()][-1] = Obj.GetValue()
        
    def OnConstDel(event):
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        del(data[name][Id])
        OnPageChanged(event)        
        
    def OnConstEdit(event):
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        const = data[name][Id][-2]
        items = data[name][Id][:-2]+[[const,'= constant'],[]]
        print items
        dlg = self.SumDialog(self,'Constraint','Enter value for each term in constraint','',items)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                print result
                data[name][Id][:-2] = result[:-1]
                data[name][Id][-2] = result[-2][0]
        finally:
            dlg.Destroy()            
        OnPageChanged(event)        
             
    
    def UpdateHAPConstr():
        HAPConstr.DestroyChildren()
        pageDisplay = wx.Panel(HAPConstr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(ConstSizer('HAP',pageDisplay))
        pageDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        pageDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)
        
    def UpdateHistConstr():
        HistConstr.DestroyChildren()
        pageDisplay = wx.Panel(HistConstr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)        
        mainSizer.Add(ConstSizer('Hist',pageDisplay))
        pageDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        pageDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)
        
    def UpdatePhaseConstr():
        PhaseConstr.DestroyChildren()
        pageDisplay = wx.Panel(PhaseConstr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)        
        mainSizer.Add(ConstSizer('Phase',pageDisplay))
        pageDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        pageDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)
    
    def OnPageChanged(event):
        page = event.GetSelection()
        self.Page = page
        text = self.dataDisplay.GetPageText(page)
        if text == 'Histogram/Phase constraints':
            UpdateHAPConstr()
        elif text == 'Histogram constraints':
            UpdateHistConstr()
        elif text == 'Phase constraints':
            UpdatePhaseConstr()

    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.ConstraintMenu)
    self.dataFrame.SetLabel('Constraints')
    self.dataFrame.CreateStatusBar()
    self.dataFrame.SetMenuBar(self.dataFrame.ConstraintMenu)
    self.dataFrame.Bind(wx.EVT_MENU, OnAddConstraint, id=wxID_CONSTRAINTADD)
    self.dataDisplay = GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())
    
    PhaseConstr = wx.ScrolledWindow(self.dataDisplay)
    self.dataDisplay.AddPage(PhaseConstr,'Phase constraints')
    HAPConstr = wx.ScrolledWindow(self.dataDisplay)
    self.dataDisplay.AddPage(HAPConstr,'Histogram/Phase constraints')
    HistConstr = wx.ScrolledWindow(self.dataDisplay)
    self.dataDisplay.AddPage(HistConstr,'Histogram constraints')
    UpdatePhaseConstr()

    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    
    
def UpdateRestraints(self,data):

    def OnAddRestraint(event):
        page = self.dataDisplay.GetSelection()
        print self.dataDisplay.GetPageText(page)

    def UpdateAtomRestr():
        AtomRestr.DestroyChildren()
        dataDisplay = wx.Panel(AtomRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Atom restraint data:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)


        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)
        
    def UpdatePhaseRestr():
        PhaseRestr.DestroyChildren()
        dataDisplay = wx.Panel(PhaseRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Phase restraint data:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)


        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)
    
    def OnPageChanged(event):
        page = event.GetSelection()
        text = self.dataDisplay.GetPageText(page)
        if text == 'Atom restraints':
            self.dataFrame.SetMenuBar(self.dataFrame.RestraintMenu)
            UpdateAtomRestr()
        elif text == 'Phase restraints':
            UpdatePhaseRestr()
            self.dataFrame.SetMenuBar(self.dataFrame.RestraintMenu)
        event.Skip()

    if self.dataDisplay:
        self.dataDisplay.Destroy()
    self.dataFrame.SetMenuBar(self.dataFrame.RestraintMenu)
    self.dataFrame.SetLabel('restraints')
    self.dataFrame.CreateStatusBar()
    self.dataFrame.Bind(wx.EVT_MENU, OnAddRestraint, id=wxID_RESTRAINTADD)
    self.dataDisplay = GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())
    
    PhaseRestr = wx.ScrolledWindow(self.dataDisplay)
    self.dataDisplay.AddPage(PhaseRestr,'Phase restraints')
    AtomRestr = wx.ScrolledWindow(self.dataDisplay)
    self.dataDisplay.AddPage(AtomRestr,'Atom restraints')
    UpdatePhaseRestr()
#    AtomRestrData = data['AtomRestr']

    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)        
             
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

def GetPatternTreeDataNames(self,dataTypes):
    names = []
    item, cookie = self.PatternTree.GetFirstChild(self.root)        
    while item:
        name = self.PatternTree.GetItemText(item)
        if name[:4] in dataTypes:
            names.append(name)
        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
    return names
                          
def GetPatternTreeItemId(self, parentId, itemText):
    item, cookie = self.PatternTree.GetFirstChild(parentId)
    while item:
        if self.PatternTree.GetItemText(item) == itemText:
            return item
        item, cookie = self.PatternTree.GetNextChild(parentId, cookie)
    return 0                

def MovePatternTreeToGrid(self,item):
    
#    print self.PatternTree.GetItemText(item)
    
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
            if not data:           #fill in defaults
                data = {
                    #least squares controls
                    'deriv type':'analytic','min dM/M':0.0001,'shift factor':1.0,
                    #Fourier controls
                    'mapType':'Fobs','d-max':100.,'d-min':0.2,'histograms':[],
                    'stepSize':[0.5,0.5,0.5],'minX':[0.,0.,0.],'maxX':[1.0,1.0,1.0],
                    #distance/angle controls
                    'distMax':0.0,'angleMax':0.0,'useMapPeaks':False}
                self.PatternTree.SetItemPyData(item,data)                             
            self.Refine.Enable(True)
            self.SeqRefine.Enable(True)
            UpdateControls(self,data)
        elif self.PatternTree.GetItemText(item) == 'Sequental results':
            data = self.PatternTree.GetItemPyData(item)
            UpdateSeqResults(self,data)            
        elif self.PatternTree.GetItemText(item) == 'Covariance':
            data = self.PatternTree.GetItemPyData(item)
            G2plt.PlotCovariance(self)
        elif self.PatternTree.GetItemText(item) == 'Constraints':
            data = self.PatternTree.GetItemPyData(item)
            UpdateConstraints(self,data)
        elif self.PatternTree.GetItemText(item) == 'Restraints':
            data = self.PatternTree.GetItemPyData(item)
            UpdateRestraints(self,data)
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
        elif 'PDF' in self.PatternTree.GetItemText(item):
            self.PatternId = item
            self.ExportPDF.Enable(True)
            G2plt.PlotISFG(self,type='S(Q)')
            
    elif 'I(Q)' in self.PatternTree.GetItemText(item):
        self.PickId = item
        self.PatternId = self.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(self,type='I(Q)',newPlot=True)
    elif 'S(Q)' in self.PatternTree.GetItemText(item):
        self.PickId = item
        self.PatternId = self.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(self,type='S(Q)',newPlot=True)
    elif 'F(Q)' in self.PatternTree.GetItemText(item):
        self.PickId = item
        self.PatternId = self.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(self,type='F(Q)',newPlot=True)
    elif 'G(R)' in self.PatternTree.GetItemText(item):
        self.PickId = item
        self.PatternId = self.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(self,type='G(R)',newPlot=True)            
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
    elif self.PatternTree.GetItemText(item) == 'PDF Controls':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.ExportPDF.Enable(True)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(self,data)
        G2plt.PlotISFG(self,type='I(Q)')
        G2plt.PlotISFG(self,type='S(Q)')
        G2plt.PlotISFG(self,type='F(Q)')
        G2plt.PlotISFG(self,type='G(R)')
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
                'DisplaceY':[0.0,False],'Diffuse':[],'Temperature':300.,'Pressure':1.0,'Humidity':0.0,'Voltage':0.0,
                'Force':0.0,'Gonio. radius':200.0}
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
            data.append([0,0.0,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, zero value, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            self.PatternTree.SetItemPyData(item,data)                             
        G2pdG.UpdateUnitCellsGrid(self,data)
        if 'PKS' in self.PatternTree.GetItemText(self.PatternId):
            G2plt.PlotPowderLines(self)
        else:
            G2plt.PlotPatterns(self)
    elif self.PatternTree.GetItemText(item) == 'Reflection Lists':
        self.PatternId = self.PatternTree.GetItemParent(item)
        self.PickId = item
        data = self.PatternTree.GetItemPyData(item)
        self.RefList = ''
        if len(data):
            self.RefList = data.keys()[0]
        G2pdG.UpdateReflectionGrid(self,data)
        G2plt.PlotPatterns(self)
     
