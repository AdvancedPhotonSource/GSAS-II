#GSASII plotting routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import math
import time
import copy
import os.path
import numpy as np
import numpy.linalg as nl
import wx
import wx.aui
import wx.glcanvas
import matplotlib as mpl
import mpl_toolkits.mplot3d.axes3d as mp3d
import GSASIIpath
import GSASIIgrid as G2gd
import GSASIIimage as G2img
import GSASIIpwd as G2pwd
import GSASIIIO as G2IO
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import pytexture as ptx
from  OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GLE import *
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar

# useful degree trig functions
sind = lambda x: math.sin(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
acosd = lambda x: 180.*math.acos(x)/math.pi
atan2d = lambda x,y: 180.*math.atan2(y,x)/math.pi
atand = lambda x: 180.*math.atan(x)/math.pi
# numpy versions
npsind = lambda x: np.sin(x*np.pi/180.)
npcosd = lambda x: np.cos(x*np.pi/180.)
npacosd = lambda x: 180.*np.arccos(x)/np.pi
npasind = lambda x: 180.*np.arcsin(x)/np.pi
npatan2d = lambda x,y: 180.*np.arctan2(x,y)/np.pi
    
class G2PlotMpl(wx.Panel):    
    def __init__(self,parent,id=-1,dpi=None,**kwargs):
        wx.Panel.__init__(self,parent,id=id,**kwargs)
        mpl.rcParams['legend.fontsize'] = 8
        self.figure = mpl.figure.Figure(dpi=dpi,figsize=(5,7))
        self.canvas = Canvas(self,-1,self.figure)
        self.toolbar = Toolbar(self.canvas)

        self.toolbar.Realize()
        
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar,0,wx.LEFT|wx.EXPAND)
        self.SetSizer(sizer)
        
class G2PlotOgl(wx.Panel):
    def __init__(self,parent,id=-1,dpi=None,**kwargs):
        self.figure = wx.Panel.__init__(self,parent,id=id,**kwargs)
        self.canvas = wx.glcanvas.GLCanvas(self,-1,**kwargs)
        self.camera = {}
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        self.SetSizer(sizer)
        
class G2Plot3D(wx.Panel):
    def __init__(self,parent,id=-1,dpi=None,**kwargs):
        wx.Panel.__init__(self,parent,id=id,**kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi,figsize=(6,6))
        self.canvas = Canvas(self,-1,self.figure)
        self.toolbar = Toolbar(self.canvas)

        self.toolbar.Realize()
        
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar,0,wx.LEFT|wx.EXPAND)
        self.SetSizer(sizer)
                              
class G2PlotNoteBook(wx.Panel):
    def __init__(self,parent,id=-1):
        wx.Panel.__init__(self,parent,id=id)
        #so one can't delete a plot page!!
        self.nb = wx.aui.AuiNotebook(self, \
            style=wx.aui.AUI_NB_DEFAULT_STYLE ^ wx.aui.AUI_NB_CLOSE_ON_ACTIVE_TAB)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb,1,wx.EXPAND)
        self.SetSizer(sizer)
        self.status = parent.CreateStatusBar()
        self.status.SetFieldsCount(2)
        self.status.SetStatusWidths([150,-1])
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        
        self.plotList = []
            
    def addMpl(self,name=""):
        page = G2PlotMpl(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
    def add3D(self,name=""):
        page = G2Plot3D(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
    def addOgl(self,name=""):
        page = G2PlotOgl(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
    def Delete(self,name):
        try:
            item = self.plotList.index(name)
            del self.plotList[item]
            self.nb.DeletePage(item)
        except ValueError:          #no plot of this name - do nothing
            return      
                
    def clear(self):
        while self.nb.GetPageCount():
            self.nb.DeletePage(0)
        self.plotList = []
        self.status.DestroyChildren()
        
    def Rename(self,oldName,newName):
        try:
            item = self.plotList.index(oldName)
            self.plotList[item] = newName
            self.nb.SetPageText(item,newName)
        except ValueError:          #no plot of this name - do nothing
            return      
        
    def OnPageChanged(self,event):        
        if self.plotList:
            self.status.SetStatusText('Better to select this from GSAS-II data tree',1)
        self.status.DestroyChildren()                           #get rid of special stuff on status bar
        
def PlotSngl(self,newPlot=False):
    '''Single crystal structure factor plotting package - displays zone of reflections as rings proportional
        to F, F**2, etc. as requested
    '''
    from matplotlib.patches import Circle
    global HKL,HKLF

    def OnSCMotion(event):
        xpos = event.xdata
        if xpos:
            xpos = round(xpos)                                        #avoid out of frame mouse position
            ypos = round(event.ydata)
            zpos = Data['Layer']
            if '100' in Data['Zone']:
                HKLtxt = '(%3d,%3d,%3d)'%(zpos,xpos,ypos)
            elif '010' in Data['Zone']:
                HKLtxt = '(%3d,%3d,%3d)'%(xpos,zpos,ypos)
            elif '001' in Data['Zone']:
                HKLtxt = '(%3d,%3d,%3d)'%(xpos,ypos,zpos)
            Page.canvas.SetToolTipString(HKLtxt)
            self.G2plotNB.status.SetFields(['HKL = '+HKLtxt,''])
                
    def OnSCPick(event):
        zpos = Data['Layer']
        pos = event.artist.center
        if '100' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(zpos,pos[0],pos[1]))
            hkl = np.array([zpos,pos[0],pos[1]])
        elif '010' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],zpos,pos[1]))
            hkl = np.array([pos[0],zpos,pos[1]])
        elif '001' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],pos[1],zpos))
            hkl = np.array([pos[0],pos[1],zpos])
        h,k,l = hkl
        hklf = HKLF[np.where(np.all(HKL-hkl == [0,0,0],axis=1))]
        if len(hklf):
            Fosq,sig,Fcsq = hklf[0]
            HKLtxt = '(%3d,%3d,%3d %.2f %.3f %.2f %.2f)'%(h,k,l,Fosq,sig,Fcsq,(Fosq-Fcsq)/(scale*sig))
            self.G2plotNB.status.SetFields(['','HKL, Fosq, sig, Fcsq, delFsq/sig = '+HKLtxt])
                                 
    def OnSCKeyPress(event):
        print event.key

    try:
        plotNum = self.G2plotNB.plotList.index('Structure Factors')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
    except ValueError:
        Plot = self.G2plotNB.addMpl('Structure Factors').gca()
        plotNum = self.G2plotNB.plotList.index('Structure Factors')
        Page = self.G2plotNB.nb.GetPage(plotNum)
#        Page.canvas.mpl_connect('key_press_event', OnSCKeyPress)
        Page.canvas.mpl_connect('pick_event', OnSCPick)
        Page.canvas.mpl_connect('motion_notify_event', OnSCMotion)
    Page.SetFocus()
    
    Plot.set_aspect(aspect='equal')
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
    Plot.set_title(self.PatternTree.GetItemText(self.Sngl)[5:])
    HKL = []
    HKLF = []
    for H,Fosq,sig,Fcsq,x,x,x in HKLref:
        HKL.append(H)
        HKLF.append([Fosq,sig,Fcsq])
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
                Plot.add_artist(Circle(xy,radius=A,ec='g',fc='w',picker=3))
            if B:
                Plot.add_artist(Circle(xy,radius=B,ec='b',fc='w'))
                radius = C
                if radius > 0:
                    if A > B:
                        Plot.add_artist(Circle(xy,radius=radius,ec='g',fc='g'))
                    else:                    
                        Plot.add_artist(Circle(xy,radius=radius,ec='r',fc='r'))
    HKL = np.array(HKL,dtype=np.int)
    HKLF = np.array(HKLF)
    Plot.set_xlabel(xlabel[izone]+str(Data['Layer']),fontsize=12)
    Plot.set_ylabel(ylabel[izone],fontsize=12)
    Plot.set_xlim((HKLmin[pzone[izone][0]],HKLmax[pzone[izone][0]]))
    Plot.set_ylim((HKLmin[pzone[izone][1]],HKLmax[pzone[izone][1]]))
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
       
def PlotPatterns(self,newPlot=False):
    '''Powder pattern plotting package - displays single or multiple powder patterns as intensity vs
    2-theta or q (future TOF). Can display multiple patterns as "waterfall plots" or contour plots. Log I 
    plotting available.
    '''
    global HKL
    
    def OnPick(event):
        if self.itemPicked is not None: return
        PatternId = self.PatternId
        try:
            Values,Names = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1::2]
        except TypeError:
            return
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        PickId = self.PickId
        pick = event.artist
        mouse = event.mouseevent       
        xpos = pick.get_xdata()
        ypos = pick.get_ydata()
        ind = event.ind
        xy = list(zip(np.take(xpos,ind),np.take(ypos,ind))[0])
        if self.PatternTree.GetItemText(PickId) == 'Peak List':
            if ind.all() != [0]:                                    #picked a data point
                if 'C' in Parms['Type']:                            #CW data - TOF later in an elif
                    ins = [Parms[x] for x in ['U','V','W','X','Y']]
                    if self.qPlot:                              #qplot - convert back to 2-theta
                        xy[0] = 2.0*asind(xy[0]*wave/(4*math.pi))
                    sig = ins[0]*tand(xy[0]/2.0)**2+ins[1]*tand(xy[0]/2.0)+ins[2]
                    gam = ins[3]/cosd(xy[0]/2.0)+ins[4]*tand(xy[0]/2.0)           
                    data = self.PatternTree.GetItemPyData(self.PickId)
                    XY = [xy[0],0, xy[1],1, sig,0, gam,0]       #default refine intensity 1st
                data.append(XY)
                G2pdG.UpdatePeakGrid(self,data)
                PlotPatterns(self)
            else:                                                   #picked a peak list line
                self.itemPicked = pick
        elif self.PatternTree.GetItemText(PickId) == 'Limits':
            if ind.all() != [0]:                                    #picked a data point
                LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
                data = self.PatternTree.GetItemPyData(LimitId)
                if 'C' in Parms['Type']:                            #CW data - TOF later in an elif
                    if self.qPlot:                              #qplot - convert back to 2-theta
                        xy[0] = 2.0*asind(xy[0]*wave/(4*math.pi))
                if mouse.button==1:
                    data[1][0] = min(xy[0],data[1][1])
                if mouse.button==3:
                    data[1][1] = max(xy[0],data[1][0])
                self.PatternTree.SetItemPyData(LimitId,data)
                G2pdG.UpdateLimitsGrid(self,data)
                PlotPatterns(self)
            else:                                                   #picked a limit line
                self.itemPicked = pick
        
    def OnPlotKeyPress(event):
        newPlot = False
        if event.key == 'w':
            if self.Weight:
                self.Weight = False
            else:
                self.Weight = True
            print 'plot weighting:',self.Weight
        elif event.key == 'n':
            if self.Contour:
                pass
            else:
                if self.logPlot:
                    self.logPlot = False
                else:
                    self.Offset[0] = 0
                    self.logPlot = True
        elif event.key == 'u':
            if self.Contour:
                self.Cmax = min(1.0,self.Cmax*1.2)
            elif self.logPlot:
                pass
            elif self.Offset[0] < 100.:
                self.Offset[0] += 1.
        elif event.key == 'd':
            if self.Contour:
                self.Cmax = max(0.0,self.Cmax*0.8)
            elif self.logPlot:
                pass
            elif self.Offset[0] > 0.:
                self.Offset[0] -= 1.
        elif event.key == 'l':
            self.Offset[1] -= 1.
        elif event.key == 'r':
            self.Offset[1] += 1.
        elif event.key == 'o':
            self.Offset = [0,0]
        elif event.key == 'c':
            newPlot = True
            if self.Contour:
                self.Contour = False
            else:
                self.Contour = True
                self.SinglePlot = False
                self.Offset = [0.,0.]
        elif event.key == 'q':
            newPlot = True
            if self.qPlot:
                self.qPlot = False
            else:
                self.qPlot = True
        elif event.key == 's':
            if self.Contour:
                choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
                choice.sort()
                dlg = wx.SingleChoiceDialog(self,'Select','Color scheme',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    self.ContourColor = choice[sel]
                else:
                    self.ContourColor = 'Paired'
                dlg.Destroy()
            else:                
                if self.SinglePlot:
                    self.SinglePlot = False
                else:
                    self.SinglePlot = True
        elif event.key == '+':
            if self.PickId:
                self.PickId = False
        elif event.key == 'i':                  #for smoothing contour plot
            choice = ['nearest','bilinear','bicubic','spline16','spline36','hanning',
               'hamming','hermite','kaiser','quadric','catrom','gaussian','bessel',
               'mitchell','sinc','lanczos']
            dlg = wx.SingleChoiceDialog(self,'Select','Interpolation',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                self.Interpolate = choice[sel]
            else:
                self.Interpolate = 'nearest'
            dlg.Destroy()
            
        PlotPatterns(self,newPlot=newPlot)
        
    def OnKeyBox(event):
        if self.G2plotNB.nb.GetSelection() == self.G2plotNB.plotList.index('Powder Patterns'):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            OnPlotKeyPress(event)
                        
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                Values,Names = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1::2]
                Parms = dict(zip(Names,Values))
                try:
                    wave = Parms['Lam']
                except KeyError:
                    wave = Parms['Lam1']
                if self.qPlot:
                    xpos = 2.0*asind(xpos*wave/(4*math.pi))
                dsp = 0.0
                if abs(xpos) > 0.:                  #avoid possible singularity at beam center
                    dsp = wave/(2.*sind(abs(xpos)/2.0))
                if self.Contour:
                    self.G2plotNB.status.SetStatusText('2-theta =%9.3f d =%9.5f pattern ID =%5d'%(xpos,dsp,int(ypos)),1)
                else:
                    self.G2plotNB.status.SetStatusText('2-theta =%9.3f d =%9.5f Intensity =%9.1f'%(xpos,dsp,ypos),1)
                if self.itemPicked:
                    Page.canvas.SetToolTipString('%9.3f'%(xpos))
                if self.PickId and self.PatternTree.GetItemText(self.PickId) in ['Index Peak List','Unit Cells List']:
                    found = []
                    if len(HKL):
                        view = Page.toolbar._views.forward()[0][:2]
                        wid = view[1]-view[0]
                        found = HKL[np.where(np.fabs(HKL.T[5]-xpos) < 0.002*wid)]
                    if len(found):
                        h,k,l = found[0][:3] 
                        Page.canvas.SetToolTipString('%d,%d,%d'%(int(h),int(k),int(l)))
                    else:
                        Page.canvas.SetToolTipString('')

            except TypeError:
                self.G2plotNB.status.SetStatusText('Select PWDR powder pattern first',1)
                                                   
    def OnRelease(event):
        if self.itemPicked is None: return
        Values,Names = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,self.PatternId, 'Instrument Parameters'))[1::2]
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            lines = []
            for line in self.Lines: lines.append(line.get_xdata()[0])
            lineNo = lines.index(self.itemPicked.get_xdata()[0])
            if  lineNo in [0,1]:
                LimitId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Limits')
                data = self.PatternTree.GetItemPyData(LimitId)
#                print 'limits',xpos
                if self.qPlot:
                    data[1][lineNo] = 2.0*asind(wave*xpos/(4*math.pi))
                else:
                    data[1][lineNo] = xpos
                self.PatternTree.SetItemPyData(LimitId,data)
                if self.PatternTree.GetItemText(self.PickId) == 'Limits':
                    G2pdG.UpdateLimitsGrid(self,data)
            else:
                PeakId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List')
                data = self.PatternTree.GetItemPyData(PeakId)
#                print 'peaks',xpos
                if event.button == 3:
                    del data[lineNo-2]
                else:
                    if self.qPlot:
                        data[lineNo-2][0] = 2.0*asind(wave*xpos/(4*math.pi))
                    else:
                        data[lineNo-2][0] = xpos
                self.PatternTree.SetItemPyData(PeakId,data)
                G2pdG.UpdatePeakGrid(self,data)
        PlotPatterns(self)
        self.itemPicked = None    

    xylim = []
    try:
        plotNum = self.G2plotNB.plotList.index('Powder Patterns')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
    except ValueError:
        newPlot = True
        self.Cmax = 1.0
        Plot = self.G2plotNB.addMpl('Powder Patterns').gca()
        plotNum = self.G2plotNB.plotList.index('Powder Patterns')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.canvas.mpl_connect('pick_event', OnPick)
        Page.canvas.mpl_connect('button_release_event', OnRelease)
    Page.SetFocus()
    self.G2plotNB.status.DestroyChildren()
    if self.Contour:
        Choice = (' key press','d: lower contour max','u: raise contour max',
            'i: interpolation method','s: color scheme','c: contour off')
    else:
        if self.logPlot:
            Choice = (' key press','n: log(I) off','l: offset left','r: offset right',
                'c: contour on','q: toggle q plot','s: toggle single plot','+: no selection')
        else:
            Choice = (' key press','l: offset left','r: offset right','d: offset down',
                'u: offset up','o: reset offset','n: log(I) on','c: contour on',
                'q: toggle q plot','s: toggle single plot','+: no selection')
    cb = wx.ComboBox(self.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=Choice)
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' key press')
    
    PickId = self.PickId
    PatternId = self.PatternId
    colors=['b','g','r','c','m','k']
    Lines = []
    if self.SinglePlot:
        Pattern = self.PatternTree.GetItemPyData(PatternId)
        Pattern.append(self.PatternTree.GetItemText(PatternId))
        PlotList = [Pattern,]
        ParmList = [self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,
            self.PatternId, 'Instrument Parameters'))[1],]
    else:        
        PlotList = []
        ParmList = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            if 'PWDR' in self.PatternTree.GetItemText(item):
                Pattern = self.PatternTree.GetItemPyData(item)
                if len(Pattern) < 3:                    # put name on end if needed
                    Pattern.append(self.PatternTree.GetItemText(item))
                PlotList.append(Pattern)
                ParmList.append(self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,
                    item,'Instrument Parameters'))[1])
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
    Ymax = 1.0
    lenX = 0
    HKL = np.array(self.HKL)
    for Pattern in PlotList:
        xye = Pattern[1]
        Ymax = max(Ymax,max(xye[1]))
    offset = self.Offset[0]*Ymax/100.0
    Title = 'Powder Patterns: '+os.path.split(self.GSASprojectfile)[1]
    if self.logPlot:
        Title = 'log('+Title+')'
    Plot.set_title(Title)
    if self.qPlot:
        Plot.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
    else:        
        Plot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    Plot.set_ylabel('Intensity',fontsize=12)
    if self.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    for N,Pattern in enumerate(PlotList):
        Parms = ParmList[N]
        ifpicked = False
        LimitId = 0
        xye = np.array(Pattern[1])
        if PickId:
            ifpicked = Pattern[2] == self.PatternTree.GetItemText(PatternId)
            LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
        if self.qPlot:
            Id = G2gd.GetPatternTreeItemId(self,self.root, Pattern[2])
            X = 4*np.pi*npsind(xye[0]/2.0)/Parms[1]
        else:
            X = xye[0]
        if not lenX:
            lenX = len(X)           
        Y = xye[1]+offset*N
        if LimitId:
            limits = np.array(self.PatternTree.GetItemPyData(LimitId))
            if self.qPlot:
                limits = 4*np.pi*npsind(limits/2.0)/Parms[1]
            Lines.append(Plot.axvline(limits[1][0],color='g',dashes=(5,5),picker=3.))    
            Lines.append(Plot.axvline(limits[1][1],color='r',dashes=(5,5),picker=3.))                    
        if self.Contour:
            
            if lenX == len(X):
                ContourY.append(N)
                ContourZ.append(Y)
                ContourX = X
                Nseq += 1
                Plot.set_ylabel('Data sequence',fontsize=12)
        else:
            X += self.Offset[1]*.005*N
            if ifpicked:
                Z = xye[3]+offset*N
                W = xye[4]+offset*N
                D = xye[5]+offset*N-Ymax*.02
                if self.Weight:
                    W2 = np.sqrt(xye[2])
                    D *= W2-Ymax*.02
                if self.logPlot:
                    Plot.semilogy(X,Y,colors[N%6]+'+',picker=3.,clip_on=False,nonposy='mask')
                    Plot.semilogy(X,Z,colors[(N+1)%6],picker=False,nonposy='mask')
                    Plot.semilogy(X,W,colors[(N+2)%6],picker=False,nonposy='mask')
                else:
                    Plot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                    Plot.plot(X,Z,colors[(N+1)%6],picker=False)
                    Plot.plot(X,W,colors[(N+2)%6],picker=False)
                    Plot.plot(X,D,colors[(N+3)%6],picker=False)
                    Plot.axhline(0.,color=wx.BLACK)
                Page.canvas.SetToolTipString('')
                if self.PatternTree.GetItemText(PickId) == 'Peak List':
                    tip = 'On data point: Pick peak - L or R MB. On line: L-move, R-delete'
                    Page.canvas.SetToolTipString(tip)
                    data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
                    for item in data:
                        if self.qPlot:
                            Lines.append(Plot.axvline(4*math.pi*sind(item[0]/2.)/Parms[1],color=colors[N%6],picker=2.))
                        else:
                            Lines.append(Plot.axvline(item[0],color=colors[N%6],picker=2.))
                if self.PatternTree.GetItemText(PickId) == 'Limits':
                    tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                    Page.canvas.SetToolTipString(tip)
                    data = self.LimitsTable.GetData()
            else:
                if self.logPlot:
                    Plot.semilogy(X,Y,colors[N%6],picker=False,nonposy='clip')
                else:
                    Plot.plot(X,Y,colors[N%6],picker=False)
    if PickId and self.PatternTree.GetItemText(PickId) in ['Index Peak List','Unit Cells List']:
        Values,Names = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))[1::2]
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        peaks = np.array((self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))))
        for peak in peaks:
            if self.qPlot:
                Plot.axvline(4*np.pi*sind(peak[0]/2.0)/wave,color='b')
            else:
                Plot.axvline(peak[0],color='b')
        for hkl in self.HKL:
            if self.qPlot:
                Plot.axvline(4*np.pi*sind(hkl[5]/2.0)/wave,color='r',dashes=(5,5))
            else:
                Plot.axvline(hkl[5],color='r',dashes=(5,5))
    if self.Contour:
        acolor = mpl.cm.get_cmap(self.ContourColor)
        Img = Plot.imshow(ContourZ,cmap=acolor,vmin=0,vmax=Ymax*self.Cmax,interpolation=self.Interpolate, 
            extent=[ContourX[0],ContourX[-1],ContourY[0],ContourY[-1]],aspect='auto',origin='lower')
        Page.figure.colorbar(Img)
    else:
        self.Lines = Lines
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
    self.Pwdr = True
    
def PlotISFG(self,newPlot=False,type=''):
    ''' PLotting package for PDF analysis; displays I(q), S(q), F(q) and G(r) as single 
    or multiple plots with waterfall and contour plots as options
    '''
    if not type:
        type = self.G2plotNB.plotList[self.G2plotNB.nb.GetSelection()]
    if type not in ['I(Q)','S(Q)','F(Q)','G(R)']:
        return
    superMinusOne = unichr(0xaf)+unichr(0xb9)
    
    def OnPlotKeyPress(event):
        newPlot = False
        if event.key == 'u':
            if self.Contour:
                self.Cmax = min(1.0,self.Cmax*1.2)
            elif self.Offset[0] < 100.:
                self.Offset[0] += 1.
        elif event.key == 'd':
            if self.Contour:
                self.Cmax = max(0.0,self.Cmax*0.8)
            elif self.Offset[0] > 0.:
                self.Offset[0] -= 1.
        elif event.key == 'l':
            self.Offset[1] -= 1.
        elif event.key == 'r':
            self.Offset[1] += 1.
        elif event.key == 'o':
            self.Offset = [0,0]
        elif event.key == 'c':
            newPlot = True
            if self.Contour:
                self.Contour = False
            else:
                self.Contour = True
                self.SinglePlot = False
                self.Offset = [0.,0.]
        elif event.key == 's':
            if self.Contour:
                choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
                choice.sort()
                dlg = wx.SingleChoiceDialog(self,'Select','Color scheme',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    self.ContourColor = choice[sel]
                else:
                    self.ContourColor = 'Paired'
                dlg.Destroy()
            else:                
                if self.SinglePlot:
                    self.SinglePlot = False
                else:
                    self.SinglePlot = True
        elif event.key == 'i':                  #for smoothing contour plot
            choice = ['nearest','bilinear','bicubic','spline16','spline36','hanning',
               'hamming','hermite','kaiser','quadric','catrom','gaussian','bessel',
               'mitchell','sinc','lanczos']
            dlg = wx.SingleChoiceDialog(self,'Select','Interpolation',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                self.Interpolate = choice[sel]
            else:
                self.Interpolate = 'nearest'
            dlg.Destroy()
        elif event.key == 't' and not self.Contour:
            if self.Legend:
                self.Legend = False
            else:
                self.Legend = True
            
            
        PlotISFG(self,newPlot=newPlot,type=type)
        
    def OnKeyBox(event):
        if self.G2plotNB.nb.GetSelection() == self.G2plotNB.plotList.index(type):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            OnPlotKeyPress(event)
                        
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                if self.Contour:
                    self.G2plotNB.status.SetStatusText('R =%.3fA pattern ID =%5d'%(xpos,int(ypos)),1)
                else:
                    self.G2plotNB.status.SetStatusText('R =%.3fA %s =%.2f'%(xpos,type,ypos),1)                   
            except TypeError:
                self.G2plotNB.status.SetStatusText('Select '+type+' pattern first',1)
    
    xylim = []
    try:
        plotNum = self.G2plotNB.plotList.index(type)
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        newPlot = True
        self.Cmax = 1.0
        Plot = self.G2plotNB.addMpl(type).gca()
        plotNum = self.G2plotNB.plotList.index(type)
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
    
    Page.SetFocus()
    self.G2plotNB.status.DestroyChildren()
    if self.Contour:
        Choice = (' key press','d: lower contour max','u: raise contour max',
            'i: interpolation method','s: color scheme','c: contour off')
    else:
        Choice = (' key press','l: offset left','r: offset right','d: offset down','u: offset up',
            'o: reset offset','t: toggle legend','c: contour on','s: toggle single plot')
    cb = wx.ComboBox(self.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=Choice)
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' key press')
    PatternId = self.PatternId
    PickId = self.PickId
    Plot.set_title(type)
    if type == 'G(R)':
        Plot.set_xlabel(r'$R,\AA$',fontsize=14)
    else:
        Plot.set_xlabel(r'$Q,\AA$'+superMinusOne,fontsize=14)
    Plot.set_ylabel(r''+type,fontsize=14)
    colors=['b','g','r','c','m','k']
    name = self.PatternTree.GetItemText(PatternId)[4:]
    Pattern = []    
    if self.SinglePlot:
        name = self.PatternTree.GetItemText(PatternId)
        name = type+name[4:]
        Id = G2gd.GetPatternTreeItemId(self,PatternId,name)
        Pattern = self.PatternTree.GetItemPyData(Id)
        if Pattern:
            Pattern.append(name)
        PlotList = [Pattern,]
    else:
        PlotList = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            if 'PDF' in self.PatternTree.GetItemText(item):
                name = type+self.PatternTree.GetItemText(item)[4:]
                Id = G2gd.GetPatternTreeItemId(self,item,name)
                Pattern = self.PatternTree.GetItemPyData(Id)
                if Pattern:
                    Pattern.append(name)
                    PlotList.append(Pattern)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
    PDFdata = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'PDF Controls'))
    numbDen = G2pwd.GetNumDensity(PDFdata['ElList'],PDFdata['Form Vol'])
    Xb = [0.,10.]
    Yb = [0.,-40.*np.pi*numbDen]
    Ymax = 1.0
    lenX = 0
    for Pattern in PlotList:
        xye = Pattern[1]
        Ymax = max(Ymax,max(xye[1]))
    offset = self.Offset[0]*Ymax/100.0
    if self.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    for N,Pattern in enumerate(PlotList):
        xye = Pattern[1]
        if PickId:
            ifpicked = Pattern[2] == self.PatternTree.GetItemText(PatternId)
        X = xye[0]
        if not lenX:
            lenX = len(X)           
        Y = xye[1]+offset*N
        if self.Contour:
            if lenX == len(X):
                ContourY.append(N)
                ContourZ.append(Y)
                ContourX = X
                Nseq += 1
                Plot.set_ylabel('Data sequence',fontsize=12)
        else:
            X = xye[0]+self.Offset[1]*.005*N
            if ifpicked:
                Plot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                Page.canvas.SetToolTipString('')
            else:
                if self.Legend:
                    Plot.plot(X,Y,colors[N%6],picker=False,label='Azm:'+Pattern[2].split('=')[1])
                else:
                    Plot.plot(X,Y,colors[N%6],picker=False)
            if type == 'G(R)':
                Plot.plot(Xb,Yb,color='k',dashes=(5,5))
            elif type == 'F(Q)':
                Plot.axhline(0.,color=wx.BLACK)
            elif type == 'S(Q)':
                Plot.axhline(1.,color=wx.BLACK)
    if self.Contour:
        acolor = mpl.cm.get_cmap(self.ContourColor)
        Img = Plot.imshow(ContourZ,cmap=acolor,vmin=0,vmax=Ymax*self.Cmax,interpolation=self.Interpolate, 
            extent=[ContourX[0],ContourX[-1],ContourY[0],ContourY[-1]],aspect='auto',origin='lower')
        Page.figure.colorbar(Img)
    elif self.Legend:
        Plot.legend(loc='best')
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
        
def PlotXY(self,XY,newPlot=False,type=''):
    '''simple plot of xy data, used for diagnostic purposes
    '''
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                self.G2plotNB.status.SetStatusText('X =%9.3f %s =%9.3f'%(xpos,type,ypos),1)                   
            except TypeError:
                self.G2plotNB.status.SetStatusText('Select '+type+' pattern first',1)

    try:
        plotNum = self.G2plotNB.plotList.index(type)
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        newPlot = True
        Plot = self.G2plotNB.addMpl(type).gca()
        plotNum = self.G2plotNB.plotList.index(type)
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
    
    Page.SetFocus()
    self.G2plotNB.status.DestroyChildren()
    Plot.set_title(type)
    Plot.set_xlabel(r'X',fontsize=14)
    Plot.set_ylabel(r''+type,fontsize=14)
    colors=['b','g','r','c','m','k']
    Ymax = 1.0
    lenX = 0
    X,Y = XY[:2]
    Ymax = max(Ymax,max(Y))
    Plot.plot(X,Y,'k',picker=False)
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()

def PlotPowderLines(self):
    ''' plotting of powder lines (i.e. no powder pattern) as sticks
    '''
    global HKL

    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            self.G2plotNB.status.SetFields(['','2-theta =%9.3f '%(xpos,)])
            if self.PickId and self.PatternTree.GetItemText(self.PickId) in ['Index Peak List','Unit Cells List']:
                found = []
                if len(HKL):
                    view = Page.toolbar._views.forward()[0][:2]
                    wid = view[1]-view[0]
                    found = HKL[np.where(np.fabs(HKL.T[5]-xpos) < 0.002*wid)]
                if len(found):
                    h,k,l = found[0][:3] 
                    Page.canvas.SetToolTipString('%d,%d,%d'%(int(h),int(k),int(l)))
                else:
                    Page.canvas.SetToolTipString('')

    try:
        plotNum = self.G2plotNB.plotList.index('Powder Lines')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        Plot = self.G2plotNB.addMpl('Powder Lines').gca()
        plotNum = self.G2plotNB.plotList.index('Powder Lines')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        
    Page.SetFocus()
    Plot.set_title('Powder Pattern Lines')
    Plot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    PickId = self.PickId
    PatternId = self.PatternId
    peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
    for peak in peaks:
        Plot.axvline(peak[0],color='b')
    HKL = np.array(self.HKL)
    for hkl in self.HKL:
        Plot.axvline(hkl[5],color='r',dashes=(5,5))
    xmin = peaks[0][0]
    xmax = peaks[-1][0]
    delt = xmax-xmin
    xlim = [max(0,xmin-delt/20.),min(180.,xmax+delt/20.)]
    Plot.set_xlim(xlim)
    Page.canvas.draw()
    Page.toolbar.push_current()

def PlotPeakWidths(self):
    ''' Plotting of instrument broadening terms as function of 2-theta (future TOF)
    Seen when "Instrument Parameters" chosen from powder pattern data tree
    '''
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
        if len(instParms[1]) == 13:
            GU,GV,GW,LX,LY = instParms[0][6:11]
        else:
            GU,GV,GW,LX,LY = instParms[0][4:9]
    peakID = G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List')
    if peakID:
        peaks = self.PatternTree.GetItemPyData(peakID)
    else:
        peaks = []
    
    try:
        plotNum = self.G2plotNB.plotList.index('Peak Widths')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        Plot = self.G2plotNB.addMpl('Peak Widths').gca()
        plotNum = self.G2plotNB.plotList.index('Peak Widths')
        Page = self.G2plotNB.nb.GetPage(plotNum)
    Page.SetFocus()
    
    Page.canvas.SetToolTipString('')
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
    V = []
    sig = lambda Th,U,V,W: 1.17741*math.sqrt(U*tand(Th)**2+V*tand(Th)+W)*math.pi/18000.
    gam = lambda Th,X,Y: (X/cosd(Th)+Y*tand(Th))*math.pi/18000.
    gamFW = lambda s,g: math.exp(math.log(s**5+2.69269*s**4*g+2.42843*s**3*g**2+4.47163*s**2*g**3+0.07842*s*g**4+g**5)/5.)
    gamFW2 = lambda s,g: math.sqrt(s**2+(0.4654996*g)**2)+.5345004*g  #Ubaldo Bafile - private communication
    for theta in thetas:
        X.append(4.0*math.pi*sind(theta)/lam)              #q
        s = sig(theta,GU,GV,GW)
        g = gam(theta,LX,LY)
        G = gamFW(g,s)
        H = gamFW2(g,s)
        Y.append(s/tand(theta))
        Z.append(g/tand(theta))
        W.append(G/tand(theta))
        V.append(H/tand(theta))
    Plot.set_title('Instrument and sample peak widths')
    Plot.set_ylabel(r'$\Delta q/q, \Delta d/d$',fontsize=14)
    Plot.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
    Plot.plot(X,Y,color='r',label='Gaussian')
    Plot.plot(X,Z,color='g',label='Lorentzian')
    Plot.plot(X,W,color='b',label='G+L')
    Plot.plot(X,V,color='k',label='G+L2')
    X = []
    Y = []
    Z = []
    W = []
    V = []
    for peak in peaks:
        X.append(4.0*math.pi*sind(peak[0]/2.0)/lam)
        s = 1.17741*math.sqrt(peak[4])*math.pi/18000.
        g = peak[6]*math.pi/18000.
        G = gamFW(g,s)
        H = gamFW2(g,s)
        Y.append(s/tand(peak[0]/2.))
        Z.append(g/tand(peak[0]/2.))
        W.append(G/tand(peak[0]/2.))
        V.append(H/tand(peak[0]/2.))
    Plot.plot(X,Y,'+',color='r',label='G peak')
    Plot.plot(X,Z,'+',color='g',label='L peak')
    Plot.plot(X,W,'+',color='b',label='G+L peak')
    Plot.plot(X,V,'+',color='k',label='G+L2 peak')
    Plot.legend(loc='best')
    Page.canvas.draw()
    
def PlotStrain(self,data):
    '''Plot 3D microstrain figure. In this instance data is for a phase
    '''
    PatternId = self.PatternId
    generalData = data['General']
    SGData = generalData['SGData']
    MuStrKeys = G2spc.MustrainNames(SGData)
    cell = generalData['Cell'][1:]
    A,B = G2lat.cell2AB(cell[:6])
    Vol = cell[6]
    useList = data['Histograms']
    numPlots = len(useList)
    
    try:
        plotNum = self.G2plotNB.plotList.index('Microstrain')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = mp3d.Axes3D(Page.figure)
    except ValueError:
        Plot = mp3d.Axes3D(self.G2plotNB.add3D('Microstrain'))
        plotNum = self.G2plotNB.plotList.index('Microstrain')
        Page = self.G2plotNB.nb.GetPage(plotNum)
    Page.SetFocus()
    self.G2plotNB.status.SetStatusText('Adjust frame size to get desired aspect ratio',1)
    
    for item in useList:
        if useList[item]['Show']:
            muStrain = useList[item]['Mustrain']
            PHI = np.linspace(0.,360.,30,True)
            PSI = np.linspace(0.,180.,30,True)
            X = np.outer(npsind(PHI),npsind(PSI))
            Y = np.outer(npcosd(PHI),npsind(PSI))
            Z = np.outer(np.ones(np.size(PHI)),npcosd(PSI))
            if muStrain[0] == 'isotropic':
                muiso = muStrain[1][0]*math.pi/0.018      #centidegrees to radians!
                X *= muiso
                Y *= muiso
                Z *= muiso                
                
            elif muStrain[0] == 'uniaxial':
                def uniaxMustrain(xyz,muiso,muaniso,axes):
                    cp = abs(np.dot(xyz,axes))
                    S = muiso+muaniso*cp
                    return S*xyz*math.pi/0.018      #centidegrees to radians!
                muiso,muaniso = muStrain[1][:2]
                axes = np.inner(A,np.array(muStrain[3]))
                axes /= nl.norm(axes)
                Shkl = np.array(muStrain[1])
                Shape = X.shape[0]
                XYZ = np.dstack((X,Y,Z))
                XYZ = np.nan_to_num(np.apply_along_axis(uniaxMustrain,2,XYZ,muiso,muaniso,axes))
                X,Y,Z = np.dsplit(XYZ,3)
                X = X[:,:,0]
                Y = Y[:,:,0]
                Z = Z[:,:,0]
                
            elif muStrain[0] == 'generalized':
                def genMustrain(xyz,SGData,A,Shkl):
                    uvw = np.inner(A.T,xyz)
                    Strm = np.array(G2spc.MustrainCoeff(uvw,SGData))
                    sum = np.sqrt(np.sum(np.multiply(Shkl,Strm)))*math.pi/0.018      #centidegrees to radians!
                    return sum*xyz
                Shkl = np.array(muStrain[4])
                if np.any(Shkl):
                    Shape = X.shape[0]
                    XYZ = np.dstack((X,Y,Z))
                    XYZ = np.nan_to_num(np.apply_along_axis(genMustrain,2,XYZ,SGData,A,Shkl))
                    X,Y,Z = np.dsplit(XYZ,3)
                    X = X[:,:,0]
                    Y = Y[:,:,0]
                    Z = Z[:,:,0]
                    
            if np.any(X) and np.any(Y) and np.any(Z):
                Plot.plot_surface(X,Y,Z,rstride=1,cstride=1,color='g')
                
            Plot.set_xlabel('X')
            Plot.set_ylabel('Y')
            Plot.set_zlabel('Z')
    Page.canvas.draw()
    
def PlotTexture(self,data,newPlot=False):
    '''Pole figure, inverse pole figure(?), 3D pole distribution and 3D inverse pole distribution(?)
    plotting; Need way to select  
    pole figure or pole distribution to be displayed - do in key enter menu
    dict generalData contains all phase info needed which is in data
    '''
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    PatternId = self.PatternId
    generalData = data['General']
    SGData = generalData['SGData']
    textureData = generalData['SH Texture']
    SHData = generalData['SH Texture']
    SHCoef = SHData['SH Coeff'][1]
    cell = generalData['Cell'][1:7]
    Start = True
               
    try:
        plotNum = self.G2plotNB.plotList.index('Texture')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
        if not Page.IsShown():
            Page.Show()
    except ValueError:
        Plot = self.G2plotNB.addMpl('Texture').gca()
        plotNum = self.G2plotNB.plotList.index('Texture')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.SetFocus()
    
    if 'Axial' in SHData['PlotType']:
        PH = np.array(SHData['PFhkl'])
        phi,beta = G2lat.CrsAng(PH,cell,SGData)
        ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
        X = np.linspace(0,90.0,26)
        Y = np.zeros_like(X)
        for i,a in enumerate(X):
            Y[i] = G2lat.polfcal(ODFln,SamSym[textureData['Model']],a,0.0)
        Plot.plot(X,Y,color='k',label=str(SHData['PFhkl']))
        Plot.legend(loc='best')
        Plot.set_title('Axial distribution for HKL='+str(SHData['PFhkl']))
        Plot.set_xlabel(r'$\psi$',fontsize=16)
        Plot.set_ylabel('MRD',fontsize=14)
        
        
    else:       
#        self.G2plotNB.status.SetStatusText('Adjust frame size to get desired aspect ratio',1)
        if 'inverse' in SHData['PlotType']:
            PX = np.array(SHData['PFxyz'])
            
        else:
            PH = np.array(SHData['PFhkl'])
            phi,beta = G2lat.CrsAng(PH,cell,SGData)
            npts = 201
            ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
            X,Y = np.meshgrid(np.linspace(1.,-1.,npts),np.linspace(-1.,1.,npts))
            R,P = np.sqrt(X**2+Y**2).flatten(),npatan2d(X,Y).flatten()
            R = np.where(R <= 1.,2.*npasind(R*0.70710678),0.0)
            Z = np.zeros_like(R)
            Z = G2lat.polfcal(ODFln,SamSym[textureData['Model']],R,P)
            Z = np.reshape(Z,(npts,npts))
            Img = Plot.imshow(Z.T,aspect='equal',cmap='binary')
            if newPlot:
                Page.figure.colorbar(Img)
                newPlot = False
            Plot.set_title('Pole figure for HKL='+str(SHData['PFhkl']))
    Page.canvas.draw()
    return

            
def PlotExposedImage(self,newPlot=False,event=None):
    '''General access module for 2D image plotting
    '''
    plotNo = self.G2plotNB.nb.GetSelection()
    if self.G2plotNB.nb.GetPageText(plotNo) == '2D Powder Image':
        PlotImage(self,newPlot,event,newImage=True)
    elif self.G2plotNB.nb.GetPageText(plotNo) == '2D Integration':
        PlotIntegration(self,newPlot,event)

def PlotImage(self,newPlot=False,event=None,newImage=True):
    '''Plot of 2D detector images as contoured plot. Also plot calibration ellipses,
    masks, etc.
    '''
    from matplotlib.patches import Ellipse,Arc,Circle,Polygon
    import numpy.ma as ma
    Dsp = lambda tth,wave: wave/(2.*sind(tth/2.))
    global Data,Masks
    colors=['b','g','r','c','m','k']
    Data = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
    Masks = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))

    def OnImMotion(event):
        Page.canvas.SetToolTipString('')
        sizexy = Data['size']
        if event.xdata and event.ydata:                 #avoid out of frame errors
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            item = self.itemPicked
            pixelSize = Data['pixelSize']
            scalex = 1000./pixelSize[0]
            scaley = 1000./pixelSize[1]
            if item and self.PatternTree.GetItemText(self.PickId) == 'Image Controls':
                if 'Text' in str(item):
                    Page.canvas.SetToolTipString('%8.3f %8.3fmm'%(event.xdata,event.ydata))
                else:
                    xcent,ycent = Data['center']
                    xpos = event.xdata-xcent
                    ypos = event.ydata-ycent
                    tth,azm = G2img.GetTthAzm(event.xdata,event.ydata,Data)
                    if 'line3' in  str(item) or 'line4' in str(item) and not Data['fullIntegrate']:
                        Page.canvas.SetToolTipString('%6d deg'%(azm))
                    elif 'line1' in  str(item) or 'line2' in str(item):
                        Page.canvas.SetToolTipString('%8.3fdeg'%(tth))                           
            else:
                xpos = event.xdata
                ypos = event.ydata
                xpix = xpos*scalex
                ypix = ypos*scaley
                Int = 0
                if (0 <= xpix <= sizexy[0]) and (0 <= ypix <= sizexy[1]):
                    Int = self.ImageZ[ypix][xpix]
                tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                Q = 2.*math.pi/dsp
                if self.setPoly:
                    self.G2plotNB.status.SetFields(['','Polygon mask pick - LB next point, RB close polygon'])
                else:
                    self.G2plotNB.status.SetFields(\
                        ['','Detector 2-th =%9.3fdeg, dsp =%9.3fA, Q = %6.5fA-1, azm = %7.2fdeg, I = %6d'%(tth,dsp,Q,azm,Int)])

    def OnImPlotKeyPress(event):
        try:
            PickName = self.PatternTree.GetItemText(self.PickId)
        except TypeError:
            return
        if PickName == 'Masks':
            Xpos = event.xdata
            if not Xpos:            #got point out of frame
                return
            Ypos = event.ydata
            if event.key == 's':
                Masks['Points'].append([Xpos,Ypos,1.])
            elif event.key == 'r':
                tth = G2img.GetTth(Xpos,Ypos,Data)
                Masks['Rings'].append([tth,0.1])
            elif event.key == 'a':
                tth,azm = G2img.GetTthAzm(Xpos,Ypos,Data)
                azm = int(azm)                
                Masks['Arcs'].append([tth,[azm-5,azm+5],0.1])
            elif event.key == 'p':
                self.setPoly = True
                Masks['Polygons'].append([])
                self.G2plotNB.status.SetFields(['','Polygon mask active - LB pick next point, RB close polygon'])
            G2imG.UpdateMasks(self,Masks)
        elif PickName == 'Image Controls':
            if event.key == 'c':
                Xpos = event.xdata
                if not Xpos:            #got point out of frame
                    return
                Ypos = event.ydata
                dlg = wx.MessageDialog(self,'Are you sure you want to change the center?',
                    'Center change',style=wx.OK|wx.CANCEL)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        print 'move center to: ',Xpos,Ypos
                        Data['center'] = [Xpos,Ypos]
                        G2imG.UpdateImageControls(self,Data,Masks)
                finally:
                    dlg.Destroy()
            elif event.key == 'l':
                if self.logPlot:
                    self.logPlot = False
                else:
                    self.logPlot = True
        PlotImage(self,newImage=True)
            
    def OnKeyBox(event):
        if self.G2plotNB.nb.GetSelection() == self.G2plotNB.plotList.index('2D Powder Image'):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            if event.key in 'l':
                OnImPlotKeyPress(event)
                        
    def OnImPick(event):
        if self.PatternTree.GetItemText(self.PickId) not in ['Image Controls','Masks']:
            return
        if self.setPoly:
            polygon = Masks['Polygons'][-1]
            xpos,ypos = event.mouseevent.xdata,event.mouseevent.ydata
            if xpos and ypos:                       #point inside image
                if len(polygon) > 2 and event.mouseevent.button == 3:
                    x0,y0 = polygon[0]
                    polygon.append([x0,y0])
                    self.setPoly = False
                    self.G2plotNB.status.SetFields(['','Polygon closed - RB drag a vertex to change shape'])
                else:
                    self.G2plotNB.status.SetFields(['','New polygon point: %.1f,%.1f'%(xpos,ypos)])
                    polygon.append([xpos,ypos])
                G2imG.UpdateMasks(self,Masks)
        else:
            if self.itemPicked is not None: return
            self.itemPicked = event.artist
            self.mousePicked = event.mouseevent
        
    def OnImRelease(event):
        try:
            PickName = self.PatternTree.GetItemText(self.PickId)
        except TypeError:
            return
        if PickName not in ['Image Controls','Masks']:
            return
        pixelSize = Data['pixelSize']
        scalex = 1000./pixelSize[0]
        scaley = 1000./pixelSize[1]
        pixLimit = Data['pixLimit']
        if self.itemPicked is None and PickName == 'Image Controls':
#            sizexy = Data['size']
            Xpos = event.xdata
            if not (Xpos and self.ifGetRing):                   #got point out of frame
                return
            Ypos = event.ydata
            if Ypos and not Page.toolbar._active:         #make sure zoom/pan not selected
                if event.button == 1:
                    Xpix = Xpos*scalex
                    Ypix = Ypos*scaley
                    xpos,ypos,I,J = G2img.ImageLocalMax(self.ImageZ,pixLimit,Xpix,Ypix)
                    if I and J:
                        xpos += .5                              #shift to pixel center
                        ypos += .5
                        xpos /= scalex                          #convert to mm
                        ypos /= scaley
                        Data['ring'].append([xpos,ypos])
                elif event.button == 3:
                    self.dataFrame.GetStatusBar().SetStatusText('Calibrating...')
                    if G2img.ImageCalibrate(self,Data):
                        self.dataFrame.GetStatusBar().SetStatusText('Calibration successful - Show ring picks to check')
                        print 'Calibration successful'
                    else:
                        self.dataFrame.GetStatusBar().SetStatusText('Calibration failed - Show ring picks to diagnose')
                        print 'Calibration failed'
                    self.ifGetRing = False
                    G2imG.UpdateImageControls(self,Data,Masks)
                    return
                PlotImage(self,newImage=False)
            return
        else:
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                ypos = event.ydata
                if self.ifGetRing:                          #delete a calibration ring pick
                    xypos = [xpos,ypos]
                    rings = Data['ring']
                    for ring in rings:
                        if np.allclose(ring,xypos,.01,0):
                            rings.remove(ring)                                                                       
                else:
                    tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                    itemPicked = str(self.itemPicked)
                    if 'Line2D' in itemPicked and PickName == 'Image Controls':
                        if 'line1' in itemPicked:
                            Data['IOtth'][0] = max(tth,0.001)
                        elif 'line2' in itemPicked:
                            Data['IOtth'][1] = tth
                        elif 'line3' in itemPicked:
                            Data['LRazimuth'][0] = int(azm)
                            if Data['fullIntegrate']:
                                Data['LRazimuth'][1] = Data['LRazimuth'][0]+360
                        elif 'line4' in itemPicked and not Data['fullIntegrate']:
                            Data['LRazimuth'][1] = int(azm)
                            
                        if Data['LRazimuth'][0] > Data['LRazimuth'][1]:
                            Data['LRazimuth'][0] -= 360
                            
                        azLim = np.array(Data['LRazimuth'])    
                        if np.any(azLim>360):
                            azLim -= 360
                            Data['LRazimuth'] = list(azLim)
                            
                        if  Data['IOtth'][0] > Data['IOtth'][1]:
                            Data['IOtth'][0],Data['IOtth'][1] = Data['IOtth'][1],Data['IOtth'][0]
                            
                        self.InnerTth.SetValue("%8.2f" % (Data['IOtth'][0]))
                        self.OuterTth.SetValue("%8.2f" % (Data['IOtth'][1]))
                        self.Lazim.SetValue("%6d" % (Data['LRazimuth'][0]))
                        self.Razim.SetValue("%6d" % (Data['LRazimuth'][1]))
                    elif 'Circle' in itemPicked and PickName == 'Masks':
                        spots = Masks['Points']
                        newPos = itemPicked.split(')')[0].split('(')[2].split(',')
                        newPos = np.array([float(newPos[0]),float(newPos[1])])
                        for spot in spots:
                            if np.allclose(np.array([spot[:2]]),newPos):
                                spot[:2] = xpos,ypos
                        G2imG.UpdateMasks(self,Masks)
                    elif 'Line2D' in itemPicked and PickName == 'Masks':
                        Obj = self.itemPicked.findobj()
                        rings = Masks['Rings']
                        arcs = Masks['Arcs']
                        polygons = Masks['Polygons']
                        for ring in self.ringList:
                            if Obj == ring[0]:
                                rN = ring[1]
                                if ring[2] == 'o':
                                    rings[rN][0] = G2img.GetTth(xpos,ypos,Data)-rings[rN][1]/2.
                                else:
                                    rings[rN][0] = G2img.GetTth(xpos,ypos,Data)+rings[rN][1]/2.
                        for arc in self.arcList:
                            if Obj == arc[0]:
                                aN = arc[1]
                                if arc[2] == 'o':
                                    arcs[aN][0] = G2img.GetTth(xpos,ypos,Data)-arcs[aN][2]/2
                                elif arc[2] == 'i':
                                    arcs[aN][0] = G2img.GetTth(xpos,ypos,Data)+arcs[aN][2]/2
                                elif arc[2] == 'l':
                                    arcs[aN][1][0] = int(G2img.GetAzm(xpos,ypos,Data))
                                else:
                                    arcs[aN][1][1] = int(G2img.GetAzm(xpos,ypos,Data))
                        for poly in self.polyList:
                            if Obj == poly[0]:
                                ind = self.itemPicked.contains(self.mousePicked)[1]['ind'][0]
                                oldPos = np.array([self.mousePicked.xdata,self.mousePicked.ydata])
                                pN = poly[1]
                                for i,xy in enumerate(polygons[pN]):
                                    if np.allclose(np.array([xy]),oldPos,atol=1.0):
                                        polygons[pN][i] = xpos,ypos
                        G2imG.UpdateMasks(self,Masks)
#                    else:                  #keep for future debugging
#                        print str(self.itemPicked),event.xdata,event.ydata,event.button
                PlotImage(self,newImage=True)
            self.itemPicked = None
            
    try:
        plotNum = self.G2plotNB.plotList.index('2D Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        if newImage:
            Page.figure.clf()
            Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = self.G2plotNB.addMpl('2D Powder Image').gca()
        plotNum = self.G2plotNB.plotList.index('2D Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnImPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnImMotion)
        Page.canvas.mpl_connect('pick_event', OnImPick)
        Page.canvas.mpl_connect('button_release_event', OnImRelease)
        xylim = []
    if not event:                       #event from GUI TextCtrl - don't want focus to change to plot!!!
        Page.SetFocus()
    Title = self.PatternTree.GetItemText(self.Image)[4:]
    self.G2plotNB.status.DestroyChildren()
    if self.logPlot:
        Title = 'log('+Title+')'
    Plot.set_title(Title)
    try:
        if self.PatternTree.GetItemText(self.PickId) in ['Image Controls',]:
            if self.logPlot:
                Choice = (' key press','l: log(I) off')
            else:
                Choice = (' key press','l: log(I) on')
            cb = wx.ComboBox(self.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
                choices=Choice)
            cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
            cb.SetValue(' key press')
    except TypeError:
        pass
    size,imagefile = self.PatternTree.GetItemPyData(self.Image)
    if imagefile != self.oldImagefile:
        imagefile = G2IO.CheckImageFile(self,imagefile)
        if not imagefile:
            self.G2plotNB.Delete('2D Powder Image')
            return
        self.PatternTree.SetItemPyData(self.Image,[size,imagefile])
        self.ImageZ = G2IO.GetImageData(self,imagefile,imageOnly=True)
#        print self.ImageZ.shape,self.ImageZ.size,Data['size'] #might be useful debugging line
        self.oldImagefile = imagefile

    imScale = 1
    if len(self.ImageZ) > 1024:
        imScale = len(self.ImageZ)/1024
    sizexy = Data['size']
    pixelSize = Data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    Xmax = sizexy[0]*pixelSize[0]/1000.
    Ymax = sizexy[1]*pixelSize[1]/1000.
    xlim = (0,Xmax)
    ylim = (Ymax,0)
    Imin,Imax = Data['range'][1]
    acolor = mpl.cm.get_cmap(Data['color'])
    xcent,ycent = Data['center']
    Plot.set_xlabel('Image x-axis, mm',fontsize=12)
    Plot.set_ylabel('Image y-axis, mm',fontsize=12)
    #do threshold mask - "real" mask - others are just bondaries
    Zlim = Masks['Thresholds'][1]
    wx.BeginBusyCursor()
    try:
            
        if newImage:                    
            MA = ma.masked_greater(ma.masked_less(self.ImageZ,Zlim[0]),Zlim[1])
            MaskA = ma.getmaskarray(MA)
            A = G2img.ImageCompress(MA,imScale)
            AM = G2img.ImageCompress(MaskA,imScale)
            if self.logPlot:
                A = np.log(A)
                AM = np.log(AM)
                Imin,Imax = [np.amin(A),np.amax(A)]
            ImgM = Plot.imshow(AM,aspect='equal',cmap='Reds',
                interpolation='nearest',vmin=0,vmax=2,extent=[0,Xmax,Xmax,0])
            Img = Plot.imshow(A,aspect='equal',cmap=acolor,
                interpolation='nearest',vmin=Imin,vmax=Imax,extent=[0,Xmax,Ymax,0])
            if self.setPoly:
                Img.set_picker(True)
    
        Plot.plot(xcent,ycent,'x')
        if Data['showLines']:
            LRAzim = Data['LRazimuth']                  #NB: integers
            Nazm = Data['outAzimuths']
            delAzm = float(LRAzim[1]-LRAzim[0])/Nazm
            AzmthOff = Data['azmthOff']
            IOtth = Data['IOtth']
            wave = Data['wavelength']
            dspI = wave/(2.0*sind(IOtth[0]/2.0))
            ellI = G2img.GetEllipse(dspI,Data)           #=False if dsp didn't yield an ellipse (ugh! a parabola or a hyperbola)
            dspO = wave/(2.0*sind(IOtth[1]/2.0))
            ellO = G2img.GetEllipse(dspO,Data)           #Ditto & more likely for outer ellipse
            Azm = np.array(range(LRAzim[0],LRAzim[1]+1))-AzmthOff
            if ellI:
                xyI = []
                for azm in Azm:
                    xyI.append(G2img.GetDetectorXY(dspI,azm-90.,Data))
                xyI = np.array(xyI)
                arcxI,arcyI = xyI.T
                Plot.plot(arcxI,arcyI,picker=3)
            if ellO:
                xyO = []
                for azm in Azm:
                    xyO.append(G2img.GetDetectorXY(dspO,azm-90.,Data))
                xyO = np.array(xyO)
                arcxO,arcyO = xyO.T
                Plot.plot(arcxO,arcyO,picker=3)
            if ellO and ellI:
                Plot.plot([arcxI[0],arcxO[0]],[arcyI[0],arcyO[0]],picker=3)
                Plot.plot([arcxI[-1],arcxO[-1]],[arcyI[-1],arcyO[-1]],picker=3)
            for i in range(Nazm):
                cake = LRAzim[0]+i*delAzm
                ind = np.searchsorted(Azm,cake)
                Plot.plot([arcxI[ind],arcxO[ind]],[arcyI[ind],arcyO[ind]],color='k',dashes=(5,5))
                    
        for xring,yring in Data['ring']:
            Plot.plot(xring,yring,'r+',picker=3)
        if Data['setRings']:
#            rings = np.concatenate((Data['rings']),axis=0)
            N = 0
            for ring in Data['rings']:
                xring,yring = np.array(ring).T[:2]
                Plot.plot(xring,yring,'+',color=colors[N%6])
                N += 1            
        for ellipse in Data['ellipses']:
            cent,phi,[width,height],col = ellipse
            Plot.add_artist(Ellipse([cent[0],cent[1]],2*width,2*height,phi,ec=col,fc='none'))
            Plot.text(cent[0],cent[1],'+',color=col,ha='center',va='center')
        #masks - mask lines numbered after integration limit lines
        spots = Masks['Points']
        rings = Masks['Rings']
        arcs = Masks['Arcs']
        polygons = Masks['Polygons']
        for x,y,d in spots:
            Plot.add_artist(Circle((x,y),radius=d/2,fc='none',ec='r',picker=3))
        self.ringList = []
        for iring,(tth,thick) in enumerate(rings):
            wave = Data['wavelength']
            x1,y1 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth+thick/2.,wave),Data))),2)
            x2,y2 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth-thick/2.,wave),Data))),2)
            self.ringList.append([Plot.plot(x1,y1,'r',picker=3),iring,'o'])            
            self.ringList.append([Plot.plot(x2,y2,'r',picker=3),iring,'i'])
        self.arcList = []
        for iarc,(tth,azm,thick) in enumerate(arcs):            
            wave = Data['wavelength']
            x1,y1 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth+thick/2.,wave),Data),azm)),2)
            x2,y2 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(max(0.01,tth-thick/2.),wave),Data),azm)),2)
            self.arcList.append([Plot.plot(x1,y1,'r',picker=3),iarc,'o'])            
            self.arcList.append([Plot.plot(x2,y2,'r',picker=3),iarc,'i'])
            self.arcList.append([Plot.plot([x1[0],x2[0]],[y1[0],y2[0]],'r',picker=3),iarc,'l'])
            self.arcList.append([Plot.plot([x1[-1],x2[-1]],[y1[-1],y2[-1]],'r',picker=3),iarc,'u'])
        self.polyList = []
        for ipoly,polygon in enumerate(polygons):
            x,y = np.hsplit(np.array(polygon),2)
            self.polyList.append([Plot.plot(x,y,'r+',picker=10),ipoly])
            Plot.plot(x,y,'r')            
        if newImage:
            colorBar = Page.figure.colorbar(Img)
        Plot.set_xlim(xlim)
        Plot.set_ylim(ylim)
        if not newPlot and xylim:
            Page.toolbar.push_current()
            Plot.set_xlim(xylim[0])
            Plot.set_ylim(xylim[1])
            xylim = []
            Page.toolbar.push_current()
            Page.toolbar.draw()
        else:
            Page.canvas.draw()
    finally:
        wx.EndBusyCursor()
        
def PlotIntegration(self,newPlot=False,event=None):
    '''Plot of 2D image after image integration with 2-theta and azimuth as coordinates
    '''
            
    def OnMotion(event):
        Page.canvas.SetToolTipString('')
        Page.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.ydata
        tth = event.xdata
        if azm and tth:
            self.G2plotNB.status.SetFields(\
                ['','Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm)])
                                
    try:
        plotNum = self.G2plotNB.plotList.index('2D Integration')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = self.G2plotNB.addMpl('2D Integration').gca()
        plotNum = self.G2plotNB.plotList.index('2D Integration')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.views = False
        view = False
    if not event:
        Page.SetFocus()
        
    Data = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
    image = self.Integrate[0]
    xsc = self.Integrate[1]
    ysc = self.Integrate[2]
    Imin,Imax = Data['range'][1]
    acolor = mpl.cm.get_cmap(Data['color'])
    Plot.set_title(self.PatternTree.GetItemText(self.Image)[4:])
    Plot.set_ylabel('azimuth',fontsize=12)
    Plot.set_xlabel('2-theta',fontsize=12)
    Img = Plot.imshow(image,cmap=acolor,vmin=Imin,vmax=Imax,interpolation='nearest', \
        extent=[ysc[0],ysc[-1],xsc[-1],xsc[0]],aspect='auto')
    colorBar = Page.figure.colorbar(Img)
    if Data['setRings'] and Data['rings']:
        rings = np.concatenate((Data['rings']),axis=0)
        for xring,yring,dsp in rings:
            x,y = G2img.GetTthAzm(xring,yring,Data)
            Plot.plot(x,y,'r+')
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2img.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2img.GetTthAzm(x,y,Data)
            azm = np.where(azm < 0.,azm+360,azm)
            Plot.plot(tth,azm,'b,')
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
                
def PlotTRImage(self,tax,tay,taz,newPlot=False):
    '''a test plot routine - not normally used
    ''' 
            
    def OnMotion(event):
        Page.canvas.SetToolTipString('')
        Page.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.xdata
        tth = event.ydata
        if azm and tth:
            self.G2plotNB.status.SetFields(\
                ['','Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm)])
                                
    try:
        plotNum = self.G2plotNB.plotList.index('2D Transformed Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = self.G2plotNB.addMpl('2D Transformed Powder Image').gca()
        plotNum = self.G2plotNB.plotList.index('2D Transformed Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.views = False
        view = False
    Page.SetFocus()
        
    Data = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
    Imin,Imax = Data['range'][1]
    step = (Imax-Imin)/5.
    V = np.arange(Imin,Imax,step)
    acolor = mpl.cm.get_cmap(Data['color'])
    Plot.set_title(self.PatternTree.GetItemText(self.Image)[4:])
    Plot.set_xlabel('azimuth',fontsize=12)
    Plot.set_ylabel('2-theta',fontsize=12)
    Plot.contour(tax,tay,taz,V,cmap=acolor)
    if Data['showLines']:
        IOtth = Data['IOtth']
        if Data['fullIntegrate']:
            LRAzim = [-180,180]
        else:
            LRAzim = Data['LRazimuth']                  #NB: integers
        Plot.plot([LRAzim[0],LRAzim[1]],[IOtth[0],IOtth[0]],picker=True)
        Plot.plot([LRAzim[0],LRAzim[1]],[IOtth[1],IOtth[1]],picker=True)
        if not Data['fullIntegrate']:
            Plot.plot([LRAzim[0],LRAzim[0]],[IOtth[0],IOtth[1]],picker=True)
            Plot.plot([LRAzim[1],LRAzim[1]],[IOtth[0],IOtth[1]],picker=True)
    if Data['setRings']:
        rings = np.concatenate((Data['rings']),axis=0)
        for xring,yring,dsp in rings:
            x,y = G2img.GetTthAzm(xring,yring,Data)
            Plot.plot(y,x,'r+')            
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2img.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2img.GetTthAzm(x,y,Data)
            Plot.plot(azm,tth,'b,')
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
        
def PlotStructure(self,data):
    '''Crystal structure plotting package. Can show structures as balls, sticks, lines,
    thermal motion ellipsoids and polyhedra
    '''
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(cell)         #Amat - crystal to cartesian, Bmat - inverse
    A4mat = np.concatenate((np.concatenate((Amat,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
    B4mat = np.concatenate((np.concatenate((Bmat,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
    Mydir = generalData['Mydir']
    atomData = data['Atoms']
    drawingData = data['Drawing']
    drawAtoms = drawingData['Atoms']
    cx,ct,cs = drawingData['atomPtrs']
    Wt = [255,255,255]
    Rd = [255,0,0]
    Gr = [0,255,0]
    Bl = [0,0,255]
    uBox = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    uEdges = np.array([
        [uBox[0],uBox[1]],[uBox[0],uBox[3]],[uBox[0],uBox[4]],[uBox[1],uBox[2]], 
        [uBox[2],uBox[3]],[uBox[1],uBox[5]],[uBox[2],uBox[6]],[uBox[3],uBox[7]], 
        [uBox[4],uBox[5]],[uBox[5],uBox[6]],[uBox[6],uBox[7]],[uBox[7],uBox[4]]])
    uColors = [Rd,Gr,Bl,Wt, Wt,Wt,Wt,Wt, Wt,Wt,Wt,Wt]
    altDown = False
    shiftDown = False
    ctrlDown = False
    
    def OnKeyBox(event):
        import Image
        Draw()                          #make sure plot is fresh!!
        mode = cb.GetValue()
        Fname = Mydir+'\\'+generalData['Name']+'.'+mode
        size = Page.canvas.GetSize()
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
        if mode in ['jpeg',]:
            Pix = glReadPixels(0,0,size[0],size[1],GL_RGBA, GL_UNSIGNED_BYTE)
            im = Image.new("RGBA", (size[0],size[1]))
        else:
            Pix = glReadPixels(0,0,size[0],size[1],GL_RGB, GL_UNSIGNED_BYTE)
            im = Image.new("RGB", (size[0],size[1]))
        im.fromstring(Pix)
        im.save(Fname,mode)
        cb.SetValue(' Save as:')
        self.G2plotNB.status.SetStatusText('Drawing saved to: '+Fname,1)
    
    def GetTruePosition(xy):
        View = glGetIntegerv(GL_VIEWPORT)
        Proj = glGetDoublev(GL_PROJECTION_MATRIX)
        Model = glGetDoublev(GL_MODELVIEW_MATRIX)
        Zmax = 1.
        for i,atom in enumerate(drawAtoms):
            x,y,z = atom[cx:cx+3]
            X,Y,Z = gluProject(x,y,z,Model,Proj,View)
            XY = [int(X),int(View[3]-Y)]
            if np.allclose(xy,XY,atol=10) and Z < Zmax:
                Zmax = Z
                SetSelectedAtoms(i)
                    
    def OnMouseDown(event):
        xy = event.GetPosition()
        if event.ShiftDown():
            GetTruePosition(xy)
        else:
            drawingData['Rotation'][3] = xy
        Draw()
        
    def OnMouseMove(event):
        newxy = event.GetPosition()
        page = getSelection()
        if event.ControlDown() and drawingData['showABC']:
            if event.LeftIsDown():
                SetTestRot(newxy)
            elif event.RightIsDown():
                SetTestPos(newxy)
            elif event.MiddleIsDown():
                SetTestRotZ(newxy)
            x,y,z = drawingData['testPos'][0]
            self.G2plotNB.status.SetStatusText('moving test point %.4f,%.4f,%.4f'%(x,y,z),1)
                
                
        if event.Dragging() and not event.ControlDown():
            if event.LeftIsDown():
                SetRotation(newxy)
                angX,angY,angZ = drawingData['Rotation'][:3]
                self.G2plotNB.status.SetStatusText('New rotation: %.2f, %.2f ,%.2f'%(angX,angY,angZ),1)
            elif event.RightIsDown():
                SetTranslation(newxy)
                Tx,Ty,Tz = drawingData['viewPoint'][0]
                self.G2plotNB.status.SetStatusText('New view point: %.4f, %.4f, %.4f'%(Tx,Ty,Tz),1)
            elif event.MiddleIsDown():
                SetRotationZ(newxy)
                angX,angY,angZ = drawingData['Rotation'][:3]
                self.G2plotNB.status.SetStatusText('New rotation: %.2f, %.2f, %.2f'%(angX,angY,angZ),1)
        Draw()
        
    def OnMouseWheel(event):
        drawingData['cameraPos'] += event.GetWheelRotation()/24
        drawingData['cameraPos'] = max(10,min(500,drawingData['cameraPos']))
        self.G2plotNB.status.SetStatusText('New camera distance: %.2f'%(drawingData['cameraPos']),1)
        page = getSelection()
        if page:
            if self.dataDisplay.GetPageText(page) == 'Draw Options':
                panel = self.dataDisplay.GetPage(page).GetChildren()[0].GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('cameraPos')].SetLabel('Camera Position: '+'%.2f'%(drawingData['cameraPos']))
                panel[names.index('cameraSlider')].SetValue(drawingData['cameraPos'])
        Draw()
        
    def getSelection():
        try:
            return self.dataDisplay.GetSelection()
        except AttributeError:
            print self.dataDisplay.GetLabel()
            self.G2plotNB.status.SetStatusText('Select this from Phase data window!')
            return 0
            
    def SetViewPointText(VP):
        page = getSelection()
        if page:
            if self.dataDisplay.GetPageText(page) == 'Draw Options':
                panel = self.dataDisplay.GetPage(page).GetChildren()[0].GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('viewPoint')].SetValue('%.3f, %.3f, %.3f'%(VP[0],VP[1],VP[2]))
            
    def ClearSelectedAtoms():
        page = getSelection()
        if page:
            if self.dataDisplay.GetPageText(page) == 'Draw Atoms':
                self.dataDisplay.GetPage(page).ClearSelection()      #this is the Atoms grid in Draw Atoms
            elif self.dataDisplay.GetPageText(page) == 'Atoms':
                self.dataDisplay.GetPage(page).ClearSelection()      #this is the Atoms grid in Atoms
                    
    def SetSelectedAtoms(ind):
        page = getSelection()
        if page:
            if self.dataDisplay.GetPageText(page) == 'Draw Atoms':
                self.dataDisplay.GetPage(page).SelectRow(ind)      #this is the Atoms grid in Draw Atoms
            elif self.dataDisplay.GetPageText(page) == 'Atoms':
                Id = drawAtoms[ind][-2]
                for i,atom in enumerate(atomData):
                    if atom[-1] == Id:
                        self.dataDisplay.GetPage(page).SelectRow(i)      #this is the Atoms grid in Atoms
                  
    def GetSelectedAtoms():
        page = getSelection()
        Ind = []
        if page:
            if self.dataDisplay.GetPageText(page) == 'Draw Atoms':
                Ind = self.dataDisplay.GetPage(page).GetSelectedRows()      #this is the Atoms grid in Draw Atoms
            elif self.dataDisplay.GetPageText(page) == 'Atoms':
                Ind = self.dataDisplay.GetPage(page).GetSelectedRows()      #this is the Atoms grid in Atoms
        return Ind
                                       
    def OnKey(event):           #on key UP!!
        keyCode = event.GetKeyCode()
        if keyCode > 255:
            keyCode = 0
        key,xyz = chr(keyCode),event.GetPosition()
        indx = drawingData['selectedAtoms']
        if key in ['c','C']:
            drawingData['viewPoint'] = [[.5,.5,.5],[0,0]]
            drawingData['testPos'] = [[-.1,-.1,-.1],[0.0,0.0,0.0],[0,0]]
            drawingData['Rotation'] = [0.0,0.0,0.0,[]]
            SetViewPointText(drawingData['viewPoint'][0])
        elif key in ['n','N']:
            drawAtoms = drawingData['Atoms']
            pI = drawingData['viewPoint'][1]
            if indx:
                pI[0] = indx[pI[1]]
                Tx,Ty,Tz = drawAtoms[pI[0]][cx:cx+3]
                pI[1] += 1
                if pI[1] >= len(indx):
                    pI[1] = 0
            else:
                Tx,Ty,Tz = drawAtoms[pI[0]][cx:cx+3]                
                pI[0] += 1
                if pI[0] >= len(drawAtoms):
                    pI[0] = 0
            drawingData['viewPoint'] = [[Tx,Ty,Tz],pI]
            SetViewPointText(drawingData['viewPoint'][0])
                
        elif key in ['p','P']:
            drawAtoms = drawingData['Atoms']
            pI = drawingData['viewPoint'][1]
            if indx:
                pI[0] = indx[pI[1]]
                Tx,Ty,Tz = drawAtoms[pI[0]][cx:cx+3]
                pI[1] -= 1
                if pI[1] < 0:
                    pI[1] = len(indx)-1
            else:
                Tx,Ty,Tz = drawAtoms[pI[0]][cx:cx+3]                
                pI[0] -= 1
                if pI[0] < 0:
                    pI[0] = len(drawAtoms)-1
            drawingData['viewPoint'] = [[Tx,Ty,Tz],pI]
            SetViewPointText(drawingData['viewPoint'][0])            
        Draw()
            
    def SetBackground():
        R,G,B,A = Page.camera['backColor']
        glClearColor(R,G,B,A)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
    def SetLights():
        glEnable(GL_DEPTH_TEST)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0)
        glLightfv(GL_LIGHT0,GL_AMBIENT,[1,1,1,.8])
        glLightfv(GL_LIGHT0,GL_DIFFUSE,[1,1,1,1])
        
    def SetTranslation(newxy):
        Tx,Ty,Tz = drawingData['viewPoint'][0]
        anglex,angley,anglez,oldxy = drawingData['Rotation']
        Rx = G2lat.rotdMat(anglex,0)
        Ry = G2lat.rotdMat(angley,1)
        Rz = G2lat.rotdMat(anglez,2)
        dxy = list(newxy-oldxy)+[0,]
        dxy = np.inner(Bmat,np.inner(Rz,np.inner(Ry,np.inner(Rx,dxy))))
        Tx -= dxy[0]*0.01
        Ty += dxy[1]*0.01
        Tz -= dxy[2]*0.01
        drawingData['Rotation'][3] = newxy
        drawingData['viewPoint'][0] =  Tx,Ty,Tz
        SetViewPointText([Tx,Ty,Tz])
        
    def SetTestPos(newxy):
        Tx,Ty,Tz = drawingData['testPos'][0]
        anglex,angley,anglez,oldxy = drawingData['Rotation']
        Rx = G2lat.rotdMat(anglex,0)
        Ry = G2lat.rotdMat(angley,1)
        Rz = G2lat.rotdMat(anglez,2)
        dxy = list(newxy-oldxy)+[0,]
        dxy = np.inner(Rz,np.inner(Ry,np.inner(Rx,dxy)))
        Tx += dxy[0]*0.001
        Ty -= dxy[1]*0.001
        Tz += dxy[2]*0.001
        drawingData['Rotation'][3] = newxy
        drawingData['testPos'][0] =  Tx,Ty,Tz
        
    def SetTestRot(newxy):
        Txyz = np.array(drawingData['testPos'][0])
        oldxy = drawingData['testPos'][2]
        Ax,Ay,Az = drawingData['testPos'][1]
        Vxyz = np.array(drawingData['viewPoint'][0])
        Dxyz = np.inner(Amat,Txyz-Vxyz)
        dxy = list(newxy-oldxy)+[0,]
        Ax += dxy[1]*0.01
        Ay += dxy[0]*0.01
        Rx = G2lat.rotdMat(Ax,0)
        Ry = G2lat.rotdMat(Ay,1)
        Dxyz = np.inner(Ry,np.inner(Rx,Dxyz))        
        Dxyz = np.inner(Bmat,Dxyz)+Vxyz
        drawingData['testPos'][1] = [Ax,Ay,Az]
        drawingData['testPos'][2] = newxy
        drawingData['testPos'][0] = Dxyz
        
    def SetTestRotZ(newxy):
        Txyz = np.array(drawingData['testPos'][0])
        oldxy = drawingData['testPos'][2]
        Ax,Ay,Az = drawingData['testPos'][1]
        Vxyz = np.array(drawingData['viewPoint'][0])
        Dxyz = np.inner(Amat,Txyz-Vxyz)        
        dxy = list(newxy-oldxy)+[0,]
        Az += (dxy[0]+dxy[1])*.01
        Rz = G2lat.rotdMat(Az,2)
        Dxyz = np.inner(Rz,Dxyz)       
        Dxyz = np.inner(Bmat,Dxyz)+Vxyz
        drawingData['testPos'][1] = [Ax,Ay,Az]
        drawingData['testPos'][2] = newxy
        drawingData['testPos'][0] = Dxyz
                              
    def SetRotation(newxy):        
        anglex,angley,anglez,oldxy = drawingData['Rotation']
        dxy = newxy-oldxy
        anglex += dxy[1]*.25
        angley += dxy[0]*.25
        oldxy = newxy
        drawingData['Rotation'] = [anglex,angley,anglez,oldxy]
        
    def SetRotationZ(newxy):                        
        anglex,angley,anglez,oldxy = drawingData['Rotation']
        dxy = newxy-oldxy
        anglez += (dxy[0]+dxy[1])*.25
        oldxy = newxy
        drawingData['Rotation'] = [anglex,angley,anglez,oldxy]
        
    def RenderBox():
        glEnable(GL_COLOR_MATERIAL)
        glLineWidth(2)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_LINE_SMOOTH)
        glBegin(GL_LINES)
        for line,color in zip(uEdges,uColors):
            glColor3ubv(color)
            glVertex3fv(line[0])
            glVertex3fv(line[1])
        glEnd()
        glColor4ubv([0,0,0,0])
        glDisable(GL_LINE_SMOOTH)
        glDisable(GL_BLEND)
        glDisable(GL_COLOR_MATERIAL)
        
    def RenderUnitVectors(x,y,z):
        xyz = np.array([x,y,z])
        glEnable(GL_COLOR_MATERIAL)
        glLineWidth(1)
        glPushMatrix()
        glTranslate(x,y,z)
        glScalef(1/cell[0],1/cell[1],1/cell[2])
        glBegin(GL_LINES)
        for line,color in zip(uEdges,uColors)[:3]:
            glColor3ubv(color)
            glVertex3fv(line[0])
            glVertex3fv(line[1])
        glEnd()
        glPopMatrix()
        glColor4ubv([0,0,0,0])
        glDisable(GL_COLOR_MATERIAL)
                
    def RenderSphere(x,y,z,radius,color):
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glPushMatrix()
        glTranslate(x,y,z)
        glMultMatrixf(B4mat.T)
        q = gluNewQuadric()
        gluSphere(q,radius,20,10)
        glPopMatrix()
        
    def RenderEllipsoid(x,y,z,ellipseProb,E,R4,color):
        s1,s2,s3 = E
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glPushMatrix()
        glTranslate(x,y,z)
        glMultMatrixf(B4mat.T)
        glMultMatrixf(R4.T)
        glEnable(GL_NORMALIZE)
        glScale(s1,s2,s3)
        q = gluNewQuadric()
        gluSphere(q,ellipseProb,20,10)
        glDisable(GL_NORMALIZE)
        glPopMatrix()
        
    def RenderBonds(x,y,z,Bonds,radius,color,slice=20):
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glPushMatrix()
        glTranslate(x,y,z)
        glMultMatrixf(B4mat.T)
        for bond in Bonds:
            glPushMatrix()
            Dx = np.inner(Amat,bond)
            Z = np.sqrt(np.sum(Dx**2))
            azm = atan2d(-Dx[1],-Dx[0])
            phi = acosd(Dx[2]/Z)
            glRotate(-azm,0,0,1)
            glRotate(phi,1,0,0)
            q = gluNewQuadric()
            gluCylinder(q,radius,radius,Z,slice,2)
            glPopMatrix()            
        glPopMatrix()
                
    def RenderLines(x,y,z,Bonds,color):
        xyz = np.array([x,y,z])
        glEnable(GL_COLOR_MATERIAL)
        glLineWidth(1)
        glColor3fv(color)
        glPushMatrix()
        glBegin(GL_LINES)
        for bond in Bonds:
            glVertex3fv(xyz)
            glVertex3fv(xyz+bond)
        glEnd()
        glColor4ubv([0,0,0,0])
        glPopMatrix()
        glDisable(GL_COLOR_MATERIAL)
        
    def RenderPolyhedra(x,y,z,Faces,color):
        glPushMatrix()
        glTranslate(x,y,z)
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glShadeModel(GL_SMOOTH)
        glMultMatrixf(B4mat.T)
        for face,norm in Faces:
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
            glFrontFace(GL_CW)
            glNormal3fv(norm)
            glBegin(GL_TRIANGLES)
            for vert in face:
                glVertex3fv(vert)
            glEnd()
        glPopMatrix()
        
    def RenderBackbone(Backbone,BackboneColor,radius):
        glPushMatrix()
        glMultMatrixf(B4mat.T)
        glEnable(GL_COLOR_MATERIAL)
        glShadeModel(GL_SMOOTH)
        gleSetJoinStyle(TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP)
        glePolyCylinder(Backbone,BackboneColor,radius)
        glPopMatrix()        
        glDisable(GL_COLOR_MATERIAL)
        
    def RenderLabel(x,y,z,label,r):       
        glPushMatrix()
        glTranslate(x,y,z)
        glMultMatrixf(B4mat.T)
        glDisable(GL_LIGHTING)
        glColor3f(0,1.,0)
        glRasterPos3f(r,r,r)
        for c in list(label):
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13,ord(c))
        glEnable(GL_LIGHTING)
        glPopMatrix()
                            
    def Draw():
        import numpy.linalg as nl
        Ind = GetSelectedAtoms()
        VS = np.array(Page.canvas.GetSize())
        aspect = float(VS[0])/float(VS[1])
        cPos = drawingData['cameraPos']
        Zclip = drawingData['Zclip']*cPos/200.
        anglex,angley,anglez = drawingData['Rotation'][:3]
        Tx,Ty,Tz = drawingData['viewPoint'][0]
        cx,ct,cs = drawingData['atomPtrs']
        bondR = drawingData['bondRadius']
        G,g = G2lat.cell2Gmat(cell)
        GS = G
        GS[0][1] = GS[1][0] = math.sqrt(GS[0][0]*GS[1][1])
        GS[0][2] = GS[2][0] = math.sqrt(GS[0][0]*GS[2][2])
        GS[1][2] = GS[2][1] = math.sqrt(GS[1][1]*GS[2][2])
        ellipseProb = G2lat.criticalEllipse(drawingData['ellipseProb']/100.)
        
        SetBackground()
        glInitNames()
        glPushName(0)
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0,0,VS[0],VS[1])
        gluPerspective(20.,aspect,cPos-Zclip,cPos+Zclip)
        gluLookAt(0,0,cPos,0,0,0,0,1,0)
        SetLights()            
            
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotate(anglez,0,0,1)
        glRotate(anglex,cosd(anglez),-sind(anglez),0)
        glRotate(angley,sind(anglez),cosd(anglez),0)
        glMultMatrixf(A4mat.T)
        glTranslate(-Tx,-Ty,-Tz)
        if drawingData['unitCellBox']:
            RenderBox()
        if drawingData['showABC']:
            x,y,z = drawingData['testPos'][0]
#            if altDown:
#                self.G2plotNB.status.SetStatusText('moving test point %.4f,%.4f,%.4f'%(x,y,z),1)
#            else:
#                self.G2plotNB.status.SetStatusText('test point %.4f,%.4f,%.4f'%(x,y,z),1)            
            RenderUnitVectors(x,y,z)
        Backbone = []
        BackboneColor = []
        time0 = time.time()
        for iat,atom in enumerate(drawingData['Atoms']):
            x,y,z = atom[cx:cx+3]
            Bonds = atom[-2]
            Faces = atom[-1]
            try:
                atNum = generalData['AtomTypes'].index(atom[ct])
            except ValueError:
                atNum = -1
            CL = atom[cs+2]
            color = np.array(CL)/255.
            if iat in Ind:
                color = np.array(Gr)/255.
            radius = 0.5
            if atom[cs] != '':
                glLoadName(atom[-3])                    
            if 'balls' in atom[cs]:
                vdwScale = drawingData['vdwScale']
                ballScale = drawingData['ballScale']
                if atNum < 0:
                    radius = 0.2
                elif 'H' == atom[ct]:
                    if drawingData['showHydrogen']:
                        if 'vdW' in atom[cs] and atNum >= 0:
                            radius = vdwScale*generalData['vdWRadii'][atNum]
                        else:
                            radius = ballScale*drawingData['sizeH']
                    else:
                        radius = 0.0
                else:
                    if 'vdW' in atom[cs]:
                        radius = vdwScale*generalData['vdWRadii'][atNum]
                    else:
                        radius = ballScale*generalData['BondRadii'][atNum]
                RenderSphere(x,y,z,radius,color)
                if 'sticks' in atom[cs]:
                    RenderBonds(x,y,z,Bonds,bondR,color)
            elif 'ellipsoids' in atom[cs]:
                RenderBonds(x,y,z,Bonds,bondR,color)
                if atom[cs+3] == 'A':                    
                    Uij = atom[cs+5:cs+11]
                    U = np.multiply(G2spc.Uij2U(Uij),GS)
                    U = np.inner(Amat,np.inner(U,Amat).T)
                    E,R = nl.eigh(U)
                    R4 = np.concatenate((np.concatenate((R,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
                    E = np.sqrt(E)
                    if atom[ct] == 'H' and not drawingData['showHydrogen']:
                        pass
                    else:
                        RenderEllipsoid(x,y,z,ellipseProb,E,R4,color)                    
                else:
                    if atom[ct] == 'H' and not drawingData['showHydrogen']:
                        pass
                    else:
                        radius = ellipseProb*math.sqrt(abs(atom[cs+4]))
                        RenderSphere(x,y,z,radius,color)
            elif 'lines' in atom[cs]:
                radius = 0.1
                RenderLines(x,y,z,Bonds,color)
#                RenderBonds(x,y,z,Bonds,0.05,color,6)
            elif atom[cs] == 'sticks':
                radius = 0.1
                RenderBonds(x,y,z,Bonds,bondR,color)
            elif atom[cs] == 'polyhedra':
                RenderPolyhedra(x,y,z,Faces,color)
            elif atom[cs] == 'backbone':
                if atom[ct-1].split()[0] in ['C','N']:
                    Backbone.append(list(np.inner(Amat,np.array([x,y,z]))))
                    BackboneColor.append(list(color))
                    
            if atom[cs+1] == 'type':
                RenderLabel(x,y,z,atom[ct],radius)
            elif atom[cs+1] == 'name':
                RenderLabel(x,y,z,atom[ct-1],radius)
            elif atom[cs+1] == 'number':
                RenderLabel(x,y,z,str(iat+1),radius)
            elif atom[cs+1] == 'residue' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-4],radius)
            elif atom[cs+1] == '1-letter' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-3],radius)
            elif atom[cs+1] == 'chain' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-2],radius)
        if Backbone:
            RenderBackbone(Backbone,BackboneColor,bondR)
#        print time.time()-time0
        Page.canvas.SwapBuffers()
       
    def OnSize(event):
        Draw()
        
    def OnFocus(event):
        Draw()
        
    try:
        plotNum = self.G2plotNB.plotList.index(generalData['Name'])
        Page = self.G2plotNB.nb.GetPage(plotNum)        
    except ValueError:
        Plot = self.G2plotNB.addOgl(generalData['Name'])
        plotNum = self.G2plotNB.plotList.index(generalData['Name'])
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.views = False
        view = False
        altDown = False
    Page.SetFocus()
    cb = wx.ComboBox(self.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=(' save as:','jpeg','tiff','bmp'))
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' save as:')
    Page.canvas.Bind(wx.EVT_MOUSEWHEEL, OnMouseWheel)
    Page.canvas.Bind(wx.EVT_LEFT_DOWN, OnMouseDown)
    Page.canvas.Bind(wx.EVT_RIGHT_DOWN, OnMouseDown)
    Page.canvas.Bind(wx.EVT_MIDDLE_DOWN, OnMouseDown)
    Page.canvas.Bind(wx.EVT_KEY_UP, OnKey)
    Page.canvas.Bind(wx.EVT_MOTION, OnMouseMove)
    Page.canvas.Bind(wx.EVT_SIZE, OnSize)
    Page.canvas.Bind(wx.EVT_SET_FOCUS, OnFocus)
    Page.camera['position'] = drawingData['cameraPos']
    Page.camera['viewPoint'] = np.inner(Amat,drawingData['viewPoint'][0])
    Page.camera['backColor'] = np.array(list(drawingData['backColor'])+[0,])/255.
    Page.canvas.SetCurrent()
    Draw()
        