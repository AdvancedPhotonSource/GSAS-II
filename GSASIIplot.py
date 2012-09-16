# -*- coding: utf-8 -*-
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
import numpy.ma as ma
import numpy.linalg as nl
import wx
import wx.aui
import wx.glcanvas
import matplotlib as mpl
import mpl_toolkits.mplot3d.axes3d as mp3d
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIgrid as G2gd
import GSASIIimage as G2img
import GSASIIpwd as G2pwd
import GSASIIIO as G2IO
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImath as G2mth
import pytexture as ptx
from  OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GLE import *
from matplotlib.backends.backend_wx import _load_bitmap
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
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda x,y: 180.*np.arctan2(x,y)/np.pi
    
class G2PlotMpl(wx.Panel):    
    def __init__(self,parent,id=-1,dpi=None,**kwargs):
        wx.Panel.__init__(self,parent,id=id,**kwargs)
        mpl.rcParams['legend.fontsize'] = 10
        self.figure = mpl.figure.Figure(dpi=dpi,figsize=(5,6))
        self.canvas = Canvas(self,-1,self.figure)
        self.toolbar = GSASIItoolbar(self.canvas)

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
        self.toolbar = GSASIItoolbar(self.canvas)

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
        
class GSASIItoolbar(Toolbar):
    ON_MPL_HELP = wx.NewId()
    def __init__(self,plotCanvas):
        Toolbar.__init__(self,plotCanvas)
        POSITION_OF_CONFIGURE_SUBPLOTS_BTN = 6
        self.DeleteToolByPos(POSITION_OF_CONFIGURE_SUBPLOTS_BTN)
        help = os.path.join(os.path.split(__file__)[0],'help.ico')
        self.AddSimpleTool(self.ON_MPL_HELP,_load_bitmap(help),'Help on','Show help on')
        wx.EVT_TOOL(self,self.ON_MPL_HELP,self.OnHelp)
    def OnHelp(self,event):
        Page = self.GetParent().GetParent()
        pageNo = Page.GetSelection()
        bookmark = Page.GetPageText(pageNo)
        bookmark = bookmark.strip(')').replace('(','_')
        G2gd.ShowHelp(bookmark,self.TopLevelParent)

################################################################################
##### PlotSngl
################################################################################
            
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
    HKLref = self.PatternTree.GetItemPyData(self.Sngl)[1]
    Data = self.PatternTree.GetItemPyData( 
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
    for refl in HKLref:
        H = np.array(refl[:3])
        sig,Fosq,Fcsq = refl[7:10]
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
       
################################################################################
##### PlotPatterns
################################################################################
            
def PlotPatterns(G2frame,newPlot=False):
    '''Powder pattern plotting package - displays single or multiple powder patterns as intensity vs
    2-theta or q (future TOF). Can display multiple patterns as "waterfall plots" or contour plots. Log I 
    plotting available.
    '''
    global HKL
    
    def OnKeyBox(event):
        if G2frame.G2plotNB.nb.GetSelection() == G2frame.G2plotNB.plotList.index('Powder Patterns'):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            wx.CallAfter(OnPlotKeyPress,event)
                        
    def OnPlotKeyPress(event):
        newPlot = False
        if event.key == 'w':
            if G2frame.Weight:
                G2frame.Weight = False
            else:
                G2frame.Weight = True
                G2frame.SinglePlot = True
            newPlot = True
        elif event.key == 'n':
            if G2frame.Contour:
                pass
            else:
                if G2frame.logPlot:
                    G2frame.logPlot = False
                else:
                    G2frame.Offset[0] = 0
                    G2frame.logPlot = True
                newPlot = True
        elif event.key == 'u':
            if G2frame.Contour:
                G2frame.Cmax = min(1.0,G2frame.Cmax*1.2)
            elif G2frame.logPlot:
                pass
            elif G2frame.Offset[0] < 100.:
                G2frame.Offset[0] += 1.
        elif event.key == 'd':
            if G2frame.Contour:
                G2frame.Cmax = max(0.0,G2frame.Cmax*0.8)
            elif G2frame.logPlot:
                pass
            elif G2frame.Offset[0] > 0.:
                G2frame.Offset[0] -= 1.
        elif event.key == 'l':
            G2frame.Offset[1] -= 1.
        elif event.key == 'r':
            G2frame.Offset[1] += 1.
        elif event.key == 'o':
            G2frame.Offset = [0,0]
        elif event.key == 'c':
            newPlot = True
            if G2frame.Contour:
                G2frame.Contour = False
            else:
                G2frame.Contour = True
                G2frame.SinglePlot = False
                G2frame.Offset = [0.,0.]
        elif event.key == 'q':
            newPlot = True
            if G2frame.qPlot:
                G2frame.qPlot = False
            else:
                G2frame.qPlot = True
        elif event.key == 's':
            if G2frame.Contour:
                choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
                choice.sort()
                dlg = wx.SingleChoiceDialog(G2frame,'Select','Color scheme',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    G2frame.ContourColor = choice[sel]
                else:
                    G2frame.ContourColor = 'Paired'
                dlg.Destroy()
            else:                
                if G2frame.SinglePlot:
                    G2frame.SinglePlot = False
                else:
                    G2frame.SinglePlot = True
            newPlot = True
        elif event.key == '+':
            if G2frame.PickId:
                G2frame.PickId = False
        elif event.key == 'i':                  #for smoothing contour plot
            choice = ['nearest','bilinear','bicubic','spline16','spline36','hanning',
               'hamming','hermite','kaiser','quadric','catrom','gaussian','bessel',
               'mitchell','sinc','lanczos']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Interpolation',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.Interpolate = choice[sel]
            else:
                G2frame.Interpolate = 'nearest'
            dlg.Destroy()
            
        PlotPatterns(G2frame,newPlot=newPlot)
        
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                Values,Names = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[1::2]
                Parms = dict(zip(Names,Values))
                try:
                    wave = Parms['Lam']
                except KeyError:
                    wave = Parms['Lam1']
                if G2frame.qPlot:
                    try:
                        xpos = 2.0*asind(xpos*wave/(4*math.pi))
                    except ValueError:      #avoid bad value in asin beyond upper limit
                        pass
                dsp = 0.0
                if abs(xpos) > 0.:                  #avoid possible singularity at beam center
                    dsp = wave/(2.*sind(abs(xpos)/2.0))
                if G2frame.Contour:
                    G2frame.G2plotNB.status.SetStatusText('2-theta =%9.3f d =%9.5f pattern ID =%5d'%(xpos,dsp,int(ypos)),1)
                else:
                    G2frame.G2plotNB.status.SetStatusText('2-theta =%9.3f d =%9.5f Intensity =%9.1f'%(xpos,dsp,ypos),1)
                if G2frame.itemPicked:
                    Page.canvas.SetToolTipString('%9.3f'%(xpos))
                if G2frame.PickId:
                    found = []
                    if G2frame.PatternTree.GetItemText(G2frame.PickId) in ['Index Peak List','Unit Cells List','Reflection Lists'] or \
                        'PWDR' in G2frame.PatternTree.GetItemText(PickId):
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
                G2frame.G2plotNB.status.SetStatusText('Select PWDR powder pattern first',1)
                                                   
    def OnPick(event):
        if G2frame.itemPicked is not None: return
        PatternId = G2frame.PatternId
        try:
            Values,Names = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[1::2]
        except TypeError:
            return
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        PickId = G2frame.PickId
        pick = event.artist
        mouse = event.mouseevent       
        xpos = pick.get_xdata()
        ypos = pick.get_ydata()
        ind = event.ind
        xy = list(zip(np.take(xpos,ind),np.take(ypos,ind))[0])
        if G2frame.PatternTree.GetItemText(PickId) == 'Peak List':
            if ind.all() != [0]:                                    #picked a data point
                if 'C' in Parms['Type']:                            #CW data - TOF later in an elif
                    ins = [Parms[x] for x in ['U','V','W','X','Y']]
                    if G2frame.qPlot:                              #qplot - convert back to 2-theta
                        xy[0] = 2.0*asind(xy[0]*wave/(4*math.pi))
                    sig = ins[0]*tand(xy[0]/2.0)**2+ins[1]*tand(xy[0]/2.0)+ins[2]
                    gam = ins[3]/cosd(xy[0]/2.0)+ins[4]*tand(xy[0]/2.0)           
                    data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
                    XY = [xy[0],0, xy[1],1, sig,0, gam,0]       #default refine intensity 1st
                data.append(XY)
                G2pdG.UpdatePeakGrid(G2frame,data)
                PlotPatterns(G2frame)
            else:                                                   #picked a peak list line
                G2frame.itemPicked = pick
        elif G2frame.PatternTree.GetItemText(PickId) == 'Limits':
            if ind.all() != [0]:                                    #picked a data point
                LimitId = G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits')
                data = G2frame.PatternTree.GetItemPyData(LimitId)
                if 'C' in Parms['Type']:                            #CW data - TOF later in an elif
                    if G2frame.qPlot:                              #qplot - convert back to 2-theta
                        xy[0] = 2.0*asind(xy[0]*wave/(4*math.pi))
                if mouse.button==1:
                    data[1][0] = min(xy[0],data[1][1])
                if mouse.button==3:
                    data[1][1] = max(xy[0],data[1][0])
                G2frame.PatternTree.SetItemPyData(LimitId,data)
                G2pdG.UpdateLimitsGrid(G2frame,data)
                PlotPatterns(G2frame)
            else:                                                   #picked a limit line
                G2frame.itemPicked = pick
        elif G2frame.PatternTree.GetItemText(PickId) == 'Reflection Lists' or \
            'PWDR' in G2frame.PatternTree.GetItemText(PickId):
            G2frame.itemPicked = pick
            pick = str(pick)
        
    def OnRelease(event):
        if G2frame.itemPicked is None: return
        Values,Names = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[1::2]
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        xpos = event.xdata
        PickId = G2frame.PickId
        if G2frame.PatternTree.GetItemText(PickId) in ['Peak List','Limits'] and xpos:
            lines = []
            for line in G2frame.Lines: 
                lines.append(line.get_xdata()[0])
#            print G2frame.itemPicked.get_xdata()
            lineNo = lines.index(G2frame.itemPicked.get_xdata()[0])
            if  lineNo in [0,1]:
                LimitId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits')
                data = G2frame.PatternTree.GetItemPyData(LimitId)
                if G2frame.qPlot:
                    data[1][lineNo] = 2.0*asind(wave*xpos/(4*math.pi))
                else:
                    data[1][lineNo] = xpos
                G2frame.PatternTree.SetItemPyData(LimitId,data)
                if G2frame.PatternTree.GetItemText(G2frame.PickId) == 'Limits':
                    G2pdG.UpdateLimitsGrid(G2frame,data)
            else:
                PeakId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List')
                data = G2frame.PatternTree.GetItemPyData(PeakId)
#                print 'peaks',xpos
                if event.button == 3:
                    del data[lineNo-2]
                else:
                    if G2frame.qPlot:
                        data[lineNo-2][0] = 2.0*asind(wave*xpos/(4*math.pi))
                    else:
                        data[lineNo-2][0] = xpos
                G2frame.PatternTree.SetItemPyData(PeakId,data)
                G2pdG.UpdatePeakGrid(G2frame,data)
        elif (G2frame.PatternTree.GetItemText(PickId) == 'Reflection Lists' or \
            'PWDR' in G2frame.PatternTree.GetItemText(PickId)) and xpos:
            Phases = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            pick = str(G2frame.itemPicked).split('(')[1].strip(')')
            if 'line' not in pick:       #avoid data points, etc.
                num = Phases.keys().index(pick)
                if num:
                    G2frame.refDelt = -(event.ydata-G2frame.refOffset)/(num*Ymax)
                else:       #1st row of refl ticks
                    G2frame.refOffset = event.ydata
        PlotPatterns(G2frame)
        G2frame.itemPicked = None    

    xylim = []
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Powder Patterns')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
    except ValueError:
        newPlot = True
        G2frame.Cmax = 1.0
        Plot = G2frame.G2plotNB.addMpl('Powder Patterns').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Powder Patterns')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.canvas.mpl_connect('pick_event', OnPick)
        Page.canvas.mpl_connect('button_release_event', OnRelease)
    Page.SetFocus()
    G2frame.G2plotNB.status.DestroyChildren()
    if G2frame.Contour:
        Choice = (' key press','d: lower contour max','u: raise contour max',
            'i: interpolation method','s: color scheme','c: contour off')
    else:
        if G2frame.logPlot:
            Choice = (' key press','n: log(I) off','l: offset left','r: offset right',
                'c: contour on','q: toggle q plot','s: toggle single plot','+: no selection')
        else:
            Choice = (' key press','l: offset left','r: offset right','d: offset down',
                'u: offset up','o: reset offset','n: log(I) on','c: contour on',
                'q: toggle q plot','s: toggle single plot','w: toggle divide by sig','+: no selection')
    cb = wx.ComboBox(G2frame.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=Choice)
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' key press')
    
    PickId = G2frame.PickId
    PatternId = G2frame.PatternId
    colors=['b','g','r','c','m','k']
    Lines = []
    if G2frame.SinglePlot:
        Pattern = G2frame.PatternTree.GetItemPyData(PatternId)
        Pattern.append(G2frame.PatternTree.GetItemText(PatternId))
        PlotList = [Pattern,]
        ParmList = [G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId, 'Instrument Parameters'))[1],]
    else:        
        PlotList = []
        ParmList = []
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while item:
            if 'PWDR' in G2frame.PatternTree.GetItemText(item):
                Pattern = G2frame.PatternTree.GetItemPyData(item)
                if len(Pattern) < 3:                    # put name on end if needed
                    Pattern.append(G2frame.PatternTree.GetItemText(item))
                PlotList.append(Pattern)
                ParmList.append(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
                    item,'Instrument Parameters'))[1])
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)                
    Ymax = 1.0
    lenX = 0
    if PickId:
        if G2frame.PatternTree.GetItemText(PickId) in ['Reflection Lists']:
            Phases = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            HKL = []
            if Phases:
                for peak in Phases[G2frame.RefList]:
                    HKL.append(peak[:6])
                HKL = np.array(HKL)
        else:
            HKL = np.array(G2frame.HKL)
    for Pattern in PlotList:
        xye = Pattern[1]
        Ymax = max(Ymax,max(xye[1]))
    offset = G2frame.Offset[0]*Ymax/100.0
    Title = 'Powder Patterns: '+os.path.split(G2frame.GSASprojectfile)[1]
    if G2frame.logPlot:
        Title = 'log('+Title+')'
    Plot.set_title(Title)
    if G2frame.qPlot:
        Plot.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
    else:        
        Plot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    if G2frame.Weight:
        Plot.set_ylabel(r'$\mathsf{I/\sigma(I)}$',fontsize=14)
    else:
        Plot.set_ylabel('Intensity',fontsize=12)
    if G2frame.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    if len(PlotList) < 2:
        G2frame.Contour = False
    for N,Pattern in enumerate(PlotList):
        Parms = ParmList[N]
        ifpicked = False
        LimitId = 0
        xye = np.array(Pattern[1])
        if PickId:
            ifpicked = Pattern[2] == G2frame.PatternTree.GetItemText(PatternId)
            LimitId = G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits')
        if G2frame.qPlot:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, Pattern[2])
            X = 4*np.pi*npsind(xye[0]/2.0)/Parms[1]
        else:
            X = xye[0]
        if not lenX:
            lenX = len(X)           
        Y = xye[1]+offset*N
        if LimitId:
            limits = np.array(G2frame.PatternTree.GetItemPyData(LimitId))
            if G2frame.qPlot:
                limits = 4*np.pi*npsind(limits/2.0)/Parms[1]
            Lines.append(Plot.axvline(limits[1][0],color='g',dashes=(5,5),picker=3.))    
            Lines.append(Plot.axvline(limits[1][1],color='r',dashes=(5,5),picker=3.))                    
        if G2frame.Contour:
            
            if lenX == len(X):
                ContourY.append(N)
                ContourZ.append(Y)
                ContourX = X
                Nseq += 1
                Plot.set_ylabel('Data sequence',fontsize=12)
        else:
            X += G2frame.Offset[1]*.005*N
            if ifpicked:
                Z = xye[3]+offset*N
                W = xye[4]+offset*N
                D = xye[5]-Ymax*G2frame.delOffset
                if G2frame.logPlot:
                    Plot.semilogy(X,Y,colors[N%6]+'+',picker=3.,clip_on=False,nonposy='mask')
                    Plot.semilogy(X,Z,colors[(N+1)%6],picker=False,nonposy='mask')
                    Plot.semilogy(X,W,colors[(N+2)%6],picker=False,nonposy='mask')
                elif G2frame.Weight:
                    DY = xye[1]*np.sqrt(xye[2])
                    DYmax = max(DY)
                    DZ = xye[3]*np.sqrt(xye[2])
                    DS = xye[5]*np.sqrt(xye[2])-DYmax*G2frame.delOffset
                    Plot.plot(X,DY,colors[N%6]+'+',picker=3.,clip_on=False)
                    Plot.plot(X,DZ,colors[(N+1)%6],picker=False)
                    Plot.plot(X,DS,colors[(N+3)%6],picker=False)
                    Plot.axhline(0.,color=wx.BLACK)
                else:
                    Plot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                    Plot.plot(X,Z,colors[(N+1)%6],picker=False)
                    Plot.plot(X,W,colors[(N+2)%6],picker=False)
                    Plot.plot(X,D,colors[(N+3)%6],picker=False)
                    Plot.axhline(0.,color=wx.BLACK)
                Page.canvas.SetToolTipString('')
                if G2frame.PatternTree.GetItemText(PickId) == 'Peak List':
                    tip = 'On data point: Pick peak - L or R MB. On line: L-move, R-delete'
                    Page.canvas.SetToolTipString(tip)
                    data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
                    for item in data:
                        if G2frame.qPlot:
                            Lines.append(Plot.axvline(4*math.pi*sind(item[0]/2.)/Parms[1],color=colors[N%6],picker=2.))
                        else:
                            Lines.append(Plot.axvline(item[0],color=colors[N%6],picker=2.))
                if G2frame.PatternTree.GetItemText(PickId) == 'Limits':
                    tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                    Page.canvas.SetToolTipString(tip)
                    data = G2frame.LimitsTable.GetData()
            else:
                if G2frame.logPlot:
                    Plot.semilogy(X,Y,colors[N%6],picker=False,nonposy='clip')
                else:
                    Plot.plot(X,Y,colors[N%6],picker=False)
    if PickId and not G2frame.Contour:
        Values,Names = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))[1::2]
        Parms = dict(zip(Names,Values))
        try:
            wave = Parms['Lam']
        except KeyError:
            wave = Parms['Lam1']
        if G2frame.PatternTree.GetItemText(PickId) in ['Index Peak List','Unit Cells List']:
            peaks = np.array((G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))))
            for peak in peaks:
                if G2frame.qPlot:
                    Plot.axvline(4*np.pi*sind(peak[0]/2.0)/wave,color='b')
                else:
                    Plot.axvline(peak[0],color='b')
            for hkl in G2frame.HKL:
                if G2frame.qPlot:
                    Plot.axvline(4*np.pi*sind(hkl[5]/2.0)/wave,color='r',dashes=(5,5))
                else:
                    Plot.axvline(hkl[5],color='r',dashes=(5,5))
        elif G2frame.PatternTree.GetItemText(PickId) in ['Reflection Lists'] or \
            'PWDR' in G2frame.PatternTree.GetItemText(PickId):
            Phases = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            for pId,phase in enumerate(Phases):
                peaks = Phases[phase]
                if not peaks:
                    continue
                peak = np.array([[peak[4],peak[5]] for peak in peaks])
                pos = G2frame.refOffset-pId*Ymax*G2frame.refDelt*np.ones_like(peak)
                if G2frame.qPlot:
                    Plot.plot(2*np.pi/peak.T[0],pos,colors[pId%6]+'|',mew=1,ms=8,picker=3.,label=phase)
                else:
                    Plot.plot(peak.T[1],pos,colors[pId%6]+'|',mew=1,ms=8,picker=3.,label=phase)
            if len(Phases):
                handles,legends = Plot.get_legend_handles_labels()  #got double entries in the legends for some reason
                if handles:
                    Plot.legend(handles[::2],legends[::2],title='Phases',loc='best')    #skip every other one
            
    if G2frame.Contour:
        acolor = mpl.cm.get_cmap(G2frame.ContourColor)
        Img = Plot.imshow(ContourZ,cmap=acolor,vmin=0,vmax=Ymax*G2frame.Cmax,interpolation=G2frame.Interpolate, 
            extent=[ContourX[0],ContourX[-1],ContourY[0],ContourY[-1]],aspect='auto',origin='lower')
        Page.figure.colorbar(Img)
    else:
        G2frame.Lines = Lines
    if not newPlot:
        Page.toolbar.push_current()
        Plot.set_xlim(xylim[0])
        Plot.set_ylim(xylim[1])
        xylim = []
        Page.toolbar.push_current()
        Page.toolbar.draw()
    else:
        Page.canvas.draw()
#    G2frame.Pwdr = True
    
################################################################################
##### PlotDeltSig
################################################################################
            
def PlotDeltSig(G2frame,kind):
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Error analysis')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
    except ValueError:
        newPlot = True
        G2frame.Cmax = 1.0
        Plot = G2frame.G2plotNB.addMpl('Error analysis').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Error analysis')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
    PatternId = G2frame.PatternId
    Pattern = G2frame.PatternTree.GetItemPyData(PatternId)
    Pattern.append(G2frame.PatternTree.GetItemText(PatternId))
    wtFactor = Pattern[0]['wtFactor']
    if kind == 'PWDR':
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        xye = np.array(Pattern[1])
        xmin = np.searchsorted(xye[0],limits[0])
        xmax = np.searchsorted(xye[0],limits[1])
        DS = xye[5][xmin:xmax]*np.sqrt(wtFactor*xye[2][xmin:xmax])
    elif kind == 'HKLF':
        refl = Pattern[1]
        DS = []
        for ref in refl:
            if ref[6] > 0.:
                DS.append((ref[5]-ref[7])/ref[6])
    Page.SetFocus()
    G2frame.G2plotNB.status.DestroyChildren()
    DS.sort()
    EDS = np.zeros_like(DS)
    DX = np.linspace(0.,1.,num=len(DS),endpoint=True)
    oldErr = np.seterr(invalid='ignore')    #avoid problem at DX==0
    T = np.sqrt(np.log(1.0/DX**2))
    top = 2.515517+0.802853*T+0.010328*T**2
    bot = 1.0+1.432788*T+0.189269*T**2+0.001308*T**3
    EDS = np.where(DX>0,-(T-top/bot),(T-top/bot))
    low1 = np.searchsorted(EDS,-1.)
    hi1 = np.searchsorted(EDS,1.)
    slp,intcp = np.polyfit(EDS[low1:hi1],DS[low1:hi1],deg=1)
    frac = 100.*(hi1-low1)/len(DS)
    G2frame.G2plotNB.status.SetStatusText(  \
        'Over range -1. to 1. :'+' slope = %.3f, intercept = %.3f for %.2f%% of the fitted data'%(slp,intcp,frac),1)
    Plot.set_title('Normal probability for '+Pattern[-1])
    Plot.set_xlabel(r'expected $\mathsf{\Delta/\sigma}$',fontsize=14)
    Plot.set_ylabel(r'observed $\mathsf{\Delta/\sigma}$',fontsize=14)
    Plot.plot(EDS,DS,'r+',label='result')
    Plot.plot([-2,2],[-2,2],'k',dashes=(5,5),label='ideal')
    Plot.legend(loc='upper left')
    np.seterr(invalid='warn')
    Page.canvas.draw()
       
################################################################################
##### PlotISFG
################################################################################
            
def PlotISFG(G2frame,newPlot=False,type=''):
    ''' PLotting package for PDF analysis; displays I(q), S(q), F(q) and G(r) as single 
    or multiple plots with waterfall and contour plots as options
    '''
    if not type:
        type = G2frame.G2plotNB.plotList[G2frame.G2plotNB.nb.GetSelection()]
    if type not in ['I(Q)','S(Q)','F(Q)','G(R)']:
        return
    superMinusOne = unichr(0xaf)+unichr(0xb9)
    
    def OnPlotKeyPress(event):
        newPlot = False
        if event.key == 'u':
            if G2frame.Contour:
                G2frame.Cmax = min(1.0,G2frame.Cmax*1.2)
            elif G2frame.Offset[0] < 100.:
                G2frame.Offset[0] += 1.
        elif event.key == 'd':
            if G2frame.Contour:
                G2frame.Cmax = max(0.0,G2frame.Cmax*0.8)
            elif G2frame.Offset[0] > 0.:
                G2frame.Offset[0] -= 1.
        elif event.key == 'l':
            G2frame.Offset[1] -= 1.
        elif event.key == 'r':
            G2frame.Offset[1] += 1.
        elif event.key == 'o':
            G2frame.Offset = [0,0]
        elif event.key == 'c':
            newPlot = True
            if G2frame.Contour:
                G2frame.Contour = False
            else:
                G2frame.Contour = True
                G2frame.SinglePlot = False
                G2frame.Offset = [0.,0.]
        elif event.key == 's':
            if G2frame.Contour:
                choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
                choice.sort()
                dlg = wx.SingleChoiceDialog(G2frame,'Select','Color scheme',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    G2frame.ContourColor = choice[sel]
                else:
                    G2frame.ContourColor = 'Paired'
                dlg.Destroy()
            else:                
                if G2frame.SinglePlot:
                    G2frame.SinglePlot = False
                else:
                    G2frame.SinglePlot = True
        elif event.key == 'i':                  #for smoothing contour plot
            choice = ['nearest','bilinear','bicubic','spline16','spline36','hanning',
               'hamming','hermite','kaiser','quadric','catrom','gaussian','bessel',
               'mitchell','sinc','lanczos']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Interpolation',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.Interpolate = choice[sel]
            else:
                G2frame.Interpolate = 'nearest'
            dlg.Destroy()
        elif event.key == 't' and not G2frame.Contour:
            if G2frame.Legend:
                G2frame.Legend = False
            else:
                G2frame.Legend = True
            
            
        PlotISFG(G2frame,newPlot=newPlot,type=type)
        
    def OnKeyBox(event):
        if G2frame.G2plotNB.nb.GetSelection() == G2frame.G2plotNB.plotList.index(type):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            wx.CallAfter(OnPlotKeyPress,event)
                        
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                if G2frame.Contour:
                    G2frame.G2plotNB.status.SetStatusText('R =%.3fA pattern ID =%5d'%(xpos,int(ypos)),1)
                else:
                    G2frame.G2plotNB.status.SetStatusText('R =%.3fA %s =%.2f'%(xpos,type,ypos),1)                   
            except TypeError:
                G2frame.G2plotNB.status.SetStatusText('Select '+type+' pattern first',1)
    
    xylim = []
    try:
        plotNum = G2frame.G2plotNB.plotList.index(type)
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        newPlot = True
        G2frame.Cmax = 1.0
        Plot = G2frame.G2plotNB.addMpl(type).gca()
        plotNum = G2frame.G2plotNB.plotList.index(type)
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
    
    Page.SetFocus()
    G2frame.G2plotNB.status.DestroyChildren()
    if G2frame.Contour:
        Choice = (' key press','d: lower contour max','u: raise contour max',
            'i: interpolation method','s: color scheme','c: contour off')
    else:
        Choice = (' key press','l: offset left','r: offset right','d: offset down','u: offset up',
            'o: reset offset','t: toggle legend','c: contour on','s: toggle single plot')
    cb = wx.ComboBox(G2frame.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=Choice)
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' key press')
    PatternId = G2frame.PatternId
    PickId = G2frame.PickId
    Plot.set_title(type)
    if type == 'G(R)':
        Plot.set_xlabel(r'$R,\AA$',fontsize=14)
    else:
        Plot.set_xlabel(r'$Q,\AA$'+superMinusOne,fontsize=14)
    Plot.set_ylabel(r''+type,fontsize=14)
    colors=['b','g','r','c','m','k']
    name = G2frame.PatternTree.GetItemText(PatternId)[4:]
    Pattern = []    
    if G2frame.SinglePlot:
        name = G2frame.PatternTree.GetItemText(PatternId)
        name = type+name[4:]
        Id = G2gd.GetPatternTreeItemId(G2frame,PatternId,name)
        Pattern = G2frame.PatternTree.GetItemPyData(Id)
        if Pattern:
            Pattern.append(name)
        PlotList = [Pattern,]
    else:
        PlotList = []
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while item:
            if 'PDF' in G2frame.PatternTree.GetItemText(item):
                name = type+G2frame.PatternTree.GetItemText(item)[4:]
                Id = G2gd.GetPatternTreeItemId(G2frame,item,name)
                Pattern = G2frame.PatternTree.GetItemPyData(Id)
                if Pattern:
                    Pattern.append(name)
                    PlotList.append(Pattern)
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    PDFdata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'PDF Controls'))
    numbDen = G2pwd.GetNumDensity(PDFdata['ElList'],PDFdata['Form Vol'])
    Xb = [0.,10.]
    Yb = [0.,-40.*np.pi*numbDen]
    Ymax = 1.0
    lenX = 0
    for Pattern in PlotList:
        xye = Pattern[1]
        Ymax = max(Ymax,max(xye[1]))
    offset = G2frame.Offset[0]*Ymax/100.0
    if G2frame.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    for N,Pattern in enumerate(PlotList):
        xye = Pattern[1]
        if PickId:
            ifpicked = Pattern[2] == G2frame.PatternTree.GetItemText(PatternId)
        X = xye[0]
        if not lenX:
            lenX = len(X)           
        Y = xye[1]+offset*N
        if G2frame.Contour:
            if lenX == len(X):
                ContourY.append(N)
                ContourZ.append(Y)
                ContourX = X
                Nseq += 1
                Plot.set_ylabel('Data sequence',fontsize=12)
        else:
            X = xye[0]+G2frame.Offset[1]*.005*N
            if ifpicked:
                Plot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                Page.canvas.SetToolTipString('')
            else:
                if G2frame.Legend:
                    Plot.plot(X,Y,colors[N%6],picker=False,label='Azm:'+Pattern[2].split('=')[1])
                else:
                    Plot.plot(X,Y,colors[N%6],picker=False)
            if type == 'G(R)':
                Plot.plot(Xb,Yb,color='k',dashes=(5,5))
            elif type == 'F(Q)':
                Plot.axhline(0.,color=wx.BLACK)
            elif type == 'S(Q)':
                Plot.axhline(1.,color=wx.BLACK)
    if G2frame.Contour:
        acolor = mpl.cm.get_cmap(G2frame.ContourColor)
        Img = Plot.imshow(ContourZ,cmap=acolor,vmin=0,vmax=Ymax*G2frame.Cmax,interpolation=G2frame.Interpolate, 
            extent=[ContourX[0],ContourX[-1],ContourY[0],ContourY[-1]],aspect='auto',origin='lower')
        Page.figure.colorbar(Img)
    elif G2frame.Legend:
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
        
################################################################################
##### PlotXY
################################################################################
            
def PlotXY(G2frame,XY,newPlot=False,type=''):
    '''simple plot of xy data, used for diagnostic purposes
    '''
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            try:
                G2frame.G2plotNB.status.SetStatusText('X =%9.3f %s =%9.3f'%(xpos,type,ypos),1)                   
            except TypeError:
                G2frame.G2plotNB.status.SetStatusText('Select '+type+' pattern first',1)

    try:
        plotNum = G2frame.G2plotNB.plotList.index(type)
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        newPlot = True
        Plot = G2frame.G2plotNB.addMpl(type).gca()
        plotNum = G2frame.G2plotNB.plotList.index(type)
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
    
    Page.SetFocus()
    G2frame.G2plotNB.status.DestroyChildren()
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

################################################################################
##### PlotPowderLines
################################################################################
            
def PlotPowderLines(G2frame):
    ''' plotting of powder lines (i.e. no powder pattern) as sticks
    '''
    global HKL

    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            G2frame.G2plotNB.status.SetFields(['','2-theta =%9.3f '%(xpos,)])
            if G2frame.PickId and G2frame.PatternTree.GetItemText(G2frame.PickId) in ['Index Peak List','Unit Cells List']:
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
        plotNum = G2frame.G2plotNB.plotList.index('Powder Lines')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('Powder Lines').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Powder Lines')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        
    Page.SetFocus()
    Plot.set_title('Powder Pattern Lines')
    Plot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    PickId = G2frame.PickId
    PatternId = G2frame.PatternId
    peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))
    for peak in peaks:
        Plot.axvline(peak[0],color='b')
    HKL = np.array(G2frame.HKL)
    for hkl in G2frame.HKL:
        Plot.axvline(hkl[5],color='r',dashes=(5,5))
    xmin = peaks[0][0]
    xmax = peaks[-1][0]
    delt = xmax-xmin
    xlim = [max(0,xmin-delt/20.),min(180.,xmax+delt/20.)]
    Plot.set_xlim(xlim)
    Page.canvas.draw()
    Page.toolbar.push_current()

################################################################################
##### PlotPeakWidths
################################################################################
            
def PlotPeakWidths(G2frame):
    ''' Plotting of instrument broadening terms as function of 2-theta (future TOF)
    Seen when "Instrument Parameters" chosen from powder pattern data tree
    '''
    PatternId = G2frame.PatternId
    limitID = G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits')
    if limitID:
        limits = G2frame.PatternTree.GetItemPyData(limitID)
    else:
        return
    instParms = G2frame.PatternTree.GetItemPyData( \
        G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
    if instParms[0][0] in ['PXC','PNC']:
        lam = instParms[1][1]
        if len(instParms[1]) == 13:
            GU,GV,GW,LX,LY = instParms[0][6:11]
        else:
            GU,GV,GW,LX,LY = instParms[0][4:9]
    peakID = G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List')
    if peakID:
        peaks = G2frame.PatternTree.GetItemPyData(peakID)
    else:
        peaks = []
    
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Peak Widths')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('Peak Widths').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Peak Widths')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
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
    sig = lambda Th,U,V,W: 1.17741*math.sqrt(U*tand(Th)**2+V*tand(Th)+W)*math.pi/18000.
    gam = lambda Th,X,Y: (X/cosd(Th)+Y*tand(Th))*math.pi/18000.
    gamFW = lambda s,g: math.exp(math.log(s**5+2.69269*s**4*g+2.42843*s**3*g**2+4.47163*s**2*g**3+0.07842*s*g**4+g**5)/5.)
#    gamFW2 = lambda s,g: math.sqrt(s**2+(0.4654996*g)**2)+.5345004*g  #Ubaldo Bafile - private communication
    try:
        for theta in thetas:
            X.append(4.0*math.pi*sind(theta)/lam)              #q
            s = sig(theta,GU,GV,GW)
            g = gam(theta,LX,LY)
            G = gamFW(g,s)
            Y.append(s/tand(theta))
            Z.append(g/tand(theta))
            W.append(G/tand(theta))
        Plot.set_title('Instrument and sample peak widths')
        Plot.set_ylabel(r'$\Delta q/q, \Delta d/d$',fontsize=14)
        Plot.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
        Plot.plot(X,Y,color='r',label='Gaussian')
        Plot.plot(X,Z,color='g',label='Lorentzian')
        Plot.plot(X,W,color='b',label='G+L')
        X = []
        Y = []
        Z = []
        W = []
        V = []
        for peak in peaks:
            X.append(4.0*math.pi*sind(peak[0]/2.0)/lam)
            try:
                s = 1.17741*math.sqrt(peak[4])*math.pi/18000.
            except ValueError:
                s = 0.01
            g = peak[6]*math.pi/18000.
            G = gamFW(g,s)
            Y.append(s/tand(peak[0]/2.))
            Z.append(g/tand(peak[0]/2.))
            W.append(G/tand(peak[0]/2.))
        Plot.plot(X,Y,'+',color='r',label='G peak')
        Plot.plot(X,Z,'+',color='g',label='L peak')
        Plot.plot(X,W,'+',color='b',label='G+L peak')
        Plot.legend(loc='best')
        Page.canvas.draw()
    except ValueError:
        print '**** ERROR - default U,V,W profile coefficients yield sqrt of negative value at 2theta =', \
            '%.3f'%(2*theta)
        G2frame.G2plotNB.Delete('Peak Widths')

    
################################################################################
##### PlotSizeStrainPO
################################################################################
            
def PlotSizeStrainPO(G2frame,data,Start=False):
    '''Plot 3D mustrain/size/preferred orientation figure. In this instance data is for a phase
    '''
    
    PatternId = G2frame.PatternId
    generalData = data['General']
    SGData = generalData['SGData']
    SGLaue = SGData['SGLaue']
    if Start:                   #initialize the spherical harmonics qlmn arrays
        ptx.pyqlmninit()
        Start = False
#    MuStrKeys = G2spc.MustrainNames(SGData)
    cell = generalData['Cell'][1:]
    A,B = G2lat.cell2AB(cell[:6])
    Vol = cell[6]
    useList = data['Histograms']
    phase = generalData['Name']
    plotType = generalData['Data plot type']
    plotDict = {'Mustrain':'Mustrain','Size':'Size','Preferred orientation':'Pref.Ori.'}
    for ptype in plotDict:
        G2frame.G2plotNB.Delete(ptype)
    if plotType in ['None']:
        return        

    for item in useList:
        if useList[item]['Show']:
            break
    else:
        return            #nothing to show!!
    
    numPlots = len(useList)

    if plotType in ['Mustrain','Size']:
        Plot = mp3d.Axes3D(G2frame.G2plotNB.add3D(plotType))
    else:
        Plot = G2frame.G2plotNB.addMpl(plotType).gca()        
    plotNum = G2frame.G2plotNB.plotList.index(plotType)
    Page = G2frame.G2plotNB.nb.GetPage(plotNum)
    Page.SetFocus()
    G2frame.G2plotNB.status.SetStatusText('',1)
    if not Page.IsShown():
        Page.Show()
    
    for item in useList:
        if useList[item]['Show']:
            PHI = np.linspace(0.,360.,30,True)
            PSI = np.linspace(0.,180.,30,True)
            X = np.outer(npsind(PHI),npsind(PSI))
            Y = np.outer(npcosd(PHI),npsind(PSI))
            Z = np.outer(np.ones(np.size(PHI)),npcosd(PSI))
            try:        #temp patch instead of 'mustrain' for old files with 'microstrain'
                coeff = useList[item][plotDict[plotType]]
            except KeyError:
                break
            if plotType in ['Mustrain','Size']:
                if coeff[0] == 'isotropic':
                    X *= coeff[1][0]
                    Y *= coeff[1][0]
                    Z *= coeff[1][0]                                
                elif coeff[0] == 'uniaxial':
                    
                    def uniaxCalc(xyz,iso,aniso,axes):
                        Z = np.array(axes)
                        cp = abs(np.dot(xyz,Z))
                        sp = np.sqrt(1.-cp**2)
                        R = iso*aniso/np.sqrt((iso*cp)**2+(aniso*sp)**2)
                        return R*xyz
                        
                    iso,aniso = coeff[1][:2]
                    axes = np.inner(A,np.array(coeff[3]))
                    axes /= nl.norm(axes)
                    Shkl = np.array(coeff[1])
                    XYZ = np.dstack((X,Y,Z))
                    XYZ = np.nan_to_num(np.apply_along_axis(uniaxCalc,2,XYZ,iso,aniso,axes))
                    X,Y,Z = np.dsplit(XYZ,3)
                    X = X[:,:,0]
                    Y = Y[:,:,0]
                    Z = Z[:,:,0]
                
                elif coeff[0] == 'ellipsoidal':
                    
                    def ellipseCalc(xyz,E,R):
                        XYZ = xyz*E.T
                        return np.inner(XYZ.T,R)
                        
                    S6 = coeff[4]
                    Sij = G2lat.U6toUij(S6)
                    E,R = nl.eigh(Sij)
                    XYZ = np.dstack((X,Y,Z))
                    XYZ = np.nan_to_num(np.apply_along_axis(ellipseCalc,2,XYZ,E,R))
                    X,Y,Z = np.dsplit(XYZ,3)
                    X = X[:,:,0]
                    Y = Y[:,:,0]
                    Z = Z[:,:,0]
                    
                elif coeff[0] == 'generalized':
                    
                    def genMustrain(xyz,SGData,A,Shkl):
                        uvw = np.inner(A.T,xyz)
                        Strm = np.array(G2spc.MustrainCoeff(uvw,SGData))
                        sum = np.sum(np.multiply(Shkl,Strm))
                        sum = np.where(sum > 0.01,sum,0.01)
                        sum = np.sqrt(sum)*math.pi/0.018      #centidegrees to radians!
                        return sum*xyz
                        
                    Shkl = np.array(coeff[4])
                    if np.any(Shkl):
                        XYZ = np.dstack((X,Y,Z))
                        XYZ = np.nan_to_num(np.apply_along_axis(genMustrain,2,XYZ,SGData,A,Shkl))
                        X,Y,Z = np.dsplit(XYZ,3)
                        X = X[:,:,0]
                        Y = Y[:,:,0]
                        Z = Z[:,:,0]
                            
                if np.any(X) and np.any(Y) and np.any(Z):
                    errFlags = np.seterr(all='ignore')
                    Plot.plot_surface(X,Y,Z,rstride=1,cstride=1,color='g',linewidth=1)
                    np.seterr(all='ignore')
                    xyzlim = np.array([Plot.get_xlim3d(),Plot.get_ylim3d(),Plot.get_zlim3d()]).T
                    XYZlim = [min(xyzlim[0]),max(xyzlim[1])]
                    Plot.set_xlim3d(XYZlim)
                    Plot.set_ylim3d(XYZlim)
                    Plot.set_zlim3d(XYZlim)
                    Plot.set_aspect('equal')
                if plotType == 'Size':
                    Plot.set_title('Crystallite size for '+phase+'\n'+coeff[0]+' model')
                    Plot.set_xlabel(r'X, $\mu$m')
                    Plot.set_ylabel(r'Y, $\mu$m')
                    Plot.set_zlabel(r'Z, $\mu$m')
                else:    
                    Plot.set_title(r'$\mu$strain for '+phase+'\n'+coeff[0]+' model')
                    Plot.set_xlabel(r'X, $\mu$strain')
                    Plot.set_ylabel(r'Y, $\mu$strain')
                    Plot.set_zlabel(r'Z, $\mu$strain')
            else:
                h,k,l = generalData['POhkl']
                if coeff[0] == 'MD':
                    print 'March-Dollase preferred orientation plot'
                
                else:
                    PH = np.array(generalData['POhkl'])
                    phi,beta = G2lat.CrsAng(PH,cell[:6],SGData)
                    SHCoef = {}
                    for item in coeff[5]:
                        L,N = eval(item.strip('C'))
                        SHCoef['C%d,0,%d'%(L,N)] = coeff[5][item]                        
                    ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
                    X = np.linspace(0,90.0,26)
                    Y = G2lat.polfcal(ODFln,'0',X,0.0)
                    Plot.plot(X,Y,color='k',label=str(PH))
                    Plot.legend(loc='best')
                    Plot.set_title('Axial distribution for HKL='+str(PH)+' in '+phase)
                    Plot.set_xlabel(r'$\psi$',fontsize=16)
                    Plot.set_ylabel('MRD',fontsize=14)
    Page.canvas.draw()
    
################################################################################
##### PlotTexture
################################################################################
            
def PlotTexture(G2frame,data,Start=False):
    '''Pole figure, inverse pole figure, 3D pole distribution and 3D inverse pole distribution
    plotting.
    dict generalData contains all phase info needed which is in data
    '''

    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    PatternId = G2frame.PatternId
    generalData = data['General']
    SGData = generalData['SGData']
    textureData = generalData['SH Texture']
    G2frame.G2plotNB.Delete('Texture')
    if not textureData['Order']:
        return                  #no plot!!
    SHData = generalData['SH Texture']
    SHCoef = SHData['SH Coeff'][1]
    cell = generalData['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(cell)
    sq2 = 1.0/math.sqrt(2.0)
    
    def rp2xyz(r,p):
        z = npcosd(r)
        xy = np.sqrt(1.-z**2)
        return xy*npsind(p),xy*npcosd(p),z
            
    def OnMotion(event):
        SHData = data['General']['SH Texture']
        if event.xdata and event.ydata:                 #avoid out of frame errors
            xpos = event.xdata
            ypos = event.ydata
            if 'Inverse' in SHData['PlotType']:
                r = xpos**2+ypos**2
                if r <= 1.0:
                    if 'equal' in G2frame.Projection: 
                        r,p = 2.*npasind(np.sqrt(r)*sq2),npatan2d(ypos,xpos)
                    else:
                        r,p = 2.*npatand(np.sqrt(r)),npatan2d(ypos,xpos)
                    ipf = G2lat.invpolfcal(IODFln,SGData,np.array([r,]),np.array([p,]))
                    xyz = np.inner(Amat.T,np.array([rp2xyz(r,p)]))
                    y,x,z = list(xyz/np.max(np.abs(xyz)))
                    
                    G2frame.G2plotNB.status.SetFields(['',
                        'psi =%9.3f, beta =%9.3f, MRD =%9.3f hkl=%5.2f,%5.2f,%5.2f'%(r,p,ipf,x,y,z)])
                                    
            elif 'Axial' in SHData['PlotType']:
                pass
                
            else:                       #ordinary pole figure
                z = xpos**2+ypos**2
                if z <= 1.0:
                    z = np.sqrt(z)
                    if 'equal' in G2frame.Projection: 
                        r,p = 2.*npasind(z*sq2),npatan2d(ypos,xpos)
                    else:
                        r,p = 2.*npatand(z),npatan2d(ypos,xpos)
                    pf = G2lat.polfcal(ODFln,SamSym[textureData['Model']],np.array([r,]),np.array([p,]))
                    G2frame.G2plotNB.status.SetFields(['','phi =%9.3f, gam =%9.3f, MRD =%9.3f'%(r,p,pf)])
    
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Texture')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
        if not Page.IsShown():
            Page.Show()
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('Texture').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Texture')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)

    Page.SetFocus()
    G2frame.G2plotNB.status.SetFields(['',''])    
    PH = np.array(SHData['PFhkl'])
    phi,beta = G2lat.CrsAng(PH,cell,SGData)
    ODFln = G2lat.Flnh(Start,SHCoef,phi,beta,SGData)
    PX = np.array(SHData['PFxyz'])
    gam = atan2d(PX[0],PX[1])
    xy = math.sqrt(PX[0]**2+PX[1]**2)
    xyz = math.sqrt(PX[0]**2+PX[1]**2+PX[2]**2)
    psi = asind(xy/xyz)
    IODFln = G2lat.Glnh(Start,SHCoef,psi,gam,SamSym[textureData['Model']])
    if 'Axial' in SHData['PlotType']:
        X = np.linspace(0,90.0,26)
        Y = G2lat.polfcal(ODFln,SamSym[textureData['Model']],X,0.0)
        Plot.plot(X,Y,color='k',label=str(SHData['PFhkl']))
        Plot.legend(loc='best')
        Plot.set_title('Axial distribution for HKL='+str(SHData['PFhkl']))
        Plot.set_xlabel(r'$\psi$',fontsize=16)
        Plot.set_ylabel('MRD',fontsize=14)
        
    else:       
        npts = 201
        if 'Inverse' in SHData['PlotType']:
            X,Y = np.meshgrid(np.linspace(1.,-1.,npts),np.linspace(-1.,1.,npts))
            R,P = np.sqrt(X**2+Y**2).flatten(),npatan2d(X,Y).flatten()
            if 'equal' in G2frame.Projection:
                R = np.where(R <= 1.,2.*npasind(R*sq2),0.0)
            else:
                R = np.where(R <= 1.,2.*npatand(R),0.0)
            Z = np.zeros_like(R)
            Z = G2lat.invpolfcal(IODFln,SGData,R,P)
            Z = np.reshape(Z,(npts,npts))
            CS = Plot.contour(Y,X,Z,aspect='equal')
            Plot.clabel(CS,fontsize=9,inline=1)
            try:
                Img = Plot.imshow(Z.T,aspect='equal',cmap=G2frame.ContourColor,extent=[-1,1,-1,1])
            except ValueError:
                pass
            Page.figure.colorbar(Img)
            Plot.set_title('Inverse pole figure for XYZ='+str(SHData['PFxyz']))
            Plot.set_xlabel(G2frame.Projection.capitalize()+' projection')
                        
        else:
            X,Y = np.meshgrid(np.linspace(1.,-1.,npts),np.linspace(-1.,1.,npts))
            R,P = np.sqrt(X**2+Y**2).flatten(),npatan2d(X,Y).flatten()
            if 'equal' in G2frame.Projection:
                R = np.where(R <= 1.,2.*npasind(R*sq2),0.0)
            else:
                R = np.where(R <= 1.,2.*npatand(R),0.0)
            Z = np.zeros_like(R)
            Z = G2lat.polfcal(ODFln,SamSym[textureData['Model']],R,P)
            Z = np.reshape(Z,(npts,npts))
            try:
                CS = Plot.contour(Y,X,Z,aspect='equal')
                Plot.clabel(CS,fontsize=9,inline=1)
            except ValueError:
                pass
            Img = Plot.imshow(Z.T,aspect='equal',cmap=G2frame.ContourColor,extent=[-1,1,-1,1])
            Page.figure.colorbar(Img)
            Plot.set_title('Pole figure for HKL='+str(SHData['PFhkl']))
            Plot.set_xlabel(G2frame.Projection.capitalize()+' projection')
    Page.canvas.draw()

################################################################################
##### PlotCovariance
################################################################################
            
def PlotCovariance(G2frame,Data={}):
    if not Data:
        Data = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Covariance'))
    if not Data:
        print 'No covariance matrix available'
        return
    varyList = Data['varyList']
    values = Data['variables']
    Xmax = len(varyList)
    covMatrix = Data['covMatrix']
    sig = np.sqrt(np.diag(covMatrix))
    xvar = np.outer(sig,np.ones_like(sig))
    covArray = np.divide(np.divide(covMatrix,xvar),xvar.T)
    title = ' for\n'+Data['title']
    newAtomDict = Data['newAtomDict']

    def OnPlotKeyPress(event):
        newPlot = False
        if event.key == 's':
            choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
            choice.sort()
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Color scheme',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.VcovColor = choice[sel]
            else:
                G2frame.VcovColor = 'RdYlGn'
            dlg.Destroy()
        PlotCovariance(G2frame)

    def OnMotion(event):
        if event.button:
            ytics = imgAx.get_yticks()
            ytics = np.where(ytics<len(varyList),ytics,-1)
            ylabs = [np.where(0<=i ,varyList[int(i)],' ') for i in ytics]
            imgAx.set_yticklabels(ylabs)
            
        if event.xdata and event.ydata:                 #avoid out of frame errors
            xpos = int(event.xdata+.5)
            ypos = int(event.ydata+.5)
            if -1 < xpos < len(varyList) and -1 < ypos < len(varyList):
                if xpos == ypos:
                    value = values[xpos]
                    name = varyList[xpos]
                    if varyList[xpos] in newAtomDict:
                        name,value = newAtomDict[name]                        
                    msg = '%s value = %.4g, esd = %.4g'%(name,value,sig[xpos])
                else:
                    msg = '%s - %s: %5.3f'%(varyList[xpos],varyList[ypos],covArray[xpos][ypos])
                Page.canvas.SetToolTipString(msg)
                G2frame.G2plotNB.status.SetFields(['Key: s to change colors',msg])
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Covariance')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
        if not Page.IsShown():
            Page.Show()
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('Covariance').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Covariance')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)

    Page.SetFocus()
    G2frame.G2plotNB.status.SetFields(['',''])    
    acolor = mpl.cm.get_cmap(G2frame.VcovColor)
    Img = Plot.imshow(covArray,aspect='equal',cmap=acolor,interpolation='nearest',origin='lower')
    imgAx = Img.get_axes()
    ytics = imgAx.get_yticks()
    ylabs = [varyList[int(i)] for i in ytics[:-1]]
    imgAx.set_yticklabels(ylabs)
    colorBar = Page.figure.colorbar(Img)
    Plot.set_title('V-Cov matrix'+title)
    Plot.set_xlabel('Variable number')
    Plot.set_ylabel('Variable name')
    Page.canvas.draw()
    
################################################################################
##### PlotSeq
################################################################################
            
def PlotSeq(G2frame,SeqData,SeqSig,SeqNames,sampleParm):
    
    def OnKeyPress(event):
        if event.key == 's' and sampleParm:
            if G2frame.xAxis:
                G2frame.xAxis = False
            else:
                G2frame.xAxis = True
            Draw(False)
    try:
        plotNum = G2frame.G2plotNB.plotList.index('Sequential refinement')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.figure.clf()
        Plot = Page.figure.gca()
        if not Page.IsShown():
            Page.Show()
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('Sequential refinement').gca()
        plotNum = G2frame.G2plotNB.plotList.index('Sequential refinement')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnKeyPress)
        G2frame.xAxis = False
        
    def Draw(newPlot):
        Page.SetFocus()
        G2frame.G2plotNB.status.SetFields(['','press s to toggle x-axis = sample environment parameter'])
        if len(SeqData):
            Plot.clear()
            if G2frame.xAxis:    
                xName = sampleParm.keys()[0]
                X = sampleParm[xName]
            else:
                X = np.arange(0,len(SeqData[0]),1)
                xName = 'Data sequence number'
            for Y,sig,name in zip(SeqData,SeqSig,SeqNames):
                Plot.errorbar(X,Y,yerr=sig,label=name)        
            Plot.legend(loc='best')
            Plot.set_ylabel('Parameter values')
            Plot.set_xlabel(xName)
            Page.canvas.draw()            
    Draw(True)
            
################################################################################
##### PlotExposedImage & PlotImage
################################################################################
            
def PlotExposedImage(G2frame,newPlot=False,event=None):
    '''General access module for 2D image plotting
    '''
    plotNo = G2frame.G2plotNB.nb.GetSelection()
    if G2frame.G2plotNB.nb.GetPageText(plotNo) == '2D Powder Image':
        PlotImage(G2frame,newPlot,event,newImage=True)
    elif G2frame.G2plotNB.nb.GetPageText(plotNo) == '2D Integration':
        PlotIntegration(G2frame,newPlot,event)

def PlotImage(G2frame,newPlot=False,event=None,newImage=True):
    '''Plot of 2D detector images as contoured plot. Also plot calibration ellipses,
    masks, etc.
    '''
    from matplotlib.patches import Ellipse,Arc,Circle,Polygon
    import numpy.ma as ma
    Dsp = lambda tth,wave: wave/(2.*sind(tth/2.))
    global Data,Masks
    colors=['b','g','r','c','m','k']
    Data = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
    Masks = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
    try:    #may be absent
        StrSta = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Stress/Strain'))
    except TypeError:   #is missing
        StrSta = {}

    def OnImMotion(event):
        Page.canvas.SetToolTipString('')
        sizexy = Data['size']
        if event.xdata and event.ydata and len(G2frame.ImageZ):                 #avoid out of frame errors
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            item = G2frame.itemPicked
            pixelSize = Data['pixelSize']
            scalex = 1000./pixelSize[0]
            scaley = 1000./pixelSize[1]
            if item and G2frame.PatternTree.GetItemText(G2frame.PickId) == 'Image Controls':
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
                    Int = G2frame.ImageZ[ypix][xpix]
                tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                Q = 2.*math.pi/dsp
                if G2frame.setPoly:
                    G2frame.G2plotNB.status.SetFields(['','Polygon mask pick - LB next point, RB close polygon'])
                else:
                    G2frame.G2plotNB.status.SetFields(\
                        ['','Detector 2-th =%9.3fdeg, dsp =%9.3fA, Q = %6.5fA-1, azm = %7.2fdeg, I = %6d'%(tth,dsp,Q,azm,Int)])

    def OnImPlotKeyPress(event):
        try:
            PickName = G2frame.PatternTree.GetItemText(G2frame.PickId)
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
                G2frame.setPoly = True
                Masks['Polygons'].append([])
                G2frame.G2plotNB.status.SetFields(['','Polygon mask active - LB pick next point, RB close polygon'])
            G2imG.UpdateMasks(G2frame,Masks)
        elif PickName == 'Image Controls':
            if event.key == 'c':
                Xpos = event.xdata
                if not Xpos:            #got point out of frame
                    return
                Ypos = event.ydata
                dlg = wx.MessageDialog(G2frame,'Are you sure you want to change the center?',
                    'Center change',style=wx.OK|wx.CANCEL)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        print 'move center to: ',Xpos,Ypos
                        Data['center'] = [Xpos,Ypos]
                        G2imG.UpdateImageControls(G2frame,Data,Masks)
                finally:
                    dlg.Destroy()
            elif event.key == 'l':
                if G2frame.logPlot:
                    G2frame.logPlot = False
                else:
                    G2frame.logPlot = True
        PlotImage(G2frame,newImage=True)
            
    def OnKeyBox(event):
        if G2frame.G2plotNB.nb.GetSelection() == G2frame.G2plotNB.plotList.index('2D Powder Image'):
            event.key = cb.GetValue()[0]
            cb.SetValue(' key press')
            if event.key in 'l':
                wx.CallAfter(OnImPlotKeyPress,event)
                        
    def OnImPick(event):
        if G2frame.PatternTree.GetItemText(G2frame.PickId) not in ['Image Controls','Masks']:
            return
        if G2frame.setPoly:
            polygon = Masks['Polygons'][-1]
            xpos,ypos = event.mouseevent.xdata,event.mouseevent.ydata
            if xpos and ypos:                       #point inside image
                if len(polygon) > 2 and event.mouseevent.button == 3:
                    x0,y0 = polygon[0]
                    polygon.append([x0,y0])
                    G2frame.setPoly = False
                    G2frame.G2plotNB.status.SetFields(['','Polygon closed - RB drag a vertex to change shape'])
                else:
                    G2frame.G2plotNB.status.SetFields(['','New polygon point: %.1f,%.1f'%(xpos,ypos)])
                    polygon.append([xpos,ypos])
                G2imG.UpdateMasks(G2frame,Masks)
        else:
            if G2frame.itemPicked is not None: return
            G2frame.itemPicked = event.artist
            G2frame.mousePicked = event.mouseevent
        
    def OnImRelease(event):
        try:
            PickName = G2frame.PatternTree.GetItemText(G2frame.PickId)
        except TypeError:
            return
        if PickName not in ['Image Controls','Masks']:
            return
        pixelSize = Data['pixelSize']
        scalex = 1000./pixelSize[0]
        scaley = 1000./pixelSize[1]
        pixLimit = Data['pixLimit']
        if G2frame.itemPicked is None and PickName == 'Image Controls' and len(G2frame.ImageZ):
#            sizexy = Data['size']
            Xpos = event.xdata
            if not (Xpos and G2frame.ifGetRing):                   #got point out of frame
                return
            Ypos = event.ydata
            if Ypos and not Page.toolbar._active:         #make sure zoom/pan not selected
                if event.button == 1:
                    Xpix = Xpos*scalex
                    Ypix = Ypos*scaley
                    xpos,ypos,I,J = G2img.ImageLocalMax(G2frame.ImageZ,pixLimit,Xpix,Ypix)
                    if I and J:
                        xpos += .5                              #shift to pixel center
                        ypos += .5
                        xpos /= scalex                          #convert to mm
                        ypos /= scaley
                        Data['ring'].append([xpos,ypos])
                elif event.button == 3:
                    G2frame.dataFrame.GetStatusBar().SetStatusText('Calibrating...')
                    if G2img.ImageCalibrate(G2frame,Data):
                        G2frame.dataFrame.GetStatusBar().SetStatusText('Calibration successful - Show ring picks to check')
                        print 'Calibration successful'
                    else:
                        G2frame.dataFrame.GetStatusBar().SetStatusText('Calibration failed - Show ring picks to diagnose')
                        print 'Calibration failed'
                    G2frame.ifGetRing = False
                    G2imG.UpdateImageControls(G2frame,Data,Masks)
                    return
                PlotImage(G2frame,newImage=False)
            return
        else:
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                ypos = event.ydata
                if G2frame.ifGetRing:                          #delete a calibration ring pick
                    xypos = [xpos,ypos]
                    rings = Data['ring']
                    for ring in rings:
                        if np.allclose(ring,xypos,.01,0):
                            rings.remove(ring)                                                                       
                else:
                    tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                    itemPicked = str(G2frame.itemPicked)
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
                            
                        G2frame.InnerTth.SetValue("%8.2f" % (Data['IOtth'][0]))
                        G2frame.OuterTth.SetValue("%8.2f" % (Data['IOtth'][1]))
                        G2frame.Lazim.SetValue("%6d" % (Data['LRazimuth'][0]))
                        G2frame.Razim.SetValue("%6d" % (Data['LRazimuth'][1]))
                    elif 'Circle' in itemPicked and PickName == 'Masks':
                        spots = Masks['Points']
                        newPos = itemPicked.split(')')[0].split('(')[2].split(',')
                        newPos = np.array([float(newPos[0]),float(newPos[1])])
                        for spot in spots:
                            if np.allclose(np.array([spot[:2]]),newPos):
                                spot[:2] = xpos,ypos
                        G2imG.UpdateMasks(G2frame,Masks)
                    elif 'Line2D' in itemPicked and PickName == 'Masks':
                        Obj = G2frame.itemPicked.findobj()
                        rings = Masks['Rings']
                        arcs = Masks['Arcs']
                        polygons = Masks['Polygons']
                        for ring in G2frame.ringList:
                            if Obj == ring[0]:
                                rN = ring[1]
                                if ring[2] == 'o':
                                    rings[rN][0] = G2img.GetTth(xpos,ypos,Data)-rings[rN][1]/2.
                                else:
                                    rings[rN][0] = G2img.GetTth(xpos,ypos,Data)+rings[rN][1]/2.
                        for arc in G2frame.arcList:
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
                        for poly in G2frame.polyList:
                            if Obj == poly[0]:
                                ind = G2frame.itemPicked.contains(G2frame.mousePicked)[1]['ind'][0]
                                oldPos = np.array([G2frame.mousePicked.xdata,G2frame.mousePicked.ydata])
                                pN = poly[1]
                                for i,xy in enumerate(polygons[pN]):
                                    if np.allclose(np.array([xy]),oldPos,atol=1.0):
                                        polygons[pN][i] = xpos,ypos
                        G2imG.UpdateMasks(G2frame,Masks)
#                    else:                  #keep for future debugging
#                        print str(G2frame.itemPicked),event.xdata,event.ydata,event.button
                PlotImage(G2frame,newImage=True)
            G2frame.itemPicked = None
            
    try:
        plotNum = G2frame.G2plotNB.plotList.index('2D Powder Image')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        if newImage:
            Page.figure.clf()
            Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('2D Powder Image').gca()
        plotNum = G2frame.G2plotNB.plotList.index('2D Powder Image')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnImPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnImMotion)
        Page.canvas.mpl_connect('pick_event', OnImPick)
        Page.canvas.mpl_connect('button_release_event', OnImRelease)
        xylim = []
    if not event:                       #event from GUI TextCtrl - don't want focus to change to plot!!!
        Page.SetFocus()
    Title = G2frame.PatternTree.GetItemText(G2frame.Image)[4:]
    G2frame.G2plotNB.status.DestroyChildren()
    if G2frame.logPlot:
        Title = 'log('+Title+')'
    Plot.set_title(Title)
    try:
        if G2frame.PatternTree.GetItemText(G2frame.PickId) in ['Image Controls',]:
            Choice = (' key press','l: log(I) on',)
            if G2frame.logPlot:
                Choice[1] = 'l: log(I) off'
            cb = wx.ComboBox(G2frame.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
                choices=Choice)
            cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
            cb.SetValue(' key press')
    except TypeError:
        pass
    size,imagefile = G2frame.PatternTree.GetItemPyData(G2frame.Image)
    if imagefile != G2frame.oldImagefile:
        imagefile = G2IO.CheckImageFile(G2frame,imagefile)
        if not imagefile:
            G2frame.G2plotNB.Delete('2D Powder Image')
            return
        G2frame.PatternTree.SetItemPyData(G2frame.Image,[size,imagefile])
        G2frame.ImageZ = G2IO.GetImageData(G2frame,imagefile,imageOnly=True)
        G2frame.oldImagefile = imagefile

    imScale = 1
    if len(G2frame.ImageZ) > 1024:
        imScale = len(G2frame.ImageZ)/1024
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
            MA = ma.masked_greater(ma.masked_less(G2frame.ImageZ,Zlim[0]),Zlim[1])
            MaskA = ma.getmaskarray(MA)
            A = G2img.ImageCompress(MA,imScale)
            AM = G2img.ImageCompress(MaskA,imScale)
            if G2frame.logPlot:
                A = np.where(A>0,np.log(A),0)
                AM = np.where(AM>0,np.log(AM),0)
                Imin,Imax = [np.amin(A),np.amax(A)]
            ImgM = Plot.imshow(AM,aspect='equal',cmap='Reds',
                interpolation='nearest',vmin=0,vmax=2,extent=[0,Xmax,Ymax,0])
            Img = Plot.imshow(A,aspect='equal',cmap=acolor,
                interpolation='nearest',vmin=Imin,vmax=Imax,extent=[0,Xmax,Ymax,0])
            if G2frame.setPoly:
                Img.set_picker(True)
    
        Plot.plot(xcent,ycent,'x')
        #G2frame.PatternTree.GetItemText(item)
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
                cake = LRAzim[0]+i*delAzm-AzmthOff
                if Data['centerAzm']:
                    cake += delAzm/2.
                ind = np.searchsorted(Azm,cake)
                Plot.plot([arcxI[ind],arcxO[ind]],[arcyI[ind],arcyO[ind]],color='k',dashes=(5,5))
                    
        if G2frame.PatternTree.GetItemText(G2frame.PickId) in 'Image Controls':
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
        if G2frame.PatternTree.GetItemText(G2frame.PickId) in 'Stress/Strain':
            print 'plot stress/strain stuff'
            for ring in StrSta['d-zero']:
                xring,yring = ring['ImxyObs']
                Plot.plot(xring,yring,'r+')
#                for xring,yring in ring['ImxyCalc'].T:
#                    Plot.add_artist(Polygon(ring['ImxyCalc'].T,ec='b',fc='none'))
#                    Plot.plot(xring,yring)
        #masks - mask lines numbered after integration limit lines
        spots = Masks['Points']
        rings = Masks['Rings']
        arcs = Masks['Arcs']
        polygons = Masks['Polygons']
        for x,y,d in spots:
            Plot.add_artist(Circle((x,y),radius=d/2,fc='none',ec='r',picker=3))
        G2frame.ringList = []
        for iring,(tth,thick) in enumerate(rings):
            wave = Data['wavelength']
            x1,y1 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth+thick/2.,wave),Data))),2)
            x2,y2 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth-thick/2.,wave),Data))),2)
            G2frame.ringList.append([Plot.plot(x1,y1,'r',picker=3),iring,'o'])            
            G2frame.ringList.append([Plot.plot(x2,y2,'r',picker=3),iring,'i'])
        G2frame.arcList = []
        for iarc,(tth,azm,thick) in enumerate(arcs):            
            wave = Data['wavelength']
            x1,y1 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(tth+thick/2.,wave),Data),azm)),2)
            x2,y2 = np.hsplit(np.array(G2img.makeIdealRing(G2img.GetEllipse(Dsp(max(0.01,tth-thick/2.),wave),Data),azm)),2)
            G2frame.arcList.append([Plot.plot(x1,y1,'r',picker=3),iarc,'o'])            
            G2frame.arcList.append([Plot.plot(x2,y2,'r',picker=3),iarc,'i'])
            G2frame.arcList.append([Plot.plot([x1[0],x2[0]],[y1[0],y2[0]],'r',picker=3),iarc,'l'])
            G2frame.arcList.append([Plot.plot([x1[-1],x2[-1]],[y1[-1],y2[-1]],'r',picker=3),iarc,'u'])
        G2frame.polyList = []
        for ipoly,polygon in enumerate(polygons):
            x,y = np.hsplit(np.array(polygon),2)
            G2frame.polyList.append([Plot.plot(x,y,'r+',picker=10),ipoly])
            Plot.plot(x,y,'r')            
        if newImage:
            colorBar = Page.figure.colorbar(Img)
        Plot.set_xlim(xlim)
        Plot.set_ylim(ylim)
        Plot.invert_yaxis()
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
        
################################################################################
##### PlotIntegration
################################################################################
            
def PlotIntegration(G2frame,newPlot=False,event=None):
    '''Plot of 2D image after image integration with 2-theta and azimuth as coordinates
    '''
            
    def OnMotion(event):
        Page.canvas.SetToolTipString('')
        Page.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.ydata
        tth = event.xdata
        if azm and tth:
            G2frame.G2plotNB.status.SetFields(\
                ['','Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm)])
                                
    try:
        plotNum = G2frame.G2plotNB.plotList.index('2D Integration')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('2D Integration').gca()
        plotNum = G2frame.G2plotNB.plotList.index('2D Integration')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.views = False
        view = False
    if not event:
        Page.SetFocus()
        
    Data = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
    image = G2frame.Integrate[0]
    xsc = G2frame.Integrate[1]
    ysc = G2frame.Integrate[2]
    Imin,Imax = Data['range'][1]
    acolor = mpl.cm.get_cmap(Data['color'])
    Plot.set_title(G2frame.PatternTree.GetItemText(G2frame.Image)[4:])
    Plot.set_ylabel('azimuth',fontsize=12)
    Plot.set_xlabel('2-theta',fontsize=12)
    Img = Plot.imshow(image,cmap=acolor,vmin=Imin,vmax=Imax,interpolation='nearest', \
        extent=[ysc[0],ysc[-1],xsc[-1],xsc[0]],aspect='auto')
    colorBar = Page.figure.colorbar(Img)
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2img.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2img.GetTthAzm(x,y,Data)
#            azm = np.where(azm < 0.,azm+360,azm)
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
                
################################################################################
##### PlotTRImage
################################################################################
            
def PlotTRImage(G2frame,tax,tay,taz,newPlot=False):
    '''a test plot routine - not normally used
    ''' 
            
    def OnMotion(event):
        Page.canvas.SetToolTipString('')
        Page.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.xdata
        tth = event.ydata
        if azm and tth:
            G2frame.G2plotNB.status.SetFields(\
                ['','Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm)])
                                
    try:
        plotNum = G2frame.G2plotNB.plotList.index('2D Transformed Powder Image')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError:
        Plot = G2frame.G2plotNB.addMpl('2D Transformed Powder Image').gca()
        plotNum = G2frame.G2plotNB.plotList.index('2D Transformed Powder Image')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.views = False
        view = False
    Page.SetFocus()
        
    Data = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.Image, 'Image Controls'))
    Imin,Imax = Data['range'][1]
    step = (Imax-Imin)/5.
    V = np.arange(Imin,Imax,step)
    acolor = mpl.cm.get_cmap(Data['color'])
    Plot.set_title(G2frame.PatternTree.GetItemText(G2frame.Image)[4:])
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
        
################################################################################
##### PlotStructure
################################################################################
            
def PlotStructure(G2frame,data):
    '''Crystal structure plotting package. Can show structures as balls, sticks, lines,
    thermal motion ellipsoids and polyhedra
    '''
    ForthirdPI = 4.0*math.pi/3.0
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    Vol = generalData['Cell'][7:8][0]
    Amat,Bmat = G2lat.cell2AB(cell)         #Amat - crystal to cartesian, Bmat - inverse
    Gmat,gmat = G2lat.cell2Gmat(cell)
    A4mat = np.concatenate((np.concatenate((Amat,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
    B4mat = np.concatenate((np.concatenate((Bmat,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
    Mydir = generalData['Mydir']
    atomData = data['Atoms']
    mapPeaks = []
    if 'Map Peaks' in data:
        mapPeaks = data['Map Peaks']
    drawingData = data['Drawing']
    try:
        drawAtoms = drawingData['Atoms']
    except KeyError:
        drawAtoms = []
    mapData = {}
    flipData = {}
    rhoXYZ = []
    if 'Map' in generalData:
        mapData = generalData['Map']
    if 'Flip' in generalData:
        flipData = generalData['Flip']                        
        flipData['mapRoll'] = [0,0,0]
    cx,ct,cs,ci = drawingData['atomPtrs']
    Wt = np.array([255,255,255])
    Rd = np.array([255,0,0])
    Gr = np.array([0,255,0])
    Bl = np.array([0,0,255])
    uBox = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    uEdges = np.array([
        [uBox[0],uBox[1]],[uBox[0],uBox[3]],[uBox[0],uBox[4]],[uBox[1],uBox[2]], 
        [uBox[2],uBox[3]],[uBox[1],uBox[5]],[uBox[2],uBox[6]],[uBox[3],uBox[7]], 
        [uBox[4],uBox[5]],[uBox[5],uBox[6]],[uBox[6],uBox[7]],[uBox[7],uBox[4]]])
    mD = 0.1
    mV = np.array([[[-mD,0,0],[mD,0,0]],[[0,-mD,0],[0,mD,0]],[[0,0,-mD],[0,0,mD]]])
    mapPeakVecs = np.inner(mV,Bmat)

    uColors = [Rd,Gr,Bl,Wt, Wt,Wt,Wt,Wt, Wt,Wt,Wt,Wt]
    altDown = False
    shiftDown = False
    ctrlDown = False
    
    def OnKeyBox(event):
        import Image
#        Draw()                          #make sure plot is fresh!!
        mode = cb.GetValue()
        if mode in ['jpeg','bmp','tiff',]:
            Fname = os.path.joint(Mydir,generalData['Name']+'.'+mode)
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
            cb.SetValue(' save as/key:')
            G2frame.G2plotNB.status.SetStatusText('Drawing saved to: '+Fname,1)
        else:
            event.key = cb.GetValue()[0]
            cb.SetValue(' save as/key:')
            wx.CallAfter(OnKey,event)

    def OnKey(event):           #on key UP!!
#        Draw()                          #make sure plot is fresh!!
        try:
            keyCode = event.GetKeyCode()
            if keyCode > 255:
                keyCode = 0
            key = chr(keyCode)
        except AttributeError:       #if from OnKeyBox above
            key = str(event.key).upper()
        indx = drawingData['selectedAtoms']
        if key in ['C']:
            drawingData['viewPoint'] = [[.5,.5,.5],[0,0]]
            drawingData['viewDir'] = [0,0,1]
            drawingData['oldxy'] = []
            V0 = np.array([0,0,1])
            V = np.inner(Amat,V0)
            V /= np.sqrt(np.sum(V**2))
            A = np.arccos(np.sum(V*V0))
            Q = G2mth.AV2Q(A,[0,1,0])
            drawingData['Quaternion'] = Q
            SetViewPointText(drawingData['viewPoint'][0])
            SetViewDirText(drawingData['viewDir'])
            Q = drawingData['Quaternion']
            G2frame.G2plotNB.status.SetStatusText('New quaternion: %.2f+, %.2fi+ ,%.2fj+, %.2fk'%(Q[0],Q[1],Q[2],Q[3]),1)
        elif key in ['N']:
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
            G2frame.G2plotNB.status.SetStatusText('View point at atom '+drawAtoms[pI[0]][ct-1]+str(pI),1)
                
        elif key in ['P']:
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
            G2frame.G2plotNB.status.SetStatusText('View point at atom '+drawAtoms[pI[0]][ct-1]+str(pI),1)
        elif key in ['U','D','L','R'] and mapData['Flip'] == True:
            dirDict = {'U':[0,1],'D':[0,-1],'L':[-1,0],'R':[1,0]}
            SetMapRoll(dirDict[key])
            SetPeakRoll(dirDict[key])
            SetMapPeaksText(mapPeaks)
        Draw()
            
    def GetTruePosition(xy,Add=False):
        View = glGetIntegerv(GL_VIEWPORT)
        Proj = glGetDoublev(GL_PROJECTION_MATRIX)
        Model = glGetDoublev(GL_MODELVIEW_MATRIX)
        Zmax = 1.
        if Add:
            Indx = GetSelectedAtoms()
        if G2frame.dataDisplay.GetPageText(getSelection()) == 'Map peaks':
            for i,peak in enumerate(mapPeaks):
                x,y,z = peak[1:4]
                X,Y,Z = gluProject(x,y,z,Model,Proj,View)
                XY = [int(X),int(View[3]-Y)]
                if np.allclose(xy,XY,atol=10) and Z < Zmax:
                    Zmax = Z
                    try:
                        Indx.remove(i)
                        ClearSelectedAtoms()
                        for id in Indx:
                            SetSelectedAtoms(id,Add)
                    except:
                        SetSelectedAtoms(i,Add)
        else:
            for i,atom in enumerate(drawAtoms):
                x,y,z = atom[cx:cx+3]
                X,Y,Z = gluProject(x,y,z,Model,Proj,View)
                XY = [int(X),int(View[3]-Y)]
                if np.allclose(xy,XY,atol=10) and Z < Zmax:
                    Zmax = Z
                    try:
                        Indx.remove(i)
                        ClearSelectedAtoms()
                        for id in Indx:
                            SetSelectedAtoms(id,Add)
                    except:
                        SetSelectedAtoms(i,Add)
                                       
    def OnMouseDown(event):
        xy = event.GetPosition()
        if event.ShiftDown():
            if event.LeftIsDown():
                GetTruePosition(xy)
            elif event.RightIsDown():
                GetTruePosition(xy,True)
        else:
            drawingData['oldxy'] = list(xy)
#        Draw()
        
    def OnMouseMove(event):
        if event.ShiftDown():           #don't want any inadvertant moves when picking
            return
        newxy = event.GetPosition()
        page = getSelection()
                                
        if event.Dragging() and not event.ControlDown():
            if event.LeftIsDown():
                SetRotation(newxy)
                Q = drawingData['Quaternion']
                G2frame.G2plotNB.status.SetStatusText('New quaternion: %.2f+, %.2fi+ ,%.2fj+, %.2fk'%(Q[0],Q[1],Q[2],Q[3]),1)
            elif event.RightIsDown():
                SetTranslation(newxy)
                Tx,Ty,Tz = drawingData['viewPoint'][0]
                G2frame.G2plotNB.status.SetStatusText('New view point: %.4f, %.4f, %.4f'%(Tx,Ty,Tz),1)
            elif event.MiddleIsDown():
                SetRotationZ(newxy)
                Q = drawingData['Quaternion']
                G2frame.G2plotNB.status.SetStatusText('New quaternion: %.2f+, %.2fi+ ,%.2fj+, %.2fk'%(Q[0],Q[1],Q[2],Q[3]),1)
            Draw()
        
    def OnMouseWheel(event):
        if event.ShiftDown():
            return
        drawingData['cameraPos'] += event.GetWheelRotation()/24
        drawingData['cameraPos'] = max(10,min(500,drawingData['cameraPos']))
        G2frame.G2plotNB.status.SetStatusText('New camera distance: %.2f'%(drawingData['cameraPos']),1)
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Options':
                panel = G2frame.dataDisplay.GetPage(page).GetChildren()[0].GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('cameraPos')].SetLabel('Camera Position: '+'%.2f'%(drawingData['cameraPos']))
                panel[names.index('cameraSlider')].SetValue(drawingData['cameraPos'])
            Draw()
        
    def getSelection():
        try:
            return G2frame.dataDisplay.GetSelection()
        except AttributeError:
            G2frame.G2plotNB.status.SetStatusText('Select this from Phase data window!')
            return 0
            
    def SetViewPointText(VP):
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Options':
                panel = G2frame.dataDisplay.GetPage(page).GetChildren()[0].GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('viewPoint')].SetValue('%.3f %.3f %.3f'%(VP[0],VP[1],VP[2]))
                
    def SetViewDirText(VD):
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Options':
                panel = G2frame.dataDisplay.GetPage(page).GetChildren()[0].GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('viewDir')].SetValue('%.3f %.3f %.3f'%(VD[0],VD[1],VD[2]))
                
    def SetMapPeaksText(mapPeaks):
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Map peaks':
                G2frame.MapPeaksTable.SetData(mapPeaks)
                panel = G2frame.dataDisplay.GetPage(page).GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('grid window')].Refresh()
            
    def ClearSelectedAtoms():
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Atoms':
                G2frame.dataDisplay.GetPage(page).ClearSelection()      #this is the Atoms grid in Draw Atoms
            elif G2frame.dataDisplay.GetPageText(page) == 'Map peaks':
                G2frame.dataDisplay.GetPage(page).ClearSelection()      #this is the Atoms grid in Atoms
            elif G2frame.dataDisplay.GetPageText(page) == 'Atoms':
                G2frame.dataDisplay.GetPage(page).ClearSelection()      #this is the Atoms grid in Atoms
                
                    
    def SetSelectedAtoms(ind,Add=False):
        page = getSelection()
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Atoms':
                G2frame.dataDisplay.GetPage(page).SelectRow(ind,Add)      #this is the Atoms grid in Draw Atoms
            elif G2frame.dataDisplay.GetPageText(page) == 'Map peaks':
                G2frame.dataDisplay.GetPage(page).SelectRow(ind,Add)                  
            elif G2frame.dataDisplay.GetPageText(page) == 'Atoms':
                Id = drawAtoms[ind][-2]
                for i,atom in enumerate(atomData):
                    if atom[-1] == Id:
                        G2frame.dataDisplay.GetPage(page).SelectRow(i)      #this is the Atoms grid in Atoms
                  
    def GetSelectedAtoms():
        page = getSelection()
        Ind = []
        if page:
            if G2frame.dataDisplay.GetPageText(page) == 'Draw Atoms':
                Ind = G2frame.dataDisplay.GetPage(page).GetSelectedRows()      #this is the Atoms grid in Draw Atoms
            elif G2frame.dataDisplay.GetPageText(page) == 'Map peaks':
                Ind = G2frame.dataDisplay.GetPage(page).GetSelectedRows()
            elif G2frame.dataDisplay.GetPageText(page) == 'Atoms':
                Ind = G2frame.dataDisplay.GetPage(page).GetSelectedRows()      #this is the Atoms grid in Atoms
        return Ind
                                       
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
        
    def GetRoll(newxy,rho):
        Q = drawingData['Quaternion']
        dxy = G2mth.prodQVQ(G2mth.invQ(Q),np.inner(Bmat,newxy+[0,]))
        dxy = np.array(dxy*rho.shape)        
        roll = np.where(dxy>0.5,1,np.where(dxy<-.5,-1,0))
        return roll
                
    def SetMapRoll(newxy):
        rho = mapData['rho']
        roll = GetRoll(newxy,rho)
        mapData['rho'] = np.roll(np.roll(np.roll(rho,roll[0],axis=0),roll[1],axis=1),roll[2],axis=2)
        drawingData['oldxy'] = list(newxy)
        
    def SetPeakRoll(newxy):
        rho = mapData['rho']
        roll = GetRoll(newxy,rho)
        steps = 1./np.array(rho.shape)
        dxy = roll*steps
        for peak in mapPeaks:
            peak[1:4] += dxy
            peak[1:4] %= 1.
            peak[4] = np.sqrt(np.sum(np.inner(Amat,peak[1:4])**2))
                
    def SetTranslation(newxy):
#first get translation vector in screen coords.       
        oldxy = drawingData['oldxy']
        if not len(oldxy): oldxy = list(newxy)
        dxy = newxy-oldxy
        drawingData['oldxy'] = list(newxy)
        V = np.array([-dxy[0],dxy[1],0.])
#then transform to rotated crystal coordinates & apply to view point        
        Q = drawingData['Quaternion']
        V = np.inner(Bmat,G2mth.prodQVQ(G2mth.invQ(Q),V))
        Tx,Ty,Tz = drawingData['viewPoint'][0]
        Tx += V[0]*0.01
        Ty += V[1]*0.01
        Tz += V[2]*0.01
        drawingData['oldxy'] = list(newxy)
        drawingData['viewPoint'][0] =  Tx,Ty,Tz
        SetViewPointText([Tx,Ty,Tz])
        
    def SetRotation(newxy):
#first get rotation vector in screen coords. & angle increment        
        oldxy = drawingData['oldxy']
        if not len(oldxy): oldxy = list(newxy)
        dxy = newxy-oldxy
        drawingData['oldxy'] = list(newxy)
        V = np.array([dxy[1],dxy[0],0.])
        A = 0.25*np.sqrt(dxy[0]**2+dxy[1]**2)
# next transform vector back to xtal coordinates via inverse quaternion
# & make new quaternion
        Q = drawingData['Quaternion']
        V = G2mth.prodQVQ(G2mth.invQ(Q),np.inner(Bmat,V))
        DQ = G2mth.AVdeg2Q(A,V)
        Q = G2mth.prodQQ(Q,DQ)
        drawingData['Quaternion'] = Q
# finally get new view vector - last row of rotation matrix
        VD = np.inner(Bmat,G2mth.Q2Mat(Q)[2])
        VD /= np.sqrt(np.sum(VD**2))
        drawingData['viewDir'] = VD
        SetViewDirText(VD)
        
    def SetRotationZ(newxy):                        
#first get rotation vector (= view vector) in screen coords. & angle increment        
        View = glGetIntegerv(GL_VIEWPORT)
        cent = [View[2]/2,View[3]/2]
        oldxy = drawingData['oldxy']
        if not len(oldxy): oldxy = list(newxy)
        dxy = newxy-oldxy
        drawingData['oldxy'] = list(newxy)
        V = drawingData['viewDir']
        A = [0,0]
        A[0] = dxy[1]*.25
        A[1] = dxy[0]*.25
        if newxy[0] > cent[0]:
            A[0] *= -1
        if newxy[1] < cent[1]:
            A[1] *= -1        
# next transform vector back to xtal coordinates & make new quaternion
        Q = drawingData['Quaternion']
        V = np.inner(Amat,V)
        Qx = G2mth.AVdeg2Q(A[0],V)
        Qy = G2mth.AVdeg2Q(A[1],V)
        Q = G2mth.prodQQ(Q,Qx)
        Q = G2mth.prodQQ(Q,Qy)
        drawingData['Quaternion'] = Q

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
            glVertex3fv(-line[1]/2.)
            glVertex3fv(line[1]/2.)
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
        
    def RenderDots(XYZ,RC):
        glEnable(GL_COLOR_MATERIAL)
        XYZ = np.array(XYZ)
        glPushMatrix()
        for xyz,rc in zip(XYZ,RC):
            x,y,z = xyz
            r,c = rc
            glColor3ubv(c)
            glPointSize(r*50)
            glBegin(GL_POINTS)
            glVertex3fv(xyz)
            glEnd()
        glPopMatrix()
        glColor4ubv([0,0,0,0])
        glDisable(GL_COLOR_MATERIAL)
        
    def RenderSmallSphere(x,y,z,radius,color):
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glPushMatrix()
        glTranslate(x,y,z)
        glMultMatrixf(B4mat.T)
        q = gluNewQuadric()
        gluSphere(q,radius,4,2)
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
            if Z:
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

    def RenderMapPeak(x,y,z,color):
        xyz = np.array([x,y,z])
        glEnable(GL_COLOR_MATERIAL)
        glLineWidth(3)
        glColor3fv(color)
        glPushMatrix()
        glBegin(GL_LINES)
        for vec in mapPeakVecs:
            glVertex3fv(vec[0]+xyz)
            glVertex3fv(vec[1]+xyz)
        glEnd()
        glColor4ubv([0,0,0,0])
        glPopMatrix()
        glDisable(GL_COLOR_MATERIAL)
        
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
        
    def RenderMap(rho,rhoXYZ,indx,Rok):
        cLevel = drawingData['contourLevel']
        XYZ = []
        RC = []
        for i,xyz in enumerate(rhoXYZ):
            if not Rok[i]:
                x,y,z = xyz
                I,J,K = indx[i]
                alpha = 1.0
                if cLevel < 1.:
                    alpha = (abs(rho[I,J,K])/mapData['rhoMax']-cLevel)/(1.-cLevel)
                if rho[I,J,K] < 0.:
                    XYZ.append(xyz)
                    RC.append([0.1*alpha,Rd])
                else:
                    XYZ.append(xyz)
                    RC.append([0.1*alpha,Gr])
        RenderDots(XYZ,RC)
                            
    def Draw():
        mapData = generalData['Map']
        pageName = ''
        page = getSelection()
        if page:
            pageName = G2frame.dataDisplay.GetPageText(page)
        rhoXYZ = []
        if len(mapData['rho']):
            VP = np.array(drawingData['viewPoint'][0])-np.array([.5,.5,.5])
            contLevel = drawingData['contourLevel']*mapData['rhoMax']
            if 'delt-F' in mapData['MapType']:
                rho = ma.array(mapData['rho'],mask=(np.abs(mapData['rho'])<contLevel))
            else:
                rho = ma.array(mapData['rho'],mask=(mapData['rho']<contLevel))
            steps = 1./np.array(rho.shape)
            incre = np.where(VP>=0,VP%steps,VP%steps-steps)
            Vsteps = -np.array(VP/steps,dtype='i')
            rho = np.roll(np.roll(np.roll(rho,Vsteps[0],axis=0),Vsteps[1],axis=1),Vsteps[2],axis=2)
            indx = np.array(ma.nonzero(rho)).T
            rhoXYZ = indx*steps+VP-incre
            Nc = len(rhoXYZ)
            rcube = 2000.*Vol/(ForthirdPI*Nc)
            rmax = math.exp(math.log(rcube)/3.)**2
            radius = min(drawingData['mapSize']**2,rmax)
            view = np.array(drawingData['viewPoint'][0])
            Rok = np.sum(np.inner(Amat,rhoXYZ-view).T**2,axis=1)>radius
        Ind = GetSelectedAtoms()
        VS = np.array(Page.canvas.GetSize())
        aspect = float(VS[0])/float(VS[1])
        cPos = drawingData['cameraPos']
        Zclip = drawingData['Zclip']*cPos/200.
        Q = drawingData['Quaternion']
        Tx,Ty,Tz = drawingData['viewPoint'][0]
        cx,ct,cs,ci = drawingData['atomPtrs']
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
        matRot = G2mth.Q2Mat(Q)
        matRot = np.concatenate((np.concatenate((matRot,[[0],[0],[0]]),axis=1),[[0,0,0,1],]),axis=0)
        glMultMatrixf(matRot.T)
        glMultMatrixf(A4mat.T)
        glTranslate(-Tx,-Ty,-Tz)
        if drawingData['unitCellBox']:
            RenderBox()
        if drawingData['showABC']:
            x,y,z = drawingData['viewPoint'][0]
            RenderUnitVectors(x,y,z)
        Backbones = {}
        BackboneColor = []
        time0 = time.time()
#        glEnable(GL_BLEND)
#        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
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
            if iat in Ind and G2frame.dataDisplay.GetPageText(getSelection()) != 'Map peaks':
                color = np.array(Gr)/255.
#            color += [.25,]
            radius = 0.5
            if atom[cs] != '':
                try:
                    glLoadName(atom[-3])
                except: #problem with old files - missing code
                    pass                    
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
                    if atom[2] not in Backbones:
                        Backbones[atom[2]] = []
                    Backbones[atom[2]].append(list(np.inner(Amat,np.array([x,y,z]))))
                    BackboneColor.append(list(color))
                    
            if atom[cs+1] == 'type':
                RenderLabel(x,y,z,atom[ct],radius)
            elif atom[cs+1] == 'name':
                RenderLabel(x,y,z,atom[ct-1],radius)
            elif atom[cs+1] == 'number':
                RenderLabel(x,y,z,str(iat),radius)
            elif atom[cs+1] == 'residue' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-4],radius)
            elif atom[cs+1] == '1-letter' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-3],radius)
            elif atom[cs+1] == 'chain' and atom[ct-1] == 'CA':
                RenderLabel(x,y,z,atom[ct-2],radius)
#        glDisable(GL_BLEND)
        if len(rhoXYZ):
            RenderMap(rho,rhoXYZ,indx,Rok)
        if len(mapPeaks):
            for ind,[mag,x,y,z,d] in enumerate(mapPeaks):
                if ind in Ind and pageName == 'Map peaks':
                    RenderMapPeak(x,y,z,Gr)
                else:
                    RenderMapPeak(x,y,z,Wt)
        if Backbones:
            for chain in Backbones:
                Backbone = Backbones[chain]
                RenderBackbone(Backbone,BackboneColor,bondR)
#        print time.time()-time0
        Page.canvas.SwapBuffers()
       
    def OnSize(event):
        Draw()
        
    def OnFocus(event):
        Draw()
        
    try:
        plotNum = G2frame.G2plotNB.plotList.index(generalData['Name'])
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)        
    except ValueError:
        Plot = G2frame.G2plotNB.addOgl(generalData['Name'])
        plotNum = G2frame.G2plotNB.plotList.index(generalData['Name'])
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        Page.views = False
        view = False
        altDown = False
    Page.SetFocus()
    if mapData['Flip']:
        choice = [' save as/key:','jpeg','tiff','bmp','u: roll up','d: roll down','l: roll left','r: roll right']
    else:
        choice = [' save as/key:','jpeg','tiff','bmp','c: center on 1/2,1/2,1/2','n: next','p: previous']
    cb = wx.ComboBox(G2frame.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,choices=choice)
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' save as/key:')
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
        
