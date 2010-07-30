import math
import time
import copy
import numpy as np
import wx
import wx.aui
import matplotlib as mpl
import GSASIIpath
import GSASIIgrid as G2gd
import GSASIIimage as G2img
import GSASIIIO as G2IO
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
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
    
class G2Plot(wx.Panel):
    
    def __init__(self,parent,id=-1,dpi=None,**kwargs):
        wx.Panel.__init__(self,parent,id=id,**kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi,figsize=(5,7))
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
        self.status.SetStatusWidths([125,-1])
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        
        self.plotList = []
            
    def add(self,name=""):
        page = G2Plot(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
    def clear(self):
        while self.nb.GetPageCount():
            self.nb.DeletePage(0)
        self.plotList = []
        self.status.DestroyChildren()
        
    def OnPageChanged(self,event):
        self.status.DestroyChildren()                           #get rid of special stuff on status bar
        
def PlotSngl(self,newPlot=False):
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
    except ValueError,error:
        Plot = self.G2plotNB.add('Structure Factors').gca()
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
    global HKL
    
    def OnPick(event):
        if self.itemPicked is not None: return
        PatternId = self.PatternId
        PickId = self.PickId
        pick = event.artist
        mouse = event.mouseevent
        xpos = pick.get_xdata()
        ypos = pick.get_ydata()
        ind = event.ind
        view = Page.toolbar._views.forward()
        if view and 'line2' in str(pick):           #apply offset only for picked powder pattern points
            ind += np.searchsorted(xye[0],view[0][0])
        xy = zip(xpos[ind],ypos[ind])[0]
        if self.PatternTree.GetItemText(PickId) == 'Peak List':
            if ind.all() != [0]:                                    #picked a data point
                inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
                if len(inst[1]) == 11:
                    ins = inst[1][4:10]
                else:
                    ins = inst[1][6:12]    
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
        if event.key == 'w':
            if self.Weight:
                self.Weight = False
            else:
                self.Weight = True
            print 'plot weighting:',self.Weight
        elif event.key == 'u' and self.Offset < 100.:
            self.Offset += 1.
        elif event.key == 'd' and self.Offset > 0.:
            self.Offset -= 1.
        elif event.key == 'c':
            if self.Contour:
                self.Contour = False
            else:
                self.Contour = True
        elif event.key == 's':
            if self.SinglePlot:
                self.SinglePlot = False
            else:
                self.SinglePlot = True
            
        PlotPatterns(self)
        
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
                wave = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self, \
                    self.PatternId, 'Instrument Parameters'))[0][1]
                dsp = 0.0
                if abs(xpos) > 0.:                  #avoid possible singularity at beam center
                    dsp = wave/(2.*sind(abs(xpos)/2.0))
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
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            lines = []
            for line in self.Lines: lines.append(line.get_xdata()[0])
            lineNo = lines.index(self.itemPicked.get_xdata()[0])
            if  lineNo in [0,1]:
                LimitId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Limits')
                data = self.PatternTree.GetItemPyData(LimitId)
#                print 'limits',xpos
                data[1][lineNo] = xpos
                self.PatternTree.SetItemPyData(LimitId,data)
                if self.PatternTree.GetItemText(self.PickId) == 'Limits':
                    G2pdG.UpdateLimitsGrid(self,data)
            else:
                PeakId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List')
                data = self.PatternTree.GetItemPyData(PeakId)
#                print 'peaks',xpos
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
    except ValueError,error:
        newPlot = True
        Plot = self.G2plotNB.add('Powder Patterns').gca()
        plotNum = self.G2plotNB.plotList.index('Powder Patterns')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.canvas.mpl_connect('pick_event', OnPick)
        Page.canvas.mpl_connect('button_release_event', OnRelease)
    Page.SetFocus()
    cb = wx.ComboBox(self.G2plotNB.status,style=wx.CB_DROPDOWN|wx.CB_READONLY,
        choices=(' key press','d: offset down','u: offset up','c: toggle contour','s: toggle single plot'))
    cb.Bind(wx.EVT_COMBOBOX, OnKeyBox)
    cb.SetValue(' key press')
    
    PickId = self.PickId
    PatternId = self.PatternId
    colors=['b','g','r','c','m','k']
    Lines = []
    if self.SinglePlot:
        Pattern = self.PatternTree.GetItemPyData(self.PatternId)
        Pattern.append(self.PatternTree.GetItemText(self.PatternId))
        PlotList = [Pattern,]
    else:        
        PlotList = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            if 'PWDR' in self.PatternTree.GetItemText(item):
                Pattern = self.PatternTree.GetItemPyData(item)
                Pattern.append(self.PatternTree.GetItemText(item))
                PlotList.append(Pattern)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
    Ymax = 1.0
    HKL = np.array(self.HKL)
    for Pattern in PlotList:
        xye = Pattern[1]
        Ymax = max(Ymax,max(xye[1]))
    offset = self.Offset*Ymax/100.0
    Plot.set_title('Powder Patterns')
    Plot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    Plot.set_ylabel('Intensity',fontsize=12)
    if self.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    for N,Pattern in enumerate(PlotList):
        ifpicked = False
        LimitId = 0
        xye = Pattern[1]
        if PickId:
            ifpicked = Pattern[2] == self.PatternTree.GetItemText(PatternId)
            LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
        X = xye[0]
        Y = xye[1]+offset*N
        if LimitId:
            limits = self.PatternTree.GetItemPyData(LimitId)
            Lines.append(Plot.axvline(limits[1][0],color='g',dashes=(5,5),picker=3.))    
            Lines.append(Plot.axvline(limits[1][1],color='r',dashes=(5,5),picker=3.))                    
        if self.Contour:
            ContourY.append(N)
            ContourZ.append(Y)
            ContourX = X
            Nseq += 1
            Plot.set_ylabel('Data sequence',fontsize=12)
        else:
            if ifpicked:
                Z = xye[3]+offset*N
                W = xye[4]+offset*N
                D = xye[5]+offset*N-Ymax*.02
                if self.Weight:
                    W2 = np.sqrt(xye[2])
                    D *= W2-Ymax*.02
                Plot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                Plot.plot(X,Z,colors[(N+1)%6],picker=False)
                Plot.plot(X,W,colors[(N+2)%6],picker=False)
                Plot.plot(X,D,colors[(N+3)%6],picker=False)
                Plot.axhline(0.,color=wx.BLACK)
                Page.canvas.SetToolTipString('')
                if self.PatternTree.GetItemText(PickId) == 'Peak List':
                    tip = 'On data point: Pick peak - L or R MB.On line: MB down to move'
                    Page.canvas.SetToolTipString(tip)
                    data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
                    for item in data:
                        Lines.append(Plot.axvline(item[0],color=colors[N%6],picker=2.))
                if self.PatternTree.GetItemText(PickId) == 'Limits':
                    tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                    Page.canvas.SetToolTipString(tip)
                    data = self.LimitsTable.GetData()
            else:
                Plot.plot(X,Y,colors[N%6],picker=False)
    if PickId and self.PatternTree.GetItemText(PickId) in ['Index Peak List','Unit Cells List']:
        peaks = np.array((self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))))
        for peak in peaks:
            Plot.axvline(peak[0],color='b')
        for hkl in self.HKL:
            Plot.axvline(hkl[5],color='r',dashes=(5,5))
    if self.Contour:
        acolor = mpl.cm.get_cmap('Paired')
        Plot.imshow(ContourZ,cmap=acolor,vmin=0,vmax=Ymax,interpolation='nearest', 
            extent=[ContourX[0],ContourX[-1],ContourY[0],ContourY[-1]],aspect='auto')
        newPlot = True
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

def PlotPowderLines(self):
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
    except ValueError,error:
        Plot = self.G2plotNB.add('Powder Lines').gca()
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

def PlotPeakWidths(self):
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
    except ValueError,error:
        Plot = self.G2plotNB.add('Peak Widths').gca()
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
            
def PlotExposedImage(self,newPlot=False,event=None):
    plotNo = self.G2plotNB.nb.GetSelection()
    if self.G2plotNB.nb.GetPageText(plotNo) == '2D Powder Image':
        PlotImage(self,newPlot,event)
    elif self.G2plotNB.nb.GetPageText(plotNo) == '2D Integration':
        PlotIntegration(self,newPlot,event)

def PlotImage(self,newPlot=False,event=None):
    from matplotlib.patches import Ellipse,Arc,Circle,Polygon
    import numpy.ma as ma
    Dsp = lambda tth,wave: wave/(2.*sind(tth/2.))

    def OnImMotion(event):
        Data = self.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
        Page.canvas.SetToolTipString('')
        size = len(self.ImageZ)
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
                    if 'line3' in  str(item) or 'line4' in str(item) and not Data['fullIntegrate']:
                        ang = int(atan2d(xpos,ypos))
                        Page.canvas.SetToolTipString('%6d deg'%(ang))
                    elif 'line1' in  str(item) or 'line2' in str(item):
                        tth = G2img.GetTth(event.xdata,event.ydata,Data)
                        Page.canvas.SetToolTipString('%8.3fdeg'%(tth))                           
            else:
                xpos = event.xdata
                ypos = event.ydata
                xpix = xpos*scalex
                ypix = ypos*scaley
                Int = 0
                if (0 <= xpix <= size) and (0 <= ypix <= size):
                    Int = self.ImageZ[ypix][xpix]
                tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                Q = 2.*math.pi/dsp
                self.G2plotNB.status.SetFields(\
                    ['','Detector 2-th =%9.2fdeg, dsp =%9.3fA, Q = %6.3fA-1, azm = %7.2fdeg, I = %6d'%(tth,dsp,Q,azm,Int)])

    def OnImPlotKeyPress(event):
        if self.PatternTree.GetItemText(self.PickId) == 'Masks':
            Data = self.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
            Masks = self.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))
            Xpos = event.xdata
            if not Xpos:            #got point out of frame
                return
            Ypos = event.ydata
            if event.key == 's':
                print 'spot mask @ ',Xpos,Ypos
                Masks['Points'].append([Xpos,Ypos,1.])
            elif event.key == 'r':
                tth = G2img.GetTth(Xpos,Ypos,Data)
                print 'ring mask @ ',Xpos,Ypos,tth
                Masks['Rings'].append([tth,0.1])
            elif event.key == 'a':
                tth,azm = G2img.GetTthAzm(Xpos,Ypos,Data)
                azm = int(azm)                
                print 'arc mask @ ', Xpos,Ypos
                Masks['Arcs'].append([tth,[azm-5,azm+5],0.1])
            elif event.key == 'p':
                self.setPoly = True
                Masks['Polygons'].append([])
                print 'Polygon mask active - pick points with mouse LB'
                print '   use RB to close when > 2 points chosen'
                print 'Vertices can be dragged with LB down after polygon closed'
            G2imG.UpdateMasks(self,Masks)
        PlotImage(self)
            
    def OnImPick(event):
        if self.PatternTree.GetItemText(self.PickId) not in ['Image Controls','Masks']:
            return
        if self.setPoly:
            Masks = self.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))
            polygon = Masks['Polygons'][-1]
            xpos,ypos = event.mouseevent.xdata,event.mouseevent.ydata
            if xpos and ypos:                       #point inside image
                if len(polygon) > 2 and event.mouseevent.button == 3:
                    x0,y0 = polygon[0]
                    polygon.append([x0,y0])
                    self.setPoly = False
                else:           
                    polygon.append([xpos,ypos])
                G2imG.UpdateMasks(self,Masks)
        else:
            if self.itemPicked is not None: return
            self.itemPicked = event.artist
            self.mousePicked = event.mouseevent
        
    def OnImRelease(event):
        PickName = self.PatternTree.GetItemText(self.PickId)
        if PickName not in ['Image Controls','Masks']:
            return
        Data = self.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
        Masks = self.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))
        pixelSize = Data['pixelSize']
        scalex = 1000./pixelSize[0]
        scaley = 1000./pixelSize[1]
        if self.itemPicked is None and PickName == 'Image Controls':
            size = len(self.ImageZ)
            Xpos = event.xdata
            if not (Xpos and self.ifGetRing):                   #got point out of frame
                return
            Ypos = event.ydata
            if Ypos and not Page.toolbar._active:         #make sure zoom/pan not selected
                if event.button == 1:
                    Xpix = Xpos*scalex
                    Ypix = Ypos*scaley
                    xpos,ypos,I,J = G2img.ImageLocalMax(self.ImageZ,20,Xpix,Ypix)
                    if I and J:
                        xpos += .5                              #shift to pixel center
                        ypos += .5
                        xpos /= scalex                          #convert to mm
                        ypos /= scaley
                        Data['ring'].append([xpos,ypos])
                PlotImage(self)
            return
        else:
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                ypos = event.ydata
                if self.ifGetRing:
                    xypos = [xpos,ypos]
                    rings = Data['ring']
                    for ring in rings:
                        if np.allclose(ring,xypos,.01,0):
                            rings.remove(ring)                                                                       
                else:
                    tth,azm,dsp = G2img.GetTthAzmDsp(xpos,ypos,Data)
                    itemPicked = str(self.itemPicked)
                    if 'Line2D' in itemPicked and PickName == 'Image Controls':
                        print int(itemPicked.split('_line')[1].strip(')'))
                        if 'line1' in itemPicked:
                            Data['IOtth'][0] = tth
                        elif 'line2' in itemPicked:
                            Data['IOtth'][1] = tth
                        elif 'line3' in itemPicked and not Data['fullIntegrate']:
                            Data['LRazimuth'][0] = int(azm)
                        elif 'line4' in itemPicked and not Data['fullIntegrate']:
                            Data['LRazimuth'][1] = int(azm)
                            
                        if Data['LRazimuth'][1] < Data['LRazimuth'][0]:
                            Data['LRazimuth'][1] += 360
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
                PlotImage(self)
            self.itemPicked = None
            
    try:
        plotNum = self.G2plotNB.plotList.index('2D Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.figure.clf()
        Plot = Page.figure.gca()          #get a fresh plot after clf()
        
    except ValueError,error:
        Plot = self.G2plotNB.add('2D Powder Image').gca()
        plotNum = self.G2plotNB.plotList.index('2D Powder Image')
        Page = self.G2plotNB.nb.GetPage(plotNum)
        Page.canvas.mpl_connect('key_press_event', OnImPlotKeyPress)
        Page.canvas.mpl_connect('motion_notify_event', OnImMotion)
        Page.canvas.mpl_connect('pick_event', OnImPick)
        Page.canvas.mpl_connect('button_release_event', OnImRelease)
    if not event:                       #event from GUI TextCtrl - don't want focus to change to plot!!!
        Page.SetFocus()
    Plot.set_title(self.PatternTree.GetItemText(self.Image)[4:])
    size,imagefile = self.PatternTree.GetItemPyData(self.Image)
    if imagefile != self.oldImagefile:
        self.ImageZ = G2IO.GetImageData(imagefile,imageOnly=True)
        self.oldImagefile = imagefile
    Data = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
    Masks = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))

    imScale = 1
    if len(self.ImageZ) > 1024:
        imScale = len(self.ImageZ)/1024
    pixelSize = Data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    xmax = len(self.ImageZ)
    Xmax = len(self.ImageZ)*pixelSize[0]/1000.
    xlim = (-0.5,Xmax-.5)
    ylim = (Xmax-.5,-0.5,)
    Imin,Imax = Data['range'][1]
    acolor = mpl.cm.get_cmap(Data['color'])
    xcent,ycent = Data['center']
    Plot.set_xlabel('Image x-axis, mm',fontsize=12)
    Plot.set_ylabel('Image y-axis, mm',fontsize=12)
    #do threshold mask - "real" mask - others are just bondaries
    Zlim = Masks['Thresholds'][1]
    wx.BeginBusyCursor()
    try:
        MA = ma.masked_greater(ma.masked_less(self.ImageZ,Zlim[0]),Zlim[1])
        MaskA = ma.getmaskarray(MA)
        A = G2img.ImageCompress(MA,imScale)
        AM = G2img.ImageCompress(MaskA,imScale)
        
        ImgM = Plot.imshow(AM,aspect='equal',cmap='Reds',
            interpolation='nearest',vmin=0,vmax=2,extent=[0,Xmax,Xmax,0])
        Img = Plot.imshow(A,aspect='equal',cmap=acolor,
            interpolation='nearest',vmin=Imin,vmax=Imax,extent=[0,Xmax,Xmax,0])
        if self.setPoly:
            Img.set_picker(True)
    
        Plot.plot(xcent,ycent,'x')
        if Data['showLines']:
            LRAzim = Data['LRazimuth']                  #NB: integers
            IOtth = Data['IOtth']
            wave = Data['wavelength']
            dspI = wave/(2.0*sind(IOtth[0]/2.0))
            ellI = G2img.GetEllipse(dspI,Data)           #=False if dsp didn't yield an ellipse (ugh! a parabola or a hyperbola)
            dspO = wave/(2.0*sind(IOtth[1]/2.0))
            ellO = G2img.GetEllipse(dspO,Data)           #Ditto & more likely for outer ellipse
            if Data['fullIntegrate']:
                Azm = np.array(range(0,361))
            else:
                Azm = np.array(range(LRAzim[0],LRAzim[1]+1))
            if ellI:
                xyI = []
                for azm in Azm:
                    xyI.append(G2img.GetDetectorXY(dspI,azm,Data))
                xyI = np.array(xyI)
                arcxI,arcyI = xyI.T
                Plot.plot(arcxI,arcyI,picker=3)
            if ellO:
                xyO = []
                for azm in Azm:
                    xyO.append(G2img.GetDetectorXY(dspO,azm,Data))
                xyO = np.array(xyO)
                arcxO,arcyO = xyO.T
                Plot.plot(arcxO,arcyO,picker=3)
            if ellO and ellI and not Data['fullIntegrate']:
                Plot.plot([arcxI[0],arcxO[0]],[arcyI[0],arcyO[0]],picker=3)
                Plot.plot([arcxI[-1],arcxO[-1]],[arcyI[-1],arcyO[-1]],picker=3)
        for xring,yring in Data['ring']:
            Plot.plot(xring,yring,'r+',picker=3)
        if Data['setRings']:
            rings = np.concatenate((Data['rings']),axis=0)
            for xring,yring,dsp in rings:
                Plot.plot(xring,yring,'r+')            
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
            self.polyList.append([Plot.plot(x,y,'r',picker=3),ipoly])            
        colorBar = Page.figure.colorbar(Img)
        Plot.set_xlim(xlim)
        Plot.set_ylim(ylim)
        if not newPlot:
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
        
    except ValueError,error:
        Plot = self.G2plotNB.add('2D Integration').gca()
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
    if Data['setRings']:
        rings = np.concatenate((Data['rings']),axis=0)
        for xring,yring,dsp in rings:
            x,y = G2img.GetTthAzm(xring,yring,Data)
            Plot.plot(x,y,'r+')
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2img.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2img.GetTthAzm(x,y,Data)
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
    #a test plot routine - not normally used 
            
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
        
    except ValueError,error:
        Plot = self.G2plotNB.add('2D Transformed Powder Image').gca()
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
        
        