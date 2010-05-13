import math
import time
import copy
import numpy as np
import wx
import wx.aui
import matplotlib as mpl
import GSASIIgrid as G2gd
import GSASIIcomp as G2cmp
import GSASIIIO as G2IO
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
#        self.Bind(wx.aui.EVT_AUI_PAGE_CHANGED, self.OnPageChanged)
        
        self.plotList = []
            
    def add(self,name=""):
        page = G2Plot(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
    def OnPageChanged(self,event):
        print 'page changed'
        
def PlotSngl(self,newPlot=False):
    from matplotlib.patches import Circle

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
            self.G2plotNB.status.SetFields(['HKL = '+HKLtxt,])
                
    def OnSCPick(event):
        zpos = Data['Layer']
        pos = event.artist.center
        if '100' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(zpos,pos[0],pos[1]))
            hkl = [zpos,pos[0],pos[1]]
        elif '010' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],zpos,pos[1]))
            hkl = [pos[0],zpos,pos[1]]
        elif '001' in Data['Zone']:
            Page.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],pos[1],zpos))
            hkl = [pos[0],pos[1],zpos]
        h,k,l = hkl
        i = HKL.all(hkl)
        print i
        HKLtxt = '(%3d,%3d,%3d %10.2f %6.3f %10.2f)'%(h,k,l,Fosq,sig,Fcsq)
        self.G2plotNB.status.SetFields(['HKL, Fosq, sig, Fcsq = '+HKLtxt,])
                         
        
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
        Page.canvas.mpl_connect('key_press_event', OnSCKeyPress)
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
    for H,Fosq,sig,Fcsq,x,x,x in HKLref:
        HKL.append(H)
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
                        Plot.add_artist(Circle(xy,radius=radius,ec='r',fc='r'))
                    else:                    
                        Plot.add_artist(Circle(xy,radius=radius,ec='g',fc='g'))
    HKL = np.array(HKL)
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
                XY = [xy[0],0, xy[1],1, sig,0, gam,0, ins[5],0]       #default refine intensity 1st   
                data.append(XY)
                G2gd.UpdatePeakGrid(self,data)
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
                G2gd.UpdateLimitsGrid(self,data)
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
            print 'contouring'
            if self.Contour:
                self.Contour = False
            else:
                self.Contour = True
        PlotPatterns(self)
                        
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
                self.G2plotNB.status.SetFields(['2-theta =%9.3f d =%9.5f Intensity =%9.1f'%(xpos,dsp,ypos),])
                if self.itemPicked:
                    Page.canvas.SetToolTipString('%9.3f'%(xpos))
            except TypeError:
                self.G2plotNB.status.SetFields(['Select PWDR powder pattern first',])
                
                   
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
                print 'limits',xpos
                data[1][lineNo] = xpos
                self.PatternTree.SetItemPyData(LimitId,data)
                if self.PatternTree.GetItemText(self.PickId) == 'Limits':
                    G2gd.UpdateLimitsGrid(self,data)
            else:
                PeakId = G2gd.GetPatternTreeItemId(self,self.PatternId, 'Peak List')
                data = self.PatternTree.GetItemPyData(PeakId)
                print 'peaks',xpos
                data[lineNo-2][0] = xpos
                self.PatternTree.SetItemPyData(PeakId,data)
                G2gd.UpdatePeakGrid(self,data)
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

    PickId = self.PickId
    PatternId = self.PatternId
    colors=['b','g','r','c','m','k']
    PlotList = []
    Lines = []
    item, cookie = self.PatternTree.GetFirstChild(self.root)
    while item:
        if 'PWDR' in self.PatternTree.GetItemText(item):
            Pattern = self.PatternTree.GetItemPyData(item)
            Pattern.append(self.PatternTree.GetItemText(item))
            PlotList.append(Pattern)
        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
    Ymax = 1.0
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
                D = xye[5]+offset*N
                if self.Weight:
                    W2 = np.sqrt(xye[2])
                    D *= W2
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
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
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

    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            Page.canvas.SetCursor(wx.CROSS_CURSOR)
            self.G2plotNB.status.SetFields(['2-theta =%9.3f '%(xpos,),])

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
    for peak in peaks:
        X.append(4.0*math.pi*sind(peak[0]/2.0)/lam)
        s = 1.17741*math.sqrt(peak[4])*math.pi/18000.
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
            
def PlotExposedImage(self,newPlot=False):
    plotNo = self.G2plotNB.nb.GetSelection()
    if self.G2plotNB.nb.GetPageText(plotNo) == '2D Powder Image':
        PlotImage(self,newPlot)
    elif self.G2plotNB.nb.GetPageText(plotNo) == '2D Integration':
        PlotIntegration(self,newPlot)

def PlotImage(self,newPlot=False):
    from matplotlib.patches import Ellipse,Arc

    def OnImMotion(event):
        Page.canvas.SetToolTipString('')
        size = len(self.ImageZ)
        if event.xdata and event.ydata:                 #avoid out of frame errors
            Data = self.PatternTree.GetItemPyData( \
                G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
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
                        tth = G2cmp.GetTth(event.xdata,event.ydata,Data)
                        Page.canvas.SetToolTipString('%8.3fdeg'%(tth))                           
            else:
                xpos = event.xdata
                ypos = event.ydata
                xpix = xpos*scalex
                ypix = ypos*scaley
                if (0 <= xpix <= size) and (0 <= ypix <= size):
                    Page.canvas.SetToolTipString('%6d'%(self.ImageZ[ypix][xpix]))
                tth,azm,dsp = G2cmp.GetTthAzmDsp(xpos,ypos,Data)
                Q = 2.*math.pi/dsp
                self.G2plotNB.status.SetFields(\
                    ['Detector 2-th =%9.2fdeg, dsp =%9.3fA, Q = %6.3fA-1, azm = %7.2fdeg'%(tth,dsp,Q,azm),])

    def OnImPlotKeyPress(event):
        if self.PatternTree.GetItemText(self.PickId) == 'Image Controls':
            Data = self.PatternTree.GetItemPyData(self.PickId)
            pixelSize = Data['pixelSize']
            size = len(self.ImageZ)
            Xpos = event.xdata
            if not Xpos:            #got point out of frame
                return
            Ypos = event.ydata
            if event.key == 'm':
                print 'mask = ',Xpos,Ypos
            
    def OnImPick(event):
        if self.PatternTree.GetItemText(self.PickId) != 'Image Controls':
            return
        if self.itemPicked is not None: return
        pick = event.artist
        self.itemPicked = pick
        
    def OnImRelease(event):
        if self.PatternTree.GetItemText(self.PickId) != 'Image Controls':
            return
        Data = self.PatternTree.GetItemPyData(self.PickId)
        pixelSize = Data['pixelSize']
        scalex = 1000./pixelSize[0]
        scaley = 1000./pixelSize[1]
        if self.itemPicked is None:
            size = len(self.ImageZ)
            Xpos = event.xdata
            if not (Xpos and self.ifGetRing):                   #got point out of frame
                return
            Ypos = event.ydata
            if Ypos and not Page.toolbar._active:         #make sure zoom/pan not selected
                if event.button == 1:
                    Xpix = Xpos*scalex
                    Ypix = Ypos*scaley
                    xpos,ypos,I,J = G2cmp.ImageLocalMax(self.ImageZ,20,Xpix,Ypix)
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
                    tth,azm,dsp = G2cmp.GetTthAzmDsp(xpos,ypos,Data)
                    if 'Line2D' in str(self.itemPicked):
                        if 'line1' in str(self.itemPicked):
                            Data['IOtth'][0] = tth
                        elif 'line2' in str(self.itemPicked):
                            Data['IOtth'][1] = tth
                        elif 'line3' in str(self.itemPicked) and not Data['fullIntegrate']:
                            Data['LRazimuth'][0] = int(azm)
                        elif 'line4' in str(self.itemPicked) and not Data['fullIntegrate']:
                            Data['LRazimuth'][1] = int(azm)
                            
                        if Data['LRazimuth'][1] < Data['LRazimuth'][0]:
                            Data['LRazimuth'][1] += 360
                        if  Data['IOtth'][0] > Data['IOtth'][1]:
                            Data['IOtth'] = G2cmp.SwapXY(Data['IOtth'][0],Data['IOtth'][1])
                            
                        self.InnerTth.SetValue("%8.2f" % (Data['IOtth'][0]))
                        self.OuterTth.SetValue("%8.2f" % (Data['IOtth'][1]))
                        self.Lazim.SetValue("%6d" % (Data['LRazimuth'][0]))
                        self.Razim.SetValue("%6d" % (Data['LRazimuth'][1]))
                    else:
                        print event.xdata,event.ydata,event.button
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
    Page.SetFocus()
        
    Plot.set_title(self.PatternTree.GetItemText(self.Image)[4:])
    size,imagefile = self.PatternTree.GetItemPyData(self.Image)
    if imagefile != self.oldImagefile:
        self.ImageZ = G2IO.GetImageData(imagefile,imageOnly=True)
        self.oldImagefile = imagefile
    Data = self.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
    try:
        Masks = self.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(self,self.Image, 'Masks'))
    except TypeError:       #missing Masks
        Masks = {}
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
    #need "applyMask" routine here
    A = G2cmp.ImageCompress(self.ImageZ,imScale)
    Img = Plot.imshow(A,aspect='equal',cmap=acolor,
        interpolation='nearest',vmin=Imin,vmax=Imax,extent=[0,Xmax,Xmax,0])

    Plot.plot(xcent,ycent,'x')
    if Data['showLines']:
        LRAzim = Data['LRazimuth']                  #NB: integers
        IOtth = Data['IOtth']
        wave = Data['wavelength']
        dspI = wave/(2.0*sind(IOtth[0]/2.0))
        ellI = G2cmp.GetEllipse(dspI,Data)           #=False if dsp didn't yield an ellipse (ugh! a parabola or a hyperbola)
        dspO = wave/(2.0*sind(IOtth[1]/2.0))
        ellO = G2cmp.GetEllipse(dspO,Data)           #Ditto & more likely for outer ellipse
        if Data['fullIntegrate']:
            Azm = np.array(range(0,361))
        else:
            Azm = np.array(range(LRAzim[0],LRAzim[1]+1))
        if ellI:
            xyI = []
            for azm in Azm:
                xyI.append(G2cmp.GetDetectorXY(dspI,azm,Data))
            xyI = np.array(xyI)
            arcxI,arcyI = xyI.T
            Plot.plot(arcxI,arcyI,picker=3)
        if ellO:
            xyO = []
            for azm in Azm:
                xyO.append(G2cmp.GetDetectorXY(dspO,azm,Data))
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
        
def PlotIntegration(self,newPlot=False):
            
    def OnMotion(event):
        Page.canvas.SetToolTipString('')
        Page.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.ydata
        tth = event.xdata
        if azm and tth:
            self.G2plotNB.status.SetFields(\
                ['Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm),])
                                
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
            x,y = G2cmp.GetTthAzm(xring,yring,Data)
            Plot.plot(x,y,'r+')
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2cmp.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2cmp.GetTthAzm(x,y,Data)
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
                ['Detector 2-th =%9.3fdeg, azm = %7.2fdeg'%(tth,azm),])
                                
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
            x,y = G2cmp.GetTthAzm(xring,yring,Data)
            Plot.plot(y,x,'r+')            
    if Data['ellipses']:            
        for ellipse in Data['ellipses']:
            ring = np.array(G2cmp.makeIdealRing(ellipse[:3])) #skip color
            x,y = np.hsplit(ring,2)
            tth,azm = G2cmp.GetTthAzm(x,y,Data)
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
        
        