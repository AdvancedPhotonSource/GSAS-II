import math
import time
import copy
import numpy as np
import wx
import wx.aui
import matplotlib as mpl
import GSASIIgrid as G2gd
import GSASIIcomp as G2cmp
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
        
        self.plotList = []
            
    def add(self,name=""):
        page = G2Plot(self.nb)
        self.nb.AddPage(page,name)
        
        self.plotList.append(name)
        
        return page.figure
        
def PlotSngl(self):
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
            snglPage.canvas.SetToolTipString(HKLtxt)
            self.G2plotNB.status.SetFields(['HKL = '+HKLtxt,])
                
    def OnSCPick(event):
        zpos = Data['Layer']
        pos = event.artist.center
        if '100' in Data['Zone']:
            snglPage.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(zpos,pos[0],pos[1]))
            hkl = [zpos,pos[0],pos[1]]
        elif '010' in Data['Zone']:
            snglPage.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],zpos,pos[1]))
            hkl = [pos[0],zpos,pos[1]]
        elif '001' in Data['Zone']:
            snglPage.canvas.SetToolTipString('(picked:(%3d,%3d,%3d))'%(pos[0],pos[1],zpos))
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
        snglPage = self.G2plotNB.nb.GetPage(plotNum)
        snglPlot = snglPage.figure.gca()
        snglPlot.cla()
    except ValueError,error:
        snglPlot = self.G2plotNB.add('Structure Factors').gca()
        plotNum = self.G2plotNB.plotList.index('Structure Factors')
        snglPage = self.G2plotNB.nb.GetPage(plotNum)
        snglPage.canvas.mpl_connect('key_press_event', OnSCKeyPress)
        snglPage.canvas.mpl_connect('pick_event', OnSCPick)
        snglPage.canvas.mpl_connect('motion_notify_event', OnSCMotion)
    snglPage.SetFocus()
    
    snglPlot.set_aspect(aspect='equal')
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
    snglPlot.set_title(self.PatternTree.GetItemText(self.Sngl)[5:])
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
                snglPlot.add_artist(Circle(xy,radius=A,ec='g',fc='w',picker=3))
            if B:
                snglPlot.add_artist(Circle(xy,radius=B,ec='b',fc='w'))
                radius = C
                if radius > 0:
                    if A > B:
                        snglPlot.add_artist(Circle(xy,radius=radius,ec='r',fc='r'))
                    else:                    
                        snglPlot.add_artist(Circle(xy,radius=radius,ec='g',fc='g'))
    HKL = np.array(HKL)
    snglPlot.set_xlabel(xlabel[izone]+str(Data['Layer']),fontsize=12)
    snglPlot.set_ylabel(ylabel[izone],fontsize=12)
    snglPlot.set_xlim((HKLmin[pzone[izone][0]],HKLmax[pzone[izone][0]]))
    snglPlot.set_ylim((HKLmin[pzone[izone][1]],HKLmax[pzone[izone][1]]))
    snglPage.canvas.draw()
       
def PlotImage(self):
    from matplotlib.patches import Ellipse,Arc

    def OnImMotion(event):
        imgPage.canvas.SetToolTipString('')
        size = len(self.ImageZ)
        if (xlim[0] < event.xdata < xlim[1]) & (ylim[0] > event.ydata > ylim[1]):
            Data = self.PatternTree.GetItemPyData( \
                G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
            imgPage.canvas.SetCursor(wx.CROSS_CURSOR)
            item = self.itemPicked
            pixelSize = Data['pixelSize']
            scalex = 1000./pixelSize[0]
            scaley = 1000./pixelSize[1]                    
            if item and self.PatternTree.GetItemText(self.PickId) == 'Image Controls':
                if 'Text' in str(item):
                    imgPage.canvas.SetToolTipString('%8.3f %8.3fmm'%(event.xdata,event.ydata))
                else:
                    xcent,ycent = Data['center']
                    xpos = event.xdata-xcent
                    ypos = event.ydata-ycent
                    if 'line3' in  str(item) or 'line4' in str(item) and not Data['fullIntegrate']:
                        ang = int(atan2d(xpos,ypos))
                        imgPage.canvas.SetToolTipString('%6d deg'%(ang))
                    elif 'line1' in  str(item) or 'line2' in str(item):
                        tth = G2cmp.GetTth(event.xdata,event.ydata,Data)
                        imgPage.canvas.SetToolTipString('%8.3fdeg'%(tth))                           
            else:
                xpos = event.xdata
                ypos = event.ydata
                xpix = xpos*scalex
                ypix = ypos*scaley
                if (0 <= xpix <= size) and (0 <= ypix <= size):
                    imgPage.canvas.SetToolTipString('%6d'%(self.ImageZ[ypix][xpix]))
                tth,azm,dsp = G2cmp.GetTthDspAzm(xpos,ypos,Data)
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
            if Ypos and not imgPage.toolbar._active:         #make sure zoom/pan not selected
                if event.button == 1:
                    Xpix = Xpos*scalex
                    Ypix = Ypos*scaley
                    xpos,ypos,I,J = G2cmp.ImageLocalMax(self.ImageZ,20,Xpix,Ypix)
                    if I and J:
                        xpos /= scalex
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
                    tth,azm,dsp = G2cmp.GetTthDspAzm(xpos,ypos,Data)
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
        imgPage = self.G2plotNB.nb.GetPage(plotNum)
        imgPage.figure.clf()
        imgPlot = imgPage.figure.gca()
        if imgPage.views:
            imgPage.toolbar._views = copy.deepcopy(imgPage.views)
        view = imgPage.toolbar._views.forward()
        
    except ValueError,error:
        imgPlot = self.G2plotNB.add('2D Powder Image').gca()
        plotNum = self.G2plotNB.plotList.index('2D Powder Image')
        imgPage = self.G2plotNB.nb.GetPage(plotNum)
        imgPage.canvas.mpl_connect('key_press_event', OnImPlotKeyPress)
        imgPage.canvas.mpl_connect('motion_notify_event', OnImMotion)
        imgPage.canvas.mpl_connect('pick_event', OnImPick)
        imgPage.canvas.mpl_connect('button_release_event', OnImRelease)
        imgPage.views = False
        view = False
    imgPage.SetFocus()
        
    imgPlot.set_title(self.PatternTree.GetItemText(self.Image)[4:])
    size,self.ImageZ = self.PatternTree.GetItemPyData(self.Image)
    Data = self.PatternTree.GetItemPyData( \
        G2gd.GetPatternTreeItemId(self,self.Image, 'Image Controls'))
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
    if self.Img:
        xlim = self.Img.axes.get_xlim()
        ylim = self.Img.axes.get_ylim()
    Imin,Imax = Data['range'][1]
    acolor = mpl.cm.get_cmap(Data['color'])
    xcent,ycent = Data['center']
    imgPlot.set_xlabel('Image x-axis, mm',fontsize=12)
    imgPlot.set_ylabel('Image y-axis, mm',fontsize=12)
    A = G2cmp.ImageCompress(self.ImageZ,imScale)
    self.Img = imgPlot.imshow(A,aspect='equal',cmap=acolor, \
        interpolation='nearest',vmin=Imin,vmax=Imax,extent=[0,Xmax,Xmax,0])
    imgPlot.plot(xcent,ycent,'x')
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
            imgPlot.plot(arcxI,arcyI,picker=3)
        if ellO:
            xyO = []
            for azm in Azm:
                xyO.append(G2cmp.GetDetectorXY(dspO,azm,Data))
            xyO = np.array(xyO)
            arcxO,arcyO = xyO.T
            imgPlot.plot(arcxO,arcyO,picker=3)
        if ellO and ellI and not Data['fullIntegrate']:
            imgPlot.plot([arcxI[0],arcxO[0]],[arcyI[0],arcyO[0]],picker=3)
            imgPlot.plot([arcxI[-1],arcxO[-1]],[arcyI[-1],arcyO[-1]],picker=3)
    for xring,yring in Data['ring']:
        imgPlot.text(xring,yring,'+',color='b',ha='center',va='center',picker=3)
    if Data['setRings']:
        rings = np.concatenate((Data['rings']),axis=0)
        for xring,yring,dsp in rings:
            imgPlot.text(xring,yring,'+',ha='center',va='center')            
    for ellipse in Data['ellipses']:
        cent,phi,[width,height],col = ellipse
        imgPlot.add_artist(Ellipse([cent[0],cent[1]],2*width,2*height,phi,ec=col,fc='none'))
        imgPlot.text(cent[0],cent[1],'+',color=col,ha='center',va='center')
    colorBar = imgPage.figure.colorbar(self.Img)
    if view:
        self.Img.axes.set_xlim(view[0][:2])
        self.Img.axes.set_ylim(view[0][2:])
    else:
        self.Img.axes.set_xlim(xlim)
        self.Img.axes.set_ylim(ylim)
    imgPage.canvas.draw()
            
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
        pkwPage = self.G2plotNB.nb.GetPage(plotNum)
        pkwPage.figure.clf()
        pkwPlot = pkwPage.figure.gca()
    except ValueError,error:
        pkwPlot = self.G2plotNB.add('Peak Widths').gca()
        plotNum = self.G2plotNB.plotList.index('Peak Widths')
        pkwPage = self.G2plotNB.nb.GetPage(plotNum)
    pkwPage.SetFocus()
    
    pkwPage.canvas.SetToolTipString('')
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
    pkwPlot.set_title('Instrument and sample peak widths')
    pkwPlot.set_ylabel(r'$\Delta q/q, \Delta d/d$',fontsize=14)
    pkwPlot.set_xlabel(r'$q, \AA^{-1}$',fontsize=14)
    pkwPlot.plot(X,Y,color='r',label='Gaussian')
    pkwPlot.plot(X,Z,color='g',label='Lorentzian')
    pkwPlot.plot(X,W,color='b',label='G+L')
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
    pkwPlot.plot(X,Y,'+',color='r',label='G peak')
    pkwPlot.plot(X,Z,'+',color='g',label='L peak')
    pkwPlot.plot(X,W,'+',color='b',label='G+L peak')
    pkwPlot.legend(loc='best')
    pkwPage.canvas.draw()

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
        view = pdrPage.toolbar._views.forward()
        if view and 'line2' in str(pick):           #apply offset only for picked powder pattern points
            ind += np.searchsorted(xye[0],view[0][0])
        xy = zip(xpos[ind],ypos[ind])[0]
        if self.PatternTree.GetItemText(PickId) == 'Peak List':
            if ind.all() != [0]:                                    #picked a data point
                inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))
                if len(inst[1]) == 10:
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
        else:
            event.Skip(True)
        PlotPatterns(self)
                        
    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            ypos = event.ydata
            wave = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Instrument Parameters'))[0][1]
            dsp = 0.0
            if abs(xpos) > 0.:
                dsp = wave/(2.*sind(abs(xpos)/2.0))
            pdrPage.canvas.SetCursor(wx.CROSS_CURSOR)
            self.G2plotNB.status.SetFields(['2-theta =%9.3f d =%9.5f Intensity =%9.1f'%(xpos,dsp,ypos),])
            if self.itemPicked:
                pdrPage.canvas.SetToolTipString('%9.3f'%(xpos))
                   
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
        pdrPage = self.G2plotNB.nb.GetPage(plotNum)
        if not newPlot:
            pdrPlot = pdrPage.figure.gca()          #get previous powder plot & get limits
            xylim = pdrPlot.get_xlim(),pdrPlot.get_ylim()
        pdrPage.figure.clf()
        pdrPlot = pdrPage.figure.gca()          #get a fresh plot after clf()
    except ValueError,error:
        newPlot = True
        pdrPlot = self.G2plotNB.add('Powder Patterns').gca()
        plotNum = self.G2plotNB.plotList.index('Powder Patterns')
        pdrPage = self.G2plotNB.nb.GetPage(plotNum)
        pdrPage.canvas.mpl_connect('key_press_event', OnPlotKeyPress)
        pdrPage.canvas.mpl_connect('motion_notify_event', OnMotion)
        pdrPage.canvas.mpl_connect('pick_event', OnPick)
        pdrPage.canvas.mpl_connect('button_release_event', OnRelease)
        
    pdrPage.SetFocus()

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
    pdrPlot.set_title('Powder Patterns')
    pdrPlot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    pdrPlot.set_ylabel('Intensity',fontsize=12)
    if self.Contour:
        ContourZ = []
        ContourY = []
        Nseq = 0
    for Pattern in PlotList:
        ifpicked = False
        LimitId = 0
        xye = Pattern[1]
        if PickId:
            ifpicked = Pattern[2] == self.PatternTree.GetItemText(PatternId)
            LimitId = G2gd.GetPatternTreeItemId(self,PatternId, 'Limits')
        N = PlotList.index(Pattern)
        X = xye[0]
        Y = xye[1]+offset*N
        if LimitId:
            limits = self.PatternTree.GetItemPyData(LimitId)
            Lines.append(pdrPlot.axvline(limits[1][0],color='g',dashes=(5,5),picker=3.))    
            Lines.append(pdrPlot.axvline(limits[1][1],color='r',dashes=(5,5),picker=3.))                    
        if self.Contour:
            ContourY.append(N)
            ContourZ.append(Y)
            ContourX = X
            Nseq += 1
            pdrPlot.set_ylabel('Data sequence',fontsize=12)
        else:
            if ifpicked:
                Z = xye[3]+offset*N
                W = xye[4]+offset*N
                D = xye[5]+offset*N
                if self.Weight:
                    W2 = np.sqrt(xye[2])
                    D *= W2
                pdrPlot.plot(X,Y,colors[N%6]+'+',picker=3.,clip_on=False)
                pdrPlot.plot(X,Z,colors[(N+1)%6],picker=False)
                pdrPlot.plot(X,W,colors[(N+2)%6],picker=False)
                pdrPlot.plot(X,D,colors[(N+3)%6],picker=False)
                pdrPlot.axhline(0.,color=wx.BLACK)
                pdrPage.canvas.SetToolTipString('')
                if self.PatternTree.GetItemText(PickId) == 'Peak List':
                    tip = 'On data point: Pick peak - L or R MB.On line: MB down to move'
                    pdrPage.canvas.SetToolTipString(tip)
                    data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Peak List'))
                    for item in data:
                        Lines.append(pdrPlot.axvline(item[0],color=colors[N%6],picker=2.))
                if self.PatternTree.GetItemText(PickId) == 'Limits':
                    tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                    pdrPage.canvas.SetToolTipString(tip)
                    data = self.LimitsTable.GetData()
            else:
                pdrPlot.plot(X,Y,colors[N%6],picker=False)
    if PickId and self.PatternTree.GetItemText(PickId) in ['Index Peak List','Unit Cells List']:
        peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
        for peak in peaks:
            pdrPlot.axvline(peak[0],color='b')
        for hkl in self.HKL:
            pdrPlot.axvline(hkl[5],color='r',dashes=(5,5))
    if self.Contour:
        acolor = mpl.cm.get_cmap('Paired')
        pdrPlot.contourf(ContourX,ContourY,ContourZ,cmap=acolor)
#        pdrPlot.set_ylim(0,Nseq-1)
    else:
        self.Lines = Lines
    if not newPlot:
        pdrPage.toolbar.push_current()
        pdrPlot.set_xlim(xylim[0])
        pdrPlot.set_ylim(xylim[1])
        xylim = []
        pdrPage.toolbar.push_current()
        pdrPage.toolbar.draw()
    else:
        pdrPage.canvas.draw()

    
    self.Pwdr = True

def PlotPowderLines(self):

    def OnMotion(event):
        xpos = event.xdata
        if xpos:                                        #avoid out of frame mouse position
            pksPage.canvas.SetCursor(wx.CROSS_CURSOR)
            self.G2plotNB.status.SetFields(['2-theta =%9.3f '%(xpos,),])

    try:
        plotNum = self.G2plotNB.plotList.index('Powder Lines')
        pksPage = self.G2plotNB.nb.GetPage(plotNum)
        pksPage.figure.clf()
        pksPlot = pksPage.figure.gca()
    except ValueError,error:
        newPlot = True
        pksPlot = self.G2plotNB.add('Powder Lines').gca()
        plotNum = self.G2plotNB.plotList.index('Powder Lines')
        pksPage = self.G2plotNB.nb.GetPage(plotNum)
        pksPage.canvas.mpl_connect('motion_notify_event', OnMotion)
        
    pksPage.SetFocus()
    pksPlot.set_title('Powder Pattern Lines')
    pksPlot.set_xlabel(r'$\mathsf{2\theta}$',fontsize=14)
    PickId = self.PickId
    PatternId = self.PatternId
    peaks = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PatternId, 'Index Peak List'))
    for peak in peaks:
        pksPlot.axvline(peak[0],color='b')
    for hkl in self.HKL:
        pksPlot.axvline(hkl[5],color='r',dashes=(5,5))
    xmin = peaks[0][0]
    xmax = peaks[-1][0]
    delt = xmax-xmin
    xlim = [max(0,xmin-delt/20.),min(180.,xmax+delt/20.)]
    pksPlot.set_xlim(xlim)
    pksPage.canvas.draw()


def PlotTRImage(self):
            
    def OnMotion(event):
        trimgPage.canvas.SetToolTipString('')
        trimgPage.canvas.SetCursor(wx.CROSS_CURSOR)
        azm = event.xdata
        tth = event.ydata
        if azm and tth:
            self.G2plotNB.status.SetFields(\
                ['Detector 2-th =%9.2fdeg, azm = %7.2fdeg'%(tth,azm),])
                    
    def OnPick(event):
        if self.PatternTree.GetItemText(self.PickId) != 'Image Controls':
            return
        if self.itemPicked is not None: return
        pick = event.artist
        self.itemPicked = pick
        
    def OnRelease(event):
        if self.PatternTree.GetItemText(self.PickId) != 'Image Controls':
            return
        Data = self.PatternTree.GetItemPyData(self.PickId)
        if self.itemPicked:
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                ypos = event.ydata
                if 'Line2D' in str(self.itemPicked):
                    if 'line0' in str(self.itemPicked):
                        Data['IOtth'][0] = ypos
                    elif 'line1' in str(self.itemPicked):
                        Data['IOtth'][1] = ypos
                    elif 'line2' in str(self.itemPicked) and not Data['fullIntegrate']:
                        Data['LRazimuth'][0] = int(xpos)
                    elif 'line3' in str(self.itemPicked) and not Data['fullIntegrate']:
                        Data['LRazimuth'][1] = int(xpos)
                        
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
                PlotTRImage(self)
            self.itemPicked = None
            
    try:
        plotNum = self.G2plotNB.plotList.index('2D Transformed Powder Image')
        trimgPage = self.G2plotNB.nb.GetPage(plotNum)
        trimgPage.figure.clf()
        trimgPlot = trimgPage.figure.gca()
        if trimgPage.views:
            trimgPage.toolbar._views = copy.deepcopy(trimgPage.views)
        view = trimgPage.toolbar._views.forward()
        
    except ValueError,error:
        trimgPlot = self.G2plotNB.add('2D Transformed Powder Image').gca()
        plotNum = self.G2plotNB.plotList.index('2D Transformed Powder Image')
        trimgPage = self.G2plotNB.nb.GetPage(plotNum)
        trimgPage.canvas.mpl_connect('motion_notify_event', OnMotion)
        trimgPage.canvas.mpl_connect('pick_event', OnPick)
        trimgPage.canvas.mpl_connect('button_release_event', OnRelease)
        trimgPage.views = False
        view = False
    trimgPage.SetFocus()
        
    data = self.PatternTree.GetItemPyData(self.PickId)
    image = self.ImageZ
    Iz = len(image)
    Imin,Imax = data['range'][1]
    step = (Imax-Imin)/5.
    V = np.arange(Imin,Imax,step)
    acolor = mpl.cm.get_cmap('Paired')
    trimgPlot.set_xlabel('azimuth',fontsize=12)
    trimgPlot.set_ylabel('2-theta',fontsize=12)
    trimgPlot.contour(self.TA[1],self.TA[0],image,V,cmap=acolor)
    if data['showLines']:
        IOtth = data['IOtth']
        LRAzim = data['LRazimuth']                  #NB: integers
        trimgPlot.plot([LRAzim[0],LRAzim[1]],[IOtth[0],IOtth[0]],picker=True)
        trimgPlot.plot([LRAzim[0],LRAzim[1]],[IOtth[1],IOtth[1]],picker=True)
        trimgPlot.plot([LRAzim[0],LRAzim[0]],[IOtth[0],IOtth[1]],picker=True)
        trimgPlot.plot([LRAzim[1],LRAzim[1]],[IOtth[0],IOtth[1]],picker=True)
    trimgPage.canvas.draw()
            
  