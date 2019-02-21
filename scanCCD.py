#scanCCD data processing
'''
*scanCCD: reduce data from scanning CCD*
========================================

Quickly prototyped routine for reduction of data from detector described in
B.H. Toby, T.J. Madden, M.R. Suchomel, J.D. Baldwin, and R.B. Von Dreele,
"A Scanning CCD Detector for Powder Diffraction Measurements".
Journal of Applied Crystallography. 46(4): p. 1058-63 (2013).

'''
import os
import os.path as ospath
import sys
import math
import time
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import wx
import matplotlib as mpl
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIimage as G2img
import GSASIIplot as G2plt

npsind = lambda x: np.sin(x*np.pi/180.)
npcosd = lambda x: np.cos(x*np.pi/180.)
npacosd = lambda x: 180.*np.arccos(x)/np.pi
npasind = lambda x: 180.*np.arcsin(x)/np.pi

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

def create(parent):
    return scanCCD(parent)
    
[wxID_FILEEXIT, wxID_FILEOPEN, wxID_INTEGRATE, wxID_OUTPUT,
] = [wx.NewId() for _init_coll_File_Items in range(4)]

def FileDlgFixExt(dlg,file):            #this is needed to fix a problem in linux wx.FileDialog
    ext = dlg.GetWildcard().split('|')[2*dlg.GetFilterIndex()+1].strip('*')
    if ext not in file:
        file += ext
    return file
    
class scanCCD(wx.Frame):

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='scanCCD', parent=parent,
            size=wx.Size(460, 250),style=wx.DEFAULT_FRAME_STYLE, title='scanCCD')
        self.scanCCDMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(help='Open scanCCD image files (*.tif)', id=wxID_FILEOPEN,
             kind=wx.ITEM_NORMAL,text='Open scanCCD files')
        self.File.Append(help='Integrate scanCCD images',id=wxID_INTEGRATE,
             kind=wx.ITEM_NORMAL,text='Integrate scanCCD images')
        self.File.Append(help='Output fxye file from integration',id=wxID_OUTPUT,
             kind=wx.ITEM_NORMAL,text='Output pattern')
        self.File.Append(help='Exit from scanCCD', id=wxID_FILEEXIT, kind=wx.ITEM_NORMAL,
            text='Exit')
        self.Bind(wx.EVT_MENU, self.OnImageRead, id=wxID_FILEOPEN)
        self.Bind(wx.EVT_MENU,self.OnImageIntegrate,id=wxID_INTEGRATE)
        self.Bind(wx.EVT_MENU,self.OnOutput,id=wxID_OUTPUT)
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.scanCCDMenu.Append(menu=self.File, title='File')
        self.SetMenuBar(self.scanCCDMenu)
        self.SCCDPanel = wx.Panel(self)        
        plotFrame = wx.Frame(None,-1,'scanCCD Plots',size=wx.Size(700,600), \
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.plotNB = G2plt.G2PlotNoteBook(plotFrame)
        plotFrame.Show()

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dirname = ''
        self.imagefiles = []
        self.dataFrame = None
        self.itemPicked = None
        self.Image = []
        self.Hxyw = []
        self.data = {'color':'Paired','range':[[0,1000],[0,1000]]}

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def OnImageRead(self,event):
        dlg = wx.FileDialog(self, 'Choose scanCCD image files', '.', '',\
        'Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|\
        All files (*.*)|*.*',
        wx.FD_OPEN | wx.FD_MULTIPLE)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            self.imagefiles = []
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetDirectory()
                self.imagefiles = dlg.GetPaths()
                self.imagefiles.sort()
                self.data['imScale'] = 8
                self.Image = []
                for imagefile in self.imagefiles:
                    Comments,Data,Npix,image = G2IO.GetImageData(self,imagefile)
                    if Comments:
                        A = G2img.ImageCompress(image,self.data['imScale'])
                        if len(self.Image):
                            self.Image = np.concatenate((self.Image,A))
                        else:
                            self.Image = A
                self.data['range'][0] = self.data['range'][1] = np.max(self.Image)
                self.data['pixel'] = Data['pixelSize']
                self.data['size'] = self.Image.shape
                self.data['Zmin'] = np.min(self.Image)
                self.data['Zmax'] = np.max(self.Image)
                self.data['Zeros'] = [-15.3096,16.4]        #for bob1-3 files
                self.data['mm/deg'] = 15.56012              #from bob1 fit to 2333 peak positions
                self.data['radius'] = 180.0*self.data['mm/deg']/math.pi
                self.data['TBlimits'] = [0,len(image)]
                self.data['LRlimits'] = [0,len(image)*len(self.imagefiles)]
                self.data['PtScan'] = len(image)/2
                self.data['2thScan'] = [0.0,90.0,0.001]
                self.data['showBlk'] = False
                self.data['skip'] = 1               #default - skip ramp block
                if len(self.imagefiles) == 1:
                    self.data['skip'] = 0
                self.UpdateControls(event)
                self.PlotImage()
        finally:
            dlg.Destroy()
            
    def OnImageIntegrate(self,event):
        
        def Make2ThetaMap(data,iLim,jLim):
            #transforms scanCCD image from x,y space to 2-theta,y space
            pixelSize = data['pixel']
            scalex = pixelSize[0]/1000.
            scaley = pixelSize[1]/1000.
            vecA = np.asfarray([1.,0.,0.],dtype=np.float32)
            
            tay,tax = np.mgrid[jLim[0]+.5:jLim[1]+.5,iLim[0]+.5:iLim[1]+.5]        #bin centers not corners
            tax = np.asfarray((tax*scalex-data['Zeros'][0])/data['mm/deg'],dtype=np.float32)  #scanCCD 2-thetas
            tay = np.asfarray(tay*scaley-data['Zeros'][1],dtype=np.float32)
            vecB = np.array([npcosd(tax)*data['radius'],npsind(tax)*data['radius'],tay])
            norm = np.sqrt(np.sum((vecB.T*vecB.T),axis=2))
            vecB /= norm.T
            tax = npacosd(np.dot(vecB.T,vecA))*tax/np.abs(tax)      #to get sign of 2-theta
            tay += data['Zeros'][1]
            return tax,tay.T           #2-theta arrays & turn y array around!!

    
        def Fill2ThetaMap(data,TA,image):
            import numpy.ma as ma
            Zmin = data['Zmin']
            Zmax = data['Zmax']
            tax,tay = TA    # 2-theta & yaxis
            taz = ma.masked_outside(image.flatten()-Zmin,0,Zmax-Zmin)
            tam = ma.getmask(taz)
            tax = ma.compressed(ma.array(tax.flatten(),mask=tam))
            tay = ma.compressed(ma.array(tay.flatten(),mask=tam))
            taz = ma.compressed(ma.array(taz.flatten(),mask=tam))
            del(tam)
            return tax,tay,taz
        
        import histosigma2d as h2d
        print('Begin image integration')
        scaley = self.data['pixel'][1]/1000.
        tthStart,tthEnd,tthStep = self.data['2thScan']
        LUtth = [tthStart,tthEnd]
        TBlim = np.array(self.data['TBlimits'],dtype=np.float32)
        nYpix = TBlim[1]-TBlim[0]
        TBlim *= scaley
        TBdelt = TBlim[1]-TBlim[0]
        nTB = 1
        numChans = (tthEnd-tthStart)/tthStep
        NST = np.zeros(shape=(numChans,1),order='F',dtype=np.float32)
        H0 = np.zeros(shape=(numChans,1),order='F',dtype=np.float32)
        Amat = np.zeros(shape=(numChans,1),order='F',dtype=np.float32)
        Qmat = np.zeros(shape=(numChans,1),order='F',dtype=np.float32)
        imSize = np.array(self.data['size'])*self.data['imScale']
        H1 = [tth for tth in np.linspace(LUtth[0],LUtth[1],numChans)]
        blkSize = 2048
        N = 4096/blkSize
        nBlk = N**2*len(self.imagefiles)       #assume 4Kx4K CCD images - done in 1Kx1K blocks
        t0 = time.time()
        dlg = wx.ProgressDialog("Elapsed time","2D image integration",nBlk,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        try:
            nBlk = 0
            iBegi = 0
            iFin = 0
            GoOn = True
            for ifile,imagefile in enumerate(self.imagefiles):
                if ifile >= self.data['skip']:
                    image = G2IO.GetImageData(self,imagefile)[3]
                    image = np.fliplr(image)
                    for iBlk in range(N):
                        iBeg = iBegi+iBlk*blkSize
                        iFin = iBeg+blkSize
                        for jBlk in range(N):
                            jBeg = jBlk*blkSize
                            jFin = jBeg+blkSize                
                            print('Process map block:',nBlk,iBlk,jBlk,' offset:',iBegi,' limits:',iBeg,iFin,jBeg,jFin)
                            TA = Make2ThetaMap(self.data,(iBeg,iFin),(jBeg,jFin))           #2-theta & Y arrays & create position mask                       
                            Block = image[iBeg-iBegi:iFin-iBegi,jBeg:jFin]
                            tax,tay,taz = Fill2ThetaMap(self.data,TA,Block)                 #and apply masks
                            NST,H0,Amat,Qmat = h2d.histosigma2d(len(tax),tax,tay,taz,numChans,nTB,LUtth,TBlim,tthStep,TBdelt,NST,H0,Amat,Qmat)
                            H0temp = np.nan_to_num(np.divide(H0,NST))
                            Qtemp = np.nan_to_num(np.divide(Qmat,NST))
                            if self.data['showBlk']:
                                self.PlotBlock(Block.T,'%s %d %d'%('Block',iBlk,jBlk))
                                OKdlg = wx.MessageDialog(self,'','Continue',wx.OK)
                                try:
                                    result = OKdlg.ShowModal()
                                finally:
                                    OKdlg.Destroy()
                            del tax,tay,taz
                            nBlk += 1
                            GoOn = dlg.Update(nBlk)[0]
                            if not GoOn:
                                break
                else:
                    print('file '+imagefile+' skipped')
                    nBlk += N*N
                    GoOn = dlg.Update(nBlk)[0]
                    iFin += N*blkSize
                if not GoOn:
                    break
                iBegi = iFin
            H0 = np.divide(H0,NST)
            H0 = np.nan_to_num(H0)
            Qmat = np.divide(Qmat,NST)
            Qmat = np.nan_to_num(Qmat)
            del NST
            t1 = time.time()
        finally:
            dlg.Destroy()
        Scale = np.sum(Qmat)/np.sum(H0)
        print('SumI: ',np.sum(H0),' SumV: ',np.sum(Qmat),' Scale:',Scale)
        print('Integration complete')
        print("Elapsed time:","%8.3f"%(t1-t0), "s")
        self.Hxyw = [H1,Scale*H0.T[0],np.sqrt(Qmat.T[0])]
        print()
        self.PlotXY(self.Hxyw,True,type='Integration result')
        
    def OnOutput(self,event):

        def powderSave(self,powderfile,Fxye=False):
            file = open(powderfile,'w')
            file.write('#%s\n'%('from scanCCD image '+self.imagefiles[0]))
            print('save powder pattern to file: ',powderfile)
            wx.BeginBusyCursor()
            try:
                x,y,e = self.Hxyw
                if Fxye:
                    file.write(powderfile+'\n')
                    file.write('BANK 1 %d %d CONS %.2f %.2f 0 0 FXYE\n'%(len(x),len(x),\
                        100.*x[0],100.*(x[1]-x[0])))                    
                XYE = zip(x,y,e)
                for X,Y,E in XYE:
                    if Fxye:
                        file.write("%15.6g %15.6g %15.6g\n" % (100.*X,Y,max(E,1.0)))                        
                    else:
                        file.write("%15.6g %15.6g %15.6g\n" % (X,Y,max(E,1.0)))
                file.close()
            finally:
                wx.EndBusyCursor()
            print('powder pattern file written')
        
        if not self.Hxyw:
            return
        dlg = wx.FileDialog(self, 'Choose output powder file name', '.', '', 
            'GSAS fxye file (*.fxye)|*.fxye|xye file (*.xye)|*.xye',
            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                print(dlg.GetFilename())
                powderfile = dlg.GetPath()
                if 'fxye' in powderfile:
                    powderSave(self,powderfile,Fxye=True)
                else:       #just xye
                    powderSave(self,powderfile)
                self.dirname = dlg.GetDirectory()
        finally:
            dlg.Destroy()

    def UpdateControls(self,event):
        ObjIndx = {}
        self.SCCDPanel.DestroyChildren()
        
        def ColorSizer(data):
            
            def OnNewColorBar(event):
                colSel = event.GetEventObject()
                data['color'] = colSel.GetValue()
                self.PlotImage()
                
            def OnShowBlk(event):
                data['showBlk'] = shoBlk.GetValue()
                
            def OnSkipFiles(event):
                data['skip'] = int(skipFile.GetValue())
                            
            colorList = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
            colorSizer = wx.FlexGridSizer(0,5,5,5)
            colorSizer.Add(wx.StaticText(self.SCCDPanel,label=' Color bar '),0,wx.ALIGN_CENTER_VERTICAL)
            colSel = wx.ComboBox(self.SCCDPanel,value=data['color'],choices=colorList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN|wx.CB_SORT)
            colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
            colorSizer.Add(colSel,0,wx.ALIGN_CENTER_VERTICAL)
            colorSizer.Add(wx.StaticText(self.SCCDPanel,label='image files to skip:'),0,wx.ALIGN_CENTER_VERTICAL)
            skipFile = wx.ComboBox(self.SCCDPanel,value=str(data['skip']),
                choices=[str(i) for i in range(len(self.imagefiles))],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            skipFile.Bind(wx.EVT_COMBOBOX,OnSkipFiles)
            colorSizer.Add(skipFile,0,wx.ALIGN_CENTER_VERTICAL)
            shoBlk = wx.CheckBox(self.SCCDPanel,-1,label='Show image blocks?')
            shoBlk.SetValue(data['showBlk'])
            shoBlk.Bind(wx.EVT_CHECKBOX, OnShowBlk)
            colorSizer.Add(shoBlk,0,wx.ALIGN_CENTER_VERTICAL)
            return colorSizer                   
            
        def MaxSizer(data):
            
            def OnMaxSlider(event):            
                imax = int(self.maxSel.GetValue())*data['range'][0]/100.
                data['range'][1] = imax
                self.maxVal.SetValue('%d'%(data['range'][1]))
                self.PlotImage()
                
            def OnMaxValue(event):
                try:
                    value = int(self.maxVal.GetValue())
                    if value < 0 or value > data['range'][0]:
                        raise ValueError
                except ValueError:
                    value = data['range'][1]
                data['range'][1] = value
                self.maxSel.SetValue(int(100*value/data['range'][0]))
                self.PlotImage()
            
            maxSizer = wx.FlexGridSizer(0,3,0,5)
            maxSizer.AddGrowableCol(1,1)
            maxSizer.SetFlexibleDirection(wx.HORIZONTAL)
            maxSizer.Add(wx.StaticText(parent=self.SCCDPanel,label=' Max intensity'),0,
                wx.ALIGN_CENTER_VERTICAL)
            self.maxSel = wx.Slider(parent=self.SCCDPanel,style=wx.SL_HORIZONTAL,
                value=int(100*data['range'][1]/data['range'][0]))
            maxSizer.Add(self.maxSel,1,wx.EXPAND)
            self.maxSel.Bind(wx.EVT_SLIDER, OnMaxSlider)
            self.maxVal = wx.TextCtrl(parent=self.SCCDPanel,value='%d'%(data['range'][1]))
            self.maxVal.Bind(wx.EVT_TEXT_ENTER,OnMaxValue)
            self.maxVal.Bind(wx.EVT_KILL_FOCUS,OnMaxValue)
            maxSizer.Add(self.maxVal,0,wx.ALIGN_CENTER_VERTICAL)    
            return maxSizer    
                    
        def ZSizer(data):
            
            def OnZValue(event):
                Obj = event.GetEventObject()
                try:
                    value = int(Obj.GetValue())
                except ValueError:
                    value = data[ObjIndx[Obj.GetId()]]
                data[ObjIndx[Obj.GetId()]] = value
                self.PlotImage()
            
            zSizer = wx.FlexGridSizer(0,4,5,5)
            zSizer.Add(wx.StaticText(self.SCCDPanel,label='Upper intensity mask:'),0,wx.ALIGN_CENTER_VERTICAL)
            zMax = wx.TextCtrl(self.SCCDPanel,value='%d'%(data['Zmax']))
            zMax.Bind(wx.EVT_TEXT_ENTER,OnZValue)
            zMax.Bind(wx.EVT_KILL_FOCUS,OnZValue)
            ObjIndx[zMax.GetId()] = 'Zmax'
            zSizer.Add(zMax)
            zSizer.Add(wx.StaticText(self.SCCDPanel,label='Intensity subtraction:'),0,wx.ALIGN_CENTER_VERTICAL)
            zMin = wx.TextCtrl(self.SCCDPanel,value='%d'%(data['Zmin']))
            ObjIndx[zMin.GetId()] = 'Zmin'
            zMin.Bind(wx.EVT_TEXT_ENTER,OnZValue)
            zMin.Bind(wx.EVT_KILL_FOCUS,OnZValue)
            zSizer.Add(zMin)
            return zSizer
                        
        def ZeroSizer(data):
            
            def OnZeroValue(event):
                Obj = event.GetEventObject()
                item = ObjIndx[Obj.GetId()]
                try:
                    value = float(Obj.GetValue())
                except ValueError:
                    value = data[item[0]][item[1]]
                data[item[0]][item[1]] = value
                Obj.SetValue('%.3f'%(value))
                self.PlotImage()
                
            def OnZpdgValue(event):
                Obj = event.GetEventObject()
                item = ObjIndx[Obj.GetId()]
                try:
                    value = float(Obj.GetValue())
                except ValueError:
                    value = self.data[item[0]]
                self.data[item[0]] = value
                self.data['radius'] = 180.0*value/math.pi
                Obj.SetValue('%.3f'%(value))
                self.PlotImage()
            
            zeroSizer = wx.FlexGridSizer(0,6,5,5)
            zeroSizer.Add(wx.StaticText(self.SCCDPanel,label='X-zero:'),0,wx.ALIGN_CENTER_VERTICAL)
            zMax = wx.TextCtrl(self.SCCDPanel,value='%.3f'%(data['Zeros'][0]))
            zMax.Bind(wx.EVT_TEXT_ENTER,OnZeroValue)
            zMax.Bind(wx.EVT_KILL_FOCUS,OnZeroValue)
            ObjIndx[zMax.GetId()] = ['Zeros',0]
            zeroSizer.Add(zMax)
            zeroSizer.Add(wx.StaticText(self.SCCDPanel,label='Y-zero:'),0,wx.ALIGN_CENTER_VERTICAL)
            zMin = wx.TextCtrl(self.SCCDPanel,value='%.3f'%(data['Zeros'][1]))
            ObjIndx[zMin.GetId()] = ['Zeros',1]
            zMin.Bind(wx.EVT_TEXT_ENTER,OnZeroValue)
            zMin.Bind(wx.EVT_KILL_FOCUS,OnZeroValue)
            zeroSizer.Add(zMin)
            zeroSizer.Add(wx.StaticText(self.SCCDPanel,label='mm per deg:'),0,wx.ALIGN_CENTER_VERTICAL)
            zpdeg = wx.TextCtrl(self.SCCDPanel,value='%.3f'%(data['mm/deg']))
            ObjIndx[zpdeg.GetId()] = ['mm/deg']
            zpdeg.Bind(wx.EVT_TEXT_ENTER,OnZpdgValue)
            zpdeg.Bind(wx.EVT_KILL_FOCUS,OnZpdgValue)
            zeroSizer.Add(zpdeg)
            return zeroSizer
                        
        def TBLRSizer(data):
            
            def OnTBLRValue(event):
                Obj = event.GetEventObject()
                item = ObjIndx[Obj.GetId()]
                try:
                    value = int(Obj.GetValue())
                except ValueError:
                    value = data[item[0]][item[1]]
                data[item[0]][item[1]] = value
                Obj.SetValue('%d'%(value))
                self.PlotImage()
            
            TBLRsizer = wx.FlexGridSizer(0,4,5,5)
            for i,item in enumerate(['Bottom','Top']): 
                TBLRsizer.Add(wx.StaticText(self.SCCDPanel,label=item+' limit, pixels:'),0,wx.ALIGN_CENTER_VERTICAL)
                TBlim = wx.TextCtrl(self.SCCDPanel,value='%d'%(data['TBlimits'][i]))
                TBlim.Bind(wx.EVT_TEXT_ENTER,OnTBLRValue)
                TBlim.Bind(wx.EVT_KILL_FOCUS,OnTBLRValue)
                ObjIndx[TBlim.GetId()] = ['TBlimits',i]
                TBLRsizer.Add(TBlim)
            return TBLRsizer
            
        def ScanSizer(data):
            
            def OnScanValue(event):
                Obj = event.GetEventObject()
                item = ObjIndx[Obj.GetId()][0]
                try:
                    value = float(Obj.GetValue())
                except ValueError:
                    value = self.data['2thScan'][item]
                data['2thScan'][item] = value
                Obj.SetValue('%.3f'%(value))
                if item in [0,1]:
                    pixel = data['pixel']
                    zero = data['Zeros'][0]
                    tthscale = data['mm/deg']
                    npixel = (value*tthscale+zero)*1000/pixel[0]
                    data['LRlimits'][item] = npixel
                self.PlotImage()
            
            scanSizer = wx.FlexGridSizer(0,6,5,5)
            for i,item in enumerate(['Lower 2-th','Upper 2-th','2-th step']):
                scanSizer.Add(wx.StaticText(self.SCCDPanel,label=item+':'),0,wx.ALIGN_CENTER_VERTICAL)
                scanParm = wx.TextCtrl(self.SCCDPanel,value='%.3f'%(data['2thScan'][i]))
                scanParm.Bind(wx.EVT_TEXT_ENTER,OnScanValue)
                scanParm.Bind(wx.EVT_KILL_FOCUS,OnScanValue)
                ObjIndx[scanParm.GetId()] = [i]
                scanSizer.Add(scanParm)
            return scanSizer                
                        
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(ColorSizer(self.data),0,wx.ALIGN_CENTER_VERTICAL)        
        mainSizer.Add(MaxSizer(self.data),0,wx.ALIGN_LEFT|wx.EXPAND)
        mainSizer.Add(ZSizer(self.data),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(ZeroSizer(self.data),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(ScanSizer(self.data),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(TBLRSizer(self.data),0,wx.ALIGN_CENTER_VERTICAL)
        self.SCCDPanel.SetSizer(mainSizer)
        mainSizer.Layout()    
        fitSize = mainSizer.Fit(self.SCCDPanel)
#        self.SCCDPanel.GetParent().SetSize(fitSize)
            
    def PlotImage(self):
        pixel = self.data['pixel']
        scalex = pixel[0]/1000.
        scaley = pixel[1]/1000.

        
        def OnMotion(event):
            Page.canvas.SetToolTipString('')
            sizexy = self.data['size']
            if event.xdata and event.ydata:                 #avoid out of frame errors
                Page.canvas.SetCursor(wx.CROSS_CURSOR)
                item = self.itemPicked
                if item:
                    if 'Line' in str(item):
                        Page.canvas.SetToolTipString('%8.3f %8.3fmm'%(event.xdata,event.ydata))
                else:
                    xpos = event.xdata
                    ypos = event.ydata
                    xpix = xpos/(self.data['imScale']*scalex)
                    ypix = sizexy[1]-ypos/(self.data['imScale']*scaley)
                    Int = 0
                    if (0 <= xpix <= sizexy[0]) and (0 <= ypix <= sizexy[1]):
                        Int = self.Image[xpix][ypix]-self.data['Zmin']
                        tth = (xpos-self.data['Zeros'][0])/self.data['mm/deg']
                        vecA = np.array([1.0,0,0])
                        vecB = np.array([npcosd(tth)*self.data['radius'],npsind(tth)*self.data['radius'],(ypos-self.data['Zeros'][1])])
                        vecB /= nl.norm(vecB)
                        tth2 = npacosd(np.dot(vecA,vecB))*tth/abs(tth)
                        self.plotNB.status.SetFields(\
                            ['','Detector x,y =%9.3fmm %9.3fmm, 2-th =%9.3f I = %6d'%(xpos,ypos,tth2,Int)])
                            
        def OnPick(event):
            if self.itemPicked is not None: return
            self.itemPicked = event.artist
            self.mousePicked = event.mouseevent
            
        def OnRelease(event):
            if self.itemPicked is None: return
            xpos,ypos = [event.xdata,event.ydata]
            if '_line0' in str(self.itemPicked):    #X-zero
                self.data['Zeros'][0] = xpos               
            elif '_line1' in str(self.itemPicked):  #Y-zero
                self.data['Zeros'][1] = ypos
            elif '_line2' in str(self.itemPicked): #Y-lower limit
                self.data['TBlimits'][0] = int(ypos/scaley)
            elif '_line3' in str(self.itemPicked): #Y-upper limit
                self.data['TBlimits'][1] = int(ypos/scaley)
            elif '_line4' in str(self.itemPicked): #X-lower limit
                self.data['LRlimits'][0] = int(xpos/scalex)
            elif '_line5' in str(self.itemPicked): #X-upper limit
                self.data['LRlimits'][1] = int(xpos/scalex)
            self.itemPicked = None    
            self.PlotImage()
            self.UpdateControls(event)
        
        try:
            plotNum = self.plotNB.plotList.index('scanCCD image')
            Page = self.plotNB.nb.GetPage(plotNum)
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
            Page.figure.clf()
            Plot = Page.figure.gca()
            if not Page.IsShown():
                Page.Show()
        except ValueError:
            Plot = self.plotNB.addMpl('scanCCD image').gca()
            plotNum = self.plotNB.plotList.index('scanCCD image')
            Page = self.plotNB.nb.GetPage(plotNum)
            Page.canvas.mpl_connect('motion_notify_event', OnMotion)
            Page.canvas.mpl_connect('pick_event', OnPick)
            Page.canvas.mpl_connect('button_release_event', OnRelease)
            xylim = []
        xlim,ylim = self.data['size']
        Imin,Imax = [0,self.data['range'][1]]
        Xmax = scalex*xlim*self.data['imScale']
        Ymax = scaley*ylim*self.data['imScale']
        TBlimit = np.array(self.data['TBlimits'])*scalex
        LRlimit = np.array(self.data['LRlimits'])*scaley

        self.plotNB.status.SetFields(['',''])
        Zmin = self.data['Zmin']
        Zmax = self.data['Zmax']
        Zeros = self.data['Zeros']
        Lines = []
            
#        MA = ma.masked_greater(ma.masked_less(self.Image.T,Zmin),Zmax)
#        MaskA = ma.getmaskarray(MA)
#        ImgM = Plot.imshow(MaskA,aspect='auto',cmap='Reds',
#            interpolation='nearest',vmin=0,vmax=2,extent=[0,Xmax,0,Ymax])
        Img = Plot.imshow(self.Image.T-Zmin,interpolation='nearest',vmin=Imin,vmax=Imax,
            aspect='auto',cmap=self.data['color'],extent=[0,Xmax,0,Ymax])
        Lines.append(Plot.axvline(Zeros[0],color='b',dashes=(5,5),picker=3.))    
        Lines.append(Plot.axhline(Zeros[1],color='b',dashes=(5,5),picker=3.))    
        Lines.append(Plot.axhline(TBlimit[0],color='g',dashes=(5,5),picker=3.))    
        Lines.append(Plot.axhline(TBlimit[1],color='r',dashes=(5,5),picker=3.))    
        Lines.append(Plot.axvline(LRlimit[0],color='g',dashes=(5,5),picker=3.))    
        Lines.append(Plot.axvline(LRlimit[1],color='r',dashes=(5,5),picker=3.))
        for blk in range(len(self.imagefiles)):
            Lines.append(Plot.axvline(blk*4096*scalex,color='k'))    
        Plot.set_title('')
        Plot.set_xlabel('Scan, mm')
        Plot.set_ylabel('Detector Y, mm')
        if xylim:
            Page.toolbar.push_current()
            Plot.set_xlim(xylim[0])
            Plot.set_ylim(xylim[1])
            xylim = []
            Page.toolbar.push_current()
            Page.toolbar.draw()
        else:
            Page.canvas.draw()
            
    def PlotBlock(self,block,title):
        try:
            plotNum = self.plotNB.plotList.index('Block')
            Page = self.plotNB.nb.GetPage(plotNum)
            Plot = Page.figure.gca()          #get previous powder plot & get limits
            xylim = Plot.get_xlim(),Plot.get_ylim()
            Page.figure.clf()
            Plot = Page.figure.gca()
            if not Page.IsShown():
                Page.Show()
        except ValueError:
            Plot = self.plotNB.addMpl('Block').gca()
            plotNum = self.plotNB.plotList.index('Block')
            Page = self.plotNB.nb.GetPage(plotNum)
            xylim = []
        Img = Plot.imshow(block,interpolation='nearest',aspect='equal',cmap=self.data['color'])
        Plot.invert_yaxis()
        Plot.set_title(title)
        if xylim:
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
                    self.plotNB.status.SetStatusText('X =%9.3f %s =%9.3f'%(xpos,type,ypos),1)                   
                except TypeError:
                    self.plotNB.status.SetStatusText('Select '+type+' pattern first',1)
    
        try:
            plotNum = self.plotNB.plotList.index(type)
            Page = self.plotNB.nb.GetPage(plotNum)
            if not newPlot:
                Plot = Page.figure.gca()
                xylim = Plot.get_xlim(),Plot.get_ylim()
            Page.figure.clf()
            Plot = Page.figure.gca()
        except ValueError:
            newPlot = True
            Plot = self.plotNB.addMpl(type).gca()
            plotNum = self.plotNB.plotList.index(type)
            Page = self.plotNB.nb.GetPage(plotNum)
            Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        
        Page.SetFocus()
        self.plotNB.status.DestroyChildren()
        Plot.set_title(type)
        if type == 'line scan':
            Plot.set_xlabel(r'image x, mm',fontsize=14)
        else:
            Plot.set_xlabel(r'2-theta, deg',fontsize=14)            
        Plot.set_ylabel(r''+type,fontsize=14)
        colors=['b','g','r','c','m','k']
        lenX = 0
        X,Y,S = XY
        Plot.plot(X,Y,'k',picker=False)
        Plot.plot(X,S,'r',picker=False)
        if not newPlot:
            Page.toolbar.push_current()
            Plot.set_xlim(xylim[0])
            Plot.set_ylim(xylim[1])
            xylim = []
            Page.toolbar.push_current()
            Page.toolbar.draw()
        else:
            Page.canvas.draw()

class scanCCDmain(wx.App):
    def OnInit(self):
        self.main = scanCCD(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    'starts main application to merge data from scanning CCD'
    application = scanCCDmain(0)
    application.MainLoop()
    
if __name__ == '__main__':
    main()
