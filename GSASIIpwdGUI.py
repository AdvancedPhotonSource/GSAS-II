# -*- coding: utf-8 -*-
#GSASIIpwdGUI - powder data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIpwdGUI: Powder Pattern GUI routines*
-------------------------------------------

Used to define GUI controls for the routines that interact
with the powder histogram (PWDR) data tree items.

'''
import sys
import os.path
import wx
import wx.grid as wg
import wx.lib.scrolledpanel as wxscroll
import numpy as np
import numpy.ma as ma
import math
import copy
import random as ran
import cPickle
import scipy.interpolate as si
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import GSASIIIO as G2IO
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIindex as G2indx
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIctrls as G2G
import GSASIIElemGUI as G2elemGUI
import GSASIIElem as G2elem
import GSASIIsasd as G2sasd
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WACV = wx.ALIGN_CENTER_VERTICAL
GkDelta = unichr(0x0394)
Pwr10 = unichr(0x0b9)+unichr(0x2070)
Pwr20 = unichr(0x0b2)+unichr(0x2070)
Pwrm1 = unichr(0x207b)+unichr(0x0b9)
Pwrm2 = unichr(0x207b)+unichr(0x0b2)
Pwrm6 = unichr(0x207b)+unichr(0x2076)
Pwrm4 = unichr(0x207b)+unichr(0x2074)
Angstr = unichr(0x00c5)   
# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
    
################################################################################
###### class definitions
################################################################################

class RDFDialog(wx.Dialog):
    def __init__(self,parent):
        wx.Dialog.__init__(self,parent,-1,'Background radial distribution function',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.result = {'UseObsCalc':True,'maxR':20.0,'Smooth':'linear'}
        
        self.Draw()
        
    def Draw(self):
        
        def OnUseOC(event):
            self.result['UseObsCalc'] = not self.result['UseObsCalc']
            
        def OnSmCombo(event):
            self.result['Smooth'] = smCombo.GetValue()
                    
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label='Background RDF controls:'),0,WACV)
        useOC = wx.CheckBox(self.panel,label=' Use obs - calc intensities?')
        useOC.SetValue(self.result['UseObsCalc'])
        useOC.Bind(wx.EVT_CHECKBOX,OnUseOC)
        mainSizer.Add(useOC,0,WACV)
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        dataSizer.Add(wx.StaticText(self.panel,label=' Smoothing type: '),0,WACV)
        smChoice = ['linear','nearest',]
        smCombo = wx.ComboBox(self.panel,value=self.result['Smooth'],choices=smChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        smCombo.Bind(wx.EVT_COMBOBOX, OnSmCombo)
        dataSizer.Add(smCombo,0,WACV)
        dataSizer.Add(wx.StaticText(self.panel,label=' Maximum radial dist.: '),0,WACV)
        maxR = G2G.ValidatedTxtCtrl(self.panel,self.result,'maxR',nDig=(10,1),min=10.,max=50.,
            typeHint=float)
        dataSizer.Add(maxR,0,WACV)
        mainSizer.Add(dataSizer,0,WACV)

        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(self.panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)
        
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        
    def GetSelection(self):
        return self.result

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)
        
    
################################################################################
##### Setup routines
################################################################################
    
def IsHistogramInAnyPhase(G2frame,histoName):
    'Needs a doc string'
    phases = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
    if phases:
        item, cookie = G2frame.PatternTree.GetFirstChild(phases)
        while item:
            data = G2frame.PatternTree.GetItemPyData(item)
            histoList = data['Histograms'].keys()
            if histoName in histoList:
                return G2frame.PatternTree.GetItemText(item)
            item, cookie = G2frame.PatternTree.GetNextChild(phases, cookie)
        return False
    else:
        return False

def SetDefaultSample():
    'Fills in default items for the Sample dictionary for Debye-Scherrer & SASD'
    return {
        'InstrName':'',
        'ranId':ran.randint(0,sys.maxint),
        'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],
        'DisplaceX':[0.0,False],'DisplaceY':[0.0,False],'Diffuse':[],
        'Temperature':300.,'Pressure':0.1,'Time':0.0,
        'FreePrm1':0.,'FreePrm2':0.,'FreePrm3':0.,
        'Gonio. radius':200.0,
        'Omega':0.0,'Chi':0.0,'Phi':0.0,'Azimuth':0.0,
#SASD items
        'Materials':[{'Name':'vacuum','VolFrac':1.0,},{'Name':'vacuum','VolFrac':0.0,}],
        'Thick':1.0,'Contrast':[0.0,0.0],       #contrast & anomalous contrast
        'Trans':1.0,                            #measured transmission
        'SlitLen':0.0,                          #Slit length - in Q(A-1)
        }
def SetupSampleLabels(histName,dataType,histType):
    '''Setup a list of labels and number formatting for use in
    labeling sample parameters.
    :param str histName: Name of histogram, ("PWDR ...")
    :param str dataType: 
    '''
    parms = []
    parms.append(['Scale','Histogram scale factor: ',[10,7]])
    if 'C' in histType:
        parms.append(['Gonio. radius','Goniometer radius (mm): ',[10,3]])
    if 'PWDR' in histName:
        if dataType == 'Debye-Scherrer':
            if 'T' in histType:
                parms += [['Absorption',u'Sample absorption (\xb5\xb7r/l): ',[10,4]],]
            else:
                parms += [['DisplaceX',u'Sample X displ. perp. to beam (\xb5m): ',[10,3]],
                    ['DisplaceY',u'Sample Y displ. || to beam (\xb5m): ',[10,3]],
                    ['Absorption',u'Sample absorption (\xb5\xb7r): ',[10,4]],]
        elif dataType == 'Bragg-Brentano':
            parms += [['Shift',u'Sample displacement(\xb5m): ',[10,4]],
                ['Transparency',u'Sample transparency(1/\xb5eff, cm): ',[10,3]],
                ['SurfRoughA','Surface roughness A: ',[10,4]],
                ['SurfRoughB','Surface roughness B: ',[10,4]]]
    elif 'SASD' in histName:
        parms.append(['Thick','Sample thickness (mm)',[10,3]])
        parms.append(['Trans','Transmission (meas)',[10,3]])
        parms.append(['SlitLen',u'Slit length (Q,\xc5'+Pwrm1+')',[10,3]])
    parms.append(['Omega','Goniometer omega:',[10,3]])
    parms.append(['Chi','Goniometer chi:',[10,3]])
    parms.append(['Phi','Goniometer phi:',[10,3]])
    parms.append(['Azimuth','Detector azimuth:',[10,3]])
    parms.append(['Time','Clock time (s):',[12,3]])
    parms.append(['Temperature','Sample temperature (K): ',[10,3]])
    parms.append(['Pressure','Sample pressure (MPa): ',[10,3]])
    return parms

def SetDefaultSASDModel():
    'Fills in default items for the SASD Models dictionary'    
    return {'Back':[0.0,False],'Size':{'MinDiam':50,'MaxDiam':10000,'Nbins':100,'logBins':True,'Method':'MaxEnt','Distribution':[],
        'Shape':['Spheroid',1.0],'MaxEnt':{'Niter':100,'Precision':0.01,'Sky':-3},
        'IPG':{'Niter':100,'Approach':0.8,'Power':-1},'Reg':{},},            
        'Particle':{'Matrix':{'Name':'vacuum','VolFrac':[0.0,False]},'Levels':[],},
        'Current':'Size dist.','BackFile':'',
        }
        
def SetDefaultREFDModel():
    '''Fills in default items for the REFD Models dictionary
    Defined as follows for each layer:
        Name: name of substance
        Thick: thickness of layer in Angstroms (not present for top & bottom layers)
        Rough: upper surface roughness for layer (not present for toplayer)
        Penetration: mixing of layer substance into layer above-is this needed?
        DenMul: multiplier for layer scattering density (default = 1.0)
    Top layer defaults to vacuum (or air/any gas); can be substituted for some other substance
    Bottom layer default: infinitely thisck Silicon; can be substituted for some other substance
    '''
    return {'Layers':[{'Name':'vacuum','DenMul':[1.0,False],},                                  #top layer
        {'Name':'vacuum','Rough':[0.,False],'Penetration':[0.,False],'DenMul':[1.0,False],}],   #bottom layer
        'Scale':[1.0,False],'FltBack':[0.0,False],'Zero':'Top','dQ type':'None',                #globals
        'Minimizer':'LMLS','Resolution':[0.,'Const dq/q'],'Recomb':0.5,'Toler':0.5,             #minimizer controls
        'DualFitFiles':['',],'DualFltBacks':[[0.0,False],],'DualScales':[[1.0,False],]}         #optional stuff for multidat fits?
        
def SetDefaultSubstances():
    'Fills in default items for the SASD Substances dictionary'
    return {'Substances':{'vacuum':{'Elements':{},'Volume':1.0,'Density':0.0,'Scatt density':0.0,'XImag density':0.0},
        'unit scatter':{'Elements':None,'Volume':None,'Density':None,'Scatt density':1.0,'XImag density':1.0}}}

def GetFileList(G2frame,fileType):
    fileList = []
    id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
    while id:
        name = G2frame.PatternTree.GetItemText(id)
        if fileType in name.split()[0]:
            fileList.append(name)
        id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    return fileList
        
def GetHistsLikeSelected(G2frame):
    '''Get the histograms that match the current selected one:
    The histogram prefix and data type (PXC etc.), the number of
    wavelengths and the instrument geometry (Debye-Scherrer etc.) 
    must all match. The current histogram is not included in the list. 

    :param wx.Frame G2frame: pointer to main GSAS-II data tree
    '''
    histList = []
    inst,inst2 = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(
            G2frame,G2frame.PatternId, 'Instrument Parameters')
        )
    hType = inst['Type'][0]
    if 'Lam1' in inst:
        hLam = 2
    elif 'Lam' in inst:
        hLam = 1
    else:
        hLam = 0
    sample = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(
            G2frame,G2frame.PatternId, 'Sample Parameters')
        )
    hGeom = sample.get('Type')
    hstName = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    hPrefix = hstName.split()[0]+' '
    # cycle through tree looking for items that match the above
    item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)        
    while item:
        name = G2frame.PatternTree.GetItemText(item)
        if name.startswith(hPrefix) and name != hstName:
            cGeom,cType,cLam, = '?','?',-1
            subitem, subcookie = G2frame.PatternTree.GetFirstChild(item)
            while subitem:
                subname = G2frame.PatternTree.GetItemText(subitem)
                if subname == 'Sample Parameters':
                    sample = G2frame.PatternTree.GetItemPyData(subitem)
                    cGeom = sample.get('Type')
                elif subname == 'Instrument Parameters':
                    inst,inst2 = G2frame.PatternTree.GetItemPyData(subitem)
                    cType = inst['Type'][0]
                    if 'Lam1' in inst:
                        cLam = 2
                    elif 'Lam' in inst:
                        cLam = 1
                    else:
                        cLam = 0
                subitem, subcookie = G2frame.PatternTree.GetNextChild(item, subcookie)
            if cLam == hLam and cType == hType and cGeom == hGeom:
                if name not in histList: histList.append(name)
        item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    return histList

def SetCopyNames(histName,dataType,addNames=[]):
    '''Determine the items in the sample parameters that should be copied,
    depending on the histogram type and the instrument type.
    '''
    copyNames = ['Scale',]
    histType = 'HKLF'
    if 'PWDR' in histName:
        histType = 'PWDR'
        if 'Debye' in dataType:
            copyNames += ['DisplaceX','DisplaceY','Absorption']
        else:       #Bragg-Brentano
            copyNames += ['Shift','Transparency','SurfRoughA','SurfRoughB']
    elif 'SASD' in histName:
        histType = 'SASD'
        copyNames += ['Materials','Thick',]
    if len(addNames):
        copyNames += addNames
    return histType,copyNames
    
def CopyPlotCtrls(G2frame):
    '''Global copy: Copy plot controls from current histogram to others.
    '''
    hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    histList = GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No other histograms match '+hst,G2frame.dataFrame)
        return
    sourceData = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)
    
    if 'Offset' not in sourceData[0]:    #patch for old data
        sourceData[0].update({'Offset':[0.0,0.0],'delOffset':0.02,'refOffset':-1.0,
            'refDelt':0.01,'qPlot':False,'dPlot':False,'sqrtPlot':False})
        G2frame.PatternTree.SetItemPyData(G2frame.PatternId,sourceData)
        
    dlg = G2G.G2MultiChoiceDialog(
        G2frame.dataFrame, 
        'Copy plot controls from\n'+str(hst[5:])+' to...',
        'Copy plot controls', histList)
    results = []
    try:
        if dlg.ShowModal() == wx.ID_OK:
            results = dlg.GetSelections()
    finally:
        dlg.Destroy()
    copyList = []
    for i in results: 
        copyList.append(histList[i])

    keys = ['Offset','delOffset','refOffset','refDelt','qPlot','dPlot','sqrtPlot']
    source = dict(zip(keys,[sourceData[0][item] for item in keys]))
    for hist in copyList:
        Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,hist)
        data = G2frame.PatternTree.GetItemPyData(Id)
        data[0].update(source)
        G2frame.PatternTree.SetItemPyData(Id,data)
    print 'Copy of plot controls successful'

def CopySelectedHistItems(G2frame):
    '''Global copy: Copy items from current histogram to others.
    '''
    hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    histList = GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No other histograms match '+hst,G2frame.dataFrame)
        return
    choices = ['Limits','Background','Instrument Parameters','Sample Parameters']
    dlg = G2G.G2MultiChoiceDialog(
        G2frame.dataFrame, 
        'Copy which histogram sections from\n'+str(hst[5:]),
        'Select copy sections', choices, filterBox=False)
    dlg.SetSelections(range(len(choices)))
    choiceList = []
    if dlg.ShowModal() == wx.ID_OK:
        choiceList = [choices[i] for i in dlg.GetSelections()]
    if not choiceList: return
    
    dlg = G2G.G2MultiChoiceDialog(
        G2frame.dataFrame, 
        'Copy parameters from\n'+str(hst[5:])+' to...',
        'Copy parameters', histList)
    results = []
    try:
        if dlg.ShowModal() == wx.ID_OK:
            results = dlg.GetSelections()
    finally:
        dlg.Destroy()
    copyList = []
    for i in results: 
        copyList.append(histList[i])

    if 'Limits' in choiceList: # Limits
        data = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Limits'))
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.SetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Limits'),
                copy.deepcopy(data))
    if 'Background' in choiceList:  # Background
        data = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Background'))
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.SetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Background'),
                copy.deepcopy(data))
    if 'Instrument Parameters' in choiceList:  # Instrument Parameters
        # for now all items in Inst. parms are copied
        data,data1 = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(
                G2frame,G2frame.PatternId,'Instrument Parameters'))
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters')
                )[0].update(copy.deepcopy(data))
            G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters')
                )[1].update(copy.deepcopy(data1))
    if 'Sample Parameters' in choiceList:  # Sample Parameters
        data = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(
                G2frame,G2frame.PatternId,'Sample Parameters'))
        # selects items to be copied
        histType,copyNames = SetCopyNames(hst,data['Type'],
            addNames = ['Omega','Chi','Phi','Gonio. radius','InstrName'])
        copyDict = {parm:data[parm] for parm in copyNames}
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters')
                ).update(copy.deepcopy(copyDict))
                         
################################################################################
#####  Powder Peaks
################################################################################           
       
def UpdatePeakGrid(G2frame, data):
    '''respond to selection of PWDR powder peaks data tree item.
    '''
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
        
    def OnAutoSearch(event):
        PatternId = G2frame.PatternId
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        profile = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        x0 = profile[0]
        iBeg = np.searchsorted(x0,limits[0])
        iFin = np.searchsorted(x0,limits[1])
        x = x0[iBeg:iFin]
        y0 = profile[1][iBeg:iFin]
        ysig = 0.5*np.std(y0)
        offset = [-1,1]
        ymask = ma.array(y0,mask=(y0<ysig))
        for off in offset:
            ymask = ma.array(ymask,mask=(ymask-np.roll(y0,off)<=0.))
        indx = ymask.nonzero()
        mags = ymask[indx]
        poss = x[indx]
        refs = zip(poss,mags)
        if 'C' in Inst['Type'][0]:    
            refs = G2mth.sortArray(refs,0,reverse=True)     #small 2-Thetas first
        else:   #'T'OF
            refs = G2mth.sortArray(refs,0,reverse=False)    #big TOFs first
        for i,ref1 in enumerate(refs):      #reject picks closer than 1 FWHM
            for ref2 in refs[i+1:]:
                if abs(ref2[0]-ref1[0]) < G2pwd.getFWHM(ref1[0],inst):
                    del(refs[i])
        if 'C' in Inst['Type'][0]:    
            refs = G2mth.sortArray(refs,1,reverse=True)
        else:   #'T'OF
            refs = G2mth.sortArray(refs,1,reverse=False)
        for pos,mag in refs:
            data['peaks'].append(G2mth.setPeakparms(inst,inst2,pos,mag))
        UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        
    def OnCopyPeaks(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy peak list from\n'+str(hst[5:])+' to...',
            'Copy peaks', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.SetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Peak List'),copy.deepcopy(data))
    
    def OnUnDo(event):
        DoUnDo()
        G2frame.dataFrame.UnDo.Enable(False)
        
    def DoUnDo():
        print 'Undo last refinement'
        file = open(G2frame.undofile,'rb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item),cPickle.load(file))
            if G2frame.dataDisplay.GetName() == item:
                if item == 'Background':
                    UpdateBackground(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
                elif item == 'Instrument Parameters':
                    UpdateInstrumentGrid(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
                elif item == 'Peak List':
                    UpdatePeakGrid(G2frame,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, item)))
            print item,' recovered'
        file.close()
        
    def SaveState():
        G2frame.undofile = os.path.join(G2frame.dirname,'GSASII.save')
        file = open(G2frame.undofile,'wb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            cPickle.dump(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,item)),file,1)
        file.close()
        G2frame.dataFrame.UnDo.Enable(True)
        
    def OnLSQPeakFit(event):
        if G2frame.dataDisplay.IsCellEditControlEnabled(): # complete any grid edits in progress
            G2frame.dataDisplay.HideCellEditControl()
            G2frame.dataDisplay.DisableCellEditControl()
        if not G2frame.GSASprojectfile:            #force a save of the gpx file so SaveState can write in the same directory
            G2frame.OnFileSaveas(event)
        wx.CallAfter(OnPeakFit,'LSQ')
        
    def OnOneCycle(event):
        if G2frame.dataDisplay.IsCellEditControlEnabled(): # complete any grid edits in progress
            G2frame.dataDisplay.HideCellEditControl()
            G2frame.dataDisplay.DisableCellEditControl()
        wx.CallAfter(OnPeakFit,'LSQ',oneCycle=True)
        
    def OnSeqPeakFit(event):
        histList = G2gd.GetPatternTreeDataNames(G2frame,['PWDR',])
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame, 'Sequential peak fits',
             'Select dataset to include',histList)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(histList[sel])
        dlg.Destroy()
        if not names:
            return
        Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequential peak fit results')
        if Id:
            SeqResult = G2frame.PatternTree.GetItemPyData(Id)
        else:
            SeqResult = {}
            Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Sequential peak fit results')
        SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
        SeqResult['histNames'] = names
        Reverse = False
        CopyForward = False
        choice = ['Reverse sequence','Copy from prev.',]
        dlg = wx.MultiChoiceDialog(G2frame.dataFrame,'Sequential controls','Select controls',choice)
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                if sel:
                    CopyForward = True
                else:
                    Reverse = True
        dlg.Destroy()
        dlg = wx.ProgressDialog('Sequential peak fit','Data set name = '+names[0],len(names), 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        controls = {'deriv type':'analytic','min dM/M':0.0001,}
        print 'Peak Fitting with '+controls['deriv type']+' derivatives:'
        oneCycle = False
        FitPgm = 'LSQ'
        prevVaryList = []
        peaks = None
        varyList = None
        if Reverse:
            names.reverse()
        try:
            for i,name in enumerate(names):
                print ' Sequential fit for ',name
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    dlg.Destroy()
                    break
                PatternId =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name)
                if i and CopyForward:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'),copy.deepcopy(peaks))
                    prevVaryList = varyList[:]
                peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
                background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
                limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
                inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
                data = G2frame.PatternTree.GetItemPyData(PatternId)[1]
                peaks['sigDict'],result,sig,Rvals,varyList,parmDict,fullvaryList,badVary = G2pwd.DoPeakFit(FitPgm,peaks['peaks'],
                    background,limits,inst,inst2,data,prevVaryList,oneCycle,controls)
                if len(result[0]) != len(fullvaryList):
                    dlg.Destroy()
                    print ' ***** Sequential peak fit stopped at '+name+' *****'
                    break
                else:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'),copy.deepcopy(peaks))
                    SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                        'covMatrix':np.eye(len(result[0])),'title':name,'parmDict':parmDict,
                        'fullVary':fullvaryList,'badVary':badVary}
            print ' ***** Sequential peak fit successful *****'
        finally:
            dlg.Destroy()
        SeqResult['histNames'] = histList
        G2frame.PatternTree.SetItemPyData(Id,SeqResult)
        G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
        G2frame.PatternTree.SelectItem(Id)
        
    def OnClearPeaks(event):
        dlg = wx.MessageDialog(G2frame,'Delete all peaks?','Clear peak list',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = {'peaks':[],'sigDict':{}}
        finally:
            dlg.Destroy()
        UpdatePeakGrid(G2frame,peaks)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        
    def OnPeakFit(FitPgm,oneCycle=False):
        SaveState()
        controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
        if not controls:
            controls = {'deriv type':'analytic','min dM/M':0.0001,}     #fill in defaults if needed
        print 'Peak Fitting with '+controls['deriv type']+' derivatives:'
        PatternId = G2frame.PatternId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        data = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        dlg = wx.ProgressDialog('Residual','Peak fit Rwp = ',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            dlg.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            peaks['sigDict'] = G2pwd.DoPeakFit(FitPgm,peaks['peaks'],background,limits,inst,inst2,data,[],oneCycle,controls,dlg)[0]
        finally:
            print 'finished'
            dlg.Destroy()
        newpeaks = copy.copy(peaks)
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'),newpeaks)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        wx.CallAfter(UpdatePeakGrid,G2frame,newpeaks)
        
    def OnResetSigGam(event):
        PatternId = G2frame.PatternId
        Inst,Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks['peaks']:
            G2frame.ErrorDialog('No peaks!','Nothing to do!')
            return
        newpeaks = {'peaks':[],'sigDict':{}}
        for peak in peaks['peaks']:
            newpeaks['peaks'].append(G2mth.setPeakparms(Inst,Inst2,peak[0],peak[2]))
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'),newpeaks)
        UpdatePeakGrid(G2frame,newpeaks)
                
    def RefreshPeakGrid(event):
        
        event.StopPropagation()
        data['peaks'] = G2frame.PeakTable.GetData()
        T = []
        for peak in data['peaks']:T.append(peak[0])
        D = dict(zip(T,data['peaks']))
        T.sort()
        X = []
        for key in T: X.append(D[key])
        data['peaks'] = X        
        
    def setBackgroundColors():
       for r in range(G2frame.dataDisplay.GetNumberRows()):
           for c in range(G2frame.dataDisplay.GetNumberCols()):
               if G2frame.dataDisplay.GetColLabelValue(c) in ['position','intensity','alpha','beta','sigma','gamma']:
                   if float(G2frame.dataDisplay.GetCellValue(r,c)) < 0.:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.RED)
                   else:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.WHITE)
                                                  
    def KeyEditPeakGrid(event):
        '''Respond to pressing a key to act on selection of a row, column or cell
        in the Peak List table
        '''
        rowList = G2frame.dataDisplay.GetSelectedRows()
        colList = G2frame.dataDisplay.GetSelectedCols()
        selectList = G2frame.dataDisplay.GetSelectedCells()
        data = G2frame.PatternTree.GetItemPyData(G2frame.PickId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif rowList and (event.GetKeyCode() == wx.WXK_DELETE or event.GetKeyCode() == 8):
            # pressing the delete key or backspace deletes selected peak(s)
            G2frame.dataDisplay.ClearSelection()
            G2frame.dataDisplay.ClearGrid()
            rowList.sort()
            rowList.reverse()
            nDel = 0
            for row in rowList:
                G2frame.PeakTable.DeleteRow(row)
                nDel += 1
            if nDel:
                msg = wg.GridTableMessage(G2frame.PeakTable, 
                    wg.GRIDTABLE_NOTIFY_ROWS_DELETED,0,nDel)
                G2frame.dataDisplay.ProcessTableMessage(msg)
            data['peaks'] = G2frame.PeakTable.GetData()[:-nDel]
            G2frame.PatternTree.SetItemPyData(G2frame.PickId,data)
            wx.CallAfter(G2frame.dataDisplay.ForceRefresh)
            setBackgroundColors()
        elif colList and (event.GetKeyCode() == 89 or event.GetKeyCode() == 78):
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.PeakTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data['peaks'][row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data['peaks'][row][col]=False
        elif selectList and (event.GetKeyCode() == 89 or event.GetKeyCode() == 78):
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for row,col in selectList:
                if G2frame.PeakTable.GetTypeName(row,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        data['peaks'][row][col]=True
                    elif key == 78:  #'N'
                        data['peaks'][row][col]=False
        else:
            event.Skip()
            return
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        #wx.CallAfter(G2frame.dataDisplay.ForceRefresh) # did not always work
        wx.CallAfter(UpdatePeakGrid,G2frame,data)
            
    def SelectVars(rows):
        '''Set or clear peak refinement variables for peaks listed in rows
        '''
        refOpts = {G2frame.dataDisplay.GetColLabelValue(i):i+1 for i in range(G2frame.dataDisplay.GetNumberCols()) if G2frame.dataDisplay.GetColLabelValue(i) != "refine"}
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Select columns to refine',
            'Refinement Selection', sorted(refOpts.keys()),
            filterBox=False,toggle=False)
        sels = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sels = [sorted(refOpts.keys())[i] for i in dlg.GetSelections()]
            else:
                return
        finally:
            dlg.Destroy()
        for r in rows:
            for lbl,c in refOpts.iteritems():
                data['peaks'][r][c] = lbl in sels
        UpdatePeakGrid(G2frame,data)
        
    def OnRefineSelected(event):
        '''set refinement flags for the selected peaks
        '''
        rows = list(set([row for row,col in G2frame.dataDisplay.GetSelectedCells()] +
                        G2frame.dataDisplay.GetSelectedRows()))
        if not rows:
            wx.MessageBox('No selected rows. You must select rows or cells before using this command',
                          caption='No selected peaks')
            return
        SelectVars(rows)

    def OnRefineAll(event):
        '''set refinement flags for all peaks
        '''
        SelectVars(range(G2frame.dataDisplay.GetNumberRows()))

    def onCellListSClick(event):
        '''Called when a peak is selected so that it can be highlighted in the plot
        '''
        event.Skip()
        c =  event.GetRow(),event.GetCol()
        if c < 0: # replot except whan a column is selected
            wx.CallAfter(G2plt.PlotPatterns,G2frame,plotType='PWDR')
        
    def onCellListDClick(event):
        '''Called after a double-click on a cell label'''
        r,c =  event.GetRow(),event.GetCol()
        grid = G2frame.dataDisplay
        if r < 0 and c < 0:
            for row in range(grid.GetNumberRows()):
                grid.SelectRow(row,True)                    
            for col in range(grid.GetNumberCols()):
                grid.SelectCol(col,True)                    
        elif r > 0:     #row label: select it and replot!
            grid.ClearSelection()
            grid.SelectRow(r,True)
            wx.CallAfter(G2frame.dataDisplay.ForceRefresh)
            wx.CallAfter(G2plt.PlotPatterns,G2frame,plotType='PWDR')
        elif c > 0:     #column label: just select it (& redisplay)
            grid.ClearSelection()
            grid.SelectCol(c,True)
            if grid.GetColLabelValue(c) != 'refine': return
            choice = ['Y - vary all','N - vary none',]
            dlg = wx.SingleChoiceDialog(G2frame.dataFrame,
                                        'Select refinement option for '+grid.GetColLabelValue(c-1),
                                        'Refinement controls',choice)
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                if sel == 0:
                    for row in range(grid.GetNumberRows()): data['peaks'][row][c]=True
                else:
                    for row in range(grid.GetNumberRows()): data['peaks'][row][c]=False
            wx.CallAfter(UpdatePeakGrid,G2frame,data)
                 
    #======================================================================
    # beginning of UpdatePeakGrid init
    #======================================================================
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PeakMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    Status.SetStatusText('Global refine: select refine column & press Y or N')
    G2frame.Bind(wx.EVT_MENU, OnAutoSearch, id=G2gd.wxID_AUTOSEARCH)
    G2frame.Bind(wx.EVT_MENU, OnCopyPeaks, id=G2gd.wxID_PEAKSCOPY)
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_UNDO)
    G2frame.Bind(wx.EVT_MENU, OnRefineSelected, id=G2frame.dataFrame.peaksSel.GetId())
    G2frame.Bind(wx.EVT_MENU, OnRefineAll, id=G2frame.dataFrame.peaksAll.GetId())
    G2frame.Bind(wx.EVT_MENU, OnLSQPeakFit, id=G2gd.wxID_LSQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnOneCycle, id=G2gd.wxID_LSQONECYCLE)
    G2frame.Bind(wx.EVT_MENU, OnSeqPeakFit, id=G2gd.wxID_SEQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnClearPeaks, id=G2gd.wxID_CLEARPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnResetSigGam, id=G2gd.wxID_RESETSIGGAM)
    if data['peaks']:
        G2frame.dataFrame.AutoSearch.Enable(False)
        G2frame.dataFrame.PeakCopy.Enable(True)
        G2frame.dataFrame.PeakFit.Enable(True)
        G2frame.dataFrame.PFOneCycle.Enable(True)
        G2frame.dataFrame.SeqPeakFit.Enable(True)
    else:
        G2frame.dataFrame.PeakFit.Enable(False)
        G2frame.dataFrame.PeakCopy.Enable(False)
        G2frame.dataFrame.PFOneCycle.Enable(False)
        G2frame.dataFrame.AutoSearch.Enable(True)
        G2frame.dataFrame.SeqPeakFit.Enable(False)
    G2frame.PickTable = []
    rowLabels = []
    PatternId = G2frame.PatternId
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))[0]
    for i in range(len(data['peaks'])): rowLabels.append(str(i+1))
    if 'C' in Inst['Type'][0]:
        colLabels = ['position','refine','intensity','refine','sigma','refine','gamma','refine']
        Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
    else:
        colLabels = ['position','refine','intensity','refine','alpha','refine',
            'beta','refine','sigma','refine','gamma','refine']
        Types = [wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
    T = []
    for peak in data['peaks']:
        T.append(peak[0])
    D = dict(zip(T,data['peaks']))
    T.sort()
    if 'T' in Inst['Type'][0]:  #want big TOF's first
        T.reverse()
    X = []
    for key in T: X.append(D[key])
    data['peaks'] = X
    G2frame.PatternTree.SetItemPyData(G2frame.PickId,data)
    G2frame.PeakTable = G2G.Table(data['peaks'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Peak List')
    G2frame.dataFrame.currentGrids = []
    G2frame.dataDisplay = G2G.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.PeakTable, True)
    setBackgroundColors()                         
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshPeakGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, onCellListSClick)
#    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, onCellListSClick)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, onCellListDClick)
#    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, onCellListDClick)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([535,350])
    #G2frame.dataFrame.SendSizeEvent()

################################################################################
#####  Background
################################################################################           
       
def UpdateBackground(G2frame,data):
    '''respond to selection of PWDR background data tree item.
    '''
    if len(data) < 2:       #add Debye diffuse & peaks scattering here
        data.append({'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]})
    if 'nPeaks' not in data[1]:
        data[1].update({'nPeaks':0,'peaksList':[]})
    G2frame.dataFrame.currentGrids = []
    
    def OnBackFlagCopy(event):
        flag = data[0][1]
        backDict = data[-1]
        if backDict['nDebye']:
            DBflags = []
            for term in backDict['debyeTerms']:
                DBflags.append(term[1::2])
        if backDict['nPeaks']:
            PKflags = []
            for term in backDict['peaksList']:
                PKflags.append(term[1::2])            
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy bkg ref. flags from\n'+str(hst[5:])+' to...',
            'Copy bkg flags', histList)
        copyList = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections(): 
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            backData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Background'))
            backData[0][1] = copy.copy(flag)
            bkDict = backData[-1]
            if bkDict['nDebye'] == backDict['nDebye']:
                for i,term in enumerate(bkDict['debyeTerms']):
                    term[1::2] = copy.copy(DBflags[i])
            if bkDict['nPeaks'] == backDict['nPeaks']:
                for i,term in enumerate(bkDict['peaksList']):
                    term[1::2] = copy.copy(PKflags[i])                    
            
    def OnBackCopy(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame,'Copy bkg params from\n'+str(hst[5:])+' to...',
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.SetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Background'),copy.copy(data))

    def OnBkgFit(event):
        
        def SetInstParms(Inst):
            dataType = Inst['Type'][0]
            insVary = []
            insNames = []
            insVals = []
            for parm in Inst:
                insNames.append(parm)
                insVals.append(Inst[parm][1])
                if parm in ['U','V','W','X','Y','SH/L','I(L2)/I(L1)','alpha',
                    'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',] and Inst[parm][2]:
                        Inst[parm][2] = False
#                        insVary.append(parm)
            instDict = dict(zip(insNames,insVals))
            instDict['X'] = max(instDict['X'],0.01)
            instDict['Y'] = max(instDict['Y'],0.01)
            if 'SH/L' in instDict:
                instDict['SH/L'] = max(instDict['SH/L'],0.002)
            return dataType,instDict,insVary
    
        PatternId = G2frame.PatternId        
        controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
        background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        # sort the points for convenience and then separate them; extend the range if needed
        background[1]['FixedPoints'] = sorted(background[1]['FixedPoints'],key=lambda pair:pair[0])        
        X = [x for x,y in background[1]['FixedPoints']]
        Y = [y for x,y in background[1]['FixedPoints']]
        if X[0] > limits[0]:
            X = [limits[0]] + X
            Y = [Y[0]] + Y
        if X[-1] < limits[1]:
            X += [limits[1]]
            Y += [Y[-1]]
        # interpolate the fixed points onto the grid of data points within limits
        pwddata = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        xBeg = np.searchsorted(pwddata[0],limits[0])
        xFin = np.searchsorted(pwddata[0],limits[1])
        xdata = pwddata[0][xBeg:xFin]
        ydata = si.interp1d(X,Y)(ma.getdata(xdata))
        #GSASIIpath.IPyBreak()
        W = [1]*len(xdata)
        Z = [0]*len(xdata)

        # load instrument and background params
        print ' NB: Any instrument parameter refinement flags will be cleared'
        dataType,insDict,insVary = SetInstParms(inst)
        bakType,bakDict,bakVary = G2pwd.SetBackgroundParms(background)
        # how many background parameters are refined?
        if len(bakVary)*1.5 > len(X):
            msg = ("You are attempting to vary "+str(len(bakVary))+
                   " background terms with only "+str(len(X))+" background points"+
                    "\nAdd more points or reduce the number of terms")
            print msg
            G2frame.ErrorDialog('Too few points',msg)
            return
        
        wx.BeginBusyCursor()
        try:
            G2pwd.DoPeakFit('LSQ',[],background,limits,inst,inst2,
                np.array((xdata,ydata,W,Z,Z,Z)),bakVary,False,controls)
        finally:
            wx.EndBusyCursor()
        # compute the background values and plot them
        parmDict = {}
        bakType,bakDict,bakVary = G2pwd.SetBackgroundParms(background)
        parmDict.update(bakDict)
        parmDict.update(insDict)
        pwddata[3][xBeg:xFin] *= 0.
        pwddata[5][xBeg:xFin] *= 0.
        pwddata[4][xBeg:xFin] = G2pwd.getBackground(
            '',parmDict,bakType,dataType,xdata)[0]
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        # show the updated background values
        wx.CallLater(100,UpdateBackground,G2frame,data)
        
    def OnBkgClear(event):
        if 'FixedPoints' not in data[1]:
            return
        else:
            del data[1]['FixedPoints']
            G2plt.PlotPatterns(G2frame,plotType='PWDR')
    
    def OnPeaksMove(event):
        if not data[1]['nPeaks']:
            G2frame.ErrorDialog('Error','No peaks to move')
            return
        Peaks = {'peaks':[],'sigDict':{}}
        for peak in data[1]['peaksList']:
            Peaks['peaks'].append([peak[0],0,peak[2],0,peak[4],0,peak[6],0])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'),Peaks)
        
    def OnMakeRDF(event):
        dlg = RDFDialog(G2frame)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                RDFcontrols = dlg.GetSelection()
            else:
                return
        finally:
            dlg.Destroy()
        PatternId = G2frame.PatternId        
        background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        pwddata = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        auxPlot = G2pwd.MakeRDF(RDFcontrols,background,inst,pwddata)
#        GSASIIpath.IPyBreak()
        superMinusOne = unichr(0xaf)+unichr(0xb9)
        for plot in auxPlot:
            XY = np.array(plot[:2])
            if plot[2] == 'D(R)':
                xlabel = r'$R, \AA$'
                ylabel = r'$D(R), arb. units$'
            else:
                xlabel = r'$Q,\AA$'+superMinusOne
                ylabel = r'$I(Q)$'
            G2plt.PlotXY(G2frame,[XY,],Title=plot[2],labelX=xlabel,labelY=ylabel,lines=True)      
        
    def BackSizer():
        
        def OnNewType(event):
            data[0][0] = bakType.GetValue()
            
        def OnBakRef(event):
            data[0][1] = bakRef.GetValue()
            
        def OnBakTerms(event):
            data[0][2] = int(bakTerms.GetValue())
            M = len(data[0])
            N = data[0][2]+3
            item = data[0]
            if N > M:       #add terms
                for i in range(M,N): 
                    item.append(0.0)
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(item[-1])
            G2frame.PatternTree.SetItemPyData(BackId,data)
            #wx.CallAfter(UpdateBackground,G2frame,data)
            wx.CallLater(100,UpdateBackground,G2frame,data)
            
        backSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background function: '),0,WACV)
        bakType = wx.ComboBox(G2frame.dataDisplay,value=data[0][0],
                choices=Choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakType.Bind(wx.EVT_COMBOBOX, OnNewType)
        topSizer.Add(bakType)
        topSizer.Add((5,0),0)
        bakRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
        bakRef.SetValue(bool(data[0][1]))
        bakRef.Bind(wx.EVT_CHECKBOX, OnBakRef)
        topSizer.Add(bakRef,0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. coeff.: '),0,WACV)
        bakTerms = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[0][2]),choices=[str(i+1) for i in range(36)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakTerms.Bind(wx.EVT_COMBOBOX,OnBakTerms)
        topSizer.Add(bakTerms,0,WACV)
        topSizer.Add((5,0),0)
        backSizer.Add(topSizer)
        backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background coefficients:'),0,WACV)
        bakSizer = wx.FlexGridSizer(0,5,5,5)
        for i,value in enumerate(data[0][3:]):
            bakVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data[0],i+3,nDig=(10,4),typeHint=float)
            bakSizer.Add(bakVal,0,WACV)
        backSizer.Add(bakSizer)
        return backSizer
        
    def DebyeSizer():
        
        def OnDebTerms(event):
            data[1]['nDebye'] = int(debTerms.GetValue())
            M = len(data[1]['debyeTerms'])
            N = data[1]['nDebye']
            if N > M:       #add terms
                for i in range(M,N): 
                    data[1]['debyeTerms'].append([1.0,False,1.0,False,0.010,False])
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(data[1]['debyeTerms'][-1])
            #wx.CallAfter(UpdateBackground,G2frame,data)
            wx.CallLater(100,UpdateBackground,G2frame,data)

        def KeyEditPeakGrid(event):
            colList = debyeGrid.GetSelectedCols()
            if event.GetKeyCode() == wx.WXK_RETURN:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_CONTROL:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_SHIFT:
                event.Skip(True)
            elif colList:
                debyeGrid.ClearSelection()
                key = event.GetKeyCode()
                for col in colList:
                    if debyeTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                        if key == 89: #'Y'
                            for row in range(debyeGrid.GetNumberRows()): data[1]['debyeTerms'][row][col]=True
                        elif key == 78:  #'N'
                            for row in range(debyeGrid.GetNumberRows()): data[1]['debyeTerms'][row][col]=False

        
        debSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Debye scattering: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. coeff.: '),0,WACV)
        debTerms = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[1]['nDebye']),choices=[str(i) for i in range(12)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        debTerms.Bind(wx.EVT_COMBOBOX,OnDebTerms)
        topSizer.Add(debTerms,0,WACV)
        topSizer.Add((5,0),0)
        debSizer.Add(topSizer)
        if data[1]['nDebye']:
            debSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Debye diffuse terms:'),0,WACV)       
            rowLabels = []
            for i in range(len(data[1]['debyeTerms'])): rowLabels.append(str(i))
            colLabels = ['A','refine','R','refine','U','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            debyeTable = G2G.Table(data[1]['debyeTerms'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            debyeGrid = G2G.GSGrid(parent=G2frame.dataDisplay)
            debyeGrid.SetTable(debyeTable, True)
            debyeGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            debyeGrid.AutoSizeColumns(False)
            debSizer.Add(debyeGrid)        
        return debSizer
      
    def PeaksSizer():

        def OnPeaks(event):
            data[1]['nPeaks'] = int(peaks.GetValue())
            M = len(data[1]['peaksList'])
            N = data[1]['nPeaks']
            if N > M:       #add terms
                for i in range(M,N): 
                    data[1]['peaksList'].append([1.0,False,1.0,False,0.10,False,0.10,False])
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(data[1]['peaksList'][-1])
            #wx.CallAfter(UpdateBackground,G2frame,data)
            wx.CallLater(100,UpdateBackground,G2frame,data)
            
        def KeyEditPeakGrid(event):
            colList = peaksGrid.GetSelectedCols()
            if event.GetKeyCode() == wx.WXK_RETURN:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_CONTROL:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_SHIFT:
                event.Skip(True)
            elif colList:
                peaksGrid.ClearSelection()
                key = event.GetKeyCode()
                for col in colList:
                    if peaksTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                        if key == 89: #'Y'
                            for row in range(peaksGrid.GetNumberRows()): data[1]['peaksList'][row][col]=True
                        elif key == 78:  #'N'
                            for row in range(peaksGrid.GetNumberRows()): data[1]['peaksList'][row][col]=False

        peaksSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Peaks in background: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' No. peaks: '),0,WACV)
        peaks = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[1]['nPeaks']),choices=[str(i) for i in range(30)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        peaks.Bind(wx.EVT_COMBOBOX,OnPeaks)
        topSizer.Add(peaks,0,WACV)
        topSizer.Add((5,0),0)
        peaksSizer.Add(topSizer)
        G2frame.dataFrame.currentGrids = []
        if data[1]['nPeaks']:
            peaksSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Peak list:'),0,WACV)       
            rowLabels = []
            for i in range(len(data[1]['peaksList'])): rowLabels.append(str(i))
            colLabels = ['pos','refine','int','refine','sig','refine','gam','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            peaksTable = G2G.Table(data[1]['peaksList'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            peaksGrid = G2G.GSGrid(parent=G2frame.dataDisplay)
            peaksGrid.SetTable(peaksTable, True)
            peaksGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            peaksGrid.AutoSizeColumns(False)
            peaksSizer.Add(peaksGrid)        
        return peaksSizer
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.BackMenu)
    G2frame.dataFrame.SetLabel('Background')
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnBackCopy,id=G2gd.wxID_BACKCOPY)
    G2frame.Bind(wx.EVT_MENU,OnBackFlagCopy,id=G2gd.wxID_BACKFLAGCOPY)
    G2frame.Bind(wx.EVT_MENU,OnPeaksMove,id=G2gd.wxID_PEAKSMOVE)
    G2frame.Bind(wx.EVT_MENU,OnMakeRDF,id=G2gd.wxID_MAKEBACKRDF)
    G2frame.Bind(wx.EVT_MENU,OnBkgFit,id=G2frame.dataFrame.wxID_BackPts['Fit'])
    G2frame.Bind(wx.EVT_MENU,OnBkgClear,id=G2frame.dataFrame.wxID_BackPts['Clear'])    
    BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Background')
    Choices = ['chebyschev','cosine','Q^2 power series','Q^-2 power series','lin interpolate','inv interpolate','log interpolate']
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(BackSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(DebyeSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(PeaksSizer())
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(Size)
        
################################################################################
#####  Limits
################################################################################           
       
def UpdateLimitsGrid(G2frame, data,plottype):
    '''respond to selection of PWDR Limits data tree item.
    '''
    def AfterChange(invalid,value,tc):
        if invalid: return
        plottype = G2frame.PatternTree.GetItemText(G2frame.PatternId)[:4]
#        print 'new plot'
        wx.CallAfter(G2plt.PlotPatterns,G2frame,newPlot=False,plotType=plottype)  #unfortunately this resets the plot width

    def LimitSizer():
        limits = wx.FlexGridSizer(0,3,0,5)
        labels = ['Tmin','Tmax']
        for i in [0,1]:
            limits.Add(wx.StaticText(G2frame.dataDisplay,label=' Original %s: %.4f'%(labels[0],data[0][i])),0,WACV)
            limits.Add(wx.StaticText(G2frame.dataDisplay,label=' New: '),0,WACV)
            limits.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data[1],i,  \
                min=data[0][0],max=data[0][1],nDig=(10,4),typeHint=float,OnLeave=AfterChange))
        return limits
        
    def ExclSizer():
        
        def OnDelExcl(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            del(data[item+2])
            G2plt.PlotPatterns(G2frame,newPlot=False,plotType=plottype)
            wx.CallAfter(UpdateLimitsGrid,G2frame,data,plottype)
        
        Indx = {}
        excl = wx.FlexGridSizer(0,3,0,5)
        excl.Add(wx.StaticText(G2frame.dataDisplay,label=' From: '),0,WACV)
        excl.Add(wx.StaticText(G2frame.dataDisplay,label=' To: '),0,WACV)
        excl.Add(wx.StaticText(G2frame.dataDisplay,label=' Delete?: '),0,WACV)
        for id,item in enumerate(data[2:]):
            for i in [0,1]:
                excl.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,item,i,  \
                    min=data[0][0],max=data[0][1],nDig=(10,4),typeHint=float,OnLeave=AfterChange))
            delExcl = wx.CheckBox(G2frame.dataDisplay,label='')
            Indx[delExcl.GetId()] = id
            delExcl.Bind(wx.EVT_CHECKBOX,OnDelExcl)
            excl.Add(delExcl,0,WACV)
        return excl
               
    def OnAddExcl(event):
        G2frame.ifGetExclude = True
        print 'Add excluded region'
        
    def OnLimitCopy(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy limits from\n'+str(hst[5:])+' to...',
            'Copy limits', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections(): 
                    item = histList[i]
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    G2frame.PatternTree.SetItemPyData(
                        G2gd.GetPatternTreeItemId(G2frame,Id,'Limits'),copy.copy(data))
        finally:
            dlg.Destroy()
            
    def Draw():
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(LimitSizer())
        if len(data)>2:
#            mainSizer.Add((0,5),0)
            mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Excluded regions:'),0,WACV)
            mainSizer.Add(ExclSizer())
        mainSizer.Layout()    
        G2frame.dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()
    G2frame.ifGetExclude = False
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.LimitMenu)
    G2frame.dataFrame.SetLabel('Limits')
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnLimitCopy,id=G2gd.wxID_LIMITCOPY)
    G2frame.Bind(wx.EVT_MENU,OnAddExcl,id=G2gd.wxID_ADDEXCLREGION)
    Draw() 
    
################################################################################
#####  Instrument parameters
################################################################################           
       
def UpdateInstrumentGrid(G2frame,data):
    '''respond to selection of PWDR/SASD/REFD Instrument Parameters
    data tree item.
    '''
    if 'Bank' not in data:  #get it from name; absent for default parms selection 
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        if 'Bank' in hst:
            bank = int(hst.split('Bank')[1].split('_')[0])
            data['Bank'] = [bank,bank,0]
        else:
            data['Bank'] = [1,1,0]

    def keycheck(keys):
        good = []
        for key in keys:
            if key in ['Type','Bank','U','V','W','X','Y','SH/L','I(L2)/I(L1)','alpha',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','Polariz.',
                'Lam','Azimuth','2-theta','fltPath','difC','difA','difB','Zero','Lam1','Lam2']:
                good.append(key)
        return good
        
    def updateData(inst,ref):
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[0]
        for item in data:
            try:
                data[item] = [data[item][0],inst[item],ref[item]]
            except KeyError:
                try:
                    data[item] = [data[item][0],inst[item]]
                except KeyError:
                    pass        #skip 'Polariz.' for N-data
    
    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
            newpeaks = []
            for peak in peaks['peaks']:
                newpeaks.append(G2mth.setPeakparms(data,Inst2,peak[0],peak[2]))
            peaks['peaks'] = newpeaks
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'),peaks)
            
    def OnCalibrate(event):
        Pattern = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)
        xye = ma.array(ma.getdata(Pattern[1]))
        cw = np.diff(xye[0])
        IndexPeaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
        if not len(IndexPeaks[0]):
            G2frame.ErrorDialog('Can not calibrate','Index Peak List empty')
            return
        if not np.any(IndexPeaks[1]):
            G2frame.ErrorDialog('Can not calibrate','Peak positions not refined')
            return False
        Ok = False
        for peak in IndexPeaks[0]:
            if peak[2] and peak[3]:
                Ok = True
        if not Ok:
            G2frame.ErrorDialog('Can not calibrate','Index Peak List not indexed')
            return            
        if G2pwd.DoCalibInst(IndexPeaks,data):
            UpdateInstrumentGrid(G2frame,data)
            XY = []
            Sigs = []
            for ip,peak in enumerate(IndexPeaks[0]):
                if peak[2] and peak[3]:
                    binwid = cw[np.searchsorted(xye[0],peak[0])]
                    XY.append([peak[-1],peak[0],binwid])
                    Sigs.append(IndexPeaks[1][ip])
            if len(XY):
                XY = np.array(XY)
                G2plt.PlotCalib(G2frame,data,XY,Sigs,newPlot=True)
        else:
            G2frame.ErrorDialog('Can not calibrate','Nothing selected for refinement')

    def OnLoad(event):
        '''Loads instrument parameters from a G2 .instprm file
        in response to the Instrument Parameters-Operations/Load Profile menu
        If instprm file has multiple banks each with header #Bank n: ..., this 
        finds matching bank no. to load - rejects nonmatches.
        
        Note that similar code is found in ReadPowderInstprm (GSASII.py)
        '''
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[0]
        bank = data['Bank'][0]
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', pth, '', 
            'instrument parameter files (*.instprm)|*.instprm',wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                newItems = []
                newVals = []
                Found = False
                while S:
                    if S[0] == '#':
                        if Found:
                            break
                        if 'Bank' in S:
                            if bank == int(S.split(':')[0].split()[1]):
                                S = File.readline()
                                continue
                            else:
                                S = File.readline()
                                while S and '#Bank' not in S:
                                    S = File.readline()
                                continue
                        else:   #a non #Bank file
                            S = File.readline()
                            continue
                    Found = True
                    [item,val] = S[:-1].split(':')
                    newItems.append(item)
                    try:
                        newVals.append(float(val))
                    except ValueError:
                        newVals.append(val)                        
                    S = File.readline()                
                File.close()
                if Found:
                    Inst,Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Instrument Parameters'))
                    if 'Bank' not in Inst:  #patch for old .instprm files - may cause faults for TOF data
                        Inst['Bank'] = [1,1,0]
                    data = G2IO.makeInstDict(newItems,newVals,len(newVals)*[False,])
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Instrument Parameters'),[data,Inst2])
                    RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
                else:
                    G2frame.ErrorDialog('No match','Bank %d not in %s'%(bank,filename),G2frame.dataFrame)
                UpdateInstrumentGrid(G2frame,data)
                G2plt.PlotPeakWidths(G2frame)
        finally:
            dlg.Destroy()
        
    def OnSave(event):
        '''Respond to the Instrument Parameters Operations/Save Profile menu
        item: writes current parameters to a .instprm file
        It does not write Bank n: on # line & thus can be used any time w/o clash of bank nos.
        '''
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', pth, '', 
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .instprm
                filename = os.path.splitext(filename)[0]+'.instprm'
                File = open(filename,'w')
                File.write("#GSAS-II instrument parameter file; do not add/delete items!\n")
                for item in data:
                    File.write(item+':'+str(data[item][1])+'\n')
                File.close()
                print 'Instrument parameters saved to: '+filename
        finally:
            dlg.Destroy()
            
    def OnSaveAll(event):
        '''Respond to the Instrument Parameters Operations/Save all Profile menu & writes 
        selected inst parms. across multiple banks into a single file
        Each block starts with #Bank n: GSAS-II instrument... where n is bank no.
        item: writes parameters from selected PWDR entries to a .instprm file
        '''
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        histList.insert(0,hst)
        saveList = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Save instrument parameters from',
            'Save instrument parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    saveList.append(histList[i])
        finally:
            dlg.Destroy()
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', pth, '', 
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .instprm
                filename = os.path.splitext(filename)[0]+'.instprm'
                File = open(filename,'w')
                for hist in saveList:
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,hist)
                    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))[0]
                    if 'Bank' not in inst:  #patch
                        bank = 1
                        if 'Bank' in hist:
                            bank = int(hist.split('Bank')[1])
                        inst['Bank'] = [bank,bank,0]
                    bank = inst['Bank'][0]                
                    File.write("#Bank %d: GSAS-II instrument parameter file; do not add/delete items!\n"%(bank))
                    for item in inst:
                        File.write(item+':'+str(inst[item][1])+'\n')                                    
                File.close()
        finally:
            dlg.Destroy()
                                                
    def OnReset(event):
        insVal.update(insDef)
        updateData(insVal,insRef)
        RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
        UpdateInstrumentGrid(G2frame,data)
        G2plt.PlotPeakWidths(G2frame)
        
    def OnInstFlagCopy(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        keys = data.keys()
        try:
            keys.remove('Source')
        except ValueError:
            pass
        flags = dict(zip(keys,[data[key][2] for key in keys]))
        instType = data['Type'][0]
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy inst ref. flags from\n'+hst[5:],
            'Copy refinement flags', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            instData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))[0]
            if 'Bank' not in instData:
                instData['Bank'] = [1,1,0]
#            GSASIIpath.IPyBreak()
            if len(data) == len(instData) and instType == instData['Type'][0]:   #don't mix data types or lam & lam1/lam2 parms!
                for item in instData:
                    if item not in ['Source',]:
                        instData[item][2] = copy.copy(flags[item])
            else:
                print item+' not copied - instrument parameters not commensurate'
        
    def OnInstCopy(event):
        #need fix for dictionary
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        instType = data['Type'][0]
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy inst params from\n'+hst,
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections(): 
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            instData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))[0]
            if 'Bank' not in instData:
                instData['Bank'] = [1,1,0]
            if len(data) == len(instData) and instType == instData['Type'][0]:  #don't mix data types or lam & lam1/lam2 parms!
                instData.update(data)
            else:
                print item+' not copied - instrument parameters not commensurate'
                         
    def AfterChange(invalid,value,tc):
        if invalid: return
        updateData(insVal,insRef)
        
    def OnItemRef(event):
        Obj = event.GetEventObject()
        item = RefObj[Obj.GetId()]
        insRef[item] = Obj.GetValue()
        updateData(insVal,insRef)

    def OnCopy1Val(event):
        '''Select one instrument parameter value to edit and copy to many histograms
        optionally allow values to be edited in a table
        '''
        updateData(insVal,insRef)
        G2G.SelectEdit1Var(G2frame,data,labelLst,elemKeysLst,dspLst,refFlgElem)
        insVal.update({key:data[key][1] for key in instkeys})
        insRef.update({key:data[key][2] for key in instkeys})
        wx.CallAfter(MakeParameterWindow)
        
    def lblWdef(lbl,dec,val):
        'Label parameter showing the default value'
        fmt = "%15."+str(dec)+"f"
        return " " + lbl + " (" + (fmt % val).strip() + "): "

    def RefineBox(item):
        'Define a refine checkbox with binding'
        wid = wx.CheckBox(G2frame.dataDisplay,label=' Refine?  ')
        wid.SetValue(bool(insRef[item]))
        RefObj[wid.GetId()] = item
        wid.Bind(wx.EVT_CHECKBOX, OnItemRef)
        return wid

    def OnLamPick(event):
        data['Source'][1] = lamType = event.GetEventObject().GetValue()
        if 'P' in insVal['Type']:
            insVal['Lam1'] = waves[lamType][0]
            insVal['Lam2'] = waves[lamType][1]
        elif 'S' in insVal['Type'] and 'synch' not in lamType:
            insVal['Lam'] = meanwaves[lamType]
        updateData(insVal,insRef)
        i,j= wx.__version__.split('.')[0:2]
        if int(i)+int(j)/10. > 2.8:
            pass # repaint crashes wxpython 2.9
            wx.CallLater(100, MakeParameterWindow)
            #wx.CallAfter(MakeParameterWindow)
        else:
            wx.CallAfter(MakeParameterWindow)

    def MakeParameterWindow():
        'Displays the Instrument parameters in the datadisplay frame'
        if G2frame.dataDisplay:
            G2frame.dataFrame.Clear()
        G2frame.dataFrame.SetLabel('Instrument Parameters')
        G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        instSizer = wx.FlexGridSizer(0,6,5,5)
        subSizer = wx.BoxSizer(wx.HORIZONTAL)
        if insVal['Bank'] == None:      #patch
            insVal['Bank'] = 1
        text = ' Histogram Type: %s  Bank: %d'%(insVal['Type'],insVal['Bank'])
        subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,text),0,WACV)
        mainSizer.Add(subSizer)
        labelLst[:],elemKeysLst[:],dspLst[:],refFlgElem[:] = [],[],[],[]
        if 'P' in insVal['Type']:                   #powder data
            if 'C' in insVal['Type']:               #constant wavelength
                labelLst.append('Azimuth angle')
                elemKeysLst.append(['Azimuth',1])
                dspLst.append([10,2])
                refFlgElem.append(None)                   
                if 'Lam1' in insVal:
                    subSizer = wx.BoxSizer(wx.HORIZONTAL)
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Azimuth: '),0,WACV)
                    txt = '%7.2f'%(insVal['Azimuth'])
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'   Ka1/Ka2: '),0,WACV)
                    txt = u'  %8.6f/%8.6f\xc5'%(insVal['Lam1'],insVal['Lam2'])
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                    waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                    waveSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Source type: '),0,WACV)
                    # PATCH?: for now at least, Source is not saved anywhere before here
                    if 'Source' not in data: data['Source'] = ['CuKa','?']
                    choice = ['TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa']
                    lamPick = wx.ComboBox(G2frame.dataDisplay,value=data['Source'][1],choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                    waveSizer.Add(lamPick,0)
                    subSizer.Add(waveSizer,0)
                    mainSizer.Add(subSizer)
                    instSizer.Add(wx.StaticText(
                        G2frame.dataDisplay,-1,
                        lblWdef('I(L2)/I(L1)',4,insDef['I(L2)/I(L1)'])),
                        0,WACV)
                    key = 'I(L2)/I(L1)'
                    labelLst.append(key)
                    elemKeysLst.append([key,1])
                    dspLst.append([10,4])
                    refFlgElem.append([key,2])                   
                    ratVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,key,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                    instSizer.Add(ratVal,0)
                    instSizer.Add(RefineBox(key),0,WACV)
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)                
                else: # single wavelength
                    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Azimuth: '),0,WACV)
                    txt = '%7.2f'%(insVal['Azimuth'])
                    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                    instSizer.Add((5,5),0)
                    key = 'Lam'
                    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef[key])),
                        0,WACV)
                    waveVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,key,nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                    labelLst.append(u'Lam (\xc5)')
                    elemKeysLst.append([key,1])
                    dspLst.append([10,6])
                    instSizer.Add(waveVal,0,WACV)
                    refFlgElem.append([key,2])                   
                    instSizer.Add(RefineBox(key),0,WACV)
#                    if ifHisto:
#                        refFlgElem.append([key,2])                   
#                        instSizer.Add(RefineBox(key),0,WACV)
#                    else:
#                        refFlgElem.append(None)                   
#                        instSizer.Add((5,5),0)
                for item in ['Zero','Polariz.']:
                    if item in insDef:
                        labelLst.append(item)
                        elemKeysLst.append([item,1])
                        dspLst.append([10,4])
                        instSizer.Add(
                            wx.StaticText(G2frame.dataDisplay,-1,lblWdef(item,4,insDef[item])),
                            0,WACV)
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,item,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                        instSizer.Add(itemVal,0,WACV)
                        refFlgElem.append([item,2])
                        instSizer.Add(RefineBox(item),0,WACV)
#                        if ifHisto:
#                            refFlgElem.append([item,2])
#                            instSizer.Add(RefineBox(item),0,WACV)
#                        else:
#                            refFlgElem.append(None)                   
#                            instSizer.Add((5,5),0)
                    else:                           #skip Polariz. for neutrons
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                for item in ['U','V','W','','X','Y','SH/L']:
                    if item == '':
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                        continue
                    nDig = (10,3)
                    if item == 'SH/L':
                        nDig = (10,5)
                    labelLst.append(item)
                    elemKeysLst.append([item,1])
                    dspLst.append(nDig)
                    refFlgElem.append([item,2])
                    instSizer.Add(
                        wx.StaticText(G2frame.dataDisplay,-1,lblWdef(item,nDig[1],insDef[item])),
                        0,WACV)
                    itemVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,item,nDig=nDig,typeHint=float,OnLeave=AfterChange)
                    instSizer.Add(itemVal,0,WACV)
                    instSizer.Add(RefineBox(item),0,WACV)
            elif 'T' in insVal['Type']:                                   #time of flight (neutrons)
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Fligth path: '),0,WACV)
                txt = '%8.3f'%(insVal['fltPath'])
                subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                labelLst.append('flight path')
                elemKeysLst.append(['fltPath',1])
                dspLst.append([10,2])
                refFlgElem.append(None)                   
                subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  2-theta: '),0,WACV)
                txt = '%7.2f'%(insVal['2-theta'])
                subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                labelLst.append('2-theta')
                elemKeysLst.append(['2-theta',1])
                dspLst.append([10,2])
                refFlgElem.append(None)                   
                if 'Pdabc' in Inst2:
                    Items = ['sig-0','sig-1','sig-2','sig-q','X','Y']
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  difC: '),0,WACV)
                    txt = '%8.2f'%(insVal['difC'])
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,txt.strip()),0,WACV)
                    labelLst.append('difC')
                    elemKeysLst.append(['difC',1])
                    dspLst.append([10,2])
                    refFlgElem.append(None)
                    subSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  alpha, beta: fixed by table'),0,WACV)
                else:
                    Items = ['difC','difA','difB','Zero','alpha','beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','X','Y']
                mainSizer.Add((5,5),0)
                mainSizer.Add(subSizer)
                mainSizer.Add((5,5),0)
                for item in Items:
                    if item == '':
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                        instSizer.Add((5,5),0)
                        continue
                    nDig = (10,3)
                    if 'beta' in item:
                        nDig = (12,6)
                    instSizer.Add(
                            wx.StaticText(G2frame.dataDisplay,-1,lblWdef(item,nDig[1],insDef[item])),
                            0,WACV)
                    itemVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,item,nDig=nDig,typeHint=float,OnLeave=AfterChange)
                    instSizer.Add(itemVal,0,WACV)
                    labelLst.append(item)
                    elemKeysLst.append([item,1])
                    dspLst.append(nDig)
                    refFlgElem.append([item,2])
                    instSizer.Add(RefineBox(item),0,WACV)
            elif 'PKS' in insVal['Type']:   #peak positions only
                key = 'Lam'
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef[key])),
                    0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,key,nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                labelLst.append(u'Lam (\xc5)')
                elemKeysLst.append([key,1])
                dspLst.append([10,6])
                instSizer.Add(waveVal,0,WACV)
                refFlgElem.append([key,2])                   
#                    instSizer.Add(RefineBox(key),0,WACV)
                for item in ['Zero',]:
                    if item in insDef:
                        labelLst.append(item)
                        elemKeysLst.append([item,1])
                        dspLst.append([10,4])
                        instSizer.Add(
                            wx.StaticText(G2frame.dataDisplay,-1,lblWdef(item,4,insDef[item])),
                            0,WACV)
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,item,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                        instSizer.Add(itemVal,0,WACV)
                        refFlgElem.append([item,2])
#                        instSizer.Add(RefineBox(item),0,WACV)
                
                
        elif 'S' in insVal['Type']:                       #single crystal data
            if 'C' in insVal['Type']:               #constant wavelength
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),
                    0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,'Lam',nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                instSizer.Add(waveVal,0,WACV)
                labelLst.append(u'Lam (\xc5)')
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Source type: '),0,WACV)
                # PATCH?: for now at least, Source is not saved anywhere before here
                if 'Source' not in data: data['Source'] = ['CuKa','?']
                choice = ['synchrotron','TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa']
                lamPick = wx.ComboBox(G2frame.dataDisplay,value=data['Source'][1],choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0,WACV)
                instSizer.Add(waveSizer,0,WACV)
                elemKeysLst.append(['Lam',1])
                dspLst.append([10,6])
                refFlgElem.append(None)
            else:                                   #time of flight (neutrons)
                pass                                #for now
        elif insVal['Type'][0] in ['L','R',]:
            if 'C' in insVal['Type']:        
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),
                    0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,insVal,'Lam',nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                instSizer.Add(waveVal,0,WACV)
                labelLst.append(u'Lam (\xc5)')
                elemKeysLst.append(['Lam',1])
                dspLst.append([10,6])
                refFlgElem.append(None)
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Azimuth: %7.2f'%(insVal['Azimuth'])),0,WACV)
                labelLst.append('Azimuth angle')
                elemKeysLst.append(['Azimuth',1])
                dspLst.append([10,2])
                refFlgElem.append(None)                   
            else:                                   #time of flight (neutrons)
                pass                                #for now

        mainSizer.Add(instSizer,0)
        mainSizer.Layout()    
        G2frame.dataDisplay.SetSizer(mainSizer)
        G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
        G2frame.dataFrame.SendSizeEvent()  # this causes a frame repaint, even if the size does not change!
        # end of MakeParameterWindow
                
    # beginning of UpdateInstrumentGrid code    
    #patch: make sure all parameter items are lists
    patched = 0
    for key in data:
        if type(data[key]) is tuple:
            data[key] = list(data[key])
            patched += 1
    if patched: print patched,' instrument parameters changed from tuples'
    #end of patch
    labelLst,elemKeysLst,dspLst,refFlgElem = [],[],[],[]
    instkeys = keycheck(data.keys())
    if 'P' in data['Type'][0]:          #powder data
        insVal = dict(zip(instkeys,[data[key][1] for key in instkeys]))
        insDef = dict(zip(instkeys,[data[key][0] for key in instkeys]))
        insRef = dict(zip(instkeys,[data[key][2] for key in instkeys]))
        if 'NC' in data['Type'][0]:
            del(insDef['Polariz.'])
            del(insVal['Polariz.'])
            del(insRef['Polariz.'])
    elif 'S' in data['Type'][0]:                               #single crystal data
        insVal = dict(zip(instkeys,[data[key][1] for key in instkeys]))
        insDef = dict(zip(instkeys,[data[key][0] for key in instkeys]))
        insRef = {}
    elif 'L' in data['Type'][0]:                               #low angle data
        insVal = dict(zip(instkeys,[data[key][1] for key in instkeys]))
        insDef = dict(zip(instkeys,[data[key][0] for key in instkeys]))
        insRef = {}
    elif 'R' in data['Type'][0]:                               #low angle data
        insVal = dict(zip(instkeys,[data[key][1] for key in instkeys]))
        insDef = dict(zip(instkeys,[data[key][0] for key in instkeys]))
        insRef = {}
    RefObj = {}
    waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
        'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'MoKa':[0.70926,0.713543],
        'AgKa':[0.559363,0.563775]}
    meanwaves = {'CuKa':1.5418,'TiKa':2.7496,'CrKa':2.2909,'FeKa':1.9373,
        'CoKa':1.7902,'MoKa':0.7107,'AgKa':0.5608}
    Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[1]        
    G2gd.SetDataMenuBar(G2frame)
    #patch
    if 'P' in insVal['Type']:                   #powder data
        if 'C' in insVal['Type']:               #constant wavelength
            if 'Azimuth' not in insVal:
                insVal['Azimuth'] = 0.0
                insDef['Azimuth'] = 0.0
                insRef['Azimuth'] = False
#        if 'T' in insVal['Type']:
#            if 'difB' not in insVal:
#                insVal['difB'] = 0.0
#                insDef['difB'] = 0.0
#                insRef['difB'] = False
    #end of patch
    if 'P' in insVal['Type']:                   #powder data menu commands
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.InstMenu)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()            
            Status.SetStatusText('NB: Azimuth is used for polarization only')
        G2frame.Bind(wx.EVT_MENU,OnCalibrate,id=G2gd.wxID_INSTCALIB)
        G2frame.Bind(wx.EVT_MENU,OnLoad,id=G2gd.wxID_INSTLOAD)
        G2frame.Bind(wx.EVT_MENU,OnSave,id=G2gd.wxID_INSTSAVE)
        G2frame.Bind(wx.EVT_MENU,OnSaveAll,id=G2gd.wxID_INSTSAVEALL)
        G2frame.Bind(wx.EVT_MENU,OnReset,id=G2gd.wxID_INSTPRMRESET)
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2gd.wxID_INSTCOPY)
        G2frame.Bind(wx.EVT_MENU,OnInstFlagCopy,id=G2gd.wxID_INSTFLAGCOPY)
        #G2frame.Bind(wx.EVT_MENU,OnWaveChange,id=G2gd.wxID_CHANGEWAVETYPE)        
        G2frame.Bind(wx.EVT_MENU,OnCopy1Val,id=G2gd.wxID_INST1VAL)
    elif 'L' in insVal['Type'] or 'R' in insVal['Type']:                   #SASD data menu commands
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.SASDInstMenu)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2gd.wxID_INSTCOPY)
    MakeParameterWindow()
        
    
################################################################################
#####  Sample parameters
################################################################################           
       
def UpdateSampleGrid(G2frame,data):
    '''respond to selection of PWDR/SASD Sample Parameters
    data tree item.
    '''

    def OnSampleSave(event):
        '''Respond to the Sample Parameters Operations/Save menu
        item: writes current parameters to a .samprm file
        '''
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II sample parameters file', pth, '', 
            'sample parameter files (*.samprm)|*.samprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .samprm
                filename = os.path.splitext(filename)[0]+'.samprm'
                File = open(filename,'w')
                File.write("#GSAS-II sample parameter file\n")
                File.write("'Type':'"+str(data['Type'])+"'\n")
                File.write("'Gonio. radius':"+str(data['Gonio. radius'])+"\n")
                if data.get('InstrName'):
                    File.write("'InstrName':'"+str(data['InstrName'])+"'\n")
                File.close()
        finally:
            dlg.Destroy()
            
    def OnSampleLoad(event):
        '''Loads sample parameters from a G2 .samprm file
        in response to the Sample Parameters-Operations/Load menu
        
        Note that similar code is found in ReadPowderInstprm (GSASII.py)
        '''
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II sample parameters file', pth, '', 
            'sample parameter files (*.samprm)|*.samprm',wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                newItems = {}
                while S:
                    if S[0] == '#':
                        S = File.readline()
                        continue
                    [item,val] = S[:-1].split(':')
                    newItems[item.strip("'")] = eval(val)
                    S = File.readline()                
                File.close()
                data.update(newItems)
                G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Sample Parameters'),data)
                UpdateSampleGrid(G2frame,data)
        finally:
            dlg.Destroy()
            
    def OnAllSampleLoad(event):
        filename = ''
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose multihistogram metadata text file', pth, '', 
            'metadata file (*.*)|*.*',wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                newItems = []
                itemNames = []
                Comments = []
                while S:
                    if S[0] == '#':
                        Comments.append(S)
                        S = File.readline()
                        continue
                    S = S.replace(',',' ').replace('\t',' ')
                    Stuff = S[:-1].split()
                    itemNames.append(Stuff[0])
                    newItems.append(Stuff[1:])
                    S = File.readline()                
                File.close()
        finally:
            dlg.Destroy()
        if not filename:
            G2frame.ErrorDialog('Nothing to do','No file selected')
            return
        dataDict = dict(zip(itemNames,newItems))
        ifany = False
        Controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Controls'))
        Names = [' ','Phi','Chi','Omega','Time','Temperature','Pressure']
        freeNames = {}
        for name in ['FreePrm1','FreePrm2','FreePrm3']:
            freeNames[Controls[name]] = name
            Names.append(Controls[name])
        dlg = G2G.G2ColumnIDDialog( G2frame,' Choose multihistogram metadata columns:',
            'Select columns',Comments,Names,np.array(newItems).T)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                colNames,newData = dlg.GetSelection()
                dataDict = dict(zip(itemNames,newData.T))
                for item in colNames:
                    if item != ' ':
                        ifany = True
        finally:
            dlg.Destroy()
        if not ifany:
            G2frame.ErrorDialog('Nothing to do','No columns identified')
            return
        histList = [G2frame.PatternTree.GetItemText(G2frame.PatternId),]
        histList += GetHistsLikeSelected(G2frame)
        colIds = {}
        for i,name in enumerate(colNames):
            if name != ' ':
                colIds[name] = i
        for hist in histList:
            name = hist.split()[1]  #this is file name
            newItems = {}
            for item in colIds:
                key = freeNames.get(item,item)
                newItems[key] = float(dataDict[name][colIds[item]])
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,hist)
            sampleData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
            sampleData.update(newItems)        
        UpdateSampleGrid(G2frame,data)        
    
    def OnSetScale(event):
        if histName[:4] in ['REFD','PWDR']:
            Scale = data['Scale'][0]
            dlg = wx.MessageDialog(G2frame,'Rescale data by %.2f?'%(Scale),'Rescale data',wx.OK|wx.CANCEL)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,histName)
                    y,w = G2frame.PatternTree.GetItemPyData(pId)[1][1:3]
                    y *= Scale
                    w /= Scale**2
                    data['Scale'][0] = 1.0
            finally:
                dlg.Destroy()
            G2plt.PlotPatterns(G2frame,plotType=histName[:4],newPlot=True)
            UpdateSampleGrid(G2frame,data)
            return
        #SASD rescaliing                
        histList = []
        item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while item:
            name = G2frame.PatternTree.GetItemText(item)
            if 'SASD' in name and name != histName:
                histList.append(name)
            item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        if not len(histList):      #nothing to copy to!
            return
        dlg = wx.SingleChoiceDialog(G2frame,'Select reference histogram for scaling',
            'Reference histogram',histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                refHist = histList[sel]
        finally:
            dlg.Destroy()
        Limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits'))
        Profile = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[1]
        Data = [Profile,Limits,data]
        refId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,refHist)
        refSample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,refId, 'Sample Parameters'))
        refLimits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,refId, 'Limits'))
        refProfile = G2frame.PatternTree.GetItemPyData(refId)[1]
        refData = [refProfile,refLimits,refSample]
        G2sasd.SetScale(Data,refData)
        G2plt.PlotPatterns(G2frame,plotType='SASD',newPlot=True)
        UpdateSampleGrid(G2frame,data)
        
    def OnRescaleAll(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        x0,y0,w0 = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[1][:3]
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        od = {'label_1':'Scaling range min','value_1':0.0,'label_2':'Scaling range max','value_2':10.}
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame,
            'Do scaling from\n'+str(hst[5:])+' to...','Rescale histograms', histList,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Xmin = od['value_1']
                Xmax = od['value_2']
                iBeg = np.searchsorted(x0,Xmin)
                iFin = np.searchsorted(x0,Xmax)
                if iBeg > iFin:
                    wx.MessageBox('Wrong order for Xmin, Xmax','Error',style=wx.ICON_EXCLAMATION)
                else:
                    sum0 = np.sum(y0[iBeg:iFin])
                    result = dlg.GetSelections()
                    for i in result: 
                        item = histList[i]
                        Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                        xi,yi,wi = G2frame.PatternTree.GetItemPyData(Id)[1][:3]
                        sumi = np.sum(yi[iBeg:iFin])
                        if sumi:
                            Scale = sum0/sumi
                            yi *= Scale
                            wi /= Scale**2
        finally:
            dlg.Destroy()
        G2plt.PlotPatterns(G2frame,plotType=histName[:4],newPlot=True)
        
    def OnSampleCopy(event):
        histType,copyNames = SetCopyNames(histName,data['Type'],
            addNames = ['Omega','Chi','Phi','Gonio. radius','InstrName'])
        copyDict = {}
        for parm in copyNames:
            copyDict[parm] = data[parm]
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame,
            'Copy sample params from\n'+str(hst[5:])+' to...',
            'Copy sample parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    item = histList[i]
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
                    sampleData.update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnSampleCopySelected(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        Controls = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        # Assemble a list of item labels
        TextTable = {key:label for key,label,dig in
                     SetupSampleLabels(hst,data.get('Type'),Inst['Type'][0])
                     }
        # get flexible labels
        TextTable.update({
            key:Controls[key] for key in Controls if key.startswith('FreePrm')
            })
        # add a few extra
        TextTable.update({
            'Type':'Diffractometer type',
            'InstrName':'Instrument Name',
            })
        # Assemble a list of dict entries that would be labeled in the Sample
        # params data window (drop ranId and items not used).
        keyList = [i for i in data.keys() if i in TextTable]
        keyText = [TextTable[i] for i in keyList]
        # sort both lists together, ordered by keyText
        keyText, keyList = zip(*sorted(zip(keyText,keyList))) # sort lists 
        selectedKeys = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame,
            'Select which sample parameters\nto copy',
            'Select sample parameters', keyText)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                selectedKeys = [keyList[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if not selectedKeys: return # nothing to copy
        copyDict = {}
        for parm in selectedKeys:
            copyDict[parm] = data[parm]
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame,
            'Copy sample params from\n'+str(hst[5:])+' to...',
            'Copy sample parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    item = histList[i]
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
                    sampleData.update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()            
        G2plt.PlotPatterns(G2frame,plotType=hst[:4],newPlot=False)

    def OnSampleFlagCopy(event):
        histType,copyNames = SetCopyNames(histName,data['Type'])
        flagDict = {}
        for parm in copyNames:
            flagDict[parm] = data[parm][1]
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy sample ref. flags from\n'+str(hst[5:])+' to...',
            'Copy sample flags', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    item = histList[i]
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
                    for name in copyNames:
                        sampleData[name][1] = copy.copy(flagDict[name])
        finally:
            dlg.Destroy()

    def OnHistoChange():
        '''Called when the histogram type is changed to refresh the window
        '''
        #wx.CallAfter(UpdateSampleGrid,G2frame,data)
        wx.CallLater(100,UpdateSampleGrid,G2frame,data)
        
    def SetNameVal():
        inst = instNameVal.GetValue()
        data['InstrName'] = inst.strip()

    def OnNameVal(event):
        event.Skip()
        wx.CallAfter(SetNameVal)
        
    def AfterChange(invalid,value,tc):
        if invalid:
            return
        if tc.key == 0 and 'SASD' in histName:          #a kluge for Scale!
            G2plt.PlotPatterns(G2frame,plotType='SASD',newPlot=True)
        elif tc.key == 'Thick':
            wx.CallAfter(UpdateSampleGrid,G2frame,data)            
            
    def OnMaterial(event):
        Obj = event.GetEventObject()
        id = Info[Obj.GetId()]
        data['Materials'][id]['Name'] = Obj.GetValue()
        wx.CallAfter(UpdateSampleGrid,G2frame,data)
        
    def OnVolFrac(invalid,value,tc):
        id = Info[tc.GetId()]
        data['Materials'][not id][key] = 1.-value
        wx.CallAfter(UpdateSampleGrid,G2frame,data)

    def OnCopy1Val(event):
        'Select one value to copy to many histograms and optionally allow values to be edited in a table'
        G2G.SelectEdit1Var(G2frame,data,labelLst,elemKeysLst,dspLst,refFlgElem)
        wx.CallAfter(UpdateSampleGrid,G2frame,data)
        
    ######## DEBUG #######################################################
    #import GSASIIpwdGUI
    #reload(GSASIIpwdGUI)
    #reload(G2gd)
    ######################################################################
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(
            G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    histName = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.SampleMenu)
    G2frame.dataFrame.SetLabel('Sample Parameters')
    G2frame.Bind(wx.EVT_MENU, OnSetScale, id=G2gd.wxID_SETSCALE)
    G2frame.Bind(wx.EVT_MENU, OnSampleCopy, id=G2gd.wxID_SAMPLECOPY)
    G2frame.Bind(wx.EVT_MENU, OnSampleCopySelected, id=G2gd.wxID_SAMPLECOPYSOME)
    G2frame.Bind(wx.EVT_MENU, OnSampleFlagCopy, id=G2gd.wxID_SAMPLEFLAGCOPY)
    G2frame.Bind(wx.EVT_MENU, OnSampleSave, id=G2gd.wxID_SAMPLESAVE)
    G2frame.Bind(wx.EVT_MENU, OnSampleLoad, id=G2gd.wxID_SAMPLELOAD)
    G2frame.Bind(wx.EVT_MENU, OnCopy1Val, id=G2gd.wxID_SAMPLE1VAL)
    G2frame.Bind(wx.EVT_MENU, OnAllSampleLoad, id=G2gd.wxID_ALLSAMPLELOAD)
    G2frame.Bind(wx.EVT_MENU, OnRescaleAll, id=G2gd.wxID_RESCALEALL)
    if histName[:4] in ['SASD','REFD','PWDR']:
        G2frame.dataFrame.SetScale.Enable(True)
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()    
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    Controls = G2frame.PatternTree.GetItemPyData(
        G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
#patch
    if 'ranId' not in data:
        data['ranId'] = ran.randint(0,sys.maxint)
    if not 'Gonio. radius' in data:
        data['Gonio. radius'] = 200.0
    if not 'Omega' in data:
        data.update({'Omega':0.0,'Chi':0.0,'Phi':0.0})
    if 'Azimuth' not in data:
        data['Azimuth'] = 0.0
    if type(data['Temperature']) is int:
        data['Temperature'] = float(data['Temperature'])
    if 'Time' not in data:
        data['Time'] = 0.0
    if 'FreePrm1' not in Controls:
        Controls['FreePrm1'] = 'Sample humidity (%)'
    if 'FreePrm2' not in Controls:
        Controls['FreePrm2'] = 'Sample voltage (V)'
    if 'FreePrm3' not in Controls:
        Controls['FreePrm3'] = 'Applied load (MN)'
    if 'FreePrm1' not in data:
        data['FreePrm1'] = 0.
    if 'FreePrm2' not in data:
        data['FreePrm2'] = 0.
    if 'FreePrm3' not in data:
        data['FreePrm3'] = 0.
    if 'SurfRoughA' not in data and 'PWDR' in histName:
        data['SurfRoughA'] = [0.,False]
        data['SurfRoughB'] = [0.,False]
    if 'Trans' not in data and 'SASD' in histName:
        data['Trans'] = 1.0
    if 'SlitLen' not in data and 'SASD' in histName:
        data['SlitLen'] = 0.0
    if 'Shift' not in data:
        data['Shift'] = [0.0,False]
    if 'Transparency' not in data:
        data['Transparency'] = [0.0,False]
    data['InstrName'] = data.get('InstrName','')
#patch end
    labelLst,elemKeysLst,dspLst,refFlgElem = [],[],[],[]
    parms = SetupSampleLabels(histName,data.get('Type'),Inst['Type'][0])
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    topSizer.Add((-1,-1),1,wx.EXPAND,1)
    topSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Sample and Experimental Parameters'))
    topSizer.Add((-1,-1),1,wx.EXPAND,1)
    mainSizer.Add(topSizer,0,wx.EXPAND,1)
    nameSizer = wx.BoxSizer(wx.HORIZONTAL)
    nameSizer.Add(wx.StaticText(G2frame.dataDisplay,wx.ID_ANY,' Instrument Name'),
                0,WACV)
    nameSizer.Add((-1,-1),1,wx.EXPAND,1)
    instNameVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,data['InstrName'],
                              size=(200,-1),style=wx.TE_PROCESS_ENTER)        
    nameSizer.Add(instNameVal)
    instNameVal.Bind(wx.EVT_CHAR,OnNameVal)
    mainSizer.Add(nameSizer,0,wx.EXPAND,1)
    mainSizer.Add((5,5),0)
    labelLst.append('Instrument Name')
    elemKeysLst.append(['InstrName'])
    dspLst.append(None)
    refFlgElem.append(None)

    if 'PWDR' in histName:
        nameSizer = wx.BoxSizer(wx.HORIZONTAL)
        nameSizer.Add(wx.StaticText(G2frame.dataDisplay,wx.ID_ANY,' Diffractometer type: '),
                    0,WACV)
        if 'T' in Inst['Type'][0]:
            choices = ['Debye-Scherrer',]
        else:
            choices = ['Debye-Scherrer','Bragg-Brentano',]
        histoType = G2G.G2ChoiceButton(G2frame.dataDisplay,choices,
                    strLoc=data,strKey='Type',
                    onChoice=OnHistoChange)
        nameSizer.Add(histoType)
        mainSizer.Add(nameSizer,0,wx.EXPAND,1)
        mainSizer.Add((5,5),0)

    parmSizer = wx.FlexGridSizer(0,2,5,0)
    for key,lbl,nDig in parms:
        labelLst.append(lbl.strip().strip(':').strip())
        dspLst.append(nDig)
        if 'list' in str(type(data[key])):
            parmRef = G2G.G2CheckBox(G2frame.dataDisplay,' '+lbl,data[key],1)
            parmSizer.Add(parmRef,0,WACV|wx.EXPAND)
            parmVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data[key],0,
                nDig=nDig,typeHint=float,OnLeave=AfterChange)
            elemKeysLst.append([key,0])
            refFlgElem.append([key,1])
        else:
            parmSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' '+lbl),
                0,WACV|wx.EXPAND)
            parmVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,key,
                typeHint=float,OnLeave=AfterChange)
            elemKeysLst.append([key])
            refFlgElem.append(None)
        parmSizer.Add(parmVal,1,wx.EXPAND)
    Info = {}
        
    for key in ('FreePrm1','FreePrm2','FreePrm3'):
        parmVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,Controls,key,typeHint=str,
                                        notBlank=False)
        parmSizer.Add(parmVal,1,wx.EXPAND)
        parmVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,key,typeHint=float)
        parmSizer.Add(parmVal,1,wx.EXPAND)
        labelLst.append(Controls[key])
        dspLst.append(None)
        elemKeysLst.append([key])
        refFlgElem.append(None)
        
    mainSizer.Add(parmSizer,1,wx.EXPAND)
    mainSizer.Add((0,5),0)    
    if histName[:4] in ['SASD',]:
        rho = [0.,0.]
        anomrho = [0.,0.]
        mu = 0.
        subSizer = wx.FlexGridSizer(0,4,5,5)
        Substances = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Substances'))
        for id,item in enumerate(data['Materials']):
            subSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Material: '),0,WACV)
            matsel = wx.ComboBox(G2frame.dataDisplay,value=item['Name'],choices=Substances['Substances'].keys(),
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Info[matsel.GetId()] = id
            matsel.Bind(wx.EVT_COMBOBOX,OnMaterial)        
            subSizer.Add(matsel,0,WACV)
            subSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Volume fraction: '),0,WACV)
            volfrac = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,item,'VolFrac',
                min=0.,max=1.,nDig=(10,3),typeHint=float,OnLeave=OnVolFrac)
            subSizer.Add(volfrac,0,WACV)
            material = Substances['Substances'][item['Name']]
            mu += item['VolFrac']*material.get('XAbsorption',0.)
            rho[id] = material['Scatt density']
            anomrho[id] = material.get('XAnom density',0.)
        data['Contrast'] = [(rho[1]-rho[0])**2,(anomrho[1]-anomrho[0])**2]
        mainSizer.Add(subSizer,0)
        conSizer = wx.BoxSizer(wx.HORIZONTAL)
        conSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Contrast: %10.2f '%(data['Contrast'][0])),0,WACV)
        conSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Anom. Contrast: %10.2f '%(data['Contrast'][1])),0,WACV)
        mut =  mu*data['Thick']
        conSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Transmission (calc): %10.3f  '%(np.exp(-mut))),0,WACV)
        mainSizer.Add(conSizer,0)
    
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(Size)
                
################################################################################
#####  Indexing Peaks
################################################################################           
       
def UpdateIndexPeaksGrid(G2frame, data):
    '''respond to selection of PWDR Index Peak List data
    tree item.
    '''
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    IndexId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List')
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    limitId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits')
    Limits = G2frame.PatternTree.GetItemPyData(limitId)
    def RefreshIndexPeaksGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        peaks = G2frame.IndexPeaksTable.GetData()
        if c == 2:
            if peaks[r][c]:
                peaks[r][c] = False
            else:
                peaks[r][c] = True
            G2frame.IndexPeaksTable.SetData(peaks)
            G2frame.PatternTree.SetItemPyData(IndexId,[peaks,data[1]])
            G2frame.dataDisplay.ForceRefresh()
            if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                G2plt.PlotPowderLines(G2frame)
            else:
                G2plt.PlotPatterns(G2frame,plotType='PWDR')
            
    def OnReload(event):
        peaks = []
        sigs = []
        Peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
        for ip,peak in enumerate(Peaks['peaks']):
            dsp = G2lat.Pos2dsp(Inst,peak[0])
            peaks.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])    #SS?
            try:
                sig = Peaks['sigDict']['pos'+str(ip)]
            except KeyError:
                sig = 0.
            sigs.append(sig)
        data = [peaks,sigs]
        G2frame.PatternTree.SetItemPyData(IndexId,data)
        UpdateIndexPeaksGrid(G2frame,data)
        
    def KeyEditPickGrid(event):
        colList = G2frame.dataDisplay.GetSelectedCols()
        data = G2frame.PatternTree.GetItemPyData(IndexId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif colList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.IndexPeaksTable.GetColLabelValue(col) in ['use',]:
                    if key == 89: #'Y'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col]=False
                    elif key == 83: # 'S'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col] = not data[0][row][col]
                        
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    if 'PWD' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.IndPeaksMenu)
        G2frame.Bind(wx.EVT_MENU, OnReload, id=G2gd.wxID_INDXRELOAD)
    G2frame.dataFrame.IndexPeaks.Enable(False)
    G2frame.IndexPeaksTable = []
    if len(data[0]):
        G2frame.dataFrame.IndexPeaks.Enable(True)
        Unit = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List'))
        if Unit:
            if len(Unit) == 4:  #patch
                Unit.append({})
            controls,bravais,cellist,dmin,ssopt = Unit
            if 'T' in Inst['Type'][0]:   #TOF - use other limit!
                dmin = G2lat.Pos2dsp(Inst,Limits[1][0])
            else:
                dmin = G2lat.Pos2dsp(Inst,Limits[1][1])
            G2frame.HKL = []
            if ssopt.get('Use',False):
                cell = controls[6:12]
                A = G2lat.cell2A(cell)
                ibrav = bravaisSymb.index(controls[5])
                spc = controls[13]
                SGData = G2spc.SpcGroup(spc)[1]
                SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
                Vec = ssopt['ModVec']
                maxH = ssopt['maxH']
                G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,Vec,maxH,A)
                G2frame.HKL = np.array(G2frame.HKL)
                data[0] = G2indx.IndexSSPeaks(data[0],G2frame.HKL)[1]
            else:        #select cell from table - no SS
                for i,cell in enumerate(cellist):
                    if cell[-2]:
                        ibrav = cell[2]
                        A = G2lat.cell2A(cell[3:9])
                        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                        for hkl in G2frame.HKL:
                            hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3]))
                        G2frame.HKL = np.array(G2frame.HKL)
                        data[0] = G2indx.IndexPeaks(data[0],G2frame.HKL)[1]
                        break
    rowLabels = []
    for i in range(len(data[0])): rowLabels.append(str(i+1))
    colLabels = ['position','intensity','use','indexed','h','k','l','d-obs','d-calc']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',]+2*[wg.GRID_VALUE_BOOL,]+ \
        3*[wg.GRID_VALUE_LONG,]+2*[wg.GRID_VALUE_FLOAT+':10,5',]
    if len(data[0]) and len(data[0][0]) > 9:
        colLabels = ['position','intensity','use','indexed','h','k','l','m','d-obs','d-calc']
        Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',]+2*[wg.GRID_VALUE_BOOL,]+ \
            4*[wg.GRID_VALUE_LONG,]+2*[wg.GRID_VALUE_FLOAT+':10,5',]
    G2frame.PatternTree.SetItemPyData(IndexId,data)
    G2frame.IndexPeaksTable = G2G.Table(data[0],rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Index Peak List')
    G2frame.dataFrame.currentGrids = []
    G2frame.dataDisplay = G2G.GSGrid(parent=G2frame.dataFrame)                
    G2frame.dataDisplay.SetTable(G2frame.IndexPeaksTable, True)
    XY = []
    Sigs = []
    for r in range(G2frame.dataDisplay.GetNumberRows()):
        for c in range(G2frame.dataDisplay.GetNumberCols()):
            if c == 2:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
        if data[0][r][2] and data[0][r][3]:
            XY.append([data[0][r][-1],data[0][r][0]])
            try:
                sig = data[1][r]
            except IndexError:
                sig = 0.
            Sigs.append(sig)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, RefreshIndexPeaksGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)                 
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([490,300])
    if len(XY):
        XY = np.array(XY)
        G2plt.PlotCalib(G2frame,Inst,XY,Sigs,newPlot=True)
    G2frame.dataFrame.SendSizeEvent()
      
################################################################################
#####  Unit cells
################################################################################           
       
def UpdateUnitCellsGrid(G2frame, data):
    '''respond to selection of PWDR Unit Cells data tree item.
    '''
    UnitCellsId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
    SPGlist = G2spc.spglist
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm','P4/mmm',
        'Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    spaceGroups = ['F m 3 m','I m 3 m','P m 3 m','R 3 m','P 6/m m m','I 4/m m m',
        'P 4/m m m','F m m m','I m m m','C m m m','P m m m','C 2/m','P 2/m','P -1']
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    Limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits'))[1]
    if 'C' in Inst['Type'][0] or 'PKS' in Inst['Type'][0]:
        wave = G2mth.getWave(Inst)
        dmin = G2lat.Pos2dsp(Inst,Limits[1])
    else:
        difC = Inst['difC'][1]
        dmin = G2lat.Pos2dsp(Inst,Limits[0])
    
    def SetLattice(controls):
        ibrav = bravaisSymb.index(controls[5])
        if controls[5] in ['Fm3m','Im3m','Pm3m']:
            controls[7] = controls[8] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
        elif controls[5] in ['R3m','P6/mmm','I4/mmm','P4/mmm']:
            controls[7] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
            if controls[5] in ['R3-H','P6/mmm']:
                controls[11] = 120.
        elif controls[5] in ['Fmmm','Immm','Cmmm','Pmmm']:
            controls[9] = controls[10] = controls[11] = 90.
        elif controls[5] in ['C2/m','P2/m']:
            controls[9] = controls[11] = 90.  # b unique
        if len(controls) < 13: controls.append(0)
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        return ibrav
        
    def OnNcNo(event):
        controls[2] = NcNo.GetValue()
        
    def OnIfX20(event):
        G2frame.ifX20 = x20.GetValue()
        
    def OnBravais(event):
        Obj = event.GetEventObject()
        bravais[bravList.index(Obj.GetId())] = Obj.GetValue()
                
    def OnZeroVar(event):
        controls[0] = zeroVar.GetValue()
        
    def OnSSopt(event):
        if controls[5] in ['Fm3m','Im3m','Pm3m']:
            SSopt.SetValue(False)
            G2frame.ErrorDialog('Cubic lattice', 'Superlattice not allowed for a cubic lattice')
            return
        ssopt['Use'] = SSopt.GetValue()
        if 'ssSymb' not in ssopt:
            ssopt.update({'ssSymb':'(abg)','ModVec':[0.1,0.1,0.1],'maxH':1})
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnSelMG(event):
        ssopt['ssSymb'] = selMG.GetValue()
        Vec = ssopt['ModVec']
        modS = G2spc.splitSSsym(ssopt['ssSymb'])[0]
        ssopt['ModVec'] = G2spc.SSGModCheck(Vec,modS)[0]
        print ' Selecting: ',controls[13],ssopt['ssSymb'], 'maxH:',ssopt['maxH']
        OnHklShow(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnModVal(invalid,value,tc):
        OnHklShow(tc.event)
        
    def OnMoveMod(event):
        Obj = event.GetEventObject()
        ObjId = Obj.GetId()
        Id,valObj = Indx[ObjId]
        move = Obj.GetValue()*0.01
        Obj.SetValue(0)
        value = min(0.98,max(-0.98,float(valObj.GetValue())+move))
        valObj.SetValue('%.4f'%(value)) 
        ssopt['ModVec'][Id] = value
        OnHklShow(event)
        
    def OnMaxMH(event):
        ssopt['maxH'] = int(maxMH.GetValue())
        print ' Selecting: ',controls[13],ssopt['ssSymb'], 'maxH:',ssopt['maxH']
        OnHklShow(event)
        
    def OnFindMV(event):
        Peaks = np.copy(peaks[0])
        print ' Trying: ',controls[13],ssopt['ssSymb'], 'maxH:',1
        dlg = wx.ProgressDialog('Elapsed time','Modulation vector search',
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        try:
            ssopt['ModVec'],result = G2indx.findMV(Peaks,controls,ssopt,Inst,dlg)
            if len(result[0]) == 2:
                G2plt.PlotXYZ(G2frame,result[2],1./result[3],labelX='a',labelY='g',
                    newPlot=True,Title='Modulation vector search')
        finally:
            dlg.Destroy()
        OnHklShow(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnBravSel(event):
        brav = bravSel.GetString(bravSel.GetSelection())
        controls[5] = brav
        controls[13] = SPGlist[brav][0]       
        wx.CallLater(100,UpdateUnitCellsGrid,G2frame,data)
        
    def OnSpcSel(event):
        controls[13] = spcSel.GetString(spcSel.GetSelection())
        G2frame.dataFrame.RefineCell.Enable(True)
        OnHklShow(event)
        
    def SetCellValue(Obj,ObjId,value):
        if controls[5] in ['Fm3m','Im3m','Pm3m']:
            controls[6] = controls[7] = controls[8] = value
            controls[9] = controls[10] = controls[11] = 90.0
            Obj.SetValue("%.5f"%(controls[6]))
        elif controls[5] in ['R3-H','P6/mmm','I4/mmm','P4/mmm']:
            if ObjId == 0:
                controls[6] = controls[7] = value
                Obj.SetValue("%.5f"%(controls[6]))
            else:
                controls[8] = value
                Obj.SetValue("%.5f"%(controls[8]))
            controls[9] = controls[10] = controls[11] = 90.0
            if controls[5] in ['R3-H','P6/mmm']:
                controls[11] = 120.
        elif controls[5] in ['Fmmm','Immm','Cmmm','Pmmm']:
            controls[6+ObjId] = value
            Obj.SetValue("%.5f"%(controls[6+ObjId]))
            controls[9] = controls[10] = controls[11] = 90.0
        elif controls[5] in ['C2/m','P2/m']:
            controls[9] = controls[11] = 90.0
            if ObjId != 3:
                controls[6+ObjId] = value
                Obj.SetValue("%.5f"%(controls[6+ObjId]))
            else:
                controls[10] = value
                Obj.SetValue("%.3f"%(controls[10]))
        else:
            controls[6+ObjId] = value
            if ObjId < 3:
                Obj.SetValue("%.5f"%(controls[6+ObjId]))
            else:
                Obj.SetValue("%.3f"%(controls[6+ObjId]))
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        volVal.SetValue("%.3f"%(controls[12]))
        
    def OnMoveCell(event):
        Obj = event.GetEventObject()
        ObjId = cellList.index(Obj.GetId())
        valObj = valDict[Obj.GetId()]
        if ObjId/2 < 3:
            move = Obj.GetValue()*0.01
        else:
            move = Obj.GetValue()*0.1
        Obj.SetValue(0)
        value = float(valObj.GetValue())+move  
        SetCellValue(valObj,ObjId/2,value)
        OnHklShow(event)
        
    def OnExportCells(event):
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose Indexing Result csv file', pth, '', 
            'indexing result file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.csv'
                File = open(filename,'w')
                names = 'M20,X20,Bravais,a,b,c,alpha,beta,gamma,volume\n'
                File.write(names)
                fmt = '%.2f,%d,%s,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%.3f\n'
                for cell in cells:
                    File.write(fmt%(cell[0],cell[1],bravaisSymb[cell[2]], cell[3],cell[4],cell[5], cell[6],cell[7],cell[8],cell[9]))
                File.close()
        finally:
            dlg.Destroy()
        
    def OnCellChange(event):
        event.Skip()
        Obj = event.GetEventObject()
        ObjId = cellList.index(Obj.GetId())
        try:
            value = max(1.0,float(Obj.GetValue()))
        except ValueError:
            if ObjId/2 < 3:               #bad cell edge - reset
                value = controls[6+ObjId/2]
            else:                       #bad angle
                value = 90.
        SetCellValue(Obj,ObjId/2,value)
        OnHklShow(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnHklShow(event):
        PatternId = G2frame.PatternId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))
        controls,bravais,cells,dminx,ssopt = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
#        ibrav = bravaisSymb.index(controls[5])
        spc = controls[13]
        SGData = G2spc.SpcGroup(spc)[1]
        if ssopt.get('Use',False):
            SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
            Vec = ssopt['ModVec']
            maxH = ssopt['maxH']
            G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,Vec,maxH,A)
            peaks = [G2indx.IndexSSPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #keep esds from peak fit
            M20,X20 = G2indx.calc_M20SS(peaks[0],G2frame.HKL)
        else:
            if len(peaks[0]):
#                dmin = peaks[0][-1][7]
                G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
                peaks = [G2indx.IndexPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #keep esds from peak fit
                M20,X20 = G2indx.calc_M20(peaks[0],G2frame.HKL)
            else:
                M20 = X20 = 0.
                G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
        G2frame.HKL = np.array(G2frame.HKL)
        if len(G2frame.HKL):
            print ' new M20,X20: %.2f %d fraction found: %.3f'%(M20,X20,float(len(peaks[0]))/len(G2frame.HKL))
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'),peaks)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
            
    def OnSortCells(event):
        controls,bravais,cells,dminx,ssopt = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        c =  event.GetCol()
        if colLabels[c] == 'M20':
            cells = G2indx.sortM20(cells)
        elif colLabels[c] in ['X20','Bravais','a','b','c','alpha','beta','gamma','Volume']:
            if c == 1:
                c += 1  #X20 before Use
            cells = G2indx.sortCells(cells,c-1)     #an extra column (Use) not in cells
        else:
            return
        data = [controls,bravais,cells,dmin,ssopt]
        G2frame.PatternTree.SetItemPyData(UnitCellsId,data)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def CopyUnitCell(event):
        controls,bravais,cells,dminx,ssopt = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        for Cell in cells:
            if Cell[-2]:
                break
        cell = Cell[2:9]
        controls[4] = 1
        controls[5] = bravaisSymb[cell[0]]
        controls[6:12] = cell[1:8]
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        controls[13] = spaceGroups[bravaisSymb.index(controls[5])]
        G2frame.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin,ssopt])
        G2frame.dataFrame.RefineCell.Enable(True)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)        
                
    def RefineCell(event):
        
        def cellPrint(ibrav,A):
            cell = G2lat.A2cell(A)
            Vol = G2lat.calc_V(A)
            if ibrav in ['Fm3m','Im3m','Pm3m']:
                print " %s%10.6f" % ('a =',cell[0])
            elif ibrav in ['R3-H','P6/mmm','I4/mmm','P4/mmm']:
                print " %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],' c =',cell[2],' volume =',Vol)
            elif ibrav in ['P4/mmm','Fmmm','Immm','Cmmm','Pmmm']:
                print " %s%10.6f %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],' volume =',Vol)
            elif ibrav in ['C2/m','P2/m']:
                print " %s%10.6f %s%10.6f %s%10.6f %s%8.3f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],'beta =',cell[4],' volume =',Vol)
            else:
                print " %s%10.6f %s%10.6f %s%10.6f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2])
                print " %s%8.3f %s%8.3f %s%8.3f %s%12.3f" % ('alpha =',cell[3],'beta =',cell[4],'gamma =',cell[5],' volume =',Vol)
                
        def vecPrint(Vec):
            print ' %s %10.5f %10.5f %10.5f'%('Modulation vector:',Vec[0],Vec[1],Vec[2])
             
        PatternId = G2frame.PatternId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))
        if not len(peaks[0]):
            G2frame.ErrorDialog('No peaks!', 'Nothing to refine!')
            return        
        print ' Refine cell'
        controls,bravais,cells,dminx,ssopt = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        SGData = G2spc.SpcGroup(controls[13])[1]
        if 'C' in Inst['Type'][0] or 'PKS' in Inst['Type'][0]:
            if ssopt.get('Use',False):
                vecFlags = [True if x in ssopt['ssSymb'] else False for x in ['a','b','g']]
                SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
                G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,ssopt['ModVec'],ssopt['maxH'],A)
                peaks = [G2indx.IndexSSPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #put peak fit esds back in peaks
                Lhkl,M20,X20,Aref,Vec,Zero = \
                    G2indx.refinePeaksZSS(peaks[0],wave,Inst,SGData,SSGData,ssopt['maxH'],ibrav,A,ssopt['ModVec'],vecFlags,controls[1],controls[0])
            else:
                G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
                peaks = [G2indx.IndexPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #put peak fit esds back in peaks
                Lhkl,M20,X20,Aref,Zero = G2indx.refinePeaksZ(peaks[0],wave,ibrav,A,controls[1],controls[0])
        else:   
            if ssopt.get('Use',False):
                vecFlags = [True if x in ssopt['ssSymb'] else False for x in ['a','b','g']]
                SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
                G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,ssopt['ModVec'],ssopt['maxH'],A)
                peaks = [G2indx.IndexSSPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #put peak fit esds back in peaks
                Lhkl,M20,X20,Aref,Vec,Zero = \
                    G2indx.refinePeaksTSS(peaks[0],difC,Inst,SGData,SSGData,ssopt['maxH'],ibrav,A,ssopt['ModVec'],vecFlags,controls[1],controls[0])
            else:
                G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
                peaks = [G2indx.IndexPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #put peak fit esds back in peaks
                Lhkl,M20,X20,Aref,Zero = G2indx.refinePeaksT(peaks[0],difC,ibrav,A,controls[1],controls[0])            
        G2frame.HKL = np.array(G2frame.HKL)
        controls[1] = Zero
        controls[6:12] = G2lat.A2cell(Aref)
        controls[12] = G2lat.calc_V(Aref)
        cells = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2]
        for cell in cells:
            cell[-2] = False
        cells.insert(0,[M20,X20,ibrav]+controls[6:13]+[True,False])
        if ssopt.get('Use',False):
            ssopt['ModVec'] = Vec
            G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,ssopt['ModVec'],ssopt['maxH'],A)
        else:
            G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
        data = [controls,bravais,cells,dmin,ssopt]
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
        print " %s%10.3f" % ('refinement M20 = ',M20)
        print ' unindexed lines = ',X20
        cellPrint(controls[5],Aref)
        ip = 4
        if ssopt.get('Use',False):
            vecPrint(Vec)
            ip = 5
        for hkl in G2frame.HKL:
            hkl[ip] = G2lat.Dsp2pos(Inst,hkl[ip-1])+controls[1]
        G2frame.HKL = np.array(G2frame.HKL)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnIndexPeaks(event):
        PatternId = G2frame.PatternId    
        print 'Peak Indexing'
        keepcells = []
        try:
            controls,bravais,cells,dminx,ssopt = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
            for cell in cells:
                if cell[11]:
                    cell[10] = False    #clear selection flag on keepers
                    keepcells.append(cell)
        except IndexError:
            pass
        except ValueError:
            G2frame.ErrorDialog('Error','Need to set controls in Unit Cell List first')
            return
        if ssopt.get('Use',False):
            G2frame.ErrorDialog('Super lattice error','Indexing not available for super lattices')
            return
        if True not in bravais:
            G2frame.ErrorDialog('Error','No Bravais lattices selected')
            return
        if not len(peaks[0]):
            G2frame.ErrorDialog('Error','Index Peak List is empty')
            return
        if len(peaks[0][0]) > 9:
            G2frame.ErrorDialog('Error','You need to reload Index Peaks List first')
            return
        G2frame.dataFrame.CopyCell.Enable(False)
        G2frame.dataFrame.RefineCell.Enable(False)
        dlg = wx.ProgressDialog("Generated reflections",'0 '+" cell search for "+bravaisNames[ibrav],101, 
#            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_SKIP|wx.PD_CAN_ABORT) #desn't work in 32 bit versions
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        try:
            OK,dmin,newcells = G2indx.DoIndexPeaks(peaks[0],controls,bravais,dlg,G2frame.ifX20)
            if not OK:
                dlg.Destroy()
        finally:
            dlg.Destroy()
        cells = keepcells+newcells
        cells = G2indx.sortM20(cells)
        if OK:
            cells[0][10] = True         #select best M20
            data = [controls,bravais,cells,dmin,ssopt]
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                G2frame.HKL = G2lat.GenHBravais(dmin,bestCell[2],G2lat.cell2A(bestCell[3:9]))
                for hkl in G2frame.HKL:
                    hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
                G2frame.HKL = np.array(G2frame.HKL)
                if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame)
                else:
                    G2plt.PlotPatterns(G2frame)
            G2frame.dataFrame.CopyCell.Enable(True)
            G2frame.dataFrame.IndexPeaks.Enable(True)
            G2frame.dataFrame.MakeNewPhase.Enable(True)
            G2frame.ifX20 = True
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
                
    def RefreshUnitCellsGrid(event):
        data = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        cells,dminx = data[2:4]
        r,c =  event.GetRow(),event.GetCol()
        if cells:
            if c == 2:
                for i in range(len(cells)):
                    cells[i][-2] = False
                    UnitCellsTable.SetValue(i,c,False)
                UnitCellsTable.SetValue(r,c,True)
                gridDisplay.ForceRefresh()
                cells[r][-2] = True
                ibrav = cells[r][2]
                A = G2lat.cell2A(cells[r][3:9])
                G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                for hkl in G2frame.HKL:
                    hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
                G2frame.HKL = np.array(G2frame.HKL)
                if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame)
                else:
                    G2plt.PlotPatterns(G2frame)
            elif c == 11:
                if UnitCellsTable.GetValue(r,c):
                    UnitCellsTable.SetValue(r,c,False)
                    cells[r][c] = False
                else:
                    cells[r][c] = True
                    UnitCellsTable.SetValue(r,c,True)
                gridDisplay.ForceRefresh()
            G2frame.PatternTree.SetItemPyData(UnitCellsId,data)
        
    def MakeNewPhase(event):
        if not G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases'):
            sub = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PhaseName = dlg.GetValue()
                cells = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2]
                for Cell in cells:
                    if Cell[-2]:
                        break
                cell = Cell[2:10]        
                sub = G2frame.PatternTree.AppendItem(parent=sub,text=PhaseName)
                E,SGData = G2spc.SpcGroup(controls[13])
                G2frame.PatternTree.SetItemPyData(sub, \
                    G2IO.SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell[1:],Super=ssopt))
                Status.SetStatusText('Change space group from '+str(controls[13])+' if needed')
        finally:
            dlg.Destroy()
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.IndexMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU, OnIndexPeaks, id=G2gd.wxID_INDEXPEAKS)
    G2frame.Bind(wx.EVT_MENU, CopyUnitCell, id=G2gd.wxID_COPYCELL)
    G2frame.Bind(wx.EVT_MENU, RefineCell, id=G2gd.wxID_REFINECELL)
    G2frame.Bind(wx.EVT_MENU, MakeNewPhase, id=G2gd.wxID_MAKENEWPHASE)
    G2frame.Bind(wx.EVT_MENU, OnExportCells, id=G2gd.wxID_EXPORTCELLS)
        
    controls,bravais,cells,dminx,ssopt = data
    if len(controls) < 13:              #add cell volume if missing
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    if len(controls) < 14:              #add space group used in indexing
        controls.append(spaceGroups[bravaisSymb.index(controls[5])])
    G2frame.PatternTree.SetItemPyData(UnitCellsId,data)            #update with volume
    bravaisNames = ['Cubic-F','Cubic-I','Cubic-P','Trigonal-R','Trigonal/Hexagonal-P',
        'Tetragonal-I','Tetragonal-P','Orthorhombic-F','Orthorhombic-I','Orthorhombic-C',
        'Orthorhombic-P','Monoclinic-C','Monoclinic-P','Triclinic']
    cellGUIlist = [[[0,1,2],4,zip([" Unit cell: a = "," Vol = "],["%.5f","%.3f"],[True,False],[0,0])],
    [[3,4,5,6],6,zip([" Unit cell: a = "," c = "," Vol = "],["%.5f","%.5f","%.3f"],[True,True,False],[0,2,0])],
    [[7,8,9,10],8,zip([" Unit cell: a = "," b = "," c = "," Vol = "],["%.5f","%.5f","%.5f","%.3f"],
        [True,True,True,False],[0,1,2,0])],
    [[11,12],10,zip([" Unit cell: a = "," b = "," c = "," beta = "," Vol = "],
        ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
    [[13,],8,zip([" Unit cell: a = "," b = "," c = "," Vol = "," alpha = "," beta = "," gamma = "],
        ["%.5f","%.5f","%.5f","%.3f","%.3f","%.3f","%.3f"],
        [True,True,True,False,True,True,True],[0,1,2,0,3,4,5])]]
    
    G2frame.dataFrame.SetLabel('Unit Cells List')
    G2frame.dataFrame.IndexPeaks.Enable(False)
    peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
    if peaks:
        G2frame.dataFrame.IndexPeaks.Enable(True)
    G2frame.dataFrame.RefineCell.Enable(False)
    if controls[12] > 1.0:                               #if a "real" volume (i.e. not default)
        G2frame.dataFrame.RefineCell.Enable(True)    
    G2frame.dataFrame.CopyCell.Enable(False)
    G2frame.dataFrame.MakeNewPhase.Enable(False)        
    G2frame.dataFrame.ExportCells.Enable(False)
    if cells:
        G2frame.dataFrame.CopyCell.Enable(True)
        G2frame.dataFrame.MakeNewPhase.Enable(True)
        G2frame.dataFrame.ExportCells.Enable(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Indexing controls: '),0,WACV)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(0,5,5,5)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Max Nc/Nobs '),0,WACV)
    NcNo = wx.SpinCtrl(G2frame.dataDisplay)
    NcNo.SetRange(2,8)
    NcNo.SetValue(controls[2])
    NcNo.Bind(wx.EVT_SPINCTRL,OnNcNo)
    littleSizer.Add(NcNo,0,WACV)
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Start Volume '),0,WACV)
    startVol = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,controls,3,typeHint=int,min=25)
    littleSizer.Add(startVol,0,WACV)
    x20 = wx.CheckBox(G2frame.dataDisplay,label='Use M20/(X20+1)?')
    x20.SetValue(G2frame.ifX20)
    x20.Bind(wx.EVT_CHECKBOX,OnIfX20)
    littleSizer.Add(x20,0,WACV)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Select Bravais Lattices for indexing: '),
        0,WACV)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(0,7,5,5)
    bravList = []
    bravs = zip(bravais,bravaisNames)
    for brav,bravName in bravs:
        bravCk = wx.CheckBox(G2frame.dataDisplay,label=bravName)
        bravList.append(bravCk.GetId())
        bravCk.SetValue(brav)
        bravCk.Bind(wx.EVT_CHECKBOX,OnBravais)
        littleSizer.Add(bravCk,0,WACV)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Cell Test & Refinement: '),0,WACV)
    mainSizer.Add((5,5),0)
    littleSizer = wx.BoxSizer(wx.HORIZONTAL)
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Bravais lattice "),0,WACV)
    bravSel = wx.Choice(G2frame.dataDisplay,choices=bravaisSymb)
    bravSel.SetSelection(bravaisSymb.index(controls[5]))
    bravSel.Bind(wx.EVT_CHOICE,OnBravSel)
    littleSizer.Add(bravSel,0,WACV)
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Space group "),0,WACV)
    spcSel = wx.Choice(G2frame.dataDisplay,choices=SPGlist[controls[5]])
    spcSel.SetSelection(SPGlist[controls[5]].index(controls[13]))
    spcSel.Bind(wx.EVT_CHOICE,OnSpcSel)
    littleSizer.Add(spcSel,0,WACV)
    if ssopt.get('Use',False):        #zero for super lattice doesn't work!
        controls[0] = False
    else:
        littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Zero offset"),0,WACV)
        zero = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,controls,1,nDig=(10,4),typeHint=float,min=-5.,max=5.)
        littleSizer.Add(zero,0,WACV)
        zeroVar = wx.CheckBox(G2frame.dataDisplay,label="Refine?")
        zeroVar.SetValue(controls[0])
        zeroVar.Bind(wx.EVT_CHECKBOX,OnZeroVar)
        littleSizer.Add(zeroVar,0,WACV)
    SSopt = wx.CheckBox(G2frame.dataDisplay,label="Super lattice?")
    SSopt.SetValue(ssopt.get('Use',False))
    SSopt.Bind(wx.EVT_CHECKBOX,OnSSopt)
    littleSizer.Add(SSopt,0,WACV)
    hklShow = wx.Button(G2frame.dataDisplay,label="Show hkl positions")
    hklShow.Bind(wx.EVT_BUTTON,OnHklShow)
    littleSizer.Add(hklShow,0,WACV)
    mainSizer.Add(littleSizer,0)
    
    mainSizer.Add((5,5),0)
    ibrav = SetLattice(controls)
    for cellGUI in cellGUIlist:
        if ibrav in cellGUI[0]:
            useGUI = cellGUI
    cellList = []
    valDict = {}
    littleSizer = wx.FlexGridSizer(0,useGUI[1],5,5)
    for txt,fmt,ifEdit,Id in useGUI[2]:
        littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=txt),0,WACV)
        if ifEdit:          #a,b,c,etc.
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
            cellVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[6+Id])),style=wx.TE_PROCESS_ENTER)
            cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
            cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
            valSizer = wx.BoxSizer(wx.HORIZONTAL)
            valSizer.Add(cellVal,0,WACV)
            cellSpin = wx.SpinButton(G2frame.dataDisplay,style=wx.SP_VERTICAL,size=wx.Size(20,20))
            cellSpin.SetValue(0)
            cellSpin.SetRange(-1,1)
            cellSpin.Bind(wx.EVT_SPIN, OnMoveCell)
            valSizer.Add(cellSpin,0,WACV)
            littleSizer.Add(valSizer,0,WACV)
            cellList.append(cellVal.GetId())
            cellList.append(cellSpin.GetId())
            valDict[cellSpin.GetId()] = cellVal
        else:               #volume
            volVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[12])),style=wx.TE_READONLY)
            volVal.SetBackgroundColour(VERY_LIGHT_GREY)
            littleSizer.Add(volVal,0,WACV)
    mainSizer.Add(littleSizer,0)
    if ssopt.get('Use',False):        #super lattice display
        indChoice = ['1','2','3','4',]
        SpSg = controls[13]
        ssChoice = G2spc.ssdict[SpSg]
        if ssopt['ssSymb'] not in ssChoice:
            ssopt['ssSymb'] = ssChoice[0]
        ssSizer = wx.BoxSizer(wx.HORIZONTAL)
        ssSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Supersymmetry space group: '+SpSg+' '),0,WACV)
        selMG = wx.ComboBox(G2frame.dataDisplay,value=ssopt['ssSymb'],
                choices=ssChoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        selMG.Bind(wx.EVT_COMBOBOX, OnSelMG)
        ssSizer.Add(selMG,0,WACV)
        ssSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Mod. vector: '),0,WACV)
        modS = G2spc.splitSSsym(ssopt['ssSymb'])[0]
        ssopt['ModVec'],ifShow = G2spc.SSGModCheck(ssopt['ModVec'],modS)
        Indx = {}
        for i,[val,show] in enumerate(zip(ssopt['ModVec'],ifShow)):
            if show:
                valSizer = wx.BoxSizer(wx.HORIZONTAL)
                modVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,ssopt['ModVec'],i,
                    min=-.98,max=.98,nDig=(10,4),typeHint=float,OnLeave=OnModVal)
                valSizer.Add(modVal,0,WACV)
                modSpin = wx.SpinButton(G2frame.dataDisplay,style=wx.SP_VERTICAL,size=wx.Size(20,20))
                modSpin.SetValue(0)
                modSpin.SetRange(-1,1)
                modSpin.Bind(wx.EVT_SPIN, OnMoveMod)
                valSizer.Add(modSpin,0,WACV)
                ssSizer.Add(valSizer,0,WACV)
                Indx[modVal.GetId()] = i
                Indx[modSpin.GetId()] = [i,modVal]
            else:
                modVal = wx.TextCtrl(G2frame.dataDisplay,value=('%.3f'%(val)),
                    size=wx.Size(50,20),style=wx.TE_READONLY)
                modVal.SetBackgroundColour(VERY_LIGHT_GREY)
                ssSizer.Add(modVal,0,WACV)
        ssSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Max. M: '),0,WACV)
        maxMH = wx.ComboBox(G2frame.dataDisplay,value=str(ssopt['maxH']),
            choices=indChoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        maxMH.Bind(wx.EVT_COMBOBOX, OnMaxMH)
        ssSizer.Add(maxMH,0,WACV)
        findMV = wx.Button(G2frame.dataDisplay,label="Find mod. vec.?")
        findMV.Bind(wx.EVT_BUTTON,OnFindMV)
        ssSizer.Add(findMV,0,WACV)
        mainSizer.Add(ssSizer,0)

    G2frame.dataFrame.currentGrids = []
    if cells:
        mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='\n Indexing Result:'),0,WACV)
        rowLabels = []
        colLabels = ['M20','X20','use','Bravais','a','b','c','alpha','beta','gamma','Volume','Keep']
        Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+ \
            [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL]
        table = []
        for cell in cells:
            rowLabels.append('')
            row = cell[0:2]+[cell[-2]]+[bravaisSymb[cell[2]]]+cell[3:10]+[cell[11],]
            if cell[-2]:
                A = G2lat.cell2A(cell[3:9])
                G2frame.HKL = G2lat.GenHBravais(dmin,cell[2],A)
                for hkl in G2frame.HKL:
                    hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
                G2frame.HKL = np.array(G2frame.HKL)
            table.append(row)
        UnitCellsTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        gridDisplay = G2G.GSGrid(G2frame.dataDisplay)
        gridDisplay.SetTable(UnitCellsTable, True)
        G2frame.dataFrame.CopyCell.Enable(True)
        gridDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,RefreshUnitCellsGrid)
        gridDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnSortCells)
        gridDisplay.SetMargins(0,0)
        gridDisplay.SetRowLabelSize(0)
        gridDisplay.AutoSizeColumns(False)
        for r in range(gridDisplay.GetNumberRows()):
            for c in range(gridDisplay.GetNumberCols()):
                if c == 2:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=True)
        mainSizer.Add(gridDisplay,0,WACV)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    Size[1] = min(600,Size[1])
    G2frame.dataFrame.setSizePosLeft(Size)    
    
################################################################################
#####  Reflection list
################################################################################           
       
def UpdateReflectionGrid(G2frame,data,HKLF=False,Name=''):
    '''respond to selection of PWDR Reflections data tree item by displaying
    a table of reflections in the data window.
    '''
    Controls = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Controls'))
    dMin = 0.05
    if 'UsrReject' in Controls:
        dMin = Controls['UsrReject'].get('MinD',0.05)
    def OnPlotHKL(event):
        '''Plots a layer of reflections
        '''
        phaseName = G2frame.RefList
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetPatternTreeItemId(G2frame,pId,phaseName)
            General = G2frame.PatternTree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
            SuperVec = General.get('SuperVec',[])
        else:
            Super = 0
            SuperVec = []       
        if 'list' in str(type(data)):   #single crystal data is 2 dict in list
            refList = data[1]['RefList']
        else:                           #powder data is a dict of dicts; each same structure as SC 2nd dict
            if 'RefList' in data[phaseName]:
                refList = np.array(data[phaseName]['RefList'])
            else:
                wx.MessageBox('No reflection list - do Refine first',caption='Reflection plotting')
                return
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        controls = {'Type' : 'Fo','ifFc' : True,'HKLmax' : Hmax,'HKLmin' : Hmin,
            'FoMax' : FoMax,'Zone' : '001','Layer' : 0,'Scale' : 1.0,'Super':Super,'SuperVec':SuperVec}
        G2plt.PlotSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
        
    def OnPlot3DHKL(event):
        '''Plots the reflections in 3D
        '''
        phaseName = G2frame.RefList
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetPatternTreeItemId(G2frame,pId,phaseName)
            General = G2frame.PatternTree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
            SuperVec = General.get('SuperVec',[])
        else:
            Super = 0
            SuperVec = []       
        if 'list' in str(type(data)):   #single crystal data is 2 dict in list
            refList = data[1]['RefList']
        else:                           #powder data is a dict of dicts; each same structure as SC 2nd dict
            if 'RefList' in data[phaseName]:
                refList = np.array(data[phaseName]['RefList'])
            else:
                wx.MessageBox('No reflection list - do Refine first',caption='Reflection plotting')
                return
        refList.T[3+Super] = np.where(refList.T[4+Super]<dMin,-refList.T[3+Super],refList.T[3+Super])
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        Vpoint = np.array([int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))])
        controls = {'Type':'Fosq','Iscale':False,'HKLmax':Hmax,'HKLmin':Hmin,'Zone':False,'viewKey':'L',
            'FoMax' : FoMax,'Scale' : 1.0,'Drawing':{'viewPoint':[Vpoint,[]],'default':Vpoint[:],
            'backColor':[0,0,0],'depthFog':False,'Zclip':10.0,'cameraPos':10.,'Zstep':0.05,'viewUp':[0,1,0],
            'Scale':1.0,'oldxy':[],'viewDir':[0,0,1]},'Super':Super,'SuperVec':SuperVec}
        G2plt.Plot3DSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
        
    def MakeReflectionTable(phaseName):
        '''Returns a wx.grid table (G2G.Table) containing a list of all reflections
        for a phase.        
        '''
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetPatternTreeItemId(G2frame,pId,phaseName)
            if not phaseId:         #phase deleted 
                return None
            General = G2frame.PatternTree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
        else:
            Super = 0
        rowLabels = []
        if HKLF:
            refList = data[1]['RefList']
            refs = refList
        else:
            if len(data) > 1:
                G2frame.dataFrame.SelectPhase.Enable(True)
            try:            #patch for old reflection lists
                if not len(data[phaseName]):
                    return None
                refList = np.array(data[phaseName]['RefList'])
                I100 = refList.T[8+Super]*refList.T[11+Super]
            except TypeError:
                refList = np.array([refl[:11+Super] for refl in data[phaseName]])
                I100 = refList.T[8+Super]*np.array([refl[11+Super] for refl in data[phaseName]])
            Imax = np.max(I100)
            if Imax:
                I100 *= 100.0/Imax
            if 'C' in Inst['Type'][0]:
                refs = np.vstack((refList.T[:15+Super],I100)).T
            elif 'T' in Inst['Type'][0]:
                refs = np.vstack((refList.T[:18+Super],I100)).T
        rowLabels = [str(i) for i in range(len(refs))]
        Types = (4+Super)*[wg.GRID_VALUE_LONG,]+4*[wg.GRID_VALUE_FLOAT+':10,4',]+ \
            2*[wg.GRID_VALUE_FLOAT+':10,2',]+[wg.GRID_VALUE_FLOAT+':10,3',]+ \
            [wg.GRID_VALUE_FLOAT+':10,3',]
        if HKLF:
            colLabels = ['H','K','L','twin','d','Fosq','sig','Fcsq','FoTsq','FcTsq','phase','ExtC',]
            if 'T' in Inst['Type'][0]:
                colLabels = ['H','K','L','twin','d','Fosq','sig','Fcsq','FoTsq','FcTsq','phase','ExtC','wave','tbar']
                Types += 2*[wg.GRID_VALUE_FLOAT+':10,3',]
            if Super:
                colLabels.insert(3,'M')
        else:
            if 'C' in Inst['Type'][0]:
                colLabels = ['H','K','L','mul','d','pos','sig','gam','Fosq','Fcsq','phase','Icorr','Prfo','Trans','ExtP','I100']
                Types += 4*[wg.GRID_VALUE_FLOAT+':10,3',]
            elif 'T' in Inst['Type'][0]:
                colLabels = ['H','K','L','mul','d','pos','sig','gam','Fosq','Fcsq','phase','Icorr','alp','bet','wave','Prfo','Abs','Ext','I100']
                Types += 7*[wg.GRID_VALUE_FLOAT+':10,3',]
            if Super:
                colLabels.insert(3,'M')
        refs.T[3+Super] = np.where(refs.T[4+Super]<dMin,-refs.T[3+Super],refs.T[3+Super])
        return G2G.Table(refs,rowLabels=rowLabels,colLabels=colLabels,types=Types)

    def ShowReflTable(phaseName):
        '''Posts a table of reflections for a phase, creating the table
        if needed using MakeReflectionTable
        '''
        def setBackgroundColors(im,it):
            for r in range(G2frame.refTable[phaseName].GetNumberRows()):
                if HKLF:
                    if float(G2frame.refTable[phaseName].GetCellValue(r,3+im)) <= 0.:
                        G2frame.refTable[phaseName].SetCellBackgroundColour(r,3+im,wx.RED)
                    Fosq = float(G2frame.refTable[phaseName].GetCellValue(r,5+im))
                    Fcsq = float(G2frame.refTable[phaseName].GetCellValue(r,7+im))
                    sig = float(G2frame.refTable[phaseName].GetCellValue(r,6+im))
                    rat = 11.
                    if sig:
                        rat = abs(Fosq-Fcsq)/sig
                    if  rat > 10.:
                        G2frame.refTable[phaseName].SetCellBackgroundColour(r,7+im,wx.RED)
                    elif rat > 3.0:
                        G2frame.refTable[phaseName].SetCellBackgroundColour(r,7+im,wx.Colour(255,255,0))
                else:   #PWDR
                    if float(G2frame.refTable[phaseName].GetCellValue(r,12+im+itof)) < 0.:
                        G2frame.refTable[phaseName].SetCellBackgroundColour(r,12+im+itof,wx.RED)
                    if float(G2frame.refTable[phaseName].GetCellValue(r,3+im)) < 0:
                        G2frame.refTable[phaseName].SetCellBackgroundColour(r,8+im,wx.RED)
                        
                                                  
        if not HKLF and not len(data[phaseName]):
            return          #deleted phase?
        G2frame.RefList = phaseName
        G2frame.dataFrame.SetLabel('Reflection List for '+phaseName)
        if HKLF:
            Status.SetStatusText('abs(DF)/sig > 10 red; > 3 yellow; twin < 0 (user rejected) red; twin=0 (sp. gp. absent) red')
        else:
            Status.SetStatusText('Prfo < 0. in red; if excluded Fosq in red & mul < 0')
        itof = 0
        if HKLF:
            im = data[1].get('Super',0)
        else:
            if 'T' in data[phaseName].get('Type',''):
                itof = 3
            im = data[phaseName].get('Super',0)
        # has this table already been displayed?
        if G2frame.refTable[phaseName].GetTable() is None:
            PeakTable = MakeReflectionTable(phaseName)
            G2frame.refTable[phaseName].SetTable(PeakTable, True)
            G2frame.refTable[phaseName].EnableEditing(False)
            G2frame.refTable[phaseName].SetMargins(0,0)
            G2frame.refTable[phaseName].AutoSizeColumns(False)
            setBackgroundColors(im,itof)
        if HKLF:
            refList = np.array([refl[:6+im] for refl in data[1]['RefList']])
        else:
            refList = np.array([refl[:6+im] for refl in data[phaseName]['RefList']])
        G2frame.HKL = np.vstack((refList.T)).T    #build for plots
        # raise the tab (needed for 1st use and from OnSelectPhase)
        for PageNum in range(G2frame.dataDisplay.GetPageCount()):
            if phaseName == G2frame.dataDisplay.GetPageText(PageNum):
                G2frame.dataDisplay.SetSelection(PageNum)
                break
        else:
            print phaseName
            print phases
            raise Exception("how did we not find a phase name?")
        
    def OnPageChanged(event):
        '''Respond to a press on a phase tab by displaying the reflections. This
        routine is needed because the reflection table may not have been created yet.
        '''
        page = event.GetSelection()
        phaseName = G2frame.dataDisplay.GetPageText(page)
        ShowReflTable(phaseName)

    def OnSelectPhase(event):
        '''For PWDR, selects a phase with a selection box. Called from menu.
        '''
        if len(phases) < 2: return
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Phase',phases)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                ShowReflTable(phases[sel])
        finally:
            dlg.Destroy()
            
    if not data:
        print 'No phases, no reflections'
        return
    if HKLF:
        G2frame.RefList = 1
        phaseName = IsHistogramInAnyPhase(G2frame,Name)
        if not phaseName:
            phaseName = 'Unknown'
        phases = [phaseName]
    else:
        phaseName = G2frame.RefList
        phases = data.keys()
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    if HKLF:
        G2gd.SetDataMenuBar(G2frame)
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ReflMenu)
        G2frame.Bind(wx.EVT_MENU, OnPlotHKL, id=G2gd.wxID_PWDHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=G2gd.wxID_PWD3DHKLPLOT)
        G2frame.dataFrame.SelectPhase.Enable(False)
    else:
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ReflMenu)
        G2frame.Bind(wx.EVT_MENU, OnSelectPhase, id=G2gd.wxID_SELECTPHASE)
        G2frame.Bind(wx.EVT_MENU, OnPlotHKL, id=G2gd.wxID_PWDHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=G2gd.wxID_PWD3DHKLPLOT)
        G2frame.dataFrame.SelectPhase.Enable(False)
            
    G2frame.dataDisplay = G2G.GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    G2frame.refTable = {}
    G2frame.dataFrame.currentGrids = []
    for tabnum,phase in enumerate(phases):
        if isinstance(data,list):           #single crystal HKLF
            G2frame.refTable[phase] = G2G.GSGrid(parent=G2frame.dataDisplay)
            G2frame.dataDisplay.AddPage(G2frame.refTable[phase],phase)
        elif len(data[phase]):              #else dict for PWDR
            G2frame.refTable[phase] = G2G.GSGrid(parent=G2frame.dataDisplay)
            G2frame.dataDisplay.AddPage(G2frame.refTable[phase],phase)
        else:       #cleanup deleted phase reflection lists
            del data[phase]
            if len(data):
                G2frame.RefList = data.keys()[0]
                phaseName = G2frame.RefList
            else:
                G2frame.RefList = ''
                phaseName = ''
    if phaseName:
        ShowReflTable(phaseName)
    G2frame.dataFrame.setSizePosLeft([550,350])
    G2frame.dataDisplay.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    
################################################################################
#####  SASD/REFD Substances 
################################################################################
           
def UpdateSubstanceGrid(G2frame,data):
    '''respond to selection of SASD/REFD Substance data tree item.
    '''
    import Substances as substFile
    
    def LoadSubstance(name):
        
        subst = substFile.Substances[name]
        ElList = subst['Elements'].keys()
        for El in ElList:
            Info = G2elem.GetAtomInfo(El.strip().capitalize())
            Info.update(subst['Elements'][El])
            data['Substances'][name]['Elements'][El] = Info
            if 'Volume' in subst:
                data['Substances'][name]['Volume'] = subst['Volume']
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
            elif 'Density' in subst:
                data['Substances'][name]['Density'] = subst['Density']
                data['Substances'][name]['Volume'] = \
                    G2mth.Den2Vol(data['Substances'][name]['Elements'],data['Substances'][name]['Density'])
            else:
                data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
            if 'X' in Inst['Type'][0]:
                data['Substances'][name]['Scatt density'] = \
                    G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            elif 'NC' in Inst['Type'][0]:
                isotopes = Info['Isotopes'].keys()
                isotopes.sort()
                data['Substances'][name]['Elements'][El]['Isotope'] = isotopes[-1]
                data['Substances'][name]['Scatt density'] = \
                    G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            data['Substances'][name]['XAnom density'] = recontrst
            data['Substances'][name]['XAbsorption'] = absorb
            data['Substances'][name]['XImag density'] = imcontrst
            
    def OnReloadSubstances(event):
        
        for name in data['Substances'].keys():
            if name not in ['vacuum','unit scatter']:
                if 'X' in Inst['Type'][0]:
                    data['Substances'][name]['Scatt density'] = \
                        G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                    recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                elif 'NC' in Inst['Type'][0]:
                    data['Substances'][name]['Scatt density'] = \
                        G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                    recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                data['Substances'][name]['XAnom density'] = recontrst
                data['Substances'][name]['XAbsorption'] = absorb
                data['Substances'][name]['XImag density'] = imcontrst
        UpdateSubstanceGrid(G2frame,data)
    
    def OnLoadSubstance(event):
        
        names = substFile.Substances.keys()
        names.sort()
        dlg = wx.SingleChoiceDialog(G2frame, 'Which substance?', 'Select substance', names, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                name = names[dlg.GetSelection()]
            else:
                return
        finally:
            dlg.Destroy()
            
        data['Substances'][name] = {'Elements':{},'Volume':1.0,'Density':1.0,
            'Scatt density':0.0,'Real density':0.0,'XAbsorption':0.0,'XImag density':0.0}
        LoadSubstance(name)            
        UpdateSubstanceGrid(G2frame,data)
        
    def OnCopySubstance(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy substances from\n'+hst[5:]+' to...',
            'Copy substances', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections(): 
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()        
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id, 'Instrument Parameters'))[0]
            wave = G2mth.getWave(Inst)
            ndata = copy.deepcopy(data)
            for name in ndata['Substances'].keys():
                if name not in ['vacuum','unit scatter']:
                    if 'X' in Inst['Type'][0]:
                        recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                    elif 'NC' in Inst['Type'][0]:
                        recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                    ndata['Substances'][name]['XAnom density'] = recontrst
                    ndata['Substances'][name]['XAbsorption'] = absorb
                    ndata['Substances'][name]['XImag density'] = imcontrst
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Substances'),ndata)
    
    def OnAddSubstance(event):
        dlg = wx.TextEntryDialog(None,'Enter a name for this substance','Substance Name Entry','New substance',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            Name = dlg.GetValue()
            data['Substances'][Name] = {'Elements':{},'Volume':1.0,'Density':1.0,
                'Scatt density':0.0,'XAnom density':0.,'XAbsorption':0.}
            AddElement(Name)
        dlg.Destroy()
        UpdateSubstanceGrid(G2frame,data)
        
    def OnDeleteSubstance(event):
        TextList = []
        for name in data['Substances']:
            if name not in ['vacuum','unit scatter']:
                TextList += [name,]
        if not TextList:
            return
        dlg = wx.SingleChoiceDialog(G2frame, 'Which substance?', 'Select substance to delete', TextList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                name = TextList[dlg.GetSelection()]
            else:
                return
        finally:
            dlg.Destroy()
        del(data['Substances'][name])
        UpdateSubstanceGrid(G2frame,data)        
                
    def OnAddElement(event):        
        TextList = []
        for name in data['Substances']:
            if name not in ['vacuum','unit scatter']:
                TextList += [name,]
        if not TextList:
            return
        dlg = wx.SingleChoiceDialog(G2frame, 'Which substance?', 'Select substance', TextList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                name = TextList[dlg.GetSelection()]
                AddElement(name)
            else:
                return
        finally:
            dlg.Destroy()
        UpdateSubstanceGrid(G2frame,data)
        
    def AddElement(name):
        ElList = data['Substances'][name]['Elements'].keys()
        dlg = G2elemGUI.PickElements(G2frame,ElList)
        if dlg.ShowModal() == wx.ID_OK:
            for El in dlg.Elem:
                El = El.strip().capitalize()
                Info = G2elem.GetAtomInfo(El)
                Info.update({'Num':1.})
                data['Substances'][name]['Elements'][El] = Info
                isotopes = Info['Isotopes'].keys()
                isotopes.sort()
                data['Substances'][name]['Elements'][El]['Isotope'] = isotopes[-1]
            data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
            data['Substances'][name]['Density'] = \
                G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
            if 'X' in Inst['Type'][0]:
                data['Substances'][name]['Scatt density'] = \
                    G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            elif 'NC' in Inst['Type'][0]:
                data['Substances'][name]['Scatt density'] = \
                    G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            data['Substances'][name]['XAnom density'] = recontrst
            data['Substances'][name]['XAbsorption'] = absorb
            data['Substances'][name]['XImag density'] = imcontrst
        dlg.Destroy()
        
    def OnDeleteElement(event):
        TextList = []
        for name in data['Substances']:
            if name not in ['vacuum','unit scatter']:
                TextList += [name,]
        if not TextList:
            return
        dlg = wx.SingleChoiceDialog(G2frame, 'Which substance?', 'Select substance', TextList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                name = TextList[dlg.GetSelection()]
            else:
                return
        finally:
            dlg.Destroy()
        ElList = data['Substances'][name]['Elements'].keys()
        if len(ElList):
            DE = G2elemGUI.DeleteElement(G2frame,ElList)
            if DE.ShowModal() == wx.ID_OK:
                El = DE.GetDeleteElement().strip().upper()
                del(data['Substances'][name]['Elements'][El])
                data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
                if 'X' in Inst['Type'][0]:
                    data['Substances'][name]['Scatt density'] = \
                        G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                    recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                elif 'NC' in Inst['Type'][0]:
                    data['Substances'][name]['Scatt density'] = \
                        G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                    recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
                data['Substances'][name]['XAnom density'] = recontrst
                data['Substances'][name]['XAbsorption'] = absorb
                data['Substances'][name]['XImag density'] = imcontrst
        UpdateSubstanceGrid(G2frame,data)
                
    def SubstSizer():
        
        def OnNum(invalid,value,tc):
            if invalid: return
            name,El,keyId = Indx[tc.GetId()]
            data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
            data['Substances'][name]['Density'] = \
                G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
            if 'X' in Inst['Type'][0]:
                recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            elif 'NC' in Inst['Type'][0]:
                recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            data['Substances'][name]['XAnom density'] = recontrst
            data['Substances'][name]['XAbsorption'] = absorb
            data['Substances'][name]['XImag density'] = imcontrst
            wx.CallAfter(UpdateSubstanceGrid,G2frame,data)
            
        def OnVolDen(invalid,value,tc):
            if invalid: return
            name,keyId = Indx[tc.GetId()]
            if keyId in 'Volume':
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],value)
            elif keyId in 'Density':
                data['Substances'][name]['Volume'] = \
                    G2mth.Den2Vol(data['Substances'][name]['Elements'],value)
            if 'X' in Inst['Type'][0]:
                data['Substances'][name]['Scatt density'] = \
                    G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            elif 'NC' in Inst['Type'][0]:
                data['Substances'][name]['Scatt density'] = \
                    G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            data['Substances'][name]['XAnom density'] = recontrst
            data['Substances'][name]['XAbsorption'] = absorb
            data['Substances'][name]['XImag density'] = imcontrst
            wx.CallAfter(UpdateSubstanceGrid,G2frame,data)
            
        def OnIsotope(event):
            Obj = event.GetEventObject()
            El,name = Indx[Obj.GetId()]
            data['Substances'][name]['Elements'][El]['Isotope'] = Obj.GetValue()
            recontrst,absorb,imcontrst = G2mth.NCScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)
            data['Substances'][name]['XAnom density'] = recontrst
            data['Substances'][name]['XAbsorption'] = absorb
            data['Substances'][name]['XImag density'] = imcontrst
            wx.CallAfter(UpdateSubstanceGrid,G2frame,data)
        
        Indx = {}
        substSizer = wx.BoxSizer(wx.VERTICAL)
        substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Substance list: wavelength: %.5fA'%(wave)),
            0,WACV)
        for name in data['Substances']:
            G2G.HorizontalLine(substSizer,G2frame.dataDisplay)    
            substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Data for '+name+':'),
                0,WACV)
            if name == 'vacuum':
                substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='        Not applicable'),
                    0,WACV)
            elif name == 'unit scatter':
                substSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Scattering density,f: %.3f *10%scm%s'%(data['Substances'][name]['Scatt density'],Pwr10,Pwrm2)),0,WACV)
            else:    
                elSizer = wx.FlexGridSizer(0,8,5,5)
                Substance = data['Substances'][name]
                Elems = Substance['Elements']
                for El in Elems:    #do elements as pull downs for isotopes for neutrons
                    elSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' '+El+': '),
                        0,WACV)
                    num = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Substances'][name]['Elements'][El],'Num',
                        nDig=(10,2,'f'),typeHint=float,OnLeave=OnNum)
                    Indx[num.GetId()] = [name,El,'Num']
                    elSizer.Add(num,0,WACV)
                    if 'N' in Inst['Type'][0]:
                        elSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Isotope: '),0,WACV)
                        isotopes = Elems[El]['Isotopes'].keys()
                        isotope = wx.ComboBox(G2frame.dataDisplay,choices=isotopes,value=Elems[El].get('Isotope','Nat. Abund.'),
                            style=wx.CB_READONLY|wx.CB_DROPDOWN)
                        Indx[isotope.GetId()] = [El,name]
                        isotope.Bind(wx.EVT_COMBOBOX,OnIsotope)
                        elSizer.Add(isotope,0,WACV)
                substSizer.Add(elSizer,0)
                vdsSizer = wx.FlexGridSizer(0,4,5,5)
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Volume: '),
                    0,WACV)
                vol = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Substances'][name],'Volume',nDig=(10,2),typeHint=float,OnLeave=OnVolDen)
                Indx[vol.GetId()] = [name,'Volume']
                vdsSizer.Add(vol,0,WACV)                
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Density: '),
                    0,WACV)
                den = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Substances'][name],'Density',nDig=(10,2),typeHint=float,OnLeave=OnVolDen)
                Indx[den.GetId()] = [name,'Density']
                vdsSizer.Add(den,0,WACV)
                substSizer.Add(vdsSizer,0)
                denSizer = wx.FlexGridSizer(0,2,0,0)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Scattering density,f'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=': %.3f *10%scm%s'%(Substance['Scatt density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Real density,f+f'"),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=': %.3f *10%scm%s'%(Substance['XAnom density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Imaginary density,f"'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=': %.3g *10%scm%s'%(Substance['XImag density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Absorption'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataDisplay,label=': %.3g cm%s'%(Substance['XAbsorption'],Pwrm1)),0,WACV)
                substSizer.Add(denSizer)
        return substSizer
            
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    wave = G2mth.getWave(Inst)
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()  # is this a ScrolledWindow? If so, bad!
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.SubstanceMenu)
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2frame.dataFrame.SetLabel('Substances')
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadSubstance, id=G2gd.wxID_LOADSUBSTANCE)    
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnReloadSubstances, id=G2gd.wxID_RELOADSUBSTANCES)    
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddSubstance, id=G2gd.wxID_ADDSUBSTANCE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopySubstance, id=G2gd.wxID_COPYSUBSTANCE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteSubstance, id=G2gd.wxID_DELETESUBSTANCE)    
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddElement, id=G2gd.wxID_ELEMENTADD)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2gd.wxID_ELEMENTDELETE)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(SubstSizer(),0)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    G2frame.dataFrame.setSizePosLeft(Size)    
       
################################################################################
#####  SASD Models 
################################################################################           
       
def UpdateModelsGrid(G2frame,data):
    '''respond to selection of SASD Models data tree item.
    '''
    #patches
    if 'Current' not in data:
        data['Current'] = 'Size dist.'
    if 'logBins' not in data['Size']:
        data['Size']['logBins'] = True
    if 'MinMaxDiam' in data['Size']:
        data['Size']['MinDiam'] = 50.
        data['Size']['MaxDiam'] = 10000.
        del data['Size']['MinMaxDiam']
    if isinstance(data['Size']['MaxEnt']['Sky'],float):
        data['Size']['MaxEnt']['Sky'] = -3
    if 'Power' not in data['Size']['IPG']:
        data['Size']['IPG']['Power'] = -1
    if 'Matrix' not in data['Particle']:
        data['Particle']['Matrix'] = {'Name':'vacuum','VolFrac':[0.0,False]}
    if 'BackFile' not in data:
        data['BackFile'] = ''
    #end patches
    
    def RefreshPlots(newPlot=False):
        PlotText = G2frame.G2plotNB.nb.GetPageText(G2frame.G2plotNB.nb.GetSelection())
        if 'Powder' in PlotText:
            G2plt.PlotPatterns(G2frame,plotType='SASD',newPlot=newPlot)
        elif 'Size' in PlotText:
            G2plt.PlotSASDSizeDist(G2frame)
                
    def OnAddModel(event):
        if data['Current'] == 'Particle fit':
            material = 'vacuum'
            if len(data['Particle']['Levels']):
                material = data['Particle']['Levels'][-1]['Controls']['Material']
            data['Particle']['Levels'].append({
                'Controls':{'FormFact':'Sphere','DistType':'LogNormal','Material':material,
                    'FFargs':{},'SFargs':{},'NumPoints':50,'Cutoff':0.01,'Contrast':0.0,
                    'SlitSmear':[0.0,False],'StrFact':'Dilute'},    #last 2 not used - future?
                'LogNormal':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[0.5,False],'MinSize':[10.,False],},
                'Gaussian':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[300.,False],},
                'LSW':{'Volume':[0.05,False],'Mean':[1000.0,False],},
                'Schulz-Zimm':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[300.,False],},
                'Unified':{'G':[1.e3,False],'Rg':[100,False],'B':[1.e-5,False],'P':[4,False],'Cutoff':[1e-5,False],},
                'Porod':{'B':[1.e-4,False],'P':[4,False],'Cutoff':[1e-5,False],},
                'Monodisperse':{'Volume':[0.05,False],'Radius':[100,False],},   #OK for spheres
                'Bragg':{'PkInt':[100,False],'PkPos':[0.2,False],
                    'PkSig':[10,False],'PkGam':[10,False],},        #reasonable 31A peak
                })
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)
                    
        wx.CallAfter(UpdateModelsGrid,G2frame,data)
        
    def OnCopyModel(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy models from\n'+hst[5:]+' to...',
            'Copy models', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections(): 
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()        
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            newdata = copy.deepcopy(data)
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Models'),newdata)
            if newdata['BackFile']:
                Profile = G2frame.PatternTree.GetItemPyData(Id)[1]
                BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,newdata['BackFile'])
                BackSample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,BackId, 'Sample Parameters'))
                Profile[5] = BackSample['Scale'][0]*G2frame.PatternTree.GetItemPyData(BackId)[1][1]
        UpdateModelsGrid(G2frame,newdata)  
        wx.CallAfter(UpdateModelsGrid,G2frame,data)
        RefreshPlots(True)
                
    def OnCopyFlags(event):
        thisModel = copy.deepcopy(data)
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        dlg = G2G.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy sample ref. flags from\n'+str(hst[5:])+' to...',
            'Copy sample flags', histList)
        distChoice = ['LogNormal','Gaussian','LSW','Schulz-Zimm','Bragg','Unified',
            'Porod','Monodisperse',]
        parmOrder = ['Volume','Radius','Mean','StdDev','G','Rg','B','P',
            'Cutoff','PkInt','PkPos','PkSig','PkGam','VolFr','Dist',]
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: 
                    item = histList[i]
                    Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    newModel = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Models'))
                    newModel['Back'][1] = copy.copy(thisModel['Back'][1])
                    for ilev,level in enumerate(newModel['Particle']['Levels']):
                        for form in level:
                            if form in distChoice:
                                thisForm = thisModel['Particle']['Levels'][ilev][form]                               
                                for item in parmOrder:
                                    if item in thisForm:
                                       level[form][item][1] = copy.copy(thisForm[item][1])
                            elif form == 'Controls':
                                thisForm = thisModel['Particle']['Levels'][ilev][form]['SFargs']
                                for item in parmOrder:
                                    if item in thisForm:
                                        level[form]['SFargs'][item][1] = copy.copy(thisForm[item][1])
        finally:
            dlg.Destroy()
                
    def OnFitModelAll(event):
        choices = G2gd.GetPatternTreeDataNames(G2frame,['SASD',])
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame, 'Sequential SASD refinement',
             'Select dataset to include',choices)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(choices[sel])
            Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequential SASD results')
            if Id:
                SeqResult = G2frame.PatternTree.GetItemPyData(Id)
            else:
                SeqResult = {}
                Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Sequential SASD results')
            SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
            SeqResult['histNames'] = names
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        Reverse = False
        CopyForward = False
        choice = ['Reverse sequence','Copy from prev.']
        dlg = wx.MultiChoiceDialog(G2frame.dataFrame,'Sequential controls','Select controls',choice)
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                if sel:
                    CopyForward = True
                else:
                    Reverse = True
        dlg.Destroy()
        dlg = wx.ProgressDialog('SASD Sequential fit','Data set name = '+names[0],len(names), 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        wx.BeginBusyCursor()
        if Reverse:
            names.reverse()
        JModel = None
        try:
            for i,name in enumerate(names):
                print ' Sequential fit for ',name
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    break
                sId =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name)
                if i and CopyForward:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'),JModel)
                IProfDict,IProfile = G2frame.PatternTree.GetItemPyData(sId)[:2]
                IModel = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'))
                ISample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Sample Parameters'))
                ILimits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Limits'))
                IfOK,result,varyList,sig,Rvals,covMatrix,parmDict,Msg = G2sasd.ModelFit(IProfile,IProfDict,ILimits,ISample,IModel)
                JModel = copy.deepcopy(IModel)
                if not IfOK:
                    G2frame.ErrorDialog('Failed sequential refinement for data '+name,
                        ' Msg: '+Msg+'\nYou need to rethink your selection of parameters\n'+    \
                        ' Model restored to previous version for'+name)
                    SeqResult['histNames'] = names[:i]
                    dlg.Destroy()
                    break
                else:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'),copy.deepcopy(IModel))
                
                G2sasd.ModelFxn(IProfile,IProfDict,ILimits,ISample,IModel)
                SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                    'covMatrix':covMatrix,'title':name,'parmDict':parmDict}
            else:
                dlg.Destroy()
                print ' ***** Small angle sequential refinement successful *****'
        finally:
            wx.EndBusyCursor()    
        G2frame.PatternTree.SetItemPyData(Id,SeqResult)
        G2frame.PatternTree.SelectItem(Id)
        
    def OnFitModel(event):
        if data['Current'] == 'Size dist.':
            if not any(Sample['Contrast']):
                G2frame.ErrorDialog('No contrast; your sample is a vacuum!',
                    'You need to define a scattering substance!\n'+    \
                    ' Do Substances and then Sample parameters')
                return
            G2sasd.SizeDistribution(Profile,ProfDict,Limits,Sample,data)
            G2plt.PlotSASDSizeDist(G2frame)
            RefreshPlots(True)
            
        elif data['Current'] == 'Particle fit':
            SaveState()
            Results = G2sasd.ModelFit(Profile,ProfDict,Limits,Sample,data)
            if not Results[0]:
                    G2frame.ErrorDialog('Failed refinement',
                        ' Msg: '+Results[-1]+'\nYou need to rethink your selection of parameters\n'+    \
                        ' Model restored to previous version')
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            
    def OnUnDo(event):
        DoUnDo()
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Models'))
        G2frame.dataFrame.SasdUndo.Enable(False)
        UpdateModelsGrid(G2frame,data)
        G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
        RefreshPlots(True)

    def DoUnDo():
        print 'Undo last refinement'
        file = open(G2frame.undosasd,'rb')
        PatternId = G2frame.PatternId
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Models'),cPickle.load(file))
        print ' Models recovered'
        file.close()
        
    def SaveState():
        G2frame.undosasd = os.path.join(G2frame.dirname,'GSASIIsasd.save')
        file = open(G2frame.undosasd,'wb')
        PatternId = G2frame.PatternId
        for item in ['Models']:
            cPickle.dump(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,item)),file,1)
        file.close()
        G2frame.dataFrame.SasdUndo.Enable(True)
        
    def OnSelectFit(event):
        data['Current'] = fitSel.GetValue()
        wx.CallAfter(UpdateModelsGrid,G2frame,data)
        
    def OnCheckBox(event):
        Obj = event.GetEventObject()
        item,ind = Indx[Obj.GetId()]
        item[ind] = Obj.GetValue()
        
    def OnIntVal(event):
        event.Skip()
        Obj = event.GetEventObject()
        item,ind,minVal = Indx[Obj.GetId()]
        try:
            value = int(Obj.GetValue())
            if value <= minVal:
                raise ValueError
        except ValueError:
            value = item[ind]
        Obj.SetValue(str(value))
        item[ind] = value

    def SizeSizer():
        
        def OnShape(event):
            data['Size']['Shape'][0] = partsh.GetValue()
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            
        def OnMethod(event):
            data['Size']['Method'] = method.GetValue()
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            
        sizeSizer = wx.BoxSizer(wx.VERTICAL)
        sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Size distribution parameters: '),0,WACV)
        binSizer = wx.FlexGridSizer(0,7,5,5)
        binSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' No. size bins: '),0,WACV)
        bins = ['50','100','150','200']
        nbins = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['Nbins']),choices=bins,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[nbins.GetId()] = [data['Size'],'Nbins',0]
        nbins.Bind(wx.EVT_COMBOBOX,OnIntVal)        
        binSizer.Add(nbins,0,WACV)
        binSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min diam.: '),0,WACV)
        minDias = ['10','25','50','100','150','200']
        mindiam = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['MinDiam']),choices=minDias,
            style=wx.CB_DROPDOWN)
        mindiam.Bind(wx.EVT_LEAVE_WINDOW,OnIntVal)
        mindiam.Bind(wx.EVT_TEXT_ENTER,OnIntVal)        
        mindiam.Bind(wx.EVT_KILL_FOCUS,OnIntVal)
        Indx[mindiam.GetId()] = [data['Size'],'MinDiam',0]
        binSizer.Add(mindiam,0,WACV)
        binSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Max diam.: '),0,WACV)
        maxDias = [str(1000*(i+1)) for i in range(10)]
        maxdiam = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['MaxDiam']),choices=maxDias,
            style=wx.CB_DROPDOWN)
        maxdiam.Bind(wx.EVT_LEAVE_WINDOW,OnIntVal)
        maxdiam.Bind(wx.EVT_TEXT_ENTER,OnIntVal)        
        maxdiam.Bind(wx.EVT_KILL_FOCUS,OnIntVal)
        Indx[maxdiam.GetId()] = [data['Size'],'MaxDiam',0]
        binSizer.Add(maxdiam,0,WACV)
        logbins = wx.CheckBox(G2frame.dataDisplay,label='Log bins?')
        Indx[logbins.GetId()] = [data['Size'],'logBins']
        logbins.SetValue(data['Size']['logBins'])
        logbins.Bind(wx.EVT_CHECKBOX, OnCheckBox)
        binSizer.Add(logbins,0,WACV)
        sizeSizer.Add(binSizer,0)
        sizeSizer.Add((5,5),0)
        partSizer = wx.BoxSizer(wx.HORIZONTAL)
        partSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Particle description: '),0,WACV)
        shapes = {'Spheroid':' Aspect ratio: ','Cylinder':' Diameter ','Cylinder AR':' Aspect ratio: ',
            'Unified sphere':'','Unified rod':' Diameter: ','Unified rod AR':' Aspect ratio: ',
            'Unified disk':' Thickness: ', 'Spherical shell': ' Shell thickness'}
        partsh = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['Shape'][0]),choices=shapes.keys(),
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        partsh.Bind(wx.EVT_COMBOBOX,OnShape)        
        partSizer.Add(partsh,0,WACV)
        if data['Size']['Shape'][0] not in ['Unified sphere',]:
            partSizer.Add(wx.StaticText(G2frame.dataDisplay,label=shapes[data['Size']['Shape'][0]]),0,WACV)
            partprm = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Size']['Shape'],1,
                nDig=(10,3),typeHint=float,min=0.)
            partSizer.Add(partprm,0,WACV)
        sizeSizer.Add(partSizer,0)
        sizeSizer.Add((5,5),0)
        fitSizer = wx.BoxSizer(wx.HORIZONTAL)
        methods = ['MaxEnt','IPG',]
        fitSizer.Add(wx.StaticText(G2frame.dataDisplay,label='Fitting method: '),0,WACV)
        method = wx.ComboBox(G2frame.dataDisplay,value=data['Size']['Method'],choices=methods,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        method.Bind(wx.EVT_COMBOBOX,OnMethod)
        fitSizer.Add(method,0,WACV)
        iters = ['10','25','50','100','150','200']        
        fitSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' No. iterations: '),0,WACV)
        Method = data['Size']['Method']
        iter = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size'][Method]['Niter']),choices=iters,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[iter.GetId()] = [data['Size'][Method],'Niter',0]
        iter.Bind(wx.EVT_COMBOBOX,OnIntVal)
        fitSizer.Add(iter,0,WACV)
        if 'MaxEnt' in data['Size']['Method']:
            fitSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Log floor factor: '),0,WACV)
            floors = [str(-i) for i in range(9)]
            floor = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['MaxEnt']['Sky']),choices=floors,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[floor.GetId()] = [data['Size']['MaxEnt'],'Sky',-10]
            floor.Bind(wx.EVT_COMBOBOX,OnIntVal)
            fitSizer.Add(floor,0,WACV)
        elif 'IPG' in data['Size']['Method']:
            fitSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Q power weight (-1 for sigma): '),0,WACV)
            choices = ['-1','0','1','2','3','4']
            power = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['IPG']['Power']),choices=choices,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[power.GetId()] = [data['Size']['IPG'],'Power',-2]
            power.Bind(wx.EVT_COMBOBOX,OnIntVal)
            fitSizer.Add(power,0,WACV)
        sizeSizer.Add(fitSizer,0)

        return sizeSizer
        
    def PartSizer():
        
        FormFactors = {'Sphere':{},'Spheroid':{'Aspect ratio':[1.0,False]},
            'Cylinder':{'Length':[100.,False]},'Cylinder diam':{'Diameter':[100.,False]},
            'Cylinder AR':{'Aspect ratio':[1.0,False]},'Unified sphere':{},
            'Unified rod':{'Length':[100.,False]},'Unified rod AR':{'Aspect ratio':[1.0,False]},
            'Unified disk':{'Thickness':[100.,False]},
            'Unified tube':{'Length':[100.,False],'Thickness':[10.,False]},
            'Spherical shell':{'Shell thickness':[1.5,False] }, }
                
        StructureFactors = {'Dilute':{},'Hard sphere':{'VolFr':[0.1,False],'Dist':[100.,False]},
            'Sticky hard sphere':{'VolFr':[0.1,False],'Dist':[100.,False],'epis':[0.05,False],'Sticky':[0.2,False]},
            'Square well':{'VolFr':[0.1,False],'Dist':[100.,False],'Depth':[0.1,False],'Width':[1.,False]},
            'InterPrecipitate':{'VolFr':[0.1,False],'Dist':[100.,False]},}
                
        ffDistChoices =  ['Sphere','Spheroid','Cylinder','Cylinder diam',
            'Cylinder AR','Unified sphere','Unified rod','Unified rod AR',
            'Unified disk','Unified tube','Spherical shell',]
                
        ffMonoChoices = ['Sphere','Spheroid','Cylinder','Cylinder AR',]
        
        sfChoices = ['Dilute','Hard sphere','Sticky hard sphere','Square well','InterPrecipitate',]
            
        slMult = 1000.
                  
        def OnValue(event):
            event.Skip()
            Obj = event.GetEventObject()
            item,key,sldrObj = Indx[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
                if value <= 0.:
                    raise ValueError
            except ValueError:
                value = item[key][0]
            item[key][0] = value
            Obj.SetValue('%.3g'%(value))
            if key in ['P','epis','Sticky','Depth','Width','VolFr','Dist']:
                sldrObj.SetValue(slMult*value)
            else:
                logv = np.log10(value)
                valMinMax = [logv-1,logv+1]
                sldrObj.SetRange(slMult*valMinMax[0],slMult*valMinMax[1])
                sldrObj.SetValue(slMult*logv)
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots()
            
        def OnSelect(event):
            Obj = event.GetEventObject()
            item,key = Indx[Obj.GetId()]
            item[key] = Obj.GetValue()
            if 'Refine' not in Obj.GetLabel():
                if 'FormFact' in key :
                    item['FFargs'] = FormFactors[Obj.GetValue()]
                elif 'StrFact' in key:
                    item['SFargs'] = StructureFactors[Obj.GetValue()]
                wx.CallAfter(UpdateModelsGrid,G2frame,data)
                G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
                RefreshPlots()
                
        def OnDelLevel(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            del data['Particle']['Levels'][item]
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots()
            
        def OnParmSlider(event):
            Obj = event.GetEventObject()
            item,key,pvObj = Indx[Obj.GetId()]
            slide = Obj.GetValue()
            if key in ['P','epis','Sticky','Depth','Width','VolFr','Dist']:
                value = float(slide/slMult)
            else:
                value = 10.**float(slide/slMult)
            item[key][0] = value
            pvObj.SetValue('%.3g'%(item[key][0]))
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots()
            
        def SizeSizer():
            sizeSizer = wx.FlexGridSizer(0,4,5,5)
            sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Distribution: '),0,WACV)
            Distchoice = ['LogNormal','Gaussian','LSW','Schulz-Zimm','Bragg','Unified','Porod','Monodisperse',]
            distChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['DistType'],choices=Distchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[distChoice.GetId()] = [level['Controls'],'DistType']
            distChoice.Bind(wx.EVT_COMBOBOX,OnSelect)
            sizeSizer.Add(distChoice,0,WACV)    #put structure factor choices here
            if level['Controls']['DistType'] not in ['Bragg','Unified','Porod',]:
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Form Factor: '),0,WACV)
                if 'Mono' in level['Controls']['DistType']:
                    ffChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['FormFact'],choices=ffMonoChoices,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                else:
                    ffChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['FormFact'],choices=ffDistChoices,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[ffChoice.GetId()] = [level['Controls'],'FormFact']
                ffChoice.Bind(wx.EVT_COMBOBOX,OnSelect)
                sizeSizer.Add(ffChoice,0,WACV)
                
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Material: '),0,WACV)
                matSel = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['Material'],
                    choices=Substances['Substances'].keys(),style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[matSel.GetId()] = [level['Controls'],'Material']
                matSel.Bind(wx.EVT_COMBOBOX,OnSelect)        
                sizeSizer.Add(matSel,0,WACV) #do neutron test here?
                rho = Substances['Substances'][level['Controls']['Material']].get('XAnom density',0.0)
                level['Controls']['Contrast'] = contrast = (rho-rhoMat)**2                 
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Resonant X-ray contrast: '),0,WACV)
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label='  %.2f 10%scm%s'%(contrast,Pwr20,Pwrm4)),0,WACV)
                if 'Mono' not in level['Controls']['DistType']:
                    sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Num. radii: '),0,WACV)
                    radii = ['25','50','75','100','200']
                    nRadii = wx.ComboBox(G2frame.dataDisplay,value=str(level['Controls']['NumPoints']),choices=radii,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    Indx[nRadii.GetId()] = [level['Controls'],'NumPoints']
                    nRadii.Bind(wx.EVT_COMBOBOX,OnSelect)
                    sizeSizer.Add(nRadii,0,WACV)
                    sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' R dist. cutoff: '),0,WACV)
                    rCutoff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,level['Controls'],'Cutoff',
                        min=0.001,max=0.1,typeHint=float)
                    sizeSizer.Add(rCutoff,0,WACV)
            elif level['Controls']['DistType']  in ['Unified',]:
                Parms = level['Unified']
                Best = G2sasd.Bestimate(Parms['G'][0],Parms['Rg'][0],Parms['P'][0])
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Estimated Dist B: %12.4g'%(Best)),0,WACV)
            return sizeSizer
            
        def ParmSizer():
            parmSizer = wx.FlexGridSizer(0,3,5,5)
            parmSizer.AddGrowableCol(2,1)
            parmSizer.SetFlexibleDirection(wx.HORIZONTAL)
            Parms = level[level['Controls']['DistType']]
            FFargs = level['Controls']['FFargs']
            SFargs = level['Controls'].get('SFargs',{})
            parmOrder = ['Volume','Radius','Mean','StdDev','MinSize','G','Rg','B','P','Cutoff',
                'PkInt','PkPos','PkSig','PkGam',]
            for parm in parmOrder:
                if parm in Parms:
                    if parm == 'MinSize':
                        parmSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Dist '+parm),0,wx.ALIGN_CENTER)
                    else:
                        parmVar = wx.CheckBox(G2frame.dataDisplay,label='Refine? Dist '+parm) 
                        parmVar.SetValue(Parms[parm][1])
                        parmVar.Bind(wx.EVT_CHECKBOX, OnSelect)
                        parmSizer.Add(parmVar,0,WACV)
                        Indx[parmVar.GetId()] = [Parms[parm],1]
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
                    parmValue = wx.TextCtrl(G2frame.dataDisplay,value='%.3g'%(Parms[parm][0]),
                        style=wx.TE_PROCESS_ENTER)
                    parmValue.Bind(wx.EVT_TEXT_ENTER,OnValue)        
                    parmValue.Bind(wx.EVT_KILL_FOCUS,OnValue)
                    parmSizer.Add(parmValue,0,WACV)
                    if parm == 'P':
                        value = Parms[parm][0]
                        valMinMax = [0.1,4.2]
                    else:
                        value = np.log10(Parms[parm][0])
                        valMinMax = [value-1,value+1]
                    parmSldr = wx.Slider(G2frame.dataDisplay,minValue=slMult*valMinMax[0],
                        maxValue=slMult*valMinMax[1],value=slMult*value)
                    Indx[parmValue.GetId()] = [Parms,parm,parmSldr]
                    Indx[parmSldr.GetId()] = [Parms,parm,parmValue]
                    parmSldr.Bind(wx.EVT_SLIDER,OnParmSlider)
                    parmSizer.Add(parmSldr,1,wx.EXPAND)
            if level['Controls']['DistType'] not in ['Bragg']:
                parmOrder = ['Aspect ratio','Length','Diameter','Thickness','VolFr','Dist','epis','Sticky','Depth','Width','Shell thickness',]
                fTypes = ['FF ','SF ']
                for iarg,Args in enumerate([FFargs,SFargs]):
                    for parm in parmOrder:
                        if parm in Args:
                            parmVar = wx.CheckBox(G2frame.dataDisplay,label='Refine? '+fTypes[iarg]+parm) 
                            parmVar.SetValue(Args[parm][1])
                            Indx[parmVar.GetId()] = [Args[parm],1]
                            parmVar.Bind(wx.EVT_CHECKBOX, OnSelect)
                            parmSizer.Add(parmVar,0,WACV)
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
                            parmValue = wx.TextCtrl(G2frame.dataDisplay,value='%.3g'%(Args[parm][0]),
                                style=wx.TE_PROCESS_ENTER)
                            parmValue.Bind(wx.EVT_TEXT_ENTER,OnValue)        
                            parmValue.Bind(wx.EVT_KILL_FOCUS,OnValue)
                            parmSizer.Add(parmValue,0,WACV)
                            value = Args[parm][0]
                            if parm == 'epis':
                                valMinMax = [0,.1]
                            elif parm in ['Sticky','Width',]:
                                valMinMax = [0,1.]
                            elif parm == 'Depth':
                                valMinMax = [-2.,2.]
                            elif parm == 'Dist':
                                valMinMax = [100.,1000.]
                            elif parm == 'VolFr':
                                valMinMax = [1.e-4,1.]
                            else:
                                value = np.log10(Args[parm][0])
                                valMinMax = [value-1,value+1]
                            parmSldr = wx.Slider(G2frame.dataDisplay,minValue=slMult*valMinMax[0],
                                maxValue=slMult*valMinMax[1],value=slMult*value)
                            Indx[parmVar.GetId()] = [Args[parm],1]
                            Indx[parmValue.GetId()] = [Args,parm,parmSldr]
                            Indx[parmSldr.GetId()] = [Args,parm,parmValue]
                            parmSldr.Bind(wx.EVT_SLIDER,OnParmSlider)
                            parmSizer.Add(parmSldr,1,wx.EXPAND)
            return parmSizer                
            
        Indx = {}
        partSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Particle fit parameters: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Matrix: '),0,WACV)
        matsel = wx.ComboBox(G2frame.dataDisplay,value=data['Particle']['Matrix']['Name'],
            choices=Substances['Substances'].keys(),style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[matsel.GetId()] = [data['Particle']['Matrix'],'Name'] 
        matsel.Bind(wx.EVT_COMBOBOX,OnSelect) #Do neutron test here?
        rhoMat = Substances['Substances'][data['Particle']['Matrix']['Name']].get('XAnom density',0.0)        
        topSizer.Add(matsel,0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Volume fraction: '),0,WACV)
        volfrac = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Particle']['Matrix']['VolFrac'],0,
                typeHint=float)
        topSizer.Add(volfrac,0,WACV)
        volVar = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
        volVar.SetValue(data['Particle']['Matrix']['VolFrac'][1])
        Indx[volVar.GetId()] = [data['Particle']['Matrix']['VolFrac'],1]
        volVar.Bind(wx.EVT_CHECKBOX, OnSelect)
        topSizer.Add(volVar,0,WACV)
        partSizer.Add(topSizer,0,)
        for ilev,level in enumerate(data['Particle']['Levels']):
            G2G.HorizontalLine(partSizer,G2frame.dataDisplay)
            topLevel = wx.BoxSizer(wx.HORIZONTAL)
            topLevel.Add(wx.StaticText(G2frame.dataDisplay,label=' Model component %d: '%(ilev)),0,WACV)
            delBtn = wx.Button(G2frame.dataDisplay,label=' Delete?')
            Indx[delBtn.GetId()] = ilev
            delBtn.Bind(wx.EVT_BUTTON,OnDelLevel)
            topLevel.Add(delBtn,0,WACV)
            partSizer.Add(topLevel,0)
            partSizer.Add(SizeSizer())
            if level['Controls']['DistType'] not in ['Bragg','Unified','Porod',]:
                topLevel.Add(wx.StaticText(G2frame.dataDisplay,label=' Structure factor: '),0,WACV)
                strfctr = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['StrFact'],
                    choices=sfChoices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[strfctr.GetId()] = [level['Controls'],'StrFact']
                strfctr.Bind(wx.EVT_COMBOBOX,OnSelect)
                topLevel.Add(strfctr,0,WACV)
            partSizer.Add(ParmSizer(),0,wx.EXPAND)
        return partSizer
        
    def OnEsdScale(event):
        event.Skip()
        try:
            value = float(esdScale.GetValue())
            if value <= 0.:
                raise ValueError
        except ValueError:
            value = 1./np.sqrt(ProfDict['wtFactor'])
        ProfDict['wtFactor'] = 1./value**2
        esdScale.SetValue('%.3f'%(value))
        RefreshPlots(True)
        
    def OnBackChange(invalid,value,tc):
        Profile[4][:] = value
        RefreshPlots()
        
    def OnBackFile(event):  #multiple backgrounds?
        data['BackFile'] = backFile.GetValue()
        if data['BackFile']:
            BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['BackFile'])
            BackSample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,BackId, 'Sample Parameters'))
            Profile[5] = BackSample['Scale'][0]*G2frame.PatternTree.GetItemPyData(BackId)[1][1]
        else:
            Profile[5] = np.zeros(len(Profile[5]))
        RefreshPlots(True)
            
    Sample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Sample Parameters'))
    Limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits'))
    Substances = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Substances'))
    ProfDict,Profile = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[:2]
    if data['BackFile']:
        BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['BackFile'])
        BackSample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,BackId, 'Sample Parameters'))
        Profile[5] = BackSample['Scale'][0]*G2frame.PatternTree.GetItemPyData(BackId)[1][1]
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()   # is this a ScrolledWindow? If so, bad!
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ModelMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.SetLabel('Modelling')
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyModel, id=G2gd.wxID_MODELCOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyFlags, id=G2gd.wxID_MODELCOPYFLAGS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitModel, id=G2gd.wxID_MODELFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitModelAll, id=G2gd.wxID_MODELFITALL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_MODELUNDO)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddModel, id=G2gd.wxID_MODELADD)
    Indx = {}
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    models = ['Size dist.','Particle fit']
    topSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Modeling by: '),0,WACV)
    fitSel = wx.ComboBox(G2frame.dataDisplay,value=data['Current'],choices=models,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    fitSel.Bind(wx.EVT_COMBOBOX,OnSelectFit)        
    topSizer.Add(fitSel,0,WACV)
    topSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Error multiplier: '),0,WACV)
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
    esdScale = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(1./np.sqrt(ProfDict['wtFactor'])),style=wx.TE_PROCESS_ENTER)
    esdScale.Bind(wx.EVT_TEXT_ENTER,OnEsdScale)        
    esdScale.Bind(wx.EVT_KILL_FOCUS,OnEsdScale)
    topSizer.Add(esdScale,0,WACV)
    mainSizer.Add(topSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)
    if 'Size' in data['Current']:
        if 'MaxEnt' in data['Size']['Method']:
            Status.SetStatusText('Size distribution by Maximum entropy')
        elif 'IPG' in data['Size']['Method']:
            Status.SetStatusText('Size distribution by Interior-Point Gradient')
        mainSizer.Add(SizeSizer())        
    elif 'Particle' in data['Current']:
        mainSizer.Add(PartSizer(),1,wx.ALIGN_LEFT|wx.EXPAND)
    G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)    
    backSizer = wx.BoxSizer(wx.HORIZONTAL)
    backSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Background:'),0,WACV)
    backVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Back'],0,
        nDig=(10,3),typeHint=float,OnLeave=OnBackChange)
    backSizer.Add(backVal,0,WACV)
    backVar = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
    Indx[backVar.GetId()] = [data['Back'],1]
    backVar.SetValue(data['Back'][1])
    backVar.Bind(wx.EVT_CHECKBOX, OnCheckBox)
    backSizer.Add(backVar,0,WACV)
    #multiple background files?
    backSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Background file: '),0,WACV)
    Choices = ['',]+G2gd.GetPatternTreeDataNames(G2frame,['SASD',])
    backFile = wx.ComboBox(parent=G2frame.dataDisplay,value=data['BackFile'],choices=Choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    backFile.Bind(wx.EVT_COMBOBOX,OnBackFile)
    backSizer.Add(backFile)    
    mainSizer.Add(backSizer)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    G2frame.dataFrame.setSizePosLeft(Size)    
    
################################################################################
#####  REFD Models 
################################################################################           
       
def UpdateREFDModelsGrid(G2frame,data):
    '''respond to selection of REFD Models data tree item.
    '''
    def OnCopyModel(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame,'Copy reflectivity models from\n'+str(hst[5:])+' to...',
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternTree.SetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,Id,'Models'),copy.deepcopy(data))
        
    def OnFitModel(event):
        
        SaveState()
        G2pwd.REFDRefine(Profile,ProfDict,Inst,Limits,Substances,data)
        x,xr,y = G2pwd.makeSLDprofile(data,Substances)
        ModelPlot(data,x,xr,y)
        G2plt.PlotPatterns(G2frame,plotType='REFD')
        wx.CallLater(100,UpdateREFDModelsGrid,G2frame,data)
        
    def OnModelPlot(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetFileList(G2frame,'REFD')
#        histList = [hst,]
#        histList += GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        plotList = []
        od = {'label_1':'Zero at substrate','value_1':False,'label_2':'Show layer transitions','value_2':True}
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame,'Plot reflectivity models for:',
            'Plot SLD models', histList,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    plotList.append(histList[i])
            else:
                dlg.Destroy()
                return
        finally:
            dlg.Destroy()
        XY = []
        LinePos = []
        for item in plotList:
            mId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
            model = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,mId,'Models'))
            Substances = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,mId,'Substances'))['Substances']       
            x,xr,y = G2pwd.makeSLDprofile(model,Substances)
            if od['value_1']:
                XY.append([xr,y])
                disLabel = r'$Distance\ from\ substrate,\ \AA$'
            else:
                XY.append([x,y])
                disLabel = r'$Distance\ from\ top\ surface,\ \AA$'
            if od['value_2']:
                laySeq = model['Layer Seq'].split()
                nLines = len(laySeq)+1
                linePos = np.zeros(nLines)
                for ilay,lay in enumerate(np.fromstring(data['Layer Seq'],dtype=int,sep=' ')):
                    linePos[ilay+1:] += model['Layers'][lay].get('Thick',[0.,False])[0]
                if od['value_1']:
                    linePos = linePos[-1]-linePos
                LinePos.append(linePos)
        G2plt.PlotXY(G2frame,XY,labelX=disLabel,labelY=r'$SLD,\ 10^{10}cm^{-2}$',newPlot=True,   
                      Title='Scattering length density',lines=True,names=[],vertLines=LinePos)
        
    def OnFitModelAll(event):
        choices = G2gd.GetPatternTreeDataNames(G2frame,['REFD',])
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame, 'Sequential REFD refinement',
             'Select dataset to include',choices)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(choices[sel])
            Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequential REFD results')
            if Id:
                SeqResult = G2frame.PatternTree.GetItemPyData(Id)
            else:
                SeqResult = {}
                Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Sequential REFD results')
            SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
            SeqResult['histNames'] = names
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        Reverse = False
        CopyForward = False
        choice = ['Reverse sequence','Copy from prev.']
        dlg = wx.MultiChoiceDialog(G2frame.dataFrame,'Sequential controls','Select controls',choice)
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                if sel:
                    CopyForward = True
                else:
                    Reverse = True
        dlg.Destroy()
        dlg = wx.ProgressDialog('REFD Sequential fit','Data set name = '+names[0],len(names), 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        wx.BeginBusyCursor()
        if Reverse:
            names.reverse()
        JModel = None
        try:
            for i,name in enumerate(names):
                print ' Sequential fit for ',name
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    break
                sId =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name)
                if i and CopyForward:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'),JModel)
                IProfDict,IProfile = G2frame.PatternTree.GetItemPyData(sId)[:2]
                IModel = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'))
                ISubstances = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Substances'))['Substances']
                ILimits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Limits'))
                IfOK,result,varyList,sig,Rvals,covMatrix,parmDict,Msg = G2pwd.REFDRefine(IProfile,IProfDict,Inst,ILimits,ISubstances,IModel)
                JModel = copy.deepcopy(IModel)
                if not IfOK:
                    G2frame.ErrorDialog('Failed sequential refinement for data '+name,
                        ' Msg: '+Msg+'\nYou need to rethink your selection of parameters\n'+    \
                        ' Model restored to previous version for'+name)
                    SeqResult['histNames'] = names[:i]
                    dlg.Destroy()
                    break
                else:
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,sId, 'Models'),copy.deepcopy(IModel))
                
                SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                    'covMatrix':covMatrix,'title':name,'parmDict':parmDict}
            else:
                dlg.Destroy()
                print ' ***** Small angle sequential refinement successful *****'
        finally:
            wx.EndBusyCursor()    
        G2frame.PatternTree.SetItemPyData(Id,SeqResult)
        G2frame.PatternTree.SelectItem(Id)
                
    def ModelPlot(data,x,xr,y):
        laySeq = data['Layer Seq'].split()
        nLines = len(laySeq)+1
        linePos = np.zeros(nLines)
        for ilay,lay in enumerate(np.fromstring(data['Layer Seq'],dtype=int,sep=' ')):
            linePos[ilay+1:] += data['Layers'][lay].get('Thick',[0.,False])[0]
        if data['Zero'] == 'Top':
            XY = [[x,y],]
            disLabel = r'$Distance\ from\ top\ surface,\ \AA$'
        else:
            XY = [[xr,y],]
            linePos = linePos[-1]-linePos
            disLabel = r'$Distance\ from\ substrate,\ \AA$'
        G2plt.PlotXY(G2frame,XY,labelX=disLabel,labelY=r'$SLD,\ 10^{10}cm^{-2}$',newPlot=True,
            Title='Scattering length density',lines=True,names=[],vertLines=[linePos,])
        
    def OnUnDo(event):
        DoUnDo()
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Models'))
        G2frame.dataFrame.REFDUndo.Enable(False)
        G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
        x,xr,y = G2pwd.makeSLDprofile(data,Substances)
        ModelPlot(data,x,xr,y)
        G2plt.PlotPatterns(G2frame,plotType='REFD')
        wx.CallLater(100,UpdateREFDModelsGrid,G2frame,data)

    def DoUnDo():
        print 'Undo last refinement'
        file = open(G2frame.undorefd,'rb')
        PatternId = G2frame.PatternId
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Models'),cPickle.load(file))
        print ' Model recovered'
        file.close()
        
    def SaveState():
        G2frame.undorefd = os.path.join(G2frame.dirname,'GSASIIrefd.save')
        file = open(G2frame.undorefd,'wb')
        PatternId = G2frame.PatternId
        cPickle.dump(G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Models')),file,1)
        file.close()
        G2frame.dataFrame.REFDUndo.Enable(True)
    
    def ControlSizer():
        
        def OnRefPos(event):
            data['Zero'] = refpos.GetValue()
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            
        def OnMinSel(event):
            data['Minimizer'] = minSel.GetValue()
            
        def OnWeight(event):
            data['2% weight'] = weight.GetValue()
            
        def OnSLDplot(event):
            sld.SetValue(False)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            
        def OnQ4fftplot(event):
            q4fft.SetValue(False)
            R,F = G2pwd.makeRefdFFT(Limits,Profile)
            XY = [[R[:2500],F[:2500]],]
            G2plt.PlotXY(G2frame,XY,labelX='thickness',labelY='F(R)',newPlot=True,
                Title='Fourier transform',lines=True)
            
        def OndQSel(event):
            data['dQ type'] = dQSel.GetStringSelection()
            
        controlSizer = wx.BoxSizer(wx.VERTICAL)
        resol = wx.BoxSizer(wx.HORIZONTAL)
        choice = ['None','const '+GkDelta+'Q/Q',]
        if ProfDict['ifDQ']:
            choice += [GkDelta+'Q/Q in data']
        dQSel = wx.RadioBox(G2frame.dataDisplay,wx.ID_ANY,'Instrument resolution type:',choices=choice,
            majorDimension=0,style=wx.RA_SPECIFY_COLS)
        dQSel.SetStringSelection(data['dQ type'])
        dQSel.Bind(wx.EVT_RADIOBOX,OndQSel)
        resol.Add(dQSel,0,WACV)
        resol.Add(wx.StaticText(G2frame.dataDisplay,label=' (FWHM %): '),0,WACV)
        resol.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Resolution'],0,nDig=(10,3),min=0.,max=5.),0,WACV)
        controlSizer.Add(resol,0,WACV)
        minimiz = wx.BoxSizer(wx.HORIZONTAL)
        minimiz.Add(wx.StaticText(G2frame.dataDisplay,label=' Minimizer: '),0,WACV)
        minlist = ['LMLS','Basin Hopping','MC/SA Anneal','L-BFGS-B',]
        minSel = wx.ComboBox(G2frame.dataDisplay,value=data['Minimizer'],choices=minlist,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        minSel.Bind(wx.EVT_COMBOBOX, OnMinSel)
        minimiz.Add(minSel,0,WACV)
        minimiz.Add(wx.StaticText(G2frame.dataDisplay,label=' Bounds factor: '),0,WACV)
        minimiz.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Toler',nDig=(10,2),max=0.99,min=0.1),0,WACV)
        weight = wx.CheckBox(G2frame.dataDisplay,label='Use 2% sig. weights')
        weight.SetValue(data.get('2% weight',False))
        weight.Bind(wx.EVT_CHECKBOX, OnWeight)
        minimiz.Add(weight,0,WACV)
        controlSizer.Add(minimiz,0,WACV)
        plotSizer = wx.BoxSizer(wx.HORIZONTAL)
        plotSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Plot controls: '),0,WACV)
        sld = wx.CheckBox(G2frame.dataDisplay,label='Plot SLD?')
        sld.Bind(wx.EVT_CHECKBOX, OnSLDplot)
        plotSizer.Add(sld,0,WACV)
        plotSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Zero position location: '),0,WACV)
        poslist = ['Top','Bottom']
        refpos = wx.ComboBox(G2frame.dataDisplay,value=data['Zero'],choices=poslist,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        refpos.Bind(wx.EVT_COMBOBOX, OnRefPos)
        plotSizer.Add(refpos,0,WACV)
#        q4fft = wx.CheckBox(G2frame.dataDisplay,label='Plot fft?')
#        q4fft.Bind(wx.EVT_CHECKBOX, OnQ4fftplot)
#        plotSizer.Add(q4fft,0,WACV)
        controlSizer.Add(plotSizer,0,WACV)
        return controlSizer
    
    def OverallSizer():
#'DualFitFile':'', 'DualFltBack':[0.0,False],'DualScale':[1.0,False] future for neutrons - more than one?
        
        def OnScaleRef(event):
            data['Scale'][1] = scaleref.GetValue()
            
        def OnBackRef(event):
            data['FltBack'][1] = backref.GetValue()
            
        def Recalculate(invalid,value,tc):
            if invalid:
                return
            
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2plt.PlotPatterns(G2frame,plotType='REFD')

        overall = wx.BoxSizer(wx.HORIZONTAL)
        overall.Add(wx.StaticText(G2frame.dataDisplay,label=' Scale: '),0,WACV)
        overall.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Scale'],0,
            nDig=(10,2),typeHint=float,OnLeave=Recalculate),0,WACV)
        scaleref = wx.CheckBox(G2frame.dataDisplay,label=' Refine?  ')
        scaleref.SetValue(data['Scale'][1])
        scaleref.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        overall.Add(scaleref,0,WACV)
        overall.Add(wx.StaticText(G2frame.dataDisplay,label=' Flat bkg.: '),0,WACV)
        overall.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['FltBack'],0,
            nDig=(10,2,'g'),typeHint=float,OnLeave=Recalculate),0,WACV)
        backref = wx.CheckBox(G2frame.dataDisplay,label=' Refine?  ')
        backref.SetValue(data['FltBack'][1])
        backref.Bind(wx.EVT_CHECKBOX, OnBackRef)
        overall.Add(backref,0,WACV)        
        return overall
        
    def LayerSizer():
    #'Penetration':[0.,False]?

        def OnSelect(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            Name = Obj.GetValue()
            data['Layers'][item]['Name'] = Name
            if 'Rough' not in data['Layers'][item]:
                data['Layers'][item]['Rough'] = [0.,False]
            if 'Thick' not in data['Layers'][item]:
                data['Layers'][item]['Thick'] = [10.,False]
            if 'N' in Inst['Type'][0]:
                data['Layers'][item]['Mag SLD'] = [0.,False]
            if Name == 'unit scatter':
                data['Layers'][item]['iDenMul'] = [0.,False]
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            G2plt.PlotPatterns(G2frame,plotType='REFD')
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)
            
        def OnCheckBox(event):
            Obj = event.GetEventObject()
            item,parm = Indx[Obj.GetId()]
            data['Layers'][item][parm][1] = Obj.GetValue()
            
        def OnInsertLayer(event):
            Obj = event.GetEventObject()
            ind = Indx[Obj.GetId()]
            data['Layers'].insert(ind+1,{'Name':'vacuum','DenMul':[1.0,False],})
            data['Layer Seq'] = ' '.join([str(i+1) for i in range(len(data['Layers'])-2)])
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            G2plt.PlotPatterns(G2frame,plotType='REFD')
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)
            
        def OnDeleteLayer(event):
            Obj = event.GetEventObject()
            ind = Indx[Obj.GetId()]
            del data['Layers'][ind]
            data['Layer Seq'] = ' '.join([str(i+1) for i in range(len(data['Layers'])-2)])
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            G2plt.PlotPatterns(G2frame,plotType='REFD')
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data) 

        def Recalculate(invalid,value,tc):
            if invalid:
                return
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2plt.PlotPatterns(G2frame,plotType='REFD')
            wx.CallLater(100,UpdateREFDModelsGrid,G2frame,data) 

        Indx = {}                       
        layerSizer = wx.BoxSizer(wx.VERTICAL)
        
        for ilay,layer in enumerate(data['Layers']):
            if not ilay:
                layerSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Top layer (superphase):'),0,WACV)
            elif ilay < len(data['Layers'])-1:
                layerSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Layer no. %d'%(ilay)),0,WACV)
            else:
                layerSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Bottom layer (substrate):'),0,WACV)
            midlayer = wx.BoxSizer(wx.HORIZONTAL)
            midlayer.Add(wx.StaticText(G2frame.dataDisplay,label=' Substance: '),0,WACV)
            midName = data['Layers'][ilay]['Name']
            midSel = wx.ComboBox(G2frame.dataDisplay,value=midName,
                choices=Substances.keys(),style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[midSel.GetId()] = ilay
            midSel.Bind(wx.EVT_COMBOBOX,OnSelect)
            midlayer.Add(midSel,0,WACV)
            if midName != 'vacuum':
                if midName != 'unit scatter':
                    midlayer.Add(wx.StaticText(G2frame.dataDisplay,label=' Den. Mult.: '),0,WACV)
                    midlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Layers'][ilay]['DenMul'],0,
                        nDig=(10,4),OnLeave=Recalculate),0,WACV)
                    varBox = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
                    Indx[varBox.GetId()] = [ilay,'DenMul']
                    varBox.SetValue(data['Layers'][ilay]['DenMul'][1])
                    varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                    midlayer.Add(varBox,0,WACV)
                    realScatt = data['Layers'][ilay]['DenMul'][0]*Substances[midName]['Scatt density']
                    midlayer.Add(wx.StaticText(G2frame.dataDisplay,
                        label=' Real scat. den.: %.4g'%(realScatt)),0,WACV)
                    imagScatt = data['Layers'][ilay]['DenMul'][0]*Substances[midName]['XImag density']
                    midlayer.Add(wx.StaticText(G2frame.dataDisplay,
                        label=' Imag scat. den.: %.4g'%(imagScatt)),0,WACV)
                else:
                    realScatt = data['Layers'][ilay]['DenMul'][0]
                    midlayer.Add(wx.StaticText(G2frame.dataDisplay,
                        label=' Real scat. den.: %.4g'%(realScatt)),0,WACV)
                    imagScatt = data['Layers'][ilay]['iDenMul'][0]
                    midlayer.Add(wx.StaticText(G2frame.dataDisplay,
                        label=' Imag scat. den.: %.4g'%(imagScatt)),0,WACV)                    
            else:
                midlayer.Add(wx.StaticText(G2frame.dataDisplay,label=', air or gas'),0,WACV)
            layerSizer.Add(midlayer)
            if midName == 'unit scatter':
                nxtlayer = wx.BoxSizer(wx.HORIZONTAL)
                nxtlayer.Add(wx.StaticText(G2frame.dataDisplay,label=' Real Den. : '),0,WACV)                
                nxtlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Layers'][ilay]['DenMul'],0,
                    nDig=(10,4),OnLeave=Recalculate),0,WACV)
                varBox = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
                Indx[varBox.GetId()] = [ilay,'DenMul']
                varBox.SetValue(data['Layers'][ilay]['DenMul'][1])
                varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                nxtlayer.Add(varBox,0,WACV)
                nxtlayer.Add(wx.StaticText(G2frame.dataDisplay,label=' Imag Den. : '),0,WACV)                
                nxtlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Layers'][ilay]['iDenMul'],0,
                    nDig=(10,4),OnLeave=Recalculate),0,WACV)
                varBox = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
                Indx[varBox.GetId()] = [ilay,'iDenMul']
                varBox.SetValue(data['Layers'][ilay]['iDenMul'][1])
                varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                nxtlayer.Add(varBox,0,WACV)
                layerSizer.Add(nxtlayer)
            if midName != 'vacuum':
                if 'N' in Inst['Type'][0] and midName not in ['vacuum','unit scatter']:
                    magLayer = wx.BoxSizer(wx.HORIZONTAL)
                    magLayer.Add(wx.StaticText(G2frame.dataDisplay,label=' Magnetic SLD: '),0,WACV)
                    magLayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Layers'][ilay]['Mag SLD'],0,
                        nDig=(10,4),OnLeave=Recalculate),0,WACV)
                    varBox = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
                    Indx[varBox.GetId()] = [ilay,'Mag SLD']
                    varBox.SetValue(data['Layers'][ilay]['Mag SLD'][1])
                    varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                    magLayer.Add(varBox,0,WACV)
                    magLayer.Add(wx.StaticText(G2frame.dataDisplay,
                        label=' Real+mag scat. den.: %.4g'%(realScatt+data['Layers'][ilay]['Mag SLD'][0])),0,WACV)
                    layerSizer.Add(magLayer)
                if ilay:
                    names = {'Rough':'Upper surface Roughness, '+Angstr,'Thick':'Layer Thickness, '+Angstr}
                    parmsline = wx.BoxSizer(wx.HORIZONTAL)
                    parms= ['Rough','Thick']
                    if ilay == len(data['Layers'])-1:
                        parms = ['Rough',]
                    for parm in parms:
                        parmsline.Add(wx.StaticText(G2frame.dataDisplay,label=' %s: '%(names[parm])),0,WACV)
                        parmsline.Add(G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['Layers'][ilay][parm],0,
                            nDig=(10,2),OnLeave=Recalculate),0,WACV)
                        varBox = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
                        Indx[varBox.GetId()] = [ilay,parm]
                        varBox.SetValue(data['Layers'][ilay][parm][1])
                        varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                        parmsline.Add(varBox,0,WACV)
                    layerSizer.Add(parmsline)
            if ilay < len(data['Layers'])-1:
                newlayer = wx.BoxSizer(wx.HORIZONTAL)
                insert = wx.Button(G2frame.dataDisplay,label='Insert')
                Indx[insert.GetId()] = ilay
                insert.Bind(wx.EVT_BUTTON,OnInsertLayer)
                newlayer.Add(insert)
                delet = wx.Button(G2frame.dataDisplay,label='Delete')
                Indx[delet.GetId()] = ilay
                delet.Bind(wx.EVT_BUTTON,OnDeleteLayer)
                newlayer.Add(delet)
                layerSizer.Add(newlayer)
                G2G.HorizontalLine(layerSizer,G2frame.dataDisplay)   
        
        return layerSizer
        
    def OnRepSeq(event):
        event.Skip()
        stack = repseq.GetValue()
        nstar = stack.count('*')
        if nstar:
            try:
                newstack = ''
                Istar = 0
                for star in range(nstar):
                    Istar = stack.index('*',Istar+1)
                    iB = stack[:Istar].rfind(' ')
                    if iB == -1:
                        mult = int(stack[:Istar])
                    else:
                        mult = int(stack[iB:Istar])
                    pattern = stack[Istar+2:stack.index(')',Istar)]+' '
                    newstack += mult*pattern
                stack = newstack
            except ValueError:
                stack += ' Error in string'
                wx.MessageBox(stack,'Error',style=wx.ICON_EXCLAMATION)
                repseq.SetValue(data['Layer Seq'])
                return
        try:
            Slist = np.array(stack.split(),dtype=int)
        except ValueError:
            stack += ' Error in string'
            repseq.SetValue(data['Layer Seq'])
            wx.MessageBox(stack,'Error',style=wx.ICON_EXCLAMATION)
            return
        if len(Slist) < 1:
            stack += ' Error in sequence - too short!'
        Stest = np.arange(1,Nlayers-1)
        if not np.all(np.array([item in Stest for item in Slist])):
            stack += ' Error: invalid layer selection'
        elif not np.all(np.ediff1d(Slist)):
            stack += ' Error: Improbable sequence or bad string'
        if 'Error' in stack:
            repseq.SetValue(data['Layer Seq'])
            wx.MessageBox(stack,'Error',style=wx.ICON_EXCLAMATION)
            return
        else:
            data['Layer Seq'] = stack
            repseq.SetValue(stack)
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2plt.PlotPatterns(G2frame,plotType='REFD')
    
    Substances = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Substances'))['Substances']
    ProfDict,Profile,Name = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[:3]
    if 'ifDQ' not in ProfDict:
        ProfDict['ifDQ'] = np.any(Profile[5])
        data['dQ type'] = 'None'
    Limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits'))
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()   # is this a ScrolledWindow? If so, bad!
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.REFDModelMenu)
    if not G2frame.dataFrame.GetStatusBar():
        G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.SetLabel('Modelling')
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyModel, id=G2gd.wxID_MODELCOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnModelPlot, id=G2gd.wxID_MODELPLOT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitModel, id=G2gd.wxID_MODELFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitModelAll, id=G2gd.wxID_MODELFITALL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_MODELUNDO)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Reflectometry fitting for: '+Name),0,WACV)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Controls:'),0,WACV)
    mainSizer.Add(ControlSizer())
    G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)   
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Global parameters:'),0,WACV)
    mainSizer.Add(OverallSizer()) 
    G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)
    Nlayers = len(data['Layers'])
    if Nlayers > 2:
        if 'Layer Seq' not in data:
            data['Layer Seq'] = ' '.join([str(i+1) for i in range(Nlayers-2)])
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Layer sequence: '),0,WACV)
        repseq = wx.TextCtrl(G2frame.dataDisplay,value = data['Layer Seq'],style=wx.TE_PROCESS_ENTER,size=(500,25))
        repseq.Bind(wx.EVT_TEXT_ENTER,OnRepSeq)        
        repseq.Bind(wx.EVT_KILL_FOCUS,OnRepSeq)
        lineSizer.Add(repseq,0,WACV)
        mainSizer.Add(lineSizer)
        Str = ' Use sequence nos. from:'
        for ilay,layer in enumerate(data['Layers'][1:-1]):
            Str += ' %d: %s'%(ilay+1,layer['Name'])
        mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=Str),0,WACV)
        mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' NB: Repeat sequence by e.g. 6*(1 2) '),0,WACV)
    G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)    
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Layers: scatt. densities are 10%scm%s = 10%s%s%s'%(Pwr10,Pwrm2,Pwrm6,Angstr,Pwrm2)),0,WACV)
    mainSizer.Add(LayerSizer())
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    Size[1] = min(700,Size[1])
    G2frame.dataFrame.setSizePosLeft(Size)
#    x,xr,y = G2pwd.makeSLDprofile(data,Substances)
#    ModelPlot(data,x,xr,y)
    
    
################################################################################
#####  PDF controls
################################################################################           
def computePDF(G2frame,data):
    '''Calls :func:`GSASIIpwd.CalcPDF` to compute the PDF and put into the data tree array.
    Called from OnComputePDF and OnComputeAllPDF and OnComputeAllPDF in
    GSASIIimgGUI.py
    '''

    xydata = {}
    for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
        name = data[key]['Name']
        if name:
            xydata[key] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name))
    powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['Sample']['Name'])
    limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Limits'))[1]
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Instrument Parameters'))[0]
    auxPlot = G2pwd.CalcPDF(data,inst,limits,xydata)
    data['I(Q)'] = xydata['IofQ']
    data['S(Q)'] = xydata['SofQ']
    data['F(Q)'] = xydata['FofQ']
    data['G(R)'] = xydata['GofR']
    return auxPlot

def OptimizePDF(G2frame,data,showFit=True,maxCycles=5):
    '''Optimize the PDF to minimize the difference between G(r) and the expected value for
    low r (-4 pi r #density). 
    '''
    import scipy.optimize as opt
    Min,Init,Done = SetupPDFEval(G2frame,data)
    xstart = Init()
    bakMul = data['Sample Bkg.']['Mult']
    if showFit:
        rms = Min(xstart)
        print('  Optimizing corrections to improve G(r) at low r')
        if data['Sample Bkg.'].get('Refine',False):
#            data['Flat Bkg'] = 0.
            print('  start: Ruland={:.3f}, Sample Bkg mult={:.3f} (RMS:{:.4f})'.format(
                data['Ruland'],data['Sample Bkg.']['Mult'],rms))
        else:
            print('  start: Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f} (RMS:{:.4f})'.format(
                data['Flat Bkg'],data['BackRatio'],data['Ruland'],rms))
    if data['Sample Bkg.'].get('Refine',False):
        res = opt.minimize(Min,xstart,bounds=([0.01,1],[1.2*bakMul,0.8*bakMul]),
                    method='L-BFGS-B',options={'maxiter':maxCycles},tol=0.001)
    else:
        res = opt.minimize(Min,xstart,bounds=([0,None],[0,1],[0.01,1]),
                    method='L-BFGS-B',options={'maxiter':maxCycles},tol=0.001)
    Done(res['x'])
    if showFit:
        if res['success']:
            msg = 'Converged'
        else:
            msg = 'Not Converged'
        if data['Sample Bkg.'].get('Refine',False):
            print('  end:   Ruland={:.3f}, Sample Bkg mult={:.3f} (RMS:{:.4f}) *** {} ***\n'.format(
                data['Ruland'],data['Sample Bkg.']['Mult'],res['fun'],msg))
        else:
            print('  end:   Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f}) *** {} ***\n'.format(
                data['Flat Bkg'],data['BackRatio'],data['Ruland'],res['fun'],msg))
    return res['success']

def SetupPDFEval(G2frame,data):
    '''Create functions needed to optimize the PDF at low r
    '''
    Data = copy.deepcopy(data)
    xydata = {}
    for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
        name = Data[key]['Name']
        if name:
            xydata[key] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name))
    powName = Data['Sample']['Name']
    powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,powName)
    limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Limits'))[1]
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Instrument Parameters'))[0]
    numbDen = G2pwd.GetNumDensity(Data['ElList'],Data['Form Vol'])
    BkgMax = 1.
    def EvalLowPDF(arg):
        '''Objective routine -- evaluates the RMS deviations in G(r)
        from -4(pi)*#density*r for for r<Rmin
        arguments are ['Flat Bkg','BackRatio','Ruland'] scaled so that
        the min & max values are between 0 and 1. 
        '''
        if Data['Sample Bkg.'].get('Refine',False):
            R,S = arg
            Data['Sample Bkg.']['Mult'] = S
        else:
            F,B,R = arg
            Data['Flat Bkg'] = F*BkgMax
            Data['BackRatio'] = B
        Data['Ruland'] = R/10.
        G2pwd.CalcPDF(Data,inst,limits,xydata)
        # test low r computation
        g = xydata['GofR'][1][1]
        r = xydata['GofR'][1][0]
        g0 = g[r < Data['Rmin']] + 4*np.pi*r[r < Data['Rmin']]*numbDen
        M = sum(g0**2)/len(g0)
        return M
    def GetCurrentVals():
        '''Get the current ['Flat Bkg','BackRatio','Ruland'] with scaling
        '''
        if data['Sample Bkg.'].get('Refine',False):
                return [max(10*data['Ruland'],.05),data['Sample']['Mult']]
        try:
            F = data['Flat Bkg']/BkgMax
        except:
            F = 0
        return [F,data['BackRatio'],max(10*data['Ruland'],.05)]
    def SetFinalVals(arg):
        '''Set the 'Flat Bkg', 'BackRatio' & 'Ruland' values from the
        scaled, refined values and plot corrected region of G(r) 
        '''
        if data['Sample Bkg.'].get('Refine',False):
            R,S = arg
            data['Sample Bkg.']['Mult'] = S
        else:
            F,B,R = arg
            data['Flat Bkg'] = F*BkgMax
            data['BackRatio'] = B
        data['Ruland'] = R/10.
        G2pwd.CalcPDF(data,inst,limits,xydata)
    EvalLowPDF(GetCurrentVals())
    BkgMax = max(xydata['IofQ'][1][1])/50.
    return EvalLowPDF,GetCurrentVals,SetFinalVals
    
def UpdatePDFGrid(G2frame,data):
    '''respond to selection of PWDR PDF data tree item.
    '''    
    
    def PDFFileSizer():
        
        def FillFileSizer(fileSizer,key):
            #fileSizer is a FlexGridSizer(3,4)
            
            def OnSelectFile(event):
                Obj = event.GetEventObject()
                fileKey,itemKey,fmt = itemDict[Obj.GetId()]
                if itemKey == 'Name':
                    value = Obj.GetValue()
                Obj.SetValue(fmt%(value))
                data[fileKey][itemKey] = value
                data[fileKey]['Mult'] = GetExposure(value)
                ResetFlatBkg()
                wx.CallLater(100,UpdatePDFGrid,G2frame,data)
                wx.CallAfter(OnComputePDF,None)
                
            def OnMoveMult(event):
                data[key]['Mult'] += multSpin.GetValue()*0.01
                mult.SetValue(data[key]['Mult'])
                multSpin.SetValue(0)
#                ResetFlatBkg()
                wx.CallLater(100,UpdatePDFGrid,G2frame,data)
                wx.CallAfter(OnComputePDF,None)
                            
            def OnMult(invalid,value,tc):
                if invalid: return
#                ResetFlatBkg()
                wx.CallLater(100,UpdatePDFGrid,G2frame,data)
                wx.CallAfter(OnComputePDF,None)
                
            def OnRefMult(event):
                item['Refine'] = refMult.GetValue()
                if item['Refine']:
                    Status.SetStatusText('Be sure Mult is close to anticipated value. '+   \
                        'Suggest setting Flat Bkg. to 0 before Optimize Mult')
            
            def GetExposure(backFile):
                dataId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'PWDR'+dataFile[4:])
                dataComments = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,dataId,'Comments'))
                if not backFile:
                    return -1.
                backId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,backFile)
                backComments = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,backId,'Comments'))
                expTime = 1.
                sumExp = 1.
                for item in dataComments:
                    if 'exposureTime' in item:
                        expTime = float(item.split('=')[1])
                    if 'summedExposures' in item:
                        sumExp = float(item.split('=')[1])
                dataExp = expTime*sumExp
                expTime = 1.
                sumExp = 1.
                for item in backComments:
                    if 'exposureTime' in item:
                        expTime = float(item.split('=')[1])
                    if 'summedExposures' in item:
                        sumExp = float(item.split('=')[1])
                backExp = expTime*sumExp
                return -dataExp/backExp
            
            item = data[key]
            fileList = [''] + GetFileList(G2frame,'PWDR')
            fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' '+key+' file:'),0,WACV)
            fileName = wx.ComboBox(G2frame.dataDisplay,value=item['Name'],choices=fileList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            itemDict[fileName.GetId()] = [key,'Name','%s']
            fileName.Bind(wx.EVT_COMBOBOX,OnSelectFile)        
            fileSizer.Add(fileName,0,)
            fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Multiplier:'),0,WACV)
            mulBox = wx.BoxSizer(wx.HORIZONTAL)
            mult = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,item,'Mult',nDig=(10,3),
                typeHint=float,OnLeave=OnMult)
            mulBox.Add(mult,0,)
            multSpin = wx.SpinButton(G2frame.dataDisplay,style=wx.SP_VERTICAL,size=wx.Size(20,25))
            multSpin.SetRange(-1,1)
            multSpin.SetValue(0)
            multSpin.Bind(wx.EVT_SPIN, OnMoveMult)
            mulBox.Add(multSpin,0,WACV)
            fileSizer.Add(mulBox,0,WACV)
            if 'Refine' in item and item['Name']:
                refMult = wx.CheckBox(parent=G2frame.dataDisplay,label='Refine?')
                refMult.SetValue(item['Refine'])
                refMult.Bind(wx.EVT_CHECKBOX, OnRefMult)
                fileSizer.Add(refMult,0,WACV)
            else:
                fileSizer.Add((5,5),0)
                            
        def ResetFlatBkg():
            Smin = np.min(G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'PWDR'+dataFile[4:]))[1][1])
            Bmin = 0; Cmin = 0.; Cmul = 0.; CBmin = 0.
            if data['Sample Bkg.']['Name']:
                Bmin = np.min(G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['Sample Bkg.']['Name']))[1][1])
                Smin += Bmin*data['Sample Bkg.']['Mult']
            if data['Container']['Name']:
                Cmin = np.min(G2frame.PatternTree.GetItemPyData(
                G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['Container']['Name']))[1][1])
                Cmul = data['Container']['Mult']
                if data['Container Bkg.']['Name']:
                    CBmin = np.min(G2frame.PatternTree.GetItemPyData(
                        G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['Container Bkg.']['Name']))[1][1])
                    Cmin += CBmin*data['Container Bkg.']['Mult']
                Smin += Cmul*Cmin
            data['Flat Bkg'] = max(0,Smin)
                            
        PDFfileSizer = wx.BoxSizer(wx.VERTICAL)
        PDFfileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' PDF data files: '),0,WACV)
        PDFfileSizer.Add((5,5),0)    
        if 'C' in inst['Type'][0]:
            str = ' Sample file: PWDR%s   Wavelength, A: %.5f  Energy, keV: %.3f  Polariz.: %.2f '%(dataFile[4:],wave,keV,polariz)
            PDFfileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=str),0,WACV)
        PDFfileSizer.Add((5,5),0)
        fileSizer = wx.FlexGridSizer(0,5,5,1)
        select = ['Sample Bkg.','Container']
        if data['Container']['Name']:
            select.append('Container Bkg.')
        for key in select:
            FillFileSizer(fileSizer,key)
        PDFfileSizer.Add(fileSizer,0)
        return PDFfileSizer
        
    def SampleSizer():
    
        def FillElemSizer(elemSizer,ElData):
            
            def AfterChange(invalid,value,tc):
                if invalid: return
                data['Form Vol'] = max(10.0,SumElementVolumes())
                wx.CallAfter(UpdatePDFGrid,G2frame,data)
                wx.CallAfter(OnComputePDF,None)
                    
            elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
                label=' Element: '+'%2s'%(ElData['Symbol'])+' * '),0,WACV)
            num = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,ElData,'FormulaNo',nDig=(10,3),min=0.0,
                typeHint=float,OnLeave=AfterChange)
            elemSizer.Add(num,0,WACV)
            elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
                label="f': %.3f"%(ElData['fp'])+' f": %.3f'%(ElData['fpp'])+' mu: %.2f barns'%(ElData['mu']) ),
                0,WACV)
            
        def AfterChange(invalid,value,tc):
            if invalid: return
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            wx.CallAfter(OnComputePDF,None)
        
        def OnGeometry(event):
            data['Geometry'] = geometry.GetValue()
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            #UpdatePDFGrid(G2frame,data)
            wx.CallAfter(OnComputePDF,None)
        
        sampleSizer = wx.BoxSizer(wx.VERTICAL)
        if not ElList:
            sampleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample information: fill in this 1st'),0,WACV)
        else:
            sampleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample information: '),0,WACV)
        sampleSizer.Add((5,5),0)    
        Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
        Trans = G2pwd.Transmission(data['Geometry'],Abs*data['Pack'],data['Diam'])
        elemSizer = wx.FlexGridSizer(0,3,5,1)
        for El in ElList:
            FillElemSizer(elemSizer,ElList[El])
        sampleSizer.Add(elemSizer,0)
        sampleSizer.Add((5,5),0)    
        midSizer = wx.BoxSizer(wx.HORIZONTAL)
        midSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Formula volume: '),0,WACV)
        formVol = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Form Vol',nDig=(10,3),min=10.0,
            typeHint=float,OnLeave=AfterChange)
        midSizer.Add(formVol,0)
        midSizer.Add(wx.StaticText(G2frame.dataDisplay,
            label=' Theoretical absorption: %.4f cm-1 Sample absorption: %.4f cm-1'%(Abs,Abs*data['Pack'])),
            0,WACV)
        sampleSizer.Add(midSizer,0)
        sampleSizer.Add((5,5),0)
        geoBox = wx.BoxSizer(wx.HORIZONTAL)
        geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample geometry: '),0,WACV)
        choice = ['Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate']
        geometry = wx.ComboBox(G2frame.dataDisplay,value=data['Geometry'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        geometry.Bind(wx.EVT_COMBOBOX, OnGeometry)
        geoBox.Add(geometry,0)
        geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample diameter/thickness, mm: '),0,WACV)
        diam = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Diam',nDig=(10,3),min=0.01,
            typeHint=float,OnLeave=AfterChange)
        geoBox.Add(diam,0)
        sampleSizer.Add(geoBox,0)
        sampleSizer.Add((5,5),0)    
        geoBox = wx.BoxSizer(wx.HORIZONTAL)
        geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Packing: '),0,WACV)
        pack = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Pack',nDig=(10,2),min=0.01,
            typeHint=float,OnLeave=AfterChange)
        geoBox.Add(pack,0)
        geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample transmission: %.3f %%'%(Trans)),0,WACV)    
        sampleSizer.Add(geoBox,0)
        return sampleSizer
        
    def SFGctrlSizer():
        
        def OnOptimizePDF(event):
            '''Optimize Flat Bkg, BackRatio & Ruland corrections to remove spurious
            "intensity" from portion of G(r) with r<Rmin.
            Invoked by Optimize PDF button and from menu command.
            '''
            if not data['ElList']:
                G2frame.ErrorDialog('PDF error','Chemical formula not defined')
                return
            Status.SetStatusText('')
            wx.BeginBusyCursor()
            try:
                OptimizePDF(G2frame,data)
            finally:
                wx.EndBusyCursor()
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            OnComputePDF(event)        
                        
        def AfterChangeNoRefresh(invalid,value,tc):
            if invalid: return
            wx.CallAfter(OnComputePDF,None)
        
        def OnDetType(event):
            data['DetType'] = detType.GetValue()
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            wx.CallAfter(OnComputePDF,None)
        
        def OnFlatSpin(event):
            data['Flat Bkg'] += flatSpin.GetValue()*0.01*data['IofQmin']
            flatBkg.SetValue(data['Flat Bkg'])
            flatSpin.SetValue(0)        
            wx.CallAfter(OnComputePDF,None)
                
        def OnBackSlider(event):
            value = int(backSldr.GetValue())/100.
            data['BackRatio'] = value
            backVal.SetValue(data['BackRatio'])
            wx.CallAfter(OnComputePDF,None)
        
        def OnRulSlider(event):
            value = int(rulandSldr.GetValue())/1000.
            data['Ruland'] = max(0.001,value)
            rulandWdt.SetValue(data['Ruland'])
            wx.CallAfter(OnComputePDF,None)
        
        def NewQmax(invalid,value,tc):
            if invalid: return
            data['QScaleLim'][0] = 0.9*value
            SQmin.SetValue(data['QScaleLim'][0])
            wx.CallAfter(OnComputePDF,None)
        
        def OnResetQ(event):
            data['QScaleLim'][1] = qLimits[1]
            SQmax.SetValue(data['QScaleLim'][1])
            data['QScaleLim'][0] = 0.9*qLimits[1]
            SQmin.SetValue(data['QScaleLim'][0])
            wx.CallAfter(OnComputePDF,None)
            
        def OnLorch(event):
            data['Lorch'] = lorch.GetValue()
            wx.CallAfter(OnComputePDF,None)
                            
        def OnNoRing(event):
            data['noRing'] = not data['noRing']
            wx.CallAfter(OnComputePDF,None)

        sfgSizer = wx.BoxSizer(wx.VERTICAL)         
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' S(Q)->F(Q)->G(r) controls: '),0,WACV)
        sqBox.Add((1,1),1,wx.EXPAND,1)
        optB = wx.Button(G2frame.dataDisplay,label='Optimize PDF',style=wx.BU_EXACTFIT)
        optB.Bind(wx.EVT_BUTTON, OnOptimizePDF)
        sqBox.Add(optB,0,WACV|wx.ALIGN_RIGHT)
        sfgSizer.Add(sqBox,0,wx.EXPAND)
        
        sfgSizer.Add((5,5),0)
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Detector type: '),0,WACV)
        choice = ['Area detector','Point detector']
        detType = wx.ComboBox(G2frame.dataDisplay,value=data['DetType'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        detType.Bind(wx.EVT_COMBOBOX, OnDetType)
        sqBox.Add(detType,0)
        if data['DetType'] == 'Area detector':
            sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' IP transmission coeff.: '),0,WACV)
            obliqCoeff = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'ObliqCoeff',nDig=(10,3),min=0.0,max=1.0,
                typeHint=float,OnLeave=AfterChangeNoRefresh)
            sqBox.Add(obliqCoeff,0)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Flat Bkg.: '),0,WACV)
        flatBkg = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Flat Bkg',nDig=(10,0),min=0,
                typeHint=float,OnLeave=AfterChangeNoRefresh)
        sqBox.Add(flatBkg,0)
        flatSpin = wx.SpinButton(G2frame.dataDisplay,style=wx.SP_VERTICAL,size=wx.Size(20,25))
        flatSpin.SetRange(-1,1)
        flatSpin.SetValue(0)
        flatSpin.Bind(wx.EVT_SPIN, OnFlatSpin)
        sqBox.Add(flatSpin,0,WACV)
        sqBox.Add((1,1),1,wx.EXPAND,1)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label='Rmin: '),0,WACV|wx.ALIGN_RIGHT)
        rmin = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Rmin',nDig=(5,1),
                typeHint=float,size=wx.Size(50,20))
        sqBox.Add(rmin,0,WACV|wx.ALIGN_RIGHT)
        sfgSizer.Add(sqBox,0,wx.EXPAND)
            
        bkBox = wx.BoxSizer(wx.HORIZONTAL)
        bkBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Background ratio: '),0,WACV)    
        backSldr = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(100*data['BackRatio']))
        bkBox.Add(backSldr,1,wx.EXPAND)
        backSldr.Bind(wx.EVT_SLIDER, OnBackSlider)
        backVal = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'BackRatio',nDig=(10,3),min=0.0,max=1.0,
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        bkBox.Add(backVal,0,WACV)    
        sfgSizer.Add(bkBox,0,wx.ALIGN_LEFT|wx.EXPAND)
    
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Ruland width: '),0,WACV)    
        rulandSldr = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
            value=int(1000*data['Ruland']))
        sqBox.Add(rulandSldr,1,wx.EXPAND)
        rulandSldr.Bind(wx.EVT_SLIDER, OnRulSlider)
        rulandWdt = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Ruland',nDig=(10,3),min=0.001,max=1.0,
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        sqBox.Add(rulandWdt,0,WACV)    
        sfgSizer.Add(sqBox,0,wx.ALIGN_LEFT|wx.EXPAND)
        
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Scaling Q-range: '),0,WACV)
        SQmin = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['QScaleLim'],0,nDig=(10,3),
            min=qLimits[0],max=.95*data['QScaleLim'][1],
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        sqBox.Add(SQmin,0,WACV)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' to Qmax '),0,WACV)
        SQmax = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data['QScaleLim'],1,nDig=(10,3),
            min=qLimits[0],max=qLimits[1],typeHint=float,OnLeave=NewQmax)
        sqBox.Add(SQmax,0,WACV)
        resetQ = wx.Button(G2frame.dataDisplay,label='Reset?',style=wx.BU_EXACTFIT)
        sqBox.Add(resetQ,0,WACV)
        resetQ.Bind(wx.EVT_BUTTON, OnResetQ)
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Rmax: '),0,WACV)
        rmax = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'Rmax',nDig=(10,1),min=10.,max=200.,
            typeHint=float,OnLeave=AfterChangeNoRefresh,size=wx.Size(50,20))
        sqBox.Add(rmax,0,WACV)
        lorch = wx.CheckBox(parent=G2frame.dataDisplay,label='Lorch damping?')
        lorch.SetValue(data['Lorch'])
        lorch.Bind(wx.EVT_CHECKBOX, OnLorch)
        sqBox.Add(lorch,0,WACV)
        noRing = wx.CheckBox(parent=G2frame.dataDisplay,label='Suppress G(0) ringing?')
        noRing.SetValue(data['noRing'])
        noRing.Bind(wx.EVT_CHECKBOX, OnNoRing)
        sqBox.Add(noRing,0,WACV)
        sfgSizer.Add(sqBox,0)
        return sfgSizer
        
    def DiffSizer():
        
        def OnSelectGR(event):
            newName = grName.GetValue()
            if newName:
                data['delt-G(R)'] = copy.deepcopy(data['G(R)'])
                id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,newName)
                pId = G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls')
                subData = G2frame.PatternTree.GetItemPyData(pId)['G(R)']
                if subData[1][0][-1] != data['G(R)'][1][0][-1]:
                    G2frame.ErrorDialog('delt-G(R) Error',' G(R) for '+newName+' not same R range')
                    grName.SetValue(data['diffGRname'])
                    return
                data['diffGRname'] = newName
                data['delt-G(R)'][1] = np.array([subData[1][0],data['G(R)'][1][1]-subData[1][1]])
                data['delt-G(R)'][2] += ('-\n'+subData[2])
                G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='delt-G(R)')
                wx.CallAfter(UpdatePDFGrid,G2frame,data)
                
        def OnMult(invalid,value,tc):
            if invalid: return
            id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,data['diffGRname'])
            pId = G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls')
            subData = G2frame.PatternTree.GetItemPyData(pId)['G(R)']
            data['delt-G(R)'][1] = np.array([subData[1][0],data['G(R)'][1][1]-data['diffMult']*subData[1][1]])
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='delt-G(R)')
        
        diffSizer = wx.BoxSizer(wx.HORIZONTAL)
        fileList = [''] + GetFileList(G2frame,'PDF')
        diffSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Subtract G(R) for: '),0,WACV)
        grName = wx.ComboBox(G2frame.dataDisplay,value=data['diffGRname'],choices=fileList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        grName.Bind(wx.EVT_COMBOBOX,OnSelectGR)        
        diffSizer.Add(grName,0,WACV)
        if data['diffGRname']:
            diffSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Mult: '),0,WACV)
            mult = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,data,'diffMult',nDig=(10,3),
                    typeHint=float,OnLeave=OnMult)
            diffSizer.Add(mult,0,WACV)
            OnMult(False,None,None)
        return diffSizer
            
    def SumElementVolumes():
        sumVol = 0.
        ElList = data['ElList']
        for El in ElList:
            Avol = (4.*math.pi/3.)*ElList[El]['Drad']**3
            sumVol += Avol*ElList[El]['FormulaNo']
        return sumVol
        wx.CallAfter(OnComputePDF,None)
        
    def OnCopyPDFControls(event):
        import copy
        TextList = GetFileList(G2frame,'PDF')
        Source = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "PDF" pattern')
            return
        od = {'label_1':'Only copy flags','value_1':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy PDF controls','Copy controls from '+Source+' to:',TextList,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PDFlist = [TextList[i] for i in dlg.GetSelections()]
                for item in PDFlist:
                    id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    olddata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'))
                    if od['value_1']:
                        olddata['Sample Bkg.']['Refine'] = data['Sample Bkg.']['Refine']    #only one flag
                    else:
                        sample = olddata['Sample']
                        olddata.update(copy.deepcopy(data))
                        olddata['Sample'] = sample
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'),olddata)
                Status.SetStatusText('PDF controls copied')
        finally:
            dlg.Destroy()
                
    def OnSavePDFControls(event):
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II PDF controls file', pth, '', 
            'PDF controls files (*.pdfprm)|*.pdfprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .pdfprm
                filename = os.path.splitext(filename)[0]+'.pdfprm'
                File = open(filename,'w')
                File.write("#GSAS-II PDF controls file; do not add/delete items!\n")
                for item in data:
                    if item[:] not in ['I(Q)','S(Q)','F(Q)','G(R)']:
                        File.write(item+':'+unicode(data[item])+'\n')
                File.close()
                print 'PDF controls saved to: '+filename
        finally:
            dlg.Destroy()
                
    def OnLoadPDFControls(event):
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II PDF controls file', pth, '', 
            'PDF controls files (*.pdfprm)|*.pdfprm',wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                newdata = {}
                S = File.readline()
                while S:
                    if '#' in S:
                        S = File.readline()
                        continue
                    key,val = S.split(':',1)
                    try:
                        newdata[key] = eval(val)
                    #except SyntaxError:
                    except:
                        newdata[key] = val
                    S = File.readline()
                File.close()
                data.update(newdata)
        finally:
            dlg.Destroy()
        OnComputePDF(event)                
        wx.CallAfter(UpdatePDFGrid,G2frame,data)
        
    def OnAddElement(event):
        ElList = data['ElList']
        choice = ElList.keys()
        PE = G2elemGUI.PickElements(G2frame,choice)
        if PE.ShowModal() == wx.ID_OK:
            for El in PE.Elem:
                if El not in ElList:
                    try:
                        data['ElList'][El] = G2elem.GetElInfo(El,inst)
                        data['ElList'][El]['FormulaNo'] = 1.0
                    except IndexError: # happens with element Q
                        pass
            data['Form Vol'] = max(10.0,SumElementVolumes())
        PE.Destroy()
        wx.CallAfter(UpdatePDFGrid,G2frame,data)
        
    def OnDeleteElement(event):
        ElList = data['ElList']
        choice = ElList.keys()
        dlg = G2elemGUI.DeleteElement(G2frame,choice=choice)
        if dlg.ShowModal() == wx.ID_OK:
            del ElList[dlg.GetDeleteElement()]
        dlg.Destroy()
        wx.CallAfter(UpdatePDFGrid,G2frame,data)

    def OnComputePDF(event):
        '''Compute and plot PDF, in response to a menu command or a change to a
        computation parameter.
        '''
        if not data['ElList']:
            G2frame.ErrorDialog('PDF error','Chemical formula not defined')
            return
        auxPlot = computePDF(G2frame,data)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()
            Status.SetStatusText('PDF computed')
        for plot in auxPlot:
            XY = np.array(plot[:2])
            G2plt.PlotXY(G2frame,[XY,],Title=plot[2])
        if event is not None:
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='I(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='S(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='F(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='G(R)')
        else:
            G2plt.PlotISFG(G2frame,data,newPlot=False)
        
    def OnComputeAllPDF(event):
        print('Calculating PDFs...')
        choices = []
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                Name = G2frame.PatternTree.GetItemText(id)
                if Name.startswith('PDF '):
                    Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls'))
                    if not Data['ElList']:
                        print('  No chemical formula for {}'.format(Name))
                    else:
                        choices.append(Name)
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        if not choices:
            print('  No PDFs to compute\n')
            return
        od = {'label_1':'Optimize PDFs','value_1':True}
        dlg = G2G.G2MultiChoiceDialog(G2frame.dataFrame, 'Select PDFs to compute','Select PDFs',
            choices,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                results = dlg.GetSelections()
            else:
                return
        finally:
            dlg.Destroy()
        if not results:
            print('  No PDFs to compute\n')
            return
        Names = [choices[i] for i in results]
        pgbar = wx.ProgressDialog('Compute PDF','PDFs done: 0',len(Names)+1, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        notConverged = 0
        id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        N = 0
        try:
            while id:
                Name = G2frame.PatternTree.GetItemText(id)
                if Name in Names:
                    N += 1
                    msg = 'PDFs done: {} of {}'.format(N-1,len(Names))
                    if not pgbar.Update(N,msg)[0]:
                        pgbar.Destroy()
                        break
                    pId = G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls')
                    Data = G2frame.PatternTree.GetItemPyData(pId)
                    print('  Computing {}'.format(Name))
                    computePDF(G2frame,Data)
                    if od['value_1']:
                        notConverged += not OptimizePDF(G2frame,Data,maxCycles=10)
                    computePDF(G2frame,Data)
                    G2frame.PatternTree.SetItemPyData(pId,Data)
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        finally:
            pgbar.Destroy()
        if od['value_1']:
            msg = '{}/{} PDFs computed; {} unconverged'.format(N,len(Names),notConverged)
        else:
            msg = '{}/{} PDFs computed'.format(N,len(Names))
        G2frame.dataFrame.GetStatusBar().SetStatusText(msg)
        print(msg)
        # what item is being plotted? -- might be better to select from tree
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='I(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='S(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='F(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='G(R)')

    # Routine UpdatePDFGrid starts here
    global inst
    tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
    tof2q = lambda t,C:2.0*math.pi*C/t
    dataFile = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    powName = 'PWDR'+dataFile[4:]
    powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, powName)
    if powId: # skip if no matching PWDR entry
        fullLimits,limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Limits'))[:2]
        inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Instrument Parameters'))[0]
        if 'C' in inst['Type'][0]:
            wave = G2mth.getWave(inst)
            keV = 12.397639/wave
            qLimits = [tth2q(fullLimits[0],wave),tth2q(fullLimits[1],wave)]
            polariz = inst['Polariz.'][1]
        else:   #'T'of
            qLimits = [tof2q(fullLimits[1],inst['difC'][1]),tof2q(fullLimits[0],inst['difC'][1])]
            polariz = 1.0
        data['QScaleLim'][1] = min(qLimits[1],data['QScaleLim'][1])
        if data['QScaleLim'][0]:
            data['QScaleLim'][0] = max(qLimits[0],data['QScaleLim'][0])
        else:                                #initial setting at 90% of max Q
            data['QScaleLim'][0] = 0.90*data['QScaleLim'][1]
        itemDict = {}
        #patch
        if 'BackRatio' not in data:
            data['BackRatio'] = 0.
        if 'noRing' not in data:
            data['noRing'] = False
        if 'Rmax' not in data:
            data['Rmax'] = 100.
        if 'Flat Bkg' not in data:
            data['Flat Bkg'] = 0.
        if 'IofQmin' not in data:
            data['IofQmin'] = 1.0        
        if 'Rmin' not in data:
            data['Rmin'] = 1.5
        if data['DetType'] == 'Image plate':
            data['DetType'] = 'Area detector'
        if 'Refine' not in data['Sample Bkg.']:
            data['Sample Bkg.']['Refine'] = False
        if 'diffGRname' not in data:
            data['diffGRname'] = ''
        if 'diffMult' not in data:
            data['diffMult'] = 1.0
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PDFMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    if powId:
        G2frame.dataFrame.PDFMenu.EnableTop(0,enable=True)
    else:
        G2frame.dataFrame.PDFMenu.EnableTop(0,enable=False)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyPDFControls, id=G2gd.wxID_PDFCOPYCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSavePDFControls, id=G2gd.wxID_PDFSAVECONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadPDFControls, id=G2gd.wxID_PDFLOADCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddElement, id=G2gd.wxID_PDFADDELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2gd.wxID_PDFDELELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputePDF, id=G2gd.wxID_PDFCOMPUTE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputeAllPDF, id=G2gd.wxID_PDFCOMPUTEALL)

    mainSizer = wx.BoxSizer(wx.VERTICAL)
    if powId:
        ElList = data['ElList']
        mainSizer.Add(PDFFileSizer(),0,WACV)
        G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)
        mainSizer.Add(SampleSizer(),0,WACV)
        G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)
        mainSizer.Add(SFGctrlSizer(),0,WACV)
        G2G.HorizontalLine(mainSizer,G2frame.dataDisplay)
    mainSizer.Add(DiffSizer(),0,WACV)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(Size)
###############################################################################################################
#UpdatePDFPeaks: peaks in G(r)
###############################################################################################################
def UpdatePDFPeaks(G2frame,peaks,data):
    
    def limitSizer():
        
        def NewLim(invalid,value,tc):
            if invalid:
                return
            G2plt.PlotISFG(G2frame,data,newPlot=False,plotType='G(R)',peaks=peaks)
                        
        limitBox = wx.BoxSizer(wx.HORIZONTAL)
        limitBox.Add(wx.StaticText(G2frame.dataDisplay,label=' PDF Limits: '),0,WACV)
        lowLim = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,peaks['Limits'],0,nDig=(10,3),
            min=0.,max=10.,typeHint=float,OnLeave=NewLim)
        limitBox.Add(lowLim,0,WACV)
        highLim = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,peaks['Limits'],1,nDig=(10,3),
            min=peaks['Limits'][0],max=10.,typeHint=float,OnLeave=NewLim)
        limitBox.Add(highLim,0,WACV)
        return limitBox
        
    def backSizer():
        
        def NewBack(invalid,value,tc):
            if invalid:
                return
            G2plt.PlotISFG(G2frame,data,newPlot=False,plotType='G(R)',peaks=peaks)
            
        def OnRefBack(event):
            peaks['Background'][2] = refbk.GetValue()
        
        backBox = wx.wx.BoxSizer(wx.HORIZONTAL)
        backBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Background slope: '),0,WACV)
        slope = G2G.ValidatedTxtCtrl(G2frame.dataDisplay,peaks['Background'][1],1,nDig=(10,3),
            min=-4.*np.pi,max=0.,typeHint=float,OnLeave=NewBack)
        backBox.Add(slope,0,WACV)
        refbk = wx.CheckBox(parent=G2frame.dataDisplay,label=' Refine?')
        refbk.SetValue(peaks['Background'][2])
        refbk.Bind(wx.EVT_CHECKBOX, OnRefBack)
        backBox.Add(refbk,0,WACV)
        return backBox
        
    def peakSizer():
        
        def PeaksRefine(event):
            c =  event.GetCol()
            if PDFPeaks.GetColLabelValue(c) == 'refine':
                choice = ['P - position','M - magnitude','S - standrd deviation']
                dlg = wx.MultiChoiceDialog(G2frame,'Select','Refinement controls',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelections()
                    parms = ''
                    for x in sel:
                        parms += choice[x][0]
                    for peak in peaks['Peaks']:
                        peak[3] = parms
                dlg.Destroy()
                wx.CallAfter(UpdatePDFPeaks,G2frame,peaks,data)
                
        def ElTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if 'Atom' in PDFPeaks.GetColLabelValue(c):
                PE = G2elemGUI.PickElement(G2frame)
                if PE.ShowModal() == wx.ID_OK:
                    el = PE.Elem.strip()
                    peaks['Peaks'][r][c] = el
                    PDFPeaks.SetCellValue(r,c,el)
                PE.Destroy()                
        
        colLabels = ['position','magnitude','sig','refine','Atom A','Atom B','Bond No.']
        Types = 3*[wg.GRID_VALUE_FLOAT+':10,3',]+[wg.GRID_VALUE_CHOICE+': ,P,M,S,PM,PS,MS,PMS',]+     \
            2*[wg.GRID_VALUE_STRING,]+[wg.GRID_VALUE_FLOAT+':10,3',]
        rowLabels = range(len(peaks['Peaks']))
        peakTable = G2G.Table(peaks['Peaks'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        PDFPeaks = G2G.GSGrid(G2frame.dataDisplay)
        PDFPeaks.SetTable(peakTable,False)
        PDFPeaks.SetMargins(0,0)
        PDFPeaks.SetRowLabelSize(40)
        PDFPeaks.AutoSizeColumns(False)
        PDFPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, PeaksRefine)
        PDFPeaks.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, ElTypeSelect)

        peakBox = wx.BoxSizer(wx.VERTICAL)
        peakBox.Add(wx.StaticText(G2frame.dataDisplay,label=' PDF Peaks:'),0,WACV)
        peakBox.Add(PDFPeaks,0,WACV)
        
        return peakBox
        
    def OnCopyPDFPeaks(event):
        import copy
        TextList = GetFileList(G2frame,'PDF')
        Source = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy PDF peaks to','There must be more than one "PDF" pattern')
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy PDF peaks','Copy peaks from '+Source+' to:',TextList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PDFlist = [TextList[i] for i in dlg.GetSelections()]
                for item in PDFlist:
                    id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,item)
                    olddata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Peaks'))
                    olddata.update(copy.deepcopy(peaks))
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Peaks'),olddata)
                Status.SetStatusText('PDF peaks copied')
        finally:
            dlg.Destroy()
        
    def OnFitPDFpeaks(event):
        PatternId = G2frame.PatternId
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'PDF Controls'))
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'PDF Peaks'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        newpeaks = G2pwd.PDFPeakFit(peaks,data['G(R)'])[0]
        print 'PDF peak fit finished'
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'PDF Peaks'),newpeaks)
        G2plt.PlotISFG(G2frame,data,peaks=newpeaks,newPlot=False)
        wx.CallAfter(UpdatePDFPeaks,G2frame,newpeaks,data)
        
    def OnFitAllPDFpeaks(event):
        Names = G2gd.GetPatternTreeDataNames(G2frame,['PDF ',])
        dlg = G2G.G2MultiChoiceDialog(G2frame,'PDF peak fitting','Select PDFs to fit:',Names)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id =  G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequential PDF peak fit results')
                if Id:
                    SeqResult = G2frame.PatternTree.GetItemPyData(Id)
                else:
                    SeqResult = {}
                    Id = G2frame.PatternTree.AppendItem(parent=G2frame.root,text='Sequential PDF peak fit results')
                SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
                items = dlg.GetSelections()
                G2frame.EnablePlot = False
                for item in items:
                    name = Names[item]
                    pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name)
                    data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,pId, 'PDF Controls'))
                    peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,pId,'PDF Peaks'))
                    newpeaks,vals,varyList,sigList,parmDict,Rvals = G2pwd.PDFPeakFit(peaks,data['G(R)'])
                    if vals is None:
                        print 'Nothing varied!'
                        dlg.Destroy()
                        return
                    SeqResult[name] = {'variables':vals,'varyList':varyList,'sig':sigList,'Rvals':Rvals,
                        'covMatrix':np.eye(len(varyList)),'title':name,'parmDict':parmDict}
                    G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,pId, 'PDF Peaks'),newpeaks)
                SeqResult['histNames'] = Names
                G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
        finally:
            dlg.Destroy()
        G2plt.PlotISFG(G2frame,data,peaks=newpeaks,newPlot=False)
        G2frame.PatternTree.SetItemPyData(Id,SeqResult)
        G2frame.PatternTree.SelectItem(Id)
        print 'All PDFs peak fitted - results in Sequential PDF peak fit results'
        
    def OnClearPDFpeaks(event):
        peaks['Peaks'] = []
        G2plt.PlotISFG(G2frame,data,peaks=peaks,newPlot=False)
        wx.CallAfter(UpdatePDFPeaks,G2frame,peaks,data)

    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PDFPksMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyPDFPeaks, id=G2gd.wxID_PDFCOPYPEAKS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitPDFpeaks, id=G2gd.wxID_PDFPKSFIT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitAllPDFpeaks, id=G2gd.wxID_PDFPKSFITALL)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnClearPDFpeaks, id=G2gd.wxID_CLEARPDFPEAKS)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0) 
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' PDF peak fit controls:'),0,WACV)
    mainSizer.Add((5,5),0) 
    mainSizer.Add(limitSizer(),0,WACV) 
    mainSizer.Add((5,5),0) 
    mainSizer.Add(backSizer())
    if len(peaks['Peaks']):
        mainSizer.Add((5,5),0)
        mainSizer.Add(peakSizer())
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataFrame.setSizePosLeft(Size)
    
