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
import time
import copy
import random as ran
import cPickle
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
import GSASIIElemGUI as G2elemGUI
import GSASIIElem as G2elem
import GSASIIsasd as G2sasd
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WACV = wx.ALIGN_CENTER_VERTICAL
Pwr10 = unichr(0x0b9)+unichr(0x0b0)
Pwr20 = unichr(0x0b2)+unichr(0x0b0)
Pwrm1 = unichr(0x207b)+unichr(0x0b9)
Pwrm2 = unichr(0x207b)+unichr(0x0b2)
Pwrm4 = unichr(0x207b)+unichr(0x2074)   #really -d but looks like -4 as a superscript
# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
    
def IsHistogramInAnyPhase(G2frame,histoName):
    'Needs a doc string'
    phases = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Phases')
    if phases:
        item, cookie = G2frame.PatternTree.GetFirstChild(phases)
        while item:
            data = G2frame.PatternTree.GetItemPyData(item)
            histoList = data['Histograms'].keys()
            if histoName in histoList:
                return True
            item, cookie = G2frame.PatternTree.GetNextChild(phases, cookie)
        return False
    else:
        return False

def SetDefaultSample():
    'Fills in default items for the Sample dictionary'
    return {
        'InstrName':'',
        'ranId':ran.randint(0,sys.maxint),
        'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],
        'DisplaceX':[0.0,False],'DisplaceY':[0.0,False],'Diffuse':[],
        'Temperature':300.,'Pressure':0.1,
        'FreePrm1':0.,'FreePrm2':0.,'FreePrm3':0.,
        'Gonio. radius':200.0,
        'Omega':0.0,'Chi':0.0,'Phi':0.0,
#SASD items
        'Materials':[{'Name':'vacuum','VolFrac':1.0,},{'Name':'vacuum','VolFrac':0.0,}],
        'Thick':1.0,'Contrast':[0.0,0.0],       #contrast & anomalous contrast
        'Trans':1.0,                            #measured transmission
        }
def SetupSampleLabels(histName,dataType):
    '''Setup a list of labels and number formatting for use in
    labeling sample parameters.
    :param str histName: Name of histogram, ("PWDR ...")
    :param str dataType: 
    '''
    parms = []
    parms.append(['Scale','Histogram scale factor: ',[10,4]])
    parms.append(['Gonio. radius','Goniometer radius (mm): ',[10,3]])
    if 'PWDR' in histName:
        if dataType == 'Debye-Scherrer':
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
    parms.append(['Omega','Goniometer omega:',[10,3]])
    parms.append(['Chi','Goniometer chi:',[10,3]])
    parms.append(['Phi','Goniometer phi:',[10,3]])
    parms.append(['Temperature','Sample temperature (K): ',[10,3]])
    parms.append(['Pressure','Sample pressure (MPa): ',[10,3]])
    return parms

def SetDefaultSASDModel():
    'Fills in default items for the SASD Models dictionary'    
    return {'Back':[0.0,False],'Size':{'MinDiam':50,'MaxDiam':10000,'Nbins':100,
        'logBins':True,'Method':'MaxEnt','Distribution':[],
        'Shape':['Spheroid',1.0],'MaxEnt':{'Niter':100,'Precision':0.01,'Sky':-3},
        'IPG':{'Niter':100,'Approach':0.8,'Power':-1},'Reg':{},},            
        'Particle':{'Matrix':{'Name':'vacuum','VolFrac':[0.0,False]},'Levels':[],},
        'Current':'Size dist.',
        }
        
def SetDefaultSubstances():
    'Fills in default items for the SASD Substances dictionary'
    return {'Substances':{'vacuum':{'Elements':{},'Volume':1.0,'Density':0.0,'Scatt density':0.0}}}

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

def CopySelectedHistItems(G2frame):
    '''Global copy: Copy items from current histogram to others.
    This is called from the menubar and is available only when the top histogram tree entry
    is selected.
    '''
    hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    histList = GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
        return
    choices = ['Limits','Background','Instrument Parameters','Sample Parameters']
    dlg = G2gd.G2MultiChoiceDialog(
        G2frame.dataFrame, 
        'Copy which histogram sections from\n'+str(hst[5:]),
        'Select copy sections', choices, filterBox=False)
    dlg.SetSelections(range(len(choices)))
    choiceList = []
    if dlg.ShowModal() == wx.ID_OK:
        choiceList = [choices[i] for i in dlg.GetSelections()]
    if not choiceList: return
    
    dlg = G2gd.G2MultiChoiceDialog(
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
        PickId = G2frame.PickId
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        profile = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        x0 = profile[0]
        iBeg = np.searchsorted(x0,limits[0])
        iFin = np.searchsorted(x0,limits[1])
        x = x0[iBeg:iFin]
        y0 = profile[1][iBeg:iFin]
        y1 = copy.copy(y0)
        ysig = np.std(y1)
        offset = [-1,1]
        ymask = ma.array(y0,mask=(y0<ysig))
        for off in offset:
            ymask = ma.array(ymask,mask=(ymask-np.roll(y0,off)<=0.))
        indx = ymask.nonzero()
        mags = ymask[indx]
        poss = x[indx]
        refs = zip(poss,mags)
        refs = G2mth.sortArray(refs,0,reverse=True)
        for i,ref1 in enumerate(refs):
            for ref2 in refs[i+1:]:
                if abs(ref2[0]-ref1[0]) < 0.1*G2pwd.getFWHM(ref1[0],inst):
                    del(refs[i])    
        refs = G2mth.sortArray(refs,1,reverse=True)
        for pos,mag in refs:
            data.append(G2mth.setPeakparms(inst,inst2,pos,mag))
        UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')        
    
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
        if not G2frame.GSASprojectfile:            #force a save of the gpx file so SaveState can wirte in the same directory
            G2frame.OnFileSaveas(event)
        OnPeakFit('LSQ')
        
    def OnOneCycle(event):
        OnPeakFit('LSQ',oneCycle=True)
        
    def OnClearPeaks(event):
        dlg = wx.MessageDialog(G2frame,'Delete all peaks?','Clear peak list',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = []
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
        PickId = G2frame.PickId
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        background = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        data = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        wx.BeginBusyCursor()
        dlg = wx.ProgressDialog('Residual','Peak fit Rwp = ',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            G2pwd.DoPeakFit(FitPgm,peaks,background,limits,inst,inst2,data,oneCycle,controls,dlg)
        finally:
            wx.EndBusyCursor()    
        UpdatePeakGrid(G2frame,peaks)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        print 'finished'
        return
        
    def OnResetSigGam(event):
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId
        Inst,Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to do!')
            return
        newpeaks = []
        for peak in peaks:
            newpeaks.append(G2mth.setPeakparms(Inst,Inst2,peak[0],peak[2]))
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Peak List'),newpeaks)
        UpdatePeakGrid(G2frame,newpeaks)
                
    def RefreshPeakGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        
        event.StopPropagation()
        data = G2frame.PeakTable.GetData()
        T = []
        for peak in data:T.append(peak[0])
        D = dict(zip(T,data))
        T.sort()
        X = []
        for key in T: X.append(D[key])
        data = X        
        
    def setBackgroundColors():
       for r in range(G2frame.dataDisplay.GetNumberRows()):
           for c in range(G2frame.dataDisplay.GetNumberCols()):
               if G2frame.dataDisplay.GetColLabelValue(c) in ['position','intensity','alpha','beta','sigma','gamma']:
                   if float(G2frame.dataDisplay.GetCellValue(r,c)) < 0.:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.RED)
                   else:
                       G2frame.dataDisplay.SetCellBackgroundColour(r,c,wx.WHITE)
                                                  
    def KeyEditPeakGrid(event):
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
        elif rowList:
            G2frame.dataDisplay.ClearSelection()
            if event.GetKeyCode() == wx.WXK_DELETE:
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
                data = G2frame.PeakTable.GetData()
                G2frame.PatternTree.SetItemPyData(G2frame.PickId,data[:-nDel])
                G2frame.dataDisplay.ForceRefresh()
                setBackgroundColors()
                if not len(G2frame.PatternTree.GetItemPyData(G2frame.PickId)): 
                    G2frame.dataFrame.PeakFit.Enable(False)
                        
        elif colList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.PeakTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data[row][col]=False
        elif selectList:
            G2frame.dataDisplay.ClearSelection()
            key = event.GetKeyCode()
            for row,col in selectList:
                if G2frame.PeakTable.GetTypeName(row,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        data[row][col]=True
                    elif key == 78:  #'N'
                        data[row][col]=False
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
            
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PeakMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    Status.SetStatusText('Global refine: select refine column & press Y or N')
    G2frame.Bind(wx.EVT_MENU, OnAutoSearch, id=G2gd.wxID_AUTOSEARCH)
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2gd.wxID_UNDO)
    G2frame.Bind(wx.EVT_MENU, OnLSQPeakFit, id=G2gd.wxID_LSQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnOneCycle, id=G2gd.wxID_LSQONECYCLE)
    G2frame.Bind(wx.EVT_MENU, OnClearPeaks, id=G2gd.wxID_CLEARPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnResetSigGam, id=G2gd.wxID_RESETSIGGAM)
    G2frame.dataFrame.PeakFit.Enable(False)
    if data:
        G2frame.dataFrame.PeakFit.Enable(True)
        G2frame.dataFrame.PFOneCycle.Enable(True)
    G2frame.PickTable = []
    rowLabels = []
    PatternId = G2frame.PatternId
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Instrument Parameters'))[0]
    for i in range(len(data)): rowLabels.append(str(i+1))
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
    for peak in data:
        T.append(peak[0])
    D = dict(zip(T,data))
    T.sort()
    X = []
    for key in T: X.append(D[key])
    data = X
    G2frame.PatternTree.SetItemPyData(G2frame.PickId,data)
    G2frame.PeakTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Peak List')
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.PeakTable, True)
    setBackgroundColors()                         
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshPeakGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([535,350])
        
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
    ValObj = {}
    
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
        dlg = G2gd.G2MultiChoiceDialog(
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
        dlg = G2gd.G2MultiChoiceDialog(
            G2frame.dataFrame, 
            'Copy bkg params from\n'+str(hst[5:])+' to...',
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
            wx.CallAfter(UpdateBackground,G2frame,data)
            
        def OnBakVal(event):
            Obj = event.GetEventObject()
            item = ValObj[Obj.GetId()][0]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = data[0][item]
            data[0][item] = value
            Obj.SetValue('%10.4f'%(value))
        
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
        bakSizer = wx.FlexGridSizer(1,5,5,5)
        for i,value in enumerate(data[0][3:]):
            bakVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.4f'%(value),style=wx.TE_PROCESS_ENTER)
            bakSizer.Add(bakVal,0,WACV)
            ValObj[bakVal.GetId()] = [i+3]
            bakVal.Bind(wx.EVT_TEXT_ENTER,OnBakVal)
            bakVal.Bind(wx.EVT_KILL_FOCUS,OnBakVal)
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
            wx.CallAfter(UpdateBackground,G2frame,data)
            
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
            debyeTable = G2gd.Table(data[1]['debyeTerms'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            debyeGrid = G2gd.GSGrid(parent=G2frame.dataDisplay)
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
            wx.CallAfter(UpdateBackground,G2frame,data)
            
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
        peaks = wx.ComboBox(G2frame.dataDisplay,-1,value=str(data[1]['nPeaks']),choices=[str(i) for i in range(12)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        peaks.Bind(wx.EVT_COMBOBOX,OnPeaks)
        topSizer.Add(peaks,0,WACV)
        topSizer.Add((5,0),0)
        peaksSizer.Add(topSizer)
        if data[1]['nPeaks']:
            peaksSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Peak list:'),0,WACV)       
            rowLabels = []
            for i in range(len(data[1]['peaksList'])): rowLabels.append(str(i))
            colLabels = ['pos','refine','int','refine','sig','refine','gam','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            peaksTable = G2gd.Table(data[1]['peaksList'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            peaksGrid = G2gd.GSGrid(parent=G2frame.dataDisplay)
            peaksGrid.SetTable(peaksTable, True)
            peaksGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            peaksGrid.AutoSizeColumns(False)
            peaksSizer.Add(peaksGrid)        
        return peaksSizer
                
    if G2frame.dataDisplay:
#        G2frame.dataFrame.Clear()
        G2frame.dataFrame.DestroyChildren()
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.BackMenu)
    G2frame.dataFrame.SetLabel('Background')
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnBackCopy,id=G2gd.wxID_BACKCOPY)
    G2frame.Bind(wx.EVT_MENU,OnBackFlagCopy,id=G2gd.wxID_BACKFLAGCOPY)
    BackId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Background')
    Choices = ['chebyschev','cosine','lin interpolate','inv interpolate','log interpolate']
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(BackSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(DebyeSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(PeaksSizer())
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
        
################################################################################
#####  Limits
################################################################################           
       
def UpdateLimitsGrid(G2frame, data,plottype):
    '''respond to selection of PWDR Limits data tree item.
    '''
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2frame.ifGetExclude = False
        
    def KeyEditPeakGrid(event):
        if event.GetKeyCode() == wx.WXK_DELETE:
            row = G2frame.dataDisplay.GetSelectedRows()[0]
            if row > 1: #can't delete limits!
                del(data[row])
                wx.CallAfter(UpdateLimitsGrid,G2frame,data,plottype)
                G2plt.PlotPatterns(G2frame,plotType=plottype)
                        
    def RefreshLimitsGrid(event):
        event.StopPropagation()
        data = G2frame.LimitsTable.GetData()
        old = data[0]
        new = data[1]
        new[0] = max(old[0],new[0])
        new[1] = max(new[0],min(old[1],new[1]))
        excl = []
        if len(data) > 2:
            excl = data[2:]
            for item in excl:
                item[0] = max(old[0],item[0])
                item[1] = max(item[0],min(old[1],item[1]))
        data = [old,new]+excl
        G2frame.LimitsTable.SetData(data)
        G2plt.PlotPatterns(G2frame,plotType=plottype)
        
    def OnLimitCopy(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2gd.G2MultiChoiceDialog(
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
            
    def OnAddExcl(event):
        G2frame.ifGetExclude = True
        print 'Add excluded region'
        
    G2frame.LimitsTable = []
    colLabels = ['Tmin','Tmax']
    rowLabels = ['original','changed']
    for i in range(len(data)-2):
        rowLabels.append('exclude')
    Types = 2*[wg.GRID_VALUE_FLOAT+':12,5',]
    G2frame.LimitsTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Limits')
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.LimitMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU,OnLimitCopy,id=G2gd.wxID_LIMITCOPY)
    G2frame.Bind(wx.EVT_MENU,OnAddExcl,id=G2gd.wxID_ADDEXCLREGION)    
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.LimitsTable, True)    
    G2frame.dataDisplay.SetCellStyle(0,0,VERY_LIGHT_GREY,True)
    G2frame.dataDisplay.SetCellStyle(0,1,VERY_LIGHT_GREY,True)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshLimitsGrid)                
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([230,260])
                                
    
################################################################################
#####  Instrument parameters
################################################################################           
       
def UpdateInstrumentGrid(G2frame,data):
    '''respond to selection of PWDR/SASD Instrument Parameters
    data tree item.
    '''
    def keycheck(keys):
        good = []
        for key in keys:
            if key in ['Type','U','V','W','X','Y','SH/L','I(L2)/I(L1)','alpha',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-q','Polariz.',
                'Lam','Azimuth','2-theta','difC','difA','Zero','Lam1','Lam2']:
                good.append(key)
        return good
        
    keys = keycheck(data.keys())
    if 'P' in data['Type'][0]:          #powder data
        insVal = dict(zip(keys,[data[key][1] for key in keys]))
        insDef = dict(zip(keys,[data[key][0] for key in keys]))
        insRef = dict(zip(keys,[data[key][2] for key in keys]))
        if 'NC' in data['Type'][0]:
            del(insDef['Polariz.'])
            del(insVal['Polariz.'])
            del(insRef['Polariz.'])
    elif 'S' in data['Type'][0]:                               #single crystal data
        insVal = dict(zip(keys,[data[key][1] for key in keys]))
        insDef = dict(zip(keys,[data[key][0] for key in keys]))
        insRef = {}
    elif 'L' in data['Type'][0]:                               #low angle data
        insVal = dict(zip(keys,[data[key][1] for key in keys]))
        insDef = dict(zip(keys,[data[key][0] for key in keys]))
        insRef = {}
    ValObj = {}
    RefObj = {}
    waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
        'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'MoKa':[0.70926,0.713543],
        'AgKa':[0.559363,0.563775]}
    Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[1]
        
    def inst2data(inst,ref,data):
        for item in data:
            try:
                data[item] = [data[item][0],inst[item],ref[item]]
            except KeyError:
                pass        #skip 'Polariz.' for N-data
        return data
        
    def updateData(inst,ref):
        return inst2data(inst,ref,G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[0])        
    
    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
            newpeaks = []
            for peak in peaks:
                newpeaks.append(G2mth.setPeakparms(data,Inst2,peak[0],peak[2]))
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'),newpeaks)
                    
    def OnLoad(event):
        '''Loads instrument parameters from a G2 .instprm file
        in response to the Instrument Parameters-Operations/Load Profile menu
        
        Note that similar code is found in ReadPowderInstprm (GSASII.py)
        '''
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', '.', '', 
            'instrument parameter files (*.instprm)|*.instprm',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                newItems = []
                newVals = []
                while S:
                    if S[0] == '#':
                        S = File.readline()
                        continue
                    [item,val] = S[:-1].split(':')
                    newItems.append(item)
                    try:
                        newVals.append(float(val))
                    except ValueError:
                        newVals.append(val)                        
                    S = File.readline()                
                File.close()
                Inst,Inst2 = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Instrument Parameters'))
                inst = G2IO.makeInstDict(newItems,newVals,len(newVals)*[False,])
                G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId,'Instrument Parameters'),[inst,Inst2])
                RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
                UpdateInstrumentGrid(G2frame,data)
                G2plt.PlotPeakWidths(G2frame)
        finally:
            dlg.Destroy()
        
    def OnSave(event):
        '''Respond to the Instrument Parameters Operations/Save Profile menu
        item: writes current parameters to a .instprm file
        '''
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', '.', '', 
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .instprm
                filename = os.path.splitext(filename)[0]+'.instprm'
                File = open(filename,'w')
                File.write("#GSAS-II instrument parameter file; do not add/delete or change order of items!\n")
                for item in data:
                    File.write(item+':'+str(data[item][1])+'\n')
                File.close()
        finally:
            dlg.Destroy()
                                                
    def OnReset(event):
        insVal.update(insDef)
        data = updateData(insVal,insRef)
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
        flags = dict(zip(keys,[data[key][2] for key in keys]))
        instType = data['Type'][0]
        copyList = []
        dlg = G2gd.G2MultiChoiceDialog(
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
            if len(data) == len(instData) and instType == instData['Type'][0]:   #don't mix data types or lam & lam1/lam2 parms!
                for item in instData:
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
        dlg = G2gd.G2MultiChoiceDialog(
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
            if len(data) == len(instData) and instType == instData['Type'][0]:  #don't mix data types or lam & lam1/lam2 parms!
                instData.update(data)
            else:
                print item+' not copied - instrument parameters not commensurate'
        
    def OnWaveChange(event):
        if 'Lam' in insVal:
            data['Lam1'] = [waves['CuKa'][0],waves['CuKa'][0],0]
            data['Lam2'] = [waves['CuKa'][1],waves['CuKa'][1],0]
            data['I(L2)/I(L1)'] = [0.5,0.5,0]
            del(data['Lam'])
        else:
            data['Lam'] = [data['Lam1'][0],data['Lam1'][0],0]
            del(data['Lam1'])
        wx.CallAfter(UpdateInstrumentGrid,data)
                
    def OnLamPick(event):
        lamType = lamPick.GetValue()
        insVal['Lam1'] = waves[lamType][0]
        insVal['Lam2'] = waves[lamType][1]
        data = updateData(insVal,insRef)
        UpdateInstrumentGrid(G2frame,data)
                 
    def OnRatValue(event):
        try:
            value = float(ratVal.GetValue())
            if value < 0:
                raise ValueError
        except ValueError:
            value = insVal['I(L2)/I(L1)']
        insVal['I(L2)/I(L1)'] = value
        ratVal.SetValue('%10.4f'%(value))
        data = updateData(insVal,insRef)
        
    def OnRatRef(event):
        insRef['I(L2)/I(L1)'] = ratRef.GetValue()
        data = updateData(insVal,insRef)
        
    def OnWaveValue(event):
        try:
            value = float(waveVal.GetValue())
            if value < 0:
                raise ValueError
        except ValueError:
            value = insVal['Lam']
        insVal['Lam'] = value
        waveVal.SetValue('%10.6f'%(value))
        data = updateData(insVal,insRef)
        
    def OnWaveRef(event):
        insRef['Lam'] = waveRef.GetValue()
        data = updateData(insVal,insRef)
        
    def OnItemValue(event):
        Obj = event.GetEventObject()
        item,fmt = ValObj[Obj.GetId()]
        try:
            value = float(Obj.GetValue())
        except ValueError:
            value = insVal[item]
        insVal[item] = value
        Obj.SetValue(fmt%(value))
        data = updateData(insVal,insRef)
        G2plt.PlotPeakWidths(G2frame)
        
    def OnItemRef(event):
        Obj = event.GetEventObject()
        item = RefObj[Obj.GetId()]
        insRef[item] = Obj.GetValue()
        data = updateData(insVal,insRef)
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    try:
        histoName = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[-1]
        ifHisto = IsHistogramInAnyPhase(G2frame,histoName)
    except TypeError:       #PKS data never used in a phase as data
        ifhisto = False
    G2gd.SetDataMenuBar(G2frame)
    G2frame.dataFrame.SetLabel('Instrument Parameters')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    topSizer = wx.FlexGridSizer(1,6,5,5)
    instSizer = wx.FlexGridSizer(2,6,5,5)
    topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Histogram Type: '+insVal['Type']),0,WACV)
#    topSizer.Add((5,5),0)
    if 'P' in insVal['Type']:                   #powder data
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.InstMenu)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()
        G2frame.Bind(wx.EVT_MENU,OnLoad,id=G2gd.wxID_INSTLOAD)
        G2frame.Bind(wx.EVT_MENU,OnSave,id=G2gd.wxID_INSTSAVE)
        G2frame.Bind(wx.EVT_MENU,OnReset,id=G2gd.wxID_INSTPRMRESET)
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2gd.wxID_INSTCOPY)
        G2frame.Bind(wx.EVT_MENU,OnInstFlagCopy,id=G2gd.wxID_INSTFLAGCOPY)
        G2frame.Bind(wx.EVT_MENU,OnWaveChange,id=G2gd.wxID_CHANGEWAVETYPE)        
        if 'C' in insVal['Type']:               #constant wavelength
            #patch
            if 'Azimuth' not in insVal:
                insVal['Azimuth'] = 0.0
                insDef['Azimuth'] = 0.0
                insRef['Azimuth'] = False
            #end of patch
            topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Azimuth: %7.2f'%(insVal['Azimuth'])),0,WACV)
            if 'Lam1' in insVal:
                topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Ka1/Ka2:'),
                        0,WACV)
                topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u'  %8.6f/%8.6f\xc5'%(insVal['Lam1'],insVal['Lam2'])),
                        0,WACV)
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Select: '),0,WACV)
                choice = ['TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa']
                lamPick = wx.ComboBox(G2frame.dataDisplay,value=' ',choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0)
                topSizer.Add(waveSizer,0)
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' I(L2)/I(L1): (%10.4f)'%(insDef['I(L2)/I(L1)'])),
                        0,WACV)
                ratVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.4f'%(insVal['I(L2)/I(L1)']),style=wx.TE_PROCESS_ENTER)
                ratVal.Bind(wx.EVT_TEXT_ENTER,OnRatValue)
                ratVal.Bind(wx.EVT_KILL_FOCUS,OnRatValue)
                instSizer.Add(ratVal,0)
                ratRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
                ratRef.SetValue(bool(insRef['I(L2)/I(L1)']))
                ratRef.Bind(wx.EVT_CHECKBOX, OnRatRef)
                instSizer.Add(ratRef,0,WACV)
                
            else:
                topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),
                    0,WACV)
                waveVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.6f'%(insVal['Lam']),style=wx.TE_PROCESS_ENTER)
                waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveValue)
                waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveValue)
                topSizer.Add(waveVal,0,WACV)
                if ifHisto:
                    waveRef = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
                    waveRef.SetValue(bool(insRef['Lam']))
                    waveRef.Bind(wx.EVT_CHECKBOX, OnWaveRef)
                    topSizer.Add(waveRef,0,WACV)
            for item in ['Zero','Polariz.']:
                fmt = '%10.4f'
                Fmt = ' %s: ('+fmt+')'
                if item in insDef:
                    instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,Fmt%(item,insDef[item])),
                            0,WACV)
                    itemVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                    ValObj[itemVal.GetId()] = [item,fmt]
                    itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                    itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                    instSizer.Add(itemVal,0,WACV)
                    if ifHisto:
                        itemRef = wx.CheckBox(G2frame.dataDisplay,wx.ID_ANY,label=' Refine?')
                        itemRef.SetValue(bool(insRef[item]))
                        RefObj[itemRef.GetId()] = item
                        itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                        instSizer.Add(itemRef,0,WACV)
                    else:
                        instSizer.Add((5,5),0)
                else:                           #skip Polariz. for neutrons
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)
                    instSizer.Add((5,5),0)
            for item in ['U','V','W','X','Y','SH/L']:
                fmt = '%10.3f'
                if item == 'SH/L':
                    fmt = '%10.5f'
                Fmt = ' %s: ('+fmt+')'
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,Fmt%(item,insDef[item])),
                        0,WACV)
                itemVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                ValObj[itemVal.GetId()] = [item,fmt]
                itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                instSizer.Add(itemVal,0,WACV)
                itemRef = wx.CheckBox(G2frame.dataDisplay,wx.ID_ANY,label=' Refine?')
                itemRef.SetValue(bool(insRef[item]))
                RefObj[itemRef.GetId()] = item
                itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                instSizer.Add(itemRef,0,WACV)
        else:                                   #time of flight (neutrons)
            topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Azimuth: %7.2f'%(insVal['Azimuth'])),0,WACV)
            topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' 2-theta: %7.2f'%(insVal['2-theta'])),0,WACV)
            if 'Pdabc' in Inst2:
                Items = ['sig-0','sig-1','X','Y']
                topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' difC: %8.2f'%(insVal['difC'])),0,WACV)
                topSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' alpha, beta: fixed by table'),0,WACV)
            else:
                Items = ['difC','difA','Zero','alpha','beta-0','beta-1','beta-q','sig-0','sig-1','sig-q','X','Y']
            for item in Items:
                fmt = '%10.3f'
                if 'beta' in item:
                    fmt = '%12.4g'
                Fmt = ' %s: ('+fmt+')'
                instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,Fmt%(item,insDef[item])),
                        0,WACV)
                itemVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,fmt%(insVal[item]),style=wx.TE_PROCESS_ENTER)
                ValObj[itemVal.GetId()] = [item,fmt]
                itemVal.Bind(wx.EVT_TEXT_ENTER,OnItemValue)
                itemVal.Bind(wx.EVT_KILL_FOCUS,OnItemValue)
                instSizer.Add(itemVal,0,WACV)
                if not ifHisto and item in ['difC','difA','Zero',]:
                    instSizer.Add((5,5),0)
                else:
                    itemRef = wx.CheckBox(G2frame.dataDisplay,wx.ID_ANY,label=' Refine?')
                    itemRef.SetValue(bool(insRef[item]))
                    RefObj[itemRef.GetId()] = item
                    itemRef.Bind(wx.EVT_CHECKBOX, OnItemRef)
                    instSizer.Add(itemRef,0,WACV)
        
    elif 'S' in insVal['Type']:                       #single crystal data
        if 'C' in insVal['Type']:               #constant wavelength
            instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),
                0,WACV)
            waveVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.6f'%(insVal['Lam']),style=wx.TE_PROCESS_ENTER)
            waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveValue)
            waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveValue)
            instSizer.Add(waveVal,0,WACV)
        else:                                   #time of flight (neutrons)
            pass                                #for now
    elif 'L' in insVal['Type']:
        if 'C' in insVal['Type']:        
            instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),
                0,WACV)
            waveVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,'%10.6f'%(insVal['Lam']),style=wx.TE_PROCESS_ENTER)
            waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveValue)
            waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveValue)
            instSizer.Add(waveVal,0,WACV)
            instSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,'  Azimuth: %7.2f'%(insVal['Azimuth'])),0,WACV)
        else:                                   #time of flight (neutrons)
            pass                                #for now
        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(topSizer,0)
    mainSizer.Add(instSizer,0)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
    
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
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II sample parameters file', '.', '', 
            'sample parameter files (*.samprm)|*.samprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
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
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II sample parameters file', '.', '', 
            'sample parameter files (*.samprm)|*.samprm',wx.OPEN|wx.CHANGE_DIR)
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
    
    def OnSetScale(event):
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
        UpdateSampleGrid(G2frame,data)       
        G2plt.PlotPatterns(G2frame,plotType='SASD',newPlot=False)
        
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
        dlg = G2gd.G2MultiChoiceDialog(
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
                     SetupSampleLabels(hst,data.get('Type'))
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
        dlg = G2gd.G2MultiChoiceDialog(
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
        dlg = G2gd.G2MultiChoiceDialog(
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
        dlg = G2gd.G2MultiChoiceDialog(
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

    def OnHistoType(event):
        Obj = event.GetEventObject()
        data['Type'] = Obj.GetValue()
        if data['Type'] == 'Bragg-Brentano' and 'Shift' not in data:    #set up defaults for new type(s)
            data['Shift'] = [0.0,False]
            data['Transparency'] = [0.0,False]
        wx.CallAfter(UpdateSampleGrid,G2frame,data)
        
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
        id,key = Info[Obj.GetId()]
        if key == 'Name':
            data['Materials'][id][key] = Obj.GetValue()
        elif key == 'VolFrac':
            try:
                value = min(max(0.,float(Obj.GetValue())),1.)
            except ValueError:
                value = data['Materials'][id][key]
            data['Materials'][id][key] = value
            data['Materials'][not id][key] = 1.-value
        wx.CallAfter(UpdateSampleGrid,G2frame,data)

    ######## DEBUG #######################################################
    #import GSASIIpwdGUI
    #reload(GSASIIpwdGUI)
    #reload(G2gd)
    ######################################################################
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
    if 'SASD' in histName:
        G2frame.dataFrame.SetScale.Enable(True)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
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
    if type(data['Temperature']) is int:
        data['Temperature'] = float(data['Temperature'])
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
#patch end
    
    parms = SetupSampleLabels(histName,data.get('Type'))
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
    instNameVal = wx.TextCtrl(G2frame.dataDisplay,wx.ID_ANY,data.get('InstrName',''),
                              size=(200,-1),style=wx.TE_PROCESS_ENTER)        
    nameSizer.Add(instNameVal)
    instNameVal.Bind(wx.EVT_CHAR,OnNameVal)
    mainSizer.Add(nameSizer,0,wx.EXPAND,1)
    mainSizer.Add((5,5),0)

    if 'PWDR' in histName:
        nameSizer = wx.BoxSizer(wx.HORIZONTAL)
        nameSizer.Add(wx.StaticText(G2frame.dataDisplay,wx.ID_ANY,' Diffractometer type: '),
                    0,WACV)
        choices = ['Debye-Scherrer','Bragg-Brentano',]
        histoType = wx.ComboBox(G2frame.dataDisplay,wx.ID_ANY,value=data['Type'],choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        histoType.Bind(wx.EVT_COMBOBOX, OnHistoType)
        nameSizer.Add(histoType)
        mainSizer.Add(nameSizer,0,wx.EXPAND,1)
        mainSizer.Add((5,5),0)

    parmSizer = wx.FlexGridSizer(10,2,5,0)
    for key,lbl,nDig in parms:
        if 'list' in str(type(data[key])):
            parmRef = G2gd.G2CheckBox(G2frame.dataDisplay,' '+lbl,data[key],1)
            parmSizer.Add(parmRef,0,WACV|wx.EXPAND)
            parmVal = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,data[key],0,
                nDig=nDig,typeHint=float,OnLeave=AfterChange)
        else:
            parmSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' '+lbl),
                0,WACV|wx.EXPAND)
            parmVal = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,data,key,
                typeHint=float,OnLeave=AfterChange)
        parmSizer.Add(parmVal,1,wx.EXPAND)
    Info = {}

        
    for key in ('FreePrm1','FreePrm2','FreePrm3'):
        parmVal = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,Controls,key,typeHint=str,
                                        notBlank=False)
        parmSizer.Add(parmVal,1,wx.EXPAND)
        parmVal = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,data,key,typeHint=float)
        parmSizer.Add(parmVal,1,wx.EXPAND)
    mainSizer.Add(parmSizer,1,wx.EXPAND)
    mainSizer.Add((0,5),0)    
    if 'SASD' in histName:
        rho = [0.,0.]
        anomrho = [0.,0.]
        mu = 0.
        subSizer = wx.FlexGridSizer(1,4,5,5)
        Substances = G2frame.PatternTree.GetItemPyData(
            G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Substances'))
        for id,item in enumerate(data['Materials']):
            subSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Material: '),0,WACV)
            matsel = wx.ComboBox(G2frame.dataDisplay,value=item['Name'],choices=Substances['Substances'].keys(),
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Info[matsel.GetId()] = [id,'Name']
            matsel.Bind(wx.EVT_COMBOBOX,OnMaterial)        
            subSizer.Add(matsel,0,WACV)
            subSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Volume fraction: '),0,WACV)
            volfrac = wx.TextCtrl(G2frame.dataDisplay,value=str('%.3f'%(item['VolFrac'])),style=wx.TE_PROCESS_ENTER)
            Info[volfrac.GetId()] = [id,'VolFrac']
            volfrac.Bind(wx.EVT_TEXT_ENTER,OnMaterial)
            volfrac.Bind(wx.EVT_KILL_FOCUS,OnMaterial)
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
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)
                
################################################################################
#####  Indexing Peaks
################################################################################           
       
def UpdateIndexPeaksGrid(G2frame, data):
    '''respond to selection of PWDR Index Peak List data
    tree item.
    '''
    IndexId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List')
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    wave = G2mth.getWave(Inst)
    
    def RefreshIndexPeaksGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        data = G2frame.IndexPeaksTable.GetData()
        if c == 2:
            if data[r][c]:
                data[r][c] = False
            else:
                data[r][c] = True
            G2frame.IndexPeaksTable.SetData(data)
            G2frame.PatternTree.SetItemPyData(IndexId,data)
            G2frame.dataDisplay.ForceRefresh()
            
    def OnReload(event):
        data = []
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
        for peak in peaks:
            dsp = wave/(2.0*sind((peak[0]-Inst['Zero'][1])/2.0))
            data.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])
        G2frame.PatternTree.SetItemPyData(IndexId,data)
        UpdateIndexPeaksGrid(G2frame,data)
        
    def KeyEditPickGrid(event):
        colList = G2frame.dataDisplay.GetSelectedCols()
        rowList = G2frame.dataDisplay.GetSelectedRows()
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
                if G2frame.IndexPeaksTable.GetColLabelValue(col) in ['use','refine']:
                    if key == 89: #'Y'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[row][col]=False
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    if 'PWD' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.IndPeaksMenu)
        G2frame.Bind(wx.EVT_MENU, OnReload, id=G2gd.wxID_INDXRELOAD)
    G2frame.dataFrame.IndexPeaks.Enable(False)
    G2frame.IndexPeaksTable = []
    if data:
        G2frame.dataFrame.IndexPeaks.Enable(True)
        cells = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List'))
        if cells:
            cellist = cells[2]
            dmin = cells[3]
            G2frame.HKL = []
            for i,cell in enumerate(cellist):
                if cell[-1]:
                    ibrav = cell[2]
                    A = G2lat.cell2A(cell[3:9])
                    G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A)
                    G2indx.IndexPeaks(data,G2frame.HKL)
                    for hkl in G2frame.HKL:
                        hkl.append(2.0*asind(wave/(2.*hkl[3]))+Inst['Zero'][1])             
    rowLabels = []
    for i in range(len(data)): rowLabels.append(str(i+1))
    colLabels = ['position','intensity','use','indexed','h','k','l','d-obs','d-calc']
    Types = [wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,1',wg.GRID_VALUE_BOOL,
        wg.GRID_VALUE_BOOL,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,wg.GRID_VALUE_LONG,
        wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5']
    G2frame.PatternTree.SetItemPyData(IndexId,data)
    G2frame.IndexPeaksTable = G2gd.Table(data,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Index Peak List')
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)                
    G2frame.dataDisplay.SetTable(G2frame.IndexPeaksTable, True)
    for r in range(G2frame.dataDisplay.GetNumberRows()):
        for c in range(G2frame.dataDisplay.GetNumberCols()):
            if c == 2:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, RefreshIndexPeaksGrid)
    G2frame.dataDisplay.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)                 
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataFrame.setSizePosLeft([490,300])
  
################################################################################
#####  Unit cells
################################################################################           
       
def UpdateUnitCellsGrid(G2frame, data):
    '''respond to selection of PWDR Unit Cells data tree item.
    '''
    UnitCellsId = G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
    SPGlist = G2spc.spglist
    bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm',
        'P4/mmm','Fmmm','Immm','Cmmm','Pmmm','C2/m','P2/m','P1']
    spaceGroups = ['F m 3 m','I m 3 m','P m 3 m','R 3 m','P 6/m m m','I 4/m m m',
        'P 4/m m m','F m m m','I m m m','C m m m','P m m m','C 2/m','P 2/m','P -1']
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    wave = G2mth.getWave(Inst)
    
    def SetLattice(controls):
        ibrav = bravaisSymb.index(controls[5])
        if ibrav in [0,1,2]:
            controls[7] = controls[8] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
        elif ibrav in [3,4,5,6]:
            controls[7] = controls[6]
            controls[9] = controls[10] = controls[11] = 90.
            if ibrav in [3,4]:
                controls[11] = 120.
        elif ibrav in [7,8,9,10]:
            controls[9] = controls[10] = controls[11] = 90.
        elif ibrav in [11,12]:
            controls[9] = controls[11] = 90.  # b unique
        if len(controls) < 13: controls.append(0)
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        return ibrav
        
    def OnNcNo(event):
        controls[2] = NcNo.GetValue()
        
    def OnStartVol(event):
        try:
            stVol = int(float(startVol.GetValue()))
            if stVol < 25:
                raise ValueError
        except ValueError:
            stVol = 25
        controls[3] = stVol
        startVol.SetValue("%d"%(stVol))
        
    def OnBravais(event):
        Obj = event.GetEventObject()
        bravais[bravList.index(Obj.GetId())] = Obj.GetValue()
        
    def OnZero(event):
        try:
            Zero = min(5.0,max(-5.0,float(zero.GetValue())))
        except ValueError:
            Zero = 0.0
        controls[1] = Zero
        zero.SetValue("%.4f"%(Zero))
        
    def OnZeroVar(event):
        controls[0] = zeroVar.GetValue()
        
    def OnBravSel(event):
        brav = bravSel.GetString(bravSel.GetSelection())
        controls[5] = brav
        controls[13] = SPGlist[brav][0]       
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def OnSpcSel(event):
        controls[13] = spcSel.GetString(spcSel.GetSelection())       
        
    def OnCellChange(event):
        ibrav = bravaisSymb.index(controls[5])
        Obj = event.GetEventObject()
        ObjId = cellList.index(Obj.GetId())
        try:
            value = max(1.0,float(Obj.GetValue()))
        except ValueError:
            if ObjId < 3:               #bad cell edge - reset
                value = controls[6+ObjId]
            else:                       #bad angle
                value = 90.
        if ibrav in [0,1,2]:
            controls[6] = controls[7] = controls[8] = value
            controls[9] = controls[10] = controls[11] = 90.0
            Obj.SetValue("%.5f"%(controls[6]))
        elif ibrav in [3,4,5,6]:
            if ObjId == 0:
                controls[6] = controls[7] = value
                Obj.SetValue("%.5f"%(controls[6]))
            else:
                controls[8] = value
                Obj.SetValue("%.5f"%(controls[8]))
            controls[9] = controls[10] = controls[11] = 90.0
            if ibrav in [3,4]:
                controls[11] = 120.
        elif ibrav in [7,8,9,10]:
            controls[6+ObjId] = value
            Obj.SetValue("%.5f"%(controls[6+ObjId]))
            controls[9] = controls[10] = controls[11] = 90.0
        elif ibrav in [11,12]:
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
        
    def OnHklShow(event):
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId    
        limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Limits'))[1]
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        spc = controls[13]
        SGData = G2spc.SpcGroup(spc)[1]
        dmin = wave/(2.0*sind(limits[1]/2.0))
        G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A)
        for hkl in G2frame.HKL:
            hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'][1])             
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
            
    def OnSortCells(event):
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        c =  event.GetCol()
        if colLabels[c] == 'M20':
            cells = G2indx.sortM20(cells)
        elif colLabels[c] == 'Volume':
            cells = G2indx.sortVolume(cells)
        else:
            return
        data = [controls,bravais,cells,dmin]
        G2frame.PatternTree.SetItemPyData(UnitCellsId,data)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def CopyUnitCell(event):
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(UnitCellsId)
        for Cell in cells:
            if Cell[-2]:
                break
        cell = Cell[2:9]
        controls[4] = 1
        controls[5] = bravaisSymb[cell[0]]
        controls[6:12] = cell[1:8]
        controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
        controls[13] = spaceGroups[bravaisSymb.index(controls[5])]
        G2frame.PatternTree.SetItemPyData(UnitCellsId,[controls,bravais,cells,dmin])
        G2frame.dataFrame.RefineCell.Enable(True)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)        
                
    def RefineCell(event):
        def cellPrint(ibrav,A):
            cell = G2lat.A2cell(A)
            Vol = G2lat.calc_V(A)
            if ibrav in [0,1,2]:
                print "%s%10.6f" % ('a =',cell[0])
            elif ibrav in [3,4,5,6]:
                print "%s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],' c =',cell[2],' volume =',Vol)
            elif ibrav in [7,8,9,10]:
                print "%s%10.6f %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],' volume =',Vol)
            elif ibrav in [11,12]:
                print "%s%10.6f %s%10.6f %s%10.6f %s%8.3f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],'beta =',cell[4],' volume =',Vol)
            else:
                print "%s%10.6f %s%10.6f %s%10.6f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2])
                print "%s%8.3f %s%8.3f %s%8.3f %s%12.3f" % ('alpha =',cell[3],'beta =',cell[4],'gamma =',cell[5],' volume =',Vol)
             
        PatternId = G2frame.PatternId
        PickId = G2frame.PickId    
        peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Index Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!', 'Nothing to refine!')
            return        
        print 'Refine cell'
        controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        ibrav = bravaisSymb.index(controls[5])
        SGData = G2spc.SpcGroup(controls[13])[1]
        dmin = G2indx.getDmin(peaks)-0.005
        G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A)
        G2indx.IndexPeaks(peaks,G2frame.HKL)
        Lhkl,M20,X20,Aref,Zero = G2indx.refinePeaksZ(peaks,wave,ibrav,A,controls[1],controls[0])            
        controls[1] = Zero
        controls[6:12] = G2lat.A2cell(Aref)
        controls[12] = G2lat.calc_V(Aref)
        data = [controls,bravais,cells,dmin]
        cells = G2frame.PatternTree.GetItemPyData(UnitCellsId)[2]
        for cell in cells:
            cell[-2] = False
        cells.insert(0,[M20,X20,ibrav]+controls[6:13]+[True,False])
        data[2] = cells
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
        print "%s%10.3f" % ('refinement M20 = ',M20)
        print 'unindexed lines = ',X20
        cellPrint(ibrav,Aref)
        for hkl in G2frame.HKL:
            hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'][1])             
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
        
    def IndexPeaks(event):
        PatternId = G2frame.PatternId    
        print 'Peak Indexing'
        keepcells = []
        try:
            controls,bravais,cells,dmin = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'))
            for cell in cells:
                if cell[11]:
                    keepcells.append(cell)
        except IndexError:
            pass
        except ValueError:
            G2frame.ErrorDialog('Error','Need to set controls in Unit Cell List first')
            return
        if True not in bravais:
            G2frame.ErrorDialog('Error','No Bravais lattices selected')
            return
        G2frame.dataFrame.CopyCell.Enable(False)
        G2frame.dataFrame.RefineCell.Enable(False)
        OK,dmin,newcells = G2indx.DoIndexPeaks(peaks,wave,controls,bravais)
        cells = keepcells+newcells
        cells = G2indx.sortM20(cells)
        cells[0][10] = True
        if OK:
            data = [controls,bravais,cells,dmin]
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                G2frame.HKL = G2lat.GenHBravais(dmin,bestCell[2],G2lat.cell2A(bestCell[3:9]))
                for hkl in G2frame.HKL:
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'][1])             
                if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame)
                else:
                    G2plt.PlotPatterns(G2frame)
            G2frame.dataFrame.CopyCell.Enable(True)
            G2frame.dataFrame.IndexPeaks.Enable(True)
            G2frame.dataFrame.MakeNewPhase.Enable(True)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)
                
    def RefreshUnitCellsGrid(event):
        data =G2frame.PatternTree.GetItemPyData(UnitCellsId)
        cells,dmin = data[2:]
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
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'][1])             
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
                    G2IO.SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell[1:]))
                Status.SetStatusText('Change space group from '+str(controls[13])+' if needed')
        finally:
            dlg.Destroy()
            
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.IndexMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.Bind(wx.EVT_MENU, IndexPeaks, id=G2gd.wxID_INDEXPEAKS)
    G2frame.Bind(wx.EVT_MENU, CopyUnitCell, id=G2gd.wxID_COPYCELL)
    G2frame.Bind(wx.EVT_MENU, RefineCell, id=G2gd.wxID_REFINECELL)
    G2frame.Bind(wx.EVT_MENU, MakeNewPhase, id=G2gd.wxID_MAKENEWPHASE)
    
    controls,bravais,cells,dmin = data
    if len(controls) < 13:              #add cell volume if missing
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    if len(controls) < 14:              #add space gropu used in indexing
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
    G2frame.sp = wx.SplitterWindow(G2frame.dataFrame)
    G2frame.dataDisplay = wx.Panel(G2frame.sp, style=wx.SUNKEN_BORDER)
    G2frame.dataFrame.IndexPeaks.Enable(False)
    peaks = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
    if peaks:
        G2frame.dataFrame.IndexPeaks.Enable(True)
    G2frame.dataFrame.RefineCell.Enable(False)
    if controls[12] > 1.0:                               #if a "real" volume (i.e. not default)
        G2frame.dataFrame.RefineCell.Enable(True)    
    G2frame.dataFrame.CopyCell.Enable(False)
    G2frame.dataFrame.MakeNewPhase.Enable(False)        
    if cells:
        G2frame.bottom = wx.Panel(G2frame.sp, style=wx.SUNKEN_BORDER)
        G2frame.sp.SplitHorizontally(G2frame.dataDisplay,G2frame.bottom,0)
        G2frame.dataFrame.CopyCell.Enable(True)
        G2frame.dataFrame.MakeNewPhase.Enable(True)        
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Indexing controls: '),0,WACV)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,5,5,5)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Max Nc/Nobs '),0,WACV)
    NcNo = wx.SpinCtrl(G2frame.dataDisplay)
    NcNo.SetRange(1,6)
    NcNo.SetValue(controls[2])
    NcNo.Bind(wx.EVT_SPINCTRL,OnNcNo)
    littleSizer.Add(NcNo,0,WACV)
    littleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Start Volume '),0,WACV)
    startVol = wx.TextCtrl(G2frame.dataDisplay,value=str('%d'%(controls[3])),style=wx.TE_PROCESS_ENTER)
    startVol.Bind(wx.EVT_TEXT_ENTER,OnStartVol)
    startVol.Bind(wx.EVT_KILL_FOCUS,OnStartVol)
    littleSizer.Add(startVol,0,WACV)
    mainSizer.Add(littleSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Select Bravais Lattices for indexing: '),
        0,WACV)
    mainSizer.Add((5,5),0)
    littleSizer = wx.FlexGridSizer(2,7,5,5)
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
    
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Cell Refinement: '),0,WACV)
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
    littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=" Zero offset"),0,WACV)
    zero = wx.TextCtrl(G2frame.dataDisplay,value="%.4f"%(controls[1]),style=wx.TE_PROCESS_ENTER)
    zero.Bind(wx.EVT_TEXT_ENTER,OnZero)
    zero.Bind(wx.EVT_KILL_FOCUS,OnZero)
    littleSizer.Add(zero,0,WACV)
    zeroVar = wx.CheckBox(G2frame.dataDisplay,label="Refine?")
    zeroVar.SetValue(controls[0])
    zeroVar.Bind(wx.EVT_CHECKBOX,OnZeroVar)
    littleSizer.Add(zeroVar,0,WACV)
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
    littleSizer = wx.FlexGridSizer(2,useGUI[1],5,5)
    for txt,fmt,ifEdit,Id in useGUI[2]:
        littleSizer.Add(wx.StaticText(G2frame.dataDisplay,label=txt),0,WACV)
        if ifEdit:          #a,b,c,etc.
            cellVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[6+Id])),style=wx.TE_PROCESS_ENTER)
            cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
            cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
            littleSizer.Add(cellVal,0,WACV)
            cellList.append(cellVal.GetId())
        else:               #volume
            volVal = wx.TextCtrl(G2frame.dataDisplay,value=(fmt%(controls[12])),style=wx.TE_READONLY)
            volVal.SetBackgroundColour(VERY_LIGHT_GREY)
            littleSizer.Add(volVal,0,WACV)
    mainSizer.Add(littleSizer,0)
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    topSize = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataDisplay.SetSize(topSize)
    if cells:
        if ibrav == 13:
            topSize[1] += 230
        else:
            topSize[1] += 200
    G2frame.dataFrame.setSizePosLeft(topSize)    
    
    if cells:
        bottomSize = topSize        #screwy but bottom doesn't have a size in linux!
        bottomSize[0] -= 20         #to reveal slider
        if ibrav == 13:
            bottomSize[1] -= 240
        else:
            bottomSize[1] -= 210
        wx.StaticText(parent=G2frame.bottom,label=' Indexing Result ')
        rowLabels = []
        colLabels = ['M20','X20','use','Bravais','a','b','c','alpha','beta','gamma','Volume','Keep']
        Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+ \
            [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL]
        numRows = len(cells)
        table = []
        for cell in cells:
            rowLabels.append('')
            row = cell[0:2]+[cell[-2]]+[bravaisSymb[cell[2]]]+cell[3:10]+[cell[11],]
            if cell[-2]:
                A = G2lat.cell2A(cell[3:9])
                G2frame.HKL = G2lat.GenHBravais(dmin,cell[2],A)
                for hkl in G2frame.HKL:
                    hkl.append(2.0*asind(wave/(2.*hkl[3]))+controls[1]+Inst['Zero'][1])             
            table.append(row)
        UnitCellsTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        gridDisplay = G2gd.GSGrid(G2frame.bottom)
        gridDisplay.SetPosition(wx.Point(0,20))                
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
        gridDisplay.SetSize(bottomSize)

################################################################################
#####  Reflection list
################################################################################           
       
def UpdateReflectionGrid(G2frame,data,HKLF=False,Name=''):
    '''respond to selection of PWDR Reflections data tree item.
    '''
    if not data:
        print 'No phases, no reflections'
        return
    if HKLF:
        G2frame.RefList = 1
        phaseName = Name
    else:
        phaseName = G2frame.RefList
        phases = data.keys()
    
        def OnSelectPhase(event):
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Phase',phases)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    G2frame.RefList = phases[sel]
                    UpdateReflectionGrid(G2frame,data)
            finally:
                dlg.Destroy()
            G2plt.PlotPatterns(G2frame)
            
        def OnPlotHKL(event):
            PlotSngl(self,newPlot=True)           #this isn't going to work; no HKL plot controls! 
        
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    rowLabels = []
    if HKLF:
        G2gd.SetDataMenuBar(G2frame)
        refs = data[1]['RefList']
    else:        
        G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ReflMenu)
        if not G2frame.dataFrame.GetStatusBar():
            Status = G2frame.dataFrame.CreateStatusBar()    
        G2frame.Bind(wx.EVT_MENU, OnSelectPhase, id=G2gd.wxID_SELECTPHASE)
        G2frame.Bind(wx.EVT_MENU, OnPlotHKL, id=G2gd.wxID_PWDHKLPLOT)   #not work yet
        G2frame.dataFrame.SelectPhase.Enable(False)
        if len(data) > 1:
            G2frame.dataFrame.SelectPhase.Enable(True)
        try:            #patch for old reflection lists
            refList = np.array(data[G2frame.RefList]['RefList'])
            I100 = refList.T[8]*refList.T[11]
        except TypeError:
            refList = np.array([refl[:11] for refl in data[G2frame.RefList]])
            I100 = refList.T[8]*np.array([refl[13] for refl in data[G2frame.RefList]])
        Imax = np.max(I100)
        if Imax:
            I100 *= 100.0/Imax
        refs = np.vstack((refList.T[:11],I100)).T
    for i in range(len(refs)): rowLabels.append(str(i))
    if HKLF:
        colLabels = ['H','K','L','mul','d','Fosq','sig','Fcsq','FoTsq','FcTsq','phase',]
    else:
        colLabels = ['H','K','L','mul','d','pos','sig','gam','Fosq','Fcsq','phase','I100',]
    Types = 4*[wg.GRID_VALUE_LONG,]+4*[wg.GRID_VALUE_FLOAT+':10,4',]+ \
        2*[wg.GRID_VALUE_FLOAT+':10,2',]+[wg.GRID_VALUE_FLOAT+':10,3',]+ \
        [wg.GRID_VALUE_FLOAT+':10,2',]
    G2frame.PeakTable = G2gd.Table(refs,rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataFrame.SetLabel('Reflection List for '+phaseName)
    G2frame.dataDisplay = G2gd.GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.PeakTable, True)
    G2frame.dataDisplay.EnableEditing(False)
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    G2frame.dataDisplay.Fit()
    size = G2frame.dataDisplay.GetSize()
    G2frame.dataFrame.setSizePosLeft([size[0]+32,350])
    
################################################################################
#####  SASD Substances 
################################################################################
           
def UpdateSubstanceGrid(G2frame,data):
    '''respond to selection of SASD Substance data tree item.
    '''
    import Substances as substFile
    
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
            'Scatt density':0.0,'XAnom density':0.0,'XAbsorption':0.0}
        subst = substFile.Substances[name]
        ElList = subst['Elements'].keys()
        for El in ElList:
            Info = G2elem.GetAtomInfo(El.strip().upper())
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
            data['Substances'][name]['Scatt density'] = \
                G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
            contrst,absorb = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)         
            data['Substances'][name]['XAnom density'] = contrst
            data['Substances'][name]['XAbsorption'] = absorb
                         
        UpdateSubstanceGrid(G2frame,data)
        
    def OnCopySubstance(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2gd.G2MultiChoiceDialog(
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
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Substances'),
                copy.copy(data))
    
    def OnAddSubstance(event):
        dlg = wx.TextEntryDialog(None,'Enter a name for this substance','Substance Name Entry','New substance',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            Name = dlg.GetValue()
            data['Substances'][Name] = {'Elements':{},'Volume':1.0,'Density':1.0,
                'Scatt density':0.0,'XAnom density':0.,'XAbsorption':0.}
        dlg.Destroy()
        AddElement(Name)
        UpdateSubstanceGrid(G2frame,data)
        
    def OnDeleteSubstance(event):
        TextList = []
        for name in data['Substances']:
            if name != 'vacuum':
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
            if name != 'vacuum':
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
        AddElement(name)
        UpdateSubstanceGrid(G2frame,data)
        
    def AddElement(name):
        ElList = data['Substances'][name]['Elements'].keys()
        dlg = G2elemGUI.PickElements(G2frame,ElList)
        if dlg.ShowModal() == wx.ID_OK:
            for El in dlg.Elem:
                El = El.strip().upper()
                Info = G2elem.GetAtomInfo(El)
                Info.update({'Num':1})
                data['Substances'][name]['Elements'][El] = Info
                data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
                data['Substances'][name]['Scatt density'] = \
                    G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                contrst,absorb = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)         
                data['Substances'][name]['XAnom density'] = contrst
                data['Substances'][name]['XAbsorption'] = absorb
        dlg.Destroy()
        
        
    def OnDeleteElement(event):
        TextList = []
        for name in data['Substances']:
            if name != 'vacuum':
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
                data['Substances'][name]['Scatt density'] = \
                    G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
                contrst,absorb = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)         
                data['Substances'][name]['XAnom density'] = contrst
                data['Substances'][name]['XAbsorption'] = absorb
        UpdateSubstanceGrid(G2frame,data)
                
    def SubstSizer():
        
        def OnValueChange(event):
            Obj = event.GetEventObject()
            if len(Indx[Obj.GetId()]) == 3:
                name,El,keyId = Indx[Obj.GetId()]
                try:
                    value = max(0,float(Obj.GetValue()))
                except ValueError:
                    value = 0
                    Obj.SetValue('%.2f'%(value))
                data['Substances'][name]['Elements'][El][keyId] = value
                data['Substances'][name]['Volume'] = G2mth.El2EstVol(data['Substances'][name]['Elements'])
                data['Substances'][name]['Density'] = \
                    G2mth.Vol2Den(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])
            else:
                name,keyId = Indx[Obj.GetId()]
                try:
                    value = max(0,float(Obj.GetValue()))
                except ValueError:
                    value = 1.0
                data['Substances'][name][keyId] = value
                if keyId in 'Volume':
                    data['Substances'][name]['Density'] = \
                        G2mth.Vol2Den(data['Substances'][name]['Elements'],value)
                elif keyId in 'Density':
                    data['Substances'][name]['Volume'] = \
                        G2mth.Den2Vol(data['Substances'][name]['Elements'],value)
            data['Substances'][name]['Scatt density'] = \
                G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'])[0]
            contrst,absorb = G2mth.XScattDen(data['Substances'][name]['Elements'],data['Substances'][name]['Volume'],wave)         
            data['Substances'][name]['XAnom density'] = contrst
            data['Substances'][name]['XAbsorption'] = absorb
            UpdateSubstanceGrid(G2frame,data)
        
        Indx = {}
        substSizer = wx.BoxSizer(wx.VERTICAL)
        substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Substance list: wavelength: %.5fA'%(wave)),
            0,WACV)
        for name in data['Substances']:
            G2gd.HorizontalLine(substSizer,G2frame.dataDisplay)    
            substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Data for '+name+':'),
                0,WACV)
            if name == 'vacuum':
                substSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='        Not applicable'),
                    0,WACV)
            else:    
                elSizer = wx.FlexGridSizer(1,6,5,5)
                Substance = data['Substances'][name]
                Elems = Substance['Elements']
                for El in Elems:    #do elements as pull downs for isotopes for neutrons
                    elSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' '+El+': '),
                        0,WACV)
                    num = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(Elems[El]['Num']),style=wx.TE_PROCESS_ENTER)
                    Indx[num.GetId()] = [name,El,'Num']
                    num.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
                    num.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
                    elSizer.Add(num,0,WACV)
                substSizer.Add(elSizer,0)
                vdsSizer = wx.FlexGridSizer(1,4,5,5)
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Volume: '),
                    0,WACV)
                vol = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(Substance['Volume']),style=wx.TE_PROCESS_ENTER)
                Indx[vol.GetId()] = [name,'Volume']
                vol.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
                vol.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
                vdsSizer.Add(vol,0,WACV)                
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Density: '),
                    0,WACV)
                den = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(Substance['Density']),style=wx.TE_PROCESS_ENTER)
                Indx[den.GetId()] = [name,'Density']
                den.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
                den.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
                vdsSizer.Add(den,0,WACV)
                substSizer.Add(vdsSizer,0)
                substSizer.Add(wx.StaticText(G2frame.dataDisplay,
                    label=' Scattering density  : %.2f *10%scm%s'%(Substance['Scatt density'],Pwr10,Pwrm2)),
                    0,WACV)                
                substSizer.Add(wx.StaticText(G2frame.dataDisplay,       #allow neutrons here into NAnom density & NAbsorption
                    label=' Anomalous density : %.2f *10%scm%s'%(Substance['XAnom density'],Pwr10,Pwrm2)),
                    0,WACV)                
                substSizer.Add(wx.StaticText(G2frame.dataDisplay,
                    label=' X-ray absorption   : %.2f cm%s'%(Substance['XAbsorption'],Pwrm1)),
                    0,WACV)                
        return substSizer
            
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    wave = G2mth.getWave(Inst)
    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.SubstanceMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2frame.dataFrame.SetLabel('Substances')
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadSubstance, id=G2gd.wxID_LOADSUBSTANCE)    
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
    G2frame.dataDisplay.SetSize(Size)
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
                    'FFargs':{},'NumPoints':50,'Cutoff':0.01,'AutoDist':True,'logR':True,
                    'StrFact':'Dilute'},
                'LogNormal':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[0.5,False],'MinSize':[10.,False],},
                'Gaussian':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[300.,False],},
                'LSW':{'Volume':[0.05,False],'Mean':[1000.0,False],},
                'Schulz-Zimm':{'Volume':[0.05,False],'Mean':[1000.,False],'StdDev':[300.,False],},
                'Unified':{'G':[100,False],'Rg':[100,False],'B':[1.e-4,False],'P':[4,False],'Cutoff':[1e-5,False],},
                'Porod':{'B':[1.e-4,False],'P':[4,False],'Cutoff':[1e-5,False],},
                'Monodisperse':{'Volume':[0.05,False],'Radius':[100,False],},   #OK for spheres
                'Bragg':{'PkInt':[100,False],'PkPos':[0.2,False],
                    'PkSig':[10,False],'PkGam':[10,False],},        #reeasonable 31A peak
                })
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots(True)
                    
        wx.CallAfter(UpdateModelsGrid,G2frame,data)
        
    def OnCopyModel(event):
        hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
            return
        copyList = []
        dlg = G2gd.G2MultiChoiceDialog(
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
            G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Models'),
                copy.deepcopy(data))
        
    def OnFitModel(event):
        if data['Current'] == 'Size dist.':
            if not any(Sample['Contrast']):
                G2frame.ErrorDialog('No contrast; your sample is a vacuum!',
                    'You need to define a scattering substance!\n'+    \
                    ' Do Substances and then Sample parameters')
                return
            G2sasd.SizeDistribution(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots(True)
            
        elif data['Current'] == 'Particle fit':
            SaveState()
            if not G2sasd.ModelFit(Profile,ProfDict,Limits,Substances,Sample,data):
                G2frame.ErrorDialog('Failed refinement',
                    'You need to rethink your selection of parameters\n'+    \
                    ' Model restored to previous version')
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots(True)
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            
    def OnUnDo(event):
        DoUnDo()
        data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,
            G2frame.PatternId,'Models'))
        G2frame.dataFrame.SasdUndo.Enable(False)
        UpdateModelsGrid(G2frame,data)
        G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
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
        
    def OnValueChange(event):
        Obj = event.GetEventObject()
        itemKey,ind,fmt = Indx[Obj.GetId()]
        try:
            value = float(Obj.GetValue())
        except ValueError:
            value = 0.0
        Obj.SetValue(fmt%(value))
        data[itemKey][ind] = value
        if itemKey == 'Back':
            Profile[4][:] = value
        RefreshPlots()
        
    def OnCheckBox(event):
        Obj = event.GetEventObject()
        item,ind = Indx[Obj.GetId()]
        item[ind] = Obj.GetValue()
        
    def OnIntVal(event):
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
            
        def OnPartVal(event):
            try:
                val = max(0.0,float(partprm.GetValue()))
            except ValueError:
                val = 1
            data['Size']['Shape'][1] = val
            partprm.SetValue('%.3f'%(val))
            
        sizeSizer = wx.BoxSizer(wx.VERTICAL)
        sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Size distribution parameters: '),0,WACV)
        binSizer = wx.FlexGridSizer(1,7,5,5)
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
            'Unified disk':' Thickness: '}
        partsh = wx.ComboBox(G2frame.dataDisplay,value=str(data['Size']['Shape'][0]),choices=shapes.keys(),
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        partsh.Bind(wx.EVT_COMBOBOX,OnShape)        
        partSizer.Add(partsh,0,WACV)
        if data['Size']['Shape'][0] not in ['Unified sphere',]:
            partSizer.Add(wx.StaticText(G2frame.dataDisplay,label=shapes[data['Size']['Shape'][0]]),0,WACV)
            partprm = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['Size']['Shape'][1]),
                style=wx.TE_PROCESS_ENTER)
            partprm.Bind(wx.EVT_TEXT_ENTER,OnPartVal)        
            partprm.Bind(wx.EVT_KILL_FOCUS,OnPartVal)
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
            'Unified tube':{'Length':[100.,False],'Thickness':[10.,False]},}
                
        ffDistChoices =  ['Sphere','Spheroid','Cylinder','Cylinder diam',
            'Cylinder AR','Unified sphere','Unified rod','Unified rod AR',
            'Unified disk','Unified tube',]
                
        ffMonoChoices = ['Sphere','Spheroid','Cylinder','Cylinder AR',
            ]
                  
        def OnValue(event):
            Obj = event.GetEventObject()
            item,key,id,sldrObj = Indx[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
                if value <= 0.:
                    raise ValueError
            except ValueError:
                value = item[key][id]
            item[key][id] = value
            Obj.SetValue('%.3g'%(value))
            if key == 'P':
                sldrObj.SetValue(1000.*value)
            else:
                sldrObj.SetRange(1000.*(np.log10(value)-2),1000.*(np.log10(value)+2))
                sldrObj.SetValue(1000.*np.log10(value))
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots()
            
        def OnSelect(event):
            Obj = event.GetEventObject()
            item,key = Indx[Obj.GetId()]
            item[key] = Obj.GetValue()
            if 'Refine' not in Obj.GetLabel():
                if 'FormFact' in key :
                    item['FFargs'] = FormFactors[Obj.GetValue()]
                wx.CallAfter(UpdateModelsGrid,G2frame,data)
                G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
                RefreshPlots()
                
        def OnDelLevel(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            del data['Particle']['Levels'][item]
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots()
            
        def OnParmSlider(event):
            Obj = event.GetEventObject()
            item,key,id,pvObj = Indx[Obj.GetId()]
            slide = Obj.GetValue()
            if key == 'P':
                value = float(slide/1000.)
            else:
                value = 10.**float(slide/1000.)
            item[key][id] = value
            pvObj.SetValue('%.3g'%(item[key][id]))
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Substances,Sample,data)
            RefreshPlots()
            
        def OnCenterSlider(event):
            Obj = event.GetEventObject()
            Obj.SetValue(False)
            sldrObj,value = Indx[Obj.GetId()]
            sldrObj.SetRange(1000*(value-1),1000*(value+1))
            sldrObj.SetValue(1000*value)
            
        def SizeSizer():
            sizeSizer = wx.FlexGridSizer(0,4,5,5)
            sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Distribution: '),0,WACV)
            Distchoice = ['LogNormal','Gaussian','LSW','Schulz-Zimm','Bragg','Unified','Porod','Monodisperse',]
            distChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['DistType'],choices=Distchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[distChoice.GetId()] = [level['Controls'],'DistType']
            distChoice.Bind(wx.EVT_COMBOBOX,OnSelect)
            sizeSizer.Add(distChoice,0,WACV)
            if level['Controls']['DistType'] not in ['Bragg','Unified','Porod',]:
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Form Factor: '),0,WACV)
                if 'Mono' not in level['Controls']['DistType']:
                    ffChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['FormFact'],choices=ffDistChoices,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                else:
                    ffChoice = wx.ComboBox(G2frame.dataDisplay,value=level['Controls']['FormFact'],choices=ffMonoChoices,
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
                contrast = rho**2-rhoMat**2
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Resonant X-ray contrast: '),0,WACV)
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label='  %.2f 10%scm%s'%(contrast,Pwr20,Pwrm4)),0,WACV)
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Num. radii: '),0,WACV)
                radii = ['25','50','75','100','200']
                nRadii = wx.ComboBox(G2frame.dataDisplay,value=str(level['Controls']['NumPoints']),choices=radii,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[nRadii.GetId()] = [level['Controls'],'NumPoints']
                nRadii.Bind(wx.EVT_COMBOBOX,OnSelect)
                sizeSizer.Add(nRadii,0,WACV)
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' R dist. cutoff: '),0,WACV)
                rCutoff = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,level['Controls'],'Cutoff',
                    min=0.001,max=0.1,typeHint=float)
                sizeSizer.Add(rCutoff,0,WACV)
            elif level['Controls']['DistType']  in ['Unified',]:
                Parms = level['Unified']
                Best = G2sasd.Bestimate(Parms['G'][0],Parms['Rg'][0],Parms['P'][0])
                sizeSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Estimated Dist B: %12.4g'%(Best)),0,WACV)
            return sizeSizer
            
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
        volfrac = G2gd.ValidatedTxtCtrl(G2frame.dataDisplay,data['Particle']['Matrix']['VolFrac'],0,
                typeHint=float)
        topSizer.Add(volfrac,0,WACV)
        volVar = wx.CheckBox(G2frame.dataDisplay,label=' Refine?')
        volVar.SetValue(data['Particle']['Matrix']['VolFrac'][1])
        Indx[volVar.GetId()] = [data['Particle']['Matrix']['VolFrac'],1]
        volVar.Bind(wx.EVT_CHECKBOX, OnSelect)
        topSizer.Add(volVar,0,WACV)
        partSizer.Add(topSizer,0,)
        for ilev,level in enumerate(data['Particle']['Levels']):
            G2gd.HorizontalLine(partSizer,G2frame.dataDisplay)
            topLevel = wx.BoxSizer(wx.HORIZONTAL)
            topLevel.Add(wx.StaticText(G2frame.dataDisplay,label=' Model component %d: '%(ilev)),0,WACV)
            delBtn = wx.Button(G2frame.dataDisplay,label=' Delete?')
            Indx[delBtn.GetId()] = ilev
            delBtn.Bind(wx.EVT_BUTTON,OnDelLevel)
            topLevel.Add(delBtn,0,WACV)
            partSizer.Add(topLevel,0)
            partSizer.Add(SizeSizer())
            
            lvlSizer = wx.BoxSizer(wx.HORIZONTAL)
            parmSizer = wx.FlexGridSizer(1,4,5,5)
            parmSizer.AddGrowableCol(2,1)
            parmSizer.SetFlexibleDirection(wx.HORIZONTAL)
            Parms = level[level['Controls']['DistType']]
            FFargs = level['Controls']['FFargs']
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
                    parmSldr = wx.Slider(G2frame.dataDisplay,minValue=1000.*valMinMax[0],
                        maxValue=1000.*valMinMax[1],value=1000.*value)
                    Indx[parmValue.GetId()] = [Parms,parm,0,parmSldr]
                    Indx[parmSldr.GetId()] = [Parms,parm,0,parmValue]
                    parmSldr.Bind(wx.EVT_SLIDER,OnParmSlider)
                    parmSizer.Add(parmSldr,1,wx.EXPAND)
                    if parm == 'P':
                        parmSizer.Add((5,5),)
                    else:
                        center = wx.CheckBox(G2frame.dataDisplay,label='Center? ')
                        Indx[center.GetId()] = [parmSldr,value]
                        center.Bind(wx.EVT_CHECKBOX,OnCenterSlider)
                        parmSizer.Add(center,0,WACV)
            if level['Controls']['DistType'] not in ['Bragg']:
                for parm in list(FFargs):
                    parmVar = wx.CheckBox(G2frame.dataDisplay,label='Refine? FF '+parm) 
                    parmVar.SetValue(FFargs[parm][1])
                    Indx[parmVar.GetId()] = [FFargs[parm],1]
                    parmVar.Bind(wx.EVT_CHECKBOX, OnSelect)
                    parmSizer.Add(parmVar,0,WACV)
                    parmValue = wx.TextCtrl(G2frame.dataDisplay,value='%.3g'%(FFargs[parm][0]),
                        style=wx.TE_PROCESS_ENTER)
                    parmValue.Bind(wx.EVT_TEXT_ENTER,OnValue)        
                    parmValue.Bind(wx.EVT_KILL_FOCUS,OnValue)
                    parmSizer.Add(parmValue,0,WACV)
                    value = np.log10(FFargs[parm][0])
                    valMinMax = [value-1,value+1]
                    parmSldr = wx.Slider(G2frame.dataDisplay,minValue=1000.*valMinMax[0],
                        maxValue=1000.*valMinMax[1],value=1000.*value)
                    Indx[parmVar.GetId()] = [FFargs[parm],1]
                    Indx[parmValue.GetId()] = [FFargs[parm],0,parmSldr]
                    Indx[parmSldr.GetId()] = [FFargs,parm,0,parmValue]
                    parmSldr.Bind(wx.EVT_SLIDER,OnParmSlider)
                    parmSizer.Add(parmSldr,1,wx.EXPAND)
                    center = wx.CheckBox(G2frame.dataDisplay,label='Center? ')
                    Indx[center.GetId()] = [parmSldr,value]
                    center.Bind(wx.EVT_CHECKBOX,OnCenterSlider)
                    parmSizer.Add(center,0,WACV)
            lvlSizer.Add(parmSizer,1,wx.EXPAND)
            partSizer.Add(lvlSizer,1,wx.EXPAND)
        return partSizer
        
    def OnEsdScale(event):
        try:
            value = float(esdScale.GetValue())
            if value <= 0.:
                raise ValueError
        except ValueError:
            value = 1./np.sqrt(ProfDict['wtFactor'])
        ProfDict['wtFactor'] = 1./value**2
        esdScale.SetValue('%.3f'%(value))
        RefreshPlots(True)
        
    Sample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Sample Parameters'))
    Limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Limits'))
    Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
    Substances = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.PatternId, 'Substances'))
    ProfDict,Profile = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)[:2]

    if G2frame.dataDisplay:
        G2frame.dataFrame.DestroyChildren()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.ModelMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.SetLabel('Modelling')
    G2frame.dataDisplay = wxscroll.ScrolledPanel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyModel, id=G2gd.wxID_MODELCOPY)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFitModel, id=G2gd.wxID_MODELFIT)
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
    esdScale = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(1./np.sqrt(ProfDict['wtFactor'])),style=wx.TE_PROCESS_ENTER)
    esdScale.Bind(wx.EVT_TEXT_ENTER,OnEsdScale)        
    esdScale.Bind(wx.EVT_KILL_FOCUS,OnEsdScale)
    topSizer.Add(esdScale,0,WACV)
    mainSizer.Add(topSizer)
    G2gd.HorizontalLine(mainSizer,G2frame.dataDisplay)
    if 'Size' in data['Current']:
        if 'MaxEnt' in data['Size']['Method']:
            Status.SetStatusText('Size distribution by Maximum entropy')
        elif 'IPG' in data['Size']['Method']:
            Status.SetStatusText('Size distribution by Interior-Point Gradient')
        mainSizer.Add(SizeSizer())        
    elif 'Particle' in data['Current']:
        mainSizer.Add(PartSizer(),1,wx.ALIGN_LEFT|wx.EXPAND)
    G2gd.HorizontalLine(mainSizer,G2frame.dataDisplay)    
    backSizer = wx.BoxSizer(wx.HORIZONTAL)
    backSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Background:'),0,WACV)
    backVal = wx.TextCtrl(G2frame.dataDisplay,value='%.3g'%(data['Back'][0]),style=wx.TE_PROCESS_ENTER)
    Indx[backVal.GetId()] = ['Back',0,'%.3g']
    backVal.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
    backVal.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
    backSizer.Add(backVal,0,WACV)
    backVar = wx.CheckBox(G2frame.dataDisplay,label='Refine?')
    Indx[backVar.GetId()] = [data['Back'],1]
    backVar.SetValue(data['Back'][1])
    backVar.Bind(wx.EVT_CHECKBOX, OnCheckBox)
    backSizer.Add(backVar,0,WACV)
    mainSizer.Add(backSizer)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetAutoLayout(1)
    G2frame.dataDisplay.SetupScrolling()
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[0] += 25
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)    
        
    
################################################################################
#####  PDF controls
################################################################################           
       
def UpdatePDFGrid(G2frame,data):
    '''respond to selection of PWDR PDF data tree item.
    '''
    global inst
    tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
    dataFile = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    powName = 'PWDR'+dataFile[4:]
    powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, powName)
    fullLimits,limits = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Limits'))[:2]
    inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId, 'Instrument Parameters'))[0]
    if 'Lam' in inst:
        keV = 12.397639/inst['Lam'][1]
    else:
        keV = 12.397639/inst['Lam1'][0]
    wave = 12.397639/keV
    qLimits = [tth2q(fullLimits[0],wave),tth2q(fullLimits[1],wave)]
    data['QScaleLim'][1] = min(qLimits[1],data['QScaleLim'][1])
    if data['QScaleLim'][0]:
        data['QScaleLim'][0] = max(qLimits[0],data['QScaleLim'][0])
    else:                                #initial setting at 90% of max Q
        data['QScaleLim'][0] = 0.90*data['QScaleLim'][1]
    polariz = inst['Polariz.'][1]
    azimuth = inst['Azimuth'][1]
    itemDict = {}
    
    def FillFileSizer(fileSizer,key):
        #fileSizer is a FlexGridSizer(3,6)
        
        def OnSelectFile(event):
            Obj = event.GetEventObject()
            fileKey,itemKey,fmt = itemDict[Obj.GetId()]
            if itemKey == 'Name':
                value = Obj.GetValue()
            Obj.SetValue(fmt%(value))
            data[fileKey][itemKey] = value
            UpdatePDFGrid(G2frame,data)
        
        def OnValueChange(event):
            Obj = event.GetEventObject()
            fileKey,itemKey,fmt = itemDict[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = -1.0
            Obj.SetValue(fmt%(value))
            data[fileKey][itemKey] = value
            auxPlot = ComputePDF(data)
            G2plt.PlotISFG(G2frame,newPlot=True)
                        
        item = data[key]
        fileList = np.array(GetFileList('PWDR')).T[1]
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' '+key+' file:'),0,WACV)
        fileName = wx.ComboBox(G2frame.dataDisplay,value=item['Name'],choices=fileList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        itemDict[fileName.GetId()] = [key,'Name','%s']
        fileName.Bind(wx.EVT_COMBOBOX,OnSelectFile)        
        fileSizer.Add(fileName,0,)
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Multiplier:'),0,WACV)
        mult = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(item['Mult']),style=wx.TE_PROCESS_ENTER)
        itemDict[mult.GetId()] = [key,'Mult','%.3f']
        mult.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
        mult.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
        fileSizer.Add(mult,0,)
        fileSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Add:'),0,WACV)
        add = wx.TextCtrl(G2frame.dataDisplay,value='%.0f'%(item['Add']),style=wx.TE_PROCESS_ENTER)
        itemDict[add.GetId()] = [key,'Add','%.0f']
        add.Bind(wx.EVT_TEXT_ENTER,OnValueChange)        
        add.Bind(wx.EVT_KILL_FOCUS,OnValueChange)
        fileSizer.Add(add,0,)
        
    def SumElementVolumes():
        sumVol = 0.
        ElList = data['ElList']
        for El in ElList:
            Avol = (4.*math.pi/3.)*ElList[El]['Drad']**3
            sumVol += Avol*ElList[El]['FormulaNo']
        return sumVol
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def FillElemSizer(elemSizer,ElData):
        
        def OnFractionChange(event):
            try:
                value = max(0.0,float(num.GetValue()))
            except ValueError:
                value = 0.0
            num.SetValue('%.3f'%(value))
            ElData['FormulaNo'] = value
            data['Form Vol'] = max(10.0,SumElementVolumes())
            formVol.SetValue('%.2f'%(data['Form Vol']))
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            auxPlot = ComputePDF(data)
            G2plt.PlotISFG(G2frame,newPlot=True)        
        
        elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
            label=' Element: '+'%2s'%(ElData['Symbol'])+' * '),0,WACV)
        num = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(ElData['FormulaNo']),style=wx.TE_PROCESS_ENTER)
        num.Bind(wx.EVT_TEXT_ENTER,OnFractionChange)        
        num.Bind(wx.EVT_KILL_FOCUS,OnFractionChange)
        elemSizer.Add(num,0,WACV)
        elemSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,
            label="f': %.3f"%(ElData['fp'])+' f": %.3f'%(ElData['fpp'])+' mu: %.2f barns'%(ElData['mu']) ),
            0,WACV)
            
    def OnGeometry(event):
        data['Geometry'] = geometry.GetValue()
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def OnDetType(event):
        data['DetType'] = detType.GetValue()
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def OnFormVol(event):
        try:
            value = float(formVol.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Form Vol']
        data['Form Vol'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
        
    def OnDiameter(event):
        try:
            value = float(diam.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Diam']
        data['Diam'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnPolaVal(event):
        try:
            value = float(polaVal.GetValue())
            if not (0.0 <= value <= 1.0):
                raise ValueError
        except ValueError:
            value = inst['Polariz.'][1]
        inst['Polariz.'][1] = value
        polaVal.SetValue('%.2f'%(inst['Polariz.'][1]))
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
                
    def OnAzimVal(event):
        try:
            value = float(azimVal.GetValue())
            if not (0. <= value <= 360.):
                raise ValueError
        except ValueError:
            value = inst['Azimuth'][1]
        inst['Azimuth'][1] = value
        azimVal.SetValue('%.1f'%(inst['Azimuth'][1]))
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
                        
    def OnObliqCoeff(event):
        try:
            value = float(obliqCoeff.GetValue())
            if value < 0.0:
                raise ValueError
            elif value > 1.0:
                value = 1.0
        except ValueError:
            value = data['ObliqCoeff']
        data['ObliqCoeff'] = value
        obliqCoeff.SetValue('%.3f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnRulandWdt(event):
        try:
            value = float(rulandWdt.GetValue())
            if value <= 0.001:
                raise ValueError
            elif value > 1.0:
                value = 1.0
        except ValueError:
            value = data['Ruland']
        data['Ruland'] = value
        rulandWdt.SetValue('%.3f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnRulSlider(event):
        value = int(rulandSldr.GetValue())/1000.
        data['Ruland'] = max(0.001,value)
        rulandWdt.SetValue('%.3f'%(data['Ruland']))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)
        
    def OnLorch(event):
        data['Lorch'] = lorch.GetValue()
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
                        
    def OnPacking(event):
        try:
            value = float(pack.GetValue())
            if value <= 0.0:
                raise ValueError
        except ValueError:
            value = data['Pack']
        data['Pack'] = value
        UpdatePDFGrid(G2frame,data)
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=False)        
                
    def OnSQmin(event):
        try:
            value = float(SQmin.GetValue())
            if value < qLimits[0]:
                raise ValueError
        except ValueError:
            value = max(qLimits[0],data['QScaleLim'][0])
        data['QScaleLim'][0] = value
        SQmin.SetValue('%.1f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        
        
    def OnSQmax(event):
        try:
            value = float(SQmax.GetValue())
            if value > qLimits[1]:
                raise ValueError
        except ValueError:
            value = min(qLimits[1],data['QScaleLim'][1])
        data['QScaleLim'][1] = value
        if value < data['QScaleLim'][0]:
            data['QScaleLim'][0] = 0.90*value
            SQmin.SetValue('%.1f'%(data['QScaleLim'][0]))
        SQmax.SetValue('%.1f'%(value))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)
        
    def OnResetQ(event):
        resetQ.SetValue(False)
        data['QScaleLim'][1] = qLimits[1]
        SQmax.SetValue('%.1f'%(data['QScaleLim'][1]))
        data['QScaleLim'][0] = 0.9*qLimits[1]
        SQmin.SetValue('%.1f'%(data['QScaleLim'][0]))
        auxPlot = ComputePDF(data)
        G2plt.PlotISFG(G2frame,newPlot=True)        

    def GetFileList(fileType,skip=None):
        fileList = [[False,'',0]]
        Source = ''
        id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
        while id:
            name = G2frame.PatternTree.GetItemText(id)
            if fileType in name:
                if id == skip:
                    Source = name
                else:
                    fileList.append([False,name,id])
            id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
        
    def OnCopyPDFControls(event):
        import copy
        TextList,Source = GetFileList('PDF',skip=G2frame.PatternId)
        TextList[0] = [False,'All PDF',0]
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "PDF" pattern')
            return
        dlg = G2frame.CopyDialog(G2frame,'Copy PDF controls','Copy controls from '+Source+' to:',TextList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                if result[0][0]:
                    result = TextList[1:]
                    for item in result: item[0] = True
                for i,item in enumerate(result):
                    ifcopy,name,id = item
                    if ifcopy:
                        olddata = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'))
                        sample = olddata['Sample']
                        olddata.update(copy.deepcopy(data))
                        olddata['Sample'] = sample
                        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id, 'PDF Controls'),olddata)
                Status.SetStatusText('PDF controls copied')
        finally:
            dlg.Destroy()
                
    def OnSavePDFControls(event):
        print 'save PDF controls?'
        
    def OnLoadPDFControls(event):
        print 'Load PDF controls?'
        
    def OnAddElement(event):
        ElList = data['ElList']
        PE = G2elemGUI.PickElement(G2frame,oneOnly=True)
        if PE.ShowModal() == wx.ID_OK:
            El = PE.Elem
            if El not in ElList and El != 'None':
                ElemSym = El.strip().upper()                
                FpMu = G2elem.FPcalc(G2elem.GetXsectionCoeff(ElemSym), keV)
                ElData = G2elem.GetFormFactorCoeff(ElemSym)[0]
                ElData['FormulaNo'] = 0.0
                ElData.update(G2elem.GetAtomInfo(ElemSym))
                ElData.update(dict(zip(['fp','fpp','mu'],FpMu)))
                ElData.update(G2elem.GetFFC5(El))
                data['ElList'][El] = ElData
            data['Form Vol'] = max(10.0,SumElementVolumes())
        PE.Destroy()
        UpdatePDFGrid(G2frame,data)
        
    def OnDeleteElement(event):
        ElList = data['ElList']
        choice = ElList.keys()
        dlg = G2elemGUI.DeleteElement(G2frame,choice=choice)
        if dlg.ShowModal() == wx.ID_OK:
            del ElList[dlg.GetDeleteElement()]
        dlg.Destroy()
        UpdatePDFGrid(G2frame,data)
                
    def ComputePDF(Data):
        xydata = {}
        for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
            name = Data[key]['Name']
            if name:
                xydata[key] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root,name))
                PDFname = name
        powName = xydata['Sample'][2]
        powId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,powName)
        inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,powId,'Instrument Parameters'))[0]
        auxPlot = G2pwd.CalcPDF(Data,inst,xydata)
        PDFId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'PDF '+powName[4:])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'I(Q)'+powName[4:]),xydata['IofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'S(Q)'+powName[4:]),xydata['SofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'F(Q)'+powName[4:]),xydata['FofQ'])
        G2frame.PatternTree.SetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PDFId,'G(R)'+powName[4:]),xydata['GofR'])
        return auxPlot
        
    def OnComputePDF(event):
        print 'Calculating PDF:'
        auxPlot = ComputePDF(data)
        print 'Done calculating PDF:'
        Status.SetStatusText('PDF computed')
        for plot in auxPlot:
            G2plt.PlotXY(G2frame,plot[:2],type=plot[2])
        
        G2plt.PlotISFG(G2frame,newPlot=True,type='I(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='S(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='F(Q)')
        G2plt.PlotISFG(G2frame,newPlot=True,type='G(R)')
        
    def OnComputeAllPDF(event):
        print 'Calculating PDFs:'
        if G2frame.PatternTree.GetCount():
            id, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while id:
                Name = G2frame.PatternTree.GetItemText(id)
                if 'PDF' in Name:
                    Data = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,id,'PDF Controls'))
                    auxPlot = ComputePDF(Data)                    
                id, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            Status.SetStatusText('All PDFs computed')
            G2plt.PlotISFG(G2frame,newPlot=True,type='G(R)')
            print ' Done calculating PDFs:'
        
    def OnShowTip(G2frame,tip):
        print tip    
                
    if G2frame.dataDisplay:
        G2frame.dataFrame.Clear()
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PDFMenu)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()    
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnCopyPDFControls, id=G2gd.wxID_PDFCOPYCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSavePDFControls, id=G2gd.wxID_PDFSAVECONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnLoadPDFControls, id=G2gd.wxID_PDFLOADCONTROLS)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddElement, id=G2gd.wxID_PDFADDELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2gd.wxID_PDFDELELEMENT)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputePDF, id=G2gd.wxID_PDFCOMPUTE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnComputeAllPDF, id=G2gd.wxID_PDFCOMPUTEALL)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' PDF data files: '),0,WACV)
    mainSizer.Add((5,5),0)
    str = ' Sample file: PWDR %s   Wavelength, A: %.5f  Energy, keV: %.3f  Polariz.: %.2f '%(dataFile[3:],wave,keV,polariz)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=str),0,WACV)
#    dataSizer = wx.BoxSizer(wx.HORIZONTAL)
#    dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Azimuth'),0,WACV)
#    azimVal = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(inst['Azimuth']))
#    azimVal.Bind(wx.EVT_TEXT_ENTER,OnAzimVal)        
#    azimVal.Bind(wx.EVT_KILL_FOCUS,OnAzimVal)
#    dataSizer.Add(azimVal,0)    
#    dataSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label='Polarization'),0,WACV)
#    polaVal = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(inst['Polariz.']))
#    polaVal.Bind(wx.EVT_TEXT_ENTER,OnPolaVal)        
#    polaVal.Bind(wx.EVT_KILL_FOCUS,OnPolaVal)
#    dataSizer.Add(polaVal,0)    
#    mainSizer.Add(dataSizer,0)
    mainSizer.Add((5,5),0)
    fileSizer = wx.FlexGridSizer(3,6,5,1)
    select = ['Sample Bkg.','Container']
    if data['Container']['Name']:
        select.append('Container Bkg.')
    for key in select:
        FillFileSizer(fileSizer,key)
    mainSizer.Add(fileSizer,0)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample information: '),0,WACV)
    mainSizer.Add((5,5),0)    

    ElList = data['ElList']
    Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
    Trans = G2pwd.Transmission(data['Geometry'],Abs*data['Pack'],data['Diam'])
    elemSizer = wx.FlexGridSizer(3,3,5,1)
    for El in ElList:
        FillElemSizer(elemSizer,ElList[El])
    mainSizer.Add(elemSizer,0)
    mainSizer.Add((5,5),0)    
    midSizer = wx.BoxSizer(wx.HORIZONTAL)
    midSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Formula volume: '),0,WACV)
    formVol = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(data['Form Vol']))
    formVol.Bind(wx.EVT_TEXT_ENTER,OnFormVol)        
    formVol.Bind(wx.EVT_KILL_FOCUS,OnFormVol)
    midSizer.Add(formVol,0)
    midSizer.Add(wx.StaticText(G2frame.dataDisplay,
        label=' Theoretical absorption: %.4f cm-1 Sample absorption: %.4f cm-1'%(Abs,Abs*data['Pack'])),
        0,WACV)
    mainSizer.Add(midSizer,0)
    mainSizer.Add((5,5),0)    

    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample geometry: '),0,WACV)
    choice = ['Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate']
    geometry = wx.ComboBox(G2frame.dataDisplay,value=data['Geometry'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    geometry.Bind(wx.EVT_COMBOBOX, OnGeometry)
    geoBox.Add(geometry,0)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample diameter/thickness, mm: '),0,WACV)
    diam = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['Diam']))
    diam.Bind(wx.EVT_TEXT_ENTER,OnDiameter)        
    diam.Bind(wx.EVT_KILL_FOCUS,OnDiameter)
#    diam.Bind(wx.EVT_SET_FOCUS,OnShowTip(G2frame,'tip')) #this doesn't work - what would????
    geoBox.Add(diam,0)
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
    geoBox = wx.BoxSizer(wx.HORIZONTAL)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Packing: '),0,WACV)
    pack = wx.TextCtrl(G2frame.dataDisplay,value='%.2f'%(data['Pack']))
    pack.Bind(wx.EVT_TEXT_ENTER,OnPacking)        
    pack.Bind(wx.EVT_KILL_FOCUS,OnPacking)
    geoBox.Add(pack,0)
    geoBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Sample transmission: %.3f %%'%(Trans)),0,WACV)    
    mainSizer.Add(geoBox,0)
    mainSizer.Add((5,5),0)    
        
    mainSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' S(Q)->F(Q)->G(R) controls: '),0,WACV)
    mainSizer.Add((5,5),0)
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Detector type: '),0,WACV)
    choice = ['Image plate','Point detector']
    detType = wx.ComboBox(G2frame.dataDisplay,value=data['DetType'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
    detType.Bind(wx.EVT_COMBOBOX, OnDetType)
    sqBox.Add(detType,0)
    if data['DetType'] == 'Image plate':
        sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' IP transmission coeff.: '),0,WACV)
        obliqCoeff = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['ObliqCoeff']))
        obliqCoeff.Bind(wx.EVT_TEXT_ENTER,OnObliqCoeff)        
        obliqCoeff.Bind(wx.EVT_KILL_FOCUS,OnObliqCoeff)
        sqBox.Add(obliqCoeff,0)
    mainSizer.Add(sqBox,0)
        
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Ruland width: '),0,WACV)    
    rulandSldr = wx.Slider(parent=G2frame.dataDisplay,style=wx.SL_HORIZONTAL,
        value=int(1000*data['Ruland']))
    sqBox.Add(rulandSldr,1,wx.EXPAND)
    rulandSldr.Bind(wx.EVT_SLIDER, OnRulSlider)
    rulandWdt = wx.TextCtrl(G2frame.dataDisplay,value='%.3f'%(data['Ruland']))
    rulandWdt.Bind(wx.EVT_TEXT_ENTER,OnRulandWdt)        
    rulandWdt.Bind(wx.EVT_KILL_FOCUS,OnRulandWdt)
    sqBox.Add(rulandWdt,0,WACV)    
    mainSizer.Add(sqBox,0,wx.ALIGN_LEFT|wx.EXPAND)
    
    sqBox = wx.BoxSizer(wx.HORIZONTAL)
    lorch = wx.CheckBox(parent=G2frame.dataDisplay,label='Lorch damping?')
    lorch.SetValue(data['Lorch'])
    lorch.Bind(wx.EVT_CHECKBOX, OnLorch)
    sqBox.Add(lorch,0,WACV)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' Scaling q-range: '),0,WACV)
    SQmin = wx.TextCtrl(G2frame.dataDisplay,value='%.1f'%(data['QScaleLim'][0]))
    SQmin.Bind(wx.EVT_TEXT_ENTER,OnSQmin)        
    SQmin.Bind(wx.EVT_KILL_FOCUS,OnSQmin)    
    sqBox.Add(SQmin,0)
    sqBox.Add(wx.StaticText(G2frame.dataDisplay,label=' to '),0,WACV)
    SQmax = wx.TextCtrl(G2frame.dataDisplay,value='%.1f'%(data['QScaleLim'][1]))
    SQmax.Bind(wx.EVT_TEXT_ENTER,OnSQmax)        
    SQmax.Bind(wx.EVT_KILL_FOCUS,OnSQmax)
    sqBox.Add(SQmax,0)
    resetQ = wx.CheckBox(parent=G2frame.dataDisplay,label='Reset?')
    sqBox.Add(resetQ,0)
    resetQ.Bind(wx.EVT_CHECKBOX, OnResetQ)
    
    mainSizer.Add(sqBox,0)

    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    G2frame.dataDisplay.SetSize(Size)
    G2frame.dataFrame.setSizePosLeft(Size)
    
