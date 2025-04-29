# -*- coding: utf-8 -*-
#GSASIIpwdGUI - powder data display routines
'''GUI routines for PWDR datadree subitems follow.
'''
from __future__ import division, print_function
import platform
import sys
import os.path
# Don't depend on graphics for scriptable
try:
    import wx
    import wx.grid as wg
except ImportError:
    pass
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import math
import copy
import random as ran
import pickle
import scipy.interpolate as si
from . import GSASIIpath
from . import GSASIImath as G2mth
from . import GSASIIpwd as G2pwd
from . import GSASIIfiles as G2fil
from . import GSASIIobj as G2obj
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIindex as G2indx
from . import GSASIIplot as G2plt
from . import GSASIIpwdplot as G2pwpl
from . import GSASIIdataGUI as G2gd
from . import GSASIIphsGUI as G2phsG
from . import GSASIIctrlGUI as G2G
from . import GSASIIElemGUI as G2elemGUI
from . import GSASIIElem as G2elem
from . import GSASIIsasd as G2sasd
from . import G2shapes
from . import SUBGROUPS as kSUB
from . import k_vector_search as kvs
try:
    VERY_LIGHT_GREY = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
    WACV = wx.ALIGN_CENTER_VERTICAL
except:
    pass
if '2' not in platform.python_version_tuple()[0]:
    unichr = chr
GkDelta = unichr(0x0394)
GkSigma = unichr(0x03a3)
GkTheta = unichr(0x03f4)
Gklambda = unichr(0x03bb)
Pwr10 = unichr(0x0b9)+unichr(0x2070)
Pwr20 = unichr(0x0b2)+unichr(0x2070)
Pwrm1 = unichr(0x207b)+unichr(0x0b9)
Pwrm2 = unichr(0x207b)+unichr(0x0b2)
Pwrm6 = unichr(0x207b)+unichr(0x2076)
Pwrm4 = unichr(0x207b)+unichr(0x2074)
Angstr = unichr(0x00c5)
superMinusOne = unichr(0xaf)+unichr(0xb9)
notEq0 = unichr(0x2260)+'0'
# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
npsind = lambda x: np.sin(x*np.pi/180.)
npasind = lambda x: 180.*np.arcsin(x)/math.pi
npcosd = lambda x: np.cos(x*math.pi/180.)

cellDisplayOpts = {'showExtinct':False}
cellGUIlist = [
    [[0,1,2],4,([" a = "," Vol = "],[(10,5),"%.3f"],[True,False],[0,0])],
    [[3,4,5,6],6,([" a = "," c = "," Vol = "],[(10,5),(10,5),"%.3f"],[True,True,False],[0,2,0])],
    [[7,8,9,10,11,12],8,([" a = "," b = "," c = "," Vol = "],[(10,5),(10,5),(10,5),"%.3f"],
        [True,True,True,False],[0,1,2,0])],
    [[13,14,15,16],10,([" a = "," b = "," c = ",u'\u03B2 = '," Vol = "],
        [(10,5),(10,5),(10,5),(10,3),"%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
    [[17,18],8,([" a = "," b = "," c = ",u'\u03B1 = ',u'\u03B2 = ',u'\u03B3 = '," Vol = "],
        [(10,5),(10,5),(10,5),(10,3),(10,3),(10,3),"%.3f"],
        [True,True,True,True,True,True,False],[0,1,2,3,4,5,0])]]
bravaisSymb = ['Fm3m','Im3m','Pm3m','R3-H','P6/mmm','I4/mmm','P4/mmm',
        'Fmmm','Immm','Ammm','Bmmm','Cmmm','Pmmm','I2/m','A2/m','C2/m',
        'P2/m','P1','C1']

def SetLattice(controls):
    '''impose constraints on lattice constaints and determine the
    Bravias lattice index (ibrav) as used in cellGUIlist
    '''
    ibrav = bravaisSymb.index(controls[5])
    if controls[5] in ['Fm3m','Im3m','Pm3m']:
        controls[7] = controls[8] = controls[6]
        controls[9] = controls[10] = controls[11] = 90.
    elif controls[5] in ['R3m','P6/mmm','I4/mmm','P4/mmm']:
        controls[7] = controls[6]
        controls[9] = controls[10] = controls[11] = 90.
        if controls[5] in ['R3-H','P6/mmm']:
            controls[11] = 120.
    elif controls[5] in ['Fmmm','Immm','Ammm','Bmmm','Cmmm','Pmmm']:
        controls[9] = controls[10] = controls[11] = 90.
    elif controls[5] in ['A2/m','C2/m','P2/m','I2/m']:
        controls[9] = controls[11] = 90.  # b unique
    controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
    return ibrav

################################################################################
###### class definitions
################################################################################

class SubCellsDialog(wx.Dialog):
    'Display magnetic subcell space group information from selection in Unit Cells table of results from k-SUBGROUPSMAG'
    def __init__(self,parent,title,controls,SGData,items,phaseDict,ifPick=False):
        wx.Dialog.__init__(self,parent,-1,title,
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = None
        self.controls = controls
        self.SGData = SGData         #for parent phase
        self.items = items
        self.phaseDict = phaseDict
        self.Pick = 0
        self.ifPick = ifPick

        self.Draw()

    def Draw(self):

        def RefreshGrid(event):
            r,c =  event.GetRow(),event.GetCol()
            br = self.items[r]
            phase = self.phaseDict[br]
            rLab = magDisplay.GetRowLabelValue(r)
            pname = '(%s) %s'%(rLab,phase['Name'])
            if c == 0:
                mSGData = phase['SGData']
                text,table = G2spc.SGPrint(mSGData,AddInv=True)
                if 'magAtms' in phase:
                    msg = 'Magnetic space group information'
                    text[0] = ' Magnetic Space Group: '+mSGData['MagSpGrp']
                    text[3] = ' The magnetic lattice point group is '+mSGData['MagPtGp']
                    OprNames,SpnFlp = G2spc.GenMagOps(mSGData)
                    G2G.SGMagSpinBox(self.panel,msg,text,table,mSGData['SGCen'],OprNames,
                        mSGData['SpnFlp'],False).Show()
                else:
                    msg = 'Space Group Information'
                    G2G.SGMessageBox(self.panel,msg,text,table).Show()
            elif c == 1:
                self.Pick = r
                self.Draw()
            elif c == 2:
                maxequiv = phase['maxequiv']
                mSGData = phase['SGData']
                Uvec = phase['Uvec']
                Trans = phase['Trans']
                ifMag = False
                if 'magAtms' in phase:
                    ifMag = True
                    allmom = phase.get('allmom',False)
                    magAtms = phase.get('magAtms','')
                    mAtoms = TestMagAtoms(phase,magAtms,self.SGData,Uvec,Trans,allmom,maxequiv)
                else:
                    mAtoms = TestAtoms(phase,self.controls[15],self.SGData,Uvec,Trans,maxequiv)
                Atms = []
                AtCods = []
                atMxyz = []
                for ia,atom in enumerate(mAtoms):
                    atom[0] += '_%d'%ia
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[2:5],mSGData)
                    Atms.append(atom[:2]+['',]+atom[2:5])
                    AtCods.append('1')
                    if 'magAtms' in phase:
                        MagSytSym = G2spc.MagSytSym(SytSym,dupDir,mSGData)
                        CSI = G2spc.GetCSpqinel(mSGData['SpnFlp'],dupDir)
                        atMxyz.append([MagSytSym,CSI[0]])
                    else:
                        CSI = G2spc.GetCSxinel(SytSym)
                        atMxyz.append([SytSym,CSI[0]])
                G2phsG.UseMagAtomDialog(self.panel,pname,Atms,AtCods,atMxyz,ifMag=ifMag,ifOK=True).ShowModal()
            elif c in [3,4]:
                if c == 3:
                    title = 'Conjugacy list for '+pname
                    items = phase['altList']

                elif c == 4:
                    title = 'Super groups list list for '+pname
                    items = phase['supList']
                    if not items[0]:
                        wx.MessageBox(pname+' is a maximal subgroup',caption='Super group is parent',style=wx.ICON_INFORMATION)
                        return
                SubCellsDialog(self.panel,title,self.controls,self.SGData,items,self.phaseDict).Show()

        if self.panel: self.panel.Destroy()
        self.panel = wx.Panel(self)
        rowLabels = [str(i+1) for i in range(len(self.items))]
        if self.ifPick:
            colLabels = ['Space Gp','Pick','Uniq','nConj','nSup','Trans','Vec','a','b','c','\u03B1','\u03B2','\u03B3','Volume']
            Types = [wg.GRID_VALUE_STRING,]+[wg.GRID_VALUE_BOOL,]+3*[wg.GRID_VALUE_LONG,]+2*[wg.GRID_VALUE_STRING,]+ \
                3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+[wg.GRID_VALUE_FLOAT+':10,2']
        else:
            colLabels = ['Space Gp','Uniq','nConj','nSup','Trans','Vec','a','b','c','\u03B1','\u03B2','\u03B3','Volume']
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_LONG,]+2*[wg.GRID_VALUE_STRING,]+ \
                3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+[wg.GRID_VALUE_FLOAT+':10,2']
        table = []
        for kp,ip in enumerate(self.items):
            phase = self.phaseDict[ip]
            natms = phase.get('nAtoms',1)
            pick = False
            if kp == self.Pick:
                pick = True
            try:
                nConj = len(phase['altList'])
                nSup = len(phase['supList'])
            except KeyError:
                nConj = 0
                nSup = 0
            cell  = list(phase['Cell'])
            trans = G2spc.Trans2Text(phase['Trans'])
            vec = G2spc.Latt2text([phase['Uvec'],])
            if self.ifPick:
                row = [phase['Name'],pick,natms,nConj,nSup,trans,vec]+cell
            else:
                row = [phase['Name'],natms,nConj,nSup,trans,vec]+cell
            table.append(row)
        CellsTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        magDisplay = G2G.GSGrid(self.panel)
        magDisplay.SetTable(CellsTable, True)
        magDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,RefreshGrid)
        magDisplay.AutoSizeColumns(False)
        mainSizer.Add(magDisplay,0)

        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def GetSelection(self):
        return self.Pick

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

class RDFDialog(wx.Dialog):
    'Display controls for generating RDF plot in Background'
    def __init__(self,parent):
        wx.Dialog.__init__(self,parent,-1,'Background radial distribution function',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = None
        self.result = {'UseObsCalc':'obs-calc','maxR':20.0,'Smooth':'linear'}

        self.Draw()

    def Draw(self):

        def OnUseOC(event):
            self.result['UseObsCalc'] = useOC.GetValue()

        def OnSmCombo(event):
            self.result['Smooth'] = smCombo.GetValue()

        if self.panel: self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label='Background RDF controls:'),0)
        plotType = wx.BoxSizer(wx.HORIZONTAL)
        plotType.Add(wx.StaticText(self.panel,label=' Select plot type:'),0,WACV)
        Choices = ['obs-back','calc-back','obs-calc','auto-back']
        useOC = wx.ComboBox(self.panel,value=Choices[2],choices=Choices,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        useOC.SetValue(self.result['UseObsCalc'])
        useOC.Bind(wx.EVT_COMBOBOX,OnUseOC)
        plotType.Add(useOC,0,WACV)
        mainSizer.Add(plotType,0)
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        dataSizer.Add(wx.StaticText(self.panel,label=' Smoothing type: '),0,WACV)
        smChoice = ['linear','nearest',]
        smCombo = wx.ComboBox(self.panel,value=self.result['Smooth'],choices=smChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        smCombo.Bind(wx.EVT_COMBOBOX, OnSmCombo)
        dataSizer.Add(smCombo,0,WACV)
        dataSizer.Add(wx.StaticText(self.panel,label=' Maximum radial dist.: '),0,WACV)
        maxR = G2G.ValidatedTxtCtrl(self.panel,self.result,'maxR',nDig=(10,1),xmin=10.,xmax=50.,
            typeHint=float)
        dataSizer.Add(maxR,0,WACV)
        mainSizer.Add(dataSizer,0)

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

def GetFileBackground(G2frame,xye,background,scale=True):
    ''' Select a background file to subtract from PWDR pattern
    param: xye list [npts,6] of PWDR pattern
    param: background PWDR file to be used as background
    param: scale bool:=True if scale mult included in background & apply it
    returns: list background to subtract
    '''
    bxye = np.zeros(len(xye[1]))
    mult = 1.0
    if 'background PWDR' in background[1]:
        backfile,mult = background[1]['background PWDR'][:2]
        if backfile:
            bId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,backfile)
            if bId:
                bxye = G2frame.GPXtree.GetItemPyData(bId)[1][1]
            else:
                print('Error: background PWDR {} not found'.format(backfile))
                background[1]['background PWDR'] = ['',1.0,False]
    if scale:
        return bxye*mult
    else:
        return bxye

def IsHistogramInAnyPhase(G2frame,histoName):
    '''Tests a Histogram to see if it is linked to any phases.
    Returns the name of the first phase where the histogram is used.
    '''
    phases = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
    if phases:
        item, cookie = G2frame.GPXtree.GetFirstChild(phases)
        while item:
            data = G2frame.GPXtree.GetItemPyData(item)
            histoList = data['Histograms'].keys()
            if histoName in histoList:
                return G2frame.GPXtree.GetItemText(item)
            item, cookie = G2frame.GPXtree.GetNextChild(phases, cookie)
        return False
    else:
        return False

def GetPhasesforHistogram(G2frame,histoName):
    '''Returns phases (if any) associated with provided Histogram
    Returns a list of phase dicts
    '''
    phList = []
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    for ph in Phases:
        if histoName in Phases[ph]['Histograms']:
            phList.append(Phases[ph])
    return phList

def SetupSampleLabels(histName,dataType,histType):
    '''Setup a list of labels and number formatting for use in
    labeling sample parameters.
    :param str histName: Name of histogram, ("PWDR ...")
    :param str dataType:
    '''
    parms = []
    parms.append(['Scale','Histogram scale factor: ',[10,7]])
    if histType[2] in ['A','B','C']:
        parms.append(['Gonio. radius','Goniometer radius (mm): ',[10,3]])
    if 'PWDR' in histName:
        if dataType == 'Debye-Scherrer':
            if 'T' in histType:
                parms += [['Absorption',u'Sample absorption (\xb5r/'+Gklambda+'): ',[10,4]],]
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

def GetFileList(G2frame,fileType):
    ''' Get list of file names containing a particular string
    param: fileType str: any string within a file name
    returns: list of file names from GSAS-II tree
    Note routine of same name in GSASIIdataGUI; it has a skip option
    '''
    fileList = []
    Id, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
    while Id:
        name = G2frame.GPXtree.GetItemText(Id)
        if fileType in name.split()[0]:
            fileList.append(name)
        Id, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
    return fileList

def GetHistsLikeSelected(G2frame):
    '''Get the histograms that match the current selected one:
    The histogram prefix and data type (PXC etc.), the number of
    wavelengths and the instrument geometry (Debye-Scherrer etc.)
    must all match. The current histogram is not included in the list.

    :param wx.Frame G2frame: pointer to main GSAS-II data tree
    '''
    histList = []
    inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
    hType = inst['Type'][0]
    if 'Lam1' in inst:
        hLam = 2
    elif 'Lam' in inst:
        hLam = 1
    else:
        hLam = 0
#    sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Sample Parameters'))
#    hGeom = sample.get('Type')
    hstName = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    hPrefix = hstName.split()[0]+' '
    # cycle through tree looking for items that match the above
    item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
    while item:
        name = G2frame.GPXtree.GetItemText(item)
        if name.startswith(hPrefix) and name != hstName:
            cType,cLam, = '?',-1
            subitem, subcookie = G2frame.GPXtree.GetFirstChild(item)
            while subitem:
                subname = G2frame.GPXtree.GetItemText(subitem)
                if subname == 'Sample Parameters':
                    # sample = G2frame.GPXtree.GetItemPyData(subitem)
                    # cGeom = sample.get('Type')
                    pass
                elif subname == 'Instrument Parameters':
                    inst,inst2 = G2frame.GPXtree.GetItemPyData(subitem)
                    cType = inst['Type'][0]
                    if 'Lam1' in inst:
                        cLam = 2
                    elif 'Lam' in inst:
                        cLam = 1
                    else:
                        cLam = 0
                subitem, subcookie = G2frame.GPXtree.GetNextChild(item, subcookie)
            if cLam == hLam and cType == hType: # and cGeom == hGeom:
                if name not in histList: histList.append(name)
        item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
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
    hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    histList = GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No other histograms match '+hst,G2frame)
        return
    sourceData = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)

    if 'Offset' not in sourceData[0]:    #patch for old data
        sourceData[0].update({'Offset':[0.0,0.0],'delOffset':0.02,'refOffset':-1.0,
            'refDelt':0.01,})
        G2frame.GPXtree.SetItemPyData(G2frame.PatternId,sourceData)

    dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy plot controls from\n'+str(hst[5:])+' to...',
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

    keys = ['Offset','delOffset','refOffset','refDelt']
    source = dict(zip(keys,[sourceData[0][item] for item in keys]))
    for hist in copyList:
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
        data = G2frame.GPXtree.GetItemPyData(Id)
        data[0].update(source)
        G2frame.GPXtree.SetItemPyData(Id,data)
    print ('Copy of plot controls successful')

def CopySelectedHistItems(G2frame):
    '''Global copy: Copy items from current histogram to others.
    '''
    hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    histList = GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No other histograms match '+hst,G2frame)
        return
    choices = ['Limits','Background','Instrument Parameters','Sample Parameters']
    dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy which histogram sections from\n'+str(hst[5:]),
        'Select copy sections', choices, filterBox=False)
    dlg.SetSelections(range(len(choices)))
    choiceList = []
    if dlg.ShowModal() == wx.ID_OK:
        choiceList = [choices[i] for i in dlg.GetSelections()]
    if not choiceList: return

    dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy parameters from\n'+str(hst[5:])+' to...',
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
        data = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Limits'))
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.SetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Limits'),
                copy.deepcopy(data))
    if 'Background' in choiceList:  # Background
        data = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Background'))
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.SetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Background'),
                copy.deepcopy(data))
    if 'Instrument Parameters' in choiceList:  # Instrument Parameters
        # for now all items in Inst. parms are copied
        data,data1 = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(
                G2frame,G2frame.PatternId,'Instrument Parameters'))
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters')
                )[0].update(copy.deepcopy(data))
            G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters')
                )[1].update(copy.deepcopy(data1))
    if 'Sample Parameters' in choiceList:  # Sample Parameters
        data = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(
                G2frame,G2frame.PatternId,'Sample Parameters'))
        # selects items to be copied
        histType,copyNames = SetCopyNames(hst,data['Type'],
            addNames = ['Omega','Chi','Phi','Gonio. radius','InstrName'])
        copyDict = {parm:data[parm] for parm in copyNames}
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters')
                ).update(copy.deepcopy(copyDict))

def TestMagAtoms(phase,magAtms,SGData,Uvec,Trans,allmom,maxequiv=100,maximal=False):
    ''' Tests substructure magnetic atoms for magnetic site symmetry
    param: phase GSAS-II phase object
    param: magAtms list:magnetic atom objects
    param: SGData dict: GSAS-II space group object
    param: Uvec array: Translation U vector
    param: Trans array: Transformation matrix
    param: allmom bool: =True if all atoms must have moments allowed
    params: maxequiv int:maximum number of atoms with moments to consider
    param: maximal bool:=True if maximal subgroups only are allowed
    returns: unique magnetic atoms (if any)
    '''
    found = False
    anymom = False
    phase['Keep'] = False
    if not magAtms:
        phase['Keep'] = True
        return []
    invTrans = nl.inv(Trans)
    atCodes = []
    Phase = {'General':{'AtomPtrs':[2,1],'SGData':copy.deepcopy(phase['SGData'])},'Atoms':[]}
    for matm in magAtms:
        XYZ = G2spc.GenAtom(matm[3:6],SGData,False,Move=True)
        xyzs = [xyz[0] for xyz in XYZ]
        atCodes += len(xyzs)*['1',]
        xyzs,atCodes = G2lat.ExpandCell(xyzs,atCodes,0,Trans)
        for ix,x in enumerate(xyzs):
            xyz = G2lat.TransformXYZ(x-Uvec,invTrans.T,np.zeros(3))%1.
            Phase['Atoms'].append(matm[:2]+list(xyz))
            SytSym,Mul,Nop,dupDir = G2spc.SytSym(xyz,phase['SGData'])
            CSI = G2spc.GetCSpqinel(phase['SGData']['SpnFlp'],dupDir)
            if any(CSI[0]):
                anymom = True
            if allmom:
                if not any(CSI[0]):
                    phase['Keep'] = False
                    found = True
    uAtms = G2lat.GetUnique(Phase,atCodes)[0]
    natm = len(uAtms)
    if anymom and natm <= maxequiv and not found:
        phase['Keep'] = True
        if maximal and phase['supList'][0]:
            phase['Keep'] = False
    return uAtms

def TestAtoms(phase,magAtms,SGData,Uvec,Trans,maxequiv=100,maximal=False):
    '''Tests atoms for substructure equivalents
    param: phase GSAS-II phase object
    param: magAtms list: atom objects
    param: SGData dict: GSAS-II space group object
    param: Uvec array: Translation U vector
    param: Trans array: Transformation matrix
    params: maxequiv int:maximum number of atoms with moments to consider
    param: maximal bool:=True if maximal subgroups only are allowed
    returns: unique atoms (if any)
    '''
    phase['Keep'] = True
    invTrans = nl.inv(Trans)
    atCodes = []
    Phase = {'General':{'AtomPtrs':[2,1],'SGData':copy.deepcopy(phase['SGData'])},'Atoms':[]}
    for matm in magAtms:
        XYZ = G2spc.GenAtom(matm[3:6],SGData,False,Move=True)
        xyzs = [xyz[0] for xyz in XYZ]
        atCodes += len(xyzs)*['1',]
        xyzs,atCodes = G2lat.ExpandCell(xyzs,atCodes,0,Trans)
        for ix,x in enumerate(xyzs):
            xyz = G2lat.TransformXYZ(x-Uvec,invTrans.T,np.zeros(3))%1.
            Phase['Atoms'].append(matm[:2]+list(xyz))
    uAtms = G2lat.GetUnique(Phase,atCodes)[0]
    natm = len(uAtms)
    if natm > maxequiv: #too many allowed atoms found
        phase['Keep'] = False
    if maximal and phase['supList'][0]:
        phase['Keep'] = False
    return uAtms

def RefineCell(G2frame):
    '''Refine a unit cell, called by OnRefine'''
    def cellPrint(ibrav,A):
        cell = G2lat.A2cell(A)
        Vol = G2lat.calc_V(A)
        if ibrav in ['Fm3m','Im3m','Pm3m']:
            print (" %s%10.6f" % ('a =',cell[0]))
        elif ibrav in ['R3-H','P6/mmm','I4/mmm','P4/mmm']:
            print (" %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],' c =',cell[2],' volume =',Vol))
        elif ibrav in ['P4/mmm','Fmmm','Immm','Ammm','Bmmm','Cmmm','Pmmm']:
            print (" %s%10.6f %s%10.6f %s%10.6f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],' volume =',Vol))
        elif ibrav in ['I2/m','A2/m','C2/m','P2/m']:
            print (" %s%10.6f %s%10.6f %s%10.6f %s%8.3f %s%12.3f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2],'beta =',cell[4],' volume =',Vol))
        else:
            print (" %s%10.6f %s%10.6f %s%10.6f" % ('a =',cell[0],'b =',cell[1],'c =',cell[2]))
            print (" %s%8.3f %s%8.3f %s%8.3f %s%12.3f" % ('alpha =',cell[3],'beta =',cell[4],'gamma =',cell[5],' volume =',Vol))

    def vecPrint(Vec):
        print (' %s %10.5f %10.5f %10.5f'%('Modulation vector:',Vec[0],Vec[1],Vec[2]))

    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))[1]
    if 'T' in Inst['Type'][0]:
        difC = Inst['difC'][1]
        dmin = G2lat.Pos2dsp(Inst,Limits[0])
    elif 'E' in Inst['Type'][0]:
        TTh = Inst['2-theta'][1]
        dmin = G2lat.Pos2dsp(Inst,Limits[1])
    else:   #'C', 'B', or 'PKS'
        wave = G2mth.getWave(Inst)
        dmin = G2lat.Pos2dsp(Inst,Limits[1])
    PatternId = G2frame.PatternId
    peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Index Peak List'))
    controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List'))
    cell = controls[6:12]
    A = G2lat.cell2A(cell)
    ibrav = bravaisSymb.index(controls[5])
    SGData = G2spc.SpcGroup(controls[13])[1]
    if 'T' in Inst['Type'][0]:
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
    elif 'E' in Inst['Type'][0]:        #no super lattice stuff for EDX data - resolution too low
        G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
        peaks = [G2indx.IndexPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #put peak fit esds back in peaks
        Lhkl,M20,X20,Aref = G2indx.refinePeaksE(peaks[0],TTh,ibrav,A)
        Zero = 0.0
    else:
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
    controls[1] = Zero
    controls[6:12] = G2lat.A2cell(Aref)
    controls[12] = G2lat.calc_V(Aref)
    if ssopt.get('Use',False):
        ssopt['ModVec'] = Vec
        G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,ssopt['ModVec'],ssopt['maxH'],A)
    else:
        G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
        newcell = [M20,X20,ibrav]+controls[6:13]+[False,False]
        cells.append(newcell)
        cells = G2indx.sortM20(cells)
    G2frame.HKL = np.array(G2frame.HKL)
    data = [controls,bravais,cells,dmin,ssopt,magcells]
    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
    print (" %s%10.3f" % ('refinement M20 = ',M20))
    print (' unindexed lines = %d'%X20)
    cellPrint(controls[5],Aref)
    ip = 4
    if ssopt.get('Use',False):
        vecPrint(Vec)
        ip = 5
    for hkl in G2frame.HKL:
        hkl[ip] = G2lat.Dsp2pos(Inst,hkl[ip-1])+controls[1]
    if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
        G2plt.PlotPowderLines(G2frame,indexFrom='Indexing from refine cell, M20=%.3f'%M20)
    else:
        G2pwpl.PlotPatterns(G2frame,indexFrom='Indexing from refine cell, M20=%.3f'%M20)
    return data

################################################################################
#####  Powder Peaks
################################################################################

def UpdatePeakGrid(G2frame, data):
    '''respond to selection of PWDR powder peaks data tree item.
    '''
    def OnAutoSearch(event):
        'Search pattern for possible peak positions'
        PatternId = G2frame.PatternId
        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
        background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
        inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        Pattern = G2frame.GPXtree.GetItemPyData(PatternId)
        profile = Pattern[1]
        bxye = GetFileBackground(G2frame,profile,background)
        x0 = profile[0]
        iBeg = np.searchsorted(x0,limits[0])
        iFin = np.searchsorted(x0,limits[1])
        x = x0[iBeg:iFin]
        y0 = (profile[1]-bxye)[iBeg:iFin]
        ysig = 1.0*np.std(y0)
        offset = [-1,1]
        ymask = ma.array(y0,mask=(y0<ysig))
        for off in offset:
            ymask = ma.array(ymask,mask=(ymask-np.roll(y0,off)<=0.))
        indx = ymask.nonzero()
        mags = ymask[indx]
        poss = x[indx]
        refs = list(zip(poss,mags))
        if 'T' in Inst['Type'][0]:
            refs = G2mth.sortArray(refs,0,reverse=True)     #big TOFs first
        else:   #'C', 'E' or 'B'
            refs = G2mth.sortArray(refs,0,reverse=False)    #small 2-Thetas or energies first
        for i,ref1 in enumerate(refs):      #reject picks closer than 1 FWHM
            for ref2 in refs[i+1:]:
                if abs(ref2[0]-ref1[0]) < 2.*G2pwd.getFWHM(ref1[0],inst):
                    del(refs[i])
        if 'T' in Inst['Type'][0]:
            refs = G2mth.sortArray(refs,1,reverse=False)
        else:   #'C', 'E' or 'B'
            refs = G2mth.sortArray(refs,1,reverse=True)
        for pos,mag in refs:
            data['peaks'].append(G2mth.setPeakparms(inst,inst2,pos,mag))
        UpdatePeakGrid(G2frame,data)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnCopyPeaks(event):
        'Copy peaks to other histograms'
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy peak list from\n'+str(hst[5:])+' to...',
            'Copy peaks', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.SetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Peak List'),copy.deepcopy(data))

    def OnLoadPeaks(event):
        'Load peak list from file'
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II PWDR peaks list file', pth, '',
            'PWDR peak list files (*.pkslst)|*.pkslst',wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = []
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                while S:
                    if '#' in S:
                        S = File.readline()
                        continue
                    try:
                        peaks.append(eval(S))
                    except:
                        break
                    S = File.readline()
                File.close()
        finally:
            dlg.Destroy()
        data = {'peaks':peaks,'sigDict':{}}
        UpdatePeakGrid(G2frame,data)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnSavePeaks(event):
        'Save peak to file suitable for OnLoadPeaks'
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II PWDR peaks list file', pth, '',
            'PWDR peak list files (*.pkslst)|*.pkslst',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .pkslst
                filename = os.path.splitext(filename)[0]+'.pkslst'
                File = open(filename,'w')
                File.write("#GSAS-II PWDR peaks list file; do not add/delete items!\n")
                for item in data:
                    if item == 'peaks':
                        for pk in data[item]:
                            File.write(str(pk)+'\n')
                File.close()
                print ('PWDR peaks list saved to: '+filename)
        finally:
            dlg.Destroy()

    def OnUnDo(event):
        'Undo a peak fit - reads a saved file from PeakFit'
        file = open(G2frame.undofile,'rb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            Id = G2gd.GetGPXtreeItemId(G2frame,PatternId, item)
            oldvals = pickle.load(file)
            G2frame.GPXtree.SetItemPyData(Id,oldvals)
            if item == 'Peak List':
                data.update(G2frame.GPXtree.GetItemPyData(Id))
            print (item+' recovered')
        file.close()
        G2frame.dataWindow.UnDo.Enable(False)
        wx.CallAfter(UpdatePeakGrid,G2frame,data)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def SaveState():
        'Saves result of a peaak fit for possible UnDo'
        G2frame.undofile = os.path.join(G2frame.dirname,'GSASII.save')
        file = open(G2frame.undofile,'wb')
        PatternId = G2frame.PatternId
        for item in ['Background','Instrument Parameters','Peak List']:
            pickle.dump(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,item)),file,1)
        file.close()
        G2frame.dataWindow.UnDo.Enable(True)

    def OnLSQPeakFit(event):
        'Do a peak fit'
        if reflGrid.IsCellEditControlEnabled(): # complete any grid edits in progress
            reflGrid.HideCellEditControl()
            reflGrid.DisableCellEditControl()
        if not G2frame.GSASprojectfile:            #force a save of the gpx file so SaveState can write in the same directory
            G2frame.OnFileSaveas(event)
        wx.CallAfter(OnPeakFit)

    def OnOneCycle(event):
        'Do a single cycle of peak fit'
        if reflGrid.IsCellEditControlEnabled(): # complete any grid edits in progress
            reflGrid.HideCellEditControl()
            reflGrid.DisableCellEditControl()
        wx.CallAfter(OnPeakFit,oneCycle=True)

    def OnSeqPeakFit(event):
        '''Do a sequential peak fit across multiple histograms - peaks must be present in all.
        results saved in Sequential peak fit results'''
        histList = G2gd.GetGPXtreeDataNames(G2frame,['PWDR',])
        od = {'label_1':'Copy to next','value_1':False,'label_2':'Reverse order','value_2':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Sequential peak fits',
             'Select dataset to include',histList,extraOpts=od)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(histList[sel])
        dlg.Destroy()
        if not names:
            return
        Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential peak fit results')
        if Id:
            SeqResult = G2frame.GPXtree.GetItemPyData(Id)
        else:
            SeqResult = {}
            Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential peak fit results')
        SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
        SeqResult['histNames'] = names
        dlg = wx.ProgressDialog('Sequential peak fit','Data set name = '+names[0],len(names),
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        controls = {'deriv type':'analytic','min dM/M':0.001,}
        print ('Peak Fitting with '+controls['deriv type']+' derivatives:')
        oneCycle = False
        prevVaryList = []
        peaks = None
        varyList = None
        if od['value_2']:
            names.reverse()
        try:
            for i,name in enumerate(names):
                print (' Sequential fit for '+name)
                dlg.Raise()
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    dlg.Destroy()
                    break
                PatternId =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                if i and od['value_1']:
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'),copy.deepcopy(peaks))
                    prevVaryList = varyList[:]
                peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'))
                background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
                limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
                inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
                Pattern = G2frame.GPXtree.GetItemPyData(PatternId)
                data = Pattern[1]
                fixback = GetFileBackground(G2frame,data,background,scale=False)
                peaks['sigDict'],result,sig,Rvals,varyList,parmDict,fullvaryList,badVary = G2pwd.DoPeakFit(None,peaks['peaks'],
                    background,limits,inst,inst2,data,fixback,prevVaryList,oneCycle,controls)   #needs wtFactor after controls?
                if len(result[0]) != len(fullvaryList):
                    dlg.Destroy()
                    print (' ***** Sequential peak fit stopped at '+name+' *****')
                    break
                else:
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'),copy.deepcopy(peaks))
                    SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                        'covMatrix':np.eye(len(result[0])),'title':name,'parmDict':parmDict,
                        'fullVary':fullvaryList,'badVary':badVary}
            print (' ***** Sequential peak fit successful *****')
        finally:
            dlg.Destroy()
        SeqResult['histNames'] = histList
        G2frame.GPXtree.SetItemPyData(Id,SeqResult)
        G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
        G2frame.GPXtree.SelectItem(Id)

    def OnDelPeaks(event):
        'Delete selected peaks from the Peak fit table'
        if G2frame.dataWindow.XtraPeakMode.IsChecked(): # which table?
            tbl = data['xtraPeaks']
        else:
            tbl = data['peaks']
        choices = [f"{i[0]:.2f}" for i in tbl]
        if not choices: return
        sel = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select peaks to delete',
                'Delete peaks',choices)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
        finally:
            dlg.Destroy()
        for i in sorted(sel,reverse=True):
            del tbl[i]
        UpdatePeakGrid(G2frame,data)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnClearPeaks(event):
        'Clear the Peak fit table'
        dlg = wx.MessageDialog(G2frame,'Delete all peaks?','Clear peak list',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                peaks = {'peaks':[],'sigDict':{}}
        finally:
            dlg.Destroy()
        UpdatePeakGrid(G2frame,peaks)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnPeakFit(oneCycle=False,noFit=False):
        'Do peak fitting by least squares'
        SaveState()
        controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        if not controls:
            controls = {'deriv type':'analytic','min dM/M':0.001,}     #fill in defaults if needed
        #print ('Peak Fitting with '+controls['deriv type']+' derivatives:')
        PatternId = G2frame.PatternId
        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        Pattern = G2frame.GPXtree.GetItemPyData(PatternId)
        data = Pattern[1]
        wtFactor = Pattern[0]['wtFactor']
        overallInfo = {}
        peaks['LaueFringe'] = peaks.get('LaueFringe',{})
        peaks['LaueFringe']['satellites'] = []
        lines = peaks['LaueFringe'].get('Show')
        if 'LF' in inst['Type'][0]:
            wave = G2mth.getWave(inst)
            overallInfo = {'ncell':peaks['LaueFringe']['ncell'],
                            'clat': peaks['LaueFringe']['clat'],
                            'clat-ref': peaks['LaueFringe']['clat-ref'],
                            'fitRange': peaks['LaueFringe'].get('fitRange',8.0),
                            'fitPowerM': peaks['LaueFringe'].get('fitPowerM',2.0),
                            'fitPowerP': peaks['LaueFringe'].get('fitPowerP',2.0),
                               } # add overall info here
            if lines:
                for i in peaks['peaks']:
                    pks = list(range(-lines,0)) + list(range(1,lines+1))
                    peaks['LaueFringe']['satellites'].extend(
                        G2pwd.LaueSatellite(i[0],wave,peaks['LaueFringe']['clat'],peaks['LaueFringe']['ncell'],pks))

        if G2frame.dataWindow.XtraPeakMode.IsChecked(): # adding peaks to computed pattern
            histoName = G2frame.GPXtree.GetItemText(G2frame.PatternId)
# do zero cycle refinement
            from . import GSASIIstrMain as G2stMn
            # recompute current pattern for current histogram, set as fixed background
            bxye = G2stMn.DoNoFit(G2frame.GSASprojectfile,histoName)
            peaksplus = peaks['xtraPeaks'] + [{}]
            # dummy out background parameters (are not used so can't be refined)
            background = [['chebyschev-1', False, 1, 0.0], {
                'nDebye': 0,'debyeTerms': [],'nPeaks': 0,'peaksList': [],
                'background PWDR': ['', 1.0, False],
                'FixedPoints': [],'autoPrms': {}}]
            #breakpoint() ##################################################
        else:
            peaksplus = peaks['peaks'] + [overallInfo]
            bxye = GetFileBackground(G2frame,data,background,scale=False)
        if noFit:
            results = G2pwd.DoPeakFit(None,peaksplus,background,limits,inst,inst2,data,bxye,[],oneCycle,controls,wtFactor,noFit=True)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
            return
        # try:
        dlg = wx.ProgressDialog('Residual','Peak fit Rwp = ',101,parent=G2frame,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        results = G2pwd.DoPeakFit(None,peaksplus,background,limits,inst,inst2,data,bxye,[],oneCycle,controls,wtFactor,dlg)
        # finally:
        #     dlg.Destroy()
        if results is None:
            return
        peaks['sigDict'] = results[0]
        if 'LF' in inst['Type'][0] and 'clat' in overallInfo:
            peaks['LaueFringe']['clat'] = overallInfo['clat']
        text = 'Peak fit: Rwp=%.2f%% Nobs= %d Nparm= %d Npeaks= %d'%(results[3]['Rwp'],results[1][2]['fjac'].shape[1],len(results[0]),len(peaks['peaks']))
        newpeaks = copy.copy(peaks)
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'),newpeaks)
        G2frame.AddToNotebook(text,'PF')
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
        wx.CallAfter(UpdatePeakGrid,G2frame,newpeaks)

    def OnResetSigGam(event):
        'Reset sig & gam values to instrument parameter values'
        PatternId = G2frame.PatternId
        Inst,Inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'))
        if not peaks['peaks']:
            G2frame.ErrorDialog('No peaks!','Nothing to do!')
            return
        newpeaks = {'peaks':[],'sigDict':{}}
        for peak in peaks['peaks']:
            newpeaks['peaks'].append(G2mth.setPeakparms(Inst,Inst2,peak[0],peak[2]))
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'),newpeaks)
        UpdatePeakGrid(G2frame,newpeaks)

    def setBackgroundColors():
        'Set background colors in peak list table; red if negative (nonsense), white if ok'
        for r in range(reflGrid.GetNumberRows()):
            for c in range(reflGrid.GetNumberCols()):
                if reflGrid.GetColLabelValue(c) in ['position','intensity','alpha','beta','sigma\u00b2','gamma']:
                    try:
                        if float(reflGrid.GetCellValue(r,c)) < 0.:
                            reflGrid.SetCellBackgroundColour(r,c,wx.RED)
                        else:
                            reflGrid.SetCellBackgroundColour(r,c,wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                    except:
                        pass

    def KeyEditPeakGrid(event):
        '''Respond to pressing a key to act on selection of a row, column or cell
        in the Peak List table
        '''
        rowList = reflGrid.GetSelectedRows()
        colList = reflGrid.GetSelectedCols()
        selectList = reflGrid.GetSelectedCells()
        data = G2frame.GPXtree.GetItemPyData(G2frame.PickId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif rowList and (event.GetKeyCode() == wx.WXK_DELETE or event.GetKeyCode() == 8):
            # pressing the delete key or backspace deletes selected peak(s)
            reflGrid.ClearSelection()
            reflGrid.ClearGrid()
            rowList.sort()
            rowList.reverse()
            nDel = 0
            for row in rowList:
                G2frame.PeakTable.DeleteRow(row)
                nDel += 1
            if nDel:
                msg = wg.GridTableMessage(G2frame.PeakTable,
                    wg.GRIDTABLE_NOTIFY_ROWS_DELETED,0,nDel)
                reflGrid.ProcessTableMessage(msg)
            data['peaks'] = G2frame.PeakTable.GetData()[:-nDel]
            G2frame.GPXtree.SetItemPyData(G2frame.PickId,data)
            setBackgroundColors()
        elif colList and (event.GetKeyCode() == 89 or event.GetKeyCode() == 78):
            reflGrid.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.PeakTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                    if key == 89: #'Y'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data['peaks'][row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.PeakTable.GetNumberRows()): data['peaks'][row][col]=False
        elif selectList and (event.GetKeyCode() == 89 or event.GetKeyCode() == 78):
            reflGrid.ClearSelection()
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
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
        wx.CallAfter(UpdatePeakGrid,G2frame,data)

    def SelectVars(rows):
        '''Set or clear peak refinement variables for peaks listed in rows
        '''
        refOpts = {reflGrid.GetColLabelValue(i):i+1 for i in range(reflGrid.GetNumberCols()) if reflGrid.GetColLabelValue(i) != "refine"}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select columns to refine',
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
        tbl = reflGrid.GetTable().data
        for r in rows:
            for lbl,c in refOpts.items():
                tbl[r][c] = lbl in sels
        UpdatePeakGrid(G2frame,data)

    def OnRefineSelected(event):
        '''set refinement flags for the selected peaks
        '''
        rows = list(set([row for row,col in reflGrid.GetSelectedCells()] +
                        reflGrid.GetSelectedRows()))
        tbl = reflGrid.GetTable().data
        if not rows:
            choices = [f"{i[0]:.2f}" for i in tbl]
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Select peaks to refine',
                'select peaks',choices)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    rows = dlg.GetSelections()
            finally:
                dlg.Destroy()
        if not rows: return
        SelectVars(rows)

    def OnRefineAll(event):
        '''set refinement flags for all peaks
        '''
        SelectVars(range(reflGrid.GetNumberRows()))

    def onCellListDClick(event):
        '''Called after a double-click on a row/column label'''
        r,c =  event.GetRow(),event.GetCol()
        if r < 0 and c < 0:
            for row in range(reflGrid.GetNumberRows()):
                reflGrid.SelectRow(row,True)
            for col in range(reflGrid.GetNumberCols()):
                reflGrid.SelectCol(col,True)
        elif r >= 0 and c < 0:     #row label: select it and replot!
            reflGrid.ClearSelection()
            reflGrid.SelectRow(r,True)
            wx.CallAfter(G2frame.reflGrid.ForceRefresh)
            wx.CallAfter(G2pwpl.PlotPatterns,G2frame,plotType='PWDR')
        elif c > 0:     #column label: just select it (& redisplay)
            reflGrid.ClearSelection()
            reflGrid.SelectCol(c,True)
            if reflGrid.GetColLabelValue(c) != 'refine': return
            choice = ['Y - vary all','N - vary none',]
            dlg = wx.SingleChoiceDialog(G2frame,'Select refinement option for '+reflGrid.GetColLabelValue(c-1),
                'Refinement controls',choice)
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                tbl = reflGrid.GetTable().data
                if sel == 0:
                    for r in range(reflGrid.GetNumberRows()): tbl[r][c]=True
                else:
                    for r in range(reflGrid.GetNumberRows()): tbl[r][c]=False
            wx.CallAfter(UpdatePeakGrid,G2frame,data)

    def RefreshPeakGrid(event):
        'recompute & plot the peaks any time a value in the table is edited'
        if 'LF' in Inst['Type'][0]:
            for i in range(len(data['LFpeaks'])):
                data['peaks'][i][2:] = data['LFpeaks'][i]
            wx.CallAfter(UpdatePeakGrid,G2frame,data)
        if data['peaks']:
            OnPeakFit(noFit=True)

    def ToggleXtraMode(event):
        '''Switch "Extra Peak" mode in response to button'''
        G2frame.dataWindow.XtraPeakMode.Check(
            not G2frame.dataWindow.XtraPeakMode.IsChecked())
        OnXtraMode(event)

    def OnXtraMode(event):
        '''Respond to change in "Extra Peak" mode from menu command
        or button via :func:`ToggleXtraMode`.
        '''
        data['xtraMode'] = G2frame.dataWindow.XtraPeakMode.IsChecked()
        wx.CallAfter(UpdatePeakGrid,G2frame,data)
        wx.CallAfter(G2pwpl.PlotPatterns,G2frame,plotType='PWDR')

    def OnSetPeakWidMode(event):
        '''Toggle G2pwd.peakInstPrmMode mode; determines if unvaried
        sigma and gamma values are set from UVW & XY
        '''
        state = G2frame.dataWindow.setPeakMode.IsChecked()
        G2pwd.setPeakInstPrmMode(state)

    def ShiftLFc(event):
        '''Shift the Laue Fringe lattice parameter
        '''
        Obj = event.GetEventObject()
        move = Obj.GetValue()
        Obj.SetValue(0)
        data['LaueFringe']['clat'] *= 1. + move/2000.
        wx.CallAfter(RefreshPeakGrid,None)

    #### beginning of UpdatePeakGrid =================================
    G2frame.GetStatusBar().SetStatusText('Global refine: select refine column & press Y or N',1)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu)
    G2frame.Bind(wx.EVT_MENU, OnAutoSearch, id=G2G.wxID_AUTOSEARCH)
    G2frame.Bind(wx.EVT_MENU, OnCopyPeaks, id=G2G.wxID_PEAKSCOPY)
    G2frame.Bind(wx.EVT_MENU, OnSavePeaks, id=G2G.wxID_PEAKSAVE)
    G2frame.Bind(wx.EVT_MENU, OnLoadPeaks, id=G2G.wxID_PEAKLOAD)
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2G.wxID_UNDO)
    G2frame.Bind(wx.EVT_MENU, OnRefineSelected, id=G2frame.dataWindow.peaksSel.GetId())
    G2frame.Bind(wx.EVT_MENU, OnRefineAll, id=G2frame.dataWindow.peaksAll.GetId())
    G2frame.Bind(wx.EVT_MENU, OnLSQPeakFit, id=G2G.wxID_LSQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnOneCycle, id=G2G.wxID_LSQONECYCLE)
    G2frame.Bind(wx.EVT_MENU, OnSeqPeakFit, id=G2G.wxID_SEQPEAKFIT)
    G2frame.Bind(wx.EVT_MENU, OnDelPeaks, id=G2G.wxID_DELPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnClearPeaks, id=G2G.wxID_CLEARPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnResetSigGam, id=G2G.wxID_RESETSIGGAM)
    G2frame.Bind(wx.EVT_MENU, OnSetPeakWidMode, id=G2G.wxID_SETUNVARIEDWIDTHS)
    # get info in phases associated with current histogram
    data['xtraPeaks'] = data.get('xtraPeaks',[])
    histoName = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    phCount = 0
    # isHistMag = False
    for ph in GetPhasesforHistogram(G2frame,histoName):
         phCount += 1
    #     if ph['General']['Type'] == 'magnetic':
    #        isHistMag = True
    #        break
    if phCount == 0: # no phases for this histogram
        G2frame.dataWindow.XtraPeakMode.Enable(False)
    else:
        G2frame.dataWindow.XtraPeakMode.Enable(True)
        data['xtraMode'] = data.get('xtraMode',False)
        G2frame.dataWindow.XtraPeakMode.Check(data['xtraMode'])
        G2frame.Bind(wx.EVT_MENU, OnXtraMode, id=G2G.wxID_XTRAPEAKMODE)
    OnSetPeakWidMode(None)
    # can peaks be refined?
    pkmode = False
    if G2frame.dataWindow.XtraPeakMode.IsChecked():
        if data['xtraPeaks']: pkmode = True
    elif data['peaks']:
        pkmode = True
    for item in (G2frame.dataWindow.PeakCopy,
                     G2frame.dataWindow.PeakFit,
                     G2frame.dataWindow.PFOneCycle,
                     G2frame.dataWindow.SeqPeakFit):
        item.Enable(pkmode)
    G2frame.dataWindow.AutoSearch.Enable(not pkmode)
    G2frame.PickTable = []
    rowLabels = []
    PatternId = G2frame.PatternId
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))[0]
    if G2frame.dataWindow.XtraPeakMode.IsChecked():
        G2frame.dataWindow.ClearData()
        topSizer = G2frame.dataWindow.topBox
        parent = G2frame.dataWindow.topPanel
        if 'N' in Inst['Type'][0]:
            lbl = 'List of k-vector/impurity or subgroup defining peaks'
        else:
            lbl = 'List of impurity or subgroup peaks'
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        btn = wx.Button(parent, wx.ID_ANY, 'Switch to Normal peak mode')
        topSizer.Add(btn,0,wx.LEFT,8)
        btn.Bind(wx.EVT_BUTTON,ToggleXtraMode)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        mainSizer =  wx.BoxSizer(wx.VERTICAL)
        G2frame.dataWindow.SetSizer(mainSizer)
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        krowLabels = [str(i+1) for i,j in enumerate(data['xtraPeaks'])]
        kcolLabels = []
        Types = []
        for _,f,l in zip(*G2pwd.getHeaderInfo(Inst['Type'][0])):
            kcolLabels.append(l)
            Types.append(wg.GRID_VALUE_FLOAT + f.replace('%',':').replace('f','').replace('.',','))
            kcolLabels.append('refine')
            Types.append(wg.GRID_VALUE_BOOL)
        T = []
        # put k-vecs in order
        for peak in data['xtraPeaks']:
            T.append(peak[0])
        D = dict(zip(T,data['xtraPeaks']))
        T.sort()
        if 'T' in Inst['Type'][0]:  #want big TOF's first
            T.reverse()
        X = []
        for key in T: X.append(D[key])
        data['xtraPeaks'] = X
        G2frame.GPXtree.SetItemPyData(G2frame.PickId,data)
        kPeakTable = G2G.Table(data['xtraPeaks'],rowLabels=krowLabels,
            colLabels=kcolLabels,types=Types)
        G2frame.dataWindow.currentGrids = []
        reflGrid = G2G.GSGrid(parent=G2frame.dataWindow)
        reflGrid.SetRowLabelSize(45)
        reflGrid.SetTable(kPeakTable, True)
        setBackgroundColors()
        reflGrid.Bind(wg.EVT_GRID_CELL_CHANGED, RefreshPeakGrid)
        reflGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
        reflGrid.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, onCellListDClick)
        reflGrid.AutoSizeColumns(False)
        G2frame.reflGrid = reflGrid
        mainSizer.Add(reflGrid,1,wx.ALL|wx.EXPAND,3)
        G2frame.dataWindow.SetDataSize()
        return

    # Create the header on the Peak Table
    for i in range(len(data['peaks'])): rowLabels.append(str(i+1))
    colLabels = []
    Types = []
    for _,f,l in zip(*G2pwd.getHeaderInfo(Inst['Type'][0])):
        colLabels.append(l)
        if l.startswith('2'):
            Types.append(wg.GRID_VALUE_FLOAT + ':7,3')
        elif l == '00l':
            Types.append(wg.GRID_VALUE_NUMBER)
        else:
            Types.append(wg.GRID_VALUE_FLOAT + f.replace('%',':').replace('f','').replace('.',','))
            colLabels.append('refine')
            Types.append(wg.GRID_VALUE_BOOL)
    # sort peaks & save into tree
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
    G2frame.GPXtree.SetItemPyData(G2frame.PickId,data)

    G2frame.dataWindow.ClearData()
    if 'LF' in Inst['Type'][0]:
        data['LFpeaks'] = []
        for i in range(len(data['peaks'])):
            if len(data['peaks'][i]) == 8:  # need to extend the entry
                if 'Lam' in Inst:
                    lam = Inst['Lam'][0]
                else:
                    lam = (Inst['Lam1'][0] +
                            Inst['I(L2)/I(L1)'][0] * Inst['Lam2'][0]) / (
                                1 + Inst['I(L2)/I(L1)'][0])
                if i == 0:
                    pos = data['peaks'][i][0]
                    data['LaueFringe']['clat'] = data['LaueFringe']['lmin'] * 0.5 * lam / np.sin(pos*np.pi/360)
                    l = data['LaueFringe']['lmin']
                else:
                    l = data['LaueFringe']['clat'] / (0.5 * lam / np.sin(data['peaks'][i][0]*np.pi/360))
                    l = int(l + 0.5)
                    data['peaks'][i][0] = 360 * np.arcsin(0.5 * lam / (data['LaueFringe']['clat']/l)) / np.pi
                data['peaks'][i] += [0., False, 0., False, l, None]
            peakline = copy.deepcopy(data['peaks'][i][2:])
            peakline[-1] = data['peaks'][i][0]
            peakline[-2] = int(peakline[-2]+.5)
            data['LFpeaks'].append(peakline)
        G2frame.PeakTable = G2G.Table(data['LFpeaks'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
    else:
        G2frame.PeakTable = G2G.Table(data['peaks'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
    #G2frame.SetLabel(G2frame.GetLabel().split('||')[0]+' || '+'Peak List')
    G2frame.dataWindow.currentGrids = []
    reflGrid = G2G.GSGrid(parent=G2frame.dataWindow)
    reflGrid.SetRowLabelSize(45)
    reflGrid.SetTable(G2frame.PeakTable, True)
    setBackgroundColors()
    reflGrid.Bind(wg.EVT_GRID_CELL_CHANGED, RefreshPeakGrid)
    reflGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
    reflGrid.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, onCellListDClick)
    reflGrid.AutoSizeColumns(False)
    G2frame.reflGrid = reflGrid
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='List of peaks to fit individually'),0,WACV)
    btn = wx.Button(parent, wx.ID_ANY, 'Switch to Extra Peak mode')
    topSizer.Add(btn,0,wx.LEFT,8)
    btn.Bind(wx.EVT_BUTTON,ToggleXtraMode)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    if 'LF' in Inst['Type'][0]:
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Laue Fringe fitting'))
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Overall parms: '),0,WACV)
        data['LaueFringe'] = data.get('LaueFringe',{})
        data['LaueFringe']['ncell'] = data['LaueFringe'].get('ncell',20)
        data['LaueFringe']['clat'] =  data['LaueFringe'].get('clat',9.0)
        data['LaueFringe']['lmin'] =  data['LaueFringe'].get('lmin',1)
        data['LaueFringe']['clat-ref'] =  data['LaueFringe'].get('clat-ref',False)
        data['LaueFringe']['Show'] =  data['LaueFringe'].get('Show',0)
        data['LaueFringe']['fitRange'] =  data['LaueFringe'].get('fitRange',8.0)
        data['LaueFringe']['fitPowerM'] =  data['LaueFringe'].get('fitPowerM',2.0)
        data['LaueFringe']['fitPowerP'] =  data['LaueFringe'].get('fitPowerP',2.0)
        prmVSizer = wx.BoxSizer(wx.VERTICAL)
        prmSizer = wx.BoxSizer(wx.HORIZONTAL)
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' c='),0,WACV)
        cVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['LaueFringe'],'clat',
                                    typeHint=float,nDig=(10,4),size=(80,-1),
                                    OnLeave=lambda *arg,**kw:RefreshPeakGrid(None))
        prmSizer.Add(cVal,0,WACV)
        cellSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(20,20))
        cellSpin.SetValue(0)
        cellSpin.SetRange(-1,1)
        cellSpin.Bind(wx.EVT_SPIN, ShiftLFc)
        prmSizer.Add(cellSpin,0,WACV)
        cRef = G2G.G2CheckBox(G2frame.dataWindow,'ref',data['LaueFringe'],'clat-ref')
        prmSizer.Add(cRef,0,WACV)
        prmSizer.Add((15,-1))
        siz = G2G.G2SpinWidget(G2frame.dataWindow,data['LaueFringe'] ,'lmin',
                                       'l min')
        prmSizer.Add(siz,0,WACV)
        prmSizer.Add((15,-1))
        siz = G2G.G2SpinWidget(G2frame.dataWindow,data['LaueFringe'] ,'ncell',
                                       'Laue ncell',
                                       onChange=RefreshPeakGrid,onChangeArgs=[None])
        prmSizer.Add(siz,0,WACV)
        # prmSizer.Add((15,-1))
        # prmSizer.Add(wx.StaticText(G2frame.dataWindow,label='  Show '),0,WACV)
        # ch = G2G.EnumSelector(G2frame.dataWindow,data['LaueFringe'],'Show',
        #                             ['None','1','2','3','4','5','6'],list(range(7)),
        #                             OnChange=RefreshPeakGrid)
        # prmSizer.Add(ch,0,WACV)
        # prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' satellites'),0,WACV)
        prmVSizer.Add(prmSizer)

        prmSizer = wx.BoxSizer(wx.HORIZONTAL)
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' fit width'),0,WACV)
        cVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['LaueFringe'],'fitRange',typeHint=float,
            nDig=(6,1),size=(60,-1),xmin=1., xmax=20.,OnLeave=lambda *arg,**kw:RefreshPeakGrid(None))
        prmSizer.Add(cVal,0,WACV)
        prmSizer.Add((15,-1))
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' fit exponent, minus side'),0,WACV)
        cVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['LaueFringe'],'fitPowerM',typeHint=float,
            nDig=(6,1),size=(60,-1),xmin=0.5, xmax=10.,OnLeave=lambda *arg,**kw:RefreshPeakGrid(None))
        prmSizer.Add(cVal,0,WACV)
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' plus side'),0,WACV)
        cVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['LaueFringe'],'fitPowerP',typeHint=float,
            nDig=(6,1),size=(60,-1),xmin=0.5, xmax=10.,OnLeave=lambda *arg,**kw:RefreshPeakGrid(None))
        prmSizer.Add(cVal,0,WACV)
        prmVSizer.Add(prmSizer)

        prmSizer = wx.BoxSizer(wx.HORIZONTAL)
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label='  Show '),0,WACV)
        ch = G2G.EnumSelector(G2frame.dataWindow,data['LaueFringe'],'Show',
            ['None','1','2','3','4','5','6'],list(range(7)),OnChange=RefreshPeakGrid)
        prmSizer.Add(ch,0,WACV)
        prmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' satellites'),0,WACV)
        prmVSizer.Add(prmSizer)
        topSizer.Add(prmVSizer,0,WACV)
        mainSizer.Add(topSizer)
    mainSizer.Add(reflGrid,1,wx.EXPAND,1)
    G2frame.dataWindow.SetDataSize()
    #RefreshPeakGrid(None)

################################################################################
#####  Background
################################################################################

def UpdateBackground(G2frame,data):
    '''respond to selection of PWDR background data tree item.
    '''
    def OnBackFlagCopy(event):
        'Copy background refonement flags to other similar histograms'
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
        FBflag = bool(backDict['background PWDR'][2])
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy bkg ref. flags from\n'+str(hst[5:])+' to...',
            'Copy bkg flags', histList)
        copyList = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            backData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Background'))
            backData[0][1] = copy.copy(flag)
            bkDict = backData[-1]
            if bkDict['nDebye'] == backDict['nDebye']:
                for i,term in enumerate(bkDict['debyeTerms']):
                    term[1::2] = copy.copy(DBflags[i])
            if bkDict['nPeaks'] == backDict['nPeaks']:
                for i,term in enumerate(bkDict['peaksList']):
                    term[1::2] = copy.copy(PKflags[i])
            try:
                backData[1]['background PWDR'][2] = FBflag
            except:
                backData[1]['background PWDR'] = ['',-1.,False]

    def OnBackCopy(event):
        'Copy background functions/values to other similar histograms'
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy bkg params from\n'+str(hst[5:])+' to...',
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.SetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Background'),copy.deepcopy(data))
            CalcBack(Id)

    def OnBackSave(event):
        'Save background values to file'
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Set name to save GSAS-II background parameters file', pth, '',
            'background parameter files (*.pwdrbck)|*.pwdrbck',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                # make sure extension is .pwdrbck
                filename = os.path.splitext(filename)[0]+'.pwdrbck'
                File = open(filename,'w')
                File.write("#GSAS-II background parameter file; do not add/delete items!\n")
                File.write(str(data[0])+'\n')
                for item in data[1]:
                    if item in ['nPeaks','background PWDR','nDebye'] or not len(data[1][item]):
                        File.write(item+':'+str(data[1][item])+'\n')
                    else:
                        File.write(item+':\n')
                        for term in data[1][item]:
                            File.write(str(term)+'\n')
                File.close()
                print ('Background parameters saved to: '+filename)
        finally:
            dlg.Destroy()

    def OnBackLoad(event):
        'load background values from file'
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II background parameters file', pth, '',
            'background parameter files (*.pwdrbck)|*.pwdrbck',wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                newback = [[],{}]
                filename = dlg.GetPath()
                File = open(filename,'r')
                S = File.readline()
                if S[0] == '#':    #skip the heading
                    S = File.readline()     #should contain the std. bck fxn
                newback[0] = eval(S.strip())
                S = File.readline()
                while S and ':' in S:
                    item,vals = S.strip().split(':')
                    if item in ['nPeaks','nDebye']:
                        newback[1][item] = int(vals)
                    elif 'PWDR' in item:
                        newback[1][item] = eval(vals)
                    elif item in ['FixedPoints','debyeTerms','peaksList']:
                        newback[1][item] = []
                        S = File.readline()
                        while S and ':' not in S:
                            newback[1][item].append(eval(S.strip()))
                            S = File.readline()
                        else:
                            continue
                    S = File.readline()
                File.close()
                G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Background'),newback)
        finally:
            dlg.Destroy()
        CalcBack(G2frame.PatternId)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
        wx.CallLater(100,UpdateBackground,G2frame,newback)

    def OnBkgFit(event):
        'Fit background functions to fixed set of background points'

        def SetInstParms(Inst):
            dataType = Inst['Type'][0]
            insVary = []
            insNames = []
            insVals = []
            for parm in Inst:
                insNames.append(parm)
                insVals.append(Inst[parm][1])
                if parm in ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha','A','B','C',
                    'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',] and Inst[parm][2]:
                        Inst[parm][2] = False
#                        insVary.append(parm)
            instDict = dict(zip(insNames,insVals))
            if 'E' not in dataType:  #exclude EDX
                instDict['X'] = max(instDict['X'],0.01)
                instDict['Y'] = max(instDict['Y'],0.01)
            if 'SH/L' in instDict:
                instDict['SH/L'] = max(instDict['SH/L'],0.002)
            return dataType,instDict,insVary

        PatternId = G2frame.PatternId
        controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        # sort the points for convenience and then separate them; extend the range if needed
        if 'FixedPoints' not in background[1] or not len(background[1]['FixedPoints']):
            msg = ("You have not defined any fixed background points. "+
                    "Use the Fixed Points/Add menu item to define points that will be fit."+
                    '\n\nSee the "Fitting the Starting Background using Fixed Points" tutorial for more details.')
            print (msg)
            G2frame.ErrorDialog('No points',msg)
            return
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
        pwddata = G2frame.GPXtree.GetItemPyData(PatternId)[1]
        xBeg = np.searchsorted(pwddata[0],limits[0])
        xFin = np.searchsorted(pwddata[0],limits[1])
        xdata = pwddata[0][xBeg:xFin]
        ydata = si.interp1d(X,Y)(ma.getdata(xdata))
        W = [1]*len(xdata)
        Z = [0]*len(xdata)
        bxye = GetFileBackground(G2frame,data,background,scale=False)[xBeg:xFin]
        if not np.any(bxye):
            bxye = None

        # load instrument and background params
        print (' NB: Any instrument parameter refinement flags will be cleared')
        dataType,insDict,insVary = SetInstParms(inst)
        bakType,bakDict,bakVary = G2pwd.SetBackgroundParms(background)
        # how many background parameters are refined?
        if len(bakVary)*1.5 > len(X):
            msg = ("You are attempting to vary "+str(len(bakVary))+
                   " background terms with only "+str(len(X))+" background points"+
                    "\nAdd more points or reduce the number of terms")
            print (msg)
            G2frame.ErrorDialog('Too few points',msg)
            return

        wx.BeginBusyCursor()
        try:
            G2pwd.DoPeakFit('LSQ',[],background,limits,inst,inst2,
                np.array((xdata,ydata,W,Z,Z,Z)),bxye,prevVaryList=bakVary,controls=controls)
        finally:
            wx.EndBusyCursor()
        # compute the background values and plot them
        parmDict = {}
        bakType,bakDict,bakVary = G2pwd.SetBackgroundParms(background)
        parmDict.update(bakDict)
        parmDict.update(insDict)
        # Note that this generates a MaskedArrayFutureWarning, but these items are not always masked
        pwddata[3][xBeg:xFin] *= 0.
        pwddata[5][xBeg:xFin] *= 0.
        pwddata[4][xBeg:xFin] = G2pwd.getBackground('',parmDict,bakType,dataType,xdata,bxye)[0]
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
        # show the updated background values
        wx.CallLater(100,UpdateBackground,G2frame,data)

    def OnBkgClear(event):
        'Clear fixed points from background'
        if 'FixedPoints' not in data[1]:
            return
        else:
            data[1]['FixedPoints'] = []
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnPeaksMove(event):
        'Move a background peak'
        if not data[1]['nPeaks']:
            G2frame.ErrorDialog('Error','No peaks to move')
            return
        Peaks = {'peaks':[],'sigDict':{}}
        for peak in data[1]['peaksList']:
            Peaks['peaks'].append([peak[0],0,peak[2],0,peak[4],0,peak[6],0])
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List'),Peaks)

    def OnMakeRDF(event):
        'Make a Radial Distribution Function from the background - useful for selecting Debye background positions'
        dlg = RDFDialog(G2frame)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                RDFcontrols = dlg.GetSelection()
            else:
                return
        finally:
            dlg.Destroy()
        PatternId = G2frame.PatternId
        background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
        inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        pwddata = G2frame.GPXtree.GetItemPyData(PatternId)[1]
        auxPlot = G2pwd.MakeRDF(RDFcontrols,background,inst,pwddata)
        for plot in auxPlot:
            XY = np.array(plot[:2])
            if 'D(R)' in plot[2]:
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
            G2frame.GPXtree.SetItemPyData(BackId,data)
            wx.CallLater(100,UpdateBackground,G2frame,data)

        def AfterChange(invalid,value,tc):
            if invalid: return
            CalcBack(G2frame.PatternId)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

        backSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Background function: '),0,WACV)
        bakType = wx.ComboBox(G2frame.dataWindow,value=data[0][0],
                choices=Choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakType.Bind(wx.EVT_COMBOBOX, OnNewType)
        topSizer.Add(bakType)
        topSizer.Add((5,0),0)
        bakRef = wx.CheckBox(G2frame.dataWindow,label=' Refine?')
        bakRef.SetValue(bool(data[0][1]))
        bakRef.Bind(wx.EVT_CHECKBOX, OnBakRef)
        topSizer.Add(bakRef,0,WACV)
        backSizer.Add(topSizer)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Number of coeff.: '),0,WACV)
        bakTerms = wx.ComboBox(G2frame.dataWindow,-1,value=str(data[0][2]),choices=[str(i+1) for i in range(36)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        bakTerms.Bind(wx.EVT_COMBOBOX,OnBakTerms)
        topSizer.Add(bakTerms,0,WACV)
        topSizer.Add((5,0),0)
        backSizer.Add(topSizer)
        backSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Background coefficients:'),0)
        bakSizer = wx.FlexGridSizer(0,5,5,5)
        for i,value in enumerate(data[0][3:]):
            bakVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[0],i+3,nDig=(10,4),OnLeave=AfterChange)
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
            if N == 0:
                CalcBack(G2frame.PatternId)
                G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
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

        def OnCellChange(event):
            CalcBack(G2frame.PatternId)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

        debSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Debye scattering: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Number of terms: '),0,WACV)
        debTerms = wx.ComboBox(G2frame.dataWindow,-1,value=str(data[1]['nDebye']),choices=[str(i) for i in range(21)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        debTerms.Bind(wx.EVT_COMBOBOX,OnDebTerms)
        topSizer.Add(debTerms,0,WACV)
        topSizer.Add((5,0),0)
        debSizer.Add(topSizer)
        if data[1]['nDebye']:
            debSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Debye diffuse terms:'),0)
            rowLabels = []
            for i in range(len(data[1]['debyeTerms'])): rowLabels.append(str(i))
            colLabels = ['A','refine','R','refine','U','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            debyeTable = G2G.Table(data[1]['debyeTerms'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            debyeGrid = G2G.GSGrid(parent=G2frame.dataWindow)
            debyeGrid.SetTable(debyeTable, True)
            debyeGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            debyeGrid.Bind(wg.EVT_GRID_CELL_CHANGED,OnCellChange)
            debyeGrid.AutoSizeColumns(False)
            debSizer.Add(debyeGrid)
        return debSizer

    def PeaksSizer():

        def OnPeaks(event):
            'Respond to a change in the number of background peaks'
            data[1]['nPeaks'] = int(peaks.GetValue())
            M = len(data[1]['peaksList'])
            N = data[1]['nPeaks']
            if N > M:       #add terms
                for i in range(M,N):
                    data[1]['peaksList'].append([1.0,False,1.0,False,0.10,False,0.10,False])
            elif N < M:     #delete terms
                for i in range(N,M):
                    del(data[1]['peaksList'][-1])
            if N == 0:
                CalcBack(G2frame.PatternId)
                G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
            # this callback is crashing wx when there is an open
            # peaksGrid cell editor, at least on Mac. Code below
            # should fix this, but it does not.
            # https://stackoverflow.com/questions/64082199/wxpython-grid-destroy-with-open-celleditor-crashes-python-even-with-disablece
            if peaksGrid and peaksGrid.IsCellEditControlEnabled():
                # complete any grid edits in progress
                #print('closing')
                peaksGrid.HideCellEditControl()
                peaksGrid.DisableCellEditControl()
                #wx.CallLater(100,peaksGrid.Destroy) # crashes python
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

        def OnCellChange(event):
            CalcBack(G2frame.PatternId)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

        peaksSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Peaks in background: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Number of peaks: '),0,WACV)
        peaks = wx.ComboBox(G2frame.dataWindow,-1,value=str(data[1]['nPeaks']),choices=[str(i) for i in range(30)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        peaks.Bind(wx.EVT_COMBOBOX,OnPeaks)
        topSizer.Add(peaks,0,WACV)
        topSizer.Add((5,0),0)
        peaksSizer.Add(topSizer)
        G2frame.dataWindow.currentGrids = []
        peaksGrid = None
        if data[1]['nPeaks']:
            peaksSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Peak list:'),0)
            rowLabels = []
            for i in range(len(data[1]['peaksList'])): rowLabels.append(str(i))
            colLabels = ['pos','refine','int','refine','sig','refine','gam','refine']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,3',wg.GRID_VALUE_BOOL,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_BOOL]
            peaksTable = G2G.Table(data[1]['peaksList'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
            peaksGrid = G2G.GSGrid(parent=G2frame.dataWindow)
            peaksGrid.SetRowLabelSize(45)
            peaksGrid.SetTable(peaksTable, True)
            peaksGrid.Bind(wx.EVT_KEY_DOWN, KeyEditPeakGrid)
            peaksGrid.Bind(wg.EVT_GRID_CELL_CHANGED,OnCellChange)
            peaksGrid.AutoSizeColumns(False)
            peaksSizer.Add(peaksGrid)
        return peaksSizer

    def BackFileSizer():

        def OnBackPWDR(event):
            data[1]['background PWDR'][0] = back.GetValue()
            if len(data[1]['background PWDR'][0]):
                curHist = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data[1]['background PWDR'][0])
                if not Id:
                    G2G.G2MessageBox(G2frame,'Histogram not found -- how did this happen?','Missing histogram')
                    back.SetValue('')
                    data[1]['background PWDR'][0] = back.GetValue()
                    return
                bkgHist = G2frame.GPXtree.GetItemPyData(Id)
                if len(bkgHist[1][0]) != len(curHist[1][0]):
                    G2G.G2MessageBox(G2frame,'Histogram have different lengths','Mismatched histograms')
                    back.SetValue('')
                    data[1]['background PWDR'][0] = back.GetValue()
                    return
            else:
                data[1]['background PWDR'][2] = False
            CalcBack()
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
            wx.CallLater(100,UpdateBackground,G2frame,data)

        def AfterChange(invalid,value,tc):
            if invalid: return
            CalcBack()
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

        def OnBackFit(event):
            data[1]['background PWDR'][2] = not data[1]['background PWDR'][2]

        fileSizer = wx.BoxSizer(wx.VERTICAL)
        fileSizer.Add((-1,5))
        backSizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(G2frame.dataWindow, wx.ID_ANY,'Fit to fixed bkg')
        backSizer.Add(btn,0,wx.RIGHT,3)
        btn.Enable(len(data[1].get('FixedPoints',[])) > 5)
        btn.Bind(wx.EVT_BUTTON,OnBkgFit)

        btn = wx.Button(G2frame.dataWindow, wx.ID_ANY,'Compute auto background')
        backSizer.Add(btn,0,wx.RIGHT,3)
        btn.Bind(wx.EVT_BUTTON,onAutoBack)

        btn = wx.Button(G2frame.dataWindow, wx.ID_ANY,'Copy auto background')
        backSizer.Add(btn,0,wx.RIGHT,3)
        data[1]['autoPrms'] = data[1].get('autoPrms',{})
        btn.Enable(bool(data[1]['autoPrms'].get('Mode')))
        btn.Bind(wx.EVT_BUTTON,copyAutoBack)
        fileSizer.Add(backSizer)
        fileSizer.Add((-1,5))

        fileSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Fixed background histogram (for point-by-point subtraction):'),0)
        if 'background PWDR' not in data[1]:
            data[1]['background PWDR'] = ['',-1.,False]
        backSizer = wx.BoxSizer(wx.HORIZONTAL)
        Choices = ['',]+G2gd.GetGPXtreeDataNames(G2frame,['PWDR',])
        Source = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        Choices.pop(Choices.index(Source))
        back = wx.ComboBox(parent=G2frame.dataWindow,value=data[1]['background PWDR'][0],choices=Choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        back.Bind(wx.EVT_COMBOBOX,OnBackPWDR)
        backSizer.Add(back)
        backSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' multiplier'),0,WACV)
        backMult = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[1]['background PWDR'],1,nDig=(10,3),OnLeave=AfterChange)
        backSizer.Add(backMult,0,WACV)
        if len(data[1]['background PWDR'][0]):
            backFit = wx.CheckBox(G2frame.dataWindow,label=' Refine?')
            backFit.SetValue(data[1]['background PWDR'][2])
            backFit.Bind(wx.EVT_CHECKBOX, OnBackFit)
            backSizer.Add(backFit,0,WACV)
        fileSizer.Add(backSizer)
        return fileSizer

    def onAutoBack(event):
        '''Open a window for auto background computation
        '''
        bkgdict = data[1]
        xydata = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1]
        autoBackground(G2frame)
        if bkgdict['autoPrms']['Mode'] == 'fixed':
            xydata[4] = G2pwd.autoBkgCalc(bkgdict,xydata[1])
            addAutoBack(G2frame,data,xydata)
            wx.CallAfter(UpdateBackground,G2frame,data)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
        elif bkgdict['autoPrms']['Mode'] == 'fit':
            xydata[4] = G2pwd.autoBkgCalc(bkgdict,xydata[1])
            npts = len(xydata[0])
            bkgdict['FixedPoints'] = [i for i in zip(
                xydata[0].data[::npts//100],
                xydata[4].data[::npts//100])]
            OnBkgFit(event)
        else:
            wx.CallAfter(UpdateBackground,G2frame,data)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def copyAutoBack(event):
        '''reproduce the auto background computation on selected
        other histograms
        '''
        savePatternId = G2frame.PatternId
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        autoBkgDict = data[1].get('autoPrms')
        if not autoBkgDict.get('Mode'):
            G2frame.ErrorDialog('No auto bkg setting','This is unexpected, no auto bkg parms for histogram '+hst,G2frame)
            return
        elif autoBkgDict['Mode'] == 'fit':
            txt = 'set fixed points from auto bkg calc'
        else:
            txt = 'use auto bkg calc to define Fixed Bkg histogram'
        histList = [i for i in GetHistsLikeSelected(G2frame)
                        if 'Autobkg for' not in i]
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,
                        f'Select histogram(s) to {txt} based on {hst}',
                        'Compute auto bkg for...', histList)
        try:
            copyList = []
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.PatternId = Id
            xydata = G2frame.GPXtree.GetItemPyData(Id)[1]
            G2frame.GPXtree.SetItemPyData
            bkgId = G2gd.GetGPXtreeItemId(G2frame,Id,'Background')
            G2frame.GPXtree.SetItemPyData(bkgId,copy.deepcopy(data))
            itemData = G2frame.GPXtree.GetItemPyData(bkgId)
            if autoBkgDict['Mode'] == 'fixed':
                xydata[4] = G2pwd.autoBkgCalc(itemData[1],xydata[1])
                addAutoBack(G2frame,itemData,xydata)
            elif autoBkgDict['Mode'] == 'fit':
                xydata[4] = G2pwd.autoBkgCalc(itemData[1],xydata[1])
                npts = len(xydata[0])
                itemData[1]['FixedPoints'] = [i for i in zip(
                xydata[0].data[::npts//100],
                xydata[4].data[::npts//100])]
                OnBkgFit(event)
        G2frame.PatternId = savePatternId
        wx.CallAfter(UpdateBackground,G2frame,data)

    def CalcBack(PatternId=G2frame.PatternId):
        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
        inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
        backData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
        dataType = inst['Type'][0]
        insDict = {inskey:inst[inskey][1] for inskey in inst}
        parmDict = {}
        bakType,bakDict,bakVary = G2pwd.SetBackgroundParms(data)
        parmDict.update(bakDict)
        parmDict.update(insDict)
        pwddata = G2frame.GPXtree.GetItemPyData(PatternId)
        xBeg = np.searchsorted(pwddata[1][0],limits[0])
        xFin = np.searchsorted(pwddata[1][0],limits[1])
        fixBack = GetFileBackground(G2frame,pwddata[1],backData,scale=False)
        pwddata[1][4][xBeg:xFin] = G2pwd.getBackground('',parmDict,bakType,dataType,pwddata[1][0][xBeg:xFin],fixBack[xBeg:xFin])[0]

    # UpdateBackground execution starts here
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu) # needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.BackMenu)
    if len(data) < 2:       #add Debye diffuse & peaks scattering here
        data.append({'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[],'background PWDR':['',1.0,False]})
    if 'nPeaks' not in data[1]:
        data[1].update({'nPeaks':0,'peaksList':[],'background PWDR':['',1.0,False]})
    if 'background PWDR' not in data[1]:
        data[1].update({'background PWDR':['',1.0,False]})
    elif len(data[1]['background PWDR']) < 3:
        data[1]['background PWDR'].append(False)
    G2frame.dataWindow.currentGrids = []
    G2frame.Bind(wx.EVT_MENU,OnBackCopy,id=G2G.wxID_BACKCOPY)
    G2frame.Bind(wx.EVT_MENU,OnBackFlagCopy,id=G2G.wxID_BACKFLAGCOPY)
    G2frame.Bind(wx.EVT_MENU,OnBackSave,id=G2G.wxID_BACKSAVE)
    G2frame.Bind(wx.EVT_MENU,OnBackLoad,id=G2G.wxID_BACKLOAD)
    G2frame.Bind(wx.EVT_MENU,OnPeaksMove,id=G2G.wxID_BACKPEAKSMOVE)
    G2frame.Bind(wx.EVT_MENU,OnMakeRDF,id=G2G.wxID_MAKEBACKRDF)
    G2frame.Bind(wx.EVT_MENU,OnBkgFit,id=G2frame.dataWindow.wxID_BackPts['Fit'])
    G2frame.Bind(wx.EVT_MENU,OnBkgClear,id=G2frame.dataWindow.wxID_BackPts['Clear'])
    BackId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background')
    Choices = ['chebyschev','chebyschev-1','cosine','Q^2 power series','Q^-2 power series','lin interpolate','inv interpolate','log interpolate']
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Background used in refinement'),0,WACV)
    # add help button to bring up help web page - at right side of window
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add(BackSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(DebyeSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(PeaksSizer())
    mainSizer.Add((0,5),0)
    mainSizer.Add(BackFileSizer())
    G2frame.dataWindow.SetDataSize()

def addAutoBack(G2frame,data,xydata):
    '''Create a new histogram for the computed auto background and place
    as the fixed background histogram
    '''
    bkgHistName = 'PWDR Autobkg for '+G2frame.GPXtree.GetItemText(G2frame.PatternId)[5:]
    # if histogram exists we should probably reuse it, but for now, just create a new one
    bkgHistName = G2obj.MakeUniqueLabel(bkgHistName,G2frame.GetHistogramNames('PWDR'))

    Ymin = min(xydata[4])
    Ymax = max(xydata[4])
    d = copy.deepcopy(G2frame.GPXtree.GetItemPyData(G2frame.PatternId))
    d[0] = {'wtFactor': 1.0, 'Dummy': False,
                  'ranId': ran.randint(0, sys.maxsize),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'Yminmax': [Ymin, Ymax]}
    d[1][1] = xydata[4]
    d[1][2] = np.ones_like(xydata[4])
    d[1][3] = np.zeros_like(xydata[4])
    d[1][4] = np.zeros_like(xydata[4])
    d[1][5] = np.zeros_like(xydata[4])
    NewId = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=bkgHistName)
    G2frame.GPXtree.SetItemPyData(NewId,d)

    item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.PatternId)
    while item:
        nam = G2frame.GPXtree.GetItemText(item)
        if nam == 'Comments':
            d = [' # background generated with Autobkg']
        elif nam == 'Background':
            d = [['chebyschev-1',True,3,1.0,0.0,0.0],
                     {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[],
                          'background PWDR':['',1.0,False],'FixedPoints':[]}]
        elif nam == 'Peak List':
            d = {'sigDict':{},'peaks':[]}
        elif nam == 'Index Peak List':
            d = [[], []]
        elif nam == 'Unit Cells List':
            d = []
        elif nam == 'Reflection Lists':
            d = {}
        else:
            d = copy.deepcopy(G2frame.GPXtree.GetItemPyData(item))
        G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(parent=NewId,text=nam),d)
        item, cookie = G2frame.GPXtree.GetNextChild(G2frame.PatternId, cookie)


    bId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Background')
    data = G2frame.GPXtree.GetItemPyData(bId)
    # set fixed bkg & turn off computed background
    data[1]['background PWDR'] = [bkgHistName, 1.0, False]
    data[0][1] = False
    data[0][3:] = data[0][2]*[0.]
    for p in data[1]['peaksList']:
        p[2] = 0.
        p[1::2] = 4*[False]
    for p in data[1]['debyeTerms']:
        p[0] = 0
        p[1::2] = 3*[False]
    G2frame.GetUsedHistogramsAndPhasesfromTree() # reindex

# Autobackground Dialog
class autoBackground(wx.Dialog):
    '''Create a file selection widget for setting background with
    pybaselines, as requested by James Feng.

    :param wx.Frame G2frame: reference to the main GSAS-II frame.

    '''
    def __init__(self,G2frame,*args,**kwargs):
        self.G2frame = G2frame
        bId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Background')
        data = G2frame.GPXtree.GetItemPyData(bId)
        self.bkgdict = data[1]
        self.xydata = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1]
        npts = len(self.xydata[0])
        # add auto bkg to background prms dict
        self.bkgdict['autoPrms'] = self.bkgdict.get('autoPrms',{})
        self.bkgdict['autoPrms']['opt'] = self.bkgdict['autoPrms'].get('opt',0)
        logLam =  min(10,float(int(10*np.log10(npts)**1.5)-9.5)/10.)
        self.bkgdict['autoPrms']['logLam'] = self.bkgdict['autoPrms'].get('logLam',logLam)
        self.bkgdict['autoPrms']['Mode'] = None
        maxLam = min(15.,1.*int(3*self.bkgdict['autoPrms']['logLam']+0.9))
        # save starting point info
        self.startingBackground = copy.deepcopy(self.xydata[4])
        # start process
        wx.Dialog.__init__(self, parent=G2frame,
                                 style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.CenterOnParent()

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self,label=' Compute autobackground'),0)
        choices = ['arpls','iarpls']
        subSiz = wx.BoxSizer(wx.HORIZONTAL)
        subSiz.Add(wx.StaticText(self,label='Computation option'))
        for w in G2G.G2RadioButtons(self,self.bkgdict['autoPrms'],'opt',choices,
                                    OnChange=self._calcBkg):
            subSiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL,0)
        mainSizer.Add(subSiz)
        siz = G2G.G2SliderWidget(self,self.bkgdict['autoPrms'],
                                'logLam','log(Lambda)',
                                1.,maxLam,100,self._calcBkg)
        mainSizer.Add(siz)

        subSiz = wx.BoxSizer(wx.HORIZONTAL)
        subSiz.Add((-1,-1),1,wx.EXPAND,1)
        btn = wx.Button(self, wx.ID_CLOSE, label='Set Fixed\nPoints && Fit')
        btn.Bind(wx.EVT_BUTTON,lambda event: self.EndModal(wx.ID_CLOSE))
        subSiz.Add(btn)
        subSiz.Add((5,-1))
        btn = wx.Button(self, wx.ID_OK, label='Define Fixed\nBkg histogram')
        btn.Bind(wx.EVT_BUTTON,lambda event: self.EndModal(wx.ID_OK))
        subSiz.Add(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        subSiz.Add((5,-1))
        subSiz.Add(btn,0,wx.CENTER)
        subSiz.Add((-1,-1),1,wx.EXPAND,1)
        mainSizer.Add((-1,5))
        mainSizer.Add(subSiz,0,wx.EXPAND)
        mainSizer.Add((-1,5))
        self.SetSizer(mainSizer)
        mainSizer.Fit(self)
        self._calcBkg()
        res = self.ShowModal()
        if res == wx.ID_CLOSE:
            self.bkgdict['autoPrms']['Mode'] = 'fit'
        elif res == wx.ID_OK:
            self.bkgdict['autoPrms']['Mode'] = 'fixed'
        else:
            # restore the background to the starting values
            self.xydata[4] = self.startingBackground
            self.bkgdict['autoPrms']['Mode'] = None

    def _calcBkg(self,event=None):
        '''respond to a change in the background parameters by recomputing
        the auto background
        '''
        self.xydata[4] = G2pwd.autoBkgCalc(self.bkgdict,self.xydata[1].data)
        G2pwpl.PlotPatterns(self.G2frame,plotType='PWDR')

################################################################################
#####  Limits
################################################################################

def UpdateLimitsGrid(G2frame, data,datatype):
    '''respond to selection of PWDR Limits data tree item.
    Allows setting of limits and excluded regions in a PWDR data set
    '''
    def AfterChange(invalid,value,tc):
        if invalid: return
        datatype = G2frame.GPXtree.GetItemText(G2frame.PatternId)[:4]
        wx.CallAfter(G2pwpl.PlotPatterns,G2frame,newPlot=False,plotType=datatype)  #unfortunately this resets the plot width

    def LimitSizer():
        limits = wx.FlexGridSizer(0,3,0,5)
        labels = ['Tmin','Tmax']
        for i in [0,1]:
            limits.Add(wx.StaticText(G2frame.dataWindow,
                label=' Original {} {:.4f}'.format(labels[i],data[0][i])),0,WACV)
            limits.Add(wx.StaticText(G2frame.dataWindow,label=' New: '),0,WACV)
            limits.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[1],i,  \
                xmin=data[0][0],xmax=data[0][1],nDig=(10,4),typeHint=float,OnLeave=AfterChange))
        return limits

    def ExclSizer():

        def OnDelExcl(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            del(data[item+2])
            G2pwpl.PlotPatterns(G2frame,newPlot=False,plotType=datatype)
            wx.CallAfter(UpdateLimitsGrid,G2frame,data,datatype)

        Indx = {}
        excl = wx.FlexGridSizer(0,3,0,5)
        excl.Add(wx.StaticText(G2frame.dataWindow,label=' From: '),0,WACV)
        excl.Add(wx.StaticText(G2frame.dataWindow,label=' To: '),0,WACV)
        excl.Add(wx.StaticText(G2frame.dataWindow,label=' Delete?: '),0,WACV)
        for Id,item in enumerate(data[2:]):
            for i in [0,1]:
                excl.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,item,i,  \
                    xmin=data[0][0],xmax=data[0][1],nDig=(10,4),typeHint=float,OnLeave=AfterChange))
            delExcl = wx.CheckBox(G2frame.dataWindow,label='')
            Indx[delExcl.GetId()] = Id
            delExcl.Bind(wx.EVT_CHECKBOX,OnDelExcl)
            excl.Add(delExcl,0,WACV)
        return excl

    def onSetLimExcl(event):
        txt = event.EventObject.FindItemById(event.GetId()).GetItemLabel() # gets menu command text
        if 'lower' in txt:
            G2frame.ifSetLimitsMode = 1
            G2frame.CancelSetLimitsMode.Enable(True)
        elif 'upper' in txt:
            G2frame.ifSetLimitsMode = 2
            G2frame.CancelSetLimitsMode.Enable(True)
        elif 'excl' in txt:
            G2frame.ifSetLimitsMode = 3
            G2frame.CancelSetLimitsMode.Enable(True)
        else:
            G2frame.ifSetLimitsMode = 0
            G2frame.CancelSetLimitsMode.Enable(False)
        G2frame.plotFrame.Raise()
        G2pwpl.PlotPatterns(G2frame,newPlot=False,plotType=datatype)

    def OnLimitCopy(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy limits from\n'+str(hst[5:])+' to...',
            'Copy limits', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    item = histList[i]
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    G2frame.GPXtree.SetItemPyData(
                        G2gd.GetGPXtreeItemId(G2frame,Id,'Limits'),copy.deepcopy(data))
        finally:
            dlg.Destroy()

    def Draw():
        G2frame.dataWindow.ClearData()
        topSizer = G2frame.dataWindow.topBox
        parent = G2frame.dataWindow.topPanel
        topSizer.Add(wx.StaticText(parent,label=' Data range to be used in fits'),0,WACV)
        # add help button to bring up help web page - at right side of window
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        mainSizer =  wx.BoxSizer(wx.VERTICAL)
        G2frame.dataWindow.SetSizer(mainSizer)
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        mainSizer.Add((5,5))
        mainSizer.Add(LimitSizer())
        if len(data)>2:
            mainSizer.Add((0,5),0)
            mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Excluded regions:'))
            mainSizer.Add(ExclSizer())
        G2frame.dataWindow.SetDataSize()

    # start of UpdateLimitsGrid
    G2frame.ifSetLimitsMode = 0
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu) # needed below
    if 'P' in datatype:                   #powder data menu commands
        G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.LimitMenu)
        G2frame.CancelSetLimitsMode.Enable(False)
        G2frame.Bind(wx.EVT_MENU,OnLimitCopy,id=G2G.wxID_LIMITCOPY)
        for n in (G2G.wxID_ADDEXCLREGION, G2G.wxID_SETLOWLIMIT,
                  G2G.wxID_SETTOPLIMIT, G2G.wxID_STOPSETLIMIT):
            G2frame.Bind(wx.EVT_MENU,onSetLimExcl,id=n)
    elif 'A' in datatype or 'R' in datatype:                   #SASD & REFD data menu commands
        G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SASDLimitMenu)
        G2frame.CancelSetLimitsMode.Enable(False)
        G2frame.Bind(wx.EVT_MENU,OnLimitCopy,id=G2G.wxID_SASDLIMITCOPY)
    Draw()

################################################################################
#####  Instrument parameters
################################################################################

def UpdateInstrumentGrid(G2frame,data):
    '''respond to selection of PWDR/SASD/REFD Instrument Parameters
    data tree item.
    '''
    if 'Bank' not in data:  #get it from name; absent for default parms selection
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        if 'Bank' in hst:
            bank = int(hst.split('Bank')[1].split('_')[0])
            data['Bank'] = [bank,bank,0]
        else:
            data['Bank'] = [1,1,0]

    def keycheck(keys):
        good = []
        for key in keys:
            if key in ['Type','Bank','U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha','A','B','C',
                'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','Polariz.','alpha-0','alpha-1',
                'Lam','Azimuth','2-theta','fltPath','difC','difA','difB','Zero','Lam1','Lam2','XE','YE','ZE','WE']:
                good.append(key)
        return good

    def updateData(inst,ref):
        Data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[0]
        for item in Data:
            try:
                Data[item] = [Data[item][0],inst[item],ref[item]]
            except KeyError:
                try:
                    Data[item] = [Data[item][0],inst[item]]
                except KeyError:
                    pass        #skip 'Polariz.' for N-data

    def RefreshInstrumentGrid(event,doAnyway=False):
        if doAnyway or event.GetRow() == 1:
            peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
            newpeaks = []
            for peak in peaks['peaks']:
                newpeaks.append(G2mth.setPeakparms(data,Inst2,peak[0],peak[2]))
            peaks['peaks'] = newpeaks
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List'),peaks)

    def OnCalibrate(event):
        Pattern = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
        xye = ma.array(ma.getdata(Pattern[1]))
        cw = np.diff(xye[0])
        IndexPeaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
        Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Sample Parameters'))
        if 'Debye' not in Sample['Type']:
            G2frame.ErrorDialog('Cannot calibrate','Only apropriate for Debye-Scherrer geometry')
            return
        if not len(IndexPeaks[0]):
            G2frame.ErrorDialog('Cannot calibrate','Index Peak List empty')
            return
        if not np.any(IndexPeaks[1]):
            G2frame.ErrorDialog('Cannot calibrate','Peak positions not refined')
            return False
        Ok = False
        for peak in IndexPeaks[0]:
            if peak[2] and peak[3]:
                Ok = True
        if not Ok:
            G2frame.ErrorDialog('Cannot calibrate','Index Peak List not indexed')
            return
        if G2pwd.DoCalibInst(IndexPeaks,data,Sample):
            UpdateInstrumentGrid(G2frame,data)
            const = 0.0
            if 'C' in data['Type'][0] or 'B' in data['Type'][0]:
                const = 18.e-2/(np.pi*Sample['Gonio. radius'])
                # const = 10**-3/Sample['Gonio. radius']
            XY = []
            Sigs = []
            for ip,peak in enumerate(IndexPeaks[0]):
                shft = 0.0
                if peak[2] and peak[3]:
                    binwid = cw[np.searchsorted(xye[0],peak[0])]
                    if const:
                        shft = -const*(Sample['DisplaceX'][0]*npcosd(peak[0])+Sample['DisplaceY'][0]*npsind(peak[0]))
                    XY.append([peak[-1],peak[0]-shft,binwid])
                    Sigs.append(IndexPeaks[1][ip])
            if len(XY):
                XY = np.array(XY)
                G2plt.PlotCalib(G2frame,data,XY,Sigs,newPlot=True)
        else:
            G2frame.ErrorDialog('Cannot calibrate','Nothing selected for refinement or refinement failed')

    def OnLoad(event):
        '''Loads instrument parameters from a G2 .instprm file
        in response to the Instrument Parameters-Operations/Load Profile menu
        If instprm file has multiple banks each with header #Bank n: ..., this
        finds matching bank no. to load - rejects nonmatches.

        Note uses ReadPowderInstprm (GSASIIdataGUI.py) to read .instprm fil
        '''

        def GetDefaultParms(rd):
            '''Solicits from user a default set of parameters & returns Inst parm dict
            param: rd: importer data structure
            returns: dict: Instrument parameter dictionary
            '''
            import defaultIparms as dI
            sind = lambda x: math.sin(x*math.pi/180.)
            tand = lambda x: math.tan(x*math.pi/180.)
            while True: # loop until we get a choice
                choices = []
                head = 'Select from default instrument parameters'

                for l in dI.defaultIparm_lbl:
                    choices.append('Defaults for '+l)
                res = G2G.BlockSelector(choices,ParentFrame=G2frame,title=head,
                    header='Select default inst parms',useCancel=True)
                if res is None: return None
                if 'Generic TOF' in choices[res]:
                    dlg = G2G.MultiDataDialog(G2frame,title='Generic TOF detector bank',
                        prompts=['Total FP','2-theta',],values=[25.0,150.,],
                            limits=[[6.,200.],[5.,175.],],formats=['%6.2f','%6.1f',])
                    if dlg.ShowModal() == wx.ID_OK: #strictly empirical approx.
                        FP,tth = dlg.GetValues()
                        difC = 505.632*FP*sind(tth/2.)
                        sig1 = 50.+2.5e-6*(difC/tand(tth/2.))**2
                        bet1 = .00226+7.76e+11/difC**4
                        Inst = G2frame.ReadPowderInstprm(dI.defaultIparms[res],bank,rd)
                        Inst[0]['difC'] = [difC,difC,0]
                        Inst[0]['sig-1'] = [sig1,sig1,0]
                        Inst[0]['beta-1'] = [bet1,bet1,0]
                        return Inst    #this is [Inst1,Inst2] a pair of dicts
                    dlg.Destroy()
                else:
                    inst1,inst2 = G2frame.ReadPowderInstprm(dI.defaultIparms[res],bank,rd)
                    return [inst1,inst2]
                if 'lab data' in choices[res]:
                    rd.Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
                        'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
                else:
                    rd.Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],
                        'DisplaceY':[0.,False]})

        data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[0]
        bank = data['Bank'][0]
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II instrument parameters file', pth, '',
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Found = False
                filename = dlg.GetPath()
                try:
                    File = open(dlg.GetPath(),'r')
                    instLines = File.readlines()
                    File.close()
                    rd = G2obj.ImportPowderData('Dummy')
                    rd.Sample = G2frame.GPXtree.GetItemPyData(
                        G2gd.GetGPXtreeItemId(
                            G2frame,G2frame.PatternId,'Sample Parameters'))
                    instvals = G2frame.ReadPowderInstprm(instLines, bank, rd)
                    Found = True
                except:
                    print('instprm read failed')
                if Found:
                    Inst,Inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Instrument Parameters'))
                    if 'Bank' not in Inst:  #patch for old .instprm files - may cause faults for TOF data
                        Inst['Bank'] = [1,1,0]
                    Inst.update(instvals[0])
                    if 'B' in Inst['Type'][1] and 'SH/L' in Inst:
                        del Inst['SH/L']
                    RefreshInstrumentGrid(event,doAnyway=True)          #to get peaks updated
                    UpdateInstrumentGrid(G2frame,data)
                else:
                    G2frame.ErrorDialog('No match','Bank %d not in %s'%(bank,filename),G2frame)
            else:
                rd = G2obj.ImportPowderData('Dummy')
                rd.Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Sample Parameters'))
                try:
                    data = GetDefaultParms(rd)[0]
                except TypeError:   #Cancel - got None
                    pass
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
        dlg = wx.FileDialog(G2frame, 'Set name to save GSAS-II instrument parameters file', pth, '',
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Sample = G2frame.GPXtree.GetItemPyData(
                    G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId,
                                              'Sample Parameters'))
                filename = dlg.GetPath()
                # make sure extension is .instprm
                filename = os.path.splitext(filename)[0]+'.instprm'
                File = open(filename,'w')
                G2fil.WriteInstprm(File, data, Sample)
                File.close()
                print ('Instrument parameters saved to: '+filename)
        finally:
            dlg.Destroy()

    def OnSaveAll(event):
        '''Respond to the Instrument Parameters Operations/Save all Profile menu & writes
        selected inst parms. across multiple banks into a single file
        Each block starts with #Bank n: GSAS-II instrument... where n is bank no.
        item: writes parameters from selected PWDR entries to a .instprm file
        '''
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        histList.insert(0,hst)
        saveList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Save instrument parameters from',
            'Save instrument parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                saveList = [histList[i] for i in dlg.GetSelections()]
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
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
                    inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
                    Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                    if 'Bank' not in inst:  #patch
                        bank = 1
                        if 'Bank' in hist:
                            bank = int(hist.split('Bank')[1])
                        inst['Bank'] = [bank,bank,0]
                    G2fil.WriteInstprm(File, inst, Sample, inst['Bank'][0])
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
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        keys = list(data.keys())
        try:
            keys.remove('Source')
        except ValueError:
            pass
        flags = dict(zip(keys,[data[key][2] for key in keys]))
        instType = data['Type'][0]
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy inst ref. flags from\n'+hst[5:],
            'Copy refinement flags', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            instData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
            if 'Bank' not in instData:
                instData['Bank'] = [1,1,0]
            if len(data) == len(instData) and instType == instData['Type'][0]:   #don't mix data types or lam & lam1/lam2 parms!
                for item in instData:
                    if item not in ['Source',]:
                        instData[item][2] = copy.copy(flags[item])
            else:
                print (item+' not copied - instrument parameters not commensurate')

    def OnInstCopy(event):
        #need fix for dictionary
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        copyData = copy.deepcopy(data)
        if 'E' not in data['Type'][0]:
            del copyData['Azimuth'] #not to be copied!
        instType = data['Type'][0]
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy inst params from\n'+hst,
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            instData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
            if 'Bank' not in data:
                data['Bank'] = [1,1,0]
            if 'Source' not in data:
                data['Source'] = ['','']
            if 'Bank' not in instData:
                instData['Bank'] = [1,1,0]
            if 'Source' not in instData:
                instData['Source'] = ['','']
            if 'Z' not in instData:
                instData['Z'] = [0.,0.,False]
            if len(data) == len(instData) and instType == instData['Type'][0]:  #don't mix data types or lam & lam1/lam2 parms!
                instData.update(copyData)
            else:
                if len(data) != len(instData):
                    print (item+' not copied - %d instrument parameters do not match source # %d'%(len(instData),len(data)))
                else:
                    print (item+' not copied - instrument type %s does not match source type %s'%(instData['Type'][0],instType))

    def AfterChange(invalid,value,tc):
        if invalid: return
        updateData(insVal,insRef)
        G2plt.PlotPeakWidths(G2frame)

    def AfterChangeEC(invalid,value,tc):
        '''for SEC data only; converts electrn energy in keV to wavelength
        '''
        if invalid: return
        if value > 10.:
            value *= 1000.      #keV -> eV
            value = 12.2639/np.sqrt(value+value**2*0.97845e-6)
            tc.SetValue(value)
            insVal.update({'Lam':value})
        updateData(insVal,insRef)

    def NewProfile(invalid,value,tc):
        if invalid: return
        G2plt.PlotPeakWidths(G2frame)
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
        wx.CallAfter(UpdateInstrumentGrid,G2frame,data)

    def OnInstMult(event):
        'If checked or unchecked, redisplay window'
        wx.CallAfter(UpdateInstrumentGrid,G2frame,data)

    def lblWdef(lbl,dec,val):
        'Label parameter showing the default value'
        fmt = "%15."+str(dec)+"f"
        return " " + lbl + " (" + (fmt % val).strip() + "): "

    def RefineBox(item):
        'Define a refine checkbox with binding'
        wid = wx.CheckBox(G2frame.dataWindow,label='')
        wid.SetValue(bool(insRef[item]))
        RefObj[wid.GetId()] = item
        wid.Bind(wx.EVT_CHECKBOX, OnItemRef)
        return wid

    def OnLamPick(event):
        'After selection of lab. x-ray source type'
        data['Source'][1] = lamType = event.GetEventObject().GetValue()
        if 'P' in insVal['Type']:
            insVal['Lam1'] = G2elem.waves[lamType][0]
            insVal['Lam2'] = G2elem.waves[lamType][1]
        elif 'S' in insVal['Type']: #and
            try:
                insVal['Lam'] = G2elem.meanwaves[lamType]
                data['Type'][0] = 'SXC'
                insVal['Type'] = 'SXC'
            except KeyError:
                if 'synch' in lamType:
                    insVal['Lam'] = 1.0  #typical?
                    data['Type'][0] = 'SXC'
                    insVal['Type'] = 'SXC'
                elif 'micro' in lamType:
                    insVal['Lam'] = 0.0251 # @200keV
                    data['Type'][0] = 'SEC'
                    insVal['Type'] = 'SEC'      #change to electron diffraction
        updateData(insVal,insRef)
        wx.CallAfter(UpdateInstrumentGrid,G2frame,data)

    def MakeParameterWindow():
        'Displays the Instrument parameters in the dataWindow frame'

        def MakeLamSizer():
            if 'Lam1' in insVal:
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Azimuth: '),0,WACV)
                txt = '%7.2f'%(insVal['Azimuth'])
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'   Ka1/Ka2: '),0,WACV)
                txt = u'  %8.6f/%8.6f\xc5'%(insVal['Lam1'],insVal['Lam2'])
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  Source type: '),0,WACV)
                # PATCH?: for now at least, Source is not saved anywhere before here
                if 'Source' not in data: data['Source'] = ['CuKa','?']
                choice = list(G2elem.waves.keys())
                # ['TiKa','CrKa','FeKa','CoKa','CuKa','GaKa','MoKa','AgKa','InKa']
                # patch if Lam1/2 & only element is specified
                if 'Lam1' in data and data['Source'][1] not in choice:
                    indxs = [i for i,t in enumerate(choice) if
                             t.lower().startswith(data['Source'][1].lower())]
                    if len(indxs) == 1: data['Source'][1] = choice[indxs[0]]
                lamPick = wx.ComboBox(G2frame.dataWindow,value=data['Source'][1],choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0)
                subSizer.Add(waveSizer,0)
                mainSizer.Add(subSizer)
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,lblWdef('I(L2)/I(L1)',4,insDef['I(L2)/I(L1)'])),0,WACV)
                key = 'I(L2)/I(L1)'
                labelLst.append(key)
                elemKeysLst.append([key,1])
                dspLst.append([10,4])
                refFlgElem.append([key,2])
                ratVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,key,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                instSizer.Add(ratVal,0)
                instSizer.Add(RefineBox(key),0,WACV)
            else: # single wavelength
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Azimuth: '),0,WACV)
                txt = '%7.2f'%(insVal['Azimuth'])
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                instSizer.Add((5,5),0)
                key = 'Lam'
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' Lam (\xc5): (%10.6f)'%(insDef[key])),0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,key,nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                labelLst.append(u'Lam (\xc5)')
                elemKeysLst.append([key,1])
                dspLst.append([10,6])
                instSizer.Add(waveVal,0,WACV)
                refFlgElem.append([key,2])
                instSizer.Add(RefineBox(key),0,WACV)

        G2frame.dataWindow.ClearData()
        topSizer = G2frame.dataWindow.topBox
        parent = G2frame.dataWindow.topPanel
        topSizer.Add(wx.StaticText(parent,label=' Instrument settings'),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        mainSizer =  wx.BoxSizer(wx.VERTICAL)
        G2frame.dataWindow.SetSizer(mainSizer)
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        if insVal['Bank'] == None:      #patch
            insVal['Bank'] = 1
        text = ' Histogram Type: %s  Bank: %d'%(insVal['Type'],insVal['Bank'])
        if 'SEC' in insVal['Type']:
            text += ' (NB: Enter microscope voltage in keV to get wavelength)'
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,text))
        instSizer = wx.FlexGridSizer(0,3,5,5)
        labelLst[:],elemKeysLst[:],dspLst[:],refFlgElem[:] = [],[],[],[]
        if 'P' in insVal['Type']:                   #powder data
            [instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt),0,WACV) for txt in [' Name (default)',' Value','Refine?']]
            Reference = "Reference?"
            if insVal['Type'][2] in ['A','B','C']:               #constant wavelength
                labelLst.append('Azimuth angle')
                elemKeysLst.append(['Azimuth',1])
                dspLst.append([10,2])
                refFlgElem.append(None)
                MakeLamSizer()
                for item in ['Zero','Polariz.']:
                    if item in insDef:
                        labelLst.append(item)
                        elemKeysLst.append([item,1])
                        dspLst.append([10,4])
                        instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,lblWdef(item,4,insDef[item])),0,WACV)
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                        instSizer.Add(itemVal,0,WACV)
                        refFlgElem.append([item,2])
                        instSizer.Add(RefineBox(item),0,WACV)
                if 'C' in insVal['Type']:
                    itemList = ['U','V','W','X','Y','Z','SH/L']
                    Reference = """References:
    Thompson, P., Cox, D.E. & Hastings, J.B. (1987). J. Appl. Cryst. 20,79-83.
    Finger, L. W., Cox, D. E. & Jephcoat, A. P. (1994). J. Appl. Cryst. 27, 892-900.
                        """
                elif 'B' in insVal['Type']:
                    itemList = ['U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1']
                    if 'X' in data['Type']:
                        Reference = "Reference: Von Dreele, R.B., Clarke, S.M. & Walsh, J.P.S. (2021). J. Appl. Cryst., 54, 3-6."
                    else:
                        Reference = "Reference: R.B. Von Dreele (2024). J. Appl. Cryst. 57, 1588-1597."
                else: #'A'
                    itemList = ['U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1','SH/L']
                    Reference = """References:
    Thompson, P., Cox, D.E. & Hastings, J.B. (1987). J. Appl. Cryst. 20,79-83.
    Finger, L. W., Cox, D. E. & Jephcoat, A. P. (1994). J. Appl. Cryst. 27, 892-900.
    Von Dreele, R.B., Clarke, S.M. & Walsh, J.P.S. (2021). J. Appl. Cryst., 54, 3-6.
                        """
                for item in itemList:
                    nDig = (10,3)
                    if item == 'SH/L':
                        nDig = (10,5)
                    labelLst.append(item)
                    elemKeysLst.append([item,1])
                    dspLst.append(nDig)
                    refFlgElem.append([item,2])
                    instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,lblWdef(item,nDig[1],insDef[item])),0,WACV)
                    if item == 'SH/L':
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=nDig,typeHint=float,OnLeave=NewProfile,xmin=0.002)
                    else:
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=nDig,typeHint=float,OnLeave=NewProfile)
                    instSizer.Add(itemVal,0,WACV)
                    instSizer.Add(RefineBox(item),0,WACV)
            elif 'E' in insVal['Type']:
                key = '2-theta'
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' 2-theta (%10.6f):'%(insDef[key])),0,WACV)
                tthVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,key,nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                labelLst.append(u'2-theta')
                elemKeysLst.append([key,1])
                dspLst.append([10,3])
                instSizer.Add(tthVal,0,WACV)
                refFlgElem.append([key,2])
                instSizer.Add(RefineBox(key),0,WACV)
                for item in ['XE','YE','ZE','WE','A','B','C','X','Y','Z']:
                    nDig = (10,6,'g')
                    labelLst.append(item)
                    elemKeysLst.append([item,1])
                    dspLst.append(nDig)
                    refFlgElem.append([item,2])
                    instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,lblWdef(item,nDig[1],insDef[item])),0,WACV)
                    if item in ['XE','YE','ZE','WE']:
                        instSizer.Add(wx.StaticText(G2frame.dataWindow,label='%10.6g'%insVal[item]))
                        instSizer.Add((5,5),0)
                    else:
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=nDig,typeHint=float,OnLeave=NewProfile)
                        instSizer.Add(itemVal,0,WACV)
                        instSizer.Add(RefineBox(item),0,WACV)
            elif 'T' in insVal['Type']:                                   #time of flight (neutrons)
                Reference = """References:
    Von Dreele, R., Jorgensen, J. D. & Windsor, C. G. (1982) J. Appl. Cryst. 15, 581-589.
    Huq, A., Kirkham, M., Peterson, P.F., Hodges, J.P. Whitfield, P.S., Page, K., Hugle, T.,
        Iverson, E.B., Parizzia, A. & Rennich, G. (2019). J. Appl. Cryst. 52, 1189–1201.
                """
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Flight path: '),0,WACV)
                txt = '%8.3f'%(insVal['fltPath'])
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                labelLst.append('flight path')
                elemKeysLst.append(['fltPath',1])
                dspLst.append([10,2])
                refFlgElem.append(None)
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  2-theta: '),0,WACV)
                txt = '%7.2f'%(insVal['2-theta'])
                subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                labelLst.append('2-theta')
                elemKeysLst.append(['2-theta',1])
                dspLst.append([10,2])
                refFlgElem.append(None)
                if 'Pdabc' in Inst2:
                    Items = ['sig-0','sig-1','sig-2','sig-q','X','Y','Z']
                    subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  difC: '),0,WACV)
                    txt = '%8.2f'%(insVal['difC'])
                    subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,txt.strip()),0,WACV)
                    labelLst.append('difC')
                    elemKeysLst.append(['difC',1])
                    dspLst.append([10,2])
                    refFlgElem.append(None)
                    subSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  alpha, beta: fixed by table'),0,WACV)
                else:
                    Items = ['difC','difA','difB','Zero','alpha','beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','X','Y','Z']
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
                    instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,lblWdef(item,nDig[1],insDef[item])),0,WACV)
                    itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=nDig,typeHint=float,OnLeave=AfterChange)
                    instSizer.Add(itemVal,0,WACV)
                    labelLst.append(item)
                    elemKeysLst.append([item,1])
                    dspLst.append(nDig)
                    refFlgElem.append([item,2])
                    instSizer.Add(RefineBox(item),0,WACV)
            elif 'PKS' in insVal['Type']:   #peak positions only
                Reference = ''
                key = 'Lam'
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' Lam (\xc5): (%10.6f)'%(insDef[key])),0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,key,nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                labelLst.append(u'Lam (\xc5)')
                elemKeysLst.append([key,1])
                dspLst.append([10,6])
                instSizer.Add(waveVal,0,WACV)
                refFlgElem.append([key,2])
                instSizer.Add((5,5),0)
                for item in ['Zero',]:
                    if item in insDef:
                        labelLst.append(item)
                        elemKeysLst.append([item,1])
                        dspLst.append([10,4])
                        instSizer.Add(
                            wx.StaticText(G2frame.dataWindow,-1,lblWdef(item,4,insDef[item])),
                            0,WACV)
                        itemVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,item,nDig=(10,4),typeHint=float,OnLeave=AfterChange)
                        instSizer.Add(itemVal,0,WACV)
                        refFlgElem.append([item,2])
                        instSizer.Add((5,5),0)

        elif 'S' in insVal['Type']:                       #single crystal data
            Reference = ''
            if 'C' in insVal['Type']:               #constant wavelength
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),0,WACV)
                if 'EC' in insVal['Type']:
                    waveVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,'Lam',nDig=(10,6),typeHint=float,OnLeave=AfterChangeEC)
                else:
                    waveVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,'Lam',nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                instSizer.Add(waveVal,0,WACV)
                labelLst.append(u'Lam (\xc5)')
                waveSizer = wx.BoxSizer(wx.HORIZONTAL)
                waveSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  Source type: '),0,WACV)
                # PATCH?: for now at least, Source is not saved anywhere before here
                if 'Source' not in data: data['Source'] = ['CuKa','?']
                choice = ['synchrotron','TiKa','CrKa','FeKa','CoKa','CuKa','MoKa','AgKa','micro-ED']
                lamPick = wx.ComboBox(G2frame.dataWindow,value=data['Source'][1],choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                lamPick.Bind(wx.EVT_COMBOBOX, OnLamPick)
                waveSizer.Add(lamPick,0,WACV)
                instSizer.Add(waveSizer,0,WACV)
                elemKeysLst.append(['Lam',1])
                dspLst.append([10,6])
                refFlgElem.append(None)
            else:                                   #time of flight (neutrons)
                pass                                #for now
        elif insVal['Type'][0] in ['L','R',]:
            Reference = ''
            if 'C' in insVal['Type']:
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' Lam (\xc5): (%10.6f)'%(insDef['Lam'])),0,WACV)
                waveVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,insVal,'Lam',nDig=(10,6),typeHint=float,OnLeave=AfterChange)
                instSizer.Add(waveVal,0,WACV)
                labelLst.append(u'Lam (\xc5)')
                elemKeysLst.append(['Lam',1])
                dspLst.append([10,6])
                refFlgElem.append(None)
                instSizer.Add(wx.StaticText(G2frame.dataWindow,-1,'  Azimuth: %7.2f'%(insVal['Azimuth'])),0,WACV)
                labelLst.append('Azimuth angle')
                elemKeysLst.append(['Azimuth',1])
                dspLst.append([10,2])
                refFlgElem.append(None)
            else:                                   #time of flight (neutrons)
                pass                                #for now

        mainSizer.Add(instSizer,0)
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=Reference))
        G2frame.dataWindow.SetDataSize()
        # end of MakeParameterWindow

    plotYsel = {}   # selected Y items
    def MakeMultiParameterWindow(selected=None):
        '''Displays the Instrument parameters for multiple histograms
        in the dataWindow panel
        '''
        plotIndex = {'plotX':0} # index for param name => plotTbl index
                                # plotX is selected X axis
        plotTbl = []   # table of values for each param
        plotLabel = []   # table of values for each param
        def onSelectHists(event):
            'select histograms to show'
            dlg = G2G.G2MultiChoiceDialog(G2frame,
                        'Select histograms to show of type '+data['Type'][1],
                        'Select histograms',hlist)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    selected = dlg.GetSelections()
                else:
                    return
            finally:
                dlg.Destroy()
            wx.CallAfter(MakeMultiParameterWindow,selected)
            return
        def onSelectX(event):
            'respond to change in plotting x axis; save and plot (if y selected)'
            plotIndex['plotX'] = event.GetEventObject().rbindex
            onPrmPlot(event)

        def onPrmPlot(event):
            '''Callback after a change to X or Y plot contents
            plots multiple instrument param values vs selected X value.
            If no Y values are selected, any previous plot is deleted.
            '''
            xvals = plotTbl[plotIndex['plotX']]
            xlbl  = plotLabel[plotIndex['plotX']]
            XY = []
            keys = ''
            for k in plotYsel:
                if k == 'plotX': continue
                if not plotYsel[k]: continue
                yvals = plotTbl[plotIndex[k]]
                if keys: keys += ', '
                keys += plotLabel[plotIndex[k]]
                XY.append((xvals,yvals))
            if not XY:
                G2frame.G2plotNB.Delete('Parameter values')
                return
            G2plt.PlotXY(G2frame,XY,labelX=xlbl,labelY=keys,Title='Parameter values',newPlot=True)

        # gather histograms matching the currently selected histogram
        histoList,histIdList = G2frame.GetHistogramNamesID(['PWDR',])
        hlist = GetHistsLikeSelected(G2frame)
        if selected is None and len(hlist) > 10:  # on initial call this is none
            onSelectHists(None) # lots of histograms, give user a chance to choose
            return
        elif selected is None:  # select all, not so many
            selected = range(len(hlist))
            lbl = 'Instrument parameters for all matching histograms'
        else:
            lbl = 'Instrument parameters for selected histograms'
        wx.BeginBusyCursor()
        h = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histnames = [h]
        histdict = {h:data}
        histnum = {h:G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[0]['hId']}
        for i in selected:
            h = hlist[i]
            if h not in histoList: # unexpected
                print('hist from GetHistsLikeSelected not in GetHistogramNamesID',h)
                continue
            hid = histIdList[h]
            inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,hid, 'Instrument Parameters'))
            histnames.append(h)
            histdict[h] = inst
            histnum[h] = G2frame.GPXtree.GetItemPyData(hid)[0]['hId']

        # start posting info into window
        G2frame.dataWindow.ClearData()
        topSizer = G2frame.dataWindow.topBox
        parent = G2frame.dataWindow.topPanel
        topSizer.Add(wx.StaticText(parent,wx.ID_ANY,lbl))
        if hlist:
            btn = wx.Button(parent, wx.ID_ANY,'Select\nHistograms')
            topSizer.Add(btn,0,wx.LEFT|wx.RIGHT,15)
            btn.Bind(wx.EVT_BUTTON,onSelectHists)
        topSizer.Add((20,-1))
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        mainSizer =  wx.BoxSizer(wx.VERTICAL)
        G2frame.dataWindow.SetSizer(mainSizer)

        # create table w/headers
        sdlg = G2frame.dataWindow
        fgs = wx.FlexGridSizer(0,len(histnames)+3,0,0)
        fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,'plot\nas X'),0,wx.LEFT|wx.RIGHT,1)
        fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,'plot\nas Y'),0,wx.LEFT|wx.RIGHT,1)
        fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,'Histogram  '),0,WACV|wx.LEFT,14)
        for i,h in enumerate(histnames):
            if len(h[:5].strip()) > 20:
                fgs.Add(G2G.ScrolledStaticText(sdlg,label=h[5:],dots=False,lbllen=20),
                    0,WACV|wx.LEFT|wx.RIGHT,5)
            else:
                fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,h[5:]),0,WACV|wx.LEFT|wx.RIGHT,5)
        firstRadio = wx.RB_GROUP
        # put non-editable values at top of table (plot as x but not y)
        keylist = ['num']
        lbllist = ['#']
        if 'T' in data['Type'][1]:
            keylist += ['2-theta']
            lbllist += ['2-theta']
        for key,lbl in zip(keylist,lbllist):
            rb = wx.RadioButton(sdlg,wx.ID_ANY,'',style=firstRadio)
            rb.rbindex = len(plotTbl)
            rb.Bind(wx.EVT_RADIOBUTTON,onSelectX)
            plotLabel.append(lbl)
            fgs.Add(rb,0,wx.ALIGN_CENTER|WACV)
            if firstRadio:
                rb.SetValue(True)
                firstRadio = 0
            fgs.Add((-1,-1)) # skip y checkbutton
            fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,lbl),0,WACV|wx.LEFT,14)
            plotvals = []
            for h in histnames:
                if key == 'num':
                    val = histnum[h]
                else:
                    val = histdict[h][key][1]
                fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,str(val)),
                    0,wx.ALIGN_CENTER|WACV,0)
                plotvals.append(val)
            plotTbl.append(plotvals)

        # determine what items will be shown based on histogram type
        Items = []
        if data['Type'][1][2] in ['A','B','C']:               #constant wavelength
            if 'Lam1' in data:
                Items = ['Lam1','Lam2','I(L2)/I(L1)']
            else:
                Items = ['Lam','Zero','Polariz.']
            if 'C' in data['Type'][1]:
                Items += ['Azimuth','U','V','W','X','Y','Z','SH/L']
            elif 'B' in data['Type'][1]:
                Items += ['Azimuth','U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1']
            else: #'A'
                Items += ['Azimuth','U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1','SH/L']
        elif 'E' in data['Type'][1]:
            Items = ['2-theta','XE','YE','ZE','WE','A','B','C','X','Y','Z']
        elif 'T' in data['Type'][1]:            # TOF
            Items = ['difC','difA','difB','Zero','alpha','beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q','X','Y','Z']
        # display the items in the table
        for k in Items:
            plotYsel[k] = plotYsel.get(k,False)
            #if not l: l = k
            l = k
            rb = wx.RadioButton(sdlg,wx.ID_ANY,'',style=firstRadio)
            rb.rbindex = len(plotTbl)
            rb.Bind(wx.EVT_RADIOBUTTON,onSelectX)
            plotLabel.append(l)
            fgs.Add(rb,0,wx.ALIGN_CENTER|WACV)
            if firstRadio:
                rb.SetValue(True)
                firstRadio = 0
            fgs.Add(G2G.G2CheckBox(sdlg,'',plotYsel,k,OnChange=onPrmPlot),0,wx.ALIGN_CENTER|WACV)
            plotIndex[k] = rb.rbindex
            fgs.Add(wx.StaticText(sdlg,wx.ID_ANY,l),0,WACV|wx.LEFT,14)
            plotvals = []
            for h in histnames:
                miniSizer = wx.BoxSizer(wx.HORIZONTAL)
                itemVal = G2G.ValidatedTxtCtrl(sdlg,histdict[h][k],1,nDig=(10,4),typeHint=float)
                plotvals.append(histdict[h][k][1])
                miniSizer.Add(itemVal)
                miniSizer.Add((2,-1))
                miniSizer.Add(G2G.G2CheckBox(sdlg,'',histdict[h][k],2),0,WACV|wx.RIGHT,15)
                fgs.Add(miniSizer,0,wx.ALIGN_CENTER)
            plotTbl.append(plotvals)

        mainSizer.Add(fgs)
        G2frame.dataWindow.SetDataSize()
        wx.EndBusyCursor()
        # end of MakeMultiParameterWindow

    #### beginning of UpdateInstrumentGrid code
    #patch: make sure all parameter items are lists
    patched = 0
    for key in data:
        if type(data[key]) is tuple:
            data[key] = list(data[key])
            patched += 1
    if patched: print (patched,' instrument parameters changed from tuples')
    if 'E' not in data['Type'][0] and 'Z' not in data:
        data['Z'] = [0.0,0.0,False]
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
    elif 'R' in data['Type'][0]:                               #Reflectometry data
        insVal = dict(zip(instkeys,[data[key][1] for key in instkeys]))
        insDef = dict(zip(instkeys,[data[key][0] for key in instkeys]))
        insRef = {}
    RefObj = {}
    Inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))[1]
    G2gd.SetDataMenuBar(G2frame)
    #patch
    if 'P' in insVal['Type']:                   #powder data
        if 'C' in insVal['Type']:               #constant wavelength
            if 'Azimuth' not in insVal:
                insVal['Azimuth'] = 0.0
                insDef['Azimuth'] = 0.0
                insRef['Azimuth'] = False
        elif 'E' in insVal['Type']:
            if 'WE' not in insVal:
                    insVal['WE'] = 0.0
                    insDef['WE'] = 0.0
                    insRef['WE'] = False
            if 'X' not in insVal:
                for item in ['X','Y','Z']:
                    insVal[item] = 0.0
                    insDef[item] = 0.0
                    insRef[item] = False
    #end of patch
    if 'P' in insVal['Type']:                   #powder data menu commands
        G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.InstMenu)
        G2frame.GetStatusBar().SetStatusText('NB: Azimuth is used for polarization only',1)
        G2frame.Bind(wx.EVT_MENU,OnCalibrate,id=G2G.wxID_INSTCALIB)
        G2frame.Bind(wx.EVT_MENU,OnLoad,id=G2G.wxID_INSTLOAD)
        G2frame.Bind(wx.EVT_MENU,OnSave,id=G2G.wxID_INSTSAVE)
        G2frame.Bind(wx.EVT_MENU,OnSaveAll,id=G2G.wxID_INSTSAVEALL)
        G2frame.Bind(wx.EVT_MENU,OnReset,id=G2G.wxID_INSTPRMRESET)
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2G.wxID_INSTCOPY)
        G2frame.Bind(wx.EVT_MENU,OnInstFlagCopy,id=G2G.wxID_INSTFLAGCOPY)
        G2frame.Bind(wx.EVT_MENU,OnCopy1Val,id=G2G.wxID_INST1VAL)
        G2frame.Bind(wx.EVT_MENU,OnInstMult,id=G2G.wxID_INSTSHOWMULT)
        menuitem = G2frame.dataWindow.InstMenu.FindItemById(G2G.wxID_INSTSHOWMULT)
        if menuitem.IsChecked():
            MakeMultiParameterWindow()
        else:
            MakeParameterWindow()
    elif 'L' in insVal['Type'] or 'R' in insVal['Type']:                   #SASD & REFD data menu commands
        MakeParameterWindow()
        G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SASDInstMenu)
        G2frame.Bind(wx.EVT_MENU,OnInstCopy,id=G2G.wxID_SASDINSTCOPY)
    G2frame.dataWindow.SetDataSize()

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
        '''
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II sample parameters file', pth, '',
            'sample parameter files (*.samprm)|*.samprm',wx.FD_OPEN)
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
                G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Sample Parameters'),data)
                UpdateSampleGrid(G2frame,data)
        finally:
            dlg.Destroy()

    def OnAllSampleLoad(event):
        filename = ''
        pth = G2G.GetImportPath(G2frame)
        if not pth: pth = '.'
        dlg = wx.FileDialog(G2frame, 'Choose multihistogram metadata text file', pth, '',
            'metadata file (*.*)|*.*',wx.FD_OPEN)
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
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Controls'))
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
        histList = [G2frame.GPXtree.GetItemText(G2frame.PatternId),]
        histList += GetHistsLikeSelected(G2frame)
        colIds = {}
        for i,name in enumerate(colNames):
            if name != ' ':
                colIds[name] = i
        for ih,hist in enumerate(histList):
            name = hist.split()[1]  #this is file name
            newItems = {}
            for item in colIds:
                key = freeNames.get(item,item)
                try:
                    newItems[key] = float(dataDict[name][colIds[item]])
                except KeyError:
                    try:
                        newItems[key] = float(dataDict['{:}'.format(ih+1)])
                    except KeyError:
                        break
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
            sampleData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
            sampleData.update(newItems)
        UpdateSampleGrid(G2frame,data)

    def OnSetScale(event):
        if histName[:4] in ['REFD','PWDR']:
            Scale = data['Scale'][0]
            dlg = wx.MessageDialog(G2frame,'Rescale data by %.2f?'%(Scale),'Rescale data',wx.OK|wx.CANCEL)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histName)
                    y,w = G2frame.GPXtree.GetItemPyData(pId)[1][1:3]
                    y *= Scale
                    w /= Scale**2
                    data['Scale'][0] = 1.0
            finally:
                dlg.Destroy()
            G2pwpl.PlotPatterns(G2frame,plotType=histName[:4],newPlot=True)
            UpdateSampleGrid(G2frame,data)
            return
        #SASD rescaliing
        histList = []
        item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
        while item:
            name = G2frame.GPXtree.GetItemText(item)
            if 'SASD' in name and name != histName:
                histList.append(name)
            item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
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
        Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))
        Profile = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1]
        Data = [Profile,Limits,data]
        refId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,refHist)
        refSample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,refId, 'Sample Parameters'))
        refLimits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,refId, 'Limits'))
        refProfile = G2frame.GPXtree.GetItemPyData(refId)[1]
        refData = [refProfile,refLimits,refSample]
        G2sasd.SetScale(Data,refData)
        G2pwpl.PlotPatterns(G2frame,plotType='SASD',newPlot=True)
        UpdateSampleGrid(G2frame,data)

    def OnRescaleAll(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        x0,y0,w0 = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1][:3]
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        od = {'label_1':'Scaling range min','value_1':0.0,'label_2':'Scaling range max','value_2':10.}
        dlg = G2G.G2MultiChoiceDialog(G2frame,
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
                        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                        xi,yi,wi = G2frame.GPXtree.GetItemPyData(Id)[1][:3]
                        sumi = np.sum(yi[iBeg:iFin])
                        if sumi:
                            Scale = sum0/sumi
                            yi *= Scale
                            wi /= Scale**2
        finally:
            dlg.Destroy()
        G2pwpl.PlotPatterns(G2frame,plotType=histName[:4],newPlot=True)

    def OnSampleCopy(event):
        histType,copyNames = SetCopyNames(histName,data['Type'],
            addNames = ['Omega','Chi','Phi','Gonio. radius','InstrName'])
        copyDict = {}
        for parm in copyNames:
            copyDict[parm] = data[parm]
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy sample params from\n'+str(hst[5:])+' to...',
            'Copy sample parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result:
                    item = histList[i]
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                    sampleData.update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnSampleCopySelected(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        Controls = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        # Assemble a list of item labels
        TextTable = {key:label for key,label,dig in
            SetupSampleLabels(hst,data.get('Type'),Inst['Type'][0])}
        # get flexible labels
        TextTable.update({key:Controls[key] for key in Controls if key.startswith('FreePrm')})
        # add a few extra
        TextTable.update({'Type':'Diffractometer type','InstrName':'Instrument Name',})
        if 'SASD' in histList[0] or 'REFD' in histList[0]:
            TextTable.update({'Materials':'Materials',})
        # Assemble a list of dict entries that would be labeled in the Sample
        # params data window (drop ranId and items not used).
        keyList = [i for i in data.keys() if i in TextTable]
        keyText = [TextTable[i] for i in keyList]
        # sort both lists together, ordered by keyText
        keyText, keyList = zip(*sorted(list(zip(keyText,keyList)))) # sort lists
        selectedKeys = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select which sample parameters\nto copy',
            'Select sample parameters', keyText)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                selectedKeys = [keyList[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if not selectedKeys: return # nothing to copy
        copyDict = {}
        for parm in selectedKeys:
            copyDict[parm] = copy.deepcopy(data[parm])
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy sample params from\n'+str(hst[5:])+' to...',
            'Copy sample parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result:
                    item = histList[i]
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                    sampleData.update(copyDict)
        finally:
            dlg.Destroy()
        G2pwpl.PlotPatterns(G2frame,plotType=hst[:4],newPlot=False)

    def OnSampleFlagCopy(event):
        histType,copyNames = SetCopyNames(histName,data['Type'])
        flagDict = {}
        for parm in copyNames:
            flagDict[parm] = data[parm][1]
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy sample ref. flags from\n'+str(hst[5:])+' to...',
            'Copy sample flags', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result:
                    item = histList[i]
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    sampleData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                    for name in copyNames:
                        if name not in sampleData:
                            sampleData[name] = [0.0,False]
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
            G2pwpl.PlotPatterns(G2frame,plotType='SASD',newPlot=True)
        elif tc.key == 'Thick':
            wx.CallAfter(UpdateSampleGrid,G2frame,data)

    def OnMaterial(event):
        Obj = event.GetEventObject()
        Id = Info[Obj.GetId()]
        data['Materials'][Id]['Name'] = Obj.GetValue()
        wx.CallAfter(UpdateSampleGrid,G2frame,data)

    def OnVolFrac(invalid,value,tc):
        Id = Info[tc.GetId()]
        data['Materials'][not Id]['VolFrac'] = 1.-value
        wx.CallAfter(UpdateSampleGrid,G2frame,data)

    def OnCopy1Val(event):
        'Select one value to copy to many histograms and optionally allow values to be edited in a table'
        G2G.SelectEdit1Var(G2frame,data,labelLst,elemKeysLst,dspLst,refFlgElem)
        wx.CallAfter(UpdateSampleGrid,G2frame,data)

    def SearchAllComments(value,tc,*args,**kwargs):
        '''Called when the label for a FreePrm is changed: the comments for all PWDR
        histograms are searched for a "label=value" pair that matches the label (case
        is ignored) and the values are then set to this value, if it can be converted
        to a float.
        '''
        Id, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
        while Id:
            name = G2frame.GPXtree.GetItemText(Id)
            if 'PWDR' in name:
                Comments = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Comments'))
                Sample =   G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Sample Parameters'))
                for i,item in enumerate(Comments):
                    itemSp = item.split('=')
                    if value.lower() == itemSp[0].lower():
                        try:
                            Sample[tc.key] = float(itemSp[1])
                        except:
                            print('"{}" has an invalid value in Comments from {}'
                                  .format(item.strip(),name))
            Id, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        wx.CallLater(100,UpdateSampleGrid,G2frame,data)

    # start of UpdateSampleGrid
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu) # needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SampleMenu)
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
            G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    histName = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    #G2frame.SetLabel(G2frame.GetLabel().split('||')[0]+' || '+'Sample Parameters')
    G2frame.Bind(wx.EVT_MENU, OnSetScale, id=G2G.wxID_SETSCALE)
    G2frame.Bind(wx.EVT_MENU, OnSampleCopy, id=G2G.wxID_SAMPLECOPY)
    G2frame.Bind(wx.EVT_MENU, OnSampleCopySelected, id=G2G.wxID_SAMPLECOPYSOME)
    G2frame.Bind(wx.EVT_MENU, OnSampleFlagCopy, id=G2G.wxID_SAMPLEFLAGCOPY)
    G2frame.Bind(wx.EVT_MENU, OnSampleSave, id=G2G.wxID_SAMPLESAVE)
    G2frame.Bind(wx.EVT_MENU, OnSampleLoad, id=G2G.wxID_SAMPLELOAD)
    G2frame.Bind(wx.EVT_MENU, OnCopy1Val, id=G2G.wxID_SAMPLE1VAL)
    G2frame.Bind(wx.EVT_MENU, OnAllSampleLoad, id=G2G.wxID_ALLSAMPLELOAD)
    G2frame.Bind(wx.EVT_MENU, OnRescaleAll, id=G2G.wxID_RESCALEALL)
    if histName[:4] in ['SASD','REFD','PWDR']:
        G2frame.dataWindow.SetScale.Enable(True)
    Controls = G2frame.GPXtree.GetItemPyData(
        G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
#patch
    if 'ranId' not in data:
        data['ranId'] = ran.randint(0,sys.maxsize)
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
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label=' Sample and Experimental Parameters'),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    nameSizer = wx.BoxSizer(wx.HORIZONTAL)
    nameSizer.Add(wx.StaticText(G2frame.dataWindow,wx.ID_ANY,' Instrument Name '),0,WACV)
    nameSizer.Add((-1,-1),1,WACV)
    instNameVal = wx.TextCtrl(G2frame.dataWindow,wx.ID_ANY,data['InstrName'],
        size=(200,-1),style=wx.TE_PROCESS_ENTER)
    nameSizer.Add(instNameVal)
    instNameVal.Bind(wx.EVT_CHAR,OnNameVal)
    mainSizer.Add(nameSizer,0)
    mainSizer.Add((5,5),0)
    labelLst.append('Instrument Name')
    elemKeysLst.append(['InstrName'])
    dspLst.append(None)
    refFlgElem.append(None)

    if 'PWDR' in histName:
        nameSizer = wx.BoxSizer(wx.HORIZONTAL)
        nameSizer.Add(wx.StaticText(G2frame.dataWindow,wx.ID_ANY,' Diffractometer type: '),0,WACV)
        if 'T' in Inst['Type'][0] or 'E' in Inst['Type'][0]:
            choices = ['Debye-Scherrer',]
        else: #'[A','B','C']
            choices = ['Debye-Scherrer','Bragg-Brentano',]
        histoType = G2G.G2ChoiceButton(G2frame.dataWindow,choices,strLoc=data,strKey='Type',
            onChoice=OnHistoChange)
        nameSizer.Add(histoType)
        mainSizer.Add(nameSizer,0)
        mainSizer.Add((5,5),0)

    parmSizer = wx.FlexGridSizer(0,2,5,0)
    for key,lbl,nDig in parms:
        if 'E' in Inst['Type'][0] and key in ['DisplaceX','DisplaceY','Absorption']:
            continue
        labelLst.append(lbl.strip().strip(':').strip())
        dspLst.append(nDig)
        if 'list' in str(type(data[key])):
            parmRef = G2G.G2CheckBox(G2frame.dataWindow,' '+lbl,data[key],1)
            parmSizer.Add(parmRef,0,wx.EXPAND)
            parmVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[key],0,
                nDig=nDig,typeHint=float,OnLeave=AfterChange)
            elemKeysLst.append([key,0])
            refFlgElem.append([key,1])
        else:
            parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' '+lbl),
                0,wx.EXPAND)
            parmVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,key,
                typeHint=float,OnLeave=AfterChange)
            elemKeysLst.append([key])
            refFlgElem.append(None)
        parmSizer.Add(parmVal,0,WACV)
    Info = {}

    for key in ('FreePrm1','FreePrm2','FreePrm3'):
        parmVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,Controls,key,typeHint=str,
            notBlank=False,OnLeave=SearchAllComments)
        parmSizer.Add(parmVal,1,wx.EXPAND)
        parmVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,key,typeHint=float)
        parmSizer.Add(parmVal,0,WACV)
        labelLst.append(Controls[key])
        dspLst.append(None)
        elemKeysLst.append([key])
        refFlgElem.append(None)

    mainSizer.Add(parmSizer,0)
    mainSizer.Add((0,5),0)
    if histName[:4] in ['SASD',]:
        rho = [0.,0.]
        anomrho = [0.,0.]
        mu = 0.
        subSizer = wx.FlexGridSizer(0,4,5,5)
        Substances = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Substances'))
        for Id,item in enumerate(data['Materials']):
            subSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Material: '),0,WACV)
            matsel = wx.ComboBox(G2frame.dataWindow,value=item['Name'],choices=list(Substances['Substances'].keys()),
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Info[matsel.GetId()] = Id
            matsel.Bind(wx.EVT_COMBOBOX,OnMaterial)
            subSizer.Add(matsel,0,WACV)
            subSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Volume fraction: '),0,WACV)
            volfrac = G2G.ValidatedTxtCtrl(G2frame.dataWindow,item,'VolFrac',
                xmin=0.,xmax=1.,nDig=(10,3),typeHint=float,OnLeave=OnVolFrac)
            Info[volfrac.GetId()] = Id
            subSizer.Add(volfrac,0,WACV)
            try:
                material = Substances['Substances'][item['Name']]
            except KeyError:
                print('ERROR - missing substance: '+item['Name'])
                material = Substances['Substances']['vacuum']
            mu += item['VolFrac']*material.get('XAbsorption',0.)
            rho[Id] = material['Scatt density']
            anomrho[Id] = material.get('XAnom density',0.)
        data['Contrast'] = [(rho[1]-rho[0])**2,(anomrho[1]-anomrho[0])**2]
        mainSizer.Add(subSizer,0)
        conSizer = wx.BoxSizer(wx.HORIZONTAL)
        conSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Contrast: %10.2f '%(data['Contrast'][0])),0,WACV)
        conSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Anom. Contrast: %10.2f '%(data['Contrast'][1])),0,WACV)
        mut =  mu*data['Thick']
        conSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Transmission (calc): %10.3f  '%(np.exp(-mut))),0,WACV)
        mainSizer.Add(conSizer,0)
    G2frame.dataWindow.SetDataSize()

################################################################################
#####  Indexing Peaks
################################################################################

def UpdateIndexPeaksGrid(G2frame, data):
    '''respond to selection of PWDR Index Peak List data
    tree item.
    '''
    IndexId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Index Peak List')
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    limitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
    Limits = G2frame.GPXtree.GetItemPyData(limitId)

    def RefreshIndexPeaksGrid(event):
        r,c =  event.GetRow(),event.GetCol()
        peaks = G2frame.IndexPeaksTable.GetData()
        if c == 2:
            peaks[r][c] = not peaks[r][c]
            G2frame.IndexPeaksTable.SetData(peaks)
            G2frame.indxPeaks.ForceRefresh()
            if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
                G2plt.PlotPowderLines(G2frame)
            else:
                G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def onCellListDClick(event):
        '''Called after a double-click on a row/column label'''
        r,c =  event.GetRow(),event.GetCol()
        if r >= 0 and c < 0:     #row label: select it and replot!
            G2frame.indxPeaks.ClearSelection()
            G2frame.indxPeaks.SelectRow(r,True)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR')
            wx.CallAfter(G2frame.indxPeaks.ForceRefresh)

    def OnReload(event):
        peaks = []
        sigs = []
        Peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
        for ip,peak in enumerate(Peaks['peaks']):
            dsp = G2lat.Pos2dsp(Inst,peak[0])
            peaks.append([peak[0],peak[2],True,False,0,0,0,dsp,0.0])    #SS?
            try:
                sig = Peaks['sigDict']['pos'+str(ip)]
            except KeyError:
                sig = 0.
            sigs.append(sig)
        data = [peaks,sigs]
        G2frame.GPXtree.SetItemPyData(IndexId,data)
        UpdateIndexPeaksGrid(G2frame,data)

    def OnSave(event):
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose Index peaks csv file', pth, '',
            'indexing peaks file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.csv'
                File = open(filename,'w')
                names = 'h,k,l,position,intensity,d-Obs,d-calc\n'
                File.write(names)
                fmt = '%d,%d,%d,%.4f,%.1f,%.5f,%.5f\n'
                for refl in data[0]:
                    if refl[3]:
                        File.write(fmt%(refl[4],refl[5],refl[6],refl[0],refl[1],refl[7],refl[8]))
                File.close()
        finally:
            dlg.Destroy()

    def OnExportPreDICT(event):
        'Place 2theta positions from Index Peak List into clipboard for cut-&-paste'
        if wx.TheClipboard.Open():
            txt = ''
            c = 0
            for refl in data[0]:
                if refl[2]:
                    c += 1
                    txt += '{}\n'.format(refl[0])
            wx.TheClipboard.SetData(wx.TextDataObject(txt))
            wx.TheClipboard.Close()
        else:
            G2frame.ErrorDialog('Clipboard locked','Sorry, unable to access the clipboard, try again later. You might need to restart GSAS-II or reboot')
            return
        G2G.G2MessageBox(G2frame,
                '{} reflection positions placed in clipboard. '.format(c)+
                'In PreDICT open the DICVOL input. Update the number of '
                'reflections (& wavelength if needed) in GUI. Then'
                '  use paste to replace the reflections '
                '(starting at line 6...) in the DICVOL input.',
                'DICVOL input generated')

    def KeyEditPickGrid(event):
        colList = G2frame.indxPeaks.GetSelectedCols()
        data = G2frame.GPXtree.GetItemPyData(IndexId)
        if event.GetKeyCode() == wx.WXK_RETURN:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_CONTROL:
            event.Skip(True)
        elif event.GetKeyCode() == wx.WXK_SHIFT:
            event.Skip(True)
        elif colList:
            G2frame.indxPeaks.ClearSelection()
            key = event.GetKeyCode()
            for col in colList:
                if G2frame.IndexPeaksTable.GetColLabelValue(col) in ['use',]:
                    if key == 89: #'Y'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col]=True
                    elif key == 78:  #'N'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col]=False
                    elif key == 83: # 'S'
                        for row in range(G2frame.IndexPeaksTable.GetNumberRows()): data[0][row][col] = not data[0][row][col]

    def onRefineCell(event):
        RefineCell(G2frame)
        UpdateIndexPeaksGrid(G2frame,data)

    # start of UpdateIndexPeaksGrid
    controls = None
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.IndexMenu) # needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu) # needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.IndPeaksMenu)
    if 'PWD' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
        G2frame.Bind(wx.EVT_MENU, OnReload, id=G2G.wxID_INDXRELOAD)
        G2frame.Bind(wx.EVT_MENU, OnSave, id=G2G.wxID_INDEXSAVE)
        G2frame.Bind(wx.EVT_MENU, OnExportPreDICT, id=G2G.wxID_INDEXEXPORTDICVOL)
        G2frame.Bind(wx.EVT_MENU, onRefineCell, id=G2G.wxID_REFINECELL2)
    G2frame.dataWindow.IndexPeaks.Enable(False)
    G2frame.IndexPeaksTable = []
    if len(data[0]):
        G2frame.dataWindow.IndexPeaks.Enable(True)
        Unit = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List'))
        if Unit:
            if len(Unit) == 4:  #patch
                Unit.append({})
            if len(Unit) == 5:  #patch
                Unit.append({})
            controls,bravais,cellist,dmin,ssopt,magcells = Unit
            if 'T' in Inst['Type'][0]:   #TOF - use other limit!
                dmin = G2lat.Pos2dsp(Inst,Limits[1][0])
            else:
                dmin = G2lat.Pos2dsp(Inst,Limits[1][1])
            G2frame.HKL = np.array([])
            G2frame.Extinct = []
            if ssopt.get('Use',False):
                cell = controls[6:12]
                A = G2lat.cell2A(cell)
                spc = controls[13]
                SGData = G2spc.SpcGroup(spc)[1]
                SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
                Vec = ssopt['ModVec']
                maxH = ssopt['maxH']
                G2frame.HKL = np.array(G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,Vec,maxH,A))
                data[0] = G2indx.IndexSSPeaks(data[0],G2frame.HKL)[1]
            else:        #select cell from table - no SS
                for i,cell in enumerate(cellist):
                    if cell[-2]:
                        ibrav = cell[2]
                        A = G2lat.cell2A(cell[3:9])
                        G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A,ifList=True)
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
    G2frame.GPXtree.SetItemPyData(IndexId,data)
    G2frame.IndexPeaksTable = G2G.Table(data[0],rowLabels=rowLabels,colLabels=colLabels,types=Types)
    G2frame.dataWindow.currentGrids = []
    G2frame.indxPeaks = G2G.GSGrid(parent=G2frame.dataWindow)
    G2frame.indxPeaks.SetRowLabelSize(45)
    G2frame.indxPeaks.SetTable(G2frame.IndexPeaksTable, True)
    XY = []
    Sigs = []
    for r in range(G2frame.indxPeaks.GetNumberRows()):
        for c in range(G2frame.indxPeaks.GetNumberCols()):
            if c == 2:
                G2frame.indxPeaks.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.indxPeaks.SetReadOnly(r,c,isReadOnly=True)
        if data[0][r][2] and data[0][r][3]:
            XY.append([data[0][r][-1],data[0][r][0]])
            try:
                sig = data[1][r]
            except IndexError:
                sig = 0.
            Sigs.append(sig)
    G2frame.indxPeaks.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, RefreshIndexPeaksGrid)
    G2frame.indxPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, onCellListDClick)
    G2frame.indxPeaks.Bind(wx.EVT_KEY_DOWN, KeyEditPickGrid)
    G2frame.indxPeaks.AutoSizeColumns(False)
    if len(XY):
        XY = np.array(XY)
        G2plt.PlotCalib(G2frame,Inst,XY,Sigs,newPlot=True)
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Index peaks list'),0,WACV)
    # add help button to bring up help web page - at right side of window
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add(G2frame.indxPeaks,1,wx.EXPAND,1)
    botSizer = G2frame.dataWindow.bottomBox
    parent = G2frame.dataWindow.bottomPanel
    if controls:
        ibrav = SetLattice(controls)
        for cellGUI in cellGUIlist:
            if ibrav in cellGUI[0]:
                useGUI = cellGUI
        botSizer.Add(wx.StaticText(parent,label='Cell: ',style=wx.ALIGN_RIGHT),0,WACV)
        cellSizer = wx.FlexGridSizer(0,min(6,useGUI[1]),3,3)
        botSizer.Add(cellSizer,0,WACV)
        for txt,fmt,ifEdit,Id in zip(*useGUI[2]):
            if 'Vol' in txt:
                val = fmt % controls[12]
                botSizer.Add(wx.StaticText(parent,label=txt,style=wx.ALIGN_RIGHT),0,WACV)
                volVal = wx.TextCtrl(parent,value=val,style=wx.TE_READONLY,size=(65,-1))
                botSizer.Add(volVal,0,WACV)
            else:
                val = f'%.{fmt[1]}f' % controls[6+Id]
                cellSizer.Add(wx.StaticText(parent,label=txt,style=wx.ALIGN_RIGHT),0,wx.ALIGN_RIGHT|WACV)
                volVal = wx.TextCtrl(parent,value=val,style=wx.TE_READONLY,size=(65,-1))
                cellSizer.Add(volVal,0,WACV)
            volVal.SetBackgroundColour(VERY_LIGHT_GREY)
        if not ssopt.get('Use',False):        #zero for super lattice doesn't work!
            botSizer.Add((15,-1))
            vcSizer =  wx.BoxSizer(wx.VERTICAL)
            vcSizer.Add(wx.StaticText(parent,label="Zero offset",
                                         style=wx.ALIGN_CENTER),0,wx.EXPAND)
            hcSizer = wx.BoxSizer(wx.HORIZONTAL)
            zero = G2G.ValidatedTxtCtrl(parent,controls,1,nDig=(10,4),typeHint=float,
                    xmin=-5.,xmax=5.,size=(50,-1))
            hcSizer.Add(zero,0,WACV)
            zeroVar = G2G.G2CheckBox(parent,'Ref?',controls,0)
            hcSizer.Add(zeroVar,0,WACV|wx.LEFT,3)
            vcSizer.Add(hcSizer)
            botSizer.Add(vcSizer,0,WACV)
# TODO: get SUs from cell refinement
# TODO: implement Enable on G2frame.dataWindow.RefineCell2 to match G2frame.dataWindow.RefineCell
    G2frame.dataWindow.SetDataSize()

################################################################################
#####  Unit cells
################################################################################
def UpdateUnitCellsGrid(G2frame, data, callSeaResSelected=False,New=False,showUse=False):
    '''respond to selection of PWDR Unit Cells data tree item.

    :param wx.Frame G2frame: Main GSAS-II window
    :param dict data: contents of "Unit Cells List" data tree item
    :param bool callSeaResSelected: when True, selects first entry in
      UnitCellsTable search results table
    :param bool New:
    :param bool showUse: when showUse is False (default) the Show flag
      is cleared in all search tables. When True, and there is a True value
      for Show, the flag is set for that in the grid and the row is scrolled
      into view. This is currently implemented for indexing (Cell Search
      Results) only.
    '''
    global KeyList
    KeyList = []


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

    def OnShowGenRefls(event):
        '''Generate the reflections from the unit cell and
        display them in the console window
        '''
        OnHklShow(None,indexFrom=' Indexing from unit cell & symmetry settings')
        for r in G2frame.HKL:
            print("{0:.0f},{1:.0f},{2:.0f}   2\u03B8={4:7.3f} d={3:8.4f}".format(*r))

    def OnHklShow(event=None,Print=True,Plot=True,indexFrom=''):
        '''Compute the location of powder diffraction peaks from the
        cell in controls[6:12] and the space group in ssopt['SGData'] if
        defined, or controls[13], if not.

        Reflections are placed in G2frame.HKL

        If cellDisplayOpts['showExtinct'] is True, all reflections are computed
        and reflections that are not present in the G2frame.HKL array are
        placed in G2frame.Extinct

        :params Print: if True, the M20 value, etc. is printed on the console
        :returns: None or [Symb,False,M20,X20,Nhkl,frfnd] where
         * Symb: Space group symbol
         * M20: line position fit metric
         * X20: number of indexed lines fit metric
         * Nhkl: number of generated reflections below dmin
         * frfnd: fraction of lines indexed
        '''
        result = None
        PatternId = G2frame.PatternId
        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Index Peak List'))
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List'))
        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
        Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))[1]
        if 'T' in Inst['Type'][0]:
            dmin = G2lat.Pos2dsp(Inst,Limits[0])
        else:
            dmin = G2lat.Pos2dsp(Inst,Limits[1])
        cell = controls[6:12]
        A = G2lat.cell2A(cell)
        spc = controls[13]
        if ssopt.get('Use',False):  #modulated
            SGData = ssopt.get('SGData',G2spc.SpcGroup(spc)[1])
        else:                       #not modulated
            SGData = G2spc.SpcGroup(spc)[1]
        Symb = SGData['SpGrp']
        M20 = X20 = 0.
        if ssopt.get('Use',False) and ssopt.get('ssSymb',''):
            # modulated is set -- and a super-space group symbol is provided
            SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'])[1]
            if SSGData is None:
                SSGData = G2spc.SSpcGroup(SGData,ssopt['ssSymb'][:-1])[1]     #skip trailing 's' for mag.
            Symb = SSGData['SSpGrp']
            Vec = ssopt['ModVec']
            maxH = ssopt['maxH']
            G2frame.HKL = G2pwd.getHKLMpeak(dmin,Inst,SGData,SSGData,Vec,maxH,A)
            if len(peaks[0]):
                peaks = [G2indx.IndexSSPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #keep esds from peak fit
                M20,X20 = G2indx.calc_M20SS(peaks[0],G2frame.HKL)
        else:
            G2frame.HKL = G2pwd.getHKLpeak(dmin,SGData,A,Inst)
            G2frame.Extinct = []
            if cellDisplayOpts['showExtinct']:
                # generate a table of extinct reflections -- not all, just those not
                # close to a allowed peak
                allpeaks = G2pwd.getHKLpeak(dmin,G2spc.SpcGroup('P 1')[1],A,Inst)
                alreadyShown = G2frame.HKL[:,4].round(3)
                for peak in allpeaks: # show one reflection only if in a region with no others
                    pos = peak[4].round(3)
                    if pos in alreadyShown: continue
                    alreadyShown = np.append(alreadyShown,pos)
                    G2frame.Extinct.append(peak)
            if len(peaks[0]): # put hkl values into the Index Peak List
                peaks = [G2indx.IndexPeaks(peaks[0],G2frame.HKL)[1],peaks[1]]   #keep esds from peak fit
                M20,X20 = G2indx.calc_M20(peaks[0],G2frame.HKL)
                G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Index Peak List'),peaks)
        G2frame.HKL = np.array(G2frame.HKL)
        frfnd = 0.0
        Nhkl = len(G2frame.HKL)
        if Nhkl:
            frfnd = min(1.0,float(len(peaks[0]))/Nhkl)
            if Print:
                print (' new M20,X20: %.2f %d, fraction found: %.3f for %s'%(M20,X20,frfnd,Symb))
            result = [Symb,False,M20,X20,Nhkl,frfnd]
        if not Plot: return result
        if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame,indexFrom=indexFrom)
        else:
            G2pwpl.PlotPatterns(G2frame,indexFrom=indexFrom)
        return result

    def OnSortCells(event):
        '''Sort the display of found unit cells
        '''
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(UnitCellsId)
        c =  event.GetCol()
        if colLabels[c] == 'M20':
            cells = G2indx.sortM20(cells)
        elif colLabels[c] in ['X20','Bravais','a','b','c','alpha','beta','gamma','Volume']:
            if c == 1:
                c += 1  #X20 before Use
            cells = G2indx.sortCells(cells,c-1)     #an extra column (Use) not in cells
        else:
            return
        data = [controls,bravais,cells,dmin,ssopt,magcells]
        G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data,showUse=True)

    def CopyUnitCell(event):
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(UnitCellsId)
        controls = controls[:5]+10*[0.,]
        if len(cells):
            for Cell in cells:
                if Cell[-2]:
                    break
            cell = Cell[2:10]
            controls[4] = 1
            controls[5] = bravaisSymb[cell[0]]
            controls[6:13] = cell[1:8]
            controls[13] = spaceGroups[bravaisSymb.index(controls[5])]
            ssopt['SgResults'] = []
            # G2frame.dataWindow.RefineCell.Enable(True) # set in UpdateUnitCellsGrid
        elif magcells:
            for phase in magcells:
                if phase['Use']:
                    break
            SGData = phase['SGData']
            controls[4] = 1
            controls[5] = (SGData['SGLatt']+SGData['SGLaue']).replace('-','')
            if controls[5][1:] == 'm3': controls[5] += 'm'
            if 'P3' in controls[5] or 'P-3' in controls[5]: controls[5] = 'P6/mmm'
            if 'R' in controls[5]: controls[5] = 'R3-H'
            controls[6:13] = phase['Cell']
            controls[13] = SGData['SpGrp']
            ssopt['SGData'] = SGData
        data = [controls,bravais,cells,dminx,ssopt,magcells]
        G2frame.dataWindow.RunSubGroups.Enable(True)
        G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
        OnHklShow(None,indexFrom=' Indexing from new unit cell & symmetry settings')
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def LoadUnitCell(event):
        '''Called in response to a Load Phase menu command'''
        UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
        data = G2frame.GPXtree.GetItemPyData(UnitCellsId)
        if len(data) < 5:
            data.append({})
        controls,bravais,cells,dminx,ssopt = data[:5]
        magcells = []           #clear away old mag cells list (if any)
        controls = controls[:14]+[['0','0','0',' ',' ',' '],[],]
        data = controls,bravais,cells,dminx,ssopt,magcells
        G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
        pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Phases')
        if not pId: return
        Phases = []
        item, cookie = G2frame.GPXtree.GetFirstChild(pId)
        while item:
            pName = G2frame.GPXtree.GetItemText(item)
            Phase = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,pName))
            if not Phase['General']['SGData'].get('SGFixed',False):
                Phases.append(G2frame.GPXtree.GetItemText(item))
            item, cookie = G2frame.GPXtree.GetNextChild(pId, cookie)
        if not len(Phases):
                wx.MessageBox('NB: Magnetic phases from mcif files are not suitable for this purpose,\n because of space group symbol - operators mismatches',
                    caption='No usable space groups',style=wx.ICON_EXCLAMATION)
                return
        elif len(Phases) == 1: # don't ask questions when there is only 1 answer
            pNum = 0
        else:
            pNum = G2G.ItemSelector(Phases,G2frame,'Select phase',header='Phase')
        if pNum is None: return
        Phase = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,Phases[pNum]))
        Phase['magPhases'] = G2frame.GPXtree.GetItemText(G2frame.PatternId)    #use as reference for recovering possible phases
        Cell = Phase['General']['Cell']
        SGData.update(Phase['General']['SGData'])
        if 'SGGray' not in SGData:
            SGData['SGGray'] = False
        if Phase['General']['Type'] == 'nuclear' and 'MagSpGrp' in SGData:
            SGData.update(G2spc.SpcGroup(SGData['SpGrp'])[1])
        G2frame.dataWindow.RunSubGroups.Enable(True)
        ssopt.update({'Use':False,'ssSymb':'(abg)','ModVec':[0.1,0.1,0.1],'maxH':1,'SgResults':[],})
        if 'SuperSg' in Phase['General'] or SGData.get('SGGray',False):
            ssopt.update({'SGData':SGData,'ssSymb':Phase['General']['SuperSg'],'ModVec':Phase['General']['SuperVec'][0],'Use':True,'maxH':1})
            ssopt['ssSymb'] = ssopt['ssSymb'].replace(',','')
            ssSym = ssopt['ssSymb']
            if SGData.get('SGGray',False):
                ssSym = ssSym[:-1]
            if ssSym not in G2spc.SSChoice(SGData):
                ssSym = ssSym.split(')')[0]+')000'
                ssopt['ssSymb'] = ssSym
                wx.MessageBox('Super space group '+SGData['SpGrp']+ssopt['ssSymb']+' not valid;\n It is set to '+ssSym,
                    caption='Unusable super space group',style=wx.ICON_EXCLAMATION)
            G2frame.dataWindow.RunSubGroups.Enable(False)
        SpGrp = SGData['SpGrp']
        if 'mono' in SGData['SGSys']:
            SpGrp = G2spc.fixMono(SpGrp)
            if SpGrp == None:
                wx.MessageBox('Monoclinic '+SGData['SpGrp']+' not usable here',caption='Unusable space group',style=wx.ICON_EXCLAMATION)
                return
        controls[13] = SpGrp
        controls[4] = 1
        controls[5] = (SGData['SGLatt']+SGData['SGLaue']).replace('-','')
        if controls[5][1:] == 'm3': controls[5] += 'm'
        if controls[5] in ['P3','P3m1','P31m','P6/m']: controls[5] = 'P6/mmm'
        if controls[5] in ['P4/m',]: controls[5] = 'P4/mmm'
        if 'R' in controls[5]: controls[5] = 'R3-H'
        controls[6:13] = Cell[1:8]
        cx,ct,cs,cia = Phase['General']['AtomPtrs']
        controls[15] = [atom[:cx+3] for atom in Phase['Atoms']]
        if 'N' in Inst['Type'][0]:
            if not ssopt.get('Use',False):
                G2frame.dataWindow.RunSubGroupsMag.Enable(True)
        data = controls,bravais,cells,dminx,ssopt,magcells
        G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
        G2frame.dataWindow.RefineCell.Enable(True)
        OnHklShow(None,indexFrom=' Indexing from loaded unit cell & symmetry settings')
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    # TODO: I think this used to work, but this needs to be revisited due to
    # AttributeError: 'G2DataWindow' object has no attribute 'ReImportMenuId'
    # from: 
    #    reqrdr = G2frame.dataWindow.ReImportMenuId.get(event.GetId())
    #
    # def ImportUnitCell(event):
    #     controls,bravais,cells,dminx,ssopt = G2frame.GPXtree.GetItemPyData(UnitCellsId)[:5]
    #     reqrdr = G2frame.dataWindow.ReImportMenuId.get(event.GetId())
    #     rdlist = G2frame.OnImportGeneric(reqrdr,
    #         G2frame.ImportPhaseReaderlist,'phase')
    #     if len(rdlist) == 0: return
    #     rd = rdlist[0]
    #     Cell = rd.Phase['General']['Cell']
    #     SGData = rd.Phase['General']['SGData']
    #     if '1 1' in SGData['SpGrp']:
    #         wx.MessageBox('Unusable space group',caption='Monoclinic '+SGData['SpGrp']+' not usable here',style=wx.ICON_EXCLAMATION)
    #         return
    #     controls[4] = 1
    #     controls[5] = (SGData['SGLatt']+SGData['SGLaue']).replace('-','')
    #     if controls[5][1:] == 'm3': controls[5] += 'm'
    #     if 'P3' in controls[5] or 'P-3' in controls[5]: controls[5] = 'P6/mmm'
    #     if 'R' in controls[5]: controls[5] = 'R3-H'
    #     controls[6:13] = Cell[1:8]
    #     controls[13] = SGData['SpGrp']
    #     ssopt['SgResults'] = []
    #     G2frame.dataWindow.RefineCell.Enable(True)
    #     OnHklShow(None,indexFrom=' Indexing from imported unit cell & symmetry settings')
    #     wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def onRefineCell(event):
        data = RefineCell(G2frame)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnIndexPeaks(event):
        PatternId = G2frame.PatternId
        #print ('Peak Indexing')
        keepcells = []
        try:
            controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List'))
            ssopt['SgResults'] = []
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
        else:
            ibrav = bravais.index(True)
        if not len(peaks[0]):
            G2frame.ErrorDialog('Error','Index Peak List is empty')
            return
        if len(peaks[0][0]) > 9:
            G2frame.ErrorDialog('Error','You need to reload Index Peaks List first')
            return
        G2frame.dataWindow.CopyCell.Enable(False)
        G2frame.dataWindow.RefineCell.Enable(False)
        dlg = wx.ProgressDialog("Generated reflections",'0 '+" cell search for "+bravaisNames[ibrav],101,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_SKIP|wx.PD_CAN_ABORT) #desn't work in 32 bit versions
#            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        try:
            OK,dmin,newcells = G2indx.DoIndexPeaks(peaks[0],controls,bravais,dlg,G2frame.ifX20)
        finally:
            dlg.Destroy()
        cells = keepcells+newcells
        cells = G2indx.sortM20(cells)
        if OK:
            cells[0][10] = True         #select best M20
            data = [controls,bravais,cells,dmin,ssopt,magcells]
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List'),data)
            bestCell = cells[0]
            if bestCell[0] > 10.:
                G2frame.HKL = G2lat.GenHBravais(dmin,bestCell[2],G2lat.cell2A(bestCell[3:9]),ifList=True)
                for hkl in G2frame.HKL:
                    hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
                G2frame.HKL = np.array(G2frame.HKL)
                if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
                    G2plt.PlotPowderLines(G2frame,indexFrom=' Indexing from result #0, M20= %.3f'%cells[0][0])
                else:
                    G2pwpl.PlotPatterns(G2frame,indexFrom=' Indexing from result #0, M20= %.3f'%cells[0][0])
            G2frame.dataWindow.CopyCell.Enable(True)
            G2frame.dataWindow.IndexPeaks.Enable(True)
            G2frame.dataWindow.MakeNewPhase.Enable(True)
            G2frame.ifX20 = True
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def SeaResSelected(event=None):
        '''Responds when "show" is pressed in the UnitCellsTable (search results table)
        or at end of UpdateUnitCellsGrid (When callSeaResSelected==True; event==None).

        Generates & plots reflections
        '''
        data = G2frame.GPXtree.GetItemPyData(UnitCellsId)
        cells,dminx = data[2:4]
        if event is None:
            r,c = 0,0
            colLabel = 'show'
        else:
            r,c =  event.GetRow(),event.GetCol()
            colLabel = event.GetEventObject().GetColLabelValue(c)
        if cells:
            clearShowFlags()
            disableCellEntries()
            Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
                    G2frame,G2frame.PatternId,'Instrument Parameters'))[0]
            if cells[0][0] == '?':  # k-vector table
                # uncheck all rows then check only the one used row
                for i in range(len(cells)):
                    UnitCellsTable.SetValue(i,c,False)
                UnitCellsTable.SetValue(r,c,True)

                if G2frame.kvecSearch['mode']:
                    phase_sel = G2frame.kvecSearch['phase']
                    _, Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
                    Phase = Phases[phase_sel]

                    lat_type = Phase["General"]["SGData"]["SGLatt"]
                    lat_sym = Phase["General"]["SGData"]["SGSys"]
                    if lat_sym == "trigonal":
                        brav_sym = "hR"
                    else:
                        brav_sym = lat_sym[0] + lat_type

                    Trans = np.eye(3)
                    Uvec = np.zeros(3)
                    Vvec = np.zeros(3)

                    newPhase = copy.deepcopy(Phase)
                    newPhase['ranId'] = ran.randint(0, sys.maxsize)
                    newPhase['General']['SGData'] = G2spc.SpcGroup('P 1')[1]
                    newPhase, _ = G2lat.TransformPhase(Phase, newPhase, Trans,
                        Uvec, Vvec, False)
                    atoms_pointer = newPhase['General']['AtomPtrs']

                    atom_coords = list()
                    atom_types = list()
                    for atom in newPhase["Atoms"]:
                        coord_tmp = atom[atoms_pointer[0]:atoms_pointer[0] + 3]
                        atom_coords.append(coord_tmp)
                        type_tmp = atom[atoms_pointer[1]]
                        atom_types.append(type_tmp)

                    atom_ids = kvs.unique_id_gen(atom_types)

                    cell_params = newPhase["General"]["Cell"][1:7]
                    lat_vectors = kvs.lat_params_to_vec(cell_params)

                    # Generate (hkl) coordinate in primitive reciprocal space.
                    # 0-5 for each direction should be generating enough
                    # satellite peaks for users to check.
                    hkl_refls = list()
                    for i in range(6):
                        for j in range(6):
                            for k in range(6):
                                hkl_refls.append([i, j, k])

                    k_search = kvs.kVector(brav_sym, lat_vectors,
                        atom_coords, atom_ids,hkl_refls, [0.], 0.)
                    kpoint = cells[r][3:6]
                    kpoint = k_search.hklConvToPrim(kpoint)

                    rep_prim_latt = k_search.kpathFinder()["reciprocal_primitive_lattice"]

                    satellite_peaks = list()
                    for nucp in k_search.nucPeaks:
                        all_zero_nuc = all(v == 0 for v in nucp)
                        all_zero_k = all(v == 0 for v in kpoint)
                        if all_zero_nuc and all_zero_k:
                            continue

                        hkl_prim = np.array(nucp)
                        hkl_p_k = hkl_prim + kpoint
                        k_cart = np.matmul(hkl_p_k,rep_prim_latt)
                        d_hkl_p_k = 2. * np.pi / np.linalg.norm(k_cart)
                        satellite_peaks.append(d_hkl_p_k)

                        hkl_m_k = hkl_prim - kpoint
                        k_cart = np.matmul(hkl_m_k,rep_prim_latt)
                        d_hkl_m_k = 2. * np.pi / np.linalg.norm(k_cart)
                        satellite_peaks.append(d_hkl_m_k)

                    satellite_peaks = sorted(list(set(satellite_peaks)),reverse=True)


                    # Here we are generating a dummy HKL list to host the
                    # satellite peak positions for the selected k vector. We
                    # set (hkl) for all the peaks to (000).
                    G2frame.HKL = list()
                    for d in satellite_peaks:
                        list_tmp = [0, 0, 0,d, G2lat.Dsp2pos(Inst, d),-1]
                        G2frame.HKL.append(list_tmp)
                    G2frame.HKL = np.array(G2frame.HKL)

                    G2pwpl.PlotPatterns(G2frame)
            if colLabel == 'show':
                for i in range(len(cells)):
                    cells[i][-2] = False
                    UnitCellsTable.SetValue(i,c,False)
                UnitCellsTable.SetValue(r,c,True)
                gridDisplay.ForceRefresh()
                cells[r][-2] = True
                ibrav = cells[r][2]
                try:
                    A = G2lat.cell2A(cells[r][3:9])
                    G2frame.HKL = G2lat.GenHBravais(dmin,ibrav,A,ifList=True)
                    for hkl in G2frame.HKL:
                        hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
                    G2frame.HKL = np.array(G2frame.HKL)
                    if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
                        G2plt.PlotPowderLines(G2frame,indexFrom=' Indexing selection #%d M20= %.3f'%(r,cells[r][0]))
                    else:
                        G2pwpl.PlotPatterns(G2frame,indexFrom=' Indexing selection #%d M20= %.3f'%(r,cells[r][0]))
                except:
                    pass
            elif colLabel == 'Keep':
                if UnitCellsTable.GetValue(r,c):
                    UnitCellsTable.SetValue(r,c,False)
                    cells[r][c] = False
                else:
                    cells[r][c] = True
                    UnitCellsTable.SetValue(r,c,True)
                gridDisplay.ForceRefresh()

            G2frame.GPXtree.SetItemPyData(UnitCellsId,data)

    def MakeNewPhase(event):
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PhaseName = dlg.GetValue()
                cells = G2frame.GPXtree.GetItemPyData(UnitCellsId)[2]
                for Cell in cells:
                    if Cell[-2]:
                        break
                cell = Cell[2:10]
                sub = G2gd.FindPhaseItem(G2frame)
                sub = G2frame.GPXtree.AppendItem(parent=sub,text=PhaseName)
                E,SGData = G2spc.SpcGroup(controls[13])
                G2frame.GPXtree.SetItemPyData(sub, \
                    G2obj.SetNewPhase(Name=PhaseName,SGData=SGData,cell=cell[1:],Super=ssopt))
                G2frame.GetStatusBar().SetStatusText('Change space group from '+str(controls[13])+' if needed',1)
        finally:
            dlg.Destroy()

    def TransformUnitCell(event):
        Trans = np.eye(3)
        Uvec = np.zeros(3)
        Vvec = np.zeros(3)
        ifMag = False
        Type = 'nuclear'
        BNSlatt = ''
        E,SGData = G2spc.SpcGroup(controls[13])
        phase = {'General':{'Name':'','Type':Type,'Cell':['',]+controls[6:13],'SGData':SGData}}
        dlg = G2phsG.TransformDialog(G2frame,phase,Trans,Uvec,Vvec,ifMag,BNSlatt)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                newPhase,Trans,Uvec,Vvec,ifMag,ifConstr,Common = dlg.GetSelection()
                sgData = newPhase['General']['SGData']
                controls[5] = sgData['SGLatt']+sgData['SGLaue']
                controls[13] = sgData['SpGrp']
                ssopt['SGData'] = sgData
                controls[6:13] = newPhase['General']['Cell'][1:8]
            else:
                return
        finally:
            dlg.Destroy()
        OnHklShow(None,indexFrom=' Indexing from transformed unit cell & symmetry settings')
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnLatSym(event):
        'Run Bilbao PseudoLattice cell search'
        # look up a space group matching Bravais lattice (should not matter which one)
        bravaisSPG = {'Fm3m':225,'Im3m':229,'Pm3m':221,'R3-H':146,'P6/mmm':191,
                       'I4/mmm':139,'P4/mmm':123,'Fmmm':69,'Immm':71,
                       'Cmmm':65,'Pmmm':47,'C2/m':12,'P2/m':10,'P1':2}
        pUCid = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(pUCid)
        sgNum = bravaisSPG.get(controls[5],0)
        if sgNum < 1:
            wx.MessageBox('Sorry, only standard cell settings are allowed, please transform axes',caption='Bilbao requires standard settings',style=wx.ICON_EXCLAMATION)
            return
        cell = controls[6:12]
        tolerance = 5.
        dlg = G2G.SingleFloatDialog(G2frame,'Tolerance',
                'Enter angular tolerance for search',5.0,[.1,30.],"%.1f")
        if dlg.ShowModal() == wx.ID_OK:
            tolerance = dlg.GetValue()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        wx.BeginBusyCursor()
        wx.MessageBox(' For use of PSEUDOLATTICE, please cite:\n\n'+
                          G2G.GetCite('Bilbao: PSEUDOLATTICE'),
                          caption='Bilbao PSEUDOLATTICE',
                          style=wx.ICON_INFORMATION)
        page = kSUB.subBilbaoCheckLattice(sgNum,cell,tolerance)
        wx.EndBusyCursor()
        if not page: return
        cells.clear()
        for i,(cell,mat) in enumerate(kSUB.parseBilbaoCheckLattice(page)):
            cells.append([])
            cells[-1] += [mat,0,16]
            cells[-1] += cell
            cells[-1] += [G2lat.calc_V(G2lat.cell2A(cell)),False,False]
        G2frame.GPXtree.SetItemPyData(pUCid,data)
        G2frame.OnFileSave(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnNISTLatSym(event):
        'Run NIST*LATTICE cell search'
        pUCid = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(pUCid)
        nistInput=[0.2,1.,2,3]
        msg = G2G.NISTlatUse(True)
        dlg = wx.Dialog(G2frame,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(dlg,label='NIST*LATTICE Cell Symmetry Search Settings'),
                      0,wx.ALIGN_CENTER_HORIZONTAL,0)
        sizer.Add((-1,15))
        sizer.Add(wx.StaticText(dlg,label=msg))
        sizer.Add((-1,15))
        sizer.Add(wx.StaticText(dlg,label=(
            ('Starting cell: '+3*'{:.3f}, '+3*'{:.2f}, ').format(
                *controls[6:12]))))
        sizer.Add((-1,15))
        tableSizer = wx.FlexGridSizer(0,2,0,0)
        tableSizer.Add(wx.StaticText(dlg,label='Cell length tolerance (A) '),
            0,WACV|wx.ALIGN_LEFT)
        w = G2G.ValidatedTxtCtrl(dlg,nistInput,0,nDig=(6,2))
        tableSizer.Add(w)
        tableSizer.Add(wx.StaticText(dlg,label='Cell angle tolerance (deg) '),
            0,WACV|wx.ALIGN_LEFT)
        w = G2G.ValidatedTxtCtrl(dlg,nistInput,1,nDig=(6,1))
        tableSizer.Add(w)
        #
        # next option makes it too easy to create a really
        # long-running (infinite?) computation. Removed for now.
        #
        #tableSizer.Add(wx.StaticText(dlg,label='Cell volume range (ratio) '),
        #    0,WACV|wx.ALIGN_LEFT)
        #w = G2G.ValidatedTxtCtrl(dlg,nistInput,2)
        #tableSizer.Add(w)
        tableSizer.Add(wx.StaticText(dlg,label='Search mode: Generate '),
            0,WACV|wx.ALIGN_LEFT)
        tableSizer.Add(G2G.EnumSelector(dlg,nistInput,3,
            ['supercells', 'subcells', 'sub- and supercells'],[1,2,3]))
        sizer.Add(tableSizer,1,wx.EXPAND)
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(dlg, wx.ID_OK)
        btn.SetDefault()
        btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
        btnsizer.AddButton(btn)
        btn = wx.Button(dlg, wx.ID_CANCEL)
        btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_CANCEL))
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
        dlg.SetSizer(sizer)
        sizer.Fit(dlg)
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        tol = 3*[nistInput[0]]+3*[nistInput[1]]
        cell = controls[6:12]
        center = controls[13].strip()[0]
        delta = nistInput[2]
        mode = nistInput[3]
        wx.BeginBusyCursor()
        import nistlat
        out = nistlat.CellSymSearch(cell, center, tolerance=tol, mode=mode,deltaV=delta)
        wx.EndBusyCursor()

        if not out: return
        cells.clear()
        for o in out:
            cells.append([])
            c = o[2][0]
            # assign a Laue class
            laue = 16 # P1
            if c[0] == c[1] == c[2] and c[3] == c[4] == c[5] == 90:
                if o[2][1] == 'F':
                    laue = 0 # Fm3m
                elif o[2][1] == 'I':
                    laue = 1 # Im3m
                else:
                    laue = 2 # Pm3m
            elif o[2][1] == 'R':
                laue = 3 # R3
            elif c[0] == c[1] and c[5] == 120:
                laue = 4 # P6/mmm
            elif c[0] == c[1] and c[3] == c[4] == c[5] == 90 and o[2][1] == 'I':
                laue = 5 # I4/mmm
            elif c[0] == c[1] and c[3] == c[4] == c[5] == 90 and o[2][1] == 'P':
                laue = 6 # P4/mmm
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'F':
                laue =  7 # 'Fmmm'
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'I':
                laue =  8 # 'Immm'
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'A':
                laue =  9 # 'Ammm'
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'B':
                laue =  10 # 'Bmmm'
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'C':
                laue =  11 # 'Cmmm'
            elif c[3] == c[4] == c[5] == 90 and o[2][1] == 'P':
                laue =  12 # 'Pmmm'
            elif c[3] == c[5] == 90 and o[2][1] == 'C':
                laue =  13 # 'C2/m'
            elif c[3] == c[5] == 90 and o[2][1] == 'P':
                laue =  14 # 'P2/m'
            elif o[2][1] == 'C':
                laue =  15 # 'C1'
            cells[-1] += [o[4],0,laue]
            cells[-1] += c
            cells[-1] += [G2lat.calc_V(G2lat.cell2A(c)),False,False]
        G2frame.GPXtree.SetItemPyData(pUCid,data)
        G2frame.OnFileSave(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnRunSubs(event):
        G2frame.dataWindow.RunSubGroupsMag.Enable(False)
        pUCid = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(pUCid)
        E,SGData = G2spc.SpcGroup(controls[13])
        Kx = [' ','0','1/2','-1/2','1/3','-1/3','2/3','1']
        Ky = [' ','0','1/2','1/3','2/3','1']
        Kz = [' ','0','1/2','3/2','1/3','2/3','1']
        kvec = [['0','0','0'],[' ',' ',' '],[' ',' ',' ',' ']]
        dlg = G2G.MultiDataDialog(G2frame,title='SUBGROUPS options',prompts=[' k-vector 1',' k-vector 2',' k-vector 3', \
            ' Use whole star',' Filter by','preserve axes','max unique'],
            values=kvec+[False,'',True,100],
            limits=[[Kx[1:],Ky[1:],Kz[1:]],[Kx,Ky,Kz],[Kx,Ky,Kz],[True,False],['',' Landau transition',' Only maximal subgroups',],
                [True,False],[1,100]],
            formats=[['choice','choice','choice'],['choice','choice','choice'],['choice','choice','choice'],'bool','choice',
                    'bool','%d',])
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            magcells = []
            newVals = dlg.GetValues()
            kvec[:9] = newVals[0]+newVals[1]+newVals[2]+[' ',]
            nkvec = kvec.index(' ')
            star = newVals[3]
            filterby = newVals[4]
            keepaxes = newVals[5]
            maxequiv = newVals[6]
            if 'maximal' in filterby:
                maximal = True
                Landau = False
            elif 'Landau' in filterby:
                maximal = False
                Landau = True
            else:
                maximal = False
                Landau = False
            if nkvec not in [0,3,6,9]:
                wx.MessageBox('Error: check your propagation vector(s)',
                    caption='Bilbao SUBGROUPS setup error',style=wx.ICON_EXCLAMATION)
                return
            if nkvec in [6,9] and Landau:
                wx.MessageBox('Error, multi k-vectors & Landau not compatible',
                    caption='Bilbao SUBGROUPS setup error',style=wx.ICON_EXCLAMATION)
                return
            wx.BeginBusyCursor()
            wx.MessageBox(' For use of SUBGROUPS, please cite:\n\n'+
                              G2G.GetCite('Bilbao: k-SUBGROUPSMAG'),
                              caption='Bilbao SUBGROUPS',
                              style=wx.ICON_INFORMATION)
            SubGroups,baseList = kSUB.GetNonStdSubgroups(SGData,kvec[:9],star,Landau)
            wx.EndBusyCursor()
            if SubGroups is None:
                wx.MessageBox('Check your internet connection?',caption='Bilbao SUBGROUPS error',style=wx.ICON_EXCLAMATION)
                return
            if not SubGroups:
                if Landau:
                    wx.MessageBox('No results from SUBGROUPS, multi k-vectors & Landau not compatible',
                        caption='Bilbao SUBGROUPS error',style=wx.ICON_EXCLAMATION)
                else:
                    wx.MessageBox('No results from SUBGROUPS, check your propagation vector(s)',
                        caption='Bilbao SUBGROUPS error',style=wx.ICON_EXCLAMATION)
                return
            controls[14] = kvec[:9]
            try:
                controls[16] = baseList
            except IndexError:
                controls.append(baseList)
            dlg = wx.ProgressDialog('SUBGROUPS results','Processing '+SubGroups[0][0],len(SubGroups),
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)
            for ir,result in enumerate(SubGroups):
                dlg.Update(ir,newmsg='Processing '+result[0])
                Trans = np.array(eval(result[1][0]))
                Uvec = np.array(eval(result[1][1]))
                phase = G2lat.makeBilbaoPhase(result,Uvec,Trans)
                phase['gid'] = result[2]
                phase['altList'] = result[3]
                phase['supList'] = eval(result[4])
                RVT = None
                if keepaxes:
                    RVT = G2lat.FindNonstandard(controls,phase)
                if RVT is not None:
                    result,Uvec,Trans = RVT
                phase.update(G2lat.makeBilbaoPhase(result,Uvec,Trans))
                phase['Cell'] = G2lat.TransformCell(controls[6:12],Trans)
                phase['maxequiv'] = maxequiv
                phase['nAtoms'] = len(TestAtoms(phase,controls[15],SGData,Uvec,Trans,maxequiv,maximal))
                magcells.append(phase)
            dlg.Destroy()
            magcells[0]['Use'] = True
            SGData = magcells[0]['SGData']
            A = G2lat.cell2A(magcells[0]['Cell'][:6])
            G2frame.HKL = np.array(G2pwd.getHKLpeak(1.0,SGData,A,Inst))
            G2pwpl.PlotPatterns(G2frame,extraKeys=KeyList)
        data = [controls,bravais,cells,dmin,ssopt,magcells]
        G2frame.GPXtree.SetItemPyData(pUCid,data)
        G2frame.OnFileSave(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnRunSubsMag(event,kvec1=None):
        def strTest(text):
            if '.' in text: # no decimals
                return False
            elif text.strip() in  [' ','0','1','-1','3/2']: # specials
                return True
            elif '/' in text: #process fraction
                nums = text.split('/')
                return (0 < int(nums[1]) < 10) and (0 < abs(int(nums[0])) < int(nums[1]))
            return False

        G2frame.dataWindow.RunSubGroups.Enable(False)
        pUCid = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
        controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(pUCid)
        E,SGData = G2spc.SpcGroup(controls[13])
        try:
            atoms = list(set([atom[1] for atom in controls[15]]))
        except:
            wx.MessageBox('Error: Problem with phase. Use Load Phase 1st.',
                    caption='k-SUBGROUPSMAG setup error: Phase loaded?',style=wx.ICON_EXCLAMATION)
            return
        testAtoms = ['',]+[atom for atom in atoms if (len(G2elem.GetMFtable([atom,],[2.0,])) and atom != 'O')]   #skip "magnetic" O atoms
        kvec = [['0','0','0'],[' ',' ',' '],[' ',' ',' ',' ']]
        msg = f"Run k-SUBGROUPSMAG on space group {SGData['SpGrp']}\n"+'For k-vectors, enter rational fraction, 0, -1 or 1; no decimals'
        dlg = G2G.MultiDataDialog(G2frame,title='k-SUBGROUPSMAG options',
            prompts=[' k-vector 1',' k-vector 2',' k-vector 3',
                     ' Use whole star',' Filter by','preserve axes',
                     'test for mag. atoms','all have moment','max unique'],
            values=kvec+[False,'',True,'',False,100],
            limits=[['0','0','0'],['0','0','0'],['0','0','0'],
                    [True,False],['',' Landau transition',' Only maximal subgroups',],
                [True,False],testAtoms,[True,False],[1,100]],
            testfxns = [[strTest,strTest,strTest],
                        [strTest,strTest,strTest],
                        [strTest,strTest,strTest],
                        None,None,None,None,None,None],
            formats=[['testfxn','testfxn','testfxn'],
                     ['testfxn','testfxn','testfxn'],
                     ['testfxn','testfxn','testfxn'],
                     'bool','choice','bool','choice','bool','%d',],
            header=msg)
        if dlg.ShowModal() == wx.ID_OK:
            magcells = []
            newVals = dlg.GetValues()
            kvec[:9] = newVals[0]+newVals[1]+newVals[2]+[' ',]
            nkvec = kvec.index(' ')
            star = newVals[3]
            filterby = newVals[4]
            keepaxes = newVals[5]
            atype = newVals[6]
            allmom = newVals[7]
            maxequiv = newVals[8]
            if 'maximal' in filterby:
                maximal = True
                Landau = False
            elif 'Landau' in filterby:
                maximal = False
                Landau = True
            else:
                maximal = False
                Landau = False
            if nkvec not in [0,3,6,9]:
                wx.MessageBox('Error: check your propagation vector(s)',
                    caption='Bilbao k-SUBGROUPSMAG setup error',style=wx.ICON_EXCLAMATION)
                return
            if nkvec in [6,9] and Landau:
                wx.MessageBox('Error, multi k-vectors & Landau not compatible',
                    caption='Bilbao k-SUBGROUPSMAG setup error',style=wx.ICON_EXCLAMATION)
                return
            magAtms = [atom for atom in controls[15] if atom[1] == atype]
            wx.BeginBusyCursor()
            wx.MessageBox(
                ' For use of k-SUBGROUPSMAG in GSAS-II, please cite:\n\n'+
                G2G.GetCite('Bilbao: k-SUBGROUPSMAG')+
                '\nand\n'+
                G2G.GetCite('Bilbao+GSAS-II magnetism'),
                caption='Bilbao/GSAS-II Magnetism',
                style=wx.ICON_INFORMATION)

            MAXMAGN,baseList = kSUB.GetNonStdSubgroupsmag(SGData,kvec[:9],star,Landau)
            wx.EndBusyCursor()
            if MAXMAGN is None:
                wx.MessageBox('Check your internet connection?',caption='Bilbao k-SUBGROUPSMAG error',style=wx.ICON_EXCLAMATION)
                return
            if not MAXMAGN:
                if Landau:
                    wx.MessageBox('No results from k-SUBGROUPSMAG, multi k-vectors & Landau not compatible',
                        caption='Bilbao k-SUBGROUPSMAG error',style=wx.ICON_EXCLAMATION)
                else:
                    wx.MessageBox('No results from k-SUBGROUPSMAG, check your propagation vector(s)',
                        caption='Bilbao k-SUBGROUPSMAG error',style=wx.ICON_EXCLAMATION)
                return
            controls[14] = kvec[:9]
            try:
                controls[16] = baseList
            except IndexError:
                controls.append(baseList)
            dlg = wx.ProgressDialog('k-SUBGROUPSMAG results','Processing '+MAXMAGN[0][0],len(MAXMAGN),
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)

            for ir,result in enumerate(MAXMAGN):
                # result is SPGP,BNS,MV,itemList,altList,superList
                dlg.Update(ir,newmsg='Processing '+result[0])
                Trans = np.array(eval(result[2][0]))
                Uvec = np.array(eval(result[2][1]))
                phase = G2lat.makeBilbaoPhase(result[:2],Uvec,Trans,True)
                phase['gid'] = result[3]
                phase['altList'] = result[4]
                phase['supList'] = eval(result[5])
                RVT = None
                if keepaxes:
                    RVT = G2lat.FindNonstandard(controls,phase)
                if RVT is not None:
                    result,Uvec,Trans = RVT
                phase.update(G2lat.makeBilbaoPhase(result,Uvec,Trans,True))
                phase['Cell'] = G2lat.TransformCell(controls[6:12],Trans)
                phase['aType'] = atype
                phase['allmom'] = allmom
                phase['magAtms'] = magAtms
                phase['maxequiv'] = maxequiv
                phase['nAtoms'] = len(TestMagAtoms(phase,magAtms,SGData,Uvec,Trans,allmom,maxequiv,maximal))
                magcells.append(phase)
            dlg.Destroy()
            magcells[0]['Use'] = True
            SGData = magcells[0]['SGData']
            A = G2lat.cell2A(magcells[0]['Cell'][:6])
            G2frame.HKL = np.array(G2pwd.getHKLpeak(1.0,SGData,A,Inst))
            G2pwpl.PlotPatterns(G2frame,extraKeys=KeyList)
        data = [controls,bravais,cells,dmin,ssopt,magcells]
        G2frame.GPXtree.SetItemPyData(pUCid,data)
        G2frame.OnFileSave(event)
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

    def OnISODISTORT_kvec(phase_nam):   # needs attention from Yuanpeng
        '''Search for k-vector using the ISODISTORT web service
        '''
        def _showWebPage(event):
            'Show a web page when the user presses the "show" button'
            import tempfile
            txt = event.GetEventObject().page
            tmp = tempfile.NamedTemporaryFile(suffix='.html',delete=False)
            with open(tmp.name,'w') as fp:
                fp.write(txt.replace('<HEAD>','<head><base href="https://stokes.byu.edu/iso/">',))
            fileList.append(tmp.name)
            G2G.ShowWebPage('file://'+tmp.name,G2frame)
        def showWebtext(txt):
            import tempfile
            tmp = tempfile.NamedTemporaryFile(suffix='.html',delete=False)
            with open(tmp.name,'w') as fp:
                fp.write(txt.replace('<HEAD>','<head><base href="https://stokes.byu.edu/iso/">',))
            fileList.append(tmp.name)
            G2G.ShowWebPage('file://'+tmp.name,G2frame)

        import tempfile
        import re
        import requests
        from exports import G2export_CIF
        from . import ISODISTORT as ISO
        isoformsite = 'https://iso.byu.edu/iso/isodistortform.php'

        #isoscript='isocifform.php'
        #isoCite = ('For use of this please cite\n'+
        #    G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...',wrap=60,indent=5)+
        #    G2G.GetCite('ISODISPLACE',wrap=60,indent=5))

        #latTol,coordTol,occTol = 0.001, 0.01, 0.1
        phaseID = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
        Phase = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
            G2frame,phaseID,phase_nam))
        data = Phase
        #oacomp,occomp = G2mth.phaseContents(data)
        #ophsnam = data['General']['Name']
        fileList = []
        # write a CIF as a scratch file
        obj = G2export_CIF.ExportPhaseCIF(G2frame)
        obj.InitExport(None)
        obj.currentExportType='phase'
        obj.loadTree()
        tmp = tempfile.NamedTemporaryFile(suffix='.cif', delete=False)
        obj.dirname,obj.filename = os.path.split(tmp.name)
        obj.phasenam = data['General']['Name']
        obj.Writer('', data['General']['Name'])
        parentcif = tmp.name
        ISOparentcif = ISO.UploadCIF(parentcif)
        up2 = {'filename': ISOparentcif, 'input': 'uploadparentcif'}
        out2 = requests.post(isoformsite, up2).text
        try:
            pos = out2.index('<p><FORM')
        except ValueError:
            ISO.HandleError(out2)
            return [], []
        data = {}
        while True:
            try:
                posB = out2[pos:].index('INPUT TYPE') + pos
                posF = out2[posB:].index('>') + posB
                items = out2[posB:posF].split('=',3)
                name = items[2].split()[0].replace('"', '')
                if 'isosystem' in name:
                    break
                vals = items[3].replace('"', '')
                data[name] = vals
                pos = posF
            except ValueError:
                break

        # TODO: bring in the searched k-vector.
        # we can use the `use` check box to select the k-vector to use.
        data['input'] = 'kvector'
        data['irrepcount'] = '1'
        data['kvec1'] = ' 1 *GM, k16 (0,0,0)'
        data['nmodstar1'] = '0'
        out3 = requests.post(isoformsite, data=data).text

        try:
            pos = out3.index('irrep1')
        except ValueError:
            ISO.HandleError(out3)
            return [],[]

        def _get_opt_val(opt_name, out):
            opt_pattern = rf'NAME="{opt_name}" VALUE="(.*?)"'
            opt_match = re.search(opt_pattern, out)

            return opt_match.group(1)

        kvec1 = _get_opt_val('kvec1', out3)
        kvecnumber1 = _get_opt_val('kvecnumber1', out3)
        kparam1 = _get_opt_val('kparam1', out3)
        nmodstar1 = _get_opt_val('nmodstar1', out3)
        data["kvec1"] = kvec1
        data["kvecnumber1"] = kvecnumber1
        data["kparam1"] = kparam1
        data["nmodstar1"] = nmodstar1

        pos = out3.index("irrep1")
        pos1 = out3[pos:].index("</SELECT>")
        str_tmp = out3[pos:][:pos1]
        ir_options = re.findall(r'<OPTION VALUE="([^"]+)">([^<]+)</OPTION>',str_tmp)

        for ir_opt, _ in ir_options:
            data["input"] = "irrep"
            data['irrep1'] = ir_opt
            out4 = requests.post(isoformsite, data=data).text

            irrep1 = _get_opt_val('irrep1', out4)
            irrpointer1 = _get_opt_val('irrpointer1', out4)
            data["irrep1"] = irrep1
            data["irrpointer1"] = irrpointer1

            r_pattern = r'<input\s+type=["\']?radio["\']?\s+name=["\']?'
            r_pattern += r'orderparam["\']?\s+value=["\']?([^"\']+)["\']?'
            radio_val_pattern = re.compile(r_pattern, re.IGNORECASE)
            radio_vals = radio_val_pattern.findall(out4)
            cleaned_radio_vals = [value.strip() for value in radio_vals]

            for radio_val in cleaned_radio_vals:
                data["input"] = "distort"
                data["origintype"] = "method2"
                data["orderparam"] = radio_val + '" CHECKED'
                data["isofilename"] = ""
                requests.post(isoformsite, data=data).text

                continue

                # TODO: parse out the distortion info from out5
                # 1. download the CIF file for each search result
                # 2. load the CIF file as a phase and create a project for it
                # 3. we may want to package up all the project files and save
                #    to a user specified location for later use.

        os.unlink(tmp.name)

        return

    def updateCellsWindow(event):
        'called when k-vec mode is selected'
        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)


    def OnClearCells(event):
        'remove previous search results'
        data[2] = []
        data[5] = []
        ssopt['SgResults'] = []
        wx.CallAfter(UpdateUnitCellsGrid, G2frame, data)

    def OnISODIST(event):
        phase_sel = G2frame.kvecSearch['phase']
        if len(phase_sel.strip()) == 0:
            err_title = "Missing parent phase"
            err_msg = "Please select the parent phase from "
            err_msg += "the drop-down list."
            G2G.G2MessageBox(G2frame, err_msg, err_title)

    def _setupResultsGrid(grid):
        '''Turn on scrolling for grid and  register the grid in resizeGrids so that it will be
        adjusted when the window is resized (also when first created).

        This ensure that Search Results tables display nicely
        with scrollbars internal to the table (plays nicely with table
        headers). See _onResizeUnitCellsList for more.
        '''
        grid.SetScrollRate(10,10)
        resizeGrids.append(grid)

    def clearShowFlags():
        '''resets all the "Use" flags in all search tables
        '''
        for i,grid in enumerate(resizeGrids):
            try:
                for c in range(grid.GetNumberCols()):
                    if grid.GetColLabelValue(c) == 'show':
                        for r in range(grid.GetNumberRows()):
                            grid.GetTable().SetValue(r,c,False)
                grid.ForceRefresh()
            except:
                pass

    def disableCellEntries(mode=False):
        '''Called to disable (or enable with mode==True) the widgets
        associated with entering a unit cell into the 2nd Box
        '''
        for item in unitSizerWidgetList:
            try:
                item.Enable(mode)
            except:
                pass

    def _onResizeUnitCellsList(event):
        '''This is called to adjust the sizes of objects in
        the dataWindow display if that window is resized.

        The only task at present is to adjust the size of the
        GSGrid tables (search results) to be at most the new
        width of the window (less 15 pt to leave room for a
        scrollbar) but if the grid is smaller than that width
        add 15 pt for the scroll bar (might not be needed). For
        the table height, make sure that each results table takes
        less than 1/2 of the available vertical space.
        '''
        try:
            G2frame.dataWindow.Layout()
            wid,hgt = G2frame.dataWindow.GetSize()
            for i,grid in enumerate(resizeGrids):
                # save the initial value for BestSize, as it increases as SetMinSize
                # gets changed
                try:
                    gwid,ghgt = grid.initialGetBestSize
                except:
                    grid.initialGetBestSize = grid.GetBestSize()
                    gwid,ghgt = grid.initialGetBestSize
                grid.SetMaxSize((wid-15,int(hgt/2.5))) # leave inside on right room for scrollbar
                # add a bit to the minimum width for the scroll bar (if there is room)
                grid.SetMinSize((min(wid-15,gwid+15),min(int(hgt/2.5),ghgt)))
            event.Skip()
        except: # can fail after window is destroyed
            pass

# Display of Autoindexing controls
    def firstSizer():

        def OnNcNo(event):
            controls[2] = NcNo.GetValue()

        def OnIfX20(event):
            G2frame.ifX20 = x20.GetValue()

        def OnBravais(event):
            Obj = event.GetEventObject()
            bravais[bravList.index(Obj.GetId())] = Obj.GetValue()

        firstSizer = wx.BoxSizer(wx.VERTICAL)
        firstSizer.Add(wx.StaticText(parent=G2frame.dataWindow,
            label='Autoindexing of "Index Peak List" contents',style=wx.ALIGN_CENTER),0,wx.EXPAND)
        firstSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Bravais Lattice(s) for autoindexing trials:'),0)
        firstSizer.Add((5,5),0)
        indentSizer = wx.BoxSizer(wx.HORIZONTAL)
        indentSizer.Add((20,-1))
        littleSizer = wx.FlexGridSizer(0,4,5,5)
        bravList = []
        bravs = zip(bravais,bravaisNames)
        for brav,bravName in bravs:
            bravCk = wx.CheckBox(G2frame.dataWindow,label=bravName)
            bravList.append(bravCk.GetId())
            bravCk.SetValue(brav)
            bravCk.Bind(wx.EVT_CHECKBOX,OnBravais)
            littleSizer.Add(bravCk,0,WACV)

        indentSizer.Add(littleSizer,0)
        firstSizer.Add(indentSizer,0)
        littleSizer = wx.FlexGridSizer(0,5,5,5)
        littleSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' Max Nc/Nobs '),0,WACV)
        NcNo = wx.SpinCtrl(G2frame.dataWindow)
        NcNo.SetRange(2,8)
        NcNo.SetValue(controls[2])
        NcNo.Bind(wx.EVT_SPINCTRL,OnNcNo)
        littleSizer.Add(NcNo,0,WACV)
        littleSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Start Volume '),0,WACV)
        startVol = G2G.ValidatedTxtCtrl(G2frame.dataWindow,controls,3,typeHint=int,xmin=25)
        littleSizer.Add(startVol,0,WACV)
        x20 = wx.CheckBox(G2frame.dataWindow,label='Use M20/(X20+1)?')
        x20.SetValue(G2frame.ifX20)
        x20.Bind(wx.EVT_CHECKBOX,OnIfX20)
        littleSizer.Add(x20,0,WACV)
        firstSizer.Add(littleSizer,0)
        firstSizer.Add((5,5),0)
        return firstSizer

# Display of k-vector search controls
    def kvecSizer():

        def OnKvecSearch(event):
            'Run the k-vector search'

            try:
                import seekpath
                seekpath
            except:
                G2fil.NeededPackage({'magnetic k-vector search':['seekpath']})
                msg = 'Performing a k-vector search requires installation of the Python seekpath package. Use the Help/Add Package... to install that package.'
                dlg = wx.MessageDialog(G2frame, msg,'Install seekpath package')
                try:
                    dlg.ShowModal()
                finally:
                    dlg.Destroy()
                return

            # msg = G2G.NISTlatUse(True)
            _, _, cells, _, _, _ = data
            wx.BeginBusyCursor()

            # grab the satellite peaks. here, gsas-ii will only grab the data
            # for the histogram under selection.
            hist_name = G2frame.GPXtree.GetItemText(G2frame.PatternId)
            Id = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, hist_name)
            peakId = G2gd.GetGPXtreeItemId(G2frame, Id, 'Peak List')
            peakdata = G2frame.GPXtree.GetItemPyData(peakId)
            if len(peakdata["xtraPeaks"]) == 0:
                err_title = "Empty satellite peak list"
                err_msg = "The satellite peak list is empty. Please select the "
                err_msg += "extra peak positions before executing the k vector "
                err_msg += "search."
                G2G.G2MessageBox(G2frame, err_msg, err_title)
                wx.EndBusyCursor()
                return

            # we need to grab the instrument parameter and call gsas ii routine
            # to convert the satellite peaks into d-spacing.
            Id = G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId, 'Instrument Parameters')
            Parms, _ = G2frame.GPXtree.GetItemPyData(Id)
            xtra_peaks_d = list()
            for extra_peak in peakdata["xtraPeaks"]:
                dsp_tmp = G2lat.Pos2dsp(Parms, extra_peak[0])
                xtra_peaks_d.append(dsp_tmp)

            # select a parent phase. error message will be presented in a pop-up
            # window if an invalid selection is made.
            phase_sel = G2frame.kvecSearch['phase']
            if len(phase_sel.strip()) == 0:
                err_title = "Missing parent phase"
                err_msg = "Please select the parent phase from the dropdown."
                G2G.G2MessageBox(G2frame, err_msg, err_title)
                wx.EndBusyCursor()
                return

            # again, gsas ii will only grab the reflection lists for the
            # histogram under selection. also, only those phases that are being
            # used for the selected histogram will be included in the
            # `refDict` variable.
            refDict = G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId, 'Reflection Lists'))
            if phase_sel not in refDict.keys():
                err_title = "Phase selection error"
                err_msg = "The parent phase selected is not used in "
                err_msg += "current histogram. "
                err_msg += "Please select from one of the followings, "
                err_msg += f"{list(refDict.keys())}"
                G2G.G2MessageBox(G2frame, err_msg, err_title)
                wx.EndBusyCursor()
                return

            # grab the search option. By default, we would search over those
            # high symmetry points only.
            kvs_option_map = {"HighSymPts": 0,"HighSymPts & HighSymPaths": 1,"General": 2}
            kvs_option = G2frame.kvecSearch['soption']
            kvs_option = kvs_option_map[kvs_option]
            if kvs_option == 2:
                # grab the search range, which will be used for the option of
                # 'General' only -- see the part below for user selection of the
                # search option.
                try:
                    kx_s = float(G2frame.kvecSearch['kx_step'])
                    ky_s = float(G2frame.kvecSearch['ky_step'])
                    kz_s = float(G2frame.kvecSearch['kz_step'])
                    num_procs = int(G2frame.kvecSearch['num_procs'])
                except ValueError:
                    err_title = "Invalid k grid input"
                    err_msg = "The k gird values should all be float numbers"
                    G2G.G2MessageBox(G2frame, err_msg, err_title)
                    wx.EndBusyCursor()
                    return

                if kx_s <= 0 or ky_s <= 0 or kz_s <= 0:
                    err_title = "Invalid k grid input"
                    err_msg = "The k step is less than or equal to 0. "
                    err_msg += "Please check the input values."
                    G2G.G2MessageBox(G2frame, err_msg, err_title)
                    wx.EndBusyCursor()
                    return
                else:
                    warn_title = "Long execution time expected"
                    warn_msg = "Searching over general k points may take a while, "
                    warn_msg += "usually on the level of serveral hours. "
                    warn_msg += "Do you want to proceed?"
                    dialog = wx.MessageDialog(G2frame,warn_msg,warn_title,
                        wx.OK | wx.CANCEL | wx.ICON_INFORMATION)
                    result = dialog.ShowModal()

                    if result == wx.ID_OK:
                        pass
                    else:
                        dialog.Destroy()
                        wx.EndBusyCursor()
                        return

                kstep = [kx_s, ky_s, kz_s]
            else:
                num_procs = None
                kstep = None

            # grab the user defined tolerance for determining the optimal k vector.
            # refer to the following link for more detailed explanation about this,
            #
            # https://yr.iris-home.net/kvectordoc
            #
            try:
                tol_val = float(G2frame.kvecSearch['tolerance'])
            except ValueError:
                err_title = "Invalid tolerance input"
                err_msg = "The tolerance value should be a float number, "
                err_msg += "representing the instrument resolution "
                err_msg += ("level " + u"\u03B4" + "d/d.")
                G2G.G2MessageBox(G2frame, err_msg, err_title)
                wx.EndBusyCursor()
                return

            _, Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
            Phase = Phases[phase_sel]

            # grab the Bravais lattice type
            #
            # given the lattice type and lattice system, the Bravais lattice type
            # can be determined. see the comparison table in the Wikipedia link
            # below for the correspondence,
            #
            # https://en.wikipedia.org/wiki/Crystal_system
            #
            # here, we need to assign the Bravais lattice with specific names as
            # inputs for the `seekpath` routine.
            lat_type = Phase["General"]["SGData"]["SGLatt"]
            lat_sym = Phase["General"]["SGData"]["SGSys"]
            if lat_sym == "trigonal":
                brav_sym = "hR"
            else:
                brav_sym = lat_sym[0] + lat_type

            # grab all atomic coordinates in the P1 symmetry
            #
            # define some matrix as necessary inputs for generating the P1
            # structure.
            Trans = np.eye(3)
            Uvec = np.zeros(3)
            Vvec = np.zeros(3)

            # expand the structure to P1 symmetry
            newPhase = copy.deepcopy(Phase)
            newPhase['ranId'] = ran.randint(0, sys.maxsize)
            newPhase['General']['SGData'] = G2spc.SpcGroup('P 1')[1]
            newPhase, _ = G2lat.TransformPhase(Phase,newPhase,Trans,Uvec,Vvec,False)
            atoms_pointer = newPhase['General']['AtomPtrs']

            atom_coords = list()
            atom_types = list()
            for atom in newPhase["Atoms"]:
                coord_tmp = atom[atoms_pointer[0]:atoms_pointer[0] + 3]
                atom_coords.append(coord_tmp)
                type_tmp = atom[atoms_pointer[1]]
                atom_types.append(type_tmp)

            # this will turn each of the atom types into a unique integer number,
            # which is required by the `seekpath` routine
            atom_ids = kvs.unique_id_gen(atom_types)

            # grab the parent unit cell and construct the lattice vectors
            cell_params = newPhase["General"]["Cell"][1:7]
            lat_vectors = kvs.lat_params_to_vec(cell_params)

            hkl_refls = list()
            for i in range(6):
                for j in range(6):
                    for k in range(6):
                        hkl_refls.append([i, j, k])

            try:
                # if we choose option-2, we need to use the `kvec_general` module
                # Otherwise, the computation time would be unacceptably long.
                k_search = kvs.kVector(brav_sym,lat_vectors,atom_coords,atom_ids,
                    hkl_refls,xtra_peaks_d,tol_val,option=kvs_option,
                    kstep=kstep,processes=num_procs)
            except ModuleNotFoundError:
                err_title = "Module not found"
                err_msg = "The `kvec_general` module is not found. Please install "
                err_msg += "the module before running the k-vector search with "
                err_msg += "option-2."
                G2G.G2MessageBox(G2frame, err_msg, err_title)
                wx.EndBusyCursor()
                return

            k_opt = k_search.kOptFinder()
            k_opt_dist = k_opt[1]
            ave_dd = k_opt[2]
            max_dd = k_opt[3]
            k_opt = [list(k_search.kVecPrimToConv(k)) for k in k_opt[0]]

            wx.EndBusyCursor()
            cells.clear()

            # display the result
            select = False
            for i, k_v in enumerate(k_opt):
                select = True
                cells.append([])
                laue = 'P1'
                cells[-1] += ['?', 0, laue]
                c = k_v + [k_opt_dist[i], ave_dd[i], max_dd[i]]
                cells[-1] += c
                cells[-1] += [0, False, False]
                # G2frame.OnFileSave(event) # forces save of project
            wx.CallAfter(UpdateUnitCellsGrid, G2frame, data, select)  # refresh & select 1st result

        kvecSizer = wx.BoxSizer(wx.VERTICAL)
        kvecSizer.Add(wx.StaticText(
            parent=G2frame.dataWindow,label='k-Vector Search Mode',style=wx.ALIGN_CENTER),0,wx.EXPAND)
        Histograms, Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        littleSizer1x = wx.BoxSizer(wx.HORIZONTAL)
        littleSizer2 = wx.BoxSizer(wx.HORIZONTAL)
        G2frame.kvecSearch['phase'] = G2frame.kvecSearch.get('phase', '')
        G2frame.kvecSearch['tolerance'] = G2frame.kvecSearch.get('tolerance', 0.0)
        G2frame.kvecSearch['kx_step'] = G2frame.kvecSearch.get('kx_step', 0.01)
        G2frame.kvecSearch['ky_step'] = G2frame.kvecSearch.get('ky_step', 0.01)
        G2frame.kvecSearch['kz_step'] = G2frame.kvecSearch.get('k_step', 0.01)
        G2frame.kvecSearch['num_procs'] = G2frame.kvecSearch.get('num_procs', 8)
        G2frame.kvecSearch['soption'] = G2frame.kvecSearch.get('soption', "HighSymPts")
        if len(Phases) == 0:
            littleSizer.Add(wx.StaticText(G2frame.dataWindow,
                label='    You need to define a phase to use k-vector searching'),0,WACV)
        elif len(Phases) == 1:
            G2frame.kvecSearch['phase'] = list(Phases.keys())[0]
        else:
            littleSizer.Add(wx.StaticText(G2frame.dataWindow, label='Select phase'), 0, WACV)
            ch = G2G.EnumSelector(G2frame.dataWindow, G2frame.kvecSearch, 'phase',
                [''] + list(Phases.keys()))
            littleSizer.Add(ch, 10, WACV | wx.RIGHT, 0)

        littleSizer1x.Add(wx.StaticText(G2frame.dataWindow, label='kx step'),0,WACV)
        kx_s = G2G.ValidatedTxtCtrl(G2frame.dataWindow,G2frame.kvecSearch,'kx_step',
            nDig=[6, 3],typeHint=float,size=(50, -1))
        littleSizer1x.Add(kx_s, 0, WACV | wx.RIGHT, 10)
        littleSizer1x.Add((10, -1))

        littleSizer1x.Add(wx.StaticText(G2frame.dataWindow, label='ky step'),0,WACV)
        ky_s = G2G.ValidatedTxtCtrl(G2frame.dataWindow,G2frame.kvecSearch,'ky_step',
            nDig=[6, 3],typeHint=float,size=(50, -1))
        littleSizer1x.Add(ky_s, 0, WACV | wx.RIGHT, 10)
        littleSizer1x.Add((10, -1))

        littleSizer1x.Add(wx.StaticText(G2frame.dataWindow, label='kz step'),0,WACV)
        kz_s = G2G.ValidatedTxtCtrl(G2frame.dataWindow,G2frame.kvecSearch,'kz_step',
            nDig=[6, 3],typeHint=float,size=(50, -1))
        littleSizer1x.Add(kz_s, 0, WACV | wx.RIGHT, 10)
        littleSizer1x.Add((10, -1))

        littleSizer1x.Add(wx.StaticText(G2frame.dataWindow, label='Number of processors'),0,WACV)
        num_procs = G2G.ValidatedTxtCtrl(G2frame.dataWindow,G2frame.kvecSearch,'num_procs',
            nDig=[6, 2],typeHint=int,size=(50, -1))
        littleSizer1x.Add(num_procs, 0, WACV | wx.RIGHT, 10)

        littleSizer2.Add(wx.StaticText(G2frame.dataWindow, label='Search tolerance'),0,WACV)
        tolVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,G2frame.kvecSearch,'tolerance',
            nDig=[10, 6],typeHint=float,size=(70, -1))
        littleSizer2.Add(tolVal, 0, WACV | wx.RIGHT, 20)

        littleSizer2.Add(wx.StaticText(G2frame.dataWindow, label='Search option'), 0, WACV)
        search_opts = ["HighSymPts", "HighSymPts & HighSymPaths", "General"]
        ch1 = G2G.EnumSelector(G2frame.dataWindow,G2frame.kvecSearch,'soption',search_opts)
        littleSizer2.Add(ch1, 10, WACV | wx.RIGHT, 0)

        kvecSizer.Add(littleSizer, 0)
        kvecSizer.Add((-1, 5), 0)
        kvecSizer.Add(littleSizer1x, 0)
        kvecSizer.Add((-1, 5), 0)
        kvecSizer.Add(littleSizer2, 0)
        kvecSizer.Add((-1, 5), 0)
        btn = wx.Button(G2frame.dataWindow,wx.ID_ANY,'Start Search')
        btn.Bind(wx.EVT_BUTTON,OnKvecSearch)
        kvecSizer.Add(btn)
        return kvecSizer

# Unit cell display controls
    def unitSizer():

        def OnSSselect(event):
            if controls[5] in ['Fm3m','Im3m','Pm3m']:
                SSselect.SetValue(False)
                G2frame.ErrorDialog('Cubic lattice','Incommensurate superlattice not possible with a cubic lattice')
                return
            ssopt['Use'] = SSselect.GetValue()
            if 'ssSymb' not in ssopt:
                ssopt.update({'ssSymb':'(abg)','ModVec':[0.1,0.1,0.1],'maxH':1})
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnSelMG(event):
            ssopt['ssSymb'] = selMG.GetValue()
            Vec = ssopt['ModVec']
            modS = G2spc.splitSSsym(ssopt['ssSymb'])[0]
            ssopt['ModVec'] = G2spc.SSGModCheck(Vec,modS)[0]
            print (' Selecting: '+controls[13]+ssopt['ssSymb']+ 'maxH:'+str(ssopt['maxH']))
            OnHklShow(event,indexFrom=' Indexing from new super group %s'%(controls[13]+ssopt['ssSymb']))
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnModVal(invalid,value,tc):
            OnHklShow(tc.event,indexFrom=' Indexing from new modulation vector')

        def OnMoveMod(event):
            Obj = event.GetEventObject()
            ObjId = Obj.GetId()
            Id,valObj = Indx[ObjId]
            inc = float(shiftChoices[shiftSel.GetSelection()][:-1])
            move = Obj.GetValue()*inc/100.
            Obj.SetValue(0)
            value = min(0.98,max(-0.98,float(valObj.GetValue())+move))
            valObj.SetValue('%.4f'%(value))
            ssopt['ModVec'][Id] = value
            OnHklShow(event,indexFrom=' Indexing from changed modulation vector')

        def OnMaxMH(event):
            ssopt['maxH'] = int(maxMH.GetValue())
            print (' Selecting: '+controls[13]+ssopt['ssSymb']+'maxH:'+str(ssopt['maxH']))
            OnHklShow(event,indexFrom=' Indexing from new max. M')

        def OnButton(xpos,ypos):
            modSym = ssopt['ssSymb'].split(')')[0]+')'
            if modSym in ['(a0g)','(a1/2g)']:
                ssopt['ModVec'][0] = xpos
                ssopt['ModVec'][2] = ypos
            elif modSym in ['(0bg)','(1/2bg)']:
                ssopt['ModVec'][1] = xpos
                ssopt['ModVec'][2] = ypos
            elif modSym in ['(ab0)','(ab1/2)']:
                ssopt['ModVec'][0] = xpos
                ssopt['ModVec'][1] = ypos
            vec = ssopt['ModVec']
            print(' Trying: %s %s modulation vector = %.3f %.3f %.3f'%(controls[13],ssopt['ssSymb'],vec[0],vec[1],vec[2]))
            OnHklShow(None,indexFrom=' Indexing from selected modulation vector')
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnFindOneMV(event):
            Peaks = np.copy(peaks[0])
            print (' Trying: '+controls[13],ssopt['ssSymb']+' maxH: 1')
            dlg = wx.ProgressDialog('Elapsed time','Modulation vector search',
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
            try:
                ssopt['ModVec'],result = G2indx.findMV(Peaks,controls,ssopt,Inst,dlg)
                if len(result[0]) == 2:
                    G2plt.PlotXYZ(G2frame,result[2],1./result[3],labelX='a',labelY='g',newPlot=True,
                        Title='Modulation vector search for %s%s'%(controls[13],ssopt['ssSymb']),
                        buttonHandler=OnButton)
                elif len(result[0]) == 1:
                    G2plt.PlotXY(G2frame,[[result[2],1./result[3]],],labelX='k',labelY='fit',newPlot=True,
                        Title='Modulation vector search for %s%s'%(controls[13],ssopt['ssSymb']))

            finally:
                dlg.Destroy()
            OnHklShow(event,indexFrom=' Indexing from best modulation vector')
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnFindMV(event):

            best = 1.
            bestSS = ''
            ssopt['SSgResults'] = []
            for ssSym in ssChoice:
                Vref = [x for x in ['a','b','g'] if x in ssSym]
                if len(Vref) > 1:
                    print(' %s skipped - too many variables to search this way'%ssSym)
                    continue
                ssopt['ssSymb'] = ssSym
                Peaks = np.copy(peaks[0])
                ssopt['ModVec'] = G2spc.SSGModCheck(ssopt['ModVec'],G2spc.splitSSsym(ssSym)[0],True)[0]
                print (' Trying: '+controls[13]+ssSym+' maxH: 1')
                dlg = wx.ProgressDialog('Elapsed time','Modulation vector search with %s'%ssSym,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
                try:
                    ssopt['ModVec'],result = G2indx.findMV(Peaks,controls,ssopt,Inst,dlg)
                    ssopt['SSgResults'].append([ssSym,ssopt['ModVec'],])
                    OnHklShow(event,indexFrom='  Indexing from best modulation vector')
                finally:
                    dlg.Destroy()
                if result[1] < best:
                    bestSS = ssSym
                    best = result[1]
            if bestSS != '':
                ssopt['ssSymb'] = bestSS
                ssopt['ModVec'],result = G2indx.findMV(Peaks,controls,ssopt,Inst,dlg=None)
                G2plt.PlotXY(G2frame,[[result[2],1./result[3]],],labelX='k',labelY='fit',
                    newPlot=True,Title='Modulation vector search for %s'%bestSS)

            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnBravSel(event):
            brav = bravSel.GetString(bravSel.GetSelection())
            controls[5] = brav
            controls[13] = SPGlist[brav][0]
            ssopt['Use'] = False
            OnHklShow(event,indexFrom=' Indexing from Bravais lattice %s'%brav)
            wx.CallLater(100,UpdateUnitCellsGrid,G2frame,data)

        def OnSpcSel(event):
            controls[13] = spcSel.GetString(spcSel.GetSelection())
            ssopt['SGData'] = G2spc.SpcGroup(controls[13])[1]
            # ssopt['Use'] = False
            G2frame.dataWindow.RefineCell.Enable(True)
            OnHklShow(event,indexFrom=' Indexing from space group %s'%controls[13])
            wx.CallLater(100,UpdateUnitCellsGrid,G2frame,data)

        def OnTryAllSG(event):
            '''Computes extinctions for all possible space groups
            creates data for a table that is displayed in the subsequent
            UpdateUnitCellsGrid call.
            '''
            ssopt['SgResults'] = []
            ssopt['SgSettings'] = ''
            for controls[13] in SPGlist[controls[5]]:
                ssopt['SGData'] = G2spc.SpcGroup(controls[13])[1]
                ssopt['Use'] = False
                G2frame.dataWindow.RefineCell.Enable(True)
                res = OnHklShow(event,False,False,indexFrom=' ')
                if res:
                    ssopt['SgResults'].append(res)
            if ssopt['SgResults']:
                ssopt['SgResults'] = G2mth.sortArray(ssopt['SgResults'],2,reverse=True)
                controls[13] = ssopt['SgResults'][0][0]
                ssopt['SGData'] = G2spc.SpcGroup(controls[13])[1]
                G2frame.dataWindow.RefineCell.Enable(True)
                ssopt['SgResults'][0][1] = True  # select the first entry
                ssopt['SgSettings'] = 'Brav. Lat.: {}, cell: {:.3f} {:.3f} {:.3f} {:.2f} {:.2f} {:.2f}'.format(*controls[5:12])
                OnHklShow(event,False,False,indexFrom=' Indexing from space group %s'%controls[13])
            wx.CallLater(100,UpdateUnitCellsGrid,G2frame,data)

        def SetCellValue(Obj,ObjId,value):
            if hasattr(Obj,'ChangeValue'):
                setval = Obj.ChangeValue
            else:
                setval = Obj.SetValue
            if controls[5] in ['Fm3m','Im3m','Pm3m']:
                controls[6] = controls[7] = controls[8] = value
                controls[9] = controls[10] = controls[11] = 90.0
                setval(controls[6])
            elif controls[5] in ['R3-H','P6/mmm','I4/mmm','P4/mmm']:
                if ObjId == 0:
                    controls[6] = controls[7] = value
                    setval(controls[6])
                else:
                    controls[8] = value
                    setval(controls[8])
                controls[9] = controls[10] = controls[11] = 90.0
                if controls[5] in ['R3-H','P6/mmm']:
                    controls[11] = 120.
            elif controls[5] in ['Fmmm','Immm','Cmmm','Pmmm']:
                controls[6+ObjId] = value
                setval(controls[6+ObjId])
                controls[9] = controls[10] = controls[11] = 90.0
            elif controls[5] in ['I2/m','A2/m','C2/m','P2/m']:
                controls[9] = controls[11] = 90.0
                if ObjId != 3:
                    controls[6+ObjId] = value
                    setval(controls[6+ObjId])
                else:
                    controls[10] = value
                    setval(controls[10])
            else:
                controls[6+ObjId] = value
                if ObjId < 3:
                    setval(controls[6+ObjId])
                else:
                    setval(controls[6+ObjId])
            controls[12] = G2lat.calc_V(G2lat.cell2A(controls[6:12]))
            volVal.SetValue("%.3f"%(controls[12]))

        def OnMoveCell(event):
            Obj = event.GetEventObject()
            ObjId = cellList.index(Obj.GetId())
            valObj = valDict[Obj.GetId()]
            inc = float(shiftChoices[shiftSel.GetSelection()][:-1])
            move = Obj.GetValue()  # +1 or -1
            Obj.SetValue(0)
            value = float(valObj.GetValue()) * (1. + move*inc/100.)
            SetCellValue(valObj,ObjId//2,value)
            OnHklShow(event,indexFrom=' Indexing from new cell')

        def OnCellChange(invalid,value,tc):
            if invalid:
                return
            try: # fails when zero is updated
                SetCellValue(tc,Info[tc.GetId()],value)
            except:
                pass
            OnHklShow(tc.event,indexFrom=' Indexing from new cell ')
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnMagSel(event):
            Obj = event.GetEventObject()
            if Obj.GetValue():
                SGData['SGSpin'] = [1,]*len(SGData['SGSpin'])
                GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
                SGData['GenSym'] = GenSym
                SGData['GenFlg'] = GenFlg
                OprNames,SpnFlp = G2spc.GenMagOps(SGData)
                SGData['SpnFlp'] = SpnFlp
                SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
            else:
                del SGData['MagSpGrp']
            ssopt['SGData'] = SGData
            OnHklShow(event,indexFrom=' Indexing from new spin selection')
            wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnSpinOp(event):
            Obj = event.GetEventObject()
            isym = Indx[Obj.GetId()]+1
            spCode = {'red':-1,'black':1}
            SGData['SGSpin'][isym] = spCode[Obj.GetValue()]
            G2spc.CheckSpin(isym,SGData)
            GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
            SGData['GenSym'] = GenSym
            SGData['GenFlg'] = GenFlg
            OprNames,SpnFlp = G2spc.GenMagOps(SGData)
            SGData['SpnFlp'] = SpnFlp
            SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
            OnHklShow(event,indexFrom=' Indexing from new spin operators')

        def OnBNSlatt(event):
            Obj = event.GetEventObject()
            SGData.update(G2spc.SpcGroup(SGData['SpGrp'])[1])
            BNSlatt = Obj.GetValue()
            if '_' in BNSlatt:
                SGData['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
            else:
                SGData['BNSlattsym'] = [SGData['SGLatt'],[0.,0.,0.]]
            SGData['SGSpin'] = [1,]*len(SGData['SGSpin'])
            GenSym,GenFlg = G2spc.GetGenSym(SGData)[:2]
            SGData['GenSym'] = GenSym
            SGData['GenFlg'] = GenFlg
            SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
            G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
            OprNames,SpnFlp = G2spc.GenMagOps(SGData)
            SGData['SpnFlp'] = SpnFlp
            OnHklShow(event,indexFrom=' Indexing from new BNS centering')

        def OnShowSpins(event):
            msg = 'Magnetic space group information'
            text,table = G2spc.SGPrint(SGData,AddInv=True)
            text[0] = ' Magnetic Space Group: '+SGData['MagSpGrp']
            text[3] = ' The magnetic lattice point group is '+SGData['MagPtGp']
            G2G.SGMagSpinBox(G2frame.dataWindow,msg,text,table,SGData['SGCen'],OprNames,
                SGData['SpnFlp'],False).Show()

        def OnMakePks(event):
            msg = ('This will replace the current contents of the Peaks List'+
                   ' with peak positions generated from the current cell.'+
                   ' Are you sure you want to delete the previous Peak List'+
                   ' contents?')
            dlg = wx.MessageDialog(G2frame, msg,'Replace peaks?',wx.YES_NO|wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result != wx.ID_YES: return
            hkl = G2frame.HKL
            PatternId = G2frame.PatternId
            limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Limits'))[1]
            background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Background'))
            inst,inst2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Instrument Parameters'))
            peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'))
            Pattern = G2frame.GPXtree.GetItemPyData(PatternId)
            peaks = []
            profile = Pattern[1]
            bxye = GetFileBackground(G2frame,profile,background)
            x0 = profile[0]
            iBeg = np.searchsorted(x0,limits[0])
            iFin = np.searchsorted(x0,limits[1])
            x = x0[iBeg:iFin]
            y = (profile[1]-bxye)[iBeg:iFin]
    #        ysig = 1.0*np.std(y0)
            poss = [p for p in hkl[:,4] if (x[0] <= p <= x[-1])]
            mags = y[np.searchsorted(x,poss)]
            refs = list(zip(poss,mags))
            if 'T' in Inst['Type'][0]:
                refs = G2mth.sortArray(refs,0,reverse=True)     #big TOFs first
            else:   #'C', 'E' or 'B'
                refs = G2mth.sortArray(refs,0,reverse=False)    #small 2-Thetas or energies first
            for i,ref1 in enumerate(refs):      #reject picks closer than 1 FWHM
                for ref2 in refs[i+1:]:
                    if abs(ref2[0]-ref1[0]) < 2.*G2pwd.getFWHM(ref1[0],inst):
                        del(refs[i])
            if 'T' in Inst['Type'][0]:
                refs = G2mth.sortArray(refs,1,reverse=False)
            else:   #'C', 'E' or 'B'
                refs = G2mth.sortArray(refs,1,reverse=True)
            for pos,mag in refs:
                peaks.append(G2mth.setPeakparms(inst,inst2,pos,mag))
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Peak List'),peaks)

        def OnNewHklShow(event):
            disableCellEntries(True)
            clearShowFlags()
            OnHklShow(event,indexFrom=' Indexing from unit cell & symmetry settings')

        def OnExtHklShow(event):
            Obj = event.GetEventObject()
            if Obj.GetValue():
                OnHklShow(event,indexFrom=' Indexing including extinct hkls')
            else:
                OnHklShow(event,indexFrom=' Indexing from unit cell & symmetry settings')

        def OnAxHklShow():
            mode = Sel.GetStringSelection()
            if 'None' not in mode:
                OnHklShow(None,indexFrom=' Indexing showing only %s hkls'%mode)
            else:
                OnHklShow(None,indexFrom=' Indexing from unit cell & symmetry settings')

        unitSizer = wx.BoxSizer(wx.VERTICAL)
        unitSizer.Add(wx.StaticText(parent=G2frame.dataWindow,style=wx.ALIGN_CENTER,
            label='Unit Cell && Symmetry Settings for Reflection Display'),0,wx.EXPAND)
        ibrav = SetLattice(controls)
        for cellGUI in cellGUIlist:
            if ibrav in cellGUI[0]:
                useGUI = cellGUI
        cellList = []
        valDict = {}
        Info = {}

        bravSizer = wx.BoxSizer(wx.HORIZONTAL)
        bravSizer.Add(wx.StaticText(G2frame.dataWindow,label=" Bravais  \n lattice ",style=wx.ALIGN_CENTER),0,WACV,5)
        bravSel = wx.Choice(G2frame.dataWindow,choices=bravaisSymb,size=(75,-1))
        unitSizerWidgetList.append(bravSel)
        bravSel.SetSelection(bravaisSymb.index(controls[5]))
        bravSel.Bind(wx.EVT_CHOICE,OnBravSel)
        bravSizer.Add(bravSel,0,WACV)
        bravSizer.Add(wx.StaticText(G2frame.dataWindow,label=" Space  \n group  ",style=wx.ALIGN_CENTER),0,WACV,5)
        spcSel = wx.Choice(G2frame.dataWindow,choices=SPGlist[controls[5]],size=(100,-1))
        unitSizerWidgetList.append(spcSel)
        try:
            spcSel.SetSelection(SPGlist[controls[5]].index(controls[13]))
        except ValueError:
            pass
        spcSel.Bind(wx.EVT_CHOICE,OnSpcSel)
        bravSizer.Add(spcSel,0,WACV)
        bravSizer.Add((5,-1))
        tryAll = wx.Button(G2frame.dataWindow,label='Try all?')
        unitSizerWidgetList.append(tryAll)
        tryAll.Bind(wx.EVT_BUTTON,OnTryAllSG)
        bravSizer.Add(tryAll,0,WACV)
        if 'E' not in Inst['Type'][0]:
            SSselect = wx.CheckBox(G2frame.dataWindow,label="Modulated?")
            unitSizerWidgetList.append(SSselect)
            SSselect.SetValue(ssopt.get('Use',False))
            SSselect.Bind(wx.EVT_CHECKBOX,OnSSselect)
            bravSizer.Add(SSselect,0,WACV)
            if ssopt.get('Use',False):        #zero for super lattice doesn't work!
                controls[0] = False
        unitSizer.Add(bravSizer,0)

        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellSizer = wx.FlexGridSizer(0,min(6,useGUI[1]),3,3)
        for txt,fmt,ifEdit,Id in zip(*useGUI[2]):
            cellSizer.Add(wx.StaticText(G2frame.dataWindow,label=txt,style=wx.ALIGN_RIGHT),0,wx.ALIGN_RIGHT|WACV)
            if ifEdit:          #a,b,c,etc.
                cellVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,controls,6+Id,nDig=fmt,OnLeave=OnCellChange,size=(65,-1))
                unitSizerWidgetList.append(cellVal)
                Info[cellVal.GetId()] = Id
                valSizer = wx.BoxSizer(wx.HORIZONTAL)
                valSizer.Add(cellVal,0,WACV)
                cellSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(20,20))
                unitSizerWidgetList.append(cellSpin)
                cellSpin.SetValue(0)
                cellSpin.SetRange(-1,1)
                cellSpin.Bind(wx.EVT_SPIN, OnMoveCell)
                valSizer.Add(cellSpin,0,WACV)
                cellSizer.Add(valSizer,0,WACV)
                cellList.append(cellVal.GetId())
                cellList.append(cellSpin.GetId())
                valDict[cellSpin.GetId()] = cellVal
            else:               #volume
                volVal = wx.TextCtrl(G2frame.dataWindow,value=(fmt%(controls[12])),style=wx.TE_READONLY,size=(65,-1))
                volVal.SetBackgroundColour(VERY_LIGHT_GREY)
                cellSizer.Add(volVal,0,WACV)
        hSizer.Add(cellSizer,0,WACV)
        hSizer.Add((10,-1))
        vcSizer =  wx.BoxSizer(wx.VERTICAL)
        vcSizer.Add(wx.StaticText(G2frame.dataWindow,label='cell step',
            style=wx.ALIGN_CENTER),0,wx.EXPAND)
        shiftChoices = [ '0.01%','0.05%','0.1%','0.5%', '1.0%','2.5%','5.0%']
        shiftSel = wx.Choice(G2frame.dataWindow,choices=shiftChoices)
        unitSizerWidgetList.append(shiftSel)
        shiftSel.SetSelection(3)
        vcSizer.Add(shiftSel)
        hSizer.Add(vcSizer,0,WACV)
        if not ssopt.get('Use',False):        #zero for super lattice doesn't work!
            hSizer.Add((15,-1))
            vcSizer =  wx.BoxSizer(wx.VERTICAL)
            vcSizer.Add(wx.StaticText(G2frame.dataWindow,label="Zero offset",
                style=wx.ALIGN_CENTER),0,wx.EXPAND)
            hcSizer = wx.BoxSizer(wx.HORIZONTAL)
            zero = G2G.ValidatedTxtCtrl(G2frame.dataWindow,controls,1,nDig=(10,4),typeHint=float,
                xmin=-5.,xmax=5.,size=(50,-1),OnLeave=OnCellChange)
            unitSizerWidgetList.append(zero)
            hcSizer.Add(zero,0,WACV)
            zeroVar = G2G.G2CheckBox(G2frame.dataWindow,'Ref?',controls,0)
            unitSizerWidgetList.append(zeroVar)
            hcSizer.Add(zeroVar,0,WACV|wx.LEFT,3)
            vcSizer.Add(hcSizer)
            hSizer.Add(vcSizer,0,WACV)
        unitSizer.Add(hSizer,0)
        if ssopt.get('Use',False):        #super lattice display
            indChoice = ['1','2','3','4',]
            SpSg = SGData['SpGrp']
            if 'MagSpGrp' in SGData:    #limit to one for magnetic SS for now
                indChoice = ['1',]
                SpSg = SGData['MagSpGrp']
            ssChoice = G2spc.SSChoice(SGData)
            if ssopt['ssSymb'] not in ssChoice:
                ssopt['ssSymb'] = ssopt['ssSymb'][:-1]
            ssSizer = wx.BoxSizer(wx.HORIZONTAL)
            ssSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Supersymmetry space group: '+SpSg+' '),0,WACV)
            selMG = wx.ComboBox(G2frame.dataWindow,value=ssopt['ssSymb'],
                choices=ssChoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            unitSizerWidgetList.append(selMG)
            selMG.Bind(wx.EVT_COMBOBOX, OnSelMG)
            ssSizer.Add(selMG,0,WACV)
            unitSizer.Add(ssSizer,0)
            ssSizer = wx.BoxSizer(wx.HORIZONTAL)
            ssSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Mod. vector: '),0,WACV)
            modS = G2spc.splitSSsym(ssopt['ssSymb'])[0]
            ssopt['ModVec'],ifShow = G2spc.SSGModCheck(ssopt['ModVec'],modS)
            for i,[val,show] in enumerate(zip(ssopt['ModVec'],ifShow)):
                if show:
                    valSizer = wx.BoxSizer(wx.HORIZONTAL)
                    modVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,ssopt['ModVec'],i,
                        xmin=-.98,xmax=.98,nDig=(10,4),typeHint=float,
                        OnLeave=OnModVal,size=wx.Size(50,-1))
                    unitSizerWidgetList.append(modVal)
                    valSizer.Add(modVal,0,WACV)
                    modSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(20,20))
                    unitSizerWidgetList.append(modSpin)
                    modSpin.SetValue(0)
                    modSpin.SetRange(-1,1)
                    modSpin.Bind(wx.EVT_SPIN, OnMoveMod)
                    valSizer.Add(modSpin,0,WACV)
                    ssSizer.Add(valSizer,0,WACV)
                    Indx[modVal.GetId()] = i
                    Indx[modSpin.GetId()] = [i,modVal]
                else:
                    modVal = wx.TextCtrl(G2frame.dataWindow,value=('%.3f'%(val)),
                        size=wx.Size(50,20),style=wx.TE_READONLY)
                    modVal.SetBackgroundColour(VERY_LIGHT_GREY)
                    ssSizer.Add(modVal,0,WACV)
            ssSizer.Add(wx.StaticText(G2frame.dataWindow,label=' in steps of cell step above'),0,WACV)
            unitSizer.Add(ssSizer,0)
            ssSizer = wx.BoxSizer(wx.HORIZONTAL)
            ssSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Max. M: '),0,WACV)
            maxMH = wx.ComboBox(G2frame.dataWindow,value=str(ssopt['maxH']),
                choices=indChoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            unitSizerWidgetList.append(maxMH)
            maxMH.Bind(wx.EVT_COMBOBOX, OnMaxMH)
            ssSizer.Add(maxMH,0,WACV)
            if len(peaks[0]):
                findMV = wx.Button(G2frame.dataWindow,label="Find mod. vec.?")
                findMV.Bind(wx.EVT_BUTTON,OnFindOneMV)
                ssSizer.Add(findMV,0,WACV)
                findallMV = wx.Button(G2frame.dataWindow,label="Try all?")
                unitSizerWidgetList.append(findallMV)
                findallMV.Bind(wx.EVT_BUTTON,OnFindMV)
                ssSizer.Add(findallMV,0,WACV)
            unitSizer.Add(ssSizer,0)
        # cell display options
        littleSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellDisplayOpts['Show'] = True
        hklShow = wx.Button(G2frame.dataWindow,label="show HKL from cell")
        hklShow.Bind(wx.EVT_BUTTON,OnNewHklShow)
        littleSizer.Add(hklShow,0,WACV)
        if 'E' not in Inst['Type'][0]:
            littleSizer.Add((5,-1))
            if not ssopt.get('Use',False):  # Show Extinct not available for super lattice
                showExt = G2G.G2CheckBox(G2frame.dataWindow,'Show Extinct',
                    cellDisplayOpts,'showExtinct',OnChange=OnExtHklShow)
                unitSizerWidgetList.append(showExt)
                littleSizer.Add(showExt,0,WACV)
        littleSizer.Add(wx.StaticText(G2frame.dataWindow,label=' highlight ',style=wx.ALIGN_RIGHT),0,WACV)
        G2frame.PlotOpts['hklHighlight'] = G2frame.PlotOpts.get('hklHighlight',0)
        Sel = G2G.G2ChoiceButton(G2frame.dataWindow,[ 'None',] + [c+notEq0 for c in ('h','k','l')],
            indLoc=G2frame.PlotOpts,indKey='hklHighlight',onChoice=OnAxHklShow)
        unitSizerWidgetList.append(Sel)
        littleSizer.Add(Sel,0,WACV)

        if 'E' not in Inst['Type'][0]:
            if 'N' in Inst['Type'][0]:
                littleSizer.Add((5,-1))
                MagSel = wx.CheckBox(G2frame.dataWindow,label="Magnetic?")
                unitSizerWidgetList.append(MagSel)
                MagSel.SetValue('MagSpGrp' in SGData)
                MagSel.Bind(wx.EVT_CHECKBOX,OnMagSel)
                littleSizer.Add(MagSel,0,WACV)
            if len(G2frame.HKL) and 'PKS' not in G2frame.GPXtree.GetItemText(G2frame.PatternId):
                littleSizer.Add((5,-1))
                makePks = wx.Button(G2frame.dataWindow,label='Make new Peak list')
                unitSizerWidgetList.append(makePks)
                makePks.Bind(wx.EVT_BUTTON,OnMakePks)
                littleSizer.Add(makePks,0,WACV)
        unitSizer.Add(littleSizer,0)

        # magnetic cell options
        if 'N' in Inst['Type'][0] and 'MagSpGrp' in SGData:
            neutSizer = wx.BoxSizer(wx.HORIZONTAL)
            GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
            SGData['GenSym'] = GenSym
            SGData['SGGray'] = False
            neutSizer.Add(wx.StaticText(G2frame.dataWindow,label=' BNS lattice: '),0,WACV)
            BNSkeys = [SGData['SGLatt'],]+list(BNSsym.keys())
            BNSkeys.sort()
            try:        #this is an ugly kluge - bug in wx.ComboBox
                if SGData['BNSlattsym'][0][2] in ['a','b','c']:
                    BNSkeys.reverse()
            except:
                pass
            BNS = wx.ComboBox(G2frame.dataWindow,value=SGData['BNSlattsym'][0],
                choices=BNSkeys,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            BNS.Bind(wx.EVT_COMBOBOX,OnBNSlatt)
            neutSizer.Add(BNS,0,WACV)
            spinColor = ['black','red']
            spCode = {-1:'red',1:'black'}
            for isym,sym in enumerate(GenSym[1:]):
                neutSizer.Add(wx.StaticText(G2frame.dataWindow,label=' %s: '%(sym.strip())),0,WACV)
                spinOp = wx.ComboBox(G2frame.dataWindow,value=spCode[SGData['SGSpin'][isym+1]],choices=spinColor,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[spinOp.GetId()] = isym
                spinOp.Bind(wx.EVT_COMBOBOX,OnSpinOp)
                neutSizer.Add(spinOp,0,WACV)
            OprNames,SpnFlp = G2spc.GenMagOps(SGData)
            SGData['SpnFlp'] = SpnFlp
            showSpins = wx.Button(G2frame.dataWindow,label=' Show spins?')
            showSpins.Bind(wx.EVT_BUTTON,OnShowSpins)
            neutSizer.Add(showSpins,0,WACV)
            unitSizer.Add(neutSizer,0)
            unitSizer.Add((5,5),0)
        return unitSizer

# Space group grid
    def SpGrpGrid():

        def OnSelectSgrp(event):
            'Called when the Space Group Search Results show column is checked'
            if event is not None:
                clearShowFlags()
                r = event.GetRow()
                for i in range(len(ssopt['SgResults'])):
                    ssopt['SgResults'][i][1] = False
                    SgTable.SetValue(i,1,False)
                SgTable.SetValue(r,1,True)
                controls[13] = ssopt['SgResults'][r][0]
            else:
                for r in range(len(ssopt['SgResults'])):
                    if ssopt['SgResults'][r][1]:
                        controls[13] = ssopt['SgResults'][r][0]
                        break
            SgDisplay.ForceRefresh()
            ssopt['SGData'] = G2spc.SpcGroup(controls[13])[1]
            G2frame.dataWindow.RefineCell.Enable(True)
            OnHklShow(event,True,indexFrom=' Space group selection %s #%d'%(controls[13],r))

        SpGrpGrid = wx.BoxSizer(wx.VERTICAL)
        lbl = (' Space group search results from "Try all"'+
                   '\n  '+ssopt.get('SgSettings',''))
        SpGrpGrid.Add(wx.StaticText(parent=G2frame.dataWindow,label=lbl))
        colLabels = ['Sp Grp','show','M20','X20','Nhkl','fr. found']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,
            wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_FLOAT+':10,3']
        rowLabels = []
        table = []
        for result in ssopt['SgResults']:
            rowLabels.append('')
            row = result
            table.append(row)
        SgTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        SgDisplay = G2G.GSGrid(G2frame.dataWindow)
        SgDisplay.SetTable(SgTable, True)
        SgDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,OnSelectSgrp)
        SgDisplay.SetRowLabelSize(0)
        SgDisplay.AutoSizeColumns(False)
        for r in range(SgDisplay.GetNumberRows()):
            for c in range(SgDisplay.GetNumberCols()):
                if c == 1:
                    SgDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    SgDisplay.SetReadOnly(r,c,isReadOnly=True)
        SpGrpGrid.Add(SgDisplay)
        _setupResultsGrid(SgDisplay)
        OnSelectSgrp(None)
        return SpGrpGrid

    def SSGrid():

        def OnSelectSSG(event):
            'Called when the Super Space Group Search Results show column is checked'
            if event is not None:
                clearShowFlags()
                r = event.GetRow()
                for i in range(len(ssopt['SSgResults'])):
                    ssopt['SSgResults'][i][1] = False
                    SSgTable.SetValue(i,1,False)
                SSgTable.SetValue(r,1,True)
#                controls[13] = ssopt['SgResults'][r][0]
            else:
                for r in range(len(ssopt['SSgResults'])):
                    if ssopt['SSgResults'][r][1]:
#                        controls[13] = ssopt['SSgResults'][r][0]
                        break
            SSDisplay.ForceRefresh()
#            ssopt['SGData'] = G2spc.SpcGroup(controls[13])[1]
            # G2frame.dataWindow.RefineCell.Enable(True)
            OnHklShow(event,True,indexFrom=' Super Space group selection %s #%d'%(controls[13],r))

        SSGrpGrid = wx.BoxSizer(wx.VERTICAL)
        lbl = (' Super Space group search results from "Try all"'+
                   '\n  '+ssopt.get('SgSettings',''))
        SSGrpGrid.Add(wx.StaticText(parent=G2frame.dataWindow,label=lbl))
        colLabels = ['SSp Grp','show','M20','X20','Nhkl','fr. found']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,
            wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_FLOAT+':10,3']
        rowLabels = []
        table = []
        for result in ssopt['SSgResults']:
            rowLabels.append('')
            row = result
            table.append(row)
        SSgTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        SSDisplay = G2G.GSGrid(G2frame.dataWindow)
        SSDisplay.SetTable(SSgTable, True)
        SSDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,OnSelectSSG)
        SSDisplay.SetRowLabelSize(0)
        SSDisplay.AutoSizeColumns(False)
        for r in range(SSDisplay.GetNumberRows()):
            for c in range(SSDisplay.GetNumberCols()):
                if c == 1:
                    SSDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    SSDisplay.SetReadOnly(r,c,isReadOnly=True)
        SSGrpGrid.Add(SSDisplay)
        _setupResultsGrid(SSDisplay)
        OnSelectSSG(None)
        return SSGrpGrid

    def MagSubGrid():

        global KeyList
        def ClearCurrentShowNext():
            KeepShowNext(False)

        KeyList += [['j',ClearCurrentShowNext,'Show next Mag. Spc. Group, clear keep flag on current']]

        def KeepCurrentShowNext():
            KeepShowNext(True)

        KeyList += [['k',KeepCurrentShowNext,'Show next Mag. Spc. Group, keep current']]

        def KeepShowNext(KeepCurrent=True):
            '''Show next "keep" item in Magnetic Space Group list, possibly resetting the
            keep flag for the current displayed cell
            '''
            for i in range(len(magcells)): # find plotted setting
                if magcells[i]['Use']: break
            else:
                return # no Try is set
            if not KeepCurrent:  # clear current
                magcells[i]['Keep'] = False
                MagCellsTable.SetValue(i,2,False)
            keeps = [j for j in range(i+1,len(magcells)) if magcells[j]['Keep']]
            if not keeps:
                if not KeepCurrent: magDisplay.ForceRefresh()
                return # no remaining Keep-flagged entries
            next = keeps[0]
            # update table
            magcells[i]['Use'] = False
            MagCellsTable.SetValue(i,1,False)
            magcells[next]['Use'] = True
            MagCellsTable.SetValue(next,1,True)
            # get SG info and plot
            SGData = magcells[next]['SGData']
            A = G2lat.cell2A(magcells[next]['Cell'][:6])
            G2frame.HKL = G2pwd.getHKLpeak(1.0,SGData,A,Inst)
            G2pwpl.PlotPatterns(G2frame,extraKeys=KeyList)
            magDisplay.ForceRefresh()
            # change Scroll to display new setting
            xscroll = G2frame.dataWindow.GetScrollPos(wx.HORIZONTAL)
            yscroll = magDisplay.CellToRect(next,1)[1]/G2frame.dataWindow.GetScrollPixelsPerUnit()[1]
            G2frame.dataWindow.Scroll(xscroll,yscroll)

        def RefreshMagCellsGrid(event):
            'Display results from k-SUBGROUPSMAG in the Unit Cells tab & allow inspection of results'
            controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(UnitCellsId)
            r,c =  event.GetRow(),event.GetCol()
            rLab = magDisplay.GetRowLabelValue(r)
            br = baseList[r]
            phase = phaseDict[br]
            pname = '(%s) %s'%(rLab,phase['Name'])
            if magcells:
                if c == 0:
                    mSGData = phase['SGData']
                    text,table = G2spc.SGPrint(mSGData,AddInv=True)
                    if 'magAtms' in phase:
                        msg = 'Magnetic space group information'
                        text[0] = ' Magnetic Space Group: '+mSGData['MagSpGrp']
                        text[3] = ' The magnetic lattice point group is '+mSGData['MagPtGp']
                        OprNames,SpnFlp = G2spc.GenMagOps(mSGData)
                        G2G.SGMagSpinBox(G2frame.dataWindow,msg,text,table,mSGData['SGCen'],OprNames,
                            mSGData['SpnFlp'],False).Show()
                    else:
                        msg = 'Space Group Information'
                        G2G.SGMessageBox(G2frame.dataWindow,msg,text,table).Show()
                elif c == 1:
                    for i in range(len(magcells)):
                        magcells[i]['Use'] = False
                    for i in range(len(baseList)):
                        MagCellsTable.SetValue(i,c,False)
                    MagCellsTable.SetValue(r,c,True)
                    magDisplay.ForceRefresh()
                    phase['Use'] = True
                    mSGData = phase['SGData']
                    A = G2lat.cell2A(phase['Cell'][:6])
                    G2frame.HKL = np.array(G2pwd.getHKLpeak(1.0,mSGData,A,Inst))
                    G2pwpl.PlotPatterns(G2frame,extraKeys=KeyList)
                elif c == 2:
                    if MagCellsTable.GetValue(r,c):
                        MagCellsTable.SetValue(r,c,False)
                        phase['Keep'] = False
                    else:
                        phase['Keep'] = True
                        MagCellsTable.SetValue(r,c,True)
                    magDisplay.ForceRefresh()
                elif c ==3:
                    maxequiv = magcells[0].get('maxequiv',100)
                    mSGData = phase['SGData']
                    Uvec = phase['Uvec']
                    Trans = phase['Trans']
                    ifMag = False
                    if 'magAtms' in phase:
                        ifMag = True
                        allmom = phase.get('allmom',False)
                        magAtms = phase.get('magAtms','')
                        mAtoms = TestMagAtoms(phase,magAtms,SGData,Uvec,Trans,allmom,maxequiv)
                    else:
                        mAtoms = TestAtoms(phase,controls[15],SGData,Uvec,Trans,maxequiv)
                    Atms = []
                    AtCods = []
                    atMxyz = []
                    for ia,atom in enumerate(mAtoms):
                        atom[0] += '_%d'%ia
                        SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[2:5],mSGData)
                        Atms.append(atom[:2]+['',]+atom[2:5])
                        AtCods.append('1')
                        if 'magAtms' in phase:
                            MagSytSym = G2spc.MagSytSym(SytSym,dupDir,mSGData)
                            CSI = G2spc.GetCSpqinel(mSGData['SpnFlp'],dupDir)
                            atMxyz.append([MagSytSym,CSI[0]])
                        else:
                            CSI = G2spc.GetCSxinel(SytSym)
                            atMxyz.append([SytSym,CSI[0]])
                    G2phsG.UseMagAtomDialog(G2frame,pname,Atms,AtCods,atMxyz,ifMag=ifMag,ifOK=True).ShowModal()
                elif c in [4,5]:
                    if 'altList' not in phase: return
                    if c == 4:
                        title = 'Conjugacy list for '+pname
                        items = phase['altList']
                        ifPick = True
                    elif c == 5:
                        title = 'Super groups list for '+pname
                        items = phase['supList']
                        if not items[0]:
                            wx.MessageBox(pname+' is a maximal subgroup',caption='Super group is parent',style=wx.ICON_INFORMATION)
                            return
                        ifPick = False
                    Pick = -1
                    dlg = SubCellsDialog(G2frame,title,controls,SGData,items,phaseDict,ifPick)
                    try:
                        if dlg.ShowModal() == wx.ID_OK:
                            altList = copy.copy(phase['altList'])
                            Pick = dlg.GetSelection()
                            pickPhase = phaseDict[altList[Pick]]
                            pickPhase['altList'] = altList
                            pickPhase['Use'] = phase['Use']
                            baseList[r] = altList[Pick]
                    finally:
                        dlg.Destroy()
                    if Pick >= 0:
                        data = [controls,bravais,cells,dminx,ssopt,magcells]
                        G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
                        wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        def OnRefreshKeep(event):
            controls,bravais,cells,dminx,ssopt,magcells = G2frame.GPXtree.GetItemPyData(UnitCellsId)
            c =  event.GetCol()
            E,SGData = G2spc.SpcGroup(controls[13])
            if c == 2:
                testAtoms = ['',]+list(set([atom[1] for atom in controls[15]]))
                ifMag = False
                maxequiv = magcells[0]['maxequiv']
                maximal = False
                if 'magAtms' in magcells[0]:
                    ifMag = True
                    allmom = magcells[0]['allmom']
                    magAtms = magcells[0]['magAtms']
                    dlg = G2G.MultiDataDialog(G2frame,title='Keep options',
                        prompts=['max unique','test for mag. atoms','all have moment','only maximal subgroups',],
                        values=[maxequiv,'',allmom,False],limits=[[1,100],testAtoms,[True,False],[True,False]],
                        formats=['%d','choice','bool','bool'])
                else:
                    dlg = G2G.MultiDataDialog(G2frame,title='Keep options',
                        prompts=['max unique','only maximal subgroups',],
                        values=[maxequiv,False],limits=[[1,100],[True,False],],
                        formats=['%d','bool',])
                if dlg.ShowModal() == wx.ID_OK:
                    if ifMag:
                        maxequiv,atype,allmom,maximal = dlg.GetValues()
                        magAtms = [atom for atom in controls[15] if atom[1] == atype]
                    else:
                        maxequiv,maximal = dlg.GetValues()
                dlg = wx.ProgressDialog('Setting Keep flags','Processing '+magcells[0]['Name'],len(magcells),
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)
                for ip,phase in enumerate(magcells):
                    dlg.Update(ip,newmsg='Processing '+phase['Name'])
                    Uvec = phase['Uvec']
                    Trans = phase['Trans']
                    if ifMag:
                        phase['nAtoms'] = len(TestMagAtoms(phase,magAtms,SGData,Uvec,Trans,allmom,maxequiv,maximal))
                    else:
                        phase['nAtoms'] = len(TestAtoms(phase,controls[15],SGData,Uvec,Trans,maxequiv,maximal))
                dlg.Destroy()
                data = controls,bravais,cells,dminx,ssopt,magcells
                G2frame.GPXtree.SetItemPyData(UnitCellsId,data)
                wx.CallAfter(UpdateUnitCellsGrid,G2frame,data)

        MagSubGrid = wx.BoxSizer(wx.VERTICAL)
        itemList = [phase.get('gid',ip+1) for ip,phase in enumerate(magcells)]
        phaseDict = dict(zip(itemList,magcells))
        G2frame.dataWindow.CopyCell.Enable(False)
        kvec1 = ','.join(controls[14][:3])
        kvec2 = ','.join(controls[14][3:6])
        kvec3 = ','.join(controls[14][6:])
        baseList = controls[16]
        if 'magAtms' in magcells[0]:
            G2frame.dataWindow.RunSubGroupsMag.Enable(True)
            Label = '\n Magnetic subgroup cells from Bilbao k-SUBGROUPSMAG for %s; kvec1=(%s)'%(controls[13],kvec1)
        else:
            G2frame.dataWindow.RunSubGroups.Enable(True)
            Label = '\n Subgroup cells from Bilbao SUBGROUPS for %s; kvec1=(%s)'%(controls[13],kvec1)
        if ' ' not in kvec2:
            Label += ', kvec2=(%s)' % kvec2
        if ' ' not in kvec3:
            Label += ', kvec3=(%s)' % kvec3
        Label += ':'
        MagSubGrid.Add(wx.StaticText(parent=G2frame.dataWindow,label=Label))
        rowLabels = [str(i+1) for i in range(len(baseList))]
        colLabels = ['Space Gp','Try','Keep','Uniq','nConj','nSup','Trans','Vec','a','b','c','\u03B1','\u03B2','\u03B3','Volume']
        Types = [wg.GRID_VALUE_STRING,]+2*[wg.GRID_VALUE_BOOL,]+3*[wg.GRID_VALUE_LONG,]+2*[wg.GRID_VALUE_STRING,]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+3*[wg.GRID_VALUE_FLOAT+':10,3',]+[wg.GRID_VALUE_FLOAT+':10,2']
        table = []
        for ip in baseList:
            phase = phaseDict[ip]
            natms = phase.get('nAtoms',1)
            try:
                nConj = len(phase['altList'])
                nSup = len(phase['supList'])
            except KeyError:
                nConj = 0
                nSup = 0
            cell  = list(phase['Cell'])
            trans = G2spc.Trans2Text(phase['Trans'])
            vec = G2spc.Latt2text([phase['Uvec'],])
            row = [phase['Name'],phase['Use'],phase['Keep'],natms,nConj,nSup,trans,vec]+cell
            table.append(row)
        MagCellsTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        G2frame.GetStatusBar().SetStatusText(
                'Double click Keep to refresh Keep flags; click Space Gp to see sym. ops., Uniq to see unique atoms list; Try to trigger K & J keys on plot',1)
        magDisplay = G2G.GSGrid(G2frame.dataWindow)
        magDisplay.SetRowLabelSize(45)
        magDisplay.SetTable(MagCellsTable, True)
        magDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,RefreshMagCellsGrid)
        magDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnRefreshKeep)
        magDisplay.AutoSizeColumns(False)
        for r in range(magDisplay.GetNumberRows()):
            for c in range(magDisplay.GetNumberCols()):
                if c in [1,2]:
                    magDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    magDisplay.SetReadOnly(r,c,isReadOnly=True)
        MagSubGrid.Add(magDisplay)
        _setupResultsGrid(magDisplay)
        return MagSubGrid

    #### UpdateUnitCellsGrid code starts here
    # create all menubars accessed here
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.LimitMenu)   # Needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PeakMenu)   # Needed below
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.IndexMenu)
    resizeGrids = []     # track Grids to resize
    Indx = {}
    G2frame.ifSetLimitsMode = 0
    G2frame.CancelSetLimitsMode.Enable(False)
    UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Unit Cells List')
    SPGlist = G2spc.spglist
    spaceGroups = ['F m 3 m','I m 3 m','P m 3 m','R 3 m','P 6/m m m','I 4/m m m',
        'P 4/m m m','F m m m','I m m m','A m m m','B m m m','C m m m','P m m m','I 2/m','A 2/m','C 2/m','P 2/m','P -1','C -1']
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))[1]
    if 'T' in Inst['Type'][0]:
        #difC = Inst['difC'][1]
        dmin = G2lat.Pos2dsp(Inst,Limits[0])
    elif 'E' in Inst['Type'][0]:
        #TTh = Inst['2-theta'][1]
        dmin = G2lat.Pos2dsp(Inst,Limits[1])
    else:   #'C', 'B', or 'PKS'
        #wave = G2mth.getWave(Inst)
        dmin = G2lat.Pos2dsp(Inst,Limits[1])
    G2frame.GetStatusBar().SetStatusText('')
    G2frame.Bind(wx.EVT_MENU, OnIndexPeaks, id=G2G.wxID_INDEXPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnRunSubs, id=G2G.wxID_RUNSUB)
    G2frame.Bind(wx.EVT_MENU, OnRunSubsMag, id=G2G.wxID_RUNSUBMAG)
    G2frame.Bind(wx.EVT_MENU, OnLatSym, id=G2G.wxID_LATSYM)
    G2frame.Bind(wx.EVT_MENU, OnNISTLatSym, id=G2G.wxID_NISTLATSYM)
    G2frame.Bind(wx.EVT_MENU, CopyUnitCell, id=G2G.wxID_COPYCELL)
    G2frame.Bind(wx.EVT_MENU, LoadUnitCell, id=G2G.wxID_LOADCELL)
    #G2frame.Bind(wx.EVT_MENU, ImportUnitCell, id=G2G.wxID_IMPORTCELL)
    G2frame.Bind(wx.EVT_MENU, TransformUnitCell, id=G2G.wxID_TRANSFORMCELL)
    G2frame.Bind(wx.EVT_MENU, onRefineCell, id=G2G.wxID_REFINECELL)
    G2frame.Bind(wx.EVT_MENU, MakeNewPhase, id=G2G.wxID_MAKENEWPHASE)
    G2frame.Bind(wx.EVT_MENU, OnExportCells, id=G2G.wxID_EXPORTCELLS)
    G2frame.Bind(wx.EVT_MENU, OnShowGenRefls, id=G2G.wxID_SHOWGENHKLS)
    G2frame.Bind(wx.EVT_MENU, OnClearCells, id=G2G.wxID_CLEARCELLS)
    if len(data) < 6:
        data.append([])
    controls,bravais,cells,dminx,ssopt,magcells = data
    if len(controls) < 13:              #add cell volume if missing
        controls.append(G2lat.calc_V(G2lat.cell2A(controls[6:12])))
    if len(controls) < 14:              #add space group if missing
        controls.append(spaceGroups[bravaisSymb.index(controls[5])])
    if len(controls) < 15:
        controls.append(list(range(1,len(magcells)+1)))
    while len(bravais) < 18:
        bravais += [0,]
    SGData = ssopt.get('SGData',G2spc.SpcGroup(controls[13])[1])
    G2frame.GPXtree.SetItemPyData(UnitCellsId,data)            #update with volume
    bravaisNames = ['Cubic-F','Cubic-I','Cubic-P','Trigonal-R','Trigonal/Hexagonal-P',
        'Tetragonal-I','Tetragonal-P','Orthorhombic-F','Orthorhombic-I','Orthorhombic-A',
        'Orthorhombic-B','Orthorhombic-C','Orthorhombic-P',
        'Monoclinic-I','Monoclinic-A','Monoclinic-C','Monoclinic-P','Triclinic','Triclinic',]

    G2frame.dataWindow.IndexPeaks.Enable(False)
    peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
    if peaks:
        G2frame.dataWindow.IndexPeaks.Enable(True)
    G2frame.dataWindow.RefineCell.Enable(False)
    if controls[12] > 1.0 and len(peaks[0]):             #if a "real" volume (i.e. not default) and peaks
        G2frame.dataWindow.RefineCell.Enable(True)
    G2frame.dataWindow.CopyCell.Enable(False)
    G2frame.dataWindow.MakeNewPhase.Enable(False)
    G2frame.dataWindow.ExportCells.Enable(False)
    if cells:
        G2frame.dataWindow.CopyCell.Enable(True)
        G2frame.dataWindow.MakeNewPhase.Enable(True)
        G2frame.dataWindow.ExportCells.Enable(True)
    elif magcells:
        G2frame.dataWindow.CopyCell.Enable(True)
    if G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Phases'):
        G2frame.dataWindow.LoadCell.Enable(True)
    pkId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List')
    if pkId:
        peakList = G2frame.GPXtree.GetItemPyData(pkId)
        # in case we are loading this without visiting the Peak List first, initialize
        peakList['xtraMode'] = peakList.get('xtraMode',False)
        G2frame.dataWindow.XtraPeakMode.Check(peakList['xtraMode'])

    # GUI code
    setRow = None
    G2frame.dataWindow.ClearData()
    # setup for resizing
    G2frame.dataWindow.customResize = _onResizeUnitCellsList
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Indexing tools'),0,WACV)
    if not hasattr(G2frame,'kvecSearch'):
        G2frame.kvecSearch = {'mode':False}

    if G2frame.dataWindow.XtraPeakMode.IsChecked():
        cb = G2G.G2CheckBox(parent,'Search for k-vector',G2frame.kvecSearch,'mode',OnChange=updateCellsWindow)
        topSizer.Add(cb,0,WACV|wx.LEFT,15)
    else:
        G2frame.kvecSearch['mode'] = False
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))

    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add((-1,3),0)
    if not G2frame.kvecSearch['mode']:
        # autoindexing GUI
        mainSizer.Add(firstSizer())
    else:
        # k-vector GUI
        mainSizer.Add(kvecSizer())

    # 2nd "box": unit cell/sym info
    mainSizer.Add((-1,3),0)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add((-1, 3), 0)
    unitSizerWidgetList = []   # entries in unitSizer to disable/enable
    mainSizer.Add(unitSizer())
    # 3rd "box" search results
    mainSizer.Add((-1,3),0)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add(wx.StaticText(parent=G2frame.dataWindow,
        label='Cell Search Results',style=wx.ALIGN_CENTER),0,wx.EXPAND,3)
    G2frame.dataWindow.currentGrids = []

    # space group search results
    if len(ssopt.get('SgResults',[])):
         mainSizer.Add(SpGrpGrid())

    # cell search results
    if cells:
        mode = 0
        try: # for Cell sym, 1st entry is cell xform matrix;
             # for kvector search, 1st item is a question mark
            len(cells[0][0])
            mode = 1
            if cells[0][0] == '?': mode = 2
        except:
            pass
    # k-vector search results table
        if mode == 2:
            G2frame.kvecSearch['mode'] == True
            colLabels = ['show']
            Types = [wg.GRID_VALUE_BOOL]
            colLabels += [
                'kx', 'ky', 'kz',
                'Ave. ' + u'\u03B4' + 'd/d',
                'Ave. ' + u'\u03B4' + 'd',
                'Max. ' + u'\u03B4' + 'd'
            ]
            Types += (6 * [wg.GRID_VALUE_FLOAT + ':10, 5'])
            mainSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label='\n k-vector search results:'))
        elif mode == 1:
            mainSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label='\n Cell symmetry search:'))
            colLabels = ['show']
            Types = [wg.GRID_VALUE_BOOL]
            colLabels += ['a','b','c','\u03B1','\u03B2','\u03B3','Volume','Keep']
            Types += (3*[wg.GRID_VALUE_FLOAT+':10,5',]+
                  3*[wg.GRID_VALUE_FLOAT+':10,3',]+
                  [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL])
        else:
    # indexing search results table
            mainSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label='\n Indexing Result:'))
            colLabels = ['M20','X20','show','Bravais']
            Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_NUMBER,
                         wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING]
            colLabels += ['a','b','c','\u03B1','\u03B2','\u03B3','Volume','Keep']
            Types += (3*[wg.GRID_VALUE_FLOAT+':10,5',]+
                    3*[wg.GRID_VALUE_FLOAT+':10,3',]+
                    [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_BOOL])
        rowLabels = []
        table = []
        # are we going to keep a current "Show" flag?
        isUseSet = any([cell[-2] for cell in cells]) and showUse and mode == 0
        setRow = None
        for i,cell in enumerate(cells):
            # reset all "show" flags when table is first created, unless
            # we are going to set the Show flag in a row
            if not isUseSet:
                cell[-2] = False
            elif cell[-2]:
                setRow = i
            rowLabels.append('')
            if mode:
                row = [cell[-2]]+cell[3:10]+[cell[11],]
            else:
                row = cell[0:2]+[cell[-2]]+[bravaisSymb[cell[2]]]+cell[3:10]+[cell[11],]
            # if cell[-2]:
            #     if mode != 2:
            #         A = G2lat.cell2A(cell[3:9])
            #         G2frame.HKL = G2lat.GenHBravais(dmin,cell[2],A)
            #         for hkl in G2frame.HKL:
            #             hkl.insert(4,G2lat.Dsp2pos(Inst,hkl[3])+controls[1])
            #         G2frame.HKL = np.array(G2frame.HKL)
            #     else:
            #         # We need to fill in the todos when the mode is 2, i.e.,
            #         # the k-vector search.
            #         pass
            table.append(row)
        UnitCellsTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        gridDisplay = G2G.GSGrid(G2frame.dataWindow)
        gridDisplay.SetTable(UnitCellsTable, True)
        G2frame.dataWindow.CopyCell.Enable(True)
        gridDisplay.Bind(wg.EVT_GRID_CELL_LEFT_CLICK,SeaResSelected)
        gridDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK,OnSortCells)
        gridDisplay.SetRowLabelSize(0)
        gridDisplay.AutoSizeColumns(False)
        for r in range(gridDisplay.GetNumberRows()):
            for c in range(gridDisplay.GetNumberCols()):
                if c == 2:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=False)
                else:
                    gridDisplay.SetReadOnly(r,c,isReadOnly=True)
        if mode == 2:
            #OnISODISTORT_kvec(phase_sel) # TODO: not ready yet

            hSizer = wx.BoxSizer(wx.HORIZONTAL)
            hSizer.Add(gridDisplay)

            # TODO: Add a button to call ISODISTORT
            # ISObut = wx.Button(G2frame.dataWindow,label='Call ISODISTORT')
            # ISObut.Bind(wx.EVT_BUTTON,OnISODIST)
            # hSizer.Add(ISObut)

            mainSizer.Add(hSizer)
        else:
            mainSizer.Add(gridDisplay)
        _setupResultsGrid(gridDisplay)
        if setRow:   # we have a "Show" flag to set
            setGrid = gridDisplay
            gridDisplay.SetCellValue(setRow,2,'True')   # set flag in grid
            disableCellEntries()

    # Subgroup/magnetic s.g. search results
    if magcells and len(controls) > 16:
        mainSizer.Add(MagSubGrid())

    # GUI creation done -- finally
    G2frame.Contour = False
    if New:
       OnHklShow(None,indexFrom=' Indexing from unit cell & symmetry settings')
    G2frame.dataWindow.SetDataSize()
    if callSeaResSelected:
        SeaResSelected(None)  # select 1st item in table
    elif setRow:
        wx.CallAfter(setGrid.MakeCellVisible,setRow,0)  # scroll to line setRow

################################################################################
#####  Reflection list
################################################################################
def UpdateReflectionGrid(G2frame,data,HKLF=False,Name=''):
    '''respond to selection of PWDR or HKLF Reflections data tree
    item by displaying a table of reflections in the data window.

    Note that this is used for Single Xtal even though in pwdGUI.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    dMin = 0.05
    if not isinstance(G2frame.Hide, bool):
        G2frame.Hide = False
    if 'UsrReject' in Controls:
        dMin = Controls['UsrReject'].get('MinD',0.05)

    def OnPlot1DHKL(event):
        phaseName = G2frame.RefList
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetGPXtreeItemId(G2frame,pId,phaseName)
            General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
        else:
            Super = 0
        if 'list' in str(type(data)):   #single crystal data is 2 dict in list
            refList = data[1]['RefList']
        else:                           #powder data is a dict of dicts; each same structure as SC 2nd dict
            if 'RefList' in data[phaseName]:
                refList = np.array(data[phaseName]['RefList'])
            else:
                wx.MessageBox('No reflection list - do Refine first',caption='Reflection plotting')
                return
        G2plt.Plot1DSngl(G2frame,newPlot=True,hklRef=refList,Super=Super,Title=phaseName)

    def OnPlotHKL(event):
        '''Plots a layer of reflections
        '''
        phaseName = G2frame.RefList
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetGPXtreeItemId(G2frame,pId,phaseName)
            General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
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
        Super = 0
        SuperVec = []
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetGPXtreeItemId(G2frame,pId,phaseName)
            General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
            if General.get('Modulated',False):
                Super = 1
                SuperVec = General['SuperVec']
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

    def OnWilsonStat(event):
        ''' Show Wilson plot for PWDR and HKLF & return Wilson statistics <<E>, <E^2> & <E^2-1> to console
        '''
        phaseName = G2frame.RefList
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  G2gd.GetGPXtreeItemId(G2frame,pId,phaseName)
            General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
        else:
            Super = 0
        if 'list' in str(type(data)):   #single crystal data is 2 dict in list
            refList = data[1]['RefList']
        else:                           #powder data is a dict of dicts; each same structure as SC 2nd dict
            if 'RefList' in data[phaseName]:
                refList = np.array(data[phaseName]['RefList'])
            else:
                wx.MessageBox('No reflection list - do Refine first',caption='Reflection plotting')
                return

        PE = G2elemGUI.PickElement(G2frame,ifNone=False)
        if PE.ShowModal() == wx.ID_OK:
            normEle = PE.Elem.strip()
        PE.Destroy()
        Estat,Ehist = G2mth.DoWilsonStat(refList,Super,normEle,Inst)
        print(' Wilson statistics   : <|E|>: %.3f, <|E^2|>: %.3f, <|E^2-1|>: %.3f'%(Estat[0]/Estat[1],1.,Estat[2]/Estat[1]))
        print(' Expected: random P-1: <|E|>: 0.798, <|E^2|>: 1.000, <|E^2-1|>: 0.968')
        print('           random P1 : <|E|>: 0.886, <|E^2|>: 1.000, <|E^2-1|>: 0.736')
        XY = [[Ehist[0],Ehist[1]],[Ehist[0],Ehist[2]]]
        G2plt.PlotXY(G2frame,XY,labelX='sin$^2$%s/%s$^2$'%(GkTheta,Gklambda),
            labelY=r'ln(<|F$_o$|$^2$>/%sf$^2$)'%GkSigma,newPlot=True,Title='Wilson plot')

    def OnMakeCSV(event):
        '''Make csv file from displayed ref table.
        '''
        phaseName = G2frame.RefList
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose Reflection List csv file', pth, '',
            'Reflection List (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.csv'
                File = open(filename,'w')
                File.write('%s\n'%phaseName)
                colLabels = [G2frame.PeakTable.GetColLabelValue(i) for i in range(G2frame.PeakTable.GetNumberCols())]
                File.write('%s\n'%(','.join(colLabels)))
                nRows = G2frame.PeakTable.GetNumberRows()
                for i in range(nRows):
                    refLine = G2frame.PeakTable.GetRowValues(i)
                    strLine = ','.join([str(item) for item in refLine])
                    File.write('%s\n'%strLine)
                File.close()
        finally:
            dlg.Destroy()

    def MakeReflectionTable(phaseName):
        '''Returns a wx.grid table (G2G.Table) containing a list of all reflections
        for a phase.
        '''
        Super = 0
        if phaseName not in ['Unknown',]:
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            if pId: # phase section missing from file (unusual)
                phaseId =  G2gd.GetGPXtreeItemId(G2frame,pId,phaseName)
                if phaseId:         #is phase deleted?
                    General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
                    G,g = G2lat.cell2Gmat(General['Cell'][1:7])
                    GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
                    SGData = General['SGData']
                    if General.get('Modulated',False):
                        Super = 1
                    try:
                        Histograms = G2frame.GPXtree.GetItemPyData(phaseId)['Histograms']
                        histName = G2frame.GPXtree.GetItemText(G2frame.PatternId)
                        histData = Histograms[histName]
                    except:
                        if GSASIIpath.GetConfigValue('debug'):
                            print('Reflection table problem: histogram {} not found in phase {}'.format(histName,phaseName))
                        return
        rowLabels = []
        if HKLF:
            if G2frame.Hide:
                refList = np.array([refl for refl in data[1]['RefList'] if refl[3]])
            else:
                refList = data[1]['RefList']
            refs = refList
        else:
            muStrData = histData['Mustrain']
            sizeData = histData['Size']
            if len(data) > 1:
                G2frame.dataWindow.SelectPhase.Enable(True)
            try:            #patch for old reflection lists
                if not len(data[phaseName]):
                    return None
                refList = np.array(data[phaseName]['RefList'])
                I100 = refList.T[8+Super]*refList.T[11+Super]
            except TypeError:
                refList = np.array([refl[:11+Super] for refl in data[phaseName]])
                I100 = refList.T[8+Super]*np.array([refl[11+Super] for refl in data[phaseName]])
            MuStr = G2pwd.getMustrain(refList.T[:3],G,SGData,muStrData)
            CrSize = G2pwd.getCrSize(refList.T[:3],G,GB,sizeData)
            Imax = np.max(I100)
            if Imax:
                I100 *= 100.0/Imax
            if 'C' in Inst['Type'][0]:
                refs = np.vstack((refList.T[:15+Super],I100,MuStr,CrSize)).T
            elif 'T' in Inst['Type'][0]:
                refs = np.vstack((refList.T[:18+Super],I100,MuStr,CrSize)).T
            elif Inst['Type'][0][2] in ['A','B']:
                refs = np.vstack((refList.T[:17+Super],I100,MuStr,CrSize)).T
            elif 'E' in Inst['Type'][0]:
                refs = np.vstack((refList.T[:12+Super],I100,MuStr,CrSize)).T        #last two not shown for now
        rowLabels = [str(i) for i in range(len(refs))]
        Types = (4+Super)*[wg.GRID_VALUE_LONG,]+4*[wg.GRID_VALUE_FLOAT+':10,4',]+ \
            2*[wg.GRID_VALUE_FLOAT+':10,2',]+[wg.GRID_VALUE_FLOAT+':10,3',]+ \
            [wg.GRID_VALUE_FLOAT+':10,3',]
        if HKLF:
            colLabels = ['H','K','L','flag','d','Fosq','sig','Fcsq','FoTsq','FcTsq','phase','ExtC',]
            if 'T' in Inst['Type'][0]:
                colLabels = ['H','K','L','flag','d','Fosq','sig','Fcsq','FoTsq','FcTsq','phase','ExtC','wave','tbar']
                Types += 2*[wg.GRID_VALUE_FLOAT+':10,3',]
            if Super:
                colLabels.insert(3,'M')
        else:
            if 'C' in Inst['Type'][0]:
                colLabels = ['H','K','L','mul','d','pos','sig\u00b2','gam','Fo\u00b2','Fc\u00b2','phase','Icorr','Prfo','Trans','ExtP','I100','\u03bcstrain','Size']
                Types += 6*[wg.GRID_VALUE_FLOAT+':10,3',]
            elif 'T' in Inst['Type'][0]:
                colLabels = ['H','K','L','mul','d','pos','sig\u00b2','gam','Fo\u00b2','Fc\u00b2','phase','Icorr','alp','bet','wave','Prfo','Abs','Ext','I100','\u03bcstrain','Size']
                Types += 9*[wg.GRID_VALUE_FLOAT+':10,3',]
            elif Inst['Type'][0][2] in ['A','B']:
                colLabels = ['H','K','L','mul','d','pos','sig\u00b2','gam','Fo\u00b2','Fc\u00b2','phase','Icorr','alp','bet','Prfo','Abs','Ext','I100','\u03bcstrain','Size']
                Types += 8*[wg.GRID_VALUE_FLOAT+':10,3',]
            elif 'E' in Inst['Type'][0]:
                colLabels = ['H','K','L','mul','d','pos','sig\u00b2','gam','Fo\u00b2','Fc\u00b2','phase','Icorr','I100']
                Types += [wg.GRID_VALUE_FLOAT+':10,3',]
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
        if HKLF:
            G2frame.GetStatusBar().SetStatusText('abs(DF)/sig > 10 red; > 3 yellow; flag:>0 twin no., 0 sp.gp absent, -1 user rejected, -2 Rfree',1)
        else:
            G2frame.GetStatusBar().SetStatusText('Prfo < 0. in red; if excluded Fosq in red & mul < 0',1)
        itof = 0
        if HKLF:
            im = data[1].get('Super',0)
        else:
            if 'T' in data[phaseName].get('Type',''):
                itof = 3
            im = data[phaseName].get('Super',0)
        # has this table already been displayed?
        if G2frame.refTable[phaseName].GetTable() is None:
            G2frame.PeakTable = MakeReflectionTable(phaseName)
            if not G2frame.PeakTable: return
            G2frame.refTable[phaseName].SetTable(G2frame.PeakTable, True)
            G2frame.refTable[phaseName].EnableEditing(False)
            G2frame.refTable[phaseName].SetMargins(0,0)
            G2frame.refTable[phaseName].AutoSizeColumns(False)
            setBackgroundColors(im,itof)
        if HKLF:
            if G2frame.Hide:
                refList = np.array([refl[:6+im] for refl in data[1]['RefList'] if refl[3]])
            else:
                refList = np.array([refl[:6+im] for refl in data[1]['RefList']])
        else:
            refList = np.array([refl[:6+im] for refl in data[phaseName]['RefList']])
        G2frame.HKL = np.vstack((refList.T)).T    #build for plots
        # raise the tab (needed for 1st use and from OnSelectPhase)
        for PageNum in range(G2frame.refBook.GetPageCount()):
            if phaseName == G2frame.refBook.GetPageText(PageNum):
                G2frame.refBook.SetSelection(PageNum)
                break
        else:
            print (phaseName)
            print (phases)
            raise Exception("how did we not find a phase name?")

    def OnToggleExt(event):
        G2frame.Hide = not G2frame.Hide
        UpdateReflectionGrid(G2frame,data,HKLF=True,Name=Name)

    def OnPageChanged(event):
        '''Respond to a press on a phase tab by displaying the reflections. This
        routine is needed because the reflection table may not have been created yet.
        '''
        G2frame.refBook.SetSize(G2frame.dataWindow.GetClientSize())    #TODO -almost right
        page = event.GetSelection()
        phaseName = G2frame.refBook.GetPageText(page)
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

    # start of UpdateReflectionGrid
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.ReflMenu)
    if not data:
        print ('No phases, no reflections')
        lbl = 'No phases, no reflections'
    elif HKLF:
        lbl = 'Single crystal reflections'
    else:
        lbl = 'Powder reflections, selected by phase'
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,wx.ID_ANY,lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    if not data:
        return
    if HKLF:
        G2frame.RefList = 1
        phaseName = IsHistogramInAnyPhase(G2frame,Name)
        if not phaseName:
            phaseName = 'Unknown'
        phases = [phaseName]
    else:
        phaseName = G2frame.RefList
        phases = list(data.keys())
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    if HKLF:
        G2frame.Bind(wx.EVT_MENU, OnPlotHKL, id=G2G.wxID_PWDHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot1DHKL, id=G2G.wxID_1DHKLSTICKPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=G2G.wxID_PWD3DHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnWilsonStat, id=G2G.wxID_WILSONSTAT)
        G2frame.Bind(wx.EVT_MENU, OnToggleExt, id=G2G.wxID_SHOWHIDEEXTINCT)
        G2frame.dataWindow.SelectPhase.Enable(False)
    else:
        G2frame.Bind(wx.EVT_MENU, OnSelectPhase, id=G2G.wxID_SELECTPHASE)
        G2frame.Bind(wx.EVT_MENU, OnPlot1DHKL, id=G2G.wxID_1DHKLSTICKPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlotHKL, id=G2G.wxID_PWDHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=G2G.wxID_PWD3DHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnMakeCSV, id=G2G.wxID_CSVFROMTABLE)
        G2frame.Bind(wx.EVT_MENU, OnWilsonStat, id=G2G.wxID_WILSONSTAT)
        G2frame.dataWindow.SelectPhase.Enable(False)

    G2frame.refBook = G2G.GSNoteBook(parent=G2frame.dataWindow)
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(G2frame.refBook,1,wx.ALL|wx.EXPAND,1)
    G2frame.refTable = {}
    G2frame.dataWindow.currentGrids = []
    for tabnum,phase in enumerate(phases):
        if isinstance(data,list):           #single crystal HKLF
            G2frame.refTable[phase] = G2G.GSGrid(parent=G2frame.refBook)
            G2frame.refTable[phase].SetRowLabelSize(60)  # leave room for big numbers
            G2frame.refBook.AddPage(G2frame.refTable[phase],phase)
            G2frame.refTable[phase].SetScrollRate(10,10) # reflection grids (inside tab) need scroll bars
        elif len(data[phase]):              #else dict for PWDR
            G2frame.refTable[phase] = G2G.GSGrid(parent=G2frame.refBook)
            G2frame.refTable[phase].SetRowLabelSize(60)
            G2frame.refBook.AddPage(G2frame.refTable[phase],phase)
            G2frame.refTable[phase].SetScrollRate(10,10) # as above
        else:       #cleanup deleted phase reflection lists
            del data[phase]
            if len(data):
                G2frame.RefList = list(data.keys())[0]
                phaseName = G2frame.RefList
            else:
                G2frame.RefList = ''
                phaseName = ''
    if phaseName: ShowReflTable(phaseName)
    G2frame.refBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    G2frame.dataWindow.SetDataSize()

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
                isotopes = list(Info['Isotopes'].keys())
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

        names = list(substFile.Substances.keys())
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
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy substances from\n'+hst[5:]+' to...',
            'Copy substances', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Instrument Parameters'))[0]
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
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Substances'),ndata)

    def OnAddSubstance(event):
        dlg = wx.TextEntryDialog(None,'Enter a name for this substance','Substance Name Entry','New substance',
            style=wx.OK|wx.CANCEL)
        if dlg.ShowModal() == wx.ID_OK:
            Name = dlg.GetValue()
            data['Substances'][Name] = {'Elements':{},'Volume':1.0,'Density':1.0,
                'Scatt density':0.0,'XAnom density':0.,'XAbsorption':0.,'XImag density':0.}
            AddElement(Name)
        else:
            return
        dlg.Destroy()
        if not data['Substances'][Name]['XAbsorption']:
            del data['Substances'][Name]
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
        ElList = list(data['Substances'][name]['Elements'].keys())
        dlg = G2elemGUI.PickElements(G2frame,ElList)
        if dlg.ShowModal() == wx.ID_OK:
            for El in dlg.Elem:
                El = El.strip().capitalize()
                Info = G2elem.GetAtomInfo(El)
                Info.update({'Num':1.})
                data['Substances'][name]['Elements'][El] = Info
                isotopes = list(Info['Isotopes'].keys())
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
        else:
            return
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
        ElList = list(data['Substances'][name]['Elements'].keys())
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
        substSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' Substance list: wavelength: %.5fA'%(wave)))
        for name in data['Substances']:
            G2G.HorizontalLine(substSizer,G2frame.dataWindow)
            substSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' Data for '+name+':'),0)
            if name == 'vacuum':
                substSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label='        Not applicable'),0)
            elif name == 'unit scatter':
                substSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Scattering density,f: %.3f *10%scm%s'%(data['Substances'][name]['Scatt density'],Pwr10,Pwrm2)),0)
            else:
                elSizer = wx.FlexGridSizer(0,8,5,5)
                Substance = data['Substances'][name]
                Elems = Substance['Elements']
                for El in Elems:    #do elements as pull downs for isotopes for neutrons
                    elSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' '+El+': '),
                        0,WACV)
                    num = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Substances'][name]['Elements'][El],'Num',
                        nDig=(10,2,'f'),typeHint=float,OnLeave=OnNum)
                    Indx[num.GetId()] = [name,El,'Num']
                    elSizer.Add(num,0,WACV)
                    if 'N' in Inst['Type'][0]:
                        elSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Isotope: '),0,WACV)
                        isotopes = list(Elems[El]['Isotopes'].keys())
                        isotope = wx.ComboBox(G2frame.dataWindow,choices=isotopes,value=Elems[El].get('Isotope','Nat. Abund.'),
                            style=wx.CB_READONLY|wx.CB_DROPDOWN)
                        Indx[isotope.GetId()] = [El,name]
                        isotope.Bind(wx.EVT_COMBOBOX,OnIsotope)
                        elSizer.Add(isotope,0,WACV)
                substSizer.Add(elSizer,0)
                vdsSizer = wx.FlexGridSizer(0,4,5,5)
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' Volume: '),
                    0,WACV)
                vol = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Substances'][name],'Volume',nDig=(10,2),typeHint=float,OnLeave=OnVolDen)
                Indx[vol.GetId()] = [name,'Volume']
                vdsSizer.Add(vol,0,WACV)
                vdsSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' Density: '),
                    0,WACV)
                den = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Substances'][name],'Density',nDig=(10,2),typeHint=float,OnLeave=OnVolDen)
                Indx[den.GetId()] = [name,'Density']
                vdsSizer.Add(den,0,WACV)
                substSizer.Add(vdsSizer,0)
                denSizer = wx.FlexGridSizer(0,2,0,0)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Scattering density,f'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=': %.3f *10%scm%s'%(Substance['Scatt density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=" Real density,f+f'"),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=': %.3f *10%scm%s'%(Substance['XAnom density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Imaginary density,f"'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=': %.3g *10%scm%s'%(Substance['XImag density'],Pwr10,Pwrm2)),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Absorption'),0,WACV)
                denSizer.Add(wx.StaticText(G2frame.dataWindow,label=': %.3g cm%s'%(Substance['XAbsorption'],Pwrm1)),0,WACV)
                substSizer.Add(denSizer)
        return substSizer

    # start of UpdateSubstanceGrid
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SubstanceMenu)
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    Name = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[2]
    wave = G2mth.getWave(Inst)
    G2frame.Bind(wx.EVT_MENU, OnLoadSubstance, id=G2G.wxID_LOADSUBSTANCE)
    G2frame.Bind(wx.EVT_MENU, OnReloadSubstances, id=G2G.wxID_RELOADSUBSTANCES)
    G2frame.Bind(wx.EVT_MENU, OnAddSubstance, id=G2G.wxID_ADDSUBSTANCE)
    G2frame.Bind(wx.EVT_MENU, OnCopySubstance, id=G2G.wxID_COPYSUBSTANCE)
    G2frame.Bind(wx.EVT_MENU, OnDeleteSubstance, id=G2G.wxID_DELETESUBSTANCE)
    G2frame.Bind(wx.EVT_MENU, OnAddElement, id=G2G.wxID_ELEMENTADD)
    G2frame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2G.wxID_ELEMENTDELETE)
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Sample substances for %s:'%Name),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    #print('Substance ',G2frame.dataWindow.helpKey)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(SubstSizer(),0)
    G2frame.dataWindow.SetDataSize()

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
    if 'Pair' not in data:
        data['Pair'] = {'Method':'Moore','MaxRadius':100.,'NBins':100,'Errors':'User','Result':[],
            'Percent error':2.5,'Background':[0,False],'Distribution':[],'Moore':10,'Dist G':100.,}
    if 'Shapes' not in data:
        data['Shapes'] = {'outName':'run','NumAA':100,'Niter':1,'AAscale':1.0,'Symm':1,'bias-z':0.0,
            'inflateV':1.0,'AAglue':0.0,'pdbOut':False,'boxStep':4.0}
    if 'boxStep' not in data['Shapes']:
        data['Shapes']['boxStep'] = 4.0
    plotDefaults = {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':150.,'viewDir':[0,0,1],}

    #end patches

    def RefreshPlots(newPlot=False):
        PlotText = G2frame.G2plotNB.nb.GetPageText(G2frame.G2plotNB.nb.GetSelection())
        if 'Powder' in PlotText:
            G2pwpl.PlotPatterns(G2frame,plotType='SASD',newPlot=newPlot)
        elif 'Size' in PlotText:
            G2plt.PlotSASDSizeDist(G2frame)
        elif 'Pair' in PlotText:
            G2plt.PlotSASDPairDist(G2frame)


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
                'Unified':{'G':[1.e3,False],'Rg':[100.,False],'B':[1.e-5,False],'P':[4.,False],'Cutoff':[1.e-5,False],},
                'Porod':{'B':[1.e-4,False],'P':[4.,False],'Cutoff':[1.e-5,False],},
                'Monodisperse':{'Volume':[0.05,False],'Radius':[100.,False],},   #OK for spheres
                'Bragg':{'PkInt':[100.,False],'PkPos':[0.2,False],
                    'PkSig':[10.,False],'PkGam':[10.,False],},        #reasonable 31A peak
                })
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)

        wx.CallAfter(UpdateModelsGrid,G2frame,data)

    def OnCopyModel(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        wtFactor = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[0]['wtFactor']
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy models from\n'+hst[5:]+' to...',
            'Copy models', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            ProfDict = G2frame.GPXtree.GetItemPyData(Id)[0]
            ProfDict['wtFactor'] = wtFactor
            newdata = copy.deepcopy(data)
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Models'),newdata)
            if newdata['BackFile']:
                Profile = G2frame.GPXtree.GetItemPyData(Id)[1]
                BackId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,newdata['BackFile'])
                BackSample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,BackId, 'Sample Parameters'))
                Profile[5] = BackSample['Scale'][0]*G2frame.GPXtree.GetItemPyData(BackId)[1][1]
        UpdateModelsGrid(G2frame,newdata)
        wx.CallAfter(UpdateModelsGrid,G2frame,data)
        RefreshPlots(True)

    def OnCopyFlags(event):
        thisModel = copy.deepcopy(data)
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy sample ref. flags from\n'+str(hst[5:])+' to...',
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
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    newModel = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Models'))
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
        choices = G2gd.GetGPXtreeDataNames(G2frame,['SASD',])
        od = {'label_1':'Copy to next','value_1':False,'label_2':'Reverse order','value_2':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Sequential SASD refinement',
             'Select dataset to include',choices,extraOpts=od)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(choices[sel])
            Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential SASD results')
            if Id:
                SeqResult = G2frame.GPXtree.GetItemPyData(Id)
            else:
                SeqResult = {}
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential SASD results')
            SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
            SeqResult['histNames'] = names
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        dlg = wx.ProgressDialog('SASD Sequential fit','Data set name = '+names[0],len(names),
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        wx.BeginBusyCursor()
        if od['value_2']:
            names.reverse()
        JModel = None
        try:
            for i,name in enumerate(names):
                print (' Sequential fit for '+name)
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    break
                sId =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                if i and od['value_1']:
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'),JModel)
                IProfDict,IProfile = G2frame.GPXtree.GetItemPyData(sId)[:2]
                IModel = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'))
                ISample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Sample Parameters'))
                ILimits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Limits'))
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
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'),copy.deepcopy(IModel))

                G2sasd.ModelFxn(IProfile,IProfDict,ILimits,ISample,IModel)
                saveDict = {}
                for item in parmDict:
                    if ';' in item:
                        trm = item.split(';')[1]
                        if trm in ['FFVolume','StrFact','FormFact']:
                            continue
                    saveDict.update({item:parmDict[item]})
                SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                    'covMatrix':covMatrix,'title':name,'parmDict':saveDict}
            else:
                dlg.Destroy()
                print (' ***** Small angle sequential refinement successful *****')
        finally:
            wx.EndBusyCursor()
        G2frame.GPXtree.SetItemPyData(Id,SeqResult)
        G2frame.GPXtree.SelectItem(Id)

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

        elif data['Current'] == 'Pair distance':
            SaveState()
            G2sasd.PairDistFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)
            G2plt.PlotSASDPairDist(G2frame)
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

        elif data['Current'] == 'Shapes':
            SaveState()
            wx.MessageBox(' For use of SHAPES, please cite:\n\n'+
                G2G.GetCite('SHAPES'),
                caption='Program Shapes',style=wx.ICON_INFORMATION)
            dlg = wx.ProgressDialog('Running SHAPES','Cycle no.: 0 of 160',161,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)

            data['Pair']['Result'] = []       #clear old results (if any) for now
            data['Pair']['Result'] = G2shapes.G2shapes(Profile,ProfDict,Limits,data,dlg)
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

    def OnUnDo(event):
        DoUnDo()
        data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Models'))
        G2frame.dataWindow.SasdUndo.Enable(False)
        UpdateModelsGrid(G2frame,data)
        G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
        RefreshPlots(True)

    def DoUnDo():
        print ('Undo last refinement')
        file = open(G2frame.undosasd,'rb')
        PatternId = G2frame.PatternId
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Models'),pickle.load(file))
        print (' Models recovered')
        file.close()

    def SaveState():
        G2frame.undosasd = os.path.join(G2frame.dirname,'GSASIIsasd.save')
        file = open(G2frame.undosasd,'wb')
        PatternId = G2frame.PatternId
        for item in ['Models']:
            pickle.dump(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,item)),file,1)
        file.close()
        G2frame.dataWindow.SasdUndo.Enable(True)

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
        sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Size distribution parameters: '),0)
        binSizer = wx.FlexGridSizer(0,7,5,5)
        binSizer.Add(wx.StaticText(G2frame.dataWindow,label=' No. size bins: '),0,WACV)
        bins = ['50','100','150','200']
        nbins = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['Nbins']),choices=bins,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[nbins.GetId()] = [data['Size'],'Nbins',0]
        nbins.Bind(wx.EVT_COMBOBOX,OnIntVal)
        binSizer.Add(nbins,0,WACV)
        binSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Min diam.: '),0,WACV)
        minDias = ['10','25','50','100','150','200']
        mindiam = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['MinDiam']),choices=minDias,
            style=wx.CB_DROPDOWN)
        mindiam.Bind(wx.EVT_LEAVE_WINDOW,OnIntVal)
        mindiam.Bind(wx.EVT_TEXT_ENTER,OnIntVal)
        mindiam.Bind(wx.EVT_KILL_FOCUS,OnIntVal)
        Indx[mindiam.GetId()] = [data['Size'],'MinDiam',0]
        binSizer.Add(mindiam,0,WACV)
        binSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Max diam.: '),0,WACV)
        maxDias = [str(1000*(i+1)) for i in range(10)]
        maxdiam = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['MaxDiam']),choices=maxDias,
            style=wx.CB_DROPDOWN)
        maxdiam.Bind(wx.EVT_LEAVE_WINDOW,OnIntVal)
        maxdiam.Bind(wx.EVT_TEXT_ENTER,OnIntVal)
        maxdiam.Bind(wx.EVT_KILL_FOCUS,OnIntVal)
        Indx[maxdiam.GetId()] = [data['Size'],'MaxDiam',0]
        binSizer.Add(maxdiam,0,WACV)
        logbins = wx.CheckBox(G2frame.dataWindow,label='Log bins?')
        Indx[logbins.GetId()] = [data['Size'],'logBins']
        logbins.SetValue(data['Size']['logBins'])
        logbins.Bind(wx.EVT_CHECKBOX, OnCheckBox)
        binSizer.Add(logbins,0,WACV)
        sizeSizer.Add(binSizer,0)
        sizeSizer.Add((5,5),0)
        partSizer = wx.BoxSizer(wx.HORIZONTAL)
        partSizer.Add(wx.StaticText(G2frame.dataWindow,label='Particle description: '),0,WACV)
        shapes = {'Spheroid':' Aspect ratio: ','Cylinder':' Diameter ','Cylinder AR':' Aspect ratio: ',
            'Unified sphere':'','Unified rod':' Diameter: ','Unified rod AR':' Aspect ratio: ',
            'Unified disk':' Thickness: ', 'Spherical shell': ' Shell thickness'}
        partsh = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['Shape'][0]),choices=list(shapes.keys()),
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        partsh.Bind(wx.EVT_COMBOBOX,OnShape)
        partSizer.Add(partsh,0,WACV)
        if data['Size']['Shape'][0] not in ['Unified sphere',]:
            partSizer.Add(wx.StaticText(G2frame.dataWindow,label=shapes[data['Size']['Shape'][0]]),0,WACV)
            partprm = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Size']['Shape'],1,
                nDig=(10,3),typeHint=float,xmin=0.)
            partSizer.Add(partprm,0,WACV)
        sizeSizer.Add(partSizer,0)
        sizeSizer.Add((5,5),0)
        fitSizer = wx.BoxSizer(wx.HORIZONTAL)
        methods = ['MaxEnt','IPG',]
        fitSizer.Add(wx.StaticText(G2frame.dataWindow,label='Fitting method: '),0,WACV)
        method = wx.ComboBox(G2frame.dataWindow,value=data['Size']['Method'],choices=methods,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        method.Bind(wx.EVT_COMBOBOX,OnMethod)
        fitSizer.Add(method,0,WACV)
        iters = ['10','25','50','100','150','200']
        fitSizer.Add(wx.StaticText(G2frame.dataWindow,label=' No. iterations: '),0,WACV)
        Method = data['Size']['Method']
        iter = wx.ComboBox(G2frame.dataWindow,value=str(data['Size'][Method]['Niter']),choices=iters,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[iter.GetId()] = [data['Size'][Method],'Niter',0]
        iter.Bind(wx.EVT_COMBOBOX,OnIntVal)
        fitSizer.Add(iter,0,WACV)
        if 'MaxEnt' in data['Size']['Method']:
            fitSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Log floor factor: '),0,WACV)
            floors = [str(-i) for i in range(9)]
            floor = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['MaxEnt']['Sky']),choices=floors,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[floor.GetId()] = [data['Size']['MaxEnt'],'Sky',-10]
            floor.Bind(wx.EVT_COMBOBOX,OnIntVal)
            fitSizer.Add(floor,0,WACV)
        elif 'IPG' in data['Size']['Method']:
            fitSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Q power weight (-1 for sigma): '),0,WACV)
            choices = ['-1','0','1','2','3','4']
            power = wx.ComboBox(G2frame.dataWindow,value=str(data['Size']['IPG']['Power']),choices=choices,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[power.GetId()] = [data['Size']['IPG'],'Power',-2]
            power.Bind(wx.EVT_COMBOBOX,OnIntVal)
            fitSizer.Add(power,0,WACV)
        sizeSizer.Add(fitSizer,0)

        return sizeSizer

    def PairSizer():

        def OnMethod(event):
            data['Pair']['Method'] = method.GetValue()
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

        def OnError(event):
            data['Pair']['Errors'] = error.GetValue()
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

        def OnMaxRadEst(event):
            Results = G2sasd.RgFit(Profile,ProfDict,Limits,Sample,data)
            if not Results[0]:
                    G2frame.ErrorDialog('Failed refinement',
                        ' Msg: '+Results[-1]+'\nYou need to rethink your selection of parameters\n'+    \
                        ' Model restored to previous version')
            RefreshPlots(True)
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

        def OnMooreTerms(event):
            data['Pair']['Moore'] = int(round(Limits[1][1]*data['Pair']['MaxRadius']/np.pi))-1
            wx.CallAfter(UpdateModelsGrid,G2frame,data)

        def OnNewVal(invalid,value,tc):
            if invalid: return
            parmDict = {'Rg':data['Pair']['MaxRadius']/2.5,'G':data['Pair']['Dist G'],
                'B':data['Pair'].get('Dist B',Profile[1][-1]*Profile[0][-1]**4),
                'Back':data['Back'][0]}
            Profile[2] = G2sasd.getSASDRg(Profile[0],parmDict)
            RefreshPlots(True)

        pairSizer = wx.BoxSizer(wx.VERTICAL)
        pairSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Pair distribution parameters: '),0)
        binSizer = wx.FlexGridSizer(0,6,5,5)
        binSizer.Add(wx.StaticText(G2frame.dataWindow,label=' No. R bins: '),0,WACV)
        bins = ['50','100','150','200']
        nbins = wx.ComboBox(G2frame.dataWindow,value=str(data['Pair']['NBins']),choices=bins,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[nbins.GetId()] = [data['Pair'],'NBins',0]
        nbins.Bind(wx.EVT_COMBOBOX,OnIntVal)
        binSizer.Add(nbins,0,WACV)
        binSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Max diam.: '),0,WACV)
        maxdiam = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Pair'],'MaxRadius',xmin=10.,nDig=(10,1),OnLeave=OnNewVal)
        binSizer.Add(maxdiam,0,WACV)
        maxest = wx.Button(G2frame.dataWindow,label='Make estimate')
        maxest.Bind(wx.EVT_BUTTON,OnMaxRadEst)
        binSizer.Add(maxest,0,WACV)
        pairSizer.Add(binSizer,0)
        pairSizer.Add((5,5),0)
        fitSizer = wx.BoxSizer(wx.HORIZONTAL)
        methods = ['Moore',]            #'Regularization',
        fitSizer.Add(wx.StaticText(G2frame.dataWindow,label='Fitting method: '),0,WACV)
        method = wx.ComboBox(G2frame.dataWindow,value=data['Pair']['Method'],choices=methods,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        method.Bind(wx.EVT_COMBOBOX,OnMethod)
        fitSizer.Add(method,0,WACV)
        if data['Pair']['Method'] == 'Moore':
            fitSizer.Add(wx.StaticText(G2frame.dataWindow,label=" P.B. Moore, J. Appl. Cryst., 13, 168-175 (1980)"),0,WACV)
        else:
            fitSizer.Add(wx.StaticText(G2frame.dataWindow,label=" D.I. Svergun, J. Appl. Cryst., 24, 485-492 (1991)"),0,WACV)
        pairSizer.Add(fitSizer,0)
        if 'Moore' in data['Pair']['Method']:
            mooreSizer = wx.BoxSizer(wx.HORIZONTAL)
            mooreSizer.Add(wx.StaticText(G2frame.dataWindow,label='Number of functions: '),0,WACV)
            moore = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Pair'],'Moore',xmin=2,xmax=20)
            mooreSizer.Add(moore,0,WACV)
            mooreterms = wx.Button(G2frame.dataWindow,label = 'Auto determine?')
            mooreterms.Bind(wx.EVT_BUTTON,OnMooreTerms)
            mooreSizer.Add(mooreterms,0,WACV)
            pairSizer.Add(mooreSizer,0)
        errorSizer = wx.BoxSizer(wx.HORIZONTAL)
        errorSizer.Add(wx.StaticText(G2frame.dataWindow,label='Error method: '),0,WACV)
        errors = ['User','Sqrt','Percent']
        error = wx.ComboBox(G2frame.dataWindow,value=data['Pair']['Errors'],choices=errors,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        error.Bind(wx.EVT_COMBOBOX,OnError)
        if 'Percent' in data['Pair']['Errors']:
            percent = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Pair'],'Percent error',xmin=0.5,nDig=(10,1))
            errorSizer.Add(percent,0,WACV)
        errorSizer.Add(error,0,WACV)
        pairSizer.Add(errorSizer,0)
        return pairSizer

    def ShapesSizer():

#        def OnPDBout(event):
#            data['Shapes']['pdbOut'] = not data['Shapes']['pdbOut']

        def OnShapeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            for i in [1,2]:
                for j in range(len(Patterns)):
                    shapeTable.SetValue(j,i,False)
            shapeTable.SetValue(r,c,True)
            ShapesResult.ForceRefresh()
            Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))[1]
            ProfDict,Profile = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[:2]
            iBeg = np.searchsorted(Profile[0],Limits[0])
            iFin = np.searchsorted(Profile[0],Limits[1])
            pattern = Patterns[r]
            Profile[3][iBeg:iFin+1] = np.array(pattern[2])
            selAtoms = Atoms[2*r+(c-1)]
            prCalc = PRcalc[r][2]
            prDelt= np.diff(PRcalc[r][0])[0]
            prsum = np.sum(prCalc)
            prCalc /= prsum*prDelt
            data['Pair']['Pair Calc'] = np.array([PRcalc[r][0],prCalc]).T
            print('%s %d'%('num. beads',len(selAtoms[1])))
            print('%s %.3f'%('selected r value',pattern[-1]))
            print('%s %.3f'%('selected Delta P(r)',PRcalc[r][-1]))
            PDBtext = 'P(R) dif: %.3f r-value: %.3f Nbeads: %d'%(PRcalc[r][-1],pattern[-1],len(selAtoms[1]))
#            RefreshPlots(True)
            G2pwpl.PlotPatterns(G2frame,plotType='SASD',newPlot=True)
            G2plt.PlotSASDPairDist(G2frame)
            G2plt.PlotBeadModel(G2frame,selAtoms,plotDefaults,PDBtext)

        shapeSizer = wx.BoxSizer(wx.VERTICAL)
        shapeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Shape parameters:'),0)
        parmSizer = wx.FlexGridSizer(0,4,5,5)
#1st row
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' No. amino acids: '),0,WACV)
        numAA = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'NumAA',xmin=10)
        parmSizer.Add(numAA,0,WACV)
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Nballs=no. amino acids*'),0,WACV)
        scaleAA = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'AAscale',xmin=0.01,xmax=10.,nDig=(10,2))
        parmSizer.Add(scaleAA,0,WACV)
#2nd row
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Inflate by (1.-1.4): '),0,WACV)
        inflate = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'inflateV',xmin=1.,xmax=1.4,nDig=(10,2))
        parmSizer.Add(inflate,0,WACV)
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Axial symmetry (1-12): '),0,WACV)
        symm = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'Symm',xmin=1,xmax=12)
        parmSizer.Add(symm,0,WACV)
#3rd row
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' z-axis bias (-2 to 2): '),0,WACV)
        zaxis = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'bias-z',xmin=-2.,xmax=2.,nDig=(10,2))
        parmSizer.Add(zaxis,0,WACV)
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' elongation (0-20): '),0,WACV)
        glue = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'AAglue',xmin=0.,xmax=20.,nDig=(10,2))
        parmSizer.Add(glue,0,WACV)
#4th row
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' No. iterations (1-10): '),0,WACV)
        niter = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'Niter',xmin=1,xmax=10)
        parmSizer.Add(niter,0,WACV)
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Output name: '),0,WACV)
        name = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'outName')
        parmSizer.Add(name,0,WACV)
#last row
        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Bead separation (3.5-5): '),0,WACV)
        beadsep = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Shapes'],'boxStep',xmin=3.5,xmax=5,nDig=(10,1))
        parmSizer.Add(beadsep,0,WACV)
#        pdb = wx.CheckBox(G2frame.dataWindow,label=' Save as pdb files?: ')
#        pdb.SetValue(data['Shapes']['pdbOut'])
#        pdb.Bind(wx.EVT_CHECKBOX, OnPDBout)
#        parmSizer.Add(pdb,0,WACV)

        shapeSizer.Add(parmSizer)

        if len(data['Pair'].get('Result',[])):
            shapeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' SHAPES run results:'),0)
            Atoms,Patterns,PRcalc = data['Pair']['Result']
            colLabels = ['name','show beads','show shape','Rvalue','P(r) dif','Nbeads','Nshape']
            Types = [wg.GRID_VALUE_STRING,]+2*[wg.GRID_VALUE_BOOL,]+2*[wg.GRID_VALUE_FLOAT+':10,3',]+2*[wg.GRID_VALUE_LONG,]
            rowLabels = [str(i) for i in range(len(Patterns))]
            tableVals = []
            for i in range(len(Patterns)):
                tableVals.append([Atoms[2*i][0],False,False,Patterns[i][-1],PRcalc[i][-1],len(Atoms[2*i][1]),len(Atoms[2*i+1][1])])
            shapeTable = G2G.Table(tableVals,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            ShapesResult = G2G.GSGrid(G2frame.dataWindow)
            ShapesResult.SetRowLabelSize(45)
            ShapesResult.SetTable(shapeTable,True)
            ShapesResult.AutoSizeColumns(False)
            ShapesResult.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, OnShapeSelect)
            for r in range(len(Patterns)):
                for c in range(7):
                    if c in [1,2]:
                        ShapesResult.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        ShapesResult.SetReadOnly(r,c,isReadOnly=True)

            shapeSizer.Add(ShapesResult,0)
        return shapeSizer

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
                sldrObj.SetScaledRange(slMult*valMinMax[0],slMult*valMinMax[1])
                sldrObj.SetValue(slMult*logv)
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)

        def OnSelect(event):
            Obj = event.GetEventObject()
            item,key = Indx[Obj.GetId()]
            if key in ['NumPoints',]:
                item[key] = int(Obj.GetValue())
            else:
                item[key] = Obj.GetValue()
            if 'Refine' not in Obj.GetLabel():
                if 'FormFact' in key :
                    item['FFargs'] = FormFactors[Obj.GetValue()]
                elif 'StrFact' in key:
                    item['SFargs'] = StructureFactors[Obj.GetValue()]
                wx.CallAfter(UpdateModelsGrid,G2frame,data)
                G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
                RefreshPlots(True)

        def OnDelLevel(event):
            Obj = event.GetEventObject()
            item = Indx[Obj.GetId()]
            del data['Particle']['Levels'][item]
            wx.CallAfter(UpdateModelsGrid,G2frame,data)
            G2sasd.ModelFxn(Profile,ProfDict,Limits,Sample,data)
            RefreshPlots(True)

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
            RefreshPlots(True)

        def SizeSizer():
            sizeSizer = wx.FlexGridSizer(0,4,5,5)
            sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Distribution: '),0,WACV)
            Distchoice = ['LogNormal','Gaussian','LSW','Schulz-Zimm','Bragg','Unified','Porod','Monodisperse',]
            distChoice = wx.ComboBox(G2frame.dataWindow,value=level['Controls']['DistType'],choices=Distchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[distChoice.GetId()] = [level['Controls'],'DistType']
            distChoice.Bind(wx.EVT_COMBOBOX,OnSelect)
            sizeSizer.Add(distChoice,0,WACV)    #put structure factor choices here
            if level['Controls']['DistType'] not in ['Bragg','Unified','Porod',]:
                sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Form Factor: '),0,WACV)
                if 'Mono' in level['Controls']['DistType']:
                    ffChoice = wx.ComboBox(G2frame.dataWindow,value=level['Controls']['FormFact'],choices=ffMonoChoices,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                else:
                    ffChoice = wx.ComboBox(G2frame.dataWindow,value=level['Controls']['FormFact'],choices=ffDistChoices,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[ffChoice.GetId()] = [level['Controls'],'FormFact']
                ffChoice.Bind(wx.EVT_COMBOBOX,OnSelect)
                sizeSizer.Add(ffChoice,0,WACV)

                sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Material: '),0,WACV)
                matSel = wx.ComboBox(G2frame.dataWindow,value=level['Controls']['Material'],
                    choices=list(Substances['Substances'].keys()),style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[matSel.GetId()] = [level['Controls'],'Material']
                matSel.Bind(wx.EVT_COMBOBOX,OnSelect)
                sizeSizer.Add(matSel,0,WACV) #do neutron test here?
                rho = Substances['Substances'][level['Controls']['Material']].get('XAnom density',0.0)
                level['Controls']['Contrast'] = contrast = (rho-rhoMat)**2
                sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Resonant X-ray contrast: '),0,WACV)
                sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label='  %.2f 10%scm%s'%(contrast,Pwr20,Pwrm4)),0,WACV)
                if 'Mono' not in level['Controls']['DistType']:
                    sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Num. radii: '),0,WACV)
                    radii = ['25','50','75','100','200']
                    nRadii = wx.ComboBox(G2frame.dataWindow,value=str(level['Controls']['NumPoints']),choices=radii,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    Indx[nRadii.GetId()] = [level['Controls'],'NumPoints']
                    nRadii.Bind(wx.EVT_COMBOBOX,OnSelect)
                    sizeSizer.Add(nRadii,0,WACV)
                    sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' R dist. cutoff: '),0,WACV)
                    rCutoff = G2G.ValidatedTxtCtrl(G2frame.dataWindow,level['Controls'],'Cutoff',
                        xmin=0.001,xmax=0.1,typeHint=float)
                    sizeSizer.Add(rCutoff,0,WACV)
            elif level['Controls']['DistType']  in ['Unified',]:
                Parms = level['Unified']
                Best = G2sasd.Bestimate(Parms['G'][0],Parms['Rg'][0],Parms['P'][0])
                sizeSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Estimated Dist B: %12.4g'%(Best)),0,WACV)
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
                        parmSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Dist '+parm,style=wx.ALIGN_CENTER))
                    else:
                        parmVar = wx.CheckBox(G2frame.dataWindow,label='Refine? Dist '+parm)
                        parmVar.SetValue(Parms[parm][1])
                        parmVar.Bind(wx.EVT_CHECKBOX, OnSelect)
                        parmSizer.Add(parmVar,0,WACV)
                        Indx[parmVar.GetId()] = [Parms[parm],1]
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
                    parmValue = wx.TextCtrl(G2frame.dataWindow,value='%.3g'%(Parms[parm][0]),
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
                    parmSldr = G2G.G2Slider(G2frame.dataWindow,minValue=slMult*valMinMax[0],
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
                            parmVar = wx.CheckBox(G2frame.dataWindow,label='Refine? '+fTypes[iarg]+parm)
                            parmVar.SetValue(Args[parm][1])
                            Indx[parmVar.GetId()] = [Args[parm],1]
                            parmVar.Bind(wx.EVT_CHECKBOX, OnSelect)
                            parmSizer.Add(parmVar,0,WACV)
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
                            parmValue = wx.TextCtrl(G2frame.dataWindow,value='%.3g'%(Args[parm][0]),
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
                            parmSldr = G2G.G2Slider(G2frame.dataWindow,minValue=slMult*valMinMax[0],
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
        topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Particle fit parameters: '),0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Matrix: '),0,WACV)
        matsel = wx.ComboBox(G2frame.dataWindow,value=data['Particle']['Matrix']['Name'],
            choices=list(Substances['Substances'].keys()),style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[matsel.GetId()] = [data['Particle']['Matrix'],'Name']
        matsel.Bind(wx.EVT_COMBOBOX,OnSelect) #Do neutron test here?
        rhoMat = Substances['Substances'][data['Particle']['Matrix']['Name']].get('XAnom density',0.0)
        topSizer.Add(matsel,0,WACV)
        topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Volume fraction: '),0,WACV)
        volfrac = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Particle']['Matrix']['VolFrac'],0,typeHint=float)
        topSizer.Add(volfrac,0,WACV)
        volVar = wx.CheckBox(G2frame.dataWindow,label=' Refine?')
        volVar.SetValue(data['Particle']['Matrix']['VolFrac'][1])
        Indx[volVar.GetId()] = [data['Particle']['Matrix']['VolFrac'],1]
        volVar.Bind(wx.EVT_CHECKBOX, OnSelect)
        topSizer.Add(volVar,0,WACV)
        partSizer.Add(topSizer,0,)
        for ilev,level in enumerate(data['Particle']['Levels']):
            G2G.HorizontalLine(partSizer,G2frame.dataWindow)
            topLevel = wx.BoxSizer(wx.HORIZONTAL)
            topLevel.Add(wx.StaticText(G2frame.dataWindow,label=' Model component %d: '%(ilev)),0,WACV)
            delBtn = wx.Button(G2frame.dataWindow,label=' Delete?')
            Indx[delBtn.GetId()] = ilev
            delBtn.Bind(wx.EVT_BUTTON,OnDelLevel)
            topLevel.Add(delBtn,0,WACV)
            partSizer.Add(topLevel,0)
            partSizer.Add(SizeSizer())
            if level['Controls']['DistType'] not in ['Bragg','Unified','Porod',]:
                topLevel.Add(wx.StaticText(G2frame.dataWindow,label=' Structure factor: '),0,WACV)
                strfctr = wx.ComboBox(G2frame.dataWindow,value=level['Controls']['StrFact'],
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
            BackId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['BackFile'])
            BackSample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,BackId, 'Sample Parameters'))
            Profile[5] = BackSample['Scale'][0]*G2frame.GPXtree.GetItemPyData(BackId)[1][1]
        else:
            Profile[5] = np.zeros(len(Profile[5]))
        RefreshPlots(True)

    # start of UpdateModelsGrid
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.ModelMenu)
    Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Sample Parameters'))
    Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))
    Substances = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Substances'))
    ProfDict,Profile = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[:2]
    if data['BackFile']:
        BackId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['BackFile'])
        BackSample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,BackId, 'Sample Parameters'))
        Profile[5] = BackSample['Scale'][0]*G2frame.GPXtree.GetItemPyData(BackId)[1][1]
    G2frame.dataWindow.ClearData()
    G2frame.Bind(wx.EVT_MENU, OnCopyModel, id=G2G.wxID_MODELCOPY)
    G2frame.Bind(wx.EVT_MENU, OnCopyFlags, id=G2G.wxID_MODELCOPYFLAGS)
    G2frame.Bind(wx.EVT_MENU, OnFitModel, id=G2G.wxID_MODELFIT)
    G2frame.Bind(wx.EVT_MENU, OnFitModelAll, id=G2G.wxID_MODELFITALL)
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2G.wxID_MODELUNDO)
    G2frame.Bind(wx.EVT_MENU, OnAddModel, id=G2G.wxID_MODELADD)
    Indx = {}
    G2frame.dataWindow.ClearData()
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    models = ['Size dist.','Particle fit','Pair distance',]
    if len(data['Pair']['Distribution']):
        models += ['Shapes',]
    topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Modeling by: '),0,WACV)
    fitSel = wx.ComboBox(G2frame.dataWindow,value=data['Current'],choices=models,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    fitSel.Bind(wx.EVT_COMBOBOX,OnSelectFit)
    topSizer.Add(fitSel,0,WACV)
    topSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Error multiplier: '),0,WACV)
#        azmthOff = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'azmthOff',nDig=(10,2),typeHint=float,OnLeave=OnAzmthOff)
    esdScale = wx.TextCtrl(G2frame.dataWindow,value='%.3f'%(1./np.sqrt(ProfDict['wtFactor'])),style=wx.TE_PROCESS_ENTER)
    esdScale.Bind(wx.EVT_TEXT_ENTER,OnEsdScale)
    esdScale.Bind(wx.EVT_KILL_FOCUS,OnEsdScale)
    topSizer.Add(esdScale,0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(G2frame.dataWindow,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer.Add(topSizer,0,wx.EXPAND)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    if 'Size' in data['Current']:
        G2frame.dataWindow.SasSeqFit.Enable(False)
        if 'MaxEnt' in data['Size']['Method']:
            G2frame.GetStatusBar().SetStatusText('Size distribution by Maximum entropy',1)
        elif 'IPG' in data['Size']['Method']:
            G2frame.GetStatusBar().SetStatusText('Size distribution by Interior-Point Gradient',1)
        mainSizer.Add(SizeSizer())
    elif 'Particle' in data['Current']:
        G2frame.dataWindow.SasSeqFit.Enable(True)
        mainSizer.Add(PartSizer(),1,wx.EXPAND)
    elif 'Pair' in data['Current']:
        G2frame.dataWindow.SasSeqFit.Enable(False)
        mainSizer.Add(PairSizer(),1,wx.EXPAND)
    elif 'Shape' in data['Current']:
        G2frame.dataWindow.SasSeqFit.Enable(False)
        mainSizer.Add(ShapesSizer(),1,wx.EXPAND)
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    backSizer = wx.BoxSizer(wx.HORIZONTAL)
    backSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Background:'),0,WACV)
    backVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Back'],0,
        nDig=(10,3,'g'),OnLeave=OnBackChange)
    backSizer.Add(backVal,0,WACV)
    if 'Shape' not in data['Current']:
        backVar = wx.CheckBox(G2frame.dataWindow,label='Refine?')
        Indx[backVar.GetId()] = [data['Back'],1]
        backVar.SetValue(data['Back'][1])
        backVar.Bind(wx.EVT_CHECKBOX, OnCheckBox)
        backSizer.Add(backVar,0,WACV)
        #multiple background files?
    backSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Background file: '),0,WACV)
    Choices = ['',]+G2gd.GetGPXtreeDataNames(G2frame,['SASD',])
    backFile = wx.ComboBox(parent=G2frame.dataWindow,value=data['BackFile'],choices=Choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    backFile.Bind(wx.EVT_COMBOBOX,OnBackFile)
    backSizer.Add(backFile)
    mainSizer.Add(backSizer)
    G2frame.dataWindow.SetDataSize()

################################################################################
#####  REFD Models
################################################################################

def UpdateREFDModelsGrid(G2frame,data):
    '''respond to selection of REFD Models data tree item.
    '''
    def OnCopyModel(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        copyList = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy reflectivity models from\n'+str(hst[5:])+' to...',
            'Copy parameters', histList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for i in dlg.GetSelections():
                    copyList.append(histList[i])
        finally:
            dlg.Destroy()
        for item in copyList:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            G2frame.GPXtree.SetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,Id,'Models'),copy.deepcopy(data))

    def OnFitModel(event):

        SaveState()
        G2pwd.REFDRefine(Profile,ProfDict,Inst,Limits,Substances,data)
        x,xr,y = G2pwd.makeSLDprofile(data,Substances)
        ModelPlot(data,x,xr,y)
        G2pwpl.PlotPatterns(G2frame,plotType='REFD')
        wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)

    def OnModelPlot(event):
        hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        histList = GetFileList(G2frame,'REFD')
#        histList = [hst,]
#        histList += GetHistsLikeSelected(G2frame)
        if not histList:
            G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
            return
        plotList = []
        od = {'label_1':'Zero at substrate','value_1':False,'label_2':'Show layer transitions','value_2':True}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Plot reflectivity models for:',
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
            mId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
            model = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,mId,'Models'))
            if len(model['Layers']) == 2:   #see if any real layers defined; will be > 2
                continue
            Substances = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,mId,'Substances'))['Substances']
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
        choices = G2gd.GetGPXtreeDataNames(G2frame,['REFD',])
        od = {'label_1':'Copy to next','value_1':False,'label_2':'Reverse order','value_2':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Sequential REFD refinement',
             'Select dataset to include',choices,extraOpts=od)
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                names.append(choices[sel])
            Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential REFD results')
            if Id:
                SeqResult = G2frame.GPXtree.GetItemPyData(Id)
            else:
                SeqResult = {}
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential REFD results')
            SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
            SeqResult['histNames'] = names
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        dlg = wx.ProgressDialog('REFD Sequential fit','Data set name = '+names[0],len(names),
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
        wx.BeginBusyCursor()
        if od['value_2']:
            names.reverse()
        JModel = None
        try:
            for i,name in enumerate(names):
                print (' Sequential fit for '+name)
                GoOn = dlg.Update(i,newmsg='Data set name = '+name)[0]
                if not GoOn:
                    break
                sId =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                if i and od['value_1']:
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'),JModel)
                IProfDict,IProfile = G2frame.GPXtree.GetItemPyData(sId)[:2]
                IModel = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'))
                ISubstances = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Substances'))['Substances']
                ILimits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Limits'))
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
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,sId, 'Models'),copy.deepcopy(IModel))

                SeqResult[name] = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                    'covMatrix':covMatrix,'title':name,'parmDict':parmDict}
            else:
                dlg.Destroy()
                print (' ***** Small angle sequential refinement successful *****')
        finally:
            wx.EndBusyCursor()
        G2frame.GPXtree.SetItemPyData(Id,SeqResult)
        G2frame.GPXtree.SelectItem(Id)

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
        data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Models'))
        G2frame.dataWindow.REFDUndo.Enable(False)
        G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
        x,xr,y = G2pwd.makeSLDprofile(data,Substances)
        ModelPlot(data,x,xr,y)
        G2pwpl.PlotPatterns(G2frame,plotType='REFD')
        wx.CallLater(100,UpdateREFDModelsGrid,G2frame,data)

    def DoUnDo():
        print ('Undo last refinement')
        file = open(G2frame.undorefd,'rb')
        PatternId = G2frame.PatternId
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Models'),pickle.load(file))
        print (' Model recovered')
        file.close()

    def SaveState():
        G2frame.undorefd = os.path.join(G2frame.dirname,'GSASIIrefd.save')
        file = open(G2frame.undorefd,'wb')
        PatternId = G2frame.PatternId
        pickle.dump(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Models')),file,1)
        file.close()
        G2frame.dataWindow.REFDUndo.Enable(True)

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
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)

#        def OnQ4fftplot(event):
#            q4fft.SetValue(False)
#            R,F = G2pwd.makeRefdFFT(Limits,Profile)
#            XY = [[R[:2500],F[:2500]],]
#            G2plt.PlotXY(G2frame,XY,labelX='thickness',labelY='F(R)',newPlot=True,
#                Title='Fourier transform',lines=True)

        def OndQSel(event):
            data['dQ type'] = dQSel.GetStringSelection()
            Recalculate()

        def NewRes(invalid,value,tc):
            Recalculate()

        def Recalculate():
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')

        controlSizer = wx.BoxSizer(wx.VERTICAL)
        resol = wx.BoxSizer(wx.HORIZONTAL)
        choice = ['None','const '+GkDelta+'Q/Q',]
        if ProfDict['ifDQ']:
            choice += [GkDelta+'Q/Q in data']
        dQSel = wx.RadioBox(G2frame.dataWindow,wx.ID_ANY,'Instrument resolution type:',choices=choice,
            majorDimension=0,style=wx.RA_SPECIFY_COLS)
        dQSel.SetStringSelection(data['dQ type'])
        dQSel.Bind(wx.EVT_RADIOBOX,OndQSel)
        resol.Add(dQSel,0,WACV)
        resol.Add(wx.StaticText(G2frame.dataWindow,label=' (FWHM %): '),0,WACV)
        resol.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Resolution'],0,nDig=(10,3),xmin=0.,xmax=5.,OnLeave=NewRes),0,WACV)
        controlSizer.Add(resol,0)
        minimiz = wx.BoxSizer(wx.HORIZONTAL)
        minimiz.Add(wx.StaticText(G2frame.dataWindow,label=' Minimizer: '),0,WACV)
        minlist = ['LMLS','Basin Hopping','MC/SA Anneal','L-BFGS-B',]
        minSel = wx.ComboBox(G2frame.dataWindow,value=data['Minimizer'],choices=minlist,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        minSel.Bind(wx.EVT_COMBOBOX, OnMinSel)
        minimiz.Add(minSel,0,WACV)
        minimiz.Add(wx.StaticText(G2frame.dataWindow,label=' Bounds factor: '),0,WACV)
        minimiz.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Toler',nDig=(10,2),xmax=0.99,xmin=0.1),0,WACV)
        weight = wx.CheckBox(G2frame.dataWindow,label='Use 2% sig. weights')
        weight.SetValue(data.get('2% weight',False))
        weight.Bind(wx.EVT_CHECKBOX, OnWeight)
        minimiz.Add(weight,0,WACV)
        controlSizer.Add(minimiz,0)
        plotSizer = wx.BoxSizer(wx.HORIZONTAL)
        plotSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Plot controls: '),0,WACV)
        sld = wx.Button(G2frame.dataWindow,label='Plot SLD?')
        sld.Bind(wx.EVT_BUTTON, OnSLDplot)
        plotSizer.Add(sld,0,WACV)
        plotSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Zero position location: '),0,WACV)
        poslist = ['Top','Bottom']
        refpos = wx.ComboBox(G2frame.dataWindow,value=data['Zero'],choices=poslist,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        refpos.Bind(wx.EVT_COMBOBOX, OnRefPos)
        plotSizer.Add(refpos,0,WACV)
#        q4fft = wx.CheckBox(G2frame.dataWindow,label='Plot fft?')
#        q4fft.Bind(wx.EVT_CHECKBOX, OnQ4fftplot)
#        plotSizer.Add(q4fft,0,WACV)
        controlSizer.Add(plotSizer,0)
        return controlSizer

    def OverallSizer():
#'DualFitFile':'', 'DualFltBack':[0.0,False],'DualScale':[1.0,False] future for neutrons - more than one?

        def OnScaleRef(event):
            data['Scale'][1] = scaleref.GetValue()

        def OnBackRef(event):
            data['FltBack'][1] = backref.GetValue()

        def OnSliderMax(event):
            data['slider max'] = float(slidermax.GetValue())
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)

        def Recalculate(invalid,value,tc):
            if invalid:
                return

            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')

        overall = wx.BoxSizer(wx.HORIZONTAL)
        overall.Add(wx.StaticText(G2frame.dataWindow,label=' Scale: '),0,WACV)
        overall.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Scale'],0,
            nDig=(10,2),typeHint=float,OnLeave=Recalculate),0,WACV)
        scaleref = wx.CheckBox(G2frame.dataWindow,label=' Refine?  ')
        scaleref.SetValue(data['Scale'][1])
        scaleref.Bind(wx.EVT_CHECKBOX, OnScaleRef)
        overall.Add(scaleref,0,WACV)
        overall.Add(wx.StaticText(G2frame.dataWindow,label=' Flat bkg.: '),0,WACV)
        overall.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['FltBack'],0,
            nDig=(10,2,'g'),typeHint=float,OnLeave=Recalculate),0,WACV)
        backref = wx.CheckBox(G2frame.dataWindow,label=' Refine?  ')
        backref.SetValue(data['FltBack'][1])
        backref.Bind(wx.EVT_CHECKBOX, OnBackRef)
        overall.Add(backref,0,WACV)
        overall.Add(wx.StaticText(G2frame.dataWindow,label=' Select slider range 0-'),0,WACV)
        choice = ['1000','2000','5000','10000']
        slidermax = wx.ComboBox(G2frame.dataWindow,value='%d'%data['slider max'],choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        slidermax.Bind(wx.EVT_COMBOBOX,OnSliderMax)
        overall.Add(slidermax,0,WACV)
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
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')
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
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)

        def OnDeleteLayer(event):
            Obj = event.GetEventObject()
            ind = Indx[Obj.GetId()]
            del data['Layers'][ind]
            data['Layer Seq'] = ' '.join([str(i+1) for i in range(len(data['Layers'])-2)])
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')
            wx.CallAfter(UpdateREFDModelsGrid,G2frame,data)

        def Recalculate(invalid,value,tc):
            if invalid:
                return
            if tc and tc.GetId() in Indx:  #sets slider
                Indx[tc.GetId()].SetValue(value)
            G2pwd.REFDModelFxn(Profile,Inst,Limits,Substances,data)
            x,xr,y = G2pwd.makeSLDprofile(data,Substances)
            ModelPlot(data,x,xr,y)
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')

        def OnMoveParm(event):
            Obj = event.GetEventObject()
            ilay,parm,parmObj = Indx[Obj.GetId()]
            move = Obj.GetValue()  # +1 or -1
            Obj.SetValue(0)
            data['Layers'][ilay][parm][0] += move
            parmObj.SetValue(data['Layers'][ilay][parm][0])
            Recalculate(False,1.0,None)

        def OnParmSlider(event):
            Obj = event.GetEventObject()
            ilay,parmObj = Indx[Obj.GetId()]
            value = Obj.GetValue()
            data['Layers'][ilay]['Thick'][0] = value
            parmObj.SetValue(data['Layers'][ilay]['Thick'][0])
            Recalculate(False,1.0,None)

        Indx = {}
        layerSizer = wx.BoxSizer(wx.VERTICAL)

        for ilay,layer in enumerate(data['Layers']):
            if not ilay:
                layerSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Top layer (superphase):'),0)
            elif ilay < len(data['Layers'])-1:
                layerSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Layer no. %d'%(ilay)),0)
            else:
                layerSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Bottom layer (substrate):'),0)
            midlayer = wx.BoxSizer(wx.HORIZONTAL)
            midlayer.Add(wx.StaticText(G2frame.dataWindow,label=' Substance: '),0,WACV)
            midName = data['Layers'][ilay]['Name']
            midSel = wx.ComboBox(G2frame.dataWindow,value=midName,
                choices=list(Substances.keys()),style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[midSel.GetId()] = ilay
            midSel.Bind(wx.EVT_COMBOBOX,OnSelect)
            midlayer.Add(midSel,0,WACV)
            if midName != 'vacuum':
                if midName != 'unit scatter':
                    midlayer.Add(wx.StaticText(G2frame.dataWindow,label=' Den. Mult.: '),0,WACV)
                    midlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Layers'][ilay]['DenMul'],0,
                        nDig=(10,4),OnLeave=Recalculate),0,WACV)
                    varBox = wx.CheckBox(G2frame.dataWindow,label='Refine?')
                    Indx[varBox.GetId()] = [ilay,'DenMul']
                    varBox.SetValue(data['Layers'][ilay]['DenMul'][1])
                    varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                    midlayer.Add(varBox,0,WACV)
                    realScatt = data['Layers'][ilay]['DenMul'][0]*Substances[midName]['Scatt density']
                    midlayer.Add(wx.StaticText(G2frame.dataWindow,
                        label=' Real scat. den.: %.4g'%(realScatt)),0,WACV)
                    imagScatt = data['Layers'][ilay]['DenMul'][0]*Substances[midName]['XImag density']
                    midlayer.Add(wx.StaticText(G2frame.dataWindow,
                        label=' Imag scat. den.: %.4g'%(imagScatt)),0,WACV)
                else:
                    realScatt = data['Layers'][ilay]['DenMul'][0]
                    midlayer.Add(wx.StaticText(G2frame.dataWindow,
                        label=' Real scat. den.: %.4g'%(realScatt)),0,WACV)
                    imagScatt = data['Layers'][ilay]['iDenMul'][0]
                    midlayer.Add(wx.StaticText(G2frame.dataWindow,
                        label=' Imag scat. den.: %.4g'%(imagScatt)),0,WACV)
            else:
                midlayer.Add(wx.StaticText(G2frame.dataWindow,label=', air or gas'),0,WACV)
            layerSizer.Add(midlayer)
            if midName == 'unit scatter':
                nxtlayer = wx.BoxSizer(wx.HORIZONTAL)
                nxtlayer.Add(wx.StaticText(G2frame.dataWindow,label=' Real Den. : '),0,WACV)
                nxtlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Layers'][ilay]['DenMul'],0,
                    nDig=(10,4),OnLeave=Recalculate),0,WACV)
                varBox = wx.CheckBox(G2frame.dataWindow,label='Refine?')
                Indx[varBox.GetId()] = [ilay,'DenMul']
                varBox.SetValue(data['Layers'][ilay]['DenMul'][1])
                varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                nxtlayer.Add(varBox,0,WACV)
                nxtlayer.Add(wx.StaticText(G2frame.dataWindow,label=' Imag Den. : '),0,WACV)
                nxtlayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Layers'][ilay]['iDenMul'],0,
                    nDig=(10,4),OnLeave=Recalculate),0,WACV)
                varBox = wx.CheckBox(G2frame.dataWindow,label='Refine?')
                Indx[varBox.GetId()] = [ilay,'iDenMul']
                varBox.SetValue(data['Layers'][ilay]['iDenMul'][1])
                varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                nxtlayer.Add(varBox,0,WACV)
                layerSizer.Add(nxtlayer)
            if midName != 'vacuum':
                if 'N' in Inst['Type'][0] and midName not in ['vacuum','unit scatter']:
                    magLayer = wx.BoxSizer(wx.HORIZONTAL)
                    magLayer.Add(wx.StaticText(G2frame.dataWindow,label=' Magnetic SLD: '),0,WACV)
                    magLayer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Layers'][ilay]['Mag SLD'],0,
                        nDig=(10,4),OnLeave=Recalculate),0,WACV)
                    varBox = wx.CheckBox(G2frame.dataWindow,label='Refine?')
                    Indx[varBox.GetId()] = [ilay,'Mag SLD']
                    varBox.SetValue(data['Layers'][ilay]['Mag SLD'][1])
                    varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                    magLayer.Add(varBox,0,WACV)
                    magLayer.Add(wx.StaticText(G2frame.dataWindow,
                        label=' Real+mag scat. den.: %.4g'%(realScatt+data['Layers'][ilay]['Mag SLD'][0])),0,WACV)
                    layerSizer.Add(magLayer)
                if ilay:
                    names = {'Rough':'Upper surface Roughness, '+Angstr,'Thick':'Layer Thickness, '+Angstr}
                    parmsline = wx.BoxSizer(wx.HORIZONTAL)
                    parms= ['Rough','Thick']
                    if ilay == len(data['Layers'])-1:
                        parms = ['Rough',]
                    if len(parms) > 1:
                        slide = wx.BoxSizer(wx.HORIZONTAL)
                        slide.Add(wx.StaticText(G2frame.dataWindow,label=' Layer thickness: '),0,WACV)
                        parmSldr = G2G.G2Slider(G2frame.dataWindow,minValue=0,maxValue=data['slider max'],value=data['Layers'][ilay]['Thick'][0])
                        parmSldr.Bind(wx.EVT_SLIDER,OnParmSlider)
                        slide.Add(parmSldr,1,wx.EXPAND)
                    for parm in parms:
                        parmsline.Add(wx.StaticText(G2frame.dataWindow,label=' %s: '%(names[parm])),0,WACV)
                        parmval = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['Layers'][ilay][parm],0,nDig=(10,2),OnLeave=Recalculate)
                        if parm =='Thick':
                            Indx[parmval.GetId()] = parmSldr
                        parmsline.Add(parmval,0,WACV)
                        parmSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(25,25))
                        Indx[parmSpin.GetId()] = [ilay,parm,parmval]
                        parmSpin.SetValue(0)
                        parmSpin.SetRange(-1,1)
                        parmSpin.Bind(wx.EVT_SPIN, OnMoveParm)
                        parmsline.Add(parmSpin,0,WACV)
                        varBox = wx.CheckBox(G2frame.dataWindow,label='Refine?')
                        Indx[varBox.GetId()] = [ilay,parm]
                        varBox.SetValue(data['Layers'][ilay][parm][1])
                        varBox.Bind(wx.EVT_CHECKBOX, OnCheckBox)
                        parmsline.Add(varBox,0,WACV)
                    layerSizer.Add(parmsline)
                    if len(parms) > 1:
                        Indx[parmSldr.GetId()] = [ilay,parmval] #parmval is always for Thick
                        layerSizer.Add(slide,1,wx.EXPAND)
            if ilay < len(data['Layers'])-1:
                newlayer = wx.BoxSizer(wx.HORIZONTAL)
                insert = wx.Button(G2frame.dataWindow,label='Insert')
                Indx[insert.GetId()] = ilay
                insert.Bind(wx.EVT_BUTTON,OnInsertLayer)
                newlayer.Add(insert)
                delet = wx.Button(G2frame.dataWindow,label='Delete')
                Indx[delet.GetId()] = ilay
                delet.Bind(wx.EVT_BUTTON,OnDeleteLayer)
                newlayer.Add(delet)
                layerSizer.Add(newlayer)
                G2G.HorizontalLine(layerSizer,G2frame.dataWindow)

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
            G2pwpl.PlotPatterns(G2frame,plotType='REFD')

    # start of UpdateREFDModelsGrid
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.REFDModelMenu)
    Substances = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Substances'))['Substances']
    ProfDict,Profile,Name = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[:3]
    #patch
    if 'slider max' not in data:
        data['slider max'] = 2000.
    if 'ifDQ' not in ProfDict:
        ProfDict['ifDQ'] = np.any(Profile[5])
        data['dQ type'] = 'None'
    Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))
    Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
    G2frame.dataWindow.ClearData()
    G2frame.Bind(wx.EVT_MENU, OnCopyModel, id=G2G.wxID_MODELCOPY)
    G2frame.Bind(wx.EVT_MENU, OnModelPlot, id=G2G.wxID_MODELPLOT)
    G2frame.Bind(wx.EVT_MENU, OnFitModel, id=G2G.wxID_MODELFIT)
    G2frame.Bind(wx.EVT_MENU, OnFitModelAll, id=G2G.wxID_MODELFITALL)
    G2frame.Bind(wx.EVT_MENU, OnUnDo, id=G2G.wxID_MODELUNDO)
    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label=' Reflectometry fitting for: '+Name),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Controls:'))
    mainSizer.Add(ControlSizer())
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Global parameters:'))
    mainSizer.Add(OverallSizer())
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    Nlayers = len(data['Layers'])
    if Nlayers > 2:
        if 'Layer Seq' not in data:
            data['Layer Seq'] = ' '.join([str(i+1) for i in range(Nlayers-2)])
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Layer sequence: '),0,WACV)
        repseq = wx.TextCtrl(G2frame.dataWindow,value = data['Layer Seq'],style=wx.TE_PROCESS_ENTER,size=(500,25))
        repseq.Bind(wx.EVT_TEXT_ENTER,OnRepSeq)
        repseq.Bind(wx.EVT_KILL_FOCUS,OnRepSeq)
        lineSizer.Add(repseq,0,WACV)
        mainSizer.Add(lineSizer)
        Str = ' Use sequence nos. from:'
        for ilay,layer in enumerate(data['Layers'][1:-1]):
            Str += ' %d: %s'%(ilay+1,layer['Name'])
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=Str))
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' NB: Repeat sequence by e.g. 6*(1 2) '))
    G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
    mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Layers: scatt. densities are 10%scm%s = 10%s%s%s'%(Pwr10,Pwrm2,Pwrm6,Angstr,Pwrm2)))
    mainSizer.Add(LayerSizer())
    G2frame.dataWindow.SetDataSize()

################################################################################
#####  PDF controls
################################################################################

def computePDF(G2frame,data):
    '''Calls :func:`GSASIIpwd.CalcPDF` to compute the PDF and put into the data tree array.
    Called from OnComputePDF and OnComputeAllPDF and OnComputeAllPDF in
    GSASIIimgGUI.py
    '''

    xydata = {}
    problem = False
    for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
        name = data[key]['Name']
        if name.strip():
            pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
            if not pId:
                print(key,'Entry',name,'Not found.')
                problem = True
                continue
            xydata[key] = G2frame.GPXtree.GetItemPyData(pId)
    if problem:
        print('PDF computation aborted')
        return
    powId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['Sample']['Name'])
    limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId,'Limits'))
    Xlimits = [limits[1][0],limits[0][1]]       #use lower limit but ignore upper limit
    inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId,'Instrument Parameters'))[0]
    auxPlot = G2pwd.CalcPDF(data,inst,Xlimits,xydata)
    try:
        data['I(Q)'] = xydata['IofQ']
        data['S(Q)'] = xydata['SofQ']
        data['F(Q)'] = xydata['FofQ']
        data['G(R)'] = xydata['GofR']
        data['g(r)'] = xydata['gofr']
        return auxPlot
    except:   # PDF Calc aborted
        pass

def OptimizePDF(G2frame,data,showFit=True,maxCycles=5):
    '''Optimize the PDF to minimize the difference between G(r) and the expected value for
    low r (-4 pi r #density).
    '''
    xydata = {}
    for key in ['Sample','Sample Bkg.','Container','Container Bkg.']:
        name = data[key]['Name']
        if name:
            xydata[key] = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name))
    powName = data['Sample']['Name']
    powId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,powName)
    limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId,'Limits'))
    Xlimits = [limits[1][0],limits[0][1]]       #use lower limit but ignore upper limit
    inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId,'Instrument Parameters'))[0]
    res = G2pwd.OptimizePDF(data,xydata,Xlimits,inst,showFit,maxCycles)
    return res['success']

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
                mult.SetValue(data[fileKey]['Mult'])
                ResetFlatBkg()
                wx.CallAfter(OnComputePDF,None)

            def OnMoveMult(event):
                data[key]['Mult'] += multSpin.GetValue()*0.01
                mult.SetValue(data[key]['Mult'])
                multSpin.SetValue(0)
                wx.CallAfter(OnComputePDF,None)

            def OnMult(invalid,value,tc):
                if invalid: return
                ResetFlatBkg()
                wx.CallAfter(OnComputePDF,None)

            def OnRefMult(event):
                item['Refine'] = refMult.GetValue()
                if item['Refine']:
                    G2frame.GetStatusBar().SetStatusText('Be sure Mult is close to anticipated value. '+   \
                        'Suggest setting Flat Bkg. to 0 before Optimize Mult',1)

            def GetExposure(backFile):
                dataId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'PWDR'+dataFile[4:])
                dataComments = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,dataId,'Comments'))
                if not backFile:
                    return -1.
                backId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,backFile)
                backComments = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,backId,'Comments'))
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
            fileSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' '+key+' file:'),0,WACV)
            fileName = wx.ComboBox(G2frame.dataWindow,value=item['Name'],choices=fileList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            itemDict[fileName.GetId()] = [key,'Name','%s']
            fileName.Bind(wx.EVT_COMBOBOX,OnSelectFile)
            fileSizer.Add(fileName,0,)
            fileSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label='Multiplier:'),0,WACV)
            mulBox = wx.BoxSizer(wx.HORIZONTAL)
            mult = G2G.ValidatedTxtCtrl(G2frame.dataWindow,item,'Mult',nDig=(10,3),
                typeHint=float,OnLeave=OnMult)
            mulBox.Add(mult,0,)
            multSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(20,25))
            multSpin.SetRange(-1,1)
            multSpin.SetValue(0)
            multSpin.Bind(wx.EVT_SPIN, OnMoveMult)
            mulBox.Add(multSpin,0,WACV)
            fileSizer.Add(mulBox,0,WACV)
            if 'Refine' in item and item['Name'] and 'Sample' in key:
                refMult = wx.CheckBox(parent=G2frame.dataWindow,label='Refine?')
                refMult.SetValue(item['Refine'])
                refMult.Bind(wx.EVT_CHECKBOX, OnRefMult)
                fileSizer.Add(refMult,0,WACV)
            else:
                fileSizer.Add((5,5),0)

        def ResetFlatBkg():
            Smin = np.min(G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'PWDR'+dataFile[4:]))[1][1])
            Bmin = 0; Cmin = 0.; Cmul = 0.; CBmin = 0.
            if data['Sample Bkg.']['Name']:
                Bmin = np.min(G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['Sample Bkg.']['Name']))[1][1])
                Smin += Bmin*data['Sample Bkg.']['Mult']
            if data['Container']['Name']:
                Cmin = np.min(G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['Container']['Name']))[1][1])
                Cmul = data['Container']['Mult']
                if data['Container Bkg.']['Name']:
                    CBmin = np.min(G2frame.GPXtree.GetItemPyData(
                        G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['Container Bkg.']['Name']))[1][1])
                    Cmin += CBmin*data['Container Bkg.']['Mult']
                Smin += Cmul*Cmin
            data['Flat Bkg'] = max(0,Smin)
            G2frame.flatBkg.SetValue(data['Flat Bkg'])

        PDFfileSizer = wx.BoxSizer(wx.VERTICAL)
        PDFfileSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=' PDF data files: '),0)
        PDFfileSizer.Add((5,5),0)
        if 'C' in inst['Type'][0]:
            str = ' Sample file: PWDR%s   Wavelength, A: %.5f  Energy, keV: %.3f  Polariz.: %.2f '%(dataFile[4:],wave,keV,polariz)
            PDFfileSizer.Add(wx.StaticText(parent=G2frame.dataWindow,label=str),0)
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

            def OnElemNum(invalid,value,tc):
                if invalid: return
                data['Form Vol'] = max(10.0,SumElementVolumes())
                wx.CallAfter(UpdatePDFGrid,G2frame,data)
                wx.CallAfter(OnComputePDF,tc.event)

            elemSizer.Add(wx.StaticText(parent=G2frame.dataWindow,
                label=' Element: '+'%2s'%(ElData['Symbol'])+' * '),0,WACV)
            num = G2G.ValidatedTxtCtrl(G2frame.dataWindow,ElData,'FormulaNo',nDig=(10,3),xmin=0.0,
                typeHint=float,OnLeave=OnElemNum)
            elemSizer.Add(num,0,WACV)
            elemSizer.Add(wx.StaticText(parent=G2frame.dataWindow,
                label="f': %.3f"%(ElData['fp'])+' f": %.3f'%(ElData['fpp'])+' mu: %.2f barns'%(ElData['mu']) ),
                0,WACV)

        def AfterChange(invalid,value,tc):
            if invalid: return
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            wx.CallAfter(OnComputePDF,tc.event)

        def OnGeometry(event):
            data['Geometry'] = geometry.GetValue()
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            #UpdatePDFGrid(G2frame,data)
            wx.CallAfter(OnComputePDF,event)

        sampleSizer = wx.BoxSizer(wx.VERTICAL)
        if not ElList:
            sampleSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Sample information: fill in this 1st'),0)
        else:
            sampleSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Sample information: '),0)
        sampleSizer.Add((5,5),0)
        Abs = G2lat.CellAbsorption(ElList,data['Form Vol'])
        Trans = G2pwd.Transmission(data['Geometry'],Abs*data['Pack'],data['Diam'])
        elemSizer = wx.FlexGridSizer(0,3,5,1)
        for El in ElList:
            FillElemSizer(elemSizer,ElList[El])
        sampleSizer.Add(elemSizer,0)
        sampleSizer.Add((5,5),0)
        midSizer = wx.BoxSizer(wx.HORIZONTAL)
        midSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Formula volume: '),0,WACV)
        formVol = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Form Vol',nDig=(10,3),xmin=10.0,
            typeHint=float,OnLeave=AfterChange)
        midSizer.Add(formVol,0)
        midSizer.Add(wx.StaticText(G2frame.dataWindow,
            label=' Theoretical absorption: %.4f cm-1 Sample absorption: %.4f cm-1'%(Abs,Abs*data['Pack'])),
            0,WACV)
        sampleSizer.Add(midSizer,0)
        sampleSizer.Add((5,5),0)
        geoBox = wx.BoxSizer(wx.HORIZONTAL)
        geoBox.Add(wx.StaticText(G2frame.dataWindow,label=' Sample geometry: '),0,WACV)
        if 'C' in inst['Type'][0]:
            choice = ['Cylinder','Bragg-Brentano','Tilting flat plate in transmission','Fixed flat plate']
        else:
            choice = ['Cylinder',]
        geometry = wx.ComboBox(G2frame.dataWindow,value=data['Geometry'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        geometry.Bind(wx.EVT_COMBOBOX, OnGeometry)
        geoBox.Add(geometry,0)
        geoBox.Add(wx.StaticText(G2frame.dataWindow,label=' Sample diameter/thickness, mm: '),0,WACV)
        diam = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Diam',nDig=(10,3),xmin=0.01,
            typeHint=float,OnLeave=AfterChange)
        geoBox.Add(diam,0)
        sampleSizer.Add(geoBox,0)
        sampleSizer.Add((5,5),0)
        geoBox = wx.BoxSizer(wx.HORIZONTAL)
        geoBox.Add(wx.StaticText(G2frame.dataWindow,label=' Packing: '),0,WACV)
        pack = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Pack',nDig=(10,2),xmin=0.01,
            typeHint=float,OnLeave=AfterChange)
        geoBox.Add(pack,0)
        geoBox.Add(wx.StaticText(G2frame.dataWindow,label=' Sample transmission: %.3f %%'%(100.*Trans)),0,WACV)
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
            G2frame.GetStatusBar().SetStatusText('',1)
            wx.BeginBusyCursor()
            try:
                data['Ruland'] = 0.01       #always set small to start
                OptimizePDF(G2frame,data)
            finally:
                wx.EndBusyCursor()
            OnComputePDF(event)
            wx.CallAfter(UpdatePDFGrid,G2frame,data)

        def AfterChangeNoRefresh(invalid,value,tc):
            if invalid: return
            if tc.GetId() in Indx:
                Indx[tc.GetId()][0].SetValue(int(value*Indx[tc.GetId()][1]))
            wx.CallAfter(OnComputePDF,None)

        def OnDetType(event):
            data['DetType'] = detType.GetValue()
            wx.CallAfter(UpdatePDFGrid,G2frame,data)
            wx.CallAfter(OnComputePDF,None)

        def OnFlatSpin(event):
            data['Flat Bkg'] += flatSpin.GetValue()*0.01*data['IofQmin']
            G2frame.flatBkg.SetValue(data['Flat Bkg'])
            flatSpin.SetValue(0)
            wx.CallAfter(OnComputePDF,None)

        def OnBackSlider(event):
            value = int(backSldr.GetValue())/100.
            data['BackRatio'] = value
            backVal.SetValue(data['BackRatio'])
            wx.CallAfter(OnComputePDF,None)

        def OnRulSlider(event):
            value = int(rulandSldr.GetValue())/100.
            data['Ruland'] = max(0.001,value)
            rulandWdt.SetValue(data['Ruland'])
            wx.CallAfter(OnComputePDF,None)

        def OnGRscaleSlider(event):
            value = int(gscaleSldr.GetValue())/50.
            data['GR Scale'] = max(0.1,min(2.,value))
            gscale.SetValue(data['GR Scale'])
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

        Indx = {}
        sfgSizer = wx.BoxSizer(wx.VERTICAL)
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' S(Q)->F(Q)->G(r) controls: '),0,WACV)
        sqBox.Add((1,1),1,wx.EXPAND,1)
        optB = wx.Button(G2frame.dataWindow,label='Optimize PDF',style=wx.BU_EXACTFIT)
        optB.Bind(wx.EVT_BUTTON, OnOptimizePDF)
        sqBox.Add(optB,0,WACV)
        sfgSizer.Add(sqBox,0,wx.EXPAND)

        sfgSizer.Add((5,5),0)
        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' Detector type: '),0,WACV)
        choice = ['Area detector','Point detector']
        detType = wx.ComboBox(G2frame.dataWindow,value=data['DetType'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        detType.Bind(wx.EVT_COMBOBOX, OnDetType)
        sqBox.Add(detType,0)
        if data['DetType'] == 'Area detector':
            sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' IP transmission coeff.: '),0,WACV)
            obliqCoeff = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'ObliqCoeff',nDig=(10,3),xmin=0.0,xmax=1.0,
                typeHint=float,OnLeave=AfterChangeNoRefresh)
            sqBox.Add(obliqCoeff,0)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' Flat Bkg.: '),0,WACV)
        G2frame.flatBkg = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Flat Bkg',nDig=(10,0),
                typeHint=float,OnLeave=AfterChangeNoRefresh)
        sqBox.Add(G2frame.flatBkg,0)
        flatSpin = wx.SpinButton(G2frame.dataWindow,style=wx.SP_VERTICAL,size=wx.Size(20,25))
        flatSpin.SetRange(-1,1)
        flatSpin.SetValue(0)
        flatSpin.Bind(wx.EVT_SPIN, OnFlatSpin)
        sqBox.Add(flatSpin,0,WACV)
        sqBox.Add((1,1),1,wx.EXPAND,1)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label='Rmin: '),0,WACV)
        rmin = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Rmin',nDig=(5,1),
                xmin=0.01,typeHint=float,size=wx.Size(50,20))
        sqBox.Add(rmin,0,WACV)
        sfgSizer.Add(sqBox,0,wx.EXPAND)

        bkBox = wx.BoxSizer(wx.HORIZONTAL)
        bkBox.Add(wx.StaticText(G2frame.dataWindow,label=' Background ratio: '),0,WACV)
        backSldr = G2G.G2Slider(parent=G2frame.dataWindow,style=wx.SL_HORIZONTAL,
            value=int(100*data['BackRatio']))
        bkBox.Add(backSldr,1,wx.EXPAND)
        backSldr.Bind(wx.EVT_SLIDER, OnBackSlider)
        backVal = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'BackRatio',nDig=(10,3),xmin=0.0,xmax=1.0,
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        Indx[backVal.GetId()] = [backSldr,100]
        bkBox.Add(backVal,0,WACV)
        sfgSizer.Add(bkBox,0,wx.EXPAND)

        if 'XC' in inst['Type'][0]:
            sqBox = wx.BoxSizer(wx.HORIZONTAL)
            sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' Ruland width: '),0,WACV)
            rulandSldr = G2G.G2Slider(parent=G2frame.dataWindow,style=wx.SL_HORIZONTAL,
                value=int(100*data['Ruland']))
            sqBox.Add(rulandSldr,1,wx.EXPAND)
            rulandSldr.Bind(wx.EVT_SLIDER, OnRulSlider)
            rulandWdt = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Ruland',nDig=(10,3),xmin=0.001,xmax=1.0,
                typeHint=float,OnLeave=AfterChangeNoRefresh)
            Indx[rulandWdt.GetId()] = [rulandSldr,100]
            sqBox.Add(rulandWdt,0,WACV)
            sfgSizer.Add(sqBox,0,wx.EXPAND)

        gscaleBox = wx.BoxSizer(wx.HORIZONTAL)
        gscaleBox.Add(wx.StaticText(G2frame.dataWindow,label=' G(R) scale: '),0,WACV)
        gscaleSldr = G2G.G2Slider(parent=G2frame.dataWindow,style=wx.SL_HORIZONTAL,
            value=int(50*data['GR Scale']))
        gscaleBox.Add(gscaleSldr,1,wx.EXPAND)
        gscaleSldr.Bind(wx.EVT_SLIDER, OnGRscaleSlider)
        gscale = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'GR Scale',nDig=(10,3),xmin=0.1,xmax=2.,
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        Indx[gscale.GetId()] = [gscaleSldr,50]
        gscaleBox.Add(gscale,0,WACV)
        sfgSizer.Add(gscaleBox,0,wx.EXPAND)

        sqBox = wx.BoxSizer(wx.HORIZONTAL)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' Scaling Q-range: '),0,WACV)
        SQmin = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['QScaleLim'],0,nDig=(10,3),
            xmin=qLimits[0],xmax=.95*data['QScaleLim'][1],
            typeHint=float,OnLeave=AfterChangeNoRefresh)
        sqBox.Add(SQmin,0,WACV)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' to Qmax '),0,WACV)
        SQmax = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data['QScaleLim'],1,nDig=(10,3),
            xmin=qLimits[0],xmax=qLimits[1],typeHint=float,OnLeave=NewQmax)
        sqBox.Add(SQmax,0,WACV)
        resetQ = wx.Button(G2frame.dataWindow,label='Reset?',style=wx.BU_EXACTFIT)
        sqBox.Add(resetQ,0,WACV)
        resetQ.Bind(wx.EVT_BUTTON, OnResetQ)
        sqBox.Add(wx.StaticText(G2frame.dataWindow,label=' Plot Rmax: '),0,WACV)
        rmax = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'Rmax',nDig=(10,1),xmin=10.,xmax=200.,
            typeHint=float,OnLeave=AfterChangeNoRefresh,size=wx.Size(50,20))
        sqBox.Add(rmax,0,WACV)
        lorch = wx.CheckBox(parent=G2frame.dataWindow,label='Lorch damping?')
        lorch.SetValue(data['Lorch'])
        lorch.Bind(wx.EVT_CHECKBOX, OnLorch)
        sqBox.Add(lorch,0,WACV)
        noRing = wx.CheckBox(parent=G2frame.dataWindow,label='Suppress G(0) ringing?')
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
                Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,newName)
                pId = G2gd.GetGPXtreeItemId(G2frame,Id,'PDF Controls')
                subData = G2frame.GPXtree.GetItemPyData(pId)['G(R)']
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
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,data['diffGRname'])
            if Id == 0: return
            pId = G2gd.GetGPXtreeItemId(G2frame,Id,'PDF Controls')
            if pId == 0: return
            subData = G2frame.GPXtree.GetItemPyData(pId)['G(R)']
            data['delt-G(R)'][1] = np.array([subData[1][0],data['G(R)'][1][1]-data['diffMult']*subData[1][1]])
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='delt-G(R)')

        diffSizer = wx.BoxSizer(wx.HORIZONTAL)
        fileList = [''] + GetFileList(G2frame,'PDF')
        diffSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Subtract G(R) for: '),0,WACV)
        grName = wx.ComboBox(G2frame.dataWindow,value=data['diffGRname'],choices=fileList,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        grName.Bind(wx.EVT_COMBOBOX,OnSelectGR)
        diffSizer.Add(grName,0,WACV)
        if data['diffGRname']:
            diffSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Mult: '),0,WACV)
            mult = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'diffMult',nDig=(10,3),
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
        TextList = GetFileList(G2frame,'PDF')
        Source = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy controls to','There must be more than one "PDF" pattern')
            return
        od = {'label_1':'Only refine flag','value_1':False,'label_2':'Only Lorch flag','value_2':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy PDF controls','Copy controls from '+Source+' to:',TextList,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PDFlist = [TextList[i] for i in dlg.GetSelections()]
                for item in PDFlist:
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    olddata = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'PDF Controls'))
                    if od['value_1']:
                        olddata['Sample Bkg.']['Refine'] = data['Sample Bkg.']['Refine']    #only one flag
                    elif od['value_2']:
                        olddata['Lorch'] = data['Lorch']    #only one flag
                    else:
                        sample = olddata['Sample']
                        olddata.update(copy.deepcopy(data))
                        olddata['Sample'] = sample
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'PDF Controls'),olddata)
                G2frame.GetStatusBar().SetStatusText('PDF controls copied',1)
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
                    if item[:] not in ['Sample','I(Q)','S(Q)','F(Q)','G(R)']:
                        File.write(item+':'+data[item]+'\n')
                File.close()
                print ('PDF controls saved to: '+filename)
        finally:
            dlg.Destroy()

    def OnLoadPDFControls(event):
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II PDF controls file', pth, '',
            'PDF controls files (*.pdfprm)|*.pdfprm',wx.FD_OPEN)
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
                        newdata[key] = val.strip()
                    S = File.readline()
                File.close()
                data.update(newdata)
        finally:
            dlg.Destroy()
        OnComputePDF(event)
        wx.CallAfter(UpdatePDFGrid,G2frame,data)

    def OnAddElement(event):
        ElList = data['ElList']
        choice = list(ElList.keys())
        PE = G2elemGUI.PickElements(G2frame,choice)
        if PE.ShowModal() == wx.ID_OK:
            for El in PE.Elem:
                if El not in ElList:
                    try:
                        data['ElList'][El] = copy.deepcopy(G2elem.GetElInfo(El,inst))
                        data['ElList'][El]['FormulaNo'] = 1.0
                    except IndexError: # happens with element Q
                        pass
            data['Form Vol'] = max(10.0,SumElementVolumes())
        PE.Destroy()
        wx.CallAfter(UpdatePDFGrid,G2frame,data)

    def OnDeleteElement(event):
        ElList = data['ElList']
        choice = list(ElList.keys())
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
            OnAddElement(event)
        auxPlot = computePDF(G2frame,data)
        if auxPlot is None: return
        G2frame.GetStatusBar().SetStatusText('PDF computed',1)
        for plot in auxPlot:
            XY = np.array(plot[:2])
            G2plt.PlotXY(G2frame,[XY,],Title=plot[2])
        if event is not None:
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='I(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='S(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='F(Q)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='g(r)')
            G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='G(R)')
        else:
            G2plt.PlotISFG(G2frame,data,newPlot=False)

    def OnComputeAllPDF(event):
        print('Calculating PDFs...')
        choices = []
        if G2frame.GPXtree.GetCount():
            Id, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while Id:
                Name = G2frame.GPXtree.GetItemText(Id)
                if Name.startswith('PDF '):
                    Data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'PDF Controls'))
                    if not Data['ElList']:
                        print('  No chemical formula for {}'.format(Name))
                    else:
                        choices.append(Name)
                Id, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        if not choices:
            print('  No PDFs to compute\n')
            return
        od = {'label_1':'Optimize PDFs','value_1':True}
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select PDFs to compute','Select PDFs',
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
        Id, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
        N = 0
        try:
            while Id:
                Name = G2frame.GPXtree.GetItemText(Id)
                if Name in Names:
                    N += 1
                    msg = 'PDFs done: {} of {}'.format(N-1,len(Names))
                    if not pgbar.Update(N,msg)[0]:
                        pgbar.Destroy()
                        break
                    pId = G2gd.GetGPXtreeItemId(G2frame,Id,'PDF Controls')
                    Data = G2frame.GPXtree.GetItemPyData(pId)
                    print('  Computing {}'.format(Name))
                    computePDF(G2frame,Data)
                    if od['value_1']:
                        notConverged += not OptimizePDF(G2frame,Data,maxCycles=10)
                    computePDF(G2frame,Data)
                    G2frame.GPXtree.SetItemPyData(pId,Data)
                Id, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        finally:
            pgbar.Destroy()
        if od['value_1']:
            msg = '{}/{} PDFs computed; {} unconverged'.format(N,len(Names),notConverged)
        else:
            msg = '{}/{} PDFs computed'.format(N,len(Names))
        G2frame.GetStatusBar().SetStatusText(msg,1)
        print(msg)
        # what item is being plotted? -- might be better to select from tree
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='I(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='S(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='F(Q)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='G(R)')
        G2plt.PlotISFG(G2frame,data,newPlot=True,plotType='g(r)')

    # Routine UpdatePDFGrid starts here
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PDFMenu)
    global inst
    tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
    tof2q = lambda t,C:2.0*math.pi*C/t
    dataFile = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    powName = 'PWDR'+dataFile[4:]
    powId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, powName)
    if powId: # skip if no matching PWDR entry
        fullLimits,limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId, 'Limits'))[:2]
        inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,powId, 'Instrument Parameters'))[0]
        if 'C' in inst['Type'][0]:
            wave = G2mth.getWave(inst)
            keV = G2mth.wavekE(wave)
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
        if 'GR Scale' not in data:
            data['GR Scale'] = 1.0
    if powId:
        G2frame.dataWindow.PDFMenu.EnableTop(0,enable=True)
    else:
        G2frame.dataWindow.PDFMenu.EnableTop(0,enable=False)
    G2frame.Bind(wx.EVT_MENU, OnCopyPDFControls, id=G2G.wxID_PDFCOPYCONTROLS)
    G2frame.Bind(wx.EVT_MENU, OnSavePDFControls, id=G2G.wxID_PDFSAVECONTROLS)
    G2frame.Bind(wx.EVT_MENU, OnLoadPDFControls, id=G2G.wxID_PDFLOADCONTROLS)
    G2frame.Bind(wx.EVT_MENU, OnAddElement, id=G2G.wxID_PDFADDELEMENT)
    G2frame.Bind(wx.EVT_MENU, OnDeleteElement, id=G2G.wxID_PDFDELELEMENT)
    G2frame.Bind(wx.EVT_MENU, OnComputePDF, id=G2G.wxID_PDFCOMPUTE)
    G2frame.Bind(wx.EVT_MENU, OnComputeAllPDF, id=G2G.wxID_PDFCOMPUTEALL)

    G2frame.dataWindow.ClearData()
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    if powId:
        ElList = data['ElList']
        mainSizer.Add(PDFFileSizer())
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        mainSizer.Add(SampleSizer())
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        mainSizer.Add(SFGctrlSizer())
        G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
        mainSizer.Add(DiffSizer())
    else:
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='Controls for %s:\n'%powName))
        mainSizer.Add(DiffSizer())
    G2frame.dataWindow.SetDataSize()

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
        limitBox.Add(wx.StaticText(G2frame.dataWindow,label=' PDF Limits: '),0,WACV)
        lowLim = G2G.ValidatedTxtCtrl(G2frame.dataWindow,peaks['Limits'],0,nDig=(10,3),
            xmin=0.,xmax=10.,typeHint=float,OnLeave=NewLim)
        limitBox.Add(lowLim,0,WACV)
        highLim = G2G.ValidatedTxtCtrl(G2frame.dataWindow,peaks['Limits'],1,nDig=(10,3),
            xmin=peaks['Limits'][0],xmax=10.,typeHint=float,OnLeave=NewLim)
        limitBox.Add(highLim,0,WACV)
        return limitBox

    def backSizer():

        def NewBack(invalid,value,tc):
            if invalid:
                return
            G2plt.PlotISFG(G2frame,data,newPlot=False,plotType='G(R)',peaks=peaks)

        def OnRefBack(event):
            peaks['Background'][2] = refbk.GetValue()

        backBox = wx.BoxSizer(wx.HORIZONTAL)
        backBox.Add(wx.StaticText(G2frame.dataWindow,label=' Background slope: '),0,WACV)
        slope = G2G.ValidatedTxtCtrl(G2frame.dataWindow,peaks['Background'][1],1,nDig=(10,3),
            xmin=-4.*np.pi,xmax=0.,typeHint=float,OnLeave=NewBack)
        backBox.Add(slope,0,WACV)
        refbk = wx.CheckBox(parent=G2frame.dataWindow,label=' Refine?')
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
        rowLabels = [str(i) for i in range(len(peaks['Peaks']))]
        peakTable = G2G.Table(peaks['Peaks'],rowLabels=rowLabels,colLabels=colLabels,types=Types)
        PDFPeaks = G2G.GSGrid(G2frame.dataWindow)
        PDFPeaks.SetRowLabelSize(45)
        PDFPeaks.SetTable(peakTable,True)
        PDFPeaks.AutoSizeColumns(False)
        PDFPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, PeaksRefine)
        PDFPeaks.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, ElTypeSelect)

        peakBox = wx.BoxSizer(wx.VERTICAL)
        peakBox.Add(wx.StaticText(G2frame.dataWindow,label=' PDF Peaks:'),0)
        peakBox.Add(PDFPeaks,0)

        return peakBox

    def OnCopyPDFPeaks(event):
        TextList = GetFileList(G2frame,'PDF')
        Source = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        if len(TextList) == 1:
            G2frame.ErrorDialog('Nothing to copy PDF peaks to','There must be more than one "PDF" pattern')
            return
        od = {'label_1':'Only refine flags','value_1':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Copy PDF peaks','Copy peaks from '+Source+' to:',TextList,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                PDFlist = [TextList[i] for i in dlg.GetSelections()]
                for item in PDFlist:
                    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,item)
                    olddata = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'PDF Peaks'))
                    if od['value_1']:
                        olddata['Background'][2] = peaks['Background'][2]
                        for ip,peak in enumerate(olddata['Peaks']):
                            peak[3] = peaks['Peaks'][ip][3]
                    else:
                        olddata.update(copy.deepcopy(peaks))
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'PDF Peaks'),olddata)
                G2frame.GetStatusBar().SetStatusText('PDF peaks copied',1)
        finally:
            dlg.Destroy()

    def OnFitPDFpeaks(event):
        PatternId = G2frame.PatternId
        data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'PDF Controls'))
        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'PDF Peaks'))
        if not peaks:
            G2frame.ErrorDialog('No peaks!','Nothing to fit!')
            return
        newpeaks = G2pwd.PDFPeakFit(peaks,data['G(R)'])[0]
        print ('PDF peak fit finished')
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId, 'PDF Peaks'),newpeaks)
        G2plt.PlotISFG(G2frame,data,peaks=newpeaks,newPlot=False)
        wx.CallAfter(UpdatePDFPeaks,G2frame,newpeaks,data)

    def OnFitAllPDFpeaks(event):
        Names = G2gd.GetGPXtreeDataNames(G2frame,['PDF ',])
        od = {'label_1':'Copy to next','value_1':False,'label_2':'Reverse order','value_2':False}
        dlg = G2G.G2MultiChoiceDialog(G2frame,'PDF peak fitting','Select PDFs to fit:',Names,extraOpts=od)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential PDF peak fit results')
                if Id:
                    SeqResult = G2frame.GPXtree.GetItemPyData(Id)
                else:
                    SeqResult = {}
                    Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential PDF peak fit results')
                SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
                items = dlg.GetSelections()
                if od['value_2']:
                    items.reverse()
                newpeaks = None
                G2frame.EnablePlot = False
                for item in items:
                    name = Names[item]
                    pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                    data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId, 'PDF Controls'))
                    if od['value_1'] and newpeaks is not None:
                        peaks = copy.deepcopy(newpeaks)
                    else:
                        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'PDF Peaks'))
                    newpeaks,vals,varyList,sigList,parmDict,Rvals = G2pwd.PDFPeakFit(peaks,data['G(R)'])
                    if vals is None:
                        print ('Nothing varied!')
                        dlg.Destroy()
                        return
                    SeqResult[name] = {'variables':vals,'varyList':varyList,'sig':sigList,'Rvals':Rvals,
                        'covMatrix':np.eye(len(varyList)),'title':name,'parmDict':parmDict}
                    G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId, 'PDF Peaks'),newpeaks)
                SeqResult['histNames'] = Names
                G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
                G2plt.PlotISFG(G2frame,data,peaks=newpeaks,newPlot=False)
                G2frame.GPXtree.SetItemPyData(Id,SeqResult)
                G2frame.GPXtree.SelectItem(Id)
                print ('All PDFs peak fitted - results in Sequential PDF peak fit results')
            else:
                print ('Sequential fit cancelled')
        finally:
            dlg.Destroy()

    def OnClearPDFpeaks(event):
        peaks['Peaks'] = []
        G2plt.PlotISFG(G2frame,data,peaks=peaks,newPlot=False)
        wx.CallAfter(UpdatePDFPeaks,G2frame,peaks,data)

    # start of UpdatePDFPeaks
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PDFPksMenu)
    G2frame.Bind(wx.EVT_MENU, OnCopyPDFPeaks, id=G2G.wxID_PDFCOPYPEAKS)
    G2frame.Bind(wx.EVT_MENU, OnFitPDFpeaks, id=G2G.wxID_PDFPKSFIT)
    G2frame.Bind(wx.EVT_MENU, OnFitAllPDFpeaks, id=G2G.wxID_PDFPKSFITALL)
    G2frame.Bind(wx.EVT_MENU, OnClearPDFpeaks, id=G2G.wxID_CLEARPDFPEAKS)
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=' PDF peak fit controls:'))
    mainSizer.Add((5,5),0)
    mainSizer.Add(limitSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(backSizer())
    if len(peaks['Peaks']):
        mainSizer.Add((5,5),0)
        mainSizer.Add(peakSizer())
    G2frame.dataWindow.SetDataSize()
