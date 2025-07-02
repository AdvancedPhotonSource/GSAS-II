# -*- coding: utf-8 -*-
#GSASII - phase data display routines
'''

Routines for Phase dataframes follow.
'''
from __future__ import division, print_function
import platform
import os
import wx
import wx.grid as wg
import wx.lib.scrolledpanel as wxscroll
import matplotlib as mpl
#import math
import copy
import time
import sys
import random as ran
import subprocess as subp
#import platform
import shutil

import numpy as np
import numpy.linalg as nl
import numpy.ma as ma
import scipy.optimize as so
from . import GSASIIpath
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIElem as G2elem
from . import GSASIIElemGUI as G2elemGUI
from . import GSASIIddataGUI as G2ddG
from . import GSASIIplot as G2plt
from . import GSASIIpwdplot as G2pwpl
# if GSASIIpath.GetConfigValue('debug'):
#     print('Debug reloading',G2plt)
#     import imp
#     imp.reload(G2plt)
from . import GSASIIdataGUI as G2gd
from . import GSASIImiscGUI as G2IO
from . import GSASIIstrMain as G2stMn
from . import GSASIIstrIO as G2stIO
from . import GSASIImath as G2mth
from . import GSASIIpwd as G2pwd
from . import GSASIIobj as G2obj
from . import GSASIIctrlGUI as G2G
from . import GSASIIphsGUI as G2phsG
from . import GSASIIfiles as G2fil
from . import GSASIIconstrGUI as G2cnstG
from . import atmdata
from . import ISODISTORT as ISO
from . import SUBGROUPS

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

try:
    VERY_LIGHT_GREY = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
    WHITE = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
    BLACK = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT)
    RED = wx.Colour(255,0,0)
    WACV = wx.ALIGN_CENTER_VERTICAL
except:
    pass
mapDefault = G2elem.mapDefault
TabSelectionIdDict = {}
# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda x,y: 180.*np.arctan2(y,x)/np.pi
is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
sqt2 = np.sqrt(2.)
sqt3 = np.sqrt(3.)

# previous rigid body selections
prevResId = None
prevVecId = None
prevSpnId = None

GkDelta = chr(0x0394)
Angstr = chr(0x00c5)

RMCmisc = {}
ranDrwDict = {}
ranDrwDict['atomList'] = []
DrawStyleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']

#### Dysnomia (MEM) Data page ##############################################################################
def UpdateDysnomia(G2frame,data):
    ''' Present the controls for running Dysnomia
    '''
    def OnOptMeth(event):
        DysData['Optimize'] = OptMeth.GetValue()
        wx.CallAfter(UpdateDysnomia,G2frame,data)

    def OnZmult(event):
        DysData['Lagrange'][0] = Zmult.GetValue()
        wx.CallAfter(UpdateDysnomia,G2frame,data)

    def OnStart(event):
        DysData['DenStart'] = Start.GetValue()

    def OnPrior(event):
        DysData['prior'] = Prior.GetValue()
        if DysData['prior'] == 'last run':
            if os.path.isfile(pName+'_prior.pgrid'):
                os.remove(pName+'_prior.pgrid')
            os.rename(pName+'.pgrid',pName+'_prior.pgrid')

    def OnFileCheck(event):
        DysData['clear'] = fileCheck.GetValue()

    generalData = data['General']
    pName = generalData['Name'].replace(' ','_')
    Map = generalData['Map']
    UseList = Map['RefList']
    pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,UseList[0])       #only use 1st histogram
    if not pId:
        wx.MessageBox('You must prepare a fourier map before running Dysnomia','Dysnomia Error',
            style=wx.ICON_ERROR)
        return
    reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists'))
    reflData = reflSets[generalData['Name']]['RefList']
    refDmin = reflData[-1][4]
    mulMin = np.argmin(reflData[:][3])
    if reflData[mulMin][3] < 0:
        refDmin = reflData[mulMin-1][4]
    MEMData = G2frame.MEMData
    if MEMData.GetSizer():
        MEMData.GetSizer().Clear(True)
    DysData = data['Dysnomia']
    if 'overlap' not in DysData:
        DysData['overlap'] = 1.0
    if 'MEMdmin' not in DysData:
        DysData['MEMdmin'] = refDmin
    if 'clear' not in DysData:
        DysData['clear'] = True
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl = f"Maximum Entropy Method (Dysnomia) controls for {data['General']['Name']!r}"
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)

    mainSizer.Add(wx.StaticText(MEMData,label=
        ' For use of Dysnomia, please cite:\n'+
          G2G.GetCite('Dysnomia',wrap=60,indent=5)))
    lineSizer = wx.BoxSizer(wx.HORIZONTAL)
    lineSizer.Add(wx.StaticText(MEMData,label=' MEM Optimization method: '),0,WACV)
    OptMeth = wx.ComboBox(MEMData,-1,value=DysData['Optimize'],choices=['ZSPA','L-BFGS'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    OptMeth.Bind(wx.EVT_COMBOBOX,OnOptMeth)
    lineSizer.Add(OptMeth,0,WACV)
    lineSizer.Add(wx.StaticText(MEMData,label=' Peak overlap factor'),0,WACV)
    overlap = G2G.ValidatedTxtCtrl(MEMData,DysData,'overlap',nDig=(10,4),xmin=0.1,xmax=1.)
    lineSizer.Add(overlap,0,WACV)
    mainSizer.Add(lineSizer)
    if DysData['Optimize'] == 'ZSPA':
        Zsizer = wx.BoxSizer(wx.HORIZONTAL)
        Zsizer.Add(wx.StaticText(MEMData,label=' Initial Lagrangian multiplier: from '),0,WACV)
        Zmult = wx.ComboBox(MEMData,value=DysData['Lagrange'][0],choices=['user','Dysnomia'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Zmult.Bind(wx.EVT_COMBOBOX,OnZmult)
        Zsizer.Add(Zmult,0,WACV)
        if DysData['Lagrange'][0] == 'user':
            Zsizer.Add(wx.StaticText(MEMData,label=' value: '),0,WACV)
            lamb = G2G.ValidatedTxtCtrl(MEMData,DysData['Lagrange'],1,nDig=(10,4),xmin=0.0001,xmax=1.)
            Zsizer.Add(lamb,0,WACV)
        Zsizer.Add(wx.StaticText(MEMData,label=' Adjust by: '),0,WACV)
        dlamb = G2G.ValidatedTxtCtrl(MEMData,DysData['Lagrange'],2,nDig=(8,2),xmin=0.05,xmax=0.1)
        Zsizer.Add(dlamb,0,WACV)
        mainSizer.Add(Zsizer)

    Esizer = wx.BoxSizer(wx.HORIZONTAL)
    Esizer.Add(wx.StaticText(MEMData,label=' Weight by d-spacing**'),0,WACV)
    Efact = G2G.ValidatedTxtCtrl(MEMData,DysData,'wt pwr',xmin=0,xmax=4,size=(50,20))
    Esizer.Add(Efact,0,WACV)
    Dmin = G2G.ValidatedTxtCtrl(MEMData,DysData,'MEMdmin',xmin=0.5,xmax=refDmin,size=(50,20))
    Esizer.Add(wx.StaticText(MEMData,label=' Minimum d-spacing for generated reflections: '),0,WACV)
    Esizer.Add(Dmin,0,WACV)
    mainSizer.Add(Esizer)

    if os.path.isfile(pName+'.pgrid'):
        PriorSizer = wx.BoxSizer(wx.HORIZONTAL)
        PriorSizer.Add(wx.StaticText(MEMData,label=' Start from densities: '),0,WACV)
        Start = wx.ComboBox(MEMData,-1,value=DysData['DenStart'],choices=['uniform','last run'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Start.Bind(wx.EVT_COMBOBOX,OnStart)
        PriorSizer.Add(Start,0,WACV)
        PriorSizer.Add(wx.StaticText(MEMData,label=' Use as prior: '),0,WACV)
        Prior = wx.ComboBox(MEMData,-1,value=DysData['prior'],choices=['uniform','last run'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Prior.Bind(wx.EVT_COMBOBOX,OnPrior)
        PriorSizer.Add(Prior,0,WACV)
        mainSizer.Add(PriorSizer)
    else:
        DysData['DenStart'] = 'uniform'
        DysData['prior'] = 'uniform'

    Csizer = wx.BoxSizer(wx.HORIZONTAL)
    Csizer.Add(wx.StaticText(MEMData,label=' Maximum number of cycles: '),0,WACV)
    Cyc = G2G.ValidatedTxtCtrl(MEMData,DysData,'Ncyc',xmin=0,xmax=10000,size=(50,20))
    Csizer.Add(Cyc,0,WACV)
    fileCheck = wx.CheckBox(MEMData,label='Clear Dynsomia files? ')
    fileCheck.SetValue(DysData['clear'])
    fileCheck.Bind(wx.EVT_CHECKBOX,OnFileCheck)
    Csizer.Add(fileCheck,0,WACV)
    mainSizer.Add(Csizer)
    G2phsG.SetPhaseWindow(G2frame.MEMData,mainSizer)

def OnLoadDysnomia(event,G2frame,data):
    print('Load MEM - might not be implemented')

def OnSaveDysnomia(event,G2frame,data):
    print('Save MEM - might not be implemented')

def OnRunDysnomia(event,G2frame,data):

    path2GSAS2 = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
    DYSNOMIA = os.path.join(path2GSAS2,'Dysnomia','Dysnomia64.exe')
    DysData = data['Dysnomia']

    if not os.path.exists(DYSNOMIA):
        wx.MessageBox(''' Dysnomia is not installed. Please download it from
https://jp-minerals.org/dysnomia/en/
and install it at.'''+DYSNOMIA,
            caption='Dysnomia not installed',style=wx.ICON_ERROR)
        return

    generalData = data['General']
    Map = generalData['Map']
    UseList = Map['RefList']
    pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,UseList[0])       #only use 1st histogram
    if not pId:
        wx.MessageBox('You must prepare a Fourier map before running Dysnomia','Dysnomia Error',
            style=wx.ICON_ERROR)
        return
    reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists'))
    reflData = reflSets[generalData['Name']]['RefList']
    if 'Type' not in Map:
        wx.MessageBox('You must prepare a Fourier map before running Dysnomia','Dysnomia Error',
            style=wx.ICON_ERROR)
        return
    Type = Map['Type']
    MEMtype = 0
    if 'N' in Type:
        for el in generalData['Isotope']:
            isotope = generalData['Isotope'][el]
            if el not in generalData['Isotopes']:
                continue
            if generalData['Isotopes'][el][isotope]['SL'][0] < 0.:
                MEMtype = 1
    prfName = str(G2pwd.makePRFfile(data,MEMtype))
    if not G2pwd.makeMEMfile(data,reflData,MEMtype,DYSNOMIA):
        SpGrp = generalData['SGData']['SpGrp']
        wx.MessageBox('Non standard space group '+SpGrp+' not permitted in Dysnomia','Dysnomia Error',
            style=wx.ICON_ERROR)
        return
    wx.MessageBox(' For use of Dysnomia, please cite:\n\n'+
                      G2G.GetCite('Dysnomia'),
                      caption='Dysnomia (MEM)',style=wx.ICON_INFORMATION)

    print('Run '+DYSNOMIA)
    subp.call([DYSNOMIA,prfName])

    DysData['DenStart'] = 'uniform'
    DysData['prior'] = 'uniform'
    wx.CallAfter(UpdateDysnomia,data)

    goon,reflData = G2pwd.MEMupdateReflData(prfName,data,reflData)
    if goon:
        reflSets[generalData['Name']]['RefList'] = reflData
        G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists'),reflSets)
        G2phsG.DoFourierMaps(G2frame,data)           #auto run Fourier
        if DysData['clear']:
            os.remove(os.path.splitext(prfName)[0]+'.fba')
            os.remove(os.path.splitext(prfName)[0]+'.mem')
            os.remove(os.path.splitext(prfName)[0]+'.out')
            os.remove(os.path.splitext(prfName)[0]+'.prf')
            os.remove(os.path.splitext(prfName)[0]+'_eps.raw')
    else:
        wx.MessageBox('Dysnomia failed to make new structure factors','Dysnomia Error',
            style=wx.ICON_ERROR)

# #### RMC Data page ################################################################################
# # fullrmc stuff TODO:
# #  1) need to implement swapping in scripts
# #  2) fullrmc tutorials

# def UpdateRMC(event=None):
#     ''' Present the controls for running fullrmc, RMCProfile or PDFfit
#     '''
#     global runFile
#     def OnRMCselect(event):
#         G2frame.RMCchoice = RMCsel.GetStringSelection()
#         wx.CallLater(200,UpdateRMC)

#     def GetAtmChoice(pnl,RMCPdict):

#         Indx = {}
#         def OnAtSel(event):
#             Obj = event.GetEventObject()
#             itype = Indx[Obj.GetId()]
#             tid = RMCPdict['atSeq'].index(Obj.GetStringSelection())
#             if itype < nTypes:
#                 if itype == tid:
#                     tid += 1
#                 RMCPdict['atSeq'] = G2lat.SwapItems(RMCPdict['atSeq'],itype,tid)
#             Pairs= []
#             atSeq = RMCPdict['atSeq']
#             lenA = len(atSeq)
#             for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA)] for i in range(lenA)]:
# #                for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
#                 Pairs += pair
#             RMCPdict['Pairs'] = {pairs:[0.0,0.0,0.0] for pairs in Pairs}
#             if RMCPdict['useBVS']:
#                 BVSpairs = []
#                 for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i+1,lenA)] for i in range(lenA)]:
#                     BVSpairs += pair
#                 RMCPdict['BVS'] = {pairs:[0.0,0.0,0.0,0.0] for pairs in BVSpairs}
#             wx.CallAfter(UpdateRMC)

#         def OnValSel(event):
#             Obj = event.GetEventObject()
#             itype = Indx[Obj.GetId()]
#             RMCPdict['Oxid'][itype][0] = Obj.GetStringSelection()
#             wx.CallAfter(UpdateRMC)

#         nTypes = len(RMCPdict['aTypes'])
#         atmChoice = wx.FlexGridSizer(nTypes+1,5,5)
#         atmChoice.Add(wx.StaticText(pnl,label='atom ordering: '),0,WACV)
#         for iType in range(nTypes):
#             atChoice = RMCPdict['atSeq'][iType:]
#             atmSel = wx.ComboBox(pnl,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#             atmSel.SetStringSelection(RMCPdict['atSeq'][iType])
#             atmSel.Bind(wx.EVT_COMBOBOX,OnAtSel)
#             Indx[atmSel.GetId()] = iType
#             atmChoice.Add(atmSel,0,WACV)
#         if RMCPdict['useBVS']:
#             atmChoice.Add(wx.StaticText(pnl,label='Valence: '),0,WACV)
#             for itype in range(nTypes):
#                 valChoice = atmdata.BVSoxid[RMCPdict['atSeq'][itype]]
#                 valSel = wx.ComboBox(pnl,choices=valChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                 try:
#                     valSel.SetStringSelection(RMCPdict['Oxid'][itype][0])
#                 except IndexError:
#                     RMCPdict['Oxid'].append([RMCPdict['atSeq'][itype],0.0])
#                 valSel.Bind(wx.EVT_COMBOBOX,OnValSel)
#                 Indx[valSel.GetId()] = itype
#                 atmChoice.Add(valSel,0,WACV)
#             atmChoice.Add(wx.StaticText(pnl,label='BVS weight: '),0,WACV)
#             for itype in range(nTypes):
#                 atmChoice.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['Oxid'][itype],1,xmin=0.),0,WACV)
#         if G2frame.RMCchoice == 'RMCProfile':
#             atmChoice.Add(wx.StaticText(pnl,label='max shift: '),0,WACV)
#             for iType in range(nTypes):
#                 atId = RMCPdict['atSeq'][iType]
#                 atmChoice.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['aTypes'],atId,xmin=0.,xmax=1.),0,WACV)
#             atmChoice.Add(wx.StaticText(pnl,label='Isotope: '),0,WACV)
#             for iType in range(nTypes):
#                 atId = RMCPdict['atSeq'][iType]
#                 try:
#                     lbl = RMCPdict['Isotope'][atId]
#                 except:
#                     lbl = '?'
#                 atmChoice.Add(wx.StaticText(pnl,label=lbl),0,WACV)
#         return atmChoice

#     def GetSwapSizer(RMCPdict):

#         def OnDelSwap(event):
#             Obj = event.GetEventObject()
#             swap = Indx[Obj.GetId()]
#             del RMCPdict['Swaps'][swap]
#             wx.CallAfter(UpdateRMC)

#         Indx = {}
#         atChoice = RMCPdict['atSeq']
#         # if G2frame.RMCchoice == 'fullrmc':
#         #     atChoice = atNames
#         swapSizer = wx.FlexGridSizer(6,5,5)
#         swapLabels = [' ','Atom-A','Atom-B',' Swap prob.',' ','delete']
#         for lab in swapLabels:
#             swapSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#         for ifx,swap in enumerate(RMCPdict['Swaps']):
#             swapSizer.Add((20,-1))
#             for i in [0,1]:
#                 if swap[i] not in atChoice: swap[i] = atChoice[0]
#                 atmSel = G2G.EnumSelector(G2frame.FRMC,swap,i,atChoice)
#                 swapSizer.Add(atmSel,0,WACV)
#             swapSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,swap,2,xmin=0.01,xmax=0.5,size=(50,25)),0,WACV)
#             swapSizer.Add((20,-1))
#             delBtn = wx.Button(G2frame.FRMC,label='Del',style=wx.BU_EXACTFIT)
#             delBtn.Bind(wx.EVT_BUTTON,OnDelSwap)
#             Indx[delBtn.GetId()] = ifx
#             swapSizer.Add(delBtn,0,WACV)
#         return swapSizer

#     def GetPairSizer(pnl,RMCPdict):
#         pairSizer = wx.FlexGridSizer(len(RMCPdict['Pairs'])+1,5,5)
#         pairSizer.Add((5,5),0)
#         for pair in RMCPdict['Pairs']:
#             pairSizer.Add(wx.StaticText(pnl,label=pair),0,WACV)
#         if G2frame.RMCchoice == 'RMCProfile':
#             pairSizer.Add(wx.StaticText(pnl,label='%14s'%' Hard min: '),0,WACV)
#             for pair in RMCPdict['Pairs']:
#                 pairSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['Pairs'][pair],0,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#             pairSizer.Add(wx.StaticText(pnl,label='%14s'%' Search from: '),0,WACV)
#         elif G2frame.RMCchoice == 'fullrmc':
#             pairSizer.Add(wx.StaticText(pnl,label='%14s'%' Distance min: '),0,WACV)
#         for pair in RMCPdict['Pairs']:
#             pairSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['Pairs'][pair],
#                 1,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#         pairSizer.Add(wx.StaticText(pnl,label='%14s'%'to: '),0,WACV)
#         for pair in RMCPdict['Pairs']:
#             pairSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['Pairs'][pair],2,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#         return pairSizer

#     def GetMetaSizer(RMCPdict,metalist):
#         metaSizer = wx.FlexGridSizer(0,2,5,5)
#         for item in metalist:
#             metaSizer.Add(wx.StaticText(G2frame.FRMC,label=' Metadata item: '+item+' '),0,WACV)
#             metaSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['metadata'],item),0,WACV)
#         return metaSizer

#     def SetRestart(invalid,value,tc):
#         RMCPdict['ReStart'] = [True,True]

#     def GetSuperSizer(RMCPdict,Xmax):
#        superSizer = wx.BoxSizer(wx.HORIZONTAL)
#        axes = ['X','Y','Z']
#        for i,ax in enumerate(axes):
#            superSizer.Add(wx.StaticText(G2frame.FRMC,label=' %s-axis: '%ax),0,WACV)
#            superSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['SuperCell'],
#                i,xmin=1,xmax=Xmax,size=(50,25),OnLeave=SetRestart),0,WACV)
#        return superSizer

#     def FileSizer(RMCPdict):

#         def OnFileSel(event):
#             Obj = event.GetEventObject()
#             fil = Indx[Obj.GetId()]
#             G2frame.OnFileSave(event)
#             dlg = wx.FileDialog(G2frame.FRMC, 'Choose '+fil,G2G.GetImportPath(G2frame),
#                 style=wx.FD_OPEN ,wildcard=fil+'(*.*)|*.*')
#             if dlg.ShowModal() == wx.ID_OK:
#                 fpath,fName = os.path.split(dlg.GetPath())
#                 if os.path.exists(fName): # is there a file by this name in the current directory?
#                     RMCPdict['files'][fil][0] = fName
#                 else: # nope, copy it
#                     # TODO: is G2frame.LastGPXdir the right choice here or
#                     #       do I want the current working directory (same?)
#                     shutil.copy(dlg.GetPath(), os.path.join(G2frame.LastGPXdir,fName))
#                 if not os.path.exists(fName): # sanity check
#                     print(f'Error: file {fName} not found in .gpx directory ({G2frame.LastGPXdir})')
#                     return
#                 G2frame.LastImportDir = fpath    #set so next file is found in same place
#                 dlg.Destroy()
#                 RMCPdict['ReStart'][0] = True
#                 if G2frame.RMCchoice == 'PDFfit':
#                     start = 0
#                     XY = np.empty((1,2))
#                     while XY.shape[0] == 1:
#                         try:
#                             XY = np.loadtxt(fName,skiprows=start)
#                         except ValueError:
#                             start += 1
#                     name = 'Ndata'
#                     if 'X' in fil:
#                         name = 'Xdata'
#                     RMCPdict[name]['Datarange'][0] = np.min(XY.T[0])
#                     RMCPdict[name]['Datarange'][1] = np.max(XY.T[0])
#                     RMCPdict[name]['Fitrange'][1] = np.max(XY.T[0])
#             else:
#                 dlg.Destroy()

#             wx.CallAfter(UpdateRMC)

#         def OnFileFormat(event):
#             Obj = event.GetEventObject()
#             fil = Indx[Obj.GetId()]
#             RMCPdict['files'][fil][3] = Obj.GetStringSelection()

#         def OnPlotBtn(event):
#             Obj = event.GetEventObject()
#             fil = Indx[Obj.GetId()]
#             fileItem = RMCPdict['files'][fil]
#             start = 0
#             XY = np.empty((1,2))
#             while XY.shape[0] == 1:
#                 try:
#                     XY = np.loadtxt(fileItem[0],skiprows=start)
#                 except ValueError:
#                     start += 1
#                     if start > 500:     #absurd number of header lines!
#                         wx.MessageBox('WARNING: %s has bad data at end;\n RMCProfile may fail to read it'%fileItem[0],
#                             style=wx.ICON_ERROR)
#                         break
#             Xlab = 'Q'
#             if 'G(R)' in fileItem[2].upper():
#                 Xlab = 'R'
#             G2plt.PlotXY(G2frame,[XY.T[:2],],labelX=Xlab,
#                 labelY=fileItem[2],newPlot=True,Title=fileItem[0],
#                 lines=True)

#         def OnCorrChk(event):
#             Obj = event.GetEventObject()
#             fil = Indx[Obj.GetId()]
#             RMCPdict['files'][fil][3] = not RMCPdict['files'][fil][3]

#         def OnDelBtn(event):
#             Obj = event.GetEventObject()
#             fil = Indx[Obj.GetId()]
#             RMCPdict['files'][fil][0] = 'Select'
#             RMCPdict['ReStart'][0] = True
#             wx.CallAfter(UpdateRMC)

#         def OnRef(event):
#             Obj = event.GetEventObject()
#             name,item = Indx[Obj.GetId()]
#             RMCPdict[name][item][1] = not RMCPdict[name][item][1]

#         def OnRefSel(event):
#             RMCPdict['refinement'] = reftype.GetStringSelection()
#             wx.CallLater(100,UpdateRMC)

#         def OnDataSel(event):
#             RMCPdict['SeqDataType'] = dataType.GetStringSelection()

#         def OnSeqCopy(event):
#             RMCPdict['SeqCopy'] = not RMCPdict['SeqCopy']

#         def OnSeqReverse(event):
#             RMCPdict['SeqReverse'] = not RMCPdict['SeqReverse']

#         # --- FileSizer starts here
#         Indx = {}
#         mainSizer = wx.BoxSizer(wx.VERTICAL)
#         if G2frame.RMCchoice == 'PDFfit':
#             topSizer = wx.BoxSizer(wx.HORIZONTAL)
#             reftype = wx.RadioBox(G2frame.FRMC,label='PDFfit refinement type:',choices=['normal','sequential'])
#             reftype.SetStringSelection(RMCPdict.get('refinement','normal'))
#             reftype.Bind(wx.EVT_RADIOBOX,OnRefSel)
#             topSizer.Add(reftype)
#             if 'seq' in RMCPdict.get('refinement','normal'):
#                 dataType = wx.RadioBox(G2frame.FRMC,label='Seq data type:',choices=['X','N'])
#                 dataType.SetStringSelection(RMCPdict.get('SeqDataType','X'))
#                 dataType.Bind(wx.EVT_RADIOBOX,OnDataSel)
#                 topSizer.Add(dataType)
#                 endSizer = wx.BoxSizer(wx.VERTICAL)
#                 seqcopy = wx.CheckBox(G2frame.FRMC,label=' Copy to next')
#                 seqcopy.SetValue(RMCPdict['SeqCopy'])
#                 seqcopy.Bind(wx.EVT_CHECKBOX,OnSeqCopy)
#                 endSizer.Add(seqcopy)
#                 seqreverse = wx.CheckBox(G2frame.FRMC,label=' Reverse processing')
#                 seqreverse.SetValue(RMCPdict['SeqReverse'])
#                 seqreverse.Bind(wx.EVT_CHECKBOX,OnSeqReverse)
#                 endSizer.Add(seqreverse)
#                 topSizer.Add(endSizer,0,WACV)
#             mainSizer.Add(topSizer)
#         elif G2frame.RMCchoice == 'fullrmc':
#             topSizer = wx.BoxSizer(wx.HORIZONTAL)
#             topSizer.Add(wx.StaticText(G2frame.FRMC,label='  Select data for processing (files must be 2 columns w/headers preceeded by "#"; edit if needed)'))
#             mainSizer.Add(topSizer)
#             Heads = ['Name','File','type','Plot','Delete']
#             fileSizer = wx.FlexGridSizer(5,5,5)
#             Formats = ['RMC','GUDRUN','STOG']
#             for head in Heads:
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=head),0,WACV)
#             for fil in RMCPdict['files']:
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=fil),0,WACV)
#                 Rfile = RMCPdict['files'][fil][0]
#                 filSel = wx.Button(G2frame.FRMC,label=Rfile)
#                 filSel.Bind(wx.EVT_BUTTON,OnFileSel)
#                 Indx[filSel.GetId()] = fil
#                 fileSizer.Add(filSel,0,WACV)
#                 if Rfile and os.path.exists(Rfile): # in case .gpx file is moved away from G(R), F(Q), etc. files
#                     #fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['files'][fil],1,size=(50,25)),0,WACV)
#                     #patch
#                     if len(RMCPdict['files'][fil]) < 4:
#                         RMCPdict['files'][fil].append(0)
#                     if len(RMCPdict['files'][fil]) < 5:
#                         RMCPdict['files'][fil].append(True)
#                     #end patch
#                     if 'G(r)' in fil:
#                         choices = 'G(r)-RMCProfile','G(r)-PDFFIT','g(r)'
#                         if type(RMCPdict['files'][fil][3]) is bool: RMCPdict['files'][fil][3] = 0
#                         fmtTyp = G2G.G2ChoiceButton(G2frame.FRMC,choices,RMCPdict['files'][fil],3)
#                     elif '(Q)' in fil:
#                         choices = 'F(Q)-RMCProfile','S(Q)-PDFFIT'
#                         if type(RMCPdict['files'][fil][3]) is bool: RMCPdict['files'][fil][3] = 0
#                         fmtTyp = G2G.G2ChoiceButton(G2frame.FRMC,choices,RMCPdict['files'][fil],3)
#                     else:
#                         fmtTyp = (-1,-1)
#                     fileSizer.Add(fmtTyp,0,WACV)
#                     plotBtn = wx.Button(G2frame.FRMC,label='Plot',style=wx.BU_EXACTFIT)
#                     plotBtn.Bind(wx.EVT_BUTTON,OnPlotBtn)
#                     Indx[plotBtn.GetId()] = fil
#                     fileSizer.Add(plotBtn,0,WACV)
#                     delBtn = wx.Button(G2frame.FRMC,label='Del',style=wx.BU_EXACTFIT)
#                     delBtn.Bind(wx.EVT_BUTTON,OnDelBtn)
#                     Indx[delBtn.GetId()] = fil
#                     fileSizer.Add(delBtn,0,WACV)
#                     if '(Q)' in fil:
#                         fileSizer.Add((-1,-1),0)
#                         corrChk = wx.CheckBox(G2frame.FRMC,label='Apply sinc convolution? ')
#                         corrChk.SetValue(RMCPdict['files'][fil][4])
#                         Indx[corrChk.GetId()] = fil
#                         corrChk.Bind(wx.EVT_CHECKBOX,OnCorrChk)
#                         fileSizer.Add(corrChk,0,WACV)
#                         #fileSizer.Add((-1,-1),0)
#                         fileSizer.Add((-1,-1),0)
#                         fileSizer.Add((-1,-1),0)
#                         fileSizer.Add((-1,-1),0)
#                 elif 'Select' not in Rfile: # file specified, but must not exist
#                     RMCPdict['files'][fil][0] = 'Select' # set filSel?
#                     fileSizer.Add(wx.StaticText(G2frame.FRMC,
#                             label='Warning: file not found.\nWill be removed'),0)
#                     fileSizer.Add((-1,-1),0)
#                     fileSizer.Add((-1,-1),0)
#                 else:
#                     RMCPdict['files'][fil][0] = 'Select' # set filSel?
#                     #fileSizer.Add((-1,-1),0)
#                     fileSizer.Add((-1,-1),0)
#                     fileSizer.Add((-1,-1),0)
#                     fileSizer.Add((-1,-1),0)
#             mainSizer.Add(fileSizer,0)
#             return mainSizer

#         if G2frame.RMCchoice == 'PDFfit' and RMCPdict['refinement'] == 'sequential':

#             def OnAddPDF(event):
#                 ''' Add PDF G(r)s while maintanining original sequence
#                 '''
#                 usedList = RMCPdict['seqfiles']
#                 PDFlist = [item[1:][0] for item in G2frame.GetFileList('PDF')]
#                 PDFdict = dict([item[1:] for item in G2frame.GetFileList('PDF')])
#                 PDFnames = [item for item in PDFdict if item not in [itm[0] for itm in usedList]]
#                 dlg = G2G.G2MultiChoiceDialog(G2frame.FRMC,'Add PDF dataset',
#                     'Select G(r) data to use in seq. PDFfit',PDFnames)
#                 if dlg.ShowModal() == wx.ID_OK:
#                     PDFuse = dlg.GetSelections()
#                     for item in PDFuse:
#                         pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PDFnames[item])
#                         data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'PDF Controls'))
#                         try:
#                             insrt = PDFlist.index(PDFnames[item])-1
#                             RMCPdict['seqfiles'].insert(insrt+1,[PDFnames[item],data])
#                         except ValueError:
#                             RMCPdict['seqfiles'].append([PDFnames[item],data])
#                 dlg.Destroy()
#                 wx.CallAfter(UpdateRMC)

#             def OnDelPDF(event):
#                 usedList = [item[0] for item in RMCPdict['seqfiles']]
#                 dlg = G2G.G2MultiChoiceDialog(G2frame.FRMC,'Delete PDF dataset',
#                     'Select G(r) data to delete frpm seq. PDFfit',usedList)
#                 if dlg.ShowModal() == wx.ID_OK:
#                     PDFdel = dlg.GetSelections()
#                     PDFdel.reverse()
#                     for item in PDFdel:
#                         del RMCPdict['seqfiles'][item]
#                 dlg.Destroy()
#                 wx.CallAfter(UpdateRMC)

#             def OnSetColVal(event):
#                 parms = {'Rmin':[0.01,5.0],'Rmax':[5.,30.],'dscale':[0.5,2.0],
#                     'qdamp':[0.0,0.5],'qbroad':[0.0,0.1],'Temp':300}
#                 c =  event.GetCol()
#                 if c >= 0:
#                     if c in [3,5,7]:
#                         seqGrid.ClearSelection()
#                         seqGrid.SelectCol(c,True)
#                         if seqGrid.GetColLabelValue(c) != 'refine': return
#                         choice = ['Y - vary all','N - vary none',]
#                         dlg = wx.SingleChoiceDialog(G2frame,'Select refinement option for '+seqGrid.GetColLabelValue(c-1),
#                             'Refinement controls',choice)
#                         dlg.CenterOnParent()
#                         if dlg.ShowModal() == wx.ID_OK:
#                             sel = dlg.GetSelection()
#                             varib = colLabels[c-1]
#                             if sel == 0:
#                                 for row in range(seqGrid.GetNumberRows()): RMCPdict['seqfiles'][row][1][varib][1]=True
#                             else:
#                                 for row in range(seqGrid.GetNumberRows()): RMCPdict['seqfiles'][row][1][varib][1]=False
#                     elif c in [0,1,2,4,6,8]:
#                         seqGrid.ClearSelection()
#                         seqGrid.SelectCol(c,True)
#                         parm = colLabels[c]
#                         dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter value for '+parm,0.0,parms[parm])
#                         if dlg.ShowModal() == wx.ID_OK:
#                             value = dlg.GetValue()
#                             if c in [2,4,6]:
#                                 for row in range(seqGrid.GetNumberRows()): RMCPdict['seqfiles'][row][1][parm][0] = value
#                             elif c == 8:
#                                 for row in range(seqGrid.GetNumberRows()): RMCPdict['seqfiles'][row][1][parm] = value
#                             else:
#                                 for row in range(seqGrid.GetNumberRows()): RMCPdict['seqfiles'][row][1]['Fitrange'][c] = value
#                     wx.CallAfter(UpdateRMC)

#             def OnSetVal(event):
#                 r,c= event.GetRow(),event.GetCol()
#                 if c >= 0:
#                     if c in [3,5,7]:
#                         varib = colLabels[c-1]
#                         RMCPdict['seqfiles'][r][1][varib][1] = bool(seqGrid.GetCellValue(r,c))
#                     elif c in [0,1,2,4,6,8]:
#                         parm = colLabels[c]
#                         if c in [2,4,6]:
#                             RMCPdict['seqfiles'][r][1][parm][0] = float(seqGrid.GetCellValue(r,c))
#                         elif c == 8:
#                             RMCPdict['seqfiles'][r][1][parm] = float(seqGrid.GetCellValue(r,c))
#                         else:
#                             RMCPdict['seqfiles'][r][1]['Fitrange'][c] = float(seqGrid.GetCellValue(r,c))

#             topSizer = wx.BoxSizer(wx.HORIZONTAL)
#             topSizer.Add(wx.StaticText(G2frame.FRMC,label='  Select data for processing: '))
#             mainSizer.Add(topSizer)
#             G2frame.GetStatusBar().SetStatusText('NB: All PDFs used in sequential PDFfit must be the same type ("X" or "N") - there is no check',1)
#             if 'seqfiles' not in RMCPdict:
#                 RMCPdict['seqfiles'] = []
#             topSizer = wx.BoxSizer(wx.HORIZONTAL)
#             topSizer.Add(wx.StaticText(G2frame.FRMC,label=' Sequential data list for PDFfit:  '),0,WACV)
#             addPDF = wx.Button(G2frame.FRMC,label='Add PDF G(r) data sets')
#             addPDF.Bind(wx.EVT_BUTTON,OnAddPDF)
#             topSizer.Add(addPDF,0,WACV)
#             delPDF = wx.Button(G2frame.FRMC,label='Delete PDF G(r) data sets')
#             delPDF.Bind(wx.EVT_BUTTON,OnDelPDF)
#             topSizer.Add(delPDF,0,WACV)
#             mainSizer.Add(topSizer)
#             table = [[item[1]['Fitrange'][0],item[1]['Fitrange'][1],
#                 item[1]['dscale'][0],item[1]['dscale'][1],item[1]['qdamp'][0],item[1]['qdamp'][1],
#                 item[1]['qbroad'][0],item[1]['qbroad'][1],item[1].get('Temp',300.)] for item in RMCPdict['seqfiles']]
#             colLabels = ['Rmin','Rmax','dscale','refine','qdamp','refine','qbroad','refine','Temp']
#             rowLabels = [item[0] for item in RMCPdict['seqfiles']]
#             Types = [wg.GRID_VALUE_FLOAT+':10,2',wg.GRID_VALUE_FLOAT+':10,2',
#                      wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
#                      wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,
#                      wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,2']
#             seqTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
#             seqGrid = G2G.GSGrid(G2frame.FRMC)
#             seqGrid.SetTable(seqTable, True)
#             seqGrid.AutoSizeColumns(True)
#             seqGrid.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, OnSetColVal)
#             seqGrid.Bind(wg.EVT_GRID_CELL_CHANGED, OnSetVal)
#             mainSizer.Add(seqGrid)
#             return mainSizer

# # begin FileSizer
#         topSizer = wx.BoxSizer(wx.HORIZONTAL)
#         topSizer.Add(wx.StaticText(G2frame.FRMC,label='  Select data for processing: '))
#         mainSizer.Add(topSizer)
#         # RMCProfile & PDFfit (Normal)
#         Heads = ['Name','File','Format','Weight','Plot','Delete']
#         fileSizer = wx.FlexGridSizer(6,5,5)
#         Formats = ['RMC','GUDRUN','STOG']
#         for head in Heads:
#             fileSizer.Add(wx.StaticText(G2frame.FRMC,label=head),0,WACV)
#         for fil in RMCPdict['files']:
#             for head in Heads:
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=20*'-'),0,WACV)
#             fileSizer.Add(wx.StaticText(G2frame.FRMC,label=fil),0,WACV)
#             Rfile = RMCPdict['files'][fil][0]
#             filSel = wx.Button(G2frame.FRMC,label=Rfile)
#             filSel.Bind(wx.EVT_BUTTON,OnFileSel)
#             Indx[filSel.GetId()] = fil
#             fileSizer.Add(filSel,0,WACV)
#             nform = 3
#             Name = 'Ndata'
#             if 'Xray' in fil:
#                 nform = 1
#                 Name = 'Xdata'
#             if Rfile and os.path.exists(Rfile): #incase .gpx file is moved away from G(R), F(Q), etc. files
#                 fileFormat = wx.ComboBox(G2frame.FRMC,choices=Formats[:nform],style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                 fileFormat.SetStringSelection(RMCPdict['files'][fil][3])
#                 Indx[fileFormat.GetId()] = fil
#                 fileFormat.Bind(wx.EVT_COMBOBOX,OnFileFormat)
#                 fileSizer.Add(fileFormat,0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['files'][fil],1),0,WACV)
#                 plotBtn = wx.Button(G2frame.FRMC,label='Plot?',style=wx.BU_EXACTFIT)
#                 plotBtn.Bind(wx.EVT_BUTTON,OnPlotBtn)
#                 Indx[plotBtn.GetId()] = fil
#                 fileSizer.Add(plotBtn,0,WACV)
#                 delBtn = wx.Button(G2frame.FRMC,label='Del',style=wx.BU_EXACTFIT)
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelBtn)
#                 Indx[delBtn.GetId()] = fil
#                 fileSizer.Add(delBtn,0,WACV)
#             else:
#                 RMCPdict['files'][fil][0] = 'Select'
#                 fileSizer.Add((5,5),0)
#                 fileSizer.Add((5,5),0)
#                 fileSizer.Add((5,5),0)
#                 fileSizer.Add((5,5),0)
#             if 'Select' not in Rfile and 'PDFfit' in G2frame.RMCchoice:
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=' R-range (from/to)'),0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[Name]['Fitrange'],0,xmin=RMCPdict[Name]['Datarange'][0],xmax=3.0),0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[Name]['Fitrange'],1,xmin=10.0,xmax=RMCPdict[Name]['Datarange'][1]),0,WACV)
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=' Scale factor: '),0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[Name]['dscale'],0,xmin=0.001,xmax=20.0),0,WACV)
#                 scaleref = wx.CheckBox(G2frame.FRMC,label='refine')
#                 scaleref.SetValue(RMCPdict[Name]['dscale'][1])
#                 Indx[scaleref.GetId()] = [Name,'dscale']
#                 scaleref.Bind(wx.EVT_CHECKBOX,OnRef)
#                 fileSizer.Add(scaleref,0,WACV)
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=' Qdamp '),0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[Name]['qdamp'],0,xmin=0.0,xmax=1.0),0,WACV)
#                 qdampref = wx.CheckBox(G2frame.FRMC,label='refine')
#                 qdampref.SetValue(RMCPdict[Name]['qdamp'][1])
#                 Indx[qdampref.GetId()] = [Name,'qdamp']
#                 qdampref.Bind(wx.EVT_CHECKBOX,OnRef)
#                 fileSizer.Add(qdampref,0,WACV)
#                 fileSizer.Add(wx.StaticText(G2frame.FRMC,label=' Qbroad '),0,WACV)
#                 fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[Name]['qbroad'],0,xmin=0.0,xmax=1.0),0,WACV)
#                 qbroadref = wx.CheckBox(G2frame.FRMC,label='refine')
#                 qbroadref.SetValue(RMCPdict[Name]['qbroad'][1])
#                 Indx[qbroadref.GetId()] = [Name,'qbroad']
#                 qbroadref.Bind(wx.EVT_CHECKBOX,OnRef)
#                 fileSizer.Add(qbroadref,0,WACV)

#         mainSizer.Add(fileSizer,0)

#         return mainSizer

#     def fullrmcSizer(RMCPdict):
#         mainSizer = wx.BoxSizer(wx.VERTICAL)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=
# '''* "Atomic Stochastic Modeling & Optimization with fullrmc", B. Aoun, J. Appl. Cryst. 2022, 55(6) 1664-1676,
#  DOI: 10.1107/S1600576722008536;
# * "Fullrmc, a Rigid Body Reverse Monte Carlo Modeling Package Enabled with Machine Learning and Artificial
#    Intelligence", B. Aoun, Jour. Comp. Chem. (2016), 37, 1102-1111. DOI: 10.1002/jcc.24304;
# * www.fullrmc.com
#  '''))
#         # if G2pwd.findfullrmc() is None:
#         #     mainSizer.Add(wx.StaticText(G2frame.FRMC,
#         #         label="\nsorry, fullrmc not installed or was not located"))
#         #     return mainSizer
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC,True)
#         #mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' fullrmc big box starting pdb file preparation:'),0)

#         # initialize fullrmc dictionary if needed
#         RMCPdict = data['RMC']['fullrmc'] = data['RMC'].get('fullrmc',{})
#         # update, if atoms list has been updated
#         Atypes = [atype.split('+')[0].split('-')[0] for atype in data['General']['AtomTypes']]
#         aTypes = dict(zip(Atypes,len(Atypes)*[0.10,]))
#         if len(data['RMC']['fullrmc'].get('aTypes',{})) != len(aTypes):
#             #print('atypes has changed')
#             atSeq = list(aTypes.keys())
#             lenA = len(atSeq)
#             Pairs= []
#             for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]]
#                              for i in range(lenA) if 'Va' not in atSeq[i]]:
#                 Pairs += pair
#             Pairs = {pairs:[0.0,0.0,0.0] for pairs in Pairs}
#             RMCPdict.update({'aTypes':aTypes,'atSeq':atSeq,'Pairs':Pairs})
#         RMCPdict['files'] = RMCPdict.get('files',
#                         {'Neutron real space data; G(r): ':['Select',1.,'G(r)',0,True],
#                         'Neutron reciprocal space data; S(Q)-1: ':['Select',1.,'F(Q)',0,True],
#                         'Xray real space data; G(r): ':['Select',1.,'G(r)',0,True],
#                         'Xray reciprocal space data; S(Q)-1: ':['Select',1.,'F(Q)',0,True]})
#         if 'moleculePdb' not in RMCPdict:
#             RMCPdict.update({'moleculePdb':'Select','targetDensity':1.0,'maxRecursion':10000})
#         if 'Angles' not in RMCPdict:
#             RMCPdict.update({'Angles':[],'Angle Weight':1.e-5,'Bond Weight':1.e-5,'Torsions':[],'Torsion Weight':1.e-5})
#         for key,val in {'SuperCell':[1,1,1],'Box':[10.,10.,10.],'ReStart':[False,False],'Cycles':1,
#                 'Swaps':[],'useBVS':False,'FitScale':False,'AveCN':[],'FxCN':[],
#                 'min Contact':1.5,'periodicBound':True}.items():
#             RMCPdict[key] = RMCPdict.get(key,val)

#         def GetSuperSizer():
#             def ShowRmax(*args,**kwargs):
#                 cell = data['General']['Cell'][1:7]
#                 bigcell = np.array(cell)*np.array(RMCPdict['SuperCell']+[1,1,1])
#                 bigG = G2lat.cell2Gmat(bigcell)[0]
#                 rmax = min([0.5/np.sqrt(G2lat.calc_rDsq2(H,bigG)) for H in np.eye(3)])
#                 rmaxlbl.SetLabel('  Rmax = {:.1f}'.format(rmax))
#             superSizer = wx.BoxSizer(wx.HORIZONTAL)
#             axes = ['X','Y','Z']
#             for i,ax in enumerate(axes):
#                 superSizer.Add(wx.StaticText(G2frame.FRMC,label=' %s-axis: '%ax),0,WACV)
#                 superSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['SuperCell'],
#                     i,xmin=1,xmax=20,size=(50,25),OnLeave=ShowRmax),0,WACV)
#             rmaxlbl = wx.StaticText(G2frame.FRMC,label=' Rmax=?')
#             superSizer.Add(rmaxlbl,0,WACV)
#             ShowRmax()
#             return superSizer

#         def GetBoxSizer():
#             boxSizer = wx.BoxSizer(wx.HORIZONTAL)
#             axes = ['X','Y','Z']
#             for i,ax in enumerate(axes):
#                 boxSizer.Add(wx.StaticText(G2frame.FRMC,label=' %s-axis: '%ax),0,WACV)
#                 boxSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['Box'],
#                     i,xmin=10.,xmax=50.,size=(50,25)),0,WACV)
#             return boxSizer

#         def OnReStart(event):
#             RMCPdict['ReStart'][0] = not RMCPdict['ReStart'][0]

#         def OnAddSwap(event):
#             RMCPdict['Swaps'].append(['','',0.0,])
#             wx.CallAfter(UpdateRMC)

#         def OnPdbButton(event):
#             dlg = wx.FileDialog(G2frame.FRMC, 'Choose molecule pdb file',G2frame.LastGPXdir,
#                 style=wx.FD_OPEN ,wildcard='PDB file(*.pdb)|*.pdb')
#             if dlg.ShowModal() == wx.ID_OK:
#                 fpath,fName = os.path.split(dlg.GetPath())
#                 RMCPdict['moleculePdb'] = fName
#                 pdbButton.SetLabel(fName)

#         def OnAddAngle(event):
#             RMCPdict['Angles'].append(['','','',0.,0.,0.,0.])
#             wx.CallAfter(UpdateRMC)

#         # def OnAddTorsion(event):
#         #     RMCPdict['Torsions'].append(['','','','',0.,0.,0.,0.,0.,0.])
#         #     wx.CallAfter(UpdateRMC)

#         def GetAngleSizer():

#             def OnDelAngle(event):
#                 Obj = event.GetEventObject()
#                 angle = Indx[Obj.GetId()]
#                 del RMCPdict['Angles'][angle]
#                 wx.CallAfter(UpdateRMC)

#             # def OnAngleAtSel(event):
#             #     Obj = event.GetEventObject()
#             #     angle,i = Indx[Obj.GetId()]
#             #     RMCPdict['Angles'][angle][i] = Obj.GetStringSelection()

#             def SetRestart1(invalid,value,tc):
#                 RMCPdict['ReStart'][1] = True

#             Indx = {}
#             atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#             angleSizer = wx.GridBagSizer(0,5)
#             fxcnLabels1 = [' ',' ','Central','',None,'angle restraint values (deg)',None,'search distance (A)']
#             fxcnLabels2 = [' ','Atom-A','Atom','Atom-C','min','max','from','to']
#             for i in range(8):
#                 if fxcnLabels1[i]:
#                     cspan=1
#                     coloff = 0
#                     if fxcnLabels1[i-1] is None:
#                         cspan=2
#                         coloff = 1
#                     angleSizer.Add(wx.StaticText(G2frame.FRMC,label=fxcnLabels1[i]),
#                         (0,i-coloff),(1,cspan))
#                 if fxcnLabels2[i]:
#                     angleSizer.Add(wx.StaticText(G2frame.FRMC,wx.ID_ANY,
#                         label=fxcnLabels2[i],style=wx.CENTER),(1,i))
#             row = 1
#             for ifx,angle in enumerate(RMCPdict['Angles']):
#                 row += 1
#                 angleSizer.Add((30,-1),(row,0))
#                 for i in range(3):
#                     if angle[i] not in atChoice: angle[i] = atChoice[0]
#                     atmSel = G2G.EnumSelector(G2frame.FRMC,angle,i,atChoice)
#                     angleSizer.Add(atmSel,(row,1+i))
#                 for i in range(4):
#                     if i == 0:
#                         xmin,xmax=0.,180.
#                     elif i == 2:
#                         xmin,xmax=0.1,6.
#                     angleSizer.Add(
#                         G2G.ValidatedTxtCtrl(G2frame.FRMC,angle,3+i,xmin=xmin,xmax=xmax,
#                             OnLeave=SetRestart1,size=(50,25)),(row,4+i))
#                 delBtn = wx.Button(G2frame.FRMC,label='Del',style=wx.BU_EXACTFIT)
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelAngle)
#                 Indx[delBtn.GetId()] = ifx
#                 angleSizer.Add(delBtn,(row,9))
#             return angleSizer

#         # def GetTorsionSizer():

#         #     def OnDelTorsion(event):
#         #         Obj = event.GetEventObject()
#         #         angle = Indx[Obj.GetId()]
#         #         del RMCPdict['Torsions'][angle]
#         #         wx.CallAfter(UpdateRMC)

#         #     def OnTorsionAtSel(event):
#         #         Obj = event.GetEventObject()
#         #         torsion,i = Indx[Obj.GetId()]
#         #         RMCPdict['Torsions'][torsion][i] = Obj.GetStringSelection()

#         #     def SetRestart1(invalid,value,tc):
#         #         RMCPdict['ReStart'][1] = True

#         #     Indx = {}
#         #     atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#         #     torsionSizer = wx.FlexGridSizer(11,5,5)
#         #     fxcnLabels = [' ','Atom-A','Atom-B','Atom-C','Atom-D',' min angle1',' max angle1',' min angle2',' max angle2',' min angle3',' max angle3']
#         #     for lab in fxcnLabels:
#         #         torsionSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#         #     for ifx,torsion in enumerate(RMCPdict['Torsions']):
#         #         delBtn = wx.Button(G2frame.FRMC,label='Delete')
#         #         delBtn.Bind(wx.EVT_BUTTON,OnDelTorsion)
#         #         Indx[delBtn.GetId()] = ifx
#         #         torsionSizer.Add(delBtn,0,WACV)
#         #         for i in [0,1,2,3]:
#         #             atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#         #             atmSel.SetStringSelection(torsion[i])
#         #             atmSel.Bind(wx.EVT_COMBOBOX,OnTorsionAtSel)
#         #             Indx[atmSel.GetId()] = [ifx,i]
#         #             torsionSizer.Add(atmSel,0,WACV)
#         #         for i in  [4,5,6,7,8,9]:
#         #             torsionSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,torsion,i,xmin=0.,xmax=360.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#         #     return torsionSizer

#         generalData = data['General']
#         cx,ct,cs,cia = generalData['AtomPtrs']
#         #atomData = data['Atoms']
#         # atNames = [atom[ct-1] for atom in atomData]
#         # ifP1 = False
#         # if generalData['SGData']['SpGrp'] == 'P 1':
#         #     ifP1 = True
#         ifBox = False
#         if 'macromolecular' in generalData['Type']:
#             ifBox = True
#         lineSizer = wx.BoxSizer(wx.HORIZONTAL)
#         if ifBox:
#             lineSizer.Add(wx.StaticText(G2frame.FRMC,label=' Big box dimensions, %s:'%Angstr),0,WACV)
#             lineSizer.Add(GetBoxSizer(),0,WACV)
# #            elif ifP1:
#         lineSizer.Add(wx.StaticText(G2frame.FRMC,label=' Lattice multipliers:'),0,WACV)
#         lineSizer.Add(GetSuperSizer(),0,WACV)
#         lineSizer.Add((5,-1))
#         # Bachir suggests that w/o periodic boundaries, users are likely to use fullrmc wrong
#         #lineSizer.Add(G2G.G2CheckBox(G2frame.FRMC,'Impose periodic boundaries',RMCPdict,'periodicBound'),
#         #                  0,WACV)
#         mainSizer.Add(lineSizer,0)
#         if ifBox:
#             molecSizer = wx.BoxSizer(wx.HORIZONTAL)
#             molecSizer.Add(wx.StaticText(G2frame.FRMC,label=' Source molecule file '),0,WACV)
#             pdbButton = wx.Button(G2frame.FRMC,label=RMCPdict['moleculePdb'])
#             pdbButton.Bind(wx.EVT_BUTTON,OnPdbButton)
#             molecSizer.Add(pdbButton,0,WACV)
#             molecSizer.Add(wx.StaticText(G2frame.FRMC,label=' target density, gm/cc '),0,WACV)
#             molecSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'targetDensity',xmin=0.1,size=[60,25]),0,WACV)
#             molecSizer.Add(wx.StaticText(G2frame.FRMC,label=' max tries '),0,WACV)
#             molecSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'maxRecursion',xmin=1000,xmax=1000000,size=[60,25]),0,WACV)
#             mainSizer.Add(molecSizer,0)
#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' fullrmc run file preparation:'))
#         resLine = wx.BoxSizer(wx.HORIZONTAL)
#         resLine.Add(wx.StaticText(G2frame.FRMC,label=' Run '),0,WACV)
#         resLine.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Cycles',xmin=1,size=[60,25]))
#         resLine.Add(wx.StaticText(G2frame.FRMC,
#                     label=' computation cycles of '),0,WACV)
#         RMCPdict['Steps/cycle'] = RMCPdict.get('Steps/cycle',5000)
#         resLine.Add(G2G.EnumSelector(G2frame.FRMC,RMCPdict,'Steps/cycle',
#                     ['1K','5K','10K','50K'],[1000,5000,10000,50000]),0,WACV)
#         resLine.Add(wx.StaticText(G2frame.FRMC,
#                     label=' steps per cycle'),0,WACV)
#         mainSizer.Add(resLine,0)
#         resLine = wx.BoxSizer(wx.HORIZONTAL)
#         resLine.Add(wx.StaticText(G2frame.FRMC,
#                     label=' Restart fullrmc Engine? '),0,WACV)
#         restart = wx.CheckBox(G2frame.FRMC,label='(will clear old result!) ')
#         resLine.Add(restart,0,WACV)

#         restart.SetValue(RMCPdict['ReStart'][0])
#         restart.Bind(wx.EVT_CHECKBOX,OnReStart)
#         mainSizer.Add(resLine,0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(GetAtmChoice(G2frame.FRMC,RMCPdict),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         swapBox = wx.BoxSizer(wx.HORIZONTAL)
#         swapBox.Add(wx.StaticText(G2frame.FRMC,label='Atom swap probabiities: '),0,WACV)
#         swapAdd = wx.Button(G2frame.FRMC,label='Add',style=wx.BU_EXACTFIT)
#         swapAdd.Bind(wx.EVT_BUTTON,OnAddSwap)
#         swapBox.Add(swapAdd,0,WACV)
#         mainSizer.Add(swapBox,0)
#         if len(RMCPdict['Swaps']):
#             mainSizer.Add(GetSwapSizer(RMCPdict),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label='Geometry constraints && restraints'),0)
#         distBox = wx.BoxSizer(wx.HORIZONTAL)
#         distBox.Add(wx.StaticText(G2frame.FRMC,label='Distance constraints'),0,WACV)
#         # weights removed for now
#         #distBox.Add(wx.StaticText(G2frame.FRMC,label=', distance weight:'),0,WACV)
#         #distBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Bond Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
#         distBox.Add(wx.StaticText(G2frame.FRMC,label=' min contact dist: '),0,WACV)
#         distBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'min Contact',xmin=0.,xmax=4.,size=(50,25)),0,WACV)

#         RMCPdict['useBondConstraints'] = RMCPdict.get('useBondConstraints',True)
#         distBox.Add(wx.StaticText(G2frame.FRMC,label='  Use bond constraints? '),0,WACV)
#         distBox.Add(G2G.G2CheckBox(G2frame.FRMC,'',RMCPdict,'useBondConstraints',OnChange=UpdateRMC),
#                         0,WACV)
#         mainSizer.Add(distBox,0)

#         if RMCPdict['useBondConstraints']:
#             mainSizer.Add(GetPairSizer(G2frame.FRMC,RMCPdict),0)
#             mainSizer.Add((-1,10))
#         angBox = wx.BoxSizer(wx.HORIZONTAL)
#         angBox.Add(wx.StaticText(G2frame.FRMC,label='A-B-C angle restraints'),0,WACV)
#         # weights removed for now
#         #angBox.Add(wx.StaticText(G2frame.FRMC,label=', angle weight:'),0,WACV)
#         #angBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Angle Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
#         angBox.Add((20,-1))
#         angAdd = wx.Button(G2frame.FRMC,label='Add',style=wx.BU_EXACTFIT)
#         angAdd.Bind(wx.EVT_BUTTON,OnAddAngle)
#         angBox.Add(angAdd,0,WACV)
#         mainSizer.Add(angBox,0)
#         if len(RMCPdict['Angles']):
#             mainSizer.Add(GetAngleSizer(),0)
#         RMCPdict['Groups'] = RMCPdict.get('Groups',[])
#         def OnAddGroup(event):
#             index = len(RMCPdict['Groups'])
#             RMCPdict['Groups'].append([])
#             GroupEditor(index)
#         def OnDelGroup(event):
#             index = event.EventObject.index
#             del RMCPdict['Groups'][index]
#             wx.CallAfter(UpdateRMC)
#         def OnEdtGroup(event):
#             index = event.EventObject.index
#             GroupEditor(index)
#         def GroupEditor(index):
#             cx,ct,cs,cia = data['General']['AtomPtrs']
#             atomlbs = [a[ct-1] for a in data['Atoms']]
#             dlg = G2G.G2MultiChoiceDialog(G2frame.FRMC,'Atom Selector',
#                         'Select atoms to include in group',atomlbs,selected=RMCPdict['Groups'][index])
#             if dlg.ShowModal() == wx.ID_OK:
#                 RMCPdict['Groups'][index] = dlg.GetSelections()
#             dlg.Destroy()
#             if len(RMCPdict['Groups'][index]) == 0:
#                 del RMCPdict['Groups'][index]
#             wx.CallAfter(UpdateRMC)
#         if len(RMCPdict['Groups']) == 0:
#             grpAdd = wx.Button(G2frame.FRMC,label='Define atom group',style=wx.BU_EXACTFIT)
#             grpAdd.Bind(wx.EVT_BUTTON,OnAddGroup)
#             mainSizer.Add(grpAdd,0)
#         else:
#             grpBox = wx.BoxSizer(wx.HORIZONTAL)
#             grpBox.Add(wx.StaticText(G2frame.FRMC,label='Atom Groups:  '),0,WACV)
#             grpAdd = wx.Button(G2frame.FRMC,label='Add group',style=wx.BU_EXACTFIT)
#             grpAdd.Bind(wx.EVT_BUTTON,OnAddGroup)
#             RMCPdict['GroupMode'] = RMCPdict.get('GroupMode',0)
#             grpBox.Add(grpAdd,0,WACV)
#             grpBox.Add(wx.StaticText(G2frame.FRMC,
#                         label='  Group refinement mode: '),0,WACV)
#             grpBox.Add(G2G.EnumSelector(G2frame.FRMC,RMCPdict,'GroupMode',
#                     ('Rotate & Translate','Rotate only','Translate only'),
#                     [0,1,2]),0,WACV)
#             mainSizer.Add(grpBox,0)
#             for i,g in enumerate(RMCPdict['Groups']):
#                 grpBox = wx.BoxSizer(wx.HORIZONTAL)
#                 grpBox.Add((20,-1))
#                 grpBox.Add(wx.StaticText(G2frame.FRMC,label='Group #'+str(i+1)),0,WACV)
#                 grpBox.Add((4,-1))
#                 grpdel = wx.Button(G2frame.FRMC,label='Del',style=wx.BU_EXACTFIT)
#                 grpdel.Bind(wx.EVT_BUTTON,OnDelGroup)
#                 grpdel.index = i
#                 grpBox.Add(grpdel,0,WACV)
#                 grpadd = wx.Button(G2frame.FRMC,label='Edit',style=wx.BU_EXACTFIT)
#                 grpadd.Bind(wx.EVT_BUTTON,OnEdtGroup)
#                 grpadd.index = i
#                 grpBox.Add(grpadd,0,WACV)
#                 msg = ' Contains atoms: '
#                 for i,n in enumerate(g):
#                     if i+1 == len(g):
#                         msg += ' && '
#                     elif i > 0:
#                         msg += ', '
#                     msg += str(i)
#                 grpBox.Add(wx.StaticText(G2frame.FRMC,label=msg),0,WACV)
#                 mainSizer.Add(grpBox,0)

#         RMCPdict['addThermalBroadening'] = RMCPdict.get('addThermalBroadening',False)
#         mainSizer.Add((-1,5))
#         distBox = wx.BoxSizer(wx.HORIZONTAL)
#         distBox.Add(wx.StaticText(G2frame.FRMC,label=' Add thermal broadening? '),0,WACV)
#         distBox.Add(G2G.G2CheckBox(G2frame.FRMC,'',RMCPdict,'addThermalBroadening',OnChange=UpdateRMC),
#                         0,WACV)
#         if RMCPdict['addThermalBroadening']:
#             distBox.Add((15,-1))
#             distBox.Add(wx.StaticText(G2frame.FRMC,label='Uiso equiv.'),0,WACV)
#             RMCPdict['ThermalU'] = RMCPdict.get('ThermalU',{})
#             for atm in RMCPdict['aTypes']:
#                 RMCPdict['ThermalU'][atm] = RMCPdict['ThermalU'].get(atm,0.005)
#                 distBox.Add(wx.StaticText(G2frame.FRMC,label='  '+atm+':'),0,WACV)
#                 distBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['ThermalU'],atm,
#                                             xmin=0.0001,xmax=0.25,size=(50,25)),0,WACV)
#         mainSizer.Add(distBox,0)
#         if RMCPdict['addThermalBroadening']: mainSizer.Add((-1,5))

#         # Torsions are difficult to implement. Need to be internal to a unit cell & named with fullrmc
#         # atom labels. Leave this out, at least for now.
#         # torBox = wx.BoxSizer(wx.HORIZONTAL)
#         # torAdd = wx.Button(G2frame.FRMC,label='Add')
#         # torAdd.Bind(wx.EVT_BUTTON,OnAddTorsion)
#         # torBox.Add(torAdd,0,WACV)
#         # torBox.Add(wx.StaticText(G2frame.FRMC,label=' A-B-C-D torsion angle restraints (intracell only), weight: '),0,WACV)
#         # torBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Torsion Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
#         # mainSizer.Add(torBox,0)
#         # if len(RMCPdict['Torsions']):
#         #     mainSizer.Add(GetTorsionSizer(),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(FileSizer(RMCPdict))
#         return mainSizer

#     def RMCProfileSizer(RMCPdict):

#         def CheckAtms(Atypes):
#             newAtm = False
#             for atm in Atypes:
#                 if atm not in data['RMC']['RMCProfile'].get('aTypes',{}):
#                     newAtm = True
#                     break
#             for atm in data['RMC']['RMCProfile'].get('aTypes',{}):
#                 if atm not in Atypes:
#                     newAtm = True
#                     break
#             return newAtm

#         mainSizer = wx.BoxSizer(wx.VERTICAL)
#         subSizer = wx.BoxSizer(wx.HORIZONTAL)
#         subSizer.Add((-1,-1),1,wx.EXPAND)
#         subSizer.Add(wx.StaticText(G2frame.FRMC,label='RMCProfile setup'),0,WACV)
#         subSizer.Add((-1,-1),1,wx.EXPAND)
#         mainSizer.Add(subSizer)
#         mainSizer.Add((5,5))
#         txt = wx.StaticText(G2frame.FRMC,label=
#                 f"Please cite: {G2G.GetCite('RMCProfile')}")
#         txt.Wrap(500)
#         mainSizer.Add(txt)
#         mainSizer.Add((5,5))
#         Atypes = [atype.split('+')[0].split('-')[0] for atype in data['General']['AtomTypes']]
#         aTypes = dict(zip(Atypes,len(Atypes)*[0.10,]))
#         atSeq = list(aTypes.keys())
#         lenA = len(atSeq)
#         atOxid = [[atmdata.BVSoxid[atm][0],0.001] for atm in atSeq]
#         if CheckAtms(Atypes):
#             oldPairs = data['RMC']['RMCProfile'].get('Pairs',{})
#             Pairs = {}
# #                for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
#             for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA)] for i in range(lenA)]:
#                 for pair in pairs:
#                     if pair in oldPairs:
#                         Pairs[pair] = oldPairs[pair]
#                     else:
#                         Pairs[pair] = [0.0,0.0,0.0]
#             data['RMC']['RMCProfile'].update({'aTypes':aTypes,'atSeq':atSeq,'Pairs':Pairs,'Oxid':atOxid,})

#         if not data['RMC']['RMCProfile'] or 'metadata' not in RMCPdict:
#             Pairs = {}
# #                for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
#             for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA)] for i in range(lenA)]:
#                 for pair in pairs:
#                     Pairs[pair] = [0.0,0.0,0.0]
#             BVSpairs = []
#             if lenA > 1:
# #                    for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
#                 for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA)] for i in range(lenA)]:
#                     BVSpairs += pair
#             BVS = {pairs:[0.0,0.0,0.0,0.0] for pairs in BVSpairs}
#             files = {'Neutron real space data; G(r): ':['Select',0.05,'G(r)','RMC',],
#                       'Neutron reciprocal space data; F(Q): ':['Select',0.05,'F(Q)','RMC',],
#                       'Neutron reciprocal space data; S(Q): ':['Select',0.05,'S(Q)','RMC',],
# #                          'Xray real space data; G(r): ':['Select',0.01,'G(r)','RMC',],
#                       'Xray reciprocal space data; F(Q): ':['Select',0.01,'F(Q)','RMC',],}
#             runTimes = [10.,1.]
#             metadata = {'title':'none','owner':'no one','date':str(time.ctime()),'temperature':'300K',
#                 'material':'nothing','phase':'vacuum','comment':'none ','source':'nowhere'}
#             data['RMC']['RMCProfile'].update({'SuperCell':[1,1,1],'UseSampBrd':[True,True],'aTypes':aTypes,
#                 'histogram':['',1.0],'files':files,'metadata':metadata,'FitScale':False,'atSeq':atSeq,
#                 'runTimes':runTimes,'ReStart':[False,False],'BVS':BVS,'Oxid':atOxid,'useBVS':False,'Swaps':[],
#                 'AveCN':[],'FxCN':[],'Potentials':{'Angles':[],'Angle search':10.,'Stretch':[],'Pairs':Pairs,
#                 'Stretch search':10.,'Pot. Temp.':300.,'useGPU':False,}})

# #            data['RMC']['RMCProfile']['aTypes'] = {aTypes[atype] for atype in aTypes if atype in Atypes}
#         data['RMC']['RMCProfile']['Isotope'] = copy.copy(data['General']['Isotope'])
#         data['RMC']['RMCProfile']['Isotopes'] = copy.deepcopy(data['General']['Isotopes'])
#         data['RMC']['RMCProfile']['NoAtoms'] = copy.copy(data['General']['NoAtoms'])
#         RMCPdict = data['RMC']['RMCProfile']
# #patches
#         if 'FitScale' not in RMCPdict:
#             RMCPdict['FitScale'] = False
#         if 'useGPU' not in RMCPdict:
#             RMCPdict['useGPU'] = False

# #end patches

#         def OnHisto(event):
#             RMCPdict['histogram'][0] = histo.GetStringSelection()

#         def OnSize(event):
#             RMCPdict['UseSampBrd'][0] = samSize.GetValue()

#         def OnStrain(event):
#             RMCPdict['UseSampBrd'][1] = strain.GetValue()

#         def OnFitScale(event):
#             RMCPdict['FitScale'] = not RMCPdict['FitScale']

#         def SetRestart(invalid,value,tc):
#             RMCPdict['ReStart'] = [True,True]

#         def OnUseBVS(event):
#             RMCPdict['useBVS'] = not RMCPdict['useBVS']
#             wx.CallAfter(UpdateRMC)

#         def OnAddSwap(event):
#             RMCPdict['Swaps'].append(['','',0.0,])
#             wx.CallAfter(UpdateRMC)

#         def OnAddFxCN(event):
#             RMCPdict['FxCN'].append(['','',0.5,2.0,6,1.0,0.00001])
#             wx.CallAfter(UpdateRMC)

#         def OnAddAveCN(event):
#             RMCPdict['AveCN'].append(['','',0.5,2.0,6.,0.00001])
#             wx.CallAfter(UpdateRMC)

#         def OnAddAnglePot(event):
#             RMCPdict['Potentials']['Angles'].append(['','','',0.,0.,0.,0.])
#             wx.CallAfter(UpdateRMC)

#         def OnAddBondPot(event):
#             RMCPdict['Potentials']['Stretch'].append(['','',0.,0.])
#             wx.CallAfter(UpdateRMC)

#         def GetTimeSizer():

#             def OnUseGPU(event):
#                 RMCPdict['useGPU'] = not RMCPdict['useGPU']

#             timeSizer = wx.BoxSizer(wx.HORIZONTAL)
#             timeSizer.Add(wx.StaticText(G2frame.FRMC,label=' Total running time (min): '),0,WACV)
#             timeSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['runTimes'],0,xmin=0.,size=(70,25)),0,WACV)
#             timeSizer.Add(wx.StaticText(G2frame.FRMC,label=' Save interval time (min): '),0,WACV)
#             timeSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['runTimes'],1,xmin=0.1,xmax=20.,size=(50,25)),0,WACV)
#             usegpu = wx.CheckBox(G2frame.FRMC,label=' use GPU?')
#             usegpu.SetValue(RMCPdict['useGPU'])
#             usegpu.Bind(wx.EVT_CHECKBOX,OnUseGPU)
#             timeSizer.Add(usegpu,0,WACV)
#             return timeSizer

#         # def GetSuperSizer(Xmax):
#         #     superSizer = wx.BoxSizer(wx.HORIZONTAL)
#         #     axes = ['X','Y','Z']
#         #     for i,ax in enumerate(axes):
#         #         superSizer.Add(wx.StaticText(G2frame.FRMC,label=' %s-axis: '%ax),0,WACV)
#         #         superSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['SuperCell'],
#         #             i,xmin=1,xmax=xamx,size=(50,25),OnLeave=SetRestart),0,WACV)
#         #     return superSizer

#         def GetBvsSizer(pnl):

#             def OnResetBVS(event):
#                 Obj = event.GetEventObject()
#                 pair = Indx[Obj.GetId()]
#                 pId = [key for key in RMCPdict['BVS']].index(pair)+1
#                 nId = len(RMCPdict['BVS'])+1
#                 dist = G2elem.GetBVS(pair,RMCPdict['atSeq'],RMCPdict['Oxid'])
#                 if dist:
#                     RMCPdict['BVS'][pair] = [dist,0.37,3.0]
#                     bvsCh = bvsSizer.GetChildren()
#                     addr = 2*nId+pId
#                     bvsCh[addr].Window.SetValue('%6.3f'%dist)
#                     bvsCh[addr+nId].Window.SetValue('0.37')
#                     bvsCh[addr+2*nId].Window.SetValue('3.00')

#             bvsSizer = wx.FlexGridSizer(len(RMCPdict['BVS'])+1,5,5)
#             bvsSizer.Add((5,5),0)
#             for pair in RMCPdict['BVS']:
#                 bvsSizer.Add(wx.StaticText(pnl,label=pair),0,WACV)
#             bvsSizer.Add(wx.StaticText(pnl,label=' Reset:'),0,WACV)
#             for pair in RMCPdict['BVS']:
#                 reset = wx.Button(pnl,label='Yes')
#                 bvsSizer.Add(reset,0,WACV)
#                 reset.Bind(wx.EVT_BUTTON,OnResetBVS)
#                 Indx[reset.GetId()] = pair
#             bvsSizer.Add(wx.StaticText(pnl,label=' Bond length:'),0,WACV)
#             for pair in RMCPdict['BVS']:
#                 bvsSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['BVS'][pair],0,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#             bvsSizer.Add(wx.StaticText(pnl,label=' B constant (0.37): '),0,WACV)
#             for pair in RMCPdict['BVS']:
#                 bvsSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['BVS'][pair],1,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#             bvsSizer.Add(wx.StaticText(pnl,label=' Cut off: '),0,WACV)
#             for pair in RMCPdict['BVS']:
#                 bvsSizer.Add(G2G.ValidatedTxtCtrl(pnl,RMCPdict['BVS'][pair],2,xmin=0.,xmax=10.,size=(50,25)),0,WACV)
#             return bvsSizer

#         def GetFxcnSizer():

#             def OnDelFxCN(event):
#                 Obj = event.GetEventObject()
#                 fxCN = Indx[Obj.GetId()]
#                 del RMCPdict['FxCN'][fxCN]
#                 wx.CallAfter(UpdateRMC)

#             def OnFxcnAtSel(event):
#                 Obj = event.GetEventObject()
#                 ifxCN,i = Indx[Obj.GetId()]
#                 RMCPdict['FxCN'][ifxCN][i] = Obj.GetStringSelection()

#             fxcnSizer = wx.FlexGridSizer(8,5,5)
#             atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#             fxcnLabels = [' ','Atom-1','Atom-2','min dist','max dist','CN','fraction','weight']
#             for lab in fxcnLabels:
#                 fxcnSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#             for ifx,fxCN in enumerate(RMCPdict['FxCN']):
#                 delBtn = wx.Button(G2frame.FRMC,label='Delete')
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelFxCN)
#                 Indx[delBtn.GetId()] = ifx
#                 fxcnSizer.Add(delBtn,0,WACV)
#                 for i in [0,1]:
#                     atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                     atmSel.SetStringSelection(fxCN[i])
#                     atmSel.Bind(wx.EVT_COMBOBOX,OnFxcnAtSel)
#                     Indx[atmSel.GetId()] = [ifx,i]
#                     fxcnSizer.Add(atmSel,0,WACV)
#                 fxcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,2,xmin=0.,xmax=5.,size=(50,25)),0,WACV)
#                 fxcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,3,xmin=0.,xmax=5.,size=(50,25)),0,WACV)
#                 fxcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,4,xmin=1,xmax=12,size=(50,25)),0,WACV)
#                 fxcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,5,xmin=0.,xmax=1.,size=(50,25)),0,WACV)
#                 fxcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,6,xmin=0.,size=(50,25)),0,WACV)
#             return fxcnSizer

#         def GetAvcnSizer():

#             def OnDelAvCN(event):
#                 Obj = event.GetEventObject()
#                 fxCN = Indx[Obj.GetId()]
#                 del RMCPdict['AveCN'][fxCN]
#                 wx.CallAfter(UpdateRMC)

#             def OnAvcnAtSel(event):
#                 Obj = event.GetEventObject()
#                 ifxCN,i = Indx[Obj.GetId()]
#                 RMCPdict['AveCN'][ifxCN][i] = Obj.GetStringSelection()

#             avcnSizer = wx.FlexGridSizer(7,5,5)
#             atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#             fxcnLabels = [' ','Atom-1','Atom-2','min dist','max dist','CN','weight']
#             for lab in fxcnLabels:
#                 avcnSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#             for ifx,fxCN in enumerate(RMCPdict['AveCN']):
#                 delBtn = wx.Button(G2frame.FRMC,label='Delete')
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelAvCN)
#                 Indx[delBtn.GetId()] = ifx
#                 avcnSizer.Add(delBtn,0,WACV)
#                 for i in [0,1]:
#                     atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                     atmSel.SetStringSelection(fxCN[i])
#                     atmSel.Bind(wx.EVT_COMBOBOX,OnAvcnAtSel)
#                     Indx[atmSel.GetId()] = [ifx,i]
#                     avcnSizer.Add(atmSel,0,WACV)
#                 avcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,2,xmin=0.,xmax=5.,size=(50,25)),0,WACV)
#                 avcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,3,xmin=0.,xmax=5.,size=(50,25)),0,WACV)
#                 avcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,4,xmin=1.,xmax=12.,size=(50,25)),0,WACV)
#                 avcnSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,fxCN,5,xmin=0.,size=(50,25)),0,WACV)
#             return avcnSizer

#         def GetAngleSizer():

#             def OnDelAngle(event):
#                 Obj = event.GetEventObject()
#                 angle = Indx[Obj.GetId()]
#                 del RMCPdict['Potentials']['Angles'][angle]
#                 wx.CallAfter(UpdateRMC)

#             def OnAngleAtSel(event):
#                 Obj = event.GetEventObject()
#                 angle,i = Indx[Obj.GetId()]
#                 RMCPdict['Potentials']['Angles'][angle][i] = Obj.GetStringSelection()

#             def SetRestart1(invalid,value,tc):
#                 RMCPdict['ReStart'][1] = True

#             atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#             angleSizer = wx.FlexGridSizer(8,5,5)
#             fxcnLabels = [' ','Atom-A','Atom-B','Atom-C',' ABC angle','AB dist','BC dist','potential']
#             for lab in fxcnLabels:
#                 angleSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#             for ifx,angle in enumerate(RMCPdict['Potentials']['Angles']):
#                 delBtn = wx.Button(G2frame.FRMC,label='Delete')
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelAngle)
#                 Indx[delBtn.GetId()] = ifx
#                 angleSizer.Add(delBtn,0,WACV)
#                 for i in [0,1,2]:
#                     atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                     atmSel.SetStringSelection(angle[i])
#                     atmSel.Bind(wx.EVT_COMBOBOX,OnAngleAtSel)
#                     Indx[atmSel.GetId()] = [ifx,i]
#                     angleSizer.Add(atmSel,0,WACV)
#                 angleSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,angle,3,xmin=0.,xmax=180.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#                 angleSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,angle,4,xmin=0.5,xmax=5.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#                 angleSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,angle,5,xmin=0.5,xmax=5.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#                 angleSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,angle,6,xmin=0.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#             return angleSizer

#         def GetBondSizer():

#             def OnDelBond(event):
#                 Obj = event.GetEventObject()
#                 bond = Indx[Obj.GetId()]
#                 del RMCPdict['Potentials']['Stretch'][bond]
#                 wx.CallAfter(UpdateRMC)

#             def OnBondAtSel(event):
#                 Obj = event.GetEventObject()
#                 bond,i = Indx[Obj.GetId()]
#                 RMCPdict['Potentials']['Stretch'][bond][i] = Obj.GetStringSelection()

#             def SetRestart1(invalid,value,tc):
#                 RMCPdict['ReStart'][1] = True

#             atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
#             bondSizer = wx.FlexGridSizer(5,5,5)
#             fxcnLabels = [' ','Atom-A','Atom-B',' AB dist','potential']
#             for lab in fxcnLabels:
#                 bondSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
#             for ifx,bond in enumerate(RMCPdict['Potentials']['Stretch']):
#                 delBtn = wx.Button(G2frame.FRMC,label='Delete')
#                 delBtn.Bind(wx.EVT_BUTTON,OnDelBond)
#                 Indx[delBtn.GetId()] = ifx
#                 bondSizer.Add(delBtn,0,WACV)
#                 for i in [0,1]:
#                     atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#                     atmSel.SetStringSelection(bond[i])
#                     atmSel.Bind(wx.EVT_COMBOBOX,OnBondAtSel)
#                     Indx[atmSel.GetId()] = [ifx,i]
#                     bondSizer.Add(atmSel,0,WACV)
#                 bondSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,bond,2,xmin=0.,xmax=5.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
#                 bondSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,bond,3,xmin=0.,size=(50,25)),0,WACV)
#             return bondSizer

#         Indx = {}

#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' Enter metadata items:'),0)
#         mainSizer.Add(GetMetaSizer(RMCPdict,['title','owner','material','phase','comment','source',]),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(GetTimeSizer(),0)

#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' Lattice multipliers; if changed will force reset of atom positions:'),0)
#         mainSizer.Add(GetSuperSizer(RMCPdict,20),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)

#         mSizer = wx.BoxSizer(wx.VERTICAL)
#         mSizer.Add(wx.StaticText(G2frame.FRMC,label='Enter atom settings'),0)
#         mSizer.Add(GetAtmChoice(G2frame.FRMC,RMCPdict),0)
#         mSizer.Add(wx.StaticText(G2frame.FRMC,label=' N.B.: be sure to set cations first && anions last in atom ordering'))
#         mainSizer.Add(mSizer)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         swapBox = wx.BoxSizer(wx.HORIZONTAL)
#         swapAdd = wx.Button(G2frame.FRMC,label='Add')
#         swapAdd.Bind(wx.EVT_BUTTON,OnAddSwap)
#         swapBox.Add(swapAdd,0,WACV)
#         swapBox.Add(wx.StaticText(G2frame.FRMC,label=' Atom swap probabilities: '),0,WACV)
#         mainSizer.Add(swapBox,0)
#         if len(RMCPdict['Swaps']):
#             mainSizer.Add(GetSwapSizer(RMCPdict),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)

#         mSizer = wx.BoxSizer(wx.VERTICAL)
#         mSizer.Add(wx.StaticText(G2frame.FRMC,label='Enter constraints && restraints via minimum && maximum distances for atom pairs:'),0)
#         mSizer.Add(GetPairSizer(G2frame.FRMC,RMCPdict),0)
#         mainSizer.Add(mSizer)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         useBVS = wx.CheckBox(G2frame.FRMC,label=' Use bond valence sum restraints for (set to 0 for non-bonded ones):')
#         useBVS.SetValue(RMCPdict.get('useBVS',False))
#         useBVS.Bind(wx.EVT_CHECKBOX,OnUseBVS)
#         mainSizer.Add(useBVS,0)
#         if RMCPdict.get('useBVS',False):
#             mSizer = wx.BoxSizer(wx.VERTICAL)
#             mSizer.Add(GetBvsSizer(G2frame.FRMC),0)
#             mainSizer.Add(mSizer)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         fxcnBox = wx.BoxSizer(wx.HORIZONTAL)
#         fxcnAdd = wx.Button(G2frame.FRMC,label='Add')
#         fxcnAdd.Bind(wx.EVT_BUTTON,OnAddFxCN)
#         fxcnBox.Add(fxcnAdd,0,WACV)
#         fxcnBox.Add(wx.StaticText(G2frame.FRMC,label=' Fixed coordination number restraint: '),0,WACV)
#         mainSizer.Add(fxcnBox,0)
#         if len(RMCPdict['FxCN']):
#             mainSizer.Add(GetFxcnSizer(),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         avcnBox = wx.BoxSizer(wx.HORIZONTAL)
#         avcnAdd = wx.Button(G2frame.FRMC,label='Add')
#         avcnAdd.Bind(wx.EVT_BUTTON,OnAddAveCN)
#         avcnBox.Add(avcnAdd,0,WACV)
#         avcnBox.Add(wx.StaticText(G2frame.FRMC,label=' Average coordination number restraint: '),0,WACV)
#         mainSizer.Add(avcnBox,0)
#         if len(RMCPdict['AveCN']):
#             mainSizer.Add(GetAvcnSizer(),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         pottempBox = wx.BoxSizer(wx.HORIZONTAL)
#         pottempBox.Add(wx.StaticText(G2frame.FRMC,label=' Potential temperature (K): '),0,WACV)
#         pottempBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['Potentials'],'Pot. Temp.',xmin=0.,xmax=1000.,size=(50,25)),0,WACV)
#         mainSizer.Add(pottempBox,0)
#         bondpotBox = wx.BoxSizer(wx.HORIZONTAL)
#         bondpotAdd = wx.Button(G2frame.FRMC,label='Add')
#         bondpotAdd.Bind(wx.EVT_BUTTON,OnAddBondPot)
#         bondpotBox.Add(bondpotAdd,0,WACV)
#         bondpotBox.Add(wx.StaticText(G2frame.FRMC,label=' A-B stretch potential restraints, search range (%): '),0,WACV)
#         bondpotBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['Potentials'],'Stretch search',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
#         mainSizer.Add(bondpotBox,0)
#         if len(RMCPdict['Potentials']['Stretch']):
#             mainSizer.Add(GetBondSizer(),0)

#         angpotBox = wx.BoxSizer(wx.HORIZONTAL)
#         angpotAdd = wx.Button(G2frame.FRMC,label='Add')
#         angpotAdd.Bind(wx.EVT_BUTTON,OnAddAnglePot)
#         angpotBox.Add(angpotAdd,0,WACV)
#         angpotBox.Add(wx.StaticText(G2frame.FRMC,label=' A-B-C angle potential restraints, search range (%): '),0,WACV)
#         angpotBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['Potentials'],'Angle search',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
#         mainSizer.Add(angpotBox,0)
#         if len(RMCPdict['Potentials']['Angles']):
#             mainSizer.Add(GetAngleSizer(),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' Select data:'),0)
#         histograms = data['Histograms']
#         histNames = list(histograms.keys())
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' Select one histogram for Bragg processing:'),0)
#         histoSizer = wx.BoxSizer(wx.HORIZONTAL)
#         histo = wx.ComboBox(G2frame.FRMC,choices=histNames,style=wx.CB_DROPDOWN|wx.TE_READONLY)
#         if RMCPdict['histogram'][0] == '':
#             RMCPdict['histogram'][0] = histo.GetStringSelection()
#         histo.SetStringSelection(RMCPdict['histogram'][0])
#         histo.Bind(wx.EVT_COMBOBOX,OnHisto)
#         histoSizer.Add(histo,0,WACV)
#         histoSizer.Add(wx.StaticText(G2frame.FRMC,label=' Weight '),0,WACV)
#         histoSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['histogram'],1,xmin=0.,xmax=10000.,size=(50,25)),0,WACV)
#         mainSizer.Add(histoSizer,0)

#         samSizer = wx.BoxSizer(wx.HORIZONTAL)
#         samSize = wx.CheckBox(G2frame.FRMC,label=' Use size broadening?')
#         samSize.SetValue(RMCPdict['UseSampBrd'][0])
#         samSize.Bind(wx.EVT_CHECKBOX,OnSize)
#         strain = wx.CheckBox(G2frame.FRMC,label=' Use mustrain broadening?')
#         strain.SetValue(RMCPdict['UseSampBrd'][1])
#         strain.Bind(wx.EVT_CHECKBOX,OnStrain)
#         samSizer.Add(samSize,0,WACV)
#         samSizer.Add(strain,0,WACV)
#         fitscale = wx.CheckBox(G2frame.FRMC,label=' Fit scale factors?')
#         fitscale.SetValue(RMCPdict['FitScale'])
#         fitscale.Bind(wx.EVT_CHECKBOX,OnFitScale)
#         samSizer.Add(fitscale,0,WACV)
#         mainSizer.Add(samSizer,0)

#         mainSizer.Add(FileSizer(RMCPdict))
#         return mainSizer

#     def PDFfitSizer(data):

#         mainSizer = wx.BoxSizer(wx.VERTICAL)
#         Indx = {}
#         def PDFParmSizer():

#             def OnShape(event):
#                 RMCPdict['shape'] = shape.GetValue()
#                 wx.CallAfter(UpdateRMC)

#             parmSizer = wx.FlexGridSizer(3,6,5,5)
#             Names = ['delta1','delta2','sratio','rcut','spdiameter']
#             Names2 = ['stepcut',]
#             for name in Names:

#                 def OnRefine(event):
#                     Obj = event.GetEventObject()
#                     name = Indx[Obj.GetId()]
#                     RMCPdict[name][1] = not RMCPdict[name][1]

#                 if name == 'spdiameter' and RMCPdict.get('shape','sphere') != 'sphere':
#                     pass
#                 else:
#                     parmSizer.Add(wx.StaticText(G2frame.FRMC,label=name),0,WACV)
#                     if name == 'rcut':
#                         parmSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,name,xmin=0.,size=(70,25)),0,WACV)
#                         parmSizer.Add((5,5))
#                         continue
#                     parmSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict[name],0,xmin=0.,size=(70,25)),0,WACV)
#                     refine = wx.CheckBox(G2frame.FRMC,label='Refine')
#                     refine.SetValue(RMCPdict[name][1])
#                     refine.Bind(wx.EVT_CHECKBOX,OnRefine)
#                     Indx[refine.GetId()] = name
#                     parmSizer.Add(refine,0,WACV)
#             parmSizer.Add(wx.StaticText(G2frame.FRMC,label=' Shape'),0,WACV)
#             shape = wx.ComboBox(G2frame.FRMC,choices=['sphere','stepcut'],style=wx.CB_DROPDOWN|wx.TE_READONLY)
#             shape.SetStringSelection(RMCPdict.get('shape','sphere'))
#             shape.Bind(wx.EVT_COMBOBOX,OnShape)
#             parmSizer.Add(shape,0,WACV)
#             if RMCPdict.get('shape','sphere') == 'stepcut':
#                 for name in Names2:
#                     parmSizer.Add(wx.StaticText(G2frame.FRMC,label=name),0,WACV)
#                     parmSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,name,xmin=0.,size=(70,25)),0,WACV)

#             return parmSizer

#         def OnSpaceGroup(event):
#             # try a lookup on the user-supplied name
#             SpcGp = GetSpGrpfromUser(G2frame.FRMC,SpGrp)
#             SGErr,SGData = G2spc.SpcGroup(SpcGp)
#             if SGErr:
#                 text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
#                 SGTxt.SetLabel(RMCPdict.get('Target SpGrp','P 1'))
#                 msg = 'Space Group Error'
#                 Text = '\n'.join(text)
#                 wx.MessageBox(Text,caption=msg,style=wx.ICON_EXCLAMATION)
#             else:
#                 text,table = G2spc.SGPrint(SGData)
#                 RMCPdict['SGData'] = SGData
#                 SGTxt.SetLabel(RMCPdict['SGData']['SpGrp'])
#                 msg = 'Target Space Group Information'
#                 G2G.SGMessageBox(G2frame.FRMC,msg,text,table).Show()
#             G2spc.UpdateSytSym(data)

#         def OnCellRef(event):
#             RMCPdict['cellref'] = not RMCPdict['cellref']

#         def AtomSizer():

#             def OnSetVal(event):
#                 r,c = event.GetRow(),event.GetCol()
#                 if c > 0:
#                     strval = atmGrid.GetCellValue(r,c).strip()
#                     try:
# #                            if strval == '' or ('@' in strval and int(strval.split('@')[-1]) >= 20):
#                         if strval == '' or '@' in strval:
#                             RMCPdict['AtomConstr'][r][c+1] = strval
#                         else:
#                             raise ValueError
#                     except ValueError:
#                         atmGrid.SetCellValue(r,c,RMCPdict['AtomConstr'][r][c+1])
#                         wx.MessageBox('ERROR - atom constraints must be blank or have "@n" with n >= 20',
#                             style=wx.ICON_ERROR)
#                     wx.CallAfter(UpdateRMC)

#             def OnUisoRefine(event):
#                 RMCPdict['UisoRefine'] = uiso.GetValue()
#                 nextP = 80
#                 oType = ''
#                 oName = ''
#                 for atom in RMCPdict['AtomConstr']:
#                     if RMCPdict['UisoRefine'] == 'No':
#                         atom[6] = ''
#                     elif 'same' in RMCPdict['UisoRefine']:
#                         atom[6] = '@81'
#                         RMCPdict['AtomVar']['@81'] = 0.005
#                     elif 'type' in RMCPdict['UisoRefine']:
#                         if atom[1] != oType:
#                             oType = atom[1]
#                             nextP += 1
#                         atom[6] = '@%d'%nextP
#                         RMCPdict['AtomVar']['@%d'%nextP] = 0.005
#                     elif 'parent' in RMCPdict['UisoRefine']:
#                         if atom[0] != oName:
#                             oName = atom[0]
#                             nextP += 1
#                         atom[6] = '@%d'%nextP
#                         RMCPdict['AtomVar']['@%d'%nextP] = 0.005
#                 wx.CallAfter(UpdateRMC)

#             atmSizer = wx.BoxSizer(wx.VERTICAL)
#             atmSizer.Add(wx.StaticText(G2frame.FRMC,label=' Atom Constraints; enter as e.g. "@n" or "0.5-@n"; n>=20 && "@n" should be at end'))
#             uisoSizer = wx.BoxSizer(wx.HORIZONTAL)
#             uisoSizer.Add(wx.StaticText(G2frame.FRMC,label=' Refine Uiso? '),0,WACV)
#             uiso = wx.ComboBox(G2frame.FRMC,choices=['No','by type','by parent name','all same'],style=wx.CB_DROPDOWN|wx.TE_READONLY)
#             uiso.SetValue(RMCPdict['UisoRefine'])
#             uiso.Bind(wx.EVT_COMBOBOX,OnUisoRefine)
#             uisoSizer.Add(uiso,0,WACV)
#             atmSizer.Add(uisoSizer)

#             table = [item[1:] for item in RMCPdict['AtomConstr']]
#             addCol = False
#             if len(RMCPdict['AtomConstr'][0]) > 6:
#                 addCol = True
#             colLabels = ['Type','x constraint','y constraint','z  constraint','frac constr','Uiso constr']
#             rowLabels = [item[0] for item in RMCPdict['AtomConstr']]
#             Types = 6*[wg.GRID_VALUE_STRING,]
#             if addCol:
#                 colLabels += ['sym opr',]
#                 Types = 7*[wg.GRID_VALUE_STRING,]
#             atmTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
#             atmGrid = G2G.GSGrid(G2frame.FRMC)
#             atmGrid.SetTable(atmTable, True,useFracEdit=False)
#             atmGrid.AutoSizeColumns(True)
#             atmGrid.Bind(wg.EVT_GRID_CELL_CHANGED, OnSetVal)
#             atmSizer.Add(atmGrid)
#             return atmSizer

#         def AtomVarSizer():
#             atomVarSizer = wx.FlexGridSizer(0,8,5,5)
#             for item in RMCPdict['AtomVar']:
#                 atomVarSizer.Add(wx.StaticText(G2frame.FRMC,label=item),0,WACV)
#                 atomVarSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['AtomVar'],
#                     item,xmin=-3.,xmax=3.,size=(70,25)),0,WACV)
#             return atomVarSizer

#         txt = wx.StaticText(G2frame.FRMC,label=
#                                 'For use of PDFfit, please cite: '+
#                                 G2G.GetCite('PDFfit2'))
#         txt.Wrap(500)
#         mainSizer.Add(txt)
#         mainSizer.Add((5,5))
#         if 'PDFfit' not in data['RMC'] or not data['RMC']['PDFfit'] or 'delta1' not in data['RMC']['PDFfit']:
#             if 'PDFfit' not in data['RMC']:
#                 data['RMC']['PDFfit'] = {}
#             metadata = {'title':'none','date':str(time.ctime()),'temperature':'300K','doping':0}
#             files = {'Neutron real space data; G(r): ':['Select',0.05,'G(r)','RMC',],
#                       'Xray real space data; G(r): ':['Select',0.01,'G(r)','RMC',],}
#             data['RMC']['PDFfit'].update({'files':files,'ReStart':[False,False],'metadata':metadata,
#             'delta1':[0.,False],'delta2':[0.,False],'spdiameter':[0.,False],'refinement':'normal',
#             'sratio':[1.,False],'rcut':0.0,'stepcut':0.0,'shape':'sphere','cellref':False,
#             'SeqDataType':'X','SeqCopy':True,'SeqReverse':False,
#             'Xdata':{'dscale':[1.0,False],'Datarange':[0.,30.],'Fitrange':[0.,30.],'qdamp':[0.03,False],'qbroad':[0.,False]},
#             'Ndata':{'dscale':[1.0,False],'Datarange':[0.,30.],'Fitrange':[0.,30.],'qdamp':[0.03,False],'qbroad':[0.,False]},})

#         RMCPdict = data['RMC']['PDFfit']
# #patch
#         if 'AtomConstr' not in RMCPdict:        #keep this one
#             RMCPdict['AtomConstr'] = []
#         if 'AtomVar' not in RMCPdict:
#             RMCPdict['AtomVar'] = {}
#         if 'SGData' not in RMCPdict:
#             RMCPdict['SGData'] = G2spc.SpcGroup('P 1')[1]
#         if 'refinement' not in RMCPdict:
#             RMCPdict['refinement'] = 'normal'
#         if 'cellref' not in RMCPdict:
#             RMCPdict['cellref'] = False
#         if 'metadata' not in RMCPdict:
#             RMCPdict['metadata'] = {'title':'none','date':str(time.ctime()),'temperature':'300K','doping':0}
#         if 'SeqDataType' not in RMCPdict:
#             RMCPdict['SeqDataType'] = 'X'
#         if 'SeqCopy' not in RMCPdict:
#             RMCPdict['SeqCopy'] = False
#             RMCPdict['SeqReverse'] = False
#         if 'UisoRefine' not in RMCPdict:
#             RMCPdict['UisoRefine'] = 'No'
# #end patch
#         Atoms = data['Atoms']
#         cx,ct,cs,ci = G2mth.getAtomPtrs(data)
#         if not RMCPdict['AtomConstr']:
#             for atom in Atoms:
#                 RMCPdict['AtomConstr'].append([atom[ct-1],atom[ct],'','','','',''])
#         else:       #update name/type changes
#             for iatm,atom in enumerate(Atoms):
#                 RMCPdict['AtomConstr'][iatm][:2] = atom[ct-1:ct+1]

#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' Enter metadata items:'),0)
#         mainSizer.Add(GetMetaSizer(RMCPdict,['title','date','temperature','doping']),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         SgSizer = wx.BoxSizer(wx.HORIZONTAL)
#         SgSizer.Add(wx.StaticText(G2frame.FRMC,label=' Target space group: '),0,WACV)

#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label='PDFfit phase structure parameters:'))

#         SpGrp = RMCPdict['SGData']['SpGrp']
#         SGTxt = wx.Button(G2frame.FRMC,wx.ID_ANY,SpGrp,size=(100,-1))
#         SGTxt.Bind(wx.EVT_BUTTON,OnSpaceGroup)
#         SgSizer.Add(SGTxt,0,WACV)
#         mainSizer.Add(SgSizer)

#         cellref = wx.CheckBox(G2frame.FRMC,label=' Refine unit cell?')
#         cellref.SetValue(RMCPdict['cellref'])
#         cellref.Bind(wx.EVT_CHECKBOX,OnCellRef)
#         mainSizer.Add(cellref)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label='PDFfit atom parameters:'))
#         mainSizer.Add(AtomSizer())

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label='PDFfit starting atom variables:'))
#         G2pwd.GetPDFfitAtomVar(data,RMCPdict)
#         mainSizer.Add(AtomVarSizer())

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' PDFfit phase profile coefficients:'))
#         mainSizer.Add(PDFParmSizer(),0)

#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         mainSizer.Add(FileSizer(RMCPdict))
#         return mainSizer

# ####start of UpdateRMC
#     G2frame.GetStatusBar().SetStatusText('',1)
#     G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC,False)
#     G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC,False)
#     if G2frame.RMCchoice == 'RMCProfile':
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC,True)
#         #G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_STOPRMC,False)
#     elif G2frame.RMCchoice == 'fullrmc':
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC,False)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC,False)
#         #G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_STOPRMC,False)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC,True)
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC,True)
#     try:
#         if G2frame.FRMC.GetSizer():
#             G2frame.FRMC.GetSizer().Clear(True)
#     except: #wxAssertionError from C++
#         pass
#     mainSizer = wx.BoxSizer(wx.VERTICAL)
#     if not len(data['Atoms']):
#         mainSizer.Add(wx.StaticText(G2frame.FRMC,label='No atoms found - PDF fitting not possible'))
#     else:
#         runFile = ' '
#         choice = ['RMCProfile','fullrmc','PDFfit']
#         topSizer = wx.BoxSizer(wx.HORIZONTAL)
#         RMCsel = wx.RadioBox(G2frame.FRMC,-1,' Select RMC method:',choices=choice)
#         RMCsel.SetStringSelection(G2frame.RMCchoice)
#         RMCsel.Bind(wx.EVT_RADIOBOX, OnRMCselect)
#         topSizer.Add(RMCsel,0)
#         topSizer.Add((20,0))
#         txt = wx.StaticText(G2frame.FRMC,
#             label='NB: if you change any of the entries below, you must redo the Operations/Setup RMC step above to apply them before doing Operations/Execute')
#         txt.Wrap(250)
#         topSizer.Add(txt,0)
#         mainSizer.Add(topSizer,0)
#         RMCmisc['RMCnote'] = wx.StaticText(G2frame.FRMC)
#         mainSizer.Add(RMCmisc['RMCnote'])
#         G2G.HorizontalLine(mainSizer,G2frame.FRMC)
#         if G2frame.RMCchoice == 'fullrmc':
#             RMCPdict = data['RMC']['fullrmc']
#             mainSizer.Add(fullrmcSizer(RMCPdict))

#         elif G2frame.RMCchoice ==  'RMCProfile':
#             RMCPdict = data['RMC']['RMCProfile']
#             mainSizer.Add(RMCProfileSizer(RMCPdict))

#         else:       #PDFfit
#             mainSizer.Add(PDFfitSizer(data))

#     topSizer = G2frame.dataWindow.topBox
#     topSizer.Clear(True)
#     parent = G2frame.dataWindow.topPanel
#     lbl= f"PDF fitting options {data['General']['Name']!r}"[:60]
#     topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
#     topSizer.Add((-1,-1),1,wx.EXPAND)
#     topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
#     wx.CallAfter(G2frame.dataWindow.SetDataSize)
#     SetPhaseWindow(G2frame.FRMC,mainSizer)

#     if G2frame.RMCchoice == 'PDFfit' and not checkPDFfit(G2frame):
#         RMCmisc['RMCnote'].SetLabel('PDFfit may not be installed or operational')
#     elif G2frame.RMCchoice == 'fullrmc' and G2pwd.findfullrmc() is None:
#         msg = ('The fullrmc Python image is not found.'+
#                 ' Do you want it installed for you from '+
#     'https://github.com/bachiraoun/fullrmc/tree/master/standalones?'+
#                 '\n\n40-50 Mb (download times vary)')
#         dlg = wx.MessageDialog(G2frame,msg,'Install fullrmc',wx.YES|wx.NO)
#         try:
#             dlg.CenterOnParent()
#             result = dlg.ShowModal()
#         finally:
#             dlg.Destroy()
#         if result == wx.ID_YES:
#             wx.BeginBusyCursor()
#             G2pwd.fullrmcDownload()
#             wx.EndBusyCursor()
#         else:
#             RMCmisc['RMCnote'].SetLabel('Note that fullrmc is not installed or was not located')

# def OnSetupRMC(event):
#     written = lambda fil: print(f' {fil} written')
#     generalData = data['General']
#     if not G2frame.GSASprojectfile:     #force a project save
#         G2frame.OnFileSaveas(event)
#     dName = G2frame.LastGPXdir
#     os.chdir(dName)
#     print(f'Writing input files in directory {dName!r}')
#     if G2frame.RMCchoice == 'fullrmc':
#         RMCPdict = data['RMC']['fullrmc']
#         pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
#         pName = pName.replace(' ','_')
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
#         if RMCPdict['Swaps']:
#             wx.MessageDialog(G2frame, G2G.StripIndents(
#                     '''GSAS-II does not yet fully support use of swapping in fullrmc.
#                     Edit the script by hand before using.''',True),
#                     'No swaps yet',wx.OK).ShowModal()
#         #--------- debug stuff
#         # if GSASIIpath.GetConfigValue('debug'):
#         #     print('reloading',G2pwd)
#         #     import imp
#         #     imp.reload(G2pwd)
#         #--------- end debug stuff
#         rname = G2pwd.MakefullrmcRun(pName,data,RMCPdict)
#         print('build of fullrmc file {} completed'.format(rname))
#     elif G2frame.RMCchoice == 'RMCProfile':
#         if ' ' in generalData['Name']:
#             wx.MessageDialog(G2frame,'ERROR: Phase name has space; change phase name','Bad phase name',wx.ICON_ERROR).ShowModal()
#             G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#             return
#         dName = G2frame.LastGPXdir
#         pName = generalData['Name']
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
#         RMCPdict = data['RMC']['RMCProfile']
#         PWId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,RMCPdict['histogram'][0])
#         if PWId:
#             PWDdata = G2frame.GetPWDRdatafromTree(PWId)
#             histoName = G2frame.GPXtree.GetItemPyData(PWId)[2]
#             Size = data['Histograms'][histoName]['Size']
#             Mustrain = data['Histograms'][histoName]['Mustrain']
#             reset = False
#             written(G2pwd.MakeInst(PWDdata,pName,Size,Mustrain,RMCPdict['UseSampBrd']))
#             backfile = G2pwd.MakeBack(PWDdata,pName)
#             if backfile is None:
#                 print(' Chebyschev-1 background not used; no .back file written')
#                 wx.MessageDialog(G2frame,' Chebyschev-1 background not used; '+ \
#                     'no .back file written & RMCProfile will not run','Wrong background function',wx.OK).ShowModal()
#                 G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#                 return
#             else:
#                 written(backfile)
#             written(G2pwd.MakeBragg(PWDdata,pName,data))
#             if RMCPdict['ReStart'][0]:
#                 if os.path.isfile(pName+'.his6f'):
#                     os.remove(pName+'.his6f')
#                 RMC6f,reset = G2pwd.MakeRMC6f(PWDdata,pName,data,RMCPdict)
#                 written(RMC6f)
#             fname = G2pwd.MakeRMCPdat(PWDdata,pName,data,RMCPdict)
#             if 'Error' in fname:
#                 print(fname)
#                 wx.MessageDialog(G2frame,fname,'Missing reflection list',wx.OK).ShowModal()
#                 G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#                 return
#             written(fname)
#             print('RMCProfile file build completed')
#             RMCPdict['ReStart'] = [False,False]
#             if reset:
#                 wx.MessageDialog(G2frame,' Vacancies found & "Va" atoms added to list. '+ \
#                     'You may need to revise RMCProfile setup parameters.','Repeat Setup RMC',wx.OK).ShowModal()
#                 wx.CallAfter(UpdateRMC)
#         else:
#             print('RMCProfile file build failed - no histogram selected')
#             G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#     elif G2frame.RMCchoice == 'PDFfit':
#         if ' ' in generalData['Name']:
#             wx.MessageDialog(G2frame,'ERROR: Phase name has space; change phase name','Bad phase name',wx.ICON_ERROR).ShowModal()
#             G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
#             return
#         G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
#         RMCPdict = data['RMC']['PDFfit']
#         msg = G2pwd.MakePDFfitAtomsFile(data,RMCPdict)
#         if msg:
#             G2G.G2MessageBox(G2frame,'ERROR: '+msg,'PDFfit setup failure')
#             return
#         fname = G2pwd.MakePDFfitRunFile(data,RMCPdict)
#         if fname is None:
#             wx.MessageDialog(G2frame,'ERROR: failure to setup PDFfit; check console','PDFfit setup failure',wx.ICON_ERROR).ShowModal()
#         else:
#             written(fname)
#             print('PDFfit file build completed')

# def RunPDFfit(event):
#     generalData = data['General']
#     ISOdict = data['ISODISTORT']
#     PDFfit_exec = G2pwd.findPDFfit()  #returns location of python with PDFfit installed
#     if not PDFfit_exec:
#         wx.MessageBox(''' PDFfit2 is not currently installed for this platform.
# Please contact us for assistance''',caption='No PDFfit2',style=wx.ICON_INFORMATION)
#         return
#     RMCPdict = data['RMC']['PDFfit']
#     pName = generalData['Name'].replace(' ','_')
#     if 'sequential' in RMCPdict['refinement']:
#         rname = 'Seq_PDFfit.py'
#     else:
#         rname = pName+'-PDFfit.py'
#         if not os.path.exists(rname):
#             wx.MessageBox(f'File {rname} does not exist. Has the Operations/"Setup RMC" menu command been run?',
#                               caption='Run setup',style=wx.ICON_WARNING)
#             return
#     wx.MessageBox(' For use of PDFfit2, please cite:\n\n'+
#                       G2G.GetCite('PDFfit2'),
#                       caption='PDFfit2',style=wx.ICON_INFORMATION)
#     G2frame.OnFileSave(event)
#     print (' GSAS-II project saved')
#     if sys.platform.lower().startswith('win'):
#         batch = open('pdffit2.bat','w')
#         # Include an activate command here
#         p = os.path.split(PDFfit_exec)[0]
#         while p:
#             if os.path.exists(os.path.join(p,'Scripts','activate')):
#                 batch.write('call '+os.path.join(p,'Scripts','activate')+'\n')
#                 break
#             prevp = p
#             p = os.path.split(p)[0]
#             if prevp == p:
#                 print('Note, no activate command found')
#                 break
#         batch.write(PDFfit_exec+' '+rname+'\n')
#         # batch.write('pause')
#         if 'normal' in RMCPdict['refinement']:
#             batch.write('pause')
#         batch.close()
#     else:
#         batch = open('pdffit2.sh','w')
#         batch.write('#!/bin/bash\n')
#         # include an activate command here
#         p = os.path.split(PDFfit_exec)[0]
#         while p:
#             if os.path.exists(os.path.join(p,'bin','activate')):
#                 batch.write('source '+os.path.join(p,'Scripts','activate')+'\n')
#                 break
#             prevp = p
#             p = os.path.split(p)[0]
#             if prevp == p:
#                 print('Note, no activate command found')
#                 break

#         batch.write('cd ' + os.path.split(os.path.abspath(rname))[0] + '\n')
#         batch.write(PDFfit_exec + ' ' + os.path.abspath(rname) + '\n')
#         batch.close()
#     if 'sequential' in RMCPdict['refinement']:
#         Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential PDFfit2 results')
# #            if Id:
# #                saveSeqResult = G2frame.GPXtree.GetItemPyData(Id)
# #            else:
#         if not Id:
#             SeqResult = {}
#             Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential PDFfit2 results')
#         G2Names = [item.name for item in ISOdict['G2ModeList']]
#         SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
#         SeqResult['histNames'] = []         #this clears the previous seq. result!
#         SeqNames = []
#         for itm in range(len(RMCPdict['seqfiles'])):
#             SeqNames.append([itm,RMCPdict['seqfiles'][itm][0]])
#         if RMCPdict['SeqReverse']:
#             SeqNames.reverse()
#         nPDF = len(SeqNames)
#         pgbar = wx.ProgressDialog('Sequential PDFfit','PDF G(R) done = 0',nPDF+1,
#             style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
#         newParms = {}
#         for itm,Item in enumerate(SeqNames):
#             PDFfile = RMCPdict['seqfiles'][Item[0]]
#             pfdata = PDFfile[1]['G(R)'][1].T
# #                    pfname = PDFfile[0].replace(' ','_')
#             pfname = 'Seq_PDF.gr'
#             pfile = open(pfname,'w')
#             for dp in pfdata:
#                 pfile.write('%12.5f%12.5f\n'%(dp[0],dp[1]))
#             pfile.close()
#             rfile = open('Seq_PDFfit_template.py','r')
#             lines = rfile.readlines()       #template lines
#             rfile.close()
#             newlines = []
#             parms = {}
#             Np = 0
#             for line in lines:
#                 if '#sequential' in line:
#                     newlines += "pf.read_data('%s', '%s', 30.0, %.4f)\n"%(pfname,PDFfile[1]['Type'][0],PDFfile[1]['qdamp'][0])
#                     newlines += 'pf.setdata(1)\n'
#                     newlines += 'pf.pdfrange(1, %6.2f, %6.2f)\n'%(PDFfile[1]['Fitrange'][0],PDFfile[1]['Fitrange'][1])
#                     for item in ['dscale','qdamp','qbroad']:
#                         if PDFfile[1][item][1]:
#                             Np += 1
#                             newlines += 'pf.constrain(pf.%s(),"@%d")\n'%(item,Np)
#                             parms[item] = '%d'%Np
#                             if itm and RMCPdict['SeqCopy']:
#                                 newParms[parms[item]] = RMCPdict['Parms'][parms[item]]
#                             else:
#                                 if not itm and 'result' not in PDFfile[1]:
#                                     newParms[parms[item]] = PDFfile[1][item][0]
#                                 else:
#                                     newParms[parms[item]] = PDFfile[1]['result'][parms[item]][0]
#                 elif '#parameters' in line:
#                     startParms = RMCPdict['Parms']
#                     if newParms or RMCPdict['SeqCopy']:
#                         if newParms:
#                             startParms = newParms
#                         for iprm in startParms:
#                             if int(iprm) > 9:
#                                 break
#                             newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,startParms[iprm])
#                         print('Begin dscale: %d %.4f'%(itm,startParms['1']))
#                         for iprm in RMCPdict['Parms']:
#                             if isinstance(RMCPdict['Parms'][iprm],float):
#                                 newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,RMCPdict['Parms'][iprm])
#                             else:
#                                 newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,RMCPdict['Parms'][iprm][0])
#                     elif not RMCPdict['SeqCopy']:
#                         startParms = PDFfile[1]['result']
#                         for iprm in startParms:
#                             newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,startParms[iprm][0])
#                         print('Begin dscale: %d %.4f'%(itm,startParms['1']))
#                 else:
#                     newlines += line
#             rfile= open('Seq_PDFfit.py','w')
#             rfile.writelines(newlines)
#             rfile.close()
#             fName = 'Sequential_PDFfit'     #clean out old PDFfit output files
#             if os.path.isfile(fName+'.res'):
#                 os.remove(fName+'.res')
#             if os.path.isfile(fName+'.rstr'):
#                 os.remove(fName+'.rstr')
#             if os.path.isfile(fName+'.fgr'):
#                 os.remove(fName+'.fgr')

#             if sys.platform.lower().startswith('win'):
#                 Proc = subp.Popen('pdffit2.bat',creationflags=subp.CREATE_NEW_CONSOLE)
#                 Proc.wait()     #for it to finish before continuing on
#             else:
#                 if sys.platform == "darwin":
#                     GSASIIpath.MacRunScript(os.path.abspath('pdffit2.sh'))
#                 else:
#                     Proc = subp.Popen(['/bin/bash','pdffit2.sh'])
#                     Proc.wait()

#             newParms,Rwp =  G2pwd.UpdatePDFfit(data,RMCPdict)
#             if isinstance(newParms,str):
#                 wx.MessageBox('Singular matrix in PDFfit',caption='PDFfit2 failed',style=wx.ICON_INFORMATION)
#                 break
#             for item in ['dscale','qdamp','qbroad']:
#                 if PDFfile[1][item][1]:
#                     PDFfile[1][item][0] = newParms[parms[item]][0]
#             PDFfile[1]['result'] = copy.deepcopy(newParms)
#             parmDict = copy.deepcopy(newParms)
#             parmDict.update({'Temperature':PDFfile[1]['Temp']})
#             tempList = ['%s-%s'%(parms[item],item) for item in parms]       #these come first
#             parmkeys = [int(item) for item in RMCPdict['ParmNames']]
#             parmkeys.sort()
#             tempList += ['%s-%s'%(item,RMCPdict['ParmNames'][item]) for item in parmkeys]
#             print('result dscale: ',parmDict['1'],' Rw: ',Rwp)
#             atParms = [str(i+21) for i in range(len(G2Names))]
#             varyList = []
#             for item in tempList:
#                 pid = item.split('-')[0]
#                 if pid in atParms:
#                     item = '%s-%s'%(pid,G2Names[int(pid)-21])
#                 varyList.append(item)
#             result = np.array(list(newParms.values())).T
#             SeqResult[PDFfile[0]] = {'variables':result[0],'varyList':varyList,'sig':result[1],'Rvals':{'Rwp':Rwp,},
#                 'covMatrix':[],'title':PDFfile[0],'parmDict':parmDict}

#             pfile = open('Sequential_PDFfit.fgr')
#             XYcalc = np.loadtxt(pfile).T[:2]
#             pfile.close()
#             pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PDFfile[0])
#             PDFctrl = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'PDF Controls'))
#             XYobs = PDFctrl['G(R)'][1]
#             if XYobs.shape[0] < 4:
#                 XYobs = np.concatenate((XYobs,np.zeros_like(XYobs)),axis=0)
#             ibeg = np.searchsorted( XYobs[0],XYcalc[0][0])
#             ifin = ibeg+XYcalc.shape[1]
#             XYobs[2][ibeg:ifin] = XYcalc[1]
#             XYobs[3] = XYobs[1]-XYobs[2]
#             PDFctrl['G(R)'][1] = XYobs
#             SeqResult['histNames'].append(Item[1])
#             GoOn = pgbar.Update(itm,newmsg='PDF G(R) done = %d'%(itm))
#             if not GoOn[0]:
#                 print(' Sequential PDFfit aborted')
#                 break

#         pgbar.Destroy()
#         G2frame.GPXtree.SetItemPyData(Id,SeqResult)
#         G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
#         G2frame.GPXtree.SelectItem(Id)

#     else: #normal
#         #remove any old PDFfit output files
#         fName = generalData['Name'].replace(' ','_')+'-PDFfit'
#         if os.path.isfile(fName+'.res'):
#             os.remove(fName+'.res')
#         if os.path.isfile(fName+'.rstr'):
#             os.remove(fName+'.rstr')
#         if os.path.isfile(fName+'N.fgr'):
#             os.remove(fName+'N.fgr')
#         if os.path.isfile(fName+'X.fgr'):
#             os.remove(fName+'X.fgr')

#         if sys.platform.lower().startswith('win'):
#             Proc = subp.Popen('pdffit2.bat',creationflags=subp.CREATE_NEW_CONSOLE)
#             Proc.wait()     #for it to finish before continuing on
#         else:
#             if sys.platform == "darwin":
#                 GSASIIpath.MacRunScript(os.path.abspath('pdffit2.sh'))
#             else:
#                 Proc = subp.Popen(['/bin/bash','pdffit2.sh'])
#                 Proc.wait()     #for it to finish before continuing on
#         #update choice? here?
#         dlg = wx.MessageDialog(G2frame,'Check PDFfit console for results; do you want to update?',
#             'PDFfit run finished',wx.YES|wx.NO)
#         try:
#             dlg.CenterOnParent()
#             result = dlg.ShowModal()
#         finally:
#             dlg.Destroy()
#         if result == wx.ID_YES:
#             Error =  G2pwd.UpdatePDFfit(data,RMCPdict)
#             if Error:
#                 wx.MessageBox('PDFfit failed',caption='%s not found'%Error[0],style=wx.ICON_EXCLAMATION)
#         UpdateRMC()

# def Runfullrmc(event):
#     fullrmc_exec = G2pwd.findfullrmc()
#     if fullrmc_exec is None:
#         G2G.G2MessageBox(G2frame,'fullrmc Python not found. How did we get here?')
#         return
#     generalData = data['General']
#     pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
#     pName = pName.replace(' ','_')
#     rname = pName+'-fullrmc.py'
#     if not os.path.exists(rname):
#         G2G.G2MessageBox(G2frame,'The fullrmc script has not been created. Running setup.',
#             'Not setup')
#         OnSetupRMC(event)
#     RMCPdict = data['RMC']['fullrmc']
#     rmcname = pName+'-fullrmc.rmc'
#     if os.path.isdir(rmcname) and RMCPdict['ReStart'][0]:
#         msg = '''You have asked to start a new fullrmc run rather than
#              continue the existing {} run.
#              %%Press "Yes" to continue, deleting this
#              previous run or "No" to change the restart checkbox to
#              continue from the previous results.'''.format(rmcname)

#         dlg = wx.MessageDialog(G2frame,G2G.StripIndents(msg,True),
#             'Restart or continue',wx.YES|wx.NO)
#         try:
#             dlg.CenterOnParent()
#             result = dlg.ShowModal()
#         finally:
#             dlg.Destroy()
#         if result == wx.ID_YES:
#             shutil.rmtree(rmcname)
#         else:
#             return
#     G2G.G2MessageBox(G2frame,'For use of fullrmc, please cite:\n\n'+
#                          G2G.GetCite('fullrmc')+
#                          '\n\nNote: A more advanced version of fullrmc can be found at www.fullrmc.com',
#                          'Please cite fullrmc')
#     ilog = 0
#     while True:
#         logname = '%s_%d.log'%(pName,ilog)
#         if os.path.isfile(logname):
#             if GSASIIpath.GetConfigValue('debug'):
#                 print('removing',logname)
#             os.remove(logname)
#         else:
#             break
#         ilog += 1
#     if sys.platform.lower().startswith('win'):
#         batch = open('fullrmc.bat','w')
#         #batch.write('CALL '+sys.exec_prefix+'\\Scripts\\activate\n')
#         batch.write(fullrmc_exec+' '+rname+'\n')
#         batch.write('pause')
#         batch.close()
#         Proc = subp.Popen('fullrmc.bat',creationflags=subp.CREATE_NEW_CONSOLE)
#         Proc.wait()     #for it to finish before continuing on
#     else:
#         batch = open('fullrmc.sh','w')
#         batch.write('#!/bin/bash\n')
#         #activate = os.path.split(os.environ.get('CONDA_EXE',''))[0] +'/activate'
#         batch.write('cd ' + os.path.split(os.path.abspath(rname))[0] + '\n')
#         #if os.path.exists(activate):
#         #    batch.write('source ' + activate + ' ' +
#         #                os.environ['CONDA_DEFAULT_ENV'] +'\n')
#         #    batch.write('python ' + rname + '\n')
#         #else:
#         #    batch.write(sys.exec_prefix+'/python ' + rname + '\n')
#         batch.write(fullrmc_exec + ' ' + os.path.abspath(rname) + '\n')
#         batch.close()
#         if sys.platform == "darwin":
#             GSASIIpath.MacRunScript(os.path.abspath('fullrmc.sh'))
#         else:
#             Proc = subp.Popen(['/bin/bash','fullrmc.sh'])
# #                Proc.wait()     #for it to finish before continuing on
#     UpdateRMC()

# def RunRMCProfile(event):
#     generalData = data['General']
#     pName = generalData['Name'].replace(' ','_')
#     rmcfile = G2pwd.findrmcprofile()
    
#     wx.MessageBox(
#         ' For use of RMCProfile, please cite:\n\n'+
#         G2G.GetCite("RMCProfile"),
#         caption='RMCProfile',style=wx.ICON_INFORMATION)
#     if os.path.isfile(pName+'.his6f'):
#         os.remove(pName+'.his6f')
#     if os.path.isfile(pName+'.xray'):
#         os.remove(pName+'.xray')
#     if os.path.isfile(pName+'.neigh'):
#         os.remove(pName+'.neigh')
#     if os.path.isfile(pName+'.bonds'):
#         os.remove(pName+'.bonds')
#     if os.path.isfile(pName+'.triplets'):
#         os.remove(pName+'.triplets')
#     i = 1
#     while True:
#         if os.path.isfile(pName+'.bondodf_%d'%i):
#             os.remove(pName+'.bondodf_%d'%i)
#             os.remove(pName+'_bondplot_%d.ppm'%i)
#             i += 1
#         else:
#             break
#     i = 1
#     while True:
#         if os.path.isfile(pName+'_anglehist_%d.csv'%i):
#             os.remove(pName+'_anglehist_%d.csv'%i)
#             i += 1
#         else:
#             break

#     G2frame.OnFileSave(event)
#     print ('GSAS-II project saved')
#     pName = generalData['Name'].replace(' ','_')

#     if rmcfile is None:
#         wx.MessageBox('''RMCProfile is not correctly installed for use in GSAS-II
#     This software must be downloaded separately (from 
#     https://rmcprofile.ornl.gov/download). Install the rmcprofile or 
#     rmcprofile.exe file in a location where GSAS-II can find it 
#     (see config variable rmcprofile_exec in preferences.)''',
#         caption='RMCProfile',style=wx.ICON_INFORMATION)
#         return

#     if sys.platform == "darwin":
#         script_file = os.path.join(os.getcwd(), "runrmc.sh")
#         with open(script_file, 'w') as f:
#             f.write("#!/bin/bash\n")
#             f.write(f'cd "{os.getcwd()}"\n')
#             f.write(f'export PATH="{os.path.dirname(rmcfile)}":$PATH\n')
#             f.write(f'"{rmcfile}" "{pName}"\n')
#         os.system("chmod +x runrmc.sh")
#         ascript_file = os.path.join(os.getcwd(), "runrmc.script")
#         with open(ascript_file, 'w') as f:
#             f.write('tell application "Terminal"\n')
#             f.write(f'''  do script "echo 'Running RMCprofile'"\n''')
#             f.write(f'  do script "bash {script_file}" in window 1\n')
#             f.write("end tell\n")
#         subp.Popen(['osascript', ascript_file])
#     else:
#         script_file = os.path.join(os.getcwd(), "runrmc.bat")
#         with open(script_file,'w') as batch:
#             batch.write('Title RMCProfile\n')   # BHT: is Title a Windows command?
#             batch.write(f'"{rmcfile}" "{pName}"\n')
#             batch.write('pause\n')
#             batch.close()
#         subp.Popen(script_file,creationflags=subp.CREATE_NEW_CONSOLE)
# #        Proc.wait()     #for it to finish before continuing on
#     UpdateRMC()

# def OnRunRMC(event):
#     '''Run a previously created RMCProfile/fullrmc/PDFfit2 script
#     '''
#     if G2frame.RMCchoice == 'fullrmc':
#          Runfullrmc(event)
#     elif G2frame.RMCchoice == 'RMCProfile':
#         RunRMCProfile(event)
#     elif G2frame.RMCchoice == 'PDFfit':
#         RunPDFfit(event)

# # def OnStopRMC(event):
# #     if G2frame.RMCchoice == 'fullrmc':
# #         generalData = data['General']
# #         pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
# #         pName = pName.replace(' ','_')
# #         engineFilePath = pName+'.rmc'
# #         if not os.path.exists(engineFilePath):
# #             print('fullrmc repository {} not found'.format(engineFilePath))
# #             return
# #         try:
# #             from fullrmc import InterceptHook
# #             hook = InterceptHook(path=engineFilePath)
# #             hook.stop_engine()
# #             print('hook.stop_engine() sent to {}'.format(engineFilePath))
# #         except Exception as msg:
# #             print('failed, msg=',msg)

# def OnLoadRMC(event):
#     '''Used to load the output from fullrmc with all atoms placed in the
#     original cell
#     '''
#     fullrmcLoadPhase(super=False)
# def OnLoadRMCsuper(event):
#     '''Used to load the output from fullrmc with atoms in the simulation
#     supercell cell
#     '''
#     fullrmcLoadPhase(super=True)
# def fullrmcLoadPhase(super):
#     '''Used to load the output from fullrmc. Creates a new phase,
#     reads all atoms & converts coordinates to fractional.
#     If super is False all atoms placed in the original cell.

#     Called from :func:`OnLoadRMC` or :func:`OnLoadRMCsuper` from
#     the RMC tab Operations menu commands 'Superimpose into cell'
#     and 'Load Supercell'.
#     '''
#     if G2frame.RMCchoice != 'fullrmc':
#         print('fullrmcLoadPhase: How did this happen?')
#         return
#     # create a new phase
#     phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
#     phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
#     phaseNameList = list(usedHistograms.keys()) # phase names in use
#     PhaseName = G2obj.MakeUniqueLabel(data['General']['Name']+'_fullrmc',phaseNameList)
#     psub = G2frame.GPXtree.AppendItem(parent=phId,text=PhaseName)
#     E,SGData = G2spc.SpcGroup('P 1')
#     G2frame.GPXtree.SetItemPyData(psub,G2obj.SetNewPhase(Name=PhaseName,SGData=SGData))
#     newPhase = G2frame.GPXtree.GetItemPyData(psub)
#     # read in info from file
#     pName = (G2frame.GSASprojectfile.split('.')[0] + '-'
#                  + data['General']['Name'])
#     pName = pName.replace(' ','_')
#     try:
#         with open(pName+'-fullrmc.atoms','r') as fp:
#             cell = [float(i) for i in fp.readline().split(':')[1].split()]
#             supercell = [int(i) for i in fp.readline().split(':')[1].split()]
#             if super:
#                 for i in range(3): cell[i] *= supercell[i]
#             # set cell & volume
#             newPhase['General']['Cell'][1:7] = cell
#             newPhase['General']['Cell'][7] = G2lat.calc_V(G2lat.cell2A(cell))
#             A,B = G2lat.cell2AB(cell)
#             # add atoms
#             for line in fp:
#                 Name,El,ox,oy,oz = line.split()
#                 oxyz = [float(i) for i in (ox,oy,oz)]
#                 if super:
#                     (x,y,z) = np.inner(B,oxyz)
#                 else:
#                     (x,y,z),disp = G2spc.MoveToUnitCell(np.inner(B,oxyz))
#                 atId = ran.randint(0,sys.maxsize)
#                 newPhase['Atoms'].append([Name,El,'',x,y,z,1.,'1',1,'I',0.01,0,0,0,0,0,0,atId])
#     except:
#         G2G.G2MessageBox(G2frame,
#                 'Unable to open or read file '+pName+'-fullrmc.atoms. '
#                 'Was a fullrmc run from the current .gpx file '
#                 'and for the current phase?',
#                 'Error on read')
#         G2frame.GPXtree.Delete(psub)
#         return
#     # add a restraint tree entry for new phase
#     subr = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
#     G2frame.GPXtree.GetItemPyData(subr).update({PhaseName:{}})
#     SetupGeneral()  # index elements

#     #wx.CallAfter(G2frame.GPXtree.SelectItem,psub) # should call SelectDataT

# def OnViewRMC(event):
#     if G2frame.RMCchoice == 'fullrmc':
#         RMCPdict = data['RMC']['fullrmc']
#         generalData = data['General']
#         pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
#         pName = pName.replace(' ','_')
#         engineFilePath = pName+'-fullrmc.rmc'
#         if not os.path.exists(engineFilePath):
#             dlg = wx.FileDialog(G2frame, 'Open fullrmc directory',
#                 defaultFile='*.rmc',wildcard='*.rmc')
#             try:
#                 if dlg.ShowModal() == wx.ID_OK:
#                     engineFilePath = dlg.GetPath()
#                 else:
#                     return
#             finally:
#                 dlg.Destroy()
#             engineFilePath = os.path.splitext(engineFilePath)[0] + '.rmc'
#             if not os.path.exists(engineFilePath): return
#         choices = []
#         statFilePath = os.path.splitext(engineFilePath)[0] + '.stats'
#         plotFilePath = os.path.splitext(engineFilePath)[0] + '.plots'
#         imgDict = {}
#         if os.path.exists(statFilePath):
#             fp = open(statFilePath,'r')
#             vals = []
#             for i,line in enumerate(fp):
#                 v = line.strip().split(',')[:-1] # ends with comma, remove last empty element
#                 if i == 0:
#                     lbls = [i.strip() for i in v[1:]]
#                     continue
#                 try:
#                     vals.append([float(i) for i in v])
#                 except:
#                     print('Error reading line ',i,'in',statFilePath)
#             fp.close()
#             steps = np.array(vals)[:,0]
#             yvals = np.array(vals)[:,1:].T
#             choices = ['Constraints vs. Steps']
#         if os.path.exists(plotFilePath):
#             import pickle
#             fp = open(plotFilePath,'rb')
#             while True:
#                 try:
#                     title = pickle.load(fp)
#                     imgDict[title] = fp.tell()
#                     im = pickle.load(fp)
#                     choices += [title]
#                 except:
#                     break
#             fp.close()
#         if not choices:
#             G2G.G2MessageBox(G2frame,
#                 'Nothing to plot. '+
#                 'No results in '+statFilePath+' or '+plotFilePath,
#                 'Nothing to plot')
#             return
#         dlg = G2G.G2MultiChoiceDialog(G2frame,'Select plots to see displayed.',
#                                           'Select plots',choices)
#         try:
#             result = dlg.ShowModal()
#             if result == wx.ID_OK:
#                 selectedPlots = [choices[i] for i in dlg.GetSelections()]
#             else:
#                 return
#         finally:
#             dlg.Destroy()

#         for plt in selectedPlots:
#             if plt in imgDict:
#                 fp = open(plotFilePath,'rb')
#                 fp.seek(imgDict[plt])
#                 try:
#                     im = pickle.load(fp)
#                     G2plt.PlotRawImage(G2frame,im,plt)
#                 except:
#                     pass
#                 fp.close()
#             else:
#                 plotLbls = []
#                 plotVals = []
#                 for lbl,row in zip(lbls,yvals): # deal with <=0 costs
#                     if sum(row**2) == 0: continue # drop if all zeros
#                     if min(row) <= 0:
#                         row = np.where(row>0,row,min(row[np.where(row>0)])/10.)
#                     plotLbls.append(lbl)
#                     plotVals.append([steps,np.log10(row)])
#                 title = 'fullrmc residuals for '+pName
#                 G2plt.PlotXY(G2frame,plotVals,
#                             labelX='generated steps',
#                             labelY=r'$log_{10}$ ($\mathsf{\chi^2})$',
#                             newPlot=True,Title=title,
#                             lines=True,names=plotLbls)
#         return
#     elif G2frame.RMCchoice == 'RMCProfile':
#         generalData = data['General']
#         RMCPdict = data['RMC']['RMCProfile']
#         pName = generalData['Name'].replace(' ','_')
#         dlg = wx.FileDialog(G2frame, "Choose any RMCProfile csv results file for "+pName+":",
#             defaultDir=G2frame.LastGPXdir,style=wx.FD_CHANGE_DIR,wildcard='RMCProfile result csv files|'+pName+'*.csv')
#         if dlg.ShowModal() == wx.ID_OK:
#             path = os.path.split(dlg.GetPath())[0]
#             dlg.Destroy()
#         else:
#             dlg.Destroy()
#             return

#         ifXray = False
#         ifNeut = False
#         try:
#             datFile = open(os.path.join(path,pName+'.dat'),'r')
#             datLines = datFile.readlines()
#             datFile.close()
#             for line in datLines:
#                 if 'xray' in line:
#                     ifXray = True
#         except:
#             pass
#         files =  {'_PDF1.csv':[],'_PDF2.csv':[],'_PDFpartials.csv':[],'_SQ1.csv':[],'_SQ2.csv':[],'_XFQ1.csv':[],
#                   '_SQ1partials.csv':[],'_SQ2partials.csv':[],'_FQ1.csv':[],'_FT_XFQ1.csv':[],
#                   '_FQ1partials.csv':[],'_bragg.csv':[],'.chi2':[]}
#         for item in files:
#             if os.path.exists(os.path.join(path,pName+item)):
#                 OutFile = open(pName+item,'r')
#                 files[item] = OutFile.readlines()
#                 OutFile.close()
#                 print('RMCProfile file %s read'%(pName+item))
# #                else:
# #                    print('RMCProfile file %s not found'%(pName+item))
# #total result plots
#         Labels = {'_PDF1.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R) for '],
#             '_PDF2.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R)-2 for '],
#             '_SQ1.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q) for '],
#             '_SQ2.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q)-2 for '],
#             '_FQ1.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP x-ray S(Q) for '],
#             '_FT_XFQ1.csv':[r'$\mathsf{R,\AA}$','g(R)','RMCP x-ray g(R) for '],
#             '_bragg.csv':[r'$\mathsf{TOF,\mu s}$','Normalized Intensity','RMCP bragg for ']}
#         Ysave = []
#         for label in Labels:
#             X = []
#             Yobs = []
#             Ycalc = []
#             if len(files[label]):
#                 if 'PDF1' in label:
#                     ifNeut = True
#                 Names = files[label][0][:-1].split(',')
#                 Xmax = 100.
#                 if 'XFQ' in label:
#                     Xmax = RMCPdict.get('Rmax',100.)
#                 for line in files[label][1:]:
#                     items = line.split(',')
#                     if 'XFQ' in label and float(items[0]) > Xmax:
#                         break
#                     X.append(float(items[0]))
#                     Yobs.append(float(items[2]))
#                     Ycalc.append(float(items[1]))
#                 Yobs = np.array([X,Yobs])
#                 Ycalc = np.array([X,Ycalc])
#                 if '(R)' in Labels[label][1]:
#                     Ysave.append(Ycalc)
#                     Ymin = Ysave[0][1][0]
#                 if 'bragg' in label:
#                     Ydiff = np.array([X,(Yobs-Ycalc)[1]])
#                     Yoff = np.max(Ydiff[1])-np.min(Yobs[1])
#                     Ydiff[1] -= Yoff
#                     if ifXray:
#                         Labels[label][0] = r'$\mathsf{2\theta ,deg}$'
#                         Labels[label][1] = 'Intensity'
#                     G2plt.PlotXY(G2frame,[Yobs,Ycalc],XY2=[Ydiff,],labelX=Labels[label][0],
#                         labelY=Labels[label][1],newPlot=True,Title=Labels[label][2]+pName,
#                         lines=True,names=Names[1:])
#                 else:
#                     G2plt.PlotXY(G2frame,[Ycalc,Yobs],labelX=Labels[label][0],
#                         labelY=Labels[label][1],newPlot=True,Title=Labels[label][2]+pName,
#                         lines=True,names=Names[1:])
#                     RMCPdict[pName+label] = np.sum(Ycalc[1])/np.sum(Yobs[1])
#                     print(' %s scale Ycalc/Yobs: %.4f'%(label,RMCPdict[pName+label]))
# #partials plots
#         Labels = {'_PDFpartials.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R) partials for '],
#             '_SQ1partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q) partials for '],
#             '_SQ2partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)-2','RMCP S(Q) partials for '],
#             '_FQ1partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','F(Q)','RMCP F(Q) partials for ']}
#         for label in Labels:
#             X = []
#             Partials = []
#             if len(files[label]):
#                 Names = files[label][0][:-1].split(',')
#                 for line in files[label][1:]:
#                     items = line.split(',')[:-1]
#                     X.append(float(items[0]))
#                     Partials.append([float(item) for item in items[1:]])
#                 X = np.array(X)
#                 DX = X[1]-X[0]  #should be 0.02
#                 Partials = np.array(Partials).T
#                 if 'Q' in label:
#                     continue            #skip these partials
#                     XY = [[X.T,Y.T] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
#                 else:
#                     if ifNeut:
#                         XY = [[X.T,(DX*Y.T)] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
#                     else:
#                         XY = [[X.T,(DX*Y.T)*X.T] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
#                 Names = [name for name in Names if 'Va' not in name]
#                 ylabel = Labels[label][1]
#                 if 'G(R)' in Labels[label][1]:
#                     if ifNeut:
#                         title = 'Neutron '+Labels[label][2]+pName
#                     else:
#                         continue        #skip for now - x-ray partials are missing header record
#                         title = 'X-ray '+Labels[label][2].replace('G','g')+pName
#                         ylabel = 'g(R)'
#                     sumAtm = 0
#                     BLtables = G2elem.GetBLtable(generalData)
#                     AtNum = generalData['NoAtoms']
#                     if ifNeut:
#                         for atm in AtNum: sumAtm += AtNum[atm]
#                     else:
#                         for atm in AtNum:
#                             sumAtm += AtNum[atm]*G2elem.GetFFtable([atm,])[atm]['Z']
#                     bfac = {}
#                     bcorr = []
#                     for atm in AtNum:
#                         if ifNeut:
#                             if 'SL' in BLtables[atm][1]:
#                                 bfac[atm] = 10.*BLtables[atm][1]['SL'][0]*AtNum[atm]/sumAtm  #scale to pm
#                             else:   #resonant scatters (unlikely!)
#                                 bfac[atm] = AtNum[atm]/sumAtm
#                         else:
#                             bfac[atm] = 10.*G2elem.GetFFtable([atm,])[atm]['Z']*AtNum[atm]/sumAtm
#                     for name in Names:
#                         if '-' in name:
#                             at1,at2 = name.strip().split('-')
#                             if 'Va' in name:
#                                 bcorr.append(0.)
#                             else:
#                                 bcorr.append(bfac[at1]*bfac[at2])
#                             if at1 == at2:
#                                 bcorr[-1] /= 2.         #no double counting
#                     for ixy,xy in enumerate(XY):
#                         xy[1] *= bcorr[ixy]
#                         xy[1] += Ymin
#                     Xmax = np.searchsorted(Ysave[0][0],XY[0][0][-1])
#                     G2plt.PlotXY(G2frame,XY2=XY,XY=[Ysave[0][:,0:Xmax],],labelX=Labels[label][0],
#                         labelY=ylabel,newPlot=True,Title=title,
#                         lines=False,names=[r'   $G(R)_{calc}$',]+Names[1:])
#                 else:
#                     G2plt.PlotXY(G2frame,XY,labelX=Labels[label][0],
#                         labelY=ylabel,newPlot=True,Title=Labels[label][2]+pName,
#                         lines=True,names=Names[1:])
# #chi**2 plot
#         X = []
#         Chi = []
#         if len(files['.chi2']) > 2:
#             Names = files['.chi2'][0][:-1].split()
#             for line in files['.chi2'][1:]:
#                 items = line[:-1].split()
#                 X.append(float(items[1]))
#                 Chi.append([float(item) for item in items[3:]])
#             X = np.array(X)
#             Chi = np.array(Chi).T
#             XY = [[X.T,np.log10(Y.T)] for Y in Chi]
#             G2plt.PlotXY(G2frame,XY,labelX='no. generated',
#                 labelY=r'$log_{10}$ (reduced $\mathsf{\chi^2})$',newPlot=True,Title='RMCP Chi^2 for '+pName,
#                 lines=True,names=Names[3:])

# #get atoms from rmc6f file
#         rmc6fName = pName+'.rmc6f'
#         rmc6f = open(rmc6fName,'r')
#         rmc6fAtoms = []
#         while True:
#             line = rmc6f.readline()
#             if 'Number of atoms:' in line:
#                 Natm = int(line.split(':')[1])
#             if 'Atoms:' in line:
#                 break
#         for iAtm in range(Natm):
#             line = rmc6f.readline().split()
#             rmc6fAtoms.append([line[1],float(line[3]),float(line[4]),float(line[5])])
#         rmc6f.close()
# #alt bond histograms - from rmc6 & bond files

#         bondName = pName+'.bonds'
#         if os.path.exists(os.path.join(path,bondName)):
#             nBond = len(RMCPdict['Potentials']['Stretch'])
#             bondList = []
#             bonds = open(bondName,'r')
#             while True:
#                 line = bonds.readline()
#                 if '............' in line:
#                     break       #positions at start of 1st bond list
#             for iBond in range(nBond):
#                 bondList.append([])
#                 Bonds = RMCPdict['Potentials']['Stretch'][iBond]
#                 for iAtm in range(Natm):     #same as in rmc6f above
#                     line = bonds.readline().split('::')
#                     if Bonds[0] in line[0]:
#                         items = line[1].replace(';','').split()[:-1:2]
#                         items = np.array([int(item) for item in items])
#                         bondList[iBond].append([iAtm,items])
#             bonds.close()
#             for iBond in range(nBond):
#                 Bonds = RMCPdict['Potentials']['Stretch'][iBond]
#                 title = '%s-%s'%(Bonds[0],Bonds[1])
#                 bondDist = G2pwd.GetRMCBonds(generalData,RMCPdict,rmc6fAtoms,bondList[iBond])
#                 print('%s mean %.3f(%d)'%(title,np.mean(bondDist),int(1000*np.std(bondDist))))
#                 G2plt.PlotBarGraph(G2frame,bondDist,Xname=r'$Bond, \AA$',Title=title+' from Potential Energy Restraint',
#                     PlotName='%s Bond for %s'%(title,pName))
#                 print(' %d %s bonds found'%(len(bondDist),title))

# #alt angle histograms - from rmc6 & triplets files
#         tripName = pName+'.triplets'
#         if os.path.exists(os.path.join(path,tripName)):
#             nAng = len(RMCPdict['Potentials']['Angles'])
#             tripList = []
#             triples = open(tripName,'r')
#             while True:
#                 line = triples.readline()
#                 if '............' in line:
#                     break       #positions at start of 1st triple list
#             for iAng in range(nAng):
#                 tripList.append([])
#                 Angles = RMCPdict['Potentials']['Angles'][iAng]
#                 for iAtm in range(Natm):     #same as in rmc6f above
#                     line = triples.readline().split('::')
#                     if Angles[1] in line[0]:
#                         items = line[1].replace(';','').split()[:-1:2]
#                         items = np.array([int(item) for item in items]).reshape((-1,2))
#                         tripList[iAng].append([iAtm,items])
#                     line = triples.readline()       #to skip a line
#             triples.close()
#             for iAng in range(nAng):
#                 Angles = RMCPdict['Potentials']['Angles'][iAng]
#                 angles = G2pwd.GetRMCAngles(generalData,RMCPdict,rmc6fAtoms,tripList[iAng])
#                 title = '%s-%s-%s'%(Angles[0],Angles[1],Angles[2])
#                 print('%s mean %.2f(%d)'%(title,np.mean(angles),int(100*np.std(angles))))
#                 G2plt.PlotBarGraph(G2frame,angles,Xname=r'$Angle, \AA$',Title=title+' from Potential Energy Restraint',
#                     PlotName='%s Angle for %s'%(title,pName))
#                 print(' %d %s angles found'%(len(angles),title))

# #bond odf plots
#         nPot = len(RMCPdict['Potentials']['Stretch'])
#         for iPot in range(nPot):
#             fname = pName+'.bondodf_%d'%(iPot+1)
#             bond = RMCPdict['Potentials']['Stretch'][iPot]
#             if os.path.exists(os.path.join(path,fname)):
#                 OutFile = open(fname,'r')
#                 odfFile = OutFile.readlines()
#                 if len(odfFile) > 1:
#                     OutFile.seek(0)
#                     odfData = np.fromfile(OutFile,sep=' ')
#                     numx,numy = odfData[:2]
#                     G2plt.Plot3dXYZ(G2frame,int(numx),int(numy),odfData[2:],
#                         newPlot=False,Title='Number of %s-%s Bonds'%(bond[0],bond[1]),Centro=True)
#                 OutFile.close()
#     elif G2frame.RMCchoice == 'PDFfit':
#         generalData = data['General']
#         RMCPdict = data['RMC']['PDFfit']
#         pName = generalData['Name'].replace(' ','_')
#         Labels = [r'$\mathsf{R,\AA}$','G(R)','PDFfit2 G(R) for ']
#         files = RMCPdict['files']
#         for file in files:
#             if 'Select' not in files[file][0]:
#                 XY = np.empty((1,2))
#                 start = 0
#                 while XY.shape[0] == 1:
#                     try:
#                         XY = np.loadtxt(files[file][0],skiprows=start)
#                     except ValueError:
#                         start += 1
#                 if 'Neutron' in file:
#                     cname = pName+'-PDFfitN.fgr'
#                 else:
#                     cname = pName+'-PDFfitX.fgr'
#                 XYobs = XY.T[:2]
#                 XY = np.empty((1,2))
#                 start = 0
#                 while XY.shape[0] == 1:
#                     try:
#                         XY = np.loadtxt(cname,skiprows=start)
#                     except ValueError:
#                         start += 1
#                 XYcalc = XY.T[:2]
#                 XYdiff = np.zeros_like(XYcalc)
#                 XYdiff[0] = XYcalc[0]
#                 ibeg = np.searchsorted(XYobs[0],XYdiff[0][0])
#                 ifin = ibeg+XYcalc.shape[1]
#                 XYdiff[1] = XYobs[1][ibeg:ifin]-XYcalc[1]
#                 offset = 1.1*(np.max(XYdiff[1])-np.min(XYcalc[1]))
#                 XYdiff[1] -= offset
#                 G2plt.PlotXY(G2frame,[XYobs,],XY2=[XYcalc,XYdiff],labelX=Labels[0],
#                     labelY=Labels[1],newPlot=True,Title=Labels[2]+files[file][0],
#                     lines=False,names=['G(R) obs','G(R) calc','diff',])


# #### ISODISTORT tab ###############################################################################

# def UpdateISODISTORT(Scroll=0):
#     ''' Setup ISODISTORT and present the results. Allow selection of a distortion model for PDFfit or
#     GSAS-II structure refinement as a cif file produced by ISODISTORT. Allows manipulation of distortion
#     mode displacements selection their refinement for this new phase.
#     '''

#     def displaySetup():

#         def OnParentCif(event):
#             dlg = wx.FileDialog(ISODIST, 'Select parent cif file',G2frame.LastGPXdir,
#                 style=wx.FD_OPEN ,wildcard='cif file(*.cif)|*.cif')
#             if dlg.ShowModal() == wx.ID_OK:
#                 fName = dlg.GetFilename()
#                 fDir = dlg.GetDirectory()
#                 ISOdata['ParentCIF'] = os.path.join(fDir,fName)
#                 dlg.Destroy()
#             else:
#                 dlg.Destroy()
#             UpdateISODISTORT()

#         def OnUsePhase(event):
#             ISOdata['ParentCIF'] = 'Use this phase'
#             UpdateISODISTORT()

#         def OnMethodSel(event):
#             method = methodSel.GetSelection()+1
#             if method in [1,4]:
#                 ISOdata['ISOmethod'] = method
#             UpdateISODISTORT()

#         def OnChildCif(event):
#             dlg = wx.FileDialog(ISODIST, 'Select child cif file',G2frame.LastGPXdir,
#                 style=wx.FD_OPEN ,wildcard='cif file(*.cif)|*.cif')
#             if dlg.ShowModal() == wx.ID_OK:
#                 fName = dlg.GetFilename()
#                 fDir = dlg.GetDirectory()
#                 ISOdata['ChildCIF'] = os.path.join(fDir,fName)
#                 dlg.Destroy()
#             else:
#                 dlg.Destroy()
#             UpdateISODISTORT()

#         def OnUsePhase2(event):
#             ISOdata['ChildCIF'] = 'Use this phase'
#             UpdateISODISTORT()

#         topSizer = wx.BoxSizer(wx.VERTICAL)
#         topSizer.Add(wx.StaticText(ISODIST,label=' ISODISTORT setup controls:'))
#         parentSizer = wx.BoxSizer(wx.HORIZONTAL)
#         parentSizer.Add(wx.StaticText(ISODIST,label=' Parent cif file:'),0,WACV)
#         parentCif = wx.Button(ISODIST,label=ISOdata['ParentCIF'],size=(300,24))
#         parentCif.Bind(wx.EVT_BUTTON,OnParentCif)
#         parentSizer.Add(parentCif,0,WACV)
#         if 'Use this phase' not in ISOdata['ChildCIF'] and 'Use this phase' not in ISOdata['ParentCIF']:
#             usePhase = wx.Button(ISODIST,label=' Use this phase? ')
#             usePhase.Bind(wx.EVT_BUTTON,OnUsePhase)
#             parentSizer.Add(usePhase,0,WACV)
#         topSizer.Add(parentSizer)
#         #patch
#         if 'ISOmethod' not in ISOdata:
#             ISOdata['ISOmethod'] = 1
#         #end patch
#         choice = ['Method 1: Search over all special k points - yields only single Irrep models',
#             'Method 2: not implemented in GSAS-II',
#             'Method 3: not implemented in GSAS-II',
#             'Method 4: Mode decomposition of known child structure']
#         methodSel = wx.RadioBox(ISODIST,label='Select ISODISTORT method:',choices=choice,
#             majorDimension=1,style=wx.RA_SPECIFY_COLS)
#         methodSel.SetSelection(ISOdata['ISOmethod']-1)
#         methodSel.Bind(wx.EVT_RADIOBOX,OnMethodSel)
#         topSizer.Add(methodSel)
#         if ISOdata['ISOmethod'] == 4:
#             childSizer = wx.BoxSizer(wx.HORIZONTAL)
#             childSizer.Add(wx.StaticText(ISODIST,label=' Child cif file:'),0,WACV)
#             childCif = wx.Button(ISODIST,label=ISOdata['ChildCIF'],size=(300,24))
#             childCif.Bind(wx.EVT_BUTTON,OnChildCif)
#             childSizer.Add(childCif,0,WACV)
#             if 'Use this phase' not in ISOdata['ChildCIF'] and 'Use this phase' not in ISOdata['ParentCIF']:
#                 usePhase2 = wx.Button(ISODIST,label=' Use this phase? ')
#                 usePhase2.Bind(wx.EVT_BUTTON,OnUsePhase2)
#                 childSizer.Add(usePhase2,0,WACV)
#             topSizer.Add(childSizer)

#         return topSizer

#     def displaySubset():

#         def OnLaue(event):
#             Obj = event.GetEventObject()
#             name = Indx[Obj.GetId()]
#             ISOdata['SGselect'][name[:4]] = not ISOdata['SGselect'][name[:4]]
#             ISOdata['selection'] = None
#             UpdateISODISTORT()

#         def OnAllBtn(event):
#             for item in ISOdata['SGselect']:
#                 ISOdata['SGselect'][item] = not ISOdata['SGselect'][item]
#             ISOdata['selection'] = None
#             UpdateISODISTORT()

#         topSizer = wx.BoxSizer(wx.VERTICAL)
#         G2G.HorizontalLine(topSizer,ISODIST)
#         topSizer.Add(wx.StaticText(ISODIST,label='ISODISTORT Method 1 distortion search results:'))
#         topSizer.Add(wx.StaticText(ISODIST,label=' Subset selection if desired:'))
#         laueName = ['Cubic','Hexagonal','Trigonal','Tetragonal','Orthorhombic','Monoclinic','Triclinic']
#         littleSizer = wx.FlexGridSizer(0,8,5,5)
#         Indx = {}
#         for name in laueName:
#             laueCk = wx.CheckBox(ISODIST,label=name)
#             Indx[laueCk.GetId()] = name
#             laueCk.SetValue(ISOdata['SGselect'][name[:4]])
#             laueCk.Bind(wx.EVT_CHECKBOX,OnLaue)
#             littleSizer.Add(laueCk,0,WACV)
#         allBtn = wx.Button(ISODIST,label='Toggle all')
#         allBtn.Bind(wx.EVT_BUTTON,OnAllBtn)
#         littleSizer.Add(allBtn)
#         topSizer.Add(littleSizer)
#         return topSizer

#     def displayRadio():

#         def CheckItem(item):
#             SGnum = int(item.split()[1].split('*')[0])
#             for SGtype in ISOdata['SGselect']:
#                 if ISOdata['SGselect'][SGtype] and SGnum in SGrange[SGtype]:
#                     return True
#             return False

#         def OnSelect(event):
#            r,c = event.GetRow(),event.GetCol()
#            if c == 0:
#                ISOdata['selection'] = [r,isoTable.GetValue(r,1)]
#                for row in range(isoGrid.GetNumberRows()):
#                    isoTable.SetValue(row,c,False)
#                isoTable.SetValue(r,c,True)
#                isoGrid.ForceRefresh()

#         SGrange = {'Cubi':np.arange(195,231),'Hexa':np.arange(168,195),'Trig':np.arange(143,168),'Tetr':np.arange(75,143),
#                    'Orth':np.arange(16,75),'Mono':np.arange(3,16),'Tric':np.arange(1,3)}
#         bottomSizer = wx.BoxSizer(wx.VERTICAL)
#         colLabels = ['select',' ISODISTORT order parameter direction description']
#         colTypes = [wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,]

#         Radio = ISOdata['radio']
#         rowLabels = []
#         table = []
#         for i,item in enumerate(Radio):
#             if CheckItem(Radio[item]):
#                 if ISOdata['selection'] and ISOdata['selection'][0] == i:
#                     table.append([True,Radio[item]])
#                 else:
#                     table.append([False,Radio[item]])
#                 rowLabels.append(str(i))
#         isoTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=colTypes)
#         isoGrid = G2G.GSGrid(ISODIST)
#         isoGrid.SetTable(isoTable,True,useFracEdit=False)
#         isoGrid.AutoSizeColumns(True)
#         isoGrid.SetColLabelAlignment(wx.ALIGN_LEFT,wx.ALIGN_CENTRE)
#         bottomSizer.Add(isoGrid)
#         attr = wg.GridCellAttr()
#         attr.SetReadOnly(True)
#         attr.SetBackgroundColour(VERY_LIGHT_GREY)
#         isoGrid.SetColAttr(1,attr)
#         isoGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, OnSelect)
#         return bottomSizer

#     def displayModes():

#         def OnDispl(event):
#             '''Respond to movement of distortion mode slider'''
#             Obj = event.GetEventObject()
#             idsp,dispVal = Indx[Obj.GetId()]
#             modeDisp[idsp] = Obj.GetValue()/1000.
#             dispVal.SetValue(modeDisp[idsp])
#             err = G2mth.ApplyModeDisp(data)
#             if err:
#                 G2G.G2MessageBox(G2frame,'Do Draw atoms first')
#             FindBondsDraw(data)
#             G2plt.PlotStructure(G2frame,data)

#         def OnDispVal(invalid,value,tc):
#             '''Respond to entry of a value into a distortion mode entry widget'''
#             idsp,displ = Indx[tc.GetId()]
#             displ.SetValue(int(value*1000))
#             err = G2mth.ApplyModeDisp(data)
#             if err:
#                 G2G.G2MessageBox(G2frame,'Do Draw atoms first')
#             FindBondsDraw(data)
#             G2plt.PlotStructure(G2frame,data)

#         def OnRefDispl(event):
#             Obj = event.GetEventObject()
#             idsp,item = Indx[Obj.GetId()]
#             item[-2] = not item[-2]

#         def OnReset(event):
#             '''Reset all distortion mode values to initial values'''
#             ISOdata['modeDispl'] = copy.deepcopy(ISOdata['ISOmodeDispl'])
#             err = G2mth.ApplyModeDisp(data)
#             if err:
#                 G2G.G2MessageBox(G2frame,'Do Draw atoms first')
#             FindBondsDraw(data)
#             G2plt.PlotStructure(G2frame,data)
#             UpdateISODISTORT()

#         def OnSetZero(event):
#             '''Reset all distortion mode values to 0'''
#             ISOdata['modeDispl'] = [0.0 for i in ISOdata['ISOmodeDispl']]
#             err = G2mth.ApplyModeDisp(data)
#             if err:
#                 G2G.G2MessageBox(G2frame,'Do Draw atoms first')
#             FindBondsDraw(data)
#             G2plt.PlotStructure(G2frame,data)
#             UpdateISODISTORT()

#         def OnSaveModes(event):
#             '''Set saved distortion mode values to displayed values'''
#             dlg = wx.MessageDialog(G2frame,'Are you sure you want to replace the saved mode values?',
#                 'Confirm replace',wx.YES|wx.NO)
#             try:
#                 dlg.CenterOnParent()
#                 result = dlg.ShowModal()
#             finally:
#                 dlg.Destroy()
#             if result != wx.ID_YES: return
#             ISOdata['ISOmodeDispl'] = copy.deepcopy(ISOdata['modeDispl'])
#             G2plt.PlotStructure(G2frame,data)
#             UpdateISODISTORT()

#         #### displayModes code starts here
#         ConstrData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Constraints'))
#         pId = data['ranId']
#         mainSizer = wx.BoxSizer(wx.VERTICAL)
#         if SGLaue not in ['mmm','2/m','-1']:
#             mainSizer.Add(wx.StaticText(ISODIST,label=' NB: ISODISTORT distortion mode symmetry is too high to be used in PDFfit'))
#         txt = wx.StaticText(ISODIST,label=
#                         ' For use of ISODISTORT, please cite: '+
#                         G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...'))
#         txt.Wrap(500)
#         mainSizer.Add(txt)
#         mainSizer.Add(wx.StaticText(ISODIST,label=
# u''' The 2nd column below shows the last saved mode values. The 3rd && 4th columns will set the
#  display mode values. The positions in the Atoms and Draw Atoms tabs, as well as the atom
#  positions shown in the Plot Window are changed to reflect the display mode values. The
#  range of the slider corresponds to making a maximum atomic displacement between -2 && +2 \u212B.'''))
#         mainSizer.Add((-1,10))
#         slideSizer = wx.FlexGridSizer(0,5,0,0)
#         modeDisp = ISOdata['modeDispl']
#         idsp = 0
#         slideSizer.Add(wx.StaticText(ISODIST,label='Name'),0,wx.ALIGN_CENTER)
#         slideSizer.Add(wx.StaticText(ISODIST,label='Save value'))
#         slideSizer.Add(wx.StaticText(ISODIST,label='Value'),0,wx.ALIGN_CENTER)
#         slideSizer.Add(wx.StaticText(ISODIST,label='Refine?'),0,wx.ALIGN_CENTER)
#         slideSizer.Add(wx.StaticText(ISODIST,label='Atom displacements'),0,wx.EXPAND|wx.LEFT,15)
#         isoDict = {i.name:j for (i,j) in zip(data['ISODISTORT']['G2ModeList'],data['ISODISTORT']['IsoModeList'])}
#         for item in ConstrData['Phase']:
#             if item[-1] != 'f': continue # only want new vars
#             if item[-3] is None: continue # unnamed new var is not ISO
#             try:
#                 if pId != item[-3].phase: continue # at present only ISO modes are associated with a phase
#             except AttributeError:
#                 continue
#             if  item[-3].name not in isoDict: continue
#             isoName = item[-3].varname().split('::')[1]
#             slideSizer.Add(wx.StaticText(ISODIST,label=isoName),0,WACV)
#             slideSizer.Add(wx.StaticText(ISODIST,label=' %.5g '%ISOdata['ISOmodeDispl'][idsp],
#                 style=wx.ALIGN_CENTER_HORIZONTAL),0,WACV|wx.EXPAND)
#             lineSizer = wx.BoxSizer(wx.HORIZONTAL)
#             dispVal = G2G.ValidatedTxtCtrl(ISODIST,modeDisp,idsp,xmin=-2.,xmax=2.,size=(75,20),OnLeave=OnDispVal)
#             lineSizer.Add(dispVal,0,WACV)
#             displ = G2G.G2Slider(ISODIST,style=wx.SL_HORIZONTAL,minValue=-2000,maxValue=2000,
#                 value=int(modeDisp[idsp]*1000),size=(250,20))
#             displ.Bind(wx.EVT_SLIDER, OnDispl)
#             Indx[displ.GetId()] = [idsp,dispVal]
#             Indx[dispVal.GetId()] = [idsp,displ]
#             lineSizer.Add(displ)
#             slideSizer.Add(lineSizer)
#             refDispl = wx.CheckBox(ISODIST)
#             refDispl.SetValue(item[-2])
#             refDispl.Bind(wx.EVT_CHECKBOX,OnRefDispl)
#             Indx[refDispl.GetId()] = [idsp,item]
#             slideSizer.Add(refDispl,0,WACV|wx.EXPAND|wx.LEFT,15)
#             slideSizer.Add(wx.StaticText(ISODIST,label=', '.join(ModeDispList[idsp])),0,wx.EXPAND|wx.LEFT,15)
#             idsp += 1
#         slideSizer.SetMinSize(wx.Size(650,10))
#         mainSizer.Add(slideSizer)
#         lineSizer = wx.BoxSizer(wx.HORIZONTAL)
#         reset = wx.Button(ISODIST,label='Reset modes to save values')
#         reset.Bind(wx.EVT_BUTTON,OnReset)
#         lineSizer.Add(reset,0,WACV)
#         reset = wx.Button(ISODIST,label='Set all modes to zero')
#         reset.Bind(wx.EVT_BUTTON,OnSetZero)
#         lineSizer.Add(reset,0,wx.ALL,10)
#         reset = wx.Button(ISODIST,label='Save mode values')
#         reset.Bind(wx.EVT_BUTTON,OnSaveModes)
#         lineSizer.Add(reset,0,WACV)
#         mainSizer.Add(lineSizer,0,wx.TOP,5)
#         mainSizer.Layout()
#         SetPhaseWindow(ISODIST,mainSizer,Scroll=Scroll)

#     #### UpdateISODISTORT code starts here
#     topSizer = G2frame.dataWindow.topBox
#     topSizer.Clear(True)
#     parent = G2frame.dataWindow.topPanel
#     lbl= f"ISODISTORT distortion modes for {data['General']['Name']!r}"[:60]
#     topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
#     topSizer.Add((-1,-1),1,wx.EXPAND)
#     topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
#     wx.CallAfter(G2frame.dataWindow.SetDataSize)
#     Indx = {}
#     ISOdata = data['ISODISTORT']
#     SGLaue = data['General']['SGData']['SGLaue']
#     G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_ISODNEWPHASE,'rundata' in ISOdata)
#     G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_ISOPDFFIT,(('G2VarList' in ISOdata) and (SGLaue in ['mmm','2/m','-1'])))
#     G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_SHOWISO1,('G2VarList' in ISOdata)
#         or ('G2OccVarList' in ISOdata))
#     G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_SHOWISOMODES,('G2VarList' in ISOdata))

#     if ISODIST.GetSizer():
#         ISODIST.GetSizer().Clear(True)

#     if 'G2ModeList' in ISOdata:      #invoked only if phase is from a ISODISTORT cif file & thus contains distortion mode constraints

# # #patch
# #             if 'modeDispl' not in ISOdata:
# #                 ISOdata['modeDispl'] = np.zeros(len(ISOdata['G2ModeList']))
# # #end patch
#         ModeDispList = G2pwd.GetAtmDispList(ISOdata)
#         displayModes()
#         return

# #initialization
#     if 'ParentCIF' not in ISOdata:
#         ISOdata.update({'ParentCIF':'Select','ChildCIF':'Select','ISOmethod':4,
#             'ChildMatrix':np.eye(3),'ChildSprGp':'P 1','ChildCell':'abc',})         #these last 3 currently unused
# #end initialization

#     mainSizer = wx.BoxSizer(wx.VERTICAL)
#     txt = wx.StaticText(ISODIST,label=
#                         ' For use of ISODISTORT, please cite: '+
#                         G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...'))
#     txt.Wrap(500)
#     mainSizer.Add(txt)
#     mainSizer.Add((-1,5))
#     G2G.HorizontalLine(mainSizer,ISODIST)
#     mainSizer.Add((-1,5))
#     mainSizer.Add(displaySetup())

#     if 'radio' in ISOdata:
#         mainSizer.Add(displaySubset())
#         mainSizer.Add(displayRadio())
#     SetPhaseWindow(ISODIST,mainSizer,Scroll=Scroll)

# def OnRunISODISTORT(event):
#     ''' this needs to setup for method #3 or #4 in ISODISTORT
#     after providing parent cif:
#     #3 asks for transformation matrix & space group of child structure
#     #4 asks for cif file of child structure
#     '''

#     if not G2frame.GSASprojectfile:     #force a project save to establish location of output cif file
#         G2frame.OnFileSaveas(event)

#     radio,rundata = ISO.GetISODISTORT(data)
#     if radio:
#         data['ISODISTORT']['radio'] = radio
#         data['ISODISTORT']['rundata'] = rundata
#         data['ISODISTORT']['SGselect'] =  {'Tric':True,'Mono':True,'Orth':True,'Tetr':True,'Trig':True,'Hexa':True,'Cubi':True}
#         data['ISODISTORT']['selection'] = None
#         print('ISODISTORT run complete')
#         wx.CallAfter(UpdateISODISTORT)
#     elif data['ISODISTORT']['ISOmethod'] != 4 or radio is None:
#         G2G.G2MessageBox(G2frame,'ISODISTORT run failed - see page opened in web browser')
#     else:
#         G2G.G2MessageBox(G2frame,'ISODISTORT run complete; new cif file %s created.\n To use, import it as a new phase.'%rundata)
#         print(' ISODISTORT run complete; new cif file %s created. To use, import it as a new phase.'%rundata)

# def OnNewISOPhase(event):
#     ''' Make CIF file with ISODISTORT
#     '''
#     if 'rundata' in data['ISODISTORT'] and data['ISODISTORT']['selection'] is not None:
#         CIFfile = ISO.GetISODISTORTcif(data)
#         G2G.G2MessageBox(G2frame,'ISODISTORT generated cif file %s has been created.'%CIFfile)
#     elif 'rundata' in data['ISODISTORT']:
#         G2G.G2MessageBox(G2frame,'Need to select an ISODISTORTdistortion model first before creating a CIF')
#     else:
#         G2G.G2MessageBox(G2frame,'ERROR - need to run ISODISTORT first - see General/Compute menu')

# def OnNewPDFfitPhase(event):
#     ''' Make new phase for PDFfit using ISODISTORT mode definitions as constraints
#     '''
#     newPhase = G2pwd.ISO2PDFfit(data)
#     phaseName = newPhase['General']['Name']
#     sub = G2frame.GPXtree.AppendItem(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases'),text=phaseName)
#     G2frame.GPXtree.SetItemPyData(sub,newPhase)
#     G2frame.GPXtree.SelectItem(sub)

# #### DIFFax Layer Data page ################################################################################
# def UpdateLayerData(Scroll=0):
#     '''Present the contents of the Phase/Layers tab for stacking fault simulation
#     '''

#     laueChoice = ['-1','2/m(ab)','2/m(c)','mmm','-3','-3m','4/m','4/mmm',
#         '6/m','6/mmm','unknown']
#     colLabels = ['Name','Type','x','y','z','frac','Uiso']
#     transLabels = ['Prob','Dx','Dy','Dz','refine','plot']
#     colTypes = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,]+ \
#         3*[wg.GRID_VALUE_FLOAT+':10,5',]+2*[wg.GRID_VALUE_FLOAT+':10,4',] #x,y,z,frac,Uiso
#     transTypes = [wg.GRID_VALUE_FLOAT+':10,3',]+3*[wg.GRID_VALUE_FLOAT+':10,5',]+ \
#         [wg.GRID_VALUE_CHOICE+": ,P,Dx,Dy,Dz,Dxy,Dxz,Dyz,Dxyz",wg.GRID_VALUE_BOOL,]
#     plotDefaults = {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':30.,'viewDir':[0,0,1],
#         'viewPoint':[[0.,0.,0.],[]],}
#     Indx = {}

#     def OnLaue(event):
#         Obj = event.GetEventObject()
#         data['Layers']['Laue'] = Obj.GetValue()
#         wx.CallAfter(UpdateLayerData)

#     def OnSadpPlot(event):
#         sadpPlot.SetValue(False)
#         labels = Layers['Sadp']['Plane']
#         lmax = float(Layers['Sadp']['Lmax'])
#         XY = 2*lmax*np.mgrid[0:256:256j,0:256:256j]/256.-lmax
#         G2frame.Cmax = 1.0
#         G2plt.PlotXYZ(G2frame,XY,Layers['Sadp']['Img'].T,labelX=labels[:-1],
#             labelY=labels[-1],newPlot=False,Title=Layers['Sadp']['Plane'])

#     def OnSeqPlot(event):
#         seqPlot.SetValue(False)
#         resultXY,resultXY2,seqNames = Layers['seqResults']
#         pName = Layers['seqCodes'][0]
#         G2plt.PlotXY(G2frame,resultXY,XY2=resultXY2,labelX=r'$\mathsf{2\theta}$',
#             labelY='Intensity',newPlot=True,Title='Sequential simulations on '+pName,
#             lines=True,names=seqNames)

#     def CellSizer():

#         cellGUIlist = [
#             [['-3','-3m','6/m','6/mmm','4/m','4/mmm'],6,zip([" a = "," c = "],["%.5f","%.5f",],[True,True],[0,2])],
#             [['mmm'],8,zip([" a = "," b = "," c = "],["%.5f","%.5f","%.5f"],[True,True,True],[0,1,2,])],
#             [['2/m(ab)','2/m(c)','-1','axial','unknown'],10,zip([" a = "," b = "," c = "," gamma = "],
#                 ["%.5f","%.5f","%.5f","%.3f"],[True,True,True,True],[0,1,2,5])]]

#         def OnCellRef(event):
#             data['Layers']['Cell'][0] = cellRef.GetValue()

#         def OnCellChange(event):
#             event.Skip()
#             laue = data['Layers']['Laue']
#             cell = data['Layers']['Cell']
#             Obj = event.GetEventObject()
#             ObjId = cellList.index(Obj.GetId())
#             try:
#                 value = max(1.0,float(Obj.GetValue()))
#             except ValueError:
#                 if ObjId < 3:               #bad cell edge - reset
#                     value = cell[ObjId+1]
#                 else:                       #bad angle
#                     value = 90.
#             if laue in ['-3','-3m','6/m','6/mmm','4/m','4/mmm']:
#                 cell[4] = cell[5] = 90.
#                 cell[6] = 120.
#                 if laue in ['4/m','4/mmm']:
#                     cell[6] = 90.
#                 if ObjId == 0:
#                     cell[1] = cell[2] = value
#                     Obj.SetValue("%.5f"%(cell[1]))
#                 else:
#                     cell[3] = value
#                     Obj.SetValue("%.5f"%(cell[3]))
#             elif laue in ['mmm']:
#                 cell[ObjId+1] = value
#                 cell[4] = cell[5] = cell[6] = 90.
#                 Obj.SetValue("%.5f"%(cell[ObjId+1]))
#             elif laue in ['2/m','-1']:
#                 cell[4] = cell[5] = 90.
#                 if ObjId != 3:
#                     cell[ObjId+1] = value
#                     Obj.SetValue("%.5f"%(cell[ObjId+1]))
#                 else:
#                     cell[6] = value
#                     Obj.SetValue("%.3f"%(cell[6]))
#             cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
#             volVal.SetLabel(' Vol = %.3f'%(cell[7]))

#         cell = data['Layers']['Cell']
#         laue = data['Layers']['Laue']
#         for cellGUI in cellGUIlist:
#             if laue in cellGUI[0]:
#                 useGUI = cellGUI
#         cellSizer = wx.FlexGridSizer(0,useGUI[1]+1,5,5)
#         cellRef = wx.CheckBox(layerData,-1,label='Refine unit cell:')
#         cellSizer.Add(cellRef,0,WACV)
#         cellRef.Bind(wx.EVT_CHECKBOX, OnCellRef)
#         cellRef.SetValue(cell[0])
#         cellList = []
#         for txt,fmt,ifEdit,Id in useGUI[2]:
#             cellSizer.Add(wx.StaticText(layerData,label=txt),0,WACV)
# #            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
#             cellVal = wx.TextCtrl(layerData,value=(fmt%(cell[Id+1])),
#                 style=wx.TE_PROCESS_ENTER)
#             cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)
#             cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
#             cellSizer.Add(cellVal,0,WACV)
#             cellList.append(cellVal.GetId())
#         volVal = wx.StaticText(layerData,label=' Vol = %.3f'%(cell[7]))
#         cellSizer.Add(volVal,0,WACV)
#         return cellSizer

#     def WidthSizer():

#         def OnRefWidth(event):
#             Id = Indx[event.GetEventObject()]
#             Layers['Width'][1][Id] = not Layers['Width'][1][Id]

#         Labels = ['a','b']
#         flags = Layers['Width'][1]
#         widthSizer = wx.BoxSizer(wx.HORIZONTAL)
#         for i in range(2):
#             widthSizer.Add(wx.StaticText(layerData,label=u' layer width(%s) (<= 1\xb5m): '%(Labels[i])),0,WACV)
#             widthVal = G2G.ValidatedTxtCtrl(layerData,Layers['Width'][0],i,nDig=(10,3),xmin=0.005,xmax=1.0)
#             widthSizer.Add(widthVal,0,WACV)
#             widthRef = wx.CheckBox(layerData,label='Refine?')
#             widthRef.SetValue(flags[i])
#             Indx[widthRef] = i
#             widthRef.Bind(wx.EVT_CHECKBOX, OnRefWidth)
#             widthSizer.Add(widthRef,0,WACV)
#         return widthSizer

#     def OnNewLayer(event):
#         data['Layers']['Layers'].append({'Name':'Unk','SameAs':'','Symm':'None','Atoms':[]})
#         Trans = data['Layers']['Transitions']
#         if len(Trans):
#             Trans.append([[0.,0.,0.,0.,'',False] for trans in Trans])
#             for trans in Trans:
#                 trans.append([0.,0.,0.,0.,'',False])
#         else:
#             Trans = [[[1.,0.,0.,0.,'',False],],]
#         data['Layers']['Transitions'] = Trans
#         wx.CallLater(100,UpdateLayerData)

#     def OnDeleteLast(event):
#         del(data['Layers']['Layers'][-1])
#         del(data['Layers']['Transitions'][-1])
#         for trans in data['Layers']['Transitions']:
#             del trans[-1]
#         wx.CallAfter(UpdateLayerData)

#     def OnImportLayer(event):
#         dlg = wx.FileDialog(G2frame, 'Choose GSAS-II project file', G2G.GetImportPath(G2frame),
#             wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN| wx.FD_CHANGE_DIR)
#         try:
#             if dlg.ShowModal() == wx.ID_OK:
#                 GPXFile = dlg.GetPath()
#                 phaseNames = G2stIO.GetPhaseNames(GPXFile)
#             else:
#                 return
#         finally:
#             dlg.Destroy()
#         dlg = wx.SingleChoiceDialog(G2frame,'Phase to use for layer','Select',phaseNames)
#         if dlg.ShowModal() == wx.ID_OK:
#             sel = dlg.GetSelection()
#             PhaseName = phaseNames[sel]
#         else:
#             return
#         Phase = G2stIO.GetAllPhaseData(GPXFile,PhaseName)
#         #need cell compatibility check here
#         Layer = {'Name':Phase['General']['Name'],'SameAs':'','Symm':'None'}
#         cx,ct,cs,cia = Phase['General']['AtomPtrs']
#         atoms = Phase['Atoms']
#         Atoms = []
#         for atom in atoms:
#             x,y,z,f = atom[cx:cx+4]
#             u = atom[cia+1]
#             if not u: u = 0.01
#             Atoms.append([atom[ct-1],atom[ct],x,y,z,f,u])
#             if atom[ct] not in data['Layers']['AtInfo']:
#                 data['Layers']['AtInfo'][atom[ct]] = G2elem.GetAtomInfo(atom[ct])
#         Layer['Atoms'] = Atoms
#         data['Layers']['Layers'].append(Layer)
#         Trans = data['Layers']['Transitions']
#         if len(Trans):
#             Trans.append([[0.,0.,0.,0.,'',False] for trans in Trans])
#             for trans in Trans:
#                 trans.append([0.,0.,0.,0.,'',False])
#         else:
#             Trans = [[[1.,0.,0.,0.,'',False],],]
#         data['Layers']['Transitions'] = Trans
#         wx.CallAfter(UpdateLayerData)

#     def LayerSizer(il,Layer):

#         def OnNameChange(event):
#             event.Skip()
#             Layer['Name'] = layerName.GetValue()
#             wx.CallLater(100,UpdateLayerData)

#         def OnAddAtom(event):
#             Layer['Atoms'].append(['Unk','Unk',0.,0.,0.,1.,0.01])
#             wx.CallAfter(UpdateLayerData)

#         def OnSymm(event):
#             Layer['Symm'] = symm.GetValue()

#         def AtomTypeSelect(event):
#             r,c =  event.GetRow(),event.GetCol()
#             if atomGrid.GetColLabelValue(c) == 'Type':
#                 PE = G2elemGUI.PickElement(G2frame)
#                 if PE.ShowModal() == wx.ID_OK:
#                     if PE.Elem != 'None':
#                         atType = PE.Elem.strip()
#                         Layer['Atoms'][r][c] = atType
#                         name = Layer['Atoms'][r][c]
#                         if len(name) in [2,4]:
#                             Layer['Atoms'][r][c-1] = name[:2]+'%d'%(r+1)
#                         else:
#                             Layer['Atoms'][r][c-1] = name[:1]+'%d'%(r+1)
#                         if atType not in data['Layers']['AtInfo']:
#                             data['Layers']['AtInfo'][atType] = G2elem.GetAtomInfo(atType)
#                 PE.Destroy()
#                 wx.CallAfter(UpdateLayerData)
#             else:
#                 event.Skip()

#         def OnDrawLayer(event):
#             drawLayer.SetValue(False)
#             G2plt.PlotLayers(G2frame,Layers,[il,],plotDefaults,firstCall=True)

#         def OnSameAs(event):
#             Layer['SameAs'] = sameas.GetValue()
#             wx.CallLater(100,UpdateLayerData)

#         layerSizer = wx.BoxSizer(wx.VERTICAL)
#         nameSizer = wx.BoxSizer(wx.HORIZONTAL)
#         nameSizer.Add(wx.StaticText(layerData,label=' Layer name: '),0,WACV)
# #            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
#         layerName = wx.TextCtrl(layerData,value=Layer['Name'],style=wx.TE_PROCESS_ENTER)
#         layerName.Bind(wx.EVT_TEXT_ENTER,OnNameChange)
#         layerName.Bind(wx.EVT_KILL_FOCUS,OnNameChange)
#         layerName.Bind(wx.EVT_LEAVE_WINDOW,OnNameChange)
#         nameSizer.Add(layerName,0,WACV)
#         if il:
#             nameSizer.Add(wx.StaticText(layerData,label=' Same as: '),0,WACV)
#             sameas = wx.ComboBox(layerData,value=Layer['SameAs'],choices=['',]+layerNames[:-1],
#                 style=wx.CB_READONLY|wx.CB_DROPDOWN)
#             sameas.Bind(wx.EVT_COMBOBOX, OnSameAs)
#             nameSizer.Add(sameas,0,WACV)
#             if Layer['SameAs']:
#                 indx = layerNames.index(Layer['SameAs'])
#                 if indx < il:    #previously used : same layer
#                     layerSizer.Add(nameSizer)
#                     return layerSizer
#         nameSizer.Add(wx.StaticText(layerData,label=' Layer symmetry: '),0,WACV)
#         symmChoice = ['-1','None']
#         symm = wx.ComboBox(layerData,value=Layer['Symm'],choices=symmChoice,
#             style=wx.CB_READONLY|wx.CB_DROPDOWN)
#         symm.Bind(wx.EVT_COMBOBOX,OnSymm)
#         nameSizer.Add(symm,0,WACV)
#         addAtom = wx.CheckBox(layerData,label=' Add atom? ')
#         addAtom.Bind(wx.EVT_CHECKBOX, OnAddAtom)
#         nameSizer.Add(addAtom,0,WACV)
#         drawLayer = wx.CheckBox(layerData,label=' Draw layer? ')
#         drawLayer.Bind(wx.EVT_CHECKBOX, OnDrawLayer)
#         nameSizer.Add(drawLayer,0,WACV)
#         layerSizer.Add(nameSizer)
#         table = []
#         rowLabels = []
#         for i,atom in enumerate(Layer['Atoms']):
#             table.append(atom)
#             rowLabels.append(str(i))
#         atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=colTypes)
#         atomGrid = G2G.GSGrid(layerData)
#         atomGrid.SetTable(atomTable,True,useFracEdit=False)
# #            atomGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
#         # loop over all cols in table, set cell editor for numerical items
#         for c,t in enumerate(colTypes):
#             if not t.startswith(wg.GRID_VALUE_FLOAT): continue
#             attr = wx.grid.GridCellAttr()
#             attr.IncRef()               #fix from Jim Hester
#             attr.SetEditor(G2G.GridFractionEditor(atomGrid))
#             atomGrid.SetColAttr(c, attr)
#         for row,atom in enumerate(Layer['Atoms']):
#             atomGrid.SetReadOnly(row,1,True)
#         atomGrid.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
#         atomGrid.AutoSizeColumns(True)
#         layerSizer.Add(atomGrid)
#         return layerSizer

#     def TransSizer():

#         def PlotSelect(event):
#             Obj = event.GetEventObject()
#             Yi = Indx[Obj.GetId()]
#             Xi,c =  event.GetRow(),event.GetCol()
#             if Xi >= 0 and c == 5:   #plot column
#                 G2plt.PlotLayers(G2frame,Layers,[Yi,Xi,],plotDefaults,firstCall=True)
#             else:
#                 Psum = 0.
#                 for Xi in range(len(transArray)):
#                     Psum += transArray[Xi][Xi][0]
#                 Psum /= len(transArray)
#                 totalFault.SetLabel(' Total fault density = %.3f'%(1.-Psum))
#                 event.Skip()

#         def OnNormProb(event):
#             for Yi,Yname in enumerate(Names):
#                 Psum = 0.
#                 for Xi,Xname in enumerate(Names):
#                     Psum += transArray[Yi][Xi][0]
#                 if not Psum:
#                     transArray[Yi][0][0] = 1.0
#                     Psum = 1.0
#                 for Xi,Xname in enumerate(Names):
#                     transArray[Yi][Xi][0] /= Psum
#             wx.CallAfter(UpdateLayerData)

#         def OnSymProb(event):
#             if symprob.GetValue():
#                 Nx = len(Names)-1
#                 Layers['SymTrans'] = True
#                 for Yi,Yname in enumerate(Names):
#                     for Xi,Xname in enumerate(Names):
#                         if transArray[Nx-Yi][Nx-Xi][0] != transArray[Yi][Xi][0]:
#                             Layers['SymTrans'] = False
#                             symprob.SetValue(False)
#                             wx.MessageBox('%s-%s not equal %s-%s'%(Yname,Xname,Xname,Yname),
#                                 caption='Probability symmetry error',style=wx.ICON_EXCLAMATION)
#                             break
#             else:
#                 Layers['SymTrans'] = False

#         transSizer = wx.BoxSizer(wx.VERTICAL)
#         transSizer.Add(wx.StaticText(layerData,label=' Layer-Layer transition probabilities: '),0)
#         topSizer = wx.BoxSizer(wx.HORIZONTAL)
#         normprob = wx.CheckBox(layerData,label=' Normalize probabilities?')
#         normprob.Bind(wx.EVT_CHECKBOX,OnNormProb)
#         topSizer.Add(normprob,0,WACV)
#         symprob = wx.CheckBox(layerData,label=' Symmetric probabilities?')
#         symprob.SetValue(Layers.get('SymTrans',False))
#         symprob.Bind(wx.EVT_CHECKBOX,OnSymProb)
#         topSizer.Add(symprob,0,WACV)
#         transSizer.Add(topSizer,0)
#         Names = [layer['Name'] for layer in Layers['Layers']]
#         transArray = Layers['Transitions']
#         layerData.transGrids = []
#         if not Names or not transArray:
#             return transSizer
#         diagSum = 0.
#         for Yi,Yname in enumerate(Names):
#             transSizer.Add(wx.StaticText(layerData,label=' From %s to:'%(Yname)),0)
#             table = []
#             rowLabels = []
#             diagSum += transArray[Yi][Yi][0]
#             for Xi,Xname in enumerate(Names):
#                 table.append(transArray[Yi][Xi])
#                 rowLabels.append(Xname)
#                 if transArray[Yi][Xi][0] > 0.:
#                     Layers['allowedTrans'].append([str(Yi+1),str(Xi+1)])
#             transTable = G2G.Table(table,rowLabels=rowLabels,colLabels=transLabels,types=transTypes)
#             transGrid = G2G.GSGrid(layerData)
#             transGrid.SetTable(transTable,True,useFracEdit=False)
# #                transGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
#             Indx[transGrid.GetId()] = Yi
#             for c,t in enumerate(transTypes):
#                 if not t.startswith(wg.GRID_VALUE_FLOAT): continue
#                 attr = wx.grid.GridCellAttr()
#                 attr.IncRef()               #fix from Jim Hester
#                 attr.SetEditor(G2G.GridFractionEditor(transGrid))
#                 transGrid.SetColAttr(c, attr)
#             transGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, PlotSelect)
#             transGrid.AutoSizeColumns(True)
#             transSizer.Add(transGrid)
#             layerData.transGrids.append(transGrid)
#         if len(transArray):
#             diagSum /= len(transArray)
#             totalFault = wx.StaticText(layerData,
#                 label=' Total fault density = %.3f'%(1.-diagSum))
#             transSizer.Add(totalFault,0)
#         return transSizer

#     def PlotSizer():

#         def OnPlotSeq(event):
#             event.Skip()
#             vals = plotSeq.GetValue().split()
#             try:
#                 vals = [int(val)-1 for val in vals]
#                 if not all([0 <= val < len(Names) for val in vals]):
#                     raise ValueError
#             except ValueError:
#                 plotSeq.SetValue('Error in string '+plotSeq.GetValue())
#                 return
#             G2plt.PlotLayers(G2frame,Layers,vals,plotDefaults,firstCall=True)

#         Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
#         plotSizer = wx.BoxSizer(wx.VERTICAL)
#         Str = ' Using sequence nos. from:'
#         for name in Names:
#             Str += name
#         plotSizer.Add(wx.StaticText(layerData,label=Str[:-1]),0)
#         lineSizer = wx.BoxSizer(wx.HORIZONTAL)
#         lineSizer.Add(wx.StaticText(layerData,label=' Enter sequence of layers to plot:'),0,WACV)
# #            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
#         plotSeq = wx.TextCtrl(layerData,value = '',style=wx.TE_PROCESS_ENTER)
#         plotSeq.Bind(wx.EVT_TEXT_ENTER,OnPlotSeq)
#         plotSeq.Bind(wx.EVT_KILL_FOCUS,OnPlotSeq)
#         lineSizer.Add(plotSeq,0,WACV)
#         plotSizer.Add(lineSizer,0)
#         return plotSizer

#     def StackSizer():

#         stackChoice = ['recursive','explicit',]
#         seqChoice = ['random','list',]

#         def OnStackType(event):
#             newType = stackType.GetValue()
#             if newType == data['Layers']['Stacking'][0]:
#                 return
#             data['Layers']['Stacking'][0] = newType
#             if newType == 'recursive':
#                 data['Layers']['Stacking'][1] = 'infinite'
#             else:  #explicit
#                 data['Layers']['Stacking'][1] = 'random'
#                 data['Layers']['Stacking'][2] = '250'
#             wx.CallAfter(UpdateLayerData)

#         def OnSeqType(event):
#             newType = seqType.GetValue()
#             if newType == data['Layers']['Stacking'][1]:
#                 return
#             data['Layers']['Stacking'][1] = newType
#             if newType == 'random':
#                 data['Layers']['Stacking'][2] = '250'
#             else: #List
#                 data['Layers']['Stacking'][2] = ''
#             wx.CallAfter(UpdateLayerData)

#         def OnNumLayers(event):
#             event.Skip()
#             val = numLayers.GetValue()
#             if val == 'infinite':
#                 data['Layers']['Stacking'][1] = val
#             else:
#                 try:
#                     if 0 < int(val) < 1023:
#                         data['Layers']['Stacking'][1] = val
#                     else:
#                         data['Layers']['Stacking'][1] = 'infinite'
#                 except ValueError:
#                     pass
#             numLayers.SetValue(data['Layers']['Stacking'][1])

#         def OnNumRan(event):
#             event.Skip()
#             val = numRan.GetValue()
#             try:
#                 if 0 > int(val) > 1022:
#                     raise ValueError
#                 else:
#                     data['Layers']['Stacking'][2] = val
#             except ValueError:
#                 val = data['Layers']['Stacking'][2]
#             numRan.SetValue(val)

#         def OnStackList(event):
#             event.Skip()
#             stack = stackList.GetValue()
#             stack = stack.replace('\n',' ').strip().strip('\n')
#             nstar = stack.count('*')
#             if nstar:
#                 try:
#                     newstack = ''
#                     Istar = 0
#                     for star in range(nstar):
#                         Istar = stack.index('*',Istar+1)
#                         iB = stack[:Istar].rfind(' ')
#                         if iB == -1:
#                             mult = int(stack[:Istar])
#                         else:
#                             mult = int(stack[iB:Istar])
#                         pattern = stack[Istar+2:stack.index(')',Istar)]+' '
#                         newstack += mult*pattern
#                     stack = newstack
#                 except ValueError:
#                     stack += ' Error in string'
#             Slist = stack.split()
#             if len(Slist) < 2:
#                 stack = 'Error in sequence - too short!'
#             OKlist = [Slist[i:i+2] in Layers['allowedTrans'] for i in range(len(Slist[:-1]))]
#             if all(OKlist):
#                 data['Layers']['Stacking'][2] = stack
#             else:
#                 stack = 'Improbable sequence or bad string'
#             stackList.SetValue(stack)

#         stackSizer = wx.BoxSizer(wx.VERTICAL)
#         stackSizer.Add(wx.StaticText(layerData,label=' Layer stacking parameters:'),0)
#         if not Layers['Stacking']:
#             Layers['Stacking'] = ['recursive','infinite','']
#         topLine = wx.BoxSizer(wx.HORIZONTAL)
#         topLine.Add(wx.StaticText(layerData,label=' Stacking type: '),0,WACV)
#         stackType = wx.ComboBox(layerData,value=Layers['Stacking'][0],choices=stackChoice,
#             style=wx.CB_READONLY|wx.CB_DROPDOWN)
#         stackType.Bind(wx.EVT_COMBOBOX,OnStackType)
#         topLine.Add(stackType,0,WACV)
#         if Layers['Stacking'][0] == 'recursive':
#             topLine.Add(wx.StaticText(layerData,label=' number of layers (<1022 or "infinite"): '),0,WACV)
# #            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
#             numLayers = wx.TextCtrl(layerData,value=data['Layers']['Stacking'][1],style=wx.TE_PROCESS_ENTER)
#             numLayers.Bind(wx.EVT_TEXT_ENTER,OnNumLayers)
#             numLayers.Bind(wx.EVT_KILL_FOCUS,OnNumLayers)
#             topLine.Add(numLayers,0,WACV)
#             stackSizer.Add(topLine)
#         elif Layers['Stacking'][0] == 'explicit':
#             topLine.Add(wx.StaticText(layerData,label=' layer sequence: '),0,WACV)
#             seqType = wx.ComboBox(layerData,value=data['Layers']['Stacking'][1],choices=seqChoice,
#                 style=wx.CB_READONLY|wx.CB_DROPDOWN)
#             seqType.Bind(wx.EVT_COMBOBOX,OnSeqType)
#             topLine.Add(seqType,0,WACV)
#             if Layers['Stacking'][1] == 'list':
#                 stackSizer.Add(topLine,0)
#                 Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
#                 stackSizer.Add(wx.StaticText(layerData,label=' Explicit layer sequence; enter space delimited list of numbers:'),0)
#                 Str = ' Use sequence nos. from:'
#                 for name in Names:
#                     Str += name
#                 stackSizer.Add(wx.StaticText(layerData,label=Str[:-1]+' Repeat sequences can be used: e.g. 6*(1 2) '),0)
#                 stackSizer.Add(wx.StaticText(layerData,label=' Zero probability sequences not allowed'),0)
#                 stackList = wx.TextCtrl(layerData,value=Layers['Stacking'][2],size=(600,-1),
#                     style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER)
#                 stackList.Bind(wx.EVT_TEXT_ENTER,OnStackList)
#                 stackList.Bind(wx.EVT_KILL_FOCUS,OnStackList)
#                 stackSizer.Add(stackList,0,wx.ALL|wx.EXPAND,8)
#             else:   #random
#                 topLine.Add(wx.StaticText(layerData,label=' Length of random sequence: '),0,WACV)
# #            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
#                 numRan = wx.TextCtrl(layerData,value=Layers['Stacking'][2],style=wx.TE_PROCESS_ENTER)
#                 numRan.Bind(wx.EVT_TEXT_ENTER,OnNumRan)
#                 numRan.Bind(wx.EVT_KILL_FOCUS,OnNumRan)
#                 topLine.Add(numRan,0,WACV)
#                 stackSizer.Add(topLine,0)
#         return stackSizer

#     Layers = data['Layers']
#     layerNames = []
#     Layers['allowedTrans'] = []
#     if len(Layers['Layers']):
#         layerNames = [layer['Name'] for layer in Layers['Layers']]
#     G2frame.GetStatusBar().SetStatusText('',1)
#     layerData = G2frame.layerData
#     try:
#         if layerData.GetSizer():
#             layerData.GetSizer().Clear(True)
#     except:
#         pass
#     mainSizer = wx.BoxSizer(wx.VERTICAL)
#     topSizer = wx.BoxSizer(wx.VERTICAL)
#     bottomSizer = wx.BoxSizer(wx.VERTICAL)
#     headSizer = wx.BoxSizer(wx.HORIZONTAL)
#     headSizer.Add(wx.StaticText(layerData,label=' Global layer description:  '),0,WACV)
#     if 'Sadp' in Layers:
#         sadpPlot = wx.CheckBox(layerData,label=' Plot selected area diffraction?')
#         sadpPlot.Bind(wx.EVT_CHECKBOX,OnSadpPlot)
#         headSizer.Add(sadpPlot,0,WACV)
#     if 'seqResults' in Layers:
#         seqPlot = wx.CheckBox(layerData,label=' Plot sequential result?')
#         seqPlot.Bind(wx.EVT_CHECKBOX,OnSeqPlot)
#         headSizer.Add(seqPlot,0,WACV)
#     topSizer.Add(headSizer)
#     laueSizer = wx.BoxSizer(wx.HORIZONTAL)
#     laueSizer.Add(wx.StaticText(layerData,label=' Diffraction Laue symmetry:'),0,WACV)
#     laue = wx.ComboBox(layerData,value=Layers['Laue'],choices=laueChoice,
#         style=wx.CB_READONLY|wx.CB_DROPDOWN)
#     laue.Bind(wx.EVT_COMBOBOX,OnLaue)
#     laueSizer.Add(laue,0,WACV)
#     if Layers['Laue'] == 'unknown':
#         laueSizer.Add(wx.StaticText(layerData,label=' Diffraction symmetry tolerance: '),0,WACV)
#         toler = G2G.ValidatedTxtCtrl(layerData,Layers,'Toler',nDig=(10,3))
#         laueSizer.Add(toler,0,WACV)
#     topSizer.Add(laueSizer,0)
#     topSizer.Add(wx.StaticText(layerData,label=' Reference unit cell for all layers:'),0)
#     topSizer.Add(CellSizer(),0)
#     topSizer.Add(WidthSizer())
#     topSizer.Add(wx.StaticText(layerData,label=' NB: stacking fault refinement currently not available'),0)
#     G2G.HorizontalLine(topSizer,layerData)
#     titleSizer = wx.BoxSizer(wx.HORIZONTAL)
#     titleSizer.Add(wx.StaticText(layerData,label=' Layer descriptions: '),0,WACV)
#     newLayer = wx.Button(layerData,label='Add new layer')
#     newLayer.Bind(wx.EVT_BUTTON, OnNewLayer)
#     titleSizer.Add(newLayer,0,WACV)
#     importLayer = wx.Button(layerData,label=' Import new layer')
#     importLayer.Bind(wx.EVT_BUTTON, OnImportLayer)
#     titleSizer.Add(importLayer,0,WACV)
#     deleteLast = wx.Button(layerData,label='Delete last layer')
#     deleteLast.Bind(wx.EVT_BUTTON, OnDeleteLast)
#     titleSizer.Add(deleteLast,0,WACV)
#     topSizer.Add(titleSizer,0)
#     for il,layer in enumerate(Layers['Layers']):
#         topSizer.Add(LayerSizer(il,layer))
#     G2G.HorizontalLine(topSizer,layerData)
#     mainSizer.Add(topSizer)
#     bottomSizer.Add(TransSizer())
#     G2G.HorizontalLine(bottomSizer,layerData)
#     bottomSizer.Add(PlotSizer(),0)
#     G2G.HorizontalLine(bottomSizer,layerData)
#     bottomSizer.Add(StackSizer())
#     mainSizer.Add(bottomSizer)
#     SetPhaseWindow(G2frame.layerData,mainSizer,Scroll=Scroll)

# def OnCopyPhase(event):
#     dlg = wx.FileDialog(G2frame, 'Choose GSAS-II project file', G2G.GetImportPath(G2frame),
#         wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN| wx.FD_CHANGE_DIR)
#     try:
#         if dlg.ShowModal() == wx.ID_OK:
#             GPXFile = dlg.GetPath()
#             phaseNames = G2stIO.GetPhaseNames(GPXFile)
#         else:
#             return
#     finally:
#         dlg.Destroy()
#     dlg = wx.SingleChoiceDialog(G2frame,'Phase to use for cell data','Select',phaseNames)
#     if dlg.ShowModal() == wx.ID_OK:
#         sel = dlg.GetSelection()
#         PhaseName = phaseNames[sel]
#     else:
#         return
#     General = G2stIO.GetAllPhaseData(GPXFile,PhaseName)['General']
#     data['Layers']['Cell'] = General['Cell']
#     wx.CallAfter(UpdateLayerData)

# def OnLoadDIFFaX(event):
#     if len(data['Layers']['Layers']):
#         dlg = wx.MessageDialog(G2frame,'Do you really want to replace the Layer data?','Load DIFFaX file',
#             wx.YES_NO | wx.ICON_QUESTION)
#         try:
#             result = dlg.ShowModal()
#             if result == wx.ID_NO:
#                 return
#         finally:
#             dlg.Destroy()
#     dlg = wx.FileDialog(G2frame, 'Choose DIFFaX file name to read', G2G.GetImportPath(G2frame), '',
#         'DIFFaX file (*.*)|*.*',style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
#     try:
#         if dlg.ShowModal() == wx.ID_OK:
#             DIFFaXfile = dlg.GetPath()
#             data['Layers'] = G2IO.ReadDIFFaX(DIFFaXfile)
#     finally:
#         dlg.Destroy()
#     wx.CallAfter(UpdateLayerData)

# def OnSimulate(event):
#     debug = False       #set True to run DIFFax to compare/debug (must be in bin)
#     idebug = 0
#     if debug: idebug = 1
#     wx.MessageBox(' For use of DIFFaX, please cite:\n\n'+
#                       G2G.GetCite('DIFFaX'),
#                       caption='DIFFaX',style=wx.ICON_INFORMATION)
#     ctrls = ''
#     dlg = DIFFaXcontrols(G2frame,ctrls)
#     if dlg.ShowModal() == wx.ID_OK:
#         simCodes = dlg.GetSelection()
#     else:
#         return

#     if 'PWDR' in  simCodes[0]:    #powder pattern
#         data['Layers']['selInst'] = simCodes[1]
#         UseList = [item for item in data['Histograms'] if 'PWDR' in item]
#         if not UseList:
#             wx.MessageBox('No PWDR data for this phase to simulate',caption='Data error',style=wx.ICON_EXCLAMATION)
#             return
#         elif len(UseList) == 1: # don't ask questions when we know the answer!
#             HistName = UseList[0]
#         else:
#             dlg = wx.SingleChoiceDialog(G2frame,'Data to simulate','Select',UseList)
#             if dlg.ShowModal() == wx.ID_OK:
#                 sel = dlg.GetSelection()
#                 HistName = UseList[sel]
#             else:
#                 return
#             dlg.Destroy()
#         G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistName)
#         sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#             G2frame,G2frame.PatternId, 'Sample Parameters'))
#         scale = sample['Scale'][0]
#         background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#             G2frame,G2frame.PatternId, 'Background'))
#         limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#             G2frame,G2frame.PatternId, 'Limits'))[1]
#         inst = G2frame.GPXtree.GetItemPyData(
#             G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
#         if 'T' in inst['Type'][0]:
#             wx.MessageBox("Can't simulate neutron TOF patterns yet",caption='Data error',style=wx.ICON_EXCLAMATION)
#             return
#         profile = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1]
#         G2pwd.CalcStackingPWDR(data['Layers'],scale,background,limits,inst,profile,debug)
#         if debug:
#             ctrls = '0\n%d\n3\n'%(idebug)
#             G2pwd.StackSim(data['Layers'],ctrls,scale,background,limits,inst,profile)
#             test1 = np.copy(profile[3])
#             test1 = np.where(test1,test1,1.0)
#             G2pwd.CalcStackingPWDR(data['Layers'],scale,background,limits,inst,profile,debug)
#             test2 = np.copy(profile[3])
#             rat = test1-test2
#             XY = np.vstack((profile[0],rat))
#             G2plt.PlotXY(G2frame,[XY,],XY2=[],labelX=r'$\mathsf{2\theta}$',
#                 labelY='difference',newPlot=True,Title='DIFFaX vs GSASII',lines=True)
#         G2pwpl.PlotPatterns(G2frame,plotType='PWDR',newPlot=True)
#     else:   #selected area
#         data['Layers']['Sadp'] = {}
#         data['Layers']['Sadp']['Plane'] = simCodes[1]
#         data['Layers']['Sadp']['Lmax'] = simCodes[2]
#         if debug:
#             planeChoice = ['h0l','0kl','hhl','h-hl',]
#             lmaxChoice = [str(i+1) for i in range(6)]
#             ctrls = '0\n%d\n4\n1\n%d\n%d\n16\n1\n1\n0\nend\n'%    \
#                 (idebug,planeChoice.index(simCodes[1])+1,lmaxChoice.index(simCodes[2])+1)
#             G2pwd.StackSim(data['Layers'],ctrls)
#         G2pwd.CalcStackingSADP(data['Layers'],debug)
#     wx.MessageBox('Simulation finished',caption='Stacking fault simulation',style=wx.ICON_EXCLAMATION)
#     wx.CallAfter(UpdateLayerData)

# def OnFitLayers(event):
#     print (' fit stacking fault model TBD')
# #        import scipy.optimize as opt
#     wx.BeginBusyCursor()
#     # see pwd.SetupPDFEval() and pwd.OptimizePDF() for an example minimization
#     wx.EndBusyCursor()
#     wx.CallAfter(UpdateLayerData)
#     G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

# def OnSeqSimulate(event):

#     cellSel = ['cellA','cellB','cellC','cellG']
#     transSel = ['TransP','TransX','TransY','TransZ']
#     ctrls = ''
#     cell = data['Layers']['Cell']
#     data['Layers']['seqResults'] = []
#     data['Layers']['seqCodes'] = []
#     Parms = G2pwd.GetStackParms(data['Layers'])
#     dlg = DIFFaXcontrols(G2frame,ctrls,Parms)
#     if dlg.ShowModal() == wx.ID_OK:
#         simCodes = dlg.GetSelection()
#     else:
#         return
#     UseList = []
#     for item in data['Histograms']:
#         if 'PWDR' in item:
#             UseList.append(item)
#     if not UseList:
#         wx.MessageBox('No PWDR data for this phase to simulate',caption='Data error',style=wx.ICON_EXCLAMATION)
#         return
#     dlg = wx.SingleChoiceDialog(G2frame,'Data to simulate','Select',UseList)
#     if dlg.ShowModal() == wx.ID_OK:
#         sel = dlg.GetSelection()
#         HistName = UseList[sel]
#     else:
#         return
#     dlg.Destroy()
#     G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistName)
#     sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#         G2frame,G2frame.PatternId, 'Sample Parameters'))
#     scale = sample['Scale'][0]
#     background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#         G2frame,G2frame.PatternId, 'Background'))
#     limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
#         G2frame,G2frame.PatternId, 'Limits'))[1]
#     inst = G2frame.GPXtree.GetItemPyData(
#         G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
#     if 'T' in inst['Type'][0]:
#         wx.MessageBox("Can't simulate neutron TOF patterns yet",caption='Data error',style=wx.ICON_EXCLAMATION)
#         return
#     profile = np.copy(G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1])
#     resultXY2 = []
#     resultXY = [np.vstack((profile[0],profile[1])),]    #observed data
#     data['Layers']['selInst'] = simCodes[1]
#     data['Layers']['seqCodes'] = simCodes[2:]
#     Layers = copy.deepcopy(data['Layers'])
#     pName = simCodes[2]
#     BegFin = simCodes[3]
#     nSteps = simCodes[4]
#     laue = Layers['Laue']
#     vals = np.linspace(BegFin[0],BegFin[1],nSteps+1,True)
#     simNames = []
#     for val in vals:
#         print (' Stacking simulation step for %s = %.5f'%(pName,val))
#         simNames.append('%.3f'%(val))
#         if 'cell' in pName:
#             cellId = cellSel.index(pName)
#             cell = Layers['Cell']
#             cell[cellId+1] = val
#             if laue in ['-3','-3m','6/m','6/mmm','4/m','4/mmm']:
#                 cell[2] = cell[1]
#             cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
#             Layers['Cell'] = cell
#         elif 'Trans' in pName:
#             names = pName.split(';')
#             transId = transSel.index(names[0])
#             iY = int(names[1])
#             iX = int(names[2])
#             Trans = Layers['Transitions'][iY]
#             Nx = len(Trans)-1
#             if not transId:     #i.e. probability
#                 osum = 1.-Trans[iX][0]
#                 nsum = 1.-val
#                 for i in range(Nx+1):
#                     if i != iX:
#                         Trans[i][0] *= (nsum/osum)
#                 Trans[iX][0] = val
#                 if Layers.get('SymTrans',False):
#                     Layers['Transitions'][Nx-iX][Nx-iY][0] = val
#                     for i in range(Nx+1):
#                         Layers['Transitions'][Nx-iY][Nx-i][0] = Layers['Transitions'][iY][i][0]
#                 print (' Transition matrix:')
#                 for trans in Layers['Transitions']:
#                     line = str([' %.3f'%(item[0]) for item in trans])
#                     print (line[1:-2].replace("'",''))
#             else:
#                 Trans[iX][transId] = val
#         G2pwd.CalcStackingPWDR(Layers,scale,background,limits,inst,profile,False)
#         resultXY2.append([np.vstack((profile[0],profile[3])),][0])
#     data['Layers']['seqResults'] = [resultXY,resultXY2,simNames]
#     wx.MessageBox('Sequential simulation finished',caption='Stacking fault simulation',style=wx.ICON_EXCLAMATION)
#     wx.CallAfter(UpdateLayerData)
