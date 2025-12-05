# -*- coding: utf-8 -*-
#GSASII - phase data display routines
'''
Routines for Phase dataframes follow. Only a few Update routines are here
all others are in GSASIIphsGUI.py
'''
import os
import wx
import wx.grid as wg
import matplotlib as mpl
#import math
import copy

import numpy as np
import numpy.linalg as nl
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIElem as G2elem
from . import GSASIIElemGUI as G2elemGUI
from . import GSASIIplot as G2plt
from . import GSASIIdataGUI as G2gd
from . import GSASIIstrIO as G2stIO
from . import GSASIImath as G2mth
from . import GSASIIpwd as G2pwd
from . import GSASIIctrlGUI as G2G
from . import GSASIIphsGUI as G2phsG

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
#### UpdateDynsomia GUI
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
    treeId = G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists')
    if treeId:
        reflSets = G2frame.GPXtree.GetItemPyData(treeId)    
        reflData = reflSets[generalData['Name']]['RefList']
    else:
        wx.MessageBox('You must have PWDR reflections before running Dysnomia','Dysnomia Error',
            style=wx.ICON_ERROR)
        return
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

####  UpdateDeformation form factor routines ################################################################

def UpdateDeformation(G2frame,data,AtdId):
    
    def OnRadFxn(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        deformationData[-dId]['Radial'] = radFxn.GetStringSelection()
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
        
    def OnSSchoice(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        deformationData[-dId]['LocSS'] = SSchoice.GetStringSelection()
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
    
    def MakeUVmat(defData,U,V):
        MX = U/nl.norm(U)
        if 'A' in defData['MUV']:
            MY = V/nl.norm(V)
            MZ = np.cross(MX,MY)
            MZ /= nl.norm(MZ)
            MY = np.cross(MZ,MX)
            MY /= nl.norm(MY)
        else:
            MZ = V/nl.norm(V)
            MY = np.cross(MZ,MX)
            MY /= nl.norm(MY)
            MZ = np.cross(MX,MY)
            MZ /= nl.norm(MZ)
        return np.array([MX,MY,MZ]).T
    
    def OnDeformRef(event):
        Obj = event.GetEventObject()
        dId,oId,dkey = Indx[Obj.GetId()]
        deformationData[dId][oId][1][dkey][1] = not deformationData[dId][oId][1][dkey][1]

    def OnPlotAtm(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        atom = atomData[AtLookUp[dId]]
        neigh = G2mth.FindAllNeighbors(data,atom[ct-1],AtNames)
        deform = deformationData[dId]
        UVmat = deformationData[-dId]['UVmat']
        radial = deformationData[-dId]['Radial']
        G2plt.PlotDeform(G2frame,generalData,atom[ct-1],atom[ct],deform,UVmat,radial,neigh)            

    def OnDelAtm(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        del deformationData[dId]
        wx.CallAfter(UpdateDeformation,G2frame,data,None)

    def OnMatSel(event):
        "Cartesian axes: A: X'=U, Y'=(UxV)xU & Z'=UxV,B: X'=U, Y'=UxV & Z'=Ux(UxV)"
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        deformationData[-dId]['MUV'] = Obj.GetValue()
        U = UVvec[dId][UVchoice[dId].index(deformationData[-dId]['U'])]
        V = UVvec[dId][UVchoice[dId].index(deformationData[-dId]['V'])]
        UVmat = MakeUVmat(deformationData[-dId],U,V)
        data['Deformations'][-dId]['UVmat'] = UVmat
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)

    def OnUvec(event):
        "Cartesian axes: A: X'=U, Y'=(UxV)xU & Z'=UxV,B: X'=U, Y'=UxV & Z'=Ux(UxV)"
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        if Obj.GetValue() == deformationData[-dId]['V']:
            Obj.SetValue(deformationData[-dId]['U'])
        else:
            U = UVvec[dId][Obj.GetSelection()]
            V = UVvec[dId][UVchoice[dId].index(deformationData[-dId]['V'])]
            UVmat = MakeUVmat(deformationData[-dId],U,V)
            if np.any(np.isnan(UVmat)):
                Obj.SetValue(deformationData[-dId]['U'])
                G2G.G2MessageBox(G2frame,'ERROR: Z: U-vector zero or parallel to V','Invalid vector choice')
                return
            if nl.det(UVmat) < 0.:  #ensure right hand
                UVmat *= -1.
            deformationData[-dId]['U'] =  Obj.GetValue()
            data['Deformations'][-dId]['UVmat'] = UVmat
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
        
    def OnVvec(event):
        "Cartesian axes: A: X'=U, Y'=(UxV)xU & Z'=UxV,B: X'=U, Y'=UxV & Z'=Ux(UxV)"
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        if Obj.GetValue() == deformationData[-dId]['U']:
            Obj.SetValue(deformationData[-dId]['V'])
        else:
            U = UVvec[dId][UVchoice[dId].index(deformationData[-dId]['U'])]
            V = UVvec[dId][Obj.GetSelection()]
            UVmat = MakeUVmat(deformationData[-dId],U,V)
            if np.any(np.isnan(UVmat)):
                Obj.SetValue(deformationData[-dId]['V'])
                G2G.G2MessageBox(G2frame,'ERROR: V-vector zero or parallel to U','Invalid vector choice')
                return
            if nl.det(UVmat) < 0.:  #ensure right hand
                UVmat *= -1.
            deformationData[-dId]['V'] =  Obj.GetValue()
            data['Deformations'][-dId]['UVmat'] = UVmat
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
            
    def OnAtSel(event):
        dId = atomList[atSel.GetValue()]
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
        
    def Kappa(deformation,orbSizer,dId,orb,kname,Indx):
        orbSizer.Add(G2G.ValidatedTxtCtrl(deformation,orb[1][kname],0,nDig=(8,3),xmin=0.5,xmax=1.5))
        Tcheck = wx.CheckBox(deformation,-1,'Refine?')
        Tcheck.SetValue(orb[1][kname][1])
        Tcheck.Bind(wx.EVT_CHECKBOX,OnDeformRef)
        Indx[Tcheck.GetId()] = [dId,iorb,kname]
        orbSizer.Add(Tcheck)
        
    def NeSizer(deformation,orbSizer,dId,orb,Indx):
        orbSizer.Add(G2G.ValidatedTxtCtrl(deformation,orb[1]['Ne'],0,nDig=(8,3),xmin=0.,xmax=10.))
        Tcheck = wx.CheckBox(deformation,-1,'Refine?')
        Tcheck.SetValue(orb[1]['Ne'][1])
        Tcheck.Bind(wx.EVT_CHECKBOX,OnDeformRef)
        Indx[Tcheck.GetId()] = [dId,iorb,'Ne']
        orbSizer.Add(Tcheck)
        
    def Dsizer(deformation,orbSizer,Names,dId,orb,Indx):
        name = Names.get(item,'') #Names only go to order = 3
        orbSizer.Add(wx.StaticText(deformation,label=item+name+':'))
        orbSizer.Add(G2G.ValidatedTxtCtrl(deformation,orb[1][item],0,nDig=(8,5),xmin=-1.5,xmax=1.5))
        Tcheck = wx.CheckBox(deformation,-1,'Refine?')
        Tcheck.SetValue(orb[1][item][1])
        Tcheck.Bind(wx.EVT_CHECKBOX,OnDeformRef)
        Indx[Tcheck.GetId()] = [dId,iorb,item]
        orbSizer.Add(Tcheck)
        
    def OnNewHarm(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        atom = atomData[AtLookUp[dId]]
        sytsym = atom[cs].strip()
        rbSym = deformationData[-dId]['LocSS']
        for harm in data['Deformations'][dId]:
            if 'Sl' in harm[0]:
                Harm = harm
        Order = 1
        Hkeys = list(Harm[1].keys())
        orders = [int(item[2]) for item in Hkeys if 'D' in item]
        if len(orders):
            Order = max(orders)+1
        cofNames = []
        notFound = True
        while notFound and Order < 6:
            cofNames,cofSgns = G2lat.GenRBCoeff(sytsym,rbSym,Order)
            cofNames = [name.replace('C','D') for name in cofNames]
            for name in cofNames:
                if name not in Hkeys:   #new names found
                    notFound = False
                    Harm[1].update({name:[0.0,False]})
                    # if '0' not in name:
                    #     negname = name.replace(',',',-')
                    #     Harm[1].update({negname:[0.0,False]})
            Order += 1
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
        
    def OnDelHarm(event):
        Obj = event.GetEventObject()
        dId = Indx[Obj.GetId()]
        minL = 1  #always an "Ne"
        for harm in data['Deformations'][dId]:
            if 'Sl' in harm[0]:
                Harm = harm
                minL = 2    # & 'kappa'
        Hkeys = list(Harm[1].keys())
        if len(Hkeys) > minL:
            maxord = max([int(item[2]) for item in Hkeys if 'D' in item])
            for item in Hkeys:
                if 'D' in item and int(item[2]) == maxord:
                    del Harm[1][item]
            if len(Harm[1]) == 3:
                del Harm[1]["kappa'"]
        wx.CallAfter(UpdateDeformation,G2frame,data,dId)
        
    def OnShowDef(event):
        dId = Indx[event.GetEventObject().GetId()]
        deformationData[-dId]['showDef'] = not deformationData[-dId]['showDef']
        G2plt.PlotStructure(G2frame,data)
    
    def OnAtCol(event):
        dId = Indx[event.GetEventObject().GetId()]
        deformationData[-dId]['atColor'] = not deformationData[-dId]['atColor']
        G2plt.PlotStructure(G2frame,data)
        
    # UpdateDeformation executable code starts here
    deformation = G2frame.deformation
    alpha = ['A','B','C','D','E','F','G','H',]
    generalData = data['General']
    cx,ct,cs,cia = generalData['AtomPtrs']
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    atomData = data['Atoms']
    AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
    AtNames = [atom[ct-1] for atom in atomData]
    deformationData = data['Deformations']
    dId = AtdId
    if deformation.GetSizer():
        deformation.GetSizer().Clear(True)
    atomList = {}
    for item in deformationData:
        if item in AtLookUp:
            atom = atomData[AtLookUp[item]]
            atomList.update({atom[ct-1]:item})
    AtChoice = ' '
    if dId is not None:
        AtChoice = atomData[AtLookUp[dId]][ct-1]
    elif len(atomList):
        AtChoice = list(atomList.keys())[0]
        dId = atomList[AtChoice]
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl= f"Atomic deformation data for {data['General']['Name']!r}"[:60]
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)

    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    if dId is None:
        topSizer.Add(wx.StaticText(deformation,
            label='No atoms in deformation list. Do add atom first (neutral atoms only)'),0,WACV)
    else:
        topSizer.Add(wx.StaticText(deformation,label=' Select an atom '),0,WACV)
        atSel = wx.ComboBox(deformation,value=AtChoice,choices=list(atomList.keys()),style=wx.CB_READONLY|wx.CB_DROPDOWN)
        atSel.Bind(wx.EVT_COMBOBOX,OnAtSel)
        topSizer.Add(atSel,0,WACV)
    mainSizer.Add(topSizer,0,wx.EXPAND)
    if dId is not None:
        Indx = {}
        UVchoice = {}
        UVvec = {}
        #patch
        if 'UVmat' not in deformationData[-dId] or 'MUV' not in deformationData[-dId]:
            deformationData[-dId] = {'U':'X','V':'Y','UVmat':np.eye(3),'MUV':"A: X'=U, Y'=(UxV)xU & Z'=UxV"}
        #end patch
        atom = atomData[AtLookUp[dId]]
        neigh = G2mth.FindAllNeighbors(data,atom[ct-1],AtNames)[0]
        neigh = G2mth.sortArray(neigh,2)    #sort by dist
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(deformation,label=' For atom %s, site sym %s:'%(atom[ct-1],atom[cs])),0,WACV)
        plotAtm = wx.Button(deformation,label='Plot')
        plotAtm.Bind(wx.EVT_BUTTON,OnPlotAtm)
        Indx[plotAtm.GetId()] = dId
        lineSizer.Add(plotAtm,0,WACV)
        deformationData[-dId]['showDef'] = deformationData[-dId].get('showDef',False)
        deformationData[-dId]['atColor'] = deformationData[-dId].get('atColor',True)
        showDef = wx.CheckBox(deformation,label='show def.?')
        showDef.SetValue(deformationData[-dId]['showDef'])
        Indx[showDef.GetId()] = dId
        showDef.Bind(wx.EVT_CHECKBOX,OnShowDef)
        lineSizer.Add(showDef,0,WACV)
        atCol = wx.CheckBox(deformation,label='use atom colors?')
        atCol.SetValue(deformationData[-dId]['atColor'])
        Indx[atCol.GetId()] = dId
        atCol.Bind(wx.EVT_CHECKBOX,OnAtCol)
        lineSizer.Add(atCol,0,WACV)
        delAtm = wx.Button(deformation,label='Delete')
        delAtm.Bind(wx.EVT_BUTTON,OnDelAtm)
        Indx[delAtm.GetId()] = dId
        lineSizer.Add(delAtm,0,WACV)
        mainSizer.Add(lineSizer)
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        names = []
        if not len(neigh):
            lineSizer.Add(wx.StaticText(deformation,label=' No neighbors found; Do Set bond parms to expand search'),0,WACV)
        elif len(neigh) < 9:
            names = ['%s=%s'%(alpha[i],item[0].replace(' ','')) for i,item in enumerate(neigh)]
            lineSizer.Add(wx.StaticText(deformation,label=' Neighbors: '+str(names)),0,WACV)
        else:
            names = 'Too many neighbors - change atom radii to fix'
        UVchoice[dId] = ['X','Y','Z','X+Y','X+Y+Z',]
        UVvec[dId] = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,1.,0.]/sqt2,[1.,1.,1.]/sqt3,]
        NUVvec,NUVchoice = G2lat.SetUVvec(neigh)
        UVchoice[dId] += NUVchoice
        UVvec[dId] += NUVvec
        lineSizer.Add(wx.StaticText(deformation,label=' Local site sym:'),0,WACV)
        SSchoices = G2spc.GetSytSymChoice(atom[cs])
        deformationData[-dId]['LocSS'] = deformationData[-dId].get('LocSS',atom[cs][:])
        SSchoice = wx.ComboBox(deformation,value=deformationData[-dId]['LocSS'],
            choices=SSchoices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[SSchoice.GetId()] = dId
        SSchoice.Bind(wx.EVT_COMBOBOX,OnSSchoice)
        lineSizer.Add(SSchoice,0,WACV)
        mainSizer.Add(lineSizer)
        mainSizer.Add(wx.StaticText(deformation,
            label=" NB: Local site sym always has unique axis || Z' and second axis || X'; choose U && V carefully"))
        matSizer = wx.BoxSizer(wx.HORIZONTAL)
        Mchoice = ["A: X'=U, Y'=(UxV)xU & Z'=UxV","B: X'=U, Y'=UxV & Z'=Ux(UxV)"]
        matSizer.Add(wx.StaticText(deformation,label=' Orbital Cartesian axes:'),0,WACV)
        matSel = wx.ComboBox(deformation,choices=Mchoice,value=deformationData[-dId]['MUV'],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        matSel.Bind(wx.EVT_COMBOBOX,OnMatSel)
        Indx[matSel.GetId()] = dId
        matSizer.Add(matSel,0,WACV)        
        deformationData[-dId]['Radial'] = deformationData[-dId].get('Radial','Slater')
        topSizer.Add(wx.StaticText(deformation,label=' Select radial fxn: '),0,WACV)
        fxchoice = deformationData[-dId].get('fxchoice',['Slater',])
        radFxn = wx.ComboBox(deformation,value=deformationData[-dId]['Radial'],
            choices=fxchoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Indx[radFxn.GetId()] = dId
        radFxn.Bind(wx.EVT_COMBOBOX,OnRadFxn)
        topSizer.Add(radFxn,0,WACV)
        mainSizer.Add(matSizer)
        oriSizer = wx.BoxSizer(wx.HORIZONTAL)
        oriSizer.Add(wx.StaticText(deformation,label=' Select orbital U vector: '),0,WACV)
        Uvec = wx.ComboBox(deformation,value=deformationData[-dId]['U'],choices=UVchoice[dId],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Uvec.Bind(wx.EVT_COMBOBOX,OnUvec)
        Indx[Uvec.GetId()] = dId
        oriSizer.Add(Uvec,0,WACV)
        oriSizer.Add(wx.StaticText(deformation,label=' Select orbital V vector: '),0,WACV)
        Vvec = wx.ComboBox(deformation,value=deformationData[-dId]['V'],choices=UVchoice[dId],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Vvec.Bind(wx.EVT_COMBOBOX,OnVvec)
        Indx[Vvec.GetId()] = dId
        oriSizer.Add(Vvec,0,WACV)
        if 'Slater' in data['Deformations'][-dId]['Radial']:
            newHarm = wx.Button(deformation,label='Add harmonic')
            newHarm.Bind(wx.EVT_BUTTON,OnNewHarm)
            Indx[newHarm.GetId()] = dId
            oriSizer.Add(newHarm,0,WACV)
            delHarm = wx.Button(deformation,label='Delete highest harmonic')
            delHarm.Bind(wx.EVT_BUTTON,OnDelHarm)
            Indx[delHarm.GetId()] = dId
            oriSizer.Add(delHarm,0,WACV)
        mainSizer.Add(oriSizer)
        G2G.HorizontalLine(mainSizer,deformation)
        Names = {'D(1,-1)':'py','D(1,0)':'pz','D(1,1)':'px',
                 'D(2,-2)':'dxy','D(2,-1)':'dyz','D(2,0)':'dz2','D(2,1)':'dxz','D(2,2)':'dx2-y2',
                 'D(3,-3)':'fy(3x2-y2)','D(3,-2)':'fxyz','D(3,-1)':'fyz2','D(3,0)':'fz3',
                 'D(3,1)':'fxz2','D(3,2)':'fz(x2-y2)','D(3,3)':'fx(x2-3y2)',}
        mainSizer.Add(wx.StaticText(deformation,label=' Deformation parameters:'))
        orbSizer = wx.FlexGridSizer(0,9,2,2)
        for iorb,orb in enumerate(deformationData[dId]):
            # if deformationData[-dId]['Radial'] == 'Bessel' and 'Sl ' not in orb[0]:
            #     if '<j0>' in orb[0]:
            #         orbSizer.Add(wx.StaticText(deformation,label=orb[0]+' Ne:'))
            #         NeSizer(deformation,orbSizer,dId,orb,Indx)
            #         if 'kappa' in orb[1]:
            #             orbSizer.Add(wx.StaticText(deformation,label=' kappa:'))
            #             Kappa(deformation,orbSizer,dId,orb,'kappa',Indx)
            #         for i in range(3): orbSizer.Add((5,5),0)
            #         continue
            #     if 'kappa' in orb[1]:
            #         for i in range(3): orbSizer.Add((5,5),0)
            #         orbSizer.Add(wx.StaticText(deformation,label=orb[0]+" kappa':"))
            #         Kappa(deformation,orbSizer,dId,orb,"kappa'",Indx)
            #     if 'kappa' not in orb[1]:
            #         orbSizer.Add(wx.StaticText(deformation,label=orb[0]+':'))
            #         for i in range(2): orbSizer.Add((5,5),0)
            #     nItem = 0
            #     for item in orb[1]:
            #         if 'D' in item:                            
            #             nItem += 1
            #             Dsizer(deformation,orbSizer,Names,dId,orb,Indx)
            #             if nItem in [2,4,6,8,10]:
            #                 for i in range(3): orbSizer.Add((5,5),0)
            #     for i in range(3): orbSizer.Add((5,5),0)
            # elif deformationData[-dId]['Radial'] == 'Slater' and 'Sl ' in orb[0]: 
                orbSizer.Add(wx.StaticText(deformation,label=orb[0]+' Ne:'))
                NeSizer(deformation,orbSizer,dId,orb,Indx)
                orbSizer.Add(wx.StaticText(deformation,label=' kappa:'))
                Kappa(deformation,orbSizer,dId,orb,'kappa',Indx)
                if len(orb[1]) > 2:
                    orb[1]["kappa'"] = orb[1].get("kappa'",[1.0,False])
                    orbSizer.Add(wx.StaticText(deformation,label=" kappa':"))
                    Kappa(deformation,orbSizer,dId,orb,"kappa'",Indx)
                iD = 1
                for item in orb[1]:
                    if 'D' in item:
                        if iD < int(item[2]):
                            iD = int(item[2])
                            nItems = orbSizer.GetItemCount()%9
                            if nItems:
                                nB = 9-nItems
                                for i in range(nB): orbSizer.Add((5,5),0)
                        Dsizer(deformation,orbSizer,Names,dId,orb,Indx)
        mainSizer.Add(orbSizer)    

    G2phsG.SetPhaseWindow(deformation,mainSizer)


#### UpdateISODISTORT ###############################################################################

def UpdateISODISTORT(G2frame,data,Scroll=0):
    ''' Setup ISODISTORT and present the results. Allow selection of a distortion model for PDFfit or
    GSAS-II structure refinement as a cif file produced by ISODISTORT. Allows manipulation of distortion
    mode displacements selection their refinement for this new phase.
    '''

    def displaySetup():

        def OnParentCif(event):
            dlg = wx.FileDialog(G2frame.ISODIST, 'Select parent cif file',G2frame.LastGPXdir,
                style=wx.FD_OPEN ,wildcard='cif file(*.cif)|*.cif')
            if dlg.ShowModal() == wx.ID_OK:
                fName = dlg.GetFilename()
                fDir = dlg.GetDirectory()
                ISOdata['ParentCIF'] = os.path.join(fDir,fName)
                dlg.Destroy()
            else:
                dlg.Destroy()
            UpdateISODISTORT(G2frame,data)

        def OnUsePhase(event):
            ISOdata['ParentCIF'] = 'Use this phase'
            UpdateISODISTORT(G2frame,data)

        def OnMethodSel(event):
            method = methodSel.GetSelection()+1
            if method in [1,4]:
                ISOdata['ISOmethod'] = method
            UpdateISODISTORT(G2frame,data)

        def OnChildCif(event):
            dlg = wx.FileDialog(G2frame.ISODIST, 'Select child cif file',G2frame.LastGPXdir,
                style=wx.FD_OPEN ,wildcard='cif file(*.cif)|*.cif')
            if dlg.ShowModal() == wx.ID_OK:
                fName = dlg.GetFilename()
                fDir = dlg.GetDirectory()
                ISOdata['ChildCIF'] = os.path.join(fDir,fName)
                dlg.Destroy()
            else:
                dlg.Destroy()
            UpdateISODISTORT(G2frame,data)

        def OnUsePhase2(event):
            ISOdata['ChildCIF'] = 'Use this phase'
            UpdateISODISTORT(G2frame,data)

        topSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer.Add(wx.StaticText(G2frame.ISODIST,label=' ISODISTORT setup controls:'))
        parentSizer = wx.BoxSizer(wx.HORIZONTAL)
        parentSizer.Add(wx.StaticText(G2frame.ISODIST,label=' Parent cif file:'),0,WACV)
        parentCif = wx.Button(G2frame.ISODIST,label=ISOdata['ParentCIF'],size=(300,24))
        parentCif.Bind(wx.EVT_BUTTON,OnParentCif)
        parentSizer.Add(parentCif,0,WACV)
        if 'Use this phase' not in ISOdata['ChildCIF'] and 'Use this phase' not in ISOdata['ParentCIF']:
            usePhase = wx.Button(G2frame.ISODIST,label=' Use this phase? ')
            usePhase.Bind(wx.EVT_BUTTON,OnUsePhase)
            parentSizer.Add(usePhase,0,WACV)
        topSizer.Add(parentSizer)
        #patch
        if 'ISOmethod' not in ISOdata:
            ISOdata['ISOmethod'] = 1
        #end patch
        choice = ['Method 1: Search over all special k points - yields only single Irrep models',
            'Method 2: not implemented in GSAS-II',
            'Method 3: not implemented in GSAS-II',
            'Method 4: Mode decomposition of known child structure']
        methodSel = wx.RadioBox(G2frame.ISODIST,label='Select ISODISTORT method:',choices=choice,
            majorDimension=1,style=wx.RA_SPECIFY_COLS)
        methodSel.SetSelection(ISOdata['ISOmethod']-1)
        methodSel.Bind(wx.EVT_RADIOBOX,OnMethodSel)
        topSizer.Add(methodSel)
        if ISOdata['ISOmethod'] == 4:
            childSizer = wx.BoxSizer(wx.HORIZONTAL)
            childSizer.Add(wx.StaticText(G2frame.ISODIST,label=' Child cif file:'),0,WACV)
            childCif = wx.Button(G2frame.ISODIST,label=ISOdata['ChildCIF'],size=(300,24))
            childCif.Bind(wx.EVT_BUTTON,OnChildCif)
            childSizer.Add(childCif,0,WACV)
            if 'Use this phase' not in ISOdata['ChildCIF'] and 'Use this phase' not in ISOdata['ParentCIF']:
                usePhase2 = wx.Button(G2frame.ISODIST,label=' Use this phase? ')
                usePhase2.Bind(wx.EVT_BUTTON,OnUsePhase2)
                childSizer.Add(usePhase2,0,WACV)
            topSizer.Add(childSizer)

        return topSizer

    def displaySubset():

        def OnLaue(event):
            Obj = event.GetEventObject()
            name = Indx[Obj.GetId()]
            ISOdata['SGselect'][name[:4]] = not ISOdata['SGselect'][name[:4]]
            ISOdata['selection'] = None
            UpdateISODISTORT(G2frame,data)

        def OnAllBtn(event):
            for item in ISOdata['SGselect']:
                ISOdata['SGselect'][item] = not ISOdata['SGselect'][item]
            ISOdata['selection'] = None
            UpdateISODISTORT(G2frame,data)

        topSizer = wx.BoxSizer(wx.VERTICAL)
        G2G.HorizontalLine(topSizer,G2frame.ISODIST)
        topSizer.Add(wx.StaticText(G2frame.ISODIST,label='ISODISTORT Method 1 distortion search results:'))
        topSizer.Add(wx.StaticText(G2frame.ISODIST,label=' Subset selection if desired:'))
        laueName = ['Cubic','Hexagonal','Trigonal','Tetragonal','Orthorhombic','Monoclinic','Triclinic']
        littleSizer = wx.FlexGridSizer(0,8,5,5)
        Indx = {}
        for name in laueName:
            laueCk = wx.CheckBox(G2frame.ISODIST,label=name)
            Indx[laueCk.GetId()] = name
            laueCk.SetValue(ISOdata['SGselect'][name[:4]])
            laueCk.Bind(wx.EVT_CHECKBOX,OnLaue)
            littleSizer.Add(laueCk,0,WACV)
        allBtn = wx.Button(G2frame.ISODIST,label='Toggle all')
        allBtn.Bind(wx.EVT_BUTTON,OnAllBtn)
        littleSizer.Add(allBtn)
        topSizer.Add(littleSizer)
        return topSizer

    def displayRadio():

        def CheckItem(item):
            SGnum = int(item.split()[1].split('*')[0])
            for SGtype in ISOdata['SGselect']:
                if ISOdata['SGselect'][SGtype] and SGnum in SGrange[SGtype]:
                    return True
            return False

        def OnSelect(event):
           r,c = event.GetRow(),event.GetCol()
           if c == 0:
               ISOdata['selection'] = [r,isoTable.GetValue(r,1)]
               for row in range(isoGrid.GetNumberRows()):
                   isoTable.SetValue(row,c,False)
               isoTable.SetValue(r,c,True)
               isoGrid.ForceRefresh()

        SGrange = {'Cubi':np.arange(195,231),'Hexa':np.arange(168,195),'Trig':np.arange(143,168),'Tetr':np.arange(75,143),
                   'Orth':np.arange(16,75),'Mono':np.arange(3,16),'Tric':np.arange(1,3)}
        bottomSizer = wx.BoxSizer(wx.VERTICAL)
        colLabels = ['select',' ISODISTORT order parameter direction description']
        colTypes = [wg.GRID_VALUE_BOOL,wg.GRID_VALUE_STRING,]

        Radio = ISOdata['radio']
        rowLabels = []
        table = []
        for i,item in enumerate(Radio):
            if CheckItem(Radio[item]):
                if ISOdata['selection'] and ISOdata['selection'][0] == i:
                    table.append([True,Radio[item]])
                else:
                    table.append([False,Radio[item]])
                rowLabels.append(str(i))
        isoTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=colTypes)
        isoGrid = G2G.GSGrid(G2frame.ISODIST)
        isoGrid.SetTable(isoTable,True,useFracEdit=False)
        isoGrid.AutoSizeColumns(True)
        isoGrid.SetColLabelAlignment(wx.ALIGN_LEFT,wx.ALIGN_CENTRE)
        bottomSizer.Add(isoGrid)
        attr = wg.GridCellAttr()
        attr.SetReadOnly(True)
        attr.SetBackgroundColour(VERY_LIGHT_GREY)
        isoGrid.SetColAttr(1,attr)
        isoGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, OnSelect)
        return bottomSizer

    def displayModes():

        def OnDispl(event):
            '''Respond to movement of distortion mode slider'''
            Obj = event.GetEventObject()
            idsp,dispVal = Indx[Obj.GetId()]
            modeDisp[idsp] = Obj.GetValue()/1000.
            dispVal.SetValue(modeDisp[idsp])
            err = G2mth.ApplyModeDisp(data)
            if err:
                G2G.G2MessageBox(G2frame,'Do Draw atoms first')
            G2phsG.FindBondsDraw(data)
            G2plt.PlotStructure(G2frame,data)

        def OnDispVal(invalid,value,tc):
            '''Respond to entry of a value into a distortion mode entry widget'''
            idsp,displ = Indx[tc.GetId()]
            displ.SetValue(int(value*1000))
            err = G2mth.ApplyModeDisp(data)
            if err:
                G2G.G2MessageBox(G2frame,'Do Draw atoms first')
            G2phsG.FindBondsDraw(data)
            G2plt.PlotStructure(G2frame,data)

        def OnRefDispl(event):
            Obj = event.GetEventObject()
            idsp,item = Indx[Obj.GetId()]
            item[-2] = not item[-2]

        def OnReset(event):
            '''Reset all distortion mode values to initial values'''
            ISOdata['modeDispl'] = copy.deepcopy(ISOdata['ISOmodeDispl'])
            err = G2mth.ApplyModeDisp(data)
            if err:
                G2G.G2MessageBox(G2frame,'Do Draw atoms first')
            G2phsG.FindBondsDraw(data)
            G2plt.PlotStructure(G2frame,data)
            UpdateISODISTORT(G2frame,data)

        def OnSetZero(event):
            '''Reset all distortion mode values to 0'''
            ISOdata['modeDispl'] = [0.0 for i in ISOdata['ISOmodeDispl']]
            err = G2mth.ApplyModeDisp(data)
            if err:
                G2G.G2MessageBox(G2frame,'Do Draw atoms first')
            G2phsG.FindBondsDraw(data)
            G2plt.PlotStructure(G2frame,data)
            UpdateISODISTORT(G2frame,data)

        def OnSaveModes(event):
            '''Set saved distortion mode values to displayed values'''
            dlg = wx.MessageDialog(G2frame,'Are you sure you want to replace the saved mode values?',
                'Confirm replace',wx.YES|wx.NO)
            try:
                dlg.CenterOnParent()
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result != wx.ID_YES: return
            ISOdata['ISOmodeDispl'] = copy.deepcopy(ISOdata['modeDispl'])
            G2plt.PlotStructure(G2frame,data)
            UpdateISODISTORT(G2frame,data)

        #### displayModes code starts here
        ConstrData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Constraints'))
        pId = data['ranId']
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if SGLaue not in ['mmm','2/m','-1']:
            mainSizer.Add(wx.StaticText(G2frame.ISODIST,label=' NB: ISODISTORT distortion mode symmetry is too high to be used in PDFfit'))
        txt = wx.StaticText(G2frame.ISODIST,label=
                        ' For use of ISODISTORT, please cite: '+
                        G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...'))
        txt.Wrap(500)
        mainSizer.Add(txt)
        mainSizer.Add(wx.StaticText(G2frame.ISODIST,label=
u''' The 2nd column below shows the last saved mode values. The 3rd && 4th columns will set the
 display mode values. The positions in the Atoms and Draw Atoms tabs, as well as the atom
 positions shown in the Plot Window are changed to reflect the display mode values. The
 range of the slider corresponds to making a maximum atomic displacement between -2 && +2 \u212B.'''))
        mainSizer.Add((-1,10))
        slideSizer = wx.FlexGridSizer(0,5,0,0)
        modeDisp = ISOdata['modeDispl']
        idsp = 0
        slideSizer.Add(wx.StaticText(G2frame.ISODIST,label='Name'),0,wx.ALIGN_CENTER)
        slideSizer.Add(wx.StaticText(G2frame.ISODIST,label='Save value'))
        slideSizer.Add(wx.StaticText(G2frame.ISODIST,label='Value'),0,wx.ALIGN_CENTER)
        slideSizer.Add(wx.StaticText(G2frame.ISODIST,label='Refine?'),0,wx.ALIGN_CENTER)
        slideSizer.Add(wx.StaticText(G2frame.ISODIST,label='Atom displacements'),0,wx.EXPAND|wx.LEFT,15)
        isoDict = {i.name:j for (i,j) in zip(data['ISODISTORT']['G2ModeList'],data['ISODISTORT']['IsoModeList'])}
        for item in ConstrData['Phase']:
            if item[-1] != 'f': continue # only want new vars
            if item[-3] is None: continue # unnamed new var is not ISO
            try:
                if pId != item[-3].phase: continue # at present only ISO modes are associated with a phase
            except AttributeError:
                continue
            if  item[-3].name not in isoDict: continue
            isoName = item[-3].varname().split('::')[1]
            slideSizer.Add(wx.StaticText(G2frame.ISODIST,label=isoName),0,WACV)
            slideSizer.Add(wx.StaticText(G2frame.ISODIST,label=' %.5g '%ISOdata['ISOmodeDispl'][idsp],
                style=wx.ALIGN_CENTER_HORIZONTAL),0,WACV|wx.EXPAND)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            dispVal = G2G.ValidatedTxtCtrl(G2frame.ISODIST,modeDisp,idsp,xmin=-2.,xmax=2.,size=(75,20),OnLeave=OnDispVal)
            lineSizer.Add(dispVal,0,WACV)
            displ = G2G.G2Slider(G2frame.ISODIST,style=wx.SL_HORIZONTAL,minValue=-2000,maxValue=2000,
                value=int(modeDisp[idsp]*1000),size=(250,20))
            displ.Bind(wx.EVT_SLIDER, OnDispl)
            Indx[displ.GetId()] = [idsp,dispVal]
            Indx[dispVal.GetId()] = [idsp,displ]
            lineSizer.Add(displ)
            slideSizer.Add(lineSizer)
            refDispl = wx.CheckBox(G2frame.ISODIST)
            refDispl.SetValue(item[-2])
            refDispl.Bind(wx.EVT_CHECKBOX,OnRefDispl)
            Indx[refDispl.GetId()] = [idsp,item]
            slideSizer.Add(refDispl,0,WACV|wx.EXPAND|wx.LEFT,15)
            slideSizer.Add(wx.StaticText(G2frame.ISODIST,label=', '.join(ModeDispList[idsp])),0,wx.EXPAND|wx.LEFT,15)
            idsp += 1
        slideSizer.SetMinSize(wx.Size(650,10))
        mainSizer.Add(slideSizer)
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        reset = wx.Button(G2frame.ISODIST,label='Reset modes to save values')
        reset.Bind(wx.EVT_BUTTON,OnReset)
        lineSizer.Add(reset,0,WACV)
        reset = wx.Button(G2frame.ISODIST,label='Set all modes to zero')
        reset.Bind(wx.EVT_BUTTON,OnSetZero)
        lineSizer.Add(reset,0,wx.ALL,10)
        reset = wx.Button(G2frame.ISODIST,label='Save mode values')
        reset.Bind(wx.EVT_BUTTON,OnSaveModes)
        lineSizer.Add(reset,0,WACV)
        mainSizer.Add(lineSizer,0,wx.TOP,5)
        mainSizer.Layout()
        G2phsG.SetPhaseWindow(G2frame.ISODIST,mainSizer,Scroll=Scroll)

    #### UpdateISODISTORT code starts here
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl= f"ISODISTORT distortion modes for {data['General']['Name']!r}"[:60]
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    Indx = {}
    ISOdata = data['ISODISTORT']
    SGLaue = data['General']['SGData']['SGLaue']
    G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_ISODNEWPHASE,'rundata' in ISOdata)
    G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_ISOPDFFIT,(('G2VarList' in ISOdata) and (SGLaue in ['mmm','2/m','-1'])))
    G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_SHOWISO1,('G2VarList' in ISOdata)
        or ('G2OccVarList' in ISOdata))
    G2frame.dataWindow.ISODDataEdit.Enable(G2G.wxID_SHOWISOMODES,('G2VarList' in ISOdata))

    if G2frame.ISODIST.GetSizer():
        G2frame.ISODIST.GetSizer().Clear(True)

    if 'G2ModeList' in ISOdata:      #invoked only if phase is from a ISODISTORT cif file & thus contains distortion mode constraints

# #patch
#             if 'modeDispl' not in ISOdata:
#                 ISOdata['modeDispl'] = np.zeros(len(ISOdata['G2ModeList']))
# #end patch
        ModeDispList = G2pwd.GetAtmDispList(ISOdata)
        displayModes()
        return

#initialization
    if 'ParentCIF' not in ISOdata:
        ISOdata.update({'ParentCIF':'Select','ChildCIF':'Select','ISOmethod':4,
            'ChildMatrix':np.eye(3),'ChildSprGp':'P 1','ChildCell':'abc',})         #these last 3 currently unused
#end initialization

    mainSizer = wx.BoxSizer(wx.VERTICAL)
    txt = wx.StaticText(G2frame.ISODIST,label=
                        ' For use of ISODISTORT, please cite: '+
                        G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...'))
    txt.Wrap(500)
    mainSizer.Add(txt)
    mainSizer.Add((-1,5))
    G2G.HorizontalLine(mainSizer,G2frame.ISODIST)
    mainSizer.Add((-1,5))
    mainSizer.Add(displaySetup())

    if 'radio' in ISOdata:
        mainSizer.Add(displaySubset())
        mainSizer.Add(displayRadio())
    G2phsG.SetPhaseWindow(G2frame.ISODIST,mainSizer,Scroll=Scroll)

#### UpdateLayerData GUI for DIFFax Layer Data ################################################################################
def UpdateLayerData(G2frame,data,Scroll=0):
    '''Present the contents of the Phase/Layers tab for stacking fault simulation
    '''

    laueChoice = ['-1','2/m(ab)','2/m(c)','mmm','-3','-3m','4/m','4/mmm',
        '6/m','6/mmm','unknown']
    colLabels = ['Name','Type','x','y','z','frac','Uiso']
    transLabels = ['Prob','Dx','Dy','Dz','refine','plot']
    colTypes = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,]+ \
        3*[wg.GRID_VALUE_FLOAT+':10,5',]+2*[wg.GRID_VALUE_FLOAT+':10,4',] #x,y,z,frac,Uiso
    transTypes = [wg.GRID_VALUE_FLOAT+':10,3',]+3*[wg.GRID_VALUE_FLOAT+':10,5',]+ \
        [wg.GRID_VALUE_CHOICE+": ,P,Dx,Dy,Dz,Dxy,Dxz,Dyz,Dxyz",wg.GRID_VALUE_BOOL,]
    plotDefaults = {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':30.,'viewDir':[0,0,1],
        'viewPoint':[[0.,0.,0.],[]],}
    Indx = {}

    def OnLaue(event):
        Obj = event.GetEventObject()
        data['Layers']['Laue'] = Obj.GetValue()
        wx.CallAfter(UpdateLayerData,G2frame,data)

    def OnSadpPlot(event):
        sadpPlot.SetValue(False)
        labels = Layers['Sadp']['Plane']
        lmax = float(Layers['Sadp']['Lmax'])
        XY = 2*lmax*np.mgrid[0:256:256j,0:256:256j]/256.-lmax
        G2frame.Cmax = 1.0
        G2plt.PlotXYZ(G2frame,XY,Layers['Sadp']['Img'].T,labelX=labels[:-1],
            labelY=labels[-1],newPlot=False,Title=Layers['Sadp']['Plane'])

    def OnSeqPlot(event):
        seqPlot.SetValue(False)
        resultXY,resultXY2,seqNames = Layers['seqResults']
        pName = Layers['seqCodes'][0]
        G2plt.PlotXY(G2frame,resultXY,XY2=resultXY2,labelX=r'$\mathsf{2\theta}$',
            labelY='Intensity',newPlot=True,Title='Sequential simulations on '+pName,
            lines=True,names=seqNames)

    def CellSizer():

        cellGUIlist = [
            [['-3','-3m','6/m','6/mmm','4/m','4/mmm'],6,zip([" a = "," c = "],["%.5f","%.5f",],[True,True],[0,2])],
            [['mmm'],8,zip([" a = "," b = "," c = "],["%.5f","%.5f","%.5f"],[True,True,True],[0,1,2,])],
            [['2/m(ab)','2/m(c)','-1','axial','unknown'],10,zip([" a = "," b = "," c = "," gamma = "],
                ["%.5f","%.5f","%.5f","%.3f"],[True,True,True,True],[0,1,2,5])]]

        def OnCellRef(event):
            data['Layers']['Cell'][0] = cellRef.GetValue()

        def OnCellChange(event):
            event.Skip()
            laue = data['Layers']['Laue']
            cell = data['Layers']['Cell']
            Obj = event.GetEventObject()
            ObjId = cellList.index(Obj.GetId())
            try:
                value = max(1.0,float(Obj.GetValue()))
            except ValueError:
                if ObjId < 3:               #bad cell edge - reset
                    value = cell[ObjId+1]
                else:                       #bad angle
                    value = 90.
            if laue in ['-3','-3m','6/m','6/mmm','4/m','4/mmm']:
                cell[4] = cell[5] = 90.
                cell[6] = 120.
                if laue in ['4/m','4/mmm']:
                    cell[6] = 90.
                if ObjId == 0:
                    cell[1] = cell[2] = value
                    Obj.SetValue("%.5f"%(cell[1]))
                else:
                    cell[3] = value
                    Obj.SetValue("%.5f"%(cell[3]))
            elif laue in ['mmm']:
                cell[ObjId+1] = value
                cell[4] = cell[5] = cell[6] = 90.
                Obj.SetValue("%.5f"%(cell[ObjId+1]))
            elif laue in ['2/m','-1']:
                cell[4] = cell[5] = 90.
                if ObjId != 3:
                    cell[ObjId+1] = value
                    Obj.SetValue("%.5f"%(cell[ObjId+1]))
                else:
                    cell[6] = value
                    Obj.SetValue("%.3f"%(cell[6]))
            cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
            volVal.SetLabel(' Vol = %.3f'%(cell[7]))

        cell = data['Layers']['Cell']
        laue = data['Layers']['Laue']
        for cellGUI in cellGUIlist:
            if laue in cellGUI[0]:
                useGUI = cellGUI
        cellSizer = wx.FlexGridSizer(0,useGUI[1]+1,5,5)
        cellRef = wx.CheckBox(layerData,-1,label='Refine unit cell:')
        cellSizer.Add(cellRef,0,WACV)
        cellRef.Bind(wx.EVT_CHECKBOX, OnCellRef)
        cellRef.SetValue(cell[0])
        cellList = []
        for txt,fmt,ifEdit,Id in useGUI[2]:
            cellSizer.Add(wx.StaticText(layerData,label=txt),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            cellVal = wx.TextCtrl(layerData,value=(fmt%(cell[Id+1])),
                style=wx.TE_PROCESS_ENTER)
            cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)
            cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
            cellSizer.Add(cellVal,0,WACV)
            cellList.append(cellVal.GetId())
        volVal = wx.StaticText(layerData,label=' Vol = %.3f'%(cell[7]))
        cellSizer.Add(volVal,0,WACV)
        return cellSizer

    def WidthSizer():

        def OnRefWidth(event):
            Id = Indx[event.GetEventObject()]
            Layers['Width'][1][Id] = not Layers['Width'][1][Id]

        Labels = ['a','b']
        flags = Layers['Width'][1]
        widthSizer = wx.BoxSizer(wx.HORIZONTAL)
        for i in range(2):
            widthSizer.Add(wx.StaticText(layerData,label=u' layer width(%s) (<= 1\xb5m): '%(Labels[i])),0,WACV)
            widthVal = G2G.ValidatedTxtCtrl(layerData,Layers['Width'][0],i,nDig=(10,3),xmin=0.005,xmax=1.0)
            widthSizer.Add(widthVal,0,WACV)
            widthRef = wx.CheckBox(layerData,label='Refine?')
            widthRef.SetValue(flags[i])
            Indx[widthRef] = i
            widthRef.Bind(wx.EVT_CHECKBOX, OnRefWidth)
            widthSizer.Add(widthRef,0,WACV)
        return widthSizer

    def OnNewLayer(event):
        data['Layers']['Layers'].append({'Name':'Unk','SameAs':'','Symm':'None','Atoms':[]})
        Trans = data['Layers']['Transitions']
        if len(Trans):
            Trans.append([[0.,0.,0.,0.,'',False] for trans in Trans])
            for trans in Trans:
                trans.append([0.,0.,0.,0.,'',False])
        else:
            Trans = [[[1.,0.,0.,0.,'',False],],]
        data['Layers']['Transitions'] = Trans
        wx.CallLater(100,UpdateLayerData,G2frame,data)

    def OnDeleteLast(event):
        del(data['Layers']['Layers'][-1])
        del(data['Layers']['Transitions'][-1])
        for trans in data['Layers']['Transitions']:
            del trans[-1]
        wx.CallAfter(UpdateLayerData,G2frame,data)

    def OnImportLayer(event):
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II project file', G2G.GetImportPath(G2frame),
            wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN| wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                GPXFile = dlg.GetPath()
                phaseNames = G2stIO.GetPhaseNames(GPXFile)
            else:
                return
        finally:
            dlg.Destroy()
        dlg = wx.SingleChoiceDialog(G2frame,'Phase to use for layer','Select',phaseNames)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            PhaseName = phaseNames[sel]
        else:
            return
        Phase = G2stIO.GetAllPhaseData(GPXFile,PhaseName)
        #need cell compatibility check here
        Layer = {'Name':Phase['General']['Name'],'SameAs':'','Symm':'None'}
        cx,ct,cs,cia = Phase['General']['AtomPtrs']
        atoms = Phase['Atoms']
        Atoms = []
        for atom in atoms:
            x,y,z,f = atom[cx:cx+4]
            u = atom[cia+1]
            if not u: u = 0.01
            Atoms.append([atom[ct-1],atom[ct],x,y,z,f,u])
            if atom[ct] not in data['Layers']['AtInfo']:
                data['Layers']['AtInfo'][atom[ct]] = G2elem.GetAtomInfo(atom[ct])
        Layer['Atoms'] = Atoms
        data['Layers']['Layers'].append(Layer)
        Trans = data['Layers']['Transitions']
        if len(Trans):
            Trans.append([[0.,0.,0.,0.,'',False] for trans in Trans])
            for trans in Trans:
                trans.append([0.,0.,0.,0.,'',False])
        else:
            Trans = [[[1.,0.,0.,0.,'',False],],]
        data['Layers']['Transitions'] = Trans
        wx.CallAfter(UpdateLayerData,G2frame,data)

    def LayerSizer(il,Layer):

        # def OnNameChange(event):
        #     event.Skip()
        #     Layer['Name'] = layerName.GetValue()
        #     wx.CallLater(100,UpdateLayerData,G2frame,data)

        def OnAddAtom(event):
            Layer['Atoms'].append(['Unk','Unk',0.,0.,0.,1.,0.01])
            wx.CallAfter(UpdateLayerData,G2frame,data)

        def OnSymm(event):
            Layer['Symm'] = symm.GetValue()

        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if atomGrid.GetColLabelValue(c) == 'Type':
                PE = G2elemGUI.PickElement(G2frame)
                if PE.ShowModal() == wx.ID_OK:
                    if PE.Elem != 'None':
                        atType = PE.Elem.strip()
                        Layer['Atoms'][r][c] = atType
                        name = Layer['Atoms'][r][c]
                        if len(name) in [2,4]:
                            Layer['Atoms'][r][c-1] = name[:2]+'%d'%(r+1)
                        else:
                            Layer['Atoms'][r][c-1] = name[:1]+'%d'%(r+1)
                        if atType not in data['Layers']['AtInfo']:
                            data['Layers']['AtInfo'][atType] = G2elem.GetAtomInfo(atType)
                PE.Destroy()
                wx.CallAfter(UpdateLayerData,G2frame,data)
            else:
                event.Skip()

        def OnDrawLayer(event):
            drawLayer.SetValue(False)
            G2plt.PlotLayers(G2frame,Layers,[il,],plotDefaults,firstCall=True)

        def OnSameAs(event):
            Layer['SameAs'] = sameas.GetValue()
            wx.CallLater(100,UpdateLayerData,G2frame,data)

        layerSizer = wx.BoxSizer(wx.VERTICAL)
        nameSizer = wx.BoxSizer(wx.HORIZONTAL)
        nameSizer.Add(wx.StaticText(layerData,label=' Layer name: '),0,WACV)
        layerName = G2G.ValidatedTxtCtrl(layerData,Layer,'Name')
        # layerName = wx.TextCtrl(layerData,value=Layer['Name'],style=wx.TE_PROCESS_ENTER)
        # layerName.Bind(wx.EVT_TEXT_ENTER,OnNameChange)
        # layerName.Bind(wx.EVT_KILL_FOCUS,OnNameChange)
        # layerName.Bind(wx.EVT_LEAVE_WINDOW,OnNameChange)
        nameSizer.Add(layerName,0,WACV)
        if il:
            nameSizer.Add(wx.StaticText(layerData,label=' Same as: '),0,WACV)
            sameas = wx.ComboBox(layerData,value=Layer['SameAs'],choices=['',]+layerNames[:-1],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            sameas.Bind(wx.EVT_COMBOBOX, OnSameAs)
            nameSizer.Add(sameas,0,WACV)
            if Layer['SameAs']:
                indx = layerNames.index(Layer['SameAs'])
                if indx < il:    #previously used : same layer
                    layerSizer.Add(nameSizer)
                    return layerSizer
        nameSizer.Add(wx.StaticText(layerData,label=' Layer symmetry: '),0,WACV)
        symmChoice = ['-1','None']
        symm = wx.ComboBox(layerData,value=Layer['Symm'],choices=symmChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        symm.Bind(wx.EVT_COMBOBOX,OnSymm)
        nameSizer.Add(symm,0,WACV)
        addAtom = wx.CheckBox(layerData,label=' Add atom? ')
        addAtom.Bind(wx.EVT_CHECKBOX, OnAddAtom)
        nameSizer.Add(addAtom,0,WACV)
        drawLayer = wx.CheckBox(layerData,label=' Draw layer? ')
        drawLayer.Bind(wx.EVT_CHECKBOX, OnDrawLayer)
        nameSizer.Add(drawLayer,0,WACV)
        layerSizer.Add(nameSizer)
        table = []
        rowLabels = []
        for i,atom in enumerate(Layer['Atoms']):
            table.append(atom)
            rowLabels.append(str(i))
        atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=colTypes)
        atomGrid = G2G.GSGrid(layerData)
        atomGrid.SetTable(atomTable,True,useFracEdit=False)
#            atomGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
        # loop over all cols in table, set cell editor for numerical items
        for c,t in enumerate(colTypes):
            if not t.startswith(wg.GRID_VALUE_FLOAT): continue
            attr = wx.grid.GridCellAttr()
            attr.IncRef()               #fix from Jim Hester
            attr.SetEditor(G2G.GridFractionEditor(atomGrid))
            atomGrid.SetColAttr(c, attr)
        for row,atom in enumerate(Layer['Atoms']):
            atomGrid.SetReadOnly(row,1,True)
        atomGrid.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
        atomGrid.AutoSizeColumns(True)
        layerSizer.Add(atomGrid)
        return layerSizer

    def TransSizer():

        def PlotSelect(event):
            Obj = event.GetEventObject()
            Yi = Indx[Obj.GetId()]
            Xi,c =  event.GetRow(),event.GetCol()
            if Xi >= 0 and c == 5:   #plot column
                G2plt.PlotLayers(G2frame,Layers,[Yi,Xi,],plotDefaults,firstCall=True)
            else:
                Psum = 0.
                for Xi in range(len(transArray)):
                    Psum += transArray[Xi][Xi][0]
                Psum /= len(transArray)
                totalFault.SetLabel(' Total fault density = %.3f'%(1.-Psum))
                event.Skip()

        def OnNormProb(event):
            for Yi,Yname in enumerate(Names):
                Psum = 0.
                for Xi,Xname in enumerate(Names):
                    Psum += transArray[Yi][Xi][0]
                if not Psum:
                    transArray[Yi][0][0] = 1.0
                    Psum = 1.0
                for Xi,Xname in enumerate(Names):
                    transArray[Yi][Xi][0] /= Psum
            wx.CallAfter(UpdateLayerData,G2frame,data)

        def OnSymProb(event):
            if symprob.GetValue():
                Nx = len(Names)-1
                Layers['SymTrans'] = True
                for Yi,Yname in enumerate(Names):
                    for Xi,Xname in enumerate(Names):
                        if transArray[Nx-Yi][Nx-Xi][0] != transArray[Yi][Xi][0]:
                            Layers['SymTrans'] = False
                            symprob.SetValue(False)
                            wx.MessageBox('%s-%s not equal %s-%s'%(Yname,Xname,Xname,Yname),
                                caption='Probability symmetry error',style=wx.ICON_EXCLAMATION)
                            break
            else:
                Layers['SymTrans'] = False

        transSizer = wx.BoxSizer(wx.VERTICAL)
        transSizer.Add(wx.StaticText(layerData,label=' Layer-Layer transition probabilities: '),0)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        normprob = wx.CheckBox(layerData,label=' Normalize probabilities?')
        normprob.Bind(wx.EVT_CHECKBOX,OnNormProb)
        topSizer.Add(normprob,0,WACV)
        symprob = wx.CheckBox(layerData,label=' Symmetric probabilities?')
        symprob.SetValue(Layers.get('SymTrans',False))
        symprob.Bind(wx.EVT_CHECKBOX,OnSymProb)
        topSizer.Add(symprob,0,WACV)
        transSizer.Add(topSizer,0)
        Names = [layer['Name'] for layer in Layers['Layers']]
        transArray = Layers['Transitions']
        layerData.transGrids = []
        if not Names or not transArray:
            return transSizer
        diagSum = 0.
        for Yi,Yname in enumerate(Names):
            transSizer.Add(wx.StaticText(layerData,label=' From %s to:'%(Yname)),0)
            table = []
            rowLabels = []
            diagSum += transArray[Yi][Yi][0]
            for Xi,Xname in enumerate(Names):
                table.append(transArray[Yi][Xi])
                rowLabels.append(Xname)
                if transArray[Yi][Xi][0] > 0.:
                    Layers['allowedTrans'].append([str(Yi+1),str(Xi+1)])
            transTable = G2G.Table(table,rowLabels=rowLabels,colLabels=transLabels,types=transTypes)
            transGrid = G2G.GSGrid(layerData)
            transGrid.SetTable(transTable,True,useFracEdit=False)
#                transGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
            Indx[transGrid.GetId()] = Yi
            for c,t in enumerate(transTypes):
                if not t.startswith(wg.GRID_VALUE_FLOAT): continue
                attr = wx.grid.GridCellAttr()
                attr.IncRef()               #fix from Jim Hester
                attr.SetEditor(G2G.GridFractionEditor(transGrid))
                transGrid.SetColAttr(c, attr)
            transGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, PlotSelect)
            transGrid.AutoSizeColumns(True)
            transSizer.Add(transGrid)
            layerData.transGrids.append(transGrid)
        if len(transArray):
            diagSum /= len(transArray)
            totalFault = wx.StaticText(layerData,
                label=' Total fault density = %.3f'%(1.-diagSum))
            transSizer.Add(totalFault,0)
        return transSizer

    def PlotSizer():

        def OnPlotSeq(event):
            event.Skip()
            vals = plotSeq.GetValue().split()
            try:
                vals = [int(val)-1 for val in vals]
                if not all([0 <= val < len(Names) for val in vals]):
                    raise ValueError
            except ValueError:
                plotSeq.SetValue('Error in string '+plotSeq.GetValue())
                return
            G2plt.PlotLayers(G2frame,Layers,vals,plotDefaults,firstCall=True)

        Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
        plotSizer = wx.BoxSizer(wx.VERTICAL)
        Str = ' Using sequence nos. from:'
        for name in Names:
            Str += name
        plotSizer.Add(wx.StaticText(layerData,label=Str[:-1]),0)
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(layerData,label=' Enter sequence of layers to plot:'),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
        plotSeq = wx.TextCtrl(layerData,value = '',style=wx.TE_PROCESS_ENTER)
        plotSeq.Bind(wx.EVT_TEXT_ENTER,OnPlotSeq)
        plotSeq.Bind(wx.EVT_KILL_FOCUS,OnPlotSeq)
        lineSizer.Add(plotSeq,0,WACV)
        plotSizer.Add(lineSizer,0)
        return plotSizer

    def StackSizer():

        stackChoice = ['recursive','explicit',]
        seqChoice = ['random','list',]

        def OnStackType(event):
            newType = stackType.GetValue()
            if newType == data['Layers']['Stacking'][0]:
                return
            data['Layers']['Stacking'][0] = newType
            if newType == 'recursive':
                data['Layers']['Stacking'][1] = 'infinite'
            else:  #explicit
                data['Layers']['Stacking'][1] = 'random'
                data['Layers']['Stacking'][2] = '250'
            wx.CallAfter(UpdateLayerData,G2frame,data)

        def OnSeqType(event):
            newType = seqType.GetValue()
            if newType == data['Layers']['Stacking'][1]:
                return
            data['Layers']['Stacking'][1] = newType
            if newType == 'random':
                data['Layers']['Stacking'][2] = '250'
            else: #List
                data['Layers']['Stacking'][2] = ''
            wx.CallAfter(UpdateLayerData,G2frame,data)

        def OnNumLayers(event):
            event.Skip()
            val = numLayers.GetValue()
            if val == 'infinite':
                data['Layers']['Stacking'][1] = val
            else:
                try:
                    if 0 < int(val) < 1023:
                        data['Layers']['Stacking'][1] = val
                    else:
                        data['Layers']['Stacking'][1] = 'infinite'
                except ValueError:
                    pass
            numLayers.SetValue(data['Layers']['Stacking'][1])

        def OnNumRan(event):
            event.Skip()
            val = numRan.GetValue()
            try:
                if 0 > int(val) > 1022:
                    raise ValueError
                else:
                    data['Layers']['Stacking'][2] = val
            except ValueError:
                val = data['Layers']['Stacking'][2]
            numRan.SetValue(val)

        def OnStackList(event):
            event.Skip()
            stack = stackList.GetValue()
            stack = stack.replace('\n',' ').strip().strip('\n')
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
            Slist = stack.split()
            if len(Slist) < 2:
                stack = 'Error in sequence - too short!'
            OKlist = [Slist[i:i+2] in Layers['allowedTrans'] for i in range(len(Slist[:-1]))]
            if all(OKlist):
                data['Layers']['Stacking'][2] = stack
            else:
                stack = 'Improbable sequence or bad string'
            stackList.SetValue(stack)

        stackSizer = wx.BoxSizer(wx.VERTICAL)
        stackSizer.Add(wx.StaticText(layerData,label=' Layer stacking parameters:'),0)
        if not Layers['Stacking']:
            Layers['Stacking'] = ['recursive','infinite','']
        topLine = wx.BoxSizer(wx.HORIZONTAL)
        topLine.Add(wx.StaticText(layerData,label=' Stacking type: '),0,WACV)
        stackType = wx.ComboBox(layerData,value=Layers['Stacking'][0],choices=stackChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        stackType.Bind(wx.EVT_COMBOBOX,OnStackType)
        topLine.Add(stackType,0,WACV)
        if Layers['Stacking'][0] == 'recursive':
            topLine.Add(wx.StaticText(layerData,label=' number of layers (<1022 or "infinite"): '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            numLayers = wx.TextCtrl(layerData,value=data['Layers']['Stacking'][1],style=wx.TE_PROCESS_ENTER)
            numLayers.Bind(wx.EVT_TEXT_ENTER,OnNumLayers)
            numLayers.Bind(wx.EVT_KILL_FOCUS,OnNumLayers)
            topLine.Add(numLayers,0,WACV)
            stackSizer.Add(topLine)
        elif Layers['Stacking'][0] == 'explicit':
            topLine.Add(wx.StaticText(layerData,label=' layer sequence: '),0,WACV)
            seqType = wx.ComboBox(layerData,value=data['Layers']['Stacking'][1],choices=seqChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            seqType.Bind(wx.EVT_COMBOBOX,OnSeqType)
            topLine.Add(seqType,0,WACV)
            if Layers['Stacking'][1] == 'list':
                stackSizer.Add(topLine,0)
                Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
                stackSizer.Add(wx.StaticText(layerData,label=' Explicit layer sequence; enter space delimited list of numbers:'),0)
                Str = ' Use sequence nos. from:'
                for name in Names:
                    Str += name
                stackSizer.Add(wx.StaticText(layerData,label=Str[:-1]+' Repeat sequences can be used: e.g. 6*(1 2) '),0)
                stackSizer.Add(wx.StaticText(layerData,label=' Zero probability sequences not allowed'),0)
                stackList = wx.TextCtrl(layerData,value=Layers['Stacking'][2],size=(600,-1),
                    style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER)
                stackList.Bind(wx.EVT_TEXT_ENTER,OnStackList)
                stackList.Bind(wx.EVT_KILL_FOCUS,OnStackList)
                stackSizer.Add(stackList,0,wx.ALL|wx.EXPAND,8)
            else:   #random
                topLine.Add(wx.StaticText(layerData,label=' Length of random sequence: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                numRan = wx.TextCtrl(layerData,value=Layers['Stacking'][2],style=wx.TE_PROCESS_ENTER)
                numRan.Bind(wx.EVT_TEXT_ENTER,OnNumRan)
                numRan.Bind(wx.EVT_KILL_FOCUS,OnNumRan)
                topLine.Add(numRan,0,WACV)
                stackSizer.Add(topLine,0)
        return stackSizer

    Layers = data['Layers']
    layerNames = []
    Layers['allowedTrans'] = []
    if len(Layers['Layers']):
        layerNames = [layer['Name'] for layer in Layers['Layers']]
    G2frame.GetStatusBar().SetStatusText('',1)
    layerData = G2frame.layerData
    try:
        if layerData.GetSizer():
            layerData.GetSizer().Clear(True)
    except:
        pass
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.VERTICAL)
    bottomSizer = wx.BoxSizer(wx.VERTICAL)
    headSizer = wx.BoxSizer(wx.HORIZONTAL)
    headSizer.Add(wx.StaticText(layerData,label=' Global layer description:  '),0,WACV)
    if 'Sadp' in Layers:
        sadpPlot = wx.CheckBox(layerData,label=' Plot selected area diffraction?')
        sadpPlot.Bind(wx.EVT_CHECKBOX,OnSadpPlot)
        headSizer.Add(sadpPlot,0,WACV)
    if 'seqResults' in Layers:
        seqPlot = wx.CheckBox(layerData,label=' Plot sequential result?')
        seqPlot.Bind(wx.EVT_CHECKBOX,OnSeqPlot)
        headSizer.Add(seqPlot,0,WACV)
    topSizer.Add(headSizer)
    laueSizer = wx.BoxSizer(wx.HORIZONTAL)
    laueSizer.Add(wx.StaticText(layerData,label=' Diffraction Laue symmetry:'),0,WACV)
    laue = wx.ComboBox(layerData,value=Layers['Laue'],choices=laueChoice,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    laue.Bind(wx.EVT_COMBOBOX,OnLaue)
    laueSizer.Add(laue,0,WACV)
    if Layers['Laue'] == 'unknown':
        laueSizer.Add(wx.StaticText(layerData,label=' Diffraction symmetry tolerance: '),0,WACV)
        toler = G2G.ValidatedTxtCtrl(layerData,Layers,'Toler',nDig=(10,3))
        laueSizer.Add(toler,0,WACV)
    topSizer.Add(laueSizer,0)
    topSizer.Add(wx.StaticText(layerData,label=' Reference unit cell for all layers:'),0)
    topSizer.Add(CellSizer(),0)
    topSizer.Add(WidthSizer())
    topSizer.Add(wx.StaticText(layerData,label=' NB: stacking fault refinement currently not available'),0)
    G2G.HorizontalLine(topSizer,layerData)
    titleSizer = wx.BoxSizer(wx.HORIZONTAL)
    titleSizer.Add(wx.StaticText(layerData,label=' Layer descriptions: '),0,WACV)
    newLayer = wx.Button(layerData,label='Add new layer')
    newLayer.Bind(wx.EVT_BUTTON, OnNewLayer)
    titleSizer.Add(newLayer,0,WACV)
    importLayer = wx.Button(layerData,label=' Import new layer')
    importLayer.Bind(wx.EVT_BUTTON, OnImportLayer)
    titleSizer.Add(importLayer,0,WACV)
    deleteLast = wx.Button(layerData,label='Delete last layer')
    deleteLast.Bind(wx.EVT_BUTTON, OnDeleteLast)
    titleSizer.Add(deleteLast,0,WACV)
    topSizer.Add(titleSizer,0)
    for il,layer in enumerate(Layers['Layers']):
        topSizer.Add(LayerSizer(il,layer))
    G2G.HorizontalLine(topSizer,layerData)
    mainSizer.Add(topSizer)
    bottomSizer.Add(TransSizer())
    G2G.HorizontalLine(bottomSizer,layerData)
    bottomSizer.Add(PlotSizer(),0)
    G2G.HorizontalLine(bottomSizer,layerData)
    bottomSizer.Add(StackSizer())
    mainSizer.Add(bottomSizer)
    G2phsG.SetPhaseWindow(G2frame.layerData,mainSizer,Scroll=Scroll)
#### UpdateTexture GUI
def UpdateTexture(G2frame,data):

    def SetSHCoef():
        cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],SamSym[textureData['Model']],textureData['Order'])
        newSHCoef = dict(zip(cofNames,np.zeros(len(cofNames))))
        SHCoeff = textureData['SH Coeff'][1]
        for cofName in SHCoeff:
            if cofName in  cofNames:
                newSHCoef[cofName] = SHCoeff[cofName]
        return newSHCoef

    def OnShOrder(event):
        Obj = event.GetEventObject()
        # the Kaduk test: is Texture appropriate? Look for a 2nd histogram
        # of any type with differing setting angles from the 1st.
        instArray = {}
        if int(Obj.GetValue()) > 0:
            for h in data['Histograms']:
                PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,h)
                if not PatternId:       #skip bogus histograms
                    continue
                Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Instrument Parameters'))[0]
                Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Sample Parameters'))
                if Inst['Type'][1] in instArray:
                    if instArray[Inst['Type'][1]] != [Sample['Chi'],Sample['Phi'],Sample['Omega']]:
                        textureData['Order'] = int(Obj.GetValue())
                        break
                else:
                    instArray[Inst['Type'][1]] = [Sample['Chi'],Sample['Phi'],Sample['Omega']]
            else:
                if textureData['Model'] != 'cylindrical':
                    textureData['Order'] = 0
                    wx.MessageBox('Incorrect use of Texture. Use preferred orientation (on data tab) unless you have multiple histograms taken with different orientations',
                            caption='Texture Error',style=wx.ICON_EXCLAMATION)
                else:
                    textureData['Order'] = int(Obj.GetValue())
        else:
            textureData['Order'] = 0
        textureData['SH Coeff'][1] = SetSHCoef()
        G2frame.GPXtree.UpdateSelection()
#            wx.CallLater(100,UpdateTexture)
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnShModel(event):
        Obj = event.GetEventObject()
        textureData['Model'] = Obj.GetValue()
        textureData['SH Coeff'][1] = SetSHCoef()
        G2frame.GPXtree.UpdateSelection()
#            wx.CallLater(100,UpdateTexture)
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnSHRefine(event):
        Obj = event.GetEventObject()
        textureData['SH Coeff'][0] = Obj.GetValue()

    def OnSHShow(event):
        Obj = event.GetEventObject()
        textureData['SHShow'] = Obj.GetValue()
        G2frame.GPXtree.UpdateSelection()
#            wx.CallLater(100,UpdateTexture)

    def OnProjSel(event):
        Obj = event.GetEventObject()
        G2frame.Projection = Obj.GetValue()
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnShoDet(event):
        Obj = event.GetEventObject()
        textureData['ShoDet'] = Obj.GetValue()
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnColorSel(event):
        Obj = event.GetEventObject()
        G2frame.ContourColor = Obj.GetValue()
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnAngRef(event):
        Obj = event.GetEventObject()
        textureData[angIndx[Obj.GetId()]][0] = Obj.GetValue()

    def OnODFValue(invalid,value,tc):
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnPfType(event):
        Obj = event.GetEventObject()
        textureData['PlotType'] = Obj.GetValue()
        G2frame.GPXtree.UpdateSelection()
#            wx.CallLater(100,UpdateTexture)
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnPFValue(event):
        event.Skip()
        Obj = event.GetEventObject()
        Saxis = Obj.GetValue().split()
        if textureData['PlotType'] in ['Pole figure','Axial pole distribution','3D pole distribution']:
            try:
                hkl = [int(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                hkl = textureData['PFhkl']
            if not np.any(np.array(hkl)):       #can't be all zeros!
                hkl = textureData['PFhkl']
            Obj.SetValue('%d %d %d'%(hkl[0],hkl[1],hkl[2]))
            textureData['PFhkl'] = hkl
        else:
            try:
                xyz = [float(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                xyz = textureData['PFxyz']
            if not np.any(np.array(xyz)):       #can't be all zeros!
                xyz = textureData['PFxyz']
            Obj.SetValue('%3.1f %3.1f %3.1f'%(xyz[0],xyz[1],xyz[2]))
            textureData['PFxyz'] = xyz
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    def OnpopLA(event):
        pfName = PhaseName
        cell = generalData['Cell'][1:7]
        PH = np.array(textureData['PFhkl'])
        phi,beta = G2lat.CrsAng(PH,cell,SGData)
        SHCoef = textureData['SH Coeff'][1]
        ODFln = G2lat.Flnh(SHCoef,phi,beta,SGData)
        pfName = PhaseName+'%d%d%d.gpf'%(PH[0],PH[1],PH[2])
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose popLA pole figure file name', pth, pfName,
            'popLA file (*.gpf)|*.gpf',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                pfFile = dlg.GetPath()
        finally:
            dlg.Destroy()
        print ('popLA save '+pfFile)
        if pfFile:
            pf = open(pfFile,'w')
            pf.write(PhaseName+'\n')
            str = ' %d%d%d   5.0 90.0  5.0360.0 1 1 2 1 3  100    1'%(PH[0],PH[1],PH[2])
            pf.write(str+'\n')
            Psi,Gam = np.mgrid[0:19,0:72]
            Psi = Psi.flatten()*5.
            Gam = Gam.flatten()*5.
            Z = np.array(G2lat.polfcal(ODFln,SamSym[textureData['Model']],Psi,Gam)*100.,dtype='int')
            Z = np.where(Z>0,Z,0)
            Z = np.where(Z<9999,Z,9999)
            for i in range(76):
                iBeg = i*18
                iFin = iBeg+18
                np.savetxt(pf,Z[iBeg:iFin],fmt='%4d',newline='')
                pf.write('\n')
            pf.close()
            print (' popLA %d %d %d pole figure saved to %s'%(PH[0],PH[1],PH[2],pfFile))

    def OnCSV(event):
        pfName = PhaseName
        pfFile = ''
        cell = generalData['Cell'][1:7]
        if 'Inverse' in textureData['PlotType']:
            SHCoef = textureData['SH Coeff'][1]
            PX = np.array(textureData['PFxyz'])
            gam = atan2d(PX[0],PX[1])
            xy = np.sqrt(PX[0]**2+PX[1]**2)
            xyz = np.sqrt(PX[0]**2+PX[1]**2+PX[2]**2)
            psi = asind(xy/xyz)
            IODFln = G2lat.Glnh(SHCoef,psi,gam,SamSym[textureData['Model']])
            pfName = PhaseName+'%d%d%dIPF.csv'%(int(PX[0]),int(PX[1]),int(PX[2]))
            pth = G2G.GetExportPath(G2frame)
            dlg = wx.FileDialog(G2frame, 'Choose CSV inverse pole figure file name', pth, pfName,
                'CSV file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        else:
            PH = np.array(textureData['PFhkl'])
            phi,beta = G2lat.CrsAng(PH,cell,SGData)
            SHCoef = textureData['SH Coeff'][1]
            ODFln = G2lat.Flnh(SHCoef,phi,beta,SGData)
            pfName = PhaseName+'%d%d%dPF.csv'%(PH[0],PH[1],PH[2])
            pth = G2G.GetExportPath(G2frame)
            dlg = wx.FileDialog(G2frame, 'Choose CSV pole figure file name', pth, pfName,
                'CSV file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                pfFile = dlg.GetPath()
                print ('CSV save '+pfFile)
        finally:
            dlg.Destroy()
        if pfFile:
            pf = open(pfFile,'w')
            pf.write('"%s"\n'%(PhaseName))
            if 'Inverse' in textureData['PlotType']:
                pf.write('" %s %d %d %d inverse pole figure"\n'%(PhaseName,int(PX[0]),int(PX[1]),int(PX[2])))
                P,R = np.mgrid[0:19,0:72]
                pf.write('"phi/beta",')
                np.savetxt(pf,np.linspace(0.,90.,19,True),fmt='%10.4f,',newline='')
                pf.write('\n')
                P = P.flatten()*5.
                R = R.flatten()*5.
                Z = G2lat.invpolfcal(IODFln,SGData,P,R)
                Z = np.reshape(Z,(19,72)).T
                for i,row in enumerate(Z):
                    pf.write('%8d,  '%(i*5))
                    np.savetxt(pf,row,fmt='%10.4f,',newline='')
                    pf.write('\n')
                pf.close()
                print (' %s %d %d %d inverse pole figure saved to %s'%(PhaseName,int(PX[0]),int(PX[1]),int(PX[2]),pfFile))
            else:
                pf.write('" %s %d %d %d pole figure"\n'%(PhaseName,PH[0],PH[1],PH[2]))
                Psi,Gam = np.mgrid[0:19,0:72]
                pf.write('"psi/gam",')
                np.savetxt(pf,np.linspace(0.,90.,19,True),fmt='%10.4f,',newline='')
                pf.write('\n')
                Psi = Psi.flatten()*5.
                Gam = Gam.flatten()*5.
                Z = np.array(G2lat.polfcal(ODFln,SamSym[textureData['Model']],Psi,Gam))
                Z = np.reshape(Z,(19,72)).T
                for i,row in enumerate(Z):
                    pf.write('%8d, '%(i*5))
                    np.savetxt(pf,row,fmt='%10.4f,',newline='')
                    pf.write('\n')
                pf.close()
                print (' %s %d %d %d pole figure saved to %s'%(PhaseName,PH[0],PH[1],PH[2],pfFile))

    def SHPenalty(Penalty):

        def OnHKLList(event):
            event.Skip()
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select penalty hkls',
                'Penalty hkls',hkls,filterBox=False)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Penalty[0] = [hkls[i] for i in dlg.GetSelections()]
                    if not Penalty[0]:
                        Penalty[0] = ['',]
                else:
                    return
            finally:
                dlg.Destroy()
            G2frame.GPXtree.UpdateSelection()
#                wx.CallLater(100,UpdateTexture)

        A = G2lat.cell2A(generalData['Cell'][1:7])
        hkls = G2lat.GenPfHKLs(10,SGData,A)
        shPenalty = wx.BoxSizer(wx.HORIZONTAL)
        shPenalty.Add(wx.StaticText(Texture,wx.ID_ANY,' Negative MRD penalty list: '),0,WACV)
        shPenalty.Add(wx.ComboBox(Texture,value=Penalty[0][0],choices=Penalty[0],
            style=wx.CB_DROPDOWN|wx.CB_READONLY),0,WACV)
        hklList = wx.Button(Texture,label='Select penalty hkls')
        hklList.Bind(wx.EVT_BUTTON,OnHKLList)
        shPenalty.Add(hklList,0,WACV)
        shPenalty.Add(wx.StaticText(Texture,wx.ID_ANY,' Zero MRD tolerance: '),0,WACV)
        shToler = G2G.ValidatedTxtCtrl(Texture,Penalty,1,nDig=(10,2),xmin=0.001)
        shPenalty.Add(shToler,0,WACV)
        return shPenalty

    def OnProj(event):
        Obj = event.GetEventObject()
        generalData['3Dproj'] = Obj.GetValue()
        wx.CallAfter(G2plt.PlotTexture,G2frame,data)

    # UpdateTexture executable starts here
    Texture = G2frame.Texture
    generalData = data['General']
    PhaseName = generalData['Name']
    SGData = generalData['SGData']
    try:
        textureData = generalData['SH Texture']
    except KeyError:            #fix old files!
        textureData = generalData['SH Texture'] = {'Order':0,'Model':'cylindrical',
            'Sample omega':[False,0.0],'Sample chi':[False,0.0],'Sample phi':[False,0.0],
            'SH Coeff':[False,{}],'SHShow':False,'PFhkl':[0,0,1],'ShoDet':False,
            'PFxyz':[0,0,1.],'PlotType':'Pole figure'}
    if 'SHShow' not in textureData:
        textureData.update({'SHShow':False,'PFhkl':[0,0,1],'PFxyz':[0,0,1.],'PlotType':'Pole figure'})
    if 'PlotType' not in textureData:
        textureData.update({'PFxyz':[0,0,1.],'PlotType':'Pole figure'})
    if 'Penalty' not in textureData:
        textureData['Penalty'] = [['',],0.1,False,1.0]
    if 'ShoDet' not in textureData:
        textureData['ShoDet'] = False
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    keyList = G2frame.GetHistogramNames(['PWDR',])
    UseList = data['Histograms']
    HistsInPhase = [name for name in keyList if name in UseList]

    textureData['det Angles'] = []
    for hist in HistsInPhase:
        pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)       #only use 1st histogram
        Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Sample Parameters'))
        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Instrument Parameters'))[0]
        if 'NT' in Inst['Type'][0]:
            Tth = Inst['2-theta'][1]/2.
            Gangls = [Sample['Phi'],Sample['Chi'],Sample['Omega'],Sample['Azimuth']]
            Sangls = [textureData['Sample omega'][1],textureData['Sample chi'][1],textureData['Sample phi'][1]]
            phi,gam,x,x = G2lat.SamAng(Tth,Gangls,Sangls,False)
            if phi > 90.:
                phi = 180.-phi
                gam += 180.
            gam %= 360.
            textureData['det Angles'].append([hist,phi,gam])

    G2frame.GetStatusBar().SetStatusText('',1)
    if Texture.GetSizer():
        Texture.GetSizer().Clear(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    # sanity check: should this project be fitting texture?
    mainSizer.Add(wx.StaticText(Texture,label=
        ' NB: Normally texture model fitting generally requires multiple datasets with differing sample orientations/detector values'))
    if G2frame.testSeqRefineMode() and G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential results'):
        mainSizer.Add(wx.StaticText(Texture,label=
            " Sequential result found. Use Texture/Refine texture above. See Method B in texture tutorial."))
    h,pd = G2frame.GetUsedHistogramsAndPhasesfromTree()
    if len(h) < 3:
        mainSizer.Add(wx.StaticText(Texture,label=' Insufficient number of patterns for usual texture analysis'))
        mainSizer.Add(wx.StaticText(Texture,label=
            ' For structural fits that need preferred orientation corrections, use terms in the Data tab instead.'))
    else:
        unique = []
        for i in h:
            setting = [h[i]['Sample Parameters'][key] for key in ('Omega', 'Chi', 'Phi', 'Azimuth')]
            if setting not in unique: unique.append(setting)
        if len(unique) == 1:
            mainSizer.Add(wx.StaticText(Texture,label=
            " All your histograms have the same values for Omega, Chi, Phi and Azimuth.\n Texture fitting requires differing values."))
        elif len(unique) == 2:
            mainSizer.Add(wx.StaticText(Texture,label=
                " You have only two unique settings for Omega, Chi, Phi and Azimuth.\n Texture fitting requires more differing values."))
    if textureData['Order'] or textureData['SH Coeff'][0]:
        mainSizer.Add(wx.StaticText(Texture,wx.ID_ANY,
        '\n *** To remove this texture: Turn off refinement and set Harmonic order to zero below.\n'))
    G2G.HorizontalLine(mainSizer,Texture)

    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl= f"Spherical harmonics texture data for {data['General']['Name']!r}"[:60]
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)

    mainSizer.Add(wx.StaticText(Texture,label=' Texture Index J = %7.3f'%(G2lat.textureIndex(textureData['SH Coeff'][1]))))
    mainSizer.Add((0,5),0)
    shSizer = wx.FlexGridSizer(0,6,5,5)
    shSizer.Add(wx.StaticText(Texture,-1,' Texture model: '),0,WACV)
    shModel = wx.ComboBox(Texture,value=textureData['Model'],choices=shModels,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    shModel.Bind(wx.EVT_COMBOBOX,OnShModel)
    shSizer.Add(shModel,0,WACV)
    shSizer.Add(wx.StaticText(Texture,-1,'  Harmonic order: '),0,WACV)
    shOrder = wx.ComboBox(Texture,value=str(textureData['Order']),choices=[str(2*i) for i in range(18)],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    shOrder.Bind(wx.EVT_COMBOBOX,OnShOrder)
    shSizer.Add(shOrder,0,WACV)
    shRef = wx.CheckBox(Texture,label=' Refine texture?')
    shRef.SetValue(textureData['SH Coeff'][0])
    shRef.Bind(wx.EVT_CHECKBOX, OnSHRefine)
    shSizer.Add(shRef,0,WACV)
    shShow = wx.CheckBox(Texture,label=' Show coeff.?')
    shShow.SetValue(textureData['SHShow'])
    shShow.Bind(wx.EVT_CHECKBOX, OnSHShow)
    shSizer.Add(shShow,0,WACV)
    mainSizer.Add(shSizer,0,0)
    mainSizer.Add((0,5),0)
    PTSizer = wx.FlexGridSizer(0,5,5,5)
    PTSizer.Add(wx.StaticText(Texture,-1,' Texture plot type: '),0,WACV)
    choices = ['Axial pole distribution','Pole figure','Inverse pole figure','3D pole distribution']
    pfType = wx.ComboBox(Texture,-1,value=str(textureData['PlotType']),choices=choices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    pfType.Bind(wx.EVT_COMBOBOX,OnPfType)
    PTSizer.Add(pfType,0,WACV)
    if 'Axial' not in textureData['PlotType'] and '3D' not in textureData['PlotType']:
        PTSizer.Add(wx.StaticText(Texture,-1,' Projection type: '),0,WACV)
        projSel = wx.ComboBox(Texture,-1,value=G2frame.Projection,choices=['equal area','stereographic'],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        projSel.Bind(wx.EVT_COMBOBOX,OnProjSel)
        PTSizer.Add(projSel,0,WACV)
        if textureData['PlotType'] == 'Pole figure' and len(textureData['det Angles']):
            shoDet = wx.CheckBox(Texture,-1,label=' Show detectors?')
            shoDet.SetValue(textureData['ShoDet'])
            shoDet.Bind(wx.EVT_CHECKBOX, OnShoDet)
            PTSizer.Add(shoDet,0,WACV)
        else:
            PTSizer.Add((0,5),0)
    elif '3D' in textureData['PlotType']:
        PTSizer.Add(wx.StaticText(Texture,-1,' Show projections for: '),0,WACV)
        proj = ['','x','y','z','xy','xz','yz','xyz']
        projType = wx.ComboBox(Texture,wx.ID_ANY,value=generalData['3Dproj'],choices=proj,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        projType.Bind(wx.EVT_COMBOBOX, OnProj)
        PTSizer.Add(projType,0,WACV)
        PTSizer.Add((0,5),0)

    if textureData['PlotType'] in ['Pole figure','Axial pole distribution','3D pole distribution']:
        PTSizer.Add(wx.StaticText(Texture,-1,' Pole figure HKL: '),0,WACV)
        PH = textureData['PFhkl']
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
        pfVal = wx.TextCtrl(Texture,-1,'%d %d %d'%(PH[0],PH[1],PH[2]),style=wx.TE_PROCESS_ENTER)
    else:
        PTSizer.Add(wx.StaticText(Texture,-1,' Inverse pole figure XYZ: '),0,WACV)
        PX = textureData['PFxyz']
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
        pfVal = wx.TextCtrl(Texture,-1,'%3.1f %3.1f %3.1f'%(PX[0],PX[1],PX[2]),style=wx.TE_PROCESS_ENTER)
    pfVal.Bind(wx.EVT_TEXT_ENTER,OnPFValue)
    pfVal.Bind(wx.EVT_KILL_FOCUS,OnPFValue)
    PTSizer.Add(pfVal,0,WACV)
    if 'Axial' not in textureData['PlotType'] and '3D' not in textureData['PlotType']:
        PTSizer.Add(wx.StaticText(Texture,-1,' Color scheme'),0,WACV)
        choice = [m for m in mpl.cm.datad.keys()]+['GSPaired','GSPaired_r',]       #if not m.endswith("_r")
        choice.sort()
        colorSel = wx.ComboBox(Texture,-1,value=G2frame.ContourColor,choices=choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        colorSel.Bind(wx.EVT_COMBOBOX,OnColorSel)
        PTSizer.Add(colorSel,0,WACV)
        if 'figure' in textureData['PlotType']:
            popLA = wx.Button(Texture,-1,"Make CSV file")
            popLA.Bind(wx.EVT_BUTTON, OnCSV)
            PTSizer.Add(popLA,0,WACV)
    mainSizer.Add(PTSizer,0)
    mainSizer.Add((0,5),0)
    if textureData['SHShow']:
        mainSizer.Add(wx.StaticText(Texture,-1,' Spherical harmonic coefficients: '))
        mainSizer.Add((0,5),0)
        ODFSizer = wx.FlexGridSizer(0,8,2,2)
        ODFkeys = list(textureData['SH Coeff'][1].keys())
        ODFkeys.sort()
        for item in ODFkeys:
            ODFSizer.Add(wx.StaticText(Texture,-1,item),0,WACV)
            ODFval = G2G.ValidatedTxtCtrl(Texture,textureData['SH Coeff'][1],item,nDig=(8,3),OnLeave=OnODFValue)
            ODFSizer.Add(ODFval,0,WACV)
        mainSizer.Add(ODFSizer,0)
        mainSizer.Add((0,5),0)
    mainSizer.Add((0,5),0)
    mainSizer.Add(wx.StaticText(Texture,-1,' Sample orientation angle zeros: '),0)
    mainSizer.Add((0,5),0)
    angSizer = wx.BoxSizer(wx.HORIZONTAL)
    angIndx = {}
    for item in ['Sample omega','Sample chi','Sample phi']:
        angRef = wx.CheckBox(Texture,-1,label=item+': ')
        angRef.SetValue(textureData[item][0])
        angIndx[angRef.GetId()] = item
        angRef.Bind(wx.EVT_CHECKBOX, OnAngRef)
        angSizer.Add(angRef,0,WACV)
        angVal = G2G.ValidatedTxtCtrl(Texture,textureData[item],1,nDig=(8,2))
        angSizer.Add(angVal,0,WACV|wx.LEFT,5)
    mainSizer.Add(angSizer,0,wx.LEFT,5)
#        mainSizer.Add(SHPenalty(textureData['Penalty']),0,wx.LEFT,5)  for future
    G2phsG.SetPhaseWindow(Texture,mainSizer)

#### UpdateWavesData GUI   
def UpdateWavesData(G2frame,data,Scroll=0):

    generalData = data['General']
    cx,ct,cs,cia = generalData['AtomPtrs']
    typeNames = {'Sfrac':' Site fraction','Spos':' Position','Sadp':' Thermal motion','Smag':' Magnetic moment'}
    numVals = {'Sfrac':2,'Spos':6,'Sadp':12,'Smag':6,'ZigZag':5,'Block':5,'Crenel':2}
    posNames = ['Xsin','Ysin','Zsin','Xcos','Ycos','Zcos','Tmin','Tmax','Xmax','Ymax','Zmax']
    adpNames = ['U11sin','U22sin','U33sin','U12sin','U13sin','U23sin',
        'U11cos','U22cos','U33cos','U12cos','U13cos','U23cos']
    magNames = ['MXsin','MYsin','MZsin','MXcos','MYcos','MZcos']
    fracNames = ['Fsin','Fcos','Fzero','Fwid']
    waveTypes = {'Sfrac':['Fourier','Crenel'],'Spos':['Fourier','ZigZag','Block',],'Sadp':['Fourier',],'Smag':['Fourier',]}
    Labels = {'Spos':posNames,'Sfrac':fracNames,'Sadp':adpNames,'Smag':magNames}
    Indx = {}
    waveData = G2frame.waveData
    G2frame.GetStatusBar().SetStatusText('',1)
    generalData = data['General']
    SGData = generalData['SGData']
    SSGData = generalData['SSGData']
    cx,ct,cs,cia = generalData['AtomPtrs']
    atomData = data['Atoms']
    D4Map = generalData.get('4DmapData',{'rho':[]})
    if waveData.GetSizer():
        waveData.GetSizer().Clear(True)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    topSizer.Add(wx.StaticText(waveData,label=' Incommensurate propagation wave data: Select atom to edit: '),0,WACV)
    atNames = []
    for atm in atomData:
        atNames.append(atm[ct-1])
    if not atNames:
        return
    if G2frame.atmSel not in atNames:
        G2frame.atmSel = atNames[0]

    def OnAtmSel(event):
        Obj = event.GetEventObject()
        G2frame.atmSel = Obj.GetValue()
        RepaintAtomInfo()

    def RepaintAtomInfo(Scroll=0):
        G2frame.bottomSizer.Clear(True)
        G2frame.bottomSizer = ShowAtomInfo()
        mainSizer.Add(G2frame.bottomSizer)
        mainSizer.Layout()
        G2frame.dataWindow.Refresh()
        waveData.SetVirtualSize(mainSizer.GetMinSize())
        waveData.Scroll(0,Scroll)
        G2frame.dataWindow.SendSizeEvent()

    def ShowAtomInfo():

        global mapSel       #so it can be seen below in OnWavePlot
        def AtomSizer(atom):
            global mapSel

            def OnShowWave(event):
                Obj = event.GetEventObject()
                atom = Indx[Obj.GetId()]
                Ax = Obj.GetValue()
                G2plt.ModulationPlot(G2frame,data,atom,Ax)

            atomSizer = wx.BoxSizer(wx.HORIZONTAL)
            atomSizer.Add(wx.StaticText(waveData,label=
            ' Modulation data for atom: %s  Site sym: %s'%(atom[0],atom[cs].strip())),0,WACV)
            axchoice = ['x','y','z']
            if len(D4Map['rho']):
                atomSizer.Add(wx.StaticText(waveData,label=' Show contour map for axis: '),0,WACV)
                mapSel = wx.ComboBox(waveData,value=' ',choices=axchoice,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                mapSel.Bind(wx.EVT_COMBOBOX,OnShowWave)
                Indx[mapSel.GetId()] = atom
                atomSizer.Add(mapSel,0,WACV)
            return atomSizer

        def WaveSizer(iatm,wavedata,Stype,typeName,Names):

            def OnWaveType(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                if len(atm[-1]['SS1'][item]) <= 1:
                    atm[-1]['SS1'][item] = [0,]
                    atm[-1]['SS1'][item][0] = waveType.GetValue()
                    wx.CallAfter(RepaintAtomInfo,G2frame.waveData.GetScrollPos(wx.VERTICAL))
                else:
                    if len(waveTypes[Stype]) > 1:
                        waveType.SetValue(atm[-1]['SS1'][Stype][0])
                        G2G.G2MessageBox(G2frame,'Warning: can only change wave type if no waves','Not changed')

            def OnAddWave(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                nt = numVals[Stype]
                if not len(atm[-1]['SS1'][item]):
                    if waveTyp in ['ZigZag','Block','Crenel']:
                        nt = numVals[waveTyp]
                    atm[-1]['SS1'][item] = [0,]
                    atm[-1]['SS1'][item][0] = waveType.GetValue()
                atm[-1]['SS1'][item].append([[0.0 for i in range(nt)],False])
                wx.CallAfter(RepaintAtomInfo,G2frame.waveData.GetScrollPos(wx.VERTICAL))

            def OnRefWave(event):
                Obj = event.GetEventObject()
                item,iwave = Indx[Obj.GetId()]
                atm[-1]['SS1'][item][iwave+1][1] = not atm[-1]['SS1'][item][iwave+1][1]

            def OnDelWave(event):
                Obj = event.GetEventObject()
                item,iwave = Indx[Obj.GetId()]
                del atm[-1]['SS1'][item][iwave+1]
                if len(atm[-1]['SS1'][item]) == 1:
                    atm[-1]['SS1'][item][0] = 'Fourier'
                wx.CallAfter(RepaintAtomInfo,G2frame.waveData.GetScrollPos(wx.VERTICAL))

            def OnWavePlot(invalid,value,tc):
                if len(D4Map['rho']):
                    Ax = mapSel.GetValue()
                    if Ax:
                        G2plt.ModulationPlot(G2frame,data,atm,Ax)

            waveTyp,waveBlk = 'Fourier',[]
            if len(wavedata):
                waveTyp = wavedata[0]
                waveBlk = wavedata[1:]
            waveSizer = wx.BoxSizer(wx.VERTICAL)
            waveHead = wx.BoxSizer(wx.HORIZONTAL)
            waveHead.Add(wx.StaticText(waveData,label=typeName+' modulation parameters: '),0,WACV)
            if Stype == 'Smag' and len(waveBlk):    #only allow one magnetic wave - keeps it simple for now
                pass
            else:
                waveAdd = wx.Button(waveData,label='Add wave?')
                waveAdd.Bind(wx.EVT_BUTTON, OnAddWave)
                Indx[waveAdd.GetId()] = Stype
                waveHead.Add(waveAdd,0,WACV)
                waveHead.Add(wx.StaticText(waveData,label='   WaveType: '),0,WACV)
                waveType = wx.ComboBox(waveData,value=waveTyp,choices=waveTypes[Stype],
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[waveType.GetId()] = Stype
                waveType.Bind(wx.EVT_COMBOBOX,OnWaveType)
                waveHead.Add(waveType,0,WACV)
            waveSizer.Add(waveHead)
            if len(waveBlk):
                nx = 0
                for iwave,wave in enumerate(waveBlk):
                    if not iwave:
                        if waveTyp in ['ZigZag','Block','Crenel']:
                            nx = 1
                        CSI = G2spc.GetSSfxuinel(waveTyp,Stype,1,xyz,SGData,SSGData)[0]
                    else:
                        CSI = G2spc.GetSSfxuinel('Fourier',Stype,iwave+1-nx,xyz,SGData,SSGData)[0]
                    waveName = 'Fourier'
                    if Stype == 'Sfrac':
                        nterm = 2
                        if 'Crenel' in waveTyp and not iwave:
                            waveName = 'Crenel'
                            names = Names[2:]
                        else:
                            names = Names[:2]
                        Waves = wx.FlexGridSizer(0,4,5,5)
                    elif Stype == 'Spos':
                        nterm = 6
                        if waveTyp in ['ZigZag','Block'] and not iwave:
                            nterm = 5
                            names = Names[6:]
                            Waves = wx.FlexGridSizer(0,7,5,5)
                            waveName = waveTyp
                        else:
                            names = Names[:6]
                            Waves = wx.FlexGridSizer(0,8,5,5)
                    elif Stype == 'Sadp':
                        nterm = 12
                        names = Names
                        Waves = wx.FlexGridSizer(0,8,5,5)
                    elif Stype == 'Smag':
                        nterm = 6
                        names = Names
                        Waves = wx.FlexGridSizer(0,8,5,5)
                    waveSizer.Add(wx.StaticText(waveData,label=' %s  parameters: %s'%(waveName,str(names).rstrip(']').lstrip('[').replace("'",''))),0)
                    for ival in range(nterm):
                        val = wave[0][ival]
                        if np.any(CSI[0][ival]):
                            minmax = [-0.2,0.2]
                            if Stype == 'Smag':
                                minmax = [-20.,20.]
                            if waveTyp in ['ZigZag','Block','Crenel'] and not iwave and ival < 2:
                                if not ival:
                                    minmax = [0.,2.]
                                else:
                                    minmax = [wave[0][0],1.0+wave[0][0]]
                            waveVal = G2G.ValidatedTxtCtrl(waveData,wave[0],ival,nDig=(10,5),xmin=minmax[0],xmax=minmax[1],OnLeave=OnWavePlot)
                        else:
                            waveVal = G2G.ReadOnlyTextCtrl(waveData,value='%.5f'%(val))
                        Waves.Add(waveVal,0,WACV)
                        if len(wave[0]) > 6 and ival == 5:
                            Waves.Add((5,5),0)
                            Waves.Add((5,5),0)
                    waveRef = wx.CheckBox(waveData,label='Refine?')
                    waveRef.SetValue(wave[1])
                    Indx[waveRef.GetId()] = [Stype,iwave]
                    waveRef.Bind(wx.EVT_CHECKBOX, OnRefWave)
                    Waves.Add(waveRef,0,WACV)
                    if iwave < len(waveBlk)-1:
                        Waves.Add((5,5),0)
                    else:
                        waveDel = wx.Button(waveData,wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
                        Indx[waveDel.GetId()] = [Stype,iwave]
                        waveDel.Bind(wx.EVT_BUTTON,OnDelWave)
                        Waves.Add(waveDel,0,WACV)
                    waveSizer.Add(Waves)
            return waveSizer

        iatm = atNames.index(G2frame.atmSel)
        atm = atomData[iatm]
        xyz = atm[cx:cx+3]
        atomSizer = wx.BoxSizer(wx.VERTICAL)
        G2G.HorizontalLine(atomSizer,waveData)
        atomSizer.Add(AtomSizer(atm))
        for Stype in ['Sfrac','Spos','Sadp','Smag']:
            if atm[cia] != 'A' and Stype == 'Sadp':    #Uiso can't have modulations! (why not?)
                continue
            if generalData['Type'] != 'magnetic' and Stype == 'Smag':
                break

            atomSizer.Add(WaveSizer(iatm,atm[-1]['SS1'][Stype],Stype,typeNames[Stype],Labels[Stype]))
        return atomSizer

    atms = wx.ComboBox(waveData,value=G2frame.atmSel,choices=atNames,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    atms.Bind(wx.EVT_COMBOBOX,OnAtmSel)
    topSizer.Add(atms,0)
    mainSizer.Add(topSizer,0)
    G2frame.bottomSizer = ShowAtomInfo()
    mainSizer.Add(G2frame.bottomSizer)
    G2phsG.SetPhaseWindow(G2frame.waveData,mainSizer,Scroll=Scroll)

