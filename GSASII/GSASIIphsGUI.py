# -*- coding: utf-8 -*-
#GSASII - phase data display routines
'''
Main routine here is :func:`UpdatePhaseData`, which displays the phase information
(called from :func:`GSASIIdataGUI:SelectDataTreeItem`).

Other top-level routines are:
:func:`GetSpGrpfromUser` (called locally only);
:func:`FindBondsDraw` and :func:`FindBondsDrawCell` (called locally and in GSASIIplot);
:func:`SetPhaseWindow` (called locally and in GSASIIddataGUI and GSASIIrestrGUI, multiple locations)
to control scrolling.

Routines for Phase dataframes follow.
'''
from __future__ import division, print_function
import platform
import os
import wx
import wx.grid as wg
import wx.lib.scrolledpanel as wxscroll
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
from . import GSASIIphsGUI2 as G2phsG2
from . import GSASIIrmcGUI as G2rmcG
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
from . import GSASIIconstrGUI as G2cnstG
from . import atmdata
from . import ISODISTORT as ISO
from . import SUBGROUPS
from pathlib import Path

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

#### phase class definitions ################################################################################
class SymOpDialog(wx.Dialog):
    '''Class to select a symmetry operator
    '''
    def __init__(self,parent,SGData,New=True,ForceUnit=False):
        wx.Dialog.__init__(self,parent,-1,'Select symmetry operator',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        self.SGData = SGData
        self.New = New
        self.Force = ForceUnit
        self.OpSelected = [0,0,0,[0,0,0],False,False]
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if ForceUnit:
            choice = ['No','Yes']
            self.force = wx.RadioBox(panel,-1,'Force to unit cell?',choices=choice)
            self.force.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.force,0,wx.TOP,5)
#        if SGData['SGInv']:
        choice = ['No','Yes']
        self.inv = wx.RadioBox(panel,-1,'Choose inversion?',choices=choice)
        self.inv.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.inv,0)
        if SGData['SGLatt'] != 'P':
            LattOp = G2spc.Latt2text(SGData['SGCen']).split(';')
            self.latt = wx.RadioBox(panel,-1,'Choose cell centering?',choices=LattOp)
            self.latt.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.latt,0)
        if SGData['SGLaue'] in ['-1','2/m','mmm','4/m','4/mmm']:
            Ncol = 2
        else:
            Ncol = 3
        OpList = []
        for Opr in SGData['SGOps']:
            OpList.append(G2spc.MT2text(Opr))
        self.oprs = wx.RadioBox(panel,-1,'Choose space group operator?',choices=OpList,
            majorDimension=Ncol)
        self.oprs.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.oprs,0,wx.BOTTOM,5)
        mainSizer.Add(wx.StaticText(panel,-1,"   Choose unit cell?"),0)
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellName = ['X','Y','Z']
        self.cell = []
        for i in range(3):
            self.cell.append(wx.SpinCtrl(panel,-1,cellName[i],size=wx.Size(50,20)))
            self.cell[-1].SetRange(-3,3)
            self.cell[-1].SetValue(0)
            self.cell[-1].Bind(wx.EVT_SPINCTRL, self.OnOpSelect)
            cellSizer.Add(self.cell[-1],0)
        mainSizer.Add(cellSizer,0,wx.BOTTOM,5)
        if self.New:
            choice = ['No','Yes']
            self.new = wx.RadioBox(panel,-1,'Generate new positions?',choices=choice)
            self.new.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.new,0)

        OkBtn = wx.Button(panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()
        self.Fit()

    def OnOpSelect(self,event):
        self.OpSelected[0] = self.inv.GetSelection()
        if self.SGData['SGLatt'] != 'P':
            self.OpSelected[1] = self.latt.GetSelection()
        self.OpSelected[2] = self.oprs.GetSelection()
        for i in range(3):
            self.OpSelected[3][i] = float(self.cell[i].GetValue())
        if self.New:
            self.OpSelected[4] = self.new.GetSelection()
        if self.Force:
            self.OpSelected[5] = self.force.GetSelection()

    def GetSelection(self):
        return self.OpSelected

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)
#==============================================================================
class SphereEnclosure(wx.Dialog):
    ''' Add atoms within sphere of enclosure to drawing

    :param wx.Frame parent: reference to parent frame (or None)
    :param general: general data (includes drawing data)
    :param atoms: drawing atoms data
    :param indx: list of selected atoms (may be empty)

    '''
    def __init__(self,parent,general,drawing,indx):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Supply sphere info',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.General = general
        self.Drawing = drawing
        self.indx = indx
        self.Sphere = [3.0,]
        self.centers = []
        self.atomTypes = [[item,False] for item in self.General['AtomTypes']]
        self.CenterOnParent()
        self.Draw()

    def Draw(self):

        def OnAtomType(event):
            Obj = event.GetEventObject()
            Id = Ind[Obj.GetId()]
            self.atomTypes[Id][1] = Obj.GetValue()

        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label=' Sphere of enclosure controls:'),0)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        atoms = []
        if self.indx is None:
            pass
        elif len(self.indx):
            topSizer.Add(wx.StaticText(self.panel,label=' Sphere centered at atoms: '),0,WACV)
            cx,ct,cs = self.Drawing['atomPtrs'][:3]
            for Id in self.indx:
                if Id < len(self.Drawing['Atoms']):
                    atom = self.Drawing['Atoms'][Id]
                    self.centers.append(atom[cx:cx+3])
                    atoms.append('%s(%s)'%(atom[ct-1],atom[cs-1]))
                else:
                    self.centers.append(list(self.Drawing['viewPoint'][0]))
                    atoms.append('View point')
            topSizer.Add(wx.ComboBox(self.panel,choices=atoms,value=atoms[0],
                style=wx.CB_READONLY|wx.CB_DROPDOWN),0,WACV)
        else:
            topSizer.Add(wx.StaticText(self.panel,label=' Sphere centered at drawing view point'),0,WACV)
            self.centers.append(self.Drawing['viewPoint'][0])
        mainSizer.Add(topSizer,0)
        sphereSizer = wx.BoxSizer(wx.HORIZONTAL)
        sphereSizer.Add(wx.StaticText(self.panel,label=' Sphere radius: '),0,WACV)
        radius = G2G.ValidatedTxtCtrl(self.panel,self.Sphere,0,nDig=(10,3),size=(65,25))
        sphereSizer.Add(radius,0,WACV)
        mainSizer.Add(sphereSizer,0)
        mainSizer.Add(wx.StaticText(self.panel,label=' Target selected atoms:'),0)
        atSizer = wx.BoxSizer(wx.HORIZONTAL)
        Ind = {}
        for i,item in enumerate(self.atomTypes):
            atm = wx.CheckBox(self.panel,label=item[0])
            atm.SetValue(item[1])
            atm.Bind(wx.EVT_CHECKBOX, OnAtomType)
            Ind[atm.GetId()] = i
            atSizer.Add(atm,0,WACV)
        mainSizer.Add(atSizer,0)

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
        used = []
        for atm in self.atomTypes:
            if atm[1]:
                used.append(str(atm[0]))
        return self.centers,self.Sphere[0],used

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

#==============================================================================
class TransformDialog(wx.Dialog):
    ''' Phase transformation X' = M*(X-U)+V

    :param wx.Frame parent: reference to parent frame (or None)
    :param phase: parent phase data

    #NB: commonNames & commonTrans defined in GSASIIdataGUI = G2gd
    '''
    def __init__(self,parent,phase,Trans=np.eye(3),Uvec=np.zeros(3),Vvec=np.zeros(3),ifMag=False,BNSlatt=''):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Setup phase transformation',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.Phase = copy.deepcopy(phase)   #will be a new phase!
#        self.Super = phase['General']['Super']
#        if self.Super:
#            self.Trans = np.eye(4)
#            self.Vec = np.zeros(4)
#        else:
        self.Trans = Trans
        self.Uvec = Uvec
        self.Vvec = Vvec
        self.oldSpGrp = copy.deepcopy(phase['General']['SGData']['SpGrp'])
        self.oldSGdata = copy.deepcopy(phase['General']['SGData'])
        self.newSpGrp = self.Phase['General']['SGData']['SpGrp']
        self.SGData = G2spc.SpcGroup(self.newSpGrp)[1]
        self.oldCell = copy.deepcopy(phase['General']['Cell'][1:8])
        self.newCell = self.Phase['General']['Cell'][1:8]
        self.Common = 'abc'
        self.ifMag = ifMag
        if ifMag:
            self.BNSlatt = BNSlatt
        self.ifConstr = False
        self.Mtrans = False
        self.kvec = [0.,0.,0.]
        self.Draw()
        self.CenterOnParent()

    def Draw(self):

        def OnCommon(event):
            Obj = event.GetEventObject()
            self.Common = Obj.GetValue()
            self.Mtrans = False
            if '*' in self.Common:
                A,B = G2lat.cell2AB(self.oldCell[:6])
                self.newCell[2:5] = [A[2,2],90.,90.]
                a,b = G2lat.cell2AB(self.newCell[:6])
                self.Trans = np.inner(a,B)    #correct!
                self.ifConstr = False
                self.newSpGrp = 'P 1'
                SGErr,SGData = G2spc.SpcGroup(self.newSpGrp)
                self.Phase['General']['SGData'] = SGData
            else:
                if self.Common == G2gd.commonNames[-1]:      #change setting
                    self.Vvec = G2spc.spg2origins[self.oldSpGrp]
                    self.newSpGrp = self.oldSpGrp
                else:
                    self.Trans = G2gd.commonTrans[self.Common]
                    if 'R' == self.Common[-1]:
                        self.newSpGrp += ' r'
                        SGErr,SGData = G2spc.SpcGroup(self.newSpGrp)
                        self.Phase['General']['SGData'] = SGData
                        SGTxt.SetLabel(self.newSpGrp)
            OnTest(event)

        def OnSpaceGroup(event):
            event.Skip()
            SpcGp = GetSpGrpfromUser(self.panel,self.newSpGrp)
            if SpcGp == self.newSpGrp or SpcGp is None: #didn't change it!
                return
            # try a lookup on the user-supplied name
            SpGrpNorm = G2spc.StandardizeSpcName(SpcGp)
            if SpGrpNorm:
                SGErr,self.SGData = G2spc.SpcGroup(SpGrpNorm)
            else:
                SGErr,self.SGData = G2spc.SpcGroup(SpcGp)
            if SGErr:
                text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                SGTxt.SetLabel(self.newSpGrp)
                msg = 'Space Group Error'
                Text = '\n'.join(text)
                wx.MessageBox(Text,caption=msg,style=wx.ICON_EXCLAMATION)
            else:
                text,table = G2spc.SGPrint(self.SGData)
                self.Phase['General']['SGData'] = self.SGData
                self.newSpGrp = SpcGp
                SGTxt.SetLabel(self.Phase['General']['SGData']['SpGrp'])
                msg = 'Space Group Information'
                G2G.SGMessageBox(self.panel,msg,text,table).Show()
            if self.ifMag:
                self.BNSlatt = self.SGData['SGLatt']
                G2spc.SetMagnetic(self.SGData)
            if self.Phase['General']['Type'] == 'magnetic':
                Nops = len(self.SGData['SGOps'])*len(self.SGData['SGCen'])
                if self.SGData['SGInv']:
                    Nops *= 2
                self.SGData['SpnFlp'] = Nops*[1,]
                del self.oldSGdata['MAXMAGN']
            wx.CallAfter(self.Draw)

        def OnShowOps(event):
            text,table = G2spc.SGPrint(self.SGData,AddInv=True)
            if self.ifMag:
                msg = 'Magnetic space group information'
                OprNames,SpnFlp = G2spc.GenMagOps(self.SGData)
                text[0] = ' Magnetic Space Group: '+self.SGData['MagSpGrp']
                text[3] = ' The magnetic lattice point group is '+self.SGData['MagPtGp']
                G2G.SGMagSpinBox(self.panel,msg,text,table,self.SGData['SGCen'],OprNames,
                    self.SGData['SpnFlp'],False).Show()
            else:
                msg = 'Space group information'
                G2G.SGMessageBox(self.panel,msg,text,table).Show()

        def OnTest(event):
            if not self.TestMat():
                return
            if self.Mtrans:
                self.newCell = G2lat.TransformCell(self.oldCell[:6],self.Trans.T)
            else:
                self.newCell = G2lat.TransformCell(self.oldCell[:6],self.Trans)
            wx.CallAfter(self.Draw)

        def OnMag(event):
            self.ifMag = True
            self.BNSlatt = self.SGData['SGLatt']
            G2spc.SetMagnetic(self.SGData)
            wx.CallAfter(self.Draw)

        def OnConstr(event):
            self.ifConstr = constr.GetValue()

        def OnBNSlatt(event):
            Obj = event.GetEventObject()
            self.BNSlatt = Obj.GetValue()
            if self.BNSlatt == self.SGData['SGLatt']:
                return
            GenSym,GenFlg,BNSsym = G2spc.GetGenSym(self.SGData)
            self.SGData['BNSlattsym'] = [self.BNSlatt,BNSsym[self.BNSlatt]]
            self.SGData['SGSpin'] = [1,]*len(self.SGData['SGSpin'])
            wx.CallAfter(self.Draw)

        def OnMtrans(event):
            Obj = event.GetEventObject()
            self.Mtrans = Obj.GetValue()

        def OnSpinOp(event):
            Obj = event.GetEventObject()
            isym = Indx[Obj.GetId()]+1
            spCode = {'red':-1,'black':1}
            self.SGData['SGSpin'][isym] = spCode[Obj.GetValue()]
            G2spc.CheckSpin(isym,self.SGData)
            G2spc.SetMagnetic(self.SGData)
            wx.CallAfter(self.Draw)

        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if self.ifMag:
            if self.BNSlatt != self.SGData['SGLatt']:
                GenSym,GenFlg,BNSsym = G2spc.GetGenSym(self.SGData)
                self.SGData['BNSlattsym'] = [self.BNSlatt,BNSsym[self.BNSlatt]]
        else:
            mag = wx.Button(self.panel,label='Make new phase magnetic?')
            mag.Bind(wx.EVT_BUTTON,OnMag)
            mainSizer.Add(mag,0)
        MatSizer = wx.BoxSizer(wx.HORIZONTAL)
        transSizer = wx.BoxSizer(wx.VERTICAL)
        transSizer.Add((5,5),0)
        transSizer.Add(wx.StaticText(self.panel,label=
            " Cell transformation via g'=gM; g=metric tensor \n XYZ transformation via M*(X-U)+V = X'; M* = inv(M)"))
#        if self.Super:
#            Trmat = wx.FlexGridSizer(4,4,0,0)
#        else:
        commonSizer = wx.BoxSizer(wx.HORIZONTAL)
        commonSizer.Add(wx.StaticText(self.panel,label=' Common transformations: '),0,WACV)
        if self.oldSpGrp not in G2spc.spg2origins:
            common = wx.ComboBox(self.panel,value=self.Common,choices=G2gd.commonNames[:-1],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        else:
            common = wx.ComboBox(self.panel,value=self.Common,choices=G2gd.commonNames,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
        common.Bind(wx.EVT_COMBOBOX,OnCommon)
        commonSizer.Add(common,0,WACV)
        transSizer.Add(commonSizer)
        transSizer.Add(G2G.XformMatrix(self.panel,self.Trans,self.Uvec,self.Vvec))
        MatSizer.Add((10,0),0)
        MatSizer.Add(transSizer)
        mainSizer.Add(MatSizer)
        if self.ifMag:
            MagSizer = wx.BoxSizer(wx.HORIZONTAL)
            if not self.oldSGdata.get('MAXMAGN',[]):
                Mtrans = wx.CheckBox(self.panel,label=' Use matrix transform?')
                Mtrans.SetValue(self.Mtrans)
                Mtrans.Bind(wx.EVT_CHECKBOX,OnMtrans)
                MagSizer.Add(Mtrans,0,WACV)
            mainSizer.Add(MagSizer,0)
        mainSizer.Add(wx.StaticText(self.panel,label=' Old lattice parameters:'),0)
        mainSizer.Add(wx.StaticText(self.panel,label=
            ' a = %.5f       b = %.5f      c = %.5f'%(self.oldCell[0],self.oldCell[1],self.oldCell[2])),0)
        mainSizer.Add(wx.StaticText(self.panel,label=' alpha = %.3f beta = %.3f gamma = %.3f'%
            (self.oldCell[3],self.oldCell[4],self.oldCell[5])),0)
        mainSizer.Add(wx.StaticText(self.panel,label=' volume = %.3f'%(self.oldCell[6])),0)
        mainSizer.Add(wx.StaticText(self.panel,label=' New lattice parameters:'),0)
        mainSizer.Add(wx.StaticText(self.panel,label=
            ' a = %.5f       b = %.5f      c = %.5f'%(self.newCell[0],self.newCell[1],self.newCell[2])),0)
        mainSizer.Add(wx.StaticText(self.panel,label=' alpha = %.3f beta = %.3f gamma = %.3f'%
            (self.newCell[3],self.newCell[4],self.newCell[5])),0)
        mainSizer.Add(wx.StaticText(self.panel,label=' volume = %.3f'%(self.newCell[6])),0)
        sgSizer = wx.BoxSizer(wx.HORIZONTAL)
        sgSizer.Add(wx.StaticText(self.panel,label=' Target space group: '),0,WACV)
        SGTxt = wx.Button(self.panel,wx.ID_ANY,self.newSpGrp,size=(100,-1))
        SGTxt.Bind(wx.EVT_BUTTON,OnSpaceGroup)
        sgSizer.Add(SGTxt,0,WACV)
        showOps = wx.Button(self.panel,label=' Show operators?')
        showOps.Bind(wx.EVT_BUTTON,OnShowOps)
        sgSizer.Add(showOps,0,WACV)
        mainSizer.Add(sgSizer,0)
        if 'magnetic' not in self.Phase['General']['Type']:
            if self.ifMag:
                Indx = {}
                GenSym,GenFlg,BNSsym = G2spc.GetGenSym(self.SGData)
                BNSizer = wx.BoxSizer(wx.HORIZONTAL)
                BNSizer.Add(wx.StaticText(self.panel,label=' Select BNS lattice:'),0,WACV)
                BNSkeys = [self.SGData['SGLatt'],]+list(BNSsym.keys())
                BNSkeys.sort()
                try:        #this is an ugly kluge - bug in wx.ComboBox
                    if self.BNSlatt[2] in ['a','b','c']:
                        BNSkeys.reverse()
                except:
                    pass
                BNS = wx.ComboBox(self.panel,choices=BNSkeys,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                BNS.SetValue(self.BNSlatt)
                BNS.Bind(wx.EVT_COMBOBOX,OnBNSlatt)
                BNSizer.Add(BNS,0,WACV)

                spinColor = ['black','red']
                spCode = {-1:'red',1:'black'}
                for isym,sym in enumerate(GenSym[1:]):
                    BNSizer.Add(wx.StaticText(self.panel,label=' %s: '%(sym.strip())),0,WACV)
                    spinOp = wx.ComboBox(self.panel,value=spCode[self.SGData['SGSpin'][isym+1]],choices=spinColor,
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    Indx[spinOp.GetId()] = isym
                    spinOp.Bind(wx.EVT_COMBOBOX,OnSpinOp)
                    BNSizer.Add(spinOp,0,WACV)
                OprNames,SpnFlp = G2spc.GenMagOps(self.SGData)
                self.SGData['SpnFlp'] = SpnFlp
                mainSizer.Add(BNSizer,0)
                mainSizer.Add(wx.StaticText(self.panel,label=' Magnetic Space Group: '+self.SGData['MagSpGrp']),0)
            if self.ifMag:
                mainSizer.Add(wx.StaticText(self.panel, \
                    label=' NB: Nonmagnetic atoms will be deleted from new phase'),0)
            constr = wx.CheckBox(self.panel,label=' Make constraints between phases?')
            constr.SetValue(self.ifConstr)
            constr.Bind(wx.EVT_CHECKBOX,OnConstr)
            mainSizer.Add(constr,0)
        TestBtn = wx.Button(self.panel,-1,"Test")
        TestBtn.Bind(wx.EVT_BUTTON, OnTest)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(self.panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(TestBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def TestMat(self):
        VC = nl.det(self.Trans)
        if VC < 0.:
            wx.MessageBox('Warning - left handed transformation',caption='Transformation matrix check',
                style=wx.ICON_EXCLAMATION)
            return True
        try:
            nl.inv(self.Trans)
        except nl.LinAlgError:
            wx.MessageBox('ERROR - bad transformation matrix',caption='Transformation matrix check',
                style=wx.ICON_ERROR)
            return False
        return True

    def GetSelection(self):
        self.Phase['General']['SGData'] = self.SGData
        if self.ifMag:
            self.Phase['General']['Name'] += ' mag: '
        else:
            self.Phase['General']['Name'] += ' %s'%(self.Common)
        if not self.TestMat():
            return None
        if self.Mtrans:
            self.Phase['General']['Cell'][1:] = G2lat.TransformCell(self.oldCell[:6],self.Trans.T)
            return self.Phase,self.Trans.T,self.Uvec,self.Vvec,self.ifMag,self.ifConstr,self.Common
        else:
            self.Phase['General']['Cell'][1:] = G2lat.TransformCell(self.oldCell[:6],self.Trans)
            return self.Phase,self.Trans,self.Uvec,self.Vvec,self.ifMag,self.ifConstr,self.Common

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

#==============================================================================
class UseMagAtomDialog(wx.Dialog):
    '''Get user selected magnetic atoms after cell transformation
    '''
    def __init__(self,parent,Name,Atoms,atCodes,atMxyz,ifMag=True,ifOK=False,ifDelete=False):
        title = 'Subgroup atom list'
        if ifMag:
            title = 'Magnetic atom selection'
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,
                            pos=wx.DefaultPosition,size=(450,275),
                            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.panel = wxscroll.ScrolledPanel(self)         #just a dummy - gets destroyed in Draw!
#        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.Name = Name
        self.Atoms = Atoms
        self.atCodes = atCodes
        self.atMxyz = atMxyz
        self.ifMag = ifMag
        self.ifOK = ifOK
        self.ifDelete = ifDelete
        self.Use = len(self.Atoms)*[True,]
        self.Draw()

    def Draw(self):

        def OnUseChk(event):
            Obj = event.GetEventObject()
            iuse = Indx[Obj.GetId()]
            self.Use[iuse] = not self.Use[iuse]
            Obj.SetValue(self.Use[iuse])

        self.panel.Destroy()
        self.panel = wxscroll.ScrolledPanel(self,style = wx.DEFAULT_DIALOG_STYLE)
        Indx = {}
        Mstr = [' Mx',' My',' Mz']
        Xstr = ['X','Y','Z']
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label='For: %s'%self.Name),0)

        if self.ifMag:
            mainSizer.Add(wx.StaticText(self.panel,label='        Name, x, y, z, allowed moments, mag. site sym:'),0)
        else:
            mainSizer.Add(wx.StaticText(self.panel,label='        Name, x, y, z, allowed xyz, site sym:'),0)
        atmSizer = wx.FlexGridSizer(0,2,5,5)
        for iuse,[use,atom,mxyz] in enumerate(zip(self.Use,self.Atoms,self.atMxyz)):
            mstr = [' ---',' ---',' ---']
            for i,mx in enumerate(mxyz[1]):
                if mx:
                    if self.ifMag:
                        mstr[i] = Mstr[i]
                    else:
                        mstr[i] = Xstr[i]
            if self.ifMag:
                useChk = wx.CheckBox(self.panel,label='Use?')
                Indx[useChk.GetId()] = iuse
                useChk.SetValue(use)
                useChk.Bind(wx.EVT_CHECKBOX, OnUseChk)
                atmSizer.Add(useChk,0,WACV)
            else:
                atmSizer.Add((2,2),0)
            text = '  %5s %10.5f %10.5f %10.5f (%s,%s,%s) %s   '%(atom[0],atom[3],atom[4],atom[5],mstr[0],mstr[1],mstr[2],mxyz[0])
            atmSizer.Add(wx.StaticText(self.panel,label=text),0,WACV)
        mainSizer.Add(atmSizer)

        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        if self.ifOK:
            OKBtn = wx.Button(self.panel,-1,"OK")
            OKBtn.Bind(wx.EVT_BUTTON, self.OnNo)
            btnSizer.Add(OKBtn)
        else:
            YesBtn = wx.Button(self.panel,-1,"Yes")
            YesBtn.Bind(wx.EVT_BUTTON, self.OnYes)
            NoBtn = wx.Button(self.panel,-1,"No")
            NoBtn.Bind(wx.EVT_BUTTON, self.OnNo)
            btnSizer.Add((20,20),1)
            btnSizer.Add(YesBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(NoBtn)
            if self.ifDelete:
                DeleteBtn = wx.Button(self.panel,-1,"Delete")
                DeleteBtn.Bind(wx.EVT_BUTTON, self.OnDelete)
                btnSizer.Add((20,20),1)
                btnSizer.Add(DeleteBtn)
            btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        size = np.array(self.GetSize())
        self.panel.SetupScrolling()
        self.panel.SetAutoLayout(1)
        size = [size[0]-5,size[1]-20]       #this fiddling is needed for older wx!
        self.panel.SetSize(size)

    def GetSelection(self):
        useAtoms = []
        useatCodes = []
        for use,atom,code in zip(self.Use,self.Atoms,self.atCodes):
            if use:
                useAtoms.append(atom)
                useatCodes.append(code)
        return useAtoms,useatCodes

    def OnYes(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_YES)

    def OnNo(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_NO)

    def OnDelete(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_DELETE)


#==============================================================================
class RotationDialog(wx.Dialog):
    ''' Get Rotate & translate matrix & vector - currently not used
    needs rethinking - possible use to rotate a group of atoms about some
    vector/origin + translation

    '''
    def __init__(self,parent):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Atom group rotation/translation',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.Trans = np.eye(3)
        self.Vec = np.zeros(3)
        self.rotAngle = 0.
        self.rotVec = np.array([0.,0.,1.])
        self.Expand = ''
        self.Draw()

    def Draw(self):

        def OnExpand(event):
            self.Expand = expand.GetValue()

        def OnRotAngle(event):
            event.Skip()
            self.rotAngle = float(rotangle.GetValue())
            rotangle.SetValue('%5.3f'%(self.rotAngle))
            Q = G2mth.AVdeg2Q(self.rotAngle,self.rotVec)
            self.Trans = G2mth.Q2Mat(Q)
            self.Draw()

        def OnRotVec(event):
            event.Skip()
            vals = rotvec.GetValue()
            vals = vals.split()
            self.rotVec = np.array([float(val) for val in vals])
            rotvec.SetValue('%5.3f %5.3f %5.3f'%(self.rotVec[0],self.rotVec[1],self.rotVec[2]))
            Q = G2mth.AVdeg2Q(self.rotAngle,self.rotVec)
            self.Trans = G2mth.Q2Mat(Q)
            self.Draw()

        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        MatSizer = wx.BoxSizer(wx.HORIZONTAL)
        transSizer = wx.BoxSizer(wx.VERTICAL)
        transSizer.Add(wx.StaticText(self.panel,label=" XYZ Transformation matrix && vector: "+ \
            "\n B*M*A*(X-V)+V = X'\n A,B: Cartesian transformation matrices"))
        Trmat = wx.FlexGridSizer(3,5,0,0)
        for iy,line in enumerate(self.Trans):
            for ix,val in enumerate(line):
                item = G2G.ValidatedTxtCtrl(self.panel,self.Trans[iy],ix,nDig=(10,3),size=(65,25))
                Trmat.Add(item)
            Trmat.Add((25,0),0)
            vec = G2G.ValidatedTxtCtrl(self.panel,self.Vec,iy,nDig=(10,3),size=(65,25))
            Trmat.Add(vec)
        transSizer.Add(Trmat)
        MatSizer.Add((10,0),0)
        MatSizer.Add(transSizer)
        mainSizer.Add(MatSizer)
        rotationBox = wx.BoxSizer(wx.HORIZONTAL)
        rotationBox.Add(wx.StaticText(self.panel,label=' Rotation angle: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
        rotangle = wx.TextCtrl(self.panel,value='%5.3f'%(self.rotAngle),
            size=(50,25),style=wx.TE_PROCESS_ENTER)
        rotangle.Bind(wx.EVT_TEXT_ENTER,OnRotAngle)
        rotangle.Bind(wx.EVT_KILL_FOCUS,OnRotAngle)
        rotationBox.Add(rotangle,0,WACV)
        rotationBox.Add(wx.StaticText(self.panel,label=' about vector: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
        rotvec = wx.TextCtrl(self.panel,value='%5.3f %5.3f %5.3f'%(self.rotVec[0],self.rotVec[1],self.rotVec[2]),
            size=(100,25),style=wx.TE_PROCESS_ENTER)
        rotvec.Bind(wx.EVT_TEXT_ENTER,OnRotVec)
        rotvec.Bind(wx.EVT_KILL_FOCUS,OnRotVec)
        rotationBox.Add(rotvec,0,WACV)
        mainSizer.Add(rotationBox,0)
        expandChoice = ['','xy','xz','yz','xyz']
        expandBox = wx.BoxSizer(wx.HORIZONTAL)
        expandBox.Add(wx.StaticText(self.panel,label=' Expand -1 to +1 on: '),0,WACV)
        expand = wx.ComboBox(self.panel,value=self.Expand,choices=expandChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        expand.Bind(wx.EVT_COMBOBOX,OnExpand)
        expandBox.Add(expand,0,WACV)
        expandBox.Add(wx.StaticText(self.panel,label=' and find unique atoms '),0,WACV)
        mainSizer.Add(expandBox)

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
        return self.Trans,self.Vec,self.Expand

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

#==============================================================================
class DIFFaXcontrols(wx.Dialog):
    ''' Solicit items needed to prepare DIFFaX control.dif file
    '''
    def __init__(self,parent,ctrls,parms=None):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'DIFFaX controls',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.ctrls = ctrls
        self.calcType = 'powder pattern'
        self.plane = 'h0l'
        self.planeChoice = ['h0l','0kl','hhl','h-hl',]
        self.lmax = '2'
        self.lmaxChoice = [str(i+1) for i in range(6)]
        self.Parms = parms
        self.Parm = None
        if self.Parms != None:
            self.Parm = self.Parms[0]
        self.parmRange = [0.,1.]
        self.parmStep = 2
        self.Inst = 'Gaussian'
        self.Draw()

    def Draw(self):

        def OnCalcType(event):
            self.calcType = calcType.GetValue()
            wx.CallAfter(self.Draw)

        def OnPlane(event):
            self.plane = plane.GetValue()

        def OnMaxL(event):
            self.lmax = lmax.GetValue()

        def OnParmSel(event):
            self.Parm = parmsel.GetValue()

        def OnNumStep(event):
            self.parmStep = int(numStep.GetValue())

        def OnParmRange(event):
            event.Skip()
            vals = parmrange.GetValue().split()
            try:
                vals = [float(vals[0]),float(vals[1])]
            except ValueError:
                vals = self.parmRange
            parmrange.SetValue('%.3f %.3f'%(vals[0],vals[1]))
            self.parmRange = vals

        def OnInstSel(event):
            self.Inst = instsel.GetValue()

        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label=' Controls for DIFFaX'),0)
        if self.Parms:
            mainSizer.Add(wx.StaticText(self.panel,label=' Sequential powder pattern simulation'),0)
        else:
            calcChoice = ['powder pattern','selected area']
            calcSizer = wx.BoxSizer(wx.HORIZONTAL)
            calcSizer.Add(wx.StaticText(self.panel,label=' Select calculation type: '),0,WACV)
            calcType = wx.ComboBox(self.panel,value=self.calcType,choices=calcChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            calcType.Bind(wx.EVT_COMBOBOX,OnCalcType)
            calcSizer.Add(calcType,0,WACV)
            mainSizer.Add(calcSizer)
        if self.Parms:
            parmSel = wx.BoxSizer(wx.HORIZONTAL)
            parmSel.Add(wx.StaticText(self.panel,label=' Select parameter to vary: '),0,WACV)
            parmsel = wx.ComboBox(self.panel,value=self.Parm,choices=self.Parms,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            parmsel.Bind(wx.EVT_COMBOBOX,OnParmSel)
            parmSel.Add(parmsel,0,WACV)
            mainSizer.Add(parmSel)
            mainSizer.Add(wx.StaticText(self.panel,label=' Enter parameter range & no. steps: '))
            parmRange =  wx.BoxSizer(wx.HORIZONTAL)
            numChoice = [str(i+1) for i in range(10)]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            parmrange = wx.TextCtrl(self.panel,value='%.3f %.3f'%(self.parmRange[0],self.parmRange[1]),
                style=wx.TE_PROCESS_ENTER)
            parmrange.Bind(wx.EVT_TEXT_ENTER,OnParmRange)
            parmrange.Bind(wx.EVT_KILL_FOCUS,OnParmRange)
            parmRange.Add(parmrange,0,WACV)
            numStep = wx.ComboBox(self.panel,value=str(self.parmStep),choices=numChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            numStep.Bind(wx.EVT_COMBOBOX,OnNumStep)
            parmRange.Add(numStep,0,WACV)
            mainSizer.Add(parmRange)
        if 'selected' in self.calcType:
            planeSizer = wx.BoxSizer(wx.HORIZONTAL)
            planeSizer.Add(wx.StaticText(self.panel,label=' Select plane: '),0,WACV)
            plane = wx.ComboBox(self.panel,value=self.plane,choices=self.planeChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            plane.Bind(wx.EVT_COMBOBOX,OnPlane)
            planeSizer.Add(plane,0,WACV)
            planeSizer.Add(wx.StaticText(self.panel,label=' Max. l index: '),0,WACV)
            lmax = wx.ComboBox(self.panel,value=self.lmax,choices=self.lmaxChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            lmax.Bind(wx.EVT_COMBOBOX,OnMaxL)
            planeSizer.Add(lmax,0,WACV)
            mainSizer.Add(planeSizer)
        else:
            instChoice = ['None','Mean Gaussian','Gaussian',]
            instSizer = wx.BoxSizer(wx.HORIZONTAL)
            instSizer.Add(wx.StaticText(self.panel,label=' Select instrument broadening: '),0,WACV)
            instsel = wx.ComboBox(self.panel,value=self.Inst,choices=instChoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            instsel.Bind(wx.EVT_COMBOBOX,OnInstSel)
            instSizer.Add(instsel,0,WACV)
            mainSizer.Add(instSizer)
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
        if 'powder' in self.calcType:
            return 'PWDR',self.Inst,self.Parm,self.parmRange,self.parmStep
        elif 'selected' in self.calcType:
            return 'SADP',self.plane,self.lmax

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

#==============================================================================
class AddHatomDialog(wx.Dialog):
    '''H atom addition dialog. After :meth:`ShowModal` returns, the results
    are found in dict :attr:`self.data`, which is accessed using :meth:`GetData`.

    :param wx.Frame parent: reference to parent frame (or None)
    :param dict Neigh: a dict of atom names with list of atom name, dist pairs for neighboring atoms
    :param dict phase: a dict containing the phase as defined by
      :ref:`Phase Tree Item <Phase_table>`
    '''
    def __init__(self,parent,Neigh,phase):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'H atom add',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wxscroll.ScrolledPanel(self)         #just a dummy - gets destroyed in Draw!
        self.Neigh = Neigh
        self.phase = phase
        self.Hatoms = []
        self.Draw(self.Neigh,self.phase)

    def Draw(self,Neigh,phase):
        '''Creates the contents of the dialog. Normally called
        by :meth:`__init__`.
        '''
        def OnHSelect(event):
            Obj = event.GetEventObject()
            item,i = Indx[Obj.GetId()]
            for obj in Indx[item]:
                obj.SetValue(False)
            Obj.SetValue(True)
            self.Neigh[item][2] = i

        def OnBond(event):
            Obj = event.GetEventObject()
            inei,ibond = Indx[Obj.GetId()]
            self.Neigh[inei][1][0][ibond][2] = Obj.GetValue()

        self.panel.Destroy()
        self.panel = wxscroll.ScrolledPanel(self,style = wx.DEFAULT_DIALOG_STYLE)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,'H atom add controls for phase %s:'%(phase['General']['Name'])),
            0,wx.LEFT|wx.TOP,10)
        mainSizer.Add(wx.StaticText(self.panel,-1,'NB: Check selections as they may not be correct'),0|wx.LEFT,10)
        mainSizer.Add(wx.StaticText(self.panel,-1," Atom:  Add # H's          Use: Neighbors, dist"),0,wx.TOP|wx.LEFT,5)
        nHatms = ['0','1','2','3']
        dataSizer = wx.FlexGridSizer(0,3,0,0)
        Indx = {}
        for inei,neigh in enumerate(Neigh):
            dataSizer.Add(wx.StaticText(self.panel,-1,' %s:  '%(neigh[0])),0,WACV)
            nH = 1      #for O atom
            if 'C' in neigh[0] or 'N' in neigh[0]:
                nH = 4-len(neigh[1][0])
            checks = wx.BoxSizer(wx.HORIZONTAL)
            Ids = []
            for i in range(nH+1):
                nHs = wx.CheckBox(self.panel,-1,label=nHatms[i])
                if i == neigh[2]:
                    nHs.SetValue(True)
                Indx[nHs.GetId()] = [inei,i]
                Ids.append(nHs)
                nHs.Bind(wx.EVT_CHECKBOX, OnHSelect)
                checks.Add(nHs,0,WACV)
            Indx[inei] = Ids
            dataSizer.Add(checks,0,WACV)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            for ib,bond in enumerate(neigh[1][0]):
                Bond = wx.CheckBox(self.panel,-1,label=': %s, %.3f'%(bond[0],bond[1]))
                Bond.SetValue(bond[2])
                Indx[Bond.GetId()] = [inei,ib]
                Bond.Bind(wx.EVT_CHECKBOX,OnBond)
                lineSizer.Add(Bond,0,WACV)
            dataSizer.Add(lineSizer,0,WACV|wx.RIGHT,10)
        mainSizer.Add(dataSizer,0,wx.LEFT,5)

        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        OkBtn = wx.Button(self.panel,-1,'Ok')
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        size = np.array(self.GetSize())
        self.panel.SetupScrolling()
        self.panel.SetAutoLayout(1)
        size = [size[0]-5,size[1]-20]       #this fiddling is needed for older wx!
        self.panel.SetSize(size)

    def GetData(self):
        'Returns the values from the dialog'
        for neigh in self.Neigh:
            for ibond,bond in enumerate(neigh[1][0]):
                if not bond[2]:
                    neigh[1][1][1][ibond] = 0   #deselected bond
            neigh[1][1][1] = [a for a in  neigh[1][1][1] if a]
        return self.Neigh       #has #Hs to add for each entry

    def OnOk(self,event):
        'Called when the OK button is pressed'
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

#### Phase editing routines ################################################################################

def getPawleydRange(G2frame,data):
    'find d-space range in used histograms'
    fmtd = lambda d: '?' if d is None else '{:.5f}'.format(d)
    dmaxAll = dminAll = None
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    nhist = 0
    chist = 0
    for item in data['Histograms']:
        if 'HKLF' in item: #count single crystal histograms
            chist += 1
        if 'PWDR' not in item: continue
        if not data['Histograms'][item]['Use']: continue
        nhist += 1
        Inst = Histograms[item]['Instrument Parameters'][0]
        if 'T' in Inst['Type'][1]:
            dmin,dmax = [G2lat.Pos2dsp(Inst,t) for t in Histograms[item]['Limits'][1]]
        else:
            dmax,dmin = [G2lat.Pos2dsp(Inst,t) for t in Histograms[item]['Limits'][1]]
        if dmaxAll is None:
            dmaxAll = dmax
        else:
            dmaxAll = max(dmaxAll,dmax)
        if dminAll is None:
            dminAll = dmin
        else:
            dminAll = min(dminAll,dmin)
    # format data range
    lbl ="   d-space range {} to {} {}-1 ({} histograms)".format(
        fmtd(dminAll),fmtd(dmaxAll),Angstr,nhist)
    if dmaxAll is None: dmaxAll = 100.
    if dminAll is None: dminAll = 0.25
    return dminAll,dmaxAll,nhist,lbl,chist

def getAtomSelections(AtmTbl,cn=0,action='action',includeView=False,ask=True):
    '''get selected atoms from table or ask user if none are selected

        :param list AtmTbl: atom or draw atom table
        :param int cn: atom name position
        :param str action: description for prompt, when needed
        :param bool includeView: if True, the viewpoint is included
          as an option in the selection dialog
        :returns: indx (list) selected atoms from indices in table.
          If includeView is True, indx can contain index n (where there
          are n atoms in table). This is indicates the viewpoint.
    '''
    indx = AtmTbl.GetSelectedRows()
    indx += [row for row,col in AtmTbl.GetSelectedCells()]
    for top,bottom in zip([r for r,c in AtmTbl.GetSelectionBlockTopLeft()],
                          [r for r,c in AtmTbl.GetSelectionBlockBottomRight()]):
        indx += list(range(top,bottom+1))
    indx = list(set(indx))
    if indx or not ask: return indx
    choices = []
    for i in range(AtmTbl.GetNumberRows()):
        val = AtmTbl.GetCellValue(i,cn)
        if val in choices:
            val += '_' + str(i)
        choices.append(val)
    if not choices: return
    if includeView:
        choices.append('View point')
    dlg = G2G.G2MultiChoiceDialog(AtmTbl.GetTopLevelParent(),
        'Select atoms','Select atoms for '+action,choices)
    if dlg.ShowModal() == wx.ID_OK:
        indx = dlg.GetSelections()
    dlg.Destroy()
    return indx

def SetPhaseWindow(phasePage,mainSizer=None,Scroll=0):
    '''Finish off processing for all items going into a phase notebook page
    This connects the sizer to the Panel/ScrolledWindow that is assigned
    as the notebook's page for a tab.

    Note that a wx.ScrolledWindow is used for most tab pages, with the
    exception of Atoms, drawAtoms, G2frame.MapPeaks and G2frame.PawleyRefl,
    where a wx.Panel is used with a single Grid inside. This allows the grid
    to handle scrolling.

    When a wx.ScrolledWindows is used, scrolling is turned on here. The
    optional Scroll parameter is used to restore the scroll position to
    the previous position so that the window can be redrawn without
    disruption.
    '''
    if mainSizer is not None:
        phasePage.SetSizer(mainSizer)
    if isinstance(phasePage,wx.ScrolledWindow):
        phasePage.SetAutoLayout(True)
        phasePage.SetScrollRate(10,10)
        phasePage.SendSizeEvent()
        phasePage.Scroll(0,Scroll)

def GetSpGrpfromUser(parent,SpGrp):
    helptext = '''\t\t\tGSAS-II space group information

Space groups are entered here as given in Volume I or Volume A of the
International Tables using the short Hermann-Mauguin symbol,except that spaces
are placed between axial fields (e.g. "P 4/m m m", "F D 3 M" or "p -3 1 m").
NB: the cubic "bar" in "F d -3 m" is unnecessary, and upper/lower case is not required.

Where a centrosymmetric tetragonal or cubic space group has alternate origin settings,
Origin choice 2 (with the center of symmetry at the origin, which gives an -x,-y,-z
symmetry operator) is always used. Refer to the relevant pages in IT I or A to find
the offset in atom positions between the two choices.

For rhombohedral space groups, (R xxx) the hexagonal setting is assumed. Append a
final R to the name (R xxx R) to indicate that a rhombohedral cell should be
used (not recommended when alpha >> 120 or << 60, due to correlation.)

For standard settings of space groups, space group numbers (1-230) can alternately
be entered.

GSAS-II will accept non-standard settings of space groups. For example, space
group "P -1" can be set to include face centering, using symbol "F -1" and "P 1 1 21/a"
as a nonstandard version of "P 21/c".

Review the symmetry operators generated by GSAS-II to confirm that you have
entered the right symbol for your structure.
'''
    dlg = G2G.SingleStringDialog(parent,'Get Space Group',
        '  Input the space group with spaces between axial fields  \n  (e.g. p 21/c, P 63/m m c, P 4/m m m) or enter a space\n  group number between 1 and 230.',
        value=SpGrp,help=helptext)
    if not dlg.Show():
        dlg.Destroy()
        return SpGrp
    else:
        try:
            # has a space group number been input?
            spcnum = int(dlg.GetValue())
            if 1 <= spcnum <= 230:
                SpcGp = G2spc.spgbyNum[spcnum]
            else:
                msg = 'Space Group Error'
                wx.MessageBox('Invalid space group number',caption=msg,style=wx.ICON_EXCLAMATION)
                return
        except:
            #get rid of extra spaces between fields first
            Flds = dlg.GetValue().split()
            for fld in Flds: fld = fld.strip()
            SpcGp = ' '.join(Flds).capitalize()
        finally:
            dlg.Destroy()
    return SpcGp


def FindBondsDraw(data):
    '''Generally used routine where cell is from data
    '''
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    FindBondsDrawCell(data,cell)

def getAtomRadii(data):
    '''Get radii for atoms, using generalData['DisAglCtls']['BondRadii']
    to override generalData['BondRadii'] when present. Fix to make sure
    that all elements in generalData are present in DisAglCtls.
    '''
    generalData = data['General']
    if 'DisAglCtls' not in generalData:
        return generalData['AtomTypes'],generalData['BondRadii']
    if 'BondRadii' not in generalData['DisAglCtls']:
        return generalData['AtomTypes'],generalData['BondRadii']
    DisAglCtls = generalData['DisAglCtls']
    if len(generalData['BondRadii']) != len(DisAglCtls['BondRadii']):
        for typ,dis in zip(generalData['AtomTypes'],generalData['BondRadii']):
            if typ not in DisAglCtls['AtomTypes']:
                DisAglCtls['AtomTypes'].append(typ)
                DisAglCtls['AngleRadii'].append(dis)
                DisAglCtls['BondRadii'].append(dis)
    return DisAglCtls['AtomTypes'],DisAglCtls['BondRadii']

def FindCoordinationByLabel(data):
    '''Map out molecular connectivity by determining the atoms bonded
    to each atom, by label. The atoms bonded to each atom in the asymmetric
    unit is determined and returned in a dict. Works best
    '''
    generalData = data['General']
    cx,ct,cs,cia = generalData['AtomPtrs']
    atomTypes,radii = getAtomRadii(data)
    SGData = generalData['SGData']
    cellArray = G2lat.CellBlock(1)
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])

    neighborArray = {}
    error = ''
    for atomA in data['Atoms']:
        lblA = atomA[0]
        if lblA in neighborArray:
            if error: error += ', '
            error += lblA
        else:
            neighborArray[lblA] = []
        xyzA = np.array(atomA[cx:cx+3])
        indA = atomTypes.index(atomA[ct])
        for atomB in data['Atoms']:
            indB = atomTypes.index(atomB[ct])
            sumR = data['Drawing']['radiusFactor']*(radii[indA]+radii[indB])
            symAtms = [atm[0] for atm in  G2spc.GenAtom(np.array(atomB[cx:cx+3]),SGData,False,6*[0],True)]
            symCellAtms = np.concatenate([cellArray+i for i in symAtms])
            dists = np.sqrt(np.sum(np.inner(Amat,symCellAtms-xyzA)**2,axis=0))
            if np.any(np.logical_and(dists < sumR, dists != 0)):
                if atomB[0] not in neighborArray[lblA]:
                    neighborArray[lblA].append(atomB[0])
    if error:
        print('Warning, duplicated atom labels:',error)
    return neighborArray

def FindCoordination(ind,data,neighborArray,coordsArray,cmx=0,targets=None):
    '''Find atoms coordinating atom ind, speed-up version.
    This only searches to atoms already added to the Draw Array, though we might want
    to search to all atoms in the asymmetric unity (which would mean searching against
    atomsAll, but would also require a reformat of atom entry to match difference in
    format between atoms and drawatoms.
    '''
    generalData = data['General']
    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
    atomTypes,radii = getAtomRadii(data)
    atomData = data['Drawing']['Atoms']
    cx,ct,cs,ci = data['Drawing']['atomPtrs']
    #atomsAll = data['Atoms']
    #cxa,cta,csa,cia = generalData['AtomPtrs']
    SGData = generalData['SGData']
    cellArray = G2lat.CellBlock(1)
    newAtomList = []
    atomA = atomData[ind]
    xyzA = np.array(atomA[cx:cx+3])
    lblA = atomA[0]
    indA = atomTypes.index(atomA[ct])
    for atomB in atomData:
        if targets and atomB[ct] not in targets: continue
        if atomB[0] not in neighborArray[lblA]: continue
        indB = atomTypes.index(atomB[ct])
        sumR = data['Drawing']['radiusFactor']*(radii[indA]+radii[indB])
        xyzB = np.array(atomB[cx:cx+3])
        Uij = atomB[cs+5:cs+5+6]
        coords = []
        symMisc = []
        for item in G2spc.GenAtom(xyzB,SGData,False,Uij,True):
            coords.append(item[0])
            symMisc.append(item[1:4])
        symCoords = np.array(coords)
        dists = np.sqrt(np.sum(np.inner(Amat,
                np.array([symCoords+i-xyzA for i in cellArray])
                                           )**2,axis=0))
        for icell,isym in np.argwhere(np.logical_and(dists < sumR, dists != 0.0)):
            xyz = symCoords[isym] + cellArray[icell]
            item = [None]+list(symMisc[isym])
            atom = copy.deepcopy(atomB)
            atom[cx:cx+3] = xyz
            Opr = abs(item[2])%100
            M = SGData['SGOps'][Opr-1][0]
            if cmx:
                opNum = G2spc.GetOpNum(item[2],SGData)
                mom = np.array(atom[cmx:cmx+3])
                if SGData['SGGray']:
                    atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)
                else:
                    atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SGData['SpnFlp'][opNum-1]
            atom[cs-1] = str(item[2])+'+'
            atom[cs+5:cs+5+6] = item[1]
            # have we already found an atom at this site?
            if np.any(np.all(np.isclose(xyz,coordsArray,atol=0.0002),axis=1)): continue
            # are we going to add it already?
            if True in [np.allclose(np.array(xyz),np.array(a[cx:cx+3]),atol=0.0002) for a in newAtomList]: continue
            C = xyz-symCoords[isym]+item[3]
            atom[cs-1] += str(int(round(C[0])))+','+str(int(round(C[1])))+','+str(int(round(C[2])))
            newAtomList.append(atom)
    return newAtomList

def FindBondsDrawCell(data,cell):
    '''uses numpy & masks - very fast even for proteins!
    allows different cell as input from seq. refinements
    '''
    cx,ct,cs,ci = data['Drawing']['atomPtrs']
    hydro = data['Drawing']['showHydrogen']
    atomData = data['Drawing']['Atoms']
    generalData = data['General']
    Amat,Bmat = G2lat.cell2AB(cell)
    atomTypes,radii = getAtomRadii(data)
    try:
        indH = atomTypes.index('H')
        radii[indH] = 0.5
    except:
        pass
    for atom in atomData:
        atom[-2] = []               #clear out old bonds/polyhedra
        atom[-1] = []
    Indx = range(len(atomData))
    Atoms = []
    Styles = []
    Radii = []
    Names = []
    for atom in atomData:
        if 'Q' in atom[ct]:     #skip spinning RB atoms
            continue
        Atoms.append(np.array(atom[cx:cx+3]))
        Styles.append(atom[cs])
        Names.append(ord(atom[ct-1].ljust(4)[3]))
        try:
            if not hydro and atom[ct] == 'H':
                Radii.append(0.0)
            else:
                Radii.append(radii[atomTypes.index(atom[ct])])
        except ValueError:          #changed atom type!
            Radii.append(0.20)
    Atoms = np.array(Atoms)
    Radii = np.array(Radii)
    Names = np.array(Names)
    IASRN = zip(Indx,Atoms,Styles,Radii,Names)
    for atomA in IASRN:
        if atomA[2] in ['lines','sticks','ellipsoids','balls & sticks','polyhedra']:
            Dx = Atoms-atomA[1]
            dist = ma.masked_less(np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0)),0.5) #gets rid of G2frame & disorder "bonds" < 0.5A
            if generalData['Type'] == 'macromolecular':     #eliminate cross disorder residue bonds
                m1 = ma.getmask(dist)
                if atomA[4] in [ord('A'),]:
                    m2 = ma.getmask(ma.masked_equal(Names,ord('B')))
                    dist = ma.array(dist,mask=ma.mask_or(m1,m2))
                if atomA[4] in [ord('B'),]:
                    m2 = ma.getmask(ma.masked_equal(Names,ord('A')))
                    dist = ma.array(dist,mask=ma.mask_or(m1,m2))
            sumR = atomA[3]+Radii
            IndB = ma.nonzero(ma.masked_greater(dist-data['Drawing']['radiusFactor']*sumR,0.))                 #get indices of bonded atoms
            i = atomA[0]
            for j in IndB[0]:
                if Styles[i] == 'polyhedra':
                    atomData[i][-2].append(np.inner(Amat,Dx[j]))
                elif Styles[j] != 'polyhedra' and j > i:
                    atomData[i][-2].append(Dx[j]*Radii[i]/sumR[j])
                    atomData[j][-2].append(-Dx[j]*Radii[j]/sumR[j])
            if Styles[i] == 'polyhedra':
                Bonds = atomData[i][-2]
                Faces = []
                if len(Bonds) > 2:
                    FaceGen = G2lat.uniqueCombinations(Bonds,3)     #N.B. this is a generator
                    for face in FaceGen:
                        vol = nl.det(face)
                        if abs(vol) > 1. or len(Bonds) == 3:
                            if vol < 0.:
                                face = [face[0],face[2],face[1]]
                            face = np.array(face)
                            if not np.array([np.array(nl.det(face-bond))+0.0001 < 0 for bond in Bonds]).any():
                                norm = np.cross(face[1]-face[0],face[2]-face[0])
                                norm /= np.sqrt(np.sum(norm**2))
                                Faces.append([face,norm])
                    atomData[i][-1] = Faces

def VoidMap(data,aMax=1,bMax=1,cMax=1,gridspacing=.25,probeRadius=.5,
                aMin=0,bMin=0,cMin=0):
    '''Compute points where there are no atoms within probeRadius A.
    All atoms in the Atoms list are considered, provided their
    occupancy is non-zero.

    :param dict data: Phase data array
    :param float aMax: Maximum along the *a* direction (fractional units).
       Defaults to 1.
    :param float bMax: Maximum along the *b* direction (fractional units).
       Defaults to 1.
    :param float cMax: Maximum along the *c* direction (fractional units).
       Defaults to 1.
    :param float gridspacing=.25: Approximate spacing of points (fractional units).
       Defaults to 1.
    :param float ,probeRadius=.5:
    :param float aMin: Minimum along the *a* direction (fractional units).
       Defaults to 0.
    :param float bMin: Minimum along the *b* direction (fractional units).
       Defaults to 0.
    :param float cMin: Minimum along the *c* direction (fractional units).
       Defaults to 0.
    '''

    VDWdict = dict(zip(data['General']['AtomTypes'],data['General']['vdWRadii']))
    cell = data['General']['Cell'][1:7]
    Amat,Bmat = G2lat.cell2AB(cell) # orthogonalization matrix
    SGData = data['General']['SGData']
    surroundingCells = G2lat.CellBlock(1)

    xx,yy,zz = np.meshgrid(
        np.linspace(aMin,aMax,int(0.5+cell[0]*(aMax-aMin)/gridspacing),endpoint=False),
        np.linspace(bMin,bMax,int(0.5+cell[1]*(bMax-bMin)/gridspacing),endpoint=False),
        np.linspace(cMin,cMax,int(0.5+cell[2]*(cMax-cMin)/gridspacing),endpoint=False))
    coordGrd = np.array([xyz for xyz in zip(xx.ravel(),yy.ravel(),zz.ravel())])

    lgclArray = [True for i in xx.ravel()]

    cx,ct,cs,cia = data['General']['AtomPtrs']
    nind = len(data['Atoms'])
    pgbar = wx.ProgressDialog('Fill unit cell for %d atoms'%nind,'Atoms done=',nind+1,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
    screenSize = wx.ClientDisplayRect()
    Size = pgbar.GetSize()
    if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
        pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))

    for i,atom in enumerate(data['Atoms']):
        if atom[cx+3] <= 0: continue
        radius = VDWdict.get(atom[ct])
        cellMin = -radius/np.array(cell[0:3])
        cellMax = radius/np.array(cell[0:3]) + (aMax,bMax,cMax)
        if radius is None:
            print('Skipping atom {}, no radius'.format(atom[0]))
            continue
        radius += probeRadius
        result = G2spc.GenAtom(atom[cx:cx+3],SGData,Move=True)
        for item in result:
            for scell in surroundingCells:
                XYZ = item[0] + scell
                if np.any((XYZ < cellMin, XYZ > cellMax)): continue
                lgclArray = np.logical_and(lgclArray,np.sqrt(np.sum(np.inner(Amat,coordGrd-XYZ)**2,axis=0))>radius)
        GoOn = pgbar.Update(i,newmsg='Atoms done=%d'%(i))
        if not GoOn[0]:
            break
    pgbar.Destroy()
    print('found ',len(coordGrd[lgclArray]),'points gridspacing,probeRadius=',gridspacing,probeRadius)
    return coordGrd[lgclArray]

def SetDrawingDefaults(drawingData):
    """Add required items into data['drawing'] array if not present. This does not add
    all the items in SetupDrawingData, but it seems that this is not a problem. Perhaps the
    two routines could be combined?
    """
    defaultDrawing = {'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,
            'backColor':[0,0,0],'depthFog':False,'Zclip':50.0,'cameraPos':50.,'Zstep':0.5,
            'radiusFactor':0.85,'contourLevel':1.,'bondRadius':0.1,'ballScale':0.33,
            'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,'unitCellBox':True,
            'showABC':True,'selectedAtoms':[],'Atoms':[],'oldxy':[],'magMult':1.0,
            'bondList':{},'viewDir':[1,0,0],'Plane':[[0,0,1],False,False,0.0,[255,255,0]],
            'peakMoveView':True,'PeakDistRadius':0.0,'showVoids':False,'showMap':False,
            'atomsExpandRadius':5.,'atomsDistRadius':2.5,'Voids':[],
            'VPPeakDistRad':0.,'VPatomsExpandRad':0.,'VPatomsDistRad':0.,
        }
    for key in defaultDrawing:
        if key not in drawingData: drawingData[key] = defaultDrawing[key]

def updateAddRBorientText(G2frame,testRBObj,Bmat):
    '''Update all origin/orientation text on the Add RB panel or
    on main RB Models page in response to Alt+mouse movement
    '''
    A,V = G2mth.Q2AVdeg(testRBObj['rbObj']['Orient'][0])
    testRBObj['rbObj']['OrientVec'][0] = A
    testRBObj['rbObj']['OrientVec'][1:] = np.inner(Bmat,V)
    for i,val in enumerate(testRBObj['rbObj']['OrientVec']):
        G2frame.testRBObjSizers['OrientVecSiz'][i].ChangeValue(val)
#        G2frame.testRBObjSizers['OrientVecSiz'][i].SetValue(val)
    try:
#        G2frame.testRBObjSizers['OrientVecSiz'][4].SetValue(
        G2frame.testRBObjSizers['OrientVecSiz'][4].ChangeValue(
            int(10*testRBObj['rbObj']['OrientVec'][0]))
    except:
        pass
    for i,sizer in enumerate(G2frame.testRBObjSizers['Xsizers']):
        sizer.ChangeValue(testRBObj['rbObj']['Orig'][0][i])
#        sizer.SetValue(testRBObj['rbObj']['Orig'][0][i])
    # redraw asymmetric unit when called on an existing body
    if G2frame.testRBObjSizers.get('OnOrien') is None: return
    G2frame.testRBObjSizers['OnOrien'](mode=testRBObj['rbObj'].get('drawMode',DrawStyleChoice[4]))

def GetReflData(G2frame,phaseName,reflNames):
    ReflData = {'RefList':[],'Type':''}
    if '' in reflNames:
        return None
    for reflName in reflNames:
        if 'PWDR' in reflName:
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
            if not PatternId:       #got 0
                return None
            reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflData = reflSets[phaseName]
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
            if not PatternId:       #got 0
                return None
            reflData = G2frame.GPXtree.GetItemPyData(PatternId)[1]
            if 'Type' not in reflData:
                reflData['Type'] = 'SXC'
        if ReflData['Type'] and reflData['Type'] != ReflData['Type']:
            G2frame.ErrorDialog('Data type conflict',
                reflName+' conflicts with previous '+ReflData['Type'])
            return None
        ReflData['RefList'] += list(reflData['RefList'])
        ReflData['Type'] = reflData['Type']
    return ReflData
def UpdatePhaseData(G2frame,Item,data):
    '''Create the data display window contents when a phase is clicked on
    in the main (data tree) window.
    Called only from :meth:`GSASIIdataGUI.SelectDataTreeItem`,
    which in turn is called from :meth:`GSASIIdataGUI.GSASII.OnDataTreeSelChanged`
    when a Phase tree item is selected. This creates all tabs on the page and fills
    their contents. Routine OnPageChanged is called each time a tab is pressed
    and updates the contents of the tab's page.

    :param wx.frame G2frame: the main GSAS-II frame object
    :param wx.TreeItemId Item: the tree item that was selected
    :param dict data: all the information on the phase in a dictionary

    '''

    def SetupGeneral():
        try:
            G2elem.SetupGeneral(data,G2frame.dirname)
        except ValueError as msg:
            wx.MessageBox(msg,caption='Element symbol error')

#### General phase routines ################################################################################
    def UpdateGeneral(Scroll=0,SkipDraw=False):
        '''Draw the controls for the General phase data subpage
        '''

        """ This is the default dictionary structure for phase data
        (taken from GSASII.py)
        'General':{
            'Name':PhaseName
            'Type':'nuclear'
            'SGData':SGData
            'Cell':[False,10.,10.,10.,90.,90.,90,1000.]
            'AtomPtrs':[]
            'Pawley dmin':1.0,
            'Pawley neg wt':0.0}
        'Atoms':[]
        'Drawing':{}
        """
        def NameSizer():

            def SetDefaultSSsymbol():
                if generalData['SGData']['SGLaue'] in '-1':
                    return '(abg)'
                elif generalData['SGData']['SGLaue'] in ['2/m']:
                    if generalData['SGData']['SGUniq'] == 'a':
                        return '(a00)'
                    elif generalData['SGData']['SGUniq'] == 'b':
                        return '(0b0)'
                    elif generalData['SGData']['SGUniq'] == 'c':
                        return '(00g)'
                else:
                    return '(00g)'

            def OnPhaseName(event):
                'called when the phase name is changed in "General"'
                event.Skip()
                newName = NameTxt.GetValue().strip()
                renamePhaseName(G2frame,data,Item,generalData,newName)
                NameTxt.SetValue(newName)
                                                
            def OnPhaseType(event):
                if not len(generalData['AtomTypes']):             #can change only if no atoms!
                    generalData['Type'] = TypeTxt.GetValue()
                    pages = [G2frame.phaseDisplay.GetPageText(PageNum) for PageNum in range(G2frame.phaseDisplay.GetPageCount())]
                    if generalData['Type'] == 'faulted':
                        G2frame.Bind(wx.EVT_MENU, OnLoadDIFFaX, id=G2G.wxID_LOADDIFFAX)
                        G2frame.Bind(wx.EVT_MENU, OnSimulate, id=G2G.wxID_LAYERSIMULATE)
                        G2frame.Bind(wx.EVT_MENU, OnSeqSimulate, id=G2G.wxID_SEQUENCESIMULATE)
                        G2frame.Bind(wx.EVT_MENU, OnFitLayers, id=G2G.wxID_LAYERSFIT)
                        if 'Wave Data' in pages:
                            pass
#                            G2frame.phaseDisplay.DeletePage(pages.index('Wave Data'))
                        if 'MC/SA' in pages:
                            pass
#                            G2frame.phaseDisplay.DeletePage(pages.index('MC/SA'))
                        if 'RB Models' in pages:
                            pass
#                            G2frame.phaseDisplay.DeletePage(pages.index('RB Models'))
                        if 'Layers' not in pages:
                            if 'Layers' not in data:
                                data['Layers'] = {'Laue':'-1','Cell':[False,1.,1.,1.,90.,90.,90,1.],
                                    'Width':[[1.,1.],[False,False]],'Toler':0.01,'AtInfo':{},
                                    'Layers':[],'Stacking':[],'Transitions':[]}
                            G2frame.layerData = wx.ScrolledWindow(G2frame.phaseDisplay)
                            G2frame.phaseDisplay.InsertPage(3,G2frame.layerData,'Layers')
                            Id = wx.NewId()
                            TabSelectionIdDict[Id] = 'Layers'
                        wx.CallAfter(UpdateGeneral)
                    elif generalData['Type'] == 'magnetic':
                        generalData['AtomPtrs'] = [3,1,10,12]
                        SGData = generalData['SGData']
                        Nops = len(SGData['SGOps'])*len(SGData['SGCen'])
                        if SGData['SGInv']:
                            Nops *= 2
                        SGData['SpnFlp'] = Nops*[1,]
                    else:
                        if 'Wave Data' in pages:
                            G2frame.phaseDisplay.DeletePage(pages.index('Wave Data'))
                        if 'MC/SA' not in pages:
                            G2frame.MCSA = wx.ScrolledWindow(G2frame.phaseDisplay)
                            G2frame.phaseDisplay.InsertPage(7,G2frame.MCSA,'MC/SA')
                            Id = wx.NewId()
                            TabSelectionIdDict[Id] = 'MC/SA'
                        wx.CallAfter(UpdateGeneral)
                else:
                    G2frame.ErrorDialog('Phase type change error','Can change phase type only if there are no atoms')
                    TypeTxt.SetValue(generalData['Type'])

            def OnSpaceGroup(event):
                if generalData['SGData']['SGFixed']:
                    msg = 'Fixed cif generated magnetic space group'
                    text = 'space group can not be changed'
                    wx.MessageBox(text,caption=msg,style=wx.ICON_EXCLAMATION)
                    text,table = G2spc.SGPrint(generalData['SGData'])
                    SGTxt.SetLabel(generalData['SGData']['SpGrp'])
                    msg = 'cif based Space Group Information'
                    G2G.SGMessageBox(General,msg,text,table).Show()
                    return
                # try a lookup on the user-supplied name
                SpcGp = GetSpGrpfromUser(General,SpGrp)
                if SpcGp == SpGrp:
                    text,table = G2spc.SGPrint(generalData['SGData'])
                    SGTxt.SetLabel(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Information'
                    G2G.SGMessageBox(General,msg,text,table).Show()
                    return      #unchanged - do nothing but show info
                SpGrpNorm = G2spc.StandardizeSpcName(SpcGp)
                if SpGrpNorm:
                    SGErr,SGData = G2spc.SpcGroup(SpGrpNorm)
                else:
                    SGErr,SGData = G2spc.SpcGroup(SpcGp)
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    SGTxt.SetLabel(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Error'
                    Text = '\n'.join(text)
                    wx.MessageBox(Text,caption=msg,style=wx.ICON_EXCLAMATION)
                else:
                    if "1'" in SpcGp:
                        generalData['Type'] = 'magnetic'
                        generalData['Modulated'] = True
                    if generalData['Type'] == 'magnetic':
                        Nops = len(SGData['SGOps'])*len(SGData['SGCen'])
                        if SGData['SGInv']:
                            Nops *= 2
                        GenSym,GenFlg = G2spc.GetGenSym(SGData)[:2]
                        SGData['GenSym'] = GenSym
                        SGData['GenFlg'] = GenFlg
                        SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
                        G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
                        SGData['SpnFlp'] = Nops*[1,]
                    if generalData['Modulated']:
                        generalData['SuperSg'] = SetDefaultSSsymbol()
                        generalData['SSGData'] = G2spc.SSpcGroup(generalData['SGData'],generalData['SuperSg'])[1]
                        if SGData['SGGray']:
                            SGData['SpnFlp'] += Nops*[-1,]
                    text,table = G2spc.SGPrint(SGData)
                    generalData['SGData'] = SGData
                    SGTxt.SetLabel(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Information'
                    G2G.SGMessageBox(General,msg,text,table).Show()
                G2spc.UpdateSytSym(data)
                NShkl = len(G2spc.MustrainNames(SGData))
                NDij = len(G2spc.HStrainNames(SGData))
                for hist in data['Histograms']:
                    if 'HStrain' in data['Histograms'][hist]:       #PWDR only
                        data['Histograms'][hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
                        data['Histograms'][hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
                if data['Drawing']: data['Drawing']['Atoms'] = []
                wx.CallAfter(UpdateGeneral)

            def OnModulated(event):
                if not len(generalData['AtomTypes']):             #can change only if no atoms!
                    pages = [G2frame.phaseDisplay.GetPageText(PageNum) for PageNum in range(G2frame.phaseDisplay.GetPageCount())]
                    if generalData['Type'] in ['nuclear','magnetic']:
                        generalData['Modulated'] = modulated.GetValue()
                        if generalData['Modulated']:
                            if 'SuperSg' not in generalData:
                                generalData['SuperSg'] = SetDefaultSSsymbol()
                            generalData['SSGData'] = G2spc.SSpcGroup(generalData['SGData'],generalData['SuperSg'])[1]
                            if 'SuperVec' not in generalData:
                                generalData['Super'] = True
                                generalData['SuperVec'] = [[0.,0.,0.],False,4]
                                generalData['SSGData'] = {}
                            if '4DmapData' not in generalData:
                                generalData['4DmapData'] = mapDefault.copy()
                                generalData['4DmapData'].update({'MapType':'Fobs'})
                            if 'MC/SA' in pages:
                                pass
    #                            G2frame.phaseDisplay.DeletePage(pages.index('MC/SA'))   #this crashes!!
                            if 'Layers' in pages:
                                pass
    #                            G2frame.phaseDisplay.DeletePage(pages.index('Layers'))
                            if 'Wave Data' not in pages:
                                G2frame.waveData = wx.ScrolledWindow(G2frame.phaseDisplay)
                                G2frame.phaseDisplay.InsertPage(3,G2frame.waveData,'Wave Data')
                                Id = wx.NewId()
                                TabSelectionIdDict[Id] = 'Wave Data'
# deleting page now causes Mac crash, postpone until page is redrawn
#                        else:
#                            if 'Wave Data' in pages:
#                                G2frame.phaseDisplay.DeletePage(pages.index('Wave Data'))
                        wx.CallAfter(UpdateGeneral)
                else:
                    if generalData['Type'] == 'magnetic':
                        pages = [G2frame.phaseDisplay.GetPageText(PageNum) for PageNum in range(G2frame.phaseDisplay.GetPageCount())]
                        generalData['Modulated'] = modulated.GetValue()
                        if generalData['Modulated']:
                            if 'SuperSg' not in generalData:
                                generalData['SuperSg'] = SetDefaultSSsymbol()
                            generalData['SSGData'] = G2spc.SSpcGroup(generalData['SGData'],generalData['SuperSg'])[1]
                            if 'SuperVec' not in generalData:
                                generalData['Super'] = 1
                                generalData['SuperVec'] = [[0.,0.,0.],False,4]
                                generalData['SSGData'] = {}
                            if '4DmapData' not in generalData:
                                generalData['4DmapData'] = mapDefault.copy()
                                generalData['4DmapData'].update({'MapType':'Fobs'})
                            if 'Wave Data' not in pages:
                                G2frame.waveData = wx.ScrolledWindow(G2frame.phaseDisplay)
                                G2frame.phaseDisplay.InsertPage(3,G2frame.waveData,'Wave Data')
                                Id = wx.NewId()
                                TabSelectionIdDict[Id] = 'Wave Data'
                        Atoms = data['Atoms']
                        for atom in Atoms:
                            atom += [{'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}]
                        wx.CallAfter(UpdateGeneral)
                    else:
                        G2frame.ErrorDialog('Modulation type change error','Can change modulation only if there are no atoms')
                        modulated.SetValue(generalData['Modulated'])

            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(General,-1,' Phase name: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            NameTxt = wx.TextCtrl(General,-1,value=generalData['Name'],style=wx.TE_PROCESS_ENTER)
            NameTxt.Bind(wx.EVT_TEXT_ENTER,OnPhaseName)
            NameTxt.Bind(wx.EVT_KILL_FOCUS,OnPhaseName)
            nameSizer.Add(NameTxt,0,WACV)
            nameSizer.Add(wx.StaticText(General,-1,'  Phase type: '),0,WACV)
            TypeTxt = wx.ComboBox(General,-1,value=generalData['Type'],choices=phaseTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            TypeTxt.Bind(wx.EVT_COMBOBOX, OnPhaseType)
            nameSizer.Add(TypeTxt,0,WACV)
            nameSizer.Add(wx.StaticText(General,-1,'  Space group: '),0,WACV)
            SpGrp = generalData['SGData']['SpGrp']
            if generalData['SGData']['SGGray']: SpGrp += " 1'"
            SGTxt = wx.Button(General,wx.ID_ANY,SpGrp,size=(100,-1))
            SGTxt.Bind(wx.EVT_BUTTON,OnSpaceGroup)
            nameSizer.Add(SGTxt,0,WACV)
            if generalData['Type'] in ['nuclear','magnetic']:
                modulated = wx.CheckBox(General,label='Modulated? ')
                modulated.SetValue(generalData['Modulated'])
                modulated.Bind(wx.EVT_CHECKBOX,OnModulated)
                nameSizer.Add(modulated,0,WACV)
            return nameSizer

        def CellSizer():
            cellGUIlist = [[['m3','m3m'],4,zip([" Unit cell: a = "," Vol = "],["%.5f","%.3f"],[True,False],[0,0])],
            [['3R','3mR'],6,zip([" a = "," alpha = "," Vol = "],["%.5f","%.3f","%.3f"],[True,True,False],[0,3,0])],
            [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],6,zip([" a = "," c = "," Vol = "],["%.5f","%.5f","%.3f"],[True,True,False],[0,2,0])],
            [['mmm'],8,zip([" a = "," b = "," c = "," Vol = "],["%.5f","%.5f","%.5f","%.3f"],
                [True,True,True,False],[0,1,2,0])],
            [['2/m'+'a'],10,zip([" a = "," b = "," c = "," alpha = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,3,0])],
            [['2/m'+'b'],10,zip([" a = "," b = "," c = "," beta = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
            [['2/m'+'c'],10,zip([" a = "," b = "," c = "," gamma = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,5,0])],
            [['-1'],7,zip([" a = "," b = "," c = "," Vol = "," alpha = "," beta = "," gamma = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f","%.3f","%.3f"],
                [True,True,True,False,True,True,True],[0,1,2,0,3,4,5])]]

            def OnCellRef(event):
                generalData['Cell'][0] = cellRef.GetValue()

            def OnCellChange(invalid,value,tc):
                SGData = generalData['SGData']
                laue = SGData['SGLaue']
                if laue == '2/m':
                    laue += SGData['SGUniq']
                cell = generalData['Cell']
                Obj = tc
                ObjId = cellList.index(Obj.GetId())
                try:
                    value = max(1.0,float(tc.GetValue()))
                except ValueError:
                    if ObjId < 3:               #bad cell edge - reset
                        value = cell[ObjId+1]
                    else:                       #bad angle
                        value = 90.
                if laue in ['m3','m3m']:
                    cell[1] = cell[2] = cell[3] = value
                    cell[4] = cell[5] = cell[6] = 90.0
                    Obj.ChangeValue(cell[1])
                elif laue in ['3R','3mR']:
                    if ObjId == 0:
                        cell[1] = cell[2] = cell[3] = value
                        Obj.ChangeValue(cell[1])
                    else:
                        cell[4] = cell[5] = cell[6] = value
                        Obj.ChangeValue(cell[4])
                elif laue in ['3','3m1','31m','6/m','6/mmm','4/m','4/mmm']:
                    cell[4] = cell[5] = 90.
                    cell[6] = 120.
                    if laue in ['4/m','4/mmm']:
                        cell[6] = 90.
                    if ObjId == 0:
                        cell[1] = cell[2] = value
                        Obj.ChangeValue(cell[1])
                    else:
                        cell[3] = value
                        Obj.ChangeValue(cell[3])
                elif laue in ['mmm']:
                    cell[ObjId+1] = value
                    cell[4] = cell[5] = cell[6] = 90.
                    Obj.ChangeValue(cell[ObjId+1])
                elif laue in ['2/m'+'a']:
                    cell[5] = cell[6] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.ChangeValue(cell[ObjId+1])
                    else:
                        cell[4] = value
                        Obj.ChangeValue(cell[4])
                elif laue in ['2/m'+'b']:
                    cell[4] = cell[6] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.ChangeValue(cell[ObjId+1])
                    else:
                        cell[5] = value
                        Obj.ChangeValue(cell[5])
                elif laue in ['2/m'+'c']:
                    cell[4] = cell[5] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.ChangeValue(cell[ObjId+1])
                    else:
                        cell[6] = value
                        Obj.ChangeValue(cell[6])
                else:
                    cell[ObjId+1] = value
                    Obj.ChangeValue(cell[1+ObjId])
                cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
                volVal.ChangeValue("%.3f"%(cell[7]))
                density,mattCoeff = G2mth.getDensity(generalData)
                if denSizer:
                    denSizer[1].ChangeValue('%.3f'%(density))
                    if denSizer[2]:
                        denSizer[2].ChangeValue('%.3f'%(mattCoeff))

            cell = generalData['Cell']
            laue = generalData['SGData']['SGLaue']
            if laue == '2/m':
                laue += generalData['SGData']['SGUniq']
            for cellGUI in cellGUIlist:
                if laue in cellGUI[0]:
                    useGUI = cellGUI
            cellSizer = wx.FlexGridSizer(0,useGUI[1]+1,5,5)
            if PWDR:
                cellRef = wx.CheckBox(General,-1,label='Refine unit cell:')
                cellSizer.Add(cellRef,0,WACV)
                cellRef.Bind(wx.EVT_CHECKBOX, OnCellRef)
                cellRef.SetValue(cell[0])
            cellList = []
            for txt,fmt,ifEdit,Id in useGUI[2]:
                cellSizer.Add(wx.StaticText(General,label=txt),0,WACV)
                Fmt = (10,5)
                if '.3' in fmt:
                    Fmt = (10,3)
                if ifEdit:          #a,b,c,etc.
                    cellVal = G2G.ValidatedTxtCtrl(General,generalData['Cell'],Id+1,
                            xmin=0.1,xmax=500.,nDig=Fmt,OnLeave=OnCellChange)
                    cellSizer.Add(cellVal,0,WACV)
                    cellList.append(cellVal.GetId())
                else:               #volume
                    volVal = G2G.ReadOnlyTextCtrl(General,value=(fmt%(cell[7])))
                    cellSizer.Add(volVal,0,WACV)
            return cellSizer

        def ElemSizer():

            def OnIsotope(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                isotope = Obj.GetValue()
                nCols = len(generalData['AtomTypes'])+1
                data['General']['Isotope'][item] = isotope
                indx = generalData['AtomTypes'].index(item)
                wt = generalData['Isotopes'][item][isotope]['Mass']
                elemSizer.GetChildren()[indx+3*nCols+1].Window.ChangeValue('%.3f'%(wt))    #tricky
                data['General']['AtomMass'][indx] = wt
                density,mattCoeff = G2mth.getDensity(generalData)
                denSizer[1].ChangeValue('%.3f'%(density))
                if denSizer[2]:
                    denSizer[2].ChangeValue('%.3f'%(mattCoeff))

            def onDefColor(event):
                '''Called when a color bar in elements table is clicked on.
                Changes default color for element in all phases
                N.B. Change is not saved; will go back to original color when
                GSAS-II is restarted.
                Changes colors of matching atoms in Draw Atoms table for
                current phase -- only.
                '''
                if not hasattr(event.GetEventObject(),'atomNum'): return
                anum = event.GetEventObject().atomNum
                (R,G,B) = generalData['Color'][anum]
                dlg = wx.ColourDialog(event.GetEventObject().GetTopLevelParent())
                dlg.GetColourData().SetChooseFull(False)
                dlg.GetColourData().SetColour(wx.Colour(R,G,B))
                if dlg.ShowModal() == wx.ID_OK:
                    El = generalData['AtomTypes'][anum]
                    RGB = dlg.GetColourData().GetColour()[0:3]
                    G2elem.SetAtomColor(El, RGB)
                    cx,ct,cs,ci = data['Drawing']['atomPtrs']
                    for atom in data['Drawing']['Atoms']:
                        if atom[ct] != El: continue
                        atom[cs+2] = RGB
                    wx.CallAfter(UpdateGeneral)
                dlg.Destroy()

            elemSizer = wx.FlexGridSizer(0,len(generalData['AtomTypes'])+1,1,1)
            elemSizer.Add(wx.StaticText(General,label=' Elements'),0,WACV)
            for elem in generalData['AtomTypes']:
                typTxt = G2G.ReadOnlyTextCtrl(General,value=elem)
                elemSizer.Add(typTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Isotope'),0,WACV)
            for elem in generalData['AtomTypes']:
                choices = list(generalData['Isotopes'][elem].keys())
                isoSel = wx.ComboBox(General,-1,value=generalData['Isotope'][elem],choices=choices,
                    style=wx.CB_READONLY)
                isoSel.Bind(wx.EVT_COMBOBOX,OnIsotope)
                Indx[isoSel.GetId()] = elem
                elemSizer.Add(isoSel,1,wx.EXPAND)
            elemSizer.Add(wx.StaticText(General,label=' No. per cell'),0,WACV)
            for elem in generalData['AtomTypes']:
                numbTxt = G2G.ReadOnlyTextCtrl(General,value='%.1f'%(generalData['NoAtoms'][elem]))
                elemSizer.Add(numbTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Atom weight'),0,WACV)
            for wt in generalData['AtomMass']:
                wtTxt = G2G.ReadOnlyTextCtrl(General,value='%.3f'%(wt))
                elemSizer.Add(wtTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Bond radii'),0,WACV)
            for rad in generalData['BondRadii']:
                bondRadii = G2G.ReadOnlyTextCtrl(General,value='%.2f'%(rad))
                elemSizer.Add(bondRadii,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Angle radii'),0,WACV)
            for rad in generalData['AngleRadii']:
                elemTxt = G2G.ReadOnlyTextCtrl(General,value='%.2f'%(rad))
                elemSizer.Add(elemTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' van der Waals radii'),0,WACV)
            for rad in generalData['vdWRadii']:
                elemTxt = G2G.ReadOnlyTextCtrl(General,value='%.2f'%(rad))
                elemSizer.Add(elemTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Default color'),0,WACV)
            for i,(R,G,B) in enumerate(generalData['Color']):
                colorTxt = G2G.ReadOnlyTextCtrl(General,value='')
                colorTxt.SetBackgroundColour(wx.Colour(R,G,B))
                colorTxt.atomNum = i
                colorTxt.Bind(wx.EVT_SET_FOCUS,onDefColor)
                elemSizer.Add(colorTxt,0,WACV)
            if generalData['Type'] == 'magnetic':
                elemSizer.Add(wx.StaticText(General,label=' Lande g factor: '),0,WACV)
                for ig,elem in enumerate(generalData['AtomTypes']):
                    gfac = generalData['Lande g'][ig]
                    if gfac == None:
                        elemSizer.Add((5,0),)
                    else:
                        gfacTxt = G2G.ValidatedTxtCtrl(General,generalData['Lande g'],ig,
                            xmin=0.5,xmax=3.0,nDig=(10,2))
                        elemSizer.Add(gfacTxt,0,WACV)
            return elemSizer

        def DenSizer():

            generalData['Mass'] = G2mth.getMass(generalData)
            density,mattCoeff = G2mth.getDensity(generalData)
            denSizer = wx.BoxSizer(wx.HORIZONTAL)
            denSizer.Add(wx.StaticText(General,-1,' Density: '),0,WACV)
            denTxt = G2G.ReadOnlyTextCtrl(General,-1,'%.3f'%(density))
            denSizer.Add(denTxt,0,WACV)
            mattTxt = None
            if generalData['Type'] == 'macromolecular' and generalData['Mass'] > 0.0:
                denSizer.Add(wx.StaticText(General,-1,' Matthews coeff.: '),
                    0,WACV)
                mattTxt = G2G.ReadOnlyTextCtrl(General,-1,'%.3f'%(mattCoeff))
                denSizer.Add(mattTxt,0,WACV)
            return denSizer,denTxt,mattTxt

        def MagSizer():

            def OnSpinOp(event):
                if SGData['SGFixed']:
                    msg = 'Fixed cif generated spins'
                    text = 'Spin configuration can not be changed; it will be reset'
                    wx.MessageBox(text,caption=msg,style=wx.ICON_EXCLAMATION)
                    wx.CallAfter(UpdateGeneral)
                    return
                Obj = event.GetEventObject()
                isym = Indx[Obj.GetId()]+1
                spCode = {'red':-1,'black':1}
                SGData['SGSpin'][isym] = spCode[Obj.GetValue()]
                G2spc.CheckSpin(isym,SGData)
                wx.CallAfter(UpdateGeneral)

            def OnBNSlatt(event):
                Obj = event.GetEventObject()
                BNSlatt = Obj.GetValue()
                SGData = generalData['SGData']
                SpcGrp = SGData['SpGrp']
                if SGData['SGGray']:
                    SpcGrp += " 1'"
                SGErr,SGData = G2spc.SpcGroup(SpcGrp)
                if '_' in BNSlatt:
                    SGData['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
                else:
                    SGData['BNSlattsym'] = [SGData['SGLatt'],[0.,0.,0.]]
                SGData['SGSpin'] = [1,]*len(SGData['SGSpin'])   #set to all black
                GenSym,GenFlg = G2spc.GetGenSym(SGData)[:2]
                SGData['GenSym'] = GenSym
                SGData['GenFlg'] = GenFlg
                SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
                G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
                generalData['SGData'] = SGData
                G2spc.UpdateSytSym(data)
                wx.CallAfter(UpdateGeneral)

            def OnShowSpins(event):
                msg = 'Magnetic space group information'
                text,table = G2spc.SGPrint(SGData,AddInv=not SGData['SGFixed'])
                text[0] = ' Magnetic Space Group: '+SGData['MagSpGrp']
                text[3] = ' The magnetic lattice point group is '+SGData['MagPtGp']
                if SGData['SGGray'] and "1'" not in text[0]:
                    text[0] += " 1'"
                    text[3] += "1'"
                G2G.SGMagSpinBox(General,msg,text,table,SGData['SGCen'],OprNames,
                    SGData['SpnFlp'],SGData['SGGray']& (not SGData['SGFixed'])).Show()

            SGData = generalData['SGData']
            GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
            if 'BNSlattsym' not in SGData:
                SGData['BNSlattsym'] = [SGData['SGLatt'],[0,0,0]]
            Indx = {}
            MagSym = SGData['MagSpGrp']
            if SGData['SGGray'] and "1'" not in MagSym:
                MagSym += " 1'"
            magSizer = wx.BoxSizer(wx.VERTICAL)
            magSizer.Add(wx.StaticText(General,label=' Magnetic spin operator selection:'),0)
            spinSizer = wx.BoxSizer(wx.HORIZONTAL)
            if SGData['SGFixed']:
                SpnFlp = SGData['SpnFlp']
                spinSizer.Add(wx.StaticText(General,label=' Magnetic phase from mcif file; no change in spin inversion allowed'),0,WACV)
                OprNames = G2spc.GenMagOps(SGData)[0]
            else:
                if not len(GenSym): # or SGData['SGGray']:
                    spinSizer.Add(wx.StaticText(General,label=' No spin inversion allowed'),0,WACV)
                    OprNames,SpnFlp = G2spc.GenMagOps(SGData)
                else:
                    spinSizer.Add(wx.StaticText(General,label=' BNS lattice: '),0,WACV)
                    BNSkeys = [SGData['SGLatt'],]+list(BNSsym.keys())
                    BNSkeys.sort()
                    try:        #this is an ugly kluge - bug in wx.ComboBox
                        if SGData['BNSlattsym'][0][2] in ['a','b','c']:
                            BNSkeys.reverse()
                    except:
                        pass
                    BNS = wx.ComboBox(General,
                        choices=BNSkeys,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    BNS.SetValue(SGData['BNSlattsym'][0])
                    BNS.Bind(wx.EVT_COMBOBOX,OnBNSlatt)
                    spinSizer.Add(BNS,0,WACV)
                    spinColor = ['black','red']
                    spCode = {-1:'red',1:'black'}
                    for isym,sym in enumerate(GenSym[1:]):
                        spinSizer.Add(wx.StaticText(General,label=' %s: '%(sym.strip())),0,WACV)
                        spinOp = wx.ComboBox(General,value=spCode[SGData['SGSpin'][isym+1]],choices=spinColor,
                            style=wx.CB_READONLY|wx.CB_DROPDOWN)
                        Indx[spinOp.GetId()] = isym
                        spinOp.Bind(wx.EVT_COMBOBOX,OnSpinOp)
                        spinSizer.Add(spinOp,0,WACV)
                    OprNames,SpnFlp = G2spc.GenMagOps(SGData)
                    SGData['SpnFlp'] = SpnFlp
            SGData['OprNames'] = OprNames
            magSizer.Add(spinSizer)
            msgSizer = wx.BoxSizer(wx.HORIZONTAL)
            msgSizer.Add(wx.StaticText(General,label=' Magnetic space group: %s  '%(MagSym)),0,WACV)
            showSpins = wx.Button(General,label=' Show spins?')
            showSpins.Bind(wx.EVT_BUTTON,OnShowSpins)
            msgSizer.Add(showSpins,0,WACV)
            magSizer.Add(msgSizer)
            dminSizer = wx.BoxSizer(wx.HORIZONTAL)
            dminSizer.Add(wx.StaticText(General,label=' Magnetic reflection d-min: '),0,WACV)
            dminVal = G2G.ValidatedTxtCtrl(General,generalData,'MagDmin',nDig=(10,4),xmin=0.7)
            dminSizer.Add(dminVal,0,WACV)
            magSizer.Add(dminSizer,0)
            return magSizer

        def ModulatedSizer(name):

            def OnShowSOps(event):
                SSGData = generalData['SSGData']
                text,table = G2spc.SSGPrint(generalData['SGData'],SSGData,not SGData['SGFixed'])
                msg = 'Superspace Group Information'
                G2G.SGMessageBox(General,msg,text,table,SGData.get('SpnFlp',[])).ShowModal()

            def OnSuperGp(event):   #for HKLF needs to reject SSgps not agreeing with modVec!
                'Respond to selection of a modulation group'
                wx.BeginBusyCursor()
                Choice = G2spc.SSChoice(SGData)
                wx.EndBusyCursor()
                dlg = wx.Dialog(General,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
                sizer = wx.BoxSizer(wx.VERTICAL)
                sizer.Add(wx.StaticText(dlg,wx.ID_OK,'Select or enter a modulation group'),0,wx.ALIGN_CENTER)
                sizer.Add((10,10))
                superGp = wx.ComboBox(dlg,value=generalData['SuperSg'],choices=Choice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
                sizer.Add(superGp)
                sizer.Add((10,10))
                btnsizer = wx.StdDialogButtonSizer()
                btn = wx.Button(dlg, wx.ID_OK)
                btn.SetDefault()
                btnsizer.AddButton(btn)
                btn = wx.Button(dlg, wx.ID_CANCEL)
                btnsizer.AddButton(btn)
                btnsizer.Realize()
                sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
                dlg.SetSizer(sizer)
                sizer.Fit(dlg)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    try:
                        SSymbol = superGp.GetValue().strip()
                    except AttributeError: # not sure why needed
                        SSymbol = superGp.GetLabel().strip()
                else:
                    dlg.Destroy()
                    return
                dlg.Destroy()
                try:
                    E,SSGData = G2spc.SSpcGroup(generalData['SGData'],SSymbol)
                except:
                    E = 'Invalid modulation group'
                    SSGData = None
                if SSGData:
                    Vec = generalData['SuperVec'][0]     #(3+1) only
                    modSymb = SSGData['modSymb']
                    generalData['SuperVec'][0] = G2spc.SSGModCheck(Vec,modSymb)[0]
                    generalData['SSGData'] = SSGData
                    generalData['SuperSg'] = SSymbol
                    OnShowSOps(event)
                else:
                    # needed in case someone manually enters an invalid SSG?
                    Text = '\n'.join([E+'\nSuperspace Group entry ignored'])
                    G2G.G2MessageBox(General,Text,'Superspace Group Error')
                wx.CallAfter(UpdateGeneral)

            def OnVecRef(event):
                generalData['SuperVec'][1] = Ref.GetValue()

            def OnMax(event):
                generalData['SuperVec'][2] = int(Max.GetValue())

            Indx = {}
            ssSizer = wx.BoxSizer(wx.VERTICAL)
            modSizer = wx.BoxSizer(wx.HORIZONTAL)
            modSizer.Add(wx.StaticText(General,label=' '+name.capitalize()+' structure controls: '),0,WACV)
            SGData = generalData['SGData']
            SpGrp = SGData.get('MagSpGrp',SGData['SpGrp'])
            if SGData['SGGray']:
                SpGrp += " 1'"
            modSizer.Add(wx.StaticText(General,label=' Superspace group: %s '%SpGrp),0,WACV)
            if not SGData['SGFixed']:
                val = generalData['SuperSg']
                superGp = wx.Button(General,wx.ID_ANY,val,size=(100,-1))
                superGp.Bind(wx.EVT_BUTTON,OnSuperGp)
            else:
                superGp = wx.StaticText(General,label=generalData['SuperSg'])
            modSizer.Add(superGp,0,WACV)
            modSizer.Add((5,5),0)
            showOps = wx.Button(General,label='Show ops.')
            showOps.Bind(wx.EVT_BUTTON,OnShowSOps)
            modSizer.Add(showOps,0,WACV)
            if PWDR:
                modSizer.Add(wx.StaticText(General,label=' Max index: '),0,WACV)
                indChoice = ['1','2','3','4','5','6','7']
                if 'Magnetic' in name.capitalize():    #limit to one for now
                    indChoice = ['1',]
                Max = wx.ComboBox(General,-1,value='%d'%(generalData['SuperVec'][2]),choices=indChoice,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Max.Bind(wx.EVT_COMBOBOX,OnMax)
                modSizer.Add(Max,0,WACV)
            ssSizer.Add(modSizer,0)
            vecSizer = wx.FlexGridSizer(1,5,5,5)
            vecSizer.Add(wx.StaticText(General,label=' Modulation vector: '),0,WACV)
            modS = G2spc.splitSSsym(generalData['SuperSg'])[0]
            generalData['SuperVec'][0],ifShow = G2spc.SSGModCheck(generalData['SuperVec'][0],modS)
            for i,[val,show] in enumerate(zip(generalData['SuperVec'][0],ifShow)):
                if show:
                    modVal = G2G.ValidatedTxtCtrl(General,generalData['SuperVec'][0],i,nDig=(10,4),xmin=-1.,xmax=2.)
                    vecSizer.Add(modVal,0,WACV)
                    Indx[modVal.GetId()] = i
                else:
                    modVal = G2G.ReadOnlyTextCtrl(General,value=('%.3f'%(val)),
                        size=wx.Size(50,20))
                    vecSizer.Add(modVal,0,WACV)
            if PWDR:
                Ref = wx.CheckBox(General,label='Refine?')
                Ref.SetValue(generalData['SuperVec'][1])
                Ref.Bind(wx.EVT_CHECKBOX, OnVecRef)
                vecSizer.Add(Ref,0,WACV)
            ssSizer.Add(vecSizer)
            return ssSizer

        def PawleySizer():
            # find d-space range in used histograms
            def enablePawley(*args):
                for c in PawleyCtrlsList:
                    c.Enable(generalData['doPawley'])
                # If Pawley is on, turn off Le Bail settings (since they will
                # be hidden)
                if generalData['doPawley']:
                    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.GPXtree.root, 'Controls'))
                    Controls['newLeBail'] = False
                    for h in data['Histograms']:
                        data['Histograms'][h]['LeBail'] = False
                    G2G.G2MessageBox(G2frame,title='Note:',
                            msg='Use Pawley Create in Operations menu of Pawley'+
                            ' Reflections tab to complete the Pawley setup')

            dmin,dmax,nhist,lbl,chist = getPawleydRange(G2frame,data)

            pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley controls: '),0,WACV)
            if nhist == 0:  # no data, no Pawley
                pawleySizer.Add(wx.StaticText(General,label='   no data'),0,WACV)
                generalData['doPawley'] = False
                return pawleySizer
            # force limits on dmin & dmax
            generalData['Pawley dmax'] = min(generalData['Pawley dmax'],dmax)
            generalData['Pawley dmin'] = max(generalData['Pawley dmin'],dmin)

            PawleyCtrlsList = []
            pawlRef = G2G.G2CheckBoxFrontLbl(General,' Do Pawley refinement?',
                                         generalData,'doPawley',enablePawley)
            pawleySizer.Add(pawlRef,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label='  dmin: '),0,WACV)
            pawlMin = G2G.ValidatedTxtCtrl(General,generalData,'Pawley dmin',size=(75,-1),
                xmin=dmin,xmax=20.,nDig=(10,5))
            PawleyCtrlsList.append(pawlMin)
            pawleySizer.Add(pawlMin,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label='  dmax: '),0,WACV)
            pawlMax = G2G.ValidatedTxtCtrl(General,generalData,'Pawley dmax',size=(75,-1),
                xmin=2.0,xmax=dmax,nDig=(10,5))
            PawleyCtrlsList.append(pawlMax)
            pawleySizer.Add(pawlMax,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley neg. wt.: '),0,WACV)
            pawlNegWt = G2G.ValidatedTxtCtrl(General,generalData,'Pawley neg wt',size=(65,-1),
                xmin=0.,xmax=1.,nDig=(10,3,'g'))
            PawleyCtrlsList.append(pawlNegWt)
            pawleySizer.Add(pawlNegWt,0,WACV)
            pawleyOuter = wx.BoxSizer(wx.VERTICAL)
            pawleyOuter.Add(pawleySizer)
            pawleyOuter.Add(wx.StaticText(General,label=lbl),0,wx.LEFT,120)
            for c in PawleyCtrlsList:
                c.Enable(generalData['doPawley'])
            return pawleyOuter

        def MapSizer():

            def OnMapType(event):
                Map['MapType'] = mapType.GetValue()
                if 'delt-F' in Map['MapType']:
                    data['Drawing']['contourColor'] = 'RdYlGn'
                else:
                    data['Drawing']['contourColor'] = 'YlGnBu'

            def OnRefList(event):
                if not refsList:
                    G2G.G2MessageBox(G2frame,'No reflections')
                    return
                dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select reflection sets to use',
                    'Use data',refsList)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        if not len(dlg.GetSelections()):
                            dlg.Destroy()
                            return
                        Map['RefList'] = [refsList[i] for i in dlg.GetSelections()]
                    else:
                        dlg.Destroy()
                        return
                finally:
                    dlg.Destroy()
                wx.CallAfter(UpdateGeneral,General.GetScrollPos(wx.VERTICAL))

            def OnDysnomia(event):
                data['General']['doDysnomia'] = not data['General']['doDysnomia']
                pages = [G2frame.phaseDisplay.GetPageText(PageNum) for PageNum in range(G2frame.phaseDisplay.GetPageCount())]
                if generalData['doDysnomia']:
                    if 'Dysnomia' not in pages:
                        G2frame.MEMData = wx.ScrolledWindow(G2frame.phaseDisplay)
                        # G2frame.Bind(wx.EVT_MENU, lambda event:G2phsG2.OnLoadDysnomia(event,G2frame,data), id=G2G.wxID_LOADDYSNOMIA)
                        # G2frame.Bind(wx.EVT_MENU, lambda event:G2phsG2.OnSaveDysnomia(event,G2frame,data), id=G2G.wxID_SAVEDYSNOMIA)
                        # G2frame.Bind(wx.EVT_MENU, lambda event:G2phsG2.OnRunDysnomia(event,G2frame,data), id=G2G.wxID_RUNDYSNOMIA)
                        G2frame.Bind(wx.EVT_MENU, OnLoadDysnomia, id=G2G.wxID_LOADDYSNOMIA)
                        G2frame.Bind(wx.EVT_MENU, OnSaveDysnomia, id=G2G.wxID_SAVEDYSNOMIA)
                        G2frame.Bind(wx.EVT_MENU, OnRunDysnomia, id=G2G.wxID_RUNDYSNOMIA)
                        G2frame.phaseDisplay.InsertPage(7,G2frame.MEMData,'Dysnomia')
                        Id = wx.NewId()
                        TabSelectionIdDict[Id] = 'Dysnomia'
                        if 'Dysnomia' not in data:  #set defaults here
                            data['Dysnomia'] = {'DenStart':'uniform','Optimize':'ZSPA','Lagrange':['user',0.001,0.05],
                                'wt pwr':0,'E_factor':1.,'Ncyc':5000,'prior':'uniform','Lam frac':[1,0,0,0,0,0,0,0],
                                'overlap':0.2,'MEMdmin':1.0}
                else:
                    if 'Dysnomia' in pages:
                        G2frame.phaseDisplay.DeletePage(pages.index('Dysnomia'))

            #patch
            if 'cutOff' not in Map:
                Map['cutOff'] = 100.0
            if 'Resolution' in Map:
                Map['GridStep'] = Map['Resolution']
            mapTypes = ['Fobs','Fcalc','delt-F','2*Fo-Fc','Omit','2Fo-Fc Omit','Patterson']
            refsList = [item for item in G2gd.GetGPXtreeDataNames(G2frame,['HKLF','PWDR']) if item in data['Histograms'].keys()]
            if not generalData['AtomTypes']:
                 mapTypes = ['Patterson',]
                 Map['MapType'] = 'Patterson'
            mapSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Fourier map controls: Map type: '),0,WACV)
            mapType = wx.ComboBox(General,value=Map['MapType'],choices=mapTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            mapType.Bind(wx.EVT_COMBOBOX,OnMapType)
            lineSizer.Add(mapType,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' Reflection sets: '),0,WACV)
            if 'list' not in str(type(Map['RefList'])):     #patch
                Map['RefList'] = [Map['RefList'],]
            lineSizer.Add(wx.ComboBox(General,value=Map['RefList'][0],choices=Map['RefList'],
                style=wx.CB_DROPDOWN|wx.CB_READONLY),0,WACV)
            refList = wx.Button(General,label='Select reflection sets')
            refList.Bind(wx.EVT_BUTTON,OnRefList)
            lineSizer.Add(refList,0,WACV)
            mapSizer.Add(lineSizer,0)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Map grid step: '),0,WACV)
            mapRes = G2G.ValidatedTxtCtrl(General,Map,'GridStep',nDig=(10,2),xmin=0.1,xmax=2.)
            line2Sizer.Add(mapRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' Peak cutoff %: '),0,WACV)
            cutOff = G2G.ValidatedTxtCtrl(General,Map,'cutOff',nDig=(10,1),xmin=1.0,xmax=100.)
            line2Sizer.Add(cutOff,0,WACV)
            if len(Map['RefList']) and not generalData['Modulated']:
                if all(['PWDR' in map for map in Map['RefList']]):
                    Dysno = wx.CheckBox(General,-1,label=' Use Dysnomia?')
                    Dysno.SetValue(generalData['doDysnomia'])
                    Dysno.Bind(wx.EVT_CHECKBOX,OnDysnomia)
                    line2Sizer.Add(Dysno,0,WACV)
                    hlpText = '''Dysnomia uses the maximum entropy method
                    to compute intensities for unobserved reflections.
                    '''
                    hlp = G2G.HelpButton(General,hlpText)
                    line2Sizer.Add(hlp,0,WACV)
            mapSizer.Add(line2Sizer,0)
            return mapSizer

        def FlipSizer():
            #patches
            if 'k-Max' not in Flip: Flip['k-Max'] = 20.
            if 'MScorr' not in Flip: Flip['MScorr'] = 0.0
            if 'Resolution' in Flip:
                Flip['GridStep'] = Flip['Resolution']

            def OnRefList(event):
                dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select reflection sets to use',
                    'Use data',refsList)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        if not len(dlg.GetSelections()):
                            dlg.Destroy()
                            return
                        Flip['RefList'] = [refsList[i] for i in dlg.GetSelections()]
                    else:
                        dlg.Destroy()
                        return
                finally:
                    dlg.Destroy()
                wx.CallAfter(UpdateGeneral,General.GetScrollPos(wx.VERTICAL))

            def OnNormElem(event):
                PE = G2elemGUI.PickElement(G2frame,ifNone=True)
                if PE.ShowModal() == wx.ID_OK:
                    Flip['Norm element'] = PE.Elem.strip()
                    normElem.SetLabel(Flip['Norm element'])
                PE.Destroy()

            def OnTestHKL(event):
                event.Skip()
                Obj = event.GetEventObject()
                name = Obj.GetName()
                try:
                    vals = Obj.GetValue().split()
                    Id = int(name.split('hkl')[1])
                    HKL = [int(val) for val in vals]
                    Flip['testHKL'][Id] = HKL
                except ValueError:
                    HKL = Flip['testHKL'][Id]
                Obj.ChangeValue('%3d %3d %3d'%(HKL[0],HKL[1],HKL[2]))

            refsList = [item for item in G2gd.GetGPXtreeDataNames(G2frame,['HKLF','PWDR']) if item in data['Histograms'].keys()]
            flipSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Charge flip controls: Reflection sets: '),0,WACV)
            if 'list' not in str(type(Flip['RefList'])):     #patch
                Flip['RefList'] = [Flip['RefList'],]
            refName = Flip['RefList'][0]
            refType = ''
            if refName and refName in data['Histograms']:
                refType = data['Histograms'][refName].get('Type','')
            lineSizer.Add(wx.ComboBox(General,value=Flip['RefList'][0],choices=Flip['RefList'],
                style=wx.CB_DROPDOWN|wx.CB_READONLY),0,WACV)
            refList = wx.Button(General,label='Select reflection sets')
            refList.Bind(wx.EVT_BUTTON,OnRefList)
            lineSizer.Add(refList,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' Normalizing element: '),0,WACV)
            normElem = wx.Button(General,label=Flip['Norm element'],style=wx.TE_READONLY)
            normElem.Bind(wx.EVT_BUTTON,OnNormElem)
            lineSizer.Add(normElem,0,WACV)
            flipSizer.Add(lineSizer,0)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Map grid step: '),0,WACV)
            flipRes = G2G.ValidatedTxtCtrl(General,Flip,'GridStep',nDig=(10,2),xmin=0.10,xmax=2.)
            line2Sizer.Add(flipRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Factor (0.01-1.2): '),0,WACV)
            kFactor = G2G.ValidatedTxtCtrl(General,Flip,'k-factor',nDig=(10,3),xmin=0.01,xmax=1.2)
            line2Sizer.Add(kFactor,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Max (>=10.0): '),0,WACV)
            kMax = G2G.ValidatedTxtCtrl(General,Flip,'k-Max',nDig=(10,1),xmin=10.)
            line2Sizer.Add(kMax,0,WACV)
            if refType == 'SEC':
                line2Sizer.Add(wx.StaticText(General,label=' MScorr (0-0.1): '),0,WACV)
                line2Sizer.Add(G2G.ValidatedTxtCtrl(General,Flip,'MScorr',nDig=(10,4),xmin=0.,xmax=0.10),0,WACV)
            flipSizer.Add(line2Sizer,0)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line3Sizer.Add(wx.StaticText(General,label=' Test HKLs:'),0,WACV)
            if len(Flip['testHKL']) < 5:
                Flip['testHKL'] += [[1,1,1],[0,2,0],[1,2,3]]
            HKL = Flip['testHKL']
            for ih,hkl in enumerate(Flip['testHKL']):
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                hkl = wx.TextCtrl(General,value='%3d %3d %3d'%(HKL[ih][0],HKL[ih][1],HKL[ih][2]),
                    style=wx.TE_PROCESS_ENTER,name='hkl%d'%(ih))
                hkl.Bind(wx.EVT_TEXT_ENTER,OnTestHKL)
                hkl.Bind(wx.EVT_KILL_FOCUS,OnTestHKL)
                line3Sizer.Add(hkl,0,WACV)
            flipSizer.Add(line3Sizer)
            return flipSizer

        def MCSASizer():

            def OnRefList(event):
                MCSAdata['Data source'] = refList.GetValue()

            def OnCycles(event):
                MCSAdata['Cycles'] = int(cycles.GetValue())

            def OnAlist(event):
                MCSAdata['Algorithm'] = Alist.GetValue()
                OnShowTsched()
                wx.CallAfter(UpdateGeneral,General.GetScrollPos(wx.VERTICAL))

            def OnRanStart(event):
                MCSAdata['ranStart'] = ranStart.GetValue()

#            def OnAutoRan(event):
#                MCSAdata['autoRan'] = autoRan.GetValue()

            def OnAnneal(event):
                event.Skip()
                Obj = event.GetEventObject()
                ind,fmt = Indx[Obj.GetId()]
                if ind == 2:        #No. trials
                    try:
                        val = int(Obj.GetValue())
                        if 1 <= val:
                            MCSAdata['Annealing'][ind] = val
                    except ValueError:
                        Obj.SetValue(fmt%(MCSAdata['Annealing'][ind]))
                else:
                    try:
                        val = float(Obj.GetValue())
                        if .0 <= val:
                            MCSAdata['Annealing'][ind] = val
                        Obj.SetValue(fmt%(MCSAdata['Annealing'][ind]))
                    except ValueError:
                        MCSAdata['Annealing'][ind] = None
                        Obj.SetValue(str(MCSAdata['Annealing'][ind]))

            def ShowTsched(invalid,value,tc):
                OnShowTsched()

            def OnShowTsched():
                if MCSAdata['Algorithm'] in ['fast','log']:
                    Y = G2mth.makeTsched(MCSAdata)
                    XY = [np.arange(len(Y)),np.log10(Y)]
                    G2plt.PlotXY(G2frame,[XY,],labelX='T-step',labelY='log(T)',newPlot=True,lines=True,Title='Annealing schedule')

#            OnShowTsched()
            refList = []
            if len(data['Pawley ref']):
                refList = ['Pawley reflections',]
            refList += [item for item in G2gd.GetGPXtreeDataNames(G2frame,['HKLF','PWDR']) if item in data['Histograms'].keys()]
            mcsaSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Monte Carlo/Simulated Annealing controls: Reflection set from: '),0,WACV)
            refList = wx.ComboBox(General,-1,value=MCSAdata['Data source'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' d-min: '),0,WACV)
            dmin = G2G.ValidatedTxtCtrl(General,MCSAdata,'dmin',nDig=(10,3),xmin=1.,xmax=5.)
            lineSizer.Add(dmin,0,WACV)
            mcsaSizer.Add(lineSizer)
            mcsaSizer.Add((5,5),)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' MC/SA runs: '),0,WACV)
            Cchoice = [str(2**i) for i in range(13)]
            cycles = wx.ComboBox(General,-1,value=str(MCSAdata.get('Cycles',1)),choices=Cchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            cycles.Bind(wx.EVT_COMBOBOX,OnCycles)
            line2Sizer.Add(cycles,0,WACV)
            line2Sizer.Add((5,0),)
            ranStart = wx.CheckBox(General,-1,label=' MC/SA Refine at ')
            ranStart.Bind(wx.EVT_CHECKBOX, OnRanStart)
            ranStart.SetValue(MCSAdata.get('ranStart',False))
            line2Sizer.Add(ranStart,0,WACV)
            MCSAdata['ranRange'] = MCSAdata.get('ranRange',10.)  #patch for old gpx files
            ranRange = G2G.ValidatedTxtCtrl(General,MCSAdata,'ranRange',nDig=(10,1),xmin=1.,xmax=99.)
            line2Sizer.Add(ranRange,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label='% of ranges. '),0,WACV)
            mcsaSizer.Add(line2Sizer)
            mcsaSizer.Add((5,5),)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            Achoice = ['log','fast','Basin Hopping']                #these work not 'boltzmann','cauchy',
            line3Sizer.Add(wx.StaticText(General,label=' MC/SA schedule: '),0,WACV)
            Alist = wx.ComboBox(General,-1,value=MCSAdata['Algorithm'],choices=Achoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Alist.Bind(wx.EVT_COMBOBOX,OnAlist)
            line3Sizer.Add(Alist,0,WACV)
            if MCSAdata['Algorithm'] in ['fast',]:
                Names = [' quench: ',' c-factor: ']
                parms = 'fast parms'
                for i,name in enumerate(Names):
                    line3Sizer.Add(wx.StaticText(General,label=name),0,WACV)
                    Ajump = G2G.ValidatedTxtCtrl(General,MCSAdata[parms],i,nDig=(10,2),xmin=0.1,xmax=1.,OnLeave=ShowTsched)
                    line3Sizer.Add(Ajump,0,WACV)
            elif 'log' in MCSAdata['Algorithm']:
                line3Sizer.Add(wx.StaticText(General,label=' slope: '),0,WACV)
                slope = G2G.ValidatedTxtCtrl(General,MCSAdata,'log slope',nDig=(10,3),xmin=0.25,xmax=1.0,OnLeave=ShowTsched)
                line3Sizer.Add(slope,0,WACV)
            elif 'Basin Hopping' in MCSAdata['Algorithm']:
                pass        #TODO basinhopping controls here?
            mcsaSizer.Add(line3Sizer)
            mcsaSizer.Add((5,5),)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line3Sizer.Add(wx.StaticText(General,label=' Annealing schedule: '),0,WACV)
            if 'Basin Hopping' in MCSAdata['Algorithm']:
                line3Sizer.Add(wx.StaticText(General,label=' Test temp: '),0,WACV)
                line3Sizer.Add(G2G.ValidatedTxtCtrl(General,MCSAdata['Annealing'],0,nDig=(10,5)),0,WACV)
            else:
                line3Sizer.Add(wx.StaticText(General,label=' Start temp: '),0,WACV)
                line3Sizer.Add(G2G.ValidatedTxtCtrl(General,MCSAdata['Annealing'],0,nDig=(10,5),OnLeave=ShowTsched),0,WACV)
                line3Sizer.Add(wx.StaticText(General,label=' Final temp: '),0,WACV)
                line3Sizer.Add(G2G.ValidatedTxtCtrl(General,MCSAdata['Annealing'],1,nDig=(10,5),OnLeave=ShowTsched),0,WACV)
            line3Sizer.Add(wx.StaticText(General,label=' No. trials: '),0,WACV)
            line3Sizer.Add(G2G.ValidatedTxtCtrl(General,MCSAdata['Annealing'],2),0,WACV)
            mcsaSizer.Add(line3Sizer)
            return mcsaSizer

        def compareSizer():

            def OnOatmOsel(event):
                generalData['Compare']['Oatoms'] = oatmsel.GetStringSelection()

            def OnTatmOsel(event):
                generalData['Compare']['Tatoms'] = tatmsel.GetStringSelection()

            def OnCompPlots(event):
                pName = generalData['Name']
                Oatoms = generalData['Compare']['Oatoms']
                Tatoms = generalData['Compare']['Tatoms']
                bName = '%s-%s'%(Oatoms,Tatoms)
                try:
                    Bonds = generalData['Compare']['Bonds'][bName]
                except KeyError:
                    print('need to do Compare first for %s polyhedra plots'%bName)
                    return
                Tilts = generalData['Compare']['Tilts'][bName]
                Vects = generalData['Compare']['Vects'][bName]
                dVects = generalData['Compare']['dVects'][bName]
                if len(Bonds['Obonds']):
                    print(' Octahedra:')
                    Bonds['Obonds'] = np.array(Bonds['Obonds'])
                    Bmean = np.mean(Bonds['Obonds'])
                    Bstd = np.std(Bonds['Obonds'])
                    title = '%s-%s Octahedral bond lengths'%(Oatoms,Tatoms)
                    G2plt.PlotBarGraph(G2frame,Bonds['Obonds'],Xname=r'$Bond, \AA$',Title=title,
                        PlotName='Oct %s Bond for %s'%(bName,pName))
                    Tilts['Otilts'] = np.array(Tilts['Otilts'])
                    Tmean = np.mean(Tilts['Otilts'])
                    Tstd = np.std(Tilts['Otilts'])
                    G2plt.PlotBarGraph(G2frame,Tilts['Otilts'],Xname='Tilts, deg',
                        Title='Octahedral %s tilts'%Oatoms,PlotName='Oct %s Tilts for %s'%(bName,pName))
                    dVects['Ovec'] = np.reshape(np.array(dVects['Ovec']),(-1,3))
                    for ix,aX in enumerate(['X','Y','Z']):
                        G2plt.PlotBarGraph(G2frame,dVects['Ovec'].T[ix],Xname=r'$%s%s, \AA$'%(GkDelta,aX),
                            Title='%s Octahedral distortion'%Oatoms,PlotName='Oct %s %s-Delta for %s'%(bName,aX,pName))
                    Vects['Ovec'] = np.array(Vects['Ovec'])                    #3D plot of tilt vectors
                    X = Vects['Ovec'].T[0]
                    Y = Vects['Ovec'].T[1]
                    Z = Vects['Ovec'].T[2]
                    R = Tilts['Otilts']
                    G2plt.PlotXYZvect(G2frame,X,Y,Z,R,r'X-axis',r'Y-axis',r'Z-axis',
                        Title=r'%s Octahedral tilt vectors'%Oatoms,PlotName='Oct %s tilts for %s'%(bName,pName))
                    print(' %s-%s bond distance: %.3f(%d)'%(Oatoms,Tatoms,Bmean,Bstd*1000))
                    print(' %s tilt angle: %.2f(%d)'%(Oatoms,Tmean,Tstd*100))

                if len(Bonds['Tbonds']):
                    print('Tetrahedra:')
                    Bonds['Tbonds'] = np.array(Bonds['Tbonds'])
                    Bmean = np.mean(Bonds['Tbonds'])
                    Bstd = np.std(Bonds['Tbonds'])
                    title = '%s-%s Tetrahedral bond lengths'%(Oatoms,Tatoms)
                    G2plt.PlotBarGraph(G2frame,Bonds['Tbonds'],Xname=r'$Bond, \AA$',Title=title,
                        PlotName='Tet %s Bond for %s'%(bName,pName))
                    Tilts['Ttilts'] = np.array(Tilts['Ttilts'])
                    Tmean = np.mean(Tilts['Ttilts'])
                    Tstd = np.std(Tilts['Ttilts'])
                    G2plt.PlotBarGraph(G2frame,Tilts['Ttilts'],Xname='Tilts, deg',
                        Title='Tetrahedral %s tilts'%Oatoms,PlotName='Tet %s Tilts for %s'%(bName,pName))
                    dVects['Tvec'] = np.reshape(np.array(dVects['Tvec']),(-1,3))
                    for ix,aX in enumerate(['X','Y','Z']):
                        G2plt.PlotBarGraph(G2frame,dVects['Tvec'].T[ix],Xname=r'$%s%s, \AA$'%(GkDelta,aX),
                            Title='%s Tetrahedral distortion'%Oatoms,PlotName='Tet %s %s-Delta for %s'%(bName,aX,pName))
                    Vects['Tvec'] = np.array(Vects['Tvec'])
                    X = Vects['Tvec'].T[0]
                    Y = Vects['Tvec'].T[1]
                    Z = Vects['Tvec'].T[2]
                    R = Tilts['Ttilts']
                    G2plt.PlotXYZvect(G2frame,X,Y,Z,R,r'X-axis',r'Y-axis',r'Z-axis',
                        Title=r'%s Tetrahedral tilt vectors'%Oatoms,PlotName='Tet %s tilts for %s'%(bName,pName))
                    print(' %s-%s bond distance: %.3f(%d)'%(Oatoms,Tatoms,Bmean,Bstd*1000))
                    print(' %s tilt angle: %.2f(%d)'%(Oatoms,Tmean,Tstd*100))
                    
            def OnSaveDist(event):
                pName = generalData['Name']
                fName = os.path.abspath(os.path.splitext(G2frame.GSASprojectfile)[0]+'_TOdist.csv')
                fp = open(fName,'w')
                fp.write('Phase name:  %s,'%pName)
                Oatoms = generalData['Compare']['Oatoms']
                Tatoms = generalData['Compare']['Tatoms']
                bName = '%s-%s'%(Oatoms,Tatoms)
                Bonds = generalData['Compare']['Bonds'][bName]
                # Vects = generalData['Compare']['Vects'][bName]
                # dVects = generalData['Compare']['dVects'][bName]
                if len(Bonds['Obonds']):
                    fp.write('\n%s-%s Octahedral bond lengths,\n'%(Oatoms,Tatoms))
                    for ib,bond in enumerate(Bonds['Obonds']):
                       fp.write('%7.4f,'%bond) 
                       if not (ib+1)%10: fp.write('\n')
                    
                if len(Bonds['Tbonds']):
                    fp.write('\n%s-%s Tetrahedral bond lengths,\n'%(Oatoms,Tatoms))
                    for ib,bond in enumerate(Bonds['Tbonds']):
                       fp.write('%7.4f,'%bond) 
                       if not (ib+1)%10: fp.write('\n')
                fp.close()
                print(' %s written'%fName)
                
            Oatoms = generalData['Compare']['Oatoms']
            Tatoms = generalData['Compare']['Tatoms']
            bName = '%s-%s'%(Oatoms,Tatoms)
            atTypes = generalData['AtomTypes']
            compSizer = wx.BoxSizer(wx.VERTICAL)
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add(wx.StaticText(General,label=' Compare polyhedra to ideal octahedra/tetrahedra:'),0,WACV)
            try:
                if len(generalData['Compare']['Bonds'][bName]['Obonds']) or len(generalData['Compare']['Bonds'][bName]['Tbonds']):
                    plotBtn = wx.Button(General,label='Show plots?')
                    plotBtn.Bind(wx.EVT_BUTTON,OnCompPlots)
                    topSizer.Add(plotBtn,0,WACV)
                    saveBtn = wx.Button(General,label='Save distributions?')
                    saveBtn.Bind(wx.EVT_BUTTON,OnSaveDist)
                    topSizer.Add(saveBtn,0,WACV)
            except KeyError:
                pass

            compSizer.Add(topSizer)
            atmselSizer = wx.BoxSizer(wx.HORIZONTAL)
            atmselSizer.Add(wx.StaticText(General,label=' Select origin atom type: '),0,WACV)
            oatmsel = wx.ComboBox(General,choices=atTypes,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            oatmsel.SetStringSelection(generalData['Compare']['Oatoms'])
            oatmsel.Bind(wx.EVT_COMBOBOX,OnOatmOsel)
            atmselSizer.Add(oatmsel,0,WACV)

            atmselSizer.Add(wx.StaticText(General,label=' Select target atom type: '),0,WACV)
            tatmsel = wx.ComboBox(General,choices=atTypes,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            tatmsel.SetStringSelection(generalData['Compare']['Tatoms'])
            tatmsel.Bind(wx.EVT_COMBOBOX,OnTatmOsel)
            atmselSizer.Add(tatmsel,0,WACV)

            atmselSizer.Add(wx.StaticText(General,label=' Sampling fraction: '),0,WACV)
            atmselSizer.Add(G2G.ValidatedTxtCtrl(General,generalData['Compare'],'Sampling',nDig=(8,3),xmin=0.0,xmax=1.0,),0,WACV)

            compSizer.Add(atmselSizer,0)
            return compSizer

        # UpdateGeneral execution starts here
        phaseTypes = ['nuclear','magnetic','macromolecular','faulted']
        SetupGeneral()
        generalData = data['General']
        # remove the Wave Data tab when present and not needed
        if not generalData['Modulated']:
            pages = [G2frame.phaseDisplay.GetPageText(PageNum) for PageNum in range(G2frame.phaseDisplay.GetPageCount())]
            if 'Wave Data' in pages:
                G2frame.phaseDisplay.DeletePage(pages.index('Wave Data'))
        Map = generalData['Map']
        Flip = generalData['Flip']
        MCSAdata = generalData['MCSA controls']
        PWDR = any(['PWDR' in item for item in data['Histograms'].keys()])
#patches
        if 'Pawley dmax' not in data['General']:
            data['General']['Pawley dmax'] = 100.0
        if 'SGFixed' not in data['General']['SGData']:
            data['General']['SGData']['SGFixed'] = False
        if 'SGGray' not in data['General']['SGData']:
            data['General']['SGData']['SGGray'] = False
        if 'Pawley ref' not in data:
            data['Pawley ref'] = []
#end patches
        if General.GetSizer():
            General.GetSizer().Clear(True)
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Overall params for {generalData['Name']!r} (type={data['General']['Type']})"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(NameSizer(),0,wx.EXPAND)
        mainSizer.Add((5,5),0)
        mainSizer.Add(CellSizer(),0)
        mainSizer.Add((5,5),0)

        Indx = {}
        denSizer = None
        if len(generalData['AtomTypes']):
            denSizer = DenSizer()
            mainSizer.Add(denSizer[0])
            mainSizer.Add((5,5),0)
            mainSizer.Add(ElemSizer())
        G2G.HorizontalLine(mainSizer,General)
        if generalData['Type'] == 'magnetic':
            if not generalData['SGData']['SGFixed']:
                GenSym,GenFlg,BNSsym = G2spc.GetGenSym(generalData['SGData'])
                generalData['SGData']['GenSym'] = GenSym
                generalData['SGData']['GenFlg'] = GenFlg
                generalData['SGData']['MagSpGrp'] = G2spc.MagSGSym(generalData['SGData'])
            mainSizer.Add(MagSizer())
            G2G.HorizontalLine(mainSizer,General)

        if generalData['Modulated']:
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_SINGLEMCSA,False)
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_MULTIMCSA,False)
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_4DCHARGEFLIP,True)
            mainSizer.Add(ModulatedSizer(generalData['Type']+' modulated'))
            G2G.HorizontalLine(mainSizer,General)
        else:
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_SINGLEMCSA,True)
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_MULTIMCSA,True)
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_4DCHARGEFLIP,False)
            
        dmin,dmax,nhist,lbl,chist = getPawleydRange(G2frame,data)
        if nhist > 0:
            mainSizer.Add(PawleySizer())
            G2G.HorizontalLine(mainSizer,General)
        
        mainSizer.Add(MapSizer())
        G2G.HorizontalLine(mainSizer,General)

        mainSizer.Add(FlipSizer())
        if generalData['Type'] in ['nuclear','macromolecular','faulted',]:
            G2G.HorizontalLine(mainSizer,General)
            mainSizer.Add(MCSASizer())
        G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_COMPARESTRUCTURE,False)
        if generalData['SGData']['SpGrp'] == 'P 1':
            G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_COMPARESTRUCTURE,True)
            G2G.HorizontalLine(mainSizer,General)
            mainSizer.Add(compareSizer())
        if SkipDraw:
            mainSizer.Clear(True)
            return
        G2frame.GetStatusBar().SetStatusText('',1)
        SetPhaseWindow(General,mainSizer,Scroll=Scroll)

    def OnTransform(event):
        Trans = np.eye(3)
        Uvec = np.zeros(3)
        Vvec = np.zeros(3)
        ifMag = False
        BNSlatt = ''
        ifRB = False
        RBModels = data['RBModels']
        for item in RBModels:
            if RBModels[item]:
                ifRB = True
        if ifRB:
            msg = 'Rigid Bodies present'
            text = 'Rigid bodies must be deleted before transform can be done'
            wx.MessageBox(text,caption=msg,style=wx.ICON_EXCLAMATION)
            return
        while True:
            dlg = TransformDialog(G2frame,data,Trans,Uvec,Vvec,ifMag,BNSlatt)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelection()
                    if result is None:
                        return
                    newPhase,Trans,Uvec,Vvec,ifMag,ifConstr,Common = result
                    newPhase['ranId'] = ran.randint(0,sys.maxsize)
                    SGData = newPhase['General']['SGData']
                    if ifMag:
                        BNSlatt = SGData['BNSlattsym'][0]

                    SGData['GenSym'],SGData['GenFlg'],BNSsym = G2spc.GetGenSym(SGData)
                    if '_' in BNSlatt:
                        SGData['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
                        G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
                    if SGData['SGGray']:
                        if SGData['SGInv']:
                            SGData['SpnFlp'] = np.concatenate((SGData['SpnFlp'],SGData['SpnFlp']))
                        SGData['SpnFlp'] = np.concatenate((SGData['SpnFlp'],-1*SGData['SpnFlp']))
                    SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
                    if not '_' in BNSlatt:
                        SGData['SGSpin'] = G2spc.GetSGSpin(SGData,SGData['MagSpGrp'])
                else:
                    return
            finally:
                dlg.Destroy()
            if 'setting' in Common:         #don't make new phase, Just move atoms!
                generalData = data['General']
                cx,ct,cs,cia = generalData['AtomPtrs']
                SGData = generalData['SGData']
                if SGData['SpGrp'] in G2spc.spg2origins:
                    T = G2spc.spg2origins[SGData['SpGrp']]
                Atoms = data['Atoms']
                for atom in Atoms:
                    for i in [0,1,2]:
                        atom[cx+i] += T[i]
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData) # update symmetry & mult
                    atom[cs:cs+2] = SytSym,Mul
                data['Drawing'] = []
                #force redraw of page
                sub = G2frame.GPXtree.GetSelection()
                G2frame.GPXtree.SelectItem(G2frame.GPXtree.root)
                wx.CallAfter(G2frame.GPXtree.SelectItem,sub)
                return
            else:
                phaseName = newPhase['General']['Name']
                newPhase,atCodes = G2lat.TransformPhase(data,newPhase,Trans,Uvec,Vvec,ifMag)
                detTrans = np.abs(nl.det(Trans))
                generalData = newPhase['General']
                SGData = generalData['SGData']
                SGData['fromParent'] = [Trans,Uvec,Vvec]        #save these
                Atoms = newPhase['Atoms']
                if ifMag:
                    atMxyz = []
                    for atom in Atoms:
                        if data['General']['Super']:
                            atom += [{'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}]
                        SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData)
                        CSI = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
                        MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                        atMxyz.append([MagSytSym,CSI[0]])
                    dlg = UseMagAtomDialog(G2frame,SGData['MagSpGrp'],Atoms,atCodes,atMxyz,ifDelete=False)
                    try:
                        opt = dlg.ShowModal()
                        if  opt == wx.ID_YES:
                            newPhase['Atoms'],atCodes = dlg.GetSelection()
                            generalData['Lande g'] = len(generalData['AtomTypes'])*[2.,]
                            break
                        else:
                            return
                    finally:
                        dlg.Destroy()
                else:
                    break

        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        UseList = newPhase['Histograms']
        for hist in UseList:
            UseList[hist]['Scale'] /= detTrans      #scale by 1/volume ratio
            if 'P' in UseList[hist]['Type']:
                UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
                UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
        newPhase['General']['Map'] = mapDefault.copy()
        sub = G2frame.GPXtree.AppendItem(parent=
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases'),text=phaseName)
        G2frame.GPXtree.SetItemPyData(sub,newPhase)
        newPhase['Drawing'] = []
        if 'RMC' in data: del newPhase['RMC']
        if ifConstr:
            G2cnstG.TransConstraints(G2frame,data,newPhase,Trans,Vvec,atCodes)     #data is old phase
        G2frame.GPXtree.SelectItem(sub)

    def OnISOSearch(event):
        '''Search for a higher symmetry structure consistent with the
        current phase using the ISOCIF web service
        '''
        def _showWebPage(event):
            'Show a web page when the user presses the "show" button'
            import tempfile
            txt = event.GetEventObject().page
            tmp = tempfile.NamedTemporaryFile(suffix='.html',
                        delete=False)
            with open(tmp.name,'w') as fp:
                fp.write(txt.replace(
                '<HEAD>',
                '<head><base href="https://iso.byu.edu/">',
                ))
            fileList.append(tmp.name)
            G2G.ShowWebPage('file://'+tmp.name,G2frame)
        import tempfile
        import re
        import requests
        from GSASII.exports import G2export_CIF
        isosite="https://iso.byu.edu/"
        upscript='isocifuploadfile.php'
        isoscript='isocifform.php'
        isoSubCite = ('For use of this supergroup search, please cite:\n'+
           G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...',wrap=60,indent=5))
        latTol,coordTol,occTol = 0.001, 0.01, 0.1
        oacomp,occomp = G2mth.phaseContents(data)
        ophsnam = data['General']['Name']
        fileList = []
        # write a CIF as a scratch file
        obj = G2export_CIF.ExportPhaseCIF(G2frame)
        obj.InitExport(None)
        obj.currentExportType='phase'
        obj.loadTree()
        tmp = tempfile.NamedTemporaryFile(suffix='.cif',delete=False)
        try:
            obj.dirname,obj.filename = os.path.split(tmp.name)
            obj.phasenam = data['General']['Name']
            obj.Writer('',data['General']['Name'])
            # isocif upload
            files = {'toProcess': open(tmp.name,'rb')}
            values = {'submit':'OK','input':'uploadcif'}
            r0 = requests.post(isosite+upscript, files=files, data=values)
            # get the filename out from the upload response
            for f in r0.text.split('INPUT')[1:]:
                if 'filename' in f: break
            else:
                G2G.G2MessageBox(G2frame,
                    f'ISOCIF upload problem\nHTML output={r0.text}',
                    'Upload Error')
                return
            isofile = f.split('VALUE')[1].split('"')[1]
        finally:
            os.unlink(tmp.name)

        # CIF uploaded, now process
        values1 = {'submit':'OK','input':'uploadcif','filename':isofile}
        r1 = requests.post(isosite+isoscript, data=values1)

        # from processed page, pull out form
        try:
            form = re.split('</form',
                [i for i in re.split('<form',r1.text,flags=re.IGNORECASE)
                     if 'Detect higher' in i][0],
                flags=re.IGNORECASE)[0]
        except IndexError:
            tmp1 = tempfile.NamedTemporaryFile(suffix='.html')
            with open(tmp1.name,'w') as fp:
                fp.write(r1.text.replace(
                '<HEAD>',
                '<head><base href="https://iso.byu.edu/">',
                ))
            G2G.ShowWebPage('file://'+tmp1.name,G2frame)
            G2G.G2MessageBox(G2frame,
                    'ISOCIF CIF processing problem, showing ISOCIF web page',
                    'ISOCIF Error')
            return

        formDict = {}
        for f in form.split('<'):
            try:
                name = re.split('name=',f,flags=re.IGNORECASE)[1].split('"')[1]
                value = re.split('value=',f,flags=re.IGNORECASE)[1].split('"')[1]
                formDict[name] = value
            except IndexError:
                pass

        repeat = True
        while repeat:
            repeat = False
            dlg = G2G.MultiDataDialog(G2frame,title='ISOCIF search',
                                prompts=['lattice constants tolerance',
                                     'coordinate tolerance',
                                     'occupancy tolerance'],
                                values=[latTol,coordTol,occTol],
                                limits=3*[[0.,2.]],formats=3*['%.5g'],
                                header=isoSubCite)
            res = dlg.ShowModal()
            latTol,coordTol,occTol = dlg.GetValues()
            dlg.Destroy()
            if res != wx.ID_OK:
                for i in fileList: os.unlink(i) # cleanup tmp web pages
                return
            formDict['acclat'], formDict['accpos'], formDict['accocc'] = latTol,coordTol,occTol
            r2 = requests.post(isosite+isoscript, data=formDict)
            # now parse out structural info from r2 and insert into current phase
            structure = r2.text.split('</H1>')[1].split('<p>')[2]
            atomList = []
            change = True
            for l in structure.split('\n')[1:]:
                l = l.split('<')[0]
                if 'No change in actual' in l:
                    change = False
                elif 'Space Group:' in l:
                    sgnum,sglbl,shname = l.split(':')[1].split()
                elif 'Lattice parameters:' in l:
                    cell = [float(i.split(',')[0].split('<')[0]) for i in l.split('=')[1:]]
                elif ':' in l or not l.strip():
                    pass
                else:  # must be an atom record
                    try:
                        nam = l.split()[0]
                        exec(l.split('),')[1].replace(',',';').strip())
                        xyz = [float(i) for i in eval(l.split('(')[1].split(')')[0])]
                        atomList.append([nam,xyz])
                    except:
                        if GSASIIpath.GetConfigValue('debug'):
                            print(f'could not parse "{l}"')
            if change:
                dlg = wx.Dialog(G2frame,wx.ID_ANY,'Supergroup found',
                                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
                mainSizer = wx.BoxSizer(wx.VERTICAL)
                dlg.SetSizer(mainSizer)
                msg = f'A supergroup structure in space group {sglbl} was found  with {len(atomList)} atoms in the asymmetric unit. '
                newPhase = makeIsoNewPhase(data,cell,atomList,sglbl,sgnum)
                nacomp,nccomp = G2mth.phaseContents(newPhase)
                msg += f"Unit cell {G2mth.fmtPhaseContents(nccomp)}"
                msg += f", vol={newPhase['General']['Cell'][7]:.2f} A^3"
                msg += f", density={G2mth.getDensity(newPhase['General'])[0]:.2f} g/cm^3."
                msg += f"\nAsymmetric unit {G2mth.fmtPhaseContents(nacomp)}."

                msg += f"\n\nOriginal structure is in space group {data['General']['SGData']['SpGrp']}"
                msg += f" and has {len(data['Atoms'])} atoms. "
                msg += f"Unit cell {G2mth.fmtPhaseContents(occomp)}"
                msg += f", vol={data['General']['Cell'][7]:.2f} A^3"
                msg += f", density={G2mth.getDensity(data['General'])[0]:.2f} g/cm^3."
                msg += f"\nAsymmetric unit {G2mth.fmtPhaseContents(oacomp)}."
                txt = wx.StaticText(dlg,wx.ID_ANY,msg)
                txt.Wrap(450)
                mainSizer.Add(txt)
                mainSizer.Add((-1,10))
                showSizer = wx.BoxSizer(wx.HORIZONTAL)
                btn = wx.Button(dlg, wx.ID_ANY,label='Show')
                btn.page = r2.text
                btn.Bind(wx.EVT_BUTTON,_showWebPage)
                showSizer.Add(btn)
                showSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                                ' Web page with ISOCIF search results'))
                mainSizer.Add(showSizer)
                mainSizer.Add((-1,5))
                G2G.HorizontalLine(mainSizer,dlg)
                mainSizer.Add((-1,10))
                mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                                                'Choose an action:'))
                mainSizer.Add((-1,10))
                btnsizer = wx.BoxSizer(wx.HORIZONTAL)
                btn = wx.Button(dlg, wx.ID_CANCEL, label="Quit")
                btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
                btnsizer.Add(btn)
                btnsizer.Add((5,5))
                btn = wx.Button(dlg, wx.ID_CLOSE, label="Change tolerances")
                btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CLOSE))
                btnsizer.Add(btn)
                btnsizer.Add((5,5))
                btn = wx.Button(dlg, wx.ID_OK, label="Create GSAS-II project\nwith this supergroup")
                btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
                btn.SetDefault()
                btnsizer.Add(btn)
                mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)

                mainSizer.Fit(dlg)
                dlg.CenterOnParent()
                ans = dlg.ShowModal()
                dlg.Destroy()
                if ans == wx.ID_CANCEL:
                    for i in fileList: os.unlink(i) # cleanup tmp web pages
                    return
                if ans == wx.ID_CLOSE:
                    repeat = True
                    continue
                break
            else:
                dlg = wx.MessageDialog(G2frame,'A higher symmetry structure was not found, repeat search with different tolerances?',
                                        'Repeat?',wx.YES|wx.NO)
                try:
                    dlg.CenterOnParent()
                    repeat = wx.ID_YES == dlg.ShowModal()
                    if not repeat:
                        for i in fileList: os.unlink(i) # cleanup tmp web pages
                        return
                finally:
                    dlg.Destroy()

        # now create a new .gpx file
        G2frame.OnFileSave(None) # save project on disk to restore to this later
        orgFilName = G2frame.GSASprojectfile
        # get restraints for later use (to clear them)
        resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        Restraints = G2frame.GPXtree.GetItemPyData(resId)
        resId = G2gd.GetGPXtreeItemId(G2frame,resId,ophsnam)
        cx,ct,cs,cia = data['General']['AtomPtrs']
        if ophsnam in Restraints:
            Restraints[ophsnam]['Bond']['Bonds'] = []
            Restraints[ophsnam]['Angle']['Angles'] = []
        f = saveIsoNewPhase(G2frame,data,newPhase,orgFilName)

        wx.MessageBox(f"{isoSubCite}\n\nCreated project in space group {sglbl} as file {f}",
                          style=wx.ICON_INFORMATION,caption='File created')
        for i in fileList: os.unlink(i) # cleanup tmp web pages
        def _GetPhase():
            'After search complete reload the project from the saved .gpx file'
            G2frame.OnFileOpen(None,filename=orgFilName,askSave=False)
            wx.CallLater(100,_ShowPhase)
        def _ShowPhase():
            'After search complete and project is reloaded, reopen tree to the original phase'
            phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            G2frame.GPXtree.Expand(phId)
            phId = G2gd.GetGPXtreeItemId(G2frame,phId,ophsnam)
            G2frame.GPXtree.SelectItem(phId)
        wx.CallLater(100,_GetPhase)

    def OnSuperSearch(event):
        '''Search for a supergroup matching the current phase using the
        Bilbao Pseudosymmetry search (PSEUDO) program
        '''
        def _showWebPage(event):
            'Show a web page when the user presses the "show" button'
            import tempfile
            num = event.GetEventObject().IndexNum
            tmp = tempfile.NamedTemporaryFile(suffix='.html',
                        delete=False)
            with open(tmp.name,'w') as fp:
                fp.write(pagelist[num].replace(
                '<head>',
                '<head><base href="https://www.cryst.ehu.es/">',
                ))
            fileList.append(tmp.name)
            G2G.ShowWebPage('file://'+tmp.name,G2frame)
        def _selectSuperGroups(rowdict,csdict,msg,depth=0,key=0):
            '''Present the user with a set of supergroups and allow
            selection of ones to be tested. Used initially with higher
            symmetry cells and again later to recheck structures that
            have been cast into spacegroups.
            '''
            width = 450
            dlg = wx.Dialog(G2frame,wx.ID_ANY,'Supergroup Choices',
                                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            dlg.SetSizer(mainSizer)
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                f'''Using the Bilbao Crystallographic Server Pseudosymmetry search (PSEUDO) 
program; Please cite:
{G2G.GetCite('Bilbao: PSEUDO',wrap=70,indent=5)}'''))
            mainSizer.Add((-1,5))
            G2G.HorizontalLine(mainSizer,dlg)
            txt = wx.StaticText(dlg,wx.ID_ANY,
                    'Searched for subgroups of model '+msg)
            txt.Wrap(width)
            mainSizer.Add(txt)
            if depth > 0:
                mainSizer.Add((-1,5))
                mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                    f'This is a {depth}-level supergroup of the original parent model'))
            showSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_ANY,label='Show')
            btn.IndexNum = key
            btn.Bind(wx.EVT_BUTTON,_showWebPage)
            showSizer.Add(btn)
            showSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                                ' Web page with supergroup search results'))
            mainSizer.Add(showSizer)
            mainSizer.Add((-1,5))
            G2G.HorizontalLine(mainSizer,dlg)
            mainSizer.Add((-1,5))
            txt = wx.StaticText(dlg,wx.ID_ANY,
                    'Select supergroups below to test with coordinate transformation.'+
                    ' Unselected options have been rejected due to lattice parameter'+
                    ' mismatch/value limits, but can be included.')
            txt.Wrap(width)
            mainSizer.Add(txt)
            mainSizer.Add((-1,10))
            sSizer = wx.BoxSizer(wx.HORIZONTAL)
            sSizer.Add((20,0))  # indent results
            spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width, 200))
            txtSizer = wx.FlexGridSizer(0, 3, 2, 2)
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Use'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Space group'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Unit Cell'))
            for key in rowdict:
                if key not in csdict: continue
                txtSizer.Add(G2G.G2CheckBox(spanel,key,csdict,key))
                cell = rowdict[key][5].split('\n')
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,rowdict[key][0]))
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,cell[0]))
            sSizer.Add(txtSizer)
            spanel.SetSizer(sSizer)
            mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)
            btnsizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_CANCEL, label="Quit")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
            btnsizer.Add(btn)
            btnsizer.Add((5,5))
            btn = wx.Button(dlg, wx.ID_CLOSE, label="Continue")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
            btnsizer.Add(btn)
            mainSizer.Add((-1,10))
            mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            dlg.SetSizer(mainSizer)
            mainSizer.Fit(dlg)
            spanel.SetAutoLayout(1)
            spanel.SetupScrolling()
            dlg.CenterOnParent()
            ans = dlg.ShowModal()
            return ans
        def _selectHiSymCell(rowdict,csdict):
            '''Present the user with a set of higher symmetry cells that are
            consistent with the starting cell. Used with monoclinic and triclinic
            starting cells only.
            '''
            width = 450
            dlg = wx.Dialog(G2frame,wx.ID_ANY,'Supercell Choices',
                                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            dlg.SetSizer(mainSizer)
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                    f'''Using the Bilbao Crystallographic Server Pseudosymmetry search (PSEUDO) 
program; Please cite:
{G2G.GetCite('Bilbao: PSEUDO',wrap=70,indent=5)}'''))
            mainSizer.Add((-1,5))
            G2G.HorizontalLine(mainSizer,dlg)
            showSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_ANY,label='Show')
            btn.IndexNum = 0
            btn.Bind(wx.EVT_BUTTON,_showWebPage)
            showSizer.Add(btn)
            showSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                                ' Web page with cell search results'))
            mainSizer.Add(showSizer)
            mainSizer.Add((-1,10))
            txt = wx.StaticText(dlg,wx.ID_ANY,
                    'Below select high symmetry cells to test with coordinate transformation.')
            txt.Wrap(width)
            mainSizer.Add(txt)
            G2G.HorizontalLine(mainSizer,dlg)
            mainSizer.Add((-1,10))
            sSizer = wx.BoxSizer(wx.HORIZONTAL)
            sSizer.Add((20,0))  # indent results
            spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width, 200))
            txtSizer = wx.FlexGridSizer(0, 5, 2, 20)
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Use'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Bravais\nLattice'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Unit Cell\n(as calc/symmetrized)'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Strain'))
            txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,'Tolerance'))
            for i,row in enumerate(rowdict):
                txtSizer.Add(G2G.G2CheckBox(spanel,str(i),csdict,i))
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,row[1]))
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,row[3]+'\n'+row[2]))
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,row[5]))
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,row[6]))
            sSizer.Add(txtSizer)
            spanel.SetSizer(sSizer)
            mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)
            btnsizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_CANCEL, label="Quit")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
            btnsizer.Add(btn)
            btnsizer.Add((5,5))
            btn = wx.Button(dlg, wx.ID_CLOSE, label="Continue")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
            btnsizer.Add(btn)
            mainSizer.Add((-1,10))
            mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            dlg.SetSizer(mainSizer)
            mainSizer.Fit(dlg)
            spanel.SetAutoLayout(1)
            spanel.SetupScrolling()
            dlg.CenterOnParent()
            ans = dlg.ShowModal()
            return ans
        def _GetPhase():
            'After search complete reload the project from the saved .gpx file'
            G2frame.OnFileOpen(None,filename=orgFilName,askSave=False)
            wx.CallLater(100,_ShowPhase)
        def _ShowPhase():
            'After search complete and project is reloaded, reopen tree to the original phase'
            phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            G2frame.GPXtree.Expand(phId)
            phId = G2gd.GetGPXtreeItemId(G2frame,phId,ophsnam)
            G2frame.GPXtree.SelectItem(phId)
        def _testSuperGroups(ophsnam,rowdict,csdict,valsdict,savedcookies,pagelist):
            'Use the Bilbao site to test selected supergroups'
            print(f'*** Testing {sum(csdict.values())} transformed structures')
            pgbar = wx.ProgressDialog('Supergroup Search',
                    f'Searching for supergroup(s) consistent with phase {ophsnam}',
                    len(csdict)+2,parent=G2frame,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
            try:
                pgbar.CenterOnParent()
                structDict = SUBGROUPS.BilbaoSymSearch2(valsdict,csdict,rowdict,savedcookies,
                                                            pagelist=pagelist,dlg=pgbar,ophsnam=ophsnam)
                GoOn = pgbar.Update(len(csdict)+2,newmsg=
                        f'Searching for supergroup(s) consistent with phase {ophsnam}'+
                            '\ndone')
                if not GoOn: return
                wx.GetApp().Yield()
            finally:
                pgbar.Destroy()
            return structDict
        def showSuperResults(G2frame,msgs,pagelist,fileList,ReSearch,parentpage,msg=None):
            '''Show a summary with info from a search of supergroups in
            :func:`OnSuperSearch` (in :func:`UpdatePhaseData`)
            '''
            def _showWebPage(event):
                import tempfile
                f = event.GetEventObject().webFile
                tmp = tempfile.NamedTemporaryFile(suffix='.html',
                                delete=False)
                with open(tmp.name,'w') as fp:
                    fp.write(f.replace(
                        '<head>',
                        '<head><base href="https://www.cryst.ehu.es/">',
                        ))
                fileList.append(tmp.name)
                G2G.ShowWebPage('file://'+tmp.name,G2frame)
            width = 500
            dlg = wx.Dialog(G2frame,wx.ID_ANY,'Search results',
                                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            dlg.SetSizer(mainSizer)
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                    f'''Using the Bilbao Crystallographic Server Pseudosymmetry search (PSEUDO) 
program; Please cite:
{G2G.GetCite('Bilbao: PSEUDO',wrap=70,indent=5)}'''))
            if msg:
                txt = wx.StaticText(dlg,wx.ID_ANY,'Starting from '+msg.replace('\n',' '))
            txt.Wrap(width)
            mainSizer.Add((-1,10))
            mainSizer.Add(txt)
            mainSizer.Add((-1,5))
            showSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_ANY,label='Show')
            btn.webFile = parentpage
            btn.Bind(wx.EVT_BUTTON,_showWebPage)
            showSizer.Add(btn)
            showSizer.Add(wx.StaticText(dlg,wx.ID_ANY,' Web page with supergroup search results'))
            mainSizer.Add(showSizer)
            mainSizer.Add((-1,10))
            sSizer = wx.BoxSizer(wx.HORIZONTAL)
            sSizer.Add((30,0))  # indent results
            spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width, 200))
            txtSizer = wx.BoxSizer(wx.VERTICAL)
            smsgs = sorted([i for i in pagelist if i !=0])
            if len(smsgs) == 0:  # prevent an empty panel
                txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,
                            '*** No higher symmetry candidate structures found ***'))
                txtSizer.Add((-1,10))
            for num in smsgs:
                if num in msgs:
                    if '@' not in num: G2G.HorizontalLine(txtSizer,spanel)
                    txt = wx.StaticText(spanel,wx.ID_ANY,f'Case {num}: ' + msgs[num])
                    txt.Wrap(width-50)
                    txtSizer.Add(txt)
                    if pagelist[num] is None:
                        txtSizer.Add((-1,10))
                        continue
                elif pagelist[num] is not None:
                    txt = wx.StaticText(spanel,wx.ID_ANY,
                                        f'Processing for case {num} incomplete')
                    txt.Wrap(width-50)
                    txtSizer.Add(txt)
                elif pagelist[num] is None:
                    txtSizer.Add((-1,5))
                    txtSizer.Add(wx.StaticText(spanel,wx.ID_ANY,
                                        f'Processing of case {num} failed'))
                    txtSizer.Add((-1,10))
                    continue
                txtSizer.Add((-1,5))
                showSizer = wx.BoxSizer(wx.HORIZONTAL)
                if '@' in num: showSizer.Add((20,-1))
                btn = wx.Button(spanel, wx.ID_ANY,label='Show')
                btn.webFile = pagelist[num]
                btn.Bind(wx.EVT_BUTTON,_showWebPage)
                showSizer.Add(btn)
                showSizer.Add(wx.StaticText(spanel,wx.ID_ANY,' Web page with transform info'))
                txtSizer.Add(showSizer)
                if num in ReSearch:
                    key = 'use_' + num
                    if key not in ReSearch: ReSearch[key] = True
                    reSizer = wx.BoxSizer(wx.HORIZONTAL)
                    reSizer.Add((20,-1))
                    redo = G2G.G2CheckBoxFrontLbl(spanel,
                                    ' Search again for supergroups of this result?',
                                    ReSearch,key)
                    reSizer.Add(redo)
                    txtSizer.Add(reSizer)
                if num in msgs and '@' not in num: G2G.HorizontalLine(txtSizer,spanel)

                txtSizer.Add((-1,10))
            sSizer.Add(txtSizer)
            spanel.SetSizer(sSizer)
            mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)
            btnsizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_CLOSE, label="Continue")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
            btnsizer.Add(btn)
            mainSizer.Add((-1,10))
            mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            dlg.SetSizer(mainSizer)
            mainSizer.Fit(dlg)
            spanel.SetAutoLayout(1)
            spanel.SetupScrolling()
            dlg.CenterOnParent()
            ans = dlg.ShowModal()
            dlg.Destroy()
            return ans

        def _showSummary(G2frame,msgs,gpxList):
            '''Summarize the final results from all steps'''

            width = 500
            dlg = wx.Dialog(G2frame,wx.ID_ANY,'Final Supergroup Search Results',
                                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            dlg.SetSizer(mainSizer)
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                f'''Using the Bilbao Crystallographic Server Pseudosymmetry search (PSEUDO) 
program; Please cite:
{G2G.GetCite('Bilbao: PSEUDO',wrap=70,indent=5)}'''))
            mainSizer.Add((-1,10))
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                f'From the starting model, {len(gpxList)} possible supergroups were located.'))
            if 0 in msgs:
                txt = wx.StaticText(dlg,wx.ID_ANY,msgs[0]
                                .replace('/volume:\n',' && volume:').replace('\n',' '))
                txt.Wrap(width-50)
                mainSizer.Add(txt)
            spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width, 200))
            txtSizer = wx.BoxSizer(wx.VERTICAL)
            G2G.HorizontalLine(mainSizer,dlg)
            mainSizer.Add((-1,4))
            for m in gpxList:
                if m == 0: continue
                msg = 'Found ' + m.replace('\n  after','. After').replace('/volume:\n',' && volume:')
                txt = wx.StaticText(spanel,wx.ID_ANY,msg)
                txt.Wrap(width-50)
                txtSizer.Add(txt)
                txtSizer.Add((-1,4))
                G2G.HorizontalLine(txtSizer,spanel)
            spanel.SetSizer(txtSizer)
            mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)
            btnsizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(dlg, wx.ID_CLOSE, label="Continue")
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
            btnsizer.Add(btn)
            mainSizer.Add((-1,10))
            mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            dlg.SetSizer(mainSizer)
            mainSizer.Fit(dlg)
            spanel.SetAutoLayout(1)
            spanel.SetupScrolling()
            dlg.CenterOnParent()
            ans = dlg.ShowModal()
            dlg.Destroy()
            return ans
        def fmtCell(cell):
            s = ''
            for i in cell[0:3]: s += f"{i:.3f}, "
            for i in cell[3:5]: s += f"{i:.2f}, "
            s += f"{cell[5]:.2f}"
            return s

        #### processing for OnSuperSearch starts here ####
        fileList = []
        ReSearch = {}
        gpxList = []
        ophsnam = data['General']['Name']
        pgbar = wx.ProgressDialog('Supergroup Search',
            f'Searching for supergroup(s) consistent with phase {ophsnam}',5,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT,
                parent=G2frame)
        pgbar.CenterOnParent()
        try:
            G2frame.OnFileSave(None) # save project on disk to restore to this later
            orgFilName = G2frame.GSASprojectfile
            # get restraints for later use (to clear them)
            resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
            Restraints = G2frame.GPXtree.GetItemPyData(resId)
            resId = G2gd.GetGPXtreeItemId(G2frame,resId,ophsnam)
            cx,ct,cs,cia = data['General']['AtomPtrs']
            # get starting unit cell contents
            oacomp,occomp = G2mth.phaseContents(data)
            msgs = {}
            msgs[0] = f"initial structure: cell = {fmtCell(data['General']['Cell'][1:7])}"
            msgs[0] += f", vol={data['General']['Cell'][7]:.2f} A^3"
            msgs[0] += f", Space group {data['General']['SGData']['SpGrp']}. "
            msgs[0] += f"Before any transform, unit cell {G2mth.fmtPhaseContents(occomp)}"
            msgs[0] += f", density={G2mth.getDensity(data['General'])[0]:.2f} g/cm^3"
            msgs[0] += f". Asymmetric unit {G2mth.fmtPhaseContents(oacomp)} ({len(data['Atoms'])} atoms)."
            startSet = msgs[0]

            # fix non-standard space group settings
            GoOn = pgbar.Update(1,newmsg=
                f'Searching for supergroup(s) consistent with phase {ophsnam}'+
                '\nTesting structure for standard setting')
            if not GoOn: return
            wx.GetApp().Yield()
            # need to convert non-standard space group settings
            print('*** Checking space group setting')
            sgnum,sgsym,xmat,xoff = SUBGROUPS.GetStdSGset(data['General']['SGData'])
            newPhase = copy.deepcopy(data)
            try:
                if np.allclose(np.eye(3),xmat) and np.allclose(xoff,np.zeros_like(xoff)):
                    print('*** Structure in standard setting')
                else:
                    print('*** Transforming structure to standard setting')
                    newPhase['ranId'] = ran.randint(0,sys.maxsize),
                    newPhase['General']['SGData'] = G2spc.SpcGroup(sgsym)[1]
                    newPhase['General']['Cell'][1:] = G2lat.TransformCell(
                        newPhase['General']['Cell'][1:-1],xmat)
                    uvec = np.array(xoff)
                    vvec = np.array([0.,0.,0.])
                    newPhase['MagXform'] = (xmat,xoff,vvec)
                    newPhase,atCodes = G2lat.TransformPhase(
                        data,newPhase,xmat,uvec,vvec,False)
                    startSet = "transformed starting structure: cell = "
                    startSet += fmtCell(newPhase['General']['Cell'][1:7])
                    startSet += f". Space group {newPhase['General']['SGData']['SpGrp']}."
            except:
                G2G.G2MessageBox(G2frame,
                        'Standard setting check failed. Try again later.',
                        'Unexpected error')
                return

            # search from a standard space group setting
            GoOn = pgbar.Update(2,newmsg=
                f'Searching for supergroup(s) consistent with phase {ophsnam}'+
                '\nSearching with phase')
            if not GoOn: return
            wx.GetApp().Yield()
            pagelist = {}
            valsdict,csdict,rowdict,savedcookies = SUBGROUPS.BilbaoSymSearch1(
                sgnum,newPhase,pagelist=pagelist)
        finally:
            pgbar.Destroy()

        # process initial PSEUDO results
        if csdict is None and len(rowdict) == 0:   # this was monoclinic or triclinic
            # look for supergroups of the current cell
            pgbar = wx.ProgressDialog('Supergroup Search',
                    f'Searching for supergroup(s) consistent with phase {ophsnam}',
                    1+len(rowdict),
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE,
                    parent=G2frame)
            try:
                pgbar.CenterOnParent()
                wx.GetApp().Yield()
                valsdict,csdict,rowdict,savedcookies = SUBGROUPS.BilbaoSymSearch1(
                    sgnum,newPhase,pagelist=pagelist,keepCell=True)
            finally:
                pgbar.Destroy()
            ans = _selectSuperGroups(rowdict,csdict,'from '+startSet+
                        '\n*** Note, no higher symmetry cells found.')
            if ans == wx.ID_CANCEL: return
            structDict = _testSuperGroups(ophsnam,rowdict,csdict,valsdict,savedcookies,pagelist)
            if len(structDict) != 0: ReSearch = SUBGROUPS.find2SearchAgain(pagelist,'')
        elif csdict is None:   # this was monoclinic or triclinic
            structDict = {}
            csdict = len(rowdict)*[True]
            ans = _selectHiSymCell(rowdict,csdict)
            if ans == wx.ID_CANCEL: return
            pgbar = wx.ProgressDialog('Supergroup Search',
                    f'Searching for supergroup(s) consistent with phase {ophsnam}',
                    1+len(rowdict),
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT,
                    parent=G2frame)
            try:
                pgbar.CenterOnParent()
                wx.GetApp().Yield()

                for i,row in enumerate(rowdict):
                    if not csdict[i]: continue
                    GoOn = pgbar.Update(i,newmsg=
        f'Searching for supergroup(s) w/cell {str(row[2])}\nLattice {row[1]} -- starting')
                    wx.GetApp().Yield()
                    lbl,latticeList,vals1dict,rowList = SUBGROUPS.BilbaoLowSymSea1(
                    valsdict,row,savedcookies,pagelist=pagelist)
                    msgs[lbl] = (f'Using cell {str(row[2])} with lattice {row[1]}'
                                + f'. Checking {len([i for i in rowList if i[0]])}'
                                + f' supergroups (of {len(rowList)} found)')
                    for row1 in rowList:
                        if not row1[0]: continue
                        GoOn = pgbar.Update(i,newmsg=
                                            f'Searching for supergroup(s) w/cell {str(row[2])}'
                                            + f'\nLattice {row[1]} && spacegroup {row1[2]}'
                                            )
                        wx.GetApp().Yield()
                        lbl1,structure = SUBGROUPS.BilbaoLowSymSea2(
                                                row[0],vals1dict,row1,savedcookies,
                                                pagelist=pagelist)
                        if structure is not None:
                            structDict[lbl1] = structure
                        else:
                            msgs[lbl1] = f'Coordinates inconsistent with space group {row1[2]}'
            finally:
                pgbar.Destroy()
            if len(structDict) != 0: ReSearch = SUBGROUPS.find2SearchAgain(pagelist)
        else: # not monoclinic or triclinic
            ans = _selectSuperGroups(rowdict,csdict,'from '+startSet)
            if ans == wx.ID_CANCEL: return
            structDict = _testSuperGroups(ophsnam,rowdict,csdict,valsdict,savedcookies,pagelist)
            if len(structDict) != 0: ReSearch = SUBGROUPS.find2SearchAgain(pagelist,'')

        # searches completed.
        if len(structDict) != 0:  # were new structures generated?
            # new phases will need different restraints clear them (probably
            # should clear constraints too)
            if ophsnam in Restraints:
                Restraints[ophsnam]['Bond']['Bonds'] = []
                Restraints[ophsnam]['Angle']['Angles'] = []
            # Now generate .gpx files and show results
            for num,s in structDict.items():   # loop over supergroup settings
                f = G2IO.saveNewPhase(G2frame,data,s,num,msgs,orgFilName)
                if f: gpxList.append(msgs[num])
        ans = showSuperResults(G2frame,msgs,pagelist,fileList,ReSearch,pagelist[0],msgs[0])
        for i in fileList: os.unlink(i) # cleanup tmp web pages
        fileList = []

        # repeat search on any identified (& selected) supergroups
        repeatcount = 0
        while ReSearch:
            repeatcount += 1
            NextSearch = {}
            for key in ReSearch:
                #print(key,'in ReSearch')
                pagelist = {}
                if key.startswith('use_'): continue
                if not ReSearch.get('use_'+key,False): continue
                fromMsg = msgs[key]
                del msgs[key]
                # need a status bar here
                if GSASIIpath.GetConfigValue('debug'): print(f"processing {key}")
                pgbar = wx.ProgressDialog('Supergroup Search',
                    f'Searching for supergroup(s) from case {key}',
                    1,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE,
                    parent=G2frame)
                try:
                    pgbar.CenterOnParent()
                    wx.GetApp().Yield()
                    valsdict,csdict,rowdict,savedcookies = SUBGROUPS.BilbaoReSymSearch(
                        key,ReSearch[key],pagelist=pagelist)
                finally:
                    pgbar.Destroy()
                ans = _selectSuperGroups(rowdict,csdict,
                    f'case {key}, {fromMsg}',repeatcount,key=key)
                parentpage = pagelist[key]
                del pagelist[key]
                structDict = _testSuperGroups(ophsnam,rowdict,csdict,valsdict,savedcookies,pagelist)
                for num,s in structDict.items():   # loop over supergroup settings
                    f = G2IO.saveNewPhase(G2frame,data,s,num,msgs,orgFilName)
                    if f:
                        gpxList.append(msgs[num])
                fndStruct = SUBGROUPS.find2SearchAgain(pagelist,'')
                if not fndStruct: continue
                ans = showSuperResults(G2frame,msgs,pagelist,fileList,fndStruct,parentpage,fromMsg)
                # rename the msg & structure entry to have a reference to the parent
                for k in fndStruct:
                    if k.startswith('use_'): continue
                    if not fndStruct.get('use_'+k,False): continue
                    nkey = key + '_' + k
                    NextSearch[nkey] = fndStruct[k]
                    NextSearch['use_'+nkey] = True
                    msgs[nkey] = msgs.pop(k)
            ReSearch = NextSearch

        for i in fileList: os.unlink(i) # cleanup tmp web pages

        # show final message
        if len(gpxList):
            _showSummary(G2frame, msgs, gpxList)
            print(f'Search done, from {msgs[0]}\n{len(gpxList)} supergroups located:\n')
            for i in gpxList: print(i)
        else:
            G2G.G2MessageBox(G2frame,
                    'No possible supergroups were found to match the starting model.',
                    'Search complete')
            print('Search done, no supergroups located')

        # Restore the original saved project
        wx.CallLater(100,_GetPhase)

    def OnSubSearch(event):
        '''Search for a lower symmetry structure consistent with the
        current phase using the Bilbao SUBGROUPS web service

        Based on GSASIIpwdGUI.OnRunSubs
        '''
        generalData = data['General']
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
        if dlg.ShowModal() != wx.ID_OK: return
            
        subcells = []
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
        wx.MessageBox(' For use of SUBGROUPS, please cite:\n\n'+
                          G2G.GetCite('Bilbao: k-SUBGROUPSMAG'),
                          caption='Bilbao SUBGROUPS',
                          style=wx.ICON_INFORMATION)
        wx.BeginBusyCursor()
        SGData = generalData['SGData']
        from . import SUBGROUPS as kSUB
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

        # emulate indexing controls 
        controls = 16*[None]
        controls[6:12] = generalData['Cell'][1:7] # cell
        cx,ct,cs,cia = data['General']['AtomPtrs']
        controls[15] = [atom[:cx+3] for atom in data['Atoms']]
        
        dlg = wx.ProgressDialog('SUBGROUPS results',f'Processing {SubGroups[0][0]}',len(SubGroups),
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)
        import GSASII.GSASIIpwdGUI
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
            phase['nAtoms'] = len(GSASII.GSASIIpwdGUI.TestAtoms(phase,controls[15],SGData,Uvec,Trans,maxequiv,maximal))
            subcells.append(phase)
        dlg.Destroy()
        data['SUBGROUPS'] = (subcells,baseList)
        msg = f'''{len(subcells)} subgroup entries were generated. Use menu commands 
"Select magnetic/subgroup phase" or "Make subgroup project file(s)"
to use these entries'''
        G2G.G2MessageBox(G2frame,msg,'Sugbroups stored')
        wx.CallAfter(UpdatePhaseData,G2frame,Item,data)

    def OnTransform2Std(event):
        '''Uses the Bilbao web site to transform a space group and coordinates
        to a standard setting
        '''
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        sgnum,sgsym,xmat,xoff = SUBGROUPS.GetStdSGset(generalData['SGData'])
        if np.allclose(np.eye(3),xmat) and np.allclose(xoff,np.zeros_like(xoff)):
            msg = "Nothing to do. Structure is already set in the standard setting."
            G2G.G2MessageBox(G2frame,msg,'No change needed')
            return
        G2frame.OnFileSave(None) # save
        orgFilName = G2frame.GSASprojectfile
        ophsnam = data['General']['Name']
        # get starting unit cell contents
        ocomp = {}
        for atom in data['Atoms']:
            if atom[ct] in ocomp:
                ocomp[atom[ct]] += atom[cs+1]
            else:
                ocomp[atom[ct]] = atom[cs+1]
        ovol = data['General']['Cell'][7]
        # get restraints & clear geometrical restraints
        resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        Restraints = G2frame.GPXtree.GetItemPyData(resId)
        resId = G2gd.GetGPXtreeItemId(G2frame,resId,ophsnam)
        if ophsnam in Restraints:
            Restraints[ophsnam]['Bond']['Bonds'] = []
            Restraints[ophsnam]['Angle']['Angles'] = []
        # create a new phase
        phaseName = f'{ophsnam}-std'
        newPhase = copy.deepcopy(data)
        newPhase['ranId'] = ran.randint(0,sys.maxsize),
        if 'magPhases' in data: del newPhase['magPhases']
        generalData = newPhase['General']
        generalData['Name'] = phaseName
        generalData['SGData'] = SGData = G2spc.SpcGroup(sgsym)[1]
        generalData['Cell'][1:] = G2lat.TransformCell(generalData['Cell'][1:-1],xmat)
        uvec = np.array(xoff)
        vvec = np.array([0.,0.,0.])
        newPhase['MagXform'] = (xmat,xoff,vvec)
        newPhase,atCodes = G2lat.TransformPhase(data,newPhase,xmat,uvec,vvec,False)
        Atoms = newPhase['Atoms']
        Atms = []
        AtCods = []
        atMxyz = []
        for ia,atom in enumerate(Atoms):
            atom[0] += '_%d'%ia
            atom[2] = ''                    #clear away refinement flags
            SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData)
            Atms.append(atom)
            AtCods.append(atCodes[ia])
            CSI = G2spc.GetCSxinel(SytSym)
            atMxyz.append([SytSym,CSI[0]])
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        UseList = newPhase['Histograms']
        detTrans = np.abs(nl.det(xmat))
        newPhase['Drawing'] = []
        for hist in UseList:
            UseList[hist]['Scale'] /= detTrans      #scale by 1/volume ratio
            # reset Dij & microstrain terms where # of terms changes
            if len(UseList[hist]['Mustrain'][4]) != NShkl:
                UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
            if len(UseList[hist]['HStrain'][0]) != NDij:
                UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
            newPhase['General']['Map'] = mapDefault.copy()
        # phase name rename
        newName = generalData['Name'] = phaseName
        phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
        #phaseNameList = usedHistograms.keys() # phase names in use
        generalData['Name'] = newName
        G2frame.GPXtree.SetItemText(Item,generalData['Name'])
        # change phase name key in Reflection Lists for each histogram
        for hist in data['Histograms']:
            ht = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
            rt = G2gd.GetGPXtreeItemId(G2frame,ht,'Reflection Lists')
            if not rt: continue
            RfList = G2frame.GPXtree.GetItemPyData(rt)
            RfList[newName] = []
            if ophsnam in RfList:
                del RfList[ophsnam]
        # copy restraints w/o the geometry-based ones
        if ophsnam in Restraints:
            Restraints[newName] = Restraints[ophsnam]
            del Restraints[ophsnam]
        if resId: G2frame.GPXtree.SetItemText(resId,newName)
        data.update(newPhase)
        # save new file
        G2frame.GSASprojectfile = os.path.splitext(orgFilName
                            )[0]+'_std.gpx'
        #G2frame.OnFileSaveas(event)
        G2IO.ProjFileSave(G2frame)
        # get transformed contents
        ncomp = {}
        for atom in data['Atoms']:
            if atom[ct] in ncomp:
                ncomp[atom[ct]] += atom[cs+1]
            else:
                ncomp[atom[ct]] = atom[cs+1]
        nvol = data['General']['Cell'][7]
        o = ', '.join([f'{i}:{j}' for i,j in ocomp.items()])
        n = ', '.join([f'{i}:{j}' for i,j in ncomp.items()])
        msg = f"""New project file, {G2frame.GSASprojectfile} created.

        Before & after transform, unit cell contents/volume:

        Before {o}, {ovol:.2f} A^3
        After  {n}, {nvol:.2f} A^3
        """
        G2G.G2MessageBox(G2frame,msg,'Transform complete')

        # Restore the original saved project
        G2frame.OnFileOpen(None,filename=orgFilName,askSave=False)
        # reopen tree to the original phase
        def _ShowPhase():
            phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            G2frame.GPXtree.Expand(phId)
            phId = G2gd.GetGPXtreeItemId(G2frame,phId,ophsnam)
            G2frame.GPXtree.SelectItem(phId)
        wx.CallLater(100,_ShowPhase)

    def OnCompareCells(event):
        G2G.Load2Cells(G2frame,data)

    def OnCompare(event):
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        atNames = [atm[ct-1] for atm in data['Atoms']]
        if not generalData['Compare']['Oatoms'] or not generalData['Compare']['Tatoms']:
            G2frame.ErrorDialog('Compare atom selection error','Select atoms for polygon comparison first')
            return
        bName = '%s-%s'%(generalData['Compare']['Oatoms'],generalData['Compare']['Tatoms'])
        DisAglCtls = generalData.get('DisAglCtls',{})
        dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData)
        if dlg.ShowModal() == wx.ID_OK:
            generalData['DisAglCtls'] = dlg.GetData()
        dlg.Destroy()
        Natm = len(data['Atoms'])
        iAtm= 0
        pgbar = wx.ProgressDialog('Process polyhedron compare for %d atoms'%Natm,'Atoms done=',Natm+1,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        Tilts = generalData['Compare']['Tilts']
        Tilts.update({bName:{'Otilts':[],'Ttilts':[]}})
        Bonds = generalData['Compare']['Bonds']
        Bonds.update({bName:{'Obonds':[],'Tbonds':[]}})
        Vects = generalData['Compare']['Vects']
        Vects.update({bName:{'Ovec':[],'Tvec':[]}})
        dVects = generalData['Compare']['dVects']
        dVects.update({bName:{'Ovec':[],'Tvec':[]}})
        Oatoms = generalData['Compare']['Oatoms']
        nOct = 0
        nTet = 0
        nElse = 0
        for iatm,atom in enumerate(data['Atoms']):
            if atom[ct] == Oatoms:
                if ran.random() > generalData['Compare']['Sampling']:
                    iAtm += 1
                    continue
                results = G2mth.FindAllNeighbors(data,atom[ct-1],atNames,Orig=iatm)[0]      #slow step
                if len(results) == 4:
                    bond,std,meanDisp,stdDisp,A,V,dVec = G2mth.FindTetrahedron(results)
                    Bonds[bName]['Tbonds'].append(bond)
                    Tilts[bName]['Ttilts'].append(A)
                    Vects[bName]['Tvec'].append(V)
                    dVects[bName]['Tvec'].append(dVec)
                    nTet += 1
                elif len(results) == 6:
                    bond,std,meanDisp,stdDisp,A,V,dVec = G2mth.FindOctahedron(results)
                    Bonds[bName]['Obonds'].append(bond)
                    Tilts[bName]['Otilts'].append(A)
                    Vects[bName]['Ovec'].append(V)
                    dVects[bName]['Ovec'].append(dVec)
                    nOct += 1
                else:
                    print('%s is something else with %d vertices'%(atom[ct-1],len(results)))
                    nElse += 1
            GoOn = pgbar.Update(iAtm,newmsg='Atoms done=%d'%(iAtm))
            iAtm += 1
            if not GoOn[0]:
                print(' Compare aborted')
                break
        print(' Compare finished; nOct =  %d, nTet = %d, nOther = %d'%(nOct,nTet,nElse))
        pgbar.Destroy()
        wx.CallAfter(UpdateGeneral,General.GetScrollPos(wx.VERTICAL))

    def OnUseBilbao(event):
        '''Select and apply a transformation matrix from the Bilbao web site
        to create a new phase
        '''
        if 'SUBGROUPS' in data:
            magData,baseList = data['SUBGROUPS']
        else:
            PatternName = data['magPhases']
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PatternName)
            UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List')
            UCdata = list(G2frame.GPXtree.GetItemPyData(UnitCellsId))
            magData = UCdata[5]
            if len(UCdata[0]) < 17:     #old version of k-SUBGROUPSMAG
                baseList = range(1,len(magData)+1)
            else:
                baseList = UCdata[0][16]
        
        magKeep = []
        magIds = []
        magchoices = []
        ifMag = False
        itemList = [phase.get('gid',ip+1) for ip,phase in enumerate(magData)]
        phaseDict = dict(zip(itemList,magData))
        for im,mid in enumerate(baseList):
            magdata = phaseDict[mid]
            if magdata['Keep']:
                if 'magAtms' in magdata:
                    ifMag = True
                magdata['No.'] = im+1
                trans = G2spc.Trans2Text(magdata['Trans'])
                vec = G2spc.Latt2text([magdata['Uvec'],])
                magKeep.append(magdata)
                magIds.append(mid)
                magchoices.append('(%d) %s; (%s) + (%s)'%(im+1,magdata['Name'],trans,vec))
        if not len(magKeep):
            G2frame.ErrorDialog('Subgroup/magnetic phase selection error','No magnetic phases found; be sure to "Keep" some')
            return
        if ifMag:
            dlg = wx.SingleChoiceDialog(G2frame,'Select magnetic space group','Make new magnetic phase',magchoices)
        else:
            dlg = wx.SingleChoiceDialog(G2frame,'Select subgroup','Make new subgroup phase',magchoices)
        opt = dlg.ShowModal()
        if opt == wx.ID_OK:
            sel = dlg.GetSelection()
            magchoice = magKeep[sel]
            magId = magIds[sel]
            if ifMag:
                phaseName = '%s-mag_%d'%(data['General']['Name'],magchoice['No.'])
            else:
                phaseName = '%s-sub_%d'%(data['General']['Name'],magchoice['No.'])
            newPhase = copy.deepcopy(data)
            newPhase['ranId'] = ran.randint(0,sys.maxsize)
            # remove subgroup info on parent from new phase
            if 'magPhases' in newPhase: del newPhase['magPhases']
            if 'SUBGROUPS' in newPhase: del newPhase['SUBGROUPS']
            generalData = newPhase['General']
            generalData['Name'] = phaseName
            generalData['SGData'] = copy.deepcopy(magchoice['SGData'])
            generalData['Cell'][1:] = magchoice['Cell'][:]
            generalData['MagDmin'] = 1.0
            SGData = generalData['SGData']
            vvec = np.array([0.,0.,0.])
            newPhase['MagXform'] = (magchoice['Trans'],magchoice['Uvec'],vvec)
            newPhase,atCodes = G2lat.TransformPhase(data,newPhase,magchoice['Trans'],magchoice['Uvec'],vvec,ifMag)
            Atoms = newPhase['Atoms']
            Atms = []
            AtCods = []
            atMxyz = []
            for ia,atom in enumerate(Atoms):
                if ifMag and atom[1] == 'O':       #skip "magnetic" O atoms
                    continue
                if ifMag and not len(G2elem.GetMFtable([atom[1],],[2.0,])):
                    continue
                atom[0] += '_%d'%ia
                atom[2] = ''                    #clear away refinement flags
                SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData)
                Atms.append(atom)
                AtCods.append(atCodes[ia])
                if ifMag:
                    MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                    CSI = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
                    atMxyz.append([MagSytSym,CSI[0]])
                else:
                    CSI = G2spc.GetCSxinel(SytSym)
                    atMxyz.append([SytSym,CSI[0]])
            if ifMag:
                dlg = UseMagAtomDialog(G2frame,magchoices[sel],Atms,AtCods,atMxyz,ifMag=ifMag,ifDelete=True)
                try:
                    opt = dlg.ShowModal()
                    if  opt == wx.ID_YES:
                        G2frame.OnFileSave(event)       #saves current state of Unit Cell List
                        newPhase['Atoms'],atCodes = dlg.GetSelection()
                        generalData['Lande g'] = len(generalData['AtomTypes'])*[2.,]
                    elif opt == wx.ID_DELETE:
                        magData[magId]['Keep'] = False
                        return
                    else:   #wx.ID_NO
                        return
                finally:
                    dlg.Destroy()
        else:
            return
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        UseList = newPhase['Histograms']
        detTrans = np.abs(nl.det(magchoice['Trans']))
        for hist in UseList:
            UseList[hist]['Scale'] /= detTrans      #scale by 1/volume ratio
            UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
            UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
        newPhase['General']['Map'] = mapDefault.copy()
        sub = G2frame.GPXtree.AppendItem(parent=
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases'),text=phaseName)
        G2frame.GPXtree.SetItemPyData(sub,newPhase)
        newPhase['Drawing'] = []
        if ifMag:
            G2cnstG.TransConstraints(G2frame,data,newPhase,magchoice['Trans'],vvec,atCodes)     #data is old phase
            G2frame.newGPXfile = phaseName.replace('.','_')+'.gpx'          #'.' in file names is a bad idea
            UCdata[5] = []      #clear away other mag choices from chem phase in new project
            G2frame.GPXtree.SetItemPyData(UnitCellsId,UCdata)
            G2frame.OnFileSaveas(event)
        G2frame.GPXtree.SelectItem(sub)

    def OnApplySubgroups(event):
        '''Select and apply a transformation matrix from the Bilbao web site
        to replace current phase with a subgroup
        '''
        if 'SUBGROUPS' in data:
            subPhases,baseList = data['SUBGROUPS']
            UnitCellsId = None
        else:
            PatternName = data['magPhases']
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PatternName)
            UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List')
            UCdata = list(G2frame.GPXtree.GetItemPyData(UnitCellsId))
            if len(UCdata[0]) < 17:     #old version of k-SUBGROUPSMAG
                baseList = range(1,len(UCdata[5])+1)
            else:
                baseList = UCdata[0][16]
            subPhases = UCdata[5]
        subKeep = []
        subIds = []
        subchoices = []
        ifMag = False
        itemList = [phase.get('gid',ip+1) for ip,phase in enumerate(subPhases)]
        phaseDict = dict(zip(itemList,subPhases))
        for im,mid in enumerate(baseList):
            if phaseDict[mid]['Keep']:
                phaseDict[mid]['No.'] = im+1
                trans = G2spc.Trans2Text(phaseDict[mid]['Trans'])
                vec = G2spc.Latt2text([phaseDict[mid]['Uvec'],])
                subKeep.append(phaseDict[mid])
                subIds.append(mid)
                subchoices.append('(%d) %s; (%s) + (%s)'%(im+1,phaseDict[mid]['Name'],trans,vec))
        if not len(subKeep):
            G2frame.ErrorDialog('Subgroup phase selection error','No subgroups available; be sure to "Keep" some')
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,
                    'Make new project using subgroups','Select subgroup(s)',
                    subchoices)
        opt = dlg.ShowModal()
        sels = []
        if opt == wx.ID_OK:
            sels = dlg.GetSelections()
        if not sels: return
        #
        G2frame.OnFileSave(None) # save
        orgFilName = G2frame.GSASprojectfile
        phsnam = data['General']['Name']
        # get restraints & clear geometrical restraints
        resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        Restraints = G2frame.GPXtree.GetItemPyData(resId)
        resId = G2gd.GetGPXtreeItemId(G2frame,resId,phsnam)
        savedRestraints = None
        if phsnam in Restraints:
            Restraints[phsnam]['Bond']['Bonds'] = []
            Restraints[phsnam]['Angle']['Angles'] = []
            savedRestraints = Restraints[phsnam]
            del Restraints[phsnam]
        orgData = copy.deepcopy(data)
        for sel in sels:
            data.update(copy.deepcopy(orgData))   # get rid of prev phase
            magchoice = subKeep[sel]
            spg = magchoice['SGData']['SpGrp'].replace(' ','')
            #subId = subIds[sel]
            # generate the new phase
            newPhase = copy.deepcopy(data)
            generalData = newPhase['General']
            generalData['SGData'] = copy.deepcopy(magchoice['SGData'])
            generalData['Cell'][1:] = magchoice['Cell'][:]
            generalData['MagDmin'] = 1.0
            SGData = generalData['SGData']
            vvec = np.array([0.,0.,0.])
            newPhase['MagXform'] = (magchoice['Trans'],magchoice['Uvec'],vvec)
            newPhase,atCodes = G2lat.TransformPhase(data,newPhase,magchoice['Trans'],magchoice['Uvec'],vvec,ifMag)
            Atoms = newPhase['Atoms']
            Atms = []
            AtCods = []
            atMxyz = []
            for ia,atom in enumerate(Atoms):
                atom[0] += '_%d'%ia
                atom[2] = ''                    #clear away refinement flags
                SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData)
                Atms.append(atom)
                AtCods.append(atCodes[ia])
                CSI = G2spc.GetCSxinel(SytSym)
                atMxyz.append([SytSym,CSI[0]])
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            UseList = newPhase['Histograms']
            detTrans = np.abs(nl.det(magchoice['Trans']))
            newPhase['Drawing'] = []
            for hist in UseList:
                UseList[hist]['Scale'] /= detTrans      #scale by 1/volume ratio
                # reset Dij & microstrain terms where # of terms changes
                if len(UseList[hist]['Mustrain'][4]) != NShkl:
                    UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
                if len(UseList[hist]['HStrain'][0]) != NDij:
                    UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
            newPhase['General']['Map'] = mapDefault.copy()
            # phase name rename
            newName = generalData['Name'] = f"{phsnam}_{magchoice['No.']}_{spg}"
            phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
            #phaseNameList = usedHistograms.keys() # phase names in use
            generalData['Name'] = newName
            G2frame.GPXtree.SetItemText(Item,generalData['Name'])
            # change phase name key in Reflection Lists for each histogram
            for hist in data['Histograms']:
                ht = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
                rt = G2gd.GetGPXtreeItemId(G2frame,ht,'Reflection Lists')
                if not rt: continue
                RfList = G2frame.GPXtree.GetItemPyData(rt)
                RfList[newName] = []
                if phsnam in RfList:
                    del RfList[phsnam]
            if savedRestraints: # restore cleared restraints
                Restraints[generalData['Name']] = savedRestraints
            if resId: G2frame.GPXtree.SetItemText(resId,newName)
            data.update(newPhase)
            #clear away prev subgroup choices
            if 'magPhases' in data: del data['magPhases']
            if 'SUBGROUPS' in data: del data['SUBGROUPS']
            if UnitCellsId:
                UCdata[5] = []
                G2frame.GPXtree.SetItemPyData(UnitCellsId,UCdata)
            # save new file
            G2frame.GSASprojectfile = os.path.splitext(orgFilName
                            )[0]+'_'+spg.replace('/','_')+'.gpx'
            #G2frame.OnFileSaveas(event)
            G2IO.ProjFileSave(G2frame)

        # restore the original saved project
        G2frame.OnFileOpen(None,filename=orgFilName,askSave=False)
        # reopen tree to the original phase
        def _ShowPhase():
            phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            G2frame.GPXtree.Expand(phId)
            phId = G2gd.GetGPXtreeItemId(G2frame,phId,phsnam)
            G2frame.GPXtree.SelectItem(phId)
        wx.CallLater(100,_ShowPhase)
####  Atom routines ################################################################################
    def FillAtomsGrid(Atoms):
        '''Display the contents of the Atoms tab
        '''
        def RefreshAtomGrid(event):

            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                for row in range(Atoms.GetNumberRows()):
                    Atoms.SelectRow(row,True)
                return
            if r < 0:                          #double click on col label! Change all atoms!
                noSkip = True
                parms = ''
                if Atoms.GetColLabelValue(c) == 'refine':
                    Type = generalData['Type']
                    if Type in ['nuclear','macromolecular','faulted',]:
                        choice = ['F - site fraction','X - coordinates','U - thermal parameters']
                    elif Type in ['magnetic',]:
                        choice = ['F - site fraction','X - coordinates','U - thermal parameters','M - magnetic moment']
                    dlg = wx.MultiChoiceDialog(G2frame,'Select','Refinement controls',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelections()
                        parms = ''
                        for x in sel:
                            parms += choice[x][0]
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'I/A':
                    choice = ['Isotropic','Anisotropic']
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Thermal Motion',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel][0]
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'Type':
                    choice = generalData['AtomTypes']
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom types',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if parms == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                    dlg.Destroy()
                    SetupGeneral()
                elif Atoms.GetColLabelValue(c) == 'Name':
                    choice = []
                    for r in range(Atoms.GetNumberRows()):
                        if str(atomData[r][c]) not in choice:
                            choice.append(str(atomData[r][c]))
                    choice.sort()
                    dlg = wx.SingleChoiceDialog(G2frame,'Select atom(s) by name',
                                            'Name Selection',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if parms == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                        dlg.Destroy()
                    else:
                        dlg.Destroy()
                        return
                elif Atoms.GetColLabelValue(c) == 'residue':
                    choice = []
                    for r in range(Atoms.GetNumberRows()):
                        if str(atomData[r][c]) not in choice:
                            choice.append(str(atomData[r][c]))
                    choice.sort()
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Residue',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if parms == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'res no':
                    choice = []
                    for r in range(Atoms.GetNumberRows()):
                        if str(atomData[r][c]) not in choice:
                            choice.append(str(atomData[r][c]))
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Residue no.',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if int(parms) == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'chain':
                    choice = []
                    for r in range(Atoms.GetNumberRows()):
                        if atomData[r][c] not in choice:
                            choice.append(atomData[r][c])
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Chain',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if parms == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'Uiso':       #this needs to ask for value
                    return                                        #& then change all 'I' atoms; now do nothing
                else:
                    return
                if noSkip:
                    ui = colLabels.index('U11')
                    us = colLabels.index('Uiso')
                    ss = colLabels.index('site sym')
                    for r in range(Atoms.GetNumberRows()):
                        ID = atomData[r][ui+6]
                        if parms != atomData[r][c] and Atoms.GetColLabelValue(c) == 'I/A':
                            if parms == 'A':                #'I' --> 'A'
                                Uiso = float(Atoms.GetCellValue(r,us))
                                sytsym = atomData[r][ss]
                                CSI = G2spc.GetCSuinel(sytsym)
                                atomData[r][ui:ui+6] = Uiso*np.array(CSI[3])
                                atomData[r][us] = 0.0
                                Atoms.SetCellStyle(r,us,VERY_LIGHT_GREY,True)
                                for i in range(6):
                                    ci = ui+i
                                    Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                                    if CSI[2][i]:
                                        Atoms.SetCellStyle(r,ci,WHITE,False)
                            else:                           #'A' --> 'I'
                                Uij = atomData[r][ui:ui+6]
                                Uiso = (Uij[0]+Uij[1]+Uij[2])/3.0
                                atomData[r][us] = Uiso
                                Atoms.SetCellStyle(r,us,WHITE,False)
                                for i in range(6):
                                    ci = ui+i
                                    atomData[r][ci] = 0.0
                                    Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                        if not Atoms.IsReadOnly(r,c):
                            if Atoms.GetColLabelValue(c) == 'refine':
                                rbExcl = rbAtmDict.get(atomData[r][ui+6],'')
                                if rbExcl:
                                    for excl in rbExcl:
                                        atomData[r][c] = parms.replace(excl,'')
                                else:
                                    atomData[r][c] = parms
                            else:
                                atomData[r][c] = parms
                        if 'Atoms' in data['Drawing']:
                            G2mth.DrawAtomsReplaceByID(data,ui+6,atomData[r],ID)
                    wx.CallAfter(Paint)

        def ChangeAtomCell(event):

            def chkUij(Uij,CSI): #needs to do something!!!
                return Uij

            SGData = generalData['SGData']
            r,c =  event.GetRow(),event.GetCol()
            replot = True
            if r >= 0 and c >= 0:
                ci = colLabels.index('I/A')
                ID = atomData[r][ci+8]
                if Atoms.GetColLabelValue(c) in ['x','y','z']:
                    ci = colLabels.index('x')
                    XYZ = atomData[r][ci:ci+3]
                    if None in XYZ:
                        XYZ = [0,0,0]
                    SScol = colLabels.index('site sym')
                    Mulcol = colLabels.index('mult')
                    Sytsym,Mult = G2spc.SytSym(XYZ,SGData)[:2]
                    atomData[r][SScol] = Sytsym
                    atomData[r][Mulcol] = Mult
                    if atomData[r][colLabels.index('I/A')] == 'A':
                        ui = colLabels.index('U11')
                        CSI = G2spc.GetCSuinel(Sytsym)
                        atomData[r][ui:ui+6] = chkUij(atomData[r][ui:ui+6],Sytsym)
                        for i in range(6):
                            ci = i+ui
                            Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                            if CSI[2][i]:
                                Atoms.SetCellStyle(r,ci,WHITE,False)
                    SetupGeneral()
                elif Atoms.GetColLabelValue(c) == 'Type':
                    AtomTypeSelect(event)
                elif Atoms.GetColLabelValue(c) == 'I/A':            #note use of text color to make it vanish!
                    if atomData[r][c] == 'I':
                        Uij = atomData[r][c+2:c+8]
                        atomData[r][c+1] = (Uij[0]+Uij[1]+Uij[2])/3.0
                        Atoms.SetCellStyle(r,c+1,WHITE,False)
                        Atoms.SetCellTextColour(r,c+1,BLACK)
                        for i in range(6):
                            ci = i+colLabels.index('U11')
                            Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                            Atoms.SetCellTextColour(r,ci,VERY_LIGHT_GREY)
                            atomData[r][ci] = 0.0
                    else:
                        value = atomData[r][c+1]
                        CSI = G2spc.GetCSuinel(atomData[r][colLabels.index('site sym')])
                        Atoms.SetCellStyle(r,c+1,VERY_LIGHT_GREY,True)
                        Atoms.SetCellTextColour(r,c+1,VERY_LIGHT_GREY)
                        for i in range(6):
                            ci = i+colLabels.index('U11')
                            atomData[r][ci] = value*CSI[3][i]
                            Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                            Atoms.SetCellTextColour(r,ci,BLACK)
                            if CSI[2][i]:
                                Atoms.SetCellStyle(r,ci,WHITE,False)
                        atomData[r][c+1] =  0.0
                elif Atoms.GetColLabelValue(c) in ['U11','U22','U33','U12','U13','U23']:
                    value = atomData[r][c]
                    CSI = G2spc.GetCSuinel(atomData[r][colLabels.index('site sym')])
                    iUij = CSI[0][c-colLabels.index('U11')]
                    for i in range(6):
                        if iUij == CSI[0][i]:
                            atomData[r][i+colLabels.index('U11')] = value*CSI[1][i]
                elif Atoms.GetColLabelValue(c) == 'refine':
                    ci = colLabels.index('I/A')
                    atomData[r][c] = atomData[r][c].replace(rbAtmDict.get(atomData[r][ci+8],''),'')
                    replot = False
                elif Atoms.GetColLabelValue(c) in ['Mx','My','Mz']:
                    value = atomData[r][c]
                    cx = colLabels.index('x')
                    SpnFlp = generalData['SGData']['SpnFlp']
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(atomData[r][cx:cx+3],SGData)
                    CSI = G2spc.GetCSpqinel(SpnFlp,dupDir)
                    iM = CSI[0][c-colLabels.index('Mx')]
                    for i in range(3):
                        if iM == CSI[0][i]:
                            atomData[r][i+colLabels.index('Mx')] = value*CSI[1][i]
                if 'Atoms' in data['Drawing'] and replot:
                    ci = colLabels.index('I/A')
                    G2mth.DrawAtomsReplaceByID(data,ci+8,atomData[r],ID)
                    G2frame.GetStatusBar().SetStatusText('Structure changed: Do "Edit Atoms/Update draw atoms" to refresh structure drawing',1)
                    # G2plt.PlotStructure(G2frame,data)
                if SGData['SpGrp'] != 'P 1':    #no need to update P 1 structures!
                    wx.CallAfter(Paint)

        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if Atoms.GetColLabelValue(c) == 'Type':
                PE = G2elemGUI.PickElement(G2frame,ifMag=ifMag)
                if PE.ShowModal() == wx.ID_OK:
                    if PE.Elem != 'None':
                        atomData[r][c] = PE.Elem.strip()
                        name = atomData[r][c]
                        if len(name) in [2,4]:
                            atomData[r][c-1] = name[:2]+'%d'%(r+1)
                        else:
                            atomData[r][c-1] = name[:1]+'%d'%(r+1)
                PE.Destroy()
                SetupGeneral()
                wx.CallAfter(Paint)
                value = Atoms.GetCellValue(r,c)
                atomData[r][c] = value
                ci = colLabels.index('I/A')
                ID = atomData[r][ci+8]
                if 'Atoms' in data['Drawing']:
                    G2mth.DrawAtomsReplaceByID(data,ci+8,atomData[r],ID)
                    G2frame.GetStatusBar().SetStatusText('Structure changed: Do "Edit Atoms/Update draw atoms" to refresh structure drawing',1)
                    # G2plt.PlotStructure(G2frame,data)
                SetupGeneral()
            else:
                event.Skip()

        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if not (event.AltDown() or (event.ShiftDown() and event.ControlDown())):
                Atoms.frm = -1
            if r < 0 and c < 0:
                if Atoms.IsSelection():
                    Atoms.ClearSelection()
            elif c < 0:                   #only row clicks
                ci = colLabels.index('I/A')
                if event.ControlDown() and not event.ShiftDown():
                    if r in Atoms.GetSelectedRows():
                        Atoms.DeselectRow(r)
                    else:
                        Atoms.SelectRow(r,True)
                elif event.ShiftDown() and not event.ControlDown():
                    indxs = Atoms.GetSelectedRows()
                    Atoms.ClearSelection()
                    ibeg = 0
                    if indxs:
                        ibeg = indxs[-1]
                    for row in range(ibeg,r+1):
                        Atoms.SelectRow(row,True)
                elif event.AltDown() or (event.ShiftDown() and event.ControlDown()):
                    if atomData[r][ci+8] in rbAtmDict:
                        G2frame.ErrorDialog('Atom move error','Atoms in rigid bodies can not be moved')
                        Atoms.frm = -1
                        Atoms.ClearSelection()
                    else:
                        if Atoms.frm < 0:           #pick atom to be moved
                            Atoms.frm = r
                            Atoms.SelectRow(r,True)
                            n = colLabels.index('Name')
                            G2frame.GetStatusBar().SetStatusText('Atom '+atomData[r][n]+' is to be moved',1)
                        else:                       #move it
                            item = atomData.pop(Atoms.frm)
                            atomData.insert(r,item)
                            Atoms.frm = -1
                            G2frame.GetStatusBar().SetStatusText('',1)
                            data['Drawing']['Atoms'] = []           #clear & rebuild Draw atoms table
                            UpdateDrawAtoms(G2frame,data)
                            wx.CallAfter(Paint)
                else:
                    G2frame.GetStatusBar().SetStatusText('Use right mouse click to brng up Atom editing options',1)
                    Atoms.ClearSelection()
                    Atoms.SelectRow(r,True)
            G2frame.GetStatusBar().SetStatusText('Structure changed: Do "Edit Atoms/Update draw atoms" to refresh structure drawing',1)
            # G2plt.PlotStructure(G2frame,data)

        def ChangeSelection(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                Atoms.ClearSelection()
            if c < 0:
                if r in Atoms.GetSelectedRows():
                    Atoms.DeselectRow(r)
                else:
                    Atoms.SelectRow(r,True)
            if r < 0:
                if c in Atoms.GetSelectedCols():
                    Atoms.DeselectCol(c)
                else:
                    Atoms.SelectCol(c,True)

        def Paint(Scroll=0):
            'Place atom info into the table'
            table = []
            rowLabels = []
            for i,atom in enumerate(atomData):
                table.append(atom)
                rowLabels.append(str(i))
            atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            try:
                Atoms.SetTable(atomTable, True, useFracEdit=False)    # Paint may be called after the Grid has been deleted
            except:
                return
            Atoms.frm = -1
            colType = colLabels.index('Type')
            colR = colLabels.index('refine')
            colSS = colLabels.index('site sym')
            colF = colLabels.index('frac')
            colX = colLabels.index('x')
            colIA = colLabels.index('I/A')
            colU11 = colLabels.index('U11')
            colUiso = colLabels.index('Uiso')
            colM = 0
            if 'Mx' in colLabels:
                colM = colLabels.index('Mx')
            # loop over all cols in table, set cell editor for numerical items
            for c,t in enumerate(Types):
                if not t.startswith(wg.GRID_VALUE_FLOAT): continue
                attr = wx.grid.GridCellAttr()
                attr.IncRef()               #fix from Jim Hester
                attr.SetEditor(G2G.GridFractionEditor(Atoms))
                if c in range(colU11-1,colU11+6):
                    Atoms.SetColSize(c,50)
                    attr.SetBackgroundColour(VERY_LIGHT_GREY) # make invisible
                    attr.SetTextColour(VERY_LIGHT_GREY)
                    attr.SetReadOnly(True)
                Atoms.SetColAttr(c, attr)
            for row in range(Atoms.GetNumberRows()):    #this is slow for large numbers of atoms
                atId = atomData[row][colIA+8]
                rbExcl = rbAtmDict.get(atId,'')
                if Atoms.GetCellValue(row,colIA) == 'A':
                    try:    #patch for sytsym name changes
                        CSI = G2spc.GetCSuinel(atomData[row][colSS])
                    except KeyError:
                        Sytsym,Mult = G2spc.SytSym(atomData[row][colX:colX+3],SGData)[:2]
                        atomData[row][colSS] = Sytsym
                        atomData[row][colSS+1] = Mult
                        CSI = G2spc.GetCSuinel(Sytsym)
                    Atoms.SetCellStyle(row,colUiso,VERY_LIGHT_GREY,True)
                    Atoms.SetCellTextColour(row,colUiso,VERY_LIGHT_GREY)
                    for i in range(6):
                        cj = colU11+i
                        Atoms.SetCellTextColour(row,cj,BLACK)
                        if G2lat.Uij2Ueqv(atomData[row][colU11:colU11+6],GS,Amat)[1]:
                            Atoms.SetCellTextColour(row,cj,RED)
                        Atoms.SetCellStyle(row,cj,VERY_LIGHT_GREY,True)
                        if CSI[2][i] and 'U' not in rbExcl:
                            Atoms.SetCellStyle(row,cj,WHITE,False)
                else:
                    Atoms.SetCellStyle(row,colUiso,WHITE,False)
                    Atoms.SetCellTextColour(row,colUiso,BLACK)
                    if atomData[row][colUiso] < 0.:
                        Atoms.SetCellTextColour(row,colUiso,RED)
                    if 'U' in rbExcl:
                        Atoms.SetCellStyle(row,colUiso,VERY_LIGHT_GREY,True)
                if colM:
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(atomData[row][colX:colX+3],SGData)
                    MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                    Atoms.SetCellValue(row,colSS,MagSytSym)
                    CSI = []
                    if not SGData['SGGray']:
                        CSI = G2spc.GetCSpqinel(SpnFlp,dupDir)
                    saveCSI = 0
                    for i in range(3):
                        ci = i+colM
                        Atoms.SetCellStyle(row,ci,VERY_LIGHT_GREY,True)
                        Atoms.SetCellTextColour(row,ci,VERY_LIGHT_GREY)
                        if CSI and CSI[1][i]:       # and AtInfo and AtInfo[atomData[row][colType]]:
                            if saveCSI != CSI[0][i]:
                                Atoms.SetCellStyle(row,ci,WHITE,False)
                                saveCSI = CSI[0][i]
                            Atoms.SetCellTextColour(row,ci,BLACK)
                # if 'F' in rbExcl:
                #     Atoms.SetCellStyle(row,colF,VERY_LIGHT_GREY,True)
                if 'X' in rbExcl:
                    for c in range(0,colX+3):
                        if c != colR:
                            Atoms.SetCellStyle(row,c,VERY_LIGHT_GREY,True)
                Atoms.SetReadOnly(row,colType,True)
                Atoms.SetReadOnly(row,colSS,True)                         #site sym
                Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            oldSizer = AtomList.GetSizer()
            if oldSizer:  # 2nd+ use, clear out old entries
                for i in oldSizer.GetChildren(): # look for grids in sizer
                    if type(i.GetWindow()) is G2G.GSGrid:
                        oldSizer.Detach(i.GetWindow())  # don't delete them
                oldSizer.Clear(True)
            Atoms.AutoSizeColumns(False)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add(Atoms,1,wx.EXPAND)
            Atoms.SetScrollRate(10,10) # allow grid to scroll
            SetPhaseWindow(AtomList,mainSizer,Scroll=Atoms.GetScrollPos(wx.VERTICAL))

#### FillAtomsGrid main code
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        generalData = data['General']
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Atom parameters list for {generalData['Name']}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)
        cell = generalData['Cell'][1:7]
        Amat = G2lat.cell2AB(cell)[0]
        GS = G2lat.cell2GS(cell)
        SpnFlp = generalData['SGData'].get('SpnFlp',[])
        atomData = data['Atoms']
        resRBData = data['RBModels'].get('Residue',[])
        vecRBData = data['RBModels'].get('Vector',[])
        global rbAtmDict
        rbAtmDict = {}
        for rbObj in resRBData+vecRBData:
            exclList = ['FX' for i in range(len(rbObj['Ids']))]
            rbAtmDict.update(dict(zip(rbObj['Ids'],exclList)))
            if rbObj['ThermalMotion'][0] != 'None':
                for id in rbObj['Ids']:
                    rbAtmDict[id] += 'U'
        # exclList will be 'fx' or 'fxu' if TLS used in RB
        Items = [G2G.wxID_ATOMSEDITINSERT,G2G.wxID_ATOMSEDITDELETE,
            G2G.wxID_ATOMSMODIFY,G2G.wxID_ATOMSTRANSFORM,G2G.wxID_MAKEMOLECULE,
            G2G.wxID_ATOMVIEWINSERT,G2G.wxID_ATOMMOVE,G2G.wxID_ADDHATOM]
        if atomData:
            for item in Items:
                G2frame.dataWindow.AtomsMenu.Enable(item,True)
        else:
            for item in Items:
                G2frame.dataWindow.AtomsMenu.Enable(item,False)
        Items = [G2G.wxID_ATOMVIEWINSERT, G2G.wxID_ATOMSVIEWADD,G2G.wxID_ATOMMOVE]
        if 'showABC' in data['Drawing']:
            for item in Items:
                G2frame.dataWindow.AtomsMenu.Enable(item,True)
        else:
            for item in Items:
                G2frame.dataWindow.AtomsMenu.Enable(item,False)
        parmChoice = ': ,X,XU,U,F,FX,FXU,FU'
        if generalData['Type'] == 'magnetic':
            parmChoice += ',M,MX,MXU,MU,MF,MFX,MFXU,MFU'
        AAchoice = ": ,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE,HOH,UNK"
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+parmChoice,]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+[wg.GRID_VALUE_FLOAT+':10,4', #x,y,z,frac
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+":I,A",]
        Types += 7*[wg.GRID_VALUE_FLOAT+':10,5',]
        colLabels = ['Name','Type','refine','x','y','z','frac','site sym','mult','I/A','Uiso','U11','U22','U33','U12','U13','U23']
        ifMag = False
        if generalData['Type'] == 'macromolecular':
            colLabels = ['res no','residue','chain'] + colLabels
            Types = [wg.GRID_VALUE_STRING,
                wg.GRID_VALUE_CHOICE+AAchoice,
                wg.GRID_VALUE_STRING] + Types
        elif generalData['Type'] == 'magnetic':
            ifMag = True
            colLabels = colLabels[:7]+['Mx','My','Mz']+colLabels[7:]
            Types = Types[:7]+3*[wg.GRID_VALUE_FLOAT+':10,4',]+Types[7:]
        SGData = data['General']['SGData']
        G2frame.GetStatusBar().SetStatusText('',1)
        if SGData['SpGrp'] in G2spc.spg2origins:
            G2frame.GetStatusBar().SetStatusText('Warning: Atom positions must correspond to 2nd setting for the space group '+SGData['SpGrp'],1)
        if SGData['SGPolax']:
            G2frame.GetStatusBar().SetStatusText('Warning: The location of the origin is arbitrary in '+SGData['SGPolax'],1)
        if 'phoenix' in wx.version():
            Atoms.Unbind(wg.EVT_GRID_CELL_CHANGED)
            Atoms.Bind(wg.EVT_GRID_CELL_CHANGED, ChangeAtomCell)
        else:
            Atoms.Unbind(wg.EVT_GRID_CELL_CHANGE)
            Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeAtomCell)
        Atoms.Unbind(wg.EVT_GRID_CELL_LEFT_DCLICK)
        Atoms.Unbind(wg.EVT_GRID_LABEL_LEFT_DCLICK)
        Atoms.Unbind(wg.EVT_GRID_LABEL_LEFT_CLICK)
        Atoms.Unbind(wg.EVT_GRID_LABEL_RIGHT_CLICK)
        Atoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, ChangeSelection)
        Atoms.SetRowLabelSize(45)
        lblList = ('Set refine flags','Modify selected parameters',
                       'Insert new before selected','Transform selected',
                       'Set selected xyz to view point','Select all',
                       'Select from list','Set view point from selected',
                       'Delete')
        callList = (AtomRefine,AtomModify,OnAtomInsert,AtomTransform,
                        OnAtomMove,OnSetAll,OnSetbyList,SetAtomsViewPoint,
                        AtomDelete)
        onRightClick = Atoms.setupPopup(lblList,callList)
        Atoms.Bind(wg.EVT_GRID_CELL_RIGHT_CLICK, onRightClick)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, onRightClick)
        Atoms.SetMargins(0,0)

        Paint()

    def OnAtomAdd(event):
        Elem = 'H'
        if data['General']['Type'] == 'magnetic':
            Elem = 'Fe'
        AtomAdd(0.,0.,0.,El=Elem)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        if data['Drawing']:
            G2plt.PlotStructure(G2frame,data)

    def OnAtomViewAdd(event):
        Elem = 'H'
        if data['General']['Type'] == 'magnetic':
            Elem = 'Fe'
        try:
            drawData = data['Drawing']
            x,y,z = drawData['viewPoint'][0]
            AtomAdd(x,y,z,El=Elem)
        except:
            AtomAdd(0.,0.,0.,El=Elem)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def AtomAdd(x,y,z,El='H',Name='UNK',update=True):
        atomData = data['Atoms']
        generalData = data['General']
        atId = ran.randint(0,sys.maxsize)
        if 'Q' in El:   #dummy fill spin rb pointer
            generalData['SpnIds'][atId] = -1
        SGData = generalData['SGData']
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)[:2]
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,Name,'',Name,El,'',x,y,z,1.,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] in ['nuclear','faulted',]:
            if generalData['Modulated']:
                atomData.append([Name,El,'',x,y,z,1.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.append([Name,El,'',x,y,z,1.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'magnetic':
            if generalData['Modulated']:
                atomData.append([Name,El,'',x,y,z,1.,0.,0.,0.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.append([Name,El,'',x,y,z,1.,0.,0.,0.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])

        SetupGeneral()
        # might be better to add new atom to Draw Atoms
        data['Drawing']['Atoms'] = []
#        if 'Atoms' in data['Drawing']:
#            DrawAtomAdd(data['Drawing'],atomData[-1])
        if update:
            UpdateDrawAtoms(G2frame,data)
            G2plt.PlotStructure(G2frame,data)

    def OnAtomInsert(event):
        '''Inserts a new atom into list immediately before every selected atom
        '''
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        for a in reversed(sorted(indx)):
            AtomInsert(a,0.,0.,0.)
        event.StopPropagation()
        FillAtomsGrid(Atoms)
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def OnAtomViewInsert(event):
        if 'Drawing' in data:
            drawData = data['Drawing']
            x,y,z = drawData['viewPoint'][0]
            AtomAdd(x,y,z)
            FillAtomsGrid(Atoms)
        event.StopPropagation()

    def OnHydAtomAdd(event):
        '''Adds H atoms to fill out coordination sphere for selected atoms
        '''
        cx,ct,cs,cia = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx: return
        DisAglCtls = {}
        generalData = data['General']
        if 'DisAglCtls' in generalData:
            DisAglCtls = generalData['DisAglCtls']
        dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData,Reset=False)
        if dlg.ShowModal() == wx.ID_OK:
            DisAglCtls = dlg.GetData()
            if 'H' not in DisAglCtls['AtomTypes']:
                DisAglCtls['AtomTypes'].append('H')
                DisAglCtls['AngleRadii'].append(0.5)
                DisAglCtls['BondRadii'].append(0.5)
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        generalData['DisAglCtls'] = DisAglCtls
        atomData = data['Atoms']
        AtNames = [atom[ct-1] for atom in atomData]
        AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
        Neigh = []
        AddHydIds = []
        for ind in indx:
            atom = atomData[ind]
            if atom[ct] not in ['C','N','O']:
                continue
            neigh = [atom[ct-1],G2mth.FindNeighbors(data,atom[ct-1],AtNames),0]
            if len(neigh[1][0]) > 3 or (atom[ct] == 'O' and len(neigh[1][0]) > 1):
                continue
            nH = 1      #for O atom
            if atom[ct] in ['C','N']:
                nH = 4-len(neigh[1][0])
            bonds = {item[0]:item[1:] for item in neigh[1][0]}
            nextName = ''
            if len(bonds) == 1:
                nextName = list(bonds.keys())[0]
            for bond in bonds:
                if 'C' in atom[ct]:
                    if 'C' in bond and bonds[bond][0] < 1.42:
                        nH -= 1
                        break
                    elif 'O' in bond and bonds[bond][0] < 1.3:
                        nH -= 1
                        break
                elif 'O' in atom[ct] and 'C' in bonds and bonds[bond][0] < 1.3:
                    nH -= 1
                    break
            nextneigh = []
            if nextName:
                nextneigh = G2mth.FindNeighbors(data,nextName,AtNames,notName=neigh[0])
                if nextneigh[0]:
                    neigh[1][1].append(nextneigh[1][1][0])
            neigh[2] = max(0,nH)  #set expected no. H's needed
            if len(neigh[1][0]):
                AddHydIds.append(neigh[1][1])
                Neigh.append(neigh)
        if Neigh:
            letters = ['A','B','C']
            HydIds = {}
            mapError = False
            dlg = AddHatomDialog(G2frame,Neigh,data)
            if dlg.ShowModal() == wx.ID_OK:
                Nat = len(atomData)
                Neigh = dlg.GetData()
                mapData = generalData['Map']
                for ineigh,neigh in enumerate(Neigh):
                    AddHydIds[ineigh].append(neigh[2])
                    loc = AtLookUp[AddHydIds[ineigh][0]]+1
                    if 'O' in neigh[0] and (not len(mapData['rho']) or not 'delt-F' in mapData['MapType']):
                        mapError = True
                        continue
                    Hxyz,HU = G2mth.AddHydrogens(AtLookUp,generalData,atomData,AddHydIds[ineigh])
                    for iX,X in enumerate(Hxyz):
                        AtomInsert(loc+iX,X[0],X[1],X[2],'H','H%s'%(neigh[0][1:]+letters[iX]))
                        data['Atoms'][loc+iX][cia+1] = HU[iX]
                        Id = data['Atoms'][loc+iX][cia+8]
                        HydIds[Id] = [iX,AddHydIds[ineigh]]
                        Nat += 1
                        AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
            if mapError:
                G2frame.ErrorDialog('Add H atom error','Adding O-H atoms requires delt-F map')
            SetupGeneral()
            data['General']['HydIds'].update(HydIds)
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,True)
            data['Drawing']['Atoms'] = []
            UpdateDrawAtoms(G2frame,data)
            FillAtomsGrid(Atoms)
            dlg.Destroy()
            G2plt.PlotStructure(G2frame,data)
        else:
            wx.MessageBox('No candidates found',caption='Add H atom Error',style=wx.ICON_EXCLAMATION)

    def OnHydAtomUpdate(event):
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        atomData = data['Atoms']
        AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
        HydIds = data['General']['HydIds']
        delList = []
        for HId in HydIds:
            hydIds = HydIds[HId]
            num = hydIds[0]
            Hxyz,HU = G2mth.AddHydrogens(AtLookUp,generalData,atomData,hydIds[1])
            try:
                if data['Atoms'][AtLookUp[HId]][ct] != 'H':
                    raise KeyError
                data['Atoms'][AtLookUp[HId]][cx:cx+3] = Hxyz[num]
                data['Atoms'][AtLookUp[HId]][cia+1] = HU[num]
            except KeyError:
                delList.append(HId)
                continue
        for HId in delList: #clear out deleted H-atom pointers
            del HydIds[HId]
        if not len(HydIds):
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,False)
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        FillAtomsGrid(Atoms)
        G2plt.PlotStructure(G2frame,data)

    def OnAtomMove(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        drawData = data['Drawing']
        atomData = data['Atoms']
        x,y,z = drawData['viewPoint'][0]
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        cx = colLabels.index('x')
        ci = colLabels.index('I/A')
        if len(indx) != 1:
            G2frame.ErrorDialog('Atom move error','Only one atom can be moved')
        elif atomData[indx[0]][ci+8] in rbAtmDict:
            G2frame.ErrorDialog('Atom move error','Atoms in rigid bodies can not be moved')
        else:
            atomData[indx[0]][cx:cx+3] = [x,y,z]
            SetupGeneral()
            FillAtomsGrid(Atoms)
            ID = atomData[indx[0]][ci+8]
            G2mth.DrawAtomsReplaceByID(data,ci+8,atomData[indx[0]],ID)
            G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()

    def AtomInsert(indx,x,y,z,El='H',Name='UNK'):
        atomData = data['Atoms']
        generalData = data['General']
        SGData = generalData['SGData']
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)[:2]
        atId = ran.randint(0,sys.maxsize)
        if generalData['Type'] == 'macromolecular':
            atomData.insert(indx,[0,Name,'',Name,El,'',x,y,z,1.,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] in ['nuclear','faulted',]:
            if generalData['Modulated']:
                atomData.insert(indx,[Name,El,'',x,y,z,1,Sytsym,Mult,0.,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.insert(indx,[Name,El,'',x,y,z,1.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
            SetupGeneral()
        elif generalData['Type'] == 'magnetic':
            if generalData['Modulated']:
                atomData.insert(indx,[Name,El,'',x,y,z,1.,0.,0.,0.,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.insert(indx,[Name,El,'',x,y,z,1.,0.,0.,0.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def AtomDelete(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms)
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        delList = ''
        for i in indx:
            if delList: delList += ', '
            delList += data['Atoms'][i][0]
        dlg = wx.MessageDialog(G2frame,'Do you want to delete atom(s): {}?'.format(delList),
            'Confirm delete',wx.YES|wx.NO)
        try:
            dlg.CenterOnParent()
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result != wx.ID_YES: return
        HydIds = data['General']['HydIds']
        ci = colLabels.index('I/A')
        IDs = []
        if not indx: return
        atomData = data['Atoms']
        indx.sort()
        indx.reverse()
        msg = ''
        for ind in indx:
            atom = atomData[ind]
            if atom[ci+8] in rbAtmDict:
                if msg: msg += ', '
                msg += atom[0]
            else:
                if atom[ci+8] in HydIds:    #remove Hs from Hatom update dict
                    del HydIds[atom[ci+8]]
                IDs.append(atom[ci+8])
                del atomData[ind]
        if msg: G2frame.ErrorDialog('Atom delete error',
                    'ERROR - atom(s) in a rigid body were not deleted: '+msg)
        if 'Atoms' in data['Drawing']:
            Atoms.ClearSelection()
            DrawAtomsDeleteByIDs(IDs)
            data['Drawing']['Atoms'] = []
            UpdateDrawAtoms(G2frame,data)
            wx.CallAfter(FillAtomsGrid,Atoms)
            G2plt.PlotStructure(G2frame,data)
        SetupGeneral()
        if not len(HydIds):
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,False)
        event.StopPropagation()

    def AtomRefine(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx: return
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        c = colLabels.index('refine')
        atomData = data['Atoms']
        generalData = data['General']
        Type = generalData['Type']
        if Type in ['nuclear','macromolecular','faulted',]:
            choice = ['F - site fraction','X - coordinates','U - thermal parameters']
        elif Type in ['magnetic',]:
            choice = ['F - site fraction','X - coordinates','U - thermal parameters','M - magnetic moment']
        dlg = wx.MultiChoiceDialog(G2frame,'Select','Refinement controls',choice)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            parms = ''
            for x in sel:
                parms += choice[x][0]
            for r in indx:
                if not Atoms.IsReadOnly(r,c):
                    atomData[r][c] = parms
            Atoms.ForceRefresh()
        dlg.Destroy()

    def AtomModify(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx: return
        atomData = data['Atoms']
        generalData = data['General']
        ifMag = False
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        ci = colLabels.index('I/A')
        choices = ['Type','Name','x','y','z','frac','I/A','Uiso','Uij']
        if generalData['Type'] == 'magnetic':
            choices += ['Mx','My','Mz',]
            ifMag = True
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom parameter',choices)
        parm = ''
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            parm = choices[sel]
            cid = -1
            if parm != 'Uij':
                cid = colLabels.index(parm)
        dlg.Destroy()
        if parm in ['Type']:
            dlg = G2elemGUI.PickElement(G2frame,ifMag=ifMag)
            if dlg.ShowModal() == wx.ID_OK:
                if dlg.Elem not in ['None']:
                    El = dlg.Elem.strip()
                    for r in indx:
                        if not Atoms.IsReadOnly(r,0):   #not if in RB!
                            atomData[r][cid] = El
                            if len(El) in [2,4]:
                                atomData[r][cid-1] = El[:2]+'%d'%(r+1)
                            else:
                                atomData[r][cid-1] = El[:1]+'%d'%(r+1)
                    SetupGeneral()
                    if 'Atoms' in data['Drawing']:
                        for r in indx:
                            ID = atomData[r][ci+8]
                            G2mth.DrawAtomsReplaceByID(data,ci+8,atomData[r],ID)
                FillAtomsGrid(Atoms)
            dlg.Destroy()
        elif parm in ['Name',]:
            dlg = wx.MessageDialog(G2frame,'Do you really want to rename the selected atoms?','Rename',
                wx.YES_NO | wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    for r in indx:
                        if not Atoms.IsReadOnly(r,0):   #not if in RB!
                            El = atomData[r][cid+1]
                            if len(El) in [2,4]:
                                atomData[r][cid] = El[:2]+'%d'%(r+1)
                            else:
                                atomData[r][cid] = El[:1]+'%d'%(r+1)
                FillAtomsGrid(Atoms)
            finally:
                dlg.Destroy()

        elif parm in ['I/A']:
            choices = ['Isotropic','Anisotropic']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Thermal parameter model',choices)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parm = choices[sel][0]
                for r in indx:
                    if Atoms.IsReadOnly(r,0):   #not if in RB!
                        continue
                    atomData[r][cid] = parm
                    if parm == 'A' and not any(atomData[r][cid+2:cid+8]):
                        sytsym = atomData[r][cs]
                        CSI = G2spc.GetCSuinel(sytsym)
                        atomData[r][ci+2:ci+8] = atomData[r][ci+1]*np.array(CSI[3])
                FillAtomsGrid(Atoms)
            dlg.Destroy()
        elif parm in ['frac','Uiso']:
            limits = [-1.,1.]
            val = 1.0
            if  parm in ['Uiso']:
                limits = [-0.25,0.25]
                val = 0.01
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for '+parm,val,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in indx:
#                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][cid] = parm
                SetupGeneral()
                FillAtomsGrid(Atoms)
            dlg.Destroy()
        elif parm == 'Uij':
            limits = [0.,0.25]
            val = 0.01
            dlg = G2G.SingleFloatDialog(G2frame,'Convert Uiso to Uij','Enter default for Uiso'+parm,val,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in indx:
                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][ci] = 'A'
                        sytsym = atomData[r][cs]
                        CSI = G2spc.GetCSuinel(sytsym)
                        if atomData[r][ci+1] > 0.:
                            atomData[r][ci+2:ci+8] = atomData[r][ci+1]*np.array(CSI[3])
                        else:
                            atomData[r][ci+2:ci+8] = parm*np.array(CSI[3])
                SetupGeneral()
                FillAtomsGrid(Atoms)
            dlg.Destroy()

        elif parm in ['x','y','z']:
            limits = [-1.,1.]
            val = 0.
            dlg = G2G.SingleFloatDialog(G2frame,'Atom shift','Enter shift for '+parm,val,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in indx:
                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][cid] += parm
                SetupGeneral()
                FillAtomsGrid(Atoms)
            dlg.Destroy()
        elif parm in ['Mx','My','Mz',]:
            limits = [-10.,10.]
            val = 0.
            dlg = G2G.SingleFloatDialog(G2frame,'Atom moment','Enter new value for '+parm,val,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in indx:
                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][cid] = parm
                SetupGeneral()
                FillAtomsGrid(Atoms)
            dlg.Destroy()

        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def AtomTransform(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx: return
        generalData = data['General']
        SpnFlp = generalData['SGData'].get('SpnFlp',[])
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        cx = colLabels.index('x')
        cuia = colLabels.index('I/A')
        cuij = colLabels.index('U11')
        css = colLabels.index('site sym')
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        atomData = data['Atoms']
        SGData = generalData['SGData']
        dlg = SymOpDialog(G2frame,SGData,True,True)
        New = False
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Inv,Cent,Opr,Cell,New,Force = dlg.GetSelection()
                Cell = np.array(Cell)
                cent = SGData['SGCen'][Cent]
                M,T = SGData['SGOps'][Opr]
                for ind in indx:
                    XYZ = np.array(atomData[ind][cx:cx+3])
                    XYZ = np.inner(M,XYZ)+T
                    if Inv and not SGData['SGFixed']:
                        XYZ = -XYZ
                    XYZ = XYZ+cent+Cell
                    if Force:
                        XYZ,cell = G2spc.MoveToUnitCell(XYZ)
                        Cell += cell
                    if New:
                        atom = copy.copy(atomData[ind])
                    else:
                        atom = atomData[ind]
                    atom[cx:cx+3] = XYZ
                    atom[css:css+2] = G2spc.SytSym(XYZ,SGData)[:2]
                    OprNum = ((Opr+1)+100*Cent)*(1-2*Inv)
                    if atom[cuia] == 'A':
                        Uij = atom[cuij:cuij+6]
                        U = G2spc.Uij2U(Uij)
                        U = np.inner(np.inner(M,U),M)
                        Uij = G2spc.U2Uij(U)
                        atom[cuij:cuij+6] = Uij
                    if cmx:
                        opNum = G2spc.GetOpNum(OprNum,SGData)
                        mom = np.array(atom[cmx:cmx+3])
                        atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SpnFlp[opNum-1]
                    if New:
                        atomData.append(atom)
        finally:
            dlg.Destroy()
        Atoms.ClearSelection()
        if New:
            FillAtomsGrid(Atoms)
        else:
            Atoms.ForceRefresh()
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def CollectAtoms(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        Ind = getAtomSelections(Atoms,ct-1)
        if Ind:
            choice = ['x=0','y=0','z=0','origin','center']
            dlg = wx.SingleChoiceDialog(G2frame,'Atoms closest to:','Select',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()+1
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
            wx.BeginBusyCursor()
            data['Atoms'] = G2mth.AtomsCollect(data,Ind,sel)
            wx.EndBusyCursor()
            Atoms.ClearSelection()
            data['Drawing']['Atoms'] = []
            OnReloadDrawAtoms(event)
            FillAtomsGrid(Atoms)

    def MakeMolecule(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        DisAglCtls = {}
        if indx is not None and len(indx) == 1:
            generalData = data['General']
            if 'DisAglCtls' in generalData:
                DisAglCtls = generalData['DisAglCtls']
            dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData)
            if dlg.ShowModal() == wx.ID_OK:
                DisAglCtls = dlg.GetData()
            else:
                dlg.Destroy()
                return
            dlg.Destroy()
            generalData['DisAglCtls'] = DisAglCtls
            atomData = copy.deepcopy(data['Atoms'])
            result = G2mth.FindMolecule(indx[0],generalData,atomData)
            if 'str' in str(type(result)):
                G2frame.ErrorDialog('Assemble molecule',result)
            else:
                data['Atoms'] = result
            Atoms.ClearSelection()
            data['Drawing']['Atoms'] = []
            OnReloadDrawAtoms(event)
            FillAtomsGrid(Atoms)
#                G2frame.ErrorDialog('Distance/Angle calculation','try again but do "Reset" to fill in missing atom types')
        else:
            print ("select one atom")
            G2frame.ErrorDialog('Select one atom',"select one atom to begin molecule build then redo")

    def OnDensity(event):
        'show the density for the current phase'
        density,mattCoeff = G2mth.getDensity(data['General'])
        msg = f"Density of phase {data['General']['Name']!r} = {density:.3f} g/cc"
        msg += '\n\nCell contents:'
        for typ,num in data['General']['NoAtoms'].items():
            msg += f'\n      {typ}: {num}'
        print(msg)
        G2G.G2MessageBox(G2frame,msg,'Density')

    def OnSetAll(event):
        'set all atoms in table as selected'
        for row in range(Atoms.GetNumberRows()):
            Atoms.SelectRow(row,True)
        G2plt.PlotStructure(G2frame,data)

    def OnSetbyList(event):
        'select atoms using a filtered listbox'
        choices = [atm[0] for atm in data['Atoms']]
        dlg = G2G.G2MultiChoiceDialog(G2frame,
                    'Select atoms','Choose atoms to select',choices)
        indx = []
        if dlg.ShowModal() == wx.ID_OK:
            indx = dlg.GetSelections()
        dlg.Destroy()
        if len(indx) == 0: return
        Atoms.ClearSelection()
        for row in indx:
            Atoms.SelectRow(row,True)
        G2plt.PlotStructure(G2frame,data)

    def SetAtomsViewPoint(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx: return
        pos = np.zeros(3)
        for i in indx:
            pos += data['Atoms'][i][cx:cx+3]
        data['Drawing']['viewPoint'] = [list(pos/len(indx)),[indx[0],0]]
        G2plt.PlotStructure(G2frame,data)

    def OnDistAnglePrt(event):
        'save distances and angles to a file'
        fp = open(os.path.abspath(os.path.splitext(G2frame.GSASprojectfile)[0]+'.disagl'),'w')
        OnDistAngle(event,fp=fp)
        fp.close()

    def OnDistAngleHist(event):
        OnDistAngle(event,hist=True)
        
    def OnSaveDAHist(event):
        Bonds = data.get('Bonds',{})
        Angles = data.get('Angles',{})
        if not len(Bonds) and not len(Angles):
            print('No bond or angle histograms to save!')
            return
        generalData = data['General']
        pName = generalData['Name']
        fName = os.path.abspath(os.path.splitext(G2frame.GSASprojectfile)[0]+'_BAdist.csv')
        fp = open(fName,'w')
        fp.write('Phase name:  %s,'%pName)
        for Atypes in Bonds:
            fp.write('\nBond distances for %s, \n'%Atypes)
            for ib,bond in enumerate(Bonds[Atypes]):
               fp.write('%7.4f,'%bond) 
               if not (ib+1)%10: fp.write('\n')
        for Atypes in Angles:
            fp.write('\nBond angles about %s, \n'%Atypes)
            for ib,angle in enumerate(Angles[Atypes]):
               fp.write('%7.2f,'%angle) 
               if not (ib+1)%10: fp.write('\n')
        fp.close()
        print(' %s written'%fName)
        
    def OnClearDAHist(event):
        data['Bonds'] = {}
        data['Angles'] = {}
        
        print('Bond & Angle histograms now empty')

    def OnDistAngle(event,fp=None,hist=False):
        '''Compute distances and angles in response to a menu command
        or may be called by :func:`OnDistAnglePrt` or :func:`OnDistAngleHist`.

        This calls :func:`GSASIIstrMain.PrintDistAngle` to compute
        bond distances & angles (computed in
        :func:`GSASIIstrMain.RetDistAngle`) and
        then format the results in a convenient way or plot them
        using G2plt.PlotBarGraph.
        '''

        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if not indx:
            print ("select one or more rows of atoms")
            G2frame.ErrorDialog('Select atom',"select one or more rows of atoms then redo")
            return
        Oxyz = []
        xyz = []
        DisAglData = {}
        DisAglCtls = {}
        # add RB stuff to DisAglData
        DisAglData['RBlist'] = []          # list of atom numbers used in the RB
        for d in data['RBModels'].get('Residue',{}):
            for rId in d['Ids']:
                num = int(G2obj.LookupAtomId(ranId=rId,pId=data['pId']))
                DisAglData['RBlist'].append(num)
        DisAglData['rigidbodyDict'] = G2frame.GPXtree.GetItemPyData(
             G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        _,DisAglData['Phases'] = G2frame.GetUsedHistogramsAndPhasesfromTree()
        DisAglData['parmDict'] = G2IO.mkParmDictfromTree(G2frame)
        # other setup prep
        generalData = data['General']
        DisAglData['OrigIndx'] = indx
        if 'DisAglCtls' in generalData:
            DisAglCtls = generalData['DisAglCtls']
        # get atomic radii
        dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData)
        if dlg.ShowModal() == wx.ID_OK:
            DisAglCtls = dlg.GetData()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        generalData['DisAglCtls'] = DisAglCtls
        atomData = data['Atoms']
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        cx = colLabels.index('x')
        cn = colLabels.index('Name')
        ct = colLabels.index('Type')
        Atypes = []
        DisAglData['TargIndx'] = list(range(len(atomData)))
        for i,atom in enumerate(atomData):
            xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
            if i in indx:
                Oxyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
                Atypes.append(atom[ct])
        Atypes = set(Atypes)
        Atypes = ', '.join(Atypes)
        DisAglData['OrigAtoms'] = Oxyz
        DisAglData['TargAtoms'] = xyz
        generalData = data['General']
        DisAglData['SGData'] = generalData['SGData']
        DisAglData['Cell'] = generalData['Cell'][1:] #+ volume
        if 'pId' in data:
            DisAglData['pId'] = data['pId']
            DisAglData['covData'] = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Covariance'))
        try:
            if hist:
                Bonds = data['Bonds'] = data.get('Bonds',{})
                Angles = data['Angles'] = data.get('Angles',{})
                pgbar = wx.ProgressDialog('Distance Angle calculation','Atoms done=',len(Oxyz)+1,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
                AtomLabels,DistArray,AngArray = G2stMn.RetDistAngle(DisAglCtls,DisAglData,pgbar)
                pgbar.Destroy()
                bonds = []
                for dists in DistArray:
                    bonds += [item[3] for item in DistArray[dists]]
                G2plt.PlotBarGraph(G2frame,bonds,Xname=r'$\mathsf{Bonds,\AA}$',
                    Title='Bond distances for %s'%Atypes,PlotName='%s Bonds'%Atypes)
                print('Total number of bonds to %s is %d'%(Atypes,len(bonds)))
                Bonds[Atypes] = bonds
                angles = []
                for Item in AngArray:
                    angles += [item[2][0] for item in AngArray[Item]]
                G2plt.PlotBarGraph(G2frame,angles,Xname=r'$\mathsf{Angles,{^o}}$',
                    Title='Bond angles about %s'%Atypes,PlotName='%s Angles'%Atypes)
                print('Total number of angles about %s is %d'%(Atypes,len(angles)))
                Angles[Atypes] = angles

            elif fp:
                G2stMn.PrintDistAngle(DisAglCtls,DisAglData,fp)
            else:
                G2stMn.PrintDistAngle(DisAglCtls,DisAglData)
        except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
            G2frame.ErrorDialog('Distance/Angle calculation','try again but do "Reset" to fill in missing atom types')

    def OnFracSplit(event):
        'Split atom frac accordintg to atom type & refined site fraction'
        cx,ct,cs,ci = G2mth.getAtomPtrs(data)
        indx = getAtomSelections(Atoms,ct-1)
        if indx:
            B2 = {}
            atomData = data['Atoms']
            PE = G2elemGUI.PickElement(G2frame,oneOnly=True)
            if PE.ShowModal() == wx.ID_OK:
                Atype2 = PE.Elem.strip()
                AtomInfo2 = G2elem.GetAtomInfo(Atype2)
                FF2 = G2elem.GetFFtable([Atype2])[Atype2]
                for iso in AtomInfo2['Isotopes']:
                    B2[iso] = AtomInfo2['Isotopes'][iso]['SL'][0]
            PE.Destroy()
            SQ = 0.0
            PE2 = G2G.SingleFloatDialog(G2frame,'form factor','Enter sinth/lam for atom frac split',SQ)
            if PE2.ShowModal() == wx.ID_OK:
                SQ = PE2.GetValue()
            PE2.Destroy()    
            ff2 = G2elem.ScatFac(FF2,SQ**2)
        print(' X-ray site fractions for sin(th)/lam = %.3f'%SQ)
        for ind in indx:
            Aname1 = atomData[ind][ct-1]
            Atype1 = G2elem.FixValence(atomData[ind][ct])
            if Atype1 == 'Q':
                PE = G2elemGUI.PickElement(G2frame,oneOnly=False)
                if PE.ShowModal() == wx.ID_OK:
                    Atype1 = PE.Elem.strip()
                PE.Destroy()
            Afrac = atomData[ind][cx+3]
            Amult = float(atomData[ind][cx+5])
            AtomInfo1 = G2elem.GetAtomInfo(Atype1)
            if Atype1 == Atype2:
                print('ERROR - 2nd atom type must be different from selected atom')
                continue
            FF1 = G2elem.GetFFtable([Atype1])[Atype1]
            ff1 = G2elem.ScatFac(FF1,SQ**2)
            if ff1 != ff2:
                frac1 = (ff1*Afrac-ff2)/(ff1-ff2)
                print(' For %s: X-ray based site fractions %s = %.3f, %.3f/cell; %s = %.3f, %.3f/cell'     \
                  %(Aname1,Atype1,frac1,frac1*Amult,Atype2,(1.-frac1),(1.-frac1)*Amult))
            
            B1 = AtomInfo1['Isotopes']['Nat. Abund.']['SL'][0]
            for iso in B2:
                if B1 != B2[iso]:
                    bfrac1 = (B1*Afrac-B2[iso])/(B1-B2[iso])
                    atype2 = Atype2
                    if 'nat' not in iso:
                        atype2 += '_%s'%iso
                    print('        neutron based site fractions %s = %.3f, %.3f/cell; %s = %.3f, %.3f/cell'      \
                          %(Atype1,bfrac1,bfrac1*Amult,atype2,(1.-bfrac1),(1.-bfrac1)*Amult))
    
    def OnValidProtein(event):

        def pickHandler(resName):
            drawData = data['Drawing']
            resid = resIDs[resName]
            drawData['viewPoint'][0] = atomData[AtLookUp[resid]][cx:cx+3]
            UpdateDrawAtoms(G2frame,data)
            G2plt.PlotStructure(G2frame,data)

        atomData = data['Atoms']
        cx,ct,cs,cia = data['General']['AtomPtrs']
        AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
        resNames,Probs1,resIDs = G2mth.validProtein(data,True)         #old version
        resNames,Probs2,resIDs = G2mth.validProtein(data,False)        #new version
        print ('Plot 1 is Protein validation based on errat.f')
        print ('Ref: Colovos, C. & Yeates, T.O. Protein Science 2, 1511-1519 (1991).')
        print ('Residue error scores >6 for 5% & >8 for 1% likelihood of being correct')
        print ('NB: this calc. matches errat.f result')
        print ('Plot 2 is Protein validation based on erratv2.cpp; by D. Obukhov & T. Yeates (2002)')
        print ('Residue error scores >11.5 for 5% & >17.2 for 1% likelihood of being correct')
        print ('NB: this calc. gives a close approximate to original erratv2 result')
        G2plt.PlotAAProb(G2frame,resNames,Probs1,Probs2,Title='Error score for %s'%(data['General']['Name']),
            thresh=[[8.0,6.0],[17.191,11.527]],pickHandler=pickHandler)

    def OnShowIsoDistortCalc(event):
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        G2cnstG.ShowIsoDistortCalc(G2frame,data['General']['Name'])

    def OnShowIsoModes(event):
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        #import imp
        #imp.reload(G2cnstG)
        G2cnstG.ShowIsoModes(G2frame,data['General']['Name'])

    def OnReplacePhase(event):
        'Called to replace the current phase with a new phase from a file'
        reqrdr = G2frame.dataWindow.ReplaceMenuId.get(event.GetId())
        rdlist = G2frame.OnImportGeneric(reqrdr,
            G2frame.ImportPhaseReaderlist,'phase')
        if len(rdlist) == 0: return
        # rdlist is only expected to have one entry
        rd = rdlist[0]
        # Strict = True
        # if 'rmc6f' in rd.readfilename:
        #     Strict = False
        #     idx = -1
        phsnam = data['General']['Name']
        # clear out stuff that in general should not be
        # transferred from one phase to another
        resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        Restraints = G2frame.GPXtree.GetItemPyData(resId)
        # resId = G2gd.GetGPXtreeItemId(G2frame,resId,phsnam)
        if phsnam in Restraints:
            del Restraints[phsnam]
        consId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints')
        Constraints = G2frame.GPXtree.GetItemPyData(consId)
        # TODO: would be cleaner to just delete constraints w/phase data['pId']
        Constraints['Phase'] = []
        if 'magPhases' in data: del data['magPhases']
        data['MCSA'] = {'Models':
            [{'Type': 'MD', 'Coef': [1.0, False, [0.8, 1.2]], 'axis': [0, 0, 1]}],
            'Results': [], 'AtInfo': {}}
        data['RMC'] = {'RMCProfile': {}, 'fullrmc': {}, 'PDFfit': {}}
        data['ISODISTORT'] = {}
        data['Deformations'] = {}
        # copy over most of the data from the reader object
        # but keep original name, pId & ranId
        for key in ['General', 'Atoms', 'Drawing', 'Histograms', 'Pawley ref', 'RBModels']:
            data[key] = rd.Phase[key]
        # TODO: want to keep Histograms rather than wipe it out?
        # If so, need to think about Dij & microstrain:
            # UseList = newPhase['Histograms']
            # for hist in UseList:
            #     # reset Dij & microstrain terms where # of terms changes
            #     if len(UseList[hist]['Mustrain'][4]) != NShkl:
            #         UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
            #     if len(UseList[hist]['HStrain'][0]) != NDij:
            #         UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]

        # restore existing phase name
        newname = rd.Phase['General']['Name']
        data['General']['Name'] = phsnam
        # rename phase to new name from file
        renamePhaseName(G2frame,data,G2frame.PickId, data['General'],newname)
        # force a reload of current tree item
        G2frame.PickIdText = []
        wx.CallAfter(G2gd.SelectDataTreeItem,G2frame,G2frame.PickId)

    def OnReImport(event):
        'Called to replace the coordinates with "original" values from a file'
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        reqrdr = G2frame.dataWindow.ReImportMenuId.get(event.GetId())
        rdlist = G2frame.OnImportGeneric(reqrdr,
            G2frame.ImportPhaseReaderlist,'phase')
        if len(rdlist) == 0: return
        # rdlist is only expected to have one element
        rd = rdlist[0]
        Strict = True
        if 'rmc6f' in rd.readfilename:
            Strict = False
            idx = -1
        G2frame.OnFileSave(event)
        # rd contains all info for a phase
        PhaseName = rd.Phase['General']['Name']
        print ('Read phase '+str(PhaseName)+' from file '+str(G2frame.lastimport))
        atomData = data['Atoms']
        atomNames = []
        All = False
        for atom in atomData:
            atomNames.append(''.join(atom[:ct+1]).capitalize())  #eliminate spurious differences
        for atom in rd.Phase['Atoms']:
            try:
                if Strict:
                    idx = atomNames.index(''.join(atom[:ct+1]).capitalize())  #eliminate spurious differences
                else:
                    idx += 1
                atId = atomData[idx][cia+8]                                 #save old Id
                atomData[idx][:cia+8] = atom[:cia+8]+[atId,]
            except ValueError:
                if All:
                    atomData.append(atom)
                else:
                    dlg = wx.MessageDialog(G2frame,'Some atoms not in List; do you want to append them all',   \
                        'Unknown atom '+atom[0],wx.YES_NO|wx.ICON_QUESTION)
                    try:
                        result = dlg.ShowModal()
                        if result in [wx.ID_YES,]:
                            All = True
                            atomData.append(atom)
                        else:
                            print (atom[:ct+1]+ 'not in Atom array; not updated')
                    finally:
                        dlg.Destroy()
        SetupGeneral()
        OnReloadDrawAtoms(event)
        wx.CallAfter(FillAtomsGrid,Atoms)

#### Dysnomia (MEM) Data page ##############################################################################

    def OnLoadDysnomia(event):
        print('Load MEM - might not be implemented')

    def OnSaveDysnomia(event):
        print('Save MEM - might not be implemented')

    def OnRunDysnomia(event):

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
        treeId = G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists')
        if treeId:
            reflSets = G2frame.GPXtree.GetItemPyData(treeId)    
            reflData = reflSets[generalData['Name']]['RefList']
        else:
            wx.MessageBox('You must have PWDR reflections before running Dysnomia','Dysnomia Error',
                style=wx.ICON_ERROR)
            return
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
        wx.MessageBox(' For use of Dysnomia, please cite:\n\n'+G2G.GetCite('Dysnomia'),
            caption='Dysnomia (MEM)',style=wx.ICON_INFORMATION)

        print('Run '+DYSNOMIA)
        subp.call([DYSNOMIA,prfName])

        DysData['DenStart'] = 'uniform'
        DysData['prior'] = 'uniform'
        wx.CallAfter(G2phsG2.UpdateDysnomia,G2frame,data)

        goon,reflData = G2pwd.MEMupdateReflData(prfName,data,reflData)
        if goon:
            reflSets[generalData['Name']]['RefList'] = reflData
            G2frame.GPXtree.SetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'Reflection Lists'),reflSets)
            OnFourierMaps(event)           #auto run Fourier
            if DysData['clear']:
                os.remove(os.path.splitext(prfName)[0]+'.fba')
                os.remove(os.path.splitext(prfName)[0]+'.mem')
                os.remove(os.path.splitext(prfName)[0]+'.out')
                os.remove(os.path.splitext(prfName)[0]+'.prf')
                os.remove(os.path.splitext(prfName)[0]+'_eps.raw')
        else:
            wx.MessageBox('Dysnomia failed to make new structure factors','Dysnomia Error',
                style=wx.ICON_ERROR)

#### RMC Data page ################################################################################
# fullrmc stuff TODO:
#  1) need to implement swapping in scripts
#  2) fullrmc tutorials

    def OnSetupRMC(event):
        written = lambda fil: print(f' {fil} written')
        generalData = data['General']
        if not G2frame.GSASprojectfile:     #force a project save
            G2frame.OnFileSaveas(event)
        dName = G2frame.LastGPXdir
        os.chdir(dName)
        print(f'Writing input files in directory {dName!r}')
        if G2frame.RMCchoice == 'fullrmc':
            RMCPdict = data['RMC']['fullrmc']
            pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
            pName = pName.replace(' ','_')
            G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
            if RMCPdict['Swaps']:
                wx.MessageDialog(G2frame, G2G.StripIndents(
                        '''GSAS-II does not yet fully support use of swapping in fullrmc.
                        Edit the script by hand before using.''',True),
                        'No swaps yet',wx.OK).ShowModal()
            #--------- debug stuff
            # if GSASIIpath.GetConfigValue('debug'):
            #     print('reloading',G2pwd)
            #     import imp
            #     imp.reload(G2pwd)
            #--------- end debug stuff
            rname = G2pwd.MakefullrmcRun(pName,data,RMCPdict)
            print('build of fullrmc file {} completed'.format(rname))
        elif G2frame.RMCchoice == 'RMCProfile':
            if ' ' in generalData['Name']:
                wx.MessageDialog(G2frame,'ERROR: Phase name has space; change phase name','Bad phase name',wx.ICON_ERROR).ShowModal()
                G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
                return
            dName = G2frame.LastGPXdir
            pName = generalData['Name']
            G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
            RMCPdict = data['RMC']['RMCProfile']
            PWId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,RMCPdict['histogram'][0])
            if PWId:
                PWDdata = G2frame.GetPWDRdatafromTree(PWId)
                histoName = G2frame.GPXtree.GetItemPyData(PWId)[2]
                Size = data['Histograms'][histoName]['Size']
                Mustrain = data['Histograms'][histoName]['Mustrain']
                reset = False
                written(G2pwd.MakeInst(PWDdata,pName,Size,Mustrain,RMCPdict['UseSampBrd']))
                backfile = G2pwd.MakeBack(PWDdata,pName)
                if backfile is None:
                    print(' Chebyschev-1 background not used; no .back file written')
                    wx.MessageDialog(G2frame,' Chebyschev-1 background not used; '+ \
                        'no .back file written & RMCProfile will not run','Wrong background function',wx.OK).ShowModal()
                    G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
                    return
                else:
                    written(backfile)
                written(G2pwd.MakeBragg(PWDdata,pName,data))
                if RMCPdict['ReStart'][0]:
                    if os.path.isfile(pName+'.his6f'):
                        os.remove(pName+'.his6f')
                    RMC6f,reset = G2pwd.MakeRMC6f(PWDdata,pName,data,RMCPdict)
                    written(RMC6f)
                fname = G2pwd.MakeRMCPdat(PWDdata,pName,data,RMCPdict)
                if 'Error' in fname:
                    print(fname)
                    wx.MessageDialog(G2frame,fname,'Missing reflection list',wx.OK).ShowModal()
                    G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
                    return
                written(fname)
                print('RMCProfile file build completed')
                RMCPdict['ReStart'] = [False,False]
                if reset:
                    wx.MessageDialog(G2frame,' Vacancies found & "Va" atoms added to list. '+ \
                        'You may need to revise RMCProfile setup parameters.','Repeat Setup RMC',wx.OK).ShowModal()
                    wx.CallAfter(G2rmcG.UpdateRMC,G2frame,data)
            else:
                print('RMCProfile file build failed - no histogram selected')
                G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
        elif G2frame.RMCchoice == 'PDFfit':
            if ' ' in generalData['Name']:
                wx.MessageDialog(G2frame,'ERROR: Phase name has space; change phase name','Bad phase name',wx.ICON_ERROR).ShowModal()
                G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,False)
                return
            G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC,True)
            RMCPdict = data['RMC']['PDFfit']
            msg = G2pwd.MakePDFfitAtomsFile(data,RMCPdict)
            if msg:
                G2G.G2MessageBox(G2frame,'ERROR: '+msg,'PDFfit setup failure')
                return
            fname = G2pwd.MakePDFfitRunFile(data,RMCPdict)
            if fname is None:
                wx.MessageDialog(G2frame,'ERROR: failure to setup PDFfit; check console','PDFfit setup failure',wx.ICON_ERROR).ShowModal()
            else:
                written(fname)
                print('PDFfit file build completed')

    def RunPDFfit(event):
        generalData = data['General']
        ISOdict = data['ISODISTORT']
        PDFfit_exec = G2pwd.findPDFfit()  #returns location of python with PDFfit installed
        if not PDFfit_exec:
            wx.MessageBox(''' PDFfit2 is not currently installed for this platform.
    Please contact us for assistance''',caption='No PDFfit2',style=wx.ICON_INFORMATION)
            return
        RMCPdict = data['RMC']['PDFfit']
        pName = generalData['Name'].replace(' ','_')
        if 'sequential' in RMCPdict['refinement']:
            rname = 'Seq_PDFfit.py'
        else:
            rname = pName+'-PDFfit.py'
            if not os.path.exists(rname):
                wx.MessageBox(f'File {rname} does not exist. Has the Operations/"Setup RMC" menu command been run?',
                                  caption='Run setup',style=wx.ICON_WARNING)
                return
        wx.MessageBox(' For use of PDFfit2, please cite:\n\n'+
                          G2G.GetCite('PDFfit2'),
                          caption='PDFfit2',style=wx.ICON_INFORMATION)
        G2frame.OnFileSave(event)
        print (' GSAS-II project saved')
        if sys.platform.lower().startswith('win'):
            batch = open('pdffit2.bat','w')
            # Include an activate command here
            p = os.path.split(PDFfit_exec)[0]
            while p:
                if os.path.exists(os.path.join(p,'Scripts','activate')):
                    batch.write('call '+os.path.join(p,'Scripts','activate')+'\n')
                    break
                prevp = p
                p = os.path.split(p)[0]
                if prevp == p:
                    print('Note, no activate command found')
                    break
            batch.write(PDFfit_exec+' '+rname+'\n')
            # batch.write('pause')
            if 'normal' in RMCPdict['refinement']:
                batch.write('pause')
            batch.close()
        else:
            batch = open('pdffit2.sh','w')
            batch.write('#!/bin/bash\n')
            # include an activate command here
            p = os.path.split(PDFfit_exec)[0]
            while p:
                if os.path.exists(os.path.join(p,'bin','activate')):
                    batch.write('source '+os.path.join(p,'Scripts','activate')+'\n')
                    break
                prevp = p
                p = os.path.split(p)[0]
                if prevp == p:
                    print('Note, no activate command found')
                    break

            batch.write('cd ' + os.path.split(os.path.abspath(rname))[0] + '\n')
            batch.write(PDFfit_exec + ' ' + os.path.abspath(rname) + '\n')
            batch.close()
        if 'sequential' in RMCPdict['refinement']:
            Id =  G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Sequential PDFfit2 results')
#            if Id:
#                saveSeqResult = G2frame.GPXtree.GetItemPyData(Id)
#            else:
            if not Id:
                SeqResult = {}
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Sequential PDFfit2 results')
            G2Names = [item.name for item in ISOdict['G2ModeList']]
            SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
            SeqResult['histNames'] = []         #this clears the previous seq. result!
            SeqNames = []
            for itm in range(len(RMCPdict['seqfiles'])):
                SeqNames.append([itm,RMCPdict['seqfiles'][itm][0]])
            if RMCPdict['SeqReverse']:
                SeqNames.reverse()
            nPDF = len(SeqNames)
            pgbar = wx.ProgressDialog('Sequential PDFfit','PDF G(R) done = 0',nPDF+1,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
            newParms = {}
            for itm,Item in enumerate(SeqNames):
                PDFfile = RMCPdict['seqfiles'][Item[0]]
                pfdata = PDFfile[1]['G(R)'][1].T
#                    pfname = PDFfile[0].replace(' ','_')
                pfname = 'Seq_PDF.gr'
                pfile = open(pfname,'w')
                for dp in pfdata:
                    pfile.write('%12.5f%12.5f\n'%(dp[0],dp[1]))
                pfile.close()
                rfile = open('Seq_PDFfit_template.py','r')
                lines = rfile.readlines()       #template lines
                rfile.close()
                newlines = []
                parms = {}
                Np = 0
                for line in lines:
                    if '#sequential' in line:
                        newlines += "pf.read_data('%s', '%s', 30.0, %.4f)\n"%(pfname,PDFfile[1]['Type'][0],PDFfile[1]['qdamp'][0])
                        newlines += 'pf.setdata(1)\n'
                        newlines += 'pf.pdfrange(1, %6.2f, %6.2f)\n'%(PDFfile[1]['Fitrange'][0],PDFfile[1]['Fitrange'][1])
                        for item in ['dscale','qdamp','qbroad']:
                            if PDFfile[1][item][1]:
                                Np += 1
                                newlines += 'pf.constrain(pf.%s(),"@%d")\n'%(item,Np)
                                parms[item] = '%d'%Np
                                if itm and RMCPdict['SeqCopy']:
                                    newParms[parms[item]] = RMCPdict['Parms'][parms[item]]
                                else:
                                    if not itm and 'result' not in PDFfile[1]:
                                        newParms[parms[item]] = PDFfile[1][item][0]
                                    else:
                                        newParms[parms[item]] = PDFfile[1]['result'][parms[item]][0]
                    elif '#parameters' in line:
                        startParms = RMCPdict['Parms']
                        if newParms or RMCPdict['SeqCopy']:
                            if newParms:
                                startParms = newParms
                            for iprm in startParms:
                                if int(iprm) > 9:
                                    break
                                newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,startParms[iprm])
                            print('Begin dscale: %d %.4f'%(itm,startParms['1']))
                            for iprm in RMCPdict['Parms']:
                                if isinstance(RMCPdict['Parms'][iprm],float):
                                    newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,RMCPdict['Parms'][iprm])
                                else:
                                    newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,RMCPdict['Parms'][iprm][0])
                        elif not RMCPdict['SeqCopy']:
                            startParms = PDFfile[1]['result']
                            for iprm in startParms:
                                newlines += 'pf.setpar(%s,%.6f)\n'%(iprm,startParms[iprm][0])
                            print('Begin dscale: %d %.4f'%(itm,startParms['1']))
                    else:
                        newlines += line
                rfile= open('Seq_PDFfit.py','w')
                rfile.writelines(newlines)
                rfile.close()
                fName = 'Sequential_PDFfit'     #clean out old PDFfit output files
                if os.path.isfile(fName+'.res'):
                    os.remove(fName+'.res')
                if os.path.isfile(fName+'.rstr'):
                    os.remove(fName+'.rstr')
                if os.path.isfile(fName+'.fgr'):
                    os.remove(fName+'.fgr')

                if sys.platform.lower().startswith('win'):
                    Proc = subp.Popen('pdffit2.bat',creationflags=subp.CREATE_NEW_CONSOLE)
                    Proc.wait()     #for it to finish before continuing on
                else:
                    if sys.platform == "darwin":
                        GSASIIpath.MacRunScript(os.path.abspath('pdffit2.sh'))
                    else:
                        Proc = subp.Popen(['/bin/bash','pdffit2.sh'])
                        Proc.wait()

                newParms,Rwp =  G2pwd.UpdatePDFfit(data,RMCPdict)
                if isinstance(newParms,str):
                    wx.MessageBox('Singular matrix in PDFfit',caption='PDFfit2 failed',style=wx.ICON_INFORMATION)
                    break
                for item in ['dscale','qdamp','qbroad']:
                    if PDFfile[1][item][1]:
                        PDFfile[1][item][0] = newParms[parms[item]][0]
                PDFfile[1]['result'] = copy.deepcopy(newParms)
                parmDict = copy.deepcopy(newParms)
                parmDict.update({'Temperature':PDFfile[1]['Temp']})
                tempList = ['%s-%s'%(parms[item],item) for item in parms]       #these come first
                parmkeys = [int(item) for item in RMCPdict['ParmNames']]
                parmkeys.sort()
                tempList += ['%s-%s'%(item,RMCPdict['ParmNames'][item]) for item in parmkeys]
                print('result dscale: ',parmDict['1'],' Rw: ',Rwp)
                atParms = [str(i+21) for i in range(len(G2Names))]
                varyList = []
                for item in tempList:
                    pid = item.split('-')[0]
                    if pid in atParms:
                        item = '%s-%s'%(pid,G2Names[int(pid)-21])
                    varyList.append(item)
                result = np.array(list(newParms.values())).T
                SeqResult[PDFfile[0]] = {'variables':result[0],'varyList':varyList,'sig':result[1],'Rvals':{'Rwp':Rwp,},
                    'covMatrix':[],'title':PDFfile[0],'parmDict':parmDict}

                pfile = open('Sequential_PDFfit.fgr')
                XYcalc = np.loadtxt(pfile).T[:2]
                pfile.close()
                pId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PDFfile[0])
                PDFctrl = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,pId,'PDF Controls'))
                XYobs = PDFctrl['G(R)'][1]
                if XYobs.shape[0] < 4:
                    XYobs = np.concatenate((XYobs,np.zeros_like(XYobs)),axis=0)
                ibeg = np.searchsorted( XYobs[0],XYcalc[0][0])
                ifin = ibeg+XYcalc.shape[1]
                XYobs[2][ibeg:ifin] = XYcalc[1]
                XYobs[3] = XYobs[1]-XYobs[2]
                PDFctrl['G(R)'][1] = XYobs
                SeqResult['histNames'].append(Item[1])
                GoOn = pgbar.Update(itm,newmsg='PDF G(R) done = %d'%(itm))
                if not GoOn[0]:
                    print(' Sequential PDFfit aborted')
                    break

            pgbar.Destroy()
            G2frame.GPXtree.SetItemPyData(Id,SeqResult)
            G2frame.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
            G2frame.GPXtree.SelectItem(Id)

        else: #normal
            #remove any old PDFfit output files
            fName = generalData['Name'].replace(' ','_')+'-PDFfit'
            if os.path.isfile(fName+'.res'):
                os.remove(fName+'.res')
            if os.path.isfile(fName+'.rstr'):
                os.remove(fName+'.rstr')
            if os.path.isfile(fName+'N.fgr'):
                os.remove(fName+'N.fgr')
            if os.path.isfile(fName+'X.fgr'):
                os.remove(fName+'X.fgr')

            if sys.platform.lower().startswith('win'):
                Proc = subp.Popen('pdffit2.bat',creationflags=subp.CREATE_NEW_CONSOLE)
                Proc.wait()     #for it to finish before continuing on
            else:
                if sys.platform == "darwin":
                    GSASIIpath.MacRunScript(os.path.abspath('pdffit2.sh'))
                else:
                    Proc = subp.Popen(['/bin/bash','pdffit2.sh'])
                    Proc.wait()     #for it to finish before continuing on
            #update choice? here?
            dlg = wx.MessageDialog(G2frame,'Check PDFfit console for results; do you want to update?',
                'PDFfit run finished',wx.YES|wx.NO)
            try:
                dlg.CenterOnParent()
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result == wx.ID_YES:
                Error =  G2pwd.UpdatePDFfit(data,RMCPdict)
                if Error:
                    wx.MessageBox('PDFfit failed',caption='%s not found'%Error[0],style=wx.ICON_EXCLAMATION)
            G2rmcG.UpdateRMC(G2frame,data)

    def Runfullrmc(event):
        fullrmc_exec = G2pwd.findfullrmc()
        if fullrmc_exec is None:
            G2G.G2MessageBox(G2frame,'fullrmc Python not found. How did we get here?')
            return
        generalData = data['General']
        pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
        pName = pName.replace(' ','_')
        rname = pName+'-fullrmc.py'
        if not os.path.exists(rname):
            G2G.G2MessageBox(G2frame,'The fullrmc script has not been created. Running setup.',
                'Not setup')
            OnSetupRMC(event)
        RMCPdict = data['RMC']['fullrmc']
        rmcname = pName+'-fullrmc.rmc'
        if os.path.isdir(rmcname) and RMCPdict['ReStart'][0]:
            msg = '''You have asked to start a new fullrmc run rather than
                 continue the existing {} run.
                 %%Press "Yes" to continue, deleting this
                 previous run or "No" to change the restart checkbox to
                 continue from the previous results.'''.format(rmcname)

            dlg = wx.MessageDialog(G2frame,G2G.StripIndents(msg,True),
                'Restart or continue',wx.YES|wx.NO)
            try:
                dlg.CenterOnParent()
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result == wx.ID_YES:
                shutil.rmtree(rmcname)
            else:
                return
        G2G.G2MessageBox(G2frame,'For use of fullrmc, please cite:\n\n'+
                             G2G.GetCite('fullrmc')+
                             '\n\nNote: A more advanced version of fullrmc can be found at www.fullrmc.com',
                             'Please cite fullrmc')
        ilog = 0
        while True:
            logname = '%s_%d.log'%(pName,ilog)
            if os.path.isfile(logname):
                if GSASIIpath.GetConfigValue('debug'):
                    print('removing',logname)
                os.remove(logname)
            else:
                break
            ilog += 1
        if sys.platform.lower().startswith('win'):
            batch = open('fullrmc.bat','w')
            #batch.write('CALL '+sys.exec_prefix+'\\Scripts\\activate\n')
            batch.write(fullrmc_exec+' '+rname+'\n')
            batch.write('pause')
            batch.close()
            Proc = subp.Popen('fullrmc.bat',creationflags=subp.CREATE_NEW_CONSOLE)
            Proc.wait()     #for it to finish before continuing on
        else:
            batch = open('fullrmc.sh','w')
            batch.write('#!/bin/bash\n')
            #activate = os.path.split(os.environ.get('CONDA_EXE',''))[0] +'/activate'
            batch.write('cd ' + os.path.split(os.path.abspath(rname))[0] + '\n')
            #if os.path.exists(activate):
            #    batch.write('source ' + activate + ' ' +
            #                os.environ['CONDA_DEFAULT_ENV'] +'\n')
            #    batch.write('python ' + rname + '\n')
            #else:
            #    batch.write(sys.exec_prefix+'/python ' + rname + '\n')
            batch.write(fullrmc_exec + ' ' + os.path.abspath(rname) + '\n')
            batch.close()
            if sys.platform == "darwin":
                GSASIIpath.MacRunScript(os.path.abspath('fullrmc.sh'))
            else:
                Proc = subp.Popen(['/bin/bash','fullrmc.sh'])
#                Proc.wait()     #for it to finish before continuing on
        G2rmcG.UpdateRMC(G2frame,data)

    def RunRMCProfile(event):
        generalData = data['General']
        pName = generalData['Name'].replace(' ','_')
        rmcfile = G2pwd.findrmcprofile()
        
        wx.MessageBox(
            ' For use of RMCProfile, please cite:\n\n'+
            G2G.GetCite("RMCProfile"),
            caption='RMCProfile',style=wx.ICON_INFORMATION)
        if os.path.isfile(pName+'.his6f'):
            os.remove(pName+'.his6f')
        if os.path.isfile(pName+'.xray'):
            os.remove(pName+'.xray')
        if os.path.isfile(pName+'.neigh'):
            os.remove(pName+'.neigh')
        if os.path.isfile(pName+'.bonds'):
            os.remove(pName+'.bonds')
        if os.path.isfile(pName+'.triplets'):
            os.remove(pName+'.triplets')
        i = 1
        while True:
            if os.path.isfile(pName+'.bondodf_%d'%i):
                os.remove(pName+'.bondodf_%d'%i)
                os.remove(pName+'_bondplot_%d.ppm'%i)
                i += 1
            else:
                break
        i = 1
        while True:
            if os.path.isfile(pName+'_anglehist_%d.csv'%i):
                os.remove(pName+'_anglehist_%d.csv'%i)
                i += 1
            else:
                break

        G2frame.OnFileSave(event)
        print ('GSAS-II project saved')
        pName = generalData['Name'].replace(' ','_')

        if rmcfile is None:
            wx.MessageBox('''RMCProfile is not correctly installed for use in GSAS-II
        This software must be downloaded separately (from 
        https://rmcprofile.ornl.gov/download). Install the rmcprofile or 
        rmcprofile.exe file in a location where GSAS-II can find it 
        (see config variable rmcprofile_exec in preferences.)''',
            caption='RMCProfile',style=wx.ICON_INFORMATION)
            return

        if sys.platform == "darwin":
            script_file = os.path.join(os.getcwd(), "runrmc.sh")
            with open(script_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f'cd "{os.getcwd()}"\n')
                f.write(f'export PATH="{os.path.dirname(rmcfile)}":$PATH\n')
                f.write(f'"{rmcfile}" "{pName}"\n')
            os.system("chmod +x runrmc.sh")
            ascript_file = os.path.join(os.getcwd(), "runrmc.script")
            with open(ascript_file, 'w') as f:
                f.write('tell application "Terminal"\n')
                f.write(f'''  do script "echo 'Running RMCprofile'"\n''')
                f.write(f'  do script "bash {script_file}" in window 1\n')
                f.write("end tell\n")
            subp.Popen(['osascript', ascript_file])
        elif sys.platform.startswith("linux"):
            script_file = os.path.join(os.getcwd(), "runrmc.sh")
            with open(script_file, 'w') as frmc:
                frmc.write("#!/bin/bash\n")
                frmc.write(f"export LD_LIBRARY_PATH={Path(rmcfile).parent.parent}/exe/libs\n")
                frmc.write(f"export LIBRARY_PATH={Path(rmcfile).parent.parent}/exe/libs\n")
                frmc.write(f"export PATH=$PATH:{Path(rmcfile).parent.parent}/exe\n")
                frmc.write(f'"rmcprofile" "{pName}"\n')
            subp.Popen(['bash', script_file])
        else:
            script_file = os.path.join(os.getcwd(), "runrmc.bat")
            with open(script_file,'w') as batch:
                batch.write('Title RMCProfile\n')   # BHT: is Title a Windows command?
                batch.write(f'"{rmcfile}" "{pName}"\n')
                batch.write('pause\n')
                batch.close()
            subp.Popen(script_file,creationflags=subp.CREATE_NEW_CONSOLE)
#        Proc.wait()     #for it to finish before continuing on
        G2rmcG.UpdateRMC(G2frame,data)

    def OnRunRMC(event):
        '''Run a previously created RMCProfile/fullrmc/PDFfit2 script
        '''
        if G2frame.RMCchoice == 'fullrmc':
             Runfullrmc(event)
        elif G2frame.RMCchoice == 'RMCProfile':
            RunRMCProfile(event)
        elif G2frame.RMCchoice == 'PDFfit':
            RunPDFfit(event)

    # def OnStopRMC(event):
    #     if G2frame.RMCchoice == 'fullrmc':
    #         generalData = data['General']
    #         pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
    #         pName = pName.replace(' ','_')
    #         engineFilePath = pName+'.rmc'
    #         if not os.path.exists(engineFilePath):
    #             print('fullrmc repository {} not found'.format(engineFilePath))
    #             return
    #         try:
    #             from fullrmc import InterceptHook
    #             hook = InterceptHook(path=engineFilePath)
    #             hook.stop_engine()
    #             print('hook.stop_engine() sent to {}'.format(engineFilePath))
    #         except Exception as msg:
    #             print('failed, msg=',msg)

    def OnLoadRMC(event):
        '''Used to load the output from fullrmc with all atoms placed in the
        original cell
        '''
        fullrmcLoadPhase(super=False)
    def OnLoadRMCsuper(event):
        '''Used to load the output from fullrmc with atoms in the simulation
        supercell cell
        '''
        fullrmcLoadPhase(super=True)
    def fullrmcLoadPhase(super):
        '''Used to load the output from fullrmc. Creates a new phase,
        reads all atoms & converts coordinates to fractional.
        If super is False all atoms placed in the original cell.

        Called from :func:`OnLoadRMC` or :func:`OnLoadRMCsuper` from
        the RMC tab Operations menu commands 'Superimpose into cell'
        and 'Load Supercell'.
        '''
        if G2frame.RMCchoice != 'fullrmc':
            print('fullrmcLoadPhase: How did this happen?')
            return
        # create a new phase
        phId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
        phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
        phaseNameList = list(usedHistograms.keys()) # phase names in use
        PhaseName = G2obj.MakeUniqueLabel(data['General']['Name']+'_fullrmc',phaseNameList)
        psub = G2frame.GPXtree.AppendItem(parent=phId,text=PhaseName)
        E,SGData = G2spc.SpcGroup('P 1')
        G2frame.GPXtree.SetItemPyData(psub,G2obj.SetNewPhase(Name=PhaseName,SGData=SGData))
        newPhase = G2frame.GPXtree.GetItemPyData(psub)
        # read in info from file
        pName = (G2frame.GSASprojectfile.split('.')[0] + '-'
                     + data['General']['Name'])
        pName = pName.replace(' ','_')
        try:
            with open(pName+'-fullrmc.atoms','r') as fp:
                cell = [float(i) for i in fp.readline().split(':')[1].split()]
                supercell = [int(i) for i in fp.readline().split(':')[1].split()]
                if super:
                    for i in range(3): cell[i] *= supercell[i]
                # set cell & volume
                newPhase['General']['Cell'][1:7] = cell
                newPhase['General']['Cell'][7] = G2lat.calc_V(G2lat.cell2A(cell))
                A,B = G2lat.cell2AB(cell)
                # add atoms
                for line in fp:
                    Name,El,ox,oy,oz = line.split()
                    oxyz = [float(i) for i in (ox,oy,oz)]
                    if super:
                        (x,y,z) = np.inner(B,oxyz)
                    else:
                        (x,y,z),disp = G2spc.MoveToUnitCell(np.inner(B,oxyz))
                    atId = ran.randint(0,sys.maxsize)
                    newPhase['Atoms'].append([Name,El,'',x,y,z,1.,'1',1,'I',0.01,0,0,0,0,0,0,atId])
        except:
            G2G.G2MessageBox(G2frame,
                    'Unable to open or read file '+pName+'-fullrmc.atoms. '
                    'Was a fullrmc run from the current .gpx file '
                    'and for the current phase?',
                    'Error on read')
            G2frame.GPXtree.Delete(psub)
            return
        # add a restraint tree entry for new phase
        subr = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        G2frame.GPXtree.GetItemPyData(subr).update({PhaseName:{}})
        SetupGeneral()  # index elements

        #wx.CallAfter(G2frame.GPXtree.SelectItem,psub) # should call SelectDataT

    def OnViewRMC(event):
        if G2frame.RMCchoice == 'fullrmc':
            RMCPdict = data['RMC']['fullrmc']
            generalData = data['General']
            pName = G2frame.GSASprojectfile.split('.')[0] + '-' + generalData['Name']
            pName = pName.replace(' ','_')
            engineFilePath = pName+'-fullrmc.rmc'
            if not os.path.exists(engineFilePath):
                dlg = wx.FileDialog(G2frame, 'Open fullrmc directory',
                    defaultFile='*.rmc',wildcard='*.rmc')
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        engineFilePath = dlg.GetPath()
                    else:
                        return
                finally:
                    dlg.Destroy()
                engineFilePath = os.path.splitext(engineFilePath)[0] + '.rmc'
                if not os.path.exists(engineFilePath): return
            choices = []
            statFilePath = os.path.splitext(engineFilePath)[0] + '.stats'
            plotFilePath = os.path.splitext(engineFilePath)[0] + '.plots'
            imgDict = {}
            if os.path.exists(statFilePath):
                fp = open(statFilePath,'r')
                vals = []
                for i,line in enumerate(fp):
                    v = line.strip().split(',')[:-1] # ends with comma, remove last empty element
                    if i == 0:
                        lbls = [i.strip() for i in v[1:]]
                        continue
                    try:
                        vals.append([float(i) for i in v])
                    except:
                        print('Error reading line ',i,'in',statFilePath)
                fp.close()
                steps = np.array(vals)[:,0]
                yvals = np.array(vals)[:,1:].T
                choices = ['Constraints vs. Steps']
            if os.path.exists(plotFilePath):
                import pickle
                fp = open(plotFilePath,'rb')
                while True:
                    try:
                        title = pickle.load(fp)
                        imgDict[title] = fp.tell()
                        im = pickle.load(fp)
                        choices += [title]
                    except:
                        break
                fp.close()
            if not choices:
                G2G.G2MessageBox(G2frame,
                    'Nothing to plot. '+
                    'No results in '+statFilePath+' or '+plotFilePath,
                    'Nothing to plot')
                return
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Select plots to see displayed.',
                                              'Select plots',choices)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_OK:
                    selectedPlots = [choices[i] for i in dlg.GetSelections()]
                else:
                    return
            finally:
                dlg.Destroy()

            for plt in selectedPlots:
                if plt in imgDict:
                    fp = open(plotFilePath,'rb')
                    fp.seek(imgDict[plt])
                    try:
                        im = pickle.load(fp)
                        G2plt.PlotRawImage(G2frame,im,plt)
                    except:
                        pass
                    fp.close()
                else:
                    plotLbls = []
                    plotVals = []
                    for lbl,row in zip(lbls,yvals): # deal with <=0 costs
                        if sum(row**2) == 0: continue # drop if all zeros
                        if min(row) <= 0:
                            row = np.where(row>0,row,min(row[np.where(row>0)])/10.)
                        plotLbls.append(lbl)
                        plotVals.append([steps,np.log10(row)])
                    title = 'fullrmc residuals for '+pName
                    G2plt.PlotXY(G2frame,plotVals,
                                labelX='generated steps',
                                labelY=r'$log_{10}$ ($\mathsf{\chi^2})$',
                                newPlot=True,Title=title,
                                lines=True,names=plotLbls)
            return
        elif G2frame.RMCchoice == 'RMCProfile':
            generalData = data['General']
            RMCPdict = data['RMC']['RMCProfile']
            pName = generalData['Name'].replace(' ','_')
            dlg = wx.FileDialog(G2frame, "Choose any RMCProfile csv results file for "+pName+":",
                defaultDir=G2frame.LastGPXdir,style=wx.FD_CHANGE_DIR,wildcard='RMCProfile result csv files|'+pName+'*.csv')
            if dlg.ShowModal() == wx.ID_OK:
                path = os.path.split(dlg.GetPath())[0]
                dlg.Destroy()
            else:
                dlg.Destroy()
                return

            ifXray = False
            ifNeut = False
            try:
                datFile = open(os.path.join(path,pName+'.dat'),'r')
                datLines = datFile.readlines()
                datFile.close()
                for line in datLines:
                    if 'xray' in line:
                        ifXray = True
            except:
                pass
            files =  {'_PDF1.csv':[],'_PDF2.csv':[],'_PDFpartials.csv':[],'_SQ1.csv':[],'_SQ2.csv':[],'_XFQ1.csv':[],
                      '_SQ1partials.csv':[],'_SQ2partials.csv':[],'_FQ1.csv':[],'_FT_XFQ1.csv':[],
                      '_FQ1partials.csv':[],'_bragg.csv':[],'.chi2':[]}
            for item in files:
                if os.path.exists(os.path.join(path,pName+item)):
                    OutFile = open(pName+item,'r')
                    files[item] = OutFile.readlines()
                    OutFile.close()
                    print('RMCProfile file %s read'%(pName+item))
#                else:
#                    print('RMCProfile file %s not found'%(pName+item))
#total result plots
            Labels = {'_PDF1.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R) for '],
                '_PDF2.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R)-2 for '],
                '_SQ1.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q) for '],
                '_SQ2.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q)-2 for '],
                '_FQ1.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP x-ray S(Q) for '],
                '_FT_XFQ1.csv':[r'$\mathsf{R,\AA}$','g(R)','RMCP x-ray g(R) for '],
                '_bragg.csv':[r'$\mathsf{TOF,\mu s}$','Normalized Intensity','RMCP bragg for ']}
            Ysave = []
            for label in Labels:
                X = []
                Yobs = []
                Ycalc = []
                if len(files[label]):
                    if 'PDF1' in label:
                        ifNeut = True
                    Names = files[label][0][:-1].split(',')
                    Xmax = 100.
                    if 'XFQ' in label:
                        Xmax = RMCPdict.get('Rmax',100.)
                    for line in files[label][1:]:
                        items = line.split(',')
                        if 'XFQ' in label and float(items[0]) > Xmax:
                            break
                        X.append(float(items[0]))
                        Yobs.append(float(items[2]))
                        Ycalc.append(float(items[1]))
                    Yobs = np.array([X,Yobs])
                    Ycalc = np.array([X,Ycalc])
                    if '(R)' in Labels[label][1]:
                        Ysave.append(Ycalc)
                        Ymin = Ysave[0][1][0]
                    if 'bragg' in label:
                        Ydiff = np.array([X,(Yobs-Ycalc)[1]])
                        Yoff = np.max(Ydiff[1])-np.min(Yobs[1])
                        Ydiff[1] -= Yoff
                        if ifXray:
                            Labels[label][0] = r'$\mathsf{2\theta ,deg}$'
                            Labels[label][1] = 'Intensity'
                        G2plt.PlotXY(G2frame,[Yobs,Ycalc],XY2=[Ydiff,],labelX=Labels[label][0],
                            labelY=Labels[label][1],newPlot=True,Title=Labels[label][2]+pName,
                            lines=True,names=Names[1:])
                    else:
                        G2plt.PlotXY(G2frame,[Ycalc,Yobs],labelX=Labels[label][0],
                            labelY=Labels[label][1],newPlot=True,Title=Labels[label][2]+pName,
                            lines=True,names=Names[1:])
                        RMCPdict[pName+label] = np.sum(Ycalc[1])/np.sum(Yobs[1])
                        print(' %s scale Ycalc/Yobs: %.4f'%(label,RMCPdict[pName+label]))
#partials plots
            Labels = {'_PDFpartials.csv':[r'$\mathsf{R,\AA}$','G(R)','RMCP G(R) partials for '],
                '_SQ1partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)','RMCP S(Q) partials for '],
                '_SQ2partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','S(Q)-2','RMCP S(Q) partials for '],
                '_FQ1partials.csv':[r'$\mathsf{Q,\AA^{-1}}$','F(Q)','RMCP F(Q) partials for ']}
            for label in Labels:
                X = []
                Partials = []
                if len(files[label]):
                    Names = files[label][0][:-1].split(',')
                    for line in files[label][1:]:
                        items = line.split(',')[:-1]
                        X.append(float(items[0]))
                        Partials.append([float(item) for item in items[1:]])
                    X = np.array(X)
                    DX = X[1]-X[0]  #should be 0.02
                    Partials = np.array(Partials).T
                    if 'Q' in label:
                        continue            #skip these partials
                        XY = [[X.T,Y.T] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
                    else:
                        if ifNeut:
                            XY = [[X.T,(DX*Y.T)] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
                        else:
                            XY = [[X.T,(DX*Y.T)*X.T] for iy,Y in enumerate(Partials) if 'Va' not in Names[iy+1]]
                    Names = [name for name in Names if 'Va' not in name]
                    ylabel = Labels[label][1]
                    if 'G(R)' in Labels[label][1]:
                        if ifNeut:
                            title = 'Neutron '+Labels[label][2]+pName
                        else:
                            continue        #skip for now - x-ray partials are missing header record
                            title = 'X-ray '+Labels[label][2].replace('G','g')+pName
                            ylabel = 'g(R)'
                        sumAtm = 0
                        BLtables = G2elem.GetBLtable(generalData)
                        AtNum = generalData['NoAtoms']
                        if ifNeut:
                            for atm in AtNum: sumAtm += AtNum[atm]
                        else:
                            for atm in AtNum:
                                sumAtm += AtNum[atm]*G2elem.GetFFtable([atm,])[atm]['Z']
                        bfac = {}
                        bcorr = []
                        for atm in AtNum:
                            if ifNeut:
                                if 'SL' in BLtables[atm][1]:
                                    bfac[atm] = 10.*BLtables[atm][1]['SL'][0]*AtNum[atm]/sumAtm  #scale to pm
                                else:   #resonant scatters (unlikely!)
                                    bfac[atm] = AtNum[atm]/sumAtm
                            else:
                                bfac[atm] = 10.*G2elem.GetFFtable([atm,])[atm]['Z']*AtNum[atm]/sumAtm
                        for name in Names:
                            if '-' in name:
                                at1,at2 = name.strip().split('-')
                                if 'Va' in name:
                                    bcorr.append(0.)
                                else:
                                    bcorr.append(bfac[at1]*bfac[at2])
                                if at1 == at2:
                                    bcorr[-1] /= 2.         #no double counting
                        for ixy,xy in enumerate(XY):
                            xy[1] *= bcorr[ixy]
                            xy[1] += Ymin
                        Xmax = np.searchsorted(Ysave[0][0],XY[0][0][-1])
                        G2plt.PlotXY(G2frame,XY2=XY,XY=[Ysave[0][:,0:Xmax],],labelX=Labels[label][0],
                            labelY=ylabel,newPlot=True,Title=title,
                            lines=False,names=[r'   $G(R)_{calc}$',]+Names[1:])
                    else:
                        G2plt.PlotXY(G2frame,XY,labelX=Labels[label][0],
                            labelY=ylabel,newPlot=True,Title=Labels[label][2]+pName,
                            lines=True,names=Names[1:])
#chi**2 plot
            X = []
            Chi = []
            if len(files['.chi2']) > 2:
                Names = files['.chi2'][0][:-1].split()
                for line in files['.chi2'][1:]:
                    items = line[:-1].split()
                    X.append(float(items[1]))
                    Chi.append([float(item) for item in items[3:]])
                X = np.array(X)
                Chi = np.array(Chi).T
                XY = [[X.T,np.log10(Y.T)] for Y in Chi]
                G2plt.PlotXY(G2frame,XY,labelX='no. generated',
                    labelY=r'$log_{10}$ (reduced $\mathsf{\chi^2})$',newPlot=True,Title='RMCP Chi^2 for '+pName,
                    lines=True,names=Names[3:])

#get atoms from rmc6f file
            rmc6fName = pName+'.rmc6f'
            rmc6f = open(rmc6fName,'r')
            rmc6fAtoms = []
            while True:
                line = rmc6f.readline()
                if 'Number of atoms:' in line:
                    Natm = int(line.split(':')[1])
                if 'Atoms:' in line:
                    break
            for iAtm in range(Natm):
                line = rmc6f.readline().split()
                rmc6fAtoms.append([line[1],float(line[3]),float(line[4]),float(line[5])])
            rmc6f.close()
#alt bond histograms - from rmc6 & bond files

            bondName = pName+'.bonds'
            if os.path.exists(os.path.join(path,bondName)):
                nBond = len(RMCPdict['Potentials']['Stretch'])
                bondList = []
                bonds = open(bondName,'r')
                while True:
                    line = bonds.readline()
                    if '............' in line:
                        break       #positions at start of 1st bond list
                for iBond in range(nBond):
                    bondList.append([])
                    Bonds = RMCPdict['Potentials']['Stretch'][iBond]
                    for iAtm in range(Natm):     #same as in rmc6f above
                        line = bonds.readline().split('::')
                        if Bonds[0] in line[0]:
                            items = line[1].replace(';','').split()[:-1:2]
                            items = np.array([int(item) for item in items])
                            bondList[iBond].append([iAtm,items])
                bonds.close()
                for iBond in range(nBond):
                    Bonds = RMCPdict['Potentials']['Stretch'][iBond]
                    title = '%s-%s'%(Bonds[0],Bonds[1])
                    bondDist = G2pwd.GetRMCBonds(generalData,RMCPdict,rmc6fAtoms,bondList[iBond])
                    print('%s mean %.3f(%d)'%(title,np.mean(bondDist),int(1000*np.std(bondDist))))
                    G2plt.PlotBarGraph(G2frame,bondDist,Xname=r'$Bond, \AA$',Title=title+' from Potential Energy Restraint',
                        PlotName='%s Bond for %s'%(title,pName))
                    print(' %d %s bonds found'%(len(bondDist),title))

#alt angle histograms - from rmc6 & triplets files
            tripName = pName+'.triplets'
            if os.path.exists(os.path.join(path,tripName)):
                nAng = len(RMCPdict['Potentials']['Angles'])
                tripList = []
                triples = open(tripName,'r')
                while True:
                    line = triples.readline()
                    if '............' in line:
                        break       #positions at start of 1st triple list
                for iAng in range(nAng):
                    tripList.append([])
                    Angles = RMCPdict['Potentials']['Angles'][iAng]
                    for iAtm in range(Natm):     #same as in rmc6f above
                        line = triples.readline().split('::')
                        if Angles[1] in line[0]:
                            items = line[1].replace(';','').split()[:-1:2]
                            items = np.array([int(item) for item in items]).reshape((-1,2))
                            tripList[iAng].append([iAtm,items])
                        line = triples.readline()       #to skip a line
                triples.close()
                for iAng in range(nAng):
                    Angles = RMCPdict['Potentials']['Angles'][iAng]
                    angles = G2pwd.GetRMCAngles(generalData,RMCPdict,rmc6fAtoms,tripList[iAng])
                    title = '%s-%s-%s'%(Angles[0],Angles[1],Angles[2])
                    print('%s mean %.2f(%d)'%(title,np.mean(angles),int(100*np.std(angles))))
                    G2plt.PlotBarGraph(G2frame,angles,Xname=r'$Angle, \AA$',Title=title+' from Potential Energy Restraint',
                        PlotName='%s Angle for %s'%(title,pName))
                    print(' %d %s angles found'%(len(angles),title))

#bond odf plots
            nPot = len(RMCPdict['Potentials']['Stretch'])
            for iPot in range(nPot):
                fname = pName+'.bondodf_%d'%(iPot+1)
                bond = RMCPdict['Potentials']['Stretch'][iPot]
                if os.path.exists(os.path.join(path,fname)):
                    OutFile = open(fname,'r')
                    odfFile = OutFile.readlines()
                    if len(odfFile) > 1:
                        OutFile.seek(0)
                        odfData = np.fromfile(OutFile,sep=' ')
                        numx,numy = odfData[:2]
                        G2plt.Plot3dXYZ(G2frame,int(numx),int(numy),odfData[2:],
                            newPlot=False,Title='Number of %s-%s Bonds'%(bond[0],bond[1]),Centro=True)
                    OutFile.close()
        elif G2frame.RMCchoice == 'PDFfit':
            generalData = data['General']
            RMCPdict = data['RMC']['PDFfit']
            pName = generalData['Name'].replace(' ','_')
            Labels = [r'$\mathsf{R,\AA}$','G(R)','PDFfit2 G(R) for ']
            files = RMCPdict['files']
            for file in files:
                if 'Select' not in files[file][0]:
                    XY = np.empty((1,2))
                    start = 0
                    while XY.shape[0] == 1:
                        try:
                            XY = np.loadtxt(files[file][0],skiprows=start)
                        except ValueError:
                            start += 1
                    if 'Neutron' in file:
                        cname = pName+'-PDFfitN.fgr'
                    else:
                        cname = pName+'-PDFfitX.fgr'
                    XYobs = XY.T[:2]
                    XY = np.empty((1,2))
                    start = 0
                    while XY.shape[0] == 1:
                        try:
                            XY = np.loadtxt(cname,skiprows=start)
                        except ValueError:
                            start += 1
                    XYcalc = XY.T[:2]
                    XYdiff = np.zeros_like(XYcalc)
                    XYdiff[0] = XYcalc[0]
                    ibeg = np.searchsorted(XYobs[0],XYdiff[0][0])
                    ifin = ibeg+XYcalc.shape[1]
                    XYdiff[1] = XYobs[1][ibeg:ifin]-XYcalc[1]
                    offset = 1.1*(np.max(XYdiff[1])-np.min(XYcalc[1]))
                    XYdiff[1] -= offset
                    G2plt.PlotXY(G2frame,[XYobs,],XY2=[XYcalc,XYdiff],labelX=Labels[0],
                        labelY=Labels[1],newPlot=True,Title=Labels[2]+files[file][0],
                        lines=False,names=['G(R) obs','G(R) calc','diff',])

#### ISODISTORT tab ###############################################################################

    def OnRunISODISTORT(event):
        ''' this needs to setup for method #3 or #4 in ISODISTORT
        after providing parent cif:
        #3 asks for transformation matrix & space group of child structure
        #4 asks for cif file of child structure
        '''

        if not G2frame.GSASprojectfile:     #force a project save to establish location of output cif file
            G2frame.OnFileSaveas(event)

        radio,rundata = ISO.GetISODISTORT(data)
        if radio:
            data['ISODISTORT']['radio'] = radio
            data['ISODISTORT']['rundata'] = rundata
            data['ISODISTORT']['SGselect'] =  {'Tric':True,'Mono':True,'Orth':True,'Tetr':True,'Trig':True,'Hexa':True,'Cubi':True}
            data['ISODISTORT']['selection'] = None
            print('ISODISTORT run complete')
            wx.CallAfter(G2phsG2.UpdateISODISTORT,G2frame.data)
        elif data['ISODISTORT']['ISOmethod'] != 4 or radio is None:
            G2G.G2MessageBox(G2frame,'ISODISTORT run failed - see page opened in web browser')
        else:
            G2G.G2MessageBox(G2frame,'ISODISTORT run complete; new cif file %s created.\n To use, import it as a new phase.'%rundata)
            print(' ISODISTORT run complete; new cif file %s created. To use, import it as a new phase.'%rundata)

    def OnNewISOPhase(event):
        ''' Make CIF file with ISODISTORT
        '''
        if 'rundata' in data['ISODISTORT'] and data['ISODISTORT']['selection'] is not None:
            CIFfile = ISO.GetISODISTORTcif(data)
            G2G.G2MessageBox(G2frame,'ISODISTORT generated cif file %s has been created.'%CIFfile)
        elif 'rundata' in data['ISODISTORT']:
            G2G.G2MessageBox(G2frame,'Need to select an ISODISTORTdistortion model first before creating a CIF')
        else:
            G2G.G2MessageBox(G2frame,'ERROR - need to run ISODISTORT first - see General/Compute menu')

    def OnNewPDFfitPhase(event):
        ''' Make new phase for PDFfit using ISODISTORT mode definitions as constraints
        '''
        newPhase = G2pwd.ISO2PDFfit(data)
        phaseName = newPhase['General']['Name']
        sub = G2frame.GPXtree.AppendItem(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases'),text=phaseName)
        G2frame.GPXtree.SetItemPyData(sub,newPhase)
        G2frame.GPXtree.SelectItem(sub)

#### DIFFax Layer Data page ################################################################################
# Utility routines her
    def OnCopyPhase(event):
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
        dlg = wx.SingleChoiceDialog(G2frame,'Phase to use for cell data','Select',phaseNames)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            PhaseName = phaseNames[sel]
        else:
            return
        General = G2stIO.GetAllPhaseData(GPXFile,PhaseName)['General']
        data['Layers']['Cell'] = General['Cell']
        wx.CallAfter(G2phsG2.UpdateLayerData,G2frame,data)

    def OnLoadDIFFaX(event):
        if len(data['Layers']['Layers']):
            dlg = wx.MessageDialog(G2frame,'Do you really want to replace the Layer data?','Load DIFFaX file',
                wx.YES_NO | wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_NO:
                    return
            finally:
                dlg.Destroy()
        dlg = wx.FileDialog(G2frame, 'Choose DIFFaX file name to read', G2G.GetImportPath(G2frame), '',
            'DIFFaX file (*.*)|*.*',style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                DIFFaXfile = dlg.GetPath()
                data['Layers'] = G2IO.ReadDIFFaX(DIFFaXfile)
        finally:
            dlg.Destroy()
        wx.CallAfter(G2phsG2.UpdateLayerData,G2frame,data)

    def OnSimulate(event):
        debug = False       #set True to run DIFFax to compare/debug (must be in bin)
        idebug = 0
        if debug: idebug = 1
        wx.MessageBox(' For use of DIFFaX, please cite:\n\n'+
                          G2G.GetCite('DIFFaX'),
                          caption='DIFFaX',style=wx.ICON_INFORMATION)
        ctrls = ''
        dlg = DIFFaXcontrols(G2frame,ctrls)
        if dlg.ShowModal() == wx.ID_OK:
            simCodes = dlg.GetSelection()
        else:
            return

        if 'PWDR' in  simCodes[0]:    #powder pattern
            data['Layers']['selInst'] = simCodes[1]
            UseList = [item for item in data['Histograms'] if 'PWDR' in item]
            if not UseList:
                wx.MessageBox('No PWDR data for this phase to simulate',caption='Data error',style=wx.ICON_EXCLAMATION)
                return
            elif len(UseList) == 1: # don't ask questions when we know the answer!
                HistName = UseList[0]
            else:
                dlg = wx.SingleChoiceDialog(G2frame,'Data to simulate','Select',UseList)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    HistName = UseList[sel]
                else:
                    return
                dlg.Destroy()
            G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistName)
            sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
                G2frame,G2frame.PatternId, 'Sample Parameters'))
            scale = sample['Scale'][0]
            background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
                G2frame,G2frame.PatternId, 'Background'))
            limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
                G2frame,G2frame.PatternId, 'Limits'))[1]
            inst = G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
            if 'T' in inst['Type'][0]:
                wx.MessageBox("Can't simulate neutron TOF patterns yet",caption='Data error',style=wx.ICON_EXCLAMATION)
                return
            profile = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1]
            G2pwd.CalcStackingPWDR(data['Layers'],scale,background,limits,inst,profile,debug)
            if debug:
                ctrls = '0\n%d\n3\n'%(idebug)
                G2pwd.StackSim(data['Layers'],ctrls,scale,background,limits,inst,profile)
                test1 = np.copy(profile[3])
                test1 = np.where(test1,test1,1.0)
                G2pwd.CalcStackingPWDR(data['Layers'],scale,background,limits,inst,profile,debug)
                test2 = np.copy(profile[3])
                rat = test1-test2
                XY = np.vstack((profile[0],rat))
                G2plt.PlotXY(G2frame,[XY,],XY2=[],labelX=r'$\mathsf{2\theta}$',
                    labelY='difference',newPlot=True,Title='DIFFaX vs GSASII',lines=True)
            G2pwpl.PlotPatterns(G2frame,plotType='PWDR',newPlot=True)
        else:   #selected area
            data['Layers']['Sadp'] = {}
            data['Layers']['Sadp']['Plane'] = simCodes[1]
            data['Layers']['Sadp']['Lmax'] = simCodes[2]
            if debug:
                planeChoice = ['h0l','0kl','hhl','h-hl',]
                lmaxChoice = [str(i+1) for i in range(6)]
                ctrls = '0\n%d\n4\n1\n%d\n%d\n16\n1\n1\n0\nend\n'%    \
                    (idebug,planeChoice.index(simCodes[1])+1,lmaxChoice.index(simCodes[2])+1)
                G2pwd.StackSim(data['Layers'],ctrls)
            G2pwd.CalcStackingSADP(data['Layers'],debug)
        wx.MessageBox('Simulation finished',caption='Stacking fault simulation',style=wx.ICON_EXCLAMATION)
        wx.CallAfter(G2phsG2.UpdateLayerData,G2frame,data)

    def OnFitLayers(event):
        print (' fit stacking fault model TBD')
#        import scipy.optimize as opt
        wx.BeginBusyCursor()
        # see pwd.SetupPDFEval() and pwd.OptimizePDF() for an example minimization
        wx.EndBusyCursor()
        wx.CallAfter(G2phsG2.UpdateLayerData,G2frame,data)
        G2pwpl.PlotPatterns(G2frame,plotType='PWDR')

    def OnSeqSimulate(event):

        cellSel = ['cellA','cellB','cellC','cellG']
        transSel = ['TransP','TransX','TransY','TransZ']
        ctrls = ''
        cell = data['Layers']['Cell']
        data['Layers']['seqResults'] = []
        data['Layers']['seqCodes'] = []
        Parms = G2pwd.GetStackParms(data['Layers'])
        dlg = DIFFaXcontrols(G2frame,ctrls,Parms)
        if dlg.ShowModal() == wx.ID_OK:
            simCodes = dlg.GetSelection()
        else:
            return
        UseList = []
        for item in data['Histograms']:
            if 'PWDR' in item:
                UseList.append(item)
        if not UseList:
            wx.MessageBox('No PWDR data for this phase to simulate',caption='Data error',style=wx.ICON_EXCLAMATION)
            return
        dlg = wx.SingleChoiceDialog(G2frame,'Data to simulate','Select',UseList)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            HistName = UseList[sel]
        else:
            return
        dlg.Destroy()
        G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistName)
        sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
            G2frame,G2frame.PatternId, 'Sample Parameters'))
        scale = sample['Scale'][0]
        background = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
            G2frame,G2frame.PatternId, 'Background'))
        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
            G2frame,G2frame.PatternId, 'Limits'))[1]
        inst = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))[0]
        if 'T' in inst['Type'][0]:
            wx.MessageBox("Can't simulate neutron TOF patterns yet",caption='Data error',style=wx.ICON_EXCLAMATION)
            return
        profile = np.copy(G2frame.GPXtree.GetItemPyData(G2frame.PatternId)[1])
        resultXY2 = []
        resultXY = [np.vstack((profile[0],profile[1])),]    #observed data
        data['Layers']['selInst'] = simCodes[1]
        data['Layers']['seqCodes'] = simCodes[2:]
        Layers = copy.deepcopy(data['Layers'])
        pName = simCodes[2]
        BegFin = simCodes[3]
        nSteps = simCodes[4]
        laue = Layers['Laue']
        vals = np.linspace(BegFin[0],BegFin[1],nSteps+1,True)
        simNames = []
        for val in vals:
            print (' Stacking simulation step for %s = %.5f'%(pName,val))
            simNames.append('%.3f'%(val))
            if 'cell' in pName:
                cellId = cellSel.index(pName)
                cell = Layers['Cell']
                cell[cellId+1] = val
                if laue in ['-3','-3m','6/m','6/mmm','4/m','4/mmm']:
                    cell[2] = cell[1]
                cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
                Layers['Cell'] = cell
            elif 'Trans' in pName:
                names = pName.split(';')
                transId = transSel.index(names[0])
                iY = int(names[1])
                iX = int(names[2])
                Trans = Layers['Transitions'][iY]
                Nx = len(Trans)-1
                if not transId:     #i.e. probability
                    osum = 1.-Trans[iX][0]
                    nsum = 1.-val
                    for i in range(Nx+1):
                        if i != iX:
                            Trans[i][0] *= (nsum/osum)
                    Trans[iX][0] = val
                    if Layers.get('SymTrans',False):
                        Layers['Transitions'][Nx-iX][Nx-iY][0] = val
                        for i in range(Nx+1):
                            Layers['Transitions'][Nx-iY][Nx-i][0] = Layers['Transitions'][iY][i][0]
                    print (' Transition matrix:')
                    for trans in Layers['Transitions']:
                        line = str([' %.3f'%(item[0]) for item in trans])
                        print (line[1:-2].replace("'",''))
                else:
                    Trans[iX][transId] = val
            G2pwd.CalcStackingPWDR(Layers,scale,background,limits,inst,profile,False)
            resultXY2.append([np.vstack((profile[0],profile[3])),][0])
        data['Layers']['seqResults'] = [resultXY,resultXY2,simNames]
        wx.MessageBox('Sequential simulation finished',caption='Stacking fault simulation',style=wx.ICON_EXCLAMATION)
        wx.CallAfter(G2phsG2.UpdateLayerData,G2frame,data)

#### Wave Data page ################################################################################

    def OnWaveVary(event):
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        atomData = data['Atoms']
        atNames = []
        names = ['Sfrac','Spos','Sadp','Smag']
        flags = dict(zip(names,[[],[],[],[]]))
        for atom in atomData:
            atNames.append(atom[ct-1])
            waves = atom[-1]['SS1']
            for name in names:
                if waves[name]:
                    flags[name].append(True)
                else:
                    flags[name].append(False)
        dlg = G2G.FlagSetDialog(G2frame,'Wave refinement flags',['Atom',]+names,atNames,flags)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                flags = dlg.GetSelection()
                for ia,atom in enumerate(atomData):
                    for name in names:
                        for wave in atom[-1]['SS1'][name][1:]:
                            wave[1] = flags[name][ia]
        finally:
            dlg.Destroy()
        G2phsG2.UpdateWavesData(G2frame,data)

#### Structure drawing GUI stuff ################################################################################
    def SetupDrawingData():
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = data['Atoms']
        defaultDrawing = {'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,
            'backColor':[0,0,0],'depthFog':False,'Zclip':50.0,'cameraPos':50.,'Zstep':0.5,
            'radiusFactor':0.85,'contourLevel':1.,'bondRadius':0.1,'ballScale':0.33,
            'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,'unitCellBox':True,'contourMax':1.0,
            'showABC':True,'selectedAtoms':[],'Atoms':[],'oldxy':[],'magMult':1.0,'SymFade':False,
            'bondList':{},'viewDir':[1,0,0],'Plane':[[0,0,1],False,False,0.0,[255,255,0]]}
        V0 = np.array([0,0,1])
        V = np.inner(Amat,V0)
        V /= np.sqrt(np.sum(V**2))
        A = np.arccos(np.sum(V*V0))
        defaultDrawing['Quaternion'] = G2mth.AV2Q(A,[0,1,0])
        try:
            drawingData = data['Drawing']
        except KeyError:
            data['Drawing'] = {}
            drawingData = data['Drawing']
        if not drawingData:                 #fill with defaults if empty
            drawingData = defaultDrawing.copy()
        if 'Zstep' not in drawingData:
            drawingData['Zstep'] = 0.5
        if 'contourLevel' not in drawingData:
            drawingData['contourLevel'] = 1.
        if 'contourMax' not in drawingData:
            drawingData['contourMax'] = 1.
        if 'viewDir' not in drawingData:
            drawingData['viewDir'] = [0,0,1]
        if 'Quaternion' not in drawingData:
            drawingData['Quaternion'] = G2mth.AV2Q(2*np.pi,np.inner(Amat,[0,0,1]))
        if 'showRigidBodies' not in drawingData:
            drawingData['showRigidBodies'] = True
        try:  # patch of sorts; this had been set to a string; needs to be an int between 0 & 3
            int(drawingData['showSlice'])
        except:
            drawingData['showSlice'] = 0
        if 'sliceSize' not in drawingData:
            drawingData['sliceSize'] = 5.0
        if 'contourColor' not in drawingData:
            drawingData['contourColor'] = 'RdYlGn'
        if 'Plane' not in drawingData:
            drawingData['Plane'] = [[0,0,1],False,False,0.0,[255,255,0]]
        if 'magMult' not in drawingData:
            drawingData['magMult'] = 1.0
        if 'SymFade' not in drawingData:
            drawingData['SymFade'] = False
        if 'Voids' not in drawingData:
            drawingData['Voids'] = []
            drawingData['showVoids'] = False
            drawingData['showMap'] = False
        cx,ct,cs,ci = [0,0,0,0]
        if generalData['Type'] in ['nuclear','faulted',]:
            cx,ct,cs,ci = [2,1,6,17]         #x, type, style & index
        elif generalData['Type'] == 'macromolecular':
            cx,ct,cs,ci = [5,4,9,20]         #x, type, style & index
        elif generalData['Type'] == 'magnetic':
            cx,ct,cs,ci = [2,1,9,20]         #x, type, style & index
            drawingData['vdwScale'] = 0.20
        drawingData['atomPtrs'] = [cx,ct,cs,ci]
        if not drawingData.get('Atoms'):
            for atom in atomData:
                DrawAtomAdd(drawingData,atom)
            data['Drawing'] = drawingData
        if len(drawingData['Plane']) < 5:
            drawingData['Plane'].append([255,255,0])

    def DrawAtomAdd(drawingData,atom):
        drawingData['Atoms'].append(G2mth.MakeDrawAtom(data,atom))

    def OnRestraint(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        #indx = drawAtoms.GetSelectedRows()
        restData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints'))
        drawingData = data['Drawing']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = drawingData['Atoms']
        atXYZ = []
        atSymOp = []
        atIndx = []
        for item in indx:
            atXYZ.append(np.array(atomData[item][cx:cx+3]))
            atSymOp.append(atomData[item][cs-1])
            atIndx.append(atomData[item][ci])
        if event.GetId() == G2G.wxID_DRAWRESTRBOND and len(indx) == 2:
            try:
                bondData = restData[PhaseName]['Bond']
            except KeyError:
                bondData = {'wtFactor':1.0,'Bonds':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Bond'] = bondData
            bondData['Bonds'].append([atIndx,atSymOp,1.54,0.01])
        elif event.GetId() == G2G.wxID_DRAWRESTRANGLE and len(indx) == 3:
            try:
                angleData = restData[PhaseName]['Angle']
            except KeyError:
                angleData = {'wtFactor':1.0,'Angles':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Angle'] = angleData
            angleData['Angles'].append([atIndx,atSymOp,109.5,1.0])
        elif event.GetId() == G2G.wxID_DRAWRESTRPLANE and len(indx) > 3:
            try:
                planeData = restData[PhaseName]['Plane']
            except KeyError:
                planeData = {'wtFactor':1.0,'Planes':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Plane'] = planeData
            planeData['Planes'].append([atIndx,atSymOp,0.0,0.01])
        elif event.GetId() == G2G.wxID_DRAWRESTRCHIRAL and len(indx) == 4:
            try:
                chiralData = restData[PhaseName]['Chiral']
            except KeyError:
                chiralData = {'wtFactor':1.0,'Volumes':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Chiral'] = chiralData
            chiralData['Volumes'].append([atIndx,atSymOp,2.5,0.1])
        else:
            print ('**** ERROR wrong number of atoms selected for this restraint')
            return
        G2frame.GPXtree.SetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints'),restData)

    def OnDefineRB(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        indx.sort()
        RBData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        drawingData = data['Drawing']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = drawingData['Atoms']
        rbXYZ = []
        rbType = []
        atNames = []
        AtInfo = RBData['Residue']['AtInfo']
        for i,item in enumerate(indx):
            rbtype = atomData[item][ct]
            atNames.append(rbtype+str(i))
            rbType.append(rbtype)
            if rbtype not in AtInfo:
                Info = G2elem.GetAtomInfo(rbtype)
                AtInfo[rbtype] = [Info['Drad'],Info['Color']]
            rbXYZ.append(np.inner(np.array(atomData[item][cx:cx+3]),Amat))
        rbXYZ = np.array(rbXYZ)
        rbXYZ -= rbXYZ[0]
        rbId = ran.randint(0,sys.maxsize)
        rbName = 'UNKRB'
        dlg = wx.TextEntryDialog(G2frame,'Enter the name for the new rigid body',
            'Edit rigid body name',rbName ,style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            rbName = dlg.GetValue()
        dlg.Destroy()
        RBData['Residue'][rbId] = {'RBname':rbName,'rbXYZ':rbXYZ,'rbTypes':rbType,
            'atNames':atNames,'rbRef':[0,1,2,False],'rbSeq':[],'SelSeq':[0,0],'useCount':0}
        RBData['RBIds']['Residue'].append(rbId)
        G2frame.GetStatusBar().SetStatusText('New rigid body UNKRB added to set of Residue rigid bodies',1)

#### Draw Atom routines ################################################################################
    def UpdateDrawAtoms(G2frame,data,atomStyle=''):
        drawAtoms = G2frame.drawAtoms
        def RefreshDrawAtomGrid(event):
            def SetChoice(name,c,n=0):
                choice = []
                for r in range(len(atomData)):
                    if n:
                        srchStr = str(atomData[r][c][:n])
                    else:
                        srchStr = str(atomData[r][c])
                    if srchStr not in choice:
                        if n:
                            choice.append(str(atomData[r][c][:n]))
                        else:
                            choice.append(str(atomData[r][c]))
                choice.sort()
    
                dlg = wx.MultiChoiceDialog(G2frame,'Select',name,choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelections()
                    parms = []
                    for x in sel:
                        parms.append(choice[x])
                    drawAtoms.ClearSelection()
                    drawingData['selectedAtoms'] = []
                    for row in range(len(atomData)):
                        test = atomData[row][c]
                        if n:
                            test = test[:n]
                        if  test in parms:
                            drawAtoms.SelectRow(row,True)
                            drawingData['selectedAtoms'].append(row)
                    G2plt.PlotStructure(G2frame,data)
                dlg.Destroy()
    
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                for row in range(drawAtoms.GetNumberRows()):
                    drawingData['selectedAtoms'].append(row)
                    drawAtoms.SelectRow(row,True)
            elif r < 0:                          #dclick on col label
                # sel = -1
                if drawAtoms.GetColLabelValue(c) == 'Style':
                    DrawAtomStyle(event)
                elif drawAtoms.GetColLabelValue(c) == 'Label':
                    DrawAtomLabel(event)
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    DrawAtomColor(event)
                elif drawAtoms.GetColLabelValue(c) == 'Residue':
                    SetChoice('Residue',c,3)
                elif drawAtoms.GetColLabelValue(c) == '1-letter':
                    SetChoice('1-letter',c,1)
                elif drawAtoms.GetColLabelValue(c) == 'Chain':
                    SetChoice('Chain',c)
                elif drawAtoms.GetColLabelValue(c) == 'Name':
                    SetChoice('Name',c)
                elif drawAtoms.GetColLabelValue(c) == 'Sym Op':
                    SetChoice('Name',c)
                elif drawAtoms.GetColLabelValue(c) == 'Type':
                    SetChoice('Type',c)
                elif drawAtoms.GetColLabelValue(c) in ['x','y','z','I/A']:
                    drawAtoms.ClearSelection()
            else:
                if drawAtoms.GetColLabelValue(c) in ['Style','Label']:
                    atomData[r][c] = drawAtoms.GetCellValue(r,c)
                    FindBondsDraw(data)
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    colors = wx.ColourData()
                    colors.SetChooseFull(True)
                    dlg = wx.ColourDialog(G2frame.GetParent(),colors)
                    if dlg.ShowModal() == wx.ID_OK:
                        color = dlg.GetColourData().GetColour()[:3]
                        attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                        attr.SetReadOnly(True)
                        attr.SetBackgroundColour(color)
                        atomData[r][c] = color
                        drawingData['Atoms'][r][c] = color
                        drawAtoms.SetAttr(i,cs+2,attr)
                    dlg.Destroy()
                    UpdateDrawAtoms(G2frame,data)
            G2plt.PlotStructure(G2frame,data)

        def NextAtom(event):
            'respond to a tab by cycling through the atoms'
            next = 0
            for r in drawAtoms.GetSelectedRows():
                next = r + 1
                break
            if next >= drawAtoms.GetNumberRows():
                next = 0
            drawAtoms.ClearSelection()
            drawAtoms.SelectRow(next,True)
            drawAtoms.MakeCellVisible(next,0)
            drawingData['selectedAtoms'] = drawAtoms.GetSelectedRows()
            G2plt.PlotStructure(G2frame,data)
            G2frame.Raise()
    
        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                if drawAtoms.IsSelection():
                    drawAtoms.ClearSelection()
            elif c < 0:                   #only row clicks
                if event.ControlDown():
                    if r in drawAtoms.GetSelectedRows():
                        drawAtoms.DeselectRow(r)
                    else:
                        drawAtoms.SelectRow(r,True)
                elif event.ShiftDown():
                    indxs = drawAtoms.GetSelectedRows()
                    drawAtoms.ClearSelection()
                    ibeg = 0
                    if indxs:
                        ibeg = indxs[-1]
                    for row in range(ibeg,r+1):
                        drawAtoms.SelectRow(row,True)
                else:
                    G2frame.GetStatusBar().SetStatusText('Use right mouse click to brng up Draw Atom editing options',1)
                    drawAtoms.ClearSelection()
                    drawAtoms.SelectRow(r,True)
            drawingData['selectedAtoms'] = []
            drawingData['selectedAtoms'] = drawAtoms.GetSelectedRows()
            G2plt.PlotStructure(G2frame,data)
    
    #### UpdateDrawAtoms executable code starts here
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Draw Atom list for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)
        G2frame.GetStatusBar().SetStatusText('',1)
        oldSizer = G2frame.drawAtomsList.GetSizer()
        if oldSizer: # 2nd+ use, clear out old entries
            for i in oldSizer.GetChildren(): # look for grids in sizer
                if type(i.GetWindow()) is G2G.GSGrid:
                    oldSizer.Detach(i.GetWindow())  # don't delete them
            oldSizer.Clear(True)
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
        SetDrawingDefaults(drawingData)
        cx,ct,cs,ci = drawingData['atomPtrs']
        atomData = drawingData['Atoms']
        if atomStyle:
            for atom in atomData:
                atom[cs] = atomStyle
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,5',]+ \
            [wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,lines,vdW balls,sticks,balls & sticks,ellipsoids,polyhedra",
            wg.GRID_VALUE_CHOICE+": ,type,name,number",wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,]
        # styleChoice = DrawStyleChoice
        # labelChoice = [' ','type','name','number']
        colLabels = ['Name','Type','x','y','z','Sym Op','Style','Label','Color','I/A']
        if generalData['Type'] == 'macromolecular':
            colLabels = ['Residue','1-letter','Chain'] + colLabels
            Types = 3*[wg.GRID_VALUE_STRING,]+Types
            Types[8] = wg.GRID_VALUE_CHOICE+": ,lines,vdW balls,sticks,balls & sticks,ellipsoids,backbone,ribbons,schematic"
            # styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','backbone','ribbons','schematic']
            # labelChoice = [' ','type','name','number','residue','1-letter','chain']
            Types[9] = wg.GRID_VALUE_CHOICE+": ,type,name,number,residue,1-letter,chain"
        elif generalData['Type'] == 'magnetic':
            colLabels = colLabels[:5]+['Mx','My','Mz']+colLabels[5:]
            Types = Types[:5]+3*[wg.GRID_VALUE_FLOAT+':10,4',]+Types[5:]
        table = []
        rowLabels = []
        for i,atom in enumerate(drawingData['Atoms']):
            table.append(atom[:colLabels.index('I/A')+1])
            rowLabels.append(str(i))
    
        G2frame.atomTable = None
        G2frame.atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        drawAtoms.SetTable(G2frame.atomTable, True)
        drawAtoms.SetMargins(0,0)
        drawAtoms.AutoSizeColumns(True)
        drawAtoms.SetColSize(colLabels.index('Style'),80)
        drawAtoms.SetColSize(colLabels.index('Color'),50)
        drawAtoms.SetRowLabelSize(45)
        if 'phoenix' in wx.version():
            drawAtoms.Unbind(wg.EVT_GRID_CELL_CHANGED)
            drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGED, RefreshDrawAtomGrid)
        else:
            drawAtoms.Unbind(wg.EVT_GRID_CELL_CHANGE)
            drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshDrawAtomGrid)
        drawAtoms.Unbind(wg.EVT_GRID_LABEL_LEFT_DCLICK)
        drawAtoms.Unbind(wg.EVT_GRID_CELL_LEFT_DCLICK)
        drawAtoms.Unbind(wg.EVT_GRID_LABEL_LEFT_CLICK)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshDrawAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, RefreshDrawAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
    
        lblList = ('Delete','Set atom style','Set atom label',
                           'Set atom color','Set view point','Generate copy',
                           'Generate surrounding sphere','Transform atoms',
                           'Generate bonded','Select from list')
        callList = (DrawAtomsDelete,DrawAtomStyle, DrawAtomLabel,
                            DrawAtomColor,SetViewPoint,AddSymEquiv,
                            AddSphere,AddBox,TransformSymEquiv,
                            FillCoordSphere,SelDrawList)
        onRightClick = drawAtoms.setupPopup(lblList,callList)
        drawAtoms.Bind(wg.EVT_GRID_CELL_RIGHT_CLICK, onRightClick)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, onRightClick)
    
        try:
            drawAtoms.Bind(wg.EVT_GRID_TABBING, NextAtom)
        except: # patch: for pre-2.9.5 wx
            pass
        for i,atom in enumerate(drawingData['Atoms']):
            attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
            attr.SetReadOnly(True)
            attr.SetBackgroundColour(atom[cs+2])
            drawAtoms.SetAttr(i,cs+2,attr)
            drawAtoms.SetCellValue(i,cs+2,'')
        indx = drawingData['selectedAtoms']
        if indx:
            for r in range(len(atomData)):
                if r in indx:
                    drawAtoms.SelectRow(r)
        for c in range(len(colLabels)):
           attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
           attr.SetReadOnly(True)
           attr.SetBackgroundColour(VERY_LIGHT_GREY)
           if colLabels[c] not in ['Style','Label','Color']:
                drawAtoms.SetColAttr(c,attr)
    
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(drawAtoms,1,wx.EXPAND)
        drawAtoms.SetScrollRate(10,10) # allow grid to scroll
        SetPhaseWindow(G2frame.drawAtomsList,mainSizer)
    
        FindBondsDraw(data)
        drawAtoms.ClearSelection()
        G2frame.drawAtoms = drawAtoms
    
    def DrawAtomStyle(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        styleChoice = DrawStyleChoice
        if generalData['Type'] == 'macromolecular':
            styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids',
            'backbone','ribbons','schematic']
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom drawing style',styleChoice)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            parms = styleChoice[sel]
            for r in indx:
                atomData[r][cs] = parms
                G2frame.drawAtoms.SetCellValue(r,cs,parms)
        dlg.Destroy()
        FindBondsDraw(data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def DrawAtomLabel(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        styleChoice = [' ','type','name','number']
        if generalData['Type'] == 'macromolecular':
            styleChoice = [' ','type','name','number','residue','1-letter','chain']
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom label style',styleChoice)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            parms = styleChoice[sel]
            for r in indx:
                atomData[r][cs+1] = parms
                G2frame.drawAtoms.SetCellValue(r,cs+1,parms)
        dlg.Destroy()
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def DrawAtomColor(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        if len(indx) > 1:
            G2frame.GetStatusBar().SetStatusText('Select Custom Color, change color, Add to Custom Colors, then OK',1)
        else:
            G2frame.GetStatusBar().SetStatusText('Change color, Add to Custom Colors, then OK',1)
        atomData = data['Drawing']['Atoms']
        atmColors = []
        atmTypes = []
        for r in indx:
            if atomData[r][cs+2] not in atmColors:
                atmColors.append(atomData[r][cs+2])
                atmTypes.append(atomData[r][ct])
                if len(atmColors) > 16:
                    break
        colors = wx.ColourData()
        colors.SetChooseFull(True)
        dlg = wx.ColourDialog(None,colors)
        if dlg.ShowModal() == wx.ID_OK:
            for i in range(len(atmColors)):
                atmColors[i] = dlg.GetColourData().GetColour()[:3]
            colorDict = dict(zip(atmTypes,atmColors))
            for r in indx:
                color = colorDict[atomData[r][ct]]
                atomData[r][cs+2] = color
                attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                attr.SetBackgroundColour(color)
                G2frame.drawAtoms.SetAttr(r,cs+2,attr)
                data['Drawing']['Atoms'][r][cs+2] = color
        dlg.Destroy()
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def ResetAtomColors(event):
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        for atom in atomData:
            atNum = generalData['AtomTypes'].index(atom[ct])
            atom[cs+2] = list(generalData['Color'][atNum])
        UpdateDrawAtoms(G2frame,data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def OnEditAtomRadii(event):
        DisAglCtls = {}
        generalData = data['General']
        if 'DisAglCtls' in generalData:
            DisAglCtls = generalData['DisAglCtls']
        dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData,Angle=False)
        if dlg.ShowModal() == wx.ID_OK:
            DisAglCtls = dlg.GetData()
        dlg.Destroy()
        generalData['DisAglCtls'] = DisAglCtls
        FindBondsDraw(data)
        G2plt.PlotStructure(G2frame,data)

    def SetViewPoint(event):
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(G2frame.drawAtoms,ct-1)
        if not indx: return
        atomData = data['Drawing']['Atoms']
        pos = np.zeros(3)
        for i in indx:
            pos += atomData[i][cx:cx+3]
        data['Drawing']['viewPoint'] = [list(pos/len(indx)),[indx[0],0]]
        G2plt.PlotStructure(G2frame,data)

    def noDuplicate(xyz,atomData):                  #be careful where this is used - it's slow
        cx = data['Drawing']['atomPtrs'][0]
        if True in [np.allclose(np.array(xyz),np.array(atom[cx:cx+3]),atol=0.0002) for atom in atomData]:
            return False
        else:
            return True

    def AddSymEquiv(event):
        drawAtoms = G2frame.drawAtoms
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        indx = getAtomSelections(drawAtoms,ct-1)
        if not indx: return
        indx.sort()
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cuij = ci+2
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        atomData = data['Drawing']['Atoms']
        generalData = data['General']
        SGData = generalData['SGData']
        SpnFlp = SGData.get('SpnFlp',[])
        dlg = SymOpDialog(G2frame,SGData,False,True)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Inv,Cent,Opr,Cell,New,Force = dlg.GetSelection()
                Cell = np.array(Cell)
                cent = SGData['SGCen'][Cent]
                M,T = SGData['SGOps'][Opr]
                for ind in indx:
                    XYZ = np.array(atomData[ind][cx:cx+3])
                    XYZ = np.inner(M,XYZ)+T
                    if Inv and not SGData['SGFixed']:
                        XYZ = -XYZ
                    XYZ = XYZ+cent+Cell
                    if Force:
                        XYZ %= 1.       #G2spc.MoveToUnitCell(XYZ)
                    if noDuplicate(XYZ,atomData):
                        atom = copy.copy(atomData[ind])
                        atom[cx:cx+3] = XYZ
                        atomOp = atom[cs-1]
                        OprNum = ((Opr+1)+100*Cent)*(1-2*Inv)
                        newOp = str(OprNum)+'+'+ \
                            str(int(Cell[0]))+','+str(int(Cell[1]))+','+str(int(Cell[2]))
                        atom[cs-1] = G2spc.StringOpsProd(atomOp,newOp,SGData)
                        if cmx:         #magnetic moment
                            opNum = G2spc.GetOpNum(OprNum,SGData)
                            mom = np.array(atom[cmx:cmx+3])
                            if SGData['SGGray']:
                                atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)
                            else:
                                atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SpnFlp[opNum-1]
                        if atom[cuij] == 'A':
                            Uij = atom[cuij:cuij+6]
                            Uij = G2spc.U2Uij(np.inner(np.inner(M,G2spc.Uij2U(Uij)),M))
                            atom[cuij:cuij+6] = Uij
                        atomData.append(atom[:cuij+9])  #not SS stuff
        finally:
            dlg.Destroy()
        UpdateDrawAtoms(G2frame,data)
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
        
    def AddBox(event):
        print('Box of enclosure, TBD')
        

    def AddSphere(event=None,selection=None,radius=None,targets=None):
        drawAtoms = G2frame.drawAtoms
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        if selection:
            indx = selection
        else:
            indx = getAtomSelections(drawAtoms,ct-1,
                'as center of sphere addition',includeView=True)
        if not indx: return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = data['Drawing']['Atoms']
        numAtoms = len(atomData)
        cuij = cs+5
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        SGData = generalData['SGData']
        SpnFlp = SGData.get('SpnFlp',[])
        cellArray = G2lat.CellBlock(1)
        indx.sort()
        if radius is None or targets is None:
            dlg = SphereEnclosure(G2frame,data['General'],data['Drawing'],indx)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    centers,radius,targets = dlg.GetSelection()
                else:
                    return
            finally:
                dlg.Destroy()
        else:
            centers = []
            for Id in indx:
                if Id < len(data['Drawing']['Atoms']):
                    atom = data['Drawing']['Atoms'][Id]
                    centers.append(atom[cx:cx+3])

        ncent = len(centers)
        pgbar = wx.ProgressDialog('Sphere of enclosure for %d atoms'%ncent,'Centers done=',ncent+1,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        for ic,orig in enumerate(centers):
            xyzA = np.array(orig)
            for atomB in atomData[:numAtoms]:
                if atomB[ct] not in targets:
                    continue
                xyzB = np.array(atomB[cx:cx+3])
                Uij = atomB[cuij:cuij+6]
                result = G2spc.GenAtom(xyzB,SGData,False,Uij,True)
                for item in result:
                    atom = copy.copy(atomB)
                    atom[cx:cx+3] = item[0]
                    Opr = abs(item[2])%100
                    M = SGData['SGOps'][Opr-1][0]
                    if cmx:
                        opNum = G2spc.GetOpNum(item[2],SGData)
                        mom = np.array(atom[cmx:cmx+3])
                        if SGData['SGGray']:
                            atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)
                        else:
                            atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SpnFlp[opNum-1]
                    atom[cs-1] = str(item[2])+'+'
                    atom[cuij:cuij+6] = item[1]
                    for xyz in cellArray+np.array(atom[cx:cx+3]):
                        dist = np.sqrt(np.sum(np.inner(Amat,xyz-xyzA)**2))
                        if 0 < dist <= radius:
                            if noDuplicate(xyz,atomData):
                                C = xyz-atom[cx:cx+3]+item[3]
                                newAtom = atom[:]
                                newAtom[cx:cx+3] = xyz
                                newAtom[cs-1] += str(int(round(C[0])))+','+str(int(round(C[1])))+','+str(int(round(C[2])))
                                atomData.append(newAtom)
            GoOn = pgbar.Update(ic,newmsg='Centers done=%d'%(ic))
            if not GoOn[0]:
                break
        pgbar.Destroy()
        UpdateDrawAtoms(G2frame,data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def TransformSymEquiv(event):
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms)
        if not indx: return
        indx.sort()
        atomData = data['Drawing']['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        cuij = ci+2
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        atomData = data['Drawing']['Atoms']
        generalData = data['General']
        SGData = generalData['SGData']
        SpnFlp = SGData.get('SpnFlp',[])
        dlg = SymOpDialog(G2frame,SGData,False,True)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Inv,Cent,Opr,Cell,New,Force = dlg.GetSelection()
                Cell = np.array(Cell)
                cent = SGData['SGCen'][Cent]
                M,T = SGData['SGOps'][Opr]
                for ind in indx:
                    XYZ = np.array(atomData[ind][cx:cx+3])
                    XYZ = np.inner(M,XYZ)+T
                    if Inv and not SGData['SGFixed']:
                        XYZ = -XYZ
                    XYZ = XYZ+cent+Cell
                    if Force:
                        XYZ,cell = G2spc.MoveToUnitCell(XYZ)
                        Cell += cell
                    atom = atomData[ind]
                    atom[cx:cx+3] = XYZ
                    OprNum = ((Opr+1)+100*Cent)*(1-2*Inv)
                    if cmx:
                        opNum = G2spc.GetOpNum(OprNum,SGData)
                        mom = np.array(atom[cmx:cmx+3])
                        atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SpnFlp[opNum-1]
                    atomOp = atom[cs-1]
                    newOp = str(((Opr+1)+100*Cent)*(1-2*Inv))+'+'+ \
                        str(int(Cell[0]))+','+str(int(Cell[1]))+','+str(int(Cell[2]))
                    atom[cs-1] = G2spc.StringOpsProd(atomOp,newOp,SGData)
                    if atom[ci] == 'A':
                        Uij = atom[cuij:cuij+6]
                        U = G2spc.Uij2U(Uij)
                        U = np.inner(np.inner(M,U),M)
                        Uij = G2spc.U2Uij(U)
                        atom[cuij:cuij+6] = Uij
                data['Drawing']['Atoms'] = atomData
        finally:
            dlg.Destroy()
        UpdateDrawAtoms(G2frame,data)
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def FillMolecule(event):
        '''This is called by the Complete Molecule command. It adds a layer
        of bonded atoms of the selected types for all selected atoms in
        the Draw Atoms table. If the number of repetitions is greater than
        one, the added atoms (other than H atoms, which are assumed to only
        have one bond) are then searched for the next surrounding layer of
        bonded atoms.
        '''
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms)
        if not indx: return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomTypes,radii = getAtomRadii(data)

        dlg = wx.Dialog(G2frame,wx.ID_ANY,'Addition criteria',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        dlg.CenterOnParent()
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Molecular completion parameters'),0)
        mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Pause after: '),0)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add((45,-1))
        choices = [1,2,5,10,50]
        params = {'maxrep':10, 'maxatm':1000}
        topSizer.Add(G2G.EnumSelector(dlg,params,'maxrep',
                        [str(i) for i in choices],choices))
        topSizer.Add(wx.StaticText(dlg,wx.ID_ANY,' repetitions'),0,WACV)
        mainSizer.Add(topSizer,0)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add((45,-1))
        choices = [100,500,1000,5000,10000]
        topSizer.Add(G2G.EnumSelector(dlg,params,'maxatm',
                        [str(i) for i in choices],choices))
        topSizer.Add(wx.StaticText(dlg,wx.ID_ANY,' added atoms'),0,WACV)
        mainSizer.Add(topSizer,0)
        mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Atom types to add:'),0)
        atSizer = wx.BoxSizer(wx.HORIZONTAL)
        for i,item in enumerate(atomTypes):
            params[item] = True
            atm = G2G.G2CheckBox(dlg,item,params,item)
            atSizer.Add(atm,0,WACV)
        mainSizer.Add(atSizer,0)

        OkBtn = wx.Button(dlg,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
        cancelBtn = wx.Button(dlg,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_CANCEL))
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        dlg.SetSizer(mainSizer)
        dlg.Fit()
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        dlg.Destroy()

        try:
            indH = atomTypes.index('H')
            radii[indH] = 0.5
        except:
            pass
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        neighborArray = FindCoordinationByLabel(data)

        time1 = time.time()
        added = 0
        targets = [item for item in atomTypes if params[item]]
        cx,ct,cs,ci = G2mth.getAtomPtrs(data,draw=True)
        rep = 0
        allrep = 0
        while True:
            if rep == 0:
                pgbar = wx.ProgressDialog('Fill molecular coordination',
                    'Passes done=0 %0',params['maxrep']+1,
                    parent=G2frame,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
                screenSize = wx.ClientDisplayRect()
                Size = pgbar.GetSize()
                if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
                    pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
                    pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
                pgbar.Raise()
                wx.Yield()
            startlen = len(data['Drawing']['Atoms'])
            coordsArray = np.array([a[cx:cx+3] for a in data['Drawing']['Atoms']])
            addedAtoms = []
            for Ind,ind in enumerate(indx):
                addedAtoms += FindCoordination(ind,data,neighborArray,coordsArray,cmx,targets)
                GoOn = pgbar.Update(rep+1,
                    newmsg='Passes done={} atom #{} of {}'
                                        .format(rep+1,Ind+1,len(indx)))
                if not GoOn[0]: break
            if not GoOn[0]: break
            print('pass {} processed {} atoms adding {}; Search time: {:.2f}s'.format(
                allrep+1,len(indx),len(addedAtoms),time.time()-time1))
            time1 = time.time()
            rep += 1
            allrep += 1
            if len(addedAtoms) == 0: break
            added += len(addedAtoms)
            data['Drawing']['Atoms'] += addedAtoms
            # atoms to search over (omit H)
            indx = [i+startlen for i in range(len(addedAtoms)) if addedAtoms[i][ct] != 'H']
            if len(indx) == 0: break
            if added > params['maxatm']:
                msg = "Exceeded number added atoms. Continue?"
                dlg = wx.MessageDialog(G2frame,msg,caption='Continue?',style=wx.YES_NO)
                if dlg.ShowModal() != wx.ID_YES:
                    dlg.Destroy()
                    break
                dlg.Destroy()
                rep = 0
                added = 1
                pgbar.Destroy()
            if rep >= params['maxrep']:
                msg = "Exceeded number of repetitions. Continue?"
                dlg = wx.MessageDialog(G2frame,msg,caption='Continue?',style=wx.YES_NO)
                if dlg.ShowModal() != wx.ID_YES:
                    dlg.Destroy()
                    break
                dlg.Destroy()
                rep = 0
                added = 1
                pgbar.Destroy()
            # UpdateDrawAtoms(G2frame,data)
            # G2plt.PlotStructure(G2frame,data)
        pgbar.Destroy()
        UpdateDrawAtoms(G2frame,data)
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def FillCoordSphere(event=None,selection=None):
        if selection:
            indx = selection
        else:
            indx = getAtomSelections(G2frame.drawAtoms)
        if not indx: return
        time0 = time.time()
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        radii = generalData['BondRadii']
        atomTypes = generalData['AtomTypes']
        try:
            indH = atomTypes.index('H')
            radii[indH] = 0.5
        except:
            pass
        if indx:
            indx.sort()
            atomData = data['Drawing']['Atoms']
            numAtoms = len(atomData)
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
            cij = ci+2
            SGData = generalData['SGData']
            cellArray = G2lat.CellBlock(1)
            nind = len(indx)
            pgbar = wx.ProgressDialog('Fill CN sphere for %d atoms'%nind,'Atoms done=',nind+1,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
            screenSize = wx.ClientDisplayRect()
            Size = pgbar.GetSize()
            if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
                pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
                pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
            for Ind,ind in enumerate(indx):
                atomA = atomData[ind]
                xyzA = np.array(atomA[cx:cx+3])
                indA = atomTypes.index(atomA[ct])
                for atomB in atomData[:numAtoms]:
                    indB = atomTypes.index(atomB[ct])
                    sumR = radii[indA]+radii[indB]
                    xyzB = np.array(atomB[cx:cx+3])
                    for xyz in cellArray+xyzB:
                        dist = np.sqrt(np.sum(np.inner(Amat,xyz-xyzA)**2))
                        if 0 < dist <= data['Drawing']['radiusFactor']*sumR:
                            if noDuplicate(xyz,atomData):
                                oprB = atomB[cs-1]
                                C = xyz-xyzB
                                newOp = '1+'+str(int(round(C[0])))+','+str(int(round(C[1])))+','+str(int(round(C[2])))
                                newAtom = atomB[:]
                                newAtom[cx:cx+3] = xyz
                                newAtom[cs-1] = G2spc.StringOpsProd(oprB,newOp,SGData)
                                atomData.append(newAtom[:cij+9])  #not SS stuff
                GoOn = pgbar.Update(Ind,newmsg='Atoms done=%d'%(Ind))
                if not GoOn[0]:
                    break
            pgbar.Destroy()
            data['Drawing']['Atoms'] = atomData
            print('search time: %.3f'%(time.time()-time0))
            UpdateDrawAtoms(G2frame,data)
            G2frame.drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
        else:
            G2G.G2MessageBox(G2frame,'Select atoms first')

    def FillCoordSphereNew(event):
        time0 = time.time()
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms)
        if not indx: return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomTypes,radii = getAtomRadii(data)
        try:
            indH = atomTypes.index('H')
            radii[indH] = 0.5
        except:
            pass
        indx.sort()
        atomData = data['Drawing']['Atoms']
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        neighborArray = FindCoordinationByLabel(data)
        coordsArray = np.array([a[cx:cx+3] for a in data['Drawing']['Atoms']])
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        nind = len(indx)
        pgbar = wx.ProgressDialog('Fill CN sphere for %d atoms'%nind,'Atoms done=',nind+1,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        for Ind,ind in enumerate(indx):
            atomData += FindCoordination(ind,data,neighborArray,coordsArray,cmx,atomTypes)
            GoOn = pgbar.Update(Ind,newmsg='Atoms done=%d'%(Ind))
            if not GoOn[0]: break
        pgbar.Destroy()
        data['Drawing']['Atoms'] = atomData
        print('search time: %.3f'%(time.time()-time0))
        UpdateDrawAtoms(G2frame,data)
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def FillUnitCell(event,selectAll=None):
        drawAtoms = G2frame.drawAtoms
        if selectAll is not None:
            indx = list(range(drawAtoms.NumberRows))
        else:
            indx = getAtomSelections(drawAtoms)
        if not indx: return
        indx.sort()
        atomData = data['Drawing']['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        cuij = cs+5
        generalData = data['General']
        SGData = generalData['SGData']
        SpnFlp = SGData.get('SpnFlp',[])
        nind = len(indx)
        pgbar = wx.ProgressDialog('Fill unit cell for %d atoms'%nind,'Atoms done=',nind+1,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        for Ind,ind in enumerate(indx):
            atom = atomData[ind]
            XYZ = np.array(atom[cx:cx+3])
            Uij = atom[cuij:cuij+6]
            result = G2spc.GenAtom(XYZ,SGData,False,Uij,True)
            for item in result:
                atom = copy.copy(atomData[ind])
                atom[cx:cx+3] = item[0]
                if cmx:
                    Opr = abs(item[2])%100
                    M = SGData['SGOps'][Opr-1][0]
                    opNum = G2spc.GetOpNum(item[2],SGData)
                    mom = np.array(atom[cmx:cmx+3])
                    if SGData['SGGray']:
                        atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)
                    else:
                        atom[cmx:cmx+3] = np.inner(mom,M)*nl.det(M)*SpnFlp[opNum-1]
                atom[cs-1] = str(item[2])+'+' \
                    +str(item[3][0])+','+str(item[3][1])+','+str(item[3][2])
                atom[cuij:cuij+6] = item[1]
                Opp = G2spc.Opposite(item[0])
                for key in Opp:
                    if noDuplicate(Opp[key],atomData) or atom[ct] == 'Q':
                        unit = item[3]+np.array(eval(key))*1.
                        cell = '%d+%d,%d,%d'%(item[2],unit[0],unit[1],unit[2])
                        atom[cx:cx+3] = Opp[key]
                        atom[cs-1] = cell
                        atomData.append(atom[:cuij+7]+atom[-2:])  #not SS stuff
            data['Drawing']['Atoms'] = atomData
            GoOn = pgbar.Update(Ind,newmsg='Atoms done=%d'%(Ind))
            if not GoOn[0]:
                break
        pgbar.Destroy()
        UpdateDrawAtoms(G2frame,data)
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def DrawAtomsDelete(event):
        indx = getAtomSelections(G2frame.drawAtoms)
        if not indx: return
        indx.sort()
        atomData = data['Drawing']['Atoms']
        indx.reverse()
        for ind in indx:
            del atomData[ind]
        UpdateDrawAtoms(G2frame,data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()

    def SelDrawList(event):
        'select atoms using a filtered listbox'
        drawAtoms = G2frame.drawAtoms
        choices = []
        for i in range(drawAtoms.GetNumberRows()):
            val = drawAtoms.GetCellValue(i,0)
            if val in choices:
                val += '_' + str(i)
            choices.append(val)
        if not choices: return
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select atoms','Choose atoms to select',choices)
        indx = []
        if dlg.ShowModal() == wx.ID_OK:
            indx = dlg.GetSelections()
        dlg.Destroy()
        if len(indx) == 0: return
        drawAtoms.ClearSelection()
        for row in indx:
            drawAtoms.SelectRow(row,True)
        G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()

    def DrawLoadSel(event):
        '''Copy selected atoms from the atoms list into the draw atoms list, making
        sure not to duplicate any.
        '''
        choices = [atm[0] for atm in data['Atoms']]
        dlg = G2G.G2MultiChoiceDialog(G2frame,
                    'Select atoms','Choose atoms to select',choices)
        indx = []
        if dlg.ShowModal() == wx.ID_OK:
            indx = dlg.GetSelections()
        dlg.Destroy()
        if len(indx) == 0: return
        drawingData = data['Drawing']
        cxD,ctD,_,_ = data['Drawing']['atomPtrs']
        cx,ct,cs,cia = data['General']['AtomPtrs']
        atmsXYZ = [np.array(a[cx:cx+3]) for a in data['Atoms']]
        for i in indx:
            found = False
            for dA in drawingData['Atoms']:
                if (dA[ctD] == data['Atoms'][i][ct] and
                    dA[ctD-1] == data['Atoms'][i][ct-1] and
                    dA[cxD+3] == '1' and
                    np.sum((atmsXYZ[i]-dA[cxD:cxD+3])**2) < 0.001):
                    found = True
                    break
            if not found:
                DrawAtomAdd(drawingData,data['Atoms'][i])
        UpdateDrawAtoms(G2frame,data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()

    def OnReloadDrawAtoms(event=None):
        atomData = data['Atoms']
        cx,ct,cs,ci = data['General']['AtomPtrs']
        for atom in atomData:
            ID = atom[ci+8]
            G2mth.DrawAtomsReplaceByID(data,ci+8,atom,ID)
        UpdateDrawAtoms(G2frame,data)
        G2frame.drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
        if event:
            event.StopPropagation()

    def DrawAtomsDeleteByIDs(IDs):
        atomData = data['Drawing']['Atoms']
        cx,ct,cs,ci = data['General']['AtomPtrs']
        loc = ci+8
        indx = G2mth.FindAtomIndexByIDs(atomData,loc,IDs)
        indx.reverse()
        for ind in indx:
            del atomData[ind]

    def ChangeDrawAtomsByIDs(colName,IDs,value):
        atomData = data['Drawing']['Atoms']
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        if colName == 'Name':
            col = ct-1
        elif colName == 'Type':
            col = ct
        elif colName == 'I/A':
            col = cs
        indx = G2mth.FindAtomIndexByIDs(atomData,ci+8,IDs)
        for ind in indx:
            atomData[ind][col] = value

    def OnDrawPlane(event):
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms)
        if len(indx) < 4:
            G2G.G2MessageBox(G2frame,'Select four or more atoms first')
            print ('**** ERROR - need 4+ atoms for plane calculation')
            return
        PlaneData = {}
        drawingData = data['Drawing']
        atomData = drawingData['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx = colLabels.index('x')
        cn = colLabels.index('Name')
        xyz = []
        for i,atom in enumerate(atomData):
            if i in indx:
                xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
        generalData = data['General']
        PlaneData['Name'] = generalData['Name']
        PlaneData['Atoms'] = xyz
        PlaneData['Cell'] = generalData['Cell'][1:] #+ volume
        G2stMn.BestPlane(PlaneData)

    def OnDrawDistVP(event):
        # distance to view point
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms,action='distance calc')
        if not indx: return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        drawingData = data['Drawing']
        viewPt = np.array(drawingData['viewPoint'][0])
        print (' Distance from view point at %.3f %.3f %.3f to:'%(viewPt[0],viewPt[1],viewPt[2]))
        atomDData = drawingData['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx = colLabels.index('x')
        cn = colLabels.index('Name')
        cs = colLabels.index('Sym Op')
        for i in indx:
            atom = atomDData[i]
            Dx = np.array(atom[cx:cx+3])-viewPt
            dist = np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0))
            print ('Atom: %8s (%12s) distance = %.3f'%(atom[cn],atom[cs],dist))

    def OnDrawDAT(event):
        #compute distance, angle, or torsion depending on number of selections
        drawAtoms = G2frame.drawAtoms
        indx = getAtomSelections(drawAtoms)
        if len(indx) not in [2,3,4]:
            G2G.G2MessageBox(G2frame,'Select 2, 3 or 4 atoms first')
            print ('**** ERROR - wrong number of atoms for distance, angle or torsion calculation')
            return
        DATData = {}
        ocx,oct,ocs,cia = data['General']['AtomPtrs']
        drawingData = data['Drawing']
        atomData = data['Atoms']
        atomDData = drawingData['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx = colLabels.index('x')
        cn = colLabels.index('Name')
        cid = colLabels.index('I/A')+8
        xyz = []
        Oxyz = []
        DATData['Natoms'] = len(indx)
        for i in indx:
            atom = atomDData[i]
            xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+4]) #also gets Sym Op
            Id = G2mth.FindAtomIndexByIDs(atomData,cid,[atom[cid],],False)[0]
            Oxyz.append([Id,]+atomData[Id][cx+1:cx+4])
        DATData['Datoms'] = xyz
        DATData['Oatoms'] = Oxyz
        generalData = data['General']
        DATData['Name'] = generalData['Name']
        DATData['SGData'] = generalData['SGData']
        DATData['Cell'] = generalData['Cell'][1:] #+ volume
        if 'pId' in data:
            DATData['pId'] = data['pId']
            DATData['covData'] = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Covariance'))
        G2stMn.DisAglTor(DATData)

    def MapVoid(event):
        generalData = data['General']
        rMax = max(data['General']['vdWRadii'])
        cell = data['General']['Cell'][1:7]
        drawingData = data['Drawing']
        voidDlg = wx.Dialog(G2frame,wx.ID_ANY,'Void computation parameters',style=wx.DEFAULT_DIALOG_STYLE)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(voidDlg,wx.ID_ANY,
                    'Set parameters for void computation for phase '+generalData.get('Name','?')))
        # get cell ranges
        xmax = 2. - rMax/cell[0]
        voidPar = {'a':1., 'b':1., 'c':1., 'grid':.25, 'probe':0.5}
        for i in ('a', 'b', 'c'):
            mainSizer.Add(G2G.G2SliderWidget(voidDlg,voidPar,i,'Max '+i+' value: ',0.,xmax,100))
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(wx.StaticText(voidDlg,wx.ID_ANY,'Grid spacing (A)'))
        hSizer.Add(G2G.ValidatedTxtCtrl(voidDlg,voidPar,'grid',nDig=(5,2), xmin=0.1, xmax=2., typeHint=float))
        mainSizer.Add(hSizer)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(wx.StaticText(voidDlg,wx.ID_ANY,'Probe radius (A)'))
        hSizer.Add(G2G.ValidatedTxtCtrl(voidDlg,voidPar,'probe',nDig=(5,2), xmin=0.1, xmax=2., typeHint=float))
        mainSizer.Add(hSizer)

        def OnOK(event): voidDlg.EndModal(wx.ID_OK)
        mainSizer.Add([5,5])
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(voidDlg, wx.ID_OK)
        btn.Bind(wx.EVT_BUTTON, OnOK)
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(voidDlg, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)

        voidDlg.SetSizer(mainSizer)
        mainSizer.Fit(voidDlg)
        voidDlg.CenterOnParent()
        res = voidDlg.ShowModal()
        voidDlg.Destroy()
        if res != wx.ID_OK: return
        drawingData['Voids'] = VoidMap(data, voidPar['a'], voidPar['b'], voidPar['c'],
            voidPar['grid'],voidPar['probe'])
        drawingData['showVoids'] = True
        G2plt.PlotStructure(G2frame,data)

    def RandomizedAction(event):
        '''perform a selected action on a random sequence from a selected
        list of atoms. After each selection, one chooses if the
        action should be performed on the selected atom
        '''
        ranDrwDict['opt'] = 1
        ranDrwDict['optList'] = ['Delete selection',
                          'Change color of selection',
                          'Atom drawing style for selection',
                          'Add sphere of atoms around selection',
                          'Add coordination-sphere around selection']

        dlg = G2G.G2SingleChoiceDialog(G2frame,'Select option from list',
                                    'Select option',ranDrwDict['optList'])
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                ranDrwDict['opt'] = dlg.GetSelection()
                ranDrwDict['optList'][dlg.GetSelection()]
            else:
                return
        finally:
            dlg.Destroy()
        if ranDrwDict['opt'] == 1:
            colors = wx.ColourData()
            colors.SetChooseFull(True)
            try:
                dlg = wx.ColourDialog(G2frame.GetParent(),colors)
                if dlg.ShowModal() == wx.ID_OK:
                    ranDrwDict['color'] = dlg.GetColourData().GetColour()[:3]
                else:
                    return
            finally:
                dlg.Destroy()
        elif ranDrwDict['opt'] == 2:
            styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']
            if data['General']['Type'] == 'macromolecular':
                styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids',
                                   'backbone','ribbons','schematic']
            try:
                dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom drawing style',
                                        styleChoice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    ranDrwDict['style'] = styleChoice[sel]
                else:
                    return
            finally:
                dlg.Destroy()
        elif ranDrwDict['opt'] == 3:
            dlg = SphereEnclosure(G2frame,data['General'],data['Drawing'],None)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    centers,ranDrwDict['radius'],ranDrwDict['targets'] = dlg.GetSelection()
                    ranDrwDict['2call'] = AddSphere
                else:
                    return
            finally:
                dlg.Destroy()
        elif ranDrwDict['opt'] == 4:
            ranDrwDict['2call'] = FillCoordSphere

        indx = getAtomSelections(G2frame.drawAtoms)
        if not indx: return
        ranDrwDict['atomList'] = list(indx)
        ran.shuffle(ranDrwDict['atomList'])
        i = ranDrwDict['atomList'][0]
        cx,ct = data['Drawing']['atomPtrs'][:2]
        data['Drawing']['viewPoint'][0] = data['Drawing']['Atoms'][i][cx:cx+3]
        G2frame.drawAtoms.SelectRow(i)
        G2plt.PlotStructure(G2frame,data)
        msg = f"Atom #{i} selected ({data['Drawing']['Atoms'][i][ct-1]})"
        print(msg)
        ranDrwDict['delAtomsList'] = []
        ranDrwDict['msgWin'] = wx.Frame(G2frame, wx.ID_ANY, size=(300, 300),
                        style=wx.DEFAULT_FRAME_STYLE|wx.FRAME_FLOAT_ON_PARENT|
                                    wx.STAY_ON_TOP)
        ranDrwDict['msgWin'].SetTitle("Random action messages")
        siz = wx.BoxSizer(wx.VERTICAL)
        ranDrwDict['msgWin'].text1 = wx.StaticText(ranDrwDict['msgWin'],  wx.ID_ANY,
            f"For random Draw Atoms action, in plot window press:  \n\t'Y' to {ranDrwDict['optList'][ranDrwDict['opt']]}\n\t'N' to advance to the next atom\n\t'Q' to end\n"
            )
        siz.Add(ranDrwDict['msgWin'].text1)
        ranDrwDict['msgWin'].text2 = wx.StaticText(ranDrwDict['msgWin'],  wx.ID_ANY, msg)
        siz.Add(ranDrwDict['msgWin'].text2)
        ranDrwDict['msgWin'].SetSizer(siz)
        siz.Fit(ranDrwDict['msgWin'])
        ranDrwDict['msgWin'].CentreOnParent()
        ranDrwDict['msgWin'].Show()
        UpdateDrawAtoms(G2frame,data)

# #### Draw Options page ################################################################################
    def UpdateDrawOptions(G2frame,data):
        def SlopSizer():
            def OnCameraPos():
                #old code
                # drawingData['cameraPos'] = cameraPos.GetValue()
                # cameraPosTxt.SetLabel(' Camera Distance: '+'%.2f'%(drawingData['cameraPos']))
                # Zclip.SetLabel(' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
                #new code
#                drawingData['Zclip'] = min(drawingData['Zclip'],0.95*drawingData['cameraPos'])
                Zclip.SetScaledValue(drawingData['Zclip'])
                Zval.ChangeValue(drawingData['Zclip'])
                xmin=1.0    #.01*drawingData['Zclip']*drawingData['cameraPos']/100.
                xmax=2.*drawingData['cameraPos']
                Zclip.SetScaledRange(xmin,xmax)
                Zclip.SetMax(xmax)
                Zval.Validator.xmin = xmin
                Zval.Validator.xmax = xmax
                #end new code
                G2plt.PlotStructure(G2frame,data)

            def OnMoveZ(event):
                #old code
                # drawingData['Zclip'] = Zclip.GetValue()
                # Zclip.SetLabel(' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
                #new code
                move = MoveZ.GetValue()*drawingData['Zstep']
                MoveZ.SetValue(0)
                VP = np.inner(Amat,np.array(drawingData['viewPoint'][0]))
                VD = np.inner(Amat,np.array(drawingData['viewDir']))
                VD /= np.sqrt(np.sum(VD**2))
                VP += move*VD
                VP = np.inner(Bmat,VP)
                drawingData['viewPoint'][0] = VP
                panel = drawOptions.GetChildren()
                names = [child.GetName() for child in panel]
                panel[names.index('viewPoint')].SetValue('%.3f %.3f %.3f'%(VP[0],VP[1],VP[2]))
                #end new code
                G2plt.PlotStructure(G2frame,data)

            def OnRadFactor(invalid,value,tc):
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data)

            slopSizer = wx.BoxSizer(wx.HORIZONTAL)
            slideSizer = wx.FlexGridSizer(0,3,0,0)
            slideSizer.AddGrowableCol(2,1)
            valSize = (50,20)

            cameraPosTxt,cameraPos = G2G.G2SliderWidget(
                drawOptions,drawingData,'cameraPos',
                sizer=slideSizer,
                nDig=(10,1),xmin=5.,xmax=500.,size=valSize,
                label=' Camera Distance, '+Angstr+': ',iscale=5.,
                onChange=OnCameraPos)
            G2frame.phaseDisplay.cameraPosTxt = cameraPosTxt
            G2frame.phaseDisplay.cameraSlider = cameraPos

            Zval,Zclip = G2G.G2SliderWidget(
                drawOptions,drawingData,'Zclip',
                sizer=slideSizer,size=valSize,
                xmin=1.0,    #.01*drawingData['Zclip']*drawingData['cameraPos']/100.,
                xmax=2.*drawingData['cameraPos'],
                nDig=(10,2),
                label=' Z clipping, '+Angstr+': ',iscale=50.,
                onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])
            G2frame.phaseDisplay.Zval = Zval
            G2frame.phaseDisplay.Zclip = Zclip
            OnCameraPos()

            slideSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,' Z step, '+Angstr+': '),0,WACV)
            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0,size=valSize)
            slideSizer.Add(Zstep,0,WACV)
            MoveSizer = wx.BoxSizer(wx.HORIZONTAL)
            MoveSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,'   Press to step:'),0,WACV)
            MoveZ = wx.SpinButton(drawOptions,style=wx.SP_HORIZONTAL,size=valSize)
            MoveZ.SetValue(0)
            MoveZ.SetRange(-1,1)
            MoveZ.Bind(wx.EVT_SPIN, OnMoveZ)
            MoveSizer.Add(MoveZ)
            slideSizer.Add(MoveSizer,1,wx.EXPAND|wx.RIGHT)

            G2G.G2SliderWidget(
                drawOptions,drawingData,'vdwScale',
                sizer=slideSizer,size=valSize,
                xmin=0.01, xmax=1.0,
                nDig=(10,3),iscale=50.,
                label=' van der Waals scale: ',
                onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])
            G2G.G2SliderWidget(
                drawOptions,drawingData,'ellipseProb',
                sizer=slideSizer,size=valSize,
                xmin=1, xmax=99,
                nDig=(10,2),iscale=2.,
                label=' Ellipsoid probability, %: ',
                onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])
            G2G.G2SliderWidget(
                drawOptions,drawingData,'ballScale',
                sizer=slideSizer,size=valSize,
                xmin=0.01, xmax=0.99,
                nDig=(10,3),iscale=100.,
                label=' Ball scale: ',
                onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])
            G2G.G2SliderWidget(
                drawOptions,drawingData,'bondRadius',
                sizer=slideSizer,size=valSize,
                xmin=0.01, xmax=0.25,
                nDig=(10,2),iscale=500.,
                label=' Bond radius, '+Angstr+': ',
                onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])
            if generalData['Type'] == 'magnetic':
                G2G.G2SliderWidget(
                    drawOptions,drawingData,'magMult',
                    sizer=slideSizer,size=valSize,
                    xmin=0.1, xmax=1.2,
                    nDig=(10,2),iscale=100.,
                    label=' Mag. mom. mult.: ',
                    onChange=G2plt.PlotStructure,onChangeArgs=[G2frame,data])

            slideSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,' Bond search factor: '),0,WACV)
            slideSizer.Add(G2G.ValidatedTxtCtrl(drawOptions,drawingData,'radiusFactor',
                nDig=(10,2),xmin=0.1,xmax=1.2,size=valSize,OnLeave=OnRadFactor),0,WACV)
            slideSizer.Add((-1,-1))

            slopSizer.Add(slideSizer,1,wx.EXPAND|wx.RIGHT)
            slopSizer.Add((10,5),0)
            slopSizer.SetMinSize(wx.Size(350,10))
            return slopSizer

        def ShowSizer():
            def OnShowABC(event):
                drawingData['showABC'] = showABC.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnShowUnitCell(event):
                drawingData['unitCellBox'] = unitCellBox.GetValue()
                G2plt.PlotStructure(G2frame,data)
                G2frame.GPXtree.UpdateSelection()
#                wx.CallAfter(UpdateDrawOptions)

            def OnShowHyd(event):
                drawingData['showHydrogen'] = showHydrogen.GetValue()
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data)

            def OnShowRB(event):
                drawingData['showRigidBodies'] = showRB.GetValue()
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data)

            def OnSymFade(event):
                drawingData['SymFade'] = symFade.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnShowVoids(event):
                drawingData['showVoids'] = showVoids.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnViewPoint(event):
                event.Skip()
                Obj = event.GetEventObject()
                viewPt = Obj.GetValue().split()
                try:
                    VP = [float(viewPt[i]) for i in range(3)]
                except (ValueError,IndexError):
                    VP = drawingData['viewPoint'][0]
                Obj.SetValue('%.3f %.3f %.3f'%(VP[0],VP[1],VP[2]))
                drawingData['viewPoint'][0] = VP
                G2plt.PlotStructure(G2frame,data)

            def OnViewDir(event):
                event.Skip()
                Obj = event.GetEventObject()
                viewDir = Obj.GetValue().split()
                try:
                    Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
#reset view to stndard
                    drawingData['viewDir'] = [0,0,1]
                    drawingData['oldxy'] = []
                    V0 = np.array([0,0,1])
                    V = np.inner(Amat,V0)
                    V /= np.sqrt(np.sum(V**2))
                    A = np.arccos(np.sum(V*V0))
                    Q = G2mth.AV2Q(A,[0,1,0])
                    drawingData['Quaternion'] = Q
#new view made here
                    VD = np.array([float(viewDir[i]) for i in range(3)])
                    VC = np.inner(Amat,VD)
                    VC /= np.sqrt(np.sum(VC**2))
                    V = np.array(drawingData['viewDir'])
                    VB = np.inner(Amat,V)
                    VB /= np.sqrt(np.sum(VB**2))
                    VX = np.cross(VC,VB)
                    A = acosd(max((2.-np.sum((VB-VC)**2))/2.,-1.))
                    if A == 0.0 and len(viewDir) == 3:
                        raise ValueError        #avoid a no op
#                    print('\nnew view =',viewDir)
#                    print('A=%.3f, V='%A,VX)
                    QV = G2mth.AVdeg2Q(A,VX)
                    if len(viewDir) > 3:
                        Model = drawingData['modelView'][:3,:3]
                        invModel = nl.inv(Model)
                        rt2 = np.sqrt(2.)/2.
                        VX0 = np.array([-1.,0.,0.])
                        VY0 = np.array([0.,-1.,0.])
                        if 'H' == viewDir[3].upper():
                            QV = G2mth.prodQQ(np.array([rt2,0.,rt2,0.]),QV)     #rotate 90deg about +ve Y
                            VD = np.inner(invModel.T,VX0)
                        elif 'V' == viewDir[3].upper():
                            QV = G2mth.prodQQ(np.array([rt2,-rt2,0.,0.]),QV)     #rotate 90deg about -ve X
                            VD = np.inner(invModel.T,VY0)
#                        NAV = G2mth.Q2AVdeg(QV)
#                        print('dir = %s, A= %.3f, V = '%(viewDir[3],NAV[0]),NAV[1:])
                        VD /= np.sqrt(np.sum(VD**2))
#                        print('new view dir = ',VD)
                    Q = drawingData['Quaternion']
                    QN = G2mth.prodQQ(Q,QV)
                    drawingData['Quaternion'] = QN
                except (ValueError,IndexError):
                    VD = drawingData['viewDir']
                Obj.SetValue('%.3f %.3f %.3f'%(VD[0],VD[1],VD[2]))
                drawingData['viewDir'] = VD
                G2plt.PlotStructure(G2frame,data)

            showSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(drawOptions,label=' Background color:'),0,WACV)
            backColor = G2G.setColorButton(drawOptions,drawingData, 'backColor',
                                       G2plt.PlotStructure,[G2frame,data])
            lineSizer.Add(backColor,0,WACV)
            lineSizer.Add(wx.StaticText(drawOptions,-1,' View Dir.:'),0,WACV)
            VD = drawingData['viewDir']
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            viewDir = wx.TextCtrl(drawOptions,value='%.3f %.3f %.3f'%(VD[0],VD[1],VD[2]),
                style=wx.TE_PROCESS_ENTER,size=wx.Size(140,20),name='viewDir')
            viewDir.Bind(wx.EVT_TEXT_ENTER,OnViewDir)
#            viewDir.Bind(wx.EVT_KILL_FOCUS,OnViewDir)
            G2frame.phaseDisplay.viewDir = viewDir
            lineSizer.Add(viewDir,0,WACV)
            showSizer.Add(lineSizer)
            showSizer.Add((0,5),0)

            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            showABC = wx.CheckBox(drawOptions,-1,label=' Show view point?')
            showABC.Bind(wx.EVT_CHECKBOX, OnShowABC)
            showABC.SetValue(drawingData['showABC'])
            lineSizer.Add(showABC,0,WACV)
            lineSizer.Add(wx.StaticText(drawOptions,-1,' View Point:'),0,WACV)
            VP = drawingData['viewPoint'][0]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            viewPoint = wx.TextCtrl(drawOptions,value='%.3f %.3f %.3f'%(VP[0],VP[1],VP[2]),
                style=wx.TE_PROCESS_ENTER,size=wx.Size(140,20),name='viewPoint')
            G2frame.phaseDisplay.viewPoint = viewPoint
            viewPoint.Bind(wx.EVT_TEXT_ENTER,OnViewPoint)
            viewPoint.Bind(wx.EVT_KILL_FOCUS,OnViewPoint)
            lineSizer.Add(viewPoint,0,WACV)
            showSizer.Add(lineSizer)
            mapSizer = DistanceSettingSizer(drawingData,
                        'VPPeakDistRad','VPatomsExpandRad','VPatomsDistRad')
            showSizer.Add(mapSizer,0,wx.LEFT,20)
            showSizer.Add((0,5),0)

            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)

            unitCellBox = wx.CheckBox(drawOptions,-1,label=' Show unit cell?')
            unitCellBox.Bind(wx.EVT_CHECKBOX, OnShowUnitCell)
            unitCellBox.SetValue(drawingData['unitCellBox'])
            line2Sizer.Add(unitCellBox,0,WACV)

            showHydrogen = wx.CheckBox(drawOptions,-1,label=' Show hydrogens?')
            showHydrogen.Bind(wx.EVT_CHECKBOX, OnShowHyd)
            showHydrogen.SetValue(drawingData['showHydrogen'])
            line2Sizer.Add(showHydrogen,0,WACV)

            showRB = wx.CheckBox(drawOptions,-1,label=' Show Rigid Bodies?')
            showRB.Bind(wx.EVT_CHECKBOX, OnShowRB)
            showRB.SetValue(drawingData['showRigidBodies'])
            line2Sizer.Add(showRB,0,WACV)

            showSizer.Add(line2Sizer)

            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            symFade = wx.CheckBox(drawOptions,-1,label=' Fade sym equivs?')
            symFade.Bind(wx.EVT_CHECKBOX, OnSymFade)
            symFade.SetValue(drawingData['SymFade'])
            line3Sizer.Add(symFade,0,WACV)
            showVoids = wx.CheckBox(drawOptions,-1,label=' Show void map?')
            showVoids.Bind(wx.EVT_CHECKBOX, OnShowVoids)
            showVoids.SetValue(drawingData['showVoids'])
            line3Sizer.Add(showVoids,0,WACV)
            showSizer.Add(line3Sizer)

            return showSizer

        def MapSizer():

            def OnShowMap(event):
                drawingData['showMap'] = showMap.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnShowSlice(event):
                drawingData['showSlice'] = G2frame.phaseDisplay.showCS.GetSelection()
                G2frame.phaseDisplay.showCS.SetValue(slices[drawingData['showSlice']])
                G2plt.PlotStructure(G2frame,data)

            def OnSliceSize(invalid,value,tc):
                G2plt.PlotStructure(G2frame,data)

            def OnContourMax(event):
                drawingData['contourMax'] = contourMax.GetValue()/100.
                contourMaxTxt.SetLabel(' Max.: '+'%.2f'%(drawingData['contourMax']*generalData['Map']['rhoMax']))
                G2plt.PlotStructure(G2frame,data)

            mapSizer = wx.BoxSizer(wx.VERTICAL)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            slices = ['','lines','colors','lines+colors']
            line3Sizer.Add(wx.StaticText(drawOptions,label=' Show map slice as '),0,WACV)
            G2frame.phaseDisplay.showCS = wx.ComboBox(drawOptions,value=slices[drawingData['showSlice']],
                choices=slices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            G2frame.phaseDisplay.showCS.Bind(wx.EVT_COMBOBOX, OnShowSlice)
            G2frame.phaseDisplay.showCS.SetValue(slices[drawingData['showSlice']])
            line3Sizer.Add(G2frame.phaseDisplay.showCS,0,WACV)
            line3Sizer.Add(wx.StaticText(drawOptions,label=' Slice size 2X(2-20)A: '),0,WACV)
            line3Sizer.Add(G2G.ValidatedTxtCtrl(drawOptions,drawingData,'sliceSize',nDig=(10,2),xmin=2.0,xmax=20.0,OnLeave=OnSliceSize),0,WACV)
            mapSizer.Add(line3Sizer)
            line4Sizer = wx.BoxSizer(wx.HORIZONTAL)
            contourMaxTxt = wx.StaticText(drawOptions,label=' Max.: '+'%.2f'%(drawingData['contourMax']*generalData['Map']['rhoMax']))
            line4Sizer.Add(contourMaxTxt,0,WACV)
            contourMax = G2G.G2Slider(drawOptions,style=wx.SL_HORIZONTAL,size=(150,25),
                value=int(100*drawingData['contourMax']),minValue=1,maxValue=100)
            contourMax.Bind(wx.EVT_SLIDER, OnContourMax)
            line4Sizer.Add(contourMax,1,wx.EXPAND|wx.RIGHT)
            mapSizer.Add(line4Sizer)
            valSize = (50,20)
            showMap = wx.CheckBox(drawOptions,label=' Show density map?')
            showMap.Bind(wx.EVT_CHECKBOX, OnShowMap)
            showMap.SetValue(drawingData['showMap'])
            mapSizer.Add(showMap,0)
            sliders = wx.FlexGridSizer(0,3,5,5)
            G2G.G2SliderWidget(drawOptions,drawingData,'contourLevel',
                'Fraction of rho max ({:.2f}): '.format(generalData['Map']['rhoMax']),0.01,1.0,100.,
                sizer=sliders,size=valSize,onChange=G2plt.PlotStructure,onChangeArgs=(G2frame,data))
            G2G.G2SliderWidget(drawOptions,drawingData,'mapSize',
                'Visible map radius: ',0.1,10.,10.,sizer=sliders,size=valSize,
                onChange=G2plt.PlotStructure,onChangeArgs=(G2frame,data))
            mapSizer.Add(sliders)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,'On map peak selection:  '),0,WACV)
            lineSizer.Add(G2G.G2CheckBox(drawOptions,'Move view point',drawingData,'peakMoveView'))
            mapSizer.Add(lineSizer,0)
            mapSizer.Add(DistanceSettingSizer(drawingData,
                'PeakDistRadius','atomsExpandRadius','atomsDistRadius'),0,wx.LEFT,20)
            return mapSizer

        def PlaneSizer():

            def OnPlane(event):
                event.Skip()
                vals = plane.GetValue().split()
                try:
                    hkl = [float(vals[i]) for i in range(3)]
                    if not any(hkl):       #can't be all zeros!
                        raise ValueError
                except (ValueError,IndexError):
                    hkl = drawingData['Plane'][0]
                drawingData['Plane'][0] = hkl
                plane.SetValue('%5.3f %5.3f %5.3f'%(hkl[0],hkl[1],hkl[2]))
                G2plt.PlotStructure(G2frame,data)

            def OnShowPlane(event):
                drawingData['Plane'][1] = showPlane.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnShowStack(event):
                drawingData['Plane'][2] = showStack.GetValue()
                G2plt.PlotStructure(G2frame,data)

            def OnPhase(invalid,value,tc):
                G2plt.PlotStructure(G2frame,data)

            planeSizer = wx.BoxSizer(wx.VERTICAL)
            planeSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            planeSizer1.Add(wx.StaticText(drawOptions,label=' Plane: '),0,WACV)
            H = drawingData['Plane'][0]
            plane = wx.TextCtrl(drawOptions,value='%5.3f %5.3f %5.3f'%(H[0],H[1],H[2]),
                style=wx.TE_PROCESS_ENTER,size=(140,20))
            plane.Bind(wx.EVT_TEXT_ENTER,OnPlane)
            plane.Bind(wx.EVT_KILL_FOCUS,OnPlane)
            planeSizer1.Add(plane,0,WACV)
            showPlane = wx.CheckBox(drawOptions,label=' Show plane?')
            showPlane.SetValue(drawingData['Plane'][1])
            showPlane.Bind(wx.EVT_CHECKBOX, OnShowPlane)
            planeSizer1.Add(showPlane,0,WACV)
            showStack = wx.CheckBox(drawOptions,label=' As a stack?')
            showStack.SetValue(drawingData['Plane'][2])
            showStack.Bind(wx.EVT_CHECKBOX, OnShowStack)
            planeSizer1.Add(showStack,0,WACV)
            planeSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            planeSizer2.Add(wx.StaticText(drawOptions,label=' Phase shift (deg): '),0,WACV)
            phase = G2G.ValidatedTxtCtrl(drawOptions,drawingData['Plane'],3,nDig=(10,2),OnLeave=OnPhase)
            planeSizer2.Add(phase,0,WACV)
            planeSizer2.Add(wx.StaticText(drawOptions,-1,' Plane color: '),0,WACV)
            planeColor = G2G.setColorButton(drawOptions,drawingData['Plane'], 4,
                                       G2plt.PlotStructure,[G2frame,data])
            planeSizer2.Add(planeColor,0,WACV)
            planeSizer.Add(planeSizer1)
            planeSizer.Add(planeSizer2)
            return planeSizer

        def DistanceSettingSizer(var,key1,key2,key3):
            '''Sizer to get distances to show'''
            def onLeave(*args,**kwargs):
                G2plt.PlotStructure(G2frame,data)
            for key in key1,key2,key3:
                if key not in var: var[key] = 0.
            mapSizer = wx.FlexGridSizer(0,3,5,5)
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,'Show Map points within:'),0,WACV)
            mapSizer.Add(G2G.ValidatedTxtCtrl(drawOptions,var,key1,
                xmin=0.0,xmax=5.0,nDig=(10,1),size=(50,-1),
                typeHint=float,OnLeave=onLeave))
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,u"\u212B"),0,WACV)
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,'Show atoms within:'),0,WACV)
            mapSizer.Add(G2G.ValidatedTxtCtrl(drawOptions,var,key2,
                xmin=0.0,xmax=15.0,nDig=(10,1),size=(50,-1),
                typeHint=float,OnLeave=onLeave))
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,u"\u212B"),0,WACV)
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,'Label distance to atoms within:'),0,WACV)
            mapSizer.Add(G2G.ValidatedTxtCtrl(drawOptions,var,key3,
                xmin=0.0,xmax=15.0,nDig=(10,1),size=(50,-1),
                typeHint=float,OnLeave=onLeave))
            mapSizer.Add(wx.StaticText(drawOptions,wx.ID_ANY,u"\u212B"),0,WACV)
            return mapSizer

        # UpdateDrawOptions exectable code starts here
        drawOptions = G2frame.drawOptions
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Drawing controls for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        SetupDrawingData()
        drawingData = data['Drawing']
        SetDrawingDefaults(drawingData)

        G2frame.GetStatusBar().SetStatusText('Add h or v to View Dir to set vector horizontal or vertical',1)
        if drawOptions.GetSizer():
            drawOptions.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(SlopSizer(),0)
        G2G.HorizontalLine(mainSizer,drawOptions)
        mainSizer.Add(ShowSizer(),0,)
        if generalData['Map']['rhoMax'] and not generalData.get('4DmapData',{}):
            G2G.HorizontalLine(mainSizer,drawOptions)
            mainSizer.Add(MapSizer())
        G2G.HorizontalLine(mainSizer,drawOptions)
        mainSizer.Add(PlaneSizer(),0,)

        SetPhaseWindow(drawOptions,mainSizer)

####  Deformation form factor routines ################################################################

    def SetDefDist(event):
        generalData = data['General']
        DisAglCtls = {}
        if 'DisAglCtls' in generalData:
            DisAglCtls = generalData['DisAglCtls']
        dlg = G2G.DisAglDialog(G2frame,DisAglCtls,generalData,Angle=False)
        if dlg.ShowModal() == wx.ID_OK:
            generalData['DisAglCtls'] = dlg.GetData()
        G2phsG2.UpdateDeformation(G2frame,data,None)
        event.StopPropagation()

    def SelDeformAtom(event):
        'select deformation atom using a filtered listbox'
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        choices = []
        types = []
        Ids = []
        sytsyms = []
        for atom in data['Atoms']:
            if atom[ct] in atmdata.OrbFF and atom[cia+8] not in data['Deformations']:
                choices.append(atom[ct-1])
                types.append(atom[ct])
                Ids.append(atom[cia+8])
                sytsyms.append(atom[cs])
        if not choices: return      #no atoms in phase!
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select atom','Choose atom to select',choices)
        indxes = []
        if dlg.ShowModal() == wx.ID_OK:
            indxes = dlg.GetSelections()
            for indx in indxes:
                fxchoice = []
                orbs = atmdata.OrbFF[types[indx]]
                data['Deformations'][Ids[indx]] = []
                newj0 = True
                newjn = True
                radial = ''
                for orb in orbs:
                    if 'core' in orb:
                        continue        #skip core - has no parameters
                    else:
                        if 'j0' in orb:
                            if newj0:
                                data['Deformations'][Ids[indx]].append([orb,{'Ne':[float(orbs[orb]['Ne']),False],'kappa':[1.0,False]}])   #no sp. harm for j0 terms
                                if 'Bessel' not in fxchoice:
                                    fxchoice.append('Bessel')
                                radial = 'Bessel'
                                newj0 = False
                            else:
                                data['Deformations'][Ids[indx]].append([orb,{'Ne':[float(orbs[orb]['Ne']),False]}])   #no sp. harm for j0 terms; one kappa only
                        elif 'j' in orb:
                            if newjn:
                                orbDict = {'kappa':[1.0,False],}
                                newjn = False
                            else:
                                orbDict = {}
                            Order = int(orb.split('>')[0][-1])
                            cofNames,cofSgns = G2lat.GenRBCoeff(sytsyms[indx],'1',Order)
                            cofNames = [name.replace('C','D') for name in cofNames]
                            cofTerms = {name:[0.0,False] for name in cofNames if str(Order) in name}
                            for name in cofNames:
                                if str(Order) in name and '0' not in name:
                                    negname = name.replace(',',',-')
                                    cofTerms.update({negname:[0.0,False]})
                            orbDict.update(cofTerms)
                            data['Deformations'][Ids[indx]].append([orb,orbDict])
                        elif 'Sl ' in orb:
                            if 'Slater' not in fxchoice:
                                fxchoice.append('Slater')
                                if not radial:
                                    radial = 'Slater'
                            if 'Sl val' in orb: #valence; no harmonics
                                data['Deformations'][Ids[indx]].append([orb,{'Ne':[float(orbs[orb]['Ne']),False],'kappa':[1.0,False]}])
                                break
                            # else:   #p, d or f
                            #     orbDict = {}
                            #     Order = 'spdf'.index(orb[-1])
                            #     cofNames,cofSgns = G2lat.GenRBCoeff(sytsyms[indx],'1',Order)      #sytsym, RBsym = '1'
                            #     cofNames = [name.replace('C','D') for name in cofNames]
                            #     cofTerms = {name:[0.0,False] for name in cofNames if str(Order) in name}
                            #     for name in cofNames:
                            #         if str(Order) in name and '0' not in name:
                            #             negname = name.replace(',',',-')
                            #             cofTerms.update({negname:[0.0,False]})
                            #     orbDict.update(cofTerms)
                            #     orbDict.update({'Ne':[float(orbs[orb]['Ne']),False]})
                            #     data['Deformations'][Ids[indx]].append([orb,orbDict])
                data['Deformations'][-Ids[indx]] = {'U':'X','V':'Y','UVmat':np.eye(3),
                    'MUV':"A: X'=U, Y'=(UxV)xU & Z'=UxV",'Radial':radial,'fxchoice':fxchoice}
        dlg.Destroy()
        if not len(indxes):
            return
        G2frame.drawAtoms.ClearSelection()
        G2frame.drawAtoms.SelectRow(indx,True)
        G2plt.PlotStructure(G2frame,data)
        G2phsG2.UpdateDeformation(G2frame,data,None)
        event.StopPropagation()

#### DData routines - GUI stuff in GSASIIddataGUI.py. Used for Phase/data "Edit Phase" menu
    def OnHklfAdd(event):
        '''Called to link a Single Xtal (HKLF) dataset to the current phase.
        Most commonly, the histogram and phase are linked when the latter
        item is read in (routines OnImportPhase or OnImportSfact in
        func:`GSASIIdataGUI.GSASIImain`, but one can defer this or change
        the linking later using this routine.

        Note that the capability here is duplicated in routine OnHklfAdd
        inside :func:`GSASIIddataGUI.MakeHistPhaseWin`.
        '''
        result = CheckAddHKLF(G2frame,data)
        if result is None: return
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)

    def OnDataUse(event):
#        hist = G2frame.hist
        if data['Histograms']:
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Use histograms',
                'Use which histograms?',G2frame.dataWindow.HistsInPhase)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelections()
                    for Id,item in enumerate(G2frame.dataWindow.HistsInPhase):
                        if Id in sel:
                            data['Histograms'][item]['Use'] = True
                        else:
                            data['Histograms'][item]['Use'] = False
            finally:
                dlg.Destroy()
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)


    def OnDataCopy(event):
        hist = G2frame.hist
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        sourceDict = copy.deepcopy(data['Histograms'][hist])
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        copyDict = {}
        for name in copyNames:
            if name not in sourceDict: continue
            copyDict[name] = copy.deepcopy(sourceDict[name])        #force copy
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy phase/histogram parameters\nfrom '+hist[5:][:35],
                'Copy phase/hist parameters', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    data['Histograms'][keyList[sel]].update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnDataCopyFlags(event):
        hist = G2frame.hist
        sourceDict = copy.deepcopy(data['Histograms'][hist])
        copyDict = {}
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        babNames = ['BabA','BabU']
        for name in copyNames:
            if name not in sourceDict: continue
            if name in ['Scale','Extinction','HStrain','Flack','Twins','Layer Disp']:
                if name == 'Extinction' and 'HKLF' in sourceDict['Histogram']:
                    copyDict[name] = {name:[sourceDict[name][:2]]}
                    for item in ['Eg','Es','Ep']:
                        copyDict[name][item] = sourceDict[name][2][item][1]
                elif name == 'Twins':
                    copyDict[name] = sourceDict[name][0][1][1]
                else:
                    copyDict[name] = sourceDict[name][1]
            elif name in ['Size','Mustrain']:
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2],sourceDict[name][5]]
            elif name == 'Pref.Ori.':
                copyDict[name] = [sourceDict[name][0],sourceDict[name][2]]
                if sourceDict[name][0] == 'SH':
                    SHterms = sourceDict[name][5]
                    SHflags = {}
                    for item in SHterms:
                        SHflags[item] = SHterms[item]
                    copyDict[name].append(SHflags)
            elif name == 'Babinet':
                copyDict[name] = {}
                for bab in babNames:
                    copyDict[name][bab] = sourceDict[name][bab][1]
            elif name == 'Fix FXU' or name == 'FixedSeqVars':
                copyDict[name] = copy.deepcopy(sourceDict[name])
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy phase/histogram flags\nfrom '+hist[5:][:35],
                'Copy phase/hist flags', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    item = keyList[sel]
                    for name in copyNames:
                        if name not in sourceDict: continue
                        if name in ['Scale','Extinction','HStrain','Flack','Twins','Layer Disp']:
                            if name == 'Extinction' and 'HKLF' in sourceDict['Histogram']:
                                data['Histograms'][item][name][:2] = copy.deepcopy(sourceDict[name][:2])
                                for itm in ['Eg','Es','Ep']:
                                    data['Histograms'][item][name][2][itm][1] = copy.deepcopy(copyDict[name][itm])
                            elif name == 'Twins':
                                data['Histograms'][item]['Twins'][0][1][1] = copy.deepcopy(copyDict['Twins'])
                            else:
                                try:
                                    data['Histograms'][item][name][1] = copy.deepcopy(copyDict[name])
                                except KeyError:
                                    continue
                        elif name in ['Size','Mustrain']:
                            data['Histograms'][item][name][0] = copy.deepcopy(copyDict[name][0])
                            data['Histograms'][item][name][2] = copy.deepcopy(copyDict[name][1])
                            data['Histograms'][item][name][5] = copy.deepcopy(copyDict[name][2])
                        elif name == 'Pref.Ori.':
                            data['Histograms'][item][name][0] = copy.deepcopy(copyDict[name][0])
                            data['Histograms'][item][name][2] = copy.deepcopy(copyDict[name][1])
                            if sourceDict[name][0] == 'SH':
                               SHflags = copy.deepcopy(copyDict[name][2])
                               SHterms = copy.deepcopy(sourceDict[name][5])
                               data['Histograms'][item][name][6] = copy.deepcopy(sourceDict[name][6])
                               data['Histograms'][item][name][7] = copy.deepcopy(sourceDict[name][7])
                        elif name == 'Babinet':
                            for bab in babNames:
                                data['Histograms'][item][name][bab][1] = copy.deepcopy(copyDict[name][bab])
                        elif name == 'Fix FXU' or name == 'FixedSeqVars':
                            data['Histograms'][item][name] = copy.deepcopy(sourceDict[name])
        finally:
            dlg.Destroy()

    def OnSelDataCopy(event):
        '''Select HAP items to copy from one Phase/Hist to other(s)
        '''
        hist = G2frame.hist
        sourceDict = data['Histograms'][hist]
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select which parameters to copy',
            'Select phase data parameters', copyNames)
        selectedItems = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                selectedItems = [copyNames[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if not selectedItems: return # nothing to copy
        copyDict = {}
        for parm in selectedItems:
            if parm not in sourceDict: continue
            copyDict[parm] = copy.deepcopy(sourceDict[parm])
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy selected phase/histogram parameters\nfrom '+hist[5:][:35],
            'Copy selected phase/hist parameters', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    data['Histograms'][keyList[sel]].update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()

    def OnSelDataRead(event):
        '''Select HAP items to copy from another GPX file to current
        phase & hist
        '''
        sourceDict = data['Histograms'][G2frame.hist]
        try:
            pth = G2G.GetImportPath(G2frame)
            dlg = G2G.gpxFileSelector(parent=G2frame,startdir=pth)
            if wx.ID_OK == dlg.ShowModal():
                filename = dlg.Selection
            else:
                return
        finally:
            dlg.Destroy()

        import pickle
        phases = None
        phasenames = []
        try:
            fp = open(filename,'rb')
            while True:
                try:
                    d = pickle.load(fp)
                except EOFError:
                    break
                if d[0][0] == 'Phases':
                    phases = d
                    phasenames = [phases[i][0] for i in range(1,len(phases))]
        except:
            return
        finally:
            fp.close()
        if not phasenames: return
        if len(phasenames) == 1:
            phNum = 1
        else:
            dlg = wx.SingleChoiceDialog(G2frame,'Select Phase to use',
                                            'Select',phasenames)
            if dlg.ShowModal() == wx.ID_OK:
                phNum = dlg.GetSelection()+1
            else:
                return
        histograms = list(phases[phNum][1]['Histograms'].keys())
        if len(histograms) == 0:
            return
        elif len(histograms) == 1:
            histNam = histograms[0]
        else:
            dlg = wx.SingleChoiceDialog(G2frame,'Select histogram to use',
                                            'Select',histograms)
            if dlg.ShowModal() == wx.ID_OK:
                histNam = histograms[dlg.GetSelection()]
            else:
                return
        if 'HKLF' in histNam:
            copyNames = ['Extinction','Babinet','Flack','Twins']
        else:  #PWDR
            copyNames = ['Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Layer Disp']
        copyNames += ['Scale','Fix FXU','FixedSeqVars']
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select which parameters to copy',
            'Select phase data parameters', copyNames)
        selectedItems = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                selectedItems = [copyNames[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if not selectedItems: return # nothing to copy
        for i in selectedItems:
            sourceDict[i] = phases[phNum][1]['Histograms'][histNam][i]
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)

    def OnPwdrAdd(event):
        generalData = data['General']
        SGData = generalData['SGData']
        newList = []
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        keyList = data['Histograms'].keys()
        TextList = []
        if G2frame.GPXtree.GetCount():
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.GPXtree.GetItemText(item)
                if name not in keyList and 'PWDR' in name:
                    TextList.append(name)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
            if not TextList:
                G2G.G2MessageBox(G2frame,'No histograms')
                return
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select powder histograms to use',
                'Use data',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: newList.append(TextList[i])
                    if 'All PWDR' in newList:
                        newList = TextList[1:]
                    for histoName in newList:
                        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
                        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
                        data['Histograms'][histoName] = {
                            'Histogram':histoName,'Show':False,
                            'LeBail':False,
                            'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},['',],0.1],'Type':Inst['Type'][0],
                            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                                [1.,1.,1.,0.,0.,0.],6*[False,]],
                            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                            'HStrain':[NDij*[0.0,],NDij*[False,]],
                            'Layer Disp':[0.0,False],
                            'Extinction':[0.0,False],
                            'Flack':[0.0,False],
                            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
                            'Fix FXU':' ','FixedSeqVars':[]}
                        refList = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Reflection Lists'))
                        refList[generalData['Name']] = {}
                    wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
            finally:
                dlg.Destroy()

    def OnDataDelete(event):
        if G2frame.dataWindow.HistsInPhase:
            DelList = []
            extraOpts= {'label_0':'Remove from all phases','value_0':False}
            h,pd = G2frame.GetUsedHistogramsAndPhasesfromTree()
            if len(pd) > 1:
                opts = extraOpts
            else:
                opts = {}
            dlg = G2G.G2MultiChoiceDialog(G2frame,
                'Select histogram(s) to remove   \nfrom this phase:',
                'Remove histograms', G2frame.dataWindow.HistsInPhase,
                extraOpts=opts)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    DelList = [G2frame.dataWindow.HistsInPhase[i] for i in dlg.GetSelections()]
            finally:
                dlg.Destroy()
            if extraOpts['value_0']:
                for p in pd:
                    for i in DelList:
                        if i in pd[p]['Histograms']: del pd[p]['Histograms'][i]
            else:
                for i in DelList:
                    del data['Histograms'][i]
            #wx.CallLater(100,G2ddG.UpdateDData,G2frame,DData,data) #  produces error
            G2frame.GPXtree.UpdateSelection()
            # TId = G2frame.GPXtree.GetFocusedItem()
            # G2frame.GPXtree.SelectItem(G2frame.root)
            # G2frame.GPXtree.SelectItem(TId)
#            UpdatePhaseData(G2frame,Item,data)

    def OnDataApplyStrain(event):
        SGData = data['General']['SGData']
        DijVals = data['Histograms'][G2frame.hist]['HStrain'][0][:]
        # apply the Dij values to the reciprocal cell
        newA = []
        Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
        for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),['D11','D22','D33','D12','D13','D23']):
            newA.append(Aij + Dijdict.get(lbl,0.0))
        # convert back to direct cell
        data['General']['Cell'][1:7] = G2lat.A2cell(newA)
        data['General']['Cell'][7] = G2lat.calc_V(newA)
        # subtract the selected histograms Dij values from all for this phase
        for hist in data['Histograms']:
            for i,val in enumerate(DijVals):
                data['Histograms'][hist]['HStrain'][0][i] -= val
        # for hist in sorted(data['Histograms']): # list effective lattice constants applying Dij values
        #     DijVals = data['Histograms'][hist]['HStrain'][0]
        #     newA = []
        #     Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
        #     for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),
        #                     ['D11','D22','D33','D12','D13','D23']):
        #         newA.append(Aij + Dijdict.get(lbl,0.0))
        #     print(hist, G2lat.A2cell(newA)[:3], G2lat.calc_V(newA))
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)

#### Rigid bodies ################################################################################
    def FillRigidBodyGrid(refresh=True,vecId=None,resId=None,spnId=None):
        '''Fill the Rigid Body Phase information tab page.
        Note that the page is a ScrolledWindow, not a Grid
        '''
        resVarLookup = []
        def OnThermSel(event):
            Obj = event.GetEventObject()
            RBObj = Indx[Obj.GetId()]
            val = Obj.GetValue()
            Ttype = 'A'
            if val == 'Uiso':
                Ttype = 'I'
                RBObj['ThermalMotion'][0] = 'Uiso'
            elif val == 'T':
                RBObj['ThermalMotion'][0] = 'T'
            elif val == 'TL':
                RBObj['ThermalMotion'][0] = 'TL'
            elif val == 'TLS':
                RBObj['ThermalMotion'][0] = 'TLS'
            elif val == 'None':
                RBObj['ThermalMotion'][0] = 'None'
            if val != 'None':
                cia = data['General']['AtomPtrs'][3]
                for i,Id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[Id]][cia] = Ttype
            resId,vecId = None,None         #,spnId,None
            if resSelect is not None:
                resId = resSelect.GetSelection()
            if vecSelect is not None:
                vecId = vecSelect.GetSelection()
            # if spnSelect:
            #     spnId = spnSelect.GetSelection()
            wx.CallAfter(FillRigidBodyGrid,True,vecId=vecId,resId=resId)    ##,spnId=spnId
            G2plt.PlotStructure(G2frame,data)

        def ThermDataSizer(RBObj,rbType):

            def OnThermval(invalid,value,tc):
                Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[1]
                Uout = G2mth.UpdateRBUIJ(Bmat,Cart,RBObj)
                cia = data['General']['AtomPtrs'][3]
                for i,Id in enumerate(RBObj['Ids']):
                    if Uout[i][0] == 'I':
                        data['Atoms'][AtLookUp[Id]][cia+1] = Uout[i][1]
                    else:
                        data['Atoms'][AtLookUp[Id]][cia+2:cia+8] = Uout[i][2:8]
                G2plt.PlotStructure(G2frame,data)

            def OnTLSRef(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                RBObj['ThermalMotion'][2][item] = Obj.GetValue()

            thermSizer = wx.FlexGridSizer(0,9,5,5)
            model = RBObj['ThermalMotion']
            if model[0] == 'Uiso':
                names = ['Uiso',]
            elif 'T' in model[0]:
                names = ['T11','T22','T33','T12','T13','T23']
            if 'L' in model[0]:
                names += ['L11','L22','L33','L12','L13','L23']
            if 'S' in model[0]:
                names += ['S12','S13','S21','S23','S31','S32','SAA','SBB']
            for i,name in enumerate(names):
                thermSizer.Add(wx.StaticText(RigidBodies,-1,name+': '),0,WACV)
                thermVal = G2G.ValidatedTxtCtrl(RigidBodies,model[1],i,nDig=(8,4),OnLeave=OnThermval)
                thermSizer.Add(thermVal)
                Tcheck = wx.CheckBox(RigidBodies,-1,'Refine?')
                Tcheck.Bind(wx.EVT_CHECKBOX,OnTLSRef)
                Tcheck.SetValue(model[2][i])
                Indx[Tcheck.GetId()] = i
                thermSizer.Add(Tcheck,0,WACV)
            return thermSizer

        def LocationSizer(RBObj,rbType):

            def OnOrigRef(event):
                RBObj['Orig'][1] = Ocheck.GetValue()

            def OnOrienRef(event):
                RBObj['Orient'][1] = Qcheck.GetValue()

            def OnOrigX(invalid,value,tc):
                '''Called when the position info is changed (vector
                or azimuth)
                '''
                newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[0]
                Sytsym,Mult = G2spc.SytSym(RBObj['Orig'][0],SGData)[:2]
                sytsymtxt.SetLabel('Origin site symmetry: %s, multiplicity: %d '%(Sytsym,Mult))
                maxFrac = 0.0
                for Id in RBObj['Ids']:
                    maxFrac = max(maxFrac,data['Atoms'][AtLookUp[Id]][cx+3])
                for i,Id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[Id]][cx:cx+3] = newXYZ[i]
                    data['Atoms'][AtLookUp[Id]][cx+3] = maxFrac
                data['Atoms'] = G2lat.RBsymCheck(data['Atoms'],ct,cx,cs,AtLookUp,Amat,RBObj['Ids'],SGData)
                data['Drawing']['Atoms'] = []
                UpdateDrawAtoms(G2frame,data)
                G2plt.PlotStructure(G2frame,data)

            def OnOrien(*args, **kwargs):
                '''Called when the orientation info is changed (vector
                or azimuth). When called after a move of RB with alt key pressed,
                an optional keyword arg is provided with the mode to draw atoms
                (lines, ball & sticks,...)
                '''
                try:
                    orient = [float(Indx['Orien'][i].GetValue()) for i in range(4)]
                    A = orient[0]
                    V = np.inner(Amat,orient[1:]) # normalized in AVdeg2Q
                    Q = G2mth.AVdeg2Q(A,V)
                    if not any(Q):
                        raise ValueError
                    RBObj['Orient'][0] = Q
                    if rbType != 'Spin':
                        newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[0]
                        maxFrac = 0.0
                        for Id in RBObj['Ids']:
                            maxFrac = max(maxFrac,data['Atoms'][AtLookUp[Id]][cx+3])
                        for i,Id in enumerate(RBObj['Ids']):
                            data['Atoms'][AtLookUp[Id]][cx:cx+3] = newXYZ[i]
                            data['Atoms'][AtLookUp[Id]][cx+3] = maxFrac
                        data['Atoms'] = G2lat.RBsymCheck(data['Atoms'],ct,cx,cs,AtLookUp,Amat,RBObj['Ids'],SGData)
                    data['Drawing']['Atoms'] = []
                    if 'mode' in kwargs:
                        UpdateDrawAtoms(G2frame,data,kwargs['mode'])
                    else:
                        UpdateDrawAtoms(G2frame,data)
                    G2plt.PlotStructure(G2frame,data)
                except ValueError:
                    pass

            SGData = data['General']['SGData']
            rbSizer = wx.BoxSizer(wx.VERTICAL)
            topSizer = wx.FlexGridSizer(0,6,5,5)
            if rbType != 'Spin':
                if type(RBObj['Orig'][0]) is tuple:      # patch because somehow adding RB origin is becoming a tuple
                    if GSASIIpath.GetConfigValue('debug'): print('patching origin!')
                    RBObj['Orig'][0] = list(RBObj['Orig'][0])
                Sytsym,Mult = G2spc.SytSym(RBObj['Orig'][0],SGData)[:2]
                topSizer.Add(wx.StaticText(RigidBodies,-1,'Origin x,y,z (frac)'),0,WACV)
                topSizer.Add((-1,-1))
                Xsizers = []
                for ix in range(3):
                    origX = G2G.ValidatedTxtCtrl(RigidBodies,RBObj['Orig'][0],ix,nDig=(8,5),
                        typeHint=float,OnLeave=OnOrigX,xmin=-1,xmax=1.,size=(70,-1))
                    topSizer.Add(origX,0,WACV)
                    Xsizers.append(origX)
                G2frame.testRBObjSizers.update({'Xsizers':Xsizers})
                Ocheck = wx.CheckBox(RigidBodies,-1,'Refine?')
                Ocheck.Bind(wx.EVT_CHECKBOX,OnOrigRef)
                Ocheck.SetValue(RBObj['Orig'][1])
                # TODO: does spin RB need orientation vector? Does need angle & fix vector = [0,0,1]?
                topSizer.Add(Ocheck,0,WACV)
                Name = 'Origin'
                G2frame.testRBObjSizers['OnOrien'] = OnOrien
                G2frame.testRBObjSizers['FillUnitCell'] = FillUnitCell
            else:
                if 'Orig' in RBObj:     #cleanout - not using Orig for spinning RBs!
                    del RBObj['Orig']
                atId = RBObj['Ids'][0]
                Atom = data['Atoms'][AtLookUp[atId]]
                atXYZ = Atom[cx:cx+3]
                Sytsym,Mult = G2spc.SytSym(atXYZ,SGData)[:2]
                Name = Atom[ct-1]
            topSizer.Add(wx.StaticText(RigidBodies,-1,
                'Rotation angle (deg)\n&& Orient. vector (frac)'),0,WACV)
            Indx['Orien'] = {}
            Orien,OrienV = G2mth.Q2AVdeg(RBObj['Orient'][0])
            Orien = [Orien,]
            Orien.extend(np.inner(Bmat,OrienV)) # fractional coords
            dp,xmin,xmax = 2,-180.,360.
            OrientVecSiz = []
            for ix,x in enumerate(Orien):
                orien = G2G.ValidatedTxtCtrl(RigidBodies,Orien,ix,nDig=(8,dp),
                    typeHint=float,OnLeave=OnOrien,xmin=xmin,xmax=xmax,size=(70,-1))
                OrientVecSiz.append(orien)
                dp, xmin,xmax = 4,-1.,1.
                Indx['Orien'][ix] = orien
                topSizer.Add(orien,0,WACV)
            G2frame.testRBObjSizers.update({'OrientVecSiz':OrientVecSiz})
            Qchoice = [' ','A','AV','V']
            Qcheck = wx.ComboBox(RigidBodies,-1,value='',choices=Qchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Qcheck.Bind(wx.EVT_COMBOBOX,OnOrienRef)
            Qcheck.SetValue(RBObj['Orient'][1])
            topSizer.Add(Qcheck,0,WACV)
            RBObj['SytSym'] = Sytsym    #Only needed for spinning RBs
            sytsymtxt = wx.StaticText(RigidBodies,label='%s site symmetry: %s, multiplicity: %d '%(Name,Sytsym,Mult))
            rbSizer.Add(topSizer)
            rbSizer.Add(sytsymtxt)
            return rbSizer

        def SpnrbSizer(RBObj,spnIndx):
            '''Displays details for selected spinning rigid body'''

            def OnDelSpnRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBData['Spin'][RBId[0]]['useCount'] -= 1
                cia = data['General']['AtomPtrs'][3]
                atomData = data['Atoms']
                atomId = data['RBModels']['Spin'][spnIndx]['Ids'][0]
                for ia,atom in enumerate(atomData):
                    if atomId == atom[cia+8]:   #by definition only one
                        del atomData[ia]
                        break
                del data['RBModels']['Spin'][spnIndx]
                del data['General']['SpnIds'][atomId]
                data['Drawing']['Atoms'] = []
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True)

            def OnSymRadioSet(event):
                '''Set the polar axis for the sp. harm. as
                RBdata['Spin'][RBId]['symAxis']. This may never be
                set, so use RBdata['Spin'][RBId].get('symAxis') to
                access this so the default value is [0,0,1].
                '''
                Obj = event.GetEventObject()
                axis = ([1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1])[Obj.GetSelection()]
                RBObj['symAxis'] = axis
                G2plt.PlotStructure(G2frame,data)

            def OnAddShell(event):
                rbNames = []
                rbIds = {}
                for rbId in RBData['Spin']:
                    name = RBData['Spin'][rbId]['RBname']
                    rbNames.append(name)
                    rbIds[name] = rbId
                if len(rbNames) == 1:
                    selection = rbNames[0]
                else:
                    choices = sorted(rbNames)
                    dlg = G2G.G2SingleChoiceDialog(
                        G2frame,'Select rigid body to\nadd to structure',
                        'Select rigid body',choices)
                    dlg.CenterOnParent()
                    try:
                        if dlg.ShowModal() == wx.ID_OK:
                            selection = choices[dlg.GetSelection()]
                        else:
                            return
                    finally:
                        dlg.Destroy()

                rbData = RBData['Spin'][rbIds[selection]]
                RBObj['RBId'].append(rbIds[selection])
                RBObj['SHC'].append({})
                for name in ['atColor','atType','Natoms','nSH','RBname','RBsym']:
                    RBObj[name].append(rbData[name])
                RBObj['Radius'].append([1.0,False])
                RBObj['hide'].append(False)
                RBData['Spin'][rbIds[selection]]['useCount'] += 1
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True,spnId=rbId)

            def SHsizer():
                def OnSHOrder(event):
                    Obj = event.GetEventObject()
                    iSh = Indx[Obj.GetId()]
                    RBObj['nSH'][iSh] = int(Obj.GetValue())
                    RBObj['SHC'][iSh] = SetSHCoef(iSh,RBObj['nSH'][iSh])
                    G2plt.PlotStructure(G2frame,data)
                    wx.CallAfter(FillRigidBodyGrid,True,spnId=rbId)

                def SetSHCoef(iSh,Order):
                    Sytsym = RBObj['SytSym']
                    cofNames,cofSgns = G2lat.GenRBCoeff(Sytsym,RBObj['RBsym'][iSh],Order)
                    cofTerms = [[0.0,val,False] for val in cofSgns]
                    newSHcoef = dict(zip(cofNames,cofTerms))
                    SHcoef = RBObj['SHC'][iSh]
                    for cofName in SHcoef:      #transfer old values to new set
                        if cofName in newSHcoef:
                            newSHcoef[cofName] = SHcoef[cofName]
                    return newSHcoef

                def OnSchRef(event):
                    Obj = event.GetEventObject()
                    iSh,name = Indx[Obj.GetId()]
                    RBObj['SHC'][iSh][name][2] = not RBObj['SHC'][iSh][name][2]

                def OnRadRef(event):
                    Obj = event.GetEventObject()
                    iSh = Indx[Obj.GetId()]
                    RBObj['Radius'][iSh][1] = not RBObj['Radius'][iSh][1]

                def NewSHC(invalid,value,tc):
                    G2plt.PlotStructure(G2frame,data)

                def OnDelShell(event):
                    Obj = event.GetEventObject()
                    iSh = Indx[Obj.GetId()]
                    rbId = RBObj['RBId'][iSh]
                    RBData['Spin'][rbId]['useCount'] -= 1
                    RBData['Spin'][rbId]['useCount'] = max(0,RBData['Spin'][rbId]['useCount'])
                    for name in ['atColor','atType','Natoms','nSH','RBId','RBname','RBsym','SHC','Radius']:
                        del RBObj[name][iSh]
                    G2plt.PlotStructure(G2frame,data)
                    wx.CallAfter(FillRigidBodyGrid,True,spnId=rbId)

                shSizer = wx.BoxSizer(wx.VERTICAL)
                for iSh,nSh in enumerate(RBObj['nSH']):
                    #patch
                    if 'Radius' not in RBObj:
                        RBObj['Radius'] = [[1.0,False] for i in range(len(RBObj['nSH']))]
                    #end patch
                    rbId = RBObj['RBId'][iSh]
                    RBObj['atType'][iSh] = RBData['Spin'][rbId]['atType']
                    RBObj['atColor'][iSh] = G2elem.GetAtomInfo(RBObj['atType'][iSh])['Color']  #correct atom color for shell
                    if iSh:
                        subLine = wx.BoxSizer(wx.HORIZONTAL)
                        subLine.Add(wx.StaticText(RigidBodies,label='Shell %d: Name: %s   Atom type: %s RB sym: %s '  \
                            %(iSh,RBObj['RBname'][iSh],RBObj['atType'][iSh],RBObj['RBsym'][iSh])),0,WACV)
                        delShell = wx.Button(RigidBodies,label='Delete shell',style=wx.BU_EXACTFIT)
                        Indx[delShell.GetId()] = iSh
                        delShell.Bind(wx.EVT_BUTTON,OnDelShell)
                        subLine.Add(delShell,0,WACV)
                        hidesh = wx.CheckBox(RigidBodies,label='Hide shell?')
                        hidesh.SetValue(RBObj['hide'][iSh])
                        hidesh.Bind(wx.EVT_CHECKBOX,OnHideSh)
                        Indx[hidesh.GetId()] = iSh
                        subLine.Add(hidesh,0,WACV)
                        shSizer.Add(subLine)
                    shoSizer = wx.BoxSizer(wx.HORIZONTAL)
                    shoSizer.Add(wx.StaticText(RigidBodies,label=' Bessel/Harmonic order: '),0,WACV)
                    shOrder = wx.ComboBox(RigidBodies,value=str(RBObj['nSH'][iSh]),choices=[str(i) for i in range(19)],
                        style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    Indx[shOrder.GetId()] = iSh
                    shOrder.Bind(wx.EVT_COMBOBOX,OnSHOrder)
                    shoSizer.Add(shOrder,0,WACV)
                    if RBObj['nSH'][iSh]>0:
                        shoSizer.Add(wx.StaticText(RigidBodies,label=" 'c' for cubic harmonic term. "),0,WACV)
                    shoSizer.Add(wx.StaticText(RigidBodies,label=' Radius: '),0,WACV)
                    shoSizer.Add(G2G.ValidatedTxtCtrl(RigidBodies,RBObj['Radius'][iSh],0,nDig=(8,5),xmin=0.0,xmax=5.0,
                        typeHint=float,size=(70,-1),OnLeave=NewSHC),0,WACV)
                    if 'Q' in RBObj['atType'][iSh]:
                        radref = wx.StaticText(RigidBodies,label=' For drawing only')
                    else:
                        radref = wx.CheckBox(RigidBodies,label=' refine? ')
                        radref.SetValue(RBObj['Radius'][iSh][1])
                        radref.Bind(wx.EVT_CHECKBOX,OnRadRef)
                        Indx[radref.GetId()] = iSh
                    shoSizer.Add(radref,0,WACV)
                    shSizer.Add(shoSizer)
                    if not RBObj['nSH'][iSh]:
                        shSizer.Add(wx.StaticText(RigidBodies,
                            label=' Select harmonic order or try different equivalent position'))
                    elif len(RBObj['SHC'][iSh]) > 24:
                        shSizer.Add(wx.StaticText(RigidBodies,
                            label=' WARNING: More than 24 terms found; use lower harmonic order'))
                    else:
                        shcSizer = wx.FlexGridSizer(0,9,5,5)
                        for item in RBObj['SHC'][iSh]:
                            shcSizer.Add(wx.StaticText(RigidBodies,label='%s'%item.strip('+').strip('-')))
                            shcSizer.Add(G2G.ValidatedTxtCtrl(RigidBodies,RBObj['SHC'][iSh][item],0,nDig=(8,5),
                                typeHint=float,size=(70,-1),OnLeave=NewSHC))
                            schref = wx.CheckBox(RigidBodies,label=' refine? ')
                            schref.SetValue(RBObj['SHC'][iSh][item][2])
                            schref.Bind(wx.EVT_CHECKBOX,OnSchRef)
                            Indx[schref.GetId()] = iSh,item
                            shcSizer.Add(schref)
                        shSizer.Add(shcSizer)
                return shSizer

            def OnHideSh(event):
                Obj = event.GetEventObject()
                iSh = Indx[Obj.GetId()]
                RBObj['hide'][iSh] = not RBObj['hide'][iSh]
                G2plt.PlotStructure(G2frame,data)
                
            def OnAtColor(event):
                RBObj['useAtColor'] = not RBObj['useAtColor']
                G2plt.PlotStructure(G2frame,data)
                
            def OnFadeShell(event):
                RBObj['fadeSh'] = not RBObj['fadeSh']
                G2plt.PlotStructure(G2frame,data)

            RBObj['hide'] = RBObj.get('hide',[False for i in range(len(RBObj['atType']))])
            rbId = RBObj['RBId'][0]
            atId = RBObj['Ids'][0]
            atName = data['Atoms'][AtLookUp[atId]][ct-1]
            RBObj['atType'][0] = RBData['Spin'][rbId]['atType']
            sprbSizer = wx.BoxSizer(wx.VERTICAL)
            G2G.HorizontalLine(sprbSizer,RigidBodies)
            topLine = wx.BoxSizer(wx.HORIZONTAL)
            topLine.Add(wx.StaticText(RigidBodies,label='Shell 0: Name: %s Atom name: %s Atom type: %s RB sym: %s '%
                (RBObj['RBname'][0],atName,RBObj['atType'][0],RBObj['RBsym'][0])),0,WACV)
            rbId = RBObj['RBId']
            if len(RBObj['nSH']) == 1:
                delRB = wx.Button(RigidBodies,wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
                delRB.Bind(wx.EVT_BUTTON,OnDelSpnRB)
                Indx[delRB.GetId()] = rbId
                topLine.Add(delRB,0,WACV)
            addShell = wx.Button(RigidBodies,wx.ID_ANY,'Add new shell',style=wx.BU_EXACTFIT)
            addShell.Bind(wx.EVT_BUTTON,OnAddShell)
            Indx[addShell.GetId()] = rbId
            topLine.Add(addShell,0,WACV)
            hidesh = wx.CheckBox(RigidBodies,label='Hide shell?')
            hidesh.SetValue(RBObj['hide'][0])
            hidesh.Bind(wx.EVT_CHECKBOX,OnHideSh)
            Indx[hidesh.GetId()] = 0
            topLine.Add(hidesh,0,WACV)
            sprbSizer.Add(wx.StaticText(RigidBodies,label='Spinning RB orientation parameters for %s:'%RBObj['RBname'][0]))
            sprbSizer.Add(LocationSizer(RBObj,'Spin'))
            choices = [' x ',' y ',' z ','x+y','x+y+z']
            RBObj['symAxis'] = RBObj.get('symAxis',[0,0,1])   #set default as 'z'
            symax = dict(zip([str(x) for x in [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1]]],choices))[str(RBObj['symAxis'])]
            symRadioSet = wx.RadioBox(RigidBodies,choices=choices,label='Sp harm polar axis is aligned along:')
            symRadioSet.SetStringSelection(symax)
            symRadioSet.Bind(wx.EVT_RADIOBOX, OnSymRadioSet)
            Indx[symRadioSet.GetId()] = rbId
            sprbSizer.Add(symRadioSet)
            plotLine = wx.BoxSizer(wx.HORIZONTAL)
            RBObj['useAtColor'] = RBObj.get('useAtColor',True)
            atColor = wx.CheckBox(RigidBodies,label='Use atom color?')
            atColor.SetValue(RBObj['useAtColor'])
            atColor.Bind(wx.EVT_CHECKBOX,OnAtColor)
            plotLine.Add(atColor,0,WACV)
            RBObj['fadeSh'] = RBObj.get('fadeSh',True)
            fadeShell = wx.CheckBox(RigidBodies,label='Fade shells?')
            fadeShell.SetValue(RBObj['fadeSh'])
            fadeShell.Bind(wx.EVT_CHECKBOX,OnFadeShell)
            plotLine.Add(fadeShell,0,WACV)
            sprbSizer.Add(plotLine)
            G2G.HorizontalLine(sprbSizer,RigidBodies)
            sprbSizer.Add(topLine)
            sprbSizer.Add(SHsizer())
            return sprbSizer

        def ResrbSizer(RBObj,resIndx):
            '''Displays details for selected residue rigid body'''
            def OnTorsionRef(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                RBObj['Torsions'][item][1] = Obj.GetValue()

            def OnTorsion(invalid,value,tc):
                newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')[0]
                for i,Id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[Id]][cx:cx+3] = newXYZ[i]
                data['Drawing']['Atoms'] = []
                UpdateDrawAtoms(G2frame,data,atomStyle)
                G2frame.drawAtoms.ClearSelection()
                G2plt.PlotStructure(G2frame,data)

            def OnFrac(invalid,value,tc):
                for i,Id in enumerate(RBObj['Ids']):
                    if data['Atoms'][AtLookUp[Id]][cx+3]:
                        data['Atoms'][AtLookUp[Id]][cx+3] = value

            def OnRefFrac(event):
                RBObj['AtomFrac'][1] = not RBObj['AtomFrac'][1]

            def OnDelResRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBData['Residue'][RBId]['useCount'] -= 1
                del data['RBModels']['Residue'][resIndx]
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True,resId=resIndx)

            resrbSizer = wx.BoxSizer(wx.VERTICAL)
            G2G.HorizontalLine(resrbSizer,RigidBodies)
            topLine = wx.BoxSizer(wx.HORIZONTAL)
            topLine.Add(wx.StaticText(RigidBodies,-1,'Name: '+RBObj['RBname']+RBObj['numChain']+'   '),0,WACV)
            rbId = RBObj['RBId']
            delRB = wx.Button(RigidBodies,wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
            delRB.Bind(wx.EVT_BUTTON,OnDelResRB)
            Indx[delRB.GetId()] = rbId
            topLine.Add(delRB,0,WACV)
            symAxis = RBObj.get('symAxis')
            if np.any(symAxis):
                if np.all(symAxis):
                    lbl = 'x+y+z'
                elif np.all(symAxis[:2]):
                    lbl = 'x+y'
                elif symAxis[0]:
                    lbl = 'x'
                elif symAxis[1]:
                    lbl = 'y'
                else:
                    lbl = 'z'
                topLine.Add(wx.StaticText(RigidBodies,-1,
                    '   Rigid body {} axis is aligned along oriention vector'.format(lbl)),0,WACV)
            try:
                varname = str(data['pId'])+'::RBRxxx:'+resVarLookup[resIndx]
            except:  # happens when phase has no histograms
                varname = '?::RBRxxx:'+resVarLookup[resIndx]
            topLine.Add(wx.StaticText(RigidBodies,-1,
                    '  (variables '+varname+')'),0,WACV)
            resrbSizer.Add(topLine)
            resrbSizer.Add(LocationSizer(RBObj,'Residue'))
            if len(RBObj['Torsions']):
                resrbSizer.Add(wx.StaticText(RigidBodies,-1,'Torsions:'),0)
            torSizer = wx.FlexGridSizer(0,6,5,5)
            for itors,tors in enumerate(RBObj['Torsions']):
                torSizer.Add(wx.StaticText(RigidBodies,-1,'Torsion '+'%d'%(itors)),0,WACV)
                torsTxt = G2G.ValidatedTxtCtrl(RigidBodies,RBObj['Torsions'][itors],0,nDig=(10,3),OnLeave=OnTorsion)
                torSizer.Add(torsTxt)
                torCheck = wx.CheckBox(RigidBodies,-1,'Refine?')
                torCheck.Bind(wx.EVT_CHECKBOX,OnTorsionRef)
                torCheck.SetValue(tors[1])
                Indx[torCheck.GetId()] = itors
                torSizer.Add(torCheck,0,WACV)
            resrbSizer.Add(torSizer)
            members = 'Rigid body members: '
            for nId,Id in enumerate(RBObj['Ids']):
                if nId and not nId%10:
                    members += ('\n'+30*' ')
                members += data['Atoms'][AtLookUp[Id]][ct-1].strip()+', '
            resrbSizer.Add(wx.StaticText(RigidBodies,label=members[:-2]),0)
            fracSizer = wx.BoxSizer(wx.HORIZONTAL)
            fracSizer.Add(wx.StaticText(RigidBodies,label='Rigid Body atom site fraction: '))
            fracTxt = G2G.ValidatedTxtCtrl(RigidBodies,RBObj['AtomFrac'],0,nDig=(10,3),OnLeave=OnFrac)
            fracSizer.Add(fracTxt,0,WACV)
            fracRef = wx.CheckBox(RigidBodies,label='Refine?')
            fracRef.SetValue(RBObj['AtomFrac'][1])
            fracRef.Bind(wx.EVT_CHECKBOX,OnRefFrac)
            fracSizer.Add(fracRef,0,WACV)
            resrbSizer.Add(fracSizer)
            tchoice = ['None','Uiso','T','TL','TLS']
            thermSizer = wx.BoxSizer(wx.HORIZONTAL)
            thermSizer.Add(wx.StaticText(RigidBodies,-1,'Rigid body thermal motion model: '),0,WACV)
            thermSel = wx.ComboBox(RigidBodies,-1,value=RBObj['ThermalMotion'][0],choices=tchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[thermSel.GetId()] = RBObj
            thermSel.Bind(wx.EVT_COMBOBOX,OnThermSel)
            thermSizer.Add(thermSel,0,WACV)
            thermSizer.Add(wx.StaticText(RigidBodies,-1,' Units: T A^2, L deg^2, S deg-A'),0,WACV)
            resrbSizer.Add(thermSizer)
            if RBObj['ThermalMotion'][0] != 'None':
                resrbSizer.Add(ThermDataSizer(RBObj,'Residue'))
            dragSizer = wx.BoxSizer(wx.HORIZONTAL)
            dragSizer.Add(wx.StaticText(RigidBodies,wx.ID_ANY,'Draw mode after dragging the rigid body: '),0,WACV)
            RBObj['drawMode'] = RBObj.get('drawMode',DrawStyleChoice[4])
            modeOpt = G2G.G2ChoiceButton(RigidBodies, DrawStyleChoice[1:],
                                             strLoc=RBObj, strKey='drawMode')
            dragSizer.Add(modeOpt)
            modeOpt.Enable(False) # not implemented yet
            G2frame.testRBObjSizers['fillMode'] = G2frame.testRBObjSizers.get('fillMode',False)
            dragSizer.Add(G2G.G2CheckBoxFrontLbl(RigidBodies, ' Fill cell on mouse up?', G2frame.testRBObjSizers, 'fillMode'))
            resrbSizer.Add((-1,5))
            resrbSizer.Add(dragSizer)
            return resrbSizer

        def VecrbSizer(RBObj,resIndx):
            '''Displays details for selected vector rigid body'''
            def OnFrac(invalid,value,tc):
                for Id in RBObj['Ids']:
                    if data['Atoms'][AtLookUp[Id]][cx+3]:
                        data['Atoms'][AtLookUp[Id]][cx+3] = value

            def OnRefFrac(event):
                RBObj['AtomFrac'][1] = not RBObj['AtomFrac'][1]

            def OnDelVecRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBData['Vector'][RBId]['useCount'] -= 1
                del data['RBModels']['Vector'][resIndx]
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True,vecId=resIndx)

            vecrbSizer = wx.BoxSizer(wx.VERTICAL)
            G2G.HorizontalLine(vecrbSizer,RigidBodies)
            topLine = wx.BoxSizer(wx.HORIZONTAL)
            topLine.Add(wx.StaticText(RigidBodies,-1,
                'Name: '+RBObj['RBname']+'   '),0,WACV)
            rbId = RBObj['RBId']
            delRB = wx.Button(RigidBodies,wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
            delRB.Bind(wx.EVT_BUTTON,OnDelVecRB)
            Indx[delRB.GetId()] = rbId
            topLine.Add(delRB,0,WACV)
            vecrbSizer.Add(topLine)
            vecrbSizer.Add(LocationSizer(RBObj,'Vector'))
            members = 'Rigid body members: '
            for Id in RBObj['Ids']:
                members += data['Atoms'][AtLookUp[Id]][ct-1].strip()+', '
            vecrbSizer.Add(wx.StaticText(RigidBodies,label=members[:-2]))
            fracSizer = wx.BoxSizer(wx.HORIZONTAL)
            fracSizer.Add(wx.StaticText(RigidBodies,label='Rigid Body atom site fraction: '))
            fracTxt = G2G.ValidatedTxtCtrl(RigidBodies,RBObj['AtomFrac'],0,nDig=(10,3),OnLeave=OnFrac)
            fracSizer.Add(fracTxt,0,WACV)
            fracRef = wx.CheckBox(RigidBodies,label='Refine?')
            fracRef.SetValue(RBObj['AtomFrac'][1])
            fracRef.Bind(wx.EVT_CHECKBOX,OnRefFrac)
            fracSizer.Add(fracRef,0,WACV)
            vecrbSizer.Add(fracSizer)
            tchoice = ['None','Uiso','T','TL','TLS']
            thermSizer = wx.BoxSizer(wx.HORIZONTAL)
            thermSizer.Add(wx.StaticText(RigidBodies,-1,'Rigid body thermal motion model: '),0,WACV)
            thermSel = wx.ComboBox(RigidBodies,-1,value=RBObj['ThermalMotion'][0],choices=tchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[thermSel.GetId()] = RBObj
            thermSel.Bind(wx.EVT_COMBOBOX,OnThermSel)
            thermSizer.Add(thermSel,0,WACV)
            thermSizer.Add(wx.StaticText(RigidBodies,-1,' Units: T A^2, L deg^2, S deg-A'),0,WACV)
            vecrbSizer.Add(thermSizer)
            if RBObj['ThermalMotion'][0] != 'None':
                vecrbSizer.Add(ThermDataSizer(RBObj,'Vector'))
            return vecrbSizer

        def OnVecSelect(event):
            global prevVecId
            prevVecId = vecSelect.GetSelection()
            try:
                resSelect.Deselect(resSelect.GetSelection())
            except:
                pass
            try:
                spnSelect.Deselect(spnSelect.GetSelection())
            except:
                pass
            wx.CallLater(100,RepaintRBInfo,'Vector',prevVecId)

        def OnResSelect(event):
            global prevResId
            prevResId = resSelect.GetSelection()
            try:
                vecSelect.Deselect(vecSelect.GetSelection())
            except:
                pass
            try:
                spnSelect.Deselect(spnSelect.GetSelection())
            except:
                pass
            # define the parameters needed to drag the RB with the mouse
            data['testRBObj'] = {}
            rbType = 'Residue'
            data['testRBObj']['rbObj'] = copy.deepcopy(data['RBModels'][rbType][prevResId])
            rbId = data['RBModels'][rbType][prevResId]['RBId']
            RBdata = G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
            data['testRBObj']['rbData'] = RBdata
            data['testRBObj']['rbType'] = rbType
            data['testRBObj']['rbAtTypes'] = RBdata[rbType][rbId]['rbTypes']
            data['testRBObj']['AtInfo'] = RBData[rbType]['AtInfo']
            data['testRBObj']['NameLookup'] = RBData[rbType][rbId].get('atNames',[])    #only for residues
            data['testRBObj']['Sizers'] = {}
            data['testRBObj']['rbRef'] = RBData[rbType][rbId]['rbRef']

            refType = []
            for ref in data['testRBObj']['rbRef'][:3]:
                reftype = data['testRBObj']['rbAtTypes'][ref]
                refType.append(reftype)
                #refName.append(reftype+' '+str(rbRef[0]))
            atNames = [{},{},{}]
            AtNames = {}
            cx,ct,cs,cia = data['General']['AtomPtrs']
            for iatm,atom in enumerate(data['Atoms']):
                AtNames[atom[ct-1]] = iatm
                for i,reftype in enumerate(refType):
                    if atom[ct] == reftype:
                        atNames[i][atom[ct-1]] = iatm
            data['testRBObj']['atNames'] = atNames
            data['testRBObj']['AtNames'] = AtNames
            data['testRBObj']['torAtms'] = []
            # unclear why these torsion entries are being added to rbObj.
            for item in RBData[rbType][rbId].get('rbSeq',[]):
                data['testRBObj']['rbObj']['Torsions'].append([item[2],False])  # Needed?
                data['testRBObj']['torAtms'].append([-1,-1,-1])
            wx.CallLater(100,RepaintRBInfo,'Residue',prevResId)

        def OnSpnSelect(event):
            global prevSpnId
            prevSpnId = spnSelect.GetSelection()
            try:
                resSelect.Deselect(resSelect.GetSelection())
            except:
                pass
            try:
                vecSelect.Deselect(vecSelect.GetSelection())
            except:
                pass
            wx.CallLater(100,RepaintRBInfo,'Spin',prevSpnId)

        def RepaintRBInfo(rbType,rbIndx,Scroll=0):
            oldFocus = wx.Window.FindFocus()
            try:
                if 'phoenix' in wx.version():
                    G2frame.bottomSizer.Clear(True)
                else:
                    G2frame.bottomSizer.DeleteWindows()
            except:
                return
            Indx.clear()
            rbObj = data['RBModels'][rbType][rbIndx]
            Quad = rbObj['Orient'][0]
            data['Drawing']['Quaternion'] = G2mth.invQ(Quad)
            if rbType == 'Residue':
                data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
                G2frame.bottomSizer =  ResrbSizer(rbObj,rbIndx)
            elif rbType == 'Spin':
                data['Drawing']['viewPoint'][0] = data['Atoms'][AtLookUp[rbObj['Ids'][0]]][cx:cx+3]
                G2frame.bottomSizer =  SpnrbSizer(rbObj,rbIndx)
            else: #Vector
                data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
                G2frame.bottomSizer =  VecrbSizer(rbObj,rbIndx)
            mainSizer.Add(G2frame.bottomSizer)
            mainSizer.Layout()
            G2frame.dataWindow.Refresh()
            RigidBodies.SetVirtualSize(mainSizer.GetMinSize())
            RigidBodies.Scroll(0,Scroll)
            G2frame.dataWindow.SendSizeEvent()
            G2plt.PlotStructure(G2frame,data)
            if oldFocus: wx.CallAfter(oldFocus.SetFocus)

        # FillRigidBodyGrid executable code starts here
        if refresh:
            if RigidBodies.GetSizer(): RigidBodies.GetSizer().Clear(True)
        if 'testRBObj' in data: del data['testRBObj']
        general = data['General']
        cx,ct,cs,cia = general['AtomPtrs']
        AtLookUp = G2mth.FillAtomLookUp(data['Atoms'],cia+8)
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies')
        if not Id:
            G2frame.CheckNotebook()
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies')
        if not Id:
            print('Strange! Why no Rigid Bodies tree entry?')
            return
        RBData = G2frame.GPXtree.GetItemPyData(Id)
        Indx = {}
        atomStyle = 'balls & sticks'
        if 'macro' in general['Type']:
            atomStyle = 'sticks'
        G2frame.GetStatusBar().SetStatusText('',1)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Select rigid body to view for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)
        nobody = True
        resSelect = None
        vecSelect = None
        spnSelect = None
        rbSizer = wx.BoxSizer(wx.HORIZONTAL)
        if 'Residue' in data['RBModels'] and len(data['RBModels']['Residue']):
            nobody = False
            resSizer = wx.BoxSizer(wx.VERTICAL)
            resSizer.Add(wx.StaticText(RigidBodies,label='Residue RBs:'))
            RBnames = []
            resVarLookup = []
            for irb,RBObj in enumerate(data['RBModels']['Residue']):
                # patch
                if 'AtomFrac' not in RBObj:
                    RBObj['AtomFrac'] = [1.0,False]
                #end patch
                name = RBObj['RBname']+RBObj['numChain']
                RBnames.append(name)
                jrb = RBData['RBIds']['Residue'].index(RBObj['RBId'])
                resVarLookup.append(str(irb)+':'+str(jrb)) # var name suffix following G2strMath.ApplyRBModels()
            if prevResId is not None:
                resId = prevResId
            try:
                RBnames[resId]
            except:
                resId = 0
                #rbName = RBnames[resId]
            rbObj = data['RBModels']['Residue'][resId]
            data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
            data['Drawing']['Quaternion'] = rbObj['Orient'][0]
            resSelect = wx.ListBox(RigidBodies,choices=RBnames,style=wx.LB_SINGLE,size=(-1,80))
            if resId:
                resSelect.SetSelection(resId)
                OnResSelect(None)
            resSelect.Bind(wx.EVT_LISTBOX,OnResSelect)
            resSizer.Add(resSelect)
            rbSizer.Add(resSizer)
        if 'Vector' in data['RBModels'] and len(data['RBModels']['Vector']):
            vecSizer = wx.BoxSizer(wx.VERTICAL)
            vecSizer.Add(wx.StaticText(RigidBodies,label='Vector RBs:'))
            nobody = False
            RBnames = []
            for RBObj in data['RBModels']['Vector']:
                # patch
                if 'AtomFrac' not in RBObj:
                    RBObj['AtomFrac'] = [1.0,False]
                #end patch
                RBnames.append(RBObj['RBname'])
            vecId = -1
            if prevVecId is not None:
                vecId = prevVecId
            try:
                RBnames[vecId]
            except:
                vecId = 0
                #rbName = RBnames[vecId]
            rbObj = data['RBModels']['Vector'][vecId]
            data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
            data['Drawing']['Quaternion'] = rbObj['Orient'][0]
            vecSelect = wx.ListBox(RigidBodies,choices=RBnames,style=wx.LB_SINGLE,size=(-1,80))
            if vecId is not None:
                vecSelect.SetSelection(vecId)
                OnVecSelect(None)
            vecSelect.Bind(wx.EVT_LISTBOX,OnVecSelect)
            vecSizer.Add(vecSelect)
            rbSizer.Add(vecSizer)
        if 'Spin' in data['RBModels'] and len(data['RBModels']['Spin']):
            spnSizer = wx.BoxSizer(wx.VERTICAL)
            spnSizer.Add(wx.StaticText(RigidBodies,label='Spinning RBs:'))
            nobody = False
            RBnames = []
            for RBObj in data['RBModels']['Spin']:
                RBnames.append(RBObj['RBname'][0])
            spnId = -1
            if prevSpnId is not None:
                spnId = prevSpnId
                try:
                    RBnames[spnId]
                except:
                    spnId = 0
                    #rbName = RBnames[spnId]
            rbObj = data['RBModels']['Spin'][spnId]
            data['Drawing']['viewPoint'][0] = data['Atoms'][AtLookUp[RBObj['Ids'][0]]][cx:cx+3]
            data['Drawing']['Quaternion'] = rbObj['Orient'][0]
            spnSelect = wx.ListBox(RigidBodies,choices=RBnames,style=wx.LB_SINGLE,size=(-1,80))
            if spnId != -1:
                spnSelect.SetSelection(spnId)
                OnSpnSelect(None)
            spnSelect.Bind(wx.EVT_LISTBOX,OnSpnSelect)
            spnSizer.Add(spnSelect,0)
            rbSizer.Add(spnSizer)
        mainSizer.Add(rbSizer,0,wx.EXPAND)
        G2frame.bottomSizer = wx.BoxSizer(wx.VERTICAL)
        G2frame.bottomSizer.Add(wx.StaticText(RigidBodies,label=' '))
        mainSizer.Add(G2frame.bottomSizer)
        G2plt.PlotStructure(G2frame,data)
        G2plt.PlotStructure(G2frame,data) # draw twice initially for mac
        if nobody:
            msg = 'Define a rigid body with the "Rigid Bodies" tree entry before adding it to the phase here'
            if RBData.get('RBIds') is None:
                pass
            elif len(RBData['RBIds'].get('Vector',[])) + len(RBData['RBIds'].get('Residue',[])) == 0:
                pass
            else:
                msg = 'No rigid bodies defined in phase. Use "Edit Body"/"Locate & Insert..."\ncommand to add them.'
            mainSizer.Add(wx.StaticText(RigidBodies,label=msg),0)
        SetPhaseWindow(RigidBodies,mainSizer)

    def OnRBCopyParms(event):
        RBObjs = []
        for rbType in ['Vector','Residue']:
            RBObjs += data['RBModels'].get(rbType,[])
        if not len(RBObjs):
            print ('**** ERROR - no rigid bodies defined ****')
            return
        if len(RBObjs) == 1:
            print ('**** INFO - only one rigid body defined; nothing to copy to ****')
            return
        Source = []
        sourceRB = {}
        for RBObj in RBObjs:
            Source.append(RBObj['RBname'])
        dlg = wx.SingleChoiceDialog(G2frame,'Select source','Copy rigid body parameters',Source)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            for item in ['Orig','Orient','ThermalMotion','AtomFract']:
                sourceRB.update({item:RBObjs[sel][item],})
        dlg.Destroy()
        if not sourceRB:
            return
        dlg = wx.MultiChoiceDialog(G2frame,'Select targets','Copy rigid body parameters',Source)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                RBObjs[x].update(copy.copy(sourceRB))
        G2plt.PlotStructure(G2frame,data)
        wx.CallAfter(FillRigidBodyGrid,True)

    def assignAtoms(RBData,selDict={},unmatchedRBatoms=None):
        '''Find the closest RB atoms to atoms in the structure
        If selDict is specified, it overrides the assignments to specify
        atoms that should be matched.
        '''
        general = data['General']
        cx,ct = general['AtomPtrs'][:2]
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])

        rbType = data['testRBObj']['rbType']
        rbObj = data['testRBObj']['rbObj']
        rbId = rbObj['RBId']
        matchTable = []
        if rbType == 'Spin':
            rbData = RBData['Spin'][rbId]
            rbData['Orig'] = [0.,0.,0.]
            matchTable.append([rbData['RBname'],rbData['rbType']]+list(rbData['Orig']))
            return matchTable
        rbAtmTypes = RBData[rbType][rbId]['rbTypes']
        atomData = data['Atoms']
        if 'atNames' in RBData[rbType][rbId]:
            rbAtmLbs = RBData[rbType][rbId]['atNames']
        else:
            rbAtmLbs = None

        newXYZ = G2mth.UpdateRBXYZ(Bmat,rbObj,RBData,rbType)[0]
        rbUsedIds = []   # Ids of atoms in current phase used inside RBs
        for i in data['RBModels']:
            for j in data['RBModels'][i]:
                rbUsedIds += j['Ids']
        # categorize atoms by type, omitting any that are already assigned
        # in a rigid body
        atmTypes = [None if atomData[i][-1] in rbUsedIds else atomData[i][ct] for i in range(len(atomData))]
        # remove assigned atoms from search groups
        for i in selDict:
            if selDict[i] is None: continue
            atmTypes[selDict[i]] = None
        atmXYZ = G2mth.getAtomXYZ(atomData,cx)
        # separate structure's atoms by type (w/o assigned atoms)
        oXYZbyT = {}
        atmNumByT = {}
        for t in set(atmTypes):
            if t is None: continue
            oXYZbyT[t] = np.array([atmXYZ[i] for i in range(len(atmXYZ)) if atmTypes[i] == t])
            atmNumByT[t] = [i for i in range(len(atmXYZ)) if atmTypes[i] == t]
        Ids = []
        # create table of fixed and found assignments
        unmatched = []
        for i,xyz in enumerate(newXYZ):
            t = rbAtmTypes[i]
            if rbAtmLbs:
                lbl = rbAtmLbs[i]
            else:
                lbl = ''
            if i in selDict and selDict[i] is None:
                matchTable.append([t , lbl] + list(xyz))
                continue
            elif i in selDict:
                searchXYZ = [atmXYZ[selDict[i]]] #assigned
                numLookup = [selDict[i]]
            else:
                if t not in oXYZbyT:
                    unmatched.append(i)
                    matchTable.append([t , lbl] + list(xyz))
                    continue
                searchXYZ = oXYZbyT[t]
                numLookup = atmNumByT[t]
            dist = G2mth.GetXYZDist(xyz,searchXYZ,Amat)
            while True:
                pidIndx = np.argmin(dist)
                d = dist[pidIndx]
                pid = numLookup[pidIndx]
                if atomData[pid][-1] in Ids:   #duplicate - 2 atoms on same site; invalidate & look again
                    dist[pidIndx] = 100.
                    if min(dist) == 100:
                        pid = None
                        break
                else:
                    break
            if pid is not None:
                Ids.append(atomData[pid][-1])
                matchTable.append([t , lbl] + list(xyz) + [pid, atomData[pid][0]]
                                  + atomData[pid][cx:cx+3] + [d, Ids[-1]])
            else:
                unmatched.append(i)
                matchTable.append([t , lbl] + list(xyz))
        if unmatched and unmatchedRBatoms is not None:
            unmatchedRBatoms[:] = unmatched
        return matchTable

    def OnRBAssign(event):
        '''Assign RB to atoms in a phase with tools to locate the RB in the structure
        '''

        def Draw():
            '''Create the window for assigning RB to atoms'''

            def OnAddRB(event):
                'respond to RB Add button, sets RB info in phase'
                dataGeneral = data['General']
                dataGeneral['SpnIds'] = dataGeneral.get('SpnIds',{})
                cx,ct,cs,cia = dataGeneral['AtomPtrs']
                rbType = data['testRBObj']['rbType']
                atomData = data['Atoms']
                if RigidBodies.atomsGrid:
                    matchTable = UpdateTable()
                else:
                    matchTable = assignAtoms(RBData)
                dmax = 0.
                for line in matchTable:
                    if len(line) >= 11:
                        dmax = max(dmax,line[10])
                if dmax > 1.0:
                    msg = "Atoms may not be properly located. Are you sure you want to do this?"
                    dlg = wx.MessageDialog(G2frame,msg,caption='Continue?',style=wx.YES_NO|wx.ICON_EXCLAMATION)
                    if dlg.ShowModal() != wx.ID_YES:
                        dlg.Destroy()
                        return
                    dlg.Destroy()
                Ids = []
                updateNeeded = False
                for line in matchTable:
                    if len(line) < 11:
                        if rbType == 'Spin':
                            elem = line[1]
                            x,y,z = rbObj['Orig'][0]
                        else:
                            elem = line[0]
                            x,y,z = line[2:5]
                        nextNum = len(data['Atoms'])
                        lbl = 'Rb' + elem + str(nextNum)
                        AtomAdd(x,y,z,El=elem,Name=lbl,update=False)
                        if rbType == 'Spin':
                            dataGeneral['SpnIds'][atomData[nextNum][-1]] = rbObj['RBId']
                        Ids.append(atomData[nextNum][-1])
                        updateNeeded = True
                    else:
                        atomData[line[5]][cx:cx+3] = line[2:5]
                        Ids.append(line[11])
                if len(Ids):
                    AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
                    G2lat.RBsymCheck(atomData,ct,cx,cs,AtLookUp,Amat,Ids,SGData)
                if updateNeeded:
                    SetupGeneral()
                    UpdateDrawAtoms(G2frame,data)
                    G2plt.PlotStructure(G2frame,data)

                rbNames = [[]]
                for item in data['RBModels'].get(rbType,[]):
                    if rbType == 'Spin':
                        rbNames.append(item['RBname'][0])
                    else:
                        rbNames.append(item['RBname'])
                rbObj['Ids'] = Ids          #atomids
                rbId = rbObj['RBId']        #RB obj id
                rbObj['AtomFract'] = [1.0,False]
                if rbType == 'Spin':
                    del data['testRBObj']['rbData']['Spin'][rbId]['Orig']
                    rbObj.update(data['testRBObj']['rbData']['Spin'][rbId])
                    del rbObj['rbPos']
                else:
                    rbObj['ThermalMotion'] = ['None',[0. for i in range(21)],[False for i in range(21)]] #type,values,flags
                i = 0
                while True:     #find unique name
                    rbName = '%s:%d'%(rbObj['RBname'],i)
                    if rbName in rbNames:
                        i += 1
                    else:
                        break
                rbObj['RBname'] = rbName
                if not rbType in data['RBModels']:
                    data['RBModels'][rbType] = []
                if rbType == 'Spin':    #convert items to lists of shells
                    for name in ['atColor','atType','Natoms','nSH','RBId','RBname','RBsym']:
                        item = rbObj[name]
                        rbObj[name] = [rbObj[name],]
                    rbObj['Radius'] = [[1.0,False],]
                data['RBModels'][rbType].append(copy.deepcopy(rbObj))
                RBData[rbType][rbId]['useCount'] += 1
                del data['testRBObj']

                # Update the draw atoms array & recompute bonds
                for atom in atomData:
                    ID = atom[cia+8]
                    G2mth.DrawAtomsReplaceByID(data,cia+8,atom,ID)
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data,False)
                FillRigidBodyGrid(True)

            def OnCancel(event):
                del data['testRBObj']
                FillRigidBodyGrid(True)

            def OnTorAngle(invalid,value,tc):
                '''respond to a number entered into the torsion editor.
                Update the slider (number is already saved)
                recompute atom distances
                '''
                [tor,torSlide] = Indx[tc.GetId()]
                torSlide.ChangeValue(int(value*10))
                UpdateTablePlot()

            def OnTorSlide(event):
                '''respond to the slider moving. Put value into editor & save
                recompute atom distances
                '''
                Obj = event.GetEventObject()
                tor,ang = Indx[Obj.GetId()]
                Tors = data['testRBObj']['rbObj']['Torsions'][tor]
                val = float(Obj.GetValue())/10.
                Tors[0] = val
                ang.ChangeValue(val)
                UpdateTablePlot()

            def UpdateTable(event=None):
                '''get fixed atom assignments, find closest mappings &
                update displayed table
                '''
                if not RigidBodies.atomsGrid: return []
                RigidBodies.atomsGrid.completeEdits()
                # add new atoms and reassign
                added = False
                selDict = getSelectedAtoms()
                if selDict is None: return
                for i,l in enumerate(RigidBodies.atomsTable.data):
                    if l[4] == 'Create new':
                        l[1:4] = -1,'?',-1
                        rbAssignments[i] = None
                        selDict[i] = None
                matchTable = assignAtoms(RBData,selDict)
                for i,l in enumerate(matchTable):
                    if len(l) < 11: continue
                    RigidBodies.atomsTable.data[i][1:4] = l[5],l[6],l[10]
                if RigidBodies.atomsGrid:
                    RigidBodies.atomsGrid.ForceRefresh()
                if added: wx.CallLater(100,Draw)
                return matchTable

            def OnAzSlide(event):
                '''respond to the azimuth slider moving.
                Save & put value into azimuth edit widget; show Q &
                update the plot and table
                '''
                Obj = event.GetEventObject()
                rbObj['OrientVec'][0] = float(Obj.GetValue())/10.
                for i in range(4):
                    val = rbObj['OrientVec'][i]
                    G2frame.testRBObjSizers['OrientVecSiz'][i].ChangeValue(val)
                Q = G2mth.AVdeg2Q(rbObj['OrientVec'][0],
                                np.inner(Amat,rbObj['OrientVec'][1:]))
                rbObj['Orient'][0] = Q
                A,V = G2mth.Q2AVdeg(Q)
                rbObj['OrientVec'][1:] = np.inner(Bmat,V)
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                UpdateTable()

            def UpdateOrientation(*args,**kwargs):
                '''Respond to a change in the azimuth or vector via
                the edit widget; update Q display & azimuth slider
                '''
                Q = G2mth.AVdeg2Q(rbObj['OrientVec'][0],
                    np.inner(Amat,rbObj['OrientVec'][1:]))
                rbObj['Orient'][0] = Q
                try:
                    G2frame.testRBObjSizers['OrientVecSiz'][4].ChangeValue(
                        int(10*rbObj['OrientVec'][0]))
                except:
                    pass
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                UpdateTable()

            def UpdateTablePlot(*args,**kwargs):
                '''update displayed table and plot
                '''
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                UpdateTable()

            def UpdateSytSym(*args,**kwargs):

                Sytsym,Mult = G2spc.SytSym(rbObj['Orig'][0],data['General']['SGData'])[:2]
                sytsymtxt.SetLabel('Origin site symmetry: %s, multiplicity: %d '%(Sytsym,Mult))
                UpdateTablePlot(args,kwargs)

            def getSelectedAtoms():
                'Find the FB atoms that have been assigned to specific atoms in structure'
                if not RigidBodies.atomsGrid: return
                RigidBodies.atomsGrid.completeEdits()
                tbl = RigidBodies.atomsGrid.GetTable()
                selDict = {}
                dups = []
                assigned = []
                atomsLabels = [a[0] for a in data['Atoms']]
                for r in range(tbl.GetRowsCount()):
                    sel = tbl.GetValue(r,4).strip()
                    if sel == 'Create new': continue # ignore positions of new atoms
                    if sel not in atomsLabels: continue
                    atmNum = atomsLabels.index(sel)
                    if atmNum < 0: continue
                    if atmNum in assigned:
                        if sel not in dups:
                            dups.append(sel)
                    else:
                        assigned.append(atmNum)
                    selDict[r] = atmNum
                if dups:
                    msg = 'Error: The following atom(s) are assigned multiple times: '
                    for i in dups:
                        msg += i
                        msg += ', '
                    wx.MessageBox(msg[:-2],caption='Duplicated Fixed Atoms',style=wx.ICON_EXCLAMATION)
                    return
                return selDict

            def getDeltaXYZ(selDict,data,rbObj):
                '''Evaluate the RB position & return the difference between the coordinates
                and RB coordinates with origin and orientation from data['testRBObj']
                '''
                Amat,Bmat = G2lat.cell2AB(data['General']['Cell'][1:7])
                rbXYZ = G2mth.UpdateRBXYZ(Bmat,rbObj,data['testRBObj']['rbData'],data['testRBObj']['rbType'])[0]
                phaseXYZ = G2mth.getAtomXYZ(atomData,cx)
                deltaList = []
                for i in selDict:
                    if selDict[i] is None: continue
                    deltaList.append(phaseXYZ[selDict[i]]-rbXYZ[i])
                return np.array(deltaList)

            def objectiveDeltaPos(vals,selDict,data,rbObj_in):
                '''Objective function for minimization.
                Returns a list of distances between atom positions and
                located rigid body positions

                :param list vals: a 4 or 7 element array with 4 quaterian values
                   or 4 quaterian values followed by 3 origin (x,y,z) values

                :returns: the 3*n distances between the n selected atoms
                '''
                rbObj = copy.deepcopy(rbObj_in)
                rbObj['Orient'][0][:] = G2mth.normQ(vals[:4])
                if len(vals) == 7:
                    rbObj['Orig'][0][:] = vals[4:]
                #print(np.sqrt(sum(getDeltaXYZ(selDict,data,rbObj).flatten()**2)))
                return getDeltaXYZ(selDict,data,rbObj).flatten()

            def onSetOrigin(event):
                'Set Origin to best fit selected atoms'
                if rbObj['fixOrig']:
                    wx.MessageBox('Not possible with origin locked',caption='Οrigin Locked',
                                      style=wx.ICON_EXCLAMATION)
                    return
                selDict = getSelectedAtoms()
                if len(selDict) < 1:
                    wx.MessageBox('No existing atoms were selected',caption='Select Atom(s)',
                                      style=wx.ICON_EXCLAMATION)
                    return
                deltaList = getDeltaXYZ(selDict,data,rbObj)
                data['testRBObj']['rbObj']['Orig'][0] += deltaList.sum(axis=0)/len(deltaList)
                for i,item in enumerate(Xsizers):
#                    item.SetValue(data['testRBObj']['rbObj']['Orig'][0][i])
                    item.ChangeValue(data['testRBObj']['rbObj']['Orig'][0][i])
                UpdateTablePlot()

            def onFitOrientation(event):
                'Set Orientation to best fit selected atoms'
                selDict = getSelectedAtoms()
                if len(selDict) < 2:
                    wx.MessageBox('At least two existing atoms must be selected',caption='Select Atoms',
                                      style=wx.ICON_EXCLAMATION)
                    return
                vals = rbObj['Orient'][0][:] #+ rbObj['Orig'][0][:]
                out = so.leastsq(objectiveDeltaPos,vals,(selDict,data,rbObj))
                data['testRBObj']['rbObj']['Orient'][0][:] = G2mth.normQ(out[0])
                updateAddRBorientText(G2frame,data['testRBObj'],Bmat)
                UpdateTablePlot()

            def onFitBoth(event):
                'Set Orientation and origin to best fit selected atoms'
                if rbObj['fixOrig']:
                    wx.MessageBox('Not possible with origin locked',caption='Οrigin Locked',
                                      style=wx.ICON_EXCLAMATION)
                    return
                selDict = getSelectedAtoms()
                if len(selDict) < 3:
                    wx.MessageBox('At least three existing atoms must be selected',caption='Select Atoms',
                                      style=wx.ICON_EXCLAMATION)
                    return
                vals = np.concatenate((rbObj['Orient'][0], rbObj['Orig'][0]))
                out = so.leastsq(objectiveDeltaPos,vals,(selDict,data,rbObj))
                data['testRBObj']['rbObj']['Orig'][0][:] = out[0][4:]
                data['testRBObj']['rbObj']['Orient'][0][:] = G2mth.normQ(out[0][:4])
                for i,item in enumerate(Xsizers):
#                    item.SetValue(data['testRBObj']['rbObj']['Orig'][0][i])
                    item.ChangeValue(data['testRBObj']['rbObj']['Orig'][0][i])
                updateAddRBorientText(G2frame,data['testRBObj'],Bmat)
                UpdateTablePlot()

            def BallnSticks(event):
                '''Set all draw atoms in crystal structure to balls & stick
                '''
                for a in data['Drawing']['Atoms']: a[6] = "balls & sticks"
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                RigidBodies.SetFocus() # make sure tab presses go to panel

            def Sticks(event):
                '''Set all draw atoms in crystal structure to stick
                '''
                for a in data['Drawing']['Atoms']: a[6] = "sticks"
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                RigidBodies.SetFocus() # make sure tab presses go to panel

            def OnRowSelect(event):
                '''Respond to the selection of a rigid body atom.
                Highlight the atom in the body and the paired atom in the
                crystal
                '''
                event.Skip()
                cryatom = event.GetEventObject().Table.GetValue(event.GetRow(),4)
                if not cryatom:
                    cryatom = event.GetEventObject().Table.GetValue(event.GetRow(),2)
                data['testRBObj']['RBhighLight'] = event.GetRow()
                data['testRBObj']['CRYhighLight'] = [
                    i for i,a in enumerate(data['Atoms']) if a[0] == cryatom]
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)
                misc['showSelect'].setByString(cryatom)
                G2frame.Raise()
                #RigidBodies.atomsGrid.SetFocus()

            def OnSymRadioSet(event):
                '''Set the symmetry axis for the body as
                data['testRBObj']['rbObj']['symAxis']. This may never be
                set, so use data['testRBObj']['rbObj'].get('symAxis') to
                access this so the default value is None.
                '''
                axis = (None,[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,1,1]
                    )[event.GetEventObject().GetSelection()]
                if axis:
                    axis = np.array(axis)/nl.norm(axis)
                data['testRBObj']['rbObj']['symAxis'] = axis
                UpdateTablePlot()

            def OnOrigSet(event):
                data['testRBObj']['rbObj']['Orig'][0] = data['Drawing']['viewPoint'][0]
                for i,item in enumerate(Xsizers):
                    item.ChangeValue(data['testRBObj']['rbObj']['Orig'][0][i])
                UpdateSytSym()
                UpdateTablePlot()

            showAtom = [None]
            def showCryAtom(*args,**kwargs):
                '''Respond to selection of a crystal atom
                '''
                data['testRBObj']['CRYhighLight'] = [
                    i for i,a in enumerate(data['Atoms']) if a[0] == showAtom[0]]
                G2plt.PlotStructure(G2frame,data,False,UpdateTable)

            # Start of Draw()
            RBdirlbl = [' x ',' y ',' z ','x+y','x+y+z']
            if not data['testRBObj']: return
            if RigidBodies.GetSizer(): RigidBodies.GetSizer().Clear(True)
            unmatchedRBatoms = []
            matchTable = assignAtoms(RBData,unmatchedRBatoms=unmatchedRBatoms)
            if unmatchedRBatoms:
                # msg = 'There are {} atoms that will need to be added to atoms list.'.format(len(unmatchedRBatoms))
                # G2G.G2MessageBox(G2frame,msg,title='Please note')
                for i in unmatchedRBatoms:
                    rbAssignments[i] = None
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add((5,5),0)
            rbObj = data['testRBObj']['rbObj']
            data['RBModels']['SpnIds'] = data['RBModels'].get('SpnIds',{})
            rbType = data['testRBObj']['rbType']
            rbName = rbObj['RBname']
            rbId = rbObj['RBId']
            Torsions = rbObj['Torsions']
            atNames = data['testRBObj']['atNames']
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            helpText = '''
This window is used to insert a rigid body into a unit cell, determining
the initial settings for the position and orientation of the body, as well as
any internal torsion angles. The origin determines where the origin of the rigid body
will be placed in the unit cell (expressed in fractional coordinates). The
orientation is determined by a quaternion that is expressed here as vector
(in fraction coordinates) and an azimuthal rotation (in degrees) around that
quaternion vector. For systems where the rigid body is placed on a
crystallographic symmetry element, the "Rigid body symmetry axis"
("x", "y", "z", "x+y" or "x+y+z") specifies the rotation axis that aligns with
an allowed rotation for this symmetry element (NB: could be perpendicular to a
mirror). This places the selected Cartesian axis along the quaternion vector.
The quaternion vector may be set by the user to force the rotation to be about
a particular crystallographic direction (e.g. along a rotation axis or
perpendicular to a mirror). The rotation action can be tested via the slider.

%%If there are atoms in the unit cell that are of the appropriate type and
are not already assigned to rigid bodies, a table shows each atom in the rigid
body and the closest crystallographic atom, and the distance between them.
Set this pairing by using the pulldown menu in the "Assign as atom" column by
selecting a particular atom. If the selection is changed, "Update Assignments"
recomputes distance between paired atoms. If one atom is paired manually using
"Assign as", the "Set Origin" button can be used to place the rigid body origin
to best fit the paired atom(s). Likewise, with two or more atoms assigned, the
"Set Orientation" button will determine the quaternion azimuth and vector. Three
or more pairs are assignments allow use of the "Set both" button, which
sets the orientation and origin to best give the smallest distances between
the assigned atoms. NB: if apparently stuck with a poor fit, try shifting the
Orientation azimuth slider and try Set both again.
Note that when a row in the table is selected, the corresponding atoms
are highlighted in green. The tab key or the "Crystal Highlight" pulldown
can be used to highlight differing unit cell atoms. Alt-Tab highlights different
RB atoms.

%%If there are no unassigned atoms of the right type (existing RBs will have
orange sticks for bonds), then the table will show "Create new" in the
"Assign as atom" column. The proposed RB can be positioned via Alt mouse
operations and/or by entering appropriate values in the Orientation azimuth
and vector x,y,z boxes. Symmetry element issues should be attended to by
proper selections as noted above. "Add" will then add the rigid body atoms to
the Atom list.

%%The GSAS-II graphics window shows the unit cell contents (use the
"Ball & Sticks" or "Sticks" buttons to change the display of this) and
the rigid body is shown with small balls and green sticks. At the origin of the
RB the axes are indicated with red, green and blue lines (for x, y, & z).
A white line indicates the quaternion vector direction.
The mouse can also be used to position the rigid body in the plot by holding
the Alt key down while dragging with the mouse: Alt+left button to rotates
the RB in the screen x & y; Alt+middle mouse to rotate the RB around the viewing
direction; Alt+right mouse translates the RB in the screen x-y plane.
Note that dragging the mouse without the Alt button changes the view
of the crystal structure.
%%Once the rigid body has been placed in the desired position, press the "Add" button.
'''
            topSizer.Add(wx.StaticText(RigidBodies,label='Locating rigid body: '+rbName), 0, WACV)
            topSizer.Add((10,-1),0)
            topSizer.Add(wx.StaticText(RigidBodies,label='Display crystal structure as:'), 0, WACV)
            btn = wx.Button(RigidBodies,label="Ball && Sticks")
            btn.Bind(wx.EVT_BUTTON, BallnSticks)
            topSizer.Add(btn)
            btn = wx.Button(RigidBodies,label="Sticks")
            btn.Bind(wx.EVT_BUTTON, Sticks)
            topSizer.Add(btn)
            topSizer.Add((-1,-1),1,wx.EXPAND,1)
            topSizer.Add(G2G.HelpButton(RigidBodies,helpText,wrap=700),0,wx.RIGHT,5)
            mainSizer.Add(topSizer,0,wx.EXPAND,0)
            mainSizer.Add((-1,5),0)
            OriSizer = wx.BoxSizer(wx.HORIZONTAL)
            OriSizer.Add(wx.StaticText(RigidBodies,label='Origin: '),0,WACV)
            Xsizers = []
            for ix in range(3):
                OriSizer.Add(wx.StaticText(RigidBodies,label=RBdirlbl[ix]),0,WACV,4)
                origX = G2G.ValidatedTxtCtrl(RigidBodies,rbObj['Orig'][0],ix,nDig=(10,5),
                    xmin=-1.5,xmax=1.5,typeHint=float,OnLeave=UpdateSytSym)
                OriSizer.Add(origX,0,WACV)
                Xsizers.append(origX)
            try:
                rbObj['fixOrig']
            except:
                rbObj['fixOrig'] = False
            if rbType != 'Spin':
                fixOrig = G2G.G2CheckBox(RigidBodies,'Lock',rbObj,'fixOrig')
                OriSizer.Add(fixOrig,0,WACV,10)
            if not rbObj['fixOrig']:
                origSet = wx.Button(RigidBodies,label='Set to view point')
                origSet.Bind(wx.EVT_BUTTON,OnOrigSet)
                OriSizer.Add(origSet,0,WACV)
                mainSizer.Add(OriSizer)
            Sytsym,Mult = G2spc.SytSym(rbObj['Orig'][0],data['General']['SGData'])[:2]
            sytsymtxt = wx.StaticText(RigidBodies,label='Origin site symmetry: %s, multiplicity: %d '%(Sytsym,Mult))
            mainSizer.Add(sytsymtxt)
            if rbType != 'Spin':
                OriSizer1 = wx.FlexGridSizer(0,5,5,5)
                if len(atomData):
                    choice = list(atNames[0].keys())
                    choice.sort()
                    G2frame.testRBObjSizers.update({'Xsizers':Xsizers})
                mainSizer.Add(OriSizer1)
                mainSizer.Add((5,5),0)
                OriSizer2 = wx.BoxSizer(wx.HORIZONTAL)
                if 'OrientVec' not in rbObj: rbObj['OrientVec'] = [0.,0.,0.,0.]
                rbObj['OrientVec'][0],V = G2mth.Q2AVdeg(rbObj['Orient'][0])
                rbObj['OrientVec'][1:] = np.inner(Bmat,V)
                OriSizer2.Add(wx.StaticText(RigidBodies,label='Orientation azimuth: '),0,WACV)
                OrientVecSiz = []
                OrientVecSiz.append(G2G.ValidatedTxtCtrl(RigidBodies,rbObj['OrientVec'],0,nDig=(10,2),
                    xmin=0.,xmax=360.,typeHint=float,OnLeave=UpdateOrientation))
                OriSizer2.Add(OrientVecSiz[-1],0,WACV)
                azSlide = G2G.G2Slider(RigidBodies,style=wx.SL_HORIZONTAL,size=(200,25),
                    minValue=0,maxValue=3600,value=int(10*rbObj['OrientVec'][0]))
                azSlide.Bind(wx.EVT_SLIDER, OnAzSlide)
                OriSizer2.Add(azSlide,0,WACV)
                mainSizer.Add(OriSizer2)
                OriSizer3 = wx.BoxSizer(wx.HORIZONTAL)
                OriSizer3.Add(wx.StaticText(RigidBodies,label='Orientation vector'),0,WACV)
                for ix,lbl in enumerate('xyz'):
                    OriSizer3.Add(wx.StaticText(RigidBodies,label='  '+lbl+': '),0,WACV)
                    OrientVecSiz.append(G2G.ValidatedTxtCtrl(RigidBodies,rbObj['OrientVec'],ix+1,nDig=(10,4),
                        xmin=-1.,xmax=1.,typeHint=float,OnLeave=UpdateOrientation))
                    OriSizer3.Add(OrientVecSiz[-1],0,WACV)
                OrientVecSiz.append(azSlide)
                OriSizer3.Add(wx.StaticText(RigidBodies,label=' (frac coords)'),0,WACV)
                G2frame.testRBObjSizers.update({'OrientVecSiz':OrientVecSiz})
                mainSizer.Add(OriSizer3)
                mainSizer.Add((5,5),0)
                OriSizer4 = wx.BoxSizer(wx.HORIZONTAL)
                OriSizer4.Add(wx.StaticText(RigidBodies,label='Rigid body symmetry axis: '),0, WACV)
                choices = ['None']+RBdirlbl
                symRadioSet = wx.RadioBox(RigidBodies,choices=choices)
                symRadioSet.Bind(wx.EVT_RADIOBOX, OnSymRadioSet)
                OriSizer4.Add(symRadioSet)
                mainSizer.Add(OriSizer4)
                mainSizer.Add((5,5),0)
                RefSizer = wx.FlexGridSizer(0,7,5,5)
                mainSizer.Add(RefSizer)
                mainSizer.Add((5,5),0)
                if Torsions:
                    rbSeq = RBData['Residue'][rbId]['rbSeq']
                    TorSizer = wx.FlexGridSizer(0,4,5,5)
                    TorSizer.AddGrowableCol(1,1)
                    for t,[torsion,seq] in enumerate(zip(Torsions,rbSeq)):
                        torName = ''
                        for item in [seq[0],seq[1],seq[3][0]]:
                            if data['testRBObj']['NameLookup'][item]:
                                torName += data['testRBObj']['NameLookup'][item]+' '
                            else:
                                torName += data['testRBObj']['rbAtTypes'][item]+str(item)+' '
                        TorSizer.Add(wx.StaticText(RigidBodies,label='Side chain torsion for rb seq: '+torName),0,WACV)
                        torSlide = G2G.G2Slider(RigidBodies,style=wx.SL_HORIZONTAL,minValue=0,maxValue=3600,value=int(torsion[0]*10.))
                        torSlide.Bind(wx.EVT_SLIDER, OnTorSlide)
                        TorSizer.Add(torSlide,1,wx.EXPAND|wx.RIGHT)
                        TorSizer.Add(wx.StaticText(RigidBodies,label=' Angle: '),0,WACV)
                        ang = G2G.ValidatedTxtCtrl(RigidBodies,torsion,0,nDig=(8,3),typeHint=float,OnLeave=OnTorAngle)
                        Indx[torSlide.GetId()] = [t,ang]
                        Indx[ang.GetId()] = [t,torSlide]
                        TorSizer.Add(ang,0,WACV)
                    mainSizer.Add(TorSizer,0,wx.EXPAND|wx.RIGHT)
                else:
                    mainSizer.Add(wx.StaticText(RigidBodies,label='No side chain torsions'),0)
            OkBtn = wx.Button(RigidBodies,label="Add")
            OkBtn.Bind(wx.EVT_BUTTON, OnAddRB)
            CancelBtn = wx.Button(RigidBodies,label='Cancel')
            CancelBtn.Bind(wx.EVT_BUTTON, OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            if OkBtn: btnSizer.Add(OkBtn)
            btnSizer.Add(CancelBtn)
            btnSizer.Add((20,20),1)
            mainSizer.Add(btnSizer,0,wx.BOTTOM|wx.TOP, 20)

            SetPhaseWindow(RigidBodies,mainSizer)
            data['testRBObj']['RBhighLight'] = None
            data['testRBObj']['CRYhighLight'] = []
            assignable = [a[0] for a in data['Atoms'] if a[-1] not in rbUsedIds]
            data['testRBObj']['availAtoms'] = ['         '] + assignable
            if len(assignable) == 0 and rbType != 'Spin':
                mainSizer.Add(wx.StaticText(RigidBodies,label=
                    'No matching atoms between rigid body and crystal.'+
                    ' All rigid body atoms will be added to structure.'),0)
                misc['UpdateTable'] = None
                mainSizer.Layout()
#                G2plt.PlotStructure(G2frame,data,True)
                RigidBodies.atomsGrid = None
                return

            G2plt.PlotStructure(G2frame,data,True,UpdateTable)

            if rbType == 'Spin':
                mainSizer.Add(wx.StaticText(RigidBodies,label=' Spinning rigid body:'),0)
            else:
                mainSizer.Add(wx.StaticText(RigidBodies,label=
                    'Match between atoms in rigid body and crystal. Use assignments to align bodies.'),0)
            mainSizer.Add((5,5))
            gridSizer = wx.BoxSizer(wx.HORIZONTAL)
            colLabels = ['RB\ntype','Atom\n#','Atom\nlabel','delta, A','Assign as atom']
            rowLabels = [l[1] for l in matchTable]
            displayTable = []
            Actions = True
            for i,l in enumerate(matchTable):
                lbl = ''
                if i in rbAssignments:
                    if rbAssignments[i] is None:
                        displayTable.append([l[0],-1,'?',-1,'Create new'])
                    else:
                        lbl = rbAssignments[i]
                        displayTable.append([l[0],l[5],l[6],l[10],lbl])
                else:
                    if rbType == 'Spin':
                        displayTable.append([l[1],-1,l[1]+str(i),-1,'Create new'])
                    else:
                        displayTable.append([l[0],l[5],l[6],l[10],lbl])
            Types = [wg.GRID_VALUE_STRING, wg.GRID_VALUE_NUMBER,
                     wg.GRID_VALUE_STRING, wg.GRID_VALUE_FLOAT+':8,3',wg.GRID_VALUE_STRING]
            RigidBodies.atomsTable = G2G.Table(displayTable,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            RigidBodies.atomsGrid = G2G.GSGrid(RigidBodies)
            RigidBodies.atomsGrid.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            choiceeditor = wg.GridCellChoiceEditor(data['testRBObj']['availAtoms']+['Create new'], False)
            RigidBodies.atomsGrid.SetTable(RigidBodies.atomsTable,True)
            # make all grid entries read only
            attr = wg.GridCellAttr()
            attr.SetReadOnly(True)
            for i in range(len(colLabels)-1):
                attr.IncRef()
                RigidBodies.atomsGrid.SetColAttr(i, attr)
            attr = wg.GridCellAttr()
            attr.SetAlignment(wx.ALIGN_RIGHT,wx.ALIGN_CENTRE)
            RigidBodies.atomsGrid.SetColAttr(1, attr)
            attr = wg.GridCellAttr()
            attr.SetEditor(choiceeditor)
            RigidBodies.atomsGrid.SetColAttr(4, attr)

            RigidBodies.atomsGrid.AutoSizeColumns(True)
            RigidBodies.atomsGrid.SetMargins(0,0)

            gridSizer.Add(RigidBodies.atomsGrid)#,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            gridSizer.Add((5,5))

            if Actions:
                btnSizer = wx.BoxSizer(wx.VERTICAL)
                hSizer = wx.BoxSizer(wx.HORIZONTAL)
                hSizer.Add(wx.StaticText(RigidBodies,label='Crystal Highlight: '))
                misc['showSelect'] = G2G.G2ChoiceButton(RigidBodies,
                    data['testRBObj']['availAtoms'],None,None,showAtom,0,showCryAtom,size=(75,-1))
                hSizer.Add(misc['showSelect'])
                btnSizer.Add(hSizer)
                btnSizer.Add((-1,20))
                btnSizer.Add(wx.StaticText(RigidBodies,label='Actions with assigned\natom(s)...'),0,wx.ALL)
                btnSizer.Add((-1,10))
                btn = wx.Button(RigidBodies,label='Update Assignments')
                btn.Bind(wx.EVT_BUTTON,UpdateTable)
                btnSizer.Add(btn,0,wx.ALIGN_CENTER)
                btnSizer.Add((-1,10))

                btn = wx.Button(RigidBodies,label='Set Origin')
                btn.Bind(wx.EVT_BUTTON,onSetOrigin)
                btnSizer.Add(btn,0,wx.ALIGN_CENTER)
                btnSizer.Add((-1,5))
                btn = wx.Button(RigidBodies,label='Set Orientation')
                btn.Bind(wx.EVT_BUTTON,onFitOrientation)
                btnSizer.Add(btn,0,wx.ALIGN_CENTER)
                btnSizer.Add((-1,5))
                btn = wx.Button(RigidBodies,label='Set both')
                btn.Bind(wx.EVT_BUTTON,onFitBoth)
                btnSizer.Add(btn,0,wx.ALIGN_CENTER)
                gridSizer.Add(btnSizer)
            mainSizer.Add(gridSizer)

            mainSizer.Layout()
            RigidBodies.SetScrollRate(10,10)
            RigidBodies.SendSizeEvent()
            RigidBodies.Scroll(0,0)
            RigidBodies.SetFocus() # make sure tab presses go to panel
            misc['UpdateTable'] = UpdateTable

        # start of OnRBAssign(event)
        rbAssignments = {}
        rbUsedIds = []   # Ids of atoms in current phase used inside RBs
        for i in data['RBModels']:
            if i == 'SpnIds':
                continue
            for j in data['RBModels'][i]:
                rbUsedIds += j['Ids']
        G2frame.GetStatusBar().SetStatusText('',1)
        RBData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbNames = {}
        dups = []
        for rbVec in RBData['Vector']:
            if rbVec == 'AtInfo': continue
            key = RBData['Vector'][rbVec]['RBname']
            if key in rbNames:
                dups.append(key)
            else:
                rbNames[key] = ['Vector',rbVec]
        for rbRes in RBData['Residue']:
            if rbRes == 'AtInfo': continue
            key = RBData['Residue'][rbRes]['RBname']
            if key in rbNames:
                dups.append(key)
            else:
                rbNames[key] = ['Residue',rbRes]
        if data['General']['Type'] == 'nuclear': #exclude other phase types
            for rbSpn in RBData.get('Spin',{}):     #patch get for old files
                key = RBData['Spin'][rbSpn]['RBname']
                if key in rbNames:
                    dups.append(key)
                else:
                    rbNames[key] = ['Spin',rbSpn]

        if dups:
            msg = 'Two or more rigid bodies have the same name. This must be corrected before bodies can be added.'
            msg += '\n\n Duplicated name(s): '
            for d in dups:
                msg += d
                msg += ' '
            G2G.G2MessageBox(G2frame,msg,'Duplicate Rigid Body Names')
            return
        if not rbNames:
            print ('**** ERROR - no rigid bodies defined ****')
            G2G.G2MessageBox(G2frame,
                'You must define bodies in the Rigid Bodies tree item before trying to add them to a phase',
                'No Rigid Bodies')
            return
        general = data['General']
        SGData = general['SGData']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        cx,ct,cs,cia = general['AtomPtrs']
        atomData = data['Atoms']
        #AtLookUp = G2mth.FillAtomLookUp(atomData,cia+8)
        Indx = {}
        data['testRBObj'] = {}
        if len(rbNames) == 1:
            selection = list(rbNames.keys())[0]
        else:
            choices = sorted(list(rbNames.keys()))
            dlg = G2G.G2SingleChoiceDialog(
                G2frame,'Select rigid body to\nadd to structure',
                'Select rigid body',choices)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    selection = choices[dlg.GetSelection()]
                else:
                    return
            finally:
                dlg.Destroy()
            if selection not in rbNames:
                print('Invalid RB selection',selection,'How did this happen?')
                return
        rbType,rbId = rbNames[selection]
        if rbType == 'Spin':
            data['testRBObj']['rbAtTypes'] = [RBData[rbType][rbId]['rbType'],]
            data['testRBObj']['AtInfo'] = {RBData[rbType][rbId]['rbType']:[1.0,(128, 128, 255)],}
            data['testRBObj']['rbType'] = rbType
            data['testRBObj']['rbData'] = RBData
            data['testRBObj']['Sizers'] = {}
            rbRef = RBData[rbType][rbId].get('rbRef',[])    #get only for spinning RBs
            data['testRBObj']['rbRef'] = rbRef
        else:
            data['testRBObj']['rbAtTypes'] = RBData[rbType][rbId]['rbTypes']
            data['testRBObj']['AtInfo'] = RBData[rbType]['AtInfo']
            data['testRBObj']['NameLookup'] = RBData[rbType][rbId].get('atNames',[])    #only for residues
            data['testRBObj']['rbType'] = rbType
            data['testRBObj']['rbData'] = RBData
            data['testRBObj']['Sizers'] = {}
            rbRef = RBData[rbType][rbId]['rbRef']
            data['testRBObj']['rbRef'] = rbRef
        refType = []
        #refName = []
        for ref in rbRef[:3]:
            reftype = data['testRBObj']['rbAtTypes'][ref]
            refType.append(reftype)
            #refName.append(reftype+' '+str(rbRef[0]))
        atNames = [{},{},{}]
        AtNames = {}
        for iatm,atom in enumerate(atomData):
            AtNames[atom[ct-1]] = iatm
            for i,reftype in enumerate(refType):
                if atom[ct] == reftype:
                    atNames[i][atom[ct-1]] = iatm
        data['testRBObj']['atNames'] = atNames
        data['testRBObj']['AtNames'] = AtNames
        data['testRBObj']['rbObj'] = {'Orig':[[0,0,0],False],
                    'Orient':[[0.,0.,0.,1.],' '],'Ids':[],'RBId':rbId,'Torsions':[],
                    'numChain':'','RBname':RBData[rbType][rbId]['RBname']}
        # if rbType == 'Spin':
        #     data['testRBObj']['rbObj']['Radius'] = [1.0,False]
        data['testRBObj']['torAtms'] = []
        for item in RBData[rbType][rbId].get('rbSeq',[]):
            data['testRBObj']['rbObj']['Torsions'].append([item[2],False])
            data['testRBObj']['torAtms'].append([-1,-1,-1])
        wx.CallAfter(Draw)

    def OnAutoFindResRB(event):
        RBData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbKeys = list(RBData['Residue'].keys())
        rbKeys.remove('AtInfo')
        if not len(rbKeys):
            print ('**** ERROR - no residue rigid bodies are defined ****')
            return
        RBNames = [RBData['Residue'][k]['RBname'] for k in rbKeys]
        RBIds = dict(zip(RBNames,rbKeys))
        general = data['General']
        cx,ct,cs,cia = general['AtomPtrs']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        Atoms = data['Atoms']
        AtLookUp = G2mth.FillAtomLookUp(Atoms,cia+8)
        if 'macro' not in general['Type']:
            print ('**** ERROR - this phase is not a macromolecule ****')
            return
        if not len(Atoms):
            print ('**** ERROR - this phase has no atoms ****')
            return
        RBObjs = []
        cx,ct = general['AtomPtrs'][:2]
        iatm = 0
        wx.BeginBusyCursor()
        try:
            isave = 0
            while iatm < len(Atoms):
                atom = Atoms[iatm]
                res = atom[1].strip()
                numChain = ' %s %s'%(atom[0],atom[2])
                if res not in RBIds or atom[ct-1] == 'OXT':
                    iatm += 1
                    continue        #skip for OXT, water molecules, etc.
                rbRes = RBData['Residue'][RBIds[res]]
                rbRef = rbRes['rbRef']
                VAR = rbRes['rbXYZ'][rbRef[1]]-rbRes['rbXYZ'][rbRef[0]]
                VBR = rbRes['rbXYZ'][rbRef[2]]-rbRes['rbXYZ'][rbRef[0]]
                incr = 1
                isave = 0
                if 'N  A' in atom[3]:   #disordered residue - read every other atom
                    isave = iatm+1
                    incr = 2
                if 'N  B' in atom[3]:
                    incr = 2
                rbObj = {'RBname':rbRes['RBname']+':'+str(rbRes['useCount']),'numChain':numChain}
                rbAtoms = []
                rbIds = []
                for iratm in range(len(rbRes['atNames'])):
                    if res != Atoms[iatm][1]:
                        print('Short %s residue at atom %d is skipped'%(Atoms[iatm][1],iatm))
                        break
                    rbAtoms.append(np.array(Atoms[iatm][cx:cx+3]))
                    rbIds.append(Atoms[iatm][20])
                    iatm += incr    #puts this at beginning of next residue?
                else:
                    if 'N  B' in atom[3]:   #end of disorder - reset next atom position
                        iatm -= 1
                        incr = 1
                    Orig = rbAtoms[rbRef[0]]
                    rbObj['RBId'] = RBIds[res]
                    rbObj['Ids'] = rbIds
                    rbObj['Orig'] = [Orig,False]
    #                print ' residue '+rbRes['RBname']+str(atom[0]).strip()+ \
    #                    ' origin at: ','%.5f %.5f %.5f'%(Orig[0],Orig[1],Orig[2])
                    VAC = np.inner(Amat,rbAtoms[rbRef[1]]-Orig)
                    VBC = np.inner(Amat,rbAtoms[rbRef[2]]-Orig)
                    VCC = np.cross(VAR,VAC)
                    QuatA = G2mth.makeQuat(VAR,VAC,VCC)[0]
                    VAR = G2mth.prodQVQ(QuatA,VAR)
                    VBR = G2mth.prodQVQ(QuatA,VBR)
                    QuatB = G2mth.makeQuat(VBR,VBC,VAR)[0]
                    QuatC = G2mth.prodQQ(QuatB,QuatA)
                    rbObj['Orient'] = [QuatC,' ']
                    rbObj['AtomFract'] = [1.0,False]
                    rbObj['ThermalMotion'] = ['None',[0. for i in range(21)],[False for i in range(21)]] #type,values,flags
                    SXYZ = []
                    TXYZ = []
                    rbObj['Torsions'] = []
                    for i,xyz in enumerate(rbRes['rbXYZ']):
                        SXYZ.append(G2mth.prodQVQ(QuatC,xyz))
                        TXYZ.append(np.inner(Amat,rbAtoms[i]-Orig))
                    for Oatm,Patm,x,Riders in rbRes['rbSeq']:
                        VBR = SXYZ[Oatm]-SXYZ[Patm]
                        VAR = SXYZ[Riders[0]]-SXYZ[Patm]
                        VAC = TXYZ[Riders[0]]-TXYZ[Patm]
                        QuatA,D = G2mth.makeQuat(VAR,VAC,VBR)
                        ang = 180.*D/np.pi
                        rbObj['Torsions'].append([ang,False])
                        for ride in Riders:
                            SXYZ[ride] = G2mth.prodQVQ(QuatA,SXYZ[ride]-SXYZ[Patm])+SXYZ[Patm]
                    rbRes['useCount'] += 1
                    RBObjs.append(rbObj)
                if isave:
                    iatm = isave
            data['RBModels']['Residue'] = RBObjs
            for RBObj in RBObjs:
                newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')[0]
                for i,Id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[Id]][cx:cx+3] = newXYZ[i]
        finally:
            wx.EndBusyCursor()
        wx.CallAfter(FillRigidBodyGrid,True)

    def OnRBRemoveAll(event):
        data['RBModels']['Residue'] = []
        data['RBModels']['Vector'] = []
        data['RBModels']['Spin'] = []
        RBData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        for RBType in ['Vector','Residue','Spin']:
            for rbId in RBData[RBType]:
                RBData[RBType][rbId]['useCount'] = 0
        FillRigidBodyGrid(True)

    def OnGlobalResRBTherm(event):
        RBObjs = data['RBModels']['Residue']
        names = ['None','Uiso','T','TL','TLS']
        cia = data['General']['AtomPtrs'][3]
        AtLookUp = G2mth.FillAtomLookUp(data['Atoms'],cia+8)
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Residue thermal motion model',names)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            parm = names[sel]
            Ttype = 'A'
            if parm == 'Uiso':
                Ttype = 'I'
            for rbObj in RBObjs:
                rbObj['ThermalMotion'][0] = parm
                if parm != 'None':
                    for i,Id in enumerate(rbObj['Ids']):
                        data['Atoms'][AtLookUp[Id]][cia] = Ttype
        dlg.Destroy()
        wx.CallAfter(FillRigidBodyGrid,True)

    def OnGlobalResRBRef(event):
        RBObjs = data['RBModels']['Residue']
        names = ['Origin','Orient. angle','Full Orient.']
        nTor = 0
        for rbObj in RBObjs:
            nTor = max(nTor,len(rbObj['Torsions']))
        names += ['Torsion '+str(i) for i in range(nTor)]
        if np.any([rbObj['ThermalMotion'][0] == 'Uiso' for rbObj in RBObjs]):
           names += ['Uiso',]
        if np.any([rbObj['ThermalMotion'][0] == 'TLS' for rbObj in RBObjs]):
           names += ['Tii','Tij','Lii','Lij','Sij']
        elif np.any([rbObj['ThermalMotion'][0] == 'TL' for rbObj in RBObjs]):
           names += ['Tii','Tij','Lii','Lij']
        elif np.any([rbObj['ThermalMotion'][0] == 'T' for rbObj in RBObjs]):
           names += ['Tii','Tij']

        dlg = wx.MultiChoiceDialog(G2frame,'Select','Refinement controls',names)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            parms = []
            for x in sel:
                parms.append(names[x])
            wx.BeginBusyCursor()
            try:
                for rbObj in RBObjs:
                    if 'Origin' in parms:
                        rbObj['Orig'][1] = True
                    else:
                        rbObj['Orig'][1] = False
                    if 'Full Orient.' in parms:
                        rbObj['Orient'][1] = 'AV'
                    elif 'Orient. angle' in parms:
                        rbObj['Orient'][1] = 'A'
                    else:
                        rbObj['Orient'][1] = ' '
                    for i in range(len(rbObj['Torsions'])):
                        if 'Torsion '+str(i) in parms:
                            rbObj['Torsions'][i][1] = True
                        else:
                            rbObj['Torsions'][i][1] = False
                    if rbObj['ThermalMotion'][0] == 'Uiso':
                        if 'Uiso' in parms:
                           rbObj['ThermalMotion'][2][0] = True
                        else:
                           rbObj['ThermalMotion'][2][0] = False
                    if 'T' in rbObj['ThermalMotion'][0]:
                        if 'Tii' in parms:
                            rbObj['ThermalMotion'][2][0:2] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][0:2] = [False,False,False]
                        if 'Tij' in parms:
                            rbObj['ThermalMotion'][2][3:6] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][3:6] = [False,False,False]
                    if 'L' in rbObj['ThermalMotion'][0]:
                        if 'Lii' in parms:
                            rbObj['ThermalMotion'][2][6:9] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][6:9] = [False,False,False]
                        if 'Lij' in parms:
                            rbObj['ThermalMotion'][2][9:12] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][9:12] = [False,False,False]
                    if 'S' in rbObj['ThermalMotion'][0]:
                        if 'Sij' in parms:
                            rbObj['ThermalMotion'][2][12:20] = [True,True,True,True,True,True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][12:20] = [False,False,False,False,False,False,False,False]
            finally:
                wx.EndBusyCursor()
            FillRigidBodyGrid()

#### MC/SA routines ################################################################################
    def UpdateMCSA(Scroll=0):
        Indx = {}

        def OnPosRef(event):
            Obj = event.GetEventObject()
            model,item,ix = Indx[Obj.GetId()]
            model[item][1][ix] = Obj.GetValue()

        def OnPosVal(invalid,value,tc):
            G2plt.PlotStructure(G2frame,data)

        def OnPosRange(event):
            event.Skip()
            Obj = event.GetEventObject()
            model,item,ix = Indx[Obj.GetId()]
            Range = Obj.GetValue().split()
            try:
                rmin,rmax = [float(Range[i]) for i in range(2)]
                if rmin >= rmax:
                    raise ValueError
            except (ValueError,IndexError):
                rmin,rmax = model[item][2][ix]
            model[item][2][ix] = [rmin,rmax]
            Obj.SetValue('%.3f %.3f'%(rmin,rmax))

        def atomSizer(model):

            atomsizer = wx.FlexGridSizer(0,7,5,5)
            atomsizer.Add(wx.StaticText(G2frame.MCSA,-1,' Atom: '+model['name']+': '),0,WACV)
            for ix,item in enumerate(['x','y','z']):
                posRef = wx.CheckBox(G2frame.MCSA,-1,label=item+': ')
                posRef.SetValue(model['Pos'][1][ix])
                posRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                Indx[posRef.GetId()] = [model,'Pos',ix]
                atomsizer.Add(posRef,0,WACV)
                posVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,model['Pos'][0],ix,nDig=(10,4),OnLeave=OnPosVal)
                atomsizer.Add(posVal,0,WACV)
            atomsizer.Add((5,5),0)
            for ix,item in enumerate(['x','y','z']):
                atomsizer.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
                rmin,rmax = model['Pos'][2][ix]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                posRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[posRange.GetId()] = [model,'Pos',ix]
                posRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                posRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                atomsizer.Add(posRange,0,WACV)
            return atomsizer

        def rbSizer(model):

            def OnOrVar(event):
                Obj = event.GetEventObject()
                model = Indx[Obj.GetId()]
                model['Ovar'] = Obj.GetValue()

            def OnOriVal(event):
                event.Skip()
                Obj = event.GetEventObject()
                model,ix,ObjA,ObjV = Indx[Obj.GetId()]
                A = model['Ori'][0][0]
                V = model['Ori'][0][1:]
                if ix:
                    Anew = A
                    Vec = ObjV.GetValue().split()
                    try:
                        Vnew = [float(Vec[i]) for i in range(3)]
                    except ValueError:
                        Vnew = V
                else:
                    Vnew = V
                    try:
                        Anew = float(ObjA.GetValue())
                        if not Anew:    #==0.0!
                            Anew = 360.
                    except ValueError:
                        Anew = A
                Q = G2mth.AVdeg2Q(Anew,Vnew)
                A,V = G2mth.Q2AVdeg(Q)
                model['Ori'][0][0] = A
                model['Ori'][0][1:] = V
                if ix:
                    ObjV.SetValue('%.3f %.3f %.3f'%(V[0],V[1],V[2]))
                else:
                    ObjA.SetValue('%.5f'%(A))
                    ObjV.SetValue('%.3f %.3f %.3f'%(V[0],V[1],V[2]))
                G2plt.PlotStructure(G2frame,data)

            def OnMolCent(event):
                Obj = event.GetEventObject()
                model = Indx[Obj.GetId()]
                model['MolCent'][1] = Obj.GetValue()
                if model['MolCent'][1]:
                    G2mth.SetMolCent(model,RBData)
                G2plt.PlotStructure(G2frame,data)

            rbsizer = wx.BoxSizer(wx.VERTICAL)
            rbsizer1 = wx.FlexGridSizer(0,7,5,5)
            rbsizer1.Add(wx.StaticText(G2frame.MCSA,-1,model['Type']+': '+model['name']+': '),0)
            for ix,item in enumerate(['x','y','z']):
                posRef = wx.CheckBox(G2frame.MCSA,-1,label=item+': ')
                posRef.SetValue(model['Pos'][1][ix])
                posRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                Indx[posRef.GetId()] = [model,'Pos',ix]
                rbsizer1.Add(posRef,0,WACV)
                posVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,model['Pos'][0],ix,nDig=(10,4),OnLeave=OnPosVal)
                rbsizer1.Add(posVal,0,WACV)
            molcent = wx.CheckBox(G2frame.MCSA,-1,label=' Use mol. center? ')
            molcent.SetValue(model['MolCent'][1])
            molcent.Bind(wx.EVT_CHECKBOX,OnMolCent)
            Indx[molcent.GetId()] = model
            rbsizer1.Add(molcent,0,WACV)
            for ix,item in enumerate(['x','y','z']):
                rbsizer1.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
                rmin,rmax = model['Pos'][2][ix]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                posRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[posRange.GetId()] = [model,'Pos',ix]
                posRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                posRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                rbsizer1.Add(posRange,0,WACV)

            rbsizer2 = wx.FlexGridSizer(0,6,5,5)
            Ori = model['Ori'][0]
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oa: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            angVal = wx.TextCtrl(G2frame.MCSA,-1,'%.5f'%(Ori[0]),style=wx.TE_PROCESS_ENTER)
            angVal.Bind(wx.EVT_TEXT_ENTER,OnOriVal)
            angVal.Bind(wx.EVT_KILL_FOCUS,OnOriVal)
            rbsizer2.Add(angVal,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            vecVal = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f %.3f'%(Ori[1],Ori[2],Ori[3]),style=wx.TE_PROCESS_ENTER)
            vecVal.Bind(wx.EVT_TEXT_ENTER,OnOriVal)
            vecVal.Bind(wx.EVT_KILL_FOCUS,OnOriVal)
            Indx[angVal.GetId()] = [model,0,angVal,vecVal]
            Indx[vecVal.GetId()] = [model,1,angVal,vecVal]
            rbsizer2.Add(vecVal,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,' Vary? '),0,WACV)
            choice = [' ','A','AV']
            orvar = wx.ComboBox(G2frame.MCSA,-1,value=model['Ovar'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            orvar.Bind(wx.EVT_COMBOBOX, OnOrVar)
            Indx[orvar.GetId()] = model
            rbsizer2.Add(orvar,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,' Range: Oa: '),0,WACV)
            Rge = model['Ori'][2]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            angRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(Rge[0][0],Rge[0][1]),style=wx.TE_PROCESS_ENTER)
            Indx[angRange.GetId()] = [model,'Ori',0]
            angRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
            angRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
            rbsizer2.Add(angRange,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
            for io,item in enumerate(['Oi','Oj','Ok']):
                rmin,rmax = Rge[io+1]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                vecRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f '%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[vecRange.GetId()] = [model,'Ori',io+1]
                vecRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                vecRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                rbsizer2.Add(vecRange,0,WACV)
            rbsizer.Add(rbsizer1)
            rbsizer.Add(rbsizer2)
            if model['Type'] == 'Residue':
                try:
                    atNames = RBData['Residue'][model['RBId']]['atNames']
                    rbsizer.Add(wx.StaticText(G2frame.MCSA,-1,'Torsions:'),0)
                    rbsizer3 = wx.FlexGridSizer(0,8,5,5)
                    for it,tor in enumerate(model['Tor'][0]):
                        iBeg,iFin = RBData['Residue'][model['RBId']]['rbSeq'][it][:2]
                        name = atNames[iBeg]+'-'+atNames[iFin]
                        torRef = wx.CheckBox(G2frame.MCSA,label=' %s: '%(name))
                        torRef.SetValue(model['Tor'][1][it])
                        torRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                        Indx[torRef.GetId()] = [model,'Tor',it]
                        rbsizer3.Add(torRef,0,WACV)
                        torVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,model['Tor'][0],it,nDig=(10,4),OnLeave=OnPosVal)
                        rbsizer3.Add(torVal,0,WACV)
                        rbsizer3.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
                        rmin,rmax = model['Tor'][2][it]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
                        torRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                        Indx[torRange.GetId()] = [model,'Tor',it]
                        torRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                        torRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                        rbsizer3.Add(torRange,0,WACV)
                    rbsizer.Add(rbsizer3)
                except KeyError:    #Missing RB - clear all away!
                    data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
                    wx.CallAfter(UpdateMCSA)
            return rbsizer

        def MDSizer(POData):

            def OnPORef(event):
                POData['Coef'][1] = poRef.GetValue()

            def OnPORange(event):
                event.Skip()
                Range = poRange.GetValue().split()
                try:
                    rmin,rmax = [float(Range[i]) for i in range(2)]
                    if 0. < rmin < rmax:
                        pass
                    else:
                        raise ValueError
                except (ValueError,IndexError):
                    rmin,rmax = POData['Coef'][2]
                POData['Coef'][2] = [rmin,rmax]
                poRange.SetValue('%.3f %.3f'%(rmin,rmax))

            def OnPOAxis(event):
                event.Skip()
                Saxis = poAxis.GetValue().split()
                try:
                    hkl = [int(Saxis[i]) for i in range(3)]
                except (ValueError,IndexError):
                    hkl = POData['axis']
                if not np.any(np.array(hkl)):
                    hkl = POData['axis']
                POData['axis'] = hkl
                h,k,l = hkl
                poAxis.SetValue('%3d %3d %3d'%(h,k,l))

            poSizer = wx.BoxSizer(wx.HORIZONTAL)
            poRef = wx.CheckBox(G2frame.MCSA,-1,label=' March-Dollase ratio: ')
            poRef.SetValue(POData['Coef'][1])
            poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
            poSizer.Add(poRef,0,WACV)
            poVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,POData['Coef'],0,nDig=(10,3),xmin=0.)
            poSizer.Add(poVal,0,WACV)
            poSizer.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
            rmin,rmax = POData['Coef'][2]
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            poRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
            poRange.Bind(wx.EVT_TEXT_ENTER,OnPORange)
            poRange.Bind(wx.EVT_KILL_FOCUS,OnPORange)
            poSizer.Add(poRange,0,WACV)
            poSizer.Add(wx.StaticText(G2frame.MCSA,-1,' Unique axis, H K L: '),0,WACV)
            h,k,l = POData['axis']
#            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),xmin=0.01,xmax=4.0)
            poAxis = wx.TextCtrl(G2frame.MCSA,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
            poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOAxis)
            poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOAxis)
            poSizer.Add(poAxis,0,WACV)
            return poSizer

        def ResultsSizer(Results):

            def OnCellChange(event):
                r,c = event.GetRow(),event.GetCol()
                if c == 0:
                    for row in range(resultsGrid.GetNumberRows()):
                        resultsTable.SetValue(row,c,False)
                        Results[row][0] = False
                    result = Results[r]
                    Models = data['MCSA']['Models']
                    SetSolution(result,Models)
                    Results[r][0] = True
                    resultsTable.SetValue(r,0,True)
#                    resultsGrid.ForceRefresh()
                    G2plt.PlotStructure(G2frame,data)
                    wx.CallAfter(UpdateMCSA,G2frame.MCSA.GetScrollPos(wx.VERTICAL))
                elif c == 1:
                    if Results[r][1]:
                        Results[r][1] = False
                    else:
                        Results[r][1] = True
                    resultsTable.SetValue(r,c,Results[r][1])
                    resultsGrid.ForceRefresh()

            resultsSizer = wx.BoxSizer(wx.VERTICAL)
            maxVary = 0
            resultVals = []
            for result in Results:
                maxVary = max(maxVary,len(result[-1]))
                resultVals.append(result[:-1])
            rowLabels = []
            for i in range(len(Results)): rowLabels.append(str(i))
            colLabels = ['Select','Keep','Residual','Tmin',]
            for item in result[-1]: colLabels.append(item)   #from last result from for loop above
            Types = [wg.GRID_VALUE_BOOL,wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,4',
                wg.GRID_VALUE_FLOAT+':10,4',]+maxVary*[wg.GRID_VALUE_FLOAT+':10,5',]
            resultsTable = G2G.Table(resultVals,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            resultsGrid = G2G.GSGrid(G2frame.MCSA)
            resultsGrid.SetTable(resultsTable, True)
            resultsGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, OnCellChange)
            resultsGrid.AutoSizeColumns(True)
            resultsGrid.SetMargins(0,0)
            for r in range(resultsGrid.GetNumberRows()):
                for c in range(resultsGrid.GetNumberCols()):
                    if c in [0,1]:
                        resultsGrid.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        resultsGrid.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
#            resultsGrid.SetScrollRate(0,0)
            resultsSizer.Add(resultsGrid,0,wx.EXPAND)
            return resultsSizer

        def OnSelect(event):
            rbId = rbids[select.GetSelection()]
            wx.CallLater(100,RepaintRBInfo,rbId)

        def RepaintRBInfo(rbId,Scroll=0):
            oldFocus = wx.Window.FindFocus()
            if 'phoenix' in wx.version():
                G2frame.bottomSizer.Clear(True)
            else:
                G2frame.bottomSizer.DeleteWindows()
            Indx.clear()
            rbObj = data['MCSA']['Models'][rbId]
            G2frame.bottomSizer.Insert(0,rbSizer(rbObj))
            mainSizer.Layout()
            G2frame.dataWindow.Refresh()
            G2frame.dataWindow.SendSizeEvent()
            wx.CallAfter(oldFocus.SetFocus)

        def OnShoLabels(event):
            data['MCSA']['showLabels'] = not data['MCSA']['showLabels']
            G2plt.PlotStructure(G2frame,data)

        # UpdateMCSA executable code starts here
        if G2frame.MCSA.GetSizer(): G2frame.MCSA.GetSizer().Clear(True)
        #patch
        data['MCSA']['showLabels'] = data['MCSA'].get('showLabels',False)
        #end patch
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        general = data['General']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies')
        if not Id:
            return
        RBData = G2frame.GPXtree.GetItemPyData(Id)
        Indx = {}
#        atomStyle = 'balls & sticks'
#        if 'macro' in general['Type']:
#            atomStyle = 'sticks'
        G2frame.GetStatusBar().SetStatusText('',1)
        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Monte Carlo/Simulated Annealing for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if not data['MCSA']['Models']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,-1,'No MC/SA models:'),0)
            mainSizer.Add((5,5),0)
        else:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,label='MC/SA models:'))
            mainSizer.Add((5,5),0)
            rbNames = []
            rbids = []
            for im,model in enumerate(data['MCSA']['Models']):
                if model['Type'] == 'MD':
                    mainSizer.Add(MDSizer(model))
                elif model['Type'] == 'Atom':
                    Asizer = atomSizer(model)
                    mainSizer.Add(Asizer)
                else:
                    rbNames.append(model['name'])
                    rbids.append(im)
            G2G.HorizontalLine(mainSizer,G2frame.MCSA)
            if len(rbNames):
                rbName = rbNames[0]
                select = wx.ListBox(G2frame.MCSA,choices=rbNames,style=wx.LB_SINGLE,size=(-1,65))
                select.SetSelection(rbNames.index(rbName))
                select.SetFirstItem(rbNames.index(rbName))
                select.Bind(wx.EVT_LISTBOX,OnSelect)
                mainSizer.Add(select,0)
                G2frame.bottomSizer = wx.BoxSizer(wx.VERTICAL)
                G2frame.bottomSizer.Add(rbSizer(data['MCSA']['Models'][rbids[0]]))
                mainSizer.Add(G2frame.bottomSizer)

        mainSizer.Add((5,5),0)
        bottomSizer = wx.BoxSizer(wx.HORIZONTAL)
        resStr = 'MC/SA results:  '
        if not data['MCSA']['Results']:
            resStr = 'No'+resStr
        bottomSizer.Add(wx.StaticText(G2frame.MCSA,-1,resStr),0,WACV)
        shoLabels = wx.CheckBox(G2frame.MCSA,label=' Show atom labels? ')
        shoLabels.SetValue(data['MCSA']['showLabels'])
        shoLabels.Bind(wx.EVT_CHECKBOX,OnShoLabels)
        bottomSizer.Add(shoLabels,0,WACV)
        mainSizer.Add(bottomSizer)
        mainSizer.Add((5,5),0)
        if data['MCSA']['Results']:
            Results = data['MCSA']['Results']
            mainSizer.Add(ResultsSizer(Results),0,wx.EXPAND)

        SetPhaseWindow(G2frame.MCSA,mainSizer)

    def SetSolution(result,Models):
        for key,val in zip(result[-1],result[4:-1]):
            vals = key.split(':')
            nObj,name = int(vals[0]),vals[1]
            if 'A' in name:
                ind = ['Ax','Ay','Az'].index(name)
                Models[nObj]['Pos'][0][ind] = val
            elif 'Q' in name:
                ind = ['Qa','Qi','Qj','Qk'].index(name)
                Models[nObj]['Ori'][0][ind] = val
            elif 'P' in name:
                ind = ['Px','Py','Pz'].index(name)
                Models[nObj]['Pos'][0][ind] = val
            elif 'T' in name:
                tnum = int(name.split('Tor')[1])
                Models[nObj]['Tor'][0][tnum] = val
            else:       #March Dollase
                Models[0]['Coef'][0] = val

    def OnRunMultiMCSA(event):
        RunMCSA('multi')

    def OnRunSingleMCSA(event):
        RunMCSA('single')

    def RunMCSA(process):
        generalData = data['General']
        mcsaControls = generalData['MCSA controls']
        reflName = mcsaControls['Data source']
        phaseName = generalData['Name']
        MCSAdata = data['MCSA']
        saveResult = []
        for result in MCSAdata['Results']:
            if result[1]:       #keep?
                saveResult.append(result)
        MCSAdata['Results'] = saveResult
        covData = {}
        if 'PWDR' in reflName:
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Reflection Lists'))
            try:        #patch for old reflection data
                reflData = reflSets[phaseName]['RefList']
            except TypeError:
                reflData = reflSets[phaseName]
            reflType = 'PWDR'
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
            try:
                reflData = G2frame.GPXtree.GetItemPyData(PatternId)[1]['RefList']
            except TypeError:
                reflData = G2frame.GPXtree.GetItemPyData(PatternId)[1]
            reflType = 'HKLF'
        elif reflName == 'Pawley reflections':
            reflData = data['Pawley ref']
            covData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Covariance'))
            reflType = 'Pawley'
        else:
            print ('**** ERROR - No data defined for MC/SA run')
            return
        print ('MC/SA run:')
        print ('Reflection type:'+reflType+' Total No. reflections: %d'%len(reflData))
        RBdata = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        MCSAmodels = MCSAdata['Models']
        if not len(MCSAmodels):
            print ('**** ERROR - no models defined for MC/SA run****')
            return
        time1 = time.time()
        nprocs = GSASIIpath.GetConfigValue('Multiprocessing_cores',0)
        if process == 'single' or not nprocs:
            pgbar = wx.ProgressDialog('MC/SA','Residual Rcf =',101,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
            screenSize = wx.ClientDisplayRect()
            Size = pgbar.GetSize()
            if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
                pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
                pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        else:
            pgbar = None
        try:
            tsf = 0.
            nCyc = mcsaControls['Cycles']
            if process == 'single' or not nprocs:
                for i in range(nCyc):
                    pgbar.SetTitle('MC/SA run '+str(i+1)+' of '+str(nCyc))
                    Result,tsum,nsum,rcov = G2mth.mcsaSearch(data,RBdata,reflType,reflData,covData,pgbar)
                    MCSAdata['Results'].append(Result)
                    print(' MC/SA run completed: %d residual: %.3f%% SFcalc time: %.2fs Nsfcalc: %d'%(i,100*Result[2],tsum,nsum))
                    tsf += tsum
                print (' Total SF time: %.2fs'%(tsf))
                XY = np.mgrid[0:rcov.shape[0],0:rcov.shape[1]]
                G2plt.PlotXYZ(G2frame,XY,rcov,labelX='ref No.',labelY='ref No.',newPlot=False,
                    Title='Reflection covariance matrix',zrange=[-1.,1.],color='RdYlGn')
            else:
                Results,sftime,numsf = G2mth.MPmcsaSearch(nCyc,data,RBdata,reflType,reflData,covData,nprocs)
                MCSAdata['Results'] += Results   #+= to  any saved ones
                print (' Total SF time: %.2fs MC/SA run time: %.2fs Nsfcalc: %d'%(sftime,time.time()-time1,numsf))
        finally:
            if process == 'single' or not nprocs:
                pgbar.Destroy()
        MCSAdata['Results'] = G2mth.sortArray(MCSAdata['Results'],2,reverse=False)
        MCSAdata['Results'][0][0] = True
        SetSolution(MCSAdata['Results'][0],data['MCSA']['Models'])
        G2frame.phaseDisplay.SetFocus()
        Page = G2frame.phaseDisplay.FindPage('MC/SA')
        G2frame.phaseDisplay.SetSelection(Page)
        G2plt.PlotStructure(G2frame,data)
        wx.CallAfter(UpdateMCSA)

    def OnMCSAaddAtom(event):
        dlg = G2elemGUI.PickElement(G2frame)
        if dlg.ShowModal() == wx.ID_OK:
            El = dlg.Elem.strip()
            Info = G2elem.GetAtomInfo(El)
        dlg.Destroy()

        atom = {'Type':'Atom','atType':El,'Pos':[[0.,0.,0.],
            [False,False,False],[[0.,1.],[0.,1.],[0.,1.]]],
            'name':El+'('+str(len(data['MCSA']['Models']))+')'}
        data['MCSA']['Models'].append(atom)
        data['MCSA']['AtInfo'][El] = [Info['Drad'],Info['Color']]
        G2plt.PlotStructure(G2frame,data)
        UpdateMCSA()

    def OnMCSAaddRB(event):
        rbData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbNames = {}
        for rbVec in rbData['Vector']:
            if rbVec != 'AtInfo':
                rbNames[rbData['Vector'][rbVec]['RBname']] = ['Vector',rbVec]
        for rbRes in rbData['Residue']:
            if rbRes != 'AtInfo':
                rbNames[rbData['Residue'][rbRes]['RBname']] = ['Residue',rbRes]
        if not rbNames:
            print ('**** ERROR - no rigid bodies defined ****')
            return
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Rigid body',list(rbNames.keys()))
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            rbname = list(rbNames.keys())[sel]
            rbType,rbId = rbNames[rbname]
            RB = rbData[rbType][rbId]
        body = {'name':RB['RBname']+'('+str(len(data['MCSA']['Models']))+')','RBId':rbId,'Type':rbType,
            'Pos':[[0.,0.,0.],[False,False,False],[[0.,1.],[0.,1.],[0.,1.]]],'Ovar':'','MolCent':[[0.,0.,0.],False],
            'Ori':[[180.,0.,0.,1.],[False,False,False,False],[[0.,360.],[-1.,1.],[-1.,1.],[-1.,1.]]]}
        if rbType == 'Residue':
            body['Tor'] = [[],[],[]]
            for i,tor in enumerate(RB['rbSeq']):
                body['Tor'][0].append(0.0)
                body['Tor'][1].append(False)
                body['Tor'][2].append([0.,360.])
        data['MCSA']['Models'].append(body)
        data['MCSA']['rbData'] = rbData
        data['MCSA']['AtInfo'].update(rbData[rbType]['AtInfo'])
        UpdateMCSA()
        wx.CallAfter(G2plt.PlotStructure,G2frame,data)

    def OnMCSAclear(event):
        data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
        G2plt.PlotStructure(G2frame,data)
        UpdateMCSA()

    def OnMCSAmove(event):
        general = data['General']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        xyz,aTypes = G2mth.UpdateMCSAxyz(Bmat,data['MCSA'])
        for iat,atype in enumerate(aTypes):
            x,y,z = xyz[iat]
            AtomAdd(x,y,z,atype,Name=atype+'(%d)'%(iat+1))
        G2plt.PlotStructure(G2frame,data)

    def OnClearResults(event):
        data['MCSA']['Results'] = []
        UpdateMCSA()

#### Pawley routines ################################################################################
    def FillPawleyReflectionsGrid():

        def onRefineDClick(event):
            '''Called after a double-click on a cell label'''
            c =  event.GetCol()
            if c == 5:     #refine column label: just select it (& redisplay)
                G2frame.PawleyRefl.ClearSelection()
                G2frame.PawleyRefl.SelectCol(c,True)
                choice = ['Y - vary all','N - vary none',]
                dlg = wx.SingleChoiceDialog(G2frame,'Select refinement option',
                    'Refinement controls',choice)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    if sel == 0:
                        for row in range(G2frame.PawleyRefl.GetNumberRows()): PawleyPeaks[row][c]=True
                    else:
                        for row in range(G2frame.PawleyRefl.GetNumberRows()): PawleyPeaks[row][c]=False
                wx.CallAfter(FillPawleyReflectionsGrid)

        def KeyEditPawleyGrid(event):
            colList = G2frame.PawleyRefl.GetSelectedCols()
            rowList = G2frame.PawleyRefl.GetSelectedRows()
            PawleyPeaks = data['Pawley ref']
            if event.GetKeyCode() == wx.WXK_RETURN:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_CONTROL:
                event.Skip(True)
            elif event.GetKeyCode() == wx.WXK_SHIFT:
                event.Skip(True)
            elif colList:
                G2frame.PawleyRefl.ClearSelection()
                key = event.GetKeyCode()
                for col in colList:
                    if PawleyTable.GetTypeName(0,col) == wg.GRID_VALUE_BOOL:
                        if key == 89: #'Y'
                            for row in range(PawleyTable.GetNumberRows()): PawleyPeaks[row][col]=True
                        elif key == 78:  #'N'
                            for row in range(PawleyTable.GetNumberRows()): PawleyPeaks[row][col]=False
                        FillPawleyReflectionsGrid()
            elif rowList:
                if event.GetKeyCode() == wx.WXK_DELETE:
                    rowList.reverse()
                    for row in rowList:
                        del(PawleyPeaks[row])
                    FillPawleyReflectionsGrid()

        # FillPawleyReflectionsGrid executable starts here
        G2frame.GetStatusBar().SetStatusText('To delete a Pawley reflection: select row & press Delete',1)
        if G2frame.PawleyRefl in G2frame.phaseDisplay.gridList:
            G2frame.phaseDisplay.gridList.remove(G2frame.PawleyRefl)
        oldSizer = PawleyRefList.GetSizer()
        if oldSizer: oldSizer.Clear(True)
        generalData = data['General']
        PawleyPeaks = data['Pawley ref']

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        if len(PawleyPeaks) and generalData['doPawley']:
            lbl= f"Pawley reflections for {data['General']['Name']!r}"[:60]
        else:
            lbl= f"There are no Pawley reflections for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)

        rowLabels = []
        if len(PawleyPeaks) and generalData['doPawley']:
            for i in range(len(PawleyPeaks)): rowLabels.append(str(i))
            if generalData['Modulated']:
                colLabels = ['h','k','l','m','mul','d','refine','Fsq(hkl)','sig(Fsq)']
                Types = 5*[wg.GRID_VALUE_LONG,]+[wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,]+ \
                    2*[wg.GRID_VALUE_FLOAT+':10,2',]
                pos = [6,7]
            else:
                colLabels = ['h','k','l','mul','d','refine','Fsq(hkl)','sig(Fsq)']
                Types = 4*[wg.GRID_VALUE_LONG,]+[wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,]+ \
                    2*[wg.GRID_VALUE_FLOAT+':10,2',]
                pos = [5,6]
            PawleyTable = G2G.Table(PawleyPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            G2frame.PawleyRefl = G2G.GSGrid(PawleyRefList)
            G2frame.PawleyRefl.SetTable(PawleyTable, True)
            G2frame.PawleyRefl.Unbind(wx.EVT_KEY_DOWN)
            G2frame.PawleyRefl.Bind(wx.EVT_KEY_DOWN, KeyEditPawleyGrid)
            G2frame.PawleyRefl.Unbind(wg.EVT_GRID_LABEL_LEFT_DCLICK)
            G2frame.PawleyRefl.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, onRefineDClick)
            for r in range(G2frame.PawleyRefl.GetNumberRows()):
                for c in range(G2frame.PawleyRefl.GetNumberCols()):
                    if c in pos:
                        G2frame.PawleyRefl.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        G2frame.PawleyRefl.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.PawleyRefl.SetMargins(0,0)
            G2frame.PawleyRefl.AutoSizeColumns(False)
            G2frame.PawleyRefl.SetRowLabelSize(50)
            G2frame.PawleyRefl.SetScrollRate(10,10) # allow grid to scroll
            mainSizer.Add(G2frame.PawleyRefl,1,wx.EXPAND)
            for r in range(G2frame.PawleyRefl.GetNumberRows()):
                try:
                    if float(G2frame.PawleyRefl.GetCellValue(r,6)) < 0:
                        G2frame.PawleyRefl.SetCellBackgroundColour(r,6,wx.RED)
                except:
                    pass
        else:
            msg = (
'''Pawley refinement has not yet been setup. Use the
Operations->"Pawley setttings" menu command to change this.
(Or, if Pawley settings have already been set on the General
tab, use Operations->"Pawley create")''')
            mainSizer.Add(wx.StaticText(PawleyRefList,label=msg),0,wx.ALIGN_CENTER)
        SetPhaseWindow(PawleyRefList,mainSizer)

    def OnPawleySet(event):
        '''Open dialog to set Pawley parameters and optionally recompute reflections.
        This is called from the Phase/Pawley Reflections "Pawley Settings"
        menu command. These settings are also available on the Phase/General tab.
        '''
        def DisablePawleyOpts(*args):
            'dis-/enable Pawley options'
            for c in PawleyCtrlsList:
                c.Enable(generalData['doPawley'])

        PawleyCtrlsList = []
        dmin,dmax,nhist,lbl,chist = getPawleydRange(G2frame,data)
        generalData = data['General']
        prevPawleySetting = generalData['doPawley']
        generalData['doPawley'] = True  # make Pawley extraction the default if window is opened
        if nhist == 0:  # no data, no Pawley (probably can't happen here)
            generalData['doPawley'] = False
        else:
            # force limits on dmin & dmax
            generalData['Pawley dmax'] = min(generalData['Pawley dmax'],dmax)
            generalData['Pawley dmin'] = max(generalData['Pawley dmin'],dmin)
        startDmin = generalData['Pawley dmin']
        genDlg = wx.Dialog(G2frame,title='Set Pawley Parameters',
                    style=wx.DEFAULT_DIALOG_STYLE)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(genDlg,label='Set Pawley Extraction Parameters for phase '+generalData['Name']))
        mainSizer.Add([5,10])
        pawlRef = G2G.G2CheckBoxFrontLbl(genDlg,' Do Pawley refinement?: ',generalData,'doPawley',DisablePawleyOpts)
        mainSizer.Add(pawlRef)
        mainSizer.Add(wx.StaticText(genDlg,label=lbl),0,wx.BOTTOM,10)
        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label='   Pawley dmin: '),0,WACV)
        #def d2Q(*a,**kw):
        #    temp['Qmax'] = 2 * math.pi / generalData['Pawley dmin']
        #    pawlQVal.SetValue(temp['Qmax'])
        pawlVal = G2G.ValidatedTxtCtrl(genDlg,generalData,'Pawley dmin',
            xmin=dmin,xmax=20.,nDig=(10,5),typeHint=float)
#            xmin=dmin,xmax=20.,nDig=(10,5),typeHint=float,OnLeave=d2Q)
        PawleyCtrlsList.append(pawlVal)
        pawleySizer.Add(pawlVal,0,WACV)
        #pawleySizer.Add(wx.StaticText(genDlg,label='   Qmax: '),0,WACV)
        #temp = {'Qmax':2 * math.pi / generalData['Pawley dmin']}
        #def Q2D(*args,**kw):
        #    generalData['Pawley dmin'] = 2 * math.pi / temp['Qmax']
        #    pawlVal.SetValue(generalData['Pawley dmin'])
        #pawlQVal = G2G.ValidatedTxtCtrl(genDlg,temp,'Qmax',
        #    xmin=0.314,xmax=25.,nDig=(10,5),typeHint=float,OnLeave=Q2D)
        #PawleyCtrlsList.append(pawlQVal)
        #pawleySizer.Add(pawlQVal,0,WACV)
        mainSizer.Add(pawleySizer)

        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label='   Pawley dmax: '),0,WACV)
        pawlVal = G2G.ValidatedTxtCtrl(genDlg,generalData,'Pawley dmax',
            xmin=2.,xmax=dmax,nDig=(10,5),typeHint=float)
        PawleyCtrlsList.append(pawlVal)
        pawleySizer.Add(pawlVal,0,WACV)
        mainSizer.Add(pawleySizer)

        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label=' Pawley neg. wt.: '),0,WACV)
        pawlNegWt = G2G.ValidatedTxtCtrl(genDlg,generalData,'Pawley neg wt',
            xmin=0.,xmax=1.,nDig=(10,4),typeHint=float)
        pawleySizer.Add(pawlNegWt,0,WACV)
        PawleyCtrlsList.append(pawlNegWt)
        mainSizer.Add(pawleySizer)

        # make OK button
        def OnOK(event): genDlg.EndModal(wx.ID_OK)
        mainSizer.Add([5,5])
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(genDlg, wx.ID_OK)
        btn.Bind(wx.EVT_BUTTON, OnOK)
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(genDlg, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)

        genDlg.SetSizer(mainSizer)
        mainSizer.Fit(genDlg)
        genDlg.CenterOnParent()
        res = genDlg.ShowModal()
        genDlg.Destroy()
        if res == wx.ID_NO: return
        # If Pawley is on, turn off Le Bail settings (since they will
        # be hidden)
        if generalData['doPawley']:
            Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.GPXtree.root, 'Controls'))
            Controls['newLeBail'] = False
            for h in data['Histograms']:
                data['Histograms'][h]['LeBail'] = False
        # ask to generate the reflections if the extraction setting or dmin has changed
        if generalData['doPawley'] and res == wx.ID_OK and (
                not prevPawleySetting or startDmin != generalData['Pawley dmin']):
            dlg = wx.MessageDialog(G2frame,'Do you want to generate the Pawley reflections with these settings? ("Pawley create" command)','Initialize Pawley?',
                wx.YES_NO | wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_NO:
                    wx.CallAfter(FillPawleyReflectionsGrid)
                    return
            finally:
                dlg.Destroy()
            wx.CallAfter(OnPawleyLoad,event)
        else:
            wx.CallAfter(FillPawleyReflectionsGrid)

    def OnPawleyLoad(event):
        generalData = data['General']
        histograms = data['Histograms'].keys()
        cell = generalData['Cell'][1:7]
        A = G2lat.cell2A(cell)
        SGData = generalData['SGData']
        dmin = generalData['Pawley dmin']
        dmax = generalData['Pawley dmax']
        for hist in histograms:
            if 'PWDR' in hist[:4]:
                Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
                inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
                    G2frame,Id, 'Instrument Parameters'))[0]
                limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Limits'))
                Tmin = G2lat.Dsp2pos(inst,dmin)
                if 'T' in inst['Type'][0]:
                    limits[1][0] = max(limits[0][0],Tmin)
                else:
                    limits[1][1] = min(limits[0][1],Tmin)
        PawleyPeaks = []
        HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
        if generalData['Modulated']:
            Vec,x,maxH = generalData['SuperVec']
            SSGData = G2spc.SSpcGroup(SGData,generalData['SuperSg'])[1]
            wx.BeginBusyCursor()
            try:
                HKLd = G2lat.GenSSHLaue(dmin,SGData,SSGData,Vec,maxH,A)
                for h,k,l,m,d in HKLd:
                    if d > dmax:
                        continue
                    ext,mul = G2spc.GenHKLf([h,k,l],SGData)[:2]
                    if m or not ext:
                        mul *= 2        #for powder multiplicity
                        PawleyPeaks.append([h,k,l,m,mul,d,False,100.0,1.0])
                PawleyPeaks = G2mth.sortArray(PawleyPeaks,5,reverse=True)
            finally:
                wx.EndBusyCursor()
        else:
            wx.BeginBusyCursor()
            try:
                for h,k,l,d in HKLd:
                    if d > dmax:
                        continue
                    ext,mul = G2spc.GenHKLf([h,k,l],SGData)[:2]
                    if not ext:
                        mul *= 2        #for powder multiplicity
                        PawleyPeaks.append([h,k,l,mul,d,False,100.0,1.0])
                PawleyPeaks = G2mth.sortArray(PawleyPeaks,4,reverse=True)
            finally:
                wx.EndBusyCursor()
        data['Pawley ref'] = PawleyPeaks

        dlg = wx.MessageDialog(G2frame,'Do you want to initialize Pawley reflection intensities? ("Pawley estimate" command)','Initialize Pawley?',
                wx.YES_NO | wx.ICON_QUESTION)
        try:
            result = dlg.ShowModal()
            if result == wx.ID_NO:
                wx.CallAfter(FillPawleyReflectionsGrid)
                return
        finally:
                dlg.Destroy()
        wx.CallAfter(OnPawleyEstimate,event)

    def OnPawleyEstimate(event):
        #Algorithm thanks to James Hester
        try:
            Refs = data['Pawley ref']
            Histograms = data['Histograms']
        except KeyError:
            G2frame.ErrorDialog('Pawley estimate','No histograms defined for this phase')
            return
        Vst = 1.0/data['General']['Cell'][7]     #Get volume
        generalData = data['General']
        im = 0
        if generalData['Modulated']:
            im = 1
        HistoNames = list(filter(lambda a:Histograms[a]['Use']==True,list(Histograms.keys())))
        if not len(HistoNames):
            G2frame.ErrorDialog('Pawley estimate','No histograms defined for this phase')
            return
        PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistoNames[0])       #only use 1st histogram
        xdata = G2frame.GPXtree.GetItemPyData(PatternId)[1]
        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Instrument Parameters'))[0]
        Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Sample Parameters'))
        const = 9.e-2/(np.pi*Sample['Gonio. radius'])                  #shifts in microns
        gconst = 2.35482 # sqrt(8 ln 2)
        dx = (xdata[0][1]-xdata[0][0])*20.              #this ian approximation - not correct, but CW seems to be needed
        cw = np.diff(xdata[0])
        cw = np.append(cw,cw[-1])
        gconst *= dx

        wx.BeginBusyCursor()
        try:
            for ref in Refs:
                pos = G2lat.Dsp2pos(Inst,ref[4+im])
                if 'Bragg' in Sample['Type']:
                    pos -= const*(4.*Sample['Shift'][0]*cosd(pos/2.0)+ \
                        Sample['Transparency'][0]*sind(pos)*100.0)            #trans(=1/mueff) in cm
                elif 'E' not in Inst['Type'][0]:               #Debye-Scherrer - simple but maybe not right - +Layer Disp from DData?
                    pos -= const*(Sample['DisplaceX'][0]*cosd(pos)+Sample['DisplaceY'][0]*sind(pos))
                indx = np.searchsorted(xdata[0],pos)
                try:
                    FWHM = max(0.001,G2pwd.getFWHM(pos,Inst))
                    # We want to estimate Pawley F^2 as a drop-in replacement for F^2 calculated by the structural
                    # routines, which use Icorr * F^2 * peak profile, where peak profile has an area of 1.  So
                    # we multiply the observed peak height by sqrt(8 ln 2)/(FWHM*sqrt(pi)) to determine the value of Icorr*F^2
                    # then divide by Icorr to get F^2.
                    ref[6+im] = (xdata[1][indx]-xdata[4][indx])*FWHM*np.sqrt(np.pi)  #Area of Gaussian is height * FWHM * sqrt(pi)
                    if 'E' not in Inst['Type'][0]:
                        if 'C' in Inst['Type'][0]:
                            Lorenz = 1./(2.*sind(xdata[0][indx]/2.)**2*cosd(xdata[0][indx]/2.))           #Lorentz correction
                        else:
                            Lorenz = ref[4+im]**4
                        pola = 1.0
                        if 'X' in Inst['Type']:
                            pola,dIdPola = G2pwd.Polarization(Inst['Polariz.'][1],xdata[0][indx],Inst['Azimuth'][1])
                        else:
                            pola = 1.0
                    # Include histo scale and volume in calculation
                        ref[6+im] /= (Sample['Scale'][0] * Vst * Lorenz * pola * ref[3+im])
                    else:
                        ref[6+im] /= (Sample['Scale'][0] * ref[3+im]*cw[indx])
                except IndexError:
                    ref[6+im] = 0.0
        finally:
            wx.EndBusyCursor()
        FillPawleyReflectionsGrid()

    def OnPawleyUpdate(event):
        '''This is the place for any reflection modification trick
        Patterson squared, etc.
        '''
        try:
            Refs = data['Pawley ref']
            Histograms = data['Histograms']
        except KeyError:
            G2frame.ErrorDialog('Pawley update','No histograms defined for this phase')
            return
        HistoNames = list(Histograms.keys())
        PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistoNames[0])
        refData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,  \
            PatternId,'Reflection Lists'))[PhaseName]['RefList']
        im = 0
        if data['General']['Modulated']:
            im = 1
        Inv = data['General']['SGData']['SGInv']
        mult = 0.5
        if Inv:
            mult = 0.3
        wx.BeginBusyCursor()
        try:
            for iref,ref in enumerate(Refs):
                try:
                    if ref[6+im] < 0.:
                        ref[6+im] *= -mult
                        refData[iref][8+im] *= -mult
                        refData[iref][9+im] *= -mult
                        ref[5+im] = False
                        ref[7+im] = 1.0
                except IndexError:
                    print ('skipped',ref)
                    pass
        finally:
            wx.EndBusyCursor()
        wx.CallAfter(FillPawleyReflectionsGrid)

    def OnPawleySelAll(event):
        refcol = [G2frame.PawleyRefl.GetColLabelValue(c) for c in range(G2frame.PawleyRefl.GetNumberCols())].index('refine')
        for r in range(G2frame.PawleyRefl.GetNumberRows()):
            G2frame.PawleyRefl.GetTable().SetValue(r,refcol,True)

        G2frame.PawleyRefl.ForceRefresh()
    def OnPawleySelNone(event):
        refcol = [G2frame.PawleyRefl.GetColLabelValue(c) for c in range(G2frame.PawleyRefl.GetNumberCols())].index('refine')
        for r in range(G2frame.PawleyRefl.GetNumberRows()):
            G2frame.PawleyRefl.GetTable().SetValue(r,refcol,False)
        G2frame.PawleyRefl.ForceRefresh()

    def OnPawleyToggle(event):

        refcol = [G2frame.PawleyRefl.GetColLabelValue(c) for c in range(G2frame.PawleyRefl.GetNumberCols())].index('refine')
        for r in range(G2frame.PawleyRefl.GetNumberRows()):
            G2frame.PawleyRefl.GetTable().SetValue(
                r,refcol,
                not G2frame.PawleyRefl.GetTable().GetValueAsBool(r,refcol))
        G2frame.PawleyRefl.ForceRefresh()

#### Fourier routines ################################################################################
    def FillMapPeaksGrid():

        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                if G2frame.MapPeaks.IsSelection():
                    G2frame.MapPeaks.ClearSelection()
                else:
                    for row in range(G2frame.MapPeaks.GetNumberRows()):
                        G2frame.MapPeaks.SelectRow(row,True)

            elif c < 0:                   #only row clicks
                if event.ControlDown():
                    if r in getAtomSelections(G2frame.MapPeaks):
                        G2frame.MapPeaks.DeselectRow(r)
                    else:
                        G2frame.MapPeaks.SelectRow(r,True)
                elif event.ShiftDown():
                    indxs = getAtomSelections(G2frame.MapPeaks)
                    G2frame.MapPeaks.ClearSelection()
                    ibeg = 0
                    if indxs:
                        ibeg = indxs[-1]
                    for row in range(ibeg,r+1):
                        G2frame.MapPeaks.SelectRow(row,True)
                else:
                    G2frame.MapPeaks.ClearSelection()
                    G2frame.MapPeaks.SelectRow(r,True)
            elif r < 0:                 #a column pick
                mapPeaks = data['Map Peaks']
                c =  event.GetCol()
                if colLabels[c] == 'mag':   #big to small order
                    mapPeaks = G2mth.sortArray(mapPeaks,c,reverse=True)
                elif colLabels[c] in ['x','y','z','dzero','dcent']:     #small to big
                    mapPeaks = G2mth.sortArray(mapPeaks,c)
                else:
                    return
                data['Map Peaks'] = mapPeaks
                wx.CallAfter(FillMapPeaksGrid)
            G2plt.PlotStructure(G2frame,data)

        # beginning of FillMapPeaksGrid()
        G2frame.GetStatusBar().SetStatusText('',1)
        oldSizer = MapPeakList.GetSizer()
        if oldSizer: # 2nd+ use, clear out old entries
            for i in oldSizer.GetChildren(): # look for grids in sizer
                if type(i.GetWindow()) is G2G.GSGrid:
                    oldSizer.Detach(i.GetWindow())  # don't delete them
            oldSizer.Clear(True)

        topSizer = G2frame.dataWindow.topBox
        topSizer.Clear(True)
        parent = G2frame.dataWindow.topPanel
        lbl= f"Fourier map peak positions for {data['General']['Name']!r}"[:60]
        topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
        topSizer.Add((-1,-1),1,wx.EXPAND)
        topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
        wx.CallAfter(G2frame.dataWindow.SetDataSize)

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if 'Map Peaks' in data:
            G2frame.GetStatusBar().SetStatusText('Double click any column heading to sort',1)
            mapPeaks = data['Map Peaks']
            rowLabels = []
            for i in range(len(mapPeaks)): rowLabels.append(str(i))
            colLabels = ['mag','x','y','z','dzero','dcent']
            Types = 6*[wg.GRID_VALUE_FLOAT+':10,4',]
            G2frame.MapPeaksTable = G2G.Table(mapPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            G2frame.MapPeaks.Show(True)
            G2frame.MapPeaks.SetTable(G2frame.MapPeaksTable, True)
            G2frame.MapPeaks.Unbind(wg.EVT_GRID_LABEL_LEFT_CLICK)
            G2frame.MapPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
            for r in range(G2frame.MapPeaks.GetNumberRows()):
                for c in range(G2frame.MapPeaks.GetNumberCols()):
                    G2frame.MapPeaks.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.MapPeaks.SetMargins(0,0)
            G2frame.MapPeaks.AutoSizeColumns(False)
            G2frame.MapPeaks.SetRowLabelSize(40)
            mainSizer.Add(G2frame.MapPeaks,1,wx.EXPAND)
            G2frame.MapPeaks.SetScrollRate(10,10) # allow grid to scroll
        else:
            mainSizer.Add(wx.StaticText(MapPeakList,label=' Map peak list is empty'),0,wx.ALL,10)
            G2frame.MapPeaks.Show(False)
        SetPhaseWindow(MapPeakList,mainSizer)

    def OnPeaksMove(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            peakMax = np.amax(mapPeaks.T[0])
            Ind = getAtomSelections(G2frame.MapPeaks)
            pgbar = wx.ProgressDialog('Move peaks','Map peak no. 0 processed',len(Ind)+1,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
            for i,ind in enumerate(Ind):
                mag,x,y,z = mapPeaks[ind][:4]
                AtomAdd(x,y,z,'H',Name='M '+'%d'%(int(100*mag/peakMax)))
                pgbar.Update(i+1,'Map peak no. %d processed'%ind)
            pgbar.Destroy()
            G2plt.PlotStructure(G2frame,data)

    def OnPeaksClear(event):
        data['Map Peaks'] = []
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)

    def OnPeaksSave(event):
        if 'Map Peaks' in data:
            mapPeaks = data['Map Peaks']
            pfName = PhaseName+'_peaks.csv'
            pfFile = ''
            pth = G2G.GetExportPath(G2frame)
            dlg = wx.FileDialog(G2frame, 'Choose map peaks file name', pth, pfName,
                'csv (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    pfFile = dlg.GetPath()
            finally:
                dlg.Destroy()

            if pfFile:
                pf = open(pfFile,'w')
                pf.write('"%s"\n'%(PhaseName))
                pf.write(' mag, x, y, z, dzero, dcent \n')
                for peak in mapPeaks:
                    pf.write(' %.4f, %.4f, %.4f, %.4f, %.4f, %.4f \n'%(peak[0],peak[1],peak[2],peak[3],peak[4],peak[5]))
                pf.close()

    def OnPeaksDelete(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            Ind = getAtomSelections(mapPeaks)
            Ind.sort()
            Ind.reverse()
            for ind in Ind:
                mapPeaks = np.delete(mapPeaks,ind,0)
            data['Map Peaks'] = mapPeaks
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)

    def OnPeaksInvert(event):
        if 'Map Peaks' in data:
            generalData = data['General']
            mapData = generalData['Map']
            try:
                mapData['rho'] = np.flip(mapData['rho'],(0,1,2))
            except TypeError:
                mapData['rho'] = np.flip(mapData['rho'],0)
                mapData['rho'] = np.flip(mapData['rho'],1)
                mapData['rho'] = np.flip(mapData['rho'],2)
            mapData['rho'] = np.roll(np.roll(np.roll(mapData['rho'],1,axis=0),1,axis=1),1,axis=2)
            OnSearchMaps(event)
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)

    def OnRollMap(event):
        if 'Map Peaks' in data:
            half = np.array([.5,.5,.5])
            mapPeaks = data['Map Peaks']
            generalData = data['General']
            mapData = generalData['Map']
            if mapData['Flip'] != True:
                wx.MessageBox('Only valid for charge flip maps',caption='Roll map error',style=wx.ICON_EXCLAMATION)
                return
            Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
            dims = mapData['rho'].shape
            dlg = G2G.MultiDataDialog(G2frame,title='Roll map shifts',
                prompts=['delt-X (-1. to 1.)','delt-Y (-1. to 1.)',
                         'delt-Z (-1. to 1.)'],values=[0.,0.,0.,],
                    limits=[[-1.,1.],[-1.,1.],[-1.,1.]],formats=['%.6f','%.6f','%.6f'])

            if dlg.ShowModal() == wx.ID_OK:
                rollsteps = dlg.GetValues()
                dxy = np.array([float(R) for R in rollsteps])
                rollsteps = np.array([round(float(R)*dims[iR]) for iR,R in enumerate(rollsteps)])
                mapData['rho'] = np.roll(np.roll(np.roll(mapData['rho'],rollsteps[0],axis=0),rollsteps[1],axis=1),rollsteps[2],axis=2)
                for peak in mapPeaks:
                    peak[1:4] += dxy
                    peak[1:4] %= 1.
                    peak[4] = np.sqrt(np.sum(np.inner(Amat,peak[1:4])**2))
                    if len(peak)>4:
                        peak[5] = np.sqrt(np.sum(np.inner(Amat,(peak[1:4]-half)%1.)**2))
                FillMapPeaksGrid()
                G2plt.PlotStructure(G2frame,data)
            dlg.Destroy()

    def OnPeaksEquiv(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            Ind = getAtomSelections(mapPeaks)
            if Ind:
                wx.BeginBusyCursor()
                try:
                    Ind = G2mth.PeaksEquiv(data,Ind)
                    for r in range(G2frame.MapPeaks.GetNumberRows()):
                        if r in Ind:
                            G2frame.MapPeaks.SelectRow(r,addToSelected=True)
                        else:
                            G2frame.MapPeaks.DeselectRow(r)
                finally:
                    wx.EndBusyCursor()
                G2plt.PlotStructure(G2frame,data)

    def OnShowBonds(event):
        generalData = data['General']
        if generalData['Map'].get('Show bonds',False):
            generalData['Map']['Show bonds'] = False
            G2frame.dataWindow.MapPeaksEdit.SetLabel(G2G.wxID_SHOWBONDS,'Show bonds')
        else:
            generalData['Map']['Show bonds'] = True
            G2frame.dataWindow.MapPeaksEdit.SetLabel(G2G.wxID_SHOWBONDS,'Hide bonds')
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)

    def OnPeaksUnique(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            Ind = getAtomSelections(G2frame.MapPeaks)
            if Ind:
                choice = ['x=0','y=0','z=0','origin','center']
                dlg = wx.SingleChoiceDialog(G2frame,'Peaks closest to:','Select',choice)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()+1
                    dlg.Destroy()
                else:
                    dlg.Destroy()
                    return
                pgbar = wx.ProgressDialog('Unique peaks','Map peak no. 0 processed',len(Ind)+1,
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)

                Ind = G2mth.PeaksUnique(data,Ind,sel,pgbar)
                print (' No. unique peaks: %d Unique peak fraction: %.3f'%(len(Ind),float(len(Ind))/len(mapPeaks)))
                tbl = G2frame.MapPeaks.GetTable().data
                tbl[:] = [t for i,t in enumerate(tbl) if i in Ind] + [
                    t for i,t in enumerate(tbl) if i not in Ind]
                for r in range(G2frame.MapPeaks.GetNumberRows()):
                    if r < len(Ind):
                        G2frame.MapPeaks.SelectRow(r,addToSelected=True)
                    else:
                        G2frame.MapPeaks.DeselectRow(r)
                G2frame.MapPeaks.ForceRefresh()
                G2plt.PlotStructure(G2frame,data)

    def OnPeaksViewPoint(event):
        # set view point
        indx = getAtomSelections(G2frame.MapPeaks)
        if not indx:
            G2frame.ErrorDialog('Set viewpoint','No peaks selected')
            return
        mapPeaks = data['Map Peaks']
        drawingData = data['Drawing']
        drawingData['viewPoint'][0] = mapPeaks[indx[0]][1:4]
        G2plt.PlotStructure(G2frame,data)

    def OnPeaksDistVP(event):
        # distance to view point
        indx = getAtomSelections(G2frame.MapPeaks)
        if not indx:
            G2frame.ErrorDialog('Peak distance','No peaks selected')
            return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        mapPeaks = data['Map Peaks']
        drawingData = data['Drawing']
        viewPt = np.array(drawingData['viewPoint'][0])
        print (' Distance from view point at %.3f %.3f %.3f to:'%(viewPt[0],viewPt[1],viewPt[2]))
        colLabels = [G2frame.MapPeaks.GetColLabelValue(c) for c in range(G2frame.MapPeaks.GetNumberCols())]
        cx = colLabels.index('x')
        cm = colLabels.index('mag')
        for i in indx:
            peak = mapPeaks[i]
            Dx = np.array(peak[cx:cx+3])-viewPt
            dist = np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0))
            print ('Peak: %5d mag= %8.2f distance = %.3f'%(i,peak[cm],dist))

    def OnPeaksDA(event):
        #distance, angle
        indx = getAtomSelections(G2frame.MapPeaks)
        if len(indx) not in [2,3]:
            G2frame.ErrorDialog('Peak distance/angle','Wrong number of atoms for distance or angle calculation')
            return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        mapPeaks = data['Map Peaks']
        xyz = []
        for i in indx:
            xyz.append(mapPeaks[i][1:4])
        if len(indx) == 2:
            print (' distance for atoms %s = %.3f'%(str(indx),G2mth.getRestDist(xyz,Amat)))
        else:
            print (' angle for atoms %s = %.2f'%(str(indx),G2mth.getRestAngle(xyz,Amat)))

    def OnFourierMaps(event):
        generalData = data['General']
        mapData = generalData['Map']
        reflNames = mapData['RefList']
        if not generalData['Map']['MapType']:
            G2frame.ErrorDialog('Fourier map error','Fourier map type not defined')
            return
        if not reflNames[0]:
            G2frame.ErrorDialog('Fourier map error','No reflections defined for Fourier map')
            return
        phaseName = generalData['Name']
        ReflData = GetReflData(G2frame,phaseName,reflNames)
        if ReflData == None:
            G2frame.ErrorDialog('Fourier map error','No reflections defined for Fourier map')
            return
        if 'Omit' in mapData['MapType']:
            dim = '3D '
            pgbar = wx.ProgressDialog('Omit map','Blocks done',65,
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
            mapData.update(G2mth.OmitMap(data,ReflData,pgbar))
            pgbar.Destroy()
        else:
            if generalData['Modulated']:
                dim = '4D '
                G2mth.Fourier4DMap(data,ReflData)
            else:
                dim = '3D '
                G2mth.FourierMap(data,ReflData)
        mapData['Flip'] = False
        mapSig = np.std(mapData['rho'])
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        data['Drawing']['contourLevel'] = 1.
        data['Drawing']['mapSize'] = 10.
        data['Drawing']['showMap'] = True
        ftext = dim+mapData['MapType']+' computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig)
        print (ftext)
        G2frame.AddToNotebook('Fourier '+ftext,'FM')
        UpdateDrawAtoms(G2frame,data)
        G2plt.PlotStructure(G2frame,data)

    def OnFourClear(event):
        generalData = data['General']
        generalData['Map'] = mapDefault.copy()
        data['Drawing']['showMap'] = False
        G2plt.PlotStructure(G2frame,data)

# map printing for testing purposes
    def printRho(SGLaue,rho,rhoMax):
        dim = len(rho.shape)
        if dim == 2:
            ix,jy = rho.shape
            for j in range(jy):
                line = ''
                if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                    line += (jy-j)*'  '
                for i in range(ix):
                    r = int(100*rho[i,j]/rhoMax)
                    line += '%4d'%(r)
                print (line+'\n')
        else:
            ix,jy,kz = rho.shape
            for k in range(kz):
                print ('k = %d'%k)
                for j in range(jy):
                    line = ''
                    if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                        line += (jy-j)*'  '
                    for i in range(ix):
                        r = int(100*rho[i,j,k]/rhoMax)
                        line += '%4d'%(r)
                    print (line+'\n')
## keep this

    def OnSearchMaps(event):

        print (' Begin fourier map search - can take some time')
        time0 = time.time()
        generalData = data['General']
        drawingData = data['Drawing']
        mapData = generalData['Map']
        if len(mapData['rho']):
            wx.BeginBusyCursor()
            try:
                peaks,mags,dzeros,dcents = G2mth.SearchMap(generalData,drawingData)
                if mapData['MapType'] in ['delt-F',] or 'N' in mapData['Type']:
                    npeaks,nmags,ndzeros,ndcents = G2mth.SearchMap(generalData,drawingData,Neg=True)
                    peaks = np.concatenate((peaks,npeaks))
                    mags = np.concatenate((mags,nmags))
                    dzeros = np.concatenate((dzeros,ndzeros))
                    dcents = np.concatenate((dcents,ndcents))
            finally:
                wx.EndBusyCursor()
            if len(peaks):
                mapPeaks = np.concatenate((mags,peaks,dzeros,dcents),axis=1)
                data['Map Peaks'] = G2mth.sortArray(mapPeaks,0,reverse=True)
            print (' Map search finished, time = %.2fs'%(time.time()-time0))
            print (' No.peaks found: %d'%len(peaks))
            Page = G2frame.phaseDisplay.FindPage('Map peaks')
            G2frame.phaseDisplay.SetSelection(Page)
            wx.CallAfter(FillMapPeaksGrid)
            UpdateDrawAtoms(G2frame,data)
        else:
            print ('No map available')

    def On4DChargeFlip(event):
        generalData = data['General']
        mapData = generalData['Map']
        map4DData = generalData['4DmapData']
        flipData = generalData['Flip']
        reflNames = flipData['RefList']
        if not reflNames[0]:
            G2frame.ErrorDialog('Charge flip error','No reflections defined for charge flipping')
            return
        phaseName = generalData['Name']
        ReflData = GetReflData(G2frame,phaseName,reflNames)
        if ReflData == None:
            G2frame.ErrorDialog('Charge flip error','No reflections defined for charge flipping')
            return
        pgbar = wx.ProgressDialog('Charge flipping','Residual Rcf =',101,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            result = G2mth.SSChargeFlip(data,ReflData,pgbar)
        finally:
            pgbar.Destroy()
        G2frame.AddToNotebook(f'4D Charge flip: {result[2]}\n{result[3]}','CF')
        mapData.update(result[0])
        map4DData.update(result[1])
        mapData['Flip'] = True
        mapSig = np.std(mapData['rho'])
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        data['Drawing']['contourLevel'] = 1.
        data['Drawing']['mapSize'] = 10.
        print (' 4D Charge flip map computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig))
        if mapData['Rcf'] < 99.:
            OnSearchMaps(event)             #does a plot structure at end
        else:
            print ('Bad charge flip map - no peak search done')

    def OnChargeFlip(event):
        generalData = data['General']
        mapData = generalData['Map']
        flipData = generalData['Flip']
        reflNames = flipData['RefList']
        if not reflNames[0]:
            G2frame.ErrorDialog('Charge flip error','No reflections defined for charge flipping')
            return
        phaseName = generalData['Name']
        ReflData = GetReflData(G2frame,phaseName,reflNames)
        if ReflData == None:
            G2frame.ErrorDialog('Charge flip error','No reflections defined for charge flipping')
            return
        pgbar = wx.ProgressDialog('Charge flipping','Residual Rcf =',101,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        testNames = ['%3d%3d%3d'%(h,k,l) for h,k,l in flipData['testHKL']]
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            result = G2mth.ChargeFlip(data,ReflData,pgbar)
            G2frame.AddToNotebook(f'Charge flip: {result[2]}\n{result[3]}','CF')
            mapData.update(result[0])
            X = range(len(result[1]))
            Y = 180.*np.array(result[1]).T/np.pi
            XY = [[X,y] for y in Y]
            XY = np.array(XY).reshape((5,2,-1))
            G2plt.PlotXY(G2frame,XY,labelX='charge flip cycle',labelY='phase, deg',newPlot=True,
                Title='Test HKL phases',lines=True,names=testNames)
        finally:
            pgbar.Destroy()
        mapData['Flip'] = True
        mapSig = np.std(mapData['rho'])
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        data['Drawing']['contourLevel'] = 1.
        data['Drawing']['mapSize'] = 10.
        data['Drawing']['showMap'] = True
        print (' Charge flip map computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig))
        if mapData['Rcf'] < 99.:
            OnSearchMaps(event)             #does a plot structure at end
        else:
            print ('Bad charge flip map - no peak search done')

    def OnTextureRefine(event):
        General = data['General']
        phaseName = General['Name']
        keyList = G2frame.GetHistogramNames('PWDR')
        histNames = []
        refData = {}
        Gangls = {}
        for name in keyList:
            if 'PWDR' in name:
                im = 0
                it = 0
                histNames.append(name)
                Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))
                Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                Gangls[name] = copy.copy([Sample[item] for item in['Phi','Chi','Omega','Azimuth']])
                RefDict = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Reflection Lists'))[phaseName]
                Refs = RefDict['RefList'].T  #np.array!
                if RefDict['Super']: im = 1     #(3+1) offset for m
                if 'T' in RefDict['Type']:
                    it = 3  #TOF offset for alp, bet, wave
                    tth = np.ones_like(Refs[0])*Inst[0]['2-theta'][0]
                    refData[name] = np.column_stack((Refs[0],Refs[1],Refs[2],tth,Refs[8+im],Refs[12+im+it],np.zeros_like(Refs[0])))
                else:   # xray - typical caked 2D image data
                    refData[name] = np.column_stack((Refs[0],Refs[1],Refs[2],Refs[5+im],Refs[8+im],Refs[12+im+it],np.zeros_like(Refs[0])))
        pgbar = wx.ProgressDialog('Texture fit','Residual = %5.2f'%(101.0),101,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        Error = G2mth.FitTexture(General,Gangls,refData,keyList,pgbar)
        pgbar.Destroy()
        if Error:
            wx.MessageBox(Error,caption='Fit Texture Error',style=wx.ICON_EXCLAMATION)
        XY = []
        for hist in keyList:
            x = refData[hist].T[5].T
            y = refData[hist].T[6].T
            xy = [x,y]
            XY.append(np.array(xy))
        G2plt.PlotXY(G2frame,XY,labelX='POobs',labelY='POcalc',lines=False,newPlot=False,Title='Texture fit error')
        G2phsG2.UpdateTexture(G2frame,data)
        G2plt.PlotTexture(G2frame,data,Start=False)

    def OnTextureClear(event):
        print ('clear texture? - does nothing')

#### Phase page routines ###############################################################################
    def FillSelectPageMenu(TabSelectionIdDict, menuBar):
        '''Fill "Select tab" menu with menu items for each tab and assign
        bindings to the menu item to switch between phase tabs
        '''
        def OnSelectPage(event):
            'Called when an item is selected from the Select page menu'
            # lookup the menu item that called us and get its text
            tabname = TabSelectionIdDict.get(event.GetId())
            if not tabname:
                print ('Warning: menu item not in dict! id=%d'%event.GetId())
                return
            # find the matching tab
            for PageNum in range(G2frame.phaseDisplay.GetPageCount()):
                if tabname == G2frame.phaseDisplay.GetPageText(PageNum):
                    G2frame.phaseDisplay.SetSelection(PageNum)
                    return
            else:
                print ("Warning: tab "+tabname+" was not found")
        mid = menuBar.FindMenu('Select tab')
        menu = menuBar.GetMenu(mid)
        for ipage,page in enumerate(Pages):
            if menu.FindItem(page) == wx.NOT_FOUND: # is tab already in menu?
                Id = wx.NewId()
                TabSelectionIdDict[Id] = page
                menu.Append(Id,page,'')
                G2frame.Bind(wx.EVT_MENU, OnSelectPage, id=Id)

    def OnPageChanged(event):
        '''This is called every time that a Notebook tab button is pressed
        on a Phase data item window
        '''
        page = event.GetSelection()
        G2frame.phaseDisplay.SetSize(G2frame.dataWindow.GetClientSize())    #TODO -almost right
        # force edits in open grids to complete
        newlist = []
        for p in G2frame.phaseDisplay.gridList:
            if not p: continue   # skip deleted grids
            try:
                p.ClearGrid()
                newlist.append(p)
            except:
                pass
        G2frame.phaseDisplay.gridList = newlist  # remove deleted grids from lis
        # let the pending events clear before posting the new NB page
        wx.CallAfter(ChangePage,page)

    def ChangePage(page):
        text = G2frame.phaseDisplay.GetPageText(page)
        G2frame.lastSelectedPhaseTab = text
        G2frame.dataWindow.helpKey = 'Phase-'+text # use name of Phase tab for help lookup
        if text == 'General':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataGeneral)
            UpdateGeneral()
        elif text == 'Data': # only when conf 'SeparateHistPhaseTreeItem' is False
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataMenu)
            G2plt.PlotSizeStrainPO(G2frame,data,hist='')
            G2ddG.UpdateDData(G2frame,DData,data)
        elif text == 'Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.AtomsMenu)
            FillAtomsGrid(Atoms)
        elif text == 'Layers':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.LayerData)
            G2phsG2.UpdateLayerData(G2frame,data)
        elif text == 'Wave Data' and data['General']['Modulated']:
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.WavesData)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            G2phsG2.UpdateWavesData(G2frame,data)
        elif text == 'Dysnomia':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.MEMMenu)
            G2phsG2.UpdateDysnomia(G2frame,data)
        elif text == 'RMC':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.FRMCMenu)
            G2rmcG.UpdateRMC(G2frame,data)
        elif text == 'ISODISTORT':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.ISODData)
            G2phsG2.UpdateISODISTORT(G2frame,data)
        elif text == 'Draw Options':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataDrawOptions)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            UpdateDrawOptions(G2frame,data)
        elif text == 'Draw Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DrawAtomsMenu)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            UpdateDrawAtoms(G2frame,data)
        elif text == 'Deformation':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DeformationMenu)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            G2phsG2.UpdateDeformation(G2frame,data,None)
        elif text == 'RB Models':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodiesMenu)
            FillRigidBodyGrid()
        elif text == 'Map peaks':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.MapPeaksMenu)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            FillMapPeaksGrid()
        elif text == 'MC/SA':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.MCSAMenu)
            G2plt.PlotStructure(G2frame,data,firstCall=True)
            UpdateMCSA()
        elif text == 'Texture':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.TextureMenu)
            G2plt.PlotTexture(G2frame,data,Start=True)
            G2phsG2.UpdateTexture(G2frame,data)
        elif text == 'Pawley reflections':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PawleyMenu)
            FillPawleyReflectionsGrid()
        else:
            G2gd.SetDataMenuBar(G2frame)

    def FillMenus():
        '''Create the Select tab menus and bind to all menu items
        '''
        # General
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DataGeneral)
        G2frame.Bind(wx.EVT_MENU, OnFourierMaps, id=G2G.wxID_FOURCALC)
        G2frame.Bind(wx.EVT_MENU, OnSearchMaps, id=G2G.wxID_FOURSEARCH)
        G2frame.Bind(wx.EVT_MENU, OnChargeFlip, id=G2G.wxID_CHARGEFLIP)
        G2frame.Bind(wx.EVT_MENU, On4DChargeFlip, id=G2G.wxID_4DCHARGEFLIP)
        G2frame.Bind(wx.EVT_MENU, OnFourClear, id=G2G.wxID_FOURCLEAR)
        G2frame.Bind(wx.EVT_MENU, OnRunSingleMCSA, id=G2G.wxID_SINGLEMCSA)
        G2frame.Bind(wx.EVT_MENU, OnRunMultiMCSA, id=G2G.wxID_MULTIMCSA)
        G2frame.Bind(wx.EVT_MENU, OnTransform, id=G2G.wxID_TRANSFORMSTRUCTURE)
        G2frame.Bind(wx.EVT_MENU, OnTransform2Std, id=G2G.wxID_TRANSFORMSTD)
        G2frame.Bind(wx.EVT_MENU, OnSuperSearch, id=G2G.wxID_SUPERSRCH)
        G2frame.Bind(wx.EVT_MENU, OnISOSearch, id=G2G.wxID_ISOSRCH)
        G2frame.Bind(wx.EVT_MENU, OnSubSearch, id=G2G.wxID_SUBSRCH)
        G2frame.Bind(wx.EVT_MENU, OnCompare, id=G2G.wxID_COMPARESTRUCTURE)
        G2frame.Bind(wx.EVT_MENU, OnCompareCells, id=G2G.wxID_COMPARECELLS)
        G2frame.Bind(wx.EVT_MENU, OnUseBilbao, id=G2G.wxID_USEBILBAOMAG)
        G2frame.Bind(wx.EVT_MENU, OnApplySubgroups, id=G2G.wxID_USEBILBAOSUB)
        G2frame.Bind(wx.EVT_MENU, OnValidProtein, id=G2G.wxID_VALIDPROTEIN)
        for Id in G2frame.dataWindow.ReplaceMenuId:     #loop over submenu items
            G2frame.Bind(wx.EVT_MENU, OnReplacePhase, id=Id)
        # Data (unless Hist/Phase tree entry shown)
        if not GSASIIpath.GetConfigValue('SeparateHistPhaseTreeItem',False):
            FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DataMenu)
            G2frame.Bind(wx.EVT_MENU, OnDataUse, id=G2G.wxID_DATAUSE)
            G2frame.Bind(wx.EVT_MENU, OnDataCopy, id=G2G.wxID_DATACOPY)
            G2frame.Bind(wx.EVT_MENU, OnDataCopyFlags, id=G2G.wxID_DATACOPYFLAGS)
            G2frame.Bind(wx.EVT_MENU, OnSelDataCopy, id=G2G.wxID_DATASELCOPY)
            G2frame.Bind(wx.EVT_MENU, OnSelDataRead, id=G2G.wxID_DATASELREAD)
            G2frame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2G.wxID_PWDRADD)
            G2frame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2G.wxID_HKLFADD)
            G2frame.Bind(wx.EVT_MENU, OnDataDelete, id=G2G.wxID_DATADELETE)
            G2frame.Bind(wx.EVT_MENU, OnDataApplyStrain, id=G2G.wxID_DATADIJ)
        # Atoms
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.AtomsMenu)
        G2frame.Bind(wx.EVT_MENU, OnSetAll, id=G2G.wxID_ATOMSSETALL)
        G2frame.Bind(wx.EVT_MENU, OnSetbyList, id=G2G.wxID_ATOMSSETLST)
        G2frame.Bind(wx.EVT_MENU, AtomRefine, id=G2G.wxID_ATOMSSETSEL)
        G2frame.Bind(wx.EVT_MENU, AtomModify, id=G2G.wxID_ATOMSMODIFY)
        G2frame.Bind(wx.EVT_MENU, OnAtomInsert, id=G2G.wxID_ATOMSEDITINSERT)
        G2frame.Bind(wx.EVT_MENU, OnHydAtomAdd, id=G2G.wxID_ADDHATOM)
        G2frame.Bind(wx.EVT_MENU, AtomDelete, id=G2G.wxID_ATOMSEDITDELETE)
        G2frame.Bind(wx.EVT_MENU, AtomTransform, id=G2G.wxID_ATOMSTRANSFORM)
#        G2frame.Bind(wx.EVT_MENU, AtomRotate, id=G2G.wxID_ATOMSROTATE)
        G2frame.Bind(wx.EVT_MENU, SetAtomsViewPoint, id=G2G.wxID_ATOMSSETVP)

        G2frame.Bind(wx.EVT_MENU, OnAtomAdd, id=G2G.wxID_ATOMSEDITADD)
        G2frame.Bind(wx.EVT_MENU, OnAtomViewAdd, id=G2G.wxID_ATOMSVIEWADD)
        G2frame.Bind(wx.EVT_MENU, OnAtomViewInsert, id=G2G.wxID_ATOMVIEWINSERT)
        G2frame.Bind(wx.EVT_MENU, OnHydAtomUpdate, id=G2G.wxID_UPDATEHATOM)
        G2frame.Bind(wx.EVT_MENU, OnAtomMove, id=G2G.wxID_ATOMMOVE)
        G2frame.Bind(wx.EVT_MENU, MakeMolecule, id=G2G.wxID_MAKEMOLECULE)
        G2frame.Bind(wx.EVT_MENU, CollectAtoms, id=G2G.wxID_COLLECTATOMS)
        G2frame.Bind(wx.EVT_MENU, OnReloadDrawAtoms, id=G2G.wxID_RELOADDRAWATOMS)
        G2frame.Bind(wx.EVT_MENU, OnDistAngle, id=G2G.wxID_ATOMSDISAGL)
        G2frame.Bind(wx.EVT_MENU, OnDistAnglePrt, id=G2G.wxID_ATOMSPDISAGL)
        G2frame.Bind(wx.EVT_MENU, OnDistAngleHist, id=G2G.wxID_ATOMSBNDANGLHIST)
        G2frame.Bind(wx.EVT_MENU, OnSaveDAHist, id=G2G.wxID_ATOMSSAVEHIST)
        G2frame.Bind(wx.EVT_MENU, OnClearDAHist, id=G2G.wxID_ATOMSCLEARHIST)
        G2frame.Bind(wx.EVT_MENU, OnFracSplit, id=G2G.wxID_ATOMFRACSPLIT)
        G2frame.Bind(wx.EVT_MENU, OnDensity, id=G2G.wxID_ATOMSDENSITY)
        G2frame.Bind(wx.EVT_MENU, OnShowIsoDistortCalc, id=G2G.wxID_ISODISP)
        if 'HydIds' in data['General']:
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,True)
        else:
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,False)
        for Id in G2frame.dataWindow.ReImportMenuId:     #loop over submenu items
            G2frame.Bind(wx.EVT_MENU, OnReImport, id=Id)
        # Wave Data
        if data['General']['Modulated']:
            FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.WavesData)
            G2frame.Bind(wx.EVT_MENU, OnWaveVary, id=G2G.wxID_WAVEVARY)
        # Dysnomia (MEM)
        if data['General']['doDysnomia']:
            FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.MEMMenu)
            G2frame.Bind(wx.EVT_MENU,OnLoadDysnomia,id=G2G.wxID_LOADDYSNOMIA)
            G2frame.Bind(wx.EVT_MENU,OnSaveDysnomia,id=G2G.wxID_SAVEDYSNOMIA)
            G2frame.Bind(wx.EVT_MENU,OnRunDysnomia,id=G2G.wxID_RUNDYSNOMIA)
        # Stacking faults
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.LayerData)
        G2frame.Bind(wx.EVT_MENU, OnCopyPhase, id=G2G.wxID_COPYPHASE)
        G2frame.Bind(wx.EVT_MENU, OnLoadDIFFaX, id=G2G.wxID_LOADDIFFAX)
        G2frame.Bind(wx.EVT_MENU, OnSimulate, id=G2G.wxID_LAYERSIMULATE)
        G2frame.Bind(wx.EVT_MENU, OnFitLayers, id=G2G.wxID_LAYERSFIT)
        G2frame.Bind(wx.EVT_MENU, OnSeqSimulate, id=G2G.wxID_SEQUENCESIMULATE)
        # Draw Options
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DataDrawOptions)
        # Draw Atoms
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DrawAtomsMenu)
        G2frame.Bind(wx.EVT_MENU, DrawAtomStyle, id=G2G.wxID_DRAWATOMSTYLE)
        G2frame.Bind(wx.EVT_MENU, DrawAtomLabel, id=G2G.wxID_DRAWATOMLABEL)
        G2frame.Bind(wx.EVT_MENU, DrawAtomColor, id=G2G.wxID_DRAWATOMCOLOR)
        G2frame.Bind(wx.EVT_MENU, ResetAtomColors, id=G2G.wxID_DRAWATOMRESETCOLOR)
        G2frame.Bind(wx.EVT_MENU, OnEditAtomRadii, id=G2G.wxID_DRWAEDITRADII)
        G2frame.Bind(wx.EVT_MENU, SetViewPoint, id=G2G.wxID_DRAWVIEWPOINT)
        G2frame.Bind(wx.EVT_MENU, AddSymEquiv, id=G2G.wxID_DRAWADDEQUIV)
        G2frame.Bind(wx.EVT_MENU, AddSphere, id=G2G.wxID_DRAWADDSPHERE)
        G2frame.Bind(wx.EVT_MENU, AddBox, id=G2G.wxID_DRAWADDBOX)
        G2frame.Bind(wx.EVT_MENU, TransformSymEquiv, id=G2G.wxID_DRAWTRANSFORM)
        G2frame.Bind(wx.EVT_MENU, FillCoordSphere, id=G2G.wxID_DRAWFILLCOORD)
        G2frame.Bind(wx.EVT_MENU, FillUnitCell, id=G2G.wxID_DRAWFILLCELL)
        G2frame.Bind(wx.EVT_MENU, DrawAtomsDelete, id=G2G.wxID_DRAWDELETE)
        G2frame.Bind(wx.EVT_MENU, OnReloadDrawAtoms, id=G2G.wxID_RELOADATOMS)
        G2frame.Bind(wx.EVT_MENU, OnDrawDistVP, id=G2G.wxID_DRAWDISTVP)
        G2frame.Bind(wx.EVT_MENU, OnDrawDAT, id=G2G.wxID_DRAWDISAGLTOR)
        G2frame.Bind(wx.EVT_MENU, OnDrawPlane, id=G2G.wxID_DRAWPLANE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRBOND)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRANGLE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRPLANE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRCHIRAL)
        G2frame.Bind(wx.EVT_MENU, OnDefineRB, id=G2G.wxID_DRAWDEFINERB)
        G2frame.Bind(wx.EVT_MENU, FillMolecule, id=G2G.wxID_DRAWADDMOLECULE)
        G2frame.Bind(wx.EVT_MENU, MapVoid, id=G2G.wxID_DRAWVOIDMAP)
        G2frame.Bind(wx.EVT_MENU, SelDrawList, id=G2G.wxID_DRAWSETSEL)
        G2frame.Bind(wx.EVT_MENU, DrawLoadSel, id=G2G.wxID_DRAWLOADSEL)
        G2frame.Bind(wx.EVT_MENU, RandomizedAction, id=G2G.wxID_DRAWRANDOM)

        # Deformation form factors
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DeformationMenu)
        G2frame.Bind(wx.EVT_MENU, SelDeformAtom, id=G2G.wxID_DEFORMSETSEL)
        G2frame.Bind(wx.EVT_MENU, SetDefDist, id=G2G.wxID_DEFORMDISTSET)

        # RB Models
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.RigidBodiesMenu)
        G2frame.Bind(wx.EVT_MENU, OnAutoFindResRB, id=G2G.wxID_AUTOFINDRESRB)
        G2frame.Bind(wx.EVT_MENU, OnRBAssign, id=G2G.wxID_ASSIGNATMS2RB)
        G2frame.Bind(wx.EVT_MENU, OnRBCopyParms, id=G2G.wxID_COPYRBPARMS)
        G2frame.Bind(wx.EVT_MENU, OnGlobalResRBTherm, id=G2G.wxID_GLOBALTHERM)
        G2frame.Bind(wx.EVT_MENU, OnGlobalResRBRef, id=G2G.wxID_GLOBALRESREFINE)
        G2frame.Bind(wx.EVT_MENU, OnRBRemoveAll, id=G2G.wxID_RBREMOVEALL)
        # Map peaks
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.MapPeaksMenu)
        G2frame.Bind(wx.EVT_MENU, OnPeaksMove, id=G2G.wxID_PEAKSMOVE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksViewPoint, id=G2G.wxID_PEAKSVIEWPT)
        G2frame.Bind(wx.EVT_MENU, OnPeaksDistVP, id=G2G.wxID_PEAKSDISTVP)
        G2frame.Bind(wx.EVT_MENU, OnPeaksDA, id=G2G.wxID_PEAKSDA)
        G2frame.Bind(wx.EVT_MENU, OnShowBonds, id=G2G.wxID_SHOWBONDS)
        G2frame.Bind(wx.EVT_MENU, OnPeaksEquiv, id=G2G.wxID_FINDEQVPEAKS)
        G2frame.Bind(wx.EVT_MENU, OnPeaksInvert, id=G2G.wxID_INVERTPEAKS)
        G2frame.Bind(wx.EVT_MENU, OnRollMap, id=G2G.wxID_ROLLMAP)
        G2frame.Bind(wx.EVT_MENU, OnPeaksUnique, id=G2G.wxID_PEAKSUNIQUE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksSave, id=G2G.wxID_PEAKSSAVE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksDelete, id=G2G.wxID_PEAKSDELETE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksClear, id=G2G.wxID_PEAKSCLEAR)
        # RMCProfile/fullrmc/PDFfit
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.FRMCMenu)
        G2frame.Bind(wx.EVT_MENU, OnSetupRMC, id=G2G.wxID_SETUPRMC)
        G2frame.Bind(wx.EVT_MENU, OnRunRMC, id=G2G.wxID_RUNRMC)
        G2frame.Bind(wx.EVT_MENU, OnViewRMC, id=G2G.wxID_VIEWRMC)
        G2frame.Bind(wx.EVT_MENU, OnLoadRMC, id=G2G.wxID_ATOMSRMC)
        G2frame.Bind(wx.EVT_MENU, OnLoadRMCsuper, id=G2G.wxID_SUPERRMC)
        # ISODISTORT
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.ISODData)
        G2frame.Bind(wx.EVT_MENU, OnRunISODISTORT, id=G2G.wxID_ISODISTORT)
        G2frame.Bind(wx.EVT_MENU, OnNewISOPhase, id=G2G.wxID_ISODNEWPHASE)
        G2frame.Bind(wx.EVT_MENU, OnNewPDFfitPhase, id=G2G.wxID_ISOPDFFIT)
        G2frame.Bind(wx.EVT_MENU, OnShowIsoDistortCalc, id=G2G.wxID_SHOWISO1)
        G2frame.Bind(wx.EVT_MENU, OnShowIsoModes, id=G2G.wxID_SHOWISOMODES)
        # MC/SA
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.MCSAMenu)
        G2frame.Bind(wx.EVT_MENU, OnMCSAaddAtom, id=G2G.wxID_ADDMCSAATOM)
        G2frame.Bind(wx.EVT_MENU, OnMCSAaddRB, id=G2G.wxID_ADDMCSARB)
        G2frame.Bind(wx.EVT_MENU, OnMCSAclear, id=G2G.wxID_CLEARMCSARB)
        G2frame.Bind(wx.EVT_MENU, OnMCSAmove, id=G2G.wxID_MOVEMCSA)
        G2frame.Bind(wx.EVT_MENU, OnClearResults, id=G2G.wxID_MCSACLEARRESULTS)
        # Texture
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.TextureMenu)
        G2frame.Bind(wx.EVT_MENU, OnTextureRefine, id=G2G.wxID_REFINETEXTURE)
#        G2frame.Bind(wx.EVT_MENU, OnTextureClear, id=G2G.wxID_CLEARTEXTURE)
        # Pawley reflections
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.PawleyMenu)
        G2frame.Bind(wx.EVT_MENU, OnPawleySet, id=G2G.wxID_PAWLEYSET)
        G2frame.Bind(wx.EVT_MENU, OnPawleyLoad, id=G2G.wxID_PAWLEYLOAD)
        G2frame.Bind(wx.EVT_MENU, OnPawleyEstimate, id=G2G.wxID_PAWLEYESTIMATE)
        G2frame.Bind(wx.EVT_MENU, OnPawleyUpdate, id=G2G.wxID_PAWLEYUPDATE)
        G2frame.Bind(wx.EVT_MENU, OnPawleySelAll, id=G2G.wxID_PAWLEYSELALL)
        G2frame.Bind(wx.EVT_MENU, OnPawleySelNone, id=G2G.wxID_PAWLEYSELNONE)
        G2frame.Bind(wx.EVT_MENU, OnPawleyToggle, id=G2G.wxID_PAWLEYSELTOGGLE)

    def rbKeyPress(event):
        '''Respond to a Tab to highlight the next RB or crystal atom
        '''
        if 'testRBObj' not in data: return
        if not RigidBodies.atomsGrid: return
        alt = event.GetModifiers() & wx.MOD_ALT
        event.Skip()
        try: # IsKeyInCategory not in wx 2.9
            if not event.IsKeyInCategory(wx.WXK_CATEGORY_TAB): return
        except:
            return
        if alt:    # advance RB selection
            #GSASIIpath.IPyBreak()
            rows = RigidBodies.atomsGrid.GetSelectedRows()
            if len(rows) == 0:
                rows = [0]
            else:
                rows[0] += 1
            if rows[0] > RigidBodies.atomsGrid.GetNumberRows()-1:
                rows = [0]
            elif rows[0] < 0:
                rows[0] = RigidBodies.atomsGrid.GetNumberRows()-1
            RigidBodies.atomsGrid.SelectRow(rows[0])
            RigidBodies.atomsGrid.MakeCellVisible(rows[0],0)
            data['testRBObj']['RBhighLight'] = rows[0]
        else:
            Ind = data['testRBObj'].get('CRYhighLight',[])
            if len(Ind) == 0:
                I = -1
            else:
                I = Ind[0]
            wrap = False
            while True:
                I += 1
                if I >= len(data['Atoms']) and wrap:
                    print('How did this happen?',Ind,I,
                              len(data['testRBObj']['availAtoms']),len(data['Atoms']))
                    return
                elif I >= len(data['Atoms']):
                    wrap = True
                    I = 0
                if data['Atoms'][I][0] in data['testRBObj']['availAtoms']:
                    data['testRBObj']['CRYhighLight'] = [I]
                    misc['showSelect'].setByString(data['Atoms'][I][0])
                    break
        G2plt.PlotStructure(G2frame,data,False,misc['UpdateTable'])
        G2frame.Raise()
        return

    #### UpdatePhaseData execution starts here
    # make sure that the phase menu bars get created before selecting
    # any (this will only be true on the first call to UpdatePhaseData)
    if callable(G2frame.dataWindow.DataGeneral):
        wx.BeginBusyCursor()
        G2frame.dataWindow.DataGeneral()
        wx.EndBusyCursor()

    #patch
    if 'RBModels' not in data:
        data['RBModels'] = {}
    if 'MCSA' not in data:
        data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
    #if isinstance(data['MCSA']['Results'],dict):
    if 'dict' in str(type(data['MCSA']['Results'])):
        data['MCSA']['Results'] = []
    if 'Modulated' not in data['General']:
        data['General']['Modulated'] = False
    if 'doDysnomia' not in data['General']:
        data['General']['doDysnomia'] = False
    if 'modulated' in data['General']['Type']:
        data['General']['Modulated'] = True
        data['General']['Type'] = 'nuclear'
    if 'RMC' not in data:
        data['RMC'] = {'RMCProfile':{},'fullrmc':{},'PDFfit':{}}
    if 'ISODISTORT' not in data:
        data['ISODISTORT'] = {}
    if 'Deformations' not in data:
        data['Deformations'] = {}
#end patch

    global rbAtmDict
    rbAtmDict = {}
    misc = {}
    PhaseName = G2frame.GPXtree.GetItemText(Item)
    G2gd.SetDataMenuBar(G2frame)
    G2frame.phaseDisplay = G2G.GSNoteBook(parent=G2frame.dataWindow)
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(G2frame.phaseDisplay,1,wx.ALL|wx.EXPAND,1)
    G2frame.phaseDisplay.gridList = [] # list of all grids in notebook
    Pages = []
    General = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(General,'General')
    Pages.append('General')
    if not GSASIIpath.GetConfigValue('SeparateHistPhaseTreeItem',False):
        DData = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(DData,'Data')
        Pages.append('Data')
    AtomList = wx.Panel(G2frame.phaseDisplay)
    Atoms = G2G.GSGrid(AtomList)
    G2frame.phaseDisplay.gridList.append(Atoms)
    G2frame.phaseDisplay.AddPage(AtomList,'Atoms')
    Pages.append('Atoms')
    if data['General']['Modulated']:
        G2frame.waveData = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.waveData,'Wave Data')
        Pages.append('Wave Data')
    if data['General']['Type'] == 'faulted':
        G2frame.layerData = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.layerData,'Layers')
        Pages.append('Layers')
    G2frame.drawOptions = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(G2frame.drawOptions,'Draw Options')
    Pages.append('Draw Options')
    G2frame.drawAtomsList = wx.Panel(G2frame.phaseDisplay)
    G2frame.drawAtoms = G2G.GSGrid(G2frame.drawAtomsList)
    G2frame.phaseDisplay.gridList.append(G2frame.drawAtoms)
    G2frame.phaseDisplay.AddPage(G2frame.drawAtomsList,'Draw Atoms')
    Pages.append('Draw Atoms')

    dmin,dmax,nhist,lbl,chist = getPawleydRange(G2frame,data)
    dataGeneral = data['General']
    if any('X' in item for item in G2frame.GetHistogramTypes()):
        G2frame.deformation = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.deformation,'Deformation')
        Pages.append('Deformation')
        
    if dataGeneral['Type'] not in ['faulted',] and not dataGeneral['Modulated']:
        RigidBodies = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(RigidBodies,'RB Models')
        # note the bind is here so that it is only done once, but
        # TODO: might need to be on a different widget for Windows
        RigidBodies.Bind(wx.EVT_CHAR,rbKeyPress)
        Pages.append('RB Models')

    if ('Map' in dataGeneral and dataGeneral['Map']['MapType']) or 'Flip' in dataGeneral:
        MapPeakList = wx.Panel(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(MapPeakList,'Map peaks')
        # create the grid once; N.B. need to reference at this scope
        G2frame.MapPeaks = G2G.GSGrid(MapPeakList)
        G2frame.phaseDisplay.gridList.append(G2frame.MapPeaks)
        Pages.append('Map peaks')

        if data['General']['doDysnomia']:
            G2frame.MEMData = wx.ScrolledWindow(G2frame.phaseDisplay)
            G2frame.phaseDisplay.AddPage(G2frame.MEMData,'Dysnomia')
            Pages.append('Dysnomia')
        
    if dataGeneral['Type'] not in ['faulted',] and not dataGeneral['Modulated']:
        G2frame.MCSA = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.MCSA,'MC/SA')
        Pages.append('MC/SA')
        
    if nhist:   
        G2frame.FRMC = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.FRMC,'RMC')
        Pages.append('RMC')

    if dataGeneral['Type'] == 'nuclear':
        G2frame.ISODIST = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.ISODIST,'ISODISTORT')
        Pages.append('ISODISTORT')
        
    if nhist > 0:            
        G2frame.Texture = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.Texture,'Texture')
        Pages.append('Texture')
        PawleyRefList = wx.Panel(G2frame.phaseDisplay)
        G2frame.PawleyRefl = None  # grid now created when needed
        G2frame.phaseDisplay.AddPage(PawleyRefList,'Pawley reflections')
        Pages.append('Pawley reflections')
    G2frame.dataWindow.AtomCompute.Enable(G2G.wxID_ISODISP,'ISODISTORT' in data)
    G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_VALIDPROTEIN,'macro' in data['General']['Type'])
    G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_USEBILBAOMAG,'magPhases' in data)
    flag = False
    if 'magPhases' in data:
        PatternName = data['magPhases']
        PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PatternName)
        if (PatternId):
            UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List')
            UCdata = list(G2frame.GPXtree.GetItemPyData(UnitCellsId))
            if len(UCdata) >5:
                flag = not any(['magAtms' in i for i in UCdata[5]])
        else:
            del data['magPhases']
    G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_USEBILBAOSUB,flag)
    if not flag and 'SUBGROUPS' in data:
        G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_USEBILBAOSUB,True)
        G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_USEBILBAOMAG,True)
    G2frame.phaseDisplay.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    FillMenus()
    if G2frame.lastSelectedPhaseTab in Pages:
        ind = Pages.index(G2frame.lastSelectedPhaseTab)
        if ind:
            UpdateGeneral(SkipDraw=ind)
            #G2frame.phaseDisplay.SetSelection(ind)
            # on windows, need to wait for previous to finish before
            # going to selected tab
            wx.CallAfter(G2frame.phaseDisplay.SetSelection,ind)
            return
    ChangePage(0)
    
def CheckAddHKLF(G2frame,data):
    '''GUI actions associated with linking a Phase to a HKLF histogram.

    This gets called in two routines named OnHklfAdd (one inside
    :func:`GSASIIphsGUI.UpdatePhaseData` and the other inside
    :func:`GSASIIddataGUI.MakeHistPhaseWin`).
    '''
    keyList = list(data['Histograms'].keys())
    TextList = []
    if not G2frame.GPXtree.GetCount():
        return

    item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
    while item:
        name = G2frame.GPXtree.GetItemText(item)
        if name not in keyList and 'HKLF' in name:
            TextList.append(name)
        item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
    if not TextList:
        G2G.G2MessageBox(G2frame,'No HKLF histograms')
        return
    dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select HKLF reflection sets to use',
            'Use data',TextList)
    try:
        if dlg.ShowModal() == wx.ID_OK:
            result = dlg.GetSelections()
        else:
            print('Nothing selected')
            return
    finally:
        dlg.Destroy()

    # get the histograms used in other phases
    phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
    usedHKLFhists = [] # used single-crystal histograms
    for p in usedHistograms:
        for h in usedHistograms[p]:
            if h.startswith('HKLF ') and h not in usedHKLFhists:
                usedHKLFhists.append(h)
    # check that selected single crystal histograms are not already in use!
    for i in result:
        used = [TextList[i] for i in result if TextList[i] in usedHKLFhists]
        if used:
          msg = 'The following single crystal histogram(s) are already in use'
          for i in used:
              msg += '\n  '+str(i)
          msg += '\nAre you sure you want to add them to this phase? '
          msg += 'Associating a single crystal dataset to >1 histogram is usually an error, '
          msg += 'so No is suggested here.'
          if G2frame.ErrorDialog('Likely error',msg,G2frame,wtype=wx.YES_NO) != wx.ID_YES:
            return

    wx.BeginBusyCursor()
    for i in result:
        histoName = TextList[i]
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
        refDict,reflData = G2frame.GPXtree.GetItemPyData(Id)[:2]
        G2mth.UpdateHKLFvals(histoName, data, reflData)

    wx.EndBusyCursor()
    return result

def checkPDFfit(G2frame):
    '''Checks to see if PDFfit2 is available and can be imported. PDFfit2 can be installed
    in a separate Python interpreter (saved in the pdffit2_exec config variable). If this is
    defined, no attempt is made to check that it actually runs.
    Otherwise, if diffpy.PDFfit has been installed with conda/pip, it is checked if the
    install command.

    :returns: False if PDFfit2 cannot be run/accessed. True if it appears it can be run.
    '''
    # if a separate Python interpreter has been specified, just use it, no checking
    if GSASIIpath.GetConfigValue('pdffit2_exec') is not None and is_exe(
            GSASIIpath.GetConfigValue('pdffit2_exec')):
        return True

    # see if diffpy has been installed directly
    try:
        from diffpy.pdffit2 import PdfFit
        PdfFit
        return True
    except:
        pass

    # Last effort: With conda we should be able to create a separate
    # Python in a separate environment
    try:     # have conda. Can we access it programmatically?
        import conda.cli.python_api
        conda.cli.python_api
    except:
        G2G.G2MessageBox(G2frame,'You are running a directly installed Python. You will need to install PDFfit2 directly as well, preferably in a separate virtual environment.')
        return

    msg = ('Do you want to use conda to install PDFfit2 into a separate environment? '+
               '\n\nIf successful, the pdffit2_exec configuration option will be set to the '+
               'this new Python environment.')
    dlg = wx.MessageDialog(G2frame,msg,caption='Install?',
                                   style=wx.YES_NO|wx.ICON_QUESTION)
    if dlg.ShowModal() != wx.ID_YES:
        return False
    try:
        wx.BeginBusyCursor()
        print('Preparing to create a conda environment. This may take a few minutes...')
        # for now use the older diffpy version of pdffit:
        #   conda create -n pdffit2 python=3.7 conda gsl diffpy.pdffit2=1.3.4 -c conda-forge -c diffpy
        res,PDFpython = GSASIIpath.condaEnvCreate('pdffit2',
                    ['python', 'conda', 'gsl', 'diffpy.pdffit2>=1.4.3',
                         '-c', 'conda-forge']) #  not needed , '-c', 'diffpy'])
    finally:
        wx.EndBusyCursor()
    if os.path.exists(PDFpython) and is_exe(PDFpython):
        vars = G2G.GetConfigValsDocs()
        vars['pdffit2_exec'][1] = PDFpython
        GSASIIpath.SetConfigValue(vars)
        G2G.SaveConfigVars(vars)
        print('pdffit2_exec config set with ',GSASIIpath.GetConfigValue('pdffit2_exec'))
        print('\n\nSuccess: PDFfit2 installed.')
        return True
    else:
        print(f'Failed to install PDFfit2 with error:\n{PDFpython}')
        if ('PackagesNotFoundError' in PDFpython
                and 'darwin' in sys.platform
                and 'arm' in platform.machine()):
            msg = ('It appears that PDFfit2 is not yet available as a conda package for Macs with arm processors. '+
                       '\n\nYou could install PDFfit2 with x86 Python and use that '+
                       'in compatibility mode.')
        else:
            msg = ('An attempt to install PDFfit2 has failed. '+
                       'Do you have write access to where GSAS-II is installed? '+
                       'You may be able to install PDFfit2 manually.')
        msg += '\n\nIf you install PDFfit2 yourself, set the pdffit2_exec config variable to the install location'
        G2G.G2MessageBox(G2frame,
                'PDFfit2 Install failed. See console for error message\n\n'+msg,
                'PDFfit2 install error')
        return False

def makeIsoNewPhase(phData,cell,atomList,sglbl,sgnum):
    '''create a new phase from a supergroup structure generated by ISOCIF
    '''
    if len(atomList) == 0:
        print(sglbl,'empty structure')
        return
    # create a new phase
    try:
        sgnum = int(sgnum)
        sgsym = G2spc.spgbyNum[sgnum]
        #sgname = sgsym.replace(" ","")
    except:
        print(f'Problem with processing space group name {sglbl} and number {sgnum}')
        return
    newPhase = copy.deepcopy(phData)
    newPhase['ranId'] = ran.randint(0,sys.maxsize),
    if 'magPhases' in phData: del newPhase['magPhases']
    generalData = newPhase['General']
    generalData['SGData'] = SGData = G2spc.SpcGroup(sgsym)[1]
    generalData['Cell'][1:7] = cell
    generalData['Cell'][7] = G2lat.calc_V(G2lat.cell2A(generalData['Cell'][1:7]))
    cx,ct,cs,cia = generalData['AtomPtrs']
    Atoms = newPhase['Atoms'] = []
    for nam,(x,y,z) in atomList:
        try:
            atom = []
            atom.append(nam)
            if nam[1].isdigit():
                atom.append(nam[0:1])
            else:
                atom.append(nam[0:2])
            atom.append('')
            for i in x,y,z: atom.append(float(i))
            atom.append(1.0)
            SytSym,Mult = G2spc.SytSym(np.array(atom[3:6]),SGData)[:2]
            atom.append(SytSym)
            atom.append(Mult)
            atom.append('I')
            atom += [0.02,0.,0.,0.,0.,0.,0.,]
            atom.append(ran.randint(0,sys.maxsize))
            Atoms.append(atom)
        except Exception as msg:
            print(f'{msg} error in atom {nam}')
    G2elem.SetupGeneral(newPhase,generalData['Mydir'])  # fixup composition info
    return newPhase

def saveIsoNewPhase(G2frame,phData,newPhase,orgFilName):
    '''save the new phase generated by ISOCIF created in :func:`makeIsoNewPhase`
    into a GSAS-II project (.gpx) file
    '''
    import re
    phData.update(newPhase)
    # save new file
    sgname = newPhase['General']['SGData']['SpGrp'].replace(' ','')
    G2frame.GSASprojectfile = os.path.splitext(orgFilName
                            )[0]+'_super_'+sgname.replace('/','$')+'.gpx'
    while os.path.exists(G2frame.GSASprojectfile):
        s = re.split(r'_([\d]+)\.gpx',G2frame.GSASprojectfile)
        if len(s) == 1:
            G2frame.GSASprojectfile = os.path.splitext(G2frame.GSASprojectfile)[0] + '_1.gpx'
        else:
            num = 10
            try:
                num = int(s[1]) + 1
            except:
                pass
            G2frame.GSASprojectfile = f'{s[0]}_{num}.gpx'
    G2IO.ProjFileSave(G2frame)
    return G2frame.GSASprojectfile

def renamePhaseName(G2frame,data,phaseItem,generalData,newName):
        '''Called to rename the phase. Updates the tree and items that
        reference the file name. 
        '''
        oldName = generalData['Name']
        phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if newName and newName != oldName:
            newName = G2obj.MakeUniqueLabel(newName,list(phaseNameList))
            generalData['Name'] = newName
            G2frame.G2plotNB.Rename(oldName,generalData['Name'])
            G2frame.GPXtree.SetItemText(phaseItem,generalData['Name'])
            # change phase name key in Reflection Lists for each histogram
            for hist in data['Histograms']:
                ht = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hist)
                rt = G2gd.GetGPXtreeItemId(G2frame,ht,'Reflection Lists')
                if not rt: continue
                RfList = G2frame.GPXtree.GetItemPyData(rt)
                if oldName not in RfList:
                    print('Warning: '+oldName+' not in Reflection List for '+
                          hist)
                    continue
                RfList[newName] = RfList[oldName]
                del RfList[oldName]
            # rename Restraints
            resId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
            Restraints = G2frame.GPXtree.GetItemPyData(resId)
            i = G2gd.GetGPXtreeItemId(G2frame,resId,oldName)
            if i: G2frame.GPXtree.SetItemText(i,newName)
            if len(Restraints) and oldName in Restraints:
                Restraints[newName] = Restraints[oldName]
                del Restraints[oldName]
