# -*- coding: utf-8 -*-
#GSASII - phase data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIphsGUI: Phase GUI*
-------------------------

Module to create the GUI for display of phase information
in the data display window when a phase is selected.
Phase information is stored in one or more
:ref:`Phase Tree Item <Phase_table>` objects.
Note that there are functions
that respond to some tabs in the phase GUI in other modules
(such as GSASIIddata).

Main routine here is :func:`UpdatePhaseData`, which displays the phase information
(called from :func:`GSASIIdataGUI:SelectDataTreeItem`).

Other top-level routines are: 
:func:`GetSpGrpfromUser` (called locally only);
:func:`FindBondsDraw` and :func:`FindBondsDrawCell` (called locally and in GSASIIplot); 
:func:`SetPhaseWindow` (called locally and in GSASIIddataGUI and GSASIIrestrGUI, multiple locations)
to control scrolling. 
'''
from __future__ import division, print_function
import os.path
import wx
import wx.grid as wg
import wx.lib.scrolledpanel as wxscroll
import matplotlib as mpl
import math
import copy
import time
import sys
import random as ran
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIddataGUI as G2ddG
import GSASIIplot as G2plt
import GSASIIdataGUI as G2gd
import GSASIIIO as G2IO
import GSASIIstrMain as G2stMn
import GSASIIstrIO as G2strIO
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import GSASIIpy3 as G2py3
import GSASIIobj as G2obj
import GSASIIctrlGUI as G2G
import GSASIIconstrGUI as G2cnstG
import numpy as np
import numpy.linalg as nl

VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)
BLACK = wx.Colour(0,0,0)
WACV = wx.ALIGN_CENTER_VERTICAL
mapDefault = {'MapType':'','RefList':'','Resolution':0.5,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}
TabSelectionIdDict = {}
# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda x,y: 180.*np.arctan2(y,x)/np.pi

################################################################################
#### phase class definitions
################################################################################
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
            mainSizer.Add(self.force,0,WACV|wx.TOP,5)
#        if SGData['SGInv']:
        choice = ['No','Yes']
        self.inv = wx.RadioBox(panel,-1,'Choose inversion?',choices=choice)
        self.inv.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.inv,0,WACV)
        if SGData['SGLatt'] != 'P':
            LattOp = G2spc.Latt2text(SGData['SGCen']).split(';')
            self.latt = wx.RadioBox(panel,-1,'Choose cell centering?',choices=LattOp)
            self.latt.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.latt,0,WACV)
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
        mainSizer.Add(self.oprs,0,WACV|wx.BOTTOM,5)
        mainSizer.Add(wx.StaticText(panel,-1,"   Choose unit cell?"),0,WACV)
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellName = ['X','Y','Z']
        self.cell = []
        for i in range(3):
            self.cell.append(wx.SpinCtrl(panel,-1,cellName[i],size=wx.Size(50,20)))
            self.cell[-1].SetRange(-3,3)
            self.cell[-1].SetValue(0)
            self.cell[-1].Bind(wx.EVT_SPINCTRL, self.OnOpSelect)
            cellSizer.Add(self.cell[-1],0,WACV)
        mainSizer.Add(cellSizer,0,WACV|wx.BOTTOM,5)
        if self.New:
            choice = ['No','Yes']
            self.new = wx.RadioBox(panel,-1,'Generate new positions?',choices=choice)
            self.new.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.new,0,WACV)

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
#        if self.SGData['SGInv']:
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
################################################################################
class SphereEnclosure(wx.Dialog):
    ''' Add atoms within sphere of enclosure to drawing
    
    :param wx.Frame parent: reference to parent frame (or None)
    :param general: general data (includes drawing data)
    :param atoms: drawing atoms data
    :param indx: list of selected atoms (may be empty)
    
    '''
    def __init__(self,parent,general,drawing,indx):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Setup phase transformation', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.General = general
        self.Drawing = drawing
        self.indx = indx
        self.Sphere = [1.0,]
        self.centers = []
        self.atomTypes = [[item,True] for item in self.General['AtomTypes']]
        
        self.Draw()
        
    def Draw(self):
        
        def OnAtomType(event):
            Obj = event.GetEventObject()
            id = Ind[Obj.GetId()]
            self.atomTypes[id][1] = Obj.GetValue()
        
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label=' Sphere of enclosure controls:'),0,WACV)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        atoms = []
        if len(self.indx):
            topSizer.Add(wx.StaticText(self.panel,label=' Sphere centered at atoms: '),0,WACV)
            cx,ct,cs = self.Drawing['atomPtrs'][:3]
            for id in self.indx:
                atom = self.Drawing['Atoms'][id]
                self.centers.append(atom[cx:cx+3])
                atoms.append('%s(%s)'%(atom[ct-1],atom[cs-1]))
            topSizer.Add(wx.ComboBox(self.panel,choices=atoms,value=atoms[0],
                style=wx.CB_READONLY|wx.CB_DROPDOWN),0,WACV)
        else:
            topSizer.Add(wx.StaticText(self.panel,label=' Sphere centered at drawing view point'),0,WACV)
            self.centers.append(self.Drawing['viewPoint'][0])
        mainSizer.Add(topSizer,0,WACV)
        sphereSizer = wx.BoxSizer(wx.HORIZONTAL)
        sphereSizer.Add(wx.StaticText(self.panel,label=' Sphere radius: '),0,WACV)
        radius = G2G.ValidatedTxtCtrl(self.panel,self.Sphere,0,nDig=(10,3),size=(65,25))
        sphereSizer.Add(radius,0,WACV)
        mainSizer.Add(sphereSizer,0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' Target selected atoms:'),0,WACV)
        atSizer = wx.BoxSizer(wx.HORIZONTAL)
        Ind = {}
        for i,item in enumerate(self.atomTypes):
            atm = wx.CheckBox(self.panel,label=item[0])
            atm.SetValue(item[1])
            atm.Bind(wx.EVT_CHECKBOX, OnAtomType)
            Ind[atm.GetId()] = i
            atSizer.Add(atm,0,WACV)
        mainSizer.Add(atSizer,0,WACV)
        
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

################################################################################
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
        self.oldSpGrp = phase['General']['SGData']['SpGrp']
        self.oldSGdata = phase['General']['SGData']
        self.newSpGrp = self.Phase['General']['SGData']['SpGrp']
        self.SGData = G2spc.SpcGroup(self.newSpGrp)[1]
        self.oldCell = phase['General']['Cell'][1:8]
        self.newCell = self.Phase['General']['Cell'][1:8]
        self.Common = 'abc'
        self.ifMag = ifMag
        if ifMag:
            self.BNSlatt = BNSlatt
        self.ifConstr = True
        self.Mtrans = False
        self.kvec = [0.,0.,0.]
        self.Draw()

    def Draw(self):
                
        def OnCommon(event):
            Obj = event.GetEventObject()
            self.Common = Obj.GetValue()
            self.Mtrans = False
            if '*' in self.Common:
                A,B = G2lat.cell2AB(self.oldCell[:6])
                self.newCell[2:5] = [A[2,2],90.,90.]
                a,b = G2lat.cell2AB(self.newCell[:6])
                self.Trans = np.inner(a.T,B).T    #correct!
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
            if self.Phase['General']['Type'] == 'magnetic':
                Nops = len(self.SGData['SGOps'])*len(self.SGData['SGCen'])
                if self.SGData['SGInv']:
                    Nops *= 2
                self.SGData['SpnFlp'] = Nops*[1,]
                del self.oldSGdata['MAXMAGN']
            wx.CallAfter(self.Draw)

        def OnTest(event):
            if self.Mtrans:
                self.newCell = G2lat.TransformCell(self.oldCell[:6],self.Trans.T)
            else:
                self.newCell = G2lat.TransformCell(self.oldCell[:6],self.Trans)
            wx.CallAfter(self.Draw)
            
        def OnMag(event):
            self.ifMag = True
            self.BNSlatt = self.SGData['SGLatt']
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
            mainSizer.Add(mag,0,WACV)            
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
        Trmat = wx.FlexGridSizer(4,6,0,0)
        Trmat.Add((10,0),0)
        Trmat.Add(wx.StaticText(self.panel,label='      M'),wx.ALIGN_CENTER)
        Trmat.Add((10,0),0)
        Trmat.Add((10,0),0)
        Trmat.Add(wx.StaticText(self.panel,label='      U'),wx.ALIGN_CENTER)
        Trmat.Add(wx.StaticText(self.panel,label='      V'),wx.ALIGN_CENTER)
        
        for iy,line in enumerate(self.Trans):
            for ix,val in enumerate(line):
                item = G2G.ValidatedTxtCtrl(self.panel,self.Trans[iy],ix,nDig=(10,3),size=(65,25))
                Trmat.Add(item)
            Trmat.Add((25,0),0)
            vec = G2G.ValidatedTxtCtrl(self.panel,self.Uvec,iy,nDig=(10,3),size=(65,25))
            Trmat.Add(vec)
            vec = G2G.ValidatedTxtCtrl(self.panel,self.Vvec,iy,nDig=(10,3),size=(65,25))
            Trmat.Add(vec)
        transSizer.Add(Trmat)
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
            mainSizer.Add(MagSizer,0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' Old lattice parameters:'),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=
            ' a = %.5f       b = %.5f      c = %.5f'%(self.oldCell[0],self.oldCell[1],self.oldCell[2])),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' alpha = %.3f beta = %.3f gamma = %.3f'%
            (self.oldCell[3],self.oldCell[4],self.oldCell[5])),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' volume = %.3f'%(self.oldCell[6])),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' New lattice parameters:'),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=
            ' a = %.5f       b = %.5f      c = %.5f'%(self.newCell[0],self.newCell[1],self.newCell[2])),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' alpha = %.3f beta = %.3f gamma = %.3f'%
            (self.newCell[3],self.newCell[4],self.newCell[5])),0,WACV)
        mainSizer.Add(wx.StaticText(self.panel,label=' volume = %.3f'%(self.newCell[6])),0,WACV)
        sgSizer = wx.BoxSizer(wx.HORIZONTAL)
        sgSizer.Add(wx.StaticText(self.panel,label=' Target space group: '),0,WACV)
        SGTxt = wx.Button(self.panel,wx.ID_ANY,self.newSpGrp,size=(100,-1))
        SGTxt.Bind(wx.EVT_BUTTON,OnSpaceGroup)
        sgSizer.Add(SGTxt,0,WACV)
        mainSizer.Add(sgSizer,0,WACV)
        if 'magnetic' not in self.Phase['General']['Type']:
            if self.ifMag:
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
                mainSizer.Add(BNSizer,0,WACV)
            if self.ifMag:
                mainSizer.Add(wx.StaticText(self.panel, \
                    label=' NB: Nonmagnetic atoms will be deleted from new phase'),0,WACV)
            constr = wx.CheckBox(self.panel,label=' Make constraints between phases?')
            constr.SetValue(self.ifConstr)
            constr.Bind(wx.EVT_CHECKBOX,OnConstr)
            mainSizer.Add(constr,0,WACV)
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
        
    def GetSelection(self):
        self.Phase['General']['SGData'] = self.SGData
        if self.ifMag:
            self.Phase['General']['Name'] += ' mag: '
        else:
            self.Phase['General']['Name'] += ' %s'%(self.Common)
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
        
################################################################################
class UseMagAtomDialog(wx.Dialog):
    '''Get user selected magnetic atoms after cell transformation
    '''
    def __init__(self,parent,Name,Atoms,atCodes,atMxyz,ifDelete=False):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Magnetic atom selection', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.Name = Name
        self.Atoms = Atoms
        self.atCodes = atCodes
        self.atMxyz = atMxyz
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
        self.panel = wx.Panel(self)
        Indx = {}
        Mstr = [' Mx',' My',' Mz']
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,label='For: %s'%self.Name),0,WACV)
        
        mainSizer.Add(wx.StaticText(self.panel,label='        Name, x, y, z, allowed moments, mag. site sym:'),0,WACV)
        atmSizer = wx.FlexGridSizer(0,2,5,5)
        for iuse,[use,atom,mxyz] in enumerate(zip(self.Use,self.Atoms,self.atMxyz)):
            mstr = [' ---',' ---',' ---']
            for i,mx in enumerate(mxyz[1]):
                if mx:
                    mstr[i] = Mstr[i]
            useChk = wx.CheckBox(self.panel,label='Use?')
            Indx[useChk.GetId()] = iuse
            useChk.SetValue(use)
            useChk.Bind(wx.EVT_CHECKBOX, OnUseChk)
            atmSizer.Add(useChk,0,WACV)
            text = '  %5s %10.5f %10.5f %10.5f (%s,%s,%s) %s   '%(atom[0],atom[3],atom[4],atom[5],mstr[0],mstr[1],mstr[2],mxyz[0])
            atmSizer.Add(wx.StaticText(self.panel,label=text),0,WACV)
        mainSizer.Add(atmSizer)
        
        YesBtn = wx.Button(self.panel,-1,"Yes")
        YesBtn.Bind(wx.EVT_BUTTON, self.OnYes)
        NoBtn = wx.Button(self.panel,-1,"No")
        NoBtn.Bind(wx.EVT_BUTTON, self.OnNo)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
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
        self.panel.Fit()
        self.Fit()
        
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
            
                
################################################################################
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
        rotangle = wx.TextCtrl(self.panel,value='%5.3f'%(self.rotAngle),
            size=(50,25),style=wx.TE_PROCESS_ENTER)
        rotangle.Bind(wx.EVT_TEXT_ENTER,OnRotAngle)
        rotangle.Bind(wx.EVT_KILL_FOCUS,OnRotAngle)
        rotationBox.Add(rotangle,0,WACV)
        rotationBox.Add(wx.StaticText(self.panel,label=' about vector: '),0,WACV)
        rotvec = wx.TextCtrl(self.panel,value='%5.3f %5.3f %5.3f'%(self.rotVec[0],self.rotVec[1],self.rotVec[2]),
            size=(100,25),style=wx.TE_PROCESS_ENTER)
        rotvec.Bind(wx.EVT_TEXT_ENTER,OnRotVec)
        rotvec.Bind(wx.EVT_KILL_FOCUS,OnRotVec)
        rotationBox.Add(rotvec,0,WACV)
        mainSizer.Add(rotationBox,0,WACV)
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
        
################################################################################
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
        mainSizer.Add(wx.StaticText(self.panel,label=' Controls for DIFFaX'),0,WACV)
        if self.Parms:
            mainSizer.Add(wx.StaticText(self.panel,label=' Sequential powder pattern simulation'),0,WACV)
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
            mainSizer.Add(wx.StaticText(self.panel,label=' Enter parameter range & no. steps: '),0,WACV)
            parmRange =  wx.BoxSizer(wx.HORIZONTAL)
            numChoice = [str(i+1) for i in range(10)]
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

################################################################################
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
        mainSizer.Add(wx.StaticText(self.panel,-1,'NB: Check selections as they may not be correct'),0,WACV|wx.LEFT,10)
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

################################################################################
################################################################################
################################################################################
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
    
    
def SetPhaseWindow(phasePage,mainSizer=None,Scroll=0):
    if mainSizer is not None:
        phasePage.SetSizer(mainSizer)
    phasePage.SetAutoLayout(True)
    phasePage.SetScrollRate(10,10)
    phasePage.SendSizeEvent()
    phasePage.Scroll(0,Scroll)
    
def FindBondsDraw(data):
    '''Generally used routine where cell is from data
    '''
    generalData = data['General']
    cell = generalData['Cell'][1:7]
    FindBondsDrawCell(data,cell)
    
def FindBondsDrawCell(data,cell):    
    '''uses numpy & masks - very fast even for proteins!
    allows different cell as input from seq. refinements
    '''
    import numpy.ma as ma
    cx,ct,cs,ci = data['Drawing']['atomPtrs']
    hydro = data['Drawing']['showHydrogen']
    atomData = data['Drawing']['Atoms']
    generalData = data['General']
    Amat,Bmat = G2lat.cell2AB(cell)
    radii = generalData['BondRadii']
#    if generalData.get('DisAglCtls',{}):
#        radii = generalData['DisAglCtls']['BondRadii']
    atomTypes = generalData['AtomTypes']
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
                        
################################################################################
################################################################################
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
            'bondList':{},'viewDir':[1,0,0],'Plane':[[0,0,1],False,False,0.0,[255,255,0]]}
    for key in defaultDrawing:
        if key not in drawingData: drawingData[key] = defaultDrawing[key]
            
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
    def GetReflData(G2frame,phaseName,reflNames):
        ReflData = {'RefList':[],'Type':''}
        if '' in reflNames:
            return None
        for reflName in reflNames:
            if 'PWDR' in reflName:
                PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
                reflSets = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Reflection Lists'))
                reflData = reflSets[phaseName]
            elif 'HKLF' in reflName:
                PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root, reflName)
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

    def SetupGeneral():
        generalData = data['General']
        atomData = data['Atoms']
        generalData['AtomTypes'] = []
        generalData['Isotopes'] = {}
# various patches
        if 'Isotope' not in generalData:
            generalData['Isotope'] = {}
        if 'Data plot type' not in generalData:
            generalData['Data plot type'] = 'Mustrain'
        if 'POhkl' not in generalData:
            generalData['POhkl'] = [0,0,1]
        if 'Map' not in generalData:
            generalData['Map'] = mapDefault.copy()
        if 'Flip' not in generalData:
            generalData['Flip'] = {'RefList':'','Resolution':0.5,'Norm element':'None',
                'k-factor':0.1,'k-Max':20.,}
        if 'testHKL' not in generalData['Flip']:
            generalData['Flip']['testHKL'] = [[0,0,2],[2,0,0],[1,1,1],[0,2,0],[1,2,3]]
        if 'doPawley' not in generalData:
            generalData['doPawley'] = False     #ToDo: change to ''
        if 'Pawley dmin' not in generalData:
            generalData['Pawley dmin'] = 1.0
        if 'Pawley dmax' not in generalData:
            generalData['Pawley dmax'] = 100.0
        if 'Pawley neg wt' not in generalData:
            generalData['Pawley neg wt'] = 0.0
        if '3Dproj' not in generalData:
            generalData['3Dproj'] = ''
        if 'Algolrithm' in generalData.get('MCSA controls',{}) or \
            'MCSA controls' not in generalData:
            generalData['MCSA controls'] = {'Data source':'','Annealing':[0.7,0.1,250],
            'dmin':2.8,'Algorithm':'log','fast parms':[0.8,0.6],'log slope':0.9,
            'Cycles':1,'Results':[],'newDmin':True}
        if 'AtomPtrs' not in generalData:
            generalData['AtomPtrs'] = [3,1,7,9]
            if generalData['Type'] == 'macromolecular':
                generalData['AtomPtrs'] = [6,4,10,12]
            elif generalData['Type'] == 'magnetic':
                generalData['AtomPtrs'] = [3,1,10,12]
        if generalData['Modulated']:
            if 'Super' not in generalData:
                generalData['Super'] = 1
                generalData['SuperVec'] = [[0.,0.,0.],False,4]
                generalData['SSGData'] = {}
            if '4DmapData' not in generalData:
                generalData['4DmapData'] = mapDefault.copy()
                generalData['4DmapData'].update({'MapType':'Fobs'})
            atomData = data['Atoms']
            for atom in atomData:
                if 'SS1' not in atom:
                    atom += [[],[],{'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}]
                if 'waveType' in atom[-1]['SS1']:
                    waveType = atom[-1]['SS1']['waveType']
                    for parm in ['Sfrac','Spos','Sadp','Smag']:
                        if len(atom[-1]['SS1'][parm]):
                            wType = 'Fourier'
                            if parm == 'Sfrac':
                                if 'Crenel' in waveType:
                                    wType = 'Crenel'
                            elif parm == 'Spos':
                                if not 'Crenel' in waveType:
                                    wType = waveType
                            atom[-1]['SS1'][parm] = [wType,]+list(atom[-1]['SS1'][parm])
                    del atom[-1]['SS1']['waveType']
        if 'Modulated' not in generalData:
            generalData['Modulated'] = False
        if 'HydIds' not in generalData:
            generalData['HydIds'] = {}
        if generalData['Type'] == 'magnetic':
            if 'SGGray' not in generalData['SGData']:
                generalData['SGData']['SGGray'] = False
                
# end of patches
        cx,ct,cs,cia = generalData['AtomPtrs']
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['vdWRadii'] = []
        generalData['AtomMass'] = []
        generalData['Color'] = []
        if generalData['Type'] == 'magnetic':
            generalData['MagDmin'] = generalData.get('MagDmin',1.0)
            landeg = generalData.get('Lande g',[])
        generalData['Mydir'] = G2frame.dirname
        badList = {}
        for iat,atom in enumerate(atomData):
            atom[ct] = atom[ct].lower().capitalize()              #force to standard form
            if generalData['AtomTypes'].count(atom[ct]):
                generalData['NoAtoms'][atom[ct]] += atom[cx+3]*float(atom[cs+1])
            elif atom[ct] != 'UNK':
                Info = G2elem.GetAtomInfo(atom[ct])
                if not Info:
                    if atom[ct] not in badList:
                        badList[atom[ct]] = 0
                    badList[atom[ct]] += 1
                    atom[ct] = 'UNK'
                    continue
                atom[ct] = Info['Symbol'] # N.B. symbol might be changed by GetAtomInfo
                generalData['AtomTypes'].append(atom[ct])
                generalData['Z'] = Info['Z']
                generalData['Isotopes'][atom[ct]] = Info['Isotopes']
                generalData['BondRadii'].append(Info['Drad'])
                generalData['AngleRadii'].append(Info['Arad'])
                generalData['vdWRadii'].append(Info['Vdrad'])
                if atom[ct] in generalData['Isotope']:
                    if generalData['Isotope'][atom[ct]] not in generalData['Isotopes'][atom[ct]]:
                        isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                        generalData['Isotope'][atom[ct]] = isotope
                    generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]]['Mass'])
                else:
                    generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                    if 'Nat. Abund.' not in generalData['Isotopes'][atom[ct]]:
                        isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                        generalData['Isotope'][atom[ct]] = isotope
                    generalData['AtomMass'].append(Info['Mass'])
                generalData['NoAtoms'][atom[ct]] = atom[cx+3]*float(atom[cs+1])
                generalData['Color'].append(Info['Color'])
                if generalData['Type'] == 'magnetic':
                    if len(landeg) < len(generalData['AtomTypes']):
                        landeg.append(2.0)
        if generalData['Type'] == 'magnetic':
            generalData['Lande g'] = landeg[:len(generalData['AtomTypes'])]
                        
        if badList:
            msg = 'Warning: element symbol(s) not found:'
            for key in badList:
                msg += '\n\t' + key
                if badList[key] > 1:
                    msg += ' (' + str(badList[key]) + ' times)'
            wx.MessageBox(msg,caption='Element symbol error')
        F000X = 0.
        F000N = 0.
        for i,elem in enumerate(generalData['AtomTypes']):
            F000X += generalData['NoAtoms'][elem]*generalData['Z']
            isotope = generalData['Isotope'][elem]
            F000N += generalData['NoAtoms'][elem]*generalData['Isotopes'][elem][isotope]['SL'][0]
        generalData['F000X'] = F000X
        generalData['F000N'] = F000N
        generalData['Mass'] = G2mth.getMass(generalData)
       

################################################################################
##### General phase routines
################################################################################

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
                event.Skip()
                oldName = generalData['Name']
                phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
                phaseNameList = usedHistograms.keys() # phase names in use
                newName = NameTxt.GetValue().strip()
                if newName and newName != oldName:
                    newName = G2obj.MakeUniqueLabel(newName,list(phaseNameList))             
                    generalData['Name'] = newName
                    G2frame.G2plotNB.Rename(oldName,generalData['Name'])
                    G2frame.GPXtree.SetItemText(Item,generalData['Name'])
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
                NameTxt.SetValue(generalData['Name'])
                                                
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
                    return      #unchanged - do nothing
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
                                generalData['Super'] = 1
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
                        else:
                            if 'Wave Data' in pages:
                                G2frame.phaseDisplay.DeletePage(pages.index('Wave Data'))
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
                            atom += [[],[],{'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}]
                        wx.CallAfter(UpdateGeneral)
                    else:
                        G2frame.ErrorDialog('Modulation type change error','Can change modulation only if there are no atoms')
                        modulated.SetValue(generalData['Modulated'])                
                
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(General,-1,' Phase name: '),0,WACV)
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
                
            def OnCellChange(event):
                event.Skip()
                SGData = generalData['SGData']
                laue = SGData['SGLaue']
                if laue == '2/m':
                    laue += SGData['SGUniq']
                cell = generalData['Cell']
                Obj = event.GetEventObject()
                ObjId = cellList.index(Obj.GetId())
                try:
                    value = max(1.0,float(Obj.GetValue()))
                except ValueError:
                    if ObjId < 3:               #bad cell edge - reset
                        value = cell[ObjId+1]
                    else:                       #bad angle
                        value = 90.
                if laue in ['m3','m3m']:
                    cell[1] = cell[2] = cell[3] = value
                    cell[4] = cell[5] = cell[6] = 90.0
                    Obj.SetValue("%.5f"%(cell[1]))
                elif laue in ['3R','3mR']:
                    if ObjId == 0:
                        cell[1] = cell[2] = cell[3] = value
                        Obj.SetValue("%.5f"%(cell[1]))
                    else:
                        cell[4] = cell[5] = cell[6] = value
                        Obj.SetValue("%.5f"%(cell[4]))
                elif laue in ['3','3m1','31m','6/m','6/mmm','4/m','4/mmm']:                    
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
                elif laue in ['2/m'+'a']:
                    cell[5] = cell[6] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.SetValue("%.5f"%(cell[ObjId+1]))
                    else:
                        cell[4] = value
                        Obj.SetValue("%.3f"%(cell[4]))
                elif laue in ['2/m'+'b']:
                    cell[4] = cell[6] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.SetValue("%.5f"%(cell[ObjId+1]))
                    else:
                        cell[5] = value
                        Obj.SetValue("%.3f"%(cell[5]))
                elif laue in ['2/m'+'c']:
                    cell[4] = cell[5] = 90.
                    if ObjId != 3:
                        cell[ObjId+1] = value
                        Obj.SetValue("%.5f"%(cell[ObjId+1]))
                    else:
                        cell[6] = value
                        Obj.SetValue("%.3f"%(cell[6]))
                else:
                    cell[ObjId+1] = value
                    if ObjId < 3:
                        Obj.SetValue("%.5f"%(cell[1+ObjId]))
                    else:
                        Obj.SetValue("%.3f"%(cell[1+ObjId]))                        
                cell[7] = G2lat.calc_V(G2lat.cell2A(cell[1:7]))
                volVal.SetValue("%.3f"%(cell[7]))
                density,mattCoeff = G2mth.getDensity(generalData)
                if denSizer:
                    denSizer[1].SetValue('%.3f'%(density))
                    if denSizer[2]:
                        denSizer[2].SetValue('%.3f'%(mattCoeff))
            
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
                if ifEdit:          #a,b,c,etc.
                    cellVal = wx.TextCtrl(General,value=(fmt%(cell[Id+1])),
                        style=wx.TE_PROCESS_ENTER)
                    cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
                    cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
                    cellSizer.Add(cellVal,0,WACV)
                    cellList.append(cellVal.GetId())
                else:               #volume
                    volVal = wx.TextCtrl(General,value=(fmt%(cell[7])),style=wx.TE_READONLY)
                    volVal.SetBackgroundColour(VERY_LIGHT_GREY)
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
                elemSizer.GetChildren()[indx+3*nCols+1].Window.SetValue('%.3f'%(wt))    #tricky
                data['General']['AtomMass'][indx] = wt
                density,mattCoeff = G2mth.getDensity(generalData)
                denSizer[1].SetValue('%.3f'%(density))
                if denSizer[2]:
                    denSizer[2].SetValue('%.3f'%(mattCoeff))
                    
            elemSizer = wx.FlexGridSizer(0,len(generalData['AtomTypes'])+1,1,1)
            elemSizer.Add(wx.StaticText(General,label=' Elements'),0,WACV)
            for elem in generalData['AtomTypes']:
                typTxt = wx.TextCtrl(General,value=elem,style=wx.TE_READONLY)
                typTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(typTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Isotope'),0,WACV)
            for elem in generalData['AtomTypes']:
                choices = list(generalData['Isotopes'][elem].keys())
                isoSel = wx.ComboBox(General,-1,value=generalData['Isotope'][elem],choices=choices,
                    style=wx.CB_READONLY)
                isoSel.Bind(wx.EVT_COMBOBOX,OnIsotope)
                Indx[isoSel.GetId()] = elem
                elemSizer.Add(isoSel,1,WACV|wx.EXPAND)
            elemSizer.Add(wx.StaticText(General,label=' No. per cell'),0,WACV)
            for elem in generalData['AtomTypes']:
                numbTxt = wx.TextCtrl(General,value='%.1f'%(generalData['NoAtoms'][elem]),
                    style=wx.TE_READONLY)
                numbTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(numbTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Atom weight'),0,WACV)
            for wt in generalData['AtomMass']:
                wtTxt = wx.TextCtrl(General,value='%.3f'%(wt),style=wx.TE_READONLY)
                wtTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(wtTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Bond radii'),0,WACV)
            for rad in generalData['BondRadii']:
                bondRadii = wx.TextCtrl(General,value='%.2f'%(rad),style=wx.TE_READONLY)
                bondRadii.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(bondRadii,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Angle radii'),0,WACV)
            for rad in generalData['AngleRadii']:
                elemTxt = wx.TextCtrl(General,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' van der Waals radii'),0,WACV)
            for rad in generalData['vdWRadii']:
                elemTxt = wx.TextCtrl(General,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,WACV)
            elemSizer.Add(wx.StaticText(General,label=' Default color'),0,WACV)
            for R,G,B in generalData['Color']:
                colorTxt = wx.TextCtrl(General,value='',style=wx.TE_READONLY)
                colorTxt.SetBackgroundColour(wx.Colour(R,G,B))
                elemSizer.Add(colorTxt,0,WACV)
            if generalData['Type'] == 'magnetic':
                elemSizer.Add(wx.StaticText(General,label=' Lande g factor: '),0,WACV)
                for ig,elem in enumerate(generalData['AtomTypes']):
                    gfac = generalData['Lande g'][ig]
                    if gfac == None:
                        elemSizer.Add((5,0),)
                    else:
                        gfacTxt = G2G.ValidatedTxtCtrl(General,generalData['Lande g'],ig,
                            min=0.5,max=3.0,nDig=(10,2))
                        elemSizer.Add(gfacTxt,0,WACV)
            return elemSizer
        
        def DenSizer():
            
            generalData['Mass'] = G2mth.getMass(generalData)
            density,mattCoeff = G2mth.getDensity(generalData)
            denSizer = wx.BoxSizer(wx.HORIZONTAL)
            denSizer.Add(wx.StaticText(General,-1,' Density: '),0,WACV)
            denTxt = wx.TextCtrl(General,-1,'%.3f'%(density),style=wx.TE_READONLY)
            denTxt.SetBackgroundColour(VERY_LIGHT_GREY)
            denSizer.Add(denTxt,0,WACV)
            mattTxt = None        
            if generalData['Type'] == 'macromolecular' and generalData['Mass'] > 0.0:
                denSizer.Add(wx.StaticText(General,-1,' Matthews coeff.: '),
                    0,WACV)
                mattTxt = wx.TextCtrl(General,-1,'%.3f'%(mattCoeff),style=wx.TE_READONLY)
                mattTxt.SetBackgroundColour(VERY_LIGHT_GREY)
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
                SGErr,SGData = G2spc.SpcGroup(SpcGrp)
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
            magSizer.Add(wx.StaticText(General,label=' Magnetic spin operator selection:'),0,WACV)
            spinSizer = wx.BoxSizer(wx.HORIZONTAL)
            if SGData['SGFixed']:
                SpnFlp = SGData['SpnFlp']
                spinSizer.Add(wx.StaticText(General,label=' Magnetic phase from mcif file; no change in spin inversion allowed'),0,WACV)
                OprNames = G2spc.GenMagOps(SGData)[0]
            else:
                if not len(GenSym) or SGData['SGGray']:
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
            dminVal = G2G.ValidatedTxtCtrl(General,generalData,'MagDmin',nDig=(10,4),min=0.7)
            dminSizer.Add(dminVal,0,WACV)
            magSizer.Add(dminSizer,0,WACV)
            return magSizer
            
        def ModulatedSizer(name):
            
            def OnSuperGp(event):   #for HKLF needs to reject SSgps not agreeing with modVec!
                event.Skip()
                try:
                    SSymbol = superGp.GetValue()
                except AttributeError:
                    SSymbol = superGp.GetLabel()
                SSGData = generalData['SSGData']
                if not generalData['SGData']['SGFixed']:
                    E,SSGData = G2spc.SSpcGroup(generalData['SGData'],SSymbol)
                if SSGData:
                    Vec = generalData['SuperVec'][0]     #(3+1) only
                    modSymb = SSGData['modSymb']
                    generalData['SuperVec'][0] = G2spc.SSGModCheck(Vec,modSymb)[0]
                    text,table = G2spc.SSGPrint(generalData['SGData'],SSGData)
                    generalData['SSGData'] = SSGData
                    generalData['SuperSg'] = SSymbol
                    msg = 'Superspace Group Information'
                    G2G.SGMessageBox(General,msg,text,table).Show()
                else:
                    text = [E+'\nSuperspace Group set to previous']
                    superGp.SetValue(generalData['SuperSg'])
                    msg = 'Superspace Group Error'
                    Style = wx.ICON_EXCLAMATION
                    Text = '\n'.join(text)
                    wx.MessageBox(Text,caption=msg,style=Style)
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
            SpGrp = SGData['SpGrp']
            if SGData['SGGray']:
                SpGrp += " 1'"
            modSizer.Add(wx.StaticText(General,label=' Superspace group: %s '%SpGrp),0,WACV)
            Choice = []
            if not SGData['SGFixed']:
                Choice = G2spc.SSChoice(SGData)
                if SGData['SGGray']:
                    Choice = [G2spc.fixGray(SGData,item) for item in Choice]
            if len(Choice):
                superGp = wx.ComboBox(General,value=generalData['SuperSg'],choices=Choice,style=wx.CB_DROPDOWN|wx.TE_PROCESS_ENTER)
                superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)
                superGp.Bind(wx.EVT_COMBOBOX,OnSuperGp)
            else:
                superGp = wx.StaticText(General,label=generalData['SuperSg'])
            modSizer.Add(superGp,0,WACV)
            modSizer.Add((5,5),0)
            showOps = wx.Button(General,label='Show ops.')
            showOps.Bind(wx.EVT_BUTTON,OnSuperGp)
            modSizer.Add(showOps,0,WACV)
            if PWDR:
                modSizer.Add(wx.StaticText(General,label=' Max index: '),0,WACV)
                indChoice = ['1','2','3','4','5','6','7']
                Max = wx.ComboBox(General,-1,value='%d'%(generalData['SuperVec'][2]),choices=indChoice,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Max.Bind(wx.EVT_COMBOBOX,OnMax)        
                modSizer.Add(Max,0,WACV)
            ssSizer.Add(modSizer,0,WACV)
            vecSizer = wx.FlexGridSizer(1,5,5,5)
            vecSizer.Add(wx.StaticText(General,label=' Modulation vector: '),0,WACV)
            modS = G2spc.splitSSsym(generalData['SuperSg'])[0]
            generalData['SuperVec'][0],ifShow = G2spc.SSGModCheck(generalData['SuperVec'][0],modS)
            for i,[val,show] in enumerate(zip(generalData['SuperVec'][0],ifShow)):
                if show:
                    modVal = G2G.ValidatedTxtCtrl(General,generalData['SuperVec'][0],i,nDig=(10,4),min=-1.,max=2.)
                    vecSizer.Add(modVal,0,WACV)
                    Indx[modVal.GetId()] = i
                else:
                    modVal = wx.TextCtrl(General,value=('%.3f'%(val)),
                        size=wx.Size(50,20),style=wx.TE_READONLY)
                    modVal.SetBackgroundColour(VERY_LIGHT_GREY)
                    vecSizer.Add(modVal,0,WACV)
            if PWDR:
                Ref = wx.CheckBox(General,label='Refine?')
                Ref.SetValue(generalData['SuperVec'][1])
                Ref.Bind(wx.EVT_CHECKBOX, OnVecRef)
                vecSizer.Add(Ref,0,WACV)
            ssSizer.Add(vecSizer)
            return ssSizer
            
        def PawleySizer():
            
            def OnPawleyRef(event):
                generalData['doPawley'] = pawlRef.GetValue()
            
            pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley controls: '),0,WACV)
            pawlRef = wx.CheckBox(General,-1,label=' Do Pawley refinement?')
            pawlRef.SetValue(generalData['doPawley'])
            pawlRef.Bind(wx.EVT_CHECKBOX,OnPawleyRef)
            pawleySizer.Add(pawlRef,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley dmin: '),0,WACV)
            pawlMin = G2G.ValidatedTxtCtrl(General,generalData,'Pawley dmin',size=(65,25),
                min=0.25,max=20.,nDig=(10,5))
            pawleySizer.Add(pawlMin,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley dmax: '),0,WACV)
            pawlMax = G2G.ValidatedTxtCtrl(General,generalData,'Pawley dmax',size=(65,25),
                min=2.0,max=100.,nDig=(10,5))
            pawleySizer.Add(pawlMax,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley neg. wt.: '),0,WACV)
            pawlNegWt = G2G.ValidatedTxtCtrl(General,generalData,'Pawley neg wt',size=(65,25),
                min=0.,max=1.,nDig=(10,4))
            pawleySizer.Add(pawlNegWt,0,WACV)
            return pawleySizer
            
        def MapSizer():
            
            def OnMapType(event):
                Map['MapType'] = mapType.GetValue()
                
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
                            
            #patch
            if 'cutOff' not in Map:
                Map['cutOff'] = 100.0
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
            mapSizer.Add(lineSizer,0,WACV)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Resolution: '),0,WACV)
            mapRes = G2G.ValidatedTxtCtrl(General,Map,'Resolution',nDig=(10,2),min=0.25,max=20.)
            line2Sizer.Add(mapRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' Peak cutoff %: '),0,WACV)
            cutOff = G2G.ValidatedTxtCtrl(General,Map,'cutOff',nDig=(10,1),min=1.0,max=100.)
            line2Sizer.Add(cutOff,0,WACV)
            mapSizer.Add(line2Sizer,0,WACV)
            return mapSizer
                
        def FlipSizer():
            if 'k-Max' not in Flip: Flip['k-Max'] = 20.
            
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
                    id = int(name.split('hkl')[1])
                    HKL = [int(val) for val in vals]
                    Flip['testHKL'][id] = HKL
                except ValueError:
                    HKL = Flip['testHKL'][id]
                Obj.SetValue('%3d %3d %3d'%(HKL[0],HKL[1],HKL[2]))

            refsList = [item for item in G2gd.GetGPXtreeDataNames(G2frame,['HKLF','PWDR']) if item in data['Histograms'].keys()]
            flipSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Charge flip controls: Reflection sets: '),0,WACV)
            if 'list' not in str(type(Flip['RefList'])):     #patch
                Flip['RefList'] = [Flip['RefList'],]
            lineSizer.Add(wx.ComboBox(General,value=Flip['RefList'][0],choices=Flip['RefList'],
                style=wx.CB_DROPDOWN|wx.CB_READONLY),0,WACV)
            refList = wx.Button(General,label='Select reflection sets')
            refList.Bind(wx.EVT_BUTTON,OnRefList)
            lineSizer.Add(refList,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' Normalizing element: '),0,WACV)
            normElem = wx.Button(General,label=Flip['Norm element'],style=wx.TE_READONLY)
            normElem.Bind(wx.EVT_BUTTON,OnNormElem)
            lineSizer.Add(normElem,0,WACV)
            flipSizer.Add(lineSizer,0,WACV)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Resolution: '),0,WACV)
            flipRes = G2G.ValidatedTxtCtrl(General,Flip,'Resolution',nDig=(10,2),min=0.25,max=2.)
            line2Sizer.Add(flipRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Factor (0.1-1.2): '),0,WACV)
            kFactor = G2G.ValidatedTxtCtrl(General,Flip,'k-factor',nDig=(10,3),min=0.1,max=1.2)
            line2Sizer.Add(kFactor,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Max (>=10.0): '),0,WACV)
            kMax = G2G.ValidatedTxtCtrl(General,Flip,'k-Max',nDig=(10,1),min=10.)
            line2Sizer.Add(kMax,0,WACV)
            flipSizer.Add(line2Sizer,0,WACV)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line3Sizer.Add(wx.StaticText(General,label=' Test HKLs:'),0,WACV)
            if len(Flip['testHKL']) < 5:
                Flip['testHKL'] += [[1,1,1],[0,2,0],[1,2,3]]
            HKL = Flip['testHKL']
            for ih,hkl in enumerate(Flip['testHKL']):                
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
            dmin = G2G.ValidatedTxtCtrl(General,MCSAdata,'dmin',nDig=(10,3),min=1.,max=5.)
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
            ranRange = G2G.ValidatedTxtCtrl(General,MCSAdata,'ranRange',nDig=(10,1),min=1.,max=99.)
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
                    Ajump = G2G.ValidatedTxtCtrl(General,MCSAdata[parms],i,nDig=(10,2),min=0.1,max=1.,OnLeave=ShowTsched)
                    line3Sizer.Add(Ajump,0,WACV)
            elif 'log' in MCSAdata['Algorithm']:
                line3Sizer.Add(wx.StaticText(General,label=' slope: '),0,WACV)
                slope = G2G.ValidatedTxtCtrl(General,MCSAdata,'log slope',nDig=(10,3),min=0.25,max=1.0,OnLeave=ShowTsched)
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

        # UpdateGeneral execution starts here
        phaseTypes = ['nuclear','magnetic','macromolecular','faulted']
        SetupGeneral()
        generalData = data['General']
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
#end patches
        if General.GetSizer():
            General.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(NameSizer(),0)
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

        mainSizer.Add(PawleySizer())
        G2G.HorizontalLine(mainSizer,General)
        
        mainSizer.Add(MapSizer())
        G2G.HorizontalLine(mainSizer,General)
        
        mainSizer.Add(FlipSizer())
        if generalData['Type'] in ['nuclear','macromolecular','faulted',]:
            G2G.HorizontalLine(mainSizer,General)
            mainSizer.Add(MCSASizer())
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
        while True:
            dlg = TransformDialog(G2frame,data,Trans,Uvec,Vvec,ifMag,BNSlatt)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    newPhase,Trans,Uvec,Vvec,ifMag,ifConstr,Common = dlg.GetSelection()
                    newPhase['ranId'] = ran.randint(0,sys.maxsize)
                    SGData = newPhase['General']['SGData']
                    if ifMag:
                        BNSlatt = SGData['BNSlattsym'][0]
                        
                    SGData['GenSym'],SGData['GenFlg'],BNSsym = G2spc.GetGenSym(SGData)
                    if '_' in BNSlatt:
                        SGData['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
                        G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
                    SGData['SpnFlp'] = G2spc.GenMagOps(SGData)[1]
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
                data['Drawing'] = []
                break
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
            UseList[hist]['Mustrain'][4:6] = [NShkl*[0.01,],NShkl*[False,]]
            UseList[hist]['HStrain'] = [NDij*[0.0,],NDij*[False,]]
        newPhase['General']['Map'] = mapDefault.copy()
        sub = G2frame.GPXtree.AppendItem(parent=
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases'),text=phaseName)
        G2frame.GPXtree.SetItemPyData(sub,newPhase)
        newPhase['Drawing'] = []
        if ifConstr:
            G2cnstG.TransConstraints(G2frame,data,newPhase,Trans,Vvec,atCodes)     #data is old phase
        G2frame.GPXtree.SelectItem(sub)
                
    def OnUseBilbao(event):
        PatternName = data['magPhases']
        PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,PatternName)
        UnitCellsId = G2gd.GetGPXtreeItemId(G2frame,PatternId, 'Unit Cells List')
        magData = G2frame.GPXtree.GetItemPyData(UnitCellsId)[5]
        magKeep = []
        magIds = []
        magchoices = []
        for mid,magdata in enumerate(magData):
            if magdata['Keep']:
                magdata['No.'] = mid+1
                trans = G2spc.Trans2Text(magdata['Trans'])
                vec = G2spc.Latt2text([magdata['Uvec'],])
                magKeep.append(magdata)
                magIds.append(mid)
                magchoices.append('(%d) %s; (%s) + (%s)'%(mid+1,magdata['Name'],trans,vec))
        if not len(magKeep):
            G2frame.ErrorDialog('Magnetic phase selection error','No magnetic phases found; be sure to "Keep" some')
            return
        dlg = wx.SingleChoiceDialog(G2frame,'Select magnetic space group','Make new magnetic phase',magchoices)
        opt = dlg.ShowModal()
        if opt == wx.ID_OK:
            sel = dlg.GetSelection()
            magchoice = magKeep[sel]
            magId = magIds[sel]
            phaseName = '%s mag_%d'%(data['General']['Name'],magchoice['No.'])
            newPhase = copy.deepcopy(data)
            newPhase['ranId'] = ran.randint(0,sys.maxsize),
            del newPhase['magPhases']
            generalData = newPhase['General']
            generalData['Name'] = phaseName
            generalData['SGData'] = copy.deepcopy(magchoice['SGData'])            
            generalData['Cell'][1:] = magchoice['Cell'][:]
            SGData = generalData['SGData']
            vvec = np.array([0.,0.,0.])
            newPhase,atCodes = G2lat.TransformPhase(data,newPhase,magchoice['Trans'],magchoice['Uvec'],vvec,True)
            Atoms = newPhase['Atoms']
            Atms = []
            AtCods = []
            atMxyz = []
            for ia,atom in enumerate(Atoms):
                if not len(G2elem.GetMFtable([atom[1],],[2.0,])):
                    continue
                atom[0] += '_%d'%ia
                SytSym,Mul,Nop,dupDir = G2spc.SytSym(atom[3:6],SGData)
                CSI = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
                Atms.append(atom)
                AtCods.append(atCodes[ia])
                MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                atMxyz.append([MagSytSym,CSI[0]])
            dlg = UseMagAtomDialog(G2frame,magchoices[sel],Atms,AtCods,atMxyz,ifDelete=True)
            try:
                opt = dlg.ShowModal()
                if  opt == wx.ID_YES:
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
        G2cnstG.TransConstraints(G2frame,data,newPhase,magchoice['Trans'],vvec,atCodes)     #data is old phase
        G2frame.newGPXfile = phaseName+'.gpx'
        G2frame.OnFileSaveas(event)
        G2frame.GPXtree.SelectItem(sub)
        
################################################################################
#####  Atom routines
################################################################################

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
                    pass                                        #& then change all 'I' atoms
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
                            DrawAtomsReplaceByID(data['Drawing'],ui+6,atomData[r],ID)
                    wx.CallAfter(Paint)
                    
        def ChangeAtomCell(event):

            def chkUij(Uij,CSI): #needs to do something!!!
                return Uij

            r,c =  event.GetRow(),event.GetCol()
            if r >= 0 and c >= 0:
                ci = colLabels.index('I/A')
                ID = atomData[r][ci+8]
                SGData = generalData['SGData']
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
                        atomData[r][c+1] =  0.0
                        Atoms.SetCellStyle(r,c+1,VERY_LIGHT_GREY,True)
                        Atoms.SetCellTextColour(r,c+1,VERY_LIGHT_GREY)
                        for i in range(6):
                            ci = i+colLabels.index('U11')
                            atomData[r][ci] = value*CSI[3][i]
                            Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                            Atoms.SetCellTextColour(r,ci,BLACK)
                            if CSI[2][i]:
                                Atoms.SetCellStyle(r,ci,WHITE,False)
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
                if 'Atoms' in data['Drawing']:
                    ci = colLabels.index('I/A')
                    DrawAtomsReplaceByID(data['Drawing'],ci+8,atomData[r],ID)
                    G2plt.PlotStructure(G2frame,data)
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
                    DrawAtomsReplaceByID(data['Drawing'],ci+8,atomData[r],ID)
                    G2plt.PlotStructure(G2frame,data)
                SetupGeneral()
            else:
                event.Skip()

        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if not (event.AltDown() or (event.ShiftDown() and event.ControlDown())):
                Atoms.frm = -1
                G2frame.GetStatusBar().SetStatusText('',1)                    
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
                            UpdateDrawAtoms()
                            wx.CallAfter(Paint)
                else:
                    Atoms.ClearSelection()
                    Atoms.SelectRow(r,True)
                
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
                    
        def Paint():
            
            table = []
            rowLabels = []
            for i,atom in enumerate(atomData):
                table.append(atom)
                rowLabels.append(str(i))
            atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            try:
                Atoms.SetTable(atomTable, True)    # Paint may be called after the Grid has been deleted
            except:
                return
            Atoms.frm = -1            
            colType = colLabels.index('Type')
            colR = colLabels.index('refine')
            colSS = colLabels.index('site sym')
            colX = colLabels.index('x')
            colIA = colLabels.index('I/A')
            colU11 = colLabels.index('U11')
            colUiso = colLabels.index('Uiso')
            colM = 0
            if 'Mx' in colLabels:
                colM = colLabels.index('Mx')
#                atTypes = generalData['AtomTypes']
#                Lande = generalData['Lande g']
#                AtInfo = dict(zip(atTypes,Lande))
            attr = wx.grid.GridCellAttr()
            attr.IncRef()               #fix from Jim Hester
            attr.SetEditor(G2G.GridFractionEditor(Atoms))
            for c in range(colX,colX+3):
                attr = wx.grid.GridCellAttr()
                attr.IncRef()               #fix from Jim Hester
                attr.SetEditor(G2G.GridFractionEditor(Atoms))
                Atoms.SetColAttr(c, attr)
            for i in range(colU11-1,colU11+6):
                Atoms.SetColSize(i,50)            
            for row in range(Atoms.GetNumberRows()):
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
                        Atoms.SetCellStyle(row,cj,VERY_LIGHT_GREY,True)
                        if CSI[2][i] and 'U' not in rbExcl:
                            Atoms.SetCellStyle(row,cj,WHITE,False)
                else:
                    Atoms.SetCellStyle(row,colUiso,WHITE,False)
                    Atoms.SetCellTextColour(row,colUiso,BLACK)
                    if 'U' in rbExcl:
                        Atoms.SetCellStyle(row,colUiso,VERY_LIGHT_GREY,True)
                    for i in range(6):
                        cj = colU11+i
                        Atoms.SetCellStyle(row,cj,VERY_LIGHT_GREY,True)
                        Atoms.SetCellTextColour(row,cj,VERY_LIGHT_GREY)
                if colM:
                    SytSym,Mul,Nop,dupDir = G2spc.SytSym(atomData[row][colX:colX+3],SGData)
                    MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                    Atoms.SetCellValue(row,colSS,MagSytSym)
                    CSI = []
                    if not SGData['SGGray']:
                        CSI = G2spc.GetCSpqinel(SpnFlp,dupDir)
#                    print (SytSym,Nop,SpnFlp[Nop],CSI,dupDir)
#                    print('CSI:',CSI)
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
                if 'X' in rbExcl:
                    for c in range(0,colX+3):
                        if c != colR:
                            Atoms.SetCellStyle(row,c,VERY_LIGHT_GREY,True)
                Atoms.SetReadOnly(row,colType,True)
                Atoms.SetReadOnly(row,colSS,True)                         #site sym
                Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            Atoms.AutoSizeColumns(False)
            SetPhaseWindow(Atoms)

        # FillAtomsGrid executable code starts here
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        generalData = data['General']
        SpnFlp = generalData['SGData'].get('SpnFlp',[])
#        OprNames = generalData['SGData'].get('OprNames',[])
#        print OprNames
#        print SpnFlp
#        print generalData['SGData'].get('MagMom',[])
        atomData = data['Atoms']
        resRBData = data['RBModels'].get('Residue',[])
        vecRBData = data['RBModels'].get('Vector',[])
        rbAtmDict = {}
        for rbObj in resRBData+vecRBData:
            exclList = ['X' for i in range(len(rbObj['Ids']))]
            rbAtmDict.update(dict(zip(rbObj['Ids'],exclList)))
            if rbObj['ThermalMotion'][0] != 'None':
                for id in rbObj['Ids']:
                    rbAtmDict[id] += 'U'            
        # exclList will be 'x' or 'xu' if TLS used in RB
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
        if SGData['SGPolax']:
            G2frame.GetStatusBar().SetStatusText('Warning: The location of the origin is arbitrary in '+SGData['SGPolax'],1)
        if 'phoenix' in wx.version():
            Atoms.Bind(wg.EVT_GRID_CELL_CHANGED, ChangeAtomCell)
        else:
            Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeAtomCell)
        Atoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, ChangeSelection)
        Atoms.SetMargins(0,0)
        
        Paint()

    def OnAtomAdd(event):
        Elem = 'H'
        if data['General']['Type'] == 'magnetic':
            Elem = 'Fe'
        AtomAdd(0,0,0,El=Elem)
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
            AtomAdd(0,0,0,El=Elem)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)
                
    def AtomAdd(x,y,z,El='H',Name='UNK'):
        atomData = data['Atoms']
        generalData = data['General']
        atId = ran.randint(0,sys.maxsize)
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
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)
#        if 'Atoms' in data['Drawing']:            
#            DrawAtomAdd(data['Drawing'],atomData[-1])

    def OnAtomInsert(event):
        '''Inserts a new atom into list immediately before every selected atom
        '''
        indx = GetSelectedAtoms()
        for a in reversed(sorted(indx)):
            AtomInsert(a,0,0,0)
        event.StopPropagation()
        FillAtomsGrid(Atoms)
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms()
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
        indx = GetSelectedAtoms()
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
        cx,ct,cs,cia = generalData['AtomPtrs']
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
            UpdateDrawAtoms()
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
        UpdateDrawAtoms()
        FillAtomsGrid(Atoms)
        G2plt.PlotStructure(G2frame,data)
        
    def OnAtomMove(event):
        drawData = data['Drawing']
        atomData = data['Atoms']
        x,y,z = drawData['viewPoint'][0]
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        cx = colLabels.index('x')
        ci = colLabels.index('I/A')
        indx = GetSelectedAtoms()
        if len(indx) != 1:
            G2frame.ErrorDialog('Atom move error','Only one atom can be moved')
        elif atomData[indx[0]][ci+8] in rbAtmDict:
            G2frame.ErrorDialog('Atom move error','Atoms in rigid bodies can not be moved')
        else:
            atomData[indx[0]][cx:cx+3] = [x,y,z]
            SetupGeneral()
            FillAtomsGrid(Atoms)
            ID = atomData[indx[0]][ci+8]
            DrawAtomsReplaceByID(data['Drawing'],ci+8,atomData[indx[0]],ID)
            G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()
            
    def DrawAtomsReplaceByID(drawingData,loc,atom,ID):
        IDs = [ID,]
        atomData = drawingData['Atoms']
        indx = G2mth.FindAtomIndexByIDs(atomData,loc,IDs)
        for ind in indx:
            atomData[ind] = MakeDrawAtom(atom,atomData[ind])
                
    def MakeDrawAtom(atom,oldatom=None):
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        SGData = generalData['SGData']
        if generalData['Type'] in ['nuclear','faulted',]:
            if oldatom:
                opr = oldatom[5]
                if atom[9] == 'A':                    
                    X,U = G2spc.ApplyStringOps(opr,SGData,atom[3:6],atom[11:17])
                    atomInfo = [atom[:2]+list(X)+oldatom[5:9]+atom[9:11]+list(U)+oldatom[17:]][0]
                else:
                    X = G2spc.ApplyStringOps(opr,SGData,atom[3:6])
                    atomInfo = [atom[:2]+list(X)+oldatom[5:9]+atom[9:]+oldatom[17:]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+['1',]+['vdW balls',]+
                    ['',]+[[255,255,255],]+atom[9:]+[[],[]]][0]
            ct,cs = [1,8]         #type & color
        elif  generalData['Type'] == 'magnetic':
            if oldatom:
                opr = oldatom[8]
                mom = np.array(atom[7:10])
                Mom = G2spc.ApplyStringOpsMom(opr,SGData,mom)
                if atom[12] == 'A':                    
                    X,U = G2spc.ApplyStringOps(opr,SGData,atom[3:6],atom[14:20])
                    atomInfo = [atom[:2]+list(X)+list(Mom)+oldatom[8:12]+atom[12:14]+list(U)+oldatom[20:]][0]
                else:
                    X = G2spc.ApplyStringOps(opr,SGData,atom[3:6])
                    atomInfo = [atom[:2]+list(X)+list(Mom)+oldatom[8:12]+atom[12:]+oldatom[20:]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+atom[7:10]+['1',]+['vdW balls',]+
                    ['',]+[[255,255,255],]+atom[12:]+[[],[]]][0]
            ct,cs = [1,11]         #type & color
        elif generalData['Type'] == 'macromolecular':
            try:
                oneLetter = AA3letter.index(atom[1])
            except ValueError:
                oneLetter = -1
            atomInfo = [[atom[1].strip()+atom[0],]+
                [AA1letter[oneLetter]+atom[0],]+atom[2:5]+
                atom[6:9]+['1',]+['sticks',]+['',]+[[255,255,255],]+atom[12:]+[[],[]]][0]
            ct,cs = [4,11]         #type & color
        atNum = generalData['AtomTypes'].index(atom[ct])
        atomInfo[cs] = list(generalData['Color'][atNum])
        return atomInfo
        
    def AtomInsert(indx,x,y,z,El='H',Name='UNK'):
        atomData = data['Atoms']
        generalData = data['General']
        SGData = generalData['SGData']
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)[:2]
        atId = ran.randint(0,sys.maxsize)
        if generalData['Type'] == 'macromolecular':
            atomData.insert(indx,[0,Name,'',Name,El,'',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] in ['nuclear','faulted',]:
            if generalData['Modulated']:
                atomData.insert(indx,[Name,El,'',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.insert(indx,[Name,El,'',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
            SetupGeneral()
        elif generalData['Type'] == 'magnetic':
            if generalData['Modulated']:
                atomData.insert(indx,[Name,El,'',x,y,z,1,0.,0.,0.,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.insert(indx,[Name,El,'',x,y,z,1,0.,0.,0.,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)

    def AtomDelete(event):
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        HydIds = data['General']['HydIds']
        ci = colLabels.index('I/A')
        indx = GetSelectedAtoms()
        IDs = []
        if not indx: return
        atomData = data['Atoms']
        indx.sort()
        indx.reverse()
        for ind in indx:
            atom = atomData[ind]
            if atom[ci+8] in rbAtmDict:
                G2frame.GetStatusBar().SetStatusText('**** ERROR - atom is in a rigid body and can not be deleted ****',1)
            else:
                if atom[ci+8] in HydIds:    #remove Hs from Hatom update dict
                    del HydIds[atom[ci+8]]
                IDs.append(atom[ci+8])
                del atomData[ind]
        if 'Atoms' in data['Drawing']:
            Atoms.ClearSelection()
            DrawAtomsDeleteByIDs(IDs)
            data['Drawing']['Atoms'] = []
            UpdateDrawAtoms()
            wx.CallAfter(FillAtomsGrid,Atoms)
            G2plt.PlotStructure(G2frame,data)
        SetupGeneral()
        if not len(HydIds):
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,False)
        event.StopPropagation()

    def GetSelectedAtoms():
        '''Get all atoms that are selected by row or by having any cell selected.
        produce an error message if no atoms are selected.
        '''
        indx = list(set([row for row,col in Atoms.GetSelectedCells()]+Atoms.GetSelectedRows()))
        if indx:
            return indx
        else:
            G2G.G2MessageBox(G2frame,'Warning: no atoms were selected','Nothing selected')
        
    def AtomRefine(event):
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        c = colLabels.index('refine')
        indx = GetSelectedAtoms()
        if not indx: return
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
        indx = GetSelectedAtoms()
        if not indx: return
        atomData = data['Atoms']
        generalData = data['General']
        ifMag = False
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        ci = colLabels.index('I/A')
        choices = ['Type','Name','x','y','z','frac','I/A','Uiso']
        if generalData['Type'] == 'magnetic':
            choices += ['Mx','My','Mz',]
            ifMag = True
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom parameter',choices)
        parm = ''
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            parm = choices[sel]
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
                            DrawAtomsReplaceByID(data['Drawing'],ci+8,atomData[r],ID)
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
                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][cid] = parm
                FillAtomsGrid(Atoms)
            dlg.Destroy()
        elif parm in ['frac','Uiso']:
            limits = [0.,1.]
            val = 1.0
            if  parm in ['Uiso']:
                limits = [0.,0.25]
                val = 0.01
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for '+parm,val,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in indx:                        
                    if not Atoms.IsReadOnly(r,0):   #not if in RB!
                        atomData[r][cid] = parm
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
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)

    def AtomTransform(event):
        indx = GetSelectedAtoms()
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
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)
            
#    def AtomRotate(event):
#        '''Currently not used - Bind commented out below
#        '''
#        Units = {'':np.zeros(3),
#            'xy':np.array([[i,j,0] for i in range(3) for j in range(3)])-np.array([1,1,0]),
#            'xz':np.array([[i,0,j] for i in range(3) for j in range(3)])-np.array([1,1,0]),
#            'yz':np.array([[0,i,j] for i in range(3) for j in range(3)])-np.array([1,1,0]),
#            'xyz':np.array([[i,j,k] for i in range(3) for j in range(3) for k in range(3)])-np.array([1,1,1])}
#        indx = GetSelectedAtoms()
#        if indx:
#            generalData = data['General']
#            A,B = G2lat.cell2AB(generalData['Cell'][1:7])
#            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
#            cx = colLabels.index('x')
#            cuia = colLabels.index('I/A')   #need to not do aniso atoms - stop with error? or force isotropic?
#            css = colLabels.index('site sym')
#            atomData = data['Atoms']
#            SGData = generalData['SGData']
#            dlg = G2gd.RotationDialog(G2frame)
#            try:
#                if dlg.ShowModal() == wx.ID_OK:
#                    M,T,Expand = dlg.GetSelection()
#                    Unit = Units[Expand]
#                    for ind in indx:
#                        XYZ = np.array(atomData[ind][cx:cx+3])
#                        XYZS = XYZ+Unit
#                        XYZS -= T
#                        XYZS = np.inner(A,XYZS).T   #to Cartesian
#                        XYZS = np.inner(M,XYZS).T   #rotate
#                        XYZS = np.inner(B,XYZS).T+T #back to crystal & translate
#                        GSASIIpath.IPyBreak()
#                        atomData[ind][cx:cx+3] = XYZ
#                        for unit in Unit:
#                            XYZ = np.copy(np.array(atomData[ind][cx:cx+3]))
#                            XYZ += unit 
#                            XYZ -= T
#                            XYZ = np.inner(A,XYZ)   #to Cartesian
#                            XYZ = np.inner(M,XYZ)   #rotate
#                            XYZ = np.inner(B,XYZ)+T #back to crystal & translate
#                            if np.all(XYZ>=0.) and np.all(XYZ<1.0):
#                                atomData[ind][cx:cx+3] = XYZ
##                                atom[css:css+2] = G2spc.SytSym(XYZ,SGData)[:2]
#                                break
#            finally:
#                dlg.Destroy()
#            Atoms.ClearSelection()
#            Atoms.ForceRefresh()
#        else:
#            print "select one or more rows of atoms"
#            G2frame.ErrorDialog('Select atom',"select one or more atoms then redo")
                
    def MakeMolecule(event):      
        indx = GetSelectedAtoms()
        DisAglCtls = {}
        if len(indx) == 1:
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
        msg = 'Density of phase {:s} = {:.3f} g/cc'.format(data['General']['Name'],density)
        print(msg)
        G2G.G2MessageBox(G2frame,msg,'Density')
        
    def OnSetAll(event):
        'set refinement flags for all atoms in table'
        for row in range(Atoms.GetNumberRows()):
            Atoms.SelectRow(row,True)
    
    def OnDistAnglePrt(event):
        'save distances and angles to a file'    
        fp = open(os.path.abspath(os.path.splitext(G2frame.GSASprojectfile)[0]+'.disagl'),'w')
        OnDistAngle(event,fp=fp)
        fp.close()
    
    def OnDistAngle(event,fp=None):
        'Compute distances and angles'    
        indx = GetSelectedAtoms()
        Oxyz = []
        xyz = []
        DisAglData = {}
        DisAglCtls = {}
        if indx:
            generalData = data['General']
            DisAglData['OrigIndx'] = indx
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
            atomData = data['Atoms']
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            cx = colLabels.index('x')
            cn = colLabels.index('Name')
            for i,atom in enumerate(atomData):
                xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
                if i in indx:
                    Oxyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
            DisAglData['OrigAtoms'] = Oxyz
            DisAglData['TargAtoms'] = xyz
            generalData = data['General']
            DisAglData['SGData'] = generalData['SGData']
            DisAglData['Cell'] = generalData['Cell'][1:] #+ volume
            if 'pId' in data:
                DisAglData['pId'] = data['pId']
                DisAglData['covData'] = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Covariance'))
            try:
                if fp:
                    G2stMn.PrintDistAngle(DisAglCtls,DisAglData,fp)
                else:    
                    G2stMn.PrintDistAngle(DisAglCtls,DisAglData)
            except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
                G2frame.ErrorDialog('Distance/Angle calculation','try again but do "Reset" to fill in missing atom types')
        else:
            print ("select one or more rows of atoms")
            G2frame.ErrorDialog('Select atom',"select one or more rows of atoms then redo")
                        
    def OnValidProtein(event):
        
        def pickHandler(resName):
            drawData = data['Drawing']
            resid = resIDs[resName]
            drawData['viewPoint'][0] = atomData[AtLookUp[resid]][cx:cx+3]
            UpdateDrawAtoms()
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

    def OnIsoDistortCalc(event):
        '''Compute the ISODISTORT mode values from the current coordinates.
        Called in response to the (Phase/Atoms tab) AtomCompute
        "Compute ISODISTORT mode values" menu item, which should be enabled
        only when Phase['ISODISTORT'] is defined. 
        '''
        def _onClose(event):
            dlg.EndModal(wx.ID_CANCEL)
        def fmtHelp(item,fullname):
            helptext = "A new variable"
            if item[-3]:
                helptext += " named "+str(item[-3])
            helptext += " is a linear combination of the following parameters:\n"
            first = True
            for term in item[:-3]:
                line = ''
                var = str(term[1])
                m = term[0]
                if first:
                    first = False
                    line += ' = '
                else:
                    if m >= 0:
                        line += ' + '
                    else:
                        line += ' - '
                    m = abs(m)
                line += '%.3f*%s '%(m,var)
                varMean = G2obj.fmtVarDescr(var)
                helptext += "\n" + line + " ("+ varMean + ")"
            helptext += '\n\nISODISTORT full name: '+str(fullname)
            return helptext

        if 'ISODISTORT' not in data:
            raise Exception("Should not happen: 'ISODISTORT' not in data")
        if len(data.get('Histograms',[])) == 0:
            G2frame.ErrorDialog(
                'No data',
                'Sorry, this computation requires that a histogram first be added to the phase'
                )
            return
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree() # init for constraint
        # make a lookup table for constraints
        sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints') 
        Constraints = G2frame.GPXtree.GetItemPyData(sub)
        constDict = {}
        for item in Constraints:
            if item.startswith('_'): continue
            for c in Constraints[item]:
                if c[-1] != 'f' or not c[-3]: continue
                constDict[c[-3]] = c

        ISO = data['ISODISTORT']
        parmDict,varyList = G2frame.MakeLSParmDict()
            
        dlg = wx.Dialog(G2frame,wx.ID_ANY,'ISODISTORT mode values',#size=(630,400),
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                                    'ISODISTORT mode computation for cordinates in phase '+
                                    str(data['General'].get('Name'))))
        aSizer = wx.BoxSizer(wx.HORIZONTAL)
        panel1 = wxscroll.ScrolledPanel(
            dlg, wx.ID_ANY,#size=(100,200),
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        subSizer1 = wx.FlexGridSizer(cols=2,hgap=5,vgap=2)
        panel2 = wxscroll.ScrolledPanel(
            dlg, wx.ID_ANY,#size=(100,200),
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        subSizer2 = wx.FlexGridSizer(cols=3,hgap=5,vgap=2)
        subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,'Parameter name  '))
        subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,' value'),0,wx.ALIGN_RIGHT)
        subSizer2.Add((-1,-1))
        subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'Mode name  '))
        subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,' value'),0,wx.ALIGN_RIGHT)
        
        if 'G2VarList' in ISO:
            deltaList = []
            for gv,Ilbl in zip(ISO['G2VarList'],ISO['IsoVarList']):
                dvar = gv.varname()
                var = dvar.replace('::dA','::A')
                albl = Ilbl[:Ilbl.rfind('_')]
                v = Ilbl[Ilbl.rfind('_')+1:]
                pval = ISO['ParentStructure'][albl][['dx','dy','dz'].index(v)]
                if var in parmDict:
                    cval = parmDict[var][0]
                else:
                    dlg.EndModal(wx.ID_CANCEL)
                    G2frame.ErrorDialog('Atom not found',"No value found for parameter "+str(var))
                    return
                deltaList.append(cval-pval)
            modeVals = np.inner(ISO['Var2ModeMatrix'],deltaList)
            for lbl,xyz,var,val,G2var in zip(ISO['IsoVarList'],deltaList,
                                             ISO['IsoModeList'],modeVals,ISO['G2ModeList']):
                if G2var in constDict:
                    ch = G2G.HelpButton(panel2,fmtHelp(constDict[G2var],var))
                    subSizer2.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
                else:
                    subSizer2.Add((-1,-1))
                subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,str(lbl)))
                try:
                    value = G2py3.FormatSigFigs(xyz)
                except TypeError:
                    value = str(xyz)            
                subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)
                subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,str(var)))
                try:
                    value = G2py3.FormatSigFigs(val)
                except TypeError:
                    value = str(val)            
                subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)
        if 'G2OccVarList' in ISO:
            deltaList = []
            for gv,Ilbl in zip(ISO['G2OccVarList'],ISO['OccVarList']):
                var = gv.varname()
                albl = Ilbl[:Ilbl.rfind('_')]
                #v = Ilbl[Ilbl.rfind('_')+1:]
                pval = ISO['BaseOcc'][albl]
                if var in parmDict:
                    cval = parmDict[var][0]
                else:
                    dlg.EndModal(wx.ID_CANCEL)
                    G2frame.ErrorDialog('Atom not found',"No value found for parameter "+str(var))
                    return
                deltaList.append(cval-pval)
            modeVals = np.inner(ISO['Var2OccMatrix'],deltaList)
            for lbl,xyz,var,val,G2var in zip(ISO['OccVarList'],deltaList,
                                             ISO['OccModeList'],modeVals,ISO['G2OccModeList']):
                if G2var in constDict:
                    ch = G2G.HelpButton(panel2,fmtHelp(constDict[G2var],var))
                    subSizer2.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
                else:
                    subSizer2.Add((-1,-1))
                subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,str(lbl)))
                try:
                    value = G2py3.FormatSigFigs(xyz)
                except TypeError:
                    value = str(xyz)            
                subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)
                #subSizer.Add((10,-1))
                subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,str(var)))
                try:
                    value = G2py3.FormatSigFigs(val)
                except TypeError:
                    value = str(val)            
                subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)

        # finish up ScrolledPanel
        panel1.SetSizer(subSizer1)
        panel2.SetSizer(subSizer2)
        panel1.SetAutoLayout(1)
        panel1.SetupScrolling()
        panel2.SetAutoLayout(1)
        panel2.SetupScrolling()
        # Allow window to be enlarged but not made smaller
        dlg.SetSizer(mainSizer)
        w1,l1 = subSizer1.GetSize()
        w2,l2 = subSizer2.GetSize()
        panel1.SetMinSize((w1+10,200))
        panel2.SetMinSize((w2+20,200))
        aSizer.Add(panel1,1, wx.ALL|wx.EXPAND,1)
        aSizer.Add(panel2,2, wx.ALL|wx.EXPAND,1)
        mainSizer.Add(aSizer,1, wx.ALL|wx.EXPAND,1)

        # make OK button 
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(dlg, wx.ID_CLOSE) 
        btn.Bind(wx.EVT_BUTTON,_onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)

        mainSizer.Fit(dlg)
        dlg.SetMinSize(dlg.GetSize())
        dlg.ShowModal()
        dlg.Destroy()
        
    def OnReImport(event):
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        reqrdr = G2frame.dataWindow.ReImportMenuId.get(event.GetId())
        rdlist = G2frame.OnImportGeneric(reqrdr,
            G2frame.ImportPhaseReaderlist,'phase')
        if len(rdlist) == 0: return
        # rdlist is only expected to have one element
        rd = rdlist[0]
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
                idx = atomNames.index(''.join(atom[:ct+1]).capitalize())  #eliminate spurious differences
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
        wx.CallAfter(FillAtomsGrid,Atoms)
        
################################################################################
#### Layer Data page
################################################################################
        
    def UpdateLayerData(Scroll=0):
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
            wx.CallAfter(UpdateLayerData)
        
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
                id = Indx[event.GetEventObject()]
                Layers['Width'][1][id] = not Layers['Width'][1][id]
            
            Labels = ['a','b']
            flags = Layers['Width'][1]
            widthSizer = wx.BoxSizer(wx.HORIZONTAL)
            for i in range(2):
                widthSizer.Add(wx.StaticText(layerData,label=u' layer width(%s) (<= 1\xb5m): '%(Labels[i])),0,WACV)
                widthVal = G2G.ValidatedTxtCtrl(layerData,Layers['Width'][0],i,nDig=(10,3),min=0.005,max=1.0)
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
            wx.CallAfter(UpdateLayerData)
            
        def OnDeleteLast(event):
            del(data['Layers']['Layers'][-1])
            del(data['Layers']['Transitions'][-1])
            for trans in data['Layers']['Transitions']:
                del trans[-1]
            wx.CallAfter(UpdateLayerData)
                
        def OnImportLayer(event):
            dlg = wx.FileDialog(G2frame, 'Choose GSAS-II project file', 
                wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN| wx.FD_CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    GPXFile = dlg.GetPath()
                    phaseNames = G2strIO.GetPhaseNames(GPXFile)
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
            Phase = G2strIO.GetAllPhaseData(GPXFile,PhaseName)
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
            wx.CallAfter(UpdateLayerData)
            
        def LayerSizer(il,Layer):
            
            def OnNameChange(event):
                event.Skip()
                Layer['Name'] = layerName.GetValue()                
                wx.CallAfter(UpdateLayerData)
                
            def OnAddAtom(event):
                Layer['Atoms'].append(['Unk','Unk',0.,0.,0.,1.,0.01])
                wx.CallAfter(UpdateLayerData)
                
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
                    wx.CallAfter(UpdateLayerData)
                else:
                    event.Skip()
                    
            def OnDrawLayer(event):
                drawLayer.SetValue(False)
                G2plt.PlotLayers(G2frame,Layers,[il,],plotDefaults)
                
            def OnSameAs(event):
                Layer['SameAs'] = sameas.GetValue()
                wx.CallAfter(UpdateLayerData)
                    
            layerSizer = wx.BoxSizer(wx.VERTICAL)
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)            
            nameSizer.Add(wx.StaticText(layerData,label=' Layer name: '),0,WACV)
            layerName = wx.TextCtrl(layerData,value=Layer['Name'],style=wx.TE_PROCESS_ENTER)
            layerName.Bind(wx.EVT_TEXT_ENTER,OnNameChange)        
            layerName.Bind(wx.EVT_KILL_FOCUS,OnNameChange)
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
            atomGrid.SetTable(atomTable,True)
#            atomGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
            for c in range(2,5):
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
                    G2plt.PlotLayers(G2frame,Layers,[Yi,Xi,],plotDefaults)
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
                wx.CallAfter(UpdateLayerData)
                
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
            transSizer.Add(wx.StaticText(layerData,label=' Layer-Layer transition probabilities: '),0,WACV)
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            normprob = wx.CheckBox(layerData,label=' Normalize probabilities?')
            normprob.Bind(wx.EVT_CHECKBOX,OnNormProb)
            topSizer.Add(normprob,0,WACV)
            symprob = wx.CheckBox(layerData,label=' Symmetric probabilities?')
            symprob.SetValue(Layers.get('SymTrans',False))
            symprob.Bind(wx.EVT_CHECKBOX,OnSymProb)
            topSizer.Add(symprob,0,WACV)
            transSizer.Add(topSizer,0,WACV)
            Names = [layer['Name'] for layer in Layers['Layers']]
            transArray = Layers['Transitions']
            layerData.transGrids = []
            if not Names or not transArray:
                return transSizer
            diagSum = 0.
            for Yi,Yname in enumerate(Names):
                transSizer.Add(wx.StaticText(layerData,label=' From %s to:'%(Yname)),0,WACV)
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
                transGrid.SetTable(transTable,True)
#                transGrid.SetScrollRate(0,0)    #get rid of automatic scroll bars
                Indx[transGrid.GetId()] = Yi
                for c in range(0,4):
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
                transSizer.Add(totalFault,0,WACV)
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
                G2plt.PlotLayers(G2frame,Layers,vals,plotDefaults)
            
            Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
            plotSizer = wx.BoxSizer(wx.VERTICAL)
            Str = ' Using sequence nos. from:'
            for name in Names:
                Str += name
            plotSizer.Add(wx.StaticText(layerData,label=Str[:-1]),0,WACV)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(layerData,label=' Enter sequence of layers to plot:'),0,WACV)
            plotSeq = wx.TextCtrl(layerData,value = '',style=wx.TE_PROCESS_ENTER)
            plotSeq.Bind(wx.EVT_TEXT_ENTER,OnPlotSeq)        
            plotSeq.Bind(wx.EVT_KILL_FOCUS,OnPlotSeq)
            lineSizer.Add(plotSeq,0,WACV)
            plotSizer.Add(lineSizer,0,WACV)
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
                wx.CallAfter(UpdateLayerData)
                
            def OnSeqType(event):
                newType = seqType.GetValue()
                if newType == data['Layers']['Stacking'][1]:
                    return
                data['Layers']['Stacking'][1] = newType
                if newType == 'random':
                    data['Layers']['Stacking'][2] = '250'
                else: #List
                    data['Layers']['Stacking'][2] = ''
                wx.CallAfter(UpdateLayerData)
                
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
            stackSizer.Add(wx.StaticText(layerData,label=' Layer stacking parameters:'),0,WACV)
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
                    stackSizer.Add(topLine,0,WACV)
                    Names = [' %s: %d,'%(layer['Name'],iL+1) for iL,layer in enumerate(Layers['Layers'])]
                    stackSizer.Add(wx.StaticText(layerData,label=' Explicit layer sequence; enter space delimited list of numbers:'),0,WACV)
                    Str = ' Use sequence nos. from:'
                    for name in Names:
                        Str += name
                    stackSizer.Add(wx.StaticText(layerData,label=Str[:-1]+' Repeat sequences can be used: e.g. 6*(1 2) '),0,WACV)
                    stackSizer.Add(wx.StaticText(layerData,label=' Zero probability sequences not allowed'),0,WACV)    
                    stackList = wx.TextCtrl(layerData,value=Layers['Stacking'][2],size=(600,-1),
                        style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER)
                    stackList.Bind(wx.EVT_TEXT_ENTER,OnStackList)        
                    stackList.Bind(wx.EVT_KILL_FOCUS,OnStackList)
                    stackSizer.Add(stackList,0,wx.ALL|wx.EXPAND|WACV,8)
                else:   #random
                    topLine.Add(wx.StaticText(layerData,label=' Length of random sequence: '),0,WACV)
                    numRan = wx.TextCtrl(layerData,value=Layers['Stacking'][2],style=wx.TE_PROCESS_ENTER)
                    numRan.Bind(wx.EVT_TEXT_ENTER,OnNumRan)        
                    numRan.Bind(wx.EVT_KILL_FOCUS,OnNumRan)
                    topLine.Add(numRan,0,WACV)
                    stackSizer.Add(topLine,0,WACV)
            return stackSizer
            
        Layers = data['Layers']
        layerNames = []
        Layers['allowedTrans'] = []
        if len(Layers['Layers']):
            layerNames = [layer['Name'] for layer in Layers['Layers']]
        G2frame.GetStatusBar().SetStatusText('',1)
        layerData = G2frame.layerData
        if layerData.GetSizer():
            layerData.GetSizer().Clear(True)
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
        topSizer.Add(laueSizer,0,WACV)
        topSizer.Add(wx.StaticText(layerData,label=' Reference unit cell for all layers:'),0,WACV)
        topSizer.Add(CellSizer(),0,WACV)
        topSizer.Add(WidthSizer())
        topSizer.Add(wx.StaticText(layerData,label=' NB: stacking fault refinement currently not available'),0,WACV)
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
        topSizer.Add(titleSizer,0,WACV)
        for il,layer in enumerate(Layers['Layers']):
            topSizer.Add(LayerSizer(il,layer))
        G2G.HorizontalLine(topSizer,layerData)
        mainSizer.Add(topSizer)
        bottomSizer.Add(TransSizer())
        G2G.HorizontalLine(bottomSizer,layerData)
        bottomSizer.Add(PlotSizer(),0,WACV)
        G2G.HorizontalLine(bottomSizer,layerData)
        bottomSizer.Add(StackSizer())
        mainSizer.Add(bottomSizer)
        SetPhaseWindow(G2frame.layerData,mainSizer,Scroll=Scroll)
        
    def OnCopyPhase(event):
        dlg = wx.FileDialog(G2frame, 'Choose GSAS-II project file', 
            wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN| wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                GPXFile = dlg.GetPath()
                phaseNames = G2strIO.GetPhaseNames(GPXFile)
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
        General = G2strIO.GetAllPhaseData(GPXFile,PhaseName)['General']
        data['Layers']['Cell'] = General['Cell']
        wx.CallAfter(UpdateLayerData)

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
        dlg = wx.FileDialog(G2frame, 'Choose DIFFaX file name to read', '.', '',
            'DIFFaX file (*.*)|*.*',style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                DIFFaXfile = dlg.GetPath()
                data['Layers'] = G2IO.ReadDIFFaX(DIFFaXfile)
        finally:
            dlg.Destroy()
        wx.CallAfter(UpdateLayerData)
        
    def OnSimulate(event):
        debug = False       #set True to run DIFFax to compare/debug (must be in bin)
        idebug = 0
        if debug: idebug = 1
        ctrls = ''
        dlg = DIFFaXcontrols(G2frame,ctrls)
        if dlg.ShowModal() == wx.ID_OK:
            simCodes = dlg.GetSelection()
        else:
            return
        if 'PWDR' in  simCodes[0]:    #powder pattern
            data['Layers']['selInst'] = simCodes[1]
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
            G2plt.PlotPatterns(G2frame,plotType='PWDR',newPlot=True)
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
        wx.CallAfter(UpdateLayerData)
        
    def OnFitLayers(event):
        print (' fit stacking fault model TBD')
#        import scipy.optimize as opt
        wx.BeginBusyCursor()
#        Min,Init,Done = SetupPDFEval()
#        xstart = Init()
#        rms = Min(xstart)
#        print('Optimizing corrections to improve G(r) at low r')
#        print('start: Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f} (RMS:{:.2f})'.format(
#                data['Flat Bkg'],data['BackRatio'],data['Ruland'],rms))
#        
#        res = opt.minimize(Min,xstart,bounds=([0,None],[0,1],[0.01,1]),
#                           method='L-BFGS-B',options={'maxiter':5})
#        Done(res['x'])      
#        print('end:   Flat Bkg={:.1f}, BackRatio={:.3f}, Ruland={:.3f} (RMS:{:.2f})\n'.format(
#                data['Flat Bkg'],data['BackRatio'],data['Ruland'],res['fun']))
        wx.EndBusyCursor()
        wx.CallAfter(UpdateLayerData)
        G2plt.PlotPatterns(G2frame,plotType='PWDR')
        
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
        wx.CallAfter(UpdateLayerData)
        
################################################################################
#### Wave Data page
################################################################################

    def UpdateWavesData(Scroll=0):
        
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        typeNames = {'Sfrac':' Site fraction','Spos':' Position','Sadp':' Thermal motion','Smag':' Magnetic moment'}
        numVals = {'Sfrac':2,'Spos':6,'Sadp':12,'Smag':6,'ZigZag':5,'Block':5,'Crenel':2,'SawTooth':3}
        posNames = ['Xsin','Ysin','Zsin','Xcos','Ycos','Zcos','Tmin','Tmax','Xmax','Ymax','Zmax']
        adpNames = ['U11sin','U22sin','U33sin','U12sin','U13sin','U23sin',
            'U11cos','U22cos','U33cos','U12cos','U13cos','U23cos']
        magNames = ['MXsin','MYsin','MZsin','MXcos','MYcos','MZcos']
        fracNames = ['Fsin','Fcos','Fzero','Fwid']
        waveTypes = {'Sfrac':['Fourier','Crenel'],'Spos':['Fourier','ZigZag','Block','SawTooth'],'Sadp':['Fourier',],'Smag':['Fourier',]}
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
#            mainSizer.Detach(G2frame.bottomSizer)
            G2frame.bottomSizer.Clear(True)
            G2frame.bottomSizer = ShowAtomInfo()
            mainSizer.Add(G2frame.bottomSizer)
            mainSizer.Layout()
            G2frame.dataWindow.Refresh()
            waveData.SetVirtualSize(mainSizer.GetMinSize())
            waveData.Scroll(0,Scroll)
            G2frame.dataWindow.SendSizeEvent()
            
        def ShowAtomInfo():
            
            def AtomSizer(atom):
                
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
                        if waveTyp in ['ZigZag','Block','Crenel','SawTooth']:
                            nt = numVals[waveTyp]                        
                        atm[-1]['SS1'][item] = [0,]
                        atm[-1]['SS1'][item][0] = waveType.GetValue()
                    atm[-1]['SS1'][item].append([[0.0 for i in range(nt)],False])
                    wx.CallAfter(RepaintAtomInfo,G2frame.waveData.GetScrollPos(wx.VERTICAL))
                    
                def OnWaveVal(event):
                    event.Skip()
                    Obj = event.GetEventObject()
                    item,iwave,ival = Indx[Obj.GetId()]
                    try:
                        val = float(Obj.GetValue())
                        if waveTyp in ['ZigZag','Block'] and ival < 2 and not iwave:
                            if ival == 1: #Tmax
                                val = min(1.0,max(0.0,val))
                            elif ival == 0: #Tmin
                                val = max(-1.,min(val,atm[-1]['SS1'][item][1][0][ival]))
                    except ValueError:
                        val = atm[-1]['SS1'][item][iwave+1][0][ival]
                    Obj.SetValue('%.5f'%val)
                    atm[-1]['SS1'][item][iwave+1][0][ival] = val
                    
                def OnRefWave(event):
                    Obj = event.GetEventObject()
                    item,iwave = Indx[Obj.GetId()]
                    atm[-1]['SS1'][item][iwave+1][1] = not atm[-1]['SS1'][item][iwave+1][1]
                    
                def OnDelWave(event):
                    Obj = event.GetEventObject()
                    item,iwave = Indx[Obj.GetId()]
                    del atm[-1]['SS1'][item][iwave+1]
                    wx.CallAfter(RepaintAtomInfo,G2frame.waveData.GetScrollPos(wx.VERTICAL))                
                
                waveTyp,waveBlk = 'Fourier',[]
                if len(wavedata):
                    waveTyp = wavedata[0]
                    waveBlk = wavedata[1:]
                waveSizer = wx.BoxSizer(wx.VERTICAL)
                waveHead = wx.BoxSizer(wx.HORIZONTAL)
                waveHead.Add(wx.StaticText(waveData,label=typeName+' modulation parameters: '),0,WACV)
                waveAdd = wx.CheckBox(waveData,label='Add wave?   WaveType: ')
                waveAdd.Bind(wx.EVT_CHECKBOX, OnAddWave)
                Indx[waveAdd.GetId()] = Stype
                waveHead.Add(waveAdd,0,WACV)
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
                            if waveTyp in ['ZigZag','Block','SawTooth','Crenel']:
                                nx = 1
                            CSI = G2spc.GetSSfxuinel(waveTyp,Stype,1,xyz,SGData,SSGData)
                        else:
                            CSI = G2spc.GetSSfxuinel('Fourier',Stype,iwave+1-nx,xyz,SGData,SSGData)
                        waveName = 'Fourier'
                        if Stype == 'Sfrac':
                            if 'Crenel' in waveTyp and not iwave:
                                waveName = 'Crenel'
                                names = Names[2:]
                            else:
                                names = Names[:2]
                            Waves = wx.FlexGridSizer(0,4,5,5)
                        elif Stype == 'Spos':
                            if waveTyp in ['ZigZag','Block'] and not iwave:
                                names = Names[6:]
                                Waves = wx.FlexGridSizer(0,7,5,5)
                                waveName = waveTyp
                            else:
                                names = Names[:6]
                                Waves = wx.FlexGridSizer(0,8,5,5)
                        else:
                            names = Names
                            Waves = wx.FlexGridSizer(0,8,5,5)
                        waveSizer.Add(wx.StaticText(waveData,label=' %s  parameters: %s'%(waveName,str(names).rstrip(']').lstrip('[').replace("'",''))),0,WACV)
                        for ival,val in enumerate(wave[0]):
                            if np.any(CSI[0][ival]):
                                waveVal = wx.TextCtrl(waveData,value='%.5f'%(val),style=wx.TE_PROCESS_ENTER)
                                waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveVal)
                                waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveVal)
                                Indx[waveVal.GetId()] = [Stype,iwave,ival]
                            else:
                                waveVal = wx.TextCtrl(waveData,value='%.5f'%(val),style=wx.TE_READONLY)
                                waveVal.SetBackgroundColour(VERY_LIGHT_GREY)
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
                            waveDel = wx.CheckBox(waveData,label='Delete?')
                            Indx[waveDel.GetId()] = [Stype,iwave]
                            waveDel.Bind(wx.EVT_CHECKBOX, OnDelWave)
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
        topSizer.Add(atms,0,WACV)
        mainSizer.Add(topSizer,0,WACV)
        G2frame.bottomSizer = ShowAtomInfo()
        mainSizer.Add(G2frame.bottomSizer)
        SetPhaseWindow(G2frame.waveData,mainSizer,Scroll=Scroll)
    
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
        UpdateWavesData()

################################################################################
#### Structure drawing GUI stuff                
################################################################################

    def SetupDrawingData():
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = data['Atoms']
        defaultDrawing = {'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,
            'backColor':[0,0,0],'depthFog':False,'Zclip':50.0,'cameraPos':50.,'Zstep':0.5,
            'radiusFactor':0.85,'contourLevel':1.,'bondRadius':0.1,'ballScale':0.33,
            'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,'unitCellBox':True,
            'showABC':True,'selectedAtoms':[],'Atoms':[],'oldxy':[],'magMult':1.0,
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
        if 'viewDir' not in drawingData:
            drawingData['viewDir'] = [0,0,1]
        if 'Quaternion' not in drawingData:
            drawingData['Quaternion'] = G2mth.AV2Q(2*np.pi,np.inner(Amat,[0,0,1]))
        if 'showRigidBodies' not in drawingData:
            drawingData['showRigidBodies'] = True
        if 'Plane' not in drawingData:
            drawingData['Plane'] = [[0,0,1],False,False,0.0,[255,255,0]]
        if 'magMult' not in drawingData:
            drawingData['magMult'] = 1.0
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
        drawingData['Atoms'].append(MakeDrawAtom(atom))
        
    def OnRestraint(event):        
        indx = drawAtoms.GetSelectedRows()
        restData = G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints'))
        drawingData = data['Drawing']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])            
        cx,ct,cs,ci = drawingData['atomPtrs']
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
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        RBData = G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        drawingData = data['Drawing']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])            
        cx,ct,cs,ci = drawingData['atomPtrs']
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

################################################################################
##### Draw Atom routines
################################################################################
            
    def UpdateDrawAtoms(atomStyle=''):
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
                sel = -1
                if drawAtoms.GetColLabelValue(c) == 'Style':
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom drawing style',styleChoice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = styleChoice[sel]
                        for r in range(len(atomData)):
                            atomData[r][c] = parms
                            drawAtoms.SetCellValue(r,c,parms)
                        FindBondsDraw(data)
                        G2plt.PlotStructure(G2frame,data)
                    dlg.Destroy()
                elif drawAtoms.GetColLabelValue(c) == 'Label':
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom labelling style',labelChoice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = labelChoice[sel]
                        for r in range(len(atomData)):
                            atomData[r][c] = parms
                            drawAtoms.SetCellValue(r,c,parms)
                    dlg.Destroy()                    
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    dlg = wx.ColourDialog(G2frame)
                    if dlg.ShowModal() == wx.ID_OK:
                        color = dlg.GetColourData().GetColour()
                        attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                        attr.SetReadOnly(True)
                        attr.SetBackgroundColour(color)
                        for r in range(len(atomData)):
                            atomData[r][c] = color
                            drawingData['Atoms'][r][c] = color
                            drawAtoms.SetAttr(r,c,attr)
                        UpdateDrawAtoms()
                    dlg.Destroy()
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
                    dlg = wx.ColourDialog(G2frame)
                    if dlg.ShowModal() == wx.ID_OK:
                        color = dlg.GetColourData().GetColour()
                        attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                        attr.SetReadOnly(True)
                        attr.SetBackgroundColour(color)
                        atomData[r][c] = color
                        drawingData['Atoms'][r][c] = color
                        drawAtoms.SetAttr(i,cs+2,attr)
                    dlg.Destroy()
                    UpdateDrawAtoms()
            G2plt.PlotStructure(G2frame,data)
                    
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
                    drawAtoms.ClearSelection()
                    drawAtoms.SelectRow(r,True)                
            drawingData['selectedAtoms'] = []
            drawingData['selectedAtoms'] = drawAtoms.GetSelectedRows()
            G2plt.PlotStructure(G2frame,data)                    

        # UpdateDrawAtoms executable code starts here
        G2frame.GetStatusBar().SetStatusText('',1)
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
        styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']
        labelChoice = [' ','type','name','number']
        colLabels = ['Name','Type','x','y','z','Sym Op','Style','Label','Color','I/A']
        if generalData['Type'] == 'macromolecular':
            colLabels = ['Residue','1-letter','Chain'] + colLabels
            Types = 3*[wg.GRID_VALUE_STRING,]+Types
            Types[8] = wg.GRID_VALUE_CHOICE+": ,lines,vdW balls,sticks,balls & sticks,ellipsoids,backbone,ribbons,schematic"
            styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','backbone','ribbons','schematic']
            labelChoice = [' ','type','name','number','residue','1-letter','chain']
            Types[9] = wg.GRID_VALUE_CHOICE+": ,type,name,number,residue,1-letter,chain"
        elif generalData['Type'] == 'magnetic':
            colLabels = colLabels[:5]+['Mx','My','Mz']+colLabels[5:]
            Types = Types[:5]+3*[wg.GRID_VALUE_FLOAT+':10,4',]+Types[5:]
        table = []
        rowLabels = []
        for i,atom in enumerate(drawingData['Atoms']):
            table.append(atom[:colLabels.index('I/A')+1])
            rowLabels.append(str(i))

        G2frame.atomTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        drawAtoms.SetTable(G2frame.atomTable, True)
        drawAtoms.SetMargins(0,0)
        drawAtoms.AutoSizeColumns(True)
        drawAtoms.SetColSize(colLabels.index('Style'),80)
        drawAtoms.SetColSize(colLabels.index('Color'),50)
        if 'phoenix' in wx.version():
            drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGED, RefreshDrawAtomGrid)
        else:
            drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshDrawAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshDrawAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, RefreshDrawAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
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
        SetPhaseWindow(drawAtoms)

        FindBondsDraw(data)
        drawAtoms.ClearSelection()
#        G2plt.PlotStructure(G2frame,data)

    def DrawAtomStyle(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
            styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']
            if generalData['Type'] == 'macromolecular':
                styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids',
                'backbone','ribbons','schematic']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom drawing style',styleChoice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parms = styleChoice[sel]
                for r in indx:
                    atomData[r][cs] = parms
                    drawAtoms.SetCellValue(r,cs,parms)
            dlg.Destroy()
            FindBondsDraw(data)
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)

    def DrawAtomLabel(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
            styleChoice = [' ','type','name','number']
            if generalData['Type'] == 'macromolecular':
                styleChoice = [' ','type','name','number','residue','1-letter','chain']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom label style',styleChoice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parms = styleChoice[sel]
                for r in indx:
                    atomData[r][cs+1] = parms
                    drawAtoms.SetCellValue(r,cs+1,parms)
            dlg.Destroy()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def DrawAtomColor(event):

        indx = drawAtoms.GetSelectedRows()
        if indx:
            if len(indx) > 1:
                G2frame.GetStatusBar().SetStatusText('Select Custom Color, change color, Add to Custom Colors, then OK',1)
            else:
                G2frame.GetStatusBar().SetStatusText('Change color, Add to Custom Colors, then OK',1)
            atomData = data['Drawing']['Atoms']
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
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
            dlg = wx.ColourDialog(G2frame)
            if dlg.ShowModal() == wx.ID_OK:
                for i in range(len(atmColors)):                    
                    atmColors[i] = dlg.GetColourData().GetColour()
                colorDict = dict(zip(atmTypes,atmColors))
                for r in indx:
                    color = colorDict[atomData[r][ct]]
                    atomData[r][cs+2] = color
                    attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                    attr.SetBackgroundColour(color)
                    drawAtoms.SetAttr(r,cs+2,attr)
                    data['Drawing']['Atoms'][r][cs+2] = color
            drawAtoms.ClearSelection()
            dlg.Destroy()
            G2frame.GetStatusBar().SetStatusText('',1)
            G2plt.PlotStructure(G2frame,data)
            
    def ResetAtomColors(event):
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        for atom in atomData:            
            atNum = generalData['AtomTypes'].index(atom[ct])
            atom[cs+2] = list(generalData['Color'][atNum])
        UpdateDrawAtoms()
        drawAtoms.ClearSelection()
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
        indx = drawAtoms.GetSelectedRows()
        if indx:
            atomData = data['Drawing']['Atoms']
            cx = data['Drawing']['atomPtrs'][0]
            data['Drawing']['viewPoint'] = [atomData[indx[0]][cx:cx+3],[indx[0],0]]
            drawAtoms.ClearSelection()                                  #do I really want to do this?
            G2plt.PlotStructure(G2frame,data)
            
    def noDuplicate(xyz,atomData):                  #be careful where this is used - it's slow
        cx = data['Drawing']['atomPtrs'][0]
        if True in [np.allclose(np.array(xyz),np.array(atom[cx:cx+3]),atol=0.0002) for atom in atomData]:
            return False
        else:
            return True
                
    def AddSymEquiv(event):
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
            colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
            cx,ct,cs,cui = data['Drawing']['atomPtrs']
            cuij = cui+2
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
                            if atom[cui] == 'A':
                                Uij = atom[cuij:cuij+6]
                                Uij = G2spc.U2Uij(np.inner(np.inner(M,G2spc.Uij2U(Uij)),M))
                                atom[cuij:cuij+6] = Uij
                            atomData.append(atom[:cuij+9])  #not SS stuff
            finally:
                dlg.Destroy()
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def AddSphere(event):
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = data['Drawing']['Atoms']
        numAtoms = len(atomData)
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
        cuij = cs+5
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cmx = 0
        if 'Mx' in colLabels:
            cmx = colLabels.index('Mx')
        generalData = data['General']
        SGData = generalData['SGData']
        SpnFlp = SGData.get('SpnFlp',[])
        cellArray = G2lat.CellBlock(1)
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        dlg = SphereEnclosure(G2frame,data['General'],data['Drawing'],indx)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                centers,radius,targets = dlg.GetSelection()
                for orig in centers:
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
        finally:
            dlg.Destroy()
        UpdateDrawAtoms()
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
            
    def TransformSymEquiv(event):
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
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
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def FillCoordSphere(event):
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        radii = generalData['BondRadii']
        atomTypes = generalData['AtomTypes']
        try:
            indH = atomTypes.index('H')
            radii[indH] = 0.5
        except:
            pass            
        indx = drawAtoms.GetSelectedRows()
        if indx:
            indx.sort()
            atomData = data['Drawing']['Atoms']
            numAtoms = len(atomData)
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
            cij = ci+2
            SGData = generalData['SGData']
            cellArray = G2lat.CellBlock(1)
            wx.BeginBusyCursor()
            try:
                for ind in indx:
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
            finally:
                wx.EndBusyCursor()
            data['Drawing']['Atoms'] = atomData
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def FillUnitCell(event):
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
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
            wx.BeginBusyCursor()
            try:
                for ind in indx:
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
                            if noDuplicate(Opp[key],atomData):
                                unit = item[3]+np.array(eval(key))*1.
                                cell = '%d+%d,%d,%d'%(item[2],unit[0],unit[1],unit[2])
                                atom[cx:cx+3] = Opp[key]
                                atom[cs-1] = cell
                                atomData.append(atom[:cuij+9])  #not SS stuff
                    data['Drawing']['Atoms'] = atomData
            finally:
                wx.EndBusyCursor()
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def DrawAtomsDelete(event):   
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
            atomData = data['Drawing']['Atoms']
            indx.reverse()
            for ind in indx:
                del atomData[ind]
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()
        
    def OnReloadDrawAtoms(event):
        atomData = data['Atoms']
        cx,ct,cs,ci = data['General']['AtomPtrs']
        for atom in atomData:
            ID = atom[ci+8]
            DrawAtomsReplaceByID(data['Drawing'],ci+8,atom,ID)
        UpdateDrawAtoms()
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
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
        indx = drawAtoms.GetSelectedRows()
        if len(indx) < 4:
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
        indx = drawAtoms.GetSelectedRows()
        if not indx:
            print ('***** ERROR - no atoms selected')
            return
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
        #distance, angle, torsion 
        indx = drawAtoms.GetSelectedRows()
        if len(indx) not in [2,3,4]:
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
            id = G2mth.FindAtomIndexByIDs(atomData,cid,[atom[cid],],False)[0]
            Oxyz.append([id,]+atomData[id][cx+1:cx+4])
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
                        
################################################################################
#### Draw Options page
################################################################################

    def UpdateDrawOptions():
        import wx.lib.colourselect as wcs
        def SlopSizer():            
            def OnCameraPos(event):
                drawingData['cameraPos'] = cameraPos.GetValue()
                cameraPosTxt.SetLabel(' Camera Distance: '+'%.2f'%(drawingData['cameraPos']))
                ZclipTxt.SetLabel(' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
                G2plt.PlotStructure(G2frame,data)

            def OnZclip(event):
                drawingData['Zclip'] = Zclip.GetValue()
                ZclipTxt.SetLabel(' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
                G2plt.PlotStructure(G2frame,data)
                
            def OnMoveZ(event):
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
                G2plt.PlotStructure(G2frame,data)
                
            def OnVdWScale(event):
                drawingData['vdwScale'] = vdwScale.GetValue()/100.
                vdwScaleTxt.SetLabel(' van der Waals scale: '+'%.2f'%(drawingData['vdwScale']))
                G2plt.PlotStructure(G2frame,data)
    
            def OnEllipseProb(event):
                drawingData['ellipseProb'] = ellipseProb.GetValue()
                ellipseProbTxt.SetLabel(' Ellipsoid probability: '+'%d%%'%(drawingData['ellipseProb']))
                G2plt.PlotStructure(G2frame,data)
    
            def OnBallScale(event):
                drawingData['ballScale'] = ballScale.GetValue()/100.
                ballScaleTxt.SetLabel(' Ball scale: '+'%.2f'%(drawingData['ballScale']))
                G2plt.PlotStructure(G2frame,data)

            def OnBondRadius(event):
                drawingData['bondRadius'] = bondRadius.GetValue()/100.
                bondRadiusTxt.SetLabel(' Bond radius, A: '+'%.2f'%(drawingData['bondRadius']))
                G2plt.PlotStructure(G2frame,data)

            def OnMagMult(event):
                drawingData['magMult'] = magMult.GetValue()/100.
                magMultTxt.SetLabel(' Mag. mom. mult.: '+'%.2f'%(drawingData['magMult']))
                G2plt.PlotStructure(G2frame,data)
                
            def OnContourLevel(event):
                drawingData['contourLevel'] = contourLevel.GetValue()/100.
                contourLevelTxt.SetLabel(' Contour level: '+'%.2f'%(drawingData['contourLevel']*generalData['Map']['rhoMax']))
                G2plt.PlotStructure(G2frame,data)

            def OnMapSize(event):
                drawingData['mapSize'] = mapSize.GetValue()/10.
                mapSizeTxt.SetLabel(' Map radius, A: '+'%.1f'%(drawingData['mapSize']))
                G2plt.PlotStructure(G2frame,data)

            
            slopSizer = wx.BoxSizer(wx.HORIZONTAL)
            slideSizer = wx.FlexGridSizer(0,2,0,0)
            slideSizer.AddGrowableCol(1,1)
    
            cameraPosTxt = wx.StaticText(drawOptions,-1,
                ' Camera Distance: '+'%.2f'%(drawingData['cameraPos']),name='cameraPos')
            G2frame.phaseDisplay.cameraPosTxt = cameraPosTxt
            slideSizer.Add(cameraPosTxt,0,WACV)
            cameraPos = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=drawingData['cameraPos'],name='cameraSlider')
            cameraPos.SetRange(10,500)
            cameraPos.Bind(wx.EVT_SLIDER, OnCameraPos)
            G2frame.phaseDisplay.cameraSlider = cameraPos
            slideSizer.Add(cameraPos,1,wx.EXPAND|wx.RIGHT)
            
            ZclipTxt = wx.StaticText(drawOptions,-1,' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
            slideSizer.Add(ZclipTxt,0,WACV)
            Zclip = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=drawingData['Zclip'])
            Zclip.SetRange(1,99)
            Zclip.Bind(wx.EVT_SLIDER, OnZclip)
            slideSizer.Add(Zclip,1,wx.EXPAND|wx.RIGHT)
            
            ZstepSizer = wx.BoxSizer(wx.HORIZONTAL)
            ZstepSizer.Add(wx.StaticText(drawOptions,-1,' Z step:'),0,WACV)
            Zstep = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'Zstep',nDig=(10,2),min=0.01,max=1.0)
            ZstepSizer.Add(Zstep,0,WACV)
            slideSizer.Add(ZstepSizer)
            MoveSizer = wx.BoxSizer(wx.HORIZONTAL)
            MoveSizer.Add(wx.StaticText(drawOptions,-1,'   Press to step:'),0,WACV)
            MoveZ = wx.SpinButton(drawOptions,style=wx.SP_HORIZONTAL,size=wx.Size(100,20))
            MoveZ.SetValue(0)
            MoveZ.SetRange(-1,1)
            MoveZ.Bind(wx.EVT_SPIN, OnMoveZ)
            MoveSizer.Add(MoveZ)
            slideSizer.Add(MoveSizer,1,wx.EXPAND|wx.RIGHT)
            
            vdwScaleTxt = wx.StaticText(drawOptions,-1,' van der Waals scale: '+'%.2f'%(drawingData['vdwScale']))
            slideSizer.Add(vdwScaleTxt,0,WACV)
            vdwScale = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(100*drawingData['vdwScale']))
            vdwScale.Bind(wx.EVT_SLIDER, OnVdWScale)
            slideSizer.Add(vdwScale,1,wx.EXPAND|wx.RIGHT)
    
            ellipseProbTxt = wx.StaticText(drawOptions,-1,' Ellipsoid probability: '+'%d%%'%(drawingData['ellipseProb']))
            slideSizer.Add(ellipseProbTxt,0,WACV)
            ellipseProb = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=drawingData['ellipseProb'])
            ellipseProb.SetRange(1,99)
            ellipseProb.Bind(wx.EVT_SLIDER, OnEllipseProb)
            slideSizer.Add(ellipseProb,1,wx.EXPAND|wx.RIGHT)
    
            ballScaleTxt = wx.StaticText(drawOptions,-1,' Ball scale: '+'%.2f'%(drawingData['ballScale']))
            slideSizer.Add(ballScaleTxt,0,WACV)
            ballScale = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(100*drawingData['ballScale']))
            ballScale.Bind(wx.EVT_SLIDER, OnBallScale)
            slideSizer.Add(ballScale,1,wx.EXPAND|wx.RIGHT)
    
            bondRadiusTxt = wx.StaticText(drawOptions,-1,' Bond radius, A: '+'%.2f'%(drawingData['bondRadius']))
            slideSizer.Add(bondRadiusTxt,0,WACV)
            bondRadius = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(100*drawingData['bondRadius']))
            bondRadius.SetRange(1,25)
            bondRadius.Bind(wx.EVT_SLIDER, OnBondRadius)
            slideSizer.Add(bondRadius,1,wx.EXPAND|wx.RIGHT)
            
            if generalData['Type'] == 'magnetic':
                magMultTxt = wx.StaticText(drawOptions,-1,' Mag. mom. mult.: '+'%.2f'%(drawingData['magMult']))
                slideSizer.Add(magMultTxt,0,WACV)
                magMult = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(100*drawingData['magMult']))
                magMult.SetRange(10,500)
                magMult.Bind(wx.EVT_SLIDER, OnMagMult)
                slideSizer.Add(magMult,1,wx.EXPAND|wx.RIGHT)
            
            if generalData['Map']['rhoMax']:
                contourLevelTxt = wx.StaticText(drawOptions,-1,' Contour level: '+'%.2f'%(drawingData['contourLevel']*generalData['Map']['rhoMax']))
                slideSizer.Add(contourLevelTxt,0,WACV)
                contourLevel = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(100*drawingData['contourLevel']))
                contourLevel.SetRange(1,100)
                contourLevel.Bind(wx.EVT_SLIDER, OnContourLevel)
                slideSizer.Add(contourLevel,1,wx.EXPAND|wx.RIGHT)
                mapSizeTxt = wx.StaticText(drawOptions,-1,' Map radius, A: '+'%.1f'%(drawingData['mapSize']))
                slideSizer.Add(mapSizeTxt,0,WACV)
                mapSize = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=int(10*drawingData['mapSize']))
                mapSize.SetRange(1,100)
                mapSize.Bind(wx.EVT_SLIDER, OnMapSize)
                slideSizer.Add(mapSize,1,wx.EXPAND|wx.RIGHT)
            
            slopSizer.Add(slideSizer,1,wx.EXPAND|wx.RIGHT)
            slopSizer.Add((10,5),0)
            slopSizer.SetMinSize(wx.Size(350,10))
            return slopSizer
            
        def ShowSizer():
            
            def OnBackColor(event):
                drawingData['backColor'] = list(event.GetValue())[:3]
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowABC(event):
                drawingData['showABC'] = showABC.GetValue()
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowUnitCell(event):
                drawingData['unitCellBox'] = unitCellBox.GetValue()
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowHyd(event):
                drawingData['showHydrogen'] = showHydrogen.GetValue()
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data)
                
            def OnShowRB(event):
                drawingData['showRigidBodies'] = showRB.GetValue()
                FindBondsDraw(data)
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
                    VD = np.array([float(viewDir[i]) for i in range(3)])
                    VC = np.inner(Amat,VD)
                    VC /= np.sqrt(np.sum(VC**2))
                    V = np.array(drawingData['viewDir'])
                    VB = np.inner(Amat,V)
                    VB /= np.sqrt(np.sum(VB**2))
                    VX = np.cross(VC,VB)
                    A = acosd(max((2.-np.sum((VB-VC)**2))/2.,-1.))
                    QV = G2mth.AVdeg2Q(A,VX)
                    Q = drawingData['Quaternion']
                    drawingData['Quaternion'] = G2mth.prodQQ(Q,QV)
                except (ValueError,IndexError):
                    VD = drawingData['viewDir']
                Obj.SetValue('%.3f %.3f %.3f'%(VD[0],VD[1],VD[2]))
                drawingData['viewDir'] = VD
                G2plt.PlotStructure(G2frame,data)
                                
            showSizer = wx.BoxSizer(wx.VERTICAL)            
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(drawOptions,-1,' Background color:'),0,WACV)
            backColor = wcs.ColourSelect(drawOptions, -1,colour=drawingData['backColor'],size=wx.Size(25,25))
            backColor.Bind(wcs.EVT_COLOURSELECT, OnBackColor)
            lineSizer.Add(backColor,0,WACV)
            lineSizer.Add(wx.StaticText(drawOptions,-1,' View Dir.:'),0,WACV)
            VD = drawingData['viewDir']
            viewDir = wx.TextCtrl(drawOptions,value='%.3f %.3f %.3f'%(VD[0],VD[1],VD[2]),
                style=wx.TE_PROCESS_ENTER,size=wx.Size(140,20),name='viewDir')
            viewDir.Bind(wx.EVT_TEXT_ENTER,OnViewDir)
            viewDir.Bind(wx.EVT_KILL_FOCUS,OnViewDir)
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
            viewPoint = wx.TextCtrl(drawOptions,value='%.3f %.3f %.3f'%(VP[0],VP[1],VP[2]),
                style=wx.TE_PROCESS_ENTER,size=wx.Size(140,20),name='viewPoint')
            G2frame.phaseDisplay.viewPoint = viewPoint
            viewPoint.Bind(wx.EVT_TEXT_ENTER,OnViewPoint)
            viewPoint.Bind(wx.EVT_KILL_FOCUS,OnViewPoint)
            lineSizer.Add(viewPoint,0,WACV)
            showSizer.Add(lineSizer)
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
            
            showRB = wx.CheckBox(drawOptions,-1,label=' Show rigid Bodies?')
            showRB.Bind(wx.EVT_CHECKBOX, OnShowRB)
            showRB.SetValue(drawingData['showRigidBodies'])
            line2Sizer.Add(showRB,0,WACV)
            
            showSizer.Add(line2Sizer)
            return showSizer
            
        def RadSizer():
            
            def OnSizeHatoms(invalid,value,tc):
                G2plt.PlotStructure(G2frame,data)
                
            def OnRadFactor(invalid,value,tc):
                FindBondsDraw(data)
                G2plt.PlotStructure(G2frame,data)
            
            radSizer = wx.BoxSizer(wx.HORIZONTAL)
            radSizer.Add(wx.StaticText(drawOptions,-1,' Hydrogen radius, A:  '),0,WACV)
            sizeH = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'sizeH',nDig=(10,2),min=0.1,max=1.2,size=wx.Size(60,20),OnLeave=OnSizeHatoms)
            radSizer.Add(sizeH,0,WACV)
    
            radSizer.Add(wx.StaticText(drawOptions,-1,' Bond search factor:  '),0,WACV)
            radFactor = G2G.ValidatedTxtCtrl(drawOptions,drawingData,'radiusFactor',nDig=(10,2),min=0.1,max=1.2,size=wx.Size(60,20),OnLeave=OnRadFactor)
            radSizer.Add(radFactor,0,WACV)
            return radSizer
            
        def PlaneSizer():
            
            def OnPlane(event):
                event.Skip()
                vals = plane.GetValue().split()
                try:
                    hkl = [int(vals[i]) for i in range(3)]
                    if not np.any(np.array(hkl)):       #can't be all zeros!
                        raise ValueError
                except (ValueError,IndexError):
                    hkl = drawingData['Plane'][0]
                drawingData['Plane'][0] = hkl
                plane.SetValue('%3d %3d %3d'%(hkl[0],hkl[1],hkl[2]))
                G2plt.PlotStructure(G2frame,data)
                
            def OnShowPlane(event):
                drawingData['Plane'][1] = showPlane.GetValue()
                G2plt.PlotStructure(G2frame,data)
                
            def OnShowStack(event):
                drawingData['Plane'][2] = showStack.GetValue()
                G2plt.PlotStructure(G2frame,data)
                
            def OnPhase(invalid,value,tc):
                G2plt.PlotStructure(G2frame,data)
            
            def OnPlaneColor(event):
                drawingData['Plane'][4] = list(event.GetValue())[:3]
                G2plt.PlotStructure(G2frame,data)

            planeSizer = wx.BoxSizer(wx.VERTICAL)
            planeSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            planeSizer1.Add(wx.StaticText(drawOptions,label=' Plane: '),0,WACV)
            H = drawingData['Plane'][0]
            plane = wx.TextCtrl(drawOptions,value='%3d %3d %3d'%(H[0],H[1],H[2]),
                style=wx.TE_PROCESS_ENTER)
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
            planeColor = wcs.ColourSelect(drawOptions, -1,colour=drawingData['Plane'][4],size=wx.Size(25,25))
            planeColor.Bind(wcs.EVT_COLOURSELECT, OnPlaneColor)
            planeSizer2.Add(planeColor,0,WACV)
            planeSizer.Add(planeSizer1)
            planeSizer.Add(planeSizer2)
            return planeSizer
            

        # UpdateDrawOptions exectable code starts here
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        SetupDrawingData()
        drawingData = data['Drawing']
        SetDrawingDefaults(drawingData)        

        G2frame.GetStatusBar().SetStatusText('',1)
        if drawOptions.GetSizer():
            drawOptions.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(drawOptions,-1,' Drawing controls:'),0,WACV)
        mainSizer.Add((5,5),0)        
        mainSizer.Add(SlopSizer(),0)
        mainSizer.Add((5,5),0)
        mainSizer.Add(ShowSizer(),0,)
        mainSizer.Add((5,5),0)
        mainSizer.Add(RadSizer(),0,)
        mainSizer.Add((5,5),0)
        mainSizer.Add(PlaneSizer(),0,)
        SetPhaseWindow(drawOptions,mainSizer)

################################################################################
####  Texture routines
################################################################################
        
    def UpdateTexture():
                
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
            textureData['Order'] = int(Obj.GetValue())
            textureData['SH Coeff'][1] = SetSHCoef()
            wx.CallLater(100,UpdateTexture)
            wx.CallAfter(G2plt.PlotTexture,G2frame,data)
                        
        def OnShModel(event):
            Obj = event.GetEventObject()
            textureData['Model'] = Obj.GetValue()
            textureData['SH Coeff'][1] = SetSHCoef()
            wx.CallLater(100,UpdateTexture)
            wx.CallAfter(G2plt.PlotTexture,G2frame,data)
            
        def OnSHRefine(event):
            Obj = event.GetEventObject()
            textureData['SH Coeff'][0] = Obj.GetValue()
            
        def OnSHShow(event):
            Obj = event.GetEventObject()
            textureData['SHShow'] = Obj.GetValue()
            wx.CallLater(100,UpdateTexture)
            
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
            wx.CallLater(100,UpdateTexture)
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
            ODFln = G2lat.Flnh(True,SHCoef,phi,beta,SGData)
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
                IODFln = G2lat.Glnh(True,SHCoef,psi,gam,SamSym[textureData['Model']])
                pfName = PhaseName+'%d%d%dIPF.csv'%(int(PX[0]),int(PX[1]),int(PX[2]))
                pth = G2G.GetExportPath(G2frame)
                dlg = wx.FileDialog(G2frame, 'Choose CSV inverse pole figure file name', pth, pfName, 
                    'CSV file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            else:    
                PH = np.array(textureData['PFhkl'])
                phi,beta = G2lat.CrsAng(PH,cell,SGData)
                SHCoef = textureData['SH Coeff'][1]
                ODFln = G2lat.Flnh(True,SHCoef,phi,beta,SGData)
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
                wx.CallLater(100,UpdateTexture)
                
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
            shToler = G2G.ValidatedTxtCtrl(Texture,Penalty,1,nDig=(10,2),min=0.001)
            shPenalty.Add(shToler,0,WACV)
            return shPenalty
        
        def OnProj(event):
            Obj = event.GetEventObject()
            generalData['3Dproj'] = Obj.GetValue()
            wx.CallAfter(G2plt.PlotTexture,G2frame,data)            
        
        # UpdateTexture executable starts here
        generalData = data['General']        
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
        titleSizer = wx.BoxSizer(wx.HORIZONTAL)
        titleSizer.Add(wx.StaticText(Texture,-1,' Spherical harmonics texture data for '+PhaseName+':'),0,WACV)
        titleSizer.Add(wx.StaticText(Texture,-1,
            ' Texture Index J = %7.3f'%(G2lat.textureIndex(textureData['SH Coeff'][1]))),
            0,WACV)
        mainSizer.Add(titleSizer,0)
        mainSizer.Add((0,5),0)
        shSizer = wx.FlexGridSizer(0,6,5,5)
        shSizer.Add(wx.StaticText(Texture,-1,' Texture model: '),0,WACV)
        shModel = wx.ComboBox(Texture,-1,value=textureData['Model'],choices=shModels,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        shModel.Bind(wx.EVT_COMBOBOX,OnShModel)
        shSizer.Add(shModel,0,WACV)
        shSizer.Add(wx.StaticText(Texture,-1,'  Harmonic order: '),0,WACV)
        shOrder = wx.ComboBox(Texture,-1,value=str(textureData['Order']),choices=[str(2*i) for i in range(18)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        shOrder.Bind(wx.EVT_COMBOBOX,OnShOrder)
        shSizer.Add(shOrder,0,WACV)
        shRef = wx.CheckBox(Texture,-1,label=' Refine texture?')
        shRef.SetValue(textureData['SH Coeff'][0])
        shRef.Bind(wx.EVT_CHECKBOX, OnSHRefine)
        shSizer.Add(shRef,0,WACV)
        shShow = wx.CheckBox(Texture,-1,label=' Show coeff.?')
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
            pfVal = wx.TextCtrl(Texture,-1,'%d %d %d'%(PH[0],PH[1],PH[2]),style=wx.TE_PROCESS_ENTER)
        else:
            PTSizer.Add(wx.StaticText(Texture,-1,' Inverse pole figure XYZ: '),0,WACV)
            PX = textureData['PFxyz']
            pfVal = wx.TextCtrl(Texture,-1,'%3.1f %3.1f %3.1f'%(PX[0],PX[1],PX[2]),style=wx.TE_PROCESS_ENTER)
        pfVal.Bind(wx.EVT_TEXT_ENTER,OnPFValue)
        pfVal.Bind(wx.EVT_KILL_FOCUS,OnPFValue)
        PTSizer.Add(pfVal,0,WACV)
        if 'Axial' not in textureData['PlotType'] and '3D' not in textureData['PlotType']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Color scheme'),0,WACV)
            choice = [m for m in mpl.cm.datad.keys()]       #if not m.endswith("_r")
            choice.sort()
            colorSel = wx.ComboBox(Texture,-1,value=G2frame.ContourColor,choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            colorSel.Bind(wx.EVT_COMBOBOX,OnColorSel)
            PTSizer.Add(colorSel,0,WACV)
            if 'figure' in textureData['PlotType']:
                popLA = wx.Button(Texture,-1,"Make CSV file")
                popLA.Bind(wx.EVT_BUTTON, OnCSV)
                PTSizer.Add(popLA,0,WACV)
        mainSizer.Add(PTSizer,0,WACV)
        mainSizer.Add((0,5),0)
        if textureData['SHShow']:
            mainSizer.Add(wx.StaticText(Texture,-1,' Spherical harmonic coefficients: '),0,WACV)
            mainSizer.Add((0,5),0)
            ODFSizer = wx.FlexGridSizer(0,8,2,2)
            ODFkeys = list(textureData['SH Coeff'][1].keys())
            ODFkeys.sort()
            for item in ODFkeys:
                ODFSizer.Add(wx.StaticText(Texture,-1,item),0,WACV)
                ODFval = G2G.ValidatedTxtCtrl(Texture,textureData['SH Coeff'][1],item,nDig=(8,3),OnLeave=OnODFValue)
                ODFSizer.Add(ODFval,0,WACV)
            mainSizer.Add(ODFSizer,0,WACV)
            mainSizer.Add((0,5),0)
        mainSizer.Add((0,5),0)
        mainSizer.Add(wx.StaticText(Texture,-1,' Sample orientation angle zeros: '),0,WACV)
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
        mainSizer.Add(angSizer,0,WACV|wx.LEFT,5)
#        mainSizer.Add(SHPenalty(textureData['Penalty']),0,WACV|wx.LEFT,5)  for future
        SetPhaseWindow(Texture,mainSizer)

################################################################################
##### DData routines - GUI stuff in GSASIIddataGUI.py
################################################################################
        
    def OnHklfAdd(event):
        keyList = data['Histograms'].keys()
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
                G2G.G2MessageBox(G2frame,'No reflections')
                return
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select reflection sets to use',
            'Use data',TextList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
            else:
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
                if G2frame.ErrorDialog('Likely error',msg,G2frame,wtype=wx.YES_NO) != wx.ID_YES: return

        wx.BeginBusyCursor()
        for i in result:
            histoName = TextList[i]
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
            refDict,reflData = G2frame.GPXtree.GetItemPyData(Id)
            data['Histograms'][histoName] = {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
                'Extinction':['Lorentzian','None',
                {'Tbar':0.1,'Cos2TM':0.955,'Eg':[1.e-7,False],'Es':[1.e-7,False],'Ep':[1.e-7,False]},],
                'Flack':[0.0,False],'Twins':[[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,0]],]}                        
            if 'TwMax' in reflData:     #nonmerohedral twins present
                data['Histograms'][histoName]['Twins'] = []
                for iT in range(reflData['TwMax'][0]+1):
                    if iT in reflData['TwMax'][1]:
                        data['Histograms'][histoName]['Twins'].append([False,0.0])
                    else:
                        data['Histograms'][histoName]['Twins'].append([np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,reflData['TwMax'][0]]])
            else:   #no nonmerohedral twins
                data['Histograms'][histoName]['Twins'] = [[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,0]],]
            UpdateHKLFdata(histoName)
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
        wx.EndBusyCursor()
        
    def OnDataUse(event):
#        hist = G2frame.hist
        if data['Histograms']:
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Use histograms', 
                'Use which histograms?',G2frame.dataWindow.HistsInPhase)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelections()
                    for id,item in enumerate(G2frame.dataWindow.HistsInPhase):
                        if id in sel:
                            data['Histograms'][item]['Use'] = True
                        else:
                            data['Histograms'][item]['Use'] = False                        
            finally:
                dlg.Destroy()
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)

    # Note: function replaced by one with identical name, below
    # def UpdateHKLFdata(histoName):
    #     generalData = data['General']
    #     Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
    #     refDict,reflData = G2frame.GPXtree.GetItemPyData(Id)
    #     SGData = generalData['SGData']
    #     Cell = generalData['Cell'][1:7]
    #     G,g = G2lat.cell2Gmat(Cell))
    #         try:
    #             if dlg.ShowModal() == wx.ID_OK:
    #                 sel = dlg.GetSelections()
    #                 for id,item in enumerate(keyList):  # <<<<< something is likely wrong here
    #                     if id in sel:
    #                         data['Histograms'][item]['Use'] = True
    #                     else:
    #                         data['Histograms'][item]['Use'] = False                        
    #         finally:
    #             dlg.Destroy()
    #     wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
                
    def UpdateHKLFdata(histoName):
        generalData = data['General']
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,histoName)
        refDict,reflData = G2frame.GPXtree.GetItemPyData(Id)
        SGData = generalData['SGData']
        Cell = generalData['Cell'][1:7]
        G,g = G2lat.cell2Gmat(Cell)
        for iref,ref in enumerate(reflData['RefList']):
            H = list(ref[:3])
            ref[4] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
            iabsnt,ref[3],Uniq,phi = G2spc.GenHKLf(H,SGData)
        
    def OnDataCopy(event):
        hist = G2frame.hist
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        sourceDict = copy.deepcopy(data['Histograms'][hist])
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Scale','Extinction','Babinet','Flack','Twins','Fix FXU']
        else:  #PWDR  
            copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','newLeBail','Fix FXU']
        copyDict = {}
        for name in copyNames: 
            copyDict[name] = sourceDict[name]        #force copy
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
            copyNames = ['Scale','Extinction','Babinet','Flack','Twins','Fix FXU']
        else:  #PWDR  
            copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','Fix FXU']
        babNames = ['BabA','BabU']
        for name in copyNames:
            if name in ['Scale','Extinction','HStrain','Flack','Twins']:
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
            elif name == 'Fix FXU':
                copyDict[name] = sourceDict[name]                      
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
                        if name in ['Scale','Extinction','HStrain','Flack','Twins']:
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
                        elif name == 'Fix FXU':
                            data['Histograms'][item][name] = copy.deepcopy(sourceDict[name])                      
        finally:
            dlg.Destroy()
        
    def OnSelDataCopy(event):
        hist = G2frame.hist
        sourceDict = data['Histograms'][hist]
        keyList = G2frame.dataWindow.HistsInPhase[:]
        if hist in keyList: keyList.remove(hist)
        if not keyList:
            G2G.G2MessageBox(G2frame,'No histograms to copy to')
            return
        if 'HKLF' in sourceDict['Histogram']:
            copyNames = ['Scale','Extinction','Babinet','Flack','Twins','Fix FXU']
        else:  #PWDR  
            copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction','Babinet','LeBail','Fix FXU']
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
            copyDict[parm] = copy.deepcopy(sourceDict[parm])
        dlg = G2G.G2MultiChoiceDialog(G2frame,u'Copy selected phase/histogram parameters\nfrom '+hist[5:][:35],
            'Copy selected phase/hist parameters', keyList)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    data['Histograms'][keyList[sel]].update(copy.deepcopy(copyDict))
        finally:
            dlg.Destroy()            
        
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
                        data['Histograms'][histoName] = {'Histogram':histoName,'Show':False,'LeBail':False,'newLeBail':True,
                            'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},['',],0.1],
                            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                                [1.,1.,1.,0.,0.,0.],6*[False,]],
                            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                            'HStrain':[NDij*[0.0,],NDij*[False,]],                          
                            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},'Fix FXU':' '}
                        refList = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Reflection Lists'))
                        refList[generalData['Name']] = {}                       
                    wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
            finally:
                dlg.Destroy()
                
    def OnDataDelete(event):
        if G2frame.dataWindow.HistsInPhase:
            DelList = []
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Delete histogram', 
                'Which histogram to delete from this phase?',G2frame.dataWindow.HistsInPhase)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    DelList = [G2frame.dataWindow.HistsInPhase[i] for i in dlg.GetSelections()]
                    for i in DelList:
                        del data['Histograms'][i]
            finally:
                dlg.Destroy()
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
        
    def OnDataApplyStrain(event):
        SGData = data['General']['SGData']        
        DijVals = data['Histograms'][G2frame.hist]['HStrain'][0][:]
        # apply the Dij values to the reciprocal cell
        newA = []
        Dijdict = dict(zip(G2spc.HStrainNames(SGData),DijVals))
        for Aij,lbl in zip(G2lat.cell2A(data['General']['Cell'][1:7]),
                            ['D11','D22','D33','D12','D13','D23']):
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

################################################################################
##### Rigid bodies
################################################################################

    def FillRigidBodyGrid(refresh=True):
        '''Fill the Rigid Body Phase information tab page.
        Note that the page is a ScrolledWindow, not a Grid
        '''
        def OnThermSel(event):       #needs to be seen by VecRbSizer!
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
            wx.CallAfter(FillRigidBodyGrid,True)
            if val != 'None':
                cia = data['General']['AtomPtrs'][3]
                for i,id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[id]][cia] = Ttype
            G2plt.PlotStructure(G2frame,data)
            
        def ThermDataSizer(RBObj,rbType):
            
            def OnThermval(invalid,value,tc):
                Cart = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[1]
                Uout = G2mth.UpdateRBUIJ(Bmat,Cart,RBObj)
                cia = data['General']['AtomPtrs'][3]
                for i,id in enumerate(RBObj['Ids']):
                    if Uout[i][0] == 'I':
                        data['Atoms'][AtLookUp[id]][cia+1] = Uout[i][1]
                    else:
                        data['Atoms'][AtLookUp[id]][cia+2:cia+8] = Uout[i][2:8]
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
                newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[0]
                for i,id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
                data['Drawing']['Atoms'] = []
                UpdateDrawAtoms(atomStyle)
                G2plt.PlotStructure(G2frame,data)
                
            def OnOrien(event):
                event.Skip()
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                A,V = G2mth.Q2AVdeg(RBObj['Orient'][0])
                V = np.inner(Bmat,V)
                try:
                    val = float(Obj.GetValue())
                    if item:
                        V[item-1] = val
                    else:
                        A = val
                    Obj.SetValue('%8.5f'%(val))
                    V = np.inner(Amat,V)
                    Q = G2mth.AVdeg2Q(A,V)
                    if not any(Q):
                        raise ValueError
                    RBObj['Orient'][0] = Q
                    newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[0]
                    for i,id in enumerate(RBObj['Ids']):
                        data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
                    data['Drawing']['Atoms'] = []
                    UpdateDrawAtoms(atomStyle)
                    G2plt.PlotStructure(G2frame,data)
                except ValueError:
                    pass
                
            topSizer = wx.FlexGridSizer(0,6,5,5)
            Orig = RBObj['Orig'][0]
            Orien,OrienV = G2mth.Q2AVdeg(RBObj['Orient'][0])
            Orien = [Orien,]
            Orien.extend(OrienV/nl.norm(OrienV))
            topSizer.Add(wx.StaticText(RigidBodies,-1,'Origin x,y,z:'),0,WACV)
            for ix,x in enumerate(Orig):
                origX = G2G.ValidatedTxtCtrl(RigidBodies,Orig,ix,nDig=(8,5),OnLeave=OnOrigX)
                topSizer.Add(origX,0,WACV)
            topSizer.Add((5,0),)
            Ocheck = wx.CheckBox(RigidBodies,-1,'Refine?')
            Ocheck.Bind(wx.EVT_CHECKBOX,OnOrigRef)
            Ocheck.SetValue(RBObj['Orig'][1])
            topSizer.Add(Ocheck,0,WACV)
            topSizer.Add(wx.StaticText(RigidBodies,-1,'Rotation angle, vector:'),0,WACV)
            for ix,x in enumerate(Orien):
                orien = wx.TextCtrl(RigidBodies,-1,value='%8.4f'%(x),style=wx.TE_PROCESS_ENTER)
                orien.Bind(wx.EVT_TEXT_ENTER,OnOrien)
                orien.Bind(wx.EVT_KILL_FOCUS,OnOrien)
                Indx[orien.GetId()] = ix
                topSizer.Add(orien,0,WACV)
            Qcheck = wx.ComboBox(RigidBodies,-1,value='',choices=[' ','A','AV'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Qcheck.Bind(wx.EVT_COMBOBOX,OnOrienRef)
            Qcheck.SetValue(RBObj['Orient'][1])
            topSizer.Add(Qcheck)
            return topSizer
                         
        def ResrbSizer(RBObj):
            
            def OnTorsionRef(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                RBObj['Torsions'][item][1] = Obj.GetValue()                
                
            def OnTorsion(invalid,value,tc):
                newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')[0]
                for i,id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
                data['Drawing']['Atoms'] = []
                UpdateDrawAtoms(atomStyle)
                drawAtoms.ClearSelection()
                G2plt.PlotStructure(G2frame,data)
                
            def OnDelResRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBObjs = data['RBModels']['Residue']
                for rbObj in RBObjs:
                    if RBId == rbObj['RBId']:
                        RBData['Residue'][RBId]['useCount'] -= 1
                        data['RBModels']['Residue'].remove(rbObj)                 
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True)
                
            resrbSizer = wx.BoxSizer(wx.VERTICAL)
            resrbSizer.Add(wx.StaticText(RigidBodies,-1,120*'-'))
            topLine = wx.BoxSizer(wx.HORIZONTAL)
            topLine.Add(wx.StaticText(RigidBodies,-1,
                'Name: '+RBObj['RBname']+RBObj['numChain']+'   '),0,WACV)
            rbId = RBObj['RBId']
            delRB = wx.CheckBox(RigidBodies,-1,'Delete?')
            delRB.Bind(wx.EVT_CHECKBOX,OnDelResRB)
            Indx[delRB.GetId()] = rbId
            topLine.Add(delRB,0,WACV)
            resrbSizer.Add(topLine)
            resrbSizer.Add(LocationSizer(RBObj,'Residue'))
            resrbSizer.Add(wx.StaticText(RigidBodies,-1,'Torsions:'),0,WACV)
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
            return resrbSizer
            
        def VecrbSizer(RBObj):
            G2frame.GetStatusBar().SetStatusText('NB: Rotation vector is in crystallographic space',1)
                   
            def OnDelVecRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBData['Vector'][RBId]['useCount'] -= 1                
                RBObjs = data['RBModels']['Vector']
                for rbObj in RBObjs:
                    if RBId == rbObj['RBId']:
                       data['RBModels']['Vector'].remove(rbObj)                 
                G2plt.PlotStructure(G2frame,data)
                wx.CallAfter(FillRigidBodyGrid,True)
             
            vecrbSizer = wx.BoxSizer(wx.VERTICAL)
            vecrbSizer.Add(wx.StaticText(RigidBodies,-1,120*'-'))
            topLine = wx.BoxSizer(wx.HORIZONTAL)
            topLine.Add(wx.StaticText(RigidBodies,-1,
                'Name: '+RBObj['RBname']+'   '),0,WACV)
            rbId = RBObj['RBId']
            delRB = wx.CheckBox(RigidBodies,-1,'Delete?')
            delRB.Bind(wx.EVT_CHECKBOX,OnDelVecRB)
            Indx[delRB.GetId()] = rbId
            topLine.Add(delRB,0,WACV)
            vecrbSizer.Add(topLine)
            vecrbSizer.Add(LocationSizer(RBObj,'Vector'))
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
        
        def OnSelect(event):
            rbId = select.GetSelection()
            wx.CallLater(100,RepaintRBInfo,rbId)
           
        def RepaintRBInfo(rbId,Scroll=0):
            oldFocus = wx.Window.FindFocus()
            if 'phoenix' in wx.version():
                G2frame.bottomSizer.Clear(True)
            else:
                G2frame.bottomSizer.DeleteWindows()
            Indx.clear()
            rbObj = data['RBModels']['Residue'][rbId]
            data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
            Quad = rbObj['Orient'][0]
            data['Drawing']['Quaternion'] = G2mth.invQ(Quad)
            G2frame.bottomSizer =  ResrbSizer(rbObj)
            mainSizer.Add(G2frame.bottomSizer)
            mainSizer.Layout()
            G2frame.dataWindow.Refresh()
            RigidBodies.SetVirtualSize(mainSizer.GetMinSize())
            RigidBodies.Scroll(0,Scroll)
            G2frame.dataWindow.SendSizeEvent()
            G2plt.PlotStructure(G2frame,data)
            wx.CallAfter(oldFocus.SetFocus)
        
        # FillRigidBodyGrid executable code starts here
        if refresh:
            if RigidBodies.GetSizer(): RigidBodies.GetSizer().Clear(True)
        general = data['General']
        cx,ct,cs,cia = general['AtomPtrs']
        AtLookUp = G2mth.FillAtomLookUp(data['Atoms'],cia+8)
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies')
        if not Id:
            return
        RBData = G2frame.GPXtree.GetItemPyData(Id)
        Indx = {}
        atomStyle = 'balls & sticks'
        if 'macro' in general['Type']:
            atomStyle = 'sticks'
        G2frame.GetStatusBar().SetStatusText('',1)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if not data['RBModels']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'No rigid body models:'),0,WACV)
            mainSizer.Add((5,5),0)
        if 'Residue' in data['RBModels'] and len(data['RBModels']['Residue']):
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'Residue rigid bodies:'),0,WACV)
            mainSizer.Add((5,5),0)
            RBnames = []
            for RBObj in data['RBModels']['Residue']:
                RBnames.append(RBObj['RBname'].split(':')[0]+RBObj['numChain'])
            rbName = RBnames[0]
            rbObj = data['RBModels']['Residue'][0]
            data['Drawing']['viewPoint'][0] = rbObj['Orig'][0]
            data['Drawing']['Quaternion'] = rbObj['Orient'][0]
            select = wx.ListBox(RigidBodies,choices=RBnames,style=wx.LB_SINGLE,size=(-1,120))
            select.SetSelection(RBnames.index(rbName))
            select.SetFirstItem(RBnames.index(rbName))
            select.Bind(wx.EVT_LISTBOX,OnSelect)
            mainSizer.Add(select,0,WACV)
            G2frame.bottomSizer = wx.BoxSizer(wx.VERTICAL)
            G2frame.bottomSizer.Add(ResrbSizer(rbObj))
            mainSizer.Add(G2frame.bottomSizer)
            G2frame.GetStatusBar().SetStatusText('NB: Rotation vector is in crystallographic space',1)
            G2plt.PlotStructure(G2frame,data)
        if 'Vector' in data['RBModels'] and len(data['RBModels']['Vector']):
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'Vector rigid bodies:'),0,WACV)
            mainSizer.Add((5,5),0)
            for RBObj in data['RBModels']['Vector']:
                mainSizer.Add(VecrbSizer(RBObj))

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
            for item in ['Orig','Orient','ThermalMotion']: 
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
                
    def OnRBAssign(event):
        
        G2frame.GetStatusBar().SetStatusText('',1)
        RBData = G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbNames = {}
        for rbVec in RBData['Vector']:
            if rbVec != 'AtInfo':
                rbNames[RBData['Vector'][rbVec]['RBname']] =['Vector',rbVec]
        for rbRes in RBData['Residue']:
            if rbRes != 'AtInfo':
                rbNames[RBData['Residue'][rbRes]['RBname']] = ['Residue',rbRes]
        if not rbNames:
            print ('**** ERROR - no rigid bodies defined ****')
            return
        general = data['General']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        cx,ct = general['AtomPtrs'][:2]
        atomData = data['Atoms']
        Indx = {}
        atInd = [-1,-1,-1]
        data['testRBObj'] = {}
            
        def Draw():
            
            def OnOk(event):
                rbType = data['testRBObj']['rbType']
                RBObjs = data['RBModels'].get(rbType,[])
                rbObj = data['testRBObj']['rbObj']
                rbId = rbObj['RBId']
                newXYZ = G2mth.UpdateRBXYZ(Bmat,rbObj,RBData,rbType)[0]
                Ids = []
                dmax = 0.0
                oldXYZ = G2mth.getAtomXYZ(atomData,cx)
                for xyz in newXYZ:
                    dist = G2mth.GetXYZDist(xyz,oldXYZ,Amat)
                    dmax = max(dmax,np.min(dist))
                    id = np.argmin(dist)
                    Id = atomData[id][-1]
                    if Id in Ids:   #duplicate - 2 atoms on same site; invalidate & look again
                        dist[id] = 100.
                        id =  np.argmin(dist)
                        Id = atomData[id][-1]
                    Ids.append(Id)
                    atomData[id][cx:cx+3] = xyz
                if dmax > 1.0:
                    print ('**** WARNING - some atoms not found or misidentified ****')
                    print ('****           check torsion angles & try again      ****')
                    OkBtn.SetLabel('Not Ready')
                    OkBtn.Enable(False)
                    return
                rbObj['Ids'] = Ids
                rbObj['ThermalMotion'] = ['None',[0. for i in range(21)],[False for i in range(21)]] #type,values,flags
                rbObj['RBname'] += ':'+str(RBData[rbType][rbId]['useCount'])
                RBObjs.append(rbObj)
                data['RBModels'][rbType] = RBObjs
                RBData[rbType][rbId]['useCount'] += 1
                del data['testRBObj']
                G2plt.PlotStructure(G2frame,data)
                FillRigidBodyGrid(True)
                
            def OnCancel(event):
                del data['testRBObj']
                FillRigidBodyGrid(True)
                
            def OnRBSel(event):
                selection = rbSel.GetValue()
                rbType,rbId = rbNames[selection]
                data['testRBObj']['rbAtTypes'] = RBData[rbType][rbId]['rbTypes']
                data['testRBObj']['AtInfo'] = RBData[rbType]['AtInfo']
                data['testRBObj']['rbType'] = rbType
                data['testRBObj']['rbData'] = RBData
                data['testRBObj']['Sizers'] = {}
                rbRef = RBData[rbType][rbId]['rbRef']
                data['testRBObj']['rbRef'] = rbRef
                refType = []
                refName = []
                for ref in rbRef[:3]:
                    reftype = data['testRBObj']['rbAtTypes'][ref]
                    refType.append(reftype)
                    refName.append(reftype+' '+str(rbRef[0]))
                atNames,AtNames = fillAtNames(refType,atomData,ct)
                data['testRBObj']['atNames'] = atNames
                data['testRBObj']['AtNames'] = AtNames
                data['testRBObj']['rbObj'] = {'Orig':[[0,0,0],False],
                    'Orient':[[0.,0.,0.,0.],' '],'Ids':[],'RBId':rbId,'Torsions':[],
                    'numChain':'','RBname':RBData[rbType][rbId]['RBname']}
                data['testRBObj']['torAtms'] = []                
                for item in RBData[rbType][rbId].get('rbSeq',[]):
                    data['testRBObj']['rbObj']['Torsions'].append([item[2],False])
                    data['testRBObj']['torAtms'].append([-1,-1,-1])
                wx.CallAfter(Draw)
                
            def fillAtNames(refType,atomData,ct):
                atNames = [{},{},{}]
                AtNames = {}
                for iatm,atom in enumerate(atomData):
                    AtNames[atom[ct-1]] = iatm
                    for i,reftype in enumerate(refType):
                        if atom[ct] == reftype:
                            atNames[i][atom[ct-1]] = iatm
                return atNames,AtNames
                
            def OnAtOrigPick(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                atName = Obj.GetValue()
                rbType = data['testRBObj']['rbType']
                atInd[0] = atNames[item][atName]
                if 'Vector' in rbType:
                    rbObj = data['testRBObj']['rbObj']
                    rbId = rbObj['RBId']
                    rbRef = data['testRBObj']['rbRef']
                    rbXYZ = -RBData[rbType][rbId]['rbXYZ']
                    nref = atNames[item][atName]
                    Oxyz = np.inner(Bmat,np.array(rbXYZ[rbRef[0]]))
                    Nxyz = np.array(atomData[nref][cx:cx+3])
                    Orig = Nxyz-Oxyz
                    data['testRBObj']['rbObj']['Orig'][0] = Orig   
                else:
                    Orig = atomData[atNames[item][atName]][cx:cx+3]
                    data['testRBObj']['rbObj']['Orig'][0] = Orig
                for x,item in zip(Orig,Xsizers):
                    item.SetLabel('%10.5f'%(x))
                G2plt.PlotStructure(G2frame,data)
                
            def OnAtQPick(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                atName = Obj.GetValue()
                atInd[item] = atNames[item][atName]
                if any([x<0 for x in atInd]):
                    return
                OkBtn.SetLabel('OK')
                OkBtn.Enable(True)
                rbType = data['testRBObj']['rbType']
                rbObj = data['testRBObj']['rbObj']
                rbId = rbObj['RBId']
                rbRef = data['testRBObj']['rbRef']
                rbXYZ = RBData[rbType][rbId]['rbXYZ']
                rbOrig = rbXYZ[rbRef[0]]
                VAR = rbXYZ[rbRef[1]]-rbOrig
                VBR = rbXYZ[rbRef[2]]-rbOrig
                if rbType == 'Vector':
                    Orig = np.array(atomData[atInd[0]][cx:cx+3])
                else:
                    Orig = np.array(data['testRBObj']['rbObj']['Orig'][0])                
                VAC = np.inner(Amat,np.array(atomData[atInd[1]][cx:cx+3]-Orig))
                VBC = np.inner(Amat,np.array(atomData[atInd[2]][cx:cx+3]-Orig))
                VCC = np.cross(VAR,VAC)
                if nl.norm(VCC) > 1e-7:
                    QuatA = G2mth.makeQuat(VAR,VAC,VCC)[0]
                    VAR = G2mth.prodQVQ(QuatA,VAR)
                    VBR = G2mth.prodQVQ(QuatA,VBR)
                    QuatB = G2mth.makeQuat(VBR,VBC,VAR)[0]
                    QuatC = G2mth.prodQQ(QuatB,QuatA)
                else:                               #parallel/antiparallel
                    if np.dot(VAR,VAC)/(nl.norm(VAR)*nl.norm(VAC)) > 1e-7:  #no rotation
                        QuatC = G2mth.AVdeg2Q(0.,[0.,0.,1.])
                    else:
                        QuatC = G2mth.AVdeg2Q(180.,VBR+VBC)
                data['testRBObj']['rbObj']['Orient'] = [QuatC,' ']
                for x,item in zip(QuatC,Osizers):
                    item.SetLabel('%10.4f'%(x))                
                if rbType == 'Vector':
                    Oxyz = np.inner(Bmat,G2mth.prodQVQ(QuatC,rbOrig))
                    Nxyz = np.array(atomData[atInd[0]][cx:cx+3])
                    Orig = Nxyz-Oxyz
                    data['testRBObj']['rbObj']['Orig'][0] = Orig
                    for x,item in zip(Orig,Xsizers):
                        item.SetLabel('%10.5f'%(x))
                G2plt.PlotStructure(G2frame,data)
                
            def OnTorAngle(invalid,value,tc):
                OkBtn.SetLabel('OK')
                OkBtn.Enable(True)
                Obj = tc.event.GetEventObject()
                [tor,torSlide] = Indx[Obj.GetId()]
                torSlide.SetValue(int(value*10))
                G2plt.PlotStructure(G2frame,data)
                
            def OnTorSlide(event):
                OkBtn.SetLabel('OK')
                OkBtn.Enable(True)
                Obj = event.GetEventObject()
                tor,ang = Indx[Obj.GetId()]
                Tors = data['testRBObj']['rbObj']['Torsions'][tor]
                val = float(Obj.GetValue())/10.
                Tors[0] = val
                ang.SetValue(val)
                G2plt.PlotStructure(G2frame,data)

            if len(data['testRBObj']):
                G2plt.PlotStructure(G2frame,data)
                    
            if RigidBodies.GetSizer(): RigidBodies.GetSizer().Clear(True)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add((5,5),0)
            if data['testRBObj']:
                Xsizers = []
                Osizers = []
                rbObj = data['testRBObj']['rbObj']
                rbName = rbObj['RBname']
                rbId = rbObj['RBId']
                Orig = rbObj['Orig'][0]
                Orien = rbObj['Orient'][0]
                rbRef = data['testRBObj']['rbRef']
                Torsions = rbObj['Torsions']
                refName = []
                for ref in rbRef:
                    refName.append(data['testRBObj']['rbAtTypes'][ref]+str(ref))
                atNames = data['testRBObj']['atNames']
                mainSizer.Add(wx.StaticText(RigidBodies,-1,'Locate rigid body : '+rbName),
                    0,WACV)
                mainSizer.Add((5,5),0)
                OriSizer = wx.FlexGridSizer(0,5,5,5)
                OriSizer.Add(wx.StaticText(RigidBodies,-1,'Origin x,y,z: '),0,WACV)
                for ix,x in enumerate(Orig):
                    origX = wx.StaticText(RigidBodies,-1,'%10.5f'%(x))
                    OriSizer.Add(origX,0,WACV)
                    Xsizers.append(origX)
                OriSizer.Add((5,0),)
                if len(atomData):
                    choice = list(atNames[0].keys())
                    choice.sort()
                    data['testRBObj']['Sizers']['Xsizers'] = Xsizers
                OriSizer.Add(wx.StaticText(RigidBodies,-1,'Orientation quaternion: '),0,WACV)
                for ix,x in enumerate(Orien):
                    orien = wx.StaticText(RigidBodies,-1,'%10.4f'%(x))
                    OriSizer.Add(orien,0,WACV)
                    Osizers.append(orien)
                data['testRBObj']['Sizers']['Osizers'] = Osizers
                mainSizer.Add(OriSizer)
                mainSizer.Add((5,5),0)
                RefSizer = wx.FlexGridSizer(0,7,5,5)
                if len(atomData):
                    RefSizer.Add(wx.StaticText(RigidBodies,-1,'Location setting: Select match to'),0,WACV)
                    for i in [0,1,2]:
                        choice = ['',]+list(atNames[i].keys())
                        choice.sort()
                        RefSizer.Add(wx.StaticText(RigidBodies,-1,' '+refName[i]+': '),0,WACV)
                        atPick = wx.ComboBox(RigidBodies,-1,value='',
                            choices=choice[1:],style=wx.CB_READONLY|wx.CB_DROPDOWN)
                        if i:
                            atPick.Bind(wx.EVT_COMBOBOX, OnAtQPick)
                        else:
                            atPick.Bind(wx.EVT_COMBOBOX, OnAtOrigPick)                            
                        Indx[atPick.GetId()] = i
                        RefSizer.Add(atPick,0,WACV)
                mainSizer.Add(RefSizer)
                mainSizer.Add((5,5),0)
                if Torsions:                    
                    rbSeq = RBData['Residue'][rbId]['rbSeq']
                    TorSizer = wx.FlexGridSizer(0,4)
                    TorSizer.AddGrowableCol(1,1)
                    for t,[torsion,seq] in enumerate(zip(Torsions,rbSeq)):
                        torName = ''
                        for item in [seq[0],seq[1],seq[3][0]]:
                            torName += data['testRBObj']['rbAtTypes'][item]+str(item)+' '
                        TorSizer.Add(wx.StaticText(RigidBodies,-1,'Side chain torsion for rb seq: '+torName),0,WACV)
                        torSlide = wx.Slider(RigidBodies,style=wx.SL_HORIZONTAL)
                        torSlide.SetRange(0,3600)
                        torSlide.SetValue(int(torsion[0]*10.))
                        torSlide.Bind(wx.EVT_SLIDER, OnTorSlide)
                        TorSizer.Add(torSlide,1,wx.EXPAND|wx.RIGHT)
                        TorSizer.Add(wx.StaticText(RigidBodies,-1,' Angle: '),0,WACV)
                        ang = G2G.ValidatedTxtCtrl(RigidBodies,torsion,0,nDig=(8,3),typeHint=float,OnLeave=OnTorAngle)
                        Indx[torSlide.GetId()] = [t,ang]
                        Indx[ang.GetId()] = [t,torSlide]
                        TorSizer.Add(ang,0,WACV)                            
                    mainSizer.Add(TorSizer,1,wx.EXPAND|wx.RIGHT)
                else:
                    mainSizer.Add(wx.StaticText(RigidBodies,-1,'No side chain torsions'),0,WACV)
            else:
                mainSizer.Add(wx.StaticText(RigidBodies,-1,'Assign rigid body:'),0,WACV)
                mainSizer.Add((5,5),0)
                topSizer = wx.BoxSizer(wx.HORIZONTAL)
                topSizer.Add(wx.StaticText(RigidBodies,-1,'Select rigid body model'),0,WACV)
                rbSel = wx.ComboBox(RigidBodies,-1,value='',choices=list(rbNames.keys()),
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                rbSel.Bind(wx.EVT_COMBOBOX, OnRBSel)
                topSizer.Add((5,5),0)
                topSizer.Add(rbSel,0,WACV)
                mainSizer.Add(topSizer)                
                
            OkBtn = wx.Button(RigidBodies,-1,"Not ready")
            OkBtn.Bind(wx.EVT_BUTTON, OnOk)
            OkBtn.Enable(False)
            CancelBtn = wx.Button(RigidBodies,-1,'Cancel')
            CancelBtn.Bind(wx.EVT_BUTTON, OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add(CancelBtn)
            btnSizer.Add((20,20),1)
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            SetPhaseWindow(RigidBodies,mainSizer)
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
                for i,id in enumerate(RBObj['Ids']):
                    data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
        finally:
            wx.EndBusyCursor()
        wx.CallAfter(FillRigidBodyGrid,True)
        
    def OnRBRemoveAll(event):
        data['RBModels']['Residue'] = []
        data['RBModels']['Vector'] = []
        RBData = G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        for RBType in ['Vector','Residue']:
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
                    for i,id in enumerate(rbObj['Ids']):
                        data['Atoms'][AtLookUp[id]][cia] = Ttype
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
            
################################################################################
##### MC/SA routines
################################################################################

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
            rbsizer1.Add(wx.StaticText(G2frame.MCSA,-1,model['Type']+': '+model['name']+': '),0,WACV)
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
                posRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[posRange.GetId()] = [model,'Pos',ix]
                posRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                posRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                rbsizer1.Add(posRange,0,WACV)
                
            rbsizer2 = wx.FlexGridSizer(0,6,5,5)
            Ori = model['Ori'][0]
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oa: '),0,WACV)
            angVal = wx.TextCtrl(G2frame.MCSA,-1,'%.5f'%(Ori[0]),style=wx.TE_PROCESS_ENTER)
            angVal.Bind(wx.EVT_TEXT_ENTER,OnOriVal)
            angVal.Bind(wx.EVT_KILL_FOCUS,OnOriVal)
            rbsizer2.Add(angVal,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
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
            angRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(Rge[0][0],Rge[0][1]),style=wx.TE_PROCESS_ENTER)
            Indx[angRange.GetId()] = [model,'Ori',0]
            angRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
            angRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
            rbsizer2.Add(angRange,0,WACV)
            rbsizer2.Add(wx.StaticText(G2frame.MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
            for io,item in enumerate(['Oi','Oj','Ok']):
                rmin,rmax = Rge[io+1]
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
                    rbsizer.Add(wx.StaticText(G2frame.MCSA,-1,'Torsions:'),0,WACV)
                    rbsizer3 = wx.FlexGridSizer(0,8,5,5)
                    for it,tor in enumerate(model['Tor'][0]):
                        iBeg,iFin = RBData['Residue'][model['RBId']]['rbSeq'][it][:2]
                        name = atNames[iBeg]+'-'+atNames[iFin]
                        torRef = wx.CheckBox(G2frame.MCSA,-1,label=' %s: '%(name))
                        torRef.SetValue(model['Tor'][1][it])
                        torRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                        Indx[torRef.GetId()] = [model,'Tor',it]
                        rbsizer3.Add(torRef,0,WACV)
                        torVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,model['Tor'][0],it,nDig=(10,4),OnLeave=OnPosVal)
                        rbsizer3.Add(torVal,0,WACV)
                        rbsizer3.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
                        rmin,rmax = model['Tor'][2][it]
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
            poVal = G2G.ValidatedTxtCtrl(G2frame.MCSA,POData['Coef'],0,nDig=(10,3),min=0.)
            poSizer.Add(poVal,0,WACV)
            poSizer.Add(wx.StaticText(G2frame.MCSA,-1,' Range: '),0,WACV)
            rmin,rmax = POData['Coef'][2]
            poRange = wx.TextCtrl(G2frame.MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
            poRange.Bind(wx.EVT_TEXT_ENTER,OnPORange)
            poRange.Bind(wx.EVT_KILL_FOCUS,OnPORange)
            poSizer.Add(poRange,0,WACV)                       
            poSizer.Add(wx.StaticText(G2frame.MCSA,-1,' Unique axis, H K L: '),0,WACV)
            h,k,l = POData['axis']
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
        
        # UpdateMCSA executable code starts here
        if G2frame.MCSA.GetSizer(): G2frame.MCSA.GetSizer().Clear(True)
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
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if not data['MCSA']['Models']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,-1,'No MC/SA models:'),0,WACV)
            mainSizer.Add((5,5),0)
        else:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,-1,'MC/SA models:'),0,WACV)
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
                mainSizer.Add(select,0,WACV)
                G2frame.bottomSizer = wx.BoxSizer(wx.VERTICAL)
                G2frame.bottomSizer.Add(rbSizer(data['MCSA']['Models'][rbids[0]]))
                mainSizer.Add(G2frame.bottomSizer)
                
        if not data['MCSA']['Results']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,-1,'No MC/SA results:'),0,WACV)
            mainSizer.Add((5,5),0)
        else:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(G2frame.MCSA,-1,'MC/SA results:'),0,WACV)
            mainSizer.Add((5,5),0)
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
            pgbar = wx.ProgressDialog('MC/SA','Residual Rcf =',101.0, 
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
        
################################################################################
##### Pawley routines
################################################################################

    def FillPawleyReflectionsGrid():
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
        generalData = data['General']
        if 'Pawley ref' in data:
            PawleyPeaks = data['Pawley ref']                        
            rowLabels = []
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
            G2frame.PawleyRefl.SetTable(PawleyTable, True)
            G2frame.PawleyRefl.Bind(wx.EVT_KEY_DOWN, KeyEditPawleyGrid)                 
            for r in range(G2frame.PawleyRefl.GetNumberRows()):
                for c in range(G2frame.PawleyRefl.GetNumberCols()):
                    if c in pos:
                        G2frame.PawleyRefl.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        G2frame.PawleyRefl.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.PawleyRefl.SetMargins(0,0)
            G2frame.PawleyRefl.AutoSizeColumns(False)
            SetPhaseWindow(G2frame.PawleyRefl)
                    
    def OnPawleySet(event):
        '''Set Pawley parameters and optionally recompute
        '''
        #GSASIIpath.IPyBreak()
        
        def DisablePawleyOpts(*args):
            pawlVal.Enable(generalData['doPawley'])
            pawlNegWt.Enable(generalData['doPawley'])
        generalData = data['General']
        startDmin = generalData['Pawley dmin']
        genDlg = wx.Dialog(G2frame,wx.ID_ANY,'Set Pawley Parameters',
                        style=wx.DEFAULT_DIALOG_STYLE)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(genDlg,wx.ID_ANY,
                                    'Set Pawley Extraction Parameters for phase '+
                                    generalData.get('Name','?')))
        mainSizer.Add([5,10])
        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label=' Do Pawley refinement?: '),0,WACV)
        pawlRef = G2G.G2CheckBox(genDlg,'',generalData,'doPawley',
                             DisablePawleyOpts)
        pawleySizer.Add(pawlRef,0,WACV)
        mainSizer.Add(pawleySizer)
        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label=' Pawley dmin: '),0,WACV)
        def d2Q(*a,**kw):
            temp['Qmax'] = 2 * math.pi / generalData['Pawley dmin']
            pawlQVal.SetValue(temp['Qmax'])
        pawlVal = G2G.ValidatedTxtCtrl(genDlg,generalData,'Pawley dmin',
               min=0.25,max=20.,nDig=(10,5),typeHint=float,OnLeave=d2Q)
        pawleySizer.Add(pawlVal,0,WACV)
        pawleySizer.Add(wx.StaticText(genDlg,label='   Qmax: '),0,WACV)
        temp = {'Qmax':2 * math.pi / generalData['Pawley dmin']}
        def Q2D(*args,**kw):
            generalData['Pawley dmin'] = 2 * math.pi / temp['Qmax']
            pawlVal.SetValue(generalData['Pawley dmin'])        
        pawlQVal = G2G.ValidatedTxtCtrl(genDlg,temp,'Qmax',
               min=0.314,max=25.,nDig=(10,5),typeHint=float,OnLeave=Q2D)
        pawleySizer.Add(pawlQVal,0,WACV)
        mainSizer.Add(pawleySizer)
        pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
        pawleySizer.Add(wx.StaticText(genDlg,label=' Pawley neg. wt.: '),0,WACV)
        pawlNegWt = G2G.ValidatedTxtCtrl(genDlg,generalData,'Pawley neg wt',
                    min=0.,max=1.,nDig=(10,4),typeHint=float)
        pawleySizer.Add(pawlNegWt,0,WACV)
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

        if generalData['doPawley'] and res == wx.ID_OK and startDmin != generalData['Pawley dmin']:
            dlg = wx.MessageDialog(G2frame,'Do you want to initialize the Pawley reflections with the new Dmin value?','Initialize Pawley?', 
                wx.YES_NO | wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_NO:
                    return
            finally:
                dlg.Destroy()
            OnPawleyLoad(event)
            
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
        FillPawleyReflectionsGrid()
        
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
        HistoNames = filter(lambda a:Histograms[a]['Use']==True,list(Histograms.keys()))
        if not len(HistoNames):
            G2frame.ErrorDialog('Pawley estimate','No histograms defined for this phase')
            return
        PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,HistoNames[0])       #only use 1st histogram
        xdata = G2frame.GPXtree.GetItemPyData(PatternId)[1]
        Inst = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Instrument Parameters'))[0]
        Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,PatternId,'Sample Parameters'))
        wave = G2mth.getWave(Inst)
        const = 9.e-2/(np.pi*Sample['Gonio. radius'])                  #shifts in microns
        gconst = 2.35482 # sqrt(8 ln 2)
        
        wx.BeginBusyCursor()
        try:
            for ref in Refs:
                pos = 2.0*asind(wave/(2.0*ref[4+im]))
                if 'Bragg' in Sample['Type']:
                    pos -= const*(4.*Sample['Shift'][0]*cosd(pos/2.0)+ \
                        Sample['Transparency'][0]*sind(pos)*100.0)            #trans(=1/mueff) in cm
                else:               #Debye-Scherrer - simple but maybe not right
                    pos -= const*(Sample['DisplaceX'][0]*cosd(pos)+Sample['DisplaceY'][0]*sind(pos))
                indx = np.searchsorted(xdata[0],pos)
                try:
                    FWHM = max(0.001,G2pwd.getFWHM(pos,Inst))
                    # We want to estimate Pawley F^2 as a drop-in replacement for F^2 calculated by the structural 
                    # routines, which use Icorr * F^2 * peak profile, where peak profile has an area of 1.  So
                    # we multiply the observed peak height by sqrt(8 ln 2)/(FWHM*sqrt(pi)) to determine the value of Icorr*F^2 
                    # then divide by Icorr to get F^2.
                    ref[6+im] = (xdata[1][indx]-xdata[4][indx])*gconst/(FWHM*np.sqrt(np.pi))  #Area of Gaussian is height * FWHM * sqrt(pi)
                    Lorenz = 1./(2.*sind(xdata[0][indx]/2.)**2*cosd(xdata[0][indx]/2.))           #Lorentz correction
                    pola = 1.0
                    if 'X' in Inst['Type']:
                        pola,dIdPola = G2pwd.Polarization(Inst['Polariz.'][1],xdata[0][indx],Inst['Azimuth'][1])
                    else:
                        pola = 1.0
                    # Include histo scale and volume in calculation
                    ref[6+im] /= (Sample['Scale'][0] * Vst * Lorenz * pola * ref[3+im])
                except IndexError:
                    pass
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
        raise Exception        

        refcol = [G2frame.PawleyRefl.GetColLabelValue(c) for c in range(G2frame.PawleyRefl.GetNumberCols())].index('refine')
        for r in range(G2frame.PawleyRefl.GetNumberRows()):
            G2frame.PawleyRefl.GetTable().SetValue(
                r,refcol,
                not G2frame.PawleyRefl.GetTable().GetValueAsBool(r,refcol))
        G2frame.PawleyRefl.ForceRefresh()
                            
################################################################################
##### Fourier routines
################################################################################

    def FillMapPeaksGrid():
                        
        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                if MapPeaks.IsSelection():
                    MapPeaks.ClearSelection()
                else:
                    for row in range(MapPeaks.GetNumberRows()):
                        MapPeaks.SelectRow(row,True)
                    
            elif c < 0:                   #only row clicks
                if event.ControlDown():                    
                    if r in MapPeaks.GetSelectedRows():
                        MapPeaks.DeselectRow(r)
                    else:
                        MapPeaks.SelectRow(r,True)
                elif event.ShiftDown():
                    indxs = MapPeaks.GetSelectedRows()
                    MapPeaks.ClearSelection()
                    ibeg = 0
                    if indxs:
                        ibeg = indxs[-1]
                    for row in range(ibeg,r+1):
                        MapPeaks.SelectRow(row,True)
                else:
                    MapPeaks.ClearSelection()
                    MapPeaks.SelectRow(r,True)
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
            
        G2frame.GetStatusBar().SetStatusText('',1)
        if 'Map Peaks' in data:
            G2frame.GetStatusBar().SetStatusText('Double click any column heading to sort',1)
            mapPeaks = data['Map Peaks']                        
            rowLabels = []
            for i in range(len(mapPeaks)): rowLabels.append(str(i))
            colLabels = ['mag','x','y','z','dzero','dcent']
            Types = 6*[wg.GRID_VALUE_FLOAT+':10,4',]
            G2frame.MapPeaksTable = G2G.Table(mapPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            MapPeaks.SetTable(G2frame.MapPeaksTable, True)
            MapPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
            for r in range(MapPeaks.GetNumberRows()):
                for c in range(MapPeaks.GetNumberCols()):
                    MapPeaks.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            MapPeaks.SetMargins(0,0)
            MapPeaks.AutoSizeColumns(False)
            SetPhaseWindow(MapPeaks)
                    
    def OnPeaksMove(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            peakMax = np.amax(mapPeaks.T[0])
            Ind = MapPeaks.GetSelectedRows()
            for ind in Ind:
                mag,x,y,z = mapPeaks[ind][:4]
                AtomAdd(x,y,z,'H',Name='M '+'%d'%(int(100*mag/peakMax)))
            G2plt.PlotStructure(G2frame,data)
    
    def OnPeaksClear(event):
        data['Map Peaks'] = []
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)
        
    def OnPeaksDelete(event):
        if 'Map Peaks' in data:
            mapPeaks = data['Map Peaks']
            Ind = MapPeaks.GetSelectedRows()
            Ind.sort()
            Ind.reverse()
            for ind in Ind:
                mapPeaks = np.delete(mapPeaks,ind,0)
            data['Map Peaks'] = mapPeaks
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)
        
    def OnPeaksEquiv(event):
        if 'Map Peaks' in data:
            Ind = MapPeaks.GetSelectedRows()
            if Ind:
                wx.BeginBusyCursor()
                try:
                    Ind = G2mth.PeaksEquiv(data,Ind)
                    for r in range(MapPeaks.GetNumberRows()):
                        if r in Ind:
                            MapPeaks.SelectRow(r,addToSelected=True)
                        else:
                            MapPeaks.DeselectRow(r)
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
            mapPeaks = data['Map Peaks']
            Ind = MapPeaks.GetSelectedRows()
            if Ind:
                wx.BeginBusyCursor()
                try:
                    Ind = G2mth.PeaksUnique(data,Ind)
                    print (' No. unique peaks: %d Unique peak fraction: %.3f'%(len(Ind),float(len(Ind))/len(mapPeaks)))
                    for r in range(MapPeaks.GetNumberRows()):
                        if r in Ind:
                            MapPeaks.SelectRow(r,addToSelected=True)
                        else:
                            MapPeaks.DeselectRow(r)
                finally:
                    wx.EndBusyCursor()
                G2plt.PlotStructure(G2frame,data)
                
    def OnPeaksViewPoint(event):
        # set view point
        indx = MapPeaks.GetSelectedRows()
        if not indx:
            G2frame.ErrorDialog('Set viewpoint','No peaks selected')
            return
        mapPeaks = data['Map Peaks']
        drawingData = data['Drawing']
        drawingData['viewPoint'][0] = mapPeaks[indx[0]][1:4]
        G2plt.PlotStructure(G2frame,data)
    
    def OnPeaksDistVP(event):
        # distance to view point
        indx = MapPeaks.GetSelectedRows()
        if not indx:
            G2frame.ErrorDialog('Peak distance','No peaks selected')
            return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])            
        mapPeaks = data['Map Peaks']
        drawingData = data['Drawing']
        viewPt = np.array(drawingData['viewPoint'][0])
        print (' Distance from view point at %.3f %.3f %.3f to:'%(viewPt[0],viewPt[1],viewPt[2]))
        colLabels = [MapPeaks.GetColLabelValue(c) for c in range(MapPeaks.GetNumberCols())]
        cx = colLabels.index('x')
        cm = colLabels.index('mag')
        for i in indx:
            peak = mapPeaks[i]
            Dx = np.array(peak[cx:cx+3])-viewPt
            dist = np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0))
            print ('Peak: %5d mag= %8.2f distance = %.3f'%(i,peak[cm],dist))

    def OnPeaksDA(event):
        #distance, angle 
        indx = MapPeaks.GetSelectedRows()
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
        print (dim+mapData['MapType']+' computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig))
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)
        
    def OnFourClear(event):
        generalData = data['General']
        generalData['Map'] = mapDefault.copy()
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
                if 'N' in mapData['Type']:      #look for negatives in neutron maps
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
            UpdateDrawAtoms()
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
        pgbar = wx.ProgressDialog('Charge flipping','Residual Rcf =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            newMap,new4Dmap = G2mth.SSChargeFlip(data,ReflData,pgbar)
        finally:
            pgbar.Destroy()
        mapData.update(newMap)
        map4DData.update(new4Dmap)
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
        pgbar = wx.ProgressDialog('Charge flipping','Residual Rcf =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        testNames = ['%3d%3d%3d'%(h,k,l) for h,k,l in flipData['testHKL']]
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            pgbar.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        try:
            result = G2mth.ChargeFlip(data,ReflData,pgbar)
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
        pgbar = wx.ProgressDialog('Texture fit','Residual = %5.2f'%(101.0),101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        Error = G2mth.FitTexture(General,Gangls,refData,keyList,pgbar)
        pgbar.Destroy()
        if Error:
            wx.MessageBox(Error,caption='Fit Texture Error',style=wx.ICON_EXCLAMATION)
#        x = []
#        y = []
        XY = []
        for hist in keyList:
            x = refData[hist].T[5].T
            y = refData[hist].T[6].T
            xy = [x,y]
            XY.append(np.array(xy))
#        XY = np.array(XY)
        G2plt.PlotXY(G2frame,XY,XY2=[],labelX='POobs',labelY='POcalc',newPlot=False,Title='Texture fit error')
        UpdateTexture()
        G2plt.PlotTexture(G2frame,data,Start=False)            
            
    def OnTextureClear(event):
        print ('clear texture? - does nothing')

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
            if menu.FindItem(page) < 0: # is tab already in menu?
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
        ChangePage(page)
        
    def ChangePage(page):
        for p in G2frame.phaseDisplay.gridList: # clear out all grids, forcing edits in progress to complete
            p.ClearGrid()
        text = G2frame.phaseDisplay.GetPageText(page)
        G2frame.lastSelectedPhaseTab = text
        G2frame.dataWindow.helpKey = text # use name of Phase tab for help lookup
        if text == 'General':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataGeneral)
            UpdateGeneral()
        elif text == 'Data':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataMenu)
            G2ddG.UpdateDData(G2frame,DData,data)
            wx.CallAfter(G2plt.PlotSizeStrainPO,G2frame,data,hist='',Start=True)            
        elif text == 'Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.AtomsMenu)
            FillAtomsGrid(Atoms)
        elif text == 'Layers':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.LayerData)
            UpdateLayerData()
        elif text == 'Wave Data' and data['General']['Modulated']:
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.WavesData)
            UpdateWavesData()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Draw Options':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataDrawOptions)
            UpdateDrawOptions()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Draw Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DrawAtomsMenu)
            UpdateDrawAtoms()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'RB Models':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodiesMenu)
            FillRigidBodyGrid()
        elif text == 'Map peaks':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.MapPeaksMenu)
            FillMapPeaksGrid()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'MC/SA':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.MCSAMenu)
            UpdateMCSA()                        
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Texture':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.TextureMenu)
            UpdateTexture()                        
            wx.CallAfter(G2plt.PlotTexture,G2frame,data,Start=True)            
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
        G2frame.Bind(wx.EVT_MENU, OnUseBilbao, id=G2G.wxID_USEBILBAOMAG)
        G2frame.Bind(wx.EVT_MENU, OnValidProtein, id=G2G.wxID_VALIDPROTEIN)
        # Data
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.DataMenu)
        G2frame.Bind(wx.EVT_MENU, OnDataUse, id=G2G.wxID_DATAUSE)
        G2frame.Bind(wx.EVT_MENU, OnDataCopy, id=G2G.wxID_DATACOPY)
        G2frame.Bind(wx.EVT_MENU, OnDataCopyFlags, id=G2G.wxID_DATACOPYFLAGS)
        G2frame.Bind(wx.EVT_MENU, OnSelDataCopy, id=G2G.wxID_DATASELCOPY)
        G2frame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2G.wxID_PWDRADD)
        G2frame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2G.wxID_HKLFADD)
        G2frame.Bind(wx.EVT_MENU, OnDataDelete, id=G2G.wxID_DATADELETE)
        G2frame.Bind(wx.EVT_MENU, OnDataApplyStrain, id=G2G.wxID_DATADIJ)
        # Atoms
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.AtomsMenu)
        G2frame.Bind(wx.EVT_MENU, OnSetAll, id=G2G.wxID_ATOMSSETALL)
        G2frame.Bind(wx.EVT_MENU, AtomRefine, id=G2G.wxID_ATOMSSETSEL)
        G2frame.Bind(wx.EVT_MENU, AtomModify, id=G2G.wxID_ATOMSMODIFY)
        G2frame.Bind(wx.EVT_MENU, OnAtomInsert, id=G2G.wxID_ATOMSEDITINSERT)
        G2frame.Bind(wx.EVT_MENU, OnHydAtomAdd, id=G2G.wxID_ADDHATOM)
        G2frame.Bind(wx.EVT_MENU, AtomDelete, id=G2G.wxID_ATOMSEDITDELETE)
        G2frame.Bind(wx.EVT_MENU, AtomTransform, id=G2G.wxID_ATOMSTRANSFORM)
#        G2frame.Bind(wx.EVT_MENU, AtomRotate, id=G2G.wxID_ATOMSROTATE)
        
        G2frame.Bind(wx.EVT_MENU, OnAtomAdd, id=G2G.wxID_ATOMSEDITADD)
        G2frame.Bind(wx.EVT_MENU, OnAtomViewAdd, id=G2G.wxID_ATOMSVIEWADD)
        G2frame.Bind(wx.EVT_MENU, OnAtomViewInsert, id=G2G.wxID_ATOMVIEWINSERT)
        G2frame.Bind(wx.EVT_MENU, OnHydAtomUpdate, id=G2G.wxID_UPDATEHATOM)
        G2frame.Bind(wx.EVT_MENU, OnAtomMove, id=G2G.wxID_ATOMMOVE)
        G2frame.Bind(wx.EVT_MENU, MakeMolecule, id=G2G.wxID_MAKEMOLECULE)
        G2frame.Bind(wx.EVT_MENU, OnReloadDrawAtoms, id=G2G.wxID_RELOADDRAWATOMS)
        G2frame.Bind(wx.EVT_MENU, OnDistAngle, id=G2G.wxID_ATOMSDISAGL)
        G2frame.Bind(wx.EVT_MENU, OnDistAnglePrt, id=G2G.wxID_ATOMSPDISAGL)
        G2frame.Bind(wx.EVT_MENU, OnDensity, id=G2G.wxID_ATOMSDENSITY)
        G2frame.Bind(wx.EVT_MENU, OnIsoDistortCalc, id=G2G.wxID_ISODISP)
        if 'HydIds' in data['General']:
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,True)
        else:
            G2frame.dataWindow.AtomEdit.Enable(G2G.wxID_UPDATEHATOM,False)
        for id in G2frame.dataWindow.ReImportMenuId:     #loop over submenu items
            G2frame.Bind(wx.EVT_MENU, OnReImport, id=id)
        # Wave Data
        if data['General']['Modulated']:
            FillSelectPageMenu(TabSelectionIdDict, G2frame.dataWindow.WavesData)
            G2frame.Bind(wx.EVT_MENU, OnWaveVary, id=G2G.wxID_WAVEVARY)
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
        G2frame.Bind(wx.EVT_MENU, TransformSymEquiv, id=G2G.wxID_DRAWTRANSFORM)
        G2frame.Bind(wx.EVT_MENU, FillCoordSphere, id=G2G.wxID_DRAWFILLCOORD)            
        G2frame.Bind(wx.EVT_MENU, FillUnitCell, id=G2G.wxID_DRAWFILLCELL)
        G2frame.Bind(wx.EVT_MENU, DrawAtomsDelete, id=G2G.wxID_DRAWDELETE)
        G2frame.Bind(wx.EVT_MENU, OnDrawDistVP, id=G2G.wxID_DRAWDISTVP)
        G2frame.Bind(wx.EVT_MENU, OnDrawDAT, id=G2G.wxID_DRAWDISAGLTOR)
        G2frame.Bind(wx.EVT_MENU, OnDrawPlane, id=G2G.wxID_DRAWPLANE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRBOND)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRANGLE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRPLANE)
        G2frame.Bind(wx.EVT_MENU, OnRestraint, id=G2G.wxID_DRAWRESTRCHIRAL)
        G2frame.Bind(wx.EVT_MENU, OnDefineRB, id=G2G.wxID_DRAWDEFINERB)
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
        G2frame.Bind(wx.EVT_MENU, OnPeaksUnique, id=G2G.wxID_PEAKSUNIQUE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksDelete, id=G2G.wxID_PEAKSDELETE)
        G2frame.Bind(wx.EVT_MENU, OnPeaksClear, id=G2G.wxID_PEAKSCLEAR)
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
        
    # UpdatePhaseData execution starts here
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
    if 'modulated' in data['General']['Type']:
        data['General']['Modulated'] = True
        data['General']['Type'] = 'nuclear'     
#end patch    

    global rbAtmDict   
    rbAtmDict = {}
    PhaseName = G2frame.GPXtree.GetItemText(Item)
    G2gd.SetDataMenuBar(G2frame)
    # Bob: why do this differently in debug mode? Is this code to test if tabs can be moved around? #TODO - yup, flaky tho.
#    if GSASIIpath.GetConfigValue('debug'):
#        G2frame.phaseDisplay = G2G.GSNoteBook(parent=G2frame.dataWindow,size=G2frame.dataWindow.GetClientSize(),
#            style=wx.aui.AUI_NB_TOP | wx.aui.AUI_NB_TAB_SPLIT | wx.aui.AUI_NB_TAB_MOVE)
#    else:
#        G2frame.phaseDisplay = G2G.GSNoteBook(parent=G2frame.dataWindow,size=G2frame.dataWindow.GetClientSize())
    G2frame.phaseDisplay = G2G.GSNoteBook(parent=G2frame.dataWindow)
    G2frame.dataWindow.GetSizer().Add(G2frame.phaseDisplay,1,wx.ALL|wx.EXPAND,1)
    G2frame.phaseDisplay.gridList = [] # list of all grids in notebook
    Pages = []    
    General = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(General,'General')
    Pages.append('General')
    DData = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(DData,'Data')
    Pages.append('Data')
    Atoms = G2G.GSGrid(G2frame.phaseDisplay)
#    Atoms.SetScrollRate(0,0)
    G2frame.phaseDisplay.gridList.append(Atoms)
    G2frame.phaseDisplay.AddPage(Atoms,'Atoms')
    Pages.append('Atoms')
    if data['General']['Modulated']:
        G2frame.waveData = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.waveData,'Wave Data')
        Pages.append('Wave Data') 
    if data['General']['Type'] == 'faulted':
        G2frame.layerData = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.layerData,'Layers')
        Pages.append('Layers')               
    drawOptions = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(drawOptions,'Draw Options')
    Pages.append('Draw Options')
    drawAtoms = G2G.GSGrid(G2frame.phaseDisplay)
#    drawAtoms.SetScrollRate(0,0)
    G2frame.phaseDisplay.gridList.append(drawAtoms)
    G2frame.phaseDisplay.AddPage(drawAtoms,'Draw Atoms')
    Pages.append('Draw Atoms')
    if data['General']['Type'] not in ['faulted',] and not data['General']['Modulated']:
        RigidBodies = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(RigidBodies,'RB Models')
        Pages.append('RB Models')
    MapPeaks = G2G.GSGrid(G2frame.phaseDisplay)
#    MapPeaks.SetScrollRate(0,0)
    G2frame.phaseDisplay.gridList.append(MapPeaks)    
    G2frame.phaseDisplay.AddPage(MapPeaks,'Map peaks')
    Pages.append('Map peaks')
    if data['General']['Type'] not in ['faulted',] and not data['General']['Modulated']:
        G2frame.MCSA = wx.ScrolledWindow(G2frame.phaseDisplay)
        G2frame.phaseDisplay.AddPage(G2frame.MCSA,'MC/SA')
        Pages.append('MC/SA')
    Texture = wx.ScrolledWindow(G2frame.phaseDisplay)
    G2frame.phaseDisplay.AddPage(Texture,'Texture')
    Pages.append('Texture')
    G2frame.PawleyRefl = G2G.GSGrid(G2frame.phaseDisplay)
#    G2frame.PawleyRefl.SetScrollRate(0,0)
    G2frame.phaseDisplay.gridList.append(G2frame.PawleyRefl)
    G2frame.phaseDisplay.AddPage(G2frame.PawleyRefl,'Pawley reflections')
    Pages.append('Pawley reflections')
    G2frame.dataWindow.AtomCompute.Enable(G2G.wxID_ISODISP,'ISODISTORT' in data)
    G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_VALIDPROTEIN,'macro' in data['General']['Type'])
    G2frame.dataWindow.GeneralCalc.Enable(G2G.wxID_USEBILBAOMAG,'magPhases' in data)
    G2frame.phaseDisplay.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    FillMenus()
    if G2frame.lastSelectedPhaseTab in Pages:
        ind = Pages.index(G2frame.lastSelectedPhaseTab)
        if ind:
            UpdateGeneral(SkipDraw=ind)
            G2frame.phaseDisplay.SetSelection(ind)
            return
    ChangePage(0)
