# -*- coding: utf-8 -*-
#GSASII - phase data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import wx
import wx.grid as wg
import matplotlib as mpl
import math
import copy
import time
import sys
import random as ran
import cPickle
import GSASIIpath
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIIO as G2IO
import GSASIIstruct as G2str
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma

VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)
BLACK = wx.Colour(0,0,0)

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi

class SymOpDialog(wx.Dialog):
    def __init__(self,parent,SGData,New=True):
        wx.Dialog.__init__(self,parent,-1,'Select symmetry operator',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        self.SGData = SGData
        self.New = New
        self.OpSelected = [0,0,0,[0,0,0],False]
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGInv']:
            choice = ['No','Yes']
            self.inv = wx.RadioBox(panel,-1,'Choose inversion?',choices=choice)
            self.inv.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.inv,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGLatt'] != 'P':
            LattOp = G2spc.Latt2text(SGData['SGLatt']).split(';')
            self.latt = wx.RadioBox(panel,-1,'Choose cell centering?',choices=LattOp)
            self.latt.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.latt,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGLaue'] in ['-1','2/m','mmm','4/m','4/mmm']:
            Ncol = 2
        else:
            Ncol = 3
        OpList = []
        for M,T in SGData['SGOps']:
            OpList.append(G2spc.MT2text(M,T))
        self.oprs = wx.RadioBox(panel,-1,'Choose space group operator?',choices=OpList,
            majorDimension=Ncol)
        self.oprs.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.oprs,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(panel,-1,"   Choose unit cell?"),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellSizer.Add((5,0),0)
        cellName = ['X','Y','Z']
        self.cell = []
        for i in range(3):
            self.cell.append(wx.SpinCtrl(panel,-1,cellName[i],size=wx.Size(50,20)))
            self.cell[-1].SetRange(-3,3)
            self.cell[-1].SetValue(0)
            self.cell[-1].Bind(wx.EVT_SPINCTRL, self.OnOpSelect)
            cellSizer.Add(self.cell[-1],0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(cellSizer,0,)
        if self.New:
            choice = ['No','Yes']
            self.new = wx.RadioBox(panel,-1,'Generate new positions?',choices=choice)
            self.new.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.new,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)

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
        if self.SGData['SGInv']:
            self.OpSelected[0] = self.inv.GetSelection()
        if self.SGData['SGLatt'] != 'P':
            self.OpSelected[1] = self.latt.GetSelection()
        self.OpSelected[2] = self.oprs.GetSelection()
        for i in range(3):
            self.OpSelected[3][i] = float(self.cell[i].GetValue())
        if self.New:
            self.OpSelected[4] = self.new.GetSelection()

    def GetSelection(self):
        return self.OpSelected

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)
        self.Destroy()

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)
        self.Destroy()

class DisAglDialog(wx.Dialog):
    
    def __default__(self,data,default):
        if data:
            self.data = data
        else:
            self.data = {}
            self.data['Name'] = default['Name']
            self.data['Factors'] = [0.85,0.85]
            self.data['AtomTypes'] = default['AtomTypes']
            self.data['BondRadii'] = default['BondRadii']
            self.data['AngleRadii'] = default['AngleRadii']
        
    def __init__(self,parent,data,default):
        wx.Dialog.__init__(self,parent,-1,'Distance Angle Controls', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.default = default
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.__default__(data,self.default)
        self.Draw(self.data)
                
    def Draw(self,data):
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,'Controls for phase '+data['Name']),
            0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
        mainSizer.Add((10,10),1)
        
        radiiSizer = wx.FlexGridSizer(2,3,5,5)
        radiiSizer.Add(wx.StaticText(self.panel,-1,' Type'),0,wx.ALIGN_CENTER_VERTICAL)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Bond radii'),0,wx.ALIGN_CENTER_VERTICAL)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Angle radii'),0,wx.ALIGN_CENTER_VERTICAL)
        self.objList = {}
        for id,item in enumerate(self.data['AtomTypes']):
            radiiSizer.Add(wx.StaticText(self.panel,-1,' '+item),0,wx.ALIGN_CENTER_VERTICAL)
            bRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['BondRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[bRadii.GetId()] = ['BondRadii',id]
            bRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(bRadii,0,wx.ALIGN_CENTER_VERTICAL)
            aRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['AngleRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[aRadii.GetId()] = ['AngleRadii',id]
            aRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            aRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(aRadii,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(radiiSizer,0,wx.EXPAND)
        factorSizer = wx.FlexGridSizer(2,2,5,5)
        Names = ['Bond','Angle']
        for i,name in enumerate(Names):
            factorSizer.Add(wx.StaticText(self.panel,-1,name+' search factor'),0,wx.ALIGN_CENTER_VERTICAL)
            bondFact = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['Factors'][i]),style=wx.TE_PROCESS_ENTER)
            self.objList[bondFact.GetId()] = ['Factors',i]
            bondFact.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bondFact.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            factorSizer.Add(bondFact)
        mainSizer.Add(factorSizer,0,wx.EXPAND)
        
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        ResetBtn = wx.Button(self.panel,-1,'Reset')
        ResetBtn.Bind(wx.EVT_BUTTON, self.OnReset)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(ResetBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
    
    def OnRadiiVal(self,event):
        Obj = event.GetEventObject()
        item = self.objList[Obj.GetId()]
        try:
            self.data[item[0]][item[1]] = float(Obj.GetValue())
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(self.data[item[0]][item[1]]))          #reset in case of error
        
    def GetData(self):
        return self.data
        
    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        self.Destroy()
        
    def OnReset(self,event):
        data = {}
        self.__default__(data,self.default)
        self.Draw(self.data)
        
class SingleFloatDialog(wx.Dialog):
    
    def __init__(self,parent,title,prompt,value,limits=[0.,1.]):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.limits = limits
        self.value = value
        self.prompt = prompt
        self.Draw()
        
    def Draw(self):
        
        def OnValItem(event):
            try:
                val = float(valItem.GetValue())
                if val < self.limits[0] or val > self.limits[1]:
                    raise ValueError
            except ValueError:
                val = self.value
            self.value = val
            valItem.SetValue('%.5g'%(self.value))
            
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,self.prompt),0,wx.ALIGN_CENTER)
        valItem = wx.TextCtrl(self.panel,-1,value='%.5g'%(self.value),style=wx.TE_PROCESS_ENTER)
        mainSizer.Add(valItem,0,wx.ALIGN_CENTER)
        valItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
        valItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def GetValue(self):
        return self.value
        
    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        self.Destroy()
        
    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)              
        self.Destroy()
        
def UpdatePhaseData(G2frame,Item,data,oldPage):

    Atoms = []
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    PhaseName = G2frame.PatternTree.GetItemText(Item)
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
    G2frame.dataFrame.SetLabel('Phase Data for '+PhaseName)
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataDisplay = G2gd.GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())

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
            generalData['Map'] = {'MapType':'','RefList':'','Resolution':1.0,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}
        if 'Flip' not in generalData:
            generalData['Flip'] = {'RefList':'','Resolution':1.0,'Norm element':'None',
                'k-factor':1.1}
        if 'doPawley' not in generalData:
            generalData['doPawley'] = False
        if 'Pawley dmin' not in generalData:
            generalData['Pawley dmin'] = 1.0
            
#        if 'SH Texture' not in generalData:
#            generalData['SH Texture'] = data['SH Texture']
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['vdWRadii'] = []
        generalData['AtomMass'] = []
        generalData['Color'] = []
        generalData['Mydir'] = G2frame.dirname
        cx,ct,cs,cia = [3,1,7,9]
        generalData['AtomPtrs'] = [cx,ct,cs,cia]
        if generalData['Type'] =='macromolecular':
            cx,ct,cs,cia = [6,4,10,12]
            generalData['AtomPtrs'] = [cx,ct,cs,cia]
        for atom in atomData:
            atom[ct] = atom[ct].lower().capitalize()              #force to standard form
            if generalData['AtomTypes'].count(atom[ct]):
                generalData['NoAtoms'][atom[ct]] += atom[cs-1]*float(atom[cs+1])
            elif atom[ct] != 'UNK':
                Info = G2elem.GetAtomInfo(atom[ct])
                generalData['AtomTypes'].append(atom[ct])
                generalData['Z'] = Info['Z']
                generalData['Isotopes'][atom[ct]] = Info['Isotopes']
                generalData['BondRadii'].append(Info['Drad'])
                generalData['AngleRadii'].append(Info['Arad'])
                generalData['vdWRadii'].append(Info['Vdrad'])
                if atom[ct] in generalData['Isotope']:
                    generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]][0])
                else:
                    generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                    generalData['AtomMass'].append(Info['Mass'])
                generalData['NoAtoms'][atom[ct]] = atom[cs-1]*float(atom[cs+1])
                generalData['Color'].append(Info['Color'])
        F000X = 0.
        F000N = 0.
        for i,elem in enumerate(generalData['AtomTypes']):
            F000X += generalData['NoAtoms'][elem]*generalData['Z']
            isotope = generalData['Isotope'][elem]
            F000N += generalData['NoAtoms'][elem]*generalData['Isotopes'][elem][isotope][1]
        generalData['F000X'] = F000X
        generalData['F000N'] = F000N
       

################################################################################
##### General phase routines
################################################################################

    def UpdateGeneral():
        
        ''' default dictionary structure for phase data: (taken from GSASII.py)
        'General':{
            'Name':PhaseName
            'Type':'nuclear'
            'SGData':SGData
            'Cell':[False,10.,10.,10.,90.,90.,90,1000.]
            'AtomPtrs':[]
            'Histogram list':['',]
            'Pawley dmin':1.0}
        'Atoms':[]
        'Drawing':{}
        '''
        
        phaseTypes = ['nuclear','modulated','magnetic','macromolecular']
        SetupGeneral()
        generalData = data['General']
        Map = generalData['Map']
        Flip = generalData['Flip']  
        
        def NameSizer():
                   
            def OnPhaseName(event):
                oldName = generalData['Name']
                generalData['Name'] = NameTxt.GetValue()
                G2frame.G2plotNB.Rename(oldName,generalData['Name'])
                G2frame.dataFrame.SetLabel('Phase Data for '+generalData['Name'])
                G2frame.PatternTree.SetItemText(Item,generalData['Name'])
                #Hmm, need to change phase name key in Reflection Lists for each histogram
                            
            def OnPhaseType(event):
                if not generalData['AtomTypes']:             #can change only if no atoms!
                    generalData['Type'] = TypeTxt.GetValue()
                    dataDisplay.DestroyChildren()           #needed to clear away bad cellSizer, etc.
                    UpdateGeneral()         #must use this way!
                else:
                    TypeTxt.SetValue(generalData['Type'])                
                
            def OnSpaceGroup(event):
                SpcGp = SGTxt.GetValue()
                SGErr,SGData = G2spc.SpcGroup(SpcGp)
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    SGTxt.SetValue(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Error'
                    Style = wx.ICON_EXCLAMATION
                else:
                    text = G2spc.SGPrint(SGData)
                    generalData['SGData'] = SGData
                    msg = 'Space Group Information'
                    Style = wx.ICON_INFORMATION
                Text = ''
                for line in text:
                    Text += line+'\n'
                wx.MessageBox(Text,caption=msg,style=Style)
                dataDisplay.DestroyChildren()           #needed to clear away bad cellSizer, etc.
                UpdateGeneral()
                
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(dataDisplay,-1,' Phase name: '),0,wx.ALIGN_CENTER_VERTICAL)
            NameTxt = wx.TextCtrl(dataDisplay,-1,value=generalData['Name'],style=wx.TE_PROCESS_ENTER)
            NameTxt.Bind(wx.EVT_TEXT_ENTER,OnPhaseName)
            NameTxt.Bind(wx.EVT_KILL_FOCUS,OnPhaseName)
            nameSizer.Add(NameTxt,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add(wx.StaticText(dataDisplay,-1,'  Phase type: '),0,wx.ALIGN_CENTER_VERTICAL)
            if len(data['Atoms']):
                choices = phaseTypes[:-1]
            else:
                choices = phaseTypes            
            TypeTxt = wx.ComboBox(dataDisplay,-1,value=generalData['Type'],choices=choices,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            TypeTxt.Bind(wx.EVT_COMBOBOX, OnPhaseType)
            nameSizer.Add(TypeTxt,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add(wx.StaticText(dataDisplay,-1,'  Space group: '),0,wx.ALIGN_CENTER_VERTICAL)
            SGTxt = wx.TextCtrl(dataDisplay,-1,value=generalData['SGData']['SpGrp'],style=wx.TE_PROCESS_ENTER)
            SGTxt.Bind(wx.EVT_TEXT_ENTER,OnSpaceGroup)
            nameSizer.Add(SGTxt,0,wx.ALIGN_CENTER_VERTICAL)
            return nameSizer
            
        def CellSizer():
            
            cellGUIlist = [[['m3','m3m'],4,zip([" Unit cell: a = "," Vol = "],["%.5f","%.3f"],[True,False],[0,0])],
            [['3R','3mR'],6,zip([" a = "," alpha = "," Vol = "],["%.5f","%.3f","%.3f"],[True,True,False],[0,2,0])],
            [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],6,zip([" a = "," c = "," Vol = "],["%.5f","%.5f","%.3f"],[True,True,False],[0,2,0])],
            [['mmm'],8,zip([" a = "," b = "," c = "," Vol = "],["%.5f","%.5f","%.5f","%.3f"],
                [True,True,True,False],[0,1,2,0])],
            [['2/m'+'a'],10,zip([" a = "," b = "," c = "," alpha = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
            [['2/m'+'b'],10,zip([" a = "," b = "," c = "," beta = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
            [['2/m'+'c'],10,zip([" a = "," b = "," c = "," gamma = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
            [['-1'],8,zip([" a = "," b = "," c = "," Vol = "," alpha = "," beta = "," gamma = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f","%.3f","%.3f"],
                [True,True,True,False,True,True,True],[0,1,2,0,3,4,5])]]
                
            def OnCellRef(event):
                generalData['Cell'][0] = cellRef.GetValue()
                
            def OnCellChange(event):
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
                        value = controls[6+ObjId]
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
                    cell[5] = cell[6] = 90.
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
                generalData['Cell'] = cell
                dataDisplay.DestroyChildren()           #needed to clear away bad cellSizer, etc.
                UpdateGeneral()
            
            cell = generalData['Cell']
            laue = generalData['SGData']['SGLaue']
            if laue == '2/m':
                laue += generalData['SGData']['SGUniq']
            for cellGUI in cellGUIlist:
                if laue in cellGUI[0]:
                    useGUI = cellGUI
            cellSizer = wx.FlexGridSizer(2,useGUI[1]+1,5,5)
            cellRef = wx.CheckBox(dataDisplay,-1,label='Refine unit cell:')
            cellSizer.Add(cellRef,0,wx.ALIGN_CENTER_VERTICAL)
            cellRef.Bind(wx.EVT_CHECKBOX, OnCellRef)
            cellRef.SetValue(cell[0])
            cellList = []
            for txt,fmt,ifEdit,Id in useGUI[2]:
                cellSizer.Add(wx.StaticText(dataDisplay,label=txt),0,wx.ALIGN_CENTER_VERTICAL)
                if ifEdit:          #a,b,c,etc.
                    cellVal = wx.TextCtrl(dataDisplay,value=(fmt%(cell[Id+1])),
                        style=wx.TE_PROCESS_ENTER)
                    cellVal.Bind(wx.EVT_TEXT_ENTER,OnCellChange)        
                    cellVal.Bind(wx.EVT_KILL_FOCUS,OnCellChange)
                    cellSizer.Add(cellVal,0,wx.ALIGN_CENTER_VERTICAL)
                    cellList.append(cellVal.GetId())
                else:               #volume
                    volVal = wx.TextCtrl(dataDisplay,value=(fmt%(cell[7])),style=wx.TE_READONLY)
                    volVal.SetBackgroundColour(VERY_LIGHT_GREY)
                    cellSizer.Add(volVal,0,wx.ALIGN_CENTER_VERTICAL)
            return cellSizer
            
        def ElemSizer():
            
            def OnIsotope(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                isotope = Obj.GetValue()
                generalData['Isotope'][item] = isotope
                indx = generalData['AtomTypes'].index(item)
                data['General']['AtomMass'][indx] = generalData['Isotopes'][item][isotope][0]
                dataDisplay.DestroyChildren()           #needed to clear away bad cellSizer, etc.
                UpdateGeneral()
                
            elemSizer = wx.FlexGridSizer(8,len(generalData['AtomTypes'])+1,1,1)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Elements'),0,wx.ALIGN_CENTER_VERTICAL)
            for elem in generalData['AtomTypes']:
                typTxt = wx.TextCtrl(dataDisplay,value=elem,style=wx.TE_READONLY)
                typTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(typTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Isotope'),0,wx.ALIGN_CENTER_VERTICAL)
            for elem in generalData['AtomTypes']:
                choices = generalData['Isotopes'][elem].keys()
                isoSel = wx.ComboBox(dataDisplay,-1,value=generalData['Isotope'][elem],choices=choices,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                isoSel.Bind(wx.EVT_COMBOBOX,OnIsotope)
                Indx[isoSel.GetId()] = elem
                elemSizer.Add(isoSel,1,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' No. per cell'),0,wx.ALIGN_CENTER_VERTICAL)
            for elem in generalData['AtomTypes']:
                numbTxt = wx.TextCtrl(dataDisplay,value='%.1f'%(generalData['NoAtoms'][elem]),
                    style=wx.TE_READONLY)
                numbTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(numbTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Atom weight'),0,wx.ALIGN_CENTER_VERTICAL)
            for wt in generalData['AtomMass']:
                wtTxt = wx.TextCtrl(dataDisplay,value='%.3f'%(wt),style=wx.TE_READONLY)
                wtTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(wtTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Bond radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['BondRadii']:
                bondRadii = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                bondRadii.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(bondRadii,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Angle radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['AngleRadii']:
                elemTxt = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' van der Waals radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['vdWRadii']:
                elemTxt = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label=' Default color'),0,wx.ALIGN_CENTER_VERTICAL)
            for R,G,B in generalData['Color']:
                colorTxt = wx.TextCtrl(dataDisplay,value='',style=wx.TE_READONLY)
                colorTxt.SetBackgroundColour(wx.Colour(R,G,B))
                elemSizer.Add(colorTxt,0,wx.ALIGN_CENTER_VERTICAL)
            return elemSizer
            
        def DenSizer():
            
            mass = 0.
            for i,elem in enumerate(generalData['AtomTypes']):
                mass += generalData['NoAtoms'][elem]*generalData['AtomMass'][i]
            denSizer = wx.BoxSizer(wx.HORIZONTAL)
            denSizer.Add(wx.StaticText(dataDisplay,-1,' Density: '),0,wx.ALIGN_CENTER_VERTICAL)
            Volume = generalData['Cell'][7]
            density = mass/(0.6022137*Volume)
            denTxt = wx.TextCtrl(dataDisplay,-1,'%.3f'%(density),style=wx.TE_READONLY)
            denTxt.SetBackgroundColour(VERY_LIGHT_GREY)
            denSizer.Add(denTxt,0,wx.ALIGN_CENTER_VERTICAL)        
            if generalData['Type'] == 'macromolecular' and mass > 0.0:
                denSizer.Add(wx.StaticText(dataDisplay,-1,' Matthews coeff.: '),
                    0,wx.ALIGN_CENTER_VERTICAL)
                mattTxt = wx.TextCtrl(dataDisplay,-1,'%.3f'%(Volume/mass),style=wx.TE_READONLY)
                mattTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                denSizer.Add(mattTxt,0,wx.ALIGN_CENTER_VERTICAL)
            return denSizer
            
        def PawleySizer():
            
            def OnPawleyRef(event):
                generalData['doPawley'] = pawlRef.GetValue()
            
            def OnPawleyVal(event):
                try:
                    dmin = float(pawlVal.GetValue())
                    if 0.25 <= dmin <= 20.:
                        generalData['Pawley dmin'] = dmin
                except ValueError:
                    pass
                pawlVal.SetValue("%.3f"%(generalData['Pawley dmin']))          #reset in case of error                
            
            pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
            pawleySizer.Add(wx.StaticText(dataDisplay,label=' Pawley controls: '),0,wx.ALIGN_CENTER_VERTICAL)
            pawlRef = wx.CheckBox(dataDisplay,-1,label=' Do Pawley refinement?')
            pawlRef.SetValue(generalData['doPawley'])
            pawlRef.Bind(wx.EVT_CHECKBOX,OnPawleyRef)
            pawleySizer.Add(pawlRef,0,wx.ALIGN_CENTER_VERTICAL)
            pawleySizer.Add(wx.StaticText(dataDisplay,label=' Pawley dmin: '),0,wx.ALIGN_CENTER_VERTICAL)
            pawlVal = wx.TextCtrl(dataDisplay,value='%.3f'%(generalData['Pawley dmin']),style=wx.TE_PROCESS_ENTER)
            pawlVal.Bind(wx.EVT_TEXT_ENTER,OnPawleyVal)        
            pawlVal.Bind(wx.EVT_KILL_FOCUS,OnPawleyVal)
            pawleySizer.Add(pawlVal,0,wx.ALIGN_CENTER_VERTICAL)
            return pawleySizer
            
        def MapSizer():
            
            def OnMapType(event):
                Map['MapType'] = mapType.GetValue()
                
            def OnRefList(event):
                Map['RefList'] = refList.GetValue()
                
            def OnResVal(event):
                try:
                    res = float(mapRes.GetValue())
                    if 0.25 <= res <= 20.:
                        Map['Resolution'] = res
                except ValueError:
                    pass
                mapRes.SetValue("%.2f"%(Map['Resolution']))          #reset in case of error
            
            def OnCutOff(event):
                try:
                    res = float(cutOff.GetValue())
                    if 1.0 <= res <= 100.:
                        Map['cutOff'] = res
                except ValueError:
                    pass
                cutOff.SetValue("%.1f"%(Map['cutOff']))          #reset in case of error
            
            #patch
            if 'cutOff' not in Map:
                Map['cutOff'] = 100.0
            mapTypes = ['Fobs','Fcalc','delt-F','2*Fo-Fc','Patterson']
            refList = data['Histograms'].keys()
            if not generalData['AtomTypes']:
                 mapTypes = ['Patterson',]
                 Map['MapType'] = 'Patterson'
            mapSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(dataDisplay,label=' Fourier map controls: Map type: '),0,wx.ALIGN_CENTER_VERTICAL)
            mapType = wx.ComboBox(dataDisplay,-1,value=Map['MapType'],choices=mapTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            mapType.Bind(wx.EVT_COMBOBOX,OnMapType)
            lineSizer.Add(mapType,0,wx.ALIGN_CENTER_VERTICAL)
            lineSizer.Add(wx.StaticText(dataDisplay,label=' Reflection set from: '),0,wx.ALIGN_CENTER_VERTICAL)
            refList = wx.ComboBox(dataDisplay,-1,value=Map['RefList'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,wx.ALIGN_CENTER_VERTICAL)
            mapSizer.Add(lineSizer,0,wx.ALIGN_CENTER_VERTICAL)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(dataDisplay,label=' Resolution: '),0,wx.ALIGN_CENTER_VERTICAL)
            mapRes =  wx.TextCtrl(dataDisplay,value='%.2f'%(Map['Resolution']),style=wx.TE_PROCESS_ENTER)
            mapRes.Bind(wx.EVT_TEXT_ENTER,OnResVal)        
            mapRes.Bind(wx.EVT_KILL_FOCUS,OnResVal)
            line2Sizer.Add(mapRes,0,wx.ALIGN_CENTER_VERTICAL)
            line2Sizer.Add(wx.StaticText(dataDisplay,label=' Peak cutoff %: '),0,wx.ALIGN_CENTER_VERTICAL)
            cutOff =  wx.TextCtrl(dataDisplay,value='%.1f'%(Map['cutOff']),style=wx.TE_PROCESS_ENTER)
            cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)        
            cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
            line2Sizer.Add(cutOff,0,wx.ALIGN_CENTER_VERTICAL)
            mapSizer.Add(line2Sizer,0,wx.ALIGN_CENTER_VERTICAL)
            return mapSizer
                
        def FlipSizer():
            
            def OnRefList(event):
                Flip['RefList'] = refList.GetValue()
                
            def OnNormElem(event):
                PE = G2elemGUI.PickElement(G2frame,ifNone=True)
                if PE.ShowModal() == wx.ID_OK:
                    Flip['Norm element'] = PE.Elem.strip()
                    normElem.SetLabel(Flip['Norm element'])
                PE.Destroy()                
                
            def OnResVal(event):
                try:
                    res = float(flipRes.GetValue())
                    if 0.25 <= res <= 20.:
                        Flip['Resolution'] = res
                except ValueError:
                    pass
                flipRes.SetValue("%.2f"%(Flip['Resolution']))          #reset in case of error
            
            def OnkFactor(event):
                try:
                    res = float(kFactor.GetValue())
                    if 0.1 <= res <= 1.2:
                        Flip['k-factor'] = res
                except ValueError:
                    pass
                kFactor.SetValue("%.3f"%(Flip['k-factor']))          #reset in case of error
            
            refList = data['Histograms'].keys()
            flipSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(dataDisplay,label=' Charge flip controls: Reflection set from: '),0,wx.ALIGN_CENTER_VERTICAL)
            refList = wx.ComboBox(dataDisplay,-1,value=Flip['RefList'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,wx.ALIGN_CENTER_VERTICAL)
            flipSizer.Add(lineSizer,0,wx.ALIGN_CENTER_VERTICAL)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(dataDisplay,label=' Normalizing element: '),0,wx.ALIGN_CENTER_VERTICAL)
            normElem = wx.Button(dataDisplay,label=Flip['Norm element'],style=wx.TE_READONLY)
            normElem.Bind(wx.EVT_BUTTON,OnNormElem)
            line2Sizer.Add(normElem,0,wx.ALIGN_CENTER_VERTICAL)
            line2Sizer.Add(wx.StaticText(dataDisplay,label=' Resolution: '),0,wx.ALIGN_CENTER_VERTICAL)
            flipRes =  wx.TextCtrl(dataDisplay,value='%.2f'%(Flip['Resolution']),style=wx.TE_PROCESS_ENTER)
            flipRes.Bind(wx.EVT_TEXT_ENTER,OnResVal)        
            flipRes.Bind(wx.EVT_KILL_FOCUS,OnResVal)
            line2Sizer.Add(flipRes,0,wx.ALIGN_CENTER_VERTICAL)
            line2Sizer.Add(wx.StaticText(dataDisplay,label=' k-Factor (0.1-1.2): '),0,wx.ALIGN_CENTER_VERTICAL)
            kFactor =  wx.TextCtrl(dataDisplay,value='%.3f'%(Flip['k-factor']),style=wx.TE_PROCESS_ENTER)
            kFactor.Bind(wx.EVT_TEXT_ENTER,OnkFactor)        
            kFactor.Bind(wx.EVT_KILL_FOCUS,OnkFactor)
            line2Sizer.Add(kFactor,0,wx.ALIGN_CENTER_VERTICAL)
            flipSizer.Add(line2Sizer,0,wx.ALIGN_CENTER_VERTICAL)
            return flipSizer
                
        General.DestroyChildren()
        dataDisplay = wx.Panel(General)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(NameSizer(),0)
        mainSizer.Add((5,5),0)        
        mainSizer.Add(CellSizer(),0)
        mainSizer.Add((5,5),0)
        
        Indx = {}
        if len(generalData['AtomTypes']):
            mainSizer.Add(DenSizer())
            mainSizer.Add((5,5),0)            
            mainSizer.Add(ElemSizer())
            
        mainSizer.Add((5,5),0)
        mainSizer.Add(PawleySizer())

        mainSizer.Add((5,5),0)
        mainSizer.Add(MapSizer())

        mainSizer.Add((5,5),0)
        mainSizer.Add(FlipSizer())

        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)

################################################################################
#####  Atom routines
################################################################################

    def FillAtomsGrid():

        G2frame.dataFrame.setSizePosLeft([700,300])
        generalData = data['General']
        atomData = data['Atoms']
        Items = [G2gd.wxID_ATOMSEDITINSERT, G2gd.wxID_ATOMSEDITDELETE, G2gd.wxID_ATOMSREFINE, 
            G2gd.wxID_ATOMSMODIFY, G2gd.wxID_ATOMSTRANSFORM, G2gd.wxID_ATONTESTINSERT]
        if atomData:
            for item in Items:    
                G2frame.dataFrame.AtomsMenu.Enable(item,True)
        else:
            for item in Items:
                G2frame.dataFrame.AtomsMenu.Enable(item,False)            
            
        AAchoice = ": ,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE,HOH,UNK"
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,X,XU,U,F,FX,FXU,FU",]+ \
            3*[wg.GRID_VALUE_FLOAT+':10,5',]+[wg.GRID_VALUE_FLOAT+':10,4', #x,y,z,frac
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+":I,A",]
        Types += 7*[wg.GRID_VALUE_FLOAT+':10,5',]
        colLabels = ['Name','Type','refine','x','y','z','frac','site sym','mult','I/A','Uiso','U11','U22','U33','U12','U13','U23']
        if generalData['Type'] == 'magnetic':
            colLabels += ['Mx','My','Mz']
            Types[2] = wg.GRID_VALUE_CHOICE+": ,X,XU,U,M,MX,MXU,MU,F,FX,FXU,FU,FM,FMX,FMU,"
            Types += 3*[wg.GRID_VALUE_FLOAT+':10,4',]
        elif generalData['Type'] == 'macromolecular':
            colLabels = ['res no','residue','chain'] + colLabels
            Types = [wg.GRID_VALUE_STRING,
                wg.GRID_VALUE_CHOICE+AAchoice,
                wg.GRID_VALUE_STRING] + Types
        elif generalData['Type'] == 'modulated':
            Types += []
            colLabels += []

        def RefreshAtomGrid(event):

            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                for row in range(Atoms.GetNumberRows()):
                    Atoms.SelectRow(row,True)                    
            if r < 0:                          #double click on col label! Change all atoms!
                sel = -1
                noSkip = True
                if Atoms.GetColLabelValue(c) == 'refine':
                    Type = generalData['Type']
                    if Type in ['nuclear','macromolecular']:
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
                if sel >= 0 and noSkip:
                    ui = colLabels.index('U11')
                    us = colLabels.index('Uiso')
                    ss = colLabels.index('site sym')
                    for r in range(Atoms.GetNumberRows()):
                        ID = atomData[r][-1]
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
                        atomData[r][c] = parms
                        if 'Atoms' in data['Drawing']:
                            DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                    FillAtomsGrid()
                    
        def ChangeAtomCell(event):
            
            def chkUij(Uij,CSI): #needs to do something!!!
                return Uij

            r,c =  event.GetRow(),event.GetCol()
            if r >= 0 and c >= 0:
                ID = atomData[r][-1]
                if Atoms.GetColLabelValue(c) in ['x','y','z']:
                    ci = colLabels.index('x')
                    XYZ = atomData[r][ci:ci+3]
                    if None in XYZ:
                        XYZ = [0,0,0]
                    SScol = colLabels.index('site sym')
                    Mulcol = colLabels.index('mult')
                    E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
                    Sytsym,Mult = G2spc.SytSym(XYZ,SGData)
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
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                    FindBondsDraw()
                    
        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if Atoms.GetColLabelValue(c) == 'Type':
                PE = G2elemGUI.PickElement(G2frame)
                if PE.ShowModal() == wx.ID_OK:
                    if PE.Elem not in 'None':                        
                        atomData[r][c] = PE.Elem.strip()
                        name = atomData[r][c]
                        if len(name) in [2,4]:
                            atomData[r][c-1] = name[:2]+'(%d)'%(r+1)
                        else:
                            atomData[r][c-1] = name[:1]+'(%d)'%(r+1)
                PE.Destroy()
                SetupGeneral()
                FillAtomsGrid()
                value = Atoms.GetCellValue(r,c)
                atomData[r][c] = value
                ID = atomData[r][-1]
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                SetupGeneral()
            else:
                event.Skip()

        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                if Atoms.IsSelection():
                    Atoms.ClearSelection()
            elif c < 0:                   #only row clicks
                if event.ControlDown():                    
                    if r in Atoms.GetSelectedRows():
                        Atoms.DeselectRow(r)
                    else:
                        Atoms.SelectRow(r,True)
                elif event.ShiftDown():
                    for row in range(r+1):
                        Atoms.SelectRow(row,True)
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
        
        SGData = data['General']['SGData']
        if SGData['SGPolax']:
            G2frame.dataFrame.SetStatusText('Warning: The location of the origin is arbitrary in '+SGData['SGPolax'])
        table = []
        rowLabels = []
        for i,atom in enumerate(atomData):
            table.append(atom)
            rowLabels.append(str(i))
        atomTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        Atoms.SetTable(atomTable, True)
        Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeAtomCell)
        Atoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, ChangeSelection)
        Atoms.SetMargins(0,0)
        Atoms.AutoSizeColumns(False)
        colType = colLabels.index('Type')
        colSS = colLabels.index('site sym')
        colIA = colLabels.index('I/A')
        colU11 = colLabels.index('U11')
        colUiso = colLabels.index('Uiso')
        for i in range(colU11-1,colU11+6):
            Atoms.SetColSize(i,50)            
        for row in range(Atoms.GetNumberRows()):
            Atoms.SetReadOnly(row,colType,True)
            Atoms.SetReadOnly(row,colSS,True)                         #site sym
            Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            if Atoms.GetCellValue(row,colIA) == 'A':
                CSI = G2spc.GetCSuinel(atomData[row][colLabels.index('site sym')])
                Atoms.SetCellStyle(row,colUiso,VERY_LIGHT_GREY,True)
                Atoms.SetCellTextColour(row,colUiso,VERY_LIGHT_GREY)
                for i in range(6):
                    ci = colU11+i
                    Atoms.SetCellTextColour(row,ci,BLACK)
                    Atoms.SetCellStyle(row,ci,VERY_LIGHT_GREY,True)
                    if CSI[2][i]:
                        Atoms.SetCellStyle(row,ci,WHITE,False)
            else:
                Atoms.SetCellStyle(row,colUiso,WHITE,False)
                Atoms.SetCellTextColour(row,colUiso,BLACK)
                for i in range(6):
                    ci = colU11+i
                    Atoms.SetCellStyle(row,ci,VERY_LIGHT_GREY,True)
                    Atoms.SetCellTextColour(row,ci,VERY_LIGHT_GREY)

    def OnAtomAdd(event):
        AtomAdd(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
        
    def OnAtomTestAdd(event):
        try:
            drawData = data['Drawing']
            x,y,z = drawData['testPos'][0]
            AtomAdd(x,y,z)
        except:
            AtomAdd(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
                
    def AtomAdd(x,y,z,El='H'):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        atId = ran.randint(0,sys.maxint)
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,'UNK','','UNK',El,'',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'nuclear':
            atomData.append(['UNK',El,'',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'magnetic':
            atomData.append(['UNK',El,'',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0,atId])
        SetupGeneral()
        if 'Atoms' in data['Drawing']:            
            DrawAtomAdd(data['Drawing'],atomData[-1])
            G2plt.PlotStructure(G2frame,data)

    def OnAtomInsert(event):
        AtomInsert(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
        
    def OnAtomTestInsert(event):
        if 'Drawing' in data:
            drawData = data['Drawing']
            x,y,z = drawData['testPos'][0]
            AtomAdd(x,y,z)
            FillAtomsGrid()
        event.StopPropagation()
            
    def AtomInsert(x,y,z):
        indx = Atoms.GetSelectedRows()
        if indx:
            indx = indx[0]
            atomData = data['Atoms']
            generalData = data['General']
            Ncol = Atoms.GetNumberCols()
            E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
            Sytsym,Mult = G2spc.SytSym([0,0,0],SGData)
            atId = ran.randint(0,sys.maxint)
            if generalData['Type'] == 'macromolecular':
                atomData.insert(indx,[0,'UNK','','UNK','UNK','',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
            elif generalData['Type'] == 'nuclear':
                atomData.insert(indx,['UNK','UNK','',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
            elif generalData['Type'] == 'magnetic':
                atomData.insert(indx,['UNK','UNK','',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0,atId])
            SetupGeneral()

    def AtomDelete(event):
        indx = Atoms.GetSelectedRows()
        IDs = []
        if indx:
            atomData = data['Atoms']
            indx.reverse()
            for ind in indx:
                atom = atomData[ind]
                IDs.append(atom[-1])
                del atomData[ind]
            if 'Atoms' in data['Drawing']:
                DrawAtomsDeleteByIDs(IDs)
                FillAtomsGrid()
                G2plt.PlotStructure(G2frame,data)
            SetupGeneral()
        event.StopPropagation()

    def AtomRefine(event):
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        c = colLabels.index('refine')
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            generalData = data['General']
            Type = generalData['Type']
            if Type in ['nuclear','macromolecular']:
                choice = ['F - site fraction','X - coordinates','U - thermal parameters']
            elif Type == 'magnetic':
                choice = ['F - site fraction','X - coordinates','U - thermal parameters','M - magnetic moment']
            dlg = wx.MultiChoiceDialog(G2frame,'Select','Refinement controls',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                parms = ''
                for x in sel:
                    parms += choice[x][0]
                for r in indx:
                    atomData[r][c] = parms
                Atoms.ForceRefresh()
            dlg.Destroy()

    def AtomModify(event):                  #intent to implement global modifications (+,-,*,/, etc.)?
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            generalData = data['General']
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            choices = ['Type','x','y','z','frac','I/A','Uiso']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom parameter',choices)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parm = choices[sel]
                cid = colLabels.index(parm)
            dlg.Destroy()
            print parm,cid,indx
            if parm in ['Type']:
                dlg = G2elemGUI.PickElement(G2frame)
                if dlg.ShowModal() == wx.ID_OK:
                    if dlg.Elem not in ['None']:
                        El = dlg.Elem.strip()
                        print El,indx,cid
                        for r in indx:                        
                            atomData[r][cid] = El
                            if len(El) in [2,4]:
                                atomData[r][cid-1] = El[:2]+'(%d)'%(r+1)
                            else:
                                atomData[r][cid-1] = El[:1]+'(%d)'%(r+1)
                        SetupGeneral()
                        if 'Atoms' in data['Drawing']:
                            for r in indx:
                                ID = atomData[r][-1]
                                DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                    FillAtomsGrid()
                dlg.Destroy()
            elif parm in ['I/A']:
                choices = ['Isotropic','Anisotropic']
                dlg = wx.SingleChoiceDialog(G2frame,'Select','Thermal parameter model',choices)
                if dlg.ShowModal() == wx.ID_OK:
                    sel = dlg.GetSelection()
                    parm = choices[sel][0]
                    for r in indx:                        
                        atomData[r][cid] = parm
                    FillAtomsGrid()
            elif parm in ['frac','Uiso']:
                limits = [0.,1.]
                val = 1.0
                if  parm in ['Uiso']:
                    limits = [0.,0.25]
                    val = 0.01
                dlg = SingleFloatDialog(G2frame,'New value','Enter new value for '+parm,val,limits)
                if dlg.ShowModal() == wx.ID_OK:
                    parm = dlg.GetValue()
                    for r in indx:                        
                        atomData[r][cid] = parm
                    SetupGeneral()
                    FillAtomsGrid()
            elif parm in ['x','y','z']:
                limits = [-1.,1.]
                val = 0.
                dlg = SingleFloatDialog(G2frame,'Atom shift','Enter shift for '+parm,val,limits)
                if dlg.ShowModal() == wx.ID_OK:
                    parm = dlg.GetValue()
                    for r in indx:                        
                        atomData[r][cid] += parm
                    SetupGeneral()
                    FillAtomsGrid()

    def AtomTransform(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = colLabels.index('U11')
            css = colLabels.index('site sym')
            atomData = data['Atoms']
            generalData = data['General']
            SGData = generalData['SGData']
            dlg = SymOpDialog(G2frame,SGData,True)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Inv,Cent,Opr,Cell,New = dlg.GetSelection()
                    Cell = np.array(Cell)
                    cent = SGData['SGCen'][Cent]
                    M,T = SGData['SGOps'][Opr]
                    for ind in indx:
                        XYZ = np.array(atomData[ind][cx:cx+3])
                        XYZ = np.inner(M,XYZ)+T
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        if New:
                            atom = copy.copy(atomData[ind])
                        else:
                            atom = atomData[ind]
                        atom[cx:cx+3] = XYZ
                        atom[css:css+2] = G2spc.SytSym(XYZ,SGData)
                        if atom[cuia] == 'A':
                            Uij = atom[cuij:cuij+6]
                            U = G2spc.Uij2U(Uij)
                            U = np.inner(np.inner(M,U),M)
                            Uij = G2spc.U2Uij(U)
                            atom[cuij:cuij+6] = Uij
                        if New:
                            atomData.append(atom)
            finally:
                dlg.Destroy()
            Atoms.ClearSelection()
            if New:
                FillAtomsGrid()
            else:
                Atoms.ForceRefresh()

    def OnDistAngle(event):
        indx = Atoms.GetSelectedRows()
        Oxyz = []
        xyz = []
        DisAglData = {}
        DisAglCtls = {}
        if indx:
            generalData = data['General']
            DisAglData['OrigIndx'] = indx
            if 'DisAglCtls' in generalData:
                DisAglCtls = generalData['DisAglCtls']
            dlg = DisAglDialog(G2frame,DisAglCtls,generalData)
            if dlg.ShowModal() == wx.ID_OK:
                DisAglCtls = dlg.GetData()
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
                DisAglData['covData'] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Covariance'))
            G2str.DistAngle(DisAglCtls,DisAglData)
            
################################################################################
#Structure drawing GUI stuff                
################################################################################

    def SetupDrawingData():
        generalData = data['General']
        atomData = data['Atoms']
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        defaultDrawing = {'Atoms':[],'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,
            'backColor':[0,0,0],'depthFog':False,'Zclip':50.0,'cameraPos':50.,
            'radiusFactor':0.85,'contourLevel':1.,'bondRadius':0.1,'ballScale':0.33,
            'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,'unitCellBox':False,
            'showABC':True,'selectedAtoms':[],'Atoms':[],'Rotation':[0.0,0.0,0.0,[]],
            'bondList':{},'testPos':[[-.1,-.1,-.1],[0.0,0.0,0.0],[0,0]]}
        try:
            drawingData = data['Drawing']
        except KeyError:
            data['Drawing'] = {}
            drawingData = data['Drawing']
        if not drawingData:                 #fill with defaults if empty
            drawingData.update(defaultDrawing)
        if 'contourLevel' not in drawingData:
            drawingData['contourLevel'] = 1.
        cx,ct,cs = [0,0,0]
        if generalData['Type'] == 'nuclear':
            cx,ct,cs = [2,1,6]         #x, type & style
        elif generalData['Type'] == 'macromolecular':
            cx,ct,cs = [5,4,9]         #x, type & style
        elif generalData['Type'] == 'magnetic':
            cx,ct,cs = [2,1,6]         #x, type & style
#        elif generalData['Type'] == 'modulated':
#           ?????   for future
        if not drawingData.get('Atoms'):
            for atom in atomData:
                DrawAtomAdd(drawingData,atom)
            drawingData['atomPtrs'] = [cx,ct,cs]
            data['Drawing'] = drawingData
            
    def MakeDrawAtom(atom,oldatom=None):
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        generalData = data['General']
        SGData = generalData['SGData']
        if generalData['Type'] == 'nuclear':
            if oldatom:
                opr = oldatom[5]
                if atom[9] == 'A':                    
                    X,U = G2spc.ApplyStringOps(opr,SGData,atom[3:6],atom[11:17])
                    atomInfo = [atom[:2]+list(X)+oldatom[5:9]+atom[9:11]+list(U)+oldatom[17:]][0]
                else:
                    X = G2spc.ApplyStringOps(opr,SGData,atom[3:6])
                    atomInfo = [atom[:2]+list(X)+oldatom[5:9]+atom[9:]+[oldatom[-1]]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+['1',]+['vdW balls',]+
                    ['',]+[[255,255,255],]+atom[9:]+[[],[]]][0]
            ct,cs = [1,8]         #type & color
        elif generalData['Type'] == 'macromolecular':
            try:
                oneLetter = AA3letter.index(atom[1])
            except ValueError:
                oneLetter = -1
            atomInfo = [[atom[1].strip()+atom[0],]+
                [AA1letter[oneLetter]+atom[0],]+atom[2:5]+
                atom[6:9]+['1',]+['sticks',]+['',]+[[255,255,255],]+atom[12:]+[[],[]]][0]
            ct,cs = [4,11]         #type & color
        elif generalData['Type'] == 'magnetic':
            if oldatom:
                atomInfo = [atom[:2]+oldatom[3:]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+['vdW balls',]+['',]+atom[9:]+[[],[]]][0]
            ct,cs = [1,8]         #type & color
#        elif generalData['Type'] == 'modulated':
#           ?????   for future
        atNum = generalData['AtomTypes'].index(atom[ct])
        atomInfo[cs] = list(generalData['Color'][atNum])
        return atomInfo
            
    def DrawAtomAdd(drawingData,atom):
        drawingData['Atoms'].append(MakeDrawAtom(atom))
        
    def DrawAtomsReplaceByID(drawingData,atom,ID):
        IDs = [ID,]
        atomData = drawingData['Atoms']
        indx = FindAtomIndexByIDs(atomData,IDs)
        for ind in indx:
            atomData[ind] = MakeDrawAtom(atom,atomData[ind])

################################################################################
##### Atom draw routines
################################################################################
            
    def UpdateDrawAtoms():
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
        cx,ct,cs = drawingData['atomPtrs']
        atomData = drawingData['Atoms']
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
#        elif generalData['Type'] == 'modulated':
#            Types += []
#            colLabels += []

        def RefreshAtomGrid(event):

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
                    noSkip = False
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
                Parms = False
                noSkip = True
                if drawAtoms.GetColLabelValue(c) == 'Style':
                    dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom drawing style',styleChoice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = styleChoice[sel]
                        for r in range(len(atomData)):
                            atomData[r][c] = parms
                            drawAtoms.SetCellValue(r,c,parms)
                        FindBondsDraw()
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
                    FindBondsDraw()
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    color = atomData[r][c]
                    colors = wx.ColourData()
                    colors.SetChooseFull(True)
                    colors.SetCustomColour(0,color)
                    colors.SetColour(color)
                    dlg = wx.ColourDialog(G2frame,colors)
                    dlg.GetColourData().SetCustomColour(0,color)
                    if dlg.ShowModal() == wx.ID_OK:
                        color = dlg.GetColourData().GetColour()
                        attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                        attr.SetReadOnly(True)
                        attr.SetBackgroundColour(color)
                        atomData[r][c] = color
                        drawingData['Atoms'][r][c] = color
                        drawAtoms.SetAttr(i,cs+2,attr)
                    dlg.Destroy()
                    event.StopPropagation()
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
                    for row in range(r+1):
                        drawAtoms.SelectRow(row,True)
                else:
                    drawAtoms.ClearSelection()
                    drawAtoms.SelectRow(r,True)                
            drawingData['selectedAtoms'] = []
            drawingData['selectedAtoms'] = drawAtoms.GetSelectedRows()
            G2plt.PlotStructure(G2frame,data)                    
                
        table = []
        rowLabels = []
        for i,atom in enumerate(drawingData['Atoms']):
            table.append(atom[:colLabels.index('I/A')+1])
            rowLabels.append(str(i))

        atomTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        drawAtoms.SetTable(atomTable, True)
        drawAtoms.SetMargins(0,0)
        drawAtoms.AutoSizeColumns(True)
        drawAtoms.SetColSize(colLabels.index('Style'),80)
        drawAtoms.SetColSize(colLabels.index('Color'),50)
        drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, RefreshAtomGrid)
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
        G2frame.dataFrame.setSizePosLeft([600,300])
        
        FindBondsDraw()
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)

    def DrawAtomStyle(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs = data['Drawing']['atomPtrs']
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
            FindBondsDraw()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)

    def DrawAtomLabel(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs = data['Drawing']['atomPtrs']
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
                G2frame.dataFrame.SetStatusText('Select Custom Color, change color, Add to Custom Colors, then OK')
            else:
                G2frame.dataFrame.SetStatusText('Change color, Add to Custom Colors, then OK')
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs = data['Drawing']['atomPtrs']
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
            for i,color in enumerate(atmColors):
                colors.SetCustomColour(i,color)
            dlg = wx.ColourDialog(G2frame,colors)
            if dlg.ShowModal() == wx.ID_OK:
                for i in range(len(atmColors)):                    
                    atmColors[i] = dlg.GetColourData().GetCustomColour(i)
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
            G2frame.dataFrame.SetStatusText('')
            G2plt.PlotStructure(G2frame,data)
            
    def ResetAtomColors(event):
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        cx,ct,cs = data['Drawing']['atomPtrs']
        for atom in atomData:            
            atNum = generalData['AtomTypes'].index(atom[ct])
            atom[cs+2] = list(generalData['Color'][atNum])
        UpdateDrawAtoms()
        drawAtoms.ClearSelection()
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
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = cuia+2
            atomData = data['Drawing']['Atoms']
            generalData = data['General']
            SGData = generalData['SGData']
            dlg = SymOpDialog(G2frame,SGData,False)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Inv,Cent,Opr,Cell,New = dlg.GetSelection()
                    Cell = np.array(Cell)
                    cent = SGData['SGCen'][Cent]
                    M,T = SGData['SGOps'][Opr]
                    for ind in indx:
                        XYZ = np.array(atomData[ind][cx:cx+3])
                        XYZ = np.inner(M,XYZ)+T
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        if noDuplicate(XYZ,atomData):
                            atom = copy.copy(atomData[ind])
                            atom[cx:cx+3] = XYZ
                            atomOp = atom[cx+3]
                            newOp = str(((Opr+1)+100*Cent)*(1-2*Inv))+'+'+ \
                                str(int(Cell[0]))+','+str(int(Cell[1]))+','+str(int(Cell[2]))                            
                            atom[cx+3] = G2spc.StringOpsProd(atomOp,newOp,SGData)
                            if atom[cuia] == 'A':
                                Uij = atom[cuij:cuij+6]
                                U = G2spc.Uij2U(Uij)
                                U = np.inner(np.inner(M,U),M)
                                Uij = G2spc.U2Uij(U)
                                atom[cuij:cuij+6] = Uij
                            atomData.append(atom)
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
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = cuia+2
            atomData = data['Drawing']['Atoms']
            generalData = data['General']
            SGData = generalData['SGData']
            dlg = SymOpDialog(G2frame,SGData,False)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Inv,Cent,Opr,Cell,New = dlg.GetSelection()
                    Cell = np.array(Cell)
                    cent = SGData['SGCen'][Cent]
                    M,T = SGData['SGOps'][Opr]
                    for ind in indx:
                        XYZ = np.array(atomData[ind][cx:cx+3])
                        XYZ = np.inner(M,XYZ)+T
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        atom = atomData[ind]
                        atom[cx:cx+3] = XYZ
                        atomOp = atom[cx+3]
                        newOp = str(((Opr+1)+100*Cent)*(1-2*Inv))+'+'+ \
                            str(int(Cell[0]))+','+str(int(Cell[1]))+','+str(int(Cell[2]))
                        atom[cx+3] = G2spc.StringOpsProd(atomOp,newOp,SGData)
                        if atom[cuia] == 'A':
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
            cx,ct,cs = data['Drawing']['atomPtrs']
            generalData = data['General']
            SGData = generalData['SGData']
            cellArray = G2lat.CellBlock(1)
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
                                oprB = atomB[cx+3]
                                C = xyz-xyzB
                                newOp = '1+'+str(int(round(C[0])))+','+str(int(round(C[1])))+','+str(int(round(C[2])))
                                newAtom = atomB[:]
                                newAtom[cx:cx+3] = xyz
                                newAtom[cx+3] = G2spc.StringOpsProd(oprB,newOp,SGData)
                                atomData.append(newAtom)
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
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = cuia+2
            generalData = data['General']
            SGData = generalData['SGData']
            for ind in indx:
                atom = atomData[ind]
                XYZ = np.array(atom[cx:cx+3])
                if atom[cuia] == 'A':
                    Uij = atom[cuij:cuij+6]
                    result = G2spc.GenAtom(XYZ,SGData,False,Uij,False)
                    for item in result:
                        atom = copy.copy(atomData[ind])
                        atom[cx:cx+3] = item[0]
                        atom[cx+3] = str(item[2])+'+' \
                            +str(item[3][0])+','+str(item[3][1])+','+str(item[3][2])
                        atom[cuij:cuij+6] = item[1]
                        Opp = G2spc.Opposite(item[0])
                        for xyz in Opp:
                            if noDuplicate(xyz,atomData):
                                cell = np.asarray(np.rint(xyz-atom[cx:cx+3]),dtype=np.int32)
                                cell = '1'+'+'+ \
                                    str(cell[0])+','+str(cell[1])+','+str(cell[2])
                                atom[cx:cx+3] = xyz
                                atom[cx+3] = G2spc.StringOpsProd(cell,atom[cx+3],SGData)
                                atomData.append(atom[:])
                else:
                    result = G2spc.GenAtom(XYZ,SGData,False,Move=False)
                    for item in result:
                        atom = copy.copy(atomData[ind])
                        atom[cx:cx+3] = item[0]
                        atom[cx+3] = str(item[1])+'+' \
                            +str(item[2][0])+','+str(item[2][1])+','+str(item[2][2])
                        Opp = G2spc.Opposite(item[0])
                        for xyz in Opp:
                            if noDuplicate(xyz,atomData):
                                cell = np.asarray(np.rint(xyz-atom[cx:cx+3]),dtype=np.int32)
                                cell = '1'+'+'+ \
                                    str(cell[0])+','+str(cell[1])+','+str(cell[2])
                                atom[cx:cx+3] = xyz
                                atom[cx+3] = G2spc.StringOpsProd(cell,atom[cx+3],SGData)
                                atomData.append(atom[:])               
                data['Drawing']['Atoms'] = atomData
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def FindBondsToo():                         #works but slow for large structures - keep as reference
        cx,ct,cs = data['Drawing']['atomPtrs']
        atomData = data['Drawing']['Atoms']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        radii = generalData['BondRadii']
        atomTypes = generalData['AtomTypes']
        try:
            indH = atomTypes.index('H')
            radii[indH] = 0.5
        except:
            pass            
        for atom in atomData:
            atom[-1] = []
        Atoms = []
        for i,atom in enumerate(atomData):
            Atoms.append([i,np.array(atom[cx:cx+3]),atom[cs],radii[atomTypes.index(atom[ct])]])
        for atomA in Atoms:
            if atomA[2] in ['lines','sticks','ellipsoids','balls & sticks','polyhedra']:
                for atomB in Atoms:                    
                    Dx = atomB[1]-atomA[1]
                    DX = np.inner(Amat,Dx)
                    dist = np.sqrt(np.sum(DX**2))
                    sumR = atomA[3]+atomB[3]
                    if 0.5 < dist <= 0.85*sumR:
                        i = atomA[0]
                        if atomA[2] == 'polyhedra':
                            atomData[i][-1].append(DX)
                        elif atomB[1] != 'polyhedra':
                            j = atomB[0]
                            atomData[i][-1].append(Dx*atomA[3]/sumR)
                            atomData[j][-1].append(-Dx*atomB[3]/sumR)
                    
    def FindBondsDraw():                    #uses numpy & masks - very fast even for proteins!
        import numpy.ma as ma
        cx,ct,cs = data['Drawing']['atomPtrs']
        hydro = data['Drawing']['showHydrogen']
        atomData = data['Drawing']['Atoms']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        radii = generalData['BondRadii']
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
        for atom in atomData:
            Atoms.append(np.array(atom[cx:cx+3]))
            Styles.append(atom[cs])
            try:
                if not hydro and atom[ct] == 'H':
                    Radii.append(0.0)
                else:
                    Radii.append(radii[atomTypes.index(atom[ct])])
            except ValueError:          #changed atom type!
                Radii.append(0.20)
        Atoms = np.array(Atoms)
        Radii = np.array(Radii)
        IASR = zip(Indx,Atoms,Styles,Radii)
        for atomA in IASR:
            if atomA[2] in ['lines','sticks','ellipsoids','balls & sticks','polyhedra']:
                Dx = Atoms-atomA[1]
                dist = ma.masked_less(np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0)),0.5) #gets rid of G2frame & disorder "bonds" < 0.5A
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
        data['Drawing']['Atoms'] = []
        UpdateDrawAtoms()
        drawAtoms.ClearSelection()
        G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()
        
    def FindAtomIndexByIDs(atomData,IDs,Draw=True):
        indx = []
        for i,atom in enumerate(atomData):
            if Draw and atom[-3] in IDs:
                indx.append(i)
            elif atom[-1] in IDs:
                indx.append(i)
        return indx
        
    def DrawAtomsDeleteByIDs(IDs):
        atomData = data['Drawing']['Atoms']
        indx = FindAtomIndexByIDs(atomData,IDs)
        indx.reverse()
        for ind in indx:
            del atomData[ind]
            
    def ChangeDrawAtomsByIDs(colName,IDs,value):
        atomData = data['Drawing']['Atoms']
        cx,ct,cs = data['Drawing']['atomPtrs']
        if colName == 'Name':
            col = ct-1
        elif colName == 'Type':
            col = ct
        elif colName == 'I/A':
            col = cs
        indx = FindAtomIndexByIDs(atomData,IDs)
        for ind in indx:
            atomData[ind][col] = value
                
    def OnDrawPlane(event):
        indx = drawAtoms.GetSelectedRows()
        if len(indx) < 4:
            print '**** ERROR - need 4+ atoms for plane calculation'
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
        G2str.BestPlane(PlaneData)
    
    def OnDrawDAT(event):
        #distance, angle, torsion 
        indx = drawAtoms.GetSelectedRows()
        if len(indx) not in [2,3,4]:
            print '**** ERROR - wrong number of atoms for distance, angle or torsion calculation'
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
            id = FindAtomIndexByIDs(atomData,[atom[cid],],False)[0]
            Oxyz.append([id,]+atomData[id][cx+1:cx+4])
        DATData['Datoms'] = xyz
        DATData['Oatoms'] = Oxyz
        generalData = data['General']
        DATData['Name'] = generalData['Name']
        DATData['SGData'] = generalData['SGData']
        DATData['Cell'] = generalData['Cell'][1:] #+ volume
        if 'pId' in data:
            DATData['pId'] = data['pId']
            DATData['covData'] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Covariance'))
        G2str.DisAglTor(DATData)
                
################################################################################
#### Draw Options page
################################################################################

    def UpdateDrawOptions():
        import copy
        import wx.lib.colourselect as wcs
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
        if generalData['Type'] == 'nuclear':
            pickChoice = ['Atoms','Bonds','Torsions','Planes']
        elif generalData['Type'] == 'macromolecular':
            pickChoice = ['Atoms','Residues','Chains','Bonds','Torsions','Planes','phi/psi']

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
                
            def OnContourLevel(event):
                drawingData['contourLevel'] = contourLevel.GetValue()/100.
                contourLevelTxt.SetLabel(' Contour level: '+'%.2f'%(drawingData['contourLevel']*generalData['Map']['rhoMax']))
                G2plt.PlotStructure(G2frame,data)

            def OnMapSize(event):
                drawingData['mapSize'] = mapSize.GetValue()/10.
                mapSizeTxt.SetLabel(' Map radius, A: '+'%.1f'%(drawingData['mapSize']))
                G2plt.PlotStructure(G2frame,data)

            
            slopSizer = wx.BoxSizer(wx.HORIZONTAL)
            slideSizer = wx.FlexGridSizer(7,2)
            slideSizer.AddGrowableCol(1,1)
    
            cameraPosTxt = wx.StaticText(dataDisplay,-1,
                ' Camera Distance: '+'%.2f'%(drawingData['cameraPos']),name='cameraPos')
            slideSizer.Add(cameraPosTxt,0,wx.ALIGN_CENTER_VERTICAL)
            cameraPos = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['cameraPos'],name='cameraSlider')
            cameraPos.SetRange(10,500)
            cameraPos.Bind(wx.EVT_SLIDER, OnCameraPos)
            slideSizer.Add(cameraPos,1,wx.EXPAND|wx.RIGHT)
            
            ZclipTxt = wx.StaticText(dataDisplay,-1,' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
            slideSizer.Add(ZclipTxt,0,wx.ALIGN_CENTER_VERTICAL)
            Zclip = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['Zclip'])
            Zclip.SetRange(1,99)
            Zclip.Bind(wx.EVT_SLIDER, OnZclip)
            slideSizer.Add(Zclip,1,wx.EXPAND|wx.RIGHT)
            
            vdwScaleTxt = wx.StaticText(dataDisplay,-1,' van der Waals scale: '+'%.2f'%(drawingData['vdwScale']))
            slideSizer.Add(vdwScaleTxt,0,wx.ALIGN_CENTER_VERTICAL)
            vdwScale = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['vdwScale']))
            vdwScale.Bind(wx.EVT_SLIDER, OnVdWScale)
            slideSizer.Add(vdwScale,1,wx.EXPAND|wx.RIGHT)
    
            ellipseProbTxt = wx.StaticText(dataDisplay,-1,' Ellipsoid probability: '+'%d%%'%(drawingData['ellipseProb']))
            slideSizer.Add(ellipseProbTxt,0,wx.ALIGN_CENTER_VERTICAL)
            ellipseProb = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['ellipseProb'])
            ellipseProb.SetRange(1,99)
            ellipseProb.Bind(wx.EVT_SLIDER, OnEllipseProb)
            slideSizer.Add(ellipseProb,1,wx.EXPAND|wx.RIGHT)
    
            ballScaleTxt = wx.StaticText(dataDisplay,-1,' Ball scale: '+'%.2f'%(drawingData['ballScale']))
            slideSizer.Add(ballScaleTxt,0,wx.ALIGN_CENTER_VERTICAL)
            ballScale = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['ballScale']))
            ballScale.Bind(wx.EVT_SLIDER, OnBallScale)
            slideSizer.Add(ballScale,1,wx.EXPAND|wx.RIGHT)
    
            bondRadiusTxt = wx.StaticText(dataDisplay,-1,' Bond radius, A: '+'%.2f'%(drawingData['bondRadius']))
            slideSizer.Add(bondRadiusTxt,0,wx.ALIGN_CENTER_VERTICAL)
            bondRadius = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['bondRadius']))
            bondRadius.SetRange(1,25)
            bondRadius.Bind(wx.EVT_SLIDER, OnBondRadius)
            slideSizer.Add(bondRadius,1,wx.EXPAND|wx.RIGHT)
            
            if generalData['Map']['rhoMax']:
                contourLevelTxt = wx.StaticText(dataDisplay,-1,' Contour level: '+'%.2f'%(drawingData['contourLevel']*generalData['Map']['rhoMax']))
                slideSizer.Add(contourLevelTxt,0,wx.ALIGN_CENTER_VERTICAL)
                contourLevel = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['contourLevel']))
                contourLevel.SetRange(1,100)
                contourLevel.Bind(wx.EVT_SLIDER, OnContourLevel)
                slideSizer.Add(contourLevel,1,wx.EXPAND|wx.RIGHT)
                mapSizeTxt = wx.StaticText(dataDisplay,-1,' Map radius, A: '+'%.1f'%(drawingData['mapSize']))
                slideSizer.Add(mapSizeTxt,0,wx.ALIGN_CENTER_VERTICAL)
                mapSize = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(10*drawingData['mapSize']))
                mapSize.SetRange(1,100)
                mapSize.Bind(wx.EVT_SLIDER, OnMapSize)
                slideSizer.Add(mapSize,1,wx.EXPAND|wx.RIGHT)
            
            slopSizer.Add(slideSizer,1,wx.EXPAND|wx.RIGHT)
            slopSizer.Add((10,5),0)
            slopSizer.SetMinSize(wx.Size(350,10))
            return slopSizer
            
        def ShowSizer():
            
            def OnBackColor(event):
                drawingData['backColor'] = event.GetValue()
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowABC(event):
                drawingData['showABC'] = showABC.GetValue()
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowUnitCell(event):
                drawingData['unitCellBox'] = unitCellBox.GetValue()
                G2plt.PlotStructure(G2frame,data)
    
            def OnShowHyd(event):
                drawingData['showHydrogen'] = showHydrogen.GetValue()
                FindBondsDraw()
                G2plt.PlotStructure(G2frame,data)
                
            showSizer = wx.FlexGridSizer(5,3,5,0)            
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(dataDisplay,-1,' Background color:'),0,wx.ALIGN_CENTER_VERTICAL)
            backColor = wcs.ColourSelect(dataDisplay, -1,colour=drawingData['backColor'],size=wx.Size(25,25))
            backColor.Bind(wcs.EVT_COLOURSELECT, OnBackColor)
            lineSizer.Add(backColor,0,wx.ALIGN_CENTER_VERTICAL)
            showSizer.Add(lineSizer,0,)
            
            showSizer.Add(wx.StaticText(dataDisplay,-1,' View Point:  '),0,wx.ALIGN_CENTER_VERTICAL)
            VP = drawingData['viewPoint'][0]
            viewPoint = wx.TextCtrl(dataDisplay,value='%.3f, %.3f, %.3f'%(VP[0],VP[1],VP[2]),
                style=wx.TE_READONLY,size=wx.Size(120,20),name='viewPoint')
            viewPoint.SetBackgroundColour(VERY_LIGHT_GREY)
            showSizer.Add(viewPoint,0,wx.ALIGN_CENTER_VERTICAL)
            
            showABC = wx.CheckBox(dataDisplay,-1,label=' Show test point?')
            showABC.Bind(wx.EVT_CHECKBOX, OnShowABC)
            showABC.SetValue(drawingData['showABC'])
            showSizer.Add(showABC,0,wx.ALIGN_CENTER_VERTICAL)
    
            unitCellBox = wx.CheckBox(dataDisplay,-1,label=' Show unit cell?')
            unitCellBox.Bind(wx.EVT_CHECKBOX, OnShowUnitCell)
            unitCellBox.SetValue(drawingData['unitCellBox'])
            showSizer.Add(unitCellBox,0,wx.ALIGN_CENTER_VERTICAL)
    
            showHydrogen = wx.CheckBox(dataDisplay,-1,label=' Show hydrogens?')
            showHydrogen.Bind(wx.EVT_CHECKBOX, OnShowHyd)
            showHydrogen.SetValue(drawingData['showHydrogen'])
            showSizer.Add(showHydrogen,0,wx.ALIGN_CENTER_VERTICAL)
            return showSizer
            
        def RadSizer():
            
            def OnSizeHatoms(event):
                try:
                    value = max(0.1,min(1.2,float(sizeH.GetValue())))
                except ValueError:
                    value = 0.5
                drawingData['sizeH'] = value
                sizeH.SetValue("%.2f"%(value))
                G2plt.PlotStructure(G2frame,data)
                
            def OnRadFactor(event):
                try:
                    value = max(0.1,min(1.2,float(radFactor.GetValue())))
                except ValueError:
                    value = 0.85
                drawingData['radiusFactor'] = value
                radFactor.SetValue("%.2f"%(value))
                FindBondsDraw()
                G2plt.PlotStructure(G2frame,data)
            
            radSizer = wx.BoxSizer(wx.HORIZONTAL)
            radSizer.Add(wx.StaticText(dataDisplay,-1,' Hydrogen radius, A:  '),0,wx.ALIGN_CENTER_VERTICAL)
            sizeH = wx.TextCtrl(dataDisplay,-1,value='%.2f'%(drawingData['sizeH']),size=wx.Size(60,20),style=wx.TE_PROCESS_ENTER)
            sizeH.Bind(wx.EVT_TEXT_ENTER,OnSizeHatoms)
            sizeH.Bind(wx.EVT_KILL_FOCUS,OnSizeHatoms)
            radSizer.Add(sizeH,0,wx.ALIGN_CENTER_VERTICAL)
    
            radSizer.Add(wx.StaticText(dataDisplay,-1,' Bond search factor:  '),0,wx.ALIGN_CENTER_VERTICAL)
            radFactor = wx.TextCtrl(dataDisplay,value='%.2f'%(drawingData['radiusFactor']),size=wx.Size(60,20),style=wx.TE_PROCESS_ENTER)
            radFactor.Bind(wx.EVT_TEXT_ENTER,OnRadFactor)
            radFactor.Bind(wx.EVT_KILL_FOCUS,OnRadFactor)
            radSizer.Add(radFactor,0,wx.ALIGN_CENTER_VERTICAL)
            return radSizer

        drawOptions.DestroyChildren()
        dataDisplay = wx.Panel(drawOptions)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,' Drawing controls:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)        
        mainSizer.Add(SlopSizer(),0)
        mainSizer.Add((5,5),0)
        mainSizer.Add(ShowSizer(),0,)
        mainSizer.Add((5,5),0)
        mainSizer.Add(RadSizer(),0,)

        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        G2frame.dataFrame.setSizePosLeft(Size)

################################################################################
####  Texture routines
################################################################################
        
    def UpdateTexture():
        generalData = data['General']        
        SGData = generalData['SGData']
        try:
            textureData = generalData['SH Texture']
        except KeyError:            #fix old files!
            textureData = generalData['SH Texture'] = {'Order':0,'Model':'cylindrical',
                'Sample omega':[False,0.0],'Sample chi':[False,0.0],'Sample phi':[False,0.0],
                'SH Coeff':[False,{}],'SHShow':False,'PFhkl':[0,0,1],
                'PFxyz':[0,0,1.],'PlotType':'Pole figure'}
        if 'SHShow' not in textureData:     #another fix
            textureData.update({'SHShow':False,'PFhkl':[0,0,1],'PFxyz':[0,0,1.],'PlotType':'Pole figure'})
        try:                        #another fix!
            x = textureData['PlotType']
        except KeyError:
            textureData.update({'PFxyz':[0,0,1.],'PlotType':'Pole figure'})
        shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
        SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
        if generalData['doPawley'] and G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Sequental results'):
            G2frame.dataFrame.RefineTexture.Enable(True)
        shAngles = ['omega','chi','phi']
        
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
            wx.CallAfter(UpdateTexture)
            G2plt.PlotTexture(G2frame,data)
                        
        def OnShModel(event):
            Obj = event.GetEventObject()
            textureData['Model'] = Obj.GetValue()
            textureData['SH Coeff'][1] = SetSHCoef()
            wx.CallAfter(UpdateTexture)
            G2plt.PlotTexture(G2frame,data)
            
        def OnSHRefine(event):
            Obj = event.GetEventObject()
            textureData['SH Coeff'][0] = Obj.GetValue()
            
        def OnSHShow(event):
            Obj = event.GetEventObject()
            textureData['SHShow'] = Obj.GetValue()
            wx.CallAfter(UpdateTexture)
            
        def OnProjSel(event):
            Obj = event.GetEventObject()
            G2frame.Projection = Obj.GetValue()
            G2plt.PlotTexture(G2frame,data)
            
        def OnColorSel(event):
            Obj = event.GetEventObject()
            G2frame.ContourColor = Obj.GetValue()
            G2plt.PlotTexture(G2frame,data)
            
        def OnAngRef(event):
            Obj = event.GetEventObject()
            textureData[angIndx[Obj.GetId()]][0] = Obj.GetValue()
            
        def OnAngValue(event):
            Obj = event.GetEventObject()
            try:
                value =  float(Obj.GetValue())
            except ValueError:
                value = textureData[valIndx[Obj.GetId()]][1]
            Obj.SetValue('%8.2f'%(value))
            textureData[valIndx[Obj.GetId()]][1] = value
            
        def OnODFValue(event): 
            Obj = event.GetEventObject()
            try:
                value =  float(Obj.GetValue())
            except ValueError:
                value = textureData['SH Coeff'][1][ODFIndx[Obj.GetId()]]
            Obj.SetValue('%8.3f'%(value))
            textureData['SH Coeff'][1][ODFIndx[Obj.GetId()]] = value
            G2plt.PlotTexture(G2frame,data)
            
        def OnPfType(event):
            Obj = event.GetEventObject()
            textureData['PlotType'] = Obj.GetValue()
            wx.CallAfter(UpdateTexture)
            G2plt.PlotTexture(G2frame,data)
            
        def OnPFValue(event):
            Obj = event.GetEventObject()
            Saxis = Obj.GetValue().split()
            if textureData['PlotType'] in ['Pole figure','Axial pole distribution']:                
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
            G2plt.PlotTexture(G2frame,data)

        if Texture.GetSizer():
            Texture.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        titleSizer = wx.BoxSizer(wx.HORIZONTAL)
        titleSizer.Add(wx.StaticText(Texture,-1,'Spherical harmonics texture data for '+PhaseName+':'),0,wx.ALIGN_CENTER_VERTICAL)
        titleSizer.Add(wx.StaticText(Texture,-1,
            ' Texture Index J = %7.3f'%(G2lat.textureIndex(textureData['SH Coeff'][1]))),
            0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(titleSizer,0)
        mainSizer.Add((0,5),0)
        shSizer = wx.FlexGridSizer(1,6,5,5)
        shSizer.Add(wx.StaticText(Texture,-1,'Texture model: '),0,wx.ALIGN_CENTER_VERTICAL)
        shModel = wx.ComboBox(Texture,-1,value=textureData['Model'],choices=shModels,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        shModel.Bind(wx.EVT_COMBOBOX,OnShModel)
        shSizer.Add(shModel,0,wx.ALIGN_CENTER_VERTICAL)
        shSizer.Add(wx.StaticText(Texture,-1,'  Harmonic order: '),0,wx.ALIGN_CENTER_VERTICAL)
        shOrder = wx.ComboBox(Texture,-1,value=str(textureData['Order']),choices=[str(2*i) for i in range(18)],
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        shOrder.Bind(wx.EVT_COMBOBOX,OnShOrder)
        shSizer.Add(shOrder,0,wx.ALIGN_CENTER_VERTICAL)
        shRef = wx.CheckBox(Texture,-1,label=' Refine texture?')
        shRef.SetValue(textureData['SH Coeff'][0])
        shRef.Bind(wx.EVT_CHECKBOX, OnSHRefine)
        shSizer.Add(shRef,0,wx.ALIGN_CENTER_VERTICAL)
        shShow = wx.CheckBox(Texture,-1,label=' Show coeff.?')
        shShow.SetValue(textureData['SHShow'])
        shShow.Bind(wx.EVT_CHECKBOX, OnSHShow)
        shSizer.Add(shShow,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(shSizer,0,0)
        mainSizer.Add((0,5),0)
        PTSizer = wx.FlexGridSizer(2,4,5,5)
        PTSizer.Add(wx.StaticText(Texture,-1,' Texture plot type: '),0,wx.ALIGN_CENTER_VERTICAL)
        choices = ['Axial pole distribution','Pole figure','Inverse pole figure']            
        pfType = wx.ComboBox(Texture,-1,value=str(textureData['PlotType']),choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        pfType.Bind(wx.EVT_COMBOBOX,OnPfType)
        PTSizer.Add(pfType,0,wx.ALIGN_CENTER_VERTICAL)
        if 'Axial' not in textureData['PlotType']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Projection type: '),0,wx.ALIGN_CENTER_VERTICAL)
            projSel = wx.ComboBox(Texture,-1,value=G2frame.Projection,choices=['equal area','stereographic','3D display'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            projSel.Bind(wx.EVT_COMBOBOX,OnProjSel)
            PTSizer.Add(projSel,0,wx.ALIGN_CENTER_VERTICAL)
        if textureData['PlotType'] in ['Pole figure','Axial pole distribution']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Pole figure HKL: '),0,wx.ALIGN_CENTER_VERTICAL)
            PH = textureData['PFhkl']
            pfVal = wx.TextCtrl(Texture,-1,'%d %d %d'%(PH[0],PH[1],PH[2]),style=wx.TE_PROCESS_ENTER)
        else:
            PTSizer.Add(wx.StaticText(Texture,-1,' Inverse pole figure XYZ: '),0,wx.ALIGN_CENTER_VERTICAL)
            PX = textureData['PFxyz']
            pfVal = wx.TextCtrl(Texture,-1,'%3.1f %3.1f %3.1f'%(PX[0],PX[1],PX[2]),style=wx.TE_PROCESS_ENTER)
        pfVal.Bind(wx.EVT_TEXT_ENTER,OnPFValue)
        pfVal.Bind(wx.EVT_KILL_FOCUS,OnPFValue)
        PTSizer.Add(pfVal,0,wx.ALIGN_CENTER_VERTICAL)
        if 'Axial' not in textureData['PlotType']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Color scheme'),0,wx.ALIGN_CENTER_VERTICAL)
            choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
            choice.sort()
            colorSel = wx.ComboBox(Texture,-1,value=G2frame.ContourColor,choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            colorSel.Bind(wx.EVT_COMBOBOX,OnColorSel)
            PTSizer.Add(colorSel,0,wx.ALIGN_CENTER_VERTICAL)        
        mainSizer.Add(PTSizer,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((0,5),0)
        if textureData['SHShow']:
            mainSizer.Add(wx.StaticText(Texture,-1,'Spherical harmonic coefficients: '),0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((0,5),0)
            ODFSizer = wx.FlexGridSizer(2,8,2,2)
            ODFIndx = {}
            ODFkeys = textureData['SH Coeff'][1].keys()
            ODFkeys.sort()
            for item in ODFkeys:
                ODFSizer.Add(wx.StaticText(Texture,-1,item),0,wx.ALIGN_CENTER_VERTICAL)
                ODFval = wx.TextCtrl(Texture,wx.ID_ANY,'%8.3f'%(textureData['SH Coeff'][1][item]),style=wx.TE_PROCESS_ENTER)
                ODFIndx[ODFval.GetId()] = item
                ODFval.Bind(wx.EVT_TEXT_ENTER,OnODFValue)
                ODFval.Bind(wx.EVT_KILL_FOCUS,OnODFValue)
                ODFSizer.Add(ODFval,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add(ODFSizer,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((0,5),0)
        mainSizer.Add((0,5),0)
        mainSizer.Add(wx.StaticText(Texture,-1,'Sample orientation angles: '),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((0,5),0)
        angSizer = wx.BoxSizer(wx.HORIZONTAL)
        angIndx = {}
        valIndx = {}
        for item in ['Sample omega','Sample chi','Sample phi']:
            angRef = wx.CheckBox(Texture,-1,label=item+': ')
            angRef.SetValue(textureData[item][0])
            angIndx[angRef.GetId()] = item
            angRef.Bind(wx.EVT_CHECKBOX, OnAngRef)
            angSizer.Add(angRef,0,wx.ALIGN_CENTER_VERTICAL)
            angVal = wx.TextCtrl(Texture,wx.ID_ANY,'%8.2f'%(textureData[item][1]),style=wx.TE_PROCESS_ENTER)
            valIndx[angVal.GetId()] = item
            angVal.Bind(wx.EVT_TEXT_ENTER,OnAngValue)
            angVal.Bind(wx.EVT_KILL_FOCUS,OnAngValue)
            angSizer.Add(angVal,0,wx.ALIGN_CENTER_VERTICAL)
            angSizer.Add((5,0),0)
        mainSizer.Add(angSizer,0,wx.ALIGN_CENTER_VERTICAL)
        Texture.SetSizer(mainSizer,True)
        mainSizer.Fit(G2frame.dataFrame)
        Size = mainSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],250) + 20
        Texture.SetSize(Size)
        Texture.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],450)
        G2frame.dataFrame.setSizePosLeft(Size)

################################################################################
##### DData routines
################################################################################
        
    def UpdateDData():
        UseList = data['Histograms']
        if UseList:
            G2frame.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,True)
            G2frame.Refine.Enable(True)
        else:
            G2frame.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,False)
            G2frame.Refine.Enable(False)            
        generalData = data['General']        
        SGData = generalData['SGData']
        keyList = UseList.keys()
        keyList.sort()
        Indx = {}
        
        def PlotSizer():

            def OnPlotSel(event):
                Obj = event.GetEventObject()
                generalData['Data plot type'] = Obj.GetStringSelection()
                wx.CallAfter(UpdateDData)
                G2plt.PlotSizeStrainPO(G2frame,data)
                
            def OnPOhkl(event):
                Obj = event.GetEventObject()
                Saxis = Obj.GetValue().split()
                try:
                    hkl = [int(Saxis[i]) for i in range(3)]
                except (ValueError,IndexError):
                    hkl = generalData['POhkl']
                if not np.any(np.array(hkl)):
                    hkl = generalData['POhkl']
                generalData['POhkl'] = hkl
                h,k,l = hkl
                Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
                G2plt.PlotSizeStrainPO(G2frame,data)
            
            plotSizer = wx.BoxSizer(wx.VERTICAL)
            choice = ['None','Mustrain','Size','Preferred orientation']
            plotSel = wx.RadioBox(DData,-1,'Select plot type:',choices=choice,
                majorDimension=2,style=wx.RA_SPECIFY_COLS)
            plotSel.SetStringSelection(generalData['Data plot type'])
            plotSel.Bind(wx.EVT_RADIOBOX,OnPlotSel)    
            plotSizer.Add(plotSel)
            if generalData['Data plot type'] == 'Preferred orientation':
                POhklSizer = wx.BoxSizer(wx.HORIZONTAL)
                POhklSizer.Add(wx.StaticText(DData,-1,' Plot preferred orientation for H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
                h,k,l = generalData['POhkl']
                poAxis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
                poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOhkl)
                poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOhkl)
                POhklSizer.Add(poAxis,0,wx.ALIGN_CENTER_VERTICAL)
                plotSizer.Add(POhklSizer)            
            return plotSizer
           
        def ScaleSizer():
            
            def OnScaleRef(event):
                Obj = event.GetEventObject()
                UseList[Indx[Obj.GetId()]]['Scale'][1] = Obj.GetValue()
                
            def OnScaleVal(event):
                Obj = event.GetEventObject()
                try:
                    scale = float(Obj.GetValue())
                    if scale > 0:
                        UseList[Indx[Obj.GetId()]]['Scale'][0] = scale
                except ValueError:
                    pass
                Obj.SetValue("%.4f"%(UseList[Indx[Obj.GetId()]]['Scale'][0]))          #reset in case of error
                            
            scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
            scaleRef = wx.CheckBox(DData,-1,label=' Phase fraction: ')
            scaleRef.SetValue(UseList[item]['Scale'][1])
            Indx[scaleRef.GetId()] = item
            scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
            scaleSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
            scaleVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.4f'%(UseList[item]['Scale'][0]),style=wx.TE_PROCESS_ENTER)
            Indx[scaleVal.GetId()] = item
            scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
            scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
            scaleSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
            return scaleSizer
            
        def OnShowData(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['Show'] = Obj.GetValue()
            wx.CallAfter(UpdateDData)
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnCopyData(event):
            #how about HKLF data? This is only for PWDR data
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            sourceDict = UseList[hist]
            copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction']
            copyDict = {}
            for name in copyNames: 
                copyDict[name] = copy.deepcopy(sourceDict[name])        #force copy
            keyList = ['All',]+UseList.keys()
            if UseList:
                copyList = []
                dlg = wx.MultiChoiceDialog(G2frame, 
                    'Copy parameters to which histograms?', 'Copy parameters', 
                    keyList, wx.CHOICEDLG_STYLE)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        result = dlg.GetSelections()
                        for i in result: 
                            copyList.append(keyList[i])
                        if 'All' in copyList: 
                            copyList = keyList[1:]
                        for item in copyList:
                            UseList[item].update(copy.deepcopy(copyDict))
                        wx.CallAfter(UpdateDData)
                finally:
                    dlg.Destroy()
                    
        def OnCopyFlags(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            sourceDict = UseList[hist]
            copyDict = {}
            copyNames = ['Scale','Pref.Ori.','Size','Mustrain','HStrain','Extinction']
            for name in copyNames:
                if name in ['Scale','Extinction','HStrain']:
                    copyDict[name] = sourceDict[name][1]
                elif name in ['Size','Mustrain']:
                    copyDict[name] = [sourceDict[name][0],sourceDict[name][2],sourceDict[name][4]]
                elif name == 'Pref.Ori.':
                    copyDict[name] = [sourceDict[name][0],sourceDict[name][2]]
                    if sourceDict[name][0] == 'SH':
                        SHterms = sourceDict[name][5]
                        SHflags = {}
                        for item in SHterms:
                            SHflags[item] = SHterms[item][1]
                        copyDict[name].append(SHflags)                        
            keyList = ['All',]+UseList.keys()
            if UseList:
                copyList = []
                dlg = wx.MultiChoiceDialog(G2frame, 
                    'Copy parameters to which histograms?', 'Copy parameters', 
                    keyList, wx.CHOICEDLG_STYLE)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        result = dlg.GetSelections()
                        for i in result: 
                            copyList.append(keyList[i])
                        if 'All' in copyList: 
                            copyList = keyList[1:]
                        for item in copyList:
                            UseList[item]                            
                            for name in copyNames:
                                if name in ['Scale','Extinction','HStrain']:
                                    UseList[item][name][1] = copy.copy(copyDict[name])
                                elif name in ['Size','Mustrain']:
                                    UseList[item][name][0] = copy.copy(copyDict[name][0])
                                    UseList[item][name][2] = copy.copy(copyDict[name][1])
                                    UseList[item][name][4] = copy.copy(copyDict[name][2])
                                elif name == 'Pref.Ori.':
                                    UseList[item][name][0] = copy.copy(copyDict[name][0])
                                    UseList[item][name][2] = copy.copy(copyDict[name][1])
                                    if sourceDict[name][0] == 'SH':
                                        SHflags = copy.copy(copyDict[name][2])
                                        SHterms = copy.copy(sourceDict[name][5])
                                        for item in SHflags:
                                            SHterms[item][1] = copy.copy(SHflags[item])                                              
                        wx.CallAfter(UpdateDData)
                finally:
                    dlg.Destroy()
            
            
        def OnLGmixRef(event):
            Obj = event.GetEventObject()
            hist,name = Indx[Obj.GetId()]
            UseList[hist][name][2][2] = Obj.GetValue()
            
        def OnLGmixVal(event):
            Obj = event.GetEventObject()
            hist,name = Indx[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
                if 0.1 <= value <= 1:
                    UseList[hist][name][1][2] = value
                else:
                    raise ValueError
            except ValueError:
                pass
            Obj.SetValue("%.4f"%(UseList[hist][name][1][2]))          #reset in case of error

        def OnSizeType(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['Size'][0] = Obj.GetValue()
            G2plt.PlotSizeStrainPO(G2frame,data)
            wx.CallAfter(UpdateDData)
            
        def OnSizeRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            if UseList[hist]['Size'][0] == 'ellipsoidal':
                UseList[hist]['Size'][5][pid] = Obj.GetValue()                
            else:
                UseList[hist]['Size'][2][pid] = Obj.GetValue()
            
        def OnSizeVal(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            if UseList[hist]['Size'][0] == 'ellipsoidal':
                try:
                    size = float(Obj.GetValue())
                    if pid < 3 and size <= 0.001:            #10A lower limit!
                        raise ValueError
                    UseList[hist]['Size'][4][pid] = size                    
                except ValueError:
                    pass
                Obj.SetValue("%.3f"%(UseList[hist]['Size'][4][pid]))          #reset in case of error
            else:
                try:
                    size = float(Obj.GetValue())
                    if size <= 0.001:            #10A lower limit!
                        raise ValueError
                    UseList[hist]['Size'][1][pid] = size
                except ValueError:
                    pass
                Obj.SetValue("%.3f"%(UseList[hist]['Size'][1][pid]))          #reset in case of error
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnSizeAxis(event):            
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            Saxis = Obj.GetValue().split()
            try:
                hkl = [int(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                hkl = UseList[hist]['Size'][3]
            if not np.any(np.array(hkl)):
                hkl = UseList[hist]['Size'][3]
            UseList[hist]['Size'][3] = hkl
            h,k,l = hkl
            Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
                        
        def OnStrainType(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['Mustrain'][0] = Obj.GetValue()
            wx.CallAfter(UpdateDData)
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnStrainRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            if UseList[hist]['Mustrain'][0] == 'generalized':
                UseList[hist]['Mustrain'][5][pid] = Obj.GetValue()
            else:
                UseList[hist]['Mustrain'][2][pid] = Obj.GetValue()
            
        def OnStrainVal(event):
            Snames = G2spc.MustrainNames(SGData)
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            try:
                strain = float(Obj.GetValue())
                if UseList[hist]['Mustrain'][0] == 'generalized':
                    if '4' in Snames[pid] and strain < 0:
                        raise ValueError
                    UseList[hist]['Mustrain'][4][pid] = strain
                else:
                    if strain <= 0:
                        raise ValueError
                    UseList[hist]['Mustrain'][1][pid] = strain
            except ValueError:
                pass
            if UseList[hist]['Mustrain'][0] == 'generalized':
                Obj.SetValue("%.3f"%(UseList[hist]['Mustrain'][4][pid]))          #reset in case of error
            else:
                Obj.SetValue("%.1f"%(UseList[hist]['Mustrain'][1][pid]))          #reset in case of error
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnStrainAxis(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            Saxis = Obj.GetValue().split()
            try:
                hkl = [int(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                hkl = UseList[hist]['Mustrain'][3]
            if not np.any(np.array(hkl)):
                hkl = UseList[hist]['Mustrain'][3]
            UseList[hist]['Mustrain'][3] = hkl
            h,k,l = hkl
            Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
            G2plt.PlotSizeStrainPO(G2frame,data)
            
        def OnHstrainRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            UseList[hist]['HStrain'][1][pid] = Obj.GetValue()
            
        def OnHstrainVal(event):
            Snames = G2spc.HStrainNames(SGData)
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            try:
                strain = float(Obj.GetValue())
                UseList[hist]['HStrain'][0][pid] = strain
            except ValueError:
                pass
            Obj.SetValue("%.5f"%(UseList[hist]['HStrain'][0][pid]))          #reset in case of error

        def OnPOType(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            if 'March' in Obj.GetValue():
                UseList[hist]['Pref.Ori.'][0] = 'MD'
            else:
                UseList[hist]['Pref.Ori.'][0] = 'SH'
            wx.CallAfter(UpdateDData)            

        def OnPORef(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['Pref.Ori.'][2] = Obj.GetValue()
            
        def OnPOVal(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            try:
                mdVal = float(Obj.GetValue())
                if mdVal > 0:
                    UseList[hist]['Pref.Ori.'][1] = mdVal
            except ValueError:
                pass
            Obj.SetValue("%.3f"%(UseList[hist]['Pref.Ori.'][1]))          #reset in case of error
            
        def OnPOAxis(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            Saxis = Obj.GetValue().split()
            try:
                hkl = [int(Saxis[i]) for i in range(3)]
            except (ValueError,IndexError):
                hkl = UseList[hist]['Pref.Ori.'][3]
            if not np.any(np.array(hkl)):
                hkl = UseList[hist]['Pref.Ori.'][3]
            UseList[hist]['Pref.Ori.'][3] = hkl
            h,k,l = hkl
            Obj.SetValue('%3d %3d %3d'%(h,k,l)) 
            
        def OnPOOrder(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            Order = int(Obj.GetValue())
            UseList[hist]['Pref.Ori.'][4] = Order
            UseList[hist]['Pref.Ori.'][5] = SetPOCoef(Order,hist)
            wx.CallAfter(UpdateDData)

        def SetPOCoef(Order,hist):
            cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',Order,False)     #cylindrical & no M
            newPOCoef = dict(zip(cofNames,np.zeros(len(cofNames))))
            POCoeff = UseList[hist]['Pref.Ori.'][5]
            for cofName in POCoeff:
                if cofName in  cofNames:
                    newPOCoef[cofName] = POCoeff[cofName]
            return newPOCoef
        
        def OnExtRef(event):
            Obj = event.GetEventObject()
            UseList[Indx[Obj.GetId()]]['Extinction'][1] = Obj.GetValue()
            
        def OnExtVal(event):
            Obj = event.GetEventObject()
            try:
                ext = float(Obj.GetValue())
                if ext >= 0:
                    UseList[Indx[Obj.GetId()]]['Extinction'][0] = ext
            except ValueError:
                pass
            Obj.SetValue("%.2f"%(UseList[Indx[Obj.GetId()]]['Extinction'][0]))          #reset in case of error
            
        def checkAxis(axis):
            if not np.any(np.array(axis)):
                return False
            return axis
            
        def TopSizer(name,choices,parm,OnType):
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add(wx.StaticText(DData,-1,name),0,wx.ALIGN_CENTER_VERTICAL)
            sizeType = wx.ComboBox(DData,wx.ID_ANY,value=UseList[item][parm][0],choices=choices,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            sizeType.Bind(wx.EVT_COMBOBOX, OnType)
            Indx[sizeType.GetId()] = item
            topSizer.Add(sizeType)
            topSizer.Add((5,0),0)
            return topSizer
            
        def LGmixSizer(name,OnVal,OnRef):
            lgmixSizer = wx.BoxSizer(wx.HORIZONTAL)
            lgmixRef = wx.CheckBox(DData,-1,label='LGmix')
            lgmixRef.thisown = False
            lgmixRef.SetValue(UseList[item][name][2][2])
            Indx[lgmixRef.GetId()] = [item,name]
            lgmixRef.Bind(wx.EVT_CHECKBOX, OnRef)
            lgmixSizer.Add(lgmixRef,0,wx.ALIGN_CENTER_VERTICAL)
            lgmixVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.4f'%(UseList[item][name][1][2]),style=wx.TE_PROCESS_ENTER)
            Indx[lgmixVal.GetId()] = [item,name]
            lgmixVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
            lgmixVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
            lgmixSizer.Add(lgmixVal,0,wx.ALIGN_CENTER_VERTICAL)
            return lgmixSizer
                        
        def IsoSizer(name,parm,fmt,OnVal,OnRef):
            isoSizer = wx.BoxSizer(wx.HORIZONTAL)
            sizeRef = wx.CheckBox(DData,-1,label=name)
            sizeRef.thisown = False
            sizeRef.SetValue(UseList[item][parm][2][0])
            Indx[sizeRef.GetId()] = [item,0]
            sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
            isoSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
            sizeVal = wx.TextCtrl(DData,wx.ID_ANY,
                fmt%(UseList[item][parm][1][0]),style=wx.TE_PROCESS_ENTER)
            Indx[sizeVal.GetId()] = [item,0]
            sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
            sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
            isoSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
            return isoSizer
            
        def UniSizer(parm,OnAxis):
            uniSizer = wx.BoxSizer(wx.HORIZONTAL)
            uniSizer.Add(wx.StaticText(DData,-1,' Unique axis, H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
            h,k,l = UseList[item][parm][3]
            Axis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
            Indx[Axis.GetId()] = item
            Axis.Bind(wx.EVT_TEXT_ENTER,OnAxis)
            Axis.Bind(wx.EVT_KILL_FOCUS,OnAxis)
            uniSizer.Add(Axis,0,wx.ALIGN_CENTER_VERTICAL)
            return uniSizer
            
        def UniDataSizer(parmName,parm,fmt,OnVal,OnRef):
            dataSizer = wx.BoxSizer(wx.HORIZONTAL)
            parms = zip([' Equatorial '+parmName,' Axial '+parmName],
                UseList[item][parm][1],UseList[item][parm][2],range(2))
            for Pa,val,ref,id in parms:
                sizeRef = wx.CheckBox(DData,-1,label=Pa)
                sizeRef.thisown = False
                sizeRef.SetValue(ref)
                Indx[sizeRef.GetId()] = [item,id]
                sizeRef.Bind(wx.EVT_CHECKBOX, OnRef)
                dataSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
                sizeVal = wx.TextCtrl(DData,wx.ID_ANY,fmt%(val),style=wx.TE_PROCESS_ENTER)
                Indx[sizeVal.GetId()] = [item,id]
                sizeVal.Bind(wx.EVT_TEXT_ENTER,OnVal)
                sizeVal.Bind(wx.EVT_KILL_FOCUS,OnVal)
                dataSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
                dataSizer.Add((5,0),0)
            return dataSizer
            
        def EllSizeDataSizer():
            parms = zip(['S11','S22','S33','S12','S13','S23'],UseList[item]['Size'][4],
                UseList[item]['Size'][5],range(6))
            dataSizer = wx.FlexGridSizer(1,6,5,5)
            for Pa,val,ref,id in parms:
                sizeRef = wx.CheckBox(DData,-1,label=Pa)
                sizeRef.thisown = False
                sizeRef.SetValue(ref)
                Indx[sizeRef.GetId()] = [item,id]
                sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
                dataSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
                sizeVal = wx.TextCtrl(DData,wx.ID_ANY,'%.3f'%(val),style=wx.TE_PROCESS_ENTER)
                Indx[sizeVal.GetId()] = [item,id]
                sizeVal.Bind(wx.EVT_TEXT_ENTER,OnSizeVal)
                sizeVal.Bind(wx.EVT_KILL_FOCUS,OnSizeVal)
                dataSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
            return dataSizer
            
        def GenStrainDataSizer():
            Snames = G2spc.MustrainNames(SGData)
            numb = len(Snames)
            if len(UseList[item]['Mustrain'][4]) < numb:
                UseList[item]['Mustrain'][4] = numb*[0.0,]
                UseList[item]['Mustrain'][5] = numb*[False,]
            parms = zip(Snames,UseList[item]['Mustrain'][4],UseList[item]['Mustrain'][5],range(numb))
            dataSizer = wx.FlexGridSizer(1,6,5,5)
            for Pa,val,ref,id in parms:
                strainRef = wx.CheckBox(DData,-1,label=Pa)
                strainRef.thisown = False
                strainRef.SetValue(ref)
                Indx[strainRef.GetId()] = [item,id]
                strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
                dataSizer.Add(strainRef,0,wx.ALIGN_CENTER_VERTICAL)
                strainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.5f'%(val),style=wx.TE_PROCESS_ENTER)
                Indx[strainVal.GetId()] = [item,id]
                strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
                strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
                dataSizer.Add(strainVal,0,wx.ALIGN_CENTER_VERTICAL)
            return dataSizer
            
        def HstrainSizer():
            hstrainSizer = wx.FlexGridSizer(1,6,5,5)
            Hsnames = G2spc.HStrainNames(SGData)
            parms = zip(Hsnames,UseList[item]['HStrain'][0],UseList[item]['HStrain'][1],range(len(Hsnames)))
            for Pa,val,ref,id in parms:
                hstrainRef = wx.CheckBox(DData,-1,label=Pa)
                hstrainRef.thisown = False
                hstrainRef.SetValue(ref)
                Indx[hstrainRef.GetId()] = [item,id]
                hstrainRef.Bind(wx.EVT_CHECKBOX, OnHstrainRef)
                hstrainSizer.Add(hstrainRef,0,wx.ALIGN_CENTER_VERTICAL)
                hstrainVal = wx.TextCtrl(DData,wx.ID_ANY,'%.5f'%(val),style=wx.TE_PROCESS_ENTER)
                Indx[hstrainVal.GetId()] = [item,id]
                hstrainVal.Bind(wx.EVT_TEXT_ENTER,OnHstrainVal)
                hstrainVal.Bind(wx.EVT_KILL_FOCUS,OnHstrainVal)
                hstrainSizer.Add(hstrainVal,0,wx.ALIGN_CENTER_VERTICAL)
            return hstrainSizer
            
        def PoTopSizer(POData):
            poSizer = wx.FlexGridSizer(1,6,5,5)
            choice = ['March-Dollase','Spherical harmonics']
            POtype = choice[['MD','SH'].index(POData[0])]
            poSizer.Add(wx.StaticText(DData,-1,' Preferred orientation model '),0,wx.ALIGN_CENTER_VERTICAL)
            POType = wx.ComboBox(DData,wx.ID_ANY,value=POtype,choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[POType.GetId()] = item
            POType.Bind(wx.EVT_COMBOBOX, OnPOType)
            poSizer.Add(POType)
            if POData[0] == 'SH':
                poSizer.Add(wx.StaticText(DData,-1,' Harmonic order: '),0,wx.ALIGN_CENTER_VERTICAL)
                poOrder = wx.ComboBox(DData,wx.ID_ANY,value=str(POData[4]),choices=[str(2*i) for i in range(18)],
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                Indx[poOrder.GetId()] = item
                poOrder.Bind(wx.EVT_COMBOBOX,OnPOOrder)
                poSizer.Add(poOrder,0,wx.ALIGN_CENTER_VERTICAL)
                poRef = wx.CheckBox(DData,-1,label=' Refine? ')
                poRef.SetValue(POData[2])
                Indx[poRef.GetId()] = item
                poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
                poSizer.Add(poRef,0,wx.ALIGN_CENTER_VERTICAL)
            return poSizer
           
        def MDDataSizer(POData):
            poSizer = wx.BoxSizer(wx.HORIZONTAL)
            poRef = wx.CheckBox(DData,-1,label=' March-Dollase ratio: ')
            poRef.SetValue(POData[2])
            Indx[poRef.GetId()] = item
            poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
            poSizer.Add(poRef,0,wx.ALIGN_CENTER_VERTICAL)
            poVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.3f'%(POData[1]),style=wx.TE_PROCESS_ENTER)
            Indx[poVal.GetId()] = item
            poVal.Bind(wx.EVT_TEXT_ENTER,OnPOVal)
            poVal.Bind(wx.EVT_KILL_FOCUS,OnPOVal)
            poSizer.Add(poVal,0,wx.ALIGN_CENTER_VERTICAL)
            poSizer.Add(wx.StaticText(DData,-1,' Unique axis, H K L: '),0,wx.ALIGN_CENTER_VERTICAL)
            h,k,l =POData[3]
            poAxis = wx.TextCtrl(DData,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
            Indx[poAxis.GetId()] = item
            poAxis.Bind(wx.EVT_TEXT_ENTER,OnPOAxis)
            poAxis.Bind(wx.EVT_KILL_FOCUS,OnPOAxis)
            poSizer.Add(poAxis,0,wx.ALIGN_CENTER_VERTICAL)
            return poSizer
            
        def SHDataSizer(POData):
            textJ = G2lat.textureIndex(POData[5])
            mainSizer.Add(wx.StaticText(DData,-1,' Spherical harmonic coefficients: '+'Texture index: %.3f'%(textJ)),0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((0,5),0)
            ODFSizer = wx.FlexGridSizer(2,8,2,2)
            ODFIndx = {}
            ODFkeys = POData[5].keys()
            ODFkeys.sort()
            for odf in ODFkeys:
                ODFSizer.Add(wx.StaticText(DData,-1,odf),0,wx.ALIGN_CENTER_VERTICAL)
                ODFval = wx.TextCtrl(DData,wx.ID_ANY,'%8.3f'%(POData[5][odf]),style=wx.TE_PROCESS_ENTER)
                ODFIndx[ODFval.GetId()] = odf
#                ODFval.Bind(wx.EVT_TEXT_ENTER,OnODFValue)
#                ODFval.Bind(wx.EVT_KILL_FOCUS,OnODFValue)
                ODFSizer.Add(ODFval,0,wx.ALIGN_CENTER_VERTICAL)
            return ODFSizer
            
        def ExtSizer():            
            extSizer = wx.BoxSizer(wx.HORIZONTAL)
            extRef = wx.CheckBox(DData,-1,label=' Extinction: ')
            extRef.SetValue(UseList[item]['Extinction'][1])
            Indx[extRef.GetId()] = item
            extRef.Bind(wx.EVT_CHECKBOX, OnExtRef)
            extSizer.Add(extRef,0,wx.ALIGN_CENTER_VERTICAL)
            extVal = wx.TextCtrl(DData,wx.ID_ANY,
                '%.2f'%(UseList[item]['Extinction'][0]),style=wx.TE_PROCESS_ENTER)
            Indx[extVal.GetId()] = item
            extVal.Bind(wx.EVT_TEXT_ENTER,OnExtVal)
            extVal.Bind(wx.EVT_KILL_FOCUS,OnExtVal)
            extSizer.Add(extVal,0,wx.ALIGN_CENTER_VERTICAL)
            return extSizer
            
        if DData.GetSizer():
            DData.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(DData,-1,'Histogram data for '+PhaseName+':'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(PlotSizer())            
            
        for item in keyList:
            histData = UseList[item]
####### Patch to add LGmix to Size & Mustrain
#            if len(histData['Size'][1]) == 2:
#                histData['Size'][1].append(1.0)
#                histData['Size'][2].append(False)
#                histData['Mustrain'][1].append(1.0)
#                histData['Mustrain'][2].append(False)
#                UseList[item] = histData
####### end patch
            showSizer = wx.BoxSizer(wx.HORIZONTAL)
            showData = wx.CheckBox(DData,-1,label=' Show '+item)
            showData.SetValue(UseList[item]['Show'])
            Indx[showData.GetId()] = item
            showData.Bind(wx.EVT_CHECKBOX, OnShowData)
            showSizer.Add(showData,0,wx.ALIGN_CENTER_VERTICAL)
            copyData = wx.Button(DData,-1,label=' Copy?')
            Indx[copyData.GetId()] = item
            copyData.Bind(wx.EVT_BUTTON,OnCopyData)
            showSizer.Add(copyData,wx.ALIGN_CENTER_VERTICAL)
            copyFlags = wx.Button(DData,-1,label=' Copy flags?')
            Indx[copyFlags.GetId()] = item
            copyFlags.Bind(wx.EVT_BUTTON,OnCopyFlags)
            showSizer.Add(copyFlags,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((5,5),0)
            mainSizer.Add(showSizer,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((0,5),0)
            
            if UseList[item]['Show']:
                mainSizer.Add(ScaleSizer())
                mainSizer.Add((0,5),0)
                
            if item[:4] == 'PWDR' and UseList[item]['Show']:
                if UseList[item]['Size'][0] == 'isotropic':
                    isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                    isoSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                        'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                    isoSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                    mainSizer.Add(isoSizer)
                    mainSizer.Add(IsoSizer(u' Cryst. size(\xb5m): ','Size','%.3f',
                        OnSizeVal,OnSizeRef),0,wx.ALIGN_CENTER_VERTICAL)
                elif UseList[item]['Size'][0] == 'uniaxial':
                    uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                    uniSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                        'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                    uniSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                    uniSizer.Add(UniSizer('Size',OnSizeAxis),0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(uniSizer)
                    mainSizer.Add(UniDataSizer(u'size(\xb5m): ','Size','%.3f',OnSizeVal,OnSizeRef))
                elif UseList[item]['Size'][0] == 'ellipsoidal':
                    ellSizer = wx.BoxSizer(wx.HORIZONTAL)
                    ellSizer.Add(TopSizer(' Size model: ',['isotropic','uniaxial','ellipsoidal'],
                        'Size',OnSizeType),0,wx.ALIGN_CENTER_VERTICAL)
                    ellSizer.Add(LGmixSizer('Size',OnLGmixVal,OnLGmixRef))
                    mainSizer.Add(ellSizer)
                    mainSizer.Add(EllSizeDataSizer())
                mainSizer.Add((0,5),0)                    
                
                if UseList[item]['Mustrain'][0] == 'isotropic':
                    isoSizer = wx.BoxSizer(wx.HORIZONTAL)
                    isoSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                        'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                    isoSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                    mainSizer.Add(isoSizer)
                    mainSizer.Add(IsoSizer(' microstrain: ','Mustrain','%.1f',
                        OnStrainVal,OnStrainRef),0,wx.ALIGN_CENTER_VERTICAL)                   
                    mainSizer.Add((0,5),0)
                elif UseList[item]['Mustrain'][0] == 'uniaxial':
                    uniSizer = wx.BoxSizer(wx.HORIZONTAL)
                    uniSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                        'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                    uniSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                    uniSizer.Add(UniSizer('Mustrain',OnStrainAxis),0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(uniSizer)
                    mainSizer.Add(UniDataSizer('mustrain: ','Mustrain','%.1f',OnStrainVal,OnStrainRef))
                elif UseList[item]['Mustrain'][0] == 'generalized':
                    genSizer = wx.BoxSizer(wx.HORIZONTAL)
                    genSizer.Add(TopSizer(' Mustrain model: ',['isotropic','uniaxial','generalized',],
                        'Mustrain',OnStrainType),0,wx.ALIGN_CENTER_VERTICAL)
                    genSizer.Add(LGmixSizer('Mustrain',OnLGmixVal,OnLGmixRef))
                    mainSizer.Add(genSizer)
                    mainSizer.Add(GenStrainDataSizer())                        
                mainSizer.Add((0,5),0)
                
                mainSizer.Add(wx.StaticText(DData,-1,' Hydrostatic/elastic strain:'))
                mainSizer.Add(HstrainSizer())
                    
                #texture  'Pref. Ori.':['MD',1.0,False,[0,0,1],0,[]] last two for 'SH' are SHorder & coeff
                poSizer = wx.BoxSizer(wx.VERTICAL)
                POData = UseList[item]['Pref.Ori.']
                poSizer.Add(PoTopSizer(POData))
                if POData[0] == 'MD':
                    poSizer.Add(MDDataSizer(POData))
                else:           #'SH'
                    if POData[4]:       #SH order > 0
                        poSizer.Add(SHDataSizer(POData))
                        
                mainSizer.Add(poSizer)
                mainSizer.Add((0,5),0)                
                #Extinction  'Extinction':[0.0,False]
                mainSizer.Add(ExtSizer())
                mainSizer.Add((0,5),0)
            elif item[:4] == 'HKLF' and UseList[item]['Show']:
                pass
        mainSizer.Add((5,5),0)

        DData.SetSizer(mainSizer,True)
        mainSizer.FitInside(G2frame.dataFrame)
        Size = mainSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],250) + 20
        DData.SetSize(Size)
        DData.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(Size[1],450)
        G2frame.dataFrame.setSizePosLeft(Size)
        
    def OnHklfAdd(event):
        UseList = data['Histograms']
        keyList = UseList.keys()
        TextList = []
        if G2frame.PatternTree.GetCount():
            item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.PatternTree.GetItemText(item)
                if name not in keyList and 'HKLF' in name:
                    TextList.append(name)
                item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)                        
            dlg = wx.MultiChoiceDialog(G2frame, 'Which new data to use?', 'Use data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result:
                        histoName = TextList[i]
                        UseList[histoName] = {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                            'Extinction':['Lorentzian','Secondary Type I',{'Eg':[0.0,False]},]}                        
                    data['Histograms'] = UseList
                    wx.BeginBusyCursor()
                    UpdateHKLFdata(histoName)
                    wx.EndBusyCursor()
                    wx.CallAfter(UpdateDData)
            finally:
                dlg.Destroy()
                
    def UpdateHKLFdata(histoName):
        generalData = data['General']
        Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,histoName)
        reflData = G2frame.PatternTree.GetItemPyData(Id)
        SGData = generalData['SGData']
        Cell = generalData['Cell'][1:7]
        G,g = G2lat.cell2Gmat(Cell)
        for ref in reflData:
            H = ref[:3]
            ref[4] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
            iabsnt,mulp,Uniq,phi = G2spc.GenHKLf(H,SGData,Friedel=False)
            ref[3] = mulp/2             #convert from powder mulp.
            ref[11] = Uniq
            ref[12] = phi
        G2frame.PatternTree.SetItemPyData(Id,reflData)
        
    def OnPwdrAdd(event):
        generalData = data['General']
        SGData = generalData['SGData']
        UseList = data['Histograms']
        newList = []
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        keyList = UseList.keys()
        TextList = ['All PWDR']
        if G2frame.PatternTree.GetCount():
            item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.PatternTree.GetItemText(item)
                if name not in keyList and 'PWDR' in name:
                    TextList.append(name)
                item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
            dlg = wx.MultiChoiceDialog(G2frame, 'Which new data to use?', 'Use data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: newList.append(TextList[i])
                    if 'All PWDR' in newList:
                        newList = TextList[1:]
                    for histoName in newList:
                        pId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,histoName)
                        UseList[histoName] = {'Histogram':histoName,'Show':False,
                            'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{}],
                            'Size':['isotropic',[4.,4.,1.0],[False,False,False],[0,0,1],
                                [4.,4.,4.,0.,0.,0.],6*[False,]],
                            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                            'HStrain':[NDij*[0.0,],NDij*[False,]],                          
                            'Extinction':[0.0,False]}
                        refList = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,pId,'Reflection Lists'))
                        refList[generalData['Name']] = []                       
                    data['Histograms'] = UseList
                    wx.CallAfter(UpdateDData)
            finally:
                dlg.Destroy()
        
    def OnDataDelete(event):
        UseList = data['Histograms']
        keyList = ['All',]+UseList.keys()
        keyList.sort()
        DelList = []
        if UseList:
            DelList = []
            dlg = wx.MultiChoiceDialog(G2frame, 
                'Which histogram to delete from this phase?', 'Delete histogram', 
                keyList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: 
                        DelList.append(keyList[i])
                    if 'All' in DelList:
                        DelList = keyList[1:]
                    for i in DelList:
                        del UseList[i]
                    wx.CallAfter(UpdateDData)
            finally:
                dlg.Destroy()

################################################################################
##### Pawley routines
################################################################################

    def FillPawleyReflectionsGrid():
                        
        def KeyEditPawleyGrid(event):
            colList = G2frame.PawleyRefl.GetSelectedCols()
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
            
        if 'Pawley ref' in data:
            PawleyPeaks = data['Pawley ref']                        
            rowLabels = []
            for i in range(len(PawleyPeaks)): rowLabels.append(str(i))
            colLabels = ['h','k','l','mul','d','refine','Fsq(hkl)','sig(Fsq)']
            Types = 4*[wg.GRID_VALUE_LONG,]+[wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,]+ \
                2*[wg.GRID_VALUE_FLOAT+':10,2',]
            PawleyTable = G2gd.Table(PawleyPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            G2frame.PawleyRefl.SetTable(PawleyTable, True)
            G2frame.PawleyRefl.Bind(wx.EVT_KEY_DOWN, KeyEditPawleyGrid)                 
            for r in range(G2frame.PawleyRefl.GetNumberRows()):
                for c in range(G2frame.PawleyRefl.GetNumberCols()):
                    if c in [5,6]:
                        G2frame.PawleyRefl.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        G2frame.PawleyRefl.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            G2frame.PawleyRefl.SetMargins(0,0)
            G2frame.PawleyRefl.AutoSizeColumns(False)
            G2frame.dataFrame.setSizePosLeft([500,300])
                    
    def OnPawleyLoad(event):
        generalData = data['General']
        dmin = generalData['Pawley dmin']
        cell = generalData['Cell'][1:7]
        A = G2lat.cell2A(cell)
        SGData = generalData['SGData']
        HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
        PawleyPeaks = []
        wx.BeginBusyCursor()
        try:
            for h,k,l,d in HKLd:
                ext,mul = G2spc.GenHKLf([h,k,l],SGData)[:2]
                if not ext:
                    PawleyPeaks.append([h,k,l,mul,d,False,100.0,1.0])
        finally:
            wx.EndBusyCursor()
        data['Pawley ref'] = PawleyPeaks
        FillPawleyReflectionsGrid()
        
    def OnPawleyEstimate(event):
        try:
            Refs = data['Pawley ref']
            Histograms = data['Histograms']
        except KeyError:
            print '**** Error - no histograms defined for this phase ****'
            return
        HistoNames = Histograms.keys()
        PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,HistoNames[0])
        xdata = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Instrument Parameters'))
        Inst = dict(zip(Inst[3],Inst[1]))
        Sample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Sample Parameters'))
        if 'Lam' in Inst:
            wave = Inst['Lam']
        else:
            wave = Inst['Lam1']
        posCorr = Inst['Zero']
        const = 9.e-2/(np.pi*Sample['Gonio. radius'])                  #shifts in microns
        
        for ref in Refs:
            pos = 2.0*asind(wave/(2.0*ref[4]))
            if 'Bragg' in Sample['Type']:
                pos -= const*(4.*Sample['Shift'][0]*cosd(pos/2.0)+ \
                    Sample['Transparency'][0]*sind(pos)*100.0)            #trans(=1/mueff) in cm
            else:               #Debye-Scherrer - simple but maybe not right
                pos -= const*(Sample['DisplaceX'][0]*cosd(pos)+Sample['DisplaceY'][0]*sind(pos))
            indx = np.searchsorted(xdata[0],pos)
            try:
                ref[6] = xdata[1][indx]/ref[3]
                pola,dIdPola = G2pwd.Polarization(Inst['Polariz.'],xdata[0][indx],0.0)
                ref[6] /= pola
            except IndexError:
                pass
        FillPawleyReflectionsGrid()
                            
    def OnPawleyDelete(event):
        dlg = wx.MessageDialog(G2frame,'Do you really want to delete Pawley reflections?','Delete', 
            wx.YES_NO | wx.ICON_QUESTION)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_YES: 
            data['Pawley ref'] = []
            FillPawleyReflectionsGrid()

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
                    for row in range(r+1):
                        MapPeaks.SelectRow(row,True)
                else:
                    MapPeaks.ClearSelection()
                    MapPeaks.SelectRow(r,True)                
            G2plt.PlotStructure(G2frame,data)                    
            
        G2frame.dataFrame.setSizePosLeft([450,300])
        if 'Map Peaks' in data:
            mapPeaks = data['Map Peaks']                        
            rowLabels = []
            for i in range(len(mapPeaks)): rowLabels.append(str(i))
            colLabels = ['mag','x','y','z']
            Types = 4*[wg.GRID_VALUE_FLOAT+':10,4',]
            MapPeaksTable = G2gd.Table(mapPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            MapPeaks.SetTable(MapPeaksTable, True)
            MapPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
            for r in range(MapPeaks.GetNumberRows()):
                for c in range(MapPeaks.GetNumberCols()):
                    MapPeaks.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            MapPeaks.SetMargins(0,0)
            MapPeaks.AutoSizeColumns(False)
                    
    def OnPeaksMove(event):
        if 'Map Peaks' in data:
            mapPeaks = data['Map Peaks']                        
            Ind = MapPeaks.GetSelectedRows()
            for ind in Ind:
                print mapPeaks[ind]
                x,y,z = mapPeaks[ind][1:]
                AtomAdd(x,y,z,'C')
    
    def OnPeaksClear(event):
        data['Map Peaks'] = []
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)
    
    def OnFourierMaps(event):
        generalData = data['General']
        mapData = generalData['Map']
        reflName = mapData['RefList']
        phaseName = generalData['Name']
        if 'PWDR' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflData = reflSets[phaseName]
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflData = G2frame.PatternTree.GetItemPyData(PatternId)
        mapData.update(G2mth.FourierMap(data,reflData))
        mapData['Flip'] = False
        mapSig = np.std(mapData['rho'])
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        data['Drawing']['contourLevel'] = 1.
        data['Drawing']['mapSize'] = 10.
        print mapData['MapType']+' computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig)
        UpdateDrawAtoms()
        G2plt.PlotStructure(G2frame,data)
        
    def printRho(SGLaue,rho,rhoMax):                          
# map printing for testing purposes
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
                print line+'\n'
        else:
            ix,jy,kz = rho.shape
            for k in range(kz):
                print 'k = ',k
                for j in range(jy):
                    line = ''
                    if SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                        line += (jy-j)*'  '
                    for i in range(ix):
                        r = int(100*rho[i,j,k]/rhoMax)
                        line += '%4d'%(r)
                    print line+'\n'
## keep this                
    
    def OnSearchMaps(event):
        
        peaks = []
        mags = []
        print ' Begin fourier map search - can take some time'
        time0 = time.time()
        wx.BeginBusyCursor()
        try:
            peaks,mags = G2mth.SearchMap(data,keepDup=True)
        finally:
            wx.EndBusyCursor()
        sortIdx = np.argsort(mags.flatten())
        if len(peaks):
            data['Map Peaks'] = np.concatenate((mags,peaks),axis=1)
            
            print ' Map search peaks found:'
            print '  No.    Mag.      x        y        z'
            for j,i in enumerate(sortIdx[::-1]):
                print ' %3d %8.3f %8.5f %8.5f %8.5f'%(j,mags[i],peaks[i][0],peaks[i][1],peaks[i][2])
            print ' Map search finished, time = %.2fs'%(time.time()-time0)
        Page = G2frame.dataDisplay.FindPage('Map peaks')
        G2frame.dataDisplay.ChangeSelection(Page)
        G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.MapPeaksMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksMove, id=G2gd.wxID_PEAKSMOVE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksClear, id=G2gd.wxID_PEAKSCLEAR)
        UpdateDrawAtoms()
        FillMapPeaksGrid()
        G2plt.PlotStructure(G2frame,data)
        
    def OnChargeFlip(event):
        generalData = data['General']
        mapData = generalData['Map']
        flipData = generalData['Flip']
        reflName = flipData['RefList']
        phaseName = generalData['Name']
        if 'PWDR' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflData = reflSets[phaseName]
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflData = G2frame.PatternTree.GetItemPyData(PatternId)
        else:
            print '**** ERROR - No data defined for charge flipping'
            return
        pgbar = wx.ProgressDialog('Charge flipping','Residual Rcf =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = pgbar.GetSize()
        Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
        pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        pgbar.SetSize(Size)
        try:
            mapData.update(G2mth.ChargeFlip(data,reflData,pgbar))
        finally:
            pgbar.Destroy()
        mapData['Flip'] = True        
        mapSig = np.std(mapData['rho'])
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        data['Drawing']['contourLevel'] = 1.
        data['Drawing']['mapSize'] = 10.
        print ' Charge flip map computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig)
        if mapData['Rcf'] < 99.:
            OnSearchMaps(event)             #does a plot structure at end
        else:
            print 'Bad charge flip map - no peak search done'
                
    def OnTextureRefine(event):
        print 'refine texture?'
        event.Skip()        
            
    def OnTextureClear(event):
        print 'clear texture?'
        event.Skip()

    def OnPageChanged(event):
        page = event.GetSelection()
        text = G2frame.dataDisplay.GetPageText(page)
        if text == 'Atoms':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.AtomsMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomAdd, id=G2gd.wxID_ATOMSEDITADD)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomTestAdd, id=G2gd.wxID_ATOMSTESTADD)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomInsert, id=G2gd.wxID_ATOMSEDITINSERT)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomTestInsert, id=G2gd.wxID_ATONTESTINSERT)
            G2frame.dataFrame.Bind(wx.EVT_MENU, AtomDelete, id=G2gd.wxID_ATOMSEDITDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, AtomRefine, id=G2gd.wxID_ATOMSREFINE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, AtomModify, id=G2gd.wxID_ATOMSMODIFY)
            G2frame.dataFrame.Bind(wx.EVT_MENU, AtomTransform, id=G2gd.wxID_ATOMSTRANSFORM)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnReloadDrawAtoms, id=G2gd.wxID_RELOADDRAWATOMS)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDistAngle, id=G2gd.wxID_ATOMSDISAGL)
            FillAtomsGrid()
        elif text == 'General':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataGeneral)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnFourierMaps, id=G2gd.wxID_FOURCALC)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnSearchMaps, id=G2gd.wxID_FOURSEARCH)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChargeFlip, id=G2gd.wxID_CHARGEFLIP)
            UpdateGeneral()
        elif text == 'Data':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2gd.wxID_PWDRADD)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2gd.wxID_HKLFADD)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDataDelete, id=G2gd.wxID_DATADELETE)
            UpdateDData()
            G2plt.PlotSizeStrainPO(G2frame,data,Start=True)
        elif text == 'Draw Options':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataDrawOptions)
            UpdateDrawOptions()
            G2plt.PlotStructure(G2frame,data)
        elif text == 'Draw Atoms':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DrawAtomsMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, DrawAtomStyle, id=G2gd.wxID_DRAWATOMSTYLE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, DrawAtomLabel, id=G2gd.wxID_DRAWATOMLABEL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, DrawAtomColor, id=G2gd.wxID_DRAWATOMCOLOR)
            G2frame.dataFrame.Bind(wx.EVT_MENU, ResetAtomColors, id=G2gd.wxID_DRAWATOMRESETCOLOR)
            G2frame.dataFrame.Bind(wx.EVT_MENU, SetViewPoint, id=G2gd.wxID_DRAWVIEWPOINT)
            G2frame.dataFrame.Bind(wx.EVT_MENU, AddSymEquiv, id=G2gd.wxID_DRAWADDEQUIV)
            G2frame.dataFrame.Bind(wx.EVT_MENU, TransformSymEquiv, id=G2gd.wxID_DRAWTRANSFORM)
            G2frame.dataFrame.Bind(wx.EVT_MENU, FillCoordSphere, id=G2gd.wxID_DRAWFILLCOORD)            
            G2frame.dataFrame.Bind(wx.EVT_MENU, FillUnitCell, id=G2gd.wxID_DRAWFILLCELL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, DrawAtomsDelete, id=G2gd.wxID_DRAWDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDrawDAT, id=G2gd.wxID_DRAWDISAGLTOR)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDrawPlane, id=G2gd.wxID_DRAWPLANE)
            UpdateDrawAtoms()
            G2plt.PlotStructure(G2frame,data)
        elif text == 'Pawley reflections':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.PawleyMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyLoad, id=G2gd.wxID_PAWLEYLOAD)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyEstimate, id=G2gd.wxID_PAWLEYESTIMATE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyDelete, id=G2gd.wxID_PAWLEYDELETE)            
            FillPawleyReflectionsGrid()
        elif text == 'Texture':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.TextureMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnTextureRefine, id=G2gd.wxID_REFINETEXTURE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnTextureClear, id=G2gd.wxID_CLEARTEXTURE)
            UpdateTexture()                        
            G2plt.PlotTexture(G2frame,data,Start=True)
        elif text == 'Map peaks':
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.MapPeaksMenu)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksMove, id=G2gd.wxID_PEAKSMOVE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksClear, id=G2gd.wxID_PEAKSCLEAR)
            FillMapPeaksGrid()
            G2plt.PlotStructure(G2frame,data)
            
        else:
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
        event.Skip()
        
    General = wx.Window(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(General,'General')
    G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.DataGeneral)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnFourierMaps, id=G2gd.wxID_FOURCALC)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSearchMaps, id=G2gd.wxID_FOURSEARCH)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnChargeFlip, id=G2gd.wxID_CHARGEFLIP)
    SetupGeneral()
    GeneralData = data['General']
    UpdateGeneral()

    DData = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(DData,'Data')
    Atoms = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(Atoms,'Atoms')
    drawOptions = wx.Window(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(drawOptions,'Draw Options')
    drawAtoms = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(drawAtoms,'Draw Atoms')
    Texture = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(Texture,'Texture')
    MapPeaks = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(MapPeaks,'Map peaks')
    G2frame.PawleyRefl = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(G2frame.PawleyRefl,'Pawley reflections')
            
    G2frame.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    G2frame.dataDisplay.SetSelection(oldPage)
    
            
